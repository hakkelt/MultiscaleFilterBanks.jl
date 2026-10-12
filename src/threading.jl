# ── Threading Policy API ────────────────────────────────────────────────────────

"""
    ThreadingPolicy

Abstract type for configuring the threading behavior of transforms.
Available policies are:
- [`Auto`](@ref)
- [`Enabled`](@ref)
- [`Disabled`](@ref)
"""
abstract type ThreadingPolicy end

"""
    Auto <: ThreadingPolicy

Automatically selects the optimal threading strategy based on the data type.
By default, multithreading is disabled for `Real` data (which uses SIMD) and 
enabled for `Complex` data.
"""
struct Auto <: ThreadingPolicy end

"""
    Enabled <: ThreadingPolicy

Explicitly enables multithreading for loops.
"""
struct Enabled <: ThreadingPolicy end

"""
    Disabled <: ThreadingPolicy

Explicitly disables multithreading.
"""
struct Disabled <: ThreadingPolicy end

# Resolving whether to use threading based on type:
@inline _use_threading(::Enabled, ::Type{T}) where {T} = true
@inline _use_threading(::Disabled, ::Type{T}) where {T} = false
@inline function _use_threading(::Auto, ::Type{T}) where {T}
    # By default, disabled for Real, enabled for Complex.
    return !(T <: Real)
end

# ── Task partitioning ───────────────────────────────────────────────────────────

# Tasks a loop over independent FFTs can split into: one per thread on the host,
# one on a device (whose FFTs are already parallel).
_n_fft_tasks(::Array) = Threads.nthreads()
_n_fft_tasks(::AbstractArray) = 1

# Items of the `t`-th of `nt` contiguous, balanced chunks of `1:n`.
_chunk(n::Int, nt::Int, t::Int) = (div((t - 1) * n, nt) + 1):div(t * n, nt)

# Number of tasks for a loop over `n` independent FFTs on arrays like `X`.
_n_tasks(X::AbstractArray, n::Int) = max(1, min(_n_fft_tasks(X), n))

# The threaded loop over independent FFTs (a filter bank's subbands): each of
# the `nt` tasks runs a contiguous chunk with its own buffers, indexed by `t`.
# FFTW's own threads split one FFT and gain little at these sizes; this splits
# the FFTs.  Calls `f(range, t)` per chunk, threaded when `threaded`; otherwise
# `f(1:n, 1)` once.  Buffers are indexed by the chunk `t`, never by thread id, so
# they need no locks and the tasks may migrate.
#
# Two runners, chosen per call (`_batch_safe`).  Polyester's `@batch` starts its
# tasks in under a microsecond where `Threads.@threads` takes 3–10 µs, and
# allocates ~150 B instead of 2–4 KB: standalone it is 5–25% faster (NSCT 128²
# and TICT 256² on 8 threads, NS surfacelets 64³ on 4).  But its workers are
# pinned to threads 2…nthreads() and wait for those threads to be free, so
# called from inside a caller's own `@threads` loop it ran ~100× slower than
# `@threads` (nested NSCT 64²: 43 ms against 0.46 ms), and a busy second
# interactive thread stalls it outright.  So `@batch` runs only from the
# interactive thread with no second one; anywhere else the cooperative scheduler
# runs the tasks.
#
# An exception escaping a `@batch` body leaves Polyester's workers unusable, so
# every later `@batch` in the process silently runs serially; the body therefore
# catches it, and the first one is rethrown once all tasks are done.
function _foreach_chunk(f::F, n::Int, nt::Int, threaded::Bool) where {F}
    if threaded && nt > 1
        failed = Ref{Any}(nothing)
        if _batch_safe()
            @batch for t in 1:nt
                _guarded_chunk(f, failed, n, nt, t)
            end
        else
            Threads.@threads :dynamic for t in 1:nt
                _guarded_chunk(f, failed, n, nt, t)
            end
        end
        # Either way the caller sees the task's own exception, not a wrapper.
        failed[] === nothing || throw(failed[])
    else
        f(1:n, 1)
    end
    return nothing
end

# Thread budget of these loops, registered with NestedThreading as the counted pool
# `:multiscalefilterbanks` (`__init__`): a caller's `@budgeted_threads` loop or
# `with_restricted_threads` scope lowers it to 1, and the loops then run serially,
# which is the fastest nested option (a caller's 8-way loop over NSCT 64²: 0.37 ms
# serial inside against 0.46 ms with `@threads` and 43 ms with `@batch`).  A
# transform asks once per call, with `_loop_threaded`, and passes the answer to
# every loop of the call, so a budget changing mid-call cannot split it.
const _LOOP_THREADS = Threads.Atomic{Int}(1)

_get_loop_threads() = _LOOP_THREADS[]
_set_loop_threads(n::Integer) = (_LOOP_THREADS[] = max(1, Int(n)); nothing)

_register_loop_pool() = NestedThreading.register_counted_pool!(
    _get_loop_threads, _set_loop_threads; name = :multiscalefilterbanks
)

# Whether a call asked to thread (`threaded`) may: no enclosing NestedThreading
# restriction has taken its budget down to one thread.
_loop_threaded(threaded::Bool) = threaded && _LOOP_THREADS[] > 1

# Whether `@batch` may run: the caller is on the interactive pool (the main task
# from Julia 1.12, not a task of a caller's threaded loop, which run on the
# default pool), and that pool has one thread, so Polyester's workers (threads
# 2…nthreads()) are all default threads.
_batch_safe() = Threads.nthreads(:interactive) == 1 && Threads.threadpool() === :interactive

# Chunk `t`, recording the first exception in `failed` instead of letting it reach
# Polyester.  (A function of its own: `@batch` mis-captures a `catch` variable.)
function _guarded_chunk(f::F, failed::Ref{Any}, n::Int, nt::Int, t::Int) where {F}
    try
        f(_chunk(n, nt, t), t)
    catch err
        @lock _CHUNK_ERROR_LOCK (failed[] === nothing && (failed[] = err))
    end
    return nothing
end

const _CHUNK_ERROR_LOCK = ReentrantLock()

# A sum over `1:n` split like `_foreach_chunk`: `f(a, range, t)` adds its chunk's
# terms into `a`, which is `acc` for task 1 and the zeroed `partials[t - 1]` for
# the others; the partials are then added into `acc` (`_add_partials!`).
# `length(partials) + 1` is the task count.
function _chunked_sum!(f::F, acc::AbstractArray, partials::AbstractVector, n::Int, threaded::Bool) where {F}
    nt = length(partials) + 1
    _foreach_chunk(n, nt, threaded) do range, t
        a = t == 1 ? acc : fill!(partials[t - 1], zero(eltype(acc)))
        f(a, range, t)
    end
    threaded && nt > 1 && _add_partials!(acc, partials, true)
    return acc
end

# acc += Σ partials.  On the host each task takes a contiguous slice of `acc` and
# adds every partial into it a cache-sized block at a time, so `acc` is read and
# written once; adding whole partials in turn was the serial tail of the threaded
# inverses (eight tasks at 96³: ~60 ms of 180).  A slice is at least 8192 elements.
_add_partials!(acc::AbstractArray, partials::AbstractVector, threaded::Bool) = (foreach(a -> (acc .+= a), partials); acc)

function _add_partials!(acc::Array, partials::AbstractVector{<:Array}, threaded::Bool)
    isempty(partials) && return acc
    n = length(acc)
    nt = threaded ? max(1, min(Threads.nthreads(), n >> 13)) : 1
    _foreach_chunk(n, nt, threaded) do range, _
        for lo in first(range):_ADD_BLOCK:last(range)
            block = lo:min(lo + _ADD_BLOCK - 1, last(range))
            for a in partials
                @inbounds @simd for i in block
                    acc[i] += a[i]
                end
            end
        end
    end
    return acc
end

const _ADD_BLOCK = 2048

# ── FFTW thread budget ──────────────────────────────────────────────────────────

# Run `f()`, which makes FFTW plans, with FFTW on `n` threads.  An FFTW plan keeps
# the thread count it was made with, so this scopes the planning.  It goes through
# NestedThreading rather than `FFTW.set_num_threads`: the caller's count comes
# back afterwards, and a restriction opened around the call (`@budgeted_threads`,
# `with_restricted_threads`) still wins, so plans made inside a caller's threaded
# loop do not oversubscribe it.  Without FFTW loaded the scope does nothing.
#
# `with_thread_budget` itself infers as `Any` (its guards are called through
# `Function`), which would leave the plans untyped downstream, so the result goes
# through a `Ref` typed by `f`'s inferred return type.
function _with_fftw_threads(f::F, n::Integer) where {F}
    out = Ref{Base.promote_op(f)}()
    NestedThreading.with_thread_budget(n; only = (:fftw,)) do
        out[] = f()
        nothing
    end
    return out[]
end
