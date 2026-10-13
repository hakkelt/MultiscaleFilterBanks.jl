# ── Threading keyword ──────────────────────────────────────────────────────────

"""
    use_threading(threaded::Union{Bool, Nothing}, T::Type) -> Bool

Resolve a transform's `threaded` keyword for data of element type `T`: `true`
and `false` are taken as given; `nothing` (let the transform decide) threads
complex data only, since a vectorised serial loop wins for real data in the
direct-convolution transforms.  Transforms whose work is dominated by FFTs may
resolve `nothing` to `true` instead.

```jldoctest
julia> using MultiscaleFilterBanks: use_threading

julia> use_threading(nothing, Float64), use_threading(nothing, ComplexF64), use_threading(true, Float64)
(false, true, true)
```
"""
function use_threading end

@inline use_threading(threaded::Bool, ::Type) = threaded
@inline use_threading(::Nothing, ::Type{T}) where {T} = !(T <: Real)

# ── Task partitioning ───────────────────────────────────────────────────────────

# Tasks a loop over independent FFTs can split into: one per thread on the host,
# one on a device (whose FFTs are already parallel).
_n_fft_tasks(::Array) = Threads.nthreads()
_n_fft_tasks(::AbstractArray) = 1

# Items of the `t`-th of `nt` contiguous, balanced chunks of `1:n`.
_chunk(n::Int, nt::Int, t::Int) = (div((t - 1) * n, nt) + 1):div(t * n, nt)

"""
    n_tasks(X::AbstractArray, n::Int) -> Int

The number of tasks a loop over `n` independent FFTs on arrays like `X` splits
into: one per Julia thread for an `Array`, one for any other array type (a device
array, whose FFTs are already parallel), and never more than `n`.  Allocate one
set of per-task buffers for each.
"""
n_tasks(X::AbstractArray, n::Int) = max(1, min(_n_fft_tasks(X), n))

"""
    foreach_chunk(f, n::Int, nt::Int, threaded::Bool)

The threaded loop over independent work items, such as the FFTs of a filter
bank's subbands.  Calls `f(range, t)` for each of `nt` contiguous, balanced chunks
of `1:n` in parallel when `threaded`, otherwise `f(1:n, 1)` once.  Index per-task
buffers by the chunk number `t`, never by thread id: they then need no locks and
the tasks may migrate.  An exception thrown by `f` is rethrown after all chunks
finish.  Call it with `threaded = loop_threaded(...)` (see [`loop_threaded`](@ref)).
"""
function foreach_chunk end

# Splitting the FFTs scales better than FFTW's own threads at these sizes.
#
# Polyester's `@batch` starts tasks faster than `Threads.@threads`, but its
# workers are pinned to threads 2…nthreads(): called from inside a caller's own
# threaded loop it is orders of magnitude slower, and a busy second interactive
# thread stalls it.  So `@batch` runs only where `_batch_safe()`, `@threads`
# everywhere else.
#
# An exception escaping a `@batch` body leaves Polyester's workers unusable (every
# later `@batch` in the process silently runs serially), so the body catches it
# and the first one is rethrown once all tasks are done.
function foreach_chunk(f::F, n::Int, nt::Int, threaded::Bool) where {F}
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

# Thread budget of these loops, registered with NestedThreading as the counted
# pool `:multiscalefilterbanks` (`__init__`).  A caller's `@budgeted_threads` loop
# or `with_restricted_threads` scope lowers it to 1 and the loops run serially,
# the fastest option when nested.  A transform reads it once per call
# (`loop_threaded`), so a budget changing mid-call cannot split the call.
const _LOOP_THREADS = Threads.Atomic{Int}(1)

_get_loop_threads() = _LOOP_THREADS[]
_set_loop_threads(n::Integer) = (_LOOP_THREADS[] = max(1, Int(n)); nothing)

_register_loop_pool() = NestedThreading.register_counted_pool!(
    _get_loop_threads, _set_loop_threads; name = :multiscalefilterbanks
)

"""
    loop_threaded(threaded::Bool) -> Bool

Whether a call that asked to thread may: `threaded`, and no enclosing
NestedThreading.jl restriction (a caller's `@budgeted_threads` loop or
`with_restricted_threads` scope) has lowered this package's thread budget to one.
Read it once per transform call and pass the answer to every loop of the call.
"""
loop_threaded(threaded::Bool) = threaded && _LOOP_THREADS[] > 1

# Whether `@batch` may run: the caller is on the interactive pool (the main task
# from Julia 1.12, not a task of a caller's threaded loop) and that pool has one
# thread, so Polyester's workers (threads 2…nthreads()) are all default threads.
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

"""
    chunked_sum!(f, acc, partials, n::Int, threaded::Bool) -> acc

A sum over `1:n` split like [`foreach_chunk`](@ref) into `length(partials) + 1`
tasks: `f(a, range, t)` adds the terms of its chunk into `a`, which is `acc` for
the first task and the zeroed `partials[t - 1]` for the others.  The partials are
then added into `acc` ([`add_partials!`](@ref)).
"""
function chunked_sum!(f::F, acc::AbstractArray, partials::AbstractVector, n::Int, threaded::Bool) where {F}
    nt = length(partials) + 1
    foreach_chunk(n, nt, threaded) do range, t
        a = t == 1 ? acc : fill!(partials[t - 1], zero(eltype(acc)))
        f(a, range, t)
    end
    threaded && nt > 1 && add_partials!(acc, partials, true)
    return acc
end

"""
    add_partials!(acc, partials, threaded::Bool) -> acc

`acc .+= sum(partials)`.  For host arrays each task takes a contiguous slice of
`acc` and adds every partial into it one cache-sized block at a time, so `acc` is
read and written once.
"""
function add_partials! end

add_partials!(acc::AbstractArray, partials::AbstractVector, threaded::Bool) = (foreach(a -> (acc .+= a), partials); acc)

function add_partials!(acc::Array, partials::AbstractVector{<:Array}, threaded::Bool)
    isempty(partials) && return acc
    n = length(acc)
    nt = threaded ? max(1, min(Threads.nthreads(), n >> 13)) : 1
    foreach_chunk(n, nt, threaded) do range, _
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

"""
    with_fftw_threads(f, n::Integer)

Run `f()`, which makes FFTW plans, with FFTW on `n` threads, and return its
result; a plan keeps the thread count it was made with.  Unlike
`FFTW.set_num_threads`, it restores the caller's count afterwards and lets an
enclosing NestedThreading.jl restriction win, so plans made inside a caller's
threaded loop do not oversubscribe it.  Without FFTW loaded it just runs `f()`.
"""
function with_fftw_threads end

# `with_thread_budget` infers as `Any`, so the result goes through a `Ref` typed
# by `f`'s inferred return type to keep the plans typed downstream.
function with_fftw_threads(f::F, n::Integer) where {F}
    out = Ref{Base.promote_op(f)}()
    NestedThreading.with_thread_budget(n; only = (:fftw,)) do
        out[] = f()
        nothing
    end
    return out[]
end
