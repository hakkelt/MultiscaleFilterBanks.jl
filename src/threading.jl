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
function _foreach_chunk(f::F, n::Int, nt::Int, threaded::Bool) where {F}
    if threaded && nt > 1
        Threads.@threads :dynamic for t in 1:nt
            f(_chunk(n, nt, t), t)
        end
    else
        f(1:n, 1)
    end
    return nothing
end

# A sum over `1:n` split like `_foreach_chunk`: `f(a, range, t)` adds its chunk's
# terms into `a`, which is `acc` for task 1 and the zeroed `partials[t - 1]` for
# the others; the partials are then added into `acc`.  `length(partials) + 1` is
# the task count.  `foreach` runs the chunks; a package with a cheaper task
# runner (Polyester's `@batch` allocates a tenth of `Threads.@threads`) passes a
# function with `_foreach_chunk`'s signature.
function _chunked_sum!(
        f::F, acc::AbstractArray, partials::AbstractVector, n::Int, threaded::Bool;
        foreach::G = _foreach_chunk
    ) where {F, G}
    nt = length(partials) + 1
    foreach(n, nt, threaded) do range, t
        a = t == 1 ? acc : fill!(partials[t - 1], zero(eltype(acc)))
        f(a, range, t)
    end
    if threaded && nt > 1
        for a in partials
            acc .+= a
        end
    end
    return acc
end
