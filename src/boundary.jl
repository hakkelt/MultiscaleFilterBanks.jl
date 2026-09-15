# Boundary extension helpers shared by the 2-D and 3-D convolution primitives.
#
# All helpers take the boundary mode `B` as a `Val` type parameter so the compiler
# emits fully specialised code with no dispatch overhead per call.
# The symmetric helper uses an in-range fast path (skips mod for interior
# samples) followed by a branch-free abs-fold for out-of-range indices.

@inline function _clamp_idx(i::Int, n::Int, ::Val{:symmetric})::Int
    1 <= i <= n && return i
    n == 1 && return 1
    p = 2 * (n - 1)
    m = mod(i - 1, p)
    return n - abs(m - (n - 1))
end

@inline function _clamp_idx(i::Int, n::Int, ::Val{:periodic})::Int
    return mod(i - 1, n) + 1
end

@inline function _clamp_idx(i::Int, n::Int, ::Val{:zero})::Int
    return i
end
