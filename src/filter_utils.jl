# Utilities operating on FilterPair.

"""
    upsample_filter(h::AbstractVector, factor::Int) -> Vector

Insert `factor - 1` zeros between every tap of `h`.  This implements the
"à trous" (with holes) upsampling used by the nonsubsampled pyramid and the
nonsubsampled directional filter banks to avoid downsampling at coarser levels.

# Examples
```jldoctest
julia> using MultiscaleFilterBanks

julia> upsample_filter([1.0, 2.0, 3.0], 2)
5-element Vector{Float64}:
 1.0
 0.0
 2.0
 0.0
 3.0
```
"""
function upsample_filter(h::AbstractVector{T}, factor::Int)::Vector{T} where {T}
    factor >= 1 || throw(ArgumentError("factor must be ≥ 1"))
    factor == 1 && return copy(h)
    n = length(h)
    out = zeros(T, factor * (n - 1) + 1)
    @inbounds for k in 1:n
        out[(k - 1) * factor + 1] = h[k]
    end
    return out
end
