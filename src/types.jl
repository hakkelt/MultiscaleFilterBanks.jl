# Types shared across the multiscale filter-bank packages.
# FilterPair: 1-D biorthogonal filter pair (used for the Laplacian / nonsubsampled
# pyramid stages of both the 2-D contourlet and the 3-D surfacelet transforms).
#
# Element-type convention: filters are always real (`Tf <: AbstractFloat`); image
# data may be real or complex (`Td <: Number`).  The transforms thread the two
# independently — `Td = _data_eltype(image)` and `Tf = _filter_eltype(Td)` — so
# complex images are filtered by real kernels (complex·real, not complex·complex).

# Data (sample) element type: preserves complex, promotes integers to float
# (Int → Float64, Float32 → Float32, ComplexF32 → ComplexF32).
_data_eltype(x::AbstractArray) = float(eltype(x))
_data_eltype(::Type{T}) where {T <: Number} = float(T)

# Real filter precision matching a data element type (real part of the float
# type): Float64 → Float64, ComplexF64 → Float64, ComplexF32 → Float32.
_filter_eltype(::Type{Td}) where {Td <: Number} = real(float(Td))

"""
    FilterPair{T}

A 1-D biorthogonal filter pair holding an analysis (`h`) and synthesis (`g`) low-pass
filter vector.  Both filters are stored as coefficient vectors whose centre tap is at
index `ceil(Int, length/2)` (1-based Julia indexing).

# Examples
```jldoctest
julia> using MultiscaleFilterBanks

julia> fp = FilterPair([0.5, 1.0, 0.5], [0.25, 0.5, 0.25]);

julia> eltype(fp)
Float64
```
"""
struct FilterPair{T <: AbstractFloat}
    h::Vector{T}   # analysis (low-pass)
    g::Vector{T}   # synthesis (low-pass)
end

FilterPair(h::AbstractVector, g::AbstractVector) =
    FilterPair(promote_type(eltype(h), eltype(g)).(h), promote_type(eltype(h), eltype(g)).(g))

Base.eltype(::FilterPair{T}) where {T} = T
