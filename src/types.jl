"""
    data_eltype(x::AbstractArray)
    data_eltype(T::Type{<:Number})

The element type a transform computes in for data `x` (or of element type `T`):
`float(T)`, so complex stays complex and integers become `Float64`.

Filters are always real; transforms take `Td = data_eltype(x)` for the data and
`Tf = filter_eltype(Td)` for the filters, so complex data is filtered by real
kernels (complex × real, never complex × complex).  That keeps every transform
exactly real-linear.
"""
data_eltype(x::AbstractArray) = float(eltype(x))
data_eltype(::Type{T}) where {T <: Number} = float(T)

"""
    filter_eltype(Td::Type{<:Number})

The real filter precision for data of element type `Td`: `real(float(Td))`
(`ComplexF32` → `Float32`).  See [`data_eltype`](@ref).
"""
filter_eltype(::Type{Td}) where {Td <: Number} = real(float(Td))

"""
    FilterPair{T}

A 1-D biorthogonal filter pair holding an analysis (`h`) and synthesis (`g`) low-pass
filter vector.  Both filters are stored as coefficient vectors whose centre tap is at
index `ceil(Int, length/2)` (1-based Julia indexing).

The pair drives the pyramid stage of the downstream transforms: `h` produces the
coarse approximation, and `g` predicts the input from it, so that the difference
is a bandpass (the Laplacian pyramid of Burt and Adelson [1]).  Adding the
prediction back reconstructs the input for any pair.  Do and Vetterli [2] show
that the pyramid is a frame, tight for orthogonal filters, and give a projection
reconstruction that is more robust to noise than the usual one.  Its nonsubsampled (à trous)
form is the pyramid of the nonsubsampled contourlet transform [3].

# Examples
```jldoctest
julia> using MultiscaleFilterBanks: FilterPair

julia> fp = FilterPair([0.5, 1.0, 0.5], [0.25, 0.5, 0.25]);

julia> eltype(fp)
Float64
```

# References
1. P. J. Burt and E. H. Adelson, "The Laplacian Pyramid as a Compact Image Code",
   *IEEE Transactions on Communications*, 31(4):532–540, 1983.
   [doi:10.1109/TCOM.1983.1095851](https://doi.org/10.1109/TCOM.1983.1095851)
2. M. N. Do and M. Vetterli, "Framing Pyramids", *IEEE Transactions on Signal
   Processing*, 51(9):2329–2342, 2003.
   [doi:10.1109/TSP.2003.815389](https://doi.org/10.1109/TSP.2003.815389)
3. A. L. da Cunha, J. Zhou, and M. N. Do, "The Nonsubsampled Contourlet Transform:
   Theory, Design, and Applications", *IEEE Transactions on Image Processing*,
   15(10):3089–3101, 2006.
   [doi:10.1109/TIP.2006.877507](https://doi.org/10.1109/TIP.2006.877507)
"""
struct FilterPair{T <: AbstractFloat}
    h::Vector{T}   # analysis (low-pass)
    g::Vector{T}   # synthesis (low-pass)
end

FilterPair(h::AbstractVector, g::AbstractVector) =
    FilterPair(promote_type(eltype(h), eltype(g)).(h), promote_type(eltype(h), eltype(g)).(g))

Base.eltype(::FilterPair{T}) where {T} = T
