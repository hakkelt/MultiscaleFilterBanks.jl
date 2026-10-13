"""
    upsample_filter(h::AbstractVector, factor::Int) -> Vector

Insert `factor - 1` zeros between every tap of `h`.  If `h` has frequency
response ``H(ω)``, the result has response ``H(\\text{factor}\\,ω)``.

This is the filter upsampling of the *algorithme à trous* ("with holes") [1]:
level ``j`` of an undecimated pyramid filters with the level-1 filters upsampled by
``2^{j-1}`` instead of decimating the signal, which computes the undecimated
(shift-invariant) wavelet transform [2].  The nonsubsampled pyramid and the
nonsubsampled directional filter banks of the downstream packages use it in the
same way [3].

# Examples
```jldoctest
julia> using MultiscaleFilterBanks: upsample_filter

julia> upsample_filter([1.0, 2.0, 3.0], 2)
5-element Vector{Float64}:
 1.0
 0.0
 2.0
 0.0
 3.0
```

# References
1. M. Holschneider, R. Kronland-Martinet, J. Morlet, and Ph. Tchamitchian, "A
   Real-Time Algorithm for Signal Analysis with the Help of the Wavelet Transform",
   in *Wavelets: Time-Frequency Methods and Phase Space*, Springer, 1989,
   pp. 286–297. [doi:10.1007/978-3-642-75988-8_28](https://doi.org/10.1007/978-3-642-75988-8_28)
2. M. J. Shensa, "The Discrete Wavelet Transform: Wedding the À Trous and Mallat
   Algorithms", *IEEE Transactions on Signal Processing*, 40(10):2464–2482, 1992.
   [doi:10.1109/78.157290](https://doi.org/10.1109/78.157290)
3. A. L. da Cunha, J. Zhou, and M. N. Do, "The Nonsubsampled Contourlet Transform:
   Theory, Design, and Applications", *IEEE Transactions on Image Processing*,
   15(10):3089–3101, 2006.
   [doi:10.1109/TIP.2006.877507](https://doi.org/10.1109/TIP.2006.877507)
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
