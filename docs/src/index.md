# MultiscaleFilterBanks.jl

Shared building blocks for multiscale, multidirectional filter-bank transforms.

This package is the common core of two transform packages:

- [Contourlets.jl](https://github.com/hakkelt/Contourlets.jl): the 2-D
  contourlet and nonsubsampled contourlet transforms.
- [Surfacelets.jl](https://github.com/hakkelt/Surfacelets.jl): the 3-D
  nonsubsampled and decimated surfacelet transforms.

It holds only what is neither 2-D nor 3-D specific: the biorthogonal pyramid
filter pair ([`FilterPair`](@ref), [`CDF97`](@ref)), à trous filter upsampling
([`upsample_filter`](@ref)) and the threading policy types
([`ThreadingPolicy`](@ref)). Most users never load it directly — both
downstream packages re-export its public names.

It has no non-stdlib dependencies.

```jldoctest
julia> using MultiscaleFilterBanks

julia> length(CDF97.h), length(CDF97.g)
(9, 7)

julia> upsample_filter([1.0, 2.0, 3.0], 2)
5-element Vector{Float64}:
 1.0
 0.0
 2.0
 0.0
 3.0
```

## References

1. A. Cohen, I. Daubechies, and J.-C. Feauveau, "Biorthogonal Bases of Compactly
   Supported Wavelets", *Communications on Pure and Applied Mathematics*,
   45(5):485–560, 1992 — the CDF 9/7 pair, here with the JPEG 2000 Part 1
   (ISO/IEC 15444-1) coefficients.
2. M. Holschneider, R. Kronland-Martinet, J. Morlet, and Ph. Tchamitchian, "A
   Real-Time Algorithm for Signal Analysis with the Help of the Wavelet
   Transform", in *Wavelets, Time-Frequency Methods and Phase Space*, Springer,
   1989, pp. 286–297 — the à trous algorithm.
3. M. J. Shensa, "The Discrete Wavelet Transform: Wedding the À Trous and Mallat
   Algorithms", *IEEE Transactions on Signal Processing*, 40(10):2464–2482, 1992.
