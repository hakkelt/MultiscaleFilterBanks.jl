"""
    CDF97 :: FilterPair{Float64}

The Cohen–Daubechies–Feauveau 9/7 biorthogonal filter pair [1], the irreversible
wavelet filter of JPEG 2000 Part 1 [3], introduced for image coding by Antonini et
al. [2].  The 9-tap analysis low-pass `h` is normalised to unit DC gain, the 7-tap
synthesis low-pass `g` to DC gain 2, which compensates the ×2 upsampling of the
Laplacian pyramid's synthesis step.  Both are symmetric, centred at taps 5 and 4.

```jldoctest
julia> using MultiscaleFilterBanks: CDF97

julia> isapprox(sum(CDF97.h), 1.0; atol = 1e-12)
true

julia> isapprox(sum(CDF97.g), 2.0; atol = 1e-12)
true
```

# References
1. A. Cohen, I. Daubechies, and J.-C. Feauveau, "Biorthogonal Bases of Compactly
   Supported Wavelets", *Communications on Pure and Applied Mathematics*,
   45(5):485–560, 1992.
   [doi:10.1002/cpa.3160450502](https://doi.org/10.1002/cpa.3160450502)
2. M. Antonini, M. Barlaud, P. Mathieu, and I. Daubechies, "Image Coding Using
   Wavelet Transform", *IEEE Transactions on Image Processing*, 1(2):205–220, 1992.
   [doi:10.1109/83.136597](https://doi.org/10.1109/83.136597)
3. ISO/IEC 15444-1, *Information technology — JPEG 2000 image coding system: Core
   coding system*.
"""
const CDF97 = FilterPair{Float64}(
    # analysis low-pass h
    Float64[
        0.02674875741081,
        -0.016864118442875,
        -0.07822326652899,
        0.266864118442875,
        0.60294901823636,
        0.266864118442875,
        -0.07822326652899,
        -0.016864118442875,
        0.02674875741081,
    ],
    # synthesis low-pass g
    Float64[
        -0.09127176311425,
        -0.0575435262285,
        0.59127176311425,
        1.115087052456994,
        0.59127176311425,
        -0.0575435262285,
        -0.09127176311425,
    ],
)
