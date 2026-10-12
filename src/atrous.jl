# À trous (with holes) filter construction for the nonsubsampled pyramid.
#
# Level j filters with h and g upsampled by 2^(j-1) (zeros inserted between the
# taps) and does not decimate: the "algorithme à trous" of Holschneider et al.
# (1989), shown by Shensa (1992) to compute the undecimated wavelet transform.
#
# Synthesis-filter normalisation: the analysis low-pass `h` has unit DC gain, but
# the biorthogonal synthesis `g` has DC gain 2 because it is designed to
# compensate the ×2 upsampling of the *decimated* Laplacian pyramid.  The
# nonsubsampled pyramid does no decimation, so `g` is rescaled to unit DC gain
# (`g/2` per axis).  Without this the "bandpass" would be dominated by low-pass
# leakage (prediction ≈ 2× the image on smooth regions) instead of being a true
# bandpass — defeating the directional decomposition that follows.
#
# References:
#   M. Holschneider, R. Kronland-Martinet, J. Morlet & Ph. Tchamitchian, "A
#     real-time algorithm for signal analysis with the help of the wavelet
#     transform", in Wavelets, Time-Frequency Methods and Phase Space, Springer,
#     1989, pp. 286–297.
#   M. J. Shensa, "The discrete wavelet transform: wedding the à trous and Mallat
#     algorithms", IEEE Trans. Signal Process. 40(10), 1992, pp. 2464–2482.
const _NSP_SYNTH_SCALE = 0.5

# à trous analysis/synthesis filter pair for pyramid `level`, in real precision
# `Tf`.  The returned `g_j` already carries the synthesis rescale above.  Nothing
# here is dimension-specific: the 2-D pyramid applies the pair per axis, and so
# does the 3-D one.
function _nsp_atrous_filters(fp::FilterPair, level::Int, ::Type{Tf}) where {Tf}
    factor = 2^(level - 1)
    h_up = upsample_filter(fp.h, factor)
    h_j = Tf === eltype(h_up) ? h_up : Tf.(h_up)
    g_up = upsample_filter(fp.g, factor)
    g_j = Tf(_NSP_SYNTH_SCALE) .* (Tf === eltype(g_up) ? g_up : Tf.(g_up))
    return h_j, g_j
end

# DTFT at θ of a filter whose centre tap sits at `(length + 1) ÷ 2`: the
# frequency response of the periodic, centred convolution the pyramids apply.
function _centred_dtft(f::AbstractVector, θ::Real)
    c = (length(f) + 1) ÷ 2
    return sum(f[l] * cis(-θ * (l - c)) for l in eachindex(f))
end

"""
    _nsp_responses(fp, J, n; bins=1:n) -> (a, gh, s)

The `J`-level nonsubsampled pyramid on one `n`-point DFT axis, as frequency
responses at the bins `bins` (`ComplexF64` vectors).  With Hⱼ, Gⱼ the responses
of the level-j à trous pair,

    a[j] = Π_{i<j} Hᵢ  (j ≤ J + 1),   gh[j] = Gⱼ Hⱼ  (j ≤ J),   s[j] = Π_{i<j} Gᵢ  (j ≤ J + 1).

Periodic separable filtering is a pointwise product on the DFT grid, so on a
D-dimensional grid the bandpass of level j is `Π_d a_d[j] · (1 − Π_d gh_d[j])`,
the lowpass `Π_d a_d[J + 1]`, and their synthesis masks `Π_d s_d[j]`; with them
`Σⱼ synthesis · analysis = 1` telescopes.  Each mask is a broadcast of per-axis
vectors, never stored at full size.
"""
function _nsp_responses(fp::FilterPair, J::Int, n::Int; bins = 1:n)
    a = [ones(ComplexF64, length(bins))]
    s = [ones(ComplexF64, length(bins))]
    gh = Vector{ComplexF64}[]
    for j in 1:J
        h, g = _nsp_atrous_filters(fp, j, Float64)
        H = [_centred_dtft(h, 2π * (k - 1) / n) for k in bins]
        G = [_centred_dtft(g, 2π * (k - 1) / n) for k in bins]
        push!(gh, G .* H)
        push!(a, a[end] .* H)
        push!(s, s[end] .* G)
    end
    return a, gh, s
end
