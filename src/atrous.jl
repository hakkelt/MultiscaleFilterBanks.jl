# À trous filters of the nonsubsampled pyramid: level j filters with h and g
# upsampled by 2^(j-1) and does not decimate (Holschneider et al. 1989; Shensa
# 1992, see `upsample_filter`).
#
# The biorthogonal synthesis `g` has DC gain 2 to compensate the ×2 upsampling of
# the *decimated* Laplacian pyramid.  The nonsubsampled pyramid does not decimate,
# so `g` is rescaled to unit DC gain (`g/2` per axis).  Without the rescale the
# prediction is about twice the input on smooth regions, and the "bandpass" is
# dominated by lowpass leakage.  Both the 2-D and the 3-D pyramid rely on it being
# applied exactly once, here.
const _NSP_SYNTH_SCALE = 0.5

"""
    nsp_atrous_filters(fp::FilterPair, level, Tf) -> (h, g)

The à trous analysis and synthesis filters of nonsubsampled-pyramid level `level`
in precision `Tf`: `fp.h` and `fp.g` upsampled by `2^(level - 1)`
([`upsample_filter`](@ref)), with the synthesis filter halved to unit DC gain.
The pyramid does not decimate, so the gain of 2 that compensates the ×2
upsampling of the decimated Laplacian pyramid must go.  Apply the pair along
every axis.
"""
function nsp_atrous_filters(fp::FilterPair, level::Int, ::Type{Tf}) where {Tf}
    factor = 2^(level - 1)
    h_up = upsample_filter(fp.h, factor)
    h_j = Tf === eltype(h_up) ? h_up : Tf.(h_up)
    g_up = upsample_filter(fp.g, factor)
    g_j = Tf(_NSP_SYNTH_SCALE) .* (Tf === eltype(g_up) ? g_up : Tf.(g_up))
    return h_j, g_j
end

"""
    centred_dtft(f, θ) -> Complex

The discrete-time Fourier transform at `θ` of the filter `f` with its centre tap at
`(length(f) + 1) ÷ 2`: the frequency response of the centred convolution the
pyramids apply.
"""
function centred_dtft(f::AbstractVector, θ::Real)
    c = (length(f) + 1) ÷ 2
    return sum(f[l] * cis(-θ * (l - c)) for l in eachindex(f))
end

"""
    nsp_responses(fp::FilterPair, J, n; bins = 1:n) -> (a, gh, s)

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
function nsp_responses(fp::FilterPair, J::Int, n::Int; bins = 1:n)
    a = [ones(ComplexF64, length(bins))]
    s = [ones(ComplexF64, length(bins))]
    gh = Vector{ComplexF64}[]
    for j in 1:J
        h, g = nsp_atrous_filters(fp, j, Float64)
        H = [centred_dtft(h, 2π * (k - 1) / n) for k in bins]
        G = [centred_dtft(g, 2π * (k - 1) / n) for k in bins]
        push!(gh, G .* H)
        push!(a, a[end] .* H)
        push!(s, s[end] .* G)
    end
    return a, gh, s
end
