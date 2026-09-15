# À trous (with holes) filter construction for the nonsubsampled pyramid.
#
# Synthesis-filter normalisation: the analysis low-pass `h` has unit DC gain, but
# the biorthogonal synthesis `g` has DC gain 2 because it is designed to
# compensate the ×2 upsampling of the *decimated* Laplacian pyramid.  The
# nonsubsampled pyramid does no decimation, so `g` is rescaled to unit DC gain
# (`g/2` per axis).  Without this the "bandpass" would be dominated by low-pass
# leakage (prediction ≈ 2× the image on smooth regions) instead of being a true
# bandpass — defeating the directional decomposition that follows.
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
