"""
    MultiscaleFilterBanks

Shared building blocks for multiscale, multidirectional filter-bank transforms.

This package holds the pieces that are neither 2-D nor 3-D specific and are used
by both [Contourlets.jl](https://github.com/hakkelt/Contourlets.jl) (the 2-D
contourlet / nonsubsampled contourlet transform) and
[Surfacelets.jl](https://github.com/hakkelt/Surfacelets.jl) (the 3-D
nonsubsampled and decimated surfacelet transforms):

| Export | Description |
|--------|-------------|
| `FilterPair` | 1-D biorthogonal analysis/synthesis low-pass filter pair |
| `CDF97` | The CDF 9/7 filter pair used by the pyramid stages |
| `upsample_filter` | À trous ("with holes") filter upsampling |
| `ThreadingPolicy`, `Auto`, `Enabled`, `Disabled` | Threading policy types |

It also carries the internal helpers both downstream packages share: the
element-type convention (`_data_eltype` / `_filter_eltype`), the boundary
index folding (`_clamp_idx`) and the à trous pyramid filter construction
(`_nsp_atrous_filters`).

This package has no non-stdlib dependencies.
"""
module MultiscaleFilterBanks

include("threading.jl")
include("types.jl")
include("filter_utils.jl")
include("cdf97.jl")
include("boundary.jl")
include("atrous.jl")

export
    FilterPair,
    CDF97,
    upsample_filter,
    ThreadingPolicy, Auto, Enabled, Disabled

end # module
