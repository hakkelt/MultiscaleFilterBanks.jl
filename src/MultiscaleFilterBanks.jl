"""
    MultiscaleFilterBanks

Shared building blocks for multiscale, multidirectional filter-bank transforms.

This package holds the pieces that are neither 2-D nor 3-D specific and are used
by [Contourlets.jl](https://github.com/hakkelt/Contourlets.jl) (the 2-D
contourlet / nonsubsampled contourlet transform),
[Surfacelets.jl](https://github.com/hakkelt/Surfacelets.jl) (the 3-D
nonsubsampled and decimated surfacelet transforms) and
[Shearlets.jl](https://github.com/hakkelt/Shearlets.jl) (2-D and 3-D shearlets):

| Export | Description |
|--------|-------------|
| `FilterPair` | 1-D biorthogonal analysis/synthesis low-pass filter pair |
| `CDF97` | The CDF 9/7 filter pair used by the pyramid stages |
| `upsample_filter` | À trous ("with holes") filter upsampling |
| `ThreadingPolicy`, `Auto`, `Enabled`, `Disabled` | Threading policy types |

It also carries the internal helpers both downstream packages share: the
element-type convention (`_data_eltype` / `_filter_eltype`), the threaded loop
over independent FFTs (`_foreach_chunk` / `_chunked_sum!`), the scoped FFTW thread
count (`_with_fftw_threads`), the boundary index folding (`_clamp_idx`) and the
à trous pyramid filter construction (`_nsp_atrous_filters`) and its frequency
responses (`_nsp_responses`), and the mask-bank engine that runs a nonsubsampled
pyramid plus directional banks as products on the DFT grid (`_MaskBankEngine`).

Its only non-stdlib dependencies are Polyester (the task runner), NestedThreading
(the FFTW thread budget) and FastBroadcast (the engine's same-size products).
"""
module MultiscaleFilterBanks

using FastBroadcast: @..
using LinearAlgebra: mul!
using Polyester: @batch
using NestedThreading: NestedThreading

include("threading.jl")
include("types.jl")
include("filter_utils.jl")
include("cdf97.jl")
include("boundary.jl")
include("atrous.jl")
include("mask_bank.jl")

function __init__()
    _LOOP_THREADS[] = Threads.nthreads()
    _register_loop_pool()
    return nothing
end

export
    FilterPair,
    CDF97,
    upsample_filter,
    ThreadingPolicy, Auto, Enabled, Disabled

end # module
