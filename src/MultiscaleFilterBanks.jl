"""
    MultiscaleFilterBanks

Shared core of the multiscale, multidirectional filter-bank transforms
[Contourlets.jl](https://github.com/hakkelt/Contourlets.jl) (2-D contourlets),
[Surfacelets.jl](https://github.com/hakkelt/Surfacelets.jl) (3-D surfacelets) and
[Shearlets.jl](https://github.com/hakkelt/Shearlets.jl) (2-D and 3-D shearlets):
everything that is neither 2-D nor 3-D specific.

The package exports nothing; its API is declared `public`.  Load names
explicitly, e.g. `using MultiscaleFilterBanks: CDF97, upsample_filter`.  The
downstream packages re-export the filter names, so their users do not need to.

- Filters: [`FilterPair`](@ref), [`CDF97`](@ref), [`upsample_filter`](@ref).
- Nonsubsampled pyramid: [`nsp_atrous_filters`](@ref), [`nsp_responses`](@ref),
  [`centred_dtft`](@ref).
- Mask-bank engine: [`MaskBankEngine`](@ref), [`MaskItem`](@ref),
  [`mask_bank_forward!`](@ref), [`mask_bank_inverse!`](@ref),
  [`mask_bank_subband`](@ref), [`mask_bank_coarse`](@ref).
- Threading: [`use_threading`](@ref), [`loop_threaded`](@ref),
  [`foreach_chunk`](@ref), [`chunked_sum!`](@ref), [`add_partials!`](@ref),
  [`n_tasks`](@ref), [`with_fftw_threads`](@ref).
- Element types and boundaries: [`data_eltype`](@ref), [`filter_eltype`](@ref),
  [`clamp_idx`](@ref).
"""
module MultiscaleFilterBanks

using Compat: @compat
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

@compat public FilterPair, CDF97, upsample_filter,
    nsp_atrous_filters, nsp_responses, centred_dtft,
    MaskBankEngine, MaskItem, mask_bank_forward!, mask_bank_inverse!, mask_bank_subband, mask_bank_coarse,
    use_threading, loop_threaded, foreach_chunk, chunked_sum!, add_partials!, n_tasks, with_fftw_threads,
    data_eltype, filter_eltype, clamp_idx

end # module
