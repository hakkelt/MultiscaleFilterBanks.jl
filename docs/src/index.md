```@meta
CurrentModule = MultiscaleFilterBanks
```

# MultiscaleFilterBanks.jl

Shared building blocks for multiscale, multidirectional filter-bank transforms.

A wavelet transform splits an image by scale only; contourlets, surfacelets and
shearlets split every scale further by direction, so each subband holds edges of
one scale and one orientation. They share a multiscale stage, a pyramid built from
a biorthogonal low-pass filter pair, and the machinery to run it fast on CPUs and
GPUs. This package holds that shared part; the [Theory](@ref theory) page explains
it with figures.

## Downstream packages

| Package | Transforms |
|---------|------------|
| [Contourlets.jl](https://github.com/hakkelt/Contourlets.jl) | 2-D contourlet (CT), nonsubsampled contourlet (NSCT), translation-invariant (TICT) and semi-translation-invariant (STICT) contourlet transforms |
| [Surfacelets.jl](https://github.com/hakkelt/Surfacelets.jl) | 3-D surfacelet transform: nonsubsampled, decimated and semi-translation-invariant forms |
| [Shearlets.jl](https://github.com/hakkelt/Shearlets.jl) | 2-D and 3-D shearlet transforms: band-limited (FFST), nonsubsampled (NSST) and compactly supported |

All three re-export its filter types, so most users never load it directly. Load
it to write a new transform of the same kind, or to use the pyramid filters on
their own.

## What it provides

- **Filters**: the biorthogonal pyramid filter pair [`FilterPair`](@ref), its
  default [`CDF97`](@ref), and à trous upsampling with [`upsample_filter`](@ref).
- **Nonsubsampled pyramid**: the filters of each level
  ([`nsp_atrous_filters`](@ref)) and their frequency responses on the DFT grid
  ([`nsp_responses`](@ref)).
- **Mask-bank engine**: [`MaskBankEngine`](@ref) runs a nonsubsampled pyramid
  followed by directional filter banks as pointwise products on the DFT grid, in
  any dimension and on CPU or GPU arrays.
- **Threading**: the threaded loop over independent subbands
  ([`foreach_chunk`](@ref)) and the resolution of the transforms' `threaded`
  keyword ([`use_threading`](@ref)).
- **Utilities**: element-type rules for real filters on real or complex data
  ([`data_eltype`](@ref)) and boundary extension ([`clamp_idx`](@ref)).

The [Theory](@ref theory) page explains the filters, the pyramid and the engine
with figures. The [API Reference](api.md) lists everything.

The package exports nothing, so `using MultiscaleFilterBanks` cannot clash with
other packages; every name is declared `public`. Load the ones you need explicitly:

```jldoctest
julia> using MultiscaleFilterBanks: CDF97, upsample_filter

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

## Multithreading

A filter bank's subbands are independent, so the transforms built on this package
compute them in parallel across Julia threads (start Julia with `-t auto`). Every
transform takes a `threaded` keyword: `true` to thread, `false` to run serially,
or `nothing` (the default) to let the transform decide. The threaded loop splits
the subbands, not the FFTs: each task runs a contiguous chunk of subbands with its
own buffers. It uses Polyester's low-overhead `@batch` when called from the main
task (on Julia 1.12 and later, with one interactive thread) and `Threads.@threads`
from anywhere else, so a transform called inside your own threaded loop stays
fast.

## GPU support

The downstream transforms run on these GPUArrays backends:

| Backend | Array type | Devices |
|---------|------------|---------|
| [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl) | `CuArray` | NVIDIA GPUs |
| [AMDGPU.jl](https://github.com/JuliaGPU/AMDGPU.jl) | `ROCArray` | AMD GPUs |
| [Metal.jl](https://github.com/JuliaGPU/Metal.jl) | `MtlArray` | Apple silicon GPUs |
| [oneAPI.jl](https://github.com/JuliaGPU/oneAPI.jl) | `oneArray` | Intel GPUs |
| [OpenCL.jl](https://github.com/JuliaGPU/OpenCL.jl) | `CLArray` | any OpenCL device |
| [JLArrays.jl](https://github.com/JuliaGPU/GPUArrays.jl/tree/master/lib/JLArrays) | `JLArray` | CPU reference backend, for testing |

The test suites run on CUDA, OpenCL and JLArrays. Move a transform's workspace to
the device with Adapt.jl and pass device arrays, for example with Shearlets.jl:

```julia
using Shearlets, CUDA, Adapt

x = randn(Float32, 256, 256)
sys = BandlimitedShearlets(scales = 4)
ws = Adapt.adapt(CuArray, make_shearlet_workspace(Float32, size(x), sys))
c = shearlet_forward(CuArray(x), sys; workspace = ws)
```
