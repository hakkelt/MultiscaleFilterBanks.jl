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

All three re-export its filter types, so most users never load it directly. Load it to write a new transform of the same kind, or to use the
pyramid filters on their own.

## What it provides

- [`FilterPair`](@ref): a 1-D biorthogonal analysis/synthesis low-pass filter
  pair, and [`CDF97`](@ref), the CDF 9/7 pair of JPEG 2000 that the pyramids use
  by default.
- [`upsample_filter`](@ref): the à trous filter upsampling of the nonsubsampled
  (shift-invariant) pyramids and directional banks.
- [`MaskBankEngine`](@ref): the engine that runs a nonsubsampled pyramid plus
  directional banks as pointwise products on the DFT grid, in any dimension.
- The threaded loop over independent FFTs ([`foreach_chunk`](@ref)) and the
  helpers for element types and boundary extension.

The package exports nothing, so `using MultiscaleFilterBanks` cannot clash with
other packages; every name in the [API Reference](api.md) is declared `public`.
Load the ones you need explicitly:

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
or `nothing` (the default) to let the transform decide. The threaded loop splits the subbands, not the FFTs: each task runs a contiguous
chunk of subbands with its own buffers. It uses Polyester's low-overhead `@batch`
when called from the main task (on Julia 1.12 and later, with one interactive
thread) and `Threads.@threads` from anywhere else, so a transform called inside
your own threaded loop stays fast.

The loop's thread budget is registered with
[NestedThreading.jl](https://github.com/hakkelt/NestedThreading.jl). Inside a
`@budgeted_threads` loop or a `with_restricted_threads` scope the transforms run
serially instead of oversubscribing the machine, and FFTW plans made inside such a
scope are single-threaded.

## GPU support

The nonsubsampled transforms run as products of precomputed masks on the DFT
grid: broadcasts and FFT plans that the caller supplies. Nothing in that path is
specific to `Array`, so the downstream transforms accept GPU arrays such as
`CuArray` and run on the device, with one task per transform since the device FFTs
are already parallel. The path allocates nothing per call, which lets
Contourlets.jl capture its NSCT in a CUDA graph.

## Dependencies

Besides the standard library, the package depends on
[Polyester.jl](https://github.com/JuliaSIMD/Polyester.jl) (the task runner),
[NestedThreading.jl](https://github.com/hakkelt/NestedThreading.jl) (the thread
budget) and [FastBroadcast.jl](https://github.com/YingboMa/FastBroadcast.jl) (the
engine's same-size products). FFTW is deliberately not a dependency: the engine
takes its FFT plans from the caller.

## References

See the [Theory](@ref theory) page and the docstrings of [`FilterPair`](@ref),
[`CDF97`](@ref) and [`upsample_filter`](@ref).
