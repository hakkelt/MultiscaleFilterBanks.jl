# MultiscaleFilterBanks.jl

[![CI](https://github.com/hakkelt/MultiscaleFilterBanks.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/hakkelt/MultiscaleFilterBanks.jl/actions/workflows/CI.yml)
[![Lint](https://github.com/hakkelt/MultiscaleFilterBanks.jl/actions/workflows/lint.yml/badge.svg?branch=master)](https://github.com/hakkelt/MultiscaleFilterBanks.jl/actions/workflows/lint.yml)
[![Documentation (stable)](https://img.shields.io/badge/docs-stable-blue.svg)](https://hakkelt.github.io/MultiscaleFilterBanks.jl/stable/)
[![Documentation (dev)](https://img.shields.io/badge/docs-dev-blue.svg)](https://hakkelt.github.io/MultiscaleFilterBanks.jl/dev/)
[![Coverage](https://codecov.io/gh/hakkelt/MultiscaleFilterBanks.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/hakkelt/MultiscaleFilterBanks.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![JET](https://img.shields.io/badge/%F0%9F%9B%A9%EF%B8%8F_tested_with-JET.jl-233f9a)](https://github.com/aviatesk/JET.jl)
[![Runic](https://img.shields.io/badge/code%20style-runic-blue.svg)](https://github.com/fredrikekre/Runic.jl)
[![License: MIT](https://img.shields.io/badge/license-MIT-yellow.svg)](LICENSE)

The shared core of a family of multiscale, multidirectional transforms for images
and volumes:

| Package | Transforms |
|---------|------------|
| [Contourlets.jl](https://github.com/hakkelt/Contourlets.jl) | 2-D contourlet (CT), nonsubsampled contourlet (NSCT), translation-invariant (TICT) and semi-translation-invariant (STICT) contourlet transforms |
| [Surfacelets.jl](https://github.com/hakkelt/Surfacelets.jl) | 3-D surfacelet transform: nonsubsampled, decimated and semi-translation-invariant forms |
| [Shearlets.jl](https://github.com/hakkelt/Shearlets.jl) | 2-D and 3-D shearlet transforms: band-limited (FFST), nonsubsampled (NSST) and compactly supported |

Most users never load this package directly: all three re-export its filter
types. It is for writing a new transform of the same kind, or for using the
pyramid filters on their own.

## What is a multiscale filter bank?

A wavelet transform splits a signal by *scale*: each level separates the details
of one frequency octave from a coarser approximation. In two and more dimensions
the separable wavelet sees only horizontal, vertical and diagonal detail, so a
curved edge costs many coefficients at every scale. Multiscale, multidirectional
transforms fix that by splitting each octave further by *orientation*:

1. A **pyramid** (multiscale stage) splits the input into a bandpass and a
   coarser lowpass, and repeats on the lowpass. The analysis low-pass `h` produces
   the lowpass; the synthesis low-pass `g` predicts the input from it, and the
   bandpass is what the prediction misses (the Laplacian pyramid of Burt and
   Adelson). With the *à trous* variant nothing is decimated, so the transform is
   invariant to shifts.
2. A **directional filter bank** (multidirectional stage) splits each bandpass
   into wedges of orientation in the frequency plane (contourlets), in frequency
   space (surfacelets), or into sheared copies of one generator (shearlets).

Every subband then holds edges of one scale and one direction, which makes these
transforms sparser than wavelets for images and volumes with edges and surfaces,
and useful for denoising, enhancement, compression and compressed-sensing
reconstruction. This package provides the parts the three transforms share; the
[theory page](https://hakkelt.github.io/MultiscaleFilterBanks.jl/stable/theory/)
of the documentation explains them with figures.

## Installation

```julia
using Pkg
Pkg.add("MultiscaleFilterBanks")
```

## Quick example

```julia
using MultiscaleFilterBanks: FilterPair, CDF97, upsample_filter

CDF97                                  # FilterPair{Float64}: 9-tap h, 7-tap g
sum(CDF97.h), sum(CDF97.g)             # ≈ (1, 2): DC gains of analysis and synthesis
upsample_filter(CDF97.h, 4)            # level-3 à trous filter: 3 zeros between taps

fp = FilterPair([0.5, 1.0, 0.5] ./ 2, [0.5, 1.0, 0.5])   # a custom pair (linear spline)
```

The downstream packages take the pair as a parameter, e.g.
`ContourletParams(...; lp_filters = fp)` in Contourlets.jl or
`NSSTShearlets(; pyramid = fp)` in Shearlets.jl.

## What it provides

- **Filters**: the biorthogonal pyramid filter pair `FilterPair`, its default
  `CDF97`, and à trous upsampling with `upsample_filter`.
- **Nonsubsampled pyramid**: the filters of each level and their frequency
  responses on the DFT grid.
- **Mask-bank engine**: `MaskBankEngine` runs a nonsubsampled pyramid followed by
  directional filter banks as pointwise products on the DFT grid, in any
  dimension and on CPU or GPU arrays.
- **Threading**: the threaded loop over independent subbands and the resolution
  of the transforms' `threaded` keyword.
- **Utilities**: element-type rules for real filters on real or complex data, and
  boundary extension.

The [theory page](https://hakkelt.github.io/MultiscaleFilterBanks.jl/stable/theory/)
explains the filters, the pyramid and the engine with figures, and the
[API reference](https://hakkelt.github.io/MultiscaleFilterBanks.jl/stable/api/)
lists everything. The package exports nothing, so `using MultiscaleFilterBanks`
cannot clash with other packages; every name is declared `public`. Load the ones
you need explicitly, as in the example above.

## Multithreading

The transforms built on this package compute their subbands independently, so
they run them in parallel across Julia threads (start Julia with `-t auto`). Each
transform takes a `threaded` keyword: `true` to thread, `false` to run serially,
or `nothing` (the default) to let the transform decide.

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

## Element types

Filters are real (`Float32` or `Float64`); data may be real or complex. Complex
data is filtered with real kernels, so every transform is exactly real-linear:
`T(x + im*y) == T(x) + im*T(y)`.

## Development

```bash
julia --project=test test/runtests.jl          # all tests
julia --project=test test/runtests.jl :filters # filter by tag or item name
julia --project=docs docs/make.jl              # build the documentation
```

Formatting and static analysis are mandatory before committing:

```bash
runic src/
JULIA_LOAD_PATH="@relint:@stdlib" julia --startup-file=no scripts/relint.jl src/
```

## Disclaimer

This package was created by [Claude](https://claude.com/claude-code), Anthropic's
AI model, under the direction of the author.

## License

MIT; see [`LICENSE`](LICENSE).
