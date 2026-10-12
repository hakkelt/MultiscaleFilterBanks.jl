# AGENTS.md — MultiscaleFilterBanks.jl

This file guides AI agents (and human contributors) through the repository
structure and the conventions that are *specific to this package*. General Julia
practice is not repeated here.

## What this package is

The shared core of a package family:

```
MultiscaleFilterBanks.jl   ← this package
        ├── Contourlets.jl   (2-D CT / NSCT)
        ├── Surfacelets.jl   (3-D nonsubsampled + decimated surfacelet transform)
        └── Shearlets.jl     (2-D and 3-D shearlets)
```

The downstream packages depend on this one; none depends on another. Nothing
here may be dimension-specific — if a change only makes sense for matrices or only
for volumes, it belongs downstream.

## Commit Rules

- **Commit Message Formatting**: "Conventional Commits" — `type(scope): subject`,
  then a blank line and a detailed description. Always append a co-author line
  referencing yourself.
- **Committing Changes**: ONLY commit code upon explicit instruction from the user.

## Invariants

- **Dependencies only where the downstream packages need them.** Adding a dep
  here adds it to all of them. `Polyester` (the task runner of `_foreach_chunk`),
  `NestedThreading` (`_with_fftw_threads`) and `FastBroadcast` (the engine's
  same-size products) qualify; `LoopVectorization` and `FFTW` deliberately stay
  downstream — the engine takes its FFT plans from the caller.
- **The mask-bank engine (`mask_bank.jl`) is dimension-agnostic.** It knows
  nothing of coefficient containers: downstream packages add `_subband` (and,
  if the lowpass is not a `coarse` field, `_coarse`) methods. Its bank paths are
  broadcast + FFT only and allocate nothing — Contourlets' CUDA graph capture of
  the NSCT relies on it.
- **`@batch` only where `_batch_safe()`.** Polyester's workers are pinned to
  threads 2…nthreads(): from a task of a caller's threaded loop it ran ~100×
  slower than `Threads.@threads`, and a busy second interactive thread stalls it.
  `_foreach_chunk` falls back to `Threads.@threads` there; a new threaded loop
  should go through it rather than call `@batch`.
- **No `@batch` body may throw.** An exception escaping one leaves Polyester's
  workers unusable, and every later `@batch` in the process runs serially without
  a warning. `_foreach_chunk` catches and rethrows.
- **Element types: real filters, real-or-complex data.** `FilterPair` is
  `Tf <: AbstractFloat`; data may be complex (`Td <: Number`). Downstream threads
  the two independently — `Td = _data_eltype(x)`, `Tf = _filter_eltype(Td) = real(float(Td))`
  — so a complex image is filtered as complex·real, never complex·complex. The
  transforms built on this are real-linear, so `T(x + i·y) == T(x) + i·T(y)`
  bit-for-bit; downstream test suites assert that with `== 0.0`. Any change to
  `_data_eltype` / `_filter_eltype` breaks that invariant silently.
- **`_clamp_idx` takes the boundary mode as `Val{B}`**, never a runtime `Symbol`;
  the callers rely on full specialisation with no per-sample dispatch.
- **`_nsp_atrous_filters` carries the `0.5` synthesis rescale.** The biorthogonal
  synthesis `g` has DC gain 2 to compensate the ×2 upsampling of the *decimated*
  Laplacian pyramid; the nonsubsampled pyramid does no decimation, so `g` is
  rescaled per axis. Both the 2-D and the 3-D pyramid call this and depend on the
  rescale being applied exactly once. It is the only helper here with real
  behavioural coupling to downstream numerics — changing it changes reconstruction
  in both packages.
- **Formatting.** Runic is mandatory — run `runic src/` before committing, then
  ReLint:
  ```bash
  runic src/
  JULIA_LOAD_PATH="@relint:@stdlib" julia --startup-file=no scripts/relint.jl src/
  ```

## Testing Conventions

- All tests are self-contained `@testitem` blocks (TestItemRunner): each block
  imports everything it needs except `MultiscaleFilterBanks` and `Test`.
- Tags: `:filters`, `:primitives`, `:pyramid`, `:quality`, `:docs`.
- Run: `julia --project=test test/runtests.jl`, optionally with a comma-separated
  tag/name filter (`julia --project=test test/runtests.jl :filters`).
