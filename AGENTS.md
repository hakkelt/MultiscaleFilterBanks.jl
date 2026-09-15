# AGENTS.md — MultiscaleFilterBanks.jl

This file guides AI agents (and human contributors) through the repository
structure and the conventions that are *specific to this package*. General Julia
best practices live in the skills under `.agents/skills/` and are not repeated here.

## What this package is

The shared core of a three-package family:

```
MultiscaleFilterBanks.jl   ← this package (no non-stdlib deps)
        ├── Contourlets.jl   (2-D CT / NSCT)
        └── Surfacelets.jl   (3-D nonsubsampled + decimated surfacelet transform)
```

Both downstream packages depend on this one; neither depends on the other. Nothing
here may be dimension-specific — if a change only makes sense for matrices or only
for volumes, it belongs downstream.

## Repository Layout

```
MultiscaleFilterBanks/
├── src/
│   ├── MultiscaleFilterBanks.jl  # module — all exports here
│   ├── threading.jl              # ThreadingPolicy / Auto / Enabled / Disabled
│   ├── types.jl                  # FilterPair, _data_eltype, _filter_eltype
│   ├── filter_utils.jl           # upsample_filter
│   ├── cdf97.jl                  # CDF 9/7 LP filter constants
│   ├── boundary.jl               # _clamp_idx (symmetric / periodic / zero)
│   └── atrous.jl                 # _NSP_SYNTH_SCALE, _nsp_atrous_filters
├── test/
│   ├── runtests.jl               # TestItemRunner entry point
│   └── items/
└── scripts/relint.jl
```

## Commit Rules

- **Commit Message Formatting**: "Conventional Commits" — `type(scope): subject`,
  then a blank line and a detailed description. Always append a co-author line
  referencing yourself. Do **not** add a `Claude-Session:` trailer line.
- **Committing Changes**: ONLY commit code upon explicit instruction from the user.

## Invariants

- **Zero non-stdlib dependencies.** Adding a dep here adds it to both downstream
  packages. `LoopVectorization`, `Polyester`, `FFTW` deliberately stay downstream.
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
