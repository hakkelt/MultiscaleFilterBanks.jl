# MultiscaleFilterBanks.jl

Shared building blocks for multiscale, multidirectional filter-bank transforms.

This package holds the pieces that are neither 2-D nor 3-D specific and are used by
both [Contourlets.jl](https://github.com/hakkelt/Contourlets.jl) (the 2-D contourlet
and nonsubsampled contourlet transform) and
[Surfacelets.jl](https://github.com/hakkelt/Surfacelets.jl) (the 3-D nonsubsampled
and decimated surfacelet transforms).

| Export | Description |
|--------|-------------|
| `FilterPair` | 1-D biorthogonal analysis/synthesis low-pass filter pair |
| `CDF97` | The CDF 9/7 filter pair used by the pyramid stages |
| `upsample_filter` | À trous ("with holes") filter upsampling |
| `ThreadingPolicy`, `Auto`, `Enabled`, `Disabled` | Threading policy types |

It also carries the internal helpers both downstream packages share: the element-type
convention (`_data_eltype` / `_filter_eltype`), the boundary index folding
(`_clamp_idx`) and the à trous pyramid filter construction (`_nsp_atrous_filters`).

There are no non-stdlib dependencies.

## Development

```bash
julia --project=test test/runtests.jl          # all tests
julia --project=test test/runtests.jl :filters # filter by tag or item name
```

Formatting and static analysis are mandatory before committing:

```bash
runic src/
JULIA_LOAD_PATH="@relint:@stdlib" julia --startup-file=no scripts/relint.jl src/
```

## Status

This package is not registered yet. Downstream packages resolve it through
`Pkg.develop(path = ...)` (see their `dev/bootstrap.jl`). Before CI can go green on
GitHub, the repository needs `CODECOV_TOKEN` set, and dependents need either
registration or `Pkg.develop(url = ...)` in their CI setup step.
