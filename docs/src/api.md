```@meta
CurrentModule = MultiscaleFilterBanks
```

# API Reference

MultiscaleFilterBanks exports nothing; every name below is `public`. Load the ones
you need explicitly, e.g. `using MultiscaleFilterBanks: CDF97, upsample_filter`.

```@docs
MultiscaleFilterBanks
```

## Filters

```@docs
FilterPair
CDF97
upsample_filter
```

## Nonsubsampled pyramid

```@docs
nsp_atrous_filters
nsp_responses
centred_dtft
```

## Mask-bank engine

The engine behind the nonsubsampled transforms of the downstream packages: a
nonsubsampled pyramid followed by nonsubsampled directional banks, run as
pointwise products on the DFT grid in any dimension.

```@docs
MaskBankEngine
MaskItem
mask_bank_forward!
mask_bank_inverse!
mask_bank_subband
mask_bank_coarse
```

## Threading

```@docs
use_threading
loop_threaded
foreach_chunk
chunked_sum!
add_partials!
n_tasks
with_fftw_threads
```

## Element types and boundaries

```@docs
data_eltype
filter_eltype
clamp_idx
```
