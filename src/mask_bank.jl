# The mask-bank engine: a nonsubsampled pyramid followed by nonsubsampled
# directional banks, run as pointwise products on the DFT grid of the input.
# Nothing here is dimension-specific; Contourlets' NSCT and TICT workspaces are
# built on it.
#
# Every output is a pyramid mask g (from the input spectrum to a bandpass) times a
# channel of the directional bank that splits that bandpass.  With X the input
# spectrum,
#
#     subband = ifft(bank.ana[k] · (pyramid_ana[g] · X)),    lowpass = ifft(coarse_ana · X),
#
# and the synthesis is
#
#     X = Σ_g pyramid_syn[g] · (Σ_k bank.syn[k] · fft(subband)) + coarse_syn · fft(lowpass).
#
# The items of one pyramid mask g are consecutive, so the bandpass spectrum
# pyramid_ana[g] · X is formed once per group and the synthesis applies
# pyramid_syn[g] once to the group's directional sum.  That costs a pass per group
# but saves a complex mask and a complex product per item: measured at 256² (the
# NSCT), the elementwise work is ~35% (analysis) and ~45% (synthesis) cheaper than
# the fused per-item product for 16 items per group, and equal for 2.
#
# The outputs are independent, so both loops run in the threaded chunks of
# `threading.jl` (`_foreach_chunk`, `_chunked_sum!`): per task an FFT buffer and a
# group buffer and, for the synthesis from the second task on, a partial sum.
# Plans stay single-threaded.  A chunk that starts inside a group forms that
# group's bandpass again, so the analysis gives the same bits whatever the
# chunking.
#
# A downstream package supplies the plans, the masks and, by adding methods to
# `_subband`, where each output lives in its coefficient container.

# One output: where it sits in the coefficient set (`index`, see `_subband`), the
# pyramid mask that feeds it and the channel of the bank that splits that mask.
struct _MaskItem{N}
    index::NTuple{N, Int}
    pyramid::Int
    bank::Int
    channel::Int
end

"""
    _MaskBankEngine

The precomputed masks, FFT plans and buffers of a nonsubsampled pyramid followed by
nonsubsampled directional banks, all on one DFT grid (the `rfft` half-space for
real data), in any dimension.  `items` lists every output with its pyramid mask
and bank channel; `_mask_bank_forward!` / `_mask_bank_inverse!` run them.  A bank
is anything with `ana` and `syn` vectors of channel masks.
"""
struct _MaskBankEngine{C <: AbstractArray, B, N, P, Q}
    pyramid_ana::Vector{C}      # pyramid masks, input spectrum → bandpass g
    pyramid_syn::Vector{C}      # and back
    coarse_ana::C               # input spectrum → lowpass
    coarse_syn::C
    banks::Vector{B}            # directional banks (shared between groups where equal)
    items::Vector{_MaskItem{N}} # grouped by pyramid mask
    X::C                        # input spectrum / synthesis accumulator
    buffers::Vector{C}          # per task: FFT buffer (consumed by the inverse FFT)
    groups::Vector{C}           # per task: bandpass spectrum / directional sum of a group
    partials::Vector{C}         # per task ≥ 2: partial synthesis sum
    pfwd::P                     # single-threaded plans between data and spectrum
    pinv::Q
end

# Build an engine for the spectrum buffer `X` and the single-threaded plans
# `pfwd` (data → spectrum) and `pinv` (spectrum → data) between them.
function _MaskBankEngine(
        pyramid_ana::Vector{C}, pyramid_syn::Vector{C}, coarse_ana::C, coarse_syn::C,
        banks::Vector{B}, items::Vector{_MaskItem{N}}, X::C, pfwd, pinv
    ) where {C, B, N}
    issorted(items; by = it -> it.pyramid) ||
        throw(ArgumentError("the items of one pyramid mask must be consecutive"))
    nt = _n_tasks(X, length(items))
    return _MaskBankEngine(
        pyramid_ana, pyramid_syn, coarse_ana, coarse_syn, banks, items, X,
        C[similar(X) for _ in 1:nt], C[similar(X) for _ in 1:nt], C[similar(X) for _ in 2:nt], pfwd, pinv
    )
end

_n_engine_tasks(e::_MaskBankEngine) = length(e.buffers)

"""
    _subband(coeffs, index)

The array of the output `index` (an `_MaskItem`'s `index`) in the coefficient
container `coeffs`.  Downstream packages add methods for their containers.
"""
function _subband end

# The lowpass output of a coefficient container: its `coarse` field unless a
# downstream method says otherwise.
_coarse(coeffs) = coeffs.coarse

# Analysis of the dense input `x` into `coeffs` (subbands and lowpass).
function _mask_bank_forward!(coeffs, x::AbstractArray, e::_MaskBankEngine, threaded::Bool)
    threaded = _loop_threaded(threaded)     # serial under a caller's NestedThreading budget
    X = e.X
    mul!(X, e.pfwd, x)
    _foreach_chunk(length(e.items), _n_engine_tasks(e), threaded) do range, t
        _mask_bank_forward_items!(coeffs, e, range, e.groups[t], e.buffers[t])
    end
    Z, K = e.buffers[1], e.coarse_ana
    @.. Z = K * X
    mul!(_coarse(coeffs), e.pinv, Z)
    return coeffs
end

# The outputs `e.items[range]`, with the task's group buffer Y (the bandpass
# spectrum, formed when the group changes) and FFT buffer Z.
#
# The products are `@..` (FastBroadcast): Base's broadcast checks every axis for
# a singleton at run time, which keeps these same-size products from vectorising
# and makes them ~2× slower (8–13% of a single-threaded NSCT or TICT at 256²).
# Device arrays fall back to Base's broadcast.
function _mask_bank_forward_items!(coeffs, e::_MaskBankEngine, range, Y, Z)
    X = e.X
    g = 0
    for n in range
        it = e.items[n]
        if it.pyramid != g
            g = it.pyramid
            P = e.pyramid_ana[g]
            @.. Y = P * X
        end
        H = e.banks[it.bank].ana[it.channel]
        @.. Z = H * Y
        mul!(_subband(coeffs, it.index), e.pinv, Z)
    end
    return coeffs
end

# Synthesis of `coeffs` into the dense output `x`.
function _mask_bank_inverse!(x::AbstractArray, coeffs, e::_MaskBankEngine, threaded::Bool)
    threaded = _loop_threaded(threaded)
    acc = e.X
    fill!(acc, zero(eltype(acc)))
    _chunked_sum!(acc, e.partials, length(e.items), threaded) do a, range, t
        _mask_bank_inverse_items!(a, coeffs, e, range, e.buffers[t], e.groups[t])
    end
    F, K = e.buffers[1], e.coarse_syn
    mul!(F, e.pfwd, _coarse(coeffs))
    @.. acc += K * F
    mul!(x, e.pinv, acc)
    return x
end

# acc += the synthesis of the outputs `e.items[range]`, with the task's FFT buffer F
# and group buffer Y: a group's directional syntheses are summed in Y and its
# pyramid mask applied once.
function _mask_bank_inverse_items!(acc, coeffs, e::_MaskBankEngine, range, F, Y)
    g = 0
    for n in range
        it = e.items[n]
        if it.pyramid != g
            g == 0 || _add_masked!(acc, e.pyramid_syn[g], Y)
            fill!(Y, zero(eltype(Y)))
            g = it.pyramid
        end
        mul!(F, e.pfwd, _subband(coeffs, it.index))
        H = e.banks[it.bank].syn[it.channel]
        @.. Y += H * F
    end
    g == 0 || _add_masked!(acc, e.pyramid_syn[g], Y)
    return acc
end

_add_masked!(acc, P, Y) = @.. acc += P * Y
