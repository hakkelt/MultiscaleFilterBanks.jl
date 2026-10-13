# The mask-bank engine: a nonsubsampled pyramid followed by nonsubsampled
# directional banks, run as pointwise products on the DFT grid of the input, in
# any dimension.  With X the input spectrum,
#
#     subband = ifft(bank.ana[k] · (pyramid_ana[g] · X)),    lowpass = ifft(coarse_ana · X),
#
# and the synthesis is
#
#     X = Σ_g pyramid_syn[g] · (Σ_k bank.syn[k] · fft(subband)) + coarse_syn · fft(lowpass).
#
# The items of one pyramid mask g are consecutive, so the bandpass spectrum
# pyramid_ana[g] · X is formed once per group, and the synthesis applies
# pyramid_syn[g] once to the group's directional sum.  Both loops run in the
# threaded chunks of `threading.jl`.  A chunk that starts inside a group forms
# that group's bandpass again, so the result does not depend on the chunking.
#
# The bank paths are broadcasts and FFTs only and allocate nothing, which CUDA
# graph capture of a transform built on the engine relies on.  A downstream
# package supplies the plans, the masks and `mask_bank_subband` methods for its
# coefficient container.

"""
    MaskItem(index::NTuple{N, Int}, pyramid::Int, bank::Int, channel::Int)

One output of a [`MaskBankEngine`](@ref): where it sits in the coefficient
container (`index`, passed to [`mask_bank_subband`](@ref)), the pyramid mask that
feeds it, and the bank and channel that split that mask.
"""
struct MaskItem{N}
    index::NTuple{N, Int}
    pyramid::Int
    bank::Int
    channel::Int
end

"""
    MaskBankEngine(pyramid_ana, pyramid_syn, coarse_ana, coarse_syn, banks, items, X, pfwd, pinv)

The precomputed masks, FFT plans and buffers of a nonsubsampled pyramid followed by
nonsubsampled directional banks, all on one DFT grid (the `rfft` half-space for
real data), in any dimension.

- `pyramid_ana[g]`, `pyramid_syn[g]`: the analysis and synthesis masks of pyramid
  bandpass `g`; `coarse_ana`, `coarse_syn` those of the lowpass.
- `banks`: directional banks, anything with `ana` and `syn` vectors of channel
  masks.
- `items`: every output as a [`MaskItem`](@ref), with the items of one pyramid
  mask consecutive.
- `X`: the spectrum buffer; `pfwd` / `pinv`: single-threaded plans from data to
  spectrum and back.

With `X̂` the input spectrum, output `(g, b, k)` is
`ifft(banks[b].ana[k] · pyramid_ana[g] · X̂)` and the lowpass is
`ifft(coarse_ana · X̂)`; the synthesis sums the outputs back with the `syn` masks.
[`mask_bank_forward!`](@ref) and [`mask_bank_inverse!`](@ref) run the engine.
The constructor allocates one FFT buffer per task ([`n_tasks`](@ref)).
"""
struct MaskBankEngine{C <: AbstractArray, B, N, P, Q}
    pyramid_ana::Vector{C}      # pyramid masks, input spectrum → bandpass g
    pyramid_syn::Vector{C}      # and back
    coarse_ana::C               # input spectrum → lowpass
    coarse_syn::C
    banks::Vector{B}            # directional banks (shared between groups where equal)
    items::Vector{MaskItem{N}}  # grouped by pyramid mask
    X::C                        # input spectrum / synthesis accumulator
    buffers::Vector{C}          # per task: FFT buffer (consumed by the inverse FFT)
    groups::Vector{C}           # per task: bandpass spectrum / directional sum of a group
    partials::Vector{C}         # per task ≥ 2: partial synthesis sum
    pfwd::P                     # single-threaded plans between data and spectrum
    pinv::Q
end

function MaskBankEngine(
        pyramid_ana::Vector{C}, pyramid_syn::Vector{C}, coarse_ana::C, coarse_syn::C,
        banks::Vector{B}, items::Vector{MaskItem{N}}, X::C, pfwd, pinv
    ) where {C, B, N}
    issorted(items; by = it -> it.pyramid) ||
        throw(ArgumentError("the items of one pyramid mask must be consecutive"))
    nt = n_tasks(X, length(items))
    return MaskBankEngine(
        pyramid_ana, pyramid_syn, coarse_ana, coarse_syn, banks, items, X,
        C[similar(X) for _ in 1:nt], C[similar(X) for _ in 1:nt], C[similar(X) for _ in 2:nt], pfwd, pinv
    )
end

_n_engine_tasks(e::MaskBankEngine) = length(e.buffers)

"""
    mask_bank_subband(coeffs, index)

The array of the output `index` (a [`MaskItem`](@ref)'s `index`) in the
coefficient container `coeffs`.  A package built on [`MaskBankEngine`](@ref) adds
methods for its containers.
"""
function mask_bank_subband end

"""
    mask_bank_coarse(coeffs)

The lowpass array of the coefficient container `coeffs`: its `coarse` field,
unless a method for the container says otherwise.
"""
mask_bank_coarse(coeffs) = coeffs.coarse

"""
    mask_bank_forward!(coeffs, x, e::MaskBankEngine, threaded::Bool) -> coeffs

Analysis of the dense input `x` into the subbands and lowpass of `coeffs`.  With
`threaded`, the outputs are split between tasks ([`foreach_chunk`](@ref)) unless
an enclosing NestedThreading budget forbids it ([`loop_threaded`](@ref)).
Allocates nothing.
"""
function mask_bank_forward!(coeffs, x::AbstractArray, e::MaskBankEngine, threaded::Bool)
    threaded = loop_threaded(threaded)
    X = e.X
    mul!(X, e.pfwd, x)
    foreach_chunk(length(e.items), _n_engine_tasks(e), threaded) do range, t
        _mask_bank_forward_items!(coeffs, e, range, e.groups[t], e.buffers[t])
    end
    Z, K = e.buffers[1], e.coarse_ana
    @.. Z = K * X
    mul!(mask_bank_coarse(coeffs), e.pinv, Z)
    return coeffs
end

# The outputs `e.items[range]`, with the task's group buffer Y (the bandpass
# spectrum, formed when the group changes) and FFT buffer Z.
#
# The products are `@..` (FastBroadcast): Base's broadcast checks every axis for
# a singleton at run time, which keeps these same-size products from vectorising.
# Device arrays fall back to Base's broadcast.
function _mask_bank_forward_items!(coeffs, e::MaskBankEngine, range, Y, Z)
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
        mul!(mask_bank_subband(coeffs, it.index), e.pinv, Z)
    end
    return coeffs
end

"""
    mask_bank_inverse!(x, coeffs, e::MaskBankEngine, threaded::Bool) -> x

Synthesis of the subbands and lowpass of `coeffs` into the dense output `x`, the
inverse of [`mask_bank_forward!`](@ref) for masks that sum to one.  Threaded like
the analysis, with one partial sum per task ([`chunked_sum!`](@ref)).  Allocates
nothing.
"""
function mask_bank_inverse!(x::AbstractArray, coeffs, e::MaskBankEngine, threaded::Bool)
    threaded = loop_threaded(threaded)
    acc = e.X
    fill!(acc, zero(eltype(acc)))
    chunked_sum!(acc, e.partials, length(e.items), threaded) do a, range, t
        _mask_bank_inverse_items!(a, coeffs, e, range, e.buffers[t], e.groups[t])
    end
    F, K = e.buffers[1], e.coarse_syn
    mul!(F, e.pfwd, mask_bank_coarse(coeffs))
    @.. acc += K * F
    mul!(x, e.pinv, acc)
    return x
end

# acc += the synthesis of the outputs `e.items[range]`, with the task's FFT buffer F
# and group buffer Y.
function _mask_bank_inverse_items!(acc, coeffs, e::MaskBankEngine, range, F, Y)
    g = 0
    for n in range
        it = e.items[n]
        if it.pyramid != g
            g == 0 || _add_masked!(acc, e.pyramid_syn[g], Y)
            fill!(Y, zero(eltype(Y)))
            g = it.pyramid
        end
        mul!(F, e.pfwd, mask_bank_subband(coeffs, it.index))
        H = e.banks[it.bank].syn[it.channel]
        @.. Y += H * F
    end
    g == 0 || _add_masked!(acc, e.pyramid_syn[g], Y)
    return acc
end

_add_masked!(acc, P, Y) = @.. acc += P * Y
