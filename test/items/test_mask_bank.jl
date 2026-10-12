@testitem "Mask-bank engine: perfect reconstruction on a 3-D grid" tags = [:pyramid] begin
    using MultiscaleFilterBanks: _MaskItem, _MaskBankEngine, _mask_bank_forward!, _mask_bank_inverse!, _subband
    using FFTW, Random

    # A container of this test's own, wired in through `_subband`.
    struct ToyCoefficients{A}
        subbands::Vector{Vector{A}}
        coarse::A
    end
    MultiscaleFilterBanks._subband(c::ToyCoefficients, i::NTuple{2, Int}) = c.subbands[i[1]][i[2]]

    dims = (8, 10, 6)
    rng = Xoshiro(42)
    # Two pyramid levels and a coarse band whose masks telescope:
    # Σ_g syn_g ana_g + coarse_syn coarse_ana = 1.
    w = [rand(rng) for _ in CartesianIndices(dims)]
    ana = [complex.(w), complex.(1 .- w) .* 0.5]
    syn = [ones(ComplexF64, dims), ones(ComplexF64, dims)]
    coarse_ana = complex.(1 .- w) .* 0.5
    coarse_syn = ones(ComplexF64, dims)
    # Each level split by a two-channel bank that sums to one.
    m = complex.([rand(rng) for _ in CartesianIndices(dims)])
    bank = (ana = [m, 1 .- m], syn = [ones(ComplexF64, dims), ones(ComplexF64, dims)])
    items = [_MaskItem((g, k), g, 1, k) for g in 1:2 for k in 1:2]
    x = complex.(randn(rng, dims...), randn(rng, dims...))
    X = similar(x)
    engine = _MaskBankEngine(ana, syn, coarse_ana, coarse_syn, [bank], items, X, plan_fft(x), plan_ifft(X))
    coeffs = ToyCoefficients([[similar(x) for _ in 1:2] for _ in 1:2], similar(x))
    for threaded in (false, true)
        _mask_bank_forward!(coeffs, x, engine, threaded)
        @test _subband(coeffs, (2, 1)) ≈ ifft(m .* ana[2] .* fft(x))
        @test coeffs.coarse ≈ ifft(coarse_ana .* fft(x))
        y = similar(x)
        _mask_bank_inverse!(y, coeffs, engine, threaded)
        @test y ≈ x
    end
end
