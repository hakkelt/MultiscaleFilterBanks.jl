@testitem "@test_call MultiscaleFilterBanks public API" tags = [:jet] begin
    using JET
    using MultiscaleFilterBanks

    # ── Types / constructors ─────────────────────────────────────────────────
    @test_call target_modules = (MultiscaleFilterBanks,) FilterPair([1.0, 2.0], [1.0, 2.0])
    @test_call target_modules = (MultiscaleFilterBanks,) FilterPair([1.0f0, 2.0f0], [1, 2])
    @test_call target_modules = (MultiscaleFilterBanks,) eltype(CDF97)
    @test_call target_modules = (MultiscaleFilterBanks,) Auto()
    @test_call target_modules = (MultiscaleFilterBanks,) Enabled()
    @test_call target_modules = (MultiscaleFilterBanks,) Disabled()

    # ── Filters ──────────────────────────────────────────────────────────────
    @test_call target_modules = (MultiscaleFilterBanks,) upsample_filter([1.0, 2.0, 3.0], 2)
    @test_call target_modules = (MultiscaleFilterBanks,) upsample_filter(CDF97.h, 4)
    @test_call target_modules = (MultiscaleFilterBanks,) upsample_filter(Float32[1, 2, 3], 1)
end

@testitem "@test_opt MultiscaleFilterBanks public API" tags = [:jet] begin
    using JET
    using MultiscaleFilterBanks

    h = [1.0, 2.0, 3.0]
    @test_opt target_modules = (MultiscaleFilterBanks,) FilterPair(h, h)
    @test_opt target_modules = (MultiscaleFilterBanks,) eltype(CDF97)
    @test_opt target_modules = (MultiscaleFilterBanks,) upsample_filter(h, 2)
    @test_opt target_modules = (MultiscaleFilterBanks,) upsample_filter(CDF97.h, 4)
    @test_opt target_modules = (MultiscaleFilterBanks,) Auto()
    @test_opt target_modules = (MultiscaleFilterBanks,) Enabled()
    @test_opt target_modules = (MultiscaleFilterBanks,) Disabled()
end

@testitem "@test_call/@test_opt MultiscaleFilterBanks shared internals" tags = [:jet] begin
    using JET
    using MultiscaleFilterBanks
    using MultiscaleFilterBanks: _data_eltype, _filter_eltype, _nsp_atrous_filters, _nsp_responses,
        _clamp_idx, _foreach_chunk, _with_fftw_threads

    x = randn(8, 8)
    @test_call target_modules = (MultiscaleFilterBanks,) _data_eltype(x)
    @test_call target_modules = (MultiscaleFilterBanks,) _filter_eltype(ComplexF32)
    @test_opt target_modules = (MultiscaleFilterBanks,) _data_eltype(x)
    @test_opt target_modules = (MultiscaleFilterBanks,) _filter_eltype(ComplexF32)

    @test_call target_modules = (MultiscaleFilterBanks,) _nsp_atrous_filters(CDF97, 2, Float64)
    @test_opt target_modules = (MultiscaleFilterBanks,) _nsp_atrous_filters(CDF97, 2, Float64)
    @test_call target_modules = (MultiscaleFilterBanks,) _nsp_responses(CDF97, 2, 16)

    out = zeros(Int, 4)
    f = (range, _) -> (
        for i in range
            out[i] += 1
        end
    )
    @test_call target_modules = (MultiscaleFilterBanks,) _foreach_chunk(f, 4, 2, false)
    @test_call target_modules = (MultiscaleFilterBanks,) _foreach_chunk(f, 4, 2, true)
    @test_call target_modules = (MultiscaleFilterBanks,) _clamp_idx(0, 8, Val(:symmetric))
    @test_call target_modules = (MultiscaleFilterBanks,) _clamp_idx(9, 8, Val(:periodic))
    @test_opt target_modules = (MultiscaleFilterBanks,) _clamp_idx(0, 8, Val(:symmetric))
    @test_call target_modules = (MultiscaleFilterBanks,) _with_fftw_threads(() -> nothing, 1)
end
