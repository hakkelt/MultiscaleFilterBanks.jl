@testitem "@test_call MultiscaleFilterBanks public API" tags = [:jet] begin
    using MultiscaleFilterBanks: FilterPair, CDF97, upsample_filter
    using JET
    using MultiscaleFilterBanks

    # ── Types / constructors ─────────────────────────────────────────────────
    @test_call target_modules = (MultiscaleFilterBanks,) FilterPair([1.0, 2.0], [1.0, 2.0])
    @test_call target_modules = (MultiscaleFilterBanks,) FilterPair([1.0f0, 2.0f0], [1, 2])
    @test_call target_modules = (MultiscaleFilterBanks,) eltype(CDF97)

    # ── Filters ──────────────────────────────────────────────────────────────
    @test_call target_modules = (MultiscaleFilterBanks,) upsample_filter([1.0, 2.0, 3.0], 2)
    @test_call target_modules = (MultiscaleFilterBanks,) upsample_filter(CDF97.h, 4)
    @test_call target_modules = (MultiscaleFilterBanks,) upsample_filter(Float32[1, 2, 3], 1)
end

@testitem "@test_opt MultiscaleFilterBanks public API" tags = [:jet] begin
    using MultiscaleFilterBanks: FilterPair, CDF97, upsample_filter
    using JET
    using MultiscaleFilterBanks

    h = [1.0, 2.0, 3.0]
    @test_opt target_modules = (MultiscaleFilterBanks,) FilterPair(h, h)
    @test_opt target_modules = (MultiscaleFilterBanks,) eltype(CDF97)
    @test_opt target_modules = (MultiscaleFilterBanks,) upsample_filter(h, 2)
    @test_opt target_modules = (MultiscaleFilterBanks,) upsample_filter(CDF97.h, 4)
end

@testitem "@test_call/@test_opt MultiscaleFilterBanks shared internals" tags = [:jet] begin
    using MultiscaleFilterBanks: CDF97
    using JET
    using MultiscaleFilterBanks
    using MultiscaleFilterBanks: data_eltype, filter_eltype, nsp_atrous_filters, nsp_responses,
        clamp_idx, foreach_chunk, with_fftw_threads, use_threading

    x = randn(8, 8)
    for threaded in (nothing, true, false), T in (Float64, ComplexF32)
        @test_call target_modules = (MultiscaleFilterBanks,) use_threading(threaded, T)
        @test_opt target_modules = (MultiscaleFilterBanks,) use_threading(threaded, T)
    end
    @test_call target_modules = (MultiscaleFilterBanks,) data_eltype(x)
    @test_call target_modules = (MultiscaleFilterBanks,) filter_eltype(ComplexF32)
    @test_opt target_modules = (MultiscaleFilterBanks,) data_eltype(x)
    @test_opt target_modules = (MultiscaleFilterBanks,) filter_eltype(ComplexF32)

    @test_call target_modules = (MultiscaleFilterBanks,) nsp_atrous_filters(CDF97, 2, Float64)
    @test_opt target_modules = (MultiscaleFilterBanks,) nsp_atrous_filters(CDF97, 2, Float64)
    @test_call target_modules = (MultiscaleFilterBanks,) nsp_responses(CDF97, 2, 16)

    out = zeros(Int, 4)
    f = (range, _) -> (
        for i in range
            out[i] += 1
        end
    )
    @test_call target_modules = (MultiscaleFilterBanks,) foreach_chunk(f, 4, 2, false)
    @test_call target_modules = (MultiscaleFilterBanks,) foreach_chunk(f, 4, 2, true)
    @test_call target_modules = (MultiscaleFilterBanks,) clamp_idx(0, 8, Val(:symmetric))
    @test_call target_modules = (MultiscaleFilterBanks,) clamp_idx(9, 8, Val(:periodic))
    @test_opt target_modules = (MultiscaleFilterBanks,) clamp_idx(0, 8, Val(:symmetric))
    @test_call target_modules = (MultiscaleFilterBanks,) with_fftw_threads(() -> nothing, 1)
end
