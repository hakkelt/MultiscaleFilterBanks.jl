@testitem "_nsp_atrous_filters levels 1..4" tags = [:pyramid] begin
    using MultiscaleFilterBanks: _nsp_atrous_filters, _NSP_SYNTH_SCALE
    for lv in 1:4
        h_j, g_j = _nsp_atrous_filters(CDF97, lv, Float64)
        factor = 2^(lv - 1)
        @test h_j == upsample_filter(CDF97.h, factor)
        @test g_j ≈ _NSP_SYNTH_SCALE .* upsample_filter(CDF97.g, factor)
        @test length(h_j) == factor * (length(CDF97.h) - 1) + 1
        @test length(g_j) == factor * (length(CDF97.g) - 1) + 1
        # à trous upsampling preserves DC gain; the synthesis rescale halves it
        @test sum(h_j) ≈ 1.0 atol = 1.0e-12
        @test sum(g_j) ≈ 1.0 atol = 1.0e-12
    end
end

@testitem "_nsp_atrous_filters precision" tags = [:pyramid] begin
    using MultiscaleFilterBanks: _nsp_atrous_filters
    h32, g32 = _nsp_atrous_filters(CDF97, 2, Float32)
    @test eltype(h32) === Float32
    @test eltype(g32) === Float32
    h64, g64 = _nsp_atrous_filters(CDF97, 2, Float64)
    @test h32 ≈ Float32.(h64)
    @test g32 ≈ Float32.(g64)
end
