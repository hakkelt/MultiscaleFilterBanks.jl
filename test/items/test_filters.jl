@testitem "CDF97 coefficients" tags = [:filters] begin
    # DC gain of analysis LP should be ≈ 1
    @test abs(sum(CDF97.h) - 1.0) < 1.0e-12
    # Synthesis LP sums to ≈ 2 (compensates for factor-2 upsampling in LP stage)
    @test abs(sum(CDF97.g) - 2.0) < 1.0e-12
end

@testitem "FilterPair eltype" tags = [:filters] begin
    @test eltype(CDF97) == Float64
    @test eltype(FilterPair([1.0f0], [1.0f0])) == Float32
end

@testitem "FilterPair type promotion" tags = [:filters] begin
    fp = FilterPair([1.0f0, 2.0f0], [1.0, 2.0])   # Float32 vs Float64
    @test eltype(fp) == Float64
    @test fp.h == [1.0, 2.0]
    @test fp.g == [1.0, 2.0]
end

@testitem "upsample_filter" tags = [:filters] begin
    h = [1.0, 2.0, 3.0]
    @test upsample_filter(h, 1) == h
    hu = upsample_filter(h, 2)
    @test length(hu) == 5
    @test hu == [1.0, 0.0, 2.0, 0.0, 3.0]
    @test upsample_filter([1.0f0, 2.0f0], 3) == Float32[1, 0, 0, 2]
end

@testitem "upsample_filter invalid factor" tags = [:filters] begin
    @test_throws ArgumentError upsample_filter([1.0, 2.0], 0)
end

@testitem "element-type convention" tags = [:filters] begin
    using MultiscaleFilterBanks: _data_eltype, _filter_eltype
    @test _data_eltype(Int) === Float64
    @test _data_eltype(Float32) === Float32
    @test _data_eltype(ComplexF32) === ComplexF32
    @test _data_eltype(zeros(Int, 2, 2)) === Float64
    @test _data_eltype(zeros(ComplexF32, 2)) === ComplexF32

    @test _filter_eltype(Float64) === Float64
    @test _filter_eltype(Float32) === Float32
    @test _filter_eltype(ComplexF64) === Float64
    @test _filter_eltype(ComplexF32) === Float32
    @test _filter_eltype(Int) === Float64
end
