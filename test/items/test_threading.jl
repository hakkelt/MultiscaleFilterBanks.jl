@testitem "_use_threading truth table" tags = [:primitives] begin
    using MultiscaleFilterBanks: _use_threading
    for T in (Float32, Float64, ComplexF32, ComplexF64, Int)
        @test _use_threading(Enabled(), T) === true
        @test _use_threading(Disabled(), T) === false
    end
    # Auto: off for real (SIMD path), on for complex (multi-core path)
    @test _use_threading(Auto(), Float32) === false
    @test _use_threading(Auto(), Float64) === false
    @test _use_threading(Auto(), Int) === false
    @test _use_threading(Auto(), ComplexF32) === true
    @test _use_threading(Auto(), ComplexF64) === true
end

@testitem "ThreadingPolicy hierarchy" tags = [:primitives] begin
    @test Auto <: ThreadingPolicy
    @test Enabled <: ThreadingPolicy
    @test Disabled <: ThreadingPolicy
end
