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

@testitem "_foreach_chunk / _chunked_sum! cover 1:n once" tags = [:primitives] begin
    using MultiscaleFilterBanks: _foreach_chunk, _chunked_sum!, _chunk, _n_tasks
    for n in (1, 5, 17), nt in (1, 2, 3, 4), threaded in (false, true)
        seen = zeros(Int, n)
        tasks = zeros(Int, n)
        _foreach_chunk(n, nt, threaded) do range, t
            seen[range] .+= 1
            tasks[range] .= t
        end
        @test all(==(1), seen)
        @test issorted(tasks)
        acc = fill(1.0, 3)
        partials = [zeros(3) for _ in 2:nt]
        _chunked_sum!(acc, partials, n, threaded) do a, range, t
            for i in range
                a .+= i
            end
        end
        # Unthreaded, task 1 sums everything into `acc` itself.
        @test acc == fill(1.0 + n * (n + 1) ÷ 2, 3)
    end
    @test _n_tasks(zeros(2), 1) == 1
    @test _n_tasks(zeros(2), 100) == Threads.nthreads()
end
