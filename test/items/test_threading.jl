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

@testitem "_foreach_chunk rethrows and stays threaded" tags = [:primitives] begin
    using MultiscaleFilterBanks: _foreach_chunk
    nt = max(2, Threads.nthreads())
    @test_throws DimensionMismatch _foreach_chunk(nt, nt, true) do range, t
        t == 2 && throw(DimensionMismatch("chunk $t"))
    end
    # Polyester's workers survive the error: the next loop still covers 1:n once
    # and, with threads, runs on more than one of them.
    seen = zeros(Int, 4nt)
    ids = zeros(Int, nt)
    work = zeros(nt)
    _foreach_chunk(4nt, nt, true) do range, t
        seen[range] .+= 1
        ids[t] = Threads.threadid()
        for i in 1:2_000_000
            work[t] += sin(i)
        end
    end
    @test all(==(1), seen)
    Threads.nthreads() > 1 && @test length(unique(ids)) > 1
end

@testitem "_foreach_chunk avoids @batch inside a caller's threaded loop" tags = [:primitives] begin
    using MultiscaleFilterBanks: _batch_safe, _foreach_chunk
    # A task of a caller's `@threads`/`@spawn` runs on the default pool.
    @test fetch(Threads.@spawn _batch_safe()) == false
    # Nested use still covers every item once.
    seen = zeros(Int, 4, 16)
    Threads.@threads for b in 1:4
        _foreach_chunk(16, 4, true) do range, t
            seen[b, range] .+= 1
        end
    end
    @test all(==(1), seen)
end

@testitem "NestedThreading budgets turn the loops serial" tags = [:primitives] begin
    using MultiscaleFilterBanks: _loop_threaded
    using NestedThreading
    @test _loop_threaded(true) == (Threads.nthreads() > 1)
    @test _loop_threaded(false) == false
    @test with_restricted_threads(() -> _loop_threaded(true)) == false
    @test with_restricted_threads(() -> _loop_threaded(true); only = (:fftw,)) == (Threads.nthreads() > 1)
    inside = zeros(Bool, 4)
    @budgeted_threads for b in 1:4
        inside[b] = _loop_threaded(true)
    end
    Threads.nthreads() >= 4 && @test !any(inside)
    @test _loop_threaded(true) == (Threads.nthreads() > 1)
end

@testitem "_add_partials! host and generic paths" tags = [:primitives] begin
    using MultiscaleFilterBanks: _add_partials!
    using Random
    rng = Xoshiro(42)
    for n in (7, 5000, 100_000), np in (0, 1, 3), threaded in (false, true)
        acc = randn(rng, ComplexF64, n)
        parts = [randn(rng, ComplexF64, n) for _ in 1:np]
        expect = foldl(+, parts; init = copy(acc))
        # A non-`Array` accumulator takes the broadcast path.
        acc2 = view(copy(acc), :)
        @test _add_partials!(acc, parts, threaded) ≈ expect
        @test _add_partials!(acc2, parts, threaded) ≈ expect
    end
end

@testitem "_with_fftw_threads scopes the FFTW count" tags = [:primitives] begin
    using MultiscaleFilterBanks: _with_fftw_threads
    using FFTW, NestedThreading
    FFTW.set_num_threads(1)
    @test _with_fftw_threads(() -> FFTW.get_num_threads(), 3) == 3
    @test FFTW.get_num_threads() == 1
    # An enclosing restriction wins over the request.
    @test with_restricted_threads(() -> _with_fftw_threads(() -> FFTW.get_num_threads(), 3)) == 1
    @test FFTW.get_num_threads() == 1
end
