@testitem "clamp_idx symmetric" tags = [:primitives] begin
    using MultiscaleFilterBanks: clamp_idx
    n = 5
    # interior is the identity
    @test [clamp_idx(i, n, Val(:symmetric)) for i in 1:n] == collect(1:n)
    # mirror without repeating the edge sample
    @test clamp_idx(0, n, Val(:symmetric)) == 2
    @test clamp_idx(-1, n, Val(:symmetric)) == 3
    @test clamp_idx(n + 1, n, Val(:symmetric)) == n - 1
    @test clamp_idx(n + 2, n, Val(:symmetric)) == n - 2
    # degenerate axis
    @test clamp_idx(7, 1, Val(:symmetric)) == 1
    # period is 2(n-1)
    @test clamp_idx(1 + 2 * (n - 1), n, Val(:symmetric)) == 1
end

@testitem "clamp_idx periodic" tags = [:primitives] begin
    using MultiscaleFilterBanks: clamp_idx
    n = 5
    @test [clamp_idx(i, n, Val(:periodic)) for i in 1:n] == collect(1:n)
    @test clamp_idx(0, n, Val(:periodic)) == n
    @test clamp_idx(-1, n, Val(:periodic)) == n - 1
    @test clamp_idx(n + 1, n, Val(:periodic)) == 1
    @test clamp_idx(2n + 3, n, Val(:periodic)) == 3
end

@testitem "clamp_idx zero" tags = [:primitives] begin
    using MultiscaleFilterBanks: clamp_idx
    # the zero mode passes the index through; the caller tests the range
    @test clamp_idx(-3, 5, Val(:zero)) == -3
    @test clamp_idx(9, 5, Val(:zero)) == 9
end
