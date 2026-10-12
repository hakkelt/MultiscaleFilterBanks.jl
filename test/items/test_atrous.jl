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

@testitem "_nsp_responses: DTFTs of the à trous pairs, telescoping to 1" tags = [:pyramid] begin
    using MultiscaleFilterBanks: _nsp_responses, _nsp_atrous_filters, _centred_dtft
    n, J = 24, 3
    a, gh, s = _nsp_responses(CDF97, J, n)
    @test length(a) == length(s) == J + 1
    @test length(gh) == J
    # The response of a centred periodic convolution: compare with a DFT of the
    # filter wrapped onto the grid with its centre tap at index 1.
    h, g = _nsp_atrous_filters(CDF97, 2, Float64)
    c = (length(h) + 1) ÷ 2
    w = zeros(n)
    for l in eachindex(h)
        w[mod(l - c, n) + 1] += h[l]
    end
    dft = [sum(w[m] * cis(-2π * (k - 1) * (m - 1) / n) for m in 1:n) for k in 1:n]
    @test maximum(abs, a[3] ./ a[2] .- dft) < 1.0e-12
    # Σⱼ sⱼ aⱼ (1 − ghⱼ) + s_{J+1} a_{J+1} = 1, in D = 3 dimensions too.
    total = sum(s[j] .* a[j] .* (1 .- gh[j]) for j in 1:J) .+ s[J + 1] .* a[J + 1]
    @test maximum(abs, total .- 1) < 1.0e-12
    sep(v) = v .* transpose(v) .* reshape(v, 1, 1, :)
    total3 = sum(sep(s[j]) .* sep(a[j]) .* (1 .- sep(gh[j])) for j in 1:J) .+ sep(s[J + 1]) .* sep(a[J + 1])
    @test maximum(abs, total3 .- 1) < 1.0e-12
    # A subset of bins is the same responses, restricted.
    ar, _, _ = _nsp_responses(CDF97, J, n; bins = 1:(n ÷ 2 + 1))
    @test ar[J + 1] == a[J + 1][1:(n ÷ 2 + 1)]
    @test _centred_dtft([0.0, 1.0, 0.0], 0.7) ≈ 1.0
end
