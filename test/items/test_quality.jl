@testitem "Aqua quality checks" tags = [:quality] begin
    using Aqua
    Aqua.test_all(
        MultiscaleFilterBanks;
        persistent_tasks = false,  # system resource check can fail in constrained envs
    )
end

@testitem "Documentation doctests" tags = [:docs, :quality] begin
    using Documenter
    DocMeta.setdocmeta!(
        MultiscaleFilterBanks, :DocTestSetup,
        :(using MultiscaleFilterBanks, Random); recursive = true
    )
    doctest(MultiscaleFilterBanks; manual = false)
end
