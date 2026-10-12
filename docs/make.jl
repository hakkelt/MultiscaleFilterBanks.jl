using Documenter, MultiscaleFilterBanks
using Documenter: Remotes

DocMeta.setdocmeta!(
    MultiscaleFilterBanks, :DocTestSetup,
    :(using MultiscaleFilterBanks, Random); recursive = true
)

makedocs(;
    modules = [MultiscaleFilterBanks],
    sitename = "MultiscaleFilterBanks.jl",
    # The repository has no `origin` remote in a fresh local clone, and
    # Documenter derives source links from it, so name it explicitly.
    repo = Remotes.GitHub("hakkelt", "MultiscaleFilterBanks.jl"),
    authors = "Tamás Hakkel",
    format = Documenter.HTML(;
        canonical = "https://hakkelt.github.io/MultiscaleFilterBanks.jl",
    ),
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md",
    ],
    checkdocs = :exports,
)

deploydocs(;
    repo = "github.com/hakkelt/MultiscaleFilterBanks.jl",
    devbranch = "master",
)
