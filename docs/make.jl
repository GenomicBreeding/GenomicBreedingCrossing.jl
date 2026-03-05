using GenomicBreedingCrossing
using Documenter

DocMeta.setdocmeta!(
    GenomicBreedingCrossing,
    :DocTestSetup,
    :(using GenomicBreedingCrossing);
    recursive = true,
)

makedocs(;
    modules = [GenomicBreedingCrossing],
    authors = "jeffersonparil@gmail.com",
    sitename = "GenomicBreedingCrossing.jl",
    format = Documenter.HTML(;
        canonical = "https://GenomicBreeding.github.io/GenomicBreedingCrossing.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)

deploydocs(;
    repo = "github.com/GenomicBreeding/GenomicBreedingCrossing.jl",
    devbranch = "main",
)
