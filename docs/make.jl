using Documenter, StochasticDiffEq

cp("./docs/Manifest.toml", "./docs/src/assets/Manifest.toml", force = true)
cp("./docs/Project.toml", "./docs/src/assets/Project.toml", force = true)

# Keep pages.jl separate for the DiffEqDocs.jl build
include("pages.jl")

makedocs(
    sitename = "StochasticDiffEq.jl",
    authors = "Chris Rackauckas et al.",
    clean = true,
    doctest = false,
    modules = [StochasticDiffEq],
    warnonly = [:docs_block, :missing_docs, :eval_block],
    format = Documenter.HTML(
        analytics = "UA-90474609-3",
        assets = ["assets/favicon.ico"],
        canonical = "https://stochasticdiffeq.sciml.ai/stable/",
        size_threshold_ignore = []
    ),
    pages = pages
)

deploydocs(
    repo = "github.com/SciML/StochasticDiffEq.jl";
    push_preview = true
)
