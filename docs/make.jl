using Documenter
using PositionVelocityTime

makedocs(;
    modules = [PositionVelocityTime],
    authors = "Soeren Schoenbrod, Michael Niestroj, Erik Deinzer",
    sitename = "PositionVelocityTime.jl",
    format = Documenter.HTML(;
        canonical = "https://JuliaGNSS.github.io/PositionVelocityTime.jl",
    ),
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md",
    ],
    warnonly = [:missing_docs],
)

deploydocs(;
    repo = "github.com/JuliaGNSS/PositionVelocityTime.jl",
    devbranch = "master",
)
