push!(LOAD_PATH,"../src/")
push!(LOAD_PATH,"../docs/")

using Documenter, DocStringExtensions, VSPSearch
makedocs(
         sitename = "VSPSearch.jl",
         modules  = [VSPSearch],
         format = Documenter.HTML(),
         pages = ["Introduction" => "index.md",
         "API" => "api.md",
        ]
         )
deploydocs(;
    repo="github.com/RenuSolanki/VSPSearch.jl",
    branch="gh-pages"
)
