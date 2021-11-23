using PSplines
using Documenter

DocMeta.setdocmeta!(PSplines, :DocTestSetup, :(using PSplines); recursive=true)

makedocs(;
    modules=[PSplines],
    authors="Víctor Amores",
    repo="https://github.com/victorjesusamoresmedianero/PSplines.jl/blob/{commit}{path}#{line}",
    sitename="PSplines.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://victorjesusamoresmedianero.github.io/PSplines.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/victorjesusamoresmedianero/PSplines.jl",
)
