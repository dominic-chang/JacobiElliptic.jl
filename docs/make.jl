using Documenter, FElliptic

makedocs(sitename="FElliptic.jl")

makedocs(;
    modules=[FElliptic],
    authors="Dominic <dchang3419@hotmail.com> and contributors",
    repo="https://github.com/dchang10/FElliptic/blob/{commit}{path}#{line}",
    sitename="FElliptic",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://dchang10.github.io/FElliptic",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "api.md",
    ],
)

deploydocs(repo = "github.com/dchang10/FElliptic.git",)