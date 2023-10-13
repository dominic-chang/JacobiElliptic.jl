using Documenter, JacobiElliptic

makedocs(sitename="JacobiElliptic.jl")

makedocs(;
    modules=[JacobiElliptic],
    authors="Dominic <dchang3419@hotmail.com> and contributors",
    repo="https://github.com/dchang10/JacobiElliptic.jl/blob/{commit}{path}#{line}",
    sitename="JacobiElliptic.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://dchang10.github.io/JacobiElliptic.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "api.md",
    ],
)

deploydocs(
    repo = "github.com/dchang10/JacobiElliptic.jl",
    devbranch = "main",
    )