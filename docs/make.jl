using Documenter, FastElliptic

makedocs(sitename="FastElliptic.jl")

makedocs(;
    modules=[FastElliptic],
    authors="Dominic <dchang3419@hotmail.com> and contributors",
    repo="https://github.com/dchang10/FastElliptic/blob/{commit}{path}#{line}",
    sitename="FastElliptic",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://dchang10.github.io/FastElliptic",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "api.md",
    ],
)

deploydocs(repo = "github.com/dchang10/FastElliptic.git",)