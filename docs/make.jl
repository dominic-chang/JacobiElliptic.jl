using Documenter, JacobiElliptic
using Documenter: Remotes

makedocs(;
    modules = [JacobiElliptic],
    authors = "Dominic <dochang@g.harvard.edu.com> and contributors",
    repo = Remotes.GitHub("dominic-chang", "JacobiElliptic.jl"),
    sitename = "JacobiElliptic.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://dominic-chang.github.io/JacobiElliptic.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = ["Home" => "index.md", "api.md"],
)

deploydocs(repo = "github.com/dominic-chang/JacobiElliptic.jl", devbranch = "main")
