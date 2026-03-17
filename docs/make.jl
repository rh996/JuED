using Documenter
using JuED

makedocs(
    sitename="JuED",
    modules=[JuED],
    checkdocs=:exports,
    format=Documenter.HTML(
        prettyurls=false,
        edit_link=nothing,
        repolink="https://github.com/rh996/JuED",
    ),
    remotes=nothing,
    pages=[
        "Home" => "index.md",
        "API" => "api.md",
        "Internals" => "internals.md",
        "Benchmarks" => "benchmarks.md",
    ],
)
