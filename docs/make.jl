using Documenter, GeoParams

@info "Making documentation..."
makedocs(
  sitename="GeoParams.jl",
  authors = "Boris Kaus and contributors",
  modules = [GeoParams],
  format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"), # easier local build
  pages = [
        "Home" => "index.md",
        "User Guide" => Any[
        "Nondimensionalization" =>  "man/nondimensionalize.md",
        ],
        "Contributing" => "man/contributing.md",
    ],
)

deploydocs(
    repo = "github.com/JuliaGeodynamics/GeoParams.jl.git",
)
