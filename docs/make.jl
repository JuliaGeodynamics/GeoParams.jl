using Documenter, GeoParams, Plots
push!(LOAD_PATH, "../src/")

@info "Making documentation..."
makedocs(;
    sitename="GeoParams.jl",
    authors="Boris Kaus and contributors",
    modules=[GeoParams],
    format=Documenter.HTML(; prettyurls=get(ENV, "CI", nothing) == "true"), # easier local build
    pages=[
        "Home" => "index.md",
        "User Guide" => Any[
            "GeoUnit" => "man/geounit.md",
            "Nondimensionalization" => "man/nondimensionalize.md",
            "Material Parameters" => "man/materialparameters.md",
            "Density" => "man/density.md",
            "Creep laws" => "man/creeplaws.md",
            "Elasticity" => "man/elasticity.md",
            "Plasticity" => "man/plasticity.md",
            "Heat Capacity" => "man/heatcapacity.md",
            "Conductivity" => "man/conductivity.md",
            "Latent heat" => "man/latentheat.md",
            "Radioactive heat" => "man/radioactiveheating.md",
            "Shear heating" => "man/shearheating.md",
            "Gravity" => "man/gravity.md",
            "Partial Melting" => "man/melting.md",
            "TAS classification" => "man/TASclassification.md",
            "Zircon Ages" => "man/zirconages.md",
            "Phase Diagrams" => "man/phasediagrams.md",
            "Seismic Velocity" => "man/seismicvelocity.md",
        ],
        "Plotting" => "man/plotting.md",
        "List of functions" => "man/listfunctions.md",
        "Contributing" => "man/contributing.md",
    ],
)

deploydocs(; repo="github.com/JuliaGeodynamics/GeoParams.jl.git", devbranch="main")
