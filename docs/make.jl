using Documenter
using DocumenterVitepress
using GeoParams, Makie
# push!(LOAD_PATH, "../src/")

@info "Making documentation..."
makedocs(;
    sitename = "GeoParams.jl",
    authors = "Boris Kaus and contributors",
    repo = Remotes.GitHub("JuliaGeodynamics", "GeoParams.jl"),
    modules = [
        GeoParams,
        isdefined(Base, :get_extension) ?
            Base.get_extension(GeoParams, :GeoParamsMakieExt) :
            GeoParams.GeoParamsMakieExt,
    ],
    warnonly = Documenter.except(:footnote),
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "github.com/JuliaGeodynamics/GeoParams.jl",
        devbranch = "main",
        devurl = "dev",
    ),
    pages = [
        "Home" => "index.md",
        "User Guide" => Any[
            "GeoUnit" => "man/geounit.md",
            "Nondimensionalization" => "man/nondimensionalize.md",
            "Material Parameters" => Any[
                "Overview" => "man/materialparameters.md",
                "Permeability" => "man/permeability.md",
                "Heat Capacity" => "man/heatcapacity.md",
                "Conductivity" => "man/conductivity.md",
                "Latent heat" => "man/latentheat.md",
                "Radioactive heat" => "man/radioactiveheating.md",
                "Shear heating" => "man/shearheating.md",
                "Gravity" => "man/gravity.md",
                "Partial Melting" => "man/melting.md",
                "Density" => "man/density.md",
            ],
            "Constitutive Relationships" => Any[
            "Creep laws" => "man/creeplaws.md",
            "Custom rheology" => "man/customrheology.md",
            "Viscosity" => "man/viscosity.md",
            "Elasticity" => "man/elasticity.md",
            "Plasticity" => "man/plasticity.md",
            ],
            "Chemical Diffusion" => Any[
                "Computational routines" => "man/chemicaldiffusion.md",
                "Garnet" => "man/Garnet.md",
                "Melt" => "man/Melt.md",
                "Olivine" => "man/Olivine.md",
                "Rutile" => "man/Rutile.md",
            ],
            "TAS classification" => "man/TASclassification.md",
            "Zircon Ages" => "man/zirconages.md",
            "Phase Diagrams" => "man/phasediagrams.md",
            "Seismic Velocity" => "man/seismicvelocity.md",
            "1D Strength Envelope" => "man/strengthenvelope.md",
        ],
        "Plotting" => "man/plotting.md",
        "List of functions" => "man/listfunctions.md",
        "Contributing" => "man/contributing.md",
    ],
)

deploydocs(
    repo = "github.com/JuliaGeodynamics/GeoParams.jl",
    devbranch = "main",
    branch = "gh-pages",
    push_preview = true,
)
