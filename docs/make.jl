using Documenter, GeoParams, Makie
push!(LOAD_PATH, "../src/")

@info "Making documentation..."
makedocs(;
    sitename = "GeoParams.jl",
    authors = "Boris Kaus and contributors",
    # add extension modules here for the docs.
    modules = [
        GeoParams,
        isdefined(Base, :get_extension) ?
            Base.get_extension(GeoParams, :GeoParamsMakieExt) :
            GeoParams.GeoParamsMakieExt,
    ],
    format = Documenter.HTML(;
        assets = ["assets/favicon.ico"],
        prettyurls = get(ENV, "CI", nothing) == "true",
        mathengine =
            Documenter.MathJax3(
            Dict(  # use MathJax3 as engine for latex (to be able to reference equations)
                :loader => Dict("load" => ["[tex]/physics"]),
                :tex => Dict(
                    "inlineMath" => [["\$", "\$"], ["\\(", "\\)"]],
                    "tags" => "ams",
                )
            )
        )
    ), # easier local build
    pages = [
        "Home" => "index.md",
        "User Guide" => Any[
            "GeoUnit" => "man/geounit.md",
            "Nondimensionalization" => "man/nondimensionalize.md",
            "Material Parameters" => "man/materialparameters.md",
            "Density" => "man/density.md",
            "Creep laws" => "man/creeplaws.md",
            "Custom rheology" => "man/customrheology.md",
            "Viscosity" => "man/viscosity.md",
            "Elasticity" => "man/elasticity.md",
            "Plasticity" => "man/plasticity.md",
            "Chemical Diffusion" => Any[
                "Computational routines" => "man/chemicaldiffusion.md",
                "Garnet" => "man/Garnet.md",
                "Melt" => "man/Melt.md",
                "Olivine" => "man/Olivine.md",
                "Rutile" => "man/Rutile.md",
            ],
            "Permeability" => "man/permeability.md",
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
            "1D Strength Envelope" => "man/strengthenvelope.md",
        ],
        "Plotting" => "man/plotting.md",
        "List of functions" => "man/listfunctions.md",
        "Contributing" => "man/contributing.md",
    ],
)

deploydocs(; repo = "github.com/JuliaGeodynamics/GeoParams.jl.git", devbranch = "main")
