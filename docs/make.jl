using Documenter
using DocumenterVitepress
using GeoParams, Makie
# push!(LOAD_PATH, "../src/")

# Unitful defaults to unicode exponents (e.g. `m⁻³·⁰`) only on macOS unless this
# is set explicitly. The doctests in the docstrings expect the unicode form, so
# force it here to keep results identical across the OSes that build the docs.
ENV["UNITFUL_FANCY_EXPONENTS"] = "true"

DocMeta.setdocmeta!(GeoParams, :DocTestSetup, :(using GeoParams); recursive = true)

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
    checkdocs = :exports,
    warnonly = [:missing_docs],
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
                "Solubility" => "man/solubility.md",
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
        "Parameter tables" => "man/tables.md",
        "List of functions" => "man/listfunctions.md",
        "Contributing" => "man/contributing.md",
    ],
)

DocumenterVitepress.deploydocs(
    repo = "github.com/JuliaGeodynamics/GeoParams.jl",
    devbranch = "main",
    branch = "gh-pages",
    push_preview = true,
)
