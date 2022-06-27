"""
    This provides a few plotting routines, for example, for CreepLaws
"""

using LaTeXStrings
using Unitful
using Parameters
using ..Units
using ..MaterialParameters
using ..MeltingParam
using .Plots

using GeoParams: AbstractMaterialParam, AbstractMaterialParamsStruct
using .MaterialParameters.CreepLaw:
    CreepLawVariables, computeCreepLaw_TauII, AbstractCreepLaw
using .MaterialParameters.HeatCapacity: AbstractHeatCapacity, compute_heatcapacity
using .MaterialParameters.Conductivity: AbstractConductivity, compute_conductivity
using .MeltingParam: AbstractMeltingParam, compute_meltfraction

export PlotStressStrainrate_CreepLaw,
    PlotHeatCapacity,
    PlotConductivity,
    PlotMeltFraction,
    PlotPhaseDiagram,
    Plot_TAS_diagram,
    Plot_ZirconAge_PDF

"""
    PlotStressStrainrate_CreepLaw(x::AbstractCreepLaw; p::CreepLawParams=nothing, Strainrate=(1e-18,1e-12), CreatePlot::Bool=false)

Plots deviatoric stress versus deviatoric strain rate for a single creeplaw. 
    Note: if you want to create plots or use the `CreatePlot=true` option you need to install the `Plots.jl` package in julia
    which is not added as a dependency here (as it is a rather large dependency).

# Example 1    
```julia-repl
julia> x=LinearViscous()
Linear viscosity: η=1.0e20 Pa s
julia> Tau_II, Eps_II,  = PlotStressStrainrate_CreepLaw(x);
```
Next you can plot this with
```julia-repl
julia> using Plots;
julia> plot(ustrip(Eps_II),ustrip(Tau_II), xaxis=:log, yaxis=:log,xlabel="strain rate [1/s]",ylabel="Dev. Stress [MPa]")
```
Note that `ustrip` removes the units of the arrays, as many of the plotting packages don't know how to deal with that.

You could also have done:
```julia-repl
julia> using Plots;
julia> Tau_II, Eps_II, pl = PlotStressStrainrate_CreepLaw(x,CreatePlot=true);
```
which will generate the following plot
![subet1](./assets/img/Stress_Strainrate_LinearViscous.png)

The plot can be customized as 
```julia-repl
julia> plot(pl, title="Linear viscosity", linecolor=:red)
```
See the [Plots.jl](https://github.com/JuliaPlots/Plots.jl) package for more options.

"""
function PlotStressStrainrate_CreepLaw(
    x::AbstractCreepLaw;
    p=nothing,
    Strainrate=(1e-18 / s, 1e-12 / s),
    CreatePlot::Bool=false,
)
    if isnothing(p)
        p = CreepLawParams()
    end

    if isDimensional(x) == false
        error(
            "The struct with Creep Law parameters: $(typeof(x)) should be in dimensional units for plotting. You can use Dimensionalize! to do that.",
        )
    end

    # Define strainrate 
    Eps_II =
        range(ustrip(Strainrate[1]) / s; stop=ustrip(Strainrate[2]) / s, length=101) / s
    Tau_II = computeCreepLaw_TauII(Eps_II, x, p)                  # deviatoric stress

    # Transfer to GeoUnits
    Eps_II = GeoUnit(Eps_II)
    Tau_II = GeoUnit(Tau_II / 1e6)

    if CreatePlot
        try
            pl = plot(
                ustrip(Eps_II),
                ustrip(Tau_II);
                xaxis=:log,
                xlabel=L"\textrm{deviatoric strain rate  } \dot{\varepsilon}_{II} \textrm{    [1/s]}",
                yaxis=:log,
                ylabel=L"\textrm{deviatoric stress  }\tau_{II} \textrm{    [MPa]}",
                legend=false,
                show=true,
            )
        catch
            error(
                "It seems that you did not install, or did not load Plots.jl. For plotting, please add that with `add Plots` in the package manager and type `using Plots` before running this.",
            )
        end

        return Tau_II, Eps_II, pl
    else
        return Tau_II, Eps_II
    end
end

"""
    T,Cp,plt = PlotHeatCapacity(cp::AbstractHeatCapacity; T=nothing, plt=nothing, lbl=nothing)

Creates a plot of temperature `T` vs. heat capacity, as specified in cp (which can be temperature-dependent).

# Optional parameters
- T: temperature range
- plt: a previously generated plotting object
- lbl: label of the curve

# Example
```
julia> cp = T_HeatCapacity_Whittacker()
julia> T,Cp,plt = PlotHeatCapacity(cp)
```
you can now save the figure to disk with:
```
julia> using Plots
julia> savefig(plt,"Tdependent_heatcapacity.png")
```

"""
function PlotHeatCapacity(cp::AbstractHeatCapacity; T=nothing, plt=nothing, lbl=nothing)
    if isnothing(T)
        T = collect(273.0:10:1250) * K
    end

    args = (; T=ustrip.(T))
    Cp = zeros(size(T))
    compute_heatcapacity!(Cp, cp, args)
    if length(Cp) == 1
        Cp = ones(size(T)) * Cp
    end

    if isnothing(plt)
        plt = plot(ustrip(T), ustrip(Cp); label=lbl)
    else
        plt = plot!(ustrip(T), ustrip(Cp); label=lbl)
    end
    plot!(plt; xlabel="Temperature [$(unit(T[1]))]", ylabel="Cp [$(unit(Cp[1]))]")
    gui(plt)

    return T, Cp, plt
end

"""
    T,Kk,plt = PlotConductivity(cp::AbstractConductivity; T=nothing, plt=nothing, lbl=nothing)

Creates a plot of temperature `T` vs. thermal conductivity, as specified in `k` (which can be temperature-dependent).

# Optional parameters
- `T`: temperature range
- `plt`: a previously generated plotting object
- `lbl`: label of the curve

# Example
```
julia> k = T_Conductivity_Whittacker()
julia> T,KK,plt = PlotConductivity(k)
```
you can now save the figure to disk with:
```
julia> using Plots
julia> savefig(plt,"Tdependent_conductivity.png")
```

"""
function PlotConductivity(
    k::AbstractConductivity; T=nothing, P=nothing, plt=nothing, lbl=nothing
)
    if isnothing(T)
        T = collect(273.0:10:1250) * K
    end
    if isnothing(P)
        P = 1e6Pa * ones(size(T))
    end

    args = (; T=ustrip.(T))
    Cond = zeros(size(T))

    compute_conductivity!(Cond, k, args)
    if length(Cond) == 1
        Cond = ones(size(T)) * Cond
    end

    if isnothing(plt)
        plt = plot(ustrip(T), ustrip(Cond); label=lbl)
    else
        plt = plot!(ustrip(T), ustrip(Cond); label=lbl)
    end
    plot!(
        plt;
        xlabel="Temperature [$(unit(T[1]))]",
        ylabel="Thermal conductivity [$(unit(Cond[1]))]",
    )
    gui(plt)

    return T, Cond, plt
end

"""
    T,phi,plt = PlotMeltFraction(p::AbstractMeltingParam; T=nothing, plt=nothing, lbl=nothing)

Creates a plot of temperature `T` vs. melt fraction, as specified in `p`. 
The 1D curve can be evaluated at a specific pressure `P` which can be given as a scalar or as an array of the same size as `T`

# Optional parameters
- `T`: temperature range
- `P`: pressure 
- `plt`: a previously generated plotting object
- `lbl`: label of the curve

# Example
```
julia> p          =  MeltingParam_Caricchi()
julia> T,phi,dϕdT =  PlotMeltFraction(p)
```
you can now save the figure to disk with:
```
julia> using Plots
julia> savefig(plt,"MeltFraction.png")
```

"""
function PlotMeltFraction(
    p::AbstractMeltingParam; T=nothing, P=nothing, plt=nothing, lbl=nothing
)
    if isnothing(T)
        T = (500.0:10:1500.0) * K
    end
    T_C = ustrip(T) .- 273.15

    if isnothing(P)
        P = 1e6Pa * ones(size(T))
    end
    if length(P) == 1
        P = P * ones(size(T))
    end

    phi = ones(size(T))
    dϕdT = zeros(size(T))
    args = (; T=ustrip.(Vector(T)))
    compute_meltfraction!(phi, p, args)
    compute_dϕdT!(dϕdT, p, args)

    if isnothing(plt)
        plt1 = plot(T_C, ustrip(phi); label=lbl, ylabel="Melt Fraction \\Phi")
        plt2 = plot(T_C, ustrip(dϕdT); label=lbl, ylabel="d\\Phi / dT")
    else
        plt1 = plot!(T_C, ustrip(phi); label=lbl, ylabel="Melt Fraction \\Phi")
        plt2 = plot!(T_C, ustrip(phi); label=lbl, ylabel="d\\Phi / dT")
    end
    plt = plot!(plt1, plt2; xlabel="Temperature [C]", layout=(2, 1))

    gui(plt)

    return T, phi, dϕdT
end

"""
    plt, data, Tvec, Pvec = PlotPhaseDiagram(p::PhaseDiagram_LookupTable; fieldname::Symbol, Tvec=nothing, Pvec=nothing)

Plots a phase diagram as a function of `T` (x-axis) and `P` (y-axis).
We either use the default ranges of the diagram, or you can specify the temperature and pressure ranges (while specifying units).
The return arguments are the plotting object `plt` (so you can modify properties) as well as the data that is being plotted

Example
=======
```julia
julia> PD_data =  Read_LaMEM_Perple_X_Diagram("Peridotite.in")
Perple_X/LaMEM Phase Diagram Lookup Table: 
                      File    :   Peridotite.in
                      T       :   293.0 - 1573.000039
                      P       :   1.0e7 - 2.9999999944e9
                      fields  :   :meltRho, :meltRho, :meltFrac, :rockRho, :Rho, :rockVp
                                  :rockVs, :rockVpVs, :meltVp, :meltVs, :meltVpVs
                                  :Vp, :Vs, :VpVs, :cpxFrac
julia> PlotPhaseDiagram(PD_data,:meltFrac, Tvec=(100:1:1400).*C, Pvec=(.1:.1:30).*kbar )
```
This will generate the following plot
![subet2](./assets/img/PhaseDiagram.png)

You can also use the default pressure/temperature ranges in the diagrams:
```julia
julia> PlotPhaseDiagram(PD_data,:Rho)
```

"""
function PlotPhaseDiagram(
    p::PhaseDiagram_LookupTable, fieldn::Symbol; Tvec=nothing, Pvec=nothing
)
    data = getfield(p, fieldn)
    if isnothing(Tvec)
        Tvec_K = data.itp.knots[1]
        Tvec = Tvec_K
    else
        Tvec_K = Float64.(uconvert.(K, Tvec))
    end
    if isnothing(Pvec)
        Pvec_Pa = data.itp.knots[2]
        Pvec = Pvec_Pa
    else
        Pvec_Pa = Float64.(uconvert.(Pa, Pvec))
    end

    data_scalar = data(ustrip.(Tvec_K), ustrip.(Pvec_Pa))

    plt = heatmap(
        ustrip.(Tvec),
        ustrip.(Pvec),
        data_scalar';
        title=string(fieldn),
        xlabel="T [$(unit(Tvec[1]))]",
        ylabel="P [$(unit(Pvec[1]))]",
        c=:batlow,
    )

    display(plt)

    return plt, data_scalar, Tvec, Pvec
end

"""
	plt = Plot_ZirconAge_PDF(time_Ma, PDF_zircons, time_Ma_average, PDF_zircon_average)

Creates a plot of the Zircon Age probability density function from the parameters in a simulation
"""
function Plot_ZirconAge_PDF(time_Ma, PDF_zircons, time_Ma_average, PDF_zircon_average)
    plt = Plots.plot(
        time_Ma[1],
        PDF_zircons[1];
        color=:lightgray,
        linewidth=0.1,
        xlabel="Time [Ma]",
        ylabel="probability []",
        title="Zircon age probability distribution",
        legend=:none,
    )
    for i in 2:length(PDF_zircons)
        plt = Plots.plot!(time_Ma[i], PDF_zircons[i]; color=:lightgray, linewidth=0.1)
    end
    Plots.plot!(time_Ma_average, PDF_zircon_average; color=:black, linewidth=2.0)

    display(plt)

    return plt
end

"""
	plt = Plot_TAS_diagram()

Creates a TAS diagram plot
"""
function Plot_TAS_diagram(displayLabel=nothing)
    if isnothing(displayLabel)
        displayLabel = 1
    end

    # get TAS diagram data from TASclassification routine
    ClassTASdata = TASclassificationData()
    @unpack litho, n_ver, ver = ClassTASdata

    plt = Plots.plot(
        0, 0; xlabel="SiO2 [wt%]", ylabel="Na2O+K2O [wt%]", title="TAS Diagram"
    )

    n_poly = size(litho, 2)
    shift = 1
    for poly in 1:n_poly
        x = sum(ver[shift:(shift + n_ver[poly] - 1), 1]) / n_ver[poly]
        y = sum(ver[shift:(shift + n_ver[poly] - 1), 2]) / n_ver[poly]

        plt = Plots.plot!(
            Shape(
                ver[shift:(shift + n_ver[poly] - 1), 1],
                ver[shift:(shift + n_ver[poly] - 1), 2],
            );
            c=:transparent,
            xlims=(35, 100),
            xticks=35:5:100,
            ylims=(0, 16),
            yticks=0:2:16,
            legend=false,
        )
        if displayLabel == 1
            annotate!(x, y, (poly, :topleft, :blue, 8))
        end

        shift += n_ver[poly]
    end
    if displayLabel == 1
        for i in 1:n_poly
            annotate!(86, 16 - i * 3 / 4, (string(i) * ": " * litho[i], :left, :black, 6))
        end
    end
    display(plt)

    return plt
end
