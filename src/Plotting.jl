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
using .MaterialParameters.CreepLaw: CreepLawVariables, computeCreepLaw_TauII, AbstractCreepLaw
using .MaterialParameters.HeatCapacity: AbstractHeatCapacity, compute_heatcapacity
using .MaterialParameters.Conductivity: AbstractConductivity, compute_conductivity
using .MeltingParam: AbstractMeltingParam, compute_meltfraction

export 
    PlotStressStrainrate_CreepLaw,
    PlotHeatCapacity,
    PlotConductivity,
    PlotMeltFraction,
    PlotPhaseDiagram 


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
function PlotStressStrainrate_CreepLaw(x::AbstractCreepLaw; p=nothing, Strainrate=(1e-18/s,1e-12/s), CreatePlot::Bool=false)

    if isnothing(p); p = CreepLawParams();  end

    if isDimensional(x)==false
       error("The struct with Creep Law parameters: $(typeof(x)) should be in dimensional units for plotting. You can use Dimensionalize! to do that.")
    end

    # Define strainrate 
    Eps_II = range(ustrip(Strainrate[1])/s, stop=ustrip(Strainrate[2])/s, length=101)/s
    Tau_II = computeCreepLaw_TauII(Eps_II, x, p)                  # deviatoric stress

    # Transfer to GeoUnits
    Eps_II = GeoUnit(Eps_II);
    Tau_II = GeoUnit(Tau_II/1e6);

    if CreatePlot
        try 
            pl = plot(ustrip(Eps_II), ustrip(Tau_II), 
                xaxis=:log, xlabel=L"\textrm{deviatoric strain rate  } \dot{\varepsilon}_{II} \textrm{    [1/s]}", 
                yaxis=:log, ylabel=L"\textrm{deviatoric stress  }\tau_{II} \textrm{    [MPa]}",
                legend=false,show = true)
        catch
            error("It seems that you did not install, or did not load Plots.jl. For plotting, please add that with `add Plots` in the package manager and type `using Plots` before running this.")
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
        T = (273:10:1250)*K
    end

    Cp       =   compute_heatcapacity(T,cp)
    if length(Cp) == 1
        Cp = ones(size(T))*Cp
    end

    if isnothing(plt)
        plt = plot(ustrip(T), ustrip(Cp), label=lbl)
    else
        plt = plot!(ustrip(T), ustrip(Cp), label=lbl)
    end   
    plot!(plt,   xlabel="Temperature [$(unit(T[1]))]",
                 ylabel="Cp [$(unit(Cp[1]))]")
    gui(plt)

    return T,Cp, plt
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
function PlotConductivity(k::AbstractConductivity; T=nothing, P=nothing, plt=nothing, lbl=nothing)

    if isnothing(T)
        T = (273:10:1250)*K
    end
    if isnothing(P)
        P = 1e6Pa*ones(size(T))
    end

    Cond       =   compute_conductivity(P,T,k)
    if length(Cond) == 1
        Cond = ones(size(T))*Cond
    end

    if isnothing(plt)
        plt = plot(ustrip(T), ustrip(Cond), label=lbl)
    else
        plt = plot!(ustrip(T), ustrip(Cond), label=lbl)
    end   
    plot!(plt,   xlabel="Temperature [$(unit(T[1]))]",
                 ylabel="Thermal conductivity [$(unit(Cond[1]))]")
    gui(plt)

    return T,Cond, plt
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
julia> p        =  MeltingParam_Caricchi()
julia> T,phi,plt = PlotMeltFraction(p)
```
you can now save the figure to disk with:
```
julia> using Plots
julia> savefig(plt,"MeltFraction.png")
```

"""
function PlotMeltFraction(p::AbstractMeltingParam; T=nothing, P=nothing, plt=nothing, lbl=nothing)
    
    if isnothing(T)
        T = (500:10:1500)*K
    end
    T_C = ustrip(T) .- 273.15

    if isnothing(P) 
        P = 1e6Pa*ones(size(T))
    end
    if length(P) == 1
        P = P*ones(size(T))
    end

    phi       =   compute_meltfraction(P,T,p)
    if length(phi) == 1
        phi = ones(size(T))*phi
    end

    if isnothing(plt)
       plt = plot(T_C, ustrip(phi), label=lbl)
    else
       plt = plot!(T_C, ustrip(phi), label=lbl)
    end   
    plot!(plt,   xlabel="Temperature [C]",
                 ylabel="Melt Fraction \\Phi")
    gui(plt)

    return T,phi, plt
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
function PlotPhaseDiagram(p::PhaseDiagram_LookupTable, fieldn::Symbol; Tvec=nothing, Pvec=nothing)

    data = getfield( p, fieldn)
    if isnothing(Tvec)
        Tvec_K  =   data.itp.knots[1]
        Tvec    =   Tvec_K;
    else
        Tvec_K  =   Float64.(uconvert.(K,Tvec))
    end
    if isnothing(Pvec)
        Pvec_Pa =   data.itp.knots[2]
        Pvec    =   Pvec_Pa
    else
        Pvec_Pa =   Float64.(uconvert.(Pa,Pvec))
    end


    data_scalar = data(ustrip.(Tvec_K), ustrip.(Pvec_Pa) )


    plt = heatmap(ustrip.(Tvec), ustrip.(Pvec), data_scalar', 
                title  = string(fieldn),
                xlabel = "T [$(unit(Tvec[1]))]",
                ylabel = "P [$(unit(Pvec[1]))]", c=:batlow)

    display(plt)

    return  plt, data_scalar, Tvec, Pvec  
end
            
"""
    Tz_plot, Pz_plot, TauIIz_plot, T, P, z, TauII = PlotTPTauII(c::AbstractMaterialParam, d::AbstractDensity, g::AbstractMaterialParam, T::GeoUnit, z::GeoUnit, EpsII=(1e-18s^(-1.0)))
    
c=AbstractCreepLaw, d=AbstractDensity, g=AbstractGravity 

Is dependent on Plots and GeoParams.
Plots Temperature-Depth, Pressure-Depth and Deviatoric-stress-Depth diagrams.
Temperature-Depth is plotted as a function 'T' (x-axis) and 'z' (y-axis).
Pressure-Depth is plotted as a function 'P' (x-axis) and 'z' (y-axis).
Deviatoric-stress-Depth is plotted as a function 'TauII' (x-axis) and 'z' (y-axis).
The approximation of lithostatic pressure was used.
Was used with GeoUnit(range(...)) input for T and z with the same range()-length (see default values). CreepLaw (c), density (d) and gravity (g) still need to be put in.
Returns the T-z-plot, P-z-plot, TauII-z-plot, temperature (T), pressure (P), depth (z) and deviatoric stress (TauII).

EXAMPLE
=======
using Plots
using GeoParams

CharDim = GEO_units()
PhaseTest = SetMaterialParams(Name="Test", Phase=2, Density=ConstantDensity(ρ=3300kg/m^3), CreepLaws=(LinearViscous(), SetDislocationCreep("Dry Olivine | Hirth & Kohlstedt (2003)")), CharDim=CharDim)
PlotTPTauII(PhaseTest.CreepLaws[2], PhaseTest.Density[1], PhaseTest.Gravity[1])

----- Specific depth and temperature ranges can be used: 

using Plots
using GeoParams

CharDim = GEO_units()
PhaseTest = SetMaterialParams(Name="Test", Phase=2, Density=ConstantDensity(ρ=3300kg/m^3), CreepLaws=(LinearViscous(), SetDislocationCreep("Dry Olivine | Hirth & Kohlstedt (2003)")), CharDim=CharDim)
T = GeoUnit(range(303K , 1200K, length=30))
z = GeoUnit(range(1000m, 35000m, length=30))
PlotTPTauII(PhaseTest.CreepLaws[2], PhaseTest.Density[1], PhaseTest.Gravity[1], T, z)

-----

This code should produce 3 plots as describes above.
"""

"""

function PlotTPTauII(c::AbstractMaterialParam, d::AbstractDensity, g::AbstractMaterialParam, T::GeoUnit=GeoUnit(range(283K , 1200K, length=30)), z::GeoUnit=GeoUnit(range(0m, 30000m, length=30)), EpsII=(1e-18s^(-1.0)))

    P = zeros(length(z))
    CharDim = GEO_units()
    TauII = zeros(length(z))
    TauII_single = 0.0

    P[1:end] = dimensionalize(d.ρ, CharDim) * dimensionalize(g.g, CharDim) * z.val[1:end]   # simple shear approximation

    for i = 1:length(z)
        Temp = T[i]
        Pr = P[i]
        v = CreepLawVariables(T=nondimensionalize(Temp, CharDim), P=nondimensionalize(GeoUnit(Pr, Pa, true), CharDim))
        TauII_single = dimensionalize(GeoUnit(computeCreepLaw_TauII(nondimensionalize(EpsII[1], CharDim),c,v), MPa, false), CharDim)                                                     #deviatoric stress
        TauII[i] = TauII_single
    end

    try
        Tz_plot     = plot(T.val.-273.15, z.val.*1e-3, yflip=true, xmirror=true, color=:red, lab="geotherm", ylabel="Depth [km]", xlabel="Temperature [°C]", legend=:bottom)
        Pz_plot     = plot(P.*1e-6, z.val.*1e-3, yflip=true, xmirror=true, lab="geobar", ylabel="Depth [km]", xlabel="Pressure [MPa]", legend=:bottom)
        TauIIz_plot = scatter(TauII, z.val.*1e-3, yflip=true, xlim=(0,1000), xmirror=true, lab="stress", ylabel="Depth [km]", xlabel="Deviatoric stress [MPa]", legend=:bottom)
        display(plot(Tz_plot, Pz_plot, TauIIz_plot))
        return Tz_plot, Pz_plot, TauIIz_plot, T, P, z, TauII
    catch
        error("It seems that you did not install, or did not load Plots.jl. For plotting, please add that with `add Plots` in the package manager and type `using Plots` before running this.")
    end
end

"""