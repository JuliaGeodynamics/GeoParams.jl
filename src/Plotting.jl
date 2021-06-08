"""
    This provides a few plotting routines, for example, for CreepLaws
"""

module Plotting
using LaTeXStrings
using Unitful
using Parameters
using ..Units
using ..MaterialParameters

using GeoParams: AbstractMaterialParam, AbstractMaterialParamsStruct
using .MaterialParameters.CreepLaw: CreepLawParams, ComputeCreepLaw_TauII, AbstractCreepLaw


export 
    PlotStressStrainrate_CreepLaw 


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
    Tau_II = ComputeCreepLaw_TauII(Eps_II, x, p)                  # deviatoric stress

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


end


