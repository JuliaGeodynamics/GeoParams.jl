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
using .MaterialParameters.ConstitutiveRelationships
using .MaterialParameters.HeatCapacity: AbstractHeatCapacity, compute_heatcapacity
using .MaterialParameters.Conductivity: AbstractConductivity, compute_conductivity
using .MeltingParam: AbstractMeltingParam, compute_meltfraction

export 
    PlotStrainrateStress,
    PlotStressStrainrate,
    PlotStrainrateViscosity,
    PlotHeatCapacity,
    PlotConductivity,
    PlotMeltFraction,
    PlotPhaseDiagram,
    Plot_TAS_diagram,
    Plot_ZirconAge_PDF


"""
    PlotStrainrateStress(x; args=(T=1000.0, P=0.0, d=1e-3, f=1.0), Strainrate=(1e-18,1e-12), plt=nothing)

Plots deviatoric stress versus deviatoric strain rate for a single creeplaw. 
    Note: if you want to create plots or use the `CreatePlot=true` option you need to install the `Plots.jl` package in julia
    which is not added as a dependency here (as it is a rather large dependency).

# Example 1    
```julia-repl
julia> pp   = SetDiffusionCreep("Dry Anorthite | Rybacki et al. (2006)")
DiffusionCreep: Name = Dry Anorthite | Rybacki et al. (2006), n=1.0, r=0.0, p=-3.0, A=1.258925411794166e-12 m³·⁰ Pa⁻¹·⁰ s⁻¹·⁰, E=460000.0 J mol⁻¹·⁰, V=2.4e-5 m³·⁰ mol⁻¹·⁰, Apparatus=1
julia> pp1   = SetDislocationCreep("Dry Anorthite | Rybacki et al. (2006)")
DislocationCreep: Name = Dry Anorthite | Rybacki et al. (2006), n=3.0, r=0.0, A=5.011872336272715e-6 Pa⁻³·⁰ s⁻¹·⁰, E=641000.0 J mol⁻¹·⁰, V=2.4e-5 m³·⁰ mol⁻¹·⁰, Apparatus=1
```
Next you can plot this with
```julia-repl
julia> using Plots;
julia> args=((T=900.0, d=100e-6), (;T=900.0))
((T = 900.0, d = 0.0001), (T = 900.0,))
julia> plt = PlotStrainrateStress((pp,pp1), args=args, Strainrate=(1e-22,1e-12))
```

The plot can be customized as 
```julia-repl
julia> plot(plt, title="Diffusion and Dislocation Creep for Anorthite")
```

which will generate the following plot
![subet1](./assets/img/Stress_Strainrate_DislocationDiffusion_Anorthite.png)


See the [Plots.jl](https://github.com/JuliaPlots/Plots.jl) package for more options.

"""
function PlotStrainrateStress(x; args=(T=1000.0, P=0.0, d=1e-3, f=1.0), Strainrate=(1e-18,1e-12), plt=nothing)

    n = 1
    if isa(x,Tuple)
        n = length(x)
    end

    if isnothing(plt)
        plot()      # new plot
    end
    for i=1:n   
        if isa(x,Tuple)
            p = x[i]
        else
            p = x;
        end
        if isa(args,Tuple)
            args_in = args[i]
        else
            args_in = args;
        end
      
        # Define strainrate 
        Eps_II = range(ustrip(Strainrate[1]), stop=ustrip(Strainrate[2]), length=101)
        Tau_II = zeros(size(Eps_II))

        compute_τII!(Tau_II, p, Eps_II, args_in)       # Compute stress

        η = Tau_II./(2 * Eps_II)                        # effective viscosity

        Tau_II_MPa = Tau_II./1e6;

        # Create Plot    
        if isa(p,Tuple)
            Name = ""
            Type = ""
        else
            Name = String(collect(p.Name))
            
            # determine type of creeplaw 
            Type = "$(typeof(p))"           # full name of type
            id = findfirst("{", Type)
            Type = Type[1:id[1]-1]

        end
        
        
        plt = plot!(Eps_II,  Tau_II_MPa, 
                    xaxis=:log, xlabel=L"\dot{\varepsilon}_{II} \textrm{[s}^{-1}\textrm{]}", 
                    yaxis=:log, ylabel=L"\tau_{II} \textrm{    [MPa]}",
                    label="$Type: $Name $args_in",
                    title="",
                    legendfont=font(4))
            
    end

    display(plt)
    return plt
    
end

"""
    PlotStressStrainrate(x; args=(T=1000.0, P=0.0, d=1e-3, f=1.0), Stress=(1e0,1e8), plt=nothing)

Same as `PlotStrainrateStress` but versus stress (in MPa) versus strainrate instead.

"""
function PlotStressStrainrate(x; args=(T=1000.0, P=0.0, d=1e-3, f=1.0), Stress=(1e0,1e8), plt=nothing)

    n = 1
    if isa(x,Tuple)
        n = length(x)
    end

    if isnothing(plt)
        plot()      # new plot
    end
    for i=1:n   
        if isa(x,Tuple)
            p = x[i]
        else
            p = x;
        end
        if isa(args,Tuple)
            args_in = args[i]
        else
            args_in = args;
        end
        
        if isDimensional(p)==false
            error("The struct with Creep Law parameters: $(typeof(x)) should be in dimensional units for plotting. You can use Dimensionalize! to do that.")
        end

        # Define strainrate 
        Tau_II_MPa  = range(ustrip(Stress[1]), stop=ustrip(Stress[2]), length=101)
        Tau_II      = Tau_II_MPa.*1e6
        Eps_II      = zeros(size(Tau_II))

        compute_εII!(Eps_II, p, Tau_II, args_in)       # Compute strainrate
        
        η = Tau_II./(2 * Eps_II)                        # effective viscosity

        # Create Plot    
        Name = String(collect(p.Name))
        
        # determine type of creeplaw 
        Type = "$(typeof(p))"           # full name of type
        id = findfirst("{", Type)
        Type = Type[1:id[1]-1]

        plt = plot!(Tau_II_MPa,  Eps_II, 
                    xaxis=:log, ylabel=L"\tau_{II} \textrm{    [MPa]}",
                    yaxis=:log, xlabel=L"\dot{\varepsilon}_{II} \textrm{[s}^{-1}\textrm{]}", 
                    label="$Type: $Name $args_in",
                    title="",
                    legendfont=font(4))
            
    end

    display(plt)
    return plt
    
end

"""
    PlotStrainrateViscosity(x; args=(T=1000.0, P=0.0, d=1e-3, f=1.0), Strainrate=(1e-18,1e-12), plt=nothing)

Same as `PlotStrainrateStress` but versus viscosity instead of stress.

"""
function PlotStrainrateViscosity(x; args=(T=1000.0, P=0.0, d=1e-3, f=1.0), Strainrate=(1e-18,1e-12), plt=nothing)

    n = 1
    if isa(x,Tuple)
        n = length(x)
    end

    if isnothing(plt)
        plot()      # new plot
    end
    for i=1:n   
        if isa(x,Tuple)
            p = x[i]
        else
            p = x;
        end
        if isa(args,Tuple)
            args_in = args[i]
        else
            args_in = args;
        end
        
        if isDimensional(p)==false
            error("The struct with Creep Law parameters: $(typeof(x)) should be in dimensional units for plotting. You can use Dimensionalize! to do that.")
        end

        # Define strainrate 
        Eps_II = range(ustrip(Strainrate[1]), stop=ustrip(Strainrate[2]), length=101)
        Tau_II = zeros(size(Eps_II))

        compute_τII!(Tau_II, p, Eps_II, args_in)       # Compute stress

        η = Tau_II./(2 * Eps_II)                        # effective viscosity

        Tau_II_MPa = Tau_II./1e6;

        # Create Plot    
        Name = String(collect(p.Name))
        
        # determine type of creeplaw 
        Type = "$(typeof(p))"           # full name of type
        id = findfirst("{", Type)
        Type = Type[1:id[1]-1]

        plt = plot!(Eps_II,  η, 
                    xaxis=:log, xlabel=L"\dot{\varepsilon}_{II} \textrm{    [1/s]}", 
                    yaxis=:log, ylabel=L"\eta \textrm{    [Pa S]}",
                    label="$Type: $Name $args_in",
                    title="",
                    legendfont=font(4))
            
    end

    display(plt)
    return plt
    
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
        T = collect(273.:10:1250)*K
    end

    args = (;T=ustrip.(T))
    Cp = zeros(size(T))
    compute_heatcapacity!(Cp, cp, args)
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
        T = collect(273.:10:1250)*K
    end
    if isnothing(P)
        P = 1e6Pa*ones(size(T))
    end

    args = (;T=ustrip.(T))
    Cond = zeros(size(T))
    
    compute_conductivity!(Cond, k, args)
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
julia> p          =  MeltingParam_Caricchi()
julia> T,phi,dϕdT =  PlotMeltFraction(p)
```
you can now save the figure to disk with:
```
julia> using Plots
julia> savefig(plt,"MeltFraction.png")
```

"""
function PlotMeltFraction(p::AbstractMeltingParam; T=nothing, P=nothing, plt=nothing, lbl=nothing)
    
    if isnothing(T)
        T = (500.:10:1500.)*K
    end
    T_C = ustrip(T) .- 273.15

    if isnothing(P) 
        P = 1e6Pa*ones(size(T))
    end
    if length(P) == 1
        P = P*ones(size(T))
    end

    phi = ones(size(T))
    dϕdT= zeros(size(T))
    args=(;T=ustrip.(Vector(T)))
    compute_meltfraction!(phi, p, args)
    compute_dϕdT!(dϕdT, p, args)

    if isnothing(plt)
       plt1 = plot(T_C, ustrip(phi), label=lbl, ylabel="Melt Fraction \\Phi")
       plt2 = plot(T_C, ustrip(dϕdT), label=lbl, ylabel="d\\Phi / dT")
    else
       plt1 = plot!(T_C, ustrip(phi), label=lbl, ylabel="Melt Fraction \\Phi")
       plt2 = plot!(T_C, ustrip(phi), label=lbl, ylabel="d\\Phi / dT")
    end   
    plt = plot!(plt1,plt2,   xlabel="Temperature [C]", layout=(2,1))
                 
    gui(plt)

    return T,phi,dϕdT
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
	plt = Plot_ZirconAge_PDF(time_Ma, PDF_zircons, time_Ma_average, PDF_zircon_average)

Creates a plot of the Zircon Age probability density function from the parameters in a simulation
"""
function Plot_ZirconAge_PDF(time_Ma, PDF_zircons, time_Ma_average, PDF_zircon_average)

	plt = Plots.plot(time_Ma[1], PDF_zircons[1], color=:lightgray,linewidth=0.1, 
				xlabel="Time [Ma]", ylabel="probability []", title = "Zircon age probability distribution", legend=:none)
	for i in 2:length(PDF_zircons)
		plt = Plots.plot!(time_Ma[i], PDF_zircons[i], color=:lightgray,linewidth=0.1)
	end
	Plots.plot!(time_Ma_average, PDF_zircon_average, color=:black,linewidth=2.)
	
	display(plt)

	return plt
end

"""
	plt = Plot_TAS_diagram()

Creates a TAS diagram plot
"""
function Plot_TAS_diagram(displayLabel=nothing)

    if isnothing(displayLabel)
        displayLabel = 1;
    end

    # get TAS diagram data from TASclassification routine
    ClassTASdata    = TASclassificationData();
    @unpack litho, n_ver,  ver = ClassTASdata

    plt = Plots.plot(0, 0, xlabel="SiO2 [wt%]", ylabel="Na2O+K2O [wt%]", title = "TAS Diagram")

    n_poly  = size(litho,2);
    shift   = 1;
    for poly=1:n_poly

            x = sum(ver[shift:shift+n_ver[poly]-1,1])/n_ver[poly];
            y = sum(ver[shift:shift+n_ver[poly]-1,2])/n_ver[poly];

            plt = Plots.plot!(Shape(ver[shift:shift+n_ver[poly]-1,1], ver[shift:shift+n_ver[poly]-1,2]), 
            c = :transparent, xlims=(35,100),xticks=35:5:100, 
            ylims=(0,16),yticks=0:2:16, legend = false,
            )
            if displayLabel == 1
                annotate!(x,y,  (poly,:topleft,:blue,8))
            end

            shift  += n_ver[poly]
    end
    if displayLabel == 1
        for i=1:n_poly
            annotate!(86,16-i*3/4,  (string(i)*": "*litho[i],:left,:black,6))
        end
    end
	display(plt)

	return plt
end
