# This provides a few plotting routines, for example, for CreepLaws

# These functions are implemented in this file
import GeoParams: PlotStrainrateStress,
    PlotStressStrainrate,
    PlotStrainrateViscosity,
    PlotStressViscosity,
    PlotHeatCapacity,
    PlotConductivity,
    PlotMeltFraction,
    PlotPhaseDiagram,
    Plot_TAS_diagram,
    Plot_ZirconAge_PDF,
    PlotDeformationMap,
    PlotStressTime_0D,
    PlotPressureStressTime_0D,
    PlotDiffusionCoefArrhenius

# Make all `export`ed names from GeoParams.jl available
using GeoParams
# We also need the following un-`export`ed names
using GeoParams: AbstractTempStruct, DiffusionData, ptr2string
using GeoParams: AbstractMaterialParam, AbstractMaterialParamsStruct
using GeoParams.MaterialParameters.ConstitutiveRelationships
using GeoParams.MaterialParameters.HeatCapacity: AbstractHeatCapacity, compute_heatcapacity
using GeoParams.MaterialParameters.Conductivity: AbstractConductivity, compute_conductivity
using GeoParams.MeltingParam: AbstractMeltingParam, compute_meltfraction


using GeoParams.Units
using GeoParams.MaterialParameters
using GeoParams.MeltingParam

using Unitful
using Parameters

"""
    fig, ax, εII,τII = PlotStrainrateStress(x; Strainrate=(1e-18,1e-12), args =(T=1000.0, P=0.0, d=1e-3, f=1.0),
                                            linestyle=:solid, linewidth=1, color=nothing, label=nothing, title="",
                                            fig=nothing, filename=nothing, res=(1200, 1200), legendsize=15, labelsize=35)

Plots deviatoric stress versus deviatoric strain rate for a single or multiple creeplaws
    Note: if you want to create plots you need to install and load the `GLMakie.jl` package in julia.


# Example

First, we retrieve the data for anorthite creeplaws
```julia-repl
julia> import GeoParams.Diffusion, GeoParams.Dislocation
julia> pp  = SetDiffusionCreep(Diffusion.dry_anorthite_Rybacki_2006);
julia> pp1 = SetDislocationCreep(Dislocation.dry_anorthite_Rybacki_2006);
```
Next you can define each of the creeplaws individually, plus a combined diffusion & dislocation creep law:
```julia-repl
julia> v   = (pp,pp1,CompositeRheology(pp,pp1));
```
Next, define temperature to be `900K` and grainsize to be `100 μm` and create a default plot of the 3 mechanisms:
```julia-repl
julia> using GLMakie;
julia> args=(T=900.0, d=100e-6)
julia> PlotStrainrateStress(v, args=args, Strainrate=(1e-22,1e-15));
```

We have quite a few options to customize the look & feel of the plot:
```julia-repl
julia> fig,ax,εII,τII = PlotStrainrateStress(v, args=args, Strainrate=(1e-22,1e-15),
                                            color=(:red,:blue,:green), linewidth=(1,1,3), linestyle=(:dash,:dash,:solid), label=("diffusion creep","dislocation creep","diffusion+dislocation creep"),
                                            title="Dry Anorthite after Rybacki et al. (2006) for T=900K, d=100μm");
```

which will generate the following plot
![subet1](./assets/img/Stress_Strainrate_DislocationDiffusion_Anorthite.png)


See the [Makie.jl](https://makie.juliaplots.org/stable/) package for more options.

"""
function PlotStrainrateStress(
        x;
        Strainrate = (1.0e-18, 1.0e-12),
        args = (T = 1000.0, P = 0.0, d = 1.0e-3, f = 1.0),
        linestyle = :solid,
        linewidth = 1,
        color = nothing,
        label = nothing,
        title = "",
        fig = nothing,
        filename = nothing,
        res = (1200, 1200),
        legendsize = 15,
        labelsize = 35,
    )
    n = 1
    if isa(x, Tuple)
        n = length(x)
    end
    if isnothing(fig)
        fig = Figure(; fontsize = 25, size = res)
    end
    ax = Axis(
        fig[1, 1];
        yscale = log10,
        xscale = log10,
        xlabel = L"Deviatoric strain rate $\dot{ε}_{II}$ [1/s]",
        ylabel = L"Deviatoric stress $\tau_{II}$ [MPa]",
        xlabelsize = labelsize,
        ylabelsize = labelsize,
        title = title,
    )

    Eps_II = []
    Tau_II_MPa = []
    for i in 1:n

        # This allows plotting different curves on the same plot
        if isa(x, Tuple)
            p = x[i]
        else
            p = x
        end

        # This way we can set different args for every input argument in the tuple
        if isa(args, Tuple)
            args_in = args[i]
        else
            args_in = args
        end

        # Define strainrate
        Eps_II =
            exp10.(
            range(
                ustrip(log10(Strainrate[1]));
                stop = ustrip(log10(Strainrate[2])),
                length = 101,
            )
        )
        Tau_II = zeros(size(Eps_II))

        # Compute stress
        #compute_τII!(Tau_II, p, Eps_II, args_in)
        for j in eachindex(Tau_II)
            Tau_II[j] = compute_τII(p, Eps_II[j], args)
        end

        Tau_II_MPa = Tau_II ./ 1.0e6

        # Retrieve plot arguments (label, color etc.)
        plot_args = ObtainPlotArgs(i, p, args_in, linewidth, linestyle, color, label)

        # Create plot:
        li = lines!(Eps_II, Tau_II_MPa)    # plot line

        # Customize line:
        customize_plot!(li, plot_args)
    end
    axislegend(ax; labelsize = legendsize)

    if !isnothing(filename)
        save(filename, fig)
    else
        display(fig)
    end

    return fig, ax, Eps_II, Tau_II_MPa
end

# Gelper function that simplifies customising the plots
function ObtainPlotArgs(i, p, args_in, linewidth, linestyle, color, label_in)
    if isa(linewidth, Tuple)
        linewidth_in = linewidth[i]
    else
        linewidth_in = linewidth
    end

    if isa(color, Tuple)
        color_in = color[i]
    else
        color_in = color
    end

    if isa(linestyle, Tuple)
        linestyle_in = linestyle[i]
    else
        linestyle_in = linestyle
    end

    # Create a label name from the input parameters
    if isa(p, Tuple)
        # Combined creep law
        Name = ""
        Type = ""
        label = "$Type: $Name $args_in"
    else
        #if haskey(p,"Name")
        #    Name = String(collect(p.Name))
        #    # determine type of creeplaw
        #    Type = "$(typeof(p))"           # full name of type
        #else
        Name = ""
        Type = ""
        #end
        #id = findfirst("{", Type)
        #Type = Type[1:(id[1] - 1)]

        label = "$Type: $Name $args_in"
    end

    # We can manually overrule the auto-generated label
    if !isnothing(label_in)
        if isa(label_in, Tuple)
            label = label_in[i]
        else
            label = label_in
        end
    end

    # Create NamedTuple with arguments
    args = (linewidth = linewidth_in, linestyle = linestyle_in, label = label, color = color_in)

    return args
end

# Internal function that customizes the plot
function customize_plot!(li, args)

    # Customize line:
    li.label = args.label
    li.linewidth = args.linewidth
    li.linestyle = args.linestyle
    return if !isnothing(args.color)
        li.color = args.color
    end
end

"""
    fig,ax,τII,εII =  PlotStressStrainrate(x; args=(T=1000.0, P=0.0, d=1e-3, f=1.0), Stress=(1e0,1e8))

Same as `PlotStrainrateStress` but with stress (in MPa) versus strainrate (in 1/s) instead.

Example
===

```julia
julia> import GeoParams.Dislocation, GeoParams.Diffusion;
julia> a1=SetDislocationCreep(Dislocation.wet_olivine_Hirth_2003);
julia> a2=SetDiffusionCreep(Diffusion.wet_olivine_Hirth_2003);
julia> x=CompositeRheology(a1,a2);
julia> fig,ax,τII,εII =  PlotStressStrainrate(x; args=(T=1000.0, P=0.0, d=1e-3, f=1.0), Stress=(1e0,1e8));
```

"""
function PlotStressStrainrate(
        x;
        args = (T = 1000.0, P = 0.0, d = 1.0e-3, f = 1.0),
        Stress = (1.0e0, 1.0e8),
        linestyle = :solid,
        linewidth = 1,
        color = nothing,
        label = nothing,
        title = "",
        fig = nothing,
        filename = nothing,
        res = (1200, 1200),
        legendsize = 15,
        labelsize = 35,
    )
    n = 1
    if isa(x, Tuple)
        n = length(x)
    end

    if isnothing(fig)
        fig = Figure(; fontsize = 25, size = res)
    end
    ax = Axis(
        fig[1, 1];
        yscale = log10,
        xscale = log10,
        xlabel = L"Deviatoric stress $\tau_{II}$ [MPa]",
        ylabel = L"Deviatoric strain rate $\dot{ε}_{II}$ [1/s]",
        xlabelsize = labelsize,
        ylabelsize = labelsize,
        title = title,
    )

    Eps_II = []
    Tau_II_MPa = []
    for i in 1:n
        if isa(x, Tuple)
            p = x[i]
        else
            p = x
        end
        if isa(args, Tuple)
            args_in = args[i]
        else
            args_in = args
        end

        # Define strainrate
        Tau_II_MPa = range(ustrip(Stress[1]); stop = ustrip(Stress[2]), length = 101)
        Tau_II = Tau_II_MPa .* 1.0e6
        Eps_II = zeros(size(Tau_II))


        # Compute stress
        #compute_εII!(Eps_II, p, Tau_II, args_in)       # Compute strainrate
        for j in eachindex(Tau_II)
            Eps_II[j] = compute_εII(p, Tau_II[j], args)
        end

        η = Tau_II ./ (2 * Eps_II)                        # effective viscosity

        # Retrieve plot arguments (label, color etc.)
        plot_args = ObtainPlotArgs(i, p, args_in, linewidth, linestyle, color, label)

        # Create plot:
        li = lines!(Tau_II_MPa, Eps_II)    # plot line

        # Customize plot:
        customize_plot!(li, plot_args)
    end

    axislegend(ax; labelsize = legendsize)

    if !isnothing(filename)
        save(filename, fig)
    else
        display(fig)
    end

    return fig, ax, Tau_II_MPa, Eps_II
end

"""
    fig, ax, εII, η = PlotStrainrateViscosity(x; args=(T=1000.0, P=0.0, d=1e-3, f=1.0), Strainrate=(1e-18,1e-12),
                                linestyle=:solid, linewidth=1, color=nothing, label=nothing, title="",
                                fig=nothing, filename=nothing, res=(1200, 1200), legendsize=15, labelsize=35)

Same as `PlotStrainrateStress` but versus viscosity instead of stress.

"""
function PlotStrainrateViscosity(
        x;
        args = (T = 1000.0, P = 0.0, d = 1.0e-3, f = 1.0),
        Strainrate = (1.0e-18, 1.0e-12),
        linestyle = :solid,
        linewidth = 1,
        color = nothing,
        label = nothing,
        title = "",
        fig = nothing,
        filename = nothing,
        res = (1200, 1200),
        legendsize = 15,
        labelsize = 35,
    )
    n = 1
    if isa(x, Tuple)
        n = length(x)
    end

    if isnothing(fig)
        fig = Figure(; fontsize = 25, size = res)
    end
    ax = Axis(
        fig[1, 1];
        yscale = log10,
        xscale = log10,
        xlabel = L"Deviatoric strain rate $\dot{ε}_{II}$ [1/s]",
        ylabel = L"Effective viscosity $\eta$ [Pa s]",
        xlabelsize = labelsize,
        ylabelsize = labelsize,
        title = title,
    )

    Eps_II = []
    η = []
    for i in 1:n
        if isa(x, Tuple)
            p = x[i]
        else
            p = x
        end
        if isa(args, Tuple)
            args_in = args[i]
        else
            args_in = args
        end

        # Define strainrate
        Eps_II =
            exp10.(
            range(
                ustrip(log10(Strainrate[1]));
                stop = ustrip(log10(Strainrate[2])),
                length = 101,
            )
        )
        Tau_II = zeros(size(Eps_II))

        # Compute stress
        #compute_τII!(Tau_II, p, Eps_II, args_in)
        for j in eachindex(Tau_II)
            Tau_II[j] = compute_τII(p, Eps_II[j], args)
        end

        η = Tau_II ./ (2 * Eps_II)                        # effective viscosity

        if maximum(η) ≈ minimum(η)
            # if all values are the same, plotting creates an error can occur, so we perturb the values a bit
            η[1] *= (1.0 - 1.0e-5)
        end

        Tau_II_MPa = Tau_II ./ 1.0e6

        # Retrieve plot arguments (label, color etc.)
        plot_args = ObtainPlotArgs(i, p, args_in, linewidth, linestyle, color, label)

        # Create plot:
        li = lines!(Eps_II, η)    # plot line

        # Customize line:
        customize_plot!(li, plot_args)
    end

    axislegend(ax; labelsize = legendsize)

    if !isnothing(filename)
        save(filename, fig)
    else
        display(fig)
    end

    return fig, ax, Eps_II, η
end

"""
    fig,ax,τII,η =  PlotStressViscosity(x; args=(T=1000.0, P=0.0, d=1e-3, f=1.0), Stress=(1e0,1e8),
                                    linestyle=:solid, linewidth=1, color=nothing, label=nothing, title="",
                                    fig=nothing, filename=nothing, res=(1200, 1200), legendsize=15, labelsize=35)


Same as `PlotStrainrateStress` but versus stress (in MPa) and viscosity (Pas) instead.

"""
function PlotStressViscosity(
        x;
        args = (T = 1000.0, P = 0.0, d = 1.0e-3, f = 1.0),
        Stress = (1.0e0, 1.0e8),
        linestyle = :solid,
        linewidth = 1,
        color = nothing,
        label = nothing,
        title = "",
        fig = nothing,
        filename = nothing,
        res = (1200, 1200),
        legendsize = 15,
        labelsize = 35,
    )
    n = 1
    if isa(x, Tuple)
        n = length(x)
    end

    if isnothing(fig)
        fig = Figure(; fontsize = 25, size = res)
    end
    ax = Axis(
        fig[1, 1];
        yscale = log10,
        xscale = log10,
        xlabel = L"Deviatoric stress $\tau_{II}$ [MPa]",
        ylabel = L"Effective viscosity $\eta$ [Pa s]",
        xlabelsize = labelsize,
        ylabelsize = labelsize,
        title = title,
    )

    η = []
    Tau_II_MPa = []
    for i in 1:n
        if isa(x, Tuple)
            p = x[i]
        else
            p = x
        end
        @show typeof(p)
        if isa(args, Tuple)
            args_in = args[i]
        else
            args_in = args
        end

        # Define strainrate
        Tau_II_MPa = range(ustrip(Stress[1]); stop = ustrip(Stress[2]), length = 101)
        Tau_II = Tau_II_MPa .* 1.0e6
        Eps_II = zeros(size(Tau_II))

        #compute_εII!(Eps_II, p, Tau_II, args_in)       # Compute strainrate
        for j in eachindex(Tau_II)
            Eps_II[j] = compute_εII(p, Tau_II[j], args)
        end

        η = Tau_II ./ (2 * Eps_II)                        # effective viscosity

        if maximum(η) ≈ minimum(η)
            # if all values are the same, plotting creates an error can occur, so we perturb the values a bit
            η[1] *= (1.0 - 1.0e-5)
        end

        # Retrieve plot arguments (label, color etc.)
        plot_args = ObtainPlotArgs(i, p, args_in, linewidth, linestyle, color, label)

        # Create plot:
        li = lines!(Tau_II_MPa, η)    # plot line

        # Customize plot:
        customize_plot!(li, plot_args)
    end

    axislegend(ax; labelsize = legendsize)

    if !isnothing(filename)
        save(filename, fig)
    else
        display(fig)
    end

    return fig, ax, Tau_II_MPa, η
end

"""
     fig,ax,T,Cp_vec = PlotHeatCapacity(Cp::AbstractHeatCapacity; T=nothing, plt=nothing, lbl=nothing)

Creates a plot of temperature `T` vs. heat capacity, as specified in Cp (which can be temperature-dependent).

# Optional parameters
- T: temperature range
- plt: a previously generated plotting object
- lbl: label of the curve

# Example
```
julia> Cp = T_HeatCapacity_Whittacker()
julia> fig,ax,T,Cp_vec = PlotHeatCapacity(Cp)
```
"""
function PlotHeatCapacity(
        x;
        args = (;),
        Stress = (1.0e0, 1.0e8),
        linestyle = :solid,
        linewidth = 1,
        color = nothing,
        label = nothing,
        title = "",
        fig = nothing,
        filename = nothing,
        res = (1200, 1200),
        legendsize = 15,
        labelsize = 35,
        T = nothing
    )

    if isnothing(T)
        T = collect(273.0:10:1250) * K
    end

    n = 1
    if isa(x, Tuple)
        n = length(x)
    end

    if isnothing(fig)
        fig = Figure(; fontsize = 25, size = res)
    end
    ax = Axis(
        fig[1, 1];
        xlabel = "Temperature [K]",
        ylabel = "Heat Capacity [J kg⁻¹·⁰ K⁻¹·⁰]",
        xlabelsize = labelsize,
        ylabelsize = labelsize,
        title = title,
    )

    args = (; T = ustrip.(T))
    Cp1 = zeros(size(T))
    #compute_heatcapacity!(Cp1, Cp, args)
    #if length(Cp) == 1
    #    Cp1 = ones(size(T)) * Cp1
    #end

    #if isnothing(plt)
    #    plt = plot(ustrip(T), ustrip(Cp); label=lbl)
    #else
    #    plt = plot!(ustrip(T), ustrip(Cp); label=lbl)
    #end
    #lines(plt; xlabel="Temperature [$(unit(T[1]))]", ylabel="Cp [$(unit(Cp[1]))]")
    #gui(plt)

    for i in 1:n
        if isa(x, Tuple)
            p = x[i]
        else
            p = x
        end

        if isa(args, Tuple)
            args_in = args[i]
        else
            args_in = args
        end

        compute_heatcapacity!(Cp1, p, args)
        # Retrieve plot arguments (label, color etc.)
        plot_args = ObtainPlotArgs(i, p, args_in, linewidth, linestyle, color, label)

        # Create plot:
        li = lines!(ustrip.(T), Cp1)    # plot line

        # Customize plot:
        customize_plot!(li, plot_args)
    end

    #axislegend(ax; labelsize=legendsize)

    if !isnothing(filename)
        save(filename, fig)
    else
        display(fig)
    end

    return fig, ax, T, Cp1
end


# TO BE FIXED
"""
    T,Kk,plt = PlotConductivity(Cp::AbstractConductivity; T=nothing, plt=nothing, lbl=nothing)

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
        k::AbstractConductivity; T = nothing, P = nothing, plt = nothing, lbl = nothing
    )
    if isnothing(T)
        T = collect(273.0:10:1250) * K
    end
    if isnothing(P)
        P = 1.0e6Pa * ones(size(T))
    end

    args = (; T = ustrip.(T))
    Cond = zeros(size(T))

    compute_conductivity!(Cond, k, args)
    if length(Cond) == 1
        Cond = ones(size(T)) * Cond
    end

    if isnothing(plt)
        plt = plot(ustrip(T), ustrip(Cond); label = lbl)
    else
        plt = plot!(ustrip(T), ustrip(Cond); label = lbl)
    end
    plot!(
        plt;
        xlabel = "Temperature [$(unit(T[1]))]",
        ylabel = "Thermal conductivity [$(unit(Cond[1]))]",
    )
    gui(plt)

    return T, Cond, plt
end

# TO BE FIXED
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
        p::AbstractMeltingParam; T = nothing, P = nothing, plt = nothing, lbl = nothing
    )
    if isnothing(T)
        T = (873.0:10:1500.0) * K
    end
    T_C = ustrip(T) .- 273.15

    if isnothing(P)
        P = 1.0e6Pa * ones(size(T))
    end
    if length(P) == 1
        P = P * ones(size(T))
    end

    phi = ones(size(T))
    dϕdT = zeros(size(T))
    args = (; T = ustrip.(Vector(T)))
    compute_meltfraction!(phi, p, args)
    compute_dϕdT!(dϕdT, p, args)

    if isnothing(plt)
        plt1 = plot(T_C, ustrip(phi); label = lbl, ylabel = "Melt Fraction \\Phi")
        plt2 = plot(T_C, ustrip(dϕdT); label = lbl, ylabel = "d\\Phi / dT")
    else
        plt1 = plot!(T_C, ustrip(phi); label = lbl, ylabel = "Melt Fraction \\Phi")
        plt2 = plot!(T_C, ustrip(phi); label = lbl, ylabel = "d\\Phi / dT")
    end
    plt = plot!(plt1, plt2; xlabel = "Temperature [C]", layout = (2, 1))

    gui(plt)

    return T, phi, dϕdT
end

# BROKEN
"""
    plt, data, Tvec, Pvec = PlotPhaseDiagram(p::PhaseDiagram_LookupTable; fieldname::Symbol, Tvec=nothing, Pvec=nothing)

Plots a phase diagram as a function of `T` (x-axis) and `P` (y-axis).
We either use the default ranges of the diagram, or you can specify the temperature and pressure ranges (while specifying units).
The return arguments are the plotting object `plt` (so you can modify properties) as well as the data that is being plotted

Example
=======
```julia
julia> PD_Data = PerpleX_LaMEM_Diagram("./test/test_data/Peridotite.in")
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
        p::PhaseDiagram_LookupTable, fieldn::Symbol; Tvec = nothing, Pvec = nothing
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
        title = string(fieldn),
        xlabel = "T [$(unit(Tvec[1]))]",
        ylabel = "P [$(unit(Pvec[1]))]",
        c = :batlow,
    )

    display(plt)

    return plt, data_scalar, Tvec, Pvec
end

"""
	plt = Plot_ZirconAge_PDF(time_Ma, PDF_zircons, time_Ma_average, PDF_zircon_average)

Creates a plot of the Zircon Age probability density function from the parameters in a simulation
"""
function Plot_ZirconAge_PDF(time_Ma, PDF_zircons, time_Ma_average, PDF_zircon_average)
    f = Figure()
    Axis(f[1, 1], xlabel = "Age [Myr]", ylabel = "Kernel density [ ]", title = "Zircon age probability distribution")
    for i in 1:length(PDF_zircons)
        lines!(time_Ma[i] / 1.0e6, PDF_zircons[i], color = "gray66", linewidth = 0.25)
    end
    lines!(time_Ma_average / 1.0e6, PDF_zircon_average, color = "grey0", linewidth = 2.0)
    xlims!(-1.0e5 / 1.0e6, 1.5e6 / 1.0e6)

    display(f)

    return f
end


"""
	plt = Plot_TAS_diagram(point; displayLabel=true,size=(1500,1500), fontsize=18)

Creates a TAS diagram plot

## Input
- `point::AbstractArray`: point to plot on the TAS diagram (SiO2, Na2O+K2O)

## Optional parameters
- `displayLabel::Bool`: display the label of the lithology
- `size::Tuple`: size of the plot
- `fontsz::Int`: fontsize of the plot
- `colormp::Symbol`: colormap to use

## Example
```julia
julia> point = [50.0 5.0;
                61.6 6.25]
julia> Plot_TAS_diagram(point)
```
"""
function Plot_TAS_diagram(point::AbstractArray{_T}; displayLabel = true, sz = (1500, 1500), fontsize = 18, colormap = :lipari) where {_T}
    # get TAS diagram data from TASclassification routine
    ClassTASdata = TASclassificationData()
    @unpack litho, n_ver, ver = ClassTASdata

    f = Figure(size = sz, fontsize = fontsize)
    p1 = GridLayout(f[1, 1])
    ax1 = Axis(
        p1[1, 1],
        xlabel = "SiO2 [wt%]",
        ylabel = "Na2O+K2O [wt%]",
        title = "TAS Diagram",
        aspect = 1,
        xticks = 35:5:100,
        yticks = 0:2:16,
    )
    n_poly = size(litho, 2)
    shift = 1

    for poly in 1:n_poly
        shift_poly = shift:(shift + n_ver[poly] - 1)
        x = sum(ver[i, 1] for i in shift_poly) / n_ver[poly]
        y = sum(ver[i, 2] for i in shift_poly) / n_ver[poly]
        ps = [Point2f(ver[i, :]) for i in shift_poly]

        poly!(ax1, ps, color = :white, strokecolor = :black, strokewidth = 1)
        if displayLabel
            text!(ax1, (x, y), text = "$poly")
        end
        shift += n_ver[poly]
    end

    # Ensure the number of colors matches the number of points
    num_points = size(point, 1)
    colormap = cgrad(colormap, num_points; categorical = true)
    for i in 1:num_points
        scatter!(
            ax1, [point[i, 1]], [point[i, 2]], color = colormap[i], markersize = 10, strokewidth = 0.5,
            strokecolor = :grey45, label = "Composition $i: ($(point[i, 1]), $(point[i, 2]))"
        )
    end
    xlims!(ax1, 35, 100)
    ylims!(ax1, 0, 16)

    if displayLabel
        p2 = GridLayout(f[1, 2])
        ax2 = Axis(
            p2[1, 1],
            bottomspinevisible = false,
            xgridvisible = false,
            ygridvisible = false,
            rightspinevisible = false,
            leftspinevisible = false,
            topspinevisible = false
        )
        for i in 1:n_poly
            text!(ax2, 0, 16 - i * 3 / 4, text = (string(i) * ": " * litho[i]))
        end

        xlims!(ax2, 0, 2)
        hidedecorations!(ax2)
        rowsize!(p2, 1, 700)
        colsize!(p2, 1, 225)
    end

    display(f)

    return f
end

"""
    fig,ax,τII,η =  PlotStressTime_0D(x;    args=(T=1000.0, P=0.0, d=1e-3, f=1.0),
                                            εII::Union{Number, AbstractVector},
                                            τ0=0,
                                            Time=(1e0,1e8), nt=100,
                                            t_vec=nothing,
                                            linestyle=:solid, linewidth=1, color=nothing, label=nothing, title="",
                                            fig=nothing, filename=nothing, res=(1200, 1200), legendsize=15, labelsize=35, position=:rt)


Creates a plot of stress vs. time, or stress vs. strain (in)
"""
function PlotStressTime_0D(
        x;
        args = (T = 1000.0, P = 0.0, d = 1.0e-3, f = 1.0),
        εII = 1.0e-15,
        τ0 = 0, τ_scale = 1.0e6,
        Time = (1.0e0, 1.0e8), nt = 100,
        t_vec = nothing, t_scale = 3600 * 24 * 365.25 * 1.0e6,
        verbose = false,
        linestyle = :solid,
        linewidth = 1,
        color = nothing,
        label = nothing,
        title = "",
        fig = nothing,
        filename = nothing,
        res = (1200, 1200),
        legendsize = 15,
        labelsize = 35,
        position = :rt
    )
    n = 1
    if isa(x, Tuple)
        n = length(x)
    end

    if isnothing(fig)
        fig = Figure(; fontsize = 25, size = res)
    end
    ax = Axis(
        fig[1, 1];
        ylabel = L"Deviatoric stress $\tau_{II}$ [MPa]",
        xlabel = "Time [Myr]",
        xlabelsize = labelsize,
        ylabelsize = labelsize,
        title = title,
    )

    t_vec = []
    Tau_II_MPa = []
    for i in 1:n
        if isa(x, Tuple)
            p = x[i]
        else
            p = x
        end
        if isa(args, Tuple)
            args_in = args[i]
        else
            args_in = args
        end

        # Compute
        t_vec, τ_vec = time_τII_0D(p, εII, args; t = Time, nt = nt, verbose = verbose)
        Tau_II_MPa = τ_vec / τ_scale

        # Retrieve plot arguments (label, color etc.)
        plot_args = ObtainPlotArgs(i, p, args_in, linewidth, linestyle, color, label)

        # Create plot:
        li = lines!(t_vec / t_scale, Tau_II_MPa)    # plot line

        # Customize plot:
        customize_plot!(li, plot_args)
    end

    axislegend(ax; labelsize = legendsize, position = position)

    if !isnothing(filename)
        save(filename, fig)
    else
        display(fig)
    end

    return fig, ax, Tau_II_MPa, t_vec
end

"""
	fig = PlotDeformationMap(v;    args=(P=0.0, d=1e-3, f=1.0),
                                σ = (1e-2, 1e8),                # in MPa
                                T = (10, 1000),                 # in C
                                ε = (1e-22, 1e-8),              # in 1/s
                                n = 400,                        # number of points
                                rotate_axes = false,            # flip x & y axes
                                strainrate = true,              # strainrate (otherwise stress)
                                viscosity = false,              # plot viscosity instead of strainrate/stress
                                boundaries = true,              # plot deformation boundaries
                                levels = 20,                    # number of contour levels
                                colormap=:viridis,              # colormap
                                filename=nothing,               # if you want to save this to file
                                fontsize=40,                    # fontsize of labels
                                res=(1200, 900))                # resolution in pixels

Creates a deformation mechanism map (T/εII vs. stress/viscosity or T/τII vs. strainrate/viscosity) for given (composite) rheology `v`

# Example
```julia
julia> import GeoParams.Diffusion, GeoParams.Dislocation
julia> v1 = SetDiffusionCreep(Diffusion.dry_anorthite_Rybacki_2006);
julia> v2 = SetDislocationCreep(Dislocation.dry_anorthite_Rybacki_2006);
julia> v = CompositeRheology(v1,v2)
julia> PlotDeformationMap(v, levels=100, colormap=:roma)
```
Next, let's plot viscosity and flip x & y axis:
```julia
julia> PlotDeformationMap(v, viscosity=true, rotate_axes=true)
```
Instead of plotting stress vs. T and computing strainrate, we can also provide strainrate/T and compute stress:
```julia
julia> PlotDeformationMap(v, strainrate=false)
```
Or plot viscosity but only add contours in a certain range:
```julia
julia> PlotDeformationMap(v,  strainrate=false, viscosity=true, levels=Vector(18:.25:24))
```

"""
function PlotDeformationMap(
        v;
        args = (P = 0.0, T = 1250, d = 3.0e-3, f = 1.0),
        d = (1.0e-6, 1.0e-1),               # in m
        σ = (1.0e-2, 1.0e8),                # in MPa
        T = (10, 1000),                 # in C
        ε = (1.0e-22, 1.0e-8),              # in 1/s
        n = 400,                        # number of points
        rotate_axes = false,            # flip x & y axes
        strainrate = true,              # strainrate (otherwise stress)
        viscosity = false,              # plot viscosity instead of strainrate/stress
        grainsize = false,              # plot strainrate with grainsize as x-axis
        depth = false,                  # plot strainrate with depth as x-axis
        boundaries = true,              # plot deformation boundaries
        levels = 30,                    # number of contour levels
        colormap = :viridis,
        filename = nothing,
        fontsize = 40,
        res = (1200, 900),
    )

    # allocating ticks
    xtick = 0

    # Parameters
    T_vec = Vector(range(T[1], T[2], n + 1)) .+ 273.15   # in K

    # this is used to determine the main deformation mechanism (and color that later)
    n_components = 1
    if isa(v, CompositeRheology)
        n_components = length(v.elements)
    else
        n_components = length(v)
    end

    d_vec = range(d[1], d[2], n + 1)
    if strainrate
        # compute ε as a function of τ and T

        σ_vec = 10.0 .^ Vector(range(log10(σ[1] * 1.0e6), log10(σ[2] * 1.0e6), n))        # in Pa, equally spaced in log10 space
        εII = zeros(n + 1, n)
        η = zeros(n + 1, n)
        mainDef = zeros(n + 1, n)     # indicates the main components
        for i in CartesianIndices(εII)
            τlocal = σ_vec[i[2]]
            if grainsize
                Tlocal = 1250
                dlocal = d_vec[i[1]]
                args_local = merge(args, (T = Tlocal, d = dlocal))
                xtick = log10(d[1]):0.1:(log10(d[2]) .+ 3.0)
            else
                dlocal = 3.0e-3
                Tlocal = T_vec[i[1]]
                args_local = merge(args, (T = Tlocal, d = dlocal))
            end
            εII[i] = compute_εII(v, τlocal, args_local)       # compute strainrate (1/s)
            ε_components = [ compute_εII(v[i], τlocal, args_local) for i in 1:n_components]
            ε_components = ε_components ./ sum(ε_components)
            mainDef[i] = argmax(ε_components)                 # index of max. strainrate
        end
        log_σ = log10.(σ_vec ./ 1.0e6)
    else
        # compute τ as a function of ε and T

        ε_vec = 10.0 .^ Vector(range(log10(ε[1]), log10(ε[2]), n))        # in Pa, equally spaced in log10 space

        τII = zeros(n + 1, n)
        η = zeros(n + 1, n)
        mainDef = zeros(n + 1, n)     # indicates the main components
        for i in CartesianIndices(τII)
            Tlocal = T_vec[i[1]]
            dlocal = d_vec[i[1]]
            εlocal = ε_vec[i[2]]
            args_local = merge(args, (T = Tlocal, d = dlocal))

            τII[i] = compute_τII(v, εlocal, args_local)       # compute stress (Pa)
            η[i] = τII[i] / (2 * εlocal)

            τ_components = [ compute_τII(v[i], εlocal, args_local) for i in 1:n_components]
            τ_components = τ_components ./ sum(τ_components)
            mainDef[i] = argmin(τ_components)                 # index of max. strainrate
        end
        log_ε = log10.(ε_vec)
    end
    T_plot = T_vec .- 273.15

    # determine axis of plot
    if strainrate
        x = T_plot
        xlabel = "T [°C]"
        y = log_σ
        ylabel = L"\log_{10}(\tau_{II}) [MPa]"
        label = L"\log_{10}(\varepsilon_{II}) [s^{-1}]"
        data = log10.(εII)
    else
        x = T_plot
        xlabel = "T [°C]"
        y = log_ε
        ylabel = L"\log_{10}({\varepsilon}_{II}) [s^{-1}]"
        label = L"\log_{10}(\tau_{II}) [MPa]"
        data = log10.(τII / 1.0e6)
    end
    if viscosity
        label = L"\log_{10}(\eta_{eff}) [Pa s]"
        data = log10.(η)
    end

    if rotate_axes
        x, y = y, x
        xlabel, ylabel = ylabel, xlabel
        data, mainDef = data', mainDef'
    end

    # Plotting with Makie
    fig = Figure(; fontsize = fontsize, size = res)

    ax = Axis(
        fig[1, 1],
        title = "Deformation mechanism map",
        xlabel = xlabel, xlabelsize = fontsize,
        xticks = -3:0.5:2,
        xminorticks = IntervalsBetween(5),
        xminorticksvisible = true,
        ylabel = ylabel, ylabelsize = fontsize,
        yticks = -1:0.5:4,
        yminorticks = IntervalsBetween(5),
        yminorticksvisible = true
    )

    c1 = heatmap!(ax, x, y, data, colormap = colormap)

    if boundaries
        # plot boundaries between deformation regimes
        contour!(ax, x, y, mainDef, color = :red, linewidth = 2, linestyle = :solid, levels = n_components - 1)
    end

    contour!(ax, x, y, data; color = :black, levels = -20:1:-2, labels = true, labelsize = 25, labelfont = :bold, labelcolor = :black)

    Colorbar(fig[1, 2], c1, label = label, labelsize = fontsize)

    if !isnothing(filename)
        save(filename, fig)
    else
        display(fig)
    end

    return fig
end

"""
    fig, ax1, ax2, P_MPa, Tau_II_MPa, t_vec =  PlotPressureStressTime_0D(x;    args=(T=1000.0, P=0.0, d=1e-3, f=1.0),  εII::Union{Number, AbstractVector}, εvol::Union{Number, AbstractVector},
                                            τ0=0,
                                            P0=0,
                                            Time=(1e0,1e8), nt=100,
                                            t_vec=nothing,
                                            linestyle=:solid, linewidth=1, color=nothing, label=nothing, title="",
                                            fig=nothing, filename=nothing, res=(1200, 1200), legendsize=15, labelsize=35)


Creates a plot of Pressure and stress vs. time
"""
function PlotPressureStressTime_0D(
        x;
        args = (T = 1000.0, P = 0.0, d = 1.0e-3, f = 1.0),
        εII = 1.0e-15,
        εvol = -1.0e-18,
        τ0 = 0, P0 = 0, τ_scale = 1.0e6,
        Time = (1.0e0, 1.0e8), nt = 100,
        t_vec = nothing, t_scale = 3600 * 24 * 365.25 * 1.0e6,
        verbose = false,
        linestyle = :solid,
        linewidth = 1,
        color = nothing,
        label = nothing,
        title = "",
        fig = nothing,
        filename = nothing,
        res = (1200, 1200),
        legendsize = 15,
        labelsize = 35,
    )
    n = 1
    if isa(x, Tuple)
        n = length(x)
    end

    if isnothing(fig)
        fig = Figure(; fontsize = 25, size = res)
    end
    if τ_scale == 1.0
        ylabel_str = "Deviatoric stress"
    else
        ylabel_str = L"Deviatoric stress $\tau_{II}$ [MPa]"
    end
    if t_scale == 1.0
        xlabel_str = "Time"
    else
        xlabel_str = "Time [Myr]"
    end

    ax1 = Axis(
        fig[1, 1];
        ylabel = ylabel_str,
        xlabel = xlabel_str,
        xlabelsize = labelsize,
        ylabelsize = labelsize,
        title = title,
    )

    if τ_scale == 1.0
        ylabel_str = "Pressure"
    else
        ylabel_str = "Pressure [MPa]"
    end

    ax2 = Axis(
        fig[2, 1];
        ylabel = ylabel_str,
        xlabel = xlabel_str,
        xlabelsize = labelsize,
        ylabelsize = labelsize,
        title = title,
    )

    t_vec = []
    Tau_II_MPa = []
    P_MPa = []
    for i in 1:n
        if isa(x, Tuple)
            p = x[i]
        else
            p = x
        end
        if isa(args, Tuple)
            args_in = args[i]
        else
            args_in = args
        end

        # Compute
        t_vec, p_vec, τ_vec = time_p_τII_0D(p, εII, εvol, args; t = Time, nt = nt, verbose = verbose)
        Tau_II_MPa = τ_vec / τ_scale
        P_MPa = p_vec / τ_scale

        # Retrieve plot arguments (label, color etc.)
        plot_args = ObtainPlotArgs(i, p, args_in, linewidth, linestyle, color, label)

        # Create plot:
        li_1 = lines!(ax1, t_vec / t_scale, Tau_II_MPa)
        li_2 = lines!(ax2, t_vec / t_scale, P_MPa)

        # Customize plot:
        customize_plot!(li_1, plot_args)
        customize_plot!(li_2, plot_args)

    end

    axislegend(ax2; labelsize = legendsize)

    if !isnothing(filename)
        save(filename, fig)
    else
        display(fig)
    end

    return fig, ax1, ax2, P_MPa, Tau_II_MPa, t_vec
end


"""
    fig, ax = PlotDiffusionCoefPlotDiffusionCoefArrhenius(x::Union{Tuple{Vararg{AbstractChemicalDiffusion}}, NTuple{N, AbstractChemicalDiffusion} where N, AbstractChemicalDiffusion};
                                P=1u"GPa", fO2=1NoUnits, log_type=:log10, linestyle=:solid, linewidth=1, color=nothing, label=nothing,
                                title="", fig=nothing, filename=nothing, res=(1200, 1200), legend=true, legendsize=15, position=:rt,
                                labelsize=35, xlims=(nothing, nothing), ylims=(nothing, nothing))

Creates a plot of log(D) versus 10^4/T for one or a tuple of `ChemicalDiffusionData` structures.

# Optional parameters
- `P`: Pressure (default: 1 GPa)
- `fO2`: Oxygen fugacity (default: 1 NoUnits)
- `log_type`: Logarithm type for `D` (default: :log10, options: :log10, :ln)
- `linestyle`: Line style for the plot (default: :solid)
- `linewidth`: Line width for the plot (default: 1)
- `color`: Line color for the plot (default: nothing)
- `label`: Label for the plot (default: nothing)
- `title`: Title for the plot (default: "")
- `fig`: Existing figure to plot on (default: nothing)
- `filename`: Filename to save the plot (default: nothing)
- `res`: Resolution of the plot (default: (1200, 1200))
- `legend`: Whether to display the legend (default: true)
- `legendsize`: Size of the legend text (default: 15)
- `position`: Position of the legend (default: :rt)
- `labelsize`: Size of the axis labels (default: 35)
- `xlims`: Limits for the x-axis (default: (nothing, nothing))
- `ylims`: Limits for the y-axis (default: (nothing, nothing))
- `ticklabelsize`: Size of the tick labels (default: 35)
- `xlims2`: Limits for the x-axis of the second axis (temperature) (default: (nothing, nothing))

# Example

```julia
using GeoParams
using GLMakie

# obtain diffusion data
Fe_Grt = Garnet.Grt_Fe_Chakraborty1992
Fe_Grt = SetChemicalDiffusion(Fe_Grt)
Mg_Grt = Garnet.Grt_Mg_Chakraborty1992
Mg_Grt = SetChemicalDiffusion(Mg_Grt)
Mn_Grt = Garnet.Grt_Mn_Chakraborty1992
Mn_Grt = SetChemicalDiffusion(Mn_Grt)

fig, ax = PlotDiffusionCoefArrhenius((Fe_Grt, Mg_Grt, Mn_Grt), P= 1u"GPa", linewidth=3)
```
"""
function PlotDiffusionCoefArrhenius(
        x::Union{Tuple{Vararg{DiffusionData}}, NTuple{N, DiffusionData} where {N}, DiffusionData};
        P::Quantity = 1u"GPa",
        fO2 = 1NoUnits,
        X = 0NoUnits,
        log_type::Symbol = :log10,
        linestyle::Symbol = :solid,
        linewidth = 1,
        color = nothing,
        label = nothing,
        title = "",
        fig = nothing,
        filename = nothing,
        res = (1200, 1200),
        legend = true,
        legendsize = 25,
        position = :rt,
        labelsize = 35,
        xlims = (nothing, nothing),
        ylims = (nothing, nothing),
        ticklabelsize = 35,
        xlims2 = (nothing, nothing),
    )

    if isa(x, AbstractChemicalDiffusion)
        x = (x,)
    end

    n = length(x)
    name = ptr2string.(([x[i].Name for i in 1:n]))
    T_min = Value.([x[i].T_range_min for i in 1:n])
    T_max = Value.([x[i].T_range_max for i in 1:n])

    D = zeros(n * 2)
    for i in 1:n
        D[i] = compute_D(x[i], T = T_min[i], P = P, fO2 = fO2, X = X) |> upreferred |> ustrip
        D[i + n] = compute_D(x[i], T = T_max[i], P = P, fO2 = fO2, X = X) |> upreferred |> ustrip
    end

    T_inv_min = 1.0e4 ./ T_min .|> upreferred .|> ustrip
    T_inv_max = 1.0e4 ./ T_max .|> upreferred .|> ustrip

    if isnothing(fig)
        fig = Figure(; fontsize = 25, size = res)
    end

    if log_type == :log10
        D_log = log10.(D)
    elseif log_type == :ln
        D_log = log.(D)
    end

    ax1 = Axis(
        fig[1, 1];
        xlabel = "10⁴/T [K⁻¹]",
        xlabelsize = labelsize,
        ylabelsize = labelsize,
        xticklabelsize = ticklabelsize,
        yticklabelsize = ticklabelsize,
        title = title,
    )

    if log_type == :log10
        ax1.ylabel = "log₁₀ D (D in [m²/s])"
    elseif log_type == :ln
        ax1.ylabel = "ln D (D in [m²/s])"
    end

    if isnothing(label)
        label = name
    end

    for i in 1:n
        li = lines!(ax1, [T_inv_min[i], T_inv_max[i]], [D_log[i], D_log[i + n]])

        if isa(label, Array) || isa(label, Tuple)
            label_i = label[i]
        else
            label_i = label
        end

        plot_args = ObtainPlotArgs(i, "", "", linewidth, linestyle, color, label)

        args = (linewidth = plot_args.linewidth, linestyle = plot_args.linestyle, color = plot_args.color, label = label_i)

        customize_plot!(li, args)
    end

    if legend
        axislegend(ax1; labelsize = legendsize, position = position)
    end

    ax2 = Axis(
        fig[1, 1];
        xlabel = "Temperature [°C]",
        ylabel = "",
        xlabelsize = labelsize,
        ylabelsize = labelsize,
        title = "",
        xaxisposition = :top,
        xreversed = true,
        xticklabelsize = ticklabelsize,
        yticklabelsize = ticklabelsize,
    )

    hidexdecorations!(ax2; grid = true, label = false, ticklabels = false, ticks = false)

    T_min = uconvert.(u"°C", T_min) |> ustrip
    T_max = uconvert.(u"°C", T_max) |> ustrip

    for i in 1:n
        lines!(ax2, [T_min[i], T_max[i]], [D_log[i], D_log[i + n]], label = "", color = RGBAf(0, 0, 0, 0))
    end

    linkyaxes!(ax1, ax2)

    xlims!(ax1, xlims...)
    xlims!(ax2, xlims2...)
    ylims!(ax1, ylims...)
    ylims!(ax2, ylims...)

    if !isnothing(filename)
        save(filename, fig)
    else
        display(fig)
    end

    return fig, ax1
end
