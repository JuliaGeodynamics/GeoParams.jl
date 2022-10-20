"""
    This provides a few plotting routines, for example, for CreepLaws
"""

using Unitful
using Parameters
using ..Units
using ..MaterialParameters
using ..MeltingParam
using .GLMakie

using GeoParams: AbstractMaterialParam, AbstractMaterialParamsStruct
using .MaterialParameters.ConstitutiveRelationships
using .MaterialParameters.HeatCapacity: AbstractHeatCapacity, compute_heatcapacity
using .MaterialParameters.Conductivity: AbstractConductivity, compute_conductivity
using .MeltingParam: AbstractMeltingParam, compute_meltfraction

#Makie.inline!(true)

export PlotStrainrateStress,
    PlotStressStrainrate,
    PlotStrainrateViscosity,
    PlotStressViscosity,
    PlotHeatCapacity,
    PlotConductivity,
    PlotMeltFraction,
    PlotPhaseDiagram,
    Plot_TAS_diagram,
    Plot_ZirconAge_PDF,
    PlotStressTime_0D,
    PlotPressureStressTime_0D

"""
    fig, ax, εII,τII = PlotStrainrateStress(x; Strainrate=(1e-18,1e-12), args =(T=1000.0, P=0.0, d=1e-3, f=1.0), 
                                            linestyle=:solid, linewidth=1, color=nothing, label=nothing, title="", 
                                            fig=nothing, filename=nothing, res=(1200, 1200), legendsize=15, labelsize=35)
                                            
Plots deviatoric stress versus deviatoric strain rate for a single or multiple creeplaws 
    Note: if you want to create plots you need to install and load the `GLMakie.jl` package in julia.


# Example

First, we retrieve the data for anorthite creeplaws
```julia-repl
julia> pp  = SetDiffusionCreep("Dry Anorthite | Rybacki et al. (2006)");
julia> pp1 = SetDislocationCreep("Dry Anorthite | Rybacki et al. (2006)");
```
Next you can define each of the creeplaws inidvidually, plus a combined diffusion & dislocation creep law:
```julia-repl
julia> v   = (pp,pp1,(pp,pp1));   
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
    Strainrate=(1e-18, 1e-12),
    args=(T=1000.0, P=0.0, d=1e-3, f=1.0),
    linestyle=:solid,
    linewidth=1,
    color=nothing,
    label=nothing,
    title="",
    fig=nothing,
    filename=nothing,
    res=(1200, 1200),
    legendsize=15,
    labelsize=35,
)
    n = 1
    if isa(x, Tuple)
        n = length(x)
    end
    if isnothing(fig)
        fig = Figure(; fontsize=25, resolution=res)
    end
    ax = Axis(
        fig[1, 1];
        yscale=log10,
        xscale=log10,
        xlabel=L"Deviatoric strain rate $\dot{ε}_{II}$ [1/s]",
        ylabel=L"Deviatoric stress $\tau_{II}$ [MPa]",
        xlabelsize=labelsize,
        ylabelsize=labelsize,
        title=title,
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
                    stop=ustrip(log10(Strainrate[2])),
                    length=101,
                )
            )
        Tau_II = zeros(size(Eps_II))

        # Compute stress
        #compute_τII!(Tau_II, p, Eps_II, args_in)
        for j in eachindex(Tau_II)
            Tau_II[j] = compute_τII(p, Eps_II[j], args)
        end

        Tau_II_MPa = Tau_II ./ 1e6

        # Retrieve plot arguments (label, color etc.)
        plot_args = ObtainPlotArgs(i, p, args_in, linewidth, linestyle, color, label)

        # Create plot:
        li = lines!(Eps_II, Tau_II_MPa)    # plot line

        # Customize line:
        customize_plot!(li, plot_args)
    end
    axislegend(ax; labelsize=legendsize)

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
            Name = "";
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
    args = (linewidth=linewidth_in, linestyle=linestyle_in, label=label, color=color_in)

    return args
end

# Internal fucntion that customizes the plot
function customize_plot!(li, args)

    # Customize line:
    li.label = args.label
    li.linewidth = args.linewidth
    li.linestyle = args.linestyle
    if !isnothing(args.color)
        li.color = args.color
    end
end

"""
    fig,ax,τII,εII =  PlotStressStrainrate(x; args=(T=1000.0, P=0.0, d=1e-3, f=1.0), Stress=(1e0,1e8), plt=nothing)

Same as `PlotStrainrateStress` but with stress (in MPa) versus strainrate (in 1/s) instead.

"""
function PlotStressStrainrate(
    x;
    args=(T=1000.0, P=0.0, d=1e-3, f=1.0),
    Stress=(1e0, 1e8),
    linestyle=:solid,
    linewidth=1,
    color=nothing,
    label=nothing,
    title="",
    fig=nothing,
    filename=nothing,
    res=(1200, 1200),
    legendsize=15,
    labelsize=35,
)
    n = 1
    if isa(x, Tuple)
        n = length(x)
    end

    if isnothing(fig)
        fig = Figure(; fontsize=25, resolution=res)
    end
    ax = Axis(
        fig[1, 1];
        yscale=log10,
        xscale=log10,
        xlabel=L"Deviatoric stress $\tau_{II}$ [MPa]",
        ylabel=L"Deviatoric strain rate $\dot{ε}_{II}$ [1/s]",
        xlabelsize=labelsize,
        ylabelsize=labelsize,
        title=title,
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
        Tau_II_MPa = range(ustrip(Stress[1]); stop=ustrip(Stress[2]), length=101)
        Tau_II = Tau_II_MPa .* 1e6
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

    axislegend(ax; labelsize=legendsize)

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
    args=(T=1000.0, P=0.0, d=1e-3, f=1.0),
    Strainrate=(1e-18, 1e-12),
    linestyle=:solid,
    linewidth=1,
    color=nothing,
    label=nothing,
    title="",
    fig=nothing,
    filename=nothing,
    res=(1200, 1200),
    legendsize=15,
    labelsize=35,
)
    n = 1
    if isa(x, Tuple)
        n = length(x)
    end

    if isnothing(fig)
        fig = Figure(; fontsize=25, resolution=res)
    end
    ax = Axis(
        fig[1, 1];
        yscale=log10,
        xscale=log10,
        xlabel=L"Deviatoric strain rate $\dot{ε}_{II}$ [1/s]",
        ylabel=L"Effective viscosity $\eta$ [Pa s]",
        xlabelsize=labelsize,
        ylabelsize=labelsize,
        title=title,
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
                    stop=ustrip(log10(Strainrate[2])),
                    length=101,
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
            η[1] *= (1.0 - 1e-5)
        end

        Tau_II_MPa = Tau_II ./ 1e6

        # Retrieve plot arguments (label, color etc.)
        plot_args = ObtainPlotArgs(i, p, args_in, linewidth, linestyle, color, label)

        # Create plot:
        li = lines!(Eps_II, η)    # plot line

        # Customize line:
        customize_plot!(li, plot_args)
    end

    axislegend(ax; labelsize=legendsize)

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
    args=(T=1000.0, P=0.0, d=1e-3, f=1.0),
    Stress=(1e0, 1e8),
    linestyle=:solid,
    linewidth=1,
    color=nothing,
    label=nothing,
    title="",
    fig=nothing,
    filename=nothing,
    res=(1200, 1200),
    legendsize=15,
    labelsize=35,
)
    n = 1
    if isa(x, Tuple)
        n = length(x)
    end

    if isnothing(fig)
        fig = Figure(; fontsize=25, resolution=res)
    end
    ax = Axis(
        fig[1, 1];
        yscale=log10,
        xscale=log10,
        xlabel=L"Deviatoric stress $\tau_{II}$ [MPa]",
        ylabel=L"Effective viscosity $\eta$ [Pa s]",
        xlabelsize=labelsize,
        ylabelsize=labelsize,
        title=title,
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
        Tau_II_MPa = range(ustrip(Stress[1]); stop=ustrip(Stress[2]), length=101)
        Tau_II = Tau_II_MPa .* 1e6
        Eps_II = zeros(size(Tau_II))

        #compute_εII!(Eps_II, p, Tau_II, args_in)       # Compute strainrate
        for j in eachindex(Tau_II)
            Eps_II[j] = compute_εII(p, Tau_II[j], args)
        end

        η = Tau_II ./ (2 * Eps_II)                        # effective viscosity

        if maximum(η) ≈ minimum(η)
            # if all values are the same, plotting creates an error can occur, so we perturb the values a bit
            η[1] *= (1.0 - 1e-5)
        end

        # Retrieve plot arguments (label, color etc.)
        plot_args = ObtainPlotArgs(i, p, args_in, linewidth, linestyle, color, label)

        # Create plot:
        li = lines!(Tau_II_MPa, η)    # plot line

        # Customize plot:
        customize_plot!(li, plot_args)
    end

    axislegend(ax; labelsize=legendsize)

    if !isnothing(filename)
        save(filename, fig)
    else
        display(fig)
    end

    return fig, ax, Tau_II_MPa, η
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


"""
    fig,ax,τII,η =  PlotStressTime_0D(x;    args=(T=1000.0, P=0.0, d=1e-3, f=1.0),  
                                            εII::Union{Number, AbstractVector}, 
                                            τ0=0,
                                            Time=(1e0,1e8), nt=100,
                                            t_vec=nothing,
                                            linestyle=:solid, linewidth=1, color=nothing, label=nothing, title="", 
                                            fig=nothing, filename=nothing, res=(1200, 1200), legendsize=15, labelsize=35)


Creates a plot of stress vs. time, or stress vs. strain (in)
"""
function PlotStressTime_0D(
    x;
    args=(T=1000.0, P=0.0, d=1e-3, f=1.0),
    εII=1e-15,
    τ0=0,τ_scale=1e6,
    Time=(1e0, 1e8), nt=100,
    t_vec=nothing, t_scale=3600*24*365.25*1e6,
    verbose=false,
    linestyle=:solid,
    linewidth=1,
    color=nothing,
    label=nothing,
    title="",
    fig=nothing,
    filename=nothing,
    res=(1200, 1200),
    legendsize=15,
    labelsize=35,
)
    n = 1
    if isa(x, Tuple)
        n = length(x)
    end

    if isnothing(fig)
        fig = Figure(; fontsize=25, resolution=res)
    end
    ax = Axis(
        fig[1, 1];
        ylabel=L"Deviatoric stress $\tau_{II}$ [MPa]",
        xlabel="Time [Myrs]",
        xlabelsize=labelsize,
        ylabelsize=labelsize,
        title=title,
    )

    t_vec=[]
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
        t_vec, τ_vec = time_τII_0D(p, εII, args; t=Time, nt=nt, verbose=verbose)
        Tau_II_MPa = τ_vec/τ_scale;

        # Retrieve plot arguments (label, color etc.)
        plot_args = ObtainPlotArgs(i, p, args_in, linewidth, linestyle, color, label)

        # Create plot:
        li = lines!(t_vec/t_scale, Tau_II_MPa)    # plot line

        
        # Customize plot:
        customize_plot!(li, plot_args)
    end

    axislegend(ax; labelsize=legendsize)

    if !isnothing(filename)
        save(filename, fig)
    else
        display(fig)
    end

    return fig, ax, Tau_II_MPa, t_vec
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
    args=(T=1000.0, P=0.0, d=1e-3, f=1.0),
    εII=1e-15,
    εvol=-1e-18,
    τ0=0,P0=0, τ_scale=1e6,
    Time=(1e0, 1e8), nt=100,
    t_vec=nothing, t_scale=3600*24*365.25*1e6,
    verbose=false,
    linestyle=:solid,
    linewidth=1,
    color=nothing,
    label=nothing,
    title="",
    fig=nothing,
    filename=nothing,
    res=(1200, 1200),
    legendsize=15,
    labelsize=35,
)
    n = 1
    if isa(x, Tuple)
        n = length(x)
    end

    if isnothing(fig)
        fig = Figure(; fontsize=25, resolution=res)
    end
    if τ_scale == 1.0
        ylabel_str = "Deviatoric stress";
    else
        ylabel_str = L"Deviatoric stress $\tau_{II}$ [MPa]";
    end
    if t_scale == 1.0
        xlabel_str = "Time";
    else
        xlabel_str = "Time [Myrs]";
    end
    
    ax1 = Axis(
        fig[1, 1];
        ylabel=ylabel_str,
        xlabel=xlabel_str,
        xlabelsize=labelsize,
        ylabelsize=labelsize,
        title=title,
    )

    if τ_scale == 1.0
        ylabel_str = "Pressure";
    else
        ylabel_str = "Pressure [MPa]";
    end

    ax2 = Axis(
        fig[2, 1];
        ylabel=ylabel_str,
        xlabel=xlabel_str,
        xlabelsize=labelsize,
        ylabelsize=labelsize,
        title=title,
    )

    t_vec=[]
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
        t_vec, p_vec, τ_vec = time_p_τII_0D(p, εII, εvol, args; t=Time, nt=nt, verbose=verbose)
        Tau_II_MPa = τ_vec/τ_scale;
        P_MPa      = p_vec/τ_scale;
        
        # Retrieve plot arguments (label, color etc.)
        plot_args = ObtainPlotArgs(i, p, args_in, linewidth, linestyle, color, label)

        # Create plot:
        li_1 = lines!(ax1, t_vec/t_scale, Tau_II_MPa)    
        li_2 = lines!(ax2, t_vec/t_scale, P_MPa)    
        
        # Customize plot:
        customize_plot!(li_1, plot_args)
        customize_plot!(li_2, plot_args)
        
    end

    axislegend(ax2; labelsize=legendsize)
    
    if !isnothing(filename)
        save(filename, fig)
    else
        display(fig)
    end

    return fig, ax1, ax2, P_MPa, Tau_II_MPa, t_vec
end