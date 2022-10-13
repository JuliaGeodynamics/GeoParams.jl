using Parameters
using .GLMakie

export StrengthEnvelopePlot

"""
    StrengthEnvelopePlot(MatParam, Thickness; TempType, nz)

Creates a GUI that plots a 1D strength envelope. In the GUI, temperature profile and strain rate can be adjusted. The Drucker-Prager plasticity uses lithostatic pressure.

Parameters:
- MatParam:  a tuple of materials (including the following properties: Phase, Density, CreepLaws, Plasticity)
- Thickness: a vector listing the thicknesses of the respective layers (should carry units)
- TempType:  the type of temperature profile (LinearTemp=default, HalfspaceCoolingTemp, ConstantTemp)
- nz:        optional argument controlling the number of points along the profile (default = 101)

# Example:
```julia-repl
julia> using GLMakie
julia> MatParam = (SetMaterialParams(Name="UC", Phase=1, Density=ConstantDensity(ρ=2700kg/m^3), CreepLaws = SetDislocationCreep("Wet Quartzite | Ueda et al. (2008)"), Plasticity = DruckerPrager(ϕ=30.0, C=10MPa)),
                   SetMaterialParams(Name="MC", Phase=2, Density=Density=ConstantDensity(ρ=2900kg/m^3), CreepLaws = SetDislocationCreep("Plagioclase An75 | Ji and Zhao (1993)"), Plasticity = DruckerPrager(ϕ=20.0, C=10MPa)),
                   SetMaterialParams(Name="LC", Phase=3, Density=PT_Density(ρ0=2900kg/m^3, α=3e-5/K, β=1e-10/Pa), CreepLaws = SetDislocationCreep("Maryland strong diabase | Mackwell et al. (1998)"), Plasticity = DruckerPrager(ϕ=30.0, C=10MPa)));
julia> Thickness = [15,10,15]*km;

julia> StrengthEnvelopePlot(MatParam, Thickness, LinearTemp())
```
"""
function StrengthEnvelopePlot(MatParam::NTuple{N, AbstractMaterialParamsStruct}, Thickness::Vector{U}; TempType::AbstractThermalStructure=LinearTemp(), nz::Int64=101) where {N, U}

    # check temperature structure
    if typeof(TempType) == LinearTemp
        Ttype     = 1
        title     = "1D Strength Envelope (Linear T-profile, Ttop = 0C)"
        Ttop      = 0C
        # Tbot controlled by slider
    elseif typeof(TempType) == HalfspaceCoolingTemp
        Ttype     = 2
        title     = "1D Strength Envelope (Halfspace-cooling T-profile, Ttop = 0C, Tmantle = 1350C)"
        Ttop      = 0C
        Tbot      = 1350C
        Adiabat   = 0K/km
        kappa     = 1e-6m^2/s
        # Age controlled by slider
    elseif typeof(TempType) == ConstantTemp
        Ttype     = 3
        title     = "1D Strength Envelope (Constant T)"
        # Temp controlled by slider
    end

    # build Figure
    fig = Figure(resolution = (1200, 900));
    ax1 = fig[1, 1:3] = Axis(fig,
        # title
        title  = title,
        # x-axis
        xlabel = "Maximum Strength [MPa]",
        # y-axis
        ylabel = "Depth [km]",
    )
    ax2 = fig[1, 4] = Axis(fig,
        # title
        title  = "Temperature",
        # x-axis
        xlabel = "Temperature [C]",
        # y-axis
        ylabel = "Depth [km]"
    )

    # make sliders
    if Ttype == 1
        lsgrid = SliderGrid(fig[2, :],
            (label = "Tbot [C]", range = 500:1:1400, startvalue = 800),
            (label = "log10(ε) [s-1]", range=-18:0.01:-10, startvalue = -15)
        )
    elseif Ttype == 2
        lsgrid = SliderGrid(fig[2, :],
            (label = "Plate age [Myr]", range = 1:300, startvalue = 100),
            (label = "log10(ε) [s-1]", range=-18:0.01:-10, startvalue = -15)
        )
    elseif Ttype == 3
        lsgrid = SliderGrid(fig[2, :],# show figure
            (label = "Temp [C]", range = 0:2:1400, startvalue = 500),
            (label = "log10(ε) [s-1]", range=-18:0.01:-10, startvalue = -15)
        )
    end

    # create listener to sliders
    T_slider     = lsgrid.sliders[1].value
    
    exp_slider   = lsgrid.sliders[2].value
    ε            = @lift(10^$exp_slider/s)

    # build temperature structure
    if Ttype == 1
        Tbot     = @lift($T_slider*C)
        Temp     = @lift(LinearTemp(Ttop, $Tbot))
    elseif Ttype == 2
        Age      = @lift($T_slider*Myrs)
        Temp     = @lift(HalfspaceCoolingTemp(Ttop, Tbot, $Age, Adiabat, kappa))
    elseif Ttype == 3
        Tcon     = @lift($T_slider*C)
        Temp     = @lift(ConstantTemp($Tcon))
    end

    # get results from computational routine
    res          = @lift(StrengthEnvelopeComp(MatParam, Thickness, $Temp, $ε, nz))
    
    # extract components for plotting
    z_plot       = @lift(extractFromResult($res, 1, nz))
    τ_plot       = @lift(extractFromResult($res, 2, nz))
    T_plot       = @lift(extractFromResult($res, 3, nz))

    # plotting
    lines!(ax1, τ_plot, z_plot)
    lines!(ax2, T_plot, z_plot)
    ylims!(ax1, [maximum(z_plot.val), 0])
    ylims!(ax2, [maximum(z_plot.val), 0])
    xlims!(ax1, low = 0)
    xlims!(ax2, (0, 1350))
    ax1.yreversed = true
    ax2.yreversed = true

    return fig
end
