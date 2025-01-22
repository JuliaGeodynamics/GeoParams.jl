# This function is implemented in this file
import GeoParams: StrengthEnvelopePlot

# Make all `export`ed names from GeoParams.jl available
using GeoParams
# We also need the following un-`export`ed names
using GeoParams: AbstractTempStruct

using Parameters

function StrengthEnvelopePlot(MatParam::NTuple{N, AbstractMaterialParamsStruct}, Thickness::Vector{U}; TempType::AbstractTempStruct = LinTemp(), nz::Int64 = 101) where {N, U}

    # check temperature structure
    if typeof(TempType) == LinTemp
        Ttype = 1
        title = "1D Strength Envelope (Linear T-profile, Ttop = 0C)"
        Ttop = 0C
        # Tbot controlled by slider
    elseif typeof(TempType) == HalfspaceCoolTemp
        Ttype = 2
        title = "1D Strength Envelope (Halfspace-cooling T-profile, Ttop = 0C, Tmantle = 1350C)"
        Ttop = 0C
        Tbot = 1350C
        Adiabat = 0K / km
        kappa = 1.0e-6m^2 / s
        # Age controlled by slider
    elseif typeof(TempType) == ConstTemp
        Ttype = 3
        title = "1D Strength Envelope (Constant T)"
        # Temp controlled by slider
    end

    # build Figure
    fig = Figure(size = (1200, 900))
    ax1 = fig[1, 1:3] = Axis(
        fig,
        # title
        title = title,
        # x-axis
        xlabel = "Maximum Strength [MPa]",
        # y-axis
        ylabel = "Depth [km]",
    )
    ax2 = fig[1, 4] = Axis(
        fig,
        # title
        title = "Temperature",
        # x-axis
        xlabel = "Temperature [C]",
        # y-axis
        ylabel = "Depth [km]"
    )

    # make sliders
    if Ttype == 1
        lsgrid = SliderGrid(
            fig[2, :],
            (label = "Tbot [C]", range = 500:1:1400, startvalue = 800),
            (label = "log10(ε) [s-1]", range = -18:0.01:-10, startvalue = -15)
        )
    elseif Ttype == 2
        lsgrid = SliderGrid(
            fig[2, :],
            (label = "Plate age [Myr]", range = 1:300, startvalue = 100),
            (label = "log10(ε) [s-1]", range = -18:0.01:-10, startvalue = -15)
        )
    elseif Ttype == 3
        lsgrid = SliderGrid(
            fig[2, :], # show figure
            (label = "Temp [C]", range = 0:2:1400, startvalue = 500),
            (label = "log10(ε) [s-1]", range = -18:0.01:-10, startvalue = -15)
        )
    end

    # create listener to sliders
    T_slider = lsgrid.sliders[1].value

    exp_slider = lsgrid.sliders[2].value
    ε = @lift(10^$exp_slider / s)

    # build temperature structure
    if Ttype == 1
        Tbot = @lift($T_slider * C)
        Temp = @lift(LinTemp(Ttop, $Tbot))
    elseif Ttype == 2
        Age = @lift($T_slider * Myrs)
        Temp = @lift(HalfspaceCoolTemp(Ttop, Tbot, $Age, Adiabat, kappa))
    elseif Ttype == 3
        Tcon = @lift($T_slider * C)
        Temp = @lift(ConstTemp($Tcon))
    end

    # get results from computational routine
    res = @lift(StrengthEnvelopeComp(MatParam, Thickness, $Temp, $ε, nz))

    # extract components for plotting
    z_plot = @lift(extractFromResult($res, 1, nz))
    τ_plot = @lift(extractFromResult($res, 2, nz))
    T_plot = @lift(extractFromResult($res, 3, nz))

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
