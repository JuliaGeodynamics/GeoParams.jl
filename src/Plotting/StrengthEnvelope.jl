using Parameters
using SpecialFunctions: erfc
using .GLMakie

export StrengthEnvelope, ConstantTemp, LinearTemp, HalfspaceCoolingTemp, GP_Compute_ThermalStructure, LithPres

abstract type AbstractThermalStructure end  

"""
    ConstantTemp(T=1000)
    
Sets a constant temperature inside the box
Parameters
===
- T : the value
"""
@with_kw_noshow mutable struct ConstantTemp <: AbstractThermalStructure
    T = 1000
end

function GP_Compute_ThermalStructure(Z, s::ConstantTemp)
    Temp  = zeros(Float64, length(Z))
    Temp .= s.T
    return Temp
end


"""
    LinearTemp(Ttop=0, Tbot=1000)
    
Set a linear temperature structure from top to bottom
Parameters
===
- Ttop : the value @ the top
- Tbot : the value @ the bottom
"""
@with_kw_noshow mutable struct LinearTemp <: AbstractThermalStructure
    Ttop = 0
    Tbot = 1350
end

function GP_Compute_ThermalStructure(Z, s::LinearTemp)
    @unpack Ttop, Tbot  = s

    Temp  = zeros(Float64, length(Z))

    dz    = Z[end]-Z[1];
    dT    = Tbot - Ttop

    Temp .= abs.(Z./dz).*dT .+ Ttop
    return Temp
end

"""
    HalfspaceCoolingTemp(Tsurface=0, Tmantle=1350, Age=60, Adiabat=0)
    
Sets a halfspace temperature structure in plate
Parameters
========
- Tsurface : surface temperature [C]
- Tmantle : mantle temperature [C]
- Age : Thermal Age of plate [Myrs]
- Adiabat : Mantle Adiabat [K/km]
"""
@with_kw_noshow mutable struct HalfspaceCoolingTemp <: AbstractThermalStructure
    Tsurface = 0        # top T
    Tmantle  = 1350     # bottom T
    Age      = 60       # thermal age of plate [in Myrs]
    Adiabat  = 0        # Adiabatic gradient in K/km
    kappa    = 1e-6     # thermal diffusivity
end

function GP_Compute_ThermalStructure(Z, s::HalfspaceCoolingTemp)
    @unpack Tsurface, Tmantle, Age, Adiabat, kappa  = s

    Temp  = zeros(Float64, length(Z))

    MantleAdiabaticT    =   Tmantle .+ Adiabat*abs.(Z);   # Adiabatic temperature of mantle
    
    for i in eachindex(Temp)
        Temp[i] =   (Tsurface .- Tmantle)*erfc((abs.(Z[i]))./(2*sqrt(kappa*Age))) + MantleAdiabaticT[i];
    end
    return Temp
end

"""
    LithPres(MatParam, Phases, ρ, T, nz, dz, g)

Iteratively solves a 1D lithostatic pressure profile (compatible with temperature- and pressure-dependent densities)
Parameters
========
- ρ:  density vector for initial guess(can be zeros)
- T:  temperature vector
- dz: grid spacing
- g:  gravitational accelaration
"""
function LithPres(MatParam, Phases, ρ, T, dz, g)
    nz   = length(T)
    P    = zeros(Float64, nz)

    tol  = 1e-6
    n    = 0
    stop = 1
    iter = 0
    
    while stop > tol && iter < 100
        # store current norm
        nprev = n
        iter += 1

        # update densities
        args  = (T=T, P=P)
        compute_density!(ρ, MatParam, Phases, args)

        # update pressure
        for i = 2 : nz
            P[i] = P[i-1] + ρ[i-1] * g * dz
        end

        # compute new norm
        n    = (sum(P.^2)).^0.5

        # evaluate stopping criterium
        stop = abs((n-nprev) / (n+nprev))
    end

    #println("Finding lithostatic pressure took $iter iterations.")
    return P
end

function solveStress(MatParam, Phases, ε, P, T)
    # solve for stress
    nz = length(T)
    τ  = zeros(Float64, nz)
    for i = 1 : nz
        Pres = P[i]
        Temp = T[i]
        args = (T=Temp, P=Pres)
        Mat  = MatParam[Phases[i]]
        τ[i] = compute_τII(Mat.CreepLaws[1], ε, args)

        F    = compute_yieldfunction(Mat.Plasticity[1], P=Pres, τII=τ[i])
        if F > 0
            ϕ = Mat.Plasticity[1].ϕ.val
            c = Mat.Plasticity[1].C.val
            τ[i] = Pres * sind(ϕ) + cosd(ϕ) * c
        end
    end

    return τ
end

function redimensionalize(x, nz, param_dim::Unitful.FreeUnits, CharDim)
    x_val = zeros(Float64, nz)
    x_dim = dimensionalize(x, param_dim, CharDim)
    for i = 1 : nz
        x_val[i] = x_dim[i].val
    end
    return x_val
end

"""
    StrengthEnvelopeSliders(MatParam, Thickness, TempType=LinearTemp()::AbstractThermalStructure)

Creates a GUI that plots a 1D strength envelope. In the GUI, temperature profile and strain rate can be adjusted. The Drucker-Prager plasticity uses lithostatic pressure.

Parameters:
- MatParam:  a tuple of materials (including the following properties: Phase, Density, CreepLaws, Plasticity)
- Thickness: a vector listing the thicknesses of the respective layers (should carry units)
- TempType:  the type of temperature profile (LinearTemp=default, HalfspaceCoolingTemp, ConstantTemp)

# Example:
```julia-repl
julia> using GLMakie
julia> MatParam = (SetMaterialParams(Name="UC", Phase=1, Density=ConstantDensity(ρ=2700kg/m^3), CreepLaws = SetDislocationCreep("Wet Quartzite | Ueda et al. (2008)"), Plasticity = DruckerPrager(ϕ=30.0, C=10MPa)),
                   SetMaterialParams(Name="MC", Phase=2, Density=Density=ConstantDensity(ρ=2900kg/m^3), CreepLaws = SetDislocationCreep("Plagioclase An75 | Ji and Zhao (1993)"), Plasticity = DruckerPrager(ϕ=20.0, C=10MPa)),
                   SetMaterialParams(Name="LC", Phase=3, Density=PT_Density(ρ0=2900kg/m^3, α=3e-5/K, β=1e-10/Pa), CreepLaws = SetDislocationCreep("Maryland strong diabase | Mackwell et al. (1998)"), Plasticity = DruckerPrager(ϕ=30.0, C=10MPa)));
julia> Thickness = [15,10,15]*km;

julia> StrengthEnvelope(MatParam, Thickness, LinearTemp())
```
"""
function StrengthEnvelope(MatParam::NTuple{N, AbstractMaterialParamsStruct}, Thickness::Vector{U}, TempType::AbstractThermalStructure=LinearTemp(), mode::String="normal") where {N, U}

    # hardcoded input
    nz        = 101
    g         = 9.81m/s^2

    # nondimensionalize
    CharDim   = GEO_units(length=10km, temperature=1000C, stress=10MPa, viscosity=1e20Pas)
    MatParam  = nondimensionalize(MatParam, CharDim)
    Thickness = nondimensionalize(Thickness, CharDim)
    g         = nondimensionalize(9.81m/s^2, CharDim)
    
    # derived properties
    nLayer    = length(MatParam)
    Inter     = cumsum(Thickness)
    dz        = Inter[end] / (nz-1)
    z         = collect(0:dz:Inter[end])

    # distribute phases
    Phases    = ones(Int64, nz) * MatParam[1].Phase
    for i = 1 : nLayer - 1
        Phases[z .> Inter[i]] .= MatParam[i+1].Phase
    end

    # check temperature structure
    if typeof(TempType) == LinearTemp
        Ttype     = 1
        title     = "1D Strength Envelope (Linear T-profile, Ttop = 0C)"
        Ttop      = 0C
        Ttop      = nondimensionalize(Ttop, CharDim)
        # Tbot controlled by slider
    elseif typeof(TempType) == HalfspaceCoolingTemp
        Ttype     = 2
        title     = "1D Strength Envelope (Halfspace-cooling T-profile, Ttop = 0C, Tmantle = 1350C)"
        Ttop      = 0C
        Tbot      = 1350C
        Adiabat   = 0K/km
        kappa     = 1e-6m^2/s
        Ttop      = nondimensionalize(Ttop, CharDim)
        Tbot      = nondimensionalize(Tbot, CharDim)
        Adiabat   = nondimensionalize(Adiabat, CharDim)
        kappa     = nondimensionalize(kappa, CharDim)
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
        lsgrid = SliderGrid(fig[2, :],
            (label = "Temp [C]", range = 0:2:1400, startvalue = 500),
            (label = "log10(ε) [s-1]", range=-18:0.01:-10, startvalue = -15)
        )
    end

    # create listener to sliders
    T_slider     = lsgrid.sliders[1].value
    

    exp_slider   = lsgrid.sliders[2].value
    ε_dim        = @lift(10^$exp_slider/s)
    ε            = @lift(nondimensionalize($ε_dim, CharDim))

    # build temperature structure

    if Ttype == 1
        Tbot_dim  = @lift($T_slider*C)
        Tbot      = @lift(nondimensionalize($Tbot_dim, CharDim))
        T         = @lift(GP_Compute_ThermalStructure(z, LinearTemp(Ttop, $Tbot)))
    elseif Ttype == 2
        Age_dim   = @lift($T_slider*Myrs)
        Age       = @lift(nondimensionalize($Age_dim, CharDim))
        T         = @lift(GP_Compute_ThermalStructure(z, HalfspaceCoolingTemp(Ttop, Tbot, $Age, Adiabat, kappa)))
    elseif Ttype == 3
        Tcon_dim  = @lift($T_slider*C)
        Tcon      = @lift(nondimensionalize($Tcon_dim, CharDim))
        T         = @lift(GP_Compute_ThermalStructure(z, ConstantTemp($Tcon)))
    end

    # pressure and density
    ρ         = Observable(zeros(Float64, nz))
    P         = @lift(LithPres(MatParam, Phases, $ρ, $T, dz, g))

    # solve for stress
    τ         = @lift(solveStress(MatParam, Phases, $ε, $P, $T))
    
    # redimensionalize
    z_plot    = dimensionalize(GeoUnit(z, km, false), CharDim).val
    τ_plot    = @lift(redimensionalize($τ, nz, MPa, CharDim))
    T_plot    = @lift(redimensionalize($T, nz, C, CharDim))

    # plotting
    lines!(ax1, τ_plot, z_plot)
    lines!(ax2, T_plot, z_plot)
    ylims!(ax1, [maximum(z_plot), 0])
    ylims!(ax2, [maximum(z_plot), 0])
    xlims!(ax1, low = 0)
    xlims!(ax2, (0, 1350))

    # show figure
    if mode == "normal"
        return fig
    elseif mode == "test"
        return τ_plot
    end
end