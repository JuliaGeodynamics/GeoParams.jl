using Parameters

export StrengthEnvelope, StrengthEnvelopeSliders, ConstantTemp, LinearTemp, HalfspaceCoolingTemp, GP_Compute_ThermalStructure!, 
       GP_Compute_ThermalStructure, LithPres

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

function GP_Compute_ThermalStructure!(Temp, Z, s::ConstantTemp)
    Temp .= s.T
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

function GP_Compute_ThermalStructure!(Temp, Z, s::LinearTemp)
    @unpack Ttop, Tbot  = s

    dz    = Z[end]-Z[1];
    dT    = Tbot - Ttop

    Temp .= abs.(Z./dz).*dT .+ Ttop
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
    Tsurface = 0       # top T
    Tmantle = 1350     # bottom T
    Age  = 60          # thermal age of plate [in Myrs]
    Adiabat = 0        # Adiabatic gradient in K/km
end

function GP_Compute_ThermalStructure!(Temp, Z, s::HalfspaceCoolingTemp)
    @unpack Tsurface, Tmantle, Age, Adiabat  = s

    kappa       =   1e-6;
    SecYear     =   3600*24*365.25
    ThermalAge  =   Age*1e6*SecYear;

    MantleAdiabaticT    =   Tmantle .+ Adiabat*abs.(Z);   # Adiabatic temperature of mantle
    
    for i in eachindex(Temp)
        Temp[i] =   (Tsurface .- Tmantle)*erfc((abs.(Z[i])*1e3)./(2*sqrt(kappa*ThermalAge))) + MantleAdiabaticT[i];
    end
end

"""
    This creates a 1D strength envelope
"""
function StrengthEnvelope(MatParam, Thickness, ε, TempType::AbstractThermalStructure)

    # hardcoded input
    nz        = 101
    g         = 9.81m/s^2

    # nondimensionalize
    CharDim   = GEO_units(length=10km, temperature=1000C, stress=10MPa, viscosity=1e20Pas)
    MatParam  = nondimensionalize(MatParam, CharDim)
    Thickness = nondimensionalize(Thickness, CharDim)
    ε         = nondimensionalize(ε, CharDim)
    g         = nondimensionalize(9.81m/s^2, CharDim)
    
    # derived properties
    nLayer    = length(MatParam)
    Inter     = cumsum(Thickness)
    dz        = Inter[end] / (nz-1)
    z         = collect(0:dz:Inter[end])

    # build temperature structure
    T         = zeros(Float64, nz).*C
    T         = GP_Compute_ThermalStructure!(T, z, TempType)
    T         = nondimensionalize(T, CharDim)

    # distribute phases
    Phases    = ones(Int64, nz) * MatParam[1].Phase
    for i = 1 : nLayer - 1
        Phases[z .> Inter[i]] .= MatParam[i+1].Phase
    end

    # pressure and density
    ρ       = zeros(Float64, nz)
    P       = zeros(Float64, nz)
    LithPres!(ρ, P, MatParam, Phases, T, nz, dz, g)

    # solve for stress
    τ = zeros(Float64, nz)
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
    
    # redimensionalize
    z = dimensionalize(GeoUnit(z, km,     false), CharDim)
    τ = dimensionalize(GeoUnit(τ, MPa,    false), CharDim)
    ρ = dimensionalize(GeoUnit(ρ, kg/m^3, false), CharDim)
    P = dimensionalize(GeoUnit(P, MPa,    false), CharDim)
    T = dimensionalize(GeoUnit(T, C,    false), CharDim)

    # plotting
    fig = Figure()
    lines(fig[1,1:2], ustrip(Value(τ)), -ustrip(Value(z)); axis = (; xlabel="Maximum Deviatoric Stress [MPa]", ylabel="Depth [km]"))
    lines(fig[1,3], ustrip(Value(ρ)), -ustrip(Value(z)); axis = (; xlabel="Density [kg m-3]", ylabel="Depth [km]"))
    lines(fig[1,4], ustrip(Value(P)), -ustrip(Value(z)); axis = (; xlabel="Lithostatic Pressure [MPa]", ylabel="Depth [km]"))
    #lines(fig[1,4], ustrip(Value(T)), -ustrip(Value(z)); axis = (; xlabel="Temperature [C]", ylabel="Depth [km]"))
    fig
end

function LithPres!(ρ, P, MatParam, Phases, T, nz, dz, g)
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

    println("Finding lithostatic pressure took $iter iterations.")
end

function LithPres(MatParam, Phases, ρ, T, nz, dz, g)
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

function solveStress(MatParam, Phases, ε, P, T, nz)
    # solve for stress
    τ = zeros(Float64, nz)
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

function redimensionalize!(τ, ρ, P, T, CharDim)
    #τ = dimensionalize(GeoUnit(τ, MPa,    false), CharDim).val
    #ρ = dimensionalize(GeoUnit(ρ, kg/m^3, false), CharDim).val
    #P = dimensionalize(GeoUnit(P, MPa,    false), CharDim).val
    #T = dimensionalize(GeoUnit(T, C,      false), CharDim).val

    return τ
end

function StrengthEnvelopeSliders(MatParam, Thickness)

    # hardcoded input
    nz        = 101
    g         = 9.81m/s^2
    Ttop      = 0C

    # nondimensionalize
    CharDim   = GEO_units(length=10km, temperature=1000C, stress=1MPa, viscosity=1e20Pas)
    MatParam  = nondimensionalize(MatParam, CharDim)
    Thickness = nondimensionalize(Thickness, CharDim)
    g         = nondimensionalize(9.81m/s^2, CharDim)
    Ttop      = nondimensionalize(Ttop, CharDim)
    
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

    # build Figure
    fig = Figure(resolution = (1200, 900))
    ax1 = fig[1, 1] = Axis(fig,
        # title
        title  = "1D Strength Envelope",
        # x-axis
        xlabel = "Maximum Strength [MPa]",
        # y-axis
        ylabel = "Depth [km]"
    )

    # make sliders
    lsgrid = SliderGrid(fig[2, 1],
        (label = "Tbot [C]", range = 500:1:1200, startvalue = 800),
        (label = "log10(ε) [s-1]", range=-18:0.01:-10, startvalue = -15)
    )

    # create listener to sliders
    Tbot_slider  = lsgrid.sliders[1].value
    Tbot_dim     = @lift($Tbot_slider*C)

    exp_slider   = lsgrid.sliders[2].value
    ε_dim        = @lift(10^$exp_slider/s)
    ε            = @lift(nondimensionalize($ε_dim, CharDim))

    # build temperature structure
    Tbot      = @lift(nondimensionalize($Tbot_dim, CharDim))
    T         = @lift(GP_Compute_ThermalStructure(z, LinearTemp(Ttop, $Tbot)))

    # pressure and density
    ρ         = Observable(zeros(Float64, nz))
    P         = @lift(LithPres(MatParam, Phases, $ρ, $T, nz, dz, g))

    # solve for stress
    τ         = @lift(solveStress(MatParam, Phases, $ε, $P, $T, nz))
    
    # redimensionalize
    z         = dimensionalize(GeoUnit(z, km, false), CharDim).val
    τ         = @lift(redimensionalize!($τ, $ρ, $P, $T, CharDim))

    # plotting
    lines!(ax1, τ, -z)
    #lines(fig[1,1:2], ustrip(Value(τ)), -ustrip(Value(z)); axis = (; xlabel="Maximum Deviatoric Stress [MPa]", ylabel="Depth [km]"))
    #lines(fig[1,3], ustrip(Value(ρ)), -ustrip(Value(z)); axis = (; xlabel="Density [kg m-3]", ylabel="Depth [km]"))
    #lines(fig[1,4], ustrip(Value(P)), -ustrip(Value(z)); axis = (; xlabel="Lithostatic Pressure [MPa]", ylabel="Depth [km]"))
    ##lines(fig[1,4], ustrip(Value(T)), -ustrip(Value(z)); axis = (; xlabel="Temperature [C]", ylabel="Depth [km]"))
    #fig
    return fig
end