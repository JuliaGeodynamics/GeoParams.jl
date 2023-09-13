using Parameters
using SpecialFunctions: erfc

export StrengthEnvelopeComp, StrengthEnvelopePlot, ConstTemp, LinTemp, HalfspaceCoolTemp, CompTempStruct, LithPres

abstract type AbstractTempStruct end  

"""
    ConstTemp(T=1000)
    
Sets a constant temperature inside the box

Parameters:
- T : the value
"""
@with_kw_noshow struct ConstTemp <: AbstractTempStruct
    T = 1000
end


"""
    CompTempStruct(Z, s::AbstractTempStruct)

Returns a temperature vector that matches the coordinate vector and temperature structure that were passed in

Parameters:
- Z: vector with depth coordinates
- s: Temperature structure (CostTemp, LinTemp, HalfspaceCoolTemp)
"""
CompTempStruct(Z, s::ConstTemp) = fill(s.T, length(Z))


"""
    LinTemp(Ttop=0, Tbot=1000)
    
Sets a linear temperature structure from top to bottom

Parameters:
- Ttop: the value @ the top
- Tbot: the value @ the bottom
"""
@with_kw_noshow struct LinTemp <: AbstractTempStruct
    Ttop = 0
    Tbot = 1350
end

function CompTempStruct(Z, s::LinTemp)
    @unpack Ttop, Tbot  = s

    dz    = Z[end]-Z[1];
    dT    = Tbot - Ttop

    Temp  = abs.(Z./dz).*dT .+ Ttop
    return Temp
end

"""
    HalfspaceCoolTemp(Tsurface=0, Tmantle=1350, Age=60, Adiabat=0)
    
Sets a halfspace temperature structure in plate

Parameters:
- Tsurface: surface temperature [C]
- Tmantle:  mantle temperature [C]
- Age:      Thermal Age of plate [Myrs]
- Adiabat:  Mantle Adiabat [K/km]
"""
@with_kw_noshow struct HalfspaceCoolTemp <: AbstractTempStruct
    Tsurface = 0        # top T
    Tmantle  = 1350     # bottom T
    Age      = 60       # thermal age of plate [in Myrs]
    Adiabat  = 0        # Adiabatic gradient in K/km
    kappa    = 1e-6     # thermal diffusivity
end

function CompTempStruct(Z, t::HalfspaceCoolTemp)
    @unpack Tsurface, Tmantle, Age, Adiabat, kappa  = t    

    Temp  = zeros(Float64, length(Z))K

    MantleAdiabaticT    =   Tmantle .+ Adiabat*abs.(Z);   # Adiabatic temperature of mantle
    
    for i in eachindex(Temp)
        Temp[i] =   (Tsurface .- Tmantle)*erfc(NoUnits.((abs.(Z[i]))./(2*sqrt(kappa*Age)))) + MantleAdiabaticT[i];
    end
    return Temp
end

"""
    LithPres(MatParam, Phases, ρ, T, dz, g)

Iteratively solves a 1D lithostatic pressure profile (compatible with temperature- and pressure-dependent densities)

Parameters:
- MatParam: a tuple of materials (including the following properties: Phase, Density)
- Phases:   vector with the distribution of phases
- ρ:        density vector for initial guess(can be zeros)
- T:        temperature vector
- dz:       grid spacing
- g:        gravitational acceleration
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
            c = Mat.Plasticity[1].C.val
            τ[i] = Pres * Mat.Plasticity[1].sinϕ.val + Mat.Plasticity[1].cosϕ.val * c
        end
    end

    return τ
end

function extractFromResult(res, ind, nz)
    x = zeros(Float64, nz)
    for i = 1 : nz
        x[i] = res[ind][i].val
    end
    return x
end

"""
    StrengthEnvelopeComp(MatParam::NTuple{N, AbstractMaterialParamsStruct}, Thickness::Vector{U}, TempType::AbstractTempStruct=LinTemp(0C, 800C), ε=1e-15/s, nz::Int64=101) where {N, U}

Calculates a 1D strength envelope. Pressure used for Drucker Prager plasticity is lithostatic. To visualize the results in a GUI, use StrengthEnvelopePlot.

Parameters:
- MatParam:  a tuple of materials (including the following properties: Phase, Density, CreepLaws, Plasticity)
- Thickness: a vector listing the thicknesses of the respective layers (should carry units)
- TempType:  the type of temperature profile (ConstTemp(), LinTemp(), HalfspaceCoolTemp())
- ε:         background strainrate
- nz:        optional argument controlling the number of points along the profile (default = 101)

# Example:
```julia-repl
julia> using GLMakie
julia> MatParam = (SetMaterialParams(Name="UC", Phase=1, Density=ConstantDensity(ρ=2700kg/m^3), CreepLaws = SetDislocationCreep("Wet Quartzite | Ueda et al. (2008)"), Plasticity = DruckerPrager(ϕ=30.0, C=10MPa)),
                   SetMaterialParams(Name="MC", Phase=2, Density=Density=ConstantDensity(ρ=2900kg/m^3), CreepLaws = SetDislocationCreep("Plagioclase An75 | Ji and Zhao (1993)"), Plasticity = DruckerPrager(ϕ=20.0, C=10MPa)),
                   SetMaterialParams(Name="LC", Phase=3, Density=PT_Density(ρ0=2900kg/m^3, α=3e-5/K, β=1e-10/Pa), CreepLaws = SetDislocationCreep("Maryland strong diabase | Mackwell et al. (1998)"), Plasticity = DruckerPrager(ϕ=30.0, C=10MPa)));
julia> Thickness = [15,10,15]*km;

julia> StrengthEnvelopeComp(MatParam, Thickness, LinTemp(), ε=1e-15/s)
```
"""
function StrengthEnvelopeComp(MatParam::NTuple{N, AbstractMaterialParamsStruct}, Thickness::Vector{U}, TempType::AbstractTempStruct=LinTemp(0C, 800C), ε=1e-15/s, nz::Int64=101) where {N, U}

    # hardcoded input
    g         = 9.81m/s^2

    # nondimensionalize
    CharDim   = GEO_units(length=10km, temperature=1000C, stress=10MPa, viscosity=1e20Pas)
    MatParam  = nondimensionalize(MatParam,  CharDim)
    Thickness = nondimensionalize(Thickness, CharDim)
    g         = nondimensionalize(9.81m/s^2, CharDim)
    ε         = nondimensionalize(ε,         CharDim)
    
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

    # build temperature structure
    T         = CompTempStruct(dimensionalize(z, km, CharDim), TempType)
    T         = nondimensionalize(T, CharDim)

    # pressure and density
    ρ         = zeros(Float64, nz)
    P         = LithPres(MatParam, Phases, ρ, T, dz, g)

    # solve for stress
    τ         = solveStress(MatParam, Phases, ε, P, T)

    # redimensionalize
    z         = dimensionalize(z, km,  CharDim)
    τ         = dimensionalize(τ, MPa, CharDim)
    T         = dimensionalize(T, C,   CharDim)

    return z, τ, T
end


"""
    StrengthEnvelopePlot(MatParam, Thickness; TempType, nz)

Requires GLMakie

Creates a GUI that plots a 1D strength envelope. In the GUI, temperature profile and strain rate can be adjusted. The Drucker-Prager plasticity uses lithostatic pressure.

Parameters:
- MatParam:  a tuple of materials (including the following properties: Phase, Density, CreepLaws, Plasticity)
- Thickness: a vector listing the thicknesses of the respective layers (should carry units)
- TempType:  the type of temperature profile (LinTemp=default, HalfspaceCoolTemp, ConstTemp)
- nz:        optional argument controlling the number of points along the profile (default = 101)

# Example:
```julia-repl
julia> using GLMakie
julia> MatParam = (SetMaterialParams(Name="UC", Phase=1, Density=ConstantDensity(ρ=2700kg/m^3), CreepLaws = SetDislocationCreep("Wet Quartzite | Ueda et al. (2008)"), Plasticity = DruckerPrager(ϕ=30.0, C=10MPa)),
                   SetMaterialParams(Name="MC", Phase=2, Density=Density=ConstantDensity(ρ=2900kg/m^3), CreepLaws = SetDislocationCreep("Plagioclase An75 | Ji and Zhao (1993)"), Plasticity = DruckerPrager(ϕ=20.0, C=10MPa)),
                   SetMaterialParams(Name="LC", Phase=3, Density=PT_Density(ρ0=2900kg/m^3, α=3e-5/K, β=1e-10/Pa), CreepLaws = SetDislocationCreep("Maryland strong diabase | Mackwell et al. (1998)"), Plasticity = DruckerPrager(ϕ=30.0, C=10MPa)));
julia> Thickness = [15,10,15]*km;

julia> StrengthEnvelopePlot(MatParam, Thickness, LinTemp())
```
"""
function StrengthEnvelopePlot()
    println("Load GLMakie to use this function: \n julia> using GLMakie") 
end
