export MeltViscosity,
    compute_εII!,
    compute_εII,
    compute_τII!,
    compute_τII

#=----Melt Viscosity---

Defines effective shear viscosity and effective bulk viscosity as a function of melt fraction and geometry of the pores following parameterization of *Schmeling, A., Kruse, J. P., Richard, G., Effective shear and bulk viscosity of partially molten rock based on elastic moduli theory of a fluid filled poroelastic medium, Geophysical Journal International, Volume 190, Issue 3, September 2012, Pages 1571–1578, https://doi.org/10.1111/j.1365-246X.2012.05596.x*

α is the aspect ratio (between 1 and 0)
a1
a2
b1
b2
k1
k2
k3
c1
c2
c3


=#

struct MeltViscosity{T,N,U1} <: AbstractCreepLaw{T}
    Name::NTuple{N,Char}
    α::GeoUnit{T,U1} # aspect ratio
    a1::GeoUnit{T,U1} # fitting constant
    a2::GeoUnit{T,U1} # fitting constant
    b1::GeoUnit{T,U1} # fitting constant
    b2::GeoUnit{T,U1} # fitting constant
    k2::GeoUnit{T,U1} # fitting constant
    k3::GeoUnit{T,U1} # fitting constant
    c3::GeoUnit{T,U1} # fitting constant

    function MeltViscosity(;
        Name="",
        α= 1.0NoUnits,
        a1= 0.97NoUnits,
        a2= 0.8NoUnits,
        b1= 2.2455NoUnits,
        b2= 3.45NoUnits,
        k2= 1.25NoUnits,
        k3= 1.29NoUnits,
        c3= 2.4NoUnits,
    )

        # Rheology name
        Name = String(join(Name))
        N = length(Name)
        NameU = NTuple{N,Char}(collect.(Name))
        # Corrections from lab experiments
        αU = α isa GeoUnit ? α : convert(GeoUnit, α)
        a1U = a1 isa GeoUnit ? a1 : convert(GeoUnit, a1)
        a2U = a2 isa GeoUnit ? a2 : convert(GeoUnit, a2)
        b1U = b1 isa GeoUnit ? b1 : convert(GeoUnit, b1)
        b2U = b2 isa GeoUnit ? b2 : convert(GeoUnit, b2)
        k2U = k2 isa GeoUnit ? k2 : convert(GeoUnit, k2)
        k3U = k3 isa GeoUnit ? k3 : convert(GeoUnit, k3)
        c3U = c3 isa GeoUnit ? c3 : convert(GeoUnit, c3)
        # Extract struct types
        T = typeof(αU).types[1]
        U1 = typeof(αU).types[2]
        # Create struct
        return new{T,N,U1}(
            NameU, αU, a1U, a2U, b1U, b2U, k2U, k3U, c3U
        )
    end

    function MeltViscosity(Name, α, a1, a2, b1, b2, k2, k3, c3)
        return MeltViscosity(;
            Name=Name, α=α, a1=a1, a2=a2, b1=b1, b2=b2, k2=k2, k3=k3, c3=c3
        )
    end
end

@inline function compute_εII(
    a::MeltViscosity, TauII::_T; η=one(precision(a)), ϕ=one(precision(a)), kwargs...
) where {_T}
    @unpack_val α,a1,a2,b1,b2,k2,k3,c3 = a

    k1 = a1 * (a2 + α * (1 - a2))
    c1 = (b1 * α) / (1 + b2 * fastpow(α,k3))

    if ϕ < c1
        η_melt = η * (1 - ϕ / c1)^k1
    else
        η_melt = 0
    end

    @show ϕ
    @show η
    ε = TauII / (2*η_melt)

    return ε
end

@inline function compute_εII(
    a::MeltViscosity, TauII::Quantity; η=1e18Pa * s, ϕ=0.10NoUnits, args...
)
    @unpack_val α,a1,a2,b1,b2,k2,k3,c3 = a

    k1 = a1 * (a2 + α * (1 - a2))
    c1 = (b1 * α) / (1 + b2 * fastpow(α,k3))

    if ϕ < c1
        η_melt = η * (1 - ϕ / c1)^k1
    else
        η_melt = 0
    end

    ε = TauII / (2*η_melt)

    @show ϕ
    @show η
    @show c1
    @show η_melt

    return ε
end

"""
    compute_εII!(EpsII::AbstractArray{_T,N}, a, TauII::AbstractArray{_T,N}; η, ϕ, kwargs...)

Computes strainrate as a function of stress
"""
function compute_εII!(
    EpsII::AbstractArray{_T,N},
    a::MeltViscosity,
    TauII::AbstractArray{_T,N};
    η=ones(size(TauII))::AbstractArray{_T,N},
    ϕ=ones(size(TauII))::AbstractArray{_T,N},
    kwargs...,
) where {N,_T}
    @inbounds for I in eachindex(EpsII)
        EpsII[I] = compute_εII(a, TauII[I]; η=η[I], ϕ=ϕ[I])
    end

    return nothing
end

"""
    computeCreepLaw_TauII(EpsII::_T, a::DiffusionCreep; T::_T, P=zero(_T), f=one(_T), d=one(_T), kwargs...)

Returns diffusion creep stress as a function of 2nd invariant of the strain rate
"""
@inline function compute_τII(
    a::MeltViscosity, EpsII::_T;η=one(precision(a)), ϕ=zero(precision(a)), kwargs...
) where {_T}
    @unpack_val α,a1,a2,b1,b2,k2,k3,c3 = a
    k1 = a1 * (a2 + α * (1 - a2))
    c1 = (b1 * α) / (1 + b2 * fastpow(α,k3))

    if ϕ < c1
        η_melt = η * (1 - ϕ / c1)^k1
    else
        η_melt = 0
    end

    τ = 2 * η_melt * EpsII

    return τ
end

@inline function compute_τII(
    a::MeltViscosity, EpsII::Quantity; η=10e18Pa.s, ϕ=10NoUnits, kwargs...
)
    @unpack_val α,a1,a2,b1,b2,k2,k3,c3 = a
    k1 = a1 * (a2 + α * (1 - a2))
    c1 = (b1 * α) / (1 + b2 * fastpow(α,k3))

    if ϕ < c1
        η_melt = η * (1 - ϕ / c1)^k1
    else
        η_melt = 0
    end

    τ = 2 * η_melt * EpsII

    return τ
end

function compute_τII!(
    TauII::AbstractArray{_T,N},
    a::MeltViscosity,
    EpsII::AbstractArray{_T,N};
    η=ones(size(TauII))::AbstractArray{_T,N},
    ϕ=ones(size(TauII))::AbstractArray{_T,N},
    kwargs...,
) where {N,_T}
    @inbounds for I in eachindex(EpsII)
        TauII[I] = compute_τII(a, EpsII[I]; η=η[I], ϕ=ϕ[I])
    end

    return nothing
end


# Print info
function show(io::IO, g::MeltViscosity)
    return print(
        io,
        "MeltViscosity: Name = $(String(collect(g.Name))), α=$(Value(g.α))",
    )
end

