export VolumeBoundaryDiffusionCreep,
    SetVolumeBoundaryDiffusionCreep,
    Transform_VolumeBoundaryDiffusionCreep

"""
    VolumeBoundaryDiffusionCreep(; Name = "", d = 1mm, D0v = 9.1e-4m^2/s, Qv = 59.4kJ/mol, D0b = 6.0e-4m^2/s,
                        Qb = 49kJ/mol, δ = 9.04e-10m, Vm = 1.97e-5m^3/mol, Apparatus = AxialCompression)

Diffusion creep by volume (Nabarro-Herring) plus grain-boundary (Coble) diffusion, in the form of
eq. (4) of Goldsby & Kohlstedt (2001); the defaults are their values for ice,
```math
    \\dot{\\varepsilon} = \\frac{42 V_m \\sigma}{R T d^2}
    \\left( D_{0,v} e^{-Q_v/RT} + \\frac{\\pi \\delta}{d} D_{0,b} e^{-Q_b/RT} \\right)
```
where
- ``\\sigma`` is the differential stress and ``\\dot{\\varepsilon}`` the axial strain rate of the experiments
- ``d`` is the grain size ``[\\mathrm{m}]`` (default 1 mm), stored in the law
- ``D_{0,v}``, ``Q_v`` are the pre-exponential factor and activation energy of volume diffusion
- ``D_{0,b}``, ``Q_b`` are those of grain-boundary diffusion, and ``\\delta`` is the grain-boundary width
- ``V_m`` is the molar volume.

The prefactor scales with ``1/RT``, which [`DiffusionCreep`](@ref) cannot express (see #52).
The law is linear in stress. `Apparatus` converts it to second invariants as for the other creep laws;
with `AxialCompression` (the default, matching the uniaxial experiments) ``\\dot{\\varepsilon}_{II} = \\tau_{II}/(2\\eta)`` with
```math
    \\eta = \\frac{R T d^2}{126 V_m \\left( D_{0,v} e^{-Q_v/RT} + \\frac{\\pi \\delta}{d} D_{0,b} e^{-Q_b/RT} \\right)},
```
where ``126 = 2 \\cdot \\tfrac{3}{2} \\cdot 42`` and ``\\tfrac{3}{2} = F_T/F_E`` is the axial-compression correction.

The grain size used in the computations is the stored field `d`, unless a grain size is passed at
runtime, e.g. `compute_εII(law, τII; T, d = 1e-4)`, which overrides it.
"""
struct VolumeBoundaryDiffusionCreep{T, U1, U2, U3, U4, U5, U6, U7, U8} <: AbstractCreepLaw{T}
    Name::Ptr{UInt8}
    d::GeoUnit{T, U1}   # grain size
    D0v::GeoUnit{T, U2} # pre-exponential factor of volume diffusion
    Qv::GeoUnit{T, U3}  # activation energy of volume diffusion
    D0b::GeoUnit{T, U4} # pre-exponential factor of grain-boundary diffusion
    Qb::GeoUnit{T, U5}  # activation energy of grain-boundary diffusion
    δ::GeoUnit{T, U6}   # grain-boundary width
    Vm::GeoUnit{T, U7}  # molar volume
    R::GeoUnit{T, U8}   # universal gas constant
    Apparatus::Int8     # experimental apparatus: AxialCompression, SimpleShear or Invariant
    FT::T               # stress correction factor of the apparatus
    FE::T               # strain-rate correction factor of the apparatus

    function VolumeBoundaryDiffusionCreep(;
            Name = "",
            d = 1.0e-3m,
            D0v = 9.1e-4m^2 / s,
            Qv = 59.4e3J / mol,
            D0b = 6.0e-4m^2 / s,
            Qb = 49.0e3J / mol,
            δ = 9.04e-10m,
            Vm = 1.97e-5m^3 / mol,
            R = 8.3145J / mol / K,
            Apparatus = AxialCompression,
        )
        FT, FE = CorrectionFactor(Apparatus)
        vals = map(x -> convert(GeoUnit, x), (d, D0v, Qv, D0b, Qb, δ, Vm, R))
        T = typeof(first(vals)).types[1]
        U = map(x -> typeof(x).types[2], vals)
        return new{T, U...}(pointer(ptr2string(Name)), vals..., Int8(Apparatus), FT, FE)
    end
end

function VolumeBoundaryDiffusionCreep(Name, d, D0v, Qv, D0b, Qb, δ, Vm, R, Apparatus, FT, FE)
    return VolumeBoundaryDiffusionCreep(; Name, d, D0v, Qv, D0b, Qb, δ, Vm, R, Apparatus)
end

function remove_tensor_correction(s::VolumeBoundaryDiffusionCreep)
    (; d, D0v, Qv, D0b, Qb, δ, Vm, R) = s
    return VolumeBoundaryDiffusionCreep(; Name = unsafe_string(s.Name), d, D0v, Qv, D0b, Qb, δ, Vm, R, Apparatus = Invariant)
end

function param_info(s::VolumeBoundaryDiffusionCreep)
    name = unsafe_string(s.Name)
    eq = L"\dot{\varepsilon}_{II} = \frac{42 V_m F_T \tau_{II}}{F_E R T d^2}\left(D_{0,v} e^{-Q_v/RT} + \frac{\pi\delta}{d} D_{0,b} e^{-Q_b/RT}\right)"
    name == "" && return MaterialParamsInfo(; Equation = eq)
    f = find_creep_law(VolumeBoundaryDiffusion, name)
    isnothing(f) && return MaterialParamsInfo(; Equation = eq)
    inf = volume_boundary_diffusion_database_info(f)
    return MaterialParamsInfo(; Equation = eq, Comment = inf.Comment, BibTex_Reference = inf.BibTex_Reference)
end

# C in ε̇II = C τII; all arguments in consistent units.
@inline function _volume_boundary_diffusion_coefficient(D0v, Qv, D0b, Qb, δ, Vm, R, T, d, FT, FE)
    RT = R * T
    D = D0v * exp(-Qv / RT) + π * δ / d * D0b * exp(-Qb / RT)
    return 42 * Vm * D * FT / (RT * d^2 * FE)
end

@inline function _volume_boundary_diffusion_coefficient(a::VolumeBoundaryDiffusionCreep, ::Type{Tc}, T, d) where {Tc}
    @unpack_val Tc D0v, Qv, D0b, Qb, δ, Vm, R = a
    T = convert_precision(Tc, T)
    d = _grain_size(Tc, a, d)
    FT, FE = convert_precision(Tc, a.FT), convert_precision(Tc, a.FE)
    return _volume_boundary_diffusion_coefficient(D0v, Qv, D0b, Qb, δ, Vm, R, T, d, FT, FE)
end

@inline function _volume_boundary_diffusion_coefficient_units(a::VolumeBoundaryDiffusionCreep, ::Type{Tc}, T, d) where {Tc}
    @unpack_units Tc D0v, Qv, D0b, Qb, δ, Vm, R = a
    T = convert_precision(Tc, T)
    d = _grain_size_units(Tc, a, d)
    FT, FE = convert_precision(Tc, a.FT), convert_precision(Tc, a.FE)
    return _volume_boundary_diffusion_coefficient(D0v, Qv, D0b, Qb, δ, Vm, R, T, d, FT, FE)
end

"""
    compute_εII(a::VolumeBoundaryDiffusionCreep, TauII; T, d = a.d, kwargs...)

Returns the strain rate of volume and grain-boundary diffusion creep for the second invariant of the deviatoric stress `TauII`;
see [`VolumeBoundaryDiffusionCreep`](@ref).
"""
@inline function compute_εII(a::VolumeBoundaryDiffusionCreep, TauII; T = one(precision(a)), d = nothing, kwargs...)
    ε = _volume_boundary_diffusion_coefficient(a, precision_of(TauII), T, d) * TauII
    return retry_wider(ε, TauII) do τ
        compute_εII(a, τ; T, d)
    end
end

@inline function compute_εII(a::VolumeBoundaryDiffusionCreep, TauII::Quantity; T = 1K, d = nothing, kwargs...)
    return uconvert(s^-1, _volume_boundary_diffusion_coefficient_units(a, precision_of(TauII), T, d) * TauII)
end

@inline function dεII_dτII(a::VolumeBoundaryDiffusionCreep, TauII; T = one(precision(a)), d = nothing, kwargs...)
    C = _volume_boundary_diffusion_coefficient(a, precision_of(TauII), T, d)
    return retry_wider(C, one(TauII)) do τ
        dεII_dτII(a, τ; T, d)
    end
end

@inline function dεII_dτII(a::VolumeBoundaryDiffusionCreep, TauII::Quantity; T = 1K, d = nothing, kwargs...)
    return _volume_boundary_diffusion_coefficient_units(a, precision_of(TauII), T, d)
end

"""
    compute_τII(a::VolumeBoundaryDiffusionCreep, EpsII; T, d = a.d, kwargs...)

Returns the second invariant of the deviatoric stress of volume and grain-boundary diffusion creep for the strain rate `EpsII`.
"""
@inline function compute_τII(a::VolumeBoundaryDiffusionCreep, EpsII; T = one(precision(a)), d = nothing, kwargs...)
    τ = EpsII / _volume_boundary_diffusion_coefficient(a, precision_of(EpsII), T, d)
    return retry_wider(τ, EpsII) do ε
        compute_τII(a, ε; T, d)
    end
end

@inline function compute_τII(a::VolumeBoundaryDiffusionCreep, EpsII::Quantity; T = 1K, d = nothing, kwargs...)
    return uconvert(Pa, EpsII / _volume_boundary_diffusion_coefficient_units(a, precision_of(EpsII), T, d))
end

@inline function dτII_dεII(a::VolumeBoundaryDiffusionCreep, EpsII; T = one(precision(a)), d = nothing, kwargs...)
    dτ = inv(_volume_boundary_diffusion_coefficient(a, precision_of(EpsII), T, d))
    return retry_wider(dτ, one(EpsII)) do ε
        dτII_dεII(a, ε; T, d)
    end
end

@inline function dτII_dεII(a::VolumeBoundaryDiffusionCreep, EpsII::Quantity; T = 1K, d = nothing, kwargs...)
    return inv(_volume_boundary_diffusion_coefficient_units(a, precision_of(EpsII), T, d))
end

function compute_εII!(
        EpsII::AbstractArray{_T, N}, a::VolumeBoundaryDiffusionCreep, TauII::AbstractArray;
        T = one(_T), d = nothing, kwargs...,
    ) where {N, _T}
    for i in each_argument_index(EpsII, TauII, T, d)
        EpsII[i] = compute_εII(a, convert_precision(_T, TauII[i]); T = argument_at(T, i), d = argument_at(d, i))
    end
    return nothing
end

function compute_τII!(
        TauII::AbstractArray{_T, N}, a::VolumeBoundaryDiffusionCreep, EpsII::AbstractArray;
        T = one(_T), d = nothing, kwargs...,
    ) where {N, _T}
    for i in each_argument_index(TauII, EpsII, T, d)
        TauII[i] = compute_τII(a, convert_precision(_T, EpsII[i]); T = argument_at(T, i), d = argument_at(d, i))
    end
    return nothing
end

function show(io::IO, g::VolumeBoundaryDiffusionCreep)
    return print(
        io,
        "VolumeBoundaryDiffusionCreep: Name = $(unsafe_string(g.Name)), d=$(Value(g.d)), D0v=$(Value(g.D0v)), Qv=$(Value(g.Qv)), D0b=$(Value(g.D0b)), Qb=$(Value(g.Qb)), δ=$(Value(g.δ)), Vm=$(Value(g.Vm)), FT=$(g.FT), FE=$(g.FE)",
    )
end

include("Data/VolumeBoundaryDiffusionCreep.jl")
using .VolumeBoundaryDiffusion

"""
    SetVolumeBoundaryDiffusionCreep(name; d, D0v, Qv, D0b, Qb, δ, Vm)
    SetVolumeBoundaryDiffusionCreep(name, CharDim; ...)

Returns the pre-defined volume and grain-boundary diffusion creep law `name` (e.g. `GeoParams.VolumeBoundaryDiffusion.ice_Goldsby_2001`)
in SI units, nondimensionalized if `CharDim` is given. Keywords override the stored parameters.
"""
function SetVolumeBoundaryDiffusionCreep(name::F; d = nothing, D0v = nothing, Qv = nothing, D0b = nothing, Qb = nothing, δ = nothing, Vm = nothing) where {F}
    return Transform_VolumeBoundaryDiffusionCreep(name, (; d, D0v, Qv, D0b, Qb, δ, Vm))
end

function SetVolumeBoundaryDiffusionCreep(name::F, CharDim::GeoUnits{T}; kwargs...) where {F, T <: Union{GEO, SI}}
    return nondimensionalize(SetVolumeBoundaryDiffusionCreep(name; kwargs...), CharDim)
end

"""
    Transform_VolumeBoundaryDiffusionCreep(name[, kwargs::NamedTuple])

Converts the parameters of an volume and grain-boundary diffusion creep law to SI units; non-`nothing` entries of `kwargs`
replace the stored values.
"""
Transform_VolumeBoundaryDiffusionCreep(name::F) where {F} = Transform_VolumeBoundaryDiffusionCreep(volume_boundary_diffusion_database(name))
Transform_VolumeBoundaryDiffusionCreep(name::F, kwargs::NamedTuple) where {F} = Transform_VolumeBoundaryDiffusionCreep(volume_boundary_diffusion_database(name), kwargs)

function Transform_VolumeBoundaryDiffusionCreep(name::F, CharDim::GeoUnits{U}) where {F, U <: Union{GEO, SI}}
    return nondimensionalize(Transform_VolumeBoundaryDiffusionCreep(name), CharDim)
end

function Transform_VolumeBoundaryDiffusionCreep(pp::VolumeBoundaryDiffusionCreep, kwargs::NamedTuple = (;))
    si(k, u) = uconvert(u, Value(convert(GeoUnit, something(get(kwargs, k, nothing), getfield(pp, k)))))
    return VolumeBoundaryDiffusionCreep(;
        Name = unsafe_string(pp.Name),
        d = si(:d, m), D0v = si(:D0v, m^2 / s), Qv = si(:Qv, J / mol), D0b = si(:D0b, m^2 / s),
        Qb = si(:Qb, J / mol), δ = si(:δ, m), Vm = si(:Vm, m^3 / mol), Apparatus = pp.Apparatus,
    )
end
