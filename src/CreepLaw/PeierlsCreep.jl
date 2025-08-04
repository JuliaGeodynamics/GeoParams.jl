export PeierlsCreep,
    Peierls_info,
    remove_tensor_correction,
    dεII_dτII,
    dτII_dεII,
    Transform_PeierlsCreep

# Peierls Creep ------------------------------------------------
"""
    PeierlsCreep(n = 1.0NoUnits, r = 0.0NoUnits, A = 1.5MPa/s, E = 476.0kJ/mol, V = 6e-6m^3/mol, apparatus = AxialCompression )

Defines the flow law parameter of a peierls creep law.

The peierls creep law, as used by experimentalists, is given by
```math
    \\dot{\\gamma} = A \\exp{\\left( -\\frac{E}{RT} \\left( 1 - \\left( \\frac{\\sigma_d}{\\sigma_p} \\right) ^o \\right)^q\\right)}
```
where
- ``n`` is the power law exponent
- ``q`` is the stress relation exponent
- ``o`` is the stress relation exponent (normally called p)
- ``r`` is the exponent of fugacity dependence
- ``A`` is a pre-exponential factor ``[\\mathrm{MPa}^{-n}s^{-1}]`` (if manually defined, ``n`` must be either pre-defined or substituted)
- ``E`` is the activation energy ``\\mathrm{[kJ/mol]}``
- ``V`` is the activation volume ``\\mathrm{[m^3/mol]}``
- ``\\dot{\\gamma}`` is the strain rate ``\\mathrm{[1/s]}``
- ``\\sigma_\\mathrm{p}`` is the peierls stress ``\\mathrm{[MPa]}``
- ``\\sigma_\\mathrm{d}`` is the differential stress ``\\mathrm{[MPa]}`` which are converted into second invariants using the `Apparatus` variable that can be
either `AxialCompression`, `SimpleShear` or `Invariant`. If the flow law parameters are already given as a function of second invariants, choose `Apparatus=Invariant`.

# Example
```julia-repl
julia> x2 = PeierlsCreep(n=1)
PeierlsCreep: Name = , n=1.0, q=2.0, o=1.0, TauP=8.5e9 Pa, A=5.7e11 s^-1.0, E=476.0 kJ mol^-1.0, FT=1.7320508075688772, FE=1.1547005383792517, Apparatus=1
```
"""
struct PeierlsCreep{T, U1, U2, U3, U4, U5} <: AbstractCreepLaw{T}
    Name::Ptr{UInt8} # name of the creep law
    n::GeoUnit{T, U1} # power-law exponent
    q::GeoUnit{T, U1} # stress relation exponent
    o::GeoUnit{T, U1} # ... (normally called p but used as 'o' since p already exists)
    TauP::GeoUnit{T, U2} # Peierls stress
    A::GeoUnit{T, U3} # material specific rheological parameter
    E::GeoUnit{T, U4} # activation energy
    R::GeoUnit{T, U5} # universal gas constant
    Apparatus::Int8 # type of experimental apparatus, either AxialCompression, SimpleShear or Invariant
    FT::T # type of experimental apparatus, either AxialCompression, SimpleShear or Invariant
    FE::T # type of experimental apparatus, either AxialCompression, SimpleShear or Invariant

    function PeierlsCreep(;
            Name = "",
            n = 1.0NoUnits,
            q = 2.0NoUnits,
            o = 1.0NoUnits,
            TauP = 8.5e9Pa,
            A = 5.7e11s^(-1),
            E = 476.0e3J / mol,
            R = 8.3145J / mol / K,
            Apparatus = AxialCompression,
        )
        # Corrections from lab experiments
        FT, FE = CorrectionFactor(Apparatus)
        # Convert to GeoUnits
        nU = convert(GeoUnit, n)
        qU = convert(GeoUnit, q)
        oU = convert(GeoUnit, o)
        TauPU = convert(GeoUnit, TauP)
        AU = convert(GeoUnit, A)
        EU = convert(GeoUnit, E)
        RU = convert(GeoUnit, R)
        # Extract struct types
        T = typeof(nU).types[1]
        U1 = typeof(nU).types[2]
        U2 = typeof(TauPU).types[2]
        U3 = typeof(AU).types[2]
        U4 = typeof(EU).types[2]
        U5 = typeof(RU).types[2]
        name = pointer(ptr2string(Name))
        # Create struct
        return new{T, U1, U2, U3, U4, U5}(
            name, nU, qU, oU, TauPU, AU, EU, RU, Int8(Apparatus), FT, FE
        )
    end

    function PeierlsCreep(Name, n, q, o, TauP, A, E, R, Apparatus, FT, FE)
        return PeierlsCreep(;
            Name = Name, n = n, q = q, o = o, TauP = TauP, A = A, E = E, R = R, Apparatus = Apparatus
        )
    end
end

"""
    s = remove_tensor_correction(s::PeierlsCreep)

Removes the tensor correction of the creeplaw, which is useful to compare the implemented creeplaws
with the curves of the original publications, as those publications usually do not transfer their data to tensor format.
"""
function remove_tensor_correction(s::PeierlsCreep)
    name = ptr2string(s.Name)

    return PeierlsCreep(;
        Name = name, n = s.n, q = s.q, o = s.o, TauP = s.TauP, A = s.A, E = s.E, Apparatus = Invariant
    )
end

function param_info(s::PeierlsCreep)
    name = ptr2string(s.Name)
    eq = L"\tau_{ij} = 2 \eta  \dot{\varepsilon}_{ij}"
    if name == ""
        return MaterialParamsInfo(; Equation = eq)
    end
    inf = peierls_database_info(name)
    return MaterialParamsInfo(;
        Equation = eq, Comment = inf.Comment, BibTex_Reference = inf.BibTex_Reference
    )
end

# Calculation routines for linear viscous rheologies
# All inputs must be non-dimensionalized (or converted to consistent units) GeoUnits
@inline function compute_εII(
        a::PeierlsCreep, TauII::_T; T = one(precision(a)), args...
    ) where {_T}
    @unpack_val n, q, o, TauP, A, E, R = a
    FT, FE = a.FT, a.FE

    ε = @pow (A * exp(-(E / (R * T)) * ((1.0 - ((FT * TauII) / TauP)^o)^q))) / FE

    return ε
end

@inline function compute_εII(a::PeierlsCreep, TauII::Quantity; T = 1K, args...)
    @unpack_units n, q, o, TauP, A, E, R = a
    FT, FE = a.FT, a.FE

    ε = @pow (A * exp(-(E / (R * T)) * ((1 - (((FT * TauII) / TauP)^o))^q))) / FE

    return ε
end

function compute_εII!(
        EpsII::AbstractArray{_T, N},
        a::PeierlsCreep,
        TauII::AbstractArray{_T, N};
        T = ones(size(TauII))::AbstractArray{_T, N},
        kwargs...,
    ) where {N, _T}
    @inbounds for i in eachindex(EpsII)
        EpsII[i] = compute_εII(a, TauII[i]; T = T[i])
    end

    return nothing
end

function dεII_dτII(a::PeierlsCreep, TauII; args...)
    return ForwardDiff.derivative(x -> compute_εII(a, x; args...), TauII)
end

"""
    compute_τII(a::PeierlsCreep, EpsII; P, T, f, args...)

Computes the stress for a peierls creep law given a certain strain rate.

"""
@inline function compute_τII(
        a::PeierlsCreep, EpsII::_T; T = one(precision(a)), args...
    ) where {_T}
    @unpack_val n, q, o, TauP, A, E, R = a

    FT, FE = a.FT, a.FE

    q_inv = inv(q)
    o_inv = inv(o)

    τ = @pow (TauP * (1.0 - (-((R * T * log((FE * EpsII) / A)) / E))^q_inv)^o_inv) / FT

    return τ
end

@inline function compute_τII(a::PeierlsCreep, EpsII::Quantity; T = 1K, args...)
    @unpack_units n, q, o, TauP, A, E, R = a
    FT, FE = a.FT, a.FE

    q_inv = inv(q)
    o_inv = inv(o)

    τ = @pow (TauP * (1.0 - (-((R * T * log((FE * EpsII) / A)) / E))^q_inv)^o_inv) / FT

    return τ
end

"""
    compute_τII!(TauII::AbstractArray{_T,N}, a::PeierlsCreep, EpsII::AbstractArray{_T,N};
        T = ones(size(TauII))::AbstractArray{_T,N}, args...)

Computes the deviatoric stress invariant for a peierls creep law.

"""
function compute_τII!(
        TauII::AbstractArray{_T, N},
        a::PeierlsCreep,
        EpsII::AbstractArray{_T, N};
        T = ones(size(TauII))::AbstractArray{_T, N},
        kwargs...,
    ) where {N, _T}
    @inbounds for i in eachindex(TauII)
        TauII[i] = compute_τII(a, EpsII[i]; T = T[i])
    end

    return nothing
end

"""
    dτII_dεII(v::PeierlsCreep, EpsII; args...)

Computes the derivative `dτII/dεII` for a peierls creep law using automatic differentiation.

"""

function dτII_dεII(v::PeierlsCreep, EpsII; args...)
    return ForwardDiff.derivative(x -> compute_τII(v, x; args...), EpsII)
end

# Print info
function show(io::IO, g::PeierlsCreep)
    return print(
        io,
        "PeierlsCreep: Name = $(unsafe_string(g.Name)), n=$(Value(g.n)), q=$(Value(g.q)), o=$(Value(g.o)), TauP=$(Value(g.TauP)), A=$(Value(g.A)), E=$(Value(g.E)), FT=$(g.FT), FE=$(g.FE), Apparatus=$(g.Apparatus)",
    )
end
#-------------------------------------------------------------------------

# load collection of peierls creep laws
include("Data/PeierlsCreep.jl")

using .Peierls
export SetPeierlsCreep

"""
    SetPeierlsCreep["Name of peierls creep law"]
This is a dictionary with pre-defined creep laws
"""
function SetPeierlsCreep(
        name::F;
        n = nothing,
        q = nothing,
        o = nothing,
        A = nothing,
        E = nothing,
        TauP = nothing,
    ) where {F}
    kwargs = (; n, q, o, A, E, TauP)
    return Transform_PeierlsCreep(name, kwargs)
end

function SetPeierlsCreep(
        name::F,
        CharDim::GeoUnits{T};
        n = nothing,
        q = nothing,
        o = nothing,
        A = nothing,
        E = nothing,
        TauP = nothing,
    ) where {F, T <: Union{GEO, SI}}
    kwargs = (; n, q, o, A, E, TauP)
    return nondimensionalize(Transform_PeierlsCreep(name, kwargs), CharDim)
end

"""
    Transforms units from GPa, MPa, kJ etc. to basic units such as Pa, J etc.
"""
Transform_PeierlsCreep(name::F) where {F} = Transform_PeierlsCreep(peierls_database(name))
Transform_PeierlsCreep(name::F, kwargs::NamedTuple) where {F} = Transform_PeierlsCreep(peierls_database(name), kwargs)

function Transform_PeierlsCreep(name::F, CharDim::GeoUnits{U}) where {U <: Union{GEO, SI}} where {F}
    return Transform_PeierlsCreep(peierls_database(name), CharDim)
end

function Transform_PeierlsCreep(p::AbstractCreepLaw{T}, CharDim::GeoUnits{U}) where {T, U <: Union{GEO, SI}}
    return nondimensionalize(Transform_PeierlsCreep(p), CharDim)
end

function Transform_PeierlsCreep(p::AbstractCreepLaw{T}) where {T}
    n = Value(p.n)
    q = Value(p.q)
    o = Value(p.o)
    TauP = uconvert(Pa, Value(p.TauP))
    A_Pa = uconvert(s^(-1), Value(p.A))
    E_J = uconvert(J / mol, Value(p.E))

    Apparatus = p.Apparatus

    # args from database
    args = (Name = p.Name, n = n, q = q, o = o, TauP = TauP, A = A_Pa, E = E_J, Apparatus = Apparatus)

    return PeierlsCreep(; args...)
end

function Transform_PeierlsCreep(p::AbstractCreepLaw{T}, kwargs::NamedTuple) where {T}
    (; n, q, o, A, E, TauP) = kwargs

    n_new = isnothing(n) ? Value(p.n) : Value(GeoUnit(n))
    q_new = isnothing(q) ? Value(p.q) : Value(GeoUnit(q))
    o_new = isnothing(o) ? Value(p.o) : Value(GeoUnit(o))
    A_new = isnothing(A) ? p.A : GeoUnit(A)
    E_new = isnothing(E) ? Value(p.E) : Value(GeoUnit(E))
    TauP_new = isnothing(E) ? Value(p.TauP) : Value(GeoUnit(TauP))

    A_Pa = uconvert(s^(-1), Value(A_new))
    E_J = uconvert(J / mol, E_new)
    Apparatus = p.Apparatus

    # args from database
    args = (Name = p.Name, n = n_new, q = q_new, o = o_new, TauP = TauP_new, A = A_Pa, E = E_J, Apparatus = Apparatus)

    return PeierlsCreep(; args...)
end
