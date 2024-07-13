# This implements viscous creep laws and routines to compute with them
#
# Note that various simple creep laws are defined in this file; 
# more complex ones (such as DislocationCreep) are in separate files 
# in the same directory
#
# In case you want to add new creep laws, have a look at how the ones
# here are implemented. Please add tests as well!

# include("Data/DislocationCreep_data.jl")

export DislocationCreep,
    Transform_DislocationCreep,
    remove_tensor_correction,
    dεII_dτII,
    dτII_dεII

# Dislocation Creep ------------------------------------------------
"""
    DislocationCreep(n = 1.0NoUnits, r = 0.0NoUnits, A = 1.5MPa/s, E = 476.0kJ/mol, V = 6e-6m^3/mol, apparatus = AxialCompression )
    
Defines the flow law parameter of a dislocation creep law.

The (isotropic) dislocation creep law, as used by experimentalists, is given by  
```math  
     \\dot{\\gamma} = A \\sigma_\\mathrm{d}^n f_\\mathrm{H2O}^r \\exp\\left(-\\frac{E+PV}{RT}\\right)
```
where 
- ``n`` is the power law exponent  
- ``r`` is the exponent of fugacity dependence 
- ``A`` is a pre-exponential factor ``[\\mathrm{MPa}^{-n}s^{-1}]`` (if manually defined, ``n`` and ``r`` must be either pre-defined or substituted) 
- ``E`` is the activation energy ``\\mathrm{[kJ/mol]}`` 
- ``V`` is the activation volume ``\\mathrm{[m^3/mol]}`` 
- ``\\dot{\\gamma}`` is the strain rate ``\\mathrm{[1/s]}`` 
- ``\\sigma_\\mathrm{d}`` is the differential stress ``\\mathrm{[MPa]}`` which are converted into second invariants using the `Apparatus` variable that can be
either `AxialCompression`, `SimpleShear` or `Invariant`. If the flow law parameters are already given as a function of second invariants, choose `Apparatus=Invariant`.

# Example
```julia-repl 
julia> x2 = DislocationCreep(n=3)
DislocationCreep: n=3, r=0.0, A=1.5 MPa^-3 s^-1, E=476.0 kJ mol^-1, V=6.0e-6 m^3 mol^-1, Apparatus=AxialCompression
```
"""
struct DislocationCreep{T,U1,U2,U3,U4,U5,U6} <: AbstractCreepLaw{T}
    Name::Ptr{UInt8}
    n::GeoUnit{T,U1} # power-law exponent
    r::GeoUnit{T,U1} # exponent of water-fugacity
    A::GeoUnit{T,U2} # material specific rheological parameter
    E::GeoUnit{T,U3} # activation energy
    V::GeoUnit{T,U4} # activation volume
    R::GeoUnit{T,U5} # universal gas constant
    A_to_minus_n::GeoUnit{T,U6} # universal gas constant
    Apparatus::Int8 # type of experimental apparatus, either AxialCompression, SimpleShear or Invariant
    FT::T # type of experimental apparatus, either AxialCompression, SimpleShear or Invariant
    FE::T # type of experimental apparatus, either AxialCompression, SimpleShear or Invariant

    function DislocationCreep(;
        Name="",
        n=1NoUnits,
        r=0NoUnits,
        A=1.5e6Pa^(-n) / s,
        E=476.0e3J / mol,
        V=6e-6m^3 / mol,
        R=8.3145J / mol / K,
        Apparatus=AxialCompression,
    )
        # Corrections from lab experiments
        FT, FE = CorrectionFactor(Apparatus)
        # Convert to GeoUnits
        nU = convert(GeoUnit, rat2float(n))
        rU = convert(GeoUnit, r)
        AU = convert(GeoUnit, A)
        EU = convert(GeoUnit, E)
        VU = convert(GeoUnit, V)
        RU = convert(GeoUnit, R)
        u  = (Pa^(-nU.val) / s) ^ (-1/nU.val)
        A_to_minus_n = convert(GeoUnit, (AU.val^(-1/nU.val)) * u)
        # Extract struct types
        T  = typeof(nU).types[1]
        U1 = typeof(nU).types[2]
        U2 = typeof(AU).types[2]
        U3 = typeof(EU).types[2]
        U4 = typeof(VU).types[2]
        U5 = typeof(RU).types[2]
        U6 = typeof(A_to_minus_n).types[2]
        name = pointer(ptr2string(Name))
        # Create struct
        return new{T,U1,U2,U3,U4,U5,U6}(
            name, nU, rU, AU, EU, VU, RU, A_to_minus_n, Int8(Apparatus), FT, FE
        )
    end

    function DislocationCreep(Name, n, r, A, E, V, R, A_to_minus_nU, Apparatus, FT, FE)
        return DislocationCreep(;
            Name=Name, n=n, r=r, A=A, E=E, V=V, R=R, Apparatus=Apparatus
        )
    end
end

"""
    s = remove_tensor_correction(s::DiffusionCreep)

Removes the tensor correction of the creeplaw, which is useful to compare the implemented creeplaws
with the curves of the original publications, as those publications usually do not transfer their data to tensor format
"""
function remove_tensor_correction(s::DislocationCreep)
    name = ptr2string(s.Name)

    return DislocationCreep(;
        Name=name, n=s.n, r=s.r, A=s.A, E=s.E, V=s.V, Apparatus=Invariant
    )
end

# Calculation routines for linear viscous rheologies
# All inputs must be non-dimensionalized (or converted to consistent units) GeoUnits
@inline function compute_εII(
    a::DislocationCreep,
    TauII;
    T=one(precision(a)),
    P=zero(precision(a)),
    f=one(precision(a)),
    args...,
)
    @unpack_val n, r, A, E, V, R = a
    FT, FE = a.FT, a.FE

    ε = @pow A * (TauII * FT)^n * f^r * exp(-(E + P * V) / (R * T)) / FE
    return ε
end

@inline function compute_εII(
    a::DislocationCreep, TauII::Quantity; T=1K, P=0Pa, f=1NoUnits, args...
)
    @unpack_units n, r, A, E, V, R = a
    FT, FE = a.FT, a.FE

    ε = @pow A * (TauII * FT)^n * f^r * exp(-(E + P * V) / (R * T)) / FE

    return ε
end

function compute_εII!(
    EpsII::AbstractArray{_T,N},
    a::DislocationCreep,
    TauII::AbstractArray{_T,N};
    P=zero(TauII)::AbstractArray{_T,N},
    T=ones(size(TauII))::AbstractArray{_T,N},
    f=ones(size(TauII))::AbstractArray{_T,N},
    kwargs...,
) where {N,_T}
    @inbounds for i in eachindex(EpsII)
        EpsII[i] = compute_εII(a, TauII[i]; T=T[i], P=P[i], f=f[i])
    end

    return nothing
end

@inline function dεII_dτII(
    a::DislocationCreep,
    TauII;
    T=one(precision(a)),
    P=zero(precision(a)),
    f=one(precision(a)),
    args...,
)
    @unpack_val n, r, A, E, V, R = a
    FT, FE = a.FT, a.FE

    return @pow (FT * TauII)^(n - 1) *
        f^r *
        A *
        FT *
        n *
        exp(-(E + P * V) / (R * T)) *
        inv(FE)
end

@inline function dεII_dτII(
    a::DislocationCreep, TauII::Quantity; T=1K, P=0Pa, f=1NoUnits, args...
)
    @unpack_units n, r, A, E, V, R = a
    FT, FE = a.FT, a.FE

    return @pow (FT * TauII)^(n - 1) *
        f^r *
        A *
        FT *
        n *
        exp(-(E + P * V) / (R * T)) *
        inv(FE)
end

"""
    compute_τII(a::DislocationCreep, EpsII; P, T, f, args...)

Computes the stress for a Dislocation creep law given a certain strain rate

"""
@inline function compute_τII(
    a::DislocationCreep,
    EpsII;
    T=one(precision(a)),
    P=zero(precision(a)),
    f=one(precision(a)),
    args...,
)
    n, r, A, E, V, R = if EpsII isa Quantity
        @unpack_units n, r, A, E, V, R, A_to_minus_n = a
        n, r, A, E, V, R
    else
        @unpack_val n, r, A, E, V, R, A_to_minus_n = a
        n, r, A, E, V, R
    end

    FT, FE = a.FT, a.FE
    _n = inv(n)

    return @pow A_to_minus_n * (EpsII * FE)^_n * f^(-r * _n) * exp((E + P * V) / (n * R * T)) / FT
end

@inline function compute_τII(
    a::DislocationCreep, EpsII::Quantity; P=0Pa, T=1K, f=1NoUnits, args...
)
    @unpack_units n, r, A, E, V, R, A_to_minus_n = a
    FT, FE = a.FT, a.FE
    _n = inv(n)

    return @pow A_to_minus_n * f^(-r * _n) * (EpsII * FE)^_n * exp((E + P * V) / (n * R * T)) / FT
end

"""
    compute_τII!(TauII::AbstractArray{_T,N}, a::DislocationCreep, EpsII::AbstractArray{_T,N}; 
        P =       zero(TauII)::AbstractArray{_T,N}, 
        T = ones(size(TauII))::AbstractArray{_T,N}, 
        f = ones(size(TauII))::AbstractArray{_T,N})

Computes the deviatoric stress invariant for a dislocation creep law
"""
function compute_τII!(
    TauII::AbstractArray{_T,N},
    a::DislocationCreep,
    EpsII::AbstractArray{_T,N};
    T=ones(size(TauII))::AbstractArray{_T,N},
    P=zero(TauII)::AbstractArray{_T,N},
    f=ones(size(TauII))::AbstractArray{_T,N},
    kwargs...,
) where {N,_T}
    @inbounds for i in eachindex(TauII)
        TauII[i] = compute_τII(a, EpsII[i]; T=T[i], P=P[i], f=f[i])
    end

    return nothing
end

@inline function dτII_dεII(
    a::DislocationCreep,
    EpsII;
    T=one(precision(a)),
    P=zero(precision(a)),
    f=one(precision(a)),
    args...,
)
    @unpack_val n, r, A, E, V, R, A_to_minus_n = a
    FT, FE = a.FT, a.FE
    _n = inv(n)

    return @pow (
        FE * A_to_minus_n * f^(-r * _n) * (EpsII * FE)^(_n - 1) * exp((E + P * V) / (R * T * n))
    ) / (FT * n)
end

@inline function dτII_dεII(
    a::DislocationCreep, EpsII::Quantity; P=0Pa, T=1K, f=1NoUnits, args...
)
    @unpack_units n, r, A, E, V, R, A_to_minus_n = a
    FT, FE = a.FT, a.FE
    _n = inv(n)

    return @pow (
        FE * A_to_minus_n * f^(-r * _n) * (EpsII * FE)^(_n - 1) * exp((E + P * V) / (R * T * n))
    ) / (FT * n)
end

# Print info 
function show(io::IO, g::DislocationCreep)
    return print(
        io,
        "DislocationCreep: Name = $(unsafe_string(g.Name)), n=$(Value(g.n)), r=$(Value(g.r)), A=$(Value(g.A)), E=$(Value(g.E)), V=$(Value(g.V)), FT=$(g.FT), FE=$(g.FE), Apparatus=$(g.Apparatus)",
    )
end
#-------------------------------------------------------------------------

# load collection of dislocation creep laws
include("Data/DislocationCreep.jl")
include("Data_deprecated/DislocationCreep.jl")

using .Dislocation
export SetDislocationCreep

function param_info(s::DislocationCreep)
    name = unsafe_string(s.Name)
    eq = L"\tau_{ij} = 2 \eta  \dot{\varepsilon}_{ij}"
    if name == ""
        return MaterialParamsInfo(; Equation=eq)
    end
    inf = DislocationCreep_info[name][2]
    return MaterialParamsInfo(;
        Equation=eq, Comment=inf.Comment, BibTex_Reference=inf.BibTex_Reference
    )
end

function SetDislocationCreep(
    name::F;
    n = nothing,
    r = nothing,
    A = nothing,
    E = nothing,
    V = nothing,
) where F 
    kwargs = (; n, r, A, E, V)
    Transform_DislocationCreep(name, kwargs)
end

function SetDislocationCreep(
    name::F,
    CharDim::GeoUnits{T};
    n = nothing,
    r = nothing,
    A = nothing,
    E = nothing,
    V = nothing,
) where {F, T<:Union{GEO, SI}}
    kwargs = (; n, r, A, E, V)
    nondimensionalize(Transform_DislocationCreep(name, kwargs), CharDim)
end

"""
    Transforms units from MPa, kJ etc. to basic units such as Pa, J etc.
"""
Transform_DislocationCreep(name::F) where F = Transform_DislocationCreep(dislocation_database(name))
Transform_DislocationCreep(name::F, kwargs::NamedTuple) where F = Transform_DislocationCreep(dislocation_database(name), kwargs)

function Transform_DislocationCreep(name::F, CharDim::GeoUnits{U}) where {U<:Union{GEO,SI}} where F
    Transform_DislocationCreep(dislocation_database(name), CharDim)
end

function Transform_DislocationCreep(p::AbstractCreepLaw{T}, CharDim::GeoUnits{U}) where {T,U<:Union{GEO,SI}}
    nondimensionalize(Transform_DislocationCreep(p), CharDim)
end

function Transform_DislocationCreep(p::AbstractCreepLaw{T}) where T
    n         = Value(p.n)
    A_Pa      = uconvert(Pa^unit_power(p.A) / s, Value(p.A))
    E_J       = uconvert(J / mol, Value(p.E))
    V_m3      = uconvert(m^3 / mol, Value(p.V))
    Apparatus = p.Apparatus
    r         = Value(p.r)
    # args from database
    args = (Name=unsafe_string(p.Name), n=n, r=r, A=A_Pa, E=E_J, V=V_m3, Apparatus=Apparatus)

    return DislocationCreep(; args...)
end

function Transform_DislocationCreep(p::AbstractCreepLaw{T}, kwargs::NamedTuple) where T
    
    (; n, r, A, E, V) = kwargs

    n_new = isnothing(n) ? Value(p.n) : Value(GeoUnit(n))
    r_new = isnothing(r) ? Value(p.r) : Value(GeoUnit(r))
    A_new = isnothing(A) ?        p.A : GeoUnit(A)
    E_new = isnothing(E) ? Value(p.E) : Value(GeoUnit(E))
    V_new = isnothing(E) ? Value(p.V) : Value(GeoUnit(V))
   
    A_Pa      = uconvert(Pa^unit_power(A_new) / s, Value(A_new))
    E_J       = uconvert(J / mol, E_new)
    V_m3      = uconvert(m^3 / mol, V_new)
    Apparatus = p.Apparatus
    # args from database
    args = (Name=unsafe_string(p.Name), n=n_new, r=r_new, A=A_Pa, E=E_J, V=V_m3, Apparatus=Apparatus)

    return DislocationCreep(; args...)
end