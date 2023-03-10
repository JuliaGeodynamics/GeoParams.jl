export PeierlsCreep,
    Peierls_info,
    SetPeierlsCreep,
    remove_tensor_correction,
    dεII_dτII,
    dτII_dεII

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
either `AxialCompression`, `SimpleShear` or `Invariant`. If the flow law paramters are already given as a function of second invariants, choose `Apparatus=Invariant`.

# Example
```julia-repl 
julia> x2 = PeierlsCreep(n=1)
PeierlsCreep: n=1, A=1.5 MPa^-3 s^-1, E=476.0 kJ mol^-1, Apparatus=AxialCompression
```
"""
struct PeierlsCreep{T,N,U1,U2,U3,U4,U5} <: AbstractCreepLaw{T}
    Name::NTuple{N,Char}
    n::GeoUnit{T,U1} # power-law exponent
    q::GeoUnit{T,U1} # stress relation exponent
    o::GeoUnit{T,U1} # ... (normally called p but used as 'o' since p already exists)
    TauP::GeoUnit{T, U2} # Peierls stress
    A::GeoUnit{T,U3} # material specific rheological parameter
    E::GeoUnit{T,U4} # activation energy
    R::GeoUnit{T,U5} # universal gas constant
    Apparatus::Int8 # type of experimental apparatus, either AxialCompression, SimpleShear or Invariant
    FT::T # type of experimental apparatus, either AxialCompression, SimpleShear or Invariant
    FE::T # type of experimental apparatus, either AxialCompression, SimpleShear or Invariant

    function PeierlsCreep(;
        Name="",
        n=1.0NoUnits,
        q=2.0NoUnits,
        o=1.0NoUnits,
        TauP=8.5e9Pa,
        A=5.7e11s^(-1.0),
        E=476.0kJ / mol,
        R=8.3145J / mol / K,
        Apparatus=AxialCompression,
    )

        # Rheology name
        Name = String(join(Name))
        N = length(Name)
        NameU = NTuple{N,Char}(collect.(Name))
        
        # Corrections from lab experiments
        FT, FE = CorrectionFactor(Apparatus)
        # Convert to GeoUnits
        nU = n isa GeoUnit ? n : convert(GeoUnit, n)
        qU = q isa GeoUnit ? q : convert(GeoUnit, q)
        oU = o isa GeoUnit ? o : convert(GeoUnit, o)
        TauPU = TauP isa GeoUnit ? TauP : convert(GeoUnit, TauP)
        AU = A isa GeoUnit ? A : convert(GeoUnit, A)
        EU = E isa GeoUnit ? E : convert(GeoUnit, E)
        RU = R isa GeoUnit ? R : convert(GeoUnit, R)
        # Extract struct types
        T = typeof(nU).types[1]
        U1 = typeof(nU).types[2]
        U2 = typeof(TauPU).types[2]
        U3 = typeof(AU).types[2]
        U4 = typeof(EU).types[2]
        U5 = typeof(RU).types[2]
        
        # Create struct
        return new{T,N,U1,U2,U3,U4,U5}(
            NameU, nU, qU, oU, TauPU, AU, EU, RU, Int8(Apparatus), FT, FE
        )
    end

    function PeierlsCreep(Name, n, q, o, TauP, A, E, R, Apparatus, FT, FE)
        return PeierlsCreep(;
            Name=Name, n=n, q=q, o=o, TauP=TauP, A=A, E=E, R=R, Apparatus=Apparatus
        )
    end
end

"""
    Transforms units from GPa, MPa, kJ etc. to basic units such as Pa, J etc.
"""

function Transform_PeierlsCreep(name; kwargs)
    p_in = PeierlsCreep_info[name][1]

    # Take optional arguments 
    v_kwargs = values(kwargs)
    val = GeoUnit.(values(v_kwargs))
    
    args = (Name=p_in.Name, n=p_in.n, q=p_in.q, o=p_in.o, A=p_in.A, E=p_in.E, Apparatus=p_in.Apparatus)
    p = merge(args, NamedTuple{keys(v_kwargs)}(val))
    
    Name = String(collect(p.Name))
    n = Value(p.n)
    q = Value(p.q)
    o = Value(p.o)
    TauP = uconvert(Pa, Value(p.TauP))
    A_Pa = uconvert(Pa^(-NumValue(p.n)) / s, Value(p.A))
    E_J = uconvert(J / mol, Value(p.E))

    Apparatus = p.Apparatus

    # args from database
    args = (Name=Name, n=n, q=q, o=o, TauP=TauP, A=A_Pa, E=E_J, Apparatus=Apparatus)
    
    return PeierlsCreep(; args...)
end

"""
    s = remove_tensor_correction(s::PeierlsCreep)

Removes the tensor correction of the creeplaw, which is useful to compare the implemented creeplaws
with the curves of the original publications, as those publications usually do not transfer their data to tensor format
"""
function remove_tensor_correction(s::PeierlsCreep)
    name = String(collect(s.Name))

    return PeierlsCreep(;
        Name=name, n=s.n, q=s.q, o=s.o, TauP=s.TauP, A=s.A, E=s.E, Apparatus=Invariant
    )
end

function param_info(s::PeierlsCreep)
    name = String(collect(s.Name))
    eq = L"\tau_{ij} = 2 \eta  \dot{\varepsilon}_{ij}"
    if name == ""
        return MaterialParamsInfo(; Equation=eq)
    end
    inf = PeierlsCreep_info[name][2]
    return MaterialParamsInfo(;
        Equation=eq, Comment=inf.Comment, BibTex_Reference=inf.BibTex_Reference
    )
end

# Calculation routines for linear viscous rheologies
# All inputs must be non-dimensionalized (or converted to consitent units) GeoUnits
@inline function compute_εII(
    a::PeierlsCreep, TauII::_T; T=one(precision(a)), args...
) where {_T}
    @unpack_val n, q, o, TauP, A, E, R = a
    FT, FE = a.FT, a.FE

    ε = A * exp(-(E / (R * T)) * ((1 - (TauII / TauP)^o)^q)) / FE
    return ε
end

@inline function compute_εII(
    a::PeierlsCreep, TauII::Quantity; T=1K, args...
)
    @unpack_units n, q, o, TauP, A, E, R = a
    FT, FE = a.FT, a.FE

    ε = A * exp(-(E / (R * T)) * ((1 - (TauII / TauP)^o)^q)) / FE

    return ε
end

function compute_εII!(
    EpsII::AbstractArray{_T,N},
    a::PeierlsCreep,
    TauII::AbstractArray{_T,N};
    T=ones(size(TauII))::AbstractArray{_T,N},
    kwargs...,
) where {N,_T}
    @inbounds for i in eachindex(EpsII)
        EpsII[i] = compute_εII(a, TauII[i]; T=T[i])
    end

    return nothing
end

@inline function dεII_dτII(
    a::PeierlsCreep, TauII::_T; T=one(precision(a)), args...
) where {_T}
    @unpack_val n, q, o, TauP, A, E, R = a
    FT, FE = a.FT, a.FE

    return o * 
           fastpow((FT * TauII) / TauP, o) *
           fastpow(1 - fastpow((FT * TauII) / TauP, o), q) *
           A *
           E *
           q *
           FT *
           n *
           exp(-(E * fastpow(1 - fastpow((FT * TauII) / TauP, o), q) / (R * T))) *
           (1 / (FE * R * T * Tau * (1 - fastpow((FT * TauII) / TauP, o))))
end

@inline function dεII_dτII(
    a::PeierlsCreep, TauII::Quantity; T=1K, args...
)
    @unpack_units n, q, o, TauP, A, E, R = a
    FT, FE = a.FT, a.FE

    return o * 
           fastpow((FT * TauII) / TauP, o) *
           fastpow(1 - fastpow((FT * TauII) / TauP, o), q) *
           A *
           E *
           q *
           FT *
           n *
           exp(-(E * fastpow(1 - fastpow((FT * TauII) / TauP, o), q) / (R * T))) *
           (1 / (FE * R * T * Tau * (1 - fastpow((FT * TauII) / TauP, o))))
end


"""
    compute_τII(a::PeierlsCreep, EpsII; P, T, f, args...)

Computes the stress for a peierls creep law given a certain strain rate

"""
@inline function compute_τII(
    a::PeierlsCreep, EpsII::_T; T=one(precision(a)), args...
) where {_T}
    local n, q, o, TauP, A, E, R
    if EpsII isa Quantity
        @unpack_units n, q, o, TauP, A, E, R = a
    else
        @unpack_val n, q, o, TauP, A, E, R = a
    end

    FT, FE = a.FT, a.FE

    τ = TauP *
        fastpow(1 - fastpow(- (R * T * log((FE * EpsII) / A) /Q), 
        1 / q), 
        1 / o) / 
        FT

    return τ
end

@inline function compute_τII(
    a::PeierlsCreep, EpsII::Quantity; T=1K, args...
)
    @unpack_units n, q, o, TauP, A, E, R = a
    FT, FE = a.FT, a.FE

    τ = TauP *
        fastpow(1 - fastpow(- (R * T * log((FE * EpsII) / A) /Q), 
        1 / q), 
        1 / o) / 
        FT

    return τ
end

"""
    compute_τII!(TauII::AbstractArray{_T,N}, a::PeierlsCreep, EpsII::AbstractArray{_T,N}; 
        T = ones(size(TauII))::AbstractArray{_T,N}, 

Computes the deviatoric stress invariant for a peierls creep law
"""
function compute_τII!(
    TauII::AbstractArray{_T,N},
    a::PeierlsCreep,
    EpsII::AbstractArray{_T,N};
    T=ones(size(TauII))::AbstractArray{_T,N},
    kwargs...,
) where {N,_T}
    @inbounds for i in eachindex(TauII)
        TauII[i] = compute_τII(a, EpsII[i]; T=T[i])
    end

    return nothing
end
#
@inline function dτII_dεII(
    a::PeierlsCreep, EpsII::_T; T=one(precision(a)), args...
) where {_T}
    @unpack_val n, q, o, TauP, A, E, R = a
    FT, FE = a.FT, a.FE

    # derived in WolframAlpha
    return (TauP * exp(-1 / q) * 
            R * 
            T * 
            fastpow(-R * T * log((FE * EpsII) / A)), (1 / q - 1)) * 
            fastpow(1 - exp(-1 / q) * fastpow(-R * T * log((FE * EpsII) / A), (1 / q)), (1/o - 1)) / 
            (FT * o * q *  EpsII)
end

@inline function dτII_dεII(
    a::PeierlsCreep, EpsII::Quantity; T=1K, args...
)
    @unpack_units n, q, o, TauP, A, E, R = a
    FT, FE = a.FT, a.FE

    # derived in WolframAlpha
    return (TauP * exp(-1 / q) * 
            R * 
            T * 
            fastpow(-R * T * log((FE * EpsII) / A)), (1 / q - 1)) * 
            fastpow(1 - exp(-1 / q) * fastpow(-R * T * log((FE * EpsII) / A), (1 / q)), (1/o - 1)) / 
            (FT * o * q *  EpsII)
end

# Print info 
function show(io::IO, g::PeierlsCreep)
    return print(
        io,
        "PeierlsCreep: Name = $(String(collect(g.Name))), n=$(Value(g.n)), q=$(Value(g.q)), o=$(Value(g.o)), TauP=$(Value(g.TauP)), A=$(Value(g.A)), E=$(Value(g.E)), FT=$(g.FT), FE=$(g.FE), Apparatus=$(g.Apparatus)",
    )
end
#-------------------------------------------------------------------------

# load collection of peierls creep laws
include("Data/PeierlsCreep.jl")
