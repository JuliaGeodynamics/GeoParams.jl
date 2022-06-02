# This implements viscous creep laws and routines to compute with them
#
# Note that various simple creep laws are defined in this file; 
# more complex ones (such as DislocationCreep) are in separate files 
# in the same directory
#
# In case you want to add new creep laws, have a look at how the ones
# here are implemented. Please add tests as well!

export DislocationCreep,
    DislocationCreep_info,
    SetDislocationCreep,
    remove_tensor_correction,
    dεII_dτII,
    dτII_dεII

# Dislocation Creep ------------------------------------------------
"""
    DislocationCreep(n = 1.0NoUnits, r = 0.00.0NoUnits, A = 1.5MPa/s, E = 476.0kJ/mol, V = 6e-6m^3/mol, apparatus = AxialCompression )
    
Defines the flow law parameter of a dislocation creep law 
The (isotropic) dislocation creep law, as used by experimtalists, is given by  
```math  
     \\dot{\\gamma} = A \\sigma_\\mathrm{d}^n f_\\mathrm{H2O}^r \\exp(-\\frac{E+PV}{RT})
```
where ``n`` is the power law exponent,  
``r`` is the exponent of fugacity dependence, 
``A`` is a pre-exponential factor [MPa^(n+r)] (if manually defined, n and r must be either pre-defined or substituted),  
``E`` is the activation energy [kJ/mol], ``V`` is the activation volume [m^3/mol]. ``\\dot{\\gamma}`` is the ordinary strain rate [1/s], 
and ``\\sigma_\\mathrm{d}`` is the differential stress which are converted into second invariants using the apparatus type that can be
either AxialCompression, SimpleShear or Invariant.
If the flow law paramters are already given as a function of second invariants, choose apparatus = "Invariant"
```julia-repl 
julia> x2      =   DislocationCreep(n=3)
DislocationCreep: n=3, r=0.0, A=1.5 MPa^-3 s^-1, E=476.0 kJ mol^-1, V=6.0e-6 m^3 mol^-1, Apparatus=AxialCompression
```
"""
# @with_kw_noshow struct DislocationCreep{T,N,U1,U2,U3,U4,U5} <: AbstractCreepLaw{T}
#     Name::NTuple{N,Char}    =   ""               # The name is encoded as a NTuple{Char} to make it isbits    
#     n::GeoUnit{T,U1}        = 1.0NoUnits         # power-law exponent
#     r::GeoUnit{T,U1}        = 0.0NoUnits         # exponent of water-fugacity dependence
#     A::GeoUnit{T,U2}        = 1.5MPa^(-n-r)/s    # pre-exponential factor
#     E::GeoUnit{T,U3}        = 476.0kJ/mol        # activation energy
#     V::GeoUnit{T,U4}        = 6e-6m^3/mol        # activation volume
#     R::GeoUnit{T,U5}        = 8.3145J/mol/K      # Universal gas constant
#     Apparatus::Int64        = 1                  # type of experimental apparatus, either AxialCompression, SimpleShear or Invariant
# end
# DislocationCreep(args...) = DislocationCreep(NTuple{length(args[1]), Char}(
#                                 collect.(args[1])), 
#                                 convert.(GeoUnit,args[2:end-1])..., 
#                                 args[end])

struct DislocationCreep{T,N,U1,U2,U3,U4,U5} <: AbstractCreepLaw{T}
    Name::NTuple{N,Char}
    n::GeoUnit{T,U1} # power-law exponent
    r::GeoUnit{T,U1} # exponent of water-fugacity
    A::GeoUnit{T,U2} # material specific rheological parameter
    E::GeoUnit{T,U3} # activation energy
    V::GeoUnit{T,U4} # activation volume
    R::GeoUnit{T,U5}  # universal gas constant
    Apparatus::Int8 # type of experimental apparatus, either AxialCompression, SimpleShear or Invariant
    FT::T # type of experimental apparatus, either AxialCompression, SimpleShear or Invariant
    FE::T # type of experimental apparatus, either AxialCompression, SimpleShear or Invariant

    function DislocationCreep(;
        Name="",
        n=1.0NoUnits,
        r=0.0NoUnits,
        A=1.5MPa^(-n - r) / s,
        E=476.0kJ / mol,
        V=6e-6m^3 / mol,
        R=8.3145J / mol / K,
        Apparatus=AxialCompression,
    )

        # Rheology name
        N = length(Name)
        NameU = NTuple{N,Char}(collect.(Name))
        # Corrections from lab experiments
        FT, FE = CorrectionFactor(Apparatus)
        # Convert to GeoUnits
        nU = n isa GeoUnit ? n : convert(GeoUnit, n)
        rU = r isa GeoUnit ? r : convert(GeoUnit, r)
        AU = A isa GeoUnit ? A : convert(GeoUnit, A)
        EU = E isa GeoUnit ? E : convert(GeoUnit, E)
        VU = V isa GeoUnit ? V : convert(GeoUnit, V)
        RU = R isa GeoUnit ? R : convert(GeoUnit, R)
        # Extract struct types
        T = typeof(nU).types[1]
        U1 = typeof(nU).types[2]
        U2 = typeof(AU).types[2]
        U3 = typeof(EU).types[2]
        U4 = typeof(VU).types[2]
        U5 = typeof(RU).types[2]
        # Create struct
        return new{T,N,U1,U2,U3,U4,U5}(
            NameU, nU, rU, AU, EU, VU, RU, Int8(Apparatus), FT, FE
        )
    end

    function DislocationCreep(Name, n, r, A, E, V, R, Apparatus, FT, FE)

        # Rheology name
        N = length(Name)
        NameU = NTuple{N,Char}(collect.(Name))
        # Convert to GeoUnits
        nU = n isa GeoUnit ? n : convert(GeoUnit, n)
        rU = r isa GeoUnit ? r : convert(GeoUnit, r)
        AU = A isa GeoUnit ? A : convert(GeoUnit, A)
        EU = E isa GeoUnit ? E : convert(GeoUnit, E)
        VU = V isa GeoUnit ? V : convert(GeoUnit, V)
        RU = R isa GeoUnit ? R : convert(GeoUnit, R)
        # Extract struct types
        T = typeof(nU).types[1]
        U1 = typeof(nU).types[2]
        U2 = typeof(AU).types[2]
        U3 = typeof(EU).types[2]
        U4 = typeof(VU).types[2]
        U5 = typeof(RU).types[2]
        # Create struct
        return new{T,N,U1,U2,U3,U4,U5}(
            NameU, nU, rU, AU, EU, VU, RU, Int8(Apparatus), FT, FE
        )
    end
end

"""
    Transforms units from MPa, kJ etc. to basic units such as Pa, J etc.
"""

function Transform_DislocationCreep(name)
    p = DislocationCreep_info[name][1]

    Name =  String(collect(p.Name))
    n    =  Value(p.n)
    A_Pa =  Value(p.A) |> Pa^(-NumValue(p.n) - NumValue(p.r))/s
    E_J  =  Value(p.E) |> J/mol
    V_m3 =  Value(p.V) |> m^3/mol
    
    Apparatus = p.Apparatus
    r = Value(p.r)

    return DislocationCreep(;
        Name=Name, n=n, r=r, A=A_Pa, E=E_J, V=V_m3, Apparatus=Apparatus
    )
end

"""
    s = remove_tensor_correction(s::DiffusionCreep)

Removes the tensor correction of the creeplaw, which is useful to compare the implemented creeplaws
with the curves of the original publications, as those publications usually do not transfer their data to tensor format
"""
function remove_tensor_correction(s::DislocationCreep)
    name = String(collect(s.Name))

    return DislocationCreep(;
        Name=name, n=s.n, r=s.r, A=s.A, E=s.E, V=s.V, Apparatus=Invariant
    )
end

function param_info(s::DislocationCreep)
    name = String(collect(s.Name))
    eq = L"\tau_{ij} = 2 \eta  \dot{\varepsilon}_{ij}"
    if name == ""
        return MaterialParamsInfo(; Equation=eq)
    end
    inf = DislocationCreep_info[name][2]
    return MaterialParamsInfo(;
        Equation=eq, Comment=inf.Comment, BibTex_Reference=inf.BibTex_Reference
    )
end

# Calculation routines for linear viscous rheologies
# All inputs must be non-dimensionalized (or converted to consitent units) GeoUnits
function compute_εII(
    a::DislocationCreep, TauII::_T; T::_T, P::_T=zero(_T), f::_T=one(_T), args...
) where {_T}
    @unpack_val n, r, A, E, V, R = a
    FT, FE = a.FT, a.FE

    ε = A * fastpow(TauII * FT, n) * fastpow(f, r) * exp(-(E + P * V) / (R * T)) / FE
    return ε
end

function compute_εII(a::DislocationCreep, TauII::Quantity; T=1K, P=0Pa, f=1NoUnits, args...)
    @unpack_units n, r, A, E, V, R = a
    FT, FE = a.FT, a.FE

    ε = A * (TauII * FT)^n * f^r * exp(-(E + P * V) / (R * T)) / FE

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

function dεII_dτII(
    a::DislocationCreep, TauII::_T; T::_T=one(_T), P::_T=zero(_T), f::_T=one(_T), args...
) where {_T}
    @unpack_val n, r, A, E, V, R = a
    FT, FE = a.FT, a.FE

    return fastpow(FT * TauII, -1 + n) *
           fastpow(f, r) *
           A *
           FT *
           n *
           exp((-E - P * V) / (R * T)) *
           (1 / FE)
end

"""
    compute_τII(a::DislocationCreep, EpsII; P, T, f, args...)

Computes the stress for a Dislocation creep law given a certain strain rate

"""
function compute_τII(
    a::DislocationCreep, EpsII::_T; T::_T=one(_T), P::_T=zero(_T), f::_T=one(_T), args...
) where {_T}
    local n, r, A, E, V, R
    if EpsII isa Quantity
        @unpack_units n, r, A, E, V, R = a
    else
        @unpack_val n, r, A, E, V, R = a
    end

    FT, FE = a.FT, a.FE

    return fastpow(A, -1 / n) *
           fastpow(EpsII * FE, 1 / n) *
           fastpow(f, -r / n) *
           exp((E + P * V) / (n * R * T)) / FT
end

function compute_τII(
    a::DislocationCreep, EpsII::Quantity; P=0Pa, T=1K, f=1NoUnits, args...
) where {_T}
    @unpack_units n, r, A, E, V, R = a
    FT, FE = a.FT, a.FE

    τ = A^(-1 / n) * (EpsII * FE)^(1 / n) * f^(-r / n) * exp((E + P * V) / (n * R * T)) / FT

    return τ
end

"""
    compute_τII!(TauII::AbstractArray{_T,N}, a::DislocationCreep, EpsII::AbstractArray{_T,N}; 
        P =       zero(TauII)::AbstractArray{_T,N}, 
        T = ones(size(TauII))::AbstractArray{_T,N}, 
        f = ones(size(TauII))::AbstractArray{_T,N})

Computes the stress for a Dislocation creep law
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

function dτII_dεII(
    a::DislocationCreep, EpsII::_T; T::_T=one(_T), P::_T=zero(_T), f::_T=one(_T), args...
) where {_T}
    @unpack_val n, r, A, E, V, R = a
    FT, FE = a.FT, a.FE

    return (
        FE *
        (A^(-1 / n)) *
        (f^((-r) / n)) *
        ((EpsII * FE)^(1 / n - 1)) *
        exp((E + P * V) / (R * T * n))
    ) / (FT * n)
end

# Print info 
function show(io::IO, g::DislocationCreep)
    return print(
        io,
        "DislocationCreep: Name = $(String(collect(g.Name))), n=$(Value(g.n)), r=$(Value(g.r)), A=$(Value(g.A)), E=$(Value(g.E)), V=$(Value(g.V)), FT=$(g.FT), FE=$(g.FE), Apparatus=$(g.Apparatus)",
    )
end
#-------------------------------------------------------------------------

# Add pre-defined creep laws 
include("DislocationCreep_Data.jl")
