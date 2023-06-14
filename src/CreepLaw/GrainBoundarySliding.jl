export GrainBoundarySliding,
    SetGrainBoundarySliding,
    GrainBoundarySliding_info,
    remove_tensor_correction,
    compute_εII!,
    compute_εII,
    compute_τII!,
    compute_τII

#=----Grain Boundary Sliding---
Defines grain boundary sliding parameters

p is the negative defined grain size exponent (value is set negative in the calculation of σ and ε)
A is a material specific rheological parameter
E is the activation energy
V is the activation volume
R is the universal gas constant
Apparatus defines the appartus type that shall be recreated (Axial Compression, Simple Shear, Invariant)
=#
"""
    GrainBoundarySliding(p = A = 1.5MPa/s, E = 476.0kJ/mol, V = 6e-6m^3/mol, apparatus = AxialCompression )
    
Defines the flow law parameter of a grain boundary sliding.

The grain boundary sliding, as used by experimentalists, is given by  
```math  
     \\dot{\\gamma} = A \\sigma_\\mathrm{d} d^{\\mathrm{p}} \\exp\\left(-\\frac{E+PV}{RT}\\right) #ÄNDERN!!!!!!!!
```
where 
- ``p`` is the exponent of grain size
- ``A`` is a pre-exponential factor ``[\\mathrm{MPa}^{-n}s^{-1}]`` (if manually defined, ``n`` must be either pre-defined or substituted) 
- ``E`` is the activation energy ``\\mathrm{[kJ/mol]}`` 
- ``V`` is the activation volume ``\\mathrm{[m^3/mol]}`` 
- ``\\dot{\\gamma}`` is the strain rate ``\\mathrm{[1/s]}`` 
- ``\\sigma_\\mathrm{d}`` is the differential stress ``\\mathrm{[MPa]}``

The experimental paramaters are converted into second invariants using the `Apparatus` variable that can be
either `AxialCompression`, `SimpleShear` or `Invariant`. If the flow law paramters are already given as a function of second invariants, choose `Apparatus=Invariant`.

# Example
```julia-repl 
julia> x2 = GrainBoundarySliding(Name="test")
GrainBoundarySliding: Name = test, n=1.0, p=-3.0, A=1.5 m³·⁰ MPa⁻¹·⁰ s⁻¹·⁰, E=500.0 kJ mol⁻¹·⁰, V=2.4e-5 m³·⁰ mol⁻¹·⁰, FT=1.7320508075688772, FE=1.1547005383792517)
```
"""

struct GrainBoundarySliding{T,N,U1,U2,U3,U4,U5} <: AbstractCreepLaw{T}
    Name::NTuple{N,Char}
    n::GeoUnit{T,U1} # powerlaw exponent
    p::GeoUnit{T,U1} # grain size exponent
    A::GeoUnit{T,U2} # material specific rheological parameter
    E::GeoUnit{T,U3} # activation energy
    V::GeoUnit{T,U4} # activation volume
    R::GeoUnit{T,U5} # universal gas constant
    Apparatus::Int8 # type of experimental apparatus, either AxialCompression, SimpleShear or Invariant
    FT::T # type of experimental apparatus, either AxialCompression, SimpleShear or Invariant
    FE::T # type of experimental apparatus, either AxialCompression, SimpleShear or Invariant

    function GrainBoundarySliding(;
        Name="",
        n=3.5NoUnits,
        p=-2.0NoUnits,
        A=6500MPa^(-3.5) * s^(-1.0) * µm^(2.0),
        E=400kJ / mol,
        V=18e-6m^3 / mol,
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
        pU = p isa GeoUnit ? p : convert(GeoUnit, p)
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
            NameU, nU, pU, AU, EU, VU, RU, Int8(Apparatus), FT, FE
        )
    end

    function GrainBoundarySliding(Name, n, p, A, E, V, R, Apparatus, FT, FE)
        return GrainBoundarySliding(;
            Name=Name, n=n, p=p, A=A, E=E, V=V, R=R, Apparatus=Apparatus
        )
    end
end

"""
    Transform_GrainBoundarySliding(name)
Transforms units from MPa, kJ etc. to basic units such as Pa, J etc.
"""
function Transform_GrainBoundarySliding(name; kwargs)
    pp_in = GrainBoundarySliding_info[name][1]
    
    # Take optional arguments 
    v_kwargs = values(kwargs)
    val = GeoUnit.(values(v_kwargs))
    
    args = (Name=pp_in.Name, n = pp_in.n, p=pp_in.p, A=pp_in.A, E=pp_in.E, V=pp_in.V, Apparatus=pp_in.Apparatus)
    pp = merge(args, NamedTuple{keys(v_kwargs)}(val))
    
     
    Name = String(collect(pp.Name))
    n = Value(pp.n)
    p = Value(pp.p)
    A_Pa = uconvert(
        Pa^(-NumValue(pp.n)) * m^(-NumValue(p)) / s, Value(pp.A)
    )
    E_J = uconvert(J / mol, Value(pp.E))
    V_m3 = uconvert(m^3 / mol, Value(pp.V))

    Apparatus = pp.Apparatus

    args = (Name=Name, n=n, p=p, A=A_Pa, E=E_J, V=V_m3, Apparatus=Apparatus)
    
    return GrainBoundarySliding(; args...)
end

"""
    s = remove_tensor_correction(s::GrainBoundarySliding)

Removes the tensor correction of the creeplaw, which is useful to compare the implemented creeplaws
with the curves of the original publications, as those publications usually do not transfer their data to tensor format
"""
function remove_tensor_correction(s::GrainBoundarySliding)
    name = String(collect(s.Name))
    return GrainBoundarySliding(;
        Name=name, n=s.n, p=s.p, A=s.A, E=s.E, V=s.V, Apparatus=Invariant
    )
end

function param_info(s::GrainBoundarySliding)
    name = String(collect(s.Name))
    eq = L"\tau_{ij} = 2 \eta  \dot{\varepsilon}_{ij}" #ÄNDERN!!!!!!!!!!!
    if name == ""
        return MaterialParamsInfo(; Equation=eq)
    end

    inf = GrainBoundarySliding_info[name][2]
    return MaterialParamsInfo(;
        Equation=eq, Comment=inf.Comment, BibTex_Reference=inf.BibTex_Reference
    )
end

# Calculation routines for linear viscous rheologies
"""
    compute_εII(a::GrainBoundarySliding, TauII::_T; T::_T, P=one(_T), f=one(_T), d=one(_T), kwargs...)

Returns grain boundary sliding strainrate as a function of 2nd invariant of the stress tensor ``\\tau_{II}`` 
```math
    \\dot{ε}_{II} = A τ_{II}^n d^{p} \\exp \\left(- {{E + PV} \\over RT} \\right)  #ÄNDERN!!!!!!!!
```


"""
@inline function compute_εII(
    a::GrainBoundarySliding, TauII::_T; T=one(precision(a)), P=zero(precision(a)), d=one(precision(a)), args...
) where {_T}
    @unpack_val n, p, A, E, V, R = a
    FT, FE = a.FT, a.FE

    ε = A *
        fastpow(TauII * FT, n) *
        fastpow(d, p) *
        exp(-(E + P * V) / (R * T)) / FE

    return ε 
end

@inline function compute_εII(
    a::GrainBoundarySliding, TauII::Quantity; T=1K, P=0Pa, d=1e-3m, args...
)
    @unpack_units n, p, A, E, V, R = a
    FT, FE = a.FT, a.FE

    ε = A * 
        fastpow(TauII * FT, n) * 
        fastpow(d, p) * 
        exp(-(E + P * V) / (R * T)) / FE

    return ε
end

"""
    compute_εII!(EpsII::AbstractArray{_T,N}, a, TauII::AbstractArray{_T,N}; T, P, f,d,kwargs...)

Computes strainrate as a function of stress
"""
function compute_εII!(
    EpsII::AbstractArray{_T,N},
    a::GrainBoundarySliding,
    TauII::AbstractArray{_T,N};
    T=ones(size(TauII))::AbstractArray{_T,N},
    P=zero(TauII)::AbstractArray{_T,N},
    d=ones(size(TauII))::AbstractArray{_T,N},
    kwargs...,
) where {N,_T}
    @inbounds for i in eachindex(EpsII)
        EpsII[i] = compute_εII(a, TauII[i]; T=T[i], P=P[i], d=d[i])
    end

    return nothing
end

@inline function dεII_dτII(
    a::GrainBoundarySliding, TauII::_T; T=one(precision(a)), P=zero(precision(a)), d=one(precision(a)), args...
) where {_T}
    @unpack_val n, p, A, E, V, R = a
    FT, FE = a.FT, a.FE

    # computed symbolically
    return (A *
            fastpow(d, p) *
            n *
            fastpow(FT * TauII, n) *
            exp(-(E + P * V) / (R * T))) \
            (FE *
            TauII)
end

@inline function dεII_dτII(
    a::GrainBoundarySliding, TauII::Quantity; T=1K, P=0Pa, d=1e-3m, args...
)
    @unpack_units n, p, A, E, V, R = a
    FT, FE = a.FT, a.FE

    #computed symbolically
    return (A *
            fastpow(d, p) *
            n *
            fastpow(FT * TauII, n) *
            exp(-(E + P * V) / (R * T))) \
            (FE *
            TauII)
end

"""
    computeCreepLaw_TauII(EpsII::_T, a::GrainBoundarySliding; T::_T, P=zero(_T), d=one(_T), kwargs...)

Returns grain boundary sliding stress as a function of 2nd invariant of the strain rate 
"""
@inline function compute_τII(
    a::GrainBoundarySliding, EpsII::_T; T=one(precision(a)), P=zero(precision(a)), d=one(precision(a)), kwargs...
) where {_T}
    @unpack_val n, p, A, E, V, R = a
    FT, FE = a.FT, a.FE
    
    n_inv = inv(n)
    
    τ = fastpow(A, -n_inv) *
        fastpow(EpsII * FE, n_inv) *
        fastpow(d, -p * n_inv) *
        fastpow(exp((E + P * V) / (R * T)), n_inv) / FT

    return τ
end

@inline function compute_τII(
    a::GrainBoundarySliding, EpsII::Quantity; T=1K, P=0Pa, d=1m, kwargs...
)
    @unpack_units n, p, A, E, V, R = a
    FT, FE = a.FT, a.FE

    n_inv = inv(n)
    
    τ = fastpow(A, -n_inv) *
        fastpow(EpsII * FE, n_inv) *
        fastpow(d, -p * n_inv) *
        fastpow(exp((E + P * V) / (R * T)), n_inv) / FT

    return τ
end

function compute_τII!(
    TauII::AbstractArray{_T,N},
    a::GrainBoundarySliding,
    EpsII::AbstractArray{_T,N};
    T=ones(size(TauII))::AbstractArray{_T,N},
    P=zero(TauII)::AbstractArray{_T,N},
    d=ones(size(TauII))::AbstractArray{_T,N},
    kwargs...,
) where {N,_T}
    @inbounds for i in eachindex(EpsII)
        TauII[i] = compute_τII(a, EpsII[i]; T=T[i], P=P[i], d=d[i])
    end

    return nothing
end

@inline function dτII_dεII(
    a::GrainBoundarySliding, EpsII::_T; T=one(precision(a)), P=zero(precision(a)), d=one(precision(a)), args...
) where {_T}
    @unpack_val n, p, A, E, V, R = a
    FT, FE = a.FT, a.FE

    n_inv = inv(n)

    # derived symbolically
    return (fastpow(A, -n_inv) *
            fastpow(d, -p * n_inv) *
            fastpow(FE * EpsII, n_inv) *
            fastpow(exp((E + P * V) / (R * T)), n_inv) / 
            (EpsII * FT * n))
end

@inline function dτII_dεII(
    a::GrainBoundarySliding, EpsII::Quantity; T=1K, P=0Pa, d=1e-3m, args...
)
    @unpack_units n, p, A, E, V, R = a
    FT, FE = a.FT, a.FE

    n_inv = inv(n)

    # derived symbolically
    return (fastpow(A, -n_inv) *
            fastpow(d, -p * n_inv) *
            fastpow(FE * EpsII, n_inv) *
            fastpow(exp((E + P * V) / (R * T)), n_inv) / 
            (EpsII * FT * n))
end

# Print info 
function show(io::IO, g::GrainBoundarySliding)
    return print(
        io,
        "GrainBoundarySliding: Name = $(String(collect(g.Name))), n=$(Value(g.n)), p=$(Value(g.p)), A=$(Value(g.A)), E=$(Value(g.E)), V=$(Value(g.V)), FT=$(g.FT), FE=$(g.FE)",
    )
end

# load collection of grain boundary sliding laws
include("Data/GrainBoundarySliding.jl")