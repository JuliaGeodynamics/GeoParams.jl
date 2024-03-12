export NonLinearPeierlsCreep,
    NonLinearPeierlsCreep_info,
    dεII_dτII,
    remove_tensor_correction,
    Peierls_stress_iterations,
    Transform_NonLinearPeierlsCreep

# NonLinearPeierls Creep ------------------------------------------------
"""
    NonLinearPeierlsCreep(n = 1.0NoUnits, r = 0.0NoUnits, A = 1.5MPa/s, E = 476.0kJ/mol, V = 6e-6m^3/mol, apparatus = AxialCompression )
    
Defines the flow law parameter of a non linear peierls creep law.

The non linear peierls creep law, as used by experimentalists, is given by  
```math  
    \\dot{\\gamma} = A \\sigma_d \\exp{\\left( -\\frac{E}{RT} \\left( 1 - \\left( \\frac{\\sigma_d}{\\sigma_p} \\right) ^o \\right)^q\\right)}
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
julia> x2 = NonLinearPeierlsCreep(n=1)
NonLinearPeierlsCreep: n=1, A=1.5 MPa^-3 s^-1, E=476.0 kJ mol^-1, Apparatus=AxialCompression
```
"""
struct NonLinearPeierlsCreep{T,N,U1,U2,U3,U4,U5} <: AbstractCreepLaw{T}
    Name::NTuple{100,UInt8}
    n::GeoUnit{T,U1} # power-law exponent
    q::GeoUnit{T,U1} # stress relation exponent
    o::GeoUnit{T,U1} # ... (normally called p but used as 'o' since p already exists)
    TauP::GeoUnit{T,U2} # Peierls stress
    A::GeoUnit{T,U3} # material specific rheological parameter
    E::GeoUnit{T,U4} # activation energy
    R::GeoUnit{T,U5} # universal gas constant
    Apparatus::Int8 # type of experimental apparatus, either AxialCompression, SimpleShear or Invariant
    FT::T # type of experimental apparatus, either AxialCompression, SimpleShear or Invariant
    FE::T # type of experimental apparatus, either AxialCompression, SimpleShear or Invariant

    function NonLinearPeierlsCreep(;
        Name="",
        n=2NoUnits,
        q=1.0NoUnits,
        o=0.5NoUnits,
        TauP=8.5e9Pa,
        A=0.57Pa^(-2) * s^(-1),
        E=476e3J / mol,
        R=8.3145J / mol / K,
        Apparatus=AxialCompression,
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
        N = length(Name)
        name = str2tuple(Name)    
        # Create struct
        return new{T,100,U1,U2,U3,U4,U5}(
            name, nU, qU, oU, TauPU, AU, EU, RU, Int8(Apparatus), FT, FE
        )
    end

    function NonLinearPeierlsCreep(Name, n, q, o, TauP, A, E, R, Apparatus, FT, FE)
        return NonLinearPeierlsCreep(;
            Name=Name, n=n, q=q, o=o, TauP=TauP, A=A, E=E, R=R, Apparatus=Apparatus
        )
    end
end

"""
    s = remove_tensor_correction(s::NonLinearPeierlsCreep)

Removes the tensor correction of the creeplaw, which is useful to compare the implemented creeplaws
with the curves of the original publications, as those publications usually do not transfer their data to tensor format
"""
function remove_tensor_correction(s::NonLinearPeierlsCreep)
    name = uint2str(s.Name)

    return NonLinearPeierlsCreep(;
        Name=name, n=s.n, q=s.q, o=s.o, TauP=s.TauP, A=s.A, E=s.E, Apparatus=Invariant
    )
end

function param_info(s::NonLinearPeierlsCreep)
    name = uint2str(s.Name)
    eq = ""
    if name == ""
        return MaterialParamsInfo(; Equation=eq)
    end
    inf = NonLinearPeierlsCreep_info[name][2]
    return MaterialParamsInfo(;
        Equation=eq, Comment=inf.Comment, BibTex_Reference=inf.BibTex_Reference
    )
end

# Calculation routines for linear viscous rheologies
# All inputs must be non-dimensionalized (or converted to consistent units) GeoUnits
@inline function compute_εII(
    a::NonLinearPeierlsCreep, TauII::_T; T=one(precision(a)), args...
) where {_T}
    @unpack_val n, q, o, TauP, A, E, R = a
    FT, FE = a.FT, a.FE

    ε = @pow A * (TauII * FT)^n * exp(-(E / (R * T)) * (1.0 - ((TauII * FT) / TauP)^o)^q) /
        FE

    return ε
end

@inline function compute_εII(a::NonLinearPeierlsCreep, TauII::Quantity; T=1K, args...)
    @unpack_units n, q, o, TauP, A, E, R = a
    FT, FE = a.FT, a.FE

    ε = @pow A * (TauII * FT)^n * exp(-(E / (R * T)) * (1.0 - ((TauII * FT) / TauP)^o)^q) /
        FE

    return ε
end

function compute_εII!(
    EpsII::AbstractArray{_T,N},
    a::NonLinearPeierlsCreep,
    TauII::AbstractArray{_T,N};
    T=ones(size(TauII))::AbstractArray{_T,N},
    args...,
) where {N,_T}
    @inbounds for i in eachindex(EpsII)
        EpsII[i] = compute_εII(a, TauII[i]; T=T[i])
    end

    return nothing
end

function dεII_dτII(a::NonLinearPeierlsCreep, TauII; args...)
    return ForwardDiff.derivative(x -> compute_εII(a, x; args...), TauII)
end

"""
    Peierls_stress_iterations(rheo::NonLinearPeierlsCreep, Tau::Float64, EpsII::Float64, args)

Nonlinear iterations for Peierls creep stress using Newton-Raphson Iterations. Every number needs to be a child of type Real (don't use units here).
The initial stress guess Tau should be at least in the same order of magnitude as the value of TauP is in the used creep law. Example: 1.75e9 is a 
good initial guess for the preexisting "Wet Olivine | Mei et al. (2010)" creep law. Find the sweet spot in the Tau/TauP relation (initial guess/TauP)
if the stress is diverging. Maximum iterations are by default 500 but can be changed as optional argument.
"""
function PeierlsResidual(rheo::NonLinearPeierlsCreep, TauII, EpsII, args)
    return EpsII - compute_εII(rheo, TauII; args...)
end

# implement nonlinear iterations function to iterate until stable stress value
function Peierls_stress_iterations(
    rheo::NonLinearPeierlsCreep, Tau, EpsII, args; max_iter=500
)
    err = 1.0
    dfdtau = 0.0
    i = 0
    while err > 1.0e-6
        i += 1
        if i > max_iter
            print("Stress iterations did not converge.\n")
            break
        end
        Tau_old = Tau

        fTau_n, dfdtau = value_and_partial(x -> PeierlsResidual(rheo, x, EpsII, args), Tau)
        Tau = Tau - (fTau_n / dfdtau)
        err = abs(Tau_old - Tau)

        if err < 1e-6
            println("Converged in $i iterations with err = $err.")
        end
    end
    return Tau
end

# Print info 
function show(io::IO, g::NonLinearPeierlsCreep)
    return print(
        io,
        "NonLinearPeierlsCreep: Name = $(String(collect(g.Name))), n=$(Value(g.n)), q=$(Value(g.q)), o=$(Value(g.o)), TauP=$(Value(g.TauP)), A=$(Value(g.A)), E=$(Value(g.E)), FT=$(g.FT), FE=$(g.FE), Apparatus=$(g.Apparatus)",
    )
end
#-------------------------------------------------------------------------

# load collection of peierls creep laws
include("Data/NonLinearPeierlsCreep.jl")
include("Data_deprecated/NonLinearPeierlsCreep.jl")

using .NonLinearPeierls
export SetNonLinearPeierlsCreep

"""
    SetNonLinearPeierlsCreep["Name of non linear peierls creep law"]
This is a dictionary with pre-defined creep laws    
"""
function SetNonLinearPeierlsCreep(
    name::F;
    n    = nothing,
    q    = nothing,
    o    = nothing,
    A    = nothing,
    E    = nothing,
    TauP = nothing,
) where F 
    kwargs = (; n, q, o, A, E, TauP)
    Transform_NonLinearPeierlsCreep(name, kwargs)
end

function SetNonLinearPeierlsCreep(
    name::F,
    CharDim::GeoUnits{T};
    n    = nothing,
    q    = nothing,
    o    = nothing,
    A    = nothing,
    E    = nothing,
    TauP = nothing,
) where {F, T<:Union{GEO, SI}}
    kwargs = (; n, q, o, A, E, TauP)
    nondimensionalize(Transform_NonLinearPeierlsCreep(name, kwargs), CharDim)
end

"""
    Transforms units from GPa, MPa, kJ etc. to basic units such as Pa, J etc.
"""
Transform_NonLinearPeierlsCreep(name::F) where F = Transform_NonLinearPeierlsCreep(nonlinear_peierls_database(name))
Transform_NonLinearPeierlsCreep(name::F, kwargs::NamedTuple) where F = Transform_NonLinearPeierlsCreep(nonlinear_peierls_database(name), kwargs)

function Transform_NonLinearPeierlsCreep(name::F, CharDim::GeoUnits{U}) where {F, U<:Union{GEO,SI}}
    Transform_NonLinearPeierlsCreep(nonlinear_peierls_database(name), CharDim)
end

function Transform_NonLinearPeierlsCreep(p::AbstractCreepLaw{T}, CharDim::GeoUnits{U}) where {T,U<:Union{GEO,SI}}
    nondimensionalize(Transform_NonLinearPeierlsCreep(p), CharDim)
end

function Transform_NonLinearPeierlsCreep(p::AbstractCreepLaw)
    n = Value(p.n)
    q = Value(p.q)
    o = Value(p.o)
    TauP = uconvert(Pa, Value(p.TauP))
    A_Pa = uconvert(Pa^(unit_power(p.A)) / s, Value(p.A))
    E_J = uconvert(J / mol, Value(p.E))
    Apparatus = p.Apparatus

    # args from database
    args = (Name=p.Name, n=n, q=q, o=o, TauP=TauP, A=A_Pa, E=E_J, Apparatus=Apparatus)

    return NonLinearPeierlsCreep(; args...)
end

function Transform_NonLinearPeierlsCreep(p::AbstractCreepLaw, kwargs::NamedTuple)
    (; n, q, o, A, E, TauP) = kwargs

    n_new    = isnothing(n) ? Value(p.n)    : Value(GeoUnit(n))
    q_new    = isnothing(q) ? Value(p.q)    : Value(GeoUnit(q))
    o_new    = isnothing(o) ? Value(p.o)    : Value(GeoUnit(o))
    A_new    = isnothing(A) ?        p.A    : GeoUnit(A)
    E_new    = isnothing(E) ? Value(p.E)    : Value(GeoUnit(E))
    TauP_new = isnothing(E) ? Value(p.TauP) : Value(GeoUnit(TauP))

    TauP      = uconvert(Pa, TauP_new)
    A_Pa      = uconvert(Pa^(unit_power(A_new)) / s, Value(A_new))
    E_J       = uconvert(J / mol, E_new)
    Apparatus = p.Apparatus

    # args from database
    args      = (Name=p.Name, n=n_new, q=q_new, o=o_new, TauP=TauP_new, A=A_Pa, E=E_J, Apparatus=Apparatus)

    return NonLinearPeierlsCreep(; args...)
end
