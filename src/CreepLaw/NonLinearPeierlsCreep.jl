export NonLinearPeierlsCreep,
    NonLinearPeierlsCreep_info,
    SetNonLinearPeierlsCreep,
    remove_tensor_correction,
    dεII_dτII,
    Peierls_stress_iterations

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

    function NonLinearPeierlsCreep(;
        Name="",
        n=2.0NoUnits,
        q=1.0NoUnits,
        o=0.5NoUnits,
        TauP=8.5e9Pa,
        A=5.7e11MPa^(-2.0) * s^(-1.0),
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

    function NonLinearPeierlsCreep(Name, n, q, o, TauP, A, E, R, Apparatus, FT, FE)
        return NonLinearPeierlsCreep(;
            Name=Name, n=n, q=q, o=o, TauP=TauP, A=A, E=E, R=R, Apparatus=Apparatus
        )
    end
end

"""
    Transforms units from GPa, MPa, kJ etc. to basic units such as Pa, J etc.
"""

function Transform_NonLinearPeierlsCreep(name; kwargs)
    p_in = NonLinearPeierlsCreep_info[name][1]

    # Take optional arguments 
    v_kwargs = values(kwargs)
    val = GeoUnit.(values(v_kwargs))
    
    args = (Name=p_in.Name, n=p_in.n, q=p_in.q, o=p_in.o, TauP=p_in.TauP, A=p_in.A, E=p_in.E, Apparatus=p_in.Apparatus)
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
    
    return NonLinearPeierlsCreep(; args...)
end

"""
    s = remove_tensor_correction(s::NonLinearPeierlsCreep)

Removes the tensor correction of the creeplaw, which is useful to compare the implemented creeplaws
with the curves of the original publications, as those publications usually do not transfer their data to tensor format
"""
function remove_tensor_correction(s::NonLinearPeierlsCreep)
    name = String(collect(s.Name))

    return NonLinearPeierlsCreep(;
        Name=name, n=s.n, q=s.q, o=s.o, TauP=s.TauP, A=s.A, E=s.E, Apparatus=Invariant
    )
end

function param_info(s::NonLinearPeierlsCreep)
    name = String(collect(s.Name))
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

    TauII_FT_n = pow_check(TauII * FT, n)
    TauII_FT_TauP_o = pow_check((TauII * FT) / TauP, o)
    one_minus_TauII_FT_TauP_o_q = pow_check(1.0 - TauII_FT_TauP_o, q)

    ε = A * 
        TauII_FT_n * 
        exp(-(E / (R * T)) * 
        (one_minus_TauII_FT_TauP_o_q)) / 
        FE

    return ε
end

@inline function compute_εII(
    a::NonLinearPeierlsCreep, TauII::Quantity; T=1K, args...
)
    @unpack_units n, q, o, TauP, A, E, R = a
    FT, FE = a.FT, a.FE

    TauII_FT_n = pow_check(TauII * FT, n)
    TauII_FT_TauP_o = pow_check((TauII * FT) / TauP, o)
    one_minus_TauII_FT_TauP_o_q = pow_check(1.0 - TauII_FT_TauP_o, q)

    ε = A * 
        TauII_FT_n * 
        exp(-(E / (R * T)) * 
        (one_minus_TauII_FT_TauP_o_q)) / 
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

dεII_dτII(a::NonLinearPeierlsCreep, TauII; args...) = ForwardDiff.derivative(x -> compute_εII(a, x; args...), TauII)

"""
    Peierls_stress_iterations(rheo::NonLinearPeierlsCreep, Tau::Float64, EpsII::Float64, args)

Nonlinear iterations for Peierls creep stress using Newton-Raphson Iterations. Every number needs to be a child of type Real (don't use units here).
The initial stress guess Tau should be at least in the same order of magnitude as the value of TauP is in the used creep law. Example: 1.75e9 is a 
good initial guess for the preexisting "Wet Olivine | Mei et al. (2010)" creep law. Find the sweet spot in the Tau/TauP relation (initial guess/TauP)
if the stress is diverging. Maximum iterations are by default 500 but can be changed as optional argument.
"""
PeierlsResidual(rheo::NonLinearPeierlsCreep, TauII, EpsII, args) = EpsII - compute_εII(rheo, TauII; args...)

# implement nonlinear iterations function to iterate until stable stress value
function Peierls_stress_iterations(rheo::NonLinearPeierlsCreep, Tau, EpsII, args; max_iter=500)
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
        
        fTau_n, dfdtau = value_and_partial(x->PeierlsResidual(rheo, x, EpsII, args), Tau)
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
