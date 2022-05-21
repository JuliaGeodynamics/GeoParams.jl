export  DiffusionCreep,
        SetDiffusionCreep,
        DiffusionCreep_info,
        dεII_dτII,    dτII_dεII,
        compute_εII!, compute_εII,
        compute_τII!, compute_τII
        

#=----Diffusion Creep---
Defines diffusion creep law parameters

n is the power-law exponent
r is the water-fugacity exponent
p is the negative defined grain size exponent (value is set negative in the calculation of σ and ε)
A is a   material specific rheological parameter
E is the activation energy
V is the activation volume
R is the universal gas constant
Apparatus defines the appartus type that shall be recreated (Axial Compression, Simple Shear, Invariant)
=#
@with_kw_noshow mutable struct DiffusionCreep{T,N,U1,U2,U3,U4,U5} <: AbstractCreepLaw{T}
    Name::NTuple{N,Char}        = ""
    n::GeoUnit{T,U1}            = 1.0NoUnits         # power-law exponent
    r::GeoUnit{T,U1}            = 0.0NoUnits         # exponent of water-fugacity
    p::GeoUnit{T,U1}            = -3.0NoUnits         # grain size exponent
    A::GeoUnit{T,U2}            = 1.5MPa^(-n-r)*s^(-1)*m^(-p)    # material specific rheological parameter
    E::GeoUnit{T,U3}            = 500kJ/mol          # activation energy
    V::GeoUnit{T,U4}            = 24e-6m^3/mol       # activation volume
    R::GeoUnit{T,U5}            = 8.3145J/mol/K      # universal gas constant
    Apparatus::Int32            = AxialCompression   # type of experimental apparatus, either AxialCompression, SimpleShear or Invariant
end

DiffusionCreep(args...) = DiffusionCreep(NTuple{length(args[1]), Char}(collect.(args[1])), convert.(GeoUnit,args[2:end-1])..., args[end])

function param_info(s::DiffusionCreep)
    name = String(collect(s.Name))
    eq = L"\tau_{ij} = 2 \eta  \dot{\varepsilon}_{ij}"
    if name == "" 
        return MaterialParamsInfo(Equation=eq)
    end
    inf = DiffusionCreep_info[name][2]
    return MaterialParamsInfo(Equation=eq, Comment=inf.Comment, BibTex_Reference=inf.BibTex_Reference)
end

# Calculation routines for linear viscous rheologies
"""
    compute_εII(a::DiffusionCreep, TauII::_T; T::_T, P=one(_T), f=one(_T), d=one(_T), kwargs...)

Returns diffusion creep strainrate as a function of 2nd invariant of the stress tensor ``\\tau_{II}`` 
```math
    \\dot{ε}_{II} = A τ_{II}^n d^{p} f_{H_2O}^r \\exp \\left(- {{E + PV} \\over RT} \\right)
```


"""
function compute_εII(a::DiffusionCreep, TauII::_T; T::_T, P=zero(_T), f=one(_T), d=one(_T), kwargs...) where _T

    @unpack_val n,r,p,A,E,V,R = a
    
    FT, FE = CorrectionFactor(a)
   
    # ε = A τ^n f^r d^p exp( -(E + pV)/(RT) )
    # τ = τII*FT,  εII=
    return A*fastpow(TauII*FT,n)*fastpow(f,r)*fastpow(d,p)*exp(-(E + P*V)/(R*T))/FE
end
#compute_εII(s::DiffusionCreep, TauII::_T, args) where _T = compute_εII(s, TauII ; args...)

"""
    compute_εII!(EpsII::AbstractArray{_T,N}, a, TauII::AbstractArray{_T,N}; T, P, f,d,kwargs...)

Computes strainrate as a function of stress
"""
function compute_εII!(EpsII::AbstractArray{_T,N}, a::DiffusionCreep, TauII::AbstractArray{_T,N}; 
    T = ones(size(TauII))::AbstractArray{_T,N}, 
    P = zero(TauII)::AbstractArray{_T,N}, 
    f = ones(size(TauII))::AbstractArray{_T,N},
    d = ones(size(TauII))::AbstractArray{_T,N},
    kwargs...)  where {N,_T}
   
    @inbounds for i in eachindex(EpsII)
        EpsII[i] = compute_εII(a, TauII[i], T=T[i], P=P[i], f=f[i], d=d[i])
    end
  
    return nothing
end

"""
    dεII_dτII(a::DiffusionCreep, TauII::_T; T::_T, P=zero(_T), f=one(_T), d=one(_T), kwargs...)

returns the derivative of strainrate versus stress 
"""
function dεII_dτII(a::DiffusionCreep, TauII::_T; T::_T, P=zero(_T), f=one(_T), d=one(_T), kwargs...) where _T
    @unpack_val n,r,p,A,E,V,R = a
    FT, FE = CorrectionFactor(a)
    
    return fastpow(FT*TauII, -1+n)*fastpow(f,r)*fastpow(d,p)*A*FT*n*exp((-E-P*V)/(R*T))*(1/FE)

end
#dεII_dτII(s::DiffusionCreep, TauII::_T, args) where _T = dεII_dτII(s, TauII ; args...)

"""
    computeCreepLaw_TauII(EpsII::_T, a::DiffusionCreep; T::_T, P=zero(_T), f=one(_T), d=one(_T), kwargs...)

Returns dislocation creep stress as a function of 2nd invariant of the strain rate 
"""
function compute_τII(a::DiffusionCreep, EpsII::_T; T::_T, P=zero(_T), f=one(_T), d=one(_T), kwargs...) where _T
    @unpack_val n,r,p,A,E,V,R = a

    FT, FE = CorrectionFactor(a);    

    return fastpow(A,-1/n)*fastpow(EpsII*FE, 1/n)*fastpow(f,-r/n)*fastpow(d,-p/n)*exp((E + P*V)/(n * R*T))/FT
end
#compute_τII(s::DiffusionCreep, EpsII::_T, args) where _T = compute_τII(s, EpsII ; args...)

function compute_τII!(TauII::AbstractArray{_T,N}, a::DiffusionCreep, EpsII::AbstractArray{_T,N}; 
    T = ones(size(TauII))::AbstractArray{_T,N}, 
    P = zero(TauII)::AbstractArray{_T,N}, 
    f = ones(size(TauII))::AbstractArray{_T,N},
    d = ones(size(TauII))::AbstractArray{_T,N},
    kwargs...)  where {N,_T}
   
    @inbounds for i in eachindex(EpsII)
        TauII[i] = compute_τII(a, EpsII[i], T=T[i], P=P[i], f=f[i], d=d[i])
    end
  
    return nothing
end


function dτII_dεII(a::DiffusionCreep, EpsII,; T::_T, P=zero(_T), f=one(_T), d=one(_T), kwargs...) where _T
    @unpack_val n,r,p,A,E,V,R = a
    FT, FE = CorrectionFactor(a)
    return fastpow(FT*EpsII,-1+1/n)*fastpow(f,-r/n)*fastpow(d,-p/n)*fastpow(A,-1/n)*FE*n*exp((E+P*V)/(n*R*T))*(1/(FT*n))
end
#dτII_dεII(s::DiffusionCreep, EpsII::_T, args) where _T = dτII_dεII(s, EpsII ; args...)


# Print info 
function show(io::IO, g::DiffusionCreep)  
    print(io, "DiffusionCreep: Name = $(String(collect(g.Name))), n=$(g.n.val), r=$(g.r.val), p=$(g.p.val), A=$(g.A.val), E=$(g.E.val), V=$(g.V.val), Apparatus=$(g.Apparatus)" )  
end



# Add a list of pre-defined creep laws 
include("DiffusionCreep_Data.jl")