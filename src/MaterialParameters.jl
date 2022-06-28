"""
    This creates the overall material parameters structure
"""

module MaterialParameters
using Unitful: Energy
using Unitful
using Parameters, LaTeXStrings, BibTeX
using ..Units

import Base.show, Base.convert
using GeoParams: AbstractMaterialParam, AbstractMaterialParamsStruct, AbstractPhaseDiagramsStruct

# Define an "empty" Material parameter structure
struct No_MaterialParam{_T} <: AbstractMaterialParam end
No_MaterialParam() = No_MaterialParam{Float64}();

export 
    MaterialParams, SetMaterialParams, No_MaterialParam, MaterialParamsInfo    


"""
    MaterialParamsInfo

Structure that holds information (Equation, Comment, BibTex_Reference) about a given material parameter, which
can be used to create parameter tables, documentation etc.

Usually used in combination with `param_info(the_parameter_of_interest)`

"""
@with_kw_noshow struct MaterialParamsInfo
    Equation::LaTeXString   =   L"" 
    Comment::String         =   ""
    BibTex_Reference::String = ""
end

# Link the modules with various definitions:
include("./PhaseDiagrams/PhaseDiagrams.jl")
include("./CreepLaw/CreepLaw.jl")
include("./Density/Density.jl")
include("./GravitationalAcceleration/GravitationalAcceleration.jl")
include("./Energy/HeatCapacity.jl")
include("./Energy/Conductivity.jl")
include("./Energy/LatentHeat.jl")
include("./Energy/RadioactiveHeat.jl")
include("./Energy/Shearheating.jl")
include("./SeismicVelocity/SeismicVelocity.jl")

using .Density: AbstractDensity

"""
    MaterialParams
    
Structure that holds all material parameters for a given phase

"""
 @with_kw_noshow struct MaterialParams{ N,
                                        Vdensity  <: Tuple, 
                                        Vgravity  <: Tuple,
                                        Vcreep    <: Tuple,
                                        Velastic  <: Tuple,
                                        Vplastic  <: Tuple,
                                        Vcond     <: Tuple,
                                        Vheatc    <: Tuple,
                                        Vensource <: Tuple,
                                        Vmelting  <: Tuple,
                                        Vseismvel <: Tuple } <: AbstractMaterialParamsStruct
    Name::NTuple{N,Char}                            #       The name is encoded as a NTuple{Char} (to make it isbits and the whole MaterialParams isbits as well; required to use this on the GPU)
    Phase::Int64                 =   1;             #       Number of the phase (optional)
    Nondimensional::Bool         =   false;         #       Are all fields non-dimensionalized or not?
    Density::Vdensity            =   ()             #       Density equation of state
    Gravity::Vgravity            =   ()             #       Gravitational acceleration (set automatically)
    CreepLaws::Vcreep            =   ()             #       Creep laws
    Elasticity::Velastic         =   ()             #       Elastic parameters
    Plasticity::Vplastic         =   ()             #       Plasticity
    Conductivity::Vcond          =   ()             #       Parameters related to the energy equation 
    HeatCapacity::Vheatc         =   ()             #       Heat capacity 
    EnergySourceTerms::Vensource =   ()             #       Source terms in energy conservation equation (such as radioactive heat)
    Melting::Vmelting            =   ()             #       Melting model
    SeismicVelocity::Vseismvel   =   ()             #       Seismic velocity
end

"""
    SetMaterialParams(; Name::String="", Phase::Int64=1,
                        Density             =   nothing, 
                        Gravity             =   nothing,
                        CreepLaws           =   nothing, 
                        Elasticity          =   nothing, 
                        Plasticity          =   nothing, 
                        Conductivity        =   nothing, 
                        HeatCapacity        =   nothing, 
                        EnergySourceTerms   =   nothing, 
                        Melting             =   nothing,
                        SeismicVelocity     =   nothing,
                        CharDim::GeoUnits   =   nothing)

Sets material parameters for a given phase. 

If `CharDim` is specified the input parameters are non-dimensionalized.  
Note that if `Density` is specified, we also set `Gravity` even if not explicitly listed
    
# Examples

Define two viscous creep laws & constant density:
```julia-repl
julia> Phase = SetMaterialParams(Name="Viscous Matrix",
                       Density   = ConstantDensity(),
                       CreepLaws = (PowerlawViscous(), LinearViscous(η=1e21Pa*s)))
Phase 1 : Viscous Matrix
        | [dimensional units]
        | 
        |-- Density           : Constant density: ρ=2900 kg m⁻³ 
        |-- Gravity           : Gravitational acceleration: g=9.81 m s⁻² 
        |-- CreepLaws         : Powerlaw viscosity: η0=1.0e18 Pa s, n=2.0, ε0=1.0e-15 s⁻¹  
        |                       Linear viscosity: η=1.0e21 Pa s
```

Define two viscous creep laws & P/T dependent density and nondimensionalize
```julia-repl
julia> CharUnits_GEO   =   GEO_units(viscosity=1e19, length=1000km);
julia> Phase = SetMaterialParams(Name="Viscous Matrix", Phase=33,
                              Density   = PT_Density(),
                              CreepLaws = (PowerlawViscous(n=3), LinearViscous(η=1e23Pa*s)),
                              CharDim   = CharUnits_GEO)
Phase 33: Viscous Matrix
        | [non-dimensional units]
        | 
        |-- Density           : P/T-dependent density: ρ0=2.9e-16, α=0.038194500000000006, β=0.01, T0=0.21454659702313156, P0=0.0 
        |-- Gravity           : Gravitational acceleration: g=9.810000000000002e18 
        |-- CreepLaws         : Powerlaw viscosity: η0=0.1, n=3, ε0=0.001  
        |                       Linear viscosity: η=10000.0 
```

You can also create an array that holds several parameters:
```julia-repl
julia> MatParam        =   Array{MaterialParams, 1}(undef, 2);
julia> Phase           =   1;
julia> MatParam[Phase] =   SetMaterialParams(Name="Upper Crust", Phase=Phase,
                            CreepLaws= (PowerlawViscous(), LinearViscous(η=1e23Pa*s)),
                            Density  = ConstantDensity(ρ=2900kg/m^3));
julia> Phase           =   2;
julia> MatParam[Phase] =   SetMaterialParams(Name="Lower Crust", Phase=Phase,
                            CreepLaws= (PowerlawViscous(n=5), LinearViscous(η=1e21Pa*s)),
                            Density  = PT_Density(ρ0=3000kg/m^3));
julia> MatParam
2-element Vector{MaterialParams}:
 Phase 1 : Upper Crust
    | [dimensional units]
    | 
    |-- Density           : Constant density: ρ=2900 kg m⁻³ 
    |-- Gravity           : Gravitational acceleration: g=9.81 m s⁻² 
    |-- CreepLaws         : Powerlaw viscosity: η0=1.0e18 Pa s, n=2.0, ε0=1.0e-15 s⁻¹  
    |                       Linear viscosity: η=1.0e23 Pa s 
                            
 Phase 2 : Lower Crust
    | [dimensional units]
    | 
    |-- Density           : P/T-dependent density: ρ0=3000 kg m⁻³, α=3.0e-5 K⁻¹, β=1.0e-9 Pa⁻¹, T0=0 °C, P0=0 MPa 
    |-- Gravity           : Gravitational acceleration: g=9.81 m s⁻² 
    |-- CreepLaws         : Powerlaw viscosity: η0=1.0e18 Pa s, n=5, ε0=1.0e-15 s⁻¹  
    |                       Linear viscosity: η=1.0e21 Pa s 
```


"""
function SetMaterialParams(; 
            Name::String        =   "",         # this makes the struct !isbits(); as that sucks for portability 
            Phase               =   1,
            Density             =   nothing, 
            Gravity             =   nothing,
            CreepLaws           =   nothing, 
            Elasticity          =   nothing, 
            Plasticity          =   nothing, 
            Conductivity        =   nothing, 
            HeatCapacity        =   nothing, 
            EnergySourceTerms   =   nothing, 
            Melting             =   nothing,
            SeismicVelocity     =   nothing,
            CharDim             =   nothing)

        # In case density is defined and gravity not, set gravity to default value
        if ~isnothing(Density) & isnothing(Gravity)
            Gravity = GravitationalAcceleration.ConstantGravity();
        end  


        # define struct for phase, while also specifying the maximum number of definitions for every field   
        #=
        phase = MaterialParams(Name, Phase, false,
                                     ConvField(Density,             :Density,           maxAllowedFields=1),     
                                     ConvField(Gravity,             :Gravity,           maxAllowedFields=1),       
                                     ConvField(CreepLaws,           :Creeplaws),       
                                     ConvField(Elasticity,          :Elasticity,        maxAllowedFields=1), 
                                     ConvField(Plasticity,          :Plasticity),  
                                     ConvField(Conductivity,        :Conductivity,      maxAllowedFields=1),    
                                     ConvField(HeatCapacity,        :HeatCapacity,      maxAllowedFields=1), 
                                     ConvField(EnergySourceTerms,   :EnergySourceTerms), 
                                     ConvField(Melting,             :Melting,           maxAllowedFields=1),
                                     ConvField(SeismicVelocity,     :SeismicVelocity,   maxAllowedFields=1)
                                     ) 
        =#                        
        phase = MaterialParams(NTuple{length(Name), Char}(collect.(Name)),  
                                     Phase, 
                                     false,
                                     ConvField(Density,             :Density,           maxAllowedFields=1),     
                                     ConvField(Gravity,             :Gravity,           maxAllowedFields=1),       
                                     ConvField(CreepLaws,           :Creeplaws),       
                                     ConvField(Elasticity,          :Elasticity,        maxAllowedFields=1), 
                                     ConvField(Plasticity,          :Plasticity),  
                                     ConvField(Conductivity,        :Conductivity,      maxAllowedFields=1),    
                                     ConvField(HeatCapacity,        :HeatCapacity,      maxAllowedFields=1), 
                                     ConvField(EnergySourceTerms,   :EnergySourceTerms), 
                                     ConvField(Melting,             :Melting,           maxAllowedFields=1),
                                     ConvField(SeismicVelocity,     :SeismicVelocity,   maxAllowedFields=1)
                                     ) 

        # [optionally] non-dimensionalize the struct
        if ~isnothing(CharDim) 
            if typeof(CharDim) <: GeoUnits
                phase = nondimensionalize(phase, CharDim)
            else
                error("CharDim should be of type GeoUnits")
            end

        end

        return phase
end

# Helper function that converts a field to a Tuple, provided it is not nothing
# This also checks for the maximum allowed number of definitions 
# (some rheological phases may allow for an arbitrary combination per phase; others like density EoS not) 
function ConvField(field, fieldname::Symbol; maxAllowedFields=1e6)
    if ~isnothing(field)
        if typeof(field) <: AbstractMaterialParam
            field = (field, )       # transform to tuple
        end
        if typeof(field[1]) <: AbstractMaterialParam
            if length(field)>maxAllowedFields
                error("Maximum $(maxAllowedFields) field allowed for: $fieldname")
            end
        end
    else
        field = ()  # empty tuple
    end
    return field
end

# Helper that prints info about each of the material parameters
#  for this to look nice, you need to define a Base.show 
function Print_MaterialParam(io::IO, name::Symbol, Data)
    if length(Data)>0 
        if typeof(Data[1]) <: AbstractMaterialParam
            print(io, "        |-- $(rpad(name,18)):")
            for i=1:length(Data)
                if i==1
                    print(io, " $(Data[1]) \n")
                else
                    print(io, "        |  $(rpad("     ",18))   $(Data[i]) \n")
                end
            end
        end
    end
end


# Specify how the printing of the MaterialParam structure is done
function Base.show(io::IO, phase::MaterialParams)
    println(io, "Phase $(rpad(phase.Phase,2)): $(String(collect(phase.Name)))")
    if phase.Nondimensional
        println(io,"        | [non-dimensional units]")
    else
        println(io,"        | [dimensional units]")

    end
    println(io,"        | ")
        
    for param in fieldnames(typeof(phase))
        Print_MaterialParam(io, param, getfield(phase, param))
    end
    
end

# Slightly nicer printout in case we have a tuple with material parameters 
function Base.show(io::IO, phase_tuple::NTuple{N,MaterialParams}) where N
    for i=1:N
        Base.show(io, phase_tuple[i])
    end
    return
end

# Automatically fill tuples with No_MaterialParam given a length n
function fill_tup(v::NTuple{N,Tuple{Vararg{AbstractMaterialParam}}}, n) where N
    ntuple(i->ntuple(j-> j <= length(v[i]) ? v[i][j] : No_MaterialParam(), Val(n)),Val(N))
end

end


