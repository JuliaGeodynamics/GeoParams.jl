"""
    This creates the overall material parameters structure
"""

module MaterialParameters
using Unitful: Energy
using Unitful
using Parameters
using ..Units

import Base.show, Base.convert
using GeoParams: AbstractMaterialParam, AbstractMaterialParamsStruct

export 
    MaterialParams, SetMaterialParams    


# Link the modules with various definitions:
include("./CreepLaw/CreepLaw.jl")
include("./Density/Density.jl")
include("./GravitationalAcceleration/GravitationalAcceleration.jl")
include("./Energy/HeatCapacity.jl")
include("./Energy/Conductivity.jl")


"""
    MaterialParams
    
Structure that holds all material parameters for a given phase

"""
 @with_kw_noshow mutable struct MaterialParams <: AbstractMaterialParamsStruct
    # 
    Name::String         =   ""                  #       Description/name of the phase
    Phase::Int64         =   1;                  #       Number of the phase (optional)
    Nondimensional::Bool =   false;              #       Are all fields non-dimensionalized or not?
    Density              =   nothing             #       Density equation of state
    Gravity              =   nothing             #       Gravitational acceleration (set automatically)
    CreepLaws            =   nothing             #       Creep laws
    Elasticity           =   nothing             #       Elastic parameters
    Plasticity           =   nothing             #       Plasticity
    Conductivity         =   nothing             #       Parameters related to the energy equation 
    HeatCapacity         =   nothing             #       Heat capacity 
    EnergySourceTerms    =   nothing             #       Source terms in energy conservation equation (such as radioactive heat)
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
function SetMaterialParams(; Name::String="", Phase=1,
            Density             =   nothing, 
            Gravity             =   nothing,
            CreepLaws           =   nothing, 
            Elasticity          =   nothing, 
            Plasticity          =   nothing, 
            Conductivity        =   nothing, 
            HeatCapacity        =   nothing, 
            EnergySourceTerms   =   nothing, 
            CharDim             =   nothing)

        # In case density is defined and gravity not, set gravity to default value
        if ~isnothing(Density) & isnothing(Gravity)
            Gravity = GravitationalAcceleration.ConstantGravity();
        end  


        # define struct for phase, while also specifying the maximum number of definitions for every field   
        phase = MaterialParams(Name, Phase, false,
                                     ConvField(Density,             :Density,           maxAllowedFields=1),     
                                     ConvField(Gravity,             :Gravity,           maxAllowedFields=1),       
                                     ConvField(CreepLaws,           :Creeplaws),       
                                     ConvField(Elasticity,          :Elasticity,        maxAllowedFields=1), 
                                     ConvField(Plasticity,          :Plasticity),  
                                     ConvField(Conductivity,        :Conductivity,      maxAllowedFields=1),    
                                     ConvField(HeatCapacity,        :HeatCapacity,      maxAllowedFields=1), 
                                     ConvField(EnergySourceTerms,   :EnergySourceTerms) )

        # [optionally] non-dimensionalize the struct
        if ~isnothing(CharDim) 
            if typeof(CharDim) <: GeoUnits
                Nondimensionalize!(phase, CharDim)
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
    end
    return field
end

# Helper that prints info about each of the material parameters
#  for this to look nice, you need to define a Base.show 
function Print_MaterialParam(io::IO, name::Symbol, Data)
    if ~isnothing(Data) 
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
    println(io, "Phase $(rpad(phase.Phase,2)): $(phase.Name)")
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




end


