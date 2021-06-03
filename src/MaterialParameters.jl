"""
    This creates the overall material parameters structure
"""

module MaterialParameters
using Unitful: Energy
using Unitful
using Parameters
using ..Units

import Base.show
using GeoParams: AbstractMaterialParam

export 
    MaterialParams, SetMaterialParams    

"""
    MaterialParams
    
Structure that holds all material parameters for a given phase

"""
 @with_kw_noshow mutable struct MaterialParams
    # 
    Name::String        =   ""                  #       Description/name of the phase
    Density             =   nothing             #       Density equation of state
    CreepLaws::Tuple    =   nothing             #       Creep laws
    Elasticity          =   nothing             #       Elastic parameters
    Plasticity          =   nothing             #       Plasticity
    Conductivity        =   nothing             #       Parameters related to the energy equation 
    HeatCapacity        =   nothing             #        
    EnergySourceTerms   =   nothing             #       Source terms in energy conservation equation
end

"""
    SetMaterialParams(; Name::String="", 
                        Density=nothing, 
                        CreepLaws=nothing, 
                        Elasticity=nothing, 
                        Plasticity=nothing, 
                        Conductivity=nothing, 
                        HeatCapacity=nothing, 
                        EnergySourceTerms=nothing, 
                        CharDim::GeoUnits=nothing)

Sets material parameters for a given phase. 

if `CharDim` is specified the input parameters are non-dimensionalized   
    
# Examples

Define two viscous creep laws & non-dimensionalize
```julia-repl
julia> Phase = SetMaterialParams(Name="Viscous Matrix",
                CreepLaws= (PowerlawViscous(), LinearViscous(eta=1e21Pa*s)),
                CharDim=CharUnits_GEO);
```

"""
function SetMaterialParams(; Name::String="", 
            Density=nothing, 
            CreepLaws=nothing, 
            Elasticity=nothing, 
            Plasticity=nothing, 
            Conductivity=nothing, 
            HeatCapacity=nothing, 
            EnergySourceTerms=nothing, 
            CharDim=nothing)

        # define struct for phase   
        phase = MaterialParams(Name,Density, CreepLaws, Elasticity, Plasticity, Conductivity, HeatCapacity, EnergySourceTerms)

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


# Helper that prints info about each of the material parameters
#  for this to look nice, you need to define a Base.show 
function Print_MaterialParam(io::IO, name::Symbol, Data)
    if ~isnothing(Data) 
        if typeof(Data[1]) <: AbstractMaterialParam
            print(io, "     |-- $(name):")
            for i=1:length(Data)
                if i==1
                    print(io, " $(Data[1]) \n")
                else
                    print(io, "     |              $(Data[i]) \n")
                end
            end
        end
    end
end

function Base.show(io::IO, phase::MaterialParams)
    println(io, "Phase: $(phase.Name)")
    println(io, "     | ")
    
    for param in fieldnames(typeof(phase))
        Print_MaterialParam(io, param, getfield(phase, param))
    end
    
end



"""
    Nondimensionalize!(phase_mat::MaterialParams, g::GeoUnits{TYPE})

Nondimensionalizes all fields within the Material Parameters structure that contain material parameters
"""
function Nondimensionalize!(phase_mat::MaterialParams, g::GeoUnits{TYPE}) where {TYPE} 

    for param in fieldnames(typeof(phase_mat))
        fld = getfield(phase_mat, param)
        if ~isnothing(fld)
            for i=1:length(fld)
                if typeof(fld[i]) <: AbstractMaterialParam
                    Units.Nondimensionalize!(fld[i],g)
                end
            end
        end
    end
end

"""
    Dimensionalize!(phase_mat::MaterialParams, g::GeoUnits{TYPE})

Dimensionalizes all fields within the Material Parameters structure that contain material parameters
"""

function Dimensionalize!(phase_mat::MaterialParams, g::GeoUnits{TYPE}) where {TYPE} 

    for param in fieldnames(typeof(phase_mat))
        fld = getfield(phase_mat, param)
        if ~isnothing(fld)
            for i=1:length(fld)
                if typeof(fld[i]) <: AbstractMaterialParam
                    Units.Dimensionalize!(fld[i],g)
                end
            end
        end
    end
end

# Link the modules with various definitions:
include("./CreepLaw/CreepLaw.jl")


end


