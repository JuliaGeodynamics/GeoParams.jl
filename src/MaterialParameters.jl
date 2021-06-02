"""
    This creates the overall material parameters structure
"""

module MaterialParameters
using Unitful
using Parameters
using GeoParams.Units

import Base.show




export 
    MaterialParams,
    CreepLaw_Eps,
    LinearViscous


 @with_kw_noshow mutable struct MaterialParams
    # 
    Name::String        =   ""                  #       Description/name of the phase
    Density             =   nothing             #       Density equation of state
    CreepLaws           =   nothing             #       Creep laws
    Elasticity          =   nothing             #       Elastic parameters
    Plasticity          =   nothing             #       Plasticity
    Conductivity        =   nothing             #       Parameters related to the energy equation 
    HeatCapacity        =   nothing             #        
    EnergySourceTerms   =   nothing             #       Source terms in energy conservation equation
end




# just for testing: define a few creep laws

"""
    Defines a linear viscous creeplaw 

"""
@with_kw_noshow mutable struct LinearViscous
 #   latex_expr  =   "\tau_{ij} = 2*\eta*\dot(\varepsilon}_{ij}"       
    eta         =   1e20            # viscosity
    eta_units   =   Pa*s            # units of viscosity
end



"""
    CreepLaw_EpsII(TauII::Float64, s::LinearViscous)

Returns the strainrate invariant given a stress invariant for a linear viscous rheology

"""
function CreepLaw_Eps(TauII::Float64, s::LinearViscous)
    @unpack eta   = s
    
    return EpsII = TauII/(2.0*eta);
end

"""
    CreepLaw_TauII(EpsII::Float64, s::LinearViscous)

Returns the stress invariant given a strainrate invariant for a linear viscous rheology

"""
function CreepLaw_Eps(EpsII::Float64, s::LinearViscous)
    @unpack eta   = s
    
    return TauII = 2.0*eta*EpsII;
end


struct PowerlawViscous
    eta     ::  Float64
    n       ::  Float64         # powerlaw viscosity
    ε0      ::  Float64;        # 
end






end


