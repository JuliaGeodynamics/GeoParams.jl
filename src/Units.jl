"""
    This provides units and creates a non-dimensionalization object
"""
module Units
using Unitful
using Parameters

import Base.show

# Define a few additional units that are useful in geodynamics 
@unit    Myrs  "Myrs"   MillionYears    1000000u"yr"    false


function __init__()
    #merge!(Unitful.basefactors, localunits)
    Unitful.register(Units)
end


Unitful.register(Units)

# Define a number of useful units
const km    = u"km"
const m     = u"m"
const cm    = u"cm"
const Myrs  = u"Myrs"
const yr    = u"yr"
const s     = u"s"
const Pa    = u"Pa"
const MPa   = u"MPa"
const Pas   = u"Pa*s"
const K     = u"K"
const C     = u"°C"


export 
    km, m, cm, Mtrs, yr, s, MPa, Pa, Pas, K, C, 
    GeoUnits, GEO_units, SI_units, NO_units, AbstractGeoUnits, Nondimensionalize

"""
AbstractGeoUnits

Abstract supertype for geo units.
"""
abstract type AbstractGeoUnits{TYPE} end

abstract type AbstractUnitType end

struct GEO <: AbstractUnitType end
struct SI  <: AbstractUnitType end
struct NONE<: AbstractUnitType end

#abstract type GeoUnits{TYPE} <: AbstractGeoUnits{TYPE} end


@with_kw_noshow struct GeoUnits{TYPE} 
    # Selectable input parameters
    temperature     =   1               #   Characteristic temperature  [C or K]
    length          =   1               #   Characteristic length unit  [km, m or -]
    stress          =   1               #   Characteristic stress unit  [MPa, Pa or -]
    time            =   1               #   Characteristic time unit    [Myrs, s pr -]
    viscosity       =   1

    # primary characteristic units
    K               =   1;              # temperature in SI units for material parameter scaling
    s               =   1;              # time in SI units for material parameter scaling
    m               =   1;              # length in SI units for material parameter scaling
    Pa              =   1;              # stress in SI units 
    kg              =   Pa*m^2/(m/s^2)  # compute mass from pascal
    
    # Derived units
    N               =   kg*m/s^2        # Newton
    J               =   N*m             # Joule
    W               =   J/s             # Watt
    area            =   m^2             # area   in SI units for material parameter scaling
    volume          =   m^3             # volume
    velocity        =   m/s             
    density         =   kg/m^3
    acceleration    =   m/s^2
    force           =   kg*m/s^2
    
    # Helpful
    SecYear         =   3600*24*365.25
    Myrs            =   1e6
    cmYear          =   SecYear*100     # to transfer m/s -> cm/yr
    
end

"""
    Specify the characteristic values using GEO units (more convenient for geodynamic simulations)
    This is as SI units, but has length in km, time in Myrs, temperature in C, stress in MPa
    
    Usage:
        CharUnits = GEO_units(length=1000km, temperature=1000C, stress=10MPa, viscosity=1e20Pas)
    
    CharUnits contains the units that can be used to nondimensionalize or dimensionalize all parameters    
    in the simulations. It contains all basic SI units

    Examples:
        CharUnits.Pa  - Characteristic value in units of Pascal
        CharUnits.s   - Characteristic value in units of seconds
        
"""
function GEO_units(;length=1000km, temperature=1000C, stress=10MPa, viscosity=1e20Pas)
    
    if unit(temperature)==NoUnits;  temperature = temperature*C;        end
    if unit(length)==NoUnits;       length      = length*u"km";         end
    if unit(stress)==NoUnits;       stress      = stress*u"MPa";        end
    if unit(viscosity)==NoUnits;    viscosity   = viscosity*u"Pa*s";    end
    
    T       =   uconvert(C,      temperature) 
    Le      =   uconvert(km,     length);
    Sigma   =   uconvert(MPa,    stress)
    Eta     =   uconvert(Pas,    viscosity)
    
    T_SI    =   uconvert(K,      T);
    Le_SI   =   uconvert(m,      Le);
    Sigma_SI=   uconvert(Pa,     Sigma)
    Time_SI =   Eta/Sigma_SI;
    t       =   uconvert(Myrs,   Time_SI)

    GeoUnits{GEO}(length=Le, temperature=T,      stress=Sigma,  viscosity=Eta, time=t,
                  m=Le_SI,   K=T_SI,   Pa=Sigma_SI,  s=Time_SI)
end


"""
    Specify the characteristic values using SI units 
    
    Usage:
        CharUnits = SI_units(length=1000m, temperature=1000K, stress=10Pa, viscosity=1e20)
    
    CharUnits contains the units that can be used to nondimensionalize or dimensionalize all parameters    
    in the simulations. It contains all basic SI units.

    Examples:
        CharUnits.Pa  - Units of Pascal
        CharUnits.s   - Units of seconds
        

"""
function SI_units(;length=1000m, temperature=1000K, stress=10Pa, viscosity=1e20Pas)
    
    if unit(temperature)==NoUnits;  temperature = temperature*K;     end
    if unit(length)==NoUnits;       length      = length*u"m";         end
    if unit(stress)==NoUnits;       stress      = stress*u"Pa";        end
    if unit(viscosity)==NoUnits;    viscosity   = viscosity*u"Pa*s";    end
    
    T       =   uconvert(K,     temperature) 
    Le      =   uconvert(m,     length);
    Sigma   =   uconvert(Pa,    stress)
    Eta     =   uconvert(Pas,   viscosity)
    
    T_SI    =   uconvert(K,      T);
    Le_SI   =   uconvert(m,      Le);
    Sigma_SI=   uconvert(Pa,     Sigma)
    Time_SI =   Eta/Sigma_SI;
    t       =   uconvert(s,      Time_SI)

    GeoUnits{SI}(length=Le, temperature=T,      stress=Sigma,  viscosity=Eta, time=t,
                 m=Le_SI,   K=T_SI,   Pa=Sigma_SI,  s=Time_SI)
end

"""
    In case we employ nondimension al units throughout
    
    Usage:
        CharUnits = NO_units(length=1000km, temperature=1000C, stress=10MPa, viscosity=1e20Pas)
    
    CharUnits contains the units that can be used to nondimensionalize or dimensionalize all parameters    
    in the simulations. It contains all basic SI units

    Examples:
        CharUnits.Pa  - Units of Pascal
        CharUnits.s   - Units of seconds
        

"""
function NO_units(;length=1000, temperature=1000, stress=10, viscosity=1e20)
    
    if unit(temperature)!=NoUnits;  error("temperature should not have units")    end
    if unit(length)!=NoUnits;       error("length should not have units")    end
    if unit(stress)!=NoUnits;       error("stress should not have units")    end
    if unit(viscosity)!=NoUnits;    error("viscosity should not have units")    end
    
    T       =   temperature
    Le      =   length;
    Sigma   =   stress
    Eta     =   viscosity
    Time    =   Eta/Sigma;

    GeoUnits{NONE}(length=Le, temperature=T, stress=Sigma, viscosity=Eta, time=Time,
                 m=Le, K=T, Pa=Sigma, s=Time)
end

"""
    This nondimensionalizes a parameter using the characteristic values 
    listed in the GeoUnits structure

    Usage:
        param_ND = Nondimensionalize(param, g::GeoUnits{TYPE}) where {TYPE}

    Input:
        param       -   a (dimensional) parameter 
        g           -   the dimensionalization structure 

    Output:
        param_ND    -   the nondimensional parameter    


    Example:
        CharUnits   =   GEO_units();
        v           =   100m/s 
        v_ND        =   Nondimensionalize(v, CharUnits) 

"""
function Nondimensionalize(param, g::GeoUnits{TYPE}) where {TYPE}
    # to be finished; idea is to automatically determine the units of param, and use the same units     

end

# Define a view for the GEO_Units structure
function show(io::IO, g::GeoUnits{TYPE})  where {TYPE}
    print(io, "Employing $TYPE units \n",
              "Characteristic values: \n",  
              "         length:      $(g.length)\n",
              "         time:        $(round(ustrip(g.time),digits=4)) $(unit(g.time))\n",
              "         stress:      $(g.stress)\n",
              "         temperature: $(Float64(g.temperature))\n")
end


end