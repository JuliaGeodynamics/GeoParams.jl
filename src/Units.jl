"""
    This provides units and creates a non-dimensionalization object
"""
module Units
using Unitful
import Unitful: superscript
using Parameters

import Base.show

# Define additional units that are useful in geodynamics 
@unit    Myrs  "Myrs"   MillionYears    1000000u"yr"    false

function __init__()
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
const kg    = u"kg"
const Pa    = u"Pa"
const MPa   = u"MPa"
const Pas   = u"Pa*s"
const K     = u"K"
const C     = u"°C"
const mol   = u"mol"  


export 
    km, m, cm, Mtrs, yr, s, MPa, Pa, Pas, K, C, kg, mol, 
    GeoUnit, GeoUnits, GEO_units, SI_units, NO_units, AbstractGeoUnits, 
    Nondimensionalize, Nondimensionalize!, Dimensionalize, Dimensionalize!,
    superscript, upreferred

"""
AbstractGeoUnits

Abstract supertype for geo units.
"""
abstract type AbstractGeoUnits{TYPE} end

abstract type AbstractUnitType end

struct GEO <: AbstractUnitType end
struct SI  <: AbstractUnitType end
struct NONE<: AbstractUnitType end

"""
    Structure that holds a GeoUnit parameter and their dimensions

    Having that is useful, as non-dimensionalization removes the units from a number
    and we thus no longer know how to transfer it back to the correct units.

"""
mutable struct GeoUnit 
    val                 # the actual value ()
    unit :: Unitful.FreeUnits
end


GeoUnit(v::Unitful.Quantity)            =   GeoUnit(v, unit(v))     # store the units 
GeoUnit(v::Number)                      =   GeoUnit(v, NoUnits)     # in case we just have a number with no units
GeoUnit(v::Array)                       =   GeoUnit(v, NoUnits)     # in case we just have a number with no units
GeoUnit(v::Array{Unitful.Quantity})     =   GeoUnit(v, unit.(v))    # in case we just have a number with no units

Base.convert(::Type{Float64}, v::GeoUnit) = v.val

# define a few basic routines so we can easily operate with GeoUnits
Base.show(io::IO, x::GeoUnit)  = println(x.val)


Base.:*(x::GeoUnit, y::Number)  = x.val*y
Base.:+(x::GeoUnit, y::Number)  = x.val+y
Base.:/(x::GeoUnit, y::Number)  = x.val/y
Base.:-(x::GeoUnit, y::Number)  = x.val-y

Base.:*(x::GeoUnit, y::Unitful.Quantity)  = GeoUnit(x.val*y, x.unit)
Base.:+(x::GeoUnit, y::Unitful.Quantity)  = GeoUnit(x.val+y, x.unit)
Base.:/(x::GeoUnit, y::Unitful.Quantity)  = GeoUnit(x.val/y, x.unit)
Base.:-(x::GeoUnit, y::Unitful.Quantity)  = GeoUnit(x.val-y, x.unit)

Base.:*(x::GeoUnit, y::Array)   = GeoUnit(x.val*y, x.unit)
Base.:/(x::GeoUnit, y::Array)   = GeoUnit(x.val/y, x.unit)
Base.:+(x::GeoUnit, y::Array)   = GeoUnit(x.val+y, x.unit)
Base.:-(x::GeoUnit, y::Array)   = GeoUnit(x.val-y, x.unit)

Base.getindex(x::GeoUnit, i::Int64, j::Int64) = x.val[i,j]
Base.getindex(x::GeoUnit, i::Int64) = x.val[i]

Base.setindex!(x::GeoUnit, v::Any, i::Int64, j::Int64) = x.val[i,j] = v
Base.setindex!(x::GeoUnit, v::Any, i::Int64) = x.val[i] = v

"""
    GeoUnits

    Structure that holds parameters used for non-dimensionalization
"""
@with_kw_noshow struct GeoUnits{TYPE} 
    # Selectable input parameters
    temperature     =   1               #   Characteristic temperature  [C or K]
    length          =   1               #   Characteristic length unit  [km, m or -]
    stress          =   1               #   Characteristic stress unit  [MPa, Pa or -]
    time            =   1               #   Characteristic time unit    [Myrs, s pr -]
    viscosity       =   1

    # primary characteristic units
    K               =   1;                      # temperature in SI units for material parameter scaling
    s               =   1;                      # time in SI units for material parameter scaling
    m               =   1;                      # length in SI units for material parameter scaling
    Pa              =   1;                      # stress in SI units 
    kg              =   upreferred(Pa*m*s^2)    # compute mass from pascal. Note: this may result in very large values
    
    # Main SI units (used later for nondimensionalization)
    Length          =   m
    Mass            =   kg
    Time            =   s    
    Temperature     =   K
    Amount          =   1mol
    Second          =   s
    

    # Not defined, as they are not common in geodynamics:
    #Current
    #Luminosity
    
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
    strainrate      =   1/s
    
    # Helpful
    SecYear         =   3600*24*365.25
    Myrs            =   1e6
    cmYear          =   SecYear*100     # to transfer m/s -> cm/yr
    
end

"""
    GEO_units(;length=1000km, temperature=1000C, stress=10MPa, viscosity=1e20Pas)

Creates a non-dimensionalization object using GEO units.

GEO units implies that upon dimensionalization, `time` will be in `Myrs`, `length` in `km`, stress in `MPa`, etc.
which is more convenient for typical geodynamic simulations than SI units
The characteristic values given as input can be in arbitrary units (`km` or `m`), provided the unit is specified.

# Examples:
```julia-repl
julia> CharUnits = GEO_units()
Employing GeoParams.Units.GEO units 
Characteristic values: 
         length:      1000 km
         time:        0.3169 Myrs
         stress:      10 MPa
         temperature: 1000.0 °C
julia> CharUnits.velocity
1.0e-7 m s⁻¹
```
If we instead have a crustal-scale simulation, it is likely more appropriate to use a different characteristic `length`:
```julia-repl
julia> CharUnits = GEO_units(length=10km)
Employing GeoParams.Units.GEO units 
Characteristic values: 
         length:      10 km
         time:        0.3169 Myrs
         stress:      10 MPa
         temperature: 1000.0 °C
```
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
    CharUnits = SI_units(length=1000m, temperature=1000K, stress=10Pa, viscosity=1e20)

Specify the characteristic values using SI units 

# Examples:
```julia-repl
julia> CharUnits = SI_units(length=1000m)
Employing GeoParams.Units.SI units 
Characteristic values: 
         length:      1000 m
         time:        1.0e19 s
         stress:      10 Pa
         temperature: 1000.0 K
```
Note that the same can be achieved if the input is given in `km`:
```julia-repl
julia> CharUnits = SI_units(length=1km)
```
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
    CharUnits = NO_units(length=1, temperature=1, stress=1, viscosity=1)
   
Specify the characteristic values in non-dimensional units

# Examples:
```julia-repl
julia> using GeoParams;
julia> CharUnits = NO_units()
Employing GeoParams.Units.SI units 
Characteristic values: 
         length:      1
         time:        1.0 
         stress:      1
         temperature: 1.0
```
"""
function NO_units(;length=1, temperature=1, stress=1, viscosity=1)
    
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
    Nondimensionalize(param, CharUnits::GeoUnits{TYPE})

Nondimensionalizes `param` using the characteristic values specified in `CharUnits`

# Example 1
```julia-repl
julia> using GeoParams;
julia> CharUnits =   GEO_units();
julia> v         =   3cm/yr
3 cm yr⁻¹ 
julia> v_ND      =   Nondimensionalize(v, CharUnits) 
0.009506426344208684
```
# Example 2
In geodynamics one sometimes encounters more funky units
```julia-repl
julia> CharUnits =   GEO_units();
julia> A         =   6.3e-2MPa^-3.05*s^-1
0.063 MPa⁻³·⁰⁵ s⁻¹
julia> A_ND      =   Nondimensionalize(A, CharUnits) 
7.068716262102384e14
```

In case you are interested to see how the units of `A` look like in different units, use this function from the [Unitful](https://github.com/PainterQubits/Unitful.jl) package:
```julia-repl
julia> uconvert(u"Pa^-3.05*s^-1",A) 
3.157479571851836e-20 Pa⁻³·⁰⁵
```
and to see it decomposed in the basic `SI` units of length, mass and time:
```julia-repl
julia> upreferred(A)
3.1574795718518295e-20 m³·⁰⁵ s⁵·¹ kg⁻³·⁰⁵
```
"""
function Nondimensionalize(param, g::GeoUnits{TYPE}) where {TYPE}

    if unit(param)!=NoUnits
        dim         =   Unitful.dimension(param);                   # Basic SI units
        char_val    =   1.0;
        foreach((typeof(dim).parameters[1])) do y
            val = upreferred(getproperty(g, Unitful.name(y)))       # Retrieve the characteristic value from structure g
            pow = Float64(y.power)                                  # power by which it should be multiplied   
            char_val *= val^pow                                     # multiply characteristic value
        end
        param_ND = upreferred.(param)/char_val
    else
        param_ND = param # The parameter has no units, so there is no way to determine how to nondimensionize it 
    end
    return param_ND
end

"""
    Nondimensionalize!(param::GeoUnit, CharUnits::GeoUnits{TYPE})

Nondimensionalizes `param` (given as GeoUnit) using the characteristic values specified in `CharUnits` in-place

# Example 1
```julia-repl
julia> using GeoParams;
julia> CharUnits =   GEO_units();
julia> v         =   GeoUnit(3cm/yr)
3 cm yr⁻¹ 
julia> Nondimensionalize!(v, CharUnits) 
0.009506426344208684
```
# Example 2
```julia-repl
julia> CharUnits =   GEO_units();
julia> A         =   GeoUnit(6.3e-2MPa^-3.05*s^-1)
0.063 MPa⁻³·⁰⁵ s⁻¹
julia> A_ND      =   Nondimensionalize(A, CharUnits) 
7.068716262102384e14
```
"""
function Nondimensionalize!(param::GeoUnit, g::GeoUnits{TYPE}) where {TYPE}

    if unit.(param.val)!=NoUnits
        dim         =   Unitful.dimension(param.val[1]);                   # Basic SI units
        char_val    =   1.0;
        foreach((typeof(dim).parameters[1])) do y
            val = upreferred(getproperty(g, Unitful.name(y)))       # Retrieve the characteristic value from structure g
            pow = Float64(y.power)                                  # power by which it should be multiplied   
            char_val *= val^pow                                     # multiply characteristic value
        end
        param.val = upreferred.(param.val)/char_val;
    else
        param = param # The parameter has no units, so there is no way to determine how to nondimensionize it 
    end
end

"""
    Dimensionalize(param, param_dim::Unitful.FreeUnits, CharUnits::GeoUnits{TYPE})

Dimensionalizes `param` into the dimensions `param_dim` using the characteristic values specified in `CharUnits`.  

# Example
```julia-repl
julia> CharUnits =   GEO_units();
julia> v_ND      =   Nondimensionalize(3cm/yr, CharUnits) 
0.031688087814028945
julia> v_dim     =   Dimensionalize(v_ND, cm/yr, CharUnits) 
3.0 cm yr⁻¹
```

"""
function Dimensionalize(param_ND, param_dim::Unitful.FreeUnits, g::GeoUnits{TYPE}) where {TYPE}

    dim         =   Unitful.dimension(param_dim);                   # Basic SI units
    char_val    =   1.0;
    foreach((typeof(dim).parameters[1])) do y
        val = upreferred(getproperty(g, Unitful.name(y)))       # Retrieve the characteristic value from structure g
        pow = Float64(y.power)                                  # power by which it should be multiplied   
        char_val *= val^pow                                     # multiply characteristic value
    end
    param = uconvert.(param_dim, param_ND*char_val)
  
    return param
end

"""
    Dimensionalize!(param::GeoUnit, CharUnits::GeoUnits{TYPE})

Dimensionalizes `param` again to the values that it used to have using the characteristic values specified in `CharUnits`.  

# Example
```julia-repl
julia> CharUnits =   GEO_units();
julia> x = GeoUnit(3cm/yr)
julia> Nondimensionalize!(x, CharUnits)
julia> Dimensionalize!(x, CharUnits) 
3.0 cm yr⁻¹
```

"""
function Dimensionalize!(param::GeoUnit, g::GeoUnits{TYPE}) where {TYPE}
    
    dim         =   Unitful.dimension(param.unit);                   # Basic SI units
    char_val    =   1.0;
    foreach((typeof(dim).parameters[1])) do y
        val = upreferred(getproperty(g, Unitful.name(y)))       # Retrieve the characteristic value from structure g
        pow = Float64(y.power)                                  # power by which it should be multiplied   
        char_val *= val^pow                                     # multiply characteristic value
    end
    param.val = uconvert.(param.unit, param.val*char_val)
  
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


# This replaces the viewer of the Unitful package
"""
    superscript(i::Rational)
Prints exponents. 

Note that we redefine this method (from Unitful) here, as we regularly deal with exponents 
such as Pa^-4.2, which are not so nicely displayed in the Unitful package
"""
function superscript(i::Rational)
    val = Float64(i);   # the numerical value of the exponent
    v = get(ENV, "UNITFUL_FANCY_EXPONENTS", Sys.isapple() ? "true" : "false")
    t = tryparse(Bool, lowercase(v))
    k = (t === nothing) ? false : t
    if k
         return superscript_mac(Float64(i.num/i.den))
    else
        i.den == 1 ? "^" * string(val) : "^" * replace(string(val), "//" => "/")
    end
end

function superscript_mac(number::Float64)
    # this transfers the float number to a superscript string, which can be used for visualization on mac 
    x           =   trunc(Int,number)               # retrieve the part before the decimal point
    dec         =   number-x                        # decimal part
    super_str   =   superscriptnumber(x)            # transfer to superscript string
    str         =   string(number);
    ind         =   findfirst(isequal('.'), str)    
    if ~isempty(ind)  && abs(dec)>0
        st        = [super_str;'\uB7']              # dot
        for id=ind+1:length(str)
            y         = parse(Int,str[id]);
            st        = [st; superscriptnumber(y)]
        end

        super_str = join(st)                        # add decimal part
    end

    return super_str
end

function superscriptnumber(i::Int)
    if i < 0
        c = [Char(0x207B)]
    else
        c = []
    end
    for d in reverse(digits(abs(i)))
        if d == 0 push!(c, Char(0x2070)) end
        if d == 1 push!(c, Char(0x00B9)) end
        if d == 2 push!(c, Char(0x00B2)) end
        if d == 3 push!(c, Char(0x00B3)) end
        if d > 3 push!(c, Char(0x2070+d)) end
    end
    return join(c)
end


end
