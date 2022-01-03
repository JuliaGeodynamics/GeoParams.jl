"""
    This provides units and creates a non-dimensionalization object
"""
module Units
using Unitful
import Unitful: superscript
using Parameters
using Setfield # allows modifying fields in immutable struct

import Base.show, Base.isapprox, Base.isequal, Base.convert
using GeoParams: AbstractMaterialParam, AbstractMaterialParamsStruct, AbstractPhaseDiagramsStruct, PerpleX_LaMEM_Diagram


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
const mm    = u"mm"
const Myrs  = u"Myrs"
const yr    = u"yr"
const s     = u"s"
const kg    = u"kg"
const g     = u"g"
const Pa    = u"Pa"
const MPa   = u"MPa"
const kbar  = u"kbar"
const Pas   = u"Pa*s"
const K     = u"K"
const C     = u"°C"
const mol   = u"mol"  
const kJ    = u"kJ"
const J     = u"J"
const Watt  = u"W"
const μW    = u"μW"
const μm    = u"μm"

export 
    km, m, cm, mm, μm, Myrs, yr, s, MPa, Pa, kbar, Pas, K, C, g, kg, mol, J, kJ, Watt, μW, 
    GeoUnit, GeoUnits, GEO_units, SI_units, NO_units, AbstractGeoUnits, 
    Nondimensionalize, Nondimensionalize!, Dimensionalize, Dimensionalize!,
    superscript, upreferred, GEO, SI, NONE, isDimensional, Value, NumValue, Unit, UnitValue,
    isdimensional

export AbstractGeoUnit1,   GeoUnit1  

"""
AbstractGeoUnits

Abstract supertype for geo units.
"""
abstract type AbstractGeoUnit{TYPE, DIMENSIONAL} <: Number end

abstract type AbstractUnitType end

struct GEO <: AbstractUnitType end
struct SI  <: AbstractUnitType end
struct NONE<: AbstractUnitType end


# The GeoUnit struct encodes dimensional info in the type info
"""
    Structure that holds a GeoUnit parameter, and encodes the units and whether it is dimensional or not in the name.

    Having that is useful, as non-dimensionalization removes the units from a number and we thus no longer know how to transfer it back to the correct units.
    With the GeoUnit struct, this information is retained, and we can thus redimensionalize it at a later StepRange


"""
struct GeoUnit{T}
    val::T
    unit::Unitful.FreeUnits
    isdimensional::Bool
end

# Different ways of specifying the GeoUnit:
GeoUnit(val) = GeoUnit{typeof(ustrip.(val))}(  ustrip.(val), unit(val[1]),
                           isa(val[1], Union{Unitful.FreeUnits, Unitful.Quantity} ) )

GeoUnit{T}(val) where {T,U}   = GeoUnit{T}( T.(ustrip.(val)), unit(val[1]),
                                            isa(val[1], Union{Unitful.FreeUnits, Unitful.Quantity} ) )


#Base.convert(t::Type{GeoUnit{T,U}}, x::GeoUnit)      where {T,U} =  (t(x.val,x.isdimensional), unit(x.val[1]))
#Base.convert(t::Type{GeoUnit{T,U}}, x::Quantity)     where {T,U} =  GeoUnit{T,unit(x)}(T(ustrip(x)), true)
#Base.convert(t::Type{GeoUnit{T,U}}, x::GeoUnit{T})  where {T,U}   =  x
Base.convert(t::Type{GeoUnit{T}}, x::Quantity)  where {T} = GeoUnit(T.(x))
Base.convert(t::Type{GeoUnit{T}}, x::T)  where {T} = GeoUnit(x)



# Define methods to deal with cases when the input has integers 
GeoUnit(val::Union{Int64, AbstractArray{Int64}}) = GeoUnit(Float64.(val))
GeoUnit(val::Union{Int32, AbstractArray{Int32}}) = GeoUnit(Float32.(val))
GeoUnit(val::Union{Quantity{Int64}, AbstractArray{<:Quantity{<:Int64}}}) = GeoUnit(Float64.(val))
GeoUnit(val::Union{Quantity{Int32}, AbstractArray{<:Quantity{<:Int32}}}) = GeoUnit(Float32.(val))

# helper functions
Unit(v::GeoUnit{T}) where {T} = v.unit
isdimensional(v::GeoUnit{T})  where {T}     =   v.isdimensional         # is it a nondimensional number or not?
NumValue(v::GeoUnit)                        =   v.val                   # numeric value, with no units
function UnitValue(v::GeoUnit{T}) where {T,U}
    if v.isdimensional
        return v.val*v.unit             # returns value with units
    else
        return v.val
    end
end

Base.isapprox(x::GeoUnit, y::Number; kwargs...) = Base.isapprox(x.val, y; kwargs...)
Base.isequal(x::GeoUnit,  y::Number) = Base.isequal(x.val, y)
Base.isequal(x::GeoUnit,  y::AbstractArray) = Base.isequal(x.val, y)
Base.isequal(x::GeoUnit,  y::GeoUnit) = Base.isequal(x.val, y.val)

Base.convert(::Type{<:AbstractArray},   v::GeoUnit)          =   v.val
Base.convert(::Type{<:Number},          v::GeoUnit)          =   v.val
Base.convert(::Type{GeoUnit},  v::Number)           =   GeoUnit(v) 
Base.convert(::Type{GeoUnit},  v::Int32)            =   GeoUnit(Float32(v)) 
Base.convert(::Type{GeoUnit},  v::Int64)            =   GeoUnit(Float64(v)) 
Base.convert(::Type{GeoUnit},  v::Quantity)         =   GeoUnit(v) 
Base.convert(::Type{GeoUnit},  v::AbstractArray)    =   GeoUnit(v) 
Base.convert(::Type{GeoUnit},  v::AbstractArray{Int32}) =   GeoUnit(ustrip.(Float32.(v))*unit(v[1])) 
Base.convert(::Type{GeoUnit},  v::AbstractArray{Int64}) =   GeoUnit(ustrip.(Float64.(v))*unit(v[1])) 

Base.promote_rule(::Type{GeoUnit}, ::Type{Quantity}) = GeoUnit

function Base.show(io::IO, x::GeoUnit{T})  where {T} # output
    val = x.val
    if x.isdimensional == true
        println("GeoUnit{dimensional, $(x.unit)}, ")
    else    
        println("GeoUnit{nondimensional, $(x.unit)}, ")
    end
    show(io, MIME("text/plain"), val)
end

# define a few basic routines so we can easily operate with GeoUnits
Base.length(v::GeoUnit)         =   length(v.val)
Base.size(v::GeoUnit)           =   size(v.val)

# Multiply with numnber
Base.:*(x::GeoUnit, y::Number)  =   x.val*y
Base.:+(x::GeoUnit, y::Number)  =   x.val+y
Base.:/(x::GeoUnit, y::Number)  =   x.val/y
Base.:-(x::GeoUnit, y::Number)  =   x.val-y

Base.:*(x::Number, y::GeoUnit)  = y.val*x
Base.:+(x::Number, y::GeoUnit)  = y.val+x
Base.:/(x::Number, y::GeoUnit)  = x/y.val
Base.:-(x::Number, y::GeoUnit)  = x-y.val

# Multiplying a GeoUnit with another one, returns a GeoUnit
Base.:+(x::GeoUnit{T1}, y::GeoUnit{T2}) where {T1,T2} = GeoUnit(x.val*x.unit + y.val*y.unit) 
Base.:/(x::GeoUnit{T1}, y::GeoUnit{T2}) where {T1,T2} = GeoUnit( (x.val*x.unit) / (y.val*y.unit) ) 
Base.:*(x::GeoUnit{T1}, y::GeoUnit{T2}) where {T1,T2} = GeoUnit(x.val*x.unit * y.val*y.unit)  
Base.:-(x::GeoUnit{T1}, y::GeoUnit{T2}) where {T1,T2} = GeoUnit(x.val*x.unit - y.val*y.unit) 
Base.:^(x::GeoUnit{T1}, y::GeoUnit{T2}) where {T1,T2} = GeoUnit(x.val*x.unit^(y.val*y.unit)) 

Base.:*(x::GeoUnit, y::Quantity)   = UnitValue(x).*y
Base.:/(x::GeoUnit, y::Quantity)   = UnitValue(x)./y
Base.:+(x::GeoUnit, y::Quantity)   = UnitValue(x) .+ y
Base.:-(x::GeoUnit, y::Quantity)   = UnitValue(x) .- y

Base.:*(x::Quantity, y::GeoUnit)   = x.*UnitValue(y)
Base.:/(x::Quantity, y::GeoUnit)   = x./UnitValue(y)
Base.:+(x::Quantity, y::GeoUnit)   = x .+ UnitValue(y)
Base.:-(x::Quantity, y::GeoUnit)   = x .- UnitValue(y) 

# If we multiply a GeoUnit with an abstract array, we only return values, not units (use GeoUnits for that)
Base.:*(x::GeoUnit, y::AbstractArray)   = NumValue(x).*y 
Base.:/(x::GeoUnit, y::AbstractArray)   = NumValue(x) ./y 
Base.:+(x::GeoUnit, y::AbstractArray)   = NumValue(x) .+ y 
Base.:-(x::GeoUnit, y::AbstractArray)   = NumValue(x) .-y

Base.:*(y::AbstractArray, x::GeoUnit)   = NumValue(x).*y
Base.:/(y::AbstractArray, x::GeoUnit)   = y./NumValue(x)
Base.:+(y::AbstractArray, x::GeoUnit)   = NumValue(x) .+ y
Base.:-(y::AbstractArray, x::GeoUnit)   = y .- NumValue(x)

Base.getindex(x::GeoUnit{T}, i::Int64, j::Int64, k::Int64) where {T,U} = GeoUnit(x.val[i,j,k]*x.unit)
Base.getindex(x::GeoUnit{T}, i::Int64, j::Int64) where {T,U} = GeoUnit(x.val[i,j]*x.unit)
Base.getindex(x::GeoUnit{T}, i::Int64) where {T,U} = GeoUnit(x.val[i]*x.unit)

Base.setindex!(x::GeoUnit{T}, v::Any, i::Int64, j::Int64, k::Int64)  where {T} = x.val[i,j,k] = v
Base.setindex!(x::GeoUnit{T}, v::Any, i::Int64, j::Int64)  where {T} = x.val[i,j] = v
Base.setindex!(x::GeoUnit{T}, v::Any, i::Int64) where {T}  = x.val[i] = v

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
    heatcapacity    =   J/kg/K 
    conductivity    =   W/m/K
    
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
Employing GEO units 
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
Employing GEO units 
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
Employing SI units 
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
Employing NONE units 
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
function Nondimensionalize(param::GeoUnit{T}, g::GeoUnits{TYPE}) where {T,TYPE}
    if param.isdimensional
        dim         =   Unitful.dimension(param.unit);                   # Basic SI units
        char_val    =   1.0;
        foreach((typeof(dim).parameters[1])) do y
            val = upreferred(getproperty(g, Unitful.name(y)))       # Retrieve the characteristic value from structure g
            pow = Float64(y.power)                                  # power by which it should be multiplied   
            char_val *= val^pow                                     # multiply characteristic value
        end
        val_ND = upreferred.(param.val*param.unit)/char_val
        param_ND = GeoUnit{T}(val_ND, param.unit, false)                      # store new value, but keep original dimensions
    else
        param_ND = param;
    end
    return param_ND
end

# in case the parameter is already non-dimensional:
Nondimensionalize(param::String, g::GeoUnits{TYPE}) where {TYPE} = param 

# in case it is a uniful quantity
function Nondimensionalize(param::Union{Unitful.Quantity, AbstractArray{<:Quantity}}, g::GeoUnits{TYPE}) where {TYPE} 
    
    result = Nondimensionalize(GeoUnit(param),g)

    return UnitValue(result)

end

# If it is an array, but has no units we cannot know how to nondimensionalize it
Nondimensionalize(param::AbstractArray{<:Number}, g::GeoUnits{TYPE}) where {TYPE} = param

#=
function Nondimensionalize(param, g::GeoUnits{TYPE}) where {TYPE}
    if typeof(param) == String
        param_ND = param # The parameter is a string, cannot be nondimensionalized
    elseif unit(param)!=NoUnits 
        dim         =   Unitful.dimension(param);                   # Basic SI units
        char_val    =   1.0;
        foreach((typeof(dim).parameters[1])) do y
            val = upreferred(getproperty(g, Unitful.name(y)))       # Retrieve the characteristic value from structure g
            pow = Float64(y.power)                                  # power by which it should be multiplied   
            char_val *= val^pow                                     # multiply characteristic value
        end
        param_ND = upreferred.(param.val*unit(param))/char_val
    else
        param_ND = param # The parameter has no units, so there is no way to determine how to nondimensionize it 
    end
    return param_ND
end

function Nondimensionalize(param::GeoUnit, g::GeoUnits{TYPE}) where {TYPE}
    if param.dimensional
        dim         =   Unitful.dimension(param.unit);                   # Basic SI units
        char_val    =   1.0;
        foreach((typeof(dim).parameters[1])) do y
            val = upreferred(getproperty(g, Unitful.name(y)))       # Retrieve the characteristic value from structure g
            pow = Float64(y.power)                                  # power by which it should be multiplied   
            char_val *= val^pow                                     # multiply characteristic value
        end
        param_ND = upreferred.(param.val*param.unit)/char_val
    else
        param_ND = param
    end
    return param_ND

end

function Nondimensionalize(param::Array, g::GeoUnits{TYPE}) where {TYPE}
    if  param[1] !=NoUnits
        dim         =   Unitful.dimension.(param);                   # Basic SI units
        char_val    =   1.0;
        foreach((typeof(dim[1]).parameters[1])) do y
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
=#

## Debugging; testing alternative way to define a GeoUnit and nondimensionalize it


## 


#=
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

function Nondimensionalize!(param::GeoUnit{T,U,true}, g::GeoUnits{TYPE}) where {T,U,TYPE}
        dim         =   Unitful.dimension(U);                   # Basic SI units
        char_val    =   1.0;
        foreach((typeof(dim).parameters[1])) do y
            val = upreferred(getproperty(g, Unitful.name(y)))       # Retrieve the characteristic value from structure g
            pow = Float64(y.power)                                  # power by which it should be multiplied   
            char_val *= val^pow                                     # multiply characteristic value
        end
        param.val = upreferred.(param.val*U)/char_val;
end


function Nondimensionalize!(param::String, g::GeoUnits{TYPE}) where {TYPE}
    param_ND = param
    return nothing
end

function Nondimensionalize(param::String, g::GeoUnits{TYPE}) where {TYPE}
    return param
end
=#

#=
# TO BE FIXED OR REMOVED!
function Nondimensionalize!(param::Array, g::GeoUnits{TYPE}) where {TYPE}
  
    if (param.unit != NoUnits) & isdimensional(param)
        dim         =   Unitful.dimension(param.unit);                   # Basic SI units
        char_val    =   1.0;
        foreach((typeof(dim).parameters[1])) do y
            val = upreferred(getproperty(g, Unitful.name(y)))       # Retrieve the characteristic value from structure g
            pow = Float64(y.power)                                  # power by which it should be multiplied   
            char_val *= val^pow                                     # multiply characteristic value
        end
        param = upreferred.(param.val.*param.unit)/char_val;
    else
        param = param # The parameter has no units, so there is no way to determine how to nondimensionize it 
    end

end
=#

"""
    MatParam_ND = Nondimensionalize(MatParam::AbstractMaterialParam, CharUnits::GeoUnits{TYPE})

Non-dimensionalizes a material parameter structure (e.g., Density, CreepLaw)

"""
function Nondimensionalize(MatParam::AbstractMaterialParam, g::GeoUnits{TYPE}) where {TYPE} 
    for param in fieldnames(typeof(MatParam))
        if isa(getfield(MatParam, param), GeoUnit)
            z=getfield(MatParam, param)
            
            # non-dimensionalize:
            z = Nondimensionalize(z, g)  
            
            # Replace field (using Setfield package):
            MatParam = set(MatParam, Setfield.PropertyLens{param}(), z)
        
        end
    end
    return MatParam 
end
    
"""
    Nondimensionalize!(phase_mat::MaterialParams, g::GeoUnits{TYPE})

Nondimensionalizes all fields within the Material Parameters structure that contain material parameters
"""
function Nondimensionalize!(phase_mat::AbstractMaterialParamsStruct, g::GeoUnits{TYPE}) where {TYPE} 

    for param in fieldnames(typeof(phase_mat))
        fld = getfield(phase_mat, param)
        if length(fld)>0
            if typeof(fld[1]) <: AbstractPhaseDiagramsStruct
                
                # in case we employ a phase diagram 
                temp = PerpleX_LaMEM_Diagram(fld[1].Name, CharDim = g)
                fld_new = (temp,)
                setfield!(phase_mat, param, fld_new)
            else
                # otherwise non-dimensionalize all fields and create a new tuple
                id = findall(isa.(fld,AbstractMaterialParam))
                if length(id)>0
                    # Create a new tuple with non-dimensionalized fields:
                    fld_new = ntuple(i -> Units.Nondimensionalize(fld[id[i]],g), length(id) )
                    setfield!(phase_mat, param, fld_new)        # to be changed for immutable struct
                end
            end
        end
    end
    phase_mat.Nondimensional = true
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
function Dimensionalize(param_ND::GeoUnit{T}, g::GeoUnits{TYPE}) where {T,TYPE}
    if isdimensional(param_ND)==false
        dim         =   Unitful.dimension(param_ND.unit);                   # Basic SI units
        char_val    =   1.0;
        foreach((typeof(dim).parameters[1])) do y
            val = upreferred(getproperty(g, Unitful.name(y)))       # Retrieve the characteristic value from structure g
            pow = Float64(y.power)                                  # power by which it should be multiplied   
            char_val *= val^pow                                     # multiply characteristic value
        end
        val = uconvert.(param_ND.unit, param_ND.val*char_val)
        param = GeoUnit{T}(ustrip.(val), param_ND.unit, true)     # store new value, but keep original dimensions

        return param
    else
        return param_ND
    end
end

function Dimensionalize(param_ND, param_dim::Unitful.FreeUnits, g::GeoUnits{TYPE}) where {TYPE}

    dim         =   Unitful.dimension(param_dim);               # Basic SI units
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
    Dimensionalize!(MatParam::AbstractMaterialParam, CharUnits::GeoUnits{TYPE})

Dimensionalizes a material parameter structure (e.g., Density, CreepLaw)

"""
function Dimensionalize(MatParam::AbstractMaterialParam, g::GeoUnits{TYPE}) where {TYPE} 

    for param in fieldnames(typeof(MatParam))
        if isa(getfield(MatParam, param), GeoUnit)
            z=getfield(MatParam, param)
            z_dim = Dimensionalize(z, g)

            # Replace field (using Setfield package):
            MatParam = set(MatParam, Setfield.PropertyLens{param}(), z_dim)

        end
    end
    return MatParam
end


"""
    Dimensionalize!(phase_mat::MaterialParams, g::GeoUnits{TYPE})

Dimensionalizes all fields within the Material Parameters structure that contain material parameters
"""
function Dimensionalize!(phase_mat::AbstractMaterialParamsStruct, g::GeoUnits{TYPE}) where {TYPE} 

    for param in fieldnames(typeof(phase_mat))
        fld = getfield(phase_mat, param)
        if length(fld)>0
            id = findall(isa.(fld,AbstractMaterialParam))
            if length(id)>0
                # Create a new tuple with non-dimensionalized fields:
                fld_new = ntuple(i -> Units.Dimensionalize(fld[id[i]],g), length(id) )
              

                setfield!(phase_mat, param, fld_new)        # to be changed for immutable struct
               # phase_mat = set(phase_mat, Setfield.PropertyLens{param}(), fld_new)
            end
        end
    end
    phase_mat.Nondimensional = false
end


## Debugging

# This dimensionalizes the GeoUnit1 parameter


## 

#=
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
    
    if !param.dimensional
        dim         =   Unitful.dimension(param.unit);                   # Basic SI units
        char_val    =   1.0;
        foreach((typeof(dim).parameters[1])) do y
            val = upreferred(getproperty(g, Unitful.name(y)))       # Retrieve the characteristic value from structure g
            pow = Float64(y.power)                                  # power by which it should be multiplied   
            char_val *= val^pow                                     # multiply characteristic value
        end
        param.val = ustrip.(uconvert.(param.unit, param.val*char_val))
        param.dimensional = true
    end
  
end





=#

"""
    isDimensional(MatParam::AbstractMaterialParam)

`true` if MatParam is in dimensional units.    
"""
function isDimensional(MatParam::AbstractMaterialParam)
    isDim = false;
    for param in fieldnames(typeof(MatParam))
        if isa(getfield(MatParam, param), GeoUnit)
            z=getfield(MatParam, param)
            if z.isdimensional==true
                isDim=true;
            end
        end
    end
    return isDim
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




# This replaces the viewer of the Unitful package, such that the printing of units is done as floats (better)
function Unitful.superscript(i::Rational{Int64}; io=nothing) 
    string(superscript(float(i)))
    if io === nothing
        iocontext_value = nothing
    else
        iocontext_value = get(io, :fancy_exponent, nothing)
    end
    if iocontext_value isa Bool
        fancy_exponent = iocontext_value
    else
        v = get(ENV, "UNITFUL_FANCY_EXPONENTS", Sys.isapple() ? "true" : "false")
        t = tryparse(Bool, lowercase(v))
        fancy_exponent = (t === nothing) ? false : t
    end
    if fancy_exponent
        return superscript(float(i)) 
    else
        return  "^" * string(float(i)) 
    end

end

Unitful.superscript(i::Float64) = map(repr(i)) do c
    c == '-' ? '\u207b' :
    c == '1' ? '\u00b9' :
    c == '2' ? '\u00b2' :
    c == '3' ? '\u00b3' :
    c == '4' ? '\u2074' :
    c == '5' ? '\u2075' :
    c == '6' ? '\u2076' :
    c == '7' ? '\u2077' :
    c == '8' ? '\u2078' :
    c == '9' ? '\u2079' :
    c == '0' ? '\u2070' :
    c == '.' ? '\u0387' :
    error("unexpected character")
end



end
