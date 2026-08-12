---
---

# Nondimensionalization {#Nondimensionalization}

Create a nondimensionalization object in which we specify characteristic values, which are later employed to non-dimensionalize (or dimensionalize) all model parameters. Choose between `GEO`, `SI` or `NO` units:
- `SI` units assume all input and output is in `SI` units. Very general, but for typical geodynamic simulations often not so useful (as a million years has many seconds, resulting in large numbers).
  
- `GEO` units uses `SI` units throughout but assumes that input/output are in a format that is more convenient for typical geodynamic use cases, such as `Myr`,`km` or `MPa`
  
- `NO` units are nondimensional units. Note that for parameters to be correctly non-dimensionalized in this case, you still have to indicate units (such as that `velocity` is given in `m/s`).
  

A dimensional parameter can be transformed into a non-dimensional one with `nondimensionalize`.

# Specify characteristic values {#Specify-characteristic-values}

Characteristic values can be defined in 3 ways.
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.Units' href='#GeoParams.Units'><span class="jlbinding">GeoParams.Units</span></a> <Badge type="info" class="jlObjectType jlModule" text="Module" /></summary>



```julia
This provides units and creates a non-dimensionalization object
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Units.jl#L1-L3" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.Units.AbstractGeoUnit' href='#GeoParams.Units.AbstractGeoUnit'><span class="jlbinding">GeoParams.Units.AbstractGeoUnit</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



AbstractGeoUnit

Abstract supertype for geo units.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Units.jl#L101-L105" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.Units.GeoUnit' href='#GeoParams.Units.GeoUnit'><span class="jlbinding">GeoParams.Units.GeoUnit</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
Structure that holds a GeoUnit parameter, and encodes the units and whether it is dimensional or not in the name.

Having that is useful, as non-dimensionalization removes the units from a number and we thus no longer know how to transfer it back to the correct units.
With the GeoUnit struct, this information is retained, and we can thus redimensionalize it at a later StepRange
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Units.jl#L115-L122" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.Units.GeoUnits' href='#GeoParams.Units.GeoUnits'><span class="jlbinding">GeoParams.Units.GeoUnits</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
GeoUnits

Structure that holds parameters used for non-dimensionalization
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Units.jl#L324-L328" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.Units.GEO_units' href='#GeoParams.Units.GEO_units'><span class="jlbinding">GeoParams.Units.GEO_units</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
GEO_units(;length=1000km, temperature=1000C, stress=10MPa, viscosity=1e20Pas)
```


Creates a non-dimensionalization object using GEO units.

GEO units implies that upon dimensionalization, `time` will be in `Myr`, `length` in `km`, stress in `MPa`, etc. which is more convenient for typical geodynamic simulations than SI units The characteristic values given as input can be in arbitrary units (`km` or `m`), provided the unit is specified.

**Examples:**

```julia
julia> CharUnits = GEO_units()
Employing GEO units
Characteristic values:
         length:      1000 km
         time:        0.3169 Myr
         stress:      10 MPa
         temperature: 1000.0 °C
julia> CharUnits.velocity
1.0e-7 m s⁻¹
```


If we instead have a crustal-scale simulation, it is likely more appropriate to use a different characteristic `length`:

```julia
julia> CharUnits = GEO_units(length=10km)
Employing GEO units
Characteristic values:
         length:      10 km
         time:        0.3169 Myr
         stress:      10 MPa
         temperature: 1000.0 °C
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Units.jl#L377-L408" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.Units.SI_units' href='#GeoParams.Units.SI_units'><span class="jlbinding">GeoParams.Units.SI_units</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
CharUnits = SI_units(length=1000m, temperature=1000K, stress=10Pa, viscosity=1e20)
```


Specify the characteristic values using SI units

**Examples:**

```julia
julia> CharUnits = SI_units(length=1000m)
Employing SI units
Characteristic values:
         length:      1000 m
         time:        1.0e19 s
         stress:      10 Pa
         temperature: 1000.0 K
```


Note that the same can be achieved if the input is given in `km`:

```julia
julia> CharUnits = SI_units(length=1km)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Units.jl#L447-L466" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.Units.NO_units' href='#GeoParams.Units.NO_units'><span class="jlbinding">GeoParams.Units.NO_units</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
CharUnits = NO_units(length=1, temperature=1, stress=1, viscosity=1)
```


Specify the characteristic values in non-dimensional units

**Examples:**

```julia
julia> using GeoParams;
julia> CharUnits = NO_units()
Employing NONE units
Characteristic values:
         length:      1
         time:        1.0
         stress:      1
         temperature: 1.0
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Units.jl#L505-L521" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.Units.@unpack_val' href='#GeoParams.Units.@unpack_val'><span class="jlbinding">GeoParams.Units.@unpack_val</span></a> <Badge type="info" class="jlObjectType jlMacro" text="Macro" /></summary>



```julia
This unpacks the numerical value of `GeoUnit` in a structure, without the units.
All requested variables must be GeoUnits.

This is a modification of the `@unpack` macro as implemented in the `UnPack.jl` package, which can be used to retrieve the full variables.
```


**Example**

```julia
julia> struct Density{T} 
        ρ::GeoUnit{T} 
        α::GeoUnit{T} 
       end
julia> r = Density(GeoUnit(100kg/m^3),GeoUnit(4e-5/K));
julia> @unpack_val ρ,α = r
julia> α
4.0e-5
julia> typeof(α)
Float64
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/unpack.jl#L4-L23" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.Units.@unpack_units' href='#GeoParams.Units.@unpack_units'><span class="jlbinding">GeoParams.Units.@unpack_units</span></a> <Badge type="info" class="jlObjectType jlMacro" text="Macro" /></summary>



```julia
This unpacks the numerical value with units of `GeoUnit` parameters in a structure
All requested variables must be GeoUnits.

This is a modification of the `@unpack` macro as implemented in the `UnPack.jl` package, which can be used to retrieve the full variables.
```


**Example**

```julia
julia> struct Density{T} 
        ρ::GeoUnit{T} 
        α::GeoUnit{T} 
       end
julia> r = Density(GeoUnit(100kg/m^3),GeoUnit(4e-5/K));
julia> @unpack_units ρ,α = r
julia> α
4.0e-5 K⁻¹·⁰
julia> typeof(α)
Quantity{Float64, 𝚯⁻¹·⁰, Unitful.FreeUnits{(K⁻¹·⁰,), 𝚯⁻¹·⁰, nothing}}
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/unpack.jl#L46-L65" target="_blank" rel="noreferrer">source</a></Badge>

</details>


# (Non)-dimensionalize parameters {#Non-dimensionalize-parameters}

Once characteristic values have been defined, you can use them to non-dimensionalize or dimensionalize any parameter.
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.Units.nondimensionalize' href='#GeoParams.Units.nondimensionalize'><span class="jlbinding">GeoParams.Units.nondimensionalize</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
nondimensionalize(param, CharUnits::GeoUnits{TYPE})
```


Nondimensionalizes `param` using the characteristic values specified in `CharUnits`

**Example 1**

```julia
julia> using GeoParams;
julia> CharUnits =   GEO_units();
julia> v         =   3cm/yr
3 cm yr⁻¹
julia> v_ND      =   nondimensionalize(v, CharUnits)
0.009506426344208684
```


**Example 2**

In geodynamics one sometimes encounters more funky units

```julia
julia> CharUnits =   GEO_units();
julia> A         =   6.3e-2MPa^-3.05*s^-1
0.063 MPa⁻³·⁰⁵ s⁻¹
julia> A_ND      =   nondimensionalize(A, CharUnits)
7.068716262102384e14
```


In case you are interested to see how the units of `A` look like in different units, use this function from the [Unitful](https://github.com/PainterQubits/Unitful.jl) package:

```julia
julia> uconvert(u"Pa^-3.05*s^-1",A)
3.157479571851836e-20 Pa⁻³·⁰⁵
```


and to see it decomposed in the basic `SI` units of length, mass and time:

```julia
julia> upreferred(A)
3.1574795718518295e-20 m³·⁰⁵ s⁵·¹ kg⁻³·⁰⁵
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Units.jl#L555-L589" target="_blank" rel="noreferrer">source</a></Badge>



```julia
param = nondimensionalize(param::NTuple{N,Quantity}, g::GeoUnits{TYPE})
```


nondimensionalizes a tuple of parameters


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Units.jl#L665-L669" target="_blank" rel="noreferrer">source</a></Badge>



```julia
MatParam_ND = nondimensionalize(MatParam::AbstractMaterialParam, CharUnits::GeoUnits{TYPE})
```


Non-dimensionalizes a material parameter structure (e.g., Density, CreepLaw)

****


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Units.jl#L702-L707" target="_blank" rel="noreferrer">source</a></Badge>



```julia
nondimensionalize(phase_mat::MaterialParams, g::GeoUnits{TYPE})
```


nondimensionalizes all fields within the Material Parameters structure that contain material parameters


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Units.jl#L733-L737" target="_blank" rel="noreferrer">source</a></Badge>



```julia
MatParam_ND = nondimensionalize(MatParam::NTuple{N, AbstractMaterialParamsStruct}, CharUnits::GeoUnits)
```


Non-dimensionalizes a tuple of material parameter structures (e.g., Density, CreepLaw)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Units.jl#L773-L778" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.Units.dimensionalize' href='#GeoParams.Units.dimensionalize'><span class="jlbinding">GeoParams.Units.dimensionalize</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
dimensionalize(param, param_dim::Unitful.FreeUnits, CharUnits::GeoUnits{TYPE})
```


Dimensionalizes `param` into the dimensions `param_dim` using the characteristic values specified in `CharUnits`.

**Example**

```julia
julia> CharUnits =   GEO_units();
julia> v_ND      =   nondimensionalize(3cm/yr, CharUnits)
0.031688087814028945
julia> v_dim     =   dimensionalize(v_ND, cm/yr, CharUnits)
3.0 cm yr⁻¹
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Units.jl#L785-L799" target="_blank" rel="noreferrer">source</a></Badge>



```julia
dimensionalize(MatParam::AbstractMaterialParam, CharUnits::GeoUnits{TYPE})
```


Dimensionalizes a material parameter structure (e.g., Density, CreepLaw)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Units.jl#L826-L831" target="_blank" rel="noreferrer">source</a></Badge>



```julia
Dimensionalize(phase_mat::MaterialParams, g::GeoUnits{TYPE})
```


Dimensionalizes all fields within the Material Parameters structure that contain material parameters


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Units.jl#L860-L864" target="_blank" rel="noreferrer">source</a></Badge>



```julia
dimensionalize(MatParam::NTuple{N, AbstractMaterialParamsStruct}, CharUnits::GeoUnits)
```


dimensionalizes a tuple of material parameter structures (e.g., Density, CreepLaw)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Units.jl#L907-L911" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.Units.isDimensional' href='#GeoParams.Units.isDimensional'><span class="jlbinding">GeoParams.Units.isDimensional</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
isDimensional(MatParam::AbstractMaterialParam)
```


`true` if MatParam is in dimensional units.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Units.jl#L935-L939" target="_blank" rel="noreferrer">source</a></Badge>

</details>

