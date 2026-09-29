---
---

# Phase Diagrams {#Phase-Diagrams}

It is possible to employ phase diagrams as lookup tables, which can be used to compute density, melt fraction or seismic properties, for example. 

# Perple_X {#Perple_X}

A popular way to compute phase diagrams (pseudosections) for a given (fixed) chemical composition as a function of pressure and temperature is by using [Perple_X](https://www.perplex.ethz.ch). Once a stable assemblage is computed, through Gibbs energy minimisation, it allows you to output a large number of properties such as density, melt fraction or seismic velocities.

Within GeoParams, we can import this diagram (provided it is formatted in the same manner as `LaMEM` expects `Perple_X` input to be).
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.PerpleX_LaMEM_Diagram' href='#GeoParams.PerpleX_LaMEM_Diagram'><span class="jlbinding">GeoParams.PerpleX_LaMEM_Diagram</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
PD_Data = PerpleX_LaMEM_Diagram(fname::String; CharDim = nothing, type::AbstractString = "Perple_X/LaMEM")
```


Reads a precomputed phase diagram in the `Perple_X/MAGEMin/LaMEM` format (which is a phase diagram computed using `Perple_X/MAGEMin`, but formatted in a manner that is readable using LaMEM or potentially any other geodynamic code). The data is stored in the `PhaseDiagram_LookupTable` structure.

If the `CharDim` object is specified, the values of all diagrams will be non-dimensionalized. The type kwarg is a string that indicates the type of phase diagram (default is "Perple_X/LaMEM").

**Example**

```julia
julia> PD_Data = PerpleX_LaMEM_Diagram("./test_data/Peridotite.in")
Perple_X/LaMEM Phase Diagram Lookup Table:
                      File    :   ./test_data/Peridotite.in
                      T       :   293.0 - 1573.000039
                      P       :   1.0e7 - 2.9999999944e9
                      fields  :   :meltRho, :meltRho, :meltFrac, :rockRho, :Rho, :rockVp
                                  :rockVs, :rockVpVs, :meltVp, :meltVs, :meltVpVs
                                  :Vp, :Vs, :VpVs
```


Once imported, the properties on the diagram can be interpolated in a simple manner:

```
julia> PD_Data.Rho(1500,1e7)
3042.836820256982
```


This also works for vectors or arrays:

```julia
julia> T = [1500 1800; 1233 1300]
julia> P = [1e8 1e9; 1e7 1e7]
julia> rho = PD_Data.Rho.(T,P)
```


(Note the dot `.` in front of the bracket while evaluating arrays).

The fields that are available depend on what is listed in the diagram file. The units of the fields are automatically evaluated, and employed to non-dimensionalize the parameters if `CharDim` is specified.

**Algorithm**

Internally, we employ a linear 2D interpolation scheme for evaluating the phase diagram values at arbitrary (T,P) points using bilinear interpolation. Values outside the range of the diagram are set to the boundary of the diagram. The interpolation object is directly encoded in the `PhaseDiagram_LookupTable`` object.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/PhaseDiagrams/PhaseDiagrams.jl#L75-L117" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.PhaseDiagrams.PhaseDiagram_LookupTable' href='#GeoParams.MaterialParameters.PhaseDiagrams.PhaseDiagram_LookupTable'><span class="jlbinding">GeoParams.MaterialParameters.PhaseDiagrams.PhaseDiagram_LookupTable</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
Contains data of a Phase Diagram that is regularly spaced in P & T
```


**Fields**
- `Type::Ptr{UInt8}` : String pointer indicating the type of phase diagram (e.g., "Perple_X/MAGEMin/LaMEM")
  
- `Name::Ptr{UInt8}` : String pointer to the name of the phase diagram file
  
- `rockRho::Union{T, Nothing}` : Interpolation object for rock density
  
- `meltRho::Union{T, Nothing}` : Interpolation object for melt density
  
- `meltFrac::Union{T, Nothing}` : Interpolation object for melt fraction
  
- `Rho::Union{T, Nothing}` : Interpolation object for total density
  
- `rockVp::Union{T, Nothing}` : Interpolation object for rock P-wave velocity
  
- `rockVs::Union{T, Nothing}` : Interpolation object for rock S-wave velocity
  
- `rockVpVs::Union{T, Nothing}` : Interpolation object for rock Vp/Vs ratio
  
- `meltVp::Union{T, Nothing}` : Interpolation object for melt P-wave velocity
  
- `meltVs::Union{T, Nothing}` : Interpolation object for melt S-wave velocity
  
- `meltVpVs::Union{T, Nothing}` : Interpolation object for melt Vp/Vs ratio
  
- `Vp::Union{T, Nothing}` : Interpolation object for total P-wave velocity
  
- `Vs::Union{T, Nothing}` : Interpolation object for total S-wave velocity
  
- `VpVs::Union{T, Nothing}` : Interpolation object for total Vp/Vs ratio
  
- `SpecificCp::Union{T, Nothing}` : Interpolation object for specific heat capacity
  
- `solid_Vp::Union{T, Nothing}` : Interpolation object for solid P-wave velocity
  
- `solid_Vs::Union{T, Nothing}` : Interpolation object for solid S-wave velocity
  
- `melt_bulkModulus::Union{T, Nothing}` : Interpolation object for melt bulk modulus
  
- `solid_bulkModulus::Union{T, Nothing}` : Interpolation object for solid bulk modulus
  
- `solid_shearModulus::Union{T, Nothing}` : Interpolation object for solid shear modulus
  
- `Vp_uncorrected::Union{T, Nothing}` : Interpolation object for uncorrected P-wave velocity
  
- `Vs_uncorrected::Union{T, Nothing}` : Interpolation object for uncorrected S-wave velocity
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/PhaseDiagrams/PhaseDiagrams.jl#L17-L44" target="_blank" rel="noreferrer">source</a></Badge>

</details>

