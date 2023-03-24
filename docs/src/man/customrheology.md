# `User-defined rheology`

`CustomRheology` allows the user to interface with `GeoParams.jl` API for rheology calculations: 

```julia
struct CustomRheology{F1, F2, T} <: AbstractCreepLaw
    strain::F1
    stress::F2
    args::T
end
```

where `strain` and `stress` are functions that compute the strain rate and deviatoric stress tensors, respectively, and `args` is a `NamedTuple` containing any constant (i.e. **immutable**) physical parameters needed by our `CustomRheology`. 

## Example: Arrhenious-type rheology
The deviatoric stress $\tau_{ij}$ and strain rate tensors $\dot{\varepsilon}_{ij}$ are simply
    $$\tau_{ij} = 2 \eta \dot{\varepsilon}_{ij} \qquad \dot{\varepsilon}_{ij}  = {\tau_{ij}  \over 2 \eta}$$

and the viscosity $\eta$ is temperature-dependant
    $$\eta = \eta_{0} * \exp\left(\frac{E}{T+T_{o}}-\frac{E}{T_{\eta}+T_{o}}\right)$$

where $\eta_0$ and $T_{\eta}$ are the respective reference viscosity and temperature, $T_o$ is the offset temperature, $T$ is the local temperature, and $E$ is activation energy. 

Before defining the functions to compute $\tau$ and $\dot{\varepsilon}$, it is convinient to define a helper function to compute the viscosity:
```julia
@inline function custom_viscosity(
    a::CustomRheology; T = 0.0, kwargs...
)
    (; η0, E, To, Tη) = a.args
    η = η0 * exp(E / (T + To) - E / (Tη + To))
    return  η
end
```
Then we just need to define two simple functions to compue the second invariants of $\tau_{ij}$ and $\dot{\varepsilon}_{ij}$:

```julia
# function to compute deviatoric stress
@inline function custom_τII(
    a::CustomRheology, EpsII; kwargs...
)
    η = custom_viscosity(a; kwargs...)
    return 2.0 * (η * EpsII)
end

# function to compute strain rate
@inline function custom_εII(
    a::CustomRheology, TauII; kwargs...
)
    η = custom_viscosity(a; kwargs...)
    return (TauII / η) * 0.5
end
```
Note that the key word argument `kwargs...` is needed in all the above functions for compatibility with `GeoParams.jl`. 

Finally, we need a `NamedTuple` containing the physical parameters needed by `custom_τII` and `custom_εII`:
```julia
# constant parameters, these are typically wrapped into a struct (compulsory)
parameters = (;
    η0 = 1.0,
    E = 23.03,
    To = 1.0,
    Tη = 1.0,
)
```

Then the `CustomRheology` object is created
```julia
a = CustomRheology(custom_εII, custom_τII, parameters)
```

Now we can use our new rheology object with `GeoParams.jl` methods:

```julia-repl
julia> args = (; T=0.5)
(T = 0.5,)

julia> τII, εII = 1.5,  0.5
(1.5, 0.5)

julia> ε = compute_εII(v, τII, args)
0.016147090412110487

julia> τ = compute_τII(v, εII, args)
46.44799656521971

julia> dεdτ = dεII_dτII(v, τII, args)
0.010764726941406991

julia> dτdε = dτII_dεII(v, εII, args)
92.89599313043942
```
 
And also works for composite rheologies:

```julia-repl
julia> el = ConstantElasticity(; G=1.0)
Linear elasticity with shear modulus: G = 1.0, Poisson's ratio: ν = 0.5, bulk modulus: Kb = Inf and Young's module: E=NaN

julia> c = CompositeRheology(v, el)
--?????----/\/\/--

julia> args = (; dt=Inf, T=0.5)
(dt = Inf, T = 0.5)

julia> ε = compute_εII(c, τII, args)
0.016147090412110487

julia> τ = compute_τII(c, εII, args)
46.44799656521971
```
