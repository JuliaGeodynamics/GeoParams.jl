export MeltWeakening

"""
    MeltWeakening(law; α)

Creep law `law` weakened by a small melt fraction ``\\phi``,
```math
    \\dot{\\varepsilon}_{II}(\\tau_{II}, \\phi) = e^{\\alpha \\phi} \\, \\dot{\\varepsilon}_{II}^{\\mathrm{law}}(\\tau_{II}),
```
so the viscosity at a given stress drops by ``e^{-\\alpha\\phi}``. For water-saturated
olivine–basalt aggregates, Mei et al. (2002, Table 2) give ``\\alpha = 26 \\pm 2`` in diffusion
creep and ``\\alpha = 31 \\pm 4`` in dislocation creep, for melt fractions of 0.02–0.12. The same
values describe anhydrous aggregates below a melt fraction of about 0.05: the dislocation-creep data
of Hirth & Kohlstedt (1995b, Table 1, 1523 K) give ``\\alpha = 32 \\pm 3`` (``\\alpha = 38 \\pm 3`` with
``n = 3.2``), and their diffusion creep is 3× faster at ``\\phi = 0.04`` (1995a), as ``\\alpha = 26`` predicts.
Above ``\\phi \\approx 0.05`` anhydrous diffusion creep speeds up more: 25× at ``\\phi = 0.07`` (1995a),
where ``\\alpha = 26`` gives 6×. The form
describes partially molten rock below the disaggregation threshold; it does not describe
crystal-bearing magma.

`α` is dimensionless and must be non-negative. The melt fraction is passed at runtime,
e.g. `compute_εII(law, τII; T, ϕ = 0.05)`; without it ``\\phi = 0`` and the result equals that of `law`.
All other keyword arguments are passed on to `law`.

# Example
```julia-repl
julia> law0 = SetDiffusionCreep(GeoParams.Diffusion.wet_olivine_Mei_2000a);

julia> law = MeltWeakening(law0; α = 26);

julia> args = (; T = 1473.0, P = 300e6, d = 15e-6, f = 300e6);

julia> compute_εII(law, 1e6; args..., ϕ = 0.05) / compute_εII(law0, 1e6; args...)
3.6692966676192444
```

# References
- Hirth, G., & Kohlstedt, D. L. (1995a). Experimental constraints on the dynamics of the partially molten upper mantle: Deformation in the diffusion creep regime. Journal of Geophysical Research, 100(B2), 1981–2001. https://doi.org/10.1029/94JB02128
- Hirth, G., & Kohlstedt, D. L. (1995b). Experimental constraints on the dynamics of the partially molten upper mantle: 2. Deformation in the dislocation creep regime. Journal of Geophysical Research, 100(B8), 15441–15449. https://doi.org/10.1029/95JB01292
- Mei, S., Bai, W., Hiraga, T., & Kohlstedt, D. L. (2002). Influence of melt on the creep behavior of olivine–basalt aggregates under hydrous conditions. Earth and Planetary Science Letters, 201, 491–507. https://doi.org/10.1016/S0012-821X(02)00745-8
"""
struct MeltWeakening{T, L <: AbstractCreepLaw{T}} <: AbstractCreepLaw{T}
    law::L                          # creep law of the melt-free rock
    α::GeoUnit{T, typeof(NoUnits)}  # melt-weakening exponent

    function MeltWeakening(law::AbstractCreepLaw{T}, α) where {T}
        g = convert(GeoUnit, α)
        g.unit == NoUnits || throw(ArgumentError("MeltWeakening: α must be dimensionless, got $(Value(g))"))
        g.val >= 0 || throw(ArgumentError("MeltWeakening: α must be non-negative (ε̇ ∝ exp(αϕ)), got $(g.val)"))
        return new{T, typeof(law)}(law, GeoUnit(convert(T, g.val), NoUnits, g.isdimensional))
    end
end
MeltWeakening(law::AbstractCreepLaw; α) = MeltWeakening(law, α)

isDimensional(a::MeltWeakening) = isDimensional(a.law)

function param_info(a::MeltWeakening)
    inner = param_info(a.law)
    return MaterialParamsInfo(;
        Equation = L"\dot{\varepsilon}_{II} = e^{\alpha\phi}\,\dot{\varepsilon}_{II}^{\mathrm{law}}",
        Comment = inner.Comment,
        BibTex_Reference = inner.BibTex_Reference,
    )
end

# e^{αϕ} in the precision of the computation
@inline function _melt_weakening(a::MeltWeakening, ::Type{Tc}, ϕ) where {Tc}
    @unpack_val Tc α = a
    return exp(α * convert_precision(Tc, ϕ))
end

@inline function compute_εII(a::MeltWeakening, TauII; ϕ = 0, kwargs...)
    return _melt_weakening(a, precision_of(TauII), ϕ) * compute_εII(a.law, TauII; ϕ, kwargs...)
end

@inline function compute_τII(a::MeltWeakening, EpsII; ϕ = 0, kwargs...)
    return compute_τII(a.law, EpsII / _melt_weakening(a, precision_of(EpsII), ϕ); ϕ, kwargs...)
end

@inline function dεII_dτII(a::MeltWeakening, TauII; ϕ = 0, kwargs...)
    return _melt_weakening(a, precision_of(TauII), ϕ) * dεII_dτII(a.law, TauII; ϕ, kwargs...)
end

@inline function dτII_dεII(a::MeltWeakening, EpsII; ϕ = 0, kwargs...)
    f = _melt_weakening(a, precision_of(EpsII), ϕ)
    return dτII_dεII(a.law, EpsII / f; ϕ, kwargs...) / f
end

function compute_εII!(
        EpsII::AbstractArray{_T, N}, a::MeltWeakening, TauII::AbstractArray; ϕ = 0, kwargs...
    ) where {N, _T}
    for i in each_argument_index(EpsII, TauII, ϕ, values(kwargs)...)
        args_i = map(x -> argument_at(x, i), values(kwargs))
        EpsII[i] = compute_εII(a, convert_precision(_T, TauII[i]); ϕ = argument_at(ϕ, i), args_i...)
    end
    return nothing
end

function compute_τII!(
        TauII::AbstractArray{_T, N}, a::MeltWeakening, EpsII::AbstractArray; ϕ = 0, kwargs...
    ) where {N, _T}
    for i in each_argument_index(TauII, EpsII, ϕ, values(kwargs)...)
        args_i = map(x -> argument_at(x, i), values(kwargs))
        TauII[i] = compute_τII(a, convert_precision(_T, EpsII[i]); ϕ = argument_at(ϕ, i), args_i...)
    end
    return nothing
end

function show(io::IO, g::MeltWeakening)
    return print(io, "MeltWeakening: α = $(Value(g.α)), law = $(g.law)")
end
