export Suspension, Roscoe1952_Crystals, Mueller2010_Crystals, IshiiZuber1979_Bubbles, Mader2013_Bubbles, Vona2011_Crystals, compute_relative_viscosity

"""
    AbstractRelativeViscosity

Relative viscosity ``\\eta_r(\\phi) = \\eta / \\eta_\\mathrm{melt}`` of a suspension with volume
fraction ``\\phi`` of crystals or bubbles; used by [`Suspension`](@ref).
"""
abstract type AbstractRelativeViscosity <: AbstractMaterialParam end

"""
    Roscoe1952_Crystals(; ϕ_max = 1/1.35, n = 2.5)

Einstein–Roscoe relative viscosity of a suspension of rigid spheres,
```math
    \\eta_r = \\left(1 - \\frac{\\phi_x}{\\phi_\\mathrm{max}}\\right)^{-n},
```
with ``\\phi_x`` the crystal volume fraction. Roscoe (1952) uses ``n = 5/2`` and gives
``\\phi_\\mathrm{max} = 1/1.35`` for equal spheres at high concentration and ``\\phi_\\mathrm{max} = 1``
for spheres of very different sizes. The Krieger & Dougherty (1959) form has ``n = B\\,\\phi_\\mathrm{max}``
with Einstein coefficient ``B``; Keller & Suckale (2019) fit ``\\phi_\\mathrm{max} = 0.58`` and ``B = 4``
(``n = 2.32``) to crystal-bearing basalt. ``\\eta_r`` diverges at ``\\phi_x = \\phi_\\mathrm{max}`` and is `NaN` outside
``0 \\le \\phi_x \\le \\phi_\\mathrm{max}``.

# References
- Roscoe, R. (1952). The viscosity of suspensions of rigid spheres. British Journal of Applied Physics, 3(8), 267–269. https://doi.org/10.1088/0508-3443/3/8/306
- Krieger, I. M., & Dougherty, T. J. (1959). A mechanism for non-Newtonian flow in suspensions of rigid spheres. Transactions of the Society of Rheology, 3, 137–152. https://doi.org/10.1122/1.548848
- Keller, T., & Suckale, J. (2019). A continuum model of multi-phase reactive transport in igneous systems. Geophysical Journal International, 219(1), 185–222. https://doi.org/10.1093/gji/ggz287
"""
struct Roscoe1952_Crystals{T} <: AbstractRelativeViscosity
    ϕ_max::T  # maximum packing fraction
    n::T      # exponent

    function Roscoe1952_Crystals{T}(ϕ_max, n) where {T}
        0 < ϕ_max <= 1 || throw(ArgumentError("Roscoe1952_Crystals: ϕ_max must lie in (0, 1], got $ϕ_max"))
        n > 0 || throw(ArgumentError("Roscoe1952_Crystals: n must be positive, got $n"))
        return new{T}(ϕ_max, n)
    end
end
Roscoe1952_Crystals(ϕ_max, n) = Roscoe1952_Crystals{float(promote_type(typeof(ϕ_max), typeof(n)))}(ϕ_max, n)
Roscoe1952_Crystals(; ϕ_max = 1 / 1.35, n = 2.5) = Roscoe1952_Crystals(ϕ_max, n)

@inline function compute_relative_viscosity(f::Roscoe1952_Crystals, ϕ)
    ϕ = float(ϕ)
    ϕ_max = convert(typeof(ϕ), f.ϕ_max)
    0 <= ϕ <= ϕ_max || return oftype(ϕ, NaN)
    return (1 - ϕ / ϕ_max)^-convert(typeof(ϕ), f.n)
end

"""
    Mueller2010_Crystals(; aspect_ratio)

Maron–Pierce relative viscosity of a suspension of particles with aspect ratio ``r_p`` = `aspect_ratio`,
after Mueller et al. (2010, eqs 2.4 and 5.1),
```math
    \\eta_r = \\left(1 - \\frac{\\phi_x}{\\phi_m}\\right)^{-2}, \\qquad \\phi_m = \\frac{2}{0.321\\,r_p + 3.02},
```
returned as the equivalent [`Roscoe1952_Crystals`](@ref) law with ``\\phi_\\mathrm{max} = \\phi_m`` and ``n = 2``.
The fit of ``\\phi_m`` covers spheres and prolate particles, ``1 \\le r_p \\le 9.17``. Mueller et al. (2010)
describe their suspensions as Herschel–Bulkley fluids; this law is their consistency, which equals ``\\eta_r`` at a shear rate of
1 s⁻¹ only for suspensions without yield stress, and omits their shear thinning (eq. 5.2) and yield stress
(eq. 4.1).

# References
- Mueller, S., Llewellin, E. W., & Mader, H. M. (2010). The rheology of suspensions of solid particles. Proceedings of the Royal Society A, 466(2116), 1201–1228. https://doi.org/10.1098/rspa.2009.0445
"""
function Mueller2010_Crystals(; aspect_ratio)
    aspect_ratio >= 1 || throw(ArgumentError("Mueller2010_Crystals: aspect_ratio must be at least 1, got $aspect_ratio"))
    return Roscoe1952_Crystals(2 / (0.321 * aspect_ratio + 3.02), 2)
end

"""
    IshiiZuber1979_Bubbles()

Ishii–Zuber relative viscosity of a suspension of undeformable spherical bubbles (capillary
number ``\\mathrm{Ca} \\ll 1``),
```math
    \\eta_r = \\frac{1}{1 - \\phi_\\mathrm{gas}},
```
with ``\\phi_\\mathrm{gas}`` the gas volume fraction. This is the Ishii & Zuber (1979) mixture
viscosity for an inviscid dispersed phase with maximum packing 1. ``\\eta_r`` is `NaN` outside
``0 \\le \\phi_\\mathrm{gas} \\le 1``.

# References
- Ishii, M., & Zuber, N. (1979). Drag coefficient and relative velocity in bubbly, droplet or particulate flows. AIChE Journal, 25(5), 843–855. https://doi.org/10.1002/aic.690250513
"""
struct IshiiZuber1979_Bubbles <: AbstractRelativeViscosity end

@inline function compute_relative_viscosity(::IshiiZuber1979_Bubbles, ϕ)
    ϕ = float(ϕ)
    0 <= ϕ <= 1 || return oftype(ϕ, NaN)
    return inv(1 - ϕ)
end

@inline compute_relative_viscosity(::Nothing, ϕ) = one(float(ϕ))

"""
    AbstractRateDependentViscosity

Relative viscosity ``\\eta_r(\\phi, x)`` that also depends on a rate ``x``, such as the capillary number of an
[`AbstractCapillaryViscosity`](@ref GeoParams.MaterialParameters.ConstitutiveRelationships.AbstractCapillaryViscosity)
or the shear rate of [`Vona2011_Crystals`](@ref).
"""
abstract type AbstractRateDependentViscosity <: AbstractRelativeViscosity end

"""
    AbstractCapillaryViscosity

Relative viscosity ``\\eta_r(\\phi, \\mathrm{Ca})`` of a suspension of deformable bubbles, which also
depends on the capillary number
```math
    \\mathrm{Ca} = \\frac{\\eta\\,\\dot{\\gamma}\\,r}{\\Gamma} = \\frac{r\\,\\tau_{II}}{\\Gamma},
```
with bubble radius ``r``, surface tension ``\\Gamma`` and the viscosity ``\\eta`` of the liquid around the
bubbles. The shear rate ``\\dot{\\gamma} = 2\\dot{\\varepsilon}_{II}`` equals ``\\sqrt{2\\,\\dot{\\varepsilon}_{ij}\\dot{\\varepsilon}_{ij}}``
and is the engineering shear rate in simple shear, so ``\\eta\\,\\dot{\\gamma}`` is the second stress invariant
``\\tau_{II}`` of that liquid.
"""
abstract type AbstractCapillaryViscosity <: AbstractRateDependentViscosity end

"""
    Mader2013_Bubbles(; r = 1e-4m, Γ = 0.3J/m^2)

Relative viscosity of a suspension of deformable bubbles in steady flow, the Cross model of Mader et al.
(2013, eqs 27–29),
```math
    \\eta_r = \\eta_{r\\infty} + \\frac{\\eta_{r0} - \\eta_{r\\infty}}{1 + (K\\,\\mathrm{Ca})^2}, \\qquad
    \\eta_{r0} = \\frac{1}{1 - \\phi_\\mathrm{gas}}, \\qquad \\eta_{r\\infty} = (1 - \\phi_\\mathrm{gas})^{5/3}, \\qquad K = \\frac{6}{5},
```
with gas volume fraction ``\\phi_\\mathrm{gas}``, bubble radius `r` and surface tension `Γ`. ``K`` and the exponent 2
are from Llewellin et al. (2002); the limits ``\\eta_{r0}`` and ``\\eta_{r\\infty}`` are those of Llewellin & Manga
(2005, eqs 7–8), who switch between them at ``\\mathrm{Ca} = 1``. The capillary number
``\\mathrm{Ca}`` is defined in [`AbstractCapillaryViscosity`](@ref GeoParams.MaterialParameters.ConstitutiveRelationships.AbstractCapillaryViscosity).
Undeformed bubbles (``\\mathrm{Ca} \\to 0``) stiffen the melt as in [`IshiiZuber1979_Bubbles`](@ref); fully deformed bubbles
(``\\mathrm{Ca} \\to \\infty``) weaken it. For dilute suspensions the law reduces to
``\\eta_r = 1 + \\phi_\\mathrm{gas}(1 - \\tfrac{12}{5}\\mathrm{Ca}^2)/(1 + \\tfrac{36}{25}\\mathrm{Ca}^2)``, eq. 11 of
Llewellin et al. (2002) with the theoretical (Taylor) dilute value ``b = 1``; Llewellin et al. (2002) fit
``b = 9`` to their data.
``\\eta_r`` is `NaN` outside ``0 \\le \\phi_\\mathrm{gas} \\le 1``.

The default `Γ` = 0.3 N/m is the value of the worked example of Mader et al. (2013, §7.3.1, p. 155), the upper
end of the range 0.05–0.3 N/m they give for hydrous to anhydrous silicate melts (§3.1, p. 137). The default
bubble radius `r` = 0.1 mm is a placeholder; set it for each application.

# References
- Mader, H. M., Llewellin, E. W., & Mueller, S. P. (2013). The rheology of two-phase magmas: A review and analysis. Journal of Volcanology and Geothermal Research, 257, 135–158. https://doi.org/10.1016/j.jvolgeores.2013.02.014
- Llewellin, E. W., Mader, H. M., & Wilson, S. D. R. (2002). The constitutive equation and flow dynamics of bubbly magmas. Geophysical Research Letters, 29(24), 2170. https://doi.org/10.1029/2002GL015697
- Llewellin, E. W., & Manga, M. (2005). Bubble suspension rheology and implications for conduit flow. Journal of Volcanology and Geothermal Research, 143(1–3), 205–217. https://doi.org/10.1016/j.jvolgeores.2004.09.018
"""
struct Mader2013_Bubbles{T, U1, U2} <: AbstractCapillaryViscosity
    r::GeoUnit{T, U1}  # bubble radius
    Γ::GeoUnit{T, U2}  # surface tension

    function Mader2013_Bubbles(r, Γ)
        rU, ΓU = convert(GeoUnit, r), convert(GeoUnit, Γ)
        rU.val > 0 || throw(ArgumentError("bubble radius r must be positive, got $r"))
        ΓU.val > 0 || throw(ArgumentError("surface tension Γ must be positive, got $Γ"))
        T = float(promote_type(typeof(rU.val), typeof(ΓU.val)))
        rU, ΓU = GeoUnit(T(rU.val), rU.unit, rU.isdimensional), GeoUnit(T(ΓU.val), ΓU.unit, ΓU.isdimensional)
        return new{T, typeof(rU.unit), typeof(ΓU.unit)}(rU, ΓU)
    end
end
Mader2013_Bubbles(; r = 1.0e-4m, Γ = 0.3J / m^2) = Mader2013_Bubbles(r, Γ)

@inline function relative_viscosity_and_slope(::Mader2013_Bubbles, ϕ, Ca)
    ϕ, Ca = promote(float(ϕ), float(Ca))
    Tc = typeof(ϕ)
    0 <= ϕ <= 1 || return Tc(NaN), Tc(NaN)
    η0 = inv(1 - ϕ)
    η∞ = (1 - ϕ)^(Tc(5) / 3)
    y = (Tc(6) / 5 * Ca)^2
    ηr = η∞ + (η0 - η∞) / (1 + y)
    return ηr, -2 * (η0 - η∞) * y / ((1 + y)^2 * ηr)
end

"""
    Vona2011_Crystals(; aspect_ratio = 6.5, a = 0.06)

Shear-thinning relative viscosity of crystal-bearing basalt after Vona et al. (2011, eq. 23),
```math
    \\eta_r = \\left(1 - \\frac{\\phi_x}{\\phi_m}\\right)^{-2\\,[1 - a \\log_{10}(\\dot{\\gamma}/\\dot{\\gamma}_0)]},
    \\qquad \\dot{\\gamma}_0 = 1\\,\\mathrm{s}^{-1},
```
with crystal volume fraction ``\\phi_x``, shear rate ``\\dot{\\gamma}`` and the maximum packing fraction ``\\phi_m`` of
[`Mueller2010_Crystals`](@ref) for crystals of mean aspect ratio `aspect_ratio`, as Vona et al. (2011, eq. 6) compute it.
Equivalently, the suspension is a power-law fluid,
``\\eta_r = (1 - \\phi_x/\\phi_m)^{-2}\\,(\\dot{\\gamma}/\\dot{\\gamma}_0)^{n - 1}`` with flow index
``n = 1 + 2a\\log_{10}(1 - \\phi_x/\\phi_m)`` (eq. 20), whose consistency is the Maron–Pierce law. Vona et al. fit
``a = 0.06`` to Stromboli and Etna basalts with ``0.07 \\le \\phi_x \\le 0.27``, mean aspect ratios 5.6–7.3 and
shear rates 0.02–2.2 s⁻¹ in a concentric-cylinder viscometer. That shear rate of simple shear is
``\\dot{\\gamma} = 2\\dot{\\varepsilon}_{II}``; in [`Suspension`](@ref) it is the shear rate of the whole suspension.
The default aspect ratio is the representative value of their Fig. 10b. ``\\eta_r`` is `NaN` outside
``0 \\le \\phi_x \\le \\phi_m``.

# Example
```julia-repl
julia> compute_relative_viscosity(Vona2011_Crystals(aspect_ratio = 7.25), 0.272, 0.02)
17.515418342481468
```

# References
- Vona, A., Romano, C., Dingwell, D. B., & Giordano, D. (2011). The rheology of crystal-bearing basaltic magmas from Stromboli and Etna. Geochimica et Cosmochimica Acta, 75(11), 3214–3236. https://doi.org/10.1016/j.gca.2011.03.031
"""
struct Vona2011_Crystals{T, U} <: AbstractRateDependentViscosity
    ϕ_max::T           # maximum packing fraction
    a::T               # shear-thinning coefficient
    γ̇0::GeoUnit{T, U}  # reference shear rate

    function Vona2011_Crystals(ϕ_max, a, γ̇0)
        0 < ϕ_max <= 1 || throw(ArgumentError("Vona2011_Crystals: ϕ_max must lie in (0, 1], got $ϕ_max"))
        a >= 0 || throw(ArgumentError("Vona2011_Crystals: a must be non-negative, got $a"))
        γU = convert(GeoUnit, γ̇0)
        T = float(promote_type(typeof(ϕ_max), typeof(a), typeof(γU.val)))
        γU = GeoUnit(T(γU.val), γU.unit, γU.isdimensional)
        return new{T, typeof(γU.unit)}(ϕ_max, a, γU)
    end
end
Vona2011_Crystals(; aspect_ratio = 6.5, a = 0.06) = Vona2011_Crystals(Mueller2010_Crystals(; aspect_ratio).ϕ_max, a, 1.0 / s)

@inline function relative_viscosity_and_slope(f::Vona2011_Crystals, ϕ, γ̇)
    ϕ, γ̇ = promote(float(ϕ), float(γ̇))
    Tc = typeof(ϕ)
    @unpack_like γ̇ Tc γ̇0 = f
    ϕ_max, a = convert(Tc, f.ϕ_max), convert(Tc, f.a)
    0 <= ϕ <= ϕ_max || return Tc(NaN), Tc(NaN)
    x = 1 - ϕ / ϕ_max
    # the floor keeps η_r finite at γ̇ = 0, so that η_r γ̇ vanishes there
    ηr = x^(-2 * (1 - a * log10(max(γ̇ / γ̇0, floatmin(Tc)))))
    return ηr, 2 * a * log(x) / log(Tc(10))
end

"""
    compute_relative_viscosity(f::AbstractRateDependentViscosity, ϕ, x)

Relative viscosity of the law `f` at volume fraction `ϕ` and rate `x`: the capillary number for bubble laws
such as [`Mader2013_Bubbles`](@ref), and the shear rate ``\\dot{\\gamma} = 2\\dot{\\varepsilon}_{II}`` for [`Vona2011_Crystals`](@ref).

# Example
```julia-repl
julia> compute_relative_viscosity(Mader2013_Bubbles(), 0.5, 0.0)
2.0
```
"""
@inline compute_relative_viscosity(f::AbstractRateDependentViscosity, ϕ, x) = first(relative_viscosity_and_slope(f, ϕ, x))

compute_relative_viscosity(f::AbstractRateDependentViscosity, ϕ) =
    throw(ArgumentError("$(nameof(typeof(f))) depends on the strain rate; call compute_relative_viscosity(f, ϕ, x)"))

# η_r of a law that does not depend on the rate x, with slope d ln η_r / d ln x = 0
@inline function relative_viscosity_and_slope(f::Union{Nothing, AbstractRelativeViscosity}, ϕ, x)
    ηr = compute_relative_viscosity(f, ϕ)
    return ηr, zero(ηr)
end

# capillary number r τII / Γ of a bubble in liquid at stress TauII
@inline function capillary_number(f::AbstractCapillaryViscosity, TauII)
    Tc = precision_of(TauII)
    @unpack_like TauII Tc r, Γ = f
    return r * TauII / Γ
end

"""
    Suspension(law; crystals = nothing, bubbles = nothing)

Creep law `law` of the melt, stiffened by suspended crystals and bubbles. As in the three-phase model of
Truby et al. (2015, eq. 3.2), the crystals are suspended in the bubbly melt:
```math
    \\dot{\\varepsilon}_{II}(\\tau_{II}) = \\frac{\\dot{\\varepsilon}_{II}^{\\mathrm{law}}(\\tau_{II})}{\\eta_{r,x}(\\phi_x)\\,\\eta_{r,g}(\\phi_b)},
    \\qquad \\phi_b = \\frac{\\phi_\\mathrm{gas}}{1 - \\phi_x},
```
where ``\\phi_b`` is the gas fraction of the bubbly melt (eqs 3.4, 3.5 and 3.7), so a Newtonian melt viscosity
``\\eta_\\mathrm{melt}`` becomes ``\\eta = \\eta_\\mathrm{melt}\\,\\eta_{r,x}\\,\\eta_{r,g}``. Truby et al. validate this
for spherical particles in the Maron–Pierce law and undeformed bubbles ([`IshiiZuber1979_Bubbles`](@ref)) up to
``\\phi_x \\approx 0.5`` and ``\\phi_b \\approx 0.3``. `crystals` and `bubbles` are relative viscosities such as
[`Roscoe1952_Crystals`](@ref), [`Mueller2010_Crystals`](@ref) and [`IshiiZuber1979_Bubbles`](@ref); `nothing` means no effect.

The relative viscosities may also depend on the strain rate ``\\dot{\\varepsilon}_{II}`` of the suspension.
Crystal laws such as [`Vona2011_Crystals`](@ref) depend on its shear rate ``2\\dot{\\varepsilon}_{II}``. Bubble laws such as
[`Mader2013_Bubbles`](@ref) depend on the capillary number, which the stress of the melt sets at the strain rate
``\\sqrt{\\eta_{r,x}/(1 - \\phi_x)}\\,\\dot{\\varepsilon}_{II}`` of the bubbly melt, the average melt strain rate
that dissipates the same energy as the crystal suspension (Chateau et al. 2008). Then `compute_τII` stays explicit, and
`compute_εII` solves for ``\\dot{\\varepsilon}_{II}`` by Newton iteration and returns `NaN` if it does not converge. At
high gas fractions the stress of capillary bubble laws can decrease with strain rate, so the strain rate at a
given stress need not be unique.

The crystal and gas volume fractions of the whole suspension are passed at runtime as `ϕ_x` and `ϕ_gas`, the
same keywords that [`ThreePhase_Density`](@ref GeoParams.MaterialParameters.Density.ThreePhase_Density) uses, e.g. `compute_εII(law, τII; T, ϕ_x = 0.3, ϕ_gas = 0.1)`.
Both are required. All keyword arguments are also passed on to `law`.

# Example
```julia-repl
julia> law = Suspension(LinearViscous(η = 1e5Pa*s); crystals = Roscoe1952_Crystals(ϕ_max = 1), bubbles = IshiiZuber1979_Bubbles());

julia> compute_τII(law, 1.0; ϕ_x = 0.2, ϕ_gas = 0.5) / compute_τII(law.law, 1.0)
4.658474953124561
```

# References
- Chateau, X., Ovarlez, G., & Trung, K. L. (2008). Homogenization approach to the behavior of suspensions of noncolloidal particles in yield stress fluids. Journal of Rheology, 52(2), 489–506. https://doi.org/10.1122/1.2838254
- Truby, J. M., Mueller, S. P., Llewellin, E. W., & Mader, H. M. (2015). The rheology of three-phase suspensions at low bubble capillary number. Proceedings of the Royal Society A, 471(2173), 20140557. https://doi.org/10.1098/rspa.2014.0557
"""
struct Suspension{T, L <: AbstractCreepLaw{T}, X, G} <: AbstractCreepLaw{T}
    law::L       # creep law of the melt
    crystals::X  # relative viscosity of the crystal suspension, or nothing
    bubbles::G   # relative viscosity of the bubble suspension, or nothing

    function Suspension(law::AbstractCreepLaw{T}, crystals, bubbles) where {T}
        for f in (crystals, bubbles)
            f isa Union{Nothing, AbstractRelativeViscosity} ||
                throw(ArgumentError("Suspension: crystals and bubbles must be relative viscosities or nothing, got $f"))
        end
        crystals isa AbstractCapillaryViscosity &&
            throw(ArgumentError("Suspension: crystals cannot depend on the capillary number, got $crystals"))
        return new{T, typeof(law), typeof(crystals), typeof(bubbles)}(law, crystals, bubbles)
    end
end
Suspension(law::AbstractCreepLaw; crystals = nothing, bubbles = nothing) = Suspension(law, crystals, bubbles)

function param_info(a::Suspension)
    inner = param_info(a.law)
    return MaterialParamsInfo(;
        Equation = L"\eta = \eta_\mathrm{melt}\,\eta_{r,x}(\phi_x)\,\eta_{r,g}(\phi_\mathrm{gas}/(1 - \phi_x))",
        Comment = inner.Comment,
        BibTex_Reference = inner.BibTex_Reference,
    )
end

"""
    compute_relative_viscosity(f::AbstractRelativeViscosity, ϕ)
    compute_relative_viscosity(a::Suspension; ϕ_x, ϕ_gas)
    compute_relative_viscosity(a::Suspension, args)

Relative viscosity ``\\eta_r = \\eta / \\eta_\\mathrm{melt}`` of the law `f` at volume fraction `ϕ`, or of
the suspension `a`, which is ``\\eta_{r,x}(\\phi_x)\\,\\eta_{r,g}(\\phi_\\mathrm{gas}/(1 - \\phi_x))``. `args` is a
`NamedTuple` that must contain `ϕ_x` and `ϕ_gas`; its other entries are ignored.

# Example
```julia-repl
julia> compute_relative_viscosity(IshiiZuber1979_Bubbles(), 0.5)
2.0
```
"""
@inline function compute_relative_viscosity(a::Suspension; ϕ_x, ϕ_gas, kwargs...)
    return compute_relative_viscosity(a.crystals, ϕ_x) * compute_relative_viscosity(a.bubbles, ϕ_gas / (1 - ϕ_x))
end
@inline compute_relative_viscosity(a::Suspension, args) = compute_relative_viscosity(a; args...)

# η_r in the precision `Tc` of the stress or strain rate
@inline function compute_relative_viscosity(a::Suspension, ::Type{Tc}, ϕ_x, ϕ_gas) where {Tc}
    return compute_relative_viscosity(a; ϕ_x = convert_precision(Tc, ϕ_x), ϕ_gas = convert_precision(Tc, ϕ_gas))
end

# η_r = η_r,x η_r,g at strain rate EpsII and, if `slope` is true, d ln η_r / d ln ε̇ (otherwise 0). Crystals are
# sheared at 2ε̇; bubbles fill the fraction ϕ_gas / (1 - ϕ_x) of the bubbly melt, sheared at √(η_r,x / (1 - ϕ_x)) ε̇.
@inline function relative_viscosity_and_slope(a::Suspension, EpsII, ϕ_x, ϕ_gas, slope; kwargs...)
    ηrx, sx = relative_viscosity_and_slope(a.crystals, ϕ_x, 2 * EpsII)
    ηrg, sg = bubble_viscosity_and_slope(a.bubbles, a.law, sqrt(ηrx / (1 - ϕ_x)) * EpsII, ϕ_gas / (1 - ϕ_x), slope; kwargs...)
    return ηrx * ηrg, sx + sg * (1 + sx / 2)
end

# η_r,g of bubbles at fraction ϕ in melt sheared at strain rate ε and, if `slope` is true, d ln η_r,g / d ln ε
@inline function bubble_viscosity_and_slope(f, law, ε, ϕ, slope; kwargs...)
    ηr = compute_relative_viscosity(f, ϕ)
    return ηr, zero(ηr)
end

@inline function bubble_viscosity_and_slope(f::AbstractCapillaryViscosity, law, ε, ϕ, slope; kwargs...)
    τm = compute_τII(law, ε; kwargs...)
    ηr, s = relative_viscosity_and_slope(f, ϕ, capillary_number(f, τm))
    # at τm = 0 the slope of η_r,g vanishes and d ln τm / d ln ε is undefined
    (slope && !iszero(τm)) || return ηr, zero(s)
    return ηr, s * ε * dτII_dεII(law, ε; kwargs...) / τm
end

@inline function compute_τII(a::Suspension, EpsII; ϕ_x, ϕ_gas, kwargs...)
    Tc = precision_of(EpsII)
    ϕx, ϕg = convert_precision(Tc, ϕ_x), convert_precision(Tc, ϕ_gas)
    ηr, _ = relative_viscosity_and_slope(a, EpsII, ϕx, ϕg, false; ϕ_x, ϕ_gas, kwargs...)
    return compute_τII(a.law, ηr * EpsII; ϕ_x, ϕ_gas, kwargs...)
end

@inline function dτII_dεII(a::Suspension, EpsII; ϕ_x, ϕ_gas, kwargs...)
    Tc = precision_of(EpsII)
    ϕx, ϕg = convert_precision(Tc, ϕ_x), convert_precision(Tc, ϕ_gas)
    ηr, s = relative_viscosity_and_slope(a, EpsII, ϕx, ϕg, true; ϕ_x, ϕ_gas, kwargs...)
    return ηr * (1 + s) * dτII_dεII(a.law, ηr * EpsII; ϕ_x, ϕ_gas, kwargs...)
end

@inline function compute_εII(a::Suspension, TauII; ϕ_x, ϕ_gas, kwargs...)
    ηr = compute_relative_viscosity(a, precision_of(TauII), ϕ_x, ϕ_gas)
    return compute_εII(a.law, TauII; ϕ_x, ϕ_gas, kwargs...) / ηr
end

@inline function dεII_dτII(a::Suspension, TauII; ϕ_x, ϕ_gas, kwargs...)
    ηr = compute_relative_viscosity(a, precision_of(TauII), ϕ_x, ϕ_gas)
    return dεII_dτII(a.law, TauII; ϕ_x, ϕ_gas, kwargs...) / ηr
end

# The melt law fixes η_r(ε̇) ε̇ = ε̇_law(τII); Newton iteration on ln ε̇ solves it for ε̇.
# Suspensions whose relative viscosity depends on the strain rate
@inline function compute_εII(
        a::Union{
            Suspension{T, L, <:AbstractRateDependentViscosity} where {T, L <: AbstractCreepLaw{T}},
            Suspension{T, L, X, <:AbstractRateDependentViscosity} where {T, L <: AbstractCreepLaw{T}, X},
        },
        TauII; ϕ_x, ϕ_gas, kwargs...
    )
    Tc = precision_of(TauII)
    iszero(TauII) && return zero(TauII)
    TauII < 0 && return -compute_εII(a, -TauII; ϕ_x, ϕ_gas, kwargs...)   # τ(ε̇) is odd
    ϕx, ϕg = convert_precision(Tc, ϕ_x), convert_precision(Tc, ϕ_gas)
    lnε = log(compute_εII(a.law, TauII; ϕ_x, ϕ_gas, kwargs...))  # ln(η_r ε̇)
    x = lnε
    for _ in 1:50
        ηr, s = relative_viscosity_and_slope(a, exp(x), ϕx, ϕg, true; ϕ_x, ϕ_gas, kwargs...)
        dx = (log(ηr) + x - lnε) / (1 + s)
        x -= dx
        abs(dx) < sqrt(eps(Tc)) && return exp(x)
    end
    return Tc(NaN)
end

@inline function dεII_dτII(
        a::Union{
            Suspension{T, L, <:AbstractRateDependentViscosity} where {T, L <: AbstractCreepLaw{T}},
            Suspension{T, L, X, <:AbstractRateDependentViscosity} where {T, L <: AbstractCreepLaw{T}, X},
        },
        TauII; kwargs...
    )
    return inv(dτII_dεII(a, compute_εII(a, TauII; kwargs...); kwargs...))
end

function compute_εII!(
        EpsII::AbstractArray{_T, N}, a::Suspension, TauII::AbstractArray; ϕ_x, ϕ_gas, kwargs...
    ) where {N, _T}
    for i in each_argument_index(EpsII, TauII, ϕ_x, ϕ_gas, values(kwargs)...)
        args_i = map(x -> argument_at(x, i), values(kwargs))
        EpsII[i] = compute_εII(
            a, convert_precision(_T, TauII[i]); ϕ_x = argument_at(ϕ_x, i), ϕ_gas = argument_at(ϕ_gas, i), args_i...
        )
    end
    return nothing
end

function compute_τII!(
        TauII::AbstractArray{_T, N}, a::Suspension, EpsII::AbstractArray; ϕ_x, ϕ_gas, kwargs...
    ) where {N, _T}
    for i in each_argument_index(TauII, EpsII, ϕ_x, ϕ_gas, values(kwargs)...)
        args_i = map(x -> argument_at(x, i), values(kwargs))
        TauII[i] = compute_τII(
            a, convert_precision(_T, EpsII[i]); ϕ_x = argument_at(ϕ_x, i), ϕ_gas = argument_at(ϕ_gas, i), args_i...
        )
    end
    return nothing
end

function show(io::IO, g::Suspension)
    return print(io, "Suspension: crystals = $(nameof(typeof(g.crystals))), bubbles = $(nameof(typeof(g.bubbles))), law = $(g.law)")
end
