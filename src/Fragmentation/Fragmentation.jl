module Fragmentation

using Unitful
using ..Units
using GeoParams: AbstractMaterialParam
import GeoParams: DruckerPrager_regularised

export AbstractFragmentation,
    CriticalPorosityFragmentation,
    StrainRateFragmentation,
    compute_fragmentation

"""
    AbstractFragmentation <: AbstractMaterialParam

Supertype of magma fragmentation criteria. [`compute_fragmentation`](@ref) returns the
dimensionless ratio of a criterion's driving quantity to its threshold; the magma is
fragmented where this ratio is at least 1. Implemented criteria are
[`CriticalPorosityFragmentation`](@ref) and [`StrainRateFragmentation`](@ref).
"""
abstract type AbstractFragmentation <: AbstractMaterialParam end

"""
    compute_fragmentation(f::AbstractFragmentation, args)

Fragmentation ratio of criterion `f` at a point, with the inputs of `f` passed in the
`NamedTuple` `args`; entries that `f` does not use are ignored. The magma is fragmented
where the ratio is at least 1. The keyword form `f(; args...)` is equivalent. The ratio
is dimensionless and therefore the same for dimensional and nondimensionalized inputs.
"""
@inline compute_fragmentation(f::AbstractFragmentation, args) = f(; args...)

"""
    CriticalPorosityFragmentation(; ϕ_crit = 0.75)

Fragmentation at a critical gas volume fraction. [`compute_fragmentation`](@ref) returns
```math
    \\frac{\\phi_\\mathrm{gas}}{\\phi_\\mathrm{crit}},
```
with the gas volume fraction `ϕ_gas` of the magma as input. Sparks (1978, abstract and
p. 28) finds void-to-liquid volume ratios of 3:1 to 5:1 in pumice, i.e. vesicularities
of 0.75 to 0.83 just before disruption; the default `ϕ_crit = 0.75` is the lower bound,
the value Papale (1999, Table 1, p. 426) uses for his critical-porosity runs.

# References
- Sparks, R. S. J. (1978). The dynamics of bubble formation and growth in magmas: A review and analysis. Journal of Volcanology and Geothermal Research, 3(1–2), 1–37. https://doi.org/10.1016/0377-0273(78)90002-1
- Papale, P. (1999). Strain-induced magma fragmentation in explosive eruptions. Nature, 397, 425–428. https://doi.org/10.1038/17109
"""
struct CriticalPorosityFragmentation{T} <: AbstractFragmentation
    ϕ_crit::T  # critical gas volume fraction

    function CriticalPorosityFragmentation{T}(ϕ_crit) where {T}
        0 < ϕ_crit < 1 || throw(ArgumentError("CriticalPorosityFragmentation: ϕ_crit must lie in (0, 1), got $ϕ_crit"))
        return new{T}(ϕ_crit)
    end
end
CriticalPorosityFragmentation(ϕ_crit) = CriticalPorosityFragmentation{float(typeof(ϕ_crit))}(ϕ_crit)
CriticalPorosityFragmentation(; ϕ_crit = 0.75) = CriticalPorosityFragmentation(ϕ_crit)

@inline (f::CriticalPorosityFragmentation)(; ϕ_gas, kwargs...) =
    ϕ_gas / convert_precision(precision_of(ϕ_gas), f.ϕ_crit)

"""
    StrainRateFragmentation(; k = 0.01, G = 1e10Pa)

Strain-induced fragmentation of Papale (1999, eq. 1, p. 426): magma breaks when its
elongational strain rate exceeds a fraction `k` of the inverse Maxwell relaxation time
``G_\\infty/\\eta``. [`compute_fragmentation`](@ref) returns
```math
    \\frac{\\dot{\\varepsilon}_e\\,\\eta}{k\\,G_\\infty}, \\qquad
    \\dot{\\varepsilon}_e = \\frac{2}{\\sqrt{3}}\\,\\dot{\\varepsilon}_{II},
```
with inputs `εII`, the second invariant of the deviatoric strain rate, and `η`, the
Newtonian (zero-frequency) shear viscosity of the magma.

Papale states eq. 1 for the elongational strain rate ``\\partial v_z/\\partial z`` of a
pseudofluid with the mixture viscosity of melt, crystals and bubbles, and keeps `k` and
`G` at their pure-melt values; pass the suspension viscosity (e.g. that of a
[`Suspension`](@ref GeoParams.MaterialParameters.ConstitutiveRelationships.Suspension) law) as `η` to apply eq. 1 as published, or the melt viscosity to
use the melt relaxation time of Gonnermann (2015, sect. 3). ``\\dot{\\varepsilon}_e`` is
the von Mises equivalent strain rate, which equals the elongation rate in incompressible
uniaxial elongation, the geometry of the fiber-elongation experiments that calibrate `k`.
In the compressible one-dimensional flow of Papale's conduit model,
``\\dot{\\varepsilon}_{II} = (\\partial v_z/\\partial z)/\\sqrt{3}``, so
``\\dot{\\varepsilon}_e`` is 2/3 of the elongation rate there.

Papale (1999, p. 426) gives ``k \\approx 0.01`` and ``G_\\infty`` = 3–30 GPa; his
simulations use 25 GPa and change by less than 10% for a tenfold change of ``k\\,G_\\infty``.
`G` is the unrelaxed shear modulus ``G_\\infty`` (at infinite frequency). Its default `G = 1e10Pa` is the
composition-independent estimate
``\\log_{10} G_\\infty\\,(\\mathrm{Pa}) = 10 \\pm 0.5`` of Dingwell & Webb (1989, p. 510),
also used by Gonnermann (2015, sect. 3). `G` must carry stress units, so that
`nondimensionalize` scales it like `η` and `εII`.

# References
- Papale, P. (1999). Strain-induced magma fragmentation in explosive eruptions. Nature, 397, 425–428. https://doi.org/10.1038/17109
- Dingwell, D. B., & Webb, S. L. (1989). Structural relaxation in silicate melts and non-Newtonian melt rheology in geologic processes. Physics and Chemistry of Minerals, 16, 508–516. https://doi.org/10.1007/BF00197020
- Gonnermann, H. M. (2015). Magma fragmentation. Annual Review of Earth and Planetary Sciences, 43, 431–458. https://doi.org/10.1146/annurev-earth-060614-105206
"""
struct StrainRateFragmentation{T, U} <: AbstractFragmentation
    k::T                # ratio of critical strain rate to inverse relaxation time
    G::GeoUnit{T, U}    # unrelaxed shear modulus G∞

    function StrainRateFragmentation(k, G)
        GU = convert(GeoUnit, G isa Quantity ? uconvert(Pa, G) : G)
        k > 0 || throw(ArgumentError("StrainRateFragmentation: k must be positive, got $k"))
        dimension(GU.unit) == dimension(Pa) ||
            throw(ArgumentError("StrainRateFragmentation: G must carry stress units, got $G"))
        GU.val > 0 || throw(ArgumentError("StrainRateFragmentation: G must be positive, got $G"))
        T = float(promote_type(typeof(k), typeof(GU.val)))
        return new{T, typeof(GU.unit)}(k, GeoUnit(T(GU.val), GU.unit, GU.isdimensional))
    end
end
StrainRateFragmentation(; k = 0.01, G = 1.0e10Pa) = StrainRateFragmentation(k, G)

@inline function (f::StrainRateFragmentation)(; εII, η, kwargs...)
    Tc = precision_of((εII, η))
    @unpack_like η Tc G = f
    return 2 * εII * η / (sqrt(Tc(3)) * convert_precision(Tc, f.k) * G)
end

"""
    DruckerPrager_regularised(f::StrainRateFragmentation; η_vp)

Yield surface equivalent to the strain-rate criterion `f`: with ``\\tau_{II} = 2\\eta\\dot{\\varepsilon}_{II}``
the criterion reads ``\\tau_{II} \\ge \\sqrt{3}\\,k\\,G_\\infty``, a von Mises surface, i.e. a regularized
Drucker–Prager surface with `ϕ = 0`, `Ψ = 0` and cohesion `C` ``= \\sqrt{3}\\,k\\,G_\\infty``. In a viscoelastic
model whose magma rheology includes ``G_\\infty`` as elastic modulus, the solver's plastic return mapping then
detects fragmentation on the actual stress. `η_vp` regularizes the localization of the failure; it must be small
compared with the magma viscosity, or the stress overshoots the yield stress. Build it from
a dimensional `f` and nondimensionalize the result.
"""
function DruckerPrager_regularised(f::StrainRateFragmentation; η_vp)
    return DruckerPrager_regularised(; ϕ = 0, Ψ = 0, C = sqrt(3) * f.k * UnitValue(f.G), η_vp)
end

end
