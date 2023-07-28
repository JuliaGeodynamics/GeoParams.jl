






### Als nächstes versuchen die Farben in Thorstens plot nach zu bilden
### ODER die A values zurück umrechnen und gucken ob sie zu Thorstens values passen oder nicht







using Pkg
Pkg.activate(".")
using GeoParams
using GLMakie
using CairoMakie

#x1 = SetGrainBoundarySliding("Dry Olivine < 1523K | Hirth and Kohlstedt (2003)")
#x1 = SetGrainBoundarySliding("Dry Olivine >= 1523K | Hirth and Kohlstedt (2003)")
x1 = SetGrainBoundarySliding("TEST GBS")
x1 = remove_tensor_correction(x1)
#x2 = SetNonLinearPeierlsCreep("Dry Olivine | Mei et al. (2010)")
x2 = SetNonLinearPeierlsCreep("TEST PEIERLS")
x2 = remove_tensor_correction(x2)
#x3 = SetDislocationCreep("Dry Olivine | Hirth & Kohlstedt (2003)")
x3 = SetDislocationCreep("TEST DISL")
x3 = remove_tensor_correction(x3)
#x4 = SetDiffusionCreep("Dry Olivine | Hirth & Kohlstedt (2003)")
x4 = SetDiffusionCreep("TEST DIFF")
x4 = remove_tensor_correction(x4)
c2 = CompositeRheology(x1, x2, x3, x4)

# angepasste values damit der plot UNGEFÄHR passt. rausbekommen wie man von ursprünglichen values auf angepasste values kommt und warum.

args = (P=1.0e9,)
T = (0.0,1523.0)
σ = (1.0e-1, 1.0e4)
d = (1.0e-6, 1.0e-1)

p = PlotDeformationMap(c2, args=args, d=d, σ=σ, T=T, grainsize=true)

#=
function compute_εtestnoinvar(
    a::GrainBoundarySliding, TauII::_T; T=one(precision(a)), P=zero(precision(a)), d=one(precision(a)), args...
) where {_T}
    @unpack_val n, p, A, E, V, R = a
    FT, FE = a.FT, a.FE
    @show n
    @show p
    @show A
    @show E
    @show V
    @show R
    @show FT
    @show FE

    ε = fastpow(TauII / fastpow(A, -1.0 / n), n) *
        fastpow(d, p) *
        exp(-(E + P * V) / (R * T))

    return ε 
end=#