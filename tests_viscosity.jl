using SpecialFunctions
using Plots
using Pkg: Pkg;
Pkg.activate(".");
using GeoParams

yr = 365.25 * 3600 * 24
z = LinRange(0, 2890e3, 2890)
t = 100 * yr
κ = 3 / 3300 / 1200
Tm = 2500
Tp = 1900
T = @. (Tm - Tp) / 2890e3 * z + Tp + 273
@. T *= erf(z * 0.5 / sqrt(κ * t))  # half-space cooling
P = @. z * 9.81 * 3300

# plot(T, -z)
# display(gcf())

function miguel(P, T, strain; B=1e-28, n=4, E=223e3, V=0e0, R=8.3145)
    F = 1 / (2^((n - 1) / n) * 3^((n + 1) / 2 / n))
    a = F * B^(-1 / n) * strain^((1 - n) / n) * exp((E + P * V) / n / R / T)
    return a
end

# Physics
εII = 1e-14          # second invariant of effective viscoelastic strain rate, 1/s
μ = 1e11           # shear modolus, Pa
A = 10^-8.65       # pre-exponentialfactor, Pa^(-n) s^(-1) 
n = 3.3            # power-law exponent, []
E = 335e3          # activation energy, J/mol 
R = 8.3145         # universal gas constant, J/mol/K
V = 4e-6
# numerics 
dt = 1e3 * yr   # time step, s

# Define composite rheology
v_dis = DislocationCreep(; A=10^-8.65, E=335e3, n=1, V=4e-6) # wet olivine
v_dif = DiffusionCreep(; A=10^-15.81, E=480e3, n=3.5, V=10e-6)  # wet olivine
v_dis = DislocationCreep(; A=10^-28.0, E=223e3, n=4, V=0.0)  # wet quartzite

args = (P=P, T=T, f=1.0, d=4e-3)

visc_dif = similar(T)
visc_dis = similar(T)
visc_eff = similar(T)
visc_miguel_dis = similar(T)
to = TimerOutput()
for i in eachindex(T)
    argsi = (P=P[i], T=T[i], f=1.0, d=4e-3)
    @timeit to "dislocation" visc_dif[i] = computeViscosity_εII(v_dif, εII, argsi)
    visc_dis[i] = computeViscosity_εII(v_dis, εII, argsi)
    visc_eff[i] = 1 / (1 / visc_dif[i] + 1 / visc_dis[i])
    @timeit to "miguel" visc_miguel_dis[i] = miguel(P[i], T[i], εII)
end
to
plot(log10.(visc_dif), -z .* 1e-3)
plot!(log10.(visc_dis), -z .* 1e-3)
plot!(log10.(visc_eff), -z .* 1e-3)
plot!(log10.(visc_miguel_dis), -z .* 1e-3)

ProfileCanvas.@profview for _ in 1:100000
    miguel(P[i], T[i], εII)
end
ProfileCanvas.@profview for _ in 1:100000
    computeViscosity_εII(v_dif, εII, argsi)
end
@btime computeViscosity_εII($v_dif, $εII, $argsi)
@btime miguel($(P[i]), $(T[i]), $εII)

# KelvinVoigt(v_el, v_dis)
# vk = KelvinVoigt(v_el, (v_dis, v_dis))

# vk = KelvinVoigt(v_el, v_dis)
# compute_εII(vk, rand(), args)

# Local_Iterations(v, args)

# compute_τII(v, εII_vis, args)

# @btime KelvinVoigt($v_el, ($v_dis, $v_dis))
# @btime KelvinVoigt($v_el, $(v_dis, v_dis))
