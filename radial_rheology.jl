using GeoParams
using SpecialFunctions
using Plots

function viscosity(P, T, strain; B=1e-28, n=4, E=223e3, V=0e0, R=8.3145)
    F = 1 / (2^((n-1)/n)*3^((n+1)/2/n))
    η = F*B^(-1/n)*strain^((1-n)/n)*exp((E+P*V)/n/R/T)
    return max(min(η, 1e25), 1e16)
end

## TEMPERATURE AND PRESSURE RADIAL PROFILES
npoints = 600*2
yr = 365.25 * 3600 * 24
z = LinRange(0, 600e3, npoints)
t = 100 * yr
κ = 3/3300/1200
Tm = 2500
Tp = 1900
T = @. (Tm-Tp)/2890e3*z + Tp + 273
@. T *= erf(z*0.5/sqrt(κ*t))  # half-space cooling
P = @. z*9.81*3300

## RHEOLOGICAL PROPERTIES

# LITHOSPHERE
dis_litho = DislocationCreep(A=10^-28.0, E=223e3, n=4, V=0.0) # wet quartzite
lithosphere = SetMaterialParams(;
    Name="Lithosphere",
    Phase=1,
    CreepLaws=(dis_litho,),
    Density=ConstantDensity(; ρ=3300kg / m^3),
)

# MANTLE
dis_mantle = DislocationCreep(A=10^-15.96, E=530e3, n=3.5, V=13e-6) # dry olivine
dif_mantle = DiffusionCreep(A=10^-8.16, E=375e3, n=1, V=6e-6) # dry olivine
lithospheric_mantle = SetMaterialParams(;
    Name="Mantle",
    Phase=2,
    CreepLaws=(dis_mantle, dif_mantle),
    Density=ConstantDensity(; ρ=2900kg / m^3),
)

# phase tuple
material_phases = (lithosphere, lithospheric_mantle)

# define radial phase profile
Phases = fill(2, length(T))
idx = z .< 120e3
Phases[idx] .= 1

# arguments to compute the viscosities
args = (; P=P, T=T, f=1e0, d=1e0)

# compute viscosity profile
εII_bg = 1e-14 # background strain rate
εII = fill(εII_bg, npoints) # background strain rate

η = similar(T)
computeViscosity_εII!(
    η,
    material_phases,
    εII,
    Phases,
    args;
)

η_hardcoded = ones(npoints)
for i in eachindex(T)
    if z[i] > 120e3 # lithospheric mantle
        # nu_dis = viscosity(P[i], T[i], εII[i]; B=10^-15.96, n=3.5, E=530e3, V=13e-6)
        # nu_dif = viscosity(P[i], T[i], εII[i]; B=10^-8.16, n=1, E=375e3, V=6e-6)
        # η_hardcoded[i] = 1/( 1/nu_dis + 1/nu_dif)
    else
        nu_dis = viscosity(P[i], T[i], εII[i]; B=1e-28, n=4, E=223e3, V=6e-6)
        η_hardcoded[i] = nu_dis
    end
end

plot(log10.(η), -z)
plot!(log10.(η_hardcoded), -z)
