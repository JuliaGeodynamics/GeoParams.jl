import Pkg; Pkg.activate(".")
using GeoParams

function Local_Iterations(v, args)
    # Physics
    εII_ve = 1e-14          # second invariant of effective viscoelastic strain rate, 1/s
    εII_vis = εII_ve
    
    # Initial guess
    η_ve = computeViscosity_EpsII(εII_ve, v, args) # guess: it computes 1.0 / (1.0 / η + 1.0 / (μ * dt))
    τII = 2 * η_ve * εII_ve            # guess
    
    # Local Iterations
    iter = 0
    while iter < 5  # Newton
        iter = iter + 1
        f = εII_ve - strain_rate_circuit(τII, v, args)
        dfdτII = 0.0 - dεII_dτII(τII, v[1], args) - dεII_dτII(τII, v[2], args) # yet to be generalized
        τII = τII - f / dfdτII
    end

end

# Physics
εII_ve = 1e-14          # second invariant of effective viscoelastic strain rate, 1/s
μ = 1e11           # shear modolus, Pa
A = 1.67e-24       # pre-exponentialfactor, Pa^(-n) s^(-1) 
n = 3.3            # power-law exponent, []
E = 187e3          # activation energy, J/mol 
R = 8.3145         # universal gas constant, J/mol/K
T = 500 + 273.15     # temperature, Ks
# numerics 
dt = 1e3 * (365.25 * 24 * 3600)   # time step, s
nt = 5                   # number of iterations, []
# Initialiazation
εII_vis = εII_ve
iter = 0

# Define composite rheology
v_dis = DislocationCreep(A=A, E=E, n=n, V=0.0)
v_el = ConstantElasticity(G=μ*Pa)
v = (v_el, v_el) # dislocation and elasticity are in serie
args = (P=0.0, T=T, f=1.0, dt=dt)

Local_Iterations(v, args)
