"""
bulk_viscosity(η_shear, ϕ, C, R, λ, P_eff)

# η_shear: shear viscosity
# ϕ: porosity
# C: shear viscosity ratio (value, dimensionless)
# R: compaction strength ratio (value, dimensionless)
# λ: effective pressure transition zone
# P_eff: effective Pressure
"""
function bulk_viscosity(η_shear, ϕ, C, R, λ, P_eff)

    η_bulk = η_shear / ϕ * C * (1.0 + 0.5 * (inv(R) - 1.0) * (1.0 + tanh(-P_eff * inv(λ))))

    return η_bulk
end


### Georg Adjoint 2 pahse paper,
### Ludo's 2019 paper has it quite nicely also with a calc for R
