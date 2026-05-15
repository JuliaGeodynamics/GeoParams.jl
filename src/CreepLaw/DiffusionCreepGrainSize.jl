export DiffusionCreepGrainSize,
    remove_tensor_correction,
    compute_εII,
    compute_τII,
    compute_grain_size_reduction,
    compute_growth_rate

#=----Diffusion Creep---
Defines diffusion creep law parameters

r is the water-fugacity exponent
p is the negative defined grain size exponent (value is set negative in the calculation of σ and ε)
A is a   material specific rheological parameter
E is the activation energy
V is the activation volume
R is the universal gas constant
Apparatus defines the appartus type that shall be recreated (Axial Compression, Simple Shear, Invariant)
=#
"""
    DiffusionCreepGrainSize(; Name="", d0=5e-3m, n=1NoUnits, r=0NoUnits, p=-3NoUnits,
        A=1.5e6Pa^(-n-r)/s*m^(-p), Ag=1.5e6Pa^(-rg)/s*m^(-p),
        Qg=134.0e3J/mol, rg=1.38NoUnits, E=500.0e3J/mol,
        V=24.0e-6m^3/mol, R=8.3145J/mol/K, λg=0.015NoUnits,
        γd=1J/m^2, c=π*NoUnits, Apparatus=AxialCompression)

Defines the parameters of a diffusion creep law with grain-size evolution.

The diffusion creep strain-rate law is
```math
    \\dot{\\varepsilon}_{II}
        = A \\tau_{II}^{n} d^{p} f_\\mathrm{H2O}^{r}
          \\exp\\left(-\\frac{E + PV}{RT}\\right)
```
where ``\\tau_{II}`` and ``\\dot{\\varepsilon}_{II}`` are second invariants of
deviatoric stress and strain rate. Grain size is updated from the difference
between a grain-growth term and a grain-size-reduction term.

# Parameters
- `Name`: optional rheology name used by database and display utilities.
- `n`: stress exponent.
- `r`: water-fugacity exponent for diffusion creep.
- `p`: grain-size exponent. Diffusion creep laws commonly use a negative value.
- `A`: diffusion creep pre-exponential factor.
- `Ag`: grain-growth pre-exponential factor.
- `Qg`: activation energy for grain growth.
- `rg`: water-fugacity exponent for grain growth.
- `E`: activation energy for diffusion creep.
- `V`: activation volume for diffusion creep.
- `R`: universal gas constant.
- `λg`: fraction of work converted into grain-size reduction.
- `γd`: grain-boundary energy.
- `c`: geometrical factor in the grain-size-reduction law.
- `Apparatus`: experimental geometry. Use `AxialCompression`, `SimpleShear`, or
  `Invariant`. Choose `Invariant` when parameters are already expressed in
  second-invariant form.

# Example
```julia-repl
julia> x2 = DiffusionCreepGrainSize(Name="test")
DiffusionCreepGrainSize: Name = test, n=1.0, r=0.0, p=-3.0, A=1.5 m³·⁰ MPa⁻¹·⁰ s⁻¹·⁰, E=500.0 kJ mol⁻¹·⁰, V=2.4e-5 m³·⁰ mol⁻¹·⁰, FT=1.7320508075688772, FE=1.1547005383792517)
```
"""
struct DiffusionCreepGrainSize{T, U0, U1, U2, U3, U4, U5, U6, U7} <: AbstractCreepLaw{T}
    Name::Ptr{UInt8}
    d0::GeoUnit{T, U0}  # initial grain size
    n::GeoUnit{T, U1}  # powerlaw exponent
    r::GeoUnit{T, U1}  # exponent of water-fugacity
    p::GeoUnit{T, U1}  # grain size exponent
    A::GeoUnit{T, U2}  # material specific rheological parameter
    E::GeoUnit{T, U3}  # activation energy
    V::GeoUnit{T, U4}  # activation volume
    R::GeoUnit{T, U5}  # universal gas constant
    λg::GeoUnit{T, U1} # fraction of work going into grain size reduction
    γd::GeoUnit{T, U6} # [J/m^2] grain boundary energy
    c::GeoUnit{T, U1}  # geometrical constant
    Ag::GeoUnit{T, U7} # material specific rheological parameter for grain growth
    Qg::GeoUnit{T, U3} # activation energy for grain growth
    rg::GeoUnit{T, U1} # exponent of water-fugacity for grain growth
    Apparatus::Int8 # type of experimental apparatus, either AxialCompression, SimpleShear or Invariant
    FT::T # type of experimental apparatus, either AxialCompression, SimpleShear or Invariant
    FE::T # type of experimental apparatus, either AxialCompression, SimpleShear or Invariant

    function DiffusionCreepGrainSize(;
            Name = "",
            d0 = 5e-3m,
            n = 1NoUnits,
            r = 0NoUnits,
            p = -3NoUnits,
            rg = 1.38NoUnits,
            A = 1.5e6Pa^(-n - r) / s * m^(-p),
            Ag = 1.5e6Pa^(-rg) / s * m^(-p),
            Qg = 134.0e3J / mol,
            E = 500.0e3J / mol,
            V = 24.0e-6m^3 / mol,
            R = 8.3145J / mol / K,
            λg = 0.015NoUnits,
            γd = 1J / m^2,
            c = π*NoUnits,
            Apparatus = AxialCompression,
        )

        # Corrections from lab experiments
        FT, FE = CorrectionFactor(Apparatus)
        # Convert to GeoUnits
        nU = convert(GeoUnit, rat2float(n))
        d0U = convert(GeoUnit, d0)
        rU = convert(GeoUnit, r)
        pU = convert(GeoUnit, p)
        AU = convert(GeoUnit, A)
        EU = convert(GeoUnit, E)
        VU = convert(GeoUnit, V)
        RU = convert(GeoUnit, R)
        λgU = convert(GeoUnit, λg)
        γdU = convert(GeoUnit, γd)
        cU = convert(GeoUnit, c)
        AgU = convert(GeoUnit, Ag)
        QgU = convert(GeoUnit, Qg)
        rgU = convert(GeoUnit, rg)
        # Extract struct types
        T = typeof(rU).types[1]
        U0 = typeof(d0U).types[2]
        U1 = typeof(rU).types[2]
        U2 = typeof(AU).types[2]
        U3 = typeof(EU).types[2]
        U4 = typeof(VU).types[2]
        U5 = typeof(RU).types[2]
        U6 = typeof(γdU).types[2]
        U7 = typeof(AgU).types[2]
        name = pointer(ptr2string(Name))
        # Create struct
        return new{T, U0, U1, U2, U3, U4, U5, U6, U7}(
            name, d0U, nU, rU, pU, AU, EU, VU, RU, λgU, γdU, cU, AgU, QgU, rgU, Int8(Apparatus), FT, FE
        )
    end
end

function DiffusionCreepGrainSize(Name, d0, n, r, p, A, E, V, R, λg, γd, c, Ag, Qg, rg, Apparatus, FT, FE)
    return DiffusionCreepGrainSize(;
        Name = Name, d0 = d0, n = n, r = r, p = p, A = A, Ag = Ag, Qg = Qg, rg = rg, E = E, V = V, R = R, λg = λg, γd = γd, c = c, Apparatus = Apparatus
    )
end

"""
    s = remove_tensor_correction(s::DiffusionCreepGrainSize)

Return a copy of `s` with `Apparatus=Invariant`.

This removes the stress and strain-rate tensor correction factors, which is
useful when comparing the implemented creep law with experimental curves from
publications that report differential stress and axial or shear strain rate
instead of tensor invariants.
"""
function remove_tensor_correction(s::DiffusionCreepGrainSize)
    # name = String(collect(s.Name))
    return DiffusionCreepGrainSize(;
        Name = unsafe_string(s.Name), d0 = s.d0, n = s.n, r = s.r, p = s.p, A = s.A, Ag = s.Ag, Qg = s.Qg, rg = s.rg, E = s.E, V = s.V, λg = s.λg, γd = s.γd, c = s.c, Apparatus = Invariant
    )
end

function param_info(s::DiffusionCreepGrainSize)
    name = unsafe_string(s.Name)
    eq = L"\tau_{ij} = 2 \eta  \dot{\varepsilon}_{ij}"
    if name == ""
        return MaterialParamsInfo(; Equation = eq)
    end

    inf = diffusion_database_info(name)
    return MaterialParamsInfo(;
        Equation = eq, Comment = inf.Comment, BibTex_Reference = inf.BibTex_Reference
    )
end

# Calculation routines for linear viscous rheologies
"""
    compute_εII(a::DiffusionCreepGrainSize, τII; T=one(precision(a)),
        P=zero(precision(a)), f=one(precision(a)), d=NumValue(a.d),
        kwargs...)

Return the diffusion creep strain-rate invariant as a function of the stress
invariant ``\\tau_{II}``.

```math
    \\dot{\\varepsilon}_{II}
        = A (F_T\\tau_{II})^n d^p f_\\mathrm{H2O}^r
          \\exp\\left(-\\frac{E + PV}{RT}\\right) / F_E
```

For this grain-size-evolution law, `d` may be updated from the grain-growth and
grain-size-reduction rates before evaluating the creep law.
"""
@inline function compute_εII(
        a::DiffusionCreepGrainSize,
        τII;
        T = one(precision(a)),
        P = zero(precision(a)),
        f = one(precision(a)),
        d = NumValue(a.d0),
        dt = zero(precision(a)),
        kwargs...,
    )
    @unpack_val n, r, p, A, E, V, R = a
    FT, FE = a.FT, a.FE

    if !iszero(dt)
        d_red = compute_grain_size_reduction(a, τII; d = d, kwargs...)
        d_g = compute_growth_rate(a, τII; f = f, T = T, d = d, kwargs...)
        d = d + (d_g - d_red) * dt
    end

    ε = @pow A * (τII * FT)^n * f^r * d^p * exp(-(E + P * V) / (R * T)) / FE
    return ε
end

@inline function compute_εII(
        a::DiffusionCreepGrainSize, τII::Quantity; T = 1K, P = 0Pa, f = 1NoUnits, d = Value(a.d0), dt = 0s, args...
    )
    @unpack_units n, r, p, A, E, V, R = a
    FT, FE = a.FT, a.FE

    if !iszero(dt)
        d_red = compute_grain_size_reduction(a, τII; d = d, args...)
        d_g = compute_growth_rate(a, τII; f = f, T = T, d = d, args...)
        d = d + (d_g - d_red) * dt
    end

    ε = @pow A * (τII * FT)^n * f^r * d^p * exp(-(E + P * V) / (R * T)) / FE
    return ε
end

"""
    compute_τII(a::DiffusionCreepGrainSize, EpsII; T=one(precision(a)),
        P=zero(precision(a)), f=one(precision(a)), kwargs...)

Return the stress invariant corresponding to the diffusion creep strain-rate
invariant `EpsII`.
"""
@inline function compute_τII(
        a::DiffusionCreepGrainSize,
        EpsII;
        T = one(precision(a)),
        P = zero(precision(a)),
        f = one(precision(a)),
        kwargs...,
    )
    @unpack_val d0, n, r, p, A, E, V, R = a
    FT, FE = a.FT, a.FE

    n_inv = inv(n)

    τ = @pow A^-n_inv *
        (EpsII * FE)^n_inv *
        f^(-r * n_inv) *
        d0^(-p * n_inv) *
        exp((E + P * V) / (n * R * T)) / FT

    return τ
end

@inline function compute_τII(
        a::DiffusionCreepGrainSize, EpsII::Quantity; T = 1K, P = 0Pa, f = 1NoUnits, kwargs...
    )
    @unpack_units d0, n, r, p, A, E, V, R = a
    FT, FE = a.FT, a.FE

    n_inv = inv(n)

    τ = @pow A^(-n_inv) *
        EpsII *
        FE *
        f^(-r * n_inv) *
        d0^(-p * n_inv) *
        exp((E + P * V) / (n * R * T)) / FT

    return τ
end

"""
    compute_grain_size_reduction(a::DiffusionCreepGrainSize, τII; kwargs...)

Return the grain-size-reduction rate associated with deformational work.
"""
function compute_grain_size_reduction(a::DiffusionCreepGrainSize, τII; εII = one(τII), d = one(τII), args...)
    @unpack_val λg, c, γd = a
    d_reduction = τII * εII * λg * d^2 / (c * γd)
    return d_reduction
end

function compute_grain_size_reduction(a::DiffusionCreepGrainSize, τII::Quantity; εII = one(τII), d = 1m, args...)
    @unpack_units λg, c, γd = a
    d_reduction = τII * εII * λg * d^2 / (c * γd)
    return d_reduction
end

"""
    compute_growth_rate(a::DiffusionCreepGrainSize, τII; f=1e0, T=0e0,
        d=1e-3, kwargs...)

Return the grain-growth rate for fugacity `f`, temperature `T`, and grain size
`d`.
"""
function compute_growth_rate(a::DiffusionCreepGrainSize, τII; f=1e0, T=0e0, d=1e-3, args...)
    @unpack_val Ag, rg, p, Qg, R = a
    d_growth = Ag * f^rg * d^(1-p) * exp(-Qg / (R * T))
    return d_growth
end

function compute_growth_rate(a::DiffusionCreepGrainSize, τII::Quantity; f = 1NoUnits, T = 1K, d = 1m, args...)
    @unpack_units Ag, rg, p, Qg, R = a
    d_growth = Ag * f^rg * d^(1-p) * exp(-Qg / (R * T))
    return d_growth
end

# Print info
function show(io::IO, g::DiffusionCreepGrainSize)
    return print(
        io,
        "DiffusionCreepGrainSize: Name = $(unsafe_string(g.Name)), n=$(Value(g.n)), r=$(Value(g.r)), p=$(Value(g.p)), A=$(Value(g.A)), E=$(Value(g.E)), V=$(Value(g.V)), FT=$(g.FT), FE=$(g.FE)",
    )
end
