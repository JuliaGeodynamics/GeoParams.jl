"""
    module GeoParams

Typical geodynamic simulations involve a large number of material parameters that have units that are often inconvenient to be directly used in numerical models
This package has two main features that help with this:
- Create a nondimensionalization object, which can be used to transfer dimensional to non-dimensional parameters (usually better for numerical solvers)
- Create an object in which you can specify material parameters employed in the geodynamic simulations

The material parameter object is designed to be extensible and can be passed on to the solvers, such that new creep laws or features can be readily added.
We also implement some typically used creep law parameters, together with tools to plot them versus and compare our results with those of published papers (to minimize mistakes).
"""
module GeoParams

    using Parameters         # helps setting default parameters in structures
    using Unitful            # Units
    using BibTeX             # references of creep laws
    using StaticArrays
    using LinearAlgebra
    using ForwardDiff
    using MuladdMacro

    import Base: getindex

    # overload to account for cases where this is an integer
    for T in (:Real, :Symbol)
        @eval begin
            Base.getindex(val::$(T), I::Vararg{Integer, N}) where {N} = val
            Base.getindex(val::$(T), I::Integer) = val
        end
    end

    export @u_str,
        uconvert,
        unit,
        ustrip,
        NoUnits, #  Units
        GeoUnit,
        GeoUnits,
        GEO_units,
        SI_units,
        NO_units,
        AbstractGeoUnit,
        nondimensionalize,
        dimensionalize,
        dimensionalize_and_strip,
        @dimstrip,
        superscript,
        upreferred,
        GEO,
        SI,
        NONE,
        isDimensional,
        Value,
        NumValue,
        unpack_units,
        Unit,
        UnitValue,
        isdimensional,
        km,
        m,
        cm,
        mm,
        μm,
        Myr,
        yr,
        s,
        GPa,
        MPa,
        Pa,
        bar,
        kbar,
        Pas,
        K,
        C,
        g,
        kg,
        mol,
        J,
        kJ,
        Watt,
        μW,
        Quantity

    #
    """
        AbstractMaterialParam

    Supertype of every individual material property parameterization (density, elasticity,
    viscosity, conductivity, …). Concrete subtypes are callable and evaluate the property they
    describe.
    """
    abstract type AbstractMaterialParam end

    """
        AbstractMaterialParamsStruct

    Supertype of the per-phase container that bundles all `AbstractMaterialParam`s belonging to a
    single material phase (see `MaterialParams`).
    """
    abstract type AbstractMaterialParamsStruct end

    """
        AbstractPhaseDiagramsStruct <: AbstractMaterialParam

    Supertype of material parameters obtained by interpolating a precomputed phase diagram lookup
    table as a function of pressure and temperature.
    """
    abstract type AbstractPhaseDiagramsStruct <: AbstractMaterialParam end

    """
        AbstractConstitutiveLaw{T} <: AbstractMaterialParam

    Supertype of the constitutive laws that relate stress and strain rate (creep laws, elasticity,
    plasticity). `T` is the numeric element type.
    """
    abstract type AbstractConstitutiveLaw{T} <: AbstractMaterialParam end

    """
        AbstractComposite <: AbstractMaterialParam

    Supertype of composite rheologies that combine several constitutive laws (see
    [`CompositeRheology`](@ref) and [`Parallel`](@ref)).
    """
    abstract type AbstractComposite <: AbstractMaterialParam end

    function PerpleX_LaMEM_Diagram end                                         # necessary as we already use this function in Units, but only define it later in PhaseDiagrams

    """
        param_info(s::AbstractMaterialParam) -> MaterialParamsInfo

    Returns a `MaterialParamsInfo` describing the parameterization `s`: its governing equation (as a
    `LaTeXString`) and, where available, a comment and BibTeX reference. Each concrete
    material-parameter type provides its own method.
    """
    function param_info end
    export AbstractMaterialParam, AbstractMaterialParamsStruct, AbstractPhaseDiagramsStruct

    include("Utils.jl")
    export value_and_partial

    include("TensorAlgebra/TensorAlgebra.jl")
    export second_invariant, second_invariant_staggered, rotate_elastic_stress

    # note that this throws a "Method definition warning regarding superscript"; that is expected & safe
    #  as we add a nicer way to create output of superscripts. I have been unable to get rid of this warning,
    #  as I am indeed redefining a method originally defined in Unitful
    include("Units.jl")
    using .Units
    export @unpack_units, @unpack_val
    export compute_units, udim

    include("Interpolations.jl")
    using .Interpolations
    export LinearInterpolator, interpolate_field

    # Define Material Parameter structure
    include("MaterialParameters.jl")
    using .MaterialParameters
    export MaterialParams, SetMaterialParams, No_MaterialParam, MaterialParamsInfo

    # Phase Diagrams
    using .MaterialParameters.PhaseDiagrams
    export PhaseDiagram_LookupTable, PerpleX_LaMEM_Diagram, MAGEMin_Diagram, MAGEMin_LookupTable

    # Density
    using .MaterialParameters.Density
    export compute_density, # computational routines
        compute_density!,
        param_info,
        AbstractDensity,
        ConduitDensity,
        ConstantDensity,
        PT_Density,
        Compressible_Density,
        T_Density,
        Vector_Density,
        PhaseDiagram_LookupTable,
        MeltDependent_Density,
        BubbleFlow_Density,
        GasPyroclast_Density,
        RedlichKwong_Density,
        IdealGas_Density,
        ThreePhase_Density,
        Melt_DensityX,
        compute_density_ratio

    # Constitutive relationships laws
    using .MaterialParameters.ConstitutiveRelationships
    export AxialCompression, SimpleShear, Invariant

    #       Calculation routines
    export dεII_dτII,
        dτII_dεII,
        dεII_dτII_AD,
        dτII_dεII_AD,
        dεvol_dp,
        dp_dεvol,
        compute_εII!,
        compute_εII,
        compute_εII_AD,
        compute_τII!,
        compute_τII,
        compute_τII_AD,
        compute_εvol!,
        compute_εvol,
        compute_p!,
        compute_p,
        CorrectionFactor,
        remove_tensor_correction,
        isvolumetric,

        #       Viscous creep laws
        AbstractCreepLaw,
        LinearViscous,
        LinearMeltViscosity,
        ViscosityPartialMelt_Costa_etal_2009,
        GiordanoMeltViscosity,
        PowerlawViscous,
        ArrheniusType,
        HerschelBulkley,
        CustomRheology,
        DislocationCreep,
        SetDislocationCreep,
        DiffusionCreep,
        SetDiffusionCreep,
        GrainBoundarySliding,
        SetGrainBoundarySliding,
        PeierlsCreep,
        SetPeierlsCreep,
        NonLinearPeierlsCreep,
        SetNonLinearPeierlsCreep,
        Transform_DislocationCreep,
        Transform_DiffusionCreep,
        Transform_GrainBoundarySliding,
        Transform_PeierlsCreep,
        Transform_NonLinearPeierlsCreep,
        Peierls_stress_iterations,

        #       Elasticity
        AbstractElasticity,
        ConstantElasticity,
        SetConstantElasticity,
        effective_εII,
        iselastic,
        get_shearmodulus,
        get_bulkmodulus,

        #       softening
        AbstractSoftening,
        NoSoftening,
        LinearSoftening,
        NonLinearSoftening,
        DecaySoftening,

        #       Plasticity
        AbstractPlasticity,
        compute_yieldfunction,
        compute_yieldfunction!,
        compute_flowpotential,
        compute_flowpotential!,
        DruckerPrager,
        DruckerPrager_regularised,
        DruckerPragerCap,
        compute_plasticpotentialDerivative,
        ∂Q∂τ,
        ∂Q∂P, ∂Q∂τII,
        ∂F∂τII, ∂F∂P, ∂F∂λ,

        #       Composite rheologies
        AbstractConstitutiveLaw,
        AbstractComposite,
        computeViscosity_εII,
        computeViscosity_εII_AD,
        local_iterations_εII,
        local_iterations_εII_AD,
        local_iterations_τII,
        local_iterations_τII_AD,
        InverseCreepLaw,
        CompositeRheology,
        Parallel,
        create_rheology_string, print_rheology_matrix,
        compute_εII_harmonic, compute_τII_AD,
        isplastic, isvolumetricplastic,
        compute_p_τII,
        local_iterations_εvol,
        compute_p_harmonic

    include("CreepLaw/Data_deprecated/DislocationCreep.jl")
    include("CreepLaw/Data_deprecated/DiffusionCreep.jl")
    include("CreepLaw/Data_deprecated/GrainBoundarySliding.jl")
    include("CreepLaw/Data_deprecated/NonLinearPeierlsCreep.jl")
    include("CreepLaw/Data_deprecated/PeierlsCreep.jl")
    export DislocationCreep_info,
        DiffusionCreep_info,
        GrainBoundarySliding_info,
        PeierlsCreep_info,
        NonLinearPeierlsCreep_info

    # Constitutive relationships laws
    include("StressComputations/StressComputations.jl")
    export compute_τij, compute_p_τij, compute_τij_stagcenter!, compute_p_τij_stagcenter!, compute_τij!, compute_p_τij!

    include("Rheology_Utils.jl")
    export time_τII_0D, time_τII_0D!, time_p_τII_0D, time_p_τII_0D!

    include("Viscosity/Viscosity.jl")
    export compute_viscosity_εII,
        compute_viscosity_τII,
        compute_elastoviscosity,
        compute_elastoviscosity_εII,
        compute_elastoviscosity_τII,
        compute_viscosity,
        compute_elasticviscosity

    # Gravitational Acceleration
    using .MaterialParameters.GravitationalAcceleration
    export compute_gravity, # computational routines
        ConstantGravity,
        DippingGravity


    using .MaterialParameters.ChemicalDiffusion
    export AbstractChemicalDiffusion,
        DiffusionData,
        MeltMulticompDiffusionData,
        compute_D,
        compute_D!,
        compute_λ,
        compute_λ!,
        SetChemicalDiffusion,
        SetMulticompChemicalDiffusion,
        Transform_ChemicalDiffusion


    export Rutile,
        Garnet,
        Olivine,
        Melt

    # Energy parameters: Heat Capacity, Thermal conductivity, latent heat, radioactive heat
    using .MaterialParameters.HeatCapacity
    export compute_heatcapacity,
        compute_heatcapacity!, ConstantHeatCapacity, T_HeatCapacity_Whittington, Latent_HeatCapacity, Vector_HeatCapacity

    using .MaterialParameters.Conductivity
    export compute_conductivity,
        compute_conductivity!,
        ConstantConductivity,
        T_Conductivity_Whittington,
        T_Conductivity_Whittington_parameterised,
        TP_Conductivity,
        Set_TP_Conductivity

    using .MaterialParameters.LatentHeat
    export compute_latent_heat, compute_latent_heat!, ConstantLatentHeat

    using .MaterialParameters.RadioactiveHeat
    export compute_radioactive_heat,
        compute_radioactive_heat!, ConstantRadioactiveHeat, ExpDepthDependentRadioactiveHeat

    using .MaterialParameters.Shearheating
    export compute_shearheating!, compute_shearheating, ConstantShearheating

    # Add TAS classification
    include("./RockClassification/TASclassification.jl")
    using .TASclassification
    export TASclassificationData, computeTASclassification, retrieveTASrockType

    # Add zircon saturation parameterizations
    include("./ZirconAge/ZirconAges.jl")
    using .ZirconAges
    export ZirconAgeData,
        compute_zircon_age_PDF,
        compute_zircons_Ttpath,
        zircon_age_PDF,
        compute_zircons_convert_vecs2mat

    # Seismic velocities
    using .MaterialParameters.SeismicVelocity
    export compute_wave_velocity,
        compute_wave_velocity!,
        ConstantSeismicVelocity,
        anelastic_correction,
        melt_correction,
        porosity_correction,
        correct_wavevelocities_phasediagrams,
        melt_correction_Takei

    # Add melting parameterizations
    include("./MeltFraction/MeltingParameterization.jl")
    using .MeltingParam
    export compute_meltfraction,
        compute_meltfraction!, # calculation routines
        compute_meltfraction_ratio,
        compute_dϕdT,
        compute_dϕdT!,
        MeltingParam_Caricchi,
        MeltingParam_Smooth3rdOrder,
        MeltingParam_4thOrder,
        MeltingParam_5thOrder,
        MeltingParam_Quadratic,
        MeltingParam_Assimilation,
        MeltingParam_Volatile,
        MeltingParam_MaficVolatile,
        Vector_MeltingParam,
        SmoothMelting


    using .MaterialParameters.Permeability
    export compute_permeability,
        compute_permeability!,
        compute_permeability_ratio,
        param_info,
        AbstractPermeability,
        ConstantPermeability,
        HazenPermeability,
        PowerLawPermeability,
        CarmanKozenyPermeability

    using .MaterialParameters.Solubility
    export compute_dissolved,
        compute_dissolved!,
        compute_dissolved_ratio,
        ∂dissolved_∂P,
        ∂dissolved_∂T,
        ∂dissolved_∂Xco2,
        find_Xco2,
        AbstractSolubility,
        Liu2005_Solubility,
        Mafic_Solubility,
        GasMixture,
        compute_gas_heatcapacity,
        effective_molar_mass

    include("Traits/rheology.jl")
    export RheologyTrait
    export islinear, LinearRheologyTrait, NonLinearRheologyTrait
    export isviscoelastic, ElasticRheologyTrait, NonElasticRheologyTrait
    export isplasticity, PlasticRheologyTrait, NonPlasticRheologyTrait

    include("Traits/density.jl")
    export isconstant, DensityTrait, ConstantDensityTrait, NonConstantDensityTrait

    include("CreepLaw/Data/DislocationCreep.jl")
    using .Dislocation

    include("CreepLaw/Data/DiffusionCreep.jl")
    using .Diffusion

    include("CreepLaw/Data/GrainBoundarySliding.jl")
    using .GBS

    include("CreepLaw/Data/NonLinearPeierlsCreep.jl")
    using .NonLinearPeierls

    include("CreepLaw/Data/PeierlsCreep.jl")
    using .Peierls

    function creeplaw_list(m::Module)
        out = string.(names(m; all = true, imported = true))
        filter!(x -> !startswith(x, "#"), out)
        return [getfield(m, Symbol(x)) for x in out if !isnothing(tryparse(Int, string(x[end]))) || endswith(x, "a") || endswith(x, "b")]
    end

    "Returns the list of pre-defined diffusion creep laws (entries of the `Diffusion` submodule)."
    diffusion_law_list() = creeplaw_list(Diffusion)
    "Returns the list of pre-defined dislocation creep laws (entries of the `Dislocation` submodule)."
    dislocation_law_list() = creeplaw_list(Dislocation)
    "Returns the list of pre-defined grain-boundary-sliding creep laws (entries of the `GBS` submodule)."
    grainboundarysliding_law_list() = creeplaw_list(GBS)
    "Returns the list of pre-defined non-linear Peierls creep laws (entries of the `NonLinearPeierls` submodule)."
    nonlinearpeierls_law_list() = creeplaw_list(NonLinearPeierls)
    "Returns the list of pre-defined Peierls creep laws (entries of the `Peierls` submodule)."
    peierls_law_list() = creeplaw_list(Peierls)

    export diffusion_law_list,
        dislocation_law_list,
        grainboundarysliding_law_list,
        nonlinearpeierls_law_list,
        peierls_law_list


    # Define Table output functions
    include("Tables.jl")
    using .Tables
    export detachFloatfromExponent, extract_parameters_from_phases, Dict2LatexTable, extract_parameters_from_phases_md, Dict2MarkdownTable, ParameterTable

    # Add 1D Strength Envelope
    include("./StrengthEnvelope/StrengthEnvelope.jl")

    # Add plotting routines - only activated if the "GLMakie.jl" package is loaded
    #
    # Add function definitions here such that they can be exported from GeoParams.jl
    # and extended in the GeoParamsMakieExt package extension or by the
    # GLMakie-specific code loaded by Requires.jl
    # Each plotting routine is a stub here and implemented in the `GeoParamsMakieExt`
    # package extension; a `Makie.jl` backend (e.g. `GLMakie.jl`, or `CairoMakie.jl`
    # for headless use) must be loaded for the methods to become available.

    """
        PlotStrainrateStress(x; kwargs...)

    Plots deviatoric stress versus deviatoric strain rate for one or more creep laws `x`.
    """
    function PlotStrainrateStress end

    """
        PlotStressStrainrate(x; kwargs...)

    Plots deviatoric strain rate versus deviatoric stress for one or more creep laws `x`
    (the transpose of [`PlotStrainrateStress`](@ref)).
    """
    function PlotStressStrainrate end

    """
        PlotStrainrateViscosity(x; kwargs...)

    Plots effective viscosity versus deviatoric strain rate for one or more creep laws `x`.
    """
    function PlotStrainrateViscosity end

    """
        PlotStressViscosity(x; kwargs...)

    Plots effective viscosity versus deviatoric stress for one or more creep laws `x`.
    """
    function PlotStressViscosity end

    """
        PlotHeatCapacity(Cp::AbstractHeatCapacity; kwargs...)

    Plots heat capacity as a function of temperature for the parameterization `Cp`.
    """
    function PlotHeatCapacity end

    """
        PlotConductivity(k::AbstractConductivity; kwargs...)

    Plots thermal conductivity as a function of temperature for the parameterization `k`.
    """
    function PlotConductivity end

    """
        PlotMeltFraction(p::AbstractMeltingParam; kwargs...)

    Plots melt fraction and `dϕ/dT` as a function of temperature for the parameterization `p`.
    """
    function PlotMeltFraction end

    """
        PlotPhaseDiagram(p::AbstractPhaseDiagramsStruct, fieldname::Symbol; kwargs...)

    Plots the field `fieldname` of a phase diagram as a function of temperature (x-axis) and pressure (y-axis).
    """
    function PlotPhaseDiagram end

    """
        Plot_TAS_diagram(point; kwargs...)

    Plots a TAS (total-alkali versus silica) classification diagram for the given composition `point`.
    """
    function Plot_TAS_diagram end

    """
        Plot_ZirconAge_PDF(time_Ma, PDF_zircons, time_Ma_average, PDF_zircon_average)

    Plots the zircon-age probability density function computed from a simulation.
    """
    function Plot_ZirconAge_PDF end

    """
        PlotDeformationMap(v; kwargs...)

    Plots a deformation-mechanism map (deformation regime as a function of temperature and stress or strain rate) for the rheology `v`.
    """
    function PlotDeformationMap end

    """
        PlotStressTime_0D(x; εII, kwargs...)

    Plots the stress evolution over time of a 0-D visco-elasto-(plastic) model for the rheology `x`.
    """
    function PlotStressTime_0D end

    """
        PlotPressureStressTime_0D(x; εII, εvol, kwargs...)

    Plots the pressure and stress evolution over time of a 0-D visco-elasto-(plastic) model for the rheology `x`.
    """
    function PlotPressureStressTime_0D end

    function StrengthEnvelopePlot end
    function PlotDiffusionCoef end

    """
        PlotDiffusionCoefArrhenius(x; kwargs...)

    Arrhenius plot of the diffusion coefficient (`log(D)` versus `10⁴/T`) for one or more `ChemicalDiffusionData` structures `x`.
    """
    function PlotDiffusionCoefArrhenius end

    export PlotStrainrateStress,
        PlotStressStrainrate,
        PlotStrainrateViscosity,
        PlotStressViscosity,
        PlotHeatCapacity,
        PlotConductivity,
        PlotMeltFraction,
        PlotPhaseDiagram,
        Plot_TAS_diagram,
        Plot_ZirconAge_PDF,
        PlotDeformationMap,
        PlotStressTime_0D,
        PlotPressureStressTime_0D,
        StrengthEnvelopePlot,
        PlotDiffusionCoefArrhenius

    #Set functions aliases using @use
    include("aliases.jl")
    export ntuple_idx

    for modulus in (:G, :Kb)
        fun = Symbol("get_$(string(modulus))")
        @eval begin
            @inline $(fun)(a::ConstantElasticity) = a.$(modulus).val
            @inline $(fun)(c::CompositeRheology) = $(fun)(isviscoelastic(c), c)
            @inline $(fun)(::ElasticRheologyTrait, c::CompositeRheology) = mapreduce(x -> $(fun)(x), +, c.elements)
            @inline $(fun)(r::AbstractMaterialParamsStruct) = $(fun)(r.CompositeRheology[1])
            @inline $(fun)(a::NTuple{N, AbstractMaterialParamsStruct}, phase) where {N} = nphase($(fun), phase, a)
            @inline $(fun)(::NonElasticRheologyTrait, c::CompositeRheology) = 0
            @inline $(fun)(::Union{NonElasticRheologyTrait, AbstractCreepLaw, AbstractPlasticity, AbstractConstitutiveLaw}) = 0
        end
    end

    """
        get_G(r)

    Returns the elastic shear modulus `G` of the rheology `r` (a `ConstantElasticity`,
    [`CompositeRheology`](@ref), or `MaterialParams`), summing the contributions of elastic elements
    and returning `0` when `r` is non-elastic. Also available as `get_shearmodulus`.
    """
    get_G

    """
        get_Kb(r)

    Returns the elastic bulk modulus `Kb` of the rheology `r` (a `ConstantElasticity`,
    [`CompositeRheology`](@ref), or `MaterialParams`), summing the contributions of elastic elements
    and returning `0` when `r` is non-elastic. Also available as `get_bulkmodulus`.
    """
    get_Kb

    export get_G, get_Kb

    const get_shearmodulus = get_G
    const get_bulkmodulus = get_Kb

end # module GeoParams
