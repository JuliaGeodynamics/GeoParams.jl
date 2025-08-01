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
        upreffered,
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
        superscript,
        upreferred,
        GEO,
        SI,
        NONE,
        isDimensional,
        Value,
        NumValue,
        unpack_units,
        unpack_val,
        Unit,
        UnitValue,
        isdimensional,
        km,
        m,
        cm,
        mm,
        μm,
        Myrs, # to remove at some point
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

    export AbstractGeoUnit1, GeoUnit1

    #
    abstract type AbstractMaterialParam end                                    # structure that holds material parameters (density, elasticity, viscosity)
    abstract type AbstractMaterialParamsStruct end                             # will hold all info for a phase
    abstract type AbstractPhaseDiagramsStruct <: AbstractMaterialParam end    # will hold all info for phase diagrams
    abstract type AbstractConstitutiveLaw{T} <: AbstractMaterialParam end
    abstract type AbstractComposite <: AbstractMaterialParam end

    function PerpleX_LaMEM_Diagram end                                         # necessary as we already use this function in Units, but only define it later in PhaseDiagrams
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

    # Define Material Parameter structure
    include("MaterialParameters.jl")
    using .MaterialParameters
    export MaterialParams, SetMaterialParams, No_MaterialParam, MaterialParamsInfo

    # Phase Diagrams
    using .MaterialParameters.PhaseDiagrams
    export PhaseDiagram_LookupTable, PerpleX_LaMEM_Diagram

    # Density
    using .MaterialParameters.Density
    export compute_density, # computational routines
        compute_density!,
        param_info,
        AbstractDensity,
        ConduitDensity,
        No_Density,
        ConstantDensity,
        PT_Density,
        Compressible_Density,
        T_Density,
        Vector_Density,
        PhaseDiagram_LookupTable,
        Read_LaMEM_Perple_X_Diagram,
        MeltDependent_Density,
        BubbleFlow_Density,
        GasPyroclast_Density,
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
        DislocationCreep_data,
        DiffusionCreep_data,
        GrainBoundarySliding_data,
        PeierlsCreep_data,
        NonLinearPeierlsCreep_data,
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
        DruckerPrager,
        DruckerPrager_regularised,
        compute_plasticpotentialDerivative,
        ∂Q∂τ,
        ∂Q∂P, ∂Q∂τII,
        ∂F∂τII, ∂F∂P, ∂F∂λ,

        #       Composite rheologies
        AbstractConstitutiveLaw,
        AbstractComposite,
        computeViscosity_εII,
        computeViscosity_εII!,
        computeViscosity_εII_AD,
        local_iterations_εII,
        local_iterations_εII_AD,
        local_iterations_τII,
        local_iterations_τII_AD,
        computeViscosity,
        InverseCreepLaw,
        KelvinVoigt,
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
        # melt_correction,
        # porosity_correction,
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

    diffusion_law_list() = creeplaw_list(Diffusion)
    dislocation_law_list() = creeplaw_list(Dislocation)
    grainboundarysliding_law_list() = creeplaw_list(GBS)
    nonlinearpeierls_law_list() = creeplaw_list(NonLinearPeierls)
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
    function PlotStrainrateStress end
    function PlotStressStrainrate end
    function PlotStrainrateViscosity end
    function PlotStressViscosity end
    function PlotHeatCapacity end
    function PlotConductivity end
    function PlotMeltFraction end
    function PlotPhaseDiagram end
    function Plot_TAS_diagram end
    function Plot_ZirconAge_PDF end
    function PlotDeformationMap end
    function PlotStressTime_0D end
    function PlotPressureStressTime_0D end
    function StrengthEnvelopePlot end
    function PlotDiffusionCoef end
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

    export get_G, get_Kb

    const get_shearmodulus = get_G
    const get_bulkmodulus = get_Kb

end # module GeoParams
