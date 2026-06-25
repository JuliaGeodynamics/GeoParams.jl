using Test
using GeoParams
using CairoMakie
import GeoParams.Dislocation, GeoParams.Diffusion, GeoParams.Garnet

CairoMakie.activate!()

# All figures are written to pngs in a temp dir (forces the `save` path, never `display`),
# which lets these run on a headless CI runner without OpenGL/xvfb.
@testset "Plotting (headless / CairoMakie)" begin
    tmp = mktempdir()
    npng = Ref(0)
    png() = (npng[] += 1; joinpath(tmp, "fig_$(npng[]).png"))

    # rheologies used throughout
    v_disl = SetDislocationCreep(Dislocation.dry_olivine_Hirth_2003)
    v_diff = SetDiffusionCreep(Diffusion.dry_anorthite_Rybacki_2006)
    comp = CompositeRheology(v_diff, v_disl)                 # multi-component (deformation map)
    ve = CompositeRheology(LinearViscous(; η = 1.0e21Pa * s), ConstantElasticity())  # visco-elastic (time evolution)
    args = (T = 1000.0, P = 0.0, d = 1.0e-3, f = 1.0)

    @testset "stress / strainrate / viscosity" begin
        @test PlotStrainrateStress(v_disl; args = args, Strainrate = (1.0e-18, 1.0e-12), filename = png())[1] isa Figure
        @test PlotStressStrainrate(v_disl; args = args, Stress = (1.0e0, 1.0e8), filename = png())[1] isa Figure
        @test PlotStrainrateViscosity(v_disl; args = args, Strainrate = (1.0e-18, 1.0e-12), filename = png())[1] isa Figure
        @test PlotStressViscosity(v_disl; args = args, Stress = (1.0e0, 1.0e8), filename = png())[1] isa Figure
        # tuple input + per-curve args (covers the tuple/per-curve branches of each)
        vt = (v_disl, v_diff)
        @test PlotStrainrateStress(vt; args = args, filename = png())[1] isa Figure
        @test PlotStressStrainrate(vt; args = args, filename = png())[1] isa Figure
        @test PlotStrainrateViscosity(vt; args = args, filename = png())[1] isa Figure
        @test PlotStressViscosity(vt; args = args, filename = png())[1] isa Figure
        # no `filename` -> exercises the `display(fig)` branch (a no-op headless with CairoMakie)
        @test PlotStrainrateStress(v_disl; args = args)[1] isa Figure
    end

    @testset "energy & melt fraction" begin
        @test PlotConductivity(T_Conductivity_Whittington(); filename = png())[1] isa Figure
        @test PlotConductivity(ConstantConductivity(); filename = png())[1] isa Figure
        @test PlotHeatCapacity(T_HeatCapacity_Whittington(); filename = png())[1] isa Figure
        T, phi, dϕdT = PlotMeltFraction(MeltingParam_4thOrder(); filename = png())
        @test length(phi) == length(T) == length(dϕdT)
    end

    @testset "phase diagrams (Perple_X + MAGEMin)" begin
        PD = PerpleX_LaMEM_Diagram(joinpath(@__DIR__, "test_data", "Peridotite.in"))
        @test PlotPhaseDiagram(PD, :Rho; filename = png())[1] isa Figure
        @test PlotPhaseDiagram(PD, :meltFrac; filename = png())[1] isa Figure
        PDm = MAGEMin_Diagram(joinpath(@__DIR__, "test_data", "MAGEMin_Rhyolite.in"))
        @test PlotPhaseDiagram(PDm, :Rho; filename = png())[1] isa Figure
        @test PlotPhaseDiagram(PDm, :meltFrac; filename = png())[1] isa Figure
    end

    @testset "0D time evolution & deformation map" begin
        @test PlotStressTime_0D(ve; args = args, εII = 1.0e-15, Time = (1.0e0, 1.0e10), nt = 10, filename = png())[1] isa Figure
        @test PlotPressureStressTime_0D(ve; args = args, εII = 1.0e-15, εvol = -1.0e-18, Time = (1.0e0, 1.0e10), nt = 10, filename = png())[1] isa Figure
        @test PlotDeformationMap(comp; n = 20, filename = png()) isa Figure
        # option branches (degenerate-data contours now guarded, so these no longer error)
        @test PlotDeformationMap(comp; n = 20, strainrate = false, filename = png()) isa Figure
        @test PlotDeformationMap(comp; n = 20, viscosity = true, filename = png()) isa Figure
        @test PlotDeformationMap(comp; n = 20, grainsize = true, filename = png()) isa Figure
        @test PlotDeformationMap(comp; n = 20, rotate_axes = true, filename = png()) isa Figure
        @test PlotDeformationMap(comp; n = 20, annotate = true, filename = png()) isa Figure
        @test PlotDeformationMap(comp; n = 20, log_stress = false, filename = png()) isa Figure
    end

    @testset "TAS / Arrhenius / zircon" begin
        @test Plot_TAS_diagram([50.0 5.0; 61.6 6.25]) isa Figure          # uses display() -> no-op headless
        Fe = SetChemicalDiffusion(Garnet.Grt_Fe_Chakraborty1992)
        @test PlotDiffusionCoefArrhenius(Fe; filename = png())[1] isa Figure
        # synthetic zircon-age PDF data
        tM = [collect(0.0:1.0e5:1.0e6) for _ in 1:3]
        pdf = [exp.(-(tM[1] .- 5.0e5) .^ 2 ./ 1.0e11) for _ in 1:3]
        @test Plot_ZirconAge_PDF(tM, pdf, tM[1], pdf[1]) isa Figure
    end

    @testset "StrengthEnvelopePlot (all 3 temperature structures)" begin
        SE_MatParam = (
            SetMaterialParams(;
                Name = "UC", Phase = 1, Density = ConstantDensity(; ρ = 2700kg / m^3),
                CreepLaws = SetDislocationCreep(Dislocation.wet_quartzite_Ueda_2008),
                Plasticity = DruckerPrager(; ϕ = 30.0, C = 10MPa),
            ),
            SetMaterialParams(;
                Name = "LC", Phase = 2, Density = ConstantDensity(; ρ = 2900kg / m^3),
                CreepLaws = SetDislocationCreep(Dislocation.plagioclase_An75_Ji_1993),
                Plasticity = DruckerPrager(; ϕ = 20.0, C = 10MPa),
            ),
        )
        SE_thick = [20, 20] * km
        @test StrengthEnvelopePlot(SE_MatParam, SE_thick; TempType = LinTemp(), nz = 51) isa Figure
        @test StrengthEnvelopePlot(SE_MatParam, SE_thick; TempType = HalfspaceCoolTemp(0C, 1350C, 10Myr, 0K / km, 1.0e-6m^2 / s), nz = 51) isa Figure
        @test StrengthEnvelopePlot(SE_MatParam, SE_thick; TempType = ConstTemp(), nz = 51) isa Figure
    end

    # every figure that used `filename` was actually written
    @test all(i -> isfile(joinpath(tmp, "fig_$(i).png")), 1:npng[])
end
