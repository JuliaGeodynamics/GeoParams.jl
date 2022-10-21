using Test, GeoParams, StaticArrays

@testset "Density.jl" begin

    # J2 TENSOR INVARIANT: 2D TESTS 
    τ_xx, τ_yy, τ_xy = 1.0, 2.0, 3.0
    τII = √( 0.5*(τ_xx^2 + τ_yy^2) + τ_xy^2)

    # test standard 2D interfaces
    @test τII == second_invariant(τ_xx, τ_yy, τ_xy)
    @test τII == second_invariant((τ_xx, τ_yy, τ_xy))
    @test τII == second_invariant(@SVector([τ_xx, τ_yy, τ_xy]))

    # test for staggered grids
    τ_xy_ij = τ_xy_11, τ_xy_12, τ_xy_21, τ_xy_22 = 3.0, 3.0, 3.5, 3.5 # shear stress at the center of the cell around the grid nodes
    τ_xy2_av = sum(ij^2 for ij in τ_xy_ij) * 0.25
    τII = √( 0.5*(τ_xx^2 + τ_yy^2) + τ_xy2_av)
    @test τII == second_invariant_staggered(τ_xx, τ_yy, τ_xy_ij)


    # J2 TENSOR INVARIANT: 3D TESTS 
    τ_xx, τ_yy, τ_zz, τ_yz, τ_xz, τ_xy = 1.0, 2.0, 3.0, 4.0, 5.0, 6.0
    τII = √( 0.5*(τ_xx^2 + τ_yy^2 + τ_zz^2) + τ_yz^2 + τ_xz^2 + τ_xy^2)

    # test standard 3D interfaces
    @test τII == second_invariant(τ_xx, τ_yy, τ_zz, τ_yz, τ_xz, τ_xy)
    @test τII == second_invariant((τ_xx, τ_yy, τ_zz, τ_yz, τ_xz, τ_xy))
    @test τII == second_invariant(@SVector([τ_xx, τ_yy, τ_zz, τ_yz, τ_xz, τ_xy]))

    # test for staggered grids
    # shear stress at the center of the cell around the grid nodes
    τ_yz_ij = τ_yz_11, τ_yz_12, τ_yz_21, τ_yz_22 = 4.0, 4.0, 4.5, 4.5 
    τ_xz_ij = τ_xz_11, τ_xz_12, τ_xz_21, τ_xz_22 = 5.0, 5.0, 5.5, 5.5
    τ_zy_ij = τ_zy_11, τ_zy_12, τ_zy_21, τ_zy_22 = 6.0, 6.0, 6.5, 6.5
    τ_yz2_av = sum(ij^2 for ij in τ_yz_ij) * 0.25
    τ_xz2_av = sum(ij^2 for ij in τ_xz_ij) * 0.25
    τ_xy2_av = sum(ij^2 for ij in τ_xy_ij) * 0.25
    τII = √( 0.5*(τ_xx^2 + τ_yy^2 + τ_zz^2) + τ_yz2_av + τ_xz2_av + τ_xy2_av)

    @test τII == second_invariant_staggered(τ_xx, τ_yy, τ_zz, τ_yz_ij, τ_xz_ij, τ_xy_ij)
    
    # normal stress at the center of the cell around the grid nodes
    τ_xx_ij = τ_xx_11, τ_xx_12, τ_yz_21, τ_xx_22 = 1.0, 1.0, 1.5, 1.5 
    τ_yy_ij = τ_yy_11, τ_yy_12, τ_xz_21, τ_yy_22 = 2.0, 2.0, 2.5, 2.5
    τ_zz_ij = τ_zz_11, τ_zz_12, τ_zy_21, τ_zz_22 = 3.0, 3.0, 3.5, 3.5
    τ_xx2_av = sum(ij^2 for ij in τ_xx_ij) * 0.25
    τ_yy2_av = sum(ij^2 for ij in τ_yy_ij) * 0.25
    τ_zz2_av = sum(ij^2 for ij in τ_zz_ij) * 0.25
    τII = √( 0.5*(τ_xx2_av + τ_yy2_av + τ_zz2_av) + τ_yz^2 + τ_xz^2 + τ_xy^2)

    @test τII == second_invariant_staggered(τ_xx_ij, τ_yy_ij, τ_zz_ij, τ_yz, τ_xz, τ_xy)
end
