using Aqua, Test, GeoParams

@testset "Aqua" begin
    # unbound_args: several tuple methods constrain a shared element type that `N == 0` leaves unbound.
    # piracies: `getindex` on `Real` and `Symbol` lets scalars stand in for arrays in the compute routines.
    # undocumented_names: the re-exported `Unitful` symbols are documented by `Unitful`.
    # persistent_tasks: the check walks every dependency's `Project.toml`, and `InternedStrings` ships none.
    Aqua.test_all(
        GeoParams;
        unbound_args = false,
        piracies = false,
        undocumented_names = false,
        persistent_tasks = false,
    )
end
