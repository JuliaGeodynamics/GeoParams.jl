using Test
using GeoParams, LinearAlgebra
@testset "TASclassification.jl" begin
    ClassTASdata = TASclassificationData() # default data

    # list of test points S-TA
    point_list = [
        60.0 5.0
        50.0 14.0
        70.0 2.0
        41.0 2.0
        42.0 6.5
        56.0 15.0
        90.0 7.0
        57.0 5.9
        50.0 8.1
        58.0 11.7
        49.0 7.3
        54.0 4.5
        50.0 6.0
        49.0 3.4
        54.0 12.0
        60.0 7.0
        40.0 10.0
        67.0 9.4
        53.1 6.4
        70.5 12.1
        38.0 4.4
    ]

    field = [13 1 14 2 3 6 15 9 4 11 3 12 8 7 5 10 1 11 9 15 1]
    name =
        ["andesite" "foidite" "dacite" "picrobasalt" "basanite" "phonolite" "rhyolite" "basaltic trachyandesite" "phonotephrite" "trachyte" "basanite" "basaltic andesite" "trachybasalt" "basalt" "tephriphonolite" "trachyandesite" "foidite" "trachyte" "basaltic trachyandesite" "rhyolite" "foidite"]

    # add tests to check that results are consistent
    for i in 1:size(point_list, 1)
        index = computeTASclassification(point_list[i, :]; ClassTASdata = ClassTASdata)
        @test index == field[i]
        @test retrieveTASrockType(index; ClassTASdata = ClassTASdata) == name[i]
    end
end
