module TASclassification

# Julia script to retrieve  the classification of igneous rocks using a TAS diagram (Total Alkali (TA) vs Silica (S))
# Positions of the fields are taken from Le Maitre et al., 2002
# The classification is defined between [S = SiO2 = 0-100] and [TA = Na2O + K2O = 0-16]
# NOTE that the oxide composition should be given in weight%
# 15/05/2022, Nicolas Riel & Boris Kaus

import Base.Threads
using Parameters

export TASclassificationData, computeTASclassification, retrieveTASrockType

"""
    Declare function to test intersection between two vectors
"""
function testIntersection(v::AbstractArray{_T}, p::AbstractArray{_T}) where {_T}
    n = 1

    #check if l2 interset v1
    a1 = v[2, 2] - v[1, 2]
    b1 = v[1, 1] - v[2, 1]
    c1 = (v[2, 1] * v[1, 2]) - (v[1, 1] * v[2, 2])

    d1 = (a1 * p[1, 1]) + (b1 * p[1, 2]) + c1
    d2 = (a1 * p[2, 1]) + (b1 * p[2, 2]) + c1

    if (d1 > 0 && d2 > 0) || (d1 < 0 && d2 < 0)
        n = 0
    end

    #check if l1 interset v2
    a2 = p[2, 2] - p[1, 2]
    b2 = p[1, 1] - p[2, 1]
    c2 = (p[2, 1] * p[1, 2]) - (p[1, 1] * p[2, 2])

    d1 = (a2 * v[1, 1]) + (b2 * v[1, 2]) + c2
    d2 = (a2 * v[2, 1]) + (b2 * v[2, 2]) + c2

    if (d1 > 0 && d2 > 0) || (d1 < 0 && d2 < 0)
        n = 0
    end

    return n
end

"""
    TASclassificationData

    Struct that holds default parameters for the TAS diagram
"""
@with_kw_noshow struct TASclassificationData
    # name of the TAS field
    litho::Matrix{String} =
        ["foidite" "picrobasalt" "basanite" "phonotephrite" "tephriphonolite" "phonolite" "basalt" "trachybasalt" "basaltic trachyandesite" "trachyandesite" "trachyte" "basaltic andesite" "andesite" "dacite" "rhyolite"]

    # number of vertices per field (polygon)
    n_ver::Matrix{Int64} = [8 4 6 4 4 4 4 3 4 4 5 4 4 4 5]

    # next array store the vertices coordinates for polygon of all 15 fields defined in litho
    # this array should also be used to draw the fields
    ver::Matrix{Float64} = [
        [
            35.0 0.0
            41.0 0.0
            41.0 7.0
            45.0 9.4
            48.4 11.5
            52.5 14.0
            48.0 16.0
            35.0 16.0
        ]
        [41.0 0.0; 45.0 0.0; 45.0 3.0; 41.0 3.0]
        [41.0 3.0; 45.0 3.0; 45.0 5.0; 49.4 7.3; 45.0 9.4; 41.0 7.0]
        [49.4 7.3; 53.0 9.3; 48.4 11.5; 45.0 9.4]
        [53.0 9.3; 57.6 11.7; 52.5 14.0; 48.4 11.5]
        [52.5 14.0; 57.6 11.7; 65.0 16.0; 48.0 16.0]
        [45.0 0.0; 52.0 0.0; 52.0 5.0; 45.0 5.0]
        [45.0 5.0; 52.0 5.0; 49.4 7.3]
        [52.0 5.0; 57.0 5.9; 53 9.3; 49.4 7.3]
        [57.0 5.9; 63.0 7.0; 57.6 11.7; 53.0 9.3]
        [63.0 7.0; 69.0 8.0; 69.0 16.0; 65.0 16; 57.6 11.7]
        [52.0 0.0; 57.0 0.0; 57.0 5.9; 52.0 5.0]
        [57.0 0.0; 63.0 0.0; 63.0 7.0; 57.0 5.9]
        [63.0 0.0; 77.0 0.0; 69.0 8.0; 63.0 7.0]
        [77.0 0.0; 100.0 0.0; 100.0 16.0; 69.0 16.0; 69.0 8.0]
    ]

    p::Matrix{Float64} = zeros(2, 2)
    v::Matrix{Float64} = zeros(2, 2)
end

"""
retrieveTASrockType(index::Int64; ClassTASdata::TASclassificationData = TASclassificationData())

This retrieves the name of the volcanic rock-type using the computed index

Input:
====
- `index` : integer [1-15]

Output:
- `litho` : a string of the name of corresponding volcanic rock

"""
function retrieveTASrockType(
    index::Int64; ClassTASdata::TASclassificationData=TASclassificationData()
)
    @unpack litho = ClassTASdata
    return litho[index]
end

"""
classIndex computeTASclassification(chemComp::AbstractArray{_T,N}, ClassTASdata::TASclassificationData)

This compute the classification of the igneous rock using TAS diagram (Total Alkali (TA) vs Silica (S)).

Input:
====
- `chemComp` : vector rock composition in oxide wt%

Output:
- `classIndex` : an integer [0-14] corresponding to a TAS field (TASclassificationData.litho[classIndex])

This routine was developed based the TAS classification of Le Maitre et al., 2002

"""
function computeTASclassification(
    point::AbstractArray{_T}; ClassTASdata::TASclassificationData=TASclassificationData()
) where {_T}
    @unpack litho, n_ver, ver, p, v = ClassTASdata

    p[1, 1] = 0.0
    p[2, 1] = point[1]
    p[1, 2] = 0.0
    p[2, 2] = point[2]

    # set the classIndex to -1 to track for issue
    classIndex = -1
    n_poly = size(litho, 2)
    shift = 0
    for poly in 1:n_poly
        # shit the index to properly reconstruct the edge of the polygon using ver
        if poly > 1
            shift += n_ver[poly - 1]
        end
        n = 0
        for i in 1:n_ver[poly]
            # here we get the coordinates of the polygon edges to be test for intersection
            v[1, 1] = ver[shift + i, 1]
            v[1, 2] = ver[shift + i, 2]
            if (i == n_ver[poly])
                v[2, 1] = ver[shift + 1, 1]
                v[2, 2] = ver[shift + 1, 2]
            else
                v[2, 1] = ver[shift + i + 1, 1]
                v[2, 2] = ver[shift + i + 1, 2]
            end

            n += testIntersection(v, p)
        end

        # if the number of interestions is uneven then the point is within the tested polygon
        if (mod(n, 2) != 0)
            classIndex = poly
            # break the loop, there is no need to test other polygon(s) (TAS fields)
            break
        end
    end

    if classIndex == -1
        print(
            "could not find the proper TAS field. Is Na2O + K2O > 16wt%? do you have negative compositions?",
        )
    end

    return classIndex
end

end
