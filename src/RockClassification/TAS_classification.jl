module TAS_Classification

# Julia script to retrieve  the classification of igneous rocks using a TAS diagram (Total Alkali (TA) vs Silica (S))
# Positions of the fields are taken from Le Maitre et al., 2002
# The classification is defined between [S = SiO2 = 0-100] and [TA = Na2O + K2O = 0-16]
# NOTE that the oxide composition should be given in weigth%
# 15/05/2022, Nicolas Riel & Boris Kaus

import Base.Threads
using Parameters

export  TAS_ClassificationData, 
        computeClassification!
    
"""
    Declare function to test intersection between two vectors
"""
function testIntersection(v::AbstractArray{_T},p::AbstractArray{_T}) where _T
        n   = 1;

        #check if l2 interset v1
        a1  = v[2,2] - v[1,2];
        b1  = v[1,1] - v[2,1];
        c1  = (v[2,1]*v[1,2]) - (v[1,1]*v[2,2]);

        d1  = (a1 * p[1,1]) + (b1 * p[1,2]) + c1;
        d2  = (a1 * p[2,1]) + (b1 * p[2,2]) + c1;

        if (d1 > 0 && d2 > 0) || (d1 < 0 && d2 < 0)
                n = 0;
        end

        #check if l1 interset v2
        a2  = p[2,2] - p[1,2];
        b2  = p[1,1] - p[2,1];
        c2  = (p[2,1]*p[1,2]) - (p[1,1]*p[2,2]);    

        d1  = (a2 * v[1,1]) + (b2 * v[1,2]) + c2;
        d2  = (a2 * v[2,1]) + (b2 * v[2,2]) + c2;

        if (d1 > 0 && d2 > 0) || (d1 < 0 && d2 < 0)
                n = 0;
        end

        return n;
end

"""
TAS_ClassificationData

Struct that holds default parameters for the TAS diagram
"""
@with_kw_noshow struct TAS_ClassificationData
        
        # name of the TAS field
        litho::Matrix{String} = ["foidite" "picrobasalt" "basanite" "phonotephrite" "tephriphonolite" "phonolite" "basalt" "trachybasalt" "basaltic trachyandesite" "trachyandesite" "trachyte" "basaltic andesite" "andesite" "dacite" "rhyolite"];
        
        # number of vertices per field (polygon)
        n_ver::Matrix{Int64}  = [8 4 6 4 4 4 4 4 4 4 5 4 4 4 5];

        # next array store the vertices coordinates for polygon of all 15 fields defined in litho
        # this array should also be used to draw the fields
        ver::Matrix{Float64}  = [[35 0; 41 0; 41 7; 45 9.4; 48.4 11.5; 52.5 14; 48 16; 35 16];
                                [41 0; 45 0; 45 3; 41 3];
                                [41 3; 45 3; 45 5; 49.4 7.3; 45 9.4; 41 7];
                                [49.4 7.3; 53 9.3; 48.4 11.5; 45 9.4];
                                [53 9.3; 57.6 11.7; 52.5 14; 48.4 11.5];
                                [52.5 14; 57.6 11.7; 65 16; 48 16];
                                [45 0; 52 0; 52 5; 45 5];
                                [45 5; 52 5; 49.4 7.3];
                                [52 5; 57 5.9; 53 9.3; 49.4 7.3];
                                [57 5.9; 63 7; 57.6 11.7; 53 9.3];
                                [63 7; 69 8; 69 16; 65 16; 57.6 11.7];
                                [52 0; 57 0; 57 5.9; 52 5];
                                [57 0; 63 0; 63 7; 57 5.9];
                                [63 0; 77 0; 69 8; 63 7];
                                [77 0; 100 0; 100 16; 69 16; 69 8]                              ];
end


"""
classIndex computeClassification!(chemComp::AbstractArray{_T,N}, ClassTASdata::TAS_ClassificationData)

This compute the classification of the igneous rock using TAS diagram (Total Alkali (TA) vs Silica (S)).

Input:
====
- `chemComp` : vector rock composition in oxide wt%

Output:
- `classIndex` : an integer [0-14] corresponding to a TAS field (TAS_ClassificationData.litho[classIndex])


This routine was developed based the TAS classification of Le Maitre et al., 2002

"""
function computeClassification!(chemComp::AbstractArray{_T,N}, ClassTASdata::TAS_ClassificationData) where {_T,N}
        @unpack litho, n_ver,  ver = ClassTASdata

        # set the classIndex to -1 to track for issue
        classIndex = -1;

        # !!!!!!!!!!!
        # here the composition should be retrieved correctly, shall we retrieve a string array with oxide names? or just retrieve TA and S??
        # TA = chemComp[Na2O]+chemComp[K2O];
        # S  = chemComp[SiO2];
        # point = [TA,S];
        # !!!!!!!!!!!

        n_poly  = size(litho,2);
        point   = [65.0 6.1];

        shift   = 0;
        for poly=1:n_poly
                # shit the index to properly reconstruct the edge of the polygon using ver
                if poly > 1
                        shift  += n_ver[poly-1]
                end

                n      = 0;
                for i  = 1:n_ver[poly]
                      
                        # here we get the coordinates of the polygon edges to be test for intersection
                        x1 = ver[shift+i,1]; y1 = ver[shift+i,2];
                        if (i == n_ver[poly])
                            x2 = ver[shift+1,1]; 
                            y2 = ver[shift+1,2];
                        else
                            x2 = ver[shift+i+1,1]; 
                            y2 = ver[shift+i+1,2];
                        end
                        v  = [x1 y1; x2 y2];

                        # set origin of the ray casting outside the TAS diagram, to make sure we don't start from within a field!
                        p  = [0 0; point[1] point[2]];
                        n += testIntersection(v,p);

                end
                
                # if the number of interestions is uneven then the point is within the tested polygon
                if (mod(n,2) != 0)
                        classIndex = poly;
                        # break the loop, there is no need to test other polygon(s) (TAS fields)
                        break;
                end
        end

        if classIndex == -1
                print("could not find the proper TAS field. Is Na2O + K2O > 16wt%? do you have negative compositions?")
        else

        # print(litho[classIndex])



        return classIndex;

end   

