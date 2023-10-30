
module DislocationCreep
# This contains all the dislocation creep data
using GeoParams
export dislocation_database, dislocation_database_info


"""
Dislocation creep data for dry olivine after Hirth, G. & Kohlstedt (2003)
"""
function DryOlivine_HirthKohlstedt_2003()
    # after Hirth, G. & Kohlstedt (2003), D. Rheology of the upper mantle and the mantle wedge: A view from the experimentalists.
    # Inside the subduction Factory 83?105. Table 1, "dry dislocation" parameters

    data = DislocationCreep(;
        Name="Dry Olivine | Hirth & Kohlstedt (2003)",
        n=3.5NoUnits,
        r=0.0NoUnits,
        A=1.1e5MPa^(-7//2) / s,
        E=530.0kJ / mol,
        V=14e-6m^3 / mol,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Still to be verified with the original publication (BK). Values checked, plots are not reproduced (DK).",
        BibTex_Reference="
            @incollection{eiler_rheology_2003,
            address = {Washington, D. C.},
            title = {Rheology of the upper mantle and the mantle wedge: {A} view from the experimentalists},
            volume = {138},
            isbn = {978-0-87590-997-4},
            shorttitle = {Rheology of the upper mantle and the mantle wedge},
            url = {http://www.agu.org/books/gm/v138/138GM06/138GM06.shtml},
            language = {en},
            urldate = {2019-10-09},
            booktitle = {Geophysical {Monograph} {Series}},
            publisher = {American Geophysical Union},
            author = {Hirth, Greg and Kohlstedt, David},
            editor = {Eiler, John},
            year = {2003},
            doi = {10.1029/138GM06},
            pages = {83--105},
            }
    ",
    )
    return data, info
end


"""
Dislocation creep data for wet olivine (1), constant water fugacity after Hirth, G. & Kohlstedt (2003)
"""
function WetOlivine1_HirthKohlstedt_2003()
    # After Hirth, G. & Kohlstedt (2003), D. Rheology of the upper mantle and the mantle wedge: A view from the experimentalists.
    #   Inside the subduction Factory 83?105. Table 1, "wet dislocation" parameters
    #  Note that this assumes C_OH=1000

    data = DislocationCreep(;
        Name="Wet Olivine 1 | Hirth & Kohlstedt (2003)",
        n=3.5NoUnits,
        A=90MPa^(-7//2) / s,
        E=480kJ / mol,
        V=11e-6m^3 / mol,
        r=1.2NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Still to be verified with the original publication (BK). Values checked, plots are not reproduced (DK).",
        BibTex_Reference="
            @incollection{HirthKohlstedt_OlivineRheology_2003,
            address = {Washington, D. C.},
            title = {Rheology of the upper mantle and the mantle wedge: {A} view from the experimentalists},
            volume = {138},
            isbn = {978-0-87590-997-4},
            shorttitle = {Rheology of the upper mantle and the mantle wedge},
            url = {http://www.agu.org/books/gm/v138/138GM06/138GM06.shtml},
            language = {en},
            urldate = {2019-10-09},
            booktitle = {Geophysical {Monograph} {Series}},
            publisher = {American Geophysical Union},
            author = {Hirth, Greg and Kohlstedt, David},
            editor = {Eiler, John},
            year = {2003},
            doi = {10.1029/138GM06},
            pages = {83--105},
            }
    ",
    )
    return data, info
end


"""
Dislocation creep data for wet olivine (2), 5th creep law in table 1 of Hirth & Kohlstedt (2003)
"""
function WetOlivine2_HirthKohlstedt_2003()
    # After Hirth, G. & Kohlstedt (2003), D. Rheology of the upper mantle and the mantle wedge: A view from the experimentalists.
    #   Inside the subduction Factory 83?105. Table 1, "wet dislocation" parameters
    #  Note that this assumes C_OH=1000

    data =  DislocationCreep(;
        Name="Wet Olivine 2 | Hirth & Kohlstedt (2003)",
        n=3.0NoUnits,
        A=1600MPa^(-3) / s,
        E=520.0kJ / mol,
        V=22.0m^3 / mol,
        r=1.2NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Values checked (NM).",
        BibTex_Reference="
        @incollection{eiler_rheology_2003,
        address = {Washington, D. C.},
        title = {Rheology of the upper mantle and the mantle wedge: {A} view from the experimentalists},
        volume = {138},
        isbn = {978-0-87590-997-4},
        shorttitle = {Rheology of the upper mantle and the mantle wedge},
        url = {http://www.agu.org/books/gm/v138/138GM06/138GM06.shtml},
        language = {en},
        urldate = {2019-10-09},
        booktitle = {Geophysical {Monograph} {Series}},
        publisher = {American Geophysical Union},
        author = {Hirth, Greg and Kohlstedt, David},
        editor = {Eiler, John},
        year = {2003},
        doi = {10.1029/138GM06},
        pages = {83--105},
        }
    ",
    )
    return data, info
end

"""
Dislocation creep data for Quartz Diorite after Hansen (1982), 'Semibrittle creep of selected crustal rocks at 1000 MPa.' and, Hansen & Carter (1982),  'Flow properties of continental lithosphere.'
"""
function QuartzDiorite_HansenCarter_1982()
    #  After Hansen (1982), 'Semibrittle creep of selected crustal rocks at 1000 MPa.' and, Hansen & Carter (1982),
    #  'Flow properties of continental lithosphere.'
    #  Hansen (1982), Fig. 53, page 184 in PDF viewer and table 18, page 224
    #  Carter & Tsenn (1986), table 4, page 18 in PDF viewer

    data =  DislocationCreep(;
        Name="Quartz Diorite | Hansen & Carter (1982)",
        n=2.25NoUnits,
        A=3.5e-2MPa^(-9//4) / s,
        E=212kJ / mol,
        V=0m^3 / mol,
        r=0NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Verified with the original publication Hansen(1982)(NM). Values checked, plots are not reproduced (NM).",
        BibTex_Reference="
            @article{carter1982stress,
            title={Stress magnitudes in natural rock salt},
            author={Carter, Neville L and Hansen, Francis D and Senseny, Paul E},
            journal={Journal of Geophysical Research: Solid Earth},
            volume={87},
            number={B11},
            pages={9289--9300},
            year={1982},
            publisher={Wiley Online Library}
            }
    ",
    )
    return data, info
end


"""
Dislocation creep data for Diabase after Caristan (1982), 'The transition from high temperature creep to fracture in Maryland diabase.' and, Bremond (1999), 'Hydrothermalism and diapirism in the Archean: gravitational instability constraints' Bremond (1999), page 5 in text
"""
function Diabase_Caristan_1982()
    #  After Caristan (1982), 'The transition from high temperature creep to fracture in Maryland diabase.'
    #  and, Bremond (1999),
    #  'Hydrothermalism and diapirism in the Archean: gravitational instability constraints'
    #  Bremond (1999), page 5 in text

    data =  DislocationCreep(;
        Name="Diabase | Caristan (1982)",
        n=3.05NoUnits,
        A=6.0e-2MPa^(-61//20) / s,
        E=276kJ / mol,
        V=0m^3 / mol,
        r=0NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Values checked (Bremond (1999)), plots are not reproduced (NM).",
        BibTex_Reference="
        @article{caristan1982transition,
        title={The transition from high temperature creep to fracture in Maryland diabase},
        author={Caristan, Y},
        journal={Journal of Geophysical Research: Solid Earth},
        volume={87},
        number={B8},
        pages={6781--6790},
        year={1982},
        publisher={Wiley Online Library}
        }
    ",
    )
    return data, info
end

function  dislocation_database(f::F) where F 
    data, _ = f()
    return data
end

function  dislocation_database_info(f::F) where F 
    _, info = f()
    return info
end


end


