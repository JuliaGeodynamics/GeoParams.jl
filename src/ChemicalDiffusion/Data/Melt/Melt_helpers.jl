function create_Melt_Holycross2018_data(; Name, Species, Buffer, D0, log_D0_1σ, Ea, Ea_1σ, T_range_min, T_range_max, Charge = 3)
    data = DiffusionData(
        Name = Name,
        Phase = "Melt",
        Formula = "",
        Species = Species,
        Orientation = "Amorphous",
        Crystallography = "Amorphous",
        Buffer = Buffer,
        Fluid = "Hydrated",
        D0 = D0,
        log_D0_1σ = log_D0_1σ,
        Ea = Ea,
        Ea_1σ = Ea_1σ,
        Charge = Charge,
        T_range_min = T_range_min,
        T_range_max = T_range_max,
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by HD (24.01.25). Values originally reported in log10 but converted here to ln. Taken from Table 3 (least-squares fit)",
        BibTex_Reference = "
          @article{holycross2018trace,
          title={Trace element diffusion and kinetic fractionation in wet rhyolitic melt},
          author={Holycross, Megan E and Watson, E Bruce},
          journal={Geochimica et Cosmochimica Acta},
          volume={232},
          pages={14--29},
          year={2018},
          publisher={Elsevier}
          }
          ",
    )

    return data, info
end
