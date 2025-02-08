"""
    Grt_Fe_Chakraborty1992()

Diffusion data of Fe in garnet. Calibrated between 1100-1480°C and 0.14-0.43 GPa in the C-O2 system. From Chakraborty and Ganguly (1992) (https://doi.org/10.1007/BF00296579) combined with data from Loomis et al. (1985) (https://doi.org/10.1007/BF00373040).
"""
function Grt_Fe_Chakraborty1992()
    data = DiffusionData(
        Name = "Fe diffusion in Garnet (C-O2) | Chakraborty and Ganguly (1992)",
        Phase = "Garnet",  # name of the mineral
        Formula = "X3Y2(SiO4)3",  # chemical formula of the mineral
        Species = "Fe",  # element or species being diffused
        Orientation = "Isotropic",  # Crystal orientation from the diffusion experiment
        Crystallography = "Isometric",  # Crystallographic system of the mineral
        Buffer = "C-O2",  # Buffer condition (e.g., NNO) during the experiment
        Fluid = "C-H2O",  # Fluid condition (e.g., anhydrous) during the experiment
        D0 = 6.4 * 1.0e-4cm^2 / s,  # pre-exponential factor
        Ea = 65824u"cal/mol",  # activation energy
        Ea_1σ = 8721u"cal/mol",  # uncertainty at 1σ of the activation energy
        ΔV = 5.6cm^3 / mol,  # activation volume
        ΔV_1σ = 2.9cm^3 / mol,  # uncertainty at 1σ of the activation volume
        Charge = 2,  # charge of the cation
        T_range_min = 1100C,  # temperature min of the experiment
        T_range_max = 1480C,  # temperature max of the experiment
        P0 = 1.0u"bar"  # pressure of calibration
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by HD (18.01.25)",
        BibTex_Reference = "
          @article{chakraborty1992cation,
          title={Cation diffusion in aluminosilicate garnets: experimental determination in spessartine-almandine diffusion couples, evaluation of effective binary diffusion coefficients, and applications},
          author={Chakraborty, Sumit and Ganguly, Jibamitra},
          journal={Contributions to Mineralogy and petrology},
          volume={111},
          number={1},
          pages={74--86},
          year={1992},
          publisher={Springer}
          }
          ",
    )

    return data, info
end
