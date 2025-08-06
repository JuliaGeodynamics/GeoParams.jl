"""
    Grt_Mn_Chakraborty1992()

Diffusion data of Mn in garnet. Calibrated between 1100-1480°C and 0.14-0.43 GPa in the C-O2 system. From Chakraborty and Ganguly (1992) (https://doi.org/10.1007/BF00296579) combined with data from Loomis et al. (1985) (https://doi.org/10.1007/BF00373040).
"""
function Grt_Mn_Chakraborty1992()
    data = DiffusionData(
        Name = "Mn diffusion in Garnet (C-O2) | Chakraborty and Ganguly (1992)",
        Phase = "Garnet",  # name of the mineral
        Formula = "X3Y2(SiO4)3",  # chemical formula of the mineral
        Species = "Mn",  # element or species being diffused
        Orientation = "Isotropic",  # Crystal orientation from the diffusion experiment
        Crystallography = "Isometric",  # Crystallographic system of the mineral
        Buffer = "C-O2",  # Buffer condition (e.g., NNO) during the experiment
        Fluid = "C-H2O",  # Fluid condition (e.g., anhydrous) during the experiment
        D0 = 5.1 * 1.0e-4cm^2 / s,  # pre-exponential factor
        Ea = 60569u"cal/mol",  # activation energy
        Ea_1σ = 8889u"cal/mol",  # uncertainty at 1σ of the activation energy
        ΔV = 6.0cm^3 / mol,  # activation volume
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


"""
    Grt_Mn_Carlson2006()

Diffusion data of Mn in garnet calibrated using natural garnets. From Carlson (2006) (https://doi.org/10.2138/am.2006.2043).
"""
function Grt_Mn_Carlson2006()
    data = DiffusionData(
        Name = "Mn diffusion in Garnet (C-O2) | Carlson (2006)",
        Phase = "Garnet",  # name of the mineral
        Formula = "X3Y2(SiO4)3",  # chemical formula of the mineral
        Species = "Mn",  # element or species being diffused
        Orientation = "Isotropic",  # Crystal orientation from the diffusion experiment
        Crystallography = "Isometric",  # Crystallographic system of the mineral
        Buffer = "",  # Buffer condition (e.g., NNO) during the experiment
        Fluid = "",  # Fluid condition (e.g., anhydrous) during the experiment
        D0 = exp(-17.619)u"m^2/s",  # pre-exponential factor
        log_D0_1σ = 0.876NoUnits,
        Ea = 264.44u"kJ/mol",  # activation energy
        Ea_1σ = 5.99u"kJ/mol",  # uncertainty at 1σ of the activation energy
        ΔV = 9.631cm^3 / mol,  # activation volume
        ΔV_1σ = 1.110cm^3 / mol,  # uncertainty at 1σ of the activation volume
        aX = 496.9NoUnits,  # correspond to the sensitivity of the frequency factor to unit-cell dimension (unit normally of 1/nm but here adimensional)
        bX = - 1.1525NoUnits,  # unit-cell dimension of the almandine endmember (unit normally of nm but here adimensional)
        nfO2 = (1/6)NoUnits,  # exponent for the oxygen fugacity dependency
        dfO2 = (1e-25)NoUnits,  # quotient for f(O2) dependency
        Charge = 2,  # charge of the cation
        T_range_min = 500C,  # temperature min of the experiment
        T_range_max = 1500C,  # temperature max of the experiment
        P0 = 0.0u"bar"  # pressure of calibration
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by HD (06.08.25)",
        BibTex_Reference = "
            @article{carlson2006rates,
            title={Rates of Fe, Mg, Mn, and Ca diffusion in garnet},
            author={Carlson, William D},
            journal={American Mineralogist},
            volume={91},
            number={1},
            pages={1--11},
            year={2006},
            publisher={Mineralogical Society of America}
            }
          ",
    )

    return data, info
end


"""
    Grt_Mn_Chu2015()

Diffusion data of Mn in garnet calibrated using compilation of experiments data from Cygan and Lasaga (1985), Elphick et al. (1985), Loomis et al. (1985), Chakraborty and Ganguly (1992), Chakraborty and Rubie (1996), Ganguly et al. (1998), Borinski et al. (2012), Vielzeuf et al. (2007), Vielzeuf and Saul (2011) and Borinski et al. (2012).
From Chu and Ague (2015) (https://doi.org/10.1007/s00410-015-1175-y).
"""
function Grt_Mn_Chu2015()
    data = DiffusionData(
        Name = "Mn diffusion in Garnet (C-O2) | Chu and Ague (2015)",
        Phase = "Garnet",  # name of the mineral
        Formula = "X3Y2(SiO4)3",  # chemical formula of the mineral
        Species = "Mn",  # element or species being diffused
        Orientation = "Isotropic",  # Crystal orientation from the diffusion experiment
        Crystallography = "Isometric",  # Crystallographic system of the mineral
        D0 = (10^(-9.86))u"m^2/s",  # pre-exponential factor
        log_D0_1σ = 0.70NoUnits,
        Ea = 212.9u"kJ/mol",  # activation energy
        Ea_1σ = 13.4u"kJ/mol",  # uncertainty at 1σ of the activation energy
        ΔV = 0.95u"J / bar / mol",  # activation volume
        ΔV_1σ = 0.13u"J / bar / mol",  # uncertainty at 1σ of the activation volume
        aX = 686.7NoUnits,  # correspond to the sensitivity of the frequency factor to unit-cell dimension (unit normally of 1/nm but here adimensional)
        bX = - 1.1525NoUnits,  # unit-cell dimension of the almandine endmember (unit normally of nm but here adimensional)
        nfO2 = (1/6)NoUnits,  # exponent for the oxygen fugacity dependency
        dfO2 = (1e-25)NoUnits,  # quotient for f(O2) dependency
        Charge = 2,  # charge of the cation
        T_range_min = 500C,  # temperature min of the experiment
        T_range_max = 1500C,  # temperature max of the experiment
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by HD (06.08.25)",
        BibTex_Reference = "
            @article{chu2015analysis,
            title={Analysis of experimental data on divalent cation diffusion kinetics in aluminosilicate garnets with application to timescales of peak Barrovian metamorphism, Scotland},
            author={Chu, Xu and Ague, Jay J},
            journal={Contributions to Mineralogy and Petrology},
            volume={170},
            number={2},
            pages={25},
            year={2015},
            publisher={Springer}
            }
          ",
    )

    return data, info
end