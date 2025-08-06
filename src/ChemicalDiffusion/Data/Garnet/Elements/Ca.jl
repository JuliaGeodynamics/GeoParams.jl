
"""
    Grt_Ca_Carlson2006()

Diffusion data of Ca in garnet calibrated using natural garnets. From Carlson (2006) (https://doi.org/10.2138/am.2006.2043).
"""
function Grt_Ca_Carlson2006()
    data = DiffusionData(
        Name = "Ca diffusion in Garnet (C-O2) | Carlson (2006)",
        Phase = "Garnet",  # name of the mineral
        Formula = "X3Y2(SiO4)3",  # chemical formula of the mineral
        Species = "Ca",  # element or species being diffused
        Orientation = "Isotropic",  # Crystal orientation from the diffusion experiment
        Crystallography = "Isometric",  # Crystallographic system of the mineral
        Buffer = "",  # Buffer condition (e.g., NNO) during the experiment
        Fluid = "",  # Fluid condition (e.g., anhydrous) during the experiment
        D0 = exp(-22.572)u"m^2/s",  # pre-exponential factor
        log_D0_1σ = 0.716NoUnits,
        Ea = 230.56u"kJ/mol",  # activation energy
        Ea_1σ = 7.15u"kJ/mol",  # uncertainty at 1σ of the activation energy
        ΔV = 9.795cm^3 / mol,  # activation volume
        ΔV_1σ = 1.139cm^3 / mol,  # uncertainty at 1σ of the activation volume
        aX = 511.1NoUnits,  # correspond to the sensitivity of the frequency factor to unit-cell dimension (unit normally of 1/nm but here adimensional)
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
    Grt_Ca_Chu2015()

Diffusion data of Ca in garnet calibrated using compilation of experiments data from Cygan and Lasaga (1985), Elphick et al. (1985), Loomis et al. (1985), Chakraborty and Ganguly (1992), Chakraborty and Rubie (1996), Ganguly et al. (1998), Borinski et al. (2012), Vielzeuf et al. (2007), Vielzeuf and Saul (2011) and Borinski et al. (2012).
From Chu and Ague (2015) (https://doi.org/10.1007/s00410-015-1175-y).
"""
function Grt_Ca_Chu2015()
    data = DiffusionData(
        Name = "Ca diffusion in Garnet (C-O2) | Chu and Ague (2015)",
        Phase = "Garnet",  # name of the mineral
        Formula = "X3Y2(SiO4)3",  # chemical formula of the mineral
        Species = "Ca",  # element or species being diffused
        Orientation = "Isotropic",  # Crystal orientation from the diffusion experiment
        Crystallography = "Isometric",  # Crystallographic system of the mineral
        D0 = (10^(-6.36))u"m^2/s",  # pre-exponential factor
        log_D0_1σ = 0.92NoUnits,
        Ea = 299.3u"kJ/mol",  # activation energy
        Ea_1σ = 18.4u"kJ/mol",  # uncertainty at 1σ of the activation energy
        ΔV = 1.77u"J / bar / mol",  # activation volume
        ΔV_1σ = 0.14u"J / bar / mol",  # uncertainty at 1σ of the activation volume
        aX = 302.9NoUnits,  # correspond to the sensitivity of the frequency factor to unit-cell dimension (unit normally of 1/nm but here adimensional)
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