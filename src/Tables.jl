# using GeoParams, Unidecode
# was tested with:
# TestPh = (SetMaterialParams(Name="Viscous Matrix", Phase=1, Density=ConstantDensity(),CreepLaws = SetDislocationCreep("Quartz Diorite | Hansen & Carter (1982)")),
#           SetMaterialParams(Name="Viscous Sinker", Phase=2, Density= PT_Density(),CreepLaws = LinearViscous(η=1e21Pa*s)),
#           SetMaterialParams(Name="Viscous Bottom", Phase=3, Density= PT_Density(),CreepLaws = SetDislocationCreep("Diabase | Caristan (1982)")));
#
#
#=
function Phase2Dict(s)
    # Dict has Key with Fieldname and Value with Tuple(value, symbol, unit)
    fds = Dict{String, Tuple{String, String, String, String}}()
    # Descriptions for every parameter that could occur in the table and their corresponding variable name(s) that is used in GeoParams
    phasecount = length(s)
    k = 1
    # Checks all Phases 
    for i = 1:phasecount
        fieldnames = propertynames(s[i])
        # Goes through all fields of a phase that are not "Name"
        for label in fieldnames
            if !isempty(getproperty(s[i], label)) && label != :Name
                # Goes through all components of all fields
                for j = 1:length(getproperty(s[i], label))
                    a = getproperty(s[i], label)
                    varnames = propertynames(a[j])
                    # Goes through all variables in a field component
                    for var in varnames
                        b = getproperty(a[j], var)
                        # Checks if they are GeoUnit (Value and Unit exist then)
                        if isa(b, GeoUnit)
                            unit = string(getproperty(b, :unit))
                            value = string(getproperty(b, :val))
                            # Gives back LaTex format string for the corresponding variable if it is longer than 2 chars
                            if length(unidecode("$var")) > 2
                                latexvar = string("\\" * string(unidecode("$var")))
                            else
                                latexvar = string("$var")
                            end
                            # Put value, LaTex name and units in a Dict
                            fds["$var $label $i"] = (value, latexvar, unit, "$i")
                        end
                        k += 1
                    end
                end
            # Takes field "Name" and puts it in the first entry of the Tuple, takes Phasecount of all Phases and puts it in the second entry of the Dict
            elseif !isempty(getproperty(s[i], label)) && label == :Name
                phasename = join(getproperty(s[i], :Name))
                fds["$label $i"] = (phasename, "$phasecount", "$i", "")
            end
        end
    end
    return fds
end




function Dict2LatexTable(d::Dict)
    dictkeys = keys(d)
    symbs = []

    # Descriptions for every parameter that could occur in the table and their corresponding variable name(s) that is used in GeoParams
    desc = Dict("\\rho"=>"Density \$(kg/m^{3})\$","\\rho0"=>"Reference density \$(kg/m^{3})\$","g"=>"Gravity \$(m/s^{2})\$","\\eta"=>"Viscosity \$(Pa \\cdot s)\$",
    "P"=>"Pressure (MPa)","T"=>"Temperature (C)","V"=>"Volume \$(m^{3})\$" ,"d"=>"Grain size (cm)" ,"f"=>"Water fugacity (MPa)",
    "n"=>"Power-law exponent (-)","r"=>"Water fugacity exponent (-)" ,"p"=>"Grain size exponent (-)" ,"A"=>"Coefficient \$(Pa^{-n}/s)\$" ,
    "E"=>"Activation energy (kJ/mol)" ,"R"=>"Gas constant (J/mol/K)" ,"G"=>"Shear modulus (Pa)","\\nu"=>"Poisson ratio (-)",
    "K"=>"Bulk modulus (Pa)" , "Y"=>"Young's modulus (Pa)","k"=>"Thermal conductivity (W/m/K)","cp"=>"Heat capacity (J/kg/K)",
    "Q_L"=>"Latent heat (kJ/kg)","H_r"=>"Radioactive heat \$(W/m^{3})\$","H_s"=>"Shear heating (-)","ϕ"=>"Friction angle \$(^{\\circ})\$",
    "\\psi"=>"Dilation angle \$(^{\\circ})\$","C"=>"Cohesion (Pa)" ,"Vp"=>"P-wave velocity (km/s)","Vs"=>"S-wave velocity (km/s)" ,
    "T0"=>"Reference temperature (C)","P0"=>"Reference pressure (Pa)","\\beta"=>"Compressibility (1/Pa)","\\alpha"=>"Thermal expansion coeff. (1/K)");

    # Generates latex preamble
    Table  = "\\documentclass{article}\n";
    Table *= "\\usepackage{booktabs}\n";
    Table *= "\\usepackage{graphicx}\n";
    Table *= "\\usepackage[utf8]{inputenc}\n";
    Table *= "\\begin{document}\n";
    Table *= "\\begin{table}[hbt]\n";
    Table *= "\\centering\n"
    Table *= "\\resizebox{\\columnwidth}{!}{%\n"
    Table *= "\\begin{tabular}{ " * "l " * "c " ^ (parse(Int64, d["Name 1"][2])+2) * "}\n";
    Table *= "\\toprule[1pt]\n";
    Table *= "\\midrule[0.3pt]\n";
    Table *= " & ";

    # Creates table headers
    for j = 1:parse(Int64, d["Name 1"][2])
        Table *= " & " * d["Name $j"][1]
    end

    # Latex formating and comment
    Table *= "\\\\\n";
    Table *= " \\hline \n";
    Table *= "    % Table body\n";

    # Get vector with all unique symbols without phasenames
    for key in dictkeys
        if occursin("Name", key)
            continue
        else
            push!(symbs, d[key][2])
        end
    end
    symbs = unique(sort(symbs))
    # Creates a vector of Pairs (can be iterated over)
    dictpairs = sort(collect(pairs(d)))

    # Creates columnwise output for all parameters of the input phase
    for symbol in symbs
        # Sets parametername and variable
        Table *= " " * string(desc[symbol]) * " & " * "\$" * symbol * "\$";
        # Iterates over all phases
        for j = 1:parse(Int64, d["Name 1"][2])
            hit = 0
            # Iterates over all pairs
            for i = 1:length(dictpairs)
                # Checks if symbol matches the symbol in the pair and if phase matches phase of the pair
                if symbol == dictpairs[i].second[2] && j == parse(Int64, dictpairs[i].second[4])
                    # put in the matched parameter value
                    Table *= " & " * dictpairs[i].second[1];
                    hit += 1
                end
            end

            # checks if a parameter in a phase is not existing 
            if hit == 0
                Table *= " & ";
            end

        end
        # new line in latex
        Table *= " \\\\\n";
    end

    # finishes latex table, closes all open formats
    Table *= "\\toprule[1pt]\n";
    Table *= "\\midrule[0.3pt]\n";
    Table *= "\\end{tabular}}\n";
    Table *= "\\caption{Parameter table}\n"
    Table *= "\\label{tab:para_table}\n";
    Table *= "\\end{table}\n"
    Table *= "\\end{document}\n";

    # Writes result into a .tex file
    write("MaterialParameters.tex", Table);
end
=#
