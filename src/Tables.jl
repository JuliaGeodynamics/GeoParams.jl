module Tables

using Unidecode
using GeoParams: AbstractMaterialParam, param_info
using ..Units
using ..MaterialParameters: MaterialParamsInfo

export Phase2Dict, Dict2LatexTable

# was tested with:
# TestPh = (SetMaterialParams(Name="Viscous Matrix", Phase=1, Density=ConstantDensity(),CreepLaws = SetDislocationCreep("Quartz Diorite | Hansen & Carter (1982)")),
#           SetMaterialParams(Name="Viscous Sinker", Phase=2, Density= PT_Density(),CreepLaws = LinearViscous(η=1e21Pa*s)),
#           SetMaterialParams(Name="Viscous Bottom", Phase=3, Density= PT_Density(),CreepLaws = SetDislocationCreep("Diabase | Caristan (1982)")));

"""
Phase2Dict() puts all parameters of a phase in a dict.
Dict2LatexTable() writes .tex file with all parameters from Phase2Dict() output in a table.

"""

function Phase2Dict(s)
    # Dict has Key with Fieldname and Value with Tuple(value, symbol, unit)
    fds = Dict{String,Tuple{String,String,String,String}}()
    refs = Dict{String,Tuple{String,String,String}}()
    # Descriptions for every parameter that could occur in the table and their corresponding variable name(s) that is used in GeoParams
    phasecount = length(s)
    k = 1
    flowlawcount = 0
    # Checks all Phases 
    for i in 1:phasecount
        fieldnames = propertynames(s[i])
        # Goes through all fields of a phase that are not "Name"
        for label in fieldnames
            if !isempty(getproperty(s[i], label)) && label != :Name
                # Goes through all components of all fields
                for j in 1:length(getproperty(s[i], label))
                    a = getproperty(s[i], label)
                    varnames = propertynames(a[j])
                    law = string(typeof(a[1]))
                    flowlaw = ""
                    if occursin("Disl", law)
                        flowlawcount += 1
                        flowlaw = "DislCreep"
                        name = join(a[j].Name)
                        bibinfo_disl = param_info(a[j])
                        bib_disl = bibinfo_disl.BibTex_Reference
                        refs["$name"] = (bib_disl, "$flowlawcount", "$i")
                    elseif occursin("Diff", law)
                        flowlawcount += 1
                        flowlaw = "DiffCreep"
                        name = join(a[j].Name)
                        bibinfo_diff = param_info(a[j])
                        bib_diff = bibinfo_diff.BibTex_Reference
                        refs["$name"] = (bib_diff, "$flowlawcount", "$i")
                    elseif occursin("LinearViscous", law)
                        flowlawcount += 1
                        flowlaw = "LinVisc"
                    end
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
                            fds["$var $label $i"] = (value, latexvar, "$flowlaw", "$i")
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
    return fds, refs
end

function Dict2LatexTable(d::Dict, refs::Dict)
    dictkeys = keys(d)
    symbs = []

    # Creates vectors of type Pairs (can be iterated over, sorted, etc.)
    dictpairs = sort(collect(pairs(d)))
    refpair = sort(collect(pairs(refs)))

    References = ""
    InTextRef = ""

    # Descriptions for every parameter that could occur in the table and their corresponding variable name(s) that is used in GeoParams
    desc = Dict(
        "\\rho" => "Density \$(kg/m^{3})\$",
        "\\rho0" => "Reference density \$(kg/m^{3})\$",
        "g" => "Gravity \$(m/s^{2})\$",
        "\\eta" => "Viscosity \$(Pa \\cdot s)\$",
        "P" => "Pressure (MPa)",
        "T" => "Temperature (C)",
        "V" => "Volume \$(m^{3})\$",
        "d" => "Grain size (cm)",
        "f" => "Water fugacity (MPa)",
        "n" => "Power-law exponent (-)",
        "r" => "Water fugacity exponent (-)",
        "p" => "Grain size exponent (-)",
        "A" => "Coefficient \$(Pa^{-n}/s)\$",
        "E" => "Activation energy (kJ/mol)",
        "R" => "Gas constant (J/mol/K)",
        "G" => "Shear modulus (Pa)",
        "\\nu" => "Poisson ratio (-)",
        "K" => "Bulk modulus (Pa)",
        "Y" => "Young's modulus (Pa)",
        "k" => "Thermal conductivity (W/m/K)",
        "cp" => "Heat capacity (J/kg/K)",
        "Q_L" => "Latent heat (kJ/kg)",
        "H_r" => "Radioactive heat \$(W/m^{3})\$",
        "H_s" => "Shear heating (-)",
        "ϕ" => "Friction angle \$(^{\\circ})\$",
        "\\psi" => "Dilation angle \$(^{\\circ})\$",
        "C" => "Cohesion (Pa)",
        "Vp" => "P-wave velocity (km/s)",
        "Vs" => "S-wave velocity (km/s)",
        "T0" => "Reference temperature (C)",
        "P0" => "Reference pressure (Pa)",
        "\\beta" => "Compressibility (1/Pa)",
        "\\alpha" => "Thermal expansion coeff. (1/K)",
    )

    # Generates latex preamble
    Table = "\\documentclass{article}\n"
    Table *= "\\usepackage{natbib}\n"
    Table *= "\\bibliographystyle{abbrvnat}\n"
    Table *= "\\setcitestyle{authoryear,open={(},close={)}} %Citation-related commands\n"
    Table *= "\\usepackage{booktabs}\n"
    Table *= "\\usepackage{graphicx}\n"
    Table *= "\\usepackage{multirow}\n"
    Table *= "\\usepackage[round, comma, sort, ]{natbib}\n"
    Table *= "\\usepackage[utf8]{inputenc}\n"
    Table *= "\\begin{document}\n"
    Table *= "\\begin{table}[hbt]\n"
    Table *= "\\centering\n"
    Table *= "\\resizebox{\\columnwidth}{!}{%\n"
    Table *= "\\begin{tabular}{ " * "l " * "c "^(parse(Int64, d["Name 1"][2]) + 2) * "}\n"
    Table *= "\\toprule[1pt]\n"
    Table *= "\\midrule[0.3pt]\n"
    Table *= " & "

    # Creates table headers and in text citations
    counter = 1
    for i in 1:parse(Int64, d["Name 1"][2])
        Table *= " & " * d["Name $i"][1]
        for j in 1:length(refs)
            if parse(Int64, refpair[j].second[3]) == i
                Table *= "(" * "*"^counter * ")"
                currentbib = refpair[j].second[1]
                startidx = first(findfirst("{", currentbib))
                endidx = first(findfirst(",", currentbib))
                InTextRef *=
                    "(" *
                    "*"^counter *
                    ") \\cite{" *
                    currentbib[(startidx + 1):(endidx - 1)] *
                    "}"
                if counter != length(refs)
                    InTextRef *= ", "
                end
                counter += 1
            end
        end
    end

    # Latex formating and comment
    Table *= "\\\\\n"
    Table *= "\\hline \n"
    Table *= "% Table body\n"

    # Get vector with all unique symbols without phasenames
    for key in dictkeys
        if occursin("Name", key)
            continue
        else
            push!(symbs, d[key][2])
        end
    end
    symbs = unique(sort(symbs))

    # Creates columnwise output for all parameters of the input phase
    for symbol in symbs
        # Sets parametername and variable
        Table *= " " * string(desc[symbol]) * " & " * "\$" * symbol * "\$"
        # Iterates over all phases
        for j in 1:parse(Int64, d["Name 1"][2])
            hit = 0
            # Iterates over all pairs
            for i in 1:length(dictpairs)
                # Checks if symbol matches the symbol in the pair and if phase matches phase of the pair
                if symbol == dictpairs[i].second[2] &&
                    j == parse(Int64, dictpairs[i].second[4])
                    # put in the matched parameter value
                    Table *= " & " * dictpairs[i].second[1]
                    hit += 1
                end
            end

            # checks if a parameter in a phase is not existing 
            if hit == 0
                Table *= " & "
            end
        end
        # new line in latex
        Table *= " \\\\\n"
    end

    # Adds equations for flow laws underneath the parameters
    DislCreep = 0
    DiffCreep = 0
    LinVisc = 0
    for i in 1:length(dictpairs)
        if dictpairs[i].second[3] == "DislCreep" && DislCreep == 0
            Table *=
                "\\rule[-5pt]{-3pt}{20pt} Dislocation Creep: & " *
                "\\multicolumn{4}{l}{\$ \\dot{\\gamma} = A \\sigma^n f_{H2O}^r \\exp(-\\frac{E+PV}{RT}) \$}\n"
            Table *= " \\\\\n"
            DislCreep += 1
        end
        if dictpairs[i].second[3] == "DiffCreep" && DiffCreep == 0
            Table *=
                "\\rule[-5pt]{-3pt}{20pt} Diffusion Creep: & " *
                "\\multicolumn{4}{l}{\$ \\dot{\\gamma} = A \\sigma^n d^p f_{H2O}^r \\exp(-\\frac{E+PV}{RT}) \$}\n"
            Table *= " \\\\\n"
            DiffCreep += 1
        end
        if dictpairs[i].second[3] == "LinVisc" && LinVisc == 0
            Table *=
                "\\rule[-5pt]{-3pt}{20pt} Linear Viscous: & " *
                "\\multicolumn{4}{l}{\$ \\eta  = \\frac{\\tau_{II} }{ 2\\dot{\\varepsilon_{II}}} \$}\n"
            Table *= " \\\\\n"
            LinVisc += 1
        end
    end

    Table *= "\\midrule[0.3pt]\n"

    # Creates References string which is later written to References.bib 
    for i in 1:parse(Int64, d["Name 1"][2])
        for j in 1:length(refs)
            if parse(Int64, refpair[j].second[2]) == i
                References *= refpair[j].second[1]
            end
        end
    end

    # Adds in text citation for all BibTex sources to LaTEx code. InTextRef is created simultaneously with table headers
    Table *= InTextRef * "\n"

    # finishes latex table, closes all open formats and writes References beneath table
    Table *= "\\midrule[0.3pt]\n"
    Table *= "\\bottomrule[1pt]\n"
    Table *= "\\end{tabular}}\n"
    Table *= "\\caption{Parameter table}\n"
    Table *= "\\label{tab:para_table}\n"
    Table *= "\\end{table}\n"
    Table *= "\\bibliography{References}\n"
    Table *= "\\end{document}\n"

    # Writes BibTex sources in to .bib file and Table string into .tex file
    write("References.bib", References)
    print("\n")
    return write("MaterialParameters.tex", Table)
end

end
