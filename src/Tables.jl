module Tables

using Unidecode
using GeoParams: AbstractMaterialParam, param_info
using ..Units
using ..MaterialParameters: MaterialParamsInfo

export detachFloatfromExponent, Phase2Dict, Dict2LatexTable, Phase2DictMd, Dict2MarkdownTable, ParameterTable

# was tested with:
# TestPh = (SetMaterialParams(Name="Viscous Matrix", Phase=1, Density=ConstantDensity(),CreepLaws = SetDislocationCreep("Quartz Diorite | Hansen & Carter (1982)")),
#           SetMaterialParams(Name="Viscous Sinker", Phase=2, Density= PT_Density(),CreepLaws = LinearViscous(η=1e21Pa*s)),
#           SetMaterialParams(Name="Viscous Bottom", Phase=3, Density= PT_Density(),CreepLaws = SetDislocationCreep("Diabase | Caristan (1982)")))

""" 
detachFloatfromExponent() returns the number of float decimals after the comma as Integer, the float number without the exponent as string 
and the exponent after the "e" as string.
The argument output returns "1" for "ex" if the input number has no exponent.

"""

function detachFloatfromExponent(str::String)
    s = lowercase(str)
    if 'e' in s
        s, ex = split(s, "e")
    else
        ex = "1"
    end

    dig = something(findfirst('.', reverse(s)), 1) - 1
    
    return dig, s, ex

end


"""
Phase2Dict() puts all parameters of a phase in a dict.

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

"""
Dict2LatexTable() writes a .tex file with all parameters from the Phase2Dict() output in a LaTeX table.

"""

function Dict2LatexTable(d::Dict, refs::Dict; filename=nothing, rdigits=4)
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
        "T" => "Temperature (\$^{\\circ}\$C)",
        "V" => "Volume \$(m^{3})\$",
        "d" => "Grain size (cm)",
        "f" => "Water fugacity (MPa)",
        "n" => "Power-law exponent (-)",
        "r" => "Water fugacity exponent (-)",
        "p" => "Grain size exponent (-)",
        "A" => "Coefficient \$(Pa^{-n}/s)\$",
        "E" => "Activation energy (J/mol)",
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
        "T0" => "Reference temperature (\$^{\\circ}\$C)",
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

                Table *= "(" * "*" ^ counter * ")"
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
    Table *= "\\midrule \n"
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
        Table *= " " * string(desc[symbol]) * " & " * "\$" * symbol[1:end-1] * "_" * symbol[end] *"\$"
        # Iterates over all phases
        for j in 1:parse(Int64, d["Name 1"][2])
            hit = 0
            # Iterates over all pairs
            for i in 1:length(dictpairs)
                # Checks if symbol matches the symbol in the pair and if phase matches phase of the pair
                if symbol == dictpairs[i].second[2] &&
                    j == parse(Int64, dictpairs[i].second[4])
                    # put in the matched parameter value
                    dig, num, expo = detachFloatfromExponent(dictpairs[i].second[1])
                    if  dig <= rdigits && expo != "1" 
                        Table *= " & " * "$num \\times 10^{$expo}"
                    elseif dig <= rdigits && expo == "1"
                        Table *= " & " * "$num"
                    elseif dig > rdigits && expo != "1"
                        Table *= " & " * string(round(parse(Float64, num); digits=rdigits)) * " \\times 10^{$expo}"
                    elseif dig > rdigits && expo == "1"
                        Table *= " & " * string(round(parse(Float64, num); digits=rdigits))
                    end
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
                "\\multicolumn{4}{l}{\$ \\dot{\\gamma} = A \\tau^n f_{H2O}^r \\exp(-\\frac{E+PV}{RT}) \$}\n"
            Table *= " \\\\\n"
            DislCreep += 1
        end
        if dictpairs[i].second[3] == "DiffCreep" && DiffCreep == 0
            Table *=
                "\\rule[-5pt]{-3pt}{20pt} Diffusion Creep: & " *
                "\\multicolumn{4}{l}{\$ \\dot{\\gamma} = A \\tau^n d^p f_{H2O}^r \\exp(-\\frac{E+PV}{RT}) \$}\n"
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
    return write("References.bib", References), write("$filename.tex", Table)
end

"""
Phase2DictMd() puts all parameters of a phase in a dict.

"""

function Phase2DictMd(s)
    # Dict has Key with Fieldname and Value with Tuple(value, symbol, unit)
    fds = Dict{String, Tuple{String, String, String, String}}()
    refs = Dict{String, Tuple{String, String, String}}()
    # Descriptions for every parameter that could occur in the table and their corresponding variable name(s) that is used in GeoParams
    phasecount = length(s)
    k = 1
    flowlawcount = 0
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
                            # Gives back variablename 
                            latexvar = "$var"
                            # Put value, LaTex name and units in a Dict
                            fds["$var $label $i"] = (value, latexvar, "", "$i") # 3te Stelle sollte sein: "$flowlaw" statt ""
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

"""
Dict2MarkdownTable() writes a .md file with all parameters from the Phase2DictMd() output in a Markdown table.

"""

function Dict2MarkdownTable(d::Dict; filename=nothing, rdigits=4)
    dictkeys = keys(d)
    symbs = []
    #refs = Dict{}()

    # Creates vectors of type Pairs (can be iterated over, sorted, etc.)
    dictpairs = sort(collect(pairs(d)))
    #refpair = sort(collect(pairs(refs)))

    #References = ""
    #InTextRef = ""

    # Descriptions for every parameter that could occur in the table and their corresponding variable name(s) that is used in GeoParams
    desc = Dict(
        "ρ"=>"Density (kg/m^3^)",
        "ρ0"=>"Reference density (kg/m^3^)",
        "g"=>"Gravity (m/s^2^)","η"=>"Viscosity (Pa s)",
        "P"=>"Pressure (MPa)","T"=>"Temperature (°C)",
        "V"=>"Volume (m^3^)" ,"d"=>"Grain size (cm)",
        "f"=>"Water fugacity (MPa)",
        "n"=>"Power-law exponent (-)",
        "r"=>"Water fugacity exponent (-)",
        "p"=>"Grain size exponent (-)",
        "A"=>"Coefficient (Pa^-n^/s)",
        "E"=>"Activation energy (J/mol)",
        "R"=>"Gas constant (J/mol/K)",
        "G"=>"Shear modulus (Pa)",
        "ν"=>"Poisson ratio (-)",
        "K"=>"Bulk modulus (Pa)", 
        "Y"=>"Young's modulus (Pa)",
        "k"=>"Thermal conductivity (W/m/K)",
        "cp"=>"Heat capacity (J/kg/K)",
        "Q_L"=>"Latent heat (kJ/kg)",
        "H_r"=>"Radioactive heat (W/m^3^)",
        "H_s"=>"Shear heating (-)",
        "ϕ"=>"Friction angle (°)",
        "ψ"=>"Dilation angle (°)",
        "C"=>"Cohesion (Pa)",
        "Vp"=>"P-wave velocity (km/s)",
        "Vs"=>"S-wave velocity (km/s)",
        "T0"=>"Reference temperature (°C)",
        "P0"=>"Reference pressure (Pa)",
        "β"=>"Compressibility (1/Pa)",
        "α"=>"Thermal expansion coeff. (1/K)"
        )

    # Generates latex preamble
    Table = " | Denotation | Variablename | "

    # Creates table headers and in text citations
    counter = 1
    for i = 1:parse(Int64, d["Name 1"][2])
        Table *=  d["Name $i"][1] * " | "
        # Ab hier sinds die References für am Ende -> Später machen!
        #for j = 1:length(refs)
            #if parse(Int64, refpair[j].second[3]) == i
                #Table *= "(" * "*" ^ counter * ") |"

                
                #currentbib = refpair[j].second[1]
                #startidx = first(findfirst("{", currentbib))
                #endidx = first(findfirst(",", currentbib))
                #InTextRef *= "(" * "*" ^ counter * ") \\cite{" * currentbib[startidx+1:endidx-1] * "}"
                #if counter != length(refs)
                #    InTextRef *= ", "
                #end
                #counter += 1
            #else
                #Table *= " | "
            #end
        #end
    end

    Table *= "\n"
    Table *= " | ---------- | ------------ | "

    for i = 1:parse(Int64, d["Name 1"][2])
        Table *= "-" ^ length(d["Name $i"][1]) * " | "
    end

    Table *= "\n"

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
        if length(symbol) > 1
            #checks if Unicode chars in combination with a number etc. (e.g "α0") are used which apparently do not have a second index but only first and third
            if length(unidecode(symbol)) > 2
                Table *= " " * string(desc[symbol]) * " | " * symbol[1] * "~" * symbol[3] * "~"
            else
                Table *= " " * string(desc[symbol]) * " | " * symbol[1] * "~" * symbol[2] * "~"
            end
        else
            Table *= " " * string(desc[symbol]) * " | " * symbol
        end
        # Iterates over all phases
        for j = 1:parse(Int64, d["Name 1"][2])
            hit = 0
            # Iterates over all pairs
            for i = 1:length(dictpairs)
                # Checks if symbol matches the symbol in the pair and if phase matches phase of the pair
                if symbol == dictpairs[i].second[2] && j == parse(Int64, dictpairs[i].second[4])
                    # put in the matched parameter value
                    dig, num, expo = detachFloatfromExponent(dictpairs[i].second[1])
                    if  dig <= rdigits && expo != "1" 
                        Table *= " | " * "$num x 10^$expo^"
                    elseif dig <= rdigits && expo == "1"
                        Table *= " | " * "$num"
                    elseif dig > rdigits && expo != "1"
                        Table *= " | " * string(round(parse(Float64, num); digits=rdigits)) * " x 10^$expo^"
                    elseif dig > rdigits && expo == "1"
                        Table *= " | " * string(round(parse(Float64, num); digits=rdigits))
                    end
                    hit += 1
                end
            end

            # checks if a parameter in a phase is not existing 
            if hit == 0
                Table *= " | "
            end

        end
        # new line in latex
        Table *= " \n"
    end

    #=

    # Adds equations for flow laws underneath the parameters
    DislCreep = 0
    DiffCreep = 0
    LinVisc = 0
    for i = 1:length(dictpairs)
        if dictpairs[i].second[3] == "DislCreep" && DislCreep == 0
            Table *= "\\rule[-5pt]{-3pt}{20pt} Dislocation Creep: & " * "\\multicolumn{4}{l}{ \\dot{\\gamma} = A \\sigma^n f_{H2O}^r \\exp(-\\frac{E+PV}{RT}) }\n"
            Table *= " \\\\\n"
            DislCreep += 1
        end
        if dictpairs[i].second[3] == "DiffCreep" && DiffCreep == 0
            Table *= "\\rule[-5pt]{-3pt}{20pt} Diffusion Creep: & " * "\\multicolumn{4}{l}{ \\dot{\\gamma} = A \\sigma^n d^p f_{H2O}^r \\exp(-\\frac{E+PV}{RT}) }\n"
            Table *= " \\\\\n"
            DiffCreep += 1
        end
        if dictpairs[i].second[3] == "LinVisc" && LinVisc == 0
            Table *= "\\rule[-5pt]{-3pt}{20pt} Linear Viscous: & " * "\\multicolumn{4}{l}{ \\eta  = \\frac{\\tau_{II} }{ 2\\dot{\\varepsilon_{II}}} }\n"
            Table *= " \\\\\n"
            LinVisc += 1
        end
    end
    
    
    # Creates References string which is later written to References.bib 
    for i = 1:parse(Int64, d["Name 1"][2])
        for j = 1:length(refs)
            if parse(Int64, refpair[j].second[2]) == i
                References *= refpair[j].second[1]
            end
        end
    end

    # Adds in text citation for all BibTex sources to LaTEx code. InTextRef is created simultaneously with table headers
    Table *= InTextRef * "\n"
    =#

    # Writes BibTex sources in to .bib file and Table string into .m file
    #write("References.bib", References)
    return write("$filename.md", Table)
end


"""
ParameterTable() creates a table with all parameters saved in the Phase struct. It lets you choose between "latex" and "markdown" as table formats with LaTeX as default.
Creates a filename.tex or filename.md file as output. If "latex" is chosen a "Reference.bib" file will automatically be produced with all your references.
"""
function ParameterTable(
    Phase;
    filename=nothing,
    format="latex",
    rdigits=4
    )

    if format == "latex"
        d, ref = Phase2Dict(Phase)
        Dict2LatexTable(d, ref, filename=filename)
    elseif format == "markdown"
        d = Phase2DictMd(Phase)
        Dict2MarkdownTable(d)
    end
    return nothing
end

end
