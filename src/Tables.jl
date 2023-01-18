module Tables

using Unidecode
using GeoParams: AbstractMaterialParam, param_info, LinearViscous, PowerlawViscous, DislocationCreep, DiffusionCreep, CompositeRheology, Parallel
using ..Units
using ..MaterialParameters: MaterialParamsInfo

export detachFloatfromExponent, Phase2Dict, Dict2LatexTable, Phase2DictMd, Dict2MarkdownTable, ParameterTable

# was tested with:
# v1 = SetDiffusionCreep("Dry Anorthite | Rybacki et al. (2006)")
# c1 = CompositeRheology(v1, SetDislocationCreep("Diabase | Caristan (1982)"), LinearViscous())
# TestPh = (SetMaterialParams(Name="Viscous Matrix", Phase=1, Density=ConstantDensity(),CreepLaws = SetDislocationCreep("Quartz Diorite | Hansen & Carter (1982)")),
#           SetMaterialParams(Name="Viscous Sinker", Phase=2, Density= PT_Density(),CompositeRheology = c1),
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
    # Dict has Key with Fieldname and Value with Tuple(value, symbol, flowlaw keyword, number of phases)
    fds = Dict{String,Tuple{String,String,String,String,String,String}}()
    refs = Dict{String,Tuple{String,String,String}}()
    if !(typeof(s) <: Tuple)
        s = (s,)
    end
    phasecount = length(s)
    k = 1
    flowlawcount = 0
    global Diff = 0
    global Disl = 0
    global Lin = 0
    var_counter = Dict(
        "A" => 0,
        "E" => 0,
        "R" => 0,
        "V" => 0,
        "η" => 0,
        "n" => 0,
        "p" => 0,
        "r" => 0,
    )
    # Checks all Phases 
    for i in 1:phasecount
        fieldnames = propertynames(s[i])
        global num_rheologies = 0
        for label in fieldnames
            if !isempty(getproperty(s[i], label)) && label != :Name
                # Goes through all components of all fields
                for j in 1:length(getproperty(s[i], label))
                    a = getproperty(s[i], label)
                    flowlaw = ""
                    flowdisl = "DislCreep"
                    flowdiff = "DiffCreep"
                    flowlin = "LinVisc"
                    flowcomp = "CompoRheo"
                    flowpara = "Parallel"
                    # Checks what type the CreepLaw or CompositeRheology field has 
                    if typeof(a[j]) <: DislocationCreep
                        global Disl += 1
                        global num_rheologies += 1
                        flowlawcount += 1
                        flowadd = flowdisl
                        flowlaw *= flowadd
                        name = join(a[j].Name)
                        bibinfo_disl = param_info(a[j])
                        bib_disl = bibinfo_disl.BibTex_Reference
                        refs["$name"] = (bib_disl, "$flowlawcount", "$i")

                    elseif typeof(a[j]) <: DiffusionCreep
                        global Diff += 1
                        global num_rheologies += 1
                        flowlawcount += 1
                        flowadd = flowdiff
                        flowlaw *= flowadd
                        name = join(a[j].Name)
                        bibinfo_diff = param_info(a[j])
                        bib_diff = bibinfo_diff.BibTex_Reference
                        refs["$name"] = (bib_diff, "$flowlawcount", "$i")

                    elseif typeof(a[j]) <: LinearViscous
                        global Lin += 1
                        global num_rheologies += 1
                        flowlawcount += 1
                        flowadd = flowlin
                        flowlaw *= flowadd

                    elseif typeof(a[j]) <: CompositeRheology
                        compos = getproperty(a[j], :elements)
                        global num_rheologies += length(compos)
                        fdsname = "Comp "
                        flowadd = flowcomp
                        comporheo = flowadd * "("
                        for u in 1:num_rheologies
                            if typeof(a[j][u]) <: Parallel
                                num_parallel = length(a[j][u].elements)
                                flowadd = flowpara
                                parallelrheo = flowadd * "("
                                namepre = fdsname * "Para "
                                for v in 1:num_parallel
                                    if typeof(a[j][u][v]) <: DislocationCreep
                                        global Disl += 1
                                        flowlawcount += 1
                                        flowadd = flowdisl
                                        parallelrheo *= flowadd * ","
                                        name = namepre * join(a[j][u][v].Name)
                                        bibinfo_disl = param_info(a[j][u][v])
                                        bib_disl = bibinfo_disl.BibTex_Reference
                                        refs["$name"] = (bib_disl, "$flowlawcount", "$i")

                                    elseif typeof(a[j][u][v]) <: DiffusionCreep
                                        global Diff += 1
                                        flowlawcount += 1
                                        flowadd = flowdiff
                                        parallelrheo *= flowadd * ","
                                        name = namepre * join(a[j][u][v].Name)
                                        bibinfo_diff = param_info(a[j][u][v])
                                        bib_diff = bibinfo_diff.BibTex_Reference
                                        refs["$name"] = (bib_diff, "$flowlawcount", "$i")

                                    elseif typeof(a[j][u][v]) <: LinearViscous
                                        global Lin += 1
                                        flowlawcount += 1
                                        flowadd = flowlin
                                        parallelrheo *= flowadd * ","
                                        
                                    end

                                    varnames = propertynames(a[j][u][v])
                                    # Goes through all variables in a field component
                                    for var in varnames
                                        b = getproperty(a[j][u][v], var)
                                        # Checks if they are GeoUnit (Value and Unit exist then)
                                        if isa(b, GeoUnit)
                                            value = string(getproperty(b, :val))
                                            if typeof(a[j]) <: CompositeRheology
                                                var_counter["$var"] += 1
                                            end
                                            # Gives back LaTex format string for the corresponding variable if it is longer than 2 chars
                                            if length(unidecode("$var")) > 2
                                                latexvar = string("\\" * string(unidecode("$var")))
                                            else
                                                latexvar = string("$var")
                                            end
                                            # Put value, LaTex variable name, creep law pattern, phase id and current disl, diff or linvisc count in a Dict
                                            counter = var_counter["$var"]
                                            if typeof(a[j][u][v]) <: DislocationCreep
                                                fds["$var $label $flowadd $i.$u.$v"] = (value, latexvar, "$flowlaw$comporheo$parallelrheo)", "$i", "$counter", "$num_rheologies")
                                            end
                                            if typeof(a[j][u][v]) <: DiffusionCreep
                                                fds["$var $label $flowadd $i.$u.$v"] = (value, latexvar, "$flowlaw$comporheo$parallelrheo)", "$i", "$counter", "$num_rheologies")
                                            end
                                            if typeof(a[j][u][v]) <: LinearViscous
                                                fds["$var $label $flowadd $i.$u.$v"] = (value, latexvar, "$flowlaw$comporheo$parallelrheo)", "$i", "$counter", "$num_rheologies")
                                            end
                                        end
                                        k += 1
                                    end
                                end
                                comporheo *= parallelrheo * "),"
                                global num_rheologies += num_parallel - 1

                            elseif typeof(a[j][u]) <: DislocationCreep
                                global Disl += 1
                                flowlawcount += 1
                                flowadd = flowdisl
                                comporheo *= flowadd * ","
                                name =  fdsname * join(a[j][u].Name)
                                bibinfo_disl = param_info(a[j][u])
                                bib_disl = bibinfo_disl.BibTex_Reference
                                refs["$name"] = (bib_disl, "$flowlawcount", "$i")

                            elseif typeof(a[j][u]) <: DiffusionCreep
                                global Diff += 1
                                flowlawcount += 1
                                flowadd = flowdiff
                                comporheo *= flowadd * ","
                                name =  fdsname * join(a[j][u].Name)
                                bibinfo_diff = param_info(a[j][u])
                                bib_diff = bibinfo_diff.BibTex_Reference
                                refs["$name"] = (bib_diff, "$flowlawcount", "$i")

                            elseif typeof(a[j][u]) <: LinearViscous
                                global Lin += 1
                                flowlawcount += 1
                                flowadd = flowlin
                                comporheo *= flowadd * ","
                            end

                            varnames = propertynames(a[j][u])
                            # Goes through all variables in a field component
                            c = 0
                            for var in varnames
                                c += 1
                                b = getproperty(a[j][u], var)
                                # Checks if they are GeoUnit (Value and Unit exist then)
                                if isa(b, GeoUnit)
                                    value = string(getproperty(b, :val))
                                    if typeof(a[j]) <: CompositeRheology
                                        var_counter["$var"] += 1
                                    end
                                    # Gives back LaTex format string for the corresponding variable if it is longer than 2 chars
                                    if length(unidecode("$var")) > 2
                                        latexvar = string("\\" * string(unidecode("$var")))
                                    else
                                        latexvar = string("$var")
                                    end
                                    # Put value, LaTex variable name and creep law pattern in a Dict
                                    counter = var_counter["$var"]
                                    if typeof(a[j][u]) <: DislocationCreep
                                        fds["$var $label $flowadd $i.$u"] = (value, latexvar, "$flowlaw$comporheo)", "$i", "$counter", "$num_rheologies")
                                    end
                                    if typeof(a[j][u]) <: DiffusionCreep
                                        fds["$var $label $flowadd $i.$u"] = (value, latexvar, "$flowlaw$comporheo)", "$i", "$counter", "$num_rheologies")
                                    end
                                    if typeof(a[j][u]) <: LinearViscous
                                        fds["$var $label $flowadd $i.$u"] = (value, latexvar, "$flowlaw$comporheo)", "$i", "$counter", "$num_rheologies")
                                    end
                                end
                                k += 1
                            end
                        end
                        flowlaw *= comporheo * "),"

                    elseif typeof(a[j]) <: Parallel
                        num_parallel = length(a[j].elements)
                        fdsname = "Para "
                        flowadd = flowpara
                        parallelrheo = flowadd * "("
                        for q in 1:num_parallel
                            if typeof(a[j][q]) <: DislocationCreep
                                global Disl += 1
                                flowlawcount += 1
                                flowadd = flowdisl
                                parallelrheo *= flowadd * ","
                                name = fdsname * join(a[j][q].Name)
                                bibinfo_disl = param_info(a[j][q])
                                bib_disl = bibinfo_disl.BibTex_Reference
                                refs["$name"] = (bib_disl, "$flowlawcount", "$i")

                            elseif typeof(a[j][q]) <: DiffusionCreep
                                global Diff += 1
                                flowlawcount += 1
                                flowadd = flowdiff
                                parallelrheo *= flowadd * ","
                                name = fdsname * join(a[j][q].Name)
                                bibinfo_diff = param_info(a[j][q])
                                bib_diff = bibinfo_diff.BibTex_Reference
                                refs["$name"] = (bib_diff, "$flowlawcount", "$i")

                            elseif typeof(a[j][q]) <: LinearViscous
                                global Lin += 1
                                flowlawcount += 1
                                flowadd = flowlin
                                parallelrheo *= flowadd * ","

                            elseif typeof(a[j][q]) <: CompositeRheology
                                compos = getproperty(a[j], :elements)
                                global num_rheologies += length(compos)
                                namepre = fdsname * "Comp "
                                flowadd = flowcomp
                                comporheo = flowadd * "("
                                for u in 1:num_rheologies
                                    if typeof(a[j][q][u]) <: DislocationCreep
                                        global Disl += 1
                                        flowlawcount += 1
                                        flowadd = flowdisl
                                        comporheo *= flowadd * ","
                                        name = namepre * join(a[j][q][u].Name)
                                        bibinfo_disl = param_info(a[j][q][u])
                                        bib_disl = bibinfo_disl.BibTex_Reference
                                        refs["$name"] = (bib_disl, "$flowlawcount", "$i")

                                    elseif typeof(a[j][q][u]) <: DiffusionCreep
                                        global Diff += 1
                                        flowlawcount += 1
                                        flowadd = flowdiff
                                        comporheo *= flowadd * ","
                                        name = namepre * join(a[j][q][u].Name)
                                        bibinfo_diff = param_info(a[j][q][u])
                                        bib_diff = bibinfo_diff.BibTex_Reference
                                        refs["$name"] = (bib_diff, "$flowlawcount", "$i")

                                    elseif typeof(a[j][q][u]) <: LinearViscous
                                        global Lin += 1
                                        flowlawcount += 1
                                        flowadd = flowlin
                                        comporheo *= flowadd * ","
                                    end

                                    varnames = propertynames(a[j][q][u])
                                    # Goes through all variables in a field component
                                    for var in varnames
                                        b = getproperty(a[j][q][u], var)
                                        # Checks if they are GeoUnit (Value and Unit exist then)
                                        if isa(b, GeoUnit)
                                            value = string(getproperty(b, :val))
                                            if typeof(a[j]) <: Parallel
                                                var_counter["$var"] += 1
                                            end
                                            # Gives back LaTex format string for the corresponding variable if it is longer than 2 chars
                                            if length(unidecode("$var")) > 2
                                                latexvar = string("\\" * string(unidecode("$var")))
                                            else
                                                latexvar = string("$var")
                                            end
                                            # Put value, LaTex variable name and creep law pattern in a Dict
                                            counter = var_counter["$var"]
                                            if typeof(a[j][q][u]) <: DislocationCreep
                                                fds["$var $label $flowadd $i.$q.$u"] = (value, latexvar, "$flowlaw$parallelrheo$comporheo)", "$i", "$counter", "$num_rheologies")
                                            end
                                            if typeof(a[j][q][u]) <: DiffusionCreep
                                                fds["$var $label $flowadd $i.$q.$u"] = (value, latexvar, "$flowlaw$parallelrheo$comporheo)", "$i", "$counter", "$num_rheologies")
                                            end
                                            if typeof(a[j][q][u]) <: LinearViscous
                                                fds["$var $label $flowadd $i.$q.$u"] = (value, latexvar, "$flowlaw$parallelrheo$comporheo)", "$i", "$counter", "$num_rheologies")
                                            end
                                        end
                                        k += 1
                                    end

                                end
                                parallelrheo *= comporheo * "),"
                                global num_rheologies += num_parallel - 1
                            end

                            varnames = propertynames(a[j][q])
                            # Goes through all variables in a field component
                            for var in varnames
                                b = getproperty(a[j][q], var)
                                # Checks if they are GeoUnit (Value and Unit exist then)
                                if isa(b, GeoUnit)
                                    value = string(getproperty(b, :val))
                                    if typeof(a[j]) <: CompositeRheology
                                        var_counter["$var"] += 1
                                    end
                                    # Gives back LaTex format string for the corresponding variable if it is longer than 2 chars
                                    if length(unidecode("$var")) > 2
                                        latexvar = string("\\" * string(unidecode("$var")))
                                    else
                                        latexvar = string("$var")
                                    end
                                    # Put value, LaTex variable name and creep law pattern in a Dict
                                    counter = var_counter["$var"]
                                    if typeof(a[j][q]) <: DislocationCreep
                                        fds["$var $label $flowadd $i.$q"] = (value, latexvar, "$flowlaw$parallelrheo)", "$i", "$counter", "$num_rheologies")
                                    end
                                    if typeof(a[j][q]) <: DiffusionCreep
                                        fds["$var $label $flowadd $i.$q"] = (value, latexvar, "$flowlaw$parallelrheo)", "$i", "$counter", "$num_rheologies")
                                    end
                                    if typeof(a[j][q]) <: LinearViscous
                                        fds["$var $label $flowadd $i.$q"] = (value, latexvar, "$flowlaw$parallelrheo)", "$i", "$counter", "$num_rheologies")
                                    end
                                end
                                k += 1
                            end

                        end
                        flowlaw *= parallelrheo * "),"
                    end
                    varnames = propertynames(a[j])
                    # Goes through all variables in a field component
                    if !(typeof(a[j]) <: CompositeRheology || typeof(a[j]) <: Parallel)
                        for var in varnames
                            b = getproperty(a[j], var)
                            # Checks if they are GeoUnit (Value and Unit exist then)
                            if isa(b, GeoUnit)
                                value = string(getproperty(b, :val))
                                if (typeof(a[j]) <: DiffusionCreep || typeof(a[j]) <: DislocationCreep || typeof(a[j]) <: LinearViscous)
                                    var_counter["$var"] += 1
                                    counter = var_counter["$var"]
                                end
                                # Gives back LaTex format string for the corresponding variable if it is longer than 2 chars
                                if length(unidecode("$var")) > 2
                                    latexvar = string("\\" * string(unidecode("$var")))
                                else
                                    latexvar = string("$var")
                                end
                                # Put value, LaTex variable name and creep law pattern in a Dict
                                if typeof(a[j]) <: DislocationCreep
                                    fds["$var $label $i"] = (value, latexvar, "$flowlaw", "$i", "$counter", "$num_rheologies")
                                elseif typeof(a[j]) <: DiffusionCreep
                                    fds["$var $label $i"] = (value, latexvar, "$flowlaw", "$i", "$counter", "$num_rheologies")
                                elseif typeof(a[j]) <: LinearViscous
                                    fds["$var $label $i"] = (value, latexvar, "$flowlaw", "$i", "$counter", "$num_rheologies")
                                else
                                    fds["$var $label $i"] = (value, latexvar, "", "$i", "", "")
                                end
                                
                            end
                            k += 1
                        end
                    end
                end
                # Takes field "Name" and puts it in the first entry of the Tuple, takes Phasecount of all Phases and puts it in the second entry of the Dict
            elseif !isempty(getproperty(s[i], label)) && label == :Name
                phasename = join(getproperty(s[i], :Name))
                fds["$label $i"] = (phasename, "$phasecount", "$i", "", "", "")
            end
        end
    end
    return fds, refs
end

"""
Dict2LatexTable() writes a .tex file with all parameters from the Phase2Dict() output in a LaTeX table. rdigits will round numbers with more decimals than rdigits
including numbers of 10 to power of n, n being an Integer, for representation purposes. For the exact numbers use the original impemented numbers from the creeplaws of the dict in src/CreepLaw/Data/DiffusionCreep.jl
or src/CreepLaw/Data/DislocationCreep.jl.


"""

function Dict2LatexTable(d::Dict, refs::Dict; filename="ParameterTable", rdigits=4)
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
        "E" => "Activation energy (J/mol)", # can also be Elastic Young's modulus, maybe change symbol for that?
        "R" => "Gas constant (J/mol/K)",
        "G" => "Shear modulus (Pa)",
        "\\nu" => "Poisson ratio (-)",
        "K" => "Bulk modulus (Pa)",
        "Kb" => "Elastic bulk modulus (Pa)",
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
    Table *= "\\usepackage[utf8]{inputenc}\n"
    Table *= "\\usepackage{booktabs}\n"
    Table *= "\\usepackage{graphicx}\n"
    Table *= "\\usepackage{multirow}\n"
    Table *= "\\usepackage[round, comma, sort, ]{natbib}\n"
    Table *= "\\bibliographystyle{abbrvnat}\n"
    Table *= "\\setcitestyle{authoryear,open={(},close={)}} %Citation-related commands\n"
    Table *= "\\begin{document}\n"
    Table *= "\\begin{table}[hbt]\n"
    Table *= "\\centering\n"
    Table *= "\\resizebox{\\columnwidth}{!}{\n"
    Table *= "\\begin{tabular}{ " * "l " * "c "^(parse(Int64, d["Name 1"][2]) + 2) * "}\n \\\\"
    Table *= "\\toprule[1pt]\n \\\\"
    Table *= "\\midrule[0.3pt]\n"
    Table *= " & "

    # Boris nach Layout für Tabellen mit CompositeRheologies fragen. 
    # Wie soll z.B. der gleiche Parameter für allerdings verschiedene CreepLaws in einer Phase gehandhabt werden (zwei Mal Activation Energy z.B.)? 

    # Creates table headers and in-text citations
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
    Table *= " \\\\\n"
    Table *= "\\midrule \n"
    Table *= "% Table body\n"

    # Get vector with all unique symbols without phasenames
    for key in 1:length(dictpairs)
        dictpairs_key = dictpairs[key].first
        if occursin("Name", dictpairs_key)
            continue
        # Checks if symbol is from field CompositeRheology and if the symbol occurs more often than once
        elseif (occursin("R", dictpairs_key[1:3]))
            push!(symbs, dictpairs[key].second[2])
        elseif maximum(occursin.(["CompositeRheology", "CreepLaws"], dictpairs_key)) && !(occursin("R", dictpairs_key[1:3]))
            number = parse(Int64, dictpairs[key].second[5])
            push!(symbs, dictpairs[key].second[2] * "$number")
        else
            push!(symbs, dictpairs[key].second[2])
        end
    end
    symbs = unique(sort(symbs))

    # Creates columnwise output for all parameters of the input phase
    for symbol in symbs
        # Sets parametername and variable
        if maximum(endswith.(symbol,["0","1","2","3","4","5","6","7","8","9"])) || (!occursin("\\", symbol) && length(symbol) > 1) 
            if occursin("0", symbol)
                Table *= " " * string(desc[symbol]) *  " & " * "\$" * symbol[1:end-1] * "_" * symbol[end] *"\$"
            else
                Table *= " " * string(desc[symbol[1:end-1]]) *  " & " * "\$" * symbol[1:end-1] * "_" * symbol[end] *"\$"
            end
        else
            Table *= " " * string(desc[symbol]) *  " & " * "\$" * symbol *"\$"
        end
        # Iterates over all phases
        for j in 1:parse(Int64, d["Name 1"][2])
            hit = 0
            i_dictpairs = 1
            # Iterates over all pairs
            for element in dictpairs
                # Checks if symbol matches the symbol in the pair and if phase matches phase of the pair
                if maximum(endswith.(symbol,["1","2","3","4","5","6","7","8","9"]))
                    symbol = symbol[1:end-1]
                end
                #--------------------------------------------- WORK IN PROGRESS FOR COMPOSITE RHEOLOGIES-------------------------------------------------------------
                if symbol == dictpairs[i_dictpairs].second[2] && j == parse(Int64, dictpairs[i_dictpairs].second[4]) && hit == 0
                    # put in the matched parameter value
                    dig, num, expo = detachFloatfromExponent(dictpairs[i_dictpairs].second[1])
                    if  dig <= rdigits && expo != "1" 
                        Table *= " & \$" * "$num \\times 10^{$expo} \$"
                    elseif dig <= rdigits && expo == "1"
                        Table *= " & \$" * "$num\$"
                    elseif dig > rdigits && expo != "1"
                        Table *= " & \$" * string(round(parse(Float64, num); digits=rdigits)) * " \\times 10^{$expo}\$"
                    elseif dig > rdigits && expo == "1"
                        Table *= " & \$" * string(round(parse(Float64, num); digits=rdigits)) * "\$"
                    end
                    hit += 1
                    deleteat!(dictpairs, i_dictpairs)
                end
                i_dictpairs += 1
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
    Table *= InTextRef * "\\\\\n"

    # finishes latex table, closes all open formats and writes References beneath table
    Table *= "\\midrule[0.3pt] \\\\\n"
    Table *= "\\bottomrule[1pt]\n"
    Table *= "\\end{tabular}\n"
    Table *= "}\n"
    Table *= "\\caption{Parameter table}\n"
    Table *= "\\label{tab:para_table}\n"
    Table *= "\\end{table}\n"
    Table *= "\\bibliography{$filename}\n"
    Table *= "\\end{document}\n"

    # Writes BibTex sources in to .bib file and Table string into .tex file
    return write("$filename.bib", References), write("$filename.tex", Table)
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
                            mdvar = "$var"
                            # Put value, LaTex name and units in a Dict
                            fds["$var $label $i"] = (value, mdvar, "", "$i") 
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
Dict2MarkdownTable() writes a .md file with all parameters from the Phase2DictMd() output in a Markdown table. rdigits will round numbers with more decimals than rdigits
including numbers of 10 to power of n, n being an Integer, for representation purposes. For the exact numbers use the original impemented numbers from the creeplaws of the dict in src/CreepLaw/Data/DiffusionCreep.jl
or src/CreepLaw/Data/DislocationCreep.jl.

"""

function Dict2MarkdownTable(d::Dict; filename="ParameterTable", rdigits=4)
    dictkeys = keys(d)
    symbs = []

    # Creates vectors of type Pairs (can be iterated over, sorted, etc.)
    dictpairs = sort(collect(pairs(d)))


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
        "Kb"=>"Elastic bulk modulus (Pa)", 
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

    return write("$filename.md", Table)
end


"""
ParameterTable() creates a table with all parameters saved in the Phase struct. It lets you choose between "latex" and "markdown" as table formats with LaTeX as default.
Creates a filename.tex or filename.md file as output. If "latex" is chosen a "Reference.bib" file will automatically be produced with all your references.
"""
function ParameterTable(
    Phase;
    filename="ParameterTable",
    format="latex",
    rdigits=4
    )

    format = lowercase(format)

    if (format == "latex") || (format == "tex")
        d, ref = Phase2Dict(Phase)
        Dict2LatexTable(d, ref, filename=filename, rdigits=rdigits)
    elseif (format == "markdown") || (format == "md")
        d = Phase2DictMd(Phase)
        Dict2MarkdownTable(d, filename=filename, rdigits=rdigits)
    end
    return nothing
end

end
