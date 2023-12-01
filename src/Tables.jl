module Tables

using Unidecode
using GeoParams: AbstractMaterialParam, param_info, LinearViscous, PowerlawViscous, DislocationCreep, DiffusionCreep, GrainBoundarySliding, PeierlsCreep, NonLinearPeierlsCreep, CompositeRheology, Parallel, make_tuple
using ..Units
using ..MaterialParameters: MaterialParamsInfo

export detachFloatfromExponent, check_flowadd, checkRheologyType, Phase2Dict, Dict2LatexTable, Phase2DictMd, Dict2MarkdownTable, ParameterTable


""" 
detachFloatfromExponent() returns the number of float decimals after the comma as Integer, the float number without the exponent as string 
and the exponent after the "e" as string.
The argument output returns "1" for "ex" if the input number has no exponent.

"""
function detachFloatfromExponent(str::String)
    s = lowercase(str)
    if 'e' in s
        s, ex = string.(split(s, "e"))
    else
        ex = "1"
    end

    dig = something(findfirst('.', reverse(s)), 1) - 1
    
    return dig, s, ex

end

@inline check_flowadd(a::DislocationCreep) = "DislCreep"
@inline check_flowadd(a::DiffusionCreep) = "DiffCreep"
@inline check_flowadd(a::GrainBoundarySliding) = "GBS"
@inline check_flowadd(a::PeierlsCreep) = "PeiCreep"
@inline check_flowadd(a::NonLinearPeierlsCreep) = "NLP"
@inline check_flowadd(a::LinearViscous) = "LinVisc"
@inline check_flowadd(::T) where T = throw("Rheology not supported.")
#Noch nachdenken wie jetzt nach Typ rheo zu Disl, Diff, etc. dazu addiert wird

function checkRheologyType(a, rheo, flowlaw, flowlawcount, refs, i; pre="")
    rheo += 1
    flowadd = check_flowadd(a)
    flowlawcount += 1
    flowlaw *= flowadd * ","
    name = pre * join(a.Name)
    bibinfo = param_info(a)
    bib = bibinfo.BibTex_Reference
    refs["$name"] = (bib, "$flowlawcount", "$i")
    return rheo, flowlaw, flowlawcount, refs
end

"""
Phase2Dict() puts all parameters of a phase s in a dict.
"""
function Phase2Dict(s)
    # Dict has Key with Fieldname and Value with Tuple(value, symbol, flowlaw keyword, number of phases, (Disl+Diff) or Lin, Disl or Diff)
    fds = Dict{String,Tuple{String,String,String,String,String,String}}()
    refs = Dict{String,Tuple{String,String,String}}()
    # Makes sure phase is a tuple
    s = make_tuple(s)
    phasecount = length(s)
    k = 1
    flowlawcount = 0
    # Checks all Phases 
    for i in 1:phasecount
        fieldnames = propertynames(s[i])
        Diff = 0
        Disl = 0
        Lin = 0
        GBS = 0
        Pei = 0
        NLP = 0
        for label in fieldnames # labels are symbols
            if !isempty(getproperty(s[i], label)) && label != :Name
                # Goes through all components of all fields
                for j in 1:length(getproperty(s[i], label))
                    a = getproperty(s[i], label)
                    flowlaw  = ""
                    flowdisl = "DislCreep"
                    flowdiff = "DiffCreep"
                    flowlin  = "LinVisc"
                    flowgbs  = "GBS"
                    flowpei  = "PeiCreep"
                    flownlp  = "NLP"
                    flowcomp = "CompoRheo"
                    flowpara = "Parallel"
                    a_j = a[j]
                    # Checks what type the CreepLaw or CompositeRheology field has 
                    if typeof(a_j) <: DislocationCreep
                        Disl, flowlaw, flowlawcount, refs = checkRheologyType(a_j, Disl, flowlaw, flowlawcount, refs, i)

                    elseif typeof(a_j) <: DiffusionCreep
                        Diff, flowlaw, flowlawcount, refs = checkRheologyType(a_j, Diff, flowlaw, flowlawcount, refs, i)

                    elseif typeof(a_j) <: GrainBoundarySliding
                        GBS, flowlaw, flowlawcount, refs = checkRheologyType(a_j, GBS, flowlaw, flowlawcount, refs, i)

                    elseif typeof(a_j) <: PeierlsCreep
                        Pei, flowlaw, flowlawcount, refs = checkRheologyType(a_j, Pei, flowlaw, flowlawcount, refs, i)

                    elseif typeof(a_j) <: NonLinearPeierlsCreep
                        NLP, flowlaw, flowlawcount, refs = checkRheologyType(a_j, NLP, flowlaw, flowlawcount, refs, i)

                    elseif typeof(a_j) <: LinearViscous
                        Lin += 1
                        flowlawcount += 1
                        flowadd = flowlin
                        flowlaw *= flowadd

                    elseif typeof(a_j) <: CompositeRheology
                        compos = getproperty(a_j, :elements)
                        num_rheologies = length(compos)
                        fdsname = "Comp "
                        flowadd = flowcomp
                        comporheo = flowadd * "("
                        for u in 1:num_rheologies
                            a_ju = a_j[u]
                            if typeof(a_ju) <: Parallel
                                num_parallel = length(a_ju.elements)
                                flowadd = flowpara
                                parallelrheo = flowadd * "("
                                namepre = fdsname * "Para "
                                for v in 1:num_parallel
                                    a_juv = a_ju[v]
                                    if typeof(a_juv) <: DislocationCreep
                                        Disl, parallelrheo, flowlawcount, refs = checkRheologyType(a_juv, Disl, parallelrheo, flowlawcount, refs, i; pre=namepre)
                                        
                                    elseif typeof(a_juv) <: DiffusionCreep
                                        Diff, parallelrheo, flowlawcount, refs = checkRheologyType(a_juv, Diff, parallelrheo, flowlawcount, refs, i; pre=namepre)

                                    elseif typeof(a_juv) <: GrainBoundarySliding
                                        GBS, parallelrheo, flowlawcount, refs = checkRheologyType(a_juv, GBS, parallelrheo, flowlawcount, refs, i; pre=namepre)
                
                                    elseif typeof(a_juv) <: PeierlsCreep
                                        Pei, parallelrheo, flowlawcount, refs = checkRheologyType(a_juv, Pei, parallelrheo, flowlawcount, refs, i; pre=namepre)
                
                                    elseif typeof(a_juv) <: NonLinearPeierlsCreep
                                        NLP, parallelrheo, flowlawcount, refs = checkRheologyType(a_juv, NLP, parallelrheo, flowlawcount, refs, i; pre=namepre)
                                        
                                    elseif typeof(a_juv) <: LinearViscous
                                        Lin += 1
                                        flowlawcount += 1
                                        flowadd = flowlin
                                        parallelrheo *= flowadd * ","
                                        
                                    end

                                    varnames = propertynames(a_juv)
                                    # Goes through all variables in a field component
                                    for var in varnames
                                        b = getproperty(a_juv, var)
                                        # Checks if they are GeoUnit (Value and Unit exist then)
                                        if isa(b, GeoUnit)
                                            value = string(getproperty(b, :val))
                                            # Gives back LaTex format string for the corresponding variable if it is longer than 2 chars
                                            if length(unidecode("$var")) > 2
                                                latexvar = string("\\" * string(unidecode("$var")))
                                            else
                                                latexvar = string("$var")
                                            end
                                            # Put value, LaTex variable name, creep law pattern, phase id, current disl, diff or linvisc count and corresponding string for equation in a Dict
                                            if typeof(a_juv) <: DislocationCreep
                                                fds["$var $label $flowadd $i.$u.$v"] = (value, latexvar, "$flowlaw$comporheo$parallelrheo)", "$i", "$Disl", "$flowdisl")
                                            end
                                            if typeof(a_juv) <: DiffusionCreep
                                                fds["$var $label $flowadd $i.$u.$v"] = (value, latexvar, "$flowlaw$comporheo$parallelrheo)", "$i", "$Diff", "$flowdiff")
                                            end
                                            if typeof(a_juv) <: GrainBoundarySliding
                                                fds["$var $label $flowadd $i.$u.$v"] = (value, latexvar, "$flowlaw$comporheo$parallelrheo)", "$i", "$GBS", "$flowgbs")
                                            end
                                            if typeof(a_juv) <: PeierlsCreep
                                                fds["$var $label $flowadd $i.$u.$v"] = (value, latexvar, "$flowlaw$comporheo$parallelrheo)", "$i", "$Pei", "$flowpei")
                                            end
                                            if typeof(a_juv) <: NonLinearPeierlsCreep
                                                fds["$var $label $flowadd $i.$u.$v"] = (value, latexvar, "$flowlaw$comporheo$parallelrheo)", "$i", "$NLP", "$flownlp")
                                            end
                                            if typeof(a_juv) <: LinearViscous
                                                fds["$var $label $flowadd $i.$u.$v"] = (value, latexvar, "$flowlaw$comporheo$parallelrheo)", "$i", "$Lin", "$flowlin")
                                            end
                                        end
                                        k += 1
                                    end
                                end
                                comporheo *= parallelrheo * "),"

                            elseif typeof(a_ju) <: DislocationCreep
                                Disl, comporheo, flowlawcount, refs = checkRheologyType(a_ju, Disl, comporheo, flowlawcount, refs, i; pre=fdsname)

                            elseif typeof(a_ju) <: DiffusionCreep
                                Diff, comporheo, flowlawcount, refs = checkRheologyType(a_ju, Diff, comporheo, flowlawcount, refs, i; pre=fdsname)

                            elseif typeof(a_ju) <: GrainBoundarySliding
                                GBS, comporheo, flowlawcount, refs = checkRheologyType(a_ju, GBS, comporheo, flowlawcount, refs, i; pre=fdsname)

                            elseif typeof(a_ju) <: PeierlsCreep
                                Pei, comporheo, flowlawcount, refs = checkRheologyType(a_ju, Pei, comporheo, flowlawcount, refs, i; pre=fdsname)

                            elseif typeof(a_ju) <: NonLinearPeierlsCreep
                                NLP, comporheo, flowlawcount, refs = checkRheologyType(a_ju, NLP, comporheo, flowlawcount, refs, i; pre=fdsname)

                            elseif typeof(a_ju) <: LinearViscous
                                Lin += 1
                                flowlawcount += 1
                                flowadd = flowlin
                                comporheo *= flowadd * ","
                            end

                            varnames = propertynames(a_ju)
                            # Goes through all variables in a field component
                            c = 0
                            for var in varnames
                                c += 1
                                b = getproperty(a_ju, var)
                                # Checks if they are GeoUnit (Value and Unit exist then)
                                if isa(b, GeoUnit)
                                    value = string(getproperty(b, :val))
                                    # Gives back LaTex format string for the corresponding variable if it is longer than 2 chars
                                    if length(unidecode("$var")) > 2
                                        latexvar = string("\\" * string(unidecode("$var")))
                                    else
                                        latexvar = string("$var")
                                    end
                                    # Put value, LaTex variable name, creep law pattern, phase id, current disl, diff or linvisc count and corresponding string for equation in a Dict
                                    if typeof(a_ju) <: DislocationCreep
                                        fds["$var $label $flowadd $i.$u"] = (value, latexvar, "$flowlaw$comporheo)", "$i", "$Disl", "$flowdisl")
                                    end
                                    if typeof(a_ju) <: DiffusionCreep
                                        fds["$var $label $flowadd $i.$u"] = (value, latexvar, "$flowlaw$comporheo)", "$i", "$Diff", "$flowdiff")
                                    end
                                    if typeof(a_ju) <: GrainBoundarySliding
                                        fds["$var $label $flowadd $i.$u"] = (value, latexvar, "$flowlaw$comporheo)", "$i", "$GBS", "$flowgbs")
                                    end
                                    if typeof(a_ju) <: PeierlsCreep
                                        fds["$var $label $flowadd $i.$u"] = (value, latexvar, "$flowlaw$comporheo)", "$i", "$Pei", "$flowpei")
                                    end
                                    if typeof(a_ju) <: NonLinearPeierlsCreep
                                        fds["$var $label $flowadd $i.$u"] = (value, latexvar, "$flowlaw$comporheo)", "$i", "$NLP", "$flownlp")
                                    end
                                    if typeof(a_ju) <: LinearViscous
                                        fds["$var $label $flowadd $i.$u"] = (value, latexvar, "$flowlaw$comporheo)", "$i", "$Lin", "$flowlin")
                                    end
                                end
                                k += 1
                            end
                        end
                        flowlaw *= comporheo * "),"

                    elseif typeof(a_j) <: Parallel
                        num_parallel = length(a_j.elements)
                        fdsname = "Para "
                        flowadd = flowpara
                        parallelrheo = flowadd * "("
                        for q in 1:num_parallel
                            a_jq = a_j[q]
                            if typeof(a_jq) <: DislocationCreep
                                Disl, parallelrheo, flowlawcount, refs = checkRheologyType(a_jq, Disl, parallelrheo, flowlawcount, refs, i; pre=fdsname)

                            elseif typeof(a_jq) <: DiffusionCreep
                                Diff, parallelrheo, flowlawcount, refs = checkRheologyType(a_jq, Diff, parallelrheo, flowlawcount, refs, i; pre=fdsname)

                            elseif typeof(a_jq) <: GrainBoundarySliding
                                GBS, parallelrheo, flowlawcount, refs = checkRheologyType(a_jq, GBS, parallelrheo, flowlawcount, refs, i; pre=fdsname)

                            elseif typeof(a_jq) <: PeierlsCreep
                                Pei, parallelrheo, flowlawcount, refs = checkRheologyType(a_jq, Pei, parallelrheo, flowlawcount, refs, i; pre=fdsname)

                            elseif typeof(a_jq) <: NonLinearPeierlsCreep
                                NLP, parallelrheo, flowlawcount, refs = checkRheologyType(a_jq, NLP, parallelrheo, flowlawcount, refs, i; pre=fdsname)

                            elseif typeof(a_jq) <: LinearViscous
                                Lin += 1
                                flowlawcount += 1
                                flowadd = flowlin
                                parallelrheo *= flowadd * ","

                            elseif typeof(a_jq) <: CompositeRheology
                                compos = getproperty(a_j, :elements)
                                num_rheologies = length(compos)
                                namepre = fdsname * "Comp "
                                flowadd = flowcomp
                                comporheo = flowadd * "("
                                for u in 1:num_rheologies
                                    a_jqu = a_jq[u]
                                    if typeof(a_jqu) <: DislocationCreep
                                        Disl, comporheo, flowlawcount, refs = checkRheologyType(a_jqu, Disl, comporheo, flowlawcount, refs, i; pre=namepre)

                                    elseif typeof(a_jqu) <: DiffusionCreep
                                        Diff, comporheo, flowlawcount, refs = checkRheologyType(a_jqu, Diff, comporheo, flowlawcount, refs, i; pre=namepre)

                                    elseif typeof(a_jqu) <: GrainBoundarySliding
                                        GBS, comporheo, flowlawcount, refs = checkRheologyType(a_jqu, GBS, comporheo, flowlawcount, refs, i; pre=namepre)

                                    elseif typeof(a_jqu) <: PeierlsCreep
                                        Pei, comporheo, flowlawcount, refs = checkRheologyType(a_jqu, Pei, comporheo, flowlawcount, refs, i; pre=namepre)

                                    elseif typeof(a_jqu) <: NonLinearPeierlsCreep
                                        NLP, flowlaw, flowlawcount, refs = checkRheologyType(a_jqu, NLP, flowlaw, flowlawcount, refs, i; pre=namepre)

                                    elseif typeof(a_jqu) <: LinearViscous
                                        Lin += 1
                                        flowlawcount += 1
                                        flowadd = flowlin
                                        comporheo *= flowadd * ","
                                    end

                                    varnames = propertynames(a_jqu)
                                    # Goes through all variables in a field component
                                    for var in varnames
                                        b = getproperty(a_jqu, var)
                                        # Checks if they are GeoUnit (Value and Unit exist then)
                                        if isa(b, GeoUnit)
                                            value = string(getproperty(b, :val))
                                            # Gives back LaTex format string for the corresponding variable if it is longer than 2 chars
                                            if length(unidecode("$var")) > 2
                                                latexvar = string("\\" * string(unidecode("$var")))
                                            else
                                                latexvar = string("$var")
                                            end
                                            # Put value, LaTex variable name, creep law pattern, phase id, current disl, diff or linvisc count and corresponding string for equation in a Dict
                                            if typeof(a_jqu) <: DislocationCreep
                                                fds["$var $label $flowadd $i.$q.$u"] = (value, latexvar, "$flowlaw$parallelrheo$comporheo)", "$i", "$Disl", "$flowdisl")
                                            end
                                            if typeof(a_jqu) <: DiffusionCreep
                                                fds["$var $label $flowadd $i.$q.$u"] = (value, latexvar, "$flowlaw$parallelrheo$comporheo)", "$i", "$Diff", "$flowdiff")
                                            end
                                            if typeof(a_jqu) <: GrainBoundarySliding
                                                fds["$var $label $flowadd $i.$q.$u"] = (value, latexvar, "$flowlaw$parallelrheo$comporheo)", "$i", "$GBS", "$flowgbs")
                                            end
                                            if typeof(a_jqu) <: PeierlsCreep
                                                fds["$var $label $flowadd $i.$q.$u"] = (value, latexvar, "$flowlaw$parallelrheo$comporheo)", "$i", "$Pei", "$flowpei")
                                            end
                                            if typeof(a_jqu) <: NonLinearPeierlsCreep
                                                fds["$var $label $flowadd $i.$q.$u"] = (value, latexvar, "$flowlaw$parallelrheo$comporheo)", "$i", "$NLP", "$flownlp")
                                            end
                                            if typeof(a_jqu) <: LinearViscous
                                                fds["$var $label $flowadd $i.$q.$u"] = (value, latexvar, "$flowlaw$parallelrheo$comporheo)", "$i", "$Lin", "$flowlin")
                                            end
                                        end
                                        k += 1
                                    end

                                end
                                parallelrheo *= comporheo * "),"
                            end

                            varnames = propertynames(a_jq)
                            # Goes through all variables in a field component
                            for var in varnames
                                b = getproperty(a_jq, var)
                                # Checks if they are GeoUnit (Value and Unit exist then)
                                if isa(b, GeoUnit)
                                    value = string(getproperty(b, :val))
                                    # Gives back LaTex format string for the corresponding variable if it is longer than 2 chars
                                    if length(unidecode("$var")) > 2
                                        latexvar = string("\\" * string(unidecode("$var")))
                                    else
                                        latexvar = string("$var")
                                    end
                                    # Put value, LaTex variable name, creep law pattern, phase id, current disl, diff or linvisc count and corresponding string for equation in a Dict
                                    if typeof(a_jq) <: DislocationCreep
                                        fds["$var $label $flowadd $i.$q"] = (value, latexvar, "$flowlaw$parallelrheo)", "$i","$Disl", "$flowdisl")
                                    end
                                    if typeof(a_jq) <: DiffusionCreep
                                        fds["$var $label $flowadd $i.$q"] = (value, latexvar, "$flowlaw$parallelrheo)", "$i", "$Diff", "$flowdiff")
                                    end
                                    if typeof(a_jq) <: GrainBoundarySliding
                                        fds["$var $label $flowadd $i.$q"] = (value, latexvar, "$flowlaw$parallelrheo)", "$i", "$GBS", "$flowgbs")
                                    end
                                    if typeof(a_jq) <: PeierlsCreep
                                        fds["$var $label $flowadd $i.$q"] = (value, latexvar, "$flowlaw$parallelrheo)", "$i", "$Pei", "$flowpei")
                                    end
                                    if typeof(a_jq) <: NonLinearPeierlsCreep
                                        fds["$var $label $flowadd $i.$q"] = (value, latexvar, "$flowlaw$parallelrheo)", "$i", "$NLP", "$flownlp")
                                    end
                                    if typeof(a_jq) <: LinearViscous
                                        fds["$var $label $flowadd $i.$q"] = (value, latexvar, "$flowlaw$parallelrheo)", "$i", "$Lin", "$flowlin")
                                    end
                                end
                                k += 1
                            end

                        end
                        flowlaw *= parallelrheo * "),"
                    end
                    varnames = propertynames(a_j)
                    # Goes through all variables in a field component
                    if !(typeof(a_j) <: CompositeRheology || typeof(a_j) <: Parallel)
                        for var in varnames
                            b = getproperty(a_j, var)
                            # Checks if they are GeoUnit (Value and Unit exist then)
                            if isa(b, GeoUnit)
                                value = string(getproperty(b, :val))
                                # Gives back LaTex format string for the corresponding variable if it is longer than 2 chars
                                if length(unidecode("$var")) > 2
                                    latexvar = string("\\" * string(unidecode("$var")))
                                else
                                    latexvar = string("$var")
                                end
                                # Put value, LaTex variable name, creep law pattern, phase id, current disl, diff or linvisc count and corresponding string for equation in a Dict
                                if typeof(a_j) <: DislocationCreep
                                    fds["$var $label $i"] = (value, latexvar, "$flowlaw", "$i", "$Disl", "$flowdisl")
                                elseif typeof(a_j) <: DiffusionCreep
                                    fds["$var $label $i"] = (value, latexvar, "$flowlaw", "$i", "$Diff", "$flowdiff")
                                elseif typeof(a_j) <: GrainBoundarySliding
                                    fds["$var $label $i"] = (value, latexvar, "$flowlaw", "$i", "$GBS", "$flowgbs")
                                elseif typeof(a_j) <: PeierlsCreep
                                    fds["$var $label $i"] = (value, latexvar, "$flowlaw", "$i", "$Pei", "$flowpei")
                                elseif typeof(a_j) <: NonLinearPeierlsCreep
                                    fds["$var $label $i"] = (value, latexvar, "$flowlaw", "$i", "$NLP", "$flownlp")
                                elseif typeof(a_j) <: LinearViscous
                                    fds["$var $label $i"] = (value, latexvar, "$flowlaw", "$i", "$Lin", "$flowlin")
                                # If not a CreepLaw or CompositeRheology
                                else
                                    fds["$var $label $i"] = (value, latexvar, "", "$i", "", "")
                                end
                                
                            end
                            k += 1
                        end
                    end
                end
                # Takes field "Name", puts in Phasename, maximum phase count and current phase number in the Dict
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
including numbers of 10 to power of n, n being an Integer, for representation purposes. For the exact numbers use the original implemented numbers from the creeplaws of the dict in src/CreepLaw/Data/DiffusionCreep.jl
or src/CreepLaw/Data/DislocationCreep.jl.


"""
function Dict2LatexTable(d::Dict, refs::Dict; filename="ParameterTable", rdigits=4)
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
        "o" => "Stres parametrization exponent *(-)*",    # in PeierlsCreep
        "q" => "Stress relation exponent *(-)*",          # in PeierlsCreep
        "A_diff" => "Coefficient *(Pa^-n-r^ m^p^/s)*",  # DiffusionCreep
        "A_disl" => "Coefficient *(Pa^-n^/s)*",         # DislocationCreep
        "A_gbs" => "Coefficient *(Pa^-n^ m^p^/s)*",          # GrainBoundarySliding
        "A_pei" => "Coefficient *(1/-s)*",          # PeierlsCreep
        "A_nlp" => "Coefficient *(Pa^-n^/s)*",          # NonLinearPeierlsCreep
        "Tau_p" => "Peierls stress *(Pa)*",
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

    # Creates table headers and in-text citations
    counter = 1
    for i in 1:parse(Int64, d["Name 1"][2])
        Table *= " & " * d["Name $i"][1]
        for j in 1:length(refs)
            if parse(Int64, refpair[j].second[3]) == i
                if length(refpair[j].second[1]) != 0 
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
                else
                    continue
                end

                if counter != length(refs)
                    InTextRef *= ", "
                end
                counter += 1
            end
        end
    end

    # Latex formatting and comment
    Table *= " \\\\\n"
    Table *= "\\midrule \n"
    Table *= "% Table body\n"

    # Get vector with all unique symbols without phasenames
    for key in 1:length(dictpairs)
        dictpairs_key = dictpairs[key].first
        if occursin("Name", dictpairs_key)
            continue
        # Checks if symbol is R or in generel from field CompositeRheology and if the symbol occurs more often than once
        elseif (occursin("R", dictpairs_key[1:3]))
            push!(symbs, dictpairs[key].second[2])
        elseif maximum(occursin.(["CompositeRheology", "CreepLaws"], dictpairs_key)) && !(occursin("R", dictpairs_key[1:3]))
            if occursin("Disl", dictpairs[key].second[6])
                number = parse(Int64, dictpairs[key].second[5])
                push!(symbs, "A_disl" * "^$number")
            elseif occursin("Diff", dictpairs[key].second[6])
                number = parse(Int64, dictpairs[key].second[5])
                push!(symbs, "A_diff" * "^$number")
            else
                number = parse(Int64, dictpairs[key].second[5])
                push!(symbs, dictpairs[key].second[2] * "$number")
            end
        else
            push!(symbs, dictpairs[key].second[2])
        end
    end
    symbs = unique(sort(symbs))

    # Creates columnwise output for all parameters of the input phase
    for symbol in symbs
        # Sets parametername and variable (symbol)
        # Order of signs in symbol: 1. "_", if already "_" in there -> 2. "^"
        # If no underscore in symbol
        if (maximum(endswith(symbol, string(i),) for i in 0:9) && !occursin("_", symbol)) || (!occursin("\\", symbol) && length(symbol) > 1 && !occursin("_", symbol)) 
            # If "0" in symbol (0 in key for dict)
            if occursin("0", symbol)
                Table *= " " * string(desc[symbol]) *  " & " * "\$" * symbol[1:end-1] * "_" * symbol[end] *"\$"
            # If 1 to 9 in symbol (1:9 not in key for dict)
            else
                Table *= " " * string(desc[symbol[1:end-1]]) *  " & " * "\$" * symbol[1:end-1] * "_" * symbol[end] *"\$"
            end
        # If "_" AND "^" in symbol 
        elseif occursin("_", symbol) && occursin("^", symbol)
            occ_under = findfirst('_', symbol)
            occ_roof = findfirst('^', symbol)
            Table *= " " * string(desc[symbol[1:occ_roof-1]]) *  " & " * "\$" * symbol[1:occ_under] * "{" * symbol[occ_under+1:occ_roof-1] * "}" * symbol[occ_roof] * "{" * symbol[occ_roof+1:end] * "}" * "\$"
        # If only "_" in symbol
        elseif occursin("_", symbol)
            occ_under = findfirst('_', symbol)
            Table *= " " * string(desc[symbol]) *  " & " * "\$" * symbol[1:occ_under] * "{" * symbol[occ_under+1:end] * "}" * "\$"
        # If neither a number nor one of the above signs in symbol, so just a normal symbol
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
                if maximum(endswith(symbol, string(i),) for i in 1:9) && !occursin("A_", symbol)
                    symbol = symbol[1:end-1]
                # Checks if "A_" is in symbol
                elseif occursin("A_", string(symbol)) && ((occursin("diff", symbol) && (occursin("Diff", dictpairs[i_dictpairs].second[6]))) || (occursin("disl", symbol) && occursin("Disl", dictpairs[i_dictpairs].second[6])))
                    symbol_A = string(symbol[1])
                    if symbol_A == dictpairs[i_dictpairs].second[2] && j == parse(Int64, dictpairs[i_dictpairs].second[4]) && hit == 0
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
                end
                # If "A_" is not in symbol
                if symbol == dictpairs[i_dictpairs].second[2] && j == parse(Int64, dictpairs[i_dictpairs].second[4]) && hit == 0 && !(occursin("A_", string(symbol)))
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
    dictpairs = sort(collect(pairs(d)))
    DislCreep = 0
    DiffCreep = 0
    LinVisc = 0
    for i in 1:length(dictpairs)
        if dictpairs[i].second[6] == "DislCreep" && DislCreep == 0

            Table *=
                "\\rule[-5pt]{-3pt}{20pt} Dislocation Creep: & " *
                "\\multicolumn{4}{l}{\$ \\dot{\\gamma} = A \\tau^n f_{H2O}^r \\exp(-\\frac{E+PV}{RT}) \$}\n"
            Table *= " \\\\\n"
            DislCreep += 1
        end
        if dictpairs[i].second[6] == "DiffCreep" && DiffCreep == 0
            Table *=
                "\\rule[-5pt]{-3pt}{20pt} Diffusion Creep: & " *
                "\\multicolumn{4}{l}{\$ \\dot{\\gamma} = A \\tau^n d^p f_{H2O}^r \\exp(-\\frac{E+PV}{RT}) \$}\n"
            Table *= " \\\\\n"
            DiffCreep += 1
        end
        if dictpairs[i].second[6] == "LinVisc" && LinVisc == 0
            Table *=
                "\\rule[-5pt]{-3pt}{20pt} Linear viscosity: & " *
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
    Table *= "\\multicolumn{5}{l}{"
    Table *= InTextRef * "}\\\\\n"

    # Finishes latex table, closes all open formats and writes References beneath table
    Table *= "\\midrule[0.3pt] \\\\\n"
    Table *= "\\bottomrule[1pt]\n"
    Table *= "\\end{tabular}\n"
    Table *= "}\n"
    Table *= "\\caption{Parameter table}\n"
    Table *= "\\label{tab:para_table}\n"
    Table *= "\\end{table}\n"
    Table *= "\\bibliography{References}\n"
    Table *= "\\end{document}\n"

    # Writes BibTex sources into .bib file and Table string into .tex file
    return write("References.bib", References), write("$filename.tex", Table)
end


"""
Phase2DictMd() puts all parameters of a phase in a dict.

"""
function Phase2DictMd(s)

    # Dict has Key with Fieldname and Value with Tuple(value, symbol, creep law pattern, phase number, current disl, diff or linvisc count, corresponding string)
    fds = Dict{String,Tuple{String,String,String,String,String,String}}()
    refs = Dict{String,Tuple{String,String,String}}()
    s = make_tuple(s)
    phasecount = length(s)
    k = 1
    flowlawcount = 0
    # Checks all Phases 
    for i in 1:phasecount
        fieldnames = propertynames(s[i])
        Diff = 0
        Disl = 0
        Lin = 0
        for label in fieldnames
            if !isempty(getproperty(s[i], label)) && label != :Name
                # Goes through all components of all fields
                for j in 1:length(getproperty(s[i], label))
                    a = getproperty(s[i], label)
                    flowlaw  = ""
                    flowdisl = "DislCreep"
                    flowdiff = "DiffCreep"
                    flowlin  = "LinVisc"
                    flowgbs  = "GBS"
                    flowpei  = "PeiCreep"
                    flownlp  = "NLP"
                    flowcomp = "CompoRheo"
                    flowpara = "Parallel"
                    a_j = a[j]
                    # Checks what type the CreepLaw or CompositeRheology field has 
                    if typeof(a_j) <: DislocationCreep
                        Disl, flowlaw, flowlawcount, refs = checkRheologyType(a_j, Disl, flowlaw, flowlawcount, refs, i)

                    elseif typeof(a_j) <: DiffusionCreep
                        Diff, flowlaw, flowlawcount, refs = checkRheologyType(a_j, Diff, flowlaw, flowlawcount, refs, i)

                    elseif typeof(a_j) <: GrainBoundarySliding
                        GBS, flowlaw, flowlawcount, refs = checkRheologyType(a_j, GBS, flowlaw, flowlawcount, refs, i)

                    elseif typeof(a_j) <: PeierlsCreep
                        Pei, flowlaw, flowlawcount, refs = checkRheologyType(a_j, Pei, flowlaw, flowlawcount, refs, i)

                    elseif typeof(a_j) <: NonLinearPeierlsCreep
                        NLP, flowlaw, flowlawcount, refs = checkRheologyType(a_j, NLP, flowlaw, flowlawcount, refs, i)

                    elseif typeof(a_j) <: LinearViscous
                        Lin += 1
                        flowlawcount += 1
                        flowadd = flowlin
                        flowlaw *= flowadd

                    elseif typeof(a_j) <: CompositeRheology
                        compos = getproperty(a_j, :elements)
                        num_rheologies = length(compos)
                        fdsname = "Comp "
                        flowadd = flowcomp
                        comporheo = flowadd * "("
                        for u in 1:num_rheologies
                            a_ju = a_j[u]
                            if typeof(a_ju) <: Parallel
                                num_parallel = length(a_ju.elements)
                                flowadd = flowpara
                                parallelrheo = flowadd * "("
                                namepre = fdsname * "Para "
                                for v in 1:num_parallel
                                    a_juv = a_ju[v]
                                    if typeof(a_juv) <: DislocationCreep
                                        Disl, parallelrheo, flowlawcount, refs = checkRheologyType(a_juv, Disl, parallelrheo, flowlawcount, refs, i; pre=namepre)
                                        
                                    elseif typeof(a_juv) <: DiffusionCreep
                                        Diff, parallelrheo, flowlawcount, refs = checkRheologyType(a_juv, Diff, parallelrheo, flowlawcount, refs, i; pre=namepre)

                                    elseif typeof(a_juv) <: GrainBoundarySliding
                                        GBS, parallelrheo, flowlawcount, refs = checkRheologyType(a_juv, GBS, parallelrheo, flowlawcount, refs, i; pre=namepre)
                
                                    elseif typeof(a_juv) <: PeierlsCreep
                                        Pei, parallelrheo, flowlawcount, refs = checkRheologyType(a_juv, Pei, parallelrheo, flowlawcount, refs, i; pre=namepre)
                
                                    elseif typeof(a_juv) <: NonLinearPeierlsCreep
                                        NLP, parallelrheo, flowlawcount, refs = checkRheologyType(a_juv, NLP, parallelrheo, flowlawcount, refs, i; pre=namepre)
                                        
                                    elseif typeof(a_juv) <: LinearViscous
                                        Lin += 1
                                        flowlawcount += 1
                                        flowadd = flowlin
                                        parallelrheo *= flowadd * ","
                                        
                                    end

                                    varnames = propertynames(a_juv)
                                    # Goes through all variables in a field component
                                    for var in varnames
                                        b = getproperty(a_juv, var)
                                        # Checks if they are GeoUnit (Value and Unit exist then)
                                        if isa(b, GeoUnit)
                                            value = string(getproperty(b, :val))
                                            mdvar = "$var"
                                            # Put value, LaTex variable name, creep law pattern, phase number, current disl, diff or linvisc count and corresponding string in a Dict
                                            if typeof(a_juv) <: DislocationCreep
                                                fds["$var $label $flowadd $i.$u.$v"] = (value, mdvar, "$flowlaw$comporheo$parallelrheo)", "$i", "$Disl", "$flowdisl")
                                            end
                                            if typeof(a_juv) <: DiffusionCreep
                                                fds["$var $label $flowadd $i.$u.$v"] = (value, mdvar, "$flowlaw$comporheo$parallelrheo)", "$i", "$Diff", "$flowdiff")
                                            end
                                            if typeof(a_juv) <: LinearViscous
                                                fds["$var $label $flowadd $i.$u.$v"] = (value, mdvar, "$flowlaw$comporheo$parallelrheo)", "$i", "$Lin", "$flowlin")
                                            end
                                        end
                                        k += 1
                                    end
                                end
                                comporheo *= parallelrheo * "),"

                            elseif typeof(a_ju) <: DislocationCreep
                                Disl, comporheo, flowlawcount, refs = checkRheologyType(a_ju, Disl, comporheo, flowlawcount, refs, i; pre=fdsname)

                            elseif typeof(a_ju) <: DiffusionCreep
                                Diff, comporheo, flowlawcount, refs = checkRheologyType(a_ju, Diff, comporheo, flowlawcount, refs, i; pre=fdsname)

                            elseif typeof(a_ju) <: GrainBoundarySliding
                                GBS, comporheo, flowlawcount, refs = checkRheologyType(a_ju, GBS, comporheo, flowlawcount, refs, i; pre=fdsname)

                            elseif typeof(a_ju) <: PeierlsCreep
                                Pei, comporheo, flowlawcount, refs = checkRheologyType(a_ju, Pei, comporheo, flowlawcount, refs, i; pre=fdsname)

                            elseif typeof(a_ju) <: NonLinearPeierlsCreep
                                NLP, comporheo, flowlawcount, refs = checkRheologyType(a_ju, NLP, comporheo, flowlawcount, refs, i; pre=fdsname)

                            elseif typeof(a_ju) <: LinearViscous
                                Lin += 1
                                flowlawcount += 1
                                flowadd = flowlin
                                comporheo *= flowadd * ","
                            end

                            varnames = propertynames(a_ju)
                            # Goes through all variables in a field component
                            c = 0
                            for var in varnames
                                c += 1
                                b = getproperty(a_ju, var)
                                # Checks if they are GeoUnit (Value and Unit exist then)
                                if isa(b, GeoUnit)
                                    value = string(getproperty(b, :val))
                                    mdvar = "$var"
                                    # Put value, LaTex variable name, creep law pattern, phase number, current disl, diff or linvisc count and corresponding string in a Dict
                                    if typeof(a_ju) <: DislocationCreep
                                        fds["$var $label $flowadd $i.$u"] = (value, mdvar, "$flowlaw$comporheo)", "$i", "$Disl", "$flowdisl")
                                    end
                                    if typeof(a_ju) <: DiffusionCreep
                                        fds["$var $label $flowadd $i.$u"] = (value, mdvar, "$flowlaw$comporheo)", "$i", "$Diff", "$flowdiff")
                                    end
                                    if typeof(a_ju) <: LinearViscous
                                        fds["$var $label $flowadd $i.$u"] = (value, mdvar, "$flowlaw$comporheo)", "$i", "$Lin", "$flowlin")
                                    end
                                end
                                k += 1
                            end
                        end
                        flowlaw *= comporheo * "),"

                    elseif typeof(a_j) <: Parallel
                        num_parallel = length(a_j.elements)
                        fdsname = "Para "
                        flowadd = flowpara
                        parallelrheo = flowadd * "("
                        for q in 1:num_parallel
                            a_jq = a_j[q]
                            if typeof(a_jq) <: DislocationCreep
                                Disl, parallelrheo, flowlawcount, refs = checkRheologyType(a_jq, Disl, parallelrheo, flowlawcount, refs, i; pre=fdsname)

                            elseif typeof(a_jq) <: DiffusionCreep
                                Diff, parallelrheo, flowlawcount, refs = checkRheologyType(a_jq, Diff, parallelrheo, flowlawcount, refs, i; pre=fdsname)

                            elseif typeof(a_jq) <: GrainBoundarySliding
                                GBS, parallelrheo, flowlawcount, refs = checkRheologyType(a_jq, GBS, parallelrheo, flowlawcount, refs, i; pre=fdsname)

                            elseif typeof(a_jq) <: PeierlsCreep
                                Pei, parallelrheo, flowlawcount, refs = checkRheologyType(a_jq, Pei, parallelrheo, flowlawcount, refs, i; pre=fdsname)

                            elseif typeof(a_jq) <: NonLinearPeierlsCreep
                                NLP, parallelrheo, flowlawcount, refs = checkRheologyType(a_jq, NLP, parallelrheo, flowlawcount, refs, i; pre=fdsname)

                            elseif typeof(a_jq) <: LinearViscous
                                Lin += 1
                                flowlawcount += 1
                                flowadd = flowlin
                                parallelrheo *= flowadd * ","

                            elseif typeof(a_jq) <: CompositeRheology
                                compos = getproperty(a_j, :elements)
                                num_rheologies = length(compos)
                                namepre = fdsname * "Comp "
                                flowadd = flowcomp
                                comporheo = flowadd * "("
                                for u in 1:num_rheologies
                                    a_jqu = a_jq[u]
                                    if typeof(a_jqu) <: DislocationCreep
                                        Disl, comporheo, flowlawcount, refs = checkRheologyType(a_jqu, Disl, comporheo, flowlawcount, refs, i; pre=namepre)

                                    elseif typeof(a_jqu) <: DiffusionCreep
                                        Diff, comporheo, flowlawcount, refs = checkRheologyType(a_jqu, Diff, comporheo, flowlawcount, refs, i; pre=namepre)

                                    elseif typeof(a_jqu) <: GrainBoundarySliding
                                        GBS, comporheo, flowlawcount, refs = checkRheologyType(a_jqu, GBS, comporheo, flowlawcount, refs, i; pre=namepre)

                                    elseif typeof(a_jqu) <: PeierlsCreep
                                        Pei, comporheo, flowlawcount, refs = checkRheologyType(a_jqu, Pei, comporheo, flowlawcount, refs, i; pre=namepre)

                                    elseif typeof(a_jqu) <: NonLinearPeierlsCreep
                                        NLP, flowlaw, flowlawcount, refs = checkRheologyType(a_jqu, NLP, flowlaw, flowlawcount, refs, i; pre=namepre)

                                    elseif typeof(a_jqu) <: LinearViscous
                                        Lin += 1
                                        flowlawcount += 1
                                        flowadd = flowlin
                                        comporheo *= flowadd * ","
                                    end

                                    varnames = propertynames(a_jqu)
                                    # Goes through all variables in a field component
                                    for var in varnames
                                        b = getproperty(a_jqu, var)
                                        # Checks if they are GeoUnit (Value and Unit exist then)
                                        if isa(b, GeoUnit)
                                            value = string(getproperty(b, :val))
                                            mdvar = "$var"
                                            # Put value, LaTex variable name, creep law pattern, phase number, current disl, diff or linvisc count and corresponding string in a Dict
                                            if typeof(a_jqu) <: DislocationCreep
                                                fds["$var $label $flowadd $i.$q.$u"] = (value, mdvar, "$flowlaw$parallelrheo$comporheo)", "$i", "$Disl", "$flowdisl")
                                            end
                                            if typeof(a_jqu) <: DiffusionCreep
                                                fds["$var $label $flowadd $i.$q.$u"] = (value, mdvar, "$flowlaw$parallelrheo$comporheo)", "$i", "$Diff", "$flowdiff")
                                            end
                                            if typeof(a_jqu) <: GrainBoundarySliding
                                                fds["$var $label $flowadd $i.$q.$u"] = (value, mdvar, "$flowlaw$parallelrheo$comporheo)", "$i", "$GBS", "$flowgbs")
                                            end
                                            if typeof(a_jqu) <: PeierlsCreep
                                                fds["$var $label $flowadd $i.$q.$u"] = (value, mdvar, "$flowlaw$parallelrheo$comporheo)", "$i", "$Pei", "$flowpei")
                                            end
                                            if typeof(a_jqu) <: NonLinearPeierlsCreep
                                                fds["$var $label $flowadd $i.$q.$u"] = (value, mdvar, "$flowlaw$parallelrheo$comporheo)", "$i", "$NLP", "$flownlp")
                                            end
                                            if typeof(a_jqu) <: LinearViscous
                                                fds["$var $label $flowadd $i.$q.$u"] = (value, mdvar, "$flowlaw$parallelrheo$comporheo)", "$i", "$Lin", "$flowlin")
                                            end
                                        end
                                        k += 1
                                    end

                                end
                                parallelrheo *= comporheo * "),"
                            end

                            varnames = propertynames(a_jq)
                            # Goes through all variables in a field component
                            for var in varnames
                                b = getproperty(a_jq, var)
                                # Checks if they are GeoUnit (Value and Unit exist then)
                                if isa(b, GeoUnit)
                                    value = string(getproperty(b, :val))
                                    mdvar = "$var"
                                    # Put value, LaTex variable name, creep law pattern, phase number, current disl, diff or linvisc count and corresponding string in a Dict
                                    if typeof(a_jq) <: DislocationCreep
                                        fds["$var $label $flowadd $i.$q"] = (value, mdvar, "$flowlaw$parallelrheo)", "$i", "$Disl", "$flowdisl")
                                    end
                                    if typeof(a_jq) <: DiffusionCreep
                                        fds["$var $label $flowadd $i.$q"] = (value, mdvar, "$flowlaw$parallelrheo)", "$i", "$Diff", "$flowdiff")
                                    end
                                    if typeof(a_jq) <: GrainBoundarySliding
                                        fds["$var $label $flowadd $i.$q"] = (value, mdvar, "$flowlaw$parallelrheo)", "$i", "$GBS", "$flowgbs")
                                    end
                                    if typeof(a_jq) <: PeierlsCreep
                                        fds["$var $label $flowadd $i.$q"] = (value, mdvar, "$flowlaw$parallelrheo)", "$i", "$Pei", "$flowpei")
                                    end
                                    if typeof(a_jq) <: NonLinearPeierlsCreep
                                        fds["$var $label $flowadd $i.$q"] = (value, mdvar, "$flowlaw$parallelrheo)", "$i", "$NLP", "$flownlp")
                                    end
                                    if typeof(a_jq) <: LinearViscous
                                        fds["$var $label $flowadd $i.$q"] = (value, mdvar, "$flowlaw$parallelrheo)", "$i", "$Lin", "$flowlin")
                                    end
                                end
                                k += 1
                            end

                        end
                        flowlaw *= parallelrheo * "),"
                    end
                    varnames = propertynames(a_j)
                    # Goes through all variables in a field component
                    if !(typeof(a_j) <: CompositeRheology || typeof(a_j) <: Parallel)
                        for var in varnames
                            b = getproperty(a_j, var)
                            # Checks if they are GeoUnit (Value and Unit exist then)
                            if isa(b, GeoUnit)
                                value = string(getproperty(b, :val))
                                mdvar = "$var"
                                # Put value, LaTex variable name, creep law pattern, phase number, current disl, diff or linvisc count and corresponding string in a Dict
                                if typeof(a_j) <: DislocationCreep
                                    fds["$var $label $i"] = (value, mdvar, "$flowlaw", "$i", "$Disl", "$flowdisl")
                                elseif typeof(a_j) <: DiffusionCreep
                                    fds["$var $label $i"] = (value, mdvar, "$flowlaw", "$i", "$Diff", "$flowdiff")
                                elseif typeof(a_j) <: GrainBoundarySliding
                                    fds["$var $label $i"] = (value, mdvar, "$flowlaw", "$i", "$GBS", "$flowgbs")
                                elseif typeof(a_j) <: PeierlsCreep
                                    fds["$var $label $i"] = (value, mdvar, "$flowlaw", "$i", "$Pei", "$flowpei")
                                elseif typeof(a_j) <: NonLinearPeierlsCreep
                                    fds["$var $label $i"] = (value, mdvar, "$flowlaw", "$i", "$NLP", "$flownlp")
                                elseif typeof(a_j) <: LinearViscous
                                    fds["$var $label $i"] = (value, mdvar, "$flowlaw", "$i", "$Lin", "$flowlin")
                                else
                                    fds["$var $label $i"] = (value, mdvar, "", "$i", "", "")
                                end
                                
                            end
                            k += 1
                        end
                    end
                end
                # Takes field "Name" and puts phase name, maximum phase  count and current phase number in Dict 
            elseif !isempty(getproperty(s[i], label)) && label == :Name
                phasename = join(getproperty(s[i], :Name))
                fds["$label $i"] = (phasename, "$phasecount", "$i", "", "", "")
            end
        end
    end
    return fds
end


"""
Dict2MarkdownTable() writes a .md file with all parameters from the Phase2DictMd() output in a Markdown table. rdigits will round numbers with more decimals than rdigits
including numbers of 10 to power of n, n being an Integer, for representation purposes. For the exact numbers use the original implemented numbers from the creeplaws of the
dict in src/CreepLaw/Data/DiffusionCreep.jl or src/CreepLaw/Data/DislocationCreep.jl.

"""
function Dict2MarkdownTable(d::Dict; filename="ParameterTable", rdigits=4)
    symbs = []

    # Creates vectors of type Pairs (can be iterated over, sorted, etc.)
    dictpairs = sort(collect(pairs(d)))

    # Descriptions for every parameter that could occur in the table and their corresponding variable name(s) that is used in GeoParams
    desc = Dict(
        "ρ"=>"Density *(kg/m^3^)*",
        "ρ0"=>"Reference density *(kg/m^3^)*",
        "g"=>"Gravity *(m/s^2^)*",
        "η"=>"Viscosity *(Pa s)*",
        "P"=>"Pressure *(MPa)*",
        "T"=>"Temperature *(°C)*",
        "V"=>"Volume *(m^3^)*" ,
        "d"=>"Grain size *(cm)*",
        "f"=>"Water fugacity *(MPa)*",
        "n"=>"Power-law exponent *(-)*",
        "r"=>"Water fugacity exponent *(-)*",
        "p"=>"Grain size exponent *(-)*",
        "o"=>"Stres parametrization exponent *(-)*",    # in PeierlsCreep
        "q"=>"Stress relation exponent *(-)*",          # in PeierlsCreep
        "A_diff" => "Coefficient *(Pa^-n-r^ m^p^/s)*",  # DiffusionCreep
        "A_disl" => "Coefficient *(Pa^-n^/s)*",         # DislocationCreep
        "A_gbs" => "Coefficient *(Pa^-n^/s)*",          # GrainBoundarySliding
        "A_pei" => "Coefficient *(Pa^-n^/s)*",          # PeierlsCreep
        "A_nlp" => "Coefficient *(Pa^-n^/s)*",          # NonLinearPeierlsCreep
        "Tau_p" => "Peierls stress *(Pa)*",
        "E"=>"Activation energy *(J/mol)*",
        "R"=>"Gas constant *(J/mol/K)*",
        "G"=>"Shear modulus *(Pa)*",
        "ν"=>"Poisson ratio *(-)*",
        "K"=>"Bulk modulus *(Pa)*",
        "Kb"=>"Elastic bulk modulus *(Pa)*", 
        "Y"=>"Young's modulus *(Pa)*",
        "k"=>"Thermal conductivity *(W/m/K)*",
        "cp"=>"Heat capacity *(J/kg/K)*",
        "Q_L"=>"Latent heat *(kJ/kg)*",
        "H_r"=>"Radioactive heat *(W/m^3^)*",
        "H_s"=>"Shear heating *(-)*",
        "ϕ"=>"Friction angle *(°)*",
        "ψ"=>"Dilation angle *(°)*",
        "C"=>"Cohesion *(Pa)*",
        "Vp"=>"P-wave velocity *(km/s)*",
        "Vs"=>"S-wave velocity *(km/s)*",
        "T0"=>"Reference temperature *(°C)*",
        "P0"=>"Reference pressure *(Pa)*",
        "β"=>"Compressibility *(1/Pa)*",
        "α"=>"Thermal expansion coeff. *(1/K)*"
        )

    # Generates latex preamble
    Table = " | Denotation | Variablename | "

    # Creates table headers
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
    for key in 1:length(dictpairs)
        dictpairs_key = dictpairs[key].first
        if occursin("Name", dictpairs_key)
            continue
        # Checks if symbol is from field CompositeRheology or CreepLaws, if symbol is R and if the symbol occurs more often than once
        elseif (occursin("R", dictpairs_key[1:3]))
            push!(symbs, dictpairs[key].second[2])
        elseif maximum(occursin.(["CompositeRheology", "CreepLaws"], dictpairs_key)) && !(occursin("R", dictpairs_key[1:3]))
            if occursin("Disl", dictpairs[key].second[6])
                number = parse(Int64, dictpairs[key].second[5])
                push!(symbs, "A_disl" * "^$number")
            elseif occursin("Diff", dictpairs[key].second[6])
                number = parse(Int64, dictpairs[key].second[5])
                push!(symbs, "A_diff" * "^$number")
            else
                number = parse(Int64, dictpairs[key].second[5])
                push!(symbs, dictpairs[key].second[2] * "$number")
            end
        else
            push!(symbs, dictpairs[key].second[2])
        end
    end
    symbs = unique(sort(symbs))

    # Creates columnwise output for all parameters of the input phase
    for symbol in symbs
        # Sets parametername and variable
        if length(symbol) > 1
            # Checks if Unicode chars in combination with a number etc. (e.g "α0") are used which apparently do not have a second index but only first and third
            # and checks if symbol has a number in it 
            if maximum(endswith(symbol, string(i),) for i in 0:9) && !occursin("_", symbol)
                if length(unidecode(symbol)) > 2 && occursin("0", symbol)
                    Table *= " " * string(desc[symbol]) * " | " * symbol[1] * "~" * symbol[3] * "~"
                elseif length(unidecode(symbol)) <= 2 && occursin("0", symbol)
                    Table *= " " * string(desc[symbol]) * " | " * symbol[1] * "~" * symbol[2] * "~"
                elseif length(unidecode(symbol)) > 2 && maximum(occursin(string(i), symbol) for i in 1:9)
                    symbol_1 = symbol[1]
                    Table *= " " * string(desc["$symbol_1"]) * " | " * symbol[1] * "~" * symbol[3] * "~"
                elseif length(unidecode(symbol)) <= 2 && maximum(occursin(string(i), symbol) for i in 1:9)
                    symbol_1 = symbol[1:end-1]
                    Table *= " " * string(desc["$symbol_1"]) * " | " * symbol[1] * "~" * symbol[2] * "~"
                end
            # If "_" AND "^" in symbol 
            elseif occursin("_", symbol) && occursin("^", symbol)
                occ_under = findfirst('_', symbol)
                occ_roof = findfirst('^', symbol)
                Table *= " " * string(desc[symbol[1:occ_roof-1]]) *  " | " * symbol[1:occ_under-1] * symbol[occ_roof] * symbol[occ_roof+1:end] * "^" * "~" * symbol[occ_under+1:occ_roof-1] * "~"
            # If only "_" in symbol
            elseif occursin("_", symbol)
                occ_under = findfirst('_', symbol)
                Table *= " " * string(desc[symbol]) *  " | " * symbol[1:occ_under-1] * "~" * symbol[occ_under+1:end] * "~"
            end
        # If neither a number nor one of the above signs in symbol, so just a normal symbol
        else
            Table *= " " * string(desc[symbol]) * " | " * symbol
        end

        # Iterates over all phases
        for j = 1:parse(Int64, d["Name 1"][2])
            hit = 0
            i_dictpairs = 1
            # Iterates over all pairs
            for element in dictpairs
                if maximum(endswith.(string(symbol),[string(i) for i in 1:9])) && !occursin("A_", string(symbol))
                    symbol = replace(symbol, symbol[end] => "")
                elseif occursin("A_", string(symbol)) && ((occursin("diff", symbol) && (occursin("Diff", dictpairs[i_dictpairs].second[6]))) || (occursin("disl", symbol) && occursin("Disl", dictpairs[i_dictpairs].second[6])))
                    symbol_A = string(symbol[1])
                    if symbol_A == dictpairs[i_dictpairs].second[2] && j == parse(Int64, dictpairs[i_dictpairs].second[4]) && hit == 0
                        # put in the matched parameter value
                        dig, num, expo = detachFloatfromExponent(dictpairs[i_dictpairs].second[1])
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
                        deleteat!(dictpairs, i_dictpairs)
                    end
                end
                # Checks if symbol matches the symbol in the pair and if phase matches phase of the pair
                if symbol == dictpairs[i_dictpairs].second[2] && j == parse(Int64, dictpairs[i_dictpairs].second[4]) && hit == 0 && !(occursin("A_", string(symbol)))
                    # put in the matched parameter value
                    dig, num, expo = detachFloatfromExponent(dictpairs[i_dictpairs].second[1])
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
                    deleteat!(dictpairs, i_dictpairs)
                end
                i_dictpairs += 1
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
Creates a filename.tex or filename.md file as output. If "latex" is chosen a "Reference.bib" file will automatically be produced with all your references. Storage path
by default is the GeoParams package folder. A specific storage path can be given as normal path with "\\filename" in the end as a string. Use double backslash for subfolders.
There is no file extension needed to be given.
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
        print("Created $filename.tex and References.bib files. \n")
    elseif (format == "markdown") || (format == "md")
        d = Phase2DictMd(Phase)
        Dict2MarkdownTable(d, filename=filename, rdigits=rdigits)
        print("Created $filename.md file. \n")
    end
    return nothing
end

end
