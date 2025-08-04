module Tables

using Unidecode
using GeoParams: AbstractMaterialParam, param_info, LinearViscous, PowerlawViscous, DislocationCreep, DiffusionCreep, CompositeRheology, Parallel, make_tuple
using GeoParams: GrainBoundarySliding, PeierlsCreep, NonLinearPeierlsCreep, LinearMeltViscosity, ViscosityPartialMelt_Costa_etal_2009, GiordanoMeltViscosity
using ..Units
using ..MaterialParameters: MaterialParamsInfo
import GeoParams: Dislocation, Diffusion, GBS, Peierls, NonLinearPeierls
import GeoParams.Dislocation: dislocation_database_info
import GeoParams.Diffusion: diffusion_database_info
import GeoParams.GBS: GrainBoundarySliding_database_info
import GeoParams.Peierls: peierls_database_info
import GeoParams.NonLinearPeierls: nonlinear_peierls_database_info

export detachFloatfromExponent, extract_parameters_from_phases, Dict2LatexTable, extract_parameters_from_phases_md, Dict2MarkdownTable, ParameterTable

# Internal types for better organization
struct ParameterEntry
    value::String
    symbol::String
    flowlaw::String
    phase::String
    rheology_type::String
    rheology_subtype::String
end

struct ReferenceEntry
    bibtex::String
    index::String
    phase::String
end

"""
    create_latex_symbol(field_name) -> String

Create a LaTeX symbol that matches the keys in the parameter descriptions dictionary.
"""
function create_latex_symbol(field_name)
    field_str = string(field_name)

    # Define symbol mappings that match get_parameter_descriptions() keys
    symbol_map = Dict(
        "ρ" => "\\rho",
        "rho" => "\\rho",
        "ρ0" => "\\rho_0",
        "rho0" => "\\rho_0",
        "η" => "\\eta",
        "eta" => "\\eta",
        "η0" => "\\eta_0",
        "eta0" => "\\eta_0",
        "ν" => "\\nu",
        "nu" => "\\nu",
        "ϕ" => "\\phi",
        "phi" => "\\phi",
        "ψ" => "\\psi",
        "psi" => "\\psi",
        "β" => "\\beta",
        "beta" => "\\beta",
        "α" => "\\alpha",
        "alpha" => "\\alpha",
        "T0" => "T_0",
        "P0" => "P_0",
        "A_diff" => "A_{\\text{diff}}",
        "A_disl" => "A_{\\text{disl}}",
        "A" => "A",
        "B" => "B"
    )

    # Check if we have a direct mapping
    if haskey(symbol_map, field_str)
        return symbol_map[field_str]
    end

    # For common single letters that should NOT have backslashes (these are parameter names, not LaTeX symbols)
    if field_str in ["g", "P", "T", "V", "d", "f", "n", "r", "p", "E", "R", "G", "K", "Y", "k", "C"]
        return field_str
    end

    # For multi-character parameters like A_diff, A_disl, etc.
    if length(field_str) > 2 || contains(field_str, "_")
        return field_str  # No backslash for compound names
    end

    # For single letters that might be Greek symbols, add backslash
    if field_str in ["lambda", "mu", "sigma", "tau", "gamma", "delta", "epsilon"]
        return "\\" * field_str
    end

    # Default: return as-is for other parameters
    return field_str
end


# Helper function to determine rheology type
function get_rheology_info(param::AbstractMaterialParam)
    if param isa DislocationCreep
        return "DislCreep", "DislCreep"
    elseif param isa DiffusionCreep
        return "DiffCreep", "DiffCreep"
    elseif param isa LinearViscous
        return "LinVisc", "LinVisc"
    elseif param isa PowerlawViscous
        return "PowerVisc", "PowerVisc"
    elseif param isa GrainBoundarySliding
        return "GBS", "GBS"
    elseif param isa PeierlsCreep
        return "PeierlsCreep", "PeierlsCreep"
    elseif param isa NonLinearPeierlsCreep
        return "NonLinPeierls", "NonLinPeierls"
    elseif param isa LinearMeltViscosity
        return "MeltVisc", "MeltVisc"
    elseif param isa ViscosityPartialMelt_Costa_etal_2009
        return "PartialMelt", "PartialMelt"
    elseif param isa GiordanoMeltViscosity
        return "GiordanoMelt", "GiordanoMelt"
    else
        # For unknown types, try to extract a meaningful name from the type
        type_name = string(typeof(param).name.name)
        return type_name, type_name
    end
end


"""
    detachFloatfromExponent(str::String) -> (Int, String, String)

Returns the number of decimal places, the number without exponent, and the exponent.
Returns "1" for the exponent if the input number has no exponent.
"""
function detachFloatfromExponent(str::String)
    s = lowercase(str)
    if 'e' in s
        parts = split(s, "e")
        s = string(parts[1])
        ex = string(parts[2])
    else
        ex = "1"
    end

    # Count decimal places
    decimal_idx = findfirst('.', s)
    if decimal_idx !== nothing
        dig = length(s) - decimal_idx
    else
        dig = 0
    end

    return dig, s, ex
end


"""
    extract_parameters_from_material(material, phase_idx::Int, counters::Dict)

Extract parameters from a single material parameter object recursively.
"""
function extract_parameters_from_material(material::AbstractMaterialParam, phase_idx::Int, counters::Dict, params::Dict, refs::Dict, prefix::String = "")
    rheo_type, rheo_subtype = get_rheology_info(material)

    # Only treat actual rheology laws as rheology types, not basic material properties
    is_rheology = material isa Union{
        DislocationCreep, DiffusionCreep, LinearViscous, PowerlawViscous,
        GrainBoundarySliding, PeierlsCreep, NonLinearPeierlsCreep,
        LinearMeltViscosity, ViscosityPartialMelt_Costa_etal_2009, GiordanoMeltViscosity,
    }

    # Update counters if this is a rheology law
    if is_rheology
        counters[rheo_type] = get(counters, rheo_type, 0) + 1
        counters["flowlaw"] = get(counters, "flowlaw", 0) + 1

        # Add reference if available
        try
            info = get_material_reference_info(material)

            if info !== nothing && hasfield(typeof(info), :BibTex_Reference) && !isempty(info.BibTex_Reference)
                name = unsafe_string(material.Name)
                key = prefix * name
                refs[key] = (info.BibTex_Reference, string(counters["flowlaw"]), string(phase_idx))
            end
        catch
            # Skip if reference extraction fails
        end
    else
        # For non-rheology parameters, use empty strings
        rheo_type = ""
        rheo_subtype = ""
    end

    # Extract parameter values
    for field_name in propertynames(material)
        field_val = getproperty(material, field_name)

        if field_val isa GeoUnit
            value_str = string(field_val.val)
            symbol_str = string(field_name)

            # Create LaTeX symbol for display - match the keys in get_parameter_descriptions()
            latex_symbol = create_latex_symbol(field_name)

            # Generate key in expected format
            # For composite rheologies, use special formatting: "fieldname FieldType RheoType phase.element"
            if !isempty(prefix) && occursin("Comp_", prefix)
                # Extract field name and element index from prefix
                # e.g., "CompositeRheology_Comp_3_" -> field="CompositeRheology", element=3
                prefix_parts = split(strip(prefix, '_'), "_")
                field_type = prefix_parts[1]  # "CompositeRheology"
                element_idx = prefix_parts[end]  # "3"

                # Use rheology type and format as phase.element
                key = "$(field_name) $(field_type) $(rheo_subtype) $(phase_idx).$(element_idx)"
            else
                # Standard format for non-composite rheologies
                field_type = if !isempty(prefix) && endswith(prefix, "_")
                    strip(prefix, '_')
                else
                    "param"
                end
                key = "$(field_name) $(field_type) $(phase_idx)"
            end

            params[key] = (
                value_str,
                latex_symbol,
                rheo_type,
                string(phase_idx),
                isempty(rheo_type) ? "" : string(get(counters, rheo_type, 0)),
                rheo_subtype,
            )
        end
    end
    return
end

"""
    extract_parameters_recursive(obj, phase_idx::Int, counters::Dict, params::Dict, refs::Dict, prefix::String)

Recursively extract parameters from nested composite rheologies and parallel combinations.
"""
function extract_parameters_recursive(obj, phase_idx::Int, counters::Dict, params::Dict, refs::Dict, prefix::String = "")
    return if obj isa CompositeRheology
        # Build composite rheology summary
        composite_summary = build_composite_summary(obj)
        new_prefix = prefix * "Comp_"
        for (i, element) in enumerate(obj.elements)
            extract_parameters_recursive_with_composite(element, phase_idx, counters, params, refs, new_prefix * "$(i)_", composite_summary)
        end
    elseif obj isa Parallel
        new_prefix = prefix * "Para_"
        for (i, element) in enumerate(obj.elements)
            extract_parameters_recursive(element, phase_idx, counters, params, refs, new_prefix * "$(i)_")
        end
    elseif obj isa AbstractMaterialParam
        # Base case - single material parameter
        extract_parameters_from_material(obj, phase_idx, counters, params, refs, prefix)
    end
    # If obj is not a material parameter type (e.g., Int, String), skip it
end

"""
    build_composite_summary(composite_rheology) -> String

Build a summary string describing the structure of a composite rheology using generic type detection.
"""
function build_composite_summary(obj::CompositeRheology)
    element_types = String[]
    seen_types = Set{String}()

    for element in obj.elements
        # Use the existing generic get_rheology_info function
        _, rheo_subtype = get_rheology_info(element)

        # If get_rheology_info doesn't recognize it, fall back to type name
        element_type = if !isempty(rheo_subtype)
            rheo_subtype
        else
            string(typeof(element).name.name)
        end

        # Only add if we haven't seen this type before
        if !(element_type in seen_types)
            push!(element_types, element_type)
            push!(seen_types, element_type)
        end
    end

    return "CompoRheo(" * join(element_types, ",") * ",)"
end

"""
    extract_parameters_recursive_with_composite(obj, phase_idx, counters, params, refs, prefix, composite_summary)

Extract parameters with composite rheology summary information.
"""
function extract_parameters_recursive_with_composite(obj, phase_idx::Int, counters::Dict, params::Dict, refs::Dict, prefix::String, composite_summary::String)
    return if obj isa AbstractMaterialParam
        # Extract parameters with composite summary
        extract_parameters_from_material_with_composite(obj, phase_idx, counters, params, refs, prefix, composite_summary)
    end
end

"""
    extract_parameters_from_material_with_composite(material, phase_idx, counters, params, refs, prefix, composite_summary)

Extract parameters from a single material parameter object with composite rheology information.
"""
function extract_parameters_from_material_with_composite(material::AbstractMaterialParam, phase_idx::Int, counters::Dict, params::Dict, refs::Dict, prefix::String, composite_summary::String)
    rheo_type, rheo_subtype = get_rheology_info(material)

    # Only treat actual rheology laws as rheology types, not basic material properties
    is_rheology = material isa Union{
        DislocationCreep, DiffusionCreep, LinearViscous, PowerlawViscous,
        GrainBoundarySliding, PeierlsCreep, NonLinearPeierlsCreep,
        LinearMeltViscosity, ViscosityPartialMelt_Costa_etal_2009, GiordanoMeltViscosity,
    }

    # Update counters if this is a rheology law
    if is_rheology
        counters[rheo_type] = get(counters, rheo_type, 0) + 1
        counters["flowlaw"] = get(counters, "flowlaw", 0) + 1

        # Add reference if available
        try
            info = get_material_reference_info(material)

            if info !== nothing && hasfield(typeof(info), :BibTex_Reference) && !isempty(info.BibTex_Reference)
                name = unsafe_string(material.Name)
                key = prefix * name
                refs[key] = (info.BibTex_Reference, string(counters["flowlaw"]), string(phase_idx))
            end
        catch
            # Skip if reference extraction fails
        end
    else
        # For non-rheology parameters, use empty strings
        rheo_type = ""
        rheo_subtype = ""
    end

    # Extract parameter values
    for field_name in propertynames(material)
        field_val = getproperty(material, field_name)

        if field_val isa GeoUnit
            value_str = string(field_val.val)
            symbol_str = string(field_name)

            # Create LaTeX symbol for display - match the keys in get_parameter_descriptions()
            latex_symbol = create_latex_symbol(field_name)

            # Generate key in expected format
            # For composite rheologies, use special formatting: "fieldname FieldType RheoType phase.element"
            if !isempty(prefix) && occursin("Comp_", prefix)
                # Extract field name and element index from prefix
                # e.g., "CompositeRheology_Comp_3_" -> field="CompositeRheology", element=3
                prefix_parts = split(strip(prefix, '_'), "_")
                field_type = prefix_parts[1]  # "CompositeRheology"
                element_idx = prefix_parts[end]  # "3"

                # Use rheology type and format as phase.element
                key = "$(field_name) $(field_type) $(rheo_subtype) $(phase_idx).$(element_idx)"
            else
                # Standard format for non-composite rheologies
                field_type = if !isempty(prefix) && endswith(prefix, "_")
                    strip(prefix, '_')
                else
                    "param"
                end
                key = "$(field_name) $(field_type) $(phase_idx)"
            end

            params[key] = (
                value_str,
                latex_symbol,
                composite_summary,  # Use composite summary instead of individual rheo_type
                string(phase_idx),
                isempty(rheo_type) ? "" : string(get(counters, rheo_type, 0)),
                rheo_subtype,
            )
        end
    end
    return
end

"""
    get_rheology_full_name(rheo_type) -> String

Convert abbreviated rheology type to full display name.
"""
function get_rheology_full_name(rheo_type::String)
    type_names = Dict(
        "DislCreep" => "Dislocation Creep",
        "DiffCreep" => "Diffusion Creep",
        "LinVisc" => "Linear Viscous",
        "PowerVisc" => "Power-law Viscous",
        "GBS" => "Grain Boundary Sliding",
        "PeierlsCreep" => "Peierls Creep",
        "NonLinPeierls" => "Nonlinear Peierls Creep",
        "MeltVisc" => "Melt Viscosity",
        "PartialMelt" => "Partial Melt Viscosity",
        "GiordanoMelt" => "Giordano Melt Viscosity"
    )
    return get(type_names, rheo_type, rheo_type)
end

"""
    get_category_display_name(field_type) -> String

Convert field type to display name for basic material properties.
"""
function get_category_display_name(field_type::String)
    display_names = Dict(
        "Density" => "Density",
        "Elasticity" => "Elastic Properties",
        "Conductivity" => "Thermal Conductivity",
        "HeatCapacity" => "Heat Capacity",
        "Plasticity" => "Plasticity",
        "Permeability" => "Permeability"
    )
    return get(display_names, field_type, field_type)
end

"""
    get_category_sort_order(field_type) -> Int

Get sort order for basic material property categories.
"""
function get_category_sort_order(field_type::String)
    sort_orders = Dict(
        "Density" => 10,
        "Elasticity" => 20,
        "Conductivity" => 30,
        "HeatCapacity" => 40,
        "Plasticity" => 50,
        "Permeability" => 60
    )
    return get(sort_orders, field_type, 99)
end

"""
    get_material_reference_info(material::AbstractMaterialParam) -> Union{MaterialParamsInfo, Nothing}

Get MaterialParamsInfo for any type of material parameter using the appropriate database_info function.
"""
function get_material_reference_info(material::AbstractMaterialParam)
    try
        # For melt viscosity types, use param_info directly
        if material isa Union{LinearMeltViscosity, ViscosityPartialMelt_Costa_etal_2009, GiordanoMeltViscosity}
            return param_info(material)
        end

        # For creep laws with database pattern, use the database_info approach
        # Get the original function name from the material's Name field
        name_str = if isa(material.Name, Ptr)
            unsafe_string(material.Name)
        else
            string(material.Name)
        end

        # Try to find the function that created this material
        func = find_creep_law_function(material, name_str)
        if func !== nothing
            return get_database_info(material, func)
        end
    catch
        # If any step fails, return nothing
    end
    return nothing
end

"""
    get_database_info(material::AbstractMaterialParam, func::Function) -> Union{MaterialParamsInfo, Nothing}

Dispatch to the appropriate database_info function based on material type.
"""
get_database_info(material::DislocationCreep, func::Function) = dislocation_database_info(func)
get_database_info(material::DiffusionCreep, func::Function) = diffusion_database_info(func)
get_database_info(material::GrainBoundarySliding, func::Function) = GrainBoundarySliding_database_info(func)
get_database_info(material::PeierlsCreep, func::Function) = peierls_database_info(func)
get_database_info(material::NonLinearPeierlsCreep, func::Function) = nonlinear_peierls_database_info(func)
# Fallback for any other types
get_database_info(material::AbstractMaterialParam, func::Function) = nothing

"""
    find_creep_law_function(material::AbstractMaterialParam, name_str::String) -> Union{Function, Nothing}

Find the original function that created this material parameter based on its type and name.
"""
function find_creep_law_function(material::AbstractMaterialParam, name_str::String)
    # Get the appropriate module based on material type
    if material isa DislocationCreep
        return find_function_in_module(Dislocation, name_str)
    elseif material isa DiffusionCreep
        return find_function_in_module(Diffusion, name_str)
    elseif material isa GrainBoundarySliding
        return find_function_in_module(GBS, name_str)
    elseif material isa PeierlsCreep
        return find_function_in_module(Peierls, name_str)
    elseif material isa NonLinearPeierlsCreep
        return find_function_in_module(NonLinearPeierls, name_str)
    end
    return nothing
end

"""
    find_function_in_module(mod::Module, name_str::String) -> Union{Function, Nothing}

Find a function in the given module that produces a material with the given name.
"""
function find_function_in_module(mod::Module, name_str::String)
    try
        # Get all exported functions from the module
        for name in names(mod; all = true, imported = false)
            if !startswith(string(name), "_") && name != :eval && name != :include
                try
                    func = getfield(mod, name)
                    if isa(func, Function)
                        # Try calling the function to see if it produces the right name
                        result = func()
                        if isa(result, Tuple) && length(result) >= 2
                            material_obj = result[1]
                            if hasfield(typeof(material_obj), :Name)
                                obj_name = if isa(material_obj.Name, Ptr)
                                    unsafe_string(material_obj.Name)
                                else
                                    string(material_obj.Name)
                                end
                                if obj_name == name_str
                                    return func
                                end
                            end
                        end
                    end
                catch
                    # Skip functions that can't be called or don't match
                    continue
                end
            end
        end
    catch
        # If module access fails, return nothing
    end
    return nothing
end

"""
    add_parameter_reference!(param, phase_idx, references, global_counters, category_name)

Add bibliographic reference for a parameter if available using generic dispatch.
"""
function add_parameter_reference!(param::AbstractMaterialParam, phase_idx::Int, references::Dict, global_counters::Dict, category_name::String)
    return try
        # Use the new generic function to get MaterialParamsInfo
        info = get_material_reference_info(param)

        if info !== nothing && hasfield(typeof(info), :BibTex_Reference) && !isempty(info.BibTex_Reference)
            bibtex_content = info.BibTex_Reference
            bibtex_key = extract_bibtex_key(bibtex_content)

            if !isempty(bibtex_key) && !haskey(references, bibtex_key)
                references[bibtex_key] = Dict(
                    "bibtex" => bibtex_content,
                    "number" => global_counters["reference_counter"],
                    "phase" => phase_idx,
                    "category" => category_name
                )
                global_counters["reference_counter"] += 1
            end
        end
    catch
        # Skip if reference extraction fails
    end
end

"""
    extract_parameters_from_single_param!(param, phase_idx, parameters_dict)

Extract all parameter values from a single material parameter.
"""
function extract_parameters_from_single_param!(param::AbstractMaterialParam, phase_idx::Int, parameters_dict::Dict)
    for field_name in propertynames(param)
        field_val = getproperty(param, field_name)

        if field_val isa GeoUnit
            param_symbol = string(field_name)

            # Create parameter entry if it doesn't exist
            if !haskey(parameters_dict, param_symbol)
                parameters_dict[param_symbol] = Dict(
                    "symbol" => param_symbol,
                    "latex" => create_latex_symbol(field_name),
                    "description" => get_parameter_description(param_symbol),
                    "values" => String[],
                    "phases" => Int[]
                )
            end

            # Add value for this phase
            push!(parameters_dict[param_symbol]["values"], string(field_val.val))
            push!(parameters_dict[param_symbol]["phases"], phase_idx)
        end
    end
    return
end

"""
    get_parameter_description(symbol) -> String

Get parameter description with units for a given symbol.
"""
function get_parameter_description(symbol::String)
    descriptions = get_parameter_descriptions()
    latex_symbol = create_latex_symbol(symbol)
    return get(descriptions, latex_symbol, "Parameter $(symbol)")
end

"""
    extract_parameters_from_phases(phases)

Extract all parameters from material phases into dictionaries for LaTeX table generation.
Returns (parameters_dict, references_dict).

This function processes material phase definitions and extracts all parameter values,
symbols, and metadata needed to generate parameter tables. It handles nested
composite rheologies and parallel combinations recursively.

# Arguments
- `phases`: Single phase or tuple of material phase definitions

# Returns
- `parameters_dict`: Dictionary with parameter data for table generation
- `references_dict`: Dictionary with bibliographic references

# Example
```julia
params, refs = extract_parameters_from_phases(my_phases)
```
"""
function extract_parameters_from_phases(phases)
    # Initialize output dictionaries
    params = Dict{String, Tuple{String, String, String, String, String, String}}()
    refs = Dict{String, Tuple{String, String, String}}()

    # Ensure phases is a tuple
    phases = make_tuple(phases)
    phase_count = length(phases)

    # Process each phase
    for (phase_idx, phase) in enumerate(phases)
        # Initialize counters for this phase
        counters = Dict{String, Int}()

        # Store phase name
        if hasfield(typeof(phase), :Name)
            phase_name = if isa(phase.Name, Ptr)
                unsafe_string(phase.Name)
            else
                string(phase.Name)
            end
            key = "Name $(phase_idx)"
            params[key] = (phase_name, string(phase_count), string(phase_idx), "", "", "")
        end

        # Process all fields of the phase
        for field_name in propertynames(phase)
            field_val = getproperty(phase, field_name)

            # Skip Name field (already processed) and Phase field
            if field_name == :Name || field_name == :Phase
                continue
            end

            # Skip if field is nothing or an empty tuple
            if field_val === nothing || (field_val isa Tuple && length(field_val) == 0)
                continue
            end

            # Process each element in the field (fields are usually tuples)
            if field_val isa Tuple
                for element in field_val
                    if element isa AbstractMaterialParam
                        extract_parameters_recursive(element, phase_idx, counters, params, refs, string(field_name) * "_")
                    end
                end
            elseif field_val isa AbstractMaterialParam
                # Single material parameter (not in a tuple)
                extract_parameters_recursive(field_val, phase_idx, counters, params, refs, string(field_name) * "_")
            end
        end
    end

    return params, refs
end


"""
    Dict2LatexTable(d::Dict, refs::Dict; filename="ParameterTable", rdigits=4)

Creates a LaTeX table from the parameter dictionary.
"""
function Dict2LatexTable(d::Dict, refs::Dict; filename = "ParameterTable", rdigits = 4)
    # Create sorted pairs for iteration
    dictpairs = sort(collect(pairs(d)))
    refpairs = sort(collect(pairs(refs)))

    # Get phase count
    phase_count = parse(Int64, d["Name 1"][2])

    # Parameter descriptions
    descriptions = get_parameter_descriptions()

    # Generate LaTeX preamble and table setup
    latex_content = generate_latex_preamble(phase_count)

    # Add phase headers with references
    latex_content *= generate_phase_headers(d, refpairs, phase_count)

    # Extract unique parameter symbols
    symbols = extract_unique_symbols(dictpairs)

    # Generate table rows
    for symbol in symbols
        latex_content *= generate_latex_row(symbol, descriptions, dictpairs, phase_count, rdigits)
    end

    # Add flow law equations
    latex_content *= generate_flow_law_equations(dictpairs, phase_count)

    # Add references and finish table
    references_content = generate_references(refpairs)
    latex_content *= generate_latex_footer(refpairs, phase_count)

    # Write files
    write("References.bib", references_content)
    write("$(filename).tex", latex_content)

    return nothing
end

"""
    get_parameter_descriptions() -> Dict

Returns dictionary mapping parameter symbols to their descriptions (without units).
"""
function get_parameter_descriptions()
    return Dict(
        "\\rho" => "Density",
        "\\rho_0" => "Reference density",
        "\\eta" => "Viscosity",
        "\\eta_0" => "Reference viscosity",
        "g" => "Gravity",
        "P" => "Pressure",
        "T" => "Temperature",
        "T_0" => "Reference temperature",
        "P_0" => "Reference pressure",
        "V" => "Volume",
        "d" => "Grain size",
        "f" => "Water fugacity",
        "n" => "Power-law exponent",
        "r" => "Water fugacity exponent",
        "p" => "Grain size exponent",
        "A" => "Viscosity parameter",
        "B" => "Temperature parameter",
        "A_{\\text{diff}}" => "Diffusion coefficient",
        "A_{\\text{disl}}" => "Dislocation coefficient",
        "E" => "Activation energy",
        "R" => "Gas constant",
        "G" => "Shear modulus",
        "\\nu" => "Poisson ratio",
        "K" => "Bulk modulus",
        "Kb" => "Elastic bulk modulus",
        "Y" => "Young's modulus",
        "k" => "Thermal conductivity",
        "Cp" => "Heat capacity",
        "Q_L" => "Latent heat",
        "H_r" => "Radioactive heat",
        "H_s" => "Shear heating",
        "\\phi" => "Friction angle",
        "\\psi" => "Dilation angle",
        "C" => "Cohesion",
        "Vp" => "P-wave velocity",
        "Vs" => "S-wave velocity",
        "\\beta" => "Compressibility",
        "\\alpha" => "Thermal expansion coeff."
    )
end

"""
    get_parameter_units() -> Dict

Returns dictionary mapping parameter symbols to their units.
"""
function get_parameter_units()
    return Dict(
        "\\rho" => "[kg/m\$^3\$]",
        "\\rho_0" => "[kg/m\$^3\$]",
        "\\eta" => "[Pa s]",
        "\\eta_0" => "[Pa s]",
        "g" => "[m/s\$^2\$]",
        "P" => "[MPa]",
        "T" => "[\$^\\circ\$C]",
        "T_0" => "[\$^\\circ\$C]",
        "P_0" => "[Pa]",
        "V" => "[m\$^3\$]",
        "d" => "[cm]",
        "f" => "[MPa]",
        "n" => "[-]",
        "r" => "[-]",
        "p" => "[-]",
        "A" => "[-]",
        "B" => "[K]",
        "A_{\\text{diff}}" => "[Pa\$^{-n-r}\$ m\$^p\$/s]",
        "A_{\\text{disl}}" => "[Pa\$^{-n}\$/s]",
        "E" => "[J/mol]",
        "R" => "[J/mol/K]",
        "G" => "[Pa]",
        "\\nu" => "[-]",
        "K" => "[Pa]",
        "Kb" => "[Pa]",
        "Y" => "[Pa]",
        "k" => "[W/m/K]",
        "Cp" => "[J/kg/K]",
        "Q_L" => "[kJ/kg]",
        "H_r" => "[W/m\$^3\$]",
        "H_s" => "[-]",
        "\\phi" => "[\$^\\circ\$]",
        "\\psi" => "[\$^\\circ\$]",
        "C" => "[Pa]",
        "Vp" => "[km/s]",
        "Vs" => "[km/s]",
        "\\beta" => "[1/Pa]",
        "\\alpha" => "[1/K]"
    )
end

function generate_latex_preamble(phase_count::Int)
    column_spec = "l" * "c" * "l" * "c"^phase_count  # description + symbol + units + phases
    return """
    \\documentclass{article}
    \\usepackage[utf8]{inputenc}
    \\usepackage{booktabs}
    \\usepackage{graphicx}
    \\usepackage{multirow}
    \\usepackage[round, comma, sort]{natbib}
    \\bibliographystyle{abbrvnat}
    \\setcitestyle{authoryear,open={(},close={)}}
    \\begin{document}
    \\begin{table}[hbt]
    \\centering
    \\resizebox{\\columnwidth}{!}{
    \\begin{tabular}{$(column_spec)}
    \\toprule[1pt]
    """
end

function generate_phase_headers(d::Dict, refpairs, phase_count::Int)
    header = "Denotation & Variable & Units"
    in_text_ref = ""
    counter = 1

    for i in 1:phase_count
        header *= " & " * d["Name $(i)"][1]

        # Add references
        for (_, ref_data) in refpairs
            if parse(Int64, ref_data[3]) == i
                header *= "(" * "*"^counter * ")"
                bib_key = extract_bibtex_key(ref_data[1])
                if !isempty(bib_key)
                    in_text_ref *= "(" * "*"^counter * ") \\cite{$(bib_key)}"
                    if counter < length(refpairs)
                        in_text_ref *= ", "
                    end
                end
                counter += 1
            end
        end
    end

    return header * " \\\\\n\\midrule\n"
end

function extract_bibtex_key(bibtex::String)
    start_idx = findfirst("{", bibtex)
    end_idx = findfirst(",", bibtex)
    if start_idx !== nothing && end_idx !== nothing
        return bibtex[(start_idx[1] + 1):(end_idx[1] - 1)]
    end
    return ""
end

function extract_unique_symbols(dictpairs)
    symbols = String[]
    for (key, data) in dictpairs
        if !occursin("Name", key)
            push!(symbols, data[2])  # Symbol is in second position
        end
    end
    return unique(sort(symbols))
end

function generate_latex_row(symbol::String, descriptions::Dict, dictpairs, phase_count::Int, rdigits::Int)
    # Get description and units
    desc = get(descriptions, symbol, "Unknown parameter")
    units_dict = get_parameter_units()
    units = get(units_dict, symbol, "[-]")

    # Start row with description, symbol, and units
    row = "$(desc) & \$$(symbol)\$ & $(units)"

    # Add values for each phase
    for phase_idx in 1:phase_count
        value_found = false
        for (key, data) in dictpairs
            if !occursin("Name", key) && data[2] == symbol && parse(Int64, data[4]) == phase_idx
                dig, num, expo = detachFloatfromExponent(data[1])
                formatted_value = format_number(num, expo, dig, rdigits)
                row *= " & \$$(formatted_value)\$"
                value_found = true
                break
            end
        end
        if !value_found
            row *= " & "
        end
    end

    return row * " \\\\\n"
end

function format_number(num::String, expo::String, dig::Int, rdigits::Int)
    if dig <= rdigits && expo != "1"
        return "$(num) \\times 10^{$(expo)}"
    elseif dig <= rdigits && expo == "1"
        return num
    elseif dig > rdigits && expo != "1"
        rounded = round(parse(Float64, num); digits = rdigits)
        return "$(rounded) \\times 10^{$(expo)}"
    else
        rounded = round(parse(Float64, num); digits = rdigits)
        return string(rounded)
    end
end

function generate_flow_law_equations(dictpairs, phase_count::Int)
    equations = ""
    added_disl = false
    added_diff = false
    added_lin = false
    total_columns = phase_count + 3  # description + symbol + units + phases

    for (_, data) in dictpairs
        if data[6] == "DislCreep" && !added_disl
            equations *= "Dislocation Creep: & \\multicolumn{$(total_columns - 1)}{l}{\$ \\dot{\\gamma} = A \\tau^n f_{H2O}^r \\exp(-\\frac{E+PV}{RT}) \$}\\\\\n"
            added_disl = true
        elseif data[6] == "DiffCreep" && !added_diff
            equations *= "Diffusion Creep: & \\multicolumn{$(total_columns - 1)}{l}{\$ \\dot{\\gamma} = A \\tau^n d^p f_{H2O}^r \\exp(-\\frac{E+PV}{RT}) \$}\\\\\n"
            added_diff = true
        elseif data[6] == "LinVisc" && !added_lin
            equations *= "Linear viscosity: & \\multicolumn{$(total_columns - 1)}{l}{\$ \\eta = \\frac{\\tau_{II}}{2\\dot{\\varepsilon_{II}}} \$}\\\\\n"
            added_lin = true
        end
    end

    return equations
end

function generate_references(refpairs)
    references = ""
    for (_, ref_data) in refpairs
        references *= ref_data[1]  # BibTeX content
    end
    return references
end

function generate_latex_footer(refpairs, phase_count::Int)
    in_text_refs = ""
    counter = 1
    total_columns = phase_count + 3  # description + symbol + units + phases

    for (_, ref_data) in refpairs
        bib_key = extract_bibtex_key(ref_data[1])
        if !isempty(bib_key)
            in_text_refs *= "(" * "*"^counter * ") \\cite{$(bib_key)}"
            if counter < length(refpairs)
                in_text_refs *= ", "
            end
            counter += 1
        end
    end

    footer = ""
    if !isempty(in_text_refs)
        footer *= "\\midrule[0.3pt]\n"
        footer *= "\\multicolumn{$(total_columns)}{l}{$(in_text_refs)}\\\\\n"
    end

    footer *= """
    \\bottomrule[1pt]
    \\end{tabular}
    }
    \\caption{Parameter table}
    \\label{tab:para_table}
    \\end{table}
    \\bibliography{References}
    \\end{document}
    """

    return footer
end

"""
    extract_parameters_from_phases_md(phases) -> Dict

Extract all parameters from material phases into dictionary for Markdown table generation.
Returns parameters dictionary suitable for markdown output.

This function is similar to extract_parameters_from_phases() but preserves Unicode
symbols instead of converting them to LaTeX format, making it suitable for
Markdown table generation.

# Arguments
- `phases`: Single phase or tuple of material phase definitions

# Returns
- `parameters_dict`: Dictionary with parameter data for Markdown table generation

# Example
```julia
params = extract_parameters_from_phases_md(my_phases)
```
"""
function extract_parameters_from_phases_md(phases)
    # For Markdown, we want to preserve Unicode symbols instead of converting to LaTeX
    # Initialize output dictionaries
    params = Dict{String, Tuple{String, String, String, String, String, String}}()

    # Ensure phases is a tuple
    phases = make_tuple(phases)
    phase_count = length(phases)

    # Process each phase
    for (phase_idx, phase) in enumerate(phases)
        # Initialize counters for this phase
        counters = Dict{String, Int}()

        # Store phase name
        if hasfield(typeof(phase), :Name)
            phase_name = if isa(phase.Name, Ptr)
                unsafe_string(phase.Name)
            else
                string(phase.Name)
            end
            key = "Name $(phase_idx)"
            params[key] = (phase_name, string(phase_count), string(phase_idx), "", "", "")
        end

        # Process all fields of the phase
        for field_name in propertynames(phase)
            field_val = getproperty(phase, field_name)

            # Skip Name field (already processed) and Phase field
            if field_name == :Name || field_name == :Phase
                continue
            end

            # Skip if field is nothing or an empty tuple
            if field_val === nothing || (field_val isa Tuple && length(field_val) == 0)
                continue
            end

            # Process each element in the field (fields are usually tuples)
            if field_val isa Tuple
                for element in field_val
                    if element isa AbstractMaterialParam
                        extract_parameters_recursive_md(element, phase_idx, counters, params, string(field_name) * "_")
                    end
                end
            elseif field_val isa AbstractMaterialParam
                # Single material parameter (not in a tuple)
                extract_parameters_recursive_md(field_val, phase_idx, counters, params, string(field_name) * "_")
            end
        end
    end

    return params
end

"""
    extract_parameters_recursive_md(obj, phase_idx::Int, counters::Dict, params::Dict, prefix::String)

Recursively extract parameters for Markdown format (preserves Unicode symbols).
"""
function extract_parameters_recursive_md(obj, phase_idx::Int, counters::Dict, params::Dict, prefix::String = "")
    return if obj isa CompositeRheology
        # Build composite rheology summary
        composite_summary = build_composite_summary(obj)
        new_prefix = prefix * "Comp_"
        for (i, element) in enumerate(obj.elements)
            extract_parameters_recursive_md_with_composite(element, phase_idx, counters, params, new_prefix * "$(i)_", composite_summary)
        end
    elseif obj isa Parallel
        new_prefix = prefix * "Para_"
        for (i, element) in enumerate(obj.elements)
            extract_parameters_recursive_md(element, phase_idx, counters, params, new_prefix * "$(i)_")
        end
    elseif obj isa AbstractMaterialParam
        # Base case - single material parameter
        extract_parameters_from_material_md(obj, phase_idx, counters, params, prefix)
    end
end

"""
    extract_parameters_recursive_md_with_composite(obj, phase_idx, counters, params, prefix, composite_summary)

Extract parameters with composite rheology summary information for Markdown format.
"""
function extract_parameters_recursive_md_with_composite(obj, phase_idx::Int, counters::Dict, params::Dict, prefix::String, composite_summary::String)
    return if obj isa AbstractMaterialParam
        # Extract parameters with composite summary
        extract_parameters_from_material_md_with_composite(obj, phase_idx, counters, params, prefix, composite_summary)
    end
end

"""
    extract_parameters_from_material_md(material, phase_idx::Int, counters::Dict, params::Dict, prefix::String)

Extract parameters from a single material parameter for Markdown format (preserves Unicode).
"""
function extract_parameters_from_material_md(material::AbstractMaterialParam, phase_idx::Int, counters::Dict, params::Dict, prefix::String = "")
    rheo_type, rheo_subtype = get_rheology_info(material)

    # Only treat actual rheology laws as rheology types, not basic material properties
    is_rheology = material isa Union{
        DislocationCreep, DiffusionCreep, LinearViscous, PowerlawViscous,
        GrainBoundarySliding, PeierlsCreep, NonLinearPeierlsCreep,
        LinearMeltViscosity, ViscosityPartialMelt_Costa_etal_2009, GiordanoMeltViscosity,
    }

    # Update counters if this is a rheology law
    if is_rheology
        counters[rheo_type] = get(counters, rheo_type, 0) + 1
        counters["flowlaw"] = get(counters, "flowlaw", 0) + 1
    else
        # For non-rheology parameters, use empty strings
        rheo_type = ""
        rheo_subtype = ""
    end

    # Extract parameter values
    for field_name in propertynames(material)
        field_val = getproperty(material, field_name)

        if field_val isa GeoUnit
            value_str = string(field_val.val)
            # For Markdown, preserve the original Unicode symbol
            symbol_str = string(field_name)

            # Generate key in expected format
            # For composite rheologies, use special formatting: "fieldname FieldType RheoType phase.element"
            if !isempty(prefix) && occursin("Comp_", prefix)
                # Extract field name and element index from prefix
                # e.g., "CompositeRheology_Comp_3_" -> field="CompositeRheology", element=3
                prefix_parts = split(strip(prefix, '_'), "_")
                field_type = prefix_parts[1]  # "CompositeRheology"
                element_idx = prefix_parts[end]  # "3"

                # Use rheology type and format as phase.element
                key = "$(field_name) $(field_type) $(rheo_subtype) $(phase_idx).$(element_idx)"
            else
                # Standard format for non-composite rheologies
                field_type = if !isempty(prefix) && endswith(prefix, "_")
                    strip(prefix, '_')
                else
                    "param"
                end
                key = "$(field_name) $(field_type) $(phase_idx)"
            end

            params[key] = (
                value_str,
                symbol_str,  # Keep original Unicode symbol for Markdown
                rheo_type,
                string(phase_idx),
                isempty(rheo_type) ? "" : string(get(counters, rheo_type, 0)),
                rheo_subtype,
            )
        end
    end
    return
end

"""
    extract_parameters_from_material_md_with_composite(material, phase_idx, counters, params, prefix, composite_summary)

Extract parameters from a single material parameter for Markdown format with composite rheology information.
"""
function extract_parameters_from_material_md_with_composite(material::AbstractMaterialParam, phase_idx::Int, counters::Dict, params::Dict, prefix::String, composite_summary::String)
    rheo_type, rheo_subtype = get_rheology_info(material)

    # Only treat actual rheology laws as rheology types, not basic material properties
    is_rheology = material isa Union{
        DislocationCreep, DiffusionCreep, LinearViscous, PowerlawViscous,
        GrainBoundarySliding, PeierlsCreep, NonLinearPeierlsCreep,
        LinearMeltViscosity, ViscosityPartialMelt_Costa_etal_2009, GiordanoMeltViscosity,
    }

    # Update counters if this is a rheology law
    if is_rheology
        counters[rheo_type] = get(counters, rheo_type, 0) + 1
        counters["flowlaw"] = get(counters, "flowlaw", 0) + 1
    else
        # For non-rheology parameters, use empty strings
        rheo_type = ""
        rheo_subtype = ""
    end

    # Extract parameter values
    for field_name in propertynames(material)
        field_val = getproperty(material, field_name)

        if field_val isa GeoUnit
            value_str = string(field_val.val)
            # For Markdown, preserve the original Unicode symbol
            symbol_str = string(field_name)

            # Generate key in expected format
            # For composite rheologies, use special formatting: "fieldname FieldType RheoType phase.element"
            if !isempty(prefix) && occursin("Comp_", prefix)
                # Extract field name and element index from prefix
                # e.g., "CompositeRheology_Comp_3_" -> field="CompositeRheology", element=3
                prefix_parts = split(strip(prefix, '_'), "_")
                field_type = prefix_parts[1]  # "CompositeRheology"
                element_idx = prefix_parts[end]  # "3"

                # Use rheology type and format as phase.element
                key = "$(field_name) $(field_type) $(rheo_subtype) $(phase_idx).$(element_idx)"
            else
                # Standard format for non-composite rheologies
                field_type = if !isempty(prefix) && endswith(prefix, "_")
                    strip(prefix, '_')
                else
                    "param"
                end
                key = "$(field_name) $(field_type) $(phase_idx)"
            end

            params[key] = (
                value_str,
                symbol_str,  # Keep original Unicode symbol for Markdown
                composite_summary,  # Use composite summary instead of individual rheo_type
                string(phase_idx),
                isempty(rheo_type) ? "" : string(get(counters, rheo_type, 0)),
                rheo_subtype,
            )
        end
    end
    return
end

"""
    Dict2MarkdownTable(d::Dict; filename="ParameterTable", rdigits=4)

Creates a Markdown table from the parameter dictionary.
"""
function Dict2MarkdownTable(d::Dict; filename = "ParameterTable", rdigits = 4)
    # Get phase count
    if !haskey(d, "Name 1")
        error("Invalid parameter dictionary: missing phase information")
    end

    phase_count = parse(Int64, d["Name 1"][2])

    # Parameter descriptions for markdown
    descriptions = get_markdown_descriptions()
    units_dict = get_markdown_units()

    # Create header
    table_content = " | Denotation | Variable | Units | "
    for i in 1:phase_count
        table_content *= d["Name $(i)"][1] * " | "
    end
    table_content *= "\n"

    # Create separator
    table_content *= " | ---------- | -------- | ----- | "
    for i in 1:phase_count
        table_content *= "-"^length(d["Name $(i)"][1]) * " | "
    end
    table_content *= "\n"

    # Extract unique symbols
    symbols = String[]
    for (key, data) in d
        if !occursin("Name", key)
            push!(symbols, data[2])
        end
    end
    symbols = unique(sort(symbols))

    # Generate table rows
    for symbol in symbols
        desc = get(descriptions, symbol, "Unknown parameter")
        units = get(units_dict, symbol, "[-]")
        table_content *= " | $(desc) | $(symbol) | $(units) | "

        # Add values for each phase
        for phase_idx in 1:phase_count
            value_found = false
            for (key, data) in d
                if !occursin("Name", key) && data[2] == symbol && parse(Int64, data[4]) == phase_idx
                    dig, num, expo = detachFloatfromExponent(data[1])
                    formatted_value = format_markdown_number(num, expo, dig, rdigits)
                    table_content *= formatted_value * " | "
                    value_found = true
                    break
                end
            end
            if !value_found
                table_content *= " | "
            end
        end
        table_content *= "\n"
    end

    # Write file
    write("$(filename).md", table_content)
    return nothing
end

"""
    get_markdown_descriptions() -> Dict

Returns dictionary mapping parameter symbols to their markdown descriptions (without units).
"""
function get_markdown_descriptions()
    return Dict(
        "ρ" => "Density",
        "ρ0" => "Reference density",
        "g" => "Gravity",
        "η" => "Viscosity",
        "P" => "Pressure",
        "T" => "Temperature",
        "V" => "Volume",
        "d" => "Grain size",
        "f" => "Water fugacity",
        "n" => "Power-law exponent",
        "r" => "Water fugacity exponent",
        "p" => "Grain size exponent",
        "A_diff" => "Coefficient",
        "A_disl" => "Coefficient",
        "E" => "Activation energy",
        "R" => "Gas constant",
        "G" => "Shear modulus",
        "ν" => "Poisson ratio",
        "K" => "Bulk modulus",
        "Kb" => "Elastic bulk modulus",
        "Y" => "Young's modulus",
        "k" => "Thermal conductivity",
        "Cp" => "Heat capacity",
        "Q_L" => "Latent heat",
        "H_r" => "Radioactive heat",
        "H_s" => "Shear heating",
        "ϕ" => "Friction angle",
        "ψ" => "Dilation angle",
        "C" => "Cohesion",
        "Vp" => "P-wave velocity",
        "Vs" => "S-wave velocity",
        "T0" => "Reference temperature",
        "P0" => "Reference pressure",
        "β" => "Compressibility",
        "α" => "Thermal expansion coeff."
    )
end

"""
    get_markdown_units() -> Dict

Returns dictionary mapping parameter symbols to their markdown units.
"""
function get_markdown_units()
    return Dict(
        "ρ" => "[kg/m³]",
        "ρ0" => "[kg/m³]",
        "g" => "[m/s²]",
        "η" => "[Pa·s]",
        "P" => "[MPa]",
        "T" => "[°C]",
        "V" => "[m³]",
        "d" => "[cm]",
        "f" => "[MPa]",
        "n" => "[-]",
        "r" => "[-]",
        "p" => "[-]",
        "A_diff" => "[Pa⁻ⁿ⁻ʳ·mᵖ/s]",
        "A_disl" => "[Pa⁻ⁿ/s]",
        "E" => "[J/mol]",
        "R" => "[J/mol/K]",
        "G" => "[Pa]",
        "ν" => "[-]",
        "K" => "[Pa]",
        "Kb" => "[Pa]",
        "Y" => "[Pa]",
        "k" => "[W/m/K]",
        "Cp" => "[J/kg/K]",
        "Q_L" => "[kJ/kg]",
        "H_r" => "[W/m³]",
        "H_s" => "[-]",
        "ϕ" => "[°]",
        "ψ" => "[°]",
        "C" => "[Pa]",
        "Vp" => "[km/s]",
        "Vs" => "[km/s]",
        "T0" => "[°C]",
        "P0" => "[Pa]",
        "β" => "[1/Pa]",
        "α" => "[1/K]"
    )
end

function format_markdown_number(num::String, expo::String, dig::Int, rdigits::Int)
    if dig <= rdigits && expo != "1"
        return "$(num) × 10^$(expo)"
    elseif dig <= rdigits && expo == "1"
        return num
    elseif dig > rdigits && expo != "1"
        rounded = round(parse(Float64, num); digits = rdigits)
        return "$(rounded) × 10^$(expo)"
    else
        rounded = round(parse(Float64, num); digits = rdigits)
        return string(rounded)
    end
end

"""
    ParameterTable(phases; filename="ParameterTable", format="latex", rdigits=4)

Creates a table with all parameters from the material phases.
Supports both LaTeX and Markdown formats.

# Arguments
- `phases`: Material parameter phases (single phase or tuple of phases)
- `filename`: Output filename (without extension)
- `format`: Output format ("latex"/"tex" or "markdown"/"md")
- `rdigits`: Number of decimal places for rounding

# Examples
```julia
# LaTeX table
ParameterTable(MatParam; format="latex", filename="MyTable")

# Markdown table
ParameterTable(MatParam; format="markdown", filename="MyTable")
```
"""
function ParameterTable(phases; filename = "ParameterTable", format = "latex", rdigits = 4)
    format_lower = lowercase(format)

    if format_lower in ("latex", "tex")
        d, refs = extract_parameters_from_phases(phases)
        Dict2LatexTable(d, refs; filename = filename, rdigits = rdigits)
        @info "Created $(filename).tex and References.bib files."
    elseif format_lower in ("markdown", "md")
        d = extract_parameters_from_phases_md(phases)
        Dict2MarkdownTable(d; filename = filename, rdigits = rdigits)
        @info "Created $(filename).md file."
    else
        error("Unsupported format: $(format). Use 'latex' or 'markdown'.")
    end

    return nothing
end

end
