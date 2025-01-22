# Pretty printing for CompositeRheologies

# returns a matrix with strings in the right order
function print_rheology_matrix(v::Tuple)

    n = 40
    A = Matrix{String}(undef, n, n)

    i, j, i_max = 1, 0, 1
    for entry in eachindex(v)
        out = print_rheology_matrix(v[entry])
        si = size(out)
        if prod(si) == 1
            j = j + 1
            A[i, j] = out[1]
        elseif (length(si) == 1) && prod(si) != 1
            j = j + 1
            A[i:(i + si[1] - 1), j] = out
            i_max = max(i_max, si[1])
        else
            A[i:(i + si[1] - 1), j:(j + si[2] - 1)] = out
            j = j + si[2] - 1
        end
    end
    for i in eachindex(A)
        if !isassigned(A, i)
            A[i] = ""
        end
    end
    for i in eachindex(A)
        if A[i] == ""
            A[i] = print_rheology_matrix("")[1]
        end
    end

    A = A[1:i_max, 1:j]

    return A
end


function print_rheology_matrix(v::Parallel)

    n = 40
    A = Matrix{String}(undef, n, n)
    elements = v.elements
    i, j = 1, 1
    i_vec = Int64[]
    for entry in eachindex(elements)
        out = print_rheology_matrix(elements[entry])
        si = size(out)
        if prod(si) == 1
            push!(i_vec, i)
            A[i, j] = out[1]
            i = i + 1
        elseif (length(si) == 1) && prod(si) != 1
            push!(i_vec, i)
            A[i:(i + si[1] - 1), j] = out
            i = i + si[1]
        else
            push!(i_vec, i)
            A[i:(i + si[1] - 1), j:(j + si[2] - 1)] = out
            i, j = i + si[1], j + si[2] - 1
        end
    end

    # Fill every index (with empty strings)
    for i in eachindex(A)
        if !isassigned(A, i)
            A[i] = ""
        end
    end

    # Extract the relevant part of A
    A = A[1:i, 1:j]

    # Center the strings & put brackets around it
    B = create_string_vec(A)

    nel = maximum(textwidth.(B))
    for i in 1:length(B)
        if any(in.(i_vec, i))
            B[i] = cpad(B[i], nel, "-")
            B[i] = "|" * B[i] * "|"
        else
            B[i] = cpad(B[i], nel, " ")
            B[i] = "|" * B[i] * "|"
        end
    end

    return B
end

function create_string_vec(A)
    B = String[]
    for i in 1:size(A, 1)
        str1 = join(A[i, :])
        if length(str1) > 0
            push!(B, str1)
        end
    end

    return B
end

function print_rheology_matrix(v::CompositeRheology)
    n = 40
    A = Matrix{String}(undef, n, n)
    elements = v.elements
    i, j, i_max = 1, 0, 1
    for entry in eachindex(elements)
        out = print_rheology_matrix(elements[entry])
        si = size(out)
        if prod(si) == 1
            j = j + 1
            A[i, j] = out[1]
        elseif (length(si) == 1) && prod(si) != 1
            j = j + 1
            A[i:(i + si[1] - 1), j] = out
            i_max = max(i_max, si[1])
        else
            j = j + 1
            A[i:(i + si[1] - 1), j:j] = out
            i = i + si[1]

            i_max = max(i_max, si[1])
        end
    end
    # Fill every index (with empty strings)
    for i in eachindex(A)
        if !isassigned(A, i)
            A[i] = ""
        end
    end

    for i in eachindex(A)
        if A[i] == ""
            A[i] = print_rheology_matrix("")[1]
        end
    end

    A = A[1:i_max, 1:j]             # Extract the relevant part of A
    B = create_string_vec(A)        # Create strings

    return B
end

# Print the individual rheological elements in the REPL
print_rheology_matrix(v::String) = ["         "]
print_rheology_matrix(v::AbstractCreepLaw) = ["--⟦▪̲̅▫̲̅▫̲̅▫̲̅--"]
print_rheology_matrix(v::CustomRheology) = ["--?????--"]
print_rheology_matrix(v::AbstractElasticity) = ["--/\\/\\/--"]
print_rheology_matrix(v::AbstractPlasticity) = ["--▬▬▬__--"]
#print_rheology_matrix(v::DruckerPrager)      = ["-dp▬▬__--"] # we can further


function create_rheology_string(str, rheo_Comp::CompositeRheology)

    rheology = rheo_Comp.elements
    for i in eachindex(rheology)
        str = create_rheology_string(str, rheology[i])
    end

    return str
end

function create_rheology_string(str, rheo_Parallel::Parallel)
    rheology = rheo_Parallel.elements
    str = str * "{"
    for i in eachindex(rheology)
        str = create_rheology_string(str, rheology[i])
        if str[end] == 'o'
            str = str[1:(end - 1)]
        end
        str = str * ";"
    end
    str = str[1:(end - 1)] * "}"      # removes the last ";"

    return str
end

function create_rheology_string(str, rheology::Tuple)
    for i in eachindex(rheology)
        str = create_rheology_string(str, rheology[i])
        str = str * "o"
    end
    return str
end

# Print the individual rheological elements in the REPL
create_rheology_string(str, rheo_Parallel::AbstractCreepLaw) = str = str * "--⟦▪̲̅▫̲̅▫̲̅▫̲̅--"
create_rheology_string(str, rheo_Parallel::AbstractPlasticity) = str = str * "--▬▬▬__--"
create_rheology_string(str, rheo_Parallel::AbstractElasticity) = str = str * "--/\\/\\/--"


function create_parallel_str(str)
    # Print them underneath each other:
    l_start = findfirst("{", str)
    l_end = findlast("}", str)

    if !isnothing(l_start)

        # step 1: all inner Parallel objects should be left untouched
        l_st2 = findnext("{", str, l_start[1] + 1)
        l_end2 = findprev("}", str, l_end[1] - 1)

        @show l_st2, l_end2
        if !isnothing(l_st2)
            str_sub = str[(l_st2[1] + 1):(l_end2[1] - 1)]
            str_sub = replace(str_sub, ";" => "X")

            str = str[1:l_st2[1]] * str_sub * str[l_end2[1]:end]
        end


        str1 = split(str[(l_start[1] + 1):(l_end[1] - 2)], ";")
        len = maximum(textwidth.(str1))
        for i in eachindex(str1)
            str1[i] = "|" * cpad(str1[i], len, "-") * "|\n"
        end
        str_out = join(str1)     # join the vectors back together
        str_out = str_out[1:(end - 1)]
        str_out = replace(str_out, "X" => ";")

        str_out1 = str[1:(l_start[1] - 2)] * str_out * str[(l_end[1] + 2):end]

        str2 = split(str_out1, "\n")
        len = maximum(textwidth.(str2))
        for i in eachindex(str2)
            str2[i] = cpad(str2[i], len, " ") * "\n"
        end
        str_out = join(str2)

    else
        str_out = str
    end

    return str_out
end


# Center strings
cpad(s, n::Integer, p = " ") = rpad(lpad(s, div(n + textwidth(s), 2), p), n, p)
struct InverseCreepLaw{N} <: AbstractConstitutiveLaw{Float64}
    v::NTuple{N, AbstractConstitutiveLaw}

    function InverseCreepLaw(v::Vararg{AbstractConstitutiveLaw, N}) where {N}
        return new{N}(ntuple(i -> v[i], Val(N)))
    end

    function InverseCreepLaw(v::NTuple{N, AbstractConstitutiveLaw}) where {N}
        return new{N}(v)
    end
end
