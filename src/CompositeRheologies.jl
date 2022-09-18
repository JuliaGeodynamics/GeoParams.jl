# This holds structures and computational routines for compositional rheologies


export CompositeRheology, Parallel, create_rheology_string, pretty_print_rheology_string, print_rheology_matrix!
#import Base: size

"""
    Structure that holds composite rheologies (e.g., visco-elasto-viscoplastic),
    but also indicates (in the name) whether we need to perform non-linear iterations.
"""
struct CompositeRheology{T, τ_it, P_it, λ_it} <: AbstractComposite
    rheology_chain::T
end

# Defines tuples of composite rheologies, while also checking which type of iterations need to be performed
function CompositeRheology(v::T) where T
    τ_it = false;
    P_it = false;
    λ_it = false;

    return CompositeRheology{typeof(v),τ_it, P_it, λ_it}(v)
end
CompositeRheology(a,b...) = CompositeRheology( (a,b...,)) 
#CompositeRheology(v::Tuple) =  CompositeRheology(v...) 

# Print info 
function show(io::IO, g::AbstractComposite)
   # rheology = g.rheology_chain;
    println(io,"Composite rheology:   ")

    # Compose a string with rheological elements, so we have an overview in the REPL
    str = create_rheology_string("",g)

    println(io,str)
 
    return nothing
end

#=
function size(v::AbstractComposite)
    # Returns the # of parallel elements & the maximu width (in case one of them is a tuple )
    le = 0
    wi = 1
    @show le
    for elem in v.rheology_chain
        if isa(elem, Tuple)
            println("boris")
            local_le = length(elem)
            @show local_le

            local_size = (0,0)
            for elem_local in elem
                if isa(elem_local, Tuple)

                end
            end

            #=
            # If one of the elements is a Parallel element again, determine the size of that one
            # We do this recursively until 
            if any(isa.(elem, Parallel))
                local_size = (0,0)
                for elem_local in elem
                    if isa(elem_local, Parallel)
                        local_size =.+ size(elem_local) .- (1,1)
                    else
                        local_size = (0,1)
                    end

                    le = le + local_size[2]+1;
                    wi = wi + local_size[1];
                end
                @show local_size
            end
            =#

            if local_le>le
                le += local_le
            end

        elseif isa(elem, Parallel)
            local_size =.+ size(elem) .- (1,1)
            le = le + local_size[2]+1;
            if local_size[1]>=wi
                wi += local_size[1];
            end

        else
            le += 1;
            
        end
    end
    return (wi, le)
end

# Computes the local size 
function compute_rheology_size(v::Tuple)
    local_size = (0,0)
    for elem_local in elem
        local_size = local_size .+ compute_rheology_size(v)
    end
    return local_size
end

compute_rheology_size(v::Parallel) = size(v)
compute_rheology_size(v::Any) = (0,1)
=#


"""
    Put rheological elements in parallel 
"""
struct Parallel{T}
    elements::T
end
Parallel(v...) = Parallel{typeof( (v...,))}((v...,))

#=
function size(v::Parallel)
    # Returns the # of parallel elements & the maximu width (in case one of them is a tuple )
    le = length(v.elements)
    wi = 1
    for elem in v.elements
        if isa(elem, Tuple)
            local_wi = length(elem)

            # If one of the elements is a Parallel element again, determine the size of that one
            # We do this recursively until 
            if any(isa.(elem, Parallel))
                local_size = (0,0)
                for elem_local in elem
                    if isa(elem_local, Parallel)
                        local_size =.+ size(elem_local) .- (1,1)
                    end
                end
                le = le + local_size[2];
                local_wi = local_wi + local_size[1];
            end

            if local_wi>wi
                wi = local_wi
            end

        elseif isa(elem, Parallel)
            local_size =.+ size(elem) .- (1,1)
            le = le + local_size[1];
            if local_size[2]>=wi
                wi += local_size[2];
            end
        end

    end
    return (le, wi)
end
=#

function show(io::IO, a::Parallel)
    rheology = a.elements;
    println(io,"Parallel:   ")  

    # Compose a string with rheological elements, so we have an overview in the REPL
    str = create_rheology_string("",a)
    

    A = Matrix{String}(undef, 6, 6)
    print_rheology_matrix!(A,1,1,a)    


    for i=1:size(A,1)
        le = 0;
        for j=1:size(A,2)
            if  isassigned(A,i,j)
                le = j;
            end
        end
        
        # fill them with space if they are not assigned yet
        for j=1:le
            if  !isassigned(A,i,j)
                A[i,j] = "         ";
            end
        end
    end

    for i in eachindex(A)
        if  !isassigned(A,i); A[i] = ""; end
    end
    
    B = String[];
    for i=1:size(A,1)
        str1 = join(A[i,:]);
        if length(str1)>0
            push!(B, str1)
        end
    end
    nel =  maximum(textwidth.(B));
    for i=1:length(B)
        if  (B[i][1]=='|') &&  (B[i][end]=='|') 
            B[i] = "|"*cpad(B[i][2:end-1],nel,"-")*"|"
        else
            B[i] = cpad(B[i],nel)
            B[i] = "|"*B[i]*"|"
        end

    end
    for i=1:length(B)
        println(B[i])
    end


    println("")
    println(io,str)

 
    return nothing
end

function print_rheology_matrix!(A::Matrix,i::Int64,j::Int64,v::Parallel)
    i_start = i
    elem  = v.elements
    i_max = length(elem) 
    for entry in eachindex(elem)
        if !isassigned(A,i,j) A[i,j]=""  end
        A[i,j] *= "|"

        j_start = j
        i,j, i_max1 = print_rheology_matrix!(A,i,j, elem[entry]) 
        if isa(elem[entry], Tuple)

            if !isassigned(A,i,j-1) A[i,j-1]=""  end
            A[i,j-1] *= "|"
            j = j_start;

        else
            A[i,j] *= "|"
        end

        i += i_max1
    end
    i = i_start

    return i,j, i_max
end

function print_rheology_matrix!(A::Matrix,i::Int64,j::Int64,v::Tuple)

    i_max = 1;
    for entry in eachindex(v)
        i,j, i_max1 = print_rheology_matrix!(A,i,j, v[entry]) 
        j += 1
        i_max = max(i_max,i_max1)
    end

    return i,j, i_max
end

# Print the individual rheological elements in the REPL
function print_rheology_matrix!(A::Matrix,i::Int64,j::Int64, v::AbstractCreepLaw)   
    if !isassigned(A,i,j) A[i,j]=""  end
    A[i,j] *= "--⟦▪̲̅▫̲̅▫̲̅▫̲̅--"
    return i,j, 1
end
function print_rheology_matrix!(A::Matrix,i::Int64,j::Int64, v::AbstractPlasticity) 
    if !isassigned(A,i,j) A[i,j]=""  end
    A[i,j] *= "--▬▬▬__--"    
    return i,j, 1
end

function print_rheology_matrix!(A::Matrix,i::Int64,j::Int64, v::AbstractElasticity) 
    if !isassigned(A,i,j) A[i,j]=""  end
    A[i,j] *= "--/\\/\\/--"
    return i,j, 1
end


function create_rheology_string(str, rheo_Comp::CompositeRheology)
 
    rheology = rheo_Comp.rheology_chain
    for i in eachindex(rheology)
        str = create_rheology_string(str,rheology[i])
    end

    return str
end

function create_rheology_string(str, rheo_Parallel::Parallel)
    rheology = rheo_Parallel.elements
    str = str*"{"
    for i in eachindex(rheology)
        str = create_rheology_string(str,rheology[i])
        if str[end]=='o'; str=str[1:end-1] end
        str = str*";"
    end
    str = str[1:end-1]*"}"      # removes the last ";"

    return str
end

function create_rheology_string(str, rheology::Tuple)
    for i in eachindex(rheology)
        str = create_rheology_string(str,rheology[i])
        str = str*"o"
    end
    return str
end

# Print the individual rheological elements in the REPL
create_rheology_string(str, rheo_Parallel::AbstractCreepLaw)   = str = str*"--⟦▪̲̅▫̲̅▫̲̅▫̲̅--"
create_rheology_string(str, rheo_Parallel::AbstractPlasticity) = str = str*"--▬▬▬__--"    
create_rheology_string(str, rheo_Parallel::AbstractElasticity) = str = str*"--/\\/\\/--"


function create_parallel_str(str)
    # Print them underneath each other:
    l_start = findfirst("{", str)
    l_end   = findlast("}", str)

    if !isnothing(l_start)

        # step 1: all inner Parallel objects should be left untouched
        l_st2  = findnext("{", str, l_start[1]+1)
        l_end2 = findprev("}", str, l_end[1]-1)
        
        @show l_st2, l_end2
        if !isnothing(l_st2)
            str_sub = str[l_st2[1]+1:l_end2[1]-1];
            str_sub = replace(str_sub,";"=>"X")
           
            str = str[1:l_st2[1]]*str_sub*str[l_end2[1]:end]
        end


        str1 = split(str[l_start[1]+1:l_end[1]-2],";")
        len  = maximum(textwidth.(str1))
        for i = eachindex(str1)
            str1[i] = "|"*cpad(str1[i],len,"-")*"|\n"
        end 
        str_out = join(str1)     # join the vectors back together
        str_out = str_out[1:end-1]
        str_out = replace(str_out,"X"=>";")

        str_out1 =  str[1:l_start[1]-2]*str_out*str[l_end[1]+2:end]

        str2 = split(str_out1,"\n")
        len  = maximum(textwidth.(str2))
        for i = eachindex(str2)
            str2[i] = cpad(str2[i],len," ")*"\n"
        end 
        str_out = join(str2)  
        
    else
        str_out = str
    end

    return str_out
end


function pretty_print_rheology_string(str, start="")
    # This takes the rheology string & creates a string the prints it in a nicer looking manner
    # The way we do this is decompose the string into an array of strings
    #  By placing the elements in the right order

    # WORK IN PROGRESS! 

    # Put the string elements into a matrix
    str_row = split(str,";")
    n_rows  = length(str_row)
    n_cols  = 1
    for row in str_row
        str_col = split(row,"o")
        if length(str_col)>n_cols
            n_cols = length(str_col)        
        end
    end

    A = Matrix{String}(undef, n_rows, n_cols)
    start = 1;
    for i in eachindex(str_row)    
        str_col = split(str_row[i],"o")
        for j in eachindex(str_col)    
            col = str_col[j]
            @show col
            A[i,j+start-1] = str_col[j]
            if (col[1] == '{') && (i>1)
                start = start+1
            end
            if col[end] == '}'
                start -= 1
            end

            @show start

        end
    end

    return A
   
end

# Center strings
cpad(s, n::Integer, p=" ")  = rpad(lpad(s,div(n+textwidth(s),2),p),n,p)




struct InverseCreepLaw{N} <: AbstractConstitutiveLaw{Float64}
    v::NTuple{N,AbstractConstitutiveLaw}

    function InverseCreepLaw(v::Vararg{AbstractConstitutiveLaw,N}) where {N}
        return new{N}(ntuple(i -> v[i], Val(N)))
    end

    function InverseCreepLaw(v::NTuple{N,AbstractConstitutiveLaw}) where {N}
        return new{N}(v)
    end
end

"""
    struct KelvinVoigt{N, V1, V2} <: AbstractConstitutiveLaw{Float64}
        v_el::V1
        v_vis::V2
    end

    Elastic spring and viscous dashpot in parallel. 

    τ = 2Gε + 2η̇ε
"""
struct KelvinVoigt{N1,N2,V1,V2} <: AbstractConstitutiveLaw{Float64}
    spring::V1
    dashpot::V2

    function KelvinVoigt(v1::AbstractConstitutiveLaw, v2::AbstractConstitutiveLaw)
        T1 = typeof(v1)
        T2 = typeof(v2)
        return new{1,1,T1,T2}(v1, v2)
    end

    function KelvinVoigt(
        v1::NTuple{N,AbstractConstitutiveLaw}, v2::AbstractConstitutiveLaw
    ) where {N}
        T1 = typeof(v1)
        T2 = typeof(v2)
        return new{N,1,T1,T2}(v1, v2)
    end

    function KelvinVoigt(
        v1::AbstractConstitutiveLaw, v2::NTuple{N,AbstractConstitutiveLaw}
    ) where {N}
        T1 = typeof(v1)
        T2 = typeof(v2)
        return new{1,N,T1,T2}(v1, v2)
    end

    function KelvinVoigt(
        v1::NTuple{N1,AbstractConstitutiveLaw}, v2::NTuple{N2,AbstractConstitutiveLaw}
    ) where {N1,N2}
        T1 = typeof(v1)
        T2 = typeof(v2)
        return new{N1,N2,T1,T2}(v1, v2)
    end
end

# COMPUTE STRAIN RATES

"""
    τ = 2Gε + η̇ε =  2G ̇ε Δt + η̇ε
    ̇ε = τ/(2GΔt + 2η)
"""
function compute_εII(v::KelvinVoigt{1,1,V1,V2}, τII, args) where {V1,V2}
    η = computeViscosity_τII(v.dashpot, τII, args)
    G = v.spring.G
    εII = τII / (2 * G * args.dt + 2η)

    return εII
end

@generated function compute_εII(v::KelvinVoigt{N1,N2,V1,V2}, τII, args) where {N1,N2,V1,V2}
    quote
        η, G = 0.0, 0.0

        # sum dashpot viscosities
        if N2 > 1
            η += Base.Cartesian.@nexprs $N2 i ->
                computeViscosity_τII(v.dashpot[i], τII, args)
        else
            η = computeViscosity_τII(v.dashpot, τII, args)
        end

        # sum spring shear modulus
        if N1 > 1
            G += Base.Cartesian.@nexprs $N2 i -> v.spring[i].G
        else
            G += v.spring.G
        end
        εII = 0.5 * τII / (G * args.dt + η)

        return εII
    end
end

"""
    compute_εII(v::NTuple{N, AbstractConstitutiveLaw},τII,  args)
    
Compute deviatoric strain rate given deviatoric stress invariant

"""
function compute_εII(
    v::NTuple{N,AbstractConstitutiveLaw}, τII, args; tol=1e-6, verbose=false, n=1
) where {N}
    εII = local_iterations_τII(v, τII, args; tol=tol, verbose=verbose, n=n)
    return εII
end

@inline function compute_εII!(
    εII::AbstractArray{T,nDim},
    v::NTuple{N,AbstractConstitutiveLaw},
    τII::AbstractArray{T,nDim},
    args;
) where {T,nDim,N}
    for I in eachindex(εII)
        εII[I] = compute_εII(v, τII[I], (; zip(keys(args), getindex.(values(args), I))...))
    end
end

# COMPUTE DEVIATORIC STRESS

"""
    τ = 2Gε + η̇ε =  2G ̇ε Δt + η̇ε
"""
function compute_τII(v::KelvinVoigt, εII, args)
    η = computeViscosity_εII(v.dashpot, εII, args)
    G = v.spring.G
    τII = εII * (2 * G * args.dt + η)

    return τII
end

"""
    compute_τII(v::NTuple{N, AbstractConstitutiveLaw}, εII, args; tol=1e-6, verbose=false)
    
Compute deviatoric stress invariant given strain rate 2nd invariant

"""
function compute_τII(
    v::NTuple{N,AbstractConstitutiveLaw}, εII, args; tol=1e-6, verbose=false
) where {N}
    τII = local_iterations_εII(v, εII, args; tol=tol, verbose=verbose)
    return τII
end

function compute_τII(v::CompositeRheology, εII, args; tol=1e-6, verbose=false)
    return compute_τII(v.rheology_chain, εII, args; tol=1e-6, verbose=false)
end


@inline function compute_τII!(
    τII::AbstractArray{T,nDim},
    v::NTuple{N,AbstractConstitutiveLaw},
    εII::AbstractArray{T,nDim},
    args;
) where {T,nDim,N}
    for I in eachindex(τII)
        τII[I] = compute_τII(v, εII[I], (; zip(keys(args), getindex.(values(args), I))...))
    end
end



# COMPUTE VISCOSITY

@generated function computeViscosity(
    fn::F, v::NTuple{N,AbstractConstitutiveLaw}, CII::T, args::NamedTuple; n=-1
) where {F,T,N}
    quote
        Base.@_inline_meta
        η = zero(T)
        Base.Cartesian.@nexprs $N i ->
            η += if v[i] isa Tuple
                computeViscosity(fn, v[i], CII, args; n=1) # viscosities in parallel → ηeff = 1/(η1 + η2)
            else
                fn(v[i], CII, args)^n # viscosities in series → ηeff = (1/η1 + 1/η2)^-1
            end
        return 1 / η
    end
end

## BASED ON STRAIN RATE

"""
    computeViscosity_εII(v::AbstractConstitutiveLaw, εII, args)
    
Compute viscosity given strain rate 2nd invariant for a given rheological element
"""
@inline function computeViscosity_εII(
    v::AbstractConstitutiveLaw, εII, args; cutoff=(1e16, 1e25)
)
    lower_cutoff, upper_cutoff = cutoff
    τII = compute_τII(v, εII, args)
    η = 0.5 * τII / εII
    η = max(min(upper_cutoff, η), lower_cutoff)
    return η
end

# special case for linar viscous rheology
computeViscosity_εII(v::LinearViscous, args...) = v.η.val

"""
    computeViscosity_εII(v::NTuple{N, AbstractConstitutiveLaw}, εII, args)
    
Compute viscosity given strain rate 2nd invariant
"""
function computeViscosity_εII(
    v::NTuple{N,AbstractConstitutiveLaw},
    εII,
    args;
    tol=1e-6,
    verbose=false,
    cutoff=(1e16, 1e25),
) where {N}
    lower_cutoff, upper_cutoff = cutoff
    τII = local_iterations_εII(v, εII, args; tol=tol, verbose=verbose)
    η = 0.5 * τII / εII
    η = max(min(upper_cutoff, η), lower_cutoff)
    return η
end

"""
    computeViscosity_εII(v::NTuple{N, AbstractConstitutiveLaw}, εII, args)
    
Compute viscosity given strain rate 2nd invariant
"""
function computeViscosity_εII!(
    η::AbstractArray{T,nDim},
    v::NTuple{N,AbstractConstitutiveLaw},
    εII::AbstractArray{T,nDim},
    args;
    tol=1e-6,
    verbose=false,
    cutoff=(1e16, 1e25),
) where {N,nDim,T}
    lower_cutoff, upper_cutoff = cutoff
    for I in eachindex(εII)
        argsi = (; zip(keys(args), getindex.(values(args), I))...)
        τII = local_iterations_εII(v, εII[I], argsi; tol=tol, verbose=verbose)
        ηi = 0.5 * τII / εII[I]
        η[I] = max(min(lower_cutoff, ηi), upper_cutoff)
    end
    return nothing
end

# support for multiple phases
for fn in (:computeViscosity_εII, :computeViscosity_τII)
    fn! = Symbol(fn, :!)
    @eval begin
        # local version
        function $(fn)(
            MatParam::NTuple{N,AbstractMaterialParamsStruct},
            εII::Real,
            Phase::Integer,
            args::NamedTuple;
            cutoff=(1e16, 1e25),
        ) where {N}
            lower_cutoff, upper_cutoff = cutoff
            η = compute_viscosity_param($(fn), MatParam, Phase, εII, args)
            return max(min(η, upper_cutoff), lower_cutoff)
        end

        # in-place version for Arrays
        function $(fn!)(
            η::AbstractArray,
            MatParam::NTuple{N,AbstractMaterialParamsStruct},
            εII::AbstractArray,
            Phases::AbstractArray,
            args::NamedTuple;
            cutoff=(1e16, 1e25),
        ) where {N}
            @inbounds for I in eachindex(Phases)
                k, v = keys(args), getindex.(values(args), I)
                argsi = (; zip(k, v)...)
                η[I] = $(fn)(MatParam, εII[I], Phases[I], argsi; cutoff=cutoff)
            end
            return nothing
        end
    end
end

@inline @generated function compute_viscosity_param(
    fn::F,
    MatParam::NTuple{N,AbstractMaterialParamsStruct},
    Phase::Integer,
    CII::T,
    args::NamedTuple,
) where {F,N,T}
    quote
        out = zero(T)
        Base.Cartesian.@nexprs $N i ->
            out += if MatParam[i].Phase == Phase
                computeViscosity(fn, MatParam[i].CreepLaws, CII, args)
            else
                zero(T)
            end
    end
end

@inline function computeViscosity_τII!(
    η::AbstractArray{T,nDim},
    v::NTuple{N,AbstractConstitutiveLaw},
    τII::AbstractArray{T,nDim},
    args;
    cutoff=(1e16, 1e25),
) where {T,nDim,N}
    Threads.@threads for I in eachindex(τII)
        η[I] = max(
            cutoff[1],
            min(
                cutoff[2],
                computeViscosity(
                    computeViscosity_τII,
                    v,
                    τII[I],
                    (; zip(keys(args), getindex.(values(args), I))...),
                ),
            ),
        )
    end
end

## BASED ON DEVIATORIC STRESS

"""
    computeViscosity_τII(v::AbstractConstitutiveLaw, τII, args; tol=1e-6, verbose=false)
    
Compute viscosity given stress 2nd invariant for a given rheological element
"""
@inline function computeViscosity_τII(
    v::AbstractConstitutiveLaw, τII, args; cutoff=(1e16, 1e25)
)
    lower_cutoff, upper_cutoff = cutoff
    εII = compute_εII(v, τII, args)
    η = 0.5 * τII / εII
    η = max(min(η, upper_cutoff), lower_cutoff)
    return η
end

# special case for linar viscous rheology
computeViscosity_τII(v::LinearViscous, args...) = v.η.val

"""
    compute viscosity given strain rate 2nd invariant

    τ = 2ηε -> η = τ/2/ε
"""
@inline function computeViscosity_τII(v, τII, args; tol=1e-6, verbose=false, n=1)
    εII = compute_εII(v, τII, args; tol=tol, verbose=verbose)
    η = 0.5 * τII / εII
    return η
end

"""
    computeViscosity_τII(v::NTuple{N, AbstractConstitutiveLaw}, τII, args; tol=1e-6, verbose=false)
    
Compute viscosity given deviatoric stress 2nd invariant
"""
function computeViscosity_τII(
    v::NTuple{N,AbstractConstitutiveLaw}, τII, args; tol=1e-6, verbose=false
) where {N}
    εII = local_iterations_τII(v, τII, args; tol=tol, verbose=verbose)
    η = 0.5 * τII / εII
    return η
end

## LOCAL ITERATIONS TO COMPUTE VISCOSITY

"""
Performs local iterations versus stress for a given strain rate 
"""
@inline function local_iterations_εII(
    v::NTuple{N,AbstractConstitutiveLaw}, εII, args; tol=1e-6, verbose=false
) where {N}
    # Initial guess
    η_ve = computeViscosity(computeViscosity_εII, v, εII, args) # viscosity guess
    τII = 2 * η_ve * εII # deviatoric stress guess

    # Local Iterations
    iter = 0
    ϵ = 2 * tol
    τII_prev = τII
    while ϵ > tol
        iter += 1
        f = εII - strain_rate_circuit(v, τII, args)
        dfdτII = -dεII_dτII(v, τII, args)
        τII -= f / dfdτII

        ϵ = abs(τII - τII_prev) / τII
        τII_prev = τII
        if verbose
            println(" iter $(iter) $ϵ")
        end
    end
    if verbose
        println("---")
    end
    return τII
end

"""
Performs local iterations versus strain rate for a given stress
"""
@inline function local_iterations_τII(
    v::NTuple{N,AbstractConstitutiveLaw}, τII, args; tol=1e-6, verbose=false, n=1
) where {N}
    # Initial guess
    η_ve = computeViscosity(computeViscosity_τII, v, τII, args) # viscosity guess

    εII = τII / (2 * η_ve)  # deviatoric strain rate guess
    @show η_ve, εII

    # Local Iterations
    iter = 0
    ϵ = 2 * tol
    εII_prev = εII
    while ϵ > tol
        iter += 1
        f = τII - stress_circuit(v, εII, args; n=n)
        dfdεII = -dτII_dεII(v, εII, args)
        εII -= f / dfdεII

        ϵ = abs(εII - εII_prev) / εII
        εII_prev = εII
        if verbose
            println(" iter $(iter) $ϵ")
        end
    end
    if verbose
        println("---")
    end
    return εII
end

# RHEOLOGY CIRCUITS

## STRAIN RATE 

@inline @generated function strain_rate_circuit(
    v::NTuple{N,AbstractConstitutiveLaw}, TauII, args
) where {N}
    quote
        c = 0.0
        Base.Cartesian.@nexprs $N i -> c += if v[i] isa InverseCreepLaw
            strain_rate_circuit(v[i], TauII, args)
        else
            compute_εII(v[i], TauII, args)
        end
        return c
    end
end

@generated function strain_rate_circuit(
    v_ice::InverseCreepLaw{N}, τII::_T, args
) where {_T,N}
    quote
        Base.@_inline_meta
        c = zero(_T)
        Base.Cartesian.@nexprs $N i -> c += 1 / compute_εII(v_ice.v[i], τII, args)
        return 1 / c
    end
end

## DEVIATORIC STRESS 

@inline @generated function stress_circuit(
    v::NTuple{N,AbstractConstitutiveLaw}, EpsII, args; n=1
) where {N}
    quote
        c = 0.0
        Base.Cartesian.@nexprs $N i ->
            c += if v[i] isa Tuple
                1 / stress_circuit(v[i], TauII, args; n=-1)
            else
                compute_τII(v[i], EpsII, args)^n
            end
        return c
    end
end

## VISCOSITY
@inline function viscosityCircuit_τII(v, τII, args)
    εII = strain_rate_circuit(v, τII, args)
    return 0.5 * τII / εII
end

@inline function viscosityCircuit_εII(v, εII, args)
    τII = stress_circuit(v, εII, args)
    return 0.5 * τII / εII
end

# STRESS AND STRAIN RATE DERIVATIVES

@generated function dεII_dτII(
    v::NTuple{N,AbstractConstitutiveLaw}, τII::_T, args
) where {_T,N}
    quote
        Base.@_inline_meta
        val = zero(_T)
        Base.Cartesian.@nexprs $N i -> val += dεII_dτII(v[i], τII, args)
        return val
    end
end

@generated function dτII_dεII(
    v::NTuple{N,AbstractConstitutiveLaw}, εII::_T, args
) where {_T,N}
    quote
        Base.@_inline_meta
        val = zero(_T)
        Base.Cartesian.@nexprs $N i -> val += dτII_dεII(v[i], εII, args)
        return val
    end
end
