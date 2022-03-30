abstract type AbstractViscosity end

"""
    custom_creep_law!(laws; name=:Viscosity)

Creates a custom creeping law struct such as
    ```Julia
    struct myCreep{T,T1,T2}
        val::T
        law1::T1
        law2::T2
        ...
        lawₙ::Tₙ
    end
    ```
where `val` is an n-dimensional array of viscosity values and `law_i` are creeping laws parameter structures
"""
function custom_creep_law!(laws; name=:Viscosity)
    nlaws = length(laws)
    types = [Symbol("T$i") for i in 1:nlaws]
    fields = [:($(Symbol(lowercase(String(laws[i]))))::$(Symbol("T$i"))) for i in 1:nlaws]

    Base.@eval(Main, 
        begin
            struct $(name){T, $(types...)} <: AbstractViscosity
                val::T
                $(fields...)
            end
        end
    )
end

Base.getindex(η::AbstractViscosity, I::Integer...) = η.val[I...]

for op in (:+, :-, :/, :*)
    @eval begin
        Base.$(op)(η::AbstractViscosity, A::AbstractArray) = Broadcast(Base.$(op), η.val, A)
        Base.$(op)(A::AbstractArray, η::AbstractViscosity) = Broadcast(Base.$(op), A, η.val)
    end
end

for op in (:extrema, :minimum, :maximum, :unique, :length, :size, :eachindex)
    @eval begin
        Base.$(op)(η::AbstractViscosity) = Base.$(op)(η.val)
    end
end

"""
    compute_viscosity(η::T, args::NTuple{nargs,Tuple}) where {T <: AbstractViscosity, nargs}

Compute the effective viscosity of a `AbstractViscosity`
"""
@inline @generated function compute_viscosity(η::T, args) where T <: AbstractViscosity
    functors = fieldnames(T)[2:end]
    nf = length(functors)
    if nf > 1
        ex = 0.0
        for i in 1:nf
            functor = functors[i]
            ex = :(1/(η.$(functor)(args[$i]...)) + $ex)
        end
        
        return :(1/$(ex))
    else

        return :(η.$(fieldnames(T)[2])(args...))
    end
end

# case for several material phases that are identified with an integer
@inline @generated function compute_viscosity(η::T, phase::Integer, args) where T <: AbstractViscosity
    functors = fieldnames(T)[2:end]
    nf = length(functors)
    if nf > 1
        ex = 0.0
        for i in 1:nf
            functor = functors[i]
            ex = :(1/(getindex(η.$(functor), phase)(args[$i]...)) + $ex)
        end
        
        return :(1/$(ex))
    else
        
        return :(getindex(η.$(fieldnames(T)[2]), phase)(args...))
    end
end