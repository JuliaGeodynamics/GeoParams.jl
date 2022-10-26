export time_τII_0D, time_τII_0D!, time_p_τII_0D, time_p_τII_0D!

# 0D RHEOLOGY functions

"""
    t_vec, τ_vec = time_τII_0D(v::CompositeRheology, εII::Number, args; t=(0.,100.), τ0=0., nt::Int64=100)

This performs a 0D constant strainrate experiment for a composite rheology structure `v`, and a given, constant, strainrate `εII` and rheology arguments `args`.
The initial stress `τ0`, the time range `t` and the number of timesteps `nt` can be modified 
"""
function time_τII_0D(v::Union{CompositeRheology,Tuple, Parallel}, εII::Number, args; t=(0.,100.), τ0=0., nt::Int64=100, verbose=true)
    t_vec    = range(t[1], t[2], length=nt)
    τ_vec    = zero(t_vec)
    εII_vec  = zero(t_vec) .+ εII
    τ_vec[1] = τ0;

    time_τII_0D!(τ_vec, v, εII_vec, args, t_vec, verbose=verbose);

    return t_vec, τ_vec
end

"""
    time_τII_0D!(τ_vec::Vector{T}, v::CompositeRheology, εII_vec::Vector{T}, args, t_vec::AbstractVector{T}) where {T}

Computes stress-time evolution for a 0D (homogeneous) setup with given strainrate vector (which can vary with time).
"""
function time_τII_0D!(τ_vec::Vector{T}, v::Union{CompositeRheology,Tuple, Parallel}, εII_vec::Vector{T}, args, t_vec::AbstractVector{T}; verbose=false) where {T}

    nt  = length(τ_vec)
    τII = τ_vec[1]

    for i=2:nt  
        dt      = t_vec[i]-t_vec[i-1]
        args    = merge(args, (; τII_old=τII, dt=dt))
        τII,    = compute_τII(v, εII_vec[i-1], args, verbose=verbose)
        
        τ_vec[i] = τII
    end

    return nothing
end



"""
    t_vec, P_vec, τ_vec = time_p_τII_0D(v::CompositeRheology, εII::Number, εvol::Number, args; t=(0.,100.), τ0=0., nt::Int64=100)

This performs a 0D constant strainrate experiment for a composite rheology structure `v`, and a given, constant, shear strainrate `εII` and volumetric strainrate `εvol`, as well as rheology arguments `args`.
The initial stress `τ0`, the time range `t` and the number of timesteps `nt` can be modified 
"""
function time_p_τII_0D(v::Union{CompositeRheology,Tuple, Parallel}, εII::Number, εvol::Number, args; t=(0.,100.), τ0=0., nt::Int64=100, verbose=true)
    t_vec    = range(t[1], t[2], length=nt)
    τ_vec    = zero(t_vec)
    P_vec    = zero(t_vec)
    εII_vec  = zero(t_vec) .+ εII
    εvol_vec = zero(t_vec) .+ εvol
    
    τ_vec[1] = τ0;

    time_p_τII_0D!(P_vec, τ_vec, v, εII_vec, εvol_vec, args, t_vec, verbose=verbose);

    return t_vec, P_vec, τ_vec
end


"""
    time_p_τII_0D!(P_vec::Vector{T}, τ_vec::Vector{T}, v::CompositeRheology, εII_vec::Vector{T}, εvol_vec::Vector{T}, args, t_vec::AbstractVector{T}) where {T}

Computes stress-time evolution for a 0D (homogeneous) setup with given shear and volumetric strainrate vectors (which can vary with time).
"""
function time_p_τII_0D!(P_vec::Vector{T}, τ_vec::Vector{T}, v::Union{CompositeRheology,Tuple, Parallel}, εII_vec::Vector{T}, εvol_vec::Vector{T}, args, t_vec::AbstractVector{T}; verbose=false) where {T}

    nt  = length(τ_vec)
    τII = τ_vec[1]
    P   = P_vec[1]

    for i=2:nt  
        dt      = t_vec[i]-t_vec[i-1]
        args    = merge(args, (; τII_old=τII, dt=dt, P_old=P))
        P,τII   = compute_p_τII(v, εII_vec[i-1], εvol_vec[i-1], args, verbose=verbose)
        
        τ_vec[i] = τII
        P_vec[i] = P
        
    end

    return nothing
end