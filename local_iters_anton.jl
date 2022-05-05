function Local_Iterations()
    # Physics
    εII_ve = 1e-14          # second invariant of effective viscoelastic strain rate, 1/s
    μ = 1e11           # shear modolus, Pa
    A = 1.67e-24       # pre-exponentialfactor, Pa^(-n) s^(-1) 
    n = 3.3            # power-law exponent, []
    E = 187e3          # activation energy, J/mol 
    R = 8.3145         # universal gas constant, J/mol/K
    T = 500 + 273.15     # temperature, Ks
    # numerics 
    dt = 1e3 * (365.25 * 24 * 3600)   # time step, s
    nt = 5                   # number of iterations, []
    # Initialiazation
    εII_vis = εII_ve
    τII = []
    err = []
    iter = 0
    # Local Iterations
    #=while (iter < 2 || err[end] > 1e-10) && iter < nt  # Picard
       iter = iter + 1
       η       = A^(-1.0/n) * εII_vis^((1.0-n)/n)*exp(E/(n*R*T))
       η_ve    = 1.0/(1.0/η + 1.0/(μ*dt))
       if iter == 1
           push!(err,abs(- 2.0*η_ve*εII_ve))
       else
           push!(err,abs((τII[end] - 2.0*η_ve*εII_ve)/(τII[end] + 2.0*η_ve*εII_ve)))
       end
       push!(τII,2.0*η_ve*εII_ve)
       εII_vis = τII[end]/(2.0*η)
    end =#
    # Local Iterations
    η = 0.5 * A^(-1.0 / n) * εII_ve^((1.0 - n) / n) * exp(E / (n * R * T)) # guess
    η_ve = 1.0 / (1.0 / η + 1.0 / (μ * dt)) # guess
    τII = 2 * η_ve * εII_ve            # guess

    while iter < 5  # Picard
        iter = iter + 1
        f = εII_ve - A * τII^n * exp(-E / R / T) - τII / (2 * μ * dt)
        dfdτII = 0.0 - n * A * τII^(n - 1.0) * exp(-E / R / T) - 1.0 / (2 * μ * dt)
        τII = τII - f / dfdτII
        # println(f)
        # println(τII)
        #=if iter == 1
            push!(err,abs(- 2.0*η_ve*εII_ve))
        else
            push!(err,abs((τII[end] - 2.0*η_ve*εII_ve)/(τII[end] + 2.0*η_ve*εII_ve)))
        end
        push!(τII,2.0*η_ve*εII_ve)
        εII_vis = τII[end]/(2.0*η)=#
    end
    #display(err[end])
    #p1 = plot(log10.(abs.(err)))
    #p2 = plot((τII))

    #display(plot(p1,p2))

    return nothing
end
@btime Local_Iterations()
