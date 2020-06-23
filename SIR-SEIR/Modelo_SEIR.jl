module Modelo_SEIR
    
    export SEIR, SEIR_norm,random_SEIR, random_SEIR_2
    """
        SEIR
        
        Modelo SEIR
    """
    function SEIR(t;σ=1/7,β=1/5,γ=0.1,S₀=1e6,E₀=12,I₀=12,R₀=5)
        S_t = Array{Float64,1}(undef, t+1)
        E_t = Array{Float64,1}(undef, t+1)
        I_t = Array{Float64,1}(undef, t+1)
        R_t = Array{Float64,1}(undef, t+1)
        S_t[1] = S₀
        E_t[1] = E₀      
        I_t[1] = I₀
        R_t[1] = R₀
        N = S₀+E₀+I₀+R₀
        k = β/N
        for i in 1:t
            S_t[i+1] = S_t[i]*(1-k*I_t[i])
            E_t[i+1] = E_t[i]*(1-σ)+k*I_t[i]*S_t[i]
            I_t[i+1] = I_t[i]*(1-γ)+σ*E_t[i]
            R_t[i+1] = R_t[i]+γ*I_t[i]
        end
        return S_t,E_t,I_t,R_t
    end

    """
        SEIR_norm
        
        Modelo SEIR normalizado
    """
    function SEIR_norm(t;σ=1/7,β=1/5,γ=0.1,S₀=0.8,E₀=0.1,I₀=0.08,R₀=0.02)
        S_t = Array{Float64,1}(undef, t+1)
        E_t = Array{Float64,1}(undef, t+1)
        I_t = Array{Float64,1}(undef, t+1)
        R_t = Array{Float64,1}(undef, t+1)
        S_t[1] = S₀
        E_t[1] = E₀
        I_t[1] = I₀
        R_t[1] = R₀
        for i in 1:t
            S_t[i+1] = S_t[i]*(1-β*I_t[i])
            E_t[i+1] = E_t[i]*(1-σ)+β*I_t[i]*S_t[i]
            I_t[i+1] = I_t[i]*(1-γ)+σ*E_t[i]
            R_t[i+1] = R_t[i]+γ*I_t[i]
        end
        return S_t,E_t,I_t,R_t
    end

    """
        random_SEIR
        
        Se corren varios modelos SEIR con parámetros aleatorios que se encuentran entre
        un intervalo, se utiliza SEIR_norm
    """
    function random_SEIR(dias,repeticiones;I₀=30,S₀=1e7,R₀=5,E₀=90,Rmin=2.2,Rmax=3.49,
                        amax=1/5,amin=1/5,gmax=1/7,gmin=1/7)
        I = []
        S = []
        R = []
        E = []
        for i in 0:dias
            push!(I,Array{Float64,1}(undef, repeticiones))
            push!(S,Array{Float64,1}(undef, repeticiones))
            push!(R,Array{Float64,1}(undef, repeticiones))
            push!(E,Array{Float64,1}(undef, repeticiones))
        end
        N = S₀+E₀+I₀+R₀
        S₀ /= N
        E₀ /= N
        I₀ /= N
        R₀ /= N
        for i in 1:repeticiones
            σ = rand()*(amax-amin)+amin
            γ = rand()*(gmax-gmin)+gmin
            β = γ*(rand()*(Rmax-Rmin)+Rmin)
            s_r,e_r,i_r,r_r = SEIR_norm(dias;σ=σ,β=β,γ=γ,I₀=I₀,S₀=S₀,R₀=R₀,E₀=E₀)
            for j in 1:dias+1
                I[j][i] = i_r[j]
                S[j][i] = s_r[j]
                R[j][i] = r_r[j]
                E[j][i] = e_r[j]
            end
        end
        t=0:1:dias
        return t,S,E,I,R,N
    end

    """
        random_SEIR_2
        
        Se corren varios modelos SEIR con parámetros aleatorios que se encuentran entre
        un intervalo, se utiliza SEIR, ademas se varian los 4 grupos S,E,I y R
    """
    function random_SEIR_2(dias,repeticiones;S0max=1e7,S0min=1e7,
                        E0max=60,E0min=60,I0max=30,I0min=30,R0max=5,R0min=5,
                         amax=1/5,amin=1/5,gmax=1/7,gmin=1/7,Rmin=2.2,Rmax=3.49)
        I = []
        S = []
        R = []
        E = []
        for i in 0:dias
            push!(I,Array{Float64,1}(undef, repeticiones))
            push!(S,Array{Float64,1}(undef, repeticiones))
            push!(R,Array{Float64,1}(undef, repeticiones))
            push!(E,Array{Float64,1}(undef, repeticiones))
        end
        for i in 1:repeticiones
            S₀= rand()*(S0max-S0min)+S0min
            E₀= rand()*(E0max-E0min)+E0min
            I₀= rand()*(I0max-I0min)+I0min
            R₀= rand()*(R0max-R0min)+R0min
            σ = rand()*(amax-amin)+amin
            γ = rand()*(gmax-gmin)+gmin
            β = γ*(rand()*(Rmax-Rmin)+Rmin)
            s_r,e_r,i_r,r_r = SEIR(dias;σ=σ,β=β,γ=γ,I₀=I₀,S₀=S₀,R₀=R₀,E₀=E₀)
            for j in 1:dias+1
                I[j][i] = i_r[j]
                S[j][i] = s_r[j]
                R[j][i] = r_r[j]
                E[j][i] = e_r[j]
            end
        end
        t=0:1:dias
        return t,S,E,I,R
    end

end