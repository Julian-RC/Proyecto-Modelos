module Modelo_SIR
    
    export SIR, SIR_norm, random_SIR, random_SIR_2
    
    """
        SIR
        
        Modelo SIR
    """
    function SIR(t;β=0.1,γ=1/5,S₀=1e4,I₀=14,R₀=1)
        S_t = Array{Float64,1}(undef, t+1)
        I_t = Array{Float64,1}(undef, t+1)
        R_t = Array{Float64,1}(undef, t+1)
        N = I₀+S₀+R₀
        S_t[1] = S₀
        I_t[1] = I₀
        R_t[1] = R₀
        k = β/N
        @inbounds for i in 1:t
            S_t[i+1] = S_t[i]*(1-k*I_t[i])
            I_t[i+1] = I_t[i]*(1+k*S_t[i]-γ)
            R_t[i+1] = R_t[i]+γ*I_t[i]
        end
        return S_t,I_t,R_t
    end

    """
        SIR_norm
        
        Modelo SIR normalizado
    """
    function SIR_norm(t;β=0.1,γ=1/5,S₀=0.9,I₀=0.08,R₀=0.02)
        S_t = Array{Float64,1}(undef, t+1)
        I_t = Array{Float64,1}(undef, t+1)
        R_t = Array{Float64,1}(undef, t+1)
        S_t[1] = S₀
        I_t[1] = I₀
        R_t[1] = R₀
        @inbounds for i in 1:t
            S_t[i+1] = S_t[i]*(1-β*I_t[i])
            I_t[i+1] = I_t[i]*(1+β*S_t[i]-γ)
            R_t[i+1] = R_t[i]+γ*I_t[i]
        end
        return S_t,I_t,R_t
    end

    """
        random_SIR
        
        Se corren varios modelos SIR con parámetros aleatorios que se encuentran entre
        un intervalo, se utiliza SIR_norm
    """
    function random_SIR(dias,repeticiones;S₀=1e7,I₀=30,R₀=5,
                        Rmin=2.2,Rmax=3.49,gmax=1/7,gmin=1/7)
        I = []
        S = []
        R = []
        @inbounds for i in 0:dias 
            push!(I,Array{Float64,1}(undef, repeticiones))
            push!(S,Array{Float64,1}(undef, repeticiones))
            push!(R,Array{Float64,1}(undef, repeticiones))
        end
        N=S₀+I₀+R₀
        S₀/=N
        I₀/=N
        R₀/=N
        @inbounds for i in 1:repeticiones
            γ = rand()*(gmax-gmin)+gmin
            β = γ*(rand()*(Rmax-Rmin)+Rmin)
            s_r,i_r,r_r = SIR_norm(dias,β=β,γ=γ,S₀=S₀,I₀=I₀,R₀=R₀)
            @inbounds for j in 1:dias+1
                I[j][i] = i_r[j]
                S[j][i] = s_r[j]
                R[j][i] = r_r[j]
            end
        end
        t=0:1:dias
        return t,S,I,R,N
    end

    """
        random_SIR_2
        
        Se corren varios modelos SIR con parámetros aleatorios que se encuentran entre
        un intervalo, se utiliza SIR, ademas se varian los 3 grupos S,I y R
    """
    function random_SIR_2(dias,repeticiones;S0max=1e7,S0min=1e7,I0max=30,I0min=30
                        ,R0max=5,R0min=5,Rmin=2.2,Rmax=3.49,gmax=1/7,gmin=1/7)
        I = []
        S = []
        R = []
        @inbounds for i in 0:dias 
            push!(I,Array{Float64,1}(undef, repeticiones))
            push!(S,Array{Float64,1}(undef, repeticiones))
            push!(R,Array{Float64,1}(undef, repeticiones))
        end
        @inbounds for i in 1:repeticiones
            S₀= rand()*(S0max-S0min)+S0min
            I₀= rand()*(I0max-I0min)+I0min
            R₀= rand()*(R0max-R0min)+R0min
            γ = rand()*(gmax-gmin)+gmin
            β = γ*(rand()*(Rmax-Rmin)+Rmin)
            s_r,i_r,r_r = SIR(dias,β=β,γ=γ,S₀=S₀,I₀=I₀,R₀=R₀)
            @inbounds for j in 1:dias+1
                I[j][i] = i_r[j]
                S[j][i] = s_r[j]
                R[j][i] = r_r[j]
            end
        end
        t=0:1:dias
        return t,S,I,R
    end
end