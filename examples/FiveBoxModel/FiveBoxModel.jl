# Number of seconds in a year
year = 60*60*24*365.25


# Parameters of the FAMOUS model in the following order:
# [λ = 2.66,α,β,γ,η,μ,V_N,V_S,V_T,V_IP,V_B,K_N,K_S,K_IP,F_N,F_S,F_T,F_IP,S_0,T_0,T_S]
FAMOUSparams = [2.66,0.12,0.79,0.58,66.061,7.0,3.683,10.28,5.151,21.29,88.12,5.439,1.880,89.778,0.375,1.014,-0.723,-0.666,0.035,3.26,5.571]


# Initial salinities as used for the FAMOUS model in the order:
# [S_N,S_T,S_S,S_IP,S_B]
initSalinitiesFAMOUS = [35.088,35.670,34.441,34.689,34.577]


# get the respective parameter or salinity values from array
get_S_N(arr::Array{Float64}) = ((ndims(arr) == 1) ? arr[1] : arr[1,:])
get_S_T(arr::Array{Float64}) = ((ndims(arr) == 1) ? arr[2] : arr[2,:])
get_S_S(arr::Array{Float64}) = ((ndims(arr) == 1) ? arr[3] : arr[3,:])
get_S_IP(arr::Array{Float64}) = ((ndims(arr) == 1) ? arr[4] : arr[4,:])
get_S_B(arr::Array{Float64}) = ((ndims(arr) == 1) ? arr[5] : arr[5,:])
get_λ(arr::Array{Float64}) = ((ndims(arr) == 1) ? arr[6] : arr[6,:])
get_α(arr::Array{Float64}) = ((ndims(arr) == 1) ? arr[7] : arr[7,:])
get_β(arr::Array{Float64}) = ((ndims(arr) == 1) ? arr[8] : arr[8,:])
get_γ(arr::Array{Float64}) = ((ndims(arr) == 1) ? arr[9] : arr[9,:])
get_η(arr::Array{Float64}) = ((ndims(arr) == 1) ? arr[10] : arr[10,:])
get_μ(arr::Array{Float64}) = ((ndims(arr) == 1) ? arr[11] : arr[11,:])
get_V_N(arr::Array{Float64}) = ((ndims(arr) == 1) ? arr[12] : arr[12,:])
get_V_S(arr::Array{Float64}) = ((ndims(arr) == 1) ? arr[13] : arr[13,:])
get_V_T(arr::Array{Float64}) = ((ndims(arr) == 1) ? arr[14] : arr[14,:])
get_V_IP(arr::Array{Float64}) = ((ndims(arr) == 1) ? arr[15] : arr[15,:])
get_V_B(arr::Array{Float64}) = ((ndims(arr) == 1) ? arr[16] : arr[16,:])
get_K_N(arr::Array{Float64}) = ((ndims(arr) == 1) ? arr[17] : arr[17,:])
get_K_S(arr::Array{Float64}) = ((ndims(arr) == 1) ? arr[18] : arr[18,:])
get_K_IP(arr::Array{Float64}) = ((ndims(arr) == 1) ? arr[19] : arr[19,:])
get_F_N(arr::Array{Float64}) = ((ndims(arr) == 1) ? arr[20] : arr[20,:])
get_F_S(arr::Array{Float64}) = ((ndims(arr) == 1) ? arr[21] : arr[21,:])
get_F_T(arr::Array{Float64}) = ((ndims(arr) == 1) ? arr[22] : arr[22,:])
get_F_IP(arr::Array{Float64}) = ((ndims(arr) == 1) ? arr[23] : arr[23,:])
get_S_0(arr::Array{Float64}) = ((ndims(arr) == 1) ? arr[24] : arr[24,:])
get_T_0(arr::Array{Float64}) = ((ndims(arr) == 1) ? arr[25] : arr[25,:])
get_T_S(arr::Array{Float64}) = ((ndims(arr) == 1) ? arr[26] : arr[26,:])


# get the index of salinities or parameters
get_index = Dict("S_N" => 1,
                "S_T" => 2,
                "S_S" => 3,
                "S_IP" => 4,
                "S_B" => 5,
                "λ" => 6,
                "α" => 7,
                "β" => 8,
                "γ" => 9,
                "η" => 10,
                "μ" => 11,
                "V_N" => 12,
                "V_S" => 13,
                "V_T" => 14,
                "V_IP" => 15,
                "V_B" => 16,
                "K_N" => 17,
                "K_S" => 18,
                "K_IP" => 19,
                "F_N" => 20,
                "F_S" => 21,
                "F_T" => 22,
                "F_IP" => 23,
                "S_0" => 24,
                "T_0" => 25,
                "T_S" => 26,
                "parameters" => 6:26,
                "salinities" => 1:5)


"""
    q(salility_arr::Array{Float64},param_arr::Array{Float64}) 

Returns the value of the AMOC strength q of a given system. To this end the euqations q = λ[α(Tₛ-Tₙ)+β(Sₙ-Sₛ)] and
    Tₙ = μ⋅q+T₀ where solved for q
"""
function q(salinity_arr::Array{Float64},param_arr::Array{Float64}) 
    arr = vcat(salinity_arr, param_arr)
    λ = get_λ(arr)
    μ = get_μ(arr)
    α = get_α(arr)
    β = get_β(arr)
    T_S = get_T_S(arr)
    T_0 = get_T_0(arr)
    S_N = get_S_N(arr)
    S_S = get_S_S(arr)
    return (λ./(1 .+1e-1.*λ.*α.*μ)).*(α.*(T_S-T_0)+β.*(S_N-S_S)).*10
end 


# The following function return the RHS of the Differential Equations determining the dynamics of the 
# salinities.
# Note that the equations used here are not the same as the ones given by the paper since the equations
# from the paper do not respect the water volume conservation of each box
"""
    rhs_S_N(salinity_arr::Array{Float64}, param_arr::Array{Float64}; original::Bool = false)

# Arguments

- `arr::Array{Float64}`: Array of length 26 containing the current state of the system at which the
    RHS of the ODE is to be determined
- `original::Bool`: If false then the modified version of the RHS of the ODE is returned and if true
    then the original version as presented in the paper is used

# Returns

The RHS of the ODE describing the rate of change of the varialbe `S_N`
"""
function rhs_S_N(salinity_arr::Array{Float64}, param_arr::Array{Float64}; original::Bool = false)
    arr = vcat(salinity_arr, param_arr)
    Q = q(salinity_arr, param_arr)
    V_N = get_V_N(arr)
    S_0 = get_S_0(arr)
    K_N = get_K_N(arr)
    F_N = get_F_N(arr)
    F_T = get_F_T(arr)
    F_S = get_F_S(arr)
    F_IP = get_F_IP(arr)
    S_T = get_S_T(arr)
    S_N = get_S_N(arr)
    S_B = get_S_B(arr)
    if Q >= 0
        if original
            return (Q*(S_T-S_N)+K_N*(S_T-S_N)-F_N*S_0)/V_N*year*1e-10
        else
            return (Q*(S_T-S_N)+K_N*(S_T-S_N)+F_N*S_0-F_N*S_N)/V_N*year*1e-10
        end
    else
        if original
            return (abs(Q)*(S_B-S_N)+K_N*(S_T-S_N)-F_N*S_0)/V_N*year*1e-10
        else
            return (abs(Q)*(S_B-S_N)+K_N*(S_T-S_N)+F_N*S_0+(F_T+F_S+F_IP)*S_B)/V_N*year*1e-10
        end
    end
end


"""
    rhs_S_T(salinity_arr::Array{Float64}, param_arr::Array{Float64}; original::Bool = false)

# Arguments

- `arr::Array{Float64}`: Array of length 26 containing the current state of the system at which the
    RHS of the ODE is to be determined
- `original::Bool`: If false then the modified version of the RHS of the ODE is returned and if true
    then the original version as presented in the paper is used

# Returns

The RHS of the ODE describing the rate of change of the varialbe `S_T`
"""
function rhs_S_T(salinity_arr::Array{Float64}, param_arr::Array{Float64}; original::Bool = false)
    arr = vcat(salinity_arr, param_arr)
    Q = q(salinity_arr, param_arr)
    γ = get_γ(arr)
    V_T = get_V_T(arr)
    S_0 = get_S_0(arr)
    K_N = get_K_N(arr)
    K_S = get_K_S(arr)
    F_N = get_F_N(arr)
    F_T = get_F_T(arr)
    F_S = get_F_S(arr)
    F_IP = get_F_IP(arr)
    S_T = get_S_T(arr)
    S_N = get_S_N(arr)
    S_S = get_S_S(arr)
    S_IP = get_S_IP(arr)
    if Q >= 0
        if original
            return (Q*(γ*S_S+(1-γ)*S_IP-S_T)+K_S*(S_S-S_T)+K_N*(S_N-S_T)-F_T*S_0)/V_T*year*1e-10
        else
            return (Q*(γ*S_S+(1-γ)*S_IP-S_T)+K_S*(S_S-S_T)+K_N*(S_N-S_T)+F_T*S_0+F_N*(γ*S_S+(1-γ)*S_IP)+F_S*S_S+F_IP*S_IP)/V_T*year*1e-10
        end
    else
        if original
            return (abs(Q)*(S_N-S_T)+K_S*(S_S-S_T)+K_N*(S_N-S_T)-F_T*S_0)/V_T*year*1e-10
        else
            return (abs(Q)*(S_N-S_T)+K_S*(S_S-S_T)+K_N*(S_N-S_T)+F_T*S_0-F_T*S_T)/V_T*year*1e-10
        end
    end
end


"""
    rhs_S_S(salinity_arr::Array{Float64}, param_arr::Array{Float64}; original::Bool = false)

# Arguments

- `arr::Array{Float64}`: Array of length 26 containing the current state of the system at which the
    RHS of the ODE is to be determined
- `original::Bool`: If false then the modified version of the RHS of the ODE is returned and if true
    then the original version as presented in the paper is used

# Returns

The RHS of the ODE describing the rate of change of the varialbe `S_S`
"""
function rhs_S_S(salinity_arr::Array{Float64}, param_arr::Array{Float64}; original::Bool = false)
    arr = vcat(salinity_arr, param_arr)
    Q = q(salinity_arr, param_arr)
    γ = get_γ(arr)
    η = get_η(arr)
    V_S = get_V_S(arr)
    S_0 = get_S_0(arr)
    K_IP = get_K_IP(arr)
    K_S = get_K_S(arr)
    F_N = get_F_N(arr)
    F_T = get_F_T(arr)
    F_S = get_F_S(arr)
    F_IP = get_F_IP(arr)
    S_T = get_S_T(arr)
    S_B = get_S_B(arr)
    S_S = get_S_S(arr)
    S_IP = get_S_IP(arr)
    if Q >= 0
        if original
            return (Q*γ*(S_B-S_S)+K_IP*(S_IP-S_S)+K_S*(S_T-S_S)+η*(S_B-S_S)-F_S*S_0)/V_S*year*1e-10
        else
            return (Q*γ*(S_B-S_S)+K_IP*(S_IP-S_S)+K_S*(S_T-S_S)+η*(S_B-S_S)+F_S*S_0+γ*F_N*(S_B-S_S)-F_S*S_S)/V_S*year*1e-10
        end
    else
        if original
            return (abs(Q)*γ*(S_T-S_S)+K_IP*(S_IP-S_S)+K_S*(S_T-S_S)+η*(S_B-S_S)-F_S*S_0)/V_S*year*1e-10
        else
            return (abs(Q)*γ*(S_T-S_S)+K_IP*(S_IP-S_S)+K_S*(S_T-S_S)+η*(S_B-S_S)+F_S*S_0+γ*F_T*(S_T-S_S)-F_S*S_S)/V_S*year*1e-10
        end
    end
end


"""
    rhs_S_IP(salinity_arr::Array{Float64}, param_arr::Array{Float64}; original::Bool = false)

# Arguments

- `arr::Array{Float64}`: Array of length 26 containing the current state of the system at which the
    RHS of the ODE is to be determined
- `original::Bool`: If false then the modified version of the RHS of the ODE is returned and if true
    then the original version as presented in the paper is used

# Returns

The RHS of the ODE describing the rate of change of the varialbe `S_IP`
"""
function rhs_S_IP(salinity_arr::Array{Float64}, param_arr::Array{Float64}; original::Bool = false)
    arr = vcat(salinity_arr, param_arr)
    Q = q(salinity_arr, param_arr)
    γ = get_γ(arr)
    V_IP = get_V_IP(arr)
    S_0 = get_S_0(arr)
    K_IP = get_K_IP(arr)
    F_N = get_F_N(arr)
    F_T = get_F_T(arr)
    F_S = get_F_S(arr)
    F_IP = get_F_IP(arr)
    S_T = get_S_T(arr)
    S_B = get_S_B(arr)
    S_S = get_S_S(arr)
    S_IP = get_S_IP(arr)
    if Q >= 0
        if original
            return (Q*(1-γ)*(S_B-S_IP)+K_IP*(S_S-S_IP)-F_IP*S_0)/V_IP*year*1e-10
        else
            return (Q*(1-γ)*(S_B-S_IP)+K_IP*(S_S-S_IP)+F_IP*S_0+(1-γ)*F_N*(S_B-S_IP)-F_IP*S_IP)/V_IP*year*1e-10
        end
    else
        if original
            return (abs(Q)*(1-γ)*(S_T-S_IP)+K_IP*(S_S-S_IP)-F_IP*S_0)/V_IP*year*1e-10
        else
            return (abs(Q)*(1-γ)*(S_T-S_IP)+K_IP*(S_S-S_IP)+F_IP*S_0+(1-γ)*F_T*(S_T-S_IP)-F_IP*S_IP)/V_IP*year*1e-10
        end
    end
end


"""
    rhs_S_B(salinity_arr::Array{Float64}, param_arr::Array{Float64}; original::Bool = false)

# Arguments

- `arr::Array{Float64}`: Array of length 26 containing the current state of the system at which the
    RHS of the ODE is to be determined
- `original::Bool`: If false then the modified version of the RHS of the ODE is returned and if true
    then the original version as presented in the paper is used

# Returns

The RHS of the ODE describing the rate of change of the varialbe `S_B`
"""
function rhs_S_B(salinity_arr::Array{Float64}, param_arr::Array{Float64}; original::Bool = false)
    arr = vcat(salinity_arr, param_arr)
    Q = q(salinity_arr, param_arr)
    η = get_η(arr)
    γ = get_γ(arr)
    V_B = get_V_B(arr)
    S_N = get_S_N(arr)
    S_B = get_S_B(arr)
    S_S = get_S_S(arr)
    S_IP = get_S_IP(arr)
    F_N = get_F_N(arr)
    F_T = get_F_T(arr)
    F_S = get_F_S(arr)
    F_IP = get_F_IP(arr)
    if Q >= 0
        if original
            return (Q*(S_N-S_B)+η*(S_S-S_B))/V_B*year*1e-10
        else
            return (Q*(S_N-S_B)+η*(S_S-S_B)+F_N*(S_N-S_B))/V_B*year*1e-10
        end
    else
        if original
            return (abs(Q)*γ*S_S+(1-γ)*abs(Q)*S_IP-abs(Q)*S_B+η*(S_S-S_B))/V_B*year*1e-10
        else
            return (abs(Q)*γ*S_S+(1-γ)*abs(Q)*S_IP-abs(Q)*S_B+η*(S_S-S_B)+F_T*(γ*S_S+(1-γ)*S_IP-S_B)+F_S*(S_S-S_B)+F_IP*(S_IP-S_B))/V_B*year*1e-10
        end
    end
end
