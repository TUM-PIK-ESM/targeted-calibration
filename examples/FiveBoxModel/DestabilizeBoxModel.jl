import Pkg
Pkg.activate("examples")

using TargetedCalibrationInComplexSystemModels
using DelimitedFiles
using Statistics
using SparseArrays
using LinearAlgebra

include("CalibrateParameterNoise.jl")
include("FiveBoxModel.jl")


dt = 0.1                        # time increment for integration
noise_amp = 1e-4                # noise amplitude used to push the system out of equilibrium
steps = 500                     # Determines how often the destabilization procedure is applied


# This function chooses Δp such that it optimizes v*Δp - η_1 * (Δp',-1,0,...,0)Σ_nuv (Δp',-1,0,...,0)'
# - η_2 * (0,...,0,Δp')Σ_nuv (0,...,0,Δp')' - η_n * (Δp'*n - u)^2 - η_σ * (Δp ./ σ)' * (Δp ./ σ)
# where the first penalty term corrects the error in how the hyperplane is reached, the second penalty term controlls 
# the error of v*Δp, the third term controlls the distance from the hyperplane (so η_n should be chosen large)
# and the last term controlls the length of Δp since we want to fix the value of (Δp ./ σ)' * (Δp ./ σ) to l^2
# η_σ is increased stepwise until (Δp ./ σ)' * (Δp ./ σ) is smaller than l^2 (cf. Appendix C)
function param_update(v::Vector, σ::Vector, n::Vector, u::Number, l::Number, Σ_nuv::Matrix, η_1::Number, η_2::Number, η_n::Number, η_σ_min::Number)
    η_σ = η_σ_min
    d_p = length(v)
    Δp = zeros(d_p)
    for i in 1:100
        RHS = -v - 2*η_1 * Σ_nuv[1:d_p,d_p+1] - 2 * η_n * u * n# - 2 * l * η_σ .* σ
        matrix = -2 * η_1 * Σ_nuv[1:d_p,1:d_p] - 2 * η_2 * Σ_nuv[d_p+2:end,d_p+2:end] - 2 * η_n * (n * n') - 2 * η_σ * diagm(σ)^(-2)
        Δp[:] = matrix\RHS
        if (Δp ./ σ)' * (Δp ./ σ) < l^2
            println("i = ",i)
            break
        else
            η_σ *= 2
        end
        if i == 100
            Δp[:] = zeros(d_p)
        end
    end
    println("Δp = ", Δp)
    println("length Δp = ", (Δp ./ σ)' * (Δp ./ σ))
    println("v*Δp = ", v' * Δp)
    println("v error penalty = ", √(abs([Δp' -1 zeros(d_p)'] * Σ_nuv * [Δp; -1; zeros(d_p)])))
    println("Hyperplane penalty = ", √(abs([zeros(d_p + 1)' Δp'] * Σ_nuv * [zeros(d_p + 1); Δp])))
    println("find hyperplane penalty = ", √(abs((Δp' * n - u)^2)))
    println("length penalty = ", √(abs(sum((Δp - l*σ).^2))))
    return Δp
end

# Set initial salinity and parameter values
x_init = initSalinitiesFAMOUS
p_init = FAMOUSparams

# Decide which parameters to vary (-5 so we obtain the indicies of the parameter array, i.e. without salinities)
param_idx = [get_index["γ"],get_index["η"],get_index["K_N"],get_index["K_S"],get_index["K_IP"]] .- 5
d_p = length(param_idx)


# Progressor function (the variable x in this case represents the salinities)
function f_p(x,p,r)
    p_arr = p_init
    p_arr[param_idx] = p
    arr = [x;p_arr]
    V_N = get_V_N(arr)
    V_S = get_V_S(arr)
    V_T = get_V_T(arr)
    V_IP = get_V_IP(arr)
    V_B = get_V_B(arr)
    # Compute the additional additive noise (cf. section III.B of the paper)
    uN = √(dt)*noise_amp*r[1]/V_N
    uT = √(dt)*noise_amp*r[2]/V_T
    uS = √(dt)*noise_amp*r[3]/V_S
    uIP = √(dt)*noise_amp*r[4]/V_IP
    uB = -√(dt)*noise_amp*(r[1]+r[2]+r[3]+r[4])/V_B

    S_N_new = x[1] + dt*(rhs_S_N(x, p_arr)) + uN
    S_T_new = x[2] + dt*(rhs_S_T(x, p_arr)) + uT
    S_S_new = x[3] + dt*(rhs_S_S(x, p_arr)) + uS
    S_IP_new = x[4] + dt*(rhs_S_IP(x, p_arr)) + uIP
    S_B_new = x[5] + dt*(rhs_S_B(x, p_arr)) + uB

    return [S_N_new, S_T_new, S_S_new, S_IP_new, S_B_new]
end


# Observable function mapping to the respective AMOC strength
function f_o(x,p)
    p_arr = p_init
    p_arr[param_idx] = p
    return [q(x,p_arr)]
end


# Object of the package TargetedCalibrationInComplexSystemModels which contains the information about the system
DS = DynamicalSystem(f_p, f_o, x_init, p_init[param_idx], random_vec_length=4)

# Run system for transient time to let the system reach its equilibrium
x_tr,_,_,_ = integrateTraj(DS,100000)
DS.x_init = x_tr[:,end]

# Determine the time scale of the system
tseries = timeScale(DS,T_max=1000000, T_init = 999999)
tScale = tseries[1]
tScale_err = tseries[2]
ts = tScale - tScale_err

# compute q_mean (mean AMOC strength)
_,_,q_tr,_ = integrateTraj(DS,100000)
q_mean = mean(q_tr, dims = 2)[1]

# array to store the parameter values at each interation step
p_arr = zeros(length(param_idx), steps + 1)
p_arr[:, 1] = DS.p_init

# array to store the mean of q at each iteration step
q_mean_arr = zeros(steps + 1)
q_mean_arr[1] = q_mean

# Container for storing the noise amplitudes used for the parameter trajectory
param_noise = zeros(d_p)

# array for saving the times scales and their errors for eache iteration step
ts_arr = zeros(steps + 1)
ts_err_arr = zeros(steps + 1)
ts_arr[1] = tScale
ts_err_arr[1] = tScale_err

for i in 1:steps

    println("Iteration: $i")

    # Run system for transient time so it can equilibrate after the parameter change
    global x_tr,p_tr,_,_ = integrateTraj(DS,100000)
    DS.x_init = x_tr[:,end]

    # Fit the model for fixed parameter values (reference model) (cf. step 1 of the procudure outlined
    # in sec II of the paper)
    global x_tr,_,q_tr,noise = integrateTraj(DS,1000000)
    q_mean_arr[i+1] = mean(q_tr, dims = 2)[1]
    println("q_mean = ", q_mean_arr[i+1])
    ref_model = fitVARmodel(q_tr)

    # use the reference model to determine the time scale of the system (cf. step 2 of the procudure outlined
    # in sec II of the paper)
    global ts = timeScale(ref_model, include_error = true)
    ts_arr[i+1] = ts[1]
    println("lambda = ", 1/ts[1])
    ts_err_arr[i+1] = ts[2]

    # If the error is sufficiently small then we subtract it twice from the determined time scale which
    # ensures that the time scale ts is with hight probability larger than the real time scale of the system
    if ts[1] - 1 > 4*ts[2]
        ts = ts[1] - 2 * ts[2]
    else 
        ts = ts[1]
    end

    # Find suitable noise amplitudes for the parameters using the methods form CalibrateParameterNoise.jl
    if (i-1) % 10 == 0 
        σ_x_arr = zeros(40,d_p)
        noise_arr = zeros(40,d_p)
        for i in 1:d_p
            na, sa = pre_calibration(ts, DS, 10., 100000, min(1e-5,DS.p_init[i]/4), 1e-1, i)
            println("na = ", na)
            println("sa = ", sa)
            σ_x_arr[:,i] = [sa ; zeros(40 - length(sa))]
            noise_arr[:,i] = [na ; zeros(40 - length(na))]
        end
        param_noise[:] = calibrate_noise(σ_x_arr, noise_arr, DS.p_init/10) .* √(1-(1/ts)^2)
        println(param_noise ./ √(1-(1/ts)^2))
    end

    # Simulate the parameter series (cf. step 3 of the procudure outlined in sec II of the paper)
    p_tr = parameterSeriesGenerator(ts,DS.p_init,param_noise,n=1000000, keep_const = 1)

    # Integrate the system with the varying parameters (cf. step 4 of the procudure outlined in 
    # sec II of the paper)
    _, _, q_arr, _, q_ref_arr = integrateTraj(DS, p_tr; include_reference_traj=true)

    # Fit the parameter dependent VAR model
    B, Σ_ε, Σ_β = fitVARmodel(q_arr,q_ref_arr,p_tr,DS.p_init,ref_model)

    μ = B[1,1]
    A_q = B[1,2]
    A_p = B[1,3:d_p+2]
    A_qp = B[1,d_p+3:end]

    # Use the information obtained above to update the parameter values as described in Appendix C  
    # (cf. step 5 of the procudure outlined in sec II of the paper)

    # normal vector of the hyperplane with u = n*Δp
    n = A_p + q_mean .* A_qp
    u = q_mean * (1 - A_q) - μ - n' * DS.p_init
    # Gradient of the jacobian with respect to the parameters
    v = vec(A_qp)

    # Compute Σ_nuv which is the covariance matrix of (n', u, v')' (cf. Appendix C of the paper)
    mat_n = sparse([zeros(2,d_p); [1, q_mean] ⊗ Id(d_p)])
    mat_u = [-1; -q_mean; -DS.p_init; -q_mean * DS.p_init]
    mat_v = sparse([zeros(2+d_p,d_p); Id(d_p)])
    mat_nuv = sparse([mat_n mat_u mat_v])
    Σ_nuv = mat_nuv' * Σ_β * mat_nuv

    # since we do not want to parameters to drop below 0 and we choose l = 100 in 'param_update' we ensure by this that 
    # no parameter drops below 0 and on the same time the range in which the parameters can be adjusted depends
    # in the noise amplitudes i.e. on the size of the region explored in parameter space.
    σ = min.(param_noise, 1e-2 * DS.p_init)
    
    # update the parameters
    DS.p_init += param_update(v, σ, n, u, 100, Σ_nuv, 1e3, 1, 1e12, 1e-12)
    p_arr[:, i+1] = DS.p_init 
end


# Conpute λ values and their errors (λ moving towards one implies destabilization of the system)
λ_arr = 1 ./ ts_arr
λ_err_arr = ts_err_arr ./ ts_arr.^2


# Save data 
isdir("examples/FiveBoxModel/Data") || mkdir("examples/FiveBoxModel/Data")

# save the results in txt files
open("examples/FiveBoxModel/Data/p_arr.txt", "w") do io
    writedlm(io, p_arr)
end
open("examples/FiveBoxModel/Data/q_mean_arr.txt", "w") do io
    writedlm(io, q_mean_arr)
end
open("examples/FiveBoxModel/Data/lambda_arr.txt", "w") do io
    writedlm(io, λ_arr)
end
open("examples/FiveBoxModel/Data/lambda_err_arr.txt", "w") do io
    writedlm(io, λ_err_arr)
end