# The goal is to determine how long the recorded trajectory has to be to get an equally good estimate of the
# parametrised VAR model when multiple parameters are varied simultateously compared to the case where only one
# parameter is varied at a time. As for the computation of the gradient on how to update the system only the
# coefficients of the terms x*p are relevant we compare the errors of these coefficients. The error comparison 
# is done by computing the error of the sum of all the respecitve coefficients to the sum of errors of the
# coefficients where only one parameter is varied at a time. The error for the case where only one coefficient
# is varied at a time is computed by averaging over 100 trajectories of length 1000. The error for the case
# where multiple parameters are varied simultateously is computed by averaging over 20 trajectories of increasing
# length until the condition is met that the error of the sum of the coefficients is smaller than the sum of the
# errors of the coefficients where only one parameter is varied at a time. The length of the trajectory is 
# increased in steps of 1000 starting from 6000. The computation is done for all possible combinations of 
# parameters.

using TargetedCalibrationInComplexSystemModels
using Combinatorics
using LinearAlgebra
using Base
using Kronecker
using SparseArrays
using DelimitedFiles
using Statistics


include("../CalibrateParameterNoise.jl")
include("../FiveBoxModel.jl")



dt = 0.1                # time increment for integration
noise_amp = 1e-4        # Noise amplitude of the system


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

# Need an arbitrary ref model in order to compute the errors of interest
_,_,q_tr,noise = integrateTraj(DS,10)
ref_model = fitVARmodel(q_tr)

# The possible trajectory lengths

#n_iteration = collect(1000:1000:1000000)                   # Is used when working with optimised parameter noise (see below)
n_iteration = vcat([10^i .* [1,2,5] for i in 3:8]...)       # Is used when working with non-optimised parameter noise (see below)

# p_errors contains the errors where only one parameter is varied at a time
p_errors = zeros(5)

# n_arr contains the trajectory lengths where the error of the sum of the coefficients is smaller than the sum
# of the errors of the coefficients where only one parameter is varied at a time. The first entry corresponds to
# the case where only one parameter is varied at a time, the second entry corresponds to the case where two
# parameters are varied at a time and so on such that n_arr is an array of arrays
n_arr = [[],[],[],[],[]]
push!(n_arr[1], n_iteration[1])

# Loop over all possible combinations of parameters
for indices in collect(powerset(collect(1:5)))[2:end]

    # set up the dynamical system where only the parameters in indices are varied
    function current_f_p(x,p,r)
        complement = setdiff(collect(1:5), indices)
        params = zeros(5)
        params[complement] = p_init[param_idx][complement]
        params[indices] = p
        return f_p(x,params,r)
    end

    function current_f_o(x,p)
        complement = setdiff(collect(1:5), indices)
        params = zeros(5)
        params[complement] = p_init[param_idx][complement]
        params[indices] = p
        return f_o(x,params)
    end

    DS_curr = DynamicalSystem(current_f_p, current_f_o, x_init, p_init[param_idx][indices], random_vec_length=4)
    
    # We used two sets of parameter noise. The first one is optimised (via trial and error) such that the errors of the coefficients
    # are sufficiently smaller when variing multiple parameters simultateously compared to the case where only
    # one parameter is varied at a time. The second one is not optimised and determined using the same method
    # as introduced in the paper.

    #global param_noise = [0.16384*4, 2.62144/16, 0.32768/2, 0.16384, 2.62144/8] .* √(1-(1/ts)^2)      # oprimised parameter noise
    global param_noise = [0.02048, 5.24288, 0.32768, 0.08192, 2.62144] .* √(1-(1/ts)^2)                # non-optimised parameter noise

    # Fill the p_error array with the errors of the coefficients where only one parameter is varied at a time
    if length(indices) == 1
        l = 100
        p_average = zeros(l)
        for i in 1:l
            p_traj = parameterSeriesGenerator(ts,DS_curr.p_init,param_noise[indices],n=n_iteration[1], keep_const = 1)
            _, _, x_arr, _, x_ref_arr = integrateTraj(DS_curr, p_traj; include_reference_traj=true)
            B, Σ_ε, Σ_β = fitVARmodel(x_arr,x_ref_arr,p_traj,DS_curr.p_init,ref_model)

            p_length = length(DS_curr.p_init)
            p_average[i] = Σ_β[p_length+3,p_length+3]
        end
        p_errors[indices[1]] = mean(p_average)
        println(p_errors)

    # Else fill the n_arr array with the trajectory lengths where the error of the sum of the coefficients is
    # smaller than the sum of the errors of the coefficients where only one parameter is varied at a time
    # The length of the trajectory is increased according to n_arr starting from 1000
    # The final result is an average over 20 runs
    else
        n_average = zeros(20)
        println("indices = ", indices)
        for i in 1:length(n_average)
            println("i = ", i)
            for n in n_iteration
                println("n = ", n)
                p_traj = parameterSeriesGenerator(ts,DS_curr.p_init,param_noise[indices],n=n, keep_const = 1)
                _, _, x_arr, _, x_ref_arr = integrateTraj(DS_curr, p_traj; include_reference_traj=true)
                B, Σ_ε, Σ_β = fitVARmodel(x_arr,x_ref_arr,p_traj,DS_curr.p_init,ref_model)
                p_length = length(DS_curr.p_init) 
                # Here the condition is checked
                if sum(Σ_β[p_length+3:p_length+2+length(indices),p_length+3:p_length+2+length(indices)]) < sum(p_errors[indices])
                    n_average[i] = n
                    break
                end
                # message if the condition is not met for the largest trajectory length
                if n == n_iteration[end]
                    sum([(Σ_β[p_length+2+i,p_length+2+i]) for i in 1:length(indices)])
                    println("it didnt work")
                end
            end
        end
        n = mean(n_average)
        push!(n_arr[length(indices)], n)
        println("n = ", n)
    end
end


# Save the results
isdir("examples/FiveBoxModel/Data") || mkdir("examples/FiveBoxModel/Data")

open("examples/FiveBoxModel/Data/ComputationalCost.txt", "w") do io
    writedlm(io, n_arr)
end

# Use the following for the optimised parameter noise
#open("examples/FiveBoxModel/Data/ComputationalCost_optimised.txt", "w") do io
#    writedlm(io, n_arr)
#end