# This code produces the data that can reproduce the figure 4 a from the paper.

# The goal is to determine how long the recorded trajectory has to be to get an equally good estimate of the
# parametrised VAR model when multiple parameters are varied simultateously compared to the case where only one
# parameter is varied at a time. As for the computation of the gradient on how to update the system only the
# coefficients of the terms x*p are relevant we compare the errors of these coefficients. The error comparison 
# is done by computing the error of the sum of all the respecitve coefficients to the sum of errors of the
# coefficients where only one parameter is varied at a time. The error for the case where only one coefficient
# is varied at a time is computed by averaging over 100 trajectories of length 6000. The error for the case
# where multiple parameters are varied simultateously is computed by averaging over 20 trajectories of increasing
# length until the condition is met that the error of the sum of the coefficients is smaller than the sum of the
# errors of the coefficients where only one parameter is varied at a time. The length of the trajectory is 
# increased in steps of 1000 starting from 6000. The computation is done for all possible combinations of 
# parameters.

using TargetedCalibrationInComplexSystemModels
using Combinatorics
using LinearAlgebra
using DelimitedFiles


include("../CalibrateParameterNoise.jl")


dt = 0.001
noise_amp = 1e-4

# progressor function
function f_p(x,p,r)
    return [x[1] + dt*(p[1]*x[1]^3 + p[2]*x[1]^2 + p[3]*x[1] + p[4]) + √(dt)*noise_amp*r[1]]
end

# observable function
f_o(v,p) = v 

# initial parameters
p_init = [-1., 0., 1., 0.]
# initial state variable
x_init = [1.]

DS = DynamicalSystem(f_p, f_o, x_init, p_init, random_vec_length=1)
# Determine the time scale of the system
ts = timeScale(DS,T_max=1000000, T_init = 999999)
# This ensures that the time scale ts is definitely larger than the time scale of the system
ts = ts[1] - 2*ts[2]

# Need an arbitrary ref model in order to compute the errors of interest
_,_,x_tr,noise = integrateTraj(DS,10)
ref_model = fitVARmodel(x_tr)

# The possible trajectory lengths
n_iteration = collect(6000:1000:2000000)
n_1 = 6000

# p_errors contains the errors where only one parameter is varied at a time
p_errors = zeros(4)

# n_arr contains the minimal trajectory lengths at which the error of the sum of the coefficients is smaller than the sum
# of the errors of the coefficients where only one parameter is varied at a time. The first entry corresponds to
# the case where only one parameter is varied at a time, the second entry corresponds to the case where two
# parameters are varied at a time and so on such that n_arr is an array of arrays
n_arr = [[],[],[],[]]
push!(n_arr[1], n_1)

# Loop over all possible combinations of parameters (empty set is excluded)
for indices in collect(powerset(collect(1:4)))[2:end]

    # set up the dynamical system where only the parameters in indices are varied

    # Define the correct progressor function corresponding to the currently varied parameters
    function current_f_p(x,p,r)
        complement = setdiff(collect(1:4), indices)
        params = zeros(4)
        params[complement] = p_init[complement]
        params[indices] = p
        return f_p(x,params,r)
    end
    DS_curr = DynamicalSystem(current_f_p, f_o, x_init, p_init[indices], random_vec_length=1)

    # Use same parameter noise as in the previous analysis of the Doublewell system
    global param_noise = 0.05*ones(length(indices)) .* √(1-(1/ts)^2)

    # Fill the p_error array with the errors of the coefficients where only one parameter is varied at a time
    # resulting in a parameter-trajectory of length n_1
    # To obtain an accurate results for the error the error of interest is averaged over 100 independent runs
    if length(indices) == 1
        p_average = zeros(100)
        for i in 1:length(p_average)
            p_traj = parameterSeriesGenerator(ts,DS_curr.p_init,param_noise,n=n_1, keep_const = 1)
            _, _, x_arr, _, x_ref_arr = integrateTraj(DS_curr, p_traj; include_reference_traj=true)
            B, Σ_ε, Σ_β = fitVARmodel(x_arr,x_ref_arr,p_traj,DS_curr.p_init,ref_model)

            p_length = length(DS_curr.p_init)
            p_average[i] = Σ_β[p_length+3,p_length+3]
        p_errors[indices[1]] = mean(p_average)
        end
        println(p_errors)

    # Else fill the n_arr array with the trajectory lengths where the error of the sum of the coefficients is
    # smaller than the sum of the errors of the coefficients where only one parameter is varied at a time
    # The length of the trajectory is increased in steps of 1000 starting from 6000
    # The final result is an average over 20 runs
    else
        n_average = zeros(20)
        println("indices = ", indices)
        for i in 1:length(n_average)
            println("i = ", i)
            for n in n_iteration
                p_traj = parameterSeriesGenerator(ts,DS_curr.p_init,param_noise,n=n, keep_const = 1)
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
isdir("examples/DoubleWell/Data") || mkdir("examples/DoubleWell/Data")

open("examples/DoubleWell/Data/ComputationalCost.txt", "w") do io
    writedlm(io, n_arr)
end