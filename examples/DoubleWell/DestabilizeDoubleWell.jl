import Pkg
Pkg.activate("examples")

using TargetedCalibrationInComplexSystemModels
using DelimitedFiles
using Statistics

include("CalibrateParameterNoise.jl")

dt = 0.001                  # time increment for integration
noise_amp = 1e-4            # noise amplitude used to push the system out of equilibrium
steps = 18                  # Determines how often the destabilization procedure is applied

# progressor function
function f_p(x,p,r)
    return [x[1] + dt*(p[1]*x[1]^3 + p[2]*x[1]^2 + p[3]*x[1] + p[4]) + √(dt)*noise_amp*r[1]]
end
# observable function
f_o(x,p) = x 

# initial parameters
p_init = [-1., 0., 1., 0.]
# initial state variable
x_init = [1.]

# Initialize the corresponding dynamical system
DS = DynamicalSystem(f_p, f_o, x_init, p_init, random_vec_length=1)

# Array to store the time scales where the time scale ts is defined as 1/Jacobian of the system with 
# Jacobian being the derivative of the VAR model w.r.t. the observable x.
ts_arr = zeros(steps + 1)
ts_err_arr = zeros(steps + 1)

# Determine the time scale of the system
ts = timeScale(DS,T_max=1000000, T_init = 999999)
ts_arr[1] = ts[1]
ts_err_arr[1] = ts[2]
# This ensures that the time scale ts is with hight probability larger than the time scale of the system
ts = ts[1] - 2*ts[2]

# compute x_mean for the purpose of analyzing how the mean changes over the course of the destabilization
# procedure. That is, this is not a part of the destabilization process but for purely analytic reasons.
# The mean values are saved in the array `x_mean_arr`
x_mean_arr = zeros(steps + 1)
_,_,x_tr,_ = integrateTraj(DS,100000)
x_mean_arr[1] = mean(x_tr)

# Array to track the evolution of the parameters
p_arr = zeros(4, steps + 1)
p_arr[:, 1] = DS.p_init

for i in 1:steps

    println("Iteration: $i")
    
    # Determine the time scale of the system (steps 1 and 2 in the destibilization process presented in section
    # II of the paper)
    global ts = timeScale(DS,T_max=1000000, T_init = 999999)
    ts_arr[i+1] = ts[1]
    ts_err_arr[i+1] = ts[2]
    ts = ts[1] - 2*ts[2]

    # compute x_mean (angain only to analyze the model and not as a part of the destabilization process)
    global _,_,x_tr,_ = integrateTraj(DS,100000)
    x_mean_arr[i+1] = mean(x_tr)

    # Calibrate the noise amplitude vor the parameter series (step 3 in the destibilization process presented in 
    # section II of the paper). The function calib_param_noise is defined in "CalibrateParameterNoise.jl"
    # and is implemented as described in section III.A
    param_noise = calib_param_noise(ts, DS, ones(4), 100, ones(4) * 1e-3, 1e-3) .* √(1-(1/ts)^2)

    # Record a trajectory with parameter values varying around their initial values of length n=100
    p_traj = parameterSeriesGenerator(ts,DS.p_init,param_noise,n=100, keep_const = 1)

    # Integrate a trajectory for 1000 time steps while recording the additive noise and fixed parameter values
    global _,_,x_tr,noise = integrateTraj(DS,1000)

    # Using the trajectory above we fit the corresponding parameter-independent VAR model which is needed for
    # the advanced method of determining the parameter dependent model as shown in step 4 in section II. (In
    # the paper the same model that was utilized to determine the time scale was used for this purpose but we 
    # use a new model here to show that a model derived from a much shorter trajectory is also sufficient for 
    # this job)
    # Assuming that we want to utilize our knowledge on how the noise influences the dynamics we use:
    y_tr = copy(x_tr)
    x_tr = x_tr - √(dt)*noise_amp .* noise
    global ref_model = auxFitVARmodel(Z(x_tr),Y(y_tr))

    # Assuming that the noise trajectory is unknown, i.e. intrinsic to the progressor function,
    # we would instead directly fit the model as follows (note that this will initially show the same
    # tendencies as the procedure above but is more unstable as the tipping point is approached, which can 
    # be fixed by in turn adjusting the hyperparameters):
    # global ref_model = fitVARmodel(x_tr)

    # Simulate the system with the varying parameter series and record the reference trajectory as explained
    # in step 4 of section II. The length of the trajectory is determined by the length of p_traj, i.e. in 
    # this case n=100
    _, _, x_arr, _, x_ref_arr = integrateTraj(DS, p_traj; include_reference_traj=true)
    
    # fit the parameter dependent VAR model
    B, Σ_ε, Σ_β = fitVARmodel(x_arr,x_ref_arr,p_traj,DS.p_init,ref_model)

    # Determine the VAR estimates as described in the paper
    p_length = length(DS.p_init)
    μ = B[1,1]
    A_x = B[1,2]
    A_p = B[1,3:p_length+2]
    A_xp = B[1,p_length+3:end]

    # In the following we update the parameter values (cf. step 5 in section II) by the particular process
    # introduces in section III.A

    # normal vector of the hyperplane
    n = A_p + A_xp

    # projection of A_xp onto the direction of (A_p + x_mean .* A_xp)
    proj = n' * A_xp / (n' * n) .* n

    # Correction the average of x_mean
    corr = - (μ + n' * DS.p_init - (1 - A_x)) / (n' * n) .* n

    add = A_xp - proj

    # A_x + A_xp * (p_init + corr + λ * add) = 1
    λ = (1 - A_x - A_xp' * (DS.p_init + corr)) / (A_xp' * add)

    # Since the noise amplitudes always take values of the order 0.05 we do not want to change any of the parameters more
    # than that. By not making the step size noise dependent but keep the largest step size at a value of 0.05 we ensure
    # that the parameter converge with a steady rate. 
    # Also since we do not want the system to bifurcate we do not allow any steps larger than λ/2.
    factor = min(λ/2, 0.05/max(abs.(add)...))

    # update the parameters
    DS.p_init += corr + factor .* add

    # track parameter changes
    p_arr[:, i+1] = DS.p_init
end

# Conpute λ values and their errors
λ_arr = 1 ./ts_arr
λ_err_arr = ts_err_arr ./ ts_arr.^2



# Save data 
isdir("examples/DoubleWell/Data") || mkdir("examples/DoubleWell/Data")

open("examples/DoubleWell/Data/p_arr.txt", "w") do io
    writedlm(io, p_arr)
end
open("examples/DoubleWell/Data/x_mean_arr.txt", "w") do io
    writedlm(io, x_mean_arr)
end
open("examples/DoubleWell/Data/lambda_arr.txt", "w") do io
    writedlm(io, λ_arr)
end
open("examples/DoubleWell/Data/lambda_err_arr.txt", "w") do io
    writedlm(io, λ_err_arr)
end