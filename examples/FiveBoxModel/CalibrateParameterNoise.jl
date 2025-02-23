using TargetedCalibrationInComplexSystemModels
using Kronecker


# Tests if the model significantly depends on the parameter with index `param_idx`.
function test_param_dependency(B::Matrix{Float64}, Σ_β::Matrix{Float64}, param_idx::Int64)
    # Determine how many parameters the model is using in total (this is specific for the five box model
    # where we have only one observalbe q, i.e. d_o = 1)
    numb_of_params = Int64((size(B,2) - 2)/2)
    # Find the collumns within the matrix B that display the parameter dependence of the VAR model on the
    # parameter of interest.
    param_loc = param_idx .+ [2, 2 + numb_of_params]
    return any(abs.(B[param_loc]) .> 200*sqrt.([Σ_β[i,i] for i in param_loc]))
end


# Only the parameter with index `param_idx` is varied and all other parameters are held constant.
# The noise amplitude (standard deviation: as the noise is later still multiplied by a factor of √(1-(1/ts)^2))
# is increased from `σ_init` to `σ_max` in steps of factor 2. The noise amplitudes in each step are saved in the
# array `noise_arr`. Further the variances of the observables are stored in the respactive place in the array 
# `σ_x_arr` if the parameter is causal for the observables. If this is not the case then 0 is stored in the
# respective place. If the variance of an observable exceeds `σ_x_max^2` then the loop is stopped leaving the
# remaining places in `noise_arr` and `σ_x_arr` filled with zeros.
# The idea is then to give the both these arrays to the function `calibrate_noise` for each of the parameters and
# the value of `σ_max` should be set really high (e.g. 100) here for reasons explained below.
function pre_calibration(ts::Float64, DS::DynamicalSystem, σ_max::Float64, traj_length::Int64, σ_init::Float64, σ_x_max::Float64, param_idx::Int64)
    max_iterations = Int64(round(log(σ_max / σ_init)/log(2)))+1
    param_noise = copy(σ_init)
    noise_arr = zeros(Float64, max_iterations)
    σ_x_arr = zeros(Float64, max_iterations)
    for i in 1:max_iterations
        noise_arr[i] = param_noise
        aux_param_noise = [param_noise  .* √(1-(1/ts)^2)]
        p_tr = parameterSeriesGenerator(ts,[DS.p_init[param_idx]],aux_param_noise,n=traj_length, keep_const = 1)
        _,_,x_tr,noise = integrateTraj(DS,10)
        ref_model = fitVARmodel(x_tr)
        p_traj = Array(DS.p_init ⊗ ones(1,traj_length))
        p_traj[param_idx,:] = p_tr
        _, _, x_arr, _, x_ref_arr = integrateTraj(DS, p_traj; include_reference_traj=true)
        B, Σ_ε, Σ_β = fitVARmodel(x_arr,x_ref_arr,p_tr,[DS.p_init[param_idx]],ref_model)
        x_mean = mean(x_arr, dims = 2)[1]
        if sum((x_arr[1,:] .- x_mean).*(x_arr[1,:] .- x_mean))/size(x_arr,2) > σ_x_max^2
            break
        end
        if test_param_dependency(B, Σ_β, 1)
            σ_x_arr[i] = sum((x_arr[1,:] .- x_mean).*(x_arr[1,:] .- x_mean))/size(x_arr,2)
        end
        if 2*param_noise <= σ_max
            param_noise *= 2
        end
    end
    return noise_arr, σ_x_arr
end


# The idea is that the want to choose noise amplitudes which influence the underlying dynamical in a similar
# way for each parameter. At the same time we want the noise amplitudes to be as large as possible without
# disturbing the system too much in order to obtain a good estimate for the VAR model. Also the noise amplitudes
# should not exceed the value of `σ_max` as this is an artificial bound set in order to prevent the parameter to
# take on unreasonable values. The arrays `σ_x_arr` and `noise_arr` are respective concatinations of the outputs of
# the pre_calibration function for each parameter. This function then chooses the smallest variance of the observables
# out of all of the parameters and sets this as the reference value. Then for each parameter the function looks for
# the noise amplitude which yields a variance of the observables closest to the reference value. If this value is
# larger than `2 * reference value` then the function looks for the noise amplitude which yields a variance of the
# observables closest to the reference value but smaller than the reference value. If this value is zero then the
# then the minumun noise value from the respective column of `noise_arr` is chosen as noise for this parameter.
# Otherwise the respective noise amplitude is chosen as noise for this parameter. If this noise amplitude is larger
# than `σ_max` then `σ_max` is chosen as noise for this parameter.
# The reson why we had to choose the value of `σ_max` really high in the function `pre_calibration` is that otherwise
# it could happen that the reference noise is choosen extemely small (e.g. in the case where one of the parameters
# has a values very close to the boundary so we dont want it to vary much more) and this would lead to the noise
# amplitude for the other parameters to be chosen unneccessarily small as well. 
function calibrate_noise(σ_x_arr, noise_arr, σ_max::Vector{Float64})
    numb_of_params = size(σ_x_arr,2)
    allowed_columns = findall(x -> x != 0, sum(σ_x_arr, dims = 1)[1,:])
    min_column = allowed_columns[argmin([max(σ_x_arr[:,i]...) for i in allowed_columns])]
    ref_value = max(σ_x_arr[:,min_column]...)
    param_noise = zeros(Float64, numb_of_params)
    for i in 1:numb_of_params
        potential_candidate = min(σ_x_arr[:,i][findall(x -> x >= ref_value, σ_x_arr[:,i])]...,100)
        if potential_candidate > 2 * ref_value
            potential_candidate = max(σ_x_arr[:,i][findall(x -> x <= ref_value, σ_x_arr[:,i])]...,0)
        end
        if potential_candidate == 0
            param_noise[i] = min(noise_arr[:,i][findall(x -> x != 0,noise_arr[:,i])]...)
        else
            param_noise[i] = noise_arr[argmin(abs.(σ_x_arr[:,i] .- potential_candidate)),i]
        end
        if param_noise[i] > σ_max[i]
            param_noise[i] = σ_max[i]
        end
    end
    return param_noise
end