using TargetedCalibrationInComplexSystemModels
using Kronecker
using Statistics

# Tests if the model significantly depends on the parameter with index `param_idx`.
function test_param_dependency(B::Matrix{Float64}, Σ_β::Matrix{Float64}, param_idx::Int64)
    numb_of_params = Int64((size(B,2) - 2)/2)
    param_loc = param_idx .+ [2, 2 + numb_of_params]
    return any(abs.(B[param_loc]) .> 200*sqrt.([Σ_β[i,i] for i in param_loc])), any(abs.(B[param_loc]) .> 100*sqrt.([Σ_β[i,i] for i in param_loc]))
end

# Function to calibrate the parameter noise as described in section III.A
# Increases all the noise values simultaneously until the variance of the observables exceeds `σ_x_max^2`. or
# the noise amplitude exceeds `σ_max`. The noise amplitudes are increased in steps of factor 2.
function calib_param_noise(ts::Float64, DS::DynamicalSystem, σ_max::Vector{Float64}, traj_length::Int64, σ_init::Vector{Float64}, σ_x_max::Float64)
    remaining_steps = Int64.(ceil.(log.(σ_max ./ σ_init)/log(2)))
    max_iterations = max(remaining_steps...)
    numb_of_params = length(σ_init)
    param_noise = copy(σ_init)
    prev_noise = copy(param_noise)
    reached_max = zeros(Int64, numb_of_params)
    for i in 1:max_iterations
        aux_param_noise = param_noise  .* √(1-(1/ts)^2)
        p_traj = parameterSeriesGenerator(ts,DS.p_init,aux_param_noise,n=traj_length, keep_const = 1)
        _,_,x_tr,noise = integrateTraj(DS,10)
        ref_model = fitVARmodel(x_tr)
        _, _, x_arr, _, x_ref_arr = integrateTraj(DS, p_traj; include_reference_traj=true)
        B, Σ_ε, Σ_β = fitVARmodel(x_arr,x_ref_arr,p_traj,DS.p_init,ref_model)
        x_mean = mean(x_arr)
        if sum((x_arr[1,:] .- x_mean).*(x_arr[1,:] .- x_mean))/size(x_arr,2) > σ_x_max^2
            break
        end
        for j in 1:numb_of_params
            if test_param_dependency(B, Σ_β, j)[1]
                reached_max[j] = 1
            end
            if !(test_param_dependency(B, Σ_β, j)[2])
                if reached_max[j] == 1
                    break
                end
            end
        end
        prev_noise = copy(param_noise)
        for j in 1:numb_of_params
            if max_iterations - i + 1 == remaining_steps[j]
                param_noise[j] *= 2
                remaining_steps[j] -= 1
            end
        end
    end
    prev_noise[findall(x -> x == 0, reached_max)] = σ_init[findall(x -> x == 0, reached_max)]
    return prev_noise
end

