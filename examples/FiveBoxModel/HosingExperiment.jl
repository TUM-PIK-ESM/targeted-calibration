# This code produces the data shown in figures 3 c and d, i.e. the results of the hosing experiments with
# original and destabilized parameters. This code can only be run after executing the 'DestabilizeBoxModel.jl' 
# which finds a parameter configuration of a destabilized model.

using DelimitedFiles

include("FiveBoxModel.jl")

# time in yrs up to when the hosing is first positive (t_decrease) then time where hosing starts to increase
# again (t_increase) and overall runtime (t_max)
t_max = 24000.
t_decrease = 4000.
t_increase = 16000.

# Set time interval (in yrs) used for integrating the model
dt = 1.

# Create array with the times to consider
t_arr = collect(0:dt:t_max)

# Set salinities to initial values and same for the initial parameters for the hosing experiment using the 
# FAMOUS parameters
salinity_arr = initSalinitiesFAMOUS
FAMOUS_param_arr = FAMOUSparams

# Initialize the parameters for the destabilized system based on the destabilization process that must have
# been run prior to the hosing experiment
destabilized_param_arr = copy(FAMOUSparams)
param_idx = [get_index["γ"],get_index["η"],get_index["K_N"],get_index["K_S"],get_index["K_IP"]] .- 5
param_evol = readdlm("examples/FiveBoxModel/Data/p_arr.txt", '\t', Float64, '\n')
destabilized_param_arr[param_idx] = param_evol[:,end]


# Update both the salinities and the parameters (the parameters now change in time due to the hosing)
function update_arr(arr,t)
    H_1 = 5e-4

    sal_arr = arr[1:5]
    p_arr = arr[6:end]

    arr[1] += dt*(rhs_S_N(sal_arr, p_arr)) 
    arr[2] += dt*(rhs_S_T(sal_arr, p_arr)) 
    arr[3] += dt*(rhs_S_S(sal_arr, p_arr)) 
    arr[4] += dt*(rhs_S_IP(sal_arr, p_arr))
    arr[5] += dt*(rhs_S_B(sal_arr, p_arr)) 

    if (t <= t_decrease) || (t >= t_increase)
        sign = 1
    else
        sign = -1
    end
    ## The numbers in the following four lines are parameter values AN, AS, AT and AIP from Wood et al., 2019.
    arr[get_index["F_N"]] += dt*sign*H_1*0.194
    arr[get_index["F_S"]] += -dt*sign*H_1*0.226
    arr[get_index["F_T"]] += dt*sign*H_1*0.597
    arr[get_index["F_IP"]] += -dt*sign*H_1*0.565

    return arr
end


# Update the salinities while keeping the parameter values constant
function update_arr_no_hosing(arr)
    sal_arr = arr[1:5]
    p_arr = arr[6:end]

    arr[1] += dt*(rhs_S_N(sal_arr, p_arr)) 
    arr[2] += dt*(rhs_S_T(sal_arr, p_arr)) 
    arr[3] += dt*(rhs_S_S(sal_arr, p_arr)) 
    arr[4] += dt*(rhs_S_IP(sal_arr, p_arr))
    arr[5] += dt*(rhs_S_B(sal_arr, p_arr)) 
    return arr
end


# Compute the hosing strength as a function of time t
function hosing_traj(t)
    H_1 = 5e-4
    H = zeros(length(t))
    dt = t[2]-t[1]
    for (i,val) in enumerate(t[2:end])
        if (val <= t_decrease) || (val >= t_increase)
            H[i+1] = H[i] + H_1*dt
        else
            H[i+1] = H[i] - H_1*dt
        end
    end
    return H
end

# Trajectory of the hosing strength
Hosing = hosing_traj(t_arr)

# Initialize arrays to track and save the AMOC strength during the hosing experiment
q_traj_FAMOUS = zeros(length(t_arr))
q_traj_destabilized = zeros(length(t_arr))

# Array to keep track of the salinities and parameters during the hosing experiments of both the FAMOUS and
# the destabilized system
FAMOUS_arr = [salinity_arr;FAMOUS_param_arr]
destabilized_arr = [salinity_arr;destabilized_param_arr]

# run the system with destabilized parameter values for a transient time so the salinities can adapt to the 
# new parameter values
for i in 1:round(1000/dt)
    global destabilized_arr = update_arr_no_hosing(destabilized_arr)
end

# Run the hosing experiment
for (idx,t) in enumerate(t_arr)
    FAMOUS_arr .= update_arr(FAMOUS_arr,t)
    destabilized_arr .= update_arr(destabilized_arr,t)
    q_traj_FAMOUS[idx] = q(FAMOUS_arr[1:5], FAMOUS_arr[6:end])
    q_traj_destabilized[idx] = q(destabilized_arr[1:5], destabilized_arr[6:end])
end


# Save data
open("examples/FiveBoxModel/Data/Hosing.txt", "w") do io
    writedlm(io, Hosing)
end

open("examples/FiveBoxModel/Data/q_FAMOUS_during_hosing.txt", "w") do io
    writedlm(io, q_traj_FAMOUS)
end

open("examples/FiveBoxModel/Data/q_destabilized_during_hosing.txt", "w") do io
    writedlm(io, q_traj_destabilized)
end
