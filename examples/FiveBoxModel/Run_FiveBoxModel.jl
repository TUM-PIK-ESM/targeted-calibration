# This code produces the data shown in figure 3 e, i.e. the time series when running the box model with fixed parameters (original versus destabilized parameters). 
## Added during revision of paper, November 2025, by Sebastian Bathiany.

using DelimitedFiles

include("FiveBoxModel.jl")

# overall runtime (t_max)
t_max = 10000.

noise_amp = 1e-4                # state noise amplitude used to push the system out of equilibrium

# Set time interval (in yrs) used for integrating the model
dt = 1.

# Create array with the times to consider
t_arr = collect(0:dt:t_max)

# Set salinities to initial values and same for the initial parameters for the hosing experiment using the 
# FAMOUS parameters
salinity_arr = initSalinitiesFAMOUS
FAMOUS_param_arr = FAMOUSparams

#destabilized_param_arr = FAMOUSparams 

# Initialize arrays to track and save the AMOC strength during the experiment
q_traj_FAMOUS = zeros(length(t_arr))
q_traj_destabilized = zeros(length(t_arr))

# Array to keep track of the salinities during the experiment of both the FAMOUS and
# the destabilized system
FAMOUS_arr = [salinity_arr;FAMOUS_param_arr]


# Initialize the parameters for the destabilized system based on the destabilization process that must have
# been run prior to this simulation
#destabilized_arr[get_index["γ"]]=1.533081238011955e-5  # init: 0.58
#destabilized_arr[get_index["η"]]=0.00027786696622338335   # init: 66.061
#destabilized_arr[get_index["K_N"]]=5.047929734218396e-5    # init: 5.439
#destabilized_arr[get_index["K_S"]]=18.407763040282468   # init: 1.88
#destabilized_arr[get_index["K_IP"]]=5.665062887958071  # init: 89.778

# instead, read from output:
destabilized_param_arr = copy(FAMOUSparams)
param_idx = [get_index["γ"],get_index["η"],get_index["K_N"],get_index["K_S"],get_index["K_IP"]] .- 5
param_evol = readdlm("examples/FiveBoxModel/Data/p_arr.txt", '\t', Float64, '\n')
destabilized_param_arr[param_idx] = param_evol[:,end]

##print(param_evol[:,end])  # result: [1.533081238011955e-5, 0.00027786696622338335, 5.047929734218396e-5, 18.407763040282468, 5.665062887958071]

destabilized_arr = [salinity_arr;destabilized_param_arr]

## need different initial condition in order to be on the "on" branch of AMOC also for the new parameters:
#destabilized_arr[get_index["S_N"]] = 50   # massive salinity in North Atlantic to kickstart AMOC



# Hosing strength in Sv
H = 0
 
## The numbers in the following four lines are parameter values AN, AS, AT and AIP from Wood et al., 2019.
FAMOUS_arr[get_index["F_N"]] += H*0.194
FAMOUS_arr[get_index["F_S"]] += -H*0.226
FAMOUS_arr[get_index["F_T"]] += H*0.597
FAMOUS_arr[get_index["F_IP"]] += -H*0.565

destabilized_arr[get_index["F_N"]] += H*0.194
destabilized_arr[get_index["F_S"]] += -H*0.226
destabilized_arr[get_index["F_T"]] += H*0.597
destabilized_arr[get_index["F_IP"]] += -H*0.565

## parameters needed for noise level
V_N = FAMOUS_arr[get_index["V_N"]]
V_S = FAMOUS_arr[get_index["V_S"]]
V_T = FAMOUS_arr[get_index["V_T"]]
V_IP = FAMOUS_arr[get_index["V_IP"]]
V_B = FAMOUS_arr[get_index["V_B"]]    


# Update the salinities while keeping the parameter values constant
function update_arr_constant_hosing(arr)
    sal_arr = arr[1:5]
    p_arr = arr[6:end]

    r = randn(4)

    uN = √(dt)*noise_amp*r[1]/V_N
    uT = √(dt)*noise_amp*r[2]/V_T
    uS = √(dt)*noise_amp*r[3]/V_S
    uIP = √(dt)*noise_amp*r[4]/V_IP
    uB = -√(dt)*noise_amp*(r[1]+r[2]+r[3]+r[4])/V_B
                
    arr[1] += dt*(rhs_S_N(sal_arr, p_arr)) + uN
    arr[2] += dt*(rhs_S_T(sal_arr, p_arr)) + uT
    arr[3] += dt*(rhs_S_S(sal_arr, p_arr)) + uS
    arr[4] += dt*(rhs_S_IP(sal_arr, p_arr)) + uIP
    arr[5] += dt*(rhs_S_B(sal_arr, p_arr)) + uB
    return arr
end


# Spinup: run the system with destabilized parameter values for a transient time so the salinities can adapt to the 
# parameter values
for i in 1:round(10000/dt)
    global FAMOUS_arr = update_arr_constant_hosing(FAMOUS_arr)
    global destabilized_arr = update_arr_constant_hosing(destabilized_arr)
end

# Run the experiment
for (idx,t) in enumerate(t_arr)
    FAMOUS_arr .= update_arr_constant_hosing(FAMOUS_arr)
    destabilized_arr .= update_arr_constant_hosing(destabilized_arr)
 
    q_traj_FAMOUS[idx] = q(FAMOUS_arr[1:5], FAMOUS_arr[6:end])
    q_traj_destabilized[idx] = q(destabilized_arr[1:5], destabilized_arr[6:end])
end


# Save data
open("examples/FiveBoxModel/Data/q_FAMOUS_stationary.txt", "w") do io
    writedlm(io, q_traj_FAMOUS)
end

open("examples/FiveBoxModel/Data/q_destabilized_stationary.txt", "w") do io
    writedlm(io, q_traj_destabilized)
end
