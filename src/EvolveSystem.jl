# In this file we define methods that allow us to integrate a given System dependent on the dynamics of the
# parameters

using Kronecker

include("DynamicalSystem.jl")


"""
    integrateTraj(DS::DynamicalSystem, p_traj::Matrix{Float64}; include_reference_traj::Bool=false)

This function iteratively computes the evolution of the state variables x and the observables o of a system forced by
    a given trajectory of the parameters p. The trajectories for x, p and o are always returned as the first there
    return values. If `include_reference_traj` is set to true, the function also returns additional trajectories of
    x and o for fixed parameters which in each timestep are computed from the x values of the trajectories with 
    varying parameters, allowing for a comparison of the trajectories with varying and fixed parameters. If the 
    field `random_vec_length` of the `DynamicalSystem` is not `nothing`, then the same noise is used for both the
    trajectory with variing parameters and the one with fixed parameters and the function additionaly returns the 
    noise trajectory. This is also the case if the field `random_vec_length` is not `nothing` while 
    `include_reference_traj` is set to false.
"""
function integrateTraj(DS::DynamicalSystem, p_traj::Matrix{Float64}; include_reference_traj::Bool=false)

    
    T = size(p_traj,2)
    x_arr = zeros(Float64, DS.d_x, T)
    o_arr = zeros(Float64, DS.d_o, T)
    x_arr[:,1] = DS.x_init
    o_arr[:,1] = DS.f_o(x_arr[:,1],p_traj[:,1])

    if include_reference_traj
        x_ref_arr = copy(x_arr)
        o_ref_arr = copy(o_arr)
        if DS.random_vec_length !== nothing
            noise = zeros(Float64, DS.random_vec_length, T)
            for t in 2:T
                r = randn(DS.random_vec_length)
                x_arr[:,t] = DS.f_p(x_arr[:,t-1],p_traj[:,t-1],r)
                o_arr[:,t] = DS.f_o(x_arr[:,t],p_traj[:,t])
                x_ref_arr[:,t] = DS.f_p(x_arr[:,t-1],DS.p_init,r)
                o_ref_arr[:,t] = DS.f_o(x_ref_arr[:,t],DS.p_init)
                noise[:,t] = r
            end
            return x_arr, p_traj, o_arr, x_ref_arr, o_ref_arr, noise
        else
            @warn "You are currently including a refenence trajectory with fixed parameters although the system does not allow any controll over the noise, i.e. the field `random_vec_length` of the `DynamicalSystem` is `nothing`."
            for t in 2:T
                x_arr[:,t] = DS.f_p(x_arr[:,t-1],p_traj[:,t-1])
                o_arr[:,t] = DS.f_o(x_arr[:,t],p_traj[:,t])
                x_ref_arr[:,t] = DS.f_p(x_arr[:,t-1],DS.p_init)
                o_ref_arr[:,t] = DS.f_o(x_ref_arr[:,t],DS.p_init)
            end
            return x_arr, p_traj, o_arr, x_ref_arr, o_ref_arr
        end
    else
        if DS.random_vec_length !== nothing
            noise = zeros(Float64, DS.random_vec_length, T)
            for t in 2:T
                r = randn(DS.random_vec_length)
                x_arr[:,t] = DS.f_p(x_arr[:,t-1],p_traj[:,t-1],r)
                o_arr[:,t] = DS.f_o(x_arr[:,t],p_traj[:,t])
                noise[:,t] = r
            end
            return x_arr, p_traj, o_arr, noise
        else
            for t in 2:T
                x_arr[:,t] = DS.f_p(x_arr[:,t-1],p_traj[:,t-1])
                o_arr[:,t] = DS.f_o(x_arr[:,t],p_traj[:,t])
            end
            return x_arr, p_traj, o_arr
        end
    end
end



"""
    integrateTraj(DS::DynamicalSystem, g::Function, T::Int64; include_reference_traj::Bool=false)

This function iteratively computes the evolution of the state variables x, the parameters p and the
    observables.

# Arguments
- `DS::DynamicalSystem`: Element of the type `DynamicalSystem` which encodes the dynamics of the system under
    considerantion
- `g::Function`: Function determining the one-step evolution of the parameter values. I.e. g must be such that
    it takes one Argument `p::Vector{Float64}` containing the current parameter values and returns the new 
    parameter values again in the form of `Vector{Float64}`
- `T::Int64`: Number of time steps
- `include_reference_traj::Bool=false`: If true, the reference trajectory is included in the output

# Returns
- `x_arr': Vector of dimension (DS.length_v,T) containing the time evolution of the state space variable
- `p_arr': Vector of dimension (DS.length_p,T) containing the time evolution of the parameters
- `o_arr': Vector of dimension (DS.length_x,T) containing the time evolution of the observables
- `x_ref_arr': This is only returned if `include_reference_traj==true`.
    Vector of dimension (DS.length_v,T) containing the v values which are obtained by computing
    the next timestep using the past v value as given in `v_arr` but the initial parameter values `DS.p_init`.
    If `DS.random_vec_length !== nothing` then the same noise is used for both integrations.
- `o_ref_arr': This is only returned if `include_reference_traj==true` and is obtained by applying the 
    `observable` function of the `DynamicalSystem` object to `v_ref_arr` and `DS.p_init`.
- `noise': This is only returned if `DS.random_vec_length !== nothing`.
    Vector of dimension (DS.random_vec_length,T) containing the noise values which are used for the integration
    of the reference trajectory.
"""
function integrateTraj(DS::DynamicalSystem, g::Function, T::Int64; include_reference_traj::Bool=false)
    if T <= 1
        throw(ArgumentError("The number of time steps must be larger than 1."))
    end
    try 
        g(DS.p_init)
    catch
        throw(ArgumentError("The function g must take one argument of type Vector{Float64} and return a Vector{Float64}."))
    end
    if typeof(g(DS.p_init)) != Vector{Float64}
        throw(ArgumentError("The function g must take one argument of type Vector{Float64} and return a Vector{Float64}."))
    end
    if length(g(DS.p_init)) != DS.d_p
        throw(ArgumentError("The function g must take one argument of type Vector{Float64} and return a Vector{Float64} of the same length."))
    end

    p_arr = Array(ones(Float64, 1, T) ⊗ reshape(DS.p_init, DS.d_p, 1))
    for t in 2:T
        p_arr[:,t] = g(p_arr[:,t-1])
    end
    return integrateTraj(DS, p_arr, include_reference_traj=include_reference_traj)
end


"""
    integrateTraj(DS::DynamicalSystem, T::Int64)

Integrates the trajectory of a system for T time steps with constant parameters and returns the trajectory of
    the state variables x, the parameters p and the observables. If the field `random_vec_length` of the 
    `DynamicalSystem` is not `nothing`, the function also returns the noise vector containing the noise trajectory
    used to integrate the system.
"""
function integrateTraj(DS::DynamicalSystem, T::Int64)
    g(p) = p
    return integrateTraj(DS, g, T)
end