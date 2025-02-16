# In this file functions are introduced, which makes the system we want to analyse acceceble for the analysis.

using .FormatTests


"""
    DynamicalSystem

An instance of the struct `DynamicalSystem` contains all the relevant information of the stochastic dynamical 
    system which we want to destabilize.

# Fields
- `f_p::Function`: In the case of a deterministic dynamical system of a stochastic dynamical system where
    we have no controll over the stochasticity,
    the Argument `f_p` should contains a function that takes two arguments, each of type
    `Vector{Float64}`, representing the current state variables `x` and the current parameter values `p` and returns
    the state variables for the next timestep as a `Vector{Float64}`. I.e. this function describes the progression
    of the dynamical system for one time step. In the case of a stochastic dynamical system, i.e. if
    `random_vec_length !== nothing` the function takes an additional argument `random_vec::Vector{Float64}` which 
    is a random vector of the type `Vector{Float64}` of length `random_vec_length`. So it makes sense to plug
    a vector of the form `randn(random_vec_length)` into the third argument of the function `f_p`.
- `f_o::Function`: This function should take the current state variables x and the current parameter values p,
    both as `Vector{Float64}`, as input and return an output of type `Vector{Float64}` containing all the observables
    of the system which are relevant for further analysis.
- `x_init::Vector{Float64}`: Initial values for the state space variables.
- `p_init::Vector{Float64}`: Initial parameter values.
- `p_bounds::Union{Vector{Tuple},Nothing}`: In case the parameter values are bounded to certain intervalls, the
    bounds of each of these intervalls can be specified in this argument. I.e. the i-th tuple in `p_bounds` 
    represent the interval bounds for the i-th parameter. If `p_bounds = nothing` then it is assumes that there
    are no restictions to the parameter values
- `d_x::Int64`: Number of state space variables
- `d_p::Int64`: Number of parameters
- `d_o::Int64`: Number of observables
- `random_vec_length::Union{Int64,Nothing}`: If the dynamical system is stochastic, i.e. if `random_vec_length !== nothing`
    then this argument specifies the length of the random vector which is used in the function `progressor` to
    describe the stochasticity of the system. If `random_vec_length = nothing` then it is assumed that the system
    is deterministic, or that we have no controll over the stochasticity of the system.
"""
mutable struct DynamicalSystem
    f_p::Function
    f_o::Function
    x_init::Vector{Float64}
    p_init::Vector{Float64}
    p_bounds::Union{Vector{Tuple},Nothing}
    d_x::Int64
    d_p::Int64
    d_o::Int64
    random_vec_length::Union{Int64,Nothing}
    function DynamicalSystem(f_p,f_o,x_init,p_init,p_bounds,d_x,d_p,d_o,random_vec_length)
        check_DynamicalSystem(f_p,f_o,x_init,p_init,p_bounds,d_x,d_p,d_o,random_vec_length)
        new(f_p,f_o,x_init,p_init,p_bounds,d_x,d_p,d_o,random_vec_length)
    end
end


"""
    DynamicalSystem(f_p,f_o,x_init,p_init;p_bounds=nothing)

Alternative initialisation of an element of type `DynamicalSystem` where `d_x`, `d_p` and 
    `d_o` are automatically computed.
"""
function DynamicalSystem(f_p::Function,
                            f_o::Function,
                            x_init::Vector{Float64},
                            p_init::Vector{Float64};
                            p_bounds::Union{Vector{Tuple},Nothing}=nothing,
                            random_vec_length::Union{Int64,Nothing}=nothing)
    check_DynamicalSystem(f_p,f_o,x_init,p_init,p_bounds,random_vec_length)
    d_x = length(x_init)
    d_p = length(p_init)
    d_o = length(f_o(x_init,p_init))
    return DynamicalSystem(f_p,f_o,x_init,p_init,p_bounds,d_x,d_p,d_o,random_vec_length)
end