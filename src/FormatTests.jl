# Here is a collection of functions that can be used to test wether certain objects have the correct format and
# produce errors if this is not the case.

module FormatTests

export check_traj
export check_DynamicalSystem
export check_op


"""
    check_traj(traj)

Checks if `traj` has the format of a trajectory i.e. `traj` is a two dimentional matrix with all entries
    given by numbers. Further it is always assumed that each column of trajectory represents a datapoint
    of the time series, i.e. each row denotes the time evolution of a single variable.
"""
function check_traj(traj)
    if !(typeof(traj) <: Matrix)
        error("You are currently treating an object as an trajectory which is not a matrix")
    end
    if any(x -> !(typeof(x) <: Number), traj)
        error("You are currently treating an object as an trajectory which does not only contain numbers")
    end
end

"""
    check_DynamicalSystem(f_p,f_o,x_init,p_init,p_bounds,x_length,p_length,o_length,random_vec_length)

Compatibility test for initialisation of an object of `DynamicalSystem`. `p_bounds` contains possible restrictions
    on the parameter values.
"""
function check_DynamicalSystem(f_p,f_o,x_init,p_init,p_bounds,x_length,p_length,o_length,random_vec_length)
    # Check if length of x_init is compatible with x_length
    if length(x_init) != x_length
        error("`x_init` is incompatible with `x_length`, as `length(x_init)` = $(length(x_init)) and `x_length` = $x_length")
    end

    # Check if length of p_init is compatible with p_length
    if length(p_init) != p_length
        error("`p_init` is incompatible with `p_length`, as `length(p_init)` = $(length(p_init)) and `x_length` = $p_length")
    end

    # Test if the function 'f_p' has the correct form
    if random_vec_length === nothing
        try 
            f_p(x_init,p_init)
        catch
            error("The function `f_p` is not of the correct form. It is supposed to take two vectors representing x and p as inputs and return x of type `Vector{Float64}` for the next timestep")
        end
        x_1 = f_p(x_init,p_init)
    else
        r = randn(random_vec_length)
        try 
            f_p(x_init,p_init,r)
        catch
            error("The function `f_p` is not of the correct form. It is supposed to take three vectors representing x, p and a random vector as inputs and return x of type `Vector{Float64}` for the next timestep")
        end
        x_1 = f_p(x_init,p_init,r)
    end
    if typeof(x_1) != Vector{Float64}
        error("The output of the function `f_p` must be of the type `Vector{Float64}`")
    end
    if length(x_1) != x_length
        error("The output of the function in the argument `f_p` must have length given by `x_length` = $x_length")
    end

    # Test if the function 'observable' has the correct form
    try
        f_o(x_init,p_init)
    catch
        error("The function `f_o` is not of the correct form. It is supposed to take two vectors representing x and p as inputs and return the resulting observables for the current timestep as `Vector{Float64}`")
    end
    o = f_o(x_init,p_init)
    if typeof(o) != Vector{Float64}
        error("The output of the function `f_o` must be of type Vector{Float64}")
    end
    if length(o) != o_length
        error("The output of the function in the argument `f_o` must have length given by `o_length` = $o_length")
    end

    # Check that p_bounds has the correct length
    if p_bounds !== nothing
        if length(p_bounds) != p_length
            error("`p_bounds` does not have the correct length. It is expected to have length given by `p_length` = $p_length")
        end
        for (i,p) in enumerate(p_init)
            if (p < p_bounds[i][1]) || (p > p_bounds[i][2])
                error("All initial parameter values given in `p_init` must respect the restrictions given by `p_bounds`")
            end
        end
    end
end


"""
    check_DynamicalSystem(f_p,f_o,x_init,p_init,p_bounds,random_vec_length)

Compatibility test for initialisation of an object of `DynamicalSystem`.
"""
function check_DynamicalSystem(f_p,f_o,x_init,p_init,p_bounds,random_vec_length)
    # Test if the function 'f_p' has the correct form
    if random_vec_length === nothing
        try 
            f_p(x_init,p_init)
        catch
            error("The function `f_p` is not of the correct form. It is supposed to take two vectors representing x and p as inputs and return x of type `Vector{Float64}` for the next timestep")
        end
        x_1 = f_p(x_init,p_init)
    else
        r = randn(random_vec_length)
        try 
            f_p(x_init,p_init,r)
        catch
            error("The function `f_p` is not of the correct form. It is supposed to take three vectors representing x, p and a random vector as inputs and return x of type `Vector{Float64}` for the next timestep")
        end
        x_1 = f_p(x_init,p_init,r)
    end
    if typeof(x_1) != Vector{Float64}
        error("The output of the function `f_p` must be of the type `Vector{Float64}`")
    end
    if length(x_1) != length(x_init)
        error("The output of the function in the argument `f_p` must have length given by `x_length` = $(length(x_init))")
    end

    # Test if the function 'f_o' has the correct form
    try
        f_o(x_init,p_init)
    catch
        error("The function `f_o` is not of the correct form. It is supposed to take two vectors representing x and p as inputs and return the resulting observables for the current timestep as `Vector{Float64}`")
    end
    o = f_o(x_init,p_init)
    if typeof(o) != Vector{Float64}
        error("The output of the function `f_o` must be of type Vector{Float64}")
    end

    # Check that p_bounds has the correct length
    if p_bounds !== nothing
        if length(p_bounds) != length(p_init)
            error("`p_bounds` does not have the correct length. It is expected to have length given by `length(p_init)` = $(length(p_init))")
        end
        for (i,p) in enumerate(p_init)
            if (p < p_bounds[i][1]) || (p > p_bounds[i][2])
                error("All initial parameter values given in `p_init` must respect the restrictions given by `p_bounds`")
            end
        end
    end
end


"""
    check_op(x_traj::Matrix,p_traj::Union{Matrix,Nothing}=nothing)

Compatibility test for simultaneously checking the trajectory `o_traj` and the parameter trajectory `p_traj`. 
    `o_traj` and `p_traj` are supposed to satify their individual compatibility tests 
    and in addition, the number of datapoints in `o_traj` and `p_traj` must be the same. 
    Besides checking the compatibility of `o_traj` and `p_traj`, the function also returns the dimensions
    of the trajectory `o_traj` and the parameter trajectory `p_traj` at the first datapoint 
    as `d_o` and `d_p`, respectively.
"""
function check_op(o_traj::Matrix,p_traj::Union{Matrix,Nothing})
    check_traj(o_traj)
    d_o = size(o_traj,1)
    if p_traj !== nothing
        check_traj(p_traj)
        d_p = size(p_traj,1)
        if size(o_traj,2) != size(p_traj,2)
            error("The number of datapoints in the trajectory `o_traj` and the parameter trajectory `p_traj` must be the same")
        end
    else
        d_p = 0
    end
    return d_o,d_p
end


end