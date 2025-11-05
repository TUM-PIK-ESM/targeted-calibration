# Here we define some matrices which are often used in order to compute the parameter estimates

using .FormatTests
using LinearAlgebra
using SparseArrays
using Kronecker

# Simple auxilary function to check whether there exists a natrual number n s.t. x = n*y 
isDivisible(x,y) = (x == (length(collect(0:y:x))-1)*y)

"""
    slice_traj(data::Matrix; T::Int64=100)

Slices out the last `T+1` columns from the `data` argument. Here p=1 denotes the number of presamples
and `T` encodes the number of actual samples. 
"""
function slice_traj(data::Matrix; T::Union{Int64,Nothing}=nothing)
    check_traj(data)
    data = convert(Matrix{Float64}, data)
    (T === nothing) && (T = size(traj)[2]-1)
    if size(data)[2] < T+1
        error("The time series given in the input data is to small (length of dataseries = $(size(data)[2])) in order to slice a trajectory with T+1=$(T+1) datapoints out of it")
    else
        return data[:,end-T:end]
    end
end


"""
    Z(traj::Matrix; T::Union{Int64,Nothing}=nothing)

This function turns a trajectory `traj`, which is a matrix with K rows, into a matrix Z as defined
    in Lütkepohls book in equation (3.2.1) (where it is called Y) or in equation (A2) of Appendix A
    by slicing the last `T` colomns from the `traj` argument. If `T===nothing` then `T` is set to the number of
    columns of `traj` minus 1 i.e. the `traj` argument is interpreted in such a way that it
    consists of all the presamples plus all the actual samples. If on the other hand `T` is given
    and `T+1` is not equal to the length of the trajectory, i.e. `traj` has less than `T+1` columns,
    an error is throuwn.
    As always we assume p=1 presamples and `T` encodes the number of actual samples. 
"""
function Z(traj::Matrix; T::Union{Int64,Nothing}=nothing)
    check_traj(traj)
    traj = convert(Matrix{Float64}, traj)
    (T === nothing) && (T = size(traj)[2]-1)
    if size(traj)[2] != T+1
        error("Your trajcetory `traj` should contain T+1=$(T+1) columns but it instead contains $(size(traj)[2]) columns")
    end
    return traj[:,end-T+1:end]
end


"""
    Y_t(traj::Matrix, t)

Returns the vector Yₜ as defined e.g. in equation (3.2.1) in Lütkepohls book (where it is called Zₜ) and in
    equation (A3) in Appendix A. Here it is assumed that the trajectory `traj` includes the presample 
    datapoints i.e. yₜ from the book corresponds to the column of `traj` with index `t+p`.
    Here we consider p=1 presamples.
"""
function Y_t(traj::Matrix, t)
    check_traj(traj)
    traj = convert(Matrix{Float64}, traj)
    if size(traj)[2] < t+1
        error("The time series encoded in the `traj` argument is not long enought (only includes $(size(traj)[2]) datapoints) as it must include at least t+1=$(t+1) datapoints")
    end
    return vcat([1], traj[:,t+1])
end


"""
    Y(traj; p::Int64=1, T::Union{Int64,Nothing}=nothing)

Returns Y as defined in Lütkepohls book in equation (3.2.1) (where it is called Z) or in equation (A4) in
    Appendix A. If the number of columns of the trajectory `traj` does not equal T+p then an error is throuwn. 
    Here we consider p=1 presamples and `T` encodes the number of actual samples. 
"""
function Y(traj; T::Union{Int64,Nothing}=nothing)
    check_traj(traj)
    (T === nothing) && (T = size(traj)[2]-1)
    if size(traj)[2] != T+1
        error("Your trajcetory `traj` should contain T+1=$(T+1) columns but it instead contains $(size(traj)[2]) columns")
    end
    d_y = size(traj)[1]
    arr = ones(Float64, d_y+1, T)
    for (idx,t) in enumerate(0:T-1)
        arr[:,idx] = Y_t(traj, t)
    end
    return arr
end



"""
    Id(d)

Returns a `d` by `d` sparse identity matrix
"""
Id(d) = sparse(Matrix{Float64}(I,d,d))



"""
    create_y_traj(o_traj::Matrix{Float64};p_traj::Union{Matrix{Float64},Nothing}=nothing)

This function creates a trajectory `y_traj` from `o_traj` and `p_traj` s.t. the i-th column of `y_traj` is given by
    `(o_traj[:,i]', p_traj[:,i]', (o_traj[:,i] ⊗ p_traj[:,i])')'`.
"""
function create_y_traj(o_traj::Matrix;p_traj::Union{Matrix,Nothing}=nothing)
    check_traj(o_traj)
    o_traj = convert(Matrix{Float64}, o_traj)
    d_o = size(o_traj)[1]
    if p_traj === nothing
        y_traj = o_traj
    else
        check_traj(p_traj)
        p_traj = convert(Matrix{Float64}, p_traj)
        d_p = size(p_traj)[1]
        y_traj = zeros(Float64, d_p + d_o + d_p * d_o, size(o_traj)[2])
        y_traj[1:d_o,:] = o_traj
        y_traj[d_o+1:d_o+d_p,:] = p_traj
        h(o::Vector,p::Vector) = reshape(Array(o ⊗ p), length(o)*length(p))
        for i in 1:size(o_traj)[2]
            y_traj[d_o+d_p+1:end,i] = h(o_traj[:,i],p_traj[:,i])
        end
    end
    return y_traj
end


"""
    cut_traj(traj::Matrix{Float64}, keep_const::Int64, cut_time_steps::Int64)

This function devides a trajectory `traj` into pieces of length `keep_const` and removes `cut_time_steps` time 
    steps from the beginning of each piece. The resulting trajectory is returned as a matrix.
"""
function cut_traj(traj::Matrix{Float64}, keep_const::Int64, cut_time_steps::Int64)
    check_traj(traj)
    traj = convert(Matrix{Float64}, traj)
    n = size(traj)[2]
    if !(isDivisible(n, keep_const))
        throw(ArgumentError("The number of time steps must be divisible by keep_const"))
    end
    if cut_time_steps >= keep_const
        throw(ArgumentError("The number of time steps to cut must be smaller than keep_const"))
    end

    # this is the result of n/keep_const
    quot = length(collect(0:keep_const:n))-1
    new_traj = zeros(Float64, size(traj)[1], n - quot*cut_time_steps)
    for i in 1:quot
        new_traj[:,(i-1)*(keep_const-cut_time_steps)+1:i*(keep_const-cut_time_steps)] = traj[:,(i-1)*keep_const+cut_time_steps+1:i*keep_const]
    end
    return new_traj
end
