# The methods presented here can be used to determine the timescale of a given VAR model. These also indicate
# the stability of the VAR model.

using .VARmodel_module
using LinearAlgebra
using Polynomials

"""
    toPolynomial(f::Function; deg::Int64 = 1, tolerance::Float64 = 1e-4)

Fit a polynomial of degree `deg` to the function `f`. The polynomial is fitted using the least squares method. The
    accuracy of the polynomial is checked by comparing the function `f` to the polynomial on 100 equidistant points
    in an appropriately chosen interval. If the maximum error is larger than `tolerance`, a warning is issued.
    Finally the polynomial is returned.
"""
function toPolynomial(f::Function; deg::Int64 = 1, tolerance::Float64 = 1e-4)
    # fit a polynomial of degree deg to the function f
    x = Array(LinRange(-1,1,deg+1))
    X = ones(deg+1,deg+1)
    for i in 1:deg+1
        X[:,i] = x.^(i-1)
    end
    y = X^-1 * f.(x)
    poly = Polynomial(y)

    # determine the interval relevant for checking the accuracy of the polynomial

    # At first we find the degrees of the coeficients which are larger than the coeficient of the highest degree
    idx = findall(x -> abs(x) > y[end], y[1:end-1])

    # if all coeficients are smaller than the coeficient of the highest degree, we use the interval [-3,3] since then the
    # highest order dominates the polynomial by a factor of at least 3.
    if idx == []
        test_interval = (-3,3) 

    # otherwise we determine an interval that includes [-3,3] for which the highest order dominates the polynomial at
    # the boundaries of the interval by a factor of at least 3.
    else
        x_bound = max(3, (abs(3*y[idx[end]]/y[end]))^(1/(deg-idx[end])))
        test_interval = (-x_bound,x_bound)
    end

    # check that the polynomial is accurate enough
    warn = false
    crit = zeros(2)
    for i in Array(LinRange(test_interval[1],test_interval[2],100))
        if abs(f(i) - poly(i)) > tolerance
            warn = true
            if abs(f(i) - poly(i)) > crit[2]
                crit = [i, abs(f(i) - poly(i))]
            end
        end
    end
    if warn
        @warn "Polynomial approximation is not accurate enough. Maximum error is $(crit[2]) at $(crit[1]) in [$(test_interval[1]),$(test_interval[2])]."
    end

    return poly
end


"""
    timeScale(model::VARmodel; include_error = false)

Determines a proxy of the slowest timescale of the VAR model `model` by returning the absolute value of the smallest
    root of the characteristic polynomial of the VAR model. This is done by fitting a polynomial of degree `deg` to
    the function `f(z) = det(Id - sum([z^i * B_i for i in 1:p]))` where `B_i` is the `i`-th block of the matrix `B`
    of the VAR model `model`. The root of this polynomial with the smallest absolute value (which should be larger than
    one for the model to be stable) is returned.
    Note that this function can only be applied to VAR models without parameter variables and without hidden variables.
"""
function timeScale(model::VARmodel; include_error = false)

    # check that the model is suitable for this function
    if model.d_p != 0
        error("This function can only be applied to VAR models without parameter variables.")
    end

    if include_error == false
        return timeScale(model.B)
    else
        ts = timeScale(model.B)
        sqared_error_sum = 0.
        λ, v = eigen(model.Σ_β)
        l = length(λ)
        for i in 1:l
            B_plus = reshape(vec(model.B) + √(λ[i]) .* v[:,i], size(model.B))
            B_minus = reshape(vec(model.B) - √(λ[i]) .* v[:,i], size(model.B))
            ts_plus = timeScale(B_plus)
            ts_minus = timeScale(B_minus)
            sqared_error_sum += max((ts_plus - ts)^2, (ts_minus - ts)^2)
        end
        return ts, sqrt(sqared_error_sum)
    end
end


"""
    timeScale(B::Matrix{Float64})

Auxiliary function that returns the timescale for a given matix `B`
"""
function timeScale(B::Matrix{Float64})
    
    # check that model has no parameter variables and no hidden variables
    dx = size(B,1)
    p = (size(B,2)-1)/dx
    if p != floor(p)
        error("The input matrix B does not have the right dimensions.")
    end
    p = Int64(p)

    function f(z)
        return det(Id(dx) - sum([z^i .* B[:,2+dx*(i-1):dx*i+1] for i in 1:p]))
    end

    poly = toPolynomial(f,deg=dx*p)

    rts = roots(poly)

    return min(abs.(rts)...)
end


"""
    timeScale(f_p::Function, 
                x_init::Vector{Float64}, 
                p_init::Vector{Float64}; 
                T_init::Int64 = 200,
                T_max::Int64 = 10000,
                rel_error::Float64 = 1e-1)

Determines a proxy of the slowest timescale of the dynamical system defined by the progressor function
    `f_p` with given initial state variables `x_init` and parameter variables `p_init`. 
    This is done by integrating the dynamical system for a time `T_init` and then doubling the time until
    the error on the time scale approximation is smaller than `rel_error * abs(ts - 1)`. The function then 
    returns the time scale approximation and the error on the approximation.
"""
function timeScale(f_p::Function, 
                    x_init::Vector{Float64}, 
                    p_init::Vector{Float64};
                    random_vec_length::Union{Int64,Nothing}=nothing, 
                    T_init::Int64 = 200,
                    T_max::Int64 = 10000,
                    rel_error::Float64 = 1e-1)
    
    f_o(x,p) = x
    DS = DynamicalSystem(f_p, f_o, x_init, p_init, random_vec_length = random_vec_length)
    return timeScale(DS, T_init = T_init, T_max = T_max, rel_error = rel_error)
end


"""
    timeScale(DS::DynamicalSystem; 
                    T_init::Int64 = 200,
                    T_max::Int64 = 10000,
                    rel_error::Float64 = 1e-1)

Determines a proxy of the slowest timescale of the observable dynamics of the the dynamical systen `DS`. 
    This is done by integrating the dynamical system for a time `T_init` and then doubling the time until
    the error on the time scale approximation is smaller than `rel_error * abs(ts - 1)`. The function then 
    returns the time scale approximation and the error on the approximation.
"""
function timeScale(DS::DynamicalSystem; 
                    T_init::Int64 = 200,
                    T_max::Int64 = 10000,
                    rel_error::Float64 = 1e-1)

    # check the input arguments
    if T_init <= 0
        error("T_init must be positive.")
    end
    if T_max <= T_init
        error("T_max must be larger than T_init.")
    end

    max_iter = Int64(ceil(log(T_max/T_init)/log(2)))+1
    T = T_init
    ts, err = 0., 0.
    for i in 1:max_iter
        o_tr = integrateTraj(DS, T)[3]
        model = fitVARmodel(o_tr)
        ts, err = timeScale(model, include_error = true)
        if err < abs(ts - 1) * rel_error
            return ts, err
        end
        T = 2*T
    end
    @warn "Could not determine timescale with relative error smaller than $(rel_error). For a higher accuracy, increase T_max."
    return ts, err
end