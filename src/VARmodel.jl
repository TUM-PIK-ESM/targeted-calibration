# An element of the struct VARmodel contains all the information needed to define a VAR model

module VARmodel_module

export VARmodel

"""
    VARmodel

An element of the struct VARmodel contains all the information needed to define a VAR model.

# Fields
- `d_o::Int64`: dimension of the observable variable
- `d_p::Int64`: dimension of the parameter variable
- `B::Matrix{Float64}`: estimated coefficient matrix
- `Σ_ε::Matrix{Float64}`: estimated covariance matrix of the residuals (unbiased estimator)
- `Σ_β::Matrix{Float64}`: estimated covariance matrix of the coefficient matrix (unbiased estimator)
- `Σ_o1::Matrix{Float64}`: estimated MSE matrix of the one step forecast (unbiased estimator)
"""
struct VARmodel
    d_o::Int64                      # dimension of the observable variable
    d_p::Int64                      # dimension of the parameter variable
    B::Matrix{Float64}              # estimated coefficient matrix
    Σ_ε::Matrix{Float64}            # estimated covariance matrix of the residuals (unbiased estimator)
    Σ_β::Matrix{Float64}            # estimated covariance matrix of the coefficient matrix (unbiased estimator)
    Σ_o1::Matrix{Float64}           # estimated MSE matrix of the one step forecast (unbiased estimator)
end

end