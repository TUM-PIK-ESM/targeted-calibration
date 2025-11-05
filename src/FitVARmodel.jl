# Here we introduce the functions which can be used to fit a VAR model to a trajectory

using .FormatTests
using Kronecker
using .VARmodel_module
using SparseArrays




"""
    fitVARmodel(o_traj::Matrix; p_traj::Union{Matrix,Nothing}=nothing)

Fits a VAR(p) model of the form `o_traj_t = ν + A (o_traj_{t-1}', p_traj_{t-1}', (o_traj_{t-1} ⊗ p_traj_{t-1})')'
    to the given input data. It returns a `VARmodel` object which contains the estimated parameters and the 
    estimated covariance matrices.
"""
function fitVARmodel(o_traj::Matrix; p_traj::Union{Matrix,Nothing}=nothing)

    # check input data
    d_o,d_p = check_op(o_traj,p_traj)
    
    o_traj = convert(Matrix{Float64}, o_traj) 
    (p_traj === nothing) || (p_traj = convert(Matrix{Float64}, p_traj))

    # Create y_traj from o_traj and p_traj
    y_traj = create_y_traj(o_traj,p_traj=p_traj)
    d_y = size(y_traj,1)

    # set up the matrices needed for estimating the VAR model
    Z_arr = Z(o_traj)
    Y_arr = Y(y_traj)

    auxVARmodel = auxFitVARmodel(Z_arr,Y_arr)
    return VARmodel(d_o,d_p,auxVARmodel.B,auxVARmodel.Σ_ε,auxVARmodel.Σ_β,auxVARmodel.Σ_o1)
end


"""
    fitVARmodel(o_traj::Matrix, 
                o_ref_traj::Matrix, 
                p_traj::Matrix, 
                p_init::Vector, 
                fixed_param_model::VARmodel)

Fits a suitable VAR model of order 1 to the data using the method described in step 4. That is the model is 
    computed by comparing trajectories with varying and fixed parameters.
"""
function fitVARmodel(o_traj::Matrix, 
                        o_ref_traj::Matrix, 
                        p_traj::Matrix, 
                        p_init::Vector, 
                        fixed_param_model::VARmodel)
    h(o::Vector,p::Vector) = reshape(Array(o ⊗ p), length(o)*length(p))

    # check input data
    d_o,d_p = check_op(o_traj,p_traj)
    d_y_final = d_o + d_p + d_o * d_p
    check_op(o_ref_traj,p_traj)
    
    o_traj = convert(Matrix{Float64}, o_traj) 
    o_ref_traj = convert(Matrix{Float64}, o_ref_traj)
    p_traj = convert(Matrix{Float64}, p_traj)

    T = size(o_traj,2) - 1
    p_init_traj = repeat(p_init,1,T + 1)

    # Create y_traj for the model evaluating the difference between changing and non changing parameters
    y_traj = create_y_traj(o_traj,p_traj=p_traj-p_init_traj)[d_o+1:end,:]
    d_y = size(y_traj,1)

    Z_arr = Z(o_traj-o_ref_traj, T = T)
    Y_arr = Y(y_traj, T = T)[2:end,:]        # remove the ones from the first row

    diff_VARmodel = auxFitVARmodel(Z_arr,Y_arr)

    B = zeros(Float64, d_o, 1 + d_y_final)
    B[:,1] = fixed_param_model.B[:,1] - diff_VARmodel.B[:,1:d_p] * p_init
    B[:,2:1+d_o] = fixed_param_model.B[:,2:1+d_o] - diff_VARmodel.B[:,d_p+1:end] * sparse(Id(d_o) ⊗ p_init)
    B[:,2+d_o:end] = diff_VARmodel.B

    Σ_ε = fixed_param_model.Σ_ε + diff_VARmodel.Σ_ε

    Σ_β = zeros(Float64, length(B), length(B))
    Σ_β[1:(d_o+d_o^2),1:(d_o+d_o^2)] = fixed_param_model.Σ_β
    mat = sparse([-Array(Id(d_o+1)⊗p_init⊗Id(d_o)) Id(d_p*(d_o+d_o^2))])
    Σ_β += mat' * diff_VARmodel.Σ_β * mat

    return B, Σ_ε, Σ_β
end


"""
    auxFitVARmodel(Z_arr::Matrix{Float64},Y_arr::Matrix{Float64})

Fits a simple VAR model of order `p` to the data, cf. Appendix A
"""
function auxFitVARmodel(Z_arr::Matrix{Float64},Y_arr::Matrix{Float64})

    d_o = size(Z_arr,1)
    T = size(Z_arr,2)
    d_y = size(Y_arr,1) - 1

    # estimate the VAR model
    B = Z_arr * Y_arr' * (Y_arr * Y_arr')^(-1)
    Σ_ε = (Z_arr - B * Y_arr) * (Z_arr - B * Y_arr)' / (T - d_y - 1)
    Σ_β = Matrix((Y_arr * Y_arr' / T)^(-1) ⊗ Σ_ε) / T
    Σ_o1 = (T + d_y + 1) / T * Σ_ε

    return VARmodel(d_o,0,B,Σ_ε,Σ_β,Σ_o1)
end
