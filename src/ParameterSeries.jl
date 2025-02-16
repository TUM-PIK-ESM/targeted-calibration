# Here a function is presented that allows to generate the parameter variations given the time scale of the system,
# the initial parameter values as well as the noise amplitudes.

using Random


"""
    parameterSeriesGenerator(ts::Float64, p_init::Vector{Float64}, noise::Vector{Float64}; n::Union{Int64,Nothing} = nothing)

Generate a parameter series given the time scale of the system, the initial parameter values as well as the noise amplitudes.
    The series has the form of a VAR(1) process with the given time scale and the given noise amplitudes varying around the
    given initial parameter values. The length of the parameter series is determined by the time scale if not specified.

# Arguments
- `ts::Float64`: The time scale of the system which for instance can be determined using the `timeScale` function.
- `p_init::Vector{Float64}`: The initial parameter values.
- `noise::Vector{Float64}`: The noise amplitudes.
- `n::Union{Int64,Nothing} = nothing`: The length of the parameter series. If not specified the length is determined by the time scale.
- `keep_const::Int64 = 1`: If this argument is larger than 1, the length of the parameter series is increased by a 
    factor of `keep_const` as each parameter value is kept constant for `keep_const` time steps. So the resulting
    parameter series are of the form `[p1,...,p1,p2,...,p2,...,pn,...,pn]`
"""
function parameterSeriesGenerator(ts::Float64, 
                                    p_init::Vector{Float64}, 
                                    noise::Vector{Float64}; 
                                    n::Union{Int64,Nothing} = nothing,
                                    keep_const::Int64 = 1)

    # check that the input is valid
    if length(p_init) != length(noise)
        throw(ArgumentError("The length of p_init and noise must be the same."))
    end
    if ts <= 1
        throw(ArgumentError("The time scale must be larger than 1."))
    end
    if ts > 1.5
        @warn "The time scale is larger than 1.5. Maybe the discritisation is too coarse."
    end

    # The length of the parameter series is determined by the time scale if not specified otherwise. To this end we assume that
    # for a time scale of 1/0.9 a parameter series of length 100 is sufficient for the analysis.
    if n === nothing
        n = Int(round(-log(0.9)/log(ts) * 100))
    end

    # generate the parameter series
    p_series = zeros(length(p_init),n)
    p_series[:,1] = p_init
    for i in 2:n
        p_series[:,i] = (ts - 1) / ts * p_init + 1 / ts * p_series[:,i-1] + noise .* randn(length(p_init))
    end

    if keep_const != 1
        if keep_const < 1
            throw(ArgumentError("keep_const must be larger than 1."))
        end
        p_series_new = zeros(length(p_init),Int(n*keep_const))
        for i in 1:n
            p_series_new[:,(i-1)*keep_const+1:i*keep_const] = repeat(p_series[:,i],1,keep_const)
        end
        p_series = p_series_new
    end
    return p_series
end