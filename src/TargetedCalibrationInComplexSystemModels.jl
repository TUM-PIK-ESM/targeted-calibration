module TargetedCalibrationInComplexSystemModels

export DynamicalSystem

export integrateTraj

export fitVARmodel
export auxFitVARmodel

export check_traj
export check_DynamicalSystem
export check_op

export slice_traj
export Z
export Y_t
export Y
export Id
export create_y_traj
export cut_traj

export parameterSeriesGenerator

export timeScale

export VARmodel


include("VARmodel.jl")
include("FormatTests.jl")
include("EvolveSystem.jl")
include("ParameterSeries.jl")
include("TimeScales.jl")
include("FitVARmodel.jl")
include("Matrices.jl")

end
