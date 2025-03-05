# TargetedCalibrationInComplexSystemModels

Code to accomplish the paper 'Targeted Calibration to Adjust Stability Biases in Non-differentiable Complex System Models'.

## Repo structure
- The folder 'src' contains multiple files that together make up the package 'TargetedCalibrationInComplexSystemModels' which provides a framework for mapping the dynamical system under investigation into the form as presented in the paper (cf. section II) and for estimating the VAR models needed in order to follow the destabilization (optimization) procedure as described in section II. The package further provides some usefull problem unspecific tools such as function for determining the stability of a given system (cf. step 2. in section II).
- The folder 'examples' displays two possibilities on how the method presented in the paper can be applied to actual systems. These systems are the same ones also discussed in section III of the paper, i.e. a Doulbe well system and the AMOC five box model. These systems were destabilized utilizing methods from the package 'TargetedCalibrationInComplexSystemModels'.

## Recreating the plots from the paper
The code from the folder 'examples' can be used in order to recreate the plots displayed in the paper. To this end it is necessary to activate the 'examples' julia-environment.

In a first step the double well system and the AMOC box model system can be destabilized by the procedure discussed in section II by running 'examples/DoubleWell/DestabilizeDoulbeWell.jl' and 'examples/FiveBoxModel/DestabilizeBoxModel.jl' respectively (the file 'examples/FiveBoxModel/FiveBoxModel.jl' contains some auxiliary functions and definitions simplifying the destabilization of this model without being essential to the machanics of the destabilization procedure). This will induce the creation of data-folders in which the relevant data is saved. For a detailed explanation on how the destabilization procedure was implemented we refer to section II and III of the paper and the documentation within the respective code.

In the case of the double well system the saved data can directly be plotted by running 'examples/DoubleWell/Plotting.py'. For the box model it is necessary to run 'examples/FiveBoxModel/HosingExperiment.jl' in order to generate the hystereses-data needed for the plots shown in figure 3 c and d before generating plots by running 'examples/FiveBoxModel/Plotting.py'. All plots are saved in respective folders that are automatically created.

### Analyzing the computational cost of the models
The computational cost of the models can be analyzed by running 'examples/DoubleWell/ComputationalCostAnalysis/ComputationalCost.jl' and 'examples/FiveBoxModel/ComputationalCostAnalysis/ComputationalCost.jl' respectively. This data is again saved in the respective Data folders and can be displayed in a Plot by running 'examples/DoubleWell/ComputationalCostAnalysis/PlotComputationalCost.py' and 'examples/FiveBoxModel/ComputationalCostAnalysis/PlotComputationalCost.py' respecively.
