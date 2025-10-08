# TargetedCalibrationInComplexSystemModels

Code to accomplish the paper 'Targeted Calibration to Adjust Stability Biases in Non-differentiable Complex System Models'.

## Repo structure
- The folder 'src' contains multiple files that together make up the package 'TargetedCalibrationInComplexSystemModels' which provides a framework for mapping the dynamical system under investigation into the form as presented in the paper (cf. section II) and for estimating the VAR models needed in order to follow the destabilization (optimization) procedure as described in section II. The package further provides some usefull problem unspecific tools such as function for determining the stability of a given system (cf. step 2. in section II).
- The folder 'examples' displays two possibilities on how the method presented in the paper can be applied to actual systems. These systems are the same ones also discussed in section III of the paper, i.e. a Double well system and the AMOC five box model. These systems were destabilized utilizing methods from the package 'TargetedCalibrationInComplexSystemModels'.

## Recreating the plots from the paper
The code from the folder 'examples' can be used in order to recreate the plots displayed in the paper. To this end, it is necessary to activate the 'examples' julia-environment.

In a first step, the double well system and the AMOC box model system can be destabilized by the procedure discussed in section II by running 'examples/DoubleWell/DestabilizeDoubleWell.jl' and 'examples/FiveBoxModel/DestabilizeBoxModel.jl' respectively (the file 'examples/FiveBoxModel/FiveBoxModel.jl' contains some auxiliary functions and definitions simplifying the destabilization of this model without being essential to the machanics of the destabilization procedure). This will induce the creation of data-folders in which the relevant data is saved. For a detailed explanation on how the destabilization procedure was implemented we refer to section II and III of the paper and the documentation within the respective code.

In the case of the double well system the saved data can directly be plotted by running 'examples/DoubleWell/Plotting.py'. For the box model it is necessary to run 'examples/FiveBoxModel/HosingExperiment.jl' in order to generate the hystereses-data needed for the plots shown in figure 3 c and d before generating plots by running 'examples/FiveBoxModel/Plotting.py'. All plots are saved in respective folders that are automatically created.

### Analyzing the computational cost of the models
The computational cost of the models can be analyzed by running 'examples/DoubleWell/ComputationalCostAnalysis/ComputationalCost.jl' and 'examples/FiveBoxModel/ComputationalCostAnalysis/ComputationalCost.jl' respectively. This data is again saved in the respective Data folders and can be displayed in a Plot by running 'examples/DoubleWell/ComputationalCostAnalysis/PlotComputationalCost.py' and 'examples/FiveBoxModel/ComputationalCostAnalysis/PlotComputationalCost.py' respecively.

## Applying the method to custom dynamical system
- Install julia
- Copy this repository to your local device
- In the folder 'examples' create a new folder for your custom model
- Within this new folder create a '.jl' file that will be used to destabilaze your model. Within this file activate the environment 'examples' and load the library 'TargetedCalibrationInComplexSystemModels' by starting the file with the lines:
  ```
  import Pkg
  Pkg.activate("examples")
  
  using TargetedCalibrationInComplexSystemModels
  ```
- Further libraries that are needed can be added to the 'examples' environment and called within the '.jl' file.
- To set up the destabilization procedure we must:
    - define a progressor function which takes two variables if the system is deterministic (first argument describing the state of the system and second arguemt the current parameter values) and three if it is stochastic (last argument is assumed to be a vector where each component is normally distributed). Since this function produces discrete updates continuous time models must be discretised in time by introducing a finite time step 'dt'.
      ```
      # progressor function
      function f_p(x,p,r)
          ...
      end
      ```
    - define an observable function that extracts the observable of interest whos dynamics one is interested in destabilizing (optimizing) from the current state of the system and the current paramter values.
      ```
      # observable function
      function f_o(x,p)
          ...
      end
      ```
    - Initialize the corresponding dynamics system using the progressor function and the observable function defined earlier and specifying the initial state of the system as well as the initial paramter values. The keyword argument 'random_vec_length' must match the length of the third argument of 'f_p' if it exists and is neglected otherwise.
      ```
      # Initialize the dynamical system
      DS = DynamicalSystem(f_p, f_o, x_init, p_init, random_vec_length=...)
      ```
    - Implement a function 'calib_param_noise' which chooses suitable noise amplitudes used to vary the parameter values. For more detailes refer to the paper or the given examples (Double well system and the AMOC five box model).
- After the setup one can find paramter values which destabilize the system by iteratively applying the following procedure:
    - Optionally one can run the model for a transient time for it to equilibrate with respect to the current parameter values by running
      ```
      integrateTraj(DS,transient_steps)
      ```
    - Record a trajectory with fixed parameter values of length 'steps' to generate a reference model (cf. paper). For higher accuracy it makes sense to explicitely track the noise trajectory in case one is considering a stochastic model
      ```
      _,_,x_tr,noise = integrateTraj(DS,steps)
      y_tr = copy(x_tr)
      x_tr = x_tr - noise_term(noise)              # noise_term must be chosen such that 'x_tr - noise_term(noise)' equals the update given the previous value of the systems trajectory where the noise level in the last step equals zero.
      ref_model = auxFitVARmodel(Z(x_tr),Y(y_tr))
      ```
    - Determine the slowest time-scale at the current fixed point:
      ```
      ts = timeScale(DS)
      ```
    - Find suitable noise amplitudes (may depend on ts) using the function 'calib_param_noise' defined earlier:
      ```
      param_noise = calib_param_noise(...)
      ```
    - Generate a trajectory of varying parameters of length 'len':
      ```
      p_traj = parameterSeriesGenerator(ts,DS.p_init,param_noise,n=len)
      ```
    - We can now estimate the full VAR model by reintegrating the system with varying parameter values and estimating the corresponding model ('B' contains all coefficients for the first order VAR model and how it depends on the parameter values up to first order in deviations around the current parameter values. 'Σ_ε' describes the noise amplitudes of the system when modeled as first order VAR process and 'Σ_β' can be used to estimate errors on 'B' as described in the paper in full detail):
      ```
      _, _, x_arr, _, x_ref_arr = integrateTraj(DS, p_traj; include_reference_traj=true)
      # fit the parameter dependent VAR model
      B, Σ_ε, Σ_β = fitVARmodel(x_arr,x_ref_arr,p_traj,DS.p_init,ref_model)
      ```
    - Finally, using our knowlegde of how the paramters influence the Jacobian of the systems at the current fixed point up to first order in deviations from the current parameter values, we can update the parameters accordingly and repeat the process described above.
