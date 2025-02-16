In order to reproduce run the procedure concerning the double well system as described in the paper and 
reproduce respective plots, first run 'DestabilizeDoubleWell.jl' to generate the data. By consulting the 
comments within the code one can relate important lines of code to the respective step in the destibilization
procedure as outlined in the paper. After generating the data which is automatically stored, running
'Plotting.py' will automotically generate the respective plots similar to the ones shown in the paper. The 
script 'CalibrateParameterNoise.jl' is called within 'DestabilizeDoubleWell.jl' and is used to determine the
noise amplitudes for the parameter series by the procedure described in the paper.