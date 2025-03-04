# Recreate the plot shown in figure 4 a after generating the data by running 'ComputationalCost.jl'

import matplotlib.pyplot as plt
import matplotlib
import numpy as np

matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

fig, ax = plt.subplots(1,1, figsize=(4,2.5))
ax.grid(alpha=0.3)
ax.set_xlabel(r'log(Number of parameters)')
ax.set_ylabel(r'log(Number of iterations)')

# Read in data
file = open(f'examples\\DoubleWell\\Data\\ComputationalCost.txt', "r")
str = file.read()
n_arr = [np.array(line.split('\t')).astype(float) for line in str.split('\n')[0:-1]]
param_numb = [(i+1)*np.ones(len(n_arr[i])) for i in range(len(n_arr))]
n_average = [np.mean(n_arr[i]) for i in range(len(n_arr))]
param_numb_average = np.arange(1, len(n_arr)+1, 1)

# Fit linear function to log log data neglecting the case where only one parameter is varied
coeffs = np.polyfit(np.log(param_numb_average)[1:], np.log(n_average)[1:], 1)

x = np.linspace(np.log(param_numb_average[0]), np.log(param_numb_average[-1]), 100)
ax.plot(x, coeffs[0]*x+coeffs[1], color='black', label=f'slope = {coeffs[0]:.2f}', zorder=0)
ax.scatter(np.log(np.concatenate(param_numb)), np.log(np.concatenate(n_arr)), color='black', alpha=0.2, s=10., zorder=-1)
ax.scatter(np.log(param_numb_average), np.log(n_average), s=30., color='black', zorder=1)
file.close()

fig.tight_layout()
ax.legend(loc='best', ncol=1, handlelength=1., columnspacing=0.3,frameon=False)

# save figures
fig.savefig('examples\\DoubleWell\\Plots\\Computational_cost_Doublewell.pdf', bbox_inches='tight')

plt.show()