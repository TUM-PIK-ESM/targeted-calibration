import matplotlib.pyplot as plt
import matplotlib
import numpy as np

matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

filenames = ['ComputationalCost', 'ComputationalCost_optimised']
colors = ['black', 'red']

fig, ax1 = plt.subplots(1,1, figsize=(4,2.5))
ax2 = ax1.twinx()  
ax1.grid(alpha=0.3)
ax1.set_xlabel(r'log(Number of parameters)')
ax1.set_ylabel(r'log(Number of iterations)')

color = 'tab:red'
ax2.set_ylabel(r'log(Number of iterations)', color=color)
ax2.tick_params(axis='y', labelcolor=color)
ax2.set_ylim(6.4, 9.3)

ax_arr = [ax1, ax2]

for i, filename in enumerate(filenames):
    # Read in data
    file = open(f'examples\\FiveBoxModel\\Data\\{filename}.txt', "r")
    str = file.read()
    n_arr = [np.array(line.split('\t')).astype(float) for line in str.split('\n')[0:-1]]
    param_numb = [(i+1)*np.ones(len(n_arr[i])) for i in range(len(n_arr))]
    n_average = [np.mean(n_arr[i]) for i in range(len(n_arr))]
    param_numb_average = np.arange(1, len(n_arr)+1, 1)

    coeffs = np.polyfit(np.log(param_numb_average)[1:], np.log(n_average)[1:], 1)
    #coeffs = np.polyfit(np.log(param_numb_average), np.log(n_average), 1)

    x = np.linspace(np.log(param_numb_average[0]), np.log(param_numb_average[-1]), 100)
    ax_arr[i].plot(x, coeffs[0]*x+coeffs[1], color=colors[i], label=f'slope = {coeffs[0]:.2f}', zorder=0)
    ax_arr[i].scatter(np.log(np.concatenate(param_numb)), np.log(np.concatenate(n_arr)), color=colors[i], alpha=0.2, s=10., zorder=-1)
    ax_arr[i].scatter(np.log(param_numb_average), np.log(n_average), s=30., color=colors[i], zorder=1)
    file.close()

ax1.legend(loc='upper left', ncol=1, handlelength=1., columnspacing=0.3,frameon=False)
ax2.legend(loc='lower right', handlelength=1., columnspacing=0.3,frameon=False)
fig.tight_layout()

# save figures
fig.savefig('examples\\FiveBoxModel\\Plots\\Computational_cost_BoxModel.pdf', bbox_inches='tight')

plt.show()