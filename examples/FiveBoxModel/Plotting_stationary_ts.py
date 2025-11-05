import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import os

## This python script plots the time series generated with Run_FiveBoxModel, used for Fig. 4 e.
## Added during revision of paper, November 2025, by Sebastian Bathiany.


matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'STIXGeneral'


# Read data
file = open('examples/FiveBoxModel/Data/q_FAMOUS_stationary.txt', "r")
q_FAMOUS_hosing = file.read().split('\n')[:-1]
q_FAMOUS_hosing = np.array(q_FAMOUS_hosing)
q_FAMOUS_hosing = q_FAMOUS_hosing.astype(float)
file.close()

file = open('examples/FiveBoxModel/Data/q_destabilized_stationary.txt', "r")
q_destab_hosing = file.read().split('\n')[:-1]
q_destab_hosing = np.array(q_destab_hosing)
q_destab_hosing = q_destab_hosing.astype(float)
file.close()

plt.close('all')
# Create figure 4 e from the paper

## time axis
t_max = 10000.
dt = 1.
t_arr = np.arange(0, t_max + dt, dt)

lw = 1.
fig1, ax1 = plt.subplots(1,1, figsize=(8,2))

ax1.grid(alpha=0.3)
ax1.plot(t_arr,q_FAMOUS_hosing[:len(t_arr)], color='red', linewidth=lw, label=r'using $\vec p_{init}$',zorder=0.)
ax1.plot(t_arr,q_destab_hosing[:len(t_arr)], color='black', linewidth=lw, label=r'using $\vec p_{final}$',zorder=-1)

ax1.set_xlabel('Year')
ax1.set_ylabel(r'$q$ /Sv')
#ax1.set_ylim(16, 18)
ax1.legend(loc='lower left', bbox_to_anchor=(0.0, -0.1), ncol=1, handlelength=1., columnspacing=0.3,frameon=False)

fig1.tight_layout()


# save figure
fig1.savefig('examples/FiveBoxModel/Plots/Box_model_stationary_ts.pdf', bbox_inches='tight')


plt.show()
