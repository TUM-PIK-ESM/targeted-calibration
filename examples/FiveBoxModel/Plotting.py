import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import os

matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'STIXGeneral'


# Read data
file = open('examples\\FiveBoxModel\\Data\\p_arr.txt', "r")
p_arr = file.read().split('\n')[0:5]
p_arr = [line.split('\t') for line in p_arr]
p_arr = np.array(p_arr)
p_arr = p_arr.astype(float)
file.close()

file = open('examples\\FiveBoxModel\\Data\\lambda_arr.txt', "r")
lambda_arr = file.read().split('\n')[:-1]
lambda_arr = np.array(lambda_arr)
lambda_arr = lambda_arr.astype(float)
file.close()

file = open('examples\\FiveBoxModel\\Data\\lambda_err_arr.txt', "r")
lambda_err_arr = file.read().split('\n')[:-1]
lambda_err_arr = np.array(lambda_err_arr)
lambda_err_arr = lambda_err_arr.astype(float)
file.close()

file = open('examples\\FiveBoxModel\\Data\\q_mean_arr.txt', "r")
q_mean_arr = file.read().split('\n')[:-1]
q_mean_arr = np.array(q_mean_arr)
q_mean_arr = q_mean_arr.astype(float)
file.close()

file = open('examples\\FiveBoxModel\\Data\\Hosing.txt', "r")
hosing = file.read().split('\n')[:-1]
hosing = np.array(hosing)
hosing = hosing.astype(float)
file.close()

file = open('examples\\FiveBoxModel\\Data\\q_destabilized_during_hosing.txt', "r")
q_destab_hosing = file.read().split('\n')[:-1]
q_destab_hosing = np.array(q_destab_hosing)
q_destab_hosing = q_destab_hosing.astype(float)
file.close()

file = open('examples\\FiveBoxModel\\Data\\q_FAMOUS_during_hosing.txt', "r")
q_FAMOUS_hosing = file.read().split('\n')[:-1]
q_FAMOUS_hosing = np.array(q_FAMOUS_hosing)
q_FAMOUS_hosing = q_FAMOUS_hosing.astype(float)
file.close()


plt.close('all')
# Create figure 3 a and b from the paper
fig1, (ax1, ax2) = plt.subplots(2,1, figsize=(4,4), sharex=True, gridspec_kw={'height_ratios': [1.5, 1]})
ax3 = ax2.twinx()  

fig1.subplots_adjust(hspace=0.1)
for ax in [ax1, ax2]:
    ax.grid(alpha=0.3)
iter_steps = np.arange(0, len(p_arr[0]), 1)
ax1.set_ylim(-0.1, 2.3)
lw = 1.
# plot the time evolution of the parameter values that where varied during the destabilization process
ax1.plot(iter_steps, p_arr[0]*3, label=r'$3 \cdot \gamma$', linewidth=lw, color='black')
ax1.plot(iter_steps, p_arr[1]/100 * 3, label=r'$0.03 \cdot \eta$ /Sv', linewidth=lw, color='orange')
ax1.plot(iter_steps, p_arr[0], label=" ", alpha = 0.) # invisible line for legend arrangement
ax1.plot(iter_steps, p_arr[2]/10 *2, label=r'$0.2 \cdot K_N$ /Sv', linewidth=lw, color='green')
ax1.plot(iter_steps, p_arr[3]/10, label=r'$ 0.1 \cdot K_S$ /Sv', linewidth=lw, color='blue')
ax1.plot(iter_steps, p_arr[4]/100, label=r'$ 0.01 \cdot K_{IP}$ /Sv', linewidth=lw, color='red')
ax1.set_ylabel('Parameter values')
ax1.legend(loc='upper right', ncol=2, handlelength=1., columnspacing=1.,frameon=False)

# Plot the evolution of the lambda values during the destabilization process
ax2.plot(iter_steps, lambda_arr, color='black', linewidth=lw)
ax2.fill_between(iter_steps, lambda_arr-1.95*lambda_err_arr, lambda_arr+1.95*lambda_err_arr, color = 'grey', alpha=0.5, label=r'$95\%$ confidence interval of $\lambda$')
ax2.set_ylim(0.9984, 1.0001)
ax2.set_ylabel(r'$\lambda$')

# Plot the evolution of the AMOC strength during the destabilization process
color = 'tab:red'
ax3.set_ylabel(r'$\bar{q}$ /Sv', color=color)
ax3.tick_params(axis='y', labelcolor=color)
ax3.plot(iter_steps, q_mean_arr, color=color, linewidth=lw)
ax3.plot(iter_steps, q_mean_arr[0]*np.ones(len(iter_steps)), color=color, linewidth=lw, ls='--', label=r'Target value of $\bar{q}$ /Sv')
ax3.legend(bbox_to_anchor=(0.5, 0., 0.5, 0.5), handlelength=1., columnspacing=0.3,frameon=False)
ax3.set_ylim(15.8, 17.4)

ax2.set_xlabel('Iteration steps')
ax2.legend(loc='lower right', ncol=2, handlelength=1.5, columnspacing=0.3,frameon=False)
fig1.tight_layout()


# Create figures 3 c and d from the paper
fig2, (ax4, ax5) = plt.subplots(2,1, figsize=(4,4))

max_idx = np.min(np.argmin(np.abs(hosing - 0.5)))

for ax in [ax4, ax5]:
    ax.grid(alpha=0.3)
ax4.plot(hosing[:max_idx],q_FAMOUS_hosing[:max_idx], color='red', linewidth=lw, label=r'using $\vec p_{init}$',zorder=0.)
ax4.plot(hosing[:max_idx],q_destab_hosing[:max_idx], color='black', linewidth=lw, label=r'using $\vec p_{final}$',zorder=-1)
ax4.scatter(hosing[0:1], q_FAMOUS_hosing[0:1], color='black', marker='o', s=20, label='Initial state', zorder=1.)
ax4.plot(hosing[:max_idx], np.zeros(len(hosing[:max_idx])), color='black', linewidth=lw, ls='--',alpha=0.3)
ax4.set_xlabel('Hosing strength /Sv')
ax4.set_ylabel(r'$q$ /Sv')
ax4.set_ylim(-11, 21)
ax4.legend(loc='lower left', ncol=1, handlelength=1., columnspacing=0.3,frameon=False)

ax5.plot(hosing,q_FAMOUS_hosing, color='red', linewidth=lw, label=r'using $\vec p_{init}$',zorder=0.)
ax5.plot(hosing,q_destab_hosing, color='black', linewidth=lw, label=r'using $\vec p_{final}$',zorder=-1)
ax5.scatter(hosing[0:1], q_FAMOUS_hosing[0:1], color='black', marker='o', s=20, label='Initial state', zorder=1.)
ax5.set_ylabel(r'$q$ /Sv')
ax5.set_xlabel('Hosing strength /Sv')
ax5.legend(loc='lower left', ncol=1, handlelength=1., columnspacing=0.3,frameon=False)

fig2.tight_layout()

# save figures
if not os.path.isdir("examples\\FiveBoxModel\\Plots"):
    os.mkdir("examples\\FiveBoxModel\\Plots")
fig1.savefig('examples\\FiveBoxModel\\Plots\\Box_model_iterations.pdf', bbox_inches='tight')
fig2.savefig('examples\\FiveBoxModel\\Plots\\Box_model_hosing.pdf', bbox_inches='tight')

plt.show()