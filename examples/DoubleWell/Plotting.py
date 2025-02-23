import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import os

matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

# Read data
file = open('examples\\DoubleWell\\Data\\p_arr.txt', "r")
p_arr = file.read().split('\n')[0:4]
p_arr = [line.split('\t') for line in p_arr]
p_arr = np.array(p_arr)
p_arr = p_arr.astype(float)
file.close()

file = open('examples\\DoubleWell\\Data\\lambda_arr.txt', "r")
lambda_arr = file.read().split('\n')[:-1]
lambda_arr = np.array(lambda_arr)
lambda_arr = lambda_arr.astype(float)
file.close()

file = open('examples\\DoubleWell\\Data\\lambda_err_arr.txt', "r")
lambda_err_arr = file.read().split('\n')[:-1]
lambda_err_arr = np.array(lambda_err_arr)
lambda_err_arr = lambda_err_arr.astype(float)
file.close()

file = open('examples\\DoubleWell\\Data\\x_mean_arr.txt', "r")
x_mean_arr = file.read().split('\n')[:-1]
x_mean_arr = np.array(x_mean_arr)
x_mean_arr = x_mean_arr.astype(float)
file.close()


plt.close('all')
# create figs. 2 a and b from the paper
fig1, (ax1, ax2) = plt.subplots(2,1, figsize=(4,3.5), sharex=True, gridspec_kw={'height_ratios': [1, 1]})
ax3 = ax2.twinx()  

fig1.subplots_adjust(hspace=0.1)
for ax in [ax1, ax2]:
    ax.grid(alpha=0.3)
iter_steps = np.arange(0, len(p_arr[0]), 1)
ax1.set_ylim(-1.05, 1.05)
lw = 1.

# Plots the evolution of the parameters
ax1.plot(iter_steps, p_arr[0], label=r'$p_1$', linewidth=lw, color='black')
ax1.plot(iter_steps, p_arr[1], label=r'$p_2$', linewidth=lw, color='red')
ax1.plot(iter_steps, p_arr[2], label=r'$p_3$', linewidth=lw, color='green')
ax1.plot(iter_steps, p_arr[3], label=r'$p_4$', linewidth=lw, color='blue')
ax1.set_ylabel('Parameter values')
ax1.legend(bbox_to_anchor=(0., 0.5, 0.5, 0.5), ncol=2, handlelength=1., columnspacing=1.,frameon=False)

# Plot evolution of lambda
ax2.plot(iter_steps, lambda_arr, color='black', linewidth=lw)
ax2.fill_between(iter_steps, lambda_arr-1.95*lambda_err_arr, lambda_arr+1.95*lambda_err_arr, color = 'grey', alpha=0.5, label=r'$95\%$ confidence interval of $\lambda$')
ax2.set_ylim(0.9975, 1.0001)
ax2.set_ylabel(r'$\lambda$')

# Plot the evolution of the mean value of x
color = 'tab:red'
ax3.set_ylabel(r'$\bar{x}$', color=color)
ax3.tick_params(axis='y', labelcolor=color)
ax3.plot(iter_steps, x_mean_arr, color=color, linewidth=lw)
ax3.plot(iter_steps, x_mean_arr[0]*np.ones(len(iter_steps)), color=color, linewidth=lw, ls='--', label=r'Target value of $\bar{x}$')
ax3.legend(bbox_to_anchor=(0.5, 0., 0.5, 0.5), handlelength=1., columnspacing=0.3,frameon=False)
ax3.set_ylim(0.983, 1.007)
ax3.grid(alpha=0.3)

ax2.set_xlabel('Iteration steps')
ax2.legend(loc='lower right', ncol=2, handlelength=1.5, columnspacing=0.3,frameon=False)
fig1.tight_layout()


# Create figure 2 c from the paper
fig2, ax4 = plt.subplots(1,1, figsize=(4,2.5))

x = np.linspace(-1.5, 1.5, 1000)
ax4.grid(alpha=0.3)
ax4.plot(x,p_arr[0,-1]*x**3 + p_arr[1,-1]*x**2 + p_arr[2,-1]*x + p_arr[3,-1], color='black', linewidth=lw, label=r'$p_1^{final} \, x^3 + p_2^{final} \, x^2 + p_3^{final} \, x + p_4^{final}$',zorder=0.)
ax4.plot(x,p_arr[0,0]*x**3 + p_arr[1,0]*x**2 + p_arr[2,0]*x + p_arr[3,0], color='red', linewidth=lw, label=r'$p_1^{init} \, x^3 + p_2^{init} \, x^2 + p_3^{init} \, x + p_4^{init}$',zorder=-1)
ax4.plot([1.,1.],[-2.,2.], color='black', linewidth=lw, ls='--',zorder=-2, alpha = 0.3)
ax4.set_ylim(-2.05, 2.05)
ax4.set_ylabel('RHS of ODE')
ax4.set_xlabel(r'$x$')
ax4.legend(loc='upper right', ncol=1, handlelength=1., columnspacing=0.3,frameon=False)

fig2.tight_layout()

# save figures
if not os.path.isdir("examples\\DoubleWell\\Plots"):
    os.mkdir("examples\\DoubleWell\\Plots")
fig1.savefig('examples\\DoubleWell\\Plots\\Double_well_iterations.pdf', bbox_inches='tight')
fig2.savefig('examples\\DoubleWell\\Plots\\Double_well_RHS.pdf', bbox_inches='tight')

plt.show()