import matplotlib.pyplot as plt
import numpy as np

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fontsize = 18

exc_ue = f'r_smooth_pulse/L8_[18.0]_[0.2]_UE.npy'
exc_lc = f'r_smooth_pulse/L8_[18.0]_[0.2]_LC.npy'
exc_smooth = f'r_smooth_pulse/L8_[18.0]_[0.2]_smooth.npy'
exc_tech = f'r_smooth_pulse/L8_[18.0]_[0.2]_tech.npy'

exc_ue = np.load(exc_ue, allow_pickle = 'true').item()
exc_lc = np.load(exc_lc, allow_pickle = 'true').item()
exc_smooth = np.load(exc_smooth, allow_pickle = 'true').item()
exc_tech = np.load(exc_tech, allow_pickle = 'true').item()

'''print(list(exc_ue.keys())); print(list(exc_lc.keys())); print(list(exc_smooth.keys())); print(list(exc_tech.keys()))'''

'''
['phi', 'Eta_sq', 'derivative']
['phi', 'Eta_sq', 'derivative', 'phi_double', 'Eta_sq_double', 'derivative_double']
['phi_1', 'Eta_sq_1', 'derivative_1', ... 'phi_6', 'Eta_sq_6', 'derivative_6']
['times', 'T', 't_center', 't_shift']
'''

times = exc_tech['times']
t_center = exc_tech['t_center']
T = exc_tech['T']
t_shift = exc_tech['t_shift']

phi_0 = exc_lc['phi_double']
phi_2 = exc_smooth['phi_2']
phi_4 = exc_smooth['phi_4']
phi_6 = exc_smooth['phi_6']
eta_0 = exc_lc['Eta_sq_double']
eta_2 = exc_smooth['Eta_sq_2']
eta_4 = exc_smooth['Eta_sq_4']
eta_6 = exc_smooth['Eta_sq_6']

fig, axes = plt.subplots(1, 2, figsize=(8, 4), constrained_layout=True)

t_start = t_center - T
t_end = t_center + T
start_idx = np.searchsorted(times, t_start, side='left')
end_idx = np.searchsorted(times, t_end, side="right")

axes[0].plot(times[start_idx:end_idx], phi_0[start_idx:end_idx], label = r'no filter')
axes[0].plot(times[start_idx:end_idx], phi_2[start_idx:end_idx], label = r'$\sigma = 2$')
axes[0].plot(times[start_idx:end_idx], phi_4[start_idx:end_idx], label = r'$\sigma = 4$')
axes[0].plot(times[start_idx:end_idx], phi_6[start_idx:end_idx], label = r'$\sigma = 6$')
axes[0].axvline(t_center - t_shift, color = 'black', linestyle = '--')
axes[0].axvline(t_center + t_shift, color = 'black', linestyle = '--')
axes[0].set_xlabel(r'$t$ [$t_h^{-1}$]', fontsize = fontsize)
axes[0].set_ylabel(r'$\Phi_{\rm{L}}(t)$', fontsize = fontsize)
axes[0].tick_params(axis = 'both', which = 'both', labelsize = fontsize, direction = 'in')
axes[0].xaxis.set_ticks_position('both')
axes[0].yaxis.set_ticks_position('both') 
axes[0].legend(fontsize = 16, frameon=False, loc='lower center', bbox_to_anchor=(0.5, -0.05), fancybox = False, handlelength = 0.75, handletextpad = 0.5)
axes[0].text(0.98, 0.02, r'(a)', transform=axes[0].transAxes, fontsize=16, verticalalignment='bottom', horizontalalignment='right')

axes[1].plot(times, eta_0/8, label = r'no filter', linewidth = 7)
axes[1].plot(times, eta_2/8, label = r'$\sigma = 2$', linewidth = 5)
axes[1].plot(times, eta_4/8, label = r'$\sigma = 4$', linewidth = 3)
axes[1].plot(times, eta_6/8, label = r'$\sigma = 6$', linewidth = 1)
axes[1].set_xlabel(r'$t$ [$t_h^{-1}$]', fontsize = fontsize)
axes[1].set_ylabel(r'$\langle \hat {\eta}^{2} \rangle / L$', fontsize = fontsize)
axes[1].tick_params(axis = 'both', which = 'both', labelsize = fontsize, direction = 'in')
axes[1].xaxis.set_ticks_position('both')
axes[1].yaxis.set_ticks_position('both') 
axes[1].legend(fontsize = 16, frameon=False, fancybox = False, handlelength = 0.75, handletextpad = 0.5)
axes[1].text(0.98, 0.02, r'(b)', transform=axes[1].transAxes, fontsize=16, verticalalignment='bottom', horizontalalignment='right')

import os
path = './analysis/r_plots'
os.makedirs(path, exist_ok=True)
fig.savefig(path + f'/exc_smooth.pdf', format='pdf', bbox_inches='tight')

plt.show()