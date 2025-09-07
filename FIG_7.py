import numpy as np
import matplotlib.pyplot as plt

L = 8
omega = 19.1
amplitude = 0.2
step = 50
right = 25

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fontsize = 18

exc_ue = f'r_customisation/L{L}_[{omega}]_[{amplitude}]_[{step}]_UE_[{right}].npy'
exc_lc = f'r_customisation/L{L}_[{omega}]_[{amplitude}]_[{step}]_LC_[{right}].npy'
exc_smooth = f'r_customisation/L{L}_[{omega}]_[{amplitude}]_[{step}]_smooth_[{right}].npy'
exc_tech = f'r_customisation/L{L}_[{omega}]_[{amplitude}]_[{step}]_tech_[{right}].npy'
eta_dict = f'r_customisation/L{L}_[{omega}]_[{amplitude}]_ETA-2.npy'

exc_ue = np.load(exc_ue, allow_pickle = 'true').item()
exc_lc = np.load(exc_lc, allow_pickle = 'true').item()
exc_smooth = np.load(exc_smooth, allow_pickle = 'true').item()
exc_tech = np.load(exc_tech, allow_pickle = 'true').item()
eta_dict = np.load(eta_dict, allow_pickle = 'true').item()

print(list(eta_dict.keys()))

eta_ue = exc_ue['Eta_sq']

times = exc_tech['times']
eta = exc_smooth['Eta_sq_3']
eta_4 = eta_dict["eta_4"]
eta_5 = eta_dict["eta_5"]
eta_6 = eta_dict["eta_6"]
eta_7 = eta_dict["eta_7"]
eta_8 = eta_dict["eta_8"]

phi = exc_smooth['phi_3']
phi_4 = eta_dict['phi_4']
phi_5 = eta_dict['phi_5']
phi_6 = eta_dict['phi_6']
phi_7 = eta_dict['phi_7']
phi_8 = eta_dict['phi_8']

start_4 = times[-1] - 25

def mask_4(t):
    right_4 = times[-1] - 20
    width_4 = right_4 - start_4
    T = t - start_4
    phi_1 = np.cos(np.pi*T/(width_4*2))**2
    heaviside_l = np.heaviside(T, 1)
    heaviside_r = np.heaviside(T - width_4, 1)
    phi_2 = heaviside_l - heaviside_r
    phi_3 = np.heaviside(t,1) - heaviside_l
    return phi_1 * phi_2 + phi_3

fig, axes = plt.subplots(1, 2, figsize=(8, 4), constrained_layout=True)
axes[0].plot(times, phi, color = 'black', label = r'no filter')
axes[0].plot(times, phi_8, label = r'$[-5;0]$')
axes[0].plot(times, phi_7, label = r'$[-10;-5]$')
axes[0].plot(times, phi_6, label = r'$[-15;-10]$')
axes[0].plot(times, phi_5, label = r'$[-20;-15]$')
axes[0].plot(times, phi_4, label = r'$[-25;-20]$')
axes[0].axvline(x=times[-1] - 25, color = 'black', linestyle = '--')
axes[0].axvline(x=times[-1] - 20, color = 'black', linestyle = '--')
#axes[0].plot(times, 0.25*mask_4(times), linestyle = '--', color = 'black')
                
axes[0].set_xlabel(r'$t$ [$t_h^{-1}$]', fontsize = fontsize)
axes[0].set_ylabel(r'$\Phi_{\rm{L}}(t)$', fontsize = fontsize)
axes[0].tick_params(axis = 'both', which = 'both', labelsize = fontsize, direction = 'in')
axes[0].xaxis.set_ticks_position('both')
axes[0].yaxis.set_ticks_position('both') 
axes[0].text(0.98, 0.02, r'(a)', transform=axes[0].transAxes, fontsize=16, verticalalignment='bottom', horizontalalignment='right')

start_idx = np.searchsorted(times, times[-1] - 25, side='left')

axes[1].plot(times[start_idx:],eta[start_idx:]/L,label = r'no filter', color = 'black')
axes[1].plot(times[start_idx:],eta_8[start_idx:]/L,label = r'$[-5;0]$')
axes[1].plot(times[start_idx:],eta_7[start_idx:]/L,label = r'$[-10;-5]$')
axes[1].plot(times[start_idx:],eta_6[start_idx:]/L,label = r'$[-15;-10]$')
axes[1].plot(times[start_idx:],eta_5[start_idx:]/L,label = r'$[-20;-15]$')
axes[1].plot(times[start_idx:],eta_4[start_idx:]/L,label = r'$[-25;-20]$')
axes[1].axvline(x=times[-1] - 25, color = 'black', linestyle = '--')
axes[1].axvline(x=times[-1] - 20, color = 'black', linestyle = '--')
axes[1].set_xlabel(r'$t$ [$t_h^{-1}$]', fontsize = fontsize)
axes[1].set_ylabel(r'$\langle \hat {\eta}^{2} \rangle / L$', fontsize = fontsize)
axes[1].tick_params(axis = 'both', which = 'both', labelsize = fontsize, direction = 'in')
axes[1].xaxis.set_ticks_position('both')
axes[1].yaxis.set_ticks_position('both') 
axes[1].legend(fontsize = 16, loc='lower center', bbox_to_anchor=(0.5, -0.05), labelspacing = 0.4, frameon=False, fancybox = False, handlelength = 0.75, handletextpad = 0.5)
axes[1].text(0.98, 0.02, r'(b)', transform=axes[1].transAxes, fontsize=16, verticalalignment='bottom', horizontalalignment='right')

import os
path = './analysis/r_plots'
os.makedirs(path, exist_ok=True)
fig.savefig(path + f'/full_stop_L{L}_[{omega}]_[{amplitude}].pdf', format='pdf', bbox_inches='tight')

plt.show()
