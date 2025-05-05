from parameters import L
fontsize = 18

import matplotlib.pyplot as plt
import numpy as np
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

read_dictionary = np.load(f'UE_vs_LC_dict_L{L}_17.0-0.4_18.5-0.4_t~o~17.0.npy',allow_pickle='TRUE').item()
read_dictionary_in = np.load(f'UE_vs_LC_dict_L{L}_19.1-0.2_18.0-0.2_t~o~19.1.npy',allow_pickle='TRUE').item()

eta_1_FE = read_dictionary['1_FE_eta']
eta_1_LC = read_dictionary['1_LC_eta']
phi_1_FE = read_dictionary['1_FE_phi']
eta_1_LC_DOUBLE = read_dictionary['1_LC_eta_DOUBLE']

eta_1_FE_in = read_dictionary_in['1_FE_eta']
eta_1_LC_in = read_dictionary_in['1_LC_eta']
phi_1_FE_in = read_dictionary_in['1_FE_phi']
eta_1_LC_DOUBLE_in = read_dictionary_in['1_LC_eta_DOUBLE']

read_dictionary_non = np.load(f'Nonloc_L{L}_17.0-0.4.npy',allow_pickle='TRUE').item()
read_dictionary_non_in = np.load(f'Nonloc_L{L}_19.1-0.2.npy',allow_pickle='TRUE').item()

times = read_dictionary_non['times']
eta_final = read_dictionary_non['eta_final']

times_in = read_dictionary_non_in['times']
eta_final_in = read_dictionary_non_in['eta_final']

fig, ax1 = plt.subplots(1, 1, figsize=(8, 4))
ax2 = ax1.twinx()
ax2.yaxis.set_label_position("left")
ax2.yaxis.tick_left()
ax1.yaxis.set_label_position("right")
ax1.yaxis.tick_right()

ax1.set_xlabel(r'$t$, $t^{\rm{(act.)}}$ [$t_h^{-1}$]', fontsize = fontsize, color = 'red')
ax2.set_ylabel(r'$\langle \hat {\eta}^{2} \rangle / L$, $\langle \hat {\eta}^{2} \rangle ^ {\rm{(sat.)}} / L $', fontsize = fontsize, color = 'red')
ax2.tick_params(axis = 'both', which = 'both', labelsize = fontsize, direction = 'in')
ax1.xaxis.set_ticks_position('both')

ax1.tick_params(axis = 'both', which = 'both', labelsize = fontsize, direction = 'in')

ph, = ax2.plot(times, eta_1_FE/L, label = r'$\omega_p = 17.0$'+'\n'+r'$\Phi_0 = 0.4$')
loc_max = np.max(eta_1_LC/L)

ax2.plot(times, eta_1_LC/L, color = ph.get_color())
ax2.plot(times, eta_1_LC_DOUBLE/L, label = 'LQC', linestyle = (0, (2, 3)), linewidth = 3, color = 'black')

ax2.plot(times, eta_final/L, color = 'red', label = r'$\rm{LQC^{(sat.)}}$')
ax2.legend(fontsize = 16, frameon=False, loc = 'upper left', fancybox = False, handlelength = 0.75, handletextpad = 0.25, labelcolor = 'red')
ax1.plot(times, phi_1_FE, color = ph.get_color(), alpha = 1, linewidth = 0.5)

ax1.set_ylabel(r'$\Phi(t)$', fontsize = fontsize, c = ph.get_color())
ax1.tick_params(axis = 'y', labelcolor=ph.get_color())


ins_ax = ax1.inset_axes([.68, .67, .3, .3])
ins_ax.xaxis.set_ticks_position('both')
ins_ax.yaxis.set_ticks_position('both')
ins_ax.tick_params(axis = 'both', which = 'both', labelsize = fontsize, 
    direction = 'in', labelbottom = False, labeltop = False, labelleft = False, labelright = False)
ph_in, = ins_ax.plot(times_in, eta_1_FE_in/L, label = r'$\omega_p = 19.1$'+'\n'+r'$\Phi_0 = 0.2$')
loc_max_in = np.max(eta_1_LC_in/L)

ins_ax.plot(times_in, eta_1_LC_in/L, color = ph.get_color())
ins_ax.plot(times_in, eta_1_LC_DOUBLE_in/L, label = 'LQC', color = 'black', linestyle = 'dashed')
ins_ax.plot(times_in, eta_final_in/L, color = 'red', label = r'$\rm{LQC^{(sat.)}}$')
ins_ax.set_xlabel(r'$t$, $t^{\rm{(act.)}}$ [$t_h^{-1}$]', fontsize = 16, loc = 'right', color = 'red')

ins_ax.text(0.51, 0.52, r'$\omega_p = 19.1$'+'\n'+r'$\Phi_0 = 0.2$', transform=ins_ax.transAxes, fontsize=16, 
    verticalalignment='top', horizontalalignment='left', color='black')

import os
path = './analysis/plots'
os.makedirs(path, exist_ok=True)
dpi=600
fig.savefig(path + f'/fig-4-a_L{L}_{dpi}.png', format='png', bbox_inches='tight', dpi=dpi)


plt.show()