from parameters import L

import matplotlib.pyplot as plt
import numpy as np
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fontsize = 18

dictionary_D = np.load(f'UE_vs_LC_dict_L{L}_18.0-0.2_19.1-0.2_t~o~18.0.npy',allow_pickle='TRUE').item()
dictionary_E = np.load(f'UE_vs_AC_dict_L{L}_18.0-0.2_19.1-0.2_t~o~18.0.npy',allow_pickle='TRUE').item()

times_1 = dictionary_D['1_times']
eta_1_FE = dictionary_D['1_FE_eta']
eta_1_LC = dictionary_D['1_LC_eta']
eta_1_LC_DOUBLE = dictionary_D['1_LC_eta_DOUBLE']
times_1_AC = dictionary_E['1_times']
eta_1_AC = dictionary_E['1_AC_eta']
eta_1_AC_DOUBLE = dictionary_E['1_AC_eta_DOUBLE']
phi_1_LC = dictionary_D['1_LC_phi']
phi_1_AC = dictionary_E['1_AC_phi']

times_2 = dictionary_D['2_times']
eta_2_FE = dictionary_D['2_FE_eta']
eta_2_LC = dictionary_D['2_LC_eta']
eta_2_LC_DOUBLE = dictionary_D['2_LC_eta_DOUBLE']
times_2_AC = dictionary_E['2_times']
eta_2_AC = dictionary_E['2_AC_eta']
eta_2_AC_DOUBLE = dictionary_E['2_AC_eta_DOUBLE']

fig, ax1 = plt.subplots(1, 1, figsize=(8, 4))
ax2 = ax1.twinx()

ax2.yaxis.set_label_position("left")
ax2.yaxis.tick_left()
ax1.yaxis.set_label_position("right")
ax1.yaxis.tick_right()

ph, = ax2.plot(times_1, eta_1_FE/L, label = r'$\omega_p = 18.0$, $\Phi_0 = 0.2$')
ax2.plot(times_1, eta_1_LC/L, color = ph.get_color())
ax2.plot(times_1_AC, eta_1_AC/L, color = ph.get_color())

ax2.plot(times_1[::400], eta_1_LC_DOUBLE[::400]/L, linestyle='', marker = 'o', markersize = 5 , color = ph.get_color())

ax2.plot(times_1_AC[100::400], eta_1_AC[100::400]/L, linestyle='', marker = '*', markersize = 10 , color = ph.get_color())


phh, = ax2.plot(times_2, eta_2_FE/L, label = r'$\omega_p = 19.1$, $\Phi_0 = 0.2$')
ax2.plot(times_2, eta_2_LC/L, color = phh.get_color())
ax2.plot(times_2_AC, eta_2_AC/L, color = phh.get_color())

ax2.plot(times_2[200::400], eta_2_LC_DOUBLE[200::400]/L, marker = 'o', markersize = 5 , color = phh.get_color(), linestyle = '', label = r'LQC')

ax2.plot(times_2_AC[300::400], eta_2_AC[300::400]/L, marker = '*', markersize = 10 , color = phh.get_color(), linestyle = '', label = r'AQC')

ax1.plot(times_1_AC, phi_1_LC, color = ph.get_color(), alpha = 1, linewidth = 0.5, label = r'$\Phi_{L}(t)$')

ax1.set_xlabel(r'$t$ [$t_h^{-1}$]', fontsize = fontsize)
ax2.set_ylabel(r'$\langle \hat {\eta}^{2} \rangle / L$', fontsize = fontsize)
ax2.tick_params(axis = 'both', which = 'both', labelsize = fontsize, direction = 'in')
ax1.xaxis.set_ticks_position('both')
ax1.set_ylabel(r'$\Phi_{\rm{L}}(t)$', fontsize = fontsize, c = ph.get_color())
ax1.tick_params(axis = 'both', which = 'both', labelsize = fontsize, direction = 'in')
ax1.tick_params(axis='y', labelcolor = ph.get_color())

ax2.legend(fontsize = 16, frameon=False, loc = 'upper left', fancybox = False, handlelength = 0.75, handletextpad = 0.25)
ax1.text(1, 0.98, r'(a)', transform=ax1.transAxes, fontsize=fontsize, 
    verticalalignment='top', horizontalalignment='right', color='black')

ins_ax = ax1.inset_axes([.01, .08, .3, .3])
ins_ax.set_title(r'$\Phi_{\rm{A}}(t)$', fontdict={'fontsize':fontsize}, loc = 'center', color = ph.get_color())
ins_ax.xaxis.set_ticks_position('both')
ins_ax.yaxis.set_ticks_position('both')
ins_ax.tick_params(axis = 'both', which = 'both', labelsize = fontsize, 
    direction = 'in', labelbottom = False, labeltop = False, labelleft = False, labelright = False)
ins_ax.plot(times_1_AC, phi_1_AC, color = ph.get_color(), alpha = 1, linewidth = 0.5, label = r'$\Phi_{A}(t)$')

import os
path = './analysis/plots'
os.makedirs(path, exist_ok=True)
dpi=600
fig.savefig(path + f'/fig-3-a_L{L}_{dpi}.png', format='png', bbox_inches='tight', dpi=dpi)

plt.show()