# gives envelopes of pulses and evolution of the order parameter (UE)

L = 8

import numpy as np
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

dictionary_1 = f'UE_dict_L{L}_o19.1_Phi0.2_t~o~18.0.npy'
read_dictionary_1 = np.load(dictionary_1, allow_pickle='TRUE').item()
dictionary_2 = f'UE_dict_L{L}_o18.0_Phi0.2_t~o~18.0.npy'
read_dictionary_2 = np.load(dictionary_2, allow_pickle='TRUE').item()
dictionary_3 = f'UE_dict_L{L}_o19.1_Phi0.4_t~o~18.0.npy'
read_dictionary_3 = np.load(dictionary_3, allow_pickle='TRUE').item()

fontsize = 18

fig_b, (ax1) = plt.subplots(1, 1, figsize=(8, 4))
ax2 = ax1.twinx()
ax2.yaxis.set_label_position("left")
ax2.yaxis.tick_left()
ax1.yaxis.set_label_position("right")
ax1.yaxis.tick_right()

ax2.set_ylabel(r'$\langle \hat {\eta}^{2} \rangle / L$', fontsize = fontsize)
ax2.tick_params(axis = 'both', which = 'both', labelsize = fontsize, direction = 'in')
ax1.set_xlabel(r'$t$ [$t_h^{-1}$]', fontsize = fontsize)
ax1.xaxis.set_ticks_position('both')
ax1.set_ylabel(r'$\Phi(t)$', fontsize = fontsize)
ax1.tick_params(axis = 'both', which = 'both', labelsize = fontsize, direction = 'in')

ax1.text(1, 0.98, r'(b)', transform=ax1.transAxes, fontsize=fontsize, 
    verticalalignment='top', horizontalalignment='right', color='black')

times1 = read_dictionary_1['times']
times2 = read_dictionary_2['times']
times3 = read_dictionary_3['times']
y1 = read_dictionary_1['Eta_sq_evolution']/L
y2 = read_dictionary_2['Eta_sq_evolution']/L
y3 = read_dictionary_3['Eta_sq_evolution']/L

L = 8
s_up = L // 2 + L % 2
s_dn = L // 2
from supporting import tools
from parameters import a, pbc, t_h, U, N_p, t_l, asymt_goal
op = tools.Operators(L, s_up, s_dn, a, pbc, t_h, U, 
        'sin_sq_envelope', 19.1, 0.2, 0, 0, N_p, t_l, asymt_goal)
phi1 = op.phi(times1)
op = tools.Operators(L, s_up, s_dn, a, pbc, t_h, U, 
        'sin_sq_envelope', 18.0, 0.2, 0, 0, N_p, t_l, asymt_goal)
phi2 = op.phi(times2)
op = tools.Operators(L, s_up, s_dn, a, pbc, t_h, U, 
        'sin_sq_envelope', 19.1, 0.4, 0, 0, N_p, t_l, asymt_goal)
phi3 = op.phi(times3)

ax2.plot(times2, y2, label = r'$\omega_p = 18.0$, $\Phi_0 = 0.2$')
ax2.plot(times1, y1, label = r'$\omega_p = 19.1$, $\Phi_0 = 0.2$')
ax2.plot(times3, y3, label = r'$\omega_p = 19.1$'+'\n'+'$\Phi_0 = 0.4$')

en2, = ax1.plot(times2, phi2, alpha = 1, linestyle = 'dotted')
en1, = ax1.plot(times1, phi1, alpha = 1, linestyle = 'dotted')
en3, = ax1.plot(times3, phi3, alpha = 1, linestyle = 'dotted')
ax1.plot(times2, -phi2, alpha = 1, linestyle = 'dotted', color = en2.get_color())
ax1.plot(times1, -phi1, alpha = 1, linestyle = 'dotted', color = en1.get_color())
ax1.plot(times3, -phi3, alpha = 1, linestyle = 'dotted', color = en3.get_color())

plt.subplots_adjust(left = 0.1, bottom = 0.15, right = 0.9, top = 0.975)

ax2.legend(fontsize = 16, frameon=False, loc = 'upper left', fancybox = False, handlelength = 0.75, handletextpad = 0.25)

import os
path = './analysis/plots'
os.makedirs(path, exist_ok=True)
dpi=600
fig_b.savefig(path + f'/fig-1-b_L{L}_{dpi}.png', format='png', bbox_inches='tight', dpi=dpi)

plt.show()