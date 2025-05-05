from parameters import L
omega_p, Phi_0 = 19.1, 0.2

import matplotlib.pyplot as plt
import numpy as np
from math import ceil
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fontsize = 18

dictionary = np.load(f'nonloc_sup_new_L{L}_omega{omega_p}_phi{Phi_0}.npy',allow_pickle='TRUE').item()

times_1 = dictionary['times']
times_2 = dictionary['times'] + times_1[-1]
phi_evolution = dictionary['phi_FE']
eta_FE_before = dictionary['eta_FE_before']
print(dictionary['eta_final'] == dictionary['eta_min'])
eta_final = dictionary['eta_final']
eta_FE_after = dictionary['eta_FE_after']
eta_LC = dictionary['eta_LC']
phi_LC = dictionary['phi_LC']
eta_LC_star = dictionary['eta_LC_star']
phi_LC_star = dictionary['phi_LC_star']
the_right_moment = dictionary['the_right_moment']

joint_times = np.concatenate((times_1,times_2[1:]))
joint_1 = np.concatenate((eta_FE_before,eta_FE_after[1:]))
phi_joint_1 = np.concatenate((phi_evolution,phi_evolution[1:]))
joint_2 = np.concatenate((eta_FE_before,eta_LC[1:]))
phi_joint_2 = np.concatenate((phi_evolution,phi_LC[1:]))
joint_3 = np.concatenate((eta_FE_before,eta_LC_star[1:]))
phi_joint_3 = np.concatenate((phi_evolution,phi_LC_star[1:]))

fig, ax1 = plt.subplots(1, 1, figsize=(8, 4))
ax2 = ax1.twinx()
ax2.yaxis.set_label_position("left")
ax2.yaxis.tick_left()
ax1.yaxis.set_label_position("right")
ax1.yaxis.tick_right()

#ax1.set_title(r'$\omega_p = 19.1$, $\Phi_0 = 0.2$, $N_p = 54$, $t_l = 5$', fontdict={'fontsize':fontsize}, loc = 'left')


pc_1, = ax2.plot(joint_times, joint_1/L, label = r'UE')
pc_2, = ax2.plot(joint_times, joint_2/L, label = r'$\rm{LQC^{(\downarrow)}}$')
pc_3, = ax2.plot(joint_times, joint_3/L, label = r'$\rm{LQC^{(\downarrow*)}}$')
ax2.legend(fontsize = 16, frameon=False, loc = 'upper left', fancybox = False, handlelength = 0.75, handletextpad = 0.25)
ax2.set_ylabel(r'$\langle \hat {\eta}^{2} \rangle / L$', fontsize = fontsize)
ax2.tick_params(axis = 'both', which = 'both', labelsize = fontsize, direction = 'in')

ax1.plot(joint_times, phi_joint_3, color = pc_3.get_color(), alpha = 1, linewidth = 0.5, label = r'$\Phi_{LQC^{*}}$')
ax1.set_xlabel(r'$t$ [$t_h^{-1}$]', fontsize = fontsize)
ax1.xaxis.set_ticks_position('both')
#ax1.set_ylabel(r'$\Phi(t)\cup\Phi_{\rm{L^{\downarrow *}}}(t)$', fontsize = fontsize, c = pc_3.get_color())
ax1.set_ylabel(r'$\Phi^{\rm{LQC^{(\downarrow *)}}}(t)$', fontsize = fontsize, c = pc_3.get_color())
ax1.tick_params(axis = 'both', which = 'both', labelsize = fontsize, direction = 'in')
ax1.tick_params(axis = 'y', labelcolor=pc_3.get_color())


ax1.text(0.03, 0.06, r'$t_l = 5$'+'\n'+r'$t_r = 5$'+'\n'+r'$N_p = 54$'+'\n'+r'$\Phi_0 = 0.2$'+'\n'+r'$\omega_p = 19.1$', transform=ax1.transAxes, fontsize=16, 
    verticalalignment='bottom', horizontalalignment='left', color='black')

time_length = len(joint_times)
time_step = time_length/48
aa = ceil(31.8*time_step)
bb = ceil(33*time_step)
phi_ax = ax1.inset_axes([.72, .67, .26, .3])
phi_ax.xaxis.set_ticks_position('both')
phi_ax.yaxis.set_ticks_position('both')
phi_ax.tick_params(axis = 'both', which = 'both', labelsize = fontsize, 
    direction = 'in', labelbottom = True, labeltop = False, labelleft = False, labelright = False)
phi_ax.plot(joint_times[aa:bb], phi_joint_3[aa:bb], alpha = 1, linewidth = 1, color = pc_3.get_color())
phi_ax.text(1, 0.925, r'(a)', transform=phi_ax.transAxes, fontsize=fontsize, 
    verticalalignment='top', horizontalalignment='right', color='black')
'''ax1.text(1, 0.98, r'(b)', transform=ax1.transAxes, fontsize=fontsize, 
    verticalalignment='top', horizontalalignment='right', color='black')'''

import os
path = './analysis/plots'
os.makedirs(path, exist_ok=True)
dpi=600
fig.savefig(path + f'/fig-5-a_L{L}_{dpi}.png', format='png', bbox_inches='tight', dpi=dpi)

figg, axx1 = plt.subplots(1, 1, figsize=(8, 4))
axx2 = axx1.twinx()
axx2.yaxis.set_label_position("left")
axx2.yaxis.tick_left()
axx1.yaxis.set_label_position("right")
axx1.yaxis.tick_right()
axx1.xaxis.set_ticks_position('both')


axx1.set_xlabel(r'$t^{\rm{(act.)}}$, $t$ [$t_h^{-1}$]', fontsize = fontsize, c = pc_1.get_color())
axx1.tick_params(axis = 'both', which = 'both', labelsize = fontsize, direction = 'in')


axx2.set_ylabel(r'$\langle \hat {\eta}^{2} \rangle ^ {\rm{(sat.)}} / L $', fontsize = fontsize)
axx2.tick_params(axis = 'both', which = 'both', labelsize = fontsize, direction = 'in')

#pcc, = axx2.plot(times_2, eta_FE_after/L, label = r'$\omega_p = 19.1$'+'\n'+r'$\Phi_0 = 0.2$')
axx2.plot(times_2, eta_final/L, color = 'red', label = r'$\rm{LQC^{(\downarrow)(sat.)}}$')
print(f'{eta_FE_after[0]/L:.2f}')
print(f'{eta_LC[0]/L:.2f}')
print(f'{eta_LC_star[0]/L:.2f}')
loc_min_1 = np.real(eta_FE_after[-1]/L)
print(f'{loc_min_1:.2f}')
loc_min_1_min = np.real(np.min(eta_FE_after/L))
print(f'({loc_min_1_min:.2f})')
loc_min_2 = np.real(np.min(eta_LC/L))
print(f'{loc_min_2:.2f}')
loc_min_3 = np.real(np.min(eta_LC_star/L))
print(f'{loc_min_3:.2f}')
axx2.axhline(loc_min_1, color = pc_1.get_color(), linestyle = (0, (2, 3)), linewidth = 3, label = r'UE')
axx2.axhline(loc_min_2, color = pc_2.get_color(), linestyle = (0, (2, 3)), linewidth = 3, label = r'$\rm{LQC^{(\downarrow)}}$')
axx2.axhline(loc_min_3, color = pc_3.get_color(), linestyle = (0, (2, 3)), linewidth = 3, label = r'$\rm{LQC^{(\downarrow*)}}$')
#axx2.axhline(loc_min, color = 'black', linestyle = 'dashed', linewidth = 0.5)
#axx2.plot(times_2, eta_LC/L, color = pcc.get_color())
#axx2.plot(times_2, eta_LC/L, label = r'$\rm{LQC^{(\downarrow)}}$', color = 'black', linestyle = 'dashed')
#axx2.axvline(times_2[the_right_moment], color = 'black', linestyle = 'dashed', linewidth = 0.5)

phcol, = axx1.plot(times_2, phi_evolution, alpha = 1, linewidth = 0.5)
axx1.set_ylabel(r'$\Phi(t)$', fontsize = fontsize, c = phcol.get_color())
axx1.tick_params(axis = 'y', labelcolor = phcol.get_color())
axx2.legend(fontsize = 16, frameon=False, loc = 'lower left', fancybox = False, handlelength = 0.75, handletextpad = 0.25)

axx1.text(1, 0.98, r'(b)', transform=axx1.transAxes, fontsize=fontsize, 
    verticalalignment='top', horizontalalignment='right', color='black')

dpi=600
figg.savefig(path + f'/fig-5-b_L{L}_{dpi}.png', format='png', bbox_inches='tight', dpi=dpi)

plt.show()