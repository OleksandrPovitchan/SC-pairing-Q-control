# gives spectrum after UE excitement


from parameters import L

import numpy as np
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
import matplotlib.colors as mcolors

fig_a, (ax3) = plt.subplots(1, 1, figsize=(8, 4))

read_dictionary = np.load(f'UE-spectral-weights_L{L}.npy',allow_pickle='TRUE').item()
color_energy = read_dictionary['color_energy']
color_eta = read_dictionary['color_eta']
color_weight = read_dictionary['color_weight']
initial_energy = read_dictionary['initial_energy']
initial_eta = read_dictionary['initial_eta']
black_energy = read_dictionary['black_energy']
initial_weight = read_dictionary['initial_weight']
black_eta = read_dictionary['black_eta']
precision = read_dictionary['precision']
energy_final = np.real(read_dictionary['energy_final'])
eta_sq_final = np.real(read_dictionary['eta_sq_final'])
e = energy_final
etl = eta_sq_final/L

norm = mcolors.LogNorm(vmin=1e-12, vmax=1)

scatter1 = ax3.scatter(color_energy, color_eta / L, c=color_weight, s=100, norm=norm, marker="d", alpha=1)
if initial_weight < precision:
    scatter2 = ax3.scatter(initial_energy, initial_eta / L, c='black', s=100, marker="d", alpha=1)
else:
    scatter2 = ax3.scatter(initial_energy, initial_eta / L, facecolors='none', edgecolors='black', s=100, marker="d", alpha=1)
scatter3 = ax3.scatter(black_energy, black_eta / L, c='black', s=2.5, marker="o", alpha=1)

plt.subplots_adjust(left = 0.1, bottom = 0.08, right = 0.825, top = 0.975)

cax = plt.axes((0.85, 0.08, 0.03, 0.895))
cbar = fig_a.colorbar(scatter1, cax=cax)

fontsize = 18

ax3.text(1, 0.98, r'(a)', transform=ax3.transAxes, fontsize=fontsize, 
    verticalalignment='top', horizontalalignment='right', color='black')

ax3.set_xlabel(r'$\varepsilon_{m}$, [$t_h$]', fontsize = fontsize)
ax3.set_ylabel(r'$\langle \hat{\eta}^2 \rangle_{m} / L$', fontsize = fontsize)

ax3.tick_params(axis = 'both', which = 'both', labelsize = fontsize, direction = 'in')
ax3.xaxis.set_ticks_position('both')
ax3.yaxis.set_ticks_position('both')


cbar.set_label(r'$w_m$', fontsize=fontsize)
cbar.ax.tick_params(labelsize=fontsize, direction = 'in')



ent = len(color_energy)
ax3.text(0.02, 0.785, r'PBC$_{}$'+'\n'+r'$L_{} = 8$'+'\n'+r'$U_{} = 20$', ha='left', va='bottom', transform=ax3.transAxes, fontsize=16, color='black')
ax3.text(0.20, 0.780, r'$N_p = 54$' +'\n'+ r'$\Phi_0 = 0.2$' +'\n'+ r'$\omega_p = 19.1$', ha='left', va='bottom', transform=ax3.transAxes, fontsize=16, color='black')
ax3.text(0.40, 0.775, r'$\varepsilon \approx 47.48 $' + '\n' + r'$\langle \hat{\eta}^2 \rangle / L \approx 1.08 $' +'\n'+ f'{ent}/4900 states', ha='left', va='bottom', transform=ax3.transAxes, fontsize=16, color='black')

ax3.annotate(r'$\langle\hat{\eta}^2\rangle/L=0$'+'\n'+r'$\varepsilon\approx-1.12$', (initial_energy, initial_eta / L), 
             textcoords="offset points", xytext=(-5, 12), ha='left', fontsize=16, color='blue')

w_w =np.max(color_weight)
w_index = np.where(color_weight == w_w)[0][0]
w_x = np.real(color_energy[w_index])
w_y = np.real(color_eta[w_index])
print(f'{w_w:.2f}')
print(f'{w_x:.2f}')
print(f'{w_y/L:.2f}')

ax3.annotate(r'$w \approx 0.61$' + '\n' + r'$\langle\hat{\eta}^2\rangle/L = 1.5$'+'\n'+r'$\varepsilon \approx 56.28$', (w_x, w_y/L), 
             textcoords="offset points", xytext=(-75, 5), ha='left', va = 'top', fontsize=16, color='black')

index = np.where(black_eta == np.max(black_eta))[0][0]
x = black_energy[index]
y = black_eta[index] / L
ax3.annotate(r'$\langle\hat{\eta}^2\rangle/L=2.5$'+'\n'+r'$\varepsilon=80$', (x, y), 
             textcoords="offset points", xytext=(-85, -35), ha='left', fontsize=16, color='red')

import os
path = './analysis/plots'
os.makedirs(path, exist_ok=True)
dpi=600
fig_a.savefig(path + f'/fig-1-a_L{L}_{dpi}.png', format='png', bbox_inches='tight', dpi=dpi)

zero_current_E = read_dictionary['zero_current_E']
zero_current_eta = read_dictionary['zero_current_eta']
zero_current_weight = read_dictionary['zero_current_weight']

print(zero_current_weight)

plt.show()