import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fontsize = 18
import numpy as np

L = 8
eta_norm = L/2*(L/2 + 1)

h_start = 0.0
h_end = 60.0
h_step = 10.0

h_ranges = [f"{h_start:.1f}-{(h_start+h_step):.1f}"] 
h_ranges += [f"{(i + 0.1):.1f}-{i + h_step}" for i in np.arange(h_step, h_end, h_step)]

v_start = 0.0
v_end = 3.0
v_step = 1.5

v_ranges = [f"{v_start:.1f}-{(v_start+v_step):.1f}"] 
v_ranges += [f"{(i + 0.1):.1f}-{i + v_step}" for i in np.arange(v_step, v_end, v_step)]

mode = 'final'

i = 0
for h_tag in h_ranges:
    for v_tag in v_ranges:
        file = np.loadtxt(f'./all/L8-LQC-eta_sq_{mode}-{h_tag}-(0.1)_{v_tag}-(0.1).txt')
        if i%2 == 0: 
            v_array = file
        else: 
            a = np.concatenate((v_array,file), axis = 0)
            v_array = a
        i += 1
    if i == 2: 
        h_array = v_array
    else: 
        b = np.hstack((h_array,v_array))
        h_array = b
rows, columns = h_array.shape
h_array = np.real(h_array)
h_array = h_array/L

omega_array = np.linspace(h_start, h_end, columns)
Phi0_array = np.linspace(v_start, v_end, rows)
X, Y = np.meshgrid(omega_array, Phi0_array, indexing='xy')
eta_max = np.max(h_array)

dictionary = {
    'X':X, 'Y':Y, 'data_final':h_array, 'omega_array': omega_array, 'Phi_0': Phi0_array
    }

import os
output_path = './analysis/numerics/'
os.makedirs(output_path, exist_ok=True)
np.save(output_path + f'global_LQC_dict.npy', dictionary) 

for i in range(h_array.shape[0]):
    for j in range(h_array.shape[1]):
        if h_array[i, j] < 0:
            h_array[i, j] = 0

fig, (ax1) = plt.subplots(1, 1, figsize=(8, 4))

contour1 = ax1.contourf(X, Y, h_array)

plt.subplots_adjust(left = 0.1, bottom = 0.08, right = 0.825, top = 0.975)

cax = plt.axes((0.85, 0.08, 0.03, 0.895))
cbar = fig.colorbar(contour1, cax=cax)

annotate_position = eta_max
cbar.ax.annotate(f'{eta_max:.3f}', xy=(1, annotate_position), xytext=(2.5, annotate_position),
                 arrowprops=dict(arrowstyle="-", color='black'), color='black', va='center', fontsize=fontsize)
cbar.ax.tick_params(direction='in') 

ax1.set_xlabel(r'$\omega_p$, [$t_h$]', fontsize = fontsize)
ax1.set_ylabel(r'$\Phi_0$', fontsize = fontsize)
ax1.tick_params(top = True, right = True, axis='both', which='both', labelsize=fontsize, direction='in', color = 'white')

cbar.set_label(r'$\langle \hat{\eta}^2 \rangle / L$', fontsize=fontsize)
cbar.ax.tick_params(labelsize=fontsize)

plt.show()