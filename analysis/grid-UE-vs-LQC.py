import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fontsize = 18
import numpy as np

modes = ['max', 'final']

dictionary_FE = np.load('./numerics/global_UE_dict.npy',allow_pickle='TRUE').item()
dictionary_QC = np.load('./numerics/global_LQC_dict.npy',allow_pickle='TRUE').item()

data_FE_max = np.real(dictionary_FE['data_max'])
data_FE_final = np.real(dictionary_FE['data_final'])
dataLCQ_final = np.real(dictionary_QC['data_final'])
X = dictionary_FE['X']
Y = dictionary_FE['Y']

for i in range(dataLCQ_final.shape[0]):
    for j in range(dataLCQ_final.shape[1]):
        if dataLCQ_final[i, j] <= 0:
            if dataLCQ_final[i, j] >= -1e-5:
                dataLCQ_final[i, j] = 1e-3

eta_max_1 = np.max(data_FE_max)
eta_max_2 = np.max(data_FE_final)
eta_max_3 = np.max(dataLCQ_final)

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 12))

levels = np.linspace(0,1.4,8)
contour1 = ax1.contourf(X, Y, data_FE_max, levels = levels)
contour2 = ax2.contourf(X, Y, data_FE_final, levels = levels)
contour3 = ax3.contourf(X, Y, dataLCQ_final, levels = levels)

plt.subplots_adjust(left = 0.1, bottom = 0.08, right = 0.825, top = 0.975, hspace = 0.2)
cax = plt.axes((0.85, 0.08, 0.03, 0.895))
cbar = fig.colorbar(contour3, cax=cax)
cbar.ax.tick_params(direction='in')
cbar.set_label(r'$\langle \hat{\eta}^2 \rangle / L$', fontsize=fontsize)
cbar.ax.tick_params(labelsize=fontsize)

cbar.ax.annotate(f'{eta_max_1:.2f}'+r'$^{\rm{(a)}}$', xy=(1, eta_max_1), xytext=(2, eta_max_1),
                 arrowprops=dict(arrowstyle="-", color='black'), color='black', va='center', fontsize=fontsize)
cbar.ax.annotate(f'{eta_max_2:.2f}'+r'$^{\rm{(b)}}$', xy=(1, eta_max_2), xytext=(2, eta_max_2),
                 arrowprops=dict(arrowstyle="-", color='black'), color='black', va='center', fontsize=fontsize)
cbar.ax.annotate(f'{eta_max_3:.2f}'+r'$^{\rm{(c)}}$', xy=(1, eta_max_3), xytext=(2, eta_max_3),
                 arrowprops=dict(arrowstyle="-", color='black'), color='black', va='center', fontsize=fontsize)

ax1.set_xlabel(r'$\omega_p$, [$t_h$]', fontsize = fontsize)
ax1.set_ylabel(r'$\Phi_0$', fontsize = fontsize)
ax1.tick_params(top = True, right = True, axis='both', which='both', labelsize=fontsize, direction='in', color = 'white')
ax2.set_xlabel(r'$\omega_p$, [$t_h$]', fontsize = fontsize)
ax2.set_ylabel(r'$\Phi_0$', fontsize = fontsize)
ax2.tick_params(top = True, right = True, axis='both', which='both', labelsize=fontsize, direction='in', color = 'white')
ax3.set_xlabel(r'$\omega_p$, [$t_h$]', fontsize = fontsize)
ax3.set_ylabel(r'$\Phi_0$', fontsize = fontsize)
ax3.tick_params(top = True, right = True, axis='both', which='both', labelsize=fontsize, direction='in', color = 'white')

ax1.text(1, 0.98, r'\bf{(a)}', transform=ax1.transAxes, fontsize=fontsize, 
    verticalalignment='top', horizontalalignment='right', color='white')
ax2.text(1, 0.98, r'\bf{(b)}', transform=ax2.transAxes, fontsize=fontsize, 
    verticalalignment='top', horizontalalignment='right', color='white')
ax3.text(1, 0.98, r'\bf{(c)}', transform=ax3.transAxes, fontsize=fontsize, 
    verticalalignment='top', horizontalalignment='right', color='white')

import os
output_path = './plots/'
os.makedirs(output_path, exist_ok=True)
dpi=600
fig.savefig(output_path + f'L8_FE_2_grids_{dpi}+.png', format='png', bbox_inches='tight', dpi=dpi)

plt.show()