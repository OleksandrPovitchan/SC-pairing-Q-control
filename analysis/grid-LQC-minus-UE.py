import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fontsize = 18
import numpy as np

dictionary_FE = np.load('./numerics/global_UE_dict.npy',allow_pickle='TRUE').item()
dictionary_QC = np.load('./numerics/global_LQC_dict.npy',allow_pickle='TRUE').item()

data_FE_final = dictionary_FE['data_final']
X = dictionary_FE['X']
Y = dictionary_FE['Y']

dataLCQ_final = dictionary_QC['data_final']
omega_array = dictionary_QC['omega_array']
Phi0_array = dictionary_QC['Phi_0']

h_array = dataLCQ_final - data_FE_final

h_array = np.real(h_array)

for i in range(h_array.shape[0]):
    for j in range(h_array.shape[1]):
        if h_array[i, j] <= 0:
            if h_array[i, j] >= -1e-5:
                h_array[i, j] = 1e-3

eta_max = np.max(h_array)


fig, (ax1) = plt.subplots(1, 1, figsize=(8, 4))

contour1 = ax1.contourf(X, Y, h_array)

plt.subplots_adjust(left = 0.1, bottom = 0.08, right = 0.825, top = 0.975)


cax = plt.axes((0.85, 0.08, 0.03, 0.895))
cbar = fig.colorbar(contour1, cax=cax)

annotate_position = eta_max
cbar.ax.annotate(f'{eta_max:.3f}', xy=(1, annotate_position), xytext=(2.5, annotate_position),
                 arrowprops=dict(arrowstyle="-", color='black'), color='black', va='center', fontsize=fontsize)
cbar.ax.tick_params(direction='in')

x = []
y = []
for i in range(h_array.shape[0]):
    for j in range(h_array.shape[1]):
        if h_array[i, j] <= -0.05:
            x.append(omega_array[j])
            y.append(omega_array[i])
ax1.scatter(x, y, color='deepskyblue', label = r'$\leq -0.05$', s = 15)

x = []
y = []
for i in range(h_array.shape[0]):
    for j in range(h_array.shape[1]):
        if h_array[i, j] <= 0.2:
            if h_array[i, j] >= 0.05:
                x.append(omega_array[j])
                y.append(omega_array[i])
ax1.scatter(x, y, color='magenta', label = r'$\in [0.05, 0.20]$', s = 2)

ax1.legend(fontsize = 16, frameon=False, loc = 'lower right', fancybox = False, labelcolor='white', handlelength = 0.75, handletextpad = 0.25)

ax1.set_xlabel(r'$\omega_p$, [$t_h$]', fontsize = fontsize)
ax1.set_ylabel(r'$\Phi_0$', fontsize = fontsize)
ax1.tick_params(top = True, right = True, axis='both', which='both', labelsize=fontsize, direction='in', color = 'white')

cbar.set_label(r'$\langle \hat{\eta}^2 \rangle / L$', fontsize=fontsize)
cbar.ax.tick_params(labelsize=fontsize)

ax1.text(0.99, 0.98, r'$(b)$', transform=ax1.transAxes, fontsize=fontsize, 
    verticalalignment='top', horizontalalignment='right', color='white')

plt.show()