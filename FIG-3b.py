import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import ShortTimeFFT
from scipy.signal.windows import gaussian
from scipy.ndimage import gaussian_filter1d

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fontsize = 18

#from matplotlib.colors import LogNorm

exc_ue = f'r_smooth_pulse/L8_[18.0]_[0.2]_[250]_UE.npy'
exc_lc = f'r_smooth_pulse/L8_[18.0]_[0.2]_[250]_LC.npy'
exc_smooth = f'r_smooth_pulse/L8_[18.0]_[0.2]_[250]_smooth.npy'
exc_tech = f'r_smooth_pulse/L8_[18.0]_[0.2]_[250]_tech.npy'

exc_ue = np.load(exc_ue, allow_pickle = 'true').item()
exc_lc = np.load(exc_lc, allow_pickle = 'true').item()
exc_smooth = np.load(exc_smooth, allow_pickle = 'true').item()
exc_tech = np.load(exc_tech, allow_pickle = 'true').item()

'''
['phi', 'Eta_sq', 'derivative']
['phi', 'Eta_sq', 'derivative', 'phi_double', 'Eta_sq_double', 'derivative_double']
['phi_1', 'Eta_sq_1', 'derivative_1', ... 'phi_6', 'Eta_sq_6', 'derivative_6']
['times', 'T', 't_center', 't_shift']
'''

times = exc_tech['times']
print(len(times))
t_center = exc_tech['t_center']
T = exc_tech['T']
t_shift = exc_tech['t_shift']
''' print(T/2, t_shift) '''

phi_free = exc_ue['phi']
phi_0 = exc_lc['phi_double']
phi_2 = exc_smooth['phi_2']
phi_4 = exc_smooth['phi_4']
phi_6 = exc_smooth['phi_6']
Eta_sq_6 = exc_smooth['Eta_sq_6']

t_x, N = times, len(times)
T_x = t_x[1] - t_x[0]

#########################################################

window = 4000

g_std = 1000  # standard deviation for Gaussian window in samples
w = gaussian(window, std=g_std, sym=True)  # symmetric Gaussian window
SFT = ShortTimeFFT(w, hop=50, fs=1/T_x, mfft=window*100, scale_to='magnitude')
Sx = SFT.stft(phi_6)  # perform the STFT

import numpy as np

magnitude = abs(Sx)
max_idx = magnitude.argmax(axis=0)  

o_min, o_max = 12.5, 25.0
f_min, f_max = o_min/(2*np.pi), o_max/(2*np.pi)
t_lo, t_hi, f_lo, f_hi = SFT.extent(N)
# get frequency vector 
freqs = np.linspace(f_lo, f_hi, magnitude.shape[0])  
ridge_smooth = 2*np.pi*freqs[max_idx]
time_vector = np.linspace(t_lo, t_hi, magnitude.shape[1])

mask = (time_vector >= SFT.lower_border_end[0] * SFT.T) & \
       (time_vector <= SFT.upper_border_begin(N)[0] * SFT.T)

'''plt.figure()
plt.plot(time_vector[mask], ridge_smooth[mask])
plt.axhline(y=19.1, color='r', linestyle='--', linewidth=1.5)

plt.figure()
plt.plot(t_x, phi_6)'''

fig1, ax1 = plt.subplots(figsize=(8, 2))  # enlarge plot a bit
x0=0.125
y0=0.10999999999999999
x1=0.9
y1=0.88
width = x1 - x0
height = y1 - y0
ax1.set_position([x0, y0, width, height])
t_lo, t_hi, f_lo, f_hi = SFT.extent(N)
extent=[t_lo, t_hi, 2*np.pi*f_lo, 2*np.pi*f_hi]

ax1.set_xlabel(r'$t$ [$t_h^{-1}$]', fontsize=fontsize)
ax1.set_ylabel(r"$\omega$ $[t_h]$", fontsize=fontsize)
ax1.set_ylim(o_min, o_max)
ax1.set_xlim(-1.442477796076938, 30.292033717615695)
ax1.tick_params(top = True, right = True, axis='both', which='both', labelsize=fontsize, direction='in', color = 'white')
ax1.text(1, 0.98, r'{(b)}', transform=ax1.transAxes, fontsize=fontsize, 
    verticalalignment='top', horizontalalignment='right', color='white')

im1 = ax1.imshow(abs(Sx), origin='lower', aspect='auto',
                 extent=extent, cmap='viridis')

cbar = fig1.colorbar(im1)
cbar.ax.tick_params(direction='in')
cbar.set_label(r"$|S_{\Phi_{\rm{L}}}(t, \omega)|$", fontsize=fontsize)
cbar.ax.tick_params(labelsize=fontsize)

# mark areas where window slices stick out to the side:
ax1.axvline(t_center, color='white', linestyle='--', alpha=1)
ax1.axvline(SFT.lower_border_end[0] * SFT.T, color='white', linestyle='solid', alpha=1)
ax1.axvline(SFT.upper_border_begin(N)[0] * SFT.T, color='white', linestyle='solid', alpha=1)

fig1.tight_layout()

dpi=600
import os
path = './analysis/r_plots'
os.makedirs(path, exist_ok=True)
fig1.savefig(path + f'/fig_3_b_[STFT].png', format='png', bbox_inches='tight', dpi=dpi)
fig1.savefig(path + f'/fig_3_b_[STFT].pdf', format='pdf', bbox_inches='tight')

plt.show()