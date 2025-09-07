import numpy as np
import matplotlib.pyplot as plt
import os

threads = 8
os.environ['OMP_NUM_THREADS'] = '{}'.format(threads)
os.environ['MKL_NUM_THREADS'] = '{}'.format(threads)
os.environ['NUMEXPR_NUM_THREADS'] = '{}'.format(threads)

### | 18.0 | 0.2 | ### | 19.1 | 0.2 | ###

L = 8
omega = 18.0
amplitude = 0.2
step = 50
right = 25

exc_ue = f'r_customisation/L{L}_[{omega}]_[{amplitude}]_[{step}]_UE_[{right}].npy'
exc_lc = f'r_customisation/L{L}_[{omega}]_[{amplitude}]_[{step}]_LC_[{right}].npy'
exc_smooth = f'r_customisation/L{L}_[{omega}]_[{amplitude}]_[{step}]_smooth_[{right}].npy'
exc_tech = f'r_customisation/L{L}_[{omega}]_[{amplitude}]_[{step}]_tech_[{right}].npy'

exc_ue = np.load(exc_ue, allow_pickle = 'true').item()
exc_lc = np.load(exc_lc, allow_pickle = 'true').item()
exc_smooth = np.load(exc_smooth, allow_pickle = 'true').item()
exc_tech = np.load(exc_tech, allow_pickle = 'true').item()

times = exc_tech['times']
phi_3 = exc_smooth['phi_3']
derivative_3 = exc_smooth['derivative_3']
Eta_sq_3 = exc_smooth['Eta_sq_3']

phi_free = exc_ue['phi']

plt.figure()
plt.plot(times, derivative_3)
plt.plot(times, Eta_sq_3/L)

t_l = 5
t_r = 25
width = times[-1] - t_l - t_r

def mask_1(t):
    T = t - t_l
    phi_1 = np.sin(np.pi*T/width)**2 # 2 pi / 2 width so that only one arc (half a period)
    heaviside_l = np.heaviside(T, 0.5)
    heaviside_r = np.heaviside(T - width, 0.5)
    phi_2 = heaviside_l - heaviside_r
    return phi_1 * phi_2

def mask_2(t):
    T = t - t_l
    width_2 = width + 5
    phi_1 = np.sin(np.pi*T/width_2)**2 # 2 pi / 2 width so that only one arc (half a period)
    heaviside_l = np.heaviside(T, 0.5)
    heaviside_r = np.heaviside(T - width_2, 0.5)
    phi_2 = heaviside_l - heaviside_r
    return phi_1 * phi_2

def mask_3(t):
    T = t
    width_3 = width + 10
    phi_1 = np.sin(np.pi*T/width_3)**2 # 2 pi / 2 width so that only one arc (half a period)
    heaviside_l = np.heaviside(T, 0.5)
    heaviside_r = np.heaviside(T - width_3, 0.5)
    phi_2 = heaviside_l - heaviside_r
    return phi_1 * phi_2

instance = t_l + width/2

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

start_5 = times[-1] - 20

def mask_5(t):
    right_5 = times[-1] - 15
    width_5 = right_5 - start_5
    T = t - start_5
    phi_1 = np.cos(np.pi*T/(width_5*2))**2
    heaviside_l = np.heaviside(T, 1)
    heaviside_r = np.heaviside(T - width_5, 1)
    phi_2 = heaviside_l - heaviside_r
    phi_3 = np.heaviside(t,1) - heaviside_l
    return phi_1 * phi_2 + phi_3

start_6 = times[-1] - 15

def mask_6(t):
    right_6 = times[-1] - 10
    width_6 = right_6 - start_6
    T = t - start_6
    phi_1 = np.cos(np.pi*T/(width_6*2))**2
    heaviside_l = np.heaviside(T, 1)
    heaviside_r = np.heaviside(T - width_6, 1)
    phi_2 = heaviside_l - heaviside_r
    phi_3 = np.heaviside(t,1) - heaviside_l
    return phi_1 * phi_2 + phi_3

start_7 = times[-1] - 10

def mask_7(t):
    right_7 = times[-1] - 5
    width_7 = right_7 - start_7
    T = t - start_7
    phi_1 = np.cos(np.pi*T/(width_7*2))**2
    heaviside_l = np.heaviside(T, 1)
    heaviside_r = np.heaviside(T - width_7, 1)
    phi_2 = heaviside_l - heaviside_r
    phi_3 = np.heaviside(t,1) - heaviside_l
    return phi_1 * phi_2 + phi_3

start_8 = times[-1] - 5

def mask_8(t):
    #right_8 = times[-1] - 20
    #return np.heaviside(-t + right_8, 1)
    right_7 = times[-1] - 0
    width_7 = right_7 - start_8
    T = t - start_8
    phi_1 = np.cos(np.pi*T/(width_7*2))**2
    heaviside_l = np.heaviside(T, 1)
    heaviside_r = np.heaviside(T - width_7, 1)
    phi_2 = heaviside_l - heaviside_r
    phi_3 = np.heaviside(t,1) - heaviside_l
    return phi_1 * phi_2 + phi_3

phi_A = phi_3*mask_1(times)
phi_B = phi_3*mask_2(times)
phi_C = phi_3*mask_3(times)
phi_D = phi_3*mask_4(times)
phi_E = phi_3*mask_5(times)
phi_F = phi_3*mask_6(times)
phi_G = phi_3*mask_7(times)
phi_H = phi_3*mask_8(times)

plt.figure()
plt.plot(times, phi_3, color = 'black')
plt.plot(times, 0.25*mask_4(times), label='mask_4')
plt.plot(times, phi_D)

plt.axvline(x=start_4, color="r", linestyle="--")
plt.axvline(x=times[-1] - 20, color="r", linestyle="--")

plt.legend()

plt.figure()
plt.plot(times, phi_3, color = 'black')
plt.plot(times, 0.25*mask_5(times), label='mask_5')
plt.plot(times, phi_E)

plt.axvline(x=start_5, color="r", linestyle="--")
plt.axvline(x=times[-1] - 15, color="r", linestyle="--")

plt.legend()

plt.figure()
plt.plot(times, phi_3, color = 'black')
plt.plot(times, 0.25*mask_6(times), label='mask_6')
plt.plot(times, phi_F)

plt.axvline(x=start_6, color="r", linestyle="--")
plt.axvline(x=times[-1] - 10, color="r", linestyle="--")

plt.legend()

plt.figure()
plt.plot(times, phi_3, color = 'black')
plt.plot(times, 0.25*mask_7(times), label='mask_7')
plt.plot(times, phi_G)

plt.axvline(x=start_7, color="r", linestyle="--")
plt.axvline(x=times[-1] - 5, color="r", linestyle="--")

plt.legend()

plt.figure()
plt.plot(times, phi_3, color = 'black')
plt.plot(times, 0.25*mask_8(times), label='mask_8')
plt.plot(times, phi_H)

plt.axvline(x=start_8, color="r", linestyle="--")
plt.axvline(x=times[-1] - 0, color="r", linestyle="--")

plt.legend()

output_dict = {
    "times": times,
    "phi_A": phi_A,
    "phi_B": phi_B,
    "phi_C": phi_C,
    "phi_D": phi_D,
    "phi_E": phi_E,
    "phi_F": phi_F,
    "phi_G": phi_G,
    "phi_H": phi_H,
}

os.makedirs("r_customisation", exist_ok=True)
np.save(f'r_customisation/L{L}_[{omega}]_[{amplitude}]_PHI-2.npy', output_dict)

plt.show()