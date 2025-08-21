L = 8; omega_p_num = 18.0; A_0 = 0.2 # potentially varying parameters

s_up = L // 2 + L % 2; s_dn = L // 2
a = 1.0; pbc = True; t_h = 1.0; U_numb = 20; U = U_numb*t_h

pulse_type = 'sin_sq'; omega_p = omega_p_num*t_h
Phi_0 = a*A_0; sigma_p = 2/t_h; t_p = 10/t_h; N_p = 54; t_l_numb = 5; t_l = t_l_numb/t_h	

asymt_goal = L/2*(L/2+1)

t_s = 0.0; t_f_numb = 30; t_r_numb = 5; step_divisor = 250

import os
threads = 8
os.environ['OMP_NUM_THREADS'] = '{}'.format(threads)
os.environ['MKL_NUM_THREADS'] = '{}'.format(threads)
os.environ['NUMEXPR_NUM_THREADS'] = '{}'.format(threads)
os.makedirs("r_smooth_pulse", exist_ok=True)

from supporting import tools
from quspin.tools.measurements import obs_vs_time
from quspin.tools.evolution import evolve
from quspin.operators import hamiltonian
import numpy as np

# evolution no control #
print('no control - start')

first_dict = {}

op = tools.Operators(L, s_up, s_dn, a, pbc, t_h, U,
    pulse_type, omega_p, Phi_0, sigma_p, t_p, N_p, t_l,
    asymt_goal)
phi = op.phi
H = op.Hamiltonian('dynamic')
onsite = op.Hamiltonian('onsite')
hop_left = op.Hamiltonian('hop_left')
hop_right = op.Hamiltonian('hop_right')
Eta_sq = op.Eta_sq('full')
Q = op.Q('full')
op_dict = {'Eta_sq': Eta_sq, 'Q': Q, 'H':H}
from math import ceil
Q_norm = ceil(np.max(np.abs(op.spectrum()[0][0])))

Time = tools.Time(t_h, pulse_type, omega_p, Phi_0, N_p, t_l, t_s, t_f_numb, t_r_numb, step_divisor)
times, dt = Time.times()

psi_gs = op.spectrum()[0][1][0]
psi_t = H.evolve(psi_gs, times[0], times)
op_dict_evolution = obs_vs_time(psi_t, times, op_dict)
phi_evolution = phi(times)
Eta_sq_evolution = op_dict_evolution['Eta_sq']
Q_evolution = op_dict_evolution['Q']
derivative_evolution = t_h*np.sin(phi_evolution)*Q_evolution

last_state = np.transpose(psi_t)[-1]

first_dict['phi'] = phi_evolution
first_dict['Eta_sq'] = Eta_sq_evolution
first_dict['derivative'] = derivative_evolution
first_dict['last_state'] = last_state

np.save(f'r_smooth_pulse/L{L}_[{omega_p}]_[{Phi_0}]_[{step_divisor}]_UE.npy', first_dict) 

print('no control - end')
# evolution lyapunov control discrete #
print('lyapunov control discrete - start')

second_dict = {}

def lyapunov_evolution(t, psi):
    Q_ev = Q.expt_value(psi)
    l_ev_phi = np.arcsin(Q_ev/Q_norm)
    l_ev_one = onsite.dot(psi)
    l_ev_l = np.exp(1j * l_ev_phi)*hop_left.dot(psi)
    l_ev_r = np.exp(-1j * l_ev_phi)*hop_right.dot(psi)
    return -1j*(l_ev_one + l_ev_l + l_ev_r)

beg_in_sum = 0
end_in_sum = int(round(step_divisor,0))
derivative_sum = np.sum(derivative_evolution[beg_in_sum:end_in_sum])

psi_t_ref = np.transpose(psi_t)

time_length = len(times)

while end_in_sum < time_length and derivative_sum >= 0:
    derivative_sum -= derivative_evolution[beg_in_sum]
    derivative_sum += derivative_evolution[end_in_sum]
    beg_in_sum += 1
    end_in_sum += 1
end_in_sum -= 1

psi_t = psi_t_ref[:end_in_sum]
psi_next = np.transpose(evolve(psi_t_ref[end_in_sum], times[end_in_sum], times[end_in_sum:], lyapunov_evolution))
psi_t = np.vstack((psi_t, psi_next))

psi_t = np.transpose(psi_t)

op_dict_L = obs_vs_time(psi_t, times, op_dict)
Eta_sq_L = op_dict_L['Eta_sq']
Q_L = op_dict_L['Q']
phi_L_1 = phi_evolution[:end_in_sum]
phi_L_2 = np.arcsin(Q_L/Q_norm)[end_in_sum:]
phi_L = np.concatenate((phi_L_1,phi_L_2))
derivative_L = t_h*np.sin(phi_L)*Q_L

second_dict['phi'] = phi_L
second_dict['Eta_sq'] = Eta_sq_L
second_dict['derivative'] = derivative_L

np.save(f'r_smooth_pulse/L{L}_[{omega_p}]_[{Phi_0}]_[{step_divisor}]_LC.npy', second_dict) 

print('lyapunov control discrete - end')
# evolution lyapunov control spline #
print('lyapunov control spline - start')

from scipy.interpolate import UnivariateSpline

phi_L_spline = UnivariateSpline(times, phi_L, k=3, s=0)

def exp_phi_to_left(t):
    return np.exp(1j * phi_L_spline(t))
def exp_phi_to_right(t):
    return np.exp(-1j * phi_L_spline(t))

if L == 2: pbc = False
transition_to_left = [[-t_h,i,(i+1)] for i in range(L-1)]
transition_to_right = [[t_h,i,(i+1)] for i in range(L-1)]
if pbc:
    transition_to_left.append([-t_h,L-1,0])
    transition_to_right.append([t_h,L-1,0])
interaction = [[U,i,i] for i in range(L)]

basis = op.basis
no_checks = op.no_checks

static = [['n|n', interaction]]
dynamic = [
['+-|', transition_to_left, exp_phi_to_left, []],
['-+|', transition_to_right, exp_phi_to_right, []],
['|+-', transition_to_left, exp_phi_to_left, []],
['|-+', transition_to_right, exp_phi_to_right, []]]
H_dynamic = hamiltonian(static, dynamic, basis = basis, **no_checks)

psi_t = H_dynamic.evolve(psi_gs, times[0], times)
op_spline_evolution = obs_vs_time(psi_t, times, op_dict)
Eta_spline_evolution = op_spline_evolution['Eta_sq']
Q_spline_evolution = op_spline_evolution['Q']
phi_spline_evolution = phi_L_spline(times)
derivative_spline_evolution = t_h*np.sin(phi_spline_evolution) * Q_spline_evolution

second_dict['phi_double'] = phi_spline_evolution
second_dict['Eta_sq_double'] = Eta_spline_evolution
second_dict['derivative_double'] = derivative_spline_evolution

np.save(f'r_smooth_pulse/L{L}_[{omega_p}]_[{Phi_0}]_[{step_divisor}]_LC.npy', second_dict) 

print('lyapunov control spline - end')
# evolution lyapunov control smooth #
print('lyapunov control smooth - start')

third_dict = {}

from scipy.ndimage import gaussian_filter1d

def double_smoothing(Phi, t, t_center, t_shift, sigma, truncate=4):

    # Find indices corresponding to time boundaries
    t_start = t_center - t_shift
    t_end = t_center + t_shift
    start_idx = np.searchsorted(t, t_start, side='left')
    end_idx = np.searchsorted(t, t_end, side="right")
    
    # Calculate padding size to avoid edge effects
    cutoff = int(truncate * sigma)
    
    # Apply Gaussian filter with padding on both sides
    smooth_segment = gaussian_filter1d(
        Phi[start_idx - cutoff:end_idx + cutoff], sigma, truncate=truncate
    )
    
    middle_part = smooth_segment[cutoff:-cutoff] if cutoff else smooth_segment  # Remove padding
    filter_ = 1 - np.abs((t - t_center) / t_shift)
    filter_sup = np.abs((t - t_center) / t_shift)

    filter_ = filter_[start_idx:end_idx]
    filter_sup = filter_sup[start_idx:end_idx]
    middle_part_filter = middle_part*filter_ + Phi[start_idx:end_idx]*filter_sup

    # Reconstruct full signal: original | smoothed | original
    Phi_smooth = np.hstack(
        [
            Phi[:start_idx],  # Keep original before smoothing region
            middle_part_filter,
            Phi[end_idx:]  # Keep original after smoothing region
        ]
    )
    return Phi_smooth

t_point = times[end_in_sum]
T = 2*np.pi/omega_p
t_shift = T/2

i = 1
while i <= 6:
    phi_smooth_array = double_smoothing(phi_L, times, t_point, t_shift, sigma = i)
    phi_smooth = UnivariateSpline(times, phi_smooth_array, k=3, s=0)

    def exp_phi_to_left(t):
        return np.exp(1j * phi_smooth(t))
    def exp_phi_to_right(t):
        return np.exp(-1j * phi_smooth(t))
    
    if L == 2: pbc = False
    transition_to_left = [[-t_h,i,(i+1)] for i in range(L-1)]
    transition_to_right = [[t_h,i,(i+1)] for i in range(L-1)]
    if pbc:
        transition_to_left.append([-t_h,L-1,0])
        transition_to_right.append([t_h,L-1,0])
    interaction = [[U,i,i] for i in range(L)]

    static = [['n|n', interaction]]
    dynamic = [
    ['+-|', transition_to_left, exp_phi_to_left, []],
    ['-+|', transition_to_right, exp_phi_to_right, []],
    ['|+-', transition_to_left, exp_phi_to_left, []],
    ['|-+', transition_to_right, exp_phi_to_right, []]]
    H_dynamic = hamiltonian(static, dynamic, basis = basis, **no_checks)

    psi_t = H_dynamic.evolve(psi_gs, times[0], times)
    op_smooth_evolution = obs_vs_time(psi_t, times, op_dict)
    Eta_smooth_evolution = op_smooth_evolution['Eta_sq']
    Q_smooth_evolution = op_smooth_evolution['Q']
    phi_smooth_evolution = phi_smooth(times)
    derivative_smooth_evolution = t_h*np.sin(phi_smooth_evolution) * Q_smooth_evolution
    third_dict[f'phi_{i}'] = phi_smooth_evolution
    third_dict[f'Eta_sq_{i}'] = Eta_smooth_evolution
    third_dict[f'derivative_{i}'] = derivative_smooth_evolution

    print(f'lyapunov control smooth - {i}/6')
    
    i += 1


np.save(f'r_smooth_pulse/L{L}_[{omega_p}]_[{Phi_0}]_[{step_divisor}]_smooth.npy', third_dict)

print('lyapunov control smooth - end')
# some technical details #
print('technical - start')

fourth_dict = {}

fourth_dict['times'] = times
fourth_dict['T'] = T
fourth_dict['t_center'] = t_point
fourth_dict['t_shift'] = t_shift

np.save(f'r_smooth_pulse/L{L}_[{omega_p}]_[{Phi_0}]_[{step_divisor}]_tech.npy', fourth_dict)

print('technical - end')