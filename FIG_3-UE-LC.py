# WARNING!
# {LINE 28 & LINE 152} - angular frequency and amplitude are introduced in the code!


import numpy as np

from parameters import L, s_up, s_dn, a, pbc, t_h, U
from parameters import pulse_type, sigma_p, t_p, N_p, t_l
from parameters import asymt_goal
from supporting import tools
from quspin.tools.measurements import obs_vs_time
from quspin.tools.evolution import evolve
from scipy.interpolate import UnivariateSpline
from quspin.operators import hamiltonian

output_dict = {}

import os
threads = 8
os.environ['OMP_NUM_THREADS'] = '{}'.format(threads)
os.environ['MKL_NUM_THREADS'] = '{}'.format(threads)
os.environ['NUMEXPR_NUM_THREADS'] = '{}'.format(threads)

# 
# start of the curve 1 
# 

omega_p, Phi_0 = 18.0, 0.2

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
op_double = {'Eta_sq': Eta_sq}
from math import ceil
Q_norm = ceil(np.max(op.spectrum()[0][0]))

from parameters import t_s, t_f_numb, t_r_numb, step_divisor
Time = tools.Time(t_h, pulse_type, omega_p, Phi_0, N_p, t_l, t_s, t_f_numb, t_r_numb, step_divisor)
t_s = Time.t_s
t_f = Time.t_f
times, dt = Time.times()
time_length = len(times)

output_dict['1_times'] = times

psi_gs = op.spectrum()[0][1][0]
psi_t = H.evolve(psi_gs, times[0], times)
op_dict_evolution = obs_vs_time(psi_t, times, op_dict)
phi_evolution = phi(times)
Eta_sq_evolution = op_dict_evolution['Eta_sq']
Q_evolution = op_dict_evolution['Q']
H_evolution = op_dict_evolution['H']
derivative_evolution = t_h*np.sin(phi_evolution)*Q_evolution

output_dict['1_FE_phi'] = phi_evolution
output_dict['1_FE_H'] = H_evolution
output_dict['1_FE_eta'] = Eta_sq_evolution
output_dict['1_FE_Q'] = Q_evolution
output_dict['1_FE_derivative'] = derivative_evolution

print('1st 5 are written')

beg_in_sum = 0
end_in_sum = int(round(step_divisor,0))
derivative_sum = np.sum(derivative_evolution[beg_in_sum:end_in_sum])

def lyapunov_evolution(t, psi):
    Q_ev = Q.expt_value(psi)
    l_ev_phi = np.arcsin(Q_ev/Q_norm)
    l_ev_one = onsite.dot(psi)
    l_ev_l = np.exp(1j * l_ev_phi)*hop_left.dot(psi)
    l_ev_r = np.exp(-1j * l_ev_phi)*hop_right.dot(psi)
    return -1j*(l_ev_one + l_ev_l + l_ev_r)

psi_t_ref = np.transpose(psi_t)

while end_in_sum < time_length and derivative_sum >= 0:
    derivative_sum -= derivative_evolution[beg_in_sum]
    derivative_sum += derivative_evolution[end_in_sum]
    beg_in_sum += 1
    end_in_sum += 1
end_in_sum -= 1
psi_t = psi_t_ref[:end_in_sum]
psi_next = np.transpose(evolve(psi_t_ref[end_in_sum], times[end_in_sum], times[end_in_sum:], lyapunov_evolution))
psi_t = np.vstack((psi_t, psi_next))
output_dict['1_LC_psi_last'] = psi_t[-1]
psi_t = np.transpose(psi_t)

op_dict_L = obs_vs_time(psi_t, times, op_dict)
Eta_sq_L = op_dict_L['Eta_sq']
Q_L = op_dict_L['Q']
phi_L_1 = phi_evolution[:end_in_sum]
phi_L_2 = np.arcsin(Q_L/Q_norm)[end_in_sum:]
phi_L = np.concatenate((phi_L_1,phi_L_2))
derivative_L = t_h*np.sin(phi_L)*Q_L
H_L = op_dict_L['H']

phi_L = UnivariateSpline(times, phi_L, k=3, s=0)

output_dict['1_break'] = end_in_sum
output_dict['1_LC_phi'] = phi_L(times)
output_dict['1_LC_eta'] = Eta_sq_L
output_dict['1_LC_Q'] = Q_L
output_dict['1_LC_derivative'] = derivative_L
output_dict['1_LC_H'] = H_L

print('2nd 5 are written')

def exp_phi_to_left(t):
    return np.exp(1j * phi_L(t))
def exp_phi_to_right(t):
    return np.exp(-1j * phi_L(t))

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
op_double_evolution = obs_vs_time(psi_t, times, op_double)
Eta_double_evolution = op_double_evolution['Eta_sq']

output_dict['1_LC_eta_DOUBLE'] = Eta_double_evolution

print('2.5th 1 is written')
# 
# start of the curve 2 
# 

omega_p, Phi_0 = 19.1, 0.2

op = tools.Operators(L, s_up, s_dn, a,pbc, t_h, U,
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
Q_norm = ceil(np.max(op.spectrum()[0][0]))

Time = tools.Time_manual_sin_sq(omega_p, t_s, t_f, step_divisor)
times, dt = Time.times()
time_length = len(times)

output_dict['2_times'] = times

psi_t = H.evolve(psi_gs, times[0], times)
op_dict_evolution = obs_vs_time(psi_t, times, op_dict)
phi_evolution = phi(times)
Eta_sq_evolution = op_dict_evolution['Eta_sq']
Q_evolution = op_dict_evolution['Q']
H_evolution = op_dict_evolution['H']
derivative_evolution = t_h*np.sin(phi_evolution)*Q_evolution

output_dict['2_FE_phi'] = phi_evolution
output_dict['2_FE_H'] = H_evolution
output_dict['2_FE_eta'] = Eta_sq_evolution
output_dict['2_FE_Q'] = Q_evolution
output_dict['2_FE_derivative'] = derivative_evolution

print('3rd 5 are written')

beg_in_sum = 0
end_in_sum = int(round(step_divisor,0))
derivative_sum = np.sum(derivative_evolution[beg_in_sum:end_in_sum])

def lyapunov_evolution(t, psi):
    Q_ev = Q.expt_value(psi)
    l_ev_phi = np.arcsin(Q_ev/Q_norm)
    l_ev_one = onsite.dot(psi)
    l_ev_l = np.exp(1j * l_ev_phi)*hop_left.dot(psi)
    l_ev_r = np.exp(-1j * l_ev_phi)*hop_right.dot(psi)
    return -1j*(l_ev_one + l_ev_l + l_ev_r)

psi_t_ref = np.transpose(psi_t)

while end_in_sum < time_length and derivative_sum >= 0:
    derivative_sum -= derivative_evolution[beg_in_sum]
    derivative_sum += derivative_evolution[end_in_sum]
    beg_in_sum += 1
    end_in_sum += 1
end_in_sum -= 1
psi_t = psi_t_ref[:end_in_sum]
psi_next = np.transpose(evolve(psi_t_ref[end_in_sum], times[end_in_sum], times[end_in_sum:], lyapunov_evolution))
psi_t = np.vstack((psi_t, psi_next))
output_dict['2_LC_psi_last'] = psi_t[-1]
psi_t = np.transpose(psi_t)

op_dict_L = obs_vs_time(psi_t, times, op_dict)
Eta_sq_L = op_dict_L['Eta_sq']
Q_L = op_dict_L['Q']
phi_L_1 = phi_evolution[:end_in_sum]
phi_L_2 = np.arcsin(Q_L/Q_norm)[end_in_sum:]
phi_L = np.concatenate((phi_L_1,phi_L_2))
derivative_L = t_h*np.sin(phi_L)*Q_L
H_L = op_dict_L['H']

phi_L = UnivariateSpline(times, phi_L, k=3, s=0)

output_dict['2_break'] = end_in_sum
output_dict['2_LC_phi'] = phi_L(times)
output_dict['2_LC_eta'] = Eta_sq_L
output_dict['2_LC_Q'] = Q_L
output_dict['2_LC_derivative'] = derivative_L
output_dict['2_LC_H'] = H_L

print('4th 5 are written')

def exp_phi_to_left(t):
    return np.exp(1j * phi_L(t))
def exp_phi_to_right(t):
    return np.exp(-1j * phi_L(t))

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
op_double_evolution = obs_vs_time(psi_t, times, op_double)
Eta_double_evolution = op_double_evolution['Eta_sq']

output_dict['2_LC_eta_DOUBLE'] = Eta_double_evolution

print('4.5th 1 is written')

np.save(f'UE_vs_LC_dict_L{L}_18.0-0.2_19.1-0.2_t~o~18.0.npy', output_dict)