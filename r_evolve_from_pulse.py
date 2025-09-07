import numpy as np
import os

threads = 8
os.environ['OMP_NUM_THREADS'] = '{}'.format(threads)
os.environ['MKL_NUM_THREADS'] = '{}'.format(threads)
os.environ['NUMEXPR_NUM_THREADS'] = '{}'.format(threads)

### | 18.0 | 0.2 | ### | 19.1 | 0.2 | ###

L = 8
omega = 18.0
amplitude = 0.2

A_0 = amplitude; omega_p_num = omega

s_up = L // 2 + L % 2; s_dn = L // 2
a = 1.0; pbc = True; t_h = 1.0; U_numb = 20; U = U_numb*t_h

pulse_type = 'sin_sq'; omega_p = omega_p_num*t_h
Phi_0 = a*A_0; sigma_p = 2/t_h; t_p = 10/t_h; N_p = 54; t_l_numb = 5; t_l = t_l_numb/t_h	

asymt_goal = L/2*(L/2+1)

phi_dict = f'r_customisation/L{L}_[{omega}]_[{amplitude}]_PHI-2.npy'

phi_dict = np.load(phi_dict, allow_pickle = 'true').item()

times = phi_dict['times']
phi_4 = phi_dict['phi_D']
phi_5 = phi_dict['phi_E']
phi_6 = phi_dict['phi_F']
phi_7 = phi_dict['phi_G']
phi_8 = phi_dict['phi_H']

from supporting import tools

op = tools.Operators(L, s_up, s_dn, a, pbc, t_h, U,
    pulse_type, omega_p, Phi_0, sigma_p, t_p, N_p, t_l,
    asymt_goal)

psi_gs = op.spectrum()[0][1][0]
Eta_sq = op.Eta_sq('full')
Q = op.Q('full')
op_dict = {'Eta_sq': Eta_sq, 'Q': Q}

from quspin.tools.measurements import obs_vs_time
from quspin.tools.evolution import evolve
from quspin.operators import hamiltonian

from scipy.interpolate import UnivariateSpline
phi_4_spline = UnivariateSpline(times, phi_4, k=3, s=0)

def exp_phi_to_left(t):
    return np.exp(1j * phi_4_spline(t))
def exp_phi_to_right(t):
    return np.exp(-1j * phi_4_spline(t))

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
phi_spline_evolution = phi_4_spline(times)
derivative_spline_evolution = t_h*np.sin(phi_spline_evolution) * Q_spline_evolution

output_dict = {
    "times": times,
    "phi_4": phi_4,
    "eta_4": Eta_spline_evolution,
    "derivative_4": derivative_spline_evolution
}

import os
os.makedirs("r_customisation", exist_ok=True)
np.save(f'r_customisation/L{L}_[{omega}]_[{amplitude}]_ETA-2.npy', output_dict)

print('done')

   ##### ##### ##### ##### ##### ##### ##### ##### #####
##### ##### ##### ##### ##### ##### ##### ##### ##### #####
   ##### ##### ##### ##### ##### ##### ##### ##### #####

from scipy.interpolate import UnivariateSpline
phi_5_spline = UnivariateSpline(times, phi_5, k=3, s=0)

def exp_phi_to_left(t):
    return np.exp(1j * phi_5_spline(t))
def exp_phi_to_right(t):
    return np.exp(-1j * phi_5_spline(t))

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
phi_spline_evolution = phi_5_spline(times)
derivative_spline_evolution = t_h*np.sin(phi_spline_evolution) * Q_spline_evolution

output_dict["phi_5"] = phi_5
output_dict["eta_5"] = Eta_spline_evolution
output_dict["derivative_5"] = derivative_spline_evolution

import os
os.makedirs("r_customisation", exist_ok=True)
np.save(f'r_customisation/L{L}_[{omega}]_[{amplitude}]_ETA-2.npy', output_dict)

print('done')

   ##### ##### ##### ##### ##### ##### ##### ##### #####
##### ##### ##### ##### ##### ##### ##### ##### ##### #####
   ##### ##### ##### ##### ##### ##### ##### ##### #####

from scipy.interpolate import UnivariateSpline
phi_6_spline = UnivariateSpline(times, phi_6, k=3, s=0)

def exp_phi_to_left(t):
    return np.exp(1j * phi_6_spline(t))
def exp_phi_to_right(t):
    return np.exp(-1j * phi_6_spline(t))

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
phi_spline_evolution = phi_6_spline(times)
derivative_spline_evolution = t_h*np.sin(phi_spline_evolution) * Q_spline_evolution

output_dict["phi_6"] = phi_6
output_dict["eta_6"] = Eta_spline_evolution
output_dict["derivative_6"] = derivative_spline_evolution

import os
os.makedirs("r_customisation", exist_ok=True)
np.save(f'r_customisation/L{L}_[{omega}]_[{amplitude}]_ETA-2.npy', output_dict)

print('done')

   ##### ##### ##### ##### ##### ##### ##### ##### #####
##### ##### ##### ##### ##### ##### ##### ##### ##### #####
   ##### ##### ##### ##### ##### ##### ##### ##### #####

from scipy.interpolate import UnivariateSpline
phi_7_spline = UnivariateSpline(times, phi_7, k=3, s=0)

def exp_phi_to_left(t):
    return np.exp(1j * phi_7_spline(t))
def exp_phi_to_right(t):
    return np.exp(-1j * phi_7_spline(t))

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
phi_spline_evolution = phi_7_spline(times)
derivative_spline_evolution = t_h*np.sin(phi_spline_evolution) * Q_spline_evolution

output_dict["phi_7"] = phi_7
output_dict["eta_7"] = Eta_spline_evolution
output_dict["derivative_7"] = derivative_spline_evolution

import os
os.makedirs("r_customisation", exist_ok=True)
np.save(f'r_customisation/L{L}_[{omega}]_[{amplitude}]_ETA-2.npy', output_dict)

print('done')

   ##### ##### ##### ##### ##### ##### ##### ##### #####
##### ##### ##### ##### ##### ##### ##### ##### ##### #####
   ##### ##### ##### ##### ##### ##### ##### ##### #####

from scipy.interpolate import UnivariateSpline
phi_8_spline = UnivariateSpline(times, phi_8, k=3, s=0)

def exp_phi_to_left(t):
    return np.exp(1j * phi_8_spline(t))
def exp_phi_to_right(t):
    return np.exp(-1j * phi_8_spline(t))

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
phi_spline_evolution = phi_8_spline(times)
derivative_spline_evolution = t_h*np.sin(phi_spline_evolution) * Q_spline_evolution

output_dict["phi_8"] = phi_8
output_dict["eta_8"] = Eta_spline_evolution
output_dict["derivative_8"] = derivative_spline_evolution

import os
os.makedirs("r_customisation", exist_ok=True)
np.save(f'r_customisation/L{L}_[{omega}]_[{amplitude}]_ETA-2.npy', output_dict)

print('done')