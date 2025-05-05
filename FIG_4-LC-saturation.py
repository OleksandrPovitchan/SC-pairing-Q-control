# Warning: angular frequency and amplitude are introduced by hands! LINES 20-21, 98-99

from supporting import tools
from parameters import L, s_up, s_dn, a, pbc, t_h, U
from parameters import pulse_type, sigma_p, t_p, N_p, t_l, asymt_goal
from quspin.tools.evolution import evolve
from quspin.tools.measurements import obs_vs_time

from tqdm import tqdm
import numpy as np

import os
threads = 8
os.environ['OMP_NUM_THREADS'] = '{}'.format(threads)
os.environ['MKL_NUM_THREADS'] = '{}'.format(threads)
os.environ['NUMEXPR_NUM_THREADS'] = '{}'.format(threads)

output_dict = {}

omega_p, Phi_0 = 18.0, 0.2
#omega_p, Phi_0 = 19.1, 0.2
op = tools.Operators(L, s_up, s_dn, a, pbc, t_h, U, 
    pulse_type, omega_p, Phi_0, sigma_p, t_p, N_p, t_l, 
    asymt_goal)
H = op.Hamiltonian('dynamic')
psi_gs = op.spectrum()[0][1][0]
 
onsite = op.Hamiltonian('onsite')
hop_left = op.Hamiltonian('hop_left')
hop_right = op.Hamiltonian('hop_right')
Eta_sq = op.Eta_sq('full')
Q = op.Q('full')
op_dict = {'Eta_sq': Eta_sq}
from math import ceil
Q_norm = ceil(np.max(op.spectrum()[0][0]))

def lyapunov_evolution(t, psi):
    Q_ev = Q.expt_value(psi)
    l_ev_phi = np.arcsin(Q_ev/Q_norm)
    l_ev_one = onsite.dot(psi)
    l_ev_l = np.exp(1j * l_ev_phi)*hop_left.dot(psi)
    l_ev_r = np.exp(-1j * l_ev_phi)*hop_right.dot(psi)
    return -1j*(l_ev_one + l_ev_l + l_ev_r)

from parameters import t_s, t_f_numb, t_r_numb, step_divisor
Time = tools.Time(t_h, pulse_type, omega_p, Phi_0, N_p, t_l, 
    t_s, t_f_numb, t_r_numb, step_divisor)
times, dt = Time.times()
time_length = len(times)

output_dict['times'] = times
output_dict['time_length'] = time_length

i = 0


eta_max = []
eta_final = []

with tqdm(total = time_length) as pbar:
    while i < time_length:

        psi_A = H.evolve(psi_gs, times[0], times[:i+1])
        psi_A = np.transpose(psi_A)

        psi_B = evolve(psi_A[-1], times[i], times[i:], lyapunov_evolution)
        psi_B = np.transpose(psi_B)

        psi_t = np.concatenate((psi_A[:-1], psi_B))
        psi_t = np.transpose(psi_t)

        op_dict_evolution = obs_vs_time(psi_t, times, op_dict)
        eta_evolution = np.real(op_dict_evolution['Eta_sq'])
        eta_m = np.max(eta_evolution)
        eta_f = eta_evolution[-1]
        eta_max.append(eta_m)
        eta_final.append(eta_f)
        

        print(i)

        if (i+1) % 10 == 0:
            '''output_dict['eta_excited'] = eta_excited'''
            output_dict['eta_max'] = np.array(eta_max)
            output_dict['eta_final'] = np.array(eta_final)
            print(f'{i+1}/{time_length}')

            np.save(f'Nonloc_L{L}_18.0-0.2.npy', output_dict)
            #np.save(f'Nonloc_L{L}_19.1-0.2.npy', output_dict)

        pbar.update(1)

        i += 1

output_dict['eta_max'] = np.array(eta_max)
output_dict['eta_final'] = np.array(eta_final)

np.save(f'Nonloc_L{L}_18.0-0.2.npy', output_dict)
#np.save(f'Nonloc_L{L}_19.1-0.2.npy', output_dict)