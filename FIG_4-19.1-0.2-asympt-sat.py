# calculates evolution under unperturbed pulse 19.1-0.2
# calculate evolution under asymptotic control ('derivative activation condition')
# does no double-check simulation (evolution under precalculated perturbed field)
# checks maximal and saturated values depending on the activation time
# exported variables are: 'phi_NC', 'eta_NC', 'eta_AC', 'eta_max', 'eta_final'

# (remember that in asymptotic unitary we take - 1j times RHS,
# i.e., it's not enough to just copy arcsine from the definition
# of the function to get the field 'phi_AC')

import numpy as np
from supporting import tools
from quspin.tools.evolution import evolve
from quspin.tools.measurements import obs_vs_time
output_dict = {}

# for smoother calculations on the cluster

import os

threads = 8
os.environ['OMP_NUM_THREADS'] = '{}'.format(threads)
os.environ['MKL_NUM_THREADS'] = '{}'.format(threads)
os.environ['NUMEXPR_NUM_THREADS'] = '{}'.format(threads)

# settings & time

from parameters import(L, s_up, s_dn, a, pbc, t_h, U,
    pulse_type, sigma_p, t_p, N_p, t_l,
    asymt_goal)
omega_p, Phi_0 = 19.1, 0.2
asymt_goal = L/2*(L/2+1)

from parameters import t_s, t_f_numb, t_r_numb, step_divisor
Time = tools.Time(t_h, pulse_type, omega_p, Phi_0, N_p, t_l,
    t_s, t_f_numb, t_r_numb, step_divisor)

times, dt = Time.times()
time_length = len(times)

output_dict['x'] = times

# import aparata

op = tools.Operators(L, s_up, s_dn, a, pbc, t_h, U,
    pulse_type, omega_p, Phi_0, sigma_p, t_p, N_p, t_l,
    asymt_goal)

H = op.Hamiltonian('dynamic')
Eta_sq = op.Eta_sq('full')
Q = op.Q('full')

phi = op.phi

psi_gs = op.spectrum()[0][1][0]

op_dict = {'Eta_sq': Eta_sq, 'Q': Q}

# asymptotic unitary

onsite = op.Hamiltonian('onsite')
hop_left = op.Hamiltonian('hop_left')
hop_right = op.Hamiltonian('hop_right')
eta_sq_max = L/2*(L/2+1)
from math import ceil
Q_norm = ceil(np.max(op.spectrum()[0][0]))

def asymptotic_evolution(t, psi):
    Q_ev = Q.expt_value(psi)
    eta_sq_evolution = Eta_sq.expt_value(psi)
    l_ev_phi = np.arcsin(Q_ev*(eta_sq_evolution-asymt_goal)/(eta_sq_max*Q_norm))
    l_ev_one = onsite.dot(psi)
    l_ev_l = np.exp(1j * l_ev_phi)*hop_left.dot(psi)
    l_ev_r = np.exp(-1j * l_ev_phi)*hop_right.dot(psi)
    return 1j*(l_ev_one + l_ev_l + l_ev_r)

# algorithm (no control)

psi_t = H.evolve(psi_gs, times[0], times)
psi_t_ref = np.transpose(psi_t)
phi_evolution = phi(times)
op_dict_evolution = obs_vs_time(psi_t, times, op_dict)
Eta_sq_evolution = op_dict_evolution['Eta_sq']
Q_evolution = op_dict_evolution['Q']

output_dict['phi_NC'] = phi(times)
output_dict['eta_NC'] = Eta_sq_evolution
derivative_evolution = t_h*np.sin(phi_evolution)*Q_evolution
np.save(f'{omega_p}-{Phi_0}-asympt-sat-(L{L}).npy', output_dict)

# algorithm (asymptotic with usual condition)

beg_in_sum = 0
end_in_sum = int(round(step_divisor,0))
derivative_sum = np.sum(derivative_evolution[beg_in_sum:end_in_sum])

while end_in_sum < time_length and derivative_sum >= 0:
    derivative_sum -= derivative_evolution[beg_in_sum]
    derivative_sum += derivative_evolution[end_in_sum]
    beg_in_sum += 1
    end_in_sum += 1
end_in_sum -= 1
psi_t = psi_t_ref[:end_in_sum]
psi_next = np.transpose(evolve(psi_t_ref[end_in_sum], times[end_in_sum], times[end_in_sum:], asymptotic_evolution))
psi_t = np.vstack((psi_t, psi_next))
psi_t = np.transpose(psi_t)

op_dict_evolution = obs_vs_time(psi_t, times, op_dict)
Eta_sq_evolution = op_dict_evolution['Eta_sq']

output_dict['eta_AC'] = Eta_sq_evolution
np.save(f'{omega_p}-{Phi_0}-asympt-sat-(L{L}).npy', output_dict)

# algorithm (saturated value)

i = 0

eta_max = []
eta_final = []

from tqdm import tqdm

with tqdm(total = time_length) as pbar:
    while i < time_length:

        psi_A = H.evolve(psi_gs, times[0], times[:i+1])
        psi_A = np.transpose(psi_A)

        psi_B = evolve(psi_A[-1], times[i], times[i:], asymptotic_evolution)
        psi_B = np.transpose(psi_B)

        psi_t = np.concatenate((psi_A[:-1], psi_B))
        psi_t = np.transpose(psi_t)

        op_dict_evolution = obs_vs_time(psi_t, times, op_dict)
        eta_evolution = np.real(op_dict_evolution['Eta_sq'])
        eta_m = np.max(eta_evolution)
        eta_f = eta_evolution[-1]
        eta_max.append(eta_m)
        eta_final.append(eta_f)

        print(f'{i+1}/{time_length}')

        if (i+1) % 10 == 0:
            output_dict['eta_max'] = np.array(eta_max)
            output_dict['eta_final'] = np.array(eta_final)
            np.save(f'{omega_p}-{Phi_0}-asympt-sat-(L{L}).npy', output_dict)

        pbar.update(1)

        i += 1

output_dict['eta_max'] = np.array(eta_max)
output_dict['eta_final'] = np.array(eta_final)

np.save(f'{omega_p}-{Phi_0}-asympt-sat-(L{L}).npy', output_dict)