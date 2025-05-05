# data from the double-pulse simulation


from supporting import tools
from parameters import L, s_up, s_dn, a, pbc, t_h, U
from parameters import pulse_type, sigma_p, t_p, N_p, t_l
from math import ceil
from quspin.tools.evolution import evolve
from quspin.tools.measurements import obs_vs_time
import numpy as np

import os
threads = 8
os.environ['OMP_NUM_THREADS'] = '{}'.format(threads)
os.environ['MKL_NUM_THREADS'] = '{}'.format(threads)
os.environ['NUMEXPR_NUM_THREADS'] = '{}'.format(threads)

asymt_goal = 0
omega_p, Phi_0 = 19.1, 0.2
export_dict = {}

from parameters import t_s, t_f_numb, t_r_numb, step_divisor
Time = tools.Time(t_h, pulse_type, omega_p, Phi_0, N_p, t_l, t_s, t_f_numb, t_r_numb, step_divisor)
times, dt = Time.times()
export_dict['times'] = times
time_length = len(times)

op = tools.Operators(L, s_up, s_dn, a, pbc, t_h, U,
    pulse_type, omega_p, Phi_0, sigma_p, t_p, N_p, t_l,
    asymt_goal)

phi = op.phi
H = op.Hamiltonian('dynamic')
Eta_sq = op.Eta_sq('full')

eta_dict = {'Eta_sq': Eta_sq}

psi_gs = op.spectrum()[0][1][0]
psi_t_before = H.evolve(psi_gs, times[0], times)                            # evolution that excites system by unperturbed pulse
phi_evolution = phi(times)
export_dict['phi_FE'] = phi_evolution
eta_dict_evolution_before = obs_vs_time(psi_t_before, times, eta_dict)
eta_evolution_before = np.real(eta_dict_evolution_before['Eta_sq'])
export_dict['eta_FE_before'] = eta_evolution_before

psi_t_before = np.transpose(psi_t_before)
psi_f = psi_t_before[-1]                                                    # the last state of the uncontrolled evolution

print('pre-evolution: DONE')

Q = op.Q('full')
Q_norm = ceil(np.max(op.spectrum()[0][0]))
onsite = op.Hamiltonian('onsite')
hop_left = op.Hamiltonian('hop_left')
hop_right = op.Hamiltonian('hop_right')
def lyapunov_evolution(t, psi):                                             # defines Lyapunov suppression
    Q_ev = Q.expt_value(psi)
    l_ev_phi = - np.arcsin(Q_ev/Q_norm)
    l_ev_one = onsite.dot(psi)
    l_ev_l = np.exp(1j * l_ev_phi)*hop_left.dot(psi)
    l_ev_r = np.exp(-1j * l_ev_phi)*hop_right.dot(psi)
    return -1j*(l_ev_one + l_ev_l + l_ev_r)

i = time_length - 10

eta_min = []
eta_final = []
from tqdm import tqdm                                                       # saturated and max values depending on the activation time
with tqdm(total = time_length) as pbar:
    while i < time_length:

        psi_A = H.evolve(psi_f, times[0], times[:i+1])
        psi_A = np.transpose(psi_A)

        psi_B = evolve(psi_A[-1], times[i], times[i:], lyapunov_evolution)
        psi_B = np.transpose(psi_B)

        psi_t_steps = np.concatenate((psi_A[:-1], psi_B))
        psi_t_steps = np.transpose(psi_t_steps)

        eta_dict_evolution_steps = obs_vs_time(psi_t_steps, times, eta_dict)
        eta_evolution_steps = np.real(eta_dict_evolution_steps['Eta_sq'])
        eta_m = np.min(eta_evolution_steps)
        eta_f = eta_evolution_steps[-1]
        eta_min.append(eta_m)
        eta_final.append(eta_f)

        if (i+1) % 10 == 0:
            export_dict['eta_min'] = np.array(eta_min)
            export_dict['eta_final'] = np.array(eta_final)
            print(f'{i+1}/{time_length}')
            np.save(f'nonloc_sup_new_L{L}_omega{omega_p}_phi{Phi_0}.npy', export_dict)

        pbar.update(1)

        i += 1

export_dict['eta_min'] = np.array(eta_min)
export_dict['eta_final'] = np.array(eta_final)

np.save(f'nonloc_sup_new_L{L}_omega{omega_p}_phi{Phi_0}.npy', export_dict)
print('SAVED the majority')

dictop = {'Q': Q, 'Eta_sq': Eta_sq}
psi_t_after = H.evolve(psi_f, times[0], times)                              # suppresses by an unperturbed pulse
psi_t_ref = np.transpose(psi_t_after)
dictop_evolution = obs_vs_time(psi_t_after, times, dictop)
Eta_sq_evolution_after = dictop_evolution['Eta_sq']
export_dict['eta_FE_after'] = Eta_sq_evolution_after
Q_evolution_after = dictop_evolution['Q']
derivative_evolution_after = t_h*np.sin(phi_evolution)*Q_evolution_after

print('SAVED FE after')
beg_in_sum = 0
end_in_sum = int(round(step_divisor,0))
derivative_sum = np.sum(derivative_evolution_after[beg_in_sum:end_in_sum])

while end_in_sum < time_length and derivative_sum <= 0:                     # searches positive mean derivative
    derivative_sum -= derivative_evolution_after[beg_in_sum]
    derivative_sum += derivative_evolution_after[end_in_sum]
    beg_in_sum += 1
    end_in_sum += 1
end_in_sum -= 1

psi_t = psi_t_ref[:end_in_sum]
psi_next = np.transpose(evolve(psi_t_ref[end_in_sum], times[end_in_sum], times[end_in_sum:], lyapunov_evolution))
psi_t = np.vstack((psi_t, psi_next))
psi_t = np.transpose(psi_t)                                                 # concatenated evolution
dictop_L = obs_vs_time(psi_t, times, dictop)

Eta_sq_evolution_L = dictop_L['Eta_sq']                                     # eta_sq - derivative activation condition
Q_L = dictop_L['Q']

phi_L_1 = phi_evolution[:end_in_sum]
phi_L_2 = - np.arcsin(Q_L/Q_norm)[end_in_sum:]
phi_LC = np.concatenate((phi_L_1,phi_L_2))

export_dict['eta_LC'] = Eta_sq_evolution_L
export_dict['phi_LC'] = phi_LC

np.save(f'nonloc_sup_new_L{L}_omega{omega_p}_phi{Phi_0}.npy', export_dict)

print('SAVED LC')

ETA_FINAL = np.array(eta_final)
the_min_sup = np.min(ETA_FINAL)
the_right_moment = np.where(ETA_FINAL == the_min_sup)[0][0]                 # searching the most suitable activating moment

export_dict['the_min_sup'] = the_min_sup
export_dict['the_right_moment'] = the_right_moment

psi_t = psi_t_ref[:the_right_moment]
psi_next = np.transpose(evolve(psi_t_ref[the_right_moment], times[the_right_moment], times[the_right_moment:], lyapunov_evolution))

psi_t = np.vstack((psi_t, psi_next))
psi_t = np.transpose(psi_t)
dictop_L_star = obs_vs_time(psi_t, times, dictop)
Eta_sq_evolution_L_star = dictop_L_star['Eta_sq']
Q_L_star = dictop_L_star['Q']
phi_L_1 = phi_evolution[:the_right_moment]
phi_L_2 = - np.arcsin(Q_L/Q_norm)[the_right_moment:]
phi_LC = np.concatenate((phi_L_1,phi_L_2))

export_dict['eta_LC_star'] = Eta_sq_evolution_L_star
export_dict['phi_LC_star'] = phi_LC

np.save(f'nonloc_sup_new_L{L}_omega{omega_p}_phi{Phi_0}.npy', export_dict)

print('SAVED LC_star')