from parameters import omega_points, Phi0_points
from parameters import omega_min_numb, omega_max_numb, omega_step_numb, A0_min, A0_max, A0_step

from parameters import L, s_up, s_dn, a, pbc, t_h, U
from parameters import pulse_type, sigma_p, t_p, N_p, t_l
from parameters import t_s, t_f_numb, t_r_numb, step_divisor

from parameters import asymt_goal

import numpy as np
import os
from tqdm import tqdm
from supporting import tools
from quspin.tools.measurements import obs_vs_time
from quspin.tools.evolution import evolve
from math import ceil

from parameters import start_from_exp_7_8 as start_from
from parameters import threads_exp_7_8 as threads

os.environ['OMP_NUM_THREADS'] = '{}'.format(threads)
os.environ['MKL_NUM_THREADS'] = '{}'.format(threads)
os.environ['NUMEXPR_NUM_THREADS'] = '{}'.format(threads)

input_path = './grid/{}-{}-({})_{}-{}-({})'.format(omega_min_numb, omega_max_numb, omega_step_numb, A0_min, A0_max, A0_step)
output_path = './exp_8/L{}-AQC-{}-{}-({})_{}-{}-({})'.format(L, omega_min_numb, omega_max_numb, omega_step_numb, A0_min, A0_max, A0_step)

os.makedirs(output_path, exist_ok=True)

eta_sq_max = L/2*(L/2+1)

i = start_from
progress_bar = tqdm(total=omega_points*Phi0_points, desc="Processing")
while i < omega_points*Phi0_points:

    omega_p = np.loadtxt(input_path + '/{}-{}-({})_{}-{}-({})_{}.txt'.format(omega_min_numb, omega_max_numb, omega_step_numb, A0_min, A0_max, A0_step, i))[0]
    Phi_0 = np.loadtxt(input_path + '/{}-{}-({})_{}-{}-({})_{}.txt'.format(omega_min_numb, omega_max_numb, omega_step_numb, A0_min, A0_max, A0_step, i))[1]

    Time = tools.Time(t_h, pulse_type, omega_p, Phi_0, N_p, t_l,
        t_s, t_f_numb, t_r_numb, step_divisor)
    times, dt = Time.times() # ready for obs_vs_time
    time_length = len(times)

    op = tools.Operators(L, s_up, s_dn, a, pbc, t_h, U, pulse_type, 
        omega_p, Phi_0, sigma_p, t_p, N_p, t_l, 
        asymt_goal)
    H = op.Hamiltonian('dynamic')
    onsite = op.Hamiltonian('onsite')
    hop_left = op.Hamiltonian('hop_left')
    hop_right = op.Hamiltonian('hop_right')
    Eta_sq = op.Eta_sq('full')
    op_dict = {'Eta_sq': Eta_sq}
    Q = op.Q('full')
    Q_norm = ceil(np.max(op.spectrum()[0][0]))
    q_dict = {'Q': Q}

    phi = op.phi
    
    psi_gs = op.spectrum()[0][1][0]
    psi_t_ref = H.evolve(psi_gs, times[0], times) # ready for obs_vs_time
    psi_t_ref = np.transpose(psi_t_ref)

    # ASYMPTOT # # ASYMPTOT # # ASYMPTOT # # ASYMPTOT # # ASYMPTOT #

    def asymptotic_evolution(t, psi):
        Q_ev = Q.expt_value(psi)
        eta_sq_evolution = Eta_sq.expt_value(psi)
        l_ev_phi = np.arcsin(Q_ev*(eta_sq_evolution-asymt_goal)/(eta_sq_max*Q_norm))
        l_ev_one = onsite.dot(psi)
        l_ev_l = np.exp(1j * l_ev_phi)*hop_left.dot(psi)
        l_ev_r = np.exp(-1j * l_ev_phi)*hop_right.dot(psi)
        return 1j*(l_ev_one + l_ev_l + l_ev_r)
    
    # # #  # # # # # #  # # # # # #  # # # # # #  # # # # # #  # # #

    derivative_ref = t_h*np.sin(phi(times))*obs_vs_time(np.transpose(psi_t_ref), times, q_dict)['Q']
    beg_in_sum = 0
    end_in_sum = int(round(step_divisor,0))

    derivative_sum = np.sum(derivative_ref[beg_in_sum:end_in_sum])
    while end_in_sum < time_length and derivative_sum >= 0:
        derivative_sum -= derivative_ref[beg_in_sum]
        derivative_sum += derivative_ref[end_in_sum]
        beg_in_sum += 1
        end_in_sum += 1
    end_in_sum -= 1 # brings us to the step back when the sum is positive
    psi_t = psi_t_ref[:end_in_sum] # index beg_in_sum-1 must NOT be included

    psi_next = np.transpose(evolve(psi_t_ref[end_in_sum], times[end_in_sum], times[end_in_sum:], asymptotic_evolution))
    psi_t = np.vstack((psi_t, psi_next))
    psi_t = np.transpose(psi_t)
    Eta_sq_asympt = obs_vs_time(psi_t, times, op_dict)['Eta_sq']

    np.savetxt(output_path + '/L{}-AQC-eta_sq-{}-{}-({})_{}-{}-({})_{}.txt'.format(L ,omega_min_numb, omega_max_numb, omega_step_numb, A0_min, A0_max, A0_step, i), Eta_sq_asympt)

    i += 1
    progress_bar.update(1)

print('the process is completed')

progress_bar.close()