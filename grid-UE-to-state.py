from parameters import omega_points, Phi0_points, start_from_exp_2
from parameters import omega_min_numb, omega_max_numb, omega_step_numb, A0_min, A0_max, A0_step

from parameters import L, s_up, s_dn, a, pbc, t_h, U
from parameters import pulse_type, sigma_p, t_p, N_p, t_l
from parameters import t_s, t_f_numb, t_r_numb, step_divisor

from parameters import asymt_goal

import numpy as np
import os
from tqdm import tqdm
from supporting import tools

from parameters import threads_exp_2 as threads

os.environ['OMP_NUM_THREADS'] = '{}'.format(threads)
os.environ['MKL_NUM_THREADS'] = '{}'.format(threads)
os.environ['NUMEXPR_NUM_THREADS'] = '{}'.format(threads)

input_path = './grid/{}-{}-({})_{}-{}-({})'.format(omega_min_numb, omega_max_numb, omega_step_numb, A0_min, A0_max, A0_step)
output_path = './exp_2/{}-{}-({})_{}-{}-({})'.format(omega_min_numb, omega_max_numb, omega_step_numb, A0_min, A0_max, A0_step)

os.makedirs(output_path, exist_ok=True)

i = start_from_exp_2
progress_bar = tqdm(total=omega_points*Phi0_points, desc="Processing")
while i < omega_points*Phi0_points:

    omega_p = np.loadtxt(input_path + '/{}-{}-({})_{}-{}-({})_{}.txt'.format(omega_min_numb, omega_max_numb, omega_step_numb, A0_min, A0_max, A0_step, i))[0]
    Phi_0 = np.loadtxt(input_path + '/{}-{}-({})_{}-{}-({})_{}.txt'.format(omega_min_numb, omega_max_numb, omega_step_numb, A0_min, A0_max, A0_step, i))[1]

    op = tools.Operators(L, s_up, s_dn, a, pbc, t_h, U, pulse_type, 
        omega_p, Phi_0, sigma_p, t_p, N_p, t_l, 
        asymt_goal)

    H = op.Hamiltonian('dynamic')

    time = tools.Time(t_h, pulse_type, omega_p, Phi_0, N_p, t_l,
        t_s, t_f_numb, t_r_numb, step_divisor)

    times, dt = time.times()                                        # --->
    np.savetxt(output_path + '/{}_time_{}-{}-({})_{}-{}-({}).txt'.format(i,omega_min_numb, omega_max_numb, omega_step_numb, A0_min, A0_max, A0_step), times)

    psi_gs = op.spectrum()[0][1][0]

    psi_t = np.transpose(H.evolve(psi_gs, times[0], times))         # ---> 
    np.savetxt(output_path + '/{}_psi_{}-{}-({})_{}-{}-({}).txt'.format(i,omega_min_numb, omega_max_numb, omega_step_numb, A0_min, A0_max, A0_step), psi_t)
    
    i += 1
    progress_bar.update(1)

progress_bar.close()