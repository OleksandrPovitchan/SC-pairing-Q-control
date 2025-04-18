from parameters import L, s_up, s_dn, a, pbc, t_h, U
from parameters import pulse_type, omega_p, Phi_0, sigma_p, t_p, N_p, t_l

from parameters import omega_points, Phi0_points, start_from_exp_3
from parameters import omega_min_numb, omega_max_numb, omega_step_numb, A0_min, A0_max, A0_step

from parameters import asymt_goal

import numpy as np
from quspin.tools.measurements import obs_vs_time
import os
from supporting import tools
from tqdm import tqdm

from parameters import threads_exp_3 as threads

os.environ['OMP_NUM_THREADS'] = '{}'.format(threads)
os.environ['MKL_NUM_THREADS'] = '{}'.format(threads)
os.environ['NUMEXPR_NUM_THREADS'] = '{}'.format(threads)

input_path = './exp_2/{}-{}-({})_{}-{}-({})'.format(omega_min_numb, omega_max_numb, omega_step_numb, A0_min, A0_max, A0_step)
output_path = './grid-UE/L{}-{}-{}-({})_{}-{}-({})'.format(L, omega_min_numb, omega_max_numb, omega_step_numb, A0_min, A0_max, A0_step)

os.makedirs(output_path, exist_ok=True)

op = tools.Operators(L, s_up, s_dn, a, pbc, t_h, U, pulse_type, 
        omega_p, Phi_0, sigma_p, t_p, N_p, t_l, 
        asymt_goal)

Eta_pm, Eta_mp, Eta_z = op.Eta_sq('components')

op_dict = {'Eta_pm': Eta_pm, 'Eta_mp': Eta_mp, 'Eta_z': Eta_z}

i = start_from_exp_3
progress_bar = tqdm(total=omega_points*Phi0_points, desc="Processing")
while i < omega_points*Phi0_points:

    psi_t = np.transpose(np.genfromtxt(input_path + '/{}_psi_{}-{}-({})_{}-{}-({}).txt'.format(i,omega_min_numb, omega_max_numb, omega_step_numb, A0_min, A0_max, A0_step), dtype=complex))
    times = np.loadtxt(input_path + '/{}_time_{}-{}-({})_{}-{}-({}).txt'.format(i,omega_min_numb, omega_max_numb, omega_step_numb, A0_min, A0_max, A0_step))

    op_dict_evolution = obs_vs_time(psi_t, times, op_dict)
    
    Eta_pm_evolution = op_dict_evolution['Eta_pm']
    Eta_mp_evolution = op_dict_evolution['Eta_mp']
    Eta_z_evolution = op_dict_evolution['Eta_z']

    Eta_sq_evolution = Eta_pm_evolution + Eta_mp_evolution + Eta_z_evolution**2
    np.savetxt(output_path + '/L{}-eta_sq-{}-{}-({})_{}-{}-({})_{}.txt'.format(L ,omega_min_numb, omega_max_numb, omega_step_numb, A0_min, A0_max, A0_step, i), Eta_sq_evolution)

    i += 1
    progress_bar.update(1)

print('the process is completed')

progress_bar.close()