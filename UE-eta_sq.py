# (19.1 & 0.2)
# (19.1 & 0.4)
# (18.0 & 0.2)

# time the same as for 18.0 (max) LINE 34
# [as long as pyplot can put several t-axis in the same plot, these changes and remarks are not essential]
# phi function is invariant under the time array transformations

from parameters import L, s_up, s_dn, a, pbc, t_h, U
from parameters import pulse_type, omega_p, Phi_0, sigma_p, t_p, N_p, t_l
from parameters import t_s, t_f_numb, t_r_numb, step_divisor
from parameters import asymt_goal

import numpy as np
from quspin.tools.measurements import obs_vs_time
from supporting import tools
import os

threads = 8

os.environ['OMP_NUM_THREADS'] = '{}'.format(threads)
os.environ['MKL_NUM_THREADS'] = '{}'.format(threads)
os.environ['NUMEXPR_NUM_THREADS'] = '{}'.format(threads)

op = tools.Operators(L, s_up, s_dn, a, pbc, t_h, U, pulse_type, 
        omega_p, Phi_0, sigma_p, t_p, N_p, t_l, 
        asymt_goal)

H = op.Hamiltonian('dynamic')
phi = op.phi
Eta_pm, Eta_mp, Eta_z = op.Eta_sq('components')
op_dict = {'Eta_pm': Eta_pm, 'Eta_mp': Eta_mp, 'Eta_z': Eta_z, 'H_dyn': H}

time_0 = tools.Time(t_h, pulse_type, 18.0, Phi_0, N_p, t_l,
        t_s, t_f_numb, t_r_numb, step_divisor)

t_f = time_0.t_f
t_s = time_0.t_s

time = tools.Time_manual_sin_sq(omega_p, t_s, t_f, step_divisor)

times, dt = time.times()

psi_gs = op.spectrum()[0][1][0]
psi_t = H.evolve(psi_gs, times[0], times) # ready for obs_vs_time

op_dict_evolution = obs_vs_time(psi_t, times, op_dict)

Eta_pm_evolution = op_dict_evolution['Eta_pm']
Eta_mp_evolution = op_dict_evolution['Eta_mp']
Eta_z_evolution = op_dict_evolution['Eta_z']
Energy_evolution = op_dict_evolution['H_dyn']

Eta_sq_evolution = Eta_pm_evolution + Eta_mp_evolution + Eta_z_evolution**2

exp_dict = {
    'Eta_sq_evolution': Eta_sq_evolution,
    'Energy_evolution': Energy_evolution,
    'psi_last': np.transpose(psi_t)[-1],
    'times': times,
    'dt': dt,
    'phi': phi(times)
}

np.save(f'UE_dict_L{L}_o{omega_p}_Phi{Phi_0}_t~o~{18.0}.npy', exp_dict) 