import numpy as np
import os
from tqdm import tqdm

from parameters import omega_points, Phi0_points, L
from parameters import omega_min_numb, omega_max_numb, omega_step_numb, A0_min, A0_max, A0_step

from parameters import threads_exp_6_processing as threads

os.environ['OMP_NUM_THREADS'] = '{}'.format(threads)
os.environ['MKL_NUM_THREADS'] = '{}'.format(threads)
os.environ['NUMEXPR_NUM_THREADS'] = '{}'.format(threads)

input_path = './grid-LQC/L{}-LQC-{}-{}-({})_{}-{}-({})'.format(L, omega_min_numb, omega_max_numb, omega_step_numb, A0_min, A0_max, A0_step)
output_path = './grid-LQC/processing'

os.makedirs(output_path, exist_ok=True)

Eta_sq_evolution_max = []
Eta_sq_evolution_final = []

i = 0

progress_bar = tqdm(total=omega_points*Phi0_points, desc="Processing")

while i < omega_points*Phi0_points:
    data_i = np.genfromtxt(input_path + '/L{}-LQC-eta_sq-{}-{}-({})_{}-{}-({})_{}.txt'.format(L ,omega_min_numb, omega_max_numb, omega_step_numb, A0_min, A0_max, A0_step, i), dtype=complex)
    data_i = np.real(data_i)
    Eta_sq_evolution_max.append(np.max(data_i))
    Eta_sq_evolution_final.append(data_i[-1])
    i += 1
    progress_bar.update(1)

progress_bar.close()

print('the process is completed')

Eta_sq_evolution_max_reshape = np.array(Eta_sq_evolution_max).reshape(Phi0_points,omega_points)
Eta_sq_evolution_final_reshape = np.array(Eta_sq_evolution_final).reshape(Phi0_points,omega_points)

np.savetxt(output_path + '/L{}-LQC-eta_sq_max-{}-{}-({})_{}-{}-({}).txt'.format(L ,omega_min_numb, omega_max_numb, omega_step_numb, A0_min, A0_max, A0_step), Eta_sq_evolution_max_reshape)
np.savetxt(output_path + '/L{}-LQC-eta_sq_final-{}-{}-({})_{}-{}-({}).txt'.format(L ,omega_min_numb, omega_max_numb, omega_step_numb, A0_min, A0_max, A0_step), Eta_sq_evolution_final_reshape)