from parameters import omega_min, omega_max, omega_points, Phi0_min, Phi0_max, Phi0_points
from parameters import omega_min_numb, omega_max_numb, omega_step_numb, A0_min, A0_max, A0_step
output = './grid/{}-{}-({})_{}-{}-({})'.format(omega_min_numb, omega_max_numb, omega_step_numb, A0_min, A0_max, A0_step)

import os
import numpy as np
from tqdm import tqdm

from parameters import threads_grid as threads

os.environ['OMP_NUM_THREADS'] = '{}'.format(threads)
os.environ['MKL_NUM_THREADS'] = '{}'.format(threads)
os.environ['NUMEXPR_NUM_THREADS'] = '{}'.format(threads)

os.makedirs(output, exist_ok=True)

omega_array = np.linspace(omega_min, omega_max, omega_points)
Phi0_array = np.linspace(Phi0_min, Phi0_max, Phi0_points)

x, y = np.meshgrid(omega_array, Phi0_array, indexing='xy')
result = np.stack((x, y), axis=-1)
result_linearized = result.reshape(-1, 2)

for i in tqdm(range(omega_points*Phi0_points)):
    np.savetxt(output + '/{}-{}-({})_{}-{}-({})_{}.txt'.format(omega_min_numb, omega_max_numb, omega_step_numb, A0_min, A0_max, A0_step, i), result_linearized[i])

np.savetxt(output + '/{}-{}-({})_{}-{}-({})_omega.txt'.format(omega_min_numb, omega_max_numb, omega_step_numb, A0_min, A0_max, A0_step), omega_array)
np.savetxt(output + '/{}-{}-({})_{}-{}-({})_Phi0.txt'.format(omega_min_numb, omega_max_numb, omega_step_numb, A0_min, A0_max, A0_step), Phi0_array)