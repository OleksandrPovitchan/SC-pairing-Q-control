from parameters import L, s_up, s_dn

import numpy as np
dictionary_1 = f'UE_vs_LC_dict_L{L}_18.0-0.2_19.1-0.2_t~o~18.0.npy'
read_dictionary_1 = np.load(dictionary_1, allow_pickle='TRUE').item()
psi_state = read_dictionary_1['2_LC_psi_last']

from parameters import a, pbc, t_h, U
from parameters import pulse_type, omega_p, Phi_0, sigma_p, t_p, N_p, t_l
from parameters import asymt_goal

from supporting import tools
op = tools.Operators(L, s_up, s_dn, a, pbc, t_h, U, pulse_type, 
        omega_p, Phi_0, sigma_p, t_p, N_p, t_l, asymt_goal)

spectrum = op.spectrum()[0]
E_val = spectrum[0]
Psi_val = spectrum[1]
Eta_sq_val = spectrum[2]


weights = [np.square(np.abs(np.dot(psi_state,np.conjugate(Psi_val[i,:])))) for i in range(len(Psi_val))]
precision =  1e-12

black_energy = []
black_eta = []
color_energy = []
color_eta = []
color_weight = []

initial_energy = E_val[0]
initial_eta = Eta_sq_val[0]
initial_weight = weights[0]

i = 0
while i < len(weights):
    if weights[i] < precision:
        weights[i] = 0
        black_energy.append(E_val[i])
        black_eta.append(Eta_sq_val[i])
    else:
        color_energy.append(E_val[i])
        color_eta.append(Eta_sq_val[i])
        color_weight.append(weights[i])
    i += 1

black_eta = np.array(black_eta)
color_eta = np.array(color_eta)

energy_final = read_dictionary_1['2_LC_H'][-1]
print(energy_final)
eta_sq_final = read_dictionary_1['2_LC_eta'][-1]
print(eta_sq_final/L)

exp_dict = {
    'color_energy': color_energy,
    'color_eta': color_eta,
    'color_weight': color_weight,
    'initial_energy': initial_energy,
    'initial_eta': initial_eta,
    'black_energy': black_energy, 
    'black_eta': black_eta,
    'precision': precision,
    'initial_weight': initial_weight,
    'energy_final': energy_final,
    'eta_sq_final': eta_sq_final
}

np.save(f'excited_picture_LC_L{L}.npy', exp_dict) 