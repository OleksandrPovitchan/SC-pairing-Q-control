from quspin.operators import hamiltonian
from quspin.basis import spinful_fermion_basis_1d
import numpy as np
from math import ceil  

class Operators:

    def __init__(self, L, s_up, s_dn, a, pbc, t_h, U, pulse_type, 
        omega_p, Phi_0, sigma_p, t_p, N_p, t_l, asymt_goal):

        if L == 2: pbc = False
        transition_to_left = [[-t_h,i,(i+1)] for i in range(L-1)]
        transition_to_right = [[t_h,i,(i+1)] for i in range(L-1)]
        if pbc:
            transition_to_left.append([-t_h,L-1,0])
            transition_to_right.append([t_h,L-1,0])

        self.eta_sq_max = L/2*(L/2+1)

        self.L = L
        self.Nf = (s_up, s_dn)
        self.a = a
        self.pbc = pbc
        self.t_h = t_h
        self.U = U
        self.asymt_goal = asymt_goal
        self.no_checks = dict(check_pcon=False,check_symm=False,check_herm=False)
        self.basis = spinful_fermion_basis_1d(L = self.L, Nf = self.Nf)

        self.interaction = [[U,i,i] for i in range(L)]
        self.transition_to_left = transition_to_left
        self.transition_to_right = transition_to_right

        self.omega_p = omega_p
        self.Phi_0 = Phi_0

        if omega_p == 0 or Phi_0 == 0: self.pulse_type = 'none'
        else: self.pulse_type = pulse_type

        pulse_type = self.pulse_type

        if pulse_type == 'gaussian':
            self.sigma_p = sigma_p
            self.t_p = t_p
        elif pulse_type == 'sin_sq' or 'sin_sq_envelope':
            self.N_p = N_p
            self.t_l = t_l
    
    def phi(self, t):
        omega_p = self.omega_p
        Phi_0 = self.Phi_0
        if self.pulse_type == 'gaussian':
            sigma_p = self.sigma_p
            t_p = self.t_p
            phi_1 = np.exp( -(t-t_p)**2 / (2*sigma_p**2) )
            phi_2 = np.cos( omega_p*(t-t_p) )
            return Phi_0 * phi_1 * phi_2
        elif self.pulse_type == 'sin_sq':
            N_p = self.N_p
            t_l = self.t_l
            T = (t - t_l)
            phi_1 = np.sin( omega_p * T )
            width = 2*np.pi/(omega_p)*N_p
            phi_2 = (np.sin( omega_p * T / (2 * N_p) ))**2
            heaviside_l = np.heaviside(T, 0.5)
            heaviside_r = np.heaviside(T - width, 0.5)
            phi_3 = heaviside_l - heaviside_r
            return Phi_0 * phi_1 * phi_2 * phi_3
        elif self.pulse_type == 'sin_sq_envelope':
            N_p = self.N_p
            t_l = self.t_l
            T = (t - t_l)
            width = 2*np.pi/(omega_p)*N_p
            phi_2 = (np.sin( omega_p * T / (2 * N_p) ))**2
            heaviside_l = np.heaviside(T, 0.5)
            heaviside_r = np.heaviside(T - width, 0.5)
            phi_3 = heaviside_l - heaviside_r
            return Phi_0 * phi_2 * phi_3
        elif self.pulse_type == 'none':
            return 0*t

    def exp_phi_to_left(self, t):
        return np.exp(1j * self.phi(t))
    def exp_phi_to_right(self, t):
        return np.exp(-1j * self.phi(t))

    def Hamiltonian(self, type):

        interaction = self.interaction
        transition_to_left = self.transition_to_left
        transition_to_right = self.transition_to_right

        basis = self.basis
        no_checks = self.no_checks

        pulse_type = self.pulse_type

        if pulse_type == 'none': type = 'static'

        if type == 'static':
            static = [
            ['n|n', interaction],
            ['+-|', transition_to_left],
            ['-+|', transition_to_right],
            ['|+-', transition_to_left],
            ['|-+', transition_to_right]]
            H_static = hamiltonian(static, [], basis = basis, **no_checks)
            return H_static
        elif type == 'dynamic':
            static = [['n|n', interaction]]
            exp_phi_to_left = self.exp_phi_to_left
            exp_phi_to_right = self.exp_phi_to_right
            dynamic = [
            ['+-|', transition_to_left, exp_phi_to_left, []],
            ['-+|', transition_to_right, exp_phi_to_right, []],
            ['|+-', transition_to_left, exp_phi_to_left, []],
            ['|-+', transition_to_right, exp_phi_to_right, []]]
            H_dynamic = hamiltonian(static, dynamic, basis = basis, **no_checks)
            return H_dynamic
        elif type == 'onsite':
            onsite = hamiltonian([['n|n', interaction]], [], basis = basis, **no_checks)
            return onsite
        elif type == 'hop_left':
            hop_left = hamiltonian([['+-|', transition_to_left], ['|+-', transition_to_left]], [], basis = basis, **no_checks)
            return hop_left
        elif type == 'hop_right':
            hop_right = hamiltonian([['-+|', transition_to_right], ['|-+', transition_to_right]], [], basis = basis, **no_checks)
            return hop_right

    def Eta_sq(self, type):
        L = self.L
        basis = self.basis
        no_checks = self.no_checks
        eta_pm=[[0.5*(-1)**(i+k),i,k,i,k] for i in range(L) for k in range(L)]
        eta_mp=[[0.5*(-1)**(i+k),k,i,k,i] for k in range(L) for i in range(L)]
        if type == 'components':
            static_eta_pm=[['+-|+-', eta_pm]]
            static_eta_mp=[['-+|-+', eta_mp]]
            eta_z_1=[[0.5,i] for i in range(L)]
            eta_z_2=[[0.5,i] for i in range(L)]
            eta_z_3=[[-0.5,i] for i in range(L)]
            static_eta_z=[
                ['n|', eta_z_1],
                ['|n', eta_z_2],
                ['I|', eta_z_3]]
            Eta_pm = hamiltonian(static_eta_pm,[],basis=basis, **no_checks)
            Eta_mp = hamiltonian(static_eta_mp,[],basis=basis, **no_checks)
            Eta_z = hamiltonian(static_eta_z,[],basis=basis, **no_checks)
            return Eta_pm, Eta_mp, Eta_z
        elif type == 'full':
            z_sq_1=[[0.25,i,k] for i in range(L) for k in range(L)]
            z_sq_2=[[0.25,i,k] for i in range(L) for k in range(L)]
            z_sq_3=[[-0.25,i,k] for i in range(L) for k in range(L)]
            z_sq_4=[[0.25,k,i] for i in range(L) for k in range(L)]
            z_sq_5=[[0.25,k,i] for i in range(L) for k in range(L)]
            z_sq_6=[[-0.25,k,i] for i in range(L) for k in range(L)]
            z_sq_7=[[-0.25,i,k] for i in range(L) for k in range(L)]
            z_sq_8=[[-0.25,i,k] for i in range(L) for k in range(L)]
            z_sq_9=[[0.25,i,k] for i in range(L) for k in range(L)]
            static_eta_sq=[
            ['+-|+-', eta_pm],
            ['-+|-+', eta_mp],
            ['nn|',z_sq_1],
            ['n|n',z_sq_2],
            ['n|I',z_sq_3],
            ['n|n',z_sq_4],
            ['|nn',z_sq_5],
            ['I|n',z_sq_6],
            ['In|',z_sq_7],
            ['I|n',z_sq_8],
            ['I|I',z_sq_9]]
            Eta_sq_full = hamiltonian(static_eta_sq,[],basis=basis, **no_checks)
            return Eta_sq_full
    
    def Q(self, type):
        pbc = self.pbc
        L = self.L
        basis = self.basis
        no_checks = self.no_checks
        if type == 'half':
            if pbc:
                q_1=[[(-1.0)**(i+k+1),i,k,(i+1)%L,k] for i in range(L) for k in range(L)]
                q_2=[[(-1.0)**(i+k+1),(i+1)%L,k,i,k] for i in range(L) for k in range(L)]
                q_3=[[(-1.0)**(i+k+1),k,i,k,(i+1)%L] for i in range(L) for k in range(L)]
                q_4=[[(-1.0)**(i+k+1),k,(i+1)%L,k,i] for i in range(L) for k in range(L)]
            else:
                q_1=[[(-1.0)**(i+k+1),i,k,(i+1),k] for i in range(L-1) for k in range(L)]
                q_2=[[(-1.0)**(i+k+1),(i+1),k,i,k] for i in range(L-1) for k in range(L)]
                q_3=[[(-1.0)**(i+k+1),k,i,k,(i+1)] for i in range(L-1) for k in range(L)]
                q_4=[[(-1.0)**(i+k+1),k,(i+1),k,i] for i in range(L-1) for k in range(L)]
            static_q=[
                ['+-|+-',q_1],
                ['+-|+-',q_2],
                ['-+|-+',q_3],
                ['-+|-+',q_4]
                ]
            Q_half = hamiltonian(static_q,[],basis=basis, **no_checks)
            return Q_half
        elif type == 'full':
            if pbc:
                q_1=[[(-1.0)**(i+k+1),i,k,(i+1)%L,k] for i in range(L) for k in range(L)]
                q_2=[[(-1.0)**(i+k+1),(i+1)%L,k,i,k] for i in range(L) for k in range(L)]
                q_3=[[(-1.0)**(i+k+1),k,i,k,(i+1)%L] for i in range(L) for k in range(L)]
                q_4=[[(-1.0)**(i+k+1),k,(i+1)%L,k,i] for i in range(L) for k in range(L)]
                q_c_1=[[(-1.0)**(i+j+1),i,j,i,(j+1)%L] for i in range(L) for j in range(L)]
                q_c_2=[[(-1.0)**(i+j+1),i,(j+1)%L,i,j] for i in range(L) for j in range(L)]
                q_c_3=[[(-1.0)**(i+j+1),j,i,(j+1)%L,i] for i in range(L) for j in range(L)]
                q_c_4=[[(-1.0)**(i+j+1),(j+1)%L,i,j,i] for i in range(L) for j in range(L)]
            else:
                q_1=[[(-1.0)**(i+k+1),i,k,(i+1),k] for i in range(L-1) for k in range(L)]
                q_2=[[(-1.0)**(i+k+1),(i+1),k,i,k] for i in range(L-1) for k in range(L)]
                q_3=[[(-1.0)**(i+k+1),k,i,k,(i+1)] for i in range(L-1) for k in range(L)]
                q_4=[[(-1.0)**(i+k+1),k,(i+1),k,i] for i in range(L-1) for k in range(L)]
                q_c_1=[[(-1.0)**(i+j+1),i,j,i,(j+1)] for i in range(L) for j in range(L-1)]
                q_c_2=[[(-1.0)**(i+j+1),i,(j+1),i,j] for i in range(L) for j in range(L-1)]
                q_c_3=[[(-1.0)**(i+j+1),j,i,(j+1),i] for i in range(L) for j in range(L-1)]
                q_c_4=[[(-1.0)**(i+j+1),(j+1),i,j,i] for i in range(L) for j in range(L-1)]
            static_r=[
                ['+-|+-',q_1],
                ['+-|+-',q_2],
                ['-+|-+',q_3],
                ['-+|-+',q_4],
                ['+-|+-',q_c_1],
                ['+-|+-',q_c_2],
                ['-+|-+',q_c_3],
                ['-+|-+',q_c_4]
                ]
            Q_full = hamiltonian(static_r,[],basis=basis, **no_checks)
            return Q_full
        
    def J_half(self, type):
        basis = self.basis
        no_checks = self.no_checks
        transition_to_left = self.transition_to_left
        pulse_type = self.pulse_type

        if pulse_type == 'none': type = 'static'

        if type == 'static':
            j_half_static = [
            ['+-|', transition_to_left],
            ['|+-', transition_to_left]]
            J_half_static = hamiltonian(j_half_static, [], basis = basis, **no_checks)
            return J_half_static
        elif type == 'dynamic':
            exp_phi_to_left = self.exp_phi_to_left
            j_half_dynamic = [
            ['+-|', transition_to_left, exp_phi_to_left, []],
            ['|+-', transition_to_left, exp_phi_to_left, []]]
            J_half_dynamic = hamiltonian([], j_half_dynamic, basis = basis, **no_checks)
            return J_half_dynamic
        
    def S_sq(self, type):
        L = self.L
        basis = self.basis
        no_checks = self.no_checks
        s_pm = [[0.5,i,j,i,j] for i in range(L) for j in range(L)]
        s_mp = [[0.5,i,j,i,j] for i in range(L) for j in range(L)]
        if type == 'components':
            static_s_pm = [['+-|-+', s_pm]]
            static_s_mp = [['-+|+-', s_mp]]
            s_z_1=[[0.5,i] for i in range(L)]
            s_z_2=[[-0.5,i] for i in range(L)]
            static_s_z=[
                ['n|', s_z_1], 
                ['|n', s_z_2]]
            S_pm = hamiltonian(static_s_pm,[],basis=basis, **no_checks)
            S_mp = hamiltonian(static_s_mp,[],basis=basis, **no_checks)
            S_z = hamiltonian(static_s_z,[],basis=basis, **no_checks)
            return S_pm, S_mp, S_z
        elif type == 'full':
            s_z_sq_1 = [[0.25,i,i,j,j] for i in range(L) for j in range(L)]
            s_z_sq_2 = [[-0.25,i,i,j,j] for i in range(L) for j in range(L)]
            s_z_sq_3 = [[-0.25,j,j,i,i] for i in range(L) for j in range(L)]
            s_z_sq_4 = [[0.25,i,i,j,j] for i in range(L) for j in range(L)]
            static_S_sq = [
            ['+-|-+', s_pm],
            ['-+|+-', s_mp],
            ['+-+-|', s_z_sq_1],
            ['+-|+-', s_z_sq_2],
            ['+-|+-', s_z_sq_3],
            ['|+-+-', s_z_sq_4]]
            S_sq_full = hamiltonian(static_S_sq,[],basis=basis, **no_checks)
            return S_sq_full
        
    def S_c(self, type):
        pbc = self.pbc
        L = self.L
        basis = self.basis
        no_checks = self.no_checks
        if type == 'full':
            s_pm_c = [[0.5,i,i+1,i,i+1] for i in range(L-1)]
            s_mp_c = [[0.5,i,i+1,i,i+1] for i in range(L-1)]
            s_z_sq_c_1 = [[0.25,i,i,i+1,i+1] for i in range(L-1)]
            s_z_sq_c_2 = [[-0.25,i,i,i+1,i+1] for i in range(L-1)]
            s_z_sq_c_3 = [[-0.25,i+1,i+1,i,i] for i in range(L-1)]
            s_z_sq_c_4 = [[0.25,i,i,i+1,i+1] for i in range(L-1)]
            if pbc:
                s_pm_c.append([0.5,L-1,0,L-1,0])
                s_mp_c.append([0.5,L-1,0,L-1,0])
                s_z_sq_c_1.append([0.25,L-1,L-1,0,0])
                s_z_sq_c_2.append([-0.25,L-1,L-1,0,0])
                s_z_sq_c_3.append([-0.25,0,0,L-1,L-1])
                s_z_sq_c_4.append([0.25,L-1,L-1,0,0])
            static_S_c = [
            ['+-|-+', s_pm_c],
            ['-+|+-', s_mp_c],
            ['+-+-|', s_z_sq_c_1],
            ['+-|+-', s_z_sq_c_2],
            ['+-|+-', s_z_sq_c_3],
            ['|+-+-', s_z_sq_c_4]
            ]
            S_c_full = hamiltonian(static_S_c,[],basis=basis, **no_checks)
            return S_c_full
    
    def spectrum(self):
        H_static = self.Hamiltonian('static')
        Eta_sq_full = self.Eta_sq('full')
        S_sq_full = self.S_sq('full')
        S_c_full = self.S_c('full')
        j = self.J_half('static')
        eigenenergy = H_static.eigh()[0]
        eigenstate = np.transpose(H_static.eigh()[1])
        dim = len(eigenenergy)
        
        E_val = [H_static.expt_value(eigenstate[i]) for i in range(dim)]
        Psi_val = eigenstate
        Eta_sq_val = [Eta_sq_full.expt_value(eigenstate[i]) for i in range(dim)]
        S_sq_val = [S_sq_full.expt_value(eigenstate[i]) for i in range(dim)]
        S_c_val = [S_c_full.expt_value(eigenstate[i]) for i in range(dim)]
        j_half = [j.expt_value(eigenstate[i]) for i in range(dim)]
        J_val = 1j*self.a*(j_half - np.conjugate(j_half))

        spectrum = [E_val, Psi_val, Eta_sq_val, S_sq_val, S_c_val, J_val]
        return spectrum, dim
    
    # The following three functions here are written as an example #
    # Their import to code is extremely inefficient #

    def linear_evolution(self, t, psi):
        onsite = self.Hamiltonian('onsite')
        phi = self.phi
        hop_left = self.Hamiltonian('hop_left')
        hop_right = self.Hamiltonian('hop_right')

        ev_1 = onsite.dot(psi)
        ev_2 = np.exp(1j * phi(t))*hop_left.dot(psi)
        ev_3 = np.exp(-1j * phi(t))*hop_right.dot(psi)
        return - 1j * (ev_1 + ev_2 + ev_3)
    
    def lyapunov_evolution(self, t, psi):
        onsite = self.Hamiltonian('onsite')
        hop_left = self.Hamiltonian('hop_left')
        hop_right = self.Hamiltonian('hop_right')

        Q = self.Q('full')
        Q_norm = ceil(np.max(self.spectrum()[0][0]))

        Q_ev = Q.expt_value(psi)
        l_ev_phi = np.arcsin(Q_ev/Q_norm)
        l_ev_one = onsite.dot(psi)
        l_ev_l = np.exp(1j * l_ev_phi)*hop_left.dot(psi)
        l_ev_r = np.exp(-1j * l_ev_phi)*hop_right.dot(psi)
        return -1j*(l_ev_one + l_ev_l + l_ev_r)
    
    def asymptotic_evolution(self, t, psi):
        onsite = self.Hamiltonian('onsite')
        hop_left = self.Hamiltonian('hop_left')
        hop_right = self.Hamiltonian('hop_right')

        Q = self.Q('full')
        Q_norm = ceil(np.max(self.spectrum()[0][0]))
        Eta_sq = self.Eta_sq('full')
        eta_sq_max = self.eta_sq_max
        eta_sq_goal = self.asymt_goal

        Q_ev = Q.expt_value(psi)
        eta_sq_evolution = Eta_sq.expt_value(psi)
        l_ev_phi = np.arcsin(Q_ev*(eta_sq_evolution-eta_sq_goal)/(eta_sq_max*Q_norm))
        l_ev_one = onsite.dot(psi)
        l_ev_l = np.exp(1j * l_ev_phi)*hop_left.dot(psi)
        l_ev_r = np.exp(-1j * l_ev_phi)*hop_right.dot(psi)
        return 1j*(l_ev_one + l_ev_l + l_ev_r)

############################################################################################

class Time:

    def __init__(self, t_h, pulse_type, omega_p, Phi_0, N_p, t_l,
        t_s, t_f_numb, t_r_numb, step_divisor):

        self.omega_p = omega_p
        self.t_s = t_s
        self.step_divisor = step_divisor

        if omega_p == 0 or Phi_0 == 0: self.pulse_type = 'none'
        else: self.pulse_type = pulse_type

        pulse_type = self.pulse_type

        if pulse_type == 'gaussian':
            self.t_f = t_f_numb/t_h
        elif pulse_type == 'sin_sq':
            width = 2*np.pi/(omega_p)*N_p
            self.t_f = t_l + t_r_numb/t_h + width
        elif pulse_type == 'none':
            self.t_f = max(t_f_numb/t_h, t_l + t_r_numb/t_h)

    def times(self):

        t_s = self.t_s
        t_f = self.t_f        
        pulse_type = self.pulse_type
        step_divisor = self.step_divisor

        if pulse_type == 'none': n_steps = step_divisor
        else:
            omega_p = self.omega_p
            T = 2*np.pi/omega_p
            time_length = t_f - t_s
            step_length = T/step_divisor
            n_steps = ceil(time_length/step_length)

        n_points = n_steps + 1

        times, dt = np.linspace(t_s, t_f, n_points, endpoint=True, retstep=True)
        return times, dt
    
class Time_manual_sin_sq:

    def __init__(self, omega_p, t_s, t_f, step_divisor):
        self.t_s = t_s
        self.t_f = t_f
        self.omega_p = omega_p
        self.t_s = t_s
        self.step_divisor = step_divisor

    def times(self):

        t_s = self.t_s
        t_f = self.t_f
        step_divisor = self.step_divisor
        omega_p = self.omega_p

        T = 2*np.pi/omega_p
        time_length = t_f - t_s
        step_length = T/step_divisor
        n_steps = ceil(time_length/step_length)

        n_points = n_steps + 1

        times, dt = np.linspace(t_s, t_f, n_points, endpoint=True, retstep=True)
        return times, dt
