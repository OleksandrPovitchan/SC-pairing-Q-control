#################################################################
# SYSTEM PARAMETERS
#################################################################

L = 8
s_up = L // 2 + L % 2	    ###
s_dn = L // 2		        ###
a = 1.0 			        ###
pbc = True 
t_h = 1.0 			        ###
U_numb = 20
U = U_numb*t_h 		        ###

#################################################################
# PULSE PARAMETERS
#################################################################

pulse_type = 'sin_sq'	    ###
omega_p_num = 18.68 
omega_p = omega_p_num*t_h	###
A_0 = 0.3 
Phi_0 = a*A_0		        ###
sigma_p = 2/t_h 		    ###
t_p = 10/t_h 			    ###
N_p = 54 
t_l_numb = 5
t_l = t_l_numb/t_h		    ###

#################################################################
# TIME PARAMETERS
#################################################################

t_s = 0.0                   ###
t_f_numb = 30               ###
t_r_numb = 5
step_divisor = 50

#################################################################
# GRID PARAMETERS
#################################################################

omega_min_numb = 0.0
omega_max_numb = 10.0
omega_step_numb = 0.1
omega_points = (omega_max_numb-omega_min_numb)/omega_step_numb + 1  ###
omega_points = int(round(omega_points,0))                           ### --->
omega_min = omega_min_numb*t_h                                      ### --->
omega_max = omega_max_numb*t_h                                      ### --->

A0_min = 0.0
A0_max = 1.5
A0_step = 0.1
A0_points = (A0_max - A0_min)/A0_step + 1                           ###
Phi0_points = int(round(A0_points,0))                               ### --->
Phi0_min = A0_min*a                                                 ### --->
Phi0_max = A0_max*a                                                 ### --->

#################################################################
# EXPERIMENTS
#################################################################

threads_grid = 8

start_from_exp_2 = 0
threads_exp_2 = 8

start_from_exp_2_3 = 0
threads_exp_2_3 = 8

start_from_exp_3 = 0
threads_exp_3 = 8

threads_exp_3_processing = 8

threads_exp_4 = 8

#################################################################
# Lyapunov control
#################################################################

start_from_exp_5_6 = 0
threads_exp_5_6 = 8

threads_exp_6_processing = 8

#################################################################
# Asymptotic control
#################################################################

start_from_exp_7_8 = 0
threads_exp_7_8 = 8

threads_exp_8_processing = 8

asymt_goal = L/2*(L/2+1)
