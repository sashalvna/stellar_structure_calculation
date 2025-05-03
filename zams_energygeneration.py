import numpy as np

"""
Functions for the energy generation rate
"""

def energy_gen_ppchain(rho, T, X, Y, Z):
    #Get PP chain energy generation rate given density, T, and composition

    #Calculate things for the epsilon_pp equation
    T9 = T/10**9
    T7 = T/10**7
    g11 = 1 + 3.82*T9 + 1.51*T9**2 + 0.144*T9**3 - 0.0114*T9**4
    f11 = np.exp(5.92e-3 * (rho / T7**3)**0.5) #weak screening
    psi = 1 #just use 1

    #if 1.5 >= T7:
    #   psi = 1
    #if (1.5 < T7) and (2.5 > T7):
    #    psi = 2
    #if 2.5 <= T7:
    #    psi = 1.5

    e_pp = 2.57e4 * psi * f11 * g11 * rho * X**2 * T9**(-2/3) * np.exp(-3.381/T9**(1/3))
    return(e_pp)

def energy_gen_cno(rho, T, X, Y, Z):
    #Get CNO cycle energy generation rate given density, T, and composition

    T9 = T/10**9
    g14_1 = 1 - 2*T9 + 3.41*T9**2 - 2.43*T9**3
    X_CNO = (0.173 + 0.053 + 0.482) * Z #C, N, O mass fractions

    e_cno = 8.24e25 * g14_1 * X_CNO * X * rho * T9**(-2/3) * np.exp(-15.231*T9**(-1/3) - (T9/0.8)**2)
    return(e_cno)

def get_epsilon(rho, T, X, Y, Z) :
    #Get total energy generation rate from the PP chain and CNO rates

    epsilon = energy_gen_ppchain(rho, T, X, Y, Z) + energy_gen_cno(rho, T, X, Y, Z)
    return epsilon