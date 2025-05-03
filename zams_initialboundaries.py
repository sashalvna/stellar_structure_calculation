import numpy as np
from zams_opacity import *
from zams_eos import *
from zams_energygeneration import *
import constants as const
from scipy.integrate import solve_ivp

"""
Functions for setting up the initial conditions
"""


def get_del_rad(L_r, P_r, T_r, M_r, X, Y, Z):
    #Get radiative temperature gradient given luminosity, pressure, temperature, mass, composition, and the opacity table filename
    rho = get_density(P_r, T_r, X)
    opacity = get_opacity('opacity_x%sy%sz%s.txt'%(X, Y, Z), T_r, rho)
    del_rad =  (3/(16*np.pi)) * P_r * opacity * L_r / (const.a * const.c * const.G * M_r * T_r**4)
    return del_rad

def get_del_ad(P_r, T_r, X):
    #Get adiabatic temperature gradient given pressure, temperature, and composition 
    rho = get_density(P_r, T_r, X)
    beta = get_beta(P_r, T_r, X, rho)
    del_ad = 2*(4 - 3*beta) / (32 - 24*beta - 3*beta**2)
    return del_ad

def get_actual_del(L, P, T, M, X, Y, Z): 
    #Get actual temperature gradient given luminosity, pressure, temperature, mass, composition
    #Calculate and compare radiative and adiabatic gradients
    del_rad =  get_del_rad(L, P, T, M, X, Y, Z)
    del_ad = get_del_ad(P, T, X)

    if del_rad <= del_ad: #if radiative gradient is smaller than adiabatic, radiative
        actual_del = del_rad
    else: #if radiative gradient is smaller than adiabatic, convective
        actual_del =  del_ad
    
    return(actual_del)

def load1(M, M_r, P_c, T_c, X, Y, Z): 
    #Inner boundary conditions (very close to center but with some small enclosed mass)
    #Inputs: total mass, enclosed mass, central pressure (guess), central temperature (guess), composition

    #Get density to calculate radius
    rho_c = get_density(P_c, T_c, X)
    r = (3 / (4*np.pi * rho_c))**(1/3) * M_r

    #Get energy generation rate to get luminosity
    epsilon_c = get_epsilon(rho_c, T_c, X, Y, Z)
    l = epsilon_c * M_r

    #Get pressure
    P = P_c - (3 * const.G / (8*np.pi)) * (4*np.pi/3 * rho_c)**(4/3) * M_r**(2/3) 

    #Get temperature depending on if core is convective or radiative
    if M/const.Ms < 1.5: #radiative core
        T = (T_c**4 - (1/(2*const.a*const.c)) * (3/(4*np.pi))**(2/3) * get_opacity('opacity_x%sy%sz%s.txt'%(X, Y, Z), T_c, rho_c) * epsilon_c * rho_c**(4/3) * M_r**(2/3))**(1/4)
    elif M/const.Ms >= 1.5: #convective core
        T = np.log(T_c) - (np.pi/6)**(1/3) * const.G * get_del_ad(P_c, T_c, X) * rho_c**(4/3)/P_c * M_r**(2/3) 

    return [l, P, r, T]

def load2(M, R, L, X, Y, Z):
    #Outer boundary conditions
    #Inputs: total mass, total radius, luminosity (guess), composition

    #Get the temperature from effective temperature using the Eddington approx.
    T_eff = (L / (4 * np.pi * R**2 * const.sb))**(1/4)
    def get_T(tau):
        T = T_eff * (3/4*(tau + 2/3))**(1/4)
        return T
    
    def get_dP(tau, P):
        T = get_T(tau)
        rho = get_density(P[0], T, X)
        opacity = get_opacity('opacity_x%sy%sz%s.txt'%(X, Y, Z), T, rho)
        dP = const.G * M / (R**2) * 2/3 / opacity
        return dP
    
    #Integrate dP to get P
    def get_P(tau):
        sol = solve_ivp(get_dP, [0, tau], y0=[10])
        P = sol.y[0][-1]
        return P

    tau = 2/3 #photosphere optical depth
    T = get_T(tau)
    P = get_P(tau)

    return [L, P, R, T]