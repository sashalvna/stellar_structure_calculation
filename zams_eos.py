import constants as const

"""
Functions for the equation of state
"""

def get_mu(X):
    #Calculate mu given X, assuming fully ionized 
    mu = 4 / (3 + 5*X) #Eqn 1.55
    return mu

def get_density(P, T, X):
    #Calculate density given use EoS for gas pressure and radiation pressure for nondegenerate ideal gas
    rho = (P - (const.a * T**4 / 3)) * get_mu(X) * const.mp / (const.k * T) #Eqn 3.17
    return rho

def get_beta(P, T, X, rho):
    #Calculate beta from gas pressure and total pressure 
    beta = rho * const.k * T / get_mu(X) / const.mp / P # Eqn 3.106
    return beta
