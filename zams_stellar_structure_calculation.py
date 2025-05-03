import numpy as np
from zams_opacity import *
from zams_eos import *
from zams_energygeneration import *
from zams_initialboundaries import *
import constants as const
from scipy.integrate import solve_ivp
from scipy import optimize
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import mesa_reader as mr
import pandas as pd

#Plot setttings
plt.rc('font', family='serif')
import matplotlib
plt.rc('font', family='serif', weight='bold')
plt.rc('text', usetex=True)
matplotlib.rcParams['font.weight']= 'bold'
matplotlib.rcParams.update({'font.weight': 'bold'})

def derivs(m, x, X, Y, Z):
    #Four coupled ODEs of stellar structure at a given mass, structure parameters, and composition
    l, P, r, T = x

    #Get actual temperature gradient and density
    actual_del = get_actual_del(l, P, T, m, X, Y, Z)
    rho = get_density(P, T, X)

    dPdm = - const.G * m / (4*np.pi * r**4)
    drdm = 1/(4*np.pi * r**2 * rho)
    dldm = get_epsilon(rho, T, X, Y, Z)
    dTdm = - const.G * m * T / (4*np.pi* r**4 * P) * actual_del

    return [dldm, dPdm, drdm, dTdm]

def integrate(initial_params, M_r, M, X, Y, Z, xf, n_steps):
    #Integrate the four coupled ODEs of stellar structure outward from center and inward from surface
    #Inputs: initial conditions for L, P_c, R, T_c (as an array), the total central very small mass, total mass, composition,
    #central point where inward/outward meet, and number of steps for integration

    #Set up initial parameters and mass per step of integration
    L, P_c, R, T_c = initial_params
    step = M/n_steps

    #Mass ranges for integration
    M_range_out = [M_r, xf] 
    M_range_in = [M, xf]
    M_eval_out = np.linspace(M_r, xf, int(0.5*n_steps)) #which masses to evaluate stellar structure parameters at (going outward)
    M_eval_in = np.linspace(M, xf, int(0.5*n_steps)) #same thing, going inward

    #Integrate outwards from center and inwards from envelope given initial conditions
    outward = solve_ivp(derivs, M_range_out, load1(M, M_r, P_c, T_c, X, Y, Z), args=(X, Y, Z), t_eval=M_eval_out)
    inward = solve_ivp(derivs, M_range_in, load2(M, R, L, X, Y, Z), args=(X, Y, Z), t_eval= M_eval_in)

    return outward, inward


def newt(initial_guess, M_r, M, X, Y, Z, xf, n_steps):
    #Call integrate and check how close the values at the central point xf match up
    #Inputs: initial guess for initial conditions + the inputs to integrate

    #First make sure there are no negative values in the parameters!!! Will break the code because unphysical
    if np.any(initial_guess < 0):
        return([9e99, 9e99, 9e99, 9e99]) #just give the optimizer a high value to make it avoid more negatives

    outward, inward =  integrate(initial_guess, M_r, M, X, Y, Z, xf, n_steps)

    #Check how close the inward and outward calculations are to each other for each input variable
    check = []
    for var in range(4):
        diff = (outward.y[var, -1] - inward.y[var, -1])/(inward.y[var, -1]) #use percent error to compare values (but allow to be negative)
        check.append(diff)

    return(check)

def plotter(outward, inward, ylabels, initial_guess = [], showplot=False):
    #Plot the stellar structure parameters L, P, R, T over mass for the stellar structure calculation results
    #Inputs: outward calculation, inward calculation, y axis labels for plots (array)
    #Optional: initial guess (format as [outward, inward]) to compare to final solution; option to show+save plots or just save them

    fig, axes = plt.subplots(2, 2, figsize = (15, 11))

    for i, ax in enumerate(axes.flat):

        if len(initial_guess) > 0:
            if i==0:
                ax.plot(initial_guess[0].t/const.Ms, initial_guess[0].y[i]/const.Ls, ls='--', lw=2, c='gray', label=r'Initial guess')
                ax.plot(initial_guess[1].t/const.Ms, initial_guess[1].y[i]/const.Ls, ls='--', lw=2, c='gray')
                ax.plot(outward.t/const.Ms, outward.y[i]/const.Ls, lw=3, label=r'Converged, outward')
                ax.plot(inward.t/const.Ms, inward.y[i]/const.Ls, lw=3, label=r'Converged, inward')
            elif i==2:
                ax.plot(initial_guess[0].t/const.Ms, initial_guess[0].y[i]/const.Rs, lw=2, ls='--', c='gray')
                ax.plot(initial_guess[1].t/const.Ms, initial_guess[1].y[i]/const.Rs, lw=2, ls='--', c='gray')
                ax.plot(outward.t/const.Ms, outward.y[i]/const.Rs, lw=3)
                ax.plot(inward.t/const.Ms, inward.y[i]/const.Rs, lw=3)
            else:
                ax.plot(initial_guess[0].t/const.Ms, initial_guess[0].y[i], lw=2, ls='--', c='gray')
                ax.plot(initial_guess[1].t/const.Ms, initial_guess[1].y[i], lw=2, ls='--', c='gray')
                ax.plot(outward.t/const.Ms, outward.y[i], lw=3)
                ax.plot(inward.t/const.Ms, inward.y[i], lw=3)
        else:
            if i==0:
                ax.plot(outward.t/const.Ms, outward.y[i]/const.Ls, lw=3, c='tab:blue')
                ax.plot(inward.t/const.Ms, inward.y[i]/const.Ls, lw=3, c='tab:blue')
            elif i==2:
                ax.plot(outward.t/const.Ms, outward.y[i]/const.Rs, lw=3, c='tab:blue')
                ax.plot(inward.t/const.Ms, inward.y[i]/const.Rs, lw=3, c='tab:blue')
            else:
                ax.plot(outward.t/const.Ms, outward.y[i], lw=3, c='tab:blue')
                ax.plot(inward.t/const.Ms, inward.y[i], lw=3, c='tab:blue')
            
        ax.set_xlim(0, max(inward.t/const.Ms))
        ax.set_ylim(bottom=0)
        ax.set_xlabel(r'$M/M_\odot$', fontsize=20)
        ax.set_ylabel(ylabels[i], fontsize=20)
        ax.tick_params(axis='both', which='major', labelsize=20)
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.05))
        ax.tick_params(length=10, width=2, which='major')
        ax.tick_params(length=5, width=1, which='minor')

    if len(initial_guess) > 0:
        fig.legend(bbox_to_anchor=(0.85, 0.97), fontsize=20, ncol = 3)
        fig.savefig('figures/converged_withinitial.png', bbox_inches='tight', dpi=300)
    else:
        fig.savefig('figures/converged.png', bbox_inches='tight', dpi=300)

    if showplot==True:
        plt.show()

def mesa_plotter(mesapath, mesafile, ylabels, calc=[], showplot=False):
    #Plot the stellar structure parameters P, R, T over mass for the MESA ZAMS structure calculation
    #Inputs: path to MESA LOGS directory, name of MESA data file, y axis labels for plots (array)
    #Optional: stellar structure calculation (format as [outward, inward]) to compare to MESA solution; option to show+save plots or just save them

    fig, axes = plt.subplots(2, 2, figsize = (15, 11))
    mesa_profile = mr.MesaData(mesapath + mesafile) #read MESA data using mesa_reader module

    #For each subplot
    for i, ax in enumerate(axes.flat):
        if i==3: #hide one of the subplots (to only have three - MESA doesn't give luminosity (maybe change this?))
            ax.set_visible(False)
            continue

        #Plot both MESA and converged model
        if len(calc) > 0:
            if i==0:
                ax.plot(calc[0].t/const.Ms, calc[0].y[i+1], lw=3,c='tab:blue', label=r'Converged model')
                ax.plot(calc[1].t/const.Ms, calc[1].y[i+1], lw=3, c='tab:blue')
                ax.plot(mesa_profile.mass, mesa_profile.P, lw=3, c='tab:green', label=r'MESA')
            elif i==1:
                ax.plot(calc[0].t/const.Ms, calc[0].y[i+1]/const.Rs, lw=3,c='tab:blue')
                ax.plot(calc[1].t/const.Ms, calc[1].y[i+1]/const.Rs, lw=3, c='tab:blue')
                ax.plot(mesa_profile.mass, mesa_profile.R, lw=3, c='tab:green')
            elif i==2:
                ax.plot(calc[0].t/const.Ms, calc[0].y[i+1], lw=3,c='tab:blue')
                ax.plot(calc[1].t/const.Ms, calc[1].y[i+1], lw=3, c='tab:blue')
                ax.plot(mesa_profile.mass, mesa_profile.T, lw=3, c='tab:green')
            ax.set_xlim(0, max(calc[1].t/const.Ms))
        #Plot just the MESA result
        else:
            if i==0:
                ax.plot(mesa_profile.mass, mesa_profile.P, lw=3, c='tab:green', label=r'MESA')
            elif i==1:
                ax.plot(mesa_profile.mass, mesa_profile.R, lw=3, c='tab:green')
            elif i==2:
                ax.plot(mesa_profile.mass, mesa_profile.T, lw=3, c='tab:green')
            ax.set_xlim(left=0) 

        ax.set_ylim(bottom=0)
        ax.set_xlabel(r'$M/M_\odot$', fontsize=20)
        ax.set_ylabel(ylabels[i+1], fontsize=20)
        ax.tick_params(axis='both', which='major', labelsize=20)
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.05))
        ax.tick_params(length=10, width=2, which='major')
        ax.tick_params(length=5, width=1, which='minor')

    #Save and (optionally) show plots
    if len(calc) > 0:
        fig.legend(bbox_to_anchor=(0.65, 0.97), fontsize=20, ncol = 2)
        fig.savefig('figures/mesa_withmodel.png', bbox_inches='tight', dpi=300)
    else:
        fig.savefig('figures/mesa.png', bbox_inches='tight', dpi=300)
    if showplot==True:
        plt.show()

def write_tocsv(outward_result, inward_result, X, Y, Z, show_table=False):
    #Save results of stellar structure calculation to a csv
    
    #Make dataframe to store data
    mass_array = np.concatenate((outward_result.t, np.flip(inward_result.t, axis=0)), axis=0) #flip inward to be in correct order (small-> large M)
    full_result = np.concatenate((outward_result.y, np.flip(inward_result.y, axis=1)), axis=1).T #reformat to match masses
    full_result_with_masses = np.column_stack((mass_array.T, full_result))
    table = pd.DataFrame(data = full_result_with_masses, columns=['M', 'L', 'P', 'R', 'T'])

    #Now also need to add density, energy generation rate, opacity, adiabatic temp gradient, actual temp gradient, convective/radiative structure of shell
    table['rho'] = get_density(table['P'], table['T'], X)
    table['epsilon'] = get_epsilon(table['rho'], table['T'], X, Y, Z)
    table['kappa'] = get_opacity('opacity_x%sy%sz%s.txt'%(X, Y, Z), table['T'], table['rho'])
    table['del_ad'] = get_del_ad(table['P'], table['T'], X)
    table['del_rad'] = get_del_rad(table['L'], table['P'], table['T'], table['M'], X, Y, Z)
    table['del_actual'] = [table['del_rad'][i] if table['del_rad'][i] <= table['del_ad'][i] else table['del_ad'][i] for i in range(len(table['del_rad']))]
    table['convective'] = [0 if table['del_rad'][i] <= table['del_ad'][i] else 1 for i in range(len(table['del_rad']))]
    table['radiative'] = [1 if table['del_rad'][i] <= table['del_ad'][i] else 0 for i in range(len(table['del_rad']))]

    if show_table==True:
        print(table)

    table.to_csv('stellar_structure_calculation_result.csv')

if __name__ == "__main__":

    #Set up initial parameters - mass and composition
    M_Msun = 1.0 #mass of star in M_sun
    M = M_Msun * const.Ms #mass of star in g
    X = 0.7 #H mass fraction
    Y = 0.28 #He mass fraction
    Z = 0.02 #metallicity

    #Inner boundary guess
    M_r = 1e-10 #very small nonzero mass
    P_c = 2.5e17 #solar central density
    T_c = 1.5e7 #solar central temperature
    #Outer boundary guess
    R = 1*const.Rs
    L = 1*const.Ls
    #Central point to integrate to
    xf = 0.5*M #just make it halfway
    n_steps = 1000

    initial_guess = np.array([L, P_c, R, T_c])
    print("Initial guess: L=", initial_guess[0]/const.Ls, "L_sun, P_c=", initial_guess[1], "dyn/cm^2, R=", initial_guess[2]/const.Rs, "R_sun, T_c=", initial_guess[3], "K")

    #Set up things for plotting and reading in MESA calculation of same star
    ylabels = [r'$L/L_\odot$', r'$P$ (dyn/cm$^2$)', r'$R/R_\odot$', r'$T$ (K)']
    mesapath = './mesa_zams/LOGS/'
    mesafile = '1_0_Msun_Z_0_02.data'

    #Run optimization to find best initial parameters
    print("Running optimizer to find best initial parameters")
    #sol = optimize.root(newt, initial_guess, args=(M_r, M, X, Y, Z, xf, n_steps), tol=1e-3)
    #print(sol)
    #final_solution = sol.x
    final_solution = np.array([1.72422256e+33, 1.19114608e+17, 9.12182509e+10, 1.25068403e+07])
    print("Final solution: L=", final_solution[0]/const.Ls, "L_sun, P_c=", final_solution[1], "dyn/cm^2, R=", final_solution[2]/const.Rs, "R_sun, T_c=", final_solution[3], "K")

    #Run the stellar structure calculation with the best (final) initial parameters
    integration_result = integrate(final_solution, M_r, M, X, Y, Z, xf, n_steps)
    outward_result = integration_result[0]
    inward_result = integration_result[1]

    #Run the stellar structure calculation with the initial guess for the initial parameters (for plotting)
    initial_guess_result = integrate(initial_guess, M_r, M, X, Y, Z, xf, n_steps)
    outward_guess_result = initial_guess_result[0]
    inward_guess_result = initial_guess_result[1]
    
    #Plot results for the final stellar structure calculation, comparing it with the initial guess and the MESA result for the same star
    plotter(outward_result, inward_result, ylabels, initial_guess = [outward_guess_result, inward_guess_result], showplot=True)
    plotter(outward_result, inward_result, ylabels, showplot=True)
    mesa_plotter(mesapath, mesafile, ylabels, calc=[outward_result, inward_result], showplot=True)
    mesa_plotter(mesapath, mesafile, ylabels, showplot=True)

    #Write data to a table
    write_tocsv(outward_result, inward_result, X, Y, Z, show_table=True)

    print("Stellar structure calculation is done!")

    



