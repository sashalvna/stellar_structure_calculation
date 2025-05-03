import numpy as np
from astropy.io import ascii
from scipy import interpolate

"""
Functions for reading in and interpolating an opacity table (from OPAL - https://cds.unistra.fr/topbase/OpacityTables.html)
"""

def read_opacity_table(filename, data_start=3, header_start=2, show_table=False):
    #Read OPAL opacity table given the name of the file (name includes composition)
    #Options: line number for data and header; option to show table

    data = ascii.read(filename, data_start=data_start, header_start=header_start)
    if show_table==True:
        print(data)
    return data

def interpolate_opacity_table(table, T, R):
    #Interpolate the opacity table given the table, the temperature, and R (R = density/T6^3)

    log_R = np.array([float(n) for n in list(table.columns[1:])]) #convert logR from column names into floats
    log_T = np.array([n for n in table['logT']]) #put log(temperature) values into array

    #Create array of opacities in rows/cols of interest
    opacities_tointerp = np.lib.recfunctions.structured_to_unstructured(table.as_array()) 
    opacities_tointerp = np.delete(opacities_tointerp, 0, axis=1)

    #Interpolate the table to get the opacities at T, R
    f_interp = interpolate.interp2d(log_T, log_R, opacities_tointerp.T, kind='cubic')
    opacity = f_interp(T,R)

    return(opacity[0])

def get_opacity(filename, T, rho):
    #Get opacity given an opacity table filename (includes composition in name), the temperature, and density

    #Get logT and logR from T and rho; account for minimum T and R so that values exist on table
    T6 = 1e-6 * T
    T = np.maximum(10**3.75, T)
    R = np.maximum(1e-8, (rho/T6**3)) #convert from density to R for using OPAL tables
    logT = np.log10(T)
    logR = np.log10(R)

    #Get log opacity from table, convert to opacity
    tab = read_opacity_table(filename)
    log_opacity = interpolate_opacity_table(tab, logT, logR)
    opacity = 10**log_opacity

    return opacity
    