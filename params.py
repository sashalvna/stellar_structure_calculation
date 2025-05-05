import constants as const

#Initial parameters
M_Msun = 1.0 #mass of star in M_sun
M = M_Msun * const.Ms #mass of star in g (don't change this)
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

#MESA files
mesapath = './mesa_zams/LOGS/'
mesafile = '1_0_Msun_Z_0_02.data'