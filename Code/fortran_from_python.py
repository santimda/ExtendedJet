#!/usr/bin/env python

# First compile with make and then create the wrapper as: 
# f2py -c *f90 -m JetCode
# f2py -c controls.f90 constants.f90 system_parameters.f90 globals.f90 initialize.f90 trajectory.f90 evolution.f90 losses.f90 edistribution.f90 synrad.f90 icrad.f90 sscrad.f90 main.f90 -m JetCode

import JetCode
import numpy as np

#Constants
mJy = 1.0e-26
dist = 16.8*3.0856e24  # Distance to the source = 16.8 Mpc
h = 6.6260755e-27     # erg.seg
eV = 1.602176e-12     # 1 eV = 1.602176e-12 erg
a = 4.0*np.pi*dist**2    
M_BH = 6.5e9             #Black hole mass [Msun] (EHT+2018?)
R_g = 1.477e5*M_BH       #Gravitational radius (M_BH [Msun])

# Define the energy range for the SED
E_min = 1.e-7 * eV #Aprox 20 MHz
E_max = 1.e10 * eV
n_points = 10
E_data = np.geomspace(E_min, E_max, n_points)
z_data = np.geomspace(3*R_g, 1e3*R_g, n_points)

# Correr el codigo completo
JetCode.run_jet(True)

#Syn + SSC:
fB = 1.e-2
Lj = 1.e43
z0 = R_g * 10.0 
a_inj = 1.5
S_data = JetCode.sed_model(fB, Lj, z0, a_inj, E_data, z_data)
print(S_data)


# Solo syn:
fB = 1.e-2
Lj = 1.e43
z0 = R_g * 10.0 
a_inj = 1.5
S_data2 = JetCode.syn_model(fB, Lj, z0, a_inj, E_data)
print(S_data2)