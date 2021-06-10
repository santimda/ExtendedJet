#!/usr/bin/env python

# First compile with make and then create the wrapper as: 
# f2py -c *f90 -m JetCode
# f2py -c controls.f90 constants.f90 system_parameters.f90 globals.f90 initialize.f90 trajectory.f90 evolution.f90 losses.f90 edistribution.f90 synrad.f90 icrad.f90 sscrad.f90 main.f90 -m JetCode

import JetCode
import numpy as np

#JetCode.run_jet(True)

fB = 1.e-2
Lj = 1.e48
E_data = np.array([1e-19, 1e-16, 1e-13, 1.e-10, 1.e-7, 1e-4])
S_data = E_data # only needs to have the same dimension
S_data = JetCode.sed_model(fB, Lj, E_data)
print(S_data)

