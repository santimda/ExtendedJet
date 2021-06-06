#!/usr/bin/env python

# First compile with make and then create the wrapper as: 
# f2py -c constants.f90 global.f90 dist.f90 radsyn.f90 radssc.f90 main.f90 -m jet

import jet
import numpy as np

B = 1.e3
U = 1.e5*1.6e-12

E_data = np.array([1e-19, 1e-14, 1.e-10, 2.e-8, 1.e-5, 1.e-3])

S_data = jet.sed_model(B, U, E_data)

print(S_data)
