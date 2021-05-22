#!/usr/bin/env python

# First compile with make and then create the wrapper as: 
# f2py -c global.f90 dist.f90 radsyn.f90 main.f90 -m jet

import jet
import numpy as np

B = 1.e3

E_data = np.array([1e-15,1e-14,1e-13])

S_data = jet.sed_model(B, E_data)

