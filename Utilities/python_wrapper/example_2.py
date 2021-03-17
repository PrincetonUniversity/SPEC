#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 12:03:45 2021

@author: jonathan
"""

import numpy as np
import matplotlib.pyplot as plt
from mpi4py import MPI

# you can also use sys.path.append
# import sys
# sys.path.append("./SPEC/Utilities/python_wrapper/")
from core import SPEC

# need to copy this from the InputFiles/TestCases folder into the python_wrapper folder for now
ext = "G3V02L1Fi.001.sp"

comm = MPI.COMM_WORLD
rank = comm.rank
spec = SPEC(input_file=ext, verbose=True, comm=comm)
# read input namelist
spec.read()
# change variables freely
spec.inputlist.nppts = 0

mu0 = 4.0e-7*np.pi

Nvol=spec.inputlist.nvol
pressure = spec.inputlist.pressure[:Nvol]
print("pressure from input: ", pressure)

N = 5
scales = np.linspace(1.0, 2.0, N)

centralPressures = []
plasmaEnergies = []

for i in range(N):
    spec.inputlist.pressure[:Nvol] = np.multiply(scales[i], pressure)
    print("pressure step %d/%d: "%(i+1, N), spec.inputlist.pressure[:Nvol])
    
    centralPressures.append(spec.inputlist.pressure[0])
    
    spec.run()
    
    energy = spec.allglobal.energy
    print("resulting plasma energy: ",energy)
    plasmaEnergies.append(energy)


plt.figure()
plt.plot(centralPressures, plasmaEnergies)
plt.xlabel("central pressure")
plt.ylabel("plasma energy")
plt.grid(True)
plt.tight_layout()

plt.show()