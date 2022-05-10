#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is an example of running SPEC via the python_wrapper to scan inputs
and re-converge after each change.
"""

# Below is a hack to get example working from the python_wrapper folder directly.
#
# You don't have to do this if you properly install the python_wrapper
# using `pip install --user .` in the python_wrapper directory.
#
import os.path
python_wrapper_path = os.path.join(os.path.dirname(__file__), '..')
import sys
if not python_wrapper_path in sys.path:
  sys.path.append(python_wrapper_path)
# end of hack...


import numpy as np
import matplotlib.pyplot as plt
from mpi4py import MPI

from spec.core import SPEC

input_file = "G3V02L1Fi.001.sp"

# The test case chosen for this example has `Lfindzero = 2`,
# which tries to save the derivative matrix in between iterations for re-using it.
# The path resolution logic is somewhat broken at the moment for this file.
# Therefore, we need to make a local copy of the input file
# and access it via a path that does not include `../../` etc.
import shutil
input_file_original = os.path.join(python_wrapper_path, "../../InputFiles/TestCases", input_file)
shutil.copy2(input_file_original, input_file)


comm = MPI.COMM_WORLD
rank = comm.rank
spec = SPEC(input_file=input_file, verbose=True, comm=comm)

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
magneticEnergies = []

for i in range(N):
    spec.inputlist.pressure[:Nvol] = np.multiply(scales[i], pressure)
    print("pressure step %d/%d: "%(i+1, N), spec.inputlist.pressure[:Nvol])

    centralPressures.append(spec.inputlist.pressure[0])

    spec.run()

    energy = spec.allglobal.energy
    print("resulting magnetic energy: ",energy)
    magneticEnergies.append(energy)


plt.figure()
plt.plot(centralPressures, magneticEnergies, '.-')
plt.xlabel("central pressure")
plt.ylabel("magnetic energy")
plt.grid(True)
plt.title("pressure scan for case %s"%(input_file,))
plt.tight_layout()

plt.show()
