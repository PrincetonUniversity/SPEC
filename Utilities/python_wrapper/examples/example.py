"""
This is an example of how to run the SPEC python wrapper for optimization.
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

import os.path # again, since hack above does not actually belong to the example...
input_file = os.path.join(python_wrapper_path, "../../InputFiles/TestCases/G3V01L0Fi.001.sp")
# input_file = "G3V01L0Fi.001.sp"

comm = MPI.COMM_WORLD
rank = comm.rank
fun = SPEC(input_file=input_file, verbose=True, comm=comm)

# read input namelist
fun.read()

# change variables freely
fun.inputlist.nppts = 0

# main execution
fun.run()
#fun.run(save_output=True)

# target volume
target = 30.0

# cost function
def func(r):
    fun.allglobal.irbc[9, 1] = r[0]
    fun.lib.volume(1, 0)
    # fun.run()
    return (fun.allglobal.vvolume[0] - target) ** 2

# print function
def callback(r):
    vol = fun.allglobal.vvolume[0]
    if rank == 0:
        print(
            "rbc(1,0)=",
            fun.allglobal.irbc[9, 1],
            "; volume diff: ",
            np.abs(vol - target),
        )
    return vol

# for some reason, this line cannot be put ahead SPEC calls on my MacBook (czhu)
from scipy.optimize import minimize

if rank == 0:
    print("--- Begin to optimize rbc(1,0) to adjust the volume ---")
tol = 1e-8
callback(1.0)
res = minimize(func, [1.0], method="Nelder-Mead", callback=callback, tol=tol)
assert np.isclose(
    fun.allglobal.vvolume[0], target, atol=tol
), "Optimization not successed!"

if rank == 0:
    print(res)
