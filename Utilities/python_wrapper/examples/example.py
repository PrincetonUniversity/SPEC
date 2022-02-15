import numpy as np
import matplotlib.pyplot as plt
from mpi4py import MPI

# you can also use sys.path.append
# import sys
# sys.path.append("./SPEC/Utilities/python_wrapper/")
from spec.core import SPEC

ext = "../../../InputFiles/TestCases/G3V01L0Fi.001.sp"
# ext = "G3V01L0Fi.001.sp"

comm = MPI.COMM_WORLD
rank = comm.rank
fun = SPEC(input_file=ext, verbose=True, comm=comm)
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
