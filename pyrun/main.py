# Load spec f90 wrapped
try:
    import spec.spec_f90wrapped as spec
except ImportError as e:
    spec = None
    logger.debug(str(e))

# Load py_spec
try:
    import py_spec
except ImportError as e:
    py_spec = None
    logger.debug(str(e))

# Check that packages are available
if spec is None:
    raise RuntimeError(
        "Using Spec requires spec python wrapper to be installed.")
if py_spec is None:
    raise RuntimeError(
        "Using Spec requires py_spec to be installed.")

from mpi4py import MPI

import numpy as np


# filename
filename = "test.sp"

# Sret up MPI
spec.allglobal.set_mpi_comm(MPI.COMM_WORLD.py2f())

# Prepare input 
spec.inputlist.initialize_inputs()
spec.allglobal.ext = filename[:-3]
spec.allglobal.read_inputlists_from_file()
spec.allglobal.check_inputs()
spec.allglobal.broadcast_inputs()
spec.preset()

# Initialize output
spec.sphdf5.init_outfile()
spec.sphdf5.mirror_input_to_outfile()
spec.allglobal.wrtend()
spec.sphdf5.init_convergence_output()

# Run spec
#spec.spec()
ngdof = spec.allglobal.ngdof
lgdof = spec.allglobal.lgdof
print(lgdof)
mvol = spec.allglobal.mvol
mn = spec.allglobal.mn
irbc = spec.allglobal.irbc
irbs = spec.allglobal.irbs
izbc = spec.allglobal.izbc
izbs = spec.allglobal.izbs
position = np.zeros( (1+ngdof,))
spec.packxi( ngdof, position, mvol, mn, irbc, izbs, irbs, izbc, 'P', False, True)

pressure = spec.inputlist.pressure

for vvol in range(1,mvol):
    vflag = 0
    spec.volume( vvol, vflag ) 
    spec.inputlist.adiabatic[vvol] = pressure[vvol] * spec.allglobal.vvolume[vvol]**spec.inputlist.gamma

ifail = 1

# Expansion of Newton
# spec.newton( ngdof, position )
Ldfjac = ngdof
lr = ngdof * (ngdof+1) / 2
mode = 0
diag = np.ones( (ngdof,) )
factor = spec.inputlist.c05factor
nfcalss = 0
ndcalls = 0



# Unpack 
spec.packxi( ngdof, position, mvol, mn, irbc, izbs, irbs, izbc, 'U', False, True)

spec.allglobal.irbc = irbc
spec.allglobal.irbs = irbs
spec.allglobal.izbc = izbc
spec.allglobal.izbs = izbs

Lcomputederivatives = False
Lcomputeaxis = True
force = np.zeros( (1+ngdof,) )
spec.dforce( ngdof, position, force, Lcomputederivatives, Lcomputeaxis )

if spec.allglobal.forceerr>spec.inputlist.forcetol:
    irevcm = 0
    ihybrid = 1
    fjac = np.zeros((ngdof,ngdof))
    RR = np.zeros( (np.floor(ngdof*(ngdof+1)/2).astype(int),) )
    
    dffdrz = np.zeros( (lgdof,2,lgdof,2,mvol) )
    spec.allglobal.dffdrz = dffdrz
    spec.allglobal.dbbdmp = np.zeros( (lgdof, mvol, 2, 2) )
    if spec.allglobal.localconstraint:
        spec.allglobal.dmupfdx = np.zeros( mvol, 1, 2, lgdof, 2 )
    else:
        spec.allglobal.dmupfdx = np.zeros( mvol, mvol-1, 2, lgdof, 1)

    spec.allglobal.hessian = np.zeros( ngdof, ngdof )
    spec.allglobal.dessian = np.zeros( ngdof, lgdof )
    

    







# End of Newton - some post processing
gradient = np.zeros( (1+ngdof,) )
spec.dforce( ngdof, position, gradient, False, True )

spec.inputlist.tflux = spec.inputlist.tflux * 2 * np.pi / spec.inputlist.phiedge
spec.inputlist.pflux = spec.inputlist.pflux * 2 * np.pi / spec.inputlist.phiedge

spec.ra00aa( 'W' )
spec.allglobal.wrtend()

# Diagnostics
spec.final_diagnostics()

# Write output
spec.sphdf5.write_grid()
spec.allglobal.wrtend()
spec.sphdf5.hdfint()
spec.sphdf5.finish_outfile()

# Terminate
spec.ending()

