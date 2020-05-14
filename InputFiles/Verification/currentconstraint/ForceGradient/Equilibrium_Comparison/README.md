## Comparison of Lfindzero=1 and Lfindzero=2

In this directory we propose to compare the equilibrium found when using the force computed numerically (`>>> ls *001*`) to the one computed analytically (`>>> ls *002*`). Three cases are provided:

* *G1V05L3Fi*: Slab geometry, five plasma volumes, fixed boundary. Execution time: less than 1s 
* *G2V03L3Fi*: Cylindrical geometry, 3 plasma volumes, fixed boundary. Execution time: less than 1s 
* *G3V08L3Fr*: Toroidal geometry, 8 plasma volumes and 1 vacuum region, free boundary. Execution time: ~3000s for .001, ~1120s for .002

Run each case either in sequential or in parallel with MPI. Then compare the output with the routine of your choice (either in Python with `python3 -m py\_spec.ci.test file1.sp.h5 file2.sp.h5` or in Matlab with `compare_spec_outputs('fname1.sp.h5', 'fname2.sp.h5')`). Equilibria should match, up to machine precision.


