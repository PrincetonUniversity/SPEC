
TESTCASES COMPARISON
--------------------

In this directory, the output of SPEC is compared between a local constraint (either L0 or L1) and the current constraint (L3) for some testcases proposed in SPEC Inputfile directory. We proceed as follows:

	1. Compile the master version of SPEC and copy it in this directory with filename `xspec_master`
```bash
git checkout master
make
cp xspec InputFiles/Verification/currentconstraint/TestCases_Comparison/xspec_master
```

	2. Compile the current constraint branch of SPEC and copy it in this directory
```bash
git checkout CurrentConstraint
make
cp xspec InputFiles/Verification/currentconstraint/TestCases_Comparison/xspec
```

	3. Run TestCases with `run_local` script (~20 minutes). This will run all L0 and L1 input files with xspec.
```bash
./run_local
```

	4. Run L3 files with `run_global` script (~20 minutes). This will run all L3 files with xspec
```bash
./run_global
```

	5. Run TestCases with `run_master` script (~20 minutes). This will set the reference to which files will be compared to.
```bash
./run_master
```

	5. Compare both outputs using Matlab script
```bash
matlab -nodesktop
compare_tests
```

For each testcase, outputs (from point 3., 4. and 5.) should be the same, since the same physical system is described
(once constrained by (mu, psip) or (iota, oita), and once by the toroidal current).
The matlab script `compare_tests.m` will evaluate the force difference between two cases.
If the estimates for df are small enough (below 1E-14), outputs are considered to be equal.

	6. Test parallelization
```bash
mpirun -np 4 ./xspec_master G3V08L1Fi.master.001
mpirun -np 4 ./xspec G3V08L1Fi.001
mpirun -np 4 ./xspec G3V08L3Fi.001
matlab -nodesktop
compare_spec_outputs('G3V08L3Fi.001.sp.h5', 'G3V08L1Fi.master.001.sp.h5')
compare_spec_outputs('G3V08L1Fi.001.sp.h5', 'G3V08L1Fi.master.001.sp.h5')
exit
```

Again, `compare_spec_outputs.m` will evaluate the force difference between both runs.
If the estimate for `df` is small enough (below 1E-14), outputs are considered to be the same.

If you have the time, repeat point 6. on other input files.

