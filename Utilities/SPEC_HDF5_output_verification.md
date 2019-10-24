Verification of the new HDF5 output module in SPEC
======

J. Schilling (jonathan.schilling@ipp.mpg.de), 2019-10-23

Objective: Verify that the newly-developed HDF5 output module actually saves identical information
           as the previous implementation does

Approach:  Run SPEC with both the previous output module and the HDF5 output module on a (number of) test case(s)
           and check whether the outputs are identical

Detailed list of steps:

0. These operations are listed for execution on the COBRA cluster at IPP: https://www.mpcdf.mpg.de/services/computing/cobra
```
> ssh cobra-i.mpcdf.mpg.de
> module purge
> module load intel impi mkl git hdf5-serial fftw-mpi matlab
> module list
Currently Loaded Modulefiles:
  1) intel/19.0.4         2) impi/2019.4          3) mkl/2019.4           4) git/2.16             5) hdf5-serial/1.8.21   6) fftw-mpi/3.3.8       7) matlab/R2019a
```

1. The source code for SPEC is going into a folder called 'src':
```
> mkdir src
```
  
2. build the current state of SPEC on the "issue68" branch:
```
> git clone git@github.com:PrincetonUniversity/SPEC.git -b issue68 src/SPEC_issue68
> pushd src/SPEC_issue68
> make CC=intel_ipp xspec
> popd
```

3. build SPEC from the master branch at the commit from where "issue68" was branched off:
```
> git clone git@github.com:PrincetonUniversity/SPEC.git src/SPEC_master
> pushd src/SPEC_master
> git checkout 32e19c3
> patch Makefile < ../SPEC_issue68/Utilities/adjust_spec_makefile_for_IPP.diff
> make CC=intel_ipp xspec
> popd
```

4. The intermediate files for the test runs are going into a separate folder 'analysis/SPEC_output_comparison'
```
> mkdir -p analysis/SPEC_output_comparison
```

The folder structure should be like this now:

| folder                                  | contents |
| ----------------------------------------|----------|
| /u/jons/src/SPEC_master                 | master branch of SPEC at commit 32e19c3 (right before the branch to issue68) |
| /u/jons/src/SPEC_issue68                | latest state of the issue68 branch |
| /u/jons/analysis/SPEC_output_comparison | outputs of the two SPEC versions for a number of input files |

---
### The general setup for comparing the two SPEC versions is completed.
### Now we are going to run the two versions of SPEC on a number of input files and compare the outputs.
---


## First test case: G3V01L0Fi.002 (W7-X OP1.1)











3. do the test run for the master branch
```
> cp src/SPEC_master/
> pushd G3V01L0Fi.002_master
> ln -s ../../SPEC_master/InputFiles/TestCases/G3V01L0Fi.002.sp .
> ln -s ../../SPEC_master/xspec .
> ln -s ../slurm_spec .
> sbatch slurm_spec ./xspec G3V01L0Fi.002
> ls -lah
-rw-r--r-- 1 jons ipg 7,9K 24. Jun 16:01 8704446.err
-rw-r--r-- 1 jons ipg  982 24. Jun 16:01 8704446.out
lrwxrwxrwx 1 jons ipg   55 24. Jun 14:53 G3V01L0Fi.002.sp -> ../../SPEC_master/InputFiles/TestCases/G3V01L0Fi.002.sp
-rw-r--r-- 1 jons ipg  58K 24. Jun 15:46 .G3V01L0Fi.002.sp.A
-rw-r--r-- 1 jons ipg 118K 24. Jun 16:01 G3V01L0Fi.002.sp.end
-rw-r--r-- 1 jons ipg 2,1M 24. Jun 16:01 .G3V01L0Fi.002.sp.grid
-rw-r--r-- 1 jons ipg  59K 24. Jun 16:01 G3V01L0Fi.002.sp.h5
-rw-r--r-- 1 jons ipg 1,8K 24. Jun 16:01 .G3V01L0Fi.002.sp.its
-rw-r--r-- 1 jons ipg  37M 24. Jun 16:01 .G3V01L0Fi.002.sp.P.0001.dat
-rw-r--r-- 1 jons ipg  524 24. Jun 16:01 .G3V01L0Fi.002.sp.t.0001.dat
lrwxrwxrwx 1 jons ipg   13 24. Jun 15:40 slurm_spec -> ../slurm_spec
lrwxrwxrwx 1 jons ipg   23 24. Jun 13:35 xspec -> ../../SPEC_master/xspec
> popd
```

6. do the test run for the issue68 branch
```
> mkdir G3V01L0Fi.002_issue68
> pushd G3V01L0Fi.002_issue68
> ln -s ../../SPEC_master/InputFiles/TestCases/G3V01L0Fi.002.sp .
> ln -s ../../SPEC_issue68/xspec xspec
> ln -s ../slurm_spec .
> sbatch slurm_spec ./xspec G3V01L0Fi.002.sp
> ls -lah
-rw-r--r-- 1 jons ipg 7,9K 24. Jun 17:44 8705540.err
-rw-r--r-- 1 jons ipg  984 24. Jun 17:44 8705540.out
-rw-r--r-- 1 jons ipg  39M 24. Jun 17:44 G3V01L0Fi.002.h5
lrwxrwxrwx 1 jons ipg   55 24. Jun 17:15 G3V01L0Fi.002.sp -> ../../SPEC_master/InputFiles/TestCases/G3V01L0Fi.002.sp
-rw-r--r-- 1 jons ipg 118K 24. Jun 17:44 G3V01L0Fi.002.sp.end
lrwxrwxrwx 1 jons ipg   13 24. Jun 17:15 slurm_spec -> ../slurm_spec
lrwxrwxrwx 1 jons ipg   24 24. Jun 17:15 xspec -> ../../SPEC_issue68/xspec
> popd
```

7. run comparison routine
```
> module load matlab/R2018b
> matlab -nodesktop
>> addpath('/u/jons/src/SPEC_issue68/Utilities/matlabtools')
>> fdata = read_spec_field('/u/jons/src/SPEC_master/InputFiles/TestCases/G3V01L0Fi.002.sp.h5');
>> gdata = read_spec_grid('/u/jons/src/SPEC_master/InputFiles/TestCases/G3V01L0Fi.002.sp.h5'); 
>> idata = read_spec_iota('/u/jons/src/SPEC_master/InputFiles/TestCases/G3V01L0Fi.002.sp.h5');
>> pdata = read_spec_poincare('/u/jons/src/SPEC_master/InputFiles/TestCases/G3V01L0Fi.002.sp.h5');
>> data = read_spec('/u/jons/src/SPEC_issue68/InputFiles/TestCases/G3V01L0Fi.002.h5');          
>> specheck(fdata, gdata, idata, pdata, data);
... all output quantities
Matching :)
```

11. Repeat steps 5 to 7 for any other input files you wish to test the new output writing routine on

e.g.
```
>> fdata = read_spec_field('/u/jons/src/SPEC_output_comparison/G3V02L1Fi.001_master/G3V02L1Fi.001.sp.h5');
>> gdata = read_spec_grid('/u/jons/src/SPEC_output_comparison/G3V02L1Fi.001_master/G3V02L1Fi.001.sp.h5');
>> idata = read_spec_iota('/u/jons/src/SPEC_output_comparison/G3V02L1Fi.001_master/G3V02L1Fi.001.sp.h5');
>> pdata = read_spec_poincare('/u/jons/src/SPEC_output_comparison/G3V02L1Fi.001_master/G3V02L1Fi.001.sp.h5');
>> data = read_spec('/u/jons/src/SPEC_output_comparison/G3V02L1Fi.001_issue68/G3V02L1Fi.001.h5');
>> specheck(fdata, gdata, idata, pdata, data);
... all output quantities
Matching :)
```

