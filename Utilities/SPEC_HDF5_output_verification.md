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

5a. create a (generic) SLURM batch script for running SPEC on 1 CPU
```
> pushd analysis/SPEC_output_comparison
> cat > slurm_spec_1 << EOF
#!/bin/bash -l
# Standard output and error:
#SBATCH -o ./spec_stdout.%j
#SBATCH -e ./spec_stderr.%j
# Initial working directory:
#SBATCH -D ./
# Job Name:
#SBATCH -J SPEC
# Queue (Partition):
#SBATCH --partition=medium
# Number of nodes and MPI tasks per node:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16384 # 16GB could be enough for the beginning...
#
#SBATCH --mail-type=none
#SBATCH --mail-user=<userid>@rzg.mpg.de
#
# Wall clock limit:
#SBATCH --time=04:00:00

# Run the program:
srun \$@
EOF
```

5b. enqueue the SPEC run from the master branch
```
> mkdir G3V01L0Fi.002_master
> pushd G3V01L0Fi.002_master
> ln -s ../../../src/SPEC_master/InputFiles/TestCases/G3V01L0Fi.002.sp .
> ln -s ../../../src/SPEC_master/xspec .
> ln -s ../slurm_spec_1 .
> ls -lh # check that everything is there for running SPEC
lrwxrwxrwx 1 jons ipg   62 Oct 24 23:06 G3V01L0Fi.002.sp -> ../../../src/SPEC_master/InputFiles/TestCases/G3V01L0Fi.002.sp
lrwxrwxrwx 1 jons ipg   15 Oct 24 23:06 slurm_spec_1 -> ../slurm_spec_1
lrwxrwxrwx 1 jons ipg   30 Oct 24 22:38 xspec -> ../../../src/SPEC_master/xspec
> sbatch ../slurm_spec_1 ./xspec G3V01L0Fi.002
```

Check the output files ```spec_stdout.<jobid>``` and ```spec_stderr.<jobid>```.
In the beginning, the flags should include ```-O0```:
```
xspech :            : 
       :  compiled  : date    = Thu Oct 24 22:31:32 CEST 2019 ; 
       :            : dir     = /u/jons/src/SPEC_master ; 
       :            : macros  = macros ; 
       :            : f90     =  ; 
       :            : flags   =  -r8 -O0 -ip -no-prec-div -xHost -fPIC ; 
                                     ^^^
```

Finally, leave this run until it is finished and return to the ```SPEC_output_comparison``` directory:

```
> popd
```

5c. enqueue SPEC from the issue68 branch on the input file from the master branch
```
> mkdir G3V01L0Fi.002_issue68
> pushd G3V01L0Fi.002_issue68
> ln -s ../../../src/SPEC_master/InputFiles/TestCases/G3V01L0Fi.002.sp .
> ln -s ../../../src/SPEC_issue68/xspec .
> ln -s ../slurm_spec_1 .
> ls -lh # check that everything is there for running SPEC
lrwxrwxrwx 1 jons ipg   62 Oct 24 23:02 G3V01L0Fi.002.sp -> ../../../src/SPEC_master/InputFiles/TestCases/G3V01L0Fi.002.sp
lrwxrwxrwx 1 jons ipg   15 Oct 24 23:01 slurm_spec_1 -> ../slurm_spec_1
lrwxrwxrwx 1 jons ipg   31 Oct 24 23:03 xspec -> ../../../src/SPEC_issue68/xspec
> sbatch slurm_spec_1 ./xspec G3V01L0Fi.002.sp
```

Check the output files ```spec_stdout.<jobid>``` and ```spec_stderr.<jobid>```.
In the beginning, the flags should include ```-O0```:
```
xspech :            : version =  1.90
       :  compiled  : date    = Thu Oct 24 22:30:49 CEST 2019 ; 
       :            : srcdir  = /u/jons/src/SPEC_issue68 ; 
       :            : macros  = macros ; 
       :            : fc      = mpiifort ; 
       :            : flags   =  -r8 -O0 -ip -no-prec-div -xHost -fPIC ;
                                     ^^^
```

Finally, leave this run until it is finished and return to the ```SPEC_output_comparison``` directory:

```
> popd
```

5d. Compare the outputs of the two runs in Matlab






## Second test case: G2V32L1Fi.001 (Tokamak with 32 volumes)

6a. create a (generic) SLURM batch script for running SPEC on 32 CPUs
```
> mkdir analysis/SPEC_output_comparison/G2V32L1Fi.001_master
> pushd analysis/SPEC_output_comparison/G2V32L1Fi.001_master
> cat > slurm_spec_32 << EOF
#!/bin/bash -l
# Standard output and error:
#SBATCH -o ./spec_stdout.%j
#SBATCH -e ./spec_stderr.%j
# Initial working directory:
#SBATCH -D ./
# Job Name:
#SBATCH -J SPEC
# Queue (Partition):
#SBATCH --partition=medium
# Number of nodes and MPI tasks per node:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
# #SBATCH --mem=16384 # 16GB could be enough for the beginning...
#
#SBATCH --mail-type=none
#SBATCH --mail-user=<userid>@rzg.mpg.de
#
# Wall clock limit:
#SBATCH --time=04:00:00

# Run the program:
srun \$@
EOF
```











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
