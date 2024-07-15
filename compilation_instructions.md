## SPEC complition instructions

Guide for installing SPEC including the python wrappers, relies on usign cmake

### Compilation

Get the repository and install the necessary compilers and libraries
```
git clone git@github.com:PrincetonUniversity/SPEC.git 
conda create -n "spec_wrapper" python=3.11
conda activate spec_wrapper
conda install gcc_linux-64 gxx_linux-64 gfortran_linux-64
conda install hdf5 openblas libopenblas fftw scalapack openmpi cmake ninja
conda install h5py matplotlib f90nml scipy scikit-build mpi4py ipython
pip install f90wrap
```

Link to the correct `f2py`
```
ln -s ~/anaconda3/envs/spec_wrapper/bin/f2py ~/anaconda3/envs/spec_wrapper/bin/f2py3
```

If necessary, unset HDF5 and FFTW environmental variables
```
unset HDF5, HDF5_ROOT, HDF5_HOME, FFTW, FFTW_DIR
```

Finally, install SPEC and the wrapper (logs will be in `compile.log`)
```
pip install -v . 2>&1 | tee compile.log
```

Install the `py_spec` python library
```
cd Utilities/pythontools/
pip install -e .
```


### Testing your SPEC installation

First, verify that the stand-alone executable is usable.
A few test cases are provided in `InputFiles/TestCases`.

Create a new directory for SPEC runs and change into it

```bash
mkdir ~/SPEC_runs
cd ~/SPEC_runs
```

Copy a demo input file into the current working directory:

```bash
cp ~/SPEC/InputFiles/TestCases/G3V01L0Fi.001.sp .
```

Call SPEC with an input file (`*.sp`) as argument on the command line:

```bash
xspec G3V01L0Fi.001.sp
```

You should see the screen output of the SPEC run.
Among the last lines should be something similar to this:

```
ending :       0.88 : myid=  0 ; completion ; time=      0.88s =     0.01m =   0.00h =  0.00d ; date= 2022/02/17 ; time= 17:35:33 ; ext = G1V02L0Fi.001                                               
ending :            : 
xspech :            :
xspech :       0.88 : myid=  0 : time=    0.01m =   0.00h =  0.00d ;
```

This indicates that the stand-alone executable is usable.

Next, the python wrapper is tested.

1. Check that the SPEC version can be found:
    
    ```bash
    python -c "from spec import spec_f90wrapped as spec; print('SPEC version: {:}'.format(spec.constants.version))"
    ```
    
    This should print a message like "SPEC version: 3.1" on the screen.
    
2. Check that the Python wrapper can be used as a stand-alone code:
    
    ```bash
    OMP_NUM_THREADS=1 python ~/SPEC/Utilities/python_wrapper/spec/core.py G3V01L0Fi.001.sp
    ```
    
    This should conclude with the message `SPEC called from python finished!`.
    
3. Run the optimization example code:
    
    ```bash
    OMP_NUM_THREADS=1 python ~/SPEC/Utilities/python_wrapper/examples/example.py
    ```
    
    This should run a basic optimization problem,
    where the SPEC inputs are controlled via `scipy.optimize`.
    
4. Run the interactive re-convergence example code:
    
    ```bash
    OMP_NUM_THREADS=1 python ~/SPEC/Utilities/python_wrapper/examples/example_2.py
    ```
    
    This should compute a SPEC equilibrium, then change the central pressure,
    re-converge SPEC, etc. for a set of five values of the central pressure
    in a two-volume classical Stellarator case.
    After the pressure scan with re-convergence,
    a plot of the MHD energy vs. the central pressure is shown.