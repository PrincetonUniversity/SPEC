# SPEC compilation instructions

The default installation method for SPEC uses CMake and installs
the python wrappers and an xspec executable. 

## Installation using Anaconda

We recommend you use Anaconda to create a coherent build environment and prevent 
dependency conflicts.

Control over the installation can be had by editing `cmake_config.json`, to guide 
CMake to the right compilers etc. 
Configurations for different machines are stored in `${SPEC_ROOT}/cmake_machines`,
to use these, link them to cmake_config.json: `ln -s cmake_config cmake_machines/<config_file.json>`

>[!TIP]
>install as much as possible in your environment using the `conda` command, 
>only use 'pip' at the very end for the last packages. 
>if you have not added the `conda-forge` channel do so by
>`conda config --add channels conda-forge`

Get the repository and install the necessary compilers and libraries
```bash
git clone git@github.com:PrincetonUniversity/SPEC.git 
conda create -n "spec_wrapper" python=3.11 # create your environment for SPEC
conda activate spec_wrapper
conda install gcc_linux-64 gxx_linux-64 gfortran_linux-64 # or macOS versions, see note below
conda install hdf5 openblas libopenblas fftw scalapack openmpi cmake ninja
conda install h5py matplotlib f90nml scipy scikit-build mpi4py ipython
pip install f90wrap
```

>[!NOTE]
> for macOS users use the respective compiler packages; 
> `conda install clang_osx-64 clangxx_osx-64 gfortran_osx-64`


Finally, install SPEC and the wrapper (logs will be in `compile.log`)
```
pip install -v . 2>&1 | tee compile.log
```

Install the `py_spec` python library
```
cd Utilities/pythontools/
pip install -e .
```

### Troubleshooting Anaconda install
If using a newer version of python, `f2py3` is no longer shipped. If your system contains an old python install (for example from your OS), CMake can find its `f2py3` and give try to use it to compile the wrappers instead of your environments `f2py`. 
Test this by looking if you have an `f2py3` in your path: `$which f2py3`. 
The easiest workaround is to create a link called f2py3 that links to f2py so it is found first. 
```
ln -s ~/anaconda3/envs/spec_wrapper/bin/f2py ~/anaconda3/envs/spec_wrapper/bin/f2py3
```

You might have HDF5 or FFTW environment variables set (for example for a VMEC install). This can throw off CMake, which we want to use only anaconda. 
```
unset HDF5, HDF5_ROOT, HDF5_HOME, FFTW, FFTW_DIR
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


## Other legacy installations
It is still possible to compile SPEC using `make` or `cmake` directly, and bypass the wrapper installation. 

### CMake installation
Spec can be installed using CMake to find the relevant libraries to link against. 
You can control 
in the root directory of SPEC do the following: 
```bash
mkdir build
cd build
cmake ..
make
```
This will compile SPEC (not the wrappers). The `xspec` executable is found in ${SPEC_ROOT}/build/build/bin/xspec

### Make installation
SPEC can also be installed using the `make` command in the root directory. 

The `make` install is controlled by the `BUILD_ENV` environment variable. 
Available options are found in the `SPECfile`
where different link and compile flags for many machines are found. 

If you cannot find your machine in the list, copy a similar machine and adapt as needed. 
Then compile by running the command

```bash
BUILD_ENV=<machine_name> make
```

The `make` process creates files in the SPEC_ROOT directory, and creates the `xspec` executable there. 


## Build process
the source files are found in the `${SPEC_ROOT}/src/ directory`.
The `.f90` files contain macros that are expanded during the make process using the `m4` command. 

Depending on the build type, the macro-expanded code is either found in `build/src/`, in the root directory, or in the `_skbuild` folder. 

>[!TIP]
>The line numbers in error messages correspond to the macro-expanded code

The macros are defined in `src/macros`
