# Compilation hints for SPEC

This document tries to summarize the steps necessary to setup SPEC on your machine.
Two approaches are discussed.
The first one is the CMake setup.
The second one is the more classical `Makefile` setup.

## CMake and Anaconda

The Anaconda system provides an ecosystem of compilers, precompiled libraries and python packages ready to be used.
The main goal is to decouple the conda environment from the host system,
so that you can use modern software and tools also on machines with outdated local software.

This guide was written while testing the commands on a Debian 9 x86_64 Linux system.
This should be deemed old enough to demonstrate the weirdest errors if something is not under control,
hence facilitating the correctness of these instructions.

### Install Anaconda
The Anaconda installer is available from [the Anaconda website](https://www.anaconda.com/products/individual).
At the time of writing, the URL to the actual file is: https://repo.anaconda.com/archive/Anaconda3-2021.11-Linux-x86_64.sh

Go into folder where the anaconda installer will be downloaded to
and download the Anaconda installer:

```bash
cd ~/Downloads
wget https://repo.anaconda.com/archive/Anaconda3-2021.11-Linux-x86_64.sh
```

Launch the Anaconda installer:

```bash
bash Anaconda3-2021.11-Linux-x86_64.sh
```

 * press ENTER to contine
 * Accept License Agreement: `yes`
 * Confirm the default installation location (here: `/home/IPP-HGW/jons/anaconda3` == `~/anaconda3`)
 * Allow the installer to run `conda init`: `yes`

At this point, the Anaconda installer will modify your `~/.bashrc`.
In order to make these changes take effect, you can logout and log back in,
close and re-open your Terminal window or do:

```bash
source ~/.bashrc
```

Disable the activation of the Anaconda base environment on login:

```bash
conda config --set auto_activate_base False
```

Anaconda works in so-called virtual environments, where environment variables get managed by Anaconda
to setup paths to compilers, libraries and include directories semi-automagically.

For SPEC, it is suggested to create a new conda environment.
A setup script to perform these actions is provided in the SPEC repository.
Thus, we need to clone the SPEC repository now.

### Clone the SPEC repository

In this guide, we assume that your copy of SPEC will be located at `~/SPEC`.

If you want to upload your changes to the SPEC repository,
you should setup your SSH key in your GitHub Settings and use the following
URL for the repository instead: `git@github.com:PrincetonUniversity/SPEC.git`

Clone the repository from GitHub into a folder `SPEC` in your home directory::

```bash
cd ~
git clone https://github.com/PrincetonUniversity/SPEC.git
```

### Setup a Conda Environment for SPEC

The conda environment setup needed to compile and run SPEC is in `setup_conda.sh`.
This scripts takes a specification of the conda packages to install from `spec_conda_env.yml`
and created a conda environment called `spec_env`.
Also, two scripts are created in `etc/conda/activate.d` and `etc/conda/deactivate.d`
of the `spec_env` environment to mask the system's `LD_LIBRARY_PATH`
when entering the conda environment and to restore it to its previous value
when leaving the `spec_env` environment.
Also, the environment variable `FFTW_ROOT` is set to the conda environment
to tell SPEC to use the conda-provided version of FFTW.

Change into the freshly-cloned SPEC repository and run this script:

```bash
cd ~/SPEC
./setup_conda.sh
```

This might take a while, but at the end you should end up with a message similar to this:

```
... lots of stuff above here ...
done
#
# To activate this environment, use
#
#     $ conda activate spec_env
#
# To deactivate an active environment, use
#
#     $ conda deactivate

~/anaconda3/envs/spec_env ~/SPEC
~/SPEC
```

Activate the conda environment for SPEC:

```bash
conda activate spec_env
```

A forked version of `f90wrap` is required to build the SPEC Python wrapper (for now).
Install that next (inside the `spec_env` conda environment !!!):

```bash
pip install -U git+https://github.com/zhucaoxiang/f90wrap
```

The CMake setup is controlled via the `setup.py` Python script from the SPEC repository.
It parses the `CMAKE_ARGS` environment variable provided by conda.

Additional machine-dependent CMake options are loaded from `cmake_config.json`.
This is a soft-link to a file in `cmake_machines`.
Anaconda provides all libraries required to build SPEC,
but we still need to make sure that the `cmake_config.json` link points to 
`cmake_machines/conda_debian.json`.
Note that `conda.json` in `cmake_machines` is outdated, as it includes machine-dependent options (`-DHDF5_ROOT=~/opt/miniconda3/envs/simsopt/`).

Now force the `cmake_config.json` link to point to the correct file:

```bash
ln -sf cmake_machines/conda_debian.json cmake_config.json
```

This concludes the preliminary setup steps and we can progress by starting the build process:

```bash
python setup.py bdist_wheel
```

This step also takes quite a while.
At the end, the SPEC python package (`spec-0.0.1-cp310-cp310-linux_x86_64.whl` or similar) should be available in `dist`.

Install it now:

```bash
pip install dist/*.whl
```

Note that the Python package you just installed also contains the regular stand-alone SPEC executable `xspec`,
which gets installed into `bin/xspec` of your conda environment.
Verify this by calling `which xspec`.
You should get a message similar to:

```
/home/IPP-HGW/jons/anaconda3/envs/spec_env/bin/xspec
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

### Acknowledgements
Here are a few links that proved useful in setting this up:
* https://docs.anaconda.com/anaconda/install/linux/
* https://www.anaconda.com/products/individual
* https://www.rosehosting.com/blog/how-to-install-anaconda-python-on-debian-9/
* https://github.com/RcppCore/Rcpp/issues/770#issuecomment-346716808
* https://computing.docs.ligo.org/conda/compiling/
* https://stackoverflow.com/a/64253999
* https://stackoverflow.com/a/46833531
* https://stackoverflow.com/a/49238956
* https://conda.io/projects/conda-build/en/latest/resources/use-shared-libraries.html









































# BELOW INSTRUCTIONS ARE OUTDATED

# BELOW INSTRUCTIONS ARE OUTDATED

# BELOW INSTRUCTIONS ARE OUTDATED



In order to run SPEC, you need a copy of the HDF5 libraries installed, which has the Fortran interface enabled.

## Installation with CMake

Using CMake, SPEC can be built as a stand-alone executable and as a python extension,
where SPEC can be run directly from python, with all variables passed directly in memory.


Download the package from git. And change to the root directory of SPEC source code by running
```bash
cd <SPEC_ROOT>
```

## Stand-alone Executable Compiling

Compiling SPEC requires MPI, HDF5, and numerical libraries such as BLAS, LAPACK, FFTW. For numerical libraries, you could use system supplied libraries or you could use intel math kernel library (MKL).

Machine-specific settings when building the python wrapper are put into separate `json` files in the `cmake_machines` directory.
For building the regular SPEC executable, the default settings should work.

In order to select a machine-specific settings file, create a soft link to the indented file in `cmake_machines`:

```bash
ln -sf cmake_machines/gfortran_ubuntu.json cmake_config.json
```

###  CentOS
Here instructions are given for CentOS 7

#### Dependencies
Install OpenBLAS, FFTW3, and hdf5 using the command
```bash
yum install -y  gcc-gfortran openmpi openmpi-devel hdf5 hdf5-devel fftw3 fftw3-devel openblas openblas-devel python3 python3-devel cmake ninja-build
```
If you don't have the latest version of cmake avaialable on your system, you can create a python virtual environment ([instructions are here](https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/)), activate it, and then install cmake in that virtual environment using pip
```bash
pip install cmake ninja
```

#### Configure
When using cmake to build SPEC, the first step is to configure compilers and the locations of libraries. Cmake can detect compilers and libraries at standard locations easily but needs hand-holding when the required libraries are non-standard locations.

The following command was used to configure cmake build setup for SPEC on Centos
```bash
cmake -S. -Bbuild -GNinja -DCMAKE_Fortran_COMPILER=mpifort -DBLA_VENDOR=OpenBLAS -DHDF5_NO_FIND_PACKAGE_CONFIG_FILE=TRUE -DHDF5_PREFER_PARALLEL=TRUE -DCMAKE_INSTALL_PREFIX=${SPEC_ROOT}/install --trace-source=CMakeLists.txt 2>&1 | tee log
```
There are few points to note on the above command
  - All the build related files will be in build folder.
  - Ninja build system is used. If your system doesn't have ninja installed, remove the -G option. The default is the standard `make` tool. 
  - We ae using OpenBLAS for BLAS and LAPACK and MPI fortran compiler
  - Since most of the libraries are in standard location, we don't have to specify them. We are giving couple of options related to HDF5 libraries.
  - The installation path is install subfolder location within SPEC folder.
  - We are interested in a verbose output and also want to store the output in `log` file. 

#### Build
After successful completion of cmake configuration step, building is trivial
```bash
cmake --build build
````
#### Install
The last step is to install the executable by running
```bash
cmake --install build
```

That's it! If all the above steps completed without errors, you have the SPEC executable `xspec`  installed at  `install/bin` folder

## Python Extension Compiling
Building the SPEC python extension will also build the SPEC executable.
In the SPEC root folder, edit the `cmake_config.json` as necessary for your system. Few example `.json` files are provided in the `cmake_machines` folder. 

### Dependencies
It is strongly suggested to use a python virtual environment either conda or python venv. After virtual environment is installed and activated, install the python related dependencies. Please note that these are in addition to the dependencies listed earlier in stand-alone installation steps. If you are using conda virtual environment try installing the dependencies using `conda install` command
```bash 
conda  install -n <your_venv> numpy f90nml scikit-build cmake ninja
```

If you are using venv virtual environment, run
```bash
pip install numpy f90nml scikit-build cmake ninja
```
Now install f90wrap. Please keep in mind that numpy has to be installed before installing f90wrap.
```bash
pip install -U git+https://github.com/zhucaoxiang/f90wrap
```

Now install the SPEC extension by running the setup.py script present in the SPEC root folder. 
```bash
python setup.py bdist_wheel; cd dist/;  pip install *.whl
``` 
in succession. At this point, you should be able to import the `spec` module in python. To test this, you can try the following command from the shell:
```bash
python -c "import spec; print('success')"
```

If you want editable install, run
```
python setup.py develop 
```


## Stellar cluster at PPPL

### Python wrapper
Below are the steps to build python wrappers for SPEC on stellar.

1. Needed modules are 
   > i. hdf5/gcc/1.10.6 
   > ii. intel-mkl/2021.1.1 
   > iii. openmpi/gcc/4.1.0 
   > iv. anaconda3/2021.5. 
   ---
   **Note**

   FFTW is supplied as part of Intel MKL and we just need to link against MKL.

   ---
   Load the modules by running
   ```
   module load hdf5/gcc/1.10.6 intel-mkl/2021.1.1 openmpi/gcc/4.1.0 anaconda3/2021.5
   ```
2. Create conda virtual environment.
   ```
   conda create -n spec_ve python=3.8
   ```
   You have to press enter twice. Here a conda virtual environment named `spec_ve` is created with python version 3.8 and lot of packages are installed. Activate by running
   ```
   conda activate spec_ve
   ```
3. Install `cmake`, `ninja`, `scikit-build`, `numpy` using either conda or pip. 
   ```
   conda install cmake ninja scikit-build numpy
   ``` 
   or
   ```
   pip install cmake ninja scikit-build numpy
   ```
4. Install `f90wrap` by running 
   ```
   pip install git+https://github.com/zhucaoxiang/f90wrap.git
   ```
5. Clone the spec repo from github
   ```
   git clone https://github.com/PrincetonUniversity/SPEC.git
   ```
   Change the working directory by running `cd SPEC`.
6. Edit the cmake_config.json to populate correct cmake_flags. For stellar, cmake_config.json should look like 
   ```
   {
     "cmake_args": [
       "-DCMAKE_C_COMPILER=mpicc",
       "-DCMAKE_CXX_COMPILER=mpicxx",
       "-DCMAKE_Fortran_COMPILER=mpifort",
       "-DBLA_VENDOR=Intel10_64lp",
       "-DHDF5_ROOT=/usr/local/hdf5/gcc/1.10.6",
       "-DHDF5_PREFER_PARALLEL=False"]
   }
   ``` 
7. Then build the python wheel for SPEC wrapper using
   ```
   python setup.py bdist_wheel
   ```
   The resulting wheel is located in `dist` folder. Install SPEC python wrapper by running 
   ```
   pip install dist/spec*.whl
   ```

8. Install mpi4py using pip/conda. If using pip, don't forget to use `--no-cache-dir` flag
   ```
   pip install --no-cache-dir mpi4py
   ```
   or
   ```
   conda install mpi4py
   ```

### SPEC executable.
The python wrapper builds spec executable but it gets installed at an obscure location. If you mainly want SPEC executable `xspec`, the steps are similar.
1. Load the required modules. Refer to the first step in the python wrapper instructions.
2. Clone the SPEC repo and make SPEC as working directory. Refer to the 5th step above.
3. Run the cmake configuration by running
   ```
   cmake -Bbuild -S . -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_Fortran_COMPILER=mpifort -DBLA_VENDOR=Intel10_64lp \
   -DHDF5_ROOT=/usr/local/hdf5/gcc/1.10.6 -DHDF5_PREFER_PARALLEL=False -DCMAKE_INSTALL_PREFIX=<SPEC_install_location>
   ```
   Please note SPEC gets installed at `<SPEC_install_location>/bin`, where `<SPEC_install_location>` is the folder of your choice. Building of SPEC library will be done in the folder `build`, where all the intermediary compilation files will be located.
4. Compile the code by running
   ```
   cmake --build build
   ```
   This command will invoke `make` build generator. Alternatively, you can switch to build folder and run make utility manually.
   ```
   cd build
   make
   ```
5. Install the SPEC executable by running
   ```
   cmake --install build
   ```
   SPEC library gets installed `<SPEC_install_location>/lib` and SPEC executable get installed at `<SPEC_install_location>/bin`

## Mac

Here is how to build the HDF5 library :
1. download `hdf5-1.10.5.tar.gz` from https://www.hdfgroup.org/downloads/hdf5/source-code/
2. extract: `tar xzf hdf5-1.10.5.tar.gz`
3. cd into source folder: `cd hdf5-1.10.5`
4. make a build folder: `mkdir build`
5. cd into build folder: `cd build`
6. run cmake with options for the Fortran interface: `cmake -DHDF5_BUILD_FORTRAN:BOOL=ON ..`
7. actually build the HDF5 library: `make`

This should leave you with a file "hdf5-1.10.5.dmg" or similar, which you can install just as any other Mac application.
See e.g. this document for more detailed instructions:
https://support.hdfgroup.org/ftp/HDF5/current/src/unpacked/release_docs/INSTALL_CMake.txt

The compilation of SPEC itself then proceeds as usual.
You then only need to specify the HDF5 folder in the Makefile, which will likely be
`/Applications/HDF_Group/HDF5/1.10.5`.
