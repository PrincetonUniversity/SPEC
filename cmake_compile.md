# Installation with CMake

Using cmake, SPEC can be built as a stand-alone executable or as a python extension,
where SPEC can be run directly from python, with all variables passed directly in memory.

Download the package from git. And change to the root directory of SPEC source code by running
```bash
cd <SPEC_ROOT>
```

## Stand-alone Executable Compiling

Compiling SPEC requires MPI, HDF5, and numerical libraries such as BLAS, LAPACK, FFTW. For numerical libraries, you could use system supplied libraries or you could use intel math kernel library (MKL).


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
