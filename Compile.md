# Compilation hints for SPEC

In order to run SPEC, you need a copy of the HDF5 libraries installed which has
both the Fortran interface and the parallel (MPI I/O) enabled.

# Installation with CMake

Using cmake, SPEC can be built as a stand-alone executable or as a python extension,
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
