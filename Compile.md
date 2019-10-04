# Compilation hints for SPEC

In order to run SPEC, you need a copy of the HDF5 libraries installed which has
both the Fortran interface and the parallel (MPI I/O) enabled.

## Mac

See e.g. this document for more detailed instructions:
https://support.hdfgroup.org/ftp/HDF5/current/src/unpacked/release_docs/INSTALL_CMake.txt

In short:
1. download `hdf5-1.10.5.tar.gz` from https://www.hdfgroup.org/downloads/hdf5/source-code/

2. extract

`tar xzf hdf5-1.10.5.tar.gz`

3. cd into source folder

`cd hdf5-1.10.5`

4. make a build folder

`mkdir build`

5. cd into build folder

`cd build`

6. run cmake with options for parallel support and Fortran interface (parallel support and C++ interface are not compatible;
so we have to disable the C++ interface)

`cmake -DHDF5_BUILD_FORTRAN:BOOL=ON -DHDF5_ENABLE_PARALLEL:BOOL=ON -DHDF5_BUILD_CPP_LIB:BOLL=OFF ..`

7. actually build the HDF5 library

`make`

This should leave you with a file "hdf5-1.10.5.dmg" or similar, which you can install similar to any other Mac application.
During the build process of SPEC, you then only need to specify the HDF5 folder in the Makefile, which will likely be
`/Applications/HDF_Group/HDF5_1.10.5/...` or similar.
