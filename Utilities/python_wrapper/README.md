# Python wrapper for SPEC

This is the python wrapper for SPEC using the package `f90wrap`.

## Prerequisites
You have to install a customized version of `f90wrap` modified by Caoxiang Zhu and Jonathan Schilling.
The modifications were made to address the bugs with derived arrays and Fortran variables conflicting with python keywords.
It also fixes a problem with missing `setjmpvalue` variables due to changes in more recent numpy versions.

To install the customized version, please do
```
pip install -U git+https://github.com/jonathanschilling/f90wrap
```

## Wrapped files
The modules and subroutines in `global.f90` and `${SPECFILES}` will be wrapped.
All the input variables and output quantities together with Fortran subroutines will be directly accessable in python.

## Modifications
Here are the modifications we have made in the original SPEC sources.

1. global.f90: `CHARACTER :: ext*100` -> `CHARACTER(100)     :: ext`
2. macros:
   - `CHARACTER` replacement is commented
   - REAL is now replaced with `real(8)`

## Compile
To compile the python package, you should type
```
make BUILD_ENV=gfortran all
```
The `BUILD_ENV=gfortran` part will be the same compiler option as you use to compile SPEC.

You can check compiler options by using
```
make BUILD_ENV=gfortran compile_test
```

You can also clean the temporary files by using
```
make f90wrap_clean
```
Or clean everything including SPEC objects by using
```
make all_clean
```

## Usage
The main SPEC python class is in `core.py`.
You can now use it directly via `mpiexec python core.py ext`, where `ext` is the SPEC extension.
The more advanced way is to use it as an imported module.
All the modules and subroutines that have been interfaced are now available in `SPEC.lib`.
For convenience, the global modules can be directly accessed via `SPEC.allglobal`, `SPEC.inputlist`, etc.
There is an example importing the SPEC python class and optimize the volume inside an interface at `example.py`.
A second example where the in-memory modification and reconvergence of SPEC is demonstrated can be found in `example2.py`.

*So far, I have only tested GCC compiler, while Intel compiler should also work after updates in the Makefile*
