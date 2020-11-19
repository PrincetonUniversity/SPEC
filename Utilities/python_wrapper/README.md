# Python wrapper for SPEC

This is the python wrapper for SPEC using the package `f90wrap`.

## Prerequisites
You have to install a customized version of `f90wrap` modified by Caoxiang Zhu.
The modifications were made to address the bugs with derived arrays and Fortran variables conflicting with python keywords.

To install the customized version, please do
```
pip install -U git+https://github.com/zhucaoxiang/f90wrap
```

## Wrapped files
The modules and subroutines in `global.f90` and `${SPECFILES}` will be wrapped.
All the input variables and output quantities together with Fortran subroutines will be directly accessable in python.

## Modifications
Here are the modifications we have made in the original SPEC sources.

1. global.f90: `CHARACTER :: ext*100` -> `CHARACTER(100)     :: ext`

## Compile
To compile the python package, you should type
```
make CC=gfortran all
```
The `CC=gfortran` part will be the same compiler option as you use to compile SPEC.

You can check compiler options by using
```
make CC=gfortran compile_test
```

You can also clean the temporary files by using
```
make f90wrap_clean
```
Or clean everything including SPEC objects by using
```
make all_clean
```

*So far, I have only tested GCC compiler, while Intel compiler should also work after updates in the Makefile*


# Compile with CMake in **cmake** branch
0. Go to the folder *SPEC/Utilities/python_wrapper*
```bash
cd <SPEC_ROOT>/Utilities/python_wrapper
```
1. Edit the *cmake_config.json* file in *SPEC/Utilities/python_wrapper* to update the cmake arguments to match with your system.
2. Install scikit-buid, numpy, cmake, ninja using pip:
```
pip install scikit-build, numpy, cmake, ninja
```
3. Install Caoxiang's customized f90wrap:
```
pip install -U git+https://github.com/zhucaoxiang/f90wrap
```
4. Run the following commands in the order
```bash
python setup.py build_ext
python setup.py install
```
5. Test it by running:
```bash
python -c "import spec; print(dir(spec))
```
You should see lot of SPEC subroutines.
