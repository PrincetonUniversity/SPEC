#.rst:
#
# The purpose of the F90Wrap –Fortran to Python interface generator– project is to provide a
# connection between Python and Fortran languages.
#
# F90Wrap is a Python package (with a command line tool f90wrap and ) that facilitates
# creating/building Python C/API extension modules that make it possible to call Fortran 77/90/95
# external subroutines and Fortran 90/95 module subroutines as well as C functions; to access Fortran
# 77 COMMON blocks and Fortran 90/95 module data, including allocatable arrays from Python.
#
#
# The following variables are defined:
#
# ::
#
#   F90Wrap_EXECUTABLE      - absolute path to the F90Wrap executable
#
# ::
#
#   F90Wrap_VERSION_STRING  - the version of F90Wrap found
#   F90Wrap_VERSION_MAJOR   - the F90Wrap major version
#   F90Wrap_VERSION_MINOR   - the F90Wrap minor version
#   F90Wrap_VERSION_PATCH   - the F90Wrap patch version
#   F90Wrap_F2PY_EXECUTABLE - absolute path to the F90Wrap supplied f2py
#
#
# .. note::
#
#

# .. warning::
#
#

find_program(F90Wrap_EXECUTABLE NAMES f90wrap)
find_program(F90Wrap_F2PY_EXECUTABLE NAMES f2py-f90wrap)

set(F90Wrap_INCLUDE_DIRS "${F2PY_INCLUDE_DIR}" "${NumPy_INCLUDE_DIRS}")

if(F90Wrap_EXECUTABLE)
  # extract the version string
  execute_process(COMMAND "${F90Wrap_EXECUTABLE}" -V
                  OUTPUT_VARIABLE F90Wrap_VERSION_STRING
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  #if("${F90Wrap_VERSION_STRING}" MATCHES "^([0-9]+)(.([0-9+]))?(.([0-9+]))?$")
  #  set(F90Wrap_VERSION_MAJOR ${CMAKE_MATCH_1})
  #  set(F90Wrap_VERSION_MINOR "${CMAKE_MATCH_3}")
  #  set(F90Wrap_VERSION_PATCH "${CMAKE_MATCH_5}")
  #endif()
endif()

# handle the QUIETLY and REQUIRED arguments and set F90Wrap_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(F90Wrap
  REQUIRED_VARS F90Wrap_EXECUTABLE F90Wrap_F2PY_EXECUTABLE 
  #VERSION_VAR F90Wrap_VERSION_STRING
  )

mark_as_advanced(F90Wrap_EXECUTABLE)
mark_as_advanced(F90Wrap_F2PY_EXECUTABLE)

#include(UseF90Wrap)
#include(FunctionsF90Wrap)
