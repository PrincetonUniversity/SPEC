cmake_minimum_required(VERSION 3.20.0)

project(spec LANGUAGES C Fortran)

set(CMAKE_VERBOSE_MAKEFILE ON CACHE BOOL "ON")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Find OpenMP
find_package(OpenMP)

# Find Threads
set (THREADS_PREFER_PTHREAD_FLAG ON) # Prefer pthreads
find_package(Threads)

# Find libm
set(POW_LIBS "")
include(CheckLibraryExists)
check_library_exists(m pow "" LIBM)
if(LIBM)
    list(APPEND POW_LIBS "m")
endif()


# Find & Configure MPI
# If using the compiler wrapper, there is no need to find the MPI libraries.
get_filename_component (compiler ${CMAKE_Fortran_COMPILER} NAME)
message(STATUS "compiler is ${compiler}")
if (${compiler} STREQUAL mpiifort OR ${compiler} STREQUAL mpifort)
else ()
    find_package (MPI REQUIRED)
    if (MPI_Fortran_FOUND)
        set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_Fortran_COMPILE_FLAGS}")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_Fortran_LINK_FLAGS}")
        include_directories (${MPI_Fortran_INCLUDE_PATH})

#  Some of the tests need to be oversubscribbed to run. There is a problem where
#  depending on the current platform or implementation of MPI the oversubscribe
#  flag changes. MPIEXEC may not write to stdout correctly. Add the ability to
#  set this variable manually.
        if (NOT DEFINED ${MPI_OVERSUBSCRIBE_FLAG}) 
            execute_process(COMMAND ${MPIEXEC} --version OUTPUT_VARIABLE MPI_VERSION)
            if (${MPI_VERSION} MATCHES "Open MPI" OR
                ${MPI_VERSION} MATCHES "OpenRTE"  OR
                ${MPI_VERSION} MATCHES "slurm"    OR
                ${MPI_VERSION} MATCHES "aprun")
                set (MPI_OVERSUBSCRIBE_FLAG "--oversubscribe")
            endif ()
        endif ()
    endif ()
endif ()


if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DIFORT")
endif()

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER "10.0")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
    endif()
endif()

find_program(AWK awk mawk gawk)

# Find LAPACK
# For LAPAPCK, it is suggested to use Intel MKL even for gcc compiler suite. To use Intel MKL, INTELMKL 
# has to be defined as environment variable. If Intel compilers are loaded, MKLROOT is auto-populated. 
# Otherwise define MKLROOT to point to Intel MKL.
# Ex: For PPPL gcc suite of compilers, after loading gcc, on command line set MKLROOT as
# "export MKLROOT=/usr/pppl/intel/2020.u1/compilers_and_libraries_2020/linux/mkl"
# If Intel MKL is not desired, change BLA_VENDOR variabale during configuration
#set(BLA_VENDOR Intel10_64lp CACHE STRING 
#    "Define BLAS vendor. Some of the popular options are Intel10_64lp OpenBLAS Apple Generic")
# find_package (BLAS REQUIRED)
# message(STATUS "BLAS vendor is ${BLA_VENDOR}")
# if (NOT ${BLA_VENDOR} MATCHES "^Intel" AND NOT ${BLA_VENDOR} MATCHES "^OpenBLAS")
find_package (LAPACK REQUIRED)
#endif()

# Find FFTW3.
# Intel oneAPI has FFT3 available as part of MKL. Just linking against MKL gives FFTW functionality.
if (NOT ${BLA_VENDOR} MATCHES "^Intel")
  find_package(FFTW3 REQUIRED)

  # will find it, but FFTW_LIBRARIES is set to 'fftw3' simple string, which F2py does not like
  
  # Hack to find the proper .so looking first in conda, then in other 
  set(FFTW3_LIBRARY_PATHS) # Initialize empty list

  if(DEFINED CONDA_PREFIX)
    list(APPEND FFTW3_LIBRARY_PATHS "${CONDA_PREFIX}/lib") # Add Conda first
  endif()

  list(APPEND FFTW3_LIBRARY_PATHS
    /usr/lib
    /usr/local/lib
    ${FFTW3_ROOT_DIR}/lib
  )
  # also add a path relative to the found include dirs, which is a true path
  if(FFTW3_INCLUDE_DIRS)
    string(REPLACE "/include" "/lib" FFTW3_LIBRARY_DIR_FROM_INCLUDE "${FFTW3_INCLUDE_DIRS}")
    list(APPEND FFTW3_LIBRARY_PATHS ${FFTW3_LIBRARY_DIR_FROM_INCLUDE})
  endif()

  # Find the library using the constructed list of paths. Should find the path to the .so
  find_library(FFTW3_LIBRARY_TRUE
    NAMES fftw3
    PATHS ${FFTW3_LIBRARY_PATHS}
    NO_DEFAULT_PATH
  )

  if(FFTW3_LIBRARY_TRUE)
    message(STATUS "FFTW3 Library with full path found: ${FFTW3_LIBRARY_TRUE}")
    set(FFTW3_FULL_PATH ${FFTW3_LIBRARY})
  else()
    message(FATAL_ERROR "Could not find FFTW3 library with find_library.")
  endif()


else()   # Not required for OneAPI, but for older  versions
    set(FFTW3_INCLUDE_DIRS $ENV{MKLROOT}/include/fftw)
endif()

# Find HDF5
# To guide cmake finding HDF5 set HDF5_ROOT. If no cmake configuration file is available, also set
# HDF5_NO_FIND_PACKAGE_CONFIG_FILE
# Ex: For PPPL gcc suite of compilers, after loading HDF5, invoke cmake as 
# cmake ... -DHDF5_ROOT=$HDF5_HOME -DHDF5_NO_FIND_PACKAGE_CONFIG_FILE=FALSE
find_package(HDF5 COMPONENTS C Fortran HL REQUIRED)
message(STATUS "HDF5 libraries found: ${HDF5_FOUND}")
message(STATUS "HDF5 C include directors : ${HDF5_C_INCLUDE_DIRS}")
message(STATUS "HDF5 Fortran include directors : ${HDF5_Fortran_INCLUDE_DIRS}")

# Direct the location of build files
set (EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/build/bin)
set (LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/build/lib)

# Fortran specific things
set (CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/build/modules)

# Build SPEC
#add_library(spec_lib "")
#add_executable(xspec "")

add_subdirectory(src)
if(SKBUILD)
    add_subdirectory("Utilities/python_wrapper")
endif()

#tartget_link_libraries(xspec spec_lib)
# Add dependencies to other libraries
#target_link_libraries(spec_lib PUBLIC
#                      BLAS LAPACK FFTW HDF5::C HDF5::Fortran)


