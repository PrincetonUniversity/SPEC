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
    # For standard FFTW library
    # We are using findfftw module from github/egpbos. 
    # PPPL provides fftwconfig.cmake. We could try that first in future

    # To guide cmake to the location of FFT, set FFTW_ROOT to the path of FFTW
    # Ex: For PPPL gcc suite of compilers, 
    # We are using FFTW find module from github.com/egpbos
    # 
    configure_file(downloadFindFFTW.cmake.in findFFTW-download/CMakeLists.txt)
    execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-download )
    if(result)
        message(FATAL_ERROR "CMake step for findFFTW failed: ${result}")
    else()
        message("CMake step for findFFTW completed (${result}).")
    endif()
    execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-download )
    if(result)
        message(FATAL_ERROR "Build step for findFFTW failed: ${result}")
    endif()

    set(findFFTW_DIR ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-src)
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${findFFTW_DIR}")

    find_package(FFTW REQUIRED)
else()   # Not required for OneAPI, but for older  versions
    set(FFTW_INCLUDE_DIRS $ENV{MKLROOT}/include/fftw)
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


