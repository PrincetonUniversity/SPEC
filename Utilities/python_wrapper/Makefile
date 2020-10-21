SOURCES = ../../

# include the original Makefile
include ${SOURCES}/Makefile

#=======================================================================
#                   define the compiler names
#=======================================================================

CC       = gcc
F90      = gfortran
#F90      = ifort
#F90      =  /opt/intel/composer_xe_2015.3.187/bin/intel64/ifort
PYTHON   = python

#=======================================================================
#                     additional flags
#=======================================================================

ifeq ($(F90),gfortran)
	FPP      = $(FC) -E
	FPP_F90FLAGS = -x f95-cpp-input -fPIC
	F90FLAGS = $(RFLAGS)
    FCOMP    = $(FC) #gfortran
    LIBS     =
	CFLAGS  := -fdefault-real-8 -fPIC
	LINKS   := -L$(LAPACK_HOME) -llapack \
			   -L$(HDF5_HOME)/lib -lhdf5hl_fortran -lhdf5 -lhdf5_fortran -lhdf5 -lpthread -lz -lm \
			   -L$(FFTW_HOME)/lib -lfftw3 \
			   -L$(OPENBLAS_HOME)/lib -lopenblas
	LIBS    := -I$(HDF5_HOME)/include -I$(FFTW_HOME)/include
endif

ifeq ($(F90),ifort)

	FPP      = gfortran -E # gfortran f90wrap temp files only. not compilation
	FPP_F90FLAGS = -x f95-cpp-input -fPIC
	F90FLAGS = -fpscomp logicals -fPIC # use 1 and 0 for True and False
    FCOMP    = intelem # for f2py
    LIBS =
endif

# CFLAGS := ${CFLAGS} -fPIC #     ==> universal for ifort, gfortran, pgi

#=======================================================================
#=======================================================================

UNAME = $(shell uname)

ifeq (${UNAME}, Darwin)
  LIBTOOL = libtool -static -o
else
  LIBTOOL = ar src
endif

# ======================================================================
# PROJECT CONFIG, do not put spaced behind the variables
# ======================================================================
# Python module name
PYTHON_MODN = spec
# mapping between Fortran and C types
KIND_MAP = kind_map
# static library name
LIB_NAME = lib$(PYTHON_MODN).a

#=======================================================================
#       List all source files that require a Python interface
#=======================================================================

# names (without suffix), f90 sources
TMP_FILES := global ${SPECFILES}
TMP_FILES := $(filter-out sphdf5,$(TMP_FILES))
LIBSRC_WRAP_SOURCES := $(TMP_FILES)

# file names
LIBSRC_WRAP_FILES = $(addsuffix .f90,${LIBSRC_WRAP_SOURCES})

# macro_explained files
LIBSRC_WRAP_MFILES = $(addsuffix _m.F90,${LIBSRC_WRAP_SOURCES})

# object files
LIBSRC_WRAP_OBJECTS = $(addsuffix _r.o,${LIBSRC_WRAP_SOURCES})

# fpp files
LIBSRC_WRAP_FPP_FILES = $(addsuffix .fpp,${LIBSRC_WRAP_SOURCES})

#=======================================================================
#                 Relevant suffixes
#=======================================================================

.SUFFIXES: .f90 .fpp

#=======================================================================
#
#=======================================================================

.PHONY: all all_clean

all: _${PYTHON_MODN}.so test

all_clean: f90wrap_clean
	make -C ${SOURCES} clean

f90wrap_clean:
	-rm -f ${LIBSRC_OBJECTS} ${LIBSRC_FPP_FILES} $(LIB_NAME) _${PYTHON_MODN}*.so \
	 *.mod *.fpp f90wrap*.f90 f90wrap*.o *.o ${PYTHON_MODN}.py
	-rm -rf ${PYTHON_MODN}_pkg
	-rm -rf src.*/ .f2py_f2cmap .libs/ __pycache__/
	-rm -rf $(LIB_NAME)

$(LIBSRC_WRAP_FPP_FILES): %.fpp: ${SOURCES}/%_m.F90
	${FPP} ${FPP_F90FLAGS} $<  -o $@

$(LIB_NAME): 
	make -C ${SOURCES} xspec CC=gfortran CFLAGS="${CFLAGS}"
	cd ${SOURCES}; ${LIBTOOL} Utilities/python_wrapper/$@ $(ROBJS) $(NUMOBJ)
	ranlib $@

_${PYTHON_MODN}.so: $(LIB_NAME) ${LIBSRC_WRAP_FPP_FILES} 
	f90wrap -m ${PYTHON_MODN} ${LIBSRC_WRAP_FPP_FILES} -k ${KIND_MAP} -v
	f2py-f90wrap --f90exec=${FC} --build-dir . -c -m _${PYTHON_MODN} -L. -l$(PYTHON_MODN) -I${SOURCES} ${LIBS} $(LINKS) f90wrap*.f90

compile_test:
	@echo "SOURCES = ${SOURCES}"
	@echo "FC = ${FC}"
	@echo "LIBTOOL = ${LIBTOOL}"
	@echo "LIBSRC_WRAP_SOURCES = ${LIBSRC_WRAP_SOURCES}"
	@echo "LINKS = ${LINKS}"

test:
	@python -c "import spec; print('SPEC version: {:}'.format(spec.allglobal.version))"