SOURCES = ../..

# include the original Makefile
include ${SOURCES}/SPECfile

#=======================================================================
#                   define the compiler names
#=======================================================================

CC       = gcc

# Fortran compiler
F90      = gfortran
#F90      = ifort # use this if you are using Intel compiler

PYTHON   = python3

#=======================================================================
#                     additional flags
#=======================================================================

ifeq ($(F90),gfortran)
  FPP      = $(F90) -E
  FPP_F90FLAGS = -x f95-cpp-input -fPIC	
endif

ifeq ($(F90),ifort)
  # gfortran f90wrap temp files only. not compilation
  FPP      = gfortran -E
  FPP_F90FLAGS = -x f95-cpp-input -fPIC
endif

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
PYTHON_MODN = spec_f90wrapped
# mapping between Fortran and C types
KIND_MAP = kind_map
# static library name
LIB_NAME = lib$(PYTHON_MODN).a

#=======================================================================
#       List all source files that require a Python interface
#=======================================================================

# names (without suffix), f90 sources
TMP_FILES := ${ALLSPEC} xspech
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
	echo "Make xspec..."
	make -C ${SOURCES} BUILD_ENV=$(BUILD_ENV) xspec
	cd ${SOURCES}; ${LIBTOOL} Utilities/python_wrapper/$@ $(addsuffix _r.o,$(ALLFILES))
	ranlib $@

_${PYTHON_MODN}.so: $(LIB_NAME) ${LIBSRC_WRAP_FPP_FILES}
	f90wrap -m ${PYTHON_MODN} ${LIBSRC_WRAP_FPP_FILES} -k ${KIND_MAP} -v
	mv ${PYTHON_MODN}.py spec/
	f2py-f90wrap --f90exec=${FC} --fcompiler=gnu95 --build-dir . ${LINKS} -L. -l$(PYTHON_MODN) -I${SOURCES} ${LIBS} $(LINKS) -c f90wrap_*.f90 -m _${PYTHON_MODN}
	mv _${PYTHON_MODN}*.so spec/

compile_test:
	make -C ${SOURCES} BUILD_ENV=$(BUILD_ENV) show_makeflags
	@echo "SOURCES = ${SOURCES}"
	@echo "BUILD_ENV = ${BUILD_ENV}"
	@echo "CC = ${CC}"
	@echo "FC = ${FC}"
	@echo "LIBTOOL = ${LIBTOOL}"
	@echo "LIBSRC_WRAP_SOURCES = ${LIBSRC_WRAP_SOURCES}"
	@echo "CFLAGS = ${CFLAGS}"
	@echo "RFLAGS = ${RFLAGS}"
	@echo "FLAGS = ${FLAGS}"
	@echo "LINKS = ${LINKS}"

test:
	@${PYTHON} -c "import spec.spec_f90wrapped as spec; print('SPEC version: {:}'.format(spec.constants.version))"
