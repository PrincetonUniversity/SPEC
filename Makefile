#!/bin/sh

###############################################################################################################################################################

# This is the "classic" Makefile for SPEC.

# The hierarchy among the source tree is roughly as follows:
# 1. inputlist contains the input namelists for SPEC and a routine to initialize their default values.
# 2. global    contains the global workspace ("state") of the SPEC computation as well as the readin() and preset() routines
#                       required to initialize the workspace based on input variables.
# 3. sphdf5    contains the HDF5 output file writing routines. It depends on global and inputlist and therefore needs to be separately listed.
# 4. afiles ..
#    .. ffiles contains the main physics wisdom of SPEC.
# 5. sfiles    contains externally provided libraries, e.g. for quadrature, minimization, etc.
# 6. xspech    contains the stand-alone main SPEC program.

###############################################################################################################################################################

 # basis of SPEC: input variables, global workspace, HDF5 output file writing
 # these are split off since they require special treatment (needed by all others and/or special macros)
 BASEFILES=inputlist global
 IOFILES=sphdf5

 # (most of) physics part of SPEC
 afiles=preset manual rzaxis packxi volume coords basefn memory
 bfiles=metrix ma00aa matrix spsmat spsint mp00ac ma02aa packab tr00ab curent df00ab lforce intghs mtrxhs lbpol
 cfiles=brcast dfp100 dfp200 dforce newton
 dfiles=casing bnorml
 efiles=jo00aa pp00aa pp00ab bfield stzxyz
 ffiles=hesian ra00aa numrec

 # externally provided libraries
 # below assumes the .f files are double precision; the CFLAGS = -r8 option is not required;
 sfiles=dcuhre minpack iqpack rksuite i1mach d1mach ilut iters

###############################################################################################################################################################

 # all of SPEC except BASEFILES
 SPECFILES=$(afiles) $(bfiles) $(cfiles) $(dfiles) $(efiles) $(ffiles)

 # all of "our" (vs. contributed) files needed for the "core" of SPEC
 ALLSPEC=$(BASEFILES) $(IOFILES) $(SPECFILES)

 # *ALL* files needed for the main SPEC executable
 ALLFILES=$(sfiles) $(ALLSPEC) xspech

 # to be preprocessed by m4
 PREPROC=$(ALLSPEC:=_m.F90)

 ROBJS=$(SPECFILES:=_r.o)
 DOBJS=$(SPECFILES:=_d.o)

 ROBJS_IO=$(IOFILES:=_r.o)
 DOBJS_IO=$(IOFILES:=_d.o)

###############################################################################################################################################################

 MACROS=macros

 # if want to use gfortran: make BUILD_ENV=gfortran (x/d)spec
 # default: use Intel compiler
 BUILD_ENV?=intel

 # to enable OpenMP acceleration within volume, set OMP=yes, otherwise set OMP=no
 OMP=yes

ifeq ($(BUILD_ENV),intel)
 # Intel Defaults
 # At PPPL, you can use the following commands
 # module use /p/focus/modules
 # module load spec
 FC=mpif90
 CFLAGS=-r8 -DIFORT
 RFLAGS=-mcmodel=large -O3 -m64 -unroll0 -fno-alias -ip -traceback
 DFLAGS=-O0 -g -traceback -check bounds -check format -check output_conversion -check pointers -check uninit -debug full -D DEBUG
 LIBS=-I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include  # MKL include
 LIBS+=-I$(HDF5_HOME)/include # HDF5 include
 LINKS=-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lhdf5hl_fortran # HDF5 link
 LIBS+=-I$(FFTW_HOME)/include # FFTW include
 LINKS+=-L$(FFTW_HOME)/lib -lfftw3 # FFTW link
 LINKS+=${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a \
     -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a \
     ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl # MKL link
endif

ifeq ($(BUILD_ENV),gfortran)
 # At PPPL, you can use the following commands
 # module use /p/focus/modules
 # module load spec/gcc
 FC=mpif90
 FLAGS=-fPIC
 RFLAGS=-O3 -w -ffree-line-length-none -fexternal-blas # -fallow-argument-mismatch # only used for GCC-10
 DFLAGS=-O0 -g -w -ffree-line-length-none -Wextra -Wtarget-lifetime -fbacktrace -fbounds-check -fexternal-blas \
     -fcheck=all -DDEBUG #-ffpe-trap=invalid,zero,overflow,underflow,inexact # for some reason this will cause crash
 CFLAGS=-fdefault-real-8
 LINKS=-L$(LAPACK_HOME) -llapack #-lgfortran
 LIBS=-I$(HDF5_HOME)/include
 LINKS+=-L$(HDF5_HOME)/lib -lhdf5hl_fortran -lhdf5 -lhdf5_fortran -lhdf5 -lpthread -lz -lm
 LIBS+=-I$(FFTW_HOME)/include
 LINKS+=-L$(FFTW_HOME)/lib -lfftw3
 LINKS+=-L$(OPENBLAS_HOME)/lib -lopenblas
endif

ifeq ($(BUILD_ENV),gfortran_ubuntu)
 # You should install the following packages
 # sudo apt install gfortran
 # sudo apt install libopenmpi-dev
 # sudo apt install liblapack-dev
 # sudo apt install m4
 # sudo apt install libfftw3-dev
 # sudo apt install libhdf5-dev
 FC=mpif90
 FLAGS=-fPIC
 CFLAGS=-fdefault-real-8
 LINKS=-Wl,-rpath -Wl,/usr/lib/lapack -llapack -lblas
 LIBS=-I/usr/include/hdf5/serial
 LINKS+=-L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_fortran -lhdf5 -lpthread -lz -lm
 LIBS+=-I/usr/include
 LINKS+=-lfftw3
 RFLAGS=-O2 -ffixed-line-length-none -ffree-line-length-none -fexternal-blas
 DFLAGS=-g -fbacktrace -fbounds-check -ffree-line-length-none -fexternal-blas -DDEBUG
endif

ifeq ($(BUILD_ENV),gfortran_debian)
 # You should install the following packages
 # sudo apt install gfortran
 # sudo apt install libopenmpi-dev
 # sudo apt install liblapack-dev
 # sudo apt install m4
 # sudo apt install libfftw3-dev
 # sudo apt install libhdf5-dev
 FC=mpif90
 FLAGS=-fPIC
 CFLAGS=-fdefault-real-8
 RFLAGS=-O2 -ffixed-line-length-none -ffree-line-length-none -fexternal-blas
 DFLAGS=-g -fbacktrace -fbounds-check -ffree-line-length-none -fexternal-blas -DDEBUG
 LINKS=-llapack -lblas
 LIBS=-I/usr/include/hdf5/serial
 LINKS+=-L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_fortran -lhdf5 -lpthread -lz -lm
 LIBS+=-I/usr/include
 LINKS+=-lfftw3
endif

ifeq ($(BUILD_ENV),gfortran_arch)
 # configuration for Arch Linux
 # traps for numerical errors: https://stackoverflow.com/a/29827243
 FC=mpif90
 FLAGS=-fPIC
 CFLAGS=-fdefault-real-8 -fallow-argument-mismatch
 LINKS=-llapack -lblas
 LIBS=
 LINKS+=-lhdf5_fortran -lhdf5 -lpthread -lz -lm
 LINKS+=-lfftw3
 RFLAGS=-O3 -ffixed-line-length-none -ffree-line-length-none -fexternal-blas
 DFLAGS=-g -fbacktrace -fbounds-check -ffree-line-length-none -fexternal-blas -DDEBUG -ffpe-trap=invalid,zero,overflow,underflow
endif

ifeq ($(BUILD_ENV),gfortran_mac)
 # works on Ksenia's laptop
 FC=mpif90
 FLAGS=-fPIC
 CFLAGS=-fdefault-real-8
 LINKS=-L/usr/local/Cellar/lapack/3.8.0_1/lib -llapack -Wl,-rpath -Wl,/usr/local/Cellar/lapack/3.8.0_1/lib -lblas -L/usr/local/Cellar/openblas/0.3.7 -lgfortran -Wl,-rpath -Wl,/usr/local/Cellar/openblas/0.3.7/lib
 LIBS=-I/Applications/HDF_Group/HDF5/1.10.5/include/static
 LINKS+=-L/Applications/HDF_Group/HDF5/1.10.5/lib -lhdf5_f90cstub -lhdf5_fortran -lhdf5 -lpthread -lz -lm -Wl,-rpath -Wl,/Applications/HDF_Group/HDF5/1.10.5/lib
 LIBS+=-I/usr/local/Cellar/fftw/3.3.8/include
 LINKS+=-L/usr/local/Cellar/fftw/3.3.8/lib -lfftw3 -Wl,-rpath -Wl,/usr/local/Cellar/fftw/3.3.8/lib
 RFLAGS=-O2 -ffixed-line-length-none -ffree-line-length-none -fexternal-blas
 DFLAGS=-g -fbacktrace -fbounds-check -ffree-line-length-none -fexternal-blas -DDEBUG
endif

ifeq ($(BUILD_ENV),lff95)
 # LF95 SAL
 # Not checked
 CFLAGS=--dbl
 RFLAGS=--ap -O -I.
 DFLAGS=-Cpp -DDEBUG
 LINKS=-L$(LINKS_ROOT) -lnag -L$(LAPACKHOME) -llapack -L$(BLASHOME) -lblas
 LIBS=-I$(FFTWHOME)/include
 LINKS+=-L$(FFTWHOME)/lib -lfftw3
endif

ifeq ($(BUILD_ENV),intel_spc)
 FC=mpif90
 CFLAGS=-r8 -DIFORT
 RFLAGS=-O2 -ip -no-prec-div -xHost -fPIC
 DFLAGS=-traceback -D DEBUG -g
 LINKS=-L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
 LIBS=-I$(FFTW_DIR)/include
 LINKS+=-L$(FFTW_DIR)/lib -lfftw3
 LIBS+=-I$(HDF5_HOME)/include
 LINKS+=-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lpthread -lz -lm -Wl,-rpath -Wl,$(HDF5_HOME)/lib
endif

ifeq ($(BUILD_ENV),gfort_spc)
 FC=mpif90
 FLAGS=-fPIC
 CFLAGS=-fdefault-real-8
 RFLAGS=-O2 -ffixed-line-length-none -ffree-line-length-none -fexternal-blas
 DFLAGS=-g -fbacktrace -fbounds-check -ffree-line-length-none -fexternal-blas -DDEBUG
 LINKS=-L/lib -lmkl_rt -lpthread -lm -ldl -Wl,-rpath
 LINKS+=-L$(FFTW_DIR)/lib -lfftw3
 LINKS+=-L/usr/local/hdf5-1.8.18-gcc6.3/lib -lhdf5_fortran -lhdf5 -lz
 LIBS=-I$(FFTW_DIR)/include
 LIBS+=-I/usr/local/hdf5-1.8.18-gcc6.3/include
endif


ifeq ($(BUILD_ENV),intel_ipp)
 # tested on draco with the following modules:
 # intel/18.0.3 impi/2018.3 mkl/2018.3 hdf5-serial/1.8.21 fftw-mpi/3.3.8
 # and on cobra with the following modules:
 # intel/19.0.4 impi/2019.4 mkl/2019.4 hdf5-serial/1.8.21 fftw-mpi/3.3.8
 FC=mpiifort
 CFLAGS=-r8 -DIFORT
 RFLAGS=-O2 -ip -no-prec-div -xHost -fPIC
 DFLAGS=-traceback -D DEBUG
 LINKS=-L${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -lm -ldl -Wl,-rpath -Wl,${MKLROOT}/lib/intel64
 LIBS=-I$(HDF5_HOME)/include
 LINKS+=-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lpthread -lz -lm -Wl,-rpath -Wl,$(HDF5_HOME)/lib
 LIBS+=-I$(FFTW_HOME)/include
 LINKS+=-L$(FFTW_HOME)/lib -lfftw3 -Wl,-rpath -Wl,$(FFTW_HOME)/lib
endif

ifeq ($(BUILD_ENV),gfortran_ipp)
 # tested on draco with the following modules:
 # gcc/8 impi/2018.3 mkl/2018.3 hdf5-mpi/1.10.5 fftw-mpi/3.3.8
 FC=mpif90
 FLAGS=-fPIC
 CFLAGS=-fdefault-real-8
 RFLAGS=-O2 -fPIC -ffree-line-length-none
 DFLAGS=-g -fbacktrace -fbounds-check -DDEBUG -ffree-line-length-none
 LINKS=-L${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -lm -ldl -Wl,-rpath -Wl,${MKLROOT}/lib/intel64
 LIBS=-I$(HDF5_HOME)/include
 LINKS+=-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lpthread -lz -lm -Wl,-rpath -Wl,$(HDF5_HOME)/lib
 LIBS+=-I$(FFTW_HOME)/include
 LINKS+=-L$(FFTW_HOME)/lib -lfftw3 -Wl,-rpath -Wl,$(FFTW_HOME)/lib
endif

ifeq ($(BUILD_ENV),intel_gadi)
#module load intel-compiler/2019.3.199
#module load intel-mkl/2019.3.199
#module load openmpi/3.1.4
#module load fftw3-mkl/2019.3.199
#module load hdf5/1.10.5
 FC=mpif90
 CFLAGS=-r8 -DIFORT
 LINKS=-L${MKLROOT}/lib/intel64 -mkl=parallel -liomp5
 LIBS=-I$(HDF5_BASE)/include
 LINKS+=-L$(HDF5_BASE)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lpthread -lz -lm
 RFLAGS=-mcmodel=large -O3 -m64 -unroll0 -fno-alias -ip -traceback -fPIC
 DFLAGS=-check bounds -check format -check output_conversion -check pointers -check uninit -debug full -D DEBUG
endif

ifeq ($(BUILD_ENV),intel_stellar)
 # PPPL stellar cluster intel compiler
 # module use /home/caoxiang/module
 # module load spec
 FC=mpiifort
 CFLAGS=-r8 -DIFORT
 RFLAGS=-mcmodel=large -O3 -m64 -unroll0 -fno-alias -ip -traceback -fPIC
 DFLAGS=-O0 -g -traceback -check bounds -check format -check output_conversion -check pointers -check uninit -debug full -D DEBUG
 LIBS=-I${MKLROOT}/include/intel64/lp64  # MKL include
 LIBS+=-I$(HDF5DIR)/include # HDF5 include
 LINKS=-L$(HDF5DIR)/lib -lhdf5_fortran -lhdf5 # HDF5 link
 LIBS+=-I$(FFTW3DIR)/include # FFTW include
 LINKS+=-L$(FFTW3DIR)/lib -lfftw3 # FFTW link
 LINKS+= ${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl

endif

ifeq ($(OMP),yes)
 RFLAGS+=-DOPENMP -fopenmp
 DFLAGS+=-DOPENMP -fopenmp
 LINKS+=-lgomp
endif

###############################################################################################################################################################

 date:=$(shell date)
 text:=$(shell date +%F)

###############################################################################################################################################################

xspec: $(addsuffix _r.o,$(ALLFILES)) $(MACROS) Makefile
	$(FC) $(FLAGS) $(CFLAGS) $(RFLAGS) -o xspec $(addsuffix _r.o,$(ALLFILES)) $(LINKS)

dspec: $(addsuffix _d.o,$(ALLFILES)) $(MACROS) Makefile
	$(FC) $(FLAGS) $(CFLAGS) $(DFLAGS) -o dspec $(addsuffix _d.o,$(ALLFILES)) $(LINKS)

###############################################################################################################################################################
# inputlist needs special handling: expansion of DSCREENLIST and NSCREENLIST using awk

inputlist_r.o: %_r.o: inputlist.f90 $(MACROS)
	@awk -v allfiles='$(ALLFILES)' 'BEGIN{nfiles=split(allfiles,files," ")} \
	{if($$2=="DSCREENLIST") {for (i=1;i<=nfiles;i++) print "  LOGICAL :: W"files[i]" = .false. "}}\
	{if($$2=="NSCREENLIST") {for (i=1;i<=nfiles;i++) print "  W"files[i]" , &"}}\
	{print}' inputlist.f90 > mnputlist.f90
	m4 -P $(MACROS) mnputlist.f90 > inputlist_m.F90
	@rm -f mnputlist.f90
	$(FC) $(FLAGS) $(CFLAGS) $(RFLAGS) -o inputlist_r.o -c inputlist_m.F90 $(LIBS)
	@wc -l -L -w inputlist_m.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

inputlist_d.o: %_d.o: inputlist.f90 $(MACROS)
	@awk -v allfiles='$(ALLFILES)' 'BEGIN{nfiles=split(allfiles,files," ")} \
	{if($$2=="DSCREENLIST") {for (i=1;i<=nfiles;i++) print "  LOGICAL :: W"files[i]" = .false. "}}\
	{if($$2=="NSCREENLIST") {for (i=1;i<=nfiles;i++) print "  W"files[i]" , &"}}\
	{print}' inputlist.f90 > mnputlist.f90
	m4 -P $(MACROS) mnputlist.f90 > inputlist_m.F90
	@rm -f mnputlist.f90
	$(FC) $(FLAGS) $(CFLAGS) $(DFLAGS) -o inputlist_d.o -c inputlist_m.F90 $(LIBS)
	@wc -l -L -w inputlist_m.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

###############################################################################################################################################################
# global needs special handling: expansion of CPUVARIABLE, BSCREENLIST and WSCREENLIST using awk

global_r.o: %_r.o: inputlist_r.o global.f90 $(MACROS)
	@awk -v allfiles='$(ALLFILES)' 'BEGIN{nfiles=split(allfiles,files," ")} \
	{if($$2=="CPUVARIABLE") {for (i=1;i<=nfiles;i++) print "  REAL    :: T"files[i]" = 0.0, "files[i]"T = 0.0"}}\
	{if($$2=="BSCREENLIST") {for (i=1;i<=nfiles;i++) print "  LlBCAST(W"files[i]",1,0)"}}\
	{if($$2=="WSCREENLIST") {s="'"'"'" ; d="'"\\\""'" ; for (i=1;i<=nfiles;i++) print "  if( W"files[i]" ) write(iunit,"s"("d" W"files[i]" = "d"L1)"s")W"files[i]}}\
	{print}' global.f90 > mlobal.f90
	m4 -P $(MACROS) mlobal.f90 > global_m.F90
	@rm -f mlobal.f90
	$(FC) $(FLAGS) $(CFLAGS) $(RFLAGS) -o global_r.o -c global_m.F90 $(LIBS)
	@wc -l -L -w global_m.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

global_d.o: %_d.o: inputlist_d.o global.f90 $(MACROS)
	@awk -v allfiles='$(ALLFILES)' 'BEGIN{nfiles=split(allfiles,files," ")} \
	{if($$2=="CPUVARIABLE") {for (i=1;i<=nfiles;i++) print "  REAL    :: T"files[i]" = 0.0, "files[i]"T = 0.0"}}\
	{if($$2=="BSCREENLIST") {for (i=1;i<=nfiles;i++) print "  LlBCAST(W"files[i]",1,0)"}}\
	{if($$2=="WSCREENLIST") {s="'"'"'" ; d="'"\\\""'" ; for (i=1;i<=nfiles;i++) print "  if( W"files[i]" ) write(iunit,"s"("d" W"files[i]" = "d"L1)"s")W"files[i]}}\
	{print}' global.f90 > mlobal.f90
	m4 -P $(MACROS) mlobal.f90 > global_m.F90
	@rm -f mlobal.f90
	$(FC) $(FLAGS) $(CFLAGS) $(DFLAGS) -o global_d.o -c global_m.F90 $(LIBS)
	@wc -l -L -w global_m.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

###############################################################################################################################################################

%_r.o: %.f
	$(FC) $(FLAGS) $(RFLAGS) -o $*_r.o -c $*.f
	@wc -l -L -w $*.f | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

%_d.o: %.f
	$(FC) $(FLAGS) $(DFLAGS) -o $*_d.o -c $*.f
	@wc -l -L -w $*.f | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

###############################################################################################################################################################

$(PREPROC): %_m.F90: %.f90 $(MACROS)
	@awk -v file=$*.f90 '{ gsub("__LINE__", NR); gsub("__FILE__",file); print }' $*.f90 > $*_p.f90
	m4 -P $(MACROS) $*_p.f90 > $*_m.F90


$(ROBJS_IO): %_r.o: %_m.F90 $(addsuffix _r.o,$(BASEFILES)) $(MACROS)
	$(FC) $(FLAGS) $(CFLAGS) $(RFLAGS) -o $*_r.o -c $*_m.F90 $(LIBS)
	@wc -l -L -w $*_m.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

$(DOBJS_IO): %_d.o: %_m.F90 $(addsuffix _d.o,$(BASEFILES)) $(MACROS)
	$(FC) $(FLAGS) $(CFLAGS) $(DFLAGS) -o $*_d.o -c $*_m.F90 $(LIBS)
	@wc -l -L -w $*_m.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''


$(ROBJS): %_r.o: %_m.F90 $(addsuffix _r.o,$(BASEFILES)) $(addsuffix _r.o,$(IOFILES)) $(MACROS)
	$(FC) $(FLAGS) $(CFLAGS) $(RFLAGS) -o $*_r.o -c $*_m.F90 $(LIBS)
	@wc -l -L -w $*_m.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

$(DOBJS): %_d.o: %_m.F90 $(addsuffix _d.o,$(BASEFILES)) $(addsuffix _d.o,$(IOFILES)) $(MACROS)
	$(FC) $(FLAGS) $(CFLAGS) $(DFLAGS) -o $*_d.o -c $*_m.F90 $(LIBS)
	@wc -l -L -w $*_m.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

###############################################################################################################################################################

xspech_r.o: xspech.f90 $(addsuffix _r.o,$(ALLSPEC)) $(MACROS)
	@awk -v date='$(date)' -v pwd='$(PWD)' -v macros='$(MACROS)' -v fc='$(FC)' -v flags='$(FLAGS) $(CFLAGS) $(RFLAGS)' -v allfiles='$(ALLFILES)' \
	'BEGIN{nfiles=split(allfiles,files," ")} \
	{if($$2=="COMPILATION") {print "    write(ounit,*)\"      :  compiled  : date    = "date" ; \"" ; \
	                         print "    write(ounit,*)\"      :            : srcdir  = "pwd" ; \"" ; \
	                         print "    write(ounit,*)\"      :            : macros  = "macros" ; \"" ; \
	                         print "    write(ounit,*)\"      :            : fc      = "fc" ; \"" ; \
	                         print "    write(ounit,*)\"      :            : flags   = "flags" ; \"" }} \
	 {if($$2=="SUMTIME") {for (i=1;i<=nfiles;i++) print "   SUMTIME("files[i]")"}}\
	 {if($$2=="PRTTIME") {for (i=1;i<=nfiles;i++) print "   PRTTIME("files[i]")"}}\
	 {print}' xspech.f90 > mspech.f90
	m4 -P $(MACROS) mspech.f90 > xspech_m.F90
	@rm -f mspech.f90
	$(FC) $(FLAGS) $(CFLAGS) $(RFLAGS) -o xspech_r.o -c xspech_m.F90 $(LIBS)
	@wc -l -L -w xspech_m.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

xspech_d.o: xspech.f90 $(addsuffix _d.o,$(ALLSPEC)) $(MACROS)
	@awk -v date='$(date)' -v pwd='$(PWD)' -v macros='$(MACROS)' -v fc='$(FC)' -v flags='$(FLAGS) $(CFLAGS) $(DFLAGS)' -v allfiles='$(ALLFILES)' \
	'BEGIN{nfiles=split(allfiles,files," ")} \
	{if($$2=="COMPILATION") {print "    write(ounit,*)\"      :  compiled  : date    = "date" ; \"" ; \
	                         print "    write(ounit,*)\"      :            : srcdir  = "pwd" ; \"" ; \
	                         print "    write(ounit,*)\"      :            : macros  = "macros" ; \"" ; \
	                         print "    write(ounit,*)\"      :            : fc      = "fc" ; \"" ; \
	                         print "    write(ounit,*)\"      :            : flags   = "flags" ; \"" }} \
	 {if($$2=="SUMTIME") {for (i=1;i<=nfiles;i++) print "   SUMTIME("files[i]")"}}\
	 {if($$2=="PRTTIME") {for (i=1;i<=nfiles;i++) print "   PRTTIME("files[i]")"}}\
	 {print}' xspech.f90 > mspech.f90
	m4 -P $(MACROS) mspech.f90 > xspech_m.F90
	@rm -f mspech.f90
	$(FC) $(FLAGS) $(CFLAGS) $(DFLAGS) -o xspech_d.o -c xspech_m.F90 $(LIBS)
	@wc -l -L -w xspech_m.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

###############################################################################################################################################################

clean:
	rm -f *.o *.mod *_p.f90 *_m.F90 .*.h *.pdf *.dvi *.out *.bbl *.toc .*.date ; rm -rf ./docs/

###############################################################################################################################################################

show_makeflags:
	@echo " "
	@echo " "
	@echo "----------------------------------------- "
	@echo " "
	@echo " "
	@echo "Fortran preprocessing"
	@echo "============================"
	@echo "$(MACROS) $(ALLFILES)"
	@echo " "
	@echo "$(BUILD_ENV) compile"
	@echo "============================="
	@echo "$(FC) $(FLAGS) $(CFLAGS) $(RFLAGS)"
	@echo " "
	@echo "$(FC) $(FLAGS) $(CFLAGS) $(DFLAGS)"
	@echo " "
	@echo "Include libraries"
	@echo "============================="
	@echo "$(LIBS)"
	@echo " "
	@echo "LINKING FLAGS"
	@echo "============================="
	@echo "$(LINKS)"

###############################################################################################################################################################
