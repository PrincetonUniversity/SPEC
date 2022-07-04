#!/bin/sh

############################################
# This is the "classic" Makefile for SPEC. #
############################################

# The files that make up SPEC are registered in SPECfile
# and machine-specific build options (linker flags etc.)
# are defined there as well.
include SPECfile

ROBJS=$(SPECFILES:=_r.o)
DOBJS=$(SPECFILES:=_d.o)

ROBJS_LEVEL_0=$(LEVEL_0:=_r.o)
DOBJS_LEVEL_0=$(LEVEL_0:=_d.o)

ROBJS_LEVEL_1=$(LEVEL_1:=_r.o)
DOBJS_LEVEL_1=$(LEVEL_1:=_d.o)

ROBJS_LEVEL_2=$(LEVEL_2:=_r.o)
DOBJS_LEVEL_2=$(LEVEL_2:=_d.o)

ROBJS_LEVEL_3=$(LEVEL_3:=_r.o)
DOBJS_LEVEL_3=$(LEVEL_3:=_d.o)

ROBJS_LEVEL_4=$(LEVEL_4:=_r.o)
DOBJS_LEVEL_4=$(LEVEL_4:=_d.o)

ROBJS_BASE=$(BASEFILES:=_r.o)
DOBJS_BASE=$(BASEFILES:=_d.o)

ROBJS_IO=$(IOFILES:=_r.o)
DOBJS_IO=$(IOFILES:=_d.o)

###############################################################################################################################################################

date:=$(shell date)
text:=$(shell date +%F)

###############################################################################################################################################################

xspec: $(addsuffix _r.o,$(ALLFILES)) Makefile
	$(FC) $(FLAGS) $(CFLAGS) $(RFLAGS) -o xspec $(addsuffix _r.o,$(ALLFILES)) $(LINKS)

dspec: $(addsuffix _d.o,$(ALLFILES)) Makefile
	$(FC) $(FLAGS) $(CFLAGS) $(DFLAGS) -o dspec $(addsuffix _d.o,$(ALLFILES)) $(LINKS)

###############################################################################################################################################################

%_r.o: src/contrib/%.f
	$(FC) $(FLAGS) $(RFLAGS) -o $*_r.o -c src/contrib/$*.f

%_d.o: src/contrib/%.f
	$(FC) $(FLAGS) $(DFLAGS) -o $*_d.o -c src/contrib/$*.f

###############################################################################################################################################################

$(ROBJS_LEVEL_0): %_r.o: src/%.F90
	$(FC) $(FLAGS) $(CFLAGS) $(RFLAGS) -o $*_r.o -c src/$*.F90 $(LIBS)

$(DOBJS_LEVEL_0): %_d.o: src/%.F90
	$(FC) $(FLAGS) $(CFLAGS) $(DFLAGS) -o $*_d.o -c src/$*.F90 $(LIBS)

$(ROBJS_LEVEL_1): %_r.o: src/%.F90 $(addsuffix _r.o,$(LEVEL_0))
	$(FC) $(FLAGS) $(CFLAGS) $(RFLAGS) -o $*_r.o -c src/$*.F90 $(LIBS)

$(DOBJS_LEVEL_1): %_d.o: src/%.F90 $(addsuffix _d.o,$(LEVEL_0))
	$(FC) $(FLAGS) $(CFLAGS) $(DFLAGS) -o $*_d.o -c src/$*.F90 $(LIBS)

$(ROBJS_LEVEL_2): %_r.o: src/%.F90 mod_kinds_r.o $(addsuffix _r.o,$(LEVEL_0)) $(addsuffix _r.o,$(LEVEL_1))
	$(FC) $(FLAGS) $(CFLAGS) $(RFLAGS) -o $*_r.o -c src/$*.F90 $(LIBS)

$(DOBJS_LEVEL_2): %_d.o: src/%.F90 mod_kinds_d.o $(addsuffix _d.o,$(LEVEL_0)) $(addsuffix _d.o,$(LEVEL_1))
	$(FC) $(FLAGS) $(CFLAGS) $(DFLAGS) -o $*_d.o -c src/$*.F90 $(LIBS)

$(ROBJS_LEVEL_3): %_r.o: src/%.F90 mod_kinds_r.o $(addsuffix _r.o,$(LEVEL_0)) $(addsuffix _r.o,$(LEVEL_1)) $(addsuffix _r.o,$(LEVEL_2))
	$(FC) $(FLAGS) $(CFLAGS) $(RFLAGS) -o $*_r.o -c src/$*.F90 $(LIBS)

$(DOBJS_LEVEL_3): %_d.o: src/%.F90 mod_kinds_d.o $(addsuffix _d.o,$(LEVEL_0)) $(addsuffix _d.o,$(LEVEL_1)) $(addsuffix _d.o,$(LEVEL_2))
	$(FC) $(FLAGS) $(CFLAGS) $(DFLAGS) -o $*_d.o -c src/$*.F90 $(LIBS)

$(ROBJS_LEVEL_4): %_r.o: src/%.F90 mod_kinds_r.o $(addsuffix _r.o,$(LEVEL_0)) $(addsuffix _r.o,$(LEVEL_1)) $(addsuffix _r.o,$(LEVEL_2)) $(addsuffix _r.o,$(LEVEL_3))
	$(FC) $(FLAGS) $(CFLAGS) $(RFLAGS) -o $*_r.o -c src/$*.F90 $(LIBS)

$(DOBJS_LEVEL_4): %_d.o: src/%.F90 mod_kinds_d.o $(addsuffix _d.o,$(LEVEL_0)) $(addsuffix _d.o,$(LEVEL_1)) $(addsuffix _d.o,$(LEVEL_2)) $(addsuffix _d.o,$(LEVEL_3))
	$(FC) $(FLAGS) $(CFLAGS) $(DFLAGS) -o $*_d.o -c src/$*.F90 $(LIBS)

###############################################################################################################################################################

$(ROBJS_IO): %_r.o: src/%.F90 $(addsuffix _r.o,$(BASEFILES))
	$(FC) $(FLAGS) $(CFLAGS) $(RFLAGS) -o $*_r.o -c src/$*.F90 $(LIBS)

$(DOBJS_IO): %_d.o: src/%.F90 $(addsuffix _d.o,$(BASEFILES))
	$(FC) $(FLAGS) $(CFLAGS) $(DFLAGS) -o $*_d.o -c src/$*.F90 $(LIBS)

$(ROBJS): %_r.o: src/%.F90 $(addsuffix _r.o,$(BASEFILES)) $(addsuffix _r.o,$(IOFILES))
	$(FC) $(FLAGS) $(CFLAGS) $(RFLAGS) -o $*_r.o -c src/$*.F90 $(LIBS)

$(DOBJS): %_d.o: src/%.F90 $(addsuffix _d.o,$(BASEFILES)) $(addsuffix _d.o,$(IOFILES))
	$(FC) $(FLAGS) $(CFLAGS) $(DFLAGS) -o $*_d.o -c src/$*.F90 $(LIBS)

###############################################################################################################################################################

xspech_r.o: src/xspech.F90 global_r.o sphdf5_r.o $(addsuffix _r.o,$(files))
	$(FC) $(FLAGS) $(CFLAGS) $(RFLAGS) -o xspech_r.o -c src/xspech.F90 $(LIBS)

xspech_d.o: src/xspech.F90 global_d.o sphdf5_d.o $(addsuffix _d.o,$(files))
	$(FC) $(FLAGS) $(CFLAGS) $(DFLAGS) -o xspech_d.o -c src/xspech.F90 $(LIBS)

###############################################################################################################################################################

clean:
	rm -f *.o *.mod *.pdf *.dvi *.out *.bbl *.toc .*.date *_m.F90
	rm -rf ./docs/html ./docs/latex

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
