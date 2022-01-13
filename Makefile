#!/bin/sh

############################################
# This is the "classic" Makefile for SPEC. #
############################################

# The files that make up SPEC are registered in SPECfile
# and machine-specific build options (linker flags etc.)
# are defined there as well.
include SPECfile

# files to be preprocessed by m4
PREPROC=$(ALLSPEC:=_m.F90)

ROBJS=$(SPECFILES:=_r.o)
DOBJS=$(SPECFILES:=_d.o)

ROBJS_IO=$(IOFILES:=_r.o)
DOBJS_IO=$(IOFILES:=_d.o)

###############################################################################################################################################################

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
# inputlist needs special handling: expansion of DSCREENLIST and NSCREENLIST using awk (not anymore !!!)

inputlist_r.o: %_r.o: src/inputlist.f90 $(MACROS)
	m4 -P $(MACROS) src/inputlist.f90 > inputlist_m.F90
	$(FC) $(FLAGS) $(CFLAGS) $(RFLAGS) -o inputlist_r.o -c inputlist_m.F90 $(LIBS)
	@wc -l -L -w inputlist_m.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

inputlist_d.o: %_d.o: src/inputlist.f90 $(MACROS)
	m4 -P $(MACROS) src/inputlist.f90 > inputlist_m.F90
	$(FC) $(FLAGS) $(CFLAGS) $(DFLAGS) -o inputlist_d.o -c inputlist_m.F90 $(LIBS)
	@wc -l -L -w inputlist_m.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

###############################################################################################################################################################
# global needs special handling: expansion of CPUVARIABLE, BSCREENLIST and WSCREENLIST using awk (not anymore !!!)

global_r.o: %_r.o: inputlist_r.o src/global.f90 $(MACROS)
	m4 -P $(MACROS) src/global.f90 > global_m.F90
	$(FC) $(FLAGS) $(CFLAGS) $(RFLAGS) -o global_r.o -c global_m.F90 $(LIBS)
	@wc -l -L -w global_m.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

global_d.o: %_d.o: inputlist_d.o src/global.f90 $(MACROS)
	m4 -P $(MACROS) src/global.f90 > global_m.F90
	$(FC) $(FLAGS) $(CFLAGS) $(DFLAGS) -o global_d.o -c global_m.F90 $(LIBS)
	@wc -l -L -w global_m.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

###############################################################################################################################################################

%_r.o: src/%.f
	$(FC) $(FLAGS) $(RFLAGS) -o $*_r.o -c src/$*.f
	@wc -l -L -w src/$*.f | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

%_d.o: src/%.f
	$(FC) $(FLAGS) $(DFLAGS) -o $*_d.o -c src/$*.f
	@wc -l -L -w src/$*.f | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

###############################################################################################################################################################

$(PREPROC): %_m.F90: src/%.f90 $(MACROS)
	@awk -v file=$*.f90 '{ gsub("__LINE__", NR); gsub("__FILE__",file); print }' src/$*.f90 > $*_p.f90
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

xspech_r.o: src/xspech.f90 global_r.o sphdf5_r.o $(addsuffix _r.o,$(files)) $(MACROS)
	@awk -v date='$(date)' -v pwd='$(PWD)' -v macros='$(MACROS)' -v fc='$(FC)' -v flags='$(FLAGS) $(CFLAGS) $(RFLAGS)' -v allfiles='$(ALLFILES)' \
	'BEGIN{nfiles=split(allfiles,files," ")} \
	{if($$2=="COMPILATION") {print "    write(ounit,*)\"      :  compiled  : date    = "date" ; \"" ; \
	                         print "    write(ounit,*)\"      :            : srcdir  = "pwd" ; \"" ; \
	                         print "    write(ounit,*)\"      :            : macros  = "macros" ; \"" ; \
	                         print "    write(ounit,*)\"      :            : fc      = "fc" ; \"" ; \
	                         print "    write(ounit,*)\"      :            : flags   = "flags" ; \"" }} \
	 {if($$2=="SUMTIME") {for (i=1;i<=nfiles;i++) print "   SUMTIME("files[i]")"}}\
	 {if($$2=="PRTTIME") {for (i=1;i<=nfiles;i++) print "   PRTTIME("files[i]")"}}\
	 {print}' src/xspech.f90 > mspech.f90
	m4 -P $(MACROS) mspech.f90 > xspech_m.F90
	@rm -f mspech.f90
	$(FC) $(FLAGS) $(CFLAGS) $(RFLAGS) -o xspech_r.o -c xspech_m.F90 $(LIBS)
	@wc -l -L -w xspech_m.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

xspech_d.o: src/xspech.f90 global_d.o sphdf5_d.o $(addsuffix _d.o,$(files)) $(MACROS)
	@awk -v date='$(date)' -v pwd='$(PWD)' -v macros='$(MACROS)' -v fc='$(FC)' -v flags='$(FLAGS) $(CFLAGS) $(DFLAGS)' -v allfiles='$(ALLFILES)' \
	'BEGIN{nfiles=split(allfiles,files," ")} \
	{if($$2=="COMPILATION") {print "    write(ounit,*)\"      :  compiled  : date    = "date" ; \"" ; \
	                         print "    write(ounit,*)\"      :            : srcdir  = "pwd" ; \"" ; \
	                         print "    write(ounit,*)\"      :            : macros  = "macros" ; \"" ; \
	                         print "    write(ounit,*)\"      :            : fc      = "fc" ; \"" ; \
	                         print "    write(ounit,*)\"      :            : flags   = "flags" ; \"" }} \
	 {if($$2=="SUMTIME") {for (i=1;i<=nfiles;i++) print "   SUMTIME("files[i]")"}}\
	 {if($$2=="PRTTIME") {for (i=1;i<=nfiles;i++) print "   PRTTIME("files[i]")"}}\
	 {print}' src/xspech.f90 > mspech.f90
	m4 -P $(MACROS) mspech.f90 > xspech_m.F90
	@rm -f mspech.f90
	$(FC) $(FLAGS) $(CFLAGS) $(DFLAGS) -o xspech_d.o -c xspech_m.F90 $(LIBS)
	@wc -l -L -w xspech_m.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

###############################################################################################################################################################

clean:
	rm -f *.o *.mod *_p.f90 *_m.F90 .*.h *.pdf *.dvi *.out *.bbl *.toc .*.date ; rm -rf ./docs/html ./docs/latex

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
