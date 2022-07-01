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

%_r.o: src/%.f
	$(FC) $(FLAGS) $(RFLAGS) -o $*_r.o -c src/$*.f
	@wc -l -L -w src/$*.f | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

%_d.o: src/%.f
	$(FC) $(FLAGS) $(DFLAGS) -o $*_d.o -c src/$*.f
	@wc -l -L -w src/$*.f | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

###############################################################################################################################################################

mod_kinds_r.o : %_r.o: src/mod_kinds.F90
	$(FC) $(FLAGS) $(CFLAGS) $(RFLAGS) -o mod_kinds_r.o -c src/mod_kinds.F90 $(LIBS)
	@wc -l -L -w src/mod_kinds.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

mod_kinds_d.o : %_d.o: src/mod_kinds.F90
	$(FC) $(FLAGS) $(CFLAGS) $(DFLAGS) -o mod_kinds_d.o -c src/mod_kinds.F90 $(LIBS)
	@wc -l -L -w src/mod_kinds.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

inputlist_r.o: %_r.o: mod_kinds_r.o src/inputlist.F90
	$(FC) $(FLAGS) $(CFLAGS) $(RFLAGS) -o inputlist_r.o -c src/inputlist.F90 $(LIBS)
	@wc -l -L -w src/inputlist.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

inputlist_d.o: %_d.o: mod_kinds_d.o src/inputlist.F90
	$(FC) $(FLAGS) $(CFLAGS) $(DFLAGS) -o inputlist_d.o -c src/inputlist.F90 $(LIBS)
	@wc -l -L -w src/inputlist.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

###############################################################################################################################################################
# global needs special handling: expansion of CPUVARIABLE, BSCREENLIST and WSCREENLIST using awk (not anymore !!!)

global_r.o: %_r.o: inputlist_r.o src/global.F90
	$(FC) $(FLAGS) $(CFLAGS) $(RFLAGS) -o global_r.o -c src/global.F90 $(LIBS)
	@wc -l -L -w src/global.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

global_d.o: %_d.o: inputlist_d.o src/global.F90
	$(FC) $(FLAGS) $(CFLAGS) $(DFLAGS) -o global_d.o -c src/global.F90 $(LIBS)
	@wc -l -L -w src/global.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

###############################################################################################################################################################

$(ROBJS_IO): %_r.o: src/%.F90 $(addsuffix _r.o,$(BASEFILES))
	$(FC) $(FLAGS) $(CFLAGS) $(RFLAGS) -o $*_r.o -c src/$*.F90 $(LIBS)
	@wc -l -L -w src/$*.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

$(DOBJS_IO): %_d.o: src/%.F90 $(addsuffix _d.o,$(BASEFILES))
	$(FC) $(FLAGS) $(CFLAGS) $(DFLAGS) -o $*_d.o -c src/$*.F90 $(LIBS)
	@wc -l -L -w src/$*.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''


$(ROBJS): %_r.o: src/%.F90 $(addsuffix _r.o,$(BASEFILES)) $(addsuffix _r.o,$(IOFILES))
	$(FC) $(FLAGS) $(CFLAGS) $(RFLAGS) -o $*_r.o -c src/$*.F90 $(LIBS)
	@wc -l -L -w src/$*.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

$(DOBJS): %_d.o: src/%.F90 $(addsuffix _d.o,$(BASEFILES)) $(addsuffix _d.o,$(IOFILES))
	$(FC) $(FLAGS) $(CFLAGS) $(DFLAGS) -o $*_d.o -c src/$*.F90 $(LIBS)
	@wc -l -L -w src/$*.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

###############################################################################################################################################################

xspech_r.o: src/xspech.F90 global_r.o sphdf5_r.o $(addsuffix _r.o,$(files))
	@awk -v date='$(date)' -v pwd='$(PWD)' -v fc='$(FC)' -v flags='$(FLAGS) $(CFLAGS) $(RFLAGS)' -v allfiles='$(ALLFILES)' \
	'BEGIN{nfiles=split(allfiles,files," ")} \
	{if($$2=="COMPILATION") {print "    write(ounit,*)\"      :  compiled  : date    = "date" ; \"" ; \
	                         print "    write(ounit,*)\"      :            : srcdir  = "pwd" ; \"" ; \
	                         print "    write(ounit,*)\"      :            : fc      = "fc" ; \"" ; \
	                         print "    write(ounit,*)\"      :            : flags   = "flags" ; \"" }} \
	 {print}' src/xspech.F90 > xspech_m.F90
	$(FC) $(FLAGS) $(CFLAGS) $(RFLAGS) -o xspech_r.o -c xspech_m.F90 $(LIBS)
	@wc -l -L -w xspech_m.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

xspech_d.o: src/xspech.F90 global_d.o sphdf5_d.o $(addsuffix _d.o,$(files))
	@awk -v date='$(date)' -v pwd='$(PWD)' -v fc='$(FC)' -v flags='$(FLAGS) $(CFLAGS) $(DFLAGS)' -v allfiles='$(ALLFILES)' \
	'BEGIN{nfiles=split(allfiles,files," ")} \
	{if($$2=="COMPILATION") {print "    write(ounit,*)\"      :  compiled  : date    = "date" ; \"" ; \
	                         print "    write(ounit,*)\"      :            : srcdir  = "pwd" ; \"" ; \
	                         print "    write(ounit,*)\"      :            : fc      = "fc" ; \"" ; \
	                         print "    write(ounit,*)\"      :            : flags   = "flags" ; \"" }} \
	 {print}' src/xspech.F90 > xspech_m.F90
	$(FC) $(FLAGS) $(CFLAGS) $(DFLAGS) -o xspech_d.o -c xspech_m.F90 $(LIBS)
	@wc -l -L -w xspech_m.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

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
