#!/bin/sh

###############################################################################################################################################################

 afiles=manual rzaxis packxi volume coords basefn
 bfiles=metrix ma00aa matrix spsmat spsint mp00ac ma02aa packab tr00ab curent df00ab lforce intghs mtrxhs
#cfiles=bc00aa fc02aa jk03aa pc00aa pc00ab
 cfiles=brcast dforce newton 
 dfiles=casing bnorml 
 efiles=jo00aa pp00aa pp00ab bfield stzxyz sc00aa
 ffiles=hesian ra00aa numrec
 sfiles=dcuhre minpack iqpack rksuite i1mach d1mach ilut iters # below assumes the .f files are double precision; the CFLAGS = -r8 option is not required;

###############################################################################################################################################################

 SPECFILES=$(afiles) $(bfiles) $(cfiles) $(dfiles) $(efiles) $(ffiles)
 ALLFILES=global $(SPECFILES) $(sfiles) xspech hdfint preset
#F77FILES=$(sfiles:=.f)
 F90FILES=$(SPECFILES:=.F90)
 HFILES=global preset $(SPECFILES) hdfint xspech

 ROBJS=$(SPECFILES:=_r.o)
 DOBJS=$(SPECFILES:=_d.o)
 
###############################################################################################################################################################
 
 MACROS=macros
 
 CC=intel
 # if want to use gfortran; make CC=gfortran; otherwise using Intel
 FC=mpif90
 OMP=yes
 # to enable OpenMP acceleration within volume, set OMP=yes, otherwise set OMP=no
 
 # Intel Defaults
 # At PPPL
 # module load intel/2017
 # module load hdf5
 # module load openmpi
 # module load fftw
 # module load lapack/3.5.0rhel6
 CFLAGS=-r8
 RFLAGS=-mcmodel=large -O3 -m64 -unroll0 -fno-alias -ip -traceback
 DFLAGS=-check bounds -check format -check output_conversion -check pointers -check uninit -debug full -D DEBUG
 #Note: on the PPPL clusters, use module lapack/3.5.0rhel6 only
 NAG=-L$(LAPACKHOME) -llapack -lblas -L$(BLASHOME) -lgfortran #-lrefblas -lgfortran
 NETCDF=-L$(NETCDFHOME)/lib -lnetcdf
 HDF5compile=-I$(HDF5_HOME)/include
 HDF5link=-L$(HDF5_HOME)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lpthread -lz -lm
 FFTWcompile=-I$(FFTWHOME)/include
 FFTWlink=-L$(FFTWHOME)/lib -lfftw3

ifeq ($(CC),gfortran)
 # Not checked
 CFLAGS=-fdefault-real-8
 NAG=-L$(ACML_HOME)/lib -lacml
 NETCDF=-L$(NETCDFHOME)/lib -lnetcdf
 HDF5compile=-I$(HDF5_HOME)/include
 HDF5link=-L$(HDF5_HOME)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lpthread -lz -lm
 FFTWcompile=-I$(FFTWHOME)/include
 FFTWlink=-L$(FFTWHOME)/lib -lfftw3
 RFLAGS=-O2 -ffixed-line-length-none -ffree-line-length-none -fexternal-blas
 DFLAGS=-g -fbacktrace -fbounds-check -ffree-line-length-none -fexternal-blas -DDEBUG
endif

ifeq ($(CC),gfortran_ubuntu)
 # You should install the following packages
 # sudo apt install gfortran
 # sudo apt install libopenmpi-dev
 # sudo apt install liblapack-dev
 # sudo apt install m4
 # sudo apt install libnetcdf-dev
 # sudo apt install libfftw3-dev
 # sudo apt install libhdf5-openmpi-dev
 CFLAGS=-fdefault-real-8
 NAG=-llapack -lblas
 NETCDF=-lnetcdf
 HDF5compile=-I/usr/include/hdf5/openmpi
 HDF5link=-L/usr/lib/x86_64-linux-gnu/hdf5/openmpi -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lpthread -lz -lm
 FFTWcompile=-I/usr/include
 FFTWlink=-lfftw3
 RFLAGS=-O2 -ffixed-line-length-none -ffree-line-length-none -fexternal-blas
 DFLAGS=-g -fbacktrace -fbounds-check -ffree-line-length-none -fexternal-blas -DDEBUG
endif 

ifeq ($(CC),lff95)
 # LF95 SAL
 # Not checked
 CFLAGS=--dbl
 RFLAGS=--ap -O -I.
 DFLAGS=-Cpp -DDEBUG
 NAG=-L$(NAG_ROOT) -lnag -L$(LAPACKHOME) -llapack -L$(BLASHOME) -lblas
 FFTWcompile=-I$(FFTWHOME)/include
 FFTWlink=-L$(FFTWHOME)/lib -lfftw3
endif


ifeq ($(CC),intel_spc)
 CFLAGS=-r8
 RFLAGS=-O2 -ip -no-prec-div -xHost -fPIC
 DFLAGS=-traceback -D DEBUG
 NAG=-L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
 NETCDF=-L$(NETCDF_HOME)/lib -lnetcdf
 FFTWcompile=-I$(FFTW_DIR)/include
 FFTWlink=-L$(FFTW_DIR)/lib -lfftw3
endif


ifeq ($(CC),intel_ipp)
 CFLAGS=-r8
 RFLAGS=-O2 -ip -no-prec-div -xHost -fPIC
 DFLAGS=-traceback -D DEBUG
# NAG=-L$(NAGFLIB_HOME)/lib -lnag_nag 
 NAG=-L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl 
 NETCDF=-L$(NETCDF_HOME)/lib -lnetcdf
 FFTWcompile=-I$(FFTW_DIR)/include
 FFTWlink=-L$(FFTW_DIR)/lib -lfftw3
endif

ifeq ($(CC),gfortran_ipp)
 CFLAGS=-fdefault-real-8
 RFLAGS=-O2 -fPIC -ffree-line-length-none
 DFLAGS=-g -fbacktrace -fbounds-check -DDEBUG -ffree-line-length-none
# NAG=-L$(NAGFLIB_HOME)/lib -lnag_nag 
 NAG=-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl 
 FFTWcompile=-I$(FFTW_DIR)/include
 FFTWlink=-L$(FFTW_DIR)/lib -lfftw3
endif

ifeq ($(CC),intel_raijin)
 # One needs to load the following modules
 # module load intel-fc/2018.1.163
 # module load intel-cc/2018.1.163
 # module load intel-mkl/2018.1.163
 # module load openmpi
 # module load fftw3-mkl/2018.1.163
 # module load netcdf
 # module load hdf5
 CFLAGS=-r8
 NAG=-L${MKLROOT}/lib/intel64 -mkl=parallel -liomp5 -lpthread -lm -ldl
 NETCDF=
 HDF5compile=-I$(HDF5_BASE)/include
 HDF5link=-L$(HDF5_BASE)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lpthread -lz -lm
 FFTWcompile=
 FFTWlink=
 RFLAGS=-mcmodel=large -O3 -m64 -unroll0 -fno-alias -ip -traceback -fPIC
 DFLAGS=-check bounds -check format -check output_conversion -check pointers -check uninit -debug full -D DEBUG
endif

ifeq ($(OMP),yes)
 RFLAGS+=-DOPENMP -fopenmp
endif

###############################################################################################################################################################

 WEBDIR=$(HOME)/w3_html

###############################################################################################################################################################

 date:=$(shell date)
 text:=$(shell date +%F)

###############################################################################################################################################################

xspec: $(addsuffix _r.o,$(ALLFILES)) $(MACROS) Makefile
	$(FC) $(FLAGS) $(CFLAGS) $(RFLAGS) -o xspec $(addsuffix _r.o,$(ALLFILES)) $(NAG) $(HDF5compile) $(HDF5link) $(NETCDF) $(FFTWlink)

dspec: $(addsuffix _d.o,$(ALLFILES)) $(MACROS) Makefile
	$(FC) $(FLAGS) $(CFLAGS) $(DFLAGS) -o dspec $(addsuffix _d.o,$(ALLFILES)) $(NAG) $(HDF5compile) $(HDF5link) $(NETCDF) $(FFTWlink)

###############################################################################################################################################################

global_r.o: %_r.o: global.h $(MACROS) Makefile
	@awk -v allfiles='$(ALLFILES)' 'BEGIN{nfiles=split(allfiles,files," ")} \
	{if($$2=="CPUVARIABLE") {for (i=1;i<=nfiles;i++) print "  REAL    :: T"files[i]" = 0.0, "files[i]"T = 0.0"}}\
	{if($$2=="DSCREENLIST") {for (i=1;i<=nfiles;i++) print "  LOGICAL :: W"files[i]" = .false. "}}\
	{if($$2=="NSCREENLIST") {for (i=1;i<=nfiles;i++) print "  W"files[i]" , &"}}\
	{if($$2=="BSCREENLIST") {for (i=1;i<=nfiles;i++) print "  LlBCAST(W"files[i]",1,0)"}}\
	{if($$2=="WSCREENLIST") {s="'"'"'" ; d="'"\\\""'" ; for (i=1;i<=nfiles;i++) print "  if( W"files[i]" ) write(iunit,"s"("d" W"files[i]" = "d"L1)"s")W"files[i]}}\
	{print}' global.h > mlobal.h
	m4 -P $(MACROS) mlobal.h > global.F90
	@rm -f mlobal.h
	$(FC) $(FLAGS) $(CFLAGS) $(RFLAGS) -o global_r.o -c global.F90 $(FFTWcompile)
	@wc -l -L -w global.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

global_d.o: %_d.o: global.h $(MACROS) Makefile
	@awk -v allfiles='$(ALLFILES)' 'BEGIN{nfiles=split(allfiles,files," ")} \
	{if($$2=="CPUVARIABLE") {for (i=1;i<=nfiles;i++) print "  REAL    :: T"files[i]" = 0.0, "files[i]"T = 0.0"}}\
	{if($$2=="DSCREENLIST") {for (i=1;i<=nfiles;i++) print "  LOGICAL :: W"files[i]" = .false. "}}\
	{if($$2=="NSCREENLIST") {for (i=1;i<=nfiles;i++) print "  W"files[i]" , &"}}\
	{if($$2=="BSCREENLIST") {for (i=1;i<=nfiles;i++) print "  LlBCAST(W"files[i]",1,0)"}}\
	{if($$2=="WSCREENLIST") {s="'"'"'" ; d="'"\\\""'" ; for (i=1;i<=nfiles;i++) print "  if( W"files[i]" ) write(iunit,"s"("d" W"files[i]" = "d"L1)"s")W"files[i]}}\
	{print}' global.h > mlobal.h
	m4 -P $(MACROS) mlobal.h > global.F90
	@rm -f mlobal.h
	$(FC) $(FLAGS) $(CFLAGS) $(DFLAGS) -o global_d.o -c global.F90 $(FFTWcompile)
	@wc -l -L -w global.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

###############################################################################################################################################################

hdfint_r.o: hdfint.h global_r.o $(MACROS) Makefile
	m4 -P $(MACROS) hdfint.h > $*.F90
	$(FC) $(FLAGS) $(CFLAGS) $(RFLAGS) -o hdfint_r.o -c $*.F90 $(HDF5compile)
	@wc -l -L -w hdfint_r.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

hdfint_d.o: hdfint.h global_d.o $(MACROS) Makefile
	m4 -P $(MACROS) hdfint.h > $*.F90
	$(FC) $(FLAGS) $(CFLAGS) $(DFLAGS) -o hdfint_d.o -c $*.F90 $(HDF5compile)
	@wc -l -L -w hdfint_d.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

###############################################################################################################################################################

preset_r.o: preset.h global_r.o $(MACROS) Makefile
	m4 -P $(MACROS) preset.h > $*.F90
	$(FC) $(FLAGS) $(CFLAGS) $(RFLAGS) -o preset_r.o -c $*.F90 $(FFTWcompile)
	@wc -l -L -w hdfint_r.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

preset_d.o: preset.h global_d.o $(MACROS) Makefile
	m4 -P $(MACROS) preset.h > $*.F90
	$(FC) $(FLAGS) $(CFLAGS) $(DFLAGS) -o preset_d.o -c $*.F90 $(FFTWcompile)
	@wc -l -L -w hdfint_d.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

###############################################################################################################################################################

%_r.o: %.f Makefile
	$(FC) $(FLAGS)           $(RFLAGS) -o $*_r.o -c $*.f
	@wc -l -L -w $*.f | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

%_d.o: %.f Makefile
	$(FC) $(FLAGS)           $(DFLAGS) -o $*_d.o -c $*.f
	@wc -l -L -w $*.f | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

###############################################################################################################################################################

$(ROBJS): %_r.o: %.F90 global_r.o $(MACROS) Makefile
	$(FC) $(FLAGS) $(CFLAGS) $(RFLAGS) -o $*_r.o -c $*.F90 $(FFTWcompile)
	@wc -l -L -w $*.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

$(DOBJS): %_d.o: %.F90 global_d.o $(MACROS) Makefile
	$(FC) $(FLAGS) $(CFLAGS) $(DFLAGS) -o $*_d.o -c $*.F90 $(FFTWcompile)
	@wc -l -L -w $*.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

$(F90FILES): %.F90: %.h 
	m4 -P $(MACROS) $*.h > $*.F90

###############################################################################################################################################################


###############################################################################################################################################################

xspech_r.o: xspech.h global_r.o $(addsuffix _r.o,$(files)) $(MACROS) Makefile
	@awk -v date='$(date)' -v pwd='$(PWD)' -v macros='$(MACROS)' -v f90='$(F90)' -v flags='$(FLAGS) $(CFLAGS) $(RFLAGS)' -v allfiles='$(ALLFILES)' \
	'BEGIN{nfiles=split(allfiles,files," ")} \
	{if($$2=="COMPILATION") {print "    write(ounit,*)\"      :  compiled  : date    = "date" ; \"" ; \
	                         print "    write(ounit,*)\"      :            : dir     = "pwd" ; \"" ; \
	                         print "    write(ounit,*)\"      :            : macros  = "macros" ; \"" ; \
	                         print "    write(ounit,*)\"      :            : f90     = "f90" ; \"" ; \
	                         print "    write(ounit,*)\"      :            : flags   = "flags" ; \"" }} \
	 {if($$2=="SUMTIME") {for (i=1;i<=nfiles;i++) print "   SUMTIME("files[i]")"}}\
	 {if($$2=="PRTTIME") {for (i=1;i<=nfiles;i++) print "   PRTTIME("files[i]")"}}\
	 {print}' xspech.h > mspech.h
	m4 -P $(MACROS) mspech.h > xspech.F90
	@rm -f mspech.h
	$(FC) $(FLAGS) $(CFLAGS) $(RFLAGS) -o xspech_r.o -c xspech.F90
	@wc -l -L -w xspech.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

xspech_d.o: xspech.h global_d.o $(addsuffix _d.o,$(files)) $(MACROS) Makefile
	@awk -v date='$(date)' -v pwd='$(PWD)' -v macros='$(MACROS)' -v f90='$(F90)' -v flags='$(FLAGS) $(CFLAGS) $(DFLAGS)' -v allfiles='$(ALLFILES)' \
	'BEGIN{nfiles=split(allfiles,files," ")} \
	{if($$2=="COMPILATION") {print "    write(ounit,*)\"      :  compiled  : date    = "date" ; \"" ; \
	                         print "    write(ounit,*)\"      :            : dir     = "pwd" ; \"" ; \
	                         print "    write(ounit,*)\"      :            : macros  = "macros" ; \"" ; \
	                         print "    write(ounit,*)\"      :            : f90     = "f90" ; \"" ; \
	                         print "    write(ounit,*)\"      :            : flags   = "flags" ; \"" }} \
	 {if($$2=="SUMTIME") {for (i=1;i<=nfiles;i++) print "   SUMTIME("files[i]")"}}\
	 {if($$2=="PRTTIME") {for (i=1;i<=nfiles;i++) print "   PRTTIME("files[i]")"}}\
	 {print}' xspech.h > mspech.h
	m4 -P $(MACROS) mspech.h > xspech.F90
	@rm -f mspech.h
	$(FC) $(FLAGS) $(CFLAGS) $(DFLAGS) -o xspech_d.o -c xspech.F90
	@wc -l -L -w xspech.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

###############################################################################################################################################################

clean:
	rm -f *.o ; rm -f *.mod ; rm -f *.F90 ; rm -f .*.h ; rm -f *.pdf ; rm -f *.dvi ; rm -f *.out ; rm -f *.bbl ; rm -f *.toc ; rm -f .*.date ; rm -rf ./docs/

###############################################################################################################################################################

./docs/%.pdf: %.h head.tex end.tex
	#emacs -r -fn 7x14 -g 160x80+280 $*.h
	mkdir -p ./docs/
	cd ./docs/ ; \
	pwd ; \
	@ls --full-time ../$*.h | cut -c 35-53 > .$*.date ; \
	awk -v file=$* -v date=.$*.date 'BEGIN{getline cdate < date ; FS="!latex" ; print "\\input{../head} \\code{"file"}"} \
	{if(NF>1) print $$2} \
	END{print "\\hrule \\vspace{1mm} \\footnotesize $*.h last modified on "cdate";" ; print "\\input{../end}"}' ../$*.h > $*.tex ; \
	latex $* ; latex $* ; latex $* ; dvips -P pdf -o $*.ps $*.dvi ; ps2pdf $*.ps ; \
	rm -f $*.tex $*.aux $*.blg $*.log $*.ps 

###############################################################################################################################################################

pdfs: $(addprefix ./docs/, $(addsuffix .pdf,$(HFILES))) head.html
	mkdir -p ./docs/
	cat head.html > ./docs/subroutines.html
	for file in $(HFILES) ; do grep "!title" $${file}.h | cut -c 7- | \
	                           awk -v file=$${file} -F!\
	                            '{print "<tr><td><a href="file".pdf\>"file"</a></td><td>"$$1"</td><td>"$$2"</td></tr>"}' \
	                            >> ./docs/subroutines.html ; \
	                          done
	echo "</table></body></html>" >> ./docs/subroutines.html
	@echo "Please view the pdfs in ./docs/ directory."

publish: $(addprefix ./docs/, $(addsuffix .pdf,$(HFILES))) ./docs/subroutines.html
# push local documentations online
	git stash
	git checkout gh-pages
	git pull origin gh-pages
	cp ./docs/* .
	git add *.pdf
	git add subroutines.html
	git commit -am "update documentations"
	git push origin gh-pages
	@echo "--------------------------------------------------------------------"
	@echo "Published the updated documentations."
	@echo "You are now in gh-pages branch."
	@echo "Please checkout back to your working branch by"
	@echo "$  git checkout <master>"
	@echo "If you have stashed local changes, you can recover them by "
	@echo "$ git stash pop "
	@echo "--------------------------------------------------------------------"

###############################################################################################################################################################

help:
	#
	# make 			: identical to make xspec ;
	# make xspec 		: expands macros (*.h --> *.F90) ; compiles xspec executable ;
	# make dspec 		: expands macros (*.h --> *.F90) ; compiles dspec executable ;
	# make clean 		: clean up compilation directory : rm -f *.o ; rm -f *.mod ; rm -f *.F90 ; rm -f *.pdf ; rm -f *.dvi
	# make pdfs		: create source documentation dvi, pdf files ; user should have directory $(HOME)/w3_html/Spec ;
	#
	# Compiler Control 	: CC=lff95; CC=intel_ipp; CC=gfortran_ipp
	# ---------------
	# macro expansion 	: m4 -P MACROS files.h > files.F90
	# ---------------
	# compilation		: FC FLAGS -o files.o -c files.F90 ; FC FLAGS -o xspec *files.o -LNAG -lnag -Lhdf5
	# -----------
	# defaults
	# --------
	# FC			= $(FC)
	# FLAGS			= $(FLAGS) $(CFLAGS) $(DFLAGS)
	# NAG			= $(NAG_ROOT)
	# MACROS		= $(MACROS)
	#

###############################################################################################################################################################
