#!/bin/sh

###############################################################################################################################################################

 afiles=manual rzaxis packxi volume coords
 bfiles=metrix ma00aa        matrix        mp00ac ma02aa packab tr00ab curent df00ab lforce
#cfiles=bc00aa fc02aa jk03aa pc00aa pc00ab
 cfiles=brcast dforce newton 
 dfiles=casing bnorml 
 efiles=jo00aa pp00aa pp00ab bfield stzxyz sc00aa
 ffiles=hesian ra00aa numrec
 sfiles=dcuhre minpack iqpack rksuite i1mach d1mach # below assumes the .f files are double precision; the CFLAGS = -r8 option is not required;

###############################################################################################################################################################

 SPECFILES=sphdf5 preset $(afiles) $(bfiles) $(cfiles) $(dfiles) $(efiles) $(ffiles)
 ALLFILES=global $(SPECFILES) $(sfiles) xspech
#F77FILES=$(sfiles:=.f)
 PREPROC=$(SPECFILES:=_m.F90) # preprocessed by m4
 RAWSOURCE=global $(SPECFILES) xspech # "raw" code, with macros not expanded yet

 ROBJS=$(SPECFILES:=_r.o)
 DOBJS=$(SPECFILES:=_d.o)
 
###############################################################################################################################################################
 
 MACROS=macros
 
 CC=intel
 # if want to use gfortran; make CC=gfortran xspec; otherwise using Intel
 FC=mpif90 # at PPPL, mpifort will cause parallel HDF5 hang
 
 # Intel Defaults
 # At PPPL, you can use the following commands
 # module use /p/focus/modules
 # module load spec
 CFLAGS=-r8
 RFLAGS=-mcmodel=large -O3 -m64 -unroll0 -fno-alias -ip -traceback
 DFLAGS=-O0 -g -traceback -check bounds -check format -check output_conversion -check pointers -check uninit -debug full -D DEBUG
 LIBS=-I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include  # MKL include
 LIBS+=-I$(HDF5_HOME)/include # HDF5 include
 LINKS=-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lpthread -lz -lm # HDF5 link
 LIBS+=-I$(FFTW_HOME)/include # FFTW include
 LINKS+=-L$(FFTW_HOME)/lib -lfftw3 # FFTW link
 LINKS+=${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a \
     -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a \
     ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl # MKL link


ifeq ($(CC),gfortran)
 # At PPPL, you can use the following commands
 # module use /p/focus/modules
 # module load spec/gcc
 RFLAGS=-O3 -w -ffree-line-length-none -fexternal-blas
 DFLAGS=-O0 -g -w -ffree-line-length-none -Wextra -Wtarget-lifetime -fbacktrace -fbounds-check -fexternal-blas \
     -fcheck=all -DDEBUG #-ffpe-trap=invalid,zero,overflow,underflow,inexact # for some reason this will cause crash
 CFLAGS=-fdefault-real-8
 LINKS=-L$(LAPACK_HOME) -llapack #-lgfortran
 LIBS=-I$(HDF5_HOME)/include
 LINKS+=-L$(HDF5_HOME)/lib -lhdf5hl_fortran -lhdf5 -lhdf5_fortran -lhdf5 -lpthread -lz -lm
 LIBS+=-I$(FFTW_HOME)/include
 LINKS+=-L$(FFTW_HOME)/lib -lfftw3 -lopenblas
endif

ifeq ($(CC),gfortran_ubuntu)
 # You should install the following packages
 # sudo apt install gfortran
 # sudo apt install libopenmpi-dev
 # sudo apt install liblapack-dev
 # sudo apt install m4
 # sudo apt install libfftw3-dev
 # sudo apt install libhdf5-openmpi-dev
 CFLAGS=-fdefault-real-8
 LINKS=-Wl,-rpath -Wl,/usr/lib/lapack -llapack -lblas
 LIBS=-I/usr/include/hdf5/openmpi
 LINKS+=-L/usr/lib/x86_64-linux-gnu/hdf5/openmpi -lhdf5_fortran -lhdf5 -lpthread -lz -lm
 LIBS+=-I/usr/include
 LINKS+=-lfftw3
 RFLAGS=-O2 -ffixed-line-length-none -ffree-line-length-none -fexternal-blas
 DFLAGS=-g -fbacktrace -fbounds-check -ffree-line-length-none -fexternal-blas -DDEBUG
endif 

ifeq ($(CC),gfortran_arch)
 # configuration for Arch Linux
 FC=mpif90
 CFLAGS=-fdefault-real-8
 LINKS=-llapack -lblas
 LIBS=
 LINKS+=-lhdf5_fortran -lhdf5 -lpthread -lz -lm
 LINKS+=-lfftw3
 RFLAGS=-O3 -ffixed-line-length-none -ffree-line-length-none -fexternal-blas
 DFLAGS=-g -fbacktrace -fbounds-check -ffree-line-length-none -fexternal-blas -DDEBUG
endif 

ifeq ($(CC),gfortran_mac)
 # works on Ksenia's laptop
 FC=mpif90
 CFLAGS=-fdefault-real-8
 LINKS=-L/usr/local/Cellar/lapack/3.8.0_1/lib -llapack -Wl,-rpath -Wl,/usr/local/Cellar/lapack/3.8.0_1/lib -lblas -L/usr/local/Cellar/openblas/0.3.7 -lgfortran -Wl,-rpath -Wl,/usr/local/Cellar/openblas/0.3.7/lib
 LIBS=-I/Applications/HDF_Group/HDF5/1.10.5/include/static
 LINKS+=-L/Applications/HDF_Group/HDF5/1.10.5/lib -lhdf5_f90cstub -lhdf5_fortran -lhdf5 -lpthread -lz -lm -Wl,-rpath -Wl,/Applications/HDF_Group/HDF5/1.10.5/lib
 LIBS+=-I/usr/local/Cellar/fftw/3.3.8/include
 LINKS+=-L/usr/local/Cellar/fftw/3.3.8/lib -lfftw3 -Wl,-rpath -Wl,/usr/local/Cellar/fftw/3.3.8/lib
 RFLAGS=-O2 -ffixed-line-length-none -ffree-line-length-none -fexternal-blas
 DFLAGS=-g -fbacktrace -fbounds-check -ffree-line-length-none -fexternal-blas -DDEBUG
endif

ifeq ($(CC),lff95)
 # LF95 SAL
 # Not checked
 CFLAGS=--dbl
 RFLAGS=--ap -O -I.
 DFLAGS=-Cpp -DDEBUG
 LINKS=-L$(LINKS_ROOT) -lnag -L$(LAPACKHOME) -llapack -L$(BLASHOME) -lblas
 LIBS=-I$(FFTWHOME)/include
 LINKS+=-L$(FFTWHOME)/lib -lfftw3
endif

ifeq ($(CC),intel_spc)
 CFLAGS=-r8
 RFLAGS=-O2 -ip -no-prec-div -xHost -fPIC
 DFLAGS=-traceback -D DEBUG -g
 LINKS=-L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
 LIBS=-I$(FFTW_DIR)/include
 LINKS+=-L$(FFTW_DIR)/lib -lfftw3
 LIBS+=-I$(HDF5_HOME)/include
 LINKS+=-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lpthread -lz -lm -Wl,-rpath -Wl,$(HDF5_HOME)/lib
endif

ifeq ($(CC),intel_ipp)
 # tested on draco with the following modules:
 # intel/18.0.3 impi/2018.3 mkl/2018.3 hdf5-serial/1.8.21 fftw-mpi/3.3.8
 # and on cobra with the following modules:
 # intel/19.0.4 impi/2019.4 mkl/2019.4 hdf5-serial/1.8.21 fftw-mpi/3.3.8
 FC=mpiifort
 CFLAGS=-r8
 RFLAGS=-O0 -ip -no-prec-div -xHost -fPIC
 DFLAGS=-traceback -D DEBUG
 LINKS=-L${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -lm -ldl -Wl,-rpath -Wl,${MKLROOT}/lib/intel64
 LIBS=-I$(HDF5_HOME)/include
 LINKS+=-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lpthread -lz -lm -Wl,-rpath -Wl,$(HDF5_HOME)/lib
 LIBS+=-I$(FFTW_HOME)/include
 LINKS+=-L$(FFTW_HOME)/lib -lfftw3 -Wl,-rpath -Wl,$(FFTW_HOME)/lib
endif

ifeq ($(CC),gfortran_ipp)
 # tested on draco with the following modules:
 # gcc/8 impi/2018.3 mkl/2018.3 hdf5-mpi/1.10.5 fftw-mpi/3.3.8
 FC=mpif90
 CFLAGS=-fdefault-real-8
 RFLAGS=-O2 -fPIC -ffree-line-length-none
 DFLAGS=-g -fbacktrace -fbounds-check -DDEBUG -ffree-line-length-none
 LINKS=-L${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -lm -ldl -Wl,-rpath -Wl,${MKLROOT}/lib/intel64
 LIBS=-I$(HDF5_HOME)/include
 LINKS+=-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lpthread -lz -lm -Wl,-rpath -Wl,$(HDF5_HOME)/lib
 LIBS+=-I$(FFTW_HOME)/include
 LINKS+=-L$(FFTW_HOME)/lib -lfftw3 -Wl,-rpath -Wl,$(FFTW_HOME)/lib
endif

ifeq ($(CC),intel_raijin)
 # One needs to load the following modules
 # module load intel-fc/2018.1.163
 # module load intel-cc/2018.1.163
 # module load intel-mkl/2018.1.163
 # module load openmpi
 # module load fftw3-mkl/2018.1.163
 # module load hdf5
 CFLAGS=-r8
 LINKS=-L${MKLROOT}/lib/intel64 -mkl 
 LIBS=-I$(HDF5_BASE)/include
 LINKS+=-L$(HDF5_BASE)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lpthread -lz -lm
 RFLAGS=-mcmodel=large -O3 -m64 -unroll0 -fno-alias -ip -traceback -fPIC
 DFLAGS=-check bounds -check format -check output_conversion -check pointers -check uninit -debug full -D DEBUG
endif

###############################################################################################################################################################

 WEBDIR=$(HOME)/w3_html

###############################################################################################################################################################

 date:=$(shell date)
 text:=$(shell date +%F)

###############################################################################################################################################################

xspec: $(addsuffix _r.o,$(ALLFILES)) $(MACROS) Makefile
	$(FC) $(FLAGS) $(CFLAGS) $(RFLAGS) -o xspec $(addsuffix _r.o,$(ALLFILES)) $(LINKS) 

dspec: $(addsuffix _d.o,$(ALLFILES)) $(MACROS) Makefile
	$(FC) $(FLAGS) $(CFLAGS) $(DFLAGS) -o dspec $(addsuffix _d.o,$(ALLFILES)) $(LINKS) 

###############################################################################################################################################################

global_r.o: %_r.o: global.f90 $(MACROS) 
	@awk -v allfiles='$(ALLFILES)' 'BEGIN{nfiles=split(allfiles,files," ")} \
	{if($$2=="CPUVARIABLE") {for (i=1;i<=nfiles;i++) print "  REAL    :: T"files[i]" = 0.0, "files[i]"T = 0.0"}}\
	{if($$2=="DSCREENLIST") {for (i=1;i<=nfiles;i++) print "  LOGICAL :: W"files[i]" = .false. "}}\
	{if($$2=="NSCREENLIST") {for (i=1;i<=nfiles;i++) print "  W"files[i]" , &"}}\
	{if($$2=="BSCREENLIST") {for (i=1;i<=nfiles;i++) print "  LlBCAST(W"files[i]",1,0)"}}\
	{if($$2=="WSCREENLIST") {s="'"'"'" ; d="'"\\\""'" ; for (i=1;i<=nfiles;i++) print "  if( W"files[i]" ) write(iunit,"s"("d" W"files[i]" = "d"L1)"s")W"files[i]}}\
	{print}' global.f90 > mlobal.f90
	m4 -P $(MACROS) mlobal.f90 > global_m.F90
	@rm -f mlobal.f90
	$(FC) $(FLAGS) $(CFLAGS) $(RFLAGS) -o global_r.o -c global_m.F90 $(LIBS)
	@wc -l -L -w global_m.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

global_d.o: %_d.o: global.f90 $(MACROS) 
	@awk -v allfiles='$(ALLFILES)' 'BEGIN{nfiles=split(allfiles,files," ")} \
	{if($$2=="CPUVARIABLE") {for (i=1;i<=nfiles;i++) print "  REAL    :: T"files[i]" = 0.0, "files[i]"T = 0.0"}}\
	{if($$2=="DSCREENLIST") {for (i=1;i<=nfiles;i++) print "  LOGICAL :: W"files[i]" = .false. "}}\
	{if($$2=="NSCREENLIST") {for (i=1;i<=nfiles;i++) print "  W"files[i]" , &"}}\
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
	$(FC) $(FLAGS)           $(RFLAGS) -o $*_r.o -c $*.f
	@wc -l -L -w $*.f | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

%_d.o: %.f 
	$(FC) $(FLAGS)           $(DFLAGS) -o $*_d.o -c $*.f
	@wc -l -L -w $*.f | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

###############################################################################################################################################################

$(ROBJS): %_r.o: %_m.F90 global_r.o $(MACROS) 
	$(FC) $(FLAGS) $(CFLAGS) $(RFLAGS) -o $*_r.o -c $*_m.F90 $(LIBS)
	@wc -l -L -w $*_m.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

$(DOBJS): %_d.o: %_m.F90 global_d.o $(MACROS) 
	$(FC) $(FLAGS) $(CFLAGS) $(DFLAGS) -o $*_d.o -c $*_m.F90 $(LIBS)
	@wc -l -L -w $*_m.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

$(PREPROC): %_m.F90: %.f90 $(MACROS)
	@awk -v file=$*.f90 '{ gsub("__LINE__", NR); gsub("__FILE__",file); print }' $*.f90 > $*_p.f90
	m4 -P $(MACROS) $*_p.f90 > $*_m.F90

###############################################################################################################################################################


###############################################################################################################################################################

xspech_r.o: xspech.f90 global_r.o $(addsuffix _r.o,$(files)) $(MACROS) 
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

xspech_d.o: xspech.f90 global_d.o $(addsuffix _d.o,$(files)) $(MACROS) 
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

./docs/%.pdf: %.f90 head.tex end.tex
	#emacs -r -fn 7x14 -g 160x80+280 $*.f90
	mkdir -p ./docs/
	cd ./docs/ ; \
	pwd ; \
	ls --full-time ../$*.f90 | awk 'IFS=" " {print $$6 " " substr($$7,0,8);}' > .$*.date ; \
	awk -v file=$* -v date=.$*.date 'BEGIN{getline cdate < date ; FS="!latex" ; print "\\input{../head} \\code{"file"}"} \
	{if(NF>1) print $$2} \
	END{print "\\hrule \\vspace{1mm} \\footnotesize $*.f90 last modified on "cdate";" ; print "\\input{../end}"}' ../$*_m.F90 > $*.tex ; \
	latex $* ; latex $* ; latex $* ; dvips -P pdf -o $*.ps $*.dvi ; ps2pdf $*.ps ; \
	rm -f $*.tex $*.aux $*.blg $*.log $*.ps 

###############################################################################################################################################################

pdfs: $(addprefix ./docs/, $(addsuffix .pdf,$(RAWSOURCE))) head.html
	mkdir -p ./docs/
	cat head.html > ./docs/subroutines.html
	for file in $(RAWSOURCE) ; do grep "!title" $${file}.f90 | cut -c 7- | \
	                           awk -v file=$${file} -F!\
	                            '{print "<tr><td><a href="file".pdf>"file"</a></td><td>"$$1"</td><td>"$$2"</td></tr>"}' \
	                            >> ./docs/subroutines.html ; \
	                          done
	echo "</table></body></html>" >> ./docs/subroutines.html
	@echo "Please view the pdfs in ./docs/ directory."

publish: $(addprefix ./docs/, $(addsuffix .pdf,$(RAWSOURCE))) ./docs/subroutines.html
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

show_makeflags:
	@echo " "
	@echo " "
	@echo "----------------------------------------- "
	@echo " "
	@echo " "
	@echo "Fortran preprocessing"
	@echo "=============="
	@echo "$(MACROS) $(ALLFILES)"
	@echo " "
	@echo "$(CC) compile"
	@echo "==============="
	@echo "$(FC) $(CFLAGS) $(RFLAGS)"
	@echo " "
	@echo "$(FC) $(CFLAGS) $(DFLAGS)"
	@echo " "
	@echo "Include libraries"
	@echo "==============="
	@echo "$(LIBS)"
	@echo " "
	@echo "LINKING FLAGS"
	@echo "==============="
	@echo "$(LINKS)"

###############################################################################################################################################################
