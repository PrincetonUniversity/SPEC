#!/bin/sh

###############################################################################################################################################################

 afiles=manual preset rzaxis packxi volume coords
 bfiles=metrix ma00aa        matrix        mp00ac ma02aa packab tr00ab curent df00ab lforce
#cfiles=bc00aa fc02aa jk03aa pc00aa pc00ab
 cfiles=brcast dforce newton 
 dfiles=casing bnorml 
 efiles=jo00aa pp00aa pp00ab bfield stzxyz sc00aa
 ffiles=hesian ra00aa numrec
 sfiles=dcuhre minpack

###############################################################################################################################################################

 SPECFILES=$(afiles) $(bfiles) $(cfiles) $(dfiles) $(efiles) $(ffiles)
 ALLFILES=global $(SPECFILES) $(sfiles) xspech hdfint
 PDFFILES=global $(SPECFILES) xspech hdfint
 F77FILES=$(sfiles:=.f)
 F90FILES= $(SPECFILES:=.F90)
 HFILES = $(SPECFILES:=.h)
 ROBJS=$(SPECFILES:=_r.o)
 DOBJS=$(SPECFILES:=_d.o)

###############################################################################################################################################################

 MACROS=macros
 CC=intel # if want to use gfortran; make CC=gfortran xfocus; otherwise using Intel
 FC=mpif90

 # Intel Defaults
 RFLAGS=-O3 -r8 -vec-report0 -fp-model strict -ip
 DFLAGS=-g  -r8 -traceback   -fp-model strict -check bounds -check format -check output_conversion -check pointers -check uninit -debug full -D DEBUG
 NAG=-L$(NAG_ROOT)/lib -lnag_nag 
 NETCDF=-L$(NETCDFHOME)/lib -lnetcdf
 HDF5compile=-I$(HDF5_HOME)/include
 HDF5link=-L$(HDF5_HOME)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lpthread -lz -lm

ifeq ($(CC),gfortran)
 # Not checked
 RFLAGS=-O2 -fdefault-real-8 -ffixed-line-length-none -ffree-line-length-none -fexternal-blas
 DFLAGS=-g3 -Wextra -Wtarget-lifetime -fbacktrace -fbounds-check -ffpe-trap=zero -fcheck=all -DDEBUG
endif

ifeq ($(CC),lff95)
 # LF95 SAL
 RFLAGS=-O2 --ap --dbl
 DFLAGS=-g  --ap --dbl -Cpp -DDEBUG
 NAG=-L$(NAG_ROOT) -lnag -L$(LAPACKHOME) -llapack -L$(BLASHOME) -lblas
endif

ifeq ($(CC),intel_prof)
 # LF95 SAL
 RFLAGS=-O2 -r8 -vec-report0 -fp-model strict -ip -p 
 DFLAGS=-g -r8 -traceback -fp-model strict -check bounds -check format -check output_conversion -check pointers -check uninit -debug full -D DEBUG
 NAG=-L$(NAG_ROOT)/lib -lnag_nag 
endif

ifeq ($(CC),intel_ipp)
 RFLAGS=-r8 -O2 -ip -no-prec-div -xHost -fPIC
 DFLAGS=-r8 -g -traceback -D DEBUG
 NAG=-L$(NAGFLIB_HOME)/lib -lnag_nag 
 NETCDF=-L$(NETCDF_HOME)/lib -lnetcdf
endif

ifeq ($(CC),intel_ipp_prof)
 RFLAGS=-r8 -O2 -ip -no-prec-div -xHost -fPIC -p
 DFLAGS=-r8 -g -traceback -D DEBUG -p
 NAG=-L$(NAGFLIB_HOME)/lib -lnag_nag 
 NETCDF=-L$(NETCDF_HOME)/lib -lnetcdf
endif

ifeq ($(CC),gfortran_ipp)
 RFLAGS=-fdefault-real-8 -O2 -fPIC -ffree-line-length-none
 DFLAGS=-fdefault-real-8 -g -fbacktrace -fbounds-check -DDEBUG -ffree-line-length-none
 NAG=-L$(NAGFLIB_HOME)/lib -lnag_nag 
endif

ifeq ($(CC),gfortran_ipp_prof)
 RFLAGS=-fdefault-real-8 -O2 -fPIC -ffree-line-length-none -p
 DFLAGS=-fdefault-real-8 -g -fbacktrace -fbounds-check -DDEBUG -ffree-line-length-none -p
 NAG=-L$(NAGFLIB_HOME)/lib -lnag_nag 
endif

###############################################################################################################################################################

 WEBDIR=$(HOME)/w3_html

###############################################################################################################################################################

 date:=$(shell date)
 text:=$(shell date +%F)

###############################################################################################################################################################

xspec: $(addsuffix _r.o,$(ALLFILES)) $(MACROS) Makefile
	$(FC) $(RFLAGS)           -o xspec $(addsuffix _r.o,$(ALLFILES)) $(NAG) $(HDF5compile) $(HDF5link) $(NETCDF)
	date
	/bin/echo -e "\a"

dspec: $(addsuffix _d.o,$(ALLFILES)) $(MACROS) Makefile
	$(FC) $(RFLAGS) $(DFLAGS) -o dspec $(addsuffix _d.o,$(ALLFILES)) $(NAG) $(HDF5compile) $(HDF5link) $(NETCDF)
	date
	/bin/echo -e "\a"

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
	$(FC) $(RFLAGS)           -o global_r.o -c global.F90
	@wc -l -L -w global.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

global_d.o: global.h $(MACROS) Makefile
	@awk -v allfiles='$(ALLFILES)' 'BEGIN{nfiles=split(allfiles,files," ")} \
	{if($$2=="CPUVARIABLE") {for (i=1;i<=nfiles;i++) print "  REAL    :: T"files[i]" = 0.0, "files[i]"T = 0.0"}}\
	{if($$2=="DSCREENLIST") {for (i=1;i<=nfiles;i++) print "  LOGICAL :: W"files[i]" = .false. "}}\
	{if($$2=="NSCREENLIST") {for (i=1;i<=nfiles;i++) print "  W"files[i]" , &"}}\
	{if($$2=="BSCREENLIST") {for (i=1;i<=nfiles;i++) print "  LlBCAST(W"files[i]",1,0)"}}\
	{if($$2=="WSCREENLIST") {s="'"'"'" ; d="'"\\\""'" ; for (i=1;i<=nfiles;i++) print "  if( W"files[i]" ) write(iunit,"s"("d" W"files[i]" = "d"L1)"s")W"files[i]}}\
	{print}' global.h > mlobal.h
	m4 -P $(MACROS) mlobal.h > global.F90
	@rm -f mlobal.h
	$(FC) $(RFLAGS) $(DFLAGS) -o global_d.o -c global.F90
	@wc -l -L -w global.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

###############################################################################################################################################################

hdfint_r.o: hdfint.h global_r.o $(MACROS) Makefile
	m4 -P $(MACROS) hdfint.h > $*.F90
	$(FC) $(RFLAGS)           -o hdfint_r.o -c $*.F90 $(HDF5compile)
	@wc -l -L -w hdfint_r.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

hdfint_d.o: hdfint.h global_d.o $(MACROS) Makefile
	m4 -P $(MACROS) hdfint.h > $*.F90
	$(FC) $(RFLAGS) $(DFLAGS) -o hdfint_d.o -c $*.F90 $(HDF5compile)
	@wc -l -L -w hdfint_d.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

###############################################################################################################################################################

%_r.o: %.f Makefile
	$(FC) $(RFLAGS)           -o $*_r.o -c $*.f
	@wc -l -L -w $*.f | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

%_d.o: %.f Makefile
	$(FC) $(RFLAGS) $(DFLAGS) -o $*_d.o -c $*.f
	@wc -l -L -w $*.f | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

###############################################################################################################################################################

$(ROBJS): %_r.o: %.F90 global_r.o $(MACROS) Makefile
	$(FC) $(RFLAGS)           -o $*_r.o -c $*.F90
	@wc -l -L -w $*.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

$(DOBJS): %_d.o: %.F90 $(MACROS) Makefile
	$(FC) $(RFLAGS) $(DFLAGS) -o $*_d.o -c $*.F90
	@wc -l -L -w $*.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

$(F90FILES): %.F90: %.h 
	m4 -P $(MACROS) $*.h > $*.F90

###############################################################################################################################################################


###############################################################################################################################################################

xspech_r.o: xspech.h global_r.o $(addsuffix .o,$(files)) $(MACROS) Makefile
	@awk -v date='$(date)' -v pwd='$(PWD)' -v macros='$(MACROS)' -v f90='$(F90)' -v flags='$(RFLAGS)'           -v allfiles='$(ALLFILES)' \
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
	$(FC) $(RFLAGS)           -o xspech_r.o -c xspech.F90
	@wc -l -L -w xspech.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

xspech_d.o: xspech.h global_d.o $(addsuffix .o,$(files)) $(MACROS) Makefile
	@awk -v date='$(date)' -v pwd='$(PWD)' -v macros='$(MACROS)' -v f90='$(F90)' -v flags='$(RFLAGS) $(DFLAGS)' -v allfiles='$(ALLFILES)' \
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
	$(FC) $(RFLAGS) $(DFLAGS) -o xspech_d.o -c xspech.F90
	@wc -l -L -w xspech.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

###############################################################################################################################################################

clean:
	rm -f *.o ; rm -f *.mod ; rm -f *.F90 ; rm -f .*.h ; rm -f *.pdf ; rm -f *.dvi ; rm -f *.out ; rm -f *.bbl ; rm -f *.toc ; rm -f .*.date

###############################################################################################################################################################

tar:
	tar -cf spec.$(text).tar *.h macros Makefile
	gzip --best spec.$(text).tar ; mv spec.$(text).tar.gz Tarfiles/.

###############################################################################################################################################################

%.pdf: %.h head.tex end.tex Makefile
	#emacs -r -fn 7x14 -g 160x80+280 $*.h
	@ls --full-time $*.h | cut -c 35-53 > .$*.date
	awk -v file=$* -v date=.$*.date 'BEGIN{getline cdate < date ; FS="!latex" ; print "\\input{head} \\code{"file"}"} \
	{if(NF>1) print $$2} \
	END{print "\\hrule \\vspace{1mm} \\footnotesize $*.h last modified on "cdate";" ; print "\\input{end}"}' $*.h > $*.tex
	@echo "-------------------------------------------------------------------------------------------------------------------------------"
	@echo $*
	@echo "-------------------------------------------------------------------------------------------------------------------------------"
	latex $* ; latex $* ; latex $*
	dvips -P pdf -o $*.ps $*.dvi ; ps2pdf $*.ps
	rm -f $*.tex $*.aux $*.blg $*.log $*.ps

###############################################################################################################################################################

pdfs: $(addsuffix .pdf,$(PDFFILES)) head.html
ifeq ($(USER),shudson)

	cat head.html > $(WEBDIR)/Spec/subroutines.html

	for file in $(PDFFILES) ; do cp $${file}.pdf $(WEBDIR)/Spec/. ; grep "!title" $${file}.h | cut -c 7- | \
	                           awk -v file=$${file} -F!\
	                            '{print "<tr><td><a href="file".pdf\">"file"</a></td><td>"$$1"</td><td>"$$2"</td></tr>"}' \
	                            >> $(WEBDIR)/Spec/subroutines.html ; \
	                          done

	echo "</table></body></html>" >> $(WEBDIR)/Spec/subroutines.html
else
	@echo "-------------------------------------------------------------------------------------------------------------------------------"
	@echo "Please read pdfs at www.ppl.gov/~shudson/Spec/subroutines.html or on GitHub pages."
	@echo "-------------------------------------------------------------------------------------------------------------------------------"

endif

###############################################################################################################################################################

help:
	#
	# make 			: identical to make xspec ;
	# make xspec 		: expands macros (*.h --> *.F90) ; compiles xspec executable ;
	# make dspec 		: expands macros (*.h --> *.F90) ; compiles dspec executable ;
	# make clean 		: clean up compilation directory : rm -f *.o ; rm -f *.mod ; rm -f *.F90 ; rm -f *.pdf ; rm -f *.dvi
	# make tar		: create archive (tar) of essential files ; save in Tarfiles/ ;
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
	# FLAGS			= $(RFLAGS) $(DFLAGS)
	# NAG			= $(NAG_ROOT)
	# MACROS		= $(MACROS)
	#

###############################################################################################################################################################
