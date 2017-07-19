#!/bin/sh

###############################################################################################################################################################

 afiles=manual global xspech        preset rzaxis packxi volume coords
 bfiles=metrix ma00aa        matrix        mp00ac ma02aa packab tr00ab curent df00ab lforce
#cfiles=bc00aa fc02aa jk03aa pc00aa pc00ab
 cfiles=brcast dforce newton 
 dfiles=casing bnorml 
 efiles=jo00aa pp00aa pp00ab bfield stzxyz sc00aa
 ffiles=hesian hdfint ra00aa numrec

###############################################################################################################################################################

 allfiles=$(afiles) $(bfiles) $(cfiles) $(dfiles) $(efiles) $(ffiles)

###############################################################################################################################################################

 MACROS=macros
 CC=intel # if want to use gfortran; make CC=gfortran xfocus; otherwise using Intel
 FC=mpif90

 # Intel Defaults
 RFLAGS=-r8 -mcmodel=large -O3 -m64 -unroll0 -fno-alias -ip -traceback 
 DFLAGS=-check bounds -check format -check output_conversion -check pointers -check uninit -debug full -D DEBUG
 NAG=-L$(NAG_ROOT)/lib -lnag_nag 

ifeq ($(CC),gfortran)
 # Not checked
 RFLAGS=-O2 -fdefault-real-8 -ffixed-line-length-none -ffree-line-length-none -fexternal-blas
 DFLAGS=-g3 -Wextra -Wtarget-lifetime -fbacktrace -fbounds-check -ffpe-trap=zero -fcheck=all -DDEBUG
endif

ifeq ($(CC),lff95)
 # LF95 SAL
 FLAGS=--ap --dbl -O
 DFLAGS=
 NAG=-L$(NAG_ROOT) -lnag
endif


 NETCDF=-L$(NETCDFHOME)/lib -lnetcdf

 HDF5compile=-I$(HDF5_HOME)/include
 HDF5link=-L$(HDF5_HOME)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lpthread -lz -lm

###############################################################################################################################################################

 WEBDIR=$(HOME)/w3_html

###############################################################################################################################################################

 date:=$(shell date)
 text:=$(shell date +%F)

###############################################################################################################################################################

xspec: $(addsuffix .o,$(allfiles)) $(MACROS) Makefile
	$(FC) $(FLAGS) $(DFLAGS) -o xspec $(addsuffix .o,$(allfiles)) $(NAG) $(HDF5compile) $(HDF5link) $(NETCDF)
	date
	/bin/echo -e "\a"

dspec: $(addsuffix .o,$(allfiles)) $(MACROS) Makefile
	$(FC) $(FLAGS) $(DFLAGS) -o dspec $(addsuffix .o,$(allfiles)) $(NAG) $(HDF5compile) $(HDF5link) $(NETCDF)
	date
	/bin/echo -e "\a"

###############################################################################################################################################################

global.o: global.h $(MACROS) Makefile
	@awk -v allfiles='$(allfiles)' 'BEGIN{nfiles=split(allfiles,files," ")} \
	{if($$2=="CPUVARIABLE") {for (i=1;i<=nfiles;i++) print "  REAL    :: T"files[i]" = 0.0, "files[i]"T = 0.0"}}\
	{if($$2=="DSCREENLIST") {for (i=1;i<=nfiles;i++) print "  LOGICAL :: W"files[i]" = .false. "}}\
	{if($$2=="NSCREENLIST") {for (i=1;i<=nfiles;i++) print "  W"files[i]" , &"}}\
	{if($$2=="BSCREENLIST") {for (i=1;i<=nfiles;i++) print "  LlBCAST(W"files[i]",1,0)"}}\
	{if($$2=="WSCREENLIST") {s="'"'"'" ; d="'"\\\""'" ; for (i=1;i<=nfiles;i++) print "  if( W"files[i]" ) write(iunit,"s"("d" W"files[i]" = "d"L1)"s")W"files[i]}}\
	{print}' global.h > mlobal.h
	m4 -P $(MACROS) mlobal.h > global.F90
	@rm -f mlobal.h
	$(FC) $(FLAGS) $(DFLAGS) -o global.o -c global.F90
	@wc -l -L -w global.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

###############################################################################################################################################################

hdfint.o: hdfint.h global.o $(MACROS) Makefile
	m4 -P $(MACROS) hdfint.h > $*.F90
	$(FC) $(FLAGS) $(DFLAGS) -o hdfint.o -c $*.F90 $(HDF5compile)
	@wc -l -L -w hdfint.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

###############################################################################################################################################################

%.o: %.h global.o $(MACROS) Makefile
	m4 -P $(MACROS) $*.h > $*.F90
	$(FC) $(FLAGS) $(DFLAGS) -o $*.o -c $*.F90
	@wc -l -L -w $*.F90 | awk '{print $$4" has "$$1" lines, "$$2" words, and the longest line is "$$3" characters ;"}'
	@echo ''

###############################################################################################################################################################

xspech.o: xspech.h global.o $(addsuffix .o,$(files)) $(MACROS) Makefile
	@awk -v date='$(date)' -v pwd='$(PWD)' -v macros='$(MACROS)' -v f90='$(F90)' -v flags='$(FLAGS) $(DFLAGS)' -v allfiles='$(allfiles)' \
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
	$(FC) $(FLAGS) $(DFLAGS) -o xspech.o -c xspech.F90
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

pdfs: $(addsuffix .pdf,$(allfiles)) head.html
ifeq ($(USER),shudson)

	cat head.html > $(WEBDIR)/Spec/subroutines.html

	for file in $(allfiles) ; do cp $${file}.pdf $(WEBDIR)/Spec/. ; grep "!title" $${file}.h | cut -c 7- | \
	                           awk -v file=$${file} -F!\
	                            '{print "<tr><td><a href="file".pdf\">"file"</a></td><td>"$$1"</td><td>"$$2"</td></tr>"}' \
	                            >> $(WEBDIR)/Spec/subroutines.html ; \
	                          done

	echo "</table></body></html>" >> $(WEBDIR)/Spec/subroutines.html
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
	# macro expansion 	: m4 -P MACROS files.h > files.F90
	# ---------------
	# compilation		: FC FLAGS -o files.o -c files.F90 ; FC FLAGS -o xspec *files.o -LNAG -lnag -Lhdf5
	# -----------
	# defaults
	# --------
	# FC			= $(FC)
	# FLAGS			= $(FLAGS) $(DFLAGS)
	# NAG			= $(NAG_ROOT)
	# MACROS		= $(MACROS)
	#

###############################################################################################################################################################
