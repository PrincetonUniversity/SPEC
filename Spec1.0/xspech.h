!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Main program.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! uppercase generally indicates macros;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

program xspech

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  use constants, only : zero, one, two, half, ten, pi, pi2

  use numerical, only : vsmall, sqrtmachprec, &
                        logtolerance

  use fileunits, only : ounit, zunit, wunit

  use inputlist, only : Wmacros, Wxspech, Wsc00aa, ext, &
                        Nfp, Mpol, Ntor, Igeometry, Nvol, &
                        pressure, pscale, Ladiabatic, adiabatic, gamma, phiedge, &
                        Rwc, Zws, Rws, Zwc, helicity, pflux, &
                        LBeltrami, Lconstraint, &
                        Lrad, &
                        Lfreebound, curtor, curpol, norblend, maxfbits, Lwall, phiwall, &
                        Linitialize, Iswmin, &
                        Lminimize, Lfindzero, &
                        ForceErr, &
                        odetol, nPpts, nPtrj, Mpqits, npq, &
!                       Lwrpj, &
                        LHevalues, LHevectors, Lperturbed, Lcheck, &
                        Rbc, Zbs, Rbs, Zbc, &
                        epsr, &
                        pqs, pqt

  use cputiming, only : Txspech

  use allglobal, only : readin, writin, ncpu, myid, cpus, pi2nfp, &
                        Mvol, &
                        YESstellsym, NOTstellsym, &
                        Iquad, &
                        mn, im, in, Nt, Nz, Ntz, &
                        Ltangent, ivol, &
                        lgeometricaldof, Ngeometricaldof, &
                        DifferentiateGeometry, &
                        iRbc, iZbs, iRbs, iZbc, &
                        iBns, iBnc, &
                        BBe, IIo, BBo, IIe, &
                        MBpsi, MEpsi, &
                        Bsubtemn, Bsubzemn, Bsubtomn, Bsubzomn, &
!                       Iwrpj, &
                        vvolume, &
                        Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, &
                        Lcontinueiterations, nFreeiterations, &
                        dpflux, &
                        Lmgridexist, &
                        ImagneticOK, &
                        pqorbit, &
                        Llatex

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS

  LOGICAL              :: LComputeDerivatives

  INTEGER              :: vvol, llmodnp, ifail, wflag, iflag, mm, nn, imn, iwa00aa, Lcurvature, lquad, mi, ni, innout, ii, jj, kk, jk, vflag
  REAL                 :: rflag, lcpu, lss, lastcpu, teta, zeta, st(1:6), Bst(1:6), normalerror(0:1)
  REAL,    allocatable :: position(:), gradient(:)
  REAL,    allocatable :: original(:)
  REAL,    allocatable :: oRbc(:,:), oZbs(:,:), oRbs(:,:), oZbc(:,:)
  CHARACTER            :: packorunpack

  LOGICAL              :: prime

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  MPISTART

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( myid.eq.0 ) then
   write(ounit,'("xspech : ", 10x ," : ")')
! COMPILATION ! this line is replaced during compilation with a write statement identifying date, time, compilation flags, etc.; see Makefile;
  endif
  
#ifdef CHECKNAG
  if( myid.eq.0 ) call A00AAF() ! check NAG version; ! 23 Oct 12;
#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  cpus = GETTIME ! set initial time; 04 Dec 14;

  cpuo = cpus
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( myid.eq.0 ) then
   write(ounit,'("xspech : ", 10x ," : ")')
   write(ounit,'("xspech : ",f10.2," : begin execution ; ncpu=",i3," ; calling readin ;")') cpus-cpus, ncpu
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \subsubsection{workflow}

!latex \item The input namelists and geometry are read in via a call to \verb+readin+. A full description of the required input is given in \verb+globals+.
  
  WCALL(xspech, readin )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item The global memory is allocated in \verb+al00aa+.

  WCALL(xspech, al00aa )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! looperr = one ! perhaps it is better to initialize this in al00aa; 23 Oct 12;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( myid.eq.0 ) then ! prepare ``convergence evolution'' output; save restart file;

   open( zunit, file = "."//trim(ext)//".iterations", status = "unknown", form = "unformatted" ) ! this file is written to in globals/writin; 11 Aug 14;
   write(zunit) mn, Mvol, Nfp
   write(zunit) im(1:mn)
   write(zunit) in(1:mn)

   wflag = 0 ; iflag = 0 ; rflag = zero
   WCALL(xspech, writin,( wflag, iflag, rflag )) ! write restart file etc. ! 14 Jan 13;

  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \item The boundary is adapted to the wall, if \verb+Lwall.eq.1+.
  
!  if( Lfreebound.gt.0 ) then ! prepare vacuum field calculation; 22 Apr 13;
!   
!   
!   if( Lwall.ne.0 ) then ! shall construct a suitable computational boundary; 01 May 13;
!    
!    FATALMESS(xspech, phiwall.le.zero .or. phiwall.ge.one, please choose phiwall between zero and one )
!
!    if( myid.eq.0 ) then
!     cput = GETTIME
!     write(ounit,'("xspech : ", 10x ," : ")') 
!     write(ounit,'("xspech : ",f10.2," : constructing computational boundary ;")') cput-cpus
!     WCALL(xspech,wa00aa,(iwa00aa)) ! construct iso-surface of Laplace's equation between plasma boundary and wall; wall.dat needs to be provided by user; 24 Oct 12;
!    endif ! end of if( myid.eq.0 );
!    
!    IlBCAST(iwa00aa,1,0) ! this is an error flag that indicates successful construction of isosurface = computational boundary; 15 May 13;
!    
!    FATALMESS(xspech, iwa00aa.ne.0, failed to construct computational boundary )
!    
!    if( iwa00aa.eq.0 ) then ! computational boundary successfully constructed; 15 May 13;
!     RlBCAST(iRbc(1:mn,Mvol), mn, 0 ) ! this was updated in wa00aa/wa01aa/. . .
!     RlBCAST(iZbs(1:mn,Mvol), mn, 0 )
!     if( NOTstellsym ) then
!     RlBCAST(iRbs(1:mn,Mvol), mn, 0 ) ! this was updated in wa00aa/wa01aa/. . .
!     RlBCAST(iZbc(1:mn,Mvol), mn, 0 )
!     endif
!     do imn = 1, mn ; mm = im(imn) ; nn = in(imn) / Nfp ; Rwc(nn,mm) = iRbc(imn,Mvol) ; Zws(nn,mm) = iZbs(imn,Mvol) ! updated boundary written to restart file;
!      if( NOTstellsym ) then ;                          ; Rws(nn,mm) = iRbs(imn,Mvol) ; Zwc(nn,mm) = iZbc(imn,Mvol)
!      endif
!     enddo
!    !iBns(1:mn) = zero ! if computational boundary is changed, the given normal field is presumably incorrect; perhaps safest to start from B.n=zero; 10 Apr 13;
!    !if( NOTstellsym ) then 
!    !iBnc(1:mn) = zero
!    !endif
!    endif ! end of if( iwa00aa.eq.0 ); computational boundary has been constructed between plasma boundary and wall; 15 May 13;
!    
!   endif ! end of if( Lwall.ne.0 );
!
!   wflag = 0 ; iflag = 0 ; rflag = zero
!   WCALL(xspech, writin,( wflag, iflag, rflag )) ! write restart file etc. ! 14 Jan 13;   
!   
!
!   cput = GETTIME
!   if( myid.eq.0 ) write(ounit,'("xspech : ", 10x ," : ")') 
!   if( maxfbits.gt.0 ) then
!    if( myid.eq.0 ) write(ounit,'("xspech : ",f10.2," : Lfreebound.gt.0 and maxfbits.gt.0 => calling mg00aa to construct EZspline interpolation of mgrid;")') cput-cpus
!    WCALL(xspech, mg00aa ) ! read in mgrid data and interpolate using ezspline; this will set Lmgridexist = .true. if successful;
!   else
!    if( myid.eq.0 ) write(ounit,'("xspech : ",f10.2," : Lfreebound.gt.0 but maxfbits.le.0 => no need to construct EZspline interpolation of mgrid;")') cput-cpus
!   endif
!   if( myid.eq.0 ) write(ounit,'("xspech : ", 10x ," : ")') 
!   
!   
!  endif ! end of if( Lfreebound ) ; 15 May 13;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!  if( Igeometry.eq.3 .and. Iswmin.gt.0 ) then 
!   
!   if( myid.eq.0 ) then ! screen output; 17 Apr 13;
!    cput = GETTIME
!    write(ounit,'("xspech : ", 10x ," : ")') 
!    write(ounit,'("xspech : ",f10.2," : constructing spectrally reduced angle ; Iswmin="i2" ;")') cput-cpus, Iswmin
!   endif
!
!   do vvol = 1, Mvol ! loop over ALL interfaces (including computational boundary if relevant); coordinate axis is excluded; 15 May 13;
!    if( myid.ne.modulo(vvol-1,ncpu) ) cycle
!    WCALL(xspech,sw00ac,( vvol, mn, Ntz, iRbc(1:mn,vvol), iZbs(1:mn,vvol), iRbs(1:mn,vvol), iZbc(1:mn,vvol) )) ! construct spectrally reduced Fourier representation; 15 May 13;
!   enddo
!   
!   do vvol = 1, Mvol ! broadcast spectrally condensed representation of interfaces; 17 Apr 13; ! coordinate axis is excluded; 13 Sep 13;
!    llmodnp = modulo(vvol-1,ncpu)
!    RlBCAST(iRbc(1:mn,vvol),mn,llmodnp)
!    RlBCAST(iZbs(1:mn,vvol),mn,llmodnp)
!    if( NOTstellsym ) then
!    RlBCAST(iRbs(1:mn,vvol),mn,llmodnp)
!    RlBCAST(iZbc(1:mn,vvol),mn,llmodnp)
!    endif
!   enddo
!   
!   if( Iswmin.eq.2 ) then ! write spectrally condensed plasma/computational boundary to restart file; 15 May 13;
!    do ii = 1, mn ; mm = im(ii) ; nn = in(ii) / Nfp ; Rbc(nn,mm) = iRbc(ii,Nvol) ; Zbs(nn,mm) = iZbs(ii,Nvol) ! updated boundary written to restart file;
!     if( NOTstellsym ) then ;                       ; Rbs(nn,mm) = iRbs(ii,Nvol) ; Zbc(nn,mm) = iZbc(ii,Nvol)
!     endif
!    enddo
!    if( Lfreebound.gt.0 ) then
!    do ii = 1, mn ; mm = im(ii) ; nn = in(ii) / Nfp ; Rwc(nn,mm) = iRbc(ii,Mvol) ; Zws(nn,mm) = iZbs(ii,Mvol) ! updated boundary written to restart file;
!     if( NOTstellsym ) then ;                       ; Rws(nn,mm) = iRbs(ii,Mvol) ; Zwc(nn,mm) = iZbc(ii,Mvol)
!     endif
!    enddo
!    endif
!   endif ! end of if( Iswmin.eq.2) ; 15 May 13;
!
!   wflag = 0 ; iflag = 0 ; rflag = zero
!   WCALL(xspech, writin,( wflag, iflag, rflag )) ! write restart file etc. ! 14 Jan 13;
!   
!   if( Iswmin.ge.3 ) goto 9999
!
!  endif ! end of if( Igeometry.eq.3 .and. Iswmin.gt.0 ); 15 May 13;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  FATALMESS(xspech, Ngeometricaldof.lt.0, counting error)

  RALLOCATE(position,(0:Ngeometricaldof)) ! position ; Ngeometricaldof was computed in al00aa;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! FREE BOUNDARY ITERATION LOOP; ! this will eventually be redundant; 03 Apr 13;
  
  nFreeIterations = -1
  
9000 nFreeIterations = nFreeIterations + 1
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Ngeometricaldof.gt.0 ) then ! pack geometry into vector; 14 Jan 13;

   packorunpack = 'P'
   WCALL(xspech,gf00aa,( Ngeometricaldof, position(0:Ngeometricaldof), Mvol, mn, &
iRbc(1:mn,0:Mvol), iZbs(1:mn,0:Mvol), iRbs(1:mn,0:Mvol), iZbc(1:mn,0:Mvol), packorunpack ))

  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \item Initialize adiabatic constants.
  
  do vvol = 1, Mvol

   vflag = 0 ! this flag instructs vo00aa to terminate if the volume is invalid; 04 Dec 14;
   WCALL(xspech,vo00aa,( vvol, vflag )) ! compute volume;

   if( Ladiabatic.eq.0 ) adiabatic(vvol) = pressure(vvol) * vvolume(vvol)**gamma ! initialize adiabatic constants using supplied pressure profile;

  enddo ! end of do vvol = 1, Nvol;
  
  if( Mvol.gt.Nvol ) then ; adiabatic(Mvol) = zero ; pressure(Mvol) = zero ! these are never used; 15 May 13;
  endif
  
  if( Wxspech .and. myid.eq.0 ) then
   cput = GETTIME
   write(ounit,'("xspech : ",f10.2," : myid=",i3," ; adiabatic constants = "999es13.5)')cput-cpus, myid, adiabatic(1:Mvol)
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( Ngeometricaldof.gt.0 ) then !latex \item The following routines are called only if geometrical degrees of freedom exist; \begin{itemize}
   
!  Iwrpj = .false. ! controls construction of pressure-jump Hamiltonian; 20 Jun 14;
   
!latex \item if \verb+Lminimize.eq.1+, call \verb+pc00aa+ to find minimum of energy functional
!latex using quasi-Newton, preconditioned conjugate gradient method, \verb+E04DGF+;

    if( Lminimize.eq.1 ) then

     ifail = 1 ! this is probably not required; 26 Feb 13;

     WCALL(xspech,pc00aa,( Ngeometricaldof, position(1:Ngeometricaldof), Mvol, mn, ifail ))

    endif
   
!latex \item if \verb+Lminimize.eq.2+, call \verb+pc01aa+ to find minimum of energy functional using ``home-spun'' non-linear conjugate gradient method;
!   if( Lminimize.eq.2 ) then
!    WCALL(xspech,pc01aa,( Nvol, mn, Ngeometricaldof, position(1:Ngeometricaldof) )) ! SHOULD INCLUDE IFAIL ARGUMENT FOR CONSISTENCY;
!   endif
   
!latex \item if \verb+Lminimize.eq.4+, call \verb+pc02aa+ to find minimum of energy functional using using modified Newton method, \verb+E04LYF+;
!latex under construction;
!   if( Lminimize.eq.4 ) then
!    WCALL(xspech,pc02aa,( Nvol, mn, Ngeometricaldof, position(0:Ngeometricaldof), ifail ))
!   endif
   
!latex \item If \verb+Lfindzero.gt.0+, call \verb+jk03aa+ to find extremum of energy function using a Newton method, \verb+C05PDF+;
   
   if( Lfindzero.gt.0 ) then

    ifail = 1 ! this is probably not required; 26 Feb 13;
    WCALL(xspech,jk03aa,( Ngeometricaldof, position(0:Ngeometricaldof), ifail ))

   endif
   
   packorunpack = 'U' ! unpack geometrical degrees of freedom; 13 Sep 13;
   WCALL(xspech,gf00aa,( Ngeometricaldof, position(0:Ngeometricaldof), Mvol, mn, &
iRbc(1:mn,0:Mvol), iZbs(1:mn,0:Mvol), iRbs(1:mn,0:Mvol), iZbc(1:mn,0:Mvol), packorunpack ))

  endif !latex \end{itemize}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! set/reset input variables;
  
! if( Lconstraint.lt.2 ) helicity(1:Nvol) = lABintegral(1:Nvol) ! updated ``input'' quantity;

#ifdef DEBUG
  do vvol = 1, Mvol
   FATALMESS(xspech, vvolume(vvol).lt.vsmall, error dividing adiabatic by volume)
  enddo
#endif
  
  pressure(1:Mvol) = adiabatic(1:Mvol) / vvolume(1:Mvol)**gamma ! this matches construction of adiabatic above;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! the following are diagnostics; set appropriate flags;
  
! Iwrpj = Lwrpj ! cb00aa (which is called by fc02aa, which is called by jk03aa) will write the pressure-jump files, depending on input flag Lwrpj;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \item The Beltrami fields in each annulus are calculated by a call to \verb+fc02aa+.
  
  RALLOCATE(gradient,(0:Ngeometricaldof))
  
  lastcpu = GETTIME
  
  LComputeDerivatives = .false.

  WCALL( xspech, fc02aa, ( Ngeometricaldof, position(0:Ngeometricaldof), gradient(0:Ngeometricaldof), LComputeDerivatives ) ) ! (re-)calculate Beltrami fields;
  
  DEALLOCATE(gradient)
  
#ifdef DEBUG
  do vvol = 1, Mvol-1
   ; FATALMESS(xspech, BBe(vvol).lt.logtolerance, underflow)
   if( Igeometry.eq.3 .or. Igeometry.eq.4 ) then ! include spectral constraints; 04 Dec 14;
    ;FATALMESS(xspech, IIo(vvol).lt.logtolerance, underflow)
   endif
   if( NOTstellsym ) then
    ;FATALMESS(xspech, BBo(vvol).lt.logtolerance, underflow)
    if( Igeometry.eq.3 .or. Igeometry.eq.4 ) then ! include spectral constraints; 04 Dec 14;
     FATALMESS(xspech, IIe(vvol).lt.logtolerance, underflow)
    endif
   endif
  enddo
#endif

  if( myid.eq.0 ) then
   cput = GETTIME
   write(ounit,1000) cput-cpus, myid,                             ForceErr,  cput-lastcpu, "|BB|e", alog10(BBe(1:min(Mvol-1,28)))
   if( Igeometry.ge.3 .or. Igeometry.eq.4 ) then ! include spectral constraints; 04 Dec 14;
   write(ounit,1001)            myid,                                                      "|II|o", alog10(IIo(1:min(Mvol-1,28)))
   endif
   if( NOTstellsym ) then
   write(ounit,1001)            myid,                                                      "|BB|o", alog10(BBo(1:min(Mvol-1,28)))
   if( Igeometry.ge.3 .or. Igeometry.eq.4 ) then ! include spectral constraints; 04 Dec 14;
   write(ounit,1001)            myid,                                                      "|II|e", alog10(IIe(1:min(Mvol-1,28)))
   endif
   endif
  endif

1000 format("xspech : ",f10.2," : myid=",i3," ; ":"|f|="es12.5" ; ":"time=",f10.2,"s ;":" log"a5,:"="28f6.2" ...")
1001 format("xspech : ", 10x ," : myid=",i3," ; ":"    "  12x "   ":"     ", 10x ,"  ;":" log"a5,:"="28f6.2" ...")

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! if( myid.eq.0 ) then ; cput = GETTIME ; write(ounit,'("xspech : ", 10x ," : ")') ; write(ounit,1000) cput-cpus, ForceErr, cput-lastcpu
! endif
  
!1000 format("xspech : ",f10.2," : "6x,3x" ":"|f|="es12.5" ; ":"time=",f10.2,"s ;")
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! compute loop current and plasma current;
  
!  do vvol = 1, Mvol
!   
!   if( myid.ne.modulo(vvol-1,ncpu) ) cycle
!   
!   if( .not.ImagneticOK(vvol) ) then
!    cput = GETTIME
!    write(ounit,1002) cput-cpus
!    write(ounit,1002) cput-cpus, myid, vvol, ImagneticOK(vvol)
!    cycle
!   endif
!   
!   if( Igeometry.eq.1 .or.  vvol.gt.1 ) Lcoordinatesingularity = .false. ! either Cartesian geometry or annular volume; used elsewhere;
!   if( Igeometry.gt.1 .and. vvol.eq.1 ) Lcoordinatesingularity = .true.
!   
!   if( vvol.le.Nvol ) Lplasmaregion = .true.
!   if( vvol.gt.Nvol ) Lplasmaregion = .false.
!   
!   Lvacuumregion = .not. Lplasmaregion
!   
!   WCALL(xspech, cu00aa,(vvol, Ntz, mn) )
!
!  enddo
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
  if( Lcheck.eq.5 .or. LHevalues .or. LHevectors .or. Lperturbed.eq.1 ) then ! check construction of Hessian; 01 Jul 14;
   
   if( myid.eq.0 ) then
    cput = GETTIME
    write(ounit,'("xspech : ", 10x ," : ")')
    write(ounit,'("xspech : ",f10.2," : myid=",i3," ; calling he01aa ; see .ext.hessian.myid ;")') cput-cpus, myid
   endif
   
   WCALL(xspech, he01aa,( Ngeometricaldof, position(0:Ngeometricaldof), Mvol, mn, lgeometricaldof ))
   
  endif ! end of if( Lcheck.eq.5 ) ; 01 Jul 14;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( Igeometry.eq.1 ) pflux(1) = dpflux(1) ! Cartesian; only in this case is poloidal flux in innermost volume defined; 04 Dec 14;
  if( Igeometry.gt.1 ) pflux(1) = zero
  
  do vvol = 2, Nvol
   
#ifdef DEBUG
   FATALMESS(xspech, vvol.gt.Nvol, dpflux is not defined)
#endif
   
   pflux(vvol) = pflux(vvol-1) + dpflux(vvol) ! 01 Jul 14;
   
  enddo
  
  FATALMESS(xspech, abs(phiedge).lt.vsmall, phiedge is not valid for normalization)
  
  pflux(1:Nvol) = pflux(1:Nvol) * pi2 / phiedge ! normalize; 02 Sep 14;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! WCALL(xspech,ra00aa,('W')) ! this writes vector potential to file;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  Lcontinueiterations = .false. ! this will be set in bn00aa; 18 Oct 12;
  
  if( Lfreebound.gt.0 .and. maxfbits.gt.0 ) then
 !if( Lfreebound.gt.0 ) then
   
   WCALL(xspech,bn00aa,( mn, Ntz )) ! compute normal field etc. on computational boundary; ! this will set Lcontinueiterations; 03 Apr 13;
   
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  WCALL(xspech,ra00aa,('W')) ! this writes vector potential to file;

  if( myid.eq.0 ) then ! WRITE RESTART FILE; note that this is inside free-boundary iteration loop; 11 Aug 14;

   wflag=0 ; iflag=0 ; rflag=zero

   WCALL(xspech,writin,( wflag, iflag, rflag )) ! write restart file; save initial input;

  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Lcontinueiterations .and. Lfindzero.gt.0 ) goto 9000 

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! FREE-BOUNDARY ITERATIONS HAVE FINISHED;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( myid.eq.0 ) then ! this is just screen diagnostics; 20 Jun 14;

   cput = GETTIME
   
   if( nPpts.gt.0 ) then
    write(ounit,'("xspech : ", 10x ," :")')
    write(ounit,'("xspech : ",f10.2," : myid=",i3," ; Poincare plot ; odetol="es8.1" ; nPpts="i7" ;":" nPtrj="24(i5",")" ...")') &
  cput-cpus, myid, odetol, nPpts, nPtrj(1:min(Mvol,24))
   endif
   
   if( Lcheck.eq.1 ) then
    write(ounit,'("xspech : ", 10x ," :")')
    write(ounit,'("xspech : ",f10.2," : myid=",i3," ; calling jo00aa; computing error in field ; epsr="es8.1" ;")') cput-cpus, myid, epsr
   endif

  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  do vvol = 1, Mvol ! various diagnostic calculations are performed in parallel; 03 Apr 13;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
      
   if( myid.ne.modulo(vvol-1,ncpu) ) cycle ! the following is in parallel; 20 Jun 14;
   
#ifdef DEBUG
   FATALMESS(xspech, .not.allocated(ImagneticOK), error)
#endif

   if( .not.ImagneticOK(vvol) ) then
    cput = GETTIME
    write(ounit,1002) cput-cpus
    write(ounit,1002) cput-cpus, myid, vvol, ImagneticOK(vvol)
    cycle
   endif
   
1002 format("xspech : ",f10.2," :":" myid=",i3," ; vvol=",i3," ; IBeltrami="L2" ; construction of Beltrami field failed ;")
   
  !if( Igeometry.eq.1 .or.  vvol.gt.1 ) Lcoordinatesingularity = .false. ! either Cartesian geometry or annular volume; used elsewhere;
  !if( Igeometry.gt.1 .and. vvol.eq.1 ) Lcoordinatesingularity = .true.
   
   if( Igeometry.eq.1 .or. vvol.gt.1 ) then ; Lcoordinatesingularity = .false. ! 14 Jan 15;
   else                                     ; Lcoordinatesingularity = .true.  ! 14 Jan 15;
   endif
   
   if( vvol.le.Nvol ) Lplasmaregion = .true.
   if( vvol.gt.Nvol ) Lplasmaregion = .false.
   
   Lvacuumregion = .not. Lplasmaregion


   WCALL(xspech,sc00aa,( vvol, Ntz )) ! compute covariant field at interfaces; related to singular currents; 20 Jun 14; need to broadcast Bsubtmn, Bsubzmn;

   
   if( Lcheck.eq.1 ) then
    lquad = Iquad(vvol)*2 ! the minimal/optimal quadrature resolution needs revision; 20 Feb 13;
    WCALL(xspech,jo00aa,( vvol, Ntz, Lrad(vvol), lquad, mn ))    
   endif
   
   
   if( Mpqits.gt.0 .and. npq(vvol).gt.0 ) then ! locate periodic orbits in each volume; 03 Apr 13;
    WCALL(xspech,pq01aa,( vvol ))
    FATALMESS(xspech, pqorbit(vvol,1)%ok.ne.1, periodic orbit failed)
    pqt(vvol) = pqorbit(vvol,1)%to
    pqs(vvol) = pqorbit(vvol,1)%so
   endif
    
   if( Lvacuumregion ) then ! this will be deleted after construction of scalar potential has been fully de-bugged etc.; 03 Apr 13;

    normalerror(0:1) = zero ! initialize summation; 17 Apr 13;
    ivol = vvol ; Ltangent = 0 ! these must be passed through to bf00aa; 17 Apr 13;
    do innout = 0, 0
     lss = two * innout - one
     do kk = 0, Nz-1 ; zeta = kk * pi2nfp / Nz 
      do jj = 0, Nt-1 ; teta = jj * pi2    / Nt ; jk = 1 + jj + kk*Nt
       st(1:6) = (/ lss, teta, one, zero, zero, one /)
       WCALL(xspech, bf00aa,( zeta, st(1:6), Bst(1:6) ))
       normalerror(innout) = normalerror(innout) + abs(Bst(1)) ! confirm that interfaces are flux surfaces; 03 Apr 13;
      enddo ! end of do kk; 24 Apr 13;
     enddo ! end of do jj; 24 Apr 13;
    enddo ! end of do innout; 24 Apr 13;
    cput = GETTIME
    write(ounit,'("xspech : ", 10x ," : ")')
    write(ounit,'("xspech : ",f10.2," : myid=",i3," ; vvol=",i3," ; normal error="2es13.5" ;")') cput-cpus, myid, vvol, normalerror(0:1)/(Nt*Nz)

   endif
    
   if( nPpts.gt.0 ) then ! construct Poincare plot in each volume; 03 Apr 13;
    WCALL(xspech,pp00aa,( vvol )) ! Poincare plots in each volume
   endif ! end of if( nPpts.gt.0 );
   
!  if( Lwrpj ) then ! pressure-jump Hamiltonian; 11 Aug 13;
!    WCALL(xspech,ph01aa,( vvol ))
!   endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  enddo ! end of do vvol = 1, Mvol; ! end of parallel diagnostics loop; 03 Apr 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  do vvol = 1, Mvol ! some belated broadcasts; 11 Aug 14;
   
   llmodnp = modulo(vvol-1,ncpu)
   
   RlBCAST(pqs(vvol),1,llmodnp)
   RlBCAST(pqt(vvol),1,llmodnp)   
   
#ifdef DEBUG
   FATALMESS(xspech, .not.allocated(Bsubtemn), error)
   FATALMESS(xspech, .not.allocated(Bsubzemn), error)
   FATALMESS(xspech, .not.allocated(Bsubtomn), error)
   FATALMESS(xspech, .not.allocated(Bsubzomn), error)
#endif
   
   RlBCAST(Bsubtemn(1:mn,0:1,vvol),mn*2,llmodnp)
   RlBCAST(Bsubzemn(1:mn,0:1,vvol),mn*2,llmodnp)
   RlBCAST(Bsubtomn(1:mn,0:1,vvol),mn*2,llmodnp)
   RlBCAST(Bsubzomn(1:mn,0:1,vvol),mn*2,llmodnp)
   
  enddo ! end of do vvol = 1, Mvol; 01 Jul 14;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( myid.eq.0 ) then ! WRITE RESTART FILE;

   wflag = 0 ; iflag = 0 ; rflag = zero

   WCALL(xspech,writin,( wflag, iflag, rflag )) ! write restart file; save initial input;
   
   close(zunit) ! this file is written to in globals/writin; 11 Aug 14;
   
   cput = GETTIME
   write(ounit,'("xspech : ", 10x ," :")')
   write(ounit,'("xspech : ",f10.2," : myid=",i3," : restart file closed ; time="f8.2"m = "f6.2"h = "f5.2"d ;")') &
cput-cpus, myid, (cput-cpus) / (/ 60, 60*60, 24*60*60 /)

  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
9999 continue

  if( Llatex ) then
  write(wunit,'("\input{{end}}")')
  close(wunit)
  endif

  WCALL(xspech,finish)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  stop
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end program xspech

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Closes output files, writes screen summary.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine finish
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  use constants, only : zero

  use fileunits, only : ounit

  use inputlist, only : Wmacros, Wxspech, ext, Ltiming, Nvol

  use cputiming

  use allglobal, only : myid, cpus, mn

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  REAL      :: Ttotal, dcpu, ecpu
  CHARACTER :: date*8, time*10
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  cpui = GETTIME ; cpuo = cpui ! see macro expansion for begin; 11 Aug 14;

#ifdef DEBUG
  if( Wxspech ) write(ounit,'("finish : ",f10.2," : myid=",i3," ; start  ;")') cpui-cpus, myid
#endif
  
  cput = GETTIME
  
! write(ounit,'("finish : ", 10x ," : before sumtime ;")')

! SUMTIME !  this is expanded by Makefile, and then again by macros; do not remove;

! write(ounit,'("finish : ", 10x ," : after  sumtime ;")')

  cput = GETTIME ; dcpu = cput-cpus
  
  if( Ltiming .and. myid.eq.0 ) then

   Ttotal = zero
   
! PRTTIME ! this is expanded by Makefile, and then again by macros; do not remove;
   write(ounit,'("finish : ",f10.2," : time spent in writin =",f10.2," ;")') dcpu, Twritin ; Ttotal = Ttotal + Twritin
   write(ounit,'("finish : ",f10.2," : time spent in readin =",f10.2," ;")') dcpu, Treadin ; Ttotal = Ttotal + Treadin

   ecpu = Ttotal-dcpu ! error in actual cpu time and calculated cpu time;  7 Mar 13; 

   write(ounit,'("finish : ",f10.2," : Ttotal =",f10.2," s = "f8.2" m = "f6.2" h ; Timing Error = ",f10.2,"s = ",f10.2,"%")') &
dcpu, Ttotal / (/ 1, 60, 3600 /), ecpu, 100*ecpu/dcpu

  endif ! end of if( Ltiming .and. myid.eq.0 ) then; 01 Jul 14;
  
  if( myid.eq.0 ) then
   
   call date_and_time(date,time)
   write(ounit,'("finish : ", 10x ," : ")')
   write(ounit,1000) dcpu, myid, dcpu / (/ 1, 60, 60*60, 24*60*60 /), date(1:4), date(5:6), date(7:8), time(1:2), time(3:4), time(5:6), ext
   write(ounit,'("finish : ", 10x ," : ")')

  !write(ounit,'("finish : ", 10x ," : calling hdfint ;")')
   if( myid.eq.0 ) then
   call hdfint ! 18 Jul 14;
   endif
  !write(ounit,'("finish : ", 10x ," : called  hdfint ;")')

  endif ! end of if( myid.eq.0 ) ; 14 Jan 15;

  MPIFINALIZE
  
1000 format("finish : ",f10.2," : myid=",i3," ; completion ; time=",f10.2,"s = "f8.2"m = "f6.2"h = "f5.2"d ; date= "&
  a4"/"a2"/"a2" ; time= "a2":"a2":"a2" ; ext = "a60)

  stop

end subroutine finish

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
