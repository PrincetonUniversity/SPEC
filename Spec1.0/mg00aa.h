!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Reads mgrid file, in binary format, and constructs ezspline interpolation.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine mg00aa
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half, pi2

  use numerical, only : vsmall, small

  use fileunits, only : ounit, gunit

  use inputlist, only : Wmacros, Wmg00aa, maxgroups, extcur, ext

  use cputiming, only : Tmg00aa

  use allglobal, only : myid, ncpu, cpus, Lmgridexist, &
                        nextcur, Rmin, Zmin, Rmax, Zmax, &
!                       oRZp, &
                        Lmgridhasbeensplined ! 11 Oct 12; 

#ifdef NOEZSPLINE
#else

  use ezmgrid ! only required for free-boundary; interpolation of mgrid; 28 Nov 12;

#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS

  LOGICAL              :: exist

  INTEGER              :: NR, NZ, NP, NPp, isHermite, ii, jj, kk, bcx(1:2), bcy(1:2), bcz(1:2), iextcur
  REAL                 :: dR, dphi, dZ, rzp(1:3), xyz(1:3), Bxyz(1:3), dBxyzdxyz(1:3,1:3), RR, ZZ, phi
  REAL                 :: dRBR, dZBR, dpBR, dRBZ, dZBZ, dpBZ, dRBp, dZBp, dpBp ! 11 Oct 12; for computing divergence and curl of mgrid field;
  REAL, allocatable    :: Brzp(:,:,:,:), TBrzp(:,:,:,:)
  REAL, allocatable    :: Jrzp(:,:,:,:) ! 11 Oct 12; for computing curl of mgrid field;
  CHARACTER*30         :: curlabel(maxgroups) ! ! 22 Apr 13;

  INTEGER              :: NN, LDFJAC, LWA, ic05pbf
  REAL                 :: Xsu(1:2), FVEC(1:2), FJAC(1:2,1:2), XTOL, WA(2*(2+13)/2)

! external             :: FCN

  BEGIN(mg00aa)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  Lmgridhasbeensplined = .false.

  Rmin = zero ; Rmax = zero ; Zmin = zero ; Zmax = zero !! default values; 03 Apr 13;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef NOEZSPLINE
#else

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  inquire(file="mgrid",exist=exist) ! check if file exists;

  if( .not.exist .and. myid.eq.0 ) then 
   cput = GETTIME
   write(ounit,'("mg00aa : ", 10x ," : ")') 
   write(ounit,'("mg00aa : ",f10.2," : mgrid file does not exist ; Rmin, Zmin, Rmax, & Zmax left as default ;")') cput-cpus
  endif
  
  if( .not.exist ) then ; Lmgridexist = .false. ; goto 9999
  else                  ; Lmgridexist = .true.
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item The read format of the (binary) mgrid file is as follows:
!latex       \begin{verbatim}
!latex       open(unit,file="mgrid",status='old',action='read',form='unformatted',iostat=ios)
!latex       read(unit,iostat=ios) NR, NZ, NP, nfp, nextcur ! integers;
!latex       read(unit,iostat=ios) Rmin, Zmin, Rmax, Zmax   ! doubles;
!latex       read(unit,iostat=ios) curlabel(1:nextcur)      ! character*30
!latex       do iextcur = 1, nextcur
!latex        read(unit,iostat=ios)Brzp(1:3,1:NR,1:NZ,1:NP) ! doubles;
!latex       enddo
!latex       close(unit)
!latex       \end{verbatim}
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( myid.eq.0 ) then ! open mgrid file;
   
   open(gunit,file="mgrid",status='old',action='read',form='unformatted',iostat=ios)
   FATALMESS(mg00aa,ios.ne.0,error opening mgrid)
   
   read(gunit,iostat=ios) NR, NZ, NP, mgridnfp, nextcur
   FATALMESS(mg00aa,ios.ne.0,error reading NR NZ NP mgridnfp nextcur)

   read(gunit,iostat=ios) Rmin, Zmin, Rmax, Zmax
   FATALMESS(mg00aa,ios.ne.0,error reading Rmin Zmin Rmax Zmax)
   
   cput = GETTIME
   write(ounit,'("mg00aa : ", 10x ," : ")')
   write(ounit,1000) cput-cpus, NR, NZ, NP, mgridnfp, nextcur !! 17 Apr 13;
   write(ounit,1001) cput-cpus, Rmin, Rmax, Zmin, Zmax !! 17 Apr 13;
   
1000 format("mg00aa : ",f10.2," : NR="i5" , NZ="i5" , NP="i4" , mgridnfp=",i3," , nextcur=",i3," ;")
1001 format("mg00aa : ",f10.2," : [Rmin,Rmax]=["es13.5" ,"es13.5" ] ; [Zmin,Zmax]=["es13.5" ,"es13.5" ] ;")
   
  endif ! end of if( myid.eq.0 ) then;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  IlBCAST(NR,1,0)
  IlBCAST(NZ,1,0)
  IlBCAST(NP,1,0)
  IlBCAST(mgridnfp,1,0)
  IlBCAST(nextcur,1,0)
  
  FATALMESS(mg00aa, nextcur.gt.maxgroups, need to increase maxgroups )
  
  NPp = NP + 1 ; pi2mgridnfp = pi2 / mgridnfp ! shorthand;
  
  RlBCAST(Rmin,1,0)
  RlBCAST(Zmin,1,0)
  RlBCAST(Rmax,1,0)
  RlBCAST(Zmax,1,0)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RALLOCATE(TBrzp,(1:3,1:NR,1:NZ,1:NPp)) ! 11 Oct 12; total magnetic field; allocated on each myid;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( myid.eq.0 ) then ! read field due to each independent current; sum to determine total field;
   
   read(gunit,iostat=ios) curlabel(1:nextcur)
   FATALMESS(mg00aa,ios.ne.0,error reading curlabel)
   
   RALLOCATE(Brzp,(1:3,1:NR,1:NZ,1:NP)) ! 11 Oct 12; temporary array for reading from mgrid file; only allocated on myid=0;
   
   do iextcur = 1, nextcur
    
    read(gunit,iostat=ios) Brzp(1:3,1:NR,1:NZ,1:NP)
    FATALMESS(mg00aa,ios.ne.0,error reading Brzp)
    
    if( Wmg00aa ) then ; cput = GETTIME ; write(ounit,1010)cput-cpus, myid, iextcur, extcur(iextcur), iextcur, curlabel(iextcur)
    endif
    
1010 format("mg00aa : ",f10.2," : myid=",i3," ; extcur[",i3," ] ="es23.15" ; curlabel[",i3," ] = "a30)

    TBrzp(1:3,1:NR,1:NZ,1:NP) = TBrzp(1:3,1:NR,1:NZ,1:NP) + extcur(iextcur) * Brzp(1:3,1:NR,1:NZ,1:NP) ! total magnetic field due to external coils;

   enddo ! end of do iextcur = 1,nextcur loop;
   
!  DEALLOCATE(Brzp) ! 11 Oct 12; shall use this later

   close(gunit)

  endif ! end of if( myid.eq.0 ) then;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  FATALMESS(mg00aa, sum(abs(extcur(1:nextcur))).lt.small, supplied currents are trivial ) ! 11 Oct 12; 

!  if( sum(abs(extcur(1:nextcur))).lt.small ) then ! currents are trivial; equivalent to no externally supplied field; ! 11 Oct 12; 

!   Lmgridexist = .false. ! 11 Oct 12; 
   
!   if( myid.eq.0 ) write(ounit,'("mg00aa : ", 10x ," : myid=",i3," ; no external currents supplied ; herafter assume Lmgridexist=F ;")')myid ! 11 Oct 12; 
!   if( myid.eq.0 ) then ! 11 Oct 12; 
!    DEALLOCATE(Brzp) ! 11 Oct 12; 
!   endif ! 11 Oct 12; 
!   DEALLOCATE(TBrzp) ! 11 Oct 12; 
!   goto 9999 ! 11 Oct 12; 
!    ! 11 Oct 12; 
!  endif ! 11 Oct 12; 

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RlBCAST(TBrzp(1:3,1:NR,1:NZ,1:NP),3*NR*NZ*NP,0)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  dR = (Rmax-Rmin) / (NR-1) ; dphi = pi2mgridnfp / NP ; dZ = (Zmax-Zmin) / (NZ-1)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! DON'T DELETE THE BELOW: IT WILL EVENTUALLY BE REQUIRED TO INVERT GIVEN [R,Z,PHI] TO CONSTRUCT [\S,\T,\Z];
!
!   ic05pbf = -1 ; NN = 2 ; LDFJAC = NN ; XTOL = vsmall ; LWA = NN * (NN+13)/2
!   
!     rzp(1:3) = . . . ! selected point;
!     
!     if( ic05pbf.ne.0 ) Xsu(1:2) = (/ half, zero /) ! initial guess;
!     
!     ic05pbf = 1
!     call C05PBF( FCN, NN, Xsu(1:NN), FVEC(1:NN), FJAC(1:LDFJAC,1:NN), LDFJAC, XTOL, WA(1:LWA), LWA, ic05pbf ) ! determine whether (R,Z,p) point lies inside plasma boundary;
!
!     cput = GETTIME !                                                               1         2  
!     select case( ic05pbf ) !                                                  1234567890123456789
!     case(   0 )  ;               write(ounit,1020) myid, oRZp(1:3), ic05pbf, "success            " ,Xsu(1:2)
!     case( :-1 )  ;               write(ounit,1020) myid, oRZp(1:3), ic05pbf, "user terminated    "!,Xsu(1:2)
!     case(   1 )  ;               write(ounit,1020) myid, oRZp(1:3), ic05pbf, "input error        "!,Xsu(1:2)
!     case(   2 )  ;               write(ounit,1020) myid, oRZp(1:3), ic05pbf, "too many iterations"!,Xsu(1:2)
!     case(   3 )  ;               write(ounit,1020) myid, oRZp(1:3), ic05pbf, "xtol too small     "!,Xsu(1:2)
!     case(   4 )  ;               write(ounit,1020) myid, oRZp(1:3), ic05pbf, "no solution        "!,Xsu(1:2)
!     case default ;               write(ounit,1020) myid, oRZp(1:3), ic05pbf, "illegal ifail      "!,Xsu(1:2)
!     end select
!
!1020 format("mg00aa : "10x" : myid=",i3," ; (R,Z,p)=("es13.5" ,"es13.5" ,"es13.5" ) ; ic05pbf="i2" ; "a19" ;":" (s,t)=("f15.10" ,"f15.10" ) ;")
!     
! DON'T DELETE THE ABOVE: IT WILL EVENTUALLY BE REQUIRED TO INVERT GIVEN [R,Z,PHI] TO CONSTRUCT [\S,\T,\Z];

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  TBrzp(1:3,1:NR,1:NZ,NPp) = TBrzp(1:3,1:NR,1:NZ,1) ! include additional plane for spline completeness;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Wmg00aa .and. myid.eq.0 ) then ; cput = GETTIME ; write(ounit,'("mg00aa : ",f10.2," : myid=",i3," ; constructing splines ;")')cput-cpus,myid; cpui = cput
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  bcx(1:2)= (/ 0, 0/) ; bcy(1:2)= (/ 0, 0/) ; bcz(1:2)= (/-1,-1/) ; ierr=0 ! set not-a-knot, periodic, periodic boundary conditions
  
  call ezspline_init( sTBr, NR, NZ, NPp, bcx, bcy, bcz, ierr ) ; call ezspline_error(ierr)
  FATALMESS(mg00aa,ierr.ne.0,spline error init Br)

  call ezspline_init( sTBz, NR, NZ, NPp, bcx, bcy, bcz, ierr ) ; call ezspline_error(ierr)
  FATALMESS(mg00aa,ierr.ne.0,spline error init Bz)

  call ezspline_init( sTBp, NR, NZ, NPp, bcx, bcy, bcz, ierr ) ; call ezspline_error(ierr)
  FATALMESS(mg00aa,ierr.ne.0,spline error init Bp)
  
  isHermite=0 ; sTBr%isHermite=isHermite ; sTBz%isHermite=isHermite ; sTBp%isHermite=isHermite
  
  do ii = 1, NR  ; sTBr%x1(ii) = Rmin + (ii-1)*dR   ; sTBz%x1(ii) = Rmin + (ii-1)*dR   ; sTBp%x1(ii) = Rmin + (ii-1)*dR   ! need to reset  R  range to non-standard domain;
  enddo
  do jj = 1, NZ  ; sTBr%x2(jj) = Zmin + (jj-1)*dZ   ; sTBz%x2(jj) = Zmin + (jj-1)*dZ   ; sTBp%x2(jj) = Zmin + (jj-1)*dZ   ! need to reset  Z  range to non-standard domain;
  enddo
  do kk = 1, NPp ; sTBr%x3(kk) = zero + (kk-1)*dphi ; sTBz%x3(kk) = zero + (kk-1)*dphi ; sTBp%x3(kk) = zero + (kk-1)*dphi ! need to reset phi range to non-standard domain;
  enddo
  
  call ezspline_setup( sTBr, TBrzp(1,:,:,:), ierr ) ; call ezspline_error(ierr)
  FATALMESS(mg00aa,ierr.ne.0,spline error setup Br)

  call ezspline_setup( sTBz, TBrzp(2,:,:,:), ierr ) ; call ezspline_error(ierr)
  FATALMESS(mg00aa,ierr.ne.0,spline error setup Bz)

  call ezspline_setup( sTBp, TBrzp(3,:,:,:), ierr ) ; call ezspline_error(ierr)
  FATALMESS(mg00aa,ierr.ne.0,spline error setup Bp)
  
  if( Wmg00aa .and. myid.eq.0 ) then ; cput = GETTIME ; write(ounit,'("mg00aa : ",f10.2," : myid=",i3," ; constructed  splines ; took ",f10.2,"s ;")')cput-cpus,myid,cput-cpui
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  Lmgridhasbeensplined = .true. ! 11 Oct 12; 

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( myid.eq.0 ) then ! 11 Oct 12; evaluate grid divergence; ! 11 Oct 12; calculate curl of vacuum field;

   RALLOCATE(Jrzp,(1:3,1:NR,1:NZ,1:NP)) ! 11 Oct 12; temporary array for reading from mgrid file; only allocated on myid=0;
   
   Brzp(1:3,1:NR,1:NZ,1:NP) = zero ! 11 Oct 12; recall that Brzp is only allocated on myid=0;

   do ii = 1, NR ; RR  = Rmin + (ii-1) * dR
    do jj = 1, NZ ; ZZ  = Zmin + (jj-1) * dZ
     do kk = 1, NP ; phi = zero + (kk-1) * dphi
      
      call ezspline_derivative( sTBr, 1, 0, 0, RR, ZZ, phi, dRBR, ierr )
      FATALMESS(mg00aa,ierr.ne.0,spline error R derivative Br)
      call ezspline_derivative( sTBr, 0, 1, 0, RR, ZZ, phi, dZBR, ierr )
      FATALMESS(mg00aa,ierr.ne.0,spline error Z derivative Br)
      call ezspline_derivative( sTBr, 0, 0, 1, RR, ZZ, phi, dpBR, ierr )
      FATALMESS(mg00aa,ierr.ne.0,spline error p derivative Br)

      call ezspline_derivative( sTBz, 1, 0, 0, RR, ZZ, phi, dRBZ, ierr )
      FATALMESS(mg00aa,ierr.ne.0,spline error R derivative Bz)
      call ezspline_derivative( sTBz, 0, 1, 0, RR, ZZ, phi, dZBZ, ierr )
      FATALMESS(mg00aa,ierr.ne.0,spline error Z derivative Bz)
      call ezspline_derivative( sTBz, 0, 0, 1, RR, ZZ, phi, dpBZ, ierr )
      FATALMESS(mg00aa,ierr.ne.0,spline error p derivative Bz)

      call ezspline_derivative( sTBp, 1, 0, 0, RR, ZZ, phi, dRBp, ierr )
      FATALMESS(mg00aa,ierr.ne.0,spline error R derivative Bp)
      call ezspline_derivative( sTBp, 0, 1, 0, RR, ZZ, phi, dZBp, ierr )
      FATALMESS(mg00aa,ierr.ne.0,spline error Z derivative Bp)
      call ezspline_derivative( sTBp, 0, 0, 1, RR, ZZ, phi, dpBp, ierr )
      FATALMESS(mg00aa,ierr.ne.0,spline error p derivative Bp)
      
      Brzp(1,ii,jj,kk) = abs( RR*dRBR + RR*dZBZ + dpBp ) / ( TBrzp(1,ii,jj,kk)**2+TBrzp(1,ii,jj,kk)**2+TBrzp(1,ii,jj,kk)**2 ) ! 11 Oct 12; normalized divergence error;

      Jrzp(1,ii,jj,kk) = (                        dpBZ - RR * dZBp ) / RR ! 11 Oct 12; curl B . \hat R
      Jrzp(2,ii,jj,kk) = ( TBrzp(3,ii,jj,kk) + RR*dRBp -      dpBR ) / RR ! 11 Oct 12; curl B . \hat Z
      Jrzp(3,ii,jj,kk) = (                        dZBR -      dRBZ )      ! 11 Oct 12; curl B . \hat p
      
     enddo
    enddo
   enddo
   
   Brzp(1,1:NR,1:NZ,1:NP) = Brzp(1,1:NR,1:NZ,1:NP) / (NR*NZ*NP)
      
   open(gunit+myid,file="."//trim(ext)//".mgriddivergence",status="unknown",form="unformatted") ! 11 Oct 12; 
   write(gunit+myid) Rmin, Rmax, Zmin, Zmax
   write(gunit+myid) NR, NZ, NP
   write(gunit+myid) Brzp(1,1:NR,1:NZ,1:NP)
   write(gunit+myid) Jrzp
   close(gunit+myid)                                                                            ! 11 Oct 12; 

   DEALLOCATE(Jrzp) ! 11 Oct 12; 

   DEALLOCATE(Brzp)
   
!   if( Wmg00aa ) then ; cput = GETTIME ; write(ounit,'("mg00aa : ",f10.2," : myid=",i3," ; constructed divergence error, pseudo-mgrid ;")')cput-cpus,myid ! 11 Oct 12; 
!   endif ! 11 Oct 12; 
   
  endif ! end of if( myid.eq.0 ) ;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  DEALLOCATE(TBrzp) ! 11 Oct 12; only the spline objects sTBR, sTBz and sTBp are required hereafter;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(mg00aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine mg00aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine mgridfield( rzp, Brzp, dBrzpdrzp, imgridfail ) ! the argument list is arbitrary; this routine is NOT called by NAG;
  
  use constants, only : zero, one

  use numerical, only : small

  use fileunits, only : ounit

  use inputlist, only : Wmg00aa

  use cputiming, only : Tmg00aa

  use allglobal, only : myid, ncpu, cpus, Rmin, Zmin, Rmax, Zmax, Ltangent, Lmgridhasbeensplined

#ifdef NOEZSPLINE
#else

  use ezmgrid

#endif
  
  LOCALS

  REAL   , intent(in)  ::  rzp(1:3)
  REAL   , intent(out) :: Brzp(1:3), dBrzpdrzp(1:3,1:3)
  INTEGER, intent(out) :: imgridfail

  REAL              :: phi

  BEGIN(mg00aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef NOEZSPLINE
#else
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  Brzp(1:3) = (/ zero, zero, one /) ; dBrzpdrzp(1:3,1:3) = zero ! default intent out;

#ifdef DEBUG
  if( .not.Lmgridhasbeensplined                                                  ) then ; imgridfail = 2 ; goto 9999 ! 11 Oct 12; 
  endif
  if( rzp(1).lt.Rmin .or. rzp(1).gt.Rmax .or. rzp(2).lt.Zmin .or. rzp(2).gt.Zmax ) then ; imgridfail = 1 ; goto 9999
  endif
#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  FATALMESS(mg00aa, abs(pi2mgridnfp).lt.small, this is a problem ) !! can place this inside #ifdef DEBUG; 28 Apr 13;

  phi = modulo( rzp(3), pi2mgridnfp )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  call ezspline_interp( sTBr, rzp(1), rzp(2), phi, Brzp(1), ierr ) ; call ezspline_error(ierr)
  FATALMESS(mg00aa,ierr.ne.0,spline error interp Br)

  call ezspline_interp( sTBz, rzp(1), rzp(2), phi, Brzp(2), ierr ) ; call ezspline_error(ierr)
  FATALMESS(mg00aa,ierr.ne.0,spline error interp Bz)

  call ezspline_interp( sTBp, rzp(1), rzp(2), phi, Brzp(3), ierr ) ; call ezspline_error(ierr)
  FATALMESS(mg00aa,ierr.ne.0,spline error interp Bp)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \item Note that ${\bf B} \cdot \hat \phi \; \hat \phi = {\bf B} \cdot \hat \phi \; R^{-1} \; {\bf e}_\phi$ implies that $B^\phi = {\bf B} \cdot \hat \phi / R$.

  Brzp(3) = Brzp(3) / rzp(1) ! cylindrical contravariant field; 12 Oct 12;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  select case( Ltangent ) ! determines whether derivatives are required;
   
  case( 0 )
   
   dBrzpdrzp(1:3,1:3) = zero 
   
  case( 1 )

   call ezspline_derivative( sTBr, 1, 0, 0, rzp(1), rzp(2), phi, dBrzpdrzp(1,1), ierr ) ; call ezspline_error(ierr)
   FATALMESS(mg00aa,ierr.ne.0,spline error derivative Br)

   call ezspline_derivative( sTBr, 0, 1, 0, rzp(1), rzp(2), phi, dBrzpdrzp(1,2), ierr ) ; call ezspline_error(ierr)
   FATALMESS(mg00aa,ierr.ne.0,spline error derivative Br)

   call ezspline_derivative( sTBr, 0, 0, 1, rzp(1), rzp(2), phi, dBrzpdrzp(1,3), ierr ) ; call ezspline_error(ierr)
   FATALMESS(mg00aa,ierr.ne.0,spline error derivative Br)


   call ezspline_derivative( sTBz, 1, 0, 0, rzp(1), rzp(2), phi, dBrzpdrzp(2,1), ierr ) ; call ezspline_error(ierr)
   FATALMESS(mg00aa,ierr.ne.0,spline error derivative Bz)

   call ezspline_derivative( sTBz, 0, 1, 0, rzp(1), rzp(2), phi, dBrzpdrzp(2,2), ierr ) ; call ezspline_error(ierr)
   FATALMESS(mg00aa,ierr.ne.0,spline error derivative Bz)

   call ezspline_derivative( sTBz, 0, 0, 1, rzp(1), rzp(2), phi, dBrzpdrzp(2,3), ierr ) ; call ezspline_error(ierr)
   FATALMESS(mg00aa,ierr.ne.0,spline error derivative Bz)


   call ezspline_derivative( sTBp, 1, 0, 0, rzp(1), rzp(2), phi, dBrzpdrzp(3,1), ierr ) ; call ezspline_error(ierr)
   FATALMESS(mg00aa,ierr.ne.0,spline error derivative Bp)

   call ezspline_derivative( sTBp, 0, 1, 0, rzp(1), rzp(2), phi, dBrzpdrzp(3,2), ierr ) ; call ezspline_error(ierr)
   FATALMESS(mg00aa,ierr.ne.0,spline error derivative Bp)

   call ezspline_derivative( sTBp, 0, 0, 1, rzp(1), rzp(2), phi, dBrzpdrzp(3,3), ierr ) ; call ezspline_error(ierr)
   FATALMESS(mg00aa,ierr.ne.0,spline error derivative Bp)

   
  case default
   
   FATALMESS(fb00aa,.true.,invalid Ltangent)
   
  end select
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  imgridfail = 0

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(mg00aa)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine mgridfield

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! DON'T DELETE THE BELOW: IT WILL EVENTUALLY BE REQUIRED TO INVERT GIVEN [R,Z,PHI] TO CONSTRUCT [\S,\T,\Z];
!
!subroutine FCN( NN, Xsu, FVEC, FJAC, LDFJAC, IFLAG )
!
!  use constants, only : zero, half, one, pi2
!  use numerical, only : small
!  use fileunits, only : ounit, gunit
!  use inputlist, only : Wmg00aa, Nvol
!  use cputiming, only : Tmg00aa
!  use allglobal, only : myid, ncpu, cpus, mn, im, in, interfacelabel, oRZp
!
!  LOCALS
!  
!  INTEGER :: NN, LDFJAC, IFLAG
!  REAL    :: Xsu(1:NN), FVEC(1:NN), FJAC(1:LDFJAC,1:NN)
!  
!  INTEGER :: imn, ivol, lvol
!  REAL    :: arg, carg, sarg, stz(1:3), RpZ(1:3), dR(0:3), dZ(0:3), jacobian, guv(1:3,1:3)
!  
!  BEGIN(mg00aa)
!
!  stz(1:3) = (/ min( max(Xsu(1),small), one ), Xsu(2), oRZp(3) /)
!
!  do ivol = 1, Nvol
!   if( Xsu(1).ge.interfacelabel(ivol-1) .and. Xsu(1).le.interfacelabel(ivol) ) lvol = ivol
!  enddo
!  
!  call co00aa( lvol, stz(1:3), RpZ(1:3), dR(0:3), dZ(0:3), jacobian, guv(1:3,1:3) )
!  
!  select case( iflag )
!  case( 1 )    ; FVEC(1:2) = (/ RpZ(1), RpZ(3) /) - oRZp(1:2) ! update FVEC;
!  case( 2 )    ; FJAC(1,1:2) = dR(1:2)
!                 FJAC(2,1:2) = dZ(1:2)
!  case default ; FATALMESS(mg00aa,.true.,invalid iflag)
!  end select
! 
!  RETURN(mg00aa)
!
!end subroutine FCN
!
! DON'T DELETE THE ABOVE: IT WILL EVENTUALLY BE REQUIRED TO INVERT GIVEN [R,Z,PHI] TO CONSTRUCT [\S,\T,\Z];

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
