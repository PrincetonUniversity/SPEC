!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Constructs smooth approximation to wall.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \end{enumerate} \subsubsection{solution of Laplace's equation in two-dimensions} \begin{enumerate}

!latex \item The wall is given by a discrete set of points.
!latex \item The points must go anti-clockwise.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine wa00aa( iwa00aa )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  use constants, only : zero, half, one, ten, pi2

  use numerical, only : vsmall

  use fileunits, only : ounit, gunit

  use inputlist, only : Wmacros, Wwa00aa, Nvol, Mpol, Ntor, phiwall, odetol

  use cputiming, only : Twa00aa
  
  use allglobal, only : ncpu, myid, cpus, &
                        Mvol, &
                        mn, im, in, iRbc, iZbs, iRbs, iZbc, &
                        Nt, Nz, Ntz, Rij, Zij, trigwk, trigm, trign, isr, &
                        Lcoordinatesingularity, &
                        YESstellsym
  
  use laplaces

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS

   INTEGER, parameter   :: Nconstraints = 1, lengthwork = Nconstraints * ( 3*Nconstraints + 13 ) / 2 ! required for C05NBF;
  
   INTEGER              :: iwa00aa, Lcurvature, Nwall, iwall, ii, ifail

   REAL                 :: lss, lRZ(1:2), px, py, Rmin, Rmax, Zmin, Zmax, xtol
   REAL                 :: rho(1:Nconstraints), fvec(1:Nconstraints), realwork(1:lengthwork)

   REAL, allocatable    :: RZwall(:,:)

   external             :: VacuumPhi
  
  BEGIN(wa00aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  FATALMESS(wa00aa, myid.ne.0, error )
  FATALMESS(wa00aa, Ntor.gt.0, presently axisymmetry is assumed but this can easily be generalized )
#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  lss = one ; Lcurvature = 0 ; Lcoordinatesingularity = .false.
  
  WCALL(wa00aa,co01aa,( Nvol, lss, Lcurvature, Ntz, mn )) ! get plasma boundary, which serves as inner boundary; 10 Apr 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  open(gunit, file="wall.dat", status='old', action='read', iostat=ios ) ! read polygon, which serves as outer boundary; 10 Apr 13;
  FATALMESS(wa00aa, ios.ne.0, error opening wall.dat )
  
  Nwall = 0
  do
   read(gunit,*,iostat=ios)lRZ(1:2)
   if( ios.ne.0 ) exit
   Nwall = Nwall + 1
  enddo
  close(gunit)
  
  RALLOCATE(RZwall,(1:2,1:Nwall))
  
  open(gunit,file="wall.dat",status='old',action='read',iostat=ios)
  FATALMESS(wa00aa,ios.ne.0,error opening wall.dat)   
  
  read(gunit,*,iostat=ios) RZwall(1:2,1:Nwall) ! MUST GO ANTI-CLOCKWISE; LAST-POINT = FIRST POINT;
  FATALMESS(wa00aa,ios.ne.0,error reading RZwall from wall.dat)
  
  close(gunit)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! Rmin = minval( RZwall(1,1:Nwall) )
  Rmax = maxval( RZwall(1,1:Nwall) )
! Zmin = minval( RZwall(2,1:Nwall) )
! Zmax = maxval( RZwall(2,1:Nwall) )
!
  Rmid = half * ( maxval(Rij(1:Ntz,0,0) ) + minval(Rij(1:Ntz,0,0) ) ) ! defined by plasma boundary; 10 Apr 13;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!                                          boundary points +    nodes  + repeated point + inner boundary + inner nodes  
  Nintervals = Nwall-1 + Ntz ; Nsegments =       Nwall     + (Nwall-1) +         1      +    (Ntz+1)     +     Ntz
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RALLOCATE(xpoly,(1:Nsegments))
  RALLOCATE(ypoly,(1:Nsegments))
  
  RALLOCATE(phi,(1:Nintervals)) ! must set to boundary value of phi;
  RALLOCATE(phid,(1:Nintervals)) ! can leave this as zero;
  
  IC = Nintervals + 1
  
  NP4 = Nintervals + 4
  NP1 = Nintervals + 1
  
  RALLOCATE(CC,(1:IC,1:NP4))
  
  RALLOCATE(ICINT,(1:NP1))
  
  iwall = 0
  
! outer boundary must go anti-clockwise; 24 Oct 12;

   do ii = 1, Nwall-1
    iwall = iwall + 1 ; xpoly(iwall) =          RZwall(1,ii)                    ; ypoly(iwall) =          RZwall(2,ii)                    ! vertices; 24 Oct 12;
    iwall = iwall + 1 ; xpoly(iwall) = half * ( RZwall(1,ii) + RZwall(1,ii+1) ) ; ypoly(iwall) = half * ( RZwall(2,ii) + RZwall(2,ii+1) ) ! nodes   ; 24 Oct 12;
   enddo
      ii =    Nwall
    iwall = iwall + 1 ; xpoly(iwall) =          RZwall(1,ii)                    ; ypoly(iwall) =          RZwall(2,ii)                    ! last point=first; 24 Oct 12;
   
! repeated point indicates end of outer boundary; 24 Oct 12;

    iwall = iwall + 1 ; xpoly(iwall) =          RZwall(1,ii)                    ; ypoly(iwall) =          RZwall(2,ii)

! inner boundary must go clockwise; 24 Oct 12;
   
   do ii = 1, Ntz-1
    iwall = iwall + 1 ; xpoly(iwall) =          Rij( ii,0,0)                  ; ypoly(iwall) =          Zij( ii,0,0)                  ! vertices; 24 Oct 12;
    iwall = iwall + 1 ; xpoly(iwall) = half * ( Rij( ii,0,0)+Rij(ii+1,0,0 ) ) ; ypoly(iwall) = half * ( Zij( ii,0,0)+Zij(ii+1,0,0 ) ) ! nodes   ;
   enddo
    iwall = iwall + 1 ; xpoly(iwall) =          Rij(Ntz,0,0)                  ; ypoly(iwall) =          Zij(Ntz,0,0)                  ! vertices; 24 Oct 12;
    iwall = iwall + 1 ; xpoly(iwall) = half * ( Rij(Ntz,0,0)+Rij(   1,0,0) )  ; ypoly(iwall) = half * ( Zij(Ntz,0,0)+Zij(   1,0,0) )  ! nodes   ;
    iwall = iwall + 1 ; xpoly(iwall) =          Rij( 1,0,0)                   ; ypoly(iwall) =          Zij( 1,0,0)                   ! last point=first; 24 Oct 12;
    
    phi(      1:Nwall-1   ) = zero ! scalar potential on outer boundary; 24 Oct 12;
    phi(Nwall  :Nintervals) = one  ! scalar potential on inner boundary; 24 Oct 12;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( Wwa00aa ) then
   cput = GETTIME
   do ii = 1, iwall ; write(ounit,'("wa00aa : ",f10.2," : Nsegments="i9" ; ii="i9" ; R="f15.9" ; Z="f15.9" ;")') cput-cpus, Nsegments, ii, xpoly(ii), ypoly(ii)
   enddo
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  stage1   = .true.  ! preparation;
 !exterior = .true.  ! unbounded domain;
  exterior = .false. ! interior  domain;
  dorm     = .true.  ! Dirichlet or mixed boundary value problem;
  
  px = xpoly(1)      ! only required for stage two
  py = ypoly(1)      ! only required for stage two
  
  alpha = zero       ! alpha need not be set if dorm = .true. ;
  
  ifail = 1

  call D03EAF( stage1, exterior, dorm, Nintervals, px, py, xpoly, ypoly, Nsegments, phi, phid, alpha, CC, IC, NP4, ICINT, NP1, ifail )
  
  cput = GETTIME
  ;         ; write(ounit,'("wa00aa : ", 10x ," : ")')
  select case( ifail )
  case( 0 ) ; write(ounit,'("wa00aa : ",f10.2," : prepared vacuum calculation; stage1="L2" ; ifail=",i3," ; success ;          ")') cput-cpus, stage1, ifail
  case( 1 ) ; write(ounit,'("wa00aa : ",f10.2," : prepared vacuum calculation; stage1="L2" ; ifail=",i3," ; invalid tolerance ;")') cput-cpus, stage1, ifail
  case( 2 ) ; write(ounit,'("wa00aa : ",f10.2," : prepared vacuum calculation; stage1="L2" ; ifail=",i3," ; incorrect rank ;   ")') cput-cpus, stage1, ifail
  case default
   FATALMESS(wa00aa,.true.,invalid ifail returned by D03EAF)
  end select
  
  stage1 = .false. ; originalalpha = alpha

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! the intialization of Laplaces solution is complete; 24 Oct 12;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  DEALLOCATE( RZwall )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  rho(1) = (Rmax + Rij(1,0,0)) * half - Rmid ! initial guess; 10 Apr 13;
  
  do iangle = 1, Ntz
   
   niterations = 0 ; xtol = odetol
   
   ifail = 1
   
   call C05NBF( VacuumPhi, Nconstraints, rho(1:Nconstraints), fvec(1:Nconstraints), xtol, realwork(1:lengthwork), lengthwork, ifail )
   
   cput = GETTIME
   select case( ifail )
   case( :-1 ) ;               write(ounit,'("wa00aa : ",f10.2," : iangle="i6" ; ifail=",i3," ; outside domain ; FAIL ;   ")') cput-cpus, iangle, ifail
   case(   0 ) ; if( Wwa00aa ) write(ounit,'("wa00aa : ",f10.2," : iangle="i6" ; ifail=",i3," ; success ;                 ")') cput-cpus, iangle, ifail
   case(   1 ) ;               write(ounit,'("wa00aa : ",f10.2," : iangle="i6" ; ifail=",i3," ; input error ; FAIL ;      ")') cput-cpus, iangle, ifail
   case(   2 ) ;               write(ounit,'("wa00aa : ",f10.2," : iangle="i6" ; ifail=",i3," ; consider restart ; FAIL ; ")') cput-cpus, iangle, ifail
   case(   3 ) ;               write(ounit,'("wa00aa : ",f10.2," : iangle="i6" ; ifail=",i3," ; xtol is too small ; FAIL ;")') cput-cpus, iangle, ifail
   case(   4 ) ;               write(ounit,'("wa00aa : ",f10.2," : iangle="i6" ; ifail=",i3," ; bad progress ; FAIL ;     ")') cput-cpus, iangle, ifail
   case default
    FATALMESS(wa00aa,.true.,invalid ifail returned by C05NBF)
   end select

   if( ifail.ne.0 ) exit

  enddo ! end of do iangle = 1, Ntz loop;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  DEALLOCATE(xpoly)
  DEALLOCATE(ypoly)
  DEALLOCATE(phi)
  DEALLOCATE(phid)
  DEALLOCATE(CC)
  DEALLOCATE(ICINT)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  iwa00aa = ifail

  cput = GETTIME
  write(ounit,'("wa00aa : ",f10.2," : constructed outer boundary; phiwall="f6.3" ; iwa00aa=",i3," ;")') cput-cpus, phiwall, iwa00aa ! 24 Oct 12;

  if( ifail.ne.0 ) goto 9999
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  call tfft( Nt, Nz, Rij(1:Ntz,0,0), Zij(1:Ntz,0,0), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz), &
             mn, im(1:mn), in(1:mn), iRbc(1:mn,Mvol), iRbs(1:mn,Mvol), iZbc(1:mn,Mvol), iZbs(1:mn,Mvol), ifail )

  if( YESstellsym ) then ; iRbs(1:mn,Mvol)= zero ; iZbc(1:mn,Mvol) = zero
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(wa00aa)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine wa00aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine VacuumPhi( Nconstraints, rho, fvec, iflag )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  use constants, only : zero, one, pi2
  
  use inputlist, only : Wmacros, Wwa00aa, phiwall

  use fileunits, only : ounit
  
  use allglobal, only : ncpu, myid, cpus, Ntz, Rij, Zij
  
  use laplaces

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS
  
  INTEGER :: Nconstraints, iflag
  REAL    :: rho(1:Nconstraints), fvec(1:Nconstraints), angle, px, py

  INTEGER :: ifail

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  niterations = niterations + 1

  alpha = originalalpha ! this needs to be reset before every call to D03EAF; 24 Oct 12;
  
  angle = - (iangle-1) * pi2 / Ntz ; px = Rmid + rho(1) * cos( angle )
  ;                                  py =        rho(1) * sin( angle )
  
  Rij(iangle,0,0) = px 
  Zij(iangle,0,0) = py ! global information; passed back to wa00aa; used in FFT to construct new boundary;
  
  ifail = 1
  
  call D03EAF( stage1, exterior, dorm, Nintervals, px, py, xpoly, ypoly, Nsegments, phi, phid, alpha, CC, IC, NP4, ICINT, NP1, ifail )
  
  cput = GETTIME

  select case( ifail )
  case( 0 ) ; if( Wwa00aa ) write(ounit,'("wa00aa : ",f10.2," : stage1="L2" ; ifail=",i3," ; success ;          ")') cput-cpus, stage1, ifail
  case( 1 ) ;               write(ounit,'("wa00aa : ",f10.2," : stage1="L2" ; ifail=",i3," ; invalid tolerance ;")') cput-cpus, stage1, ifail
  case( 2 ) ;               write(ounit,'("wa00aa : ",f10.2," : stage1="L2" ; ifail=",i3," ; incorrect rank ;   ")') cput-cpus, stage1, ifail
  case default
   FATALMESS(wa00aa,.true.,invalid ifail returned by D03EAF)
  end select
  
  if( rho(1).lt.zero ) iflag = -1 ! could also check that R, Z are within domain;

  fvec(1) = alpha - phiwall

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
 
  if( Wwa00aa ) then
   write(ounit,1000)iangle, niterations, rho(1), px, py, alpha, fvec(1)
  endif

1000 format("wa00aa : ", 10x ," : iangle="i6" ; nits=",i3," ; rho="es13.5" ; R,Z="2es23.15" ; alpha="es13.5" ; F="es13.5" ;")

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  return

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine VacuumPhi

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
