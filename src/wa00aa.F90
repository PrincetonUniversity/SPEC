!> \defgroup grp_smooth_boundary Smooth boundary
!>
!> \file
!> \brief Constructs smooth approximation to wall.

!> \brief ...todo...
module laplaces
  use mod_kinds, only: wp => dp
  LOGICAL              :: stage1        !< what is this ?
  LOGICAL              :: exterior      !< what is this ?
  LOGICAL              :: dorm          !< what is this ?
  integer              :: Nintervals    !< what is this ?
  integer              :: Nsegments     !< what is this ?
  integer              :: IC            !< what is this ?
  integer              :: NP4           !< what is this ?
  integer              :: NP1           !< what is this ?
  integer, allocatable :: icint(:)      !< what is this ?
  real(wp)                 :: originalalpha !< what is this ?
  real(wp), allocatable    :: xpoly(:)      !< what is this ?
  real(wp), allocatable    :: ypoly(:)      !< what is this ?
  real(wp), allocatable    :: phi(:)        !< what is this ?
  real(wp), allocatable    :: phid(:)       !< what is this ?
  real(wp), allocatable    :: CC(:,:)       !< what is this ?

  integer              :: ilength       !< what is this ?
  real(wp)                 :: totallength   !< what is this ?

  integer              :: niterations   !< counter; eventually redundant; 24 Oct 12;

  integer              :: iangle        !< angle ; eventually redundant; 24 Oct 12;

  real(wp)                 :: Rmid          !< used to define local polar coordinate; eventually redundant; 24 Oct 12;
  real(wp)                 :: Zmid          !< used to define local polar coordinate; eventually redundant; 24 Oct 12;

  real(wp)                 :: alpha         !< eventually redundant; 24 Oct 12;

end module laplaces

!> \brief Constructs smooth approximation to wall.
!> \ingroup grp_smooth_boundary
!>
!> **solution of Laplace's equation in two-dimensions**
!>
!> <ul>
!> <li> The wall is given by a discrete set of points. </li>
!> <li> The points must go anti-clockwise. </li>
!> </ul>
!>
subroutine wa00aa( iwa00aa )
  use mod_kinds, only: wp => dp
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  use constants, only : zero, half, one, ten, pi2

  use numerical, only : vsmall

  use fileunits, only : ounit, gunit

  use inputlist, only : Wmacros, Wwa00aa, Nvol, Mpol, Ntor, odetol

  use cputiming, only : Twa00aa

  use allglobal, only : ncpu, myid, cpus, &
                        Mvol, &
                        mn, im, in, iRbc, iZbs, iRbs, iZbc, &
                        Nt, Nz, Ntz, Rij, Zij, &
                        Lcoordinatesingularity, &
                        YESstellsym

  use laplaces

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


#ifdef OPENMP
  USE OMP_LIB
#endif
  use mpi
  implicit none
  integer   :: ierr, astat, ios, nthreads, ithread
  real(wp)      :: cput, cpui, cpuo=0 ! cpu time; cpu initial; cpu old; 31 Jan 13;


   integer, parameter   :: Nconstraints = 1, lengthwork = Nconstraints * ( 3*Nconstraints + 13 ) / 2 ! required for C05NBF;

   integer              :: iwa00aa, Lcurvature, Nwall, iwall, ii, ifail

   real(wp)                 :: lss, lRZ(1:2), px, py, Rmin, Rmax, Zmin, Zmax, xtol
   real(wp)                 :: rho(1:Nconstraints), fvec(1:Nconstraints), realwork(1:lengthwork)

   real(wp), allocatable    :: RZwall(:,:)

  !REAL                 :: phiwall

   external             :: VacuumPhi


  cpui = MPI_WTIME()
  cpuo = cpui
#ifdef OPENMP
  nthreads = omp_get_max_threads()
#else
  nthreads = 1
#endif


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef DEBUG

   if( myid.ne.0 ) then
     write(6,'("wa00aa :      fatal : myid=",i3," ; myid.ne.0 ; error ;")') myid
     call MPI_ABORT( MPI_COMM_SPEC, 1, ierr )
     stop "wa00aa : myid.ne.0 : error  ;"
    endif


   if( Ntor.gt.0 ) then
     write(6,'("wa00aa :      fatal : myid=",i3," ; Ntor.gt.0 ; presently axisymmetry is assumed but this can easily be generalized ;")') myid
     call MPI_ABORT( MPI_COMM_SPEC, 1, ierr )
     stop "wa00aa : Ntor.gt.0 : presently axisymmetry is assumed but this can easily be generalized  ;"
    endif

#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  lss = one ; Lcurvature = 0 ; Lcoordinatesingularity = .false.


   cput = MPI_WTIME()
   Twa00aa = Twa00aa + ( cput-cpuo )
   call co01aa( Nvol, lss, Lcurvature, Ntz, mn )
   cpuo = MPI_WTIME()
 ! get plasma boundary, which serves as inner boundary; 10 Apr 13;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  open(gunit, file="wall.dat", status='old', action='read', iostat=ios ) ! read polygon, which serves as outer boundary; 10 Apr 13;

   if( ios.ne.0 ) then
     write(6,'("wa00aa :      fatal : myid=",i3," ; ios.ne.0 ; error opening wall.dat ;")') myid
     call MPI_ABORT( MPI_COMM_SPEC, 1, ierr )
     stop "wa00aa : ios.ne.0 : error opening wall.dat  ;"
    endif


  Nwall = 0
  do
   read(gunit,*,iostat=ios)lRZ(1:2)
   if( ios.ne.0 ) exit
   Nwall = Nwall + 1
  enddo
  close(gunit)


   allocate( RZwall(1:2,1:Nwall), stat=astat )
   RZwall(1:2,1:Nwall) = zero


  open(gunit,file="wall.dat",status='old',action='read',iostat=ios)

   if( ios.ne.0 ) then
     write(6,'("wa00aa :      fatal : myid=",i3," ; ios.ne.0 ; error opening wall.dat ;")') myid
     call MPI_ABORT( MPI_COMM_SPEC, 1, ierr )
     stop "wa00aa : ios.ne.0 : error opening wall.dat  ;"
    endif


  read(gunit,*,iostat=ios) RZwall(1:2,1:Nwall) ! MUST GO ANTI-CLOCKWISE; LAST-POINT = FIRST POINT;

   if( ios.ne.0 ) then
     write(6,'("wa00aa :      fatal : myid=",i3," ; ios.ne.0 ; error reading RZwall from wall.dat ;")') myid
     call MPI_ABORT( MPI_COMM_SPEC, 1, ierr )
     stop "wa00aa : ios.ne.0 : error reading RZwall from wall.dat  ;"
    endif


  close(gunit)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! Rmin = minval( RZwall(1,1:Nwall) )
  Rmax = maxval( RZwall(1,1:Nwall) )
! Zmin = minval( RZwall(2,1:Nwall) )
! Zmax = maxval( RZwall(2,1:Nwall) )
!
  Rmid = half * ( maxval(Rij(1:Ntz,0,0) ) + minval(Rij(1:Ntz,0,0) ) ) ! defined by plasma boundary; 10 Apr 13;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!                                          boundary points +    nodes  + repeated point + inner boundary + inner nodes
  Nintervals = Nwall-1 + Ntz ; Nsegments =       Nwall     + (Nwall-1) +         1      +    (Ntz+1)     +     Ntz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


   allocate( xpoly(1:Nsegments), stat=astat )
   xpoly(1:Nsegments) = zero


   allocate( ypoly(1:Nsegments), stat=astat )
   ypoly(1:Nsegments) = zero



   allocate( phi (1:Nintervals), stat=astat )
   phi (1:Nintervals) = zero
 ! must set to boundary value of phi;

   allocate( phid(1:Nintervals), stat=astat )
   phid(1:Nintervals) = zero
 ! can leave this as zero;

  IC = Nintervals + 1

  NP4 = Nintervals + 4
  NP1 = Nintervals + 1


   allocate( CC   (1:IC,1:NP4), stat=astat )
   CC   (1:IC,1:NP4) = zero


   allocate( ICINT(     1:NP1), stat=astat )
   ICINT(     1:NP1) = zero


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

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( Wwa00aa ) then
   cput = MPI_WTIME()
   do ii = 1, iwall ; write(ounit,'("wa00aa : ",f10.2," : Nsegments="i9" ; ii="i9" ; R="f15.9" ; Z="f15.9" ;")') cput-cpus, Nsegments, ii, xpoly(ii), ypoly(ii)
   enddo
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  stage1   = .true.  ! preparation;
 !exterior = .true.  ! unbounded domain;
  exterior = .false. ! interior  domain;
  dorm     = .true.  ! Dirichlet or mixed boundary value problem;

  px = xpoly(1)      ! only required for stage two
  py = ypoly(1)      ! only required for stage two

  alpha = zero       ! alpha need not be set if dorm = .true. ;

  ifail = 1

  call D03EAF( stage1, exterior, dorm, Nintervals, px, py, xpoly, ypoly, Nsegments, phi, phid, alpha, CC, IC, NP4, ICINT, NP1, ifail )

  cput = MPI_WTIME()
  ;         ; write(ounit,'("wa00aa : ", 10x ," : ")')
  select case( ifail )
  case( 0 ) ; write(ounit,'("wa00aa : ",f10.2," : prepared vacuum calculation; stage1="L2" ; ifail=",i3," ; success ;          ")') cput-cpus, stage1, ifail
  case( 1 ) ; write(ounit,'("wa00aa : ",f10.2," : prepared vacuum calculation; stage1="L2" ; ifail=",i3," ; invalid tolerance ;")') cput-cpus, stage1, ifail
  case( 2 ) ; write(ounit,'("wa00aa : ",f10.2," : prepared vacuum calculation; stage1="L2" ; ifail=",i3," ; incorrect rank ;   ")') cput-cpus, stage1, ifail
  case default

   if( .true. ) then
     write(6,'("wa00aa :      fatal : myid=",i3," ; .true. ; invalid ifail returned by D03EAF ;")') myid
     call MPI_ABORT( MPI_COMM_SPEC, 1, ierr )
     stop "wa00aa : .true. : invalid ifail returned by D03EAF  ;"
    endif

  end select

  stage1 = .false. ; originalalpha = alpha

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! the intialization of Laplaces solution is complete; 24 Oct 12;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  DEALLOCATE( RZwall )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  rho(1) = (Rmax + Rij(1,0,0)) * half - Rmid ! initial guess; 10 Apr 13;

  do iangle = 1, Ntz

   niterations = 0 ; xtol = odetol

   ifail = 1

   call C05NBF( VacuumPhi, Nconstraints, rho(1:Nconstraints), fvec(1:Nconstraints), xtol, realwork(1:lengthwork), lengthwork, ifail )

   cput = MPI_WTIME()
   select case( ifail )
   case( :-1 ) ;               write(ounit,'("wa00aa : ",f10.2," : iangle="i6" ; ifail=",i3," ; outside domain ; FAIL ;   ")') cput-cpus, iangle, ifail
   case(   0 ) ; if( Wwa00aa ) write(ounit,'("wa00aa : ",f10.2," : iangle="i6" ; ifail=",i3," ; success ;                 ")') cput-cpus, iangle, ifail
   case(   1 ) ;               write(ounit,'("wa00aa : ",f10.2," : iangle="i6" ; ifail=",i3," ; input error ; FAIL ;      ")') cput-cpus, iangle, ifail
   case(   2 ) ;               write(ounit,'("wa00aa : ",f10.2," : iangle="i6" ; ifail=",i3," ; consider restart ; FAIL ; ")') cput-cpus, iangle, ifail
   case(   3 ) ;               write(ounit,'("wa00aa : ",f10.2," : iangle="i6" ; ifail=",i3," ; xtol is too small ; FAIL ;")') cput-cpus, iangle, ifail
   case(   4 ) ;               write(ounit,'("wa00aa : ",f10.2," : iangle="i6" ; ifail=",i3," ; bad progress ; FAIL ;     ")') cput-cpus, iangle, ifail
   case default

   if( .true. ) then
     write(6,'("wa00aa :      fatal : myid=",i3," ; .true. ; invalid ifail returned by C05NBF ;")') myid
     call MPI_ABORT( MPI_COMM_SPEC, 1, ierr )
     stop "wa00aa : .true. : invalid ifail returned by C05NBF  ;"
    endif

   end select

   if( ifail.ne.0 ) exit

  enddo ! end of do iangle = 1, Ntz loop;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  DEALLOCATE(xpoly)
  DEALLOCATE(ypoly)
  DEALLOCATE(phi)
  DEALLOCATE(phid)
  DEALLOCATE(CC)
  DEALLOCATE(ICINT)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  iwa00aa = ifail

  cput = MPI_WTIME()
  write(ounit,'("wa00aa : ",f10.2," : constructed outer boundary; iwa00aa=",i3," ;")') cput-cpus, iwa00aa ! 24 Oct 12;

  if( ifail.ne.0 ) goto 9999

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  call tfft( Nt, Nz, Rij(1:Ntz,0,0), Zij(1:Ntz,0,0), &
             mn, im(1:mn), in(1:mn), iRbc(1:mn,Mvol), iRbs(1:mn,Mvol), iZbc(1:mn,Mvol), iZbs(1:mn,Mvol), ifail )

  if( YESstellsym ) then ; iRbs(1:mn,Mvol)= zero ; iZbc(1:mn,Mvol) = zero
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


9999 continue
  cput = MPI_WTIME()
  Twa00aa = Twa00aa + ( cput-cpuo )
  return


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine wa00aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> \brief Compute vacuum magnetic scalar potential (?)
!> \ingroup grp_smooth_boundary
!>
!> @param Nconstraints
!> @param rho
!> @param fvec
!> @param iflag
subroutine VacuumPhi( Nconstraints, rho, fvec, iflag )
  use mod_kinds, only: wp => dp
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  use constants, only : zero, half, one, pi2

  use inputlist, only : Wmacros, Wwa00aa

  use fileunits, only : ounit

  use allglobal, only : ncpu, myid, cpus, Ntz, Rij, Zij

  use laplaces

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


#ifdef OPENMP
  USE OMP_LIB
#endif
  use mpi
  implicit none
  integer   :: ierr, astat, ios, nthreads, ithread
  real(wp)      :: cput, cpui, cpuo=0 ! cpu time; cpu initial; cpu old; 31 Jan 13;


  integer :: Nconstraints, iflag
  real(wp)    :: rho(1:Nconstraints), fvec(1:Nconstraints), angle, px, py

 !REAL    :: phiwall

  integer :: ifail

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  niterations = niterations + 1

  alpha = originalalpha ! this needs to be reset before every call to D03EAF; 24 Oct 12;

  angle = - (iangle-1) * pi2 / Ntz ; px = Rmid + rho(1) * cos( angle )
  ;                                  py =        rho(1) * sin( angle )

  Rij(iangle,0,0) = px
  Zij(iangle,0,0) = py ! global information; passed back to wa00aa; used in FFT to construct new boundary;

  ifail = 1

  call D03EAF( stage1, exterior, dorm, Nintervals, px, py, xpoly, ypoly, Nsegments, phi, phid, alpha, CC, IC, NP4, ICINT, NP1, ifail )

  cput = MPI_WTIME()

  select case( ifail )
  case( 0 ) ; if( Wwa00aa ) write(ounit,'("wa00aa : ",f10.2," : stage1="L2" ; ifail=",i3," ; success ;          ")') cput-cpus, stage1, ifail
  case( 1 ) ;               write(ounit,'("wa00aa : ",f10.2," : stage1="L2" ; ifail=",i3," ; invalid tolerance ;")') cput-cpus, stage1, ifail
  case( 2 ) ;               write(ounit,'("wa00aa : ",f10.2," : stage1="L2" ; ifail=",i3," ; incorrect rank ;   ")') cput-cpus, stage1, ifail
  case default

   if( .true. ) then
     write(6,'("wa00aa :      fatal : myid=",i3," ; .true. ; invalid ifail returned by D03EAF ;")') myid
     call MPI_ABORT( MPI_COMM_SPEC, 1, ierr )
     stop "wa00aa : .true. : invalid ifail returned by D03EAF  ;"
    endif

  end select

  if( rho(1).lt.zero ) iflag = -1 ! could also check that R, Z are within domain;

 !fvec(1) = alpha - phiwall
  fvec(1) = alpha - half

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( Wwa00aa ) then
   write(ounit,1000)iangle, niterations, rho(1), px, py, alpha, fvec(1)
  endif

1000 format("wa00aa : ", 10x ," : iangle="i6" ; nits=",i3," ; rho="es13.5" ; R,Z="2es23.15" ; alpha="es13.5" ; F="es13.5" ;")

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  return

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine VacuumPhi

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
