!> \file
!> \brief Follows magnetic fieldline using ode-integration routine from rksuite.f .

!> \brief Constructs Poincaré plot and "approximate" rotational-transform (for single field line).
!> \ingroup grp_diagnostics
!>
!> **relevant input variables**
!>
!> <ul>
!> <li> The resolution of Poincaré plot is controlled by
!>      <ul>
!>      <li> \c nPpts iterations per trajectory; </li>
!>      <li> \c odetol o.d.e. integration tolerance; </li>
!>      </ul>
!>      The magnetic field is given by bfield() . </li>
!> </ul>
!>
!> **rotational-transform**
!>
!> <ul>
!> <li> The approximate rotational transform is determined by field line integration.
!>      This is constructed by fitting a least squares fit to the field line trajectory. </li>
!> </ul>
!>
!> @param[in]  lvol
!> @param      sti
!> @param[in]  Nz
!> @param[in]  nPpts
!> @param      poincaredata
!> @param      fittedtransform
!> @param[out] utflag
subroutine pp00ab( lvol, sti, Nz, nPpts, poincaredata, fittedtransform, utflag )
  use mod_kinds, only: wp => dp
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  use constants, only : zero, one, two, pi2

  use numerical, only : small

  use fileunits, only : ounit

  use inputlist, only : Wpp00ab, Nvol, odetol

  use cputiming, only : Tpp00ab

  use allglobal, only : myid, ncpu, cpus, MPI_COMM_SPEC, pi2nfp, ivol, Mvol, Node

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


#ifdef OPENMP
  USE OMP_LIB
#endif
  use mpi
  implicit none
  integer   :: ierr, astat, ios, nthreads, ithread
  real(wp)      :: cput, cpui, cpuo=0 ! cpu time; cpu initial; cpu old; 31 Jan 13;


  integer, intent(in)  :: lvol, Nz, nPpts
  integer, intent(out) :: utflag
  real(wp)                 :: sti(1:2), poincaredata(1:4,0:Nz-1,1:nPpts), fittedtransform(1:2), dzeta

  integer              :: jj, kk
  real(wp)                 :: ppt(1:4)

  integer, parameter   :: Lrwork = 20*Node
  real(wp)                 :: zst, zend, st(1:Node), rwork(1:Lrwork), tol, stz(1:3), RpZ(1:3), leastfit(1:5)
  character            :: RA

  integer, parameter   :: Lenwrk = 32*Node
  integer              :: rkmethod, outch
  real(wp)                 :: hstart, thres(1:Node), rkwork(1:Lenwrk), mchpes, dwarf
  real(wp)                 :: zgot, ygot(1:Node), ypgot(1:Node), ymax(1:Node)
  character            :: rktask
  LOGICAL              :: errass, mesage

  external             :: bfield
  external             :: SETUP, UT, ENVIRN


  cpui = MPI_WTIME()
  cpuo = cpui
#ifdef OPENMP
  nthreads = omp_get_max_threads()
#else
  nthreads = 1
#endif


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  ivol = lvol ! required to pass through to bfield;

#ifdef DEBUG

   if( lvol.lt.1.or.lvol.gt.Mvol ) then
     write(6,'("pp00ab :      fatal : myid=",i3," ; lvol.lt.1.or.lvol.gt.Mvol ; invalid volume ;")') myid
     call MPI_ABORT( MPI_COMM_SPEC, 1, ierr )
     stop "pp00ab : lvol.lt.1.or.lvol.gt.Mvol : invalid volume  ;"
    endif


   if( abs(sti(1)).gt.one ) then
     write(6,'("pp00ab :      fatal : myid=",i3," ; abs(sti(1)).gt.one ; illegal radial coordinate ;")') myid
     call MPI_ABORT( MPI_COMM_SPEC, 1, ierr )
     stop "pp00ab : abs(sti(1)).gt.one : illegal radial coordinate  ;"
    endif

#endif

  dzeta = pi2nfp

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RA = 'D' ; tol = odetol

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  st(1:Node) = (/ sti(1), sti(2) /) ! , one, zero, zero, one /)!, zero, zero /)

  ppt(1:2) = (/ sti(2), sti(1) /) ! seems redundant; just used for packing into poincaredata it seems;

  leastfit(1:5) = (/ zero, zero, zero, ppt(1), one /) ! initialize summation for least squares fit;

  utflag = 0 ; poincaredata(1:4,0:Nz-1,1:nPpts) = zero ; fittedtransform(1:2) = - two ! provide dummy defaults;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  do jj = 1, nPpts ! loop over iterations;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   zst = zero ! starting Poincare plane;

   call ENVIRN(outch,mchpes,dwarf) ! only dwarf is used to set thres=sqrt(dwarf); thres could be set larger but not smaller
   thres(1:Node) = sqrt(dwarf); rkmethod = 3; rktask = 'U'; errass = .FALSE. ; hstart = 0.0D0
#ifdef DEBUG
   mesage = .TRUE.
#else
   mesage = .FALSE.
#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   do kk = 0, Nz-1 ! loop over toroidal Poincare cross sections;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    stz(1:3) = (/ ppt(2), mod(ppt(1),pi2), zst /) ! toroidal coordinates;

    if( abs(stz(1)).gt.one ) then ;                                        ; utflag = 0 ; exit ! exit do kk loop; 22 Apr 13; ! 28 Feb 17;

    endif

1002 format("pp00ab : ", 10x ," : myid=",i3," ; lvol=",i3," ; "3x" : (s,t)=("f21.17" ,"f21.17" ) ;         "3x" ; outside domain ;")


   cput = MPI_WTIME()
   Tpp00ab = Tpp00ab + ( cput-cpuo )
   call stzxyz( lvol, stz(1:3), RpZ(1:3) )
   cpuo = MPI_WTIME()
 ! map to cylindrical;

    ppt(3:4)=(/ RpZ(1), RpZ(3) /) ! cylindrical coordinates;

    poincaredata(1:4,kk,jj) = (/ mod(ppt(1),pi2), ppt(2), ppt(3), ppt(4) /) ! save Poincare information to array;

    zend = zst + (pi2nfp/Nz)

    utflag = 0

    call SETUP(Node, zst, st(1:Node), zend, tol, thres(1:Node), rkmethod, rktask, errass, hstart, rkwork(1:Lenwrk), Lenwrk, mesage)


   cput = MPI_WTIME()
   Tpp00ab = Tpp00ab + ( cput-cpuo )
   call UT(bfield, zend, zgot, ygot(1:Node), ypgot(1:Node), ymax(1:Node), rkwork(1:Lenwrk), utflag)
   cpuo = MPI_WTIME()
 ! integrate to next plane;

    zst = zend

    st(1:Node) = ygot(1:Node)

    cput = MPI_WTIME()
    select case( utflag )                                                !         1         2         3         4         5         6
    case( 1 ) ; ! give screen output if error is encountered;             !123456789012345678901234567890123456789012345678901234567890123
    case( 2 ) ; write(ounit,2001) cput-cpus, myid, lvol, jj, kk, utflag, "step size too small (try RK method 2)                          "
    case( 3 ) ; write(ounit,2001) cput-cpus, myid, lvol, jj, kk, utflag, "integration interrupted (more than 5000 calls)                 "
    case( 4 ) ; write(ounit,2001) cput-cpus, myid, lvol, jj, kk, utflag, "problem is stiff (UT is not efficient)                         "
    case( 5 ) ; write(ounit,2001) cput-cpus, myid, lvol, jj, kk, utflag, "odetol/thres too small or RK method too low                    "
    case( 6 ) ; write(ounit,2001) cput-cpus, myid, lvol, jj, kk, utflag, "integration interrupted (error assessment not possible         "
    case default

   if( .true. ) then
     write(6,'("pp00ab :      fatal : myid=",i3," ; .true. ; illegal value of ifail returned from UT;")') myid
     call MPI_ABORT( MPI_COMM_SPEC, 1, ierr )
     stop "pp00ab : .true. : illegal value of ifail returned from UT ;"
    endif

    end select

2001 format("pp00ab : ",f10.2," : myid=",i3," ; lvol=",i3," ; (jj,kk)=("i4" ,"i4" ); ifail="i2" ; "a63)

    if( utflag.ne.1 ) exit ! an integration error was encountered; exit do kk loop;

    ppt(1:2) = (/ st(2), st(1) /) ! again, seems redundant; I think this is only required for packing into poincaredata;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   enddo ! end of do kk;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   if( utflag.ne.1 ) exit ! an integration error was encountered; exit do jj loop;

   leastfit(1:5) = leastfit(1:5) + (/ (jj*dzeta)**2, (jj*dzeta), (jj*dzeta)*ppt(1), ppt(1), one /) ! least squares fit summation;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  enddo ! end of do jj = 1, nPpts;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( utflag.eq.1 ) then
   fittedtransform(1:2) = (/ sti(1), ( leastfit(5)*leastfit(3)-leastfit(2)*leastfit(4) ) / ( leastfit(5)*leastfit(1)-leastfit(2)*leastfit(2) ) /)
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


9999 continue
  cput = MPI_WTIME()
  Tpp00ab = Tpp00ab + ( cput-cpuo )
  return


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine pp00ab

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
