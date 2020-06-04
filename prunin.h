!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (diagnostic) ! Interface Pruning Algorithm

!latex \briefly{Interface Pruning Algorithm (determine if an interface should be removed or not)}

!latex \calledby{\link{}}
!latex \calls{\link{}}

!latex \tableofcontents

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \subsubsection{list of algorithm used} 

!latex \begin{enumerate}
  
!latex \item Compute the maximal Lyapunov exponent on the interface and just next to the interface
!latex \item Compute the Greene's average residue for periodic orbits approching the interface

!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine prunin( lvol )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  use constants, only : zero, one

  use fileunits, only : ounit
  
  use inputlist, only : Wprunin, nPpts, Lfreebound 
  
  use cputiming, only : Tprunin

  use allglobal, only : myid, ncpu, cpus, &
                        Lcoordinatesingularity, Mvol

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS

  INTEGER, parameter :: n_expansion = 8

  INTEGER, intent(in) :: lvol

  REAL :: s, ds, st(6), lya(2), residue

  INTEGER :: innout, pi, qi, isuccess

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  BEGIN(prunin)  

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!  

  ! part one, compute the lyapunov exponent

  ds = 0.01
  
  cput = GETTIME

  if (.not. Lcoordinatesingularity) then

    s = -one
    st = (/ s, zero, one, zero, zero, one/)   ! inner boundary

    ! compute the maximal Lyapunov exponent on the boundary
    call lyapunov(lya, lvol, st(1:4), nPpts)
    write(ounit,1001) cput-cpus, myid, lvol, s, nPpts/2, nPpts, lya(1), lya(2)

    ! compute the Lyapunov exponent for s = -0.99
    s = -one + ds
    st = (/ s, zero, one, zero, zero, one/)  ! slightly away from inner boundary
    call lyapunov(lya, lvol, st(1:4), nPpts)
    write(ounit,1001) cput-cpus, myid, lvol, s, nPpts/2, nPpts, lya(1), lya(2)

    ! compute the Greene's residue for periodic orbits approaching the interface
    innout = 0
    call greenes_residue(residue, lvol, innout, n_expansion, pi, qi, isuccess)
    write(ounit,1002) cput-cpus, myid, lvol, innout, pi, qi, residue

  endif

  if (lvol .ne. Mvol .or. Lfreebound .eq. 0) then

    ! compute the maximal Lyapunov exponent on the boundary
    s = one
    st = (/ s, zero, -one, zero, zero, one/)   ! outer boundary
    call lyapunov(lya, lvol, st(1:4), nPpts)
    write(ounit,1001) cput-cpus, myid, lvol, s, nPpts/2, nPpts, lya(1), lya(2)

    ! compute the Lyapunov exponent for s = 0.99
    s = one - ds
    st = (/ s, zero, -one, zero, zero, one/)  ! slightly away from outer boundary
    call lyapunov(lya, lvol, st(1:4), nPpts)
    write(ounit,1001) cput-cpus, myid, lvol, s, nPpts/2, nPpts, lya(1), lya(2)

    ! compute the Greene's residue for periodic orbits approaching the interface
    innout = 1
    call greenes_residue(residue, lvol, innout, n_expansion, pi, qi, isuccess)
    write(ounit,1002) cput-cpus, myid, lvol, innout, pi, qi, residue
  endif

1001 format("prunin : ",f10.2," : myid=",i3," ; lvol=",i3," ; s =",es8.1," ; steps=", i4, ' ', i4 ," Lyapunov exponent = ", es13.5, " ", es13.5," ;")
1002 format("prunin : ",f10.2," : myid=",i3," ; lvol=",i3," ; innout =",i2," ; p/q=", i4,'/',i4 ,"; Greene's residue = ", es13.5," ;")

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!  
  
  RETURN(prunin)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine prunin

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine lyapunov(lya, lvol, sti, nPpts) 
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, one, two, pi2
  
  use numerical, only : small
  
  use fileunits, only : ounit
  
  use inputlist, only : Wpp00ab, Nvol, odetol
  
  use cputiming, only : Tprunin
  
  use allglobal, only : myid, ncpu, cpus, pi2nfp, ivol, Mvol, Node
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in)  :: lvol, nPpts
  REAL, intent(out) :: lya(1:2)
  REAL, intent(in)  :: sti(1:4)

  INTEGER              :: jj, kk
  
  INTEGER, parameter   :: Lrwork = 60*Node
  REAL                 :: zst, zend, st(1:3*Node), rwork(1:Lrwork), di
  CHARACTER            :: RA
  
  INTEGER, parameter   :: Lenwrk = 3*32*Node
  INTEGER              :: rkmethod, outch, utflag
  REAL                 :: hstart, thres(1:3*Node), rkwork(1:Lenwrk), mchpes, dwarf, dzeta, tol
  REAL                 :: zgot, ygot(1:3*Node), ypgot(1:3*Node), ymax(1:3*Node)
  CHARACTER            :: rktask
  LOGICAL              :: errass, mesage
 
  external             :: bfield_tangent
  external             :: SETUP, UT, ENVIRN

  BEGIN(prunin)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  ivol = lvol ! required to pass through to bfield;
  
  dzeta = pi2nfp
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RA = 'D' ; tol = odetol
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  st(1:4) = sti
  st(5) = zero
  st(6) = zero

  utflag = 0 ; 

  lya = zero
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  do jj = 1, nPpts ! loop over iterations;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
    zst = zero ! starting Poincare plane;
   
    call ENVIRN(outch,mchpes,dwarf) ! only dwarf is used to set thres=sqrt(dwarf); thres could be set larger but not smaller
    thres(1:3*Node) = sqrt(dwarf); rkmethod = 3; rktask = 'U'; errass = .FALSE. ; hstart = 0.0D0
#ifdef DEBUG
    mesage = .TRUE.
#else
    mesage = .FALSE.
#endif
   

    zend = zst + dzeta
        
    utflag = 0
    
    call SETUP(3*Node, zst, st, zend, tol, thres, rkmethod, rktask, errass, hstart, rkwork, Lenwrk, mesage) 
    
    CALL( prunin, UT, (bfield_tangent, zend, zgot, ygot, ypgot, ymax, rkwork, utflag) ) ! integrate to next plane;
    
    zst = zend
    
    st(1:2) = ygot(1:2)

    di = sqrt(ygot(3)**2 + ygot(4)**2)
    
    if (jj .le. nPpts / 2) then
      lya(1) = lya(1) + log(di)
    endif

    lya(2) = lya(2) + log(di)

    st(3:4) = ygot(3:4) / di

    cput = GETTIME
    select case( utflag )                                                !         1         2         3         4         5         6
    case( 1 ) ; ! give screen output if error is encountered;             !123456789012345678901234567890123456789012345678901234567890123
    case( 2 ) ; write(ounit,2001) cput-cpus, myid, lvol, jj, kk, utflag, "step size too small (try RK method 2)                          "
    case( 3 ) ; write(ounit,2001) cput-cpus, myid, lvol, jj, kk, utflag, "integration interrupted (more than 5000 calls)                 "
    case( 4 ) ; write(ounit,2001) cput-cpus, myid, lvol, jj, kk, utflag, "problem is stiff (UT is not efficient)                         "
    case( 5 ) ; write(ounit,2001) cput-cpus, myid, lvol, jj, kk, utflag, "odetol/thres too small or RK method too low                    "
    case( 6 ) ; write(ounit,2001) cput-cpus, myid, lvol, jj, kk, utflag, "integration interrupted (error assessment not possible         "
    case default
     FATAL(pp00ab,.true.,illegal value of ifail returned from UT)
    end select    
    
2001 format("prunin : ",f10.2," : myid=",i3," ; lvol=",i3," ; (jj,kk)=("i4" ,"i4" ); ifail="i2" ; "a63)
    
    if( utflag.ne.1 ) exit ! an integration error was encountered; exit do kk loop;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
 
  enddo ! end of do jj = 1, nPpts;
  
  lya(2) = lya(2) / float(nPpts) / dzeta

  lya(1) = lya(1) / float(nPpts / 2) / dzeta

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(prunin)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine lyapunov

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
subroutine greenes_residue(residue, lvol, innout, nexp, piout, qiout, isuccess) 
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, one, two, four, pi2
  
  use numerical, only : small
  
  use fileunits, only : ounit
  
  use inputlist, only : Wpp00ab, Nvol, odetol, iota, oita, Nfp, Igeometry
  
  use cputiming, only : Tprunin
  
  use allglobal, only : myid, ncpu, cpus, pi2nfp, ivol, Mvol, Node
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, parameter ::  max_ai = 1000, maxiter = 10 ! total number of iteration
  REAL, parameter :: fptol = 1e-9 ! the tolerance of finding the periodic orbit
  REAL, parameter :: thetastart = zero ! the angle to start the search
  REAL, parameter :: zetastart = zero ! the angle to start the search

  INTEGER, intent(in)  :: lvol, nexp, innout
  REAL, intent(out) :: residue
  INTEGER, intent(out) :: piout, qiout, isuccess

  INTEGER              :: ii, jj, kk
  
  INTEGER, parameter   :: Lrwork = 60*Node
  REAL                 :: zst, zend, st(1:3*Node), stini(1:3*Node), rwork(1:Lrwork), di, sumlndi(2)
  CHARACTER            :: RA
  
  INTEGER, parameter   :: Lenwrk = 3*32*Node
  INTEGER              :: rkmethod, outch, utflag
  REAL                 :: hstart, thres(1:3*Node), rkwork(1:Lenwrk), mchpes, dwarf, dzeta, tol
  REAL                 :: zgot, ygot(1:3*Node), ypgot(1:3*Node), ymax(1:3*Node)
  CHARACTER            :: rktask
  LOGICAL              :: errass, mesage
 
  external             :: bfield_tangent
  external             :: SETUP, UT, ENVIRN

  INTEGER :: ai(nexp), pi(nexp), qi(nexp), nterms
  REAL    :: iniota, outiota, transform, poverq, ds, df

  BEGIN(prunin)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  ! the iota of the inner surface and outer surface of lvol
  outiota = iota(lvol) / float(Nfp)
  iniota = oita(lvol-1) / float(Nfp)

  if (lvol .eq. 1 .and. Igeometry.ne.1) iniota = 0

  ds = 1e-6

  ! choose the interface according to innout
  if (innout .eq. 0) then
    transform = iniota
    stini = (/-one+ds, thetastart, one, zero, zero, one/)
  else
    transform = outiota
    stini = (/one-ds, thetastart, one, zero, zero, one/)
  endif

  ! obtain its continued fraction expansion
  call continued_fraction_expansion(transform, nexp, max_ai, nterms, ai, pi, qi)

  poverq = float(pi(nterms))/float(qi(nterms))

  ! check if the approximation is outside the volume, if yes, change pi and qi
  if ((iniota - poverq) * (outiota - poverq) .ge. zero) then
    qiout = qi(nterms-1)
    piout = pi(nterms-1)
  else
    qiout = qi(nterms)
    piout = pi(nterms)
  endif 
  poverq = float(piout) / float(qiout)

  ivol = lvol ! required to pass through to bfield;
  
  ! the length of integration, should be 2pi * p / NFP
  dzeta = pi2nfp
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RA = 'D' ; tol = odetol

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  utflag = 0 ; 

  residue = zero
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  do ii = 1, maxiter
   
    zst = zetastart ! starting Poincare plane;
    st = stini

    do jj = 1, qiout
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
      call ENVIRN(outch,mchpes,dwarf) ! only dwarf is used to set thres=sqrt(dwarf); thres could be set larger but not smaller
      thres(1:3*Node) = sqrt(dwarf); rkmethod = 3; rktask = 'U'; errass = .FALSE. ; hstart = 0.0D0

      zend = zst + dzeta
          
      utflag = 0
      
      call SETUP(3*Node, zst, st, zend, tol, thres, rkmethod, rktask, errass, hstart, rkwork, Lenwrk, mesage) 
      
      CALL( prunin, UT, (bfield_tangent, zend, zgot, ygot, ypgot, ymax, rkwork, utflag) ) ! integrate to next plane;
      
      zst = zend
      st = ygot

      cput = GETTIME
      select case( utflag )                                                !         1         2         3         4         5         6
      case( 1 ) ; ! give screen output if error is encountered;             !123456789012345678901234567890123456789012345678901234567890123
      case( 2 ) ; write(ounit,2001) cput-cpus, myid, lvol, jj, kk, utflag, "step size too small (try RK method 2)                          "
      case( 3 ) ; write(ounit,2001) cput-cpus, myid, lvol, jj, kk, utflag, "integration interrupted (more than 5000 calls)                 "
      case( 4 ) ; write(ounit,2001) cput-cpus, myid, lvol, jj, kk, utflag, "problem is stiff (UT is not efficient)                         "
      case( 5 ) ; write(ounit,2001) cput-cpus, myid, lvol, jj, kk, utflag, "odetol/thres too small or RK method too low                    "
      case( 6 ) ; write(ounit,2001) cput-cpus, myid, lvol, jj, kk, utflag, "integration interrupted (error assessment not possible         "
      case default
      FATAL(pp00ab,.true.,illegal value of ifail returned from UT)
      end select    
      
  2001 format("prunin : ",f10.2," : myid=",i3," ; lvol=",i3," ; (jj,kk)=("i4" ,"i4" ); ifail="i2" ; "a63)
      
      if( utflag.ne.1 ) exit ! an integration error was encountered; exit do kk loop;
    
    end do ! do jj = 1, piout
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
    df = st(2) - pi2 * float(piout) - stini(2)

    stini(1) = stini(1) - df / st(4)

    if (abs(df) < fptol) exit

  enddo ! end of do ii = 1, maxiter;
  

  if (abs(df) < fptol) then
    isuccess = 1
    residue = abs(2 - st(3) - st(6))**(1/float(qiout))
  else
    isuccess = 0
    residue = 0
  endif
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(prunin)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine greenes_residue

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine continued_fraction_expansion(irrational, n, thres_ai, nterms, ai, pi, qi)

! compute the continued fraction expansion for the first n terms, also generating p_i and q_i

  use constants, only : one
  use numerical, only : small

  implicit none

  REAL, intent(in) :: irrational  ! the irrational number to expand
  INTEGER, intent(in) :: n        ! the number of terms in the expansion
  INTEGER, intent(in) :: thres_ai ! the threshold for ai, if ai > thres_ai, the expansion will be considered as terminated
  INTEGER, intent(out) :: nterms   ! the final number of terms, if < n, the expansion terminates before n terms
  INTEGER, intent(out) :: ai(n), pi(n), qi(n) ! the ai terms in the expansion, and the corresponding pi/qi

  INTEGER :: ii, jj, numer, denom, tmpi
  REAL :: residue

  ai = 0; pi = 0; qi = 0

  residue = abs(irrational)

  do ii = 1, n
    nterms = ii

    ai(ii) = FLOOR(residue)

    if (ai(ii) > thres_ai) then
      ! ai is beyond threshold, the expansion has terminated
      ai(ii) = 0
      nterms = ii - 1
      exit
    end if

    residue = residue - ai(ii)

    if (abs(residue) < small) then
      ! residue is too small, the expansion has terminated
      exit
    end if
    residue = one / residue 

  end do

  ! construct pi / qi from ai

  do ii = 1, nterms
    
    numer = ai(ii)
    denom = 1

    do jj = ii - 1, 1, -1

      ! flip denom and numer
      tmpi = denom
      denom = numer
      numer = tmpi

      numer = numer + ai(jj) * denom
    enddo

    pi(ii) = numer
    qi(ii) = denom

  enddo

end subroutine continued_fraction_expansion
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
