!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (diagnostic) ! Follows magnetic fieldline using NAG ode-integration routine, D02BJF.

!latex \briefly{Constructs \Poincare plot and ``approximate'' rotational-transform (for single field line).}

!latex \calledby{\link{}}
!latex \calls{\link{}}

!latex \tableofcontents

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \subsubsection{relevant input variables} 

!latex \begin{enumerate}
  
!latex \item The resolution of \Poincare plot is controlled by 
!latex \begin{itemize} 
!latex \item \verb+nPpts+ iterations per trajectory;
!latex \item \verb+odetol+ o.d.e. integration tolerance;
!latex \end{itemize}
!latex The magnetic field is given by \verb+bfield+. 
!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{rotational-transform} 

!latex \begin{enumerate}
  
!latex \item The approximate rotational transform is determined by field line integration.
!latex       This is constructed by fitting a least squares fit to the field line trajectory.
!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine pp00ab( lvol, sti, Nz, nPpts, poincaredata, fittedtransform, utflag ) 
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, one, two, pi2
  
  use numerical, only : small
  
  use fileunits, only : ounit
  
  use inputlist, only : Wpp00ab, Nvol, odetol
  
  use cputiming, only : Tpp00ab
  
  use allglobal, only : myid, ncpu, cpus, pi2nfp, ivol, Mvol, Node
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in)  :: lvol, Nz, nPpts
  INTEGER, intent(out) :: utflag
  REAL                 :: sti(1:2), poincaredata(1:4,0:Nz-1,1:nPpts), fittedtransform(1:2), dzeta
  
  INTEGER              :: jj, kk
  REAL                 :: ppt(1:4)
  
  INTEGER, parameter   :: Lrwork = 20*Node
  REAL                 :: zst, zend, st(1:Node), rwork(1:Lrwork), tol, stz(1:3), RpZ(1:3), leastfit(1:5)
  CHARACTER            :: RA
  
  INTEGER, parameter   :: Lenwrk = 32*Node
  INTEGER              :: rkmethod, outch
  REAL                 :: hstart, thres(1:Node), rkwork(1:Lenwrk), mchpes, dwarf
  REAL                 :: zgot, ygot(1:Node), ypgot(1:Node), ymax(1:Node)
  CHARACTER            :: rktask
  LOGICAL              :: errass, mesage
 
  external             :: bfield
  external             :: SETUP, UT, ENVIRN
  
  BEGIN(pp00ab)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  ivol = lvol ! required to pass through to bfield;
  
#ifdef DEBUG
  FATAL(pp00ab, lvol.lt.1.or.lvol.gt.Mvol, invalid volume )
  FATAL(pp00ab, abs(sti(1)).gt.one, illegal radial coordinate )
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
    
    CALL( pp00ab, stzxyz, ( lvol, stz(1:3), RpZ(1:3) ) ) ! map to cylindrical;
    
    ppt(3:4)=(/ RpZ(1), RpZ(3) /) ! cylindrical coordinates;
    
    poincaredata(1:4,kk,jj) = (/ mod(ppt(1),pi2), ppt(2), ppt(3), ppt(4) /) ! save Poincare information to array;
    
    zend = zst + (pi2nfp/Nz)
        
    utflag = 0
    
    call SETUP(Node, zst, st(1:Node), zend, tol, thres(1:Node), rkmethod, rktask, errass, hstart, rkwork(1:Lenwrk), Lenwrk, mesage) 
    
    CALL( pp00ab, UT, (bfield, zend, zgot, ygot(1:Node), ypgot(1:Node), ymax(1:Node), rkwork(1:Lenwrk), utflag) ) ! integrate to next plane;
    
    zst = zend
    
    st(1:Node) = ygot(1:Node)
    
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
  
  RETURN(pp00ab)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine pp00ab

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
