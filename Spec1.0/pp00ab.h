!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Constructs \Poincare plot and `approximate' rotational-transform for single field line.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine pp00ab( lvol, sti, Nz, nPpts, poincaredata, fittedtransform, id02bjf ) 
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, one, pi2
  
  use numerical, only : small
  
  use fileunits, only : ounit
  
  use inputlist, only : Wpp00ab, Nvol, pknot, odetol
  
  use cputiming, only : Tpp00ab
  
  use allglobal, only : myid, ncpu, cpus, pi2nfp, ivol, Ltangent, Lfieldlinedirection, Mvol
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in)  :: lvol, Nz, nPpts
  INTEGER, intent(out) :: id02bjf
  REAL                 :: sti(1:2), poincaredata(1:4,0:Nz-1,1:nPpts), fittedtransform(1:2), dzeta
  
  INTEGER              :: jj, kk
  REAL                 :: ppt(1:4)
  
  INTEGER              :: Node
  REAL                 :: zst, zend, st(1:6), realwork(1:20*6), tol, stz(1:3), RpZ(1:3), jacobian, dR(0:3), dZ(0:3), leastfit(1:5), guv(1:3,1:3)
  CHARACTER            :: RA
  
  external             :: bf00aa!, bf00aa_out, bf00aa_end
  external             :: D02BJX, D02BJW
  
  BEGIN(pp00ab)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  ivol = lvol ; Ltangent = 0 ! required to pass through to bf00aa;
  
#ifdef DEBUG
  FATALMESS(pp00ab, lvol.lt.1.or.lvol.gt.Mvol, invalid volume )
  FATALMESS(pp00ab, abs(sti(1)).gt.one, illegal radial coordinate )
#endif
  
  dzeta = pi2nfp * pknot

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \end{enumerate} \subsubsection{relevant input variables} \begin{enumerate}
  
!latex \item The resolution of \Poincare plot is controlled by 
!latex \begin{itemize} 
!latex \item \verb+nPpts+ iterations per trajectory;
!latex \item \verb+odetol+ o.d.e. integration tolerance;
!latex \end{itemize}
!latex The magnetic field is given by \verb+bf00aa+.
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  Node = 6 ; RA = 'D' ; tol = odetol
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  st(1:6) = (/ sti(1), sti(2), one, zero, zero, one /)!, zero, zero /) 
  
  ppt(1:2) = (/ sti(2), sti(1) /) ! seems redundant; just used for packing into poincaredata it seems;
  
  leastfit(1:5) = (/ zero, zero, zero, ppt(1), one /) ! initialize summation for least squares fit;
  
  id02bjf = 1 ; poincaredata(1:4,0:Nz-1,1:nPpts) = zero ; fittedtransform(1:2) = zero ! provide dummy defaults;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  do jj = 1, nPpts ! loop over iterations;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   zst = zero ! starting Poincare plane;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   do kk = 0, Nz-1 ! loop over toroidal Poincare cross sections;
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
    stz(1:3) = (/ ppt(2), mod(ppt(1),pi2), zst /) ! toroidal coordinates;
    
    if( abs(stz(1)).gt.one ) then
    !if( Wpp00ab ) then
      write(ounit,1002) myid, lvol, stz(1:2)
    !endif
     id02bjf = 1
     exit ! ! exit do kk loop; 22 Apr 13;
    endif

1002 format("pp00ab : ", 10x ," : myid=",i3," ; lvol=",i3," ; "3x" : (s,t)=("f21.17" ,"f21.17" ) ;         "3x" ; outside domain ;")
    
    CALL(pp00ab, co00aa,( lvol, stz(1:3), RpZ(1:3), dR(0:3), dZ(0:3), jacobian, guv(1:3,1:3) )) ! map to cylindrical;
    
    ppt(3:4)=(/ RpZ(1), RpZ(3) /) ! cylindrical coordinates;
    
    poincaredata(1:4,kk,jj) = (/ mod(ppt(1),pi2), ppt(2), ppt(3), ppt(4) /) ! save Poincare information to array;
    
    zend = zst + (pi2nfp/Nz) * Lfieldlinedirection
    
    id02bjf = 1
    CALL(pp00ab, D02BJF,( zst, zend, Node, st(1:Node), bf00aa, tol, RA, D02BJX, D02BJW, realwork, id02bjf )) ! integrate to next toroidal subdivision;
    
    cput = GETTIME
    select case( id02bjf )                                                !         1         2         3         4         5         6
    case( 0 ) ; ! give screen output if error is encountered;             !123456789012345678901234567890123456789012345678901234567890123
    case( 1 ) ; write(ounit,2001) cput-cpus, myid, lvol, jj, kk, id02bjf, "input error                                                    "
    case( 2 ) ; write(ounit,2001) cput-cpus, myid, lvol, jj, kk, id02bjf, "error integrating field                                        "
    case( 3 ) ; write(ounit,2001) cput-cpus, myid, lvol, jj, kk, id02bjf, "tol is too small to take initial step                          "
    case( 4 ) ; write(ounit,2001) cput-cpus, myid, lvol, jj, kk, id02bjf, "xsol not reset or xsol is behind x after initial call to output"
    case( 5 ) ; write(ounit,2001) cput-cpus, myid, lvol, jj, kk, id02bjf, "xsol not reset or xsol is behind last xsol                     "
    case( 6 ) ; write(ounit,2001) cput-cpus, myid, lvol, jj, kk, id02bjf, "termination function did not change sign                       "
    case( 7 ) ; write(ounit,2001) cput-cpus, myid, lvol, jj, kk, id02bjf, "serious error                                                  "
    case default
     FATALMESS(pp00ab,.true.,illegal value of ifail returned from D02BJF)
    end select
    
2001 format("pp00ab : ",f10.2," : myid=",i3," ; lvol=",i3," ; (jj,kk)=("i4" ,"i4" ); ifail="i2" ; "a63)
    
    if( id02bjf.ne.0 ) exit ! an integration error was encountered; exit do kk loop;
    
    ppt(1:2) = (/ st(2), st(1) /) ! again, seems redundant; I think this is only required for packing into poincaredata;
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
   enddo ! end of do kk;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   if( id02bjf.ne.0 ) exit ! an integration error was encountered; exit do jj loop;
   
! this is probably incorrect if Lfieldlinedirection = -1; perhaps better to use zend rather than jj*pi2nfp;
   
   leastfit(1:5) = leastfit(1:5) + (/ (jj*dzeta)**2, (jj*dzeta), (jj*dzeta)*ppt(1), ppt(1), one /) ! least squares fit summation;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
  enddo ! end of do jj = 1, nPpts;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \end{enumerate} \subsubsection{rotational-transform} \begin{enumerate}
  
!latex \item The approximate rotational transform is determined by field line integration.
!latex       This is constructed by fitting a least squares fit to the field line trajectory.
  
  if( id02bjf.eq.0 ) fittedtransform(1:2) = (/ sti(1), ( leastfit(5)*leastfit(3)-leastfit(2)*leastfit(4) ) / ( leastfit(5)*leastfit(1)-leastfit(2)*leastfit(2) ) /)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  Ltangent = 1
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(pp00ab)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine pp00ab

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
