!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Locates invariant, irrational surfaces (i.e. KAM surfaces) of the magnetic field, and writes input for pressure-jump Hamiltonian analysis code.

!latex \subsubsection{construction of invariant surface}

!latex \item The method for constructing an irrational surface is described in [Hudson, 2004] \cite{Hudson_04}.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

module curveflow
  INTEGER              :: jj, kk, Nk 
  REAL   , allocatable :: sap(:,:), tap(:,:) ! full invariant surface on regular (alpha,zeta) grid;
end module curveflow

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine invarianterr( lvol, Nc, Nxdof, xdof, err, derr, ljc, iuser, user ) ! NAG is no longer used; free to change argument list;
  
  use constants, only : zero, one, pi2

  use fileunits, only : ounit

  use inputlist, only : Wir00aa, odetol, irrMpol, irrNtor

  use cputiming, only : Tir00aa

  use allglobal, only : myid, cpus, irrsurf, pi2nfp, ivol, Ltangent
  
  use curveflow
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS

  INTEGER   , intent(in)  :: lvol
  
  INTEGER                 :: Nc, Nxdof, ljc, iuser(1)
  REAL                    :: xdof(0:Nxdof-1), err(0:Nc-1), derr(0:ljc-1,0:Nxdof-1), user(*)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  INTEGER                 :: mm, Node, ifail, Npts
  REAL                    :: zst, zend, tol, Tst(6), realwork(1:20*6)
  REAL      , allocatable :: ts(:), si(:), ti(:), tsb(:), sb(:), tb(:), st(:), tt(:), dst(:,:,:)
  CHARACTER               :: RA

  EXTERNAL                :: bf00aa, curveflowout, bf00aa_end

  BEGIN(ir00aa)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  ivol = lvol ; Ltangent = 1 ! required to pass through to bf00aa; and required for sslower; 01 Dec 12;

  FATALMESS(ir00aa,(Nxdof-1)/2 .ne. irrMpol, Fourier resolution problem) ! 01 Dec 12;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  Npts = irrsurf%Npts ! shorthand;

  Nk = max( 2*irrNtor, 1 ) ! discrete resolution in toroidal direction determined by Fourier resolution;
 
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
 
! straight-field line angle;
  ALLOCATE(ts,(0:Npts-1))
  ts(0:Npts-1) = (/ ( jj, jj = 0, Npts-1 ) /) * pi2/Npts ! straight-field line angle;
  
! rigid-shifted straight-field line angle;
  ALLOCATE(tsb,(0:Npts-1))
  tsb = ts + pi2nfp * irrsurf%iota ! rigid-shifted straight-field line angle;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!  mapped s, t;
  ALLOCATE(st,(0:Npts-1))
  ALLOCATE(tt,(0:Npts-1))
  
! tangent map;
  ALLOCATE(dst,(0:Npts-1,2,2))
  
! rigid shift s, t;
  ALLOCATE(sb,(0:Npts-1)) 
  ALLOCATE(tb,(0:Npts-1))
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \item A trial curve, represented by a discrete set of points $(s_i,\theta_i)$,
!latex \begin{equation}\begin{array}{ccccc}      s_i&=&                      & &\sum_{m=0}^{M}      s_m \cos(m \alpha_i ), \\
!latex                                      \theta_i&=&  \alpha              &+&\sum_{m=1}^{M} \theta_m \sin(m \alpha_i ),
!latex \end{array}\end{equation} 
!latex where $\alpha$ is a straight field line angle, and $\alpha_i$ are regularly spaced, and $M\equiv$\verb!irrMpol!,

  irrsurf%si(0:Npts-1) = xdof(0)
  irrsurf%ti(0:Npts-1) = ts(0:Npts-1) ! initial (s,t) ;

  do mm = 1, irrMpol
   irrsurf%si(0:Npts-1) = irrsurf%si(0:Npts-1) + xdof(        mm) * cos(mm*ts(0:Npts-1))
   irrsurf%ti(0:Npts-1) = irrsurf%ti(0:Npts-1) + xdof(irrMpol+mm) * sin(mm*ts(0:Npts-1))
  enddo

  irrsurf%si(Npts) = irrsurf%si(0)       ! include periodicity; required for interpolation below;
  irrsurf%ti(Npts) = irrsurf%ti(0) + pi2
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex is mapped under the field line flow $(s,\theta)\mapsto(\tilde s,\tilde \theta)$, 
!latex with the field (and tangent map) given by \verb+bf00ac+,
  
  do jj = 0, Npts-1 ! j is a counter provided to curveflowout below;

   kk = 0 ! k is a counter provided to curveflowout below;

   Tst(1:6) = (/ irrsurf%si(jj), irrsurf%ti(jj), one, zero, zero, one /)

   zst=zero ; zend=pi2nfp ; Node=6 ; tol=odetol ; RA='D' ; realwork=zero

   ifail=0
   call D02BJF( zst, zend, Node, Tst, bf00aa, tol, RA, curveflowout, bf00aa_end, realwork, ifail ) ! field-line Hamiltonian;
  !call D02BJF( zst, zend, Node, Tst, ph00aa, tol, RA, curveflowout, bfend, realwork, ifail ) ! pressure-jump Hamiltonian;

   st(jj)=Tst(1) ; tt(jj)=Tst(2) ; dst(jj,1,1)=Tst(3) ; dst(jj,1,2)=Tst(4) ; dst(jj,2,1)=Tst(5) ; dst(jj,2,2)=Tst(6)

  enddo
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
 
!latex and compared to the rigidly shifted curve, 
!latex \begin{equation}\begin{array}{ccccc} \bar      s&=&             & &\sum      s_m \cos(m \bar \alpha ), \\
!latex                                      \bar \theta&=& \bar \alpha &+&\sum \theta_m \sin(m \bar \alpha ),
!latex \end{array}\end{equation}
!latex where $\bar \alpha = \alpha+2\pi\iotabar$.
  
  sb(0:Npts-1) = xdof(0)
  tb(0:Npts-1) = tsb(0:Npts-1) ! rigid shift (s,t) ;

  do mm = 1,irrMpol
   sb(0:Npts-1) = sb(0:Npts-1) + xdof(        mm) * cos(mm*tsb(0:Npts-1))
   tb(0:Npts-1) = tb(0:Npts-1) + xdof(irrMpol+mm) * sin(mm*tsb(0:Npts-1))
  enddo
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \item The `function' we seek a zero of is the error vector, given ${\bf e}=(\bar s - \tilde s, \bar \theta - \tilde \theta)^T$.
  
  err(   0:  Npts-1) = sb(0:Npts-1)-st(0:Npts-1) ! error in s;
  err(Npts:2*Npts-1) = tb(0:Npts-1)-tt(0:Npts-1) ! error in t;
  
  irrsurf%err=sqrt( sum(err(0:2*Npts-1)**2)/(2*Npts) )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
 
!latex The correction to the harmonics, $s_m$ and $\theta_m$, that describe the curve are obtained using a Newton method, for which the derivative matrix is required:
!latex \begin{equation} \begin{array}{ccccc}

!latex \partial(\bar      s - \tilde      s)/\partial      s_m &=
!latex & \cos(m \bar \alpha) &-& \cos(m \alpha) \partial \tilde      s/\partial      s, \\ 
!latex \partial(\bar      s - \tilde      s)/\partial \theta_m &=
!latex &                     &-& \sin(m \alpha) \partial \tilde      s/\partial \theta, \\ 
!latex \partial(\bar \theta - \tilde \theta)/\partial      s_m &=
!latex &                     &-& \cos(m \alpha) \partial \tilde \theta/\partial      s, \\ 
!latex \partial(\bar \theta - \tilde \theta)/\partial \theta_m &=
!latex & \sin(m \bar \alpha) &-& \sin(m \alpha) \partial \tilde \theta/\partial \theta,
!latex \end{array} \end{equation}
!latex where the $\partial \tilde x / \partial y$ are elements of the tangent map.

  derr(0:2*Npts-1,0:2*irrMpol)=zero ! error derivative matrix;
  
   derr(   0:  Npts-1,0)         = cos(0*tsb) - dst(0:Npts-1,1,1)!* cos(0*ts)
   derr(Npts:2*Npts-1,0)         =            - dst(0:Npts-1,2,1)!* cos(0*ts)
  
  do mm = 1,irrMpol
   derr(   0:  Npts-1,mm        ) = cos(mm*tsb) - dst(0:Npts-1,1,1) * cos(mm*ts)
   derr(   0:  Npts-1,mm+irrMpol) =             - dst(0:Npts-1,1,2) * sin(mm*ts)
   derr(Npts:2*Npts-1,mm        ) =             - dst(0:Npts-1,2,1) * cos(mm*ts)
   derr(Npts:2*Npts-1,mm+irrMpol) = sin(mm*tsb) - dst(0:Npts-1,2,2) * sin(mm*ts)
  enddo
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  DEALLOCATE(ts)
  DEALLOCATE(tsb)
  
  DEALLOCATE(st)
  DEALLOCATE(tt)
  
  DEALLOCATE(dst)
  
  DEALLOCATE(sb)
  DEALLOCATE(tb)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  iuser(1) = iuser(1) + 1 ! number of function evaluations;

  return
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine invarianterr

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine ir00aa( lvol )
  
  use constants, only : zero, half, one, two, pi2

  use fileunits, only : ounit, wunit

  use inputlist, only : Wir00aa, ext, Nfp, irrMpol, irrNtor, odetol, irrsvdcut, Mirrits, irrtol, p1, q1, p2, q2, tflux, absacc

  use cputiming

  use allglobal, only : myid, cpus, pi2nfp, irrsurf,  irrmn, irrim, irrin, Bzeta, ivol

  use curveflow

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS 

  INTEGER, intent(in)  :: lvol

  INTEGER              :: jk, mm, nn, Nc, Nxdof, ifail, its, Npts, imn, iuser(1), ljc!, Lw
  REAL                 :: zeta, ths, cs(2), arg, gtt, gtp, gpp, osqrtNptsNk, user(1)!, sumsq, cyl_err
  REAL                 :: dRsp(0:3), dZsp(0:3)
  REAL,  allocatable   :: xdof(:), err(:), derr(:, :), dm(:), wsvd(:)!, MWk(:)
  REAL,  allocatable   :: sjk(:), tjk(:), Rjk(:, :), Zjk(:, :), BR(:), Bp(:), BZ(:), Blow(:, :), irrBtmn(:), irrBpmn(:)
  REAL,  allocatable   :: E(:), trm(:), trn(:), Fwk(:) ! do not confuse these with the global arrays; 

  REAL                 :: stz(6), Bstz(6) ! really just position,  but bf00aa also returns tangent map so input array is length 6;
  REAL                 :: RpZ(3), dR(0:3), dZ(0:3), jacobian, guv(1:3,1:3)

  CHARACTER            :: isr, suff*19

  LOGICAL              :: Lsvd=.true. !, Lnag=.false. ! it seems that homespun svd works better than nag;
 !external             :: invarianterr ! only required if Lnag=.true.

  INTEGER              :: icount,irank
  INTEGER, allocatable :: rank(:)

  INTEGER              :: nfunctioncalls
  REAL                 :: thetalower, thetaupper, sslower, ssupper, differentialflux, absaccuracyreqd, estimatedintegral, estimatedtflux
  external             :: sslower, ssupper, differentialflux 

  CHARACTER            :: message*6

  BEGIN(ir00aa)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Mirrits.le.0 ) goto 9999

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  ivol = lvol ! required to pass through to sslower; 01 Dec 12;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  Nxdof = 2*irrMpol + 1 ! degrees of freedom in stellarator symmetric invariant curve;

  Npts = irrsurf%Npts ! irrsurf%Npts is defined in al00aa; 

  Nc = 2*Npts ; Nk = max(2*irrNtor,1) ! discrete resolution determined by Fourier resolution;
  osqrtNptsNk = one/sqrt(one*Npts*Nk)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! full surface s, t on regular (alpha,zeta) grid; including periodic curve;
  ALLOCATE(sap,(0:Npts,0:Nk-1)) 
  ALLOCATE(tap,(0:Npts,0:Nk-1))

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! independent variable, error vector and derivative matrix; provided by tagent map;
  ALLOCATE(xdof,(0:2*irrMpol))
  ALLOCATE(err,(0:2*Npts-1))
  ALLOCATE(derr,(0:2*Npts-1,0:2*irrMpol))

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( Lsvd ) then
   
! singular values and correction;
   ALLOCATE(wsvd,(1:Nxdof)) 
   ALLOCATE(dm,(0:2*irrMpol))
   
   iuser(1) = 0 ; ifail = 1
   
!latex \item Given an initial guess for the irrational invariant {\em curve}, 
   
   do its = 1, Mirrits ! until error is satisfactory;
    
    xdof(        0:  irrMpol) = irrsurf%sm(0:irrMpol)
    xdof(irrMpol+1:2*irrMpol) = irrsurf%tm(1:irrMpol)
    
!latex the routine \verb+invarianterr+ is called iteratively (maximum iterations given by \verb+Mirrits+)
    
    ljc = Nc
    CALL(invarianterr,( lvol, Nc, Nxdof, xdof, err, derr, ljc, iuser, user ))
    
!latex until the root mean square error is less than \verb+irrtol+.
    
    irrsurf%err = sqrt(sum(err(0:2*Npts-1)**2)/(2*Npts))
    
!   if( irrsurf%err.le.irrtol ) then ; ifail=0 ; exit ! always perform max. iterations allowed, as set by Mirrits;
!   endif
    
!latex Corrections are calculated by singular value decomposition inversion (see \verb+singvalues+) of the derivative matrix.
    
    call singvalues( Nc, Nxdof, derr( 0:2*Npts-1, 0:2*irrMpol ), err(0:2*Npts-1), dm(0:2*irrMpol), irrsvdcut, wsvd )
    
    irrsurf%sm(0:irrMpol) = irrsurf%sm(0:irrMpol) - dm(        0:  irrMpol) ! update s-Fourier harmonics of curve;
    irrsurf%tm(1:irrMpol) = irrsurf%tm(1:irrMpol) - dm(irrMpol+1:2*irrMpol) ! update t-Fourier harmonics of curve;
    
    cput = GETTIME
    write(ounit,2001) cput-cpus, myid, lvol, irrsurf%iota, irrMpol, irrsurf%sm(0), irrsurf%err, iuser(1), ifail, cput-cpuo
    cpuo=cput

   enddo ! end of do its = 1, Mirrits; until error is satisfactory;
   
   irrsurf%its = iuser(1)
   
   DEALLOCATE(wsvd)
   DEALLOCATE(dm)
   
  endif

2001 format("ir00aa : ",f10.2," : myid=",i3," ; lvol=",i3," ; iota="f15.12" ; M="i4" ; sm(0)="f9.6" ; ":"rms="es12.5" ; its="i4" ; ":"ifail="i2" ; time=",f10.2,"s ;")

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item (There is an option to solve sum-of-squares using NAG routine \verb+E04GYF+.
!latex This routine does not perform as well as the SVD routine above. 
!latex I think this is probably due to integration errors, which may confuse the NAG routine.)
!  if(Lnag) then 
!
!   xdof(                 0:  irrMpol) = irrsurf%sm(0:irrMpol)
!   xdof(irrMpol+1:2*irrMpol) = irrsurf%tm(1:irrMpol)
!   
!   Lw = 8*Nxdof + 2*Nxdof*Nxdof + 2*Nc*Nxdof + 3*Nc
!   
!   ALL!OCATE(MWk,(Lw))
!   MWk=zero
!   iuser=0 ; user=zero
!   
!   ifail=1 ; call E04GYF( Nc, Nxdof, invarianterr, xdof, sumsq, MWk, Lw, iuser, user, ifail ) ; irrsurf%err=sqrt(sumsq/(2*Npts))
!   
!   DEALL!OCATE(MWk)
!   
!   irrsurf%sm(0:irrMpol)=xdof(          0 :  irrMpol)
!   irrsurf%tm(1:irrMpol)=xdof(irrMpol+1:2*irrMpol)
!
!   irrsurf%its=iuser(1)
!
!   cpu(1)=GETTIME 
!   if(WF(5)) write(ounit,2001)cpu(1)-cpustart,Hamiltonian,iota,irrMpol,irrsurf%sm(0),irrsurf%err,iuser(1),ifail ; cpuo=cput
!
!  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  DEALLOCATE(xdof)
  DEALLOCATE(err)
  DEALLOCATE(derr)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item NOTE: If the root mean square error tolerance is greater than \verb+irrtol+, a Fourier representation of the invariant {\em surface} will {\bf not} be constructed.

  cput = GETTIME
  write(ounit,'("ir00aa : ",f10.2," : myid=",i3," ; lvol=",i3," ; irrational rms error ="es12.5" ; irrtol ="es9.2" ;")') cput-cpus, myid, lvol, irrsurf%err, irrtol ; cpuo=cput

  if( irrsurf%err.gt.irrtol ) then
   DEALLOCATE(sap) ; DEALLOCATE(tap)
   goto 9999 ! need to check below if file output will be corrupted by this abrupt return;
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item The invariant {\em surface} is obtained by allowing the invariant {\em curve} to `flow' along the field. 
!latex Note: the surface thus constructed is initially obtained on a regular $(\alpha,\zeta)$ grid, where $\alpha$ labels field lines.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  sap(Npts,0:Nk-1) = sap(0,0:Nk-1)       ! full surface on complete domain; regular alpha grid; destroyed by Fourier decomposition;
  tap(Npts,0:Nk-1) = tap(0,0:Nk-1) + pi2 ! full surface on complete domain; regular alpha grid; destroyed by Fourier decomposition;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! full surface s, t on regular (theta_s,zeta) grid;
  ALLOCATE(sjk,(Npts*Nk))
  ALLOCATE(tjk,(Npts*Nk))

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex To obtain a regular grid in $(\theta_s,\zeta)$, where $\theta_s$ is the straight-field line angle consistent with $\zeta$, it is required to interpolate.
!latex It is (probably) more accurate to interpolate to get the background coordinates $(s,\theta,\zeta)$, which probably vary slightly, 
!latex and then perform the (exact) transformation to cylindrical coordinates, $(R,\phi,Z)$.
  
  do kk = 0, Nk-1 ; zeta = kk * pi2nfp/Nk

   NN = Npts ; ifail = 1 ; call slowft(sap(0:Npts-1,kk),NN,ifail)

   sap(0:Npts-1,kk) = sap(0:Npts-1,kk) / sqrt(one*Npts) ; sap(1:Npts-1,kk) = sap(1:Npts-1,kk)*two

   tap(0:Npts,kk) = tap(0:Npts,kk) - ( (/(jj,jj=0,Npts)/)*pi2/Npts + irrsurf%iota*zeta ) ! lambda;

   NN = Npts ; ifail = 1 ; call slowft(tap(0:Npts-1,kk),NN,ifail)

   tap(0:Npts-1,kk) = tap(0:Npts-1,kk) / sqrt(one*Npts) ; tap(1:Npts-1,kk) = tap(1:Npts-1,kk)*two

   do jj = 0, Npts-1 ; jk = 1 + jj + Npts*kk ; ths = jj * pi2/Npts ! ; alpha = ths - iota*zeta;

    sjk(jk) = sap(0,kk) ; tjk(jk) = ths + tap(0,kk)

    do mm = 1, irrMpol ; cs = (/ cos( mm*(ths-irrsurf%iota*zeta) ), sin( mm*(ths-irrsurf%iota*zeta) ) /) ! grid is irregular;
     sjk(jk) = sjk(jk) + sap(mm,kk)*cs(1) - sap(Npts-mm,kk)*cs(2)
     tjk(jk) = tjk(jk) + tap(mm,kk)*cs(1) - tap(Npts-mm,kk)*cs(2)
    enddo

   enddo
  enddo
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! full surface R, Z on regular (theta_s,zeta) grid;
  ALLOCATE(Rjk,(2,Npts*Nk)) 
  ALLOCATE(Zjk,(2,Npts*Nk))

! covariant cylindrical B^R, B^phi, B^Z on regular (theta_s,zeta) grid;
  ALLOCATE(BR,(Npts*Nk)) 
  ALLOCATE(Bp,(Npts*Nk))
  ALLOCATE(BZ,(Npts*Nk))
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  Rjk(1,:) = zero ; Zjk(1,:) = zero
  
  do jj = 0, Npts-1
   do kk = 0, Nk-1 ; jk = 1 + jj + Npts*kk ; zeta = kk * pi2nfp/Nk

    stz(1:3)=(/ sjk(jk), tjk(jk), zeta /)

#ifdef DEBUG
   FATALMESS(ir00aa, abs(stz(1)).gt.one, illegal value)
#endif

    WCALL(ir00aa,co00aa,( lvol, stz(1:3), RpZ(1:3), dR(0:3), dZ(0:3), jacobian, guv(1:3,1:3) ))
    
    Rjk(1,jk) = RpZ(1) ; Zjk(1,jk) = RpZ(3)

! contravariant field in background coordinates; NOTE: missing Jacobian factors;
    stz(3:6) = (/ one, zero, zero, one /)
    CALL(bf00aa,( zeta, stz(1:6), Bstz(1:6) ))

!latex \item The contravariant field is transformed to cylindrical coordinates via
!latex \be 
!latex \left ( \begin{array}{c} B^R \\ B^\phi \\ B^Z \end{array} \right)  = 
!latex \left ( \begin{array}{ccc}   R_s  &   R_\theta  &   R_\zeta  \\ 
!latex                               0   &      0      &     -1     \\ 
!latex                              Z_s  &   Z_\theta  &   Z_\zeta  \end{array} \right)
!latex \left ( \begin{array}{c} B^s \\ B^\theta \\ B^\zeta \end{array} \right).
!latex \ee
!latex Recall that $B_\phi = R^2 B^\phi$.
    
    Bstz(1:3) = (/ Bstz(1), Bstz(2), one /) * (Bzeta/jacobian) ! Bzeta=B^\zeta is returned through global memory;
    
    BR(jk) = dR(1)*Bstz(1) + dR(2)*Bstz(2) + dR(3)*Bstz(3)
    Bp(jk) =                                     - Bstz(3) * Rjk(1,jk)**2 ! convert contravariant to covariant;
    BZ(jk) = dZ(1)*Bstz(1) + dZ(2)*Bstz(2) + dZ(3)*Bstz(3)
    
   enddo
  enddo

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  DEALLOCATE(sjk)
  DEALLOCATE(tjk)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! these are local trigonometric arrays; do not confuse with global arrays;
  ALLOCATE(E,(Npts*Nk))     
  ALLOCATE(trm,(2*Npts))
  ALLOCATE(trn,(2*Nk))
  ALLOCATE(Fwk,(2*Npts*Nk))

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \item The surface is Fourier decomposed in the straight field line angle, 
!latex with harmonics idenified by \verb+irrim(1:irrmn)+ and \verb+irrin(1:irrmn)+, which are constructed in \verb+readin+. 
!latex \item The surface geometry harmonics are stored in the global arrays \verb+irrsurf%Rbc(:)+, \verb+irrsurf%Zbs(:)+.

  Rjk(2,:) = Rjk(1,:) ; Zjk(2,:) = Zjk(1,:) ! save original data for Fourier reconstruction check;

  E=zero ; isr='I' ; ifail=0 ; call C06FUF( Npts, Nk, Rjk(2,1:Npts*Nk), E, isr, trm, trn, Fwk, ifail )

  Rjk(2,:) = Rjk(2,:)*osqrtNptsNk ; Rjk(2,2:) = Rjk(2,2:)*two

  E=zero ; isr='S' ; ifail=0 ; call C06FUF( Npts, Nk, Zjk(2,1:Npts*Nk), E, isr, trm, trn, Fwk, ifail )

  E = E*osqrtNptsNk ; E(2:)= E(2:)*two ! Fourier decomposition;

  do imn = 1, irrmn ; mm = irrim(imn) ; nn = irrin(imn)/Nfp ; jj = mm ; kk = modulo(Nk-nn,Nk)

   irrsurf%Rbc(imn) = Rjk(2,1+jj+kk*Npts) ; irrsurf%Zbs(imn) = - E(1+jj+kk*Npts)

  enddo

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  DEALLOCATE(Rjk)
  DEALLOCATE(Zjk)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  DEALLOCATE(sap)
  DEALLOCATE(tap)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! tangential contravariant B^s, and covariant components B_{\theta_s}, B_\zeta, where s is irrational surface label;
  ALLOCATE(Blow,(1:5,Npts*Nk)) 
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \item The vector transformation to straight-field-line coordinates should exploit fast Fourier transforms.

  do jj = 0, Npts-1 ; ths = jj * pi2/Npts ! construct straight-field-line coordinate derivatives, perform vector transformation;
   do kk = 0, Nk-1 ; zeta =kk * pi2nfp/Nk ; jk = 1 + jj + Npts*kk
  
    dRsp(0:3) = zero ; dZsp(0:3) = zero
  
    do imn = 1, irrmn ; arg = irrim(imn)*ths - irrin(imn)*zeta ; cs=(/ cos(arg), sin(arg) /)
     
     dRsp(0) = dRsp(0) + irrsurf%Rbc(imn)*cs(1)
    !dRsp(1) = ?? no radial coordinate is constructed for discrete irrational surfaces;
     dRsp(2) = dRsp(2) + irrsurf%Rbc(imn)*cs(2)*(-irrim(imn) )
     dRsp(3) = dRsp(3) + irrsurf%Rbc(imn)*cs(2)*( irrin(imn) )
    !dZsp(0) = dZsp(0) + irrsurf%Zbs(imn)*cs(2) ; not required;            
    !dZsp(1) = ?? no radial coordinate is constructed for discrete irrational surfaces;
     dZsp(2) = dZsp(2) + irrsurf%Zbs(imn)*cs(1)*( irrim(imn) )     
     dZsp(3) = dZsp(3) + irrsurf%Zbs(imn)*cs(1)*(-irrin(imn) )
     
    enddo

!latex \item The covariant field is transformed 
!latex \be 
!latex \left ( \begin{array}{c} B_s \\ B_\theta \\ B_\zeta \end{array} \right)  = 
!latex \left ( \begin{array}{ccc}   R_s          &   0  &   Z_s           \\ 
!latex                              R_{\theta_s} &   0  &   Z_{\theta_s}  \\ 
!latex                              R_\zeta      &  -1  &   Z_\zeta       \end{array} \right)
!latex \left ( \begin{array}{c} B_R \\ B_\phi \\ B_Z \end{array} \right).
!latex \ee

    Blow(2,jk) = dRsp(2) * BR(jk)          + dZsp(2) * BZ(jk)
    Blow(3,jk) = dRsp(3) * BR(jk) - Bp(jk) + dZsp(3) * BZ(jk)
    
! the following is for checking the field normal to the invariant surface is zero, and the transform is correct, . . .;
    
!latex \item The normal field is 
!latex \be {\bf B} \cdot {\bf e}_{\theta_s} \times {\bf e}_{\zeta} = 
!latex Z_{\theta_s} R B_R + ( Z_{\theta_s} R_\zeta - R_{\theta_s}Z_\zeta ) B_\phi R^{-1} - R_{\theta_s} R B_z.
!latex \ee
    
! B^s; omitting jacobian factor; debugging;
    Blow(1,jk) = dZsp(2)*dRsp(0)*BR(jk) + ( dZsp(2)*dRsp(3)-dRsp(2)*dZsp(3) ) * Bp(jk) / dRsp(0) - dRsp(2)*dRsp(0)*BZ(jk) 
    
!latex \item The transform is given 
!latex \be \iotabar = \frac{B^{\theta_s}}{B^\zeta} = \frac{ B_{\theta_s} g_{\zeta\zeta} - B_\zeta g_{\theta_s\zeta}}
!latex                     { -B_{\theta_s} g_{\theta_s\zeta} + B_\zeta g_{\zeta\zeta}}.
!latex \ee
    
    gtt = dRsp(2)*dRsp(2) + dZsp(2)*dZsp(2)
    gtp = dRsp(2)*dRsp(3) + dZsp(2)*dZsp(3)
    gpp = dRsp(3)*dRsp(3) + dZsp(3)*dZsp(3) + dRsp(0)**2
    
    Blow(4,jk) =  Blow(2,jk)*gpp - Blow(3,jk)*gtp ! B^{\theta_s} (well, not really -- missing some normalization factors);
    Blow(5,jk) = -Blow(2,jk)*gtp + Blow(3,jk)*gtt ! B^{\zeta}    (well, not really -- missing same normalization factors);
    
   enddo
  enddo

  
  if( Wir00aa ) then
   cput = GETTIME
   write(ounit,2005) cput-cpus, myid, lvol, sqrt((/sum(Blow(1,:)**2), sum((irrsurf%iota-Blow(4,:)/Blow(5,:))**2)/)/(Npts*Nk)), cput-cpuo
   cpuo=cput
  endif

2005 format("ir00aa : ",f10.2," :":" myid=",i3," ; lvol=",i3," ; covariant transformation ; normal field error ="es10.2" ; transform error ="es10.2" ; time=",f10.2,"s ;")
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  DEALLOCATE(BR)
  DEALLOCATE(Bp)
  DEALLOCATE(BZ)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! Fourier harmonics of irrational surfaces : tangential covariant field;

  ALLOCATE(irrBtmn,(1:irrmn))
  ALLOCATE(irrBpmn,(1:irrmn))

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  Blow(4,:)=Blow(2,:) ; Blow(5,:)=Blow(3,:) ! save for Fourier reconstruction check;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  E = zero ; isr = 'I' ; ifail = 0 ; call C06FUF( Npts, Nk, Blow(4,1:Npts*Nk), E, isr, trm, trn, Fwk, ifail )
  Blow(4,:) = Blow(4,:) * osqrtNptsNk ; Blow(4,2:) = Blow(4,2:)*two
  E = zero ; isr = 'S' ; ifail = 0 ; call C06FUF( Npts, Nk, Blow(5,1:Npts*Nk), E, isr, trm, trn, Fwk, ifail )
  Blow(5,:) = Blow(5,:) * osqrtNptsNk ; Blow(5,2:) = Blow(5,2:)*two
  do imn = 1, irrmn ; mm=irrim(imn) ; nn=irrin(imn) ; jj = mm ; kk = modulo(Nk-nn,Nk) ; irrBtmn(imn)= Blow(4,1+jj+kk*Npts) ; irrBpmn(imn)= Blow(5,1+jj+kk*Npts)
  enddo

  cput = GETTIME

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  DEALLOCATE(Blow)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! these are local trigonometric arrays; do not confuse with global arrays;
  DEALLOCATE(E)
  DEALLOCATE(trm)
  DEALLOCATE(trn)
  DEALLOCATE(Fwk)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item The surface potential $f$, where $B_{\theta_s}=\partial_{\theta_s} f$ and $B_{\zeta}=\partial_{\zeta} f$, is calculated and written to \verb!ext.irr.**.pj!. 

!l!tex \item The surface potential $f$, where $B_{\theta_s}=\partial_{\theta_s} f$ and $B_{\zeta}=\partial_{\zeta} f$, \verb+irrsurf%fmn(:)+.
!l!tex The secular terms in the surface potential, ie. $I$ and $G$, are by the `constant' harmonics, and saved in \verb+irrsurf%Itor+, \verb+irrsurf%Gpol+.

  irrsurf%Itor=zero ; irrsurf%Gpol=zero ; irrsurf%fmn(:)=zero

  do imn=1,irrmn ! integrate to get surface potential from covariant components of field;
   if( irrim(imn).eq.0 .and. irrin(imn).eq.0 ) irrsurf%Itor     = irrBtmn(imn)
   if( irrim(imn).eq.0 .and. irrin(imn).eq.0 ) irrsurf%Gpol     = irrBpmn(imn)
   if( irrim(imn).eq.0 .and. irrin(imn).ne.0 ) irrsurf%fmn(imn) =                             - irrBpmn(imn) / (irrin(imn)/Nfp)
   if( irrim(imn).ne.0 .and. irrin(imn).eq.0 ) irrsurf%fmn(imn) =   irrBtmn(imn) / irrim(imn)   
   if( irrim(imn).ne.0 .and. irrin(imn).ne.0 ) irrsurf%fmn(imn) = ( irrBtmn(imn) / irrim(imn) - irrBpmn(imn) / (irrin(imn)/Nfp) ) * half ! should be equal ?;
  enddo

  ALLOCATE(rank,(irrmn))
  ifail=0 ; call M01DAF(abs(irrsurf%fmn),1,irrmn,'D',rank,ifail) ! rank the surface potential harmonics;
  
  if( Wir00aa ) then 
   icount = 0
   do irank = 1, irrmn ! loop over all the ranks;
    do imn = 1, irrmn ! loop over all the harmonics;
     if( rank(imn).eq.irank ) then
      if( irrim(imn).ne.0 .and. irrin(imn).ne.0 ) then ; icount=icount+1 ! can only compare the m.ne.0 & n.ne.0 harmonics;
       write(ounit,2008) irrim(imn), irrin(imn), irrBtmn(imn) / irrim(imn), - irrBpmn(imn) / (irrin(imn)/Nfp), ( irrBtmn(imn) / irrim(imn) + irrBpmn(imn) / (irrin(imn)/Nfp) )
      endif
     endif
    enddo
    if( icount.ge.5 ) exit ! only show the top harmonics (to reduce screen clutter);
   enddo
  endif

  DEALLOCATE(rank)
  
2008 format("ir00aa : "10x" : surface potential mode (m,n)=(",i3,","i4" ) ; ( Btmn/m, -Bpmn/n )=("es11.3","es11.3" ) ;":" err="es9.1" ; ") ! should compare . . . .

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  DEALLOCATE(irrBtmn)
  DEALLOCATE(irrBpmn)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \subsubsection{input file for pressure-jump Hamiltonian analysis code}

!latex \item The information required for the pressure-jump Hamiltonian analysis is written to \verb+ext.p1:q1:p2:q2:M:N.pjh+,
!latex where the transform of the irrational surface is given by $\iotabar = (p_1 + \gamma p_2)/(q_1 + \gamma q_2)$, where $\gamma=( 1+ \sqrt 5 ) / 2$ is the golden mean,
!latex and $M$ and $N$ describe the poloidal and toroidal Fourier resolution.

!latex \item The format of this file is as follows:
!latex \begin{verbatim}
!latex Mpol Ntor Nfp             
!latex Itor Gpol p1 q1 p2 q2 iota
!latex m n fmn Rmn Zmn
!latex . . ... ... ...
!latex \end{verbatim}

!latex \item The tangential magnetic field on the surface is given by a surface potential, $f=I \t + G \z + \sum f_{m,n}\sin(m\t-n\z)$,
!latex via $B_\t = \partial_\t f$ and $B_\z = \partial_\z f$.
!latex \item The surface geometry is given by $R=\sum R_{m,n}\cos(m\t-n\z)$ and $Z=\sum Z_{m,n}\sin(m\t-n\z)$.

  write(suff,'(i2.2":"i2.2":"i2.2":"i2.2":",i3,.3":",i3,.3)')p1(lvol),q1(lvol),p2(lvol),q2(lvol),irrMpol,irrNtor
  
  if( Wir00aa ) then
   cput = GETTIME
   write(ounit,'("ir00aa : ",f10.2," : myid=",i3," ; lvol=",i3," ; writing ext."a19".pjh ; time=",f10.2,"s ;")') cput-cpus, myid, lvol, suff, cput-cpuo ; cpuo=cput
  endif
  
  open(wunit+myid,file=trim(ext)//"."//suff//".pjh",status="unknown")
  write(wunit+myid,'(3i6)')irrMpol,irrNtor,Nfp
  write(wunit+myid,'(2es23.15,4i6,es23.15)')irrsurf%Itor,irrsurf%Gpol,p1(lvol),q1(lvol),p2(lvol),q2(lvol),irrsurf%iota!=(p1*goldenmean+p2)/(q1*goldenmean+q2)
  do imn = 1,irrmn ; write(wunit+myid,'(2i6,3es23.15)')irrim(imn),irrin(imn),irrsurf%fmn(imn),irrsurf%Rbc(imn),irrsurf%Zbs(imn)
  enddo
  close(wunit+myid)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{calculation of enclosed flux}

!latex \item The toroidal flux between the irrational surface and the inner interface is calculated as a double integral using the NAG routine \verb+E01BEF+.

!latex \item The toroidal flux is given by $\psi\equiv \int_{\cal S} {\bf B}\cdot d{\bf s}$,
!latex where ${\cal S}$ is the surface on the $\z=0$ plane bounded by the $l-1$ interface and the irrational surface, $\bar s \equiv \bar s(\t)$,
!latex and ${\bf B}=\nabla \times (A_\t \nabla \t + A_\z \nabla \z)$.
!latex \be \psi = \int_{0}^{2\pi}\!\!\!\!\!\!d\t \int_{s_{l-1}}^{\bar s} \!\!\!\!\!\!ds \; \partial_s A_\t.
!latex \ee

!latex \item (The toroidal flux between the outer interface and the irrational interface should also be constructed for debugging and accuracy confirmation.)
  
!latex \item The accuracy of the calculation is controlled by the input parameter \verb+absacc+.

  ifail = 0 ; call E01BEF( Npts + 1, irrsurf%ti(0:Npts), irrsurf%si(0:Npts), irrsurf%di(0:Npts), ifail ) ! interpolation of upper integration boundary;

  thetalower = zero ! lower limit of outer integral; required on entry;
  thetaupper = pi2  ! upper limit of outer integral; required on entry;

  absaccuracyreqd = absacc ! absolute accuracy required; given on exit;

  if( Wir00aa ) then ; cput = GETTIME ; write(ounit,1002) cput-cpus, myid, lvol, absaccuracyreqd, "calling"
  endif
  
  ifail=1 ; call D01DAF( thetalower, thetaupper, sslower, ssupper, differentialflux, absaccuracyreqd, estimatedintegral, nfunctioncalls, ifail )
  
 !if( ivol.eq.1 ) then ; irrsurf%tflux = tflux(lvol  ) + estimatedintegral/pi2 ! WHEN ARE YOU GOING TO FIX THIS FACTOR OF PI2 ??? ! maybe I have already; 01 Jul 14;
 !else                 ; irrsurf%tflux = tflux(lvol-1) + estimatedintegral/pi2 ! WHEN ARE YOU GOING TO FIX THIS FACTOR OF PI2 ???
 !endif

  if( ivol.eq.1 ) then ; irrsurf%tflux = tflux(lvol  ) + estimatedintegral ! 18 Jul 14;
  else                 ; irrsurf%tflux = tflux(lvol-1) + estimatedintegral ! 18 Jul 14;
  endif

  if( ifail.ne.0 ) then ; message="FAILED"
  else                  ; message="      "
  endif

  cput = GETTIME ; write(ounit,1002) cput-cpus, myid, lvol, absaccuracyreqd, "called ", estimatedintegral, irrsurf%tflux, nfunctioncalls, ifail, message

1002 format("ir00aa : ",f10.2," : myid=",i3," ; lvol=",i3," : accuracy required ="es8.1" ; "a7" D01DAF :":" estimate ="2es23.15" ; #calls="i6" ; ifail="i6" ; ",a6)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(ir00aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine ir00aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

REAL function sslower( teta )
  
  use inputlist, only : Nvol
  
  use allglobal, only : myid, cpus, ivol, interfacelabel
  
  LOCALS
  
  REAL, intent(in) :: teta

  FATALMESS(ir00aa,ivol.lt.1.or.ivol.gt.Nvol,invalid ivol)
  
  if( ivol.eq.1 ) then                   ; sslower = interfacelabel(ivol  ) ! radial coordinate, ss, integration lower limit;
  endif
  if( ivol.gt.1 .and. ivol.le.Nvol) then ; sslower = interfacelabel(ivol-1) ! radial coordinate, ss, integration lower limit;
  endif

  return
  
end function sslower

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

REAL function ssupper( teta )
  
  use constants, only : zero, one
  use inputlist, only : irrMpol
  use allglobal, only : myid, cpus, ivol, irrsurf
  
  LOCALS
  
  REAL, intent(in) :: teta
  
  INTEGER          :: Npts, Mpts, ifail!, imn
  REAL             :: PX(1:1), PF(1:1) !zeta, arblabel, ss, tt, arg, carg, sarg, stz(1:8), Bstz(1:8)  
  
  Npts = irrsurf%Npts ! shorthand;
  
  Mpts = 1 ; PX(1:1) = teta 
  
  ifail = 0 ; call E01BFF( Npts + 1, irrsurf%ti(0:Npts), irrsurf%si(0:Npts), irrsurf%di(0:Npts), Mpts, PX(1:1), PF(1:1), ifail )
  
  ssupper = PF(1) ! radial coordinate, ss, integration upper limit;
  
  return

end function ssupper

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

REAL function differentialflux( lss, teta ) ! sqrt g * B^\zeta ;

  use constants, only : zero
  use inputlist, only : Wir00aa, Lpoincare
  use allglobal, only : myid, cpus, Bzeta

  LOCALS

  REAL, intent(in) :: lss, teta

  REAL             :: arblabel, stz(1:8), Bstz(1:8)

  FATALMESS(ir00aa,Lpoincare.ne.0,needs attention)

  arblabel = zero ; stz(1:8)=(/ lss, teta, zero, zero, zero, zero, zero, zero /)
  
  CALL(bf00aa,( arblabel, stz(1:8), Bstz(1:8) ) )
  
  differentialflux = Bzeta

  return

end function differentialflux

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine curveflowout(iz,Tst) ! used for mapping invariant curve through to invariant surface;

  use allglobal, only : pi2nfp
  use curveflow

  LOCALS

  REAL, intent(inout) :: iz ! toroidal (integration) variable;
  REAL, intent(in)    :: Tst(6)

  sap(jj,kk) = Tst(1) ; tap(jj,kk) = Tst(2) ! fill in surface array;

  kk = kk+1 ; iz = kk * pi2nfp / Nk

  return

end subroutine curveflowout

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
