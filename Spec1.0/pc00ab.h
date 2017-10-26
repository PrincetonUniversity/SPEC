!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Returns the energy functional, and its derivatives with respect to geometry, in a format suitable for NAG routine E04DGF.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{Energy functional}

!latex \item The energy functional is
!latex       \be F \equiv \sum_{l=1}^{N} \int_{\cal V} \left( \frac{p}{\gamma-1} + \frac{B^2}{2} \right) dv,
!latex       \label{eq:energyfunctional}
!latex       \ee
!latex       where $N \equiv$\verb+Nvol+ is the number of interfaces.

!latex \item Assuming that the toroidal and poloidal fluxes, $\psi_t$ and $\psi_p$, the helicity, ${\cal K}$, the helicity multiplier, $\mu$,
!latex       and/or the interface rotational-transforms, $\iotabar$, are appropriately constrained,
!latex       the Beltrami fields in each volume depend only the geometry of the adjacent interfaces.
!latex       So, the energy functional is assumed to be a function of ``position'', i.e. $F = F(R_{l,j},Z_{l,j})$.

!latex \item Introducing a ficitious time, $t$, the position may be advanced according to 
!latex       \be \begin{array}{cccccccccccccccccccccccc}
!latex           \ds \frac{\partial R_j}{\partial t} & \equiv & \ds - \frac{\partial }{\partial R_j} \sum_{l=1}^{N} \int \left( \frac{p}{\gamma-1} + \frac{B^2}{2} \right) dv,\\
!latex           \ds \frac{\partial Z_j}{\partial t} & \equiv & \ds - \frac{\partial }{\partial Z_j} \sum_{l=1}^{N} \int \left( \frac{p}{\gamma-1} + \frac{B^2}{2} \right) dv.
!latex       \end{array} \label{eq:descent} \ee

!latex \item There remain degrees of freedom in the angle representation of the interfaces.

!latex \subsubsection{Spectral energy minimization}

!latex \item Consider variations which do not affect the geometry of the surfaces,
!latex       \be \delta R &=& R_\t \; u,\\
!latex           \delta Z &=& Z_\t \; u,
!latex       \ee
!latex       where $u$ is a angle variation.
!latex \item The corresponding variation in each of the Fourier harmonics is
!latex       \be \delta R_j &\equiv& \ooint R_\t \; u \; \cos \alpha_j,\\
!latex           \delta Z_j &\equiv& \ooint Z_\t \; u \; \sin \alpha_j,
!latex       \ee
!latex \item Following Hirshman et al., introducing the normalized spectral width
!latex       \be M \equiv \frac{\sum_j ( m_j^p + n_j^q ) ( R_{l,j}^2+Z_{l,j}^2 )}{\sum_j ( R_{l,j}^2+Z_{l,j}^2 )},
!latex       \ee
!latex \item Using the notation
!latex       \be N &\equiv& \sum_j \lambda_j ( R_{l,j}^2+Z_{l,j}^2 ), \\
!latex           D &\equiv& \sum_j           ( R_{l,j}^2+Z_{l,j}^2 ),
!latex       \ee
!latex       where $\lambda_j \equiv m_j^p + n_j^q$,
!latex       the variation in the normalized spectral width is 
!latex       \be \delta M = (\delta N - M \delta D)/D.
!latex       \ee
!latex \item For tangential variations, 
!latex       \be \delta N &=& 2 \ooint d\t d\z \; u \left(R_\t \sum_j \lambda_j R_j \cos \alpha_j + Z_\t \sum_j \lambda_j Z_j \sin \alpha_j \right),\\
!latex           \delta D &=& 2 \ooint d\t d\z \; u \left(R_\t \sum_j           R_j \cos \alpha_j + Z_\t \sum_j           Z_j \sin \alpha_j \right).
!latex       \ee

!latex \item The ``tangential spectral-width descent direction'' is thus 
!latex       \be \frac{\partial u}{\partial t} &=& -\left[ R_\t \sum_j(\lambda_j - M) R_j \cos \alpha_j / D + Z_\t \sum_j(\lambda_j - M)Z_j \sin \alpha_j / D \right].
!latex       \ee

!latex \item This suggests that position should be advanced according to
!latex       \be \frac{\partial R_j}{\partial t} & \equiv &
!latex       - \frac{\partial }{\partial R_j} \sum_{l=1}^{N} \int \left( \frac{p}{\gamma-1} + \frac{B^2}{2} \right) dv - [R_\t (R_\t X + Z_\t Y)]_j,\\
!latex           \frac{\partial Z_j}{\partial t} & \equiv &
!latex       - \frac{\partial }{\partial Z_j} \sum_{l=1}^{N} \int \left( \frac{p}{\gamma-1} + \frac{B^2}{2} \right) dv - [Z_\t (R_\t X + Z_\t Y)]_j,
!latex       \ee
!latex       where $X \equiv \sum_j (\lambda_j - M)R_j \cos\alpha_j / D$ and $Y \equiv \sum_j (\lambda_j - M)Z_j \sin\alpha_j / D$.

!latex \subsubsection{numerical implementation}

!latex \item The spectral condensation terms,
!latex \be R_\t (R_\t X + Z_\t Y) &=& \sum_{j,k,l} m_j m_k (\lambda_l-M) R_j ( + R_k R_l \sin\alpha_j\sin\alpha_k\cos\alpha_l - Z_k Z_l \sin\alpha_j\cos\alpha_k\sin\alpha_l)/D,\\
!latex     Z_\t (R_\t X + Z_\t Y) &=& \sum_{j,k,l} m_j m_k (\lambda_l-M) Z_j ( - R_k R_l \cos\alpha_j\sin\alpha_k\cos\alpha_l + Z_k Z_l \cos\alpha_j\cos\alpha_k\sin\alpha_l)/D,
!latex \ee
!latex are calculated using triple angle expressions:
!latex IT IS VERY LIKELY THAT FFTs WOULD BE FASTER!!!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine pc00ab( mode, Ngeometricaldof, Position, Energy, Gradient, nstate, iuser, ruser ) ! argument fixed by NAG; see pc00aa;
  
  use constants, only : zero, half, one
  
  use numerical, only :
  
  use fileunits, only : ounit
  
  use inputlist, only : Wpc00ab, Igeometry, Nvol, epsilon, maxiter, ForceErr, forcetol
  
  use cputiming, only : Tpc00ab
  
  use allglobal, only : writin, myid, cpus, YESstellsym, mn, lBBintegral, dBBdRZ, dIIdRZ
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER :: mode, Ngeometricaldof, nstate, iuser(1:2)
  REAL    :: Position(1:Ngeometricaldof), Energy, Gradient(1:Ngeometricaldof), ruser(1:1)
  
  LOGICAL :: LComputeDerivatives
  INTEGER :: ii, vvol, irz, issym, totaldof, localdof, wflag, iflag!, mi, ni !idof, imn, irz, totaldof, localdof, jj, kk, ll, mi, ni, mj, nj, mk, nk, ml, nl, mjmk
  REAL    :: force(0:Ngeometricaldof), gradienterror, rflag
  
  BEGIN(pc00ab)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  iuser(1) = iuser(1) + 1 ! iteration counter;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LComputeDerivatives = .false.
  WCALL(pc00ab,fc02aa,( Ngeometricaldof, Position(1:Ngeometricaldof), force(0:Ngeometricaldof), LComputeDerivatives ))

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  Energy = sum( lBBintegral(1:Nvol) ) ! Energy must always be assigned; 26 Feb 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( nstate.eq.1 ) ruser(1) = Energy ! this is the first call; 26 Feb 13;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  select case( mode ) 
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
  case( 0 ) ! need to assign Energy; as Energy must always be assigned, it is assigned above select case; 26 Feb 13;

   gradienterror = zero
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
  case( 2 ) ! need to assign Energy & Gradient;
   
   iuser(2) = iuser(2) + 1 ! iteration counter;
   
   totaldof = 0 ! total degree-of-freedom counter;
   
   do vvol = 1, Nvol-1
    
    localdof = 0
    
    do ii = 1, mn !; mi = im(ii) ; ni = in(ii) ! loop over Fourier harmonics; 26 Feb 13;
     
     do irz = 0, 1
      
      if( irz.eq.1 .and. Igeometry.lt.3 ) cycle ! Z is not required for Cartesian and standard-cylindrical; 04 Dec 14;
      
      do issym = 0, 1
       
       if( issym.eq.1 .and. YESstellsym ) cycle
       
       if( ii.eq.1 .and. irz.eq.1 .and. issym.eq.0 ) cycle ! no dependence on Z_{m,n} \sin( m\t - n\z ) for m=0, n=0;
       if( ii.eq.1 .and. irz.eq.0 .and. issym.eq.1 ) cycle ! no dependence on Z_{m,n} \sin( m\t - n\z ) for m=0, n=0;
       
       localdof = localdof + 1
       
       totaldof = totaldof + 1

       if( Igeometry.lt.3 ) then ; Gradient(totaldof) = ( dBBdRZ(vvol,1,localdof) + dBBdRZ(vvol+1,0,localdof) ) ! no spectral constraints; 04 Dec 14;
       else                      ; Gradient(totaldof) = ( dBBdRZ(vvol,1,localdof) + dBBdRZ(vvol+1,0,localdof) ) - epsilon * dIIdRZ(vvol,localdof) ! computed in fc02aa; 26 Feb 13;
       endif
       
      enddo ! end of do issym; 26 Feb 13;
      
     enddo ! end of do irz; 26 Feb 13;
     
    enddo ! end of do ii; 26 Feb 13;
    
   enddo ! end of do vvol; 26 Feb 13;
   
   FATALMESS(pc00ab, totaldof.ne.Ngeometricaldof, counting error )
   
   gradienterror = sum( abs( Gradient(1:Ngeometricaldof) ) ) / Ngeometricaldof ! only used for screen output; 26 Feb 13;
   
   wflag = 1 ; iflag = 0 ; rflag = gradienterror
   WCALL(pc00ab,writin,( wflag, iflag, rflag)) ! write restart file etc.;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
  case( 3 ) ! second derivatives; required for E04LYF, which is called by pc02aa;
   
   FATALMESS(pc00ab, .true., have not yet computed second derivatives )
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
  case default
    
   FATALMESS(pc00ab, .true., invalid mode provided to pc00ab )
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
  end select

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  cput = GETTIME
  if( myid.eq.0 ) write(ounit,1000) cput-cpus, iuser(1:2), mode, Energy, gradienterror, ForceErr
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( iuser(2).ge.maxiter .and. maxiter.gt.0 ) mode = -1 ! iteration limit reached; E04DGF will terminate with ifail=mode;
  if( ForceErr.lt.abs(forcetol)              ) mode = -2 ! force balance satisfied; E04DGF will terminate with ifail=mode;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  Energy                      = Energy                      / ruser(1) ! normalize to initial energy; 26 Feb 13;
  Gradient(1:Ngeometricaldof) = Gradient(1:Ngeometricaldof) / ruser(1)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(pc00ab)
  
1000 format("pc00ab : ",f10.2," : iterations="2i8" ; mode=",i3," ; Energy="es23.15" ; |DF|="es13.5" ; ForceErr="es23.15" ;")
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine pc00ab

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
