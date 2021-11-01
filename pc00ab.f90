!> \file
!> \brief Returns the energy functional and it's derivatives with respect to geometry.

!> \brief Returns the energy functional and it's derivatives with respect to geometry.
!> \ingroup grp_conjugate_gradient
!>
!> **Energy functional**
!>
!> <ul>
!> <li> The energy functional is
!>      \f{eqnarray}{ F \equiv \sum_{l=1}^{N} \int_{\cal V} \left( \frac{p}{\gamma-1} + \frac{B^2}{2} \right) dv,
!>      \label{eq:energyfunctional_pc00ab}
!>      \f}
!>      where \f$N \equiv\,\f$\c Nvol is the number of interfaces. </li>
!> <li> Assuming that the toroidal and poloidal fluxes, \f$\psi_t\f$ and \f$\psi_p\f$, the helicity, \f${\cal K}\f$, the helicity multiplier, \f$\mu\f$,
!>      and/or the interface rotational-transforms, \f${{\,\iota\!\!\!}-}\f$, are appropriately constrained,
!>      the Beltrami fields in each volume depend only the geometry of the adjacent interfaces.
!>      So, the energy functional is assumed to be a function of "position", i.e. \f$F = F(R_{l,j},Z_{l,j})\f$. </li>
!> <li> Introducing a ficitious time, \f$t\f$, the position may be advanced according to
!>       \f{eqnarray}{ \begin{array}{cccccccccccccccccccccccc}
!>                     \displaystyle \frac{\partial R_j}{\partial t} & \equiv & \displaystyle
!>                     - \frac{\partial }{\partial R_j} \sum_{l=1}^{N} \int \left( \frac{p}{\gamma-1} + \frac{B^2}{2} \right) dv,\\
!>                     \displaystyle \frac{\partial Z_j}{\partial t} & \equiv & \displaystyle
!>                     - \frac{\partial }{\partial Z_j} \sum_{l=1}^{N} \int \left( \frac{p}{\gamma-1} + \frac{B^2}{2} \right) dv.
!>       \end{array} \label{eq:descent_pc00ab} \f} </li>
!> <li> There remain degrees of freedom in the angle representation of the interfaces. </li>
!> </ul>
!>
!> **Spectral energy minimization**
!>
!> <ul>
!> <li> Consider variations which do not affect the geometry of the surfaces,
!>       \f{eqnarray}{ \delta R &=& R_\theta \; u,\\
!>                     \delta Z &=& Z_\theta \; u,
!>       \f}
!>       where \f$u\f$ is a angle variation. </li>
!> <li> The corresponding variation in each of the Fourier harmonics is
!>       \f{eqnarray}{ \delta R_j &\equiv& \oint\!\!\!\oint \!d\theta d\zeta \,\,\, R_\theta \; u \; \cos \alpha_j,\\
!>                     \delta Z_j &\equiv& \oint\!\!\!\oint \!d\theta d\zeta \,\,\, Z_\theta \; u \; \sin \alpha_j,
!>       \f} </li>
!> <li> Following Hirshman et al., introducing the normalized spectral width
!>       \f{eqnarray}{ M \equiv \frac{\sum_j ( m_j^p + n_j^q ) ( R_{l,j}^2+Z_{l,j}^2 )}{\sum_j ( R_{l,j}^2+Z_{l,j}^2 )},
!>       \f} </li>
!> <li> Using the notation
!>       \f{eqnarray}{ N &\equiv& \sum_j \lambda_j ( R_{l,j}^2+Z_{l,j}^2 ), \\
!>                     D &\equiv& \sum_j           ( R_{l,j}^2+Z_{l,j}^2 ),
!>       \f}
!>       where \f$\lambda_j \equiv m_j^p + n_j^q\f$,
!>       the variation in the normalized spectral width is
!>       \f{eqnarray}{ \delta M = (\delta N - M \delta D)/D.
!>       \f} </li>
!> <li> For tangential variations,
!>       \f{eqnarray}{ \delta N &=& 2 \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \; u \left(R_\theta \sum_j \lambda_j R_j \cos \alpha_j + Z_\theta \sum_j \lambda_j Z_j \sin \alpha_j \right),\\
!>                     \delta D &=& 2 \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \; u \left(R_\theta \sum_j           R_j \cos \alpha_j + Z_\theta \sum_j           Z_j \sin \alpha_j \right).
!>       \f} </li>
!> <li> The "tangential spectral-width descent direction" is thus
!>       \f{eqnarray}{ \frac{\partial u}{\partial t} &=& -\left[ R_\theta \sum_j(\lambda_j - M) R_j \cos \alpha_j / D + Z_\theta \sum_j(\lambda_j - M)Z_j \sin \alpha_j / D \right].
!>       \f} </li>
!> <li> This suggests that position should be advanced according to
!>       \f{eqnarray}{ \frac{\partial R_j}{\partial t} & \equiv &
!>       - \frac{\partial }{\partial R_j} \sum_{l=1}^{N} \int \left( \frac{p}{\gamma-1} + \frac{B^2}{2} \right) dv - [R_\theta (R_\theta X + Z_\theta Y)]_j,\\
!>           \frac{\partial Z_j}{\partial t} & \equiv &
!>       - \frac{\partial }{\partial Z_j} \sum_{l=1}^{N} \int \left( \frac{p}{\gamma-1} + \frac{B^2}{2} \right) dv - [Z_\theta (R_\theta X + Z_\theta Y)]_j,
!>       \f}
!>       where \f$X \equiv \sum_j (\lambda_j - M)R_j \cos\alpha_j / D\f$ and \f$Y \equiv \sum_j (\lambda_j - M)Z_j \sin\alpha_j / D\f$. </li>
!> </ul>
!>
!> **numerical implementation**
!>
!> <ul>
!> <li> The spectral condensation terms,
!> \f{eqnarray}{ R_\theta (R_\theta X + Z_\theta Y) &=& \sum_{j,k,l} m_j m_k (\lambda_l-M) R_j ( + R_k R_l \sin\alpha_j\sin\alpha_k\cos\alpha_l - Z_k Z_l \sin\alpha_j\cos\alpha_k\sin\alpha_l)/D,\\
!>               Z_\theta (R_\theta X + Z_\theta Y) &=& \sum_{j,k,l} m_j m_k (\lambda_l-M) Z_j ( - R_k R_l \cos\alpha_j\sin\alpha_k\cos\alpha_l + Z_k Z_l \cos\alpha_j\cos\alpha_k\sin\alpha_l)/D,
!> \f}
!> are calculated using triple angle expressions...
!> \todo IT IS VERY LIKELY THAT FFTs WOULD BE FASTER!!!
!>
!> </li>
!> </ul>
!>
subroutine pc00ab( mode, NGdof, Position, Energy, Gradient, nstate, iuser, ruser ) ! argument fixed by NAG; see pc00aa;

  use constants, only : zero, half, one

  use numerical, only :

  use fileunits, only : ounit

  use inputlist, only : Wpc00ab, Igeometry, Nvol, epsilon, maxiter, forcetol

  use cputiming, only : Tpc00ab

  use allglobal, only : writin, myid, cpus, YESstellsym, mn_field, lBBintegral, dBBdRZ, dIIdRZ, ForceErr

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS

  INTEGER :: mode, NGdof, nstate, iuser(1:2)
  REAL    :: Position(1:NGdof), Energy, Gradient(1:NGdof), ruser(1:1)

  LOGICAL :: LComputeDerivatives, LComputeAxis
  INTEGER :: ii, vvol, irz, issym, totaldof, localdof, wflag, iflag!, mi, ni !idof, imn, irz, totaldof, localdof, jj, kk, ll, mi, ni, mj, nj, mk, nk, ml, nl, mjmk
  REAL    :: force(0:NGdof), gradienterror, rflag

  BEGIN(pc00ab)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  iuser(1) = iuser(1) + 1 ! iteration counter;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LComputeDerivatives = .false.
  LComputeAxis = .true.
  WCALL(pc00ab,dforce,( NGdof, Position(1:NGdof), force(0:NGdof), LComputeDerivatives, LComputeAxis ))

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  Energy = sum( lBBintegral(1:Nvol) ) ! Energy must always be assigned; 26 Feb 13;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( nstate.eq.1 ) ruser(1) = Energy ! this is the first call; 26 Feb 13;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  select case( mode )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  case( 0 ) ! need to assign Energy; as Energy must always be assigned, it is assigned above select case; 26 Feb 13;

   gradienterror = zero

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  case( 2 ) ! need to assign Energy & Gradient;

   iuser(2) = iuser(2) + 1 ! iteration counter;

   totaldof = 0 ! total degree-of-freedom counter;

   do vvol = 1, Nvol-1

    localdof = 0

    do ii = 1, mn_field !; mi = im(ii) ; ni = in(ii) ! loop over Fourier harmonics; 26 Feb 13;

     do irz = 0, 1

      if( irz.eq.1 .and. Igeometry.lt.3 ) cycle ! Z is not required for Cartesian and standard-cylindrical; 04 Dec 14;

      do issym = 0, 1

       if( issym.eq.1 .and. YESstellsym ) cycle

       if( ii.eq.1 .and. irz.eq.1 .and. issym.eq.0 ) cycle ! no dependence on Z_{m,n} \sin( m\t - n\z ) for m=0, n=0;
       if( ii.eq.1 .and. irz.eq.0 .and. issym.eq.1 ) cycle ! no dependence on Z_{m,n} \sin( m\t - n\z ) for m=0, n=0;

       localdof = localdof + 1

       totaldof = totaldof + 1

       if( Igeometry.lt.3 ) then ; Gradient(totaldof) = ( dBBdRZ(vvol,1,localdof) + dBBdRZ(vvol+1,0,localdof) ) ! no spectral constraints; 04 Dec 14;
       else                      ; Gradient(totaldof) = ( dBBdRZ(vvol,1,localdof) + dBBdRZ(vvol+1,0,localdof) ) - epsilon * dIIdRZ(vvol,localdof) ! computed in dforce; 26 Feb 13;
       endif

      enddo ! end of do issym; 26 Feb 13;

     enddo ! end of do irz; 26 Feb 13;

    enddo ! end of do ii; 26 Feb 13;

   enddo ! end of do vvol; 26 Feb 13;

   FATAL(pc00ab, totaldof.ne.NGdof, counting error )

   gradienterror = sum( abs( Gradient(1:NGdof) ) ) / NGdof ! only used for screen output; 26 Feb 13;

   wflag = 1 ; iflag = 0 ; rflag = gradienterror
   WCALL(pc00ab,writin,( wflag, iflag, rflag)) ! write restart file etc.;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  case( 3 ) ! second derivatives; required for E04LYF, which is called by pc02aa;

   FATAL(pc00ab, .true., have not yet computed second derivatives )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  case default

   FATAL(pc00ab, .true., invalid mode provided to pc00ab )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  end select

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  cput = GETTIME
  if( myid.eq.0 ) write(ounit,1000) cput-cpus, iuser(1:2), mode, Energy, gradienterror, ForceErr

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( iuser(2).ge.maxiter .and. maxiter.gt.0 ) mode = -1 ! iteration limit reached; E04DGF will terminate with ifail=mode;
  if( ForceErr.lt.abs(forcetol)              ) mode = -2 ! force balance satisfied; E04DGF will terminate with ifail=mode;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  Energy                      = Energy                      / ruser(1) ! normalize to initial energy; 26 Feb 13;
  Gradient(1:NGdof) = Gradient(1:NGdof) / ruser(1)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(pc00ab)

1000 format("pc00ab : ",f10.2," : iterations="2i8" ; mode=",i3," ; Energy="es23.15" ; |DF|="es13.5" ; ForceErr="es23.15" ;")

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine pc00ab

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!latex \end{enumerate}

