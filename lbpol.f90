!> \file
!> \brief Computes \f$B_{\theta,e,0,0}\f$ at the interface

!> \brief Computes \f$B_{\theta,e,0,0}\f$ at the interface
!> @param[in] lvol
!> @param[in, out] Bt00
!> @param[in] ideriv
!> @param[in] iocons
subroutine lbpol(lvol, Bt00, ideriv, iocons)

  use constants, only : mu0, pi, pi2, two, one, half, zero

  use allglobal, only : Ate, Aze, Ato, Azo, TT, &
                        YESstellsym, NOTstellsym, &
                        im, in, mne, ime, ine, Mvol, mn, &
                        sg, guvij, &
                        Ntz, Lcoordinatesingularity, &
                        efmn, ofmn, cfmn, sfmn, evmn, odmn, comn, simn, &
                        Nt, Nz, &
                        cpus, myid, dBdX, &
                        build_vector_potential

  use inputlist, only : Lrad, Wlbpol, Igeometry, Lcheck

  use fileunits, only : ounit

  use cputiming, only : Tlbpol

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS
! ------

  INTEGER, intent(in)    :: ideriv      !< lbpol will return \f$B_{\theta,e,0,0}\f$ (0) or its derivative with respect to the geometry (-1), mu (1) or the poloidal flux (2). \f$\mbox{ideriv}\in\{-1,\ldots,2\} \f$
  INTEGER, intent(in)    :: lvol        !< Volume index. \f$\mbox{lvol}\in\{1,\ldots,\mbox{Mvol}\} \f$
  INTEGER, intent(in)    :: iocons      !< \f$B_{\theta,e,0,0}\f$ is evaluated on the inner (iocons=0) or outer (iocons=1) volume boundary. \f$\mbox{iocons}\in\{0,1\} \f$
  REAL,    intent(inout) :: Bt00(1:Mvol, 0:1, -1:2) !< \f$B_{\theta,e,0,0}\f$, with indices Bt00( lvol, iocons, ideriv ).
  INTEGER   :: Lcurvature  !< Input to coords(). Specify what geometrical quantity is returned.
  INTEGER   :: ifail       !< Dummy variable used when tfft() is called.
  REAL      :: dAt(1:Ntz)  !< Poloidal vector potential \f$ A_\theta \f$, size 1:Ntz. If \f$ideriv\neq 0\f$, derivative of \f$ A_\theta \f$.
  REAL      :: dAz(1:Ntz)  !< Toroidal vector potential \f$ A_\phi   \f$, size 1:Ntz. If \f$ideriv\neq 0\f$, derivative of \f$ A_\phi \f$.
  REAL      :: Bt(1:Ntz)   !< Poloidal magnetic field \f$ B_\theta \f$, size 1:Ntz. If \f$ideriv\neq 0\f$, derivative of \f$ B_\theta \f$.
  REAL      :: Bz(1:Ntz)   !< Toroidal magnetic field \f$ B_\phi \f$, size 1:Ntz. If \f$ideriv\neq 0\f$, derivative of \f$ B_\phi \f$.
  REAL      :: dAt0(1:Ntz) !< Poloidal vector potential \f$ A_\theta \f$, size 1:Ntz. Only used if \f$ ideriv=-1 \f$.
  REAL      :: dAz0(1:Ntz) !< Toroidal vector potential \f$ A_\phi \f$, size 1:Ntz. Only used if \f$ ideriv=-1 \f$.
  REAL      :: lss         !< s-coordinate

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> Computes \f$B_{\theta,e,0,0}\f$ at the volume interfaces. This is used by dfp100 to evaluate the toroidal current at the volume
!> interfaces, and by dfp200 to construct the force gradient when the current constraint (Lconstraint=3) is used. This is also used
!> by xspech to compute the toroidal current at the volume interfaces, written in the output.
!> <ol>

  BEGIN(lbpol)

! TODO: This subroutine is very similar to curent.f90 - maybe merge both in a single subroutine to simplify?


!> <li> Call coords() to compute the metric coefficients and the jacobian. </li>
  lss = two * iocons - one ! Build s coordinate

  Lcurvature = 1
  WCALL( lbpol, coords, (lvol, lss, Lcurvature, Ntz, mn ) )


!> <li> Build coefficients efmn, ofmn, cfmn, sfmn from the field vector potential Ate, Ato,
!>      Aze and Azo, and radial derivatives of the polynomial basis TT(ll,innout,1). These variables
!>      are the derivatives with respect to s of the magnetic field vector potential in Fourier space.
!>      If ideriv\f$\neq 0\f$, construct the relevant derivatives of the vector potential. </li>
  call build_vector_potential(lvol, iocons, ideriv, 1)

!> <li> Take the inverse Fourier transform of efmn, ofmn, cfmn, sfmn. These are the covariant components of \f$\frac{\partial A}{\partial s}\f$,
!>      <em>i.e.</em> the contravariant components of \f$\mathbf{B}\f$. </li>
  call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, dAt(1:Ntz), dAz(1:Ntz) )

!> <li> Build covariant components of the field using the metric coefficients guvij and the jacobian sg. </li>
  Bt(1:Ntz) = ( - dAz(1:Ntz ) * guvij(1:Ntz,2,2, 0) + dAt(1:Ntz ) * guvij(1:Ntz,2,3, 0) )/ sg(1:Ntz,0)
  Bz(1:Ntz) = ( - dAz(1:Ntz ) * guvij(1:Ntz,2,3, 0) + dAt(1:Ntz ) * guvij(1:Ntz,3,3, 0) )/ sg(1:Ntz,0)


!> <li> If ideriv=-1 (derivatives with respect to the geometry), need to add derivatives relative to the metric elements </li>
  if( ideriv.eq.-1 ) then

!> <ol>
!> <li> Get derivatives of metric element by calling coords() </li>
    Lcurvature = 3
    WCALL( lbpol, coords, (lvol, lss, Lcurvature, Ntz, mn ) ) ! get sg times d/dx (g_mu,nu / sg)

!> <li> Compute vector potential without taking any derivatives </li>
    call build_vector_potential(lvol, iocons, 0, 1)
    call invfft( mn, im, in, efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, dAt0(1:Ntz), dAz0(1:Ntz) ) ! get covariant component of dA without derivatives

!> <li> Add to \f$\frac{\partial B_\theta}{\partial x_i}\f$ the contributions from \f$ \frac{\partial}{\partial x_i} \frac{g_{\mu\nu}}{\sqrt{g}} \f$
    Bt(1:Ntz) = Bt(1:Ntz) + ( - dAz0(1:Ntz ) * guvij(1:Ntz,2,2, 1) + dAt0(1:Ntz ) * guvij(1:Ntz,2,3, 1) ) / sg(1:Ntz, 0)
    Bz(1:Ntz) = Bz(1:Ntz) + ( - dAz0(1:Ntz ) * guvij(1:Ntz,2,3, 1) + dAt0(1:Ntz ) * guvij(1:Ntz,3,3, 1) ) / sg(1:Ntz, 0)

!> </ol>
  endif

!> <li> Fourier transform the field and store it in the variables efmn, ofmn, cfmn and sfmn. </li>
  ifail = 0
  call tfft( Nt, Nz, Bt(1:Ntz), Bz(1:Ntz), mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), ifail )

!> <li> Save first even fourier mode into Bt00( lvol, iocons, ideriv ) </li>
  Bt00(lvol, iocons, ideriv) = efmn(1)


!> </ol>
  ! End subroutine
  RETURN(lbpol)

end subroutine lbpol






