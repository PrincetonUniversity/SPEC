!> \file
!> \brief Polynomials evaluation

!> \brief Get the Chebyshev polynomials with zeroth, first derivatives
!>
!> The Chebyshev polynomial has been recombined and rescaled.
!> By doing so, the Chebyshev polynomial satisfy the zero Dirichlet boundary condition on the inner surface of the annulus with reduced ill-conditioning problem.
!>
!> Let \f$T_{l}\f$ be the Chebyshev polynomial of the first kind with degree \f$l\f$.
!> This subroutine computes
!> \f[\bar{T}_0 = 1, \f] and
!> \f[\bar{T}_l = \frac{T_l - (-1)^{l}}{l+1} . \f]
!>
!> \f$ T_l \f$ are computed iteratively.
!> \f[ T_0(s) = 1, \f]
!> \f[ T_1(s) = s, \f]
!> \f[ T_{l+1}(s) = 2 s T_l(s) - T_{l-1}(s). \f]
!>
!> @param[in] lss coordinate input lss
!> @param[in] lrad radial resolution
!> @param[out] cheby the value, first derivative of Chebyshev polynomial
subroutine get_cheby(lss, lrad, cheby)

  use constants, only : zero, one, two

  implicit none

  REAL,intent(in) :: lss
  INTEGER, intent(in) :: lrad
  REAL, intent(inout) :: cheby(0:lrad,0:1)

  integer :: ll

  cheby = zero

    ;                 cheby( 0,0:1) = (/ one                                       , zero                                                            /)
    ;                 cheby( 1,0:1) = (/ lss                                       , one                                                             /)
    do ll = 2, lrad ; cheby(ll,0:1) = (/ two * lss * cheby(ll-1,0) - cheby(ll-2,0) , two * cheby(ll-1,0) + two * lss * cheby(ll-1,1) - cheby(ll-2,1) /)
    enddo

  ! basis recombination
  do ll = 1, lrad
    cheby(ll, 0) = cheby(ll, 0)  - (-1)**ll
  enddo

  do ll = 0, lrad
    cheby(ll, 0:1) = cheby(ll, 0:1) / real(ll+1) ! scale for better conditioning
  enddo

  return
end subroutine get_cheby

!> \brief Get the Chebyshev polynomials with zeroth, first and second derivatives
!> The Chebyshev polynomial has been recombined and rescaled.
!> See get_cheby for more detail.
!> @param[in] lss coordinate input lss
!> @param[in] lrad radial resolution
!> @param[out] cheby the value, first and second derivative of Chebyshev polynomial
subroutine get_cheby_d2(lss, lrad, cheby)

  use constants, only : zero, one, two

  implicit none

  REAL,intent(in) :: lss
  INTEGER, intent(in) :: lrad
  REAL, intent(inout) :: cheby(0:lrad,0:2)

  integer :: ll

  cheby = zero

  ;                     ; cheby( 0,0:2) = (/ one, zero, zero /) ! T_0: Chebyshev initialization; function, 1st-derivative, 2nd-derivative;
  ;                     ; cheby( 1,0:2) = (/ lss,  one, zero /) ! T_1: Chebyshev initialization; function, 1st-derivative, 2nd-derivative;
  do ll = 2, lrad
    cheby(ll,0:2) = (/ two * lss * cheby(ll-1,0)                                                         - cheby(ll-2,0) , &
                       two       * cheby(ll-1,0) + two * lss * cheby(ll-1,1)                             - cheby(ll-2,1) , &
                       two       * cheby(ll-1,1) + two       * cheby(ll-1,1) + two * lss * cheby(ll-1,2) - cheby(ll-2,2) /)
  enddo

  do ll = 1, lrad
    cheby(ll, 0) = cheby(ll, 0)  - (-1)**ll
  enddo

  do ll = 0, lrad
    cheby(ll, 0:2) = cheby(ll, 0:2) / real(ll+1) ! scale for better conditioning
  enddo

  return
end subroutine get_cheby_d2

!> \brief Get the Zernike polynomials \f$\hat{R}^{m}_{l}\f$ with zeroth, first derivatives
!>
!> The original Zernike polynomial is defined by
!> The Zernike polynomials take the form
!> \f{eqnarray*}{
!>    Z^{-m}_{l}(s,\theta) &= R^{m}_{l}(s) \sin m\theta,\\
!>    Z^{m}_{l}(s,\theta) &= R^{m}_{l}(s) \cos m\theta,
!> \f}
!> where \f$R^{m}_{l}(s)\f$ is a \f$l\f$-th order polynomial given by
!> \f{eqnarray*}{
!>    R_l^m(s) = \sum^{\frac{l-m}{2}}_{k=0} \frac{(-1)^k (l-k)!}
!>    {k!\left[\frac{1}{2}(l+m) - k \right]! \left[\frac{1}{2}(l-m) - k \right]!} s^{l-2k},
!> \f}
!> and is only non-zero for \f$l \ge m\f$ and even \f$l-m\f$.
!>
!> In this subroutine, \f$R^{m}_{l}(s) \f$ is computed using the iterative relationship
!> \f{eqnarray*}{
!> R_l^m(s) = \frac{2(l-1)(2l(l-2)s^2 -m^2 -l(l-2))R_{l-2}^m(s) - l (l+m-2) (l-m-2)R_{l-4}^m(s) }{(l+m)(l-m)(l-2)}
!> \f}
!>
!> For \f$ m=0 \f$ and \f$ m=1 \f$, a basis recombination method is used by defining new radial basis functions as
!> \f{eqnarray*}{
!>    \hat{R}_{0}^{0} &= 1,
!>    \hat{R}^{0}_{l} &= \frac{1}{l+1}R^{0}_{l} - \frac{(-1)^{l/2}}{l+1},
!>    \\
!>    \hat{R}_{1}^{1} &= s,
!>    \hat{R}^{1}_{l} &= \frac{1}{l+1}R^{1}_{l} - \frac{(-1)^{(l-1)/2}}{2} s.
!> \f}
!> so that the basis scales as \f$s^{m+2}\f$ except for \f$\hat{R}_{0}^{0}\f$ and \f$\hat{R}_{1}^{1}\f$, which are excluded from the representation of \f${A}_{\theta,m,n}\f$.
!> For \f$m\ge2\f$, the radial basis functions are only rescaled as
!> \f[
!>   \hat{R}^{m}_{l} = \frac{1}{l+1}R^{m}_{l}.
!> \f]
!>
!> @param[in] r coordinate input, note that this is normalized to \f$[0, 1]\f$
!> @param[in] lrad radial resolution
!> @param[in] mpol poloidal resolution
!> @param[out] zernike the value, first derivative of Zernike polynomial
subroutine get_zernike(r, lrad, mpol, zernike)

  use constants, only : zero, one, two

  implicit none

  REAL,intent(in) :: r
  INTEGER, intent(in) :: lrad, mpol
  REAL, intent(inout) :: zernike(0:lrad,0:mpol,0:1)

  REAL ::    rm, rm1  ! r to the power of m'th and m-1'th
  REAL ::    factor1, factor2, factor3, factor4
  INTEGER :: m, n  ! Zernike R^m_n

  rm = one  ! r to the power of m'th
  rm1 = zero ! r to the power of m-1'th
  zernike(:,:,:) = zero
  do m = 0, mpol
    if (lrad >= m) then
      zernike(m,m,0:1) = (/ rm, real(m)*rm1 /)
    endif

    if (lrad >= m+2) then
      zernike(m+2,m,0) = real(m+2)*rm*r**2 - real(m+1)*rm
      zernike(m+2,m,1) = real((m+2)**2)*rm*r - real((m+1)*m)*rm1
    endif

    do n = m+4, lrad, 2
      factor1 = real(n)/real(n**2 - m**2)
      factor2 = real(4 * (n-1))
      factor3 = real((n-2+m)**2)/real(n-2) + real((n-m)**2)/real(n)
      factor4 = real((n-2)**2-m**2) / real(n-2)

      zernike(n, m, 0) = factor1 * ((factor2*r**2 - factor3)*zernike(n-2,m,0) - factor4*zernike(n-4,m,0))
      zernike(n, m, 1) = factor1 * (two*factor2*r*zernike(n-2,m,0) + (factor2*r**2 - factor3)*zernike(n-2,m,1) - factor4*zernike(n-4,m,1))
    enddo

    rm1 = rm
    rm = rm * r

  enddo

  do n = 2, lrad, 2
    zernike(n,0,0) = zernike(n,0,0) - (-1)**(n/2)
  enddo

  if (mpol >= 1) then
    do n = 3, lrad, 2
      zernike(n,1,0) = zernike(n,1,0) - (-1)**((n-1)/2) * real((n+1)/2) * r
      zernike(n,1,1) = zernike(n,1,1) - (-1)**((n-1)/2) * real((n+1)/2)
    enddo
  end if

  do m = 0, mpol
    do n = m, lrad, 2
      zernike(n,m,:) = zernike(n,m,:) / real(n+1)
    end do
  end do
end subroutine get_zernike

!> \brief Get the Zernike polynomials  \f$\hat{R}^{m}_{l}\f$ with zeroth, first, second derivatives
!>
!> See get_zernike for more detail.
!>
!> @param[in] r coordinate input, note that this is normalized to \f$[0, 1]\f$
!> @param[in] lrad radial resolution
!> @param[in] mpol poloidal resolution
!> @param[out] zernike the value, first/second derivative of Zernike polynomial
subroutine get_zernike_d2(r, lrad, mpol, zernike)

  use constants, only : zero, one, two

  implicit none

  REAL,intent(in) :: r
  INTEGER, intent(in) :: lrad, mpol
  REAL, intent(inout) :: zernike(0:lrad,0:mpol,0:2)

  REAL ::    rm, rm1, rm2  ! r to the power of m'th, m-1'th and m-2'th
  REAL ::    factor1, factor2, factor3, factor4
  INTEGER :: m, n  ! Zernike R^m_n

  rm = one  ! r to the power of m'th
  rm1 = zero ! r to the power of m-1'th
  rm2 = zero ! r to the power of m-2'th
  zernike(:,:,:) = zero
  do m = 0, mpol
    if (lrad >= m) then
      zernike(m,m,0:2) = (/ rm, real(m)*rm1, real(m*(m-1))*rm2 /)
      !write(0, *) m, m, r, zernike(m,m,:)
    endif

    if (lrad >= m+2) then
      zernike(m+2,m,0) = real(m+2)*rm*r**2 - real(m+1)*rm
      zernike(m+2,m,1) = real((m+2)**2)*rm*r - real((m+1)*m)*rm1
      zernike(m+2,m,2) = real((m+2)**2*(m+1))*rm - real((m+1)*m*(m-1))*rm2
      !write(0, *) m+2, m, r, zernike(m+2,m,:)
    endif

    do n = m+4, lrad, 2
      factor1 = real(n)/real(n**2 - m**2)
      factor2 = real(4 * (n-1))
      factor3 = real((n-2+m)**2)/real(n-2) + real((n-m)**2)/real(n)
      factor4 = real((n-2)**2-m**2) / real(n-2)

      zernike(n, m, 0) = factor1 * ((factor2*r**2 - factor3)*zernike(n-2,m,0) - factor4*zernike(n-4,m,0))
      zernike(n, m, 1) = factor1 * (two*factor2*r*zernike(n-2,m,0) + (factor2*r**2 - factor3)*zernike(n-2,m,1) - factor4*zernike(n-4,m,1))
      zernike(n, m, 2) = factor1 * (two*factor2*(two*r*zernike(n-2,m,1) + zernike(n-2,m,0)) &
                        +(factor2*r**2 - factor3)*zernike(n-2,m,2) - factor4*zernike(n-4,m,2))
      !write(0, *) n, m, r, zernike(n,m,:)
    enddo

    rm2 = rm1
    rm1 = rm
    rm = rm * r

  enddo
  do n = 2, lrad, 2
    zernike(n,0,0) = zernike(n,0,0) - (-1)**(n/2)
  enddo
  if (mpol >= 1) then
    do n = 3, lrad, 2
      zernike(n,1,0) = zernike(n,1,0) - (-1)**((n-1)/2) * real((n+1)/2) * r
      zernike(n,1,1) = zernike(n,1,1) - (-1)**((n-1)/2) * real((n+1)/2)
    enddo
  end if

  do m = 0, mpol
    do n = m, lrad, 2
      zernike(n,m,:) = zernike(n,m,:) / real(n+1)
    end do
  end do
end subroutine get_zernike_d2

!> \brief Get the Zernike polynomials \f$\hat{R}^{m}_{l}/r^m\f$
!>
!> See get_zernike for more detail.
!>
!> @param[in] r coordinate input, note that this is normalized to \f$[0, 1]\f$
!> @param[in] lrad radial resolution
!> @param[in] mpol poloidal resolution
!> @param[out] zernike the value
subroutine get_zernike_rm(r, lrad, mpol, zernike)

  use constants, only : zero, one, two

  implicit none

  REAL,intent(in) :: r
  INTEGER, intent(in) :: lrad, mpol
  REAL, intent(inout) :: zernike(0:lrad,0:mpol)

  REAL ::    factor1, factor2, factor3, factor4
  INTEGER :: m, n  ! Zernike R^m_n

  zernike(:,:) = zero
  do m = 0, mpol
    if (lrad >= m) then
      zernike(m,m) = one
    endif

    if (lrad >= m+2) then
      zernike(m+2,m) = real(m+2)*r**2 - real(m+1)
    endif

    do n = m+4, lrad, 2
      factor1 = real(n)/real(n**2 - m**2)
      factor2 = real(4 * (n-1))
      factor3 = real((n-2+m)**2)/real(n-2) + real((n-m)**2)/real(n)
      factor4 = real((n-2)**2-m**2) / real(n-2)

      zernike(n, m) = factor1 * ((factor2*r**2 - factor3)*zernike(n-2,m) - factor4*zernike(n-4,m))
    enddo

  enddo
  do n = 2, lrad, 2
    zernike(n,0) = zernike(n,0) - (-1)**(n/2)
  enddo
  if (mpol >= 1) then
    do n = 3, lrad, 2
      zernike(n,1) = zernike(n,1) - (-1)**((n-1)/2) * real((n+1)/2)
    enddo
  end if

  do m = 0, mpol
    do n = m, lrad, 2
      zernike(n,m) = zernike(n,m) / real(n+1)
    end do
  end do
end subroutine get_zernike_rm
