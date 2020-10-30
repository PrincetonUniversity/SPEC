subroutine get_cheby(lss, lrad, cheby)
  ! Get the Chebyshev polynomials with zeroth, first derivatives
  ! Inputs:
  ! lss - REAL, coordinate input lss
  ! lrad - INTEGER, radial resolution
  !
  ! Returns:
  ! cheby - REAL(0:Mrad,0:1), the value, first derivative of Chebyshev polynomial

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

subroutine get_cheby_d2(lss, lrad, cheby)
  ! Get the Chebyshev polynomials with zeroth, first and second derivatives
  ! Inputs:
  ! lss - REAL, coordinate input lss
  ! lrad - INTEGER, radial resolution
  !
  ! Returns:
  ! cheby - REAL(0:lrad,0:2), the value, first and second derivative of Chebyshev polynomial

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


subroutine get_zernike(r, lrad, mpol, zernike)
  ! Get the Zernike polynomials with zeroth, first derivatives
  ! Inputs:
  ! r - REAL, coordinate input r
  ! lrad - INTEGER, radial resolution
  ! mpol - INTEGER, poloidal resolution
  !
  ! Returns:
  ! zernike - REAL(0:lrad,0:mpol,0:1), the value, first derivative of Zernike polynomial

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

subroutine get_zernike_d2(r, lrad, mpol, zernike)
  ! Get the Zernike polynomials with zeroth, first, second derivatives
  ! Inputs:
  ! r - REAL, coordinate input r
  ! lrad - INTEGER, radial resolution
  ! mpol - INTEGER, poloidal resolution
  !
  ! Returns:
  ! zernike - REAL(0:lrad,0:mpol,0:2), the value, first/second derivative of Zernike polynomial
  
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

subroutine get_zernike_rm(r, lrad, mpol, zernike)
  ! Get the Zernike polynomials Z(n,m,r)/r^m
  ! Inputs:
  ! r - REAL, coordinate input r
  ! lrad - INTEGER, radial resolution
  ! mpol - INTEGER, poloidal resolution
  !
  ! Returns:
  ! zernike - REAL(0:lrad,0:mpol), the value

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