module zernik

contains

  subroutine get_zernike(r, lrad, mpol, zernike)
    ! Get the Zernike polynomials with zeroth, first derivatives
    ! Inputs:
    ! r - REAL, coordinate input r
    ! lrad - INTEGER, radial resolution
    ! mpol - INTEGER, poloidal resolution
    !
    ! Returns:
    ! zernike - REAL(0:lrad,0:mpol,0:2), the value, first derivative of Zernike polynomial

    use constants, only : zero, one, two

    implicit none

    REAL,intent(in) :: r
    INTEGER, intent(in) :: lrad, mpol
    REAL, intent(inout) :: zernike(0:lrad,0:mpol,0:2)

    REAL ::    rm  ! r to the power of m'th
    REAL ::    factor1, factor2, factor3, factor4
    INTEGER :: m, n  ! Zernike R^m_n
    
    rm = 1 ! r to the power of m'th
    zernike(:,:,:) = zero
    do m = 0, mpol
      if (lrad >= m) then
        zernike(m,m,0:1) = (/ rm, real(m)*rm/r /)
      endif

      if (lrad >= m+2) then
        zernike(m+2,m,0) = real(m+2)*rm*r**2 - real(m+1)*rm
        zernike(m+2,m,1) = real((m+2)**2)*rm*r - real((m+1)*m)*rm/r
      endif

      do n = m+4, lrad, 2
        factor1 = real(n)/real(n**2 - m**2)
        factor2 = real(4 * (n-1))
        factor3 = real((n-2+m)**2)/real(n-2) + real((n-m)**2)/real(n)
        factor4 = real((n-2)**2-m**2) / real(n-2)

        zernike(n, m, 0) = factor1 * ((factor2*r**2 - factor3)*zernike(n-2,m,0) - factor4*zernike(n-4,m,0))
        zernike(n, m, 1) = factor1 * (two*factor2*r*zernike(n-2,m,0) + (factor2*r**2 - factor3)*zernike(n-2,m,1) - factor4*zernike(n-4,m,1))
      enddo

      rm = rm * r

    enddo
  end subroutine

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

    REAL ::    rm  ! r to the power of m'th
    REAL ::    factor1, factor2, factor3, factor4
    INTEGER :: m, n  ! Zernike R^m_n
    
    rm = 1 ! r to the power of m'th
    zernike(:,:,:) = zero
    do m = 0, mpol
      if (lrad >= m) then
        zernike(m,m,0:2) = (/ rm, real(m)*rm/r, real(m*(m-1))*rm/r/r /)
      endif

      if (lrad >= m+2) then
        zernike(m+2,m,0) = real(m+2)*rm*r**2 - real(m+1)*rm
        zernike(m+2,m,1) = real((m+2)**2)*rm*r - real((m+1)*m)*rm/r 
        zernike(m+2,m,2) = real((m+2)**2*(m+1))*rm - real((m+1)*m*(m-1))*rm/r/r 
      endif

      do n = m+4, lrad, 2
        factor1 = real(n)/real(n**2 - m**2)
        factor2 = real(4 * (n-1))
        factor3 = real((n-2+m)**2)/real(n-2) + real((n-m)**2)/real(n)
        factor4 = real((n-2)**2-m**2) / real(n-2)

        zernike(n, m, 0) = factor1 * ((factor2*r**2 - factor3)*zernike(n-2,m,0) - factor4*zernike(n-4,m,0))
        zernike(n, m, 1) = factor1 * (two*factor2*r*zernike(n-2,m,0) + (factor2*r**2 - factor3)*zernike(n-2,m,1) - factor4*zernike(n-4,m,1))
        zernike(n, m, 2) = factor1 * (two*factor2*(r*zernike(n-2,m,1)+zernike(n-2,m,0)) + (factor2*r**2 - factor3)*zernike(n-2,m,2) - factor4*zernike(n-4,m,2))
      enddo

      rm = rm * r

    enddo

  end subroutine

end module zernik