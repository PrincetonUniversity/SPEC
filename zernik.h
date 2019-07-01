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
        zernike(m,m,0:1) = (/ r**m, real(m)*rm1 /)
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
      endif

      if (lrad >= m+2) then
        zernike(m+2,m,0) = real(m+2)*rm*r**2 - real(m+1)*rm
        zernike(m+2,m,1) = real((m+2)**2)*rm*r - real((m+1)*m)*rm1 
        zernike(m+2,m,2) = real((m+2)**2*(m+1))*rm - real((m+1)*m*(m-1))*rm2 
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
      enddo

      rm2 = rm1
      rm1 = rm
      rm = rm * r

    enddo

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

    REAL ::    rm, rm1  ! r to the power of m'th and m-1'th
    REAL ::    factor1, factor2, factor3, factor4
    INTEGER :: m, n  ! Zernike R^m_n
    
    rm = one  ! r to the power of m'th / r^m

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

      rm = rm * r
    enddo
  end subroutine get_zernike_rm

end module zernik