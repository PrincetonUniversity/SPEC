!> \defgroup grp_numerics Some miscellaneous numerical routines

!> \file numrec.f90
!> \brief Various miscellaneous "numerical" routines.

!> \brief Assign Fourier mode labels
!> \ingroup grp_numerics
!>
!> <ul>
!>
!> <li> This routine assigns the Fourier mode labels that converts a double-sum into a single sum; i.e., the \f$m_j\f$ and \f$n_j\f$ are assigned where
!> \f{eqnarray}{ f(\theta,\zeta) & = & \sum_{n=0}^{N} f_{0,n}\cos(-n \, N_P \, \zeta) + \sum_{m=1}^{M} \sum_{n=-N}^{N} f_{m,n}\cos(m\theta-n \, N_P \, \zeta) \\
!>                               & = & \sum_j f_j \cos(m_j\theta-n_j\zeta), \label{eq:condensedFourierrepresentation_numrec}
!> \f}
!> where \f$N\equiv\,\f$\c Ntor and \f$M\equiv\,\f$\c Mpol are given on input, and \f$N_P \equiv\,\f$\c Nfp is the field periodicity. </li>
!> </ul>
subroutine gi00ab( Mpol, Ntor, Nfp, mn, im, in )
  
  implicit none
  
  INTEGER, intent(in)  :: Mpol, Ntor, Nfp, mn
  INTEGER, intent(out) :: im(mn), in(mn)
  
  INTEGER              :: imn, mm, nn
  
  imn = 0
  
  ;  mm = 0  
  ;do nn = 0, Ntor
  ; imn = imn+1 ; im(imn) = mm ; in(imn) = nn*Nfp
  ;enddo
  ;

  do mm = 1, Mpol
   do nn = -Ntor, Ntor
    imn = imn+1 ; im(imn) = mm ; in(imn) = nn*Nfp
   enddo
  enddo
  
  return
  
end subroutine gi00ab


!> \brief Forward Fourier transform (fftw wrapper)
!> \ingroup grp_numerics
!> 
!> <ul>
!> <li> This constructs the "forward" Fourier transform. </li>
!> <li> Given a set of data, \f$(f_{i},g_{i})\f$ for \f$i = 1, \dots N_\theta N_\zeta\f$, on a regular two-dimensional angle grid, 
!>      where \f$\theta_j = 2 \pi j / N_\theta\f$ for \f$j = 0, N_\theta-1\f$, and
!>            \f$\zeta_k  = 2 \pi k / N_\zeta \f$ for \f$k = 0, N_\zeta -1\f$.
!>      The "packing" is governed by \f$i = 1 + j + k N_\theta\f$.
!>      The "discrete" resolution is \f$N_\theta \equiv\,\f$\c Nt, \f$N_\zeta \equiv\,\f$\c Nz and \c Ntz \f$=\f$ \c Nt \f$\times\f$ \c Nz , 
!>      which are set in preset() . </li>
!> <li> The Fourier harmonics consistent with Eqn.\f$(\ref{eq:condensedFourierrepresentation_numrec})\f$ are constructed.
!>      The mode identification labels appearing in Eqn.\f$(\ref{eq:condensedFourierrepresentation_numrec})\f$ are \f$m_j \equiv\,\f$\c im(j) and \f$n_j \equiv\,\f$\c in(j) ,
!>      which are set in readin() via a call to gi00ab() . </li>
!> </ul>
!>
!> @param  Nt
!> @param  Nz
!> @param  ijreal
!> @param  ijimag
!> @param  mn
!> @param  im
!> @param  in
!> @param  efmn
!> @param  ofmn
!> @param  cfmn
!> @param  sfmn
!> @param  ifail 
subroutine tfft( Nt, Nz, ijreal, ijimag, mn, im, in, efmn, ofmn, cfmn, sfmn, ifail )

  use constants, only : half, zero, pi2
  
  use fileunits, only : ounit

  use inputlist, only : Nfp
  use allglobal, only : pi2nfp

  use fftw_interface

  implicit none

  intrinsic aimag
  
  INTEGER :: Nt, Nz, mn, im(1:mn), in(1:mn), Ntz, imn, ifail, mm, nn
  REAL    :: ijreal(1:Nt*Nz), ijimag(1:Nt*Nz), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn)
  
  LOGICAL :: Lcheck = .false.
  INTEGER :: jj, kk
  REAL    :: jireal(1:Nt*Nz), jiimag(1:Nt*Nz), arg, ca, sa
  COMPLEX(C_DOUBLE_COMPLEX) :: z1, z2, z3

  if( Lcheck ) then ; jireal = ijreal ; jiimag = ijimag
  endif

  do jj = 1, Nz ; cplxin(:,jj) = CMPLX( ijreal((jj-1)*Nt+1:jj*Nt), ijimag((jj-1)*Nt+1:jj*Nt), KIND=C_DOUBLE_COMPLEX )
  enddo
   
  call fftw_execute_dft( planf, cplxin, cplxout ) ! Forward transform
  Ntz = Nt * Nz
  cplxout = cplxout / Ntz
  cplxout(1,1) = half*cplxout(1,1)

  do imn = 1, mn
   mm = im(imn);  nn = in(imn) / Nfp

   z1 = cplxout(1 + MOD(Nt - mm, Nt), 1 + MOD(Nz + nn, Nz))
   z2 = cplxout(1 +          mm,      1 + MOD(Nz - nn, Nz))

   z3 = z1 + z2
   efmn(imn) =  real(z3);  cfmn(imn) = aimag(z3)

   z3 = z1 - z2
   ofmn(imn) = aimag(z3);  sfmn(imn) = -real(z3)
  enddo

  if( .not.Lcheck ) return

  ijreal(1:Ntz) = zero ; ijimag(1:Ntz) = zero
  
  do jj = 0, Nt-1

   do kk = 0, Nz-1

    do imn = 1, mn ; arg = im(imn) * jj * pi2 / Nt - in(imn) * kk * pi2nfp / Nz ; ca = cos(arg) ; sa = sin(arg)
     
     ijreal(1+jj+kk*Nt) = ijreal(1+jj+kk*Nt) + efmn(imn) * ca + ofmn(imn) * sa
     ijimag(1+jj+kk*Nt) = ijimag(1+jj+kk*Nt) + cfmn(imn) * ca + sfmn(imn) * sa
     
    enddo
   enddo
  enddo
  
  write(ounit,'("tfft   : ",10x," : Fourier reconstruction error =",2es15.5," ;")') sqrt(sum((ijreal-jireal)**2)/Ntz), sqrt(sum((ijimag-jiimag)**2)/Ntz)

  return

end subroutine tfft

!> \brief Inverse Fourier transform (fftw wrapper)
!> \ingroup grp_numerics
!> 
!> <ul>
!> <li> Given the Fourier harmonics, the data on a regular angular grid are constructed. </li>
!> <li> This is the inverse routine to tfft() . </li>
!> </ul>
!>
!> @param[in]  mn
!> @param[in]  im
!> @param[in]  in
!> @param[in]  efmn
!> @param[in]  ofmn
!> @param[in]  cfmn
!> @param[in]  sfmn
!> @param[in]  Nt
!> @param[in]  Nz
!> @param[out] ijreal
!> @param[out] ijimag
subroutine invfft( mn, im, in, efmn, ofmn, cfmn, sfmn, Nt, Nz, ijreal, ijimag )

  use constants, only : zero, two, half
  use inputlist, only : Nfp
  use fftw_interface
  implicit none
  
  INTEGER, intent(in)  :: mn, im(mn), in(mn)
  REAL   , intent(in)  :: efmn(mn), ofmn(mn), cfmn(mn), sfmn(mn)
  INTEGER, intent(in)  :: Nt, Nz
  REAL   , intent(out) :: ijreal(Nt*Nz), ijimag(Nt*Nz) ! output real space;
  
  INTEGER              :: imn, jj, mm, nn

  cplxin = zero

  !Copy real arrays to complex
  do imn = 1,mn ; mm = im(imn) ; nn = in(imn) / Nfp
     cplxin(1 + MOD(Nt - mm, Nt), 1 + MOD(Nz + nn, Nz)) = &
          half * CMPLX(efmn(imn) - sfmn(imn), cfmn(imn) + ofmn(imn), KIND=C_DOUBLE_COMPLEX)
     cplxin(1 +          mm,      1 + MOD(Nz - nn, Nz)) = &
          half * CMPLX(efmn(imn) + sfmn(imn), cfmn(imn) - ofmn(imn), KIND=C_DOUBLE_COMPLEX)
  enddo
  cplxin(1,1) = two*cplxin(1,1)

  call fftw_execute_dft(planb, cplxin, cplxout) !Inverse transform

  !Copy complex result back to real arrays
  do jj=1,Nz
     ijreal((jj-1)*Nt+1:jj*Nt) =  real(cplxout(:,jj))
     ijimag((jj-1)*Nt+1:jj*Nt) = aimag(cplxout(:,jj))
  enddo

  return

end subroutine invfft


!> \brief Gauss-Legendre weights and abscissae
!> \ingroup grp_numerics
!> 
!> <ul>
!> <li> Compute Gaussian integration weights and abscissae. </li>
!> <li> From Numerical Recipes. </li>
!> </ul>
!>
!> @param[in]  n
!> @param[out] weight
!> @param[out] abscis
!> @param[out] ifail
subroutine gauleg( n, weight, abscis, ifail )
  
  use constants, only : zero, one, two, pi
  
  implicit none
  
  intrinsic abs, cos, epsilon
  
  INTEGER,            intent(in)  :: n
  REAL, dimension(n), intent(out) :: weight, abscis
  INTEGER,            intent(out) :: ifail

  INTEGER, parameter :: maxiter=16
  INTEGER            :: m, j, i, irefl, iter
  REAL               :: z1,z,pp,p3,p2,p1
  REAL, parameter    :: eps = epsilon(z)

  !Error checking
  if( n < 1 ) then ; ifail = 2 ;  return
  endif

  m = (n + 1)/2  !Roots are symmetric in interval, so we only need half
  do i=1,m       !Loop over desired roots
     irefl = n + 1 - i
     if (i .ne. irefl) then
        z = cos(pi*(i - 0.25)/(n + 0.5))  ! Approximate ith root
     else        !For an odd number of abscissae, the center must be at zero by symmetry.
        z = 0.0
     endif

     !Refine by Newton method
     do iter=1,maxiter
        p1 = one;  p2 = zero           ! Initialize recurrence relation

        do j=1,n  !Recurrence relation to get P(x)
           p3 = p2;  p2 = p1
           p1 = ((two*j - one)*z*p2 - (j - one)*p3)/j
        enddo !j

        pp = n*(z*p1 - p2)/(z*z - one) !Derivative of P(x)
        z1 = z;  z = z1 - p1/pp        !Newton iteration
        if (abs(z - z1) .le. eps) exit !Convergence test
     enddo !iter
     if (iter > maxiter) then
        ifail = 1;  return
     endif

     abscis(i) = -z;  abscis(irefl) = z
     weight(i) = two/((one - z*z)*pp*pp)
     weight(irefl) = weight(i)
  enddo !i

  ifail = 0
end subroutine gauleg
