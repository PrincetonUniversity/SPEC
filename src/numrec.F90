!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (numerics) ! Some miscellaneous numerical routines.

!latex \briefly{miscellaneous ``numerical'' routines}

!l tex \calledby{\link{}}
!l tex \calls{\link{}}

!latex \tableofcontents

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{Outline}

!latex This file contains various miscellaneous ``numerical'' routines as described below.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!l tex \begin{itemize}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!l tex \item \type{gi00aa}

!subroutine gi00aa( ii, jj, ig ) ! not used; SRH: 27 Feb 18;
!
!  implicit none
!
!  INTEGER, intent(in)  :: ii,jj
!  INTEGER, intent(out) :: ig
!
!  if( ( ii.eq.1 .and. jj.eq.1 )                                ) ig = 1
!  if( ( ii.eq.1 .and. jj.eq.2 ) .or. ( ii.eq.2 .and. jj.eq.1 ) ) ig = 2
!  if( ( ii.eq.1 .and. jj.eq.3 ) .or. ( ii.eq.3 .and. jj.eq.1 ) ) ig = 3
!  if( ( ii.eq.2 .and. jj.eq.2 )                                ) ig = 4
!  if( ( ii.eq.2 .and. jj.eq.3 ) .or. ( ii.eq.3 .and. jj.eq.2 ) ) ig = 5
!  if( ( ii.eq.3 .and. jj.eq.3 )                                ) ig = 6
!
!  return
!
!end subroutine gi00aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{\type{gi00ab}}

!latex \begin{enumerate}

!latex \item This routine assigns the Fourier mode labels that converts a double-sum into a single sum; i.e., the $m_j$ and $n_j$ are assigned where
!latex \be f(\t,\z) & = & \sum_{n=0}^{N} f_{0,n}\cos(-n \, N_P \, \z)
!latex + \sum_{m=1}^{M} \sum_{n=-N}^{N} f_{m,n}\cos(m\t-n \, N_P \, \z) \\
!latex              & = & \sum_j f_j \cos(m_j\t-n_j\z), \label{eq:condensedFourierrepresentation}
!latex \ee
!latex where $N\equiv $ \type{Ntor} and $M\equiv $ \type{Mpol} are given on input, and $N_P \equiv $ \type{Nfp} is the field periodicity.

!latex \end{enumerate}

subroutine gi00ab(Mpol, Ntor, Nfp, mn, im, in)
    use mod_kinds, only: wp => dp
    implicit none

    integer, intent(in) :: Mpol, Ntor, Nfp, mn
    integer, intent(out) :: im(mn), in(mn)

    integer :: imn, mm, nn

    imn = 0

    ; mm = 0
    ; do nn = 0, Ntor
        ; imn = imn + 1; im(imn) = mm; in(imn) = nn*Nfp
        ; end do
    ; 
    do mm = 1, Mpol
        do nn = -Ntor, Ntor
            imn = imn + 1; im(imn) = mm; in(imn) = nn*Nfp
        end do
    end do

    return

end subroutine gi00ab

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine getimn(Mpol, Ntor, Nfp, mi, ni, idx)
    use mod_kinds, only: wp => dp
    ! convert m and n to index
    implicit none
    integer, intent(in) :: Mpol, Ntor, Nfp, mi, ni
    integer, intent(out) :: idx

    if (mi .gt. Mpol .or. mi .lt. 0 .or. ni .gt. Ntor*Nfp .or. ni .lt. -Ntor*Nfp) then
        idx = 0
    elseif (mi .eq. 0) then
        idx = 1 + ni/Nfp
    else
        idx = 1 + Ntor + (2*Ntor + 1)*(mi - 1) + (ni/Nfp + Ntor + 1)
    end if

end subroutine getimn

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{\type{tfft}}

!latex \begin{enumerate}

!latex \item This constructs the ``forward'' Fourier transform.

!latex \item Given a set of data, $(f_{i},g_{i})$ for $i = 1, \dots N_\theta N_\zeta$, on a regular two-dimensional angle grid,
!latex       where $\theta_j = 2 \pi j / N_\theta$ for $j = 0, N_\theta-1$, and
!latex             $\zeta_k  = 2 \pi k / N_\zeta $ for $k = 0, N_\zeta -1$.
!latex       The ``packing'' is governed by $i = 1 + j + k N_\theta$.
!latex       The ``discrete'' resolution is $N_\theta \equiv $ \type{Nt}, $N_\zeta \equiv $ \type{Nz} and \type{Ntz} $=$ \type{Nt} $\times$ \type{Nz},
!latex       which are set in \link{preset}.
!latex \item The Fourier harmonics consistent with \Eqn{condensedFourierrepresentation} are constructed.
!latex       The mode identification labels appearing in \Eqn{condensedFourierrepresentation} are $m_j \equiv $ \type{im(j)} and $n_j \equiv $ \type{in(j)},
!latex       which are set in \link{global} via a call to \type{gi00ab}.

!latex \end{enumerate}

subroutine tfft(Nt, Nz, ijreal, ijimag, mn, im, in, efmn, ofmn, cfmn, sfmn, ifail)
    use mod_kinds, only: wp => dp
    use constants, only: half, zero, pi2

    use fileunits, only: ounit

    use inputlist, only: Nfp
    use allglobal, only: pi2nfp

    use fftw_interface
#ifdef OPENMP
    use OMP_LIB
#endif
    implicit none

    intrinsic aimag

    integer :: Nt, Nz, mn, im(1:mn), in(1:mn), Ntz, imn, ifail, mm, nn
    real(wp) :: ijreal(1:Nt*Nz), ijimag(1:Nt*Nz), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn)

    LOGICAL :: Lcheck = .false.
    integer :: jj, kk, ithread
    !REAL    :: jireal(1:Nt*Nz), jiimag(1:Nt*Nz), arg, ca, sa
    real(wp) :: arg, ca, sa
    COMPLEX(C_DOUBLE_COMPLEX) :: z1, z2, z3

#ifdef OPENMP
    ithread = omp_get_thread_num() + 1
#else
    ithread = 1
#endif

    !if( Lcheck ) then ; jireal = ijreal ; jiimag = ijimag
    !endif

    do jj = 1, Nz; cplxin(:, jj, ithread) = CMPLX(ijreal((jj - 1)*Nt + 1:jj*Nt), ijimag((jj - 1)*Nt + 1:jj*Nt), KIND=C_DOUBLE_COMPLEX)
    end do

    call fftw_execute_dft(planf, cplxin(:, :, ithread), cplxout(:, :, ithread)) !Forward transform
    Ntz = Nt*Nz
    cplxout(:, :, ithread) = cplxout(:, :, ithread)/Ntz
    cplxout(1, 1, ithread) = half*cplxout(1, 1, ithread)

    do imn = 1, mn
        mm = im(imn); nn = in(imn)/Nfp

        z1 = cplxout(1 + MOD(Nt - mm, Nt), 1 + MOD(Nz + nn, Nz), ithread)
        z2 = cplxout(1 + mm, 1 + MOD(Nz - nn, Nz), ithread)

        z3 = z1 + z2
        efmn(imn) = real(z3); cfmn(imn) = aimag(z3)

        z3 = z1 - z2
        ofmn(imn) = aimag(z3); sfmn(imn) = -real(z3)
    end do

    if (.not. Lcheck) return

    ijreal(1:Ntz) = zero; ijimag(1:Ntz) = zero

    do jj = 0, Nt - 1

        do kk = 0, Nz - 1

            do imn = 1, mn; arg = im(imn)*jj*pi2/Nt - in(imn)*kk*pi2nfp/Nz; ca = cos(arg); sa = sin(arg)

                ijreal(1 + jj + kk*Nt) = ijreal(1 + jj + kk*Nt) + efmn(imn)*ca + ofmn(imn)*sa
                ijimag(1 + jj + kk*Nt) = ijimag(1 + jj + kk*Nt) + cfmn(imn)*ca + sfmn(imn)*sa

            end do
        end do
    end do

    !write(ounit,'("tfft   : ",10x," : Fourier reconstruction error =",2es15.5," ;")') sqrt(sum((ijreal-jireal)**2)/Ntz), sqrt(sum((ijimag-jiimag)**2)/Ntz)

    return

end subroutine tfft

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{\type{invfft}}

!latex \begin{enumerate}

!latex \item Given the Fourier harmonics, the data on a regular angular grid are constructed.

!latex \item This is the inverse routine to \type{tfft}.

!latex \end{enumerate}

subroutine invfft(mn, im, in, efmn, ofmn, cfmn, sfmn, Nt, Nz, ijreal, ijimag)
    use mod_kinds, only: wp => dp
    use constants, only: zero, two, half
    use inputlist, only: Nfp
    use fftw_interface
#ifdef OPENMP
    use OMP_LIB
#endif

    implicit none

    integer, intent(in) :: mn, im(mn), in(mn)
    real(wp), intent(in) :: efmn(mn), ofmn(mn), cfmn(mn), sfmn(mn)
    integer, intent(in) :: Nt, Nz
    real(wp), intent(out) :: ijreal(Nt*Nz), ijimag(Nt*Nz) ! output real space;

    integer :: imn, jj, mm, nn, ithread

#ifdef OPENMP
    ithread = omp_get_thread_num() + 1
#else
    ithread = 1
#endif

    cplxin(:, :, ithread) = zero

    !Copy real arrays to complex
    do imn = 1, mn; mm = im(imn); nn = in(imn)/Nfp
        cplxin(1 + MOD(Nt - mm, Nt), 1 + MOD(Nz + nn, Nz), ithread) = &
            half*CMPLX(efmn(imn) - sfmn(imn), cfmn(imn) + ofmn(imn), KIND=C_DOUBLE_COMPLEX)
        cplxin(1 + mm, 1 + MOD(Nz - nn, Nz), ithread) = &
            half*CMPLX(efmn(imn) + sfmn(imn), cfmn(imn) - ofmn(imn), KIND=C_DOUBLE_COMPLEX)
    end do
    cplxin(1, 1, ithread) = two*cplxin(1, 1, ithread)

    call fftw_execute_dft(planb, cplxin(:, :, ithread), cplxout(:, :, ithread)) !Inverse transform

    !Copy complex result back to real arrays
    do jj = 1, Nz
        ijreal((jj - 1)*Nt + 1:jj*Nt) = real(cplxout(:, jj, ithread))
        ijimag((jj - 1)*Nt + 1:jj*Nt) = aimag(cplxout(:, jj, ithread))
    end do

    return

end subroutine invfft

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{\type{gauleg}}

!latex \begin{enumerate}

!latex \item Compute Gaussian integration weights and abscissae.

!latex \item From Numerical Recipes.

!latex \end{enumerate}

subroutine gauleg(n, weight, abscis, ifail)
    use mod_kinds, only: wp => dp
    use constants, only: zero, one, two, pi

    implicit none

    intrinsic abs, cos, epsilon

    integer, intent(in) :: n
    real(wp), dimension(n), intent(out) :: weight, abscis
    integer, intent(out) :: ifail

    integer, parameter :: maxiter = 16
    integer :: m, j, i, irefl, iter
    real(wp) :: z1, z, pp, p3, p2, p1
    real(wp), parameter :: eps = epsilon(z)

    !Error checking
    if (n < 1) then; ifail = 2; return
    end if

    m = (n + 1)/2  !Roots are symmetric in interval, so we only need half
    do i = 1, m       !Loop over desired roots
        irefl = n + 1 - i
        if (i .ne. irefl) then
            z = cos(pi*(i - 0.25)/(n + 0.5))  ! Approximate ith root
        else        !For an odd number of abscissae, the center must be at zero by symmetry.
            z = 0.0
        end if

        !Refine by Newton method
        do iter = 1, maxiter
            p1 = one; p2 = zero           ! Initialize recurrence relation

            do j = 1, n  !Recurrence relation to get P(x)
                p3 = p2; p2 = p1
                p1 = ((two*j - one)*z*p2 - (j - one)*p3)/j
            end do !j

            pp = n*(z*p1 - p2)/(z*z - one) !Derivative of P(x)
            z1 = z; z = z1 - p1/pp        !Newton iteration
            if (abs(z - z1) .le. eps) exit !Convergence test
        end do !iter
        if (iter > maxiter) then
            ifail = 1; return
        end if

        abscis(i) = -z; abscis(irefl) = z
        weight(i) = two/((one - z*z)*pp*pp)
        weight(irefl) = weight(i)
    end do !i

    ifail = 0
end subroutine gauleg
