!> \file
!> \brief Measures error in Beltrami equation, \f$\nabla \times {\bf B} - \mu {\bf B}\f$.

!> \brief Measures error in Beltrami equation, \f$\nabla \times {\bf B} - \mu {\bf B}\f$.
!> \ingroup grp_diagnostics
!>
!> This routine is called by xspech() as a post diagnostic and only if \c Lcheck==1.
!>
!> **construction of current,** \f${\bf j} \equiv \nabla \times \nabla \times {\bf A}\f$
!> <ul>
!> <li> The components of the vector potential, \f${\bf A}=A_\theta \nabla + A_\zeta \nabla \zeta\f$, are
!>      \f{eqnarray}{
!>        A_\theta(s,\theta,\zeta) &=& \sum_{i,l} {\color{red}  A_{\theta,e,i,l}} \; {\overline T}_{l,i}(s) \cos\alpha_i + \sum_{i,l} {\color{Orange}  A_{\theta,o,i,l}} \; {\overline T}_{l,i}(s) \sin\alpha_i, \label{eq:At_jo00aa} \\
!>        A_\zeta (s,\theta,\zeta) &=& \sum_{i,l} {\color{blue} A_{\zeta, e,i,l}} \; {\overline T}_{l,i}(s) \cos\alpha_i + \sum_{i,l} {\color{Cerulean}A_{\zeta ,o,i,l}} \; {\overline T}_{l,i}(s) \sin\alpha_i, \label{eq:Az_jo00aa}
!>      \f}
!>      where \f${\overline T}_{l,i}(s) \equiv \bar s^{m_i/2} \, T_l(s)\f$, \f$T_l(s)\f$ is the Chebyshev polynomial, and \f$\alpha_j \equiv m_j\theta-n_j\zeta\f$.
!>      The regularity factor, \f$\bar s^{m_i/2}\f$, where \f$\bar s \equiv (1+s)/2\f$, is only included if there is a coordinate singularity in the domain
!>      (i.e. only in the innermost volume, and only in cylindrical and toroidal geometry.) </li>
!> <li> The magnetic field, \f$\sqrt g \, {\bf B} = \sqrt g B^s {\bf e}_s + \sqrt g B^\theta {\bf e}_\theta + \sqrt g B^\zeta {\bf e}_\zeta\f$, is
!>      \f{eqnarray}{
!>        \begin{array}{ccccrcrcrcrcccccccccccccccccccccccccccccccccccccccccccccccccccc}
!>        \sqrt g \, {\bf B} & = & {\bf e}_s      & \sum_{i,l} [ ( & - m_i {\color{blue} A_{\zeta, e,i,l}} & - & n_i {\color{red}  A_{\theta,e,i,l}} & ) {\overline T}_{l,i}        \sin\alpha_i + ( & + m_i {\color{Cerulean}A_{\zeta ,o,i,l}} & + & n_i {\color{Orange}  A_{\theta,o,i,l}} & ) {\overline T}_{l,i}        \cos\alpha_i ] \\
!>                           & + & {\bf e}_\theta & \sum_{i,l} [ ( &                                       & - &     {\color{blue} A_{\zeta, e,i,l}} & ) {\overline T}_{l,i}^\prime \cos\alpha_i + ( &                                          & - &     {\color{Cerulean}A_{\zeta ,o,i,l}} & ) {\overline T}_{l,i}^\prime \sin\alpha_i ] \\
!>                           & + & {\bf e}_\zeta  & \sum_{i,l} [ ( &       {\color{red}  A_{\theta,e,i,l}} &   &                                     & ) {\overline T}_{l,i}^\prime \cos\alpha_i + ( &       {\color{Orange}  A_{\theta,o,i,l}} &   &                                        & ) {\overline T}_{l,i}^\prime \sin\alpha_i ]
!>        \end{array}
!>      \f} </li>
!> <li> The current is
!>      \f{eqnarray}{ \sqrt g \, {\bf j} = ( \partial_\theta B_\zeta  - \partial_\zeta  B_\theta) \; {\bf e}_s
!>                                       + ( \partial_\zeta  B_s      - \partial_s      B_\zeta ) \; {\bf e}_\theta
!>                                       + ( \partial_s      B_\theta - \partial_\theta B_s     ) \; {\bf e}_\zeta ,
!>      \f}
!>      where (for computational convenience) the covariant components of \f${\bf B}\f$ are computed as
!>      \f{eqnarray}{
!>        B_s      & = & (\sqrt g B^s) \, g_{ss     } / \sqrt g + (\sqrt g B^\theta) \, g_{s     \theta} / \sqrt g + (\sqrt g B^\zeta) \, g_{s     \zeta} / \sqrt g, \\
!>        B_\theta & = & (\sqrt g B^s) \, g_{s\theta} / \sqrt g + (\sqrt g B^\theta) \, g_{\theta\theta} / \sqrt g + (\sqrt g B^\zeta) \, g_{\theta\zeta} / \sqrt g, \\
!>        B_\zeta  & = & (\sqrt g B^s) \, g_{s\zeta } / \sqrt g + (\sqrt g B^\theta) \, g_{\theta\zeta } / \sqrt g + (\sqrt g B^\zeta) \, g_{\zeta \zeta} / \sqrt g.
!>      \f} </li>
!> </ul>
!>
!> **quantification of the error**
!> <ul>
!> <li> The measures of the error are
!>      \f{eqnarray}{
!>          ||\left( {\bf j}-\mu {\bf B}\right)\cdot\nabla s      || &  \equiv  &
!>        \int \!\! ds \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \left| \sqrt g \, {\bf j}\cdot\nabla s     -\mu \; \sqrt g \, {\bf B} \cdot\nabla s      \right|, \label{eq:Es_jo00aa} \\
!>          ||\left( {\bf j}-\mu {\bf B}\right)\cdot\nabla \theta || &  \equiv  &
!>        \int \!\! ds \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \left| \sqrt g \, {\bf j}\cdot\nabla \theta-\mu \; \sqrt g \, {\bf B} \cdot\nabla \theta \right|, \label{eq:Et_jo00aa} \\
!>          ||\left( {\bf j}-\mu {\bf B}\right)\cdot\nabla \zeta  || &  \equiv  &
!>        \int \!\! ds \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \left| \sqrt g \, {\bf j}\cdot\nabla \zeta -\mu \; \sqrt g \, {\bf B} \cdot\nabla \zeta  \right|. \label{eq:Ez_jo00aa}
!>      \f} </li>
!> </ul>
!>
!> **comments**
!> <ul>
!> <li>  Is there a better definition and quantification of the error? For example, should we employ an error measure that is dimensionless? </li>
!> <li>  If the coordinate singularity is in the domain, then \f$|\nabla\theta|\rightarrow\infty\f$ at the coordinate origin.
!>       What then happens to \f$||\left( {\bf j}-\mu {\bf B}\right)\cdot\nabla \theta ||\f$ as defined in Eqn.\f$(\ref{eq:Et_jo00aa})\f$? </li>
!> <li>  What is the predicted scaling of the error in the Chebyshev-Fourier representation scale with numerical resolution?
!>       Note that the predicted error scaling for \f$E^s\f$, \f$E^\theta\f$ and \f$E^\zeta\f$ may not be standard,
!>       as various radial derivatives are taken to compute the components of \f${\bf j}\f$.
!>       (See for example the discussion in Sec.IV.C in Hudson et al. (2011) \cite y2011_hudson ,
!>       where the expected scaling of the error for a finite-element implementation is confirmed numerically.) </li>
!> <li>  Instead of using Gaussian integration to compute the integral over \f$s\f$, an adaptive quadrature algorithm may be preferable. </li>
!> </ul>
!>
!> @param[in] lvol  in which volume should the Beltrami error be computed
!> @param[in] Ntz   number of grid points in \f$\theta\f$ and \f$\zeta\f$
!> @param[in] lquad degree of Gaussian quadrature
!> @param[in] mn    number of Fourier harmonics
subroutine jo00aa(lvol, Ntz, lquad, mn)
    use mod_kinds, only: wp => dp
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    use constants, only: zero, half, one, two, pi2

    use fileunits, only: ounit

    use inputlist, only: Wmacros, Wjo00aa, Nvol, Lrad, mu, mpol, Igeometry, Nfp, Lerrortype

    use cputiming, only: Tjo00aa

    use allglobal, only: myid, cpus, MPI_COMM_SPEC, ext, ivol, &
                         im, in, regumm, &
                         Mvol, &
                         cheby, zernike, &
                         Ate, Aze, Ato, Azo, &
                         sg, guvij, Rij, Zij, &
                         Nt, Nz, efmn, ofmn, cfmn, sfmn, &
                         NOTstellsym, &
                         Lcoordinatesingularity, &
                         beltramierror, Rij, Zij, gBzeta, Node, &
                         pi2nfp, ivol, RTT, TT, dtflux, dpflux

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef OPENMP
    USE OMP_LIB
#endif
    use mpi
    implicit none
    integer :: ierr, astat, ios, nthreads, ithread
    real(wp) :: cput, cpui, cpuo = 0 ! cpu time; cpu initial; cpu old; 31 Jan 13;

!                        these are really global, but are included in argument list to remove allocations
    integer, intent(in) :: lvol, Ntz, lquad, mn

    integer :: jquad, Lcurvature, ll, ii, jj, kk, uu, ideriv, twolquad, mm, jk

    real(wp) :: lss, sbar, sbarhim(0:2), gBu(1:Ntz, 1:3, 0:3), gJu(1:Ntz, 1:3), jerror(1:3), jerrormax(1:3), intvol
    real(wp) :: B_cartesian(1:Ntz, 1:3), J_cartesian(1:Ntz, 1:3)

    real(wp) :: Atemn(1:mn, 0:2), Azemn(1:mn, 0:2), Atomn(1:mn, 0:2), Azomn(1:mn, 0:2)

    integer :: itype, icdgqf
    real(wp) :: aa, bb, cc, dd, weight(1:lquad + 1), abscis(1:lquad), workfield(1:2*lquad)

    real(wp) :: zeta, teta, st(2), Bst(2)

    cpui = MPI_WTIME()
    cpuo = cpui
#ifdef OPENMP
    nthreads = omp_get_max_threads()
#else
    nthreads = 1
#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef DEBUG

    if (lquad .lt. 1) then
        write (6, '("jo00aa :      fatal : myid=",i3," ; lquad.lt.1 ; invalid lquad supplied to jo00aa ;")') myid
        call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
        stop "jo00aa : lquad.lt.1 : invalid lquad supplied to jo00aa  ;"
    end if

    if (lvol .lt. 1 .or. lvol .gt. Mvol) then
        write (6, '("jo00aa :      fatal : myid=",i3," ; lvol.lt.1 .or. lvol.gt.Mvol ; invalid interface label ;")') myid
        call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
        stop "jo00aa : lvol.lt.1 .or. lvol.gt.Mvol : invalid interface label  ;"
    end if

#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **details of the numerics**
!> <ul>
!> <li> The integration over \f$s\f$ is performed using Gaussian integration, e.g., \f$\displaystyle \int \!\! f(s) ds \approx \sum_k \omega_k f(s_k)\f$;
!>      with the abscissae, \f$s_k\f$, and the weights, \f$\omega_k\f$, for \f$k=1,\f$ \c Iquad\f$_v\f$, determined by \c CDGQF.
!>      The resolution, \c N \f$ \equiv \f$ \c Iquad\f$_v\f$, is determined by \c Nquad  (see global.f90 and preset() ).
!>      A fatal error is enforced by jo00aa() if \c CDGQF returns an \c ifail \f$\ne 0\f$. </li>
    itype = 1; aa = -one; bb = +one; cc = zero; dd = zero; twolquad = 2*lquad

    call CDGQF(lquad, abscis(1:lquad), weight(1:lquad), itype, aa, bb, twolquad, workfield(1:twolquad), icdgqf) ! prepare Gaussian quadrature;
    weight(lquad + 1) = zero

! write(ounit,'("jo00aa :  WARNING ! : THE ERROR FLAGS RETURNED BY CDGQF SEEM DIFFERENT TO NAG:D01BCF (may be trivial, but please revise); 2018/01/10;")')

    cput = MPI_WTIME()
    select case (icdgqf) !                                                         123456789012345
    case (0); if (Wjo00aa) write (ounit, 1000) cput - cpus, myid, lvol, icdgqf, "success        ", abscis(1:lquad)
    case (1); write (ounit, 1000) cput - cpus, myid, lvol, icdgqf, "failed         ", abscis(1:lquad)
    case (2); write (ounit, 1000) cput - cpus, myid, lvol, icdgqf, "input error    ", abscis(1:lquad)
    case (3); write (ounit, 1000) cput - cpus, myid, lvol, icdgqf, "input error    ", abscis(1:lquad)
    case (4); write (ounit, 1000) cput - cpus, myid, lvol, icdgqf, "weight overflow", abscis(1:lquad)
    case (5); write (ounit, 1000) cput - cpus, myid, lvol, icdgqf, "weight zero    ", abscis(1:lquad)
    case (6); write (ounit, 1000) cput - cpus, myid, lvol, icdgqf, "failed         ", abscis(1:lquad)
    case default; write (ounit, 1000) cput - cpus, myid, lvol, icdgqf, "weird          ", abscis(1:lquad)
    end select

    if (Wjo00aa) write (ounit, 1001) weight(1:lquad)

1000 format("jo00aa : ", f10.2, " : myid=", i3, " ; lvol=", i3, " ; icdgqf=", i3, " ; "a15" ;":" abscissae ="99f10.6)
1001 format("jo00aa : ", 10x, " :       "3x"          "3x"            "3x"    "15x" ;":" weights   ="99f10.6)

    if (icdgqf .ne. 0) then
        write (6, '("jo00aa :      fatal : myid=",i3," ; icdgqf.ne.0 ; failed to construct Gaussian integration abscisae and weights ;")') myid
        call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
        stop "jo00aa : icdgqf.ne.0 : failed to construct Gaussian integration abscisae and weights  ;"
    end if

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    jerror(1:3) = zero; ideriv = 0 ! three components of the error in \curl B - mu B; initialize summation;
    jerrormax(1:3) = zero; intvol = zero

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> <li>  Inside the Gaussian quadrature loop, i.e. for each \f$s_k\f$,
!>       <ul>

    do jquad = 1, lquad + 1 ! loop over radial sub-sub-grid (numerical quadrature);

        if (jquad .eq. lquad + 1) then
            lss = one; sbar = one
        else
            lss = abscis(jquad); sbar = (lss + one)*half
        end if

!>       <li>  The metric elements, \f$g_{\mu,\nu} \equiv \f$ \c gij(1:6,0,1:Ntz), and the Jacobian, \f$\sqrt g \equiv \f$ \c sg(0,1:Ntz),
!>             are calculated on a regular angular grid, \f$(\theta_i,\zeta_j)\f$, in coords().
!>             The derivatives \f$\partial_i g_{\mu,\nu} \equiv\f$ \c gij(1:6,i,1:Ntz) and \f$\partial_i \sqrt g \equiv \f$ \c sg(i,1:Ntz),
!>             with respect to \f$ i \in \{ s,\theta,\zeta \}\f$ are also returned. </li>

        Lcurvature = 2

        cput = MPI_WTIME()
        Tjo00aa = Tjo00aa + (cput - cpuo)
        call coords(lvol, lss, Lcurvature, Ntz, mn)
        cpuo = MPI_WTIME()
        ! returns coordinates, metrics, . . .

        if (Lcoordinatesingularity) then ! Zernike 1 Jul 2019
            call get_zernike_d2(sbar, Lrad(lvol), mpol, zernike)
        else
            call get_cheby_d2(lss, Lrad(lvol), cheby(0:Lrad(lvol), 0:2))
        end if

        Atemn(1:mn, 0:2) = zero ! initialize summation over Chebyshev/Zernike polynomials;
        Azemn(1:mn, 0:2) = zero
        if (NOTstellsym) then
            Atomn(1:mn, 0:2) = zero
            Azomn(1:mn, 0:2) = zero
        else
            Atomn(1:mn, 0:2) = zero ! these are used below;
            Azomn(1:mn, 0:2) = zero
        end if

!>       <li>  The Fourier components of the vector potential given in Eqn.\f$(\ref{eq:At_jo00aa})\f$ and Eqn.\f$(\ref{eq:Az_jo00aa})\f$, and their first and second radial derivatives, are summed. </li>

        if (Lcoordinatesingularity) then
            do ll = 0, Lrad(lvol) ! radial (Chebyshev) resolution of magnetic vector potential;

                do ii = 1, mn  ! Fourier resolution of magnetic vector potential;
                    mm = im(ii)
                    if (ll < mm) cycle
                    if (mod(ll + mm, 2) > 0) cycle

                    ; Atemn(ii, 0) = Atemn(ii, 0) + Ate(lvol, ideriv, ii)%s(ll)*zernike(ll, mm, 0)
                    ; Atemn(ii, 1) = Atemn(ii, 1) + Ate(lvol, ideriv, ii)%s(ll)*zernike(ll, mm, 1)*half
                    ; Atemn(ii, 2) = Atemn(ii, 2) + Ate(lvol, ideriv, ii)%s(ll)*zernike(ll, mm, 2)*half*half

                    ; Azemn(ii, 0) = Azemn(ii, 0) + Aze(lvol, ideriv, ii)%s(ll)*zernike(ll, mm, 0)
                    ; Azemn(ii, 1) = Azemn(ii, 1) + Aze(lvol, ideriv, ii)%s(ll)*zernike(ll, mm, 1)*half
                    ; Azemn(ii, 2) = Azemn(ii, 2) + Aze(lvol, ideriv, ii)%s(ll)*zernike(ll, mm, 2)*half*half

                    if (NOTstellsym) then

                        Atomn(ii, 0) = Atomn(ii, 0) + Ato(lvol, ideriv, ii)%s(ll)*zernike(ll, mm, 0)
                        Atomn(ii, 1) = Atomn(ii, 1) + Ato(lvol, ideriv, ii)%s(ll)*zernike(ll, mm, 1)*half
                        Atomn(ii, 2) = Atomn(ii, 2) + Ato(lvol, ideriv, ii)%s(ll)*zernike(ll, mm, 2)*half*half

                        Azomn(ii, 0) = Azomn(ii, 0) + Azo(lvol, ideriv, ii)%s(ll)*zernike(ll, mm, 0)
                        Azomn(ii, 1) = Azomn(ii, 1) + Azo(lvol, ideriv, ii)%s(ll)*zernike(ll, mm, 1)*half
                        Azomn(ii, 2) = Azomn(ii, 2) + Azo(lvol, ideriv, ii)%s(ll)*zernike(ll, mm, 2)*half*half

                    end if ! end of if( NOTstellsym) ; 20 Jun 14;

                end do ! end of do ii;

            end do ! end of do ll;

        else

            do ll = 0, Lrad(lvol) ! radial (Chebyshev) resolution of magnetic vector potential;

!>       <li>  The quantities \f$\sqrt g B^s\f$, \f$\sqrt g B^\theta\f$ and \f$\sqrt g B^\zeta\f$, and their first and second derivatives with respect to \f$(s,\theta,\zeta)\f$,
!>             are computed on the regular angular grid (using FFTs). </li>

                do ii = 1, mn  ! Fourier resolution of magnetic vector potential;

                    ; Atemn(ii, 0) = Atemn(ii, 0) + Ate(lvol, ideriv, ii)%s(ll)*cheby(ll, 0)
                    ; Atemn(ii, 1) = Atemn(ii, 1) + Ate(lvol, ideriv, ii)%s(ll)*cheby(ll, 1)
                    ; Atemn(ii, 2) = Atemn(ii, 2) + Ate(lvol, ideriv, ii)%s(ll)*cheby(ll, 2)

                    ; Azemn(ii, 0) = Azemn(ii, 0) + Aze(lvol, ideriv, ii)%s(ll)*cheby(ll, 0)
                    ; Azemn(ii, 1) = Azemn(ii, 1) + Aze(lvol, ideriv, ii)%s(ll)*cheby(ll, 1)
                    ; Azemn(ii, 2) = Azemn(ii, 2) + Aze(lvol, ideriv, ii)%s(ll)*cheby(ll, 2)

                    if (NOTstellsym) then

                        Atomn(ii, 0) = Atomn(ii, 0) + Ato(lvol, ideriv, ii)%s(ll)*cheby(ll, 0)
                        Atomn(ii, 1) = Atomn(ii, 1) + Ato(lvol, ideriv, ii)%s(ll)*cheby(ll, 1)
                        Atomn(ii, 2) = Atomn(ii, 2) + Ato(lvol, ideriv, ii)%s(ll)*cheby(ll, 2)

                        Azomn(ii, 0) = Azomn(ii, 0) + Azo(lvol, ideriv, ii)%s(ll)*cheby(ll, 0)
                        Azomn(ii, 1) = Azomn(ii, 1) + Azo(lvol, ideriv, ii)%s(ll)*cheby(ll, 1)
                        Azomn(ii, 2) = Azomn(ii, 2) + Azo(lvol, ideriv, ii)%s(ll)*cheby(ll, 2)

                    end if ! end of if( NOTstellsym) ; 20 Jun 14;

                end do ! end of do ii;

            end do ! end of do ll;
        end if

        ofmn(1:mn) = -im(1:mn)*Azemn(1:mn, 0) - in(1:mn)*Atemn(1:mn, 0)
        efmn(1:mn) = +im(1:mn)*Azomn(1:mn, 0) + in(1:mn)*Atomn(1:mn, 0)
        sfmn(1:mn) = -im(1:mn)*Azemn(1:mn, 1) - in(1:mn)*Atemn(1:mn, 1)
        cfmn(1:mn) = +im(1:mn)*Azomn(1:mn, 1) + in(1:mn)*Atomn(1:mn, 1)

        call invfft(mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, gBu(1:Ntz, 1, 0), gBu(1:Ntz, 1, 1)) !  (gB^s)   , d(gB^s)/ds;

        efmn(1:mn) = -im(1:mn)*im(1:mn)*Azemn(1:mn, 0) - im(1:mn)*in(1:mn)*Atemn(1:mn, 0)
        ofmn(1:mn) = -im(1:mn)*im(1:mn)*Azomn(1:mn, 0) - im(1:mn)*in(1:mn)*Atomn(1:mn, 0)
        cfmn(1:mn) = +in(1:mn)*im(1:mn)*Azemn(1:mn, 0) + in(1:mn)*in(1:mn)*Atemn(1:mn, 0)
        sfmn(1:mn) = +in(1:mn)*im(1:mn)*Azomn(1:mn, 0) + in(1:mn)*in(1:mn)*Atomn(1:mn, 0)

        call invfft(mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, gBu(1:Ntz, 1, 2), gBu(1:Ntz, 1, 3)) ! d(gB^s)/dt, d(gB^s)/dz;

        efmn(1:mn) = -Azemn(1:mn, 1)
        ofmn(1:mn) = -Azomn(1:mn, 1)
        cfmn(1:mn) = -Azemn(1:mn, 2)
        sfmn(1:mn) = -Azomn(1:mn, 2)

        call invfft(mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, gBu(1:Ntz, 2, 0), gBu(1:Ntz, 2, 1)) !  (gB^t)   , d(gB^t)/ds;

        ofmn(1:mn) = +im(1:mn)*Azemn(1:mn, 1)
        efmn(1:mn) = -im(1:mn)*Azomn(1:mn, 1)
        sfmn(1:mn) = -in(1:mn)*Azemn(1:mn, 1)
        cfmn(1:mn) = +in(1:mn)*Azomn(1:mn, 1)

        call invfft(mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, gBu(1:Ntz, 2, 2), gBu(1:Ntz, 2, 3)) ! d(gB^t)/dt, d(gB^t)/dz;

        efmn(1:mn) = +Atemn(1:mn, 1)
        ofmn(1:mn) = +Atomn(1:mn, 1)
        cfmn(1:mn) = +Atemn(1:mn, 2)
        sfmn(1:mn) = +Atomn(1:mn, 2)

        call invfft(mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, gBu(1:Ntz, 3, 0), gBu(1:Ntz, 3, 1)) !  (gB^z)   , d(gB^z)/ds;

        ofmn(1:mn) = -im(1:mn)*Atemn(1:mn, 1)
        efmn(1:mn) = +im(1:mn)*Atomn(1:mn, 1)
        sfmn(1:mn) = +in(1:mn)*Atemn(1:mn, 1)
        cfmn(1:mn) = -in(1:mn)*Atomn(1:mn, 1)

        call invfft(mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, gBu(1:Ntz, 3, 2), gBu(1:Ntz, 3, 3)) ! d(gB^z)/dt, d(gB^z)/dz;

!>       <li>  The following quantities are then computed on the regular angular grid
!>             \f{eqnarray}{ \sqrt g j^s      & = & \sum_u \left[
!>                           \partial_\theta (\sqrt g B^u) \; g_{u,\zeta } + (\sqrt g B^u) \; \partial_\theta g_{u,\zeta } - (\sqrt g B^u) g_{u,\zeta } \; \partial_\theta \sqrt g / \sqrt g
!>                                                    \right] / \sqrt g \nonumber \\
!>                                       & - & \sum_u \left[
!>                           \partial_\zeta  (\sqrt g B^u) \; g_{u,\theta} + (\sqrt g B^u) \; \partial_\zeta  g_{u,\theta} - (\sqrt g B^u) g_{u,\theta} \; \partial_\zeta  \sqrt g / \sqrt g
!>                                                    \right] / \sqrt g, \\
!>                           \sqrt g j^\theta & = & \sum_u \left[
!>                           \partial_\zeta  (\sqrt g B^u) \; g_{u,s     } + (\sqrt g B^u) \; \partial_\zeta  g_{u,s     } - (\sqrt g B^u) g_{u,s     } \; \partial_\zeta  \sqrt g / \sqrt g
!>                                                    \right] / \sqrt g \nonumber \\
!>                                       & - & \sum_u \left[
!>                           \partial_s      (\sqrt g B^u) \; g_{u,\zeta } + (\sqrt g B^u) \; \partial_s      g_{u,\zeta } - (\sqrt g B^u) g_{u,\zeta } \; \partial_s      \sqrt g / \sqrt g
!>                                                    \right] / \sqrt g, \\
!>                           \sqrt g j^\zeta  & = & \sum_u \left[
!>                           \partial_s      (\sqrt g B^u) \; g_{u,\theta} + (\sqrt g B^u) \; \partial_s      g_{u,\theta} - (\sqrt g B^u) g_{u,\theta} \; \partial_s      \sqrt g / \sqrt g
!>                                                    \right] / \sqrt g \nonumber \\
!>                                       & - & \sum_u \left[
!>                           \partial_\theta (\sqrt g B^u) \; g_{u,s     } + (\sqrt g B^u) \; \partial_\theta g_{u,s     } - (\sqrt g B^u) g_{u,s     } \; \partial_\theta \sqrt g / \sqrt g
!>                                                    \right] / \sqrt g.
!>             \f} </li>
!>       </ul>

        do ii = 1, 3

            select case (ii)
            case (1); jj = 2; kk = 3
            case (2); jj = 3; kk = 1
            case (3); jj = 1; kk = 2
            end select

            gJu(1:Ntz, ii) = zero

            do uu = 1, 3 ! summation over uu;

                gJu(1:Ntz, ii) = gJu(1:Ntz, ii) &
                                 + (gBu(1:Ntz, uu, jj)*guvij(1:Ntz, uu, kk, 0) &
                                    + gBu(1:Ntz, uu, 0)*guvij(1:Ntz, uu, kk, jj) &
                                    - gBu(1:Ntz, uu, 0)*guvij(1:Ntz, uu, kk, 0)*sg(1:Ntz, jj)/sg(1:Ntz, 0) &
                                    ) &
                                 - (gBu(1:Ntz, uu, kk)*guvij(1:Ntz, uu, jj, 0) &
                                    + gBu(1:Ntz, uu, 0)*guvij(1:Ntz, uu, jj, kk) &
                                    - gBu(1:Ntz, uu, 0)*guvij(1:Ntz, uu, jj, 0)*sg(1:Ntz, kk)/sg(1:Ntz, 0) &
                                    )
            end do ! end of do uu;

            gJu(1:Ntz, ii) = gJu(1:Ntz, ii)/sg(1:Ntz, 0)

        end do ! end of do ii;

        select case (Igeometry)
        case (1)
        case (2)
        case (3)
            B_cartesian(:, 1) = (gBu(1:Ntz, 1, 0)*Rij(1:Ntz, 1, 0) + gBu(1:Ntz, 2, 0)*Rij(1:Ntz, 2, 0) + gBu(1:Ntz, 3, 0)*Rij(1:Ntz, 3, 0))/sg(1:Ntz, 0)
            B_cartesian(:, 2) = (gBu(1:Ntz, 1, 0)*Zij(1:Ntz, 1, 0) + gBu(1:Ntz, 2, 0)*Zij(1:Ntz, 2, 0) + gBu(1:Ntz, 3, 0)*Zij(1:Ntz, 3, 0))/sg(1:Ntz, 0)
            B_cartesian(:, 3) = (gBu(1:Ntz, 3, 0)*Rij(1:Ntz, 0, 0))/sg(1:Ntz, 0)

            J_cartesian(:, 1) = (gJu(1:Ntz, 1)*Rij(1:Ntz, 1, 0) + gJu(1:Ntz, 2)*Rij(1:Ntz, 2, 0) + gJu(1:Ntz, 3)*Rij(1:Ntz, 3, 0))/sg(1:Ntz, 0)
            J_cartesian(:, 2) = (gJu(1:Ntz, 1)*Zij(1:Ntz, 1, 0) + gJu(1:Ntz, 2)*Zij(1:Ntz, 2, 0) + gJu(1:Ntz, 3)*Zij(1:Ntz, 3, 0))/sg(1:Ntz, 0)
            J_cartesian(:, 3) = (gJu(1:Ntz, 3)*Rij(1:Ntz, 0, 0))/sg(1:Ntz, 0)
        end select

        if (Lerrortype .eq. 1 .and. Igeometry .eq. 3) then
            do ii = 1, 3; jerror(ii) = jerror(ii) + weight(jquad)*sum(sg(1:Ntz, 0)*abs(J_cartesian(1:Ntz, ii) - mu(lvol)*B_cartesian(1:Ntz, ii)))
                !if (maxval(abs(  J_cartesian(1:Ntz,ii) - mu(lvol) * B_cartesian(1:Ntz,ii)  )) > jerrormax(ii)) write(ounit,*) ii, lss,maxval(abs(  J_cartesian(1:Ntz,ii) - mu(lvol) * B_cartesian(1:Ntz,ii)  ))
                ; ; jerrormax(ii) = max(jerrormax(ii), maxval(abs(J_cartesian(1:Ntz, ii) - mu(lvol)*B_cartesian(1:Ntz, ii))))
            end do
        else
            do ii = 1, 3; jerror(ii) = jerror(ii) + weight(jquad)*sum(abs(gJu(1:Ntz, ii) - mu(lvol)*gBu(1:Ntz, ii, 0)))
                ; ; jerrormax(ii) = max(jerrormax(ii), maxval(abs(gJu(1:Ntz, ii) - mu(lvol)*gBu(1:Ntz, ii, 0))/sg(1:Ntz, 0)))
            end do
        end if
        intvol = intvol + weight(jquad)*sum(sg(1:Ntz, 0))

    end do ! end of do jquad;

    beltramierror(lvol, 1:3) = jerror(1:3)/Ntz  ! the 'tranditional' SPEC error
    jerror(1:3) = jerror(1:3)/intvol
    beltramierror(lvol, 4:6) = jerror(1:3)        ! the volume average error
    beltramierror(lvol, 7:9) = jerrormax(1:3)     ! the max error

    if (Lerrortype .eq. 1 .and. Igeometry .eq. 3) then
        cput = MPI_WTIME(); write (ounit, 1002) cput - cpus, myid, lvol, Lrad(lvol), jerror(1:3), cput - cpui ! write error to screen;
        ; ; write (ounit, 1003) cput - cpus, myid, lvol, Lrad(lvol), jerrormax(1:3), cput - cpui ! write error to screen;
    else
        cput = MPI_WTIME(); write (ounit, 1004) cput - cpus, myid, lvol, Lrad(lvol), jerror(1:3), cput - cpui ! write error to screen;
        ; ; write (ounit, 1005) cput - cpus, myid, lvol, Lrad(lvol), jerrormax(1:3), cput - cpui ! write error to screen;
    end if

!> <li>  The error is stored into an array called \c beltramierror which is then written to the HDF5 file in hdfint(). </li>

1002 format("jo00aa : ", f10.2, " : myid=", i3, " ; lvol =", i3, " ; lrad =", i3, " ; AVG E^\R="es23.15" , E^\Z="es23.15" , E^\phi="es23.15" ; time="f8.2"s ;")
1003 format("jo00aa : ", f10.2, " : myid=", i3, " ; lvol =", i3, " ; lrad =", i3, " ; MAX E^\R="es23.15" , E^\Z="es23.15" , E^\phi="es23.15" ; time="f8.2"s ;")
1004 format("jo00aa : ", f10.2, " : myid=", i3, " ; lvol =", i3, " ; lrad =", i3, " ; AVG E^\s="es23.15" , E^\t="es23.15" , E^\z="es23.15" ; time="f8.2"s ;")
1005 format("jo00aa : ", f10.2, " : myid=", i3, " ; lvol =", i3, " ; lrad =", i3, " ; MAX E^\s="es23.15" , E^\t="es23.15" , E^\z="es23.15" ; time="f8.2"s ;")

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! check boundary condition
    jerror = zero
    ivol = lvol

    do kk = 0, Nz - 1; zeta = kk*pi2nfp/Nz
        do jj = 0, Nt - 1; teta = jj*pi2/Nt; jk = 1 + jj + kk*Nt; st(1:2) = (/one, teta/)

            cput = MPI_WTIME()
            Tjo00aa = Tjo00aa + (cput - cpuo)
            call bfield(zeta, st(1:Node), Bst(1:Node))
            cpuo = MPI_WTIME()

            jerror(2) = max(jerror(2), abs(Bst(1)*gBzeta))

        end do
    end do

    if (.not. Lcoordinatesingularity) then
        do kk = 0, Nz - 1; zeta = kk*pi2nfp/Nz
            do jj = 0, Nt - 1; teta = jj*pi2/Nt; jk = 1 + jj + kk*Nt; st(1:2) = (/-one, teta/)

                cput = MPI_WTIME()
                Tjo00aa = Tjo00aa + (cput - cpuo)
                call bfield(zeta, st(1:Node), Bst(1:Node))
                cpuo = MPI_WTIME()

                jerror(1) = max(jerror(1), abs(Bst(1)*gBzeta))

            end do
        end do
    end if
    cput = MPI_WTIME(); write (ounit, 1006) cput - cpus, myid, lvol, Lrad(lvol), jerror(1:2), cput - cpui ! write error to screen;

    ! check fluxes
    Bst = zero

    if (Lcoordinatesingularity) then
        do ll = 0, Lrad(lvol)
            Bst(1) = Bst(1) + Ate(lvol, 0, 1)%s(ll)*RTT(ll, 0, 1, 0)
        end do
        Bst(1) = abs(Bst(1) - dtflux(lvol))
    else
        do ll = 0, Lrad(lvol)
            Bst(1) = Bst(1) + Ate(lvol, 0, 1)%s(ll)*TT(ll, 1, 0)
        end do
        Bst(1) = abs(Bst(1) - dtflux(lvol))
        do ll = 0, Lrad(lvol)
            Bst(2) = Bst(2) - Aze(lvol, 0, 1)%s(ll)*TT(ll, 1, 0)
        end do
        Bst(2) = abs(Bst(2) - dpflux(lvol))
    end if

    cput = MPI_WTIME(); write (ounit, 1007) cput - cpus, myid, lvol, Lrad(lvol), Bst(1:2), cput - cpui ! write error to screen;

1006 format("jo00aa : ", f10.2, " : myid=", i3, " ; lvol =", i3, " ; lrad =", i3, " ; MAX gB^s(-1)="es23.15" , gB^s(+1) ="es23.15" ; time="f8.2"s ;")
1007 format("jo00aa : ", f10.2, " : myid=", i3, " ; lvol =", i3, " ; lrad =", i3, " ; dtfluxERR   ="es23.15" , dpfluxERR="es23.15" ; time="f8.2"s ;")

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

9999 continue
    cput = MPI_WTIME()
    Tjo00aa = Tjo00aa + (cput - cpuo)
    return

!> </ul>

end subroutine jo00aa
