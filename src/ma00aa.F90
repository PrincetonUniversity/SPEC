!> \defgroup grp_integrals Integrals
!>
!> \file
!> \brief Calculates volume integrals of Chebyshev polynomials and metric element products.

!> \brief Calculates volume integrals of Chebyshev polynomials and metric element products.
!> \ingroup grp_integrals
!>
!> **Chebyshev-metric information**
!> <ul>
!> <li> The following quantities are calculated:
!>
!>       \f{eqnarray}{ \verb+DToocc(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}^\prime \; {\overline T}_{p,j}        \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \cos\alpha_i \cos\alpha_j                  \\
!>                     \verb+DToocs(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}^\prime \; {\overline T}_{p,j}        \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \cos\alpha_i \sin\alpha_j                  \\
!>                     \verb+DToosc(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}^\prime \; {\overline T}_{p,j}        \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \sin\alpha_i \cos\alpha_j                  \\
!>                     \verb+DTooss(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}^\prime \; {\overline T}_{p,j}        \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \sin\alpha_i \sin\alpha_j
!>       \f}
!>
!>       \f{eqnarray}{ \verb+TTsscc(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}        \; {\overline T}_{p,j}        \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \cos\alpha_i \cos\alpha_j \; \bar g_{ss} \\
!>                     \verb+TTsscs(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}        \; {\overline T}_{p,j}        \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \cos\alpha_i \sin\alpha_j \; \bar g_{ss} \\
!>                     \verb+TTsssc(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}        \; {\overline T}_{p,j}        \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \sin\alpha_i \cos\alpha_j \; \bar g_{ss} \\
!>                     \verb+TTssss(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}        \; {\overline T}_{p,j}        \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \sin\alpha_i \sin\alpha_j \; \bar g_{ss}
!>       \f}
!>
!>       \f{eqnarray}{ \verb+TDstcc(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}        \; {\overline T}_{p,j}^\prime \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \cos\alpha_i \cos\alpha_j \; \bar g_{s\theta} \\
!>                     \verb+TDstcs(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}        \; {\overline T}_{p,j}^\prime \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \cos\alpha_i \sin\alpha_j \; \bar g_{s\theta} \\
!>                     \verb+TDstsc(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}        \; {\overline T}_{p,j}^\prime \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \sin\alpha_i \cos\alpha_j \; \bar g_{s\theta} \\
!>                     \verb+TDstss(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}        \; {\overline T}_{p,j}^\prime \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \sin\alpha_i \sin\alpha_j \; \bar g_{s\theta}
!>       \f}
!>
!>       \f{eqnarray}{ \verb+TDstcc(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}        \; {\overline T}_{p,j}^\prime \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \cos\alpha_i \cos\alpha_j \; \bar g_{s\zeta} \\
!>                     \verb+TDstcs(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}        \; {\overline T}_{p,j}^\prime \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \cos\alpha_i \sin\alpha_j \; \bar g_{s\zeta} \\
!>                     \verb+TDstsc(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}        \; {\overline T}_{p,j}^\prime \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \sin\alpha_i \cos\alpha_j \; \bar g_{s\zeta} \\
!>                     \verb+TDstss(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}        \; {\overline T}_{p,j}^\prime \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \sin\alpha_i \sin\alpha_j \; \bar g_{s\zeta}
!>       \f}
!>
!>       \f{eqnarray}{ \verb+DDstcc(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}^\prime \; {\overline T}_{p,j}^\prime \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \cos\alpha_i \cos\alpha_j \; \bar g_{\theta\theta} \\
!>                     \verb+DDstcs(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}^\prime \; {\overline T}_{p,j}^\prime \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \cos\alpha_i \sin\alpha_j \; \bar g_{\theta\theta} \\
!>                     \verb+DDstsc(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}^\prime \; {\overline T}_{p,j}^\prime \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \sin\alpha_i \cos\alpha_j \; \bar g_{\theta\theta} \\
!>                     \verb+DDstss(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}^\prime \; {\overline T}_{p,j}^\prime \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \sin\alpha_i \sin\alpha_j \; \bar g_{\theta\theta}
!>       \f}
!>
!>       \f{eqnarray}{ \verb+DDstcc(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}^\prime \; {\overline T}_{p,j}^\prime \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \cos\alpha_i \cos\alpha_j \; \bar g_{\theta\zeta} \\
!>                     \verb+DDstcs(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}^\prime \; {\overline T}_{p,j}^\prime \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \cos\alpha_i \sin\alpha_j \; \bar g_{\theta\zeta} \\
!>                     \verb+DDstsc(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}^\prime \; {\overline T}_{p,j}^\prime \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \sin\alpha_i \cos\alpha_j \; \bar g_{\theta\zeta} \\
!>                     \verb+DDstss(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}^\prime \; {\overline T}_{p,j}^\prime \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \sin\alpha_i \sin\alpha_j \; \bar g_{\theta\zeta}
!>       \f}
!>
!>       \f{eqnarray}{ \verb+DDstcc(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}^\prime \; {\overline T}_{p,j}^\prime \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \cos\alpha_i \cos\alpha_j \; \bar g_{\zeta\zeta} \\
!>                     \verb+DDstcs(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}^\prime \; {\overline T}_{p,j}^\prime \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \cos\alpha_i \sin\alpha_j \; \bar g_{\zeta\zeta} \\
!>                     \verb+DDstsc(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}^\prime \; {\overline T}_{p,j}^\prime \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \sin\alpha_i \cos\alpha_j \; \bar g_{\zeta\zeta} \\
!>                     \verb+DDstss(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}^\prime \; {\overline T}_{p,j}^\prime \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \sin\alpha_i \sin\alpha_j \; \bar g_{\zeta\zeta}
!>       \f}
!>
!>       where \f${\overline T}_{l,i}\equiv T_l \, \bar s^{m_i/2}\f$ if the domain includes the coordinate singularity, and \f${\overline T}_{l,i}\equiv T_l\f$ if not;
!>       and \f$\bar g_{\mu\nu} \equiv g_{\mu\nu} / \sqrt g\f$. </li>
!>
!> <li> The double-angle formulae are used to reduce the above expressions to the Fourier harmonics of \f$\bar g_{\mu\nu}\f$:
!>       see \c kija and \c kijs, which are defined in preset.f90 . </li>
!>
!> </ul>
!>
!> @param[in] lquad degree of quadrature
!> @param[in] mn    number of Fourier harmonics
!> @param[in] lvol  index of nested volume
!> @param[in] lrad  order of Chebychev polynomials
subroutine ma00aa(lquad, mn, lvol, lrad)
    use mod_kinds, only: wp => dp
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    use constants, only: zero, half, one, two, pi, pi2

    use fileunits, only: ounit

    use inputlist, only: mpol, Wma00aa, Wmacros

    use cputiming, only: Tma00aa

    use allglobal, only: myid, ncpu, cpus, MPI_COMM_SPEC, &
                         Mvol, im, in, mne, &
                         YESstellsym, NOTstellsym, &
                         gaussianweight, gaussianabscissae, &
                         DToocc, DToocs, DToosc, DTooss, &
                         TTsscc, TTsscs, TTsssc, TTssss, &
                         TDstcc, TDstcs, TDstsc, TDstss, &
                         TDszcc, TDszcs, TDszsc, TDszss, &
                         DDttcc, DDttcs, DDttsc, DDttss, &
                         DDtzcc, DDtzcs, DDtzsc, DDtzss, &
                         DDzzcc, DDzzcs, DDzzsc, DDzzss, &
                         ki, kija, kijs, &
                         goomne, goomno, &
                         gssmne, gssmno, &
                         gstmne, gstmno, &
                         gszmne, gszmno, &
                         gttmne, gttmno, &
                         gtzmne, gtzmno, &
                         gzzmne, gzzmno, &
                         Lcoordinatesingularity, regumm, &
                         pi2pi2nfp, pi2pi2nfphalf, Lsavedguvij, &
                         dBdX

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef OPENMP
    USE OMP_LIB
#endif
    use mpi
    implicit none
    integer :: ierr, astat, ios, nthreads, ithread
    real(wp) :: cput, cpui, cpuo = 0 ! cpu time; cpu initial; cpu old; 31 Jan 13;

    integer, intent(in) :: lquad, mn, lvol, lrad

    integer :: jquad, ll, pp, ll1, pp1, uv, ii, jj, io, mn2, lp2, mn2_max, lp2_max, nele

    integer :: kk, kd, kka, kks, kda, kds, Lcurvature, ideriv

    real(wp) :: lss, jthweight, fee, feo, foe, foo, Tl, Dl, Tp, Dp, TlTp, TlDp, DlTp, DlDp, ikda, ikds, imn2, ilrad, lssm

    real(wp) :: foocc, foocs, foosc, fooss
    real(wp) :: fsscc, fsscs, fsssc, fssss
    real(wp) :: fstcc, fstcs, fstsc, fstss
    real(wp) :: fszcc, fszcs, fszsc, fszss
    real(wp) :: fttcc, fttcs, fttsc, fttss
    real(wp) :: ftzcc, ftzcs, ftzsc, ftzss
    real(wp) :: fzzcc, fzzcs, fzzsc, fzzss

    real(wp) :: sbar
    real(wp), allocatable :: basis(:, :, :, :)

    cpui = MPI_WTIME()
    cpuo = cpui
#ifdef OPENMP
    nthreads = omp_get_max_threads()
#else
    nthreads = 1
#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef DEBUG

    if (lvol .lt. 1 .or. lvol .gt. Mvol) then
        write (6, '("ma00aa :      fatal : myid=",i3," ; lvol.lt.1 .or. lvol.gt.Mvol ; illegal volume label ;")') myid
        call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
        stop "ma00aa : lvol.lt.1 .or. lvol.gt.Mvol : illegal volume label  ;"
    end if

#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    mn2_max = mn*mn
    lp2_max = (lrad + 1)*(lrad + 1)
    imn2 = one/real(mn)
    ilrad = one/real(lrad + 1)

    DToocc = zero
    TTssss = zero
    TDstsc = zero
    TDszsc = zero
    DDttcc = zero
    DDtzcc = zero
    DDzzcc = zero

    if (NOTstellsym) then
        DToocs = zero
        DToosc = zero
        DTooss = zero

        TTsscc = zero
        TTsscs = zero
        TTsssc = zero

        TDstcc = zero
        TDstcs = zero
        TDstss = zero

        TDszcc = zero
        TDszcs = zero
        TDszss = zero

        DDttcs = zero
        DDttsc = zero
        DDttss = zero

        DDtzcs = zero
        DDtzsc = zero
        DDtzss = zero

        DDzzcs = zero
        DDzzsc = zero
        DDzzss = zero
    end if !NOTstellsym

    allocate (basis(0:lrad, 0:mpol, 0:1, lquad), stat=astat)
    basis(0:lrad, 0:mpol, 0:1, lquad) = zero

    if (dBdX%L) then; Lcurvature = 3; ideriv = 1
    else; Lcurvature = 1; ideriv = 0
    end if

    if (.not. Lsavedguvij) then

        cput = MPI_WTIME()
        Tma00aa = Tma00aa + (cput - cpuo)
        call compute_guvijsave(lquad, lvol, ideriv, Lcurvature)
        cpuo = MPI_WTIME()

    end if

    cput = MPI_WTIME()
    Tma00aa = Tma00aa + (cput - cpuo)
    call metrix(lquad, lvol)
    cpuo = MPI_WTIME()
    ! compute metric elements; 16 Jan 13;

    do jquad = 1, lquad
        lss = gaussianabscissae(jquad, lvol); jthweight = gaussianweight(jquad, lvol)
        sbar = (lss + one)*half
        if (Lcoordinatesingularity) then
            call get_zernike(sbar, lrad, mpol, basis(:, :, 0:1, jquad)) ! use Zernike polynomials 29 Jun 19;
        else
            call get_cheby(lss, lrad, basis(:, 0, 0:1, jquad))
        end if
    end do

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!$OMP PARALLEL DO SHARED(lquad,lrad,lvol,mn,basis,mn2_max,lp2_max) PRIVATE(jquad,lss,jthweight,sbar,mn2,ii,jj,kka,kks,ikds,ikda,lp2,ll,pp,ll1,pp1,Tl,Tp,Dl,Dp,TlTP,Tldp,DlTp,DlDp,foocc,fssss,fstsc,fszsc,fttcc,ftzcc,fzzcc,foocs,foosc,fooss,fsscc,fsscs,fsssc,fstcc,fstcs,fstss,fszcc,fszcs,fszss,fttcs,fttsc,fttss,ftzcs,ftzsc,ftzss,fzzcs,fzzsc,fzzss)
    do mn2 = 1, mn2_max
        ii = mod(mn2 - 1, mn) + 1
        jj = (mn2 - ii)/mn + 1

        do jquad = 1, lquad ! Gaussian quadrature loop;

            lss = gaussianabscissae(jquad, lvol); jthweight = gaussianweight(jquad, lvol)
            sbar = (lss + one)*half

            kks = kijs(ii, jj, 0) !; kds = kijs(ii,jj,1)
            kka = kija(ii, jj, 0) !; kda = kija(ii,jj,1)
            ikds = jthweight/kijs(ii, jj, 1)
            ikda = jthweight/kija(ii, jj, 1)

            foocc = (+goomne(kks, jquad)*abs(ikds) + goomne(kka, jquad)*abs(ikda))
            fssss = (+gssmne(kks, jquad)*abs(ikds) - gssmne(kka, jquad)*abs(ikda))
            fstsc = (+gstmno(kks, jquad)*ikds + gstmno(kka, jquad)*ikda)
            fszsc = (+gszmno(kks, jquad)*ikds + gszmno(kka, jquad)*ikda)
            fttcc = (+gttmne(kks, jquad)*abs(ikds) + gttmne(kka, jquad)*abs(ikda))
            ftzcc = (+gtzmne(kks, jquad)*abs(ikds) + gtzmne(kka, jquad)*abs(ikda))
            fzzcc = (+gzzmne(kks, jquad)*abs(ikds) + gzzmne(kka, jquad)*abs(ikda))

            if (NOTstellsym) then
                foocs = (-goomno(kks, jquad)*ikds + goomno(kka, jquad)*ikda)
                foosc = (+goomno(kks, jquad)*ikds + goomno(kka, jquad)*ikda)
                fooss = (+goomne(kks, jquad)*abs(ikds) - goomne(kka, jquad)*abs(ikda))

                fsscc = (+gssmne(kks, jquad)*abs(ikds) + gssmne(kka, jquad)*abs(ikda))
                fsscs = (-gssmno(kks, jquad)*ikds + gssmno(kka, jquad)*ikda)
                fsssc = (+gssmno(kks, jquad)*ikds + gssmno(kka, jquad)*ikda)

                fstcc = (+gstmne(kks, jquad)*abs(ikds) + gstmne(kka, jquad)*abs(ikda))
                fstcs = (-gstmno(kks, jquad)*ikds + gstmno(kka, jquad)*ikda)
                fstss = (+gstmne(kks, jquad)*abs(ikds) - gstmne(kka, jquad)*abs(ikda))

                fszcc = (+gszmne(kks, jquad)*abs(ikds) + gszmne(kka, jquad)*abs(ikda))
                fszcs = (-gszmno(kks, jquad)*ikds + gszmno(kka, jquad)*ikda)
                fszss = (+gszmne(kks, jquad)*abs(ikds) - gszmne(kka, jquad)*abs(ikda))

                fttcs = (-gttmno(kks, jquad)*ikds + gttmno(kka, jquad)*ikda)
                fttsc = (+gttmno(kks, jquad)*ikds + gttmno(kka, jquad)*ikda)
                fttss = (+gttmne(kks, jquad)*abs(ikds) - gttmne(kka, jquad)*abs(ikda))

                ftzcs = (-gtzmno(kks, jquad)*ikds + gtzmno(kka, jquad)*ikda)
                ftzsc = (+gtzmno(kks, jquad)*ikds + gtzmno(kka, jquad)*ikda)
                ftzss = (+gtzmne(kks, jquad)*abs(ikds) - gtzmne(kka, jquad)*abs(ikda))

                fzzcs = (-gzzmno(kks, jquad)*ikds + gzzmno(kka, jquad)*ikda)
                fzzsc = (+gzzmno(kks, jquad)*ikds + gzzmno(kka, jquad)*ikda)
                fzzss = (+gzzmne(kks, jquad)*abs(ikds) - gzzmne(kka, jquad)*abs(ikda))
            end if !NOTstellsym

            do lp2 = 1, lp2_max
                ll = mod(lp2 - 1, lrad + 1)
                pp = (lp2 - ll - 1)/(lrad + 1)

                if (Lcoordinatesingularity) then

                    ll1 = (ll - mod(ll, 2))/2 ! shrinked dof for Zernike; 02 Jul 19
                    pp1 = (pp - mod(pp, 2))/2 ! shrinked dof for Zernike; 02 Jul 19

                    if (ll < im(ii)) cycle ! zernike only non-zero for ll>=ii
                    if (pp < im(jj)) cycle ! zernike only non-zero for pp>=jj
                    if (mod(ll + im(ii), 2) /= 0) cycle ! zernike only non-zero if ll and ii have the same parity
                    if (mod(pp + im(jj), 2) /= 0) cycle ! zernike only non-zero if pp and jj have the same parity

                    Tl = basis(ll, im(ii), 0, jquad)         ! use Zernike polynomials 29 Jun 19;
                    Dl = basis(ll, im(ii), 1, jquad)*half  ! use Zernike polynomials 29 Jun 19;

                    Tp = basis(pp, im(jj), 0, jquad)         ! use Zernike polynomials 29 Jun 19;
                    Dp = basis(pp, im(jj), 1, jquad)*half  ! use Zernike polynomials 29 Jun 19;

                else

                    ll1 = ll
                    pp1 = pp

                    Tl = basis(ll, 0, 0, jquad)
                    Dl = basis(ll, 0, 1, jquad)

                    Tp = basis(pp, 0, 0, jquad)
                    Dp = basis(pp, 0, 1, jquad)

                end if ! Lcoordinatesingularity

                TlTp = Tl*Tp
                TlDp = Tl*Dp
                DlTp = Dl*Tp
                DlDp = Dl*Dp

                DToocc(ll1, pp1, ii, jj) = DToocc(ll1, pp1, ii, jj) + DlTp*foocc
                TTssss(ll1, pp1, ii, jj) = TTssss(ll1, pp1, ii, jj) + TlTp*fssss
                TDstsc(ll1, pp1, ii, jj) = TDstsc(ll1, pp1, ii, jj) + TlDp*fstsc
                TDszsc(ll1, pp1, ii, jj) = TDszsc(ll1, pp1, ii, jj) + TlDp*fszsc
                DDttcc(ll1, pp1, ii, jj) = DDttcc(ll1, pp1, ii, jj) + DlDp*fttcc
                DDtzcc(ll1, pp1, ii, jj) = DDtzcc(ll1, pp1, ii, jj) + DlDp*ftzcc
                DDzzcc(ll1, pp1, ii, jj) = DDzzcc(ll1, pp1, ii, jj) + DlDp*fzzcc

                if (NOTstellsym) then

                    DToocs(ll1, pp1, ii, jj) = DToocs(ll1, pp1, ii, jj) + DlTp*foocs
                    DToosc(ll1, pp1, ii, jj) = DToosc(ll1, pp1, ii, jj) + DlTp*foosc
                    DTooss(ll1, pp1, ii, jj) = DTooss(ll1, pp1, ii, jj) + DlTp*fooss

                    TTsscc(ll1, pp1, ii, jj) = TTsscc(ll1, pp1, ii, jj) + TlTp*fsscc
                    TTsscs(ll1, pp1, ii, jj) = TTsscs(ll1, pp1, ii, jj) + TlTp*fsscs
                    TTsssc(ll1, pp1, ii, jj) = TTsssc(ll1, pp1, ii, jj) + TlTp*fsssc

                    TDstcc(ll1, pp1, ii, jj) = TDstcc(ll1, pp1, ii, jj) + TlDp*fstcc
                    TDstcs(ll1, pp1, ii, jj) = TDstcs(ll1, pp1, ii, jj) + TlDp*fstcs
                    TDstss(ll1, pp1, ii, jj) = TDstss(ll1, pp1, ii, jj) + TlDp*fstss

                    TDszcc(ll1, pp1, ii, jj) = TDszcc(ll1, pp1, ii, jj) + TlDp*fszcc
                    TDszcs(ll1, pp1, ii, jj) = TDszcs(ll1, pp1, ii, jj) + TlDp*fszcs
                    TDszss(ll1, pp1, ii, jj) = TDszss(ll1, pp1, ii, jj) + TlDp*fszss

                    DDttcs(ll1, pp1, ii, jj) = DDttcs(ll1, pp1, ii, jj) + DlDp*fttcs
                    DDttsc(ll1, pp1, ii, jj) = DDttsc(ll1, pp1, ii, jj) + DlDp*fttsc
                    DDttss(ll1, pp1, ii, jj) = DDttss(ll1, pp1, ii, jj) + DlDp*fttss

                    DDtzcs(ll1, pp1, ii, jj) = DDtzcs(ll1, pp1, ii, jj) + DlDp*ftzcs
                    DDtzsc(ll1, pp1, ii, jj) = DDtzsc(ll1, pp1, ii, jj) + DlDp*ftzsc
                    DDtzss(ll1, pp1, ii, jj) = DDtzss(ll1, pp1, ii, jj) + DlDp*ftzss

                    DDzzcs(ll1, pp1, ii, jj) = DDzzcs(ll1, pp1, ii, jj) + DlDp*fzzcs
                    DDzzsc(ll1, pp1, ii, jj) = DDzzsc(ll1, pp1, ii, jj) + DlDp*fzzsc
                    DDzzss(ll1, pp1, ii, jj) = DDzzss(ll1, pp1, ii, jj) + DlDp*fzzss
                end if !NOTstellsym

            end do ! end of do lp2; 08 Feb 16;

        end do ! end of do jquad

    end do ! end of do mn
!$OMP END PARALLEL DO

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    deallocate (basis, stat=astat)

    DToocc = DToocc*pi2pi2nfphalf
    TTssss = TTssss*pi2pi2nfphalf
    TDstsc = TDstsc*pi2pi2nfphalf
    TDszsc = TDszsc*pi2pi2nfphalf
    DDttcc = DDttcc*pi2pi2nfphalf
    DDtzcc = DDtzcc*pi2pi2nfphalf
    DDzzcc = DDzzcc*pi2pi2nfphalf

    if (NOTstellsym) then
        DToocs = DToocs*pi2pi2nfphalf
        DToosc = DToosc*pi2pi2nfphalf
        DTooss = DTooss*pi2pi2nfphalf

        TTsscc = TTsscc*pi2pi2nfphalf
        TTsscs = TTsscs*pi2pi2nfphalf
        TTsssc = TTsssc*pi2pi2nfphalf

        TDstcc = TDstcc*pi2pi2nfphalf
        TDstcs = TDstcs*pi2pi2nfphalf
        TDstss = TDstss*pi2pi2nfphalf

        TDszcc = TDszcc*pi2pi2nfphalf
        TDszcs = TDszcs*pi2pi2nfphalf
        TDszss = TDszss*pi2pi2nfphalf

        DDttcs = DDttcs*pi2pi2nfphalf
        DDttsc = DDttsc*pi2pi2nfphalf
        DDttss = DDttss*pi2pi2nfphalf

        DDtzcs = DDtzcs*pi2pi2nfphalf
        DDtzsc = DDtzsc*pi2pi2nfphalf
        DDtzss = DDtzss*pi2pi2nfphalf

        DDzzcs = DDzzcs*pi2pi2nfphalf
        DDzzsc = DDzzsc*pi2pi2nfphalf
        DDzzss = DDzzss*pi2pi2nfphalf

    end if !NOTstellsym

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

9999 continue
    cput = MPI_WTIME()
    Tma00aa = Tma00aa + (cput - cpuo)
    return

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine ma00aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
