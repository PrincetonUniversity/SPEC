!> \file
!> \brief Constructs matrices that represent the Beltrami linear system, matrix-free.

!> \brief Constructs matrices that represent the Beltrami linear system, matrix-free.
!> \ingroup grp_build_matrices
!>
!> @param lvol
!> @param mn
!> @param lrad
!> @param resultA
!> @param resultD
!> @param idx
subroutine mtrxhs(lvol, mn, lrad, resultA, resultD, idx)
    use mod_kinds, only: wp => dp
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    use constants, only: zero, one, two, half

    use numerical, only: small

    use fileunits, only: ounit

    use inputlist, only: Wmacros, Wmtrxhs, Mpol

    use cputiming, only: Tmtrxhs

    use allglobal, only: ncpu, myid, cpus, MPI_COMM_SPEC, &
                         YESstellsym, NOTstellsym, &
                         im, in, &
                         NAdof, &
                         Ate, Ato, Aze, Azo, &
                         Lcoordinatesingularity, TT, RTT, RTM, &
                         Tss, Tsc, Dts, Dtc, Dzs, Dzc, &
                         Ttc, Tts, Tzc, Tzs, &
                         dBdX, &
                         Lma, Lmb, Lmc, Lmd, Lme, Lmf, Lmg, Lmh, &
                         Lmavalue, Lmbvalue, Lmcvalue, Lmdvalue, &
                         Lmevalue, Lmfvalue, Lmgvalue, Lmhvalue, &
                         TT, RTT, RTM

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef OPENMP
    USE OMP_LIB
#endif
    use mpi
    implicit none
    integer :: ierr, astat, ios, nthreads, ithread
    real(wp) :: cput, cpui, cpuo = 0 ! cpu time; cpu initial; cpu old; 31 Jan 13;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    integer, intent(in) :: lvol, mn, lrad, idx

    real(wp), intent(out) :: resultA(0:NAdof(lvol)), resultD(0:NAdof(lvol))

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    integer :: NN, ii, ll, jj, ll1, mi, ni, id, jd, kk

    real(wp) :: Wte, Wto, Wze, Wzo, Hte, Hto, Hze, Hzo

    real(wp), allocatable :: TTdata(:, :, :), TTMdata(:, :)

    cpui = MPI_WTIME()
    cpuo = cpui
#ifdef OPENMP
    nthreads = omp_get_max_threads()
#else
    nthreads = 1
#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    NN = NAdof(lvol) ! shorthand;

    ! making use of only the zeroth component of dMA
    resultA(0:NN) = zero
    resultD(0:NN) = zero

    allocate (TTdata(0:lrad, 0:mpol, 0:1), stat=astat)
    TTdata(0:lrad, 0:mpol, 0:1) = zero

    allocate (TTMdata(0:lrad, 0:mpol), stat=astat)
    TTMdata(0:lrad, 0:mpol) = zero

    ! fill in Zernike/Chebyshev polynomials depending on Lcooridnatesingularity
    if (Lcoordinatesingularity) then
        TTdata(0:lrad, 0:mpol, 0:1) = RTT(0:lrad, 0:mpol, 0:1, 0)
        TTMdata(0:lrad, 0:mpol) = RTM(0:lrad, 0:mpol)
    else
        do ii = 0, mpol
            TTdata(0:lrad, ii, 0:1) = TT(0:lrad, 0:1, 0)
            TTMdata(0:lrad, ii) = TT(0:lrad, 0, 0)
        end do
    end if

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    if (YESstellsym) then
!$OMP PARALLEL DO PRIVATE(ii,mi,ni,ll,ll1,kk,Wte,Wze,Hte,Hze,id,jd) SHARED(mn,lrad,resultA,resultD,TTMdata,TTdata)
        do ii = 1, mn; mi = im(ii); ni = in(ii)

            do ll = 0, lrad

                if (Lcoordinatesingularity) then
                    if (ll < mi) cycle ! rule out zero components of Zernike;
                    if (mod(ll + mi, 2) > 0) cycle ! rule out zero components of Zernike;
                    ll1 = (ll - mod(ll, 2))/2 ! shrinked dof for Zernike; 02 Jul 19
                else
                    ll1 = ll
                end if

                Wte = -two*ni*Tss(ll1, ii) + two*Dzc(ll1, ii)
                Wze = -two*mi*Tss(ll1, ii) - two*Dtc(ll1, ii)

                id = Ate(lvol, 0, ii)%i(ll); resultA(id) = Wte
                id = Aze(lvol, 0, ii)%i(ll); resultA(id) = Wze
                if (dBdX%L) cycle ! the derivatives of Lagrange multipliers and helicity w.r.t. interface is zero

                Hte = Ttc(ll1, ii)
                Hze = Tzc(ll1, ii)

                id = Ate(lvol, 0, ii)%i(ll); resultD(id) = Hte
                id = Aze(lvol, 0, ii)%i(ll); resultD(id) = Hze

                if (Lcoordinatesingularity .and. ii .eq. 1) then; kk = 1
                else; kk = 0
                end if

                ! add Langrange multipliers
                ; ; id = Ate(lvol, 0, ii)%i(ll); resultA(id) = resultA(id) + Lmavalue(lvol, ii)*TTMdata(ll, mi)
                ; ; id = Aze(lvol, 0, ii)%i(ll); resultA(id) = resultA(id) + Lmbvalue(lvol, ii)*TTdata(ll, mi, kk)
                if (ii .gt. 1) then; id = Ate(lvol, 0, ii)%i(ll); resultA(id) = resultA(id) + Lmevalue(lvol, ii)*(-ni*TTdata(ll, mi, 1))
                    ; ; id = Aze(lvol, 0, ii)%i(ll); resultA(id) = resultA(id) + Lmevalue(lvol, ii)*(-mi*TTdata(ll, mi, 1))
                else; id = Ate(lvol, 0, ii)%i(ll); resultA(id) = resultA(id) + Lmgvalue(lvol, ii)*TTdata(ll, mi, 1)
                    ; ; id = Aze(lvol, 0, ii)%i(ll); resultA(id) = resultA(id) + Lmhvalue(lvol, ii)*TTdata(ll, mi, 1)
                end if

                ; ; id = Lma(lvol, ii); resultA(id) = resultA(id) + Ate(lvol, idx, ii)%s(ll)*TTMdata(ll, mi)
                ; ; id = Lmb(lvol, ii); resultA(id) = resultA(id) + Aze(lvol, idx, ii)%s(ll)*TTdata(ll, mi, kk)
                if (ii .gt. 1) then; id = Lme(lvol, ii); resultA(id) = resultA(id) + Ate(lvol, idx, ii)%s(ll)*(-ni*TTdata(ll, mi, 1))
                    ; ; id = Lme(lvol, ii); resultA(id) = resultA(id) + Aze(lvol, idx, ii)%s(ll)*(-mi*TTdata(ll, mi, 1))
                else; id = Lmg(lvol, ii); resultA(id) = resultA(id) + Ate(lvol, idx, ii)%s(ll)*TTdata(ll, mi, 1)
                    ; ; id = Lmh(lvol, ii); resultA(id) = resultA(id) + Aze(lvol, idx, ii)%s(ll)*TTdata(ll, mi, 1)
                end if

            end do ! end of do ll ;

        end do ! end of do ii ;
!$OMP END PARALLEL DO

    else ! NOTstellsym ;
!$OMP PARALLEL DO PRIVATE(ii,mi,ni,ll,ll1,kk,Wte,Wze,Wzo,Wto,Hte,Hze,Hto,Hzo,id,jd) SHARED(mn,lrad,resultA,resultD,TTMdata,TTdata)
        do ii = 1, mn; mi = im(ii); ni = in(ii)

            do ll = 0, lrad

                if (Lcoordinatesingularity) then
                    if (ll < mi) cycle ! rule out zero components of Zernike;
                    if (mod(ll + mi, 2) > 0) cycle ! rule out zero components of Zernike;
                    ll1 = (ll - mod(ll, 2))/2 ! shrinked dof for Zernike; 02 Jul 19
                else
                    ll1 = ll
                end if

                Wte = -two*ni*Tss(ll1, ii) + two*Dzc(ll1, ii)
                Wze = -two*mi*Tss(ll1, ii) - two*Dtc(ll1, ii)
                Wto = +two*ni*Tsc(ll1, ii) + two*Dzs(ll1, ii)
                Wzo = +two*mi*Tsc(ll1, ii) - two*Dts(ll1, ii)

                id = Ate(lvol, 0, ii)%i(ll); resultA(id) = Wte
                id = Aze(lvol, 0, ii)%i(ll); resultA(id) = Wze
                id = Ato(lvol, 0, ii)%i(ll); resultA(id) = Wto
                id = Azo(lvol, 0, ii)%i(ll); resultA(id) = Wzo

                if (dBdX%L) cycle ! the derivatives of Lagrange multipliers and helicity w.r.t. interface is zero

                Hte = Ttc(ll1, ii)
                Hze = Tzc(ll1, ii)
                Hto = Tts(ll1, ii)
                Hzo = Tzs(ll1, ii)

                id = Ate(lvol, 0, ii)%i(ll); resultD(id) = Hte
                id = Aze(lvol, 0, ii)%i(ll); resultD(id) = Hze
                id = Ato(lvol, 0, ii)%i(ll); resultD(id) = Hto
                id = Azo(lvol, 0, ii)%i(ll); resultD(id) = Hzo

                if (Lcoordinatesingularity .and. ii .eq. 1) then; kk = 1
                else; kk = 0
                end if

                ! add Langrange multipliers
                ; ; id = Ate(lvol, 0, ii)%i(ll); resultA(id) = resultA(id) + Lmavalue(lvol, ii)*TTMdata(ll, mi)
                ; ; id = Aze(lvol, 0, ii)%i(ll); resultA(id) = resultA(id) + Lmbvalue(lvol, ii)*TTdata(ll, mi, kk)
                ; ; id = Ato(lvol, 0, ii)%i(ll); resultA(id) = resultA(id) + Lmcvalue(lvol, ii)*TTMdata(ll, mi)
                ; ; id = Azo(lvol, 0, ii)%i(ll); resultA(id) = resultA(id) + Lmdvalue(lvol, ii)*TTdata(ll, mi, kk)
                if (ii .gt. 1) then; id = Ate(lvol, 0, ii)%i(ll); resultA(id) = resultA(id) + Lmevalue(lvol, ii)*(-ni*TTdata(ll, mi, 1))
                    ; ; id = Aze(lvol, 0, ii)%i(ll); resultA(id) = resultA(id) + Lmevalue(lvol, ii)*(-mi*TTdata(ll, mi, 1))
                    ; ; id = Ato(lvol, 0, ii)%i(ll); resultA(id) = resultA(id) + Lmfvalue(lvol, ii)*(+ni*TTdata(ll, mi, 1))
                    ; ; id = Azo(lvol, 0, ii)%i(ll); resultA(id) = resultA(id) + Lmfvalue(lvol, ii)*(+mi*TTdata(ll, mi, 1))
                else; id = Ate(lvol, 0, ii)%i(ll); resultA(id) = resultA(id) + Lmgvalue(lvol, ii)*TTdata(ll, mi, 1)
                    ; ; id = Aze(lvol, 0, ii)%i(ll); resultA(id) = resultA(id) + Lmhvalue(lvol, ii)*TTdata(ll, mi, 1)
                end if

                ; ; id = Lma(lvol, ii); resultA(id) = resultA(id) + Ate(lvol, idx, ii)%s(ll)*TTMdata(ll, mi)
                ; ; id = Lmb(lvol, ii); resultA(id) = resultA(id) + Aze(lvol, idx, ii)%s(ll)*TTdata(ll, mi, kk)
                ; ; id = Lmc(lvol, ii); resultA(id) = resultA(id) + Ato(lvol, idx, ii)%s(ll)*TTMdata(ll, mi)
                ; ; id = Lmd(lvol, ii); resultA(id) = resultA(id) + Azo(lvol, idx, ii)%s(ll)*TTdata(ll, mi, kk)
                if (ii .gt. 1) then; id = Lme(lvol, ii); resultA(id) = resultA(id) + Ate(lvol, idx, ii)%s(ll)*(-ni*TTdata(ll, mi, 1))
                    ; ; id = Lme(lvol, ii); resultA(id) = resultA(id) + Aze(lvol, idx, ii)%s(ll)*(-mi*TTdata(ll, mi, 1))
                    ; ; id = Lmf(lvol, ii); resultA(id) = resultA(id) + Ato(lvol, idx, ii)%s(ll)*(+ni*TTdata(ll, mi, 1))
                    ; ; id = Lmf(lvol, ii); resultA(id) = resultA(id) + Azo(lvol, idx, ii)%s(ll)*(+mi*TTdata(ll, mi, 1))
                else; id = Lmg(lvol, ii); resultA(id) = resultA(id) + Ate(lvol, idx, ii)%s(ll)*TTdata(ll, mi, 1)
                    ; ; id = Lmh(lvol, ii); resultA(id) = resultA(id) + Aze(lvol, idx, ii)%s(ll)*TTdata(ll, mi, 1)
                end if

            end do ! end of do ll ;

        end do ! end of do ii ;
!$OMP END PARALLEL DO
    end if ! end of if( YESstellsym ) ;

    deallocate (TTdata, stat=astat)

    deallocate (TTMdata, stat=astat)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

9999 continue
    cput = MPI_WTIME()
    Tmtrxhs = Tmtrxhs + (cput - cpuo)
    return

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine mtrxhs

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
