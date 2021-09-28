!> \file
!> \brief Calculates volume integrals of Chebyshev-polynomials and metric elements for preconditioner.

!> \brief Calculates volume integrals of Chebyshev-polynomials and metric elements for preconditioner.
!> \ingroup grp_integrals
!>
!> Computes the integrals needed for spsmat.f90. Same as ma00aa.f90, but only compute the relevant terms that are non-zero.
!>
!> @param lquad
!> @param mn
!> @param lvol
!> @param lrad
subroutine spsint( lquad, mn, lvol, lrad )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  use constants, only : zero, half, one, two, pi, pi2

  use numerical, only : vsmall, small, sqrtmachprec

  use fileunits, only : ounit

  use inputlist, only : mpol, Wspsint, Wmacros

  use cputiming, only : Tspsint

  use allglobal, only : myid, ncpu, cpus, MPI_COMM_SPEC, &
                        Mvol, im, in, mne, Ntz, &
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
                        Lcoordinatesingularity, regumm, &
                        pi2pi2nfp, pi2pi2nfphalf, &
                        guvijsave


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS

  INTEGER, intent(in) :: lquad, mn, lvol, lrad

  INTEGER             :: jquad, ll, pp, ll1, pp1, uv, ii, jj, io, mn2, lp2, mn2_max, lp2_max, nele, mi

  INTEGER             :: kk, kd, kka, kks, kda, kds

  REAL                :: lss, jthweight, fee, feo, foe, foo, Tl, Dl, Tp, Dp, TlTp, TlDp, DlTp, DlDp, ikda, ikds, imn2, ilrad, lssm

  REAL                :: foocc, fooss
  REAL                :: fsscc, fssss
  REAL                :: fstcc, fstss
  REAL                :: fszcc, fszss
  REAL                :: fttcc, fttss
  REAL                :: ftzcc, ftzss
  REAL                :: fzzcc, fzzss

  REAL                :: goomne, gssmne, gstmne, gszmne, gttmne, gtzmne, gzzmne

  REAL                :: sbar

  REAL, allocatable   :: basis(:,:,:,:)

  BEGIN( spsint )
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  mn2_max = mn*mn
  lp2_max = (lrad+1)*(lrad+1)
  imn2    =  one/real(mn)
  ilrad = one/real(lrad+1)

  DToocc = zero
  TTssss = zero
  DDttcc = zero
  DDtzcc = zero
  DDzzcc = zero

  if (NOTstellsym) then
    DTooss = zero
    TTsscc = zero
    TDstcc = zero
    TDstss = zero
    TDszcc = zero
    TDszss = zero
    DDttss = zero
    DDtzss = zero
    DDzzss = zero
  endif !NOTstellsym

  SALLOCATE(basis,     (0:lrad,0:mpol,0:1,lquad), zero)
  do jquad = 1, lquad
    lss = gaussianabscissae(jquad,lvol) ; jthweight = gaussianweight(jquad,lvol)
    sbar = (lss + one) * half
    if (Lcoordinatesingularity) then
      call get_zernike(sbar, lrad, mpol, basis(:,:,0:1,jquad)) ! use Zernike polynomials 29 Jun 19;
    else
      call get_cheby(lss, lrad, basis(:,0,0:1,jquad))
    endif
  enddo
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!$OMP PARALLEL DO PRIVATE(TlTp,TlDp,DlTp,DlDp,Tl,Dl,Tp,Dp,ll1,pp1,ll,pp,lss,jthweight,sbar,goomne,gssmne,gstmne,gszmne,gttmne,gtzmne,gzzmne,foocc,fooss,fsscc,fssss,fstcc,fstss,fszcc,fszss,fttcc,fttss,ftzcc,ftzss,fzzcc,fzzss) SHARED(lquad,lp2_max,lrad,basis)
  do jquad = 1, lquad ! Gaussian quadrature loop;

    lss = gaussianabscissae(jquad,lvol) ; jthweight = gaussianweight(jquad,lvol)

    ! compute metric

    goomne = one
    gssmne = SUM(guvijsave(1:Ntz,1,1,jquad)) / real(Ntz)
    gstmne = SUM(guvijsave(1:Ntz,1,2,jquad)) / real(Ntz)
    gszmne = SUM(guvijsave(1:Ntz,1,3,jquad)) / real(Ntz)
    gttmne = SUM(guvijsave(1:Ntz,2,2,jquad)) / real(Ntz)
    gtzmne = SUM(guvijsave(1:Ntz,2,3,jquad)) / real(Ntz)
    gzzmne = SUM(guvijsave(1:Ntz,3,3,jquad)) / real(Ntz)

    foocc = goomne * jthweight
    fssss = gssmne * jthweight
    fttcc = gttmne * jthweight
    ftzcc = gtzmne * jthweight
    fzzcc = gzzmne * jthweight


    if (NOTstellsym) then
      fooss = goomne * jthweight
      fsscc = gssmne * jthweight
      fstcc = gstmne * jthweight
      fstss = gstmne * jthweight
      fszcc = gszmne * jthweight
      fszss = gszmne * jthweight
      fttss = gttmne * jthweight
      ftzss = gtzmne * jthweight
      fzzss = gzzmne * jthweight
    end if !NOTstellsym

    do mi = 0, mpol

      do lp2 = 1, lp2_max
        ll = mod(lp2-1,lrad+1)
        pp = (lp2-ll-1)/(lrad+1)

        if (Lcoordinatesingularity) then

          ll1 = (ll - mod(ll,2))/2 ! shrinked dof for Zernike; 02 Jul 19
          pp1 = (pp - mod(pp,2))/2 ! shrinked dof for Zernike; 02 Jul 19

          if (ll < mi) cycle ! zernike only non-zero for ll>=ii
          if (pp < mi) cycle ! zernike only non-zero for pp>=jj
          if (mod(ll+mi,2)/=0) cycle ! zernike only non-zero if ll and ii have the same parity
          if (mod(pp+mi,2)/=0) cycle ! zernike only non-zero if pp and jj have the same parity

          Tl = basis(ll, mi, 0, jquad)         ! use Zernike polynomials 29 Jun 19;
          Dl = basis(ll, mi, 1, jquad) * half  ! use Zernike polynomials 29 Jun 19;

          Tp = basis(pp, mi, 0, jquad)         ! use Zernike polynomials 29 Jun 19;
          Dp = basis(pp, mi, 1, jquad) * half  ! use Zernike polynomials 29 Jun 19;

        else

          if (mi .gt. 0) cycle ! we don't nee to compute more than m>0

          ll1 = ll
          pp1 = pp

          Tl = basis(ll, 0, 0, jquad) ! Cheby
          Dl = basis(ll, 0, 1, jquad) ! Cheby

          Tp = basis(pp, 0, 0, jquad) ! Cheby
          Dp = basis(pp, 0, 1, jquad) ! Cheby

        end if ! if (Lcoordinatesingularity)

        TlTp = Tl * Tp
        TlDp = Tl * Dp
        DlTp = Dl * Tp
        DlDp = Dl * Dp

!$OMP ATOMIC UPDATE
        DToocc( ll1, pp1, mi+1, 1 ) = DToocc( ll1, pp1, mi+1, 1 ) + DlTp * foocc
!$OMP ATOMIC UPDATE
        TTssss( ll1, pp1, mi+1, 1 ) = TTssss( ll1, pp1, mi+1, 1 ) + TlTp * fssss
!$OMP ATOMIC UPDATE
        DDttcc( ll1, pp1, mi+1, 1 ) = DDttcc( ll1, pp1, mi+1, 1 ) + DlDp * fttcc
!$OMP ATOMIC UPDATE
        DDtzcc( ll1, pp1, mi+1, 1 ) = DDtzcc( ll1, pp1, mi+1, 1 ) + DlDp * ftzcc
!$OMP ATOMIC UPDATE
        DDzzcc( ll1, pp1, mi+1, 1 ) = DDzzcc( ll1, pp1, mi+1, 1 ) + DlDp * fzzcc

        if (NOTstellsym) then
!$OMP ATOMIC UPDATE
          DTooss( ll1, pp1, mi+1, 1 ) = DTooss( ll1, pp1, mi+1, 1 ) + DlTp * fooss
!$OMP ATOMIC UPDATE
          TTsscc( ll1, pp1, mi+1, 1 ) = TTsscc( ll1, pp1, mi+1, 1 ) + TlTp * fsscc
!$OMP ATOMIC UPDATE
          TDstcc( ll1, pp1, mi+1, 1 ) = TDstcc( ll1, pp1, mi+1, 1 ) + TlDp * fstcc
!$OMP ATOMIC UPDATE
          TDszcc( ll1, pp1, mi+1, 1 ) = TDszcc( ll1, pp1, mi+1, 1 ) + TlDp * fszcc
!$OMP ATOMIC UPDATE
          TDszss( ll1, pp1, mi+1, 1 ) = TDszss( ll1, pp1, mi+1, 1 ) + TlDp * fszss
!$OMP ATOMIC UPDATE
          DDttss( ll1, pp1, mi+1, 1 ) = DDttss( ll1, pp1, mi+1, 1 ) + DlDp * fttss
!$OMP ATOMIC UPDATE
          DDtzss( ll1, pp1, mi+1, 1 ) = DDtzss( ll1, pp1, mi+1, 1 ) + DlDp * ftzss
!$OMP ATOMIC UPDATE
          DDzzss( ll1, pp1, mi+1, 1 ) = DDzzss( ll1, pp1, mi+1, 1 ) + DlDp * fzzss
        end if !NOTstellsym

      enddo ! end of do lp2; 08 Feb 16;

    enddo ! end of do mi

  enddo ! end of do jquad; ! 16 Jan 13;
!$OMP END PARALLEL DO
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  nele = SIZE(TTssss)

  call DSCAL(nele, pi2pi2nfphalf, DToocc, 1)
  call DSCAL(nele, pi2pi2nfphalf, TTssss, 1)
  call DSCAL(nele, pi2pi2nfphalf, DDttcc, 1)
  call DSCAL(nele, pi2pi2nfphalf, DDtzcc, 1)
  call DSCAL(nele, pi2pi2nfphalf, DDzzcc, 1)

  if (NOTstellsym) then

    call DSCAL(nele, pi2pi2nfphalf, DTooss, 1)
    call DSCAL(nele, pi2pi2nfphalf, TTsscc, 1)
    call DSCAL(nele, pi2pi2nfphalf, TDstcc, 1)
    call DSCAL(nele, pi2pi2nfphalf, TDstss, 1)
    call DSCAL(nele, pi2pi2nfphalf, TDszcc, 1)
    call DSCAL(nele, pi2pi2nfphalf, TDszss, 1)
    call DSCAL(nele, pi2pi2nfphalf, DDttss, 1)
    call DSCAL(nele, pi2pi2nfphalf, DDtzss, 1)
    call DSCAL(nele, pi2pi2nfphalf, DDzzss, 1)

  end if !NOTstellsym

  DALLOCATE(basis)
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN( spsint )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine spsint

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
