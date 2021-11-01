!> \defgroup grp_metrics Metric quantities
!>
!> \file
!> \brief Calculates the metric quantities, \f$\sqrt g \, g^{\mu\nu}\f$, which are required for the energy and helicity integrals.

!> \brief Calculates the metric quantities, \f$\sqrt g \, g^{\mu\nu}\f$, which are required for the energy and helicity integrals.
!> \ingroup grp_metrics
!>
!> **metrics**
!>
!> <ul>
!> <li> The Jacobian, \f$\sqrt g\f$, and the "lower" metric elements, \f$g_{\mu\nu}\f$, are calculated by coords(),
!>      and are provided on a regular grid in "real-space", i.e. \f$(\theta,\zeta)\f$, at a given radial location, i.e. where \f$s\f$ is input. </li>
!> </ul>
!>
!> **plasma region**
!>
!> <ul>
!> <li>  In the plasma region, the required terms are \f$\bar g_{\mu\nu} \equiv g_{\mu\nu}/\sqrt g\f$.
!>
!>       \f{eqnarray}{ \begin{array}{ccccccccccccccccccccccccc}
!>       \sqrt g \; g^{s     s     } & = & \left( g_{\theta\theta} g_{\zeta \zeta } - g_{\theta\zeta } g_{\theta\zeta } \right) / \sqrt g \\
!>       \sqrt g \; g^{s     \theta} & = & \left( g_{\theta\zeta } g_{s     \zeta } - g_{s     \theta} g_{\zeta \zeta } \right) / \sqrt g \\
!>       \sqrt g \; g^{s     \zeta } & = & \left( g_{s     \theta} g_{\theta\zeta } - g_{\theta\theta} g_{s     \zeta } \right) / \sqrt g \\
!>       \sqrt g \; g^{\theta\theta} & = & \left( g_{\zeta \zeta } g_{s     s     } - g_{s     \zeta } g_{s     \zeta } \right) / \sqrt g \\
!>       \sqrt g \; g^{\theta\zeta } & = & \left( g_{s     \zeta } g_{s     \theta} - g_{\theta\zeta } g_{s     s     } \right) / \sqrt g \\
!>       \sqrt g \; g^{\zeta \zeta } & = & \left( g_{s     s     } g_{\theta\theta} - g_{s     \theta} g_{s     \theta} \right) / \sqrt g
!>       \end{array}
!>       \f} </li>
!> </ul>
!>
!> **FFTs**
!>
!> <ul>
!> <li> After constructing the required quantities in real space, FFTs provided the required Fourier harmonics, which are returned through global.f90 .
!>      (The "extended" Fourier resolution is used.) </li>
!> </ul>
subroutine metrix( lquad, lvol )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  use constants, only : zero, one

  use numerical, only : small

  use fileunits, only : ounit

  use inputlist, only : Wmetrix

  use cputiming, only : Tmetrix

  use allglobal, only : myid, ncpu, cpus, &
                        dBdX, &
                        mn_field, im_field, in_field, mne, ime, ine, &
                        Nt, Nz, Ntz, efmn, ofmn, cfmn, sfmn, &   ! 10 Dec 15;
                        ijreal, &                                ! workspace;
                        sg, guvij, &                             ! calculated in coords;
                        gvuij, &                                 ! this is workspace: nowhere used outside of this routine;
                        goomne, goomno, &
                        gssmne, gssmno, &
                        gstmne, gstmno, &
                        gszmne, gszmno, &
                        gttmne, gttmno, &
                        gtzmne, gtzmno, &
                        gzzmne, gzzmno, &
                        guvijsave

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS

  INTEGER, intent(in) :: lvol, lquad

  INTEGER             :: Lcurvature, ifail, ideriv, jquad

  BEGIN( metrix )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  do jquad = 1, lquad

    gvuij(1:Ntz,0,0) =   one

    gvuij(1:Ntz,1,1) =   guvijsave(1:Ntz,1,1,jquad)
    gvuij(1:Ntz,1,2) =   guvijsave(1:Ntz,1,2,jquad)
    gvuij(1:Ntz,1,3) =   guvijsave(1:Ntz,1,3,jquad)
    gvuij(1:Ntz,2,2) =   guvijsave(1:Ntz,2,2,jquad)
    gvuij(1:Ntz,2,3) =   guvijsave(1:Ntz,2,3,jquad)
    gvuij(1:Ntz,3,3) =   guvijsave(1:Ntz,3,3,jquad)

    ifail = 0
    call tfft( Nt, Nz, gvuij(1:Ntz,0,0), ijreal(1:Ntz) , &
              mne, ime(1:mne), ine(1:mne), goomne(1:mne,jquad), goomno(1:mne,jquad), cfmn(1:mne)    , sfmn(1:mne)    , ifail )
    goomne(0,jquad) = zero ; goomno(0,jquad) = zero

    ifail = 0
    call tfft( Nt, Nz, gvuij(1:Ntz,1,1), gvuij(1:Ntz,1,2), &
              mne, ime(1:mne), ine(1:mne), gssmne(1:mne,jquad), gssmno(1:mne,jquad), gstmne(1:mne,jquad), gstmno(1:mne,jquad), ifail )
    gssmne(0,jquad) = zero ; gssmno(0,jquad) = zero
    gstmne(0,jquad) = zero ; gstmno(0,jquad) = zero

    ifail = 0
    call tfft( Nt, Nz, gvuij(1:Ntz,1,3), gvuij(1:Ntz,2,2), &
              mne, ime(1:mne), ine(1:mne), gszmne(1:mne,jquad), gszmno(1:mne,jquad), gttmne(1:mne,jquad), gttmno(1:mne,jquad), ifail )
    gszmne(0,jquad) = zero ; gszmno(0,jquad) = zero
    gttmne(0,jquad) = zero ; gttmno(0,jquad) = zero

    ifail = 0
    call tfft( Nt, Nz, gvuij(1:Ntz,2,3), gvuij(1:Ntz,3,3), &
              mne, ime(1:mne), ine(1:mne), gtzmne(1:mne,jquad), gtzmno(1:mne,jquad), gzzmne(1:mne,jquad), gzzmno(1:mne,jquad), ifail )
    gtzmne(0,jquad) = zero ; gtzmno(0,jquad) = zero
    gzzmne(0,jquad) = zero ; gzzmno(0,jquad) = zero

  enddo
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN( metrix )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine metrix

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> \brief compute guvijsave
!>
!> @param lquad
!> @param vvol
!> @param ideriv
!> @param Lcurvature
subroutine compute_guvijsave(lquad, vvol, ideriv, Lcurvature)

  use allglobal, only : gaussianabscissae, Ntz, mn_field, guvij, guvijsave, &
                        sg

  implicit none

  INTEGER, intent(in):: vvol, lquad, ideriv, Lcurvature
  INTEGER            :: jquad, ii, jj
  REAL               :: lss

  ! we need to compute guvij and save it in guvijsave
  do jquad = 1, lquad
    lss = gaussianabscissae(jquad,vvol)
    call coords( vvol, lss, Lcurvature, Ntz, mn_field )
    guvijsave(1:Ntz,1:3,1:3,jquad) = guvij(1:Ntz,1:3,1:3,ideriv)
    do ii = 1, 3
      do jj = 1, 3
        guvijsave(1:Ntz,jj,ii,jquad) = guvijsave(1:Ntz,jj,ii,jquad) / sg(1:Ntz, 0)
      enddo
    enddo
  enddo

end subroutine compute_guvijsave
