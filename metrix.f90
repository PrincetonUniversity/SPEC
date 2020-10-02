!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (metrics) ! Calculates the metric quantities, $\sqrt g \, g^{\mu\nu}$, which are required for the energy and helicity integrals.

!latex \briefly{Calculates the metric quantities required for volume integrals.}

!latex \calledby{\link{ma00aa}}
!latex \calls{\link{coords}}

!latex \tableofcontents

!latex \subsection{metrics}

!latex \begin{enumerate}
!latex \item The Jacobian, $\sqrt g$, and the `lower' metric elements, $g_{\mu\nu}$, are calculated by \link{coords},
!latex       and are provided on a regular grid in `real-space', i.e. $(\t,\z)$, at a given radial location, i.e. where $s$ is input.
!latex \end{enumerate}

!latex \subsection{plasma region}

!latex \begin{enumerate}
!latex \item In the plasma region, the required terms are $\bar g_{\mu\nu} \equiv g_{\mu\nu}/\sqrt g$.

!latex       \be \begin{array}{ccccccccccccccccccccccccc}
!latex       \sqrt g \; g^{\s\s} & = & \left( g_{\t\t} g_{\z\z} - g_{\t\z} g_{\t\z} \right) / \sqrt g \\
!latex       \sqrt g \; g^{\s\t} & = & \left( g_{\t\z} g_{\s\z} - g_{\s\t} g_{\z\z} \right) / \sqrt g \\
!latex       \sqrt g \; g^{\s\z} & = & \left( g_{\s\t} g_{\t\z} - g_{\t\t} g_{\s\z} \right) / \sqrt g \\
!latex       \sqrt g \; g^{\t\t} & = & \left( g_{\z\z} g_{\s\s} - g_{\s\z} g_{\s\z} \right) / \sqrt g \\
!latex       \sqrt g \; g^{\t\z} & = & \left( g_{\s\z} g_{\s\t} - g_{\t\z} g_{\s\s} \right) / \sqrt g \\
!latex       \sqrt g \; g^{\z\z} & = & \left( g_{\s\s} g_{\t\t} - g_{\s\t} g_{\s\t} \right) / \sqrt g
!latex       \end{array}
!latex       \ee
!latex \end{enumerate}

!latex \subsection{FFTs}

!latex \begin{enumerate}
!latex \item After constructing the required quantities in real space, FFTs provided the required Fourier harmonics, which are returned through global.
!latex       (The ``extended'' Fourier resolution is used.)
!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine metrix( lquad, lvol ) 
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, one
  
  use numerical, only : small
  
  use fileunits, only : ounit
  
  use inputlist, only : Wmetrix
  
  use cputiming, only : Tmetrix
  
  use allglobal, only : myid, ncpu, cpus, &
                        dBdX, &
                        mn, im, in, mne, ime, ine, &
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

subroutine compute_guvijsave(lquad, vvol, ideriv, Lcurvature)

  use allglobal, only : gaussianabscissae, Ntz, mn, guvij, guvijsave, &
                        sg

  implicit none

  INTEGER, intent(in):: vvol, lquad, ideriv, Lcurvature
  INTEGER            :: jquad, ii, jj
  REAL               :: lss

  ! we need to compute guvij and save it in guvijsave
  do jquad = 1, lquad
    lss = gaussianabscissae(jquad,vvol)
    call coords( vvol, lss, Lcurvature, Ntz, mn )
    guvijsave(1:Ntz,1:3,1:3,jquad) = guvij(1:Ntz,1:3,1:3,ideriv)
    do ii = 1, 3
      do jj = 1, 3
        guvijsave(1:Ntz,jj,ii,jquad) = guvijsave(1:Ntz,jj,ii,jquad) / sg(1:Ntz, 0)
      enddo
    enddo
  enddo

end subroutine compute_guvijsave