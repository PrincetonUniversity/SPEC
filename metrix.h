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
!latex       (The `extended' Fourier resolution is used.)
!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine metrix( lvol, lss ) 
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, one
  
  use numerical, only : small
  
  use fileunits, only : ounit
  
  use inputlist, only : Wmetrix
  
  use cputiming, only : Tmetrix
  
  use allglobal, only : myid, ncpu, cpus, &
                        dBdX, &
                        mn, im, in, mne, ime, ine, &
                        isr, Nt, Nz, Ntz, trigm, trign, trigwk, efmn, ofmn, cfmn, sfmn, & ! 10 Dec 15;
                        ijreal, &                                                         ! workspace;
                        sg, guvij, &                                                      ! calculated in coords;
                        gvuij, &                                                          ! this is workspace: nowhere used outside of this routine;
                        goomne, goomno, &
                        gssmne, gssmno, &
                        gstmne, gstmno, &
                        gszmne, gszmno, &
                        gttmne, gttmno, &
                        gtzmne, gtzmno, &
                        gzzmne, gzzmno
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in) :: lvol
  REAL   , intent(in) :: lss
  
  INTEGER             :: Lcurvature, ifail, ideriv
  
  BEGIN( metrix )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  FATAL( metrix, abs(lss).gt.one, invalid local radial coordinate )
#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( dBdX%L ) then ; Lcurvature = 3 ; ideriv = 1
  else              ; Lcurvature = 1 ; ideriv = 0
  endif
  
  WCALL( metrix, coords, ( lvol, lss, Lcurvature, Ntz, mn ) ) ! this returns guvij \equiv g_{\mu\nu}; 17 Apr 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  ijreal(1:Ntz) = one / sg(1:Ntz,0) ! SRH; 01 Aug 17;
  
  gvuij(1:Ntz,0,0) =   guvij(1:Ntz,0,0,ideriv)                 ! required for helicity calculation; 17 Dec 15;
  gvuij(1:Ntz,1,1) =   guvij(1:Ntz,1,1,ideriv) * ijreal(1:Ntz) ! 10 Dec 15;
  gvuij(1:Ntz,1,2) =   guvij(1:Ntz,1,2,ideriv) * ijreal(1:Ntz)
  gvuij(1:Ntz,1,3) =   guvij(1:Ntz,1,3,ideriv) * ijreal(1:Ntz)
  gvuij(1:Ntz,2,2) =   guvij(1:Ntz,2,2,ideriv) * ijreal(1:Ntz)
  gvuij(1:Ntz,2,3) =   guvij(1:Ntz,2,3,ideriv) * ijreal(1:Ntz)
  gvuij(1:Ntz,3,3) =   guvij(1:Ntz,3,3,ideriv) * ijreal(1:Ntz)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  ijreal(1:Ntz) = zero

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  ifail = 0
  call tfft( Nt, Nz, gvuij(1:Ntz,0,0), ijreal(1:Ntz)   , isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz), &
             mne, ime(1:mne), ine(1:mne), goomne(1:mne), goomno(1:mne), cfmn(1:mne)    , sfmn(1:mne)    , ifail )
  goomne(0) = zero ; goomno(0) = zero

  ifail = 0
  call tfft( Nt, Nz, gvuij(1:Ntz,1,1), gvuij(1:Ntz,1,2), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz), &
             mne, ime(1:mne), ine(1:mne), gssmne(1:mne), gssmno(1:mne), gstmne(1:mne), gstmno(1:mne), ifail )
  gssmne(0) = zero ; gssmno(0) = zero
  gstmne(0) = zero ; gstmno(0) = zero

  ifail = 0
  call tfft( Nt, Nz, gvuij(1:Ntz,1,3), gvuij(1:Ntz,2,2), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz), &
             mne, ime(1:mne), ine(1:mne), gszmne(1:mne), gszmno(1:mne), gttmne(1:mne), gttmno(1:mne), ifail )
  gszmne(0) = zero ; gszmno(0) = zero
  gttmne(0) = zero ; gttmno(0) = zero

  ifail = 0
  call tfft( Nt, Nz, gvuij(1:Ntz,2,3), gvuij(1:Ntz,3,3), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz), &
             mne, ime(1:mne), ine(1:mne), gtzmne(1:mne), gtzmno(1:mne), gzzmne(1:mne), gzzmno(1:mne), ifail )
  gtzmne(0) = zero ; gtzmno(0) = zero
  gzzmne(0) = zero ; gzzmno(0) = zero
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN( metrix )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine metrix

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
