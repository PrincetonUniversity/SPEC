!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (integrals) ! Calculates volume integrals of Chebyshev-polynomials and covariant field for Hessian computation.

!latex \briefly{Calculates volume integrals of Chebyshev polynomials and covariant field products.}

!latex \calledby{\link{dforce}}
!latex \calls{\link{getbco}}

!latex \tableofcontents

!l tex \newcommand{\bT}[1]{{\overline T}_{#1}}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{chebyshev-metric information}

!latex \begin{enumerate}

!latex \item The following quantities are calculated:

!latex       \be \verb+DToocc(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}  \; \ooint \cos\a_i \cos\a_j                  \\
!latex           \verb+DToocs(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}  \; \ooint \cos\a_i \sin\a_j                  \\
!latex           \verb+DToosc(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}  \; \ooint \sin\a_i \cos\a_j                  \\
!latex           \verb+DTooss(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}  \; \ooint \sin\a_i \sin\a_j                  
!latex       \ee

!latex       \be \verb+TTsscc(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}  \; \bT{p,j}  \; \ooint \cos\a_i \cos\a_j \; \bar g_{\s\s} \\
!latex           \verb+TTsscs(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}  \; \bT{p,j}  \; \ooint \cos\a_i \sin\a_j \; \bar g_{\s\s} \\
!latex           \verb+TTsssc(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}  \; \bT{p,j}  \; \ooint \sin\a_i \cos\a_j \; \bar g_{\s\s} \\
!latex           \verb+TTssss(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}  \; \bT{p,j}  \; \ooint \sin\a_i \sin\a_j \; \bar g_{\s\s}
!latex       \ee

!latex       \be \verb+TDstcc(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}  \; \bT{p,j}' \; \ooint \cos\a_i \cos\a_j \; \bar g_{\s\t} \\
!latex           \verb+TDstcs(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}  \; \bT{p,j}' \; \ooint \cos\a_i \sin\a_j \; \bar g_{\s\t} \\
!latex           \verb+TDsTsc(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}  \; \bT{p,j}' \; \ooint \sin\a_i \cos\a_j \; \bar g_{\s\t} \\
!latex           \verb+TDsTss(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}  \; \bT{p,j}' \; \ooint \sin\a_i \sin\a_j \; \bar g_{\s\t}
!latex       \ee

!latex       \be \verb+TDstcc(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}  \; \bT{p,j}' \; \ooint \cos\a_i \cos\a_j \; \bar g_{\s\z} \\
!latex           \verb+TDstcs(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}  \; \bT{p,j}' \; \ooint \cos\a_i \sin\a_j \; \bar g_{\s\z} \\
!latex           \verb+TDsTsc(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}  \; \bT{p,j}' \; \ooint \sin\a_i \cos\a_j \; \bar g_{\s\z} \\
!latex           \verb+TDsTss(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}  \; \bT{p,j}' \; \ooint \sin\a_i \sin\a_j \; \bar g_{\s\z}
!latex       \ee

!latex       \be \verb+DDstcc(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}' \; \ooint \cos\a_i \cos\a_j \; \bar g_{\t\t} \\
!latex           \verb+DDstcs(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}' \; \ooint \cos\a_i \sin\a_j \; \bar g_{\t\t} \\
!latex           \verb+DDsTsc(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}' \; \ooint \sin\a_i \cos\a_j \; \bar g_{\t\t} \\
!latex           \verb+DDsTss(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}' \; \ooint \sin\a_i \sin\a_j \; \bar g_{\t\t}
!latex       \ee

!latex       \be \verb+DDstcc(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}' \; \ooint \cos\a_i \cos\a_j \; \bar g_{\t\z} \\
!latex           \verb+DDstcs(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}' \; \ooint \cos\a_i \sin\a_j \; \bar g_{\t\z} \\
!latex           \verb+DDsTsc(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}' \; \ooint \sin\a_i \cos\a_j \; \bar g_{\t\z} \\
!latex           \verb+DDsTss(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}' \; \ooint \sin\a_i \sin\a_j \; \bar g_{\t\z}
!latex       \ee

!latex       \be \verb+DDstcc(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}' \; \ooint \cos\a_i \cos\a_j \; \bar g_{\z\z} \\
!latex           \verb+DDstcs(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}' \; \ooint \cos\a_i \sin\a_j \; \bar g_{\z\z} \\
!latex           \verb+DDsTsc(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}' \; \ooint \sin\a_i \cos\a_j \; \bar g_{\z\z} \\
!latex           \verb+DDsTss(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}' \; \ooint \sin\a_i \sin\a_j \; \bar g_{\z\z}
!latex       \ee

!latex       where $\bT{l,i}\equiv T_l \, \bar s^{m_i/2}$ if the domain includes the coordinate singularity, and $\bT{l,i}\equiv T_l$ if not;
!latex       and $\bar g_{\mu\nu} \equiv g_{\mu\nu} / \sqrt g$.

!latex \item The double-angle formulae are used to reduce the above expressions to the Fourier harmonics of $\bar g_{\mu\nu}$:
!latex       see \internal{kija} and \internal{kijs}, which are defined in \link{preset}.

!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine intghs( lquad, mn, lvol, lrad, idx )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half, one, two, pi, pi2
  
  use numerical, only : vsmall, small, sqrtmachprec
  
  use fileunits, only : ounit
  
  use inputlist, only : mpol, Wintghs
  
  use cputiming, only : Tintghs
  
  use allglobal, only : myid, ncpu, cpus, &
                        Mvol, im, in, mne, Ntz, &
                        YESstellsym, NOTstellsym, &
                        gaussianweight, gaussianabscissae, &
                        efmn, ofmn, cfmn, sfmn, evmn, odmn, comn, simn, &
                        Tsc, Tss, Dtc, Dts, Dzc, Dzs, &
                        Ttc, Tts, Tzc, Tzs, &
                        Bloweremn, Bloweromn, &
                        cheby, zernike, &
                        Lcoordinatesingularity, &
                        pi2pi2nfp, pi2pi2nfphalf, dBdX, &
                        ijreal, jireal, jkreal, kjreal
                        
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in) :: lquad, mn, lvol, lrad, idx
  
  INTEGER             :: jquad, ll, pp, ll1, pp1, uv, ii, jj, io, mn2, lp2, mn2_max, lp2_max, nele 
  
  REAL                :: lss, jthweight, Tl, Dl, ik

  BEGIN( intghs )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  Tss = zero
  Dtc = zero
  Dzc = zero

  if (.not.dBdX%L) then
    Ttc = zero
    Tzc = zero
  endif

  if (NOTstellsym) then
    Tsc = zero
    Dts = zero
    Dzs = zero
    if (.not.dBdX%L) then
      Tts = zero
      Tzs = zero
    endif
    
  endif !NOTstellsym

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  if (Lcoordinatesingularity) then

    do jquad = 1, lquad

      lss = gaussianabscissae(jquad,lvol) ; jthweight = gaussianweight(jquad,lvol)
      WCALL( intghs, getbco, (lvol, Ntz, lss, jquad, idx) )

!$OMP PARALLEL DO PRIVATE(ii,ik,ll,ll1,Tl,Dl) SHARED(mn,jthweight,lrad)
      do ii = 1, mn

        if (ii==1) then ;ik = jthweight * two
        else            ;ik = jthweight
        endif

        do ll = im(ii), lrad, 2

          ll1 = (ll - mod(ll,2))/2 ! shrinked dof for Zernike; 02 Jul 19
          
          ! get the basis functions, they are generated in getbco
          Tl = zernike(ll, im(ii), 0)
          Dl = zernike(ll, im(ii), 1) * half

          Tss(ll1,ii) = Tss(ll1,ii) + Tl * Bloweremn(ii,1) * ik
          Dtc(ll1,ii) = Dtc(ll1,ii) + Dl * Bloweremn(ii,2) * ik
          Dzc(ll1,ii) = Dzc(ll1,ii) + Dl * Bloweremn(ii,3) * ik
        
          if (NOTstellsym) then
            Tsc(ll1,ii) = Tsc(ll1,ii) + Tl * Bloweromn(ii,1) * ik
            Dts(ll1,ii) = Dts(ll1,ii) + Dl * Bloweromn(ii,2) * ik
            Dzs(ll1,ii) = Dzs(ll1,ii) + Dl * Bloweromn(ii,3) * ik
          endif

          if (dBdX%L) cycle ! dMD matrix does not depend on geometry

          Ttc(ll1,ii) = Ttc(ll1,ii) + (Tl * cfmn(ii) + Dl * jkreal(ii)) * ik
          Tzc(ll1,ii) = Tzc(ll1,ii) + (Tl * efmn(ii) - Dl * ijreal(ii)) * ik

          if (NOTstellsym) then
            Tts(ll1,ii) = Tts(ll1,ii) + (Tl * sfmn(ii) + Dl * kjreal(ii)) * ik
            Tzs(ll1,ii) = Tzs(ll1,ii) + (Tl * ofmn(ii) - Dl * ijreal(ii)) * ik
          endif

        enddo !ll
      enddo !ii
!$OMP END PARALLEL DO
    enddo !jquad
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  else ! .not.Lcoordinatesingularity; 

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
    do jquad = 1, lquad

      lss = gaussianabscissae(jquad,lvol) ; jthweight = gaussianweight(jquad,lvol)
      WCALL( intghs, getbco, (lvol, Ntz, lss, jquad, idx) )
!$OMP PARALLEL DO PRIVATE(ii,ik,ll,ll1,Tl,Dl) SHARED(mn,jthweight,lrad)
      do ii = 1, mn

        if (ii==1) then ;ik = jthweight * two
        else            ;ik = jthweight
        endif

        do ll = 0, lrad
          
          ! get the basis functions, they are generated in getbco
          Tl = cheby(ll, 0)
          Dl = cheby(ll, 1)

          ll1 = ll

          Tss(ll1,ii) = Tss(ll1,ii) + Tl * Bloweremn(ii,1) * ik
          Dtc(ll1,ii) = Dtc(ll1,ii) + Dl * Bloweremn(ii,2) * ik
          Dzc(ll1,ii) = Dzc(ll1,ii) + Dl * Bloweremn(ii,3) * ik
        
          if (NOTstellsym) then
            Tsc(ll1,ii) = Tsc(ll1,ii) + Tl * Bloweromn(ii,1) * ik
            Dts(ll1,ii) = Dts(ll1,ii) + Dl * Bloweromn(ii,2) * ik
            Dzs(ll1,ii) = Dzs(ll1,ii) + Dl * Bloweromn(ii,3) * ik
          endif

          if (dBdX%L) cycle ! dMD matrix does not depend on geometry

          Ttc(ll1,ii) = Ttc(ll1,ii) + (Tl * cfmn(ii) + Dl * jkreal(ii)) * ik
          Tzc(ll1,ii) = Tzc(ll1,ii) + (Tl * efmn(ii) - Dl * ijreal(ii)) * ik

          if (NOTstellsym) then
            Tts(ll1,ii) = Tts(ll1,ii) + (Tl * sfmn(ii) + Dl * kjreal(ii)) * ik
            Tzs(ll1,ii) = Tzs(ll1,ii) + (Tl * ofmn(ii) - Dl * ijreal(ii)) * ik
          endif
          
        enddo !ll
      enddo !ii
!OMP END PARALLEL DO
    enddo !jquad

  endif  ! Lcoordinatesingularity;
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  Tss = Tss * pi2pi2nfphalf
  Dtc = Dtc * pi2pi2nfphalf
  Dzc = Dzc * pi2pi2nfphalf

  if (.not.dBdX%L) then
    Ttc = Ttc * pi2pi2nfphalf
    Tzc = Tzc * pi2pi2nfphalf
  endif

  if (NOTstellsym) then
    
    Tsc = Tsc * pi2pi2nfphalf
    Dts = Dts * pi2pi2nfphalf
    Dzc = Dzc * pi2pi2nfphalf

    if (.not.dBdX%L) then
      Tts = Tts * pi2pi2nfphalf
      Tzs = Tzs * pi2pi2nfphalf
    endif
  
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN( intghs )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine intghs

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
