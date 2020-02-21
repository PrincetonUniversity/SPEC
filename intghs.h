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
                        Mvol, im, in, mne, Ntz, Nt, Nz, &
                        YESstellsym, NOTstellsym, &
                        gaussianweight, gaussianabscissae, &
                        Tsc, Tss, Dtc, Dts, Dzc, Dzs, &
                        Ttc, Tts, Tzc, Tzs, &
                        Lcoordinatesingularity, &
                        pi2pi2nfp, pi2pi2nfphalf, dBdX, &
                        Ntz, NOTstellsym, dBdX, Lsavedguvij, &
                        Ate, Ato, Aze, Azo, &
                        sg, guvij, guvijsave
                        
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in) :: lquad, mn, lvol, lrad, idx
  
  INTEGER             :: jquad, ll, pp, ll1, pp1, uv, ii, jj, io, mn2, lp2, mn2_max, lp2_max, nele, ideriv, ifail, Lcurvature
  INTEGER             :: mi, ni
  
  REAL                :: lss, jthweight, Tl, Dl, ik, sbar

  REAL, allocatable   :: gBupper(:,:), Blower(:,:), basis(:,:,:)

  REAL, allocatable   :: efmn(:), ofmn(:), cfmn(:), sfmn(:)
  REAL, allocatable   :: evmn(:), odmn(:)

  REAL, allocatable   :: ijreal(:), jireal(:), jkreal(:), kjreal(:)
  REAL, allocatable   :: Bloweremn(:,:), Bloweromn(:,:)

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

  if( dBdX%L ) then ; Lcurvature = 3 ; ideriv = 1
  else              ; Lcurvature = 1 ; ideriv = 0
  endif
  
  if (.not. Lsavedguvij) then
    do jquad = 1, lquad
      lss = gaussianabscissae(jquad,lvol)
      WCALL( intghs, coords, ( lvol, lss, Lcurvature, Ntz, mn ) )
      guvijsave(1:Ntz,1:3,1:3,jquad) = guvij(1:Ntz,1:3,1:3,ideriv)
      do ii = 1, 3
        do jj = 1, 3
          guvijsave(1:Ntz,jj,ii,jquad) = guvijsave(1:Ntz,jj,ii,jquad) / sg(1:Ntz, 0)
        enddo
      enddo
    enddo
  endif

  SALLOCATE(basis,     (0:lrad,0:mpol,0:1), zero)
  SALLOCATE(gBupper,   (1:Ntz,3), zero)
  SALLOCATE(Blower,    (1:Ntz,3), zero)
  SALLOCATE(Bloweremn, (1:mn,3), zero)
  SALLOCATE(Bloweromn, (1:mn,3), zero)
  SALLOCATE(efmn,      (1:mn), zero)
  SALLOCATE(ofmn,      (1:mn), zero)
  SALLOCATE(evmn,      (1:mn), zero)
  SALLOCATE(odmn,      (1:mn), zero)
  SALLOCATE(cfmn,      (1:mn), zero)
  SALLOCATE(sfmn,      (1:mn), zero)
  SALLOCATE(ijreal,    (1:mn), zero)
  SALLOCATE(jkreal,    (1:mn), zero)
  SALLOCATE(jireal,    (1:mn), zero)
  SALLOCATE(kjreal,    (1:mn), zero)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!$OMP PARALLEL PRIVATE(jquad,ll1,lss,jthweight,sbar,Tl,Dl,ii,jj,ll,mi,ni,ik,gBupper,Blower,efmn,ofmn,cfmn,sfmn,evmn,odmn,ijreal,jireal,jkreal,kjreal,Bloweremn,Bloweromn,basis) SHARED(lquad,mn,lvol,lrad,idx)
!$OMP DO
  do jquad = 1, lquad
    
    GETTHREAD

    lss = gaussianabscissae(jquad,lvol) ; jthweight = gaussianweight(jquad,lvol)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    if (Lcoordinatesingularity) then
      sbar = (lss + one) * half
      call get_zernike(sbar, lrad, mpol, basis(:,:,0:1))
    else
      call get_cheby(lss, lrad, basis(:,0,0:1))
    endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
    efmn(1:mn) = zero ; ofmn(1:mn) = zero ; cfmn(1:mn) = zero ; sfmn(1:mn) = zero
    evmn(1:mn) = zero ; odmn(1:mn) = zero ;

    ijreal(1:mn) = zero ; jireal(1:mn) = zero ; jkreal(1:mn) = zero ; kjreal(1:mn) = zero
  
    gBupper(:,:) = zero; Blower(:,:) = zero;

    ! A_t  : ijreal, jireal
    ! A_z  : jkreal, kjreal
    ! gB^s : evmn, odmn
    ! gB^t : cfmn, sfmn
    ! gB^z : efmn, ofmn
    
    do ii = 1, mn ; mi = im(ii) ; ni = in(ii) ! loop over Fourier harmonics;
      
      if (Lcoordinatesingularity) then
        do ll = mi, lrad, 2 ! loop over Zernike polynomials; Lrad is the radial resolution;
        ;                      ; efmn(ii) = efmn(ii) + Ate(lvol,idx,ii)%s(ll) * ( basis(ll,mi,1)*half) 
        ;                      ; cfmn(ii) = cfmn(ii) - Aze(lvol,idx,ii)%s(ll) * ( basis(ll,mi,1)*half)
        ;                      ; odmn(ii) = odmn(ii) - Ate(lvol,idx,ii)%s(ll) * ( basis(ll,mi,0)) * ni & 
                                                    - Aze(lvol,idx,ii)%s(ll) * ( basis(ll,mi,0)) * mi
        ;                      ; ijreal(ii) = ijreal(ii) + Ate(lvol,idx,ii)%s(ll) * basis(ll,mi,0)
        ;                      ; jkreal(ii) = jkreal(ii) + Aze(lvol,idx,ii)%s(ll) * basis(ll,mi,0)
        if( NOTstellsym ) then ; ofmn(ii) = ofmn(ii) + Ato(lvol,idx,ii)%s(ll) * ( basis(ll,mi,1)*half)
          ;                    ; sfmn(ii) = sfmn(ii) - Azo(lvol,idx,ii)%s(ll) * ( basis(ll,mi,1)*half)
          ;                    ; evmn(ii) = evmn(ii) + Ato(lvol,idx,ii)%s(ll) * ( basis(ll,mi,1)) * ni & 
                                                    + Azo(lvol,idx,ii)%s(ll) * ( basis(ll,mi,1)) * mi
          ;                    ; jireal(ii) = jireal(ii) + Ato(lvol,idx,ii)%s(ll) * basis(ll,mi,0)
          ;                    ; kjreal(ii) = kjreal(ii) + Azo(lvol,idx,ii)%s(ll) * basis(ll,mi,0)
        endif
        enddo ! end of do ll; 20 Feb 13;
      else
        do ll = 0, lrad ! loop over Chebyshev polynomials; Lrad is the radial resolution;
        ;                      ; efmn(ii) = efmn(ii) + Ate(lvol,idx,ii)%s(ll) * ( basis(ll,0,1))
        ;                      ; cfmn(ii) = cfmn(ii) - Aze(lvol,idx,ii)%s(ll) * ( basis(ll,0,1))
        ;                      ; odmn(ii) = odmn(ii) - Ate(lvol,idx,ii)%s(ll) * ( basis(ll,0,0)) * ni & 
                                                    - Aze(lvol,idx,ii)%s(ll) * ( basis(ll,0,0)) * mi
        ;                      ; ijreal(ii) = ijreal(ii) + Ate(lvol,idx,ii)%s(ll) * basis(ll,0,0)
        ;                      ; jkreal(ii) = jkreal(ii) + Aze(lvol,idx,ii)%s(ll) * basis(ll,0,0)
        if( NOTstellsym ) then ; ofmn(ii) = ofmn(ii) + Ato(lvol,idx,ii)%s(ll) * ( basis(ll,0,1))
          ;                    ; sfmn(ii) = sfmn(ii) - Azo(lvol,idx,ii)%s(ll) * ( basis(ll,0,1))
          ;                    ; evmn(ii) = evmn(ii) + Ato(lvol,idx,ii)%s(ll) * ( basis(ll,0,0)) * ni & 
                                                    + Azo(lvol,idx,ii)%s(ll) * ( basis(ll,0,0)) * mi
          ;                    ; jireal(ii) = jireal(ii) + Ato(lvol,idx,ii)%s(ll) * basis(ll,0,0)
          ;                    ; kjreal(ii) = kjreal(ii) + Azo(lvol,idx,ii)%s(ll) * basis(ll,0,0)
        endif
        enddo ! end of do ll; 20 Feb 13;
      end if ! Lcoordinatesingularity; 01 Jul 19
    enddo ! end of do ii; 20 Feb 13;
    
    call invfft( mn, im, in, efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, gBupper(1:Ntz,3), gBupper(1:Ntz,2) )
    call invfft( mn, im, in, evmn(1:mn), odmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, gBupper(1:Ntz,1), gBupper(1:Ntz,2) )

    do ii = 1, 3
      do jj = 1, 3
        Blower(:,ii) = Blower(:,ii) + gBupper(:,jj) * guvijsave(1:Ntz,jj,ii,jquad)
      enddo
    enddo

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
    ifail = 0
   
    call tfft( Nt, Nz, Blower(1:Ntz,2), Blower(1:Ntz,3), &
               mn, im(1:mn), in(1:mn), Bloweremn(1:mn,2), Bloweromn(1:mn,2), Bloweremn(1:mn,3), Bloweromn(1:mn,3), ifail )
    call tfft( Nt, Nz, Blower(1:Ntz,1), Blower(1:Ntz,3), &
               mn, im(1:mn), in(1:mn), Bloweromn(1:mn,1), Bloweremn(1:mn,1), Bloweremn(1:mn,3), Bloweromn(1:mn,3), ifail )
    Bloweremn(1,1  ) = zero
    Bloweromn(1,2:3) = zero
    
    if (Lcoordinatesingularity) then
      do ii = 1, mn

        if (ii==1) then ;ik = jthweight * two
        else            ;ik = jthweight
        endif

        do ll = im(ii), lrad, 2

          ll1 = (ll - mod(ll,2))/2 ! shrinked dof for Zernike; 02 Jul 19
          
          ! get the basis functions, they are generated in getbco
          Tl = basis(ll, im(ii), 0)
          Dl = basis(ll, im(ii), 1) * half
          
!$OMP ATOMIC UPDATE
          Tss(ll1,ii) = Tss(ll1,ii) + Tl * Bloweremn(ii,1) * ik
!$OMP ATOMIC UPDATE
          Dtc(ll1,ii) = Dtc(ll1,ii) + Dl * Bloweremn(ii,2) * ik
!$OMP ATOMIC UPDATE
          Dzc(ll1,ii) = Dzc(ll1,ii) + Dl * Bloweremn(ii,3) * ik
        
          if (NOTstellsym) then
!$OMP ATOMIC UPDATE
            Tsc(ll1,ii) = Tsc(ll1,ii) + Tl * Bloweromn(ii,1) * ik
!$OMP ATOMIC UPDATE
            Dts(ll1,ii) = Dts(ll1,ii) + Dl * Bloweromn(ii,2) * ik
!$OMP ATOMIC UPDATE
            Dzs(ll1,ii) = Dzs(ll1,ii) + Dl * Bloweromn(ii,3) * ik
          endif

          if (dBdX%L) cycle ! dMD matrix does not depend on geometry
!$OMP ATOMIC UPDATE
          Ttc(ll1,ii) = Ttc(ll1,ii) + (Tl * cfmn(ii) + Dl * jkreal(ii)) * ik
!$OMP ATOMIC UPDATE          
          Tzc(ll1,ii) = Tzc(ll1,ii) + (Tl * efmn(ii) - Dl * ijreal(ii)) * ik

          if (NOTstellsym) then
!$OMP ATOMIC UPDATE          
            Tts(ll1,ii) = Tts(ll1,ii) + (Tl * sfmn(ii) + Dl * kjreal(ii)) * ik
!$OMP ATOMIC UPDATE            
            Tzs(ll1,ii) = Tzs(ll1,ii) + (Tl * ofmn(ii) - Dl * ijreal(ii)) * ik
          endif

        enddo !ll
      enddo !ii
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
    else ! .not.Lcoordinatesingularity; 

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

      do ii = 1, mn

        if (ii==1) then ;ik = jthweight * two
        else            ;ik = jthweight
        endif

        do ll = 0, lrad
          
          ! get the basis functions, they are generated in getbco
          Tl = basis(ll,0,0)
          Dl = basis(ll,0,1)

          ll1 = ll
!$OMP ATOMIC UPDATE
          Tss(ll1,ii) = Tss(ll1,ii) + Tl * Bloweremn(ii,1) * ik
!$OMP ATOMIC UPDATE          
          Dtc(ll1,ii) = Dtc(ll1,ii) + Dl * Bloweremn(ii,2) * ik
!$OMP ATOMIC UPDATE          
          Dzc(ll1,ii) = Dzc(ll1,ii) + Dl * Bloweremn(ii,3) * ik
        
          if (NOTstellsym) then
!$OMP ATOMIC UPDATE          
            Tsc(ll1,ii) = Tsc(ll1,ii) + Tl * Bloweromn(ii,1) * ik
!$OMP ATOMIC UPDATE            
            Dts(ll1,ii) = Dts(ll1,ii) + Dl * Bloweromn(ii,2) * ik
!$OMP ATOMIC UPDATE            
            Dzs(ll1,ii) = Dzs(ll1,ii) + Dl * Bloweromn(ii,3) * ik
          endif

          if (dBdX%L) cycle ! dMD matrix does not depend on geometry
!$OMP ATOMIC UPDATE
          Ttc(ll1,ii) = Ttc(ll1,ii) + (Tl * cfmn(ii) + Dl * jkreal(ii)) * ik
!$OMP ATOMIC UPDATE          
          Tzc(ll1,ii) = Tzc(ll1,ii) + (Tl * efmn(ii) - Dl * ijreal(ii)) * ik

          if (NOTstellsym) then
!$OMP ATOMIC UPDATE          
            Tts(ll1,ii) = Tts(ll1,ii) + (Tl * sfmn(ii) + Dl * kjreal(ii)) * ik
!$OMP ATOMIC UPDATE            
            Tzs(ll1,ii) = Tzs(ll1,ii) + (Tl * ofmn(ii) - Dl * ijreal(ii)) * ik
          endif
          
        enddo !ll
      enddo !ii
      
    endif  ! Lcoordinatesingularity;
    
  enddo !jquad
!$OMP END PARALLEL
  
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
  
  DALLOCATE(basis)
  DALLOCATE(gBupper)
  DALLOCATE(Blower)
  DALLOCATE(Bloweremn)
  DALLOCATE(Bloweromn)
  DALLOCATE(efmn)
  DALLOCATE(ofmn)
  DALLOCATE(evmn)
  DALLOCATE(odmn)
  DALLOCATE(cfmn)
  DALLOCATE(sfmn)
  DALLOCATE(ijreal)
  DALLOCATE(jkreal)
  DALLOCATE(jireal)
  DALLOCATE(kjreal)
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN( intghs )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine intghs

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
