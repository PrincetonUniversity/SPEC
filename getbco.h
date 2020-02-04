!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (getbco) ! Returns &ldquo;covariant&rdquo; magnetic field, ${\bf B} = B_s \nabla s + B_\theta \nabla \theta + B_\zeta \nabla \zeta$ $.
!latex \briefly{This magnetic field is used to assist the computation of Hessian}

!latex \calledby{\link{dforce}}
!latex \calls{\link{coords}}

!latex \tableofcontents
 
!latex \subsection{covariant representation}

!latex \begin{enumerate}
    
!latex \item \Ais
!latex \item \Bis
!latex \item The covariant representation for the field is ${\bf B} = B_\s \nabla\s + B_\t \nabla\t + B_\z\nabla\z$, where
!latex       \be B_\s &=& B^\s g_{\s\s} + B^\t g_{\s\t} + B^\z g_{\s\z} \nonumber \\
!latex           B_\t &=& B^\s g_{\s\t} + B^\t g_{\t\t} + B^\z g_{\t\z} \nonumber \\
!latex           B_\z &=& B^\s g_{\s\z} + B^\t g_{\t\z} + B^\z g_{\z\z} \nonumber
!latex       \ee
!latex       where $g_{\a\b}$ are the metric elements (computed by \link{coords}).
!latex \item On the interfaces, $B^\s=0$ by construction.

!latex \end{enumerate}
 
!latex \subsection{output data}

!latex Depending on dBdX, the metrics $g_{ij}$ will be computed either in their original value, or their derivative with respect to the interface geometry.

!latex \begin{enumerate}
    
!latex \item The Fourier harmonics of the even-and-odd, covariant components of the magnetic field, $B_\s$, $B_\t$ and $B_\z$, are saved in 
!latex       \begin{itemize}
!latex       \item[] \internal{Bloweremn(1:mn,0:3)},
!latex       \item[] \internal{Bloweromn(1:mn,0:3)},
!latex       \end{itemize}

!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine getbco( lvol, Ntz, lss )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half, one, two, pi2
  
  use numerical, only : vsmall
  
  use fileunits, only : ounit
  
  use inputlist, only : Wgetbco, Mpol, Nvol, Lrad
  
  use cputiming, only : Tgetbco
  
  use allglobal, only : ncpu, myid, cpus, pi2nfp, &
                        Lcoordinatesingularity, Lvacuumregion, Mvol, &
                        NOTstellsym, &
                        Bloweremn, Bloweromn, &
                        mn, im, in, regumm, &
                        efmn, ofmn, cfmn, sfmn, evmn, odmn, comn, simn, &
                        Nt, Nz, &
                        Ate, Aze, Ato, Azo, &
                        cheby, zernike,&
                        sg, guvij, Rij, Zij, &
                        dBdX
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in)  :: lvol, Ntz
  REAL   , intent(in)  :: lss
  REAL                 :: gBupper(1:Ntz,3), Blower(1:Ntz,3)
  
  INTEGER              :: innout, Lcurvature, ii, jj, kk, ll, ifail, ideriv, mi, ni
  REAL                 :: mfactor, sbar
  
  BEGIN(getbco)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( dBdX%L ) then ; Lcurvature = 3 ; ideriv = 1
  else              ; Lcurvature = 1 ; ideriv = 0
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if (Lcoordinatesingularity) then
    sbar = (lss + one) * half
    call get_zernike(sbar, Lrad(lvol), mpol, zernike(:,:,0:1))
  else
    call get_cheby(lss, Lrad(lvol), cheby(0:Lrad(lvol),0:1))
  endif
   
  WCALL( getbco, coords,( lvol, lss, Lcurvature, Ntz, mn ) ) ! get coordinates and derivatives wrt Rj, Zj, at specific radial location;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
  efmn(1:mn) = zero ; ofmn(1:mn) = zero ; cfmn(1:mn) = zero ; sfmn(1:mn) = zero
  evmn(1:mn) = zero ; odmn(1:mn) = zero ; comn(1:mn) = zero ; simn(1:mn) = zero
  gBupper = zero; Blower = zero;

!$OMP PARALLEL DO PRIVATE(ii,mi,ni,ll) SHARED(mn,lvol)   
  do ii = 1, mn ; mi = im(ii) ; ni = in(ii) ! loop over Fourier harmonics;
    
    if (Lcoordinatesingularity) then
      do ll = mi, Lrad(lvol), 2 ! loop over Zernike polynomials; Lrad is the radial resolution;
      ;                      ; efmn(ii) = efmn(ii) + Ate(lvol,0,ii)%s(ll) * ( zernike(ll,mi,1)*half) 
      ;                      ; cfmn(ii) = cfmn(ii) - Aze(lvol,0,ii)%s(ll) * ( zernike(ll,mi,1)*half)
      ;                      ; odmn(ii) = odmn(ii) - Ate(lvol,0,ii)%s(ll) * ( zernike(ll,mi,0)) * ni & 
                                                   - Aze(lvol,0,ii)%s(ll) * ( zernike(ll,mi,0)) * mi
      if( NOTstellsym ) then ; ofmn(ii) = ofmn(ii) + Ato(lvol,0,ii)%s(ll) * ( zernike(ll,mi,1)*half)
        ;                    ; sfmn(ii) = sfmn(ii) - Azo(lvol,0,ii)%s(ll) * ( zernike(ll,mi,1)*half)
        ;                    ; evmn(ii) = evmn(ii) + Ato(lvol,0,ii)%s(ll) * ( zernike(ll,mi,1)) * ni & 
                                                   + Azo(lvol,0,ii)%s(ll) * ( zernike(ll,mi,1)) * mi
      endif
      enddo ! end of do ll; 20 Feb 13;
    else
      do ll = 0, Lrad(lvol) ! loop over Chebyshev polynomials; Lrad is the radial resolution;
      ;                      ; efmn(ii) = efmn(ii) + Ate(lvol,0,ii)%s(ll) * ( cheby(ll,1))
      ;                      ; cfmn(ii) = cfmn(ii) - Aze(lvol,0,ii)%s(ll) * ( cheby(ll,1))
      ;                      ; odmn(ii) = odmn(ii) - Ate(lvol,0,ii)%s(ll) * ( cheby(ll,0)) * ni & 
                                                   - Aze(lvol,0,ii)%s(ll) * ( cheby(ll,0)) * mi
      if( NOTstellsym ) then ; ofmn(ii) = ofmn(ii) + Ato(lvol,0,ii)%s(ll) * ( cheby(ll,1))
        ;                    ; sfmn(ii) = sfmn(ii) - Azo(lvol,0,ii)%s(ll) * ( cheby(ll,1))
        ;                    ; evmn(ii) = evmn(ii) + Ato(lvol,0,ii)%s(ll) * ( cheby(ll,0)) * ni & 
                                                   + Azo(lvol,0,ii)%s(ll) * ( cheby(ll,0)) * mi
      endif
      enddo ! end of do ll; 20 Feb 13;
    end if ! Lcoordinatesingularity; 01 Jul 19
  enddo ! end of do ii; 20 Feb 13;
!$OMP END PARALLEL DO

  call invfft( mn, im, in, efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, gBupper(1:Ntz,3), gBupper(1:Ntz,2) )
  call invfft( mn, im, in, evmn(1:mn), odmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, gBupper(1:Ntz,1), gBupper(1:Ntz,2) )
   
  do ii = 1, 3
    do jj = 1, 3
      Blower(:,ii) = Blower(:,ii) + gBupper(:,jj) * guvij(1:Ntz,jj,ii,ideriv) / sg(1:Ntz,0)
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

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
  RETURN(getbco)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
 end subroutine getbco
 
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
