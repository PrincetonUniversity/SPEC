!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Computes covariant components of the magnetic field on interfaces.

!latex \end{enumerate} \subsection{covariant representation} \begin{enumerate}
    
!latex \item From ${\bf B}\equiv\nabla\times\left[A_\t\nabla\t+A_\z\nabla\z\right]$, the magnetic field is
!latex       \be {\bf B} & = & 
!latex           (\sqrt g)^{-1}(\partial_\t A_\z - \partial_\z A_\t ) {\bf e}_\s - (\sqrt g)^{-1}\partial_\s A_\z {\bf e}_\t + (\sqrt g)^{-1}\partial_\s A_\t {\bf e}_\z \\ 
!latex                   & = & B^s {\bf e}_\s + B^\t {\bf e}_\t + B^\z{\bf e}_\z.
!latex       \ee

!latex \item The covariant representation for the field is ${\bf B} = B_\s \nabla\s + B_\t \nabla\t + B_\z\nabla\z$, where
!latex       \be B_\s &=& B^\s g_{\s\s} +  B^\t g_{\s\t} +  B^\z g_{\s\z} \nonumber \\
!latex           B_\t &=& B^\s g_{\s\t} +  B^\t g_{\t\t} +  B^\z g_{\t\z} \nonumber \\
!latex           B_\z &=& B^\s g_{\s\z} +  B^\t g_{\t\z} +  B^\z g_{\z\z} \nonumber
!latex       \ee
!latex       where $g_{\a\b}$ are the metric elements (computed by \verb+co01aa+).

!latex \item On the interfaces, $B^\s=0$ by construction.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine sc00aa( lvol, Ntz )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half, one, two, pi2
  
  use numerical, only : vsmall
  
  use fileunits, only : ounit
  
  use inputlist, only : Wsc00aa, Nvol, Lrad
  
  use cputiming, only : Tsc00aa
  
  use allglobal, only : ncpu, myid, cpus, pi2nfp, &
                        Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, Mvol, &
                        NOTstellsym, &
                        Bsubtemn, Bsubzemn, Bsubtomn, Bsubzomn, &
                        mn, im, in, halfmm, &
                        efmn, ofmn, cfmn, sfmn, evmn, odmn, comn, simn, &
                        trigm, trign, trigwk, isr, Nt, Nz, &
                        Ate, Aze, Ato, Azo, &
                        TTll, &
                        sg, guvij, Rij, Zij
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in)  :: lvol, Ntz
  REAL                 :: dAt(1:Ntz), dAz(1:Ntz), dBs(1:Ntz), Bt(1:Ntz), Bz(1:Ntz)
  
  INTEGER              :: innout, Lcurvature, ii, jj, kk, ll, ifail, ideriv, mi, ni
  REAL                 :: lss, mfactor
  
  BEGIN(sc00aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  FATALMESS(sc00aa, lvol.lt.1 .or. lvol.gt.Mvol, illegal lvol)
#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  ideriv = 0

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  do innout = 0, 1
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   if( lvol.eq.   1 .and. innout.eq.0 ) cycle
   if( lvol.eq.Mvol .and. innout.eq.1 ) cycle

   if( ( Lcoordinatesingularity .and. innout.eq.0 ) .or. ( Lvacuumregion .and. innout.eq.1 ) ) cycle

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   lss = two * innout - one ! recall that innout is effective local radial coordinate; 24 Apr 13;
   
   Lcurvature = 1 ! this will instruct sc00aa to compute the metric elements; 20 Jun 14;
   
   WCALL(sc00aa, co01aa,( lvol, lss, Lcurvature, Ntz, mn )) ! get coordinates and derivatives wrt Rj, Zj, at specific radial location;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   if( Lplasmaregion ) then
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
    efmn(1:mn) = zero ; ofmn(1:mn) = zero ; cfmn(1:mn) = zero ; sfmn(1:mn) = zero
    evmn(1:mn) = zero ; odmn(1:mn) = zero ; comn(1:mn) = zero ; simn(1:mn) = zero
    
    do ii = 1, mn ; mi = im(ii) ; ni = in(ii)! loop over Fourier harmonics; 13 Sep 13;
     
     if( Lcoordinatesingularity ) then ; mfactor = halfmm(ii) * half ! regularity factor; 01 Jul 14;
     else                              ; mfactor = zero
     endif
     
     do ll = 0, Lrad(lvol) ! loop over Chebyshev polynomials; Lrad is the radial resolution;
      ;                      ; efmn(ii) = efmn(ii) + Ate(lvol,ideriv,ii)%s(ll) * ( TTll(ll,innout,1) + mfactor ) ! B^\t; 01 Jul 14;
      ;                      ; cfmn(ii) = cfmn(ii) + Aze(lvol,ideriv,ii)%s(ll) * ( TTll(ll,innout,1) + mfactor ) ! B^\z; 01 Jul 14;
      ;                      ; odmn(ii) = odmn(ii) + Ate(lvol,ideriv,ii)%s(ll) * ( TTll(ll,innout,1) + mfactor ) * ni & 
                                                   + Aze(lvol,ideriv,ii)%s(ll) * ( TTll(ll,innout,1) + mfactor ) * mi
      if( NOTstellsym ) then ; ofmn(ii) = ofmn(ii) + Ato(lvol,ideriv,ii)%s(ll) * ( TTll(ll,innout,1) + mfactor )
       ;                     ; sfmn(ii) = sfmn(ii) + Azo(lvol,ideriv,ii)%s(ll) * ( TTll(ll,innout,1) + mfactor )
      ;                      ; evmn(ii) = evmn(ii) - Ato(lvol,ideriv,ii)%s(ll) * ( TTll(ll,innout,1) + mfactor ) * ni & 
                                                   - Azo(lvol,ideriv,ii)%s(ll) * ( TTll(ll,innout,1) + mfactor ) * mi
     !else                   ; ofmn(ii) = zero
     ! ;                     ; sfmn(ii) = zero
      endif
     enddo ! end of do ll; 20 Feb 13;
    
    enddo ! end of do ii; 20 Feb 13;
    
    call invfft( mn, im, in, efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, dAt(1:Ntz), dAz(1:Ntz), isr, trigm, trign, trigwk ) ! map to real space; 20 Feb 13;
    
    Bt(1:Ntz) = ( - dAz(1:Ntz) * guvij(1:Ntz,2,2,0) + dAt(1:Ntz) * guvij(1:Ntz,2,3,0) ) / sg(1:Ntz,0)
    Bz(1:Ntz) = ( - dAz(1:Ntz) * guvij(1:Ntz,2,3,0) + dAt(1:Ntz) * guvij(1:Ntz,3,3,0) ) / sg(1:Ntz,0)

!   call invfft( mn, im, in, evmn(1:mn), odmn(1:mn), comn(1:mn), simn(1:mn), Nt, Nz, dBs(1:Ntz), dAz(1:Ntz), isr, trigm, trign, trigwk ) ! map to real space; 20 Feb 13;
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
   else ! matches if( Lplasmavolume ) ; this is the vacuum region ; 24 Apr 13;
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
!latex \end{enumerate} \subsection{exterior vacuum volume} \begin{enumerate}
   
!latex \item With ${\bf B}=\nabla \Phi = I \nabla \t + G \nabla \z + \phi_\s \nabla \s + \phi_\t \nabla \t + \phi_\z \nabla \z$,
!latex       we have for the contravariant components of the magnetic field
!latex       \be B^\s \equiv {\bf B} \cdot \nabla \s & = & I \;g^{\t\s} + G \;g^{\z\s} + \phi_\s \;g^{\s\s} + \phi_\t \;g^{\t\s} + \phi_\z \;g^{\z\s} \\
!latex           B^\t \equiv {\bf B} \cdot \nabla \t & = & I \;g^{\t\t} + G \;g^{\z\t} + \phi_\s \;g^{\s\t} + \phi_\t \;g^{\t\t} + \phi_\z \;g^{\z\t} \\
!latex           B^\z \equiv {\bf B} \cdot \nabla \z & = & I \;g^{\t\z} + G \;g^{\z\z} + \phi_\s \;g^{\s\z} + \phi_\t \;g^{\t\z} + \phi_\z \;g^{\z\z}
!latex       \ee
!latex       and for the covariant components 
!latex       \be B_\s \equiv {\bf B} \cdot {\bf e}_\s & = & \phi_\s     \\
!latex           B_\t \equiv {\bf B} \cdot {\bf e}_\t & = & \phi_\t + I \\
!latex           B_\z \equiv {\bf B} \cdot {\bf e}_\z & = & \phi_\z + G
!latex       \ee

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
    FATALMESS(sc00aa, .true., vacuum region under construction)
   
!   efmn(1:mn) = zero ; ofmn(1:mn) = zero ; cfmn(1:mn) = zero ; sfmn(1:mn) = zero ! initialize summation; 24 Apr 13; 
!   evmn(1:mn) = zero ; odmn(1:mn) = zero ; comn(1:mn) = zero ; simn(1:mn) = zero
!   
!   do ii = 1, mn
!    
!    do ll = 0, Lrad(lvol)
!    !;                      ; comn(ii) = 
!     ;                      ; odmn(ii) = odmn(ii) +          Ate(lvol,ideriv,ii)%s(ll) *   TTll(ll,innout,1) ! the o-e convention is inverted for scalar potential;
!     ;                      ; efmn(ii) = efmn(ii) + im(ii) * Ate(lvol,ideriv,ii)%s(ll) *   TTll(ll,innout,0)
!     ;                      ; cfmn(ii) = cfmn(ii) - in(ii) * Ate(lvol,ideriv,ii)%s(ll) *   TTll(ll,innout,0)
!     if( NOTstellsym ) then
!     !;                     ; simn(ii) = 
!      ;                     ; evmn(ii) = evmn(ii) +          Ato(lvol,ideriv,ii)%s(ll) *   TTll(ll,innout,1) ! the o-e convention is inverted for scalar potential;
!      ;                     ; ofmn(ii) = ofmn(ii) - im(ii) * Ato(lvol,ideriv,ii)%s(ll) *   TTll(ll,innout,0)
!      ;                     ; sfmn(ii) = sfmn(ii) + in(ii) * Ato(lvol,ideriv,ii)%s(ll) *   TTll(ll,innout,0)
!     else                   
!      ;                     ; simn(ii) = zero
!      ;                     ; evmn(ii) = zero
!      ;                     ; ofmn(ii) = zero
!      ;                     ; sfmn(ii) = zero
!     endif
!    enddo ! end of do ll; 20 Feb 13;
!    
!   enddo ! end of do ii; 20 Feb 13;
!   
!   call invfft( mn, im, in, comn(1:mn), simn(1:mn), evmn(1:mn), odmn(1:mn), Nt, Nz, jireal(1:Ntz), dAs(1:Ntz), isr, trigm, trign, trigwk ) ! map to real space; 20 Feb 13;
!   call invfft( mn, im, in, efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz,    dAt(1:Ntz), dAz(1:Ntz), isr, trigm, trign, trigwk ) ! map to real space; 20 Feb 13;
   

!   do ii = 1, 3
!    
!    select case( ii )
!    case( 1 ) ; jj = 2 ; kk = 3
!    case( 2 ) ; jj = 3 ; kk = 1
!    case( 3 ) ; jj = 1 ; kk = 2
!    end select
!    
!    lgvuij(1:Ntz,ii,1) = ( guvij(1:Ntz,jj,2,0)*guvij(1:Ntz,kk,3,0) - guvij(1:Ntz,jj,3,0)*guvij(1:Ntz,kk,2,0) ) / sg(1:Ntz,0)**2 ! "raising" metric elements; 24 Apr 13;
!    lgvuij(1:Ntz,ii,2) = ( guvij(1:Ntz,jj,3,0)*guvij(1:Ntz,kk,1,0) - guvij(1:Ntz,jj,1,0)*guvij(1:Ntz,kk,3,0) ) / sg(1:Ntz,0)**2
!    lgvuij(1:Ntz,ii,3) = ( guvij(1:Ntz,jj,1,0)*guvij(1:Ntz,kk,2,0) - guvij(1:Ntz,jj,2,0)*guvij(1:Ntz,kk,1,0) ) / sg(1:Ntz,0)**2
!    
!   enddo ! end of do ii; 13 May 13;
!   
!   ijimag(1:Ntz) = curtor*lgvuij(1:Ntz,2,1) + curpol*lgvuij(1:Ntz,3,1) + dAs(1:Ntz)*lgvuij(1:Ntz,1,1) + dAt(1:Ntz)*lgvuij(1:Ntz,2,1) + dAz(1:Ntz)*lgvuij(1:Ntz,3,1) ! B^s;
!   jireal(1:Ntz) = curtor*lgvuij(1:Ntz,2,2) + curpol*lgvuij(1:Ntz,3,2) + dAs(1:Ntz)*lgvuij(1:Ntz,1,2) + dAt(1:Ntz)*lgvuij(1:Ntz,2,2) + dAz(1:Ntz)*lgvuij(1:Ntz,3,2) ! B^t;
!   jiimag(1:Ntz) = curtor*lgvuij(1:Ntz,2,3) + curpol*lgvuij(1:Ntz,3,3) + dAs(1:Ntz)*lgvuij(1:Ntz,1,3) + dAt(1:Ntz)*lgvuij(1:Ntz,2,3) + dAz(1:Ntz)*lgvuij(1:Ntz,3,3) ! B^z;
!   
!   dBB(1:Ntz) = jireal(1:Ntz) * ( curtor + dAt(1:Ntz) ) + jiimag(1:Ntz) * ( curpol + dAz(1:Ntz) )                              ! NOT INCLUDING RADIAL FIELD; 17 Apr 13;
!  !dBB(1:Ntz) = jireal(1:Ntz) * ( curtor + dAt(1:Ntz) ) + jiimag(1:Ntz) * ( curpol + dAz(1:Ntz) ) + ijimag(1:Ntz) * dAs(1:Ntz) !     INCLUDING RADIAL FIELD; 17 Apr 13;
!   
!   ijreal(1:Ntz) =                                                   half * dBB(1:Ntz) ! p + B^2/2; 13 Sep 13; RECALL THAT THIS IS THE VACUUM REGION; volume is not defined;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   endif ! end of if( Lplasmaregion ) ; 24 Apr 13;
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   ifail = 0
  !call tfft( Nt, Nz, Bt(1:Ntz), Bz(1:Ntz), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz), & ! compute force-imbalance and spectral constraints; 20 Feb 13;
  !           mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), ifail )
   call tfft( Nt, Nz, Bt(1:Ntz), Bz(1:Ntz), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz), & ! compute force-imbalance and spectral constraints; 20 Feb 13;
              mn, im(1:mn), in(1:mn), Bsubtemn(1:mn,innout,lvol), Bsubtomn(1:mn,innout,lvol), Bsubzemn(1:mn,innout,lvol), Bsubzomn(1:mn,innout,lvol), ifail )
   
   if( Wsc00aa ) then
    write(ounit,'("sc00aa : "10x " : myid=",i3," ; lvol=",i3," ; innout="i2" ; dB^smn=",99es18.10)') myid, lvol, innout, odmn(1:mn)
    write(ounit,'("sc00aa : "10x " : myid=",i3," ; lvol=",i3," ; innout="i2" ;  B_tmn=",99es18.10)') myid, lvol, innout, Bsubtemn(1:mn,innout,lvol)
    write(ounit,'("sc00aa : "10x " : myid=",i3," ; lvol=",i3," ; innout="i2" ;  B_zmn=",99es18.10)') myid, lvol, innout, Bsubzemn(1:mn,innout,lvol)
    if( NOTstellsym ) then
     FATALMESS(sc00aa, .true., under construction)
    endif
   endif
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
  enddo ! end of do innout = 0, 1 ; 20 Jun 14;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(sc00aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
 end subroutine sc00aa
 
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
