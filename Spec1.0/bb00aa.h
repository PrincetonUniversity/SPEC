!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Computes B.B and the spectral condensation constraints on the interfaces.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine bb00aa( lvol, iocons, ideriv, Ntz, dAt, dAz, XX, YY, length, DDl, MMl, iflag ) !dRR, dZZ, iflag )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half, one, two, pi2
  
  use numerical, only : vsmall
  
  use fileunits, only : ounit
  
  use inputlist, only : Wbb00aa, Igeometry, Nvol, Ntor, Lrad, gamma, pscale, adiabatic, curtor, curpol
  
  use cputiming, only : Tbb00aa
  
  use allglobal, only : ncpu, myid, cpus, pi2nfp, &
                        Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, Mvol, &
                        iRbc, iZbs, iRbs, iZbc, &
!                       dnjn, &
                        YESstellsym, NOTstellsym, &
                        mn, im, in, halfmm, &
                        ijreal, ijimag, jireal, jiimag, &
                        efmn, ofmn, cfmn, sfmn, evmn, odmn, comn, simn, &
                        trigm, trign, trigwk, isr, Nt, Nz, &
                        Ate, Aze, Ato, Azo, &
                        TTll, &
                        sg, guvij, iRij, iZij, dRij, dZij, tRij, tZij, &
                        mpnq, & !nq, & ! 18 Jul 14;
                        Bemn, Bomn, Iomn, Iemn, Somn, Semn, &
                        Pomn, Pemn, &
                        vvolume
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in)  :: lvol, iocons, ideriv, Ntz, iflag
  REAL                 :: dAt(1:Ntz), dAz(1:Ntz), XX(1:Ntz), YY(1:Ntz), dRR(1:Ntz,-1:1), dZZ(1:Ntz,-1:1), DDl, MMl

  REAL                 :: IIl(1:Ntz), length(1:Ntz), dLL(1:Ntz), dPP(1:Ntz,0:1)
  
  INTEGER              :: Lcurvature, ii, jj, kk, ll, ifail, ivol, lnn!, oicons
  REAL                 :: dBB(1:Ntz), lss, mfactor
  
  REAL                 :: dAs(1:Ntz)!, dRdt(-1:1,0:1), dZdt(-1:1,0:1)
  REAL                 :: lgvuij(1:Ntz,1:3,1:3) ! local workspace; 13 Sep 13;
  
  BEGIN(bb00aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  FATALMESS(bb00aa, lvol.lt.1 .or. lvol.gt.Mvol, illegal lvol)
  FATALMESS(bb00aa, lvol.eq.1 .and. iocons.eq.0, illegal combination)
  FATALMESS(bb00aa, lvol.eq.Mvol .and. iocons.eq.1, illegal combination)
  FATALMESS(bb00aa, iflag.lt.0 .or. iflag.gt.1, illegal iflag)
#endif
 
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  dAt(1:Ntz) = zero ; dAz(1:Ntz) = zero ! initialize intent out; 01 Jul 14;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  lss = two * iocons - one ! recall that iocons is effective local radial coordinate; 24 Apr 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  Lcurvature = 1

  WCALL(bb00aa, co01aa,( lvol, lss, Lcurvature, Ntz, mn )) ! get coordinates and derivatives wrt Rj, Zj, at specific radial location;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Lplasmaregion ) then
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
!latex \end{enumerate} \subsection{interior ``Beltrami'' volume} \begin{enumerate}

!latex \item The field strength is given by $B^2 = B^\s B_\s + B^\t B_\t + B^\z B_\z$, and on the interfaces $B^\s=0$ by construction.

!latex \item The magnetic field is
!latex       $\sqrt g \; {\bf B} = (\partial_\t A_\z - \partial_\z A_\t ) {\bf e}_\s - \partial_\s A_\z {\bf e}_\t + \partial_\s A_\t {\bf e}_\z$.

!latex \item The covariant components of the field are computed via $B_\t = B^\t g_{\t\t} + B^\z g_{\t\z}$ and $B_\z = B^\t g_{\t\z} + B^\z g_{\z\z}$.
  
!latex \item The expression for $B^2$ is
!latex       \be
!latex       (\sqrt g)^2 B^2 = A_\z^\prime \; A_\z^\prime \; g_{\t\t} - 2 \; A_\z^\prime \; A_\t^\prime \; g_{\t\z} + A_\t^\prime \; A_\t^\prime \; g_{\z\z},
!latex       \ee
!latex       where the ``$\prime$'' denotes derivative with respect to $\s$.
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   efmn(1:mn) = zero ; sfmn(1:mn) = zero ; cfmn(1:mn) = zero ; ofmn(1:mn) = zero
   
   do ii = 1, mn ! loop over Fourier harmonics; 13 Sep 13;
    
    if( Lcoordinatesingularity ) then ; mfactor = halfmm(ii) * half ! only required at outer interface, where \bar s = 1; 15 Jan 15;
    else                              ; mfactor = zero
    endif
    
    do ll = 0, Lrad(lvol) ! loop over Chebyshev polynomials; Lrad is the radial resolution;
     ;                      ; efmn(ii) = efmn(ii) +          Ate(lvol,ideriv,ii)%s(ll) * ( TTll(ll,iocons,1) + mfactor ) ! ideriv labels deriv. wrt mu, pflux; 
     ;                      ; cfmn(ii) = cfmn(ii) +          Aze(lvol,ideriv,ii)%s(ll) * ( TTll(ll,iocons,1) + mfactor )
     if( NOTstellsym ) then ; ofmn(ii) = ofmn(ii) +          Ato(lvol,ideriv,ii)%s(ll) * ( TTll(ll,iocons,1) + mfactor )
      ;                     ; sfmn(ii) = sfmn(ii) +          Azo(lvol,ideriv,ii)%s(ll) * ( TTll(ll,iocons,1) + mfactor )
     endif
    enddo ! end of do ll; 20 Feb 13;
    
   enddo ! end of do ii; 20 Feb 13;

   call invfft( mn, im, in, efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, dAt(1:Ntz), dAz(1:Ntz), isr, trigm, trign, trigwk ) ! map to real space;
   
   dBB(1:Ntz) = half * & ! included factor of half here to be consistent with fc02aa; 01 Jul 14;
 ( dAz(1:Ntz   )*dAz(1:Ntz   )*guvij(1:Ntz,2,2,0) - two * dAz(1:Ntz   )*dAt(1:Ntz   )*guvij(1:Ntz,2,3,0) + dAt(1:Ntz   )*dAt(1:Ntz   )*guvij(1:Ntz,3,3,0) ) &
 / sg(1:Ntz,0)**2
   
   ijreal(1:Ntz) = adiabatic(lvol) * pscale / vvolume(lvol)**gamma + dBB(1:Ntz) ! p + B^2/2; 13 Sep 13;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
  else ! matches if( Lplasmavolume ) ; this is the vacuum region ; 24 Apr 13;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
!latex \end{enumerate} \subsection{exterior vacuum volume} \begin{enumerate}
   
!latex \item With ${\bf B}=\nabla \Phi = I \nabla \t + G \nabla \z + \phi_\s \nabla \s + \phi_\t \nabla \t + \phi_\z \nabla \z$,
!latex       we have for the contravariant components of the magnetic field
!latex       \be \begin{array}{ccccccccccc}
!latex           B^\s \equiv {\bf B} \cdot \nabla \s & = & I \;g^{\t\s} + G \;g^{\z\s} + \phi_\s \;g^{\s\s} + \phi_\t \;g^{\t\s} + \phi_\z \;g^{\z\s} \\
!latex           B^\t \equiv {\bf B} \cdot \nabla \t & = & I \;g^{\t\t} + G \;g^{\z\t} + \phi_\s \;g^{\s\t} + \phi_\t \;g^{\t\t} + \phi_\z \;g^{\z\t} \\
!latex           B^\z \equiv {\bf B} \cdot \nabla \z & = & I \;g^{\t\z} + G \;g^{\z\z} + \phi_\s \;g^{\s\z} + \phi_\t \;g^{\t\z} + \phi_\z \;g^{\z\z}
!latex       \end{array} \ee
!latex       and for the covariant components 
!latex       \be \begin{array}{cccclcccccc}
!latex           B_\s & \equiv & {\bf B} \cdot {\bf e}_\s & = & \phi_\s     \\
!latex           B_\t & \equiv & {\bf B} \cdot {\bf e}_\t & = & \phi_\t + I \\
!latex           B_\z & \equiv & {\bf B} \cdot {\bf e}_\z & = & \phi_\z + G
!latex       \end{array} \ee

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   efmn(1:mn) = zero ; ofmn(1:mn) = zero ; cfmn(1:mn) = zero ; sfmn(1:mn) = zero ! initialize summation; 24 Apr 13; 
   evmn(1:mn) = zero ; odmn(1:mn) = zero ; comn(1:mn) = zero ; simn(1:mn) = zero
   
   do ii = 1, mn
    
    do ll = 0, Lrad(lvol)
!   !;                      ; comn(ii) = 
     ;                      ; odmn(ii) = odmn(ii) +          Ate(lvol,ideriv,ii)%s(ll) *   TTll(ll,iocons,1) ! the o-e convention is inverted for scalar;
     ;                      ; efmn(ii) = efmn(ii) + im(ii) * Ate(lvol,ideriv,ii)%s(ll) *   TTll(ll,iocons,0)
     ;                      ; cfmn(ii) = cfmn(ii) - in(ii) * Ate(lvol,ideriv,ii)%s(ll) *   TTll(ll,iocons,0)
     if( NOTstellsym ) then
!    !;                     ; simn(ii) = 
      ;                     ; evmn(ii) = evmn(ii) +          Ato(lvol,ideriv,ii)%s(ll) *   TTll(ll,iocons,1) ! the o-e convention is inverted for scalar;
      ;                     ; ofmn(ii) = ofmn(ii) - im(ii) * Ato(lvol,ideriv,ii)%s(ll) *   TTll(ll,iocons,0)
      ;                     ; sfmn(ii) = sfmn(ii) + in(ii) * Ato(lvol,ideriv,ii)%s(ll) *   TTll(ll,iocons,0)
     else                   
      ;                     ; simn(ii) = zero
      ;                     ; evmn(ii) = zero
      ;                     ; ofmn(ii) = zero
      ;                     ; sfmn(ii) = zero
     endif
    enddo ! end of do ll; 20 Feb 13;
    
   enddo ! end of do ii; 20 Feb 13;
   
   call invfft( mn, im, in, comn(1:mn), simn(1:mn), evmn(1:mn), odmn(1:mn), Nt, Nz, jireal(1:Ntz), dAs(1:Ntz), isr, trigm, trign, trigwk ) ! map to real space;
   call invfft( mn, im, in, efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz,    dAt(1:Ntz), dAz(1:Ntz), isr, trigm, trign, trigwk ) ! map to real space;
   
   
   do ii = 1, 3
    
    select case( ii )
    case( 1 ) ; jj = 2 ; kk = 3
    case( 2 ) ; jj = 3 ; kk = 1
    case( 3 ) ; jj = 1 ; kk = 2
    end select
    
    lgvuij(1:Ntz,ii,1) = ( guvij(1:Ntz,jj,2,0)*guvij(1:Ntz,kk,3,0) - guvij(1:Ntz,jj,3,0)*guvij(1:Ntz,kk,2,0) ) / sg(1:Ntz,0)**2 ! "raising" metric elements;
    lgvuij(1:Ntz,ii,2) = ( guvij(1:Ntz,jj,3,0)*guvij(1:Ntz,kk,1,0) - guvij(1:Ntz,jj,1,0)*guvij(1:Ntz,kk,3,0) ) / sg(1:Ntz,0)**2
    lgvuij(1:Ntz,ii,3) = ( guvij(1:Ntz,jj,1,0)*guvij(1:Ntz,kk,2,0) - guvij(1:Ntz,jj,2,0)*guvij(1:Ntz,kk,1,0) ) / sg(1:Ntz,0)**2
    
   enddo ! end of do ii; 13 May 13;
   
! !ijimag(1:Ntz) = curtor*lgvuij(1:Ntz,2,1) + curpol*lgvuij(1:Ntz,3,1) + dAs(1:Ntz)*lgvuij(1:Ntz,1,1) + dAt(1:Ntz)*lgvuij(1:Ntz,2,1) + dAz(1:Ntz)*lgvuij(1:Ntz,3,1) ! B^s;
   jireal(1:Ntz) = curtor*lgvuij(1:Ntz,2,2) + curpol*lgvuij(1:Ntz,3,2) + dAs(1:Ntz)*lgvuij(1:Ntz,1,2) + dAt(1:Ntz)*lgvuij(1:Ntz,2,2) + dAz(1:Ntz)*lgvuij(1:Ntz,3,2) ! B^t;
   jiimag(1:Ntz) = curtor*lgvuij(1:Ntz,2,3) + curpol*lgvuij(1:Ntz,3,3) + dAs(1:Ntz)*lgvuij(1:Ntz,1,3) + dAt(1:Ntz)*lgvuij(1:Ntz,2,3) + dAz(1:Ntz)*lgvuij(1:Ntz,3,3) ! B^z;
   
   dBB(1:Ntz) = half * ( jireal(1:Ntz) * ( curtor + dAt(1:Ntz) ) + jiimag(1:Ntz) * ( curpol + dAz(1:Ntz) ) )                              ! NOT including radial field;
! !dBB(1:Ntz) = half * ( jireal(1:Ntz) * ( curtor + dAt(1:Ntz) ) + jiimag(1:Ntz) * ( curpol + dAz(1:Ntz) ) + ijimag(1:Ntz) * dAs(1:Ntz) ) !     including radial field;
   
   ijreal(1:Ntz) =                                                          dBB(1:Ntz) ! p + B^2/2; 13 Sep 13; RECALL THAT THIS IS THE VACUUM REGION; volume is not defined;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
  endif ! end of if( Lplasmaregion ) ; 24 Apr 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( iflag .eq. 1 ) goto 9999 ! iflag = 1 indicates the derivatives of the force are to be calculated; derivatives of magnetic field calculated above;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Igeometry.ge.3 ) then ! include additional spectral constraints; 08 Nov 13;
   
   do ivol = 0, 1
    
    call invfft( mn, im(1:mn), in(1:mn),            iRbc(1:mn,lvol-1+ivol),              iRbs(1:mn,lvol-1+ivol), &
                                                    iZbc(1:mn,lvol-1+ivol),              iZbs(1:mn,lvol-1+ivol), & 
                                         Nt, Nz, iRij(1:Ntz,lvol-1+ivol), iZij(1:Ntz,lvol-1+ivol), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz) )

    call invfft( mn, im(1:mn), in(1:mn), im(1:mn) * iRbs(1:mn,lvol-1+ivol), - im(1:mn) * iRbc(1:mn,lvol-1+ivol), &
                                         im(1:mn) * iZbs(1:mn,lvol-1+ivol), - im(1:mn) * iZbc(1:mn,lvol-1+ivol), &
                                         Nt, Nz, tRij(1:Ntz,lvol-1+ivol), tZij(1:Ntz,lvol-1+ivol), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz) )
   enddo ! end of do ivol = 0, 1 ; 18 Jul 14;
   
   dRij(1:Ntz,lvol) = iRij(1:Ntz,lvol) - iRij(1:Ntz,lvol-1)
   dZij(1:Ntz,lvol) = iZij(1:Ntz,lvol) - iZij(1:Ntz,lvol-1)
   
   length(1:Ntz) = sqrt( dRij(1:Ntz,lvol)**2 + dZij(1:Ntz,lvol)**2 )
   
! THE FOLLOWING IS INCORRECT; the following two lines cannot possibly be identical; there is a sign factor that is mysteriously absent; 04 Dec 14;
   
   if( iocons.eq.0 ) then ; dLL(1:Ntz) = + ( dRij(1:Ntz,lvol) * tRij(1:Ntz,lvol-1+iocons) + dZij(1:Ntz,lvol) * tZij(1:Ntz,lvol-1+iocons) ) / length(1:Ntz)
   else                   ; dLL(1:Ntz) = + ( dRij(1:Ntz,lvol) * tRij(1:Ntz,lvol-1+iocons) + dZij(1:Ntz,lvol) * tZij(1:Ntz,lvol-1+iocons) ) / length(1:Ntz)
   endif

! THE   ABOVE   IS INCORRECT; the   above   two lines cannot possibly be identical; there is a sign factor that is mysteriously absent; 04 Dec 14;
   
   if( iocons.eq.1 ) then ! include spectral condensation constraints; local to interface, i.e. no tri-diagonal structure;
    ;                      ; efmn(1:mn) = ( mpnq(1:mn)            ) * iRbc(1:mn,lvol)
    ;                      ; sfmn(1:mn) = ( mpnq(1:mn)            ) * iZbs(1:mn,lvol)
    if( NOTstellsym ) then ; ofmn(1:mn) = ( mpnq(1:mn)            ) * iRbs(1:mn,lvol)
     ;                     ; cfmn(1:mn) = ( mpnq(1:mn)            ) * iZbc(1:mn,lvol)
    else                   ; ofmn(1:mn) = zero
     ;                     ; cfmn(1:mn) = zero
    endif ! end of if( NOTstellsym ) ; 20 Feb 13;
    
    call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), &
 Nt, Nz, XX(1:Ntz), YY(1:Ntz), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz) )
    
    if( YESstellsym ) then ; DDl = sum(              ( iRbc(1:mn,lvol)**2 + iZbs(1:mn,lvol)**2                                           ) )
     ;                     ; MMl = sum( mpnq(1:mn) * ( iRbc(1:mn,lvol)**2 + iZbs(1:mn,lvol)**2                                           ) ) / DDl
    else                   ; DDl = sum(              ( iRbc(1:mn,lvol)**2 + iZbs(1:mn,lvol)**2 + iRbs(1:mn,lvol)**2 + iZbc(1:mn,lvol)**2 ) )
     ;                     ; MMl = sum( mpnq(1:mn) * ( iRbc(1:mn,lvol)**2 + iZbs(1:mn,lvol)**2 + iRbs(1:mn,lvol)**2 + iZbc(1:mn,lvol)**2 ) ) / DDl
    endif
    
    IIl(1:Ntz) = tRij(1:Ntz,lvol) * ( XX(1:Ntz) - MMl * iRij(1:Ntz,lvol) ) &
               + tZij(1:Ntz,lvol) * ( YY(1:Ntz) - MMl * iZij(1:Ntz,lvol) ) 
    
   else ! matches if( iocons.eq.1 ) ; 11 Aug 14;

    IIl(1:Ntz) = zero ! placeholder; 11 Aug 14;

   endif ! end of if( iocons.eq.1 ) ; 20 Feb 13;
   
  else ! matches if( Igeometry.ge.3 ) ; 11 Aug 14;
   
   IIl(1:Ntz) = zero ! placeholder; 11 Aug 14;

  endif ! end of if( Igeometry.ge.3 ) ; 18 Jul 14;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \end{enumerate} \subsection{spectral constraints} \begin{enumerate}

!latex \item In addition to the physical-force-balance constraints, namely that $[[p+B^2/2]]=0$ across the interfaces,
!latex       additional angle constraints are required to obtain a unique Fourier representation of the interface geometry.

!latex \item Introducing the angle functional: a weighted combination of the ``polar'' constraint;
!latex       the normalized, poloidal, spectral width \cite{Hirshman_Meier_85,Hirshman_Breslau_98};
!latex       the poloidal-angle origin constraint;
!latex       and the ``length'' of the angle curves
!latex       \be F \!\equiv\! \sum_{i=1}^{N-1} \alpha_i \ooint \!\! \frac{1}{\Theta_{i,\t}}
!latex          +     \sum_{i=1}^{N-1} \beta_i \, M_i
!latex          +     \sum_{i=1}^{N-1} \gamma_i \int_{0}^{2\pi} \! \frac{1}{2}\left[Z_i(0,\z)-Z_{i,0}\right]^2 d\z
!latex          +     \ooint \!\! \sum_{i=1}^{N} \delta_i \, L_i
!latex       \ee
!latex       where $i$ labels the interfaces, and
!latex       \be \Theta_{i,\t} & \equiv & \frac{x \, y_\t - x_\t \, y}{x^2+y^2}, \\
!latex       M_i          & \equiv & \frac{\sum_j m_j^p ( R_{j,i}^2 + Z_{j,i}^2 )}{\sum_j       ( R_{j,i}^2 + Z_{j,i}^2 )} ,\\
!latex       L_i          & \equiv & \sqrt{      [R_{i}(\t,\z)-R_{i-1}(\t,\z)]^2+[Z_{i}(\t,\z)-Z_{i-1}(\t,\z)]^2},
!latex       \ee
!latex       and where $j$ labels the Fourier harmonics.
!latex       The $\alpha_i$, $\beta_i$, $\gamma_i$ and $\delta_i$ are user-supplied weight factors.

!latex \item The polar constraint is derived from defining $\tan \Theta \equiv y/x$, where
!latex       \be x(\t,\z) &  \equiv &  R_{i}(\t,\z)-R_{i,0}(\z), \\
!latex           y(\t,\z) &  \equiv &  Z_{i}(\t,\z)-Z_{i,0}(\z),
!latex       \ee
!latex       and where the geometric center of each interface is given by the arc-length weighted integrals,
!latex       \be
!latex           R_{i,0} & \equiv & \int_0^{2\pi} \!\!\!\! d\t \; R_i(\t,\z) \sqrt{R_{i,\t}(\t,\z)^2+Z_{i,\t}(\t,\z)^2}, \\
!latex           Z_{i,0} & \equiv & \int_0^{2\pi} \!\!\!\! d\t \; Z_i(\t,\z) \sqrt{R_{i,\t}(\t,\z)^2+Z_{i,\t}(\t,\z)^2},
!latex       \ee
!latex       and $\cos \Theta = x / \sqrt{x^2+y^2}$ has been used to simplify the expressions and to avoid divide-by-zero.

!latex \item Only tangential variations will be allowed to find the extremum of $F$, 
!latex       \be 
!latex           \delta R_i(\t,\z) & \equiv & R_{i,\t}(\t,\z) \, \delta u(\t,\z),\\
!latex           \delta Z_i(\t,\z) & \equiv & Z_{i,\t}(\t,\z) \, \delta u(\t,\z),
!latex       \ee
!latex       from which it follows that the variation in each Fourier harmonic is
!latex       \be
!latex           \delta R_{j,i} = \ds \ooint R_{i,\t}(\t,\z) \, \delta u(\t,\z) \, \cos(m_j\t-n_j\z), \\
!latex           \delta Z_{j,i} = \ds \ooint Z_{i,\t}(\t,\z) \, \delta u(\t,\z) \, \sin(m_j\t-n_j\z),
!latex       \ee
!latex       and
!latex       \be 
!latex           \delta R_{i,\t}(\t,\z) & \equiv & R_{i,\t\t}(\t,\z) \, \delta u(\t,\z) + R_{i,\t}(\t,\z) \, \delta u_\t (\t,\z) \\
!latex           \delta Z_{i,\t}(\t,\z) & \equiv & Z_{i,\t\t}(\t,\z) \, \delta u(\t,\z) + Z_{i,\t}(\t,\z) \, \delta u_\t (\t,\z)
!latex       \ee

!latex \item The variation in $F$ is 
!latex       \be \delta F 
!latex       & = & \sum_{i=1}^{N-1} \alpha_i \;\;\;\ooint \!\! \left( \frac{-2\Theta_{i,\t\t}}{\Theta_{i,\t}^2} \right) \delta u_i \nonumber \\
!latex       & + & \sum_{i=1}^{N-1} \beta_i  \;\;\;\ooint \!\! \left( R_{i,\t} X_i + Z_{i,\t} Y_i \right) \delta u_i \nonumber \\
!latex       & + & \sum_{i=1}^{N-1} \gamma_i  \;\;\; \int \! d\z \;   \left( Z_{i}(0,\z)-Z_{i,0} \right) Z_{i,\t} \; \delta u_i \nonumber \\
!latex       & + & \sum_{i=1}^{N-1} \delta_i \;\;\;\ooint \!\! \left( \frac{\Delta R_{i  } R_{i,\t} + \Delta Z_{i  } Z_{i,\t}}{L_{i  }} \right) \delta u_i
!latex       \nonumber \\
!latex       & + & \sum_{i=1}^{N-1} \delta_{i+1}   \ooint \!\! \left( \frac{\Delta R_{i+1} R_{i,\t} + \Delta Z_{i+1} Z_{i,\t}}{L_{i+1}} \right) \delta u_i
!latex       \label{eq:firstvariation}
!latex       \ee
!latex       where, for the stellarator symmetric case,
!latex       \be
!latex           X_i & \equiv & \ds \sum\nolimits_j ( m_j^p - M_i ) \, R_{j,i} \cos(m_j\t-n_j\z), \\
!latex           Y_i & \equiv & \ds \sum\nolimits_j ( m_j^p - M_i ) \, Z_{j,i} \sin(m_j\t-n_j\z),
!latex       \ee
!latex       and 
!latex       \be \Delta R_{i} & \equiv & R_{i}(\t,\z)-R_{i-1}(\t,\z),\\
!latex           \Delta Z_{i} & \equiv & Z_{i}(\t,\z)-Z_{i-1}(\t,\z),
!latex       \ee

!latex \item The spectral constraints derived from \Eqn{firstvariation} are
!latex       \be I_i & \equiv & - 2 \alpha_i \frac{\Theta_{i,\t\t}}{\Theta_{i,\t}^{2}}
!latex               + \beta_i \left( R_{i,\t} X_i + Z_{i,\t} Y_i \right)
!latex               + \gamma_i \left( Z_{i}(0,\z) - Z_{i,0} \right) Z_{i,\t}(0,\z) \nonumber \\
!latex             & + & \delta_{i  } \frac{\Delta R_{i  } R_{i,\t} + \Delta Z_{i  } Z_{i,\t}}{L_i}
!latex               -   \delta_{i+1} \frac{\Delta R_{i+1} R_{i,\t} + \Delta Z_{i+1} Z_{i,\t}}{L_{i+1}}
!latex       \label{eq:spectralconstraints}
!latex       \ee

!l!tex \item Note that choosing $p=2$ gives $X=-R_{\t\t}$ and $Y=-Z_{\t\t}$, and the spectrally condensed angle constraint, $R_\t X + Z_\t Y=0$, 
!l!tex       becomes $\partial_\t (R_\t^2+Z_\t^2)=0$, 
!l!tex       which defines the equal arc length angle.

!latex \item The poloidal-angle origin term, namely $\gamma_i \left( Z_{i}(0,\z) - Z_{i,0} \right) Z_{i,\t}(0,\z)$
!latex       is only used to constrain the $m_j=0$ harmonics.

!latex \item The construction of the angle functional was influenced by the following considerations:
!latex (i)   The minimal spectral width constraint is very desirable as it reduces the required Fourier resolution, 
!latex       but it does not constrain the $m=0$ harmonics
!latex       and the minimizing spectral-width poloidal-angle may not be consistent with the poloidal angle used on adjacent interfaces.
!latex (ii)  The regularization of the vector potential and the coordinate interpolation near the coordinate origin (see elsewhere)
!latex       assumes that the poloidal angle is the polar angle.
!latex (iii) The user will provide the Fourier harmonics of the boundary, and thus the user will implicitly define the poloidal angle used on the boundary.
!latex (iv)  Minimizing the length term will ensure that the poloidal angle used on each interface
!latex       is smoothly connected to the poloidal angle used on adjacent interfaces.

!latex \item A suitable choice of the weight factors, $\alpha_i$, $\beta_i$, $\gamma_i$ and $\delta_i$, will ensure
!latex       that the polar constraint dominates for the innermost surfaces and that this constraint rapidly becomes insignificant away from the origin;
!latex       that the minimal spectral constraint dominates in the ``middle'';
!latex       and that the minimizing length constraint will be significant near the origin and dominant near the edge,
!latex       so that the minimizing spectral width angle will be continuously connected to the polar angle on the innermost surfaces
!latex       and the user-implied angle at the plasma boundary.
!latex       The length constraint should not be insignificant where the spectral constraint is dominant (so that the $m=0$ harmonics are constrained).

!latex \item The polar constraint does not need normalization.
!latex       The spectral width constraint has already been normalized.
!latex       The length constraint is not yet normalized, but perhaps it should be.

!latex \item The spectral constraints given in \Eqn{spectralconstraints} need to be differentiated
!latex       with respect to the interface Fourier harmonics, $R_{j,i}$ and $Z_{j,i}$.
!latex       The first and second terms lead to a block diagonal hessian, and the length term leads to a block tri-diagonal hessian.

!latex \item Including the poloidal-angle origin constraint means that the polar angle constraint can probably be ignored, i.e. $\alpha_i=0$.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! if( iflag .eq.0 ) goto 9999 ! will not corrupt Bemn or Bomn; the spectral terms are not required; 08 Nov 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! XX(1:Ntz) = zero ; YY(1:Ntz) = zero ; IIl(1:Ntz) = zero ; dLL(1:Ntz) = zero ; dRR(1:Ntz,-1:1) = zero ; dZZ(1:Ntz,-1:1) = zero

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! if( Igeometry.ge.3 ) then ! include additional spectral constraints; 08 Nov 13;
   
!  do ivol = 0, 1
    
!   call invfft( mn, im(1:mn), in(1:mn),            iRbc(1:mn,lvol-1+ivol),              iRbs(1:mn,lvol-1+ivol), &
!                                                   iZbc(1:mn,lvol-1+ivol),              iZbs(1:mn,lvol-1+ivol), & 
!                                                   Nt, Nz, iRij(1:Ntz,lvol-1+ivol), iZij(1:Ntz,lvol-1+ivol), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz) )

!   call invfft( mn, im(1:mn), in(1:mn), im(1:mn) * iRbs(1:mn,lvol-1+ivol), - im(1:mn) * iRbc(1:mn,lvol-1+ivol), &
!                                        im(1:mn) * iZbs(1:mn,lvol-1+ivol), - im(1:mn) * iZbc(1:mn,lvol-1+ivol), &
!                                                   Nt, Nz, dRij(1:Ntz,lvol-1+ivol), dZij(1:Ntz,lvol-1+ivol), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz) )
!  enddo

!  if( iocons.eq.1 ) then ! include spectral condensation constraints; local to interface, i.e. no tri-diagonal structure;
    
!   call invfft( mn, im(1:mn), in(1:mn), im(1:mn) * iRbs(1:mn,lvol), - im(1:mn) * iRbc(1:mn,lvol), im(1:mn) * iZbs(1:mn,lvol), - im(1:mn) * iZbc(1:mn,lvol), &
!                Nt, Nz, dRij(1:Ntz,lvol), dZij(1:Ntz,lvol), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz) )
     
!   ;                      ; efmn(1:mn) = ( mp(1:mn) + nq(1:mn) ) * iRbc(1:mn,lvol)
!   ;                      ; sfmn(1:mn) = ( mp(1:mn) + nq(1:mn) ) * iZbs(1:mn,lvol)
!   if( NOTstellsym ) then ; ofmn(1:mn) = ( mp(1:mn) + nq(1:mn) ) * iRbs(1:mn,lvol)
!    ;                     ; cfmn(1:mn) = ( mp(1:mn) + nq(1:mn) ) * iZbc(1:mn,lvol)
!   else                   ; ofmn(1:mn) = zero
!    ;                     ; cfmn(1:mn) = zero
!   endif ! end of if( NOTstellsym ) ; 20 Feb 13;
    
!   call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, XX(1:Ntz), YY(1:Ntz), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz) )
    
!   IIl(1:Ntz) = tRij(1:Ntz,lvol) * XX(1:Ntz) + tZij(1:Ntz,lvol) * YY(1:Ntz) ! ( poloidal ) spectral minimization condition; 20 Feb 13;
    
!  endif ! end of if( iocons.eq.1 ) ; 20 Feb 13;
      
!  call invfft( mn, im(1:mn), in(1:mn), iRbc(1:mn,lvol-1), iRbs(1:mn,lvol-1), iZbc(1:mn,lvol+0), iZbs(1:mn,lvol+0), & 
!               Nt, Nz, iRij(1:Ntz,lvol-1), iZij(1:Ntz,lvol), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz) )
!  
!  call invfft( mn, im(1:mn), in(1:mn), im(1:mn) * iRbs(1:mn,lvol-1), - im(1:mn) * iRbc(1:mn,lvol-1), im(1:mn) * iZbs(1:mn,lvol-1), - im(1:mn) * iZbc(1:mn,lvol-1), &
!               Nt, Nz, dRij(1:Ntz,lvol-1), dZij(1:Ntz,lvol-1), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz) )

!  ;efmn(1:mn) = iRbc(1:mn,lvol) - iRbc(1:mn,lvol-1)
!  ;sfmn(1:mn) = iZbs(1:mn,lvol) - iZbs(1:mn,lvol-1)
!  if( NOTstellsym ) then
!   ofmn(1:mn) = iRbs(1:mn,lvol) - iRbs(1:mn,lvol-1)
!   cfmn(1:mn) = iZbc(1:mn,lvol) - iZbc(1:mn,lvol-1)
!  else
!   ofmn(1:mn) = zero
!   cfmn(1:mn) = zero
!  endif
   
!  call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), & 
!               Nt, Nz, dRR(1:Ntz,0), dZZ(1:Ntz,0), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz) )
   
!  ivol = lvol - 1 + iocons
   
!  ;ofmn(1:mn) = - im(1:mn)*iRbc(1:mn,ivol)
!  ;cfmn(1:mn) = + im(1:mn)*iZbs(1:mn,ivol)
!  if( NOTstellsym ) then
!   efmn(1:mn) = + im(1:mn)*iRbs(1:mn,ivol)
!   sfmn(1:mn) = - im(1:mn)*iZbc(1:mn,ivol)
!  else  
!   efmn(1:mn) =   zero
!   sfmn(1:mn) =   zero
!  endif
   
!  call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), & 
!               Nt, Nz, dRR(1:Ntz,1), dZZ(1:Ntz,1), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz) )

!  dLL(1:Ntz) = ( dRR(1:Ntz,0) * dRR(1:Ntz,1) + dZZ(1:Ntz,0) * dZZ(1:Ntz,1) ) / sqrt( dRR(1:Ntz,0)**2 + dZZ(1:Ntz,0)**2 ) ! 18 Jul 14;

!  dLL(1:Ntz) = ( (iRij(1:Ntz,lvol)-iRij(1:Ntz,lvol-1)) * tRij(1:Ntz,lvol-1+iocons) &
!               + (iZij(1:Ntz,lvol)-iZij(1:Ntz,lvol-1)) * tZij(1:Ntz,lvol-1+iocons) ) &
!               / sqrt( (iRij(1:Ntz,lvol)-iRij(1:Ntz,lvol-1))**2 + (iZij(1:Ntz,lvol)-iZij(1:Ntz,lvol-1))**2 )

!   if( lvol.lt.Mvol ) then ! include constraint on poloidal angle origin; 18 Jul 14;
!    
!    efmn(1:mn) = zero ; sfmn(1:mn) = zero ; ofmn(1:mn) = zero ; cfmn(1:mn) = zero
!    
!    do ivol = 0, 1
!     
!     do lnn = 0, Ntor
!      
!      ;efmn(1+lnn) = sum( ( iRbc(1:mn,lvol-0+ivol) - iRbc(1:mn,lvol-1+ivol) ) * abs(dnjn(1:mn,lnn)) )
!      ;sfmn(1+lnn) = sum( ( iZbs(1:mn,lvol-0+ivol) - iZbs(1:mn,lvol-1+ivol) ) *     dnjn(1:mn,lnn)  )
!      if( NOTstellsym ) then
!       ofmn(1+lnn) = sum( ( iRbs(1:mn,lvol-0+ivol) - iRbs(1:mn,lvol-1+ivol) ) *     dnjn(1:mn,lnn)  )
!       cfmn(1+lnn) = sum( ( iZbc(1:mn,lvol-0+ivol) - iZbc(1:mn,lvol-1+ivol) ) * abs(dnjn(1:mn,lnn)) )
!      else
!       ofmn(1+lnn) = zero
!       cfmn(1+lnn) = zero
!      endif
!      
!     enddo ! end of do lnn; 18 Jul 14;
!
!     call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), & 
!                Nt, Nz, dRR(1:Ntz,ivol), dZZ(1:Ntz,ivol), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz) )
!     
!    enddo ! end of do ivol = 0, 1 ; 18 Jul 14;
!
!    dPP(1:Ntz,0) = dZZ(1:Ntz,1) * dRR(1:Ntz,0) - dZZ(1:Ntz,0) * dRR(1:Ntz,1) ! should use different variable to avoid confusion; 18 Jul 14;
!
!   endif ! end of if( lvol.lt.Mvol ) ; 18 Jul 14;

   
!   efmn(1:mn) = zero ; sfmn(1:mn) = zero ; ofmn(1:mn) = zero ; cfmn(1:mn) = zero
!   
!   do lnn = 0, Ntor
!    
!    ;efmn(1+lnn) = sum( ( iRbc(1:mn,lvol) - iRbc(1:mn,lvol-1) ) * abs(dnjn(1:mn,lnn)) )
!    ;sfmn(1+lnn) = sum( ( iZbs(1:mn,lvol) - iZbs(1:mn,lvol-1) ) *     dnjn(1:mn,lnn)  )
!    if( NOTstellsym ) then
!     ofmn(1+lnn) = sum( ( iRbs(1:mn,lvol) - iRbs(1:mn,lvol-1) ) *     dnjn(1:mn,lnn)  )
!     cfmn(1+lnn) = sum( ( iZbc(1:mn,lvol) - iZbc(1:mn,lvol-1) ) * abs(dnjn(1:mn,lnn)) )
!    else
!     ofmn(1+lnn) = zero
!     cfmn(1+lnn) = zero
!    endif
!    
!   enddo ! end of do lnn; 18 Jul 14;
!
!   call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), & 
!                Nt, Nz, dRR(1:Ntz,0), dZZ(1:Ntz,0), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz) )
!
!   do ivol = 0, 1
!   
!    efmn(1:mn) = zero ; sfmn(1:mn) = zero ; ofmn(1:mn) = zero ; cfmn(1:mn) = zero
!    
!    do lnn = 0, Ntor
!     
!     ;ofmn(1+lnn) = sum( - im(1:mn) * iRbc(1:mn,lvol-1+ivol) *     dnjn(1:mn,lnn)  )
!     ;cfmn(1+lnn) = sum( + im(1:mn) * iZbs(1:mn,lvol-1+ivol) * abs(dnjn(1:mn,lnn)) )
!     if( NOTstellsym ) then
!      efmn(1+lnn) = sum( + im(1:mn) * iRbs(1:mn,lvol-1+ivol) * abs(dnjn(1:mn,lnn)) )
!      sfmn(1+lnn) = sum( - im(1:mn) * iZbc(1:mn,lvol-1+ivol) *     dnjn(1:mn,lnn)  )
!     else
!      efmn(1+lnn) = zero
!      sfmn(1+lnn) = zero
!     endif
!     
!    enddo ! end of do lnn; 18 Jul 14;
!    
!    call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), & 
!                 Nt, Nz, dRR(1:Ntz,2*ivol-1), dZZ(1:Ntz,2*ivol-1), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz) )
!    
!   enddo ! end of do ivol = 0, 1 ; 18 Jul 14;
!
!   dPP(1:Ntz,1) = ( dRR(1:Ntz,0) * dRR(1:Ntz,2*iocons-1) + dZZ(1:Ntz,0) * dZZ(1:Ntz,2*iocons-1) ) / sqrt( dRR(1:Ntz,0)**2 + dZZ(1:Ntz,0)**2 )
!
!  endif ! end of if( Igeometry.ge.3 ) ;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!  ijreal(1:Ntz) contains the pressure + magnetic energy term;
  
#ifdef DEBUG
  FATALMESS(bb00aa, iocons.lt.0 .or. iocons.gt.2, error)
#endif
  
  ;ifail = 0
  ;call tfft( Nt, Nz, ijreal(1:Ntz), IIl(1:Ntz  ), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz), & ! compute force-imbalance and spectral constraints;
              mn, im(1:mn), in(1:mn), Bemn(1:mn,lvol,iocons), Bomn(1:mn,lvol,iocons), Iemn(1:mn,lvol       ), Iomn(1:mn,lvol       ), ifail )
  
  if( Igeometry.ge.3 ) then ! add minimal length constraint; 18 Jul 14;
   
!  ifail = 0
!  call tfft( Nt, Nz, dLL(1:Ntz)   , dPP(1:Ntz,0), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz), & !                             spectral constraints; 18 Jul 14;
!             mn, im(1:mn), in(1:mn), Semn(1:mn,lvol,iocons), Somn(1:mn,lvol,iocons), Pemn(1:mn,lvol,2), Pomn(1:mn,lvol,2), ifail )

   ifail = 0 ; ijimag(1:Ntz) = zero

   call tfft( Nt, Nz, dLL(1:Ntz), ijimag(1:Ntz), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz), & !                             spectral constraints; 18 Jul 14;
              mn, im(1:mn), in(1:mn), Semn(1:mn,lvol,iocons), Somn(1:mn,lvol,iocons), Pemn(1:mn,lvol,iocons), Pomn(1:mn,lvol,iocons), ifail )

#ifdef DEBUG
   if( Wbb00aa ) then
    write(ounit,'("bb00aa : ", 10x ," : lvol=",i3," ; iocons="i2" ; Somn="999es13.5)') lvol, iocons, Somn(1:mn,lvol,iocons)
    write(ounit,'("bb00aa : ", 10x ," : lvol=",i3," ; iocons="i2" ; Semn="999es13.5)') lvol, iocons, Semn(1:mn,lvol,iocons)
   endif
#endif
   
  endif ! end of if( Igeometry.eq.3 ) ; 01 Jul 14;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(bb00aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
 end subroutine bb00aa
 
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
