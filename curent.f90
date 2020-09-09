!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (currents) ! Computes the plasma current, $I \equiv \int B_\theta \, d\theta$, and the &ldquo;linking&rdquo; current, $G \equiv \int B_\zeta \, d\zeta$.

!latex \briefly{briefly}

!latex \calledby{\link{mp00ac}, \link{dforce}}
!latex \calls{\link{coords}}

!latex \tableofcontents

!latex \subsection{enclosed currents} 

!latex \begin{enumerate}
!latex \item In the vacuum region, the enclosed currents are given by either surface integrals of the current density or line integrals of the magnetic field,
!latex       \be \int_{\cal S} {\bf j}\cdot d{\bf s} = \int_{\partial \cal S} {\bf B}\cdot d{\bf l},
!latex       \ee
!latex       and line integrals are usually easier to compute than surface integrals . . .
!latex \item The magnetic field is given by the curl of the magnetic vector potential, as described in e.g. \link{bfield}.
!latex \item The toroidal, plasma current is obtained by taking a ``poloidal'' loop, $d{\bf l}={\bf e}_\t \, d\t$, 
!latex       on the plasma boundary, where $B^s=0$, to obtain
!latex       \be I   \equiv   \int_{0}^{2\pi} \! {\bf B} \cdot {\bf e}_\t \, d\t   
!latex                  =     \int_{0}^{2\pi} \! \left( - \partial_s A_\z \; \bar g_{\t\t} + \partial_s A_\t \; \bar g_{\t\z} \right) \, d\t,
!latex       \ee
!latex       where $\bar g_{\mu\nu} \equiv g_{\mu\nu} / \sqrt g$.
!latex \item The poloidal, ``linking'' current through the torus is obtained by taking a ``toroidal'' loop, $d{\bf l}={\bf e}_\z \, d\z$, 
!latex       on the plasma boundary to obtain
!latex       \be G   \equiv   \int_{0}^{2\pi} \! {\bf B}\cdot {\bf e}_\z \, d\z   
!latex                  =     \int_{0}^{2\pi} \! \left( - \partial_s A_\z \; \bar g_{\t\z} + \partial_s A_\t \; \bar g_{\z\z} \right) \, d\z.
!latex       \ee
!latex \end{enumerate}

!latex \subsection{``Fourier integration''} 

!latex \begin{enumerate}
!l tex \item Using
!l tex       \be \begin{array}{cccccccccccccccccccccccc}\cos(\alpha+\beta) &=& \cos\alpha \cos \beta - \sin\alpha \sin \beta, \\
!l tex                                                  \cos(\alpha-\beta) &=& \cos\alpha \cos \beta + \sin\alpha \sin \beta, \end{array}
!l tex       \ee
!latex \item Using $f \equiv - \partial_s A_\z \; \bar g_{\t\t} + \partial_s A_\t \; \bar g_{\t\z}$,
!latex       the integral for the plasma current is
!latex       \be I = \sum_i^\prime f_i \cos(n_i \z) 2\pi, \label{eq:plasmacurrent}
!latex       \ee
!latex       where $\sum^\prime$ includes only the $m_i=0$ harmonics.
!latex \item Using $g \equiv - \partial_s A_\z \; \bar g_{\t\z} + \partial_s A_\t \; \bar g_{\z\z}$,
!latex       the integral for the linking current is
!latex       \be G = \sum_i^\prime g_i \cos(m_i \z) 2\pi, \label{eq:linkingcurrent}
!latex       \ee
!latex       where $\sum^\prime$ includes only the $n_i=0$ harmonics.
!latex \item     The plasma  current, \Eqn{plasmacurrent}, should be independent of $\z$,
!latex       and the linking current, \Eqn{linkingcurrent}, should be independent of $\t$.
!latex       (Perhaps this can be proved analytically; in any case it should be confirmed numerically.)
!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine curent( lvol, mn, Nt, Nz, iflag, ldItGp )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, one, two, pi2, half
  
  use numerical, only : 
  
  use fileunits, only : ounit
  
  use inputlist, only : Wmacros, Wcurent, Lrad, Lconstraint

  use cputiming, only : Tcurent

  use allglobal, only : ncpu, cpus, myid, &
                        Mvol, im, in, mne, ime, ine, &
                        YESstellsym, NOTstellsym, &
                        sg, guvij, &
                        Ntz, ijreal, ijimag, jireal, jiimag, &
                        efmn, ofmn, cfmn, sfmn, evmn, odmn, comn, simn, &
                        Ate, Aze, Ato, Azo, TT, Lcoordinatesingularity, regumm
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS 
  
  INTEGER, intent(in)  :: mn, Nt, Nz, iflag
  REAL   , intent(out) :: ldItGp(0:1,-1:2)

  INTEGER              :: innout, lvol, ideriv, ii, ll, Lcurvature, ifail
  REAL                 :: lss, mfactor
  REAL                 :: lAte(1:mn,-1:2), lAze(1:mn,-1:2), lAto(1:mn,-1:2), lAzo(1:mn,-1:2), tAze(1:mn,-1:2), zAte(1:mn,-1:2)
  REAL                 :: Bsupt(1:Nt*Nz,-1:2), Bsupz(1:Nt*Nz,-1:2), Bsups(1:Nt*Nz,-1:2), Bsups_2(1:Nt*Nz,-1:2)
  
  BEGIN(curent)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef DEBUG
  FATAL( curent, lvol.ne.Mvol, this is only defined in the vacuum region )
#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
 !ldItGp(0:1,-1:2) = zero ! initialize intent out; 12 Sep 16;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  innout = 0 ; lss = two * innout - one ! lvol = Mvol ! internal shorthand/variables; 08 Feb 16; note that lvol = Mvol is hardwired; 12 Sep 16;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  lAte(1:mn,-1:2) = zero ! radial derivatives of vector potential evaluated at interfaces; 20 Apr 13;
  lAze(1:mn,-1:2) = zero
  zAte(1:mn,-1:2) = zero
  tAze(1:mn,-1:2) = zero
! if( NOTstellsym ) then
  lAto(1:mn,-1:2) = zero ! this is used below and needs to be assigned a (trivial) value ; 26 Jan 16;
  lAzo(1:mn,-1:2) = zero ! this is used below and needs to be assigned a (trivial) value ; 26 Jan 16;
! endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if (lconstraint .eq. -2) then
    innout = 1.0
    lss = 1.0
  end if

  do ideriv = -1, 2 ! labels derivative of magnetic field wrt enclosed fluxes; 20 Apr 13;
   
   if( iflag.eq. 1 .and. ideriv.ne.0 ) cycle ! derivatives of currents                                                   are not required; 20 Jun 14;
   if( iflag.eq. 2 .and. ideriv.lt.0 ) cycle ! derivatives of currents  wrt geometry                                     is  not required; 20 Jun 14;
   if( iflag.eq.-1 .and. ideriv.gt.0 ) cycle ! derivatives of currents  wrt enclosed toroidal and enclosed poloidal flux are not required; 20 Jun 14;
   
   if( YESstellsym ) then
   
    do ii = 1, mn ! loop over Fourier harmonics; 20 Apr 13;

      if( Lcoordinatesingularity ) then ; mfactor = regumm(ii) * half ! include radial regularization factor near coordinate origin; 21 Apr 13;
        else                              ; mfactor = zero
      endif
     
     do ll = 0, Lrad(lvol) ! loop over Chebyshev polynomials; 20 Apr 13;
      
      lAte(ii,ideriv) = lAte(ii,ideriv) + Ate(lvol,ideriv,ii)%s(ll) * (TT(ll,innout,1) + mfactor) ! compute radial derivative of vector potential; 20 Apr 13;
      lAze(ii,ideriv) = lAze(ii,ideriv) - Aze(lvol,ideriv,ii)%s(ll) * (TT(ll,innout,1) + mfactor) ! note      inclusion of sign factor           ; 26 Jan 16;
      tAze(ii,ideriv) = tAze(ii,ideriv) - im(ii)*Aze(lvol,ideriv,ii)%s(ll) * TT(ll,innout,0) ! theta derivative of Az
      ZAte(ii,ideriv) = zAte(ii,ideriv) - in(ii)*Ate(lvol,ideriv,ii)%s(ll) * TT(ll,innout,0) ! zeta derivative of At
      
     enddo ! end of do ll;
     
    enddo ! end of do ii; Fourier harmonics, and their derivatives, have been constructed; 20 Apr 13;
    
   else ! if( NOTstellsym ) ; 15 Sep 16;
    
    do ii = 1, mn ! loop over Fourier harmonics; 20 Apr 13;
     
     do ll = 0, Lrad(lvol) ! loop over Chebyshev polynomials; 20 Apr 13;
      
      lAte(ii,ideriv) = lAte(ii,ideriv) + Ate(lvol,ideriv,ii)%s(ll) * TT(ll,innout,1) ! compute radial derivative of vector potential; 20 Apr 13;
      lAze(ii,ideriv) = lAze(ii,ideriv) - Aze(lvol,ideriv,ii)%s(ll) * TT(ll,innout,1) ! note      inclusion of sign factor           ; 26 Jan 16;

      lAto(ii,ideriv) = lAto(ii,ideriv) + Ato(lvol,ideriv,ii)%s(ll) * TT(ll,innout,1)
      lAzo(ii,ideriv) = lAzo(ii,ideriv) - Azo(lvol,ideriv,ii)%s(ll) * TT(ll,innout,1) ! note      inclusion of sign factor           ; 26 Jan 16;
      
     enddo ! end of do ll;
     
    enddo ! end of do ii; Fourier harmonics, and their derivatives, have been constructed; 20 Apr 13;
    
   endif ! end of if( YESstellsym ) ; 15 Sep 16;
   
   call invfft( mn, im(1:mn), in(1:mn), lAte(1:mn,ideriv), lAto(1:mn,ideriv), lAze(1:mn,ideriv), lAzo(1:mn,ideriv), &
                Nt, Nz, Bsupz(1:Ntz,ideriv), Bsupt(1:Ntz,ideriv) ) ! map to real space;

   if (Lconstraint .eq. -2) then
        call invfft( mn, im(1:mn), in(1:mn), lAto(1:mn,ideriv), tAze(1:mn,ideriv), lAzo(1:mn,ideriv), zAte(1:mn,ideriv), &
                Nt, Nz, Bsups(1:Ntz,ideriv), Bsups_2(1:Ntz,ideriv))
        Bsups(1:Ntz,ideriv) = Bsups(1:Ntz,ideriv) + Bsups_2(1:ntz,ideriv)
   else
        Bsups = zero
   end if

  enddo ! end of do ideriv; 31 Jan 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  select case( iflag ) 
  case( -1 ) ; Lcurvature = 3 ! will need derivatives of metrics w.r.t. interface deformation; 15 Sep 16;
  case(  1 ) ; Lcurvature = 1
  case(  2 ) ; Lcurvature = 1
  end select
  
  WCALL( curent, coords,( lvol, lss, Lcurvature, Ntz, mn ) ) ! get "lower" metric elements evaluated on innout interface;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  do ideriv = -1, 2 ! labels derivative of magnetic field wrt enclosed fluxes; 20 Apr 13;
   
   if( iflag.eq. 1 .and. ideriv.ne.0 ) cycle ! derivatives of currents                                                   are not required; 20 Jun 14;
   if( iflag.eq. 2 .and. ideriv.lt.0 ) cycle ! derivatives of currents  wrt geometry                                     is  not required; 20 Jun 14;
   if( iflag.eq.-1 .and. ideriv.gt.0 ) cycle ! derivatives of currents  wrt enclosed toroidal and enclosed poloidal flux are not required; 20 Jun 14;

   if (Lconstraint .eq. -2) then
   ijreal(1:Ntz) =  (Bsups(1:Ntz,ideriv) * guvij(1:Ntz,2,1,0) + Bsupt(1:Ntz,ideriv) * guvij(1:Ntz,2,2,0) + Bsupz(1:Ntz,ideriv) * guvij(1:Ntz,2,3,0) ) / sg(1:Ntz,0)
   ijimag(1:Ntz) =  (Bsups(1:Ntz,ideriv) * guvij(1:Ntz,1,3,0) + Bsupt(1:Ntz,ideriv) * guvij(1:Ntz,2,3,0) + Bsupz(1:Ntz,ideriv) * guvij(1:Ntz,3,3,0) ) / sg(1:Ntz,0)
   else
   ijreal(1:Ntz) =  (Bsupt(1:Ntz,ideriv) * guvij(1:Ntz,2,2,0) + Bsupz(1:Ntz,ideriv) * guvij(1:Ntz,2,3,0) ) / sg(1:Ntz,0)
   ijimag(1:Ntz) =  (Bsupt(1:Ntz,ideriv) * guvij(1:Ntz,2,3,0) + Bsupz(1:Ntz,ideriv) * guvij(1:Ntz,3,3,0) ) / sg(1:Ntz,0)
   end if 
   if( ideriv.eq.-1 ) then ! add derivatives of metrics with respect to interface geometry; 15 Sep 16;
   ijreal(1:Ntz) = ijreal(1:Ntz) + ( Bsupt(1:Ntz,     0) * guvij(1:Ntz,2,2,1) + Bsupz(1:Ntz,     0) * guvij(1:Ntz,2,3,1) ) / sg(1:Ntz,0)
   ijimag(1:Ntz) = ijimag(1:Ntz) + ( Bsupt(1:Ntz,     0) * guvij(1:Ntz,2,3,1) + Bsupz(1:Ntz,     0) * guvij(1:Ntz,3,3,1) ) / sg(1:Ntz,0)
   endif

   ifail = 0
   call tfft( Nt, Nz, ijreal(1:Ntz), ijimag(1:Ntz), &
              mne, ime(1:mne), ine(1:mne), efmn(1:mne), ofmn(1:mne), cfmn(1:mne), sfmn(1:mne), ifail )
   
   ldItGp(0,ideriv) = efmn(1) * pi2 ! "toroidal", plasma  current; 12 Sep 16;
   ldItGp(1,ideriv) = cfmn(1) * pi2 ! "poloidal", linking current; 12 Sep 16;
   
#ifdef DEBUG
   if( Wcurent ) then
    write(ounit,'("curent : ", 10x ," : myid=",i3," ; dItGp(0:1,",i2,")=",es23.15,",",es23.15," ;")') myid, ideriv, ldItGp(0:1,ideriv)
   endif
#endif
   
  enddo ! end of do ideriv; 31 Jan 13;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!#ifdef DEBUG
!
!  if( Wcurent ) then
!   
!   select case( ideriv )
!   case( 0 )
!    do ii = 1, mne
!     if( ime(ii).eq.0 ) write(ounit,1000) myid, ideriv, ime(ii), ine(ii), efmn(ii)
!    enddo
!    do ii = 1, mne
!     if( ine(ii).eq.0 ) write(ounit,1010) myid, ideriv, ime(ii), ine(ii), cfmn(ii)
!    enddo
!   end select
!   
!  endif ! end of if( Wcurent ) ; 05 Feb 16;
!  
!1000 format("curent : " 10x " : myid="i3" ; ideriv ="i2" ; ("i3","i3" ) : I_n ="es13.5" ;")
!1010 format("curent : " 10x " : myid="i3" ; ideriv ="i2" ; ("i3","i3" ) : G_m ="es13.5" ;")
!
!#endif
     
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(curent)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine curent

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
