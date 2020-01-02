!> \file sc00aa.f90
!> \brief The covariant components of the tangential magnetic field are computed from the singular currents at the interfaces.
!>
!> \latexonly
!> \definecolor{Orange}{rgb}{1.0,0.5,0.0}
!> \definecolor{Cerulean}{rgb}{0.0,0.5,1.0}
!> \endlatexonly

!> \brief The covariant components of the tangential magnetic field are computed from the singular currents at the interfaces.
!> \ingroup grp_diagnostics
!> 
!> Returns "covariant" magnetic field, \f${\bf B} = B_s \nabla s + B_\theta \nabla \theta + B_\zeta \nabla \zeta\f$, on \f${\cal I}\f$.
!>
!> **covariant representation**
!>
!> <ul>
!> <li> The components of the vector potential, \f${\bf A}=A_\theta \nabla + A_\zeta \nabla \zeta\f$, are
!>      \f{eqnarray}{
!>        A_\theta(s,\theta,\zeta) &=& \sum_{i,l} {\color{red}  A_{\theta,e,i,l}} \; {\overline T}_{l,i}(s) \cos\alpha_i + \sum_{i,l} {\color{Orange}  A_{\theta,o,i,l}} \; {\overline T}_{l,i}(s) \sin\alpha_i, \label{eq:At} \\
!>        A_\zeta( s,\theta,\zeta) &=& \sum_{i,l} {\color{blue} A_{\zeta, e,i,l}} \; {\overline T}_{l,i}(s) \cos\alpha_i + \sum_{i,l} {\color{Cerulean}A_{\zeta ,o,i,l}} \; {\overline T}_{l,i}(s) \sin\alpha_i, \label{eq:Az}
!>      \f}
!>      where \f${\overline T}_{l,i}(s) \equiv \bar s^{m_i/2} \, T_l(s)\f$, \f$T_l(s)\f$ is the Chebyshev polynomial, and \f$\alpha_j \equiv m_j\theta-n_j\zeta\f$.
!>      The regularity factor, \f$\bar s^{m_i/2}\f$, where \f$\bar s \equiv (1+s)/2\f$, is only included if there is a coordinate singularity in the domain
!>      (i.e. only in the innermost volume, and only in cylindrical and toroidal geometry.) </li>
!> <li> The magnetic field, \f$\sqrt g \, {\bf B} = \sqrt g B^s {\bf e}_s + \sqrt g B^\theta {\bf e}_\theta + \sqrt g B^\zeta {\bf e}_\zeta\f$, is
!>      \f{eqnarray}{
!>        \begin{array}{ccccrcrcrcrcccccccccccccccccccccccccccccccccccccccccccccccccccc}
!>        \sqrt g \, {\bf B} & = & {\bf e}_s      & \sum_{i,l} [ ( & - m_i {\color{blue}{A_{\zeta, e,i,l}}} & - & n_i {\color{red} {A_{\theta,e,i,l}}} & ) {\overline T}_{l,i}        \sin\alpha_i + ( & + m_i {\color{Cerulean}A_{\zeta ,o,i,l}} & + & n_i {\color{Orange}  A_{\theta,o,i,l}} & ) {\overline T}_{l,i}        \cos\alpha_i ] \\
!>                           & + & {\bf e}_\theta & \sum_{i,l} [ ( &                                        & - &     {\color{blue}{A_{\zeta, e,i,l}}} & ) {\overline T}_{l,i}^\prime \cos\alpha_i + ( &                                           & - &     {\color{Cerulean}A_{\zeta ,o,i,l}} & ) {\overline T}_{l,i}^\prime \sin\alpha_i ] \\
!>                           & + & {\bf e}_\zeta  & \sum_{i,l} [ ( &       {\color{red} {A_{\theta,e,i,l}}} &   &                                      & ) {\overline T}_{l,i}^\prime \cos\alpha_i + ( &       {\color{Orange}  A_{\theta,o,i,l}} &   &                                         & ) {\overline T}_{l,i}^\prime \sin\alpha_i ]
!>        \end{array}
!>      \f}
!> </li>
!> <li> On the interfaces, \f$B^s=0\f$ by construction. </li>
!> </ul>
!>
!> **output data**
!>
!> <ul>
!> <li> The Fourier harmonics of the even-and-odd, covariant components of the magnetic field, \f$B_s\f$, \f$B_\theta\f$ and \f$B_\zeta\f$, are saved in 
!>       <ul>
!>       <li>\c Btemn(1:mn,0:1,1:Mvol) , </li>
!>       <li>\c Bzemn(1:mn,0:1,1:Mvol) , </li>
!>       <li>\c Btomn(1:mn,0:1,1:Mvol) , </li>
!>       <li>\c Bzomn(1:mn,0:1,1:Mvol) ; </li>
!>       </ul>
!>       and these are written to \c ext.sp.h5 by hdfint() . </li>
!> </ul>
!> 
!> @param[in] lvol
!> @param[in] Ntz
subroutine sc00aa( lvol, Ntz )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half, one, two, pi2
  
  use numerical, only : vsmall
  
  use fileunits, only : ounit
  
  use inputlist, only : Wsc00aa, Nvol, Lrad
  
  use cputiming, only : Tsc00aa
  
  use allglobal, only : ncpu, myid, cpus, pi2nfp, &
                        Lcoordinatesingularity, Lvacuumregion, Mvol, &
                        NOTstellsym, &
                        Btemn, Bzemn, Btomn, Bzomn, &
                        mn, im, in, regumm, &
                        efmn, ofmn, cfmn, sfmn, evmn, odmn, comn, simn, &
                        Nt, Nz, &
                        Ate, Aze, Ato, Azo, &
                        TT, &
                        sg, guvij, Rij, Zij
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in)  :: lvol, Ntz
  REAL                 :: dAt(1:Ntz), dAz(1:Ntz), dBs(1:Ntz), Bt(1:Ntz), Bz(1:Ntz)
  
  INTEGER              :: innout, Lcurvature, ii, jj, kk, ll, ifail, ideriv, mi, ni
  REAL                 :: lss, mfactor
  
  BEGIN(sc00aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  FATAL( sc00aa, lvol.lt.1 .or. lvol.gt.Mvol, illegal lvol )
#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  ideriv = 0

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  do innout = 0, 1
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
  !if( lvol.eq.   1 .and. innout.eq.0 ) cycle ! REDUNDANT; 24 Nov 16;
  !if( lvol.gt.Nvol .and. innout.eq.1 ) cycle ! REDUNDANT; 24 Nov 16;

   if( ( Lcoordinatesingularity .and. innout.eq.0 ) .or. ( Lvacuumregion .and. innout.eq.1 ) ) cycle ! ignore coordinate axis and computational boundary;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   lss = two * innout - one ! recall that innout is effective local radial coordinate; 24 Apr 13;
   
   Lcurvature = 1 ! this will instruct coords to compute the metric elements; 20 Jun 14;
   
   WCALL( sc00aa, coords,( lvol, lss, Lcurvature, Ntz, mn ) ) ! get coordinates and derivatives wrt Rj, Zj, at specific radial location;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   efmn(1:mn) = zero ; ofmn(1:mn) = zero ; cfmn(1:mn) = zero ; sfmn(1:mn) = zero
   evmn(1:mn) = zero ; odmn(1:mn) = zero ; comn(1:mn) = zero ; simn(1:mn) = zero
   
   do ii = 1, mn ; mi = im(ii) ; ni = in(ii) ! loop over Fourier harmonics; 13 Sep 13;
    
    if( Lcoordinatesingularity ) then ; mfactor = regumm(ii) * half ! regularity factor; 01 Jul 14;
    else                              ; mfactor = zero
    endif
    
    do ll = 0, Lrad(lvol) ! loop over Chebyshev polynomials; Lrad is the radial resolution;
     ;                      ; efmn(ii) = efmn(ii) + Ate(lvol,ideriv,ii)%s(ll) * ( TT(ll,innout,1) + mfactor ) ! B^\t; 01 Jul 14;
     ;                      ; cfmn(ii) = cfmn(ii) + Aze(lvol,ideriv,ii)%s(ll) * ( TT(ll,innout,1) + mfactor ) ! B^\z; 01 Jul 14;
     ;                      ; odmn(ii) = odmn(ii) + Ate(lvol,ideriv,ii)%s(ll) * ( TT(ll,innout,1) + mfactor ) * ni & 
                                                  + Aze(lvol,ideriv,ii)%s(ll) * ( TT(ll,innout,1) + mfactor ) * mi
     if( NOTstellsym ) then ; ofmn(ii) = ofmn(ii) + Ato(lvol,ideriv,ii)%s(ll) * ( TT(ll,innout,1) + mfactor )
      ;                     ; sfmn(ii) = sfmn(ii) + Azo(lvol,ideriv,ii)%s(ll) * ( TT(ll,innout,1) + mfactor )
      ;                     ; evmn(ii) = evmn(ii) - Ato(lvol,ideriv,ii)%s(ll) * ( TT(ll,innout,1) + mfactor ) * ni & 
                                                  - Azo(lvol,ideriv,ii)%s(ll) * ( TT(ll,innout,1) + mfactor ) * mi
    !else                   ; ofmn(ii) = zero ! defaulted to zero above; 11 Mar 16;
    ! ;                     ; sfmn(ii) = zero
    ! ;                     ; evmn(ii) = zero
     endif
    enddo ! end of do ll; 20 Feb 13;
    
   enddo ! end of do ii; 20 Feb 13;
   
   call invfft( mn, im, in, efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, dAt(1:Ntz), dAz(1:Ntz) )
   
   Bt(1:Ntz) = ( - dAz(1:Ntz) * guvij(1:Ntz,2,2,0) + dAt(1:Ntz) * guvij(1:Ntz,2,3,0) ) / sg(1:Ntz,0)
   Bz(1:Ntz) = ( - dAz(1:Ntz) * guvij(1:Ntz,2,3,0) + dAt(1:Ntz) * guvij(1:Ntz,3,3,0) ) / sg(1:Ntz,0)
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
   ifail = 0
   
   call tfft( Nt, Nz, Bt(1:Ntz), Bz(1:Ntz), &
              mn, im(1:mn), in(1:mn), Btemn(1:mn,innout,lvol), Btomn(1:mn,innout,lvol), Bzemn(1:mn,innout,lvol), Bzomn(1:mn,innout,lvol), ifail )
   
   if( Wsc00aa                    ) then
    write(ounit,'("sc00aa : "10x " : myid=",i3," ; lvol=",i3," ; innout="i2" ; dB^somn =",99es13.05)') myid, lvol, innout, odmn(1:mn)
    write(ounit,'("sc00aa : "10x " : myid=",i3," ; lvol=",i3," ; innout="i2" ;  B_temn =",99es13.05)') myid, lvol, innout, Btemn(1:mn,innout,lvol)
    write(ounit,'("sc00aa : "10x " : myid=",i3," ; lvol=",i3," ; innout="i2" ;  B_zemn =",99es13.05)') myid, lvol, innout, Bzemn(1:mn,innout,lvol)
   endif
   if( Wsc00aa .and. NOTstellsym ) then
    write(ounit,'("sc00aa : "10x " : myid=",i3," ; lvol=",i3," ; innout="i2" ;  B_tomn =",99es13.05)') myid, lvol, innout, Btomn(1:mn,innout,lvol)
    write(ounit,'("sc00aa : "10x " : myid=",i3," ; lvol=",i3," ; innout="i2" ;  B_zomn =",99es13.05)') myid, lvol, innout, Bzomn(1:mn,innout,lvol)
   endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
  enddo ! end of do innout = 0, 1 ; 20 Jun 14;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(sc00aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
 end subroutine sc00aa
 
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
