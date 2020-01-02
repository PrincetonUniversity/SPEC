!> \file jo00aa.f90
!> \brief Measures error in Beltrami equation, \f$\nabla \times {\bf B} - \mu {\bf B}\f$.

!> \brief Measures error in Beltrami equation, \f$\nabla \times {\bf B} - \mu {\bf B}\f$.
!> \ingroup grp_diagnostics
!>
!> This routine is called by xspech() as a post diagnostic and only if \c Lcheck==1.
!>
!> **construction of current,** \f${\bf j} \equiv \nabla \times \nabla \times {\bf A}\f$
!> <ul>
!> <li> The components of the vector potential, \f${\bf A}=A_\theta \nabla + A_\zeta \nabla \zeta\f$, are
!>      \f{eqnarray}{
!>        A_\theta(s,\theta,\zeta) &=& \sum_{i,l} {\color{red} {A_{\theta,e,i,l}}} \; {\overline T}_{l,i}(s) \cos\alpha_i + \sum_{i,l} {\color{magenta}{A_{\theta,o,i,l}}} \; {\overline T}_{l,i}(s) \sin\alpha_i, \label{eq:At} \\
!>        A_\zeta (s,\theta,\zeta) &=& \sum_{i,l} {\color{blue}{A_{\zeta, e,i,l}}} \; {\overline T}_{l,i}(s) \cos\alpha_i + \sum_{i,l} {\color{cyan}   {A_{\zeta ,o,i,l}}} \; {\overline T}_{l,i}(s) \sin\alpha_i, \label{eq:Az}
!>      \f}
!>      where \f${\overline T}_{l,i}(s) \equiv \bar s^{m_i/2} \, T_l(s)\f$, \f$T_l(s)\f$ is the Chebyshev polynomial, and \f$\alpha_j \equiv m_j\theta-n_j\zeta\f$.
!>      The regularity factor, \f$\bar s^{m_i/2}\f$, where \f$\bar s \equiv (1+s)/2\f$, is only included if there is a coordinate singularity in the domain
!>      (i.e. only in the innermost volume, and only in cylindrical and toroidal geometry.) </li>
!> <li> The magnetic field, \f$\sqrt g \, {\bf B} = \sqrt g B^s {\bf e}_s + \sqrt g B^\theta {\bf e}_\theta + \sqrt g B^\zeta {\bf e}_\zeta\f$, is
!>      \f{eqnarray}{
!>        \begin{array}{ccccrcrcrcrcccccccccccccccccccccccccccccccccccccccccccccccccccc}
!>        \sqrt g \, {\bf B} & = & {\bf e}_s      & \sum_{i,l} [ ( & - m_i {\color{blue}{A_{\zeta, e,i,l}}} & - & n_i {\color{red} {A_{\theta,e,i,l}}} & ) {\overline T}_{l,i}        \sin\alpha_i + ( & + m_i {\color{cyan}   {A_{\zeta ,o,i,l}}} & + & n_i {\color{magenta}{A_{\theta,o,i,l}}} & ) {\overline T}_{l,i}        \cos\alpha_i ] \\
!>                           & + & {\bf e}_\theta & \sum_{i,l} [ ( &                                        & - &     {\color{blue}{A_{\zeta, e,i,l}}} & ) {\overline T}_{l,i}^\prime \cos\alpha_i + ( &                                           & - &     {\color{cyan}   {A_{\zeta ,o,i,l}}} & ) {\overline T}_{l,i}^\prime \sin\alpha_i ] \\
!>                           & + & {\bf e}_\zeta  & \sum_{i,l} [ ( &       {\color{red} {A_{\theta,e,i,l}}} &   &                                      & ) {\overline T}_{l,i}^\prime \cos\alpha_i + ( &       {\color{magenta}{A_{\theta,o,i,l}}} &   &                                         & ) {\overline T}_{l,i}^\prime \sin\alpha_i ]
!>        \end{array}
!>      \f} </li>
!> <li> The current is
!>      \f{eqnarray}{ \sqrt g \, {\bf j} = ( \partial_\theta B_\zeta  - \partial_\zeta  B_\theta) \; {\bf e}_s
!>                                       + ( \partial_\zeta  B_s      - \partial_s      B_\zeta ) \; {\bf e}_\theta
!>                                       + ( \partial_s      B_\theta - \partial_\theta B_s     ) \; {\bf e}_\zeta ,
!>      \f}
!>      where (for computational convenience) the covariant components of \f${\bf B}\f$ are computed as
!>      \f{eqnarray}{
!>        B_s      & = & (\sqrt g B^s) \, g_{ss     } / \sqrt g + (\sqrt g B^\theta) \, g_{s     \theta} / \sqrt g + (\sqrt g B^\zeta) \, g_{s     \zeta} / \sqrt g, \\
!>        B_\theta & = & (\sqrt g B^s) \, g_{s\theta} / \sqrt g + (\sqrt g B^\theta) \, g_{\theta\theta} / \sqrt g + (\sqrt g B^\zeta) \, g_{\theta\zeta} / \sqrt g, \\
!>        B_\zeta  & = & (\sqrt g B^s) \, g_{s\zeta } / \sqrt g + (\sqrt g B^\theta) \, g_{\theta\zeta } / \sqrt g + (\sqrt g B^\zeta) \, g_{\zeta \zeta} / \sqrt g. 
!>      \f} </li>
!> </ul> 
!>
!> **quantification of the error**
!> <ul>
!> <li> The measures of the error are
!>      \f{eqnarray}{
!>          ||\left( {\bf j}-\mu {\bf B}\right)\cdot\nabla s      || &  \equiv  & 
!>        \int \!\! ds \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \left| \sqrt g \, {\bf j}\cdot\nabla s     -\mu \; \sqrt g \, {\bf B} \cdot\nabla s      \right|, \label{eq:Es} \\
!>          ||\left( {\bf j}-\mu {\bf B}\right)\cdot\nabla \theta || &  \equiv  & 
!>        \int \!\! ds \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \left| \sqrt g \, {\bf j}\cdot\nabla \theta-\mu \; \sqrt g \, {\bf B} \cdot\nabla \theta \right|, \label{eq:Et} \\
!>          ||\left( {\bf j}-\mu {\bf B}\right)\cdot\nabla \zeta  || &  \equiv  & 
!>        \int \!\! ds \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \left| \sqrt g \, {\bf j}\cdot\nabla \zeta -\mu \; \sqrt g \, {\bf B} \cdot\nabla \zeta  \right|. \label{eq:Ez} 
!>      \f} </li>
!> </ul>
!>
!> **comments**
!> <ul>
!> <li>  Is there a better definition and quantification of the error? For example, should we employ an error measure that is dimensionless? </li>
!> <li>  If the coordinate singularity is in the domain, then \f$|\nabla\theta|\rightarrow\infty\f$ at the coordinate origin. 
!>       What then happens to \f$||\left( {\bf j}-\mu {\bf B}\right)\cdot\nabla \theta ||\f$ as defined in Eqn.\f$(\ref{eq:Et})\f$? </li>
!> <li>  What is the predicted scaling of the error in the Chebyshev-Fourier representation scale with numerical resolution?
!>       Note that the predicted error scaling for \f$E^s\f$, \f$E^\theta\f$ and \f$E^\zeta\f$ may not be standard, 
!>       as various radial derivatives are taken to compute the components of \f${\bf j}\f$.
!>       (See for example the discussion in Sec.IV.C in Hudson et al. (2012) \cite y2012_hudson ,
!>       where the expected scaling of the error for a finite-element implementation is confirmed numerically.) </li>
!> <li>  Instead of using Gaussian integration to compute the integral over \f$s\f$, an adaptive quadrature algorithm may be preferable. </li>
!> </ul>
!>
!> @param[in] lvol  in which volume should the Beltrami error be computed
!> @param[in] Ntz   number of grid points in \f$\theta\f$ and \f$\zeta\f$
!> @param[in] lquad degree of Gaussian quadrature
!> @param[in] mn    number of Fourier harmonics
subroutine jo00aa( lvol, Ntz, lquad, mn )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  use constants, only : zero, half, one, two

  use fileunits, only : ounit

  use inputlist, only : Wmacros, Wjo00aa, ext, Nvol, Lrad, mu

  use cputiming, only : Tjo00aa

  use allglobal, only : myid, cpus, ivol, &
                        im, in, regumm, &
                        Mvol, &
                        cheby, &
                        Ate, Aze, Ato, Azo, &
                        sg, guvij, Rij, Zij, &
                        Nt, Nz, efmn, ofmn, cfmn, sfmn, &
                        NOTstellsym, &
                        Lcoordinatesingularity, &
                        beltramierror  

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS
  
!                        these are really global, but are included in argument list to remove allocations
  INTEGER, intent(in) :: lvol, Ntz, lquad, mn 
  
  INTEGER             :: jquad, Lcurvature, ll, ii, jj, kk, uu, ideriv, twolquad
  
  REAL                :: lss, sbar, sbarhim(0:2), gBu(1:Ntz,1:3,0:3), gJu(1:Ntz,1:3), jerror(1:3)
  
  REAL                :: Atemn(1:mn,0:2), Azemn(1:mn,0:2), Atomn(1:mn,0:2), Azomn(1:mn,0:2)
  
  INTEGER             :: itype, icdgqf
  REAL                :: aa, bb, cc, dd, weight(1:lquad), abscis(1:lquad), workfield(1:2*lquad)
  
  BEGIN(jo00aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  FATAL( jo00aa, lquad.lt.1, invalid lquad supplied to jo00aa )
  FATAL( jo00aa, lvol.lt.1 .or. lvol.gt.Mvol, invalid interface label )
#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **details of the numerics**
!> <ul>
!> <li> The integration over \f$s\f$ is performed using Gaussian integration, e.g., \f$\displaystyle \int \!\! f(s) ds \approx \sum_k \omega_k f(s_k)\f$;
!>      with the abscissae, \f$s_k\f$, and the weights, \f$\omega_k\f$, for \f$k=1,\f$ \c Iquad\f$_v\f$, determined by \c CDGQF.
!>      The resolution, \c N \f$ \equiv \f$ \c Iquad\f$_v\f$, is determined by \c Nquad  (see global.f90 and preset.f90).
!>      A fatal error is enforced by jo00aa() if \c CDGQF returns an \c ifail \f$\ne 0\f$. </li>
  
  itype = 1 ; aa = -one ; bb = +one ; cc = zero ; dd = zero ; twolquad = 2 * lquad
  
  call CDGQF( lquad, abscis(1:lquad), weight(1:lquad), itype, aa, bb, twolquad, workfield(1:twolquad), icdgqf ) ! prepare Gaussian quadrature;
  
! write(ounit,'("jo00aa :  WARNING ! : THE ERROR FLAGS RETURNED BY CDGQF SEEM DIFFERENT TO NAG:D01BCF (may be trivial, but please revise); 2018/01/10;")')

  cput= GETTIME
  select case( icdgqf ) !                                                         123456789012345
  case( 0 )    ;  if( Wjo00aa ) write(ounit,1000) cput-cpus, myid, lvol, icdgqf, "success        ", abscis(1:lquad)
  case( 1 )    ;                write(ounit,1000) cput-cpus, myid, lvol, icdgqf, "failed         ", abscis(1:lquad)
  case( 2 )    ;                write(ounit,1000) cput-cpus, myid, lvol, icdgqf, "input error    ", abscis(1:lquad)
  case( 3 )    ;                write(ounit,1000) cput-cpus, myid, lvol, icdgqf, "input error    ", abscis(1:lquad)
  case( 4 )    ;                write(ounit,1000) cput-cpus, myid, lvol, icdgqf, "weight overflow", abscis(1:lquad)
  case( 5 )    ;                write(ounit,1000) cput-cpus, myid, lvol, icdgqf, "weight zero    ", abscis(1:lquad)
  case( 6 )    ;                write(ounit,1000) cput-cpus, myid, lvol, icdgqf, "failed         ", abscis(1:lquad)
  case default ;                write(ounit,1000) cput-cpus, myid, lvol, icdgqf, "weird          ", abscis(1:lquad)
  end select
  
  if( Wjo00aa )                 write(ounit,1001)                                                   weight(1:lquad)
  
1000 format("jo00aa : ",f10.2," : myid=",i3," ; lvol=",i3," ; icdgqf=",i3," ; "a15" ;":" abscissae ="99f10.6)
1001 format("jo00aa : ", 10x ," :       "3x"          "3x"            "3x"    "15x" ;":" weights   ="99f10.6)
  
  FATAL( jo00aa, icdgqf.ne.0, failed to construct Gaussian integration abscisae and weights )
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  jerror(1:3) = zero ; ideriv = 0 ! three components of the error in \curl B - mu B; initialize summation;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!> <li>  Inside the Gaussian quadrature loop, i.e. for each \f$s_k\f$, 
!>       <ul>
  
  do jquad = 1, lquad ! loop over radial sub-sub-grid (numerical quadrature);
   
   lss = abscis(jquad) ; sbar = ( lss + one ) * half
   
!>       <li>  The metric elements, \f$g_{\mu,\nu} \equiv \f$ \c gij(1:6,0,1:Ntz), and the Jacobian, \f$\sqrt g \equiv \f$ \c sg(0,1:Ntz),
!>             are calculated on a regular angular grid, \f$(\theta_i,\zeta_j)\f$, in coords().
!>             The derivatives \f$\partial_i g_{\mu,\nu} \equiv\f$ \c gij(1:6,i,1:Ntz) and \f$\partial_i \sqrt g \equiv \f$ \c sg(i,1:Ntz),
!>             with respect to \f$ i \in \{ s,\theta,\zeta \}\f$ are also returned. </li>
   
   Lcurvature = 2
   
   WCALL( jo00aa, coords, ( lvol, lss, Lcurvature, Ntz, mn ) ) ! returns coordinates, metrics, . . .
   
   ;                     ; cheby( 0,0:2) = (/ one, zero, zero /) ! T_0: Chebyshev initialization; function, 1st-derivative, 2nd-derivative;
   ;                     ; cheby( 1,0:2) = (/ lss,  one, zero /) ! T_1: Chebyshev initialization; function, 1st-derivative, 2nd-derivative;
   do ll = 2, Lrad(lvol) ; cheby(ll,0:2) = (/ two * lss * cheby(ll-1,0)                                                         - cheby(ll-2,0) , &
                                              two       * cheby(ll-1,0) + two * lss * cheby(ll-1,1)                             - cheby(ll-2,1) , &
                                              two       * cheby(ll-1,1) + two       * cheby(ll-1,1) + two * lss * cheby(ll-1,2) - cheby(ll-2,2) /)
   enddo ! end of do ll; 20 Jun 14;
    
   Atemn(1:mn,0:2) = zero ! initialize summation over Chebyshev polynomials;
   Azemn(1:mn,0:2) = zero
   if( NOTstellsym ) then
   Atomn(1:mn,0:2) = zero
   Azomn(1:mn,0:2) = zero
   else
   Atomn(1:mn,0:2) = zero ! these are used below;
   Azomn(1:mn,0:2) = zero
   endif
   
   do ll = 0, Lrad(lvol) ! radial (Chebyshev) resolution of magnetic vector potential;
   
    do ii = 1, mn  ! Fourier resolution of magnetic vector potential;
     
!>          <li>  The Fourier components of the vector potential given in Eqn.\f$(\ref{eq:At})\f$ and Eqn.\f$(\ref{eq:Az})\f$, and their first and second radial derivatives, are summed. </li>

     if( Lcoordinatesingularity ) then ! compute regularization factor; 10 Jan 2018;
      sbarhim(0) = sbar**regumm(ii) ; sbarhim(1) = half * regumm(ii) * sbarhim(0) / sbar ; sbarhim(2) = half * ( regumm(ii)-one ) * sbarhim(1) / sbar
     else                              
      sbarhim(0) = one              ; sbarhim(1) = zero                                  ; sbarhim(2) = zero
     endif
     
     ;Atemn(ii,0) = Atemn(ii,0) + Ate(lvol,ideriv,ii)%s(ll) * ( cheby(ll,0) * sbarhim(0)                                                             )
     ;Atemn(ii,1) = Atemn(ii,1) + Ate(lvol,ideriv,ii)%s(ll) * ( cheby(ll,1) * sbarhim(0)                                  + cheby(ll,0) * sbarhim(1) )
     ;Atemn(ii,2) = Atemn(ii,2) + Ate(lvol,ideriv,ii)%s(ll) * ( cheby(ll,2) * sbarhim(0) + two * cheby(ll,1) * sbarhim(1) + cheby(ll,0) * sbarhim(2) )

     ;Azemn(ii,0) = Azemn(ii,0) + Aze(lvol,ideriv,ii)%s(ll) * ( cheby(ll,0) * sbarhim(0)                                                             )
     ;Azemn(ii,1) = Azemn(ii,1) + Aze(lvol,ideriv,ii)%s(ll) * ( cheby(ll,1) * sbarhim(0)                                  + cheby(ll,0) * sbarhim(1) )
     ;Azemn(ii,2) = Azemn(ii,2) + Aze(lvol,ideriv,ii)%s(ll) * ( cheby(ll,2) * sbarhim(0) + two * cheby(ll,1) * sbarhim(1) + cheby(ll,0) * sbarhim(2) )

     if( NOTstellsym ) then

      Atomn(ii,0) = Atomn(ii,0) + Ato(lvol,ideriv,ii)%s(ll) * ( cheby(ll,0) * sbarhim(0)                                                             )
      Atomn(ii,1) = Atomn(ii,1) + Ato(lvol,ideriv,ii)%s(ll) * ( cheby(ll,1) * sbarhim(0)                                  + cheby(ll,0) * sbarhim(1) )
      Atomn(ii,2) = Atomn(ii,2) + Ato(lvol,ideriv,ii)%s(ll) * ( cheby(ll,2) * sbarhim(0) + two * cheby(ll,1) * sbarhim(1) + cheby(ll,0) * sbarhim(2) )

      Azomn(ii,0) = Azomn(ii,0) + Azo(lvol,ideriv,ii)%s(ll) * ( cheby(ll,0) * sbarhim(0)                                                             )
      Azomn(ii,1) = Azomn(ii,1) + Azo(lvol,ideriv,ii)%s(ll) * ( cheby(ll,1) * sbarhim(0)                                  + cheby(ll,0) * sbarhim(1) )
      Azomn(ii,2) = Azomn(ii,2) + Azo(lvol,ideriv,ii)%s(ll) * ( cheby(ll,2) * sbarhim(0) + two * cheby(ll,1) * sbarhim(1) + cheby(ll,0) * sbarhim(2) )

     endif ! end of if( NOTstellsym) ; 20 Jun 14;
      
    enddo ! end of do ii;
    
   enddo ! end of do ll;
    
!>           <li>  The quantities \f$\sqrt g B^s\f$, \f$\sqrt g B^\theta\f$ and \f$\sqrt g B^\zeta\f$, and their first and second derivatives with respect to \f$(s,\theta,\zeta)\f$, 
!>                 are computed on the regular angular grid (using FFTs). </li>

   ofmn(1:mn) = -          im(1:mn)*Azemn(1:mn,0) -          in(1:mn)*Atemn(1:mn,0)
   efmn(1:mn) = +          im(1:mn)*Azomn(1:mn,0) +          in(1:mn)*Atomn(1:mn,0)
   sfmn(1:mn) = -          im(1:mn)*Azemn(1:mn,1) -          in(1:mn)*Atemn(1:mn,1)
   cfmn(1:mn) = +          im(1:mn)*Azomn(1:mn,1) +          in(1:mn)*Atomn(1:mn,1)
   
   call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, gBu(1:Ntz,1,0), gBu(1:Ntz,1,1) ) !  (gB^s)   , d(gB^s)/ds;
   
   efmn(1:mn) = - im(1:mn)*im(1:mn)*Azemn(1:mn,0) - im(1:mn)*in(1:mn)*Atemn(1:mn,0)
   ofmn(1:mn) = - im(1:mn)*im(1:mn)*Azomn(1:mn,0) - im(1:mn)*in(1:mn)*Atomn(1:mn,0)
   cfmn(1:mn) = + in(1:mn)*im(1:mn)*Azemn(1:mn,0) + in(1:mn)*in(1:mn)*Atemn(1:mn,0)
   sfmn(1:mn) = + in(1:mn)*im(1:mn)*Azomn(1:mn,0) + in(1:mn)*in(1:mn)*Atomn(1:mn,0)
   
   call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, gBu(1:Ntz,1,2), gBu(1:Ntz,1,3) ) ! d(gB^s)/dt, d(gB^s)/dz;
   
   efmn(1:mn) =                                   -                   Azemn(1:mn,1)
   ofmn(1:mn) =                                   -                   Azomn(1:mn,1)
   cfmn(1:mn) =                                   -                   Azemn(1:mn,2)
   sfmn(1:mn) =                                   -                   Azomn(1:mn,2)
   
   call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, gBu(1:Ntz,2,0), gBu(1:Ntz,2,1) ) !  (gB^t)   , d(gB^t)/ds;

   ofmn(1:mn) =                                   +          im(1:mn)*Azemn(1:mn,1)
   efmn(1:mn) =                                   -          im(1:mn)*Azomn(1:mn,1)
   sfmn(1:mn) =                                   -          in(1:mn)*Azemn(1:mn,1)
   cfmn(1:mn) =                                   +          in(1:mn)*Azomn(1:mn,1)
   
   call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, gBu(1:Ntz,2,2), gBu(1:Ntz,2,3) ) ! d(gB^t)/dt, d(gB^t)/dz;
   
   efmn(1:mn) = +                   Atemn(1:mn,1)
   ofmn(1:mn) = +                   Atomn(1:mn,1)
   cfmn(1:mn) = +                   Atemn(1:mn,2)
   sfmn(1:mn) = +                   Atomn(1:mn,2)
   
   call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, gBu(1:Ntz,3,0), gBu(1:Ntz,3,1) ) !  (gB^z)   , d(gB^z)/ds;

   ofmn(1:mn) = -          im(1:mn)*Atemn(1:mn,1)
   efmn(1:mn) = +          im(1:mn)*Atomn(1:mn,1)
   sfmn(1:mn) = +          in(1:mn)*Atemn(1:mn,1)
   cfmn(1:mn) = -          in(1:mn)*Atomn(1:mn,1)
   
   call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, gBu(1:Ntz,3,2), gBu(1:Ntz,3,3) ) ! d(gB^z)/dt, d(gB^z)/dz;
     
!>            <li>  The following quantities are then computed on the regular angular grid
!>                  \f{eqnarray}{ \sqrt g j^s      & = & \sum_u \left[
!>                                \partial_\theta (\sqrt g B^u) \; g_{u,\zeta } + (\sqrt g B^u) \; \partial_\theta g_{u,\zeta } - (\sqrt g B^u) g_{u,\zeta } \; \partial_\theta \sqrt g / \sqrt g
!>                                                         \right] / \sqrt g \nonumber \\
!>                                            & - & \sum_u \left[
!>                                \partial_\zeta  (\sqrt g B^u) \; g_{u,\theta} + (\sqrt g B^u) \; \partial_\zeta  g_{u,\theta} - (\sqrt g B^u) g_{u,\theta} \; \partial_\zeta  \sqrt g / \sqrt g
!>                                                         \right] / \sqrt g, \\
!>                                \sqrt g j^\theta & = & \sum_u \left[
!>                                \partial_\zeta  (\sqrt g B^u) \; g_{u,s     } + (\sqrt g B^u) \; \partial_\zeta  g_{u,s     } - (\sqrt g B^u) g_{u,s     } \; \partial_\zeta  \sqrt g / \sqrt g
!>                                                         \right] / \sqrt g \nonumber \\
!>                                            & - & \sum_u \left[
!>                                \partial_s      (\sqrt g B^u) \; g_{u,\zeta } + (\sqrt g B^u) \; \partial_s      g_{u,\zeta } - (\sqrt g B^u) g_{u,\zeta } \; \partial_s      \sqrt g / \sqrt g
!>                                                         \right] / \sqrt g, \\
!>                                \sqrt g j^\zeta  & = & \sum_u \left[
!>                                \partial_s      (\sqrt g B^u) \; g_{u,\theta} + (\sqrt g B^u) \; \partial_s      g_{u,\theta} - (\sqrt g B^u) g_{u,\theta} \; \partial_s      \sqrt g / \sqrt g
!>                                                         \right] / \sqrt g \nonumber \\
!>                                            & - & \sum_u \left[
!>                                \partial_\theta (\sqrt g B^u) \; g_{u,s     } + (\sqrt g B^u) \; \partial_\theta g_{u,s     } - (\sqrt g B^u) g_{u,s     } \; \partial_\theta \sqrt g / \sqrt g
!>                                                         \right] / \sqrt g.
!>                  \f} </li>

   do ii = 1, 3
    
    select case( ii )
    case( 1 ) ; jj = 2 ; kk = 3
    case( 2 ) ; jj = 3 ; kk = 1
    case( 3 ) ; jj = 1 ; kk = 2
    end select
    
    gJu(1:Ntz,ii) = zero
    
    do uu = 1, 3 ! summation over uu;
     
     gJu(1:Ntz,ii) = gJu(1:Ntz,ii) &
                   + ( gBu(1:Ntz,uu,jj) * guvij(1:Ntz,uu,kk, 0) &
                     + gBu(1:Ntz,uu, 0) * guvij(1:Ntz,uu,kk,jj) &
                     - gBu(1:Ntz,uu, 0) * guvij(1:Ntz,uu,kk, 0) * sg(1:Ntz,jj) / sg(1:Ntz,0) &
                     ) &
                   - ( gBu(1:Ntz,uu,kk) * guvij(1:Ntz,uu,jj, 0) &
                     + gBu(1:Ntz,uu, 0) * guvij(1:Ntz,uu,jj,kk) &
                     - gBu(1:Ntz,uu, 0) * guvij(1:Ntz,uu,jj, 0) * sg(1:Ntz,kk) / sg(1:Ntz,0) &
                     )
    enddo ! end of do uu;
     
    gJu(1:Ntz,ii) = gJu(1:Ntz,ii) / sg(1:Ntz,0)
    
   enddo ! end of do ii;

!>             </ul>
!> </li>
!> <li> The final calculation of the error, which is written to screen, is a sum over the angular grid:
!>      \f{eqnarray}{ E^s      & \equiv & \frac{1}{N} \sum_k \omega_k \sum_{i,j} | \sqrt g j^s      - \mu \sqrt g B^s      |, \\
!>                    E^\theta & \equiv & \frac{1}{N} \sum_k \omega_k \sum_{i,j} | \sqrt g j^\theta - \mu \sqrt g B^\theta |, \\
!>                    E^\zeta  & \equiv & \frac{1}{N} \sum_k \omega_k \sum_{i,j} | \sqrt g j^\zeta  - \mu \sqrt g B^\zeta  |,
!>      \f}
!>      where \f$N\equiv \sum_{i,j}1\f$. </li>

   do ii = 1, 3 ; jerror(ii) = jerror(ii) + weight(jquad) * sum( abs(  gJu(1:Ntz,ii) - mu(lvol) * gBu(1:Ntz,ii,0)  ) )
   enddo
   
  enddo ! end of do jquad;
   

  jerror(1:3) = jerror(1:3) / Ntz
   
  cput = GETTIME ; write(ounit,1002) cput-cpus, myid, lvol, Lrad(lvol), jerror(1:3), cput-cpui ! write error to screen;

!> <li>  The error is stored into an array called \c beltramierror which is then written to the HDF5 file in hdfint(). </li>

  beltramierror(lvol,1:3) = jerror(1:3)   
   
1002 format("jo00aa : ",f10.2," : myid=",i3," ; lvol =",i3," ; lrad =",i3," ; E^\s="es23.15" , E^\t="es23.15" , E^\z="es23.15" ; time="f8.2"s ;")
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(jo00aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> </ul>

end subroutine jo00aa
