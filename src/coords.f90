!> \defgroup grp_geometry Geometry
!>
!> \file
!> \brief Calculates coordinates, \f${\bf x}(s,\theta,\zeta) \equiv R \, {\bf e}_R + Z \, {\bf e}_Z\f$, and metrics, using FFTs.

!> \brief Calculates coordinates, \f${\bf x}(s,\theta,\zeta) \equiv R \, {\bf e}_R + Z \, {\bf e}_Z\f$, and metrics, using FFTs.
!> \ingroup grp_geometry
!>
!> **Coordinates**
!> <ul>
!> <li>  We work in coordinates, \f$(s,\theta,\zeta)\f$, which are be defined _inversely_ via a transformation _to_ Cartesian coordinates, \f$(x,y,z)\f$. </li>
!> <li>  The toroidal angle, \f$\zeta\f$, is identical to the cylindrical angle, \f$\zeta\equiv\phi\f$. </li>
!> <li>  The radial coordinate, \f$s\f$, is _not_ a global variable: it only needs to be defined in each volume, and in each volume \f$s \in [-1,1]\f$. </li>
!> <li>  The choice of poloidal angle, \f$\theta\f$, does not affect the following. </li>
!> </ul>
!>
!> **Geometry**
!> <ul>
!> <li>  The geometry of the "ideal"-interfaces, \f${\bf x}_v(\theta,\zeta)\f$, is given by \f$R(\theta,\zeta)\f$ and \f$Z(\theta,\zeta)\f$ as follows:
!>       <ul>
!>       <li>  \c Igeometry=1 : Cartesian
!>       \f{eqnarray}{ {\bf x} & \equiv & r_{pol}\theta \; {\bf \hat i} + r_{tor}\zeta \; {\bf \hat j}+ R \; {\bf \hat k}
!>       \f}
!>       where \f$r_{pol}\f$ and \f$r_{tor}\f$ are inputs and \f$r_{pol}=r_{tor}=1\f$ by default. </li>
!>       <li>  \c Igeometry=2 : Cylindrical
!>       \f{eqnarray}{ {\bf x} & = & R \; \cos\theta \; {\bf \hat i} + R \; \sin\theta \; {\bf \hat j} + \zeta \; {\bf \hat k}
!>       \f} </li>
!>       <li>  \c Igeometry=3 : Toroidal
!>       \f{eqnarray}{ {\bf x} & \equiv & R \; {\bf \hat r} + Z \; {\bf \hat k}
!>       \f}
!>       where \f${\bf \hat    r} \equiv   \cos \phi \; {\bf \hat i} + \sin \phi \; {\bf \hat j}\f$ and
!>             \f${\bf \hat \phi} \equiv - \sin \phi \; {\bf \hat i} + \cos \phi \; {\bf \hat j}\f$. </li>
!>       </ul> </li>
!>
!> <li>  The geometry of the ideal interfaces is given as Fourier summation: e.g., for stellarator-symmetry
!>       \f{eqnarray}{ R_v(\theta,\zeta) & \equiv & \sum_j R_{j,v} \cos\alpha_j, \\
!>           Z_v(\theta,\zeta) & \equiv & \sum_j Z_{j,v} \sin\alpha_j,
!>       \f}
!>       where \f$\alpha_j \equiv m_j \theta - n_j \zeta\f$. </li>
!>
!> </ul>
!>       **interpolation between interfaces**
!> <ul>
!> <li>  The "coordinate" functions, \f$R(s,\theta,\zeta)\f$ and \f$Z(s,\theta,\zeta)\f$, are constructed by radially interpolating the Fourier representations of the
!>       ideal-interfaces. </li>
!> <li>  The \f$v\f$-th volume is bounded by \f${\bf x}_{v-1}\f$ and \f${\bf x}_{v}\f$. </li>
!> <li>  In each _annular_ volume, the coordinates are constructed by linear interpolation:
!>       \f{eqnarray}{ \begin{array}{cccccccccccccccccccccccccccc}
!>           R(s,\theta,\zeta) & \equiv & \displaystyle \sum_j & \displaystyle \left[ \; \frac{(1-s)}{2} \; R_{j,v-1} + \frac{(1+s)}{2} \; R_{j,v}\; \right] \; \cos\alpha_j,\\
!>           Z(s,\theta,\zeta) & \equiv & \displaystyle \sum_j & \displaystyle \left[ \; \frac{(1-s)}{2} \; Z_{j,v-1} + \frac{(1+s)}{2} \; Z_{j,v}\; \right] \; \sin\alpha_j,
!>           \end{array}
!>       \f} </li>
!> </ul>
!>       **coordinate singularity: regularized extrapolation**
!> <ul>
!>
!> <li>  For cylindrical or toroidal geometry, in the innermost, "simple-torus" volume, the coordinates are constructed by an interpolation that
!>       "encourages" the interpolated coordinate surfaces to not intersect. </li>
!> <li>  Introduce \f$\bar s \equiv (s+1)/2\f$, so that in each volume \f$\bar s \in [0,1]\f$, then
!>       \f{eqnarray}{ R_{j}(s) & = & R_{j,0} + (R_{j,1} - R_{j,0} ) f_j, \\
!>           Z_{j}(s) & = & Z_{j,0} + (Z_{j,1} - Z_{j,0} ) f_j,
!>       \f}
!>       where, in toroidal geometry,
!>       \f{eqnarray}{
!>       f_j \equiv \left\{
!>       \begin{array}{llcccccccccccccc} \bar s       & , & \mbox{ for } m_j=0, \\
!>                                       \bar s^{m_j} & , & \mbox{ otherwise.}
!>       \end{array}\right\}.
!>       \f} </li>
!> <li>  Note: The location of the coordinate axis, i.e. the \f$R_{j,0}\f$ and \f$Z_{j,0}\f$,
!>       is set in the coordinate "packing" and "unpacking" routine, packxi(). </li>
!> </ul>
!>       **Jacobian**
!> <ul>
!> <li>  The coordinate Jacobian (and some other metric information) is given by
!>       <ul>
!>       <li> \c Igeometry=1 : Cartesian
!>       \f{eqnarray}{     {\bf e}_\theta \times {\bf e}_\zeta & = & -r_{tor}R_\theta \; \hat {\bf i} -  r_{pol}R_\zeta \; \hat {\bf j} + r_{pol}r_{tor}\hat {\bf k} \\
!>           \boldsymbol{\xi} \cdot {\bf e}_\theta \times {\bf e}_\zeta & = &  \delta R \\
!>           \sqrt g                                           & = &         R_s \ r_{pol} \ r_{tor}
!>       \f} </li>
!>       <li>  \c Igeometry=2 : Cylindrical
!>       \f{eqnarray}{    {\bf e}_\theta \times {\bf e}_\zeta & = & (R_\theta \sin \theta + R \cos\theta ) \; {\bf \hat i} + (R    \sin \theta - R_\theta \cos\theta ) \; {\bf \hat j} - R R_\zeta \; {\bf \hat k} \\
!>           \boldsymbol{\xi}\cdot {\bf e}_\theta \times {\bf e}_\zeta & = & \delta R \; R \\
!>           \sqrt g                                          & = & R_s \; R
!>       \f} </li>
!>       <li>  \c Igeometry=3 : Toroidal
!>       \f{eqnarray}{    {\bf e}_\theta \times {\bf e}_\zeta & = & - R \, Z_\theta \, \hat r + (Z_\theta \,R_\zeta - R_\theta \,Z_\zeta) \hat \phi + R \,R_\theta \,\hat z\\
!>           \boldsymbol{\xi}\cdot {\bf e}_\theta \times {\bf e}_\zeta & = & R ( \delta Z \; R_\theta - \delta R \; Z_\theta ) \\
!>            \sqrt g                                         & = & R ( Z_s      \; R_\theta - R_s      \; Z_\theta )
!>       \f} </li>
!>       </ul>
!>       </li>
!>
!> </ul>
!>       **cartesian metrics**
!> <ul>
!> <li>  The cartesian metrics are
!>       \f{eqnarray}{ \begin{array}{cccccccccccccccccccccccccccccccccccccccccccccccc}
!>           g_{s     s     } =  R_s      R_s                  ,&
!>           g_{s     \theta} =  R_s      R_\theta             ,&
!>           g_{s     \zeta } =  R_s      R_\zeta              ,&
!>           g_{\theta\theta} =  R_\theta R_\theta + r_{pol}^2 ,&
!>           g_{\theta\zeta } =  R_\theta R_\zeta              ,&
!>           g_{\zeta \zeta } =  R_\zeta  R_\zeta  + r_{tor}^2
!>           \end{array}
!>       \f} </li>
!>
!> </ul>
!>       **cylindrical metrics**
!> <ul>
!> <li>  The cylindrical metrics are
!>       \f{eqnarray}{ \begin{array}{cccccccccccccccccccccccccccccccccccccccccccccccc}
!>           g_{s     s     }  =  R_s      R_s            ,&
!>           g_{s     \theta}  =  R_s      R_\theta       ,&
!>           g_{s     \zeta }  =  R_s      R_\zeta        ,&
!>           g_{\theta\theta}  =  R_\theta R_\theta + R^2 ,&
!>           g_{\theta\zeta }  =  R_\theta R_\zeta        ,&
!>           g_{\zeta \zeta }  =  R_\zeta  R_\zeta  + 1
!>           \end{array}
!>       \f} </li>
!>
!> </ul>
!>       **logical control**
!> <ul>
!> <li>  The logical control is provided by \c Lcurvature as follows:
!>       <ul>
!>       <li> \c Lcurvature=0 : only the coordinate transformation is computed, i.e. only \f$R\f$ and \f$Z\f$ are calculated,          e.g. global() </li>
!>       <li> \c Lcurvature=1 : the Jacobian, \f$\sqrt g \f$, and "lower" metrics, \f$g_{\mu,\nu}\f$, are calculated ,                 e.g. bnorml(), lforce(), curent(), metrix(), sc00aa() </li>
!>       <li> \c Lcurvature=2 : the "curvature" terms are calculated, by which I mean the second derivatives of the position vector;
!>                              this information is required for computing the current, \f${\bf j}=\nabla\times\nabla\times{\bf A}\f$, e.g. jo00aa() </li>
!>       <li> \c Lcurvature=3 : the derivative of the \f$g_{\mu,\nu}/\sqrt g\f$ w.r.t. the interface boundary geometry is calculated,  e.g. metrix(), curent() </li>
!>       <li> \c Lcurvature=4 : the derivative of the \f$g_{\mu,\nu}\f$ w.r.t. the interface boundary geometry is calculated,          e.g. dforce() </li>
!>       <li> \c Lcurvature=5 : the derivative of \f$\sqrt{g}\f$ w.r.t. the interface boundary geometry is calculated,                 e.g. rzaxis() </li>
!>       </ul> </li>
!> </ul>
!>
!> @param[in] lvol specified in which volume to compute coordinates
!> @param[in] lss radial coordinate \f$s\f$
!> @param[in] Lcurvature logical control flag
!> @param[in] Ntz number of points in \f$\theta\f$ and \f$\zeta\f$
!> @param[in] mn number of Fourier harmonics
subroutine coords( lvol, lss, Lcurvature, Ntz, mn )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  use constants, only : zero, half, one, two, pi2

  use numerical, only : vsmall, small

  use fileunits, only : ounit

  use inputlist, only : Wcoords, Igeometry, Ntor, rpol, rtor

  use cputiming, only : Tcoords

  use allglobal, only : myid, cpus, pi2nfp, MPI_COMM_SPEC, &
                        Mvol, im, in, halfmm, &
                        iRbc, iZbs, iRbs, iZbc, &
                        NOTstellsym, Lcoordinatesingularity, &
                        Nt, Nz, &
                        Rij, Zij, &
                        cosi, sini, &
                        sg, guvij, &
                        dBdX, &
                        dRodR, dRodZ, dZodR, dZodZ

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS

  INTEGER, intent(in) :: lvol, Lcurvature, Ntz, mn
  REAL   , intent(in) :: lss

  INTEGER             :: ii, jj, kk, irz, innout, issym, signlss, mi, ni, imn
  REAL                :: Remn(1:mn,0:2), Zomn(1:mn,0:2), Romn(1:mn,0:2), Zemn(1:mn,0:2), alss, blss, sbar, sbarhim(1:mn), fj(1:mn,0:2)

  REAL                :: Dij(1:Ntz,0:3), dguvij(1:Ntz,1:3,1:3), DRxij(1:Ntz,0:3), DZxij(1:Ntz,0:3)

  BEGIN(coords)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef DEBUG
  FATAL( coords, lvol.lt.1 .or. lvol.gt.Mvol, invalid volume label )
  FATAL( coords, abs(lss).gt.one, invalid radial coordinate )
  FATAL( coords, Lcurvature.lt.0 .or. Lcurvature.gt.5, invalid input value for Lcurvature )
#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  Rij(1:Ntz,0:3,0:3) = zero ; sg(1:Ntz,0:3) = zero ; guvij(1:Ntz,1:3,1:3,0:3) = zero ! provide trivial default for output; 16 Jan 13;
  Zij(1:Ntz,0:3,0:3) = zero                                                          ! provide trivial default for output; 16 Jan 13;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  Remn(1:mn,0:2) = zero ! interpolated coordinate harmonics; 6 Feb 13;
  Zomn(1:mn,0:2) = zero
  Romn(1:mn,0:2) = zero
  Zemn(1:mn,0:2) = zero

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( Lcoordinatesingularity ) then

   sbar = ( lss + one ) * half

#ifdef DEBUG
   FATAL( coords, sbar.lt.zero .or. sbar.gt.one, invalid sbar )
#endif
   select case( Igeometry )
   case( 2   )  ; fj(     1:Ntor+1,0) = sbar                    ! these are the mj.eq.0 harmonics; 11 Aug 14; switch to sbar=r; 29 Jun 19
   ;            ; fj(Ntor+2:mn    ,0) = sbar**(im(Ntor+2:mn)+1) ! these are the me.ne.0 harmonics; 11 Aug 14; switch to sbar=r; 29 Jun 19
   case( 3   )  ; fj(     1:Ntor+1,0) = sbar**2                 ! switch to sbar=r; 29 Jun 19
   ;            ; fj(Ntor+2:mn    ,0) = sbar**im(Ntor+2:mn)     ! these are the me.ne.0 harmonics; 11 Aug 14; switch to sbar=r; 29 Jun 19
   case default ; FATAL( coords, .true., invalid Igeometry for Lcoordinatesingularity=T )
   end select


   Remn(1:mn,0) = iRbc(1:mn,0) + ( iRbc(1:mn,1) - iRbc(1:mn,0) ) * fj(1:mn,0)
   if( NOTstellsym ) then
   Romn(1:mn,0) = iRbs(1:mn,0) + ( iRbs(1:mn,1) - iRbs(1:mn,0) ) * fj(1:mn,0)
   endif
   if( Igeometry.eq.3 ) then
   Zomn(1:mn,0) = iZbs(1:mn,0) + ( iZbs(1:mn,1) - iZbs(1:mn,0) ) * fj(1:mn,0)
   if( NOTstellsym ) then
   Zemn(1:mn,0) = iZbc(1:mn,0) + ( iZbc(1:mn,1) - iZbc(1:mn,0) ) * fj(1:mn,0)
   endif
   endif

  else ! matches if( Lcoordinatesingularity ) ; 22 Apr 13;

   alss = half * ( one - lss ) ; blss = half * ( one + lss )

   Remn(1:mn,0) = alss * iRbc(1:mn,lvol-1) + blss * iRbc(1:mn,lvol)
   if( NOTstellsym ) then
   Romn(1:mn,0) = alss * iRbs(1:mn,lvol-1) + blss * iRbs(1:mn,lvol)
   endif
   if( Igeometry.eq.3 ) then
   Zomn(1:mn,0) = alss * iZbs(1:mn,lvol-1) + blss * iZbs(1:mn,lvol)
   if( NOTstellsym ) then
   Zemn(1:mn,0) = alss * iZbc(1:mn,lvol-1) + blss * iZbc(1:mn,lvol)
   endif
   endif ! end of if( Igeometry.eq.3 ) ; 01 Feb 13;

  endif ! end of if( Lcoordinatesingularity ); 01 Feb 13;

  call invfft( mn, im(1:mn), in(1:mn), Remn(1:mn,0), Romn(1:mn,0), Zemn(1:mn,0), Zomn(1:mn,0), &
               Nt, Nz, Rij(1:Ntz,0,0), Zij(1:Ntz,0,0) ) ! maps to real space;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( Lcurvature.eq.0 ) goto 9999 ! only the coordinates are required;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( Lcoordinatesingularity ) then

#ifdef DEBUG
   FATAL( coords, sbar.lt.small, small denominator )
#endif
   select case( Igeometry )
   case( 2   )  ; fj(     1:Ntor+1,1) = half                                                  ! these are the mj.eq.0 harmonics; 11 Aug 14; switch to sbar=r; 29 Jun 19
   ;            ; fj(Ntor+2:mn    ,1) = half*(im(Ntor+2:mn)+one) * fj(Ntor+2:mn    ,0) / sbar ! these are the me.ne.0 harmonics; 11 Aug 14; switch to sbar=r; 29 Jun 19
   case( 3   )  ; fj(     1:Ntor+1,1) = sbar                                                  ! these are the mj.eq.0 harmonics; 11 Aug 14; switch to sbar=r; 29 Jun 19
   ;            ; fj(Ntor+2:mn    ,1) = half * im(Ntor+2:mn) * fj(Ntor+2:mn    ,0) / sbar     ! these are the me.ne.0 harmonics; 11 Aug 14; switch to sbar=r; 29 Jun 19
   case default ; FATAL( coords, .true., invalid Igeometry for Lcoordinatesingularity=T and Lcurvature.ne.0 )
   end select

   Remn(1:mn,1) =                       ( iRbc(1:mn,1) - iRbc(1:mn,0) ) * fj(1:mn,1)
   if( NOTstellsym ) then
   Romn(1:mn,1) =                       ( iRbs(1:mn,1) - iRbs(1:mn,0) ) * fj(1:mn,1)
   endif
   if( Igeometry.eq.3 ) then
   Zomn(1:mn,1) =                       ( iZbs(1:mn,1) - iZbs(1:mn,0) ) * fj(1:mn,1)
   if( NOTstellsym ) then
   Zemn(1:mn,1) =                       ( iZbc(1:mn,1) - iZbc(1:mn,0) ) * fj(1:mn,1)
   endif
   endif

  else ! matches if( Lcoordinatesingularity ) ; 22 Apr 13;

   Remn(1:mn,1) = (      - iRbc(1:mn,lvol-1) +        iRbc(1:mn,lvol) ) * half
   if( NOTstellsym ) then
   Romn(1:mn,1) = (      - iRbs(1:mn,lvol-1) +        iRbs(1:mn,lvol) ) * half
   endif
   if( Igeometry.eq.3 ) then
   Zomn(1:mn,1) = (      - iZbs(1:mn,lvol-1) +        iZbs(1:mn,lvol) ) * half
   if( NOTstellsym ) then
   Zemn(1:mn,1) = (      - iZbc(1:mn,lvol-1) +        iZbc(1:mn,lvol) ) * half
   endif
   endif ! end of if( Igeometry.eq.3 ) ; 01 Feb 13;

  endif ! end of if( Lcoordinatesingularity ); 01 Feb 13;

  ! Derivative w.r.t s
  call invfft( mn, im(1:mn), in(1:mn),           Remn(1:mn,1),           Romn(1:mn,1),           Zemn(1:mn,1),           Zomn(1:mn,1), &
               Nt, Nz, Rij(1:Ntz,1,0), Zij(1:Ntz,1,0) ) ! maps to real space;

  ! Derivative w.r.t theta
  call invfft( mn, im(1:mn), in(1:mn),  im(1:mn)*Romn(1:mn,0), -im(1:mn)*Remn(1:mn,0),  im(1:mn)*Zomn(1:mn,0), -im(1:mn)*Zemn(1:mn,0), &
               Nt, Nz, Rij(1:Ntz,2,0), Zij(1:Ntz,2,0) ) ! maps to real space;

  ! Derivative w.r.t zeta
  call invfft( mn, im(1:mn), in(1:mn), -in(1:mn)*Romn(1:mn,0),  in(1:mn)*Remn(1:mn,0), -in(1:mn)*Zomn(1:mn,0),  in(1:mn)*Zemn(1:mn,0), &
               Nt, Nz, Rij(1:Ntz,3,0), Zij(1:Ntz,3,0) ) ! maps to real space;


  do ii = 1, 3 ; Rij(1:Ntz,0,ii) = Rij(1:Ntz,ii,0) ! just to complete workspace arrays; 22 Apr 13;
   ;           ; Zij(1:Ntz,0,ii) = Zij(1:Ntz,ii,0) ! just to complete workspace arrays; 22 Apr 13;
  enddo


  guvij(1:Ntz, 0, 0, 0) = one ! this is (only) required for the helicity integral; 22 Apr 13;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  select case( Igeometry )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  case( 1 ) ! Igeometry=1; Cartesian;

   sg(1:Ntz,0) = Rij(1:Ntz,1,0)*rpol*rtor

   do ii = 1, 3
    do jj = ii, 3 ; guvij(1:Ntz,ii,jj,0) = Rij(1:Ntz,ii,0) * Rij(1:Ntz,jj,0)
    enddo
   enddo

   guvij(1:Ntz, 2, 2,0) = guvij(1:Ntz, 2, 2,0) + rpol*rpol
   guvij(1:Ntz, 3, 3,0) = guvij(1:Ntz, 3, 3,0) + rtor*rtor

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  case( 2 ) ! Igeometry=2; cylindrical;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   sg(1:Ntz,0) = Rij(1:Ntz, 1,0) * Rij(1:Ntz, 0,0)

   do ii = 1, 3
    do jj = ii, 3 ; guvij(1:Ntz,ii,jj,0) = Rij(1:Ntz,ii,0) * Rij(1:Ntz,jj,0)
    enddo
   enddo

   guvij(1:Ntz, 2, 2,0) = guvij(1:Ntz, 2, 2,0) + Rij(1:Ntz, 0, 0) * Rij(1:Ntz, 0, 0)
   guvij(1:Ntz, 3, 3,0) = guvij(1:Ntz, 3, 3,0) + one

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  case( 3 ) ! Igeometry=3; toroidal;

   sg(1:Ntz,0) = Rij(1:Ntz,0,0) * ( Zij(1:Ntz,1,0)*Rij(1:Ntz,2,0) - Rij(1:Ntz,1,0)*Zij(1:Ntz,2,0) )

   do ii = 1, 3
    do jj = ii, 3 ; guvij(1:Ntz,ii,jj,0) = Rij(1:Ntz,ii,0) * Rij(1:Ntz,jj,0) + Zij(1:Ntz,ii,0) * Zij(1:Ntz,jj,0)
    enddo
   enddo

   guvij(1:Ntz, 3, 3,0) = guvij(1:Ntz, 3, 3,0) + ( Rij(1:Ntz,0,0) )**2

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  case default ! Igeometry; 09 Mar 17;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   FATAL( coords, .true., selected Igeometry not supported )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  end select ! end of select case( Igeometry ) ; 15 Sep 16;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  do ii = 2, 3
   do jj = 1, ii-1 ; guvij(1:Ntz,ii,jj,0) = guvij(1:Ntz,jj,ii,0) ! complete metric array; 20 Apr 13;
   enddo
  enddo

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( Lcurvature.le.1 ) goto 9999

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  select case( Lcurvature )

  case( 2 ) ! Lcurvature=2; get second derivatives of position wrt \s, \t & \z; 19 Sep 13;

   if( Lcoordinatesingularity ) then

#ifdef DEBUG
    FATAL( coords, sbar.lt.small, small denominator )
#endif
    select case( Igeometry )
    case( 2 )    ; fj(     1:Ntor+1,2) = zero                                                            ! these are the mj.eq.0 harmonics; 11 Aug 14; switch to sbar=r; 29 Jun 19
     ;           ; fj(Ntor+2:mn    ,2) = half * ( im(Ntor+2:mn)       ) * fj(Ntor+2:mn    ,1) / sbar     ! these are the me.ne.0 harmonics; 11 Aug 14; switch to sbar=r; 29 Jun 19
    case( 3 )    ; fj(     1:Ntor+1,2) = half                                                            ! these are the mj.eq.0 harmonics; 11 Aug 14; switch to sbar=r; 29 Jun 19
     ;           ; fj(Ntor+2:mn    ,2) = half * ( im(Ntor+2:mn) - one ) * fj(Ntor+2:mn    ,1) / sbar     ! these are the me.ne.0 harmonics; 11 Aug 14; switch to sbar=r; 29 Jun 19
    case default ;
     ;           ; FATAL( coords, .true., invalid Igeometry for Lcoordinatesingularity=T and Lcurvature=2 )
    end select   ;

    Remn(1:mn,2) =                       ( iRbc(1:mn,1) - iRbc(1:mn,0) ) * fj(1:mn,2)
    if( NOTstellsym ) then
    Romn(1:mn,2) =                       ( iRbs(1:mn,1) - iRbs(1:mn,0) ) * fj(1:mn,2)
    endif
    if( Igeometry.eq.3 ) then
    Zomn(1:mn,2) =                       ( iZbs(1:mn,1) - iZbs(1:mn,0) ) * fj(1:mn,2)
    if( NOTstellsym ) then
    Zemn(1:mn,2) =                       ( iZbc(1:mn,1) - iZbc(1:mn,0) ) * fj(1:mn,2)
    endif
    endif

   else ! if( .not.Lcoordinatesingularity ) ;

    Remn(1:mn,2) = zero
    if( NOTstellsym ) then
    Romn(1:mn,2) = zero
    endif
    if( Igeometry.eq.3 ) then
    Zomn(1:mn,2) = zero
    if( NOTstellsym ) then
    Zemn(1:mn,2) = zero
    endif
    endif ! end of if( Igeometry.eq.3 ) ; 01 Feb 13;

   endif ! end of if( Lcoordinatesingularity ); 01 Feb 13;

   call invfft( mn, im(1:mn), in(1:mn),&
                   Remn(1:mn,2),                   Romn(1:mn,2),                   Zemn(1:mn,2),                   Zomn(1:mn,2), &
                Nt, Nz, Rij(1:Ntz,1,1), Zij(1:Ntz,1,1) ) ! maps to real space;

   call invfft( mn, im(1:mn), in(1:mn),&
+         im(1:mn)*Romn(1:mn,1),         -im(1:mn)*Remn(1:mn,1),          im(1:mn)*Zomn(1:mn,1),         -im(1:mn)*Zemn(1:mn,1), &
Nt, Nz, Rij(1:Ntz,1,2), Zij(1:Ntz,1,2) ) ! maps to real space;

   call invfft( mn, im(1:mn), in(1:mn),&
-         in(1:mn)*Romn(1:mn,1),          in(1:mn)*Remn(1:mn,1),         -in(1:mn)*Zomn(1:mn,1),          in(1:mn)*Zemn(1:mn,1), &
Nt, Nz, Rij(1:Ntz,1,3), Zij(1:Ntz,1,3) ) ! maps to real space;

   call invfft( mn, im(1:mn), in(1:mn),&
-im(1:mn)*im(1:mn)*Remn(1:mn,0),-im(1:mn)*im(1:mn)*Romn(1:mn,0),-im(1:mn)*im(1:mn)*Zemn(1:mn,0),-im(1:mn)*im(1:mn)*Zomn(1:mn,0), &
Nt, Nz, Rij(1:Ntz,2,2), Zij(1:Ntz,2,2) ) ! maps to real space;

   call invfft( mn, im(1:mn), in(1:mn),&
+im(1:mn)*in(1:mn)*Remn(1:mn,0), im(1:mn)*in(1:mn)*Romn(1:mn,0), im(1:mn)*in(1:mn)*Zemn(1:mn,0), im(1:mn)*in(1:mn)*Zomn(1:mn,0), &
Nt, Nz, Rij(1:Ntz,2,3), Zij(1:Ntz,2,3) ) ! maps to real space;

   call invfft( mn, im(1:mn), in(1:mn),&
-in(1:mn)*in(1:mn)*Remn(1:mn,0),-in(1:mn)*in(1:mn)*Romn(1:mn,0),-in(1:mn)*in(1:mn)*Zemn(1:mn,0),-in(1:mn)*in(1:mn)*Zomn(1:mn,0), &
Nt, Nz, Rij(1:Ntz,3,3), Zij(1:Ntz,3,3) ) ! maps to real space;


   do ii = 2, 3
    do jj = 1, ii-1 ; Rij(1:Ntz,ii,jj) = Rij(1:Ntz,jj,ii) ; Zij(1:Ntz,ii,jj) = Zij(1:Ntz,jj,ii)
    enddo
   enddo


   select case( Igeometry )

   case( 1 ) ! Lcurvature=2; Igeometry=1 ; Cartesian;

    do kk = 1, 3 ! kk labels derivative; 13 Sep 13;

!    sg(1:Ntz, 0) = Rij(1:Ntz,1, 0)*rpol*rtor
     sg(1:Ntz,kk) = Rij(1:Ntz,1,kk)*rpol*rtor

     do ii = 1, 3
      do jj = ii, 3 ; guvij(1:Ntz,ii,jj,kk) = Rij(1:Ntz,ii,kk) * Rij(1:Ntz,jj, 0) + Rij(1:Ntz,ii, 0) * Rij(1:Ntz,jj,kk)
      enddo
     enddo

    enddo ! 06 Feb 13;

   case( 2 ) ! Lcurvature=2; Igeometry=2 ; cylindrical;

    do kk = 1, 3 ! kk labels derivative; 13 Sep 13;

!    sg(1:Ntz, 0) = Rij(1:Ntz, 1, 0) * Rij(1:Ntz, 0,0)
     sg(1:Ntz,kk) = Rij(1:Ntz, 1,kk) * Rij(1:Ntz, 0,0) + Rij(1:Ntz, 1,0) * Rij(1:Ntz, 0,kk)

     do ii = 1, 3
      do jj = ii, 3 ; guvij(1:Ntz,ii,jj,kk) = Rij(1:Ntz,ii,kk) * Rij(1:Ntz,jj, 0) + Rij(1:Ntz,ii, 0) * Rij(1:Ntz,jj,kk)
      enddo
     enddo

     guvij(1:Ntz,2,2,kk) = guvij(1:Ntz,2,2,kk) + Rij(1:Ntz,0,kk) * Rij(1:Ntz,0,0) + Rij(1:Ntz,0,0) * Rij(1:Ntz,0,kk) ! 20 Jan 15;
    !guvij(1:Ntz,2,2,kk) = guvij(1:Ntz,2,2,kk) + Rij(1:Ntz,kk,0) * Rij(1:Ntz,0,0) + Rij(1:Ntz,0,0) * Rij(1:Ntz,kk,0) ! 20 Jan 15;
    !guvij(1:Ntz,3,3,kk) = additional term is unity;

    enddo ! 06 Feb 13;

   case( 3 ) ! Lcurvature=2; Igeometry=3 ; toroidal;

    do kk = 1 , 3 ! kk labels derivative; 13 Sep 13;

!    sg(1:Ntz, 0) = Rij(1:Ntz, 0,0) * ( Zij(1:Ntz,1, 0)*Rij(1:Ntz,2, 0) - Rij(1:Ntz,1, 0)*Zij(1:Ntz,2, 0) )
     sg(1:Ntz,kk) = Rij(1:Ntz,kk,0) * ( Zij(1:Ntz,1, 0)*Rij(1:Ntz,2, 0) - Rij(1:Ntz,1, 0)*Zij(1:Ntz,2, 0) ) &
                  + Rij(1:Ntz, 0,0) * ( Zij(1:Ntz,1,kk)*Rij(1:Ntz,2, 0) - Rij(1:Ntz,1,kk)*Zij(1:Ntz,2, 0) ) &
                  + Rij(1:Ntz, 0,0) * ( Zij(1:Ntz,1, 0)*Rij(1:Ntz,2,kk) - Rij(1:Ntz,1, 0)*Zij(1:Ntz,2,kk) )

     sg(1:Ntz,kk) = sg(1:Ntz,kk)

     do ii = 1, 3
      do jj = ii, 3
       guvij(1:Ntz,ii,jj,kk) = Rij(1:Ntz,ii,kk) * Rij(1:Ntz,jj, 0) &
                             + Rij(1:Ntz,ii, 0) * Rij(1:Ntz,jj,kk) &
                             + Zij(1:Ntz,ii,kk) * Zij(1:Ntz,jj, 0) &
                             + Zij(1:Ntz,ii, 0) * Zij(1:Ntz,jj,kk)
      enddo
     enddo

     guvij(1:Ntz,3,3,kk) = guvij(1:Ntz,3,3,kk) + ( Rij(1:Ntz,0,kk) * Rij(1:Ntz,0,0) + Rij(1:Ntz,0,0) * Rij(1:Ntz,0,kk) )

    enddo ! end of do kk;  5 Feb 13;

   case default

    FATAL( coords, .true., selected Igeometry not supported for Lcurvature.eq.2 )

   end select ! end of select case( Igeometry ) ; 15 Sep 16;

   do ii = 2, 3
    do jj = 1, ii-1 ; guvij(1:Ntz,ii,jj,1:3) = guvij(1:Ntz,jj,ii,1:3)
    enddo
   enddo


  case( 3,4,5 ) ! Lcurvature=3,4,5 ; get derivatives wrt R_j and Z_j; 19 Sep 13;


   ii = dBdX%ii ; innout = dBdX%innout ; irz = dBdX%irz ; issym = dBdX%issym ! shorthand;

   if( Lcoordinatesingularity ) then

#ifdef DEBUG
    FATAL( coords, innout.eq.0, cannot differentiate metric elements wrt coordinate singularity )
#endif

#ifdef DEBUG
    FATAL( coords, Igeometry.eq.1, Cartesian does not need regularization factor )
#endif

    if( ( irz.eq.0 .and. issym.eq.0 ) .or. ( irz.eq.1 .and. issym.eq.1 ) ) then     ! cosine harmonics; 13 Sep 13;

     Dij(1:Ntz,0) = fj(ii,0) * cosi(1:Ntz,ii)
     Dij(1:Ntz,1) = fj(ii,1) * cosi(1:Ntz,ii)
     Dij(1:Ntz,2) = fj(ii,0) * sini(1:Ntz,ii) * ( - im(ii) )
     Dij(1:Ntz,3) = fj(ii,0) * sini(1:Ntz,ii) * ( + in(ii) )

    else                                                                            !   sine harmonics; 13 Sep 13;

     Dij(1:Ntz,0) = fj(ii,0) * sini(1:Ntz,ii)
     Dij(1:Ntz,1) = fj(ii,1) * sini(1:Ntz,ii)
     Dij(1:Ntz,2) = fj(ii,0) * cosi(1:Ntz,ii) * ( + im(ii) )
     Dij(1:Ntz,3) = fj(ii,0) * cosi(1:Ntz,ii) * ( - in(ii) )

    endif ! if( ( irz.eq.0 .and. issym.eq.1 ) .or. ... ; 11 Aug 14;

    if (Lcurvature .eq. 3 .or. Lcurvature .eq. 4) then
      if (Igeometry .eq. 3) then
        ! add the terms due to the moving coordinate axis
        if ( irz.eq.0 ) then
          DRxij(1:Ntz,0) = (one-sbar**2) * dRodR(1:Ntz,issym,ii)      ! dRx/dR
          DRxij(1:Ntz,1) = - sbar * dRodR(1:Ntz,issym,ii)
          DRxij(1:Ntz,2) = zero
          DRxij(1:Ntz,3) = (one-sbar**2) * dRodR(1:Ntz,issym+2,ii)    ! dRx/dR, zeta derivative

          DZxij(1:Ntz,0) = (one-sbar**2) * dZodR(1:Ntz,issym,ii)      ! dZx/dR
          DZxij(1:Ntz,1) = - sbar * dZodR(1:Ntz,issym,ii)
          DZxij(1:Ntz,2) = zero
          DZxij(1:Ntz,3) = (one-sbar**2) * dZodR(1:Ntz,issym+2,ii)    ! dZx/dR, zeta derivative
        else
          DRxij(1:Ntz,0) = (one-sbar**2) * dRodZ(1:Ntz,1-issym,ii)    ! dRx/dZ
          DRxij(1:Ntz,1) = - sbar  * dRodZ(1:Ntz,1-issym,ii)
          DRxij(1:Ntz,2) = zero
          DRxij(1:Ntz,3) = (one-sbar**2) * dRodZ(1:Ntz,1-issym+2,ii)  ! dRx/dZ, zeta derivative

          DZxij(1:Ntz,0) = (one-sbar**2) * dZodZ(1:Ntz,1-issym,ii)    ! dZx/dZ
          DZxij(1:Ntz,1) = - sbar  * dZodZ(1:Ntz,1-issym,ii)
          DZxij(1:Ntz,2) = zero
          DZxij(1:Ntz,3) = (one-sbar**2) * dZodZ(1:Ntz,1-issym+2,ii)  ! dZx/dZ, zeta derivative
        endif
      endif
    endif
    !DRxij(:,:) = zero
    !DZxij(:,:) = zero

   else ! matches if( Lcoordinatesingularity ) ; 10 Mar 13;

    if( innout.eq.0 ) signlss = - 1
    if( innout.eq.1 ) signlss = + 1

    if( ( irz.eq.0 .and. issym.eq.0 ) .or. ( irz.eq.1 .and. issym.eq.1 ) ) then ! cosine; 02 Sep 14;
     Dij(1:Ntz,0) = ( one + signlss * lss ) * half * cosi(1:Ntz,ii)
     Dij(1:Ntz,1) = (       signlss       ) * half * cosi(1:Ntz,ii)
     Dij(1:Ntz,2) = ( one + signlss * lss ) * half * sini(1:Ntz,ii) * ( - im(ii) )
     Dij(1:Ntz,3) = ( one + signlss * lss ) * half * sini(1:Ntz,ii) * ( + in(ii) )
    else                                                                        !   sine; 02 Sep 14;
     Dij(1:Ntz,0) = ( one + signlss * lss ) * half * sini(1:Ntz,ii)
     Dij(1:Ntz,1) = (       signlss       ) * half * sini(1:Ntz,ii)
     Dij(1:Ntz,2) = ( one + signlss * lss ) * half * cosi(1:Ntz,ii) * ( + im(ii) )
     Dij(1:Ntz,3) = ( one + signlss * lss ) * half * cosi(1:Ntz,ii) * ( - in(ii) )
    endif

   endif ! end of if( Lcoordinatesingularity ) ;  7 Mar 13;

   if (Lcurvature .eq. 5) then ! we only need the 2D Jacobian
    if (Igeometry .eq. 3) then ! only works for toroidal

      if( irz.eq.0 ) sg(1:Ntz,1) = ( Zij(1:Ntz,1,0)*Dij(1:Ntz,2  ) - Dij(1:Ntz,1  )*Zij(1:Ntz,2,0) )
      if( irz.eq.1 ) sg(1:Ntz,1) = ( Dij(1:Ntz,1  )*Rij(1:Ntz,2,0) - Rij(1:Ntz,1,0)*Dij(1:Ntz,2  ) )

    else
      FATAL( coords, Lcurvature.eq.5 .and. Igeometry.ne.3, Lcurvature.eq.5 can only be combined with Igeometry.ne.3 )
    end if ! if (Igeometry .eq. 3) ; 13 Jan 20

   else ! we need more for Lcurvature=3,4 ; 13 Jan 20

    select case( Igeometry )

    case( 1 ) ! Lcurvature=3,4 ; Igeometry=1 ; Cartesian; 04 Dec 14;

#ifdef DEBUG
    FATAL( coords, irz.eq.1, there is no dependence on Zbs or Zbc )
#endif

!                  sg(1:Ntz,0) = Rij(1:Ntz,1,0)*rpol*rtor
                   sg(1:Ntz,1) = Dij(1:Ntz,1  )*rpol*rtor ! 20 Jun 14; 29 Apr 19;
!   if( irz.eq.1 ) sg(1:Ntz,1) =

    do ii = 1, 3 ! careful: ii was used with a different definition above; 13 Sep 13;
     do jj = ii, 3
                     dguvij(1:Ntz,ii,jj) = Dij(1:Ntz,ii) * Rij(1:Ntz,jj,0) + Rij(1:Ntz,ii,0) * Dij(1:Ntz,jj)
!     if( irz.eq.1 ) dguvij(1:Ntz,ii,jj) =
     enddo
    enddo

   case( 2 ) ! Lcurvature=3,4,5 ; Igeometry=2 ; cylindrical;

#ifdef DEBUG
      FATAL( coords, irz.eq.1, there is no dependence on Zbs or Zbc )
#endif

  !                  sg(1:Ntz,0) = Rij(1:Ntz,1,0) * Rij(1:Ntz,0,0)
      if( irz.eq.0 ) sg(1:Ntz,1) = Dij(1:Ntz,1  ) * Rij(1:Ntz,0,0) &
                                + Rij(1:Ntz,1,0) * Dij(1:Ntz,0  )

    do ii = 1, 3 ! careful: ii was used with a different definition above; 13 Sep 13;
     do jj = ii, 3
      if( irz.eq.0 ) dguvij(1:Ntz,ii,jj) = Dij(1:Ntz,ii) * Rij(1:Ntz,jj,0) + Rij(1:Ntz,ii,0) * Dij(1:Ntz,jj)
      if( irz.eq.1 ) then
         FATAL(coords, .true., No Z-geometrical degree of freedom when Igeometry=2)!dguvij(1:Ntz,ii,jj) = Dij(1:Ntz,ii) * Zij(1:Ntz,jj,0) + Zij(1:Ntz,ii,0) * Dij(1:Ntz,jj) ! TODO REMOVE
      endif
     enddo
    enddo

    dguvij(1:Ntz,2,2) = dguvij(1:Ntz,2,2) + two * Dij(1:Ntz,0) * Rij(1:Ntz,0,0)
   !dguvij(1:Ntz,3,3) = additional term is unity;

    case( 3 ) ! Lcurvature=3,4,5 ; Igeometry=3 ; toroidal; 04 Dec 14;

      if (LcoordinateSingularity) then
    !                  sg(1:Ntz,0) = Rij(1:Ntz,0,0) * ( Zij(1:Ntz,1,0)*Rij(1:Ntz,2,0) - Rij(1:Ntz,1,0)*Zij(1:Ntz,2,0) )
        if( irz.eq.0 ) sg(1:Ntz,1) = (Dij(1:Ntz,0  )+ DRxij(1:Ntz, 0)) * ( Zij(1:Ntz,1,0)*Rij(1:Ntz,2,0) - Rij(1:Ntz,1,0)*Zij(1:Ntz,2,0) ) &
                                  + Rij(1:Ntz,0,0) * ( Zij(1:Ntz,1,0)*Dij(1:Ntz,2  ) - (Dij(1:Ntz,1  )+DRxij(1:Ntz,1)) *Zij(1:Ntz,2,0) ) &
                                  + Rij(1:Ntz,0,0) * ( DZxij(1:Ntz,1)*Rij(1:Ntz,2,0) )
        if( irz.eq.1 ) sg(1:Ntz,1) = DRxij(1:Ntz, 0) * ( Zij(1:Ntz,1,0)*Rij(1:Ntz,2,0) - Rij(1:Ntz,1,0)*Zij(1:Ntz,2,0) ) &
                                  + Rij(1:Ntz,0,0) * ( (Dij(1:Ntz,1  )+DZxij(1:Ntz,1)) *Rij(1:Ntz,2,0) - Rij(1:Ntz,1,0)*Dij(1:Ntz,2  ) ) &
                                  + Rij(1:Ntz,0,0) * (-DRxij(1:Ntz,1)*Zij(1:Ntz,2,0))

        do ii = 1, 3 ! careful: ii was used with a different definition above; 13 Sep 13;
          do jj = ii, 3
            if( irz.eq.0 ) dguvij(1:Ntz,ii,jj) = (Dij(1:Ntz,ii)+DRxij(1:Ntz,ii)) * Rij(1:Ntz,jj,0) + Rij(1:Ntz,ii,0) * (Dij(1:Ntz,jj)+DRxij(1:Ntz,jj)) &
                                              + DZxij(1:Ntz,ii) * Zij(1:Ntz,jj,0) + Zij(1:Ntz,ii,0) * DZxij(1:Ntz,jj)
            if( irz.eq.1 ) dguvij(1:Ntz,ii,jj) = (Dij(1:Ntz,ii)+DZxij(1:Ntz,ii)) * Zij(1:Ntz,jj,0) + Zij(1:Ntz,ii,0) * (Dij(1:Ntz,jj)+DZxij(1:Ntz,jj)) &
                                              + DRxij(1:Ntz,ii) * Rij(1:Ntz,jj,0) + Rij(1:Ntz,ii,0) * DRxij(1:Ntz,jj)
          enddo
        enddo

        if( irz.eq.0 ) dguvij(1:Ntz,3,3) = dguvij(1:Ntz,3,3) + two * (Dij(1:Ntz,0)+DRxij(1:Ntz,0)) * Rij(1:Ntz,0,0)
        if( irz.eq.1 ) dguvij(1:Ntz,3,3) = dguvij(1:Ntz,3,3) + two * (DRxij(1:Ntz,0)) * Rij(1:Ntz,0,0)

      else

        if( irz.eq.0 ) sg(1:Ntz,1) = Dij(1:Ntz,0  ) * ( Zij(1:Ntz,1,0)*Rij(1:Ntz,2,0) - Rij(1:Ntz,1,0)*Zij(1:Ntz,2,0) ) &
                                  + Rij(1:Ntz,0,0) * ( Zij(1:Ntz,1,0)*Dij(1:Ntz,2  ) - Dij(1:Ntz,1  )*Zij(1:Ntz,2,0) )
        if( irz.eq.1 ) sg(1:Ntz,1) = Rij(1:Ntz,0,0) * ( Dij(1:Ntz,1  )*Rij(1:Ntz,2,0) - Rij(1:Ntz,1,0)*Dij(1:Ntz,2  ) )


        do ii = 1, 3 ! careful: ii was used with a different definition above; 13 Sep 13;
          do jj = ii, 3
            if( irz.eq.0 ) dguvij(1:Ntz,ii,jj) = Dij(1:Ntz,ii) * Rij(1:Ntz,jj,0) + Rij(1:Ntz,ii,0) * Dij(1:Ntz,jj)
            if( irz.eq.1 ) dguvij(1:Ntz,ii,jj) = Dij(1:Ntz,ii) * Zij(1:Ntz,jj,0) + Zij(1:Ntz,ii,0) * Dij(1:Ntz,jj)
          enddo
        enddo

        if( irz.eq.0 ) dguvij(1:Ntz,3,3) = dguvij(1:Ntz,3,3) + two * Dij(1:Ntz,0) * Rij(1:Ntz,0,0)
      endif

    case default

      FATAL( coords, .true., supplied Igeometry is not yet supported for Lcurvature.eq.3 or Lcurvature.eq.4 )

    end select ! end of select case( Igeometry );  7 Mar 13;

    do ii = 2, 3
      do jj = 1, ii-1 ; dguvij(1:Ntz,ii,jj) = dguvij(1:Ntz,jj,ii) ! symmetry of metrics; 13 Sep 13;
      enddo
    enddo

    guvij(1:Ntz,0,0,1) = zero ! this "metric" does not depend on geometry; helicity matrix does not depend on geometry; 10 Mar 13;

    if( Lcurvature.eq.3 ) then

      do ii = 1, 3
      do jj = 1, 3 ; guvij(1:Ntz,ii,jj,1) = dguvij(1:Ntz,ii,jj) - guvij(1:Ntz,ii,jj,0) * sg(1:Ntz,1) / sg(1:Ntz,0) ! differentiated metric elements; 7 Mar 13;
      enddo
      enddo

    else ! if( Lcurvature.eq.4 ) ;

      do ii = 1, 3
      do jj = 1, 3 ; guvij(1:Ntz,ii,jj,1) = dguvij(1:Ntz,ii,jj)                                                    ! differentiated metric elements; 7 Mar 13;
      enddo
      enddo

    endif ! end of if( Lcurvature.eq.3 ) ; 15 Sep 16;

   endif ! end of if( Lcurvature.eq.5 ) ; 13 Jan 20;
  end select ! matches select case( Lcurvature ) ; 10 Mar 13;


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(coords)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine coords

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
