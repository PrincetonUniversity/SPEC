!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (geometry) ! Calculates coordinates, ${\bf x}(s,\theta,\zeta) \equiv R \, {\bf e}_R + Z \, {\bf e}_Z$, and metrics, using FFTs.

!latex \briefly{Calculates coordinate transformation, and metric elements and curvatures if required, using FFTs.}

!latex \calledby{\link{global}, \link{bnorml}, \link{lforce}, \link{dforce}, \link{curent}, \link{jo00aa}, \link{metrix}, \link{sc00aa}}
!latex \calls{\link{}}

!latex \tableofcontents

!latex \subsection{coordinates} 

!latex \begin{enumerate}

!latex \item We work in coordinates, $(\s,\t,\z)$, which are be defined {\em inversely} via a transformation {\em to} Cartesian coordinates, $(x,y,z)$.

!latex \item The toroidal angle, $\z$, is identical to the cylindrical angle, $\z\equiv\phi$.

!latex \item The radial coordinate, $s$, is {\em not} a global variable: it only needs to be defined in each volume, and in each volume $\s \in [-1,1]$.

!latex \item The choice of poloidal angle, $\t$, does not affect the following.

!latex \end{enumerate} \subsection{geometry} \begin{enumerate}

!latex \item The geometry of the ``ideal''-interfaces, ${\bf x}_v(\t,\z)$, is given by $R(\t,\z)$ and $Z(\t,\z)$ as follows:
!latex       \begin{itemize}
!latex       \item \verb+Igeometry=1+ : Cartesian
!latex       \be {\bf x} & \equiv & \t \; {\bf \hat i} + \z \; {\bf \hat j}+ R \; {\bf \hat k}
!latex       \ee
!latex       \item \verb+Igeometry=2+ : Cylindrical
!latex       \be {\bf x} & = & R \; \cos\t \; {\bf \hat i} + R \; \sin\t \; {\bf \hat j} + \z \; {\bf \hat k}
!latex       \ee
!latex       \item \verb+Igeometry=3+ : Toroidal
!latex       \be {\bf x} & \equiv & R \; {\bf \hat r} + Z \; {\bf \hat k}
!latex       \ee
!latex       where ${\bf \hat r}\equiv \cos \phi \; {\bf \hat i} + \sin \phi \; {\bf \hat j}$ and
!latex             ${\bf \hat \phi} \equiv - \sin \phi \; {\bf \hat i} + \cos \phi \; {\bf \hat j}$.
!latex       \end{itemize}

!latex \item The geometry of the ideal interfaces is given as Fourier summation: e.g., for stellarator-symmetry
!latex       \be R_v(\t,\z) & \equiv & \sum_j R_{j,v} \cos\alpha_j, \\
!latex           Z_v(\t,\z) & \equiv & \sum_j Z_{j,v} \sin\alpha_j,
!latex       \ee
!latex       where $\alpha_j \equiv m_j \t - n_j \z$.

!latex \end{enumerate} \subsection{interpolation between interfaces} \begin{enumerate}

!latex \item The ``coordinate" functions, $R(\s,\t,\z)$ and $Z(\s,\t,\z)$, are constructed by radially interpolating the Fourier representations of the 
!latex       ideal-interfaces.

!latex \item The $v$-th volume is bounded by ${\bf x}_{v-1}$ and ${\bf x}_{v}$.

!latex \item In each {\em annular} volume, the coordinates are constructed by linear interpolation:
!latex       \be \begin{array}{cccccccccccccccccccccccccccc}
!latex           R(\s,\t,\z) & \equiv & \ds \sum_j & \ds \left[ \; \frac{(1-s)}{2} \; R_{j,v-1} + \frac{(1+s)}{2} \; R_{j,v}\; \right] \; \cos\alpha_j,\\
!latex           Z(\s,\t,\z) & \equiv & \ds \sum_j & \ds \left[ \; \frac{(1-s)}{2} \; Z_{j,v-1} + \frac{(1+s)}{2} \; Z_{j,v}\; \right] \; \sin\alpha_j,
!latex           \end{array}
!latex       \ee

!latex \end{enumerate} \subsubsection{coordinate singularity: regularized extrapolation} \begin{enumerate}

!latex \item For cylindrical or toroidal geometry, in the innermost, ``simple-torus'' volume, the coordinates are constructed by an interpolation that
!latex       ``encourages'' the interpolated coordinate surfaces to not intersect.
!latex \item Introduce $\bar s \equiv (s+1)/2$, so that in each volume $\bar s \in [0,1]$, then
!latex       \be R_{j}(s) & = & R_{j,0} + (R_{j,1} - R_{j,0} ) f_j, \\
!latex           Z_{j}(s) & = & Z_{j,0} + (Z_{j,1} - Z_{j,0} ) f_j,
!latex       \ee
!latex       where, in toroidal geometry, 
!latex       \be
!latex       f_j \equiv \left\{ 
!latex       \begin{array}{llcccccccccccccc} \bar s        & , & \mbox{\rm for } m_j=0, \\
!latex                                       \bar s^{m_j/2}& , & \mbox{\rm otherwise.} 
!latex       \end{array}\right. 
!latex       \ee

!latex \item Note: The location of the coordinate axis, i.e. the $R_{j,0}$ and $Z_{j,0}$,
!latex       is set in the coordinate ``packing'' and ``unpacking'' routine, \link{packxi}.

!latex \end{enumerate} \subsection{Jacobian} \begin{enumerate}

!latex \item The coordinate Jacobian (and some other metric information) is given by
!latex       \begin{itemize}
!latex       \item \verb+Igeometry=1+ : Cartesian
!latex       \be {\bf e}_\theta \times {\bf e}_\zeta & = & -R_\t \; \hat {\bf i} -  R_\z \; \hat {\bf j} + \hat {\bf k} \\
!latex           \boldxi \cdot {\bf e}_\theta \times {\bf e}_\zeta & =&  \delta R \\
!latex           \sqrt g                                           & =&         R_s
!latex       \ee
!latex       \item \verb+Igeometry=2+ : Cylindrical
!latex       \be {\bf e}_\theta \times {\bf e}_\zeta & = & 
!latex       (R_\t \sin \t + R \cos\t ) \; {\bf \hat i} + (R    \sin \t - R_\t \cos\t ) \; {\bf \hat j} - R R_\z \; {\bf \hat k} \\
!latex           \boldxi\cdot {\bf e}_\theta \times {\bf e}_\zeta & = & \delta R \; R \\
!latex           \sqrt g                                          & = & R_s \; R
!latex       \ee
!latex       \item \verb+Igeometry=3+ : Toroidal
!latex       \be {\bf e}_\theta \times {\bf e}_\zeta & = & 
!latex       - R \, Z_\theta \, \hat r + (Z_\theta \,R_\zeta - R_\theta \,Z_\zeta) \hat \phi + R \,R_\theta \,\hat z\\
!latex           \boldxi\cdot {\bf e}_\theta \times {\bf e}_\zeta & = & R ( \delta Z \; R_\t - \delta R \; Z_\t ) \\
!latex           \sqrt g                                         & = & R ( Z_s      \; R_\t - R_s      \; Z_\t )
!latex       \ee
!latex       \end{itemize}

!latex \end{enumerate} \subsubsection{cylindrical metrics} \begin{enumerate}
!latex \item The cylindrical metrics and Jacobian are
!latex       \be \begin{array}{cccccccccccccccccccccccccccccccccccccccccccccccc}
!latex           \sqrt g   =  R_\s R          ,&
!latex           g_{\s\s}  =  R_\s R_\s       ,&
!latex           g_{\s\t}  =  R_\s R_\t       ,&
!latex           g_{\s\z}  =  R_\s R_\z       ,&
!latex           g_{\t\t}  =  R_\t R_\t + R^2 ,&
!latex           g_{\t\z}  =  R_\t R_\z       ,&
!latex           g_{\z\z}  =  R_\z R_\z + 1
!latex           \end{array}
!latex       \ee

!latex \end{enumerate} \subsection{logical control} \begin{enumerate}
!latex \item The logical control is provided by \type{Lcurvature} as follows:
!latex       \bi
!latex       \item[] \type{Lcurvature=0} : only the coordinate transformation is computed, i.e. only $R$ and $Z$ are calculated \\
!latex                                     e.g. \link{global}
!latex       \item[] \type{Lcurvature=1} : the Jacobian, $\sqrt g $, and ``lower'' metrics, $g_{\mu,\nu}$, are calculated \\
!latex                                     e.g. \link{bnorml}, \link{lforce}, \link{curent}, \link{metrix}, \link{sc00aa}
!latex       \item[] \type{Lcurvature=2} : the ``curvature'' terms are calculated, by which I mean the second derivatives of the position vector;
!latex                                     this information is required for computing the current, ${\bf j}=\nabla\times\nabla\times{\bf A}$ \\
!latex                                     e.g. \link{jo00aa}
!latex       \item[] \type{Lcurvature=3} : the derivative of the $g_{\mu,\nu}/\sqrt g$ w.r.t. the interface boundary geometry is calculated \\
!latex                                     e.g. \link{metrix}, \link{curent}
!latex       \item[] \type{Lcurvature=4} : the derivative of the $g_{\mu,\nu}$ w.r.t. the interface boundary geometry is calculated \\
!latex                                     e.g. \link{dforce}
!latex       \ei
!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine coords( lvol, lss, Lcurvature, Ntz, mn )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half, one, two, pi2
  
  use numerical, only : vsmall, small
  
  use fileunits, only : ounit
  
  use inputlist, only : Wcoords, Igeometry, Ntor
  
  use cputiming, only : Tcoords
  
  use allglobal, only : myid, cpus, pi2nfp, &
                        Mvol, im, in, halfmm, &
                        iRbc, iZbs, iRbs, iZbc, &
                        NOTstellsym, Lcoordinatesingularity, &
                        Nt, Nz, isr, trigm, trign, trigwk, &
                        Rij, Zij, &
                        cosi, sini, &
                        sg, guvij, &
                        dBdX
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in) :: lvol, Lcurvature, Ntz, mn
  REAL   , intent(in) :: lss
  
  INTEGER             :: ii, jj, kk, irz, innout, issym, signlss, mi, ni, imn
  REAL                :: Remn(1:mn,0:2), Zomn(1:mn,0:2), Romn(1:mn,0:2), Zemn(1:mn,0:2), alss, blss, sbar, sbarhim(1:mn), fj(1:mn,0:2)
  
  REAL                :: Dij(1:Ntz,0:3), dguvij(1:Ntz,1:3,1:3)
  
  BEGIN(coords)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  FATAL( coords, lvol.lt.1 .or. lvol.gt.Mvol, invalid volume label )
  FATAL( coords, abs(lss).gt.one, invalid radial coordinate )
  FATAL( coords, Lcurvature.lt.0 .or. Lcurvature.gt.4, invalid input value for Lcurvature )
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
   case( 2   )  ; fj(     1:Ntor+1,0) = sbar**half              ! these are the mj.eq.0 harmonics; 11 Aug 14;
   case( 3   )  ; fj(     1:Ntor+1,0) = sbar
   case default ; FATAL( coords, .true., invalid Igeometry for Lcoordinatesingularity=T )
   end select
   ;            ; fj(Ntor+2:mn    ,0) = sbar**halfmm(Ntor+2:mn) ! these are the me.ne.0 harmonics; 11 Aug 14;
   
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
               Nt, Nz, Rij(1:Ntz,0,0), Zij(1:Ntz,0,0), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz) ) ! maps to real space;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( Lcurvature.eq.0 ) goto 9999 ! only the coordinates are required;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Lcoordinatesingularity ) then
   
#ifdef DEBUG
   FATAL( coords, sbar.lt.small, small denominator )
#endif
   
   select case( Igeometry )
   case( 2   )  ; fj(     1:Ntor+1,1) = half * half              * fj(     1:Ntor+1,0) / sbar ! these are the mj.eq.0 harmonics; 11 Aug 14;
   case( 3   )  ; fj(     1:Ntor+1,1) = half                                                  ! these are the mj.eq.0 harmonics; 11 Aug 14;
   case default ; FATAL( coords, .true., invalid Igeometry for Lcoordinatesingularity=T and Lcurvature.ne.0 )
   end select
   ;            ; fj(Ntor+2:mn    ,1) = half * halfmm(Ntor+2:mn) * fj(Ntor+2:mn    ,0) / sbar ! these are the me.ne.0 harmonics; 11 Aug 14;

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

  
  call invfft( mn, im(1:mn), in(1:mn),           Remn(1:mn,1),           Romn(1:mn,1),           Zemn(1:mn,1),           Zomn(1:mn,1), &
               Nt, Nz, Rij(1:Ntz,1,0), Zij(1:Ntz,1,0), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz) ) ! maps to real space;

  call invfft( mn, im(1:mn), in(1:mn),  im(1:mn)*Romn(1:mn,0), -im(1:mn)*Remn(1:mn,0),  im(1:mn)*Zomn(1:mn,0), -im(1:mn)*Zemn(1:mn,0), &
               Nt, Nz, Rij(1:Ntz,2,0), Zij(1:Ntz,2,0), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz) ) ! maps to real space;

  call invfft( mn, im(1:mn), in(1:mn), -in(1:mn)*Romn(1:mn,0),  in(1:mn)*Remn(1:mn,0), -in(1:mn)*Zomn(1:mn,0),  in(1:mn)*Zemn(1:mn,0), &
               Nt, Nz, Rij(1:Ntz,3,0), Zij(1:Ntz,3,0), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz) ) ! maps to real space;


  do ii = 1, 3 ; Rij(1:Ntz,0,ii) = Rij(1:Ntz,ii,0) ! just to complete workspace arrays; 22 Apr 13;
   ;           ; Zij(1:Ntz,0,ii) = Zij(1:Ntz,ii,0) ! just to complete workspace arrays; 22 Apr 13;
  enddo


  guvij(1:Ntz, 0, 0, 0) = one ! this is (only) required for the helicity integral; 22 Apr 13; REDUNDANT; see metrix; SRH; 01 Aug 17;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  select case( Igeometry )
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  case( 1 ) ! Igeometry=1; Cartesian;
   
   sg(1:Ntz,0) = Rij(1:Ntz,1,0)
   
   do ii = 1, 3
    do jj = ii, 3 ; guvij(1:Ntz,ii,jj,0) = Rij(1:Ntz,ii,0) * Rij(1:Ntz,jj,0)
    enddo
   enddo

   guvij(1:Ntz, 2, 2,0) = guvij(1:Ntz, 2, 2,0) + one
   guvij(1:Ntz, 3, 3,0) = guvij(1:Ntz, 3, 3,0) + one

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
    case( 2 )    ; fj(     1:Ntor+1,2) = half * ( half              - one ) * fj(     1:Ntor+1,1) / sbar ! these are the mj.eq.0 harmonics; 11 Aug 14;
    case( 3 )    ; fj(     1:Ntor+1,2) = zero                                                            ! these are the mj.eq.0 harmonics; 11 Aug 14;
    case default ; 
     ;           ; FATAL( coords, .true., invalid Igeometry for Lcoordinatesingularity=T and Lcurvature=2 )
    end select   ;
    ;            ; fj(Ntor+2:mn    ,2) = half * ( halfmm(Ntor+2:mn) - one ) * fj(Ntor+2:mn    ,1) / sbar ! these are the me.ne.0 harmonics; 11 Aug 14;
    
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
                Nt, Nz, Rij(1:Ntz,1,1), Zij(1:Ntz,1,1), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz) ) ! maps to real space;

   call invfft( mn, im(1:mn), in(1:mn),&
+         im(1:mn)*Romn(1:mn,1),         -im(1:mn)*Remn(1:mn,1),          im(1:mn)*Zomn(1:mn,1),         -im(1:mn)*Zemn(1:mn,1), &
Nt, Nz, Rij(1:Ntz,1,2), Zij(1:Ntz,1,2), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz) ) ! maps to real space;

   call invfft( mn, im(1:mn), in(1:mn),&
-         in(1:mn)*Romn(1:mn,1),          in(1:mn)*Remn(1:mn,1),         -in(1:mn)*Zomn(1:mn,1),          in(1:mn)*Zemn(1:mn,1), &
Nt, Nz, Rij(1:Ntz,1,3), Zij(1:Ntz,1,3), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz) ) ! maps to real space;
  
   call invfft( mn, im(1:mn), in(1:mn),&
-im(1:mn)*im(1:mn)*Remn(1:mn,0),-im(1:mn)*im(1:mn)*Romn(1:mn,0),-im(1:mn)*im(1:mn)*Zemn(1:mn,0),-im(1:mn)*im(1:mn)*Zomn(1:mn,0), &
Nt, Nz, Rij(1:Ntz,2,2), Zij(1:Ntz,2,2), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz) ) ! maps to real space;

   call invfft( mn, im(1:mn), in(1:mn),&
+im(1:mn)*in(1:mn)*Remn(1:mn,0), im(1:mn)*in(1:mn)*Romn(1:mn,0), im(1:mn)*in(1:mn)*Zemn(1:mn,0), im(1:mn)*in(1:mn)*Zomn(1:mn,0), &
Nt, Nz, Rij(1:Ntz,2,3), Zij(1:Ntz,2,3), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz) ) ! maps to real space;

   call invfft( mn, im(1:mn), in(1:mn),&
-in(1:mn)*in(1:mn)*Remn(1:mn,0),-in(1:mn)*in(1:mn)*Romn(1:mn,0),-in(1:mn)*in(1:mn)*Zemn(1:mn,0),-in(1:mn)*in(1:mn)*Zomn(1:mn,0), &
Nt, Nz, Rij(1:Ntz,3,3), Zij(1:Ntz,3,3), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz) ) ! maps to real space;

   
   do ii = 2, 3
    do jj = 1, ii-1 ; Rij(1:Ntz,ii,jj) = Rij(1:Ntz,jj,ii) ; Zij(1:Ntz,ii,jj) = Zij(1:Ntz,jj,ii)
    enddo
   enddo

   
   select case( Igeometry )
    
   case( 1 ) ! Lcurvature=2; Igeometry=1 ; Cartesian;
    
    do kk = 1, 3 ! kk labels derivative; 13 Sep 13;
     
!    sg(1:Ntz, 0) = Rij(1:Ntz,1, 0)
     sg(1:Ntz,kk) = Rij(1:Ntz,1,kk)
     
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
   
   
  case( 3,4 ) ! Lcurvature=3,4 ; get derivatives wrt R_j and Z_j; 19 Sep 13;
   
   
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
   
   
   select case( Igeometry )
    
   case( 1 ) ! Lcurvature=3,4 ; Igeometry=1 ; Cartesian; 04 Dec 14;

#ifdef DEBUG
    FATAL( coords, irz.eq.1, there is no dependence on Zbs or Zbc )
#endif
   
!                  sg(1:Ntz,0) = Rij(1:Ntz,1,0)
                   sg(1:Ntz,1) = Dij(1:Ntz,1  ) ! 20 Jun 14;
!   if( irz.eq.1 ) sg(1:Ntz,1) = 
    
    do ii = 1, 3 ! careful: ii was used with a different definition above; 13 Sep 13;
     do jj = ii, 3
                     dguvij(1:Ntz,ii,jj) = Dij(1:Ntz,ii) * Rij(1:Ntz,jj,0) + Rij(1:Ntz,ii,0) * Dij(1:Ntz,jj)
!     if( irz.eq.1 ) dguvij(1:Ntz,ii,jj) = 
     enddo
    enddo

   case( 2 ) ! Lcurvature=3,4 ; Igeometry=2 ; cylindrical;

#ifdef DEBUG
    FATAL( coords, irz.eq.1, there is no dependence on Zbs or Zbc )
#endif

!                  sg(1:Ntz,0) = Rij(1:Ntz,1,0) * Rij(1:Ntz,0,0)
    if( irz.eq.0 ) sg(1:Ntz,1) = Dij(1:Ntz,1  ) * Rij(1:Ntz,0,0) &
                               + Rij(1:Ntz,1,0) * Dij(1:Ntz,0  )

    do ii = 1, 3 ! careful: ii was used with a different definition above; 13 Sep 13;
     do jj = ii, 3
      if( irz.eq.0 ) dguvij(1:Ntz,ii,jj) = Dij(1:Ntz,ii) * Rij(1:Ntz,jj,0) + Rij(1:Ntz,ii,0) * Dij(1:Ntz,jj)
      if( irz.eq.1 ) dguvij(1:Ntz,ii,jj) = Dij(1:Ntz,ii) * Zij(1:Ntz,jj,0) + Zij(1:Ntz,ii,0) * Dij(1:Ntz,jj)
     enddo
    enddo
    
    dguvij(1:Ntz,2,2) = dguvij(1:Ntz,2,2) + two * Dij(1:Ntz,0) * Rij(1:Ntz,0,0)
   !dguvij(1:Ntz,3,3) = additional term is unity;

   case( 3 ) ! Lcurvature=3,4 ; Igeometry=3 ; toroidal; 04 Dec 14;

!                  sg(1:Ntz,0) = Rij(1:Ntz,0,0) * ( Zij(1:Ntz,1,0)*Rij(1:Ntz,2,0) - Rij(1:Ntz,1,0)*Zij(1:Ntz,2,0) )
    if( irz.eq.0 ) sg(1:Ntz,1) = Dij(1:Ntz,0  ) * ( Zij(1:Ntz,1,0)*Rij(1:Ntz,2,0) - Rij(1:Ntz,1,0)*Zij(1:Ntz,2,0) ) & 
                               + Rij(1:Ntz,0,0) * ( Zij(1:Ntz,1,0)*Dij(1:Ntz,2  ) - Dij(1:Ntz,1  )*Zij(1:Ntz,2,0) ) 
    if( irz.eq.1 ) sg(1:Ntz,1) = Rij(1:Ntz,0,0) * ( Dij(1:Ntz,1  )*Rij(1:Ntz,2,0) - Rij(1:Ntz,1,0)*Dij(1:Ntz,2  ) )

    sg(1:Ntz,1) = sg(1:Ntz,1)
    
    do ii = 1, 3 ! careful: ii was used with a different definition above; 13 Sep 13;
     do jj = ii, 3
      if( irz.eq.0 ) dguvij(1:Ntz,ii,jj) = Dij(1:Ntz,ii) * Rij(1:Ntz,jj,0) + Rij(1:Ntz,ii,0) * Dij(1:Ntz,jj)
      if( irz.eq.1 ) dguvij(1:Ntz,ii,jj) = Dij(1:Ntz,ii) * Zij(1:Ntz,jj,0) + Zij(1:Ntz,ii,0) * Dij(1:Ntz,jj)
     enddo
    enddo
    
    if( irz.eq.0 ) dguvij(1:Ntz,3,3) = dguvij(1:Ntz,3,3) + two * Dij(1:Ntz,0) * Rij(1:Ntz,0,0)

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
   
  end select ! matches select case( Lcurvature ) ; 10 Mar 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(coords)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine coords

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
