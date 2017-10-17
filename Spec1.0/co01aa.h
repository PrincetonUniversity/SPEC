!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Calculates coordinates, metric elements and curvatures using fast Fourier transform on regular grid.

!latex \end{enumerate} \subsection{Coordinates} \begin{enumerate}

!latex \item We shall work in coordinates, $(\s,\t,\z)$, which are be defined inversely via a transformation to Cartesian coordinates, $(x,y,z)$.
!latex       The toroidal angle, $\z$, is identical to the cylindrical angle, $\z\equiv\phi$.
!latex       The radial coordinate, $s$, is {\em not} a global variable: it only needs to be defined in each volume.
!latex       The poloidal angle, $\t$, remains arbitrary.

!latex \item The geometry of the interfaces, ${\bf x}_v(\t,\z)$, is given by
!latex       \begin{itemize}
!latex       \item \verb+Igeometry=1+ : Cartesian
!latex       \be {\bf x} & \equiv &  \t \; {\bf \hat i} + \z \; {\bf \hat j}+ R \; {\bf \hat k}
!latex       \ee
!latex       \item \verb+Igeometry=2+ : Cylindrical
!latex       \be {\bf x} & = & R \; \cos\t \; {\bf \hat i} + R \sin\t \; {\bf \hat j} + \z \; {\bf \hat k}
!latex       \ee
!latex       \item \verb+Igeometry=3+ : Toroidal
!latex       \be {\bf x} & \equiv & R \; {\bf \hat r} + Z \; {\bf \hat k} 
!latex       \ee
!latex       where ${\bf \hat r}\equiv \cos \phi \; {\bf \hat i} + \sin \phi \; {\bf \hat j}$ and
!latex             ${\bf \hat \phi} \equiv - \sin \phi \; {\bf \hat i} + \cos \phi \; {\bf \hat j}$.
!latex       \item \verb+Igeometry=4+ : Extended cylindrical;
!latex       \be {\bf x} & \equiv & R \; {\bf \hat r} + Z \; {\bf \hat k} 
!latex       \ee
!latex       where ${\bf \hat r}\equiv \cos \phi \; {\bf \hat i} + \sin \phi \; {\bf \hat j}$ and
!latex             ${\bf \hat \phi} \equiv - \sin \phi \; {\bf \hat i} + \cos \phi \; {\bf \hat j}$.
!latex       \end{itemize}

!latex \end{enumerate} \subsection{interpolation and coordinate regularization} \begin{enumerate}

!latex \item The $v$-th volume is bounded by the ${\bf x}_{v-1}$ and ${\bf x}_{v}$.

!latex \item In each volume, the coordinates are constructed by linear interpolation
!latex       \be R(\s,\t,\z) &\equiv& \left[ \; (1-s) \; R_{v-1} + (1+s) \; R_{v}\; \right] \; / \; 2.\\
!latex           Z(\s,\t,\z) &\equiv& \left[ \; (1-s) \; Z_{v-1} + (1+s) \; Z_{v}\; \right] \; / \; 2.
!latex       \ee

!latex \item Note that $\s \in [-1,1]$ in each volume.

!latex \item For the innermost volume, for the cylindrical and toriodal case, there is a singularity in the coordinates.
!latex       A regular coordinate system can be enforced by including regularization factors:
!latex       \be R(\s,\t,\z) & \equiv & \sum_j R_{1,e,j} \bar s^{m_j/2} \cos \a_j + \sum_j R_{1,o,j} \bar s^{m_j/2} \sin \a_j \\
!latex           Z(\s,\t,\z) & \equiv & \sum_j Z_{1,e,j} \bar s^{m_j/2} \cos \a_j + \sum_j Z_{1,o,j} \bar s^{m_j/2} \sin \a_j,
!latex       \ee
!latex       where $\bar s \equiv (s+1)/2$.

!latex \item The location of the coordinate axis is set in the ``unpacking'' routine, \verb+gf00aa+.

!latex \end{enumerate} \subsection{logic overview}

!latex The input parameter \verb+Lcurvature+ controls the operation:

!latex \begin{enumerate}

!latex \item if \verb+Lcurvature=0+: The coordinates are returned.
!latex \item if \verb+Lcurvature=1+: The coordinates and the metric elements are returned. This is required for the computation of volume integrals.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine co01aa( lvol, lss, Lcurvature, Ntz, mn )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half, one, two, pi2
  
  use numerical, only : vsmall, small
  
  use fileunits, only : ounit
  
  use inputlist, only : Wco01aa, Igeometry, pknot, Ntor
  
  use cputiming, only : Tco01aa
  
  use allglobal, only : myid, cpus, pi2nfp, &
                        Mvol, im, in, halfmm, &
                        iRbc, iZbs, iRbs, iZbc, &
                        NOTstellsym, Lcoordinatesingularity, &
                        Nt, Nz, isr, trigm, trign, trigwk, &
                        Rij, Zij, &
                        cosi, sini, &
                        sg, guvij, &
                        DifferentiateGeometry
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in) :: lvol, Lcurvature, Ntz, mn
  REAL   , intent(in) :: lss
  
  INTEGER             :: ii, jj, kk, irz, innout, issym, signlss, mi, ni, imn
  REAL                :: Remn(1:mn,0:2), Zomn(1:mn,0:2), Romn(1:mn,0:2), Zemn(1:mn,0:2), alss, blss, sbar, sbarhim(1:mn), fj(1:mn,0:2)
  
  REAL                :: Dij(1:Ntz,0:3), dguvij(1:Ntz,1:3,1:3)
  
  BEGIN(co01aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  FATALMESS(co01aa, lvol.lt.1 .or. lvol.gt.Mvol, invalid volume label)
  FATALMESS(co01aa, abs(lss).gt.one, invalid radial coordinate)
  FATALMESS(co01aa, Lcurvature.lt.0 .or. Lcurvature.gt.4, invalid input value for Lcurvature)
#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  Rij(1:Ntz,0:3,0:3) = zero ; sg(1:Ntz,0:3) = zero ; guvij(1:Ntz,1:3,1:3,0:3) = zero ! provide trivial default for output; 16 Jan 13;
  Zij(1:Ntz,0:3,0:3) = zero                                                          ! provide trivial default for output; 16 Jan 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  Remn(1:mn,0:2) = zero ! interpolated coordinate harmonics; 6 Feb 13;
  Zomn(1:mn,0:2) = zero
! if( NOTstellsym ) then
  Romn(1:mn,0:2) = zero
  Zemn(1:mn,0:2) = zero
! endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  
  select case( Igeometry )

  case( 1:2 )

   do imn = 1, mn
   !if( im(imn).ne.0 ) then ! 14 Jan 15;
     ; FATALMESS(co01aa, abs(iRbc(imn,0)).gt.vsmall, coordinate axis cannot have angle dependence)
     if( NOTstellsym ) then
      ;FATALMESS(co01aa, abs(iRbs(imn,0)).gt.vsmall, coordinate axis cannot have angle dependence)
     endif
     if( Igeometry.ge.3 ) then
      ;FATALMESS(co01aa, abs(iZbs(imn,0)).gt.vsmall, coordinate axis cannot have angle dependence)
      if( NOTstellsym ) then
       FATALMESS(co01aa, abs(iZbc(imn,0)).gt.vsmall, coordinate axis cannot have angle dependence)
      endif
     endif
   !endif ! end of if( im(imn).ne.0 ) ;
   enddo

  case( 3:4 )

   do imn = 1, mn
   !if(    imn).eq.0 ) then
    if( im(imn).ne.0 ) then ! 14 Jan 15;
     ; FATALMESS(co01aa, abs(iRbc(imn,0)).gt.vsmall, coordinate axis cannot have angle dependence)
     if( NOTstellsym ) then
      ;FATALMESS(co01aa, abs(iRbs(imn,0)).gt.vsmall, coordinate axis cannot have angle dependence)
     endif
     if( Igeometry.ge.3 ) then
      ;FATALMESS(co01aa, abs(iZbs(imn,0)).gt.vsmall, coordinate axis cannot have angle dependence)
      if( NOTstellsym ) then
       FATALMESS(co01aa, abs(iZbc(imn,0)).gt.vsmall, coordinate axis cannot have angle dependence)
      endif
     endif
    endif ! end of if( imn.eq.0 ) ; 04 Dec 14;
   enddo

  case default
    
   FATALMESS(co01aa, .true., invalid Igeometry for Lcoordinatesingularity=T)

  end select
  
#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Lcoordinatesingularity ) then
   
   sbar = ( lss + one ) * half

#ifdef DEBUG
   FATALMESS(co01aa,sbar.lt.zero .or. sbar.gt.one, invalid sbar)
#endif
   
   select case( Igeometry )
   case( 2   )  ; fj(     1:Ntor+1,0) = sbar**half              ! these are the mj.eq.0 harmonics; 11 Aug 14;
   case( 3:4 )  ; fj(     1:Ntor+1,0) = sbar
   case default ; FATALMESS(co01aa, .true., invalid Igeometry for Lcoordinatesingularity=T)
   end select
   ;            ; fj(Ntor+2:mn    ,0) = sbar**halfmm(Ntor+2:mn) ! these are the me.ne.0 harmonics; 11 Aug 14;
   
   Remn(1:mn,0) = iRbc(1:mn,0) + ( iRbc(1:mn,1) - iRbc(1:mn,0) ) * fj(1:mn,0)
   if( NOTstellsym ) then
   Romn(1:mn,0) = iRbs(1:mn,0) + ( iRbs(1:mn,1) - iRbs(1:mn,0) ) * fj(1:mn,0)
   endif
   if( Igeometry.ge.3 ) then
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
   if( Igeometry.ge.3 ) then
   Zomn(1:mn,0) = alss * iZbs(1:mn,lvol-1) + blss * iZbs(1:mn,lvol)
   if( NOTstellsym ) then
   Zemn(1:mn,0) = alss * iZbc(1:mn,lvol-1) + blss * iZbc(1:mn,lvol)
   endif
   endif ! end of if( Igeometry.ge.3 ) ;  1 Feb 13;
    
  endif ! end of if( Lcoordinatesingularity );  1 Feb 13;


  call invfft( mn, im(1:mn), in(1:mn), Remn(1:mn,0), Romn(1:mn,0), Zemn(1:mn,0), Zomn(1:mn,0), &
               Nt, Nz, Rij(1:Ntz,0,0), Zij(1:Ntz,0,0), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz) ) ! maps to real space;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( Lcurvature.eq.0 ) goto 9999 ! only the coordinates are required;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Lcoordinatesingularity ) then
   
#ifdef DEBUG
   FATALMESS(co01aa, sbar.lt.small, small denominator)
#endif
   
   select case( Igeometry )
   case( 2   )  ; fj(     1:Ntor+1,1) = half * half              * fj(     1:Ntor+1,0) / sbar ! these are the mj.eq.0 harmonics; 11 Aug 14;
   case( 3:4 )  ; fj(     1:Ntor+1,1) = half                                                  ! these are the mj.eq.0 harmonics; 11 Aug 14;
   case default ; FATALMESS(co01aa, .true., invalid Igeometry for Lcoordinatesingularity=T and Lcurvature.ne.0)
   end select
   ;            ; fj(Ntor+2:mn    ,1) = half * halfmm(Ntor+2:mn) * fj(Ntor+2:mn    ,0) / sbar ! these are the me.ne.0 harmonics; 11 Aug 14;

   Remn(1:mn,1) =                       ( iRbc(1:mn,1) - iRbc(1:mn,0) ) * fj(1:mn,1)
   if( NOTstellsym ) then
   Romn(1:mn,1) =                       ( iRbs(1:mn,1) - iRbs(1:mn,0) ) * fj(1:mn,1)
   endif
   if( Igeometry.ge.3 ) then
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
   if( Igeometry.ge.3 ) then
   Zomn(1:mn,1) = (      - iZbs(1:mn,lvol-1) +        iZbs(1:mn,lvol) ) * half
   if( NOTstellsym ) then
   Zemn(1:mn,1) = (      - iZbc(1:mn,lvol-1) +        iZbc(1:mn,lvol) ) * half
   endif     
   endif ! end of if( Igeometry.ge.3 ) ;  1 Feb 13;
    
  endif ! end of if( Lcoordinatesingularity );  1 Feb 13;

  
  call invfft( mn, im(1:mn), in(1:mn),           Remn(1:mn,1),           Romn(1:mn,1),           Zemn(1:mn,1),           Zomn(1:mn,1), &
               Nt, Nz, Rij(1:Ntz,1,0), Zij(1:Ntz,1,0), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz) ) ! maps to real space;

  call invfft( mn, im(1:mn), in(1:mn),  im(1:mn)*Romn(1:mn,0), -im(1:mn)*Remn(1:mn,0),  im(1:mn)*Zomn(1:mn,0), -im(1:mn)*Zemn(1:mn,0), &
               Nt, Nz, Rij(1:Ntz,2,0), Zij(1:Ntz,2,0), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz) ) ! maps to real space;

  call invfft( mn, im(1:mn), in(1:mn), -in(1:mn)*Romn(1:mn,0),  in(1:mn)*Remn(1:mn,0), -in(1:mn)*Zomn(1:mn,0),  in(1:mn)*Zemn(1:mn,0), &
               Nt, Nz, Rij(1:Ntz,3,0), Zij(1:Ntz,3,0), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz) ) ! maps to real space;


  do ii = 1, 3 ; Rij(1:Ntz,0,ii) = Rij(1:Ntz,ii,0) ! just to complete workspace arrays; 22 Apr 13;
   ;           ; Zij(1:Ntz,0,ii) = Zij(1:Ntz,ii,0) ! just to complete workspace arrays; 22 Apr 13;
  enddo


  guvij(1:Ntz, 0, 0, 0) = one ! this is (only) required for the helicity integral; 22 Apr 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  select case( Igeometry )
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  case( 1 ) ! Cartesian;
   
   sg(1:Ntz,0) = Rij(1:Ntz,1,0)
   
   do ii = 1, 3
    do jj = ii, 3 ; guvij(1:Ntz,ii,jj,0) = Rij(1:Ntz,ii,0) * Rij(1:Ntz,jj,0)
    enddo
   enddo

   guvij(1:Ntz, 2, 2,0) = guvij(1:Ntz, 2, 2,0) + one
   guvij(1:Ntz, 3, 3,0) = guvij(1:Ntz, 3, 3,0) + one

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  case( 2 ) ! cylindrical - standard;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \end{enumerate} \subsection{cylindrical metrics} \begin{enumerate}
!latex \item The cylindrical metrics and Jacobian are
!latex       \be \sqrt g  &=& R_\s R          \\
!latex           g_{\s\s} &=& R_\s R_\s       \\
!latex           g_{\s\t} &=& R_\s R_\t       \\
!latex           g_{\s\z} &=& R_\s R_\z       \\
!latex           g_{\t\t} &=& R_\t R_\t + R^2 \\
!latex           g_{\t\z} &=& R_\t R_\z       \\
!latex           g_{\z\z} &=& R_\z R_\z + 1
!latex       \ee

   sg(1:Ntz,0) = Rij(1:Ntz, 1,0) * Rij(1:Ntz, 0,0)
   
   do ii = 1, 3
    do jj = ii, 3 ; guvij(1:Ntz,ii,jj,0) = Rij(1:Ntz,ii,0) * Rij(1:Ntz,jj,0)
    enddo
   enddo

   guvij(1:Ntz, 2, 2,0) = guvij(1:Ntz, 2, 2,0) + Rij(1:Ntz, 0, 0) * Rij(1:Ntz, 0, 0)
   guvij(1:Ntz, 3, 3,0) = guvij(1:Ntz, 3, 3,0) + one

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  case( 3 ) ! toroidal;
   
   sg(1:Ntz,0) = pknot * Rij(1:Ntz,0,0) * ( Zij(1:Ntz,1,0)*Rij(1:Ntz,2,0) - Rij(1:Ntz,1,0)*Zij(1:Ntz,2,0) )
   
   do ii = 1, 3
    do jj = ii, 3 ; guvij(1:Ntz,ii,jj,0) = Rij(1:Ntz,ii,0) * Rij(1:Ntz,jj,0) + Zij(1:Ntz,ii,0) * Zij(1:Ntz,jj,0)
    enddo
   enddo

   guvij(1:Ntz, 3, 3,0) = guvij(1:Ntz, 3, 3,0) + ( pknot * Rij(1:Ntz,0,0) )**2

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  case( 4 ) ! cylindrical - extended;
   
   sg(1:Ntz,0) = Rij(1:Ntz,1,0)*Zij(1:Ntz,2,0) - Zij(1:Ntz,1,0)*Rij(1:Ntz,2,0)
   
   do ii = 1, 3
    do jj = ii, 3 ; guvij(1:Ntz,ii,jj,0) = Rij(1:Ntz,ii,0) * Rij(1:Ntz,jj,0) + Zij(1:Ntz,ii,0) * Zij(1:Ntz,jj,0)
    enddo
   enddo

   guvij(1:Ntz, 3, 3,0) = guvij(1:Ntz, 3, 3,0) + one

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  case default

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   FATALMESS(co01aa, .true., selected Igeometry not supported)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
      
  end select

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  do ii = 2, 3
   do jj = 1, ii-1 ; guvij(1:Ntz,ii,jj,0) = guvij(1:Ntz,jj,ii,0) ! complete metric array; 20 Apr 13;
   enddo
  enddo
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Lcurvature.le.1 ) goto 9999
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  select case( Lcurvature )
   
  case( 2 ) ! get second derivatives of position wrt \s, \t & \z; 19 Sep 13;
   
   if( Lcoordinatesingularity ) then
    
#ifdef DEBUG
    FATALMESS(co01aa, sbar.lt.small, small denominator)
#endif

    select case( Igeometry )
    case( 2   )  ; fj(     1:Ntor+1,2) = half * ( half              - one ) * fj(     1:Ntor+1,1) / sbar ! these are the mj.eq.0 harmonics; 11 Aug 14;
    case( 3:4 )  ; fj(     1:Ntor+1,2) = zero                                                            ! these are the mj.eq.0 harmonics; 11 Aug 14;
    case default ; FATALMESS(co01aa, .true., invalid Igeometry for Lcoordinatesingularity=T and Lcurvature=2)
    end select
    ;            ; fj(Ntor+2:mn    ,2) = half * ( halfmm(Ntor+2:mn) - one ) * fj(Ntor+2:mn    ,1) / sbar ! these are the me.ne.0 harmonics; 11 Aug 14;

    Remn(1:mn,2) =                       ( iRbc(1:mn,1) - iRbc(1:mn,0) ) * fj(1:mn,2)
    if( NOTstellsym ) then
    Romn(1:mn,2) =                       ( iRbs(1:mn,1) - iRbs(1:mn,0) ) * fj(1:mn,2)
    endif
    if( Igeometry.ge.3 ) then
    Zomn(1:mn,2) =                       ( iZbs(1:mn,1) - iZbs(1:mn,0) ) * fj(1:mn,2)
    if( NOTstellsym ) then
    Zemn(1:mn,2) =                       ( iZbc(1:mn,1) - iZbc(1:mn,0) ) * fj(1:mn,2)
    endif
    endif
    
   else
    
    Remn(1:mn,2) = zero
    if( NOTstellsym ) then
    Romn(1:mn,2) = zero
    endif
    if( Igeometry.ge.3 ) then
    Zomn(1:mn,2) = zero
    if( NOTstellsym ) then
    Zemn(1:mn,2) = zero
    endif
    endif ! end of if( Igeometry.ge.3 ) ;  1 Feb 13;
    
   endif ! end of if( Lcoordinatesingularity );  1 Feb 13;
   

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
    
   case( 1 ) ! Cartesian;
    
    do kk = 1, 3 ! kk labels derivative; 13 Sep 13;
     
!    sg(1:Ntz, 0) = Rij(1:Ntz,1, 0)
     sg(1:Ntz,kk) = Rij(1:Ntz,1,kk)
     
     do ii = 1, 3
      do jj = ii, 3 ; guvij(1:Ntz,ii,jj,kk) = Rij(1:Ntz,ii,kk) * Rij(1:Ntz,jj, 0) + Rij(1:Ntz,ii, 0) * Rij(1:Ntz,jj,kk)
      enddo
     enddo
     
    enddo !  6 Feb 13;
    
   case( 2 ) ! cylindrical - standard;
    
    do kk = 1, 3 ! kk labels derivative; 13 Sep 13;
     
!    sg(1:Ntz, 0) = Rij(1:Ntz, 1, 0) * Rij(1:Ntz, 0,0)
     sg(1:Ntz,kk) = Rij(1:Ntz, 1,kk) * Rij(1:Ntz, 0,0) + Rij(1:Ntz, 1,0) * Rij(1:Ntz, 0,kk)
     
     do ii = 1, 3
      do jj = ii, 3 ; guvij(1:Ntz,ii,jj,kk) = Rij(1:Ntz,ii,kk) * Rij(1:Ntz,jj, 0) + Rij(1:Ntz,ii, 0) * Rij(1:Ntz,jj,kk)
      enddo
     enddo

     guvij(1:Ntz,2,2,kk) = guvij(1:Ntz,2,2,kk) + Rij(1:Ntz,0,kk) * Rij(1:Ntz,0,0) + Rij(1:Ntz,0,0) * Rij(1:Ntz,0,kk) ! 20 Jan 15;
    !guvij(1:Ntz,2,2,kk) = guvij(1:Ntz,2,2,kk) + Rij(1:Ntz,kk,0) * Rij(1:Ntz,0,0) + Rij(1:Ntz,0,0) * Rij(1:Ntz,kk,0) ! 20 Jan 15;
    !guvij(1:Ntz,3,3,kk) = additional term is unity; !  Jan 15;
     
    enddo !  6 Feb 13;
    
   case( 3 ) ! toroidal;

    do kk =  1 , 3 ! kk labels derivative; 13 Sep 13;
     
!    sg(1:Ntz, 0) = Rij(1:Ntz, 0,0) * ( Zij(1:Ntz,1, 0)*Rij(1:Ntz,2, 0) - Rij(1:Ntz,1, 0)*Zij(1:Ntz,2, 0) )
     sg(1:Ntz,kk) = Rij(1:Ntz,kk,0) * ( Zij(1:Ntz,1, 0)*Rij(1:Ntz,2, 0) - Rij(1:Ntz,1, 0)*Zij(1:Ntz,2, 0) ) & 
                  + Rij(1:Ntz, 0,0) * ( Zij(1:Ntz,1,kk)*Rij(1:Ntz,2, 0) - Rij(1:Ntz,1,kk)*Zij(1:Ntz,2, 0) ) & 
                  + Rij(1:Ntz, 0,0) * ( Zij(1:Ntz,1, 0)*Rij(1:Ntz,2,kk) - Rij(1:Ntz,1, 0)*Zij(1:Ntz,2,kk) )

     sg(1:Ntz,kk) = pknot * sg(1:Ntz,kk)
          
     do ii = 1, 3
      do jj = ii, 3
       guvij(1:Ntz,ii,jj,kk) = Rij(1:Ntz,ii,kk) * Rij(1:Ntz,jj, 0) &
                             + Rij(1:Ntz,ii, 0) * Rij(1:Ntz,jj,kk) &
                             + Zij(1:Ntz,ii,kk) * Zij(1:Ntz,jj, 0) &
                             + Zij(1:Ntz,ii, 0) * Zij(1:Ntz,jj,kk)
      enddo
     enddo
     
     guvij(1:Ntz,3,3,kk) = guvij(1:Ntz,3,3,kk) + pknot * ( Rij(1:Ntz,0,kk) * Rij(1:Ntz,0,0) + Rij(1:Ntz,0,0) * Rij(1:Ntz,0,kk) )
     
    enddo ! end of do kk;  5 Feb 13;
    
   case( 4 ) ! cylindrical - extended;

    FATALMESS(co01aa,.true., Lcurvature=2 and Igeometry=4 is under construction)

    do kk =  1 , 3 ! kk labels derivative; 13 Sep 13;
   
!    sg(1:Ntz, 0) = Rij(1:Ntz,1, 0)*Zij(1:Ntz,2, 0) - Zij(1:Ntz,1, 0)*Rij(1:Ntz,2, 0)
     sg(1:Ntz,kk) = Rij(1:Ntz,1,kk)*Zij(1:Ntz,2, 0) - Zij(1:Ntz,1,kk)*Rij(1:Ntz,2, 0) &
                  + Rij(1:Ntz,1, 0)*Zij(1:Ntz,2,kk) - Zij(1:Ntz,1, 0)*Rij(1:Ntz,2,kk)
     
     do ii = 1, 3
      do jj = ii, 3
       guvij(1:Ntz,ii,jj,kk) = Rij(1:Ntz,ii,kk) * Rij(1:Ntz,jj, 0) &
                             + Rij(1:Ntz,ii, 0) * Rij(1:Ntz,jj,kk) &
                             + Zij(1:Ntz,ii,kk) * Zij(1:Ntz,jj, 0) &
                             + Zij(1:Ntz,ii, 0) * Zij(1:Ntz,jj,kk)
      enddo
     enddo
     
    !guvij(1:Ntz,3,3,kk) = guvij(1:Ntz,3,3,kk) + zero
     
    enddo ! end of do kk;  5 Feb 13;
    
   case default
    
    FATALMESS(co01aa, .true., selected Igeometry not supported for Lcurvature.eq.2)
    
   end select
   
   do ii = 2, 3
    do jj = 1, ii-1 ; guvij(1:Ntz,ii,jj,1:3) = guvij(1:Ntz,jj,ii,1:3)
    enddo
   enddo
   
   
  case( 3,4 ) ! matches select case( Lcurvature ) ! get derivatives wrt R_j and Z_j; 19 Sep 13;
   
   
   ii = DifferentiateGeometry%ii ; innout = DifferentiateGeometry%innout ; irz = DifferentiateGeometry%irz ; issym = DifferentiateGeometry%issym ! shorthand;
   
   if( Lcoordinatesingularity ) then
    
#ifdef DEBUG
    FATALMESS(co01aa, innout.eq.0, cannot differentiate metric elements wrt coordinate singularity)
#endif
    
#ifdef DEBUG
    FATALMESS( co01aa, Igeometry.eq.1, Cartesian does not need regularization factor) ! 29 Apr 14;
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
    
   case( 1 ) ! Cartesian; 04 Dec 14;

#ifdef DEBUG
    FATALMESS(co01aa, irz.eq.1, there is no dependence on Zbs or Zbc)
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

   case( 2 ) ! cylindrical - standard; 04 Dec 14;

#ifdef DEBUG
    FATALMESS(co01aa, irz.eq.1, there is no dependence on Zbs or Zbc)
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

   case( 3 ) ! toroidal; 04 Dec 14;

!                  sg(1:Ntz,0) = Rij(1:Ntz,0,0) * ( Zij(1:Ntz,1,0)*Rij(1:Ntz,2,0) - Rij(1:Ntz,1,0)*Zij(1:Ntz,2,0) )
    if( irz.eq.0 ) sg(1:Ntz,1) = Dij(1:Ntz,0  ) * ( Zij(1:Ntz,1,0)*Rij(1:Ntz,2,0) - Rij(1:Ntz,1,0)*Zij(1:Ntz,2,0) ) & 
                               + Rij(1:Ntz,0,0) * ( Zij(1:Ntz,1,0)*Dij(1:Ntz,2  ) - Dij(1:Ntz,1  )*Zij(1:Ntz,2,0) ) 
    if( irz.eq.1 ) sg(1:Ntz,1) = Rij(1:Ntz,0,0) * ( Dij(1:Ntz,1  )*Rij(1:Ntz,2,0) - Rij(1:Ntz,1,0)*Dij(1:Ntz,2  ) )

    sg(1:Ntz,1) = pknot * sg(1:Ntz,1)
    
    do ii = 1, 3 ! careful: ii was used with a different definition above; 13 Sep 13;
     do jj = ii, 3
      if( irz.eq.0 ) dguvij(1:Ntz,ii,jj) = Dij(1:Ntz,ii) * Rij(1:Ntz,jj,0) + Rij(1:Ntz,ii,0) * Dij(1:Ntz,jj)
      if( irz.eq.1 ) dguvij(1:Ntz,ii,jj) = Dij(1:Ntz,ii) * Zij(1:Ntz,jj,0) + Zij(1:Ntz,ii,0) * Dij(1:Ntz,jj)
     enddo
    enddo
    
    if( irz.eq.0 ) dguvij(1:Ntz,3,3) = dguvij(1:Ntz,3,3) + pknot * two * Dij(1:Ntz,0) * Rij(1:Ntz,0,0)

   case( 4 ) ! cylindrical - extended;

!                  sg(1:Ntz,0) = Rij(1:Ntz,1, 0)*Zij(1:Ntz,2, 0) - Zij(1:Ntz,1, 0)*Rij(1:Ntz,2, 0)
    if( irz.eq.0 ) sg(1:Ntz,1) = Dij(1:Ntz,1   )*Zij(1:Ntz,2, 0) - Zij(1:Ntz,1, 0)*Dij(1:Ntz,2   )
    if( irz.eq.1 ) sg(1:Ntz,1) = Rij(1:Ntz,1, 0)*Dij(1:Ntz,2   ) - Dij(1:Ntz,1   )*Rij(1:Ntz,2, 0)

    do ii = 1, 3 ! careful: ii was used with a different definition above; 13 Sep 13;
     do jj = ii, 3
      if( irz.eq.0 ) dguvij(1:Ntz,ii,jj) = Dij(1:Ntz,ii) * Rij(1:Ntz,jj,0) + Rij(1:Ntz,ii,0) * Dij(1:Ntz,jj)
      if( irz.eq.1 ) dguvij(1:Ntz,ii,jj) = Dij(1:Ntz,ii) * Zij(1:Ntz,jj,0) + Zij(1:Ntz,ii,0) * Dij(1:Ntz,jj)
     enddo
    enddo
    
!   if( irz.eq.0 ) dguvij(1:Ntz,3,3) = dguvij(1:Ntz,3,3) + zero
    
   case default
    
    FATALMESS(co01aa, .true., supplied Igeometry is not yet supported for Lcurvature.eq.3 or Lcurvature.eq.4)
    
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
    
   else
    
    do ii = 1, 3
     do jj = 1, 3 ; guvij(1:Ntz,ii,jj,1) = dguvij(1:Ntz,ii,jj)                                                    ! differentiated metric elements; 7 Mar 13; 
     enddo
    enddo
    
   endif
   
  end select ! matches select case( Lcurvature ) ; 10 Mar 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(co01aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine co01aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
