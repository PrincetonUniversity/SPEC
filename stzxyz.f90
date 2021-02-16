!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (diagnostic) ! Calculates coordinates, ${\bf x}(s,\theta,\zeta) \equiv R \, {\bf e}_R + Z \, {\bf e}_Z$, and metrics, at given $(s,\theta,\zeta)$.

!latex \briefly{briefly}

!latex \calledby{\link{pp00ab}}
!      \calls{\link{}}

!latex \tableofcontents

!latex \begin{enumerate}
!latex \item This routine is a ``copy'' of \link{co01aa},
!latex       which calculates the coordinate information on a regular, discrete grid in $\t$ and $\z$ at given $\s$;
!latex       whereas \link{stzxyz} calculates the coordinate information at a single point $(\s,\t,\z)$.
!latex \item Please see \link{co01aa} for documentation.
!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine stzxyz( lvol , stz , RpZ )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, one, half
  use numerical, only : vsmall
  use fileunits, only : ounit
  use inputlist, only : Wstzxyz, Igeometry, Nvol, Ntor
  use cputiming, only : Tstzxyz
  use allglobal, only : myid, cpus, mn, im, in, halfmm, MPI_COMM_SPEC, &
                        iRbc, iZbs, iRbs, iZbc, &
                        Lcoordinatesingularity, &
                        NOTstellsym, &
                        Mvol
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS   
  
  INTEGER, intent(in)  :: lvol
  
  REAL,    intent(in)  :: stz(1:3)
  REAL,    intent(out) :: RpZ(1:3)
  
  INTEGER              :: ii, mi, ni
  REAL                 :: Remn, Zomn, Romn, Zemn, RR, phi, ZZ, arg, carg, sarg, lss, alss, blss, sbar, sbarhim, fj
  
  BEGIN(stzxyz)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  FATAL(stzxyz, lvol.lt.1 .or. lvol.gt.Mvol, invalid interface label )
  FATAL(stzxyz, abs(stz(1)).gt.one, invalid radial coordinate )
#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RpZ(1:3) = zero ; RR = zero ; phi = stz(3) ; ZZ = zero ! initialize intent(out), summations ; 17 Dec 15;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  lss = stz(1) ; alss = half - lss * half ; blss = half + lss * half ! shorthand; 17 Dec 15;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  do ii = 1, mn ; mi = im(ii) ; ni = in(ii) ; arg = mi * stz(2) - ni * phi ; carg = cos(arg) ; sarg = sin(arg) ! Fourier sum; 17 Dec 15;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   Remn = zero ; Zomn = zero
!  if( NOTstellsym ) then
   Romn = zero ; Zemn = zero
!  endif
   
   if( Lcoordinatesingularity ) then
    
    sbar = ( lss + one ) * half
    if( mi.eq.0 ) then
     if( Igeometry.eq.2 ) then ; fj = sbar
     else                      ; fj = sbar**2
     endif
    else                       ; fj = sbar**im(ii)
    endif
    
    Remn = iRbc(ii,0) + ( iRbc(ii,1) - iRbc(ii,0) ) * fj
    if( NOTstellsym ) then
    Romn = iRbs(ii,0) + ( iRbs(ii,1) - iRbs(ii,0) ) * fj
    endif
    if( Igeometry.eq.3 ) then ! recall that for cylindrical geometry there is no need for Z; 20 Apr 13;
    Zomn = iZbs(ii,0) + ( iZbs(ii,1) - iZbs(ii,0) ) * fj
    if( NOTstellsym ) then
    Zemn = iZbc(ii,0) + ( iZbc(ii,1) - iZbc(ii,0) ) * fj
    endif
    endif

   else ! matches if( Lcoordinatesingularity) ; 20 Apr 13;
    
    Remn =   alss * iRbc(ii,lvol-1) + blss * iRbc(ii,lvol)
    if( NOTstellsym ) then
    Romn =   alss * iRbs(ii,lvol-1) + blss * iRbs(ii,lvol)
    else
    Romn = zero
    endif ! end of if( NOTstellsym ) ; 22 Apr 13;
    if( Igeometry.eq.3 ) then
    Zomn =   alss * iZbs(ii,lvol-1) + blss * iZbs(ii,lvol)
    if( NOTstellsym ) then
    Zemn =   alss * iZbc(ii,lvol-1) + blss * iZbc(ii,lvol)
    else
    Zemn = zero
    endif ! end of if( NOTstellsym ) ; 22 Apr 13; 
    else    
    Zomn = zero
    Zemn = zero
    endif ! end of if( Igeometry.eq.3 ) ; 22 Apr 13;
    
   endif ! end of if( Lcoordinatesingularity ) ; 20 Feb 13;
   
   RR = RR + Remn * carg + Romn * sarg
   if( Igeometry.eq.3 ) then
   ZZ = ZZ + Zemn * carg + Zomn * sarg
   endif
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
  enddo ! end of do ii = 1, mn ; Fourier summation loop;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RpZ(1:3) = (/ RR, phi, ZZ /) ! return coordinate functions;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(stzxyz)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine stzxyz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
