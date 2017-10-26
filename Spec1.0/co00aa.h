!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Calculates coordinates, and metrics, at arbitrary position.

!latex \end{enumerate} \subsection{Coordinates} \begin{enumerate}

!latex \item This routine is essentially a copy of \verb+co01aa+, which calculates the coordinates and metrics on a discrete grid in $\t$, $\z$.

!latex \item The location of the coordinate axis is set in the ``unpacking'' routine, \verb+gf00aa+.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine co00aa( lvol , stz , RpZ , dR , dZ , jacobian, guv )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, one, half
  use numerical, only : vsmall
  use fileunits, only : ounit
  use inputlist, only : Wco00aa, Igeometry, pknot, Nvol, Ntor
  use cputiming, only : Tco00aa
  use allglobal, only : myid, cpus, mn, im, in, halfmm, &
                        iRbc, iZbs, iRbs, iZbc, &
                        Lcoordinatesingularity, NOTstellsym, Lplasmaregion, Lvacuumregion, &
                        Mvol
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS   
  
  INTEGER, intent(in)  :: lvol
  
  REAL,    intent(in)  :: stz(1:3)
  REAL,    intent(out) :: RpZ(1:3), dR(0:3), dZ(0:3), jacobian, guv(1:3,1:3) ! coordinates R,\phi,Z;
  
  INTEGER              :: imn, ii, jj, mi, ni
  REAL                 :: Remn(0:1), Zomn(0:1), Romn(0:1), Zemn(0:1), arg, carg, sarg, lss, alss, blss, sbar, sbarhim, fj
  
  BEGIN(co00aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  FATALMESS(co00aa, lvol.lt.1 .or. lvol.gt.Mvol, invalid interface label )
  FATALMESS(co00aa, abs(stz(1)).gt.one, invalid radial coordinate )
#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RpZ(1:3) = zero ; jacobian = zero ! initialize intent(out) quantities;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  dR(0:3) = zero ! initialize summation quantities to zero; ! 28 Jan 13;
  dZ(0:3) = zero ! initialize summation quantities to zero; ! 28 Jan 13; ! not required if Igeometry < 3; 01 May 13;
  
  lss = stz(1) ; alss = one - lss ; blss = one + lss ! sbar = ( lss + one ) * half ! shorthand; 16 Jan 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  do imn = 1, mn ; mi = im(imn) ; ni = in(imn) ; arg = mi * stz(2) - ni * stz(3) ; carg = cos(arg) ; sarg = sin(arg)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   Remn(0:1) = zero
   Zomn(0:1) = zero
  !if( NOTstellsym ) then
   Romn(0:1) = zero
   Zemn(0:1) = zero
  !endif

   if( Lcoordinatesingularity ) then ! WHY ARE Remn(1) etc. NOT REQUIRED; 20 Apr 13; ! the Remn(1) etc are only required in the vacuum; 01 May 13;
    
    sbar = ( lss + one ) * half
    
    if( mi.eq.0 ) then
     if( Igeometry.eq.2 ) then ; fj = sbar**half
     else                      ; fj = sbar
     endif
    else                       ; fj = sbar**halfmm(imn)
    endif
    
    Remn(0) = iRbc(imn,0) + ( iRbc(imn,1) - iRbc(imn,0) ) * fj
    if( NOTstellsym ) then
    Romn(0) = iRbs(imn,0) + ( iRbs(imn,1) - iRbs(imn,0) ) * fj
    endif
    if( Igeometry.ge.3 ) then ! recall that for cylindrical geometry there is no need for Z; 20 Apr 13;
    Zomn(0) = iZbs(imn,0) + ( iZbs(imn,1) - iZbs(imn,0) ) * fj
    if( NOTstellsym ) then
    Zemn(0) = iZbc(imn,0) + ( iZbc(imn,1) - iZbc(imn,0) ) * fj
    endif
    endif

   else ! matches if( Lcoordinatesingularity) ; 20 Apr 13;
    
    Remn(0) = ( alss * iRbc(imn,lvol-1) + blss * iRbc(imn,lvol) ) * half
    Remn(1) = ( -one * iRbc(imn,lvol-1) +  one * iRbc(imn,lvol) ) * half
    if( NOTstellsym ) then
    Romn(0) = ( alss * iRbs(imn,lvol-1) + blss * iRbs(imn,lvol) ) * half
    Romn(1) = ( -one * iRbs(imn,lvol-1) +  one * iRbs(imn,lvol) ) * half
    else
    Romn(0) = zero
    Romn(1) = zero
    endif ! ! end of if( NOTstellsym ) ; 22 Apr 13;
    if( Igeometry.ge.3 ) then
    Zomn(0) = ( alss * iZbs(imn,lvol-1) + blss * iZbs(imn,lvol) ) * half
    Zomn(1) = ( -one * iZbs(imn,lvol-1) +  one * iZbs(imn,lvol) ) * half
    if( NOTstellsym ) then
    Zemn(0) = ( alss * iZbc(imn,lvol-1) + blss * iZbc(imn,lvol) ) * half
    Zemn(1) = ( -one * iZbc(imn,lvol-1) +  one * iZbc(imn,lvol) ) * half
    else
    Zemn(0) = zero
    Zemn(1) = zero
    endif ! ! end of if( NOTstellsym ) ; 22 Apr 13; 
    else    
    Zomn(0) = zero
    Zomn(1) = zero
    Zemn(0) = zero
    Zemn(1) = zero
    endif ! ! end of if( Igeometry.ge.3 ) ; 22 Apr 13;
    
   endif ! end of if( Lcoordinatesingularity ) ; 20 Feb 13;
   
   dR(0) = dR(0) + Remn(0) * carg           + Romn(0) * sarg
   if( Lvacuumregion ) then ! ! only in this case are the coordinate derivatives required; 22 Apr 13;
   dR(1) = dR(1) + Remn(1) * carg           + Romn(1) * sarg               
   dR(2) = dR(2) + Remn(0) * sarg * ( -mi ) + Romn(0) * carg * ( +mi )
   dR(3) = dR(3) + Remn(0) * sarg * ( +ni ) + Romn(0) * carg * ( -ni )
   endif
   if( Igeometry.ge.3 ) then
   dZ(0) = dZ(0) + Zomn(0) * sarg           + Zemn(0) * carg
   if( Lvacuumregion ) then ! ! only in this case are the coordinate derivatives required; 22 Apr 13;
   dZ(1) = dZ(1) + Zomn(1) * sarg           + Zemn(1) * carg               
   dZ(2) = dZ(2) + Zomn(0) * carg * ( +mi ) + Zemn(0) * sarg * ( -mi )
   dZ(3) = dZ(3) + Zomn(0) * carg * ( -ni ) + Zemn(0) * sarg * ( +ni )
   endif
   endif
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  enddo ! end of do imn = 1, mn ; Fourier summation loop;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RpZ(1:3) = (/ dR(0) ,  stz(3) , dZ(0) /) ! return coordinate functions;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( Lplasmaregion ) goto 9999 ! ! metric information is not needed; this routine is called ONLY by Poincare plot; 22 Apr 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! the following is only required for the vacuum region; 20 Apr 13;
  
#ifdef DEBUG
  FATALMESS(co00aa, Lplasmaregion, the coordinate derivatives are only required for the vacuum region )
#endif
  
  select case( Igeometry )
   
  case( 2 ) ! Cylindrical; 20 Apr 13;
   
   do ii = 1, 3
    do jj = ii, 3 ; guv(ii,jj) = dR(ii) * dR(jj)
    enddo
   enddo
   
   guv(2,2) = guv(2,2) + dR(0) * dR(0)
   guv(3,3) = guv(3,3) + one

  case( 3 ) ! toroidal; 20 Apr 13;
   
   do ii = 1, 3
    do jj = ii, 3 ; guv(ii,jj) = dR(ii) * dR(jj) + dZ(ii) * dZ(jj)
    enddo
   enddo
   
   guv(3,3) = guv(3,3) + ( pknot * dR(0) )**2
   
  case default
   
   FATALMESS(co00aa, .true., selected Igeometry not supported )
   
  end select
  
  do ii = 2, 3
   do jj = 1, ii-1 ; guv(ii,jj) = guv(jj,ii)
   enddo
  enddo
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(co00aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine co00aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
