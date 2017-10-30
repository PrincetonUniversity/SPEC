!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Packs, and unpacks, Beltrami field solution vector into vector potential.

!latex \end{enumerate} \subsubsection{construction of ``vector'' of independent degrees of freedom} \begin{enumerate}

!latex \item Numerical routines for solving linear equations, etc. typically require the unknown, independent degrees of freedom,
!latex which in this case are the covariant components of the vector potential, $A_{\t,e,i,l}$ and $A_{\z,e,i,l}$,
!latex to be ``packed'' into a vector, ${\bf x}$.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine up00aa( packorunpack, lvol, NN, solution, dpsi, ideriv ) 

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  use constants, only : zero, half, one
  
  use numerical, only : vsmall, small
  
  use fileunits, only : ounit
  
  use inputlist, only : Wup00aa, Nvol, Lrad
  
  use cputiming, only : Tup00aa

  use allglobal, only : myid, ncpu, cpus, mn, im, in, Mvol, &
                        Ate, Aze, Ato, Azo, Fso, Fse, &
                        NOTstellsym, &
                        Lcoordinatesingularity, Lplasmaregion, Lvacuumregion

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS
  
  CHARACTER, intent(in) :: packorunpack
  INTEGER  , intent(in) :: lvol, NN, ideriv
  REAL                  :: solution(1:NN), dpsi(1:2)
  
  INTEGER               :: ii, ll, mi, ni, idof, pmin, pmax, pskip, llmin, lLrad
  
  REAL                  :: dtf, dpf, Atom, Azon, defaultvalue
  
  BEGIN(up00aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  FATALMESS(up00aa, lvol.lt.1 .or. lvol.gt.Mvol, illegal lvol)
  FATALMESS(up00aa, ideriv.lt.-1 .or. ideriv.gt.2, illegal ideriv)
#endif
  
  defaultvalue = -9.999E+09 ; lLrad = Lrad(lvol)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  select case( packorunpack ) 
   
  case( 'U' )
   
! first, provide "default" value for output; will check below that all harmonics have been assigned; 26 Feb 13;
   
   if( Lplasmaregion ) then
    
#ifdef DEBUG
    do ii = 1, mn
     ;Ate(lvol,ideriv,ii)%s(0:lLrad) = defaultvalue
     ;Aze(lvol,ideriv,ii)%s(0:lLrad) = defaultvalue
     if( NOTstellsym ) then
      Ato(lvol,ideriv,ii)%s(0:lLrad) = defaultvalue
      Azo(lvol,ideriv,ii)%s(0:lLrad) = defaultvalue
     endif
    enddo
#endif
    
    do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
     
     do ll = 0, lLrad - 2
      
      ; idof = Ate(lvol,0,ii)%i(ll)
#ifdef DEBUG
      FATALMESS(up00aa, idof.lt.1 .or. idof.gt.NN, subscript error)
#endif
      ; Ate(lvol,ideriv,ii)%s(ll) = solution(idof)
      
      ; idof = Aze(lvol,0,ii)%i(ll)
#ifdef DEBUG
      FATALMESS(up00aa, idof.lt.1 .or. idof.gt.NN, subscript error)
#endif
      ; Aze(lvol,ideriv,ii)%s(ll) = solution(idof)
      
      if( NOTstellsym ) then
       
       if( ii.eq.1 ) then ! this is the (m,n)=(0,0) sine harmonic; 21 Apr 13;
        Ato(lvol,ideriv,ii)%s(ll) = zero
        Azo(lvol,ideriv,ii)%s(ll) = zero
       else
        
        idof = Ato(lvol,0,ii)%i(ll)
#ifdef DEBUG
        FATALMESS(up00aa, idof.lt.1 .or. idof.gt.NN, subscript error)
#endif
        Ato(lvol,ideriv,ii)%s(ll) = solution(idof)
        
        idof = Azo(lvol,0,ii)%i(ll)
#ifdef DEBUG
        FATALMESS(up00aa, idof.lt.1 .or. idof.gt.NN, subscript error)
#endif
        Azo(lvol,ideriv,ii)%s(ll) = solution(idof)
        
       endif ! end of if( ii.eq.1 ) ; 15 Jan 15;
       
      endif ! end of if( NOTstellsym ) ; 15 Jan 15;
      
     enddo ! end of do ll; 25 Jan 13;
     
    !if( Lcoordinatesingularity .and. ( mi.ne.0 .or. ni.eq.0 ) ) then ! 15 Jan 15; ! additional freedom; 15 Jan 15;
     if( Lcoordinatesingularity .and. ( mi.ne.0 .or. ii.eq.1 ) ) then ! 15 Jan 15; ! additional freedom; 15 Jan 15;
      
      ll = lLrad - 1
      
      ; idof = Aze(lvol,0,ii)%i(ll) ; Aze(lvol,ideriv,ii)%s(ll) = solution(idof)
      if( NOTstellsym ) then
       if( ii.eq.1 ) then
        FATALMESS(up00aa, .true., have not set idof) ! 16 Jan 15;
        ;                           ; Azo(lvol,ideriv,ii)%s(ll) = solution(idof) ! IT SEEMS THAT idof HAS NOT BEEN SET; 15 Jan 15;
       else
        idof = Azo(lvol,0,ii)%i(ll) ; Azo(lvol,ideriv,ii)%s(ll) = solution(idof)
       endif ! end of if( ii.eq.1 ) ; 16 Jan 15;
      endif ! end of if( NOTstellsym ) ; 16 Jan 15;
      
     endif ! end of if( Lcoordinatesingularity ) ; 20 Feb 13;
     
    enddo ! end of do ii; 28 Jan 13;
    
   endif ! ! end of if( Lplasmaregion ) ; 22 Apr 13;
   
   
   
   if( Lvacuumregion ) then
    
    ;                 Ate(lvol,ideriv, 1)%s(0:lLrad) = zero ! sin( 0 t - 0 z ) is irrelevant; 10 Apr 13;
    if( NOTstellsym ) Ato(lvol,ideriv, 1)%s(0:lLrad) = zero ! sin( 0 t - 0 z ) is irrelevant; 10 Apr 13;
    
    do ii = 2, mn ; mi = im(ii) ; ni = in(ii) ! WHY IS THIS SUMMATION FROM ii = 2 ; NON-STELLARATOR-SYMMETRIC TERM NEEDS ii = 1; 15 Jan 15;
     
     do ll = 0, lLrad
      
      ;idof = Ate(lvol,0,ii)%i(ll) ; Ate(lvol,ideriv,ii)%s(ll) = solution(idof)
      if( NOTstellsym ) then
       idof = Ato(lvol,0,ii)%i(ll) ; Ato(lvol,ideriv,ii)%s(ll) = solution(idof)
      endif
      
     enddo ! ! end of do ll; 22 Apr 13;
     
    enddo ! ! end of do ii; 22 Apr 13;
    
   endif ! end of if( Lvacuumregion ) ; 01 May 13;
   
   
  case( 'P' ) ! packing vector-potential into vector; 26 Feb 13;
   

! first, provide "default" value for output; will check below that all harmonics have been assigned; 26 Feb 13;
   
#ifdef DEBUG
   FATALMESS(up00aa, ideriv.ne.0, only the vector potential is ever packed into one-dimensional array) ! 20 Jun 14;
   solution(1:NN) = defaultvalue
#endif
   
   FATALMESS(up00aa, Lvacuumregion, packing vacuum field under construction)
   
   do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
    
    do ll = 0, lLrad - 2
     
     ; idof = Ate(lvol,0,ii)%i(ll) ; solution(idof) = Ate(lvol,ideriv,ii)%s(ll)
     ; idof = Aze(lvol,0,ii)%i(ll) ; solution(idof) = Aze(lvol,ideriv,ii)%s(ll)
     if( NOTstellsym ) then
      if( ii.gt.1 ) then
       idof = Ato(lvol,0,ii)%i(ll) ; solution(idof) = Ato(lvol,ideriv,ii)%s(ll)
       idof = Azo(lvol,0,ii)%i(ll) ; solution(idof) = Azo(lvol,ideriv,ii)%s(ll)
      endif
     endif
     
    enddo ! end of do ll; 25 Jan 13;
    
   !if( Lcoordinatesingularity .and. ( mi.ne.0 .or. ni.eq.0 ) ) then ! 15 Jan 15; ! 15 Jan 15; ! additional freedom; 15 Jan 15;
    if( Lcoordinatesingularity .and. ( mi.ne.0 .or. ii.eq.1 ) ) then ! 15 Jan 15; ! 15 Jan 15; ! additional freedom; 15 Jan 15;
     
     ll = lLrad - 1
     
     ; idof = Aze(lvol,0,ii)%i(ll) ; solution(idof) = Aze(lvol,ideriv,ii)%s(ll)
     if( NOTstellsym ) then
      if( ii.gt.1 ) then
       idof = Azo(lvol,0,ii)%i(ll) ; solution(idof) = Azo(lvol,ideriv,ii)%s(ll)
      endif
     endif
     
    endif ! end of if( Lcoordinatesingularity ) ; 20 Feb 13;
    
   enddo ! end of do ii; 28 Jan 13;
   
   do ii = 2, mn ; mi = im(ii) ; ni = in(ii) ! "pack" surface potential; 26 Feb 13; ! WHY IS THIS SUMMATION FROM ii = 2 ; NON-STELLARATOR-SYMMETRIC TERM NEEDS ii = 1; 15 Jan 15;
    
    if( mi.ne.0 .and. ni.eq.0 ) then ; idof = Fso(lvol,ii) ; solution(idof) =   sum(Ate(lvol,ideriv,ii)%s(0:lLrad)) / mi
    endif
    if( mi.eq.0 .and. ni.ne.0 ) then ; idof = Fso(lvol,ii) ; solution(idof) = - sum(Aze(lvol,ideriv,ii)%s(0:lLrad)) / ni
    endif
    if( mi.ne.0 .and. ni.ne.0 ) then ; Atom =   sum(Ate(lvol,ideriv,ii)%s(0:lLrad)) / mi
     ;                               ; Azon = - sum(Aze(lvol,ideriv,ii)%s(0:lLrad)) / ni
     ;                               ; idof = Fso(lvol,ii) ; solution(idof) = ( Atom + Azon ) * half
    endif
    
    if( NOTstellsym ) then
     
     if( mi.ne.0 .and. ni.eq.0 ) then ; idof = Fse(lvol,ii) ; solution(idof) = - sum(Ato(lvol,ideriv,ii)%s(0:lLrad)) / mi
     endif
     if( mi.eq.0 .and. ni.ne.0 ) then ; idof = Fse(lvol,ii) ; solution(idof) =   sum(Azo(lvol,ideriv,ii)%s(0:lLrad)) / ni
     endif
     if( mi.ne.0 .and. ni.ne.0 ) then ; Atom = - sum(Ato(lvol,ideriv,ii)%s(0:lLrad)) / mi
      ;                               ; Azon =   sum(Azo(lvol,ideriv,ii)%s(0:lLrad)) / ni
      ;                               ; idof = Fse(lvol,ii) ; solution(idof) = ( Atom + Azon ) * half
     endif
    endif
    
   enddo
   
#ifdef DEBUG
   do idof = 1, NN
    if( abs(solution(idof)-defaultvalue).lt.vsmall ) write(ounit,'("up00aa : ", 10x ," : myid=",i3," ; lvol=",i3," ; xi("i9" ) not assigned ;")') myid, lvol, idof
   enddo
#endif
   
  case default
   
   FATALMESS(up00aa, .true., selected packorunpack not supported)
   
  end select ! matches select case( packorunpack ) ; 21 Apr 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Lvacuumregion ) goto 9999 ! no further unpacking is required; (i.e. no "dependent" degrees-of-freedom required for boundary constraint); 21 Apr 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  select case( ideriv )
  case( -1 )    ; dtf = zero    ; dpf = zero
  case(  0 )    ; dtf = dpsi(1) ; dpf = dpsi(2)
  case(  1 )    ; dtf = zero    ; dpf = zero    
  case(  2 )    ; dtf = zero    ; dpf = one 
  case default ; FATALMESS(up00aa, .true., illegal ideriv)
  end select
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  select case( packorunpack)
   
   
  case( 'U' ) ! need to include "dependent" Chebyshev polynomials that enforce the boundary conditions; 28 Jan 13;
   
   
   if( Lcoordinatesingularity ) then
    
    ii = 1
    
    ;  Ate(lvol,ideriv,ii)%s(lLrad-1) =   half *      dtf         - sum( Ate(lvol,ideriv,ii)%s(1:lLrad-3:2) )
    ;  Ate(lvol,ideriv,ii)%s(lLrad  ) =   half *      dtf         - sum( Ate(lvol,ideriv,ii)%s(0:lLrad-2:2) )
    ;  Aze(lvol,ideriv,ii)%s(lLrad  ) =               dpf         - sum( Aze(lvol,ideriv,ii)%s(0:lLrad-1:1) ) ! should there be factor of half?; 04 Dec 14;
    
    if( NOTstellsym) then
     ; Ato(lvol,ideriv,ii)%s(lLrad-1) =   zero
     ; Ato(lvol,ideriv,ii)%s(lLrad  ) =   zero
     ; Azo(lvol,ideriv,ii)%s(lLrad  ) =   zero
    endif
    
    do ii = 2, mn ; mi = im(ii) ; ni = in(ii) ! ii = 1 is treated immediately above;
     
     ; Ate(lvol,ideriv,ii)%s(lLrad-1) =   half * mi * solution(Fso(lvol,ii)) - sum( Ate(lvol,ideriv,ii)%s(1:lLrad-3:2) )
     ; Ate(lvol,ideriv,ii)%s(lLrad  ) =   half * mi * solution(Fso(lvol,ii)) - sum( Ate(lvol,ideriv,ii)%s(0:lLrad-2:2) )
    !if( mi.ne.0 .or. ni.eq.0 ) then ! 15 Jan 15; ! additional freedom; 15 Jan 15;
     if( mi.ne.0 .or. ii.eq.1 ) then ! 15 Jan 15; ! additional freedom; 15 Jan 15;
      ;Aze(lvol,ideriv,ii)%s(lLrad  ) = -        ni * solution(Fso(lvol,ii)) - sum( Aze(lvol,ideriv,ii)%s(0:lLrad-1:1) )
     else
      ;Aze(lvol,ideriv,ii)%s(lLrad-1) = - half * ni * solution(Fso(lvol,ii)) - sum( Aze(lvol,ideriv,ii)%s(1:lLrad-3:2) )
      ;Aze(lvol,ideriv,ii)%s(lLrad  ) = - half * ni * solution(Fso(lvol,ii)) - sum( Aze(lvol,ideriv,ii)%s(0:lLrad-2:2) )
     endif
     
     if( NOTstellsym) then
      ;Ato(lvol,ideriv,ii)%s(lLrad-1) = - half * mi * solution(Fse(lvol,ii)) - sum( Ato(lvol,ideriv,ii)%s(1:lLrad-3:2) )
      ;Ato(lvol,ideriv,ii)%s(lLrad  ) = - half * mi * solution(Fse(lvol,ii)) - sum( Ato(lvol,ideriv,ii)%s(0:lLrad-2:2) )
     !if( mi.ne.0 .or. ni.eq.0 ) then ! 15 Jan 15; ! additional freedom; 15 Jan 15;
      if( mi.ne.0 .or. ii.eq.1 ) then ! 15 Jan 15; ! additional freedom; 15 Jan 15;
       Azo(lvol,ideriv,ii)%s(lLrad  ) =          ni * solution(Fse(lvol,ii)) - sum( Azo(lvol,ideriv,ii)%s(0:lLrad-1:1) )
      else
       Azo(lvol,ideriv,ii)%s(lLrad-1) =   half * ni * solution(Fse(lvol,ii)) - sum( Azo(lvol,ideriv,ii)%s(1:lLrad-3:2) )
       Azo(lvol,ideriv,ii)%s(lLrad  ) =   half * ni * solution(Fse(lvol,ii)) - sum( Azo(lvol,ideriv,ii)%s(0:lLrad-2:2) )
      endif
     endif
     
    enddo ! end of do ii; 20 Feb 13;
    
   else ! matches if( Lcoordinatesingularity ) ; 26 Feb 13;
    
    do ll = lLrad - 1, lLrad   
     
     if( ll .eq. lLrad - 1 ) then ; pmin = 1 ; pmax = lLrad-3 ; pskip = 2
     else                         ; pmin = 0 ; pmax = lLrad-2 ; pskip = 2
     endif
     
     ii = 1
     ; Ate(lvol,ideriv,ii)%s(ll) =   half *      dtf           - sum( Ate(lvol,ideriv,ii)%s(pmin:pmax:pskip) )
     ; Aze(lvol,ideriv,ii)%s(ll) =   half *      dpf           - sum( Aze(lvol,ideriv,ii)%s(pmin:pmax:pskip) )
     if( NOTstellsym) then
      ;Ato(lvol,ideriv,ii)%s(ll) =   zero
      ;Azo(lvol,ideriv,ii)%s(ll) =   zero
     endif
     
     do ii = 2, mn ; mi = im(ii) ; ni = in(ii) ! ii = 1 is treated immediately above;
      ;Ate(lvol,ideriv,ii)%s(ll) =   half * mi * solution(Fso(lvol,ii)) - sum( Ate(lvol,ideriv,ii)%s(pmin:pmax:pskip) )
      ;Aze(lvol,ideriv,ii)%s(ll) = - half * ni * solution(Fso(lvol,ii)) - sum( Aze(lvol,ideriv,ii)%s(pmin:pmax:pskip) )
      if( NOTstellsym) then
       Ato(lvol,ideriv,ii)%s(ll) = - half * mi * solution(Fse(lvol,ii)) - sum( Ato(lvol,ideriv,ii)%s(pmin:pmax:pskip) )
       Azo(lvol,ideriv,ii)%s(ll) =   half * ni * solution(Fse(lvol,ii)) - sum( Azo(lvol,ideriv,ii)%s(pmin:pmax:pskip) )
      endif
     enddo ! end of do ii; 25 Jan 13;
     
    enddo ! end of do ll; 25 Jan 13;
    
   endif ! end of if( Lcoordinatesingularity ) ; 26 Feb 13;
   
#ifdef DEBUG
   do ii = 1, mn
    do ll = 0, lLrad
     ;FATALMESS(up00aa, Ate(lvol,ideriv,ii)%s(ll) - defaultvalue .lt. small, error unpacking)
     ;FATALMESS(up00aa, Aze(lvol,ideriv,ii)%s(ll) - defaultvalue .lt. small, error unpacking)
     if( NOTstellsym ) then
      FATALMESS(up00aa, Ato(lvol,ideriv,ii)%s(ll) - defaultvalue .lt. small, error unpacking)
      FATALMESS(up00aa, Azo(lvol,ideriv,ii)%s(ll) - defaultvalue .lt. small, error unpacking)
     endif
    enddo
   enddo
#endif
   
  end select ! end of select case( packorunpack ) ; 26 Feb 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(up00aa)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine up00aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
