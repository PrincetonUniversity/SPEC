!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Allocates, and deallocates, geometric information.

!latex \end{enumerate} \subsection{Integrals of basis functions} \begin{enumerate}

!latex \item The integrals of the Chebyshev polynomials and the Fourier harmonics are given in 
!latex
!latex \verb+TTee(0:6,0:Lrad(lvol),0:Lrad(lvol),1:mn,1:mn)+
!latex
!latex \verb+TTeo(0:6,0:Lrad(lvol),0:Lrad(lvol),1:mn,1:mn)+
!latex
!latex \verb+TToe(0:6,0:Lrad(lvol),0:Lrad(lvol),1:mn,1:mn)+
!latex
!latex \verb+TToo(0:6,0:Lrad(lvol),0:Lrad(lvol),1:mn,1:mn)+

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine ma00ab( AD, lvol )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  use constants, only : zero

  use numerical, only :

  use fileunits, only : ounit

  use inputlist, only : Wmacros, Wma00ab, Lrad

  use cputiming, only : Tma00ab

  use allglobal, only : myid, cpus, mn, TTee, TTeo, TToe, TToo, Te, To, Lvacuumregion

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS
  
  CHARACTER, intent(in) :: AD ! this flag indicates whether the derivative array is to be allocated aswell; 16 Jan 13;
  INTEGER  , intent(in) :: lvol
  
  BEGIN(ma00ab)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  select case( AD )
   
  case( 'A' )
   
   RALLOCATE(TTee,(0:6,0:Lrad(lvol),0:Lrad(lvol),1:mn,1:mn))
   RALLOCATE(TTeo,(0:6,0:Lrad(lvol),0:Lrad(lvol),1:mn,1:mn))
   RALLOCATE(TToe,(0:6,0:Lrad(lvol),0:Lrad(lvol),1:mn,1:mn))
   RALLOCATE(TToo,(0:6,0:Lrad(lvol),0:Lrad(lvol),1:mn,1:mn))

   if( Lvacuumregion ) then
   RALLOCATE(Te,(2:6,0:Lrad(lvol),1:mn))
   RALLOCATE(To,(2:6,0:Lrad(lvol),1:mn))
   endif

  case( 'D' )
   
   DEALLOCATE(TTee)
   DEALLOCATE(TTeo)
   DEALLOCATE(TToe)
   DEALLOCATE(TToo)

   if( Lvacuumregion ) then
   DEALLOCATE(Te)
   DEALLOCATE(To)
   endif

  case default
   
   FATALMESS(ma00ab,.true.,invalid allocate or deallocate flag)
   
  end select
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(ma00ab)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine ma00ab

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
