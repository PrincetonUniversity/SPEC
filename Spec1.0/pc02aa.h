!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Use modified Newton algorithm for minimizing function.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine pc02aa( Nvol, mn, Ngeometricaldof, position, ie04lyf ) ! argument list is optional;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero
  use numerical, only :
  use fileunits, only : ounit
  use inputlist, only : Wpc02aa
  use cputiming, only : Tpc02aa
  use allglobal, only : ncpu, myid, cpus
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS

  INTEGER, intent(in)    :: Nvol, mn, Ngeometricaldof
  REAL   , intent(inout) :: position(0:Ngeometricaldof)
  INTEGER                :: ie04lyf

  INTEGER                :: Ibound, intwork(1:Ngeometricaldof), Lintwork, Lrealwork, Iuser(1:4)
  REAL                   :: Energy, Gradient(1:Ngeometricaldof), BL(1:Ngeometricaldof), BU(1:Ngeometricaldof), Ruser(0:2)
  
  REAL, allocatable      :: realwork(:)

  EXTERNAL               :: pc02ab, pc02ac
  
  BEGIN(pc02aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  Iuser(1:4) = (/ mn , Nvol, 0, 0 /)
  
  Ibound = 1 ! there are no bounds;
  
  BL(1:Ngeometricaldof) = zero ; BU(1:Ngeometricaldof) = zero ! provide dummy values;
  
  Lintwork = Ngeometricaldof+2 ; Lrealwork = max(Ngeometricaldof*(Ngeometricaldof+7),10)
  
  RALLOCATE(realwork,(1:Lrealwork))
  
  ie04lyf = 1
  
  call E04LYF( Ngeometricaldof, Ibound, pc02ab, pc02ac, BL, BU, position(1:Ngeometricaldof), Energy, Gradient, &
intwork, Lintwork, realwork, Lrealwork, Iuser(1:4), Ruser(0:2), ie04lyf )
  
  DEALLOCATE(realwork)
  
  cput = GETTIME
  
  select case( ie04lyf )    
  case(    0 )    ; if( myid.eq.0 ) write(ounit,'("pc02aa : ",f10.2," : success                          ; ie04lyf=",i3," ;")')cput-cpus,ie04lyf
  case(    1 )    ; if( myid.eq.0 ) write(ounit,'("pc02aa : ",f10.2," : input error                      ; ie04lyf=",i3," ;")')cput-cpus,ie04lyf
  case(    2 )    ; if( myid.eq.0 ) write(ounit,'("pc02aa : ",f10.2," : 50 x N calculations; restart ?   ; ie04lyf=",i3," ;")')cput-cpus,ie04lyf
  case(    3 )    ; if( myid.eq.0 ) write(ounit,'("pc02aa : ",f10.2," : minimum conditions not satisfied ; ie04lyf=",i3," ;")')cput-cpus,ie04lyf
  case( 5: 8 )    ; if( myid.eq.0 ) write(ounit,'("pc02aa : ",f10.2," : maybe a minimum                  ; ie04lyf=",i3," ;")')cput-cpus,ie04lyf
  case(    9 )    ; if( myid.eq.0 ) write(ounit,'("pc02aa : ",f10.2," : large variable                   ; ie04lyf=",i3," ;")')cput-cpus,ie04lyf
  case(   10 )    ; if( myid.eq.0 ) write(ounit,'("pc02aa : ",f10.2," : first derivatives error          ; ie04lyf=",i3," ;")')cput-cpus,ie04lyf
  case(   11 )    ; if( myid.eq.0 ) write(ounit,'("pc02aa : ",f10.2," : second derivatives error         ; ie04lyf=",i3," ;")')cput-cpus,ie04lyf
  case default    ; FATALMESS(pc02aa,.true.,E04LYF ifail error)
  end select
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(pc02aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine pc02aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine pc02ab( Ngeometricaldof, position, lEnergy, Gradient, Iuser, Ruser ) ! argument list is constrained by NAG;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half
  use numerical, only :
  use fileunits, only : ounit
  use inputlist, only : Wpc02aa
  use cputiming, only : Tpc02aa
  use allglobal, only : ncpu, myid, cpus
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER                :: Ngeometricaldof, Iuser(1:4)
  REAL                   :: position(1:Ngeometricaldof), lEnergy, Gradient(1:Ngeometricaldof), Ruser(0:2)
  
  INTEGER                :: mode, nstate
  
  BEGIN(pc02aa)
  
  mode = 2 ; nstate = 1 ! I think these are dummy values; required for consistency; pc00ab is also called by pc00aa via NAG routine E04DGF;

  WCALL(pc02aa,pc00ab,( mode, Ngeometricaldof, position(1:Ngeometricaldof), lEnergy, Gradient(1:Ngeometricaldof), nstate, Iuser(1:4), Ruser(0:2) ))
  
  RETURN(pc02aa)
  
end subroutine pc02ab

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine pc02ac( Ngeometricaldof, position, LowerTriangle, LH, Diagonal, Iuser, Ruser )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero
  use numerical, only :
  use fileunits, only : ounit
  use inputlist, only : Wpc02aa
  use cputiming, only : Tpc02aa
  use allglobal, only : ncpu, myid, cpus
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER                :: Ngeometricaldof, LH, Iuser(1:4)
  REAL                   :: position(1:Ngeometricaldof), LowerTriangle(1:LH), Diagonal(1:Ngeometricaldof), Ruser(0:2)
  
  INTEGER                :: mode, nstate
  REAL                   :: lEnergy, Gradient(1:Ngeometricaldof)

  BEGIN(pc02aa)

  mode = 3

  WCALL(pc02aa,pc00ab,( mode, Ngeometricaldof, position(1:Ngeometricaldof), lEnergy, Gradient(1:Ngeometricaldof), nstate, Iuser(1:4), Ruser(0:2) )) ! matrix must be returned through global;
    
  RETURN(pc02aa)
  
end subroutine pc02ac

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
