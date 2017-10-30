!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Use preconditioned conjugate gradient method (proviced by the NAG routine E04DGF) to find minimum of energy functional.

!latex \item The energy functional is described in \verb+pc00ab+.

!latex \end{enumerate} \subsection{relevant input variables} \begin{enumerate}

!latex \item The following input variables control the operation of \verb+E04DGF+
!latex \begin{itemize}
!latex \item \verb+epsilon+ : weighting of ``spectral energy'';
!latex                        see \verb+pc00ab+;
!latex \item \verb+maxstep+ : this is given to \verb+E04DGF+ for the \verb+Maximum Step Length+;
!latex                       for further details see the NAG documentation for \verb+E04DGF+;
!latex \item \verb+maxiter+ : upper limit on derivative calculations used in the conjugate gradient iterations;
!latex \item \verb+verify+ : if \verb+verify+$=1$, then \verb+E04DGF+ will confirm user supplied gradients (provided by \verb+pc00ab+) are correct;
!latex                       this is primarily used for debugging;
!latex                       for further details see the NAG documentation for \verb+E04DGF+;
!latex \end{itemize}

!latex \item Unfortunately, \verb+E04DGF+ seems to compute approximately $3 N$ function evaluations before proceeding to minimize the energy functional,
!latex where there are $N$ degrees of freedom.
!latex I don't know how to turn this off!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine pc00aa( Ngeometricaldof, position, Nvol, mn, ie04dgf ) ! argument list is optional;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, ten

  use numerical, only :

  use fileunits, only : ounit

  use inputlist, only : Wpc00aa, Lconstraint, verify, maxstep, maxiter, forcetol, Energy, ForceErr

  use cputiming, only : Tpc00aa

  use allglobal, only : ncpu, myid, cpus
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS

  INTEGER, intent(in)    :: Nvol, mn, Ngeometricaldof
  REAL   , intent(inout) :: position(0:Ngeometricaldof)
  INTEGER                :: ie04dgf

  LOGICAL                :: LComputeDerivatives, Lexit = .true.
  INTEGER                :: niterations, Iwork(1:Ngeometricaldof+1), iuser(1:2)
  REAL                   :: lEnergy, Gradient(0:Ngeometricaldof), work(1:13*Ngeometricaldof), ruser(1:1)
  CHARACTER              :: smaxstep*34

  external               :: pc00ab

  BEGIN(pc00aa)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( myid.eq.0 ) then
   cput = GETTIME
   write(ounit,'("pc00aa : ", 10x ," : ")')
   write(ounit,1000) cput-cpus, myid, Ngeometricaldof, maxstep, maxiter, verify
  endif
  
1000 format("pc00aa : ",f10.2," : myid=",i3," ; calling E04DGF : Ngeometricaldof="i6" ; maxstep="es10.2" ; maxiter="i9" ; verify=",i3," ;")
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  iuser(1:2) = (/ 0, 0 /) ! this is only used to pass information through to pc00ab; some resolution settings & iteration counter;
  ruser(1:1) = zero       ! this will be assigned in first call to pc00ab;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!  if( Lexit ) then
!   
!   LComputeDerivatives= .false.
!   WCALL(pc00aa,fc02aa,( Ngeometricaldof, position(0:Ngeometricaldof), Gradient(0:Ngeometricaldof), LComputeDerivatives ))
!   
!   if( myid.eq.0 ) then
!    cput = GETTIME
!    write(ounit,'("pc00aa : ", 10x ," : ")')
!    write(ounit,'("pc00aa : ",f10.2," : iterations="2i8" ;      "3x" ; Energy="es23.15" ;      " 13x  " ; ForceErr="es23.15" ;")') cput-cpus, iuser(1:2), Energy, ForceErr
!   endif
!   
!   if( ForceErr.lt.abs(forcetol) ) then ; ie04dgf=0 ; goto 9999 ! force-balance is satisfied;  
!   endif
!   
!  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  select case( verify ) ! supply optional parameters to E04DGF;
   
  case( -1 ) ! no checks; no screen output;
   call E04DKF('Nolist')
   call E04DKF('Print Level = 0')
   call E04DKF('Verify = -1')
  case(  0 ) ! simple check;
  !call E04DKF('Nolist')
  !call E04DKF('Print Level = 0')
   call E04DKF('Verify =  0') ! simple check
  case(  1 ) ! extensive test;
   call E04DKF('Verify =  1') ! extensive test;   
  case default
   FATALMESS(pc00aa, .true., invalid verify supplied on input : please set verify = -1 for no checks; verify = 0 for simple check; or verify = +1 for extensive test;)
  end select
  
  call E04DKF('Iteration Limit = 99999999')
  write(smaxstep,'("Maximum Step Length ="es13.5)')maxstep
  call E04DKF(smaxstep)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  ie04dgf = 1
  
  call E04DGF( Ngeometricaldof, pc00ab, niterations, lEnergy, Gradient(1:Ngeometricaldof), position(1:Ngeometricaldof), &
               Iwork(1:Ngeometricaldof+1), work(1:13*Ngeometricaldof), iuser(1:2), ruser(1:1), ie04dgf )
  
  cput = GETTIME
  
  select case( ie04dgf )    
  case(:-1)    ; if( myid.eq.0 ) write(ounit,'("pc00aa : ",f10.2," : user requested termination      ; ie04dgf=",i3," ;")')cput-cpus,ie04dgf
  case(  0)    ; if( myid.eq.0 ) write(ounit,'("pc00aa : ",f10.2," : success                         ; ie04dgf=",i3," ;")')cput-cpus,ie04dgf
  case(  3)    ; if( myid.eq.0 ) write(ounit,'("pc00aa : ",f10.2," : iteration limit reached         ; ie04dgf=",i3," ;")')cput-cpus,ie04dgf
  case(  4)    ; if( myid.eq.0 ) write(ounit,'("pc00aa : ",f10.2," : step length too small           ; ie04dgf=",i3," ;")')cput-cpus,ie04dgf
  case(  6)    ; if( myid.eq.0 ) write(ounit,'("pc00aa : ",f10.2," : not all minimum conditions met  ; ie04dgf=",i3," ;")')cput-cpus,ie04dgf
  case(  7)    ; if( myid.eq.0 ) write(ounit,'("pc00aa : ",f10.2," : user supplied derivatives error ; ie04dgf=",i3," ;")')cput-cpus,ie04dgf
  case(  8)    ; if( myid.eq.0 ) write(ounit,'("pc00aa : ",f10.2," : initial gradient too small      ; ie04dgf=",i3," ;")')cput-cpus,ie04dgf
  case(  9)    ; if( myid.eq.0 ) write(ounit,'("pc00aa : ",f10.2," : input error                     ; ie04dgf=",i3," ;")')cput-cpus,ie04dgf
  case default ; FATALMESS(pc00aa, .true., E04DGF ifail error)
  end select
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(pc00aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine pc00aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
