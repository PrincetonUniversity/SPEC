!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Use accelerated steepest descent algorithm to find minimum of energy functional.

!latex \subsubsection{nonlinear conjugate gradient method}
!latex \item The basic algorithm, as described by \cite{Hirshman_Breslau_98}, is 
!latex \be {\bf v}_{n+1} &=& {\bf v}_{n}(1-b) + \Delta t \; {\bf F}({\bf x}_n),\\
!latex     {\bf x}_{n+1} &=& {\bf x}_{n} + \Delta t \; {\bf v}_{n+1},
!latex \ee
!latex where ${\bf v}$ is the ``velocity'' vector, $\Delta t$ is a ``time step'' parameter, and $b$ is a small ``viscous-damping'' parameter, chosen
!latex \be 1 - b \approx \frac{F_n^2}{F_{n-1}^2}.
!latex \ee

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine pc01aa( Nvol, mn, Ngeometricaldof, position ) ! argument list is optional;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half, one
  use numerical, only :
  use fileunits, only : ounit
  use inputlist, only : Wpc01aa, maxiter, maxstep, forceerr, forcetol
  use cputiming, only : Tpc01aa
  use allglobal, only : ncpu, myid, cpus
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS

  INTEGER, intent(in)    :: Nvol, mn, Ngeometricaldof
  REAL   , intent(inout) :: position(1:Ngeometricaldof)
  
  REAL                   :: velocity(1:Ngeometricaldof), Gradient(1:Ngeometricaldof)

  INTEGER                :: mode, nstate, Iuser(1:4)
  REAL                   :: Energy, Ruser(0:2)

  INTEGER                :: itime
  REAL                   :: dtime, damping, FFo, FFn

  REAL                   :: oEnergy, oforceerr, lEnergy, lforceerr

  BEGIN(pc01aa)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{implementation details}

  damping = zero !latex \item The viscous-damping parameter is initialized $b=0$.

  dtime = maxstep !latex \item The time-step is given on input, $\Delta t \equiv $ \verb+maxstep+.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  mode = 2

  Iuser(1:4) = (/ mn, Nvol, 0, 0 /) ; Ruser(0:2) = zero

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  nstate = 1 ! indicates initial call;

  velocity(1:Ngeometricaldof) = zero ! initial velocity;

  FFn = one ; lEnergy = 9.9E+09

  if( myid.eq.0 ) then 
   cput = GETTIME ; write(ounit,'("pc01aa : ",f10.2," : myid=",i3," : maxiter="i6" ; maxstep="es13.5" ;")')cput-cpus,myid,maxiter,maxstep
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  itime = 0

  do ; itime = itime + 1 !latex \item The maximum iterations is given by the input variable \verb+maxiter+. Termination is controlled by \verb+pc00ab+.

   WCALL(pc01aa,pc00ab,( mode, Ngeometricaldof, position(1:Ngeometricaldof), Energy, Gradient(1:Ngeometricaldof), nstate, Iuser(1:4), Ruser(0:2) ))
   
   if( nstate.eq.1 ) then ; oEnergy = Energy ; oforceerr = forceerr ; nstate = 0
   endif

   if( mode.lt.0 ) exit ! pc00ab has returned a signal; either convergence is reached or maximum iterations is reached;

   FFo = FFn
   FFn = sum( Gradient(1:Ngeometricaldof)*Gradient(1:Ngeometricaldof) )

   velocity(1:Ngeometricaldof) = velocity(1:Ngeometricaldof) * ( one-damping ) - dtime * Gradient(1:Ngeometricaldof)

   position(1:Ngeometricaldof) = position(1:Ngeometricaldof) + dtime * velocity(1:Ngeometricaldof)

   if( itime.gt.2 ) damping = one - FFn / FFo

   if(  Energy.gt.lEnergy ) dtime = dtime * half

   if( myid.eq.0 ) then
    write(ounit,'("pc01aa : "10x" ; myid=",i3," ; "i9" ; Energy="es23.15" ; dtime="es23.15" ; damping="es23.15" ;")')myid,itime,Energy,dtime,damping
   endif

   lEnergy = Energy

   if( dtime.lt.abs(forcetol) ) exit

  enddo

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( myid.eq.0 ) then 
   cput = GETTIME ; write(ounit,1000)cput-cpus,myid,mode,oEnergy,Energy,Energy-oEnergy,oforceerr,forceerr,forceerr-oforceerr
  endif
  
1000 format("pc01aa : ",f10.2," : myid=",i3," : mode=",i3," ; Energy="es13.5" -->"es13.5", "es13.5" ; forceerr="es13.5" -->"es13.5", "es13.5" ;")

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(pc01aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine pc01aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
