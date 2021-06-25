!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (force-driver) ! Solves ${\bf F}({\bf \xi})=0$, where ${\bf F} \equiv \{ [[p+B^2/2]]_{i,l}, I_{i,l} \}$ and ${\bf \xi} \equiv \{ R_{i,l},Z_{i,l} \}$.

!latex \briefly{Employs Descent method to find ${\bf F}({\bf x})=0$, where ${\bf x}\equiv\{\mbox{\rm geometry}\}$ and ${\bf F}$ is defined in \link{dforce}.}

!latex \calledby{\link{xspech}}
!latex \calls{\link{dforce}, \link{packxi} and \link{global}:wrtend}

!latex \tableofcontents

!latex \subsubsection{iterative, reverse communication loop}

!latex \begin{enumerate}
!latex \item The iterative, Newton search to find ${\bf x} \equiv \{ \mbox{\rm geometry} \} \equiv \{ R_{i,l}, Z_{i,l} \}$ such that ${\bf F}({\bf x})=0$,
!latex       where ${\bf F}$ and its derivatives, $\nabla_{{\bf x}} {\bf F}$, are calculated by \link{dforce}, is provided by either
!latex \begin{itemize} 
!latex \item[i.]  \nag{www.nag.co.uk/numeric/FL/manual19/pdf/C05/c05ndf_fl19.pdf}{C05NDF} if \inputvar{Lfindzero=1},
!latex            which only uses function values; or
!latex \item[ii.] \nag{www.nag.co.uk/numeric/FL/manual19/pdf/C05/c05pdf_fl19.pdf}{C05PDF} if \inputvar{Lfindzero=2}, 
!latex            which uses user-provided derivatives.
!latex \end{itemize}
!latex \item The iterative search will terminate when the solution is within \inputvar{c05xtol} of the true solution (see NAG documentation).
!latex \item The input variable \inputvar{c05factor} is provided to determine the initial step bound (see NAG documentation).
!latex \end{enumerate}

!latex \subsubsection{logic, writing/reading from file}

!latex \begin{enumerate}
!latex \item Before proceeding with iterative search, \link{dforce} is called to determine the magnitude of the initial force imbalance,
!latex       and if this is less than \inputvar{forcetol} then the iterative search will not be performed.
!latex \item As the iterations proceed, \link{global}:wrtend will be called to save itermediate information (also see \link{xspech}).
!latex \item If the derivative matrix, $\nabla_{{\bf x}} {\bf F}$, is required, i.e. if \inputvar{Lfindzero=2}, and if \inputvar{LreadGF=.true.}
!latex       then the derivative matrix will initially be read from \verb+.ext.sp.DF+, if it exists, or from \verb+.sp.DF+.
!latex \item As the iterations proceed, the derivative matrix will be written to \verb+.ext.sp.DF+. 
!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

module descnttime 
  
  INTEGER :: nFcalls, nDcalls
  REAL    :: lastcpu
  
end module descnttime

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine descnt( NGdof, position, ihybrd )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, one, two, ten

  use numerical, only : sqrtmachprec

  use fileunits, only : ounit

  use inputlist, only : Wmacros, Wdescnt, ext, &
                        Igeometry, & ! only for screen output; 
                        Nvol,                    &
                        Lfindzero, forcetol, c05xmax, c05xtol, c05factor, LreadGF, &
                        Lcheck

  use cputiming, only : Tdescnt

  use allglobal, only : myid, ncpu, cpus, &
                        NOTstellsym, &
                        ForceErr, Energy, &
                        mn, im, in, iRbc, iZbs, iRbs, iZbc, Mvol, &
                        BBe, IIo, BBo, IIe, &
                        LGdof, dFFdRZ, dBBdmp, dmupfdx, hessian, dessian, Lhessianallocated , &
                        Ldescent, &
                        nfreeboundaryiterations, &
						            LocalConstraint
  
  use descnttime

  use sphdf5, only: write_convergence_output

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in)    :: NGdof
  REAL   , intent(inout) :: position(0:NGdof)
  INTEGER, intent(out)   :: ihybrd
  
  LOGICAL                :: LComputeDerivatives
  INTEGER                :: wflag, iflag, idof, jdof, ijdof, ireadhessian, igdof, lvol, ii, imn, ierr2
  REAL                   :: rflag
  CHARACTER              :: pack

  INTEGER                :: irevcm, mode, Ldfjac, LR
  REAL                   :: xtol, epsfcn, factor
  REAL                   :: diag(1:NGdof), QTF(1:NGdof), workspace(1:NGdof,1:4)
  REAL                   :: sdxtol=1e-8, sdgtol=1e-20, sdftol=1e-16
  INTEGER                :: sditmax=200, sdiprint=0, sdiflag=0


  REAL                   :: force(0:NGdof)
  REAL, allocatable      :: fjac(:,:), RR(:), work(:,:)
  
  INTEGER                :: ML, MU ! required for only Lc05ndf;
  
  LOGICAL                :: Lexit = .true. ! perhaps this could be made user input;
  LOGICAL                :: LComputeAxis

  INTEGER                :: nprint = 1, nfev, njev, niter

  INTEGER, parameter     :: maxfev = 5000 ! maximum calls per iteration;

  external               :: fcn1, fcn2!, fcnval, fcngrad
    interface
    function fcnval(n,x)
      integer,intent(in):: n
      real(8),intent(in):: x(n)
      real(8):: fcnval
    end function fcnval
    function fcngrad(n,x)
      integer,intent(in):: n
      real(8),intent(in):: x(n)
      real(8):: fcngrad(n)
    end function fcngrad
  end interface

  BEGIN(descnt)
 
  Ldescent = .TRUE.
 
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  nFcalls = 0 ; nDcalls= 0 ! counters;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  lastcpu = GETTIME
  
  if( Lexit ) then ! will call initial force, and if ForceErr.lt.forcetol will immediately exit; 

   LComputeDerivatives= .false.
   LComputeAxis = .true.
   WCALL( descnt, dforce, ( NGdof, position(0:NGdof), force(0:NGdof), LComputeDerivatives, LComputeAxis) ) ! calculate the force-imbalance;
   
   if( myid.eq.0 ) then ! screen output; 
    cput = GETTIME
    ; write(ounit,1000) cput-cpus, nFcalls, nDcalls, ForceErr,  cput-lastcpu, "|BB|e", alog10(BBe(1:min(Mvol-1,28)))
    if( Igeometry.ge.3 ) then ! include spectral constraints; 
     ;write(ounit,1001)                                                                       "|II|o", alog10(IIo(1:min(Mvol-1,28)))
    endif
    if( NOTstellsym ) then
     ;write(ounit,1001)                                                                       "|BB|o", alog10(BBo(1:min(Mvol-1,28)))
     if( Igeometry.ge.3 ) then ! include spectral constraints; 
      write(ounit,1001)                                                                       "|II|e", alog10(IIe(1:min(Mvol-1,28)))
     endif
    endif
   endif

   if( ForceErr.lt.forcetol ) then ; ihybrd = 0 ; goto 9999 ! force-balance is satisfied;  
   endif
   
  endif ! end of if( Lexit ) ; 

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

1000 format("descnt : ",f10.2," : "i9,i3," ; ":"|f|="es12.5" ; ":"time=",f10.2,"s ;":" log"a5"="28f6.2" ...")
1001 format("descnt : ", 10x ," : "9x,3x" ; ":"    "  12x "   ":"     ", 10x ,"  ;":" log"a5"="28f6.2" ...")
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
  irevcm = 0 ; ihybrd = 1 ! required for initial entry; herefater unchanged by user;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  SALLOCATE( fjac, (1:NGdof, 1:NGdof), zero)
  SALLOCATE( RR, (1:NGdof*(NGdof+1)/2), zero)

  
  Lhessianallocated = .false.

if(Lfindzero.eq.3) then
    SALLOCATE( dFFdRZ, (1:LGdof,0:1,1:LGdof,0:1,1:Mvol), zero )
    SALLOCATE( dBBdmp, (1:LGdof,1:Mvol,0:1,1:2), zero )
   if( LocalConstraint ) then
        SALLOCATE( dmupfdx, (1:Mvol,    1:1,1:2,1:LGdof,0:1), zero )
   else
        SALLOCATE( dmupfdx, (1:Mvol, 1:Mvol-1,1:2,1:LGdof,1), zero ) ! TODO change the format to put vvol in last index position...
   endif
    SALLOCATE( hessian, (1:NGdof,1:NGdof), zero )
    SALLOCATE( dessian, (1:NGdof,1:LGdof), zero )
    Lhessianallocated = .true.
endif

  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  select case( Lfindzero )

  case( 3 ) ! use only function values to find f(x)=0 by pushing the interfaces with f(x) (force-descent) using VMEC's alogrithm

   write(*,*) "-------------------- Under construction --------------------"

   WCALL(descnt, fcndescent, (position(1:NGdof),NGdof))

   write(*,*) "-------------------- Under construction --------------------"

  case( 4 ) ! use function values to find f(x)=0 by pushing the interfaces with f(x) (force-descent) using Anderson acceleration

   write(*,*) "-------------------- Under construction --------------------"
!    WCALL(descnt, frcg, (NGdof,position(1:NGdof),Energy,sdxtol,sdgtol,sdftol,sditmax &
!      ,sdiprint,sdiflag,fcnval,fcngrad))
   WCALL(descnt, fcnanderson, (position(1:NGdof),NGdof))

   write(*,*) "-------------------- Under construction --------------------"
  
  case default
   
   FATAL( descnt, .true., value of Lfindzero not supported )
   
  end select
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
	call MPI_BARRIER( MPI_COMM_WORLD, ierr2)

  DALLOCATE( fjac )
  DALLOCATE( RR )

if(Lfindzero.eq.3) then
    DALLOCATE( dFFdRZ )
    DALLOCATE( dBBdmp )
   if( LocalConstraint ) then
        DALLOCATE( dmupfdx )
   else   
        DALLOCATE( dmupfdx )
   endif
    DALLOCATE( hessian )
    DALLOCATE( dessian )
    Lhessianallocated = .false.
endif
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  RETURN(descnt)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine descnt

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

REAL function fcnval(n,xx)
!subroutine fcnval (f, xx, n)
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, one, two, ten

  use numerical, only : sqrtmachprec

  use fileunits, only : ounit

  use inputlist, only : Wmacros, Wdescnt, ext, &
                        Igeometry, & ! only for screen output; 
                        Nvol,     &
                        Lfindzero, forcetol, c05xmax, c05xtol, c05factor, LreadGF, &
                        Lcheck, &
                        Lconstraint, mu, epsilon, opsilon

  use cputiming, only : Tdescnt

  use allglobal, only : wrtend, myid, ncpu, cpus, &
                        NOTstellsym, &
                        ForceErr, Energy, &
                        mn, im, in, iRbc, iZbs, iRbs, iZbc, Mvol, &
                        BBe, IIo, BBo, IIe, &
                        LGdof, dFFdRZ, dBBdmp, dmupfdx, hessian, dessian, Lhessianallocated, &
                        nfreeboundaryiterations, &
                        lABintegral, lLLl, lMMl, sweight, pi2, pi2nfp
  
  use descnttime

  use sphdf5, only : write_convergence_output

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS

  !REAL                   :: fcnval 
  INTEGER, intent(in)    :: n
  REAL   , intent(in)    :: xx(1:n)
  !REAL   , intent(out)   :: f

  REAL                   :: position(0:n), force(0:n)

  LOGICAL                :: LComputeDerivatives, Lonlysolution, LComputeAxis
  INTEGER                :: idof, jdof, ijdof, ireadhessian, igdof, lvol, ii, imn
  CHARACTER              :: pack
    

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  position = zero ; force = zero ; position(1:n) = xx(1:n)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
    LComputeDerivatives = .false.
    LComputeAxis = .true.
    nFcalls = nFcalls + 1
    WCALL( descnt, dforce, ( n, position(0:n), force(0:n), LComputeDerivatives, LComputeAxis ) ) ! calculate the force-imbalance;

   ! if(Lconstraint.eq.0) then
   !  fcnval = Energy - sum( mu(1:Nvol)*lABintegral(1:Nvol) )   
   ! else
    !fcnval = Energy
    fcnval = opsilon * Energy / (pi2 * pi2nfp) !+ epsilon * sum(lMMl) !+ sum(lLLl * sweight)
    !write(ounit,*) f
   ! endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

 end function fcnval
 !end subroutine fcnval

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

function fcngrad(n, xx)
!subroutine fcngrad (g, xx, n)
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, one, two, ten

  use numerical, only : sqrtmachprec

  use fileunits, only : ounit

  use inputlist, only : Wmacros, Wdescnt, ext, &
                        Igeometry, & ! only for screen output; 
                        Nvol,                    &
                        Lfindzero, forcetol, c05xmax, c05xtol, c05factor, LreadGF, &
                        Lcheck

  use cputiming, only : Tdescnt

  use allglobal, only : wrtend, myid, ncpu, cpus, &
                        NOTstellsym, &
                        ForceErr, Energy, &
                        mn, im, in, iRbc, iZbs, iRbs, iZbc, Mvol, &
                        BBe, IIo, BBo, IIe, &
                        LGdof, dFFdRZ, dBBdmp, dmupfdx, hessian, dessian, Lhessianallocated, &
                        nfreeboundaryiterations, pi2nfp, pi2, &
                        lMMl, lLLl, sweight
  
  use descnttime

  use sphdf5, only : write_convergence_output

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS

  REAL                   :: fcngrad(1:n) 
  INTEGER, intent(in)    :: n
  REAL   , intent(in)    :: xx(1:n)
  !REAL   , intent(out)   :: g(1:n)

  REAL                   :: position(0:n), force(0:n)

  LOGICAL                :: LComputeDerivatives, Lonlysolution, LComputeAxis
  INTEGER                :: idof, jdof, ijdof, ireadhessian, igdof, lvol, ii, imn
  CHARACTER              :: pack
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  position = zero ; force = zero ; position(1:n) = xx(1:n)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
    LComputeDerivatives = .false.
    LComputeAxis = .true.
    WCALL( descnt, dforce, ( n, position(0:n), force(0:n), LComputeDerivatives, LComputeAxis ) ) ! calculate the force-imbalance;
    nDcalls = nDcalls + 1
    
    fcngrad(1:n) = force(1:n) 
    !g(1:n) = force(1:n) 
    !write(ounit,*) xx(n), force(n)

    if( myid.eq.0 ) then
     
     cput = GETTIME
     
     ; write(ounit,1000) cput-cpus, nFcalls, nDcalls, ForceErr, cput-lastcpu, "|BB|e", alog10(BBe(1:min(Mvol-1,28)))
     if( Igeometry.ge.3 ) then ! include spectral constraints; 
      ;write(ounit,1001)                                                                      "|II|o", alog10(IIo(1:min(Mvol-1,28)))
     endif
     if( NOTstellsym ) then
      ;write(ounit,1001)                                                                      "|BB|o", alog10(BBo(1:min(Mvol-1,28)))
      if( Igeometry.ge.3 ) then ! include spectral constraints; 
       write(ounit,1001)                                                                      "|II|e", alog10(IIe(1:min(Mvol-1,28)))
      endif
     endif
     lastcpu = GETTIME
     WCALL( descnt, wrtend ) ! write restart file; save geometry to ext.end;

    endif ! end of if( myid.eq.0 );

1000 format("fcngrad: ",f10.2," : "i9,i3," ; ":"|f|="es12.5" ; ":"time=",f10.2,"s ;":" log"a5"="28f6.2" ...")
1001 format("fcngrad: ", 10x ," : "9x,3x" ; ":"    "  12x "   ":"     ", 10x ,"  ;":" log"a5"="28f6.2" ...")
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

 end function fcngrad
 !end subroutine fcngrad
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine fcndescent(xx, NGdof)

 use constants, only  : zero, one

 use fileunits, only : ounit

 use inputlist, only  : epsilon, Igeometry, Wdescnt, &
                        dxdesc, ftoldesc, maxitdesc, Lwritedesc, nwritedesc, &
                        Wmacros

 use fileunits, only : ounit

 use allglobal, only  : wrtend, ForceErr, lMMl, lLLl, Energy, sweight, myid, &
                        BBe, BBo, IIe, IIo, NOTstellsym, ncpu, cpus, Mvol, &
                        hessian

 use cputiming, only : Tdescnt

 use descnttime

 use sphdf5, only : write_convergence_output

 LOCALS 

 INTEGER, INTENT(in)  :: NGdof
 REAL , INTENT(inout) :: xx(1:NGdof)

 REAL                 :: time_step=1.0E-3
 INTEGER,parameter    :: ndamp=15

 REAL                 :: position(0:NGdof), force(0:NGdof), precond(1:NGdof)

 REAL                 :: otau(1:ndamp), xdot(1:NGdof)
 REAL                 :: fsq, fsq1, otav, b1, fac, dtau

 INTEGER              :: it, idesc, lvol, idof
 LOGICAL              :: LComputeDerivatives, LComputeAxis
 CHARACTER            :: pack

 position = zero ; force = zero ; position(1:NGdof) = xx(1:NGdof); precond(1:NGdof) = one 

 LComputeDerivatives = .false.; LComputeAxis = .true.

 time_step = dxdesc

 it = 0
 otau(:ndamp) = 0.15/time_step
 fsq = one
 xdot = zero

 ! The following time stepping algorithm is adapted from VMEC

 do while(it < maxitdesc)
 
  it = it + 1
 
  if(it.eq.1) then
   LComputeDerivatives = .true.
  else
   LComputeDerivatives = .false.
  endif

  call dforce(NGdof, position(0:NGdof), force(0:NGdof), LComputeDerivatives, LComputeAxis)
   
  do idof=1,NGdof
  precond(idof) = 1/ABS(hessian(idof,idof))
  enddo

  fsq1 = sum(force(1:NGdof)**2)

  dtau = MIN(ABS(LOG(fsq1/fsq)), 0.15)

  ! update backup copy of fsq1
  fsq = fsq1

  ! shift array for averaging to the left to make space at end for new entry
  otau(1:ndamp-1) = otau(2:ndamp)
  otau(ndamp) = dtau/time_step
    
  ! averaging over ndamp entries : 1/ndamp*sum(otau)
  otav = SUM(otau(:ndamp))/ndamp

  dtau = time_step*otav/2

  b1  = one - dtau
  fac = one/(one + dtau)

!     THIS IS THE TIME-STEP ALGORITHM. IT IS ESSENTIALLY A CONJUGATE
!     GRADIENT METHOD, WITHOUT THE LINE SEARCHES (FLETCHER-REEVES),
!     BASED ON A METHOD GIVEN BY P. GARABEDIAN
!

  xdot = fac*(b1*xdot - time_step*precond(1:NGdof)*force(1:NGdof))                 ! update velocity
  position(1:NGdof) = position(1:NGdof) + time_step*xdot          ! advance xc by velocity given in xcdot


  if(ForceErr<ftoldesc) then
   if(myid.eq.0) then
    write(*,*) "FORCE BELOW TOLERANCE"
   endif
   exit
  endif

  if(it .eq. maxitdesc .and. myid.eq.0) then
   write(*,*) "EXCEEDED MAX NUMBER OF ITERATIONS " , ForceErr
  endif

  if (myid .eq. 0) then
    if (mod(it, nwritedesc) .eq. 0) then

     cput = GETTIME
     
     ; write(ounit,1000) cput-cpus, it, 0, ForceErr, cput-lastcpu, "|BB|e", alog10(BBe(1:min(Mvol-1,28)))
     if( Igeometry.ge.3 ) then ! include spectral constraints; 
      ;write(ounit,1001)                                                                      "|II|o", alog10(IIo(1:min(Mvol-1,28)))
     endif
     if( NOTstellsym ) then
      ;write(ounit,1001)                                                                      "|BB|o", alog10(BBo(1:min(Mvol-1,28)))
      if( Igeometry.ge.3 ) then ! include spectral constraints; 
       write(ounit,1001)                                                                      "|II|e", alog10(IIe(1:min(Mvol-1,28)))
      endif
     endif

     if ( Igeometry .eq. 3) then
      write(ounit,1003)  Energy,  epsilon *  sum(lMMl), sum(lLLl * sweight)
     else
      write(ounit,1002) Energy
     endif

     lastcpu = GETTIME
     WCALL( descnt, wrtend ) ! write restart file; save geometry to ext.end;

     if (Lwritedesc .ge. 1) then
      WCALL( descnt, write_convergence_output, ( nDcalls, ForceErr ) ) ! save iRbc, iZbs consistent with position;
     endif

    endif
  endif
 enddo
 
 xx = position(1:NGdof)

1000 format("descnt : ",f10.2," : "i9,i3," ; ":"|f|="es12.5" ; ":"time=",f10.2,"s ;":" log"a5"="28f6.2" ...")
1001 format("descnt : ", 10x ," : "9x,3x" ; ":"    "  12x "   ":"     ", 10x ,"  ;":" log"a5"="28f6.2" ...")
1002 format("descnt :            : Energy ", es23.15, "  ;")
1003 format("descnt :            : Energy ", es23.15, "  SEnergy  ", es23.15, "  LEnergy  ", es23.15)

end subroutine

!-----------------------------------

subroutine fcnanderson(xx, NGdof)

 use constants, only  : zero, one

 use fileunits, only : ounit

 use inputlist, only  : epsilon, Igeometry, Wdescnt, &
                        dxdesc, ftoldesc, maxitdesc, Lwritedesc, nwritedesc, Manderson, &
                        Wmacros

 use fileunits, only : ounit

 use allglobal, only  : wrtend, ForceErr, lMMl, lLLl, Energy, sweight, myid, &
                        BBe, BBo, IIe, IIo, NOTstellsym, ncpu, cpus, Mvol

 use cputiming, only : Tdescnt

 use descnttime

 use sphdf5, only : write_convergence_output

 LOCALS 

 INTEGER, INTENT(in)  :: NGdof
 REAL , INTENT(inout) :: xx(1:NGdof)

 REAL                 :: position(0:NGdof), force(0:NGdof)

 REAL                 :: W(1:Manderson+1,1:Manderson)
 REAL                 :: G(1:NGdof,1:Manderson+1), Gbar(1:NGdof,1:Manderson)
 REAL                 :: Q(1:NGdof,1:NGdof), R(1:Manderson,1:Manderson)
 REAL                 :: alpha(1:Manderson+1), alphabar(1:Manderson)
 REAL                 :: posref(1:NGdof,1:Manderson+1), minusQg(1:NGdof)
 REAL                 :: reflectors(1:min(NGdof,Manderson)), dummy(1)
 REAL                 :: Sval(1:min(NGdof,Manderson)), Am(1:NGdof,1:Manderson), Bm(1:max(NGdof,Manderson))

 REAL, ALLOCATABLE    :: Wo(:,:), Go(:,:), Gbaro(:,:), Ro(:,:)
 REAL, ALLOCATABLE    :: Svalo(:), Amo(:,:), Bmo(:)
 REAL, ALLOCATABLE    :: alphao(:), alphabaro(:)
 REAL, ALLOCATABLE    :: wk(:)

 INTEGER              :: it, idesc, lvol, idof, j, ii
 INTEGER              :: idginfo, iwa(1), grank, ldb
 REAL                 :: grcond
 LOGICAL              :: LComputeDerivatives, LComputeAxis
 CHARACTER            :: pack

 LComputeDerivatives = .false.; LComputeAxis = .true.

 position = zero ; force = zero ; position(1:NGdof) = xx(1:NGdof) 
  
 !check if force-descent is run with brute force or with anderson acceleration

 if(Manderson .le. 0) then 
! ----- RUN WITHOUT ANDERSON ACCELERATION -----

 it = 0

 do while(it < maxitdesc)

  it = it + 1
 
  call dforce(NGdof, position(0:NGdof), force(0:NGdof), LComputeDerivatives, LComputeAxis)

  position(1:NGdof) = position(1:NGdof) - dxdesc*force(1:NGdof)!/ForceErr

  if(ForceErr<ftoldesc) then
   if(myid.eq.0) then
    write(*,*) "FORCE BELOW TOLERANCE"
   endif
   exit
  endif

  if(it .eq. maxitdesc .and. myid.eq.0) then
   write(*,*) "EXCEEDED MAX NUMBER OF ITERATIONS " , ForceErr
  endif
  
  if (myid .eq. 0) then
  if (mod(it, nwritedesc) .eq. 0) then
    cput = GETTIME

    ; write(ounit,1000) cput-cpus, it, 0, ForceErr, cput-lastcpu, "|BB|e", alog10(BBe(1:min(Mvol-1,28)))
    if( Igeometry.ge.3 ) then ! include spectral constraints; 
    ;write(ounit,1001)                                                                      "|II|o", alog10(IIo(1:min(Mvol-1,28)))
    endif
    if( NOTstellsym ) then
     ;write(ounit,1001)                                                                      "|BB|o", alog10(BBo(1:min(Mvol-1,28)))
     if( Igeometry.ge.3 ) then ! include spectral constraints; 
      write(ounit,1001)                                                                      "|II|e", alog10(IIe(1:min(Mvol-1,28)))
     endif
    endif

    if ( Igeometry .eq. 3) then
     write(ounit,1003)  Energy,  epsilon *  sum(lMMl), sum(lLLl * sweight)
    else
     write(ounit,1002) Energy
    endif

    lastcpu = GETTIME
    WCALL( descnt, wrtend ) ! write restart file; save geometry to ext.end;

    if (Lwritedesc .ge. 1) then
     WCALL( descnt, write_convergence_output, ( nDcalls, ForceErr ) ) ! save iRbc, iZbs consistent with position;
    endif

   endif
 endif
 enddo

 else
 ! ----- RUN WITH ANDERSON ACCELERATION -----

 posref = zero; G = zero; grcond = 1.0E-8; 

 alpha(1:Manderson) = zero; alpha(Manderson+1) = 1; alphabar = zero;

 W = zero; W(1,1) = -1; W(Manderson+1,Manderson) = 1;

 do j=2,Manderson
  W(j,j)   = -1
  W(j,j-1) =  1
 enddo  

 call dforce(NGdof, position(0:NGdof), force(0:NGdof), LComputeDerivatives, LComputeAxis)
  
 posref(1:NGdof,Manderson+1) = position(1:NGdof) - dxdesc*force(1:NGdof)!/ForceErr

 G(1:NGdof,Manderson+1)      = posref(1:NGdof,Manderson+1) - position(1:NGdof)
 
 position(1:NGdof)           = posref(1:NGdof,Manderson+1)
 
 it = 0

 do while(it < maxitdesc)
 
  it = it + 1

  if(it < Manderson) then
   SALLOCATE( Wo,        (1:it+1, 1:it),    zero)
   SALLOCATE( Go,        (1:NGdof, 1:it+1), zero)
   SALLOCATE( Gbaro,     (1:NGdof, 1:it),   zero)
   SALLOCATE( Ro,        (1:NGdof, 1:it),   zero)
   SALLOCATE( Amo,       (1:NGdof, 1:it),   zero)
   SALLOCATE( Bmo,       (1:max(NGdof,it)), zero)
   SALLOCATE( Svalo,     (1:min(NGdof,it)), zero)
   SALLOCATE( alphao,    (1:it+1),          zero)
   SALLOCATE( alphabaro, (1:it),            zero)
   alphao(1:it)=zero; alphao(it+1)=1;
   Wo = zero; Wo(1,1) = -1; Wo(it+1,it)=1;
   do j=2,it
    Wo(j,j)   = -1
    Wo(j,j-1) =  1
   enddo
  endif
  
  call dforce(NGdof, position(0:NGdof), force(0:NGdof), LComputeDerivatives, LComputeAxis)
  
  posref(1:NGdof,1:Manderson) = posref(1:NGdof,2:Manderson+1)

  posref(1:NGdof,Manderson+1) = position(1:NGdof) - dxdesc*force(1:NGdof)!/ForceErr

  G(1:NGdof,1:Manderson)      = G(1:NGdof,2:Manderson+1)

  G(1:NGdof,Manderson+1)      = posref(1:NGdof,Manderson+1) - position(1:NGdof)
 
  if(it >= Manderson) then
   Gbar(1:NGdof,1:Manderson)   = matmul(G(1:NGdof,1:Manderson+1),W(1:Manderson+1,1:Manderson))
   !solve least square linear problem for alphabar
   Am(1:NGdof,1:Manderson) = Gbar(1:NGdof,1:Manderson)  
   Bm(1:NGdof)             = -G(1:NGdof,Manderson+1)
   iwa(1)=-1; idginfo = 0;  ldb = max(NGdof,Manderson);
   call dgelss(NGdof,Manderson,1,Am,NGdof,Bm,ldb,Sval,grcond,grank,dummy,iwa(1),idginfo)    
   iwa(1)=int(dummy(1)); allocate(wk(1:iwa(1)))
   call dgelss(NGdof,Manderson,1,Am,NGdof,Bm,ldb,Sval,grcond,grank,wk,iwa(1),idginfo)  
   deallocate(wk)
   alphabar(1:Manderson) = Bm(1:Manderson)
   !RQ-decomposition of GbarT
   !iwa(1)=-1; idgeqrf = 0
   !call dgerqf(Manderson,NGdof,GbarT,ldgbar,reflectors,dummy,iwa(1),idgeqrf)
   !iwa(1)=int(dummy(1)); allocate(wk(1:iwa(1)))
   !call dgerqf(Manderson,NGdof,GbarT,ldgbar,reflectors,wk,iwa(1),idgeqrf)
   !deallocate(wk)
   !R(1:Manderson,1:Manderson) = zero
   !need to fill R coefficients depending on wether M>N or M<N
   !if(Manderson .le. NGdof) then
   ! R(1:Manderson,1:Manderson) = GbarT(1:Manderson,NGdof-Manderson+1:NGdof)
   !else
   !endif
   !do j=1,Manderson
   ! R(j,j:Manderson) = Gbar(j,j:Manderson) 
   !enddo
   !need to construct the Q matrix
   !Q(1:NGdof,1:NGdof)      = zero
   !mQ                      = min(Manderson, NGdof)
   !do j=1,mQ
   ! Q(1:NGdof,j) = Gbar(1:NGdof,j)    
   !enddo
   !ii=-1; idgeqrf = 0
   !call dorgqr(NGdof,NGdof,mQ,Q,ldgbar,reflectors,dummy,ii,idgeqrf)
   !ii=int(dummy(1)); allocate(wk(1:ii))
   !call dorgqr(NGdof,NGdof,mQ,Q,ldgbar,reflectors,wk,ii,idgeqrf)
   !deallocate(wk)
   !end of RQ-decomposition
   !solve linear system for alphabar
   !minusQg(1:NGdof)    = -matmul(transpose(Q(1:NGdof,1:NGdof)),G(1:NGdof,Manderson+1))
   !alphabar(Manderson) = minusQg(mQ)/R(mQ,Manderson) ! DON'T UNDERSTAND HOW THIS WORKS...
   !do j=Manderson-1,Manderson-mQ,-1
   ! alphabar(j) = (minusQg(j)-dot_product(R(j,j+1:mQ),alphabar(j+1:mQ)))/R(j,j)
   !enddo
   !solve for alpha
   alpha(1:Manderson+1)        = matmul(W(1:Manderson+1,1:Manderson),alphabar(1:Manderson))
   alpha(Manderson+1)          = alpha(Manderson+1) + 1
   !update position according to anderson's rule
   do j=1,NGdof 
    position(j)           = dot_product(posref(j,1:Manderson+1),alpha(1:Manderson+1))
   enddo 
  else
   Go(1:NGdof,1:it+1)  = G(1:NGdof,Manderson+1-it:Manderson+1)  
   Gbaro(1:NGdof,1:it) = matmul(Go(1:NGdof,1:it+1),Wo(1:it+1,1:it))
   !solve least square linear problem for alphabar
   Amo(1:NGdof,1:it) = Gbaro(1:NGdof,1:it) 
   Bmo(1:NGdof)      = -Go(1:NGdof,it+1)
   iwa(1)=-1; idginfo = 0;  ldb = max(NGdof,it);
   call dgelss(NGdof,it,1,Amo,NGdof,Bmo,ldb,Svalo,grcond,grank,dummy,iwa(1),idginfo)
   iwa(1)=int(dummy(1)); allocate(wk(1:iwa(1)))
   call dgelss(NGdof,it,1,Amo,NGdof,Bmo,ldb,Svalo,grcond,grank,wk,iwa(1),idginfo)
   deallocate(wk)
   alphabaro(1:it) = Bmo(1:it)
   !QR-decomposition of Gbaro
   !iwa(1)=-1; idgeqrf = 0
   !call dgeqrf(NGdof,it,Gbaro,ldgbar,wa1,wa2,iwa(1),idgeqrf)
   !iwa(1)=int(wa2(1)); allocate(wk(1:iwa(1)))
   !call dgeqrf(NGdof,it,Gbaro,ldgbar,reflectors,wk,iwa(1),idgeqrf)
   !deallocate(wk)
   !Ro(1:NGdof,1:it) = zero
   !do j=1,it
   ! Ro(j,j:it) = Gbaro(j,j:it)
   !enddo
   !Q(1:NGdof,1:NGdof)      = zero
   !do j=1,it
   ! Q(1:NGdof,j) = Gbaro(1:NGdof,j)
   !enddo
   !ii=-1; idgeqrf = 0
   !call dorgqr(NGdof,NGdof,it,Q,ldgbar,reflectors,dummy,ii,idgeqrf)
   !ii=int(dummy(1)); allocate(wk(1:ii))
   !call dorgqr(NGdof,NGdof,it,Q,ldgbar,reflectors,wk,ii,idgeqrf)
   !deallocate(wk)
   !end of QR decomposition
   !solve linear system for alphabaro
   !minusQg(1:NGdof)    = -matmul(transpose(Q(1:NGdof,1:NGdof)),Go(1:NGdof,it+1))
   !alphabaro(it)       = minusQg(it)/Ro(it,it)
   !do j=it-1,1,-1
   !alphabaro(j) = (minusQg(j)-dot_product(Ro(j,j+1:it),alphabaro(j+1:it)))/Ro(j,j)
   !enddo
   !solve for alphao
   alphao(1:it+1)   = matmul(Wo(1:it+1,1:it),alphabaro(1:it))
   alphao(it+1) = alphao(it+1)+1
   !update position according to anderson's rule
   do j=1,NGdof
   position(j)           = dot_product(posref(j,Manderson+1-it:Manderson+1),alphao(1:it+1))
   enddo
  endif

  if(ForceErr<ftoldesc) then
   if(myid.eq.0) then
    write(*,*) "FORCE BELOW TOLERANCE"
   endif
   exit
  endif

  if(it .eq. maxitdesc .and. myid.eq.0) then
   write(*,*) "EXCEEDED MAX NUMBER OF ITERATIONS, Force = " , ForceErr
  endif

  if (myid .eq. 0) then
    if (mod(it, nwritedesc) .eq. 0) then

     cput = GETTIME
     
     ; write(ounit,1000) cput-cpus, it, 0, ForceErr, cput-lastcpu, "|BB|e", alog10(BBe(1:min(Mvol-1,28)))
     if( Igeometry.ge.3 ) then ! include spectral constraints; 
      ;write(ounit,1001)                                                                      "|II|o", alog10(IIo(1:min(Mvol-1,28)))
     endif
     if( NOTstellsym ) then
      ;write(ounit,1001)                                                                      "|BB|o", alog10(BBo(1:min(Mvol-1,28)))
      if( Igeometry.ge.3 ) then ! include spectral constraints; 
       write(ounit,1001)                                                                      "|II|e", alog10(IIe(1:min(Mvol-1,28)))
      endif
     endif

     if ( Igeometry .eq. 3) then
      write(ounit,1003)  Energy,  epsilon *  sum(lMMl), sum(lLLl * sweight)
     else
      write(ounit,1002) Energy
     endif

     lastcpu = GETTIME
     WCALL( descnt, wrtend ) ! write restart file; save geometry to ext.end;

     if (Lwritedesc .ge. 1) then
      WCALL( descnt, write_convergence_output, ( nDcalls, ForceErr ) ) ! save iRbc, iZbs consistent with position;
     endif

    endif
  endif

  if(it < Manderson) then
   DALLOCATE( Wo )
   DALLOCATE( Go )
   DALLOCATE( Gbaro )
   DALLOCATE( Ro )
   DALLOCATE( Bmo )
   DALLOCATE( Amo )
   DALLOCATE( Svalo )
   DALLOCATE( alphao )
   DALLOCATE( alphabaro )
  endif

 enddo

 endif !end of if(Manderson .le. 0)

 xx = position(1:NGdof)
1000 format("descnt : ",f10.2," : "i9,i3," ; ":"|f|="es12.5" ; ":"time=",f10.2,"s ;":" log"a5"="28f6.2" ...")
1001 format("descnt : ", 10x ," : "9x,3x" ; ":"    "  12x "   ":"     ", 10x ,"  ;":" log"a5"="28f6.2" ...")
1002 format("descnt :            : Energy ", es23.15, "  ;")
1003 format("descnt :            : Energy ", es23.15, "  SEnergy  ", es23.15, "  LEnergy  ", es23.15)

end subroutine
