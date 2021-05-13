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
! if there is no PETSC, skip the whole file
#ifdef PETSC

#include "petsc/finclude/petscsnes.h"

subroutine descnt_petsc( NGdof, position, ihybrd )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, one, two, ten

  use numerical, only : sqrtmachprec

  use fileunits, only : ounit

  use inputlist, only : Wmacros, Wdescnt, ext, &
                        Igeometry, & ! only for screen output; 
                        Nvol,                    &
                        Lfindzero, forcetol, c05xmax, c05xtol, c05factor, LreadGF, &
                        Lcheck, dxdesc, ftoldesc, maxitdesc, Lwritedesc, nwritedesc

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

  INTEGER                :: irevcm, mode, Ldfjac, LR, dummy(2)
  REAL                   :: xtol, epsfcn, factor


  REAL                   :: force(0:NGdof)
  REAL, allocatable      :: fjac(:,:), RR(:), work(:,:)
  
  INTEGER                :: ML, MU ! required for only Lc05ndf;
  
  LOGICAL                :: Lexit = .true. ! perhaps this could be made user input;
  LOGICAL                :: LComputeAxis

  INTEGER                :: nprint = 1, nfev, njev, niter

  INTEGER, parameter     :: maxfev = 5000 ! maximum calls per iteration;

  SNES                   :: snes  ! the nonlinear solver
  !PetscViewerAndFormat   :: vf
  Vec                    :: x,r   ! the solution "position" and residual "force"
  Mat                    :: J     ! jacobian matrix "fjac"
  PetscInt               :: its

  external FormFunction, FormMonitor, FormJacobian

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
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  Lhessianallocated = .false.


    ! setup PETSC
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Create nonlinear solver context
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  call SNESCreate(PETSC_COMM_WORLD,snes,ierr);CHKERRA(ierr)
  call SNESSetType(snes,"anderson", ierr);CHKERRA(ierr)
  call SNESSetFromOptions(snes,ierr);CHKERRQ(ierr);

  ! create global variables 
  !call DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,max(NGdof,ncpu),1,1,PETSC_NULL_INTEGER,da,ierr)
  !call DMSetUp(da,ierr);CHKERRA(ierr)

  !call DMGetGlobalVector(da,x,ierr);CHKERRA(ierr)

  call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,NGdof,x,ierr);CHKERRA(ierr)
  call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,NGdof,r,ierr);CHKERRA(ierr)

  call SNESSetFunction(snes,r,FormFunction,dummy,ierr);CHKERRA(ierr)
  call SNESSetJacobian(snes,PETSC_NULL_MAT,PETSC_NULL_MAT,SNESComputeJacobianDefault,0,ierr);CHKERRA(ierr)

  ! setup the diagnostic subroutine
  !call PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT,vf,ierr);CHKERRQ(ierr);
  call SNESMonitorSet(snes,FormMonitor,nwritedesc,PETSC_NULL_FUNCTION,ierr)

  ! fill in with the initial guess
  call put_global_vector(NGdof, position(1:NGdof), x)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  select case( Lfindzero )

  case( 5 ) ! use only function values to find f(x)=0 by pushing the interfaces with f(x) (force-descent)

   !write(*,*) "-------------------- Under construction --------------------"
   call SNESSolve(snes,PETSC_NULL_VEC,x,ierr);CHKERRA(ierr)
   call SNESGetIterationNumber(snes,its,ierr);CHKERRA(ierr)
   !write(*,*) "-------------------- Under construction --------------------"
  
  case default
   
   FATAL( descnt, .true., value of Lfindzero not supported )
   
  end select
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  call get_global_vector(NGdof, x, position(1:NGdof))
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  ! clean up the SNES space
  call VecDestroy(x,ierr);CHKERRA(ierr)
  call VecDestroy(r,ierr);CHKERRA(ierr)
  call SNESDestroy(snes,ierr);CHKERRA(ierr)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
	call MPI_BARRIER( MPI_COMM_WORLD, ierr2)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  RETURN(descnt)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine descnt_petsc

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine FormFunction(snes,x,f,dummy,ierr)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  use constants, only : zero, one, two, ten

  use numerical, only : sqrtmachprec

  use fileunits, only : ounit

  use inputlist, only : Wmacros, Wdescnt

  use cputiming, only : Tdescnt

  use allglobal, only : wrtend, myid, ncpu, cpus

  use petscsnes
  use descnttime

  LOCALS

  REAL, allocatable      :: pos(:), force(:)
  LOGICAL                :: LComputeDerivatives, LComputeAxis

  INTEGER                :: ii, NGdof, istart, iend

  SNES                   :: snes
  Vec                    :: x,f
  INTEGER                :: dummy(2)

  BEGIN(descnt)

  nFcalls = nFcalls + 1

  LComputeDerivatives = .false.
  LComputeAxis = .true.

  ! get the degree of freedom
  call VecGetSize(x,NGdof,ierr);CHKERRA(ierr)
  
  ! allocate the working vector
  SALLOCATE( pos  , (0:NGdof), zero )
  SALLOCATE( force, (0:NGdof), zero )

  call get_global_vector(NGdof, x, pos(1:NGdof))

  WCALL( descnt, dforce, ( NGdof, pos(0:NGdof), force(0:NGdof), LComputeDerivatives, LComputeAxis ) ) ! calculate the force-imbalance;
  ! return the value of force to the solver
  force(1:NGdof) = - force(1:NGdof)
  call put_global_vector(NGdof, force(1:NGdof), f)
  ! deallocate
  DALLOCATE( pos   )
  DALLOCATE( force )
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(descnt)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine FormFunction

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine FormJacobian(snes,x,jac,jac_prec,dummy,ierr)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  use constants, only : zero, one, two, ten

  use numerical, only : sqrtmachprec

  use fileunits, only : ounit

  use inputlist, only : Wmacros, Wdescnt

  use cputiming, only : Tdescnt

  use allglobal, only : wrtend, myid, ncpu, cpus, hessian

  use petscsnes
  use descnttime

  LOCALS

  REAL, allocatable      :: pos(:), force(:)
  LOGICAL                :: LComputeDerivatives

  INTEGER                :: ii, jj, NGdof

  SNES                   :: snes
  Vec                    :: x
  Mat                    :: jac, jac_prec
  INTEGER                   dummy(*)

  BEGIN(descnt)

  nDcalls = nDcalls + 1

  LComputeDerivatives = .true.

  ! get the degree of freedom
  call VecGetSize(x,NGdof,ierr);CHKERRA(ierr)
  
  ! allocate the working vector
  SALLOCATE( pos  , (0:NGdof), zero )
  SALLOCATE( force, (0:NGdof), zero )

  call get_global_vector(NGdof, x, pos(1:NGdof))

  WCALL( descnt, dforce, ( NGdof, pos(0:NGdof), force(0:NGdof), LComputeDerivatives ) ) ! calculate the force-imbalance;
  
  if (myid .eq. 0) then
    do ii = 1, NGdof
      do jj = 1, NGdof
        call MatSetValue(jac,ii-1,jj-1,hessian(ii,jj),INSERT_VALUES,ierr);CHKERRA(ierr)
      enddo
    enddo
  endif

  call MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)
  call MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)

  ! deallocate
  DALLOCATE( pos   )
  DALLOCATE( force )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(descnt)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine FormJacobian

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine FormMonitor(snes,its,norm,dummy)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, one, two, ten

  use numerical, only : sqrtmachprec

  use fileunits, only : ounit

  use inputlist, only : Wmacros, Wdescnt, ext, &
                        Igeometry, & ! only for screen output; 
                        Nvol, epsilon,              &
                        Lfindzero, forcetol, Lwritedesc

  use cputiming, only : Tdescnt

  use allglobal, only : wrtend, myid, ncpu, cpus, &
                        ForceErr, Energy, &
                        mn, im, in, iRbc, iZbs, iRbs, iZbc, Mvol, &
                        BBe, IIo, BBo, IIe, NOTstellsym, &
                        lMMl, lLLl,  sweight

                        
  use sphdf5, only : write_convergence_output

  use descnttime

  LOCALS
  
  SNES       :: snes
  PetscInt   :: its
  PetscReal  :: norm

  INTEGER    :: dummy(2)

  Vec        :: x, r

  INTEGER    :: NGdof, wflag, iflag
  REAL       :: rflag
  REAL,allocatable :: pos(:), force(:) 

  CHARACTER  :: pack
  LOGICAL    :: LComputeDerivatives, LComputeAxis

  BEGIN(descnt)
  
  call SNESGetSolution(snes,x,ierr);CHKERRQ(ierr)

  call VecGetSize(x,NGdof,ierr);CHKERRA(ierr)
  
  ! allocate the working vector
  SALLOCATE( pos  , (0:NGdof), zero )

  call get_global_vector(NGdof, x, pos(1:NGdof))
  pack = 'U' ! unpack geometrical degrees of freedom;
  LComputeDerivatives = .false.
  LComputeAxis = .true.
  WCALL( descnt, packxi, ( NGdof, pos(0:NGdof), Mvol, mn, iRbc(1:mn,0:Mvol), iZbs(1:mn,0:Mvol), iRbs(1:mn,0:Mvol), iZbc(1:mn,0:Mvol), pack, LComputeDerivatives, LComputeAxis ) )
  
  if (myid .eq. 0) then

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

  ! deallocate
  DALLOCATE( pos   )

  RETURN(descnt)
  !do nothing

1000 format("descnt : ",f10.2," : "i9,i3," ; ":"|f|="es12.5" ; ":"time=",f10.2,"s ;":" log"a5"="28f6.2" ...")
1001 format("descnt : ", 10x ," : "9x,3x" ; ":"    "  12x "   ":"     ", 10x ,"  ;":" log"a5"="28f6.2" ...")
1002 format("descnt :            : Energy ", es23.15, "  ;")
1003 format("descnt :            : Energy ", es23.15, "  SEnergy  ", es23.15, "  LEnergy  ", es23.15)
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine FormMonitor

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! helper functions, operating on global vectors
subroutine get_global_vector(n, x_in, x_out)

  use inputlist, only : Wdescnt
  use cputiming, only : Tdescnt
  use allglobal, only : ncpu, myid, cpus
  use fileunits, only : ounit
  use petscsnes

  LOCALS

  Vec      :: x_in
  REAL     :: x_out(n)
  INTEGER  :: n

  PetscScalar,pointer    :: lx(:)

  INTEGER  :: ncount(ncpu), stride(ncpu), nrange(0:ncpu)
  INTEGER  :: ii

  BEGIN(descnt)

  call VecGetArrayReadF90(x_in,lx,ierr);CHKERRQ(ierr)
  ! see how much data is stored locally
  call VecGetOwnershipRanges(x_in,nrange,ierr);CHKERRQ(ierr)
  ncount(1:ncpu) = nrange(1:ncpu) - nrange(0:ncpu-1)
  stride(1:ncpu) = nrange(0:ncpu)
  ! use mpi_allgatherv to collect the vector and give it to all nodes
  call MPI_AllGatherv( lx, ncount(myid+1), MPI_DOUBLE, x_out, ncount, stride, MPI_DOUBLE, MPI_COMM_WORLD, ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(x_in,lx,ierr);CHKERRQ(ierr)

  RETURN(descnt)

end subroutine get_global_vector

subroutine put_global_vector(n, x_in, x_out)

  use inputlist, only : Wdescnt
  use cputiming, only : Tdescnt
  use allglobal, only : ncpu, myid, cpus
  use fileunits, only : ounit
  use petscsnes

  LOCALS

  Vec      :: x_out
  REAL     :: x_in(n)
  INTEGER  :: n

  INTEGER  :: ii, istart, iend

  BEGIN(descnt)

  ! see how much data is stored locally
  call VecGetOwnershipRange(x_out,istart,iend,ierr);CHKERRQ(ierr)
  do ii = istart+1, iend
    call VecSetValues(x_out,1,ii-1,x_in(ii),INSERT_VALUES,ierr);CHKERRA(ierr)
  enddo
  call VecAssemblyBegin(x_out,ierr);CHKERRA(ierr)
  call VecAssemblyEnd(x_out,ierr);CHKERRA(ierr)

  RETURN(descnt)

end subroutine put_global_vector

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
#endif