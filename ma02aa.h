!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (solver/driver) ! Constructs Beltrami field in given volume consistent with constraints.

!latex \briefly{Constructs Beltrami field in given volume consistent with flux, helicity, rotational-transform and/or parallel-current constraints.}

!latex \calledby{\link{dforce}}
!latex \calls{\link{packab}, \link{df00ab}, \link{mp00ac}}

!latex \tableofcontents

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine ma02aa( lvol, NN )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half, one, ten
  
  use numerical, only : vsmall, small

  use fileunits, only : ounit

  use inputlist, only : Wmacros, Wma02aa, &
                        Lconstraint, mu, helicity, &
                        mupftol, mupfits, Lrad, Lcheck

  use cputiming
  
  use allglobal, only : ncpu, myid, cpus, Mvol, mn, im, in, &
                        LBlinear, LBnewton, LBsequad, &
                        dMA, dMB, dMD, solution, &
                        MBpsi, &
                        ImagneticOK, &
                        lBBintegral, lABintegral, &
                        ivol, &
                        dtflux, dpflux, &
                        xoffset, &
                        Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, &
						IndMatrixArray
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in)  :: lvol, NN ! NN is the number of degrees of freedom in the (packed format) vector potential;
  
  
  INTEGER              :: ideriv
  REAL                 :: tol, dpsi(1:2), lastcpu, ind_matrix
  CHARACTER            :: packorunpack
  
  INTEGER              :: Nxdof, Ndof, Ldfjac, iflag, maxfev, mode, LRR, nfev, njev, nprint, ihybrj
  REAL                 :: Xdof(1:2), Fdof(1:2), Ddof(1:2,1:2), oDdof(1:2,1:2)
  REAL                 :: factor, diag(1:2), RR(1:2*(2+1)/2), QTF(1:2), wk(1:2,1:4)
  
  INTEGER              :: irevcm
  
  INTEGER              :: pNN
  
  REAL                 :: xi(0:NN), Fxi(0:NN), xo(0:NN), Mxi(1:NN)
  
  external             :: mp00ac
  
#ifdef DEBUG
  INTEGER              :: ixx, jxx, jfinite ! computing finite-difference derivatives of \iota wrt \mu and \Delta \psi_p;
  REAL                 :: lfdiff, dFdof(-1:1,-1:1,1:2)
  REAL, allocatable    :: dsolution(:,:,:,:)
#endif
  
!required for hybrj1;
  INTEGER              :: ihybrj1, Ldfmuaa, lengthwork
  REAL                 :: DFxi(0:NN,0:NN), work(1:(1+NN)*(1+NN+13)/2), NewtonError
  external             :: df00ab
  
! required for E04UFF;
  INTEGER              :: NLinearConstraints, NNonLinearConstraints, LDA, LDCJ, LDR, iterations, LIWk, LRWk, ie04uff
  INTEGER, allocatable :: Istate(:), NEEDC(:), IWk(:)
  REAL                 :: objectivefunction
  REAL   , allocatable :: LinearConstraintMatrix(:,:), LowerBound(:), UpperBound(:)
  REAL   , allocatable :: constraintfunction(:), constraintgradient(:,:), multipliers(:), objectivegradient(:), RS(:,:), RWk(:)
  CHARACTER            :: optionalparameter*33
  
  BEGIN(ma02aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  FATAL( ma02aa, lvol.lt.1 .or. lvol.gt.Mvol, illegal lvol )
#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  ivol = lvol ! various subroutines (e.g. mp00ac, df00ab) that may be called below require volume identification, but the argument list is fixed by NAG;
  ind_matrix = IndMatrixArray(ivol, 2)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \subsection{seqeuntial quadratic programming}
!latex \begin{enumerate}
!latex \item Only relevant if \internal{LBsequad=T}. See \inputvar{LBeltrami} for details.
!latex \item Documentation on the implementation of \nag{www.nag.co.uk/numeric/FL/manual19/pdf/E04/e04uff_fl19.pdf}{E04UFF} is under construction.
!latex \end{enumerate}

  if( LBsequad ) then ! sequential quadratic programming (SQP); construct minimum energy with constrained helicity;
   lastcpu = GETTIME
   
   NLinearConstraints = 0 ! no linear constraints;
   
   NNonLinearConstraints = 1 ! single non-linear constraint = conserved helicity;
   
   LDA = max(1,NLinearConstraints)
   
   LDCJ = max(1,NNonLinearConstraints)
   
   LDR = NN
   
   SALLOCATE( LinearConstraintMatrix, (1:LDA,1:1), zero ) ! linear constraint matrix;
   
   SALLOCATE( LowerBound, (1:NN+NLinearConstraints+NNonLinearConstraints), zero ) ! lower bounds on variables, linear constraints and non-linear constraints;
   SALLOCATE( UpperBound, (1:NN+NLinearConstraints+NNonLinearConstraints), zero ) ! upper bounds on variables, linear constraints and non-linear constraints;

   LowerBound(                       1 : NN                                          ) = -1.0E+21       !   variable constraints; no constraint;
   UpperBound(                       1 : NN                                          ) = +1.0E+21       !
   LowerBound( NN+                   1 : NN+NLinearConstraints                       ) = -1.0E+21       !     linear constraints; no constraint;
   UpperBound( NN+                   1 : NN+NLinearConstraints                       ) = +1.0E+21       !
   LowerBound( NN+NLinearConstraints+1 : NN+NLinearConstraints+NNonLinearConstraints ) = helicity(lvol) ! non-linear constraints; enforce helicity constraint;
   UpperBound( NN+NLinearConstraints+1 : NN+NLinearConstraints+NNonLinearConstraints ) = helicity(lvol) !
   
   iterations = 0 ! iteration counter;
   
   SALLOCATE( Istate, (1:NN+NLinearConstraints+NNonLinearConstraints), 0 )
   
   SALLOCATE( constraintfunction, (1:NNonLinearConstraints), zero ) ! constraint functions;
   
   SALLOCATE( constraintgradient, (1:LDCJ,1:NN), zero ) ! derivatives of constraint functions;
   
   SALLOCATE( multipliers, (1:NN+NLinearConstraints+NNonLinearConstraints), zero ) ! Lagrange multipliers ?;
   
   objectivefunction = zero ! objective function;
   
   SALLOCATE( objectivegradient, (1:NN), zero ) ! derivatives of objective function;
   
   SALLOCATE( RS, (1:LDR,1:NN), zero )
   ideriv = 0 ; dpsi(1:2) = (/ dtflux(lvol), dpflux(lvol) /) ! these are also used below;
   
   packorunpack = 'P'

   CALL( ma02aa, packab, ( packorunpack, lvol, NN, xi(1:NN), ideriv ) )
   
   SALLOCATE( NEEDC, (1:NNonLinearConstraints), 0 )
   
   LIWk = 3*NN + NLinearConstraints + 2*NNonLinearConstraints ! workspace;
   SALLOCATE( IWk, (1:LIWk), 0 )       ! workspace;
   
   LRWk = 2*NN**2 + NN * NLinearConstraints + 2 * NN * NNonLinearConstraints + 21 * NN + 11 * NLinearConstraints + 22 * NNonLinearConstraints + 1 ! workspace;
   SALLOCATE( RWk, (1:LRWk), zero )                                                              ! workspace;
   
   irevcm = 0 ; ie04uff = 1 ! reverse communication loop control; ifail error flag;
   
! supply optional parameters to E04UFF; NAG calls commented out (this part of the code so far not used); 17 Nov 17
   
!   call E04UEF('Nolist')               ! turn of screen output;
!   call E04UEF('Print Level = 0')      ! turn of screen output;
!   call E04UEF('Derivative Level = 3') ! assume all derivatives are provided by user;
!   call E04UEF('Verify Level = -1')    ! do not verify derivatives using finite-differences; default is Verify Level = 0, which does verify gradients;
!   write(optionalparameter,'("Major Iteration Limit = "i9)') 2**2 * max( 50, 3 * ( NN + NLinearConstraints ) + 10 * NNonLinearConstraints )
!   call E04UEF(optionalparameter)
   
! pre-calculate some matrix vector products;
   
   MBpsi(ind_matrix)%arr(1:NN) =                         matmul( dMB(ind_matrix)%mat(1:NN,1: 2), dpsi(1:2) )
!  MEpsi(1:NN) = zero !                  matmul( dME(1:NN,1: 2), dpsi(1:2) )
   
!  psiMCpsi    = zero ! half * sum( dpsi(1:2) * matmul( dMC(1: 2,1: 2), dpsi(1:2) ) )
!  psiMFpsi    = zero ! half * sum( dpsi(1:2) * matmul( dMF(1: 2,1: 2), dpsi(1:2) ) )
   
   
!   do ! reverse communication loop; NAG calls commented out (this part of the code so far not used); 17 Nov 17
!    
!    
!    call E04UFF( irevcm, &
!                 NN, NLinearConstraints, NNonLinearConstraints, LDA, LDCJ, LDR, &
!                 LinearConstraintMatrix(1:LDA,1:1), &
!                 LowerBound(1:NN+NLinearConstraints+NNonLinearConstraints), UpperBound(1:NN+NLinearConstraints+NNonLinearConstraints), &
!                 iterations, Istate(1:NN+NLinearConstraints+NNonLinearConstraints), &
!                constraintfunction(1:NNonLinearConstraints), constraintgradient(1:LDCJ,1:NN), &
!                 multipliers(1:NN+NLinearConstraints+NNonLinearConstraints), &
!                 objectivefunction, objectivegradient(1:NN), &
!                 RS(1:LDR,1:NN), &
!                 xi(1:NN), &
!                 NEEDC(1:NNonLinearConstraints), IWk(1:LIWk), LIWk, RWk(1:LRWk), LRWk, ie04uff )
!
!    if( irevcm.eq.1 .or. irevcm.eq.2 .or. irevcm.eq.3 ) Mxi(1:NN) = matmul( dMA(ind_matrix)%mat(1:NN,1:NN), xi(1:NN) ) ! calculate objective  functional and/or gradient;
!    if( irevcm.eq.4 .or. irevcm.eq.5 .or. irevcm.eq.6 ) Mxi(1:NN) = matmul( dMD(ind_matrix)%mat(1:NN,1:NN), xi(1:NN) ) ! calculate constraint functional and/or gradient;
!    
!    if( irevcm.eq.1 .or. irevcm.eq.3 ) objectivefunction       = half * sum( xi(1:NN) * Mxi(1:NN) ) + sum( xi(1:NN) * MBpsi(ind_matrix)%arr(1:NN) ) + psiMCpsi
!    if( irevcm.eq.2 .or. irevcm.eq.3 ) objectivegradient(1:NN) =                        Mxi(1:NN)   +                 MBpsi(ind_matrix)%arr(1:NN)
!    
!   if( irevcm.eq.4 .or. irevcm.eq.6 .and. NEEDC(1).gt.0 ) then
!    constraintfunction(1     ) = half * sum( xi(1:NN) * Mxi(1:NN) ) + sum( xi(1:NN) * MEpsi(1:NN) ) + psiMFpsi
!   endif
!    
!    if( irevcm.eq.5 .or. irevcm.eq.6 .and. NEEDC(1).gt.0 ) then
!     constraintgradient(1,1:NN) =                        Mxi(1:NN)   +                 MEpsi(1:NN) 
!    endif
!    
!    if( irevcm.eq.0 ) then ! final exit;
!     
!     cput = GETTIME
!     
!     select case(ie04uff)
!     case( :-1 )  
!      write(ounit,1010) cput-cpus, myid, lvol, ie04uff, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "user enforced termination ;     "
!     case(   0 )  
!     if( Wma02aa ) write(ounit,1010) cput-cpus, myid, lvol, ie04uff, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "success ;                       "
!     case(   1 )  
!      if( Wma02aa ) write(ounit,1010) cput-cpus, myid, lvol, ie04uff, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "not converged;                  "
!     case(   2 )  
!      write(ounit,1010) cput-cpus, myid, lvol, ie04uff, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "infeasible (linear) ;           "
!     case(   3 ) 
!      write(ounit,1010) cput-cpus, myid, lvol, ie04uff, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "infeasible (nonlinear) ;        "
!     case(   4 )   
!      write(ounit,1010) cput-cpus, myid, lvol, ie04uff, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "major iteration limit reached ; "
!     case(   6 )  
!     if( Wma02aa ) write(ounit,1010) cput-cpus, myid, lvol, ie04uff, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "could not be improved ;         "
!     case(   7 )  
!     write(ounit,1010) cput-cpus, myid, lvol, ie04uff, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "derivatives appear incorrect ;  "
!    case(   9 )  
!     write(ounit,1010) cput-cpus, myid, lvol, ie04uff, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "input error ;                   "
!     case default
!     FATAL( ma02aa, .true., illegal ifail returned by E04UFF )
!    end select
     
!    if( irevcm.eq.0 .or. irevcm.eq.1 .or. irevcm.eq.6 ) ImagneticOK(lvol) = .true. ! set error flag; used elsewhere;
!    
!     mu(lvol) = multipliers( NN + NLinearConstraints + NNonLinearConstraints ) ! helicity multiplier, or so it seems: NAG document is not completely clear;
!     
!     exit ! sequential quadratic programming method of constructing Beltrami field is finished;
!     
!    endif ! end of if( irevcm.eq.0 ) then;
!    
!   enddo ! end of do ! reverse communication loop;
   
   
   DALLOCATE(RWk)
   DALLOCATE(IWk)
   DALLOCATE(NEEDC)
   DALLOCATE(RS)
   DALLOCATE(objectivegradient)
   DALLOCATE(multipliers)
   DALLOCATE(constraintgradient)
   DALLOCATE(constraintfunction)
   DALLOCATE(Istate)
   DALLOCATE(LowerBound)
   DALLOCATE(UpperBound)
   DALLOCATE(LinearConstraintMatrix)
   
   
   packorunpack = 'U'
   CALL( ma02aa, packab ( packorunpack, lvol, NN, xi(1:NN), ideriv ) )
   
   lBBintegral(lvol) = half * sum( xi(1:NN) * matmul( dMA(ind_matrix)%mat(1:NN,1:NN), xi(1:NN) ) ) + sum( xi(1:NN) * MBpsi(ind_matrix)%arr(1:NN) ) ! + psiMCpsi
   lABintegral(lvol) = half * sum( xi(1:NN) * matmul( dMD(ind_matrix)%mat(1:NN,1:NN), xi(1:NN) ) ) ! + sum( xi(1:NN) * MEpsi(1:NN) ) ! + psiMFpsi
   
   solution(ivol)%mat(1:NN,0) = xi(1:NN)
   
   
  endif ! end of if( LBsequad ) then;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \subsection{Newton method}
!latex \begin{enumerate}
!latex \item Only relevant if \internal{LBnewton=T}. See \inputvar{LBeltrami} for details.
!latex \end{enumerate}

  if( LBnewton ) then
   
   lastcpu = GETTIME
   
   xi(0) = mu(lvol) ! initialize; helicity multiplier is treated as an independent degree-of-freedom;
   
   ideriv = 0 ; dpsi(1:2) = (/ dtflux(lvol), dpflux(lvol) /) ! these are also used below;
   
   packorunpack = 'P'
   CALL( ma02aa, packab ( packorunpack, lvol, NN, xi(1:NN), dpsi(1:2), ideriv ) )
   
   pNN = NN + 1 ; Ldfmuaa = pNN ; tol = mupftol ; lengthwork = pNN * ( pNN+13 ) / 2
   
! pre-calculate some matrix vector products; these are used in df00ab;
   
   MBpsi(ind_matrix)%arr(1:NN) =                         matmul( dMB(ind_matrix)%mat(1:NN,1: 2), dpsi(1:2) )
!  MEpsi(1:NN) = zero !                  matmul( dME(1:NN,1: 2), dpsi(1:2) )
!  psiMCpsi    = zero ! half * sum( dpsi(1:2) * matmul( dMC(1: 2,1: 2), dpsi(1:2) ) )
!  psiMFpsi    = zero ! half * sum( dpsi(1:2) * matmul( dMF(1: 2,1: 2), dpsi(1:2) ) )
   
   call hybrj1( df00ab, pNN, xi(0:NN), Fxi(0:NN), DFxi(0:NN,0:NN), Ldfmuaa, tol, ihybrj1, work(1:lengthwork), lengthwork )

   NewtonError = maxval( abs( Fxi(0:NN) ) )
   
   mu(lvol) = xi(0)
   
   packorunpack = 'U' ; ideriv = 0

   CALL( ma02aa, packab( packorunpack, lvol, NN, xi(1:NN), ideriv ) )
   
   cput = GETTIME

   select case( ihybrj1 )
   case( :-1 )    
    ;             write(ounit,1020) cput-cpus, myid, lvol, ihybrj1, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, NewtonError, "user termination ;"
   case( 1 )    
    if( Wma02aa ) write(ounit,1020) cput-cpus, myid, lvol, ihybrj1, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, NewtonError, "success ;         "
   case( 0 )    
    ;             write(ounit,1020) cput-cpus, myid, lvol, ihybrj1, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, NewtonError, "input error ;     "
   case( 2 )    
    ;             write(ounit,1020) cput-cpus, myid, lvol, ihybrj1, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, NewtonError, "max. evaluations ;"
   case( 3 )    
    ;             write(ounit,1020) cput-cpus, myid, lvol, ihybrj1, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, NewtonError, "xtol too small ;  "
   case( 4 )    
    ;             write(ounit,1020) cput-cpus, myid, lvol, ihybrj1, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, NewtonError, "bad progress ;    "
   case default 
    FATAL( ma02aa, .true., illegal ifail returned by hybrj1 )
   end select
   
#ifdef DEBUG
   xo(1:NN) = xi(1:NN) ! save original for comparison;
   packorunpack = 'P' ; ideriv = 0
   CALL( ma02aa, packab( packorunpack, lvol, NN, xi(1:NN), ideriv ) )
   FATAL( ma02aa, sum(abs(xi(1:NN)-xo(1:NN)))/NN.gt.vsmall, un/packing routine is incorrect )
#endif
   
!if( NewtonError.lt.mupftol ) then
   ImagneticOK(lvol) = .true. 
!endif
   
   lBBintegral(lvol) = half * sum( xi(1:NN) * matmul( dMA(ind_matrix)%mat(1:NN,1:NN), xi(1:NN) ) ) + sum( xi(1:NN) * MBpsi(ind_matrix)%arr(1:NN) ) ! + psiMCpsi
   lABintegral(lvol) = half * sum( xi(1:NN) * matmul( dMD(ind_matrix)%mat(1:NN,1:NN), xi(1:NN) ) ) ! + sum( xi(1:NN) * MEpsi(1:NN) ) ! + psiMFpsi
   
   solution(ivol)%mat(1:NN,0) = xi(1:NN)
   
  endif ! end of if( LBnewton ) then
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \subsection{``linear'' method}

!latex \begin{enumerate}
!latex \item Only relevant if \internal{LBlinear=T}. See \inputvar{LBeltrami} for details.
!latex \item The quantity $\mu$ is {\em not} not treated as a ``magnetic'' degree-of-freedom
!latex       equivalent to in the degrees-of-freedom in the magnetic vector potential
!latex       (as it strictly should be, because it is a Lagrange multiplier introduced to enforce the helicity constraint).
!latex \item In this case, the Beltrami equation, $\nabla \times {\bf B} = \mu {\bf B}$, is {\em linear} in the magnetic degrees-of-freedom.
!latex \item The algorithm proceeds as follows:

!latex \subsubsection{plasma volumes} 

!latex \begin{enumerate}

!latex \item In addition to the enclosed toroidal flux, $\Delta \psi_t$, which is held constant in the plasma volumes,
!latex       the Beltrami field in a given volume is assumed to be parameterized by $\mu$ and $\Delta \psi_p$.
!latex       (Note that $\Delta \psi_p$ is not defined in a torus.)
!latex \item These are ``packed'' into an array, e.g. $\boldmu \equiv (\mu, \Delta\psi_p)^T$, so that standard library routines ,
!latex       e.g. \nag{www.nag.co.uk/numeric/FL/manual19/pdf/C05/c05pcf_fl19.pdf}{C05PCF},
!latex       can be used to (iteratively) find the appropriately-constrained Beltrami solution, i.e. ${\bf f}(\boldmu)=0$.
!latex \item The function ${\bf f}(\boldmu)$, which is computed by \link{mp00ac}, is defined by the input parameter \inputvar{Lconstraint}:
!latex \bi
!latex \item[i.] If \inputvar{Lconstraint = -1, 0}, then $\boldmu$ is {\em not} varied and \internal{Nxdof=0}.
!latex \item[ii.] If \inputvar{Lconstraint = 1},    then $\boldmu$ is           varied to satisfy the transform constraints;
!latex            and \internal{Nxdof=1} in the simple torus and \internal{Nxdof=2} in the annular regions.
!latex            (Note that in the ``simple-torus'' region, the enclosed poloidal flux $\Delta\psi_p$ is not well-defined, 
!latex            and only $\mu=\boldmu_1$ is varied in order to satisfy the transform constraint on the ``outer'' interface of that volume.)
!latex \item[iii.] If \inputvar{Lconstraint = 2}, then $\mu=\boldmu_1$ is           varied in order to satisfy the helicity constraint, 
!latex             and $\Delta\psi_p=\boldmu_2$ is {\em not} varied, and \internal{Nxdof=1}.
!latex             {\bf (under re-construction)}
!latex \ei

!latex \end{enumerate} 

!latex \subsubsection{vacuum volume } 

!latex \begin{enumerate}

!latex \item In the vacuum, $\mu=0$, and the enclosed fluxes, $\Delta \psi_t$ and $\Delta \psi_p$, are considered to parameterize the family of solutions.
!latex       (These quantities may not be well-defined if ${\bf B}\cdot{\bf n}\ne 0$ on the computational boundary.)
!latex \item These are ``packed'' into an array, $\boldmu \equiv (\Delta\psi_t,\Delta\psi_p)^T$, so that, as above, standard routines can be used
!latex       to iteratively find the appropriately constrained solution, i.e. ${\bf f}(\boldmu)=0$.
!latex \item The function ${\bf f}(\boldmu)$, which is computed by \link{mp00ac}, is defined by the input parameter \inputvar{Lconstraint}:
!latex \bi
!latex \item[i.] If \inputvar{Lconstraint = -1}, then $\boldmu$ is {\em not} varied and \internal{Nxdof=0}.
!latex \item[ii.] If \inputvar{Lconstraint = 0,2}, then $\boldmu$ is varied to satisfy the enclosed current constraints, and \internal{Nxdof=2}.
!latex \item[iii.] If \inputvar{Lconstraint = 1}, then $\boldmu$ is varied to satisfy 
!latex             the constraint on the transform on the inner boundary $\equiv$ plasma boundary and the ``linking'' current, and \internal{Nxdof=2}. 
!latex \ei

!latex \end{enumerate} 

!latex \item The Beltrami fields, and the rotational-transform and helicity etc. as required to determine the function ${\bf f}(\boldmu)$
!latex       are calculated in \link{mp00ac}.
!latex \item This routine, \link{mp00ac}, is called iteratively if \inputvar{Nxdof > 1} via
!latex       \nag{www.nag.co.uk/numeric/FL/manual19/pdf/C05/c05pcf_fl19.pdf}{C05PCF} to
!latex       determine the appropriately constrained Beltrami field, ${\bf B}_{\boldmu}$, so that ${\bf f}(\boldmu)=0$.
!latex \item The input variables \inputvar{mupftol} and \inputvar{mupfits} control the required accuracy and maximum number of iterations.
!latex \item If \inputvar{Nxdof = 1}, then \link{mp00ac} is called only once to provide the Beltrami fields with the given value of $\boldmu$.

!latex \end{enumerate}
  

  if( LBlinear ) then ! assume Beltrami field is parameterized by helicity multiplier (and poloidal flux);
   
   lastcpu = GETTIME
   
   if( Lplasmaregion ) then
    
    Xdof(1:2) = xoffset + (/     mu(lvol), dpflux(lvol) /) ! initial guess for degrees of freedom; offset from zero so that relative error is small;
    
    select case( Lconstraint )
    case( -1 )    ;                                   ; Nxdof = 0 ! multiplier & poloidal flux NOT varied                               ;
    case(  0 )    ;                                   ; Nxdof = 0 ! multiplier & poloidal flux NOT varied                               ;
    case(  1 )    ; if( Lcoordinatesingularity ) then ; Nxdof = 1 ! multiplier                 IS  varied to match       outer transform;
     ;              else                              ; Nxdof = 2 ! multiplier & poloidal flux ARE varied to match inner/outer transform;
     ;              endif                                         
    case(  2 )    ;                                     Nxdof = 1 ! multiplier                 IS  varied to match             helicity ;
    case(  3 )    ; if( Lcoordinatesingularity ) then ; Nxdof = 0 ! multiplier & poloidal flux NOT varied                               ;
     ;              else                              ; Nxdof = 0 ! Global constraint, no dof locally
     ;              endif
    end select
    
   else ! Lvacuumregion ;

    Xdof(1:2) = xoffset + (/ dtflux(lvol), dpflux(lvol) /) ! initial guess for degrees of freedom; offset from zero so that relative error is small;
    
    select case( Lconstraint )
    case( -1 )    ;                                   ; Nxdof = 0 ! poloidal   & toroidal flux NOT varied to match linking current and plasma current;
    case(  0 )    ;                                   ; Nxdof = 2 ! poloidal   & toroidal flux ARE varied to match linking current and plasma current;
    case(  1 )    ;                                   ; Nxdof = 2 ! poloidal   & toroidal flux ARE varied to match linking current and transform-constraint;
    case(  2 )    ;                                   ; Nxdof = 2 ! poloidal   & toroidal flux ARE varied to match linking current and plasma current;
    case(  3 )    ;                                   ; Nxdof = 0 ! Global constraint
    end select

   endif ! end of if( Lplasmaregion) ;
   
   select case( Nxdof )
    
   case( 0   ) ! need only call mp00ac once, to calculate Beltrami field for given helicity multiplier and enclosed fluxes;
    
    iflag = 1 ; Ndof = 1     ; Ldfjac = Ndof ; nfev = 1 ; njev = 0 ; ihybrj = 1;  ! provide dummy values for consistency;
    
    WCALL( ma02aa, mp00ac, ( Ndof, Xdof(1:Ndof), Fdof(1:Ndof), Ddof(1:Ldfjac,1:Ndof), Ldfjac, iflag ) )
    
    helicity(lvol) = lABintegral(lvol) ! this was computed in mp00ac;
    
   case( 1:2 ) ! will iteratively call mp00ac, to calculate Beltrami field that satisfies constraints;
    
    ;         ; Ndof = Nxdof ; Ldfjac = Ndof ; nfev = 0 ; njev = 0 ; ihybrj = 0;

    tol = mupftol ; LRR = Ndof * ( Ndof+1 ) / 2 ; mode = 0 ; diag(1:2) = zero ; factor = one ; maxfev = mupfits ; nprint = 0
    
    FATAL( ma02aa, Ndof.gt.2, illegal )

    WCALL( ma02aa, hybrj2, ( mp00ac, Ndof, Xdof(1:Ndof), Fdof(1:Ndof), Ddof(1:Ldfjac,1:Ndof), Ldfjac, tol, &
                             maxfev, diag(1:Ndof), mode, factor, nprint, ihybrj, nfev, njev, RR(1:LRR), LRR, QTF(1:Ndof), &
			     WK(1:Ndof,1), WK(1:Ndof,2), WK(1:Ndof,3), WK(1:Ndof,4) ) )

    if( Lplasmaregion ) then
     
     select case( ihybrj )
     case( 0: ) ;     mu(lvol) = Xdof(1)      - xoffset
      ;         ; dpflux(lvol) = Xdof(2)      - xoffset
     case( :-1) ;      Xdof(1) = mu(lvol)     + xoffset ! mu    and dpflux have been updated in mp00ac; early termination;
      ;         ;      Xdof(2) = dpflux(lvol) + xoffset ! mu    and dpflux have been updated in mp00ac; early termination;
     end select
     
    else ! Lvacuumregion;
     
     select case( ihybrj )
     case( 0: ) ; dtflux(lvol) = Xdof(1)      - xoffset
      ;         ; dpflux(lvol) = Xdof(2)      - xoffset
     case( :-1) ; Xdof(1)      = dtflux(lvol) + xoffset ! dtflux and dpflux have been updated in mp00ac; early termination;
      ;         ; Xdof(2)      = dpflux(lvol) + xoffset ! dtflux and dpflux have been updated in mp00ac; early termination;
     end select
     
    endif ! end of if( Lplasmaregion ) ;
        
    helicity(lvol) = lABintegral(lvol) ! this was computed in mp00ac;
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    if( Lconstraint.eq.1 .or. ( Lvacuumregion .and. Lconstraint.eq.0 ) ) then
     
     iflag = 2 ; Ldfjac = Ndof ! call mp00ac: tr00ab/curent to ensure the derivatives of B, transform, currents, wrt mu/dtflux & dpflux are calculated;

     WCALL( ma02aa, mp00ac, ( Ndof, Xdof(1:Ndof), Fdof(1:Ndof), Ddof(1:Ldfjac,1:Ndof), Ldfjac, iflag ) )

    endif ! end of if( Lconstraint.eq.1 .or. ( Lvacuumregion .and. Lconstraint.eq.0 ) ) ;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
   end select ! end of select case( Nxdof ) ;
   

   cput = GETTIME

   select case(ihybrj) ! this screen output may not be correct for Lvacuumregion;
   case(    1   )  
    if( Wma02aa ) write(ounit,1040) cput-cpus, myid, lvol, ihybrj, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "success         ", Fdof(1:Ndof)
   case(   -2   )  
    if( Wma02aa ) write(ounit,1040) cput-cpus, myid, lvol, ihybrj, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "|F| < mupftol   ", Fdof(1:Ndof)
   case(   -1   )  
    ;             write(ounit,1040) cput-cpus, myid, lvol, ihybrj, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "Beltrami fail   ", Fdof(1:Ndof)
   case(    0   )  
    ;             write(ounit,1040) cput-cpus, myid, lvol, ihybrj, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "input error     ", Fdof(1:Ndof)
   case(    2   )  
    ;             write(ounit,1040) cput-cpus, myid, lvol, ihybrj, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "consider restart", Fdof(1:Ndof)
   case(    3   )  
    ;             write(ounit,1040) cput-cpus, myid, lvol, ihybrj, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "xtol too small  ", Fdof(1:Ndof)
   case(    4:5 )  
    ;             write(ounit,1040) cput-cpus, myid, lvol, ihybrj, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "bad progress    ", Fdof(1:Ndof)
   case default    
    ;             write(ounit,1040) cput-cpus, myid, lvol, ihybrj, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "illegal ifail   ", Fdof(1:Ndof)
    FATAL( ma02aa, .true., illegal ifail returned by hybrj )
   end select   

  endif ! end of if( LBlinear ) then;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
 
!latex \subsection{debugging: finite-difference confirmation of the derivatives of the rotational-transform}

!latex \begin{enumerate}

!latex \item Note that the rotational-transform (if required) is calculated by \link{tr00ab}, which is called by \link{mp00ac}.
!latex \item If \inputvar{Lconstraint=1}, then \link{mp00ac} will ask \link{tr00ab} to compute the derivatives of the transform 
!latex       with respect to variations in the helicity-multiplier, $\mu$, and the enclosed poloidal-flux, $\Delta\psi_p$, so that
!latex       \nag{www.nag.co.uk/numeric/FL/manual19/pdf/C05/c05pcf_fl19.pdf}{C05PCF} may more efficiently find the solution.
!latex \item The required derivatives are
!latex       \be \frac{\partial \iotabar}{\partial \mu}\\
!latex           \frac{\partial \iotabar}{\partial \Delta \psi_p}
!latex       \ee
!latex       to improve the efficiency of the iterative search.
!latex       A finite difference estimate of these derivatives is available; need \type{DEBUG}, \inputvar{Lcheck=2} and \inputvar{Lconstraint=1}.

!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  
  if( Lcheck.eq.2 ) then ! perform finite difference check on derivative calculation;
   
   if( Lconstraint.eq.1 .or. Lvacuumregion ) then ! only in this case are the derivatives calculated;
    
    SALLOCATE( dsolution, (1:NN,0:2,-1:1,-1:1), zero ) ! packed vector potential;
    
    if( Lplasmaregion ) then
     Xdof(1:2) = xoffset + (/     mu(lvol), dpflux(lvol) /) ! initial guess for degrees of freedom; offset from zero so that relative error is small;      
     if( Lcoordinatesingularity ) then ; Ndof = 1 ! multiplier                  IS  varied to match       transform      ;
     else                              ; Ndof = 2 ! multiplier &  poloidal flux ARE varied to match inner/outer transform;
     endif
    else ! Lvacuumregion;
     Xdof(1:2) = xoffset + (/ dtflux(lvol), dpflux(lvol) /) ! initial guess for degrees of freedom; offset from zero so that relative error is small;
     ;                                 ; Ndof = 2
    endif ! end of if( Lplasmaregion) ;
    
    Ldfjac = Ndof ; dFdof(-1:1,-1:1,1:2) = zero
    
    ixx = 0 ; jxx = 0
    
    iflag = 2 ! iflag controls derivative calculation in mp00ac; analytic derivatives of rotational-transform are required;
    
    CALL( ma02aa, mp00ac( Ndof, Xdof(1:Ndof), dFdof(ixx,jxx,1:Ndof), oDdof(1:Ldfjac,1:Ndof), Ldfjac, iflag ) ) ! compute "exact" derivatives;
    
    dsolution(1:NN,1,0,0) = solution(ivol)%mat(1:NN,1) ! packed vector potential; derivative wrt mu    ;
    dsolution(1:NN,2,0,0) = solution(ivol)%mat(1:NN,2) ! packed vector potential; derivative wrt dpflux;
    
    jfinite = 0
    cput = GETTIME 
    write(ounit,2000) cput-cpus, myid, lvol, jfinite, "derivative", oDdof(1:Ldfjac,1:Ndof)
    
    do jfinite = -4,-2,+1 ; lfdiff = ten**jfinite
     
     do ixx = -1, +1 ! centered differences;
      do jxx = -1, +1 ! centered differences;
       
       if( Lcoordinatesingularity .and. jxx.ne.0 ) cycle ! Beltrami field in volume with coordinate singularity does not depend on enclosed poloidal flux;
       
       if( abs(ixx)+abs(jxx).ne.1 ) cycle ! centered finite differences;
       
       if( Lplasmaregion ) then
        Xdof(1:2) = xoffset + (/     mu(lvol), dpflux(lvol) /) + (/ ixx, jxx /) * lfdiff * half
       else ! Lvacuumregion;
        Xdof(1:2) = xoffset + (/ dtflux(lvol), dpflux(lvol) /) + (/ ixx, jxx /) * lfdiff * half
       endif ! end of if( Lplasmaregion) ;

       iflag = 1 ! analytic derivatives of rotational-transform are not required;
       
       CALL( ma02aa, mp00ac( Ndof, Xdof(1:Ndof), dFdof(ixx,jxx,1:Ndof), Ddof(1:Ldfjac,1:Ndof), Ldfjac, iflag ) ) ! compute function values only;
       
       dsolution(1:NN,0,ixx,jxx) = solution(ivol)%mat(1:NN,0)
       
      enddo ! end of do jxx;
     enddo ! end of do ixx;
     
     ;               Ddof(1:Ndof,1) = ( dFdof( 1, 0, 1:Ndof) - dFdof(-1, 0, 1:Ndof) ) / lfdiff ! derivative wrt helicity multiplier   ;
     if( Ndof.eq.2 ) Ddof(1:Ndof,2) = ( dFdof( 0, 1, 1:Ndof) - dFdof( 0,-1, 1:Ndof) ) / lfdiff ! derivative wrt enclosed poloidal flux;
     
     cput = GETTIME 
    !write(ounit,2000) cput-cpus, myid, lvol, jfinite, " error    ", Ddof(1:Ldfjac,1:Ndof) - oDdof(1:Ldfjac,1:Ndof)
     write(ounit,2000) cput-cpus, myid, lvol, jfinite, " estimate ", Ddof(1:Ldfjac,1:Ndof)
     
    enddo ! end of do jfinite;
    
    cput = GETTIME 
    write(ounit,2000) cput-cpus
    
    DALLOCATE(dsolution)
    
2000 format("ma02aa : ":,f10.2," :":" myid=",i3," : lvol=",i3," ; jj=",i3," ; "a10" : DF=["es23.15" ,"es23.15" ,"es23.15" ,"es23.15" ] ;")
    
   endif ! end of if( Lconstraint.eq.1 .or. Lvacuumregion ) ;
   
  endif ! end of if( Lcheck.eq.2 ) ;

#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(ma02aa)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
1010 format("ma02aa : ",f10.2," : myid=",i3," ; lvol=",i3," ; SQP    : ie04uff=",i3," hel="es12.4" mu="es12.4" dpflux="es12.4" time="f9.1" ":,a36)
1020 format("ma02aa : ",f10.2," : myid=",i3," ; lvol=",i3," ; Newton : ihybrj1=",i3," hel="es12.4" mu="es12.4" dpflux="es12.4" time="f9.1" ; "&
  "error="es7.0" ; ":,a18)
1040 format("ma02aa : ",f10.2," : myid=",i3," ; lvol=",i3," ; Linear : ihybrj =",i3," hel="es12.4" mu="es12.4" dpflux="es12.4" time="f9.1" ; "&
  :,a16" ; F="2es08.0)
!050 format("ma02aa : ",f10.2," : myid=",i3," ; lvol=",i3," ; Linear : ihybrj =",i3,"     "  12x " I ="es12.4"        "  12x " time="f9.1" ; "&
! :,a16" ; F="2es08.0)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine ma02aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
