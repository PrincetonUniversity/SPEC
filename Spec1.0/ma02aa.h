!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Constructs Beltrami field in given volume consistent with constraints.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine ma02aa( lvol, NN )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half, one, two, ten, goldenmean

  use numerical, only : sqrtmachprec, vsmall, small

  use fileunits, only : ounit

  use inputlist, only : Wmacros, Wma02aa, &
                        Lconstraint, mu, helicity, curtor, curpol, &
                        mupftol, mupfits, helicity, Lrad, Lcheck

  use cputiming

  use allglobal, only : ncpu, myid, cpus, Mvol, mn, im, in, &
                        LBeltramiLinear, LBeltramiNewton, LBeltramiSeQuad, &
                        dMA, dMB, dMC, dMD, dME, dMF, solution, &
                        MBpsi, MEpsi, psiMCpsi, psiMFpsi, &
                        ImagneticOK, &
                        lBBintegral, lABintegral, &
                        ivol, &
                        dtflux, dpflux, &
                        xoffset, &
                        Lcoordinatesingularity, Lplasmaregion, Lvacuumregion

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS

  INTEGER, intent(in)  :: lvol, NN ! NN is the number of degrees of freedom in the (packed format) vector potential;


  INTEGER              :: ideriv
  REAL                 :: tol, dpsi(1:2), lastcpu
  CHARACTER            :: packorunpack

  INTEGER              :: Nmudp, Ldfjac, iflag, maxfev, mode, LRR, nfev, njev, nprint, ic05pcf
  REAL                 :: Xmudp(1:2), Fmudp(1:2), Dmudp(1:2,1:2), oDmudp(1:2,1:2)
  REAL                 :: factor, diag(1:2), RR(1:2*(2+1)/2), QTF(1:2), wk(1:2,1:4)
  
  INTEGER              :: irevcm

  INTEGER              :: pNN

  REAL                 :: xi(0:NN), Fxi(0:NN), xo(0:NN), Mxi(1:NN)

  external             :: mp00ac
  
#ifdef DEBUG
  INTEGER              :: imu, idp, jfinite, icurtor ! computing finite-difference derivatives of \iota wrt \mu and \Delta \psi_p; 04 Dec 14;
  REAL                 :: lfdiff, dFmudp(-1:1,-1:1,1:2)
  REAL, allocatable    :: dsolution(:,:,:,:)
#endif

!required for C05PBF;
  INTEGER              :: ic05pbf, Ldfmuaa, lengthwork
  REAL                 :: DFxi(0:NN,0:NN), work(1:(1+NN)*(1+NN+13)/2), NewtonError
  external             :: df00aa

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
  FATALMESS(ma02aa, lvol.lt.1 .or. lvol.gt.Mvol, illegal lvol) ! 28 Jan 13;
#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

 !cput = GETTIME
 !if( myid.eq.0 ) write(ounit,'("ma02aa : ",f10.2," : myid=",i3," ; lvol="i4" ;")') cput-cpus, myid, lvol

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  ivol = lvol ! various subroutines (e.g. mp00ac, df00aa) that may be called below require volume identification, but the argument list is fixed by NAG;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Lplasmaregion .and. LBeltramiSeQuad ) then ! sequential quadratic programming (SQP); construct minimum energy with constrained helicity;
   
   lastcpu = GETTIME

   NLinearConstraints = 0 ! no linear constraints;
   
   NNonLinearConstraints = 1 ! single non-linear constraint = conserved helicity;
   
   LDA = max(1,NLinearConstraints)
   
   LDCJ = max(1,NNonLinearConstraints)
   
   LDR = NN
   
   RALLOCATE(LinearConstraintMatrix,(1:LDA,1:1)) ! linear constraint matrix;
   
   RALLOCATE(LowerBound,(1:NN+NLinearConstraints+NNonLinearConstraints)) ! lower bounds on variables, linear constraints and non-linear constraints;
   RALLOCATE(UpperBound,(1:NN+NLinearConstraints+NNonLinearConstraints)) ! upper bounds on variables, linear constraints and non-linear constraints;
   
   LowerBound(                       1 : NN                                          ) = -1.0E+21       !   variable constraints; no constraint;
   UpperBound(                       1 : NN                                          ) = +1.0E+21       !
   LowerBound( NN+                   1 : NN+NLinearConstraints                       ) = -1.0E+21       !     linear constraints; no constraint;
   UpperBound( NN+                   1 : NN+NLinearConstraints                       ) = +1.0E+21       !
   LowerBound( NN+NLinearConstraints+1 : NN+NLinearConstraints+NNonLinearConstraints ) = helicity(lvol) ! non-linear constraints; enforce helicity constraint;
   UpperBound( NN+NLinearConstraints+1 : NN+NLinearConstraints+NNonLinearConstraints ) = helicity(lvol) !
   
   iterations = 0 ! iteration counter;
   
   IALLOCATE(Istate,(1:NN+NLinearConstraints+NNonLinearConstraints))
   
   RALLOCATE(constraintfunction,(1:NNonLinearConstraints)) ! constraint functions;
   
   RALLOCATE(constraintgradient,(1:LDCJ,1:NN)) ! derivatives of constraint functions;
   
   RALLOCATE(multipliers,(1:NN+NLinearConstraints+NNonLinearConstraints)) ! Lagrange multipliers ?;
   
   objectivefunction = zero ! objective function;
   
   RALLOCATE(objectivegradient,(1:NN)) ! derivatives of objective function;
   
   RALLOCATE(RS,(1:LDR,1:NN))
   
   ideriv = 0 ; dpsi(1:2) = (/ dtflux(lvol), dpflux(lvol) /) ! these are also used below; 26 Feb 13;

   packorunpack = 'P'
   CALL(ma02aa,up00aa,( packorunpack, lvol, NN, xi(1:NN), dpsi(1:2), ideriv ))

   IALLOCATE(NEEDC,(1:NNonLinearConstraints))
   
   LIWk = 3*NN + NLinearConstraints + 2*NNonLinearConstraints ! workspace;
   IALLOCATE(IWk,(1:LIWk))       ! workspace;
   
   LRWk = 2*NN**2 + NN * NLinearConstraints + 2 * NN * NNonLinearConstraints + 21 * NN + 11 * NLinearConstraints + 22 * NNonLinearConstraints + 1 ! workspace;
   RALLOCATE(RWk,(1:LRWk))                                                              ! workspace;
   
   irevcm = 0 ; ie04uff = 1 ! reverse communication loop control; ifail error flag;
   
! supply optional parameters to E04UFF;
   
   call E04UEF('Nolist')               ! turn of screen output;
   call E04UEF('Print Level = 0')      ! turn of screen output;
   call E04UEF('Derivative Level = 3') ! assume all derivatives are provided by user;
   call E04UEF('Verify Level = -1')    ! do not verify derivatives using finite-differences; default is Verify Level = 0, which does verify gradients;
   write(optionalparameter,'("Major Iteration Limit = "i9)') 2**2 * max( 50, 3 * ( NN + NLinearConstraints ) + 10 * NNonLinearConstraints )
   call E04UEF(optionalparameter)
   
! pre-calculate some matrix vector products;

   MBpsi(1:NN) =                         matmul( dMB(1:NN,1: 2), dpsi(1:2) )
   MEpsi(1:NN) =                         matmul( dME(1:NN,1: 2), dpsi(1:2) )

   psiMCpsi    = half * sum( dpsi(1:2) * matmul( dMC(1: 2,1: 2), dpsi(1:2) ) )
   psiMFpsi    = half * sum( dpsi(1:2) * matmul( dMF(1: 2,1: 2), dpsi(1:2) ) )


   do ! reverse communication loop;

    
    call E04UFF( irevcm, &
                 NN, NLinearConstraints, NNonLinearConstraints, LDA, LDCJ, LDR, &
                 LinearConstraintMatrix(1:LDA,1:1), &
                 LowerBound(1:NN+NLinearConstraints+NNonLinearConstraints), UpperBound(1:NN+NLinearConstraints+NNonLinearConstraints), &
                 iterations, Istate(1:NN+NLinearConstraints+NNonLinearConstraints), &
                 constraintfunction(1:NNonLinearConstraints), constraintgradient(1:LDCJ,1:NN), &
                 multipliers(1:NN+NLinearConstraints+NNonLinearConstraints), &
                 objectivefunction, objectivegradient(1:NN), &
                 RS(1:LDR,1:NN), &
                 xi(1:NN), &
                 NEEDC(1:NNonLinearConstraints), IWk(1:LIWk), LIWk, RWk(1:LRWk), LRWk, ie04uff )

    if( irevcm.eq.1 .or. irevcm.eq.2 .or. irevcm.eq.3 ) Mxi(1:NN) = matmul( dMA(1:NN,1:NN), xi(1:NN) ) ! calculate objective  functional and/or gradient;
    if( irevcm.eq.4 .or. irevcm.eq.5 .or. irevcm.eq.6 ) Mxi(1:NN) = matmul( dMD(1:NN,1:NN), xi(1:NN) ) ! calculate constraint functional and/or gradient;
    
    if( irevcm.eq.1 .or. irevcm.eq.3 ) objectivefunction       = half * sum( xi(1:NN) * Mxi(1:NN) ) + sum( xi(1:NN) * MBpsi(1:NN) ) + psiMCpsi
    if( irevcm.eq.2 .or. irevcm.eq.3 ) objectivegradient(1:NN) =                        Mxi(1:NN)   +                 MBpsi(1:NN)
    
    if( irevcm.eq.4 .or. irevcm.eq.6 .and. NEEDC(1).gt.0 ) then
     constraintfunction(1     ) = half * sum( xi(1:NN) * Mxi(1:NN) ) + sum( xi(1:NN) * MEpsi(1:NN) ) + psiMFpsi
    endif
     
    if( irevcm.eq.5 .or. irevcm.eq.6 .and. NEEDC(1).gt.0 ) then
     constraintgradient(1,1:NN) =                        Mxi(1:NN)   +                 MEpsi(1:NN) 
    endif
    
    if( irevcm.eq.0 ) then ! final exit;

     cput = GETTIME

     select case(ie04uff)
     case( :-1 )  
                    write(ounit,1010) cput-cpus, myid, lvol, ie04uff, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "user enforced termination ;     "
     case(   0 )  
      if( Wma02aa ) write(ounit,1010) cput-cpus, myid, lvol, ie04uff, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "success ;                       "
     case(   1 )  
      if( Wma02aa ) write(ounit,1010) cput-cpus, myid, lvol, ie04uff, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "not converged;                  "
     case(   2 )  
                    write(ounit,1010) cput-cpus, myid, lvol, ie04uff, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "infeasible (linear) ;           "
     case(   3 ) 
                    write(ounit,1010) cput-cpus, myid, lvol, ie04uff, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "infeasible (nonlinear) ;        "
     case(   4 )   
                    write(ounit,1010) cput-cpus, myid, lvol, ie04uff, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "major iteration limit reached ; "
     case(   6 )  
      if( Wma02aa ) write(ounit,1010) cput-cpus, myid, lvol, ie04uff, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "could not be improved ;         "
     case(   7 )  
                    write(ounit,1010) cput-cpus, myid, lvol, ie04uff, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "derivatives appear incorrect ;  "
     case(   9 )  
                    write(ounit,1010) cput-cpus, myid, lvol, ie04uff, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "input error ;                   "
     case default ; FATALMESS(ma02aa,.true.,illegal ifail returned by E04UFF)
     end select

     if( irevcm.eq.0 .or. irevcm.eq.1 .or. irevcm.eq.6 ) ImagneticOK(lvol) = .true. ! set error flag; used elsewhere; 26 Feb 13;

     mu(lvol) = multipliers( NN + NLinearConstraints + NNonLinearConstraints ) ! helicity multiplier, or so it seems: NAG document is not completely clear;

     exit ! sequential quadratic programming method of constructing Beltrami field is finished;

    endif ! end of if( irevcm.eq.0 ) then;

   enddo ! end of do ! reverse communication loop;

   
   DEALLOCATE(RWk)
   DEALLOCATE(IWk)
   DEALLOCATE(NEEDC)
   DEALLOCATE(RS)
   DEALLOCATE(objectivegradient)
   DEALLOCATE(multipliers)
   DEALLOCATE(constraintgradient)
   DEALLOCATE(constraintfunction)
   DEALLOCATE(Istate)
   DEALLOCATE(LowerBound)
   DEALLOCATE(UpperBound)
   DEALLOCATE(LinearConstraintMatrix)


   packorunpack = 'U'
   CALL(ma02aa,up00aa( packorunpack, lvol, NN, xi(1:NN), dpsi(1:2), ideriv ))
   
   lBBintegral(lvol) = half * sum( xi(1:NN) * matmul( dMA(1:NN,1:NN), xi(1:NN) ) ) + sum( xi(1:NN) * MBpsi(1:NN) ) + psiMCpsi
   lABintegral(lvol) = half * sum( xi(1:NN) * matmul( dMD(1:NN,1:NN), xi(1:NN) ) ) + sum( xi(1:NN) * MEpsi(1:NN) ) + psiMFpsi

   
   solution(1:NN,0) = xi(1:NN)


  endif ! end of if( LBeltramiSeQuad ) then;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Lplasmaregion .and. LBeltramiNewton ) then
   
   lastcpu = GETTIME
   
   xi(0) = mu(lvol) ! initialize; helicity multiplier is treated as an independent degree-of-freedom;
   
   ideriv = 0 ; dpsi(1:2) = (/ dtflux(lvol), dpflux(lvol) /) ! these are also used below; 26 Feb 13;
   
   packorunpack = 'P'
   CALL(ma02aa,up00aa( packorunpack, lvol, NN, xi(1:NN), dpsi(1:2), ideriv ))
   
   pNN = NN + 1 ; Ldfmuaa = pNN ; tol = mupftol ; lengthwork = pNN * ( pNN+13 ) / 2

! pre-calculate some matrix vector products; these are used in df00aa;

   MBpsi(1:NN) =                         matmul( dMB(1:NN,1: 2), dpsi(1:2) )
   MEpsi(1:NN) =                         matmul( dME(1:NN,1: 2), dpsi(1:2) )
   psiMCpsi    = half * sum( dpsi(1:2) * matmul( dMC(1: 2,1: 2), dpsi(1:2) ) )
   psiMFpsi    = half * sum( dpsi(1:2) * matmul( dMF(1: 2,1: 2), dpsi(1:2) ) )
   
   ic05pbf = 1
   call C05PBF( df00aa, pNN, xi(0:NN), Fxi(0:NN), DFxi(0:NN,0:NN), Ldfmuaa, tol, work(1:lengthwork), lengthwork, ic05pbf )
   
   NewtonError = maxval( abs( Fxi(0:NN) ) )

   mu(lvol) = xi(0)

   packorunpack = 'U' ; ideriv = 0
   CALL(ma02aa,up00aa( packorunpack, lvol, NN, xi(1:NN), dpsi(1:2), ideriv ) )

   cput = GETTIME
   select case( ic05pbf )
   case(-1 )    
                  write(ounit,1020) cput-cpus, myid, lvol, ic05pbf, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, NewtonError, "user termination ;"
   case( 0 )    
    if( Wma02aa ) write(ounit,1020) cput-cpus, myid, lvol, ic05pbf, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, NewtonError, "success ;         "
   case( 1 )    
                  write(ounit,1020) cput-cpus, myid, lvol, ic05pbf, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, NewtonError, "input error ;     "
   case( 2 )    
                  write(ounit,1020) cput-cpus, myid, lvol, ic05pbf, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, NewtonError, "max. evaluations ;"
   case( 3 )    
                  write(ounit,1020) cput-cpus, myid, lvol, ic05pbf, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, NewtonError, "xtol too small ;  "
   case( 4 )    
                  write(ounit,1020) cput-cpus, myid, lvol, ic05pbf, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, NewtonError, "bad progress ;    "
   case default 
    FATALMESS(ma02aa, .true., illegal ifail returned byC05PBF)
   end select

#ifdef DEBUG
   xo(1:NN) = xi(1:NN) ! save original for comparison; 26 Feb 13;
   packorunpack = 'P' ; ideriv = 0
   CALL(ma02aa,up00aa( packorunpack, lvol, NN, xi(1:NN), dpsi(1:2), ideriv ))
   FATALMESS(ma02aa, sum(abs(xi(1:NN)-xo(1:NN)))/NN.gt.vsmall, un/packing routine is incorrect)
#endif

  !if( NewtonError.lt.mupftol ) then
    ImagneticOK(lvol) = .true. 
  !endif
   
   lBBintegral(lvol) = half * sum( xi(1:NN) * matmul( dMA(1:NN,1:NN), xi(1:NN) ) ) + sum( xi(1:NN) * MBpsi(1:NN) ) + psiMCpsi
   lABintegral(lvol) = half * sum( xi(1:NN) * matmul( dMD(1:NN,1:NN), xi(1:NN) ) ) + sum( xi(1:NN) * MEpsi(1:NN) ) + psiMFpsi

   solution(1:NN,0) = xi(1:NN)

   endif ! end of if( LBeltramiNewton ) then
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Lplasmaregion .and. LBeltramiLinear ) then ! assume Beltrami field is parameterized by helicity multiplier (and poloidal flux);
   
   lastcpu = GETTIME

   Xmudp(1:2) = xoffset + (/ mu(lvol), dpflux(lvol) /) ! initial guess for degrees of freedom; offset from zero so that relative error is small;

   select case( Lconstraint ) ! need to set Nmudp = degrees of freedom;

   case( 0 )    ;                                     Nmudp = 1 ! multiplier &  poloidal flux NOT varied                               ; dummy value;
   case( 1 )    ; if( Lcoordinatesingularity ) then ; Nmudp = 1 ! multiplier                  IS  varied to match       outer transform;
    ;             else                              ; Nmudp = 2 ! multiplier &  poloidal flux ARE varied to match inner/outer transform;
    ;             endif                                         
   case( 2 )    ;                                     Nmudp = 1 ! multiplier                  IS  varied to match             helicity ;
   case default ; FATALMESS(ma02aa, .true., supplied Lconstraint is not supported)
   end select

   select case( Lconstraint ) 
    
   case( 0 ) ! need only call mp00ac once, to calculate Beltrami field for given helicity multiplier and enclosed fluxes; 28 Jan 13;
    
    iflag = 1 ; Ldfjac = Nmudp
    
    WCALL(ma02aa,mp00ac,( Nmudp, Xmudp(1:Nmudp), Fmudp(1:Nmudp), Dmudp(1:Ldfjac,1:Nmudp), Ldfjac, iflag )) ! calculate Beltrami field;
   
    ic05pcf = 0 ; nfev = 1 ; njev = 0

    helicity(lvol) = lABintegral(lvol) ! this was computed in mp00ac; 26 Feb 13;
 
   case( 1 ) ! will iteratively call mp00ac, to calculate Beltrami field with prescribed interface rotational transform; 28 Jan 13;
    
   !tol = mupftol ; Ldfjac = Nmudp ; LRR = Nmudp * ( Nmudp+1 ) / 2 ; mode = 0 ; diag(1:2) = zero ; factor = 1.0e-02 ! PLEASE CHECK WHAT FACTOR IS; 18 Feb 13;
   !tol = mupftol ; Ldfjac = Nmudp ; LRR = Nmudp * ( Nmudp+1 ) / 2 ; mode = 0 ; diag(1:2) = zero ; factor = 1.0e+02 ! PLEASE CHECK WHAT FACTOR IS; 20 Jun 14;
    tol = mupftol ; Ldfjac = Nmudp ; LRR = Nmudp * ( Nmudp+1 ) / 2 ; mode = 0 ; diag(1:2) = zero ; factor = 1.0e+00 ! PLEASE CHECK WHAT FACTOR IS; 04 Dec 14;
    
    maxfev = mupfits ; nprint = 0 ; nfev = 0 ; njev = 0
    
    if( Wma02aa ) then 
     cput = GETTIME
     write(ounit,'("ma02aa : ",f10.2," : myid=",i3," ; calling C05PCF, which calls mp00ac ;")') cput-cpus, myid
    endif
    
    ic05pcf = 1
    WCALL(ma02aa,C05PCF,( mp00ac, Nmudp, Xmudp(1:Nmudp), Fmudp(1:Nmudp), Dmudp(1:Ldfjac,1:Nmudp), Ldfjac, tol, maxfev, diag(1:Nmudp), mode, factor, nprint, &
                 nfev, njev, RR(1:LRR), LRR, QTF(1:Nmudp), WK(1:Nmudp,1:4), ic05pcf ))

    select case( ic05pcf )
    case( 0: ) ;     mu(lvol) = Xmudp(1) - xoffset
     ;         ; dpflux(lvol) = Xmudp(2) - xoffset
    case( -2 ) ! mu and dpflux have been updated in mp00ac; early termination; 18 Apr 13;
    end select

    helicity(lvol) = lABintegral(lvol) ! this was computed in mp00ac; 26 Feb 13;

    if( Lconstraint.eq.1 ) then
     iflag = 2 ; Ldfjac = Nmudp ! call tr00ab once more to ensure that the derivatives of the field & transform wrt mu and dpflux are calculated; 20 Jun 14;
     WCALL(ma02aa,mp00ac,( Nmudp, Xmudp(1:Nmudp), Fmudp(1:Nmudp), Dmudp(1:Ldfjac,1:Nmudp), Ldfjac, iflag )) ! calculate Beltrami field;
    endif
    
   case default
    
    FATALMESS(ma02aa, .true., selected Lconstraint not supported)
    
   end select
     
   cput = GETTIME
   select case(ic05pcf)
   case(    0   )  
    if( Wma02aa ) write(ounit,1040) cput-cpus, myid, lvol, ic05pcf, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "success         ", Fmudp(1:Nmudp)
   case(   -2   )  
    if( Wma02aa ) write(ounit,1040) cput-cpus, myid, lvol, ic05pcf, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "|F| < mupftol   ", Fmudp(1:Nmudp)
   case(   -1   )  
                  write(ounit,1040) cput-cpus, myid, lvol, ic05pcf, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "Beltrami fail   ", Fmudp(1:Nmudp)
   case(    1   )  
                  write(ounit,1040) cput-cpus, myid, lvol, ic05pcf, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "input error     ", Fmudp(1:Nmudp)
   case(    2   )  
                  write(ounit,1040) cput-cpus, myid, lvol, ic05pcf, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "consider restart", Fmudp(1:Nmudp)
   case(    3   )  
                  write(ounit,1040) cput-cpus, myid, lvol, ic05pcf, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "xtol too small  ", Fmudp(1:Nmudp)
   case(    4:5 )  
                  write(ounit,1040) cput-cpus, myid, lvol, ic05pcf, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "bad progress    ", Fmudp(1:Nmudp)
   case default    
                  write(ounit,1040) cput-cpus, myid, lvol, ic05pcf, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "illegal ifail   ", Fmudp(1:Nmudp)
    FATALMESS(ma02aa, .true., illegal ifail returned by C05PCF)
   end select

   endif ! end of if( LBeltramiLinear ) then;

 !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

   if( Lvacuumregion ) then ! vacuum region; ! vacuum scalar potential can always be solved by linear method; 21 Apr 13;

    lastcpu = GETTIME

    Xmudp(1:2) = xoffset + (/ curtor, curpol /) ! initial guess for degrees of freedom; offset from zero so that relative error is small;

    select case( Lconstraint ) ! 03 Apr 13;

    case( 0 ) ! need only call mp00ac once, to calculate Beltrami field for given helicity multiplier and enclosed fluxes; 28 Jan 13; ! TEMPORARY; 03 Apr 13;

     Nmudp = 1 ! ! 19 Apr 13;

     iflag = 1 ; Ldfjac = Nmudp

     WCALL(ma02aa,mp00ac,( Nmudp, Xmudp(1:Nmudp), Fmudp(1:Nmudp), Dmudp(1:Ldfjac,1:Nmudp), Ldfjac, iflag )) ! calculate Beltrami field;

     ic05pcf = 0 ; nfev =1 ; njev = 0 ! ! just provide some dummy values for screen output below; 19 Apr 13;

    case( 1 ) ! will iterative call mp00ac, to calculate Beltrami field with prescribed interface rotational transform; 28 Jan 13;

     Nmudp = 1

     tol = mupftol ; Ldfjac = Nmudp ; LRR = Nmudp * (Nmudp+1) / 2 ; mode = 0 ; diag(1:2) = zero ; factor = 1.0e+02 ! factor = 100 is suggested by NAG; 02 Sep 14;

     maxfev = mupfits ; nprint = 0 ; nfev = 0 ; njev = 0

     ic05pcf = 1
     WCALL(ma02aa,C05PCF,( mp00ac, Nmudp, Xmudp(1:Nmudp), Fmudp(1:Nmudp), Dmudp(1:Ldfjac,1:Nmudp), Ldfjac, tol, maxfev, diag(1:Nmudp), mode, factor, nprint, &
                  nfev, njev, RR(1:LRR), LRR, QTF(1:Nmudp), WK(1:Nmudp,1:4), ic05pcf ))

     select case( ic05pcf )
     case( 0: ) ; curtor = Xmudp(1) - xoffset
      ;         ;!curpol = Xmudp(2) - xoffset ! ! this is not changed; 22 Apr 13;
     case( -2 ) ! curtor has been updated in mp00ac; early termination; 21 Apr 13;
     end select

     if( Lconstraint.eq.1 ) then
      iflag = 2 ; Ldfjac = Nmudp ! call tr00ab once more to ensure that the derivatives of the field & transform wrt mu and dpflux are calculated; 20 Jun 14;
      WCALL(ma02aa,mp00ac,( Nmudp, Xmudp(1:Nmudp), Fmudp(1:Nmudp), Dmudp(1:Ldfjac,1:Nmudp), Ldfjac, iflag )) ! calculate Beltrami field;
     endif

    case default

     FATALMESS(ma02aa, .true., selected Lconstraint not supported)

    end select

    cput = GETTIME
    select case(ic05pcf)
    case(    0   ) ; if( Wma02aa ) write(ounit,1050) cput-cpus, myid, lvol, ic05pcf, curtor, cput-lastcpu, "success         ", Fmudp(1:Nmudp)
    case(   -2   ) ; if( Wma02aa ) write(ounit,1050) cput-cpus, myid, lvol, ic05pcf, curtor, cput-lastcpu, "|F| < mupftol   ", Fmudp(1:Nmudp)
    case(   -1   ) ;               write(ounit,1050) cput-cpus, myid, lvol, ic05pcf, curtor, cput-lastcpu, "Beltrami fail   ", Fmudp(1:Nmudp)
    case(    1   ) ;               write(ounit,1050) cput-cpus, myid, lvol, ic05pcf, curtor, cput-lastcpu, "input error     ", Fmudp(1:Nmudp)
    case(    2   ) ;               write(ounit,1050) cput-cpus, myid, lvol, ic05pcf, curtor, cput-lastcpu, "consider restart", Fmudp(1:Nmudp)
    case(    3   ) ;               write(ounit,1050) cput-cpus, myid, lvol, ic05pcf, curtor, cput-lastcpu, "xtol too small  ", Fmudp(1:Nmudp)
    case(    4:5 ) ;               write(ounit,1050) cput-cpus, myid, lvol, ic05pcf, curtor, cput-lastcpu, "bad progress    ", Fmudp(1:Nmudp)
    case default   ;               write(ounit,1050) cput-cpus, myid, lvol, ic05pcf, curtor, cput-lastcpu, "illegal ifail   ", Fmudp(1:Nmudp)
    end select
   endif ! end of if( Lvacuumregion );

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef DEBUG

   if( Lcheck.eq.2 .and. Lconstraint.eq.1 ) then ! perform finite difference check on derivative of rotational-transform; 01 Jul 14;

    if( Lplasmaregion ) then

     RALLOCATE(dsolution,(1:NN,0:2,-1:1,-1:1)) ! packed vector potential; 01 Jul 14;

     if( Lcoordinatesingularity ) then ; Nmudp = 1 ! multiplier                  IS  varied to match       transform      ;
     else                              ; Nmudp = 2 ! multiplier &  poloidal flux ARE varied to match inner/outer transform;
     endif

     Ldfjac = Nmudp ; dFmudp(-1:1,-1:1,1:2) = zero

     Xmudp(1:2) = xoffset + (/ mu(lvol), dpflux(lvol) /) ! initial guess for degrees of freedom; offset from zero so that relative error is small;

     imu = 0 ; idp = 0

     iflag = 2 ! iflag controls derivative calculation in mp00ac; analytic derivatives of rotational-transform are required; 01 Jul 14;

     CALL(ma02aa,mp00ac( Nmudp, Xmudp(1:Nmudp), dFmudp(imu,idp,1:Nmudp), oDmudp(1:Ldfjac,1:Nmudp), Ldfjac, iflag )) ! compute "exact" derivatives;

     dsolution(1:NN,1,0,0) = solution(1:NN,1) ! packed vector potential; derivative wrt mu    ; 20 Jun 14;
     dsolution(1:NN,2,0,0) = solution(1:NN,2) ! packed vector potential; derivative wrt dpflux; 20 Jun 14;

     jfinite = 0
     cput = GETTIME 
     write(ounit,2000) cput-cpus, myid, lvol, jfinite, "derivative", oDmudp(1:Ldfjac,1:Nmudp)

     do jfinite = -4,-2,+1 ; lfdiff = ten**jfinite

      do imu = -1, +1 ! centered differences;
       do idp = -1, +1 ! centered differences;

        if( Lcoordinatesingularity .and. idp.ne.0 ) cycle ! Beltrami field in volume with coordinate singularity does not depend on enclosed poloidal flux;

        if( abs(imu)+abs(idp).ne.1 ) cycle ! centered finite differences;

        Xmudp(1:2) = xoffset + (/ mu(lvol), dpflux(lvol) /) + (/ imu, idp /) * lfdiff * half 

       !write(ounit,'("ma02aa : " 10x " : myid="3x" : imu="i3" ; idp="i3" ; Xmudp="2es13.5" ;")') imu, idp, Xmudp(1:2) - xoffset

        iflag = 1 ! analytic derivatives of rotational-transform are not required;

        CALL(ma02aa,mp00ac( Nmudp, Xmudp(1:Nmudp), dFmudp(imu,idp,1:Nmudp), Dmudp(1:Ldfjac,1:Nmudp), Ldfjac, iflag )) ! compute function values only;

        dsolution(1:NN,0,imu,idp) = solution(1:NN,0)

       enddo ! end of do idp;
      enddo ! end of do imu;

      ;                Dmudp(1:Nmudp,1) = ( dFmudp( 1, 0, 1:Nmudp) - dFmudp(-1, 0, 1:Nmudp) ) / lfdiff ! derivative wrt helicity multiplier   ;
      if( Nmudp.eq.2 ) Dmudp(1:Nmudp,2) = ( dFmudp( 0, 1, 1:Nmudp) - dFmudp( 0,-1, 1:Nmudp) ) / lfdiff ! derivative wrt enclosed poloidal flux;

      cput = GETTIME 
     !write(ounit,2000) cput-cpus, myid, lvol, jfinite, " error    ", Dmudp(1:Ldfjac,1:Nmudp) - oDmudp(1:Ldfjac,1:Nmudp)
      write(ounit,2000) cput-cpus, myid, lvol, jfinite, " estimate ", Dmudp(1:Ldfjac,1:Nmudp) ! 04 Dec 14;

    enddo ! end of do jfinite;
    
    cput = GETTIME 
    write(ounit,2000) cput-cpus
    
    DEALLOCATE(dsolution)
    
2000 format("ma02aa : ":,f10.2," :":" myid=",i3," : lvol=",i3," ; jj=",i3," ; "a10" : DF=["es23.15" ,"es23.15" ,"es23.15" ,"es23.15" ] ;") ! 20 Jun 14;
    
   endif ! end of if( Lplasmaregion ) ; 20 Jun 14;
   
   if( Lvacuumregion ) then ! perform finite difference check on derivative of rotational-transform;
    
    RALLOCATE(dsolution,(1:NN,0:2,-1:1,-1:1))
    
    Nmudp = 1 ! toroidal plasma current, curtor, is varied to match transform constraint;
    
    Ldfjac = Nmudp ; dFmudp(-1:1,-1:1,1:2) = zero
    
    Xmudp(1:2) = xoffset + (/ curtor, curpol /) ! initial guess for degrees of freedom; offset from zero so that relative error is small;
    
    iflag = 2 ! derivatives are required;
    
    CALL(ma02aa,mp00ac( Nmudp, Xmudp(1:Nmudp), dFmudp( 0 , 0 ,1:Nmudp), Dmudp(1:Ldfjac,1:Nmudp), Ldfjac, iflag )) ! compute "exact" derivatives;
    
    dsolution(1:NN,1,0,0) = solution(1:NN,1)
    
    jfinite = 0
    cput = GETTIME 
    write(ounit,2000) cput-cpus, myid, lvol, jfinite, Dmudp(1:Ldfjac,1:Nmudp)
    
    do jfinite = -6,-2,+1
     
     lfdiff = ten**jfinite
     
     do icurtor = -1, +1, 2 ! centered differences;
      
      Xmudp(1:2) = xoffset + (/ curtor, curpol /) + (/ icurtor, 0 /) * lfdiff * half 
      
      iflag = 1 ! derivatives are not required;
      
      CALL(ma02aa,mp00ac( Nmudp, Xmudp(1:Nmudp), dFmudp(icurtor,0,1:Nmudp), Dmudp(1:Ldfjac,1:Nmudp), Ldfjac, iflag )) ! compute function values only;
      
      dsolution(1:NN,0,icurtor,0) = solution(1:NN,0)
      
     enddo ! end of do icurtor; 21 Apr 13;
     
     ;                            Dmudp(1:Nmudp,1) = ( dFmudp( 1, 0, 1:Nmudp) - dFmudp(-1, 0, 1:Nmudp) ) / lfdiff ! derivative wrt curtor;
    
     cput = GETTIME 
     write(ounit,2000) cput-cpus, myid, lvol, jfinite, Dmudp(1:Ldfjac,1:Nmudp)
     
     write(ounit,'("ma02aa : ", 10x ," : "99es13.5)')  dsolution(1:min(NN,10),1, 0,0) ! analytic derivative; 21 Apr 13;
     write(ounit,'("ma02aa : ", 10x ," : "99es13.5)') (dsolution(1:min(NN,10),0,+1,0)-dsolution(1:min(NN,10),0,-1,0)) / lfdiff
     
    enddo ! end of do jfinite;
    
   !cput = GETTIME 
    write(ounit,2000)!cput-cpus
    
    DEALLOCATE(dsolution)
    
   endif ! end of if( Lvacuumregion ) ; 20 Jun 14;
  
  endif ! end of if( Lcheck.eq.2 ) ; 01 Jul 14;

#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(ma02aa)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
1010 format("ma02aa : ",f10.2," : myid=",i3," ; lvol=",i3," ; SQP    : ie04uff=",i3," hel="es12.4" mu="es12.4" dpflux="es12.4" time="f9.1" ":,a36)
1020 format("ma02aa : ",f10.2," : myid=",i3," ; lvol=",i3," ; Newton : ic05pbf=",i3," hel="es12.4" mu="es12.4" dpflux="es12.4" time="f9.1" ; "&
  "error="es7.0" ; ":,a18)
1040 format("ma02aa : ",f10.2," : myid=",i3," ; lvol=",i3," ; Linear : ic05pcf=",i3," hel="es12.4" mu="es12.4" dpflux="es12.4" time="f9.1" ; "&
  :,a16" ; F="2es08.0)
1050 format("ma02aa : ",f10.2," : myid=",i3," ; lvol=",i3," ; Linear : ic05pcf=",i3,"     "  12x " I ="es12.4"        "  12x " time="f9.1" ; "&
  :,a16" ; F="2es08.0)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine ma02aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
