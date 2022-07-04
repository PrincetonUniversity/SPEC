!> \defgroup grp_solver_driver Solver/Driver
!>
!> \file
!> \brief Constructs Beltrami field in given volume consistent with flux, helicity, rotational-transform and/or parallel-current constraints.

!> \brief Constructs Beltrami field in given volume consistent with flux, helicity, rotational-transform and/or parallel-current constraints.
!> \ingroup grp_solver_driver
!>
!> @param[in] lvol index of nested volume for which to run this
!> @param[in] NN   number of degrees of freedom in the (packed format) vector potential;
subroutine ma02aa(lvol, NN)
    use mod_kinds, only: wp => dp
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    use constants, only: zero, half, one, ten

    use numerical, only: vsmall, small

    use fileunits, only: ounit

    use inputlist, only: Wmacros, Wma02aa, &
                         Lconstraint, mu, helicity, &
                         mupftol, mupfits, Lrad, Lcheck

    use cputiming

    use allglobal, only: ncpu, myid, cpus, MPI_COMM_SPEC, &
                         Mvol, mn, im, in, &
                         LBlinear, LBnewton, LBsequad, &
                         !                       dMA, dMB, dMC, dMD, dME, dMF, solution, &
                         dMA, dMB, dMD, solution, &
                         !                       MBpsi, MEpsi, psiMCpsi, psiMFpsi, &
                         MBpsi, Ate, &
                         ImagneticOK, &
                         lBBintegral, lABintegral, &
                         ivol, Nfielddof, &
                         dtflux, dpflux, &
                         xoffset, &
                         Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, LocalConstraint

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef OPENMP
    USE OMP_LIB
#endif
    use mpi
    implicit none
    integer :: ierr, astat, ios, nthreads, ithread
    real(wp) :: cput, cpui, cpuo = 0 ! cpu time; cpu initial; cpu old; 31 Jan 13;

    integer, intent(in) :: lvol, NN

    integer :: ideriv
    real(wp) :: tol, dpsi(1:2), lastcpu
    character :: packorunpack

    integer :: Nxdof, Ndof, Ldfjac, iflag, maxfev, mode, LRR, nfev, njev, nprint, ihybrj
    real(wp) :: Xdof(1:2), Fdof(1:2), Ddof(1:2, 1:2), oDdof(1:2, 1:2)
    real(wp) :: factor, diag(1:2), RR(1:2*(2 + 1)/2), QTF(1:2), wk(1:2, 1:4)

    integer :: irevcm

    integer :: pNN

    real(wp) :: xi(0:NN), Fxi(0:NN), xo(0:NN), Mxi(1:NN)

    external :: mp00ac

#ifdef DEBUG
    integer :: ixx, jxx, jfinite ! computing finite-difference derivatives of \iota wrt \mu and \Delta \psi_p;
    real(wp) :: lfdiff, dFdof(-1:1, -1:1, 1:2)
    real(wp), allocatable :: dsolution(:, :, :, :)
#endif

!required for hybrj1;
    integer :: ihybrj1, Ldfmuaa, lengthwork
    real(wp) :: NewtonError
    real(wp), allocatable :: DFxi(:, :), work(:)
    external :: df00ab

! required for E04UFF;
    integer :: NLinearConstraints, NNonLinearConstraints, LDA, LDCJ, LDR, iterations, LIWk, LRWk, ie04uff
    integer, allocatable :: Istate(:), NEEDC(:), IWk(:)
    real(wp) :: objectivefunction
    real(wp), allocatable :: LinearConstraintMatrix(:, :), LowerBound(:), UpperBound(:)
    real(wp), allocatable :: constraintfunction(:), constraintgradient(:, :), multipliers(:), objectivegradient(:), RS(:, :), RWk(:)
    character :: optionalparameter*33

    cpui = MPI_WTIME()
    cpuo = cpui
#ifdef OPENMP
    nthreads = omp_get_max_threads()
#else
    nthreads = 1
#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef DEBUG

    if (lvol .lt. 1 .or. lvol .gt. Mvol) then
        write (6, '("ma02aa :      fatal : myid=",i3," ; lvol.lt.1 .or. lvol.gt.Mvol ; illegal lvol ;")') myid
        call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
        stop "ma02aa : lvol.lt.1 .or. lvol.gt.Mvol : illegal lvol  ;"
    end if

#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    ivol = lvol ! various subroutines (e.g. mp00ac, df00ab) that may be called below require volume identification, but the argument list is fixed by NAG;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **sequential quadratic programming**
!> <ul>
!> <li> Only relevant if \c LBsequad=T . See \c LBeltrami for details. </li>
!> <li> Documentation on the implementation of \c E04UFF is under construction. </li>
!> </ul>

    if (LBsequad) then ! sequential quadratic programming (SQP); construct minimum energy with constrained helicity;
        lastcpu = MPI_WTIME()

        NLinearConstraints = 0 ! no linear constraints;

        NNonLinearConstraints = 1 ! single non-linear constraint = conserved helicity;

        LDA = max(1, NLinearConstraints)

        LDCJ = max(1, NNonLinearConstraints)

        LDR = NN

        allocate (LinearConstraintMatrix(1:LDA, 1:1), stat=astat)
        LinearConstraintMatrix(1:LDA, 1:1) = zero
        ! linear constraint matrix;

        allocate (LowerBound(1:NN + NLinearConstraints + NNonLinearConstraints), stat=astat)
        LowerBound(1:NN + NLinearConstraints + NNonLinearConstraints) = zero
        ! lower bounds on variables, linear constraints and non-linear constraints;

        allocate (UpperBound(1:NN + NLinearConstraints + NNonLinearConstraints), stat=astat)
        UpperBound(1:NN + NLinearConstraints + NNonLinearConstraints) = zero
        ! upper bounds on variables, linear constraints and non-linear constraints;

        LowerBound(1:NN) = -1.0E+21       !   variable constraints; no constraint;
        UpperBound(1:NN) = +1.0E+21       !
        LowerBound(NN + 1:NN + NLinearConstraints) = -1.0E+21       !     linear constraints; no constraint;
        UpperBound(NN + 1:NN + NLinearConstraints) = +1.0E+21       !
        LowerBound(NN + NLinearConstraints + 1:NN + NLinearConstraints + NNonLinearConstraints) = helicity(lvol) ! non-linear constraints; enforce helicity constraint;
        UpperBound(NN + NLinearConstraints + 1:NN + NLinearConstraints + NNonLinearConstraints) = helicity(lvol) !

        iterations = 0 ! iteration counter;

        allocate (Istate(1:NN + NLinearConstraints + NNonLinearConstraints), stat=astat)
        Istate(1:NN + NLinearConstraints + NNonLinearConstraints) = 0

        allocate (constraintfunction(1:NNonLinearConstraints), stat=astat)
        constraintfunction(1:NNonLinearConstraints) = zero
        ! constraint functions;

        allocate (constraintgradient(1:LDCJ, 1:NN), stat=astat)
        constraintgradient(1:LDCJ, 1:NN) = zero
        ! derivatives of constraint functions;

        allocate (multipliers(1:NN + NLinearConstraints + NNonLinearConstraints), stat=astat)
        multipliers(1:NN + NLinearConstraints + NNonLinearConstraints) = zero
        ! Lagrange multipliers ?;

        objectivefunction = zero ! objective function;

        allocate (objectivegradient(1:NN), stat=astat)
        objectivegradient(1:NN) = zero
        ! derivatives of objective function;

        allocate (RS(1:LDR, 1:NN), stat=astat)
        RS(1:LDR, 1:NN) = zero

        ideriv = 0; dpsi(1:2) = (/dtflux(lvol), dpflux(lvol)/) ! these are also used below;

        packorunpack = 'P'

        cput = MPI_WTIME()
        Tma02aa = Tma02aa + (cput - cpuo)
        call packab(packorunpack, lvol, NN, xi(1:NN), ideriv)
        cpuo = MPI_WTIME()

        allocate (NEEDC(1:NNonLinearConstraints), stat=astat)
        NEEDC(1:NNonLinearConstraints) = 0

        LIWk = 3*NN + NLinearConstraints + 2*NNonLinearConstraints ! workspace;

        allocate (IWk(1:LIWk), stat=astat)
        IWk(1:LIWk) = 0
        ! workspace;

        LRWk = 2*NN**2 + NN*NLinearConstraints + 2*NN*NNonLinearConstraints + 21*NN + 11*NLinearConstraints + 22*NNonLinearConstraints + 1 ! workspace;

        allocate (RWk(1:LRWk), stat=astat)
        RWk(1:LRWk) = zero
        ! workspace;

        irevcm = 0; ie04uff = 1 ! reverse communication loop control; ifail error flag;

! supply optional parameters to E04UFF; NAG calls commented out (this part of the code so far not used); 17 Nov 17

!   call E04UEF('Nolist')               ! turn of screen output;
!   call E04UEF('Print Level = 0')      ! turn of screen output;
!   call E04UEF('Derivative Level = 3') ! assume all derivatives are provided by user;
!   call E04UEF('Verify Level = -1')    ! do not verify derivatives using finite-differences; default is Verify Level = 0, which does verify gradients;
!   write(optionalparameter,'("Major Iteration Limit = "i9)') 2**2 * max( 50, 3 * ( NN + NLinearConstraints ) + 10 * NNonLinearConstraints )
!   call E04UEF(optionalparameter)

! pre-calculate some matrix vector products;

        MBpsi(1:NN) = matmul(dMB(1:NN, 1:2), dpsi(1:2))
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
!    if( irevcm.eq.1 .or. irevcm.eq.2 .or. irevcm.eq.3 ) Mxi(1:NN) = matmul( dMA(1:NN,1:NN), xi(1:NN) ) ! calculate objective  functional and/or gradient;
!    if( irevcm.eq.4 .or. irevcm.eq.5 .or. irevcm.eq.6 ) Mxi(1:NN) = matmul( dMD(1:NN,1:NN), xi(1:NN) ) ! calculate constraint functional and/or gradient;
!
!    if( irevcm.eq.1 .or. irevcm.eq.3 ) objectivefunction       = half * sum( xi(1:NN) * Mxi(1:NN) ) + sum( xi(1:NN) * MBpsi(1:NN) ) + psiMCpsi
!    if( irevcm.eq.2 .or. irevcm.eq.3 ) objectivegradient(1:NN) =                        Mxi(1:NN)   +                 MBpsi(1:NN)
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

        deallocate (RWk, stat=astat)

        deallocate (IWk, stat=astat)

        deallocate (NEEDC, stat=astat)

        deallocate (RS, stat=astat)

        deallocate (objectivegradient, stat=astat)

        deallocate (multipliers, stat=astat)

        deallocate (constraintgradient, stat=astat)

        deallocate (constraintfunction, stat=astat)

        deallocate (Istate, stat=astat)

        deallocate (LowerBound, stat=astat)

        deallocate (UpperBound, stat=astat)

        deallocate (LinearConstraintMatrix, stat=astat)

        packorunpack = 'U'

        cput = MPI_WTIME()
        Tma02aa = Tma02aa + (cput - cpuo)
        call packab(packorunpack, lvol, NN, xi(1:NN), ideriv)
        cpuo = MPI_WTIME()

        lBBintegral(lvol) = half*sum(xi(1:NN)*matmul(dMA(1:NN, 1:NN), xi(1:NN))) + sum(xi(1:NN)*MBpsi(1:NN)) ! + psiMCpsi
        lABintegral(lvol) = half*sum(xi(1:NN)*matmul(dMD(1:NN, 1:NN), xi(1:NN))) ! + sum( xi(1:NN) * MEpsi(1:NN) ) ! + psiMFpsi

        solution(1:NN, 0) = xi(1:NN)

    end if ! end of if( LBsequad ) then;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **Newton method**
!> <ul>
!> <li> Only relevant if \c LBnewton=T . See \c LBeltrami for details. </li>
!> </ul>

    if (LBnewton) then

        lastcpu = MPI_WTIME()

        allocate (DFxi(0:NN, 0:NN), stat=astat)
        DFxi(0:NN, 0:NN) = zero

        allocate (work(1:(1 + NN)*(1 + NN + 13)/2), stat=astat)
        work(1:(1 + NN)*(1 + NN + 13)/2) = zero

        xi(0) = mu(lvol) ! initialize; helicity multiplier is treated as an independent degree-of-freedom;

        ideriv = 0; dpsi(1:2) = (/dtflux(lvol), dpflux(lvol)/) ! these are also used below;

        packorunpack = 'P'
!   CALL( ma02aa, packab ( packorunpack, lvol, NN, xi(1:NN), dpsi(1:2), ideriv ) )

        cput = MPI_WTIME()
        Tma02aa = Tma02aa + (cput - cpuo)
        call packab(packorunpack, lvol, NN, xi(1:NN), ideriv)
        cpuo = MPI_WTIME()

        pNN = NN + 1; Ldfmuaa = pNN; tol = mupftol; lengthwork = pNN*(pNN + 13)/2

! pre-calculate some matrix vector products; these are used in df00ab;

        MBpsi(1:NN) = matmul(dMB(1:NN, 1:2), dpsi(1:2))
!  MEpsi(1:NN) = zero !                  matmul( dME(1:NN,1: 2), dpsi(1:2) )
!  psiMCpsi    = zero ! half * sum( dpsi(1:2) * matmul( dMC(1: 2,1: 2), dpsi(1:2) ) )
!  psiMFpsi    = zero ! half * sum( dpsi(1:2) * matmul( dMF(1: 2,1: 2), dpsi(1:2) ) )

        call hybrj1(df00ab, pNN, xi(0:NN), Fxi(0:NN), DFxi(0:NN, 0:NN), Ldfmuaa, tol, ihybrj1, work(1:lengthwork), lengthwork)

        NewtonError = maxval(abs(Fxi(0:NN)))

        mu(lvol) = xi(0)

        packorunpack = 'U'; ideriv = 0

        cput = MPI_WTIME()
        Tma02aa = Tma02aa + (cput - cpuo)
        call packab(packorunpack, lvol, NN, xi(1:NN), ideriv)
        cpuo = MPI_WTIME()

        cput = MPI_WTIME()

        select case (ihybrj1)
        case (:-1)
            ; write (ounit, 1020) cput - cpus, myid, lvol, ihybrj1, helicity(lvol), mu(lvol), dpflux(lvol), cput - lastcpu, NewtonError, "user termination ;"
        case (1)
            if (Wma02aa) write (ounit, 1020) cput - cpus, myid, lvol, ihybrj1, helicity(lvol), mu(lvol), dpflux(lvol), cput - lastcpu, NewtonError, "success ;         "
        case (0)
            ; write (ounit, 1020) cput - cpus, myid, lvol, ihybrj1, helicity(lvol), mu(lvol), dpflux(lvol), cput - lastcpu, NewtonError, "input error ;     "
        case (2)
            ; write (ounit, 1020) cput - cpus, myid, lvol, ihybrj1, helicity(lvol), mu(lvol), dpflux(lvol), cput - lastcpu, NewtonError, "max. evaluations ;"
        case (3)
            if (NewtonError > tol) then
                ; write (ounit, 1020) cput - cpus, myid, lvol, ihybrj1, helicity(lvol), mu(lvol), dpflux(lvol), cput - lastcpu, NewtonError, "xtol too small ;  "
            end if
        case (4)
            if (NewtonError > tol) then
                ; write (ounit, 1020) cput - cpus, myid, lvol, ihybrj1, helicity(lvol), mu(lvol), dpflux(lvol), cput - lastcpu, NewtonError, "bad progress ;    "
            end if
        case default

            if (.true.) then
                write (6, '("ma02aa :      fatal : myid=",i3," ; .true. ; illegal ifail returned by hybrj1 ;")') myid
                call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                stop "ma02aa : .true. : illegal ifail returned by hybrj1  ;"
            end if

        end select

#ifdef DEBUG
        xo(1:NN) = xi(1:NN) ! save original for comparison;
        packorunpack = 'P'; ideriv = 0

        cput = MPI_WTIME()
        Tma02aa = Tma02aa + (cput - cpuo)
        call packab(packorunpack, lvol, NN, xi(1:NN), ideriv)
        cpuo = MPI_WTIME()

        if (sum(abs(xi(1:Nfielddof(lvol)) - xo(1:Nfielddof(lvol))))/Nfielddof(lvol) .gt. vsmall) then
            write (6, '("ma02aa :      fatal : myid=",i3," ; sum(abs(xi(1:Nfielddof(lvol))-xo(1:Nfielddof(lvol))))/Nfielddof(lvol).gt.vsmall ; un/packing routine is incorrect ;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "ma02aa : sum(abs(xi(1:Nfielddof(lvol))-xo(1:Nfielddof(lvol))))/Nfielddof(lvol).gt.vsmall : un/packing routine is incorrect  ;"
        end if

#endif

!if( NewtonError.lt.mupftol ) then
        ImagneticOK(lvol) = .true.
!endif

        lBBintegral(lvol) = half*sum(xi(1:NN)*matmul(dMA(1:NN, 1:NN), xi(1:NN))) + sum(xi(1:NN)*MBpsi(1:NN)) ! + psiMCpsi
        lABintegral(lvol) = half*sum(xi(1:NN)*matmul(dMD(1:NN, 1:NN), xi(1:NN))) ! + sum( xi(1:NN) * MEpsi(1:NN) ) ! + psiMFpsi

        solution(1:NN, 0) = xi(1:NN)

        deallocate (DFxi, stat=astat)

        deallocate (work, stat=astat)

    end if ! end of if( LBnewton ) then

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **linear method**

!> <ul>
!> <li> Only relevant if \c LBlinear=T . See \c LBeltrami for details. </li>
!> <li> The quantity \f$\mu\f$ is *not* not treated as a "magnetic" degree-of-freedom
!>       equivalent to in the degrees-of-freedom in the magnetic vector potential
!>       (as it strictly should be, because it is a Lagrange multiplier introduced to enforce the helicity constraint). </li>
!> <li> In this case, the Beltrami equation, \f$\nabla \times {\bf B} = \mu {\bf B}\f$, is *linear* in the magnetic degrees-of-freedom. </li>
!> <li> The algorithm proceeds as follows:
!>
!> **plasma volumes**
!>
!> <ul>
!>
!> <li> In addition to the enclosed toroidal flux, \f$\Delta \psi_t\f$, which is held constant in the plasma volumes,
!>       the Beltrami field in a given volume is assumed to be parameterized by \f$\mu\f$ and \f$\Delta \psi_p\f$.
!>       (Note that \f$\Delta \psi_p\f$ is not defined in a torus.) </li>
!> <li> These are "packed" into an array, e.g. \f$\boldsymbol{\mu} \equiv (\mu, \Delta\psi_p)^T\f$, so that standard library routines ,
!>       e.g. \c C05PCF, can be used to (iteratively) find the appropriately-constrained Beltrami solution, i.e. \f${\bf f}(\boldsymbol{\mu})=0\f$. </li>
!> <li> The function \f${\bf f}(\boldsymbol{\mu})\f$, which is computed by mp00ac(), is defined by the input parameter \c Lconstraint:
!> <ul>
!> <li> If \c Lconstraint = -1, 0, then \f$\boldsymbol{\mu}\f$       is *not* varied and \c Nxdof=0. </li>
!> <li> If \c Lconstraint =  1,    then \f$\boldsymbol{\mu}\f$       is       varied to satisfy the transform constraints;
!>      and \c Nxdof=1 in the simple torus and \c Nxdof=2 in the annular regions.
!>      (Note that in the "simple-torus" region, the enclosed poloidal flux \f$\Delta\psi_p\f$ is not well-defined,
!>      and only \f$\mu=\boldsymbol{\mu}_1\f$ is varied in order to satisfy the transform constraint on the "outer" interface of that volume.) </li>
!> <li> \todo If \c Lconstraint =  2,    then \f$\mu=\boldsymbol{\mu}_1\f$ is       varied in order to satisfy the helicity constraint,
!>      and \f$\Delta\psi_p=\boldsymbol{\mu}_2\f$ is *not* varied, and \c Nxdof=1.
!>      (under re-construction)
!>
!> </li>
!> </ul> </li>
!>
!> </ul>
!>
!> **vacuum volume**
!>
!> <ul>
!>
!> <li> In the vacuum, \f$\mu=0\f$, and the enclosed fluxes, \f$\Delta \psi_t\f$ and \f$\Delta \psi_p\f$, are considered to parameterize the family of solutions.
!>      (These quantities may not be well-defined if \f${\bf B}\cdot{\bf n}\ne 0\f$ on the computational boundary.) </li>
!> <li> These are "packed" into an array, \f$\boldsymbol{\mu} \equiv (\Delta\psi_t,\Delta\psi_p)^T\f$, so that, as above, standard routines can be used
!>      to iteratively find the appropriately constrained solution, i.e. \f${\bf f}(\boldsymbol{\mu})=0\f$. </li>
!> <li> The function \f${\bf f}(\boldsymbol{\mu})\f$, which is computed by mp00ac(), is defined by the input parameter \c Lconstraint:
!> <ul>
!> <li> If \c Lconstraint = -1,   then \f$\boldsymbol{\mu}\f$ is *not* varied and \c Nxdof=0. </li>
!> <li> If \c Lconstraint =  0,2, then \f$\boldsymbol{\mu}\f$ is varied to satisfy the enclosed current constraints, and \c Nxdof=2. </li>
!> <li> If \c Lconstraint =  1,   then \f$\boldsymbol{\mu}\f$ is varied to satisfy
!>      the constraint on the transform on the inner boundary \f$\equiv\f$ plasma boundary and the "linking" current, and \c Nxdof=2.  </li>
!> </ul> </li>
!>
!> </ul>  </li>
!>
!> <li> The Beltrami fields, and the rotational-transform and helicity etc. as required to determine the function \f${\bf f}(\boldsymbol{\mu})\f$
!>      are calculated in mp00ac(). </li>
!> <li> This routine, mp00ac(), is called iteratively if \c Nxdof>1 via
!>      \c C05PCF to determine the appropriately constrained Beltrami field, \f${\bf B}_{\boldsymbol{\mu}}\f$, so that \f${\bf f}(\boldsymbol{\mu})=0\f$. </li>
!> <li> The input variables \c mupftol and \c mupfits control the required accuracy and maximum number of iterations. </li>
!> <li> If \c Nxdof=1, then mp00ac() is called only once to provide the Beltrami fields with the given value of \f$\boldsymbol{\mu}\f$. </li>
!>
!> </ul>

    if (LBlinear) then ! assume Beltrami field is parameterized by helicity multiplier (and poloidal flux);

        lastcpu = MPI_WTIME()

        if (Lplasmaregion) then

            Xdof(1:2) = xoffset + (/mu(lvol), dpflux(lvol)/) ! initial guess for degrees of freedom; offset from zero so that relative error is small;

            select case (Lconstraint)
            case (-1); ; Nxdof = 0 ! multiplier & poloidal flux NOT varied                               ;
                ; ; iflag = 1 !(we don't need derivatives)
            case (0); ; Nxdof = 0 ! multiplier & poloidal flux NOT varied
                ; ; iflag = 1 !(we don't need derivatives)                            ;
            case (1); if (Lcoordinatesingularity) then; Nxdof = 1 ! multiplier                 IS  varied to match       outer transform;
                    ; else; Nxdof = 2 ! multiplier & poloidal flux ARE varied to match inner/outer transform;
                    ; end if
            case (2); Nxdof = 1 ! multiplier                 IS  varied to match             helicity ;
            case (3); if (Lcoordinatesingularity) then; Nxdof = 0 ! multiplier & poloidal flux NOT varied                               ;
                    ; else; Nxdof = 0 ! Global constraint, no dof locally
                    ; end if
                ; ; iflag = 2 !(we still need derivatives)
            end select

        else ! Lvacuumregion ;

            Xdof(1:2) = xoffset + (/dtflux(lvol), dpflux(lvol)/) ! initial guess for degrees of freedom; offset from zero so that relative error is small;

            select case (Lconstraint)
            case (-1); ; Nxdof = 0 ! poloidal   & toroidal flux NOT varied                                                  ;
            case (0); ; Nxdof = 2 ! poloidal   & toroidal flux ARE varied to match linking current and plasma current      ;
            case (1); ; Nxdof = 2 ! poloidal   & toroidal flux ARE varied to match linking current and transform-constraint;
            case (2); ; Nxdof = 2 ! poloidal   & toroidal flux ARE varied to match linking current and plasma current      ;
            case (3); ; Nxdof = 0 ! Fluxes are determined in dforce via a linear system
                ; iflag = 2
            end select

        end if ! end of if( Lplasmaregion) ;

        select case (Nxdof)

        case (0) ! need only call mp00ac once, to calculate Beltrami field for given helicity multiplier and enclosed fluxes;

            ; ; Ndof = 1; Ldfjac = Ndof; nfev = 1; njev = 0; ihybrj = 1; ! provide dummy values for consistency;

            cput = MPI_WTIME()
            Tma02aa = Tma02aa + (cput - cpuo)
            call mp00ac(Ndof, Xdof(1:Ndof), Fdof(1:Ndof), Ddof(1:Ldfjac, 1:Ndof), Ldfjac, iflag)
            cpuo = MPI_WTIME()

            helicity(lvol) = lABintegral(lvol) ! this was computed in mp00ac;

        case (1:2) ! will iteratively call mp00ac, to calculate Beltrami field that satisfies constraints;

            ; ; Ndof = Nxdof; Ldfjac = Ndof; nfev = 0; njev = 0; ihybrj = 0; 
            tol = mupftol; LRR = Ndof*(Ndof + 1)/2; mode = 0; diag(1:2) = zero; factor = one; maxfev = mupfits; nprint = 0

            if (Ndof .gt. 2) then
                write (6, '("ma02aa :      fatal : myid=",i3," ; Ndof.gt.2 ; illegal ;")') myid
                call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                stop "ma02aa : Ndof.gt.2 : illegal  ;"
            end if

            cput = MPI_WTIME()
            Tma02aa = Tma02aa + (cput - cpuo)
            call hybrj2(mp00ac, Ndof, Xdof(1:Ndof), Fdof(1:Ndof), Ddof(1:Ldfjac, 1:Ndof), Ldfjac, tol, &
                        maxfev, diag(1:Ndof), mode, factor, nprint, ihybrj, nfev, njev, RR(1:LRR), LRR, QTF(1:Ndof), &
                        WK(1:Ndof, 1), WK(1:Ndof, 2), WK(1:Ndof, 3), WK(1:Ndof, 4))
            cpuo = MPI_WTIME()

            if (Lplasmaregion) then

                select case (ihybrj)
                case (0:); mu(lvol) = Xdof(1) - xoffset
                    ; ; dpflux(lvol) = Xdof(2) - xoffset
                case (:-1); Xdof(1) = mu(lvol) + xoffset ! mu    and dpflux have been updated in mp00ac; early termination;
                    ; ; Xdof(2) = dpflux(lvol) + xoffset ! mu    and dpflux have been updated in mp00ac; early termination;
                end select

            else ! Lvacuumregion;

                select case (ihybrj)
                case (0:); dtflux(lvol) = Xdof(1) - xoffset
                    ; ; dpflux(lvol) = Xdof(2) - xoffset
                case (:-1); Xdof(1) = dtflux(lvol) + xoffset ! dtflux and dpflux have been updated in mp00ac; early termination;
                    ; ; Xdof(2) = dpflux(lvol) + xoffset ! dtflux and dpflux have been updated in mp00ac; early termination;
                end select

            end if ! end of if( Lplasmaregion ) ;

            if (Lconstraint .ne. 2) helicity(lvol) = lABintegral(lvol) ! this was computed in mp00ac;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

            if (Lconstraint .eq. 1 .or. Lconstraint .eq. 3 .or. (Lvacuumregion .and. Lconstraint .eq. 0)) then

                iflag = 2; Ldfjac = Ndof ! call mp00ac: tr00ab/curent to ensure the derivatives of B, transform, currents, wrt mu/dtflux & dpflux are calculated;

                cput = MPI_WTIME()
                Tma02aa = Tma02aa + (cput - cpuo)
                call mp00ac(Ndof, Xdof(1:Ndof), Fdof(1:Ndof), Ddof(1:Ldfjac, 1:Ndof), Ldfjac, iflag)
                cpuo = MPI_WTIME()

            end if ! end of if( Lconstraint.eq.1 .or. ( Lvacuumregion .and. Lconstraint.eq.0 ) ) ;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

        end select ! end of select case( Nxdof ) ;

        cput = MPI_WTIME()

        select case (ihybrj) ! this screen output may not be correct for Lvacuumregion;
        case (1)
            if (Wma02aa) write (ounit, 1040) cput - cpus, myid, lvol, ihybrj, helicity(lvol), mu(lvol), dpflux(lvol), cput - lastcpu, "success         ", Fdof(1:Ndof)
        case (-2)
            if (Wma02aa) write (ounit, 1040) cput - cpus, myid, lvol, ihybrj, helicity(lvol), mu(lvol), dpflux(lvol), cput - lastcpu, "|F| < mupftol   ", Fdof(1:Ndof)
        case (-1)
            ; write (ounit, 1040) cput - cpus, myid, lvol, ihybrj, helicity(lvol), mu(lvol), dpflux(lvol), cput - lastcpu, "Beltrami fail   ", Fdof(1:Ndof)
        case (0)
            ; write (ounit, 1040) cput - cpus, myid, lvol, ihybrj, helicity(lvol), mu(lvol), dpflux(lvol), cput - lastcpu, "input error     ", Fdof(1:Ndof)
        case (2)
            ; write (ounit, 1040) cput - cpus, myid, lvol, ihybrj, helicity(lvol), mu(lvol), dpflux(lvol), cput - lastcpu, "consider restart", Fdof(1:Ndof)
        case (3)
            ; write (ounit, 1040) cput - cpus, myid, lvol, ihybrj, helicity(lvol), mu(lvol), dpflux(lvol), cput - lastcpu, "xtol too small  ", Fdof(1:Ndof)
        case (4:5)
            ; write (ounit, 1040) cput - cpus, myid, lvol, ihybrj, helicity(lvol), mu(lvol), dpflux(lvol), cput - lastcpu, "bad progress    ", Fdof(1:Ndof)
        case default
            ; write (ounit, 1040) cput - cpus, myid, lvol, ihybrj, helicity(lvol), mu(lvol), dpflux(lvol), cput - lastcpu, "illegal ifail   ", Fdof(1:Ndof)

            if (.true.) then
                write (6, '("ma02aa :      fatal : myid=",i3," ; .true. ; illegal ifail returned by hybrj ;")') myid
                call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                stop "ma02aa : .true. : illegal ifail returned by hybrj  ;"
            end if

        end select

    end if ! end of if( LBlinear ) then;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!> **debugging: finite-difference confirmation of the derivatives of the rotational-transform**
!>
!> <ul>
!>
!> <li> Note that the rotational-transform (if required) is calculated by tr00ab(), which is called by mp00ac(). </li>
!> <li> If \c Lconstraint=1, then mp00ac() will ask tr00ab() to compute the derivatives of the transform
!>      with respect to variations in the helicity-multiplier, \f$\mu\f$, and the enclosed poloidal-flux, \f$\Delta\psi_p\f$, so that
!>      \c C05PCF may more efficiently find the solution. </li>
!> <li> The required derivatives are
!>      \f{eqnarray}{ \frac{\partial {{\,\iota\!\!\!}-}}{\partial \mu}\\
!>                    \frac{\partial {{\,\iota\!\!\!}-}}{\partial \Delta \psi_p}
!>      \f}
!>      to improve the efficiency of the iterative search.
!>      A finite difference estimate of these derivatives is available; need \c DEBUG, \c Lcheck=2 and \c Lconstraint=1. </li>
!> </ul>

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef DEBUG

    if (Lcheck .eq. 2) then ! perform finite difference check on derivative calculation;

        if (Lconstraint .eq. 1 .or. Lvacuumregion) then ! only in this case are the derivatives calculated;

            allocate (dsolution(1:NN, 0:2, -1:1, -1:1), stat=astat)
            dsolution(1:NN, 0:2, -1:1, -1:1) = zero
            ! packed vector potential;

            if (Lplasmaregion) then
                Xdof(1:2) = xoffset + (/mu(lvol), dpflux(lvol)/) ! initial guess for degrees of freedom; offset from zero so that relative error is small;
                if (Lcoordinatesingularity) then; Ndof = 1 ! multiplier                  IS  varied to match       transform      ;
                else; Ndof = 2 ! multiplier &  poloidal flux ARE varied to match inner/outer transform;
                end if
            else ! Lvacuumregion;
                Xdof(1:2) = xoffset + (/dtflux(lvol), dpflux(lvol)/) ! initial guess for degrees of freedom; offset from zero so that relative error is small;
                if (LocalConstraint) then
                    ; ; Ndof = 2
                else
                    ; ; Ndof = 1
                end if
            end if ! end of if( Lplasmaregion) ;

            Ldfjac = Ndof; dFdof(-1:1, -1:1, 1:2) = zero

            ixx = 0; jxx = 0

            iflag = 2 ! iflag controls derivative calculation in mp00ac; analytic derivatives of rotational-transform are required;

            cput = MPI_WTIME()
            Tma02aa = Tma02aa + (cput - cpuo)
            call mp00ac(Ndof, Xdof(1:Ndof), dFdof(ixx, jxx, 1:Ndof), oDdof(1:Ldfjac, 1:Ndof), Ldfjac, iflag)
            cpuo = MPI_WTIME()
            ! compute "exact" derivatives;

            dsolution(1:NN, 1, 0, 0) = solution(1:NN, 1) ! packed vector potential; derivative wrt mu    ;
            dsolution(1:NN, 2, 0, 0) = solution(1:NN, 2) ! packed vector potential; derivative wrt dpflux;

            jfinite = 0
            cput = MPI_WTIME()
            write (ounit, 2000) cput - cpus, myid, lvol, jfinite, "derivative", oDdof(1:Ldfjac, 1:Ndof)

            do jfinite = -4, -2, +1; lfdiff = ten**jfinite

                do ixx = -1, +1 ! centered differences;
                    do jxx = -1, +1 ! centered differences;

                        if (Lcoordinatesingularity .and. jxx .ne. 0) cycle ! Beltrami field in volume with coordinate singularity does not depend on enclosed poloidal flux;

                        if (abs(ixx) + abs(jxx) .ne. 1) cycle ! centered finite differences;

                        if (Lplasmaregion) then
                            Xdof(1:2) = xoffset + (/mu(lvol), dpflux(lvol)/) + (/ixx, jxx/)*lfdiff*half
                        else ! Lvacuumregion;
                            Xdof(1:2) = xoffset + (/dtflux(lvol), dpflux(lvol)/) + (/ixx, jxx/)*lfdiff*half
                        end if ! end of if( Lplasmaregion) ;

                        iflag = 1 ! analytic derivatives of rotational-transform are not required;

                        cput = MPI_WTIME()
                        Tma02aa = Tma02aa + (cput - cpuo)
                        call mp00ac(Ndof, Xdof(1:Ndof), dFdof(ixx, jxx, 1:Ndof), Ddof(1:Ldfjac, 1:Ndof), Ldfjac, iflag)
                        cpuo = MPI_WTIME()
                        ! compute function values only;

                        dsolution(1:NN, 0, ixx, jxx) = solution(1:NN, 0)

                    end do ! end of do jxx;
                end do ! end of do ixx;

                ; Ddof(1:Ndof, 1) = (dFdof(1, 0, 1:Ndof) - dFdof(-1, 0, 1:Ndof))/lfdiff ! derivative wrt helicity multiplier   ;
                if (Ndof .eq. 2) Ddof(1:Ndof, 2) = (dFdof(0, 1, 1:Ndof) - dFdof(0, -1, 1:Ndof))/lfdiff ! derivative wrt enclosed poloidal flux;

                cput = MPI_WTIME()
                !write(ounit,2000) cput-cpus, myid, lvol, jfinite, " error    ", Ddof(1:Ldfjac,1:Ndof) - oDdof(1:Ldfjac,1:Ndof)
                write (ounit, 2000) cput - cpus, myid, lvol, jfinite, " estimate ", Ddof(1:Ldfjac, 1:Ndof)

            end do ! end of do jfinite;

            cput = MPI_WTIME()
            write (ounit, 2000) cput - cpus

            deallocate (dsolution, stat=astat)

2000        format("ma02aa : ":, f10.2, " :":" myid=", i3, " : lvol=", i3, " ; jj=", i3, " ; "a10" : DF=["es23.15" ,"es23.15" ,"es23.15" ,"es23.15" ] ;")

        end if ! end of if( Lconstraint.eq.1 .or. Lvacuumregion ) ;

    end if ! end of if( Lcheck.eq.2 ) ;

#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

9999 continue
    cput = MPI_WTIME()
    Tma02aa = Tma02aa + (cput - cpuo)
    return

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

1010 format("ma02aa : ", f10.2, " : myid=", i3, " ; lvol=", i3, " ; SQP    : ie04uff=", i3, " hel="es12.4" mu="es12.4" dpflux="es12.4" time="f9.1" ":, a36)
1020 format("ma02aa : ", f10.2, " : myid=", i3, " ; lvol=", i3, " ; Newton : ihybrj1=", i3, " hel="es12.4" mu="es12.4" dpflux="es12.4" time="f9.1" ; " &
           "error="es7.0" ; ":, a18)
1040 format("ma02aa : ", f10.2, " : myid=", i3, " ; lvol=", i3, " ; Linear : ihybrj =", i3, " hel="es12.4" mu="es12.4" dpflux="es12.4" time="f9.1" ; " &
           :, a16" ; F="2es08.0)
!050 format("ma02aa : ",f10.2," : myid=",i3," ; lvol=",i3," ; Linear : ihybrj =",i3,"     "  12x " I ="es12.4"        "  12x " time="f9.1" ; "&
! :,a16" ; F="2es08.0)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine ma02aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
