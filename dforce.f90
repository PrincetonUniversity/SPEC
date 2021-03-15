!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (&ldquo;global&rdquo; force) ! Given &ldquo;position&rdquo;, ${\bf \xi}$, computes ${\bf F}({\bf \xi})$ and $\nabla_{\bf \xi}{\bf F}$.

!latex \briefly{Calculates ${\bf F}({\bf x})$, where ${\bf x} \equiv \{\mbox{\rm geometry}\} \equiv \{ R_{i,v}, Z_{i,v}\}$ 
!latex          and ${\bf F}\equiv[[p+B^2/2]] + \{\mbox{\rm spectral constraints}\} $, and $\nabla {\bf F}$.}

!latex \calledby{\link{hesian}, 
!latex           \link{newton}, 
!latex           \link{pc00aa}, 
!latex           \link{pc00ab} and 
!latex           \link{xspech}} \\

!latex \calls{\link{packxi}, 
!latex        \link{ma00aa}, 
!latex        \link{matrix}, 
!latex        \link{dfp100}, 
!latex        \link{dfp200} and
!latex        \link{brcast}}


!latex \tableofcontents

!latex \subsection{unpacking}

!latex \begin{enumerate}

!latex \item The geometrical degrees of freedom are represented as a vector, ${\bf x} \equiv \{ R_{i,v}, Z_{i,v}\}$, 
!latex       where $i=1,$ \internal{mn} labels the Fourier harmonic and $v=1,$ \internal{Mvol}$-1$ is the interface label.
!latex       This vector is ``unpacked'' using \link{packxi}.
!latex       (Note that \link{packxi} also sets the coordinate axis, i.e. the $R_{i,0}$ and $Z_{i,0}$.)

!latex \end{enumerate}

!latex \subsection{Matrices computation}

!latex \begin{enumerate}
!latex       \item     the volume-integrated metric arrays, \internal{DToocc}, etc. are evaluated in each volume by calling \link{ma00aa};
!latex       \item     the energy and helicity matrices, \internal{dMA(0:NN,0:NN)}, \internal{dMB(0:NN,0:2)}, etc. are evaluated in each 
!latex                 volume by calling \link{matrix};
!latex \end{enumerate}

!latex \subsection{parallelization over volumes}

!latex Two different cases emerge: either a local constraint or a global constraint is considered. This condition is determined by the 
!latex flag \inputvar{LocalConstraint}.

!latex \subsubsection{Local constraint}
!latex In each volume, \internal{vvol = 1, Mvol}, 
!latex       \begin{enumerate}
!latex       \item The logical array \internal{ImagneticOK(vvol)} is set to \internal{.false.}
!latex          \item The MPI node associated to the volume calls \link{dfp100}. This routine calls \link{ma02aa} (and might iterate on \link{mp00ac}) and computes the
!latex                field solution in each volume consistent with the constraint.
!latex         \item The MPI node associated to the volume calls \link{dfp200}. This computes $p+B^2/2$ (and the spectral constraints if required) at the interfaces in 
!latex                each volumes, as well as the derivatives of the force-balance if \internal{LComputeDerivatives = 1};
!latex       \end{enumerate}

!latex \subsubsection{Global constraint}
!latex The MPI node $0$ minimizes the constraint with HYBRID1 by iterating on \link{dfp100} until the field matches the constraint. Other MPI nodes enter
!latex the subroutine loop\_dfp100. In loop\_dfp100, each MPI node
!latex \begin{enumerate}
!latex \item calls \link{dfp100}
!latex \item solves the field in its associated volumes
!latex \item communicates the field to the node $0$
!latex \item repeats this loop until the node $0$ sends a flag \internal{iflag=5}.
!latex \end{enumerate}

!latex \subsection{broadcasting}

!latex \begin{enumerate}
!latex \item The required quantities are broadcast by \link{brcast}.
!latex \end{enumerate}

!latex \subsection{construction of force}

!latex \begin{enumerate}
!latex \item The force vector, ${\bf F}({\bf x})$, is a combination of the pressure-imbalance Fourier harmonics, $[[p+B^2/2]]_{i,v}$,
!latex       where $i$ labels Fourier harmonic and $v$ is the interface label:
!latex       \be F_{i,v} \equiv \left[ ( p_{v+1}+B^2_{i,v+1}/2 ) - ( p_v + B^2_{i,v}/2 ) \right] \times \internal{BBweight}_i,
!latex       \ee
!latex       where \internal{BBweight\_i} is defined in \link{preset};
!latex       and the spectral condensation constraints, 
!latex       \be F_{i,v} \equiv I_{i,v} \times \inputvar{epsilon} + S_{i,v,1} \times \internal{sweight}_v - S_{i,v+1,0} \times \internal{sweight}_{v+1},
!latex       \ee
!latex       where the spectral condensation constraints, $I_{i,v}$, and the ``star-like'' poloidal angle constraints, $S_{i,v,\pm 1}$,
!latex       are calculated and defined in \link{lforce};
!latex       and the \internal{sweight}$_v$ are defined in \link{preset}. All quantities local to a volume are computed in \link{dfp200},
!latex          information is then broadcasted to the MPI node $0$ in \link{dforce} and the global force is evaluated.
!latex \end{enumerate}


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef DEBUG
recursive subroutine dforce( NGdof, position, force, LComputeDerivatives, LComputeAxis)
#else
          subroutine dforce( NGdof, position, force, LComputeDerivatives, LComputeAxis)
#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half, one, pi, pi2
  
  use numerical, only : logtolerance
  
  use fileunits, only : ounit
  
  use inputlist, only : Wmacros, Wdforce, ext, Nvol, Ntor, Lrad, Igeometry, &
                        epsilon, &
                        Lconstraint, Lcheck, dRZ, &
                        Lextrap, &
                        mupftol, &
                        Lfreebound, &
                        ext ! For outputing Lcheck = 6 test
  
  use cputiming, only : Tdforce
  
  use allglobal, only : ncpu, myid, cpus, &
                        Mvol, NAdof, &
                        Iquad, &                  ! convenience; provided to ma00aa as argument to avoid allocations;
                        iRbc, iZbs, iRbs, iZbc, & ! Fourier harmonics of geometry; vector of independent variables, position, is "unpacked" into iRbc,iZbs;
                        ImagneticOK, &
                        Energy, ForceErr, &
                        YESstellsym, NOTstellsym, &
                        Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, &
                        mn, im, in, &
                        dpflux, dtflux, sweight, &
                        Bemn, Bomn, Iomn, Iemn, Somn, Semn, &
                        BBe, IIo, BBo, IIe, & ! these are just used for screen diagnostics;
                        LGdof, dBdX, &
                        Ate, Aze, Ato, Azo, & ! only required for broadcasting
                        diotadxup, dItGpdxtp, & ! only required for broadcasting
                        lBBintegral, &
                        dFFdRZ, dBBdmp, dmupfdx, hessian, dessian, Lhessianallocated, &
                        BBweight, & ! exponential weight on force-imbalance harmonics;
                        psifactor, &
                        LocalConstraint, xoffset, &
                        solution, IPdtdPf, &
                        IsMyVolume, IsMyVolumeValue, WhichCpuID
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, parameter   :: NB = 3 ! optimal workspace block size for LAPACK:DSYSVX;

  INTEGER, intent(in)  :: NGdof               ! dimensions;
  REAL,    intent(in)  :: position(0:NGdof)   ! degrees-of-freedom = internal geometry;
  REAL,    intent(out) :: force(0:NGdof)      ! force;
  LOGICAL, intent(in)  :: LComputeDerivatives ! indicates whether derivatives are to be calculated;
  
  INTEGER              :: vvol, innout, ii, jj, irz, issym, iocons, tdoc, idoc, idof, tdof, jdof, ivol, imn, ll, ihybrd1, lwa, Ndofgl, llmodnp
  INTEGER              :: maxfev, ml, muhybr, mode, nprint, nfev, ldfjac, lr, Nbc, NN, cpu_id, ideriv
  REAL                 :: epsfcn, factor
  REAL                 :: Fdof(1:Mvol-1), Xdof(1:Mvol-1)
  INTEGER              :: ipiv(1:Mvol)
  REAL, allocatable    :: fjac(:, :), r(:), Fvec(:), dpfluxout(:)

  INTEGER              :: status(MPI_STATUS_SIZE), request_recv, request_send, cpu_send
  INTEGER              :: id
  INTEGER              :: idgesv, Lwork

  CHARACTER            :: packorunpack 
  EXTERNAL             :: dfp100, dfp200

  LOGICAL              :: LComputeAxis, dfp100_logical

#ifdef DEBUG
  INTEGER              :: isymdiff
  REAL                 :: dvol(-1:+1), evolume, imupf(1:2,-2:2), lfactor
  REAL,    allocatable :: isolution(:,:)
  REAL,   allocatable :: oRbc(:,:), oZbs(:,:), oRbs(:,:), oZbc(:,:), iforce(:,:), iposition(:,:), finitediff_hessian(:,:) ! original geometry;
#endif

  BEGIN(dforce)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! Unpack position to generate arrays iRbc, iZbs, IRbs, iZbc.

  packorunpack = 'U' ! unpack geometrical degrees-of-freedom;

#ifndef DEBUG  
  LComputeAxis = .true.
#endif

  WCALL( dforce, packxi,( NGdof, position(0:NGdof), Mvol, mn, iRbc(1:mn,0:Mvol), iZbs(1:mn,0:Mvol), &
                          iRbs(1:mn,0:Mvol), iZbc(1:mn,0:Mvol), packorunpack, LcomputeDerivatives, LComputeAxis ) )
 
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( LcomputeDerivatives ) then
#ifdef DEBUG
   FATAL( dforce, .not.allocated(dBBdmp), do not pass go )
#endif
   dBBdmp(1:LGdof,1:Mvol,0:1,1:2) = zero
  endif


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! SOLVE FIELD IN AGREEMENT WITH CONSTRAINTS AND GEOMETRY
! ------------------------------------------------------
! Two different cases, both with their own parallelization. 
! If LOCAL constraint, each process iterates on the local poloidal flux and Lagrange multiplier
! to match the local constraint. Then all information is broadcasted by the master thread
! If GLOBAL constraint, only the master thread iterates on all Lagrange multipliers and poloidal
! fluxes. Over threads are stuck in an infinite loop where they "help" the master thread compute
! each iteration. At the last iteration, master thread send IconstraintOK=.TRUE. to all threads, and
! they exit their infinite loops.


! Local constraint case - simply call dfp100 and then dfp200
  Xdof(1:Mvol-1) = zero;

  if( LocalConstraint ) then

    SALLOCATE( Fvec, (1:Mvol-1), zero)

    Ndofgl = 0; Fvec(1:Mvol-1) = 0; dfp100_logical = .FALSE.;
    Xdof(1:Mvol-1) = dpflux(2:Mvol) + xoffset
    
    ! Solve for field
	  dBdX%L = LComputeDerivatives
    WCALL(dforce, dfp100, (Ndofgl, Xdof, Fvec, dfp100_logical) )

    DALLOCATE( Fvec )

    !do vvol=1,Mvol
    !  WCALL(dforce, brcast, (vvol) )
    !enddo
! --------------------------------------------------------------------------------------------------
! Global constraint - call the master thread calls hybrd1 on dfp100, others call dfp100_loop.
  else

    IPDtdPf = zero
    Xdof(1:Mvol-1)   = dpflux(2:Mvol) + xoffset

    if( Lfreebound.eq.1 ) then
      ! Mvol-1 surface current plus 1 poloidal linking current constraints
      Ndofgl = Mvol
    else
      ! Mvol-1 surface current constraints
      Ndofgl = Mvol-1
    endif 

    SALLOCATE( Fvec, (1:Ndofgl), zero )

    dfp100_logical = .TRUE.

    WCALL(dforce, dfp100, (Ndofgl, Xdof(1:Mvol-1), Fvec(1:Ndofgl), dfp100_logical))

    SALLOCATE(dpfluxout, (1:Ndofgl), zero )
    if ( myid .eq. 0 ) then 

        dpfluxout = Fvec
        call DGESV( Ndofgl, 1, IPdtdPf, Ndofgl, ipiv, dpfluxout, Ndofgl, idgesv )

        ! one step Newton's method
        dpflux(2:Mvol) = dpflux(2:Mvol) - dpfluxout(1:Mvol-1)
        if( Lfreebound.eq.1 ) then
          dtflux(Mvol) = dtflux(Mvol  ) - dpfluxout(Mvol    )
        endif
    endif

    ! Broadcast the field and pflux
    RlBCAST(dpfluxout(1:Ndofgl), Ndofgl, 0)
    RlBCAST(dpflux(1:Mvol)   , Mvol, 0)
    if( Lfreebound.eq.1 ) then
      RlBCAST(dtflux(Mvol), 1, 0)
    endif

    do vvol = 2, Mvol
  
      WCALL(dforce, IsMyVolume, (vvol))

      if( IsMyVolumeValue .EQ. 0 ) then
          cycle
      else if( IsMyVolumeValue .EQ. -1) then
          FATAL(dforce, .true., Unassociated volume)
      endif

      NN = NAdof(vvol)

      SALLOCATE( solution, (1:NN, 0:2), zero)

      ! Pack field and its derivatives
      packorunpack = 'P'
      WCALL( dforce, packab, ( packorunpack, vvol, NN, solution(1:NN,0), 0 ) ) ! packing;
      WCALL( dforce, packab, ( packorunpack, vvol, NN, solution(1:NN,2), 2 ) ) ! packing;

      ! compute the field with renewed dpflux via single Newton method step
      if( Lfreebound.eq.1 .and.(vvol.eq.Mvol) ) then
        WCALL( dforce, packab, ( packorunpack, vvol, NN, solution(1:NN,1), 1 ) ) ! packing;
        solution(1:NN, 0) = solution(1:NN, 0) - dpfluxout(vvol-1) * solution(1:NN, 2) & ! derivative w.r.t pflux
                                              - dpfluxout(vvol  ) * solution(1:NN, 1)   ! derivative w.r.t tflux
      else
        solution(1:NN, 0) = solution(1:NN, 0) - dpfluxout(vvol-1) * solution(1:NN, 2)
      endif
      
      ! Unpack field in vector potential Fourier harmonics
      packorunpack = 'U'
      WCALL( dforce, packab, ( packorunpack, vvol, NN, solution(1:NN,0), 0 ) ) ! unpacking;                              

      DALLOCATE( solution )

    enddo ! end of do vvol = 1, Mvol



    DALLOCATE(Fvec)
    DALLOCATE(dpfluxout)



! #ifdef DEBUG
!       select case( ihybrd1 )
!         case( 1   )  ; write(ounit,'("dforce : ",f10.2," : finished ; success        ; dpflux = ", es12.5, ", its="i7";")') cput-cpus, dpflux, nfev
!         case( 0   )  ; write(ounit,'("dforce : ",f10.2," : finished ; input error    ; dpflux = ", es12.5, ", its="i7";")') cput-cpus, dpflux, nfev
!         case( 2   )  ; write(ounit,'("dforce : ",f10.2," : finished ; max. iter      ; dpflux = ", es12.5, ", its="i7";")') cput-cpus, dpflux, nfev
!         case( 3   )  ; write(ounit,'("dforce : ",f10.2," : finished ; xtol too small ; dpflux = ", es12.5, ", its="i7";")') cput-cpus, dpflux, nfev
!         case( 4:5 )  ; write(ounit,'("dforce : ",f10.2," : finished ; bad progress   ; dpflux = ", es12.5, ", its="i7";")') cput-cpus, dpflux, nfev
!         case default ; write(ounit,'("dforce : ",f10.2," : finished ; illegal ifail  ; dpflux = ", es12.5, ", its="i7";")') cput-cpus, dpflux, nfev
!       end select
! #endif

  endif !matches if( LocalConstraint ) 

! --------------------------------------------------------------------------------------------------
!                                    MPI COMMUNICATIONS

! Finally broadcast the field information to all threads from the thread which did the computation
! TODO: improve MPI communication
  do vvol = 1, Mvol
    call WhichCpuID(vvol, cpu_id)

    ! Broadcast all ImagneticOK
    !write(ounit,'("dforce : " 10x " : myid="i3"; vvol="i3"; ; ImagneticOK="999L2)') myid, vvol, ImagneticOK(1:Mvol)
    !write(ounit,'("dforce : " 10x " : cpu_id="i3"; vvol="i3"; ; ImagneticOK="999L2)') cpu_id, vvol, ImagneticOK(vvol)
    LlBCAST( ImagneticOK(vvol)         , 1, cpu_id)
  
    do ideriv=0,2
      if( (.not.LcomputeDerivatives) .and. (ideriv.ne.0) ) cycle
      do ii = 1, mn
        RlBCAST( Ate(vvol,ideriv,ii)%s(0:Lrad(vvol)), Lrad(vvol)+1, cpu_id)
        RlBCAST( Aze(vvol,ideriv,ii)%s(0:Lrad(vvol)), Lrad(vvol)+1, cpu_id)
      enddo
    enddo


    if( NOTstellsym ) then
      do ideriv=0,2
      if( (.not.LcomputeDerivatives) .and. (ideriv.ne.0) ) cycle
        do ii = 1, mn
              RlBCAST( Ato(vvol,ideriv,ii)%s(0:Lrad(vvol)), Lrad(vvol)+1, cpu_id)
              RlBCAST( Azo(vvol,ideriv,ii)%s(0:Lrad(vvol)), Lrad(vvol)+1, cpu_id)
        enddo
      enddo
    endif
  enddo

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  ! Compute local force and derivatives 
  WCALL(dforce, dfp200, ( LcomputeDerivatives, vvol) )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  if( Lcheck.eq.2 ) then
   write(ounit,'("dforce : ", 10x ," : myid=",i3," ; finished computing derivatives of rotational-transform wrt mu and dpflux ;")') myid
   stop "dforce :            : myid=    ; finished computing derivatives of rotational-transform wrt mu and dpflux ;" ! this will allow other cpus to finish;
  endif

  FATAL( dforce, Lcheck.eq.2, finished computing derivatives of rotational-transform wrt mu and dpflux )

  if( Wdforce ) write(ounit,'("dforce : " 10x " : myid="i3" ; LComputeDerivatives="L2" ; ImagneticOK="999L2)') myid, LComputeDerivatives, ImagneticOK(1:Mvol)
#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! Broadcast information to all CPUs
  do vvol = 1, Mvol

    LREGION( vvol )
    WCALL( dforce, brcast, ( vvol ) )

  enddo
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef DEBUG
  if( Wdforce ) write(ounit,'("dforce : " 10x " : myid="i3" ; LComputeDerivatives="L2" ; ImagneticOK="999L2)') myid, LComputeDerivatives, ImagneticOK(1:Mvol)
#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  lBBintegral(1:Nvol) = lBBintegral(1:Nvol) * half
  
  Energy = sum( lBBintegral(1:Nvol) ) ! should also compute beta;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


! CONSTRUCT FORCE
! ---------------
  
  ;   force(0:NGdof) = zero
  
  do vvol = 1, Mvol-1

    LREGION(vvol)
   
    tdoc = (vvol-1) * LGdof 
   
    if( ImagneticOK(vvol) .and. ImagneticOK(vvol+1) ) then ! the magnetic fields in the volumes adjacent to this interface are valid;
    
      ;  idoc = 0           ! degree-of-constraint counter; set;
      
      if( Lextrap.eq.1 .and. vvol.eq.1 ) then ! to be made redundant;
        FATAL( dforce, 2.gt.Mvol, psifactor needs attention )
        ;force(tdoc+idoc+1:tdoc+idoc+mn) = position(1:mn) - ( iRbc(1:mn,2) / psifactor(1:mn,2) )
      else
        ;force(tdoc+idoc+1:tdoc+idoc+mn    ) = ( Bemn(1:mn    ,vvol+1,0) - Bemn(1:mn    ,vvol+0,1) ) * BBweight(1:mn) ! pressure imbalance;
      endif
      
      ;  BBe(vvol) = max( sum( abs( force(tdoc+idoc+1:tdoc+idoc+mn  ) ) ) / (mn  ), logtolerance ) ! screen diagnostics;
      
      ;  idoc = idoc + mn   ! degree-of-constraint counter; increment;
      
      if( Igeometry.ge.3 ) then ! add spectral constraints;
      
        force(tdoc+idoc+1:tdoc+idoc+mn-1  ) = (                           Iomn(2:mn    ,vvol+0  ) ) * epsilon         & ! spectral constraints;
                                            + (                         + Somn(2:mn    ,vvol+0,1) ) * sweight(vvol+0) & ! poloidal length constraint;
                                            - ( Somn(2:mn    ,vvol+1,0)                           ) * sweight(vvol+1)
            
  !     if( Ntor.gt.0 ) then ! poloidal angle origin is not otherwise constrained ;
  !      force(tdoc+idoc+1:tdoc+idoc+Ntor  ) = ( Pomn(2:Ntor+1,vvol+1,0) - Pomn(2:Ntor+1,vvol+0,1) ) * apsilon ! choice of spectral constraint can be enforced;
  !     endif
      
        IIo(vvol) = max( sum( abs( force(tdoc+idoc+1:tdoc+idoc+mn-1) ) ) / (mn-1), logtolerance ) ! screen diagnostics;
      
        idoc = idoc + mn-1
      
      endif ! end of if( Igeometry.ge.3 ) ;
      
      if( NOTstellsym ) then
      
        force(tdoc+idoc+1:tdoc+idoc+mn-1  ) = ( Bomn(2:mn    ,vvol+1,0) - Bomn(2:mn    ,vvol+0,1) ) * BBweight(2:mn) ! pressure imbalance;
      
        BBo(vvol) = max( sum( abs( force(tdoc+idoc+1:tdoc+idoc+mn-1) ) ) / (mn-1), logtolerance ) ! screen diagnostics;
      
        idoc = idoc + mn-1 ! degree-of-constraint counter; increment;
      
        if( Igeometry.ge.3 ) then ! add spectral constraints;
        
          force(tdoc+idoc+1:tdoc+idoc+mn    ) = (                           Iemn(1:mn    ,vvol+0  ) ) * epsilon         & ! spectral constraints;
                                              + (                         + Semn(1:mn    ,vvol+0,1) ) * sweight(vvol+0) & ! poloidal length constraint;
                                              - ( Semn(1:mn    ,vvol+1,0)                           ) * sweight(vvol+1)
        
  !     if( Ntor.ge.0 ) then
  !      force(tdoc+idoc+1:tdoc+idoc+Ntor+1) = ( Pemn(1:Ntor+1,vvol+1,0) - Pemn(1:Ntor+1,vvol+0,1) ) * apsilon ! choice of spectral constraint can be enforced;
  !     endif
        
          IIe(vvol) = max( sum( abs( force(tdoc+idoc+1:tdoc+idoc+mn  ) ) ) / (mn  ), logtolerance ) ! screen diagnostics;
        
          idoc = idoc + mn   ! degree-of-constraint counter; increment;
        
        endif ! end of if( Igeometry.ge.3 ) ;
      
      endif ! end of if( NOTstellsym ) ;
      
#ifdef DEBUG
      FATAL( dforce, idoc.ne.LGdof, counting error ) ! this has caught bugs;
#endif
    
    else ! matches if( ImagneticOK(vvol) .and. ImagneticOK(vvol+1) );
    
      ;                       ; BBe(vvol) = 9.9E+09
      ;                       ; IIo(vvol) = 9.9E+09
      if ( NOTstellsym ) then ; BBo(vvol) = 9.9E+09
      ;                      ; IIe(vvol) = 9.9E+09
      endif
      
      ; force(tdoc+1:tdoc+LGdof) = 9.9E+09
    
    endif ! end of if( ImagneticOK(vvol) .and. ImagneticOK(vvol+1) ) ;
   
  enddo ! end of do vvol;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( NGdof.ne.0 ) then ; ForceErr = sqrt( sum( force(1:NGdof)*force(1:NGdof) ) / NGdof ) ! this includes spectral constraints;
  else                  ; ForceErr = zero
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  
  if( Wdforce .and. myid.eq.0 ) then
   
    cput = GETTIME
    ;                   ; write(ounit,4000) cput-cpus, ForceErr, cput-cpuo, "|BB|e", alog10(BBe(1:min(Mvol-1,28)))
    if( Igeometry.ge.3 )  write(ounit,4001)                                 "|II|o", alog10(IIo(1:min(Mvol-1,28)))
    if( NOTstellsym ) then
      ;                  ; write(ounit,4001)                                 "|BB|o", alog10(BBo(1:min(Mvol-1,28)))
      if( Igeometry.ge.3 ) write(ounit,4001)                                 "|II|e", alog10(IIe(1:min(Mvol-1,28)))
    endif
   
  endif ! end of if( Wdforce .and. myid.eq.0 ) ;
  
#endif

4000 format("dforce : ",f10.2," : ",6x,3x,"; ",:,"|f|=",es12.5," ; ",:,"time=",f10.2,"s ;",:," log",a5,"=",28f6.2  ," ...")
4001 format("dforce : ", 10x ," : ",6x,3x,"; ",:,"    ",  12x ,"   ",:,"     ", 10x ,"  ;",:," log",a5,"=",28f6.2  ," ...")

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! CONSTRUCT HESSIAN
! -----------------

  if( LcomputeDerivatives ) then ! construct Hessian;

#ifdef DEBUG
    FATAL( dforce, .not.Lhessianallocated, need to allocate hessian )
#endif

    hessian(1:NGdof,1:NGdof) = zero 

#ifdef DEBUG
    if( Lcheck.eq.6 ) then
      SALLOCATE( finitediff_hessian, (1:NGdof, 1:NGdof), zero )
    endif
#endif

    do vvol = 1, Mvol-1 ! loop over interior surfaces;

      if( ImagneticOK(vvol) .and. ImagneticOK(vvol+1) ) then ! the magnetic fields in the volumes adjacent to this interface are valid;

        idof = 0 ! labels degree-of-freedom = Fourier harmonic of surface geometry;

#ifdef DEBUG
        if( idof.gt.LGdof ) write(ounit,1000) myid, vvol, -1, -1, -1, idof, LGdof ! can be deleted;
#endif

        do ii = 1, mn ! loop over degrees-of-freedom;

#ifdef DEBUG
          if( idof.gt.LGdof ) write(ounit,1000) myid, vvol, ii, -1, -1, idof, LGdof ! can be deleted;
#endif

          do irz = 0, 1 ! Fourier harmonic of R, Fourier harmonic of Z;

#ifdef DEBUG
            if( idof.gt.LGdof ) write(ounit,1000) myid, vvol, ii, irz, -1, idof, LGdof ! can be deleted;
#endif

            if( irz.eq.1 .and. Igeometry.lt.3 ) cycle ! no dependence on Z;

            do issym = 0, 1 ! stellarator symmetry;

#ifdef DEBUG
              if( idof.gt.LGdof ) write(ounit,1000) myid, vvol, ii, irz, issym, idof, LGdof ! can be deleted;
#endif

              if( issym.eq.1 .and. YESstellsym ) cycle ! no dependence on the non-stellarator symmetric harmonics;

              if( ii.eq.1 .and. irz.eq.1 .and. issym.eq.0 ) cycle ! no dependence on Zbs_{m=0,n=0};
              if( ii.eq.1 .and. irz.eq.0 .and. issym.eq.1 ) cycle ! no dependence on Rbs_{m=0,n=0};
      
              idof = idof + 1 ! labels degree-of-freedom;

#ifdef DEBUG
              if( idof.gt.LGdof ) write(ounit,1000) myid, vvol, ii, irz, issym, idof, LGdof ! can be deleted;
1000 format("hforce : " 10x " : myid=",i3," ; vvol=",i3," ; ii= ",i3," ; irz="i3" ; issym="i3" ; idof="i3" ; LGdof="i3" ;")
                FATAL( hforce, idof.gt.LGdof, illegal degree-of-freedom index constructing hessian ) ! can be deleted;
#endif

              if( LocalConstraint ) then
                ! Derivative with respect to previous interface
                if( vvol.gt.1 ) then
                  tdof = (vvol-2) * LGdof + idof ! labels degree-of-freedom in internal interface geometry   ;
                  tdoc = (vvol-1) * LGdof        ! labels force-balance constraint across internal interfaces;
                  idoc = 0                       ! local  force-balance constraint across internal interface ;
                  hessian(tdoc+idoc+1:tdoc+idoc+LGdof,tdof) =                                           - dFFdRZ(idoc+1:idoc+LGdof,1,idof,0,vvol+0)
                  if( Lconstraint.eq.1 ) then ! this is a little clumsy; could include Lfreebound or something . . . ;
                    hessian(tdoc+idoc+1:tdoc+idoc+LGdof,tdof) =  hessian(tdoc+idoc+1:tdoc+idoc+LGdof,tdof)                     &
                                                                - dBBdmp(idoc+1:idoc+LGdof,vvol+0,1,1) * dmupfdx(vvol,1,1,idof,0) &
                                                                - dBBdmp(idoc+1:idoc+LGdof,vvol+0,1,2) * dmupfdx(vvol,1,2,idof,0)
                  endif ! end of if( Lconstraint.eq.1 ) ; 
                endif ! end of if( vvol.gt.1 ) ;

                ! Derivative with respect to current interface
                ;tdof = (vvol-1) * LGdof + idof
                ;tdoc = (vvol-1) * LGdof ! shorthand;
                ;idoc = 0

                if( Lextrap.eq.1 .and. vvol.eq.1 ) then
                    ;hessian(tdoc+idof                  ,tdof) = one ! diagonal elements;
                else
                    ;hessian(tdoc+idoc+1:tdoc+idoc+LGdof,tdof) = dFFdRZ(idoc+1:idoc+LGdof,0,idof,0,vvol+1) - dFFdRZ(idoc+1:idoc+LGdof,1,idof,1,vvol+0)
                    if( Lconstraint.eq.1 ) then ! this is a little clumsy;
                        hessian(tdoc+idoc+1:tdoc+idoc+LGdof,tdof) =  hessian(tdoc+idoc+1:tdoc+idoc+LGdof,tdof)                       &
                                                                    + dBBdmp(idoc+1:idoc+LGdof,vvol+1,0,1) * dmupfdx(vvol+1,1,1,idof,0) &
                                                                    + dBBdmp(idoc+1:idoc+LGdof,vvol+1,0,2) * dmupfdx(vvol+1,1,2,idof,0) &
                                                                    - dBBdmp(idoc+1:idoc+LGdof,vvol+0,1,1) * dmupfdx(vvol+0,1,1,idof,1) &
                                                                    - dBBdmp(idoc+1:idoc+LGdof,vvol+0,1,2) * dmupfdx(vvol+0,1,2,idof,1)
                    endif ! end of if( Lconstraint.eq.1 );
                endif ! end of if( Lextrap.eq.1 .and. vvol.eq.1 )


                ! Derivative with respect to next interface
                if( vvol.lt.Mvol-1 ) then
                  tdof = (vvol+0) * LGdof + idof
                  tdoc = (vvol-1) * LGdof ! shorthand;
                  idoc = 0
                  if( Lextrap.eq.1 .and. vvol.eq.1 ) then
                    if    ( im(idof).le.0                     ) then ; hessian(tdoc+idof,tdof) = - one
                    else                                             ; hessian(tdoc+idof,tdof) = - one
                    endif
                  else
                    hessian(tdoc+idoc+1:tdoc+idoc+LGdof,tdof) = dFFdRZ(idoc+1:idoc+LGdof,0,idof,1,vvol+1)
                    if( Lconstraint.eq.1 ) then ! this is a little clumsy;
                      hessian(tdoc+idoc+1:tdoc+idoc+LGdof,tdof) =  hessian(tdoc+idoc+1:tdoc+idoc+LGdof,tdof)                       &
                                                                  + dBBdmp(idoc+1:idoc+LGdof,vvol+1,0,1) * dmupfdx(vvol+1,1,1,idof,1) &
                                                                  + dBBdmp(idoc+1:idoc+LGdof,vvol+1,0,2) * dmupfdx(vvol+1,1,2,idof,1)
                    endif ! end of if( Lconstraint.eq.1 ) then;
                  endif
                endif ! end of if( vvol.lt.Mvol-1 ) ;


                ! Case of last interface in case of boundary variation ?
                if( vvol.eq.Mvol-1 ) then
                  !tdof = (vvol+0) * LGdof + idof
                  tdoc = (vvol-1) * LGdof ! shorthand ;
                  idoc = 0
                  dessian(tdoc+idoc+1:tdoc+idoc+LGdof,idof) = dFFdRZ(idoc+1:idoc+LGdof,0,idof,1,vvol+1)
                  if( Lconstraint.eq.1 ) then ! this is a little clumsy;
                    dessian(tdoc+idoc+1:tdoc+idoc+LGdof,idof) =  dessian(tdoc+idoc+1:tdoc+idoc+LGdof,idof)                       &
                                                                + dBBdmp(idoc+1:idoc+LGdof,vvol+1,0,1) * dmupfdx(vvol+1,1,1,idof,1) &
                                                                + dBBdmp(idoc+1:idoc+LGdof,vvol+1,0,2) * dmupfdx(vvol+1,1,2,idof,1)
                  endif ! end of if( Lconstraint.eq.1 ) then;

                endif ! end of if( vvol.lt.Mvol-1 ) ;
                    
              else ! Global constraint
              
                ! In the general case of global constraint, there are no zero element in the hessian. We thus loop again on all volumes

                do ivol = 1, Mvol-1
                  tdoc = (ivol-1) * LGdof ! shorthand ;
                  tdof = (vvol-1) * LGdof + idof

                  if( ivol.eq.vvol-1 ) then 
                    hessian(tdoc+1:tdoc+LGdof,tdof) =  dFFdRZ(1:LGdof,0,idof,1,ivol+1)
                  elseif( ivol.eq.vvol ) then 
                    hessian(tdoc+1:tdoc+LGdof,tdof) =  dFFdRZ(1:LGdof,0,idof,0,ivol+1) - dFFdRZ(1:LGdof,1,idof,1,ivol)
                  elseif( ivol.eq.vvol+1 ) then
                    hessian(tdoc+1:tdoc+LGdof,tdof) =                                  - dFFdRZ(1:LGdof,1,idof,0,ivol)
                  endif


                  hessian(tdoc+1:tdoc+LGdof,tdof) = hessian(tdoc+1:tdoc+LGdof,tdof)                              &
                                                    + dBBdmp(1:LGdof,ivol+1,0,1) * dmupfdx(ivol+1,vvol,1,idof,1) &
                                                    + dBBdmp(1:LGdof,ivol+1,0,2) * dmupfdx(ivol+1,vvol,2,idof,1) &
                                                    - dBBdmp(1:LGdof,ivol+0,1,1) * dmupfdx(ivol+0,vvol,1,idof,1) &
                                                    - dBBdmp(1:LGdof,ivol+0,1,2) * dmupfdx(ivol+0,vvol,2,idof,1)

                enddo



              endif ! matches if( LocalConstraint );

#ifdef DEBUG
              if( Lcheck.eq.6 ) then
                dBdX%L = .false.
                SALLOCATE( oRbc, (1:mn,0:Mvol), iRbc(1:mn,0:Mvol) ) !save unperturbed geometry
                SALLOCATE( oZbs, (1:mn,0:Mvol), iZbs(1:mn,0:Mvol) )
                SALLOCATE( oRbs, (1:mn,0:Mvol), iRbs(1:mn,0:Mvol) )
                SALLOCATE( oZbc, (1:mn,0:Mvol), iZbc(1:mn,0:Mvol) ) 
                SALLOCATE( iforce,    (-2:2, 0:NGdof), zero)
                SALLOCATE( iposition, (-2:2, 0:NGdof), zero)

                lfactor = psifactor(ii,vvol)     ! this "pre-conditions" the geometrical degrees-of-freedom;
                                
                if( ncpu.eq.1) then

                  do isymdiff = -2, 2 ! symmetric fourth-order, finite-difference used to approximate derivatives;
                    if( isymdiff.eq.0 ) cycle

                    iRbc(1:mn,0:Mvol) = oRbc(1:mn,0:Mvol)
                    iZbs(1:mn,0:Mvol) = oZbs(1:mn,0:Mvol)
                    iRbs(1:mn,0:Mvol) = oRbs(1:mn,0:Mvol)
                    iZbc(1:mn,0:Mvol) = oZbc(1:mn,0:Mvol)

                    ! Perturb geometry
                    if( issym.eq.0 .and. irz.eq.0 ) then
                      iRbc(ii,vvol) = iRbc(ii,vvol) + dRZ * isymdiff ! perturb geometry;
                    else if( issym.eq.0 .and. irz.eq.1 ) then
                      iZbs(ii,vvol) = iZbs(ii,vvol) + dRZ * isymdiff ! perturb geometry;
                    else if( issym.eq.1 .and. irz.eq.0 ) then
                      iRbs(ii,vvol) = iRbs(ii,vvol) + dRZ * isymdiff ! perturb geometry;
                    else if( issym.eq.1 .and. irz.eq.1 ) then
                      iZbc(ii,vvol) = iZbc(ii,vvol) + dRZ * isymdiff ! perturb geometry;
                    endif

                    packorunpack = 'P' ! pack geometrical degrees-of-freedom;
                    LComputeAxis = .false.

                    WCALL(dforce, packxi,( NGdof, iposition(isymdiff,0:NGdof), Mvol, mn,iRbc(1:mn,0:Mvol),iZbs(1:mn,0:Mvol),iRbs(1:mn,0:Mvol),&
                                            iZbc(1:mn,0:Mvol),packorunpack, .false., LComputeAxis ) )
                    WCALL(dforce, dforce,( NGdof, iposition(isymdiff,0:NGdof), iforce(isymdiff,0:NGdof), .false., LComputeAxis) )
                    
                  enddo

                  iforce(0, 0:NGdof)               = ( - 1 * iforce(2,0:NGdof) &
                                                      + 8 * iforce(1,0:NGdof) &
                                                      - 8 * iforce(-1,0:NGdof) &
                                                      + 1 * iforce(-2,0:NGdof))  / ( 12 * dRZ )
                  tdof = (vvol-1) * LGdof + idof
                  finitediff_hessian(1:NGdof, tdof) = iforce(0, 1:NGdof)* lfactor

                  cput = GETTIME

                  iRbc(1:mn,0:Mvol) = oRbc(1:mn,0:Mvol)
                  iZbs(1:mn,0:Mvol) = oZbs(1:mn,0:Mvol)
                  iRbs(1:mn,0:Mvol) = oRbs(1:mn,0:Mvol)
                  iZbc(1:mn,0:Mvol) = oZbc(1:mn,0:Mvol)

                endif

                DALLOCATE(oRbc)
                DALLOCATE(oZbs)
                DALLOCATE(oRbs)
                DALLOCATE(oZbc)
                DALLOCATE(iforce)
                DALLOCATE(iposition)


              endif
#endif
            enddo ! matches do issym ;

          enddo ! matches do irz ;
        enddo ! matches do ii ;

      else ! matches if( ImagneticOK(vvol) .and. ImagneticOK(vvol+1) ) ; 

        FATAL( dforce, .true., need to provide suitable values for hessian in case of field failure )

      endif ! end of if( ImagneticOK(vvol) .and. ImagneticOK(vvol+1) ) ;

    enddo ! end of do vvol;



#ifdef DEBUG

! Print hessian and finite differences estimate (if single CPU). 
    if( Lcheck.eq.6 ) then
      if(myid.eq.0) then
        open(10, file=trim(ext)//'.Lcheck6_output.txt', status='unknown')
        write(ounit,'(A)') NEW_LINE('A')
        do ii=1, NGdof
          write(ounit,1345) myid, im(ii), in(ii), hessian(ii,:)
          write(10   ,1347) hessian(ii,:)
        enddo
        close(10)
        
        write(ounit,'(A)') NEW_LINE('A')

        open(10, file=trim(ext)//'.Lcheck6_output.FiniteDiff.txt', status='unknown')
        if( ncpu.eq.1 ) then
          do ii=1, NGdof
            write(10   ,1347) hessian(ii,:)
          enddo
        endif
        close(10)

          
        write(ounit,'(A)') NEW_LINE('A')

        open(10, file=trim(ext)//'.Lcheck6_output.FiniteDiff.txt', status='unknown')
        if( ncpu.eq.1 ) then
            do ii=1, NGdof
              write(10   ,1347) finitediff_hessian(ii,:)
            enddo        
            write(ounit,'(A)') NEW_LINE('A')
        endif
        close(10)

1345       format("dforce: myid=",i3," ; (",i4,",",i4," ; Hessian            = ",512f16.10 "   ;")
1346       format("dforce: myid=",i3," ; (",i4,",",i4," ; Finite differences = ",512f16.10 "   ;")
1347       format(512F22.16, " ")

      endif
      DALLOCATE(finitediff_hessian)

    FATAL(dforce, Lcheck.eq.6, Lcheck.eq.6 test has been completed. )
    endif
#endif

  endif ! end of if( LcomputeDerivatives ) ;

!call MPI_BARRIER( MPI_COMM_WORLD, ierr )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(dforce)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine dforce

