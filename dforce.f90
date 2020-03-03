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

subroutine dforce( NGdof, position, force, LComputeDerivatives )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half, one, pi
  
  use numerical, only : logtolerance
  
  use fileunits, only : ounit
  
  use inputlist, only : Wmacros, Wdforce, ext, Nvol, Ntor, Lrad, Igeometry, &
                        epsilon, &
                        Lconstraint, Lcheck, &
                        Lextrap, &
                                                mupftol
  
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
                        dpflux, sweight, &
                        Bemn, Bomn, Iomn, Iemn, Somn, Semn, &
                        BBe, IIo, BBo, IIe, & ! these are just used for screen diagnostics;
                        LGdof, dBdX, &
                        Ate, Aze, Ato, Azo, & ! only required for broadcasting
                        diotadxup, dItGpdxtp, & ! only required for broadcasting
                        lBBintegral, &
                        dFFdRZ, dBBdmp, dmupfdx, hessian, dessian, Lhessianallocated, &
                        BBweight, & ! exponential weight on force-imbalance harmonics;
                        psifactor, &
                        DToocc, DToocs, DToosc, DTooss, &
                        TTsscc, TTsscs, TTsssc, TTssss, &
                        TDstcc, TDstcs, TDstsc, TDstss, &
                        TDszcc, TDszcs, TDszsc, TDszss, &
                        DDttcc, DDttcs, DDttsc, DDttss, &
                        DDtzcc, DDtzcs, DDtzsc, DDtzss, &
                        DDzzcc, DDzzcs, DDzzsc, DDzzss, &
                        LocalConstraint, xoffset, &
                        solution, &
                        IsMyVolume, IsMyVolumeValue, WhichCpuID
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in)  :: NGdof               ! dimensions;
  REAL,    intent(in)  :: position(0:NGdof)   ! degrees-of-freedom = internal geometry;
  REAL,    intent(out) :: force(0:NGdof)      ! force;
  LOGICAL, intent(in)  :: LComputeDerivatives ! indicates whether derivatives are to be calculated;
  
  INTEGER              :: vvol, innout, ii, jj, irz, issym, iocons, tdoc, idoc, idof, tdof, jdof, ivol, imn, ll, ihybrd1, lwa, Ndofgl, llmodnp
  INTEGER              :: maxfev, ml, muhybr, mode, nprint, nfev, ldfjac, lr, Nbc, NN, cpu_id
  DOUBLE PRECISION     :: epsfcn, factor
  DOUBLE PRECISION     :: Fdof(1:Mvol-1), Xdof(1:Mvol-1), Fvec(1:Mvol-1)
  DOUBLE PRECISION     :: diag(1:Mvol-1), qtf(1:Mvol-1), wa1(1:Mvol-1), wa2(1:Mvol-1), wa3(1:Mvol-1), wa4(1:mvol-1)
  DOUBLE PRECISION, allocatable :: fjac(:, :), r(:) 

  INTEGER                        :: status(MPI_STATUS_SIZE), request_recv, request_send, cpu_send
  INTEGER              :: id
  INTEGER              :: iflag

  CHARACTER            :: packorunpack 
  EXTERNAL                     :: dfp100, dfp200, loop_dfp100

#ifdef DEBUG
  INTEGER              :: isymdiff
  REAL                 :: dRZ = 1.0e-05, dvol(-1:+1), evolume, imupf(1:2,-2:2)
  REAL,    allocatable :: isolution(:,:)
#endif

  BEGIN(dforce)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! Unpack position to generate arrays iRbc, iZbs, IRbs, iZbc.

  packorunpack = 'U' ! unpack geometrical degrees-of-freedom;

  WCALL( dforce, packxi,( NGdof, position(0:NGdof), Mvol, mn, iRbc(1:mn,0:Mvol), iZbs(1:mn,0:Mvol), iRbs(1:mn,0:Mvol), iZbc(1:mn,0:Mvol), packorunpack ) )
 
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( LcomputeDerivatives ) then
#ifdef DEBUG
   FATAL( dforce, .not.allocated(dBBdmp), do not pass go )
#endif
   dBBdmp(1:LGdof,1:Mvol,0:1,1:2) = zero
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! COMPUTE MATRICES
! ----------------

! Here we call ma00aa to compute the geometry dependent matrices
  do vvol = 1, Mvol

   LREGION(vvol) ! assigns Lcoordinatesingularity, Lplasmaregion, etc. ;

! Determines if this volume vvol should be computed by this thread.
        call IsMyVolume(vvol)

        if( IsMyVolumeValue .EQ. 0 ) then
            cycle
        else if( IsMyVolumeValue .EQ. -1) then
            FATAL(dfp100, .true., Unassociated volume)
        endif

   ll = Lrad(vvol)

   SALLOCATE( DToocc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( DToocs, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( DToosc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( DTooss, (0:ll,0:ll,1:mn,1:mn), zero )

   SALLOCATE( TTsscc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( TTsscs, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( TTsssc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( TTssss, (0:ll,0:ll,1:mn,1:mn), zero )

   SALLOCATE( TDstcc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( TDstcs, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( TDstsc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( TDstss, (0:ll,0:ll,1:mn,1:mn), zero )

   SALLOCATE( TDszcc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( TDszcs, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( TDszsc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( TDszss, (0:ll,0:ll,1:mn,1:mn), zero )

   SALLOCATE( DDttcc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( DDttcs, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( DDttsc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( DDttss, (0:ll,0:ll,1:mn,1:mn), zero )

   SALLOCATE( DDtzcc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( DDtzcs, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( DDtzsc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( DDtzss, (0:ll,0:ll,1:mn,1:mn), zero )

   SALLOCATE( DDzzcc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( DDzzcs, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( DDzzsc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( DDzzss, (0:ll,0:ll,1:mn,1:mn), zero )

   WCALL( dforce, ma00aa, ( Iquad(vvol), mn, vvol, ll ) ) ! compute volume integrals of metric elements - evaluate TD, DT, DD, ...;
   WCALL( dforce, matrix, ( vvol, mn, ll ) )

   DALLOCATE(DToocc)
   DALLOCATE(DToocs)
   DALLOCATE(DToosc)
   DALLOCATE(DTooss)

   DALLOCATE(TTsscc)
   DALLOCATE(TTsscs)
   DALLOCATE(TTsssc)
   DALLOCATE(TTssss)

   DALLOCATE(TDstcc)
   DALLOCATE(TDstcs)
   DALLOCATE(TDstsc)
   DALLOCATE(TDstss)

   DALLOCATE(TDszcc)
   DALLOCATE(TDszcs)
   DALLOCATE(TDszsc)
   DALLOCATE(TDszss)

   DALLOCATE(DDttcc)
   DALLOCATE(DDttcs)
   DALLOCATE(DDttsc)
   DALLOCATE(DDttss)

   DALLOCATE(DDtzcc)
   DALLOCATE(DDtzcs)
   DALLOCATE(DDtzsc)
   DALLOCATE(DDtzss)

   DALLOCATE(DDzzcc)
   DALLOCATE(DDzzcs)
   DALLOCATE(DDzzsc)
   DALLOCATE(DDzzss)

  enddo

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
  if( LocalConstraint ) then

    Ndofgl = 0; Fvec(1:Mvol-1) = 0; iflag = 0;
    Xdof(1:Mvol-1) = dpflux(2:Mvol) + xoffset
    
    ! Solve for field
    WCALL(dforce, dfp100, (Ndofgl, Xdof, Fvec, iflag) )
 
    ! Get force imbalance and jacobian
    do vvol = 1, Mvol
   
        WCALL(dforce, IsMyVolume, (vvol))

        if( IsMyVolumeValue .EQ. 0 ) then
            cycle
        else if( IsMyVolumeValue .EQ. -1) then
            FATAL(dforce, .true., Unassociated volume)
        endif

        WCALL(dforce, dfp200, ( LcomputeDerivatives, vvol) )
    enddo ! end of do vvol = 1, Mvol

! --------------------------------------------------------------------------------------------------
! Global constraint - call the master thread calls hybrd1 on dfp100, others call dfp100_loop.
  else

    Ndofgl = Mvol-1; 

    if( myid.eq. 0) then
          lwa = 8 * Ndofgl * Ndofgl; maxfev = 1000; nfev=0; lr=Mvol*(Mvol-1); ldfjac=Mvol-1
          ml = Mvol-2; muhybr = Mvol-2; epsfcn=1E-16; diag=0.0; mode=1; factor=0.01; nprint=1e5;    !nprint=1e5 to force last call - used for MPI communications

          Xdof(1:Mvol-1)   = dpflux(2:Mvol) + xoffset  ! xoffset reduces the number of iterations needed by hybrd for an obscure reason...

          SALLOCATE(fjac, (1:ldfjac,1:Mvol-1), 0)
          SALLOCATE(r, (1:lr), 0)

          ! Hybrid-Powell method, iterates on all poloidal fluxes to match the global constraint
          WCALL( dforce,  hybrd1, (dfp100, Ndofgl, Xdof(1:Ndofgl), Fvec(1:Ndofgl), mupftol, maxfev, ml, muhybr, epsfcn, diag(1:Ndofgl), mode, &
                      factor, nprint, ihybrd1, nfev, fjac(1:Ndofgl,1:Ndofgl), ldfjac, r(1:lr), lr, qtf(1:Ndofgl), wa1(1:Ndofgl), &
                      wa2(1:Ndofgl), wa3(1:Ndofgl), wa4(1:Ndofgl)) ) 

          DALLOCATE(fjac)
          DALLOCATE(r)
     
          dpflux(2:Mvol) = Xdof(1:Ndofgl) - xoffset

        else

            ! Slave threads call loop_dfp100 and help the master thread computation at each iteration.
            call loop_dfp100(Ndofgl, Fvec, iflag)

    endif

! --------------------------------------------------------------------------------------------------
!                                    MPI COMMUNICATIONS

    ! Gather all ImagneticOK
    do vvol=1, Mvol 
        ! Determine which thread has info on which volume
        call WhichCpuID(vvol, cpu_send) 
            
        ! For now, use MPI_RECV and MPI_SEND. TODO: change implementationo of ImagneticOK to allow the use 
        ! of MPI_GATHER
        if( cpu_send.NE.0    ) then
            if( myid.EQ.0 ) then
                call MPI_RECV(ImagneticOK(vvol), 1, MPI_LOGICAL, cpu_send, vvol, MPI_COMM_WORLD, status, ierr)
            else if( myid.EQ.cpu_send ) then
                call MPI_SEND(ImagneticOK(vvol), 1, MPI_LOGICAL,         0, vvol, MPI_COMM_WORLD, ierr)
            endif
        endif
    enddo

    ! Now master thread broadcast the poloidal flux matching the constraint. It was obtain by iteration
    ! via hybrd1
    call MPI_Bcast( dpflux, Mvol, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! And broadcast as well the ImagneticOK flag - this determines if the computation was succesful in
    ! each volume. If not, this geometry iteration goes to the trash...
    call MPI_Bcast( ImagneticOK, Mvol, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    
    ! And finally broadcast the field information to all threads from the thread which did the computation
    do vvol = 1, Mvol
        call WhichCpuID(vvol, cpu_id)

        NN = NAdof(vvol)
        Nbc = NN * 4
        call MPI_Bcast( solution(vvol)%mat(1:NN, -1:2), Nbc, MPI_DOUBLE_PRECISION, cpu_id, MPI_COMM_WORLD, ierr)
    
        RlBCAST( diotadxup(0:1, -1:2, vvol), 8, cpu_id)
        RlBCAST( dItGpdxtp(0:1, -1:2, vvol), 8, cpu_id)
      
        do ii = 1, mn  
                 RlBCAST( Ate(vvol,0,ii)%s(0:Lrad(vvol)), Lrad(vvol)+1, cpu_id)
                 RlBCAST( Aze(vvol,0,ii)%s(0:Lrad(vvol)), Lrad(vvol)+1, cpu_id)
          enddo

        if( NOTstellsym ) then
            do ii = 1, mn    
                  RlBCAST( Ato(vvol,0,ii)%s(0:Lrad(vvol)), Lrad(vvol)+1, cpu_id)
                   RlBCAST( Azo(vvol,0,ii)%s(0:Lrad(vvol)), Lrad(vvol)+1, cpu_id)
              enddo
        endif
    enddo

    ! --------------------------------------------------------------------------------------------------
    ! Now that all the communication is over, compute the local force and its derivatives
    do vvol = 1, Mvol

        WCALL(dforce, IsMyVolume, (vvol))

        if( IsMyVolumeValue .EQ. 0 ) then
            cycle
        else if( IsMyVolumeValue .EQ. -1) then
            FATAL(dforce, .true., Unassociated volume)
        endif
                
        WCALL(dforce, dfp200, ( LcomputeDerivatives, vvol) )

    enddo

#ifdef DEBUG
      select case( ihybrd1 )
        case( 1   )  ; write(ounit,'("dforce : ",f10.2," : finished ; success        ; dpflux = ", es12.5, ", its="i7";")') cput-cpus, dpflux, nfev
        case( 0   )  ; write(ounit,'("dforce : ",f10.2," : finished ; input error    ; dpflux = ", es12.5, ", its="i7";")') cput-cpus, dpflux, nfev
        case( 2   )  ; write(ounit,'("dforce : ",f10.2," : finished ; max. iter      ; dpflux = ", es12.5, ", its="i7";")') cput-cpus, dpflux, nfev
        case( 3   )  ; write(ounit,'("dforce : ",f10.2," : finished ; xtol too small ; dpflux = ", es12.5, ", its="i7";")') cput-cpus, dpflux, nfev
        case( 4:5 )  ; write(ounit,'("dforce : ",f10.2," : finished ; bad progress   ; dpflux = ", es12.5, ", its="i7";")') cput-cpus, dpflux, nfev
        case default ; write(ounit,'("dforce : ",f10.2," : finished ; illegal ifail  ; dpflux = ", es12.5, ", its="i7";")') cput-cpus, dpflux, nfev
      end select
#endif

endif !matches if( LocalConstraint ) 

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
    
     ;force(tdoc+idoc+1:tdoc+idoc+mn-1  ) = (                           Iomn(2:mn    ,vvol+0  ) ) * epsilon         & ! spectral constraints;
                                          + (                         + Somn(2:mn    ,vvol+0,1) ) * sweight(vvol+0) & ! poloidal length constraint;
                                          - ( Somn(2:mn    ,vvol+1,0)                           ) * sweight(vvol+1)
          
!     if( Ntor.gt.0 ) then ! poloidal angle origin is not otherwise constrained ;
!      force(tdoc+idoc+1:tdoc+idoc+Ntor  ) = ( Pomn(2:Ntor+1,vvol+1,0) - Pomn(2:Ntor+1,vvol+0,1) ) * apsilon ! choice of spectral constraint can be enforced;
!     endif
     
     ;IIo(vvol) = max( sum( abs( force(tdoc+idoc+1:tdoc+idoc+mn-1) ) ) / (mn-1), logtolerance ) ! screen diagnostics;
     
     ; idoc = idoc + mn-1
     
    endif ! end of if( Igeometry.ge.3 ) ;
    
    if( NOTstellsym ) then
     
     ;force(tdoc+idoc+1:tdoc+idoc+mn-1  ) = ( Bomn(2:mn    ,vvol+1,0) - Bomn(2:mn    ,vvol+0,1) ) * BBweight(2:mn) ! pressure imbalance;
     
     ; BBo(vvol) = max( sum( abs( force(tdoc+idoc+1:tdoc+idoc+mn-1) ) ) / (mn-1), logtolerance ) ! screen diagnostics;
     
     ; idoc = idoc + mn-1 ! degree-of-constraint counter; increment;
     
     if( Igeometry.ge.3 ) then ! add spectral constraints;
      
      force(tdoc+idoc+1:tdoc+idoc+mn    ) = (                           Iemn(1:mn    ,vvol+0  ) ) * epsilon         & ! spectral constraints;
                                          + (                         + Semn(1:mn    ,vvol+0,1) ) * sweight(vvol+0) & ! poloidal length constraint;
                                          - ( Semn(1:mn    ,vvol+1,0)                           ) * sweight(vvol+1)
      
!     if( Ntor.ge.0 ) then
!      force(tdoc+idoc+1:tdoc+idoc+Ntor+1) = ( Pemn(1:Ntor+1,vvol+1,0) - Pemn(1:Ntor+1,vvol+0,1) ) * apsilon ! choice of spectral constraint can be enforced;
!     endif
      
      ;IIe(vvol) = max( sum( abs( force(tdoc+idoc+1:tdoc+idoc+mn  ) ) ) / (mn  ), logtolerance ) ! screen diagnostics;
      
      ;idoc = idoc + mn   ! degree-of-constraint counter; increment;
      
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
        
        if( vvol.gt.1 ) then
         
         tdof = (vvol-2) * LGdof + idof ! labels degree-of-freedom in internal interface geometry   ;
         tdoc = (vvol-1) * LGdof        ! labels force-balance constraint across internal interfaces;
         idoc = 0                       ! local  force-balance constraint across internal interface ;
         hessian(tdoc+idoc+1:tdoc+idoc+LGdof,tdof) =                                           - dFFdRZ(idoc+1:idoc+LGdof,vvol+0,1,idof,0)
         if( Lconstraint.eq.1 ) then ! this is a little clumsy; could include Lfreebound or something . . . ;
         hessian(tdoc+idoc+1:tdoc+idoc+LGdof,tdof) = hessian(tdoc+idoc+1:tdoc+idoc+LGdof,tdof)                     &
                                                   - dBBdmp(idoc+1:idoc+LGdof,vvol+0,1,1) * dmupfdx(vvol,1,idof,0) &
                                                   - dBBdmp(idoc+1:idoc+LGdof,vvol+0,1,2) * dmupfdx(vvol,2,idof,0)
         endif ! end of if( Lconstraint.eq.1 ) ; 
      
        endif ! end of if( vvol.gt.1 ) ;
        

        ;tdof = (vvol-1) * LGdof + idof
        ;tdoc = (vvol-1) * LGdof ! shorthand;
        ;idoc = 0
        if( Lextrap.eq.1 .and. vvol.eq.1 ) then
        ;hessian(tdoc+idof                  ,tdof) = one ! diagonal elements;
        else
        ;hessian(tdoc+idoc+1:tdoc+idoc+LGdof,tdof) = dFFdRZ(idoc+1:idoc+LGdof,vvol+1,0,idof,0) - dFFdRZ(idoc+1:idoc+LGdof,vvol+0,1,idof,1)
         if( Lconstraint.eq.1 ) then ! this is a little clumsy;
         hessian(tdoc+idoc+1:tdoc+idoc+LGdof,tdof) = hessian(tdoc+idoc+1:tdoc+idoc+LGdof,tdof)                       &
                                                   + dBBdmp(idoc+1:idoc+LGdof,vvol+1,0,1) * dmupfdx(vvol+1,1,idof,0) &
                                                   + dBBdmp(idoc+1:idoc+LGdof,vvol+1,0,2) * dmupfdx(vvol+1,2,idof,0) &
                                                   - dBBdmp(idoc+1:idoc+LGdof,vvol+0,1,1) * dmupfdx(vvol+0,1,idof,1) &
                                                   - dBBdmp(idoc+1:idoc+LGdof,vvol+0,1,2) * dmupfdx(vvol+0,2,idof,1)
         endif ! end of if( Lconstraint.eq.1 );
         endif
        
         if( vvol.lt.Mvol-1 ) then

         tdof = (vvol+0) * LGdof + idof
         tdoc = (vvol-1) * LGdof ! shorthand;
         idoc = 0
         if( Lextrap.eq.1 .and. vvol.eq.1 ) then
         if    ( im(idof).le.0                     ) then ; hessian(tdoc+idof,tdof) = - one
         else                                             ; hessian(tdoc+idof,tdof) = - one
         endif
         else
         hessian(tdoc+idoc+1:tdoc+idoc+LGdof,tdof) = dFFdRZ(idoc+1:idoc+LGdof,vvol+1,0,idof,1)
         if( Lconstraint.eq.1 ) then ! this is a little clumsy;
         hessian(tdoc+idoc+1:tdoc+idoc+LGdof,tdof) = hessian(tdoc+idoc+1:tdoc+idoc+LGdof,tdof)                       &
                                                   + dBBdmp(idoc+1:idoc+LGdof,vvol+1,0,1) * dmupfdx(vvol+1,1,idof,1) &
                                                   + dBBdmp(idoc+1:idoc+LGdof,vvol+1,0,2) * dmupfdx(vvol+1,2,idof,1)
         endif ! end of if( Lconstraint.eq.1 ) then;
         endif

        endif ! end of if( vvol.lt.Mvol-1 ) ;

        if( vvol.eq.Mvol-1 ) then
        !tdof = (vvol+0) * LGdof + idof
         tdoc = (vvol-1) * LGdof ! shorthand ;
         idoc = 0
         dessian(tdoc+idoc+1:tdoc+idoc+LGdof,idof) = dFFdRZ(idoc+1:idoc+LGdof,vvol+1,0,idof,1)
         if( Lconstraint.eq.1 ) then ! this is a little clumsy;
         dessian(tdoc+idoc+1:tdoc+idoc+LGdof,idof) = dessian(tdoc+idoc+1:tdoc+idoc+LGdof,idof)                       &
                                                   + dBBdmp(idoc+1:idoc+LGdof,vvol+1,0,1) * dmupfdx(vvol+1,1,idof,1) &
                                                   + dBBdmp(idoc+1:idoc+LGdof,vvol+1,0,2) * dmupfdx(vvol+1,2,idof,1)
         endif ! end of if( Lconstraint.eq.1 ) then;

         
        endif ! end of if( vvol.lt.Mvol-1 ) ;
        
       enddo ! matches do issym ;
       
      enddo ! matches do irz ;
      
     enddo ! matches do ii ;
     
    else ! matches if( ImagneticOK(vvol) .and. ImagneticOK(vvol+1) ) ; 
     
     FATAL( dforce, .true., need to provide suitable values for hessian in case of field failure )
     
    endif ! end of if( ImagneticOK(vvol) .and. ImagneticOK(vvol+1) ) ;
    
   enddo ! end of do vvol;
   
  endif ! end of if( LcomputeDerivatives ) ;
  
  call MPI_BARRIER( MPI_COMM_WORLD, ierr )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(dforce)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine dforce
 
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


subroutine loop_dfp100(Ndofgl, Fvec, iflag)

! LOOP_DFP100 - infinite loop for slaves helping the master thread iterating to match global 
! constraint
  
  use fileunits, only : ounit                                                                ! Unit identifier for write
  
  use inputlist, only : Wmacros, Wdforce                                        ! Flags for debugging

    use cputiming, only     :  Tdforce                                                    ! Timer
    use allglobal, only     :  Mvol, &                                                    ! Total number of volume + vacuum
                               IconstraintOK, &                                    ! Flag to exit loop
                               cpus, myid

 LOCALS
!------

    INTEGER                                              :: Ndofgl, iflag            ! Input parameters to dfp100
    DOUBLE PRECISION                                     :: Fvec(1:Mvol-1), x(1:Mvol-1) ! Input parameters to dfp100
    EXTERNAL                                                        :: dfp100                            ! Field solver

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! Initialize - for now the constraint is not matched
    IconstraintOK = .false.

! Enter the infinit loop. Master thread broadcasts IconstraintOK at each iteration - the value is
! .TRUE. when the constraint is matched
    do while (.not.IconstraintOK)

! Compute solution in every associated volumes
        WCALL(dforce, dfp100, (Ndofgl, x, Fvec, iflag) )

    end do !matches do while IconstraintOK

end subroutine loop_dfp100
