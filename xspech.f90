!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (main) ! Main program.

!latex \briefly{main program}

!latex \calledby{}

!latex \calls{\link{global}:readin,
!latex        \link{global}:wrtend,
!latex        \link{preset},
!latex        \link{packxi},
!latex        \link{volume},
!latex        \link{pc00aa},
!latex        \link{newton},
!latex        \link{dforce},
!latex        \link{hesian},
!latex        \link{ra00aa},
!latex        \link{bnorml},
!latex        \link{sc00aa},
!latex        \link{jo00aa},
!latex        \link{pp00aa} and
!latex        \link{sphdf5}}

!latex \tableofcontents

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! uppercase generally indicates macros;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
program xspech

  use allglobal, only: readin, set_mpi_comm, myid, ncpu, cpus, version
  use inputlist, only: initialize_inputs
  use fileunits, only: ounit

  LOCALS

  call MPI_INIT( ierr )

  ! set default communicator to MPI_COMM_WORLD
  call set_mpi_comm(MPI_COMM_WORLD)

  ! set initial time; 04 Dec 14;
  cpus = GETTIME
  cpuo = cpus

  ! print header: version of SPEC and compilation info
  cput = GETTIME
  if( myid.eq.0 ) then
   write(ounit,'("xspech : ", 10x ," : version = "F5.2)') version
! COMPILATION ! do not delete; this line is replaced (see Makefile) with a write statement identifying date, time, compilation flags, etc.;
  endif

  ! read command-line arguments
  call read_command_args()

  ! initialize input arrays into a default state
  call initialize_inputs()

  cput = GETTIME
  if( myid.eq.0 ) then
   write(ounit,'("xspech : ", 10x ," : ")')
   write(ounit,'("xspech : ",f10.2," : begin execution ; ncpu=",i3," ; calling global:readin ;")') cput-cpus, ncpu
  endif

  call readin() ! read & broadcast input namelist
  ! initialize; readin; default





  call spec() ! main subroutine

  call ending() ! print ending info

  call write_hdf5() ! write HDF5 output

  MPIFINALIZE

  stop
end program xspech

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{input file extension $\equiv$ command line argument}

!latex \begin{enumerate}
!latex \item The input file name, \type{ext}, is given as the first command line input, and the input file itself is \verb!ext.sp!
!latex \item Additional command line inputs recognized are:
!latex \begin{enumerate}
!latex \item \type{-help, -h} ; will give help information to user; under construction;
!latex \item \type{-readin} ; will immediately set \type{Wreadin=T}; this may be over-ruled when \type{namelist/screenlist/} is read;
!latex \end{enumerate}

!latex \end{enumerate}

! read command-line arguments; in particular, determine input file (name or extension)
subroutine read_command_args

  use fileunits, only: ounit
  use inputlist, only: Wreadin
  use allglobal, only: cpus, myid, ext, MPI_COMM_SPEC

  LOCALS

  LOGICAL              :: Lspexist
  INTEGER              :: iargc, iarg, numargs, extlen, sppos

  CHARACTER(len=100)   :: arg

  if (myid.eq.0) then

  cput = GETTIME

  call getarg( 1, arg )
  extlen = len_trim(arg)
  sppos = index(arg, ".sp", .true.) ! search for ".sp" from the back of ext
  if (sppos.eq.extlen-2) then       ! check if ext ends with ".sp";
    arg = arg(1:extlen-3)           ! if this is the case, remove ".sp" from end of ext
  endif
  ext = trim(arg)

  if( ext .eq. "" .or. ext .eq. "-h" .or. ext .eq. "-help" ) then
   ;write(ounit,'("rdcmdl : ", 10x ," : ")')
   ;write(ounit,'("rdcmdl : ", 10x ," : file extension must be given as first command line argument ; extra command line options = -help -readin ;")')
   if( ext .eq. "-h" .or. ext .eq. "-help" ) then
    write(ounit,'("rdcmdl : ", 10x ," : ")')
    write(ounit,'("rdcmdl : ", 10x ," : the input file ext.sp must contain the input namelists; see global.pdf for description ;")')
   endif
   FATAL( rdcmdl, .true., the input file does not exist) ! if not, abort;
  endif

  write(ounit,'("rdcmdl : ", 10x ," : ")')
  write(ounit,'("rdcmdl : ",f10.2," : ext = ",a100)') cput-cpus, ext

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  numargs = iargc()

  if( numargs.gt.1 ) then
   iarg = 1
   do while ( iarg < numargs )
    iarg = iarg + 1 ; call getarg( iarg, arg)
    select case( arg )
    case("-help","-h") ; write(ounit,'("rdcmdl : ",f10.2," : myid=",i3," : command line options = -readin ;")') cput-cpus, myid
    case("-readin"   ) ; Wreadin = .true.
    case("-p4pg"     ) ; iarg = iarg + 1 ; call getarg( iarg, arg) ! TODO: what is this?
    case("-p4wd"     ) ; iarg = iarg + 1 ; call getarg( iarg, arg) ! TODO: what is this?
    case default       ; write(ounit,'("rdcmdl : ",f10.2," : myid=",i3," : argument not recognized ; arg = ",a100)') cput-cpus, myid, arg
    end select
   enddo
  endif

  end if ! check for myid.eq.0

end subroutine read_command_args

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine spec

  use constants, only : zero, one, pi2, mu0

  use numerical, only : vsmall, logtolerance

  use fileunits, only : ounit, lunit

  use inputlist, only : Wmacros, Wxspech, &
                        Nfp, Igeometry, Nvol, Lrad, &
                        tflux, pflux, phiedge, pressure, pscale, helicity, Ladiabatic, adiabatic, gamma, &
                        Rbc, Zbs, Rbs, Zbc, &
                        Lconstraint, &
                        Lfreebound, mfreeits, gBntol, gBnbld, vcasingtol, LautoinitBn, &
                        Lfindzero, LautoinitBn, &
                        odetol, nPpts, nPtrj, &
                        LHevalues, LHevectors, LHmatrix, Lperturbed, Lcheck, &
                        Lzerovac, &
                        mu, Isurf, Ivolume

  use cputiming, only : Txspech

  use allglobal, only : readin, wrtend, ncpu, myid, cpus, ext, &
                        Mvol, &
                        YESstellsym, NOTstellsym, &
                        Iquad, &
                        mn, im, in, &
                        Ntz, &
                        LGdof, NGdof, &
                        iRbc, iZbs, iRbs, iZbc, &
                        BBe, IIo, BBo, IIe, &
                        Btemn, Bzemn, Btomn, Bzomn, &
                        vvolume, &
                        Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, &
                        dtflux, dpflux, &
                        ImagneticOK, &
                        ForceErr, &
                        efmn, ofmn, cfmn, sfmn, &
                        iBns, iBnc, iVns, iVnc, &
                        Ate, Aze, Ato, Azo, & ! only required for debugging; 09 Mar 17;
                        nfreeboundaryiterations, &
                        beltramierror, &
                        first_free_bound, &
                        dMA, dMB, dMD, dMG, MBpsi, solution, dtflux, IPDt, &
                        version, &
                        MPI_COMM_SPEC


   ! write _all_ output quantities into a _single_ HDF5 file
   use sphdf5,   only : init_outfile, &
                        mirror_input_to_outfile, &
                        init_convergence_output, &
                        write_grid

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS

  LOGICAL              :: LComputeDerivatives, LContinueFreeboundaryIterations, exist, LupdateBn, LComputeAxis

! INTEGER              :: nfreeboundaryiterations, imn, lmn, lNfp, lim, lin, ii, lvol ! 09 Mar 17;
  INTEGER              :: imn, lmn, lNfp, lim, lin, ii, ideriv, stat, iocons, iwait, pid, status
  INTEGER              :: vvol, llmodnp, ifail, wflag, iflag, vflag
  REAL                 :: rflag, lastcpu, bnserr, lRwc, lRws, lZwc, lZws, lItor, lGpol, lgBc, lgBs, sumI
  REAL,    allocatable :: position(:), gradient(:), Bt00(:,:,:)
  CHARACTER            :: pack, hostname
  INTEGER              :: Lfindzero_old, mfreeits_old
  REAL                 :: gBnbld_old
  INTEGER              :: lnPtrj, numTrajTotal
  INTEGER, external    :: getpid, hostnm

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!



!#ifdef DEBUG
!   iwait = 0; pid = getpid()
!   status = hostnm( hostname )
!   write(*,*) 'Process with PID: ', pid, 'ready to attach. Hostname: ', hostname
!   do while( iwait .EQ. 0 )
!     !wait for debugger
!   enddo
!#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{reading input, allocating global variables}

!latex \begin{enumerate}
!latex \item The input namelists and geometry are read in via a call to \link{global}\verb+:readin+.
!latex       A full description of the required input is given in \link{global}.
!latex \item Most internal variables, global memory etc., are allocated in \link{preset}.
!latex \item All quantities in the input file are mirrored into the output file's group \type{input}.
!latex \end{enumerate}

  WCALL( xspech, preset )

  WCALL( xspech, init_outfile ) ! initialize HDF5 library and open output file ext.h5 for writing during execution

  WCALL( xspech, mirror_input_to_outfile ) ! mirror input file contents to output file

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if ( myid .eq. 0 ) then ! save restart file;
    WCALL( xspech, wrtend ) ! write restart file ! 17 May 19
  endif

!latex \subsection{preparing output file group \type{iterations}}

!latex \begin{enumerate}

!latex \item The group \verb+iterations+ is created in the output file.
!latex       This group contains the interface geometry at each iteration, which is useful for constructing movies illustrating the convergence.
!latex       The data structure in use is an unlimited array of the following compound datatype:
!latex \begin{verbatim} DATATYPE  H5T_COMPOUND {
!latex       H5T_NATIVE_INTEGER "nDcalls";
!latex       H5T_NATIVE_DOUBLE "Energy";
!latex       H5T_NATIVE_DOUBLE "ForceErr";
!latex       H5T_ARRAY { [Mvol+1][mn] H5T_NATIVE_DOUBLE } "iRbc";
!latex       H5T_ARRAY { [Mvol+1][mn] H5T_NATIVE_DOUBLE } "iZbs";
!latex       H5T_ARRAY { [Mvol+1][mn] H5T_NATIVE_DOUBLE } "iRbs";
!latex       H5T_ARRAY { [Mvol+1][mn] H5T_NATIVE_DOUBLE } "iZbc";
!latex } \end{verbatim}
!latex \end{enumerate}

  WCALL( xspech, init_convergence_output ) ! initialize convergence output arrays ! 17 May 19

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  FATAL( xspech, NGdof.lt.0, counting error )

  SALLOCATE( position, (0:NGdof), zero ) ! position ; NGdof = #geometrical degrees-of-freedom was computed in preset;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  nfreeboundaryiterations = -1

! This is the free-boundary iteration loop (implemented using GOTO); 08 Jun 16;
9000 nfreeboundaryiterations = nfreeboundaryiterations + 1

  ! run fix_boundary for the first free_boundary iteration if LautoinitBn was set to 1
  if (Lfreebound.eq.1 .and. LautoinitBn.eq.1) then
     if (nfreeboundaryiterations.eq.0) then  ! first iteration
        first_free_bound = .true.
        !Mvol = Nvol
        gBnbld_old = gBnbld
        gBnbld = zero
        Lfindzero_old = Lfindzero
        mfreeits_old = mfreeits
        Lfindzero = 0
        mfreeits = 1
        if (myid.eq.0) write(ounit,'("xspech : ",10X," : First iteration of free boundary calculation : update Bns from plasma.")')
     else
        first_free_bound = .false.
        !Mvol = Nvol + Lfreebound
        Lfindzero = Lfindzero_old
        gBnbld = gBnbld_old
        mfreeits = mfreeits_old
     endif
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{``packing" geometrical degrees-of-freedom into vector}

!latex \begin{enumerate}
!latex \item If \internal{NGdof.gt.0}, where \internal{NGdof}$\equiv $ counts the geometrical degrees-of-freedom, i.e. the $R_{bc}$, $Z_{bs}$, etc.,
!latex       then \link{packxi} is called to ``pack" the geometrical degrees-of-freedom into \internal{position(0:NGdof)}.
!latex \end{enumerate}

  if( NGdof.gt.0 ) then ! pack geometry into vector; 14 Jan 13;

   pack = 'P'
   LComputeAxis = .true.
   WCALL( xspech, packxi, ( NGdof, position(0:NGdof), Mvol, mn, iRbc(1:mn,0:Mvol), iZbs(1:mn,0:Mvol), &
                            iRbs(1:mn,0:Mvol), iZbc(1:mn,0:Mvol), pack, .false., LComputeAxis ) )

  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{initialize adiabatic constants}

!latex \begin{enumerate}
!latex \item If \inputvar{Ladiabatic.eq.0}, then the ``adiabatic constants'' in each region, $P_v$, are calculated as
!latex       \be P_v \equiv p_v V_v^\gamma, \label{eq:adiabatic}
!latex       \ee
!latex       where $p_v\equiv$ \inputvar{pressure(vvol)}, the volume $V_v$ of each region is computed by \link{volume},
!latex       and the adiabatic index $\gamma\equiv$ \inputvar{gamma}.
!latex \end{enumerate}

  do vvol = 1, Mvol

   vflag = 0
   WCALL( xspech, volume, ( vvol, vflag ) ) ! compute volume;

   if( Ladiabatic.eq.0 ) adiabatic(vvol) = pressure(vvol) * vvolume(vvol)**gamma ! initialize adiabatic constants using supplied pressure profile;

  enddo ! end of do vvol = 1, Mvol;

  if( Mvol.gt.Nvol ) then ; adiabatic(Mvol) = zero ; pressure(Mvol) = zero ! these are never used; 15 May 13;
  endif

  if( Wxspech .and. myid.eq.0 ) then
   cput = GETTIME
   write(ounit,'("xspech : ",f10.2," : myid=",i3," ; adiabatic constants = "999es13.5)') cput-cpus, myid, adiabatic(1:Mvol)
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{solving force-balance}

!latex \begin{enumerate}
!latex \item If there are geometrical degress of freedom, i.e. if \internal{NGdof.gt.0}, then \bi
!l tex \item[] if \inputvar{Lminimize.eq.1}, call \link{pc00aa} to find minimum of energy functional
!l tex         using quasi-Newton, preconditioned conjugate gradient method, \nag{www.nag.co.uk/numeric/FL/manual19/pdf/E04/e04dgf_fl18.pdf}{E04DGF};
!latex \item[] If \inputvar{Lfindzero.gt.0}, call \link{newton} to find extremum of constrained energy functional using a Newton method,
!latex         \nag{www.nag.co.uk/numeric/FL/manual19/pdf/C05/c05pdf_fl19.pdf}{C05PDF};
!latex \ei
!latex \end{enumerate}

!    if( Lminimize.eq.1 ) then
!
!     ifail = 1 ! this is probably not required; 26 Feb 13;
!
!     WCALL(xspech,pc00aa,( NGdof, position(1:NGdof), Mvol, mn, ifail ))
!
!    endif

  if( NGdof.gt.0 ) then

   if( Lfindzero.gt.0 ) then

   ! This is the call to do one fixed-boundary iteration (by a Newton method).
    ifail = 1
    WCALL( xspech, newton, ( NGdof, position(0:NGdof), ifail ) )

   endif

   pack = 'U' ! unpack geometrical degrees of freedom; 13 Sep 13;
   LComputeAxis = .true.
   WCALL( xspech, packxi, ( NGdof, position(0:NGdof), Mvol, mn, iRbc(1:mn,0:Mvol), iZbs(1:mn,0:Mvol), &
                            iRbs(1:mn,0:Mvol), iZbc(1:mn,0:Mvol), pack, .false., LComputeAxis ) )

  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{post diagnostics}

!latex \begin{enumerate}
!latex \item The pressure is computed from the adiabatic constants from \Eqn{adiabatic}, i.e. $p=P/V^\gamma$.
!latex \item The Beltrami/vacuum fields in each region are re-calculated using \link{dforce}.
!latex \item If \inputvar{Lcheck.eq.5 .or. LHevalues .or. LHevectors .or. Lperturbed.eq.1}, then the ``Hessian'' matrix is examined using \link{hesian}.
!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! set/reset input variables;

! if( Lconstraint.lt.2 ) helicity(1:Nvol) = lABintegral(1:Nvol) ! updated ``input'' quantity;

#ifdef DEBUG
  do vvol = 1, Mvol
   FATAL( xspech, vvolume(vvol).lt.vsmall, error dividing adiabatic by volume )
  enddo
#endif

  pressure(1:Mvol) = adiabatic(1:Mvol) / vvolume(1:Mvol)**gamma ! this matches construction of adiabatic above;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  SALLOCATE( gradient, (0:NGdof), zero )

  lastcpu = GETTIME

  LComputeDerivatives = .false.
  LComputeAxis = .true.
! vvol = Mvol ; ideriv = 0 ; ii = 1
! write(ounit,'("xspech : ", 10x ," : sum(Ate(",i3,",",i2,",",i2,")%s) =",99es23.15)') vvol, ideriv, ii, sum(Ate(vvol,ideriv,ii)%s(0:Lrad(vvol)))

  WCALL( xspech, dforce, ( NGdof, position(0:NGdof), gradient(0:NGdof), LComputeDerivatives, LComputeAxis) ) ! (re-)calculate Beltrami fields;

  DALLOCATE(gradient)

#ifdef DEBUG
  do vvol = 1, Mvol-1
   ; FATAL( xspech, BBe(vvol).lt.logtolerance, underflow )
   if( Igeometry.eq.3 ) then ! include spectral constraints; 04 Dec 14;
    ;FATAL( xspech, IIo(vvol).lt.logtolerance, underflow )
   endif
   if( NOTstellsym ) then
    ;FATAL( xspech, BBo(vvol).lt.logtolerance, underflow )
    if( Igeometry.eq.3 ) then ! include spectral constraints; 04 Dec 14;
     FATAL( xspech, IIe(vvol).lt.logtolerance, underflow )
    endif
   endif
  enddo
#endif

  if( myid.eq.0 ) then
   cput = GETTIME
   write(ounit,1000) cput-cpus, nfreeboundaryiterations,          ForceErr,  cput-lastcpu, "|BB|e", alog10(BBe(1:min(Mvol-1,28)))
   if( Igeometry.ge.3 ) then ! include spectral constraints; 04 Dec 14;
   write(ounit,1001)                                                                       "|II|o", alog10(IIo(1:min(Mvol-1,28)))
   endif
   if( NOTstellsym ) then
   write(ounit,1001)                                                                       "|BB|o", alog10(BBo(1:min(Mvol-1,28)))
   if( Igeometry.ge.3 ) then ! include spectral constraints; 04 Dec 14;
   write(ounit,1001)                                                                       "|II|e", alog10(IIe(1:min(Mvol-1,28)))
   endif
   endif
  endif

1000 format("xspech : ",f10.2," : #freeits=",i3," ; ":"|f|="es12.5" ; ":"time=",f10.2,"s ;":" log"a5,:"="28f6.2" ...")
1001 format("xspech : ", 10x ," :          ",3x," ; ":"    "  12x "   ":"     ", 10x ,"  ;":" log"a5,:"="28f6.2" ...")

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( Lcheck.eq.5 .or. LHevalues .or. LHevectors .or. LHmatrix .or. Lperturbed.eq.1 ) then ! check construction of Hessian; 01 Jul 14;

   if( myid.eq.0 ) then
    cput = GETTIME
    write(ounit,'("xspech : ", 10x ," : ")')
    write(ounit,'("xspech : ",f10.2," : myid=",i3," ; calling hesian ; see .ext.hessian.myid ;")') cput-cpus, myid
   endif

   WCALL( xspech, hesian, ( NGdof, position(0:NGdof), Mvol, mn, LGdof ) )

  endif ! end of if( Lcheck.eq.5 ) ; 01 Jul 14;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  select case( Igeometry )                                  ! 08 Feb 16;
  case( 1   ) ; tflux(1) = dtflux(1) ; pflux(1) = dpflux(1) ! 08 Feb 16;
  case( 2:3 ) ; tflux(1) = dtflux(1) ; pflux(1) =   zero    ! 08 Feb 16;
  end select                                                ! 08 Feb 16;

  do vvol = 2, Mvol; tflux(vvol) = tflux(vvol-1) + dtflux(vvol) ! 01 Jul 14;
   ;                  pflux(vvol) = pflux(vvol-1) + dpflux(vvol) ! 01 Jul 14;
  enddo

  tflux(1:Mvol) = tflux(1:Mvol) * pi2 / phiedge ! this is the "inverse" operation defined in preset; 19 Jul 16;
  pflux(1:Mvol) = pflux(1:Mvol) * pi2 / phiedge

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{free-boundary: re-computing normal field}

!latex \begin{enumerate}
!latex \item If \inputvar{Lfreebound.eq.1} and \inputvar{Lfindzero.gt.0} and \inputvar{mfreeits.ne.0},
!latex       then the magnetic field at the computational boundary produced by the plasma currents is computed using \link{bnorml}.
!latex \item The ``new'' magnetic field at the computational boundary produced by the plasma currents is updated using a Picard scheme:
!latex       \be \verb+Bns+_i^{j} = \lambda \, \verb+Bns+_i^{j-1} + (1-\lambda) \verb+Bns+_i, \label{eq:blending}
!latex       \ee
!latex       where $j$ labels free-boundary iterations, the ``blending parameter'' is $\lambda\equiv $ \inputvar{gBnbld},
!latex       and \verb+Bns+$_i$ is computed by virtual casing.
!latex       The subscript ``$i$'' labels Fourier harmonics.
!latex \item If the new (unblended) normal field is {\em not} sufficiently close to the old normal field, as quantified by \inputvar{gBntol},
!latex       then the free-boundary iterations continue.
!latex       This is quantified by
!latex       \be \sum_i | \verb+Bns+_i^{j-1} - \verb+Bns+_i | / N, \label{eq:gBntol}
!latex       \ee
!latex       where $N$ is the total number of Fourier harmonics.
!latex \item There are several choices that are available:
!latex       \begin{enumerate}
!latex       \item if \inputvar{mfreeits} $= -2$ : the vacuum magnetic field
!latex             (really, the normal component of the field produced by the external currents at the computational boundary)
!latex             required to hold the given equlibrium is written to file.
!latex             This information is required as input by
!latex             \paper{FOCUS}{Caoxiang Zhu, Stuart Hudson et al.}{10.1088/1741-4326/aa8e0a}{Nucl. Fusion}{58}{016008}{2017}
!latex             for example. (This option probably needs to revised.)
!latex       \item if \inputvar{mfreeits} $= -1$ : after the plasma field is computed by virtual casing,
!latex             the vacuum magnetic field is set to exactly balance the plasma field
!latex             (again, we are really talking about the normal component at the computational boundary.)
!latex             This will ensure that the computational boundary itself if a flux surface of the total magnetic field.
!latex       \item if \inputvar{mfreeits} $= 0$ : the plasma field at the computational boundary is not updated; no ``free-boundary'' iterations take place.
!latex       \item if \inputvar{mfreeits} $> 0$ : the plasma field at the computational boundary is updated according to the above blending \Eqn{blending},
!latex             and the free-boundary iterations will continue until either the tolerance condition is met (see \inputvar{gBntol} and \Eqn{gBntol})
!latex             or the maximum number of free-boundary iterations, namely \inputvar{mfreeits}, is reached.
!latex             For this case, \inputvar{Lzerovac} is relevant:
!latex             if \inputvar{Lzerovac} $= 1$, then the vacuum field is set equal to the normal field at every iteration,
!latex             which results in the computational boundary being a flux surface.
!latex             (I am not sure if this is identical to setting \inputvar{mfreeits}$= -1$; the logic etc. needs to be revised.)
!latex       \end{enumerate}
!latex \end{enumerate}

  LContinueFreeboundaryIterations = .false.

  ;                                                              LupdateBn = .false. ! default;
!  if( Lfreebound.eq.1 .and. Lfindzero.gt.0 ) then
  if( Lfreebound.eq.1) then   ! removed Lfindzero check; Loizu Dec 18;
   if( mfreeits.gt.0 .and. nfreeboundaryiterations.lt.mfreeits ) LupdateBn = .true.
   if( mfreeits.lt.0                                           ) LupdateBn = .true.
  endif

  if( LupdateBn ) then

   Mvol = Nvol + Lfreebound

   lastcpu = GETTIME

   WCALL( xspech, bnorml, ( mn, Ntz, efmn(1:mn), ofmn(1:mn) ) ) ! compute normal field etc. on computational boundary;

   !FATAL( xspech, mn-1.le.0, divide by zero )

   if(mn.eq.1) then
     if( YESstellsym ) bnserr = 0.0 !TODO: NOT SURE, this should test if bns is actually 0
     if( NOTstellsym ) bnserr = sum( abs( iBnc(1:mn) - efmn(1:mn) ) ) / (mn  )
   else
     if( YESstellsym ) bnserr = sum( abs( iBns(2:mn) - ofmn(2:mn) ) ) / (mn-1)
     if( NOTstellsym ) bnserr = sum( abs( iBns(2:mn) - ofmn(2:mn) ) ) / (mn-1) &
                              + sum( abs( iBnc(1:mn) - efmn(1:mn) ) ) / (mn  )
   endif


   if( bnserr.gt.gBntol ) then

    LContinueFreeboundaryIterations = .true.

    select case( mfreeits )

    case( -2  ) ! mfreeits = -2 ; shall set plasma normal field at computational boundary ; 24 Nov 16;

     inquire( file=trim(ext)//".Vn", exist=exist )
     FATAL( xspech, .not.exist, ext.Vn does not exist : cannot set vacuum field)

     if( myid.eq.0 ) then ! only myid = 0 reads in the vacuum field; 04 Jan 17;

      ;                    iBns(2:mn) = iVns(2:mn) - iBns(2:mn) ! temporary storage of the total field; 07 Dec 16;
      if( NOTstellsym)     iBnc(1:mn) = iVnc(1:mn) - iBnc(1:mn)

      open(lunit, file = trim(ext)//".Vn", status="old" )
      read(lunit,*) lmn, lNfp
      read(lunit,*) ( lim, lin, lRwc, lRws, lZwc, lZws, imn = 1, lmn ) ! control surface; should confirm agreement; 04 Jan 17;
      read(lunit,*) lmn, lNfp, lItor, lGpol
      do imn = 1, lmn
       read(lunit,*) lim, lin, lgBc, lgBs
       do ii = 1, mn
        if( lim.eq.im(ii) .and. lin*lNfp.eq.in(ii) ) then
         ;                 iVns(ii) = lgBs
         if( NOTstellsym ) iVnc(ii) = lgBc
        endif
       enddo
      enddo
      close(lunit)

     endif ! end of if( myid.eq.0 ) ; 07 Dec 16;

     ;RlBCAST( iVns(1:mn), mn, 0 ) ! only required for ii > 1 ;
     if( NOTstellsym ) then
      RlBCAST( iVnc(1:mn), mn, 0 )
     endif

     ;                    iBns(2:mn) = - iBns(2:mn) - iVns(2:mn) ! updated vacuum field ; 24 Nov 16;
     if( NOTstellsym)     iBnc(1:mn) = - iBnc(1:mn) - iVnc(1:mn)

     if( myid.eq.0 ) then
     write(ounit,'("xspech : " 10x " : oBns=[",999(es11.03,","))') iBns(1:mn) ! 17 Jan 17;
     write(ounit,'("xspech : " 10x " : nBns=[",999(es11.03,","))') ofmn(1:mn) ! 17 Jan 17;
     if( NOTstellsym ) then
     write(ounit,'("xspech : " 10x " : oBnc=[",999(es11.03,","))') iBnc(1:mn) ! 17 Jan 17;
     write(ounit,'("xspech : " 10x " : nBnc=[",999(es11.03,","))') efmn(1:mn) ! 17 Jan 17;
     endif
     endif

    case( -1  ) ! mfreeits = -1 ; shall set vacuum normal field at computational boundary ; 24 Nov 16;

     ;                    iVns(2:mn) = iVns(2:mn) - iBns(2:mn)                                 ! total   normal field              ; 24 Nov 16;
     if( NOTstellsym)     iVnc(1:mn) = iVnc(1:mn) - iBnc(1:mn)
     ;                    iBns(2:mn) =                                              ofmn(2:mn) ! updated normal field due to plasma; 24 Nov 16;
     if( NOTstellsym)     iBnc(1:mn) =                                              efmn(1:mn)
     ;                    iVns(2:mn) = iVns(2:mn) + iBns(2:mn)                                 ! updated vacuum field              ; 24 Nov 16;
     if( NOTstellsym)     iVnc(1:mn) = iVnc(1:mn) + iBnc(1:mn)

    case(  0  ) ! mfreeits =  0 ; 09 Mar 17;

     FATAL( xspech, .true., illegal mfreeits logic )

    case(  1: ) ! mfreeits >  0 ; 09 Mar 17;

     select case( Lzerovac )

     case( 0 ) ! Lzerovac = 0 ; 09 Mar 17;

      ;                iBns(2:mn) = gBnbld * iBns(2:mn) + ( one - gBnbld ) * ofmn(2:mn)
      if( NOTstellsym) iBnc(1:mn) = gBnbld * iBnc(1:mn) + ( one - gBnbld ) * efmn(1:mn)

     case( 1 ) ! Lzerovac = 1 ; 09 Mar 17;

      ;                iBns(2:mn) =                                          ofmn(2:mn) ! no blend; 27 Feb 17;
      if( NOTstellsym) iBnc(1:mn) =                                          efmn(1:mn)

      ;                iVns(2:mn) =        + iBns(2:mn) ! update vacuum field to cancel plasma field on computational boundary; 27 Feb 17;
      if( NOTstellsym) iVnc(1:mn) =        + iBnc(1:mn)

     case default ! Lzerovac; 09 Mar 17;

      FATAL( xspech, .true., invalid Lzerovac )

     end select ! end select case( Lzerovac ) ; 27 Feb 17;

    end select ! end select case( mfreeits ) ; 27 Feb 17;

   else ! else of if( bnserr.gt.gBntol ) ; 03 Apr 19;

    LContinueFreeboundaryIterations = .false.

   endif ! end of if( bnserr.gt.gBntol ) ; 24 Nov 16;

   cput = GETTIME

   if( myid.eq.0 ) then ; write(ounit,1003)
    ;                   ; write(ounit,1004) cput-cpus, nfreeboundaryiterations, mfreeits, gBntol, bnserr, cput-lastcpu
    ;                   ; write(ounit,1003)
   endif ! end of if( myid.eq.0 ) ; 24 Nov 16;

1003 format("xspech : " 10x " : ")
1004 format("xspech : "f10.2" : nfreeboundaryiterations = "i6" / "i6.5" ; gBntol ="es8.1" ; bnserr =",es12.5," ; bnorml time ="f10.2"s ;")

  endif ! end of if( LupdateBn ) ;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


!latex \subsection{output files: vector potential}
!latex \begin{enumerate}
!latex \item The vector potential is written to file using \link{ra00aa}.
!latex \end{enumerate}

  WCALL( xspech, ra00aa, ('W') ) ! this writes vector potential to file;

  if( myid.eq.0 ) then ! write restart file; note that this is inside free-boundary iteration loop; 11 Aug 14;
   WCALL( xspech, wrtend ) ! write restart file; save initial input;
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( LContinueFreeboundaryIterations .and. nfreeboundaryiterations.lt.mfreeits ) goto 9000  ! removed Lfindzero check; Loizu Dec 18;
  if( Lfreebound.eq.1 .and. First_free_bound ) goto 9000  ! going back to normal free_boundary calculation; Zhu 20190701;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! FREE-BOUNDARY ITERATIONS HAVE FINISHED;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!  SALLOCATE( gradient, (0:NGdof), zero )
!
!  lastcpu = GETTIME
!
!  LComputeDerivatives = .false.
!
!  WCALL( xspech, dforce, ( NGdof, position(0:NGdof), gradient(0:NGdof), LComputeDerivatives ) ) ! (re-)calculate Beltrami fields;
!
!  DALLOCATE(gradient)
!
!#ifdef DEBUG
!  do vvol = 1, Mvol-1
!   ; FATAL( xspech, BBe(vvol).lt.logtolerance, underflow )
!   if( Igeometry.eq.3 ) then ! include spectral constraints; 04 Dec 14;
!    ;FATAL( xspech, IIo(vvol).lt.logtolerance, underflow )
!   endif
!   if( NOTstellsym ) then
!    ;FATAL( xspech, BBo(vvol).lt.logtolerance, underflow )
!    if( Igeometry.eq.3 ) then ! include spectral constraints; 04 Dec 14;
!     FATAL( xspech, IIe(vvol).lt.logtolerance, underflow )
!    endif
!   endif
!  enddo
!#endif
!
!  if( myid.eq.0 ) then
!   cput = GETTIME
!   write(ounit,1000) cput-cpus, nfreeboundaryiterations,          ForceErr,  cput-lastcpu, "|BB|e", alog10(BBe(1:min(Mvol-1,28)))
!   if( Igeometry.ge.3 ) then ! include spectral constraints; 04 Dec 14;
!   write(ounit,1001)                                                                       "|II|o", alog10(IIo(1:min(Mvol-1,28)))
!   endif
!   if( NOTstellsym ) then
!   write(ounit,1001)                                                                       "|BB|o", alog10(BBo(1:min(Mvol-1,28)))
!   if( Igeometry.ge.3 ) then ! include spectral constraints; 04 Dec 14;
!   write(ounit,1001)                                                                       "|II|e", alog10(IIe(1:min(Mvol-1,28)))
!   endif
!   endif
!  endif
!
!2000 format("finish : ",f10.2," : finished ",i3," ; ":"|f|="es12.5" ; ":"time=",f10.2,"s ;":" log"a5,:"="28f6.2" ...")
!2001 format("finish : ", 10x ," :          ",3x," ; ":"    "  12x "   ":"     ", 10x ,"  ;":" log"a5,:"="28f6.2" ...")


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( myid.eq.0 ) then ! this is just screen diagnostics; 20 Jun 14;

   cput = GETTIME

   if( nPpts.gt.0 ) then
    write(ounit,'("xspech : ", 10x ," :")')
    write(ounit,'("xspech : ",f10.2," : myid=",i3," ; Poincare plot ; odetol="es8.1" ; nPpts="i7" ;":" nPtrj="24(i5",")" ...")') &
  cput-cpus, myid, odetol, nPpts, nPtrj(1:min(Mvol,24))
   endif

   if( Lcheck.eq.1 ) then
    write(ounit,'("xspech : ", 10x ," :")')
    write(ounit,'("xspech : ",f10.2," : myid=",i3," ; calling jo00aa; computing error in field ;")') cput-cpus, myid
   endif

  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{final diagnostics}

!latex \begin{enumerate}
!latex \item \link{sc00aa} is called to compute the covariant components of the magnetic field at the interfaces;
!latex       these are related to the singular currents;
!latex \item if \inputvar{Lcheck} $= 1$, \link{jo00aa} is called to compute the error in the Beltrami equation;
!latex \item \link{pp00aa} is called to construct the \Poincare plot;
!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
! Computes the surface current at each interface for output

  SALLOCATE( Bt00, (1:Mvol, 0:1, -1:2) , zero)

  do vvol = 1, Mvol

    LREGION(vvol)

    do iocons = 0, 1
	  if( ( Lcoordinatesingularity .and. iocons.eq.0 ) .or. ( Lvacuumregion .and. iocons.eq.1 ) ) cycle
      ! Compute covariant magnetic field at interface
      WCALL(xspech, lbpol, (vvol, Bt00(1:Mvol, 0:1, -1:2), 0, iocons) )

      ! Save covariant magnetic field at interface for output - computed in lbpol
      Btemn(1:mn, iocons, vvol) = efmn(1:mn)
      Btomn(1:mn, iocons, vvol) = ofmn(1:mn)
      Bzemn(1:mn, iocons, vvol) = cfmn(1:mn)
      Bzomn(1:mn, iocons, vvol) = sfmn(1:mn)
    enddo
  enddo

  ! Evaluate surface current
  do vvol = 1, Mvol-1
    IPDt(vvol) = pi2 * (Bt00(vvol+1, 0, 0) - Bt00(vvol, 1, 0))
  enddo

  DALLOCATE( Bt00 )

  ! Evaluate volume current
  sumI = 0
  do vvol = 1, Mvol
    Ivolume(vvol) = mu(vvol) * dtflux(vvol) * pi2 + sumI    ! factor pi2 due to normalization in preset
    sumI = Ivolume(vvol)                                    ! Sum over all volumes since this is how Ivolume is defined
  enddo


  do vvol = 1, Mvol

   LREGION(vvol)

   if( myid.eq.modulo(vvol-1,ncpu) .and. myid.lt.Mvol) then ! the following is in parallel; 20 Jun 14;

    if( .not.ImagneticOK(vvol) ) then ; cput = GETTIME ; write(ounit,1002) cput-cpus ; write(ounit,1002) cput-cpus, myid, vvol, ImagneticOK(vvol) ; cycle
    endif

    ! No need for sc00aa anymore - this is done in lbpol
    !;WCALL( xspech, sc00aa, ( vvol, Ntz                  ) ) ! compute covariant field (singular currents);

    if( Lcheck.eq.1 ) then
     WCALL( xspech, jo00aa, ( vvol, Ntz, Iquad(vvol), mn ) )
    endif

   endif ! myid.eq.modulo(vvol-1,ncpu)
  enddo ! end of do vvol = 1, Mvol; ! end of parallel diagnostics loop; 03 Apr 13;

  if( nPpts .gt.0 ) then
   WCALL( xspech, pp00aa ) ! do Poincare plots in all volumes; has its own paralellization over volumes internally
  endif

1002 format("xspech : ",f10.2," :":" myid=",i3," ; vvol=",i3," ; IBeltrami="L2" ; construction of Beltrami field failed ;")

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  do vvol = 1, Mvol ; llmodnp = modulo(vvol-1,ncpu)

#ifdef DEBUG
   FATAL( xspech, .not.allocated(Btemn), error )
   FATAL( xspech, .not.allocated(Bzemn), error )
   FATAL( xspech, .not.allocated(Btomn), error )
   FATAL( xspech, .not.allocated(Bzomn), error )
#endif

   RlBCAST( Btemn(1:mn,0:1,vvol), mn*2, llmodnp ) ! this is computed in sc00aa; 07 Dec 16;
   RlBCAST( Bzemn(1:mn,0:1,vvol), mn*2, llmodnp )
   RlBCAST( Btomn(1:mn,0:1,vvol), mn*2, llmodnp )
   RlBCAST( Bzomn(1:mn,0:1,vvol), mn*2, llmodnp )

   RlBCAST( beltramierror(vvol,1:9), 9, llmodnp ) ! this is computed in jo00aa; 21 Aug 18;

  enddo ! end of do vvol = 1, Mvol; 01 Jul 14;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{restart files}

!latex \begin{enumerate}
!latex \item \link{global}:\type{wrtend} is called to write the restart files.
!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  WCALL( xspech, write_grid ) ! write grid

  if( myid.eq.0 ) then
   WCALL( xspech, wrtend ) ! write restart file; save initial input;

   cput = GETTIME
   write(ounit,'("xspech : ", 10x ," :")')
   write(ounit,'("xspech : ",f10.2," : myid=",i3," : time="f8.2"m = "f6.2"h = "f5.2"d ;")') cput-cpus, myid, (cput-cpus) / (/ 60, 60*60, 24*60*60 /)

  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

9999 continue

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine spec

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{subroutine : ending}

!latex \begin{enumerate}
!latex \item Closes output files, writes screen summary.
!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine ending

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  use constants, only : zero

  use fileunits, only : ounit

  use inputlist, only : Wmacros, Wxspech, Ltiming

  use cputiming

  use allglobal, only : myid, cpus, mn, MPI_COMM_SPEC, ext

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS

  REAL      :: Ttotal, dcpu, ecpu
  CHARACTER :: date*8, time*10

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  cpui = GETTIME ; cpuo = cpui ! see macro expansion for begin; 11 Aug 14;

#ifdef DEBUG
  if( Wxspech ) write(ounit,'("ending : ",f10.2," : myid=",i3," ; start  ;")') cpui-cpus, myid
#endif

  cput = GETTIME

! SUMTIME !  this is expanded by Makefile, and then again by macros; do not remove;

  cput = GETTIME ; dcpu = cput-cpus

  if( Ltiming .and. myid.eq.0 ) then

   Ttotal = zero

! PRTTIME ! this is expanded by Makefile, and then again by macros; do not remove;
   write(ounit,'("ending : ",f10.2," : time spent in wrtend =",f10.2," ;")') dcpu, Twrtend ; Ttotal = Ttotal + Twrtend
   write(ounit,'("ending : ",f10.2," : time spent in readin =",f10.2," ;")') dcpu, Treadin ; Ttotal = Ttotal + Treadin

   ecpu = Ttotal-dcpu ! error in actual cpu time and calculated cpu time;  7 Mar 13;

   write(ounit,'("ending : ",f10.2," : Ttotal =",f10.2," s = "f8.2" m = "f6.2" h ; Timing Error = ",f10.2,"s = ",f10.2,"%")') &
dcpu, Ttotal / (/ 1, 60, 3600 /), ecpu, 100*ecpu/dcpu

  endif ! end of if( Ltiming .and. myid.eq.0 ) then; 01 Jul 14;

  if( myid.eq.0 ) then
   call date_and_time(date,time)
   write(ounit,'("ending : ", 10x ," : ")')
   write(ounit,1000) dcpu, myid, dcpu / (/ 1, 60, 60*60, 24*60*60 /), date(1:4), date(5:6), date(7:8), time(1:2), time(3:4), time(5:6), ext
   write(ounit,'("ending : ", 10x ," : ")')
  endif ! end of if( myid.eq.0 ) ; 14 Jan 15;

  ! wait for writing to finish
  call MPI_Barrier(MPI_COMM_SPEC, ierr)

1000 format("ending : ",f10.2," : myid=",i3," ; completion ; time=",f10.2,"s = "f8.2"m = "f6.2"h = "f5.2"d ; date= "&
  a4"/"a2"/"a2" ; time= "a2":"a2":"a2" ; ext = "a60)

end subroutine ending

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine write_hdf5
  use sphdf5,    only : hdfint, finish_outfile

  LOCALS

  call hdfint ! write final outputs to HDF5 file ! 18 Jul 14;

  call finish_outfile ! close HDF5 output file

  return

  end subroutine write_hdf5

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
