!> \file
!> \brief Main program

!> \brief Main program of SPEC.
!>
!> This only calls the xpech() subroutine to do a stand-alone SPEC run.
!> @return none
program spec_main
  implicit none
  call xspech
end program spec_main

!> \brief Main subroutine of SPEC.
!>
!> This orchestrates a stand-alone SPEC run:
!> <ul>
!> <li> read the input file </li>
!> <li> solve the MRxMHD equilibrium (see spec() )</li>
!> <li> run some diagnostics on the results </li>
!> <li> write the output file(s) </li>
!> </ul>
subroutine xspech

  use numerical
  use allglobal, only: set_mpi_comm, myid, ncpu, cpus, version, MPI_COMM_SPEC, &
                       wrtend, read_inputlists_from_file, check_inputs, broadcast_inputs, skip_write, &
                       ext
  use inputlist, only: initialize_inputs, Wxspech
  use fileunits, only: ounit
  use sphdf5,    only: init_outfile, &
                       mirror_input_to_outfile, &
                       init_convergence_output, &
                       hdfint, finish_outfile, write_grid
  use cputiming, only: Txspech

  LOCALS

  CHARACTER            :: ldate*8, ltime*10, arg*100

#ifdef DEBUG
  character(len=255)   :: hostname
  integer              :: iwait, pid, status
  INTEGER, external    :: getpid, hostnm
#endif

  call MPI_INIT( ierr )

  BEGIN(xspech)

  ! set default communicator to MPI_COMM_WORLD
  call set_mpi_comm(MPI_COMM_WORLD)

  ! set initial time
  cpus = GETTIME
  cpuo = cpus

  ! explicitly enable writing of HDF5 output file
  skip_write = .false.

  ! print header: version of SPEC, compilation info, current date and time, machine precision
  cput = GETTIME
  if( myid.eq.0 ) then

    ! screen output header
    write(ounit,'("xspech : ", 10x ," : version = "F5.2)') version
! COMPILATION ! do not delete; this line is replaced (see Makefile) with a write statement identifying date, time, compilation flags, etc.;
    call date_and_time( ldate, ltime )
    write(ounit,'("xspech : ", 10x ," : ")')
    write(ounit,1000) cput-cpus, ldate(1:4), ldate(5:6), ldate(7:8), ltime(1:2), ltime(3:4), ltime(5:6), machprec, vsmall, small

    write(ounit,'("xspech : ", 10x ," : ")')
    write(ounit,'("xspech : ",f10.2," : parallelism : ncpu=",i3," ; nthreads=",i3," ;")') cput-cpus, ncpu, nthreads

    ! read command-line arguments
    call read_command_args()

    ! initialize input arrays into a default state
    call initialize_inputs()

    write(ounit,'("xspech : ", 10x ," : ")')
    write(ounit,'("xspech : ",f10.2," : begin execution ; calling global:readin ;")') cput-cpus

!> **reading input, allocating global variables**
!>
!> <ul>
!> <li> The input namelists and geometry are read in via a call to readin() .
!>       A full description of the required input is given in global.f90 . </li>
!> <li> Most internal variables, global memory etc., are allocated in preset() . </li>
!> <li> All quantities in the input file are mirrored into the output file's group \c /input . </li>
!> </ul>
    call read_inputlists_from_file()

    ! check that data from input file is within allowed ranges etc.
    call check_inputs()

  endif ! myid.eq.0

  ! broadcast input file contents
  call broadcast_inputs()

  ! initialize internal arrays based on data from input file
  call preset()

  ! initialize HDF5 library and open output file ext.h5 for writing during execution
  call init_outfile()

  ! mirror input file contents to output file
  call mirror_input_to_outfile()

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if ( myid .eq. 0 ) then ! save restart file;
    call wrtend() ! write initial restart file
  endif

!> **preparing output file group iterations**
!>
!> <ul>
!> <li> The group \c /iterations is created in the output file.
!>       This group contains the interface geometry at each iteration, which is useful for constructing movies illustrating the convergence.
!>       The data structure in use is an unlimited array of the following compound datatype:
!> ```C
!> DATATYPE  H5T_COMPOUND {
!>       H5T_NATIVE_INTEGER "nDcalls";
!>       H5T_NATIVE_DOUBLE "Energy";
!>       H5T_NATIVE_DOUBLE "ForceErr";
!>       H5T_ARRAY { [Mvol+1][mn] H5T_NATIVE_DOUBLE } "iRbc";
!>       H5T_ARRAY { [Mvol+1][mn] H5T_NATIVE_DOUBLE } "iZbs";
!>       H5T_ARRAY { [Mvol+1][mn] H5T_NATIVE_DOUBLE } "iRbs";
!>       H5T_ARRAY { [Mvol+1][mn] H5T_NATIVE_DOUBLE } "iZbc";
!> }
!> ```
!> </li>
!> </ul>
  call init_convergence_output()

!#ifdef DEBUG
!  iwait = 0; pid = getpid()
!  status = hostnm( hostname )
!  write(*,*) 'Process with PID: ', pid, 'ready to attach. Hostname: ', hostname
!  do while( iwait .EQ. 0 )
!    !wait for debugger
!  enddo
!#endif

  ! MAIN SUBROUTINE: iterate until converged or #iterations exceeds limit
  call spec()

  ! some final diagnostics: compute errors, Poincare plots, ....
  call final_diagnostics()

  ! post-processing: magnetic field evaluated on a grid
  call write_grid()

  if( myid.eq.0 ) then

!> **restart files**
!>
!> <ul>
!> <li> wrtend() is called to write the restart files. </li>
!> </ul>
    call wrtend()
  endif

  ! write final outputs to HDF5 file
  call hdfint()

  ! close HDF5 output file
  call finish_outfile()

  ! print ending info
  call ending()

  ! wait for writing to finish
  call MPI_Barrier(MPI_COMM_SPEC, ierr)

  if (myid.eq.0) then
   cput = GETTIME
   write(ounit,'("xspech : ", 10x ," :")')
   write(ounit,'("xspech : ",f10.2," : myid=",i3," : time="f8.2"m = "f6.2"h = "f5.2"d ;")') cput-cpus, myid, (cput-cpus) / (/ 60, 60*60, 24*60*60 /)
  endif

  MPIFINALIZE

  stop

1000 format("xspech : ",f10.2," : date="a4"/"a2"/"a2" , "a2":"a2":"a2" ; machine precision="es9.2" ; vsmall="es9.2" ; small="es9.2" ;")

end subroutine xspech

!> \brief Read command-line arguments; in particular, determine input file (name or extension).
!>
!> <ul>
!> <li> The input file name, \c ext , is given as the first command line input, and the input file itself is then \c ext.sp .</li>
!> <li> Alternatively, you can directly specify the input file itself as \c ext.sp .</li>
!> <li> You can also generate a template input file using \c xspec -i .</li>
!> <li> Or print help information using \c xspec -h .</li>
!> <li> Additional command line inputs recognized are:
!>      <ul>
!>      <li> \c -readin will immediately set \c Wreadin=T ; this may be over-ruled when the namelist \c screenlist is read
!>      </ul> </li>
!> </ul>
subroutine read_command_args

  use fileunits, only: ounit
  use inputlist, only: Wreadin
  use allglobal, only: cpus, myid, ext, get_hidden, MPI_COMM_SPEC, write_spec_namelist

  LOCALS

  LOGICAL              :: Lspexist
  INTEGER              :: iargc, iarg, numargs, extlen, sppos

  CHARACTER(len=100)   :: arg

  if (myid.eq.0) then

    cput = GETTIME

    ! first command-line argument is likely ext or ext.sp
    call getarg( 1, arg )

    write(ounit,'("rdcmdl : ", 10x ," : ")')
    select case (trim(arg))
    case ("", "-h", "--help")
        write(ounit,'("rdcmdl : ", 10x ," : file extension must be given as first command line argument ;")')
        write(ounit,'("rdcmdl : ", 10x ," : Usage: <mpiexec> xspec input_file <arguments>")')
        write(ounit,'("rdcmdl : ", 10x ," : Other options:")')
        write(ounit,'("rdcmdl : ", 10x ," :     -h / --help :  print help information.")')
        write(ounit,'("rdcmdl : ", 10x ," :     -i / --init :  generate a template input file.")')
        write(ounit,'("rdcmdl : ", 10x ," : Additional arguments:")')
        write(ounit,'("rdcmdl : ", 10x ," :     -readin : print debugging information during reading inputs")')
        call MPI_ABORT( MPI_COMM_SPEC, 0, ierr )
    case ("-i", "--init")
        write(ounit,'("rdcmdl : ", 10x ," : write a template input file in example.sp")')
        call write_spec_namelist()
        call MPI_ABORT( MPI_COMM_SPEC, 0, ierr )
    case default
        extlen = len_trim(arg)
        sppos = index(arg, ".sp", .true.) ! search for ".sp" from the back of ext
        if (sppos.eq.extlen-2) then       ! check if ext ends with ".sp"
            arg = arg(1:extlen-3)         ! if this is the case, remove ".sp" from end of ext
        endif
        ext = trim(arg)

        write(ounit,'("rdcmdl : ", 10x ," : ")')
        write(ounit,'("rdcmdl : ",f10.2," : ext = ",a100)') cput-cpus, ext
    end select

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    numargs = iargc()

    if( numargs.gt.1 ) then
      iarg = 1
      do while ( iarg < numargs )
        iarg = iarg + 1 ; call getarg( iarg, arg)
        select case( arg )
        case("-readin"   ) ; Wreadin = .true.
        case("-p4pg"     ) ; iarg = iarg + 1 ; call getarg( iarg, arg) ! TODO: what is this?
        case("-p4wd"     ) ; iarg = iarg + 1 ; call getarg( iarg, arg) ! TODO: what is this?
        case default       ; write(ounit,'("rdcmdl : ",f10.2," : myid=",i3," : argument not recognized ; arg = ",a100)') cput-cpus, myid, arg
        end select
      enddo
    endif

  end if ! check for myid.eq.0

end subroutine read_command_args


!> \brief This is the main "driver" for the physics part of SPEC.
!>
!> Picard iterations are performed (if in free-boundary mode)
!> and within each Picard iteration, the fixed-boundary problem
!> is solved (also iteratively).
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

  use allglobal, only : wrtend, ncpu, myid, cpus, ext, &
                        Mvol, &
                        YESstellsym, NOTstellsym, &
                        mn, im, in, &
                        Ntz, &
                        LGdof, NGdof, &
                        iRbc, iZbs, iRbs, iZbc, &
                        BBe, IIo, BBo, IIe, &
                        vvolume, &
                        Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, &
                        dtflux, dpflux, &
                        ImagneticOK, &
                        ForceErr, BnsErr,&
                        efmn, ofmn, cfmn, sfmn, &
                        iBns, iBnc, iVns, iVnc, &
                        Ate, Aze, Ato, Azo, & ! only required for debugging; 09 Mar 17;
                        nfreeboundaryiterations, &
                        beltramierror, &
                        first_free_bound, &
                        dMA, dMB, dMD, dMG, MBpsi, solution, IPDt, &
                        version, &
                        MPI_COMM_SPEC, &
                        force_final, Lhessianallocated, LocalConstraint, hessian, dBBdmp, dFFdRZ, dmupfdx, &
                        dRodR, dRodZ, dZodR, dZodZ, dRadR, dRadZ, dZadR, dZadZ, dessian, LGdof

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS

  LOGICAL              :: LComputeDerivatives, LContinueFreeboundaryIterations, exist, LupdateBn, LComputeAxis
  INTEGER              :: imn, lmn, lNfp, lim, lin, ii, ideriv, stat
  INTEGER              :: vvol, ifail, wflag, iflag, vflag
  REAL                 :: rflag, lastcpu, lRwc, lRws, lZwc, lZws, lItor, lGpol, lgBc, lgBs
  REAL,    allocatable :: position(:), gradient(:)
  CHARACTER            :: pack
  INTEGER              :: Lfindzero_old, mfreeits_old
  REAL                 :: gBnbld_old
  INTEGER              :: lnPtrj, numTrajTotal, Lfindzero_temp

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  cpuo = GETTIME

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

!> **packing geometrical degrees-of-freedom into vector**
!>
!> <ul>
!> <li> If \c NGdof.gt.0 , where \c NGdof counts the geometrical degrees-of-freedom, i.e. the \f$R_{bc}\f$, \f$Z_{bs}\f$, etc.,
!>       then packxi() is called to "pack" the geometrical degrees-of-freedom into \c position(0:NGdof) . </li>
!> </ul>

  if( NGdof.gt.0 ) then ! pack geometry into vector; 14 Jan 13;

   pack = 'P'
   LComputeAxis = .true.
   WCALL( xspech, packxi, ( NGdof, position(0:NGdof), Mvol, mn, iRbc(1:mn,0:Mvol), iZbs(1:mn,0:Mvol), &
                            iRbs(1:mn,0:Mvol), iZbc(1:mn,0:Mvol), pack, .false., LComputeAxis ) )

  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **initialize adiabatic constants**
!>
!> <ul>
!> <li> If \c Ladiabatic.eq.0 , then the "adiabatic constants" in each region, \f$P_v\f$, are calculated as
!>       \f{eqnarray}{ P_v \equiv p_v V_v^\gamma, \label{eq:adiabatic_xspech}
!>       \f}
!>       where \f$p_v\equiv\,\f$\c pressure(vvol) , the volume \f$V_v\f$ of each region is computed by volume() ,
!>       and the adiabatic index \f$\gamma\equiv\,\f$\c gamma . </li>
!> </ul>

  do vvol = 1, Mvol

   LREGION(vvol)
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

!> **solving force-balance**
!>
!> <ul>
!> <li> If there are geometrical degress of freedom, i.e. if \c NGdof.gt.0 , then
!>      <ul>
!>      <li> \todo If \c Lminimize.eq.1 , call pc00aa() to find minimum of energy functional
!>              using quasi-Newton, preconditioned conjugate gradient method, \c E04DGF
!>
!>      </li>
!>      <li> If \c Lfindzero.gt.0 , call newton() to find extremum of constrained energy functional using a Newton method, \c C05PDF . </li>
!>      </ul> </li>
!> </ul>

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

!> **post diagnostics**
!>
!> <ul>
!> <li> The pressure is computed from the adiabatic constants from Eqn.\f$(\ref{eq:adiabatic_xspech})\f$, i.e. \f$p=P/V^\gamma\f$. </li>
!> <li> The Beltrami/vacuum fields in each region are re-calculated using dforce() . </li>
!> <li> If \c Lcheck.eq.5 \c .or. \c LHevalues \c .or. \c LHevectors \c .or. \c Lperturbed.eq.1 , then the force-gradient matrix is examined using hesian() . </li>
!> </ul>

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

  if(Lcheck.eq.7) then
    SALLOCATE( force_final, (0:NGdof), zero )
    SALLOCATE( hessian, (1:NGdof,1:NGdof), zero )
    SALLOCATE( dessian, (1:NGdof,1:LGdof), zero )
    SALLOCATE( dFFdRZ, (1:LGdof,0:1,1:LGdof,0:1,1:Mvol), zero )
    SALLOCATE( dBBdmp, (1:LGdof,1:Mvol,0:1,1:2), zero )
    if( LocalConstraint ) then
      SALLOCATE( dmupfdx, (1:Mvol,    1:1,1:2,1:LGdof,0:1), zero )
    else
      SALLOCATE( dmupfdx, (1:Mvol, 1:Mvol-1,1:2,1:LGdof,1), zero ) ! TODO change the format to put vvol in last index position...
    endif
    Lhessianallocated = .true.
    SALLOCATE( dRodR, (1:Ntz,0:3,1:mn), zero ) ! calculated in rzaxis; 19 Sep 16;
    SALLOCATE( dRodZ, (1:Ntz,0:3,1:mn), zero )
    SALLOCATE( dZodR, (1:Ntz,0:3,1:mn), zero )
    SALLOCATE( dZodZ, (1:Ntz,0:3,1:mn), zero )
    SALLOCATE( dRadR, (1:mn,0:1,0:1,1:mn), zero ) ! calculated in rzaxis; 19 Sep 16;
    SALLOCATE( dRadZ, (1:mn,0:1,0:1,1:mn), zero )
    SALLOCATE( dZadR, (1:mn,0:1,0:1,1:mn), zero )
    SALLOCATE( dZadZ, (1:mn,0:1,0:1,1:mn), zero )

    LComputeDerivatives = .true.

    LComputeAxis = .true.
    Lfindzero_temp = Lfindzero
    Lfindzero = 2 ! Necessary to trigger axis recomputation in packxi
    WCALL( xspech, dforce, ( NGdof, position(0:NGdof), force_final(0:NGdof), LComputeDerivatives, LComputeAxis) )
    Lfindzero = Lfindzero_temp

  else
    SALLOCATE( force_final, (0:NGdof), zero )

    LComputeDerivatives = .false.
    LComputeAxis = .true.

    WCALL( xspech, dforce, ( NGdof, position(0:NGdof), force_final(0:NGdof), LComputeDerivatives, LComputeAxis) )
  end if
  

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
    write(ounit,'("xspech : ", 10x ," : -------------------Stability Evaluations------------------ ")')
    write(ounit,'("xspech : ",f10.2," : myid=",i3," ; calling hessian; see .ext.hessian.myid ;")') cput-cpus, myid
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

!> **free-boundary: re-computing normal field**
!>
!> <ul>
!> <li> If \c Lfreebound.eq.1 and \c Lfindzero.gt.0 and \c mfreeits.ne.0 ,
!>       then the magnetic field at the computational boundary produced by the plasma currents is computed using bnorml() . </li>
!> <li> The "new" magnetic field at the computational boundary produced by the plasma currents is updated using a Picard scheme:
!>       \f{eqnarray}{ \texttt{Bns}_i^{j} = \lambda \, \texttt{Bns}_i^{j-1} + (1-\lambda) \texttt{Bns}_i, \label{eq:blending_xspech}
!>       \f}
!>       where \f$j\f$ labels free-boundary iterations, the "blending parameter" is \f$\lambda\equiv\,\f$\c gBnbld ,
!>       and \c Bns \f$_i\f$ is computed by virtual casing.
!>       The subscript "$i$" labels Fourier harmonics.  </li>
!> <li> If the new (unblended) normal field is _not_ sufficiently close to the old normal field, as quantified by \c gBntol ,
!>       then the free-boundary iterations continue.
!>       This is quantified by
!>       \f{eqnarray}{ \sum_i | \texttt{Bns}_i^{j-1} - \texttt{Bns}_i | / N, \label{eq:gBntol_xspech}
!>       \f}
!>       where \f$N\f$ is the total number of Fourier harmonics. </li>
!> <li> There are several choices that are available:
!>       <ul>
!>       <li> if \c mfreeits=-2 : the vacuum magnetic field
!>             (really, the normal component of the field produced by the external currents at the computational boundary)
!>             required to hold the given equlibrium is written to file.
!>             This information is required as input by FOCUS \cite y2017_zhu
!>             for example. (This option probably needs to revised.) </li>
!>       <li> if \c mfreeits=-1 : after the plasma field is computed by virtual casing,
!>             the vacuum magnetic field is set to exactly balance the plasma field
!>             (again, we are really talking about the normal component at the computational boundary.)
!>             This will ensure that the computational boundary itself if a flux surface of the total magnetic field. </li>
!>       <li> if \c mfreeits=0 : the plasma field at the computational boundary is not updated; no "free-boundary" iterations take place. </li>
!>       <li> if \c mfreeits>0 : the plasma field at the computational boundary is updated according to the above blending Eqn.\f$(\ref{eq:blending_xspech})\f$,
!>             and the free-boundary iterations will continue until either the tolerance condition is met (see \c gBntol and Eqn.\f$(\ref{eq:gBntol_xspech})\f$)
!>             or the maximum number of free-boundary iterations, namely \c mfreeits , is reached.
!>             For this case, \c Lzerovac is relevant:
!>             if \c Lzerovac=1 , then the vacuum field is set equal to the normal field at every iteration,
!>             which results in the computational boundary being a flux surface.
!>             (I am not sure if this is identical to setting \c mfreeits=-1 ; the logic etc. needs to be revised.) </li>
!>       </ul> </li>
!> </ul>

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

!> **output files: vector potential**
!>
!> <ul>
!> <li> The vector potential is written to file using ra00aa() . </li>
!> </ul>

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

end subroutine spec

!> \brief Final diagnostics
!>
!> <ul>
!> <li> sc00aa() is called to compute the covariant components of the magnetic field at the interfaces;
!>      these are related to the singular currents </li>
!> <li> if \c Lcheck=1 , jo00aa() is called to compute the error in the Beltrami equation </li>
!> <li> pp00aa() is called to construct the Poincare plot by field-line following. </li>
!> </ul>
subroutine final_diagnostics

  use inputlist, only: nPtrj, nPpts, Igeometry, Lcheck, Nvol, odetol, &
                       Isurf, Ivolume, mu, Wmacros, Ltransform, Lsvdiota
  use fileunits, only: ounit
  use constants, only: zero
  use allglobal, only: pi2, myid, ncpu, MPI_COMM_SPEC, cpus, Mvol, Ntz, mn, &
                       beltramierror, Lcoordinatesingularity, &
                       Lplasmaregion, Lvacuumregion, &
                       Btemn, Bzemn, Btomn, Bzomn, &
                       efmn, ofmn, cfmn, sfmn, &
                       IPDt, ImagneticOK, dtflux, Iquad, lmns, Nt, Nz, diotadxup, &
                       IsMyVolume, IsMyVolumeValue, WhichCpuID, &
                       dlambdaout, diotadxup


  LOCALS

  integer              :: iocons, llmodnp, vvol, iflag, cpu_id
  real                 :: sumI
  REAL,    allocatable :: Bt00(:,:,:)
  REAL                 :: work(0:1,-1:2) 



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


! Evaluate rotational transform and straight field line coordinate transformation
if( Ltransform ) then

  FATAL(xspech, Lsvdiota.ne.1, Lsvdiota needs to be one for s.f.l transformation)

  do vvol=1,Mvol
    call brcast(vvol)
  enddo

  iflag = -1
  do vvol = 1, Mvol
    call IsMyVolume(vvol)
    if (IsMyVolumeValue.eq.0) then
      cycle
    elseif (IsMyVolumeValue.eq.-1) then
      FATAL( xspech, .true., Unassociated volume )
    endif

    LREGION( vvol )

    call tr00ab( vvol, mn, lmns, Nt, Nz, iflag, diotadxup(0:1,-1:2, vvol) ) ! stores lambda in a global variable.
  enddo

  ! Broadcast
  do vvol = 1, Mvol
    call WhichCpuID( vvol, cpu_id )
    RlBCAST( diotadxup(0:1,-1:2,vvol), 8, cpu_id  )
    RlBCAST( dlambdaout(1:lmns, vvol, 0:1), 2*lmns, cpu_id  )
  enddo

endif


  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
! Computes the surface current at each interface for output

  SALLOCATE( Bt00, (1:Mvol, 0:1, -1:2) , zero)

  do vvol = 1, Mvol

    LREGION(vvol)

    do iocons = 0, 1
	  if( Lcoordinatesingularity .and. iocons.eq.0 ) cycle
          if( vvol.eq.Nvol+1 .and. iocons.eq.1 ) cycle
      ! Compute covariant magnetic field at interface
      call lbpol(vvol, Bt00(1:Mvol, 0:1, -1:2), 0, iocons)

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

  ! screen info about diagnostics; 20 Jun 14;
  if (myid.eq.0) then
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

  do vvol = 1, Mvol

   LREGION(vvol)

   if( myid.eq.modulo(vvol-1,ncpu) .and. myid.lt.Mvol) then ! the following is in parallel; 20 Jun 14;

    if( .not.ImagneticOK(vvol) ) then ; cput = GETTIME ; write(ounit,1002) cput-cpus ; write(ounit,1002) cput-cpus, myid, vvol, ImagneticOK(vvol) ; cycle
    endif

    ! No need for sc00aa anymore - this is done in lbpol
    !;WCALL( xspech, sc00aa, ( vvol, Ntz                  ) ) ! compute covariant field (singular currents);

    if( Lcheck.eq.1 ) then
     call jo00aa( vvol, Ntz, Iquad(vvol), mn )
    endif

   endif ! myid.eq.modulo(vvol-1,ncpu)
  enddo ! end of do vvol = 1, Mvol; ! end of parallel diagnostics loop; 03 Apr 13;

  if( nPpts .gt.0 ) then
    call pp00aa() ! do Poincare plots in all volumes; has its own paralellization over volumes internally
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

   RlBCAST( Btemn(1:mn,0:1,vvol), mn*2, llmodnp ) ! this is computed in lbpol; 07 Dec 16;
   RlBCAST( Bzemn(1:mn,0:1,vvol), mn*2, llmodnp )
   RlBCAST( Btomn(1:mn,0:1,vvol), mn*2, llmodnp )
   RlBCAST( Bzomn(1:mn,0:1,vvol), mn*2, llmodnp )

   RlBCAST( beltramierror(vvol,1:9), 9, llmodnp ) ! this is computed in jo00aa; 21 Aug 18;

  enddo ! end of do vvol = 1, Mvol; 01 Jul 14;

end subroutine final_diagnostics

!> \brief Closes output files, writes screen summary.
!>
subroutine ending

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

1000 format("ending : ",f10.2," : myid=",i3," ; completion ; time=",f10.2,"s = "f8.2"m = "f6.2"h = "f5.2"d ; date= "&
  a4"/"a2"/"a2" ; time= "a2":"a2":"a2" ; ext = "a60)

end subroutine ending

