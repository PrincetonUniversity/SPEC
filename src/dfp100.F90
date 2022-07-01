!> \file
!> \brief Split the work between MPI nodes and evaluate the global constraint

!latex \calledby{\link{dforce}} \\

!latex \calls{\link{ma02aa} and
!latex        \link{lbpol}}


!latex \tableofcontents

!latex \subsection{Local constraint}
!latex In the case of a local constraint, each MPI node calls \link{ma02aa} to obtain the field consistent with the local constraint

!latex \subsection{Global constraint}
!latex In the case of a global constraint,

!latex \begin{enumerate}
!latex \item The MPI node $0$ broadcast the value of \internal{iflag}. \internal{iflag} is equal to $5$ if the constraint has been
!latex          matched, and all nodes exit the iteration process. Otherwise, an additional iteration occurs (point 2 and below).
!latex \item The MPI node $0$ boradcast the new values of $\psi_p$.
!latex \item Each node calls \link{ma02aa} to obtain the field for the given value of $\mu$ and $\psi_p$
!latex \item Each node calls \link{lbpol} to evalue the poloidal magnetic field Fourier coefficients $B_{\theta,e,0,0}$ at the inner and outer interface
!latex          (this is only for a toroidal current constraint. Other global constraints would require different parameters)
!latex \item Each node communicate to the MPI node $0$ the parameters necessary for the evaluation of the constraint
!latex \item The node $0$ evaluate the global constraint and sets the flag \internal{IconstraintOK}
!latex \end{enumerate}

!latex \subsection{Evaluation of constraint derivatives}
!latex Not implemented for now.


!> \brief Split the work between MPI nodes and evaluate the global constraint
!>
!> @param Ndofgl
!> @param x
!> @param Fvec
!> @param LComputeDerivatives
subroutine dfp100(Ndofgl, x, Fvec, LComputeDerivatives)
  use mod_kinds, only: wp => dp
  use constants, only : zero, half, one, two, pi2, pi, mu0

  use fileunits, only : ounit

  use inputlist, only : Wmacros, Wdfp100, Igeometry, Nvol, Lrad, Isurf, &
                        Lconstraint, Lfreebound, curpol

  use cputiming, only : Tdfp100

  use allglobal, only : ncpu, myid, cpus, MPI_COMM_SPEC, &
                        ImagneticOK, NAdof, mn, &
                        Mvol, Iquad, &
                        dBdX, &
                        Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, Localconstraint, &
                        IPDt, IPDtdPf, xoffset, dpflux, &
                        IsMyVolume, IsMyVolumeValue, WhichCpuID, &
                        IconstraintOK, &
                        DToocc, DToocs, DToosc, DTooss, &
                        TTsscc, TTsscs, TTsssc, TTssss, &
                        TDstcc, TDstcs, TDstsc, TDstss, &
                        TDszcc, TDszcs, TDszsc, TDszss, &
                        DDttcc, DDttcs, DDttsc, DDttss, &
                        DDtzcc, DDtzcs, DDtzsc, DDtzss, &
                        DDzzcc, DDzzcs, DDzzsc, DDzzss, &
                        dMA, dMB, dMD, dMG, MBpsi, solution, &
                        Nt, Nz, LILUprecond, Lsavedguvij, NOTMatrixFree, guvijsave, izbs


#ifdef OPENMP
  USE OMP_LIB
#endif
  use mpi
  implicit none
  integer   :: ierr, astat, ios, nthreads, ithread
  real(wp)      :: cput, cpui, cpuo=0 ! cpu time; cpu initial; cpu old; 31 Jan 13;

  !------
  ! vvol:                       loop index on volumes
  ! Ndofgl:                     Input parameter necessary for the use of hybrd1. Unused otherwise.
  ! iflag:                      Flag changed by hybrd1
  ! cpu_send_one, cpu_send_two: CPU IDs, used for MPI communications
  ! status:                     MPI status
  ! Fvec:                       Global constraint values
  ! x:                          Degrees of freedom of hybrd1. For now contains only the poloidal flux

  integer              :: vvol, Ndofgl, iflag, cpu_send_one, cpu_send_two, ll, NN, ideriv, iocons
  integer              :: status(MPI_STATUS_SIZE), request1, request2
  real(wp)                 :: Fvec(1:Ndofgl), x(1:Mvol-1), Bt00(1:Mvol, 0:1, -1:2), ldItGp(0:1, -1:2)
  LOGICAL              :: LComputeDerivatives
  integer              :: deriv, Lcurvature




  cpui = MPI_WTIME()
  cpuo = cpui
#ifdef OPENMP
  nthreads = omp_get_max_threads()
#else
  nthreads = 1
#endif


  dpflux(2:Mvol) = x - xoffset

  ! Now each CPU perform the calculation in its volume(s)
  do vvol = 1, Mvol


   if( Igeometry.eq.1 .or. vvol.gt.1 ) then ; Lcoordinatesingularity = .false.
   else                                   ; Lcoordinatesingularity = .true.
   endif

   if( vvol.le.Nvol ) then ; Lplasmaregion = .true.
   else                  ; Lplasmaregion = .false.
   endif

   Lvacuumregion = .not.Lplasmaregion
 ! assigns Lcoordinatesingularity, Lplasmaregion, etc. ;
    ImagneticOK(vvol) = .false.
    !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    ! Determines if this volume vvol should be computed by this thread.
    call IsMyVolume(vvol)

    if( IsMyVolumeValue .EQ. 0 ) then
      cycle
    else if( IsMyVolumeValue .EQ. -1) then

   if( .true. ) then
     write(6,'("dfp100 :      fatal : myid=",i3," ; .true. ; Unassociated volume;")') myid
     call MPI_ABORT( MPI_COMM_SPEC, 1, ierr )
     stop "dfp100 : .true. : Unassociated volume ;"
    endif

    endif

    NN = NAdof(vvol)
    ll = Lrad(vvol)

    call allocate_geometry_matrices(vvol, LcomputeDerivatives)
    call allocate_Beltrami_matrices(vvol, LcomputeDerivatives)
    call intghs_workspace_init(vvol)

    dBdX%L = .false. ! No need to take derivatives of matrices w.r.t geometry.

    ideriv = 0 ; Lcurvature = 1

   cput = MPI_WTIME()
   Tdfp100 = Tdfp100 + ( cput-cpuo )
   call compute_guvijsave(Iquad(vvol), vvol, ideriv, Lcurvature)
   cpuo = MPI_WTIME()

    Lsavedguvij = .true.

  ! we need to construct the preconditioner if needed
    if (LILUprecond) then

   cput = MPI_WTIME()
   Tdfp100 = Tdfp100 + ( cput-cpuo )
   call spsint( Iquad(vvol), mn, vvol, ll )
   cpuo = MPI_WTIME()


   cput = MPI_WTIME()
   Tdfp100 = Tdfp100 + ( cput-cpuo )
   call spsmat( vvol, mn, ll)
   cpuo = MPI_WTIME()

    endif

    if (NOTMatrixFree) then ! construct Beltrami matrix

   cput = MPI_WTIME()
   Tdfp100 = Tdfp100 + ( cput-cpuo )
   call ma00aa( Iquad(vvol), mn, vvol, ll )
   cpuo = MPI_WTIME()
 ! compute volume integrals of metric elements;

   cput = MPI_WTIME()
   Tdfp100 = Tdfp100 + ( cput-cpuo )
   call matrix( vvol, mn, ll )
   cpuo = MPI_WTIME()

    else ! matrix free, so construct something else
      ! we will still need to construct the dMB and dMG matrix

   cput = MPI_WTIME()
   Tdfp100 = Tdfp100 + ( cput-cpuo )
   call matrixBG( vvol, mn, ll )
   cpuo = MPI_WTIME()

    endif


   cput = MPI_WTIME()
   Tdfp100 = Tdfp100 + ( cput-cpuo )
   call ma02aa( vvol, NN )
   cpuo = MPI_WTIME()


    Lsavedguvij = .false.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    call intghs_workspace_destroy()
    call deallocate_Beltrami_matrices(LcomputeDerivatives)
    call deallocate_geometry_matrices(LcomputeDerivatives)

    ! Compute relevant local quantities for the evaluation of global constraints
    if( Lconstraint.EQ.3 ) then
      do iocons=0,1
        if( ( Lcoordinatesingularity .and. iocons.eq.0 ) .or. ( Lvacuumregion .and. iocons.eq.1 ) ) cycle

        ideriv = 0

   cput = MPI_WTIME()
   Tdfp100 = Tdfp100 + ( cput-cpuo )
   call lbpol(vvol, Bt00(1:Mvol, 0:1, -1:2), ideriv, iocons)
   cpuo = MPI_WTIME()
 !Compute field at interface for global constraint

        ideriv = 2

   cput = MPI_WTIME()
   Tdfp100 = Tdfp100 + ( cput-cpuo )
   call lbpol(vvol, Bt00(1:Mvol, 0:1, -1:2), ideriv, iocons)
   cpuo = MPI_WTIME()
 !Compute field at interface for global constraint, d w.r.t. pflux
      enddo

      if( Lvacuumregion ) then

#ifdef DEBUG

   if( vvol.ne.Mvol ) then
     write(6,'("dfp100 :      fatal : myid=",i3," ; vvol.ne.Mvol ; Incorrect vvol in last volume;")') myid
     call MPI_ABORT( MPI_COMM_SPEC, 1, ierr )
     stop "dfp100 : vvol.ne.Mvol : Incorrect vvol in last volume ;"
    endif

#endif

        ideriv=1 ! derivatives of Btheta w.r.t tflux
        iocons=0 ! Only need inner side of volume derivatives

   cput = MPI_WTIME()
   Tdfp100 = Tdfp100 + ( cput-cpuo )
   call lbpol(Mvol, Bt00(1:Mvol, 0:1, -1:2), ideriv, iocons)
   cpuo = MPI_WTIME()


        iflag=2 ! derivatives of poloidal linking current w.r.t geometry not required

   cput = MPI_WTIME()
   Tdfp100 = Tdfp100 + ( cput-cpuo )
   call curent(Mvol, mn, Nt, Nz, iflag, ldItGp(0:1,-1:2) )
   cpuo = MPI_WTIME()

      endif
    endif
  enddo


  ! Evaluation of global constraint and communications
  if( .not.LocalConstraint ) then

    select case (Lconstraint)

      ! Case 3: toroidal current constraint
      case( 3 )

        ! Compute IPDt on each interface.
        do vvol = 1, Mvol-1

          ! --------------------------------------------------------------------------------------------------
          !                                                                     MPI COMMUNICATIONS
          call WhichCpuID(vvol  , cpu_send_one)
          call WhichCpuID(vvol+1, cpu_send_two)

          ! Broadcast magnetic field at the interface.

   call MPI_BCAST(Bt00(vvol  , 1, 0),1,MPI_DOUBLE_PRECISION,cpu_send_one,MPI_COMM_SPEC,ierr)


   call MPI_BCAST(Bt00(vvol+1, 0, 0),1,MPI_DOUBLE_PRECISION,cpu_send_two,MPI_COMM_SPEC,ierr)


   call MPI_BCAST(Bt00(vvol  , 1, 2),1,MPI_DOUBLE_PRECISION,cpu_send_one,MPI_COMM_SPEC,ierr)


   call MPI_BCAST(Bt00(vvol+1, 0, 2),1,MPI_DOUBLE_PRECISION,cpu_send_two,MPI_COMM_SPEC,ierr)


          ! Evaluate surface current
          IPDt(vvol) = pi2 * (Bt00(vvol+1, 0, 0) - Bt00(vvol, 1, 0))

          ! their derivatives
          IPDtdPf(vvol,vvol) = pi2 * Bt00(vvol+1, 0, 2)
          if (vvol .ne. 1) IPDtdPf(vvol,vvol-1) = -pi2 * Bt00(vvol, 1, 2)
        enddo

        ! Compute the constraint and store it in Fvec.
        if( myid.EQ.0 ) then
            Fvec(1:Mvol-1) = IPDt(1:Mvol-1) - Isurf(1:Mvol-1)
        endif

        ! Compute poloidal linking current constraint as well in case of free boundary computation
        if ( Lfreebound.eq.1 ) then

          ! Communicate additional derivatives
          call WhichCpuID(Mvol, cpu_send_one)

   call MPI_BCAST(ldItGp(0:1, -1:2),8,MPI_DOUBLE_PRECISION,cpu_send_one ,MPI_COMM_SPEC,ierr)


   call MPI_BCAST(Bt00(Mvol, 0:1, 1),2,MPI_DOUBLE_PRECISION,cpu_send_one ,MPI_COMM_SPEC,ierr)


          ! Complete output: RHS
          Fvec(Mvol    ) = ldItGp(1, 0) - curpol

          ! Complete output: LHS
          IPDtdPf(Mvol-1, Mvol  ) = pi2 * Bt00(Mvol, 0, 1)
          IPDtdPf(Mvol  , Mvol-1) = ldItGp(1, 2)
          IPDtdPf(Mvol  , Mvol  ) = ldItGp(1, 1)
        endif


      case default

   if( .true. ) then
     write(6,'("dfp100 :      fatal : myid=",i3," ; .true. ; Unaccepted value for Lconstraint;")') myid
     call MPI_ABORT( MPI_COMM_SPEC, 1, ierr )
     stop "dfp100 : .true. : Unaccepted value for Lconstraint ;"
    endif

    end select
  endif


9999 continue
  cput = MPI_WTIME()
  Tdfp100 = Tdfp100 + ( cput-cpuo )
  return


end subroutine dfp100
