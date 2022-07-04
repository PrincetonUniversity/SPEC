!> \file
!> \brief memory management module

!> \brief allocate Beltrami matrices
!>
!> @param vvol
!> @param LcomputeDerivatives
subroutine allocate_Beltrami_matrices(vvol, LcomputeDerivatives)
    use mod_kinds, only: wp => dp
    use fileunits

    use inputlist, only: Wmemory, Wmacros

    use allglobal

    use cputiming

#ifdef OPENMP
    USE OMP_LIB
#endif
    use mpi
    implicit none
    integer :: ierr, astat, ios, nthreads, ithread
    real(wp) :: cput, cpui, cpuo = 0 ! cpu time; cpu initial; cpu old; 31 Jan 13;

    integer, intent(in) :: vvol
    LOGICAL, intent(in) :: LcomputeDerivatives
    integer :: NN

    cpui = MPI_WTIME()
    cpuo = cpui
#ifdef OPENMP
    nthreads = omp_get_max_threads()
#else
    nthreads = 1
#endif

    NN = NAdof(vvol) ! shorthand;

    if (NOTMatrixFree .or. LcomputeDerivatives) then

        allocate (dMA(0:NN, 0:NN), stat=astat)
        dMA(0:NN, 0:NN) = zero
        ! required for both plasma region and vacuum region;

        allocate (dMD(0:NN, 0:NN), stat=astat)
        dMD(0:NN, 0:NN) = zero

    else

        allocate (Adotx(0:NN), stat=astat)
        Adotx(0:NN) = zero

        allocate (Ddotx(0:NN), stat=astat)
        Ddotx(0:NN) = zero

    end if

    ! we will need the rest even with or without matrix-free

    allocate (dMB(0:NN, 0:2), stat=astat)
    dMB(0:NN, 0:2) = zero

    allocate (dMG(0:NN), stat=astat)
    dMG(0:NN) = zero

    allocate (solution(1:NN, -1:2), stat=astat)
    solution(1:NN, -1:2) = zero
    ! this will contain the vector potential from the linear solver and its derivatives;

    allocate (MBpsi(1:NN), stat=astat)
    MBpsi(1:NN) = zero

    if (LILUprecond) then

        allocate (dMAS(1:NdMASmax(vvol)), stat=astat)
        dMAS(1:NdMASmax(vvol)) = zero

        allocate (dMDS(1:NdMASmax(vvol)), stat=astat)
        dMDS(1:NdMASmax(vvol)) = zero

        allocate (idMAS(1:NN + 1), stat=astat)
        idMAS(1:NN + 1) = 0

        allocate (jdMAS(1:NdMASmax(vvol)), stat=astat)
        jdMAS(1:NdMASmax(vvol)) = 0

    end if ! if we use GMRES and ILU preconditioner

9999 continue
    cput = MPI_WTIME()
    Tmemory = Tmemory + (cput - cpuo)
    return

end subroutine allocate_Beltrami_matrices

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> \brief deallocate Beltrami matrices
!>
!> @param LcomputeDerivatives
subroutine deallocate_Beltrami_matrices(LcomputeDerivatives)
    use mod_kinds, only: wp => dp
    use fileunits

    use inputlist, only: Wmemory, Wmacros

    use allglobal

    use cputiming

#ifdef OPENMP
    USE OMP_LIB
#endif
    use mpi
    implicit none
    integer :: ierr, astat, ios, nthreads, ithread
    real(wp) :: cput, cpui, cpuo = 0 ! cpu time; cpu initial; cpu old; 31 Jan 13;

    LOGICAL, intent(in) :: LcomputeDerivatives

    cpui = MPI_WTIME()
    cpuo = cpui
#ifdef OPENMP
    nthreads = omp_get_max_threads()
#else
    nthreads = 1
#endif

    if (NOTMatrixFree .or. LcomputeDerivatives) then

        deallocate (dMA, stat=astat)

        deallocate (dMD, stat=astat)

    else

        deallocate (Adotx, stat=astat)

        deallocate (Ddotx, stat=astat)

    end if

    deallocate (dMB, stat=astat)

    deallocate (dMG, stat=astat)

    deallocate (solution, stat=astat)

    deallocate (MBpsi, stat=astat)

    if (LILUprecond) then

        deallocate (dMAS, stat=astat)

        deallocate (dMDS, stat=astat)

        deallocate (idMAS, stat=astat)

        deallocate (jdMAS, stat=astat)

    end if ! if we use GMRES and ILU preconditioner

9999 continue
    cput = MPI_WTIME()
    Tmemory = Tmemory + (cput - cpuo)
    return

end subroutine deallocate_Beltrami_matrices

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> \brief allocate geometry matrices
!>
!> @param vvol
!> @param LcomputeDerivatives
subroutine allocate_geometry_matrices(vvol, LcomputeDerivatives)
    use mod_kinds, only: wp => dp
! Allocate all geometry dependent matrices for a given ll

    use constants, only: zero

    use fileunits

    use inputlist, only: Wmemory, Wmacros, Mpol, Lrad

    use allglobal

    use cputiming

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef OPENMP
    USE OMP_LIB
#endif
    use mpi
    implicit none
    integer :: ierr, astat, ios, nthreads, ithread
    real(wp) :: cput, cpui, cpuo = 0 ! cpu time; cpu initial; cpu old; 31 Jan 13;

    integer :: vvol

    LOGICAL, intent(in) :: LcomputeDerivatives

    integer :: ll, lldof, jjdof, iidof

    cpui = MPI_WTIME()
    cpuo = cpui
#ifdef OPENMP
    nthreads = omp_get_max_threads()
#else
    nthreads = 1
#endif

    ll = Lrad(vvol)

    if (Lcoordinatesingularity) then ! different radial dof for Zernike; 02 Jul 19
        lldof = (Lrad(vvol) - mod(Lrad(vvol), 2))/2
        if (YESMatrixFree .and. .not. LcomputeDerivatives) then
            ! we only need a reduced number of terms to be computed for the preconditioner
            iidof = Mpol + 1
            jjdof = 1
        else
            ! we need full-size matrices
            iidof = mn
            jjdof = mn
        end if
    else
        lldof = Lrad(vvol)
        if (YESMatrixFree .and. .not. LcomputeDerivatives) then
            iidof = 1
            jjdof = 1
        else
            iidof = mn
            jjdof = mn
        end if
    end if

    allocate (guvijsave(1:Ntz, 1:3, 1:3, 1:Iquad(vvol)), stat=astat)
    guvijsave(1:Ntz, 1:3, 1:3, 1:Iquad(vvol)) = zero

    allocate (DToocc(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
    DToocc(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero

    allocate (TTssss(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
    TTssss(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero

    allocate (TDstsc(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
    TDstsc(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero

    allocate (TDszsc(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
    TDszsc(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero

    allocate (DDttcc(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
    DDttcc(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero

    allocate (DDtzcc(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
    DDtzcc(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero

    allocate (DDzzcc(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
    DDzzcc(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero

    allocate (Tss(0:lldof, 1:mn), stat=astat)
    Tss(0:lldof, 1:mn) = zero

    allocate (Dtc(0:lldof, 1:mn), stat=astat)
    Dtc(0:lldof, 1:mn) = zero

    allocate (Dzc(0:lldof, 1:mn), stat=astat)
    Dzc(0:lldof, 1:mn) = zero

    allocate (Ttc(0:lldof, 1:mn), stat=astat)
    Ttc(0:lldof, 1:mn) = zero

    allocate (Tzc(0:lldof, 1:mn), stat=astat)
    Tzc(0:lldof, 1:mn) = zero

    if (NOTstellsym) then

        allocate (DToocs(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        DToocs(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero

        allocate (DToosc(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        DToosc(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero

        allocate (DTooss(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        DTooss(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero

        allocate (TTsscc(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        TTsscc(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero

        allocate (TTsscs(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        TTsscs(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero

        allocate (TTsssc(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        TTsssc(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero

        allocate (TDstcc(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        TDstcc(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero

        allocate (TDstcs(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        TDstcs(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero

        allocate (TDstss(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        TDstss(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero

        allocate (TDszcc(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        TDszcc(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero

        allocate (TDszcs(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        TDszcs(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero

        allocate (TDszss(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        TDszss(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero

        allocate (DDttcs(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        DDttcs(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero

        allocate (DDttsc(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        DDttsc(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero

        allocate (DDttss(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        DDttss(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero

        allocate (DDtzcs(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        DDtzcs(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero

        allocate (DDtzsc(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        DDtzsc(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero

        allocate (DDtzss(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        DDtzss(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero

        allocate (DDzzcs(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        DDzzcs(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero

        allocate (DDzzsc(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        DDzzsc(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero

        allocate (DDzzss(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        DDzzss(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero

        allocate (Tsc(0:lldof, 1:mn), stat=astat)
        Tsc(0:lldof, 1:mn) = zero

        allocate (Dts(0:lldof, 1:mn), stat=astat)
        Dts(0:lldof, 1:mn) = zero

        allocate (Dzs(0:lldof, 1:mn), stat=astat)
        Dzs(0:lldof, 1:mn) = zero

        allocate (Tts(0:lldof, 1:mn), stat=astat)
        Tts(0:lldof, 1:mn) = zero

        allocate (Tzs(0:lldof, 1:mn), stat=astat)
        Tzs(0:lldof, 1:mn) = zero

    end if !NOTstellsym

9999 continue
    cput = MPI_WTIME()
    Tmemory = Tmemory + (cput - cpuo)
    return

end subroutine allocate_geometry_matrices

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> \brief deallocate geometry matrices
!>
!> @param LcomputeDerivatives
subroutine deallocate_geometry_matrices(LcomputeDerivatives)
    use mod_kinds, only: wp => dp
! Deallocate all geometry dependent matrices
    use constants, only: zero

    use fileunits

    use inputlist, only: Wmemory, Wmacros

    use allglobal

    use cputiming

#ifdef OPENMP
    USE OMP_LIB
#endif
    use mpi
    implicit none
    integer :: ierr, astat, ios, nthreads, ithread
    real(wp) :: cput, cpui, cpuo = 0 ! cpu time; cpu initial; cpu old; 31 Jan 13;

    LOGICAL, intent(in) :: LcomputeDerivatives

    cpui = MPI_WTIME()
    cpuo = cpui
#ifdef OPENMP
    nthreads = omp_get_max_threads()
#else
    nthreads = 1
#endif

    Lsavedguvij = .false.

    deallocate (guvijsave, stat=astat)

    deallocate (DToocc, stat=astat)

    deallocate (TTssss, stat=astat)

    deallocate (TDstsc, stat=astat)

    deallocate (TDszsc, stat=astat)

    deallocate (DDttcc, stat=astat)

    deallocate (DDtzcc, stat=astat)

    deallocate (DDzzcc, stat=astat)

    deallocate (Tss, stat=astat)

    deallocate (Dtc, stat=astat)

    deallocate (Dzc, stat=astat)

    deallocate (Ttc, stat=astat)

    deallocate (Tzc, stat=astat)

    if (NOTstellsym) then

        deallocate (DToocs, stat=astat)

        deallocate (DToosc, stat=astat)

        deallocate (DTooss, stat=astat)

        deallocate (TTsscc, stat=astat)

        deallocate (TTsscs, stat=astat)

        deallocate (TTsssc, stat=astat)

        deallocate (TDstcc, stat=astat)

        deallocate (TDstcs, stat=astat)

        deallocate (TDstss, stat=astat)

        deallocate (TDszcc, stat=astat)

        deallocate (TDszcs, stat=astat)

        deallocate (TDszss, stat=astat)

        deallocate (DDttcs, stat=astat)

        deallocate (DDttsc, stat=astat)

        deallocate (DDttss, stat=astat)

        deallocate (DDtzcs, stat=astat)

        deallocate (DDtzsc, stat=astat)

        deallocate (DDtzss, stat=astat)

        deallocate (DDzzcs, stat=astat)

        deallocate (DDzzsc, stat=astat)

        deallocate (DDzzss, stat=astat)

        deallocate (Tsc, stat=astat)

        deallocate (Dts, stat=astat)

        deallocate (Dzs, stat=astat)

        deallocate (Tts, stat=astat)

        deallocate (Tzs, stat=astat)

    end if

9999 continue
    cput = MPI_WTIME()
    Tmemory = Tmemory + (cput - cpuo)
    return

end subroutine deallocate_geometry_matrices
