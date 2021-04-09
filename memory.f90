!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
! memory management module
! by Zhisong Qu 07 JUL 20
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
subroutine allocate_Beltrami_matrices(vvol, LcomputeDerivatives)

  use fileunits

  use inputlist, only: Wmemory, Wmacros, Lvcvacuum

  use allglobal

  use cputiming

  LOCALS

  INTEGER, intent(in) :: vvol
  LOGICAL, intent(in) :: LcomputeDerivatives
  INTEGER     :: NN

  BEGIN(memory)

  NN = NAdof(vvol) ! shorthand;

  if (NOTMatrixFree .or. LcomputeDerivatives) then
      SALLOCATE( dMA, (0:NN,0:NN), zero ) ! required for both plasma region and vacuum region;
      SALLOCATE( dMD, (0:NN,0:NN), zero )
  else
      SALLOCATE( Adotx, (0:NN), zero)
      SALLOCATE( Ddotx, (0:NN), zero)
  endif

  ! we will need the rest even with or without matrix-free
  SALLOCATE( dMB, (0:NN,0: 2), zero )
  SALLOCATE( dMG, (0:NN     ), zero )  

  if( Lvacuumregion .and. Lvcvacuum.eq.1 ) then
  SALLOCATE( dVC, (0:NN,0:NN), zero ) ! 03/03/21 ;
  endif

  SALLOCATE( solution, (1:NN,-1:2), zero ) ! this will contain the vector potential from the linear solver and its derivatives;

  SALLOCATE( MBpsi, (1:NN), zero )

  if (LILUprecond) then
      SALLOCATE( dMAS, (1:NdMASmax(vvol)), zero)
      SALLOCATE( dMDS, (1:NdMASmax(vvol)), zero)
      SALLOCATE( idMAS, (1:NN+1), 0)
      SALLOCATE( jdMAS, (1:NdMASmax(vvol)), 0)
  endif ! if we use GMRES and ILU preconditioner

  RETURN(memory)

end subroutine allocate_Beltrami_matrices

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine deallocate_Beltrami_matrices(LcomputeDerivatives)

  use fileunits

  use inputlist, only: Wmemory, Wmacros, Lvcvacuum

  use allglobal

  use cputiming

  LOCALS

  LOGICAL, intent(in) :: LcomputeDerivatives

  BEGIN(memory)

  if (NOTMatrixFree .or. LcomputeDerivatives) then
    DALLOCATE(dMA)
    DALLOCATE(dMD)
  else
    DALLOCATE(Adotx)
    DALLOCATE(Ddotx)
  endif

  DALLOCATE(dMB)

  DALLOCATE(dMG)

  if( Lvacuumregion .and. Lvcvacuum.eq.1 ) then
  DALLOCATE(dVC) ! 03/03/21 ;
  endif
  
  DALLOCATE(solution)
   
  DALLOCATE(MBpsi)
   
  if (LILUprecond) then
    DALLOCATE(dMAS)
    DALLOCATE(dMDS)
    DALLOCATE(idMAS)
    DALLOCATE(jdMAS)
  endif ! if we use GMRES and ILU preconditioner

   RETURN(memory)

end subroutine deallocate_Beltrami_matrices

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine allocate_geometry_matrices(vvol, LcomputeDerivatives)

! Allocate all geometry dependent matrices for a given ll

  use constants, only: zero

  use fileunits

  use inputlist, only:  Wmemory, Wmacros, Mpol, Lrad

  use allglobal

  use cputiming

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS

  INTEGER         :: vvol

  LOGICAL, intent(in) :: LcomputeDerivatives

  INTEGER         :: ll, lldof, jjdof, iidof

  BEGIN(memory)

  ll = Lrad(vvol)

  if (Lcoordinatesingularity) then ! different radial dof for Zernike; 02 Jul 19
    lldof = (Lrad(vvol) - mod(Lrad(vvol),2)) / 2
    if (YESMatrixFree .and. .not. LcomputeDerivatives) then
      ! we only need a reduced number of terms to be computed for the preconditioner
      iidof = Mpol + 1
      jjdof = 1
    else
      ! we need full-size matrices
      iidof = mn
      jjdof = mn
    endif
  else
    lldof = Lrad(vvol)
    if (YESMatrixFree .and. .not. LcomputeDerivatives) then
      iidof = 1
      jjdof = 1
    else
      iidof = mn
      jjdof = mn
    endif
  end if

  SALLOCATE( guvijsave, (1:Ntz,1:3,1:3,1:Iquad(vvol)), zero)

  SALLOCATE( DToocc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
  SALLOCATE( TTssss, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
  SALLOCATE( TDstsc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
  SALLOCATE( TDszsc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
  SALLOCATE( DDttcc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
  SALLOCATE( DDtzcc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
  SALLOCATE( DDzzcc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )

  SALLOCATE( Tss, (0:lldof,1:mn), zero )
  SALLOCATE( Dtc, (0:lldof,1:mn), zero )
  SALLOCATE( Dzc, (0:lldof,1:mn), zero )
  SALLOCATE( Ttc, (0:lldof,1:mn), zero )
  SALLOCATE( Tzc, (0:lldof,1:mn), zero )

  if (NOTstellsym) then

    SALLOCATE( DToocs, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( DToosc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( DTooss, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )

    SALLOCATE( TTsscc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( TTsscs, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( TTsssc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )

    SALLOCATE( TDstcc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( TDstcs, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( TDstss, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )

    SALLOCATE( TDszcc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( TDszcs, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( TDszss, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )

    SALLOCATE( DDttcs, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( DDttsc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( DDttss, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )

    SALLOCATE( DDtzcs, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( DDtzsc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( DDtzss, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )

    SALLOCATE( DDzzcs, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( DDzzsc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( DDzzss, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )

    SALLOCATE( Tsc, (0:lldof,1:mn), zero )
    SALLOCATE( Dts, (0:lldof,1:mn), zero )
    SALLOCATE( Dzs, (0:lldof,1:mn), zero )
    SALLOCATE( Tts, (0:lldof,1:mn), zero )
    SALLOCATE( Tzs, (0:lldof,1:mn), zero )

  end if !NOTstellsym

  RETURN(memory)

end subroutine allocate_geometry_matrices

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine deallocate_geometry_matrices(LcomputeDerivatives)

! Deallocate all geometry dependent matrices
  use constants, only: zero

  use fileunits

  use inputlist, only: Wmemory, Wmacros

  use allglobal

  use cputiming

  LOCALS

  LOGICAL, intent(in) :: LcomputeDerivatives

  BEGIN(memory)

  Lsavedguvij = .false.
  DALLOCATE(guvijsave)

  DALLOCATE(DToocc)
  DALLOCATE(TTssss)
  DALLOCATE(TDstsc)
  DALLOCATE(TDszsc)
  DALLOCATE(DDttcc)
  DALLOCATE(DDtzcc)
  DALLOCATE(DDzzcc)

  DALLOCATE(Tss)
  DALLOCATE(Dtc)
  DALLOCATE(Dzc)
  DALLOCATE(Ttc)
  DALLOCATE(Tzc)

  if (NOTstellsym) then

    DALLOCATE(DToocs)
    DALLOCATE(DToosc)
    DALLOCATE(DTooss)

    DALLOCATE(TTsscc)
    DALLOCATE(TTsscs)
    DALLOCATE(TTsssc)

    DALLOCATE(TDstcc)
    DALLOCATE(TDstcs)
    DALLOCATE(TDstss)

    DALLOCATE(TDszcc)
    DALLOCATE(TDszcs)
    DALLOCATE(TDszss)

    DALLOCATE(DDttcs)
    DALLOCATE(DDttsc)
    DALLOCATE(DDttss)

    DALLOCATE(DDtzcs)
    DALLOCATE(DDtzsc)
    DALLOCATE(DDtzss)

    DALLOCATE(DDzzcs)
    DALLOCATE(DDzzsc)
    DALLOCATE(DDzzss)
       
    DALLOCATE(Tsc)
    DALLOCATE(Dts)
    DALLOCATE(Dzs)
    DALLOCATE(Tts)
    DALLOCATE(Tzs)

  endif

  RETURN(memory)

end subroutine deallocate_geometry_matrices
