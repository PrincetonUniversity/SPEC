! AUTO-GENERATED; DO NOT COMMIT CHANGES TO THIS FILE !
! auto-created by a user called 'jonathan' on a machine called 'Nebuchadnezaar' at 26/11/2019 15:27:37 UTC
module read_spec
TYPE physics
    INTEGER :: Igeometry
    INTEGER :: Istellsym
    INTEGER :: Lfreebound
    DOUBLE PRECISION :: phiedge
    DOUBLE PRECISION :: curtor
    DOUBLE PRECISION :: curpol
    DOUBLE PRECISION :: gamma
    INTEGER :: Nfp
    INTEGER :: Nvol
    INTEGER :: Mpol
    INTEGER :: Ntor
    INTEGER, ALLOCATABLE :: Lrad(:)
    INTEGER :: Lconstraint
    DOUBLE PRECISION, ALLOCATABLE :: tflux(:)
    DOUBLE PRECISION, ALLOCATABLE :: pflux(:)
    DOUBLE PRECISION, ALLOCATABLE :: helicity(:)
    DOUBLE PRECISION :: pscale
    DOUBLE PRECISION, ALLOCATABLE :: pressure(:)
    INTEGER :: Ladiabatic
    DOUBLE PRECISION, ALLOCATABLE :: adiabatic(:)
    DOUBLE PRECISION, ALLOCATABLE :: mu(:)
    INTEGER, ALLOCATABLE :: pl(:)
    INTEGER, ALLOCATABLE :: ql(:)
    INTEGER, ALLOCATABLE :: pr(:)
    INTEGER, ALLOCATABLE :: qr(:)
    DOUBLE PRECISION, ALLOCATABLE :: iota(:)
    INTEGER, ALLOCATABLE :: lp(:)
    INTEGER, ALLOCATABLE :: lq(:)
    INTEGER, ALLOCATABLE :: rp(:)
    INTEGER, ALLOCATABLE :: rq(:)
    DOUBLE PRECISION, ALLOCATABLE :: oita(:)
    DOUBLE PRECISION :: mupftol
    INTEGER :: mupfits
    DOUBLE PRECISION, ALLOCATABLE :: Rac(:)
    DOUBLE PRECISION, ALLOCATABLE :: Zas(:)
    DOUBLE PRECISION, ALLOCATABLE :: Ras(:)
    DOUBLE PRECISION, ALLOCATABLE :: Zac(:)
    DOUBLE PRECISION, ALLOCATABLE :: Rbc(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: Zbs(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: Rbs(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: Zbc(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: Rwc(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: Zws(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: Rws(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: Zwc(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: Vns(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: Bns(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: Vnc(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: Bnc(:,:)
END TYPE physics
TYPE numerics
    INTEGER :: Linitialize
    INTEGER :: Lzerovac
    INTEGER :: Ndiscrete
    INTEGER :: Nquad
    INTEGER :: iMpol
    INTEGER :: iNtor
    INTEGER :: Lsparse
    INTEGER :: Lsvdiota
    INTEGER :: imethod
    INTEGER :: iorder
    INTEGER :: iprecon
    DOUBLE PRECISION :: iotatol
    INTEGER :: Lextrap
    INTEGER :: Mregular
END TYPE numerics
TYPE local
    INTEGER :: LBeltrami
    INTEGER :: Linitgues
    DOUBLE PRECISION :: maxrndgues
    INTEGER :: Lposdef
END TYPE local
TYPE global
    INTEGER :: Lfindzero
    DOUBLE PRECISION :: escale
    DOUBLE PRECISION :: opsilon
    DOUBLE PRECISION :: pcondense
    DOUBLE PRECISION :: epsilon
    DOUBLE PRECISION :: wpoloidal
    DOUBLE PRECISION :: upsilon
    DOUBLE PRECISION :: forcetol
    DOUBLE PRECISION :: c05xmax
    DOUBLE PRECISION :: c05xtol
    DOUBLE PRECISION :: c05factor
    LOGICAL :: LreadGF
    INTEGER :: mfreeits
    DOUBLE PRECISION :: bnstol
    DOUBLE PRECISION :: bnsblend
    DOUBLE PRECISION :: gBntol
    DOUBLE PRECISION :: gBnbld
    DOUBLE PRECISION :: vcasingeps
    DOUBLE PRECISION :: vcasingtol
    INTEGER :: vcasingits
    INTEGER :: vcasingper
    INTEGER :: mcasingcal
END TYPE global
TYPE diagnostics
    DOUBLE PRECISION :: odetol
    DOUBLE PRECISION :: absreq
    DOUBLE PRECISION :: relreq
    DOUBLE PRECISION :: absacc
    DOUBLE PRECISION :: epsr
    INTEGER :: nPpts
    INTEGER, ALLOCATABLE :: nPtrj(:)
    LOGICAL :: LHevalues
    LOGICAL :: LHevectors
    LOGICAL :: LHmatrix
    INTEGER :: Lperturbed
    INTEGER :: dpp
    INTEGER :: dqq
    INTEGER :: Lcheck
    LOGICAL :: Ltiming
    DOUBLE PRECISION :: fudge
    DOUBLE PRECISION :: scaling
END TYPE diagnostics
TYPE input
    TYPE(physics) :: physics
    TYPE(numerics) :: numerics
    TYPE(local) :: local
    TYPE(global) :: global
    TYPE(diagnostics) :: diagnostics
END TYPE input
TYPE output
    DOUBLE PRECISION, ALLOCATABLE :: Vns(:)
    DOUBLE PRECISION, ALLOCATABLE :: Bns(:)
    DOUBLE PRECISION, ALLOCATABLE :: Vnc(:)
    DOUBLE PRECISION, ALLOCATABLE :: Bnc(:)
    INTEGER :: mn
    INTEGER, ALLOCATABLE :: im(:)
    INTEGER, ALLOCATABLE :: in(:)
    INTEGER :: Mvol
    DOUBLE PRECISION, ALLOCATABLE :: Rbc(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: Zbs(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: Rbs(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: Zbc(:,:)
    DOUBLE PRECISION :: ForceErr
    DOUBLE PRECISION, ALLOCATABLE :: adiabatic(:)
    DOUBLE PRECISION, ALLOCATABLE :: helicity(:)
    DOUBLE PRECISION, ALLOCATABLE :: mu(:)
    DOUBLE PRECISION, ALLOCATABLE :: tflux(:)
    DOUBLE PRECISION, ALLOCATABLE :: pflux(:)
    DOUBLE PRECISION :: volume
    INTEGER :: Mrad
    DOUBLE PRECISION, ALLOCATABLE :: TT(:,:,:)
    DOUBLE PRECISION, ALLOCATABLE :: Btemn(:,:,:)
    DOUBLE PRECISION, ALLOCATABLE :: Bzemn(:,:,:)
    DOUBLE PRECISION, ALLOCATABLE :: Btomn(:,:,:)
    DOUBLE PRECISION, ALLOCATABLE :: Bzomn(:,:,:)
    DOUBLE PRECISION :: lmns
END TYPE output
TYPE vector_potential
    DOUBLE PRECISION, ALLOCATABLE :: Ate(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: Aze(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: Ato(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: Azo(:,:)
END TYPE vector_potential
TYPE SpecOutput
    DOUBLE PRECISION :: version
    TYPE(input) :: input
    TYPE(output) :: output
    TYPE(vector_potential) :: vector_potential
END TYPE SpecOutput
contains
subroutine loadSpec(s, filename, ierr)
  use hdf5
  implicit none
  type(SpecOutput), intent(inout) :: s                 ! target datastructure
  character(len=*), intent(in)    :: filename          ! filename to load
  integer, intent(out), optional  :: ierr              ! error flag; .eq.0 if ok
  integer                         :: hdfier            ! error flag for HDF5 API calls
  integer(hid_t)                  :: file_id           ! identifier for current file
  integer(hid_t)                  :: dset_id           ! temporary dataset id
  integer(hid_t)                  :: dataspace         ! dataspace used to query dataset size
  integer(hsize_t)                :: dims_1(1) ! current dimensions of rank-1 dataset
  integer(hsize_t)                :: dims_2(2) ! current dimensions of rank-2 dataset
  integer(hsize_t)                :: dims_3(3) ! current dimensions of rank-3 dataset
  integer(hsize_t)                :: max_dims_1(1)     ! maximum dimensions of rank-1 dataset
  integer(hsize_t)                :: max_dims_2(2)     ! maximum dimensions of rank-2 dataset
  integer(hsize_t)                :: max_dims_3(3)     ! maximum dimensions of rank-3 dataset
  integer                         :: logical_tmp       ! temporary integer used to read logicals
  
  call h5open_f(hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening HDF5 library" ; goto 9999 ; endif

  call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening HDF5 file '",filename,"'" ; goto 9998 ; endif

! /input/physics/Igeometry --> s%input%physics%Igeometry; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/physics/Igeometry", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/Igeometry'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%physics%Igeometry, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/Igeometry'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/Igeometry'" ; goto 9998 ; endif

! /input/physics/Istellsym --> s%input%physics%Istellsym; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/physics/Istellsym", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/Istellsym'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%physics%Istellsym, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/Istellsym'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/Istellsym'" ; goto 9998 ; endif

! /input/physics/Lfreebound --> s%input%physics%Lfreebound; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/physics/Lfreebound", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/Lfreebound'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%physics%Lfreebound, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/Lfreebound'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/Lfreebound'" ; goto 9998 ; endif

! /input/physics/phiedge --> s%input%physics%phiedge; rank=0; h5type=H5T_NATIVE_DOUBLE
  call h5dopen_f(file_id, "/input/physics/phiedge", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/phiedge'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%physics%phiedge, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/phiedge'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/phiedge'" ; goto 9998 ; endif

! /input/physics/curtor --> s%input%physics%curtor; rank=0; h5type=H5T_NATIVE_DOUBLE
  call h5dopen_f(file_id, "/input/physics/curtor", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/curtor'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%physics%curtor, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/curtor'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/curtor'" ; goto 9998 ; endif

! /input/physics/curpol --> s%input%physics%curpol; rank=0; h5type=H5T_NATIVE_DOUBLE
  call h5dopen_f(file_id, "/input/physics/curpol", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/curpol'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%physics%curpol, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/curpol'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/curpol'" ; goto 9998 ; endif

! /input/physics/gamma --> s%input%physics%gamma; rank=0; h5type=H5T_NATIVE_DOUBLE
  call h5dopen_f(file_id, "/input/physics/gamma", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/gamma'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%physics%gamma, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/gamma'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/gamma'" ; goto 9998 ; endif

! /input/physics/Nfp --> s%input%physics%Nfp; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/physics/Nfp", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/Nfp'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%physics%Nfp, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/Nfp'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/Nfp'" ; goto 9998 ; endif

! /input/physics/Nvol --> s%input%physics%Nvol; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/physics/Nvol", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/Nvol'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%physics%Nvol, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/Nvol'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/Nvol'" ; goto 9998 ; endif

! /input/physics/Mpol --> s%input%physics%Mpol; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/physics/Mpol", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/Mpol'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%physics%Mpol, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/Mpol'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/Mpol'" ; goto 9998 ; endif

! /input/physics/Ntor --> s%input%physics%Ntor; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/physics/Ntor", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/Ntor'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%physics%Ntor, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/Ntor'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/Ntor'" ; goto 9998 ; endif

! /input/physics/Lrad --> s%input%physics%Lrad(1:dims_1(1)); rank=1
  call h5dopen_f(file_id, "/input/physics/Lrad", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/Lrad'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/input/physics/Lrad'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_1, max_dims_1, hdfier)
  if (hdfier.ne.1) then ; write(*,*) "unexpected rank of dataset '/input/physics/Lrad': ",hdfier," .ne. 1" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/input/physics/Lrad'" ; goto 9998 ; endif
  
  allocate(s%input%physics%Lrad(1:dims_1(1)))
  
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%physics%Lrad(1:dims_1(1)), dims_1, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/Lrad'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/Lrad'" ; goto 9998 ; endif

! /input/physics/Lconstraint --> s%input%physics%Lconstraint; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/physics/Lconstraint", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/Lconstraint'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%physics%Lconstraint, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/Lconstraint'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/Lconstraint'" ; goto 9998 ; endif

! /input/physics/tflux --> s%input%physics%tflux(1:dims_1(1)); rank=1
  call h5dopen_f(file_id, "/input/physics/tflux", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/tflux'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/input/physics/tflux'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_1, max_dims_1, hdfier)
  if (hdfier.ne.1) then ; write(*,*) "unexpected rank of dataset '/input/physics/tflux': ",hdfier," .ne. 1" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/input/physics/tflux'" ; goto 9998 ; endif
  
  allocate(s%input%physics%tflux(1:dims_1(1)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%physics%tflux(1:dims_1(1)), dims_1, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/tflux'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/tflux'" ; goto 9998 ; endif

! /input/physics/pflux --> s%input%physics%pflux(1:dims_1(1)); rank=1
  call h5dopen_f(file_id, "/input/physics/pflux", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/pflux'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/input/physics/pflux'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_1, max_dims_1, hdfier)
  if (hdfier.ne.1) then ; write(*,*) "unexpected rank of dataset '/input/physics/pflux': ",hdfier," .ne. 1" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/input/physics/pflux'" ; goto 9998 ; endif
  
  allocate(s%input%physics%pflux(1:dims_1(1)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%physics%pflux(1:dims_1(1)), dims_1, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/pflux'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/pflux'" ; goto 9998 ; endif

! /input/physics/helicity --> s%input%physics%helicity(1:dims_1(1)); rank=1
  call h5dopen_f(file_id, "/input/physics/helicity", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/helicity'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/input/physics/helicity'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_1, max_dims_1, hdfier)
  if (hdfier.ne.1) then ; write(*,*) "unexpected rank of dataset '/input/physics/helicity': ",hdfier," .ne. 1" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/input/physics/helicity'" ; goto 9998 ; endif
  
  allocate(s%input%physics%helicity(1:dims_1(1)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%physics%helicity(1:dims_1(1)), dims_1, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/helicity'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/helicity'" ; goto 9998 ; endif

! /input/physics/pscale --> s%input%physics%pscale; rank=0; h5type=H5T_NATIVE_DOUBLE
  call h5dopen_f(file_id, "/input/physics/pscale", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/pscale'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%physics%pscale, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/pscale'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/pscale'" ; goto 9998 ; endif

! /input/physics/pressure --> s%input%physics%pressure(1:dims_1(1)); rank=1
  call h5dopen_f(file_id, "/input/physics/pressure", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/pressure'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/input/physics/pressure'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_1, max_dims_1, hdfier)
  if (hdfier.ne.1) then ; write(*,*) "unexpected rank of dataset '/input/physics/pressure': ",hdfier," .ne. 1" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/input/physics/pressure'" ; goto 9998 ; endif
  
  allocate(s%input%physics%pressure(1:dims_1(1)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%physics%pressure(1:dims_1(1)), dims_1, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/pressure'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/pressure'" ; goto 9998 ; endif

! /input/physics/Ladiabatic --> s%input%physics%Ladiabatic; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/physics/Ladiabatic", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/Ladiabatic'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%physics%Ladiabatic, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/Ladiabatic'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/Ladiabatic'" ; goto 9998 ; endif

! /input/physics/adiabatic --> s%input%physics%adiabatic(1:dims_1(1)); rank=1
  call h5dopen_f(file_id, "/input/physics/adiabatic", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/adiabatic'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/input/physics/adiabatic'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_1, max_dims_1, hdfier)
  if (hdfier.ne.1) then ; write(*,*) "unexpected rank of dataset '/input/physics/adiabatic': ",hdfier," .ne. 1" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/input/physics/adiabatic'" ; goto 9998 ; endif
  
  allocate(s%input%physics%adiabatic(1:dims_1(1)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%physics%adiabatic(1:dims_1(1)), dims_1, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/adiabatic'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/adiabatic'" ; goto 9998 ; endif

! /input/physics/mu --> s%input%physics%mu(1:dims_1(1)); rank=1
  call h5dopen_f(file_id, "/input/physics/mu", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/mu'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/input/physics/mu'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_1, max_dims_1, hdfier)
  if (hdfier.ne.1) then ; write(*,*) "unexpected rank of dataset '/input/physics/mu': ",hdfier," .ne. 1" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/input/physics/mu'" ; goto 9998 ; endif
  
  allocate(s%input%physics%mu(1:dims_1(1)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%physics%mu(1:dims_1(1)), dims_1, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/mu'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/mu'" ; goto 9998 ; endif

! /input/physics/pl --> s%input%physics%pl(1:dims_1(1)); rank=1
  call h5dopen_f(file_id, "/input/physics/pl", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/pl'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/input/physics/pl'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_1, max_dims_1, hdfier)
  if (hdfier.ne.1) then ; write(*,*) "unexpected rank of dataset '/input/physics/pl': ",hdfier," .ne. 1" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/input/physics/pl'" ; goto 9998 ; endif
  
  allocate(s%input%physics%pl(1:dims_1(1)))
  
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%physics%pl(1:dims_1(1)), dims_1, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/pl'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/pl'" ; goto 9998 ; endif

! /input/physics/ql --> s%input%physics%ql(1:dims_1(1)); rank=1
  call h5dopen_f(file_id, "/input/physics/ql", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/ql'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/input/physics/ql'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_1, max_dims_1, hdfier)
  if (hdfier.ne.1) then ; write(*,*) "unexpected rank of dataset '/input/physics/ql': ",hdfier," .ne. 1" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/input/physics/ql'" ; goto 9998 ; endif
  
  allocate(s%input%physics%ql(1:dims_1(1)))
  
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%physics%ql(1:dims_1(1)), dims_1, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/ql'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/ql'" ; goto 9998 ; endif

! /input/physics/pr --> s%input%physics%pr(1:dims_1(1)); rank=1
  call h5dopen_f(file_id, "/input/physics/pr", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/pr'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/input/physics/pr'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_1, max_dims_1, hdfier)
  if (hdfier.ne.1) then ; write(*,*) "unexpected rank of dataset '/input/physics/pr': ",hdfier," .ne. 1" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/input/physics/pr'" ; goto 9998 ; endif
  
  allocate(s%input%physics%pr(1:dims_1(1)))
  
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%physics%pr(1:dims_1(1)), dims_1, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/pr'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/pr'" ; goto 9998 ; endif

! /input/physics/qr --> s%input%physics%qr(1:dims_1(1)); rank=1
  call h5dopen_f(file_id, "/input/physics/qr", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/qr'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/input/physics/qr'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_1, max_dims_1, hdfier)
  if (hdfier.ne.1) then ; write(*,*) "unexpected rank of dataset '/input/physics/qr': ",hdfier," .ne. 1" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/input/physics/qr'" ; goto 9998 ; endif
  
  allocate(s%input%physics%qr(1:dims_1(1)))
  
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%physics%qr(1:dims_1(1)), dims_1, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/qr'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/qr'" ; goto 9998 ; endif

! /input/physics/iota --> s%input%physics%iota(1:dims_1(1)); rank=1
  call h5dopen_f(file_id, "/input/physics/iota", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/iota'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/input/physics/iota'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_1, max_dims_1, hdfier)
  if (hdfier.ne.1) then ; write(*,*) "unexpected rank of dataset '/input/physics/iota': ",hdfier," .ne. 1" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/input/physics/iota'" ; goto 9998 ; endif
  
  allocate(s%input%physics%iota(1:dims_1(1)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%physics%iota(1:dims_1(1)), dims_1, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/iota'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/iota'" ; goto 9998 ; endif

! /input/physics/lp --> s%input%physics%lp(1:dims_1(1)); rank=1
  call h5dopen_f(file_id, "/input/physics/lp", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/lp'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/input/physics/lp'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_1, max_dims_1, hdfier)
  if (hdfier.ne.1) then ; write(*,*) "unexpected rank of dataset '/input/physics/lp': ",hdfier," .ne. 1" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/input/physics/lp'" ; goto 9998 ; endif
  
  allocate(s%input%physics%lp(1:dims_1(1)))
  
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%physics%lp(1:dims_1(1)), dims_1, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/lp'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/lp'" ; goto 9998 ; endif

! /input/physics/lq --> s%input%physics%lq(1:dims_1(1)); rank=1
  call h5dopen_f(file_id, "/input/physics/lq", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/lq'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/input/physics/lq'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_1, max_dims_1, hdfier)
  if (hdfier.ne.1) then ; write(*,*) "unexpected rank of dataset '/input/physics/lq': ",hdfier," .ne. 1" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/input/physics/lq'" ; goto 9998 ; endif
  
  allocate(s%input%physics%lq(1:dims_1(1)))
  
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%physics%lq(1:dims_1(1)), dims_1, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/lq'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/lq'" ; goto 9998 ; endif

! /input/physics/rp --> s%input%physics%rp(1:dims_1(1)); rank=1
  call h5dopen_f(file_id, "/input/physics/rp", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/rp'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/input/physics/rp'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_1, max_dims_1, hdfier)
  if (hdfier.ne.1) then ; write(*,*) "unexpected rank of dataset '/input/physics/rp': ",hdfier," .ne. 1" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/input/physics/rp'" ; goto 9998 ; endif
  
  allocate(s%input%physics%rp(1:dims_1(1)))
  
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%physics%rp(1:dims_1(1)), dims_1, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/rp'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/rp'" ; goto 9998 ; endif

! /input/physics/rq --> s%input%physics%rq(1:dims_1(1)); rank=1
  call h5dopen_f(file_id, "/input/physics/rq", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/rq'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/input/physics/rq'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_1, max_dims_1, hdfier)
  if (hdfier.ne.1) then ; write(*,*) "unexpected rank of dataset '/input/physics/rq': ",hdfier," .ne. 1" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/input/physics/rq'" ; goto 9998 ; endif
  
  allocate(s%input%physics%rq(1:dims_1(1)))
  
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%physics%rq(1:dims_1(1)), dims_1, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/rq'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/rq'" ; goto 9998 ; endif

! /input/physics/oita --> s%input%physics%oita(1:dims_1(1)); rank=1
  call h5dopen_f(file_id, "/input/physics/oita", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/oita'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/input/physics/oita'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_1, max_dims_1, hdfier)
  if (hdfier.ne.1) then ; write(*,*) "unexpected rank of dataset '/input/physics/oita': ",hdfier," .ne. 1" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/input/physics/oita'" ; goto 9998 ; endif
  
  allocate(s%input%physics%oita(1:dims_1(1)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%physics%oita(1:dims_1(1)), dims_1, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/oita'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/oita'" ; goto 9998 ; endif

! /input/physics/mupftol --> s%input%physics%mupftol; rank=0; h5type=H5T_NATIVE_DOUBLE
  call h5dopen_f(file_id, "/input/physics/mupftol", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/mupftol'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%physics%mupftol, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/mupftol'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/mupftol'" ; goto 9998 ; endif

! /input/physics/mupfits --> s%input%physics%mupfits; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/physics/mupfits", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/mupfits'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%physics%mupfits, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/mupfits'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/mupfits'" ; goto 9998 ; endif

! /input/physics/Rac --> s%input%physics%Rac(1:dims_1(1)); rank=1
  call h5dopen_f(file_id, "/input/physics/Rac", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/Rac'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/input/physics/Rac'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_1, max_dims_1, hdfier)
  if (hdfier.ne.1) then ; write(*,*) "unexpected rank of dataset '/input/physics/Rac': ",hdfier," .ne. 1" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/input/physics/Rac'" ; goto 9998 ; endif
  
  allocate(s%input%physics%Rac(1:dims_1(1)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%physics%Rac(1:dims_1(1)), dims_1, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/Rac'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/Rac'" ; goto 9998 ; endif

! /input/physics/Zas --> s%input%physics%Zas(1:dims_1(1)); rank=1
  call h5dopen_f(file_id, "/input/physics/Zas", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/Zas'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/input/physics/Zas'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_1, max_dims_1, hdfier)
  if (hdfier.ne.1) then ; write(*,*) "unexpected rank of dataset '/input/physics/Zas': ",hdfier," .ne. 1" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/input/physics/Zas'" ; goto 9998 ; endif
  
  allocate(s%input%physics%Zas(1:dims_1(1)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%physics%Zas(1:dims_1(1)), dims_1, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/Zas'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/Zas'" ; goto 9998 ; endif

! /input/physics/Ras --> s%input%physics%Ras(1:dims_1(1)); rank=1
  call h5dopen_f(file_id, "/input/physics/Ras", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/Ras'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/input/physics/Ras'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_1, max_dims_1, hdfier)
  if (hdfier.ne.1) then ; write(*,*) "unexpected rank of dataset '/input/physics/Ras': ",hdfier," .ne. 1" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/input/physics/Ras'" ; goto 9998 ; endif
  
  allocate(s%input%physics%Ras(1:dims_1(1)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%physics%Ras(1:dims_1(1)), dims_1, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/Ras'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/Ras'" ; goto 9998 ; endif

! /input/physics/Zac --> s%input%physics%Zac(1:dims_1(1)); rank=1
  call h5dopen_f(file_id, "/input/physics/Zac", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/Zac'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/input/physics/Zac'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_1, max_dims_1, hdfier)
  if (hdfier.ne.1) then ; write(*,*) "unexpected rank of dataset '/input/physics/Zac': ",hdfier," .ne. 1" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/input/physics/Zac'" ; goto 9998 ; endif
  
  allocate(s%input%physics%Zac(1:dims_1(1)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%physics%Zac(1:dims_1(1)), dims_1, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/Zac'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/Zac'" ; goto 9998 ; endif

! /input/physics/Rbc --> s%input%physics%Rbc(1:dims_2(1), 1:dims_2(2)); rank=2
  call h5dopen_f(file_id, "/input/physics/Rbc", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/Rbc'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/input/physics/Rbc'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_2, max_dims_2, hdfier)
  if (hdfier.ne.2) then ; write(*,*) "unexpected rank of dataset '/input/physics/Rbc': ",hdfier," .ne. 2" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/input/physics/Rbc'" ; goto 9998 ; endif
  
  allocate(s%input%physics%Rbc(1:dims_2(1), 1:dims_2(2)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%physics%Rbc(1:dims_2(1), 1:dims_2(2)), dims_2, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/Rbc'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/Rbc'" ; goto 9998 ; endif

! /input/physics/Zbs --> s%input%physics%Zbs(1:dims_2(1), 1:dims_2(2)); rank=2
  call h5dopen_f(file_id, "/input/physics/Zbs", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/Zbs'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/input/physics/Zbs'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_2, max_dims_2, hdfier)
  if (hdfier.ne.2) then ; write(*,*) "unexpected rank of dataset '/input/physics/Zbs': ",hdfier," .ne. 2" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/input/physics/Zbs'" ; goto 9998 ; endif
  
  allocate(s%input%physics%Zbs(1:dims_2(1), 1:dims_2(2)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%physics%Zbs(1:dims_2(1), 1:dims_2(2)), dims_2, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/Zbs'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/Zbs'" ; goto 9998 ; endif

! /input/physics/Rbs --> s%input%physics%Rbs(1:dims_2(1), 1:dims_2(2)); rank=2
  call h5dopen_f(file_id, "/input/physics/Rbs", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/Rbs'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/input/physics/Rbs'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_2, max_dims_2, hdfier)
  if (hdfier.ne.2) then ; write(*,*) "unexpected rank of dataset '/input/physics/Rbs': ",hdfier," .ne. 2" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/input/physics/Rbs'" ; goto 9998 ; endif
  
  allocate(s%input%physics%Rbs(1:dims_2(1), 1:dims_2(2)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%physics%Rbs(1:dims_2(1), 1:dims_2(2)), dims_2, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/Rbs'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/Rbs'" ; goto 9998 ; endif

! /input/physics/Zbc --> s%input%physics%Zbc(1:dims_2(1), 1:dims_2(2)); rank=2
  call h5dopen_f(file_id, "/input/physics/Zbc", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/Zbc'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/input/physics/Zbc'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_2, max_dims_2, hdfier)
  if (hdfier.ne.2) then ; write(*,*) "unexpected rank of dataset '/input/physics/Zbc': ",hdfier," .ne. 2" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/input/physics/Zbc'" ; goto 9998 ; endif
  
  allocate(s%input%physics%Zbc(1:dims_2(1), 1:dims_2(2)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%physics%Zbc(1:dims_2(1), 1:dims_2(2)), dims_2, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/Zbc'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/Zbc'" ; goto 9998 ; endif

! /input/physics/Rwc --> s%input%physics%Rwc(1:dims_2(1), 1:dims_2(2)); rank=2
  call h5dopen_f(file_id, "/input/physics/Rwc", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/Rwc'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/input/physics/Rwc'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_2, max_dims_2, hdfier)
  if (hdfier.ne.2) then ; write(*,*) "unexpected rank of dataset '/input/physics/Rwc': ",hdfier," .ne. 2" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/input/physics/Rwc'" ; goto 9998 ; endif
  
  allocate(s%input%physics%Rwc(1:dims_2(1), 1:dims_2(2)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%physics%Rwc(1:dims_2(1), 1:dims_2(2)), dims_2, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/Rwc'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/Rwc'" ; goto 9998 ; endif

! /input/physics/Zws --> s%input%physics%Zws(1:dims_2(1), 1:dims_2(2)); rank=2
  call h5dopen_f(file_id, "/input/physics/Zws", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/Zws'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/input/physics/Zws'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_2, max_dims_2, hdfier)
  if (hdfier.ne.2) then ; write(*,*) "unexpected rank of dataset '/input/physics/Zws': ",hdfier," .ne. 2" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/input/physics/Zws'" ; goto 9998 ; endif
  
  allocate(s%input%physics%Zws(1:dims_2(1), 1:dims_2(2)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%physics%Zws(1:dims_2(1), 1:dims_2(2)), dims_2, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/Zws'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/Zws'" ; goto 9998 ; endif

! /input/physics/Rws --> s%input%physics%Rws(1:dims_2(1), 1:dims_2(2)); rank=2
  call h5dopen_f(file_id, "/input/physics/Rws", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/Rws'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/input/physics/Rws'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_2, max_dims_2, hdfier)
  if (hdfier.ne.2) then ; write(*,*) "unexpected rank of dataset '/input/physics/Rws': ",hdfier," .ne. 2" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/input/physics/Rws'" ; goto 9998 ; endif
  
  allocate(s%input%physics%Rws(1:dims_2(1), 1:dims_2(2)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%physics%Rws(1:dims_2(1), 1:dims_2(2)), dims_2, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/Rws'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/Rws'" ; goto 9998 ; endif

! /input/physics/Zwc --> s%input%physics%Zwc(1:dims_2(1), 1:dims_2(2)); rank=2
  call h5dopen_f(file_id, "/input/physics/Zwc", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/Zwc'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/input/physics/Zwc'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_2, max_dims_2, hdfier)
  if (hdfier.ne.2) then ; write(*,*) "unexpected rank of dataset '/input/physics/Zwc': ",hdfier," .ne. 2" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/input/physics/Zwc'" ; goto 9998 ; endif
  
  allocate(s%input%physics%Zwc(1:dims_2(1), 1:dims_2(2)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%physics%Zwc(1:dims_2(1), 1:dims_2(2)), dims_2, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/Zwc'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/Zwc'" ; goto 9998 ; endif

! /input/physics/Vns --> s%input%physics%Vns(1:dims_2(1), 1:dims_2(2)); rank=2
  call h5dopen_f(file_id, "/input/physics/Vns", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/Vns'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/input/physics/Vns'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_2, max_dims_2, hdfier)
  if (hdfier.ne.2) then ; write(*,*) "unexpected rank of dataset '/input/physics/Vns': ",hdfier," .ne. 2" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/input/physics/Vns'" ; goto 9998 ; endif
  
  allocate(s%input%physics%Vns(1:dims_2(1), 1:dims_2(2)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%physics%Vns(1:dims_2(1), 1:dims_2(2)), dims_2, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/Vns'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/Vns'" ; goto 9998 ; endif

! /input/physics/Bns --> s%input%physics%Bns(1:dims_2(1), 1:dims_2(2)); rank=2
  call h5dopen_f(file_id, "/input/physics/Bns", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/Bns'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/input/physics/Bns'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_2, max_dims_2, hdfier)
  if (hdfier.ne.2) then ; write(*,*) "unexpected rank of dataset '/input/physics/Bns': ",hdfier," .ne. 2" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/input/physics/Bns'" ; goto 9998 ; endif
  
  allocate(s%input%physics%Bns(1:dims_2(1), 1:dims_2(2)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%physics%Bns(1:dims_2(1), 1:dims_2(2)), dims_2, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/Bns'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/Bns'" ; goto 9998 ; endif

! /input/physics/Vnc --> s%input%physics%Vnc(1:dims_2(1), 1:dims_2(2)); rank=2
  call h5dopen_f(file_id, "/input/physics/Vnc", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/Vnc'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/input/physics/Vnc'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_2, max_dims_2, hdfier)
  if (hdfier.ne.2) then ; write(*,*) "unexpected rank of dataset '/input/physics/Vnc': ",hdfier," .ne. 2" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/input/physics/Vnc'" ; goto 9998 ; endif
  
  allocate(s%input%physics%Vnc(1:dims_2(1), 1:dims_2(2)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%physics%Vnc(1:dims_2(1), 1:dims_2(2)), dims_2, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/Vnc'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/Vnc'" ; goto 9998 ; endif

! /input/physics/Bnc --> s%input%physics%Bnc(1:dims_2(1), 1:dims_2(2)); rank=2
  call h5dopen_f(file_id, "/input/physics/Bnc", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/physics/Bnc'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/input/physics/Bnc'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_2, max_dims_2, hdfier)
  if (hdfier.ne.2) then ; write(*,*) "unexpected rank of dataset '/input/physics/Bnc': ",hdfier," .ne. 2" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/input/physics/Bnc'" ; goto 9998 ; endif
  
  allocate(s%input%physics%Bnc(1:dims_2(1), 1:dims_2(2)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%physics%Bnc(1:dims_2(1), 1:dims_2(2)), dims_2, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/physics/Bnc'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/physics/Bnc'" ; goto 9998 ; endif

! /input/numerics/Linitialize --> s%input%numerics%Linitialize; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/numerics/Linitialize", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/numerics/Linitialize'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%numerics%Linitialize, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/numerics/Linitialize'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/numerics/Linitialize'" ; goto 9998 ; endif

! /input/numerics/Lzerovac --> s%input%numerics%Lzerovac; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/numerics/Lzerovac", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/numerics/Lzerovac'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%numerics%Lzerovac, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/numerics/Lzerovac'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/numerics/Lzerovac'" ; goto 9998 ; endif

! /input/numerics/Ndiscrete --> s%input%numerics%Ndiscrete; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/numerics/Ndiscrete", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/numerics/Ndiscrete'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%numerics%Ndiscrete, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/numerics/Ndiscrete'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/numerics/Ndiscrete'" ; goto 9998 ; endif

! /input/numerics/Nquad --> s%input%numerics%Nquad; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/numerics/Nquad", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/numerics/Nquad'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%numerics%Nquad, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/numerics/Nquad'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/numerics/Nquad'" ; goto 9998 ; endif

! /input/numerics/iMpol --> s%input%numerics%iMpol; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/numerics/iMpol", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/numerics/iMpol'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%numerics%iMpol, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/numerics/iMpol'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/numerics/iMpol'" ; goto 9998 ; endif

! /input/numerics/iNtor --> s%input%numerics%iNtor; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/numerics/iNtor", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/numerics/iNtor'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%numerics%iNtor, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/numerics/iNtor'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/numerics/iNtor'" ; goto 9998 ; endif

! /input/numerics/Lsparse --> s%input%numerics%Lsparse; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/numerics/Lsparse", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/numerics/Lsparse'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%numerics%Lsparse, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/numerics/Lsparse'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/numerics/Lsparse'" ; goto 9998 ; endif

! /input/numerics/Lsvdiota --> s%input%numerics%Lsvdiota; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/numerics/Lsvdiota", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/numerics/Lsvdiota'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%numerics%Lsvdiota, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/numerics/Lsvdiota'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/numerics/Lsvdiota'" ; goto 9998 ; endif

! /input/numerics/imethod --> s%input%numerics%imethod; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/numerics/imethod", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/numerics/imethod'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%numerics%imethod, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/numerics/imethod'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/numerics/imethod'" ; goto 9998 ; endif

! /input/numerics/iorder --> s%input%numerics%iorder; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/numerics/iorder", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/numerics/iorder'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%numerics%iorder, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/numerics/iorder'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/numerics/iorder'" ; goto 9998 ; endif

! /input/numerics/iprecon --> s%input%numerics%iprecon; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/numerics/iprecon", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/numerics/iprecon'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%numerics%iprecon, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/numerics/iprecon'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/numerics/iprecon'" ; goto 9998 ; endif

! /input/numerics/iotatol --> s%input%numerics%iotatol; rank=0; h5type=H5T_NATIVE_DOUBLE
  call h5dopen_f(file_id, "/input/numerics/iotatol", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/numerics/iotatol'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%numerics%iotatol, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/numerics/iotatol'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/numerics/iotatol'" ; goto 9998 ; endif

! /input/numerics/Lextrap --> s%input%numerics%Lextrap; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/numerics/Lextrap", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/numerics/Lextrap'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%numerics%Lextrap, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/numerics/Lextrap'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/numerics/Lextrap'" ; goto 9998 ; endif

! /input/numerics/Mregular --> s%input%numerics%Mregular; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/numerics/Mregular", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/numerics/Mregular'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%numerics%Mregular, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/numerics/Mregular'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/numerics/Mregular'" ; goto 9998 ; endif

! /input/local/LBeltrami --> s%input%local%LBeltrami; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/local/LBeltrami", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/local/LBeltrami'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%local%LBeltrami, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/local/LBeltrami'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/local/LBeltrami'" ; goto 9998 ; endif

! /input/local/Linitgues --> s%input%local%Linitgues; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/local/Linitgues", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/local/Linitgues'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%local%Linitgues, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/local/Linitgues'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/local/Linitgues'" ; goto 9998 ; endif

! /input/local/maxrndgues --> s%input%local%maxrndgues; rank=0; h5type=H5T_NATIVE_DOUBLE
  call h5dopen_f(file_id, "/input/local/maxrndgues", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/local/maxrndgues'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%local%maxrndgues, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/local/maxrndgues'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/local/maxrndgues'" ; goto 9998 ; endif

! /input/local/Lposdef --> s%input%local%Lposdef; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/local/Lposdef", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/local/Lposdef'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%local%Lposdef, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/local/Lposdef'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/local/Lposdef'" ; goto 9998 ; endif

! /input/global/Lfindzero --> s%input%global%Lfindzero; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/global/Lfindzero", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/global/Lfindzero'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%global%Lfindzero, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/global/Lfindzero'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/global/Lfindzero'" ; goto 9998 ; endif

! /input/global/escale --> s%input%global%escale; rank=0; h5type=H5T_NATIVE_DOUBLE
  call h5dopen_f(file_id, "/input/global/escale", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/global/escale'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%global%escale, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/global/escale'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/global/escale'" ; goto 9998 ; endif

! /input/global/opsilon --> s%input%global%opsilon; rank=0; h5type=H5T_NATIVE_DOUBLE
  call h5dopen_f(file_id, "/input/global/opsilon", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/global/opsilon'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%global%opsilon, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/global/opsilon'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/global/opsilon'" ; goto 9998 ; endif

! /input/global/pcondense --> s%input%global%pcondense; rank=0; h5type=H5T_NATIVE_DOUBLE
  call h5dopen_f(file_id, "/input/global/pcondense", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/global/pcondense'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%global%pcondense, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/global/pcondense'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/global/pcondense'" ; goto 9998 ; endif

! /input/global/epsilon --> s%input%global%epsilon; rank=0; h5type=H5T_NATIVE_DOUBLE
  call h5dopen_f(file_id, "/input/global/epsilon", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/global/epsilon'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%global%epsilon, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/global/epsilon'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/global/epsilon'" ; goto 9998 ; endif

! /input/global/wpoloidal --> s%input%global%wpoloidal; rank=0; h5type=H5T_NATIVE_DOUBLE
  call h5dopen_f(file_id, "/input/global/wpoloidal", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/global/wpoloidal'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%global%wpoloidal, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/global/wpoloidal'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/global/wpoloidal'" ; goto 9998 ; endif

! /input/global/upsilon --> s%input%global%upsilon; rank=0; h5type=H5T_NATIVE_DOUBLE
  call h5dopen_f(file_id, "/input/global/upsilon", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/global/upsilon'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%global%upsilon, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/global/upsilon'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/global/upsilon'" ; goto 9998 ; endif

! /input/global/forcetol --> s%input%global%forcetol; rank=0; h5type=H5T_NATIVE_DOUBLE
  call h5dopen_f(file_id, "/input/global/forcetol", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/global/forcetol'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%global%forcetol, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/global/forcetol'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/global/forcetol'" ; goto 9998 ; endif

! /input/global/c05xmax --> s%input%global%c05xmax; rank=0; h5type=H5T_NATIVE_DOUBLE
  call h5dopen_f(file_id, "/input/global/c05xmax", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/global/c05xmax'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%global%c05xmax, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/global/c05xmax'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/global/c05xmax'" ; goto 9998 ; endif

! /input/global/c05xtol --> s%input%global%c05xtol; rank=0; h5type=H5T_NATIVE_DOUBLE
  call h5dopen_f(file_id, "/input/global/c05xtol", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/global/c05xtol'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%global%c05xtol, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/global/c05xtol'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/global/c05xtol'" ; goto 9998 ; endif

! /input/global/c05factor --> s%input%global%c05factor; rank=0; h5type=H5T_NATIVE_DOUBLE
  call h5dopen_f(file_id, "/input/global/c05factor", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/global/c05factor'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%global%c05factor, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/global/c05factor'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/global/c05factor'" ; goto 9998 ; endif

! /input/global/LreadGF --> s%input%global%LreadGF; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/global/LreadGF", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/global/LreadGF'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, logical_tmp, int((/1/), HSIZE_T), hdfier)
  s%input%global%LreadGF = merge(.TRUE., .FALSE., logical_tmp.ne.0)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/global/LreadGF'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/global/LreadGF'" ; goto 9998 ; endif

! /input/global/mfreeits --> s%input%global%mfreeits; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/global/mfreeits", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/global/mfreeits'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%global%mfreeits, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/global/mfreeits'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/global/mfreeits'" ; goto 9998 ; endif

! /input/global/bnstol --> s%input%global%bnstol; rank=0; h5type=H5T_NATIVE_DOUBLE
  call h5dopen_f(file_id, "/input/global/bnstol", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/global/bnstol'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%global%bnstol, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/global/bnstol'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/global/bnstol'" ; goto 9998 ; endif

! /input/global/bnsblend --> s%input%global%bnsblend; rank=0; h5type=H5T_NATIVE_DOUBLE
  call h5dopen_f(file_id, "/input/global/bnsblend", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/global/bnsblend'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%global%bnsblend, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/global/bnsblend'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/global/bnsblend'" ; goto 9998 ; endif

! /input/global/gBntol --> s%input%global%gBntol; rank=0; h5type=H5T_NATIVE_DOUBLE
  call h5dopen_f(file_id, "/input/global/gBntol", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/global/gBntol'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%global%gBntol, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/global/gBntol'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/global/gBntol'" ; goto 9998 ; endif

! /input/global/gBnbld --> s%input%global%gBnbld; rank=0; h5type=H5T_NATIVE_DOUBLE
  call h5dopen_f(file_id, "/input/global/gBnbld", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/global/gBnbld'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%global%gBnbld, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/global/gBnbld'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/global/gBnbld'" ; goto 9998 ; endif

! /input/global/vcasingeps --> s%input%global%vcasingeps; rank=0; h5type=H5T_NATIVE_DOUBLE
  call h5dopen_f(file_id, "/input/global/vcasingeps", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/global/vcasingeps'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%global%vcasingeps, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/global/vcasingeps'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/global/vcasingeps'" ; goto 9998 ; endif

! /input/global/vcasingtol --> s%input%global%vcasingtol; rank=0; h5type=H5T_NATIVE_DOUBLE
  call h5dopen_f(file_id, "/input/global/vcasingtol", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/global/vcasingtol'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%global%vcasingtol, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/global/vcasingtol'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/global/vcasingtol'" ; goto 9998 ; endif

! /input/global/vcasingits --> s%input%global%vcasingits; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/global/vcasingits", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/global/vcasingits'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%global%vcasingits, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/global/vcasingits'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/global/vcasingits'" ; goto 9998 ; endif

! /input/global/vcasingper --> s%input%global%vcasingper; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/global/vcasingper", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/global/vcasingper'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%global%vcasingper, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/global/vcasingper'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/global/vcasingper'" ; goto 9998 ; endif

! /input/global/mcasingcal --> s%input%global%mcasingcal; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/global/mcasingcal", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/global/mcasingcal'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%global%mcasingcal, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/global/mcasingcal'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/global/mcasingcal'" ; goto 9998 ; endif

! /input/diagnostics/odetol --> s%input%diagnostics%odetol; rank=0; h5type=H5T_NATIVE_DOUBLE
  call h5dopen_f(file_id, "/input/diagnostics/odetol", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/diagnostics/odetol'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%diagnostics%odetol, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/diagnostics/odetol'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/diagnostics/odetol'" ; goto 9998 ; endif

! /input/diagnostics/absreq --> s%input%diagnostics%absreq; rank=0; h5type=H5T_NATIVE_DOUBLE
  call h5dopen_f(file_id, "/input/diagnostics/absreq", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/diagnostics/absreq'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%diagnostics%absreq, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/diagnostics/absreq'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/diagnostics/absreq'" ; goto 9998 ; endif

! /input/diagnostics/relreq --> s%input%diagnostics%relreq; rank=0; h5type=H5T_NATIVE_DOUBLE
  call h5dopen_f(file_id, "/input/diagnostics/relreq", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/diagnostics/relreq'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%diagnostics%relreq, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/diagnostics/relreq'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/diagnostics/relreq'" ; goto 9998 ; endif

! /input/diagnostics/absacc --> s%input%diagnostics%absacc; rank=0; h5type=H5T_NATIVE_DOUBLE
  call h5dopen_f(file_id, "/input/diagnostics/absacc", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/diagnostics/absacc'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%diagnostics%absacc, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/diagnostics/absacc'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/diagnostics/absacc'" ; goto 9998 ; endif

! /input/diagnostics/epsr --> s%input%diagnostics%epsr; rank=0; h5type=H5T_NATIVE_DOUBLE
  call h5dopen_f(file_id, "/input/diagnostics/epsr", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/diagnostics/epsr'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%diagnostics%epsr, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/diagnostics/epsr'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/diagnostics/epsr'" ; goto 9998 ; endif

! /input/diagnostics/nPpts --> s%input%diagnostics%nPpts; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/diagnostics/nPpts", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/diagnostics/nPpts'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%diagnostics%nPpts, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/diagnostics/nPpts'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/diagnostics/nPpts'" ; goto 9998 ; endif

! /input/diagnostics/nPtrj --> s%input%diagnostics%nPtrj(1:dims_1(1)); rank=1
  call h5dopen_f(file_id, "/input/diagnostics/nPtrj", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/diagnostics/nPtrj'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/input/diagnostics/nPtrj'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_1, max_dims_1, hdfier)
  if (hdfier.ne.1) then ; write(*,*) "unexpected rank of dataset '/input/diagnostics/nPtrj': ",hdfier," .ne. 1" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/input/diagnostics/nPtrj'" ; goto 9998 ; endif
  
  allocate(s%input%diagnostics%nPtrj(1:dims_1(1)))
  
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%diagnostics%nPtrj(1:dims_1(1)), dims_1, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/diagnostics/nPtrj'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/diagnostics/nPtrj'" ; goto 9998 ; endif

! /input/diagnostics/LHevalues --> s%input%diagnostics%LHevalues; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/diagnostics/LHevalues", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/diagnostics/LHevalues'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, logical_tmp, int((/1/), HSIZE_T), hdfier)
  s%input%diagnostics%LHevalues = merge(.TRUE., .FALSE., logical_tmp.ne.0)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/diagnostics/LHevalues'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/diagnostics/LHevalues'" ; goto 9998 ; endif

! /input/diagnostics/LHevectors --> s%input%diagnostics%LHevectors; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/diagnostics/LHevectors", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/diagnostics/LHevectors'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, logical_tmp, int((/1/), HSIZE_T), hdfier)
  s%input%diagnostics%LHevectors = merge(.TRUE., .FALSE., logical_tmp.ne.0)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/diagnostics/LHevectors'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/diagnostics/LHevectors'" ; goto 9998 ; endif

! /input/diagnostics/LHmatrix --> s%input%diagnostics%LHmatrix; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/diagnostics/LHmatrix", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/diagnostics/LHmatrix'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, logical_tmp, int((/1/), HSIZE_T), hdfier)
  s%input%diagnostics%LHmatrix = merge(.TRUE., .FALSE., logical_tmp.ne.0)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/diagnostics/LHmatrix'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/diagnostics/LHmatrix'" ; goto 9998 ; endif

! /input/diagnostics/Lperturbed --> s%input%diagnostics%Lperturbed; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/diagnostics/Lperturbed", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/diagnostics/Lperturbed'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%diagnostics%Lperturbed, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/diagnostics/Lperturbed'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/diagnostics/Lperturbed'" ; goto 9998 ; endif

! /input/diagnostics/dpp --> s%input%diagnostics%dpp; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/diagnostics/dpp", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/diagnostics/dpp'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%diagnostics%dpp, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/diagnostics/dpp'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/diagnostics/dpp'" ; goto 9998 ; endif

! /input/diagnostics/dqq --> s%input%diagnostics%dqq; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/diagnostics/dqq", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/diagnostics/dqq'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%diagnostics%dqq, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/diagnostics/dqq'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/diagnostics/dqq'" ; goto 9998 ; endif

! /input/diagnostics/Lcheck --> s%input%diagnostics%Lcheck; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/diagnostics/Lcheck", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/diagnostics/Lcheck'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%input%diagnostics%Lcheck, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/diagnostics/Lcheck'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/diagnostics/Lcheck'" ; goto 9998 ; endif

! /input/diagnostics/Ltiming --> s%input%diagnostics%Ltiming; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/input/diagnostics/Ltiming", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/diagnostics/Ltiming'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, logical_tmp, int((/1/), HSIZE_T), hdfier)
  s%input%diagnostics%Ltiming = merge(.TRUE., .FALSE., logical_tmp.ne.0)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/diagnostics/Ltiming'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/diagnostics/Ltiming'" ; goto 9998 ; endif

! /input/diagnostics/fudge --> s%input%diagnostics%fudge; rank=0; h5type=H5T_NATIVE_DOUBLE
  call h5dopen_f(file_id, "/input/diagnostics/fudge", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/diagnostics/fudge'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%diagnostics%fudge, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/diagnostics/fudge'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/diagnostics/fudge'" ; goto 9998 ; endif

! /input/diagnostics/scaling --> s%input%diagnostics%scaling; rank=0; h5type=H5T_NATIVE_DOUBLE
  call h5dopen_f(file_id, "/input/diagnostics/scaling", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/input/diagnostics/scaling'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%input%diagnostics%scaling, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/diagnostics/scaling'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/diagnostics/scaling'" ; goto 9998 ; endif

! /output/Vns --> s%output%Vns(1:dims_1(1)); rank=1
  call h5dopen_f(file_id, "/output/Vns", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/output/Vns'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/output/Vns'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_1, max_dims_1, hdfier)
  if (hdfier.ne.1) then ; write(*,*) "unexpected rank of dataset '/output/Vns': ",hdfier," .ne. 1" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/output/Vns'" ; goto 9998 ; endif
  
  allocate(s%output%Vns(1:dims_1(1)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%output%Vns(1:dims_1(1)), dims_1, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/output/Vns'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/output/Vns'" ; goto 9998 ; endif

! /output/Bns --> s%output%Bns(1:dims_1(1)); rank=1
  call h5dopen_f(file_id, "/output/Bns", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/output/Bns'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/output/Bns'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_1, max_dims_1, hdfier)
  if (hdfier.ne.1) then ; write(*,*) "unexpected rank of dataset '/output/Bns': ",hdfier," .ne. 1" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/output/Bns'" ; goto 9998 ; endif
  
  allocate(s%output%Bns(1:dims_1(1)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%output%Bns(1:dims_1(1)), dims_1, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/output/Bns'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/output/Bns'" ; goto 9998 ; endif

! /output/Vnc --> s%output%Vnc(1:dims_1(1)); rank=1
  call h5dopen_f(file_id, "/output/Vnc", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/output/Vnc'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/output/Vnc'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_1, max_dims_1, hdfier)
  if (hdfier.ne.1) then ; write(*,*) "unexpected rank of dataset '/output/Vnc': ",hdfier," .ne. 1" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/output/Vnc'" ; goto 9998 ; endif
  
  allocate(s%output%Vnc(1:dims_1(1)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%output%Vnc(1:dims_1(1)), dims_1, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/output/Vnc'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/output/Vnc'" ; goto 9998 ; endif

! /output/Bnc --> s%output%Bnc(1:dims_1(1)); rank=1
  call h5dopen_f(file_id, "/output/Bnc", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/output/Bnc'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/output/Bnc'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_1, max_dims_1, hdfier)
  if (hdfier.ne.1) then ; write(*,*) "unexpected rank of dataset '/output/Bnc': ",hdfier," .ne. 1" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/output/Bnc'" ; goto 9998 ; endif
  
  allocate(s%output%Bnc(1:dims_1(1)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%output%Bnc(1:dims_1(1)), dims_1, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/output/Bnc'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/output/Bnc'" ; goto 9998 ; endif

! /output/mn --> s%output%mn; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/output/mn", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/output/mn'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%output%mn, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/output/mn'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/output/mn'" ; goto 9998 ; endif

! /output/im --> s%output%im(1:dims_1(1)); rank=1
  call h5dopen_f(file_id, "/output/im", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/output/im'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/output/im'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_1, max_dims_1, hdfier)
  if (hdfier.ne.1) then ; write(*,*) "unexpected rank of dataset '/output/im': ",hdfier," .ne. 1" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/output/im'" ; goto 9998 ; endif
  
  allocate(s%output%im(1:dims_1(1)))
  
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%output%im(1:dims_1(1)), dims_1, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/output/im'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/output/im'" ; goto 9998 ; endif

! /output/in --> s%output%in(1:dims_1(1)); rank=1
  call h5dopen_f(file_id, "/output/in", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/output/in'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/output/in'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_1, max_dims_1, hdfier)
  if (hdfier.ne.1) then ; write(*,*) "unexpected rank of dataset '/output/in': ",hdfier," .ne. 1" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/output/in'" ; goto 9998 ; endif
  
  allocate(s%output%in(1:dims_1(1)))
  
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%output%in(1:dims_1(1)), dims_1, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/output/in'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/output/in'" ; goto 9998 ; endif

! /output/Mvol --> s%output%Mvol; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/output/Mvol", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/output/Mvol'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%output%Mvol, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/output/Mvol'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/output/Mvol'" ; goto 9998 ; endif

! /output/Rbc --> s%output%Rbc(1:dims_2(1), 1:dims_2(2)); rank=2
  call h5dopen_f(file_id, "/output/Rbc", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/output/Rbc'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/output/Rbc'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_2, max_dims_2, hdfier)
  if (hdfier.ne.2) then ; write(*,*) "unexpected rank of dataset '/output/Rbc': ",hdfier," .ne. 2" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/output/Rbc'" ; goto 9998 ; endif
  
  allocate(s%output%Rbc(1:dims_2(1), 1:dims_2(2)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%output%Rbc(1:dims_2(1), 1:dims_2(2)), dims_2, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/output/Rbc'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/output/Rbc'" ; goto 9998 ; endif

! /output/Zbs --> s%output%Zbs(1:dims_2(1), 1:dims_2(2)); rank=2
  call h5dopen_f(file_id, "/output/Zbs", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/output/Zbs'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/output/Zbs'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_2, max_dims_2, hdfier)
  if (hdfier.ne.2) then ; write(*,*) "unexpected rank of dataset '/output/Zbs': ",hdfier," .ne. 2" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/output/Zbs'" ; goto 9998 ; endif
  
  allocate(s%output%Zbs(1:dims_2(1), 1:dims_2(2)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%output%Zbs(1:dims_2(1), 1:dims_2(2)), dims_2, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/output/Zbs'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/output/Zbs'" ; goto 9998 ; endif

! /output/Rbs --> s%output%Rbs(1:dims_2(1), 1:dims_2(2)); rank=2
  call h5dopen_f(file_id, "/output/Rbs", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/output/Rbs'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/output/Rbs'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_2, max_dims_2, hdfier)
  if (hdfier.ne.2) then ; write(*,*) "unexpected rank of dataset '/output/Rbs': ",hdfier," .ne. 2" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/output/Rbs'" ; goto 9998 ; endif
  
  allocate(s%output%Rbs(1:dims_2(1), 1:dims_2(2)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%output%Rbs(1:dims_2(1), 1:dims_2(2)), dims_2, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/output/Rbs'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/output/Rbs'" ; goto 9998 ; endif

! /output/Zbc --> s%output%Zbc(1:dims_2(1), 1:dims_2(2)); rank=2
  call h5dopen_f(file_id, "/output/Zbc", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/output/Zbc'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/output/Zbc'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_2, max_dims_2, hdfier)
  if (hdfier.ne.2) then ; write(*,*) "unexpected rank of dataset '/output/Zbc': ",hdfier," .ne. 2" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/output/Zbc'" ; goto 9998 ; endif
  
  allocate(s%output%Zbc(1:dims_2(1), 1:dims_2(2)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%output%Zbc(1:dims_2(1), 1:dims_2(2)), dims_2, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/output/Zbc'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/output/Zbc'" ; goto 9998 ; endif

! /output/ForceErr --> s%output%ForceErr; rank=0; h5type=H5T_NATIVE_DOUBLE
  call h5dopen_f(file_id, "/output/ForceErr", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/output/ForceErr'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%output%ForceErr, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/output/ForceErr'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/output/ForceErr'" ; goto 9998 ; endif

! /output/adiabatic --> s%output%adiabatic(1:dims_1(1)); rank=1
  call h5dopen_f(file_id, "/output/adiabatic", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/output/adiabatic'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/output/adiabatic'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_1, max_dims_1, hdfier)
  if (hdfier.ne.1) then ; write(*,*) "unexpected rank of dataset '/output/adiabatic': ",hdfier," .ne. 1" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/output/adiabatic'" ; goto 9998 ; endif
  
  allocate(s%output%adiabatic(1:dims_1(1)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%output%adiabatic(1:dims_1(1)), dims_1, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/output/adiabatic'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/output/adiabatic'" ; goto 9998 ; endif

! /output/helicity --> s%output%helicity(1:dims_1(1)); rank=1
  call h5dopen_f(file_id, "/output/helicity", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/output/helicity'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/output/helicity'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_1, max_dims_1, hdfier)
  if (hdfier.ne.1) then ; write(*,*) "unexpected rank of dataset '/output/helicity': ",hdfier," .ne. 1" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/output/helicity'" ; goto 9998 ; endif
  
  allocate(s%output%helicity(1:dims_1(1)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%output%helicity(1:dims_1(1)), dims_1, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/output/helicity'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/output/helicity'" ; goto 9998 ; endif

! /output/mu --> s%output%mu(1:dims_1(1)); rank=1
  call h5dopen_f(file_id, "/output/mu", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/output/mu'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/output/mu'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_1, max_dims_1, hdfier)
  if (hdfier.ne.1) then ; write(*,*) "unexpected rank of dataset '/output/mu': ",hdfier," .ne. 1" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/output/mu'" ; goto 9998 ; endif
  
  allocate(s%output%mu(1:dims_1(1)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%output%mu(1:dims_1(1)), dims_1, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/output/mu'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/output/mu'" ; goto 9998 ; endif

! /output/tflux --> s%output%tflux(1:dims_1(1)); rank=1
  call h5dopen_f(file_id, "/output/tflux", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/output/tflux'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/output/tflux'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_1, max_dims_1, hdfier)
  if (hdfier.ne.1) then ; write(*,*) "unexpected rank of dataset '/output/tflux': ",hdfier," .ne. 1" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/output/tflux'" ; goto 9998 ; endif
  
  allocate(s%output%tflux(1:dims_1(1)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%output%tflux(1:dims_1(1)), dims_1, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/output/tflux'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/output/tflux'" ; goto 9998 ; endif

! /output/pflux --> s%output%pflux(1:dims_1(1)); rank=1
  call h5dopen_f(file_id, "/output/pflux", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/output/pflux'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/output/pflux'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_1, max_dims_1, hdfier)
  if (hdfier.ne.1) then ; write(*,*) "unexpected rank of dataset '/output/pflux': ",hdfier," .ne. 1" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/output/pflux'" ; goto 9998 ; endif
  
  allocate(s%output%pflux(1:dims_1(1)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%output%pflux(1:dims_1(1)), dims_1, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/output/pflux'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/output/pflux'" ; goto 9998 ; endif

! /output/volume --> s%output%volume; rank=0; h5type=H5T_NATIVE_DOUBLE
  call h5dopen_f(file_id, "/output/volume", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/output/volume'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%output%volume, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/output/volume'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/output/volume'" ; goto 9998 ; endif

! /output/Mrad --> s%output%Mrad; rank=0; h5type=H5T_NATIVE_INTEGER
  call h5dopen_f(file_id, "/output/Mrad", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/output/Mrad'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, s%output%Mrad, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/output/Mrad'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/output/Mrad'" ; goto 9998 ; endif

! /output/TT --> s%output%TT(1:dims_3(1), 1:dims_3(2), 1:dims_3(3)); rank=3
  call h5dopen_f(file_id, "/output/TT", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/output/TT'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/output/TT'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_3, max_dims_3, hdfier)
  if (hdfier.ne.3) then ; write(*,*) "unexpected rank of dataset '/output/TT': ",hdfier," .ne. 3" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/output/TT'" ; goto 9998 ; endif
  
  allocate(s%output%TT(1:dims_3(1), 1:dims_3(2), 1:dims_3(3)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%output%TT(1:dims_3(1), 1:dims_3(2), 1:dims_3(3)), dims_3, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/output/TT'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/output/TT'" ; goto 9998 ; endif

! /output/Btemn --> s%output%Btemn(1:dims_3(1), 1:dims_3(2), 1:dims_3(3)); rank=3
  call h5dopen_f(file_id, "/output/Btemn", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/output/Btemn'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/output/Btemn'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_3, max_dims_3, hdfier)
  if (hdfier.ne.3) then ; write(*,*) "unexpected rank of dataset '/output/Btemn': ",hdfier," .ne. 3" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/output/Btemn'" ; goto 9998 ; endif
  
  allocate(s%output%Btemn(1:dims_3(1), 1:dims_3(2), 1:dims_3(3)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%output%Btemn(1:dims_3(1), 1:dims_3(2), 1:dims_3(3)), dims_3, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/output/Btemn'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/output/Btemn'" ; goto 9998 ; endif

! /output/Bzemn --> s%output%Bzemn(1:dims_3(1), 1:dims_3(2), 1:dims_3(3)); rank=3
  call h5dopen_f(file_id, "/output/Bzemn", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/output/Bzemn'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/output/Bzemn'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_3, max_dims_3, hdfier)
  if (hdfier.ne.3) then ; write(*,*) "unexpected rank of dataset '/output/Bzemn': ",hdfier," .ne. 3" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/output/Bzemn'" ; goto 9998 ; endif
  
  allocate(s%output%Bzemn(1:dims_3(1), 1:dims_3(2), 1:dims_3(3)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%output%Bzemn(1:dims_3(1), 1:dims_3(2), 1:dims_3(3)), dims_3, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/output/Bzemn'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/output/Bzemn'" ; goto 9998 ; endif

! /output/Btomn --> s%output%Btomn(1:dims_3(1), 1:dims_3(2), 1:dims_3(3)); rank=3
  call h5dopen_f(file_id, "/output/Btomn", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/output/Btomn'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/output/Btomn'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_3, max_dims_3, hdfier)
  if (hdfier.ne.3) then ; write(*,*) "unexpected rank of dataset '/output/Btomn': ",hdfier," .ne. 3" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/output/Btomn'" ; goto 9998 ; endif
  
  allocate(s%output%Btomn(1:dims_3(1), 1:dims_3(2), 1:dims_3(3)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%output%Btomn(1:dims_3(1), 1:dims_3(2), 1:dims_3(3)), dims_3, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/output/Btomn'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/output/Btomn'" ; goto 9998 ; endif

! /output/Bzomn --> s%output%Bzomn(1:dims_3(1), 1:dims_3(2), 1:dims_3(3)); rank=3
  call h5dopen_f(file_id, "/output/Bzomn", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/output/Bzomn'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/output/Bzomn'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_3, max_dims_3, hdfier)
  if (hdfier.ne.3) then ; write(*,*) "unexpected rank of dataset '/output/Bzomn': ",hdfier," .ne. 3" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/output/Bzomn'" ; goto 9998 ; endif
  
  allocate(s%output%Bzomn(1:dims_3(1), 1:dims_3(2), 1:dims_3(3)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%output%Bzomn(1:dims_3(1), 1:dims_3(2), 1:dims_3(3)), dims_3, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/output/Bzomn'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/output/Bzomn'" ; goto 9998 ; endif

! /output/lmns --> s%output%lmns; rank=0; h5type=H5T_NATIVE_DOUBLE
  call h5dopen_f(file_id, "/output/lmns", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/output/lmns'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%output%lmns, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/output/lmns'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/output/lmns'" ; goto 9998 ; endif

! /vector_potential/Ate --> s%vector_potential%Ate(1:dims_2(1), 1:dims_2(2)); rank=2
  call h5dopen_f(file_id, "/vector_potential/Ate", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/vector_potential/Ate'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/vector_potential/Ate'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_2, max_dims_2, hdfier)
  if (hdfier.ne.2) then ; write(*,*) "unexpected rank of dataset '/vector_potential/Ate': ",hdfier," .ne. 2" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/vector_potential/Ate'" ; goto 9998 ; endif
  
  allocate(s%vector_potential%Ate(1:dims_2(1), 1:dims_2(2)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%vector_potential%Ate(1:dims_2(1), 1:dims_2(2)), dims_2, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/vector_potential/Ate'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/vector_potential/Ate'" ; goto 9998 ; endif

! /vector_potential/Aze --> s%vector_potential%Aze(1:dims_2(1), 1:dims_2(2)); rank=2
  call h5dopen_f(file_id, "/vector_potential/Aze", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/vector_potential/Aze'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/vector_potential/Aze'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_2, max_dims_2, hdfier)
  if (hdfier.ne.2) then ; write(*,*) "unexpected rank of dataset '/vector_potential/Aze': ",hdfier," .ne. 2" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/vector_potential/Aze'" ; goto 9998 ; endif
  
  allocate(s%vector_potential%Aze(1:dims_2(1), 1:dims_2(2)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%vector_potential%Aze(1:dims_2(1), 1:dims_2(2)), dims_2, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/vector_potential/Aze'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/vector_potential/Aze'" ; goto 9998 ; endif

! /vector_potential/Ato --> s%vector_potential%Ato(1:dims_2(1), 1:dims_2(2)); rank=2
  call h5dopen_f(file_id, "/vector_potential/Ato", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/vector_potential/Ato'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/vector_potential/Ato'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_2, max_dims_2, hdfier)
  if (hdfier.ne.2) then ; write(*,*) "unexpected rank of dataset '/vector_potential/Ato': ",hdfier," .ne. 2" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/vector_potential/Ato'" ; goto 9998 ; endif
  
  allocate(s%vector_potential%Ato(1:dims_2(1), 1:dims_2(2)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%vector_potential%Ato(1:dims_2(1), 1:dims_2(2)), dims_2, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/vector_potential/Ato'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/vector_potential/Ato'" ; goto 9998 ; endif

! /vector_potential/Azo --> s%vector_potential%Azo(1:dims_2(1), 1:dims_2(2)); rank=2
  call h5dopen_f(file_id, "/vector_potential/Azo", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/vector_potential/Azo'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/vector_potential/Azo'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_2, max_dims_2, hdfier)
  if (hdfier.ne.2) then ; write(*,*) "unexpected rank of dataset '/vector_potential/Azo': ",hdfier," .ne. 2" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/vector_potential/Azo'" ; goto 9998 ; endif
  
  allocate(s%vector_potential%Azo(1:dims_2(1), 1:dims_2(2)))
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%vector_potential%Azo(1:dims_2(1), 1:dims_2(2)), dims_2, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/vector_potential/Azo'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/vector_potential/Azo'" ; goto 9998 ; endif

! /version --> s%version; rank=0; h5type=H5T_NATIVE_DOUBLE
  call h5dopen_f(file_id, "/version", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '/version'" ; goto 9998 ; endif
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, s%version, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/version'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/version'" ; goto 9998 ; endif

9998 continue
  
  call h5fclose_f(file_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing HDF5 file '",filename,"'" ; ierr = hdfier ; endif

9999 continue

  call h5close_f(hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing HDF5 library" ; ierr = hdfier ; endif 
    
end subroutine loadSpec

subroutine freeSpec(s)
  implicit none
  type(SpecOutput), intent(inout) :: s
  deallocate(s%input%physics%Lrad)
end subroutine freeSpec

end module read_spec

program test_read_spec
  use read_spec
  implicit none
  type(SpecOutput) :: s
  character(*), parameter :: filename = "/home/jonathan/Uni/04_PhD/00_programs/SPEC/SPEC/InputFiles/TestCases/G3V02L1Fi.001.h5"
  
  write(*,*) "reading '",filename,"'..."
  call loadSpec(s, filename)
  write(*,*) "done"
  
  write(*,"(A,F4.2)") "SPEC version: ", s%version
  write(*,"(A,99I2)") "Lrad:", s%input%physics%Lrad
  
  call freeSpec(s)
end program test_read_spec
