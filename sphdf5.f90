!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (output) ! Writes all the output information to ext.h5.

!latex \briefly{All the input and output information is contained in \type{ext.h5}.}
!latex \calledby{\link{xspech}}
!l tex \calls{\link{}}

!latex \tableofcontents

!latex \newcommand{\pb}[1]{\parbox[t]{13cm}{#1}}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

module sphdf5

  use inputlist , only : ext
  use fileunits , only : ounit
  use allglobal , only : myid
  use hdf5

  implicit none

  integer                        :: hdfier, rank                               ! error flag for HDF5 library
  integer(hid_t)                 :: file_id, space_id, dset_id                 ! default IDs used in macros
  integer(hsize_t)               :: onedims(1:1), twodims(1:2), threedims(1:3) ! dimension specifiers used in macros
  integer(size_t)                :: obj_count                                  ! number of open HDF5 objects

  integer(hid_t)                 :: iteration_dset_id                          ! Dataset identifier for "iteration"
  integer(hid_t)                 :: dataspace                                  ! dataspace for extension by 1 iteration object
  integer(hid_t)                 :: memspace                                   ! memspace for extension by 1 iteration object
  integer(hsize_t), dimension(1) :: old_data_dims, data_dims
  integer(hid_t)                 :: plist_id                                   ! Property list identifier used to activate dataset transfer property
  integer(hid_t)                 :: dt_nDcalls_id                              ! Memory datatype identifier (for "nDcalls" field)
  integer(hid_t)                 :: dt_Energy_id                               ! Memory datatype identifier (for "Energy" field)
  integer(hid_t)                 :: dt_ForceErr_id                             ! Memory datatype identifier (for "ForceErr" field)
  integer(hid_t)                 :: dt_iRbc_id                                 ! Memory datatype identifier (for "iRbc" field)
  integer(hid_t)                 :: dt_iZbs_id                                 ! Memory datatype identifier (for "iZbs" field)
  integer(hid_t)                 :: dt_iRbs_id                                 ! Memory datatype identifier (for "iRbs" field)
  integer(hid_t)                 :: dt_iZbc_id                                 ! Memory datatype identifier (for "iZbc" field)

contains

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
subroutine init_outfile

  LOCALS
  integer(hid_t) :: plist_id      ! Property list identifier used to activate MPI I/O parallel access to HDF5 library

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  ! initialize Fortran interface to the HDF5 library;
  call h5open_f( hdfier )
  FATAL( sphdf5, hdfier.ne.0, error calling h5open_f )

  ! Create file access property list to be able to tell HDF5 about MPI I/O
  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdfier)
  FATAL( sphdf5, hdfier.ne.0, h5pcreate_f returned an error )

  ! enable MPI I/O for parallel I/O access in file access property list
  call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, hdfier)
  FATAL( sphdf5, hdfier.ne.0, h5pset_fapl_mpio_f returned an error )

  ! Create the file (collectively, since called by all processes)
  call h5fcreate_f(trim(ext)//".h5", H5F_ACC_TRUNC_F, file_id, hdfier, access_prp = plist_id)
  FATAL( sphdf5, hdfier.ne.0, h5fcreate_f returned an error )

  ! file access property list is not needed after been used to specify MPI I/O during opening of file
  call h5pclose_f(plist_id, hdfier)
  FATAL( sphdf5, hdfier.ne.0, h5pclose_f returned an error )

end subroutine init_outfile

! mirror input variables into output file
subroutine mirror_input_to_outfile

  use inputlist
  use allglobal , only : Mvol

  LOCALS

  integer(hid_t) :: grpInput
  integer(hid_t) :: grpInputPhysics, grpInputNumerics, grpInputLocal, grpInputGlobal, grpInputDiagnostics

  HDEFGRP( file_id, input, grpInput )

! the following variables constitute the namelist/physicslist/; note that all variables in namelist need to be broadcasted in readin;
! they go into ext.h5/input/physics

  HDEFGRP( grpInput, physics, grpInputPhysics )

  HWRITEIV( grpInputPhysics,           1, Igeometry         , (/ Igeometry      /))
  HWRITEIV( grpInputPhysics,           1, Istellsym         , (/ Istellsym      /))
  HWRITEIV( grpInputPhysics,           1, Lfreebound        , (/ Lfreebound     /))
  HWRITERV( grpInputPhysics,           1, phiedge           , (/ phiedge        /))
  HWRITERV( grpInputPhysics,           1, curtor            , (/ curtor         /))
  HWRITERV( grpInputPhysics,           1, curpol            , (/ curpol         /))
  HWRITERV( grpInputPhysics,           1, gamma             , (/ gamma          /))
  HWRITEIV( grpInputPhysics,           1, Nfp               , (/ Nfp            /))
  HWRITEIV( grpInputPhysics,           1, Nvol              , (/ Nvol           /))
  HWRITEIV( grpInputPhysics,           1, Mpol              , (/ Mpol           /))
  HWRITEIV( grpInputPhysics,           1, Ntor              , (/ Ntor           /))
  HWRITEIV( grpInputPhysics,        Mvol, Lrad              ,      Lrad(1:Mvol)   )
  HWRITEIV( grpInputPhysics,           1, Lconstraint       , (/ Lconstraint    /))
  HWRITERV( grpInputPhysics,        Mvol, tflux             ,     tflux(1:Mvol)   )
  HWRITERV( grpInputPhysics,        Mvol, pflux             ,     pflux(1:Mvol)   )
  HWRITERV( grpInputPhysics,        Nvol, helicity          ,  helicity(1:Nvol)   )
  HWRITERV( grpInputPhysics,           1, pscale            , (/ pscale         /))
  HWRITERV( grpInputPhysics,        Nvol, pressure          ,  pressure(1:Nvol)   )
  HWRITEIV( grpInputPhysics,           1, Ladiabatic        , (/ Ladiabatic     /))
  HWRITERV( grpInputPhysics,        Mvol, adiabatic         , adiabatic(1:Nvol)   )
  HWRITERV( grpInputPhysics,    (1+Nvol), mu                ,        mu(1:Mvol)   )
  HWRITEIV( grpInputPhysics,    (1+Mvol), pl                ,        pl(0:Nvol)   )
  HWRITEIV( grpInputPhysics,    (1+Mvol), ql                ,        ql(0:Nvol)   )
  HWRITEIV( grpInputPhysics,    (1+Mvol), pr                ,        pr(0:Nvol)   )
  HWRITEIV( grpInputPhysics,    (1+Mvol), qr                ,        qr(0:Nvol)   )
  HWRITERV( grpInputPhysics,    (1+Nvol), iota              ,      iota(0:Nvol)   )
  HWRITEIV( grpInputPhysics,    (1+Mvol), lp                ,        lp(0:Nvol)   )
  HWRITEIV( grpInputPhysics,    (1+Mvol), lq                ,        lq(0:Nvol)   )
  HWRITEIV( grpInputPhysics,    (1+Mvol), rp                ,        rp(0:Nvol)   )
  HWRITEIV( grpInputPhysics,    (1+Mvol), rq                ,        rq(0:Nvol)   )
  HWRITERV( grpInputPhysics,    (1+Nvol), oita              ,      oita(0:Nvol)   )

  HWRITERV( grpInputPhysics,    (1+Ntor), Rac               ,       Rac(0:Ntor)   ) !     stellarator symmetric coordinate axis;
  HWRITERV( grpInputPhysics,    (1+Ntor), Zas               ,       Zas(0:Ntor)   )
  HWRITERV( grpInputPhysics,    (1+Ntor), Ras               ,       Ras(0:Ntor)   ) ! non-stellarator symmetric coordinate axis;
  HWRITERV( grpInputPhysics,    (1+Ntor), Zac               ,       Zac(0:Ntor)   )

  HWRITERA( grpInputPhysics,   (2*Ntor+1), (2*Mpol+1),  Rbc, Rbc(-Ntor:Ntor,-Mpol:Mpol) ) !     stellarator symmetric boundary components;
  HWRITERA( grpInputPhysics,   (2*Ntor+1), (2*Mpol+1),  Zbs, Zbs(-Ntor:Ntor,-Mpol:Mpol) ) !     stellarator symmetric boundary components;
  HWRITERA( grpInputPhysics,   (2*Ntor+1), (2*Mpol+1),  Rbs, Rbs(-Ntor:Ntor,-Mpol:Mpol) ) ! non-stellarator symmetric boundary components;
  HWRITERA( grpInputPhysics,   (2*Ntor+1), (2*Mpol+1),  Zbc, Zbc(-Ntor:Ntor,-Mpol:Mpol) ) ! non-stellarator symmetric boundary components;

  HWRITERA( grpInputPhysics,   (2*Ntor+1), (2*Mpol+1),  Rwc, Rwc(-Ntor:Ntor,-Mpol:Mpol) ) !     stellarator symmetric boundary components of wall;
  HWRITERA( grpInputPhysics,   (2*Ntor+1), (2*Mpol+1),  Zws, Zws(-Ntor:Ntor,-Mpol:Mpol) ) !     stellarator symmetric boundary components of wall;
  HWRITERA( grpInputPhysics,   (2*Ntor+1), (2*Mpol+1),  Rws, Rws(-Ntor:Ntor,-Mpol:Mpol) ) ! non-stellarator symmetric boundary components of wall;
  HWRITERA( grpInputPhysics,   (2*Ntor+1), (2*Mpol+1),  Zwc, Zwc(-Ntor:Ntor,-Mpol:Mpol) ) ! non-stellarator symmetric boundary components of wall;

  HWRITERA( grpInputPhysics,   (2*Ntor+1), (2*Mpol+1),  Vns, Vns(-Ntor:Ntor,-Mpol:Mpol) ) !     stellarator symmetric normal field at boundary; vacuum component;
  HWRITERA( grpInputPhysics,   (2*Ntor+1), (2*Mpol+1),  Bns, Bns(-Ntor:Ntor,-Mpol:Mpol) ) !     stellarator symmetric normal field at boundary; plasma component;
  HWRITERA( grpInputPhysics,   (2*Ntor+1), (2*Mpol+1),  Vnc, Vnc(-Ntor:Ntor,-Mpol:Mpol) ) ! non-stellarator symmetric normal field at boundary; vacuum component;
  HWRITERA( grpInputPhysics,   (2*Ntor+1), (2*Mpol+1),  Bnc, Bnc(-Ntor:Ntor,-Mpol:Mpol) ) ! non-stellarator symmetric normal field at boundary; plasma component;

  HWRITERV( grpInputPhysics,           1, mupftol           , (/ mupftol        /))
  HWRITEIV( grpInputPhysics,           1, mupfits           , (/ mupfits        /))

  HCLOSEGRP( grpInputPhysics )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! the following variables constitute the namelist/numericlist/; note that all variables in namelist need to be broadcasted in readin;
! they go into ext.h5/input/numerics

  HDEFGRP( grpInput, numerics, grpInputNumerics )

  HWRITEIV( grpInputNumerics,          1, Linitialize        , (/ Linitialize /))
  HWRITEIV( grpInputNumerics,          1, Lzerovac           , (/ Lzerovac    /))
  HWRITEIV( grpInputNumerics,          1, Ndiscrete          , (/ Ndiscrete   /))
  HWRITEIV( grpInputNumerics,          1, Nquad              , (/ Nquad       /))
  HWRITEIV( grpInputNumerics,          1, iMpol              , (/ iMpol       /))
  HWRITEIV( grpInputNumerics,          1, iNtor              , (/ iNtor       /))
  HWRITEIV( grpInputNumerics,          1, Lsparse            , (/ Lsparse     /))
  HWRITEIV( grpInputNumerics,          1, Lsvdiota           , (/ Lsvdiota    /))
  HWRITEIV( grpInputNumerics,          1, imethod            , (/ imethod     /))
  HWRITEIV( grpInputNumerics,          1, iorder             , (/ iorder      /))
  HWRITEIV( grpInputNumerics,          1, iprecon            , (/ iprecon     /))
  HWRITERV( grpInputNumerics,          1, iotatol            , (/ iotatol     /))
  HWRITEIV( grpInputNumerics,          1, Lextrap            , (/ Lextrap     /))
  HWRITEIV( grpInputNumerics,          1, Mregular           , (/ Mregular    /))

  HCLOSEGRP( grpInputNumerics )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! the following variables constitute the namelist/locallist/; note that all variables in namelist need to be broadcasted in readin;
! they go into ext.h5/input/local

  HDEFGRP( grpInput, local, grpInputLocal )

  HWRITEIV( grpInputLocal,             1, LBeltrami          , (/ LBeltrami   /))
  HWRITEIV( grpInputLocal,             1, Linitgues          , (/ Linitgues   /))
  HWRITEIV( grpInputLocal,             1, Lposdef            , (/ Lposdef     /)) ! redundant;
  HWRITERV( grpInputLocal,             1, maxrndgues         , (/ maxrndgues  /))

  HCLOSEGRP( grpInputLocal )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! the following variables constitute the namelist/globallist/; note that all variables in namelist need to be broadcasted in readin;
! they go into ext.h5/input/global

  HDEFGRP( grpInput, global, grpInputGlobal )

  HWRITEIV( grpInputGlobal,            1,  Lfindzero         , (/ Lfindzero   /))
  HWRITERV( grpInputGlobal,            1,  escale            , (/ escale      /))
  HWRITERV( grpInputGlobal,            1,  opsilon           , (/ opsilon     /))
  HWRITERV( grpInputGlobal,            1,  pcondense         , (/ pcondense   /))
  HWRITERV( grpInputGlobal,            1,  epsilon           , (/ epsilon     /))
  HWRITERV( grpInputGlobal,            1,  wpoloidal         , (/ wpoloidal   /))
  HWRITERV( grpInputGlobal,            1,  upsilon           , (/ upsilon     /))
  HWRITERV( grpInputGlobal,            1,  forcetol          , (/ forcetol    /))
  HWRITERV( grpInputGlobal,            1,  c05xmax           , (/ c05xmax     /))
  HWRITERV( grpInputGlobal,            1,  c05xtol           , (/ c05xtol     /))
  HWRITERV( grpInputGlobal,            1,  c05factor         , (/ c05factor   /))
  HWRITELV( grpInputGlobal,            1,  LreadGF           , (/ LreadGF     /))
  HWRITEIV( grpInputGlobal,            1,  mfreeits          , (/ mfreeits    /))
  HWRITERV( grpInputGlobal,            1,  bnstol            , (/ bnstol      /))  ! redundant;
  HWRITERV( grpInputGlobal,            1,  bnsblend          , (/ bnsblend    /))  ! redundant;
  HWRITERV( grpInputGlobal,            1,  gBntol            , (/ gBntol      /))
  HWRITERV( grpInputGlobal,            1,  gBnbld            , (/ gBnbld      /))
  HWRITERV( grpInputGlobal,            1,  vcasingeps        , (/ vcasingeps  /))
  HWRITERV( grpInputGlobal,            1,  vcasingtol        , (/ vcasingtol  /))
  HWRITEIV( grpInputGlobal,            1,  vcasingits        , (/ vcasingits  /))
  HWRITEIV( grpInputGlobal,            1,  vcasingper        , (/ vcasingper  /))
  HWRITEIV( grpInputGlobal,            1,  mcasingcal        , (/ mcasingcal  /))  ! redundant;

  HCLOSEGRP( grpInputGlobal )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! the following variables constitute the namelist/diagnosticslist/; note that all variables in namelist need to be broadcasted in readin;
! they go into ext.h5/input/diagnostics

  HDEFGRP( grpInput, diagnostics, grpInputDiagnostics )

  HWRITERV( grpInputDiagnostics,       1,  odetol            , (/ odetol         /))
  HWRITERV( grpInputDiagnostics,       1,  absreq            , (/ absreq         /))           ! redundant;
  HWRITERV( grpInputDiagnostics,       1,  relreq            , (/ relreq         /))           ! redundant;
  HWRITERV( grpInputDiagnostics,       1,  absacc            , (/ absacc         /))           ! redundant;
  HWRITERV( grpInputDiagnostics,       1,  epsr              , (/ epsr           /))           ! redundant;
  HWRITEIV( grpInputDiagnostics,       1,  nPpts             , (/ nPpts          /))
  HWRITEIV( grpInputDiagnostics,       1,  nPtrj             ,    nPtrj(1:Nvol+1)  )
  HWRITELV( grpInputDiagnostics,       1,  LHevalues         , (/ LHevalues      /))
  HWRITELV( grpInputDiagnostics,       1,  LHevectors        , (/ LHevectors     /))
  HWRITELV( grpInputDiagnostics,       1,  LHmatrix          , (/ LHmatrix       /))
  HWRITEIV( grpInputDiagnostics,       1,  Lperturbed        , (/ Lperturbed     /))
  HWRITEIV( grpInputDiagnostics,       1,  dpp               , (/ dpp            /))
  HWRITEIV( grpInputDiagnostics,       1,  dqq               , (/ dqq            /))
  HWRITEIV( grpInputDiagnostics,       1,  Lcheck            , (/ Lcheck         /))
  HWRITELV( grpInputDiagnostics,       1,  Ltiming           , (/ Ltiming        /))
  HWRITERV( grpInputDiagnostics,       1,  fudge             , (/ fudge          /))         ! redundant;
  HWRITERV( grpInputDiagnostics,       1,  scaling           , (/ scaling        /))          ! redundant;

  HCLOSEGRP( grpInputDiagnostics )

  HCLOSEGRP( grpInput )

end subroutine mirror_input_to_outfile

! prepare ``convergence evolution'' output
subroutine init_convergence_output

  use allglobal, only : mn, Mvol

  LOCALS
  integer(hid_t)                    :: iteration_dspace_id                           ! dataspace for "iteration"
  integer(hid_t)                    :: iteration_dtype_id                            ! Compound datatype for "iteration"
  integer(hid_t)                    :: iRZbscArray_id                                ! Memory datatype identifier
  integer(size_t)                   :: iteration_dtype_size                          ! Size of the "iteration" datatype
  integer(size_t)                   :: type_size_i                                   ! Size of the integer datatype
  integer(size_t)                   :: type_size_d                                   ! Size of the double precision datatype
  integer(size_t)                   :: offset                                        ! Member's offset
  integer(hid_t)                    :: crp_list                                      ! Dataset creation property identifier
  integer, parameter                :: rank = 1                                      ! logging rank: convergence logging is one-dimensional
  integer(hsize_t), dimension(rank) :: maxdims                                       ! convergence logging maximum dimensions => will be unlimited
  integer(hsize_t), dimension(rank) :: dims = (/ 0 /)                                ! current convergence logging dimensions
  integer(hsize_t), dimension(rank) :: dimsc = (/ 1 /)                               ! chunking length ???
  integer(size_t)                   :: irbc_size_template                            ! size ofiRbc array in iterations logging
  integer(size_t)                   :: irbc_size                                     ! size ofiRbc array in iterations logging

  ! Set dataset transfer property to preserve partially initialized fields
  ! during write/read to/from dataset with compound datatype.
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdfier)
  FATAL( sphdf5, hdfier.ne.0, error calling h5pcreate_f )

  call h5pset_preserve_f(plist_id, .TRUE., hdfier)
  FATAL( sphdf5, hdfier.ne.0, error calling h5pset_preserve_f )

  maxdims = (/ H5S_UNLIMITED_F /)                                                       ! unlimited array size: "converge" until you get bored
  call h5screate_simple_f(rank, dims, iteration_dspace_id, hdfier, maxdims)              ! Create the dataspace with zero initial size and allow it to grow
  FATAL( sphdf5, hdfier.ne.0, error calling h5screate_simple_f )

  call h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, hdfier) ! dataset creation property list with chunking
  FATAL( sphdf5, hdfier.ne.0, error calling h5pcreate_f )

  call h5pset_chunk_f(crp_list, rank, dimsc, hdfier)
  FATAL( sphdf5, hdfier.ne.0, error calling h5pset_chunk_f )

  ! declare "iteration" compound datatype
  ! declare array parts
  call h5tarray_create_f(H5T_NATIVE_DOUBLE, 2, int((/mn, Mvol+1/),hsize_t), iRZbscArray_id, hdfier) ! create array datatypes for i{R,Z}b{c,s}
  call h5tget_size_f(iRZbscArray_id, irbc_size, hdfier)
  call h5tget_size_f(H5T_NATIVE_INTEGER, type_size_i, hdfier)                            ! size of "iteration" field
  call h5tget_size_f(H5T_NATIVE_DOUBLE,  type_size_d, hdfier)                            ! size of "errorValue" field
  iteration_dtype_size = 2*type_size_i + 2*type_size_d + 4*irbc_size                     ! wflag, nDcalls, Energy, ForceErr, i{R,Z}b{c,s}

  call h5tcreate_f(H5T_COMPOUND_F, iteration_dtype_size, iteration_dtype_id, hdfier)     ! create compound datatype

  offset = 0                                                                             ! offset for first field starts at 0
  call h5tinsert_f(iteration_dtype_id, "nDcalls", offset, H5T_NATIVE_INTEGER, hdfier)      ! insert "nDcalls" field in datatype
  offset = offset + type_size_i                                                          ! increment offset by size of field
  call h5tinsert_f(iteration_dtype_id, "Energy", offset, H5T_NATIVE_DOUBLE, hdfier)      ! insert "Energy" field in datatype
  offset = offset + type_size_d                                                          ! increment offset by size of field
  call h5tinsert_f(iteration_dtype_id, "ForceErr", offset, H5T_NATIVE_DOUBLE, hdfier)       ! insert "ForceErr" field in datatype
  offset = offset + type_size_d                                                          ! increment offset by size of field
  call h5tinsert_f(iteration_dtype_id, "iRbc", offset, iRZbscArray_id, hdfier)           ! insert "iRbc" field in datatype
  offset = offset + irbc_size                                                            ! increment offset by size of field
  call h5tinsert_f(iteration_dtype_id, "iZbs", offset, iRZbscArray_id, hdfier)           ! insert "iZbs" field in datatype
  offset = offset + irbc_size                                                            ! increment offset by size of field
  call h5tinsert_f(iteration_dtype_id, "iRbs", offset, iRZbscArray_id, hdfier)           ! insert "iRbs" field in datatype
  offset = offset + irbc_size                                                            ! increment offset by size of field
  call h5tinsert_f(iteration_dtype_id, "iZbc", offset, iRZbscArray_id, hdfier)           ! insert "iZbc" field in datatype
  offset = offset + irbc_size                                                            ! increment offset by size of field

  call h5dcreate_f(file_id, "iterations", iteration_dtype_id, iteration_dspace_id, &     ! create dataset with compound type
                   iteration_dset_id, hdfier, crp_list)

  call h5sclose_f(iteration_dspace_id, hdfier)                                           ! Terminate access to the data space (does not show up in obj_count below)
                                                                                         ! --> only needed for creation of dataset

  ! Create memory types. We have to create a compound datatype
  ! for each member we want to write.
  offset = 0
  call h5tcreate_f(H5T_COMPOUND_F, type_size_i, dt_nDcalls_id,  hdfier)
  call h5tcreate_f(H5T_COMPOUND_F, type_size_d, dt_Energy_id, hdfier)
  call h5tcreate_f(H5T_COMPOUND_F, type_size_d, dt_ForceErr_id,  hdfier)
  call h5tcreate_f(H5T_COMPOUND_F, irbc_size,   dt_iRbc_id,   hdfier)
  call h5tcreate_f(H5T_COMPOUND_F, irbc_size,   dt_iZbs_id,   hdfier)
  call h5tcreate_f(H5T_COMPOUND_F, irbc_size,   dt_iRbs_id,   hdfier)
  call h5tcreate_f(H5T_COMPOUND_F, irbc_size,   dt_iZbc_id,   hdfier)

  call h5tinsert_f(dt_nDcalls_id,   "nDcalls", offset, H5T_NATIVE_INTEGER, hdfier)
  call h5tinsert_f(dt_Energy_id,     "Energy", offset, H5T_NATIVE_DOUBLE,  hdfier)
  call h5tinsert_f(dt_ForceErr_id, "ForceErr", offset, H5T_NATIVE_DOUBLE,  hdfier)
  call h5tinsert_f(dt_iRbc_id,         "iRbc", offset, iRZbscArray_id,     hdfier)
  call h5tinsert_f(dt_iZbs_id,         "iZbs", offset, iRZbscArray_id,     hdfier)
  call h5tinsert_f(dt_iRbs_id,         "iRbs", offset, iRZbscArray_id,     hdfier)
  call h5tinsert_f(dt_iZbc_id,         "iZbc", offset, iRZbscArray_id,     hdfier)

  ! create memspace with size of compound object to append
  dims(1) = 1 ! only append one iteration at a time
  call h5screate_simple_f (rank, dims, memspace, hdfier)

  call h5pclose_f(crp_list, hdfier)
  call h5tclose_f(iteration_dtype_id, hdfier)                                            ! Terminate access to the datatype
  call h5tclose_f(iRZbscArray_id, hdfier)                                            ! Terminate access to the datatype

end subroutine init_convergence_output


! was in global.f90/wrtend for wflag.eq.-1 previously
subroutine write_convergence_output( nDcalls, ForceErr )

  use allglobal, only : mn, Mvol, Energy, iRbc, iZbs, iRbs, iZbc

  LOCALS
  INTEGER, intent(in)  :: nDcalls
  REAL   , intent(in)  :: ForceErr

#ifdef DEBUG
  if( Wsphdf5 ) then ; cput = GETTIME ; write(ounit,'("sphdf5 : ",f10.2," : myid=",i3," ; writing convergence ;")') cput-cpus, myid
  endif
#endif

  ! append updated values to "iterations" dataset

  ! open dataspace to current state of dataset
  call h5dget_space_f(iteration_dset_id, dataspace, hdfier)

  ! get current size of dataset
  call h5sget_simple_extent_npoints_f(dataspace, old_data_dims(1), hdfier)

  ! blow up dataset to new size
  data_dims = old_data_dims+1
  call h5dset_extent_f(iteration_dset_id, data_dims, hdfier)

  ! get dataspace slab corresponding to region which the iterations dataset was extended by
  call h5dget_space_f(iteration_dset_id, dataspace, hdfier) ! re-select dataspace to update size info in HDF5 lib
  call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, old_data_dims, (/ INT(1, HSIZE_T) /), hdfier)

  ! write next iteration object
  call h5dwrite_f(iteration_dset_id, dt_nDcalls_id, nDcalls, (/ INT(1, HSIZE_T)/), hdfier, &
    & mem_space_id=memspace, file_space_id=dataspace, xfer_prp=plist_id)
  call h5dwrite_f(iteration_dset_id, dt_Energy_id, Energy, (/ INT(1, HSIZE_T)/), hdfier, &
    & mem_space_id=memspace, file_space_id=dataspace, xfer_prp=plist_id)
  call h5dwrite_f(iteration_dset_id, dt_ForceErr_id, ForceErr, (/ INT(1, HSIZE_T)/), hdfier, &
    & mem_space_id=memspace, file_space_id=dataspace, xfer_prp=plist_id)
  call h5dwrite_f(iteration_dset_id, dt_iRbc_id, iRbc, INT((/mn,Mvol+1/), HSIZE_T), hdfier, &
    & mem_space_id=memspace, file_space_id=dataspace, xfer_prp=plist_id)
  call h5dwrite_f(iteration_dset_id, dt_iZbs_id, iZbs, INT((/mn,Mvol+1/), HSIZE_T), hdfier, &
    & mem_space_id=memspace, file_space_id=dataspace, xfer_prp=plist_id)
  call h5dwrite_f(iteration_dset_id, dt_iRbs_id, iRbs, INT((/mn,Mvol+1/), HSIZE_T), hdfier, &
    & mem_space_id=memspace, file_space_id=dataspace, xfer_prp=plist_id)
  call h5dwrite_f(iteration_dset_id, dt_iZbc_id, iZbc, INT((/mn,Mvol+1/), HSIZE_T), hdfier, &
    & mem_space_id=memspace, file_space_id=dataspace, xfer_prp=plist_id)

  ! dataspace to appended object should be closed now
  ! MAYBE we otherwise keep all the iterations in memory?
  call h5sclose_f(dataspace, hdfier)

end subroutine write_convergence_output


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
! final output
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
subroutine hdfint

  use fileunits, only : ounit
  use inputlist
  use allglobal, only : ncpu, cpus, &
                        Mvol, ForceErr, &
                        mn, im, in, iRbc, iZbs, iRbs, iZbc, &
                        dRbc, dZbs, dRbs, dZbc, &
                        vvolume, dvolume, &
                        Bsupumn, Bsupvmn, &
                        Btemn, Bzemn, Btomn, Bzomn, &
                        iVns, iBns, iVnc, iBnc, &
                        lmns, &
                        TT, &
                        beltramierror

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS

  INTEGER                        :: Mrad
  REAL                           :: tvolume

  integer(hid_t) :: grpOutput

  HDEFGRP( file_id, output, grpOutput )

  HWRITERV( grpOutput,           mn,               Vns,       iVns(1:mn)   ) !     stellarator symmetric normal field at boundary; vacuum component;
  HWRITERV( grpOutput,           mn,               Bns,       iBns(1:mn)   ) !     stellarator symmetric normal field at boundary; plasma component;
  HWRITERV( grpOutput,           mn,               Vnc,       iVnc(1:mn)   ) ! non-stellarator symmetric normal field at boundary; vacuum component;
  HWRITERV( grpOutput,           mn,               Bnc,       iBnc(1:mn)   ) ! non-stellarator symmetric normal field at boundary; plasma component;

!latex \begin{enumerate}

!latex \item In addition to the input variables, which are described in \link{global}, the following quantities are written to \type{ext.h5} :
!latex
!latex \begin{tabular}{|l|l|l|} \hline

!latex \type{variable}               & type    & \pb{description} \\ \hline

!latex \type{mn}                     & integer & \pb{number of Fourier modes} \\
  HWRITEIV( grpOutput, 1, mn, (/ mn /)  )
!latex \type{im(1:mn)}               & integer & \pb{poloidal mode numbers} \\
  HWRITEIV( grpOutput, mn, im, im(1:mn) )
!latex \type{in(1:mn)}               & integer & \pb{toroidal mode numbers} \\
  HWRITEIV( grpOutput, mn, in, in(1:mn) )
!latex \type{Mvol}                   & integer & \pb{number of interfaces = number of volumes} \\
  HWRITEIV( grpOutput, 1, Mvol, (/ Mvol /))
!latex \type{iRbc(1:mn,0:Mvol)}      & real    & \pb{Fourier harmonics, $R_{m,n}$, of interfaces} \\
  HWRITERA( grpOutput, mn, (Mvol+1), Rbc, iRbc(1:mn,0:Mvol) )
!latex \type{iZbs(1:mn,0:Mvol)}      & real    & \pb{Fourier harmonics, $Z_{m,n}$, of interfaces} \\
  HWRITERA( grpOutput, mn, (Mvol+1), Zbs, iZbs(1:mn,0:Mvol) )
!latex \type{iRbs(1:mn,0:Mvol)}      & real    & \pb{Fourier harmonics, $R_{m,n}$, of interfaces} \\
  HWRITERA( grpOutput, mn, (Mvol+1), Rbs, iRbs(1:mn,0:Mvol) )
!latex \type{iZbc(1:mn,0:Mvol)}      & real    & \pb{Fourier harmonics, $Z_{m,n}$, of interfaces} \\
  HWRITERA( grpOutput, mn, (Mvol+1), Zbc, iZbc(1:mn,0:Mvol) )
!latex \type{forcetol}               & real    & \pb{force-balance error across interfaces} \\
  HWRITERV( grpOutput, 1, forcetol, (/ forcetol /))
!latex \type{ForceErr}               & real    & \pb{force-balance error across interfaces} \\
  HWRITERV( grpOutput, 1, ForceErr, (/ ForceErr /))

  if( Lcheck.eq.1 ) then
!latex \type{beltramierror}          & real    & \pb{error in beltrami field (volume integral)} \\
   HWRITERA( grpOutput, Mvol, 3, beltramierror, beltramierror(1:Mvol,1:3) )
  endif

  if( allocated(vvolume) ) then ! why is it required to confirm that vvolume has been allocated ; 24 Nov 16;

   tvolume = sum(vvolume(1:Nvol) )
!latex \type{volume}                 & real    & \pb{total volume = $\sum V_v$} \\
   HWRITERV( grpOutput, 1, volume, (/ tvolume /))

  else

   if( Wsphdf5 ) write(ounit,'("hdfint : ", 10x ," : myid=",i3," ; vvolume is not allocated ;")') myid

  endif ! end of if( allocated(vvolume) ) ; 11 Aug 14;

  Mrad  = maxval( Lrad(1:Mvol) )
!latex \type{Mrad}                   & integer & \pb{the maximum radial (Chebyshev) resolution} \\
  HWRITEIV( grpOutput, 1, Mrad, (/ Mrad /))
!latex \type{TT(0:Mrad,0:1,0:1)}     & real    & \pb{the Chebyshev polynomials, $T_l$, and their derivatives, evaluated at $s=\pm 1$} \\
  HWRITERC( grpOutput, (Mrad+1), 2, 2, TT, TT(0:Mrad,0:1,0:1) )
!latex \type{Btemn(1:mn,0:1,1:Mvol)} & real    & \pb{the cosine harmonics of the covariant poloidal field, \\
!latex                                           i.e. $[[B_{\t,j}]]$ evaluated on the inner and outer interface in each volume} \\
  HWRITERC( grpOutput, mn, 2, Mvol, Btemn, Btemn(1:mn,0:1,1:Mvol) )
!latex \type{Bzemn(1:mn,0:1,1:Mvol)} & real    & \pb{the cosine harmonics of the covariant toroidal field, \\
!latex                                           i.e. $[[B_{\z,j}]]$ evaluated on the inner and outer interface in each volume} \\
  HWRITERC( grpOutput, mn, 2, Mvol, Bzemn, Bzemn(1:mn,0:1,1:Mvol) )
!latex \type{Btomn(1:mn,0:1,1:Mvol)} & real    & \pb{the sine harmonics of the covariant poloidal field, \\
!latex                                           i.e. $[[B_{\t,j}]]$ evaluated on the inner and outer interface in each volume} \\
  HWRITERC( grpOutput, mn, 2, Mvol, Btomn, Btomn(1:mn,0:1,1:Mvol) )
!latex \type{Bzomn(1:mn,0:1,1:Mvol)} & real    & \pb{the sine harmonics of the covariant toroidal field, \\
!latex                                           i.e. $[[B_{\z,j}]]$ evaluated on the inner and outer interface in each volume} \\
  HWRITERC( grpOutput, mn, 2, Mvol, Bzomn, Bzomn(1:mn,0:1,1:Mvol) )

  if( Lperturbed.eq.1 ) then

!latex \type{dRbc(1:mn,0:Nvol)}      & real    & \pb{Fourier harmonics, $R_{j}$, of interfaces; linearly perturbed solution} \\
  HWRITERA( grpOutput, mn, (Nvol+1), dRbc, dRbc(1:mn,0:Nvol) )
!latex \type{dZbs(1:mn,0:Nvol)}      & real    & \pb{Fourier harmonics, $Z_{j}$, of interfaces; linearly perturbed solution} \\
  HWRITERA( grpOutput, mn, (Nvol+1), dZbs, dZbs(1:mn,0:Nvol) )
!latex \type{dRbs(1:mn,0:Nvol)}      & real    & \pb{Fourier harmonics, $R_{j}$, of interfaces; linearly perturbed solution} \\
  HWRITERA( grpOutput, mn, (Nvol+1), dRbs, dRbs(1:mn,0:Nvol) )
!latex \type{dZbc(1:mn,0:Nvol)}      & real    & \pb{Fourier harmonics, $Z_{j}$, of interfaces; linearly perturbed solution} \\
  HWRITERA( grpOutput, mn, (Nvol+1), dZbc, dZbc(1:mn,0:Nvol) )

  endif

!latex \type{lmns}                   & integer & \pb{resolution of straight fieldline transformation} \\
  HWRITEIV( grpOutput, 1, lmns, (/ lmns /))

!latex \hline \end{tabular}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item All quantities marked as real should be treated as double precision.

!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  !RETURN(hdfint)
  HCLOSEGRP( grpOutput )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine hdfint

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
subroutine finish_outfile

  LOCALS

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  call h5tclose_f(dt_nDcalls_id, hdfier)
  call h5tclose_f(dt_Energy_id, hdfier)
  call h5tclose_f(dt_ForceErr_id, hdfier)
  call h5tclose_f(dt_iRbc_id, hdfier)
  call h5tclose_f(dt_iZbs_id, hdfier)
  call h5tclose_f(dt_iRbs_id, hdfier)
  call h5tclose_f(dt_iZbc_id, hdfier)

  call h5dclose_f(iteration_dset_id, hdfier)                                             ! End access to the dataset and release resources used by it.
  call h5pclose_f(plist_id, hdfier)                                                      ! close plist used for 'preserve' flag (does not show up in obj_count below)

  call h5fclose_f( file_id, hdfier ) ! terminate access on output file;
  FATAL( sphdf5, hdfier.ne.0, error calling h5fclose_f )

  call h5fget_obj_count_f(INT(H5F_OBJ_ALL_F,HID_T), H5F_OBJ_ALL_F, obj_count, hdfier)    ! check whether we forgot to close some resources
  if (obj_count.gt.0 .and. myid.eq.0) then
    write(*,'("There are still ",i3," HDF5 objects open!")') obj_count
  endif ! (obj_count.gt.0)

  call h5close_f( hdfier ) ! close Fortran interface to the HDF5 library;
  FATAL( sphdf5, hdfier.ne.0, error calling h5close_f )
end subroutine finish_outfile

end module sphdf5
