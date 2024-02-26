!> \file
!> \brief Writes all the output information to \c ext.sp.h5.
!>
!> If the output file already exists, it will be deleted and replaced
!> by an empty one, which gets filled in with the updated data.
!> All calls to the HDF5 API are filtered to only happen from MPI rank-0
!> to be able to use the serial HDF5 library.
!> Parallel HDF5 was considered in the past, but abandoned due to very
!> subtle and irreproducible errors.

!> \brief writing the HDF5 output file
!> \ingroup grp_output
module sphdf5

  use inputlist , only : Wsphdf5, Wmacros
  use fileunits , only : ounit
  use cputiming , only : Tsphdf5
  use allglobal , only : myid, cpus, MPI_COMM_SPEC, ext, skip_write
  use constants , only : version
  use hdf5

  implicit none

  logical, parameter             :: hdfDebug = .false.  !< global flag to enable verbal diarrhea commenting HDF5 operations
  integer, parameter             :: internalHdf5Msg = 0 !< 1: print internal HDF5 error messages; 0: only error messages from sphdf5

  integer                        :: hdfier              !< error flag for HDF5 library
  integer                        :: rank                !< rank of data to write using macros
  integer(hid_t)                 :: file_id             !< default file ID used in macros
  integer(hid_t)                 :: space_id            !< default dataspace ID used in macros
  integer(hid_t)                 :: dset_id             !< default dataset ID used in macros
  integer(hsize_t)               :: onedims(1:1)        !< dimension specifier for one-dimensional data used in macros
  integer(hsize_t)               :: twodims(1:2)        !< dimension specifier for two-dimensional data used in macros
  integer(hsize_t)               :: threedims(1:3)      !< dimension specifier for three-dimensional data used in macros
  logical                        :: grp_exists          !< flags used to signal if a group already exists
  logical                        :: var_exists          !< flags used to signal if a variable already exists

  integer(hid_t)                 :: iteration_dset_id   !< Dataset identifier for "iteration"
  integer(hid_t)                 :: dataspace           !< dataspace for extension by 1 iteration object
  integer(hid_t)                 :: memspace            !< memspace for extension by 1 iteration object
  integer(hsize_t), dimension(1) :: old_data_dims       !< current dimensions of "iterations" dataset
  integer(hsize_t), dimension(1) :: data_dims           !< new dimensions for "iterations" dataset
  integer(hsize_t), dimension(1) :: max_dims            !< maximum dimensions for "iterations" dataset
  integer(hid_t)                 :: plist_id            !< Property list identifier used to activate dataset transfer property
  integer(hid_t)                 :: dt_nDcalls_id       !< Memory datatype identifier (for "nDcalls"  dataset in "/iterations")
  integer(hid_t)                 :: dt_Energy_id        !< Memory datatype identifier (for "Energy"   dataset in "/iterations")
  integer(hid_t)                 :: dt_ForceErr_id      !< Memory datatype identifier (for "ForceErr" dataset in "/iterations")
  integer(hid_t)                 :: dt_iRbc_id          !< Memory datatype identifier (for "iRbc"     dataset in "/iterations")
  integer(hid_t)                 :: dt_iZbs_id          !< Memory datatype identifier (for "iZbs"     dataset in "/iterations")
  integer(hid_t)                 :: dt_iRbs_id          !< Memory datatype identifier (for "iRbs"     dataset in "/iterations")
  integer(hid_t)                 :: dt_iZbc_id          !< Memory datatype identifier (for "iZbc"     dataset in "/iterations")

  integer, parameter             :: rankP=3             !< rank of Poincare data
  integer, parameter             :: rankT=2             !< rank of rotational transform data

  integer(hid_t)                 :: grpPoincare         !< group for Poincare data
  integer(HID_T)                 :: dset_id_t           !< Dataset identifier for \f$\theta\f$ coordinate of field line following
  integer(HID_T)                 :: dset_id_s           !< Dataset identifier for \f$s\f$ coordinate of field line following
  integer(HID_T)                 :: dset_id_R           !< Dataset identifier for \f$R\f$ coordinate of field line following
  integer(HID_T)                 :: dset_id_Z           !< Dataset identifier for \f$Z\f$ coordinate of field line following
  integer(HID_T)                 :: dset_id_success     !< Dataset identifier for success flag of trajectories to follow
  integer(HID_T)                 :: filespace_t         !< Dataspace identifier in file for \f$\theta\f$ coordinate of field line following
  integer(HID_T)                 :: filespace_s         !< Dataspace identifier in file for \f$s\f$ coordinate of field line following
  integer(HID_T)                 :: filespace_R         !< Dataspace identifier in file for \f$R\f$ coordinate of field line following
  integer(HID_T)                 :: filespace_Z         !< Dataspace identifier in file for \f$Z\f$ coordinate of field line following
  integer(HID_T)                 :: filespace_success   !< Dataspace identifier in file for success flag of trajectories to follow
  integer(HID_T)                 :: memspace_t          !< Dataspace identifier in memory for \f$\theta\f$ coordinate of field line following
  integer(HID_T)                 :: memspace_s          !< Dataspace identifier in memory for \f$s\f$ coordinate of field line following
  integer(HID_T)                 :: memspace_R          !< Dataspace identifier in memory for \f$R\f$ coordinate of field line following
  integer(HID_T)                 :: memspace_Z          !< Dataspace identifier in memory for \f$Z\f$ coordinate of field line following
  integer(HID_T)                 :: memspace_success    !< Dataspace identifier in memory for success flag of trajectories to follow

  integer(hid_t)                 :: grpTransform        !< group for rotational transform data
  integer(HID_T)                 :: dset_id_diotadxup   !< Dataset identifier for diotadxup (derivative of rotational transform ?)
  integer(HID_T)                 :: dset_id_fiota       !< Dataset identifier for fiota     (              rotational transform ?)
  integer(HID_T)                 :: filespace_diotadxup !< Dataspace identifier in file for diotadxup
  integer(HID_T)                 :: filespace_fiota     !< Dataspace identifier in file for fiota
  integer(HID_T)                 :: memspace_diotadxup  !< Dataspace identifier in memory for diotadxup
  integer(HID_T)                 :: memspace_fiota      !< Dataspace identifier in memory for fiota


  character(LEN=15), parameter :: aname = "description"   !< Attribute name for descriptive info

  integer(HID_T) :: attr_id       !< Attribute identifier
  integer(HID_T) :: aspace_id     !< Attribute Dataspace identifier
  integer(HID_T) :: atype_id      !< Attribute Datatype identifier

  integer, parameter     ::   arank = 1               !< Attribure rank
  integer(HSIZE_T), dimension(arank) :: adims = (/1/) !< Attribute dimension

  integer(SIZE_T) :: attrlen    !< Length of the attribute string
  character(len=:), allocatable ::  attr_data  !< Attribute data

contains

!> \brief Initialize the interface to the HDF5 library and open the output file.
!> \ingroup grp_output
!>
subroutine init_outfile

  LOCALS

  BEGIN( sphdf5 )

 if (myid.eq.0 .and. .not.skip_write) then

  ! initialize Fortran interface to the HDF5 library;
  H5CALL( sphdf5, h5open_f, (hdfier), __FILE__, __LINE__)

  ! (en/dis)able HDF5 internal error messages; sphdf5 has its own error messages coming from the macros
  H5CALL( sphdf5, h5eset_auto_f, (internalHdf5Msg, hdfier), __FILE__, __LINE__)

  ! Create the file
  H5CALL( sphdf5, h5fcreate_f, (trim(ext)//".sp.h5", H5F_ACC_TRUNC_F, file_id, hdfier ), __FILE__, __LINE__ )

  ! write version number
  HWRITERV_LO( file_id, 1, version, (/ version /), __FILE__, __LINE__)
  H5DESCR_CDSET( /version, version of SPEC,        __FILE__, __LINE__)

 endif ! myid.eq.0

end subroutine init_outfile

!> \brief Mirror input variables into output file.
!> \ingroup grp_output
!>
!> The goal of this routine is to have an exact copy of the input file contents
!> that were used to parameterize a given SPEC run.
!> This also serves to check after the run if SPEC correctly understood the text-based input file.
subroutine mirror_input_to_outfile

  use inputlist
  use allglobal , only : Mvol

  LOCALS

  integer(hid_t) :: grpInput
  integer(hid_t) :: grpInputPhysics, grpInputNumerics, grpInputLocal, grpInputGlobal, grpInputDiagnostics

  BEGIN( sphdf5 )

 if (myid.eq.0 .and. .not.skip_write) then

  HDEFGRP( file_id, input, grpInput,                        __FILE__, __LINE__ )
  H5DESCR( grpInput, /input, group for mirrored input data, __FILE__, __LINE__ )

! the following variables constitute the namelist/physicslist/; note that all variables in namelist need to be broadcasted in readin;
! they go into ext.h5/input/physics

  ! the calls used here work as follows:
  ! step 1. HWRITEIV_LO e.g. write(s an) i(nteger) v(ariable) and l(eaves) o(pen) the dataset, so that in
  ! step 2a. an attribute with descr(iptive) information can be attached to the dataset and finally, in
  ! step 2b. the attribute is closed and also we c(lose the) d(ata)set.

  HDEFGRP( grpInput, physics, grpInputPhysics,                                                                                         __FILE__, __LINE__)
  H5DESCR( grpInputPhysics, /input/physics, physics inputs,                                                                            __FILE__, __LINE__)

  HWRITEIV_LO( grpInputPhysics,         1, Igeometry  , (/ Igeometry      /),                                                          __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/Igeometry, geometry identifier,                                                                        __FILE__, __LINE__)
  HWRITEIV_LO( grpInputPhysics,         1, Istellsym  , (/ Istellsym      /),                                                          __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/Istellsym, stellarator symmetry flag,                                                                  __FILE__, __LINE__)
  HWRITEIV_LO( grpInputPhysics,         1, Lfreebound , (/ Lfreebound     /),                                                          __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/Lfreebound, free boundary flag,                                                                        __FILE__, __LINE__)
  HWRITERV_LO( grpInputPhysics,         1, phiedge    , (/ phiedge        /),                                                          __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/phiedge, total enclosed toroidal magnetic flux,                                                        __FILE__, __LINE__)
  HWRITERV_LO( grpInputPhysics,         1, curtor     , (/ curtor         /),                                                          __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/curtor, total enclosed toroidal current,                                                               __FILE__, __LINE__)
  HWRITERV_LO( grpInputPhysics,         1, curpol     , (/ curpol         /),                                                          __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/curpol, total enclosed poloidal current,                                                               __FILE__, __LINE__)
  HWRITERV_LO( grpInputPhysics,         1, gamma      , (/ gamma          /),                                                          __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/gamma, adiabatic index,                                                                                __FILE__, __LINE__)
  HWRITEIV_LO( grpInputPhysics,         1, Nfp        , (/ Nfp            /),                                                          __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/Nfp, number of stellarator field periods,                                                              __FILE__, __LINE__)
  HWRITEIV_LO( grpInputPhysics,         1, Nvol       , (/ Nvol           /),                                                          __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/Nvol, number of volumes,                                                                               __FILE__, __LINE__)
  HWRITEIV_LO( grpInputPhysics,         1, Mpol       , (/ Mpol           /),                                                          __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/Mpol, maximum poloidal mode number,                                                                    __FILE__, __LINE__)
  HWRITEIV_LO( grpInputPhysics,         1, Ntor       , (/ Ntor           /),                                                          __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/Ntor, maximum toroidal mode number,                                                                    __FILE__, __LINE__)
  HWRITEIV_LO( grpInputPhysics,      Mvol, Lrad       ,      Lrad(1:Mvol)   ,                                                          __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/Lrad, degree of radial Chebychev polynomials,                                                          __FILE__, __LINE__)
  HWRITEIV_LO( grpInputPhysics,         1, Lconstraint, (/ Lconstraint    /),                                                          __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/Lconstraint, type of constraint to enforce,                                                            __FILE__, __LINE__)
  HWRITEIV_LO( grpInputPhysics,         1, Lreflect,    (/ Lreflect       /),                                                          __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/Lreflect, whether to reflect the perturbation on both boundaries for slab geometry                     __FILE__, __LINE__)
  HWRITERV_LO( grpInputPhysics,      Mvol, tflux      ,     tflux(1:Mvol)   ,                                                          __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/tflux, toroidal magnetic flux in volumes,                                                              __FILE__, __LINE__)
  HWRITERV_LO( grpInputPhysics,      Mvol, pflux      ,     pflux(1:Mvol)   ,                                                          __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/pflux, poloidal magnetic flux in volumes,                                                              __FILE__, __LINE__)
  HWRITERV_LO( grpInputPhysics,      Nvol, helicity   ,  helicity(1:Nvol)   ,                                                          __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/helicity, helicity profile,                                                                            __FILE__, __LINE__)
  HWRITERV_LO( grpInputPhysics,         1, pscale     , (/ pscale         /),                                                          __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/pscale, scaling factor for pressure,                                                                   __FILE__, __LINE__)
  HWRITERV_LO( grpInputPhysics,      Nvol, pressure   ,  pressure(1:Nvol)   ,                                                          __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/pressure, pressure profile,                                                                            __FILE__, __LINE__)
  HWRITEIV_LO( grpInputPhysics,         1, Ladiabatic , (/ Ladiabatic     /),                                                          __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/Ladiabatic, adiabatic flag,                                                                            __FILE__, __LINE__)
  HWRITERV_LO( grpInputPhysics,      Mvol, adiabatic  , adiabatic(1:Nvol)   ,                                                          __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/adiabatic, adiabatic profile (?),                                                                      __FILE__, __LINE__)
  HWRITERV_LO( grpInputPhysics,  (1+Nvol), mu         ,        mu(1:Mvol)   ,                                                          __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/mu, Beltrami parameter{,} parallel current profile,                                                    __FILE__, __LINE__)
  HWRITERV_LO( grpInputPhysics,  (1+Nvol), Ivolume    ,   Ivolume(1:Mvol)   ,                                                          __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/Ivolume, Volume current{,} externally driven, parallel current profile,                                __FILE__, __LINE__)
  HWRITERV_LO( grpInputPhysics,  (Mvol  ), Isurf      ,   Isurf(1:Mvol  )   ,                                                          __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/mu, Surface current{,} currents that are not volume currents (pressure driven, shielding currents) ,   __FILE__, __LINE__)
  HWRITEIV_LO( grpInputPhysics,  (1+Mvol), pl         ,        pl(0:Nvol)   ,                                                          __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/pl, pl ?,                                                                                              __FILE__, __LINE__)
  HWRITEIV_LO( grpInputPhysics,  (1+Mvol), ql         ,        ql(0:Nvol)   ,                                                          __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/ql, ql ?,                                                                                              __FILE__, __LINE__)
  HWRITEIV_LO( grpInputPhysics,  (1+Mvol), pr         ,        pr(0:Nvol)   ,                                                          __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/pr, pr ?,                                                                                              __FILE__, __LINE__)
  HWRITEIV_LO( grpInputPhysics,  (1+Mvol), qr         ,        qr(0:Nvol)   ,                                                          __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/qr, qr ?,                                                                                              __FILE__, __LINE__)
  HWRITERV_LO( grpInputPhysics,  (1+Nvol), iota       ,      iota(0:Nvol)   ,                                                          __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/iota, rotational transform profile on inside of ideal interfaces,                                      __FILE__, __LINE__)
  HWRITEIV_LO( grpInputPhysics,  (1+Mvol), lp         ,        lp(0:Nvol)   ,                                                          __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/lp, lp ?,                                                                                              __FILE__, __LINE__)
  HWRITEIV_LO( grpInputPhysics,  (1+Mvol), lq         ,        lq(0:Nvol)   ,                                                          __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/lq, lq ?,                                                                                              __FILE__, __LINE__)
  HWRITEIV_LO( grpInputPhysics,  (1+Mvol), rp         ,        rp(0:Nvol)   ,                                                          __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/rp, rp ?,                                                                                              __FILE__, __LINE__)
  HWRITEIV_LO( grpInputPhysics,  (1+Mvol), rq         ,        rq(0:Nvol)   ,                                                          __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/rq, rq ?,                                                                                              __FILE__, __LINE__)
  HWRITERV_LO( grpInputPhysics,  (1+Nvol), oita       ,      oita(0:Nvol)   ,                                                          __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/oita, rotational transform profile on outside of ideal interfaces,                                     __FILE__, __LINE__)
  HWRITERV_LO( grpInputPhysics,         1, rtor  , (/ rtor      /),                                                                    __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/rpol, for aspect ratio in slab,                                                                        __FILE__, __LINE__)
  HWRITERV_LO( grpInputPhysics,         1, rpol  , (/ rpol      /),                                                                    __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/rpol, for aspect ratio in slab,                                                                        __FILE__, __LINE__)

  HWRITERV_LO( grpInputPhysics,  (1+Ntor), Rac        ,       Rac(0:Ntor)   ,                                                          __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/Rac,     stellarator symmetric coordinate axis R cosine Fourier coefficients,                          __FILE__, __LINE__)
  HWRITERV_LO( grpInputPhysics,  (1+Ntor), Zas        ,       Zas(0:Ntor)   ,                                                          __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/Zas,     stellarator symmetric coordinate axis Z   sine Fourier coefficients,                          __FILE__, __LINE__)
  HWRITERV_LO( grpInputPhysics,  (1+Ntor), Ras        ,       Ras(0:Ntor)   ,                                                          __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/Ras, non-stellarator symmetric coordinate axis R   sine Fourier coefficients,                          __FILE__, __LINE__)
  HWRITERV_LO( grpInputPhysics,  (1+Ntor), Zac        ,       Zac(0:Ntor)   ,                                                          __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/Zac, non-stellarator symmetric coordinate axis Z cosine Fourier coefficients,                          __FILE__, __LINE__)

  HWRITERA_LO( grpInputPhysics, (2*Ntor+1), (2*Mpol+1),  Rbc, Rbc(-Ntor:Ntor,-Mpol:Mpol),                                              __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/Rbc,     stellarator symmetric boundary R cosine Fourier coefficients,                                 __FILE__, __LINE__)
  HWRITERA_LO( grpInputPhysics, (2*Ntor+1), (2*Mpol+1),  Zbs, Zbs(-Ntor:Ntor,-Mpol:Mpol),                                              __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/Zbs,     stellarator symmetric boundary Z   sine Fourier coefficients,                                 __FILE__, __LINE__)
  HWRITERA_LO( grpInputPhysics, (2*Ntor+1), (2*Mpol+1),  Rbs, Rbs(-Ntor:Ntor,-Mpol:Mpol),                                              __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/Rbs, non-stellarator symmetric boundary R   sine Fourier coefficients,                                 __FILE__, __LINE__)
  HWRITERA_LO( grpInputPhysics, (2*Ntor+1), (2*Mpol+1),  Zbc, Zbc(-Ntor:Ntor,-Mpol:Mpol),                                              __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/Zbc, non-stellarator symmetric boundary Z cosine Fourier coefficients,                                 __FILE__, __LINE__)

  HWRITERA_LO( grpInputPhysics, (2*Ntor+1), (2*Mpol+1),  Rwc, Rwc(-Ntor:Ntor,-Mpol:Mpol),                                              __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/Rwc,     stellarator symmetric boundary R cosine Fourier coefficients of wall,                         __FILE__, __LINE__)
  HWRITERA_LO( grpInputPhysics, (2*Ntor+1), (2*Mpol+1),  Zws, Zws(-Ntor:Ntor,-Mpol:Mpol),                                              __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/Zws,     stellarator symmetric boundary Z   sine Fourier coefficients of wall,                         __FILE__, __LINE__)
  HWRITERA_LO( grpInputPhysics, (2*Ntor+1), (2*Mpol+1),  Rws, Rws(-Ntor:Ntor,-Mpol:Mpol),                                              __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/Rws, non-stellarator symmetric boundary R   sine Fourier coefficients of wall,                         __FILE__, __LINE__)
  HWRITERA_LO( grpInputPhysics, (2*Ntor+1), (2*Mpol+1),  Zwc, Zwc(-Ntor:Ntor,-Mpol:Mpol),                                              __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/Zwc, non-stellarator symmetric boundary Z cosine Fourier coefficients of wall,                         __FILE__, __LINE__)

  HWRITERA_LO( grpInputPhysics, (2*Ntor+1), (2*Mpol+1),  Vns, Vns(-Ntor:Ntor,-Mpol:Mpol),                                              __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/Vns,     stellarator symmetric normal field   sine Fourier coefficients at boundary; vacuum component, __FILE__, __LINE__)
  HWRITERA_LO( grpInputPhysics, (2*Ntor+1), (2*Mpol+1),  Bns, Bns(-Ntor:Ntor,-Mpol:Mpol),                                              __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/Bns,     stellarator symmetric normal field   sine Fourier coefficients at boundary; plasma component, __FILE__, __LINE__)
  HWRITERA_LO( grpInputPhysics, (2*Ntor+1), (2*Mpol+1),  Vnc, Vnc(-Ntor:Ntor,-Mpol:Mpol),                                              __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/Vnc, non-stellarator symmetric normal field cosine Fourier coefficients at boundary; vacuum component, __FILE__, __LINE__)
  HWRITERA_LO( grpInputPhysics, (2*Ntor+1), (2*Mpol+1),  Bnc, Bnc(-Ntor:Ntor,-Mpol:Mpol),                                              __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/Bnc, non-stellarator symmetric normal field cosine Fourier coefficients at boundary; plasma component, __FILE__, __LINE__)

  HWRITERV_LO( grpInputPhysics,           1,   mupftol, (/ mupftol /),                                                                 __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/mupftol, mupftol   ,                                                                                   __FILE__, __LINE__)
  HWRITEIV_LO( grpInputPhysics,           1,   mupfits, (/ mupfits /),                                                                 __FILE__, __LINE__)
  H5DESCR_CDSET( /input/physics/mupfits, mupfits   ,                                                                                   __FILE__, __LINE__)

  HCLOSEGRP( grpInputPhysics , __FILE__, __LINE__)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! the following variables constitute the namelist/numericlist/; note that all variables in namelist need to be broadcasted in readin;
! they go into ext.h5/input/numerics

  HDEFGRP( grpInput, numerics, grpInputNumerics, __FILE__, __LINE__)

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
  HWRITEIV( grpInputNumerics,          1, Lrzaxis            , (/ Lrzaxis     /))
  HWRITEIV( grpInputNumerics,          1, Ntoraxis           , (/ Ntoraxis    /))

  HCLOSEGRP( grpInputNumerics, __FILE__, __LINE__)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! the following variables constitute the namelist/locallist/; note that all variables in namelist need to be broadcasted in readin;
! they go into ext.h5/input/local

  HDEFGRP( grpInput, local, grpInputLocal )

  HWRITEIV( grpInputLocal,             1, LBeltrami          , (/ LBeltrami   /))
  HWRITEIV( grpInputLocal,             1, Linitgues          , (/ Linitgues   /))
  HWRITEIV( grpInputLocal,             1, Lposdef            , (/ Lposdef     /)) ! redundant;
  HWRITERV( grpInputLocal,             1, maxrndgues         , (/ maxrndgues  /))
  HWRITEIV( grpInputLocal,             1, Lmatsolver         , (/ Lmatsolver  /))
  HWRITEIV( grpInputLocal,             1, LGMRESprec         , (/ LGMRESprec  /))
  HWRITERV( grpInputLocal,             1, epsGMRES           , (/ epsGMRES    /))
  HWRITERV( grpInputLocal,             1, epsILU             , (/ epsILU      /))

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
  HWRITEIV( grpInputDiagnostics,       1,   Ppts             , (/  Ppts          /))
  HWRITEIV( grpInputDiagnostics,    Mvol,  nPtrj             ,    nPtrj(1:Mvol)    )
  HWRITELV( grpInputDiagnostics,       1,  LHevalues         , (/ LHevalues      /))
  HWRITELV( grpInputDiagnostics,       1,  LHevectors        , (/ LHevectors     /))
  HWRITELV( grpInputDiagnostics,       1,  LHmatrix          , (/ LHmatrix       /))
  HWRITELV( grpInputDiagnostics,       1,  Ltransform        , (/ Ltransform     /))
  HWRITEIV( grpInputDiagnostics,       1,  Lperturbed        , (/ Lperturbed     /))
  HWRITEIV( grpInputDiagnostics,       1,  dpp               , (/ dpp            /))
  HWRITEIV( grpInputDiagnostics,       1,  dqq               , (/ dqq            /))
  HWRITEIV( grpInputDiagnostics,       1,  Lcheck            , (/ Lcheck         /))
  HWRITELV( grpInputDiagnostics,       1,  Ltiming           , (/ Ltiming        /))
  HWRITEIV( grpInputDiagnostics,       1,  Lerrortype        , (/ Lerrortype     /))
  HWRITEIV( grpInputDiagnostics,       1,  Ngrid             , (/ Ngrid          /))
  HWRITERV( grpInputDiagnostics,       1,  fudge             , (/ fudge          /))          ! redundant;
  HWRITERV( grpInputDiagnostics,       1,  scaling           , (/ scaling        /))          ! redundant;

  HCLOSEGRP( grpInputDiagnostics )

  HCLOSEGRP( grpInput )

 endif ! myid.eq.0

end subroutine mirror_input_to_outfile

!> \brief Prepare convergence evolution output.
!> \ingroup grp_output
!>
!> <ul>
!> <li> The group \c iterations is created in the output file.
!>      This group contains the interface geometry at each iteration, which is useful for constructing movies illustrating the convergence.
!>      The data structure in use is an unlimited array of the following compound datatype:
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
subroutine init_convergence_output

  use allglobal, only : mn, Mvol

  LOCALS

  integer(hid_t)                    :: iteration_dspace_id  !< dataspace for "iteration"
  integer(hid_t)                    :: iteration_dtype_id   !< Compound datatype for "iteration"
  integer(hid_t)                    :: iRZbscArray_id       !< Memory datatype identifier
  integer(size_t)                   :: iteration_dtype_size !< Size of the "iteration" datatype
  integer(size_t)                   :: type_size_i          !< Size of the integer datatype
  integer(size_t)                   :: type_size_d          !< Size of the double precision datatype
  integer(size_t)                   :: offset               !< Member's offset
  integer(hid_t)                    :: crp_list             !< Dataset creation property identifier
  integer, parameter                :: rank = 1             !< logging rank: convergence logging is one-dimensional
  integer(hsize_t), dimension(rank) :: maxdims              !< convergence logging maximum dimensions => will be unlimited
  integer(hsize_t), dimension(rank) :: dims = (/ 0 /)       !< current convergence logging dimensions
  integer(hsize_t), dimension(rank) :: dimsc = (/ 1 /)      !< chunking length ???
  integer(size_t)                   :: irbc_size_template   !< size ofiRbc array in iterations logging
  integer(size_t)                   :: irbc_size            !< size ofiRbc array in iterations logging

  BEGIN( sphdf5 )

 if (myid.eq.0 .and. .not.skip_write) then

  ! Set dataset transfer property to preserve partially initialized fields
  ! during write/read to/from dataset with compound datatype.
  H5CALL( sphdf5, h5pcreate_f, (H5P_DATASET_XFER_F, plist_id, hdfier))

  H5CALL( sphdf5, h5pset_preserve_f, (plist_id, .TRUE., hdfier) )

  maxdims = (/ H5S_UNLIMITED_F /)                                                                      ! unlimited array size: "converge" until you get bored
  H5CALL( sphdf5, h5screate_simple_f, (rank, dims, iteration_dspace_id, hdfier, maxdims) )             ! Create the dataspace with zero initial size and allow it to grow

  H5CALL( sphdf5, h5pcreate_f, (H5P_DATASET_CREATE_F, crp_list, hdfier) )                              ! dataset creation property list with chunking

  H5CALL( sphdf5, h5pset_chunk_f, (crp_list, rank, dimsc, hdfier) )

  ! declare "iteration" compound datatype
  ! declare array parts
  H5CALL( sphdf5, h5tarray_create_f, (H5T_NATIVE_DOUBLE, 2, int((/mn, Mvol+1/),hsize_t), iRZbscArray_id, hdfier) ) ! create array datatypes for i{R,Z}b{c,s}
  H5CALL( sphdf5, h5tget_size_f, (iRZbscArray_id, irbc_size, hdfier) )
  H5CALL( sphdf5, h5tget_size_f, (H5T_NATIVE_INTEGER, type_size_i, hdfier) )                           ! size of an integer field
  H5CALL( sphdf5, h5tget_size_f, (H5T_NATIVE_DOUBLE,  type_size_d, hdfier) )                           ! size of a   double field
  iteration_dtype_size = 2*type_size_i + 2*type_size_d + 4*irbc_size                                   ! wflag, nDcalls, Energy, ForceErr, i{R,Z}b{c,s}

  H5CALL( sphdf5, h5tcreate_f, (H5T_COMPOUND_F, iteration_dtype_size, iteration_dtype_id, hdfier) )    ! create compound datatype

  offset = 0                                                                                           ! offset for first field starts at 0
  H5CALL( sphdf5, h5tinsert_f, (iteration_dtype_id, "nDcalls", offset, H5T_NATIVE_INTEGER, hdfier) )   ! insert "nDcalls" field in datatype
  offset = offset + type_size_i                                                                        ! increment offset by size of field
  H5CALL( sphdf5, h5tinsert_f, (iteration_dtype_id, "Energy", offset, H5T_NATIVE_DOUBLE, hdfier) )     ! insert "Energy" field in datatype
  offset = offset + type_size_d                                                                        ! increment offset by size of field
  H5CALL( sphdf5, h5tinsert_f, (iteration_dtype_id, "ForceErr", offset, H5T_NATIVE_DOUBLE, hdfier) )   ! insert "ForceErr" field in datatype
  offset = offset + type_size_d                                                                        ! increment offset by size of field
  H5CALL( sphdf5, h5tinsert_f, (iteration_dtype_id, "iRbc", offset, iRZbscArray_id, hdfier) )          ! insert "iRbc" field in datatype
  offset = offset + irbc_size                                                                          ! increment offset by size of field
  H5CALL( sphdf5, h5tinsert_f, (iteration_dtype_id, "iZbs", offset, iRZbscArray_id, hdfier) )          ! insert "iZbs" field in datatype
  offset = offset + irbc_size                                                                          ! increment offset by size of field
  H5CALL( sphdf5, h5tinsert_f, (iteration_dtype_id, "iRbs", offset, iRZbscArray_id, hdfier) )          ! insert "iRbs" field in datatype
  offset = offset + irbc_size                                                                          ! increment offset by size of field
  H5CALL( sphdf5, h5tinsert_f, (iteration_dtype_id, "iZbc", offset, iRZbscArray_id, hdfier) )          ! insert "iZbc" field in datatype
  offset = offset + irbc_size                                                                          ! increment offset by size of field

  H5CALL( sphdf5, h5dcreate_f, (file_id, "iterations", iteration_dtype_id, iteration_dspace_id, &      ! create dataset with compound type
                  & iteration_dset_id, hdfier, crp_list) )

  H5CALL( sphdf5, h5sclose_f, (iteration_dspace_id, hdfier) )                                          ! Terminate access to the data space (does not show up in obj_count below)
                                                                                                       ! --> only needed for creation of dataset

  ! Create memory types. We have to create a compound datatype
  ! for each member we want to write.
  offset = 0
  H5CALL( sphdf5, h5tcreate_f, (H5T_COMPOUND_F, type_size_i, dt_nDcalls_id,  hdfier) )
  H5CALL( sphdf5, h5tcreate_f, (H5T_COMPOUND_F, type_size_d, dt_Energy_id,   hdfier) )
  H5CALL( sphdf5, h5tcreate_f, (H5T_COMPOUND_F, type_size_d, dt_ForceErr_id, hdfier) )
  H5CALL( sphdf5, h5tcreate_f, (H5T_COMPOUND_F, irbc_size,   dt_iRbc_id,     hdfier) )
  H5CALL( sphdf5, h5tcreate_f, (H5T_COMPOUND_F, irbc_size,   dt_iZbs_id,     hdfier) )
  H5CALL( sphdf5, h5tcreate_f, (H5T_COMPOUND_F, irbc_size,   dt_iRbs_id,     hdfier) )
  H5CALL( sphdf5, h5tcreate_f, (H5T_COMPOUND_F, irbc_size,   dt_iZbc_id,     hdfier) )

  H5CALL( sphdf5, h5tinsert_f, (dt_nDcalls_id,   "nDcalls", offset, H5T_NATIVE_INTEGER, hdfier) )
  H5CALL( sphdf5, h5tinsert_f, (dt_Energy_id,     "Energy", offset, H5T_NATIVE_DOUBLE,  hdfier) )
  H5CALL( sphdf5, h5tinsert_f, (dt_ForceErr_id, "ForceErr", offset, H5T_NATIVE_DOUBLE,  hdfier) )
  H5CALL( sphdf5, h5tinsert_f, (dt_iRbc_id,         "iRbc", offset, iRZbscArray_id,     hdfier) )
  H5CALL( sphdf5, h5tinsert_f, (dt_iZbs_id,         "iZbs", offset, iRZbscArray_id,     hdfier) )
  H5CALL( sphdf5, h5tinsert_f, (dt_iRbs_id,         "iRbs", offset, iRZbscArray_id,     hdfier) )
  H5CALL( sphdf5, h5tinsert_f, (dt_iZbc_id,         "iZbc", offset, iRZbscArray_id,     hdfier) )

  ! create memspace with size of compound object to append
  dims(1) = 1 ! only append one iteration at a time
  H5CALL( sphdf5, h5screate_simple_f, (rank, dims, memspace, hdfier) )

  H5CALL( sphdf5, h5pclose_f, (crp_list, hdfier) )
  H5CALL( sphdf5, h5tclose_f, (iteration_dtype_id, hdfier) )                                       ! Terminate access to the datatype
  H5CALL( sphdf5, h5tclose_f, (iRZbscArray_id, hdfier) )                                           ! Terminate access to the datatype

 endif ! myid.eq.0

end subroutine init_convergence_output


!> \brief Write convergence output (evolution of interface geometry, force, etc).
!> \ingroup grp_output
!>
subroutine write_convergence_output( nDcalls, ForceErr )

  use allglobal, only : myid, mn, Mvol, Energy, iRbc, iZbs, iRbs, iZbc

  LOCALS
  INTEGER, intent(in)  :: nDcalls
  REAL   , intent(in)  :: ForceErr

  BEGIN(sphdf5)

 if (myid.eq.0 .and. .not.skip_write) then

  ! append updated values to "iterations" dataset

  ! open dataspace to get current state of dataset
  H5CALL( sphdf5, h5dget_space_f, (iteration_dset_id, dataspace, hdfier), __FILE__, __LINE__)

  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, old_data_dims, max_dims, hdfier)
  FATAL( sphdf5, hdfier.ne.1, rank of convergence dataspace is not 1 )

  ! blow up dataset to new size
  data_dims = old_data_dims+1
  H5CALL( sphdf5, h5dset_extent_f, (iteration_dset_id, data_dims, hdfier), __FILE__, __LINE__)

  ! get dataspace slab corresponding to region which the iterations dataset was extended by
  H5CALL( sphdf5, h5dget_space_f, (iteration_dset_id, dataspace, hdfier), __FILE__, __LINE__)                                             ! re-select dataspace to update size info in HDF5 lib
  H5CALL( sphdf5, h5sselect_hyperslab_f, (dataspace, H5S_SELECT_SET_F, old_data_dims, (/ INT(1, HSIZE_T) /), hdfier), __FILE__, __LINE__) ! newly appended slab is at old size and 1 long

  ! write next iteration object
  H5CALL( sphdf5, h5dwrite_f, (iteration_dset_id, dt_nDcalls_id, nDcalls, INT((/1/), HSIZE_T), hdfier, &
    & mem_space_id=memspace, file_space_id=dataspace, xfer_prp=plist_id), __FILE__, __LINE__)
  H5CALL( sphdf5, h5dwrite_f, (iteration_dset_id, dt_Energy_id, Energy, INT((/1/), HSIZE_T), hdfier, &
    & mem_space_id=memspace, file_space_id=dataspace, xfer_prp=plist_id), __FILE__, __LINE__)
  H5CALL( sphdf5, h5dwrite_f, (iteration_dset_id, dt_ForceErr_id, ForceErr, INT((/1/), HSIZE_T), hdfier, &
    & mem_space_id=memspace, file_space_id=dataspace, xfer_prp=plist_id), __FILE__, __LINE__)
  H5CALL( sphdf5, h5dwrite_f, (iteration_dset_id, dt_iRbc_id, iRbc, INT((/mn,Mvol+1/), HSIZE_T), hdfier, &
    & mem_space_id=memspace, file_space_id=dataspace, xfer_prp=plist_id), __FILE__, __LINE__)
  H5CALL( sphdf5, h5dwrite_f, (iteration_dset_id, dt_iZbs_id, iZbs, INT((/mn,Mvol+1/), HSIZE_T), hdfier, &
    & mem_space_id=memspace, file_space_id=dataspace, xfer_prp=plist_id), __FILE__, __LINE__)
  H5CALL( sphdf5, h5dwrite_f, (iteration_dset_id, dt_iRbs_id, iRbs, INT((/mn,Mvol+1/), HSIZE_T), hdfier, &
    & mem_space_id=memspace, file_space_id=dataspace, xfer_prp=plist_id), __FILE__, __LINE__)
  H5CALL( sphdf5, h5dwrite_f, (iteration_dset_id, dt_iZbc_id, iZbc, INT((/mn,Mvol+1/), HSIZE_T), hdfier, &
    & mem_space_id=memspace, file_space_id=dataspace, xfer_prp=plist_id), __FILE__, __LINE__)

  ! dataspace to appended object should be closed now
  ! MAYBE we otherwise keep all the iterations in memory?
  H5CALL( sphdf5, h5sclose_f, (dataspace, hdfier), __FILE__, __LINE__)

 endif ! myid.eq.0

end subroutine write_convergence_output

!> \brief Write the magnetic field on a grid.
!> \ingroup grp_output
!>
!> The magnetic field is evaluated on a regular grid in \f$(s, \theta, \zeta)\f$
!> and the corresponding cylindrical coordinates \f$(R,Z)\f$
!> as well as the cylindrical components of the magnetic field \f$(B^R, B^\varphi, B^Z)\f$
!> are written out.
subroutine write_grid

  use constants
  use allglobal, only : myid, ijreal, ijimag, jireal, &
  &                     Nt, Nz, Ntz, Mvol, pi2nfp, ivol, mn, Node, gBzeta, &
  &                     Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, &
  &                     Rij, Zij, sg
  use inputlist, only : Lrad, Igeometry, Nvol, Ngrid, rtor, rpol
  use cputiming, only : Tsphdf5

  LOCALS
  integer(hid_t) :: grpGrid
  integer :: sumLrad, alongLrad, Ngrid_local, Ngrid_sum
  INTEGER              :: vvol, ii, jj, kk, jk, Lcurvature
  REAL                 :: lss, teta, zeta, st(1:Node), Bst(1:Node)
  REAL   , allocatable :: Rij_grid(:,:), Zij_grid(:,:), sg_grid(:,:), ijreal_grid(:,:), ijimag_grid(:,:), jireal_grid(:,:)

  BEGIN(sphdf5)

 if (myid.eq.0 .and. .not.skip_write) then

  ijreal(1:Ntz) = zero ; ijimag(1:Ntz) = zero ; jireal(1:Ntz) = zero

  HDEFGRP( file_id, grid, grpGrid )

  ! Igeometry already is in input, Mvol already is in output
  HWRITEIV( grpGrid,           1, Nt               , (/ Nt            /))
  HWRITEIV( grpGrid,           1, Nz               , (/ Nz            /))
  HWRITEIV( grpGrid,           1, Ntz              , (/ Ntz           /))
  HWRITERV( grpGrid,           1, pi2nfp           , (/ pi2nfp        /))

  ! combine all radial parts into one dimension as Lrad values can be different for different volumes
  if (Ngrid .lt. 0) then
    sumLrad = sum(Lrad(1:Mvol)+1)
  else
    sumLrad = (Ngrid + 1) * Mvol
  endif

  SALLOCATE(    Rij_grid, (1:sumLrad, 1:Ntz), zero )
  SALLOCATE(    Zij_grid, (1:sumLrad, 1:Ntz), zero )
  SALLOCATE(     sg_grid, (1:sumLrad, 1:Ntz), zero )
  SALLOCATE( ijreal_grid, (1:sumLrad, 1:Ntz), zero )
  SALLOCATE( ijimag_grid, (1:sumLrad, 1:Ntz), zero )
  SALLOCATE( jireal_grid, (1:sumLrad, 1:Ntz), zero )

  Ngrid_sum = 0

  do vvol = 1, Mvol ; ivol = vvol
   LREGION(vvol) ! sets Lcoordinatesingularity and Lplasmaregion ;

   if (Ngrid .lt. 0) then
    Ngrid_local = Lrad(vvol)  ! default
   else
    Ngrid_local = Ngrid
   endif
   if (Ngrid_local .eq. 0) cycle               ! nothing to output

   do ii = 0, Ngrid_local ! sub-grid;
    lss = ii * two / Ngrid_local - one
    if( Lcoordinatesingularity .and. ii.eq.0 ) then ; Lcurvature = 0 ! Jacobian is not defined;
    else                                            ; Lcurvature = 1 ! compute Jacobian       ;
    endif

    WCALL( sphdf5, coords, ( vvol, lss, Lcurvature, Ntz, mn ) ) ! only Rij(0,:) and Zij(0,:) are required; Rmn & Zmn are available;

    alongLrad = Ngrid_sum+ii+1

    Rij_grid(alongLrad,1:Ntz) = Rij(1:Ntz,0,0)
    Zij_grid(alongLrad,1:Ntz) = Zij(1:Ntz,0,0)
    sg_grid (alongLrad,1:Ntz) =  sg(1:Ntz,0)

    if( Lcurvature.eq.1 ) then

     select case (Igeometry)

     case (3)
      do kk = 0, Nz-1 ; zeta = kk * pi2nfp / Nz
        do jj = 0, Nt-1 ; teta = jj * pi2    / Nt ; jk = 1 + jj + kk*Nt ; st(1:2) = (/ lss, teta /)
        WCALL( sphdf5, bfield, ( zeta, st(1:Node), Bst(1:Node) ) )
        ijreal(jk) = ( Rij(jk,1,0) * Bst(1) + Rij(jk,2,0) * Bst(2) + Rij(jk,3,0) * one ) * gBzeta / sg(jk,0) ! BR;
        ijimag(jk) = (                                                             one ) * gBzeta / sg(jk,0) ! Bp;
        jireal(jk) = ( Zij(jk,1,0) * Bst(1) + Zij(jk,2,0) * Bst(2) + Zij(jk,3,0) * one ) * gBzeta / sg(jk,0) ! BZ;
        enddo
      enddo

     case (1)
      do kk = 0, Nz-1 ; zeta = kk * pi2nfp / Nz
        do jj = 0, Nt-1 ; teta = jj * pi2    / Nt ; jk = 1 + jj + kk*Nt ; st(1:2) = (/ lss, teta /)
        WCALL( sphdf5, bfield, ( zeta, st(1:Node), Bst(1:Node) ) )
        ijreal(jk) = ( Rij(jk,1,0) * Bst(1) + Rij(jk,2,0) * Bst(2) + Rij(jk,3,0) * one ) * gBzeta / sg(jk,0) ! BR;
        ijimag(jk) = (                                                            rpol ) * gBzeta / sg(jk,0) ! Bzeta;
        jireal(jk) = (                      +        rtor * Bst(2)                     ) * gBzeta / sg(jk,0) ! Btheta;
        enddo
      enddo

     case (2)
      do kk = 0, Nz-1 ; zeta = kk * pi2nfp / Nz
        do jj = 0, Nt-1 ; teta = jj * pi2    / Nt ; jk = 1 + jj + kk*Nt ; st(1:2) = (/ lss, teta /)
        WCALL( sphdf5, bfield, ( zeta, st(1:Node), Bst(1:Node) ) )
        ijreal(jk) = ( Rij(jk,1,0) * Bst(1) + Rij(jk,2,0) * Bst(2) + Rij(jk,3,0) * one ) * gBzeta / sg(jk,0) ! BR;
        ijimag(jk) = (                                                             one ) * gBzeta / sg(jk,0) ! Bp;
        jireal(jk) = (                                      Bst(2)                     ) * gBzeta / sg(jk,0) ! BZ;
        enddo
      enddo

     end select !Igeometry
    endif ! end of if( Lcurvature.eq.1 ) ;

   ijreal_grid(alongLrad,1:Ntz) = ijreal(1:Ntz)
   ijimag_grid(alongLrad,1:Ntz) = ijimag(1:Ntz)
   jireal_grid(alongLrad,1:Ntz) = jireal(1:Ntz)

   enddo ! end of do ii;

   Ngrid_sum = Ngrid_sum + Ngrid_local + 1 ! offset for storing data

  enddo ! end of do vvol;

  HWRITERA( grpGrid, sumLrad, Ntz, Rij,    Rij_grid )
  HWRITERA( grpGrid, sumLrad, Ntz, Zij,    Zij_grid )
  HWRITERA( grpGrid, sumLrad, Ntz,  sg,     sg_grid )
  HWRITERA( grpGrid, sumLrad, Ntz,  BR, ijreal_grid )
  HWRITERA( grpGrid, sumLrad, Ntz,  Bp, ijimag_grid )
  HWRITERA( grpGrid, sumLrad, Ntz,  BZ, jireal_grid )

  DALLOCATE(    Rij_grid )
  DALLOCATE(    Zij_grid )
  DALLOCATE(     sg_grid )
  DALLOCATE( ijreal_grid )
  DALLOCATE( ijimag_grid )
  DALLOCATE( jireal_grid )

  HCLOSEGRP( grpGrid )

  RETURN(sphdf5)

 endif ! myid.eq.0

end subroutine write_grid

!> \brief Initialize field line tracing output group and create array datasets.
!> \ingroup grp_output
!>
!> The field-line tracing diagnostic is parallelized over volumes,
!> where all threads/ranks produce individual output.
!> This is gathered in the output file, stacked over the radial dimension.
!> The \c success flag signals if the integrator was successful in following
!> the fieldline for the derired number of toroidal periods.
!>
!> @param[in] numTrajTotal total number of Poincare trajectories
subroutine init_flt_output( numTrajTotal )

  use allglobal, only : Nz, Mvol, lmns
  use inputlist, only : nPpts

  LOCALS
  integer, intent(in)               :: numTrajTotal ! total number of trajectories
  integer(HSIZE_T), dimension(rankP) :: dims_traj   ! Dataset dimensions.
  integer(HSIZE_T), dimension(rankP) :: length      ! Dataset dimensions.

  BEGIN( sphdf5 )

 if (myid.eq.0 .and. .not.skip_write) then

  ! create Poincare group in HDF5 file
  HDEFGRP( file_id, poincare, grpPoincare )

  dims_traj = (/ Nz, nPpts, numTrajTotal /) ! dimensions for whole Poincare dataset
  length    = (/ Nz, nPpts,            1 /) ! which is written in these slice lengths

  ! Create the data space for the  dataset.
  H5CALL( sphdf5, h5screate_simple_f, (rankP, dims_traj, filespace_t, hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5screate_simple_f, (rankP, dims_traj, filespace_s, hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5screate_simple_f, (rankP, dims_traj, filespace_R, hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5screate_simple_f, (rankP, dims_traj, filespace_Z, hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5screate_simple_f, (1, int((/ numTrajTotal /),HSIZE_T), filespace_success, hdfier), __FILE__, __LINE__ )

  ! Create the dataset with default properties.
  H5CALL( sphdf5, h5dcreate_f, (grpPoincare, "t", H5T_NATIVE_DOUBLE, filespace_t, dset_id_t, hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5dcreate_f, (grpPoincare, "s", H5T_NATIVE_DOUBLE, filespace_s, dset_id_s, hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5dcreate_f, (grpPoincare, "R", H5T_NATIVE_DOUBLE, filespace_R, dset_id_R, hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5dcreate_f, (grpPoincare, "Z", H5T_NATIVE_DOUBLE, filespace_Z, dset_id_Z, hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5dcreate_f, (grpPoincare, "success", H5T_NATIVE_INTEGER, filespace_success, dset_id_success, hdfier), __FILE__, __LINE__ )

  ! filespaces can be closed as soon as datasets are created
  H5CALL( sphdf5, h5sclose_f, (filespace_t, hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5sclose_f, (filespace_s, hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5sclose_f, (filespace_R, hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5sclose_f, (filespace_Z, hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5sclose_f, (filespace_success, hdfier), __FILE__, __LINE__ )

  ! Select hyperslab in the file.
  H5CALL( sphdf5, h5dget_space_f, (dset_id_t, filespace_t, hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5dget_space_f, (dset_id_s, filespace_s, hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5dget_space_f, (dset_id_R, filespace_R, hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5dget_space_f, (dset_id_Z, filespace_Z, hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5dget_space_f, (dset_id_success, filespace_success, hdfier), __FILE__, __LINE__ )

  ! Each process defines dataset in memory and writes it to the hyperslab in the file.
  H5CALL( sphdf5, h5screate_simple_f, (rankP, length, memspace_t, hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5screate_simple_f, (rankP, length, memspace_s, hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5screate_simple_f, (rankP, length, memspace_R, hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5screate_simple_f, (rankP, length, memspace_Z, hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5screate_simple_f, (1, int((/ 1 /),HSIZE_T), memspace_success, hdfier), __FILE__, __LINE__ )

  ! create rotational transform group in HDF5 file
  HDEFGRP( file_id, transform, grpTransform )

  ! Create the data space for the  dataset.
  H5CALL( sphdf5, h5screate_simple_f, (rankT, int((/           2,Mvol/),HSIZE_T), filespace_diotadxup, hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5screate_simple_f, (rankT, int((/numTrajTotal,   2/),HSIZE_T), filespace_fiota    , hdfier), __FILE__, __LINE__ )

  ! Create the dataset with default properties.
  H5CALL( sphdf5, h5dcreate_f, (grpTransform, "diotadxup", H5T_NATIVE_DOUBLE, filespace_diotadxup, dset_id_diotadxup, hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5dcreate_f, (grpTransform,     "fiota", H5T_NATIVE_DOUBLE, filespace_fiota    , dset_id_fiota    , hdfier), __FILE__, __LINE__ )

  ! filespaces can be closed as soon as datasets are created
  H5CALL( sphdf5, h5sclose_f, (filespace_diotadxup, hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5sclose_f, (filespace_fiota    , hdfier), __FILE__, __LINE__ )

  ! Select hyperslab in the file.
  H5CALL( sphdf5, h5dget_space_f, (dset_id_diotadxup, filespace_diotadxup, hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5dget_space_f, (dset_id_fiota    , filespace_fiota    , hdfier), __FILE__, __LINE__ )

  ! Each process defines dataset in memory and writes it to the hyperslab in the file.
  H5CALL( sphdf5, h5screate_simple_f, (rankT, int((/2,1/),HSIZE_T), memspace_diotadxup, hdfier), __FILE__, __LINE__ )

 endif ! myid.eq.0

end subroutine init_flt_output

!> \brief Write a hyperslab of Poincare data corresponding to the output of one parallel worker.
!> \ingroup grp_output
!>
!> @param offset radial offset at which the data belongs
!> @param data output from field-line tracing
!> @param success flags to indicate if integrator was successful
subroutine write_poincare( offset, data, success )

  use allglobal, only : Nz
  use inputlist, only : nPpts

  LOCALS

  integer, intent(in) :: offset, success(:)
  REAL, intent(in)    :: data(:,:,:)
  integer(hsize_t), dimension(3) :: length
  integer(HSIZE_T), dimension(2) :: dims_singleTraj ! dimensions of single trajectory data

  BEGIN( sphdf5 )

 if (myid.eq.0 .and. .not.skip_write) then

  dims_singleTraj = (/ Nz, nPpts /)
  length          = (/ Nz, nPpts, 1 /)

  ! On entry, Fortran does not know that indexing in data is from 0 to Nz-1.
  ! Hence, use default indices 1:Nz in this routine
  H5CALL( sphdf5, h5sselect_hyperslab_f, (filespace_t, H5S_SELECT_SET_F, int((/0,0,offset/),HSSIZE_T), length, hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5dwrite_f, (dset_id_t, H5T_NATIVE_DOUBLE, data(1,1:Nz,1:nPpts), dims_singleTraj, hdfier, &
  &               file_space_id=filespace_t, mem_space_id=memspace_t ), __FILE__, __LINE__ )

  H5CALL( sphdf5, h5sselect_hyperslab_f, (filespace_s, H5S_SELECT_SET_F, int((/0,0,offset/),HSSIZE_T), length, hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5dwrite_f, (dset_id_s, H5T_NATIVE_DOUBLE, data(2,1:Nz,1:nPpts), dims_singleTraj, hdfier, &
  &               file_space_id=filespace_s, mem_space_id=memspace_s ), __FILE__, __LINE__ )

  H5CALL( sphdf5, h5sselect_hyperslab_f, (filespace_R, H5S_SELECT_SET_F, int((/0,0,offset/),HSSIZE_T), length, hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5dwrite_f, (dset_id_R, H5T_NATIVE_DOUBLE, data(3,1:Nz,1:nPpts), dims_singleTraj, hdfier, &
  &               file_space_id=filespace_R, mem_space_id=memspace_R ), __FILE__, __LINE__ )

  H5CALL( sphdf5, h5sselect_hyperslab_f, (filespace_Z, H5S_SELECT_SET_F, int((/0,0,offset/),HSSIZE_T), length, hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5dwrite_f, (dset_id_Z, H5T_NATIVE_DOUBLE, data(4,1:Nz,1:nPpts), dims_singleTraj, hdfier, &
  &               file_space_id=filespace_Z, mem_space_id=memspace_Z ), __FILE__, __LINE__ )

  H5CALL( sphdf5, h5sselect_hyperslab_f, (filespace_success, H5S_SELECT_SET_F, int((/offset/),HSSIZE_T), int((/1/), HSIZE_T), hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5dwrite_f, (dset_id_success, H5T_NATIVE_INTEGER, success, int((/1/), HSIZE_T), hdfier, &
  &               file_space_id=filespace_success, mem_space_id=memspace_success ), __FILE__, __LINE__ )

 endif ! myid.eq.0

end subroutine write_poincare

!> \brief Write the rotational transform output from field line following.
!> \ingroup grp_output
!>
!> @param offset radial offset at which the data belongs
!> @param length length of dataset to write
!> @param lvol nested volume index
!> @param diotadxup derivative of rotational transform (?)
!> @param fiota rotational transform
subroutine write_transform( offset, length, lvol, diotadxup, fiota )

  LOCALS
  INTEGER, intent(in) :: offset, length, lvol
  REAL, intent(in)    :: diotadxup(:), fiota(:,:)

  BEGIN( sphdf5 )

 if (myid.eq.0 .and. .not.skip_write) then

  H5CALL( sphdf5, h5sselect_hyperslab_f, (filespace_diotadxup, H5S_SELECT_SET_F, int((/0,lvol-1/),HSSIZE_T), int((/2,1/),HSSIZE_T), hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5dwrite_f, (dset_id_diotadxup, H5T_NATIVE_DOUBLE, diotadxup, int((/2,1/),HSSIZE_T), hdfier, &
  &               file_space_id=filespace_diotadxup, mem_space_id=memspace_diotadxup ), __FILE__, __LINE__ )

  ! length of fiota piece to write here may change, so open and close memspace each time a new hyperslab is written
  H5CALL( sphdf5, h5screate_simple_f, (rankT, int((/length,2/),HSIZE_T), memspace_fiota    , hdfier), __FILE__, __LINE__ )

  H5CALL( sphdf5, h5sselect_hyperslab_f, (filespace_fiota, H5S_SELECT_SET_F, int((/offset,0/),HSSIZE_T), int((/length,2/),HSSIZE_T), hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5dwrite_f, (dset_id_fiota, H5T_NATIVE_DOUBLE, fiota(1:length,1:2), int((/length,2/),HSSIZE_T), hdfier, &
  &               file_space_id=filespace_fiota, mem_space_id=memspace_fiota ), __FILE__, __LINE__ )

  H5CALL( sphdf5, h5sclose_f, (memspace_fiota, hdfier), __FILE__, __LINE__ )

 endif ! myid.eq.0

end subroutine write_transform

!> \brief Finalize Poincare output.
!> \ingroup grp_output
!>
!> This closes the still-open datasets related to field-line tracing,
!> which had to be kept open during the tracing to be able to write
!> the outputs directly when a given worker thread is finished.
subroutine finalize_flt_output

  LOCALS

  BEGIN( sphdf5 )

 if (myid.eq.0 .and. .not.skip_write) then

  ! close filespaces
  H5CALL( sphdf5, h5sclose_f, (filespace_t,         hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5sclose_f, (filespace_s,         hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5sclose_f, (filespace_R,         hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5sclose_f, (filespace_Z,         hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5sclose_f, (filespace_success,   hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5sclose_f, (filespace_diotadxup, hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5sclose_f, (filespace_fiota,     hdfier), __FILE__, __LINE__ )

  ! close dataspaces
  H5CALL( sphdf5, h5sclose_f, (memspace_t,          hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5sclose_f, (memspace_s,          hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5sclose_f, (memspace_R,          hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5sclose_f, (memspace_Z,          hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5sclose_f, (memspace_success,    hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5sclose_f, (memspace_diotadxup,  hdfier), __FILE__, __LINE__ )
  ! memspace_fiota is re-opened/closed in each iteration (see write_transform)

  ! close datasets
  H5CALL( sphdf5, h5dclose_f, (dset_id_t,           hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5dclose_f, (dset_id_s,           hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5dclose_f, (dset_id_R,           hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5dclose_f, (dset_id_Z,           hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5dclose_f, (dset_id_success,     hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5dclose_f, (dset_id_diotadxup,   hdfier), __FILE__, __LINE__ )
  H5CALL( sphdf5, h5dclose_f, (dset_id_fiota,       hdfier), __FILE__, __LINE__ )

  ! close groups
  HCLOSEGRP( grpPoincare  )
  HCLOSEGRP( grpTransform )

 endif ! myid.eq.0

end subroutine finalize_flt_output

!> \brief Write the magnetic vector potential Fourier harmonics to the output file group \c /vector_potential .
!> \ingroup grp_output
!>
!> The data is stacked in the radial direction over \c Lrad ,
!> since \c Lrad can be different in each volume, but HDF5 only supports
!> rectangular arrays. So, one needs to split the \c sumLrad dimension
!> into chunks given by the input \c Lrad array.
!>
!> @param sumLrad total sum over \c Lrad in all nested volumes
!> @param allAte \f$A^{\theta}_\mathrm{even}\f$ for all nested volumes
!> @param allAze \f$A^{\zeta}_\mathrm{even}\f$ for all nested volumes
!> @param allAto \f$A^{\theta}_\mathrm{odd}\f$ for all nested volumes
!> @param allAzo \f$A^{\zeta}_\mathrm{odd}\f$ for all nested volumes
subroutine write_vector_potential(sumLrad, allAte, allAze, allAto, allAzo)

  use allglobal, only : mn

  LOCALS
  integer, intent(in) :: sumLrad
  REAL, intent(in)    :: allAte(:,:), allAze(:,:), allAto(:,:), allAzo(:,:)
  integer(hid_t)      :: grpVectorPotential

  BEGIN( sphdf5 )

 if (myid.eq.0 .and. .not.skip_write) then

  HDEFGRP( file_id, vector_potential, grpVectorPotential )

  HWRITERA( grpVectorPotential, sumLrad, mn, Ate, allAte(1:sumLrad,1:mn) )
  HWRITERA( grpVectorPotential, sumLrad, mn, Aze, allAze(1:sumLrad,1:mn) )
  HWRITERA( grpVectorPotential, sumLrad, mn, Ato, allAto(1:sumLrad,1:mn) )
  HWRITERA( grpVectorPotential, sumLrad, mn, Azo, allAzo(1:sumLrad,1:mn) )

  HCLOSEGRP( grpVectorPotential )

 endif ! myid.eq.0

end subroutine write_vector_potential

!> \brief Write the final state of the equilibrium to the output file.
!> \ingroup grp_output
!>
subroutine hdfint

  use fileunits, only : ounit
  use inputlist
  use allglobal, only : ncpu, cpus, &
                        Mvol, ForceErr, BnsErr,&
                        mn, im, in, iRbc, iZbs, iRbs, iZbc, &
                        mns, ims, ins, &
                        dRbc, dZbs, dRbs, dZbc, &
                        vvolume, dvolume, &
                        Bsupumn, Bsupvmn, &
                        Btemn, Bzemn, Btomn, Bzomn, &
                        iVns, iBns, iVnc, iBnc, &
                        lmns, &
                        TT, &
                        beltramierror, &
                        IPDt, dlambdaout, lmns

  LOCALS

  INTEGER                        :: Mrad
  REAL                           :: tvolume

  integer(hid_t)                 :: grpOutput

  BEGIN( sphdf5 )

 if (myid.eq.0 .and. .not.skip_write) then

  HDEFGRP( file_id, output, grpOutput )

  HWRITERV( grpOutput,           mn,               Vns,       iVns(1:mn)   ) !     stellarator symmetric normal field at boundary; vacuum component;
  HWRITERV( grpOutput,           mn,               Bns,       iBns(1:mn)   ) !     stellarator symmetric normal field at boundary; plasma component;
  HWRITERV( grpOutput,           mn,               Vnc,       iVnc(1:mn)   ) ! non-stellarator symmetric normal field at boundary; vacuum component;
  HWRITERV( grpOutput,           mn,               Bnc,       iBnc(1:mn)   ) ! non-stellarator symmetric normal field at boundary; plasma component;

!> <ul>
!> <li> In addition to the input variables, which are described in global(), the following quantities are written to \c ext.sp.h5 :
!latex
!latex \begin{tabular}{|l|l|l|} \hline

!latex \type{variable}               & type    & \pb{description} \\ \hline

!latex \type{mn}                     & integer & \pb{number of Fourier modes} \\
  HWRITEIV( grpOutput,  1, mn, (/ mn /)  )
!latex \type{im(1:mn)}               & integer & \pb{poloidal mode numbers} \\
  HWRITEIV( grpOutput, mn, im, im(1:mn) )
!latex \type{in(1:mn)}               & integer & \pb{toroidal mode numbers} \\
  HWRITEIV( grpOutput, mn, in, in(1:mn) )
!latex \type{mns}                     & integer & \pb{number of Fourier modes} \\
  HWRITEIV( grpOutput,  1, mns, (/ mns /)  )
!latex \type{ims(1:mns)}               & integer & \pb{poloidal mode numbers} \\
  HWRITEIV( grpOutput, mns, ims, ims(1:mns) )
!latex \type{ins(1:mns)}               & integer & \pb{toroidal mode numbers} \\
  HWRITEIV( grpOutput, mns, ins, ins(1:mns) )
!latex \type{Mvol}                   & integer & \pb{number of interfaces = number of volumes} \\
  HWRITEIV( grpOutput,  1, Mvol, (/ Mvol /))
!latex \type{iRbc(1:mn,0:Mvol)}      & real    & \pb{Fourier harmonics, $R_{m,n}$, of interfaces} \\
  HWRITERA( grpOutput, mn, (Mvol+1), Rbc, iRbc(1:mn,0:Mvol) )
!latex \type{iZbs(1:mn,0:Mvol)}      & real    & \pb{Fourier harmonics, $Z_{m,n}$, of interfaces} \\
  HWRITERA( grpOutput, mn, (Mvol+1), Zbs, iZbs(1:mn,0:Mvol) )
!latex \type{iRbs(1:mn,0:Mvol)}      & real    & \pb{Fourier harmonics, $R_{m,n}$, of interfaces} \\
  HWRITERA( grpOutput, mn, (Mvol+1), Rbs, iRbs(1:mn,0:Mvol) )
!latex \type{iZbc(1:mn,0:Mvol)}      & real    & \pb{Fourier harmonics, $Z_{m,n}$, of interfaces} \\
  HWRITERA( grpOutput, mn, (Mvol+1), Zbc, iZbc(1:mn,0:Mvol) )
!latex \type{BnsErr}                 & real    & \pb{error in self-consistency of field on plasma boundary (in freeboundary)} \\
  HWRITERV( grpOutput, 1, BnsErr, (/ BnsErr /)) ! already in /input/global
!latex \type{ForceErr}               & real    & \pb{force-balance error across interfaces} \\
  HWRITERV( grpOutput,  1, ForceErr, (/ ForceErr /))
!latex \type{Ivolume}                & real    & \pb{Volume current at output (parallel, externally induced)}
  HWRITERV( grpOutput, Mvol, Ivolume, Ivolume(1:Mvol))
!latex \type{IPDt}                   & real    & \pb{Surface current at output}
  HWRITERV( grpOutput, Mvol, IPDt, IPDt(1:Mvol))

  ! the following quantites can be different from input value
  HWRITERV( grpOutput,   Mvol, adiabatic         , adiabatic(1:Nvol)   )
  HWRITERV( grpOutput,   Nvol, helicity          ,  helicity(1:Nvol)   )
  HWRITERV( grpOutput,   Mvol, mu                ,        mu(1:Mvol)   )
  HWRITERV( grpOutput,   Mvol, tflux             ,     tflux(1:Mvol)   )
  HWRITERV( grpOutput,   Mvol, pflux             ,     pflux(1:Mvol)   )

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


! Write lambda_mn, Fourier harmonics or transformation to straight field line coordinates.
  HWRITERC( grpOutput, lmns, Mvol, 2, lambdamn, dlambdaout(1:lmns,1:Mvol,0:1) )

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
!> </li>
!> <li> All quantities marked as real should be treated as double precision. </li>
!> </ul>

  HCLOSEGRP( grpOutput )

 endif ! myid.eq.0

end subroutine hdfint

!> \brief Close all open HDF5 objects (we know of) and list any remaining still-open objects.
!> \ingroup grp_output
!>
subroutine finish_outfile
! Close all open HDF5 objects (we know of) and list any remaining still-open objects
! The goal should be to close all objects specifically!

  LOCALS
  integer(size_t)                          :: obj_count                                  ! number of open HDF5 objects
  integer(size_t)                          :: num_objs    ! number of still-open objects
  integer(hid_t),dimension(:),allocatable  :: obj_ids    ! still-open objects
  integer                                  :: iObj
  integer(size_t)                          :: openLength
  character(len=:),allocatable             :: openName
  integer(size_t),parameter                :: dummySize=1
  character(len=dummySize+1)               :: dummyName
  integer                                  :: typeClass

  BEGIN( sphdf5 )

 if (myid.eq.0 .and. .not.skip_write) then

  ! close objects related to convergence output
  H5CALL( sphdf5, h5tclose_f, (dt_nDcalls_id, hdfier)    , __FILE__, __LINE__)
  H5CALL( sphdf5, h5tclose_f, (dt_Energy_id, hdfier)     , __FILE__, __LINE__)
  H5CALL( sphdf5, h5tclose_f, (dt_ForceErr_id, hdfier)   , __FILE__, __LINE__)
  H5CALL( sphdf5, h5tclose_f, (dt_iRbc_id, hdfier)       , __FILE__, __LINE__)
  H5CALL( sphdf5, h5tclose_f, (dt_iZbs_id, hdfier)       , __FILE__, __LINE__)
  H5CALL( sphdf5, h5tclose_f, (dt_iRbs_id, hdfier)       , __FILE__, __LINE__)
  H5CALL( sphdf5, h5tclose_f, (dt_iZbc_id, hdfier)       , __FILE__, __LINE__)
  H5CALL( sphdf5, h5dclose_f, (iteration_dset_id, hdfier), __FILE__, __LINE__) ! End access to the dataset and release resources used by it.
  H5CALL( sphdf5, h5pclose_f, (plist_id, hdfier)         , __FILE__, __LINE__) ! close plist used for 'preserve' flag (does not show up in obj_count below)

  ! check whether we forgot to close some resources; only check for group, dataset and datatype (there is only one file and that should be still open...)
  H5CALL( sphdf5, h5fget_obj_count_f, (file_id, ior(H5F_OBJ_GROUP_F, ior(H5F_OBJ_DATASET_F, H5F_OBJ_DATATYPE_F)), obj_count, hdfier), __FILE__, __LINE__ )

  if (obj_count.gt.0) then
    write(*,'("There are still ",i3," hdf5 objects open")') obj_count
    allocate(obj_ids(1:obj_count))

    ! groups
    H5CALL( sphdf5, h5fget_obj_ids_f, (file_id, H5F_OBJ_GROUP_F, obj_count, obj_ids, hdfier, num_objs), __FILE__, __LINE__) ! get for open objects
    if (num_objs.gt.0) then
      write(*,'("There are still ",i3," HDF5 groups open:")') num_objs
      do iObj=1,num_objs
        openLength=0
        H5CALL( sphdf5, h5iget_name_f, (obj_ids(iObj), dummyName, dummySize, openLength, hdfier), __FILE__, __LINE__)
        allocate(character(len=openLength+1) :: openName)
        H5CALL( sphdf5, h5iget_name_f, (obj_ids(iObj), openName, openLength, openLength, hdfier), __FILE__, __LINE__)
        write(*,*) openName
        deallocate(openName)

        H5CALL( sphdf5, h5gclose_f, (obj_ids(iObj), hdfier), __FILE__, __LINE__)
      enddo
    endif

    ! datasets
    H5CALL( sphdf5, h5fget_obj_ids_f, (file_id, H5F_OBJ_DATASET_F, obj_count, obj_ids, hdfier, num_objs), __FILE__, __LINE__) ! get for open objects
    if (num_objs.gt.0) then
      write(*,'("There are still ",i3," HDF5 datasets open:")') num_objs
      do iObj=1,num_objs
        openLength=0
        H5CALL( sphdf5, h5iget_name_f, (obj_ids(iObj), dummyName, dummySize, openLength, hdfier), __FILE__, __LINE__)
        allocate(character(len=openLength+1) :: openName)
        H5CALL( sphdf5, h5iget_name_f, (obj_ids(iObj), openName, openLength, openLength, hdfier), __FILE__, __LINE__)
        write(*,*) openName(1:openLength)
        deallocate(openName)

        H5CALL( sphdf5, h5dclose_f, (obj_ids(iObj), hdfier), __FILE__, __LINE__)
      enddo
    endif

    ! datatypes
    H5CALL( sphdf5, h5fget_obj_ids_f, (file_id, H5F_OBJ_DATATYPE_F, obj_count, obj_ids, hdfier, num_objs), __FILE__, __LINE__) ! get for open objects
    if (num_objs.gt.0) then
      write(*,'("There are still ",i3," HDF5 datatypes open:")') num_objs
      do iObj=1,num_objs
        H5CALL( sphdf5, h5tget_class_f, (obj_ids(iObj), typeClass, hdfier), __LINE__, __FILE__) ! determine class of open datatype
        if      (typeClass.eq.H5T_NO_CLASS_F ) then ; write(*,*) "H5T_NO_CLASS_F"
        else if (typeClass.eq.H5T_INTEGER_F  ) then ; write(*,*) "H5T_INTEGER_F"
        else if (typeClass.eq.H5T_FLOAT_F    ) then ; write(*,*) "H5T_FLOAT_F"
        else if (typeClass.eq.H5T_STRING_F   ) then ; write(*,*) "H5T_STRING_F"
        else if (typeClass.eq.H5T_BITFIELD_F ) then ; write(*,*) "H5T_BITFIELD_F"
        else if (typeClass.eq.H5T_OPAQUE_F   ) then ; write(*,*) "H5T_OPAQUE_F"
        else if (typeClass.eq.H5T_COMPOUND_F ) then ; write(*,*) "H5T_COMPOUND_F"
        else if (typeClass.eq.H5T_REFERENCE_F) then ; write(*,*) "H5T_REFERENCE_F"
        else if (typeClass.eq.H5T_ENUM_F     ) then ; write(*,*) "H5T_ENUM_F"
        else if (typeClass.eq.H5T_VLEN_F     ) then ; write(*,*) "H5T_VLEN_F"
        else if (typeClass.eq.H5T_ARRAY_F    ) then ; write(*,*) "H5T_ARRAY_F"
        else ; write(*,*) "UNKNOWN TYPE!"
        endif

        H5CALL( sphdf5, h5tclose_f, (obj_ids(iObj), hdfier), __FILE__, __LINE__)
      enddo
    endif

    deallocate(obj_ids)
  endif ! (obj_count.gt.0)

  H5CALL( sphdf5, h5fclose_f, ( file_id, hdfier ), __FILE__, __LINE__ ) ! terminate access on output file;
  H5CALL( sphdf5, h5close_f,  ( hdfier ),          __FILE__, __LINE__ ) ! close Fortran interface to the HDF5 library;

 endif ! myid.eq.0

end subroutine finish_outfile

end module sphdf5
