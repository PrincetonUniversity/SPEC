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
    use mod_kinds, only: wp => dp
    use fileunits, only: ounit
    use allglobal, only: myid, cpus, MPI_COMM_SPEC, ext, skip_write
    use constants, only: version
    use hdf5
    use h5utils

    implicit none

!  integer                        :: hdfier              !< error flag for HDF5 library
!   integer                        :: rank                !< rank of data to write using macros
    integer(hid_t) :: file_id             !< default file ID used in macros
!   integer(hid_t)                 :: space_id            !< default dataspace ID used in macros
!   integer(hid_t)                 :: dset_id             !< default dataset ID used in macros
!   integer(hsize_t)               :: onedims(1:1)        !< dimension specifier for one-dimensional data used in macros
!   integer(hsize_t)               :: twodims(1:2)        !< dimension specifier for two-dimensional data used in macros
!   integer(hsize_t)               :: threedims(1:3)      !< dimension specifier for three-dimensional data used in macros
!  logical                        :: grp_exists          !< flags used to signal if a group already exists
!  logical                        :: var_exists          !< flags used to signal if a variable already exists

    integer(hid_t) :: iteration_dset_id   !< Dataset identifier for "iteration"
    integer(hid_t) :: dataspace           !< dataspace for extension by 1 iteration object
    integer(hid_t) :: memspace            !< memspace for extension by 1 iteration object
    integer(hsize_t), dimension(1) :: old_data_dims       !< current dimensions of "iterations" dataset
    integer(hsize_t), dimension(1) :: data_dims           !< new dimensions for "iterations" dataset
    integer(hsize_t), dimension(1) :: max_dims            !< maximum dimensions for "iterations" dataset
    integer(hid_t) :: plist_id            !< Property list identifier used to activate dataset transfer property
    integer(hid_t) :: dt_nDcalls_id       !< Memory datatype identifier (for "nDcalls"  dataset in "/grid")
    integer(hid_t) :: dt_Energy_id        !< Memory datatype identifier (for "Energy"   dataset in "/grid")
    integer(hid_t) :: dt_ForceErr_id      !< Memory datatype identifier (for "ForceErr" dataset in "/grid")
    integer(hid_t) :: dt_iRbc_id          !< Memory datatype identifier (for "iRbc"     dataset in "/grid")
    integer(hid_t) :: dt_iZbs_id          !< Memory datatype identifier (for "iZbs"     dataset in "/grid")
    integer(hid_t) :: dt_iRbs_id          !< Memory datatype identifier (for "iRbs"     dataset in "/grid")
    integer(hid_t) :: dt_iZbc_id          !< Memory datatype identifier (for "iZbc"     dataset in "/grid")

    integer, parameter :: rankP = 3             !< rank of Poincare data
    integer, parameter :: rankT = 2             !< rank of rotational transform data

    integer(hid_t) :: grpPoincare         !< group for Poincare data
    integer(HID_T) :: dset_id_t           !< Dataset identifier for \f$\theta\f$ coordinate of field line following
    integer(HID_T) :: dset_id_s           !< Dataset identifier for \f$s\f$ coordinate of field line following
    integer(HID_T) :: dset_id_R           !< Dataset identifier for \f$R\f$ coordinate of field line following
    integer(HID_T) :: dset_id_Z           !< Dataset identifier for \f$Z\f$ coordinate of field line following
    integer(HID_T) :: dset_id_success     !< Dataset identifier for success flag of trajectories to follow
    integer(HID_T) :: filespace_t         !< Dataspace identifier in file for \f$\theta\f$ coordinate of field line following
    integer(HID_T) :: filespace_s         !< Dataspace identifier in file for \f$s\f$ coordinate of field line following
    integer(HID_T) :: filespace_R         !< Dataspace identifier in file for \f$R\f$ coordinate of field line following
    integer(HID_T) :: filespace_Z         !< Dataspace identifier in file for \f$Z\f$ coordinate of field line following
    integer(HID_T) :: filespace_success   !< Dataspace identifier in file for success flag of trajectories to follow
    integer(HID_T) :: memspace_t          !< Dataspace identifier in memory for \f$\theta\f$ coordinate of field line following
    integer(HID_T) :: memspace_s          !< Dataspace identifier in memory for \f$s\f$ coordinate of field line following
    integer(HID_T) :: memspace_R          !< Dataspace identifier in memory for \f$R\f$ coordinate of field line following
    integer(HID_T) :: memspace_Z          !< Dataspace identifier in memory for \f$Z\f$ coordinate of field line following
    integer(HID_T) :: memspace_success    !< Dataspace identifier in memory for success flag of trajectories to follow

    integer(hid_t) :: grpTransform        !< group for rotational transform data
    integer(HID_T) :: dset_id_diotadxup   !< Dataset identifier for diotadxup (derivative of rotational transform ?)
    integer(HID_T) :: dset_id_fiota       !< Dataset identifier for fiota     (              rotational transform ?)
    integer(HID_T) :: filespace_diotadxup !< Dataspace identifier in file for diotadxup
    integer(HID_T) :: filespace_fiota     !< Dataspace identifier in file for fiota
    integer(HID_T) :: memspace_diotadxup  !< Dataspace identifier in memory for diotadxup
    integer(HID_T) :: memspace_fiota      !< Dataspace identifier in memory for fiota

    private
    public init_outfile, &
        mirror_input_to_outfile, &
        init_convergence_output, &
        write_convergence_output, &
        write_grid, &
        init_flt_output, &
        write_poincare, &
        write_transform, &
        finalize_flt_output, &
        write_vector_potential, &
        hdfint, &
        finish_outfile

contains

!> \brief Initialize the interface to the HDF5 library and open the output file.
!> \ingroup grp_output
!>
    subroutine init_outfile()
        integer :: hdfier !< error flag for HDF5 library
        integer(hid_t) :: dset_id

        if (myid .eq. 0 .and. .not. skip_write) then

            ! initialize Fortran interface to the HDF5 library
            call h5open_f(hdfier)

            ! (en/dis)able HDF5 internal error messages;
            ! sphdf5 has its own error messages coming from the macros
            call h5eset_auto_f(internalHdf5Msg, hdfier)

            ! Create the file
            call h5fcreate_f(trim(ext)//".sp.h5", H5F_ACC_TRUNC_F, file_id, hdfier)

            ! write version number
            call HWRITERV_LO(file_id, "version", 1, (/version/), dset_id)
            call H5DESCR_CDSET(dset_id, "version of SPEC")

        end if ! myid.eq.0
    end subroutine ! init_outfile

!> \brief Mirror input variables into output file.
!> \ingroup grp_output
!>
!> The goal of this routine is to have an exact copy of the input file contents
!> that were used to parameterize a given SPEC run.
!> This also serves to check after the run if SPEC correctly understood the text-based input file.
    subroutine mirror_input_to_outfile

        use inputlist
        use allglobal, only: myid, Mvol, skip_write
        use h5utils

        integer(hid_t) :: dset_id
        integer(hid_t) :: grpInput
        integer(hid_t) :: grpInputPhysics, grpInputNumerics, grpInputLocal, grpInputGlobal, grpInputDiagnostics

        if (myid .eq. 0 .and. .not. skip_write) then

            call HDEFGRP(file_id, "input", grpInput)
            call H5DESCR(grpInput, "group for mirrored input data")

! the following variables constitute the namelist/physicslist/; note that all variables in namelist need to be broadcasted in readin;
! they go into ext.h5/input/physics

            ! the calls used here work as follows:
            ! step 1. HWRITEIV_LO e.g. write(s an) i(nteger) v(ariable) and l(eaves) o(pen) the dataset, so that in
            ! step 2a. an attribute with descr(iptive) information can be attached to the dataset and finally, in
            ! step 2b. the attribute is closed and also we c(lose the) d(ata)set.

            call HDEFGRP(grpInput, "physics", grpInputPhysics)
            call H5DESCR(grpInputPhysics, "physics inputs")

            call HWRITEIV_LO(grpInputPhysics, "Igeometry", 1, (/Igeometry/), dset_id)
            call H5DESCR_CDSET(dset_id, "geometry identifier")
            call HWRITEIV_LO(grpInputPhysics, "Istellsym", 1, (/Istellsym/), dset_id)
            call H5DESCR_CDSET(dset_id, "stellarator symmetry flag")
            call HWRITEIV_LO(grpInputPhysics, "Lfreebound", 1, (/Lfreebound/), dset_id)
            call H5DESCR_CDSET(dset_id, "free boundary flag")
            call HWRITERV_LO(grpInputPhysics, "phiedge", 1, (/phiedge/), dset_id)
            call H5DESCR_CDSET(dset_id, "total enclosed toroidal magnetic flux")
            call HWRITERV_LO(grpInputPhysics, "curtor", 1, (/curtor/), dset_id)
            call H5DESCR_CDSET(dset_id, "total enclosed toroidal current")
            call HWRITERV_LO(grpInputPhysics, "curpol", 1, (/curpol/), dset_id)
            call H5DESCR_CDSET(dset_id, "total enclosed poloidal current")
            call HWRITERV_LO(grpInputPhysics, "gamma", 1, (/gamma/), dset_id)
            call H5DESCR_CDSET(dset_id, "adiabatic index")
            call HWRITEIV_LO(grpInputPhysics, "Nfp", 1, (/Nfp/), dset_id)
            call H5DESCR_CDSET(dset_id, "number of stellarator field periods")
            call HWRITEIV_LO(grpInputPhysics, "Nvol", 1, (/Nvol/), dset_id)
            call H5DESCR_CDSET(dset_id, "number of volumes")
            call HWRITEIV_LO(grpInputPhysics, "Mpol", 1, (/Mpol/), dset_id)
            call H5DESCR_CDSET(dset_id, "maximum poloidal mode number")
            call HWRITEIV_LO(grpInputPhysics, "Ntor", 1, (/Ntor/), dset_id)
            call H5DESCR_CDSET(dset_id, "maximum toroidal mode number")
            call HWRITEIV_LO(grpInputPhysics, "Lrad", Mvol, Lrad(1:Mvol), dset_id)
            call H5DESCR_CDSET(dset_id, "degree of radial Chebychev polynomials")
            call HWRITEIV_LO(grpInputPhysics, "Lconstraint", 1, (/Lconstraint/), dset_id)
            call H5DESCR_CDSET(dset_id, "type of constraint to enforce")
            call HWRITEIV_LO(grpInputPhysics, "Lreflect", 1, (/Lreflect/), dset_id)
            call H5DESCR_CDSET(dset_id, "whether to reflect the perturbation on both boundaries for slab geometry")
            call HWRITERV_LO(grpInputPhysics, "tflux", Mvol, tflux(1:Mvol), dset_id)
            call H5DESCR_CDSET(dset_id, "toroidal magnetic flux in volumes")
            call HWRITERV_LO(grpInputPhysics, "pflux", Mvol, pflux(1:Mvol), dset_id)
            call H5DESCR_CDSET(dset_id, "poloidal magnetic flux in volumes")
            call HWRITERV_LO(grpInputPhysics, "helicity", Nvol, helicity(1:Nvol), dset_id)
            call H5DESCR_CDSET(dset_id, "helicity profile")
            call HWRITERV_LO(grpInputPhysics, "pscale", 1, (/pscale/), dset_id)
            call H5DESCR_CDSET(dset_id, "scaling factor for pressure")
            call HWRITERV_LO(grpInputPhysics, "pressure", Nvol, pressure(1:Nvol), dset_id)
            call H5DESCR_CDSET(dset_id, "pressure profile")
            call HWRITEIV_LO(grpInputPhysics, "Ladiabatic", 1, (/Ladiabatic/), dset_id)
            call H5DESCR_CDSET(dset_id, "adiabatic flag")
            call HWRITERV_LO(grpInputPhysics, "adiabatic", Mvol, adiabatic(1:Nvol), dset_id)
            call H5DESCR_CDSET(dset_id, "adiabatic profile (?)")
            call HWRITERV_LO(grpInputPhysics, "mu", (1 + Nvol), mu(1:Mvol), dset_id)
            call H5DESCR_CDSET(dset_id, "Beltrami parameter, parallel current profile")
            call HWRITERV_LO(grpInputPhysics, "Ivolume", (1 + Nvol), Ivolume(1:Mvol), dset_id)
            call H5DESCR_CDSET(dset_id, "Volume current, externally driven, parallel current profile")
            call HWRITERV_LO(grpInputPhysics, "Isurf", (Mvol), Isurf(1:Mvol), dset_id)
            call H5DESCR_CDSET(dset_id, "Surface current, currents that are not volume currents (pressure driven, shielding currents)")
            call HWRITEIV_LO(grpInputPhysics, "pl", (1 + Mvol), pl(0:Nvol), dset_id)
            call H5DESCR_CDSET(dset_id, "pl ?")
            call HWRITEIV_LO(grpInputPhysics, "ql", (1 + Mvol), ql(0:Nvol), dset_id)
            call H5DESCR_CDSET(dset_id, "ql ?")
            call HWRITEIV_LO(grpInputPhysics, "pr", (1 + Mvol), pr(0:Nvol), dset_id)
            call H5DESCR_CDSET(dset_id, "pr ?")
            call HWRITEIV_LO(grpInputPhysics, "qr", (1 + Mvol), qr(0:Nvol), dset_id)
            call H5DESCR_CDSET(dset_id, "qr ?")
            call HWRITERV_LO(grpInputPhysics, "iota", (1 + Nvol), iota(0:Nvol), dset_id)
            call H5DESCR_CDSET(dset_id, "rotational transform profile on inside of ideal interfaces")
            call HWRITEIV_LO(grpInputPhysics, "lp", (1 + Mvol), lp(0:Nvol), dset_id)
            call H5DESCR_CDSET(dset_id, "lp ?")
            call HWRITEIV_LO(grpInputPhysics, "lq", (1 + Mvol), lq(0:Nvol), dset_id)
            call H5DESCR_CDSET(dset_id, "lq ?")
            call HWRITEIV_LO(grpInputPhysics, "rp", (1 + Mvol), rp(0:Nvol), dset_id)
            call H5DESCR_CDSET(dset_id, "rp ?")
            call HWRITEIV_LO(grpInputPhysics, "rq", (1 + Mvol), rq(0:Nvol), dset_id)
            call H5DESCR_CDSET(dset_id, "rq ?")
            call HWRITERV_LO(grpInputPhysics, "oita", (1 + Nvol), oita(0:Nvol), dset_id)
            call H5DESCR_CDSET(dset_id, "rotational transform profile on outside of ideal interfaces")
            call HWRITERV_LO(grpInputPhysics, "rtor", 1, (/rtor/), dset_id)
            call H5DESCR_CDSET(dset_id, "for aspect ratio in slab")
            call HWRITERV_LO(grpInputPhysics, "rpol", 1, (/rpol/), dset_id)
            call H5DESCR_CDSET(dset_id, "for aspect ratio in slab")

            call HWRITERV_LO(grpInputPhysics, "Rac", (1 + Ntor), Rac(0:Ntor), dset_id)
            call H5DESCR_CDSET(dset_id, "stellarator symmetric coordinate axis R cosine Fourier coefficients")
            call HWRITERV_LO(grpInputPhysics, "Zas", (1 + Ntor), Zas(0:Ntor), dset_id)
            call H5DESCR_CDSET(dset_id, "stellarator symmetric coordinate axis Z   sine Fourier coefficients")
            call HWRITERV_LO(grpInputPhysics, "Ras", (1 + Ntor), Ras(0:Ntor), dset_id)
            call H5DESCR_CDSET(dset_id, "non-stellarator symmetric coordinate axis R   sine Fourier coefficients")
            call HWRITERV_LO(grpInputPhysics, "Zac", (1 + Ntor), Zac(0:Ntor), dset_id)
            call H5DESCR_CDSET(dset_id, "non-stellarator symmetric coordinate axis Z cosine Fourier coefficients")

            call HWRITERA_LO(grpInputPhysics, "Rbc", (2*Ntor + 1), (2*Mpol + 1), Rbc(-Ntor:Ntor, -Mpol:Mpol), dset_id)
            call H5DESCR_CDSET(dset_id, "stellarator symmetric boundary R cosine Fourier coefficients")
            call HWRITERA_LO(grpInputPhysics, "Zbs", (2*Ntor + 1), (2*Mpol + 1), Zbs(-Ntor:Ntor, -Mpol:Mpol), dset_id)
            call H5DESCR_CDSET(dset_id, "stellarator symmetric boundary Z   sine Fourier coefficients")
            call HWRITERA_LO(grpInputPhysics, "Rbs", (2*Ntor + 1), (2*Mpol + 1), Rbs(-Ntor:Ntor, -Mpol:Mpol), dset_id)
            call H5DESCR_CDSET(dset_id, "non-stellarator symmetric boundary R   sine Fourier coefficients")
            call HWRITERA_LO(grpInputPhysics, "Zbc", (2*Ntor + 1), (2*Mpol + 1), Zbc(-Ntor:Ntor, -Mpol:Mpol), dset_id)
            call H5DESCR_CDSET(dset_id, "non-stellarator symmetric boundary Z cosine Fourier coefficients")

            call HWRITERA_LO(grpInputPhysics, "Rwc", (2*Ntor + 1), (2*Mpol + 1), Rwc(-Ntor:Ntor, -Mpol:Mpol), dset_id)
            call H5DESCR_CDSET(dset_id, "stellarator symmetric boundary R cosine Fourier coefficients of wall")
            call HWRITERA_LO(grpInputPhysics, "Zws", (2*Ntor + 1), (2*Mpol + 1), Zws(-Ntor:Ntor, -Mpol:Mpol), dset_id)
            call H5DESCR_CDSET(dset_id, "stellarator symmetric boundary Z   sine Fourier coefficients of wall")
            call HWRITERA_LO(grpInputPhysics, "Rws", (2*Ntor + 1), (2*Mpol + 1), Rws(-Ntor:Ntor, -Mpol:Mpol), dset_id)
            call H5DESCR_CDSET(dset_id, "non-stellarator symmetric boundary R   sine Fourier coefficients of wall")
            call HWRITERA_LO(grpInputPhysics, "Zwc", (2*Ntor + 1), (2*Mpol + 1), Zwc(-Ntor:Ntor, -Mpol:Mpol), dset_id)
            call H5DESCR_CDSET(dset_id, "non-stellarator symmetric boundary Z cosine Fourier coefficients of wall")

            call HWRITERA_LO(grpInputPhysics, "Vns", (2*Ntor + 1), (2*Mpol + 1), Vns(-Ntor:Ntor, -Mpol:Mpol), dset_id)
            call H5DESCR_CDSET(dset_id, "stellarator symmetric normal field   sine Fourier coefficients at boundary; vacuum component")
            call HWRITERA_LO(grpInputPhysics, "Bns", (2*Ntor + 1), (2*Mpol + 1), Bns(-Ntor:Ntor, -Mpol:Mpol), dset_id)
            call H5DESCR_CDSET(dset_id, "stellarator symmetric normal field   sine Fourier coefficients at boundary; plasma component")
            call HWRITERA_LO(grpInputPhysics, "Vnc", (2*Ntor + 1), (2*Mpol + 1), Vnc(-Ntor:Ntor, -Mpol:Mpol), dset_id)
            call H5DESCR_CDSET(dset_id, "non-stellarator symmetric normal field cosine Fourier coefficients at boundary; vacuum component")
            call HWRITERA_LO(grpInputPhysics, "Bnc", (2*Ntor + 1), (2*Mpol + 1), Bnc(-Ntor:Ntor, -Mpol:Mpol), dset_id)
            call H5DESCR_CDSET(dset_id, "non-stellarator symmetric normal field cosine Fourier coefficients at boundary; plasma component")

            call HWRITERV_LO(grpInputPhysics, "mupftol", 1, (/mupftol/), dset_id)
            call H5DESCR_CDSET(dset_id, "mupftol")
            call HWRITEIV_LO(grpInputPhysics, "mupfits", 1, (/mupfits/), dset_id)
            call H5DESCR_CDSET(dset_id, "mupfits")

            call HCLOSEGRP(grpInputPhysics)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! the following variables constitute the namelist/numericlist/; note that all variables in namelist need to be broadcasted in readin;
! they go into ext.h5/input/numerics

            call HDEFGRP(grpInput, "numerics", grpInputNumerics)

            call HWRITEIV(grpInputNumerics, "Linitialize", 1, (/Linitialize/))
            call HWRITEIV(grpInputNumerics, "Lzerovac", 1, (/Lzerovac/))
            call HWRITEIV(grpInputNumerics, "Ndiscrete", 1, (/Ndiscrete/))
            call HWRITEIV(grpInputNumerics, "Nquad", 1, (/Nquad/))
            call HWRITEIV(grpInputNumerics, "iMpol", 1, (/iMpol/))
            call HWRITEIV(grpInputNumerics, "iNtor", 1, (/iNtor/))
            call HWRITEIV(grpInputNumerics, "Lsparse", 1, (/Lsparse/))
            call HWRITEIV(grpInputNumerics, "Lsvdiota", 1, (/Lsvdiota/))
            call HWRITEIV(grpInputNumerics, "imethod", 1, (/imethod/))
            call HWRITEIV(grpInputNumerics, "iorder", 1, (/iorder/))
            call HWRITEIV(grpInputNumerics, "iprecon", 1, (/iprecon/))
            call HWRITERV(grpInputNumerics, "iotatol", 1, (/iotatol/))
            call HWRITEIV(grpInputNumerics, "Lextrap", 1, (/Lextrap/))
            call HWRITEIV(grpInputNumerics, "Mregular", 1, (/Mregular/))
            call HWRITEIV(grpInputNumerics, "Lrzaxis", 1, (/Lrzaxis/))
            call HWRITEIV(grpInputNumerics, "Ntoraxis", 1, (/Ntoraxis/))

            call HCLOSEGRP(grpInputNumerics)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! the following variables constitute the namelist/locallist/; note that all variables in namelist need to be broadcasted in readin;
! they go into ext.h5/input/local

            call HDEFGRP(grpInput, "local", grpInputLocal)

            call HWRITEIV(grpInputLocal, "LBeltrami", 1, (/LBeltrami/))
            call HWRITEIV(grpInputLocal, "Linitgues", 1, (/Linitgues/))
            call HWRITEIV(grpInputLocal, "Lposdef", 1, (/Lposdef/)) ! redundant;
            call HWRITERV(grpInputLocal, "maxrndgues", 1, (/maxrndgues/))
            call HWRITEIV(grpInputLocal, "Lmatsolver", 1, (/Lmatsolver/))
            call HWRITEIV(grpInputLocal, "LGMRESprec", 1, (/LGMRESprec/))
            call HWRITERV(grpInputLocal, "epsGMRES", 1, (/epsGMRES/))
            call HWRITERV(grpInputLocal, "epsILU", 1, (/epsILU/))

            call HCLOSEGRP(grpInputLocal)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! the following variables constitute the namelist/globallist/; note that all variables in namelist need to be broadcasted in readin;
! they go into ext.h5/input/global

            call HDEFGRP(grpInput, "global", grpInputGlobal)

            call HWRITEIV(grpInputGlobal, "Lfindzero", 1, (/Lfindzero/))
            call HWRITERV(grpInputGlobal, "escale", 1, (/escale/))
            call HWRITERV(grpInputGlobal, "opsilon", 1, (/opsilon/))
            call HWRITERV(grpInputGlobal, "pcondense", 1, (/pcondense/))
            call HWRITERV(grpInputGlobal, "epsilon", 1, (/epsilon/))
            call HWRITERV(grpInputGlobal, "wpoloidal", 1, (/wpoloidal/))
            call HWRITERV(grpInputGlobal, "upsilon", 1, (/upsilon/))
            call HWRITERV(grpInputGlobal, "forcetol", 1, (/forcetol/))
            call HWRITERV(grpInputGlobal, "c05xmax", 1, (/c05xmax/))
            call HWRITERV(grpInputGlobal, "c05xtol", 1, (/c05xtol/))
            call HWRITERV(grpInputGlobal, "c05factor", 1, (/c05factor/))
            call HWRITELV(grpInputGlobal, "LreadGF", 1, (/LreadGF/))
            call HWRITEIV(grpInputGlobal, "mfreeits", 1, (/mfreeits/))
            call HWRITERV(grpInputGlobal, "bnstol", 1, (/bnstol/))  ! redundant;
            call HWRITERV(grpInputGlobal, "bnsblend", 1, (/bnsblend/))  ! redundant;
            call HWRITERV(grpInputGlobal, "gBntol", 1, (/gBntol/))
            call HWRITERV(grpInputGlobal, "gBnbld", 1, (/gBnbld/))
            call HWRITERV(grpInputGlobal, "vcasingeps", 1, (/vcasingeps/))
            call HWRITERV(grpInputGlobal, "vcasingtol", 1, (/vcasingtol/))
            call HWRITEIV(grpInputGlobal, "vcasingits", 1, (/vcasingits/))
            call HWRITEIV(grpInputGlobal, "vcasingper", 1, (/vcasingper/))
            call HWRITEIV(grpInputGlobal, "mcasingcal", 1, (/mcasingcal/))  ! redundant;

            call HCLOSEGRP(grpInputGlobal)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! the following variables constitute the namelist/diagnosticslist/; note that all variables in namelist need to be broadcasted in readin;
! they go into ext.h5/input/diagnostics

            call HDEFGRP(grpInput, "diagnostics", grpInputDiagnostics)

            call HWRITERV(grpInputDiagnostics, "odetol", 1, (/odetol/))
            call HWRITERV(grpInputDiagnostics, "absreq", 1, (/absreq/))           ! redundant;
            call HWRITERV(grpInputDiagnostics, "relreq", 1, (/relreq/))           ! redundant;
            call HWRITERV(grpInputDiagnostics, "absacc", 1, (/absacc/))           ! redundant;
            call HWRITERV(grpInputDiagnostics, "epsr", 1, (/epsr/))           ! redundant;
            call HWRITEIV(grpInputDiagnostics, "nPpts", 1, (/nPpts/))
            call HWRITERV(grpInputDiagnostics, "Ppts", 1, (/Ppts/))
            call HWRITEIV(grpInputDiagnostics, "nPtrj", Mvol, nPtrj(1:Mvol))
            call HWRITELV(grpInputDiagnostics, "LHevalues", 1, (/LHevalues/))
            call HWRITELV(grpInputDiagnostics, "LHevectors", 1, (/LHevectors/))
            call HWRITELV(grpInputDiagnostics, "LHmatrix", 1, (/LHmatrix/))
            call HWRITEIV(grpInputDiagnostics, "Lperturbed", 1, (/Lperturbed/))
            call HWRITEIV(grpInputDiagnostics, "dpp", 1, (/dpp/))
            call HWRITEIV(grpInputDiagnostics, "dqq", 1, (/dqq/))
            call HWRITEIV(grpInputDiagnostics, "Lcheck", 1, (/Lcheck/))
            call HWRITELV(grpInputDiagnostics, "Ltiming", 1, (/Ltiming/))
            call HWRITEIV(grpInputDiagnostics, "Lerrortype", 1, (/Lerrortype/))
            call HWRITEIV(grpInputDiagnostics, "Ngrid", 1, (/Ngrid/))
            call HWRITERV(grpInputDiagnostics, "fudge", 1, (/fudge/))          ! redundant;
            call HWRITERV(grpInputDiagnostics, "scaling", 1, (/scaling/))          ! redundant;

            call HCLOSEGRP(grpInputDiagnostics)

            call HCLOSEGRP(grpInput)

        end if ! myid.eq.0

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

        use allglobal, only: mn, Mvol

#ifdef OPENMP
        USE OMP_LIB
#endif
        use mpi
        implicit none
        integer :: ierr, astat, ios, nthreads, ithread
        real(wp) :: cput, cpui, cpuo = 0 ! cpu time; cpu initial; cpu old; 31 Jan 13;

        integer(hid_t) :: iteration_dspace_id  !< dataspace for "iteration"
        integer(hid_t) :: iteration_dtype_id   !< Compound datatype for "iteration"
        integer(hid_t) :: iRZbscArray_id       !< Memory datatype identifier
        integer(size_t) :: iteration_dtype_size !< Size of the "iteration" datatype
        integer(size_t) :: type_size_i          !< Size of the integer datatype
        integer(size_t) :: type_size_d          !< Size of the double precision datatype
        integer(size_t) :: offset               !< Member's offset
        integer(hid_t) :: crp_list             !< Dataset creation property identifier
        integer, parameter :: rank = 1             !< logging rank: convergence logging is one-dimensional
        integer(hsize_t), dimension(rank) :: maxdims              !< convergence logging maximum dimensions => will be unlimited
        integer(hsize_t), dimension(rank) :: dims = (/0/)       !< current convergence logging dimensions
        integer(hsize_t), dimension(rank) :: dimsc = (/1/)      !< chunking length ???
        integer(size_t) :: irbc_size_template   !< size ofiRbc array in iterations logging
        integer(size_t) :: irbc_size            !< size ofiRbc array in iterations logging

        integer :: hdfier     !< error flag for HDF5 library

        cpui = MPI_WTIME()
        cpuo = cpui
#ifdef OPENMP
        nthreads = omp_get_max_threads()
#else
        nthreads = 1
#endif

        if (myid .eq. 0 .and. .not. skip_write) then

            ! Set dataset transfer property to preserve partially initialized fields
            ! during write/read to/from dataset with compound datatype.
            call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdfier)

            call h5pset_preserve_f(plist_id, .TRUE., hdfier)

            maxdims = (/H5S_UNLIMITED_F/)                                                                      ! unlimited array size: "converge" until you get bored
            call h5screate_simple_f(rank, dims, iteration_dspace_id, hdfier, maxdims)             ! Create the dataspace with zero initial size and allow it to grow

            call h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, hdfier)                               ! dataset creation property list with chunking

            call h5pset_chunk_f(crp_list, rank, dimsc, hdfier)

            ! declare "iteration" compound datatype
            ! declare array parts
            call h5tarray_create_f(H5T_NATIVE_DOUBLE, 2, int((/mn, Mvol + 1/), hsize_t), iRZbscArray_id, hdfier)  ! create array datatypes for i{R,Z}b{c,s}
            call h5tget_size_f(iRZbscArray_id, irbc_size, hdfier)
            call h5tget_size_f(H5T_NATIVE_INTEGER, type_size_i, hdfier)                           ! size of an integer field
            call h5tget_size_f(H5T_NATIVE_DOUBLE, type_size_d, hdfier)                           ! size of a   double field
            iteration_dtype_size = 2*type_size_i + 2*type_size_d + 4*irbc_size                                   ! wflag, nDcalls, Energy, ForceErr, i{R,Z}b{c,s}

            call h5tcreate_f(H5T_COMPOUND_F, iteration_dtype_size, iteration_dtype_id, hdfier)     ! create compound datatype

            offset = 0                                                                                           ! offset for first field starts at 0

            call h5tinsert_f(iteration_dtype_id, "nDcalls", offset, H5T_NATIVE_INTEGER, hdfier)    ! insert "nDcalls" field in datatype
            offset = offset + type_size_i                                                                        ! increment offset by size of field

            call h5tinsert_f(iteration_dtype_id, "Energy", offset, H5T_NATIVE_DOUBLE, hdfier)      ! insert "Energy" field in datatype
            offset = offset + type_size_d                                                                        ! increment offset by size of field

            call h5tinsert_f(iteration_dtype_id, "ForceErr", offset, H5T_NATIVE_DOUBLE, hdfier)    ! insert "ForceErr" field in datatype
            offset = offset + type_size_d                                                                        ! increment offset by size of field

            call h5tinsert_f(iteration_dtype_id, "iRbc", offset, iRZbscArray_id, hdfier)          ! insert "iRbc" field in datatype
            offset = offset + irbc_size                                                                          ! increment offset by size of field

            call h5tinsert_f(iteration_dtype_id, "iZbs", offset, iRZbscArray_id, hdfier)          ! insert "iZbs" field in datatype
            offset = offset + irbc_size                                                                          ! increment offset by size of field

            call h5tinsert_f(iteration_dtype_id, "iRbs", offset, iRZbscArray_id, hdfier)           ! insert "iRbs" field in datatype
            offset = offset + irbc_size                                                                          ! increment offset by size of field

            call h5tinsert_f(iteration_dtype_id, "iZbc", offset, iRZbscArray_id, hdfier)           ! insert "iZbc" field in datatype
            offset = offset + irbc_size                                                                          ! increment offset by size of field

            call h5dcreate_f(file_id, "iterations", iteration_dtype_id, iteration_dspace_id, &      ! create dataset with compound type
                            & iteration_dset_id, hdfier, crp_list)

            call h5sclose_f(iteration_dspace_id, hdfier)                                           ! Terminate access to the data space (does not show up in obj_count below)
            ! --> only needed for creation of dataset

            ! Create memory types. We have to create a compound datatype
            ! for each member we want to write.
            call h5tcreate_f(H5T_COMPOUND_F, type_size_i, dt_nDcalls_id, hdfier)
            call h5tcreate_f(H5T_COMPOUND_F, type_size_d, dt_Energy_id, hdfier)
            call h5tcreate_f(H5T_COMPOUND_F, type_size_d, dt_ForceErr_id, hdfier)
            call h5tcreate_f(H5T_COMPOUND_F, irbc_size, dt_iRbc_id, hdfier)
            call h5tcreate_f(H5T_COMPOUND_F, irbc_size, dt_iZbs_id, hdfier)
            call h5tcreate_f(H5T_COMPOUND_F, irbc_size, dt_iRbs_id, hdfier)
            call h5tcreate_f(H5T_COMPOUND_F, irbc_size, dt_iZbc_id, hdfier)

            offset = 0
            call h5tinsert_f(dt_nDcalls_id, "nDcalls", offset, H5T_NATIVE_INTEGER, hdfier)
            call h5tinsert_f(dt_Energy_id, "Energy", offset, H5T_NATIVE_DOUBLE, hdfier)
            call h5tinsert_f(dt_ForceErr_id, "ForceErr", offset, H5T_NATIVE_DOUBLE, hdfier)
            call h5tinsert_f(dt_iRbc_id, "iRbc", offset, iRZbscArray_id, hdfier)
            call h5tinsert_f(dt_iZbs_id, "iZbs", offset, iRZbscArray_id, hdfier)
            call h5tinsert_f(dt_iRbs_id, "iRbs", offset, iRZbscArray_id, hdfier)
            call h5tinsert_f(dt_iZbc_id, "iZbc", offset, iRZbscArray_id, hdfier)

            ! create memspace with size of compound object to append
            dims(1) = 1 ! only append one iteration at a time
            call h5screate_simple_f(rank, dims, memspace, hdfier)

            call h5pclose_f(crp_list, hdfier)
            call h5tclose_f(iteration_dtype_id, hdfier)                                        ! Terminate access to the datatype
            call h5tclose_f(iRZbscArray_id, hdfier)                                            ! Terminate access to the datatype

        end if ! myid.eq.0

    end subroutine init_convergence_output

!> \brief Write convergence output (evolution of interface geometry, force, etc).
!> \ingroup grp_output
!>
    subroutine write_convergence_output(nDcalls, ForceErr)

        use allglobal, only: myid, mn, Mvol, Energy, iRbc, iZbs, iRbs, iZbc, MPI_COMM_SPEC

#ifdef OPENMP
        USE OMP_LIB
#endif
        use mpi
        implicit none
        integer :: ierr, astat, ios, nthreads, ithread
        real(wp) :: cput, cpui, cpuo = 0 ! cpu time; cpu initial; cpu old; 31 Jan 13;

        integer :: hdfier     !< error flag for HDF5 library

        integer, intent(in) :: nDcalls
        real(wp), intent(in) :: ForceErr

        cpui = MPI_WTIME()
        cpuo = cpui
#ifdef OPENMP
        nthreads = omp_get_max_threads()
#else
        nthreads = 1
#endif

        if (myid .eq. 0 .and. .not. skip_write) then

            ! append updated values to "iterations" dataset

            ! open dataspace to get current state of dataset
            call h5dget_space_f(iteration_dset_id, dataspace, hdfier)

            ! get current size of dataset
            call h5sget_simple_extent_dims_f(dataspace, old_data_dims, max_dims, hdfier)

            if (hdfier .ne. 1) then
                write (6, '("sphdf5 :      fatal : myid=",i3," ; hdfier.ne.1 ; rank of convergence dataspace is not 1 ;")') myid
!     call MPI_ABORT(MPI_COMM_SPEC, 1) ! TODO: get this to compile again...
                stop "sphdf5 : hdfier.ne.1 : rank of convergence dataspace is not 1  ;"
            end if

            ! blow up dataset to new size
            data_dims = old_data_dims + 1
            call h5dset_extent_f(iteration_dset_id, data_dims, hdfier)

            ! get dataspace slab corresponding to region which the iterations dataset was extended by
            call h5dget_space_f(iteration_dset_id, dataspace, hdfier)                                     ! re-select dataspace to update size info in HDF5 lib
            call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, old_data_dims, (/INT(1, HSIZE_T)/), hdfier) ! newly appended slab is at old size and 1 long

            ! write next iteration object
            call h5dwrite_f(iteration_dset_id, dt_nDcalls_id, nDcalls, INT((/1/), HSIZE_T), hdfier, &
                            mem_space_id=memspace, file_space_id=dataspace, xfer_prp=plist_id)
            call h5dwrite_f(iteration_dset_id, dt_Energy_id, Energy, INT((/1/), HSIZE_T), hdfier, &
                            mem_space_id=memspace, file_space_id=dataspace, xfer_prp=plist_id)
            call h5dwrite_f(iteration_dset_id, dt_ForceErr_id, ForceErr, INT((/1/), HSIZE_T), hdfier, &
                            mem_space_id=memspace, file_space_id=dataspace, xfer_prp=plist_id)
            call h5dwrite_f(iteration_dset_id, dt_iRbc_id, iRbc, INT((/mn, Mvol + 1/), HSIZE_T), hdfier, &
                            mem_space_id=memspace, file_space_id=dataspace, xfer_prp=plist_id)
            call h5dwrite_f(iteration_dset_id, dt_iZbs_id, iZbs, INT((/mn, Mvol + 1/), HSIZE_T), hdfier, &
                            mem_space_id=memspace, file_space_id=dataspace, xfer_prp=plist_id)
            call h5dwrite_f(iteration_dset_id, dt_iRbs_id, iRbs, INT((/mn, Mvol + 1/), HSIZE_T), hdfier, &
                            mem_space_id=memspace, file_space_id=dataspace, xfer_prp=plist_id)
            call h5dwrite_f(iteration_dset_id, dt_iZbc_id, iZbc, INT((/mn, Mvol + 1/), HSIZE_T), hdfier, &
                            mem_space_id=memspace, file_space_id=dataspace, xfer_prp=plist_id)

            ! dataspace to appended object should be closed now
            ! MAYBE we otherwise keep all the iterations in memory?
            call h5sclose_f(dataspace, hdfier)

        end if ! myid.eq.0

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
        use allglobal, only: myid, ijreal, ijimag, jireal, &
        &                     Nt, Nz, Ntz, Mvol, pi2nfp, ivol, mn, Node, gBzeta, &
        &                     Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, &
        &                     Rij, Zij, sg
        use inputlist, only: Lrad, Igeometry, Nvol, Ngrid, rtor, rpol
        use cputiming, only: Tsphdf5

#ifdef OPENMP
        USE OMP_LIB
#endif
        use mpi
        implicit none
        integer :: ierr, astat, ios, nthreads, ithread
        real(wp) :: cput, cpui, cpuo = 0 ! cpu time; cpu initial; cpu old; 31 Jan 13;

        integer(hid_t) :: grpGrid
        integer :: sumLrad, alongLrad, Ngrid_local, Ngrid_sum
        integer :: vvol, ii, jj, kk, jk, Lcurvature
        real(wp) :: lss, teta, zeta, st(1:Node), Bst(1:Node)
        real(wp), allocatable :: Rij_grid(:, :), Zij_grid(:, :), sg_grid(:, :), ijreal_grid(:, :), ijimag_grid(:, :), jireal_grid(:, :)

        cpui = MPI_WTIME()
        cpuo = cpui
#ifdef OPENMP
        nthreads = omp_get_max_threads()
#else
        nthreads = 1
#endif

        if (myid .eq. 0 .and. .not. skip_write) then

            ijreal(1:Ntz) = zero
            ijimag(1:Ntz) = zero
            jireal(1:Ntz) = zero

            call HDEFGRP(file_id, "grid", grpGrid)

            ! Igeometry already is in input, Mvol already is in output
            call HWRITEIV(grpGrid, "Nt", 1, (/Nt/))
            call HWRITEIV(grpGrid, "Nz", 1, (/Nz/))
            call HWRITEIV(grpGrid, "Ntz", 1, (/Ntz/))
            call HWRITERV(grpGrid, "pi2nfp", 1, (/pi2nfp/))

            ! combine all radial parts into one dimension as Lrad values can be different for different volumes
            if (Ngrid .lt. 0) then
                sumLrad = sum(Lrad(1:Mvol) + 1)
            else
                sumLrad = (Ngrid + 1)*Mvol
            end if

            allocate (Rij_grid(1:sumLrad, 1:Ntz), stat=astat)
            Rij_grid(1:sumLrad, 1:Ntz) = zero

            allocate (Zij_grid(1:sumLrad, 1:Ntz), stat=astat)
            Zij_grid(1:sumLrad, 1:Ntz) = zero

            allocate (sg_grid(1:sumLrad, 1:Ntz), stat=astat)
            sg_grid(1:sumLrad, 1:Ntz) = zero

            allocate (ijreal_grid(1:sumLrad, 1:Ntz), stat=astat)
            ijreal_grid(1:sumLrad, 1:Ntz) = zero

            allocate (ijimag_grid(1:sumLrad, 1:Ntz), stat=astat)
            ijimag_grid(1:sumLrad, 1:Ntz) = zero

            allocate (jireal_grid(1:sumLrad, 1:Ntz), stat=astat)
            jireal_grid(1:sumLrad, 1:Ntz) = zero

            Ngrid_sum = 0

            do vvol = 1, Mvol
                ivol = vvol

                if (Igeometry .eq. 1 .or. vvol .gt. 1) then
                    Lcoordinatesingularity = .false.
                else
                    Lcoordinatesingularity = .true.
                end if

                if (vvol .le. Nvol) then
                    Lplasmaregion = .true.
                else
                    Lplasmaregion = .false.
                end if

                Lvacuumregion = .not. Lplasmaregion
                ! sets Lcoordinatesingularity and Lplasmaregion ;

                if (Ngrid .lt. 0) then
                    Ngrid_local = Lrad(vvol)  ! default
                else
                    Ngrid_local = Ngrid
                end if
                if (Ngrid_local .eq. 0) cycle               ! nothing to output

                do ii = 0, Ngrid_local ! sub-grid;
                    lss = ii*two/Ngrid_local - one
                    if (Lcoordinatesingularity .and. ii .eq. 0) then
                        Lcurvature = 0 ! Jacobian is not defined;
                    else
                        Lcurvature = 1 ! compute Jacobian       ;
                    end if

                    cput = MPI_WTIME()
                    Tsphdf5 = Tsphdf5 + (cput - cpuo)
                    call coords(vvol, lss, Lcurvature, Ntz, mn)
                    cpuo = MPI_WTIME()
                    ! only Rij(0,:) and Zij(0,:) are required; Rmn & Zmn are available;

                    alongLrad = Ngrid_sum + ii + 1

                    Rij_grid(alongLrad, 1:Ntz) = Rij(1:Ntz, 0, 0)
                    Zij_grid(alongLrad, 1:Ntz) = Zij(1:Ntz, 0, 0)
                    sg_grid(alongLrad, 1:Ntz) = sg(1:Ntz, 0)

                    if (Lcurvature .eq. 1) then

                        select case (Igeometry)

                        case (3)
                            do kk = 0, Nz - 1; zeta = kk*pi2nfp/Nz
                                do jj = 0, Nt - 1; teta = jj*pi2/Nt; jk = 1 + jj + kk*Nt; st(1:2) = (/lss, teta/)

                                    cput = MPI_WTIME()
                                    Tsphdf5 = Tsphdf5 + (cput - cpuo)
                                    call bfield(zeta, st(1:Node), Bst(1:Node))
                                    cpuo = MPI_WTIME()

                                    ijreal(jk) = (Rij(jk, 1, 0)*Bst(1) + Rij(jk, 2, 0)*Bst(2) + Rij(jk, 3, 0)*one)*gBzeta/sg(jk, 0) ! BR;
                                    ijimag(jk) = (one)*gBzeta/sg(jk, 0) ! Bp;
                                    jireal(jk) = (Zij(jk, 1, 0)*Bst(1) + Zij(jk, 2, 0)*Bst(2) + Zij(jk, 3, 0)*one)*gBzeta/sg(jk, 0) ! BZ;
                                end do
                            end do

                        case (1)
                            do kk = 0, Nz - 1; zeta = kk*pi2nfp/Nz
                                do jj = 0, Nt - 1; teta = jj*pi2/Nt; jk = 1 + jj + kk*Nt; st(1:2) = (/lss, teta/)

                                    cput = MPI_WTIME()
                                    Tsphdf5 = Tsphdf5 + (cput - cpuo)
                                    call bfield(zeta, st(1:Node), Bst(1:Node))
                                    cpuo = MPI_WTIME()

                                    ijreal(jk) = (Rij(jk, 1, 0)*Bst(1) + Rij(jk, 2, 0)*Bst(2) + Rij(jk, 3, 0)*one)*gBzeta/sg(jk, 0) ! BR;
                                    ijimag(jk) = (rpol)*gBzeta/sg(jk, 0) ! Bzeta;
                                    jireal(jk) = (+rtor*Bst(2))*gBzeta/sg(jk, 0) ! Btheta;
                                end do
                            end do

                        case (2)
                            do kk = 0, Nz - 1; zeta = kk*pi2nfp/Nz
                                do jj = 0, Nt - 1; teta = jj*pi2/Nt; jk = 1 + jj + kk*Nt; st(1:2) = (/lss, teta/)

                                    cput = MPI_WTIME()
                                    Tsphdf5 = Tsphdf5 + (cput - cpuo)
                                    call bfield(zeta, st(1:Node), Bst(1:Node))
                                    cpuo = MPI_WTIME()

                                    ijreal(jk) = (Rij(jk, 1, 0)*Bst(1) + Rij(jk, 2, 0)*Bst(2) + Rij(jk, 3, 0)*one)*gBzeta/sg(jk, 0) ! BR;
                                    ijimag(jk) = (one)*gBzeta/sg(jk, 0) ! Bp;
                                    jireal(jk) = (Bst(2))*gBzeta/sg(jk, 0) ! BZ;
                                end do
                            end do

                        end select !Igeometry
                    end if ! end of if( Lcurvature.eq.1 ) ;

                    ijreal_grid(alongLrad, 1:Ntz) = ijreal(1:Ntz)
                    ijimag_grid(alongLrad, 1:Ntz) = ijimag(1:Ntz)
                    jireal_grid(alongLrad, 1:Ntz) = jireal(1:Ntz)

                end do ! end of do ii;

                Ngrid_sum = Ngrid_sum + Ngrid_local + 1 ! offset for storing data

            end do ! end of do vvol;

            call HWRITERA(grpGrid, "Rij", sumLrad, Ntz, Rij_grid)
            call HWRITERA(grpGrid, "Zij", sumLrad, Ntz, Zij_grid)
            call HWRITERA(grpGrid, "sg", sumLrad, Ntz, sg_grid)
            call HWRITERA(grpGrid, "BR", sumLrad, Ntz, ijreal_grid)
            call HWRITERA(grpGrid, "Bp", sumLrad, Ntz, ijimag_grid)
            call HWRITERA(grpGrid, "BZ", sumLrad, Ntz, jireal_grid)

            deallocate (Rij_grid, stat=astat)
            deallocate (Zij_grid, stat=astat)
            deallocate (sg_grid, stat=astat)
            deallocate (ijreal_grid, stat=astat)
            deallocate (ijimag_grid, stat=astat)
            deallocate (jireal_grid, stat=astat)

            call HCLOSEGRP(grpGrid)

        end if ! myid.eq.0
    end subroutine ! write_grid

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
    subroutine init_flt_output(numTrajTotal)

        use allglobal, only: Nz, Mvol, lmns
        use inputlist, only: nPpts

#ifdef OPENMP
        USE OMP_LIB
#endif
        use mpi
        implicit none
        integer :: ierr, astat, ios, nthreads, ithread
        real(wp) :: cput, cpui, cpuo = 0 ! cpu time; cpu initial; cpu old; 31 Jan 13;

        integer, intent(in) :: numTrajTotal ! total number of trajectories
        integer(HSIZE_T), dimension(rankP) :: dims_traj   ! Dataset dimensions.
        integer(HSIZE_T), dimension(rankP) :: length      ! Dataset dimensions.

        integer :: hdfier     !< error flag for HDF5 library

        cpui = MPI_WTIME()
        cpuo = cpui
#ifdef OPENMP
        nthreads = omp_get_max_threads()
#else
        nthreads = 1
#endif

        if (myid .eq. 0 .and. .not. skip_write) then

            ! create Poincare group in HDF5 file
            call HDEFGRP(file_id, "poincare", grpPoincare)

            dims_traj = (/Nz, nPpts, numTrajTotal/) ! dimensions for whole Poincare dataset
            length = (/Nz, nPpts, 1/) ! which is written in these slice lengths

            ! Create the data space for the  dataset.
            call h5screate_simple_f(rankP, dims_traj, filespace_t, hdfier)
            call h5screate_simple_f(rankP, dims_traj, filespace_s, hdfier)
            call h5screate_simple_f(rankP, dims_traj, filespace_R, hdfier)
            call h5screate_simple_f(rankP, dims_traj, filespace_Z, hdfier)
            call h5screate_simple_f(1, int((/numTrajTotal/), HSIZE_T), filespace_success, hdfier)

            ! Create the dataset with default properties.
            call h5dcreate_f(grpPoincare, "t", H5T_NATIVE_DOUBLE, filespace_t, dset_id_t, hdfier)
            call h5dcreate_f(grpPoincare, "s", H5T_NATIVE_DOUBLE, filespace_s, dset_id_s, hdfier)
            call h5dcreate_f(grpPoincare, "R", H5T_NATIVE_DOUBLE, filespace_R, dset_id_R, hdfier)
            call h5dcreate_f(grpPoincare, "Z", H5T_NATIVE_DOUBLE, filespace_Z, dset_id_Z, hdfier)
            call h5dcreate_f(grpPoincare, "success", H5T_NATIVE_INTEGER, filespace_success, dset_id_success, hdfier)

            ! filespaces can be closed as soon as datasets are created
            call h5sclose_f(filespace_t, hdfier)
            call h5sclose_f(filespace_s, hdfier)
            call h5sclose_f(filespace_R, hdfier)
            call h5sclose_f(filespace_Z, hdfier)
            call h5sclose_f(filespace_success, hdfier)

            ! Select hyperslab in the file.
            call h5dget_space_f(dset_id_t, filespace_t, hdfier)
            call h5dget_space_f(dset_id_s, filespace_s, hdfier)
            call h5dget_space_f(dset_id_R, filespace_R, hdfier)
            call h5dget_space_f(dset_id_Z, filespace_Z, hdfier)
            call h5dget_space_f(dset_id_success, filespace_success, hdfier)

            ! Each process defines dataset in memory and writes it to the hyperslab in the file.
            call h5screate_simple_f(rankP, length, memspace_t, hdfier)
            call h5screate_simple_f(rankP, length, memspace_s, hdfier)
            call h5screate_simple_f(rankP, length, memspace_R, hdfier)
            call h5screate_simple_f(rankP, length, memspace_Z, hdfier)
            call h5screate_simple_f(1, int((/1/), HSIZE_T), memspace_success, hdfier)

            ! create rotational transform group in HDF5 file
            call HDEFGRP(file_id, "transform", grpTransform)

            ! Create the data space for the  dataset.
            call h5screate_simple_f(rankT, int((/2, Mvol/), HSIZE_T), filespace_diotadxup, hdfier)
            call h5screate_simple_f(rankT, int((/numTrajTotal, 2/), HSIZE_T), filespace_fiota, hdfier)

            ! Create the dataset with default properties.
            call h5dcreate_f(grpTransform, "diotadxup", H5T_NATIVE_DOUBLE, filespace_diotadxup, dset_id_diotadxup, hdfier)
            call h5dcreate_f(grpTransform, "fiota", H5T_NATIVE_DOUBLE, filespace_fiota, dset_id_fiota, hdfier)

            ! filespaces can be closed as soon as datasets are created
            call h5sclose_f(filespace_diotadxup, hdfier)
            call h5sclose_f(filespace_fiota, hdfier)

            ! Select hyperslab in the file.
            call h5dget_space_f(dset_id_diotadxup, filespace_diotadxup, hdfier)
            call h5dget_space_f(dset_id_fiota, filespace_fiota, hdfier)

            ! Each process defines dataset in memory and writes it to the hyperslab in the file.
            call h5screate_simple_f(rankT, int((/2, 1/), HSIZE_T), memspace_diotadxup, hdfier)

        end if ! myid.eq.0

    end subroutine init_flt_output

!> \brief Write a hyperslab of Poincare data corresponding to the output of one parallel worker.
!> \ingroup grp_output
!>
!> @param offset radial offset at which the data belongs
!> @param data output from field-line tracing
!> @param success flags to indicate if integrator was successful
    subroutine write_poincare(offset, data, success)

        use allglobal, only: Nz
        use inputlist, only: nPpts

#ifdef OPENMP
        USE OMP_LIB
#endif
        use mpi
        implicit none
        integer :: ierr, astat, ios, nthreads, ithread
        real(wp) :: cput, cpui, cpuo = 0 ! cpu time; cpu initial; cpu old; 31 Jan 13;
        integer :: hdfier     !< error flag for HDF5 library

        integer, intent(in) :: offset, success(:)
        real(wp), intent(in) :: data(:, :, :)
        integer(hsize_t), dimension(3) :: length
        integer(HSIZE_T), dimension(2) :: dims_singleTraj ! dimensions of single trajectory data

        cpui = MPI_WTIME()
        cpuo = cpui
#ifdef OPENMP
        nthreads = omp_get_max_threads()
#else
        nthreads = 1
#endif

        if (myid .eq. 0 .and. .not. skip_write) then

            dims_singleTraj = (/Nz, nPpts/)
            length = (/Nz, nPpts, 1/)

            ! On entry, Fortran does not know that indexing in data is from 0 to Nz-1.
            ! Hence, use default indices 1:Nz in this routine
            call h5sselect_hyperslab_f(filespace_t, H5S_SELECT_SET_F, int((/0, 0, offset/), HSSIZE_T), length, hdfier)
            call h5dwrite_f(dset_id_t, H5T_NATIVE_DOUBLE, data(1, 1:Nz, 1:nPpts), dims_singleTraj, hdfier, &
                            file_space_id=filespace_t, mem_space_id=memspace_t)

            call h5sselect_hyperslab_f(filespace_s, H5S_SELECT_SET_F, int((/0, 0, offset/), HSSIZE_T), length, hdfier)
            call h5dwrite_f(dset_id_s, H5T_NATIVE_DOUBLE, data(2, 1:Nz, 1:nPpts), dims_singleTraj, hdfier, &
                            file_space_id=filespace_s, mem_space_id=memspace_s)

            call h5sselect_hyperslab_f(filespace_R, H5S_SELECT_SET_F, int((/0, 0, offset/), HSSIZE_T), length, hdfier)
            call h5dwrite_f(dset_id_R, H5T_NATIVE_DOUBLE, data(3, 1:Nz, 1:nPpts), dims_singleTraj, hdfier, &
                            file_space_id=filespace_R, mem_space_id=memspace_R)

            call h5sselect_hyperslab_f(filespace_Z, H5S_SELECT_SET_F, int((/0, 0, offset/), HSSIZE_T), length, hdfier)
            call h5dwrite_f(dset_id_Z, H5T_NATIVE_DOUBLE, data(4, 1:Nz, 1:nPpts), dims_singleTraj, hdfier, &
                            file_space_id=filespace_Z, mem_space_id=memspace_Z)

            call h5sselect_hyperslab_f(filespace_success, H5S_SELECT_SET_F, int((/offset/), HSSIZE_T), int((/1/), HSIZE_T), hdfier)
            call h5dwrite_f(dset_id_success, H5T_NATIVE_INTEGER, success, int((/1/), HSIZE_T), hdfier, &
                            file_space_id=filespace_success, mem_space_id=memspace_success)

        end if ! myid.eq.0

    end subroutine write_poincare

!> \brief Write the rotational transform output from field line following.
!> \ingroup grp_output
!>
!> @param offset radial offset at which the data belongs
!> @param length length of dataset to write
!> @param lvol nested volume index
!> @param diotadxup derivative of rotational transform (?)
!> @param fiota rotational transform
    subroutine write_transform(offset, length, lvol, diotadxup, fiota)

#ifdef OPENMP
        USE OMP_LIB
#endif
        use mpi
        implicit none
        integer :: ierr, astat, ios, nthreads, ithread
        real(wp) :: cput, cpui, cpuo = 0 ! cpu time; cpu initial; cpu old; 31 Jan 13;

        integer, intent(in) :: offset, length, lvol
        real(wp), intent(in) :: diotadxup(:), fiota(:, :)
        integer :: hdfier     !< error flag for HDF5 library

        cpui = MPI_WTIME()
        cpuo = cpui
#ifdef OPENMP
        nthreads = omp_get_max_threads()
#else
        nthreads = 1
#endif

        if (myid .eq. 0 .and. .not. skip_write) then

            call h5sselect_hyperslab_f(filespace_diotadxup, H5S_SELECT_SET_F, int((/0, lvol - 1/), HSSIZE_T), int((/2, 1/), HSSIZE_T), hdfier)
            call h5dwrite_f(dset_id_diotadxup, H5T_NATIVE_DOUBLE, diotadxup, int((/2, 1/), HSSIZE_T), hdfier, &
                            file_space_id=filespace_diotadxup, mem_space_id=memspace_diotadxup)

            ! length of fiota piece to write here may change, so open and close memspace each time a new hyperslab is written
            call h5screate_simple_f(rankT, int((/length, 2/), HSIZE_T), memspace_fiota, hdfier)

            call h5sselect_hyperslab_f(filespace_fiota, H5S_SELECT_SET_F, int((/offset, 0/), HSSIZE_T), int((/length, 2/), HSSIZE_T), hdfier)
            call h5dwrite_f(dset_id_fiota, H5T_NATIVE_DOUBLE, fiota(1:length, 1:2), int((/length, 2/), HSSIZE_T), hdfier, &
                            file_space_id=filespace_fiota, mem_space_id=memspace_fiota)

            call h5sclose_f(memspace_fiota, hdfier)

        end if ! myid.eq.0

    end subroutine write_transform

!> \brief Finalize Poincare output.
!> \ingroup grp_output
!>
!> This closes the still-open datasets related to field-line tracing,
!> which had to be kept open during the tracing to be able to write
!> the outputs directly when a given worker thread is finished.
    subroutine finalize_flt_output
        use allglobal, only: myid, skip_write
        integer :: hdfier     !< error flag for HDF5 library

        if (myid .eq. 0 .and. .not. skip_write) then

            ! close filespaces
            call h5sclose_f(filespace_t, hdfier)
            call h5sclose_f(filespace_s, hdfier)
            call h5sclose_f(filespace_R, hdfier)
            call h5sclose_f(filespace_Z, hdfier)
            call h5sclose_f(filespace_success, hdfier)
            call h5sclose_f(filespace_diotadxup, hdfier)
            call h5sclose_f(filespace_fiota, hdfier)

            ! close dataspaces
            call h5sclose_f(memspace_t, hdfier)
            call h5sclose_f(memspace_s, hdfier)
            call h5sclose_f(memspace_R, hdfier)
            call h5sclose_f(memspace_Z, hdfier)
            call h5sclose_f(memspace_success, hdfier)
            call h5sclose_f(memspace_diotadxup, hdfier)
            ! memspace_fiota is re-opened/closed in each iteration (see write_transform)

            ! close datasets
            call h5dclose_f(dset_id_t, hdfier)
            call h5dclose_f(dset_id_s, hdfier)
            call h5dclose_f(dset_id_R, hdfier)
            call h5dclose_f(dset_id_Z, hdfier)
            call h5dclose_f(dset_id_success, hdfier)
            call h5dclose_f(dset_id_diotadxup, hdfier)
            call h5dclose_f(dset_id_fiota, hdfier)

            ! close groups
            call HCLOSEGRP(grpPoincare)
            call HCLOSEGRP(grpTransform)

        end if ! myid.eq.0

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
        use allglobal, only: mn, myid, skip_write

        integer, intent(in) :: sumLrad
        real(wp), dimension(:, :), intent(in) :: allAte
        real(wp), dimension(:, :), intent(in) :: allAze
        real(wp), dimension(:, :), intent(in) :: allAto
        real(wp), dimension(:, :), intent(in) :: allAzo

        integer(hid_t) :: grpVectorPotential

        if (myid .eq. 0 .and. .not. skip_write) then

            call HDEFGRP(file_id, "vector_potential", grpVectorPotential)

            call HWRITERA(grpVectorPotential, "Ate", sumLrad, mn, allAte(1:sumLrad, 1:mn))
            call HWRITERA(grpVectorPotential, "Aze", sumLrad, mn, allAze(1:sumLrad, 1:mn))
            call HWRITERA(grpVectorPotential, "Ato", sumLrad, mn, allAto(1:sumLrad, 1:mn))
            call HWRITERA(grpVectorPotential, "Azo", sumLrad, mn, allAzo(1:sumLrad, 1:mn))

            call HCLOSEGRP(grpVectorPotential)

        end if ! myid.eq.0

    end subroutine write_vector_potential

!> \brief Write the final state of the equilibrium to the output file.
!> \ingroup grp_output
!>
    subroutine hdfint

        use fileunits, only: ounit
        use inputlist
        use allglobal, only: ncpu, cpus, &
                             Mvol, ForceErr, &
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

#ifdef OPENMP
        USE OMP_LIB
#endif
        use mpi
        implicit none
        integer :: ierr, astat, ios, nthreads, ithread
        real(wp) :: cput, cpui, cpuo = 0 ! cpu time; cpu initial; cpu old; 31 Jan 13;

        integer :: Mrad
        real(wp) :: tvolume

        integer(hid_t) :: grpOutput

        cpui = MPI_WTIME()
        cpuo = cpui
#ifdef OPENMP
        nthreads = omp_get_max_threads()
#else
        nthreads = 1
#endif

        if (myid .eq. 0 .and. .not. skip_write) then

            call HDEFGRP(file_id, "output", grpOutput)

            call HWRITERV(grpOutput, "Vns", mn, iVns(1:mn)) !     stellarator symmetric normal field at boundary; vacuum component;
            call HWRITERV(grpOutput, "Bns", mn, iBns(1:mn)) !     stellarator symmetric normal field at boundary; plasma component;
            call HWRITERV(grpOutput, "Vnc", mn, iVnc(1:mn)) ! non-stellarator symmetric normal field at boundary; vacuum component;
            call HWRITERV(grpOutput, "Bnc", mn, iBnc(1:mn)) ! non-stellarator symmetric normal field at boundary; plasma component;

!> <ul>
!> <li> In addition to the input variables, which are described in global(), the following quantities are written to \c ext.sp.h5 :
!latex
!latex \begin{tabular}{|l|l|l|} \hline

!latex \type{variable}               & type    & \pb{description}  \hline

!latex \type{mn}                     & integer & \pb{number of Fourier modes}
            call HWRITEIV(grpOutput, "mn", 1, (/mn/))
!latex \type{im(1:mn)}               & integer & \pb{poloidal mode numbers}
            call HWRITEIV(grpOutput, "im", mn, im(1:mn))
!latex \type{in(1:mn)}               & integer & \pb{toroidal mode numbers}
            call HWRITEIV(grpOutput, "in", mn, in(1:mn))
!latex \type{mns}                     & integer & \pb{number of Fourier modes}
            call HWRITEIV(grpOutput, "mns", 1, (/mns/))
!latex \type{ims(1:mns)}               & integer & \pb{poloidal mode numbers}
            call HWRITEIV(grpOutput, "ims", mns, ims(1:mns))
!latex \type{ins(1:mns)}               & integer & \pb{toroidal mode numbers}
            call HWRITEIV(grpOutput, "ins", mns, ins(1:mns))
!latex \type{Mvol}                   & integer & \pb{number of interfaces = number of volumes}
            call HWRITEIV(grpOutput, "Mvol", 1, (/Mvol/))
!latex \type{iRbc(1:mn,0:Mvol)}      & real    & \pb{Fourier harmonics, $R_{m,n}$, of interfaces}
            call HWRITERA(grpOutput, "Rbc", mn, (Mvol + 1), iRbc(1:mn, 0:Mvol))
!latex \type{iZbs(1:mn,0:Mvol)}      & real    & \pb{Fourier harmonics, $Z_{m,n}$, of interfaces}
            call HWRITERA(grpOutput, "Zbs", mn, (Mvol + 1), iZbs(1:mn, 0:Mvol))
!latex \type{iRbs(1:mn,0:Mvol)}      & real    & \pb{Fourier harmonics, $R_{m,n}$, of interfaces}
            call HWRITERA(grpOutput, "Rbs", mn, (Mvol + 1), iRbs(1:mn, 0:Mvol))
!latex \type{iZbc(1:mn,0:Mvol)}      & real    & \pb{Fourier harmonics, $Z_{m,n}$, of interfaces}
            call HWRITERA(grpOutput, "Zbc", mn, (Mvol + 1), iZbc(1:mn, 0:Mvol))
!l tex \type{forcetol}               & real    & \pb{force-balance error across interfaces}
!  HWRITERV( grpOutput, 1, forcetol, (/ forcetol /)) ! already in /input/global
!latex \type{ForceErr}               & real    & \pb{force-balance error across interfaces}
            call HWRITERV(grpOutput, "ForceErr", 1, (/ForceErr/))
!latex \type{Ivolume}                & real    & \pb{Volume current at output (parallel, externally induced)}
            call HWRITERV(grpOutput, "Ivolume", Mvol, Ivolume(1:Mvol))
!latex \type{IPDt}                   & real    & \pb{Surface current at output}
            call HWRITERV(grpOutput, "IPDt", Mvol, IPDt(1:Mvol))

            ! the following quantites can be different from input value
            call HWRITERV(grpOutput, "adiabatic", Mvol, adiabatic(1:Nvol))
            call HWRITERV(grpOutput, "helicity", Nvol, helicity(1:Nvol))
            call HWRITERV(grpOutput, "mu", Mvol, mu(1:Mvol))
            call HWRITERV(grpOutput, "tflux", Mvol, tflux(1:Mvol))
            call HWRITERV(grpOutput, "pflux", Mvol, pflux(1:Mvol))

            if (Lcheck .eq. 1) then
!latex \type{beltramierror}          & real    & \pb{error in beltrami field (volume integral)}
                call HWRITERA(grpOutput, "beltramierror", Mvol, 3, beltramierror(1:Mvol, 1:3))
            end if

            if (allocated(vvolume)) then ! why is it required to confirm that vvolume has been allocated ; 24 Nov 16;

                tvolume = sum(vvolume(1:Nvol))
!latex \type{volume}                 & real    & \pb{total volume = $\sum V_v$}
                call HWRITERV(grpOutput, "volume", 1, (/tvolume/))

            else

                if (Wsphdf5) write (ounit, '("hdfint : ", 10x ," : myid=",i3," ; vvolume is not allocated ;")') myid

            end if ! end of if( allocated(vvolume) ) ; 11 Aug 14;

            Mrad = maxval(Lrad(1:Mvol))
!latex \type{Mrad}                   & integer & \pb{the maximum radial (Chebyshev) resolution}
            call HWRITEIV(grpOutput, "Mrad", 1, (/Mrad/))
!latex \type{TT(0:Mrad,0:1,0:1)}     & real    & \pb{the Chebyshev polynomials, $T_l$, and their derivatives, evaluated at $s=\pm 1$}
            call HWRITERC(grpOutput, "TT", (Mrad + 1), 2, 2, TT(0:Mrad, 0:1, 0:1))
!latex \type{Btemn(1:mn,0:1,1:Mvol)} & real    & \pb{the cosine harmonics of the covariant poloidal field,
!latex                                           i.e. $[[B_{\t,j}]]$ evaluated on the inner and outer interface in each volume}
            call HWRITERC(grpOutput, "Btemn", mn, 2, Mvol, Btemn(1:mn, 0:1, 1:Mvol))
!latex \type{Bzemn(1:mn,0:1,1:Mvol)} & real    & \pb{the cosine harmonics of the covariant toroidal field,
!latex                                           i.e. $[[B_{\z,j}]]$ evaluated on the inner and outer interface in each volume}
            call HWRITERC(grpOutput, "Bzemn", mn, 2, Mvol, Bzemn(1:mn, 0:1, 1:Mvol))
!latex \type{Btomn(1:mn,0:1,1:Mvol)} & real    & \pb{the sine harmonics of the covariant poloidal field,
!latex                                           i.e. $[[B_{\t,j}]]$ evaluated on the inner and outer interface in each volume}
            call HWRITERC(grpOutput, "Btomn", mn, 2, Mvol, Btomn(1:mn, 0:1, 1:Mvol))
!latex \type{Bzomn(1:mn,0:1,1:Mvol)} & real    & \pb{the sine harmonics of the covariant toroidal field,
!latex                                           i.e. $[[B_{\z,j}]]$ evaluated on the inner and outer interface in each volume}
            call HWRITERC(grpOutput, "Bzomn", mn, 2, Mvol, Bzomn(1:mn, 0:1, 1:Mvol))

! Write lambda_mn, Fourier harmonics or transformation to straight field line coordinates.
            call HWRITERC(grpOutput, "lambdamn", lmns, Mvol, 2, dlambdaout(1:lmns, 1:Mvol, 0:1))

            if (Lperturbed .eq. 1) then

!latex \type{dRbc(1:mn,0:Nvol)}      & real    & \pb{Fourier harmonics, $R_{j}$, of interfaces; linearly perturbed solution}
                call HWRITERA(grpOutput, "dRbc", mn, (Nvol + 1), dRbc(1:mn, 0:Nvol))
!latex \type{dZbs(1:mn,0:Nvol)}      & real    & \pb{Fourier harmonics, $Z_{j}$, of interfaces; linearly perturbed solution}
                call HWRITERA(grpOutput, "dZbs", mn, (Nvol + 1), dZbs(1:mn, 0:Nvol))
!latex \type{dRbs(1:mn,0:Nvol)}      & real    & \pb{Fourier harmonics, $R_{j}$, of interfaces; linearly perturbed solution}
                call HWRITERA(grpOutput, "dRbs", mn, (Nvol + 1), dRbs(1:mn, 0:Nvol))
!latex \type{dZbc(1:mn,0:Nvol)}      & real    & \pb{Fourier harmonics, $Z_{j}$, of interfaces; linearly perturbed solution}
                call HWRITERA(grpOutput, "dZbc", mn, (Nvol + 1), dZbc(1:mn, 0:Nvol))

            end if

!latex \type{lmns}                   & integer & \pb{resolution of straight fieldline transformation}
            call HWRITEIV(grpOutput, "lmns", 1, (/lmns/))

!latex \hline \end{tabular}
!> </li>
!> <li> All quantities marked as real should be treated as double precision. </li>
!> </ul>

            call HCLOSEGRP(grpOutput)

        end if ! myid.eq.0

    end subroutine hdfint

!> \brief Close all open HDF5 objects (we know of) and list any remaining still-open objects.
!> \ingroup grp_output
!>
    subroutine finish_outfile
! Close all open HDF5 objects (we know of) and list any remaining still-open objects
! The goal should be to close all objects specifically!

#ifdef OPENMP
        USE OMP_LIB
#endif
        use mpi
        implicit none
        integer :: ierr, astat, ios, nthreads, ithread
        real(wp) :: cput, cpui, cpuo = 0 ! cpu time; cpu initial; cpu old; 31 Jan 13;

        integer(size_t) :: obj_count                                  ! number of open HDF5 objects
        integer(size_t) :: num_objs    ! number of still-open objects
        integer(hid_t), dimension(:), allocatable :: obj_ids    ! still-open objects
        integer :: iObj
        integer(size_t) :: openLength
        character(len=:), allocatable :: openName
        integer(size_t), parameter :: dummySize = 1
        character(len=dummySize + 1) :: dummyName
        integer :: typeClass
        integer :: hdfier     !< error flag for HDF5 library

        cpui = MPI_WTIME()
        cpuo = cpui
#ifdef OPENMP
        nthreads = omp_get_max_threads()
#else
        nthreads = 1
#endif

        if (myid .eq. 0 .and. .not. skip_write) then

            ! close objects related to convergence output
            call h5tclose_f(dt_nDcalls_id, hdfier)
            call h5tclose_f(dt_Energy_id, hdfier)
            call h5tclose_f(dt_ForceErr_id, hdfier)
            call h5tclose_f(dt_iRbc_id, hdfier)
            call h5tclose_f(dt_iZbs_id, hdfier)
            call h5tclose_f(dt_iRbs_id, hdfier)
            call h5tclose_f(dt_iZbc_id, hdfier)
            call h5dclose_f(iteration_dset_id, hdfier) ! End access to the dataset and release resources used by it.
            call h5pclose_f(plist_id, hdfier)          ! close plist used for 'preserve' flag (does not show up in obj_count below)

            ! check whether we forgot to close some resources; only check for group, dataset and datatype (there is only one file and that should be still open...)
            call h5fget_obj_count_f(file_id, ior(H5F_OBJ_GROUP_F, ior(H5F_OBJ_DATASET_F, H5F_OBJ_DATATYPE_F)), obj_count, hdfier)

            if (obj_count .gt. 0) then
                write (*, '("There are still ",i3," hdf5 objects open")') obj_count
                allocate (obj_ids(1:obj_count))

                ! groups
                call h5fget_obj_ids_f(file_id, H5F_OBJ_GROUP_F, obj_count, obj_ids, hdfier, num_objs) ! get for open objects
                if (num_objs .gt. 0) then
                    write (*, '("There are still ",i3," HDF5 groups open:")') num_objs
                    do iObj = 1, num_objs
                        openLength = 0
                        call h5iget_name_f(obj_ids(iObj), dummyName, dummySize, openLength, hdfier)
                        allocate (character(len=openLength + 1) :: openName)
                        call h5iget_name_f(obj_ids(iObj), openName, openLength, openLength, hdfier)
                        write (*, *) openName
                        deallocate (openName)

                        call h5gclose_f(obj_ids(iObj), hdfier)
                    end do
                end if

                ! datasets
                call h5fget_obj_ids_f(file_id, H5F_OBJ_DATASET_F, obj_count, obj_ids, hdfier, num_objs) ! get for open objects
                if (num_objs .gt. 0) then
                    write (*, '("There are still ",i3," HDF5 datasets open:")') num_objs
                    do iObj = 1, num_objs
                        openLength = 0
                        call h5iget_name_f(obj_ids(iObj), dummyName, dummySize, openLength, hdfier)
                        allocate (character(len=openLength + 1) :: openName)
                        call h5iget_name_f(obj_ids(iObj), openName, openLength, openLength, hdfier)
                        write (*, *) openName(1:openLength)
                        deallocate (openName)

                        call h5dclose_f(obj_ids(iObj), hdfier)
                    end do
                end if

                ! datatypes
                call h5fget_obj_ids_f(file_id, H5F_OBJ_DATATYPE_F, obj_count, obj_ids, hdfier, num_objs) ! get for open objects
                if (num_objs .gt. 0) then
                    write (*, '("There are still ",i3," HDF5 datatypes open:")') num_objs
                    do iObj = 1, num_objs
                        call h5tget_class_f(obj_ids(iObj), typeClass, hdfier) ! determine class of open datatype
                        if (typeClass .eq. H5T_NO_CLASS_F) then; write (*, *) "H5T_NO_CLASS_F"
                        else if (typeClass .eq. H5T_INTEGER_F) then; write (*, *) "H5T_INTEGER_F"
                        else if (typeClass .eq. H5T_FLOAT_F) then; write (*, *) "H5T_FLOAT_F"
                        else if (typeClass .eq. H5T_STRING_F) then; write (*, *) "H5T_STRING_F"
                        else if (typeClass .eq. H5T_BITFIELD_F) then; write (*, *) "H5T_BITFIELD_F"
                        else if (typeClass .eq. H5T_OPAQUE_F) then; write (*, *) "H5T_OPAQUE_F"
                        else if (typeClass .eq. H5T_COMPOUND_F) then; write (*, *) "H5T_COMPOUND_F"
                        else if (typeClass .eq. H5T_REFERENCE_F) then; write (*, *) "H5T_REFERENCE_F"
                        else if (typeClass .eq. H5T_ENUM_F) then; write (*, *) "H5T_ENUM_F"
                        else if (typeClass .eq. H5T_VLEN_F) then; write (*, *) "H5T_VLEN_F"
                        else if (typeClass .eq. H5T_ARRAY_F) then; write (*, *) "H5T_ARRAY_F"
                        else; write (*, *) "UNKNOWN TYPE!"
                        end if

                        call h5tclose_f(obj_ids(iObj), hdfier)
                    end do
                end if

                deallocate (obj_ids)
            end if ! (obj_count.gt.0)

            call h5fclose_f(file_id, hdfier) ! terminate access on output file;
            call h5close_f(hdfier) ! close Fortran interface to the HDF5 library;

        end if ! myid.eq.0

    end subroutine finish_outfile

end module sphdf5
