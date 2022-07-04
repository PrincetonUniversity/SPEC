module h5utils
    use hdf5
    use mod_kinds, only: wp => dp
    use allglobal, only: MPI_COMM_SPEC, myid
    implicit none

    logical, parameter :: hdfDebug = .false.  !< global flag to enable verbal diarrhea commenting HDF5 operations
    integer, parameter :: internalHdf5Msg = 0 !< 1: print internal HDF5 error messages; 0: only error messages from sphdf5

    character(LEN=*), parameter :: aname = "description" !< Attribute name for descriptive info

! private
! public internalHdf5Msg, &
!        HDEFGRP, &
!        HCLOSEGRP, &
!        H5DESCR, &
!        H5DESCR_CDSET, &
!        HWRITELV, &
!        HWRITELV_LO, &
!        HWRITEIV, &
!        HWRITEIV_LO, &
!        HWRITERV, &
!        HWRITERV_LO, &
!        HWRITERA, &
!        HWRITERA_LO, &
!        HWRITERC, &
!        HWRITERC_LO

contains

!> Define a HDF5 group or opens it if it already exists.
!>
!> @param file_id file in which to define the group
!> @param name    name of the new group
!> @param group_id id of the newly-created group
    subroutine HDEFGRP(file_id, name, group_id)
        integer(hid_t), intent(in) :: file_id
        character(len=*), intent(in) :: name
        integer(hid_t), intent(out) :: group_id

        logical :: grp_exists !< flags used to signal if a group already exists
        integer :: hdfier     !< error flag for HDF5 library

        call h5lexists_f(file_id, name, grp_exists, hdfier)

        if (grp_exists) then
            ! if the group already exists, open it
            call h5gopen_f(file_id, name, group_id, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5gopen_f from hdefgrp")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5gopen_f from hdefgrp"
            end if
        else
            ! if group does not exist, create it
            call h5gcreate_f(file_id, name, group_id, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5gcreate_f from hdefgrp")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5gcreate_f from hdefgrp"
            end if
        end if
    end subroutine ! HDEFGRP

!> Close a HDF5 group.
!>
!> @param group_id HDF5 group to close
    subroutine HCLOSEGRP(group_id)
        integer(hid_t), intent(in) :: group_id

        integer :: hdfier !< error flag for HDF5 library

        call h5gclose_f(group_id, hdfier)
        if (hdfier .ne. 0) then
            write (6, '("sphdf5 : "10x" : error calling h5gclose_f from hclosegrp")')
            call MPI_ABORT(MPI_COMM_SPEC, 1)
            stop "sphdf5 : error calling h5gclose_f from hclosegrp"
        end if
    end subroutine ! HCLOSEGRP

!> Describe an already-open HDF5 object identified by group_id
!> with text given in description and leave it open
!>
!> @param group_id ID of group to describe
!> @param description description for group
    subroutine H5DESCR(group_id, description)
        integer(hid_t), intent(in) :: group_id
        character(len=*), intent(in) :: description

        integer :: hdfier !< error flag for HDF5 library

        integer(HID_T) :: attr_id       !< Attribute identifier
        integer(HID_T) :: aspace_id     !< Attribute Dataspace identifier
        integer(HID_T) :: atype_id      !< Attribute Datatype identifier

        integer, parameter :: arank = 1   !< Attribure rank
        integer(HSIZE_T), dimension(arank) :: adims = (/1/) !< Attribute dimension
        integer(SIZE_T) :: attr_len !< Length of the attribute string

        attr_len = len(description)

        !> Create scalar data space for the attribute.
        call h5screate_simple_f(arank, adims, aspace_id, hdfier)
        if (hdfier .ne. 0) then
            write (6, '("sphdf5 : "10x" : error calling h5screate_simple_f from h5descr")')
            call MPI_ABORT(MPI_COMM_SPEC, 1)
            stop "sphdf5 : error calling h5screate_simple_f from h5descr"
        end if

        !> Create datatype for the attribute.
        call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, hdfier)
        if (hdfier .ne. 0) then
            write (6, '("sphdf5 : "10x" : error calling h5tcopy_f from h5descr")')
            call MPI_ABORT(MPI_COMM_SPEC, 1)
            stop "sphdf5 : error calling h5tcopy_f from h5descr"
        end if

        ! set size of datatype for attribute
        call h5tset_size_f(atype_id, attr_len, hdfier)
        if (hdfier .ne. 0) then
            write (6, '("sphdf5 : "10x" : error calling h5tset_size_f from h5descr")')
            call MPI_ABORT(MPI_COMM_SPEC, 1)
            stop "sphdf5 : error calling h5tset_size_f from h5descr"
        end if

        !> create descriptive attribute
        call h5acreate_f(group_id, aname, atype_id, aspace_id, attr_id, hdfier)
        if (hdfier .ne. 0) then
            write (6, '("sphdf5 : "10x" : error calling h5acreate_f from h5descr")')
            call MPI_ABORT(MPI_COMM_SPEC, 1)
            stop "sphdf5 : error calling h5acreate_f from h5descr"
        end if

        !> Write the attribute data.
        call h5awrite_f(attr_id, atype_id, description, adims, hdfier)
        if (hdfier .ne. 0) then
            write (6, '("sphdf5 : "10x" : error calling h5awrite_f from h5descr")')
            call MPI_ABORT(MPI_COMM_SPEC, 1)
            stop "sphdf5 : error calling h5awrite_f from h5descr"
        end if

        !> Close the attribute.
        call h5aclose_f(attr_id, hdfier)
        if (hdfier .ne. 0) then
            write (6, '("sphdf5 : "10x" : error calling h5aclose_f from h5descr")')
            call MPI_ABORT(MPI_COMM_SPEC, 1)
            stop "sphdf5 : error calling h5aclose_f from h5descr"
        end if

        !> Close the attribute datatype.
        call h5tclose_f(atype_id, hdfier)
        if (hdfier .ne. 0) then
            write (6, '("sphdf5 : "10x" : error calling h5tclose_f from h5descr")')
            call MPI_ABORT(MPI_COMM_SPEC, 1)
            stop "sphdf5 : error calling h5tclose_f from h5descr"
        end if

        !> Terminate access to the data space.
        call h5sclose_f(aspace_id, hdfier)
        if (hdfier .ne. 0) then
            write (6, '("sphdf5 : "10x" : error calling h5sclose_f from h5descr")')
            call MPI_ABORT(MPI_COMM_SPEC, 1)
            stop "sphdf5 : error calling h5sclose_f from h5descr"
        end if
    end subroutine ! H5DESCR

!> Describe an already-open HDF5 dataset identified by dset_id
!> with text given in _2 and close it at the end
!>
!> @param dset_id     ID of dataset to describe
!> @param description description for dataset
    subroutine H5DESCR_CDSET(dset_id, description)
        integer(hid_t), intent(in) :: dset_id
        character(len=*), intent(in) :: description

        integer :: hdfier !< error flag for HDF5 library

        integer(HID_T) :: aspace_id     !< Attribute Dataspace identifier
        integer(HID_T) :: atype_id      !< Attribute Datatype identifier
        integer(HID_T) :: attr_id       !< Attribute identifier

        integer(SIZE_T) :: attr_len !< Length of the attribute string
        integer, parameter :: arank = 1   !< Attribure rank
        integer(HSIZE_T), dimension(arank) :: adims = (/1/) !< Attribute dimension

        attr_len = len(description)

        ! Create scalar data space for the attribute.
        call h5screate_simple_f(arank, adims, aspace_id, hdfier)
        if (hdfier .ne. 0) then
            write (6, '("sphdf5 : "10x" : error calling h5screate_simple_f from h5descr_cdset at $3:$4 ;")')
            call MPI_ABORT(MPI_COMM_SPEC, 1)
            stop "sphdf5 : error calling h5screate_simple_f from h5descr_cdset at $3:$4 ;"
        end if

        ! Create datatype for the attribute.
        call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, hdfier)
        if (hdfier .ne. 0) then
            write (6, '("sphdf5 : "10x" : error calling h5tcopy_f from h5descr_cdset at $3:$4 ;")')
            call MPI_ABORT(MPI_COMM_SPEC, 1)
            stop "sphdf5 : error calling h5tcopy_f from h5descr_cdset at $3:$4 ;"
        end if

        call h5tset_size_f(atype_id, attr_len, hdfier)
        if (hdfier .ne. 0) then
            write (6, '("sphdf5 : "10x" : error calling h5tset_size_f from h5descr_cdset at $3:$4 ;")')
            call MPI_ABORT(MPI_COMM_SPEC, 1)
            stop "sphdf5 : error calling h5tset_size_f from h5descr_cdset at $3:$4 ;"
        end if

        ! create descriptive attribute
        call h5acreate_f(dset_id, aname, atype_id, aspace_id, attr_id, hdfier)
        if (hdfier .ne. 0) then
            write (6, '("sphdf5 : "10x" : error calling h5acreate_f from h5descr_cdset at $3:$4 ;")')
            call MPI_ABORT(MPI_COMM_SPEC, 1)
            stop "sphdf5 : error calling h5acreate_f from h5descr_cdset at $3:$4 ;"
        end if

        ! Write the attribute data.
        call h5awrite_f(attr_id, atype_id, description, adims, hdfier)
        if (hdfier .ne. 0) then
            write (6, '("sphdf5 : "10x" : error calling h5awrite_f from h5descr_cdset at $3:$4 ;")')
            call MPI_ABORT(MPI_COMM_SPEC, 1)
            stop "sphdf5 : error calling h5awrite_f from h5descr_cdset at $3:$4 ;"
        end if

        ! Close the attribute.
        call h5aclose_f(attr_id, hdfier)
        if (hdfier .ne. 0) then
            write (6, '("sphdf5 : "10x" : error calling h5aclose_f from h5descr_cdset at $3:$4 ;")')
            call MPI_ABORT(MPI_COMM_SPEC, 1)
            stop "sphdf5 : error calling h5aclose_f from h5descr_cdset at $3:$4 ;"
        end if

        ! Close the attribute datatype.
        call h5tclose_f(atype_id, hdfier)
        if (hdfier .ne. 0) then
            write (6, '("sphdf5 : "10x" : error calling h5tclose_f from h5descr_cdset at $3:$4 ;")')
            call MPI_ABORT(MPI_COMM_SPEC, 1)
            stop "sphdf5 : error calling h5tclose_f from h5descr_cdset at $3:$4 ;"
        end if

        ! Terminate access to the data space.
        call h5sclose_f(aspace_id, hdfier)
        if (hdfier .ne. 0) then
            write (6, '("sphdf5 : "10x" : error calling h5sclose_f from h5descr_cdset at $3:$4 ;")')
            call MPI_ABORT(MPI_COMM_SPEC, 1)
            stop "sphdf5 : error calling h5sclose_f from h5descr_cdset at $3:$4 ;"
        end if

        call h5dclose_f(dset_id, hdfier)    ! terminate dataset;
        if (hdfier .ne. 0) then
            write (6, '("sphdf5 : "10x" : error calling h5dclose_f from h5descr_cdset at $3:$4 ;")')
            call MPI_ABORT(MPI_COMM_SPEC, 1)
            stop "sphdf5 : error calling h5dclose_f from h5descr_cdset at $3:$4 ;"
        end if
    end subroutine ! H5DESCR_CDSET

!> write logical variable _4 (scalar (_2=1) or rank-1 (_2=length)) into a dataset named _3 into group _1
!> example: hwritelv( grpInputGlobal, 1, LreadGF, (/ LreadGF /) ) ! scalar
!> example: hwritelv( grpInput,       5, success,  success(1:5) ) ! rank-1
    subroutine HWRITELV(group_id, name, num_data, data)
        integer(hid_t), intent(in) :: group_id
        character(len=*), intent(in) :: name
        integer, intent(in) :: num_data
        logical, dimension(:), intent(in) :: data

        integer :: hdfier     !< error flag for HDF5 library
        logical :: var_exists !< flags used to signal if a variable already exists

        integer(hsize_t) :: onedims(1) !< dimension specifier for one-dimensional data used in macros
        integer(hid_t) :: space_id   !< default dataspace ID used in macros
        integer(hid_t) :: dset_id    !< default dataset ID used in macros

        onedims(1) = num_data

        if (num_data .le. 0) then
            write (6, '("sphdf5 : "10x" : error calling hwriteiv ; ",a," : ",i1," .le. 0")') name, num_data
        else

            call h5screate_simple_f(1, onedims, space_id, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5screate_simple_f from hwritelv at $5:$6 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5screate_simple_f from hwritelv at $5:$6 ;"
            end if

            ! temporarily disable error printing to not confuse users
            call h5eset_auto_f(0, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5eset_auto_f from hwritelv at $5:$6 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5eset_auto_f from hwritelv at $5:$6 ;"
            end if

            ! check if dataset can be opened
            call h5dopen_f(group_id, name, dset_id, hdfier)
            if (hdfier .lt. 0) then
                var_exists = .false.
                if (hdfDebug .and. myid .eq. 0) then
                    write (*, *) "dataset", name, "does not exist yet, creating it"
                end if
            else
                var_exists = .true.
                if (hdfDebug .and. myid .eq. 0) then
                    write (*, *) "dataset", name, "exists already, opening it"
                end if
            end if

            ! re-establish previous state of error printing to be sensitive to "real" errors
            call h5eset_auto_f(internalHdf5Msg, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5eset_auto_f from hwritelv at $5:$6 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5eset_auto_f from hwritelv at $5:$6 ;"
            end if

            ! if the dataset does not exist already, create it. Otherwise, it should be open already
            if (.not. var_exists) then
                call h5dcreate_f(group_id, name, H5T_NATIVE_INTEGER, space_id, dset_id, hdfier)
                if (hdfier .ne. 0) then
                    write (6, '("sphdf5 : "10x" : error calling h5dcreate_f from hwritelv at $5:$6 ;")')
                    call MPI_ABORT(MPI_COMM_SPEC, 1)
                    stop "sphdf5 : error calling h5dcreate_f from hwritelv at $5:$6 ;"
                end if
            end if

            call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, merge(1, 0, data), onedims, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5dwrite_f from hwritelv at $5:$6 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5dwrite_f from hwritelv at $5:$6 ;"
            end if

            call h5dclose_f(dset_id, hdfier)    ! close dataset;
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5dclose_f from hwritelv at $5:$6 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5dclose_f from hwritelv at $5:$6 ;"
            end if

        end if
    end subroutine ! HWRITELV

! write logical variable _4 (scalar (_2=1) or rank-1 (_2=length)) into a dataset named _3 into group _1 and leave dataset open for e.g. adding an attribute; _5 and _6 should be __FILE__ and __LINE__
! example: hwritelv_lo( grpInputGlobal, 1, LreadGF, (/ LreadGF /) ) ! scalar
! example: hwritelv_lo( grpInput, 5, success, success(1:5) ) ! rank-1
! and close it using h5descr_cdset( /input/global/LreadGF, reading flag for GF )
    subroutine HWRITELV_LO(group_id, name, num_data, data, dset_id)
        integer(hid_t), intent(in) :: group_id
        character(len=*), intent(in) :: name
        integer, intent(in) :: num_data
        logical, dimension(:), intent(in) :: data
        integer(hid_t), intent(out) :: dset_id    !< default dataset ID used in macros

        integer :: hdfier     !< error flag for HDF5 library
        logical :: var_exists !< flags used to signal if a variable already exists

        integer(hsize_t) :: onedims(1) !< dimension specifier for one-dimensional data used in macros
        integer(hid_t) :: space_id   !< default dataspace ID used in macros

        onedims(1) = num_data

        if (num_data .le. 0) then
            write (6, '("sphdf5 : "10x" : error calling hwriteiv ; ",a," : ",i1," .le. 0")') name, num_data
        else

            call h5screate_simple_f(1, onedims, space_id, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5screate_simple_f from hwritelv at $5:$6 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5screate_simple_f from hwritelv at $5:$6 ;"
            end if

            ! temporarily disable error printing to not confuse users
            call h5eset_auto_f(0, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5eset_auto_f from hwritelv at $5:$6 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5eset_auto_f from hwritelv at $5:$6 ;"
            end if

            ! check if dataset can be opened
            call h5dopen_f(group_id, name, dset_id, hdfier)
            if (hdfier .lt. 0) then
                var_exists = .false.
                if (hdfDebug .and. myid .eq. 0) then
                    write (*, *) "dataset", name, "does not exist yet, creating it"
                end if
            else
                var_exists = .true.
                if (hdfDebug .and. myid .eq. 0) then
                    write (*, *) "dataset", name, "exists already, opening it"
                end if
            end if

            ! re-establish previous state of error printing to be sensitive to "real" errors
            call h5eset_auto_f(internalHdf5Msg, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5eset_auto_f from hwritelv at $5:$6 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5eset_auto_f from hwritelv at $5:$6 ;"
            end if

            ! if the dataset does not exist already, create it. Otherwise, it should be open already
            if (.not. var_exists) then
                call h5dcreate_f(group_id, name, H5T_NATIVE_INTEGER, space_id, dset_id, hdfier)
                if (hdfier .ne. 0) then
                    write (6, '("sphdf5 : "10x" : error calling h5dcreate_f from hwritelv at $5:$6 ;")')
                    call MPI_ABORT(MPI_COMM_SPEC, 1)
                    stop "sphdf5 : error calling h5dcreate_f from hwritelv at $5:$6 ;"
                end if
            end if

            call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, merge(1, 0, data), onedims, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5dwrite_f from hwritelv at $5:$6 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5dwrite_f from hwritelv at $5:$6 ;"
            end if

        end if
    end subroutine ! HWRITELV_LO

! write integer variable _4 (scalar (_2=1) or rank-1 (_2=length)) into a dataset named _3 into group _1; _5 and _6 should be __FILE__ and __LINE__
! example: hwriteiv( grpInputPhysics,    1, Igeometry, (/ Igeometry /) ) ! scalar
! example: hwriteiv( grpInputPhysics, Mvol,      Lrad,  Lrad(1:Mvol)   ) ! rank-1
    subroutine HWRITEIV(group_id, name, num_data, data)
        integer(hid_t), intent(in) :: group_id
        character(len=*), intent(in) :: name
        integer, intent(in) :: num_data
        integer, dimension(:), intent(in) :: data

        integer :: hdfier     !< error flag for HDF5 library
        logical :: var_exists !< flags used to signal if a variable already exists

        integer(hsize_t) :: onedims(1) !< dimension specifier for one-dimensional data used in macros
        integer(hid_t) :: space_id   !< default dataspace ID used in macros
        integer(hid_t) :: dset_id    !< default dataset ID used in macros

        onedims(1) = num_data

        if (num_data .le. 0) then
            write (6, '("sphdf5 : "10x" : error calling hwriteiv ; $3 : $2.le.0 at $5:$6 ;")')
        else

            call h5screate_simple_f(1, onedims, space_id, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5screate_simple_f from hwriteiv at $5:$6 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5screate_simple_f from hwriteiv at $5:$6 ;"
            end if

            ! temporarily disable error printing to not confuse users
            call h5eset_auto_f(0, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5eset_auto_f from hwriteiv at $5:$6 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5eset_auto_f from hwriteiv at $5:$6 ;"
            end if

            ! check if dataset can be opened
            call h5dopen_f(group_id, name, dset_id, hdfier)
            if (hdfier .lt. 0) then
                var_exists = .false.
                if (hdfDebug .and. myid .eq. 0) then; write (*, *) "dataset $3 does not exist yet, creating it"
                end if
            else
                var_exists = .true.
                if (hdfDebug .and. myid .eq. 0) then; write (*, *) "dataset $3 exists already, opening it"
                end if
            end if

            ! re-establish previous state of error printing to be sensitive to "real" errors
            call h5eset_auto_f(internalHdf5Msg, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5eset_auto_f from hwriteiv at $5:$6 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5eset_auto_f from hwriteiv at $5:$6 ;"
            end if

            ! if the dataset does not exist already, create it. Otherwise, it should be open already
            if (.not. var_exists) then
                call h5dcreate_f(group_id, name, H5T_NATIVE_INTEGER, space_id, dset_id, hdfier)
                if (hdfier .ne. 0) then
                    write (6, '("sphdf5 : "10x" : error calling h5dcreate_f from hwriteiv at $5:$6 ;")')
                    call MPI_ABORT(MPI_COMM_SPEC, 1)
                    stop "sphdf5 : error calling h5dcreate_f from hwriteiv at $5:$6 ;"
                end if
            end if

            call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data, onedims, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5dwrite_f from hwriteiv at $5:$6 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5dwrite_f from hwriteiv at $5:$6 ;"
            end if

            call h5dclose_f(dset_id, hdfier)    ! terminate dataset;
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5dclose_f from hwriteiv at $5:$6 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5dclose_f from hwriteiv at $5:$6 ;"
            end if

        end if
    end subroutine ! HWRITEIV

! write integer variable _4 (scalar (_2=1) or rank-1 (_2=length)) into a dataset named _3 into group _1 and leave the dataset open for e.g. adding and attribute
! example: hwriteiv( grpInputPhysics,    1, Igeometry, (/ Igeometry /) ) ! scalar
! example: hwriteiv( grpInputPhysics, Mvol,      Lrad,  Lrad(1:Mvol)   ) ! rank-1
! and close it using h5descr_cdset( /input/physics/Igeometry, geometry identifier )
    subroutine HWRITEIV_LO(group_id, name, num_data, data, dset_id)
        integer(hid_t), intent(in) :: group_id
        character(len=*), intent(in) :: name
        integer, intent(in) :: num_data
        integer, dimension(:), intent(in) :: data
        integer(hid_t), intent(out) :: dset_id    !< default dataset ID used in macros

        integer :: hdfier     !< error flag for HDF5 library
        logical :: var_exists !< flags used to signal if a variable already exists

        integer(hsize_t) :: onedims(1) !< dimension specifier for one-dimensional data used in macros
        integer(hid_t) :: space_id   !< default dataspace ID used in macros

        onedims(1) = num_data

        if (num_data .le. 0) then
            write (6, '("sphdf5 : "10x" : error calling hwriteiv ; $3 : $2.le.0 at $5:$6 ;")')
        else

            call h5screate_simple_f(1, onedims, space_id, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5screate_simple_f from hwriteiv at $5:$6 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5screate_simple_f from hwriteiv at $5:$6 ;"
            end if

            ! temporarily disable error printing to not confuse users
            call h5eset_auto_f(0, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5eset_auto_f from hwriteiv at $5:$6 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5eset_auto_f from hwriteiv at $5:$6 ;"
            end if

            ! check if dataset can be opened
            call h5dopen_f(group_id, name, dset_id, hdfier)
            if (hdfier .lt. 0) then
                var_exists = .false.
                if (hdfDebug .and. myid .eq. 0) then; write (*, *) "dataset $3 does not exist yet, creating it"
                end if
            else
                var_exists = .true.
                if (hdfDebug .and. myid .eq. 0) then; write (*, *) "dataset $3 exists already, opening it"
                end if
            end if

            ! re-establish previous state of error printing to be sensitive to "real" errors
            call h5eset_auto_f(internalHdf5Msg, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5eset_auto_f from hwriteiv at $5:$6 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5eset_auto_f from hwriteiv at $5:$6 ;"
            end if

            ! if the dataset does not exist already, create it. Otherwise, it should be open already
            if (.not. var_exists) then
                call h5dcreate_f(group_id, name, H5T_NATIVE_INTEGER, space_id, dset_id, hdfier)
                if (hdfier .ne. 0) then
                    write (6, '("sphdf5 : "10x" : error calling h5dcreate_f from hwriteiv at $5:$6 ;")')
                    call MPI_ABORT(MPI_COMM_SPEC, 1)
                    stop "sphdf5 : error calling h5dcreate_f from hwriteiv at $5:$6 ;"
                end if
            end if

            call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data, onedims, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5dwrite_f from hwriteiv at $5:$6 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5dwrite_f from hwriteiv at $5:$6 ;"
            end if

        end if
    end subroutine ! HWRITEIV_LO

! write real variable _4 (scalar (_2=1) or rank-1 (_2=length)) into a dataset named _3 into group _1; _5 and _6 should be __LINE__ and __FILE__
! example: hwriterv( grpInputPhysics,    1, phiedge, (/ phiedge /) ) ! scalar
! example: hwriterv( grpInputPhysics, Mvol,   tflux, tflux(1:Mvol) ) ! rank-1
    subroutine HWRITERV(group_id, name, num_data, data)
        integer(hid_t), intent(in) :: group_id
        character(len=*), intent(in) :: name
        integer, intent(in) :: num_data
        real(wp), dimension(:), intent(in) :: data

        integer :: hdfier     !< error flag for HDF5 library
        logical :: var_exists !< flags used to signal if a variable already exists

        integer(hsize_t) :: onedims(1) !< dimension specifier for one-dimensional data used in macros
        integer(hid_t) :: space_id   !< default dataspace ID used in macros
        integer(hid_t) :: dset_id    !< default dataset ID used in macros

        onedims(1) = num_data

        if (num_data .le. 0) then
            write (6, '("sphdf5 : "10x" : error calling hwriterv ; $3 : $2.le.0 at $5:$6 ;")')
        else

            call h5screate_simple_f(1, onedims, space_id, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5screate_simple_f from hwriterv at $5:$6 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5screate_simple_f from hwriterv at $5:$6 ;"
            end if

            ! temporarily disable error printing to not confuse users
            call h5eset_auto_f(0, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5eset_auto_f from hwriterv at $5:$6 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5eset_auto_f from hwriterv at $5:$6 ;"
            end if

            ! check if dataset can be opened
            call h5dopen_f(group_id, name, dset_id, hdfier)
            if (hdfier .lt. 0) then
                var_exists = .false.
                if (hdfDebug .and. myid .eq. 0) then; write (*, *) "dataset $3 does not exist yet, creating it"
                end if
            else
                var_exists = .true.
                if (hdfDebug .and. myid .eq. 0) then; write (*, *) "dataset $3 exists already, opening it"
                end if
            end if

            ! re-establish previous state of error printing to be sensitive to "real" errors
            call h5eset_auto_f(internalHdf5Msg, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5eset_auto_f from hwriterv at $5:$6 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5eset_auto_f from hwriterv at $5:$6 ;"
            end if

            ! if the dataset does not exist already, create it. Otherwise, it should be open already
            if (.not. var_exists) then
                call h5dcreate_f(group_id, name, H5T_NATIVE_DOUBLE, space_id, dset_id, hdfier)
                if (hdfier .ne. 0) then
                    write (6, '("sphdf5 : "10x" : error calling h5dcreate_f from hwriterv at $5:$6 ;")')
                    call MPI_ABORT(MPI_COMM_SPEC, 1)
                    stop "sphdf5 : error calling h5dcreate_f from hwriterv at $5:$6 ;"
                end if
            end if

            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, onedims, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5dwrite_f from hwriterv at $5:$6 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5dwrite_f from hwriterv at $5:$6 ;"
            end if

            call h5dclose_f(dset_id, hdfier)    ! terminate dataset;
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5dclose_f from hwriterv at $5:$6 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5dclose_f from hwriterv at $5:$6 ;"
            end if

        end if
    end subroutine ! HWRITERV

! write real variable _4 (scalar (_2=1) or rank-1 (_2=length)) into a dataset named _3 into group _1 and leave it open, e.g. for adding an attribute; _5 and _6 should be __FILE__ and __LINE__
! example: hwriterv( grpInputPhysics,    1, phiedge, (/ phiedge /) ) ! scalar
! example: hwriterv( grpInputPhysics, Mvol,   tflux, tflux(1:Mvol) ) ! rank-1
! and close it with h5descr_cdset( /input/physics/phiedge, total enclosed toroidal flux )
    subroutine HWRITERV_LO(group_id, name, num_data, data, dset_id)
        integer(hid_t), intent(in) :: group_id
        character(len=*), intent(in) :: name
        integer, intent(in) :: num_data
        real(wp), dimension(:), intent(in) :: data
        integer(hid_t), intent(out) :: dset_id    !< default dataset ID used in macros

        integer :: hdfier     !< error flag for HDF5 library
        logical :: var_exists !< flags used to signal if a variable already exists

        integer(hsize_t) :: onedims(1) !< dimension specifier for one-dimensional data used in macros
        integer(hid_t) :: space_id   !< default dataspace ID used in macros

        onedims(1) = num_data

        if (num_data .le. 0) then
            write (6, '("sphdf5 : "10x" : error calling hwriterv ; $3 : $2.le.0 at $5:$6 ;")')
        else

            call h5screate_simple_f(1, onedims, space_id, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5screate_simple_f from hwriterv at $5:$6 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5screate_simple_f from hwriterv at $5:$6 ;"
            end if

            ! temporarily disable error printing to not confuse users
            call h5eset_auto_f(0, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5eset_auto_f from hwriterv at $5:$6 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5eset_auto_f from hwriterv at $5:$6 ;"
            end if

            ! check if dataset can be opened
            call h5dopen_f(group_id, name, dset_id, hdfier)
            if (hdfier .lt. 0) then
                var_exists = .false.
                if (hdfDebug .and. myid .eq. 0) then; write (*, *) "dataset $3 does not exist yet, creating it"
                end if
            else
                var_exists = .true.
                if (hdfDebug .and. myid .eq. 0) then; write (*, *) "dataset $3 exists already, opening it"
                end if
            end if

            ! re-establish previous state of error printing to be sensitive to "real" errors
            call h5eset_auto_f(internalHdf5Msg, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5eset_auto_f from hwriterv at $5:$6 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5eset_auto_f from hwriterv at $5:$6 ;"
            end if

            ! if the dataset does not exist already, create it. Otherwise, it should be open already
            if (.not. var_exists) then
                call h5dcreate_f(group_id, name, H5T_NATIVE_DOUBLE, space_id, dset_id, hdfier)
                if (hdfier .ne. 0) then
                    write (6, '("sphdf5 : "10x" : error calling h5dcreate_f from hwriterv at $5:$6 ;")')
                    call MPI_ABORT(MPI_COMM_SPEC, 1)
                    stop "sphdf5 : error calling h5dcreate_f from hwriterv at $5:$6 ;"
                end if
            end if

            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, onedims, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5dwrite_f from hwriterv at $5:$6 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5dwrite_f from hwriterv at $5:$6 ;"
            end if

        end if
    end subroutine ! HWRITERV_LO

! write real array _5 (_2 rows, _3 columns) into a dataset named _4 into group _1; _6 and _7 should be __FILE__ and __LINE__
! example: hwritera( grpInputPhysics, (2*Ntor+1), (2*Mpol+1), Rbc, Rbc(-Ntor:Ntor,-Mpol:Mpol) )
    subroutine HWRITERA(group_id, name, num_rows, num_cols, data)
        integer(hid_t), intent(in) :: group_id
        character(len=*), intent(in) :: name
        integer, intent(in) :: num_rows, num_cols
        real(wp), dimension(:, :), intent(in) :: data

        integer :: hdfier     !< error flag for HDF5 library
        logical :: var_exists !< flags used to signal if a variable already exists

        integer(hsize_t) :: twodims(2) !< dimension specifier for two-dimensional data used in macros
        integer(hid_t) :: space_id   !< default dataspace ID used in macros
        integer(hid_t) :: dset_id    !< default dataset ID used in macros

        twodims(1:2) = (/num_rows, num_cols/)

        if (num_rows .le. 0 .or. num_cols .le. 0) then
            write (6, '("sphdf5 : "10x" : error calling hwritera ; $4 : $2.le.0 .or. $3.le.0 at $6:$7 ;")')
        else

            call h5screate_simple_f(2, twodims, space_id, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5screate_simple_f from hwritera at $6:$7 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5screate_simple_f from hwritera at $6:$7 ;"
            end if

            ! temporarily disable error printing to not confuse users
            call h5eset_auto_f(0, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5eset_auto_f from hwritera at $6:$7 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5eset_auto_f from hwritera at $6:$7 ;"
            end if

            ! check if dataset can be opened
            call h5dopen_f(group_id, name, dset_id, hdfier)
            if (hdfier .lt. 0) then
                var_exists = .false.
                if (hdfDebug .and. myid .eq. 0) then; write (*, *) "dataset $4 does not exist yet, creating it"
                end if
            else
                var_exists = .true.
                if (hdfDebug .and. myid .eq. 0) then; write (*, *) "dataset $4 exists already, opening it"
                end if
            end if

            ! re-establish previous state of error printing to be sensitive to "real" errors
            call h5eset_auto_f(internalHdf5Msg, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5eset_auto_f from hwritera at $6:$7 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5eset_auto_f from hwritera at $6:$7 ;"
            end if

            ! if the dataset does not exist already, create it. Otherwise, it should be open already
            if (.not. var_exists) then
                call h5dcreate_f(group_id, name, H5T_NATIVE_DOUBLE, space_id, dset_id, hdfier)
                if (hdfier .ne. 0) then
                    write (6, '("sphdf5 : "10x" : error calling h5dcreate_f from hwritera at $6:$7 ;")')
                    call MPI_ABORT(MPI_COMM_SPEC, 1)
                    stop "sphdf5 : error calling h5dcreate_f from hwritera at $6:$7 ;"
                end if
            end if

            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, twodims, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5dwrite_f from hwritera at $6:$7 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5dwrite_f from hwritera at $6:$7 ;"
            end if

            call h5dclose_f(dset_id, hdfier)    ! terminate dataset;
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5dclose_f from hwritera at $6:$7 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5dclose_f from hwritera at $6:$7 ;"
            end if

        end if
    end subroutine ! HWRITERA

! write real array _5 (_2 rows, _3 columns) into a dataset named _4 into group _1 and leave it open, e.g. for adding an attribute; _6 and _7 should be __FILE__ and __LINE__
! example: hwritera( grpInputPhysics, (2*Ntor+1), (2*Mpol+1), Rbc, Rbc(-Ntor:Ntor,-Mpol:Mpol) )
! and close it then via h5descr_cdset( /input/physics/Rbc, boundary R cosine Fourier coefficients )
    subroutine HWRITERA_LO(group_id, name, num_rows, num_cols, data, dset_id)
        integer(hid_t), intent(in) :: group_id
        character(len=*), intent(in) :: name
        integer, intent(in) :: num_rows, num_cols
        real(wp), dimension(:, :), intent(in) :: data
        integer(hid_t), intent(out) :: dset_id    !< default dataset ID used in macros

        integer :: hdfier     !< error flag for HDF5 library
        logical :: var_exists !< flags used to signal if a variable already exists

        integer(hsize_t) :: twodims(2) !< dimension specifier for two-dimensional data used in macros
        integer(hid_t) :: space_id   !< default dataspace ID used in macros

        twodims(1:2) = (/num_rows, num_cols/)

        if (num_rows .le. 0 .or. num_cols .le. 0) then
            write (6, '("sphdf5 : "10x" : error calling hwritera ; $4 : $2.le.0 .or. $3.le.0 at $6:$7 ;")')
        else

            call h5screate_simple_f(2, twodims, space_id, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5screate_simple_f from hwritera at $6:$7 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5screate_simple_f from hwritera at $6:$7 ;"
            end if

            ! temporarily disable error printing to not confuse users
            call h5eset_auto_f(0, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5eset_auto_f from hwritera at $6:$7 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5eset_auto_f from hwritera at $6:$7 ;"
            end if

            ! check if dataset can be opened
            call h5dopen_f(group_id, name, dset_id, hdfier)
            if (hdfier .lt. 0) then
                var_exists = .false.
                if (hdfDebug .and. myid .eq. 0) then; write (*, *) "dataset $4 does not exist yet, creating it"
                end if
            else
                var_exists = .true.
                if (hdfDebug .and. myid .eq. 0) then; write (*, *) "dataset $4 exists already, opening it"
                end if
            end if

            ! re-establish previous state of error printing to be sensitive to "real" errors
            call h5eset_auto_f(internalHdf5Msg, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5eset_auto_f from hwritera at $6:$7 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5eset_auto_f from hwritera at $6:$7 ;"
            end if

            ! if the dataset does not exist already, create it. Otherwise, it should be open already
            if (.not. var_exists) then
                call h5dcreate_f(group_id, name, H5T_NATIVE_DOUBLE, space_id, dset_id, hdfier)
                if (hdfier .ne. 0) then
                    write (6, '("sphdf5 : "10x" : error calling h5dcreate_f from hwritera at $6:$7 ;")')
                    call MPI_ABORT(MPI_COMM_SPEC, 1)
                    stop "sphdf5 : error calling h5dcreate_f from hwritera at $6:$7 ;"
                end if
            end if

            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, twodims, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5dwrite_f from hwritera at $6:$7 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5dwrite_f from hwritera at $6:$7 ;"
            end if

        end if
    end subroutine ! HWRITERA_LO

! write real cube _6 (_2 rows, _3 columns, _4 pages) into a dataset named _5 into group _1; _7 and _8 should containt __FILE__ and __LINE__
! example: hwriterc( grpOutput, (Mrad+1), 2, 2, TT, TT(0:Mrad,0:1,0:1) )
    subroutine HWRITERC(group_id, name, num_rows, num_cols, num_pages, data)
        integer(hid_t), intent(in) :: group_id
        character(len=*), intent(in) :: name
        integer, intent(in) :: num_rows, num_cols, num_pages
        real(wp), dimension(:, :, :), intent(in) :: data

        integer :: hdfier     !< error flag for HDF5 library
        logical :: var_exists !< flags used to signal if a variable already exists

        integer(hsize_t) :: threedims(3) !< dimension specifier for three-dimensional data used in macros
        integer(hid_t) :: space_id   !< default dataspace ID used in macros
        integer(hid_t) :: dset_id    !< default dataset ID used in macros

        threedims(1:3) = (/num_rows, num_cols, num_pages/)

        if (num_rows .le. 0 .or. num_cols .le. 0 .or. num_pages .le. 0) then

            write (6, '("sphdf5 : "10x" : error calling hwriterc ; $5 : $2.le.0 .or. $3.le.0 .or. $4.le.0 at $7:$8 ;")')

        else

            call h5screate_simple_f(3, threedims, space_id, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5screate_simple_f from hwriterc at $7:$8 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5screate_simple_f from hwriterc at $7:$8 ;"
            end if

            ! temporarily disable error printing to not confuse users
            call h5eset_auto_f(0, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5eset_auto_f from hwriterc at $7:$8 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5eset_auto_f from hwriterc at $7:$8 ;"
            end if

            ! check if dataset can be opened
            call h5dopen_f(group_id, name, dset_id, hdfier)
            if (hdfier .lt. 0) then
                var_exists = .false.
                if (hdfDebug .and. myid .eq. 0) then; write (*, *) "dataset $5 does not exist yet, creating it"
                end if
            else
                var_exists = .true.
                if (hdfDebug .and. myid .eq. 0) then; write (*, *) "dataset $5 exists already, opening it"
                end if
            end if

            ! re-establish previous state of error printing to be sensitive to "real" errors
            call h5eset_auto_f(internalHdf5Msg, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5eset_auto_f from hwriterc at $7:$8 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5eset_auto_f from hwriterc at $7:$8 ;"
            end if

            ! if the dataset does not exist already, create it. Otherwise, it should be open already
            if (.not. var_exists) then
                call h5dcreate_f(group_id, name, H5T_NATIVE_DOUBLE, space_id, dset_id, hdfier)
                if (hdfier .ne. 0) then
                    write (6, '("sphdf5 : "10x" : error calling h5dcreate_f from hwriterc at $7:$8 ;")')
                    call MPI_ABORT(MPI_COMM_SPEC, 1)
                    stop "sphdf5 : error calling h5dcreate_f from hwriterc at $7:$8 ;"
                end if
            end if

            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, threedims, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5dwrite_f from hwriterc at $7:$8 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5dwrite_f from hwriterc at $7:$8 ;"
            end if

            call h5dclose_f(dset_id, hdfier)    ! terminate dataset;
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5dclose_f from hwriterc at $7:$8 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5dclose_f from hwriterc at $7:$8 ;"
            end if

        end if
    end subroutine ! HWRITERC

! write real cube _6 (_2 rows, _3 columns, _4 pages) into a dataset named _5 into group _1 and leave open for e.g. adding an attribute; _7 and _8 should be __FILE__ and __LINE__
! example: hwriterc( grpOutput, (Mrad+1), 2, 2, TT, TT(0:Mrad,0:1,0:1) )
! and close it with h5descr_cdset( /output/TT, something abbreviated by TT )
    subroutine HWRITERC_LO(group_id, name, num_rows, num_cols, num_pages, data, dset_id)
        integer(hid_t), intent(in) :: group_id
        character(len=*), intent(in) :: name
        integer, intent(in) :: num_rows, num_cols, num_pages
        real(wp), dimension(:, :, :), intent(in) :: data
        integer(hid_t), intent(out) :: dset_id    !< default dataset ID used in macros

        integer :: hdfier     !< error flag for HDF5 library
        logical :: var_exists !< flags used to signal if a variable already exists

        integer(hsize_t) :: threedims(3) !< dimension specifier for three-dimensional data used in macros
        integer(hid_t) :: space_id   !< default dataspace ID used in macros

        threedims(1:3) = (/num_rows, num_cols, num_pages/)

        if (num_rows .le. 0 .or. num_cols .le. 0 .or. num_pages .le. 0) then

            write (6, '("sphdf5 : "10x" : error calling hwriterc ; $5 : $2.le.0 .or. $3.le.0 .or. $4.le.0 at $7:$8 ;")')

        else

            call h5screate_simple_f(3, threedims, space_id, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5screate_simple_f from hwriterc at $7:$8 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5screate_simple_f from hwriterc at $7:$8 ;"
            end if

            ! temporarily disable error printing to not confuse users
            call h5eset_auto_f(0, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5eset_auto_f from hwriterc at $7:$8 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5eset_auto_f from hwriterc at $7:$8 ;"
            end if

            ! check if dataset can be opened
            call h5dopen_f(group_id, name, dset_id, hdfier)
            if (hdfier .lt. 0) then
                var_exists = .false.
                if (hdfDebug .and. myid .eq. 0) then; write (*, *) "dataset $5 does not exist yet, creating it"
                end if
            else
                var_exists = .true.
                if (hdfDebug .and. myid .eq. 0) then; write (*, *) "dataset $5 exists already, opening it"
                end if
            end if

            ! re-establish previous state of error printing to be sensitive to "real" errors
            call h5eset_auto_f(internalHdf5Msg, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5eset_auto_f from hwriterc at $7:$8 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5eset_auto_f from hwriterc at $7:$8 ;"
            end if

            ! if the dataset does not exist already, create it. Otherwise, it should be open already
            if (.not. var_exists) then
                call h5dcreate_f(group_id, name, H5T_NATIVE_DOUBLE, space_id, dset_id, hdfier)
                if (hdfier .ne. 0) then
                    write (6, '("sphdf5 : "10x" : error calling h5dcreate_f from hwriterc at $7:$8 ;")')
                    call MPI_ABORT(MPI_COMM_SPEC, 1)
                    stop "sphdf5 : error calling h5dcreate_f from hwriterc at $7:$8 ;"
                end if
            end if

            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, threedims, hdfier)
            if (hdfier .ne. 0) then
                write (6, '("sphdf5 : "10x" : error calling h5dwrite_f from hwriterc at $7:$8 ;")')
                call MPI_ABORT(MPI_COMM_SPEC, 1)
                stop "sphdf5 : error calling h5dwrite_f from hwriterc at $7:$8 ;"
            end if

        end if
    end subroutine ! HWRITERC_LO

end module ! h5utils
