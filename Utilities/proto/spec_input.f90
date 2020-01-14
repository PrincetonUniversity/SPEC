!> @file spec_input.f90

module fileunits

  implicit none

  integer :: iunit = 10 !< input namelists text file
  integer :: ounit =  6 !< screen output

end module fileunits

module spec_version
  implicit none
  integer, parameter :: version(3) = (/ 2, 1, 0 /) !< version of SPEC
end module spec_version

! compare two version numbers to be less or equal, i.e. compatible
pure function version_le(query_version, ref_version)
  integer, intent(in)  :: query_version(3)
  integer, intent(in)  :: ref_version(3)
  logical              :: version_le
  version_le = (query_version(1).le.ref_version(1) .and. &
             &  query_version(2).le.ref_version(2) .and. &
             &  query_version(3).le.ref_version(3))
end function version_le

module inputlist
  
  implicit none
  
  logical :: verbose = .TRUE. ! set to true to enable verbal diarreha during reading of input file
  
  ! maximum resolution, i.e. array dimensions
  integer, parameter :: MNvol     = 256 !< maximum value of \c Nvol
  integer, parameter :: MMpol     =  32 !< maximum value of \c Mpol
  integer, parameter :: MNtor     =  16 !< maximum value of \c Ntor
  
  !> maximum length of filenames on common operating systems in the year 2020
  !> \see https://en.wikipedia.org/wiki/Comparison_of_file_systems#Limits
  character :: ext*255
  
  ! input namelist 'physicslist'
  integer :: Igeometry        = 3
  integer :: Nvol             = 1
  integer :: Lrad(1:MNvol+1)  = 4
  double precision :: phiedge = 1.0
  
  namelist /physicslist/ &
  Igeometry,&
  Nvol,&
  Lrad,&
  phiedge
  
  contains
  
  !> @brief read the input data for SPEC from a textual namelist in the file specified by filename
  subroutine readin(filename)
    use fileunits
    
    implicit none
    
    character(len=*),intent(in) :: filename ! name of input file
    
    logical :: Lspexist !< flag to indicate that the file \c ext.sp exists
    
    inquire( file=filename, exist=Lspexist ) ! check if file exists
    if (Lspexist) then
      write(ounit,*) "reading input from text file '",filename,"'"
      open( iunit, file=filename, status="old")
      read( iunit, physicslist)
      close(iunit)
    else
      write(ounit,*) "the input file '",filename,"' does not exist"
      call exit(1) ! "generic fail" return code
    endif
    
  end subroutine readin
  
  !> @brief read the input data for SPEC from a HDF5 file specified by filename_h5
  subroutine readin_h5(filename_h5)
    use hdf5
    use fileunits
    use spec_version
    
    implicit none
    
    logical :: version_le ! function to compare version numbers
    
    character(len=*),intent(in) :: filename_h5 ! name of input file
    
    integer                     :: hdfier            ! error flag for HDF5 library calls
    logical                     :: Lsph5exist        !< flag to indicate that the file \c ext.sp.h5 exists
    integer(hid_t)              :: dtype_id, dtype_id_native
    integer(hid_t)              :: dataspace         ! dataspace used to query Dataset size
    integer                     :: rank              ! rank of a dataspace, i.e. number of dimensions
    integer(hsize_t)            :: dims_1(1)         ! current dimensions of rank-1 Dataset
    integer(hsize_t)            :: dims_2(2)         ! current dimensions of rank-2 Dataset
    integer(hsize_t)            :: dims_3(3)         ! current dimensions of rank-3 Dataset
    integer(hsize_t)            :: max_dims_1(1)     ! maximum dimensions of rank-1 Dataset
    integer(hsize_t)            :: max_dims_2(2)     ! maximum dimensions of rank-2 Dataset
    integer(hsize_t)            :: max_dims_3(3)     ! maximum dimensions of rank-3 Dataset
    logical                     :: item_exists, datatypes_equal
    integer                     :: item_type, dspace_type
    
    integer(hid_t)              :: file_id
    integer(hid_t)              :: dset_version
    integer(hid_t)              :: grp_input
    integer(hid_t)              :: dset_input_Igeometry
    integer(hid_t)              :: dset_input_Nvol
    integer(hid_t)              :: dset_input_Lrad
    integer(hid_t)              :: dset_input_phiedge
    
    integer                     :: input_version(3)  ! version of input data
    
    ! check if file exists
    inquire(file=filename_h5, exist=Lsph5exist)
    if (Lsph5exist) then
      write(ounit,*) "reading input from HDF5 file '",filename_h5,"'"
      
      ! open input file
      call h5fopen_f(filename_h5, H5F_ACC_RDONLY_F, file_id, hdfier)
      if (hdfier.ne.0) then
        write(*,*) "error opening input HDF5 file"
      else
        ! successfully opened input file
        
        ! The following scheme to verify the existence, correct type and reasonable value of the Dataset
        ! is outlined in the documentation of the e.g. h5oexists_by_name_f library routine:
        ! https://portal.hdfgroup.org/display/HDF5/H5O_EXISTS_BY_NAME
        ! The procedure consists essentially of the following steps logically nested within each other:
        !  1. check if the link with the given Dataset or Group name exists
        !  2. check that the link at the given name resolves to an object
        !  3. open the object at the given name
        !  4. query the object type of the just-opened object
        !  5. get the datatype of the object
        !  6. convert the datatype into the native datatype for endianess-agnostic comparison
        !  7. compare the native datatype with the expected native datatype
        !  8. read the Dataset using the verified native datatype
        !  9. close the native datatype
        ! 10. close the Dataset object
        
        write(ounit,*) " "
        
        ! query existence of /version link
        call h5lexists_f(file_id, "version", item_exists, hdfier)
        if (hdfier.ne.0) then
          write(*,*) "error checking if /version link exists"
        else
          if (verbose) then ; write(ounit,*) "successfully checked if the /version link exists" ; endif
        
          ! check that /version link exists
          if (item_exists) then
            if (verbose) then; write(ounit,*) "successfully checked that the link /version exists" ; endif
           
            ! query existence of object at /version link
            call h5oexists_by_name_f(file_id, "version", item_exists, hdfier)
            if (hdfier.ne.0) then
              write(*,*) "error checking for presence of object at /version"
            else
              if (verbose) then ; write(ounit,*) "successfully checked for presence of an object at /version link" ; endif
            
              ! check that there exists an item at the /version link
              if (item_exists) then
                if (verbose) then ; write(ounit,*) "successfully checked that there exists an item at /version" ; endif
              
                ! try to open /version object
                call h5oopen_f(file_id, "version", dset_version, hdfier)
                if (hdfier.ne.0) then
                  write(*,*) "error opening object /version"
                else
                  if (verbose) then ; write(ounit,*) "successfully opened object at /version" ; endif
                
                  ! query type of /version object
                  call h5iget_type_f(dset_version, item_type, hdfier)
                  if (hdfier.ne.0) then
                    write(*,*) "error querying type of object /version"
                  else
                    if (verbose) then ; write(ounit,*) "successfully queried item type of /version" ; endif
                  
                    ! check if /version is a Dataset
                    if (item_type.ne.H5I_DATASET_F) then
                      write(ounit,*) "error verifying that object /version is a Dataset(",H5I_DATASET_F,"); rather it is ",item_type
                    else
                      if (verbose) then ; write(ounit,*) "successfully verified that /version is a Dataset" ; endif
                      
                      ! open dataspace of /version
                      call h5dget_space_f(dset_version, dataspace, hdfier)
                      if (hdfier.ne.0) then
                        write(ounit,*) "error getting dataspace of /version"
                      else
                      
                        ! check that the dataspace of /version is simple
                        call h5sget_simple_extent_type_f(dataspace, dspace_type, hdfier)
                        if (hdfier.ne.0) then
                          write(ounit,*) "error getting type of /version dataspace"
                        else
                        
                          ! check that the type of the /version dataspace is H5S_SIMPLE_F
                          if (dspace_type.eq.H5S_SIMPLE_F) then
                            write(ounit,*) "successfully verified that the type of the /version dataspace is H5S_SIMPLE_F"
                      
                            ! query rank of /version
                            call h5sget_simple_extent_ndims_f(dataspace, rank, hdfier)
                            if (hdfier.ne.0) then
                              write(ounit,*) "error getting rank of /version; expected ",rank,", got ",hdfier
                            else
                              if (verbose) then ; write(ounit,*) "successfully queried rank of /version" ; endif
                              
                              ! get the current and maximum dimensions of /version
                              call h5sget_simple_extent_dims_f(dataspace, dims_1, max_dims_1, hdfier)
                              if (hdfier.ne.rank) then
                                write(ounit,*) "ERROR: rank mismatch of /version; expected ",rank,", but got ",hdfier
                              else
                                write(ounit,*) "successfully queried length of /version"
                                
                                ! check that /version has the correct length
                                if (dims_1(1).ne.3) then
                                  write(ounit,*) "ERROR: expected length of /version is 3, but got ", dims_1
                                else
                                  if (verbose) then ; write(ounit,*) "successfully verified the correct length of /version" ; endif
                                  
                                  ! query datatype of /version
                                  call h5dget_type_f(dset_version, dtype_id, hdfier)
                                  if (hdfier.ne.0) then
                                    write(ounit,*) "error querying datatype of /version"
                                  else
                                    if (verbose) then ; write(ounit,*) "successfully queried datatype of /version" ; endif
                                      
                                    ! convert datatype to native array for comparison
                                    call h5tget_native_type_f(dtype_id, 1, dtype_id_native, hdfier)
                                    if (hdfier.ne.0) then
                                      write(ounit,*) "error converting datatype of /version to native datatype"
                                    else
                                      if (verbose) then ; write(ounit,*) "successfully converted datatype of /version to native datatype" ; endif
                                      
                                      ! call comparison routine for /version datatype
                                      call h5tequal_f(dtype_id_native, H5T_NATIVE_INTEGER, datatypes_equal, hdfier)
                                      if (hdfier.ne.0) then
                                        write(ounit,*) "error comparing datatype of /version to H5T_NATIVE_INTEGER"
                                      else
                                        if (verbose) then ; write(ounit,*) "successfully executed comparison of datatype of /version with H5T_NATIVE_INTEGER" ; endif
                                        
                                        ! verify correct datatype of /version
                                        if (datatypes_equal) then
                                          if (verbose) then ; write(ounit,*) "successfully checked that datatype of /version is H5T_NATIVE_INTEGER :-)" ; endif
                                          
                                          ! read /version Dataset
                                          call h5dread_f(dset_version, H5T_NATIVE_INTEGER, input_version(1:3), int((/3/), HSIZE_T), hdfier)
                                          if (hdfier.ne.0) then
                                            write(*,*) "error reading Dataset /version"
                                          else
                                            if (verbose) then ; write(ounit,'(" successfully read /version from input data: ",2(i1,"."),i1)') input_version(1:3) ; endif
                                            
                                            ! verify that version number in input file is less than or equal to the current version of SPEC
                                            if (version_le(input_version, version)) then
                                              if (verbose) then
                                                write(ounit,'(" INFO: The version of the input data (",2(i1,"."),i1,") is compatible with SPEC ",2(i1,"."),i1)') &
                                                & input_version, version
                                              endif
                                            else
                                              write(ounit,'(" WARNING: the version of the input data (",2(i1,"."),i1,") may be incompatible with SPEC ",2(i1,"."),i1)') &
                                              & input_version, version
                                            endif ! verify that version number in input file is less than or equal to the current version of SPEC
                                          endif ! read /version Dataset
                                        else
                                          write(ounit,*) "ERROR: native datatype of /version should be H5T_NATIVE_INTEGER but is ",dtype_id_native
                                        endif !verify correct datatype of /version
                                      
                                      endif ! call comparison routine for /version datatype
                                    endif ! convert datatype to native array for comparison
                                    
                                    ! close datatype of /version
                                    call h5tclose_f(dtype_id, hdfier)
                                    if (hdfier.ne.0) then
                                      write(ounit,*) "error closing datatype of /version"
                                    elseif (verbose) then
                                      write(ounit,*) "successfully closed datatype of /version"
                                    endif ! close datatype of /version
                                  endif ! query datatype of /version
                                endif ! check that /version has the correct length
                              endif ! get the current and maximum dimensions of /version
                            endif ! query rank of /version
                          else
                            write(ounit,*) "ERROR: type of dataspace /version is not H5S_SIMPLE_F but ",dspace_type
                          endif ! check that the type of the /version dataspace is H5S_SIMPLE_F
                        endif ! check that the dataspace of /version is simple
                        
                        ! close dataspace of /version
                        call h5sclose_f(dataspace, hdfier)
                        if (hdfier.ne.0) then ; write(ounit,*) "error closing dataspace of /version" ; endif
                      endif ! open dataspace of /version
                    endif ! query type of /version object
                  endif ! check if /version is a Dataset
                          
                  call h5oclose_f(dset_version, hdfier)
                  if (hdfier.ne.0) then
                    write(*,*) "error closing object /version"
                  elseif (verbose) then
                    write(ounit,*) "successfully closed object /version"
                  endif
                endif ! try to open /version object
              else
                ! /version link present but does not resolve to any object
                write(ounit,*) "/version link present but does not resolve to any object"
              endif ! check that there exists an item at the /version link
            endif ! query existence of object at /version link
          else
            ! /version link not present in input file
            write(ounit,*) " WARNING: /version not found. Using default value "
          endif ! check if /version link exists
        endif ! query existence of /version link
        
        write(ounit,*) " "
        
        ! query existence of /input link
        call h5lexists_f(file_id, "input", item_exists, hdfier)
        if (hdfier.ne.0) then
          write(*,*) "error checking if /input link exists"
        else
          if (verbose) then ; write(ounit,*) "successfully checked if the /input link exists" ; endif
        
          ! check that /input link exists
          if (item_exists) then
            if (verbose) then; write(ounit,*) "successfully checked that the link /input exists" ; endif
           
            ! query existence of object at /input link
            call h5oexists_by_name_f(file_id, "input", item_exists, hdfier)
            if (hdfier.ne.0) then
              write(*,*) "error checking for presence of object at /input"
            else
              if (verbose) then ; write(ounit,*) "successfully checked for presence of an object at /input link" ; endif
            
              ! check that there exists an item at the /input link
              if (item_exists) then
                if (verbose) then ; write(ounit,*) "successfully checked that there exists an item at /input" ; endif
              
                ! try to open /input object
                call h5oopen_f(file_id, "input", grp_input, hdfier)
                if (hdfier.ne.0) then
                  write(*,*) "error opening object /input"
                else
                  if (verbose) then ; write(ounit,*) "successfully opened object at /input" ; endif
                
                  ! query type of /input object
                  call h5iget_type_f(grp_input, item_type, hdfier)
                  if (hdfier.ne.0) then
                    write(*,*) "error querying type of object /input"
                  else
                    if (verbose) then ; write(ounit,*) "successfully queried item type of /input" ; endif
                  
                    ! check if /input is a Group
                    if (item_type.ne.H5I_GROUP_F) then
                      write(ounit,*) "error verifying that object /input is a Group(",H5I_GROUP_F,"); rather it is ",item_type
                    else
                      if (verbose) then ; write(ounit,*) "successfully verified that /input is a Group" ; endif
                      
                      ! check for contents of /input Group
                      
                      write(ounit,*) " "
                      
                      ! query existence of /input/Igeometry link
                      call h5lexists_f(grp_input, "Igeometry", item_exists, hdfier)
                      if (hdfier.ne.0) then
                        write(*,*) "error checking if /input/Igeometry link exists"
                      else
                        if (verbose) then ; write(ounit,*) "successfully checked if the /input/Igeometry link exists" ; endif
                      
                        ! check that /input/Igeometry link exists
                        if (item_exists) then
                          if (verbose) then; write(ounit,*) "successfully checked that the link /input/Igeometry exists" ; endif
                         
                          ! query existence of object at /input/Igeometry link
                          call h5oexists_by_name_f(grp_input, "Igeometry", item_exists, hdfier)
                          if (hdfier.ne.0) then
                            write(*,*) "error checking for presence of object at /input/Igeometry"
                          else
                            if (verbose) then ; write(ounit,*) "successfully checked for presence of an object at /input/Igeometry link" ; endif
                          
                            ! check that there exists an item at the /input/Igeometry link
                            if (item_exists) then
                              if (verbose) then ; write(ounit,*) "successfully checked that there exists an item at /input/Igeometry" ; endif
                            
                              ! try to open /input/Igeometry object
                              call h5oopen_f(grp_input, "Igeometry", dset_input_Igeometry, hdfier)
                              if (hdfier.ne.0) then
                                write(*,*) "error opening object /input/Igeometry"
                              else
                                if (verbose) then ; write(ounit,*) "successfully opened object at /input/Igeometry" ; endif
                              
                                ! query type of /input/Igeometry object
                                call h5iget_type_f(dset_input_Igeometry, item_type, hdfier)
                                if (hdfier.ne.0) then
                                  write(*,*) "error querying type of object /input/Igeometry"
                                else
                                  if (verbose) then ; write(ounit,*) "successfully queried item type of /input/Igeometry" ; endif
                                
                                  ! check if /input/Igeometry is a Dataset
                                  if (item_type.ne.H5I_DATASET_F) then
                                    write(ounit,*) "error verifying that object /input/Igeometry is a Dataset(",H5I_DATASET_F,"); rather it is ",item_type
                                  else
                                    if (verbose) then ; write(ounit,*) "successfully verified that /input/Igeometry is a Dataset" ; endif
                                    
                                    ! open dataspace of /input/Igeometry
                                    call h5dget_space_f(dset_input_Igeometry, dataspace, hdfier)
                                    if (hdfier.ne.0) then
                                      write(ounit,*) "error getting dataspace of /input/Igeometry"
                                    else
                                    
                                      ! check that the dataspace of /input/Igeometry is scalar
                                      call h5sget_simple_extent_type_f(dataspace, dspace_type, hdfier)
                                      if (hdfier.ne.0) then
                                        write(ounit,*) "error getting type of /input/Igeometry dataspace"
                                      else
                                      
                                        ! check that the type of the /input/Igeometry dataspace is H5S_SCALAR_F
                                        if (dspace_type.eq.H5S_SCALAR_F) then
                                          write(ounit,*) "successfully verified that the type of the /input/Igeometry dataspace is H5S_SCALAR_F"
                                    
                                          ! query datatype of /input/Igeometry
                                          call h5dget_type_f(dset_input_Igeometry, dtype_id, hdfier)
                                          if (hdfier.ne.0) then
                                            write(ounit,*) "error querying datatype of /input/Igeometry"
                                          else
                                            if (verbose) then ; write(ounit,*) "successfully queried datatype of /input/Igeometry" ; endif
                                              
                                            ! convert datatype to native array for comparison
                                            call h5tget_native_type_f(dtype_id, 1, dtype_id_native, hdfier)
                                            if (hdfier.ne.0) then
                                              write(ounit,*) "error converting datatype of /input/Igeometry to native datatype"
                                            else
                                              if (verbose) then ; write(ounit,*) "successfully converted datatype of /input/Igeometry to native datatype" ; endif
                                              
                                              ! call comparison routine for /input/Igeometry datatype
                                              call h5tequal_f(dtype_id_native, H5T_NATIVE_INTEGER, datatypes_equal, hdfier)
                                              if (hdfier.ne.0) then
                                                write(ounit,*) "error comparing datatype of /input/Igeometry to H5T_NATIVE_INTEGER"
                                              else
                                                if (verbose) then ; write(ounit,*) "successfully executed comparison of datatype of /input/Igeometry with H5T_NATIVE_INTEGER" ; endif
                                                
                                                ! verify correct datatype of /input/Igeometry
                                                if (datatypes_equal) then
                                                  if (verbose) then ; write(ounit,*) "successfully checked that datatype of /input/Igeometry is H5T_NATIVE_INTEGER :-)" ; endif
                                                  
                                                  ! read /input/Igeometry Dataset
                                                  call h5dread_f(dset_input_Igeometry, H5T_NATIVE_INTEGER, Igeometry, int((/1/), HSIZE_T), hdfier)
                                                  if (hdfier.ne.0) then
                                                    write(*,*) "error reading Dataset /input/Igeometry"
                                                  else
                                                    if (verbose) then ; write(ounit,'(" successfully read /input/Igeometry from input data: ",i2)') Igeometry ; endif
                                                  endif ! read /input/Igeometry Dataset
                                                else
                                                  write(ounit,*) "ERROR: native datatype of /input/Igeometry should be H5T_NATIVE_INTEGER but is ",dtype_id_native
                                                endif !verify correct datatype of /input/Igeometry
                                              
                                              endif ! call comparison routine for /input/Igeometry datatype
                                            endif ! convert datatype to native array for comparison
                                            
                                            ! close datatype of /input/Igeometry
                                            call h5tclose_f(dtype_id, hdfier)
                                            if (hdfier.ne.0) then
                                              write(ounit,*) "error closing datatype of /input/Igeometry"
                                            elseif (verbose) then
                                              write(ounit,*) "successfully closed datatype of /input/Igeometry"
                                            endif ! close datatype of /input/Igeometry
                                            
                                          endif ! query datatype of /input/Igeometry
                                        else
                                          write(ounit,*) "ERROR: type of dataspace /input/Igeometry is not H5S_SCALAR_F but ",dspace_type
                                        endif ! check that the type of the /input/Igeometry dataspace is H5S_SCALAR_F
                                      endif ! check that the dataspace of /input/Igeometry is scalar
                                      
                                      ! close dataspace of /input/Igeometry
                                      call h5sclose_f(dataspace, hdfier)
                                      if (hdfier.ne.0) then ; write(ounit,*) "error closing dataspace of /input/Igeometry" ; endif
                                    endif ! open dataspace of /input/Igeometry
                                  endif ! check if /input/Igeometry is a Dataset
                                endif ! query type of /input/Igeometry object
                                        
                                call h5oclose_f(dset_input_Igeometry, hdfier)
                                if (hdfier.ne.0) then
                                  write(*,*) "error closing object /input/Igeometry"
                                elseif (verbose) then
                                  write(ounit,*) "successfully closed object /input/Igeometry"
                                endif
                              endif ! try to open /input/Igeometry object
                            else
                              ! /input/Igeometry link present but does not resolve to any object
                              write(ounit,*) "/input/Igeometry link present but does not resolve to any object"
                            endif ! check that there exists an item at the /input/Igeometry link
                          endif ! query existence of object at /input/Igeometry link
                        else
                          ! /input/Igeometry link not present in input file
                          write(ounit,*) " WARNING: /input/Igeometry not found. Using default value ", Igeometry
                        endif ! check if /input/Igeometry link exists
                      endif ! query existence of /input/Igeometry link
                      
                      write(ounit,*) " "
                      
                      ! query existence of /input/Nvol link
                      call h5lexists_f(grp_input, "Nvol", item_exists, hdfier)
                      if (hdfier.ne.0) then
                        write(*,*) "error checking if /input/Nvol link exists"
                      else
                        if (verbose) then ; write(ounit,*) "successfully checked if the /input/Nvol link exists" ; endif
                      
                        ! check that /input/Nvol link exists
                        if (item_exists) then
                          if (verbose) then; write(ounit,*) "successfully checked that the link /input/Nvol exists" ; endif
                         
                          ! query existence of object at /input/Nvol link
                          call h5oexists_by_name_f(grp_input, "Nvol", item_exists, hdfier)
                          if (hdfier.ne.0) then
                            write(*,*) "error checking for presence of object at /input/Nvol"
                          else
                            if (verbose) then ; write(ounit,*) "successfully checked for presence of an object at /input/Nvol link" ; endif
                          
                            ! check that there exists an item at the /input/Nvol link
                            if (item_exists) then
                              if (verbose) then ; write(ounit,*) "successfully checked that there exists an item at /input/Nvol" ; endif
                            
                              ! try to open /input/Nvol object
                              call h5oopen_f(grp_input, "Nvol", dset_input_Nvol, hdfier)
                              if (hdfier.ne.0) then
                                write(*,*) "error opening object /input/Nvol"
                              else
                                if (verbose) then ; write(ounit,*) "successfully opened object at /input/Nvol" ; endif
                              
                                ! query type of /input/Nvol object
                                call h5iget_type_f(dset_input_Nvol, item_type, hdfier)
                                if (hdfier.ne.0) then
                                  write(*,*) "error querying type of object /input/Nvol"
                                else
                                  if (verbose) then ; write(ounit,*) "successfully queried item type of /input/Nvol" ; endif
                                
                                  ! check if /input/Nvol is a Dataset
                                  if (item_type.ne.H5I_DATASET_F) then
                                    write(ounit,*) "error verifying that object /input/Nvol is a Dataset(",H5I_DATASET_F,"); rather it is ",item_type
                                  else
                                    if (verbose) then ; write(ounit,*) "successfully verified that /input/Nvol is a Dataset" ; endif
                                    
                                    ! open dataspace of /input/Nvol
                                    call h5dget_space_f(dset_input_Nvol, dataspace, hdfier)
                                    if (hdfier.ne.0) then
                                      write(ounit,*) "error getting dataspace of /input/Nvol"
                                    else
                                    
                                      ! check that the dataspace of /input/Nvol is scalar
                                      call h5sget_simple_extent_type_f(dataspace, dspace_type, hdfier)
                                      if (hdfier.ne.0) then
                                        write(ounit,*) "error getting type of /input/Nvol dataspace"
                                      else
                                      
                                        ! check that the type of the /input/Nvol dataspace is H5S_SCALAR_F
                                        if (dspace_type.eq.H5S_SCALAR_F) then
                                          write(ounit,*) "successfully verified that the type of the /input/Nvol dataspace is H5S_SCALAR_F"
                                    
                                          ! query datatype of /input/Nvol
                                          call h5dget_type_f(dset_input_Nvol, dtype_id, hdfier)
                                          if (hdfier.ne.0) then
                                            write(ounit,*) "error querying datatype of /input/Nvol"
                                          else
                                            if (verbose) then ; write(ounit,*) "successfully queried datatype of /input/Nvol" ; endif
                                              
                                            ! convert datatype to native array for comparison
                                            call h5tget_native_type_f(dtype_id, 1, dtype_id_native, hdfier)
                                            if (hdfier.ne.0) then
                                              write(ounit,*) "error converting datatype of /input/Nvol to native datatype"
                                            else
                                              if (verbose) then ; write(ounit,*) "successfully converted datatype of /input/Nvol to native datatype" ; endif
                                              
                                              ! call comparison routine for /input/Nvol datatype
                                              call h5tequal_f(dtype_id_native, H5T_NATIVE_INTEGER, datatypes_equal, hdfier)
                                              if (hdfier.ne.0) then
                                                write(ounit,*) "error comparing datatype of /input/Nvol to H5T_NATIVE_INTEGER"
                                              else
                                                if (verbose) then ; write(ounit,*) "successfully executed comparison of datatype of /input/Nvol with H5T_NATIVE_INTEGER" ; endif
                                                
                                                ! verify correct datatype of /input/Nvol
                                                if (datatypes_equal) then
                                                  if (verbose) then ; write(ounit,*) "successfully checked that datatype of /input/Nvol is H5T_NATIVE_INTEGER :-)" ; endif
                                                  
                                                  ! read /input/Nvol Dataset
                                                  call h5dread_f(dset_input_Nvol, H5T_NATIVE_INTEGER, Nvol, int((/1/), HSIZE_T), hdfier)
                                                  if (hdfier.ne.0) then
                                                    write(*,*) "error reading Dataset /input/Nvol"
                                                  else
                                                    if (verbose) then ; write(ounit,'(" successfully read /input/Nvol from input data: ",i2)') Nvol ; endif
                                                  endif ! read /input/Nvol Dataset
                                                else
                                                  write(ounit,*) "ERROR: native datatype of /input/Nvol should be H5T_NATIVE_INTEGER but is ",dtype_id_native
                                                endif !verify correct datatype of /input/Nvol
                                              
                                              endif ! call comparison routine for /input/Nvol datatype
                                            endif ! convert datatype to native array for comparison
                                            
                                            ! close datatype of /input/Nvol
                                            call h5tclose_f(dtype_id, hdfier)
                                            if (hdfier.ne.0) then
                                              write(ounit,*) "error closing datatype of /input/Nvol"
                                            elseif (verbose) then
                                              write(ounit,*) "successfully closed datatype of /input/Nvol"
                                            endif ! close datatype of /input/Nvol
                                          endif ! query datatype of /input/Nvol
                                        else
                                          write(ounit,*) "ERROR: type of dataspace /input/Nvol is not H5S_SCALAR_F but ",dspace_type
                                        endif ! check that the type of the /input/Nvol dataspace is H5S_SCALAR_F
                                      endif ! check that the dataspace of /input/Nvol is scalar
                                      
                                      ! close dataspace of /input/Nvol
                                      call h5sclose_f(dataspace, hdfier)
                                      if (hdfier.ne.0) then ; write(ounit,*) "error closing dataspace of /input/Nvol" ; endif
                                    endif ! open dataspace of /input/Nvol
                                  endif ! check if /input/Nvol is a Dataset
                                endif ! query type of /input/Nvol object
                                        
                                call h5oclose_f(dset_input_Nvol, hdfier)
                                if (hdfier.ne.0) then
                                  write(*,*) "error closing object /input/Nvol"
                                elseif (verbose) then
                                  write(ounit,*) "successfully closed object /input/Nvol"
                                endif
                              endif ! try to open /input/Nvol object
                            else
                              ! /input/Nvol link present but does not resolve to any object
                              write(ounit,*) "/input/Nvol link present but does not resolve to any object"
                            endif ! check that there exists an item at the /input/Nvol link
                          endif ! query existence of object at /input/Nvol link
                        else
                          ! /input/Nvol link not present in input file
                          write(ounit,*) "WARNING: /input/Nvol not found. Using default value ", Nvol
                        endif ! check if /input/Nvol link exists
                      endif ! query existence of /input/Nvol link
                      
                      write(ounit,*) " "
                                            
                      ! query existence of /input/Lrad link
                      call h5lexists_f(grp_input, "Lrad", item_exists, hdfier)
                      if (hdfier.ne.0) then
                        write(*,*) "error checking if /input/Lrad link exists"
                      else
                        if (verbose) then ; write(ounit,*) "successfully checked if the /input/Lrad link exists" ; endif
                      
                        ! check that /input/Lrad link exists
                        if (item_exists) then
                          if (verbose) then; write(ounit,*) "successfully checked that the link /input/Lrad exists" ; endif
                         
                          ! query existence of object at /input/Lrad link
                          call h5oexists_by_name_f(grp_input, "Lrad", item_exists, hdfier)
                          if (hdfier.ne.0) then
                            write(*,*) "error checking for presence of object at /input/Lrad"
                          else
                            if (verbose) then ; write(ounit,*) "successfully checked for presence of an object at /input/Lrad link" ; endif
                          
                            ! check that there exists an item at the /input/Lrad link
                            if (item_exists) then
                              if (verbose) then ; write(ounit,*) "successfully checked that there exists an item at /input/Lrad" ; endif
                            
                              ! try to open /input/Lrad object
                              call h5oopen_f(grp_input, "Lrad", dset_input_Lrad, hdfier)
                              if (hdfier.ne.0) then
                                write(*,*) "error opening object /input/Lrad"
                              else
                                if (verbose) then ; write(ounit,*) "successfully opened object at /input/Lrad" ; endif
                              
                                ! query type of /input/Lrad object
                                call h5iget_type_f(dset_input_Lrad, item_type, hdfier)
                                if (hdfier.ne.0) then
                                  write(*,*) "error querying type of object /input/Lrad"
                                else
                                  if (verbose) then ; write(ounit,*) "successfully queried item type of /input/Lrad" ; endif
                                
                                  ! check if /input/Lrad is a Dataset
                                  if (item_type.ne.H5I_DATASET_F) then
                                    write(ounit,*) "error verifying that object /input/Lrad is a Dataset(",H5I_DATASET_F,"); rather it is ",item_type
                                  else
                                    if (verbose) then ; write(ounit,*) "successfully verified that /input/Lrad is a Dataset" ; endif
                                    
                                    ! open dataspace of /input/Lrad
                                    call h5dget_space_f(dset_input_Lrad, dataspace, hdfier)
                                    if (hdfier.ne.0) then
                                      write(ounit,*) "error getting dataspace of /input/Lrad"
                                    else
                                    
                                      ! check that the dataspace of /input/Lrad is simple
                                      call h5sget_simple_extent_type_f(dataspace, dspace_type, hdfier)
                                      if (hdfier.ne.0) then
                                        write(ounit,*) "error getting type of /input/Lrad dataspace"
                                      else
                                      
                                        ! check that the type of the /input/Lrad dataspace is H5S_SIMPLE_F
                                        if (dspace_type.eq.H5S_SIMPLE_F) then
                                          write(ounit,*) "successfully verified that the type of the /input/Lrad dataspace is H5S_SIMPLE_F"
                                    
                                          ! query rank of /input/Lrad
                                          call h5sget_simple_extent_ndims_f(dataspace, rank, hdfier)
                                          if (hdfier.ne.0) then
                                            write(ounit,*) "error getting rank of /input/Lrad; expected ",rank,", got ",hdfier
                                          else
                                            if (verbose) then ; write(ounit,*) "successfully queried rank of /input/Lrad" ; endif
                                            
                                            ! get the current and maximum dimensions of /input/Lrad
                                            call h5sget_simple_extent_dims_f(dataspace, dims_1, max_dims_1, hdfier)
                                            if (hdfier.ne.rank) then
                                              write(ounit,*) "ERROR: rank mismatch of /input/Lrad; expected ",rank,", but got ",hdfier
                                            else
                                              write(ounit,*) "successfully queried length of /input/Lrad"
                                              
                                              ! check that /input/Lrad has the correct length
                                              if (dims_1(1).ne.Nvol) then ! TODO: Nvol+Lfreeboundary
                                                write(ounit,*) "ERROR: expected length of /input/Lrad is Nvol+Lfreeboundary, but got ", dims_1
                                              else
                                                if (verbose) then ; write(ounit,*) "successfully verified the correct length of /input/Lrad" ; endif
                                                
                                                ! query datatype of /input/Lrad
                                                call h5dget_type_f(dset_input_Lrad, dtype_id, hdfier)
                                                if (hdfier.ne.0) then
                                                  write(ounit,*) "error querying datatype of /input/Lrad"
                                                else
                                                  if (verbose) then ; write(ounit,*) "successfully queried datatype of /input/Lrad" ; endif
                                                    
                                                  ! convert datatype to native array for comparison
                                                  call h5tget_native_type_f(dtype_id, 1, dtype_id_native, hdfier)
                                                  if (hdfier.ne.0) then
                                                    write(ounit,*) "error converting datatype of /input/Lrad to native datatype"
                                                  else
                                                    if (verbose) then ; write(ounit,*) "successfully converted datatype of /input/Lrad to native datatype" ; endif
                                                    
                                                    ! call comparison routine for /input/Lrad datatype
                                                    call h5tequal_f(dtype_id_native, H5T_NATIVE_INTEGER, datatypes_equal, hdfier)
                                                    if (hdfier.ne.0) then
                                                      write(ounit,*) "error comparing datatype of /input/Lrad to H5T_NATIVE_INTEGER"
                                                    else
                                                      if (verbose) then ; write(ounit,*) "successfully executed comparison of datatype of /input/Lrad with H5T_NATIVE_INTEGER" ; endif
                                                      
                                                      ! verify correct datatype of /input/Lrad
                                                      if (datatypes_equal) then
                                                        if (verbose) then ; write(ounit,*) "successfully checked that datatype of /input/Lrad is H5T_NATIVE_INTEGER :-)" ; endif
                                                        
                                                        ! read /input/Lrad Dataset
                                                        call h5dread_f(dset_input_Lrad, H5T_NATIVE_INTEGER, Lrad(1:Nvol), int((/Nvol/), HSIZE_T), hdfier)
                                                        if (hdfier.ne.0) then
                                                          write(*,*) "error reading Dataset /input/Lrad"
                                                        else
                                                          if (verbose) then ; write(ounit,'(" successfully read /input/Lrad from input data: ",i2)') Lrad(1:Nvol) ; endif
                                                        endif ! read /input/Lrad Dataset
                                                      else
                                                        write(ounit,*) "ERROR: native datatype of /input/Lrad should be H5T_NATIVE_INTEGER but is ",dtype_id_native
                                                      endif !verify correct datatype of /input/Lrad
                                                    
                                                    endif ! call comparison routine for /input/Lrad datatype
                                                  endif ! convert datatype to native array for comparison
                                                  
                                                  ! close datatype of /input/Lrad
                                                  call h5tclose_f(dtype_id, hdfier)
                                                  if (hdfier.ne.0) then
                                                    write(ounit,*) "error closing datatype of /input/Lrad"
                                                  elseif (verbose) then
                                                    write(ounit,*) "successfully closed datatype of /input/Lrad"
                                                  endif ! close datatype of /input/Lrad
                                                endif ! query datatype of /input/Lrad
                                              endif ! check that /input/Lrad has the correct length
                                            endif ! get the current and maximum dimensions of /input/Lrad
                                          endif ! query rank of /input/Lrad
                                        else
                                          write(ounit,*) "ERROR: type of dataspace /input/Lrad is not H5S_SIMPLE_F but ",dspace_type
                                        endif ! check that the type of the /input/Lrad dataspace is H5S_SIMPLE_F
                                      endif ! check that the dataspace of /input/Lrad is scalar
                                      
                                      ! close dataspace of /input/Lrad
                                      call h5sclose_f(dataspace, hdfier)
                                      if (hdfier.ne.0) then ; write(ounit,*) "error closing dataspace of /input/Lrad" ; endif
                                    endif ! open dataspace of /input/Lrad
                                  endif ! query type of /input/Lrad object
                                endif ! check if /input/Lrad is a Dataset
                                        
                                call h5oclose_f(dset_input_Lrad, hdfier)
                                if (hdfier.ne.0) then
                                  write(*,*) "error closing object /input/Lrad"
                                elseif (verbose) then
                                  write(ounit,*) "successfully closed object /input/Lrad"
                                endif
                              endif ! try to open /input/Lrad object
                            else
                              ! /input/Lrad link present but does not resolve to any object
                              write(ounit,*) "/input/Lrad link present but does not resolve to any object"
                            endif ! check that there exists an item at the /input/Lrad link
                          endif ! query existence of object at /input/Lrad link
                        else
                          ! /input/Lrad link not present in input file
                          write(ounit,*) "WARNING: /input/Lrad not found. Using default value ", Lrad(1:Nvol)
                        endif ! check if /input/Lrad link exists
                      endif ! query existence of /input/Lrad link
                                            
                      write(ounit,*) " "
                      
                      ! query existence of /input/phiedge link
                      call h5lexists_f(grp_input, "phiedge", item_exists, hdfier)
                      if (hdfier.ne.0) then
                        write(*,*) "error checking if /input/phiedge link exists"
                      else
                        if (verbose) then ; write(ounit,*) "successfully checked if the /input/phiedge link exists" ; endif
                      
                        ! check that /input/phiedge link exists
                        if (item_exists) then
                          if (verbose) then; write(ounit,*) "successfully checked that the link /input/phiedge exists" ; endif
                         
                          ! query existence of object at /input/phiedge link
                          call h5oexists_by_name_f(grp_input, "phiedge", item_exists, hdfier)
                          if (hdfier.ne.0) then
                            write(*,*) "error checking for presence of object at /input/phiedge"
                          else
                            if (verbose) then ; write(ounit,*) "successfully checked for presence of an object at /input/phiedge link" ; endif
                          
                            ! check that there exists an item at the /input/phiedge link
                            if (item_exists) then
                              if (verbose) then ; write(ounit,*) "successfully checked that there exists an item at /input/phiedge" ; endif
                            
                              ! try to open /input/phiedge object
                              call h5oopen_f(grp_input, "phiedge", dset_input_phiedge, hdfier)
                              if (hdfier.ne.0) then
                                write(*,*) "error opening object /input/phiedge"
                              else
                                if (verbose) then ; write(ounit,*) "successfully opened object at /input/phiedge" ; endif
                              
                                ! query type of /input/phiedge object
                                call h5iget_type_f(dset_input_phiedge, item_type, hdfier)
                                if (hdfier.ne.0) then
                                  write(*,*) "error querying type of object /input/phiedge"
                                else
                                  if (verbose) then ; write(ounit,*) "successfully queried item type of /input/phiedge" ; endif
                                
                                  ! check if /input/phiedge is a Dataset
                                  if (item_type.ne.H5I_DATASET_F) then
                                    write(ounit,*) "error verifying that object /input/phiedge is a Dataset(",H5I_DATASET_F,"); rather it is ",item_type
                                  else
                                    if (verbose) then ; write(ounit,*) "successfully verified that /input/phiedge is a Dataset" ; endif
                                    
                                    ! open dataspace of /input/phiedge
                                    call h5dget_space_f(dset_input_phiedge, dataspace, hdfier)
                                    if (hdfier.ne.0) then
                                      write(ounit,*) "error getting dataspace of /input/phiedge"
                                    else
                                    
                                      ! check that the dataspace of /input/phiedge is scalar
                                      call h5sget_simple_extent_type_f(dataspace, dspace_type, hdfier)
                                      if (hdfier.ne.0) then
                                        write(ounit,*) "error getting type of /input/phiedge dataspace"
                                      else
                                      
                                        ! check that the type of the /input/phiedge dataspace is H5S_SCALAR_F
                                        if (dspace_type.eq.H5S_SCALAR_F) then
                                          write(ounit,*) "successfully verified that the type of the /input/phiedge dataspace is H5S_SCALAR_F"
                                          
                                          ! query datatype of /input/phiedge
                                          call h5dget_type_f(dset_input_phiedge, dtype_id, hdfier)
                                          if (hdfier.ne.0) then
                                            write(ounit,*) "error querying datatype of /input/phiedge"
                                          else
                                            if (verbose) then ; write(ounit,*) "successfully queried datatype of /input/phiedge" ; endif
                                              
                                            ! convert datatype to native array for comparison
                                            call h5tget_native_type_f(dtype_id, 1, dtype_id_native, hdfier)
                                            if (hdfier.ne.0) then
                                              write(ounit,*) "error converting datatype of /input/phiedge to native datatype"
                                            else
                                              if (verbose) then ; write(ounit,*) "successfully converted datatype of /input/phiedge to native datatype" ; endif
                                              
                                              ! call comparison routine for /input/phiedge datatype
                                              call h5tequal_f(dtype_id_native, H5T_NATIVE_DOUBLE, datatypes_equal, hdfier)
                                              if (hdfier.ne.0) then
                                                write(ounit,*) "error comparing datatype of /input/phiedge to H5T_NATIVE_DOUBLE"
                                              else
                                                if (verbose) then ; write(ounit,*) "successfully executed comparison of datatype of /input/phiedge with H5T_NATIVE_DOUBLE" ; endif
                                                
                                                ! verify correct datatype of /input/phiedge
                                                if (datatypes_equal) then
                                                  if (verbose) then ; write(ounit,*) "successfully checked that datatype of /input/phiedge is H5T_NATIVE_DOUBLE :-)" ; endif
                                                  
                                                  ! read /input/phiedge Dataset
                                                  call h5dread_f(dset_input_phiedge, H5T_NATIVE_DOUBLE, phiedge, int((/1/), HSIZE_T), hdfier)
                                                  if (hdfier.ne.0) then
                                                    write(*,*) "error reading Dataset /input/phiedge"
                                                  else
                                                    if (verbose) then ; write(ounit,'(" successfully read /input/phiedge from input data: ",f4.2)') phiedge ; endif
                                                  endif ! read /input/phiedge Dataset
                                                else
                                                  write(ounit,*) "ERROR: native datatype of /input/phiedge should be H5T_NATIVE_DOUBLE but is ",dtype_id_native
                                                endif !verify correct datatype of /input/phiedge
                                              
                                              endif ! call comparison routine for /input/phiedge datatype
                                            endif ! convert datatype to native array for comparison
                                            
                                            ! close datatype of /input/phiedge
                                            call h5tclose_f(dtype_id, hdfier)
                                            if (hdfier.ne.0) then
                                              write(ounit,*) "error closing datatype of /input/phiedge"
                                            elseif (verbose) then
                                              write(ounit,*) "successfully closed datatype of /input/phiedge"
                                            endif ! close datatype of /input/phiedge
                                          endif ! query datatype of /input/phiedge
                                        else
                                          write(ounit,*) "ERROR: type of dataspace /input/phiedge is not H5S_SCALAR_F but ",dspace_type
                                        endif ! check that the type of the /input/phiedge dataspace is H5S_SCALAR_F
                                      
                                      endif ! check that the dataspace of /input/phiedge is scalar
                                      
                                      ! close dataspace of /input/phiedge
                                      call h5sclose_f(dataspace, hdfier)
                                      if (hdfier.ne.0) then ; write(ounit,*) "error closing dataspace of /input/phiedge" ; endif
                                    endif ! open dataspace of /input/phiedge
                                  endif ! check if /input/phiedge is a Dataset
                                endif ! query type of /input/phiedge object
                                
                                ! close /input/phiedge object
                                call h5oclose_f(dset_input_phiedge, hdfier)
                                if (hdfier.ne.0) then
                                  write(*,*) "error closing object /input/phiedge"
                                elseif (verbose) then
                                  write(ounit,*) "successfully closed object /input/phiedge"
                                endif ! close /input/phiedge object
                              endif ! try to open /input/phiedge object
                            else
                              ! /input/phiedge link present but does not resolve to any object
                              write(ounit,*) "/input/phiedge link present but does not resolve to any object"
                            endif ! check that there exists an item at the /input/phiedge link
                          endif ! query existence of object at /input/phiedge link
                        else
                          ! /input/phiedge link not present in input file
                          write(ounit,*) "WARNING: /input/phiedge not found. Using default value ", phiedge
                        endif ! check if /input/phiedge link exists
                      endif ! query existence of /input/phiedge link
                      
                      write(ounit,*) " "
                      
                    endif ! check if /input is a Group
                  endif ! query type of /input object
                          
                  call h5oclose_f(grp_input, hdfier)
                  if (hdfier.ne.0) then
                    write(*,*) "error closing object /input"
                  elseif (verbose) then
                    write(ounit,*) "successfully closed object /input"
                  endif
                endif ! try to open /input object
              else
                ! /input link present but does not resolve to any object
                write(ounit,*) "/input link present but does not resolve to any object"
              endif ! check that there exists an item at the /input link
            endif ! query existence of object at /input link
          else
            ! /input link not present in input file
            write(ounit,*) "WARNING: /input not found; cannot read input data from the given file"
          endif ! check if /input link exists
        endif ! query existence of /input link
        
        write(ounit,*) " "
        
        ! close input file
        call h5fclose_f(file_id, hdfier)
        if (hdfier.ne.0) then ; write(*,*) "error closing input HDF5 file" ; endif
      
      endif ! open input file
    else
      write(ounit,*) "the input file '",filename_h5,"' does not exist"
      call exit(1) ! "generic fail" return code
    endif
  end subroutine readin_h5
end module inputlist

module allglobal
  use inputlist
  use spec_version
  
  implicit none
  
end module allglobal

module sphdf5
  
  contains
  
  ! attach a description attribute to the given HDF5 item (Group or Dataset)
  subroutine attach_description(item_id, description)
    use hdf5
    use fileunits
    use allglobal
    implicit none
    integer(hid_t),   intent(in) :: item_id
    character(len=*), intent(in) :: description
    
    character(LEN=15), parameter       :: attr_name = "description" !< name of the descriptive Attribute to be attached to each Dataset and Group
    integer(HID_T)                     :: attr_id                   !< Attribute identifier
    integer(HID_T)                     :: attr_dspace               !< Attribute Dataspace identifier
    integer(HID_T)                     :: attr_dtype                !< Attribute Dataspace identifier
    integer(HSIZE_T)                   :: attr_len                  !< Length of the attribute string
    
    integer :: hdfier   !< error flag for HDF5 library
    
    ! compute length of current attribute
    attr_len = len(description)
    
    ! create dataspace for descriptive Attribute; only one description per item, so H5S_SCALAR_F
    call h5screate_f(H5S_SCALAR_F, attr_dspace, hdfier)
    if (hdfier.ne.0) then
      write(ounit,*) "error creating dataspace for attribute"
    else
      if (verbose) then ; write(ounit,*) "successfully created dataspace for attribute" ; endif
      
      ! copy base type for attribute data (String)
      call h5tcopy_f(H5T_C_S1, attr_dtype, hdfier)
      if (hdfier.ne.0) then
        write(ounit,*) "error copying H5T_C_S1 to attr_dtype"
      else
        if (verbose) then ; write(ounit,*) "successfully copyied H5T_C_S1 to attr_dtype" ; endif
      
        ! set charset used for attribute string to be UTF-8; e.g. for Poincaré
        call h5tset_cset_f(attr_dtype, H5T_CSET_UTF8_F, hdfier)
        if (hdfier.ne.0) then
          write(ounit,*) "error setting charset used for attribute string to be UTF-8"
        else
          if (verbose) then ; write(ounit,*) "successfully set charset used for attribute string to be UTF-8" ; endif
        
          ! set attribute string length
          call h5tset_size_f(attr_dtype, attr_len, hdfier)
          if (hdfier.ne.0) then
            write(ounit,*) "error setting attribute string length"
          else
            if (verbose) then ; write(ounit,*) "successfully set attribute string length" ; endif
          
            ! create Attribute 'description' attached to Group or Dataset identified by item_id
            call h5acreate_f(item_id, attr_name, attr_dtype, attr_dspace, attr_id, hdfier)
            if (hdfier.ne.0) then
              write(ounit,*) "error creating 'description' attribute"
            else
              if (verbose) then ; write(ounit,*) "successfully created 'description' attribute" ; endif
            
              ! write attribute data
              call h5awrite_f(attr_id, attr_dtype, description, int((/attr_len/),size_t), hdfier)
              if (hdfier.ne.0) then
                write(ounit,*) "error writing 'description' data"
              elseif (verbose) then ; write(ounit,*) "successfully wrote 'description' data"
              endif ! write attribute data  
            endif ! create Attribute 'description' attached to Group or Dataset identified by item_id
          endif ! set attribute string length
        endif ! set charset used for attribute string to be UTF-8
      endif ! copy base type for attribute data (String)
      
      ! close dataspace of attribute
      call h5sclose_f(attr_dspace, hdfier)
      if (hdfier.ne.0) then
        write(ounit,*) "error closing dataspace of attribute"
      elseif (verbose) then ; write(ounit,*) "successfully closed dataspace of attribute"
      endif ! close dataspace of attribute
    
    endif ! create dataspace for descriptive Attribute
    
  end subroutine attach_description
  
  subroutine mirror_input_to_outfile(filename_h5)
    use hdf5
    use fileunits
    use allglobal
    
    implicit none
    
    character(LEN=15), parameter       :: attr_name = "description" !< name of the descriptive Attribute to be attached to each Dataset and Group
    integer(HID_T)                     :: attr_id                   ! Attribute identifier
    integer(HID_T)                     :: attr_dspace               ! Attribute Dataspace identifier
    integer(HID_T)                     :: attr_dtype                ! Attribute Dataspace identifier
    integer(SIZE_T)                    :: attr_len                  ! Length of the attribute string
    character(len=:), allocatable      :: description                 ! Attribute data == description content
    
    character(len=*),intent(in) :: filename_h5 ! name of output file
    logical :: file_exists
    integer :: hdfier ! error flag for HDF5 library
    
    integer(hid_t)              :: dataspace         ! dataspace used to write stuff
    
    integer(hid_t)              :: file_id
    integer(hid_t)              :: dset_version
    integer(hid_t)              :: grp_input
    integer(hid_t)              :: dset_input_Igeometry
    integer(hid_t)              :: dset_input_Nvol
    integer(hid_t)              :: dset_input_Lrad
    integer(hid_t)              :: dset_input_phiedge
    
    inquire( file=filename_h5, exist=file_exists ) ! check if file exists
    if (file_exists) then
      write(ounit,*) "output HDF5 file '",filename_h5,"' already exists"
      call exit(1) ! TODO: Loverwrite?
    endif
    
    write(ounit,*) "writing input data into output HDF5 file '",filename_h5,"'"
    
    ! create and open the output file; will be overwritten if already existent
    call h5fcreate_f(filename_h5, H5F_ACC_TRUNC_F, file_id, hdfier)
    if (hdfier.ne.0) then
      write(ounit,*) "error creating or opening output file"
    else
      if (verbose) then ; write(ounit,*) "successfully opened output file" ; endif
      
      
      
      
      write(ounit,*) " "
      
      ! create dataspace for /version
      call h5screate_simple_f(1, int((/3/),hsize_t), dataspace, hdfier)
      if (hdfier.ne.0) then
        write(ounit,*) "error creating dataspace for /version"
      else
        if (verbose) then ; write(ounit,*) "successfully created dataspace for /version" ; endif
        
        ! create dataset for /version
        call h5dcreate_f(file_id, "version", H5T_NATIVE_INTEGER, dataspace, dset_version, hdfier)
        if (hdfier.ne.0) then
          write(ounit,*) "error creating dataspace for /version"
        else
          if (verbose) then ; write(ounit,*) "successfully created dataset for /version" ; endif
          
          ! write version number to output file
          call h5dwrite_f(dset_version, H5T_NATIVE_INTEGER, version, int((/3/),hsize_t), hdfier)
          if (hdfier.ne.0) then
            write(ounit,*) "error writing data of /version Dataset"
          else
            if (verbose) then ; write(ounit,*) "successfully wrote data of /version Dataset" ; endif
            
            ! attribute content is description of /version
            description = "version of SPEC"
            
            ! attach description to /version Dataset
            call attach_description(dset_version, description)
          
          endif ! write version number to output file
          
          ! close dataset of /version
          call h5dclose_f(dset_version, hdfier)
          if (hdfier.ne.0) then
            write(ounit,*) "error closing /version Dataset"
          elseif (verbose) then ; write(ounit,*) "successfully closed /version Datset"
          endif ! close dataset of /version
        endif ! create dataset for /version
        
        ! close dataspace for /version
        call h5sclose_f(dataspace, hdfier)
        if (hdfier.ne.0) then
          write(ounit,*) "error closing dataspace for /version"
        elseif (verbose) then ; write(ounit,*) "successfully closed dataspace for /version"
        endif ! close dataspace for /version
      endif ! create dataspace for /version
      
      write(ounit,*) " "
      
      ! create /input group
      call h5gcreate_f(file_id, "input", grp_input, hdfier)
      if (hdfier.ne.0) then
        write(ounit,*) "error creating /input Group"
      else
        if (verbose) then ; write(ounit,*) "successfully created /input Group" ; endif
        
        ! attribute content is description of /input group
        description = "group for mirrored input data"
        
        ! attach description to /input Group
        call attach_description(grp_input, description)
        
        write(ounit,*) " "
        
        ! create dataspace for /input/Igeometry
        call h5screate_f(H5S_SCALAR_F, dataspace, hdfier)
        if (hdfier.ne.0) then
          write(ounit,*) "error creating dataspace for /input/Igeometry"
        else
          if (verbose) then ; write(ounit,*) "successfully created dataspace for /input/Igeometry" ; endif
          
          ! create dataset for /input/Igeometry
          call h5dcreate_f(grp_input, "Igeometry", H5T_NATIVE_INTEGER, dataspace, dset_input_Igeometry, hdfier)
          if (hdfier.ne.0) then
            write(ounit,*) "error creating dataspace for /input/Igeometry"
          else
            if (verbose) then ; write(ounit,*) "successfully created dataset for /input/Igeometry" ; endif
            
            ! write /input/Igeometry to output file
            call h5dwrite_f(dset_input_Igeometry, H5T_NATIVE_INTEGER, Igeometry, int((/1/),hsize_t), hdfier)
            if (hdfier.ne.0) then
              write(ounit,*) "error writing data of /input/Igeometry Dataset"
            else
              if (verbose) then ; write(ounit,*) "successfully wrote data of /input/Igeometry Dataset" ; endif
              
              ! attribute content is description of /input/Igeometry
              description = "selects Cartesian, cylindrical or toroidal geometry"//new_line('A') &
              &          // "<ul>"//new_line('A') &
              &          // "<li> \c Igeometry=1 : Cartesian; geometry determined by \f$R\f$              </li>"//new_line('A') &
              &          // "<li> \c Igeometry=2 : cylindrical; geometry determined by \f$R\f$            </li>"//new_line('A') &
              &          // "<li> \c Igeometry=3 : toroidal; geometry determined by \f$R\f$ *and* \f$Z\f$ </li>"//new_line('A') &
              &          // "</ul>"
              
              ! attach description to /input/Igeometry Dataset
              call attach_description(dset_input_Igeometry, description)
            
            endif ! write /input/Igeometry to output file
            
            ! close dataset of /input/Igeometry
            call h5dclose_f(dset_input_Igeometry, hdfier)
            if (hdfier.ne.0) then
              write(ounit,*) "error closing /input/Igeometry Dataset"
            elseif (verbose) then ; write(ounit,*) "successfully closed /input/Igeometry Datset"
            endif ! close dataset of /input/Igeometry
          endif ! create dataset for /input/Igeometry
          
          ! close dataspace for /input/Igeometry
          call h5sclose_f(dataspace, hdfier)
          if (hdfier.ne.0) then
            write(ounit,*) "error closing dataspace for /input/Igeometry"
          elseif (verbose) then ; write(ounit,*) "successfully closed dataspace for /input/Igeometry"
          endif ! close dataspace for /input/Igeometry
        endif ! create dataspace for /input/Igeometry
        
        write(ounit,*) " "
        
        ! create dataspace for /input/Nvol
        call h5screate_f(H5S_SCALAR_F, dataspace, hdfier)
        if (hdfier.ne.0) then
          write(ounit,*) "error creating dataspace for /input/Nvol"
        else
          if (verbose) then ; write(ounit,*) "successfully created dataspace for /input/Nvol" ; endif
          
          ! create dataset for /input/Nvol
          call h5dcreate_f(grp_input, "Nvol", H5T_NATIVE_INTEGER, dataspace, dset_input_Nvol, hdfier)
          if (hdfier.ne.0) then
            write(ounit,*) "error creating dataspace for /input/Nvol"
          else
            if (verbose) then ; write(ounit,*) "successfully created dataset for /input/Nvol" ; endif
            
            ! write /input/Nvol to output file
            call h5dwrite_f(dset_input_Nvol, H5T_NATIVE_INTEGER, Nvol, int((/1/),hsize_t), hdfier)
            if (hdfier.ne.0) then
              write(ounit,*) "error writing data of /input/Nvol Dataset"
            else
              if (verbose) then ; write(ounit,*) "successfully wrote data of /input/Nvol Dataset" ; endif
              
              ! attribute content is description of /input/Nvol
              description = "number of volumes"//new_line('A') &
              &          // "<ul>"//new_line('A') &
              &          // "<li> each volume \f${\cal V}_l\f$ is bounded by the \f${\cal I}_{l-1}\f$ and \f${\cal I}_{l}\f$ interfaces </li>"//new_line('A') &
              &          // "<li> note that in cylindrical or toroidal geometry, \f${\cal I}_{0}\f$ is the degenerate coordinate axis   </li>"//new_line('A') &
              &          // "<li> constraint: \c Nvol<=MNvol                                                                            </li>"//new_line('A') &
              &          // "</ul>"
              
              ! attach description to /input/Nvol Dataset
              call attach_description(dset_input_Nvol, description)
            
            endif ! write /input/Nvol to output file
            
            ! close dataset of /input/Nvol
            call h5dclose_f(dset_input_Nvol, hdfier)
            if (hdfier.ne.0) then
              write(ounit,*) "error closing /input/Nvol Dataset"
            elseif (verbose) then ; write(ounit,*) "successfully closed /input/Nvol Datset"
            endif ! close dataset of /input/Nvol
          endif ! create dataset for /input/Nvol
          
          ! close dataspace for /input/Nvol
          call h5sclose_f(dataspace, hdfier)
          if (hdfier.ne.0) then
            write(ounit,*) "error closing dataspace for /input/Nvol"
          elseif (verbose) then ; write(ounit,*) "successfully closed dataspace for /input/Nvol"
          endif ! close dataspace for /input/Nvol
        endif ! create dataspace for /input/Nvol
        
        write(ounit,*) " "
        
        ! create dataspace for /input/Lrad
        call h5screate_simple_f(1, int((/Nvol/),hsize_t), dataspace, hdfier) ! TODO: Nvol+Lfreeboundary
        if (hdfier.ne.0) then
          write(ounit,*) "error creating dataspace for /input/Lrad"
        else
          if (verbose) then ; write(ounit,*) "successfully created dataspace for /input/Lrad" ; endif
          
          ! create dataset for /input/Lrad
          call h5dcreate_f(grp_input, "Lrad", H5T_NATIVE_INTEGER, dataspace, dset_input_Lrad, hdfier)
          if (hdfier.ne.0) then
            write(ounit,*) "error creating dataspace for /input/Lrad"
          else
            if (verbose) then ; write(ounit,*) "successfully created dataset for /input/Lrad" ; endif
            
            ! write /input/Lrad to output file
            call h5dwrite_f(dset_input_Lrad, H5T_NATIVE_INTEGER, Lrad(1:Nvol), int((/Nvol/),hsize_t), hdfier) ! TODO: Nvol+Lfreeboundary
            if (hdfier.ne.0) then
              write(ounit,*) "error writing data of /input/Lrad Dataset"
            else
              if (verbose) then ; write(ounit,*) "successfully wrote data of /input/Lrad Dataset" ; endif
              
              ! attribute content is description of /input/Lrad
              description = "Chebyshev resolution in each volume"//new_line('A') &
              &          // "<ul>"//new_line('A') &
              &          // "<li> constraint : \c Lrad(1:Mvol) >= 2 </li>"//new_line('A') &
              &          // "</ul>"
              
              ! attach description to /input/Lrad Dataset
              call attach_description(dset_input_Lrad, description)
            
            endif ! write /input/Lrad to output file
            
            ! close dataset of /input/Lrad
            call h5dclose_f(dset_input_Lrad, hdfier)
            if (hdfier.ne.0) then
              write(ounit,*) "error closing /input/Lrad Dataset"
            elseif (verbose) then ; write(ounit,*) "successfully closed /input/Lrad Datset"
            endif ! close dataset of /input/Lrad
          endif ! create dataset for /input/Lrad
          
          ! close dataspace for /input/Lrad
          call h5sclose_f(dataspace, hdfier)
          if (hdfier.ne.0) then
            write(ounit,*) "error closing dataspace for /input/Lrad"
          elseif (verbose) then ; write(ounit,*) "successfully closed dataspace for /input/Lrad"
          endif ! close dataspace for /input/Lrad
        endif ! create dataspace for /input/Lrad
        
        write(ounit,*) " "
        
        ! create dataspace for /input/phiedge
        call h5screate_f(H5S_SCALAR_F, dataspace, hdfier)
        if (hdfier.ne.0) then
          write(ounit,*) "error creating dataspace for /input/phiedge"
        else
          if (verbose) then ; write(ounit,*) "successfully created dataspace for /input/phiedge" ; endif
          
          ! create dataset for /input/phiedge
          call h5dcreate_f(grp_input, "phiedge", H5T_NATIVE_DOUBLE, dataspace, dset_input_phiedge, hdfier)
          if (hdfier.ne.0) then
            write(ounit,*) "error creating dataspace for /input/phiedge"
          else
            if (verbose) then ; write(ounit,*) "successfully created dataset for /input/phiedge" ; endif
            
            ! write /input/phiedge to output file
            call h5dwrite_f(dset_input_phiedge, H5T_NATIVE_DOUBLE, phiedge, int((/1/),hsize_t), hdfier)
            if (hdfier.ne.0) then
              write(ounit,*) "error writing data of /input/phiedge Dataset"
            else
              if (verbose) then ; write(ounit,*) "successfully wrote data of /input/phiedge Dataset" ; endif
              
              ! attribute content is description of /input/phiedge
              description = "total enclosed toroidal magnetic flux"
              
              ! attach description to /input/phiedge Dataset
              call attach_description(dset_input_phiedge, description)
            
            endif ! write /input/phiedge to output file
            
            ! close dataset of /input/phiedge
            call h5dclose_f(dset_input_phiedge, hdfier)
            if (hdfier.ne.0) then
              write(ounit,*) "error closing /input/phiedge Dataset"
            elseif (verbose) then ; write(ounit,*) "successfully closed /input/phiedge Datset"
            endif ! close dataset of /input/phiedge
          endif ! create dataset for /input/phiedge
          
          ! close dataspace for /input/phiedge
          call h5sclose_f(dataspace, hdfier)
          if (hdfier.ne.0) then
            write(ounit,*) "error closing dataspace for /input/phiedge"
          elseif (verbose) then ; write(ounit,*) "successfully closed dataspace for /input/phiedge"
          endif ! close dataspace for /input/phiedge
        endif ! create dataspace for /input/phiedge
        
        write(ounit,*) " "
      
      endif ! create /input group
      
      ! close /input group
      call h5gclose_f(grp_input, hdfier)
      if (hdfier.ne.0) then
        write(ounit,*) "error closing /input Group"
      elseif (verbose) then ; write(ounit,*) "successfully closed /input Group"
      endif ! close /input group
      
      write(ounit,*) " "
      
      ! finally, close the output file
      call h5fclose_f(file_id, hdfier)
      if (hdfier.ne.0) then
        write(ounit,*) "error closing output file"
      elseif (verbose) then ; write(ounit,*) "successfully closed output file"
      endif ! finally, close the output file
      
    endif ! create and open the output file
  
  end subroutine mirror_input_to_outfile

end module sphdf5
