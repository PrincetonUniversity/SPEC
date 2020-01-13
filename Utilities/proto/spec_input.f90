!> @file spec_input.f90

module fileunits

  implicit none

  integer :: iunit = 10 !< input namelists text file
  integer :: ounit =  6 !< screen output

end module fileunits

module spec_version
  implicit none
  double precision, parameter :: version = 2.10 !< version of SPEC
end module spec_version

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
  
  namelist/physicslist/ &
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
    
    character(len=*),intent(in) :: filename_h5 ! name of input file
    
    logical                     :: Lsph5exist !< flag to indicate that the file \c ext.sp.h5 exists
    integer                     :: hdfier ! error flag for HDF5 library calls
    integer(hid_t)              :: file_id, dtype_id, dtype_id_native
    integer(hid_t)              :: dataspace         ! dataspace used to query Dataset size
    integer(hsize_t)            :: dims_1(1)         ! current dimensions of rank-1 Dataset
    integer(hsize_t)            :: max_dims_1(1)     ! maximum dimensions of rank-1 Dataset
    logical                     :: item_exists, datatypes_equal
    integer                     :: item_type, dspace_type
    
    integer(hid_t)              :: dset_version
    integer(hid_t)              :: grp_input
    integer(hid_t)              :: dset_input_Igeometry
    integer(hid_t)              :: dset_input_Nvol
    integer(hid_t)              :: dset_input_Lrad
    integer(hid_t)              :: dset_input_phiedge
    
    double precision            :: input_version ! version of input data
    
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
                      
                        ! check that the dataspace of /version is scalar
                        call h5sget_simple_extent_type_f(dataspace, dspace_type, hdfier)
                        if (hdfier.ne.0) then
                          write(ounit,*) "error getting type of /version dataspace"
                        else
                        
                          ! check that the type of the /version dataspace is H5S_SCALAR_F
                          if (dspace_type.eq.H5S_SCALAR_F) then
                            write(ounit,*) "successfully verified that the type of the /version dataspace is H5S_SCALAR_F"
                            
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
                                call h5tequal_f(dtype_id_native, H5T_NATIVE_DOUBLE, datatypes_equal, hdfier)
                                if (hdfier.ne.0) then
                                  write(ounit,*) "error comparing datatype of /version to H5T_NATIVE_DOUBLE"
                                else
                                  if (verbose) then ; write(ounit,*) "successfully executed comparison of datatype of /version with H5T_NATIVE_DOUBLE" ; endif
                                  
                                  ! verify correct datatype of /version
                                  if (datatypes_equal) then
                                    if (verbose) then ; write(ounit,*) "successfully checked that datatype of /version is H5T_NATIVE_DOUBLE :-)" ; endif
                                    
                                    ! read version Dataset
                                    call h5dread_f(dset_version, H5T_NATIVE_DOUBLE, input_version, int((/1/), HSIZE_T), hdfier)
                                    if (hdfier.ne.0) then
                                      write(*,*) "error reading Dataset /version"
                                    else
                                      if (verbose) then ; write(ounit,'(" successfully read /version from input data: ",f4.2)') input_version ; endif
                                      
                                      ! verify that version number in input file is less than or equal to the current version of SPEC
                                      if (input_version.le.version) then
                                        if (verbose) then
                                          write(ounit,'(" The version of the input data (",f4.2,") is compatible with SPEC ",f4.2)') &
                                          & input_version, version
                                        endif
                                      else
                                        write(ounit,'(" WARNING: the version of the input data (",f4.2,") may be incompatible with SPEC ",f4.2)') &
                                        & input_version, version
                                      endif ! verify that version number in input file is less than or equal to the current version of SPEC
                                    endif ! read version Dataset
                                  else
                                    write(ounit,*) "ERROR: native datatype of /version should be H5T_NATIVE_DOUBLE but is ",dtype_id_native
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
                          
                          else
                            write(ounit,*) "ERROR: type of dataspace /version is not H5S_SCALAR_F but ",dspace_type
                          endif ! check that the type of the /version dataspace is H5S_SCALAR_F
                        
                        endif ! check that the dataspace of /version is scalar
                        
                        ! close dataspace of /version
                        call h5sclose_f(dataspace, hdfier)
                        if (hdfier.ne.0) then ; write(ounit,*) "error closing dataspace of /version" ; endif

                      endif ! open dataspace of /version
                      
                    endif ! check if /version is a Dataset
                    
                  endif ! query type of /version object
                          
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
            write(ounit,'(" WARNING: /version not found; cannot check if the given input file works for SPEC ",f4.2)')&
            & version
          endif ! check if /version link exists
        endif ! query existence of /version link
        
        
        
        
        
        
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
                          write(ounit,'(" WARNING: /input/Igeometry not found; cannot check if the given input file works for SPEC ",f4.2)')&
                          & version
                        endif ! check if /input/Igeometry link exists
                      endif ! query existence of /input/Igeometry link
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
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
                          write(ounit,'(" WARNING: /input/Nvol not found; cannot check if the given input file works for SPEC ",f4.2)')&
                          & version
                        endif ! check if /input/Nvol link exists
                      endif ! query existence of /input/Nvol link
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      ! /input/Lrad --> Lrad(1:dims_1(1)); rank=1
                      call h5dopen_f(grp_input, "Lrad", dset_input_Lrad, hdfier)
                      if (hdfier.ne.0) then
                        write(*,*) "error opening Dataset /input/Lrad"
                      else
                        ! open dataspace to get current state of Dataset
                        call h5dget_space_f(dset_input_Lrad, dataspace, hdfier)
                        if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for Dataset /input/Lrad" ; endif
                        
                        ! get current size of Dataset
                        call h5sget_simple_extent_dims_f(dataspace, dims_1, max_dims_1, hdfier)
                        if (hdfier.ne.1) then ; write(*,*) "unexpected rank of Dataset /input/Lrad: ",hdfier," .ne. 1" ; endif
                        
                        ! close dataspace after it has been used to query the size of the variable
                        call h5sclose_f(dataspace, hdfier)
                        if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for Dataset /input/Lrad" ; endif
                        
                        call h5dread_f(dset_input_Lrad, H5T_NATIVE_INTEGER, Lrad(1:dims_1(1)), dims_1, hdfier)
                        if (hdfier.ne.0) then ; write(*,*) "error reading Dataset /input/Lrad" ; endif
                        
                        call h5dclose_f(dset_input_Lrad, hdfier)
                        if (hdfier.ne.0) then ; write(*,*) "error closing Dataset /input/Lrad" ; endif
                      endif
                      
                      ! /input/phiedge --> phiedge
                      call h5dopen_f(grp_input, "phiedge", dset_input_phiedge, hdfier)
                      if (hdfier.ne.0) then
                        write(*,*) "error opening Dataset /input/phiedge"
                      else
                        call h5dread_f(dset_input_phiedge, H5T_NATIVE_DOUBLE, phiedge, int((/1/), HSIZE_T), hdfier)
                        if (hdfier.ne.0) then ; write(*,*) "error reading Dataset /input/phiedge" ; endif
                        call h5dclose_f(dset_input_phiedge, hdfier)
                        if (hdfier.ne.0) then ; write(*,*) "error closing Dataset /input/phiedge" ; endif
                      endif
        
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
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
  
  subroutine mirror_input_to_outfile(filename_h5)
    use hdf5
    use fileunits
    use allglobal
    
    implicit none
    
    character(len=*),intent(in) :: filename_h5 ! name of output file
    logical :: file_exists
    integer :: hdfier ! error flag for HDF5 library
    
    inquire( file=filename_h5, exist=file_exists ) ! check if file exists
    if (file_exists) then
      write(ounit,*) "output HDF5 file '",filename_h5,"' already exists"
    endif
    
    write(ounit,*) "writing input data into output HDF5 file '",filename_h5,"'"
    
    
    
    
    
    
    
  
  
  end subroutine mirror_input_to_outfile

end module sphdf5
