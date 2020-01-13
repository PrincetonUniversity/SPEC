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
  
  logical :: verbose = .FALSE. ! set to true to enable verbal diarreha during reading of input file
  
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

    integer(hid_t)              :: file_id
    integer(hid_t)              :: dset_id, dtype_id, dtype_id_native
    integer(hid_t)              :: dataspace         ! dataspace used to query dataset size
    integer(hsize_t)            :: dims_1(1)         ! current dimensions of rank-1 dataset
    integer(hsize_t)            :: max_dims_1(1)     ! maximum dimensions of rank-1 dataset
    logical                     :: item_exists, datatypes_equal
    integer                     :: item_type
    
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
        
        ! The following scheme to verify the existence, correct type and reasonable value of the dataset
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
        !  8. read the dataset using the verified native datatype
        !  9. close the native datatype
        ! 10. close the dataset object
        
        
        
        ! query existence of /version link
        call h5lexists_f(file_id, "/version", item_exists, hdfier)
        if (hdfier.ne.0) then
          write(*,*) "error checking if /version link exists"
        else
          if (verbose) then ; write(ounit,*) "successfully checked if the /version link exists" ; endif
        
          ! check that /version link exists
          if (item_exists) then
            if (verbose) then; write(ounit,*) "successfully checked that the link /version exists" ; endif
           
            ! query existence of object at /version link
            call h5oexists_by_name_f(file_id, "/version", item_exists, hdfier)
            if (hdfier.ne.0) then
              write(*,*) "error checking for presence of object at /version"
            else
              if (verbose) then ; write(ounit,*) "successfully checked for presence of an object at /version link" ; endif
            
              ! check that there exists an item at the /version link
              if (item_exists) then
                if (verbose) then ; write(ounit,*) "successfully checked that there exists an item at /version" ; endif
              
                ! try to open /version object
                call h5oopen_f(file_id, "/version", dset_id, hdfier)
                if (hdfier.ne.0) then
                  write(*,*) "error opening object /version"
                else
                  if (verbose) then ; write(ounit,*) "successfully opened object at /version" ; endif
                
                  ! query type of /version object
                  call h5iget_type_f(dset_id, item_type, hdfier)
                  if (hdfier.ne.0) then
                    write(*,*) "error querying type of object /version"
                  else
                    if (verbose) then ; write(ounit,*) "successfully queried item type of /version" ; endif
                  
                    ! check if /version is a dataset
                    if (item_type.ne.H5I_DATASET_F) then
                      write(ounit,*) "error verifying that object /version is a Dataset(",H5I_DATASET_F,"); rather it is ",item_type
                    else
                      if (verbose) then ; write(ounit,*) "successfully verified that /version is a Dataset" ; endif
                      
                      ! query datatype of /version
                      call h5dget_type_f(dset_id, dtype_id, hdfier)
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
                              
                              ! read version dataset
                              call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, input_version, int((/1/), HSIZE_T), hdfier)
                              if (hdfier.ne.0) then
                                write(*,*) "error reading dataset /version"
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
                              endif ! read version dataset
                            else
                              write(ounit,*) "native datatype of /version should be H5T_NATIVE_DOUBLE but is ",dtype_id_native
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
                    endif ! check if /version is a dataset
                    
                  endif ! query type of /version object
                          
                  call h5oclose_f(dset_id, hdfier)
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
        
        
        
        
        
        
        
        
        
        
      
        ! /input/Igeometry --> Igeometry
        call h5dopen_f(file_id, "/input/Igeometry", dset_id, hdfier)
        if (hdfier.ne.0) then
          write(*,*) "error opening dataset '/input/Igeometry'"
        else
          call h5dread_f(dset_id, H5T_NATIVE_INTEGER, Igeometry, int((/1/), HSIZE_T), hdfier)
          if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/Igeometry'" ; endif
          call h5dclose_f(dset_id, hdfier)
          if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/Igeometry'" ; endif
        endif
        
        ! /input/Nvol --> Nvol
        call h5dopen_f(file_id, "/input/Nvol", dset_id, hdfier)
        if (hdfier.ne.0) then
          write(*,*) "error opening dataset '/input/Nvol'"
        else
          call h5dread_f(dset_id, H5T_NATIVE_INTEGER, Nvol, int((/1/), HSIZE_T), hdfier)
          if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/Nvol'" ; endif
          call h5dclose_f(dset_id, hdfier)
          if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/Nvol'" ; endif
        endif
        
        ! /input/Lrad --> Lrad(1:dims_1(1)); rank=1
        call h5dopen_f(file_id, "/input/Lrad", dset_id, hdfier)
        if (hdfier.ne.0) then
          write(*,*) "error opening dataset '/input/Lrad'"
        else
          ! open dataspace to get current state of dataset
          call h5dget_space_f(dset_id, dataspace, hdfier)
          if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '/input/Lrad'" ; endif
          
          ! get current size of dataset
          call h5sget_simple_extent_dims_f(dataspace, dims_1, max_dims_1, hdfier)
          if (hdfier.ne.1) then ; write(*,*) "unexpected rank of dataset '/input/Lrad': ",hdfier," .ne. 1" ; endif
          
          ! close dataspace after it has been used to query the size of the variable
          call h5sclose_f(dataspace, hdfier)
          if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '/input/Lrad'" ; endif
          
          call h5dread_f(dset_id, H5T_NATIVE_INTEGER, Lrad(1:dims_1(1)), dims_1, hdfier)
          if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/Lrad'" ; endif
          
          call h5dclose_f(dset_id, hdfier)
          if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/Lrad'" ; endif
        endif
        
        ! /input/phiedge --> phiedge
        call h5dopen_f(file_id, "/input/phiedge", dset_id, hdfier)
        if (hdfier.ne.0) then
          write(*,*) "error opening dataset '/input/phiedge'"
        else
          call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, phiedge, int((/1/), HSIZE_T), hdfier)
          if (hdfier.ne.0) then ; write(*,*) "error reading dataset '/input/phiedge'" ; endif
          call h5dclose_f(dset_id, hdfier)
          if (hdfier.ne.0) then ; write(*,*) "error closing dataset '/input/phiedge'" ; endif
        endif
        
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
  
  implicit none
  double precision :: version = 2.1

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
