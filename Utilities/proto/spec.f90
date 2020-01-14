program spec
  use hdf5
  use fileunits
  !use inputlist
  use allglobal
  use sphdf5
  
  implicit none
  
  integer :: numargs    ! number of command-line arguments
  integer :: extlen     !< length of input filename
  integer :: sppos      !< position of ".sp"    in input filename
  integer :: sph5pos    !< position of ".sp.h5" in input filename
  integer :: infile_src !< 0: no input file found; 1: textual namelist in ext.sp; 2: HDF5 file with prepared input contents
  integer :: hdfier     ! error for HDF5 library
  
  write(*,'("SPEC version = ",f5.2)') version
  
  ! open HDF5 library
  call h5open_f(hdfier)
  if (hdfier.ne.0) then
    write(*,*) "error opening interface to HDF5 library"
    call exit(1)
  endif
  
  
  infile_src = 0 ! only ext on cmd line is the default
  
  ! check that at least one command-line argument is given
  numargs = command_argument_count()
  if (numargs.ne.1) then
    write(ounit,'(" number of command-line argument should be exactly 1 (name of the input file), but it is ",i2)') numargs
    call exit(1)
  endif
  
  ! parse first command-line argument; should be the name of the input file
  call getarg( 1, ext )
  extlen = len_trim(ext)
  
  sph5pos = index(ext, ".sp.h5", .true.)                             ! search for ".sp.h5" from the back of ext
  if (sph5pos.gt.0 .and. sph5pos.eq.extlen-5) then                   ! check if ext ends with ".sp.h5"
    ext = ext(1:extlen-6)                                            ! if this is the case, remove ".sp.h5" from end of ext
    infile_src = 2
  endif
  
  sppos = index(ext, ".sp", .true.)                                  ! search for ".sp" from the back of ext
  if (infile_src.eq.0 .and. sppos.gt.0 .and. sppos.eq.extlen-2) then ! check if ext ends with ".sp"; only if not a HDF5 file was specified instead
    ext = ext(1:extlen-3)                                            ! if this is the case, remove ".sp" from end of ext
    infile_src = 1
  endif
  
  write(ounit,*) "ext = '",trim(ext),"'"
  if (infile_src.eq.0) then
    write(ounit,*) "input file name specified by ext only  --> textual namelist"
  elseif (infile_src.eq.1) then
    write(ounit,*) "input file name specified by ext.sp    --> textual namelist"
  elseif (infile_src.eq.2) then
    write(ounit,*) "input file name specified by ext.sp.h5 --> HDF5 file for input"
  endif
  
  if (infile_src.eq.0 .or. infile_src.eq.1) then
    ! populate the data in inputlist with the namelist contents
    call readin(trim(ext)//".sp")
  elseif (infile_src.eq.2) then
    ! read input data from HDF5 file
    call readin_h5(trim(ext)//".sp.h5")
  endif
  
  ! debugging output to check that reading the input file was successful
  write(ounit,'(" Igeometry = ",i1)')      Igeometry
  write(ounit,'("      Nvol = ",i1)')      Nvol
  write(ounit,'("      Lrad = ",10(i1x))') Lrad(1:Nvol)
  write(ounit,'("   phiedge = ",f4.2)')    phiedge

  ! In case the input data was given as a namelist, mirror it into the HDF5 output file.
  ! For the other case of a HDF5 input file, the input data is already there and SPEC will simply
  ! put its output data into the same file.
  if (infile_src.eq.0 .or. infile_src.eq.1) then
    ! mirror the input file contents to the output file
    call mirror_input_to_outfile(trim(ext)//".sp.h5")
  endif 
  
  ! close HDF5 library
  call h5close_f(hdfier)
  if (hdfier.ne.0) then
    write(*,*) "error closing interface to HDF5 library"
    call exit(1)
  endif
  
  call exit(0)
end program spec
