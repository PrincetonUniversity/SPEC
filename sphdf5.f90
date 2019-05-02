!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (output) ! Writes all the output information to ext.h5.

!latex \briefly{All the output information is contained in \type{ext.h5}.}
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

  INTEGER                        :: hdfier, rank  ! error flag for HDF5 library
  integer(hid_t)                 :: file_id, space_id, dset_id
  integer(hsize_t)               :: onedims(1:1), twodims(1:2), threedims(1:3)

contains

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


!
!subroutine HWRITEIV (dim1, varname, value )
!
!  LOCALS
!
!  integer, intent(in)      :: dim1
!  character(len=*), intent(in) :: varname
!  integer, intent(in)   :: value(:)
!
!#ifdef DEBUG
!  if( Wsphdf5 ) write(ounit,'("sphdf5 : ", 10x ," : myid=",i3," ; writing ",a20," ;")') myid, varname
!#endif
!
!  rank = 1 ; onedims(1) = dim1
!
!  if( dim1.le.0 ) then
!   write(ounit,'("sphdf5 : "10x" : error calling hwriteiv ; ",a20," : ",i3,".le.0 ;")') varname, dim1
!  else
!   call h5screate_simple_f( rank, onedims, space_id, hdfier )
!   FATAL( sphdf5, hdfier.ne.0, error calling h5screate_simple_f )
!
!   call h5dcreate_f( file_id, varname, H5T_NATIVE_INTEGER, space_id, dset_id, hdfier )
!   FATAL( sphdf5, hdfier.ne.0, error calling h5dcreate_f )
!
!   call h5dwrite_f( dset_id, H5T_NATIVE_INTEGER, value, onedims, hdfier )
!   FATAL( sphdf5, hdfier.ne.0, error calling h5dwrite_f )
!
!   call h5dclose_f(dset_id, hdfier)    ! terminate dataset;
!   FATAL( sphdf5, hdfier.ne.0, error calling h5dclose_f )
!  endif
!
!end subroutine HWRITEIV
!
!subroutine HWRITERV (dim1, varname, value )
!
!  LOCALS
!
!  integer, intent(in)      :: dim1
!  character(len=*), intent(in) :: varname
!  real, intent(in)      :: value(:)
!
!#ifdef DEBUG
!  if( Wsphdf5 ) write(ounit,'("sphdf5 : ", 10x ," : myid=",i3," ; writing ",a20," ;")') myid, varname
!#endif
!
!  rank = 1 ; onedims(1) = dim1
!
!  if( dim1.le.0 ) then
!   write(ounit,'("sphdf5 : "10x" : error calling hwriterv ; ",a20," : ",i3,".le.0 ;")') varname, dim1
!  else
!   call h5screate_simple_f( rank, onedims, space_id, hdfier )
!   FATAL( sphdf5, hdfier.ne.0, error calling h5screate_simple_f )
!
!   call h5dcreate_f( file_id, varname, H5T_NATIVE_DOUBLE, space_id, dset_id, hdfier )
!   FATAL( sphdf5, hdfier.ne.0, error calling h5dcreate_f )
!
!   call h5dwrite_f( dset_id, H5T_NATIVE_DOUBLE, value, onedims, hdfier )
!   FATAL( sphdf5, hdfier.ne.0, error calling h5dwrite_f )
!
!   call h5dclose_f(dset_id, hdfier)    ! terminate dataset;
!   FATAL( sphdf5, hdfier.ne.0, error calling h5dclose_f )
!  endif
!
!end subroutine HWRITERV
!
!subroutine HWRITERA (dim1, dim2, varname, value)
!
!  LOCALS
!
!  integer, intent(in)      :: dim1, dim2
!  character(len=*), intent(in) :: varname
!  real, intent(in)    :: value(:,:)
!
!#ifdef DEBUG
!  if( Wsphdf5 ) write(ounit,'("sphdf5 : ", 10x ," : myid=",i3," ; writing ",a20," ;")') myid, varname
!#endif
!
!  rank = 2 ; twodims(1:2) = (/ dim1, dim2 /)
!
!  if( dim1.le.0 .or. dim2.le.0 ) then
!   write(ounit,'("sphdf5 : "10x" : error calling hwritera ; ",a20," : ",i3,".le.0 .or. ",i3,".le.0 ;")') varname, dim1, dim2
!  else
!   call h5screate_simple_f( rank, twodims, space_id, hdfier )
!   FATAL( sphdf5, hdfier.ne.0, error calling h5screate_simple_f )
!
!   call h5dcreate_f( file_id, varname, H5T_NATIVE_DOUBLE, space_id, dset_id, hdfier )
!   FATAL( sphdf5, hdfier.ne.0, error calling h5dcreate_f )
!
!   call h5dwrite_f( dset_id, H5T_NATIVE_DOUBLE, value, twodims, hdfier )
!   FATAL( sphdf5, hdfier.ne.0, error calling h5dwrite_f )
!
!   call h5dclose_f(dset_id, hdfier)    ! terminate dataset;
!   FATAL( sphdf5, hdfier.ne.0, error calling h5dclose_f )
!  endif
!
!end subroutine HWRITERA
!
!subroutine HWRITERC (dim1, dim2, dim3, varname, value)
!
!  LOCALS
!
!  integer, intent(in)      :: dim1, dim2, dim3
!  character(len=*), intent(in) :: varname
!  real, intent(in)  :: value(:,:,:)
!
!#ifdef DEBUG
!  if( Wsphdf5 ) write(ounit,'("sphdf5 : ", 10x ," : myid=",i3," ; writing ",a20," ;")') myid, varname
!#endif
!
!  rank = 3 ; threedims(1:3) = (/ dim1, dim2, dim3 /)
!
!  if( dim1.le.0 .or. dim2.le.0 .or. dim3.le.0 ) then
!   write(ounit,'("sphdf5 : "10x" : error calling hwriterc ; ",a20," : ",i3,".le.0 .or. ",i3,".le.0 .or. ",i3,".le.0 ;")') varname, dim1, dim2, dim3
!  else
!   call h5screate_simple_f( rank, threedims, space_id, hdfier )
!   FATAL( sphdf5, hdfier.ne.0, error calling h5screate_simple_f )
!
!   call h5dcreate_f( file_id, varname, H5T_NATIVE_DOUBLE, space_id, dset_id, hdfier )
!   FATAL( sphdf5, hdfier.ne.0, error calling h5dcreate_f )
!
!   call h5dwrite_f( dset_id, H5T_NATIVE_DOUBLE, value, threedims, hdfier )
!   FATAL( sphdf5, hdfier.ne.0, error calling h5dwrite_f )
!
!   call h5dclose_f(dset_id, hdfier)    ! terminate dataset;
!   FATAL( sphdf5, hdfier.ne.0, error calling h5dclose_f )
!  endif
!
!end subroutine HWRITERC



subroutine init_outfile

  LOCALS
  integer(HID_T) :: plist_id      ! Property list identifier

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


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
subroutine finish_outfile

  LOCALS

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!





  call h5fclose_f( file_id, hdfier ) ! terminate access on output file;
  FATAL( sphdf5, hdfier.ne.0, error calling h5fclose_f )

  call h5close_f( hdfier ) ! close Fortran interface to the HDF5 library;
  FATAL( sphdf5, hdfier.ne.0, error calling h5close_f )
end subroutine finish_outfile

end module sphdf5
