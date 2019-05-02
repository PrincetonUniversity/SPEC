!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (output) ! Writes all the output information to ext.h5.

!latex \briefly{All the output information is contained in \type{ext.h5}.}
!latex \calledby{\link{xspech}}
!l tex \calls{\link{}}

!latex \tableofcontents

!latex \newcommand{\pb}[1]{\parbox[t]{13cm}{#1}}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

module sphdf5

contains

subroutine init_outfile
#ifndef NOHDF5
  use hdf5
  use allglobal , only : myid

  LOCALS
  INTEGER :: hdfier

  call h5open_f( hdfier ) ! initialize Fortran interface;
  FATAL( hdfint, hdfier.ne.0, error calling h5open_f )


#endif
end subroutine init_outfile



subroutine finish_outfile
#ifndef NOHDF5
  use hdf5
  use allglobal , only : myid

  LOCALS
  INTEGER :: hdfier

  call h5close_f( hdfier ) ! close Fortran interface;
  FATAL( hdfint, hdfier.ne.0, error calling h5close_f )
#endif
end subroutine finish_outfile

end module sphdf5
