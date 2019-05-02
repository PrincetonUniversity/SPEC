!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (output) ! Writes all the output information to ext.h5.

!latex \briefly{All the output information is contained in \type{ext.h5}.}
!latex \calledby{\link{xspech}}
!l tex \calls{\link{}}

!latex \tableofcontents

!latex \newcommand{\pb}[1]{\parbox[t]{13cm}{#1}}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

module sphdf5
  implicit none

contains

subroutine init_outfile
#ifdef NOHDF5
#else
  use hdf5


#endif
end subroutine init_outfile

end module sphdf5
