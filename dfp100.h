subroutine dfp100(Ndofgl, x, Fvec, iflag)

use constants, only : zero, half, one, two, pi2, pi, mu0

use fileunits, only : ounit

use inputlist, only : Wmacros, Wdfp100, Igeometry, Nvol, Lrad, Isurf, &
					  Lconstraint

use cputiming, only : Tdfp100

use allglobal, only : ncpu, myid, cpus, &
                      ImagneticOK, NAdof, mn, &
                      Mvol, &
                      dBdX, &
                      Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, Localconstraint, &
                      IPDt, xoffset, dpflux, &
											Btemn, &
								      IsMyVolume, IsMyVolumeValue, WhichCpuID, &
											IconstraintOK

 LOCALS
!------

INTEGER              :: vvol, NN, ll, ii, Ndofgl, iflag, cpu_send_one, cpu_send_two, tag1, tag2
INTEGER							 :: status(MPI_STATUS_SIZE), request1, request2
REAL                 :: Fvec(1:Mvol-1), x(1:Mvol-1)



BEGIN(dfp100)

	
	if( (iflag.EQ.5) .and. (myid.EQ.0) ) then
		IconstraintOK = .true.
	else
		IconstraintOK = .false.
	endif

	call MPI_BCAST( IconstraintOK, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

	if( IconstraintOK ) then
		!write(*,*) "Exiting dfp100.h with iflag=5"
		goto 6666
	endif


	! First set up the value of dpflux for calculation

  if( myid.EQ.0 ) then
		dpflux(2:Mvol) = x - xoffset
	endif

	! We could scatter dpflux on the volumes - but this depends on how the volumes are divided between the CPUs
	! We would need to restructure the data - might be complicated for little gain.
	call MPI_Bcast( dpflux, Mvol, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

	! Now each CPU perform the calculation in its volume(s)
	do vvol = 1, Mvol

		LREGION(vvol) ! assigns Lcoordinatesingularity, Lplasmaregion, etc. ;

		ImagneticOK(vvol) = .false.

	!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

	if( .not.LocalConstraint ) then
		 
		call IsMyVolume(vvol)

		if( IsMyVolumeValue .EQ. 0 ) then
			goto 5000
		else if( IsMyVolumeValue .EQ. -1) then
			FATAL(dfp100, .true., Unassociated volume)
		endif
	endif

	!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

		NN = NAdof(vvol) ! shorthand;
		ll = Lrad(vvol)

		dBdX%L = .false. ! first, compute Beltrami fields;

	!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

		WCALL( dfp100, ma02aa, ( vvol, NN ) ) !Compute fields

	!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
		
		if( Lconstraint.EQ.3 ) then
		  WCALL( dfp100, lbpol, (vvol) )				!Compute field at interface for global constraint
		endif

		5000 continue
	enddo


	if( .not.LocalConstraint ) then

		select case (Lconstraint)

			case( 3 )
				! Compute IPDt on each interface. Eventually need to put the do
				! loop inside the surfcurent subroutine... TODO
				do ii = 1, Mvol-1
				
					call WhichCpuID(ii  , cpu_send_one)
					call WhichCpuID(ii+1, cpu_send_two)
					tag1 = 10*ii + 1
					tag2 = 10*ii + 2
										
					if( myid.EQ.cpu_send_one ) then
						call MPI_ISEND(Btemn(1, 0, ii  ), 1, MPI_DOUBLE_PRECISION,            0, tag1, MPI_COMM_WORLD, request1, ierr)
					endif
		
					if( myid.EQ.cpu_send_two ) then
						call MPI_ISEND(Btemn(1, 0, ii+1), 1, MPI_DOUBLE_PRECISION,            0, tag2, MPI_COMM_WORLD, request2, ierr)
					endif

					if( myid.EQ.0 ) then
						call MPI_RECV(Btemn(1, 0, ii+1), 1, MPI_DOUBLE_PRECISION, cpu_send_two, tag2, MPI_COMM_WORLD, status, ierr)
						call MPI_RECV(Btemn(1, 0, ii  ), 1, MPI_DOUBLE_PRECISION, cpu_send_one, tag1, MPI_COMM_WORLD, status, ierr)
          	IPDt(ii) = -pi2 * (Btemn(1, 0, ii+1) - Btemn(1, 1, ii)) / mu0
					endif

					if( myid.EQ.cpu_send_one ) then
						call MPI_WAIT(request1, status, ierr)
					endif

					if( myid.EQ.cpu_send_two ) then
						call MPI_WAIT(request2, status, ierr)
					endif

				enddo

				if( myid.EQ.0 ) then
					Fvec = IPDt - Isurf(1:Mvol-1)
					!Ddof = ??? TODO: SEE IF AN ANALYTICAL FORMULATION EXISTS...
				endif

			case default
				FATAL(dfp100, .true., Unaccepted value for Lconstraint)

		end select

	endif

6666 continue

RETURN(dfp100)

end subroutine dfp100
