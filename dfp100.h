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
								      IsMyVolume, IsMyVolumeValue

 LOCALS
!------

INTEGER              :: vvol, NN, ll, ii, Ndofgl, iflag

REAL                 :: Fvec(1:Mvol-1), x(1:Mvol-1)



BEGIN(dfp100)

	! First set up the value of dpflux for calculation

	dpflux(2:Mvol) = x - xoffset

	do vvol = 1, Mvol

		LREGION(vvol) ! assigns Lcoordinatesingularity, Lplasmaregion, etc. ;

		ImagneticOK(vvol) = .false.

	!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

		!if( myid.ne.modulo(vvol-1,ncpu) ) goto 5000 ! construct Beltrami fields in parallel;

		
		if(( myid .EQ. 0 ) .and. (.not.LocalConstraint)) then
	 	  IsMyVolumeValue = 1 ! for now, execution only by master CPU in case of global constraints
                else
		  call IsMyVolume(vvol)
		endif

		if( IsMyVolumeValue .EQ. 0 ) then
			goto 5000
		else if( IsMyVolumeValue .EQ. -1) then
			FATAL(dfp100, .true., Unassociated volume)
		endif
			

		NN = NAdof(vvol) ! shorthand;

		ll = Lrad(vvol)

	!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

		dBdX%L = .false. ! first, compute Beltrami fields;

	!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

		WCALL( dfp100, ma02aa, ( vvol, NN ) )

	!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
		
    WCALL( dfp100, lbpol, (vvol) )

		5000 continue
	enddo


	if( .not.LocalConstraint ) then

		select case (Lconstraint)

			case( 3 )
				! Compute IPDt on each interface. Eventually need to put the do
				! loop inside the surfcurent subroutine... TODO
				do ii = 1, Mvol-1
          IPDt(ii) = -pi2 * (Btemn(1, 0, ii+1) - Btemn(1, 1, ii)) / mu0
					!WCALL(dfp100, surfcurent, (ii, mn))
				enddo

				Fvec = IPDt - Isurf(1:Mvol-1)
				!Ddof = ??? TODO: SEE IF AN ANALYTICAL FORMULATION EXISTS...

			case default
				FATAL(dfp100, .true., Unaccepted value for Lconstraint)

		end select
	endif

RETURN(dfp100)

end subroutine dfp100