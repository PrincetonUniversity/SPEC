subroutine dfp100(Ndofgl, x, Fvec, iflag)

use constants, only : zero, half, one, two, pi2, pi

use fileunits, only : ounit

use inputlist, only : Wmacros, Wdfp100, Igeometry, Nvol, Lrad, Isurf, &
					  Lconstraint

use cputiming, only : Tdfp100

use allglobal, only : Ate, Aze, Ato, Azo, &
                      NOTstellsym, &
		      ncpu, myid, cpus, &
                      ImagneticOK, NAdof, mn, &
                      Mvol, &
                      dBdX, Iquad, &
                      Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, Localconstraint, &
                      IPDt, xoffset, dpflux


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

		if( myid.ne.modulo(vvol-1,ncpu) ) goto 5000 ! construct Beltrami fields in parallel;

		NN = NAdof(vvol) ! shorthand;

		ll = Lrad(vvol)

	!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

		dBdX%L = .false. ! first, compute Beltrami fields;

	!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

		WCALL( dfp100, ma02aa, ( vvol, NN ) )

	!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
		

		! Broadcast relevant information to everyone
		!do ii = 1, mn  
   		!	RlBCAST( Ate(vvol, 0, ii)%s(0:Lrad(vvol)), Lrad(vvol)+1, myid )
  	!		RlBCAST( Aze(vvol, 0, ii)%s(0:Lrad(vvol)), Lrad(vvol)+1, myid )
		!enddo
  	!	
	!	if( NOTstellsym ) then
  	!		do ii = 1, mn    
	!			RlBCAST( Ato(vvol, 0, ii)%s(0:Lrad(vvol)), Lrad(vvol)+1, myid )
   	!			RlBCAST( Azo(vvol, 0, ii)%s(0:Lrad(vvol)), Lrad(vvol)+1, myid )
	!		enddo
	!	endif

		5000 continue
	enddo

	! call MPI_Barrier(MPI_COMM_WORLD, ierr)

	if( .not.LocalConstraint ) then

		select case (Lconstraint)

			case( 3 )
				! Compute IPDt on each interface. Eventually need to put the do
				! loop inside the surfcurent subroutine... TODO
				do ii = 1, Mvol-1
					WCALL(dfp100, surfcurent, (ii, mn))
				enddo

				Fvec = IPDt - Isurf(1:Mvol-1)
				!Ddof = ??? TODO: SEE IF AN ANALYTICAL FORMULATION EXISTS...

			case default
				FATAL(dfp100, .true., Unaccepted value for Lconstraint)

		end select
	endif

RETURN(dfp100)

end subroutine dfp100
