!title (&ldquo;global&rdquo; dfp100) ! Computes global constraint and take care of MPI communications

!latex \briefly{Split the work between MPI nodes and evaluate the global constraint}

!latex \calledby{\link{dforce}} \\

!latex \calls{\link{ma02aa} and
!latex        \link{lbpol}}


!latex \tableofcontents

!latex \subsection{Local constraint}
!latex In the case of a local constraint, each MPI node calls \link{ma02aa} to obtain the field consistent with the local constraint

!latex \subsection{Global constraint}
!latex In the case of a global constraint, 

!latex \begin{enumerate}
!latex \item The MPI node $0$ broadcast the value of \internal{iflag}. \internal{iflag} is equal to $5$ if the constraint has been
!latex 		 matched, and all nodes exit the iteration process. Otherwise, an additional iteration occurs (point 2 and below).
!latex \item The MPI node $0$ boradcast the new values of $\psi_p$.
!latex \item Each node calls \link{ma02aa} to obtain the field for the given value of $\mu$ and $\psi_p$
!latex \item Each node calls \link{lbpol} to evalue the poloidal magnetic field Fourier coefficients $B_{\theta,e,0,0}$ at the inner and outer interface
!latex 		 (this is only for a toroidal current constraint. Other global constraints would require different parameters)
!latex \item Each node communicate to the MPI node $0$ the parameters necessary for the evaluation of the constraint
!latex \item The node $0$ evaluate the global constraint and sets the flag \internal{IconstraintOK}
!latex \end{enumerate}

!latex \subsection{Evaluation of constraint derivatives}
!latex Not implemented for now.



subroutine dfp100(Ndofgl, x, Fvec, iflag)

use constants, only : zero, half, one, two, pi2, pi, mu0

use fileunits, only : ounit

use inputlist, only : Wmacros, Wdfp100, Igeometry, Nvol, Lrad, Isurf, &
					  					Lconstraint, Lfindzero

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
! vvol: 											loop index on volumes
! Ndofgl:											Input parameter necessary for the use of hybrd1. Unused otherwise.
! iflag:											Flag changed by hybrd1
! cpu_send_one, cpu_send_two: CPU IDs, used for MPI communications
! tag1, tag2: 								Communication tags
! status:											MPI status
! Fvec:												Global constraint values
! x:													Degrees of freedom of hybrd1. For now contains only the poloidal flux

INTEGER              :: vvol, Ndofgl, iflag, cpu_send_one, cpu_send_two, tag1, tag2
INTEGER							 :: status(MPI_STATUS_SIZE), request1, request2
REAL                 :: Fvec(1:Mvol-1), x(1:Mvol-1)
LOGICAL				 :: LcomputeDerivatives



BEGIN(dfp100)

	if( LocalConstraint ) then
		dpflux(2:Mvol) = x - xoffset
		IconstraintOK = .false.
	else
		! In case of global constraints, dfp100 is called via hybrd1, a modified version of hybrd. hybrd1
		! calls dfp100 one additional time than hybrd1, at the end of the iteration process, with iflag=5.
		! This is the signal that the iteration is over, and the master thread can broadcast the information
		! to its slaves that they can exit the infinit loop.
		
		if( (iflag.EQ.5) .and. (myid.EQ.0) ) then
			IconstraintOK = .true.

			if( Lfindzero.eq.2 ) then! Compute derivatives of solution w.r.t mu and pflux
				LcomputeDerivatives = .true.
				WCALL( dfp100, ma02aa, ( vvol, NAdof(vvol), LcomputeDerivatives ) )
			endif
		else
			IconstraintOK = .false.
		endif

		! Master broadcast to slaves the info 
		call MPI_BCAST( IconstraintOK, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

		! If constraints are OK, no need for one additional computation, skip till the end.
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
	endif


	! Now each CPU perform the calculation in its volume(s)
	do vvol = 1, Mvol

		LREGION(vvol) ! assigns Lcoordinatesingularity, Lplasmaregion, etc. ;
		ImagneticOK(vvol) = .false.

		!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

		! Determines if this volume vvol should be computed by this thread.
		call IsMyVolume(vvol)

		if( IsMyVolumeValue .EQ. 0 ) then
			cycle
		else if( IsMyVolumeValue .EQ. -1) then
			FATAL(dfp100, .true., Unassociated volume)
		endif

		!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

		dBdX%L = .false. ! first, compute Beltrami fields;

		!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

		! Compute fields
		LcomputeDerivatives = .false.
		WCALL( dfp100, ma02aa, ( vvol, NAdof(vvol), LcomputeDerivatives ) )

		! Compute relevant local quantities for the evaluation of the constraint. Doing it like this 
		! reduces the amount of data sent to the master thread. In the case of current constraint, only two
		! doubles per volume are sent.
		if( Lconstraint.EQ.3 ) then
			WCALL( dfp100, lbpol, (vvol, 0) )				!Compute field at interface for global constraint
		endif
	enddo


	! Evaluation of global constraint and communications
	if( .not.LocalConstraint ) then

		select case (Lconstraint)

			! Case 3: toroidal current constraint
			case( 3 )

				! Compute IPDt on each interface. Eventually need to put the do
				! loop inside the surfcurent subroutine... TODO
				do vvol = 1, Mvol-1
				
					! --------------------------------------------------------------------------------------------------
					! 																	MPI COMMUNICATIONS
					call WhichCpuID(vvol  , cpu_send_one)
					call WhichCpuID(vvol+1, cpu_send_two)
					tag1 = 10*vvol + 1
					tag2 = 10*vvol + 2

					! Send poloidal magnetic field at the interface. Non blocking - what happens next depends only on
					! the master thread.
					if( myid.EQ.cpu_send_one ) then
						call MPI_ISEND(Btemn(1, 1, vvol  ), 1, MPI_DOUBLE_PRECISION,            0, tag1, MPI_COMM_WORLD, request1, ierr)
					endif
					if( myid.EQ.cpu_send_two ) then
						call MPI_ISEND(Btemn(1, 0, vvol+1), 1, MPI_DOUBLE_PRECISION,            0, tag2, MPI_COMM_WORLD, request2, ierr)
					endif

					! Master thread receives the poloidal magnetic field at each interface
					if( myid.EQ.0 ) then
						call MPI_RECV(Btemn(1, 0, vvol+1), 1, MPI_DOUBLE_PRECISION, cpu_send_two, tag2, MPI_COMM_WORLD, status, ierr)
						call MPI_RECV(Btemn(1, 1, vvol  ), 1, MPI_DOUBLE_PRECISION, cpu_send_one, tag1, MPI_COMM_WORLD, status, ierr)
          				IPDt(vvol) = pi2 * (Btemn(1, 0, vvol+1) - Btemn(1, 1, vvol))
					endif

					if( myid.EQ.cpu_send_one ) then
						call MPI_WAIT(request1, status, ierr)
					endif

					if( myid.EQ.cpu_send_two ) then
						call MPI_WAIT(request2, status, ierr)
					endif

				enddo

				! Compute the constraint and store it in Fvec. TODO: Compute analytically the constraint jacobian -
				! this would improve significaly the performances...
				if( myid.EQ.0 ) then
					Fvec = IPDt - Isurf(1:Mvol-1)
					!Ddof = ??? TODO: SEE IF AN ANALYTICAL FORMULATION EXISTS...
				endif

#ifdef DEBUG 
        		!write(ounit, '("dfp100: ", 10x ," : max(IPDt) = "es12.5)') MAXVAL(IPDt)
#endif
				! write(ounit,'("xspech : ", 10x ," : sum(Ate(",i3,",",i2,",",i2,")%s) =",99es23.15)') vvol, ideriv, ii, sum(Ate(vvol,ideriv,ii)%s(0:Lrad(vvol)))
			
			case default
				FATAL(dfp100, .true., Unaccepted value for Lconstraint)
		end select
	endif

6666 continue

RETURN(dfp100)

end subroutine dfp100