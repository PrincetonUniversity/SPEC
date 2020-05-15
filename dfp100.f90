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
!latex          matched, and all nodes exit the iteration process. Otherwise, an additional iteration occurs (point 2 and below).
!latex \item The MPI node $0$ boradcast the new values of $\psi_p$.
!latex \item Each node calls \link{ma02aa} to obtain the field for the given value of $\mu$ and $\psi_p$
!latex \item Each node calls \link{lbpol} to evalue the poloidal magnetic field Fourier coefficients $B_{\theta,e,0,0}$ at the inner and outer interface
!latex          (this is only for a toroidal current constraint. Other global constraints would require different parameters)
!latex \item Each node communicate to the MPI node $0$ the parameters necessary for the evaluation of the constraint
!latex \item The node $0$ evaluate the global constraint and sets the flag \internal{IconstraintOK}
!latex \end{enumerate}

!latex \subsection{Evaluation of constraint derivatives}
!latex Not implemented for now.



subroutine dfp100(Ndofgl, x, Fvec, LComputeDerivatives)

use constants, only : zero, half, one, two, pi2, pi, mu0

use fileunits, only : ounit

use inputlist, only : Wmacros, Wdfp100, Igeometry, Nvol, Lrad, Isurf, &
                      Lconstraint

use cputiming, only : Tdfp100

use allglobal, only : ncpu, myid, cpus, &
                      ImagneticOK, NAdof, mn, &
                      Mvol, Iquad, &
                      dBdX, &
                      Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, Localconstraint, &
                      IPDt, IPDtdPf, xoffset, dpflux, &
                      IsMyVolume, IsMyVolumeValue, WhichCpuID, &
                      IconstraintOK, &
                      DToocc, DToocs, DToosc, DTooss, &
                      TTsscc, TTsscs, TTsssc, TTssss, &
                      TDstcc, TDstcs, TDstsc, TDstss, &
                      TDszcc, TDszcs, TDszsc, TDszss, &
                      DDttcc, DDttcs, DDttsc, DDttss, &
                      DDtzcc, DDtzcs, DDtzsc, DDtzss, &
                      DDzzcc, DDzzcs, DDzzsc, DDzzss, &
                      dMA, dMB, dMD, dMG, MBpsi, solution

 LOCALS
!------
! vvol:                       loop index on volumes
! Ndofgl:                     Input parameter necessary for the use of hybrd1. Unused otherwise.
! iflag:                      Flag changed by hybrd1
! cpu_send_one, cpu_send_two: CPU IDs, used for MPI communications
! status:                     MPI status
! Fvec:                       Global constraint values
! x:                          Degrees of freedom of hybrd1. For now contains only the poloidal flux

INTEGER              :: vvol, Ndofgl, iflag, cpu_send_one, cpu_send_two, ll, NN
INTEGER              :: status(MPI_STATUS_SIZE), request1, request2
REAL                 :: Fvec(1:Mvol-1), x(1:Mvol-1), Bt00(1:Mvol, 0:1, 0:1)
LOGICAL              :: LComputeDerivatives



BEGIN(dfp100)

    dpflux(2:Mvol) = x - xoffset

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

        
        ll = Lrad(vvol)
        NN = NAdof(vvol)

        SALLOCATE( DToocc, (0:ll,0:ll,1:mn,1:mn), zero )
        SALLOCATE( DToocs, (0:ll,0:ll,1:mn,1:mn), zero )
        SALLOCATE( DToosc, (0:ll,0:ll,1:mn,1:mn), zero )
        SALLOCATE( DTooss, (0:ll,0:ll,1:mn,1:mn), zero )

        SALLOCATE( TTsscc, (0:ll,0:ll,1:mn,1:mn), zero )
        SALLOCATE( TTsscs, (0:ll,0:ll,1:mn,1:mn), zero )
        SALLOCATE( TTsssc, (0:ll,0:ll,1:mn,1:mn), zero )
        SALLOCATE( TTssss, (0:ll,0:ll,1:mn,1:mn), zero )

        SALLOCATE( TDstcc, (0:ll,0:ll,1:mn,1:mn), zero )
        SALLOCATE( TDstcs, (0:ll,0:ll,1:mn,1:mn), zero )
        SALLOCATE( TDstsc, (0:ll,0:ll,1:mn,1:mn), zero )
        SALLOCATE( TDstss, (0:ll,0:ll,1:mn,1:mn), zero )

        SALLOCATE( TDszcc, (0:ll,0:ll,1:mn,1:mn), zero )
        SALLOCATE( TDszcs, (0:ll,0:ll,1:mn,1:mn), zero )
        SALLOCATE( TDszsc, (0:ll,0:ll,1:mn,1:mn), zero )
        SALLOCATE( TDszss, (0:ll,0:ll,1:mn,1:mn), zero )

        SALLOCATE( DDttcc, (0:ll,0:ll,1:mn,1:mn), zero )
        SALLOCATE( DDttcs, (0:ll,0:ll,1:mn,1:mn), zero )
        SALLOCATE( DDttsc, (0:ll,0:ll,1:mn,1:mn), zero )
        SALLOCATE( DDttss, (0:ll,0:ll,1:mn,1:mn), zero )

        SALLOCATE( DDtzcc, (0:ll,0:ll,1:mn,1:mn), zero )
        SALLOCATE( DDtzcs, (0:ll,0:ll,1:mn,1:mn), zero )
        SALLOCATE( DDtzsc, (0:ll,0:ll,1:mn,1:mn), zero )
        SALLOCATE( DDtzss, (0:ll,0:ll,1:mn,1:mn), zero )

        SALLOCATE( DDzzcc, (0:ll,0:ll,1:mn,1:mn), zero )
        SALLOCATE( DDzzcs, (0:ll,0:ll,1:mn,1:mn), zero )
        SALLOCATE( DDzzsc, (0:ll,0:ll,1:mn,1:mn), zero )
        SALLOCATE( DDzzss, (0:ll,0:ll,1:mn,1:mn), zero )

        SALLOCATE( dMA, (0:NN, 0:NN), zero )
        SALLOCATE( dMB, (0:NN, 0: 2), zero )
        SALLOCATE( dMD, (0:NN, 0:NN), zero )
        SALLOCATE( dMG, (0:NN      ), zero )

        SALLOCATE( solution, (1:NN, -1:2), zero )
        SALLOCATE( MBpsi, (1:NN), zero )

        dBdX%L = .false. ! No need to take derivatives of matrices w.r.t geometry.

        ! Compute matrices
        WCALL( dfp100, ma00aa, ( Iquad(vvol), mn, vvol, ll ) )
        WCALL( dfp100, matrix, ( vvol, mn, ll ) )

        ! Compute fields
        WCALL( dfp100, ma02aa, ( vvol, NAdof(vvol) ) )

        !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

        ! Free memory

        DALLOCATE(MBpsi)
        DALLOCATE(solution)
        DALLOCATE(dMG)
        DALLOCATE(dMD)
        DALLOCATE(dMB)
        DALLOCATE(dMA)

        DALLOCATE(DToocc)
        DALLOCATE(DToocs)
        DALLOCATE(DToosc)
        DALLOCATE(DTooss)

        DALLOCATE(TTsscc)
        DALLOCATE(TTsscs)
        DALLOCATE(TTsssc)
        DALLOCATE(TTssss)

        DALLOCATE(TDstcc)
        DALLOCATE(TDstcs)
        DALLOCATE(TDstsc)
        DALLOCATE(TDstss)

        DALLOCATE(TDszcc)
        DALLOCATE(TDszcs)
        DALLOCATE(TDszsc)
        DALLOCATE(TDszss)

        DALLOCATE(DDttcc)
        DALLOCATE(DDttcs)
        DALLOCATE(DDttsc)
        DALLOCATE(DDttss)

        DALLOCATE(DDtzcc)
        DALLOCATE(DDtzcs)
        DALLOCATE(DDtzsc)
        DALLOCATE(DDtzss)

        DALLOCATE(DDzzcc)
        DALLOCATE(DDzzcs)
        DALLOCATE(DDzzsc)
        DALLOCATE(DDzzss)



        ! Compute relevant local quantities for the evaluation of global constraints
        if( Lconstraint.EQ.3 ) then
            WCALL( dfp100, lbpol, (vvol, Bt00(1:Mvol, 0:1, 0), 0) )                !Compute field at interface for global constraint
            WCALL( dfp100, lbpol, (vvol, Bt00(1:Mvol, 0:1, 1), 2) )                !Compute field at interface for global constraint, d w.r.t. pflux
        endif
    enddo


    ! Evaluation of global constraint and communications
    if( .not.LocalConstraint ) then

        select case (Lconstraint)

            ! Case 3: toroidal current constraint
            case( 3 )

                ! Compute IPDt on each interface.
                do vvol = 1, Mvol-1
                
                    ! --------------------------------------------------------------------------------------------------
                    !                                                                     MPI COMMUNICATIONS
                    call WhichCpuID(vvol  , cpu_send_one)
                    call WhichCpuID(vvol+1, cpu_send_two)

                    ! Broadcast magnetic field at the interface.
                    RlBCAST(Bt00(vvol  , 1, 0), 1, cpu_send_one)
                    RlBCAST(Bt00(vvol+1, 0, 0), 1, cpu_send_two)
                    RlBCAST(Bt00(vvol  , 1, 1), 1, cpu_send_one)
                    RlBCAST(Bt00(vvol+1, 0, 1), 1, cpu_send_two)

                    ! Evaluate surface current
                    IPDt(vvol) = pi2 * (Bt00(vvol+1, 0,0 ) - Bt00(vvol, 1, 0))
                    
                    ! their derivatives
                    IPDtdPf(vvol,vvol) = pi2 * Bt00(vvol+1, 0, 1)
                    if (vvol .ne. 1) IPDtdPf(vvol,vvol-1) = -pi2 * Bt00(vvol, 1, 1)
                enddo

                ! Compute the constraint and store it in Fvec. TODO: Compute analytically the constraint jacobian ?
                if( myid.EQ.0 ) then
                    Fvec = IPDt - Isurf(1:Mvol-1)
                endif

#ifdef DEBUG 
                write(ounit, '("dfp100: ", 10x ," : max(IPDt) = "es12.5)') MAXVAL(IPDt)
#endif
                ! write(ounit,'("xspech : ", 10x ," : sum(Ate(",i3,",",i2,",",i2,")%s) =",99es23.15)') vvol, ideriv, ii, sum(Ate(vvol,ideriv,ii)%s(0:Lrad(vvol)))
            
            case default
                FATAL(dfp100, .true., Unaccepted value for Lconstraint)
        end select
    endif

6666 continue

RETURN(dfp100)

end subroutine dfp100
