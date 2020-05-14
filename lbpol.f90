!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (lbpol) ! Computes Btheta at the interface

!latex \briefly{Computes Btheta at the interface - used to compute the toroidal surface current}

!latex \calledby{\link{xspech} and
!latex            \link{dfp100}}
!latex \calls{\link{coords} and 
!latex           \link{numrec}}


!latex \begin{enumerate}
!latex \item Call \link{coords} to compute the metric coefficients and the jacobian.
!latex \item Build coefficients \inputvar{efmn}, \inputvar{ofmn}, \inputvar{cfmn}, \inputvar{sfmn} from the field vector potential \inputvar{Ate}, \inputvar{Ato}, 
!latex          \inputvar{Aze} and \inputvar{Azo}, and radial derivatives of the Chebyshev polynomials \inputvar{TT(ll,innout,1)}. These variables
!latex         are the radial derivative of the Fourier coefficients of the magnetic field vector potential. 
!latex \item Take the inverse Fourier transform of \inputvar{efmn}, \inputvar{ofmn}, \inputvar{cfmn}, \inputvar{sfmn}. These are the covariant components of $dA$, 
!latex          \textit{i.e.} the contravariant components of $\mathbf{B}$.
!latex \item Build covariant components of the field using the metric coefficients \inputvar{guvij} and the jacobian \inputvar{sg}.
!latex \item Fourier transform the covariant components of the field and store them in the variables \inputvar{Btemn}, \inputvar{Btomn}, \inputvar{Bzemn} and 
!latex       \inputvar{Bzomn}.
!latex \end{enumerate}


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine lbpol(lvol, Bt00, ideriv)

  use constants, only : mu0, pi, pi2, two, one, half, zero

  use allglobal, only : Ate, Aze, Ato, Azo, TT, &
                        YESstellsym, NOTstellsym, &
                        im, in, mne, ime, ine, Mvol, mn, &
                        sg, guvij, &
                        Ntz, Lcoordinatesingularity, &
                        efmn, ofmn, cfmn, sfmn, evmn, odmn, comn, simn, &
                        Nt, Nz, &
                        regumm, &
                        cpus, myid, dBdX, &
                        build_vector_potential

  use inputlist, only : Lrad, Wlbpol, Igeometry, Lcheck

  use fileunits, only : ounit

  use cputiming, only : Tlbpol

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS
! ------
  
  INTEGER                :: Lcurvature, ideriv, ii, ll, ifail, lvol, mi, ni, iocons
  REAL                   :: lss, Bt00(1:Mvol, 0:1)
  REAL                   :: lAte(1:mn), lAze(1:mn), lAto(1:mn), lAzo(1:mn)
  REAL                   :: dAt(1:Ntz), dAz(1:Ntz), Bt(1:Ntz), Bz(1:Ntz), dAt0(1:Ntz), dAz0(1:Ntz)
  REAL                   :: dBtzero      ! Value of first B_theta mode jump
  REAL                   :: mfactor           ! Regularisation factor
  LOGICAL                :: LGeometricDerivative

! Lcurvature:             Controls what the routine coords computes.
! lint:                   Interface number
! vvol, ideriv, ii, ll:   Some iteration variables
! lAte, lAze, lAto, lAzo: Reconstructed vector potential

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  BEGIN(lbpol)

! TODO: This subroutine is very similar to curent.f90 - maybe merge both in a single subroutine to simplify?

! iocons=0 -> inner boundary of volume (s=-1) and iocons=1 -> outer boundary (s=1)

  do iocons=0,1

    lss = two * iocons - one
  
    if( Lcoordinatesingularity .and. iocons.EQ.0) then
      goto 5555; ! No need to compute at the singularity
    endif

    if( lvol.eq.Mvol .and. iocons.eq.1) then
      goto 5555;
    endif

    ! First get the metric component and jacobian
    Lcurvature = 1
    WCALL( lbpol, coords, (lvol, lss, Lcurvature, Ntz, mn ) ) ! get guvij and sg
    
    ! Then compute the vector potential and its derivatives.
    call build_vector_potential(lvol, iocons, ideriv, 1)

    ! Inverse Fourier transform to map to real space
    call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, dAt(1:Ntz), dAz(1:Ntz) ) ! get covariant component of dA / contravariant of B

    ! Construct covariant Fourier components of B
    Bt(1:Ntz) = ( - dAz(1:Ntz ) * guvij(1:Ntz,2,2, 0) + dAt(1:Ntz ) * guvij(1:Ntz,2,3, 0) )/ sg(1:Ntz,0)
    Bz(1:Ntz) = ( - dAz(1:Ntz ) * guvij(1:Ntz,2,3, 0) + dAt(1:Ntz ) * guvij(1:Ntz,3,3, 0) )/ sg(1:Ntz,0)


    if( ideriv.eq.-1 ) then

      ! Get derivatives of metric element
      Lcurvature = 5
      WCALL( lbpol, coords, (lvol, lss, Lcurvature, Ntz, mn ) ) ! get guvij over sg derivatives

      ! Compute vector potential without taking derivatives
      call build_vector_potential(lvol, iocons, 0, 1)

      ! And now add variation of metric contribution
      call invfft( mn, im, in, efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, dAt0(1:Ntz), dAz0(1:Ntz) ) ! get covariant component of dA without derivatives


      Bt(1:Ntz) = Bt(1:Ntz) + ( - dAz0(1:Ntz ) * guvij(1:Ntz,2,2, 1) + dAt0(1:Ntz ) * guvij(1:Ntz,2,3, 1) ) ! Add metric derivatives
      Bz(1:Ntz) = Bz(1:Ntz) + ( - dAz0(1:Ntz ) * guvij(1:Ntz,2,3, 1) + dAt0(1:Ntz ) * guvij(1:Ntz,3,3, 1) ) 

    endif

    ! Fourier transform, map to Fourier space
    ifail = 0
    call tfft( Nt, Nz, Bt(1:Ntz), Bz(1:Ntz), mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), ifail )

    Bt00(lvol, iocons) = efmn(1)

5555 continue
  enddo ! end of do iocons;


! Now Btemn(1, 0, vvol) and Btemn(1, 1, vvol) contain Bte00(s=-1) and Bte00(s=1) for each volume vvol.


  RETURN(lbpol)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine lbpol






