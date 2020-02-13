!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (lbpol) ! Computes Btheta at the interface

!latex \briefly{Computes Btheta at the interface - used to compute the toroidal surface current}

!latex \calledby{\link{xspech} and
!latex   		 \link{dfp100}}
!latex \calls{\link{coords} and 
!latex 		  \link{numrec}}


!latex \begin{enumerate}
!latex \item Call \link{coords} to compute the metric coefficients and the jacobian.
!latex \item Build coefficients \inputvar{efmn}, \inputvar{ofmn}, \inputvar{cfmn}, \inputvar{sfmn} from the field vector potential \inputvar{Ate}, \inputvar{Ato}, 
!latex 		 \inputvar{Aze} and \inputvar{Azo}, and radial derivatives of the Chebyshev polynomials \inputvar{TT(ll,innout,1)}. These variables
!latex		 are the radial derivative of the Fourier coefficients of the magnetic field vector potential. 
!latex \item Take the inverse Fourier transform of \inputvar{efmn}, \inputvar{ofmn}, \inputvar{cfmn}, \inputvar{sfmn}. These are the covariant components of $dA$, 
!latex 		 \textit{i.e.} the contravariant components of $\mathbf{B}$.
!latex \item Build covariant components of the field using the metric coefficients \inputvar{guvij} and the jacobian \inputvar{sg}.
!latex \item Fourier transform the covariant components of the field and store them in the variables \inputvar{Btemn}, \inputvar{Btomn}, \inputvar{Bzemn} and 
!latex  	 \inputvar{Bzomn}.
!latex \end{enumerate}


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine lbpol(lvol, ideriv)

  use constants, only : mu0, pi, pi2, two, one, half, zero

  use allglobal, only : Ate, Aze, Ato, Azo, TT, &
                        YESstellsym, NOTstellsym, &
                        im, in, mne, ime, ine, Mvol, mn, &
                        sg, guvij, &
                        Ntz, Lcoordinatesingularity, &
                        efmn, ofmn, cfmn, sfmn, evmn, odmn, comn, simn, &
                        Btemn, Bzemn, Btomn, Bzomn, &
                        Nt, Nz, &
                        regumm, &
                        cpus, myid

  use inputlist, only : Lrad, Wlbpol, Igeometry

  use fileunits, only : ounit

  use cputiming, only : Tlbpol

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS
! ------
  
  INTEGER                :: Lcurvature, ideriv, ii, ll, ifail, lvol, mi, ni
  REAL                   :: lss, innout
  REAL                   :: lAte(1:mn), lAze(1:mn), lAto(1:mn), lAzo(1:mn)
  REAL                   :: dAt(1:Ntz), dAz(1:Ntz), Bt(1:Ntz), Bz(1:Ntz)
  REAL                   :: dBtzero	      ! Value of first B_theta mode jump
  REAL                   :: mfactor           ! Regularisation factor

! Lcurvature:             Controls what the routine coords computes.
! lint:                   Interface number
! vvol, ideriv, ii, ll:   Some iteration variables
! lAte, lAze, lAto, lAzo: Reconstructed vector potential

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  BEGIN(lbpol)

! First get the metric component and jacobian

  Lcurvature = 1

! innout=0 -> inner boundary of volume (s=-1) and innout=1 -> outer boundary (s=1)

  Btemn(1:mn, 0:1, lvol) = zero
  do innout=0,1

  lss = two * innout - one
  
  if((lvol==1) .and. (Igeometry/=1)) then ; Lcoordinatesingularity = .true.;
  else; Lcoordinatesingularity = .false.;
  endif

	if( Lcoordinatesingularity .and. innout.EQ.0) then
	  goto 5555; ! No need to compute at the singularity (and crash with debug...)
	endif

  WCALL( lbpol, coords, (lvol, lss, Lcurvature, Ntz, mn ) ) ! get guvij and sg
  

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! Then compute the vector potential and its derivatives. Copied from sc00aa.h
 
   efmn(1:mn) = zero ; ofmn(1:mn) = zero ; cfmn(1:mn) = zero ; sfmn(1:mn) = zero
   do ii = 1, mn ; mi = im(ii) ; ni = in(ii) ! loop over Fourier harmonics;
    
! In case of singularity, point at sbar=0 not computed - no problem here!
! For definition of the regularisation factor, see jo00aa documentation.
   if( Lcoordinatesingularity ) then ; mfactor = regumm(ii) * half ! derivative of regularisation factor;
   else                              ; mfactor = zero
   endif
   
   do ll = 0, Lrad(lvol) ! loop over Chebyshev polynomials; Lrad is the radial resolution;
   ! Note that the minus sine is included at line 122-123
      ;                      ; efmn(ii) = efmn(ii) + Ate(lvol,ideriv,ii)%s(ll) * ( TT(ll,innout,1) + mfactor ) ! B^\t;
      ;                      ; cfmn(ii) = cfmn(ii) + Aze(lvol,ideriv,ii)%s(ll) * ( TT(ll,innout,1) + mfactor ) ! B^\z;
      if( NOTstellsym ) then ; ofmn(ii) = ofmn(ii) + Ato(lvol,ideriv,ii)%s(ll) * ( TT(ll,innout,1) + mfactor )
        ;                    ; sfmn(ii) = sfmn(ii) + Azo(lvol,ideriv,ii)%s(ll) * ( TT(ll,innout,1) + mfactor )
      endif
    enddo ! end of do ll; 20 Feb 13;
    
  enddo ! end of do ii; 20 Feb 13;

! Inverse Fourier transform to map to real space
  call invfft( mn, im, in, efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, dAt(1:Ntz), dAz(1:Ntz) ) ! get covariant component of dA / contravariant of B

  Bt(1:Ntz) = ( - dAz(1:Ntz) * guvij(1:Ntz,2,2,0) + dAt(1:Ntz) * guvij(1:Ntz,2,3,0) ) / sg(1:Ntz,0) ! Get covariant components
  Bz(1:Ntz) = ( - dAz(1:Ntz) * guvij(1:Ntz,2,3,0) + dAt(1:Ntz) * guvij(1:Ntz,3,3,0) ) / sg(1:Ntz,0)

! Fourier transform, map to Fourier space
  ifail = 0
  call tfft( Nt, Nz, Bt(1:Ntz), Bz(1:Ntz), &
              mn, im(1:mn), in(1:mn), Btemn(1:mn,innout,lvol), Btomn(1:mn,innout,lvol), Bzemn(1:mn,innout,lvol), Bzomn(1:mn,innout,lvol), ifail )

5555 continue
  enddo ! end of do innout;


! Now Btemn(1, 0, vvol) and Btemn(1, 1, vvol) contain Bet00(s=-1) and Bet00(s=1) for each volume vvol.


  RETURN(lbpol)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine lbpol






