!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (surfcurent) ! Computes pressure driven currents on each interface

!latex \briefly{briefly}

!latex \calledby{?????}
!latex \calls{?????}

!latex \tableofcontents

!latex Computes the pressure driven current integral at the interface labelled by \inputvar{lint},
!latex
!latex \be
!latex I^\mathcal{S}_\phi = \int_0^{2\pi} [[B_\theta]] d\theta.
!latex \ee
!latex
!latex The magnetic field is computed as in \link{sc00aa}. This is used only when Lconstraint=3, \textit{i.e.} when the toroidal current profile is constraint.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine surfcurent(lint, mn)

  use constants, only : mu0, pi, pi2, two, one, half, zero

  use allglobal, only : Ate, Aze, Ato, Azo, TT, &
                        YESstellsym, NOTstellsym, &
                        im, in, mne, ime, ine, Mvol, &
                        sg, guvij, &
                        Ntz, Lcoordinatesingularity, &
                        efmn, ofmn, cfmn, sfmn, evmn, odmn, comn, simn, &
                        Btemn, Bzemn, Btomn, Bzomn, &
                        Nt, Nz, &
                        regumm, &
                        IPDt, &
                        cpus, myid

  use inputlist, only : Lrad, Wsurfcurent, Igeometry

  use fileunits, only : ounit

  use cputiming, only : Tsurfcurent

! Some useful notes:
! TT:                  Derivatives of Chebyshev polynomials at the inner and outer interfaces
! Ate, Aze, Ato, Azo:  Fourier coefficients of the vector potential. 't' and 'z' stands for
!                      poloidal and toroidal component, and 'e', 'o' stands for even / odd.
!                      Ate = Ate(1:Mvol,-1:2,1:mn)%s(0:Lrad(vvol), 0)
! im, in:              Fourier modes, set in readin;               
! mne, ime, ine:       Enhanced resolution for metric element.
! Mvol:                Equal to the number of volumes (+1 if free boundary)
! Ntz:                 Discrete resolution, Ntz = Nt*Nz
! sg(0:3,Ntz) contains the jacobian and its derivatives
! guvij(0:6,0:3,1:Ntz) contains the metric elements and their derivatives
! Lrad:                Radial resolution (max order of Chebyshev polynomials)
! ounit:               Used to write some stuff in terminal

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS
! ------
  
  INTEGER, intent(in)    :: mn                ! Total number of Fourier harmonics;
  INTEGER                :: Lcurvature, lint, vvol, ideriv, ii, ll, ifail, lvol, mi, ni
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

  BEGIN(surfcurent)

! First get the metric component and jacobian

  Lcurvature = 1

! innout=0 -> inner boundary of volume (s=-1) and innout=1 -> outer boundary (s=1)

  do innout=0,1

  lss = two * innout - one

  if (innout==0) then; lvol = lint+1
  else;                  lvol = lint
  endif
  
  ! Quick dirty fix to test this routine. ATTENTION NEED TO BE REMOVED 
  if((lvol==1) .and. (Igeometry/=1)) then ; Lcoordinatesingularity = .true.;
  else; Lcoordinatesingularity = .false.;
  endif

  WCALL( surfcurent, coords, (lvol, lss, Lcurvature, Ntz, mn ) ) ! get guvij and sg
  

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! Then compute the vector potential and its derivatives. Copied from sc00aa.h

  ideriv = 0; 
   efmn(1:mn) = zero ; ofmn(1:mn) = zero ; cfmn(1:mn) = zero ; sfmn(1:mn) = zero
   do ii = 1, mn ; mi = im(ii) ; ni = in(ii) ! loop over Fourier harmonics;
    

   if( Lcoordinatesingularity ) then ; mfactor = regumm(ii) * half ! regularity factor;
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

#ifdef DEBUG
  write(ounit, '("surfcurent : ", f10.2, ", lint=",i3, ", lss=", f10.1, ", lvol=",i3,", Bt(1)=",ES23.15, ", Lcoordsing=", i3)') cput-cpus, lint, lss, lvol, Bt(1), Lcoordinatesingularity
#endif

  enddo ! end of do innout;

! Get the jump in B_theta in Fourier space for the first even mode.
  dBtzero = Btemn(1, 0, lint+1) - Btemn(1, 1, lint)


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! Then compute the line integral. Only the first even mode integral is non-zero (periodic integral)

  IPDt(lint) = -pi2 * dBtzero / mu0


  RETURN(surfcurent)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine surfcurent






