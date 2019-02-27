!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (surfcurent) ! Computes pressure driven currents on each interface

!latex \briefly{briefly}

!latex \calledby{?????}
!latex \calls{?????}

!latex \tableofcontents

!latex TO COMPLETE 

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

  use inputlist, only : Lrad, Wsurfcurent

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

  if (innout.EQ.0) then; lvol = lint+1
  else;                  lvol = lint
  endif
    
  WCALL( surfcurent, coords, (lvol, lss, Lcurvature, Ntz, mn ) ) ! get guvij and sg
  

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! Then compute the vector potential and its derivatives. Copied from sc00aa.h

  ideriv = 0; 

   do ii = 1, mn ; mi = im(ii) ; ni = in(ii) ! loop over Fourier harmonics;
    
   if( Lcoordinatesingularity ) then ; mfactor = regumm(ii) * half ! regularity factor;
   else                              ; mfactor = zero
   endif
   
   do ll = 0, Lrad(lvol) ! loop over Chebyshev polynomials; Lrad is the radial resolution;
   ! Note that the minus sine is included at line 110
      ;                      ; efmn(ii) = efmn(ii) + Ate(lvol,ideriv,ii)%s(ll) * ( TT(ll,innout,1) + mfactor ) ! B^\t;
      ;                      ; cfmn(ii) = cfmn(ii) + Aze(lvol,ideriv,ii)%s(ll) * ( TT(ll,innout,1) + mfactor ) ! B^\z;
      if( NOTstellsym ) then ; ofmn(ii) = ofmn(ii) + Ato(lvol,ideriv,ii)%s(ll) * ( TT(ll,innout,1) + mfactor )
        ;                    ; sfmn(ii) = sfmn(ii) + Azo(lvol,ideriv,ii)%s(ll) * ( TT(ll,innout,1) + mfactor )
      endif
    enddo ! end of do ll; 20 Feb 13;
    
  enddo ! end of do ii; 20 Feb 13;

  call invfft( mn, im, in, efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, dAt(1:Ntz), dAz(1:Ntz) ) ! get covariant component of dA / contravariant of B

  Bt(1:Ntz) = ( - dAz(1:Ntz) * guvij(1:Ntz,2,2,0) + dAt(1:Ntz) * guvij(1:Ntz,2,3,0) ) / sg(1:Ntz,0) ! Get covariant components

  ifail = 0
  call tfft( Nt, Nz, Bt(1:Ntz), Bz(1:Ntz), &
              mn, im(1:mn), in(1:mn), Btemn(1:mn,innout,lvol), Btomn(1:mn,innout,lvol), Bzemn(1:mn,innout,lvol), Bzomn(1:mn,innout,lvol), ifail )

!#ifdef DEBUG
!  write(ounit, '("surfcurent : ", f10.2)') cput-cpus
!#endif

  enddo ! end of do innout;

  dBtzero = Btemn(1, 0, lint+1) - Btemn(1, 1, lint)


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! Then compute the line integral

  IPDt(lint) = pi2 * dBtzero / mu0


  RETURN(surfcurent)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine surfcurent



! Attention: Lcoordinatesingularity has to have the correct value once called...





