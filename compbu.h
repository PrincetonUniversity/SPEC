!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (magnetic field) ! Calculates JB^u their derivatives

!latex \calledby{\link{dforce}}
!latex \calls{\link{}}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine compbu( lvol, lss, Ntz, mn, gBu, ideriv )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half, one, two, pi2
  
  use numerical, only : vsmall, small
  
  use fileunits, only : ounit
  
  use inputlist, only : Wcompbu, Igeometry, Ntor, Lrad, Mpol
  
  use cputiming, only : Tcompbu
  
  use allglobal, only : myid, cpus, pi2nfp, &
                        Mvol, im, in, halfmm, &
                        NOTstellsym, Lcoordinatesingularity, &
                        YESstellsym, zernike, cheby, &
                        Ate, Ato, Aze, Azo, &
                        efmn, ofmn, cfmn, sfmn, &
                        im, in, Nt, Nz, dBdX
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in) :: lvol, Ntz, mn, ideriv
  REAL, intent(in)    :: lss
  REAL, intent(out)   :: gBu(1:Ntz,1:3,0:3)

  REAL                :: Atemn(1:mn,0:2), Azemn(1:mn,0:2), Atomn(1:mn,0:2), Azomn(1:mn,0:2)
  REAL                :: sbar
  INTEGER             :: ii, jj, kk, ll, mm

  BEGIN(compbu)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  

  sbar = ( lss + one ) * half
  
  if (Lcoordinatesingularity) then ! Zernike 1 Jul 2019
    call get_zernike_d2(sbar, Lrad(lvol), mpol, zernike)
  else
    call get_cheby_d2(lss, Lrad(lvol), cheby(0:Lrad(lvol),0:2))
  endif

  Atemn(1:mn,0:2) = zero ! initialize summation over Chebyshev/Zernike polynomials;
  Azemn(1:mn,0:2) = zero
  if( NOTstellsym ) then
    Atomn(1:mn,0:2) = zero
    Azomn(1:mn,0:2) = zero
  else
    Atomn(1:mn,0:2) = zero ! these are used below;
    Azomn(1:mn,0:2) = zero
  endif
  
!latex \item The Fourier components of the vector potential given in \Eqn{At} and \Eqn{Az}, and their first and second radial derivatives, are summed.

  if (Lcoordinatesingularity) then
    do ll = 0, Lrad(lvol) ! radial (Chebyshev) resolution of magnetic vector potential;
    
      do ii = 1, mn  ! Fourier resolution of magnetic vector potential;
      mm = im(ii)
      if (ll < mm) cycle
      if (mod(ll+mm, 2) > 0) cycle
      
      ;Atemn(ii,0) = Atemn(ii,0) + Ate(lvol,ideriv,ii)%s(ll) * zernike(ll,mm,0)
      ;Atemn(ii,1) = Atemn(ii,1) + Ate(lvol,ideriv,ii)%s(ll) * zernike(ll,mm,1) * half
      ;Atemn(ii,2) = Atemn(ii,2) + Ate(lvol,ideriv,ii)%s(ll) * zernike(ll,mm,2) * half * half

      ;Azemn(ii,0) = Azemn(ii,0) + Aze(lvol,ideriv,ii)%s(ll) * zernike(ll,mm,0)
      ;Azemn(ii,1) = Azemn(ii,1) + Aze(lvol,ideriv,ii)%s(ll) * zernike(ll,mm,1) * half
      ;Azemn(ii,2) = Azemn(ii,2) + Aze(lvol,ideriv,ii)%s(ll) * zernike(ll,mm,2) * half * half

      if( NOTstellsym ) then

        Atomn(ii,0) = Atomn(ii,0) + Ato(lvol,ideriv,ii)%s(ll) * zernike(ll,mm,0)
        Atomn(ii,1) = Atomn(ii,1) + Ato(lvol,ideriv,ii)%s(ll) * zernike(ll,mm,1) * half
        Atomn(ii,2) = Atomn(ii,2) + Ato(lvol,ideriv,ii)%s(ll) * zernike(ll,mm,2) * half * half

        Azomn(ii,0) = Azomn(ii,0) + Azo(lvol,ideriv,ii)%s(ll) * zernike(ll,mm,0)
        Azomn(ii,1) = Azomn(ii,1) + Azo(lvol,ideriv,ii)%s(ll) * zernike(ll,mm,1) * half
        Azomn(ii,2) = Azomn(ii,2) + Azo(lvol,ideriv,ii)%s(ll) * zernike(ll,mm,2) * half * half

      endif ! end of if( NOTstellsym) ; 20 Jun 14;
        
      enddo ! end of do ii;
      
    enddo ! end of do ll;

  else

    do ll = 0, Lrad(lvol) ! radial (Chebyshev) resolution of magnetic vector potential;
    
      do ii = 1, mn  ! Fourier resolution of magnetic vector potential;
      
      ;Atemn(ii,0) = Atemn(ii,0) + Ate(lvol,ideriv,ii)%s(ll) * cheby(ll,0)
      ;Atemn(ii,1) = Atemn(ii,1) + Ate(lvol,ideriv,ii)%s(ll) * cheby(ll,1)
      ;Atemn(ii,2) = Atemn(ii,2) + Ate(lvol,ideriv,ii)%s(ll) * cheby(ll,2)

      ;Azemn(ii,0) = Azemn(ii,0) + Aze(lvol,ideriv,ii)%s(ll) * cheby(ll,0)
      ;Azemn(ii,1) = Azemn(ii,1) + Aze(lvol,ideriv,ii)%s(ll) * cheby(ll,1)
      ;Azemn(ii,2) = Azemn(ii,2) + Aze(lvol,ideriv,ii)%s(ll) * cheby(ll,2)

      if( NOTstellsym ) then

        Atomn(ii,0) = Atomn(ii,0) + Ato(lvol,ideriv,ii)%s(ll) * cheby(ll,0)
        Atomn(ii,1) = Atomn(ii,1) + Ato(lvol,ideriv,ii)%s(ll) * cheby(ll,1)
        Atomn(ii,2) = Atomn(ii,2) + Ato(lvol,ideriv,ii)%s(ll) * cheby(ll,2)

        Azomn(ii,0) = Azomn(ii,0) + Azo(lvol,ideriv,ii)%s(ll) * cheby(ll,0)
        Azomn(ii,1) = Azomn(ii,1) + Azo(lvol,ideriv,ii)%s(ll) * cheby(ll,1)
        Azomn(ii,2) = Azomn(ii,2) + Azo(lvol,ideriv,ii)%s(ll) * cheby(ll,2)

      endif ! end of if( NOTstellsym) ; 20 Jun 14;
        
      enddo ! end of do ii;
      
    enddo ! end of do ll;
  end if 
!latex \item The quantities $\sqrt g B^\s$, $\sqrt g B^\t$ and $\sqrt g B^\z$, and their first and second derivatives with respect to $(\s,\t,\z)$, 
!latex       are computed on the regular angular grid (using FFTs).
  ii = dbdx%ii
  ofmn(1:mn) = -          im(1:mn)*Azemn(1:mn,0) -          in(1:mn)*Atemn(1:mn,0)
  efmn(1:mn) = +          im(1:mn)*Azomn(1:mn,0) +          in(1:mn)*Atomn(1:mn,0)
  sfmn(1:mn) = -          im(1:mn)*Azemn(1:mn,1) -          in(1:mn)*Atemn(1:mn,1)
  cfmn(1:mn) = +          im(1:mn)*Azomn(1:mn,1) +          in(1:mn)*Atomn(1:mn,1)
  !if (ideriv.eq.-1) write(ounit,*) lvol, im(ii), in(ii), 'gB^s', ofmn(ii)
  call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, gBu(1:Ntz,1,0), gBu(1:Ntz,1,1) ) !  (gB^s)   , d(gB^s)/ds;
  
  efmn(1:mn) = - im(1:mn)*im(1:mn)*Azemn(1:mn,0) - im(1:mn)*in(1:mn)*Atemn(1:mn,0)
  ofmn(1:mn) = - im(1:mn)*im(1:mn)*Azomn(1:mn,0) - im(1:mn)*in(1:mn)*Atomn(1:mn,0)
  cfmn(1:mn) = + in(1:mn)*im(1:mn)*Azemn(1:mn,0) + in(1:mn)*in(1:mn)*Atemn(1:mn,0)
  sfmn(1:mn) = + in(1:mn)*im(1:mn)*Azomn(1:mn,0) + in(1:mn)*in(1:mn)*Atomn(1:mn,0)
  
  call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, gBu(1:Ntz,1,2), gBu(1:Ntz,1,3) ) ! d(gB^s)/dt, d(gB^s)/dz;
  
  efmn(1:mn) =                                   -                   Azemn(1:mn,1)
  ofmn(1:mn) =                                   -                   Azomn(1:mn,1)
  cfmn(1:mn) =                                   -                   Azemn(1:mn,2)
  sfmn(1:mn) =                                   -                   Azomn(1:mn,2)
  !if (ideriv.eq.-1) write(ounit,*) lvol, im(ii), in(ii), 'gB^t', efmn(ii)
  call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, gBu(1:Ntz,2,0), gBu(1:Ntz,2,1) ) !  (gB^t)   , d(gB^t)/ds;

  ofmn(1:mn) =                                   +          im(1:mn)*Azemn(1:mn,1)
  efmn(1:mn) =                                   -          im(1:mn)*Azomn(1:mn,1)
  sfmn(1:mn) =                                   -          in(1:mn)*Azemn(1:mn,1)
  cfmn(1:mn) =                                   +          in(1:mn)*Azomn(1:mn,1)
  
  call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, gBu(1:Ntz,2,2), gBu(1:Ntz,2,3) ) ! d(gB^t)/dt, d(gB^t)/dz;
  
  efmn(1:mn) = +                   Atemn(1:mn,1)
  ofmn(1:mn) = +                   Atomn(1:mn,1)
  cfmn(1:mn) = +                   Atemn(1:mn,2)
  sfmn(1:mn) = +                   Atomn(1:mn,2)
  !if (ideriv.eq.-1) write(ounit,*) lvol, im(ii), in(ii), 'gB^z', efmn(ii)
  call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, gBu(1:Ntz,3,0), gBu(1:Ntz,3,1) ) !  (gB^z)   , d(gB^z)/ds;

  ofmn(1:mn) = -          im(1:mn)*Atemn(1:mn,1)
  efmn(1:mn) = +          im(1:mn)*Atomn(1:mn,1)
  sfmn(1:mn) = +          in(1:mn)*Atemn(1:mn,1)
  cfmn(1:mn) = -          in(1:mn)*Atomn(1:mn,1)
  
  call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, gBu(1:Ntz,3,2), gBu(1:Ntz,3,3) ) ! d(gB^z)/dt, d(gB^z)/dz;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(compbu)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine compbu

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
