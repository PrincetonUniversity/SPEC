!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (geometry perturbation) ! Calculates \xi_s, \xi_\theta, \xi_\zeta and their derivatives on the surface of the interface

!latex \calledby{\link{dforce}}
!latex \calls{\link{}}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine compxi( lvol, Ntz, mn, dxi )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half, one, two, pi2
  
  use numerical, only : vsmall, small
  
  use fileunits, only : ounit
  
  use inputlist, only : Wcompxi, Igeometry, Ntor
  
  use cputiming, only : Tcompxi
  
  use allglobal, only : myid, cpus, pi2nfp, &
                        Mvol, im, in, halfmm, &
                        NOTstellsym, Lcoordinatesingularity, &
                        Nt, Nz, &
                        Rij, Zij, &
                        cosi, sini, &
                        sg, &
                        dBdX, &
                        YESstellsym
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in) :: lvol, Ntz, mn
  REAL, intent(out)   :: dxi(1:Ntz,1:3,0:3)


  INTEGER             :: ii, innout, irz, issym
  REAL                :: lss

  BEGIN(compxi)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  dxi(1:Ntz,1:3,0:3) = zero

  ! one must have called coords using Lcurvature=2 BEFORE calling compxi

  ii = dBdX%ii ; innout = dBdX%innout ; irz = dBdX%irz ; issym = dBdX%issym ! shorthand

  ! compute J * xi

  select case (Igeometry)

  case (1) ! slab

    if (issym .eq. 0) then
      dxi(1:Ntz,1,0) = cosi(1:Ntz,ii)
      dxi(1:Ntz,1,1) = zero                       ! not used
      dxi(1:Ntz,1,2) = - im(ii) * sini(1:Ntz,ii)
      dxi(1:Ntz,1,3) = + in(ii) * sini(1:Ntz,ii)
    else ! NOTstellsym
      dxi(1:Ntz,1,0) = sini(1:Ntz,ii)
      dxi(1:Ntz,1,1) = zero                       ! not used
      dxi(1:Ntz,1,2) = + im(ii) * cosi(1:Ntz,ii)
      dxi(1:Ntz,1,3) = - in(ii) * cosi(1:Ntz,ii)
    end if ! issym .eq. 0

  case (2) ! cylinder

    if (issym .eq. 0) then
      dxi(1:Ntz,1,0) = cosi(1:Ntz,ii) * Rij(1:Ntz,0,0)
      dxi(1:Ntz,1,1) = zero                       ! not used
      dxi(1:Ntz,1,2) = (- im(ii) * sini(1:Ntz,ii)) * Rij(1:Ntz,0,0) + cosi(1:Ntz,ii) * Rij(1:Ntz,2,0)
      dxi(1:Ntz,1,3) = (+ in(ii) * sini(1:Ntz,ii)) * Rij(1:Ntz,0,0) + cosi(1:Ntz,ii) * Rij(1:Ntz,3,0)
    else ! NOTstellsym
      dxi(1:Ntz,1,0) = sini(1:Ntz, ii) * Rij(1:Ntz,0,0)
      dxi(1:Ntz,1,1) = zero                       ! not used
      dxi(1:Ntz,1,2) = (+ im(ii) * cosi(1:Ntz,ii)) * Rij(1:Ntz,0,0) + sini(1:Ntz,ii) * Rij(1:Ntz,2,0)
      dxi(1:Ntz,1,3) = (- in(ii) * cosi(1:Ntz,ii)) * Rij(1:Ntz,0,0) + sini(1:Ntz,ii) * Rij(1:Ntz,3,0)
    end if ! issym .eq. 0

  case (3) ! toroidal

    if (issym .eq. 0) then
      if (irz .eq. 0) then ! \delta R
        dxi(1:Ntz,1,0) =- cosi(1:Ntz,ii) * Zij(1:Ntz,2,0) * Rij(1:Ntz,0,0)
        dxi(1:Ntz,1,1) = zero                      ! not used
        dxi(1:Ntz,1,2) =+ im(ii) * sini(1:Ntz,ii) * Zij(1:Ntz,2,0) * Rij(1:Ntz,0,0) &
                        - cosi(1:Ntz,ii) * Zij(1:Ntz,2,2) * Rij(1:Ntz,0,0) &
                        - cosi(1:Ntz,ii) * Zij(1:Ntz,2,0) * Rij(1:Ntz,2,0)
        dxi(1:Ntz,1,3) =- in(ii) * sini(1:Ntz,ii) * Zij(1:Ntz,2,0) * Rij(1:Ntz,0,0) & 
                        - cosi(1:Ntz,ii) * Zij(1:Ntz,2,3) * Rij(1:Ntz,0,0) &
                        - cosi(1:Ntz,ii) * Zij(1:Ntz,2,0) * Rij(1:Ntz,3,0)

        dxi(1:Ntz,2,0) =+ cosi(1:Ntz,ii) * Zij(1:Ntz,1,0) * Rij(1:Ntz,0,0)
        dxi(1:Ntz,2,1) = zero                      ! not used
        dxi(1:Ntz,2,2) =- im(ii) * sini(1:Ntz,ii) * Zij(1:Ntz,1,0) * Rij(1:Ntz,0,0) &
                        + cosi(1:Ntz,ii) * Zij(1:Ntz,1,2) * Rij(1:Ntz,0,0) &
                        + cosi(1:Ntz,ii) * Zij(1:Ntz,1,0) * Rij(1:Ntz,2,0)
        dxi(1:Ntz,2,3) =+ in(ii) * sini(1:Ntz,ii) * Zij(1:Ntz,1,0) * Rij(1:Ntz,0,0) & 
                        + cosi(1:Ntz,ii) * Zij(1:Ntz,1,3) * Rij(1:Ntz,0,0) &
                        + cosi(1:Ntz,ii) * Zij(1:Ntz,1,0) * Rij(1:Ntz,3,0)
      else ! \delta Z
        dxi(1:Ntz,1,0) =+ sini(1:Ntz,ii) * Rij(1:Ntz,2,0) * Rij(1:Ntz,0,0)
        dxi(1:Ntz,1,1) = zero                      ! not used
        dxi(1:Ntz,1,2) =+ im(ii) * cosi(1:Ntz,ii) * Rij(1:Ntz,2,0) * Rij(1:Ntz,0,0) &
                        + sini(1:Ntz,ii) * Rij(1:Ntz,2,2) * Rij(1:Ntz,0,0) &
                        + sini(1:Ntz,ii) * Rij(1:Ntz,2,0) * Rij(1:Ntz,2,0)
        dxi(1:Ntz,1,3) =- in(ii) * cosi(1:Ntz,ii) * Rij(1:Ntz,2,0) * Rij(1:Ntz,0,0) & 
                        + sini(1:Ntz,ii) * Rij(1:Ntz,2,3) * Rij(1:Ntz,0,0) &
                        + sini(1:Ntz,ii) * Rij(1:Ntz,2,0) * Rij(1:Ntz,3,0)

        dxi(1:Ntz,2,0) =- sini(1:Ntz,ii) * Rij(1:Ntz,1,0) * Rij(1:Ntz,0,0)
        dxi(1:Ntz,2,1) = zero                      ! not used
        dxi(1:Ntz,2,2) =- im(ii) * cosi(1:Ntz,ii) * Rij(1:Ntz,1,0) * Rij(1:Ntz,0,0) &
                        - sini(1:Ntz,ii) * Rij(1:Ntz,1,2) * Rij(1:Ntz,0,0) &
                        - sini(1:Ntz,ii) * Rij(1:Ntz,1,0) * Rij(1:Ntz,2,0)
        dxi(1:Ntz,2,3) =+ in(ii) * cosi(1:Ntz,ii) * Rij(1:Ntz,1,0) * Rij(1:Ntz,0,0) & 
                        - sini(1:Ntz,ii) * Rij(1:Ntz,1,3) * Rij(1:Ntz,0,0) &
                        - sini(1:Ntz,ii) * Rij(1:Ntz,1,0) * Rij(1:Ntz,3,0)
      end if ! irz .eq. 0
    else ! NOTstellsym
      if (irz .eq. 0) then ! \delta R
        dxi(1:Ntz,1,0) =- sini(1:Ntz,ii) * Zij(1:Ntz,2,0) * Rij(1:Ntz,0,0)
        dxi(1:Ntz,1,1) = zero                      ! not used
        dxi(1:Ntz,1,2) =- im(ii) * cosi(1:Ntz,ii) * Zij(1:Ntz,2,0) * Rij(1:Ntz,0,0) &
                        - sini(1:Ntz,ii) * Zij(1:Ntz,2,2) * Rij(1:Ntz,0,0) &
                        - sini(1:Ntz,ii) * Zij(1:Ntz,2,0) * Rij(1:Ntz,2,0)
        dxi(1:Ntz,1,3) =+ in(ii) * cosi(1:Ntz,ii) * Zij(1:Ntz,2,0) * Rij(1:Ntz,0,0) & 
                        - sini(1:Ntz,ii) * Zij(1:Ntz,2,3) * Rij(1:Ntz,0,0) &
                        - sini(1:Ntz,ii) * Zij(1:Ntz,2,0) * Rij(1:Ntz,3,0)
      else ! \delta Z
        dxi(1:Ntz,1,0) =+ cosi(1:Ntz,ii) * Rij(1:Ntz,2,0) * Rij(1:Ntz,0,0)
        dxi(1:Ntz,1,1) = zero                      ! not used
        dxi(1:Ntz,1,2) =- im(ii) * sini(1:Ntz,ii) * Rij(1:Ntz,2,0) * Rij(1:Ntz,0,0) &
                        + cosi(1:Ntz,ii) * Rij(1:Ntz,2,2) * Rij(1:Ntz,0,0) &
                        + cosi(1:Ntz,ii) * Rij(1:Ntz,2,0) * Rij(1:Ntz,2,0)
        dxi(1:Ntz,1,3) =+ in(ii) * sini(1:Ntz,ii) * Rij(1:Ntz,2,0) * Rij(1:Ntz,0,0) & 
                        + cosi(1:Ntz,ii) * Rij(1:Ntz,2,3) * Rij(1:Ntz,0,0) &
                        + cosi(1:Ntz,ii) * Rij(1:Ntz,2,0) * Rij(1:Ntz,3,0)
      end if ! irz .eq. 0
    end if ! issym .eq. 0

  end select ! Igeometry

  ! divide by the Jacobian
  dxi(1:Ntz,1,0) = dxi(1:Ntz,1,0) / sg(1:Ntz,0)
  dxi(1:Ntz,1,1) = zero                            ! not used
  dxi(1:Ntz,1,2) = dxi(1:Ntz,1,2) / sg(1:Ntz,0) - dxi(1:Ntz,1,0) * sg(1:Ntz,2) / sg(1:Ntz,0)
  dxi(1:Ntz,1,3) = dxi(1:Ntz,1,3) / sg(1:Ntz,0) - dxi(1:Ntz,1,0) * sg(1:Ntz,3) / sg(1:Ntz,0)

  if (Igeometry .eq. 3) then
    dxi(1:Ntz,2,0) = dxi(1:Ntz,2,0) / sg(1:Ntz,0)
    dxi(1:Ntz,2,1) = zero                            ! not used
    dxi(1:Ntz,2,2) = dxi(1:Ntz,2,2) / sg(1:Ntz,0) - dxi(1:Ntz,2,0) * sg(1:Ntz,2) / sg(1:Ntz,0)
    dxi(1:Ntz,2,3) = dxi(1:Ntz,2,3) / sg(1:Ntz,0) - dxi(1:Ntz,2,0) * sg(1:Ntz,3) / sg(1:Ntz,0)
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(compxi)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine compxi

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
