!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (perturbed RHS) ! Constructs the RHS for the boundary-perturbed Beltrami solution

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine ptbrhs(lvol, mn, lrad, dxi, gBu, rhs)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, one, two, half
  
  use numerical, only : small
  
  use fileunits, only : ounit
  
  use inputlist, only : Wmacros, Wptbrhs, Mpol
  
  use cputiming, only : Tptbrhs
  
  use allglobal, only : ncpu, myid, cpus, &
                        YESstellsym, NOTstellsym, &
                        im, in, Nt, Nz, Ntz, &
                        Lcoordinatesingularity, &
                        efmn, cfmn, ijreal, &
                        NAdof, dBdX, &
                        Lma, Lmb, Lmc, Lmd, Lme, Lmf
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  INTEGER, intent(in)  :: lvol, mn, lrad

  REAL, intent(in) :: dxi(1:Ntz,1:3,0:3), gBu(1:Ntz,1:3,0:3)

  REAL, intent(out) :: rhs(0:NAdof(lvol))
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  INTEGER              :: NN, ii, ll, jj, mi, ni, id, jd, ifail, innout, irz, issym
  
  REAL                 :: dBs(1:Ntz), dAt(1:Ntz), dAz(1:Ntz), dBns(1:mn), dBnc(1:mn), dAts(1:mn), dAtc(1:mn), dAzs(1:mn), dAzc(1:mn)
  
  BEGIN(ptbrhs)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  NN = NAdof(lvol) ! shorthand;
  
  ! clear the data
  rhs(0:NN) = zero
  
  ! one must have called coords using Lcurvature=2 then compxi BEFORE calling ptbrhs

  ii = dBdX%ii ; innout = dBdX%innout ; irz = dBdX%irz ; issym = dBdX%issym ! shorthand

  ! construct the perturbed boundary condition

  dBs(1:Ntz) = gBu(1:Ntz,2,0) * dxi(1:Ntz,1,2) + gBu(1:Ntz,3,0) * dxi(1:Ntz,1,3) - gBu(1:Ntz,1,1) * dxi(1:Ntz,1,0)

  dAt(1:Ntz) = - dxi(1:Ntz,1,0) * gBu(1:Ntz,3,0)
  dAz(1:Ntz) = + dxi(1:Ntz,1,0) * gBu(1:Ntz,2,0)

  ! ijreal, efmn, cfmn are junk
  ijreal = zero
  call tfft( Nt, Nz, dBs, ijreal, mn, im(1:mn), in(1:mn), dBnc, dBns, efmn, cfmn, ifail )
  

  call tfft( Nt, Nz, dAt, dAz, mn, im(1:mn), in(1:mn), dAtc, dAts, dAzc, dAzs, ifail )
write(ounit,*) lvol, im(ii), in(ii), 'dBns', dBns(ii), dAtc(ii), dAzc(ii)
  if (innout .eq. 0) then ! inner surface

    dBnc = -dBnc ; dBns = - dBns  ! flip the sign for inner surface
    dAtc = -dAtc ; dAts = - dAts ; dAzc = -dAzc ; dAzs = -dAzc

    if( YESstellsym ) then
 
      do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
        
          ;if( mi.ge.1 ) then   ; id = Lmb(lvol,ii)     ;                           ; rhs(id   ) = -dBns(ii) / mi
          ;elseif (ni.ge.1) then; id = Lma(lvol,ii)     ;                           ; rhs(id   ) = -dBns(ii) / ni
          ;endif
        !;                      ; id = Lma(lvol,ii)     ;                           ; rhs(id   ) = dAtc(ii)
        !;                      ; id = Lmb(lvol,ii)     ;                           ; rhs(id   ) = dAzc(ii)
        
    
      enddo ! end of do ii ;
      
    else ! NOTstellsym ;

      do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
      
        ;if( mi.ge.1 ) then   ; id = Lmb(lvol,ii)       ;                           ; rhs(id   ) = -dBns(ii) / mi
        ;                     ; id = Lmd(lvol,ii)       ;                           ; rhs(id   ) = +dBnc(ii) / mi
        ;elseif (ni.ge.1) then; id = Lma(lvol,ii)       ;                           ; rhs(id   ) = -dBns(ii) / ni
        ;                     ; id = Lmc(lvol,ii)       ;                           ; rhs(id   ) = +dBnc(ii) / ni
        ;endif
        
      enddo ! end of do ii ;
      
    endif ! end of if( YESstellsym ) ;

  else ! outer surface

    if( YESstellsym ) then
 
      do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
        
        ;if( ii.gt.1 ) then ; id = Lme(lvol,  ii)       ;                           ; rhs(id   ) =  dBns(ii)
        ;endif
    
      enddo ! end of do ii ;
      
    else ! NOTstellsym ;

      do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
      
        ;if( ii.gt.1 ) then ; id = Lme(lvol,ii)         ;                           ; rhs(id   ) =  dBns(ii)
        ;                   ; id = Lmf(lvol,ii)         ;                           ; rhs(id   ) =  dBnc(ii)
        ;endif
        
      enddo ! end of do ii ;
      
    endif ! end of if( YESstellsym ) ;

  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(ptbrhs)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine ptbrhs

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
