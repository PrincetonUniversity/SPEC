!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Calculates metric quantities.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine me00ab( lvol, lss ) 

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  use constants, only : zero, one

  use numerical, only : small

  use fileunits, only : ounit

  use inputlist, only : Wme00ab

  use cputiming, only : Tme00ab

  use allglobal, only : myid, ncpu, cpus, &
                        DifferentiateGeometry, &
                        Lplasmaregion, Lvacuumregion, Mvol, &
                        mn, im, in, mne, ime, ine, &
                        isr, Nt, Nz, Ntz, trigm, trign, trigwk, efmn, ofmn, cfmn, sfmn, &
                        sg, ijreal, &
                        guvij, gvuij, guvmne, guvmno                                          

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS

  INTEGER, intent(in) :: lvol
  REAL   , intent(in) :: lss 
  
  INTEGER             :: Lcurvature, ifail, id, jd, ii, jj, kk, mi, ni, ideriv

  BEGIN(me00ab)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef DEBUG
  FATALMESS(me00ab, lvol.lt.1 .or. lvol.gt.Mvol, invalid lvol)
  FATALMESS(me00ab, abs(lss).gt.one, invalid local radial coordinate)
#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( DifferentiateGeometry%L ) then ; Lcurvature = 3 ; ideriv = 1
  else                               ; Lcurvature = 1 ; ideriv = 0
  endif

  WCALL(me00ab,co01aa,( lvol, lss, Lcurvature, Ntz, mn )) ! this returns guvij \equiv g_{u,v}; 17 Apr 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Lplasmaregion ) then
   
   ;              ; gvuij(1:Ntz, 0, 0, 0) = guvij(1:Ntz, 0, 0, ideriv)
   do id =  1, 3
    do jd = id, 3 ; gvuij(1:Ntz,id,jd, 0) = guvij(1:Ntz,id,jd, ideriv) / sg(1:Ntz,0) 
    enddo
   enddo
   
  else ! vacuum region; 17 Apr 13;
   
   gvuij(1:Ntz, 0, 0, 0) = zero ! ! this is associated with the helicity integral and is not required in the vacuum region; 18 Apr 13;
   
   do ii = 1, 3
    
    select case( ii )
    case( 1 ) ; jj = 2 ; kk = 3
    case( 2 ) ; jj = 3 ; kk = 1
    case( 3 ) ; jj = 1 ; kk = 2
    end select
    
    gvuij(1:Ntz,ii, 1, 0) = ( guvij(1:Ntz,jj, 2, ideriv) * guvij(1:Ntz,kk, 3, ideriv) - guvij(1:Ntz,jj, 3, ideriv) * guvij(1:Ntz,kk, 2, ideriv) ) / sg(1:Ntz,0)
    gvuij(1:Ntz,ii, 2, 0) = ( guvij(1:Ntz,jj, 3, ideriv) * guvij(1:Ntz,kk, 1, ideriv) - guvij(1:Ntz,jj, 1, ideriv) * guvij(1:Ntz,kk, 3, ideriv) ) / sg(1:Ntz,0)
    gvuij(1:Ntz,ii, 3, 0) = ( guvij(1:Ntz,jj, 1, ideriv) * guvij(1:Ntz,kk, 2, ideriv) - guvij(1:Ntz,jj, 2, ideriv) * guvij(1:Ntz,kk, 1, ideriv) ) / sg(1:Ntz,0)
    
   enddo

  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  ijreal(1:Ntz) = zero

  ifail = 0
  call tfft( Nt, Nz, gvuij(1:Ntz,0,0,0), ijreal(1:Ntz), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz), & ! not required in vacuum; 03 Apr 13;
             mne, ime(1:mne), ine(1:mne), guvmne(1:mne,0), guvmno(1:mne,0), cfmn(1:mne), sfmn(1:mne), ifail )

  ifail = 0
  call tfft( Nt, Nz, gvuij(1:Ntz,1,1,0), gvuij(1:Ntz,1,2,0), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz), &
             mne, ime(1:mne), ine(1:mne), guvmne(1:mne,1), guvmno(1:mne,1), guvmne(1:mne,2), guvmno(1:mne,2), ifail )

  ifail = 0
  call tfft( Nt, Nz, gvuij(1:Ntz,1,3,0), gvuij(1:Ntz,2,2,0), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz), &
             mne, ime(1:mne), ine(1:mne), guvmne(1:mne,3), guvmno(1:mne,3), guvmne(1:mne,4), guvmno(1:mne,4), ifail )

  ifail = 0
  call tfft( Nt, Nz, gvuij(1:Ntz,2,3,0), gvuij(1:Ntz,3,3,0), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz), &
             mne, ime(1:mne), ine(1:mne), guvmne(1:mne,5), guvmno(1:mne,5), guvmne(1:mne,6), guvmno(1:mne,6), ifail )
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  if( Wme00ab ) then
   if( Lvacuumregion ) then
    do ii = 1, mne ; mi = ime(ii) ; ni = ine(ii)
     if( ni.eq.0 ) cycle ! for debugging non-axisymmetric harmonics; 03 Apr 13;
     if( sum(abs(guvmne(ii,0:6)))+sum(abs(guvmno(ii,0:6))).gt.small ) then
     write(ounit,'("me00ab : ", 10x ," : lss="f9.4" ; ii=",i3," ; mi=",i3," ; ni=",i3," ; guvmne(ii,0:6)="99es13.5" ;")') lss, ii, mi, ni, guvmne(ii,0:6)
     write(ounit,'("me00ab : ", 10x ," : lss="f9.4" ; ii=",i3," ; mi=",i3," ; ni=",i3," ; guvmno(ii,0:6)="99es13.5" ;")') lss, ii, mi, ni, guvmno(ii,0:6)
     endif
    enddo
   endif ! end of if( Lvacuumregion ) ; 03 Apr 13;
  endif ! end of if( Wme00ab ) ; 03 Apr 13;
#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(me00ab)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine me00ab

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
