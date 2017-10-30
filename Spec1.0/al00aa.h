!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Allocates and initializes global arrays.

subroutine al00aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  use constants, only : zero

  use numerical, only : sqrtmachprec, vsmall

  use fileunits, only : ounit, wunit

  use inputlist

  use cputiming

  use allglobal

  use pjhamilto, only : Itor, Gpol, spmn

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS
  
  LOGICAL   :: exist
  INTEGER   :: ieo, innout
  INTEGER   :: idof, ij, ik, ip, iq, jk, io, ll, ii, imn, jmn, kmn, lmn, llc, imnc, ipq, mnpq, mm, nn, ifail, iofe, nc, ideriv, vvol, pp, qq
  INTEGER   :: mi, ni, mj, nj, mk, nk, mimj, ninj, ksgn, mkmj, nknj
  INTEGER   :: jj, kk
  INTEGER   :: lMpol, lNtor, lNi
  INTEGER   :: madd, nadd, msub, nsub, lNZ
  REAL      :: teta, zeta, arg, lss
  INTEGER   :: llmodnp
  INTEGER   :: itype, lquad, id01bcf, maxIquad, maxLrad, jquad
  REAL      :: aa, bb, cc, dd, exponent
  
  BEGIN(al00aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! set / allocate internal variables that depend on namelist/physicslist/
  
  FATALMESS(al00aa, Nfp.eq.0, illegal division)

  pi2nfp         = pi2 / Nfp
  
  pi2pi2nfp      = pi2 * pi2nfp
  pi2pi2nfphalf  = pi2 * pi2nfp * half
  pi2pi2nfpquart = pi2 * pi2nfp * quart
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! local  number of degrees of freedom in geometry;
!                                                          Rbc  Zbs    Rbs    Zbc
  if( YESstellsym .and. Igeometry.lt.3 ) lgeometricaldof = mn      
  if( YESstellsym .and. Igeometry.ge.3 ) lgeometricaldof = mn + mn-1
  if( NOTstellsym .and. Igeometry.lt.3 ) lgeometricaldof = mn        + mn-1
  if( NOTstellsym .and. Igeometry.ge.3 ) lgeometricaldof = mn + mn-1 + mn-1 + mn ! 11 Feb 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! global number of degrees of freedom in geometry;
  
  Ngeometricaldof = ( Mvol-1 ) * lgeometricaldof
  
  if( Wal00aa ) then ; cput = GETTIME ; write(ounit,'("al00aa : ",f10.2," : myid=",i3," ; Ngeometricaldof="i9" ;")') cput-cpus, myid, Ngeometricaldof
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! rotational transform on interfaces;
  
  do vvol = 0, Nvol ! 23 Jan 13;
   if( ql(vvol).eq.0 .and. qr(vvol).eq.0 ) then ! "noble" selection is invalid; 04 Dec 14;
    ;                                    ; iota(vvol) = iota(vvol) ! provide iota directly; 04 Dec 14;
   else
    ;                                    ; iota(vvol) = ( pl(vvol) + goldenmean * pr(vvol) ) / ( ql(vvol) + goldenmean * qr(vvol) ) ! noble transform; 04 Dec 14;
   endif
   if( lq(vvol).eq.0 .and. rq(vvol).eq.0 ) then ! "noble" selection is invalid; 04 Dec 14;
    if( abs(oita(vvol)).gt.vsmall ) then ; oita(vvol) = oita(vvol) ! provide oita directly; 04 Dec 14;
    else                                 ; oita(vvol) = iota(vvol)
    endif
   else
    ;                                    ; oita(vvol) = ( lp(vvol) + goldenmean * rp(vvol) ) / ( lq(vvol) + goldenmean * rq(vvol) ) ! noble transform; 04 Dec 14;
   endif
  enddo
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! enclosed fluxes in each plasma volume;
  
  RALLOCATE(dtflux,(1:Nvol))         ! \delta \psi_{toroidal} in each annulus;
  RALLOCATE(dpflux,(1:Nvol))         ! \delta \psi_{poloidal} in each annulus;  
  ;           ; dtflux(1) = tflux(1) ! it is the differential toroidal flux that controls the transform; tflux was normalized in readin; 16 Jan 15;
  select case( Igeometry )
  case( 1   ) ; dpflux(1) = pflux(1) ! Cartesian              ; 29 Apr 14;
  case( 2:4 ) ; dpflux(1) = zero     ! cylindrical or toroidal; 29 Apr 14; ! poloidal flux in inner torus is irrelevant; 01 Jul 14;
  case default ; FATALMESS(al00aa, .true., illegal Igeometry)
  end select

 !if( Igeometry.eq.1 )  ! 16 Jan 15;
 !if( Igeometry.gt.1 )  ! 16 Jan 15;
  
  dtflux(2:Nvol) = tflux(2:Nvol) - tflux(1:Nvol-1) ! 03 Apr 13;
  dpflux(2:Nvol) = pflux(2:Nvol) - pflux(1:Nvol-1) ! it is the differential poloidal flux that controls the transform;
  
#ifdef DEBUG
  FATALMESS(al00aa, abs(pi2).lt.vsmall, nonsense) ! just checking division; 04 Dec 14;
#endif

  dtflux(1:Nvol) = dtflux(1:Nvol) * phiedge / pi2 ! 20 Jun 14;
  dpflux(1:Nvol) = dpflux(1:Nvol) * phiedge / pi2 ! 18 Jul 14;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! spectral weights; 18 Jul 14;
  
  select case( Igeometry ) 
  case( 1:2 )  ; 
  case( 3:4 )  ; RALLOCATE(sweight,(1:Mvol))
   ;           ; sweight(2:Mvol) = + upsilon * tflux(2:Mvol)**wpoloidal ! weight factor for poloidal length constraint; 04 Dec 14;
  case default ; FATALMESS(al00aa, .true., illegal Igeometry)
  end select
  
! if ( Igeometry.gt.2 ) then

!  RALLOCATE(sweight,(1:Mvol))
  
!  sweight(2:Mvol) = + upsilon * tflux(2:Mvol)**wpoloidal ! weight factor for poloidal length constraint; 04 Dec 14;
   
!#ifdef DEBUG
!   if( Wal00aa .and. myid.eq.0 ) write(ounit,'("al00aa : ", 10x ," : wpoloidal="f5.2" ; sweight="999es13.5)') wpoloidal, sweight(1:Mvol)
!#endif
   
!  endif ! end of if( Igeometry.gt.2 ) ; 14 Jan 15;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! derivatives of Chebyshev polynomials at the inner and outer interfaces;
  
!latex \item The internal array, \verb+TTll(0:Lrad,0:1,0:1)+ gives the Chebyshev polynomials, and their derivatives, evaluated at $s=-1$ and $s=+1$ and is constructed as follows:
!latex       \begin{itemize}
!latex       \item \verb+TTll(ll,0:1,0)+ $\equiv T_l       (s_i) = s^l    $ for $s_i=-1,+1$; 
!latex       \item \verb+TTll(ll,0:1,1)+ $\equiv T_l^\prime(s_i) = s^l l^2$ for $s_i=-1,+1$; 
!latex       \end{itemize}

  RALLOCATE(TTll,(0:maxval(Lrad(1:Mvol)),0:1,0:1))
  
  do innout = 0, 1
   
   lss = two * innout - one
   
   do ll = 0, maxval(Lrad(1:Mvol)) ; TTll(ll,innout,0) = lss**(ll  )         !               Chebyshev polynomial at interface; 24 Apr 13;
    ;                              ; TTll(ll,innout,1) = lss**(ll+1) * ll**2 ! derivative of Chebyshev polynomial at interface; 24 Apr 13; corrected; 20 Jun 14;
   enddo
   
!#ifdef DEBUG ! 14 Jan 15;
!   if( Wal00aa ) then ! 14 Jan 15;
!    write(ounit,'("al00aa : ", 10x ," : TT(:,"i1",0)="99f5.0)') innout, TTll(0:maxval(Lrad(1:Mvol)),innout,0) ! 14 Jan 15;
!    write(ounit,'("al00aa : ", 10x ," : TT(:,"i1",1)="99f5.0)') innout, TTll(0:maxval(Lrad(1:Mvol)),innout,1) ! 14 Jan 15;
!   endif ! 14 Jan 15;
!#endif ! 14 Jan 15;
   
  enddo ! end of do innout = 0, 1 ; 20 Jun 14;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! introduce some error flags for convenience;
  
  LALLOCATE(ImagneticOK,(1:Mvol)) ! error flag indicating successful construction of Beltrami fields in each volume; !  5 Feb 13;
  
  Lhessianallocated = .false.
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! identify Fourier harmonics; used in ma00aa; 29 Jan 13;

  IALLOCATE(guvmnks,(1:mn,1:mn,0:1))
  IALLOCATE(guvmnka,(1:mn,1:mn,0:1))
  
  IALLOCATE(guvmnk,(1:mn,0:1))
  
  do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
   
   do kk = 1, mne ; mk = ime(kk) ; nk = ine(kk)
    
    if( mk.eq.mi .and. nk.eq.ni ) then
     if( ii.eq.1 ) then ; guvmnk(ii,0:1) = (/ kk, 1 /) ! this is the (m,n)=(0,0) mode; 02 Apr 13;
     else               ; guvmnk(ii,0:1) = (/ kk, 2 /) ! denominator in Fourier extraction; 02 Apr 13;
     endif
    endif
    
   enddo ! end of do kk = 1, mne; 03 Apr 13;
   
   do jj = 1, mn ; mj = im(jj) ; nj = in(jj)
    
    mimj = mi - mj ; ninj = ni - nj
    
    do kk = 1, mne ; mk = ime(kk) ; nk = ine(kk)
     
     if( mimj.gt.0 .or. ( mimj.eq.0 .and. ninj.ge.0 ) ) then
      
      if( mk.eq. mimj .and. nk.eq. ninj ) then
       if( mk.eq.0 .and. nk.eq.0 ) then ; guvmnks(ii,jj,0:1) = (/ kk ,   1 /)
       else                             ; guvmnks(ii,jj,0:1) = (/ kk ,   2 /)
       endif
      endif
      
     else
      
      if( mk.eq.-mimj .and. nk.eq.-ninj ) then
       if( mk.eq.0 .and. nk.eq.0 ) then ; FATALMESS(al00aa,.true.,how can this happen)
       else                             ; guvmnks(ii,jj,0:1) = (/ kk , - 2 /)
       endif
      endif
      
     endif
     
    enddo ! end of do kk; 29 Jan 13;
    
    mimj = mi + mj ; ninj = ni + nj
    
    do kk = 1, mne ; mk = ime(kk) ; nk = ine(kk)
     
     if( mimj.gt.0 .or. ( mimj.eq.0 .and. ninj.ge.0 ) ) then
      
      if( mk.eq.mimj .and. nk.eq.ninj ) then
       if( mk.eq.0 .and. nk.eq.0 ) then ; guvmnka(ii,jj,0:1) = (/ kk ,   1 /)
       else                             ; guvmnka(ii,jj,0:1) = (/ kk ,   2 /) ! denominator;  8 Feb 13;
       endif
      endif
      
     else
      
      FATALMESS(al00aa, .true., how can this happen)
      
     endif ! end of if( mimj.gt.0 .or. ( mimj.eq.0 .and. ninj.ge.0 ) ) then; 20 Jun 14;
     
    enddo ! end of do kk; 29 Jan 13;
    
   enddo ! end of do jj; 29 Jan 13;
   
  enddo ! end of do ii; 29 Jan 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Igeometry.eq.2 ) then ! standard cylindrical; 04 Dec 14;
   
   IALLOCATE(djkp,(1:mn,1:mn)) ! only used in vo00aa; trignometric identities; 04 Dec 14;
   IALLOCATE(djkm,(1:mn,1:mn)) ! only used in vo00aa; trignometric identities; 04 Dec 14;
   
   do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
    do jj = 1, mn ; mj = im(jj) ; nj = in(jj)
     if( mi-mj.eq.0 .and. ni-nj.eq.0 ) djkp(ii,jj) = 1
     if( mi+mj.eq.0 .and. ni+nj.eq.0 ) djkm(ii,jj) = 1
    enddo
   enddo
   
  endif ! end of if( Igeometry.eq.2 ) ; 04 Dec 14;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  IALLOCATE(iotakkii,(1:mn      )) ! used to identify matrix elements in straight-field-line angle transformation; will be redundant after real-space sparse solution implemented;
  
  IALLOCATE(iotaksub,(1:mn,1:mns))
  IALLOCATE(iotaksgn,(1:mn,1:mns))
  IALLOCATE(iotakadd,(1:mn,1:mns))
  
  do kk = 1, mn ; mk = im(kk) ; nk = in(kk)
   
   
   do ii = 1, mns ; mi = ims(ii) ; ni = ins(ii)
    
    if( mk.eq.mi .and. nk.eq.ni ) iotakkii(kk) = ii
    
   enddo
   
   
   do jj = 1, mns ; mj = ims(jj) ; nj = ins(jj)
    
    
    mkmj = mk - mj ; nknj = nk - nj
    
    do ii = 1, mns ; mi = ims(ii) ; ni = ins(ii)
     
     if( mkmj.gt.0 .or. ( mkmj.eq.0 .and. nknj.ge.0 ) ) then
      
      if( mi.eq. mkmj .and. ni.eq. nknj ) then ; iotaksub(kk,jj) = ii ; iotaksgn(kk,jj) =  1
      endif
      
     else
      
      if( mi.eq.-mkmj .and. ni.eq.-nknj ) then ; iotaksub(kk,jj) = ii ; iotaksgn(kk,jj) = -1
      endif
      
     endif
     
    enddo ! end of do ii; 30 Jan 13;
    
    
    mkmj = mk + mj ; nknj = nk + nj
    
    do ii = 1, mns ; mi = ims(ii) ; ni = ins(ii)
     
     if( mi.eq. mkmj .and. ni.eq. nknj ) then ; iotakadd(kk,jj) = ii
     endif
     
    enddo ! end of do ii; 29 Jan 13;
    
    
   enddo ! end of do jj; 29 Jan 13;
   
   
  enddo ! end of do kk; 29 Jan 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! set / allocate internal variables that depend on namelist/numericlist/
  
  IALLOCATE(Iquad,(1:Mvol)) ! 16 Jan 13;
  
  do vvol = 1, Mvol
   
  !if( Igeometry.eq.1 .or.  vvol.gt.1 ) Lcoordinatesingularity = .false.
  !if( Igeometry.gt.1 .and. vvol.eq.1 ) Lcoordinatesingularity = .true.
   
   if( Igeometry.eq.1 .or. vvol.gt.1 ) then ; Lcoordinatesingularity = .false. ! 14 Jan 15;
   else                                     ; Lcoordinatesingularity = .true.  ! 14 Jan 15;
   endif
   
   if( Nquad.gt.0 ) then ;            Iquad(vvol) =                         Nquad
   else                 
    if(      Lcoordinatesingularity ) Iquad(vvol) = Mpol + 2 * Lrad(vvol) - Nquad ! NEED TO REVISE REGULARIZATION FACTORS; 26 Feb 13;
    if( .not.Lcoordinatesingularity ) Iquad(vvol) =        2 * Lrad(vvol) - Nquad
   endif
   
  enddo ! end of do vvol; 18 Feb 13;
  
  
  maxIquad = maxval(Iquad(1:Mvol))
  maxLrad  = maxval(Lrad(1:Mvol))
  
  RALLOCATE(gaussianweight,(1:maxIquad,1:Mvol))
  RALLOCATE(gaussianabscissae,(1:maxIquad,1:Mvol))
  
  RALLOCATE(gchebyshev,(0:maxLrad,0:2,1:maxIquad,1:Mvol))
  
  do vvol = 1, Mvol
   
   itype = 0 ; aa = -one ; bb = +one ; cc = zero ; dd = zero ; lquad = Iquad(vvol)
   
   id01bcf = 1
   call D01BCF( itype, aa, bb, cc, dd, lquad, gaussianweight(1:lquad,vvol), gaussianabscissae(1:lquad,vvol), id01bcf ) ! sets gaussian weights & abscissae;
   
   cput= GETTIME
   select case( id01bcf ) !                                                        123456789012345
   case( 0 )    ; if( Wal00aa ) write(ounit,1000) cput-cpus, myid, vvol, id01bcf, "success        ", gaussianabscissae(1:lquad,vvol)
   case( 1 )    ;               write(ounit,1000) cput-cpus, myid, vvol, id01bcf, "failed         ", gaussianabscissae(1:lquad,vvol)
   case( 2 )    ;               write(ounit,1000) cput-cpus, myid, vvol, id01bcf, "input error    ", gaussianabscissae(1:lquad,vvol)
   case( 3 )    ;               write(ounit,1000) cput-cpus, myid, vvol, id01bcf, "input error    ", gaussianabscissae(1:lquad,vvol)
   case( 4 )    ;               write(ounit,1000) cput-cpus, myid, vvol, id01bcf, "weight overflow", gaussianabscissae(1:lquad,vvol)
   case( 5 )    ;               write(ounit,1000) cput-cpus, myid, vvol, id01bcf, "weight zero    ", gaussianabscissae(1:lquad,vvol)
   case( 6 )    ;               write(ounit,1000) cput-cpus, myid, vvol, id01bcf, "failed         ", gaussianabscissae(1:lquad,vvol)
   case default ;               write(ounit,1000) cput-cpus, myid, vvol, id01bcf, "weird          ", gaussianabscissae(1:lquad,vvol)
   end select
   ;              if( Wal00aa ) write(ounit,1001)                                                    gaussianweight(1:lquad,vvol)
   
1000 format("al00aa : ",f10.2," : myid=",i3," ; lvol=",i3," ; id01bcf=",i3," ; "a15" ; abscissae ="99f10.6)
1001 format("al00aa : ", 10x ," :      "3x"        "3x"           "3x"   "15x" ; weights   ="99f10.6)
   
   
   
   do jquad = 1, lquad ! loop over radial sub-sub-grid (numerical quadrature);
    lss = gaussianabscissae(jquad,vvol)
    gchebyshev(0,0:1,jquad,vvol) = (/ one, zero /) ! Chebyshev initialization; 16 Jan 13;
    gchebyshev(1,0:1,jquad,vvol) = (/ lss,  one /) ! Chebyshev initialization; 16 Jan 13;   
    do ll = 2, Lrad(vvol)
     gchebyshev(ll,0:1,jquad,vvol) = (/ two * lss * gchebyshev(ll-1,0,jquad,vvol)                                             - gchebyshev(ll-2,0,jquad,vvol) , &
                                        two       * gchebyshev(ll-1,0,jquad,vvol) + two * lss * gchebyshev(ll-1,1,jquad,vvol) - gchebyshev(ll-2,1,jquad,vvol) /)
    enddo
   enddo
   
  enddo ! end of do vvol;  7 Mar 13; 
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! set / allocate internal variables that depend on namelist/linearlist/
  
  LBeltramiSeQuad = .false.
  LBeltramiNewton = .false.
  LBeltramiLinear = .false.
  
  if( LBeltrami.eq.1 .or. LBeltrami.eq.3 .or. LBeltrami.eq.5 .or. LBeltrami.eq.7 ) LBeltramiSeQuad = .true.
  if( LBeltrami.eq.2 .or. LBeltrami.eq.3 .or. LBeltrami.eq.6 .or. LBeltrami.eq.7 ) LBeltramiNewton = .true.
  if( LBeltrami.eq.4 .or. LBeltrami.eq.5 .or. LBeltrami.eq.6 .or. LBeltrami.eq.7 ) LBeltramiLinear = .true.
  
  if( myid.eq.0 ) then
   cput = GETTIME
   write(ounit,'("al00aa : ",f10.2," : LBeltramiSeQuad="L2" , LBeltramiNewton="L2" , LBeltramiLinear="L2" ;")')cput-cpus, LBeltramiSeQuad, LBeltramiNewton, LBeltramiLinear
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! set internal variables that depend on namelist/nonlinearlist/
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! set internal variables that depend on namelist/diagnosticslist/
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! set internal variables that depend on namelist/screenlist/
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RALLOCATE(expmmnn,(1:mn)) ! exponential weight on force-imbalance harmonics; 04 Dec 14;

  do ii = 1, mn ; mi = im(ii) ; ni = abs(in(ii)) / Nfp ; exponent = escale * ( mi*mi + ni*ni ) ; expmmnn(ii) = one / exp( exponent )
  enddo
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! set spectral condensation weight factors;
  
  RALLOCATE(mpnq,(1:mn))
  
  do ii = 1, mn ; mi= im(ii) ; ni = abs(in(ii)) / Nfp
   
   if( mi.eq.0 ) then
    
   !if( ni.eq.0 ) then ! 04 Dec 14;
     
     mpnq(ii) = zero
     
   !else ! 04 Dec 14;
     
   ! if( qcondense.gt.zero ) then ; mpnq(ii) = ni**qcondense ! 04 Dec 14;
   ! else                         ; mpnq(ii) = zero ! 04 Dec 14;
   ! endif ! 04 Dec 14;
     
   !endif ! 04 Dec 14;
    
   else ! mi.ne.0 ; 11 Aug 14;
    
    mpnq(ii) = mi**pcondense
    
   endif ! matches if( mi.eq.0 ) ; 11 Aug 14;

  enddo ! end of do ii; 08 Nov 13;
   
#ifdef DEBUG
  if( Wal00aa ) then
   if( Igeometry.ge.3 .and. myid.eq.0 ) then
    cput = GETTIME
    write(ounit,'("al00aa : ", 10x ," : ")')
    write(ounit,'("al00aa : ",f10.2," : (m,n) = ("i4" ,"i4" ) ; mpnq="f13.1" ;")')( cput-cpus, im(ii), in(ii), mpnq(ii), ii = 1, mn ) 
   endif
  endif
#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! allocate vector potential and Beltrami arrays;
  
  maxLrad = maxval(Lrad(1:Mvol))
  RALLOCATE(lchebyshev,(0:maxLrad,0:2))
  
  ALLOCATE(Ate,(1:Mvol,-1:2,1:mn)) ! recall that this is type:sub-grid; 31 Jan 13;
  ALLOCATE(Aze,(1:Mvol,-1:2,1:mn)) ! recall that this is type:sub-grid; 31 Jan 13;
  ALLOCATE(Ato,(1:Mvol,-1:2,1:mn)) ! recall that this is type:sub-grid; 31 Jan 13;
  ALLOCATE(Azo,(1:Mvol,-1:2,1:mn)) ! recall that this is type:sub-grid; 31 Jan 13;
  
  IALLOCATE(Fso,(1:Nvol,1:mn)) ! 25 Jan 13;
  IALLOCATE(Fse,(1:Nvol,1:mn)) ! 25 Jan 13;
  
!latex \item The total number of degrees of freedom in the Beltrami fields in each volume is given by
!latex       \begin{itemize}
!latex       \item if \verb+Lstelsym.eq.0+ and \verb+Lcoordinatesingularity.eq.T+
!latex
!latex             \verb+Nmagneticdof+
!latex             $ = $
!latex       \item if \verb+Lstelsym.eq.0+ and \verb+Lcoordinatesingularity.eq.F+
!latex
!latex             \verb+Nmagneticdof+
!latex             $ = \underbrace{( mn + mn-1 )}_{\textrm{\small Fourier harmonics}} \times \; 2 \; \times 
!latex             \underbrace{(L-1)}_{\textrm{\small Chebyshev}} + \; \underbrace{( mn + mn-1 )}_{\textrm{\small surface potential}}$
!latex       \item if \verb+Lstelsym.eq.1+ and \verb+Lcoordinatesingularity.eq.T+
!latex
!latex             \verb+Nmagneticdof+
!latex             $ = $
!latex       \item if \verb+Lstelsym.eq.1+ and \verb+Lcoordinatesingularity.eq.F+
!latex
!latex             \verb+Nmagneticdof+
!latex             $ = \underbrace{  mn         }_{\textrm{\small Fourier harmonics}} \times \; 2 \; \times 
!latex             \underbrace{(L-1)}_{\textrm{\small Chebyshev}} + \; \underbrace{(      mn-1 )}_{\textrm{\small surface potential}}$
!latex       \end{itemize}
  
  IALLOCATE(Nmagneticdof,(1:Mvol)) ! Beltrami degrees of freedom in each annulus;
  
  do vvol = 1, Mvol
   
  !if( Igeometry.eq.1 .or.  vvol.gt.1 ) Lcoordinatesingularity = .false.
  !if( Igeometry.gt.1 .and. vvol.eq.1 ) Lcoordinatesingularity = .true.
   
   if( Igeometry.eq.1 .or. vvol.gt.1 ) then ; Lcoordinatesingularity = .false. ! 14 Jan 15;
   else                                     ; Lcoordinatesingularity = .true.  ! 14 Jan 15;
   endif
   
   if( vvol.le.Nvol ) Lplasmaregion = .true.
   if( vvol.gt.Nvol ) Lplasmaregion = .false.
   
   Lvacuumregion = .not. Lplasmaregion
   
   
   if( Lplasmaregion ) then
    
    if( YESstellsym .and.      Lcoordinatesingularity ) Nmagneticdof(vvol) = ( mn        ) * 2 * ( Lrad(vvol)-1 ) + (        mn-1 ) + mn-Ntor
    if( YESstellsym .and. .not.Lcoordinatesingularity ) Nmagneticdof(vvol) = ( mn        ) * 2 * ( Lrad(vvol)-1 ) + (        mn-1 )
    if( NOTstellsym .and.      Lcoordinatesingularity ) Nmagneticdof(vvol) = ( mn + mn-1 ) * 2 * ( Lrad(vvol)-1 ) + ( mn-1 + mn-1 ) + mn-Ntor + mn-1-Ntor
    if( NOTstellsym .and. .not.Lcoordinatesingularity ) Nmagneticdof(vvol) = ( mn + mn-1 ) * 2 * ( Lrad(vvol)-1 ) + ( mn-1 + mn-1 )
    
   else ! Lvacuumregion;  19 Apr 13;
    
    if( YESstellsym ) Nmagneticdof(vvol) = (mn-1     ) * ( Lrad(vvol)+1 )
    if( NOTstellsym ) Nmagneticdof(vvol) = (mn-1+mn-1) * ( Lrad(vvol)+1 )
    
   endif ! end of if( Lplasmaregion ) ; 19 Apr 13;

   
   
   do ii = 1, mn ! loop over Fourier harmonics;
    
    do ideriv = -1, 2 ! loop over derivatives; 14 Jan 13;
     
     RALLOCATE(Ate(vvol,ideriv,ii)%s,(0:Lrad(vvol))) ! covariant coefficients of vector potential;
     RALLOCATE(Aze(vvol,ideriv,ii)%s,(0:Lrad(vvol))) ! covariant coefficients of vector potential;
!    if( NOTstellsym ) then ! shall always allocate these and write to output file; 19 Apr 13;
     RALLOCATE(Ato(vvol,ideriv,ii)%s,(0:Lrad(vvol))) ! covariant coefficients of vector potential;
     RALLOCATE(Azo(vvol,ideriv,ii)%s,(0:Lrad(vvol))) ! covariant coefficients of vector potential;
!    endif
     
    enddo ! end of do ideriv;
    
    IALLOCATE(Ate(vvol,0,ii)%i,(0:Lrad(vvol))) ! degree of freedom index; 17 Jan 13;
    IALLOCATE(Aze(vvol,0,ii)%i,(0:Lrad(vvol)))
    if( NOTstellsym ) then ! the degree-of-freedom identification integers are not required for the non-stellarator-symmetric terms; 19 Apr 13;
    IALLOCATE(Ato(vvol,0,ii)%i,(0:Lrad(vvol))) ! degree of freedom index;
    IALLOCATE(Azo(vvol,0,ii)%i,(0:Lrad(vvol)))
    endif
    
   enddo ! end of do ii;
   
   
#ifdef OLDNAG

   do ii = 1, mn
    
    ideriv = 0
    ifail = 0 ; call G05FAF(-sqrtmachprec, sqrtmachprec, Lrad(vvol)+1, Ate(vvol,ideriv,ii)%s(0:Lrad(vvol)), ifail ) ! initialize vector potential with small random;
    ifail = 0 ; call G05FAF(-sqrtmachprec, sqrtmachprec, Lrad(vvol)+1, Aze(vvol,ideriv,ii)%s(0:Lrad(vvol)), ifail ) ! initialize vector potential with small random;     
   !if( NOTstellsym ) then
    ifail = 0 ; call G05FAF(-sqrtmachprec, sqrtmachprec, Lrad(vvol)+1, Ato(vvol,ideriv,ii)%s(0:Lrad(vvol)), ifail ) ! initialize vector potential with small random;
    ifail = 0 ; call G05FAF(-sqrtmachprec, sqrtmachprec, Lrad(vvol)+1, Azo(vvol,ideriv,ii)%s(0:Lrad(vvol)), ifail ) ! initialize vector potential with small random;    
   !endif
     
   enddo ! end of do ii; 26 Feb 13;

#endif
   
!  WARNING: this "random" initilization for vector potential does not satisfy interface constraint;


   if( Lplasmaregion ) then

    select case( Linitgues )
    case( 0 )    ; 
    case( 1 )    ; Ate(vvol,0,1)%s(0:1) = dtflux(vvol) * half ! this is an integrable approximation; NEEDS CHECKING; 26 Feb 13;
     ;           ; Aze(vvol,0,1)%s(0:1) = dpflux(vvol) * half ! this is an integrable approximation; NEEDS CHECKING; 26 Feb 13;
    case( 2 )    ;                                            ! will call ra00aa below to read initial vector potential from file;
    case default ; 
     ;           ; FATALMESS(al00aa, .true., given value of Linitgues is not supported)
    end select
    
#ifdef DEBUG
    if( Llatex ) then
    write(wunit,'("\subsection{{degrees-of-freedom in Beltrami field : plasma region : volume="i3"; [al00aa];}}")') vvol
    write(wunit,'("The degrees-of-freedom are labelled as follows: ")')
    write(wunit,'("\begin{{enumerate}}")')
    endif
#endif
    
    idof = 0 ! degree of freedom index; reset to 0 in each annulus;
   
    do ii = 1, mn ; mi = im(ii) ; ni = in(ii) ! Fourier harmonics; 17 Jan 13;
     
     do ll = 0, Lrad(vvol)-2
     
      idof = idof + 1 ; Ate(vvol,0,ii)%i(ll) = idof
#ifdef DEBUG
      if( Llatex ) write(wunit,'("\item $(m,n)=("i3","i3")$, $l="i2"$ : $\boldxi_{{"i6"}} \equiv A_{{\t,"i3","i2"}}$")') mi, ni, ll, idof, ii, ll
#endif
      idof = idof + 1 ; Aze(vvol,0,ii)%i(ll) = idof
#ifdef DEBUG
      if( Llatex ) write(wunit,'("\item $(m,n)=("i3","i3")$, $l="i2"$ : $\boldxi_{{"i6"}} \equiv A_{{\z,"i3","i2"}}$")') mi, ni, ll, idof, ii, ll
#endif
      if( NOTstellsym ) then
      if( ii.eq.1 ) then
     !idof = idof + 0 ; Ato(vvol,0,ii)%i(ll) = 0
     !idof = idof + 0 ; Azo(vvol,0,ii)%i(ll) = 0
      else
      idof = idof + 1 ; Ato(vvol,0,ii)%i(ll) = idof
      idof = idof + 1 ; Azo(vvol,0,ii)%i(ll) = idof
      endif
      endif

     enddo ! end of do ll; 17 Jan 13;
     
     if( ii.eq.1 ) then
    !idof = idof + 0 ; Fso(vvol,ii) = 0
     else               
     idof = idof + 1 ; Fso(vvol,ii) = idof
#ifdef DEBUG
     if( Llatex ) write(wunit,'("\item $(m,n)=("i3","i3")$, $l=xx$ : $\boldxi_{{"i6"}} \equiv f_{{"i3"}}$")') mi, ni, idof, ii
#endif
     endif
     if( NOTstellsym ) then
     if( ii.eq.1 ) then
    !idof = idof + 0 ; Fse(vvol,ii) = 0 ! constant part is irrelevant; 26 Feb 13;
     else               
     idof = idof + 1 ; Fse(vvol,ii) = idof
     endif
     endif
     
    enddo ! end of do ii; 25 Jan 13;
    
    if( Lcoordinatesingularity ) then ! need to add additional degrees of freedom; 22 Apr 13;
     
     ll = Lrad(vvol) - 1
     
     do ii = 1, mn ; mi = im(ii) ; ni = in(ii) ! Fourier harmonics; 17 Jan 13;
      
     !if( mi.ne.0 .or. ni.eq.0 ) then ! 15 Jan 15; ! additional freedom; 15 Jan 15;
      if( mi.ne.0 .or. ii.eq.1 ) then ! 15 Jan 15; ! additional freedom; 15 Jan 15;
      idof = idof + 1 ; Aze(vvol,0,ii)%i(ll) = idof
#ifdef DEBUG
      if( Llatex ) write(wunit,'("\item $(m,n)=("i3","i3")$, $l="i2"$ : $\boldxi_{{"i6"}} \equiv A_{{\z,"i3","i2"}}$")') mi, ni, ll, idof, ii, ll
#endif
      if( NOTstellsym .and. ii.gt.1 ) then
      idof = idof + 1 ; Azo(vvol,0,ii)%i(ll) = idof
      endif
      endif

     enddo ! end of do ii; 20 Feb 13;
    
    endif ! end of if( Lcoordinatesingularity ) ; 20 Feb 13;

#ifdef DEBUG
    if( Llatex ) write(wunit,'("\end{{enumerate}}")')
#endif

   else ! matches if( Lplasmaregion ) ; 19 Apr 13;

    ;                 Ate(vvol,0, 1)%i(0:Lrad(vvol)) = 0 ! these are not degrees of freedom; 03 Apr 13;
    if( NOTstellsym ) Ato(vvol,0, 1)%i(0:Lrad(vvol)) = 0 ! these are not degrees of freedom; 03 Apr 13;

    idof = 0

    do ii = 2, mn ; mi = im(ii) ; ni = in(ii)
     
     do ll = 0, Lrad(vvol)
      
      idof = idof + 1 ; Ate(vvol,0,ii)%i(ll) = idof ! BEWARE: e-o convention is inverted; this is the odd  harmonic; 22 Apr 13;
      if( NOTstellsym ) then
      idof = idof + 1 ; Ato(vvol,0,ii)%i(ll) = idof ! BEWARE: e-o convention is inverted; this is the even harmonic; 22 Apr 13;
      endif

     enddo ! end of do ll; 22 Apr 13;
     
    enddo ! end of do ii; 22 Apr 13;
    
   endif ! end of if( Lplasmaregion) ; 19 Apr 13;

   FATALMESS(al00aa, idof.ne.Nmagneticdof(vvol), need to count degrees of freedom more carefully)
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  enddo ! end of do vvol = 1, Nvol loop;
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

3333 format("al00aa : ", 10x ," : myid=",i3," ; vvol=",i3," ; ii=",i3,", (m,n)=(",i3," ,",i3," ), ll="i2" : "a3" : idof="i4" ;")
3334 format("al00aa : ", 10x ," : myid=",i3," ; vvol=",i3," ; ii=",i3,", (m,n)=(",i3," ,",i3," ),    "2x" : "a3" : idof="i4" ;")

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Linitgues.eq.2 ) then ! read initial guess for Beltrami field from file; 02 Jan 15;
   WCALL( al00aa, ra00aa, ('R') ) 
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! optional screen output; 14 Jan 13;
  
  if( myid.eq.0 ) then ! 17 Oct 12;
   cput = GETTIME
   write(ounit,'("al00aa : ", 10x ," : ")')         
   write(ounit,'("al00aa : ",f10.2," : Nquad="i4" ; mn="i5" ; Ngeometricaldof="i6" ; Nmagneticdof="16(i6",")" ...")') &
 cput-cpus, Nquad, mn, Ngeometricaldof, Nmagneticdof(1:min(Mvol,16))
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! Fourier transforms;

  Nt = max( Ndiscrete*4*Mpol, 1 ) ; Nz = max( Ndiscrete*4*Ntor, 1 ) ; Ntz = Nt*Nz ; soNtz = one / sqrt( one*Ntz ) ! discrete resolution; exaggerated discrete resolution;

  ;                  ; hNt = Nt / 2
  if( Nz.gt.1 ) then ; hNz = Nz / 2
  else               ; hNz = 0
  endif
  
  if( myid.eq.0 ) then ! 17 Oct 12;
   cput = GETTIME
   write(ounit,'("al00aa : ", 10x ," : ")')         
   write(ounit,'("al00aa : ",f10.2," : Nt="i6" ; Nz="i6" ; Ntz="i9" ;")') cput-cpus, Nt, Nz, Ntz
  endif

  RALLOCATE(iRij,(1:Ntz,0:Mvol)) ! interface geometry in real space; ! 18 Jul 14;
  RALLOCATE(iZij,(1:Ntz,0:Mvol)) ! 

  RALLOCATE(dRij,(1:Ntz,1:Mvol)) ! interface geometry in real space; poloidal derivative; ! 18 Jul 14;
  RALLOCATE(dZij,(1:Ntz,1:Mvol))

  RALLOCATE(tRij,(1:Ntz,0:Mvol)) ! interface geometry in real space; poloidal derivative; ! 18 Jul 14;
  RALLOCATE(tZij,(1:Ntz,0:Mvol))

  RALLOCATE(Rij,(1:Ntz,0:3,0:3)) ! these are used for inverse fft to reconstruct real space geometry from interpolated Fourier harmonics;
  RALLOCATE(Zij,(1:Ntz,0:3,0:3))
  RALLOCATE(sg,(1:Ntz,0:3))
  RALLOCATE(guvij,(1:Ntz,0:3,0:3,0:3)) ! need this on higher resolution grid for accurate Fourier decomposition;
  RALLOCATE(gvuij,(1:Ntz,0:3,0:3,0:3)) ! need this on higher resolution grid for accurate Fourier decomposition; ! workspace; 13 Sep 13;

  RALLOCATE(guvmne,(1:mne,0:6)) ! workspace for Fourier decomposition of metric terms; see me00ab; WASTEFUL MEMORY;
  RALLOCATE(guvmno,(1:mne,0:6)) ! workspace for Fourier decomposition of metric terms; see me00ab; WASTEFUL MEMORY;

  RALLOCATE(trigm,(1:2*Nt)) ! trignometric factors required for fast Fourier transform;
  RALLOCATE(trign,(1:2*Nz))
  RALLOCATE(trigwk,(1:2*Ntz))

  RALLOCATE(ijreal,(1:Ntz)) ! real space grid;
  RALLOCATE(ijimag,(1:Ntz))

  RALLOCATE(jireal,(1:Ntz))
  RALLOCATE(jiimag,(1:Ntz))

  isr = 'I' ; ifail = 0

  WCALL(al00aa, C06FUF, ( Nt, Nz, ijreal(1:Ntz), ijimag(1:Ntz), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz), ifail ) )

  isr = 'S' ! prepare FFTs; 23 Jan 13;

  FATALMESS(al00aa, ifail.ne.0, error constructing Fourier transform)

  RALLOCATE(efmn,(1:mne)) ! Fourier harmonics workspace; 24 Apr 13;
  RALLOCATE(ofmn,(1:mne))
  RALLOCATE(cfmn,(1:mne))
  RALLOCATE(sfmn,(1:mne))

  RALLOCATE(evmn,(1:mne)) ! Fourier harmonics workspace; 24 Apr 13;
  RALLOCATE(odmn,(1:mne))
  RALLOCATE(comn,(1:mne))
  RALLOCATE(simn,(1:mne))

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RALLOCATE(Bsupumn,(1:Nvol,0:1,1:mn)) ! Fourier components of {\bf B}\cdot\nabla \theta on boundary; required for virtual casing;
  RALLOCATE(Bsupvmn,(1:Nvol,0:1,1:mn)) ! Fourier components of {\bf B}\cdot\nabla \zeta  on boundary;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! Transformation to straight field line angle;

  RALLOCATE(diota,(0:1,-1:2,1:Mvol)) ! measured rotational transform on inner/outer interfaces in each annulus;

  RALLOCATE(glambda,(1:Ntz+1,0:2,0:1,1:Mvol)) ! save initial guesses for iterative calculation of rotational-transform; 21 Apr 13;
  
!  do vvol = 1, Mvol
!   
!   if( Igeometry.eq.1 .or.  vvol.gt.1 ) Lcoordinatesingularity = .false. ! either Cartesian geometry or annular volume; used elsewhere;
!   if( Igeometry.gt.1 .and. vvol.eq.1 ) Lcoordinatesingularity = .true.
!   
!   if( vvol.le.Nvol ) Lplasmaregion = .true.
!   if( vvol.gt.Nvol ) Lplasmaregion = .false.
!   
!   Lvacuumregion = .not. Lplasmaregion
!   
!   do innout = 0, 1 ! loop over inner and outer interface of each volume; 21 Apr 13;
!    
!    if( Lcoordinatesingularity .and. innout.eq.0 ) cycle ! transform on coordinates   axis     is not required              ; 20 Apr 13;
!    if( Lvacuumregion          .and. innout.eq.1 ) cycle ! transform on computational boundary is not required (not defined); 20 Apr 13;
!    
!    glambda(Ntz+1,0,innout,vvol) = iota(vvol-1+innout) ! global lambda: real-space transformation to straight-field-line angle; 21 Apr 13;
!   !glambda(Ntz+1,?,innout,vvol) = should not be too hard to estimate change in transform wrt change in helicity/poloidal flux or plasma current;
!
!   enddo ! end of do innout; 21 Apr 13;
!
!  enddo ! end of do vvol; 21 Apr 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! Construction of `force';

  RALLOCATE( Bemn,(1:mn,1:Mvol,0:1))
  RALLOCATE( Bomn,(1:mn,1:Mvol,0:1))
  RALLOCATE( Iomn,(1:mn,1:Mvol    ))
  RALLOCATE( Iemn,(1:mn,1:Mvol    ))
  RALLOCATE( Somn,(1:mn,1:Mvol,0:1))
  RALLOCATE( Semn,(1:mn,1:Mvol,0:1))
  RALLOCATE( Pomn,(1:mn,1:Mvol,0:2))
  RALLOCATE( Pemn,(1:mn,1:Mvol,0:2))

  RALLOCATE( BBe ,(1:Mvol-1))
  RALLOCATE( IIo ,(1:Mvol-1))
  RALLOCATE( BBo ,(1:Mvol-1))
  RALLOCATE( IIe ,(1:Mvol-1))

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! Parallel construction of derivative matrix;
  
  if( Lminimize.gt.0 ) then
   RALLOCATE( dBBdRZ,(1:Mvol,0:1,1:lgeometricaldof)) ! derivatives of            energy   integrals   wrt geometric Fourier harmonics;
   RALLOCATE( dIIdRZ,(1:Mvol    ,1:lgeometricaldof)) ! derivatives of (weighted) spectral constraints wrt geometric Fourier harmonics;
  endif

! RALLOCATE(   DenergyDrz,( 1:Nvol, 0:1, 1:lgeometricaldof                         )) ! derivatives of BB                   wrt geometry;

! RALLOCATE( DanalyticDrz,( 1:Nvol, 0:1, 1:lgeometricaldof                         )) ! derivatives of BB                   wrt geometry;
! RALLOCATE( DanalyDrzDrz,( 1:Nvol, 0:1, 1:lgeometricaldof, 0:1, 1:lgeometricaldof )) ! derivatives of BB                   wrt geometry; wrt geometry;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
  RALLOCATE( Bsubtemn,(1:mn,0:1,1:Mvol))
  RALLOCATE( Bsubzemn,(1:mn,0:1,1:Mvol))
  RALLOCATE( Bsubtomn,(1:mn,0:1,1:Mvol))
  RALLOCATE( Bsubzomn,(1:mn,0:1,1:Mvol))

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! Trigonometric factors;
   
! FATALMESS(al00aa,Ntz.le.0,have not yet computed Ntz)

! if( Igeometry.eq.2 ) then ! only required for cylindrical geometry;
   
!  RALLOCATE(coszeta,(1:Ntz))
!  RALLOCATE(sinzeta,(1:Ntz))
!  RALLOCATE(costeta,(1:Ntz))
!  RALLOCATE(sinteta,(1:Ntz))
   
!  do kk = 0, Nz-1 ; zeta = kk * pi2nfp / Nz
!   do jj = 0, Nt-1 ; teta = jj * pi2 / Nt ; jk = 1 + jj + kk*Nt
!    coszeta(jk) = cos(zeta) ; sinzeta(jk) = sin(zeta)
!    costeta(jk) = cos(teta) ! sinteta(jk) = sin(teta)
!   enddo
!  enddo
  
! endif ! end of if( Igeometry.eq.3 );

  RALLOCATE(cosi,(1:Ntz,1:mn))
  RALLOCATE(sini,(1:Ntz,1:mn))
  
  FATALMESS(al00aa, Nz.eq.0, illegal division)
  FATALMESS(al00aa, Nt.eq.0, illegal division)

  do ii = 1, mn ; mi = im(ii) ; ni = in(ii) ! loop over Fourier harmonics;
   
   do kk = 0, Nz-1 ; zeta = kk * pi2nfp / Nz ! loop over real space grid; toroidal;
    do jj = 0, Nt-1 ; teta = jj * pi2    / Nt ; jk = 1 + jj + kk*Nt ! loop over real space grid; poloidal;
     
     arg = mi * teta - ni * zeta ! shorthand;
     
     cosi(jk,ii) = cos(arg)
     sini(jk,ii) = sin(arg)
     
    enddo
   enddo
   
  enddo ! end of do ii; 13 May 13;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!  IALLOCATE(ijk,(1:mn,1:mn))
!  IALLOCATE(ikj,(1:mn,1:mn))
!
!  IALLOCATE(bmask,(1:mn))
!  
!  do jj = 1, mn ; mj = im(jj) ; nj = in(jj)
!   do kk = 1, mn ; mk = im(kk) ; nk = in(kk)
!    
!    do ii = 1, mn ; mi = im(ii) ; ni = in(ii)     
!     if( mk+mj.eq.mi .and. nk+nj.eq.ni ) ijk(kk,jj) =   ii
!     if( mk-mj.eq.mi .and. nk-nj.eq.ni ) ikj(kk,jj) =   ii
!     if( mj-mk.eq.mi .and. nj-nk.eq.ni ) ikj(kk,jj) = - ii
!    enddo ! end of do ii; 11 Aug 14;
!    
!   enddo ! end of do kk; 11 Aug 14;
!   
!   if( mj.eq.0 ) bmask(jj) = 1
!
!  enddo ! end of do jj; 11 Aug 14;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! Volume integrals;
  
  RALLOCATE(vvolume    ,(1:Mvol)) ! volume integral of \sqrt g;
  RALLOCATE(lBBintegral,(1:Mvol)) ! volume integral of B.B    ;
  RALLOCATE(lABintegral,(1:Mvol)) ! volume integral of A.B    ;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! Periodic orbits;
  
! Dzeta = pi / Npl ! 10 Oct 12; ! DEBUGGING; why is there a minus sign; under construction; Lagrangian integration; Lagrangian integration to be deleted; 11 Aug 14;
  
 !write(ounit,'("al00aa : ", 10x ," : myid=",i3," ; maxval(npq)=",i3," ;")') myid, maxval(npq)

 !if( maxval(npq).gt.0 ) then

   ALLOCATE(pqorbit,(1:Mvol,1:max(1,maxval(npq)))) ! periodic orbits; this is a structure; sub-arrays etc. are not initialized;  11 Aug 14;
  
   pqorbit(1:Mvol,1:max(1,maxval(npq)))%ok = 0     ! set flag that indicates if periodic orbit has been constructed;
   
   do vvol = 1, Mvol ! prepare workspace for constructing periodic orbits;
    
#ifdef DEBUG
    if( Wal00aa ) then
     write(ounit,'("al00aa : ", 10x ," : myid=",i3," ; vvol=",i3," ; npq=",i3," ;")') myid, vvol, npq(vvol)
    endif
#endif
    
    if( npq(vvol).eq.0 ) cycle

!   if( myid.ne.modulo(vvol-1,ncpu) ) cycle ! 10 Oct 12; this cycle statement can probably be reinserted -- periodic orbits are no longer broadcast; 
    
    FATALMESS(al00aa, q1(vvol).eq.0, denominator q1 cannot be 0)
    
    if( npq(vvol).ge.1 ) then ; pqorbit(vvol,1)%pq(1:2) = (/ p1(vvol), q1(vvol) * Nfp /) ! this is given on input; 27 Nov 12;
    endif
    
    FATALMESS(al00aa, npq(vvol).gt.1 .and. q2(vvol).eq.0,denominator q2 cannot be 0)
    
    if( npq(vvol).ge.2 ) then ; pqorbit(vvol,2)%pq(1:2) = (/ p2(vvol), q2(vvol) * Nfp /) ! this is given on input; 27 Nov 12;
    endif
    
    FATALMESS(al00aa, npq(vvol).gt.2 .and. abs( p1(vvol)*q2(vvol)-q1(vvol)*p2(vvol) ).ne.1, not neighboring : cannot construct noble irrational)
    
    if( npq(vvol).ge.3 ) then
     do nc = 3, npq(vvol) ; pqorbit(vvol,nc)%pq(1:2) = pqorbit(vvol,nc-2)%pq(1:2) + pqorbit(vvol,nc-1)%pq(1:2) ! noble convergents;
     enddo
    endif
    
    if( myid.ne.modulo(vvol-1,ncpu) ) cycle ! computation is distributed; 14 Nov 12;
    
    if( npq(vvol).gt.0 ) then
     
     do nc = 1, npq(vvol) ; pp = pqorbit(vvol,nc)%pq(1) ; qq = pqorbit(vvol,nc)%pq(2)
     
      ;                        pqorbit(vvol,nc)%to = zero ! default initial guess;
      ;                        pqorbit(vvol,nc)%so = zero
      
      if( nc.eq.1 ) then
       ;                       pqorbit(vvol,nc)%to = pqt(vvol)
       if( pqs(vvol).gt.zero ) pqorbit(vvol,nc)%so = pqs(vvol)
      endif

      FATALMESS(al00aa, qq.lt.0, illegal allocation)
     
      RALLOCATE(pqorbit(vvol,nc)%s,(0:qq)) ! periodic orbit constructed from field line following, pq00aa;
      RALLOCATE(pqorbit(vvol,nc)%t,(0:qq)) ! periodic orbit constructed from field line following, pq00aa;
      RALLOCATE(pqorbit(vvol,nc)%a,(0:qq)) ! periodic orbit constructed from field line following, pq00aa;
    
      FATALMESS(al00aa, qq*Npl.lt.1, illegal allocation)
 
      RALLOCATE(pqorbit(vvol,nc)%tt,(0:2*qq*Npl)) ! periodic orbit constructed from Lagrangian integration;
      RALLOCATE(pqorbit(vvol,nc)%ss,(1:2*qq*Npl)) ! periodic orbit constructed from Lagrangian integration;
      RALLOCATE(pqorbit(vvol,nc)%dt,(1:2*qq*Npl)) ! periodic orbit constructed from Lagrangian integration;
      
#ifdef DEBUG
      if( Wal00aa ) then
       write(ounit,'("al00aa : ", 10x ," : myid=",i3," ; qq="i4" ; Npl=",i3," ;")') myid, qq, Npl
      endif
#endif
      
      pqorbit(vvol,nc)%tt(0:2*qq*Npl) = (/ ( ii, ii = 0, 2*qq*Npl ) /) * pp * pi2 / ( 2*qq*Npl )
      pqorbit(vvol,nc)%ss(1:2*qq*Npl) = pqorbit(vvol,nc)%so
      
      FATALMESS(al00aa, abs(pqorbit(vvol,nc)%tt(2*qq*Npl)-pqorbit(vvol,nc)%tt(0)-pi2*pp).gt.vsmall, incorrectly enforcing periodicity constraint)
      
      if( Wal00aa ) then
       cput = GETTIME
       write(ounit, 1002) cput-cpus, myid, vvol, nc, pqorbit(vvol,nc)%pq(1:2), pqorbit(vvol,nc)%to, pqorbit(vvol,nc)%so
      endif
      
1002  format("al00aa : ",f10.2," : myid=",i3," : vvol=",i3," ; nc=",i3," ; periodic orbit : ("i5" ,"i6" ) ; to="f15.10" ; so="f15.10" ;")
      
     enddo ! end of do nc = 1, npq(vvol) ;
    
    endif ! end of if( npq(vvol).gt.0 ) ; 11 Aug 14;

   enddo ! end of do vvol = 1, Mvol
  
 !endif ! end of if( maxval(nfp).gt.0 ) ; 11 Aug 14;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! Quadratic-flux minimizing surfaces;  

  mnpq = max(1,maxval(npq)) ! perhaps this should be maxval(npq(1:Mvol)) ?;

  ALLOCATE(qfms,(1:Mvol,1:mnpq))

  do vvol = 1, Mvol
   do ipq = 1, mnpq
    qfms(vvol,ipq)%i = 0 ; qfms(vvol,ipq)%pq(1:2) = (/ 0, 0 /) ! intialize; these are given values in pq02aa;
   enddo
  enddo

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( Lwrpj ) then ! allocate array used to hold surface potentials on inner/outer of each interface;
   RALLOCATE(Itor,(1:Mvol,0:1))
   RALLOCATE(Gpol,(1:Mvol,0:1))
   RALLOCATE(spmn,(1:Mvol,0:1,1:mn))
  endif

  
  if( YESstellsym ) lmns = 1 + (mns-1)           ! number of independent degrees of freedom in angle transformation; 30 Jan 13; 
  if( NOTstellsym ) lmns = 1 + (mns-1) + (mns-1) ! number of independent degrees of freedom in angle transformation; 30 Jan 13; 
  
! RALLOCATE(dlambda,(1:lmns,-1:2,0:1,1:Mvol)) ! NOWHERE USED; 15 Sep 15;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(al00aa)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine al00aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
