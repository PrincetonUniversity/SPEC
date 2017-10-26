!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Constructs quadratic-flux minimizing surfaces using fieldline tracing.

!latex \item The routine \verb+pq02aa+ is called from \verb+spec+, and only if \verb+Nghd.gt.0+.
!latex       The \verb+qfmin+ structure array is defined in \verb+global+ and allocated (and assigned default values) in \verb+al0aa+.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine pq02aa( lvol ) 
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, one, two, pi, pi2, goldenmean
  use numerical, only :
  use fileunits, only : ounit, qunit
  use inputlist, only : Wpq02aa, ext, Nfp, pscale, Ni, npq, Nghd, Mpqits, odetol, iota, oita
  use cputiming, only : Tpq02aa
  use allglobal, only : myid, cpus, pi2nfp, ss, Nz, pqorbit, pqorbit, qfms, jgd, kgd, actiongradient, ivol
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in) :: lvol
  
  INTEGER            :: pp, qq, nc, its, nqq, id02bjf, Node, itgt, npmodq, isfla, ii, jgdmin, ifail, fNghd
  REAL               :: low, upp, lu(2), ti, ri, acti, dr, dv, tm(2, 2), st, tt, rr, err, tol, wk(7*20), trd(7), rto(7), acto, rtd(7), rtm(7), xt, xe, fdiff=1.0e-04
  REAL               :: det, iotal, iotau, liota, ts, sflerr, xx, shear
  REAL,  allocatable :: nu(:)
  CHARACTER          :: relabs, suf*4
  EXTERNAL           :: pfield, bf00aa_out, bf00aa_end
  
  BEGIN(pq02aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  ivol = lvol

  Node = 7 ; relabs = 'D' ; tol = odetol ; fNghd = 4*Nghd
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  low = ss(lvol)%s(0) ; upp = ss(lvol)%s(Ni(lvol)) ; lu = (/ low, upp /)

 !FATALMESS(pq02aa, .true., definition of iota has changed)

  FATALMESS(pq02aa, lvol.eq.1, need to supply lower bound for rational search)

  iotal = oita(lvol-1) ; iotau = iota(lvol)
 !iotal = iota(lvol-1) ; iotau = iota(lvol) ! NEED TO EXAMINE DEFINITION OF OITA; 29 Apr 15;

  shear = ( iotau-iotal ) / ( upp-low)
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  do nc = 1, npq(lvol)

   pp = pqorbit(lvol,nc)%pq(1)
   qq = pqorbit(lvol,nc)%pq(2)*Nfp
   
   liota = pp*Nfp * one / qq ! loop over convergent rationals;

   cput = GETTIME ; write(ounit,'("pq02aa : ",f10.2," : myid=",i3," ; lvol=",i3," ; constructing qfms ; ",i3," /"i4" ="f9.6" ;")') cput-cpus, myid, lvol, pp, qq, pp*one/qq

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
!latex \item If the selected rational is within the rotational transform range of the given annulus, the flag \verb+qfms(lvol,nc)%i=1+ is set. 
!latex This will be used in \verb+pp00aa+ when constructing the \Poincare plot.
!latex If the selected rational is not within the rotational transform range of the given annulus, control will cycle.

! the following line will encounter an error if npq(lvol).gt.1 ;

   qfms(lvol,nc)%pq(1:2)=(/ pp, qq /)
   
   ALLOCATE(qfms(lvol,nc)%trvs,(1:4,0:fNghd,0:qq)) ! full family of full periodic pseudo curves;
   qfms(lvol,nc)%trvs=zero
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   if( ( liota - iotal ) * ( liota - iotau ) .ge. zero ) then ; write(ounit,1001) myid, lvol, pp, qq, liota ; cycle
   endif

1001 format("pq02aa : "10x" : myid=",i3," ; lvol=",i3," ; ",i3," /"i4" ="f9.6" ; not in range ;")
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   qfms(lvol,nc)%i = 1 ; qfms(lvol,nc)%pq(1:2)=(/ pp, qq /)
   
   do npmodq = 1, qq-1
    if( mod( npmodq * pp, qq ).eq.1 ) exit ! periodicity;
   enddo
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item The qfmin-surface is constructed iteratively (at max. \verb+Mpqits+) in order to obtain periodic pseudo curves equally spaced in the straight-field line angle.

   do isfla = 0, Mpqits
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
    if( isfla.eq.0 ) jgdmin = 0
    if( isfla.gt.0 ) jgdmin = 0 ! no need to recompute? careful with shifting!!!

    do jgd = jgdmin,fNghd-1 ! for each periodic pseudo orbit; NOTE : should be jgd=tNgd from periodicity; poloidal index;
     
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
     
     if( isfla.eq.0 ) then ! no previous estimate of straight field line angle is available;
      
      if( jgd.eq.0 ) then ; ti = zero                                                     ; ri = pqorbit(lvol,nc)%so           ; acti = zero
      else                ; ti = jgd * mod( qfms(lvol,nc)%trvs(1,0,npmodq), pi2 ) / fNghd ; ri = qfms(lvol,nc)%trvs(2,jgd-1,0) ; acti = qfms(lvol,nc)%trvs(3,jgd-1,0)
      endif
      
     else ! can estimate location of straight field line angle;
      
      ts = jgd * ( pi2/qq ) / fNghd
      
      do ii = 0,fNghd-1
       
       if( ts.ge.qfms(lvol,nc)%trvs(4,ii,0) .and. ts.le.qfms(lvol,nc)%trvs(4,ii+1,0) ) then
        
        xx = ( ts - qfms(lvol,nc)%trvs(4,ii,0) ) / ( qfms(lvol,nc)%trvs(4,ii+1,0) - qfms(lvol,nc)%trvs(4,ii,0) )
        
        ti   = qfms(lvol,nc)%trvs(1,ii,0) + ( qfms(lvol,nc)%trvs(1,ii+1,0)-qfms(lvol,nc)%trvs(1,ii,0) ) * xx
        ri   = qfms(lvol,nc)%trvs(2,ii,0) + ( qfms(lvol,nc)%trvs(2,ii+1,0)-qfms(lvol,nc)%trvs(2,ii,0) ) * xx
        acti = qfms(lvol,nc)%trvs(3,ii,0) + ( qfms(lvol,nc)%trvs(3,ii+1,0)-qfms(lvol,nc)%trvs(3,ii,0) ) * xx
        
        exit
        
       endif

      enddo
      
     endif
     
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Each periodic pseudo orbit is located by integrating the pseudo-field the full periodicity and using an iterative method.

     its = 0 ; dr = zero ; dv = zero
     
     do ; its = its+1 ! iterate to locate periodic-pseudo orbit;
      
      ri = ri + dr ; acti = acti + dv ! update guess;
      
      rto = (/ ri, ti, zero, zero, zero, one, zero /) ; acto = acti
      
      kgd = 0 ; qfms(lvol,nc)%trvs(1:3,jgd,kgd) = (/ ti, ri, acti /)
      
!latex \item The tangent mapping is constructed by finite-differences (but a better approach could construct the tangent map by integrating the differentiated equations).

      do itgt = 0, 2 ! compute tangent matrix by finite-differences;
       
       rtd(1:7) = rto(1:7) ; actiongradient = acto
       
       if( itgt.eq.1 ) rtd(1)         = rto(1) + fdiff
       if( itgt.eq.2 ) actiongradient =  acto  + fdiff

       do nqq = 1, qq ! full q-period mapping; note that for the symmetric periodic pseudo orbit = periodic true orbit, only need to follow half periodicity distance;
        
        xt = zero ; xe = xt + pi2nfp ; id02bjf = 0 ; call D02BJF( xt, xe, Node, rtd, pfield, tol, relabs, bf00aa_out, bf00aa_end, wk, id02bjf )
        
!latex \item The periodic pseudo orbits are presently only saved on one toroidal plane -- this needs to be developed if the qmfin surface is to be Fourier decomposed.
        
        if( itgt.eq.0 ) then ; kgd = kgd+1 ; qfms(lvol,nc)%trvs(1:3,jgd,kgd) = (/ rtd(2), rtd(1), acti /)
        endif
        
       enddo

       if( itgt.eq.0 ) qfms(lvol,nc)%trvs(4,jgd,0) = rtd(7) / ( pi2 * qq ) - pp * pi ! straight field line angle;
       
       if( itgt.eq.0 ) rtm(1:7) = rtd(1:7) ! original mapped;
       
       if( itgt.gt.0 ) tm(1:2,itgt) = ( rtd(1:2) - rtm(1:2) ) / fdiff ! tangent map;
       
      enddo ! end of do itgt=0,2;
      
      err = sqrt( (rto(1)-rtm(1))**2 + (rto(2)+pi2*pp-rtm(2))**2 )
      
      tm(1,1) = tm(1,1) - one ! periodicity requirement;
      
      det = ( tm(1,1)*tm(2,2) - tm(1,2)*tm(2,1) )
      
      dr = (  tm(2,2) * (rto(1)-rtm(1)) - tm(1,2) * (rto(2)+pi2*pp-rtm(2)) ) / det ! determine correction;
      dv = ( -tm(2,1) * (rto(1)-rtm(1)) + tm(1,1) * (rto(2)+pi2*pp-rtm(2)) ) / det ! determine correction;
      
!latex \item Periodic pseudo orbits are located to a tolerance of \verb+qq*odetol+, or the iterations are terminated if \verb+its.gt.Mpqits+.

      if( err.lt. qq*odetol .or. its.gt.Mpqits ) exit ! converged, or failed;
      
     enddo ! end of do ; its = its+1 loop;
     
    !cput = GETTIME ; write(ounit,1002)cput-cpus,myid,pp,qq,ti,its,rto(1),acto,err
     
    enddo ! end of do jgd=0,tNgd-1 loop;

    qfms(lvol,nc)%trvs(1,fNghd,0) = mod( qfms(lvol,nc)%trvs(1,0,npmodq), pi2 )  ! add periodic curve; theta;
    qfms(lvol,nc)%trvs(2,fNghd,0) = qfms(lvol,nc)%trvs(2,0,npmodq)              ! add periodic curve; radial;
    qfms(lvol,nc)%trvs(3,fNghd,0) = qfms(lvol,nc)%trvs(3,0,0)                   ! add periodic curve; action-gradient;
    qfms(lvol,nc)%trvs(4,fNghd,0) = qfms(lvol,nc)%trvs(4,0,0) + pi2/qq          ! add periodic curve; theta_s;

    qfms(lvol,nc)%trvs(4,0:fNghd,0) = qfms(lvol,nc)%trvs(4,0:fNghd,0)-qfms(lvol,nc)%trvs(4,0,0) ! recenter;
    
    sflerr = sum( abs( qfms(lvol,nc)%trvs(4,0:fNghd-1,0) - (/(jgd,jgd=0,fNghd-1)/)*pi2/qq/fNghd ) ) / fNghd

!latex \item The repetitive construction of qfmin surfaces to obtain periodic pseudo curves equally spaced in the straight-field-line angle is terminated 
!latex when a measure of the straight field line angle error is less than $10^{-3}$.

    if( sflerr.lt.1.0e-03 ) exit ! exit isfla iterations;

   enddo ! end of do isfla=0,1 loop;
   
  enddo ! end of do nc=1,npq(lvol) loop over convergent rationals;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item The periodic pseudo curves are saved to files, \verb+.ext.qfms:xxxx+ where \verb+xxxx+ is an annulus label.

  cput = GETTIME

  write(suf,'(i4.4)')lvol
  open(qunit+myid,file="."//trim(ext)//".qfms."//suf,status="unknown",form="unformatted") ! this will contain Poincare plot data; 
  write(qunit+myid)Nghd,npq(lvol)

  ALLOCATE(nu,(0:fNghd-1)) 

  do nc = 1, npq(lvol) ; pp = pqorbit(lvol,nc)%pq(1) ; qq = pqorbit(lvol,nc)%pq(2)

   nu(0:fNghd-1) = qfms(lvol,nc)%trvs(3,0:fNghd-1,0)
   ifail=0 ; call slowft( nu(0:fNghd-1), fNghd, ifail )

   write(ounit,'("pq02aa : ",f10.2," : myid=",i3," ; lvol=",i3," ; pscale="es13.5" ; ",i3," /"i4" : nu_m=",99es12.4)')cput-cpus,myid,lvol,pscale,pp,qq,nu(fNghd-1:fNghd-Nghd:-1)

   write(qunit+myid)pp,qq,shear
   write(qunit+myid)qfms(lvol,nc)%trvs(1,0:fNghd,0:qq-1)
   write(qunit+myid)qfms(lvol,nc)%trvs(2,0:fNghd,0:qq-1)
   write(qunit+myid)qfms(lvol,nc)%trvs(3,0:fNghd,0:qq-1)
  enddo

  DEALLOCATE(nu)

  close(qunit)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(pq02aa)

1002 format("pq02aa : ",f10.2," : myid=",i3," ; (",i3,","i4" ) : t="f6.3" : its=",i3," ; (r,v)=("f13.10,es10.2" ) ; err="es7.0" ;")
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine pq02aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine pfield( phi, rtd, Brtd ) ! ode's for pseudo field line and tangent pseudo map;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero,one
  use inputlist, only :
  use allglobal, only : myid,cpus,actiongradient
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  REAL, intent(in)  :: phi, rtd(7)
  REAL, intent(out) :: Brtd(7)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  call bf00aa( phi, rtd(1:6), Brtd(1:6) )

  Brtd(1) = Brtd(1) + actiongradient ! convert real magnetic field to pseudo magnetic field;
  
  Brtd(7) = rtd(2) ! integrate angle to get area; used for constructing straight pseudo field line angle;

! will probably need to restructure bf00aa to allow pseudo tangent to be constructed by integration;

 ! Btrd(3:6)=(/ Btrd(3) * trd(3) + Btrd(4) * trd(5), Btrd(3) * trd(4) + Btrd(4) * trd(6), Btrd(5) * trd(3) + Btrd(6) * trd(5) + one, Btrd(5) * trd(4) + Btrd(6) * trd(6)/)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  return
  
end subroutine pfield

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
