!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Constructs \Poincare plot and `approximate' rotational-transform.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \end{enumerate} \subsection{relevant input variables} \begin{enumerate}
  
!latex \item The resolution of \Poincare plot is controlled by 
!latex       \begin{itemize} 
!latex       \item \verb+nPtraj+ trajectories will be located in each volume;
!latex       \item \verb+nPpts+ iterations per trajectory;
!latex       \item \verb+odetol+ o.d.e. integration tolerance;
!latex       \end{itemize}
!latex       The magnetic field (and tangent field) is given by \verb+bf00aa+.
  
!latex \item The approximate rotational transform is determined (in \verb+pq00aa+) by field line integration.
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex  \end{enumerate} \subsection{format of output: \Poincare} \begin{enumerate}
!latex \item The \Poincare data is written to \verb+.ext.poincare:xxxx+, where \verb+xxxx+ is an integer indicating the volume.
!latex       The format of this file is as follows:
!latex \begin{verbatim}
!latex  write(svol,'(i4.4)')lvol ! lvol labels volume;
!latex  open(lunit+myid,file="."//trim(ext)//".poincare."//svol,status="unknown",form="unformatted")
!latex  do until end of file
!latex   write(lunit+myid) Nz, nPpts                                   ! integers
!latex   write(lunit+myid) poincaredata(1:4,0:Nz-1,1:nPpts)            ! doubles
!latex  enddo
!latex  close(lunit+myid)
!latex \end{verbatim}
!latex       where \begin{itemize}
!latex       \item $\t \equiv$ \verb+poincaredata(1,k,j)+ is the poloidal angle,
!latex       \item $\s \equiv$ \verb+poincaredata(2,k,j)+ is the radial coordinate, 
!latex       \item $ R \equiv$ \verb+poincaredata(3,k,j)+ is the cylindrical $R$, 
!latex       \item $ Z \equiv$ \verb+poincaredata(4,k,j)+ is the cylindrical $Z$,
!latex       \end{itemize}
!latex \item The integer \verb+k=0,Nz-1+ labels toroidal planes, so that $\phi = ( 2 \pi / $\verb+Nfp+$ ) ( k / $\verb+Nz+$)$,
!latex \item The integer \verb+j=1,nPpts+ labels toroidal iterations.
!latex \item Usually (if no field line integration errors are encountered) the number of field lines followed in volume \verb+lvol+
!latex       is given by $N+1$, where the radial resolution, $N\equiv$\verb+Ni(lvol)+, is given on input.
!latex       This will be over-ruled by if \verb+nPtrj(lvol)+, given on input, is non-negative.
!latex       The starting location for the field line integrations are equally spaced in the radial coordinate $s_i=s_{l-1}+ i (s_{l}-s_{l-1})/N$ for $i=0,N$,
!latex       along the line $\t=0$, $\z=0$.
!latex \item Additional field lines inside islands may be followed if quadratic-flux minimizing surfaces have been constructed,
!latex       and for these field lines $i_{k,j}$ is negative.
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine pp00aa( lvol ) 
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, one, two, pi
  
  use numerical, only :
  
  use fileunits, only : ounit, lunit
  
  use inputlist, only : Wmacros, Wpp00aa, Nvol, Lrad, ext, odetol, nPpts, nPtrj, Lconstraint, iota, oita, pqs, pqt, npq
  
  use cputiming, only : Tpp00aa
  
  use allglobal, only : myid, ncpu, cpus, &
                        Nz, pi2nfp, &
                        Ltangent, ivol, Mvol, &
                        Lcoordinatesingularity, Lvacuumregion, &
                        Lfieldlinedirection, diota, &
                        pqorbit
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in)  :: lvol
  INTEGER              :: lnPtrj, ioff, itrj, ith, ld02bjf, nc, pp, qq
  INTEGER, allocatable :: id02bjf(:)
  REAL                 :: sti(1:2), dth, ltransform(1:2)
  REAL, allocatable    :: poincaredata(:,:,:), fittedtransform(:,:)
  CHARACTER            :: svol*4 ! perhaps this should be global; 11 Aug 13;
  
  BEGIN(pp00aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  FATALMESS(pp00aa, lvol.lt.1 .or. lvol.gt.Mvol, invalid volume )
  FATALMESS(pp00aa, nPpts.le.0, nothing to do )
#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  ivol = lvol ; Ltangent = 0 ! tangent map is not required for Poincare plot; 28 Jan 13;
  write(svol,'(i4.4)')lvol
  
  open(lunit+myid,file="."//trim(ext)//".poincare."//svol,status="unknown",form="unformatted")
  
1001 format("pp00aa : ",f10.2," : myid=",i3," ; lvol=",i3," ; odetol=",es8.1," ; nPpts=",i8," ; lnPtrj=",i3," ;")

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  Lfieldlinedirection = +1 ! hereafter follow field lines forwards;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RALLOCATE(poincaredata,(1:4,0:Nz-1,1:nPpts)) ! for block writing to file (allows faster reading of output data files for post-processing plotting routines);
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  do nc = 1, npq(lvol) ! loop over convergents in each volume ; 11 Aug 13;
   
   pp = pqorbit(lvol,nc)%pq(1)
   qq = pqorbit(lvol,nc)%pq(2)! * Nfp
   
   if( pqorbit(lvol,nc)%ok.ne.1 ) exit
   
!if( Wpp00aa ) then
! cput = GETTIME
! write(ounit,'("pp00aa : ",f10.2," : myid=",i3," ; (",i3," ,"i4" ) : (s,t)=("f10.6" ,"f10.6" ) ;")') &
!  cput-cpus, myid, pp, qq, pqorbit(lvol,nc)%to, pqorbit(lvol,nc)%so
!endif
   
   sti(1:2) = (/ pqorbit(lvol,nc)%so, pqorbit(lvol,nc)%to /)
   
   CALL(pp00aa, pp00ab,( lvol, sti(1:2), Nz, nPpts, poincaredata(1:4,0:Nz-1,1:nPpts), ltransform(1:2), ld02bjf ) )
   
   if( Wpp00aa ) then
    cput = GETTIME
    write(ounit,1002) cput-cpus, myid, lvol, nc, sti(1:2), ld02bjf, ltransform(2)
   endif

   if( ld02bjf.eq.0 ) then ! will only write successfully followed trajectories to file; 28 Jan 13;
    write(lunit+myid) Nz, nPpts
    write(lunit+myid) poincaredata(1:4,0:Nz-1,1:nPpts)
   endif
   
  enddo ! end of do nc; 11 Aug 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
1002 format("pp00aa : ",f10.2," : myid=",i3," ; lvol=",i3," ; ",i3," : (s,t)=(",f21.17," ,",f21.17," ) ;":" id02bjf=",i3," ; transform=",es23.15,&
  " ;":" error=",es13.5," ;")
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( nPtrj(lvol).ge.0 ) lnPtrj =    nPtrj(lvol) ! selected Poincare resolution;
  if( nPtrj(lvol).lt.0 ) lnPtrj = 2 * Lrad(lvol) ! adapted  Poincare resolution;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( lnPtrj.gt.0 ) then
   
   ;                            ioff = 0 ! default is to start from inner interface ; 11 Aug 13;
   if( Lcoordinatesingularity ) ioff = 1 ! keep away from coordinate axis;
   
   IALLOCATE(id02bjf,(ioff:lnPtrj)) ! error flag that indicates if field lines successfully followed; 22 Apr 13;
   
   RALLOCATE(fittedtransform,(ioff:lnPtrj,1:2)) ! will always need fittedtransform(0,1:2);
   
   do itrj = ioff, lnPtrj ! initialize Poincare plot with trajectories regularly spaced between interfaces along \t=0;
    
    if( lvol.eq.1 ) then ; sti(1:2) = (/ - one + itrj * itrj * two / lnPtrj / lnPtrj, pqt(lvol) * pi /) ! start points along \t = given on input; 02 May 13;
    else                 ; sti(1:2) = (/ - one + itrj        * two / lnPtrj         , pqt(lvol) * pi /) ! start points along \t = given on input; 02 May 13;
    endif
    
    if( itrj.eq.lnPtrj ) sti(1:2) = (/ one, pqt(lvol) * pi /) ! avoid starting outside of domain 21 Jan 16;
    
    CALL(pp00aa, pp00ab,( lvol, sti(1:2), Nz, nPpts, poincaredata(1:4,0:Nz-1,1:nPpts), fittedtransform(itrj,1:2), id02bjf(itrj) ) )
    
    if( Wpp00aa ) then
     cput = GETTIME
     if( Lconstraint.eq.1 ) then
      if( itrj.eq.0                      ) &
   write(ounit,1002) cput-cpus, myid, lvol, itrj, sti(1:2), id02bjf(itrj), fittedtransform(itrj,2), fittedtransform(itrj,2)-oita(lvol-1)
      if( itrj.gt.0 .and. itrj.lt.lnPtrj ) &
   write(ounit,1002) cput-cpus, myid, lvol, itrj, sti(1:2), id02bjf(itrj), fittedtransform(itrj,2)
      if(                 itrj.eq.lnPtrj ) &
   write(ounit,1002) cput-cpus, myid, lvol, itrj, sti(1:2), id02bjf(itrj), fittedtransform(itrj,2), fittedtransform(itrj,2)-iota(lvol  )
     else
   write(ounit,1002) cput-cpus, myid, lvol, itrj, sti(1:2), id02bjf(itrj), fittedtransform(itrj,2)
     endif
    endif
    
    if( id02bjf(itrj).eq.0 ) then ! will only write successfully followed trajectories to file; 28 Jan 13;
     write(lunit+myid) Nz, nPpts
     write(lunit+myid) poincaredata(1:4,0:Nz-1,1:nPpts)
    endif
    
   enddo ! end of do itrj; 25 Jan 13;

  endif ! end of if( lnPtrj.gt.0 ) ; 11 Aug 14;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
  DEALLOCATE(poincaredata)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  close(lunit+myid) ! Poincare plot is finished; will re-use lunit to write rotational transform data;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \end{enumerate} \subsection{format of output: rotational-transform} \begin{enumerate}
  
!latex \item The rotational-transform data is written to \verb+.exttransform:xxxx+, where \verb+xxxx+ is an integer indicating the volume.
!latex       The format of this file is as follows:
!latex \begin{verbatim}
!latex  open(lunit+myid,file="."//trim(ext)//".transform."//svol,status="unknown",form="unformatted") ! rotational transform data;
!latex  write(lunit+myid) lnPtrj-ioff+1                                      ! integer
!latex  write(lunit+myid) diota(0:1,0,lvol)                                  ! doubles
!latex  write(lunit+myid) ( fittedtransform(itrj,1:2), itrj = ioff, lnPtrj ) ! doubles
!latex  close(lunit+myid)
!latex \end{verbatim}
  
  if( lnPtrj.gt.0 ) then

   open(lunit+myid,file="."//trim(ext)//".transform."//svol,status="unknown",form="unformatted") ! rotational transform data;
   write(lunit+myid) lnPtrj-ioff+1
   write(lunit+myid) diota(0:1,0,lvol)
   write(lunit+myid) ( fittedtransform(itrj,1:2), itrj = ioff, lnPtrj )
   close(lunit+myid)
  
   DEALLOCATE(id02bjf)
   
   DEALLOCATE(fittedtransform)
   
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(pp00aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine pp00aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!l!tex \end{enumerate} \subsection{construction of unstable manifold} \begin{enumerate}

!l!tex \item The unstable manifold is illustrated by mapping a short line segment emanating from the unstable periodic orbit in the unstable direction.

!  if( npq(lvol).gt.0 .and. Munstable.gt.0 ) then ! could also include a check on Mpqits;
    
!   RALLOCATE(poincaredata,(1:4,0:Nz-1,1:Nunstable,1:Munstable)) ! for block writing to file (allows faster reading of output data files for post-processing plotting routines);
    
!   do nc = 1, npq(lvol) ! loop over convergents;
     
!    cput = GETTIME
     
!    pp = pqorbit(lvol,nc)%pq(1) ; qq = pqorbit(lvol,nc)%pq(2) ! identify periodicity; shorthand;
     
!    if( Wpp00aa ) then 
!     write(ounit,1000)cput-cpus, myid, lvol, pp, qq, pqorbit(lvol,nc)%ok
!    endif
     
!    if( pqorbit(lvol,nc)%ok.eq.1 ) then ! periodic orbit has been successfully located;
      
!     wr(1:2) = pqorbit(lvol,nc)%wr(1:2) ; wi(1:2) = pqorbit(lvol,nc)%wi(1:2) ! eigenvalues: real and imaginary; shorthand;
      
!1000 format("pp00aa : ":,,f10.2," : myid=",i3," ; lvol=",i3," ; (",i3," ,",i3," ) ; ok="i2" ;":" (s,t)=("es23.15" ,"es23.15" ) ; residue="es13.5" ; eval="2(f9.5" +"f9.5" i ,"))
      
!     if( pqorbit(lvol,nc)%residue.lt.zero ) then ! periodic orbit is unstable; eigenvalue/eigenvector should be real;
!      
!     !write(ounit,1000) ! 14 Nov 12;
!      write(ounit,1000)cput-cpus, myid, lvol, pp, qq, pqorbit(lvol,nc)%ok, pqorbit(lvol,nc)%so, pqorbit(lvol,nc)%to, pqorbit(lvol,nc)%residue, ( wr(iev), wi(iev), iev=1,2 )
!      
!      do iev = 1, 2 ! loop over stable and unstable directions;
!       
!       if( pqorbit(lvol,nc)%wr(iev).gt.one ) Lfieldlinedirection = +1 ! this direction is un-stable; follow field lines  forwards; Lfieldlinedirection is used elsewhere;
!       if( pqorbit(lvol,nc)%wr(iev).lt.one ) Lfieldlinedirection = -1 ! this direction is    stable; follow field lines backwards; Lfieldlinedirection is used elsewhere;
!
!       if( Lfieldlinedirection.lt.0 ) cycle    ! this direction is    stable; follow field lines backwards; Lfieldlinedirection is used elsewhere;
!       
!       do linnout = -1, 1, 2 ! follow in both directions along eigenvector;
!        
!        do iupdown = -1, 1, 2 ! follow stellarator symmetric pair;
!
!         if( Lpqsym.eq.1 .and. iupdown.ne.1 ) cycle ! Lpqsym=1 enforces stellarator symmetry, so \t = 0; ! 15 Oct 12;
!
!         do ii = 1, Nunstable 
!    
!          sti(1:2) = (/ pqorbit(lvol,nc)%so , iupdown * pqorbit(lvol,nc)%to /) + ii * linnout * pqorbit(lvol,nc)%vr(1:2,iev) * dunstable / Nunstable
!          
!          call pp00ab( lvol, sti(1:2), Nz, Munstable, poincaredata(1:4,0:Nz-1,ii,1:Munstable), fittedtransform(0,1:2), id02bjf ) ! 15 Oct 12;
!          
!         !poincaredata(5,0:Nz-1,ii,1:Munstable) = -one ! labels field line; e.g. can be used to color trajectories separately in external plotting diagnostic;
!          
!         enddo ! end of do ii;
!         
!         do ii = 1, Munstable ! loop over short loop;
!          write(lunit+myid) Nz, Nunstable
!          write(lunit+myid) poincaredata(1:4,0:Nz-1,1:Nunstable,ii) ! write field line trajectory to file; block data over long loop;
!         enddo
!         
!        enddo ! end of do iupdown;
!        
!       enddo ! end of do linnout;
!       
!      enddo ! end of do iev;
!      
!     endif ! end of if( pqorbit(lvol,nc)%residue.lt.zero );
      
!    endif ! end of if( pqorbit(lvol,nc)%ok.eq.1 );
     
!   enddo ! end of do nc;
    
!   DEALLOCATE(poincaredata)
    
!  endif ! end of if( npq(lvol).gt.0 .and. Munstable.gt.0 ) then ; 20 Jun 14;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! do ii = 1 + lnPtrj, - Nghd*nqfmsok, -1 ; iigh = ii + Nghd*nqfmsok ! begin loop over trajectories;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!  Node = 6 ; RA = 'D' ; tol = odetol
   
!  if( Lpoincare.eq.2 ) Node = Node + 2 ! two extra o.d.e.s in this case
   
!  if( ii.lt.0 ) then ; st(1:6) = (/ qfmsst(iigh,1), qfmsst(iigh,2), one, zero, zero, one /) ! first quadratic-flux minimizing surface;
!  endif
   
!  if( ii.ge.0 ) then ! equally spaced (in surface label, along \t=0 \z=0) starting point for field line integrations;


!   if( Lpqsym.eq.0 .and. pqs(lvol).gt.zero ) lss = interfacelabel(lvol-1) + ii * (            pqs(lvol) - interfacelabel(lvol-1) ) / lnPtrj

!                     st(1:8) = (/ lss, zero     , one, zero, zero, one, zero, zero /) ! default is to go along \t = 0 line;
!   if( Lpqsym.eq.0 ) st(1:8) = (/ lss, pqt(lvol), one, zero, zero, one, zero, zero /)

!  endif

!  if( ii.eq.1+lnPtrj ) then ! include periodic orbit
!   st(1:8) = (/ pqorbit(lvol,1)%so, pqorbit(lvol,1)%to, one, zero, zero, one, zero, zero /)
!   st(1:2) = st(1:2) + pqorbit(lvol,1)%umanifold(1:2) * small ! perturb away from periodic orbit in direction of unstable manifold;
!  endif

!  if( ii.eq.1+lnPtrj .and. pqorbit(lvol,1)%ok.eq.0 ) cycle

!  if( lvol.eq.1 .and. ii.ge.0 .and. ii.lt.ioff ) cycle ! dont follow field lines on/near coordinate axis;
   
!  sti(1:2) = st(1:2)
   
!  ppt(1:2) = (/ st(2), st(1) /) ! seems redundant;
   
!  leastsqfit(1:5) = (/ zero, zero, zero, ppt(1), one /) ! initialize summations for least squares fit;
   
!  lerror = 0 ! used to exit loop over iterations, do jj and do kk loops, if an integration error is encountered;

!  poincaredata(1:5,0:Nz-1,1:nPpts) = zero

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!  do jj = 1, nPpts ! loop over iterations;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
!   zst = zero ! starting Poincare plane;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!   do kk = 0, Nz-1 ! loop over toroidal Poincare cross sections;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
     
!    stz(1:3) = (/ ppt(2), mod(ppt(1),pi2), zst /) ! toroidal coordinates;

!    if( Lpoincare.eq.2 ) stz(3) = st(7) ! zeta is not our integration coordinate, but one of the ODE variables

!    call co00aa( lvol, stz(1:3), RpZ(1:3), dR(0:3), dZ(0:3), jacobian, guv(1:3,1:3) )

!    ppt(3:4)=(/ RpZ(1), RpZ(3) /) ! cylindrical coordinates;

! ii labels field line; if ii lt 0 then island chain (as calculated by quadratic-flux min surface);
!    poincaredata(1:5,kk,jj) = (/ mod(ppt(1),pi2), ppt(2), ppt(3), ppt(4), ii*one /) ! save Poincare information to array;
        
!    zend = zst + pi2nfp/Nz

!    if( Lpoincare.eq.2 ) then
!     if( zst .le. 0 ) then
!      zend = 1.0/small ! choose a large number so that the first zero crossing is certain to be found
!     else
!      zend = zst + st(8)* 10 ! chose a more sensible value to save integration time
!     endif
!     st(8) = 0 ! Reset the \zeta distance from our starting position back to zero
!    endif
    
!    id02bjf = 1
!    if( Lpoincare.eq.2 ) then ; CALL(D02BJF,( zst, zend, Node, st, bf00aa, tol, RA, D02BJX    , bf00aa_end, realwork, id02bjf )) ! integrate to next toroidal subdivision;
!    else                      ; CALL(D02BJF,( zst, zend, Node, st, bf00aa, tol, RA, D02BJX    , D02BJW    , realwork, id02bjf )) ! integrate to next toroidal subdivision;
!    endif

!    cput = GETTIME
!    select case(id02bjf)                                                      !         1         2         3         4         5         6
!    case(0) ; ! will only give screen output if an error is encountered;      !123456789012345678901234567890123456789012345678901234567890123
!    case(1) ; lerror=1 ; write(ounit,2001)cput-cpus,myid,lvol,ii,jj,kk,id02bjf,"input error                                                    "
!    case(2) ; lerror=2 ; write(ounit,2001)cput-cpus,myid,lvol,ii,jj,kk,id02bjf,"error integrating field                                        "
!    case(3) ; lerror=3 ; write(ounit,2001)cput-cpus,myid,lvol,ii,jj,kk,id02bjf,"tol is too small to take initial step                          "
!    case(4) ; lerror=4 ; write(ounit,2001)cput-cpus,myid,lvol,ii,jj,kk,id02bjf,"xsol not reset or xsol is behind x after initial call to output"
!    case(5) ; lerror=5 ; write(ounit,2001)cput-cpus,myid,lvol,ii,jj,kk,id02bjf,"xsol not reset or xsol is behind last xsol                     "
!    case(6) ; lerror=6 ; write(ounit,2001)cput-cpus,myid,lvol,ii,jj,kk,id02bjf,"termination function did not change sign                       "
!    case(7) ; lerror=7 ; write(ounit,2001)cput-cpus,myid,lvol,ii,jj,kk,id02bjf,"serious error                                                  "
!   !case default
!    end select
!2001 format("pp00aa : ",f10.2," : myid=",i3," ; lvol=",i3," ; (ii,jj,kk)=("i4" ,"i4" ,"i4" ); ifail="i2" ; "a63)
 
!    if( lerror.ne.0 ) exit ! an integration error was encountered; exit do kk loop;

!    ppt(1:2) = (/ st(2), st(1) /) ! again, seems redundant;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
     
!   enddo ! end of do kk;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
!   if( lerror.ne.0 ) exit ! an integration error was encountered; exit do jj loop;
    
!   leastsqfit(1:5) = leastsqfit(1:5) + (/ (jj*pi2nfp)**2, jj*pi2nfp, jj*pi2nfp*ppt(1), ppt(1), one /) ! least squares fit summation;
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!  enddo ! end of do jj = 1,nPpts

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!  if( lerror.eq.0 ) then 
    
!   write(lunit+myid)poincaredata(1:5,0:Nz-1,1:nPpts) ! write field line trajectory to file;
    
!   if( ii.ge.0 ) then ! fit straight line, least squares fit to trajectory;
!    fittedtransform(ii,1:2) = (/ sti(1), ( leastsqfit(5)*leastsqfit(3)-leastsqfit(2)*leastsqfit(4) ) / ( leastsqfit(5)*leastsqfit(1)-leastsqfit(2)*leastsqfit(2) ) /)
!   endif

!  endif

!  if( Wpp00aa ) then
!   cput = GETTIME
!   if( ii.eq.0                  ) write(ounit,1000)cput-cpus,myid,lvol,ii,lerror,sti,fittedtransform(ii,2),diota(lvol,0,0),fittedtransform(ii,2)-diota(lvol,0,0)
!   if( ii.gt.0.and.ii.lt.lnPtrj ) write(ounit,1000)cput-cpus,myid,lvol,ii,lerror,sti,fittedtransform(ii,2)
!   if(             ii.eq.lnPtrj ) write(ounit,1000)cput-cpus,myid,lvol,ii,lerror,sti,fittedtransform(ii,2),diota(lvol,0,1),fittedtransform(ii,2)-diota(lvol,0,1)
!  endif

!1000 format("pp00aa : ",f10.2," : myid=",i3," ; lvol=",i3," , ii="i4" ; lerror="i1" ; (s,t)= ("2f9.5" ) ; iota="2es23.15" ("es10.2") ;")

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
! enddo ! end of do ii = - qfms(lvol,1)%i * ( Nghd*2 + 1 ), lnPtrj ; iigh = ii + ( Nghd*2 + 1 )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! DEALLOCATE(qfmsst)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! Ltangent = 1
