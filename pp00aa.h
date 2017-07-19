!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (diagnostic) ! Constructs Poincar&eacute; plot and &ldquo;approximate&rdquo; rotational-transform (driver).

!latex \briefly{Constructs \Poincare plot and ``approximate" rotational-transform (driver).}

!latex \calledby{\link{xspech}}
!latex \calls{\link{pp00ab}}

!latex \tableofcontents

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \subsection{relevant input variables}

!latex \begin{enumerate}
!latex \item The resolution of \Poincare plot is controlled by 
!latex       \begin{itemize}
!latex       \item[i.] \inputvar{nPtraj} trajectories will be located in each volume;
!latex       \item[ii.] \inputvar{nPpts} iterations per trajectory;
!latex       \item[iii.] \inputvar{odetol} o.d.e. integration tolerance;
!latex       \end{itemize}
!latex \item The magnetic field is given by \link{bfield}.
!latex \item The approximate rotational transform is determined, in \link{pp00ab}, by fieldline integration.
!latex  \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{format of output: \Poincare}

!latex \begin{enumerate}
!latex \item The \Poincare data is written to \type{.ext.poincare:xxxx}, where \type{xxxx} is an integer indicating the volume.
!latex       The format of this file is as follows:
!latex \begin{verbatim}
!latex  write(svol,'(i4.4)')lvol ! lvol labels volume;
!latex  open(lunit+myid,file="."//trim(ext)//".poincare."//svol,status="unknown",form="unformatted")
!latex  do until end of file
!latex   write(lunit+myid) Nz, nPpts                ! integers
!latex   write(lunit+myid) data(1:4,0:Nz-1,1:nPpts) ! doubles
!latex  enddo
!latex  close(lunit+myid)
!latex \end{verbatim}
!latex       where \begin{itemize}
!latex       \item[i.] $\t \equiv$ \type{data(1,k,j)} is the poloidal angle,
!latex       \item[ii.] $\s \equiv$ \type{data(2,k,j)} is the radial coordinate, 
!latex       \item[iii.] $ R \equiv$ \type{data(3,k,j)} is the cylindrical $R$, 
!latex       \item[iv.]   $ Z \equiv$ \type{data(4,k,j)} is the cylindrical $Z$,
!latex       \end{itemize}
!latex \item The integer \type{k=0,Nz-1} labels toroidal planes, so that $\phi = ( 2 \pi / $\inputvar{Nfp}$ ) ( k / \type{Nz})$,
!latex \item The integer \type{j=1,}\inputvar{nPpts} labels toroidal iterations.
!latex \item Usually (if no fieldline integration errors are encountered) the number of fieldlines followed in volume \type{lvol}
!latex       is given by $N+1$, where the radial resolution, $N\equiv$\type{Ni(lvol)}, is given on input.
!latex       This will be over-ruled by if \inputvar{nPtrj(lvol)}, given on input, is non-negative.
!latex \item The starting location for the fieldline integrations are equally spaced in the radial coordinate $s_i=s_{l-1}+ i (s_{l}-s_{l-1})/N$ for $i=0,N$,
!latex       along the line $\t=0$, $\z=0$.
!latex \end{enumerate} 
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \subsection{format of output: rotational-transform}

!latex \begin{enumerate}
  
!latex \item The rotational-transform data is written to \type{.exttransform:xxxx}, where \type{xxxx} is an integer indicating the volume.
!latex       The format of this file is as follows:
!latex \begin{verbatim}
!latex  open(lunit+myid,file="."//trim(ext)//".sp.t."//svol,status="unknown",form="unformatted")
!latex  write(lunit+myid) lnPtrj-ioff+1                                      ! integer
!latex  write(lunit+myid) diotadxup(0:1,0,lvol)                              ! doubles
!latex  write(lunit+myid) ( fiota(itrj,1:2), itrj = ioff, lnPtrj ) ! doubles
!latex  close(lunit+myid)
!latex \end{verbatim}
!latex \end{enumerate}
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine pp00aa( lvol ) 
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half, one, two, pi
  
  use numerical, only :
  
  use fileunits, only : ounit, lunit
  
  use inputlist, only : Wmacros, Wpp00aa, Nvol, Lrad, ext, odetol, nPpts, nPtrj, Lconstraint, iota, oita
  
  use cputiming, only : Tpp00aa
  
  use allglobal, only : myid, ncpu, cpus, &
                        Nz, pi2nfp, &
                        ivol, Mvol, &
                        Lcoordinatesingularity, &
                        diotadxup
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in)  :: lvol
  INTEGER              :: lnPtrj, ioff, itrj
  INTEGER, allocatable :: id02bjf(:)
  REAL                 :: sti(1:2), ltransform(1:2)
  REAL, allocatable    :: data(:,:,:), fiota(:,:)
  CHARACTER            :: svol*4 ! perhaps this should be global; 11 Aug 13;
  
  BEGIN(pp00aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( nPtrj(lvol).ge.0 ) then ; lnPtrj =    nPtrj(lvol) ! selected Poincare resolution;
  else                        ; lnPtrj = 2 * Lrad(lvol) ! adapted  Poincare resolution;
  endif

  if( lnPtrj.le.0 ) goto 9999

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  ivol = lvol ; write(svol,'(i4.4)') lvol
  
  open( lunit+myid, file="."//trim(ext)//".sp.P."//svol//".dat", status="unknown", form="unformatted" )
  
1001 format("pp00aa : ",f10.2," : myid=",i3," ; lvol=",i3," ; odetol=",es8.1," ; nPpts=",i8," ; lnPtrj=",i3," ;")
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  SALLOCATE( data, (1:4,0:Nz-1,1:nPpts), zero ) ! for block writing to file (allows faster reading of output data files for post-processing plotting routines);
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
1002 format("pp00aa : ",f10.2," : myid=",i3," ; lvol=",i3," ; ",i3," : (s,t)=(",f21.17," ,",f21.17," ) ;":" id02bjf=",i3," ; transform=",es23.15,&
  " ;":" error=",es13.5," ;")
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Lcoordinatesingularity ) then ; ioff = 1 ! keep away from coordinate axis;
  else                              ; ioff = 0
  endif
  
  SALLOCATE( id02bjf, (ioff:lnPtrj), 0 ) ! error flag that indicates if fieldlines successfully followed; 22 Apr 13;
  
  SALLOCATE( fiota, (ioff:lnPtrj,1:2), zero ) ! will always need fiota(0,1:2);
  
  do itrj = ioff, lnPtrj ! initialize Poincare plot with trajectories regularly spaced between interfaces along \t=0;
   
   if( Lcoordinatesingularity ) then ; sti(1:2) = (/ - one + itrj**2 * two / lnPtrj**2, zero /) ! equal increments in rr = \sqrt(ss) ; 08 Feb 16;
   else                              ; sti(1:2) = (/ - one + itrj    * two / lnPtrj   , zero /)
   endif
   
   if( itrj.eq.     0 ) sti(1) = - one ! avoid machine precision errors; 08 Feb 16;
   if( itrj.eq.lnPtrj ) sti(1) =   one ! avoid machine precision errors; 08 Feb 16;
   
   CALL( pp00aa, pp00ab, ( lvol, sti(1:2), Nz, nPpts, data(1:4,0:Nz-1,1:nPpts), fiota(itrj,1:2), id02bjf(itrj) ) )
   
   if( Wpp00aa ) then
    cput = GETTIME
    if( Lconstraint.eq.1 ) then
     if( itrj.eq.0                      ) write(ounit,1002) cput-cpus, myid, lvol, itrj, sti(1:2), id02bjf(itrj), fiota(itrj,2), fiota(itrj,2)-oita(lvol-1)
     if( itrj.gt.0 .and. itrj.lt.lnPtrj ) write(ounit,1002) cput-cpus, myid, lvol, itrj, sti(1:2), id02bjf(itrj), fiota(itrj,2)
     if(                 itrj.eq.lnPtrj ) write(ounit,1002) cput-cpus, myid, lvol, itrj, sti(1:2), id02bjf(itrj), fiota(itrj,2), fiota(itrj,2)-iota(lvol  )
    else                                ; write(ounit,1002) cput-cpus, myid, lvol, itrj, sti(1:2), id02bjf(itrj), fiota(itrj,2)
    endif
   endif
   
   if( id02bjf(itrj).eq.0 ) then ! will only write successfully followed trajectories to file; 28 Jan 13;
    write(lunit+myid) Nz, nPpts
    write(lunit+myid) data(1:4,0:Nz-1,1:nPpts)
   endif
   
  enddo ! end of do itrj; 25 Jan 13;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  DALLOCATE(data)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  close( lunit+myid ) ! Poincare plot is finished; will re-use lunit to write rotational transform data;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  open( lunit+myid, file="."//trim(ext)//".sp.t."//svol//".dat", status="unknown", form="unformatted" ) ! rotational transform data;
  write( lunit+myid ) lnPtrj-ioff+1
  write( lunit+myid ) diotadxup(0:1,0,lvol)
  write( lunit+myid ) ( fiota(itrj,1:2), itrj = ioff, lnPtrj )
  close( lunit+myid )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  DALLOCATE(id02bjf)
  
  DALLOCATE(fiota)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(pp00aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine pp00aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
