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

subroutine pp00aa
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half, one, two, pi
  
  use numerical, only :
  
  use fileunits, only : ounit
  
  use inputlist, only : Wmacros, Wpp00aa, Nvol, Lrad, ext, odetol, nPpts, nPtrj, Lconstraint, iota, oita, Igeometry
  
  use cputiming, only : Tpp00aa
  
  use allglobal, only : myid, ncpu, cpus, &
                        Nz, pi2nfp, &
                        ivol, Mvol, &
                        Lcoordinatesingularity, &
                        diotadxup, Lplasmaregion, Lvacuumregion

  use sphdf5,    only : init_flt_output, write_poincare, write_transform, finalize_flt_output
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER              :: lnPtrj, ioff, vvol, itrj, lvol
  INTEGER, allocatable :: utflag(:), numTrajs(:)
  REAL                 :: sti(1:2), ltransform(1:2)
  REAL, allocatable    :: data(:,:,:,:), fiota(:,:)
  
  integer :: id, numTraj
  integer :: status(MPI_STATUS_SIZE)

  BEGIN(pp00aa)
  
  ! count how many Poincare trajectories should be computed in total ; executed on each CPU
  allocate(numTrajs(1:Mvol))
  do vvol = 1, Mvol
    LREGION(vvol) ! sets e.g. Lcoordinatesingularity
    if( Lcoordinatesingularity ) then ; ioff = 1 ! keep away from coordinate axis;
    else                              ; ioff = 0
    endif

    if( nPtrj(vvol).ge.0 ) then ; lnPtrj =    nPtrj(vvol) ! selected Poincare resolution;
    else                        ; lnPtrj = 2 * Lrad(vvol) ! adapted  Poincare resolution;
    endif

    numTrajs(vvol) = lnPtrj - ioff + 1 ! interval includes edge
  enddo

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  ! the total number of Poincare trajectories to be saved later on is given by sum(numTrajs)
  call init_flt_output( sum(numTrajs) )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  ! start tracing
  do vvol = 1, Mvol

    ! skip loop content if in the wrong CPU
    if( myid.eq.modulo(vvol-1,ncpu) .and. myid.lt.Mvol) then ! the following is in parallel; 20 Jun 14;

      ! lower bound for radial indices
      LREGION(vvol) ! sets e.g. Lcoordinatesingularity
      if( Lcoordinatesingularity ) then ; ioff = 1 ! keep away from coordinate axis;
      else                              ; ioff = 0
      endif

      ! upper bound for radial indices
      if( nPtrj(vvol).ge.0 ) then ; lnPtrj =    nPtrj(vvol) ! selected Poincare resolution;
      else                        ; lnPtrj = 2 * Lrad(vvol) ! adapted  Poincare resolution;
      endif

      write(*,*) "CPU ",myid," works on trajectories ",ioff," to ",lnPtrj

      SALLOCATE(   data, (ioff:lnPtrj, 1:4,0:Nz-1,1:nPpts), zero ) ! for block writing to file (allows faster reading of output data files for post-processing plotting routines);
      SALLOCATE( utflag, (ioff:lnPtrj                    ),    0 ) ! error flag that indicates if fieldlines successfully followed; 22 Apr 13;
      SALLOCATE(  fiota, (ioff:lnPtrj, 1:2               ), zero ) ! will always need fiota(0,1:2);

      do itrj = ioff, lnPtrj ! initialize Poincare plot with trajectories regularly spaced between interfaces along \t=0;

        if( Lcoordinatesingularity ) then ; sti(1:2) = (/ - one + itrj**2 * two / lnPtrj**2, zero /) ! equal increments in rr = \sqrt(ss) ; 08 Feb 16;
        else                              ; sti(1:2) = (/ - one + itrj    * two / lnPtrj   , zero /)
        endif

        if( itrj.eq.     0 ) sti(1) = - one ! avoid machine precision errors; 08 Feb 16;
        if( itrj.eq.lnPtrj ) sti(1) =   one ! avoid machine precision errors; 08 Feb 16;

        ! call actual field line integration subroutine
        CALL( pp00aa, pp00ab, ( vvol, sti(1:2), Nz, nPpts, data(itrj,1:4,0:Nz-1,1:nPpts), fiota(itrj,1:2), utflag(itrj) ) )

        if( Wpp00aa ) then
          cput = GETTIME
          if( Lconstraint.eq.1 ) then
            if( itrj.eq.0                      ) write(ounit,1002) cput-cpus, myid, vvol, itrj, sti(1:2), utflag(itrj), fiota(itrj,2), fiota(itrj,2)-oita(vvol-1)
            if( itrj.gt.0 .and. itrj.lt.lnPtrj ) write(ounit,1002) cput-cpus, myid, vvol, itrj, sti(1:2), utflag(itrj), fiota(itrj,2)
            if(                 itrj.eq.lnPtrj ) write(ounit,1002) cput-cpus, myid, vvol, itrj, sti(1:2), utflag(itrj), fiota(itrj,2), fiota(itrj,2)-iota(vvol  )
          else                                 ; write(ounit,1002) cput-cpus, myid, vvol, itrj, sti(1:2), utflag(itrj), fiota(itrj,2)
          endif
        endif ! Wpp00aa

      enddo ! itrj = ioff, lnPtrj

      write(*,*) "CPU ",myid," finished field line tracing for volume ",vvol
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

      ! write data
      if (myid.eq.0) then
        write(*,*) "CPU 0 writes its Poincare data"

        do itrj = ioff, lnPtrj
          ! write utflag --> success flag vector for field line tracing
          ! write data --> actual Poincare data
          if (vvol.gt.1) then
            write(*,*) "CPU 0 writes a trajectory at offset ",sum(numTrajs(1:vvol-1))+itrj-ioff
            call write_poincare ( sum(numTrajs(1:vvol-1))+itrj-ioff, data(itrj,:,:,:), utflag )
          else
            write(*,*) "CPU 0 writes a trajectory at offset ",itrj-ioff
            call write_poincare (                       itrj-ioff, data(itrj,:,:,:), utflag )
          endif
        enddo
        call write_transform( 0, numTrajs(1), 1, diotadxup(0:1,0,1), fiota(1:numTrajs(1),1:2) )

        if (Mvol.gt.1) then
          do lvol = 2, Mvol
            write(*,*) "CPU 0 writes Poincare data for volume ",lvol

            deallocate(utflag)
            deallocate(data)
            deallocate(fiota)

            allocate(utflag(1:numTrajs(lvol)))
            allocate(  data(1:numTrajs(lvol),1:4,0:Nz-1,1:nPpts))
            allocate( fiota(1:numTrajs(lvol),1:2))

            call MPI_Recv( utflag, numTrajs(lvol)           , MPI_INTEGER         , modulo(lvol-1,ncpu), lvol, MPI_COMM_WORLD, status, ierr)
            write(*,*) "CPU 0 got utflag vector from CPU ",modulo(lvol-1,ncpu)

            call MPI_Recv(   data, numTrajs(lvol)*4*Nz*nPpts, MPI_DOUBLE_PRECISION, modulo(lvol-1,ncpu), lvol, MPI_COMM_WORLD, status, ierr)
            write(*,*) "CPU 0 got the corresponding Poincare data from CPU ",modulo(lvol-1,ncpu)

            call MPI_Recv(  fiota, numTrajs(lvol)*2         , MPI_DOUBLE_PRECISION, modulo(lvol-1,ncpu), lvol, MPI_COMM_WORLD, status, ierr)
            write(*,*) "CPU 0 got the iota profile from CPU ",modulo(lvol-1,ncpu)

            ! write utflag vector of CPU id
            ! write data of CPU id

            do itrj = 1, numTrajs(lvol)
              write(*,*) "CPU 0 writes a trajectory at offset ",sum(numTrajs(1:lvol-1))+itrj-1
              call write_poincare( sum(numTrajs(1:lvol-1))+itrj-1, data(itrj,:,:,:), utflag )
            enddo

            ! write fiota --> iota from field line tracing
            ! write diotadxup --> iota from Beltrami field(?)
            !call write_transform( sum(numTrajs(1:lvol-1)), numTrajs(lvol), lvol, diotadxup(0:1,0,lvol), fiota(1:numTrajs(lvol),1:2) )
            call write_transform( sum(numTrajs(1:lvol-1)), numTrajs(lvol), lvol, diotadxup(0:1,0,lvol), fiota(1:numTrajs(lvol),1:2) )

            ! write fiota of CPU id
          enddo ! vvol = 2, Mvol
        endif

        deallocate(numTrajs)
      else
        call MPI_Send( utflag, numTrajs(vvol)           , MPI_INTEGER         , 0, vvol, MPI_COMM_WORLD, ierr) ! success flag vector
        call MPI_Send(   data, numTrajs(vvol)*4*Nz*nPpts, MPI_DOUBLE_PRECISION, 0, vvol, MPI_COMM_WORLD, ierr) ! Poincare data
        call MPI_Send(  fiota, numTrajs(vvol)*2         , MPI_DOUBLE_PRECISION, 0, vvol, MPI_COMM_WORLD, ierr) ! rotational transform profile from field line tracing
        ! diotadxup should be available in the master already, since it is stored in global
      endif ! myid.eq.0

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

      DALLOCATE(data)
      DALLOCATE(utflag)
      DALLOCATE(fiota)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    endif ! myid.eq.modulo(vvol-1,ncpu)
  enddo ! vvol = 1, Mvol



  if (myid.eq.0) then
    call finalize_flt_output
  endif




       ! write all trajectories, but only mark successfully followed trajectories with success.eq.1; 21 May 19;
      ! WCALL( pp00aa, write_poincare, (data, numTrajTotal+itrj-ioff, utflag(itrj)) )

      ! write rotational transform data to output file
      !WCALL( pp00aa, write_transform, (numTrajTotal, lnPtrj-ioff+1, lvol, diotadxup(0:1,0,lvol), fiota(ioff:lnPtrj,1:2)) ) ! 21 May 19;


  







  RETURN(pp00aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

1001 format("pp00aa : ",f10.2," : myid=",i3," ; lvol=",i3," ; odetol=",es8.1," ; nPpts=",i8," ; lnPtrj=",i3," ;")
1002 format("pp00aa : ",f10.2," : myid=",i3," ; lvol=",i3," ; ",i3," : (s,t)=(",f21.17," ,",f21.17," ) ;":" utflag=",i3," ; transform=",es23.15,&
  " ;":" error=",es13.5," ;")

end subroutine pp00aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!