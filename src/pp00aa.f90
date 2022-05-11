!> \file
!> \brief Constructs Poincaré plot and "approximate" rotational-transform (driver).

!> \brief Constructs Poincaré plot and "approximate" rotational-transform (driver).
!> \ingroup grp_diagnostics
!>
!> **relevant input variables**
!>
!> <ul>
!> <li> The resolution of Poincaré plot is controlled by
!>       <ul>
!>       <li> \c nPtraj trajectories will be located in each volume;
!>       <li> \c nPpts  iterations per trajectory;
!>       <li> \c odetol o.d.e. integration tolerance;
!>       </ul> </li>
!> <li> The magnetic field is given by bfield() . </li>
!> <li> The approximate rotational transform is determined, in pp00ab() , by fieldline integration. </li>
!> </ul>
!>
!> **format of output: Poincaré**
!>
!> <ul>
!> <li> The Poincaré data is written to \c .ext.poincare:xxxx , where \c xxxx is an integer indicating the volume.
!>       The format of this file is as follows:
!>
!>~~~~~~~~~~~~
!> write(svol,'(i4.4)')lvol ! lvol labels volume;
!> open(lunit+myid,file="."//trim(ext)//".poincare."//svol,status="unknown",form="unformatted")
!> do until end of file
!>   write(lunit+myid) Nz, nPpts                ! integers
!>   write(lunit+myid) data(1:4,0:Nz-1,1:nPpts) ! doubles
!> enddo
!> close(lunit+myid)
!>~~~~~~~~~~~~
!> `where`
!>       <ul>
!>       <li> \f$\theta \equiv\,\f$\c data(1,k,j) is the poloidal angle,      </li>
!>       <li> \f$ s     \equiv\,\f$\c data(2,k,j) is the radial coordinate,   </li>
!>       <li> \f$ R     \equiv\,\f$\c data(3,k,j) is the cylindrical \f$R\f$, </li>
!>       <li> \f$ Z     \equiv\,\f$\c data(4,k,j) is the cylindrical \f$Z\f$, </li>
!>       </ul>
!> <li> The integer \c k=0,Nz-1 labels toroidal planes, so that \f$\phi = ( 2 \pi / \texttt{Nfp}) ( k / \texttt{Nz})\f$,
!> <li> The integer \c j=1,nPpts labels toroidal iterations.
!> <li> Usually (if no fieldline integration errors are encountered) the number of fieldlines followed in volume \c lvol
!>       is given by \f$N+1\f$, where the radial resolution, \f$N\equiv\,\f$\c Ni(lvol) , is given on input.
!>       This will be over-ruled by if \c nPtrj(lvol) , given on input, is non-negative.
!> <li> The starting location for the fieldline integrations are equally spaced in the radial coordinate \f$s_i=s_{l-1}+ i (s_{l}-s_{l-1})/N\f$ for \f$i=0,N\f$,
!>       along the line \f$\theta=0\f$, \f$\zeta=0\f$.
!> </ul>
!>
!> **format of output: rotational-transform**
!>
!> <ul>
!>
!> <li> The rotational-transform data is written to \c .ext.transform:xxxx , where \c xxxx is an integer indicating the volume.
!>       The format of this file is as follows:
!> ```
!>  open(lunit+myid,file="."//trim(ext)//".sp.t."//svol,status="unknown",form="unformatted")
!>  write(lunit+myid) lnPtrj-ioff+1                                      ! integer
!>  write(lunit+myid) diotadxup(0:1,0,lvol)                              ! doubles
!>  write(lunit+myid) ( fiota(itrj,1:2), itrj = ioff, lnPtrj ) ! doubles
!>  close(lunit+myid)
!> ```
!> </ul>
!>
subroutine pp00aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  use constants, only : zero, half, one, two, pi

  use numerical, only :

  use fileunits, only : ounit

  use inputlist, only : Wmacros, Wpp00aa, Nvol, Lrad, odetol, nPpts, Ppts, nPtrj, Lconstraint, iota, oita, Igeometry

  use cputiming, only : Tpp00aa

  use allglobal, only : myid, ncpu, cpus, MPI_COMM_SPEC, ext, &
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

  integer :: id, numTraj, recvId
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
!    write(*,*) "we have numTrajs(",vvol,")=",numTrajs(vvol)
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

      !write(*,'(ai2ai2ai2ai2)') "CPU ",myid," works on trajectories ",ioff," to ",lnPtrj," in volume ",vvol

      SALLOCATE(   data, (ioff:lnPtrj, 1:4,0:Nz-1,1:nPpts), zero ) ! for block writing to file (allows faster reading of output data files for post-processing plotting routines);
      SALLOCATE( utflag, (ioff:lnPtrj                    ),    0 ) ! error flag that indicates if fieldlines successfully followed; 22 Apr 13;
      SALLOCATE(  fiota, (ioff:lnPtrj, 1:2               ), zero ) ! will always need fiota(0,1:2);

!$OMP PARALLEL DO SHARED(lnPtrj,ioff,Wpp00aa,Nz,data,fiota,utflag,iota,oita,myid,vvol,cpus,Lconstraint,nPpts,ppts) PRIVATE(itrj,sti)
      do itrj = ioff, lnPtrj ! initialize Poincare plot with trajectories regularly spaced between interfaces along \t=0;

        ; sti(1:2) = (/ - one + itrj    * two / lnPtrj   , Ppts*pi /)

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
!$OMP END PARALLEL DO

      ! write(*,*) "CPU ",myid," finished field line tracing for volume ",vvol
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

      ! TODO: replace below logic with a single call to MPI_allgather into rank-0
      ! and the write from there at once

      ! write data
      if (myid.eq.0) then
        !write(*,*) "CPU 0 writes its own Poincare data for numTrajs(",vvol,")=",numTrajs(vvol)

        ! The Poincare trajectories are dumped into the file one after another since this makes it easier to combine the results
        ! from all the nested volumes at the moment. Feel free to change the code, but verfiy that it works afterwards as well ;-)
        do itrj = ioff, lnPtrj
          ! write utflag --> success flag vector for field line tracing
          ! write data --> actual Poincare data
          if (vvol.gt.1) then
            !write(*,*) "CPU 0 writes a trajectory at offset ",sum(numTrajs(1:vvol-1))+itrj-ioff
            call write_poincare ( sum(numTrajs(1:vvol-1))+itrj-ioff, data(itrj,:,:,:), utflag )
          else
            !write(*,*) "CPU 0 writes a trajectory at offset ",itrj-ioff
            call write_poincare (                         itrj-ioff, data(itrj,:,:,:), utflag )
          endif
        enddo

        ! The rotational transform data is written at once for a volume
        if (vvol.gt.1) then
          call write_transform( sum(numTrajs(1:vvol-1)), numTrajs(vvol), vvol, diotadxup(0:1,0,vvol), fiota(0:lnPtrj,1:2) )
        else
          call write_transform( 0, numTrajs(1), 1, diotadxup(0:1,0,1), fiota(ioff:lnPtrj,1:2) )
        endif

        if (Mvol.gt.1 .and. ncpu.gt.1) then

         ! Gather data from all other parallelly running CPUs; there are min(ncpu-1, Mvol-vvol) of these in this iteration
         ! If we have so few CPUs that all of them need to perform multiple iteration over the set of volumes in batches of ncpu,
         ! there are a number ncpu-1 CPUs apart from the master who still have data that needs to be written before they can contine.
         ! In the last iteration of the batch processing, or if we have many CPUs but just not enough to cover all volumes in one batch,
         ! there are Mvol-vvol (where vvol is the value of vvol in CPU 0) other CPUs waiting to have their data written.
         do recvId = 1, min(ncpu-1, Mvol-vvol)

            ! lvol is the volume that the parallel CPU was working on
            lvol = vvol+recvId

            !write(*,'(ai2ai2)') "CPU 0 writes Poincare data of id ",recvId," for volume ",lvol

            deallocate(utflag)
            deallocate(data)
            deallocate(fiota)

            allocate(utflag(1:numTrajs(lvol)))
            allocate(  data(1:numTrajs(lvol),1:4,0:Nz-1,1:nPpts))
            allocate( fiota(1:numTrajs(lvol),1:2))

            call MPI_Recv( utflag, numTrajs(lvol)           , MPI_INTEGER         , recvId, lvol, MPI_COMM_SPEC, status, ierr)
            !write(*,*) "CPU 0 got utflag vector from CPU ",recvId

            call MPI_Recv(   data, numTrajs(lvol)*4*Nz*nPpts, MPI_DOUBLE_PRECISION, recvId, lvol, MPI_COMM_SPEC, status, ierr)
            !write(*,*) "CPU 0 got the corresponding Poincare data from CPU ",recvId

            call MPI_Recv(  fiota, numTrajs(lvol)*2         , MPI_DOUBLE_PRECISION, recvId, lvol, MPI_COMM_SPEC, status, ierr)
!            write(*,*) "CPU 0 got the iota profile from CPU ",recvId,": sarr: "
!            write(*,*) fiota(:,1)

            ! write utflag vector of CPU id
            ! write data of CPU id
            do itrj = 1, numTrajs(lvol)
              !write(*,*) "CPU 0 writes a trajectory at offset ",sum(numTrajs(1:lvol-1))+itrj-1
              call write_poincare( sum(numTrajs(1:lvol-1))+itrj-1, data(itrj,:,:,:), utflag )
            enddo

            ! write fiota --> iota from field line tracing
            ! write diotadxup --> radial derivative of iota from Beltrami field (?)
            call write_transform( sum(numTrajs(1:lvol-1)), numTrajs(lvol), lvol, diotadxup(0:1,0,lvol), fiota(1:numTrajs(lvol),1:2) )

            !write(*,*) "wrote iota data at offset ",sum(numTrajs(1:lvol-1))," with length ",numTrajs(lvol)

            ! write fiota of CPU id
          enddo ! vvol = 2, Mvol
        endif

      else
        call MPI_Send( utflag, numTrajs(vvol)           , MPI_INTEGER         , 0, vvol, MPI_COMM_SPEC, ierr) ! success flag vector
        call MPI_Send(   data, numTrajs(vvol)*4*Nz*nPpts, MPI_DOUBLE_PRECISION, 0, vvol, MPI_COMM_SPEC, ierr) ! Poincare data
        call MPI_Send(  fiota, numTrajs(vvol)*2         , MPI_DOUBLE_PRECISION, 0, vvol, MPI_COMM_SPEC, ierr) ! rotational transform profile from field line tracing
        ! diotadxup should be available in the master already, since it is stored in global
      endif ! myid.eq.0

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

      DALLOCATE(data)
      DALLOCATE(utflag)
      DALLOCATE(fiota)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    endif ! myid.eq.modulo(vvol-1,ncpu)
  enddo ! vvol = 1, Mvol

  ! keep this until the end because we need it if we write with less CPUs than nested volumes
  deallocate(numTrajs)

  if (myid.eq.0) then
    call finalize_flt_output
  endif

  RETURN(pp00aa)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

1001 format("pp00aa : ",f10.2," : myid=",i3," ; lvol=",i3," ; odetol=",es8.1," ; nPpts=",i8," ; lnPtrj=",i3," ;")
1002 format("pp00aa : ",f10.2," : myid=",i3," ; lvol=",i3," ; ",i3," : (s,t)=(",f21.17," ,",f21.17," ) ;":" utflag=",i3," ; transform=",es23.15,&
  " ;":" error=",es13.5," ;")

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine pp00aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
