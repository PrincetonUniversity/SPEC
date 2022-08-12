!> \defgroup grp_parallel Parallelization
!>
!> \file
!> \brief Broadcasts Beltrami fields, profiles, . . .

!> \ingroup grp_parallel
!> \brief Broadcasts Beltrami fields, profiles, . . .
!>
!> **broadcasting**
!> <ul>
!> <li> The construction of the Beltrami fields is distributed on separate cpus. </li>
!> <li> All "local" information needs to be broadcast so that the "global" force vector,
!>       \f{eqnarray}{ {\bf F}_i \equiv [[p+B^2/2]]_i = (p+B^2/2)_{v,i} - (p+B^2/2)_{v-1,i}
!>       \f}
!>       can be constructed, and so that restart and output files can be saved to file. </li>
!> </ul>
!> @param[in] lvol index of nested volume
subroutine brcast( lvol )
  use mod_kinds, only: wp => dp
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  use constants, only : zero

  use numerical, only :

  use fileunits, only : ounit

  use inputlist, only : Wmacros, Wbrcast, Wcurent, MNvol, Nvol, mu, Lrad, &
                        curtor, curpol, Lconstraint, Lfindzero, helicity

  use cputiming, only : Tbrcast

  use allglobal, only : myid, cpus, ncpu, MPI_COMM_SPEC, &
                        dtflux, dpflux, Ntz, mn, Mvol, &
                        diotadxup, dItGpdxtp, &
                        Ate, Aze, Ato, Azo, &
                        Bemn, Bomn, Iomn, Iemn, Somn, Semn, Pomn, Pemn, &
                        ImagneticOK, &
                       !dBBdRZ, dIIdRZ, &
                        Lhessianallocated, LGdof, dFFdRZ, dBBdmp, dmupfdx, &
                        denergydrr,denergydzr,Lhessian3Dallocated, &
                        lBBintegral, lABintegral, &
                        vvolume, &
                        NOTstellsym, LocalConstraint, &
						IsMyVolume, IsMyVolumeValue

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


#ifdef OPENMP
  USE OMP_LIB
#endif
  use mpi
  implicit none
  integer   :: ierr, astat, ios, nthreads, ithread
  real(wp)      :: cput, cpui, cpuo=0 ! cpu time; cpu initial; cpu old; 31 Jan 13;


  integer, intent(in) :: lvol

  integer             :: llmodnp, io, iRZl, ii, ideriv, Nbc


  cpui = MPI_WTIME()
  cpuo = cpui
#ifdef OPENMP
  nthreads = omp_get_max_threads()
#else
  nthreads = 1
#endif


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! recall this routine is inside do vvol = 1, Mvol loop; see dforce;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


   if( lvol.le.0 .or. lvol.gt.Mvol ) then
     write(6,'("brcast :      fatal : myid=",i3," ; lvol.le.0 .or. lvol.gt.Mvol ; error ;")') myid
     call MPI_ABORT( MPI_COMM_SPEC, 1, ierr )
     stop "brcast : lvol.le.0 .or. lvol.gt.Mvol : error  ;"
    endif


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  llmodnp = modulo(lvol-1,ncpu) ! identify which node contains data; this must be consistent with previous looping / parallelization;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


   call MPI_BCAST(mu(lvol),1,MPI_DOUBLE_PRECISION,llmodnp ,MPI_COMM_SPEC,ierr)


   call MPI_BCAST(dtflux(lvol),1,MPI_DOUBLE_PRECISION,llmodnp ,MPI_COMM_SPEC,ierr)


   call MPI_BCAST(dpflux(lvol),1,MPI_DOUBLE_PRECISION,llmodnp ,MPI_COMM_SPEC,ierr)


   call MPI_BCAST(helicity(lvol),1,MPI_DOUBLE_PRECISION,llmodnp,MPI_COMM_SPEC,ierr)



   call MPI_BCAST(vvolume(lvol),1,MPI_DOUBLE_PRECISION,llmodnp ,MPI_COMM_SPEC,ierr)


   call MPI_BCAST(lBBintegral(lvol),1,MPI_DOUBLE_PRECISION,llmodnp ,MPI_COMM_SPEC,ierr)


   call MPI_BCAST(lABintegral(lvol),1,MPI_DOUBLE_PRECISION,llmodnp ,MPI_COMM_SPEC,ierr)



   call MPI_BCAST(diotadxup(0:1,-1:2,lvol),8,MPI_DOUBLE_PRECISION,llmodnp ,MPI_COMM_SPEC,ierr)


   call MPI_BCAST(dItGpdxtp(0:1,-1:2,lvol),8,MPI_DOUBLE_PRECISION,llmodnp ,MPI_COMM_SPEC,ierr)


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( Lhessianallocated ) then


   if( LocalConstraint ) then
 	  Nbc =             LGdof*       2*  LGdof*  2

   call MPI_BCAST(dFFdRZ(1:LGdof,0:1,1:LGdof,0:1,lvol),Nbc,MPI_DOUBLE_PRECISION,llmodnp ,MPI_COMM_SPEC,ierr)



	  Nbc =             LGdof*       2*  2

   call MPI_BCAST(dBBdmp(1:LGdof,lvol,0:1,1:2),Nbc,MPI_DOUBLE_PRECISION,llmodnp ,MPI_COMM_SPEC,ierr)


	  Nbc =                   2*  LGdof*  2

   call MPI_BCAST(dmupfdx(lvol,1:1   ,1:2,1:LGdof,0:1),Nbc,MPI_DOUBLE_PRECISION,llmodnp ,MPI_COMM_SPEC,ierr)
 ! why is this broadcast; 02 Sep 14;
   endif


  endif ! end of if( Lhessianallocated ) ; 12 Sep 16;

  if (Lhessian3Dallocated) then

      Nbc =             LGdof*       2*  LGdof*  2

   call MPI_BCAST(denergydrr(1:LGdof,lvol,0:1,1:LGdof,0:1),Nbc,MPI_DOUBLE_PRECISION,llmodnp ,MPI_COMM_SPEC,ierr)


   call MPI_BCAST(denergydzr(1:LGdof,lvol,0:1,1:LGdof,0:1),Nbc,MPI_DOUBLE_PRECISION,llmodnp ,MPI_COMM_SPEC,ierr)

  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


   call MPI_BCAST(ImagneticOK(lvol),1,MPI_LOGICAL,llmodnp ,MPI_COMM_SPEC,ierr)


  ! Commented - broadcasted in dfp200
  ! do ideriv = 0, 2
  ! 	if( (ideriv.ne.0) .and. (Lconstraint.ne.3) ) cycle
  !     do ii = 1, mn
  !       RlBCAST( Ate(lvol,ideriv,ii)%s(0:Lrad(lvol)), Lrad(lvol)+1, llmodnp )
  !       RlBCAST( Aze(lvol,ideriv,ii)%s(0:Lrad(lvol)), Lrad(lvol)+1, llmodnp )
  !     enddo
  ! enddo



   call MPI_BCAST(Bemn(1:mn,lvol,0:1),2*mn,MPI_DOUBLE_PRECISION,llmodnp ,MPI_COMM_SPEC,ierr)
 ! perhaps all these should be re-ordered; 18 Jul 14;

   call MPI_BCAST(Iomn(1:mn,lvol    ),mn,MPI_DOUBLE_PRECISION,llmodnp ,MPI_COMM_SPEC,ierr)


   call MPI_BCAST(Somn(1:mn,lvol,0:1),2*mn,MPI_DOUBLE_PRECISION,llmodnp ,MPI_COMM_SPEC,ierr)


   call MPI_BCAST(Pomn(1:mn,lvol,0:2),3*mn,MPI_DOUBLE_PRECISION,llmodnp ,MPI_COMM_SPEC,ierr)
 ! 15 Sep 15;

  if( NOTstellsym ) then
    ! do ideriv = 0, 2
    !   do ii = 1, mn
    !     RlBCAST( Ato(lvol,ideriv,ii)%s(0:Lrad(lvol)), Lrad(lvol)+1, llmodnp )
    !     RlBCAST( Azo(lvol,ideriv,ii)%s(0:Lrad(lvol)), Lrad(lvol)+1, llmodnp )
    !   enddo
    ! enddo


   call MPI_BCAST(Bomn(1:mn,lvol,0:1),2*mn,MPI_DOUBLE_PRECISION,llmodnp ,MPI_COMM_SPEC,ierr)


   call MPI_BCAST(Iemn(1:mn,lvol    ),mn,MPI_DOUBLE_PRECISION,llmodnp ,MPI_COMM_SPEC,ierr)


   call MPI_BCAST(Semn(1:mn,lvol,0:1),2*mn,MPI_DOUBLE_PRECISION,llmodnp ,MPI_COMM_SPEC,ierr)


   call MPI_BCAST(Pemn(1:mn,lvol,0:2),3*mn,MPI_DOUBLE_PRECISION,llmodnp ,MPI_COMM_SPEC,ierr)

  endif ! end of if( NOTstellsym) ; 11 Aug 14;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! if( lvol.gt.Nvol .and. Lconstraint.eq.-1 .and. Wcurent ) then ! 27 Feb 17;
  if( lvol.gt.Nvol                         .and. Wcurent ) then ! 27 Feb 17;
  !write(ounit,'("brcast : " 10x " : myid="i3" ; broadcasting : curtor="es13.5" ; curpol="es13.5" ;")') myid, curtor, curpol

   call MPI_BCAST(curtor,1,MPI_DOUBLE_PRECISION,llmodnp ,MPI_COMM_SPEC,ierr)


   call MPI_BCAST(curpol,1,MPI_DOUBLE_PRECISION,llmodnp ,MPI_COMM_SPEC,ierr)

  endif
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


9999 continue
  cput = MPI_WTIME()
  Tbrcast = Tbrcast + ( cput-cpuo )
  return


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine brcast

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
