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
                        ImagneticOK, Ldescent, &
                       !dBBdRZ, dIIdRZ, &
                        Lhessianallocated, LGdof, dFFdRZ, dBBdmp, dmupfdx, &
                        lBBintegral, lABintegral, lMMl, lLLl, &
                        vvolume, &
                        NOTstellsym, LocalConstraint, &
						IsMyVolume, IsMyVolumeValue

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS

  INTEGER, intent(in) :: lvol

  INTEGER             :: llmodnp, io, iRZl, ii, ideriv, Nbc

  BEGIN(brcast)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! recall this routine is inside do vvol = 1, Mvol loop; see dforce;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  FATAL( brcast, lvol.le.0 .or. lvol.gt.Mvol, error )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  llmodnp = modulo(lvol-1,ncpu) ! identify which node contains data; this must be consistent with previous looping / parallelization;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RlBCAST(     mu(lvol), 1, llmodnp )
  RlBCAST( dtflux(lvol), 1, llmodnp )
  RlBCAST( dpflux(lvol), 1, llmodnp )
  RlBCAST( helicity(lvol), 1, llmodnp)

  RlBCAST(     vvolume(lvol), 1, llmodnp )
  RlBCAST( lBBintegral(lvol), 1, llmodnp )
  RlBCAST( lABintegral(lvol), 1, llmodnp )
  RlBCAST( lMMl(lvol), 1, llmodnp )
  RlBCAST( lLLl(lvol), 1, llmodnp )
    

  RlBCAST( diotadxup(0:1,-1:2,lvol), 8, llmodnp )
  RlBCAST( dItGpdxtp(0:1,-1:2,lvol), 8, llmodnp )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( Lhessianallocated ) then


   if( LocalConstraint ) then
 	  Nbc =             LGdof*       2*  LGdof*  2
 	  RlBCAST( dFFdRZ(1:LGdof,0:1,1:LGdof,0:1,lvol), Nbc, llmodnp )

	  Nbc =             LGdof*       2*  2
	  RlBCAST( dBBdmp(1:LGdof,lvol,0:1,1:2), Nbc, llmodnp )

	  Nbc =                   2*  LGdof*  2
	  RlBCAST( dmupfdx(lvol,1:1   ,1:2,1:LGdof,0:1), Nbc, llmodnp ) ! why is this broadcast; 02 Sep 14;
   endif


  endif ! end of if( Lhessianallocated ) ; 12 Sep 16;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LlBCAST( ImagneticOK(lvol), 1, llmodnp )

  ! Commented - broadcasted in dfp200
  ! do ideriv = 0, 2
  ! 	if( (ideriv.ne.0) .and. (Lconstraint.ne.3) ) cycle
  !     do ii = 1, mn
  !       RlBCAST( Ate(lvol,ideriv,ii)%s(0:Lrad(lvol)), Lrad(lvol)+1, llmodnp )
  !       RlBCAST( Aze(lvol,ideriv,ii)%s(0:Lrad(lvol)), Lrad(lvol)+1, llmodnp )
  !     enddo
  ! enddo

  if (Ldescent) then
    RlBCAST( Bemn(1:mn,lvol,0:3), 4*mn, llmodnp ) ! perhaps all these should be re-ordered; 18 Jul 14;
    RlBCAST( Iomn(1:mn,lvol,0:3), 4*mn, llmodnp )
    RlBCAST( Somn(1:mn,lvol,0:3), 4*mn, llmodnp )
    RlBCAST( Pomn(1:mn,lvol,0:3), 4*mn, llmodnp ) ! 15 Sep 15;

    RlBCAST( Bomn(1:mn,lvol,0:3), 4*mn, llmodnp ) ! perhaps all these should be re-ordered; 18 Jul 14;
    RlBCAST( Iemn(1:mn,lvol,0:3), 4*mn, llmodnp )
    RlBCAST( Semn(1:mn,lvol,0:3), 4*mn, llmodnp )
    RlBCAST( Pemn(1:mn,lvol,0:3), 4*mn, llmodnp ) ! 15 Sep 15;
  else
    RlBCAST( Bemn(1:mn,lvol,0:1), 2*mn, llmodnp )
    RlBCAST( Iomn(1:mn,lvol,0:0),   mn, llmodnp )
    RlBCAST( Somn(1:mn,lvol,0:1), 2*mn, llmodnp )
    RlBCAST( Pomn(1:mn,lvol,0:2), 3*mn, llmodnp )

    if( NOTstellsym ) then
      RlBCAST( Bomn(1:mn,lvol,0:1), 2*mn, llmodnp )
      RlBCAST( Iemn(1:mn,lvol,0:0),   mn, llmodnp )
      RlBCAST( Semn(1:mn,lvol,0:1), 2*mn, llmodnp )
      RlBCAST( Pemn(1:mn,lvol,0:2), 3*mn, llmodnp )
    endif ! end of if( NOTstellsym) ; 11 Aug 14;
  endif
  


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! if( lvol.gt.Nvol .and. Lconstraint.eq.-1 .and. Wcurent ) then ! 27 Feb 17;
  if( lvol.gt.Nvol                         .and. Wcurent ) then ! 27 Feb 17;
  !write(ounit,'("brcast : " 10x " : myid="i3" ; broadcasting : curtor="es13.5" ; curpol="es13.5" ;")') myid, curtor, curpol
   RlBCAST( curtor, 1, llmodnp )
   RlBCAST( curpol, 1, llmodnp )
  endif
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(brcast)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine brcast

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
