!> \file
!> \brief Evaluates volume integrals, and their derivatives w.r.t. interface geometry, using "packed" format.

!> \brief Evaluates volume integrals, and their derivatives w.r.t. interface geometry, using "packed" format.
!> \ingroup grp_integrals
!>
!> @param[in] pNN
!> @param[in] xi
!> @param[out] Fxi
!> @param[out] DFxi
!> @param[in] Ldfjac
!> @param[in] iflag
subroutine df00ab( pNN , xi , Fxi , DFxi , Ldfjac , iflag )
  use mod_kinds, only: wp => dp
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  use constants, only : zero, half, one, two
  use numerical, only : small
  use fileunits, only : ounit
  use inputlist, only : Wdf00ab, Nvol, helicity
  use cputiming
  use allglobal, only : myid, cpus, MPI_COMM_SPEC, &
                        dMA, dMD, & ! energy and helicity matrices; 26 Feb 13;
!                       MBpsi, MEpsi, psiMCpsi, psiMFpsi, & ! pre-calculated matrix vector products; 26 Feb 13;
                        MBpsi,                            & ! pre-calculated matrix vector products; 26 Feb 13;
                        ivol


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


#ifdef OPENMP
  USE OMP_LIB
#endif
  use mpi
  implicit none
  integer   :: ierr, astat, ios, nthreads, ithread
  real(wp)      :: cput, cpui, cpuo=0 ! cpu time; cpu initial; cpu old; 31 Jan 13;


  integer, intent(in)  :: pNN, Ldfjac
  integer, intent(in)  :: iflag
  real(wp)   , intent(in)  :: xi(0:pNN-1)
  real(wp)   , intent(out) :: Fxi(0:pNN-1), DFxi(0:Ldfjac-1,0:pNN-1)

  integer             :: NN
  real(wp)                :: lmu ! , dpsi(1:2)


  cpui = MPI_WTIME()
  cpuo = cpui
#ifdef OPENMP
  nthreads = omp_get_max_threads()
#else
  nthreads = 1
#endif


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef DEBUG


   if( ivol.lt.1 .or. ivol.gt.Nvol ) then
     write(6,'("df00ab :      fatal : myid=",i3," ; ivol.lt.1 .or. ivol.gt.Nvol ; ivol invalid ;")') myid
     call MPI_ABORT( MPI_COMM_SPEC, 1, ierr )
     stop "df00ab : ivol.lt.1 .or. ivol.gt.Nvol : ivol invalid  ;"
    endif
 ! 26 Feb 13;

#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  NN = pNN-1 ; lmu = xi(0)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  select case( iflag ) ! F = B.B - mu ( A.B - helicity );

  case( 1 ) ! return d F;

! !Fxi(   0) = - half * sum( xi(1:NN) * matmul( dMD(1:NN,1:NN), xi(1:NN) ) ) - sum( xi(1:NN) * MEpsi(1:NN) ) - psiMFpsi + helicity(ivol)
   Fxi(   0) = - half * sum( xi(1:NN) * matmul( dMD(1:NN,1:NN), xi(1:NN) ) )                                            + helicity(ivol)

! !Fxi(1:NN) = matmul( dMA(1:NN,1:NN), xi(1:NN)  ) + MBpsi(1:NN) - lmu * ( matmul( dMD(1:NN,1:NN), xi(1:NN)  ) + MEpsi(1:NN) )
   Fxi(1:NN) = matmul( dMA(1:NN,1:NN), xi(1:NN)  ) + MBpsi(1:NN) - lmu * ( matmul( dMD(1:NN,1:NN), xi(1:NN)  )               )

  case( 2 ) ! return ddF;

   DFxi(   0,   0) = zero

   DFxi(   0,1:NN) = - matmul( dMD(1:NN,1:NN), xi(1:NN)  ) ! - MEpsi(1:NN)
   DFxi(1:NN,   0) = - matmul( dMD(1:NN,1:NN), xi(1:NN)  ) ! - MEpsi(1:NN)

   DFxi(1:NN,1:NN) = dMA(1:NN,1:NN) - lmu * dMD(1:NN,1:NN)

  case default


   if( .true. ) then
     write(6,'("df00ab :      fatal : myid=",i3," ; .true. ; supplied value of iflag is not supported ;")') myid
     call MPI_ABORT( MPI_COMM_SPEC, 1, ierr )
     stop "df00ab : .true. : supplied value of iflag is not supported  ;"
    endif


  end select

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


9999 continue
  cput = MPI_WTIME()
  Tdf00ab = Tdf00ab + ( cput-cpuo )
  return


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine df00ab

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
