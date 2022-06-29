!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title ! Evaluates volume integrals, and their derivatives w.r.t. interface geometry, using &ldquo;packed&rdquo; format.

!latex \briefly{briefly}
!latex \calledby{\link{}}
!latex \calls{\link{}}
!latex \tableofcontents

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine df00ab( pNN , xi , Fxi , DFxi , Ldfjac , iflag )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half, one, two
  use numerical, only : small
  use fileunits, only : ounit
  use inputlist, only : Wdf00ab, Nvol, helicity
  use cputiming
  use allglobal, only : myid, cpus, &
                        dMA, dMD, & ! energy and helicity matrices; 26 Feb 13;
!                       MBpsi, MEpsi, psiMCpsi, psiMFpsi, & ! pre-calculated matrix vector products; 26 Feb 13;
                        MBpsi,                            & ! pre-calculated matrix vector products; 26 Feb 13;
                        ivol

  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in) :: pNN, Ldfjac
  INTEGER             :: iflag
  REAL   , intent(in) :: xi(0:pNN-1)
  REAL                :: Fxi(0:pNN-1), DFxi(0:Ldfjac-1,0:pNN-1)
  
  INTEGER             :: NN
  REAL                :: lmu ! , dpsi(1:2)
  
  BEGIN(df00ab)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! write(ounit,'("df00ab : ", 10x ," : 1.1000 xi =",99(es23.15,","))') xi(0)
  
#ifdef DEBUG

  FATAL(df00ab, ivol.lt.1 .or. ivol.gt.Nvol, ivol invalid ) ! 26 Feb 13;

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

   FATAL(df00ab, .true., supplied value of iflag is not supported )

  end select
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! write(ounit,'("df00ab : ", 10x ," : 1.9000 xi =",99(es23.15,","))') xi(0)

  RETURN(df00ab)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine df00ab

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
