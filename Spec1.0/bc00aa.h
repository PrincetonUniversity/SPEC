!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Broadcasts Beltrami field, profiles . . .

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine bc00aa( lvol )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero

  use numerical, only :

  use fileunits, only : ounit

  use inputlist, only : Wmacros, Wbc00aa, MNvol, Nvol, mu, Lrad, &
!                       Lwrpj, &
                        curtor, curpol, Lminimize, Lfindzero

  use cputiming, only : Tbc00aa

  use allglobal, only : myid, cpus, ncpu, dpflux, Ntz, mn, Mvol, &
                        mns, diota, &
                        Ate, Aze, Ato, Azo, &
                        Bemn, Bomn, Iomn, Iemn, Somn, Semn, Pomn, Pemn, &
                        ImagneticOK, &
                        dBBdRZ, dIIdRZ, dFFdRZ, dBBdmp, dmupfdx, lgeometricaldof, Lhessianallocated, &
                        lBBintegral, lABintegral, &
                        vvolume, &
!                       lmns, dlambda, &
                        NOTstellsym, Lplasmaregion, Lvacuumregion
  
  use pjhamilto, only : Itor, Gpol, spmn

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS

  INTEGER, intent(in) :: lvol
  
  INTEGER             :: llmodnp, io, iRZl, ii, ideriv, Nbc
  
  BEGIN(bc00aa)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! recall this routine is inside do vvol = 1, Nvol loop;
  
  FATALMESS(bc00aa, lvol.le.0   , error)
  FATALMESS(bc00aa, lvol.gt.Mvol, error)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  llmodnp = modulo(lvol-1,ncpu) ! identify which node contains data; this must be consistent with previous looping / parallelization;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( Lplasmaregion ) then

#ifdef DEBUG
   FATALMESS(bc00aa, lvol.gt.MNvol, mu is not defined)
#endif
   RlBCAST(mu(lvol),1,llmodnp)

#ifdef DEBUG
   FATALMESS(bc00aa, lvol.gt.Nvol, dpflux is not defined)
#endif
   RlBCAST(dpflux(lvol),1,llmodnp)

#ifdef DEBUG
   FATALMESS(bc00aa, lvol.gt.Mvol, volume etc is not defined)
   FATALMESS(bc00aa, .not.allocated(vvolume), error)
   FATALMESS(bc00aa, .not.allocated(lBBintegral), error)
   FATALMESS(bc00aa, .not.allocated(lABintegral), error)
#endif
   RlBCAST(vvolume(lvol),1,llmodnp)
   RlBCAST(lBBintegral(lvol),1,llmodnp)
   RlBCAST(lABintegral(lvol),1,llmodnp)
  else ! vacuum region; 19 Sep 13;
   RlBCAST(curtor,1,llmodnp)
   RlBCAST(curpol,1,llmodnp) !! probably does not need to be broadcast; this is set by coil geometry/currents and should be set earlier; 21 Apr 13;
  endif
  
#ifdef DEBUG
   FATALMESS(bc00aa, .not.allocated(diota), error)
#endif
 !RlBCAST(diota(-1:2,0:1,lvol),8,llmodnp) ! can be deleted; 15 Sep 15;
  RlBCAST(diota(0:1,-1:2,lvol),8,llmodnp)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! the derivatives of the lvol-th energy functional with respect to the inner & outer, interface Fourier harmonics, Rj & Zj;
  
  if( Lminimize.gt.0 ) then
  RlBCAST(dBBdRZ(lvol,0:1,1:lgeometricaldof),2*lgeometricaldof,llmodnp)
  RlBCAST(dIIdRZ(lvol,    1:lgeometricaldof),  lgeometricaldof,llmodnp)
  endif
  
  if( Lhessianallocated ) then

#ifdef DEBUG
   FATALMESS(bc00aa, .not.allocated(dFFdRZ), error)
#endif
  Nbc = lgeometricaldof*2*lgeometricaldof*2
  RlBCAST(dFFdRZ(1:lgeometricaldof,lvol,0:1,1:lgeometricaldof,0:1),Nbc,llmodnp)

#ifdef DEBUG
   FATALMESS(bc00aa, .not.allocated(dBBdmp), error)
#endif
  Nbc = lgeometricaldof*2*2                
  RlBCAST(dBBdmp(1:lgeometricaldof,lvol,0:1,1:2                  ),Nbc,llmodnp)

#ifdef DEBUG
   FATALMESS(bc00aa, .not.allocated(dmupfdx), error)
#endif
  RlBCAST(dmupfdx(lvol,1:2,1:lgeometricaldof,0:1),2*lgeometricaldof*2,llmodnp) ! why is this broadcast; 02 Sep 14;

  endif

! RlBCAST(   DenergyDrz(lvol,0:1,1:lgeometricaldof),2*lgeometricaldof,llmodnp)
! RlBCAST( DanalyticDrz(lvol,0:1,1:lgeometricaldof),2*lgeometricaldof,llmodnp)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! force-balance and angle-constraints;

#ifdef DEBUG
  FATALMESS(bc00aa, .not.allocated(ImagneticOK), error)
  FATALMESS(bc00aa, .not.allocated(Bemn), error)
  FATALMESS(bc00aa, .not.allocated(Bomn), error)
  FATALMESS(bc00aa, .not.allocated(Iomn), error)
  FATALMESS(bc00aa, .not.allocated(Iemn), error)
  FATALMESS(bc00aa, .not.allocated(Somn), error)
  FATALMESS(bc00aa, .not.allocated(Semn), error)
  FATALMESS(bc00aa, .not.allocated(Pomn), error)
  FATALMESS(bc00aa, .not.allocated(Pemn), error)
#endif

  LlBCAST(ImagneticOK(lvol),1,llmodnp)
  RlBCAST(Bemn(1:mn,lvol,0:1),2*mn,llmodnp) ! perhaps all these should be re-ordered; 18 Jul 14;
  RlBCAST(Bomn(1:mn,lvol,0:1),2*mn,llmodnp)
  RlBCAST(Iomn(1:mn,lvol    ),  mn,llmodnp)
  RlBCAST(Iemn(1:mn,lvol    ),  mn,llmodnp)
  RlBCAST(Somn(1:mn,lvol,0:1),2*mn,llmodnp)
  RlBCAST(Semn(1:mn,lvol,0:1),2*mn,llmodnp)

  RlBCAST(Pomn(1:mn,lvol,0:2),3*mn,llmodnp) ! 15 Sep 15;
  RlBCAST(Pemn(1:mn,lvol,0:2),3*mn,llmodnp)

! the derivatives of force-balance and angle-constraints with respect to interface Fourier harmonics;  

! if( Lhessianallocated ) then
   
!  ;RlBCAST( DforcebalDrz(lvol,0:1,1:lgeometricaldof,0:1,1:mn),2*lgeometricaldof*2* mn   ,llmodnp)
!  ;LlBCAST( LforcebalDrz(lvol    ,1:lgeometricaldof         ),  lgeometricaldof         ,llmodnp)
!   RlBCAST( DspectralDrz(lvol,    1:lgeometricaldof    ,2:mn),  lgeometricaldof  *(mn-1),llmodnp)
   
! endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG

  FATALMESS(bc00aa, lvol.gt.Mvol, Ate Azo Ato Azo are not defined)

  FATALMESS(bc00aa, .not.allocated(Ate), error)
  FATALMESS(bc00aa, .not.allocated(Aze), error)

  if( NOTstellsym ) then
  FATALMESS(bc00aa, .not.allocated(Ato), error)
  FATALMESS(bc00aa, .not.allocated(Azo), error)
  endif

#endif
  
  ideriv = 0

  do ii = 1, mn

#ifdef DEBUG
   FATALMESS(bc00aa, .not.allocated(Ate(lvol,ideriv,ii)%s), error)
   FATALMESS(bc00aa, .not.allocated(Aze(lvol,ideriv,ii)%s), error)
#endif
   
   RlBCAST( Ate(lvol,ideriv,ii)%s(0:Lrad(lvol)), Lrad(lvol)+1, llmodnp)
   RlBCAST( Aze(lvol,ideriv,ii)%s(0:Lrad(lvol)), Lrad(lvol)+1, llmodnp) ! this is not defined in vacuum region, but it should have been allocated; 24 Apr 13;

  enddo

  if( NOTstellsym ) then

  do ii = 1, mn 

#ifdef DEBUG
   FATALMESS(bc00aa, .not.allocated(Ato(lvol,ideriv,ii)%s), error)
   FATALMESS(bc00aa, .not.allocated(Azo(lvol,ideriv,ii)%s), error)
#endif

   RlBCAST( Ato(lvol,ideriv,ii)%s(0:Lrad(lvol)), Lrad(lvol)+1, llmodnp)
   RlBCAST( Azo(lvol,ideriv,ii)%s(0:Lrad(lvol)), Lrad(lvol)+1, llmodnp) ! this is not defined in vacuum region, but it should have been allocated; 24 Apr 13;

  enddo

  endif ! end of if( NOTstellsym) ; 11 Aug 14;

! FATALMESS(bc00aa, .not.allocated(dlambda), error)
! FATALMESS(bc00aa, lmns.le.0, allocation error in dlambda)

! RlBCAST( dlambda(1:lmns,-1:2,0:1,lvol), lmns*4*2, llmodnp)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(bc00aa)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine bc00aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
