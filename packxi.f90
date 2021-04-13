!> \file packxi.f90
!> \brief Packs, and unpacks, geometrical degrees of freedom; and sets coordinate axis.

!> \ingroup grp_packing
!> \brief Packs, and unpacks, geometrical degrees of freedom; and sets coordinate axis.
!>
!> **geometrical degrees of freedom**
!> <ul>
!> <li> The geometrical degrees-of-freedom,
!>      namely the \f$R_{j,v}\f$ and \f$Z_{j,v}\f$ where \f$v\f$ labels the interface and \f$j\f$ labels the Fourier harmonic,
!>      must be "packxi", and "unpackxi", into a single vector, \f$\boldsymbol{\xi}\f$, so that standard numerical routines can be
!>      called to find solutions to force-balance, i.e. \f${\bf F}[\boldsymbol{\xi}]=0\f$. </li>
!> <li> A coordinate "pre-conditioning" factor is included: 
!>      \f{eqnarray}{ \boldsymbol{\xi}_k \equiv \frac{R_{j,v}}{\Psi_{j,v}},
!>      \f}
!>      where \f$\Psi_{j,v} \equiv\,\f$\c psifactor(j,v) , which is defined in global.f90 . </li>
!> </ul> 
!>
!> **coordinate axis**
!> <ul>
!> <li> The coordinate axis is not an independent degree-of-freedom of the geometry.
!>      It is constructed by extrapolating the geometry of the innermost interface down to a line. </li>
!> <li> Note that if the coordinate axis depends only on the geometry of the innermost interface
!>      then the block tridiagonal structure of the the force-derivative matrix is preserved. </li>
!> <li> Define the arc-length weighted averages,
!>      \f{eqnarray}{ R_0(\zeta) \equiv \frac{\int_{0}^{2\pi} R_1(\theta,\zeta) dl}{L(\zeta)}, \qquad
!>                    Z_0(\zeta) \equiv \frac{\int_{0}^{2\pi} Z_1(\theta,\zeta) dl}{L(\zeta)},
!>      \f}
!>      where \f$L(\zeta)\equiv \int_{0}^{2\pi} dl\f$ and \f$dl \equiv \sqrt{ \partial_\theta R_1(\theta,\zeta)^2 + \partial_\theta Z_1(\theta,\zeta)^2 } \, d\theta\f$. </li>
!> <li> Note that if \f$dl\f$ does not depend on \f$\theta\f$, i.e. if \f$\theta\f$ is the equal arc-length angle, then the expressions simplify. </li>
!> <li> Note that the geometry of the coordinate axis thus constructed only depends on the geometry of the innermost interface, by which I 
!>      mean that the geometry of the coordinate axis is independent of the angle parameterization. </li>
!> </ul>
!>
!> **some numerical comments**
!> <ul>
!> <li> First, the differential poloidal length, \f$dl \equiv \sqrt{ R_\theta^2 + Z_\theta^2 }\f$, is computed in real space using 
!>      an inverse FFT from the Fourier harmonics of \f$R\f$ and \f$Z\f$. </li>
!> <li> Second, the Fourier harmonics of the \f$dl\f$ are computed using an FFT.
!>      The integration over \f$\theta\f$ to construct \f$L\equiv \int dl\f$ is now trivial: just multiply the \f$m=0\f$ harmonics of \f$dl\f$ by \f$2\pi\f$.
!>      The \c ajk(1:mn) variable is used. </li>
!> <li> Next, the weighted \f$R\,dl\f$ and \f$Z\,dl\f$ are computed in real space, and the poloidal integral is similarly taken. </li>
!> <li> Lastly, the Fourier harmonics are constructed using an FFT after dividing in real space. </li>
!> </ul>
!>
!> @param[in] NGdof
!> @param position
!> @param[in] Mvol
!> @param[in] mn
!> @param iRbc
!> @param iZbs
!> @param iRbs
!> @param iZbc
!> @param packorunpack
!> @param[in] LComputeDerivatives
!> @param[in] LComputeAxis
subroutine packxi( NGdof, position, Mvol, mn, iRbc, iZbs, iRbs, iZbc, packorunpack, LComputeDerivatives, LComputeAxis )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero
  
  use numerical, only :
  
  use fileunits, only : ounit
  
  use inputlist, only : Wpackxi, Igeometry, Ntor, Nvol, Lfindzero
  
  use cputiming, only : Tpackxi
  
  use allglobal, only : ncpu, myid, cpus, im, in, &
                        YESstellsym, NOTstellsym, &
                        ajk, Nt, Nz, Ntz, iRij, iZij, tRij, tZij, &
                        ijreal, ijimag, jireal, jiimag, efmn, ofmn, cfmn, sfmn, evmn, odmn, comn, simn, &
                        psifactor, Rscale 
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS

  LOGICAL, intent(in)    :: LComputeDerivatives ! indicates whether derivatives are to be calculated;
  LOGICAL, intent(in)    :: LComputeAxis        ! if to recompute the axis

  INTEGER, intent(in)    :: NGdof, Mvol, mn
  REAL                   :: position(0:NGdof), iRbc(1:mn,0:Mvol), iZbs(1:mn,0:Mvol), iRbs(1:mn,0:Mvol), iZbc(1:mn,0:Mvol)
  CHARACTER              :: packorunpack
  
  INTEGER                :: lvol, jj, kk, irz, issym, idof, ifail, ivol
  
  BEGIN(packxi)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  idof = 0 ! initialize counter; 14 Jan 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  do lvol = 1, Mvol-1 ! loop over internal interfaces;
   
   do jj = 1, mn ! loop over Fourier harmonics;
    
    do irz = 0, 1 ! loop over R & Z;
     
     if( Igeometry.lt.3 .and. irz.eq.1 ) cycle ! no dependence on Z; 14 Jan 13;
     
     do issym = 0, 1 ! loop over even & odd;
      
      if( YESstellsym .and. issym.eq.1 ) cycle
      
      if( issym.eq.0 .and. irz.eq.1 .and. jj.eq.1 ) cycle ! no dependence on Zbs_{0,0}; 14 Jan 13;
      if( issym.eq.1 .and. irz.eq.0 .and. jj.eq.1 ) cycle ! no dependence on Rbs_{0,0}; 14 Jan 13;
      
      idof = idof + 1
      
#ifdef DEBUG
      FATAL( packxi, idof.le.0 .or. idof.gt.NGdof, out of bounds )
#endif
      
      select case( packorunpack )
       
      case( 'P' ) !   pack vector of unknowns;
       
       if( irz.eq.0 .and. issym.eq.0 ) position(idof) = iRbc(jj,lvol) / psifactor(jj,lvol)
       if( irz.eq.1 .and. issym.eq.0 ) position(idof) = iZbs(jj,lvol) / psifactor(jj,lvol)
       if( irz.eq.0 .and. issym.eq.1 ) position(idof) = iRbs(jj,lvol) / psifactor(jj,lvol)
       if( irz.eq.1 .and. issym.eq.1 ) position(idof) = iZbc(jj,lvol) / psifactor(jj,lvol)
       
      case( 'U' ) ! unpack vector of unknowns;
       
       if( irz.eq.0 .and. issym.eq.0 ) iRbc(jj,lvol) = position(idof) * psifactor(jj,lvol)
       if( irz.eq.1 .and. issym.eq.0 ) iZbs(jj,lvol) = position(idof) * psifactor(jj,lvol)
       if( irz.eq.0 .and. issym.eq.1 ) iRbs(jj,lvol) = position(idof) * psifactor(jj,lvol)
       if( irz.eq.1 .and. issym.eq.1 ) iZbc(jj,lvol) = position(idof) * psifactor(jj,lvol)
       
      end select
      
     enddo ! end of do issym;
     
    enddo ! end of do irz;
    
   enddo ! end of do jj;
   
  enddo ! end of do lvol;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  FATAL( packxi, idof.ne.NGdof, counting error )
#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( YESstellsym ) then ! iRbc(    ,0:Mvol) = zero
   ;                     ; iZbs(1   ,0:Mvol) = zero
   ;                     ; iRbs(1:mn,0:Mvol) = zero
   ;                     ; iZbc(1:mn,0:Mvol) = zero
  else                   ! iRbc(    ,0:Mvol) = zero
   ;                     ; iZbs(1   ,0:Mvol) = zero
   ;                     ; iRbs(1   ,0:Mvol) = zero
   ;                     ! iZbc(    ,0:Mvol) = zero
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  select case( packorunpack )

  case( 'P' )

  case( 'U' )

   ivol = 1 ! take care with ivol: this variable name might be a global variable, but here it is local; 19 Jul 16; 
 
   if( (Mvol .ne. 1) .and. (Lfindzero .ne. 0) ) then
    if (LComputeAxis) then
      WCALL( packxi, rzaxis, ( Mvol, mn, iRbc(1:mn,0:Mvol), iZbs(1:mn,0:Mvol), iRbs(1:mn,0:Mvol), iZbc(1:mn,0:Mvol), ivol, LComputeDerivatives ) ) ! set coordinate axis; 19 Jul 16; 
    endif
   endif

  end select
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(packxi)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine packxi

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
