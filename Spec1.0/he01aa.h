!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Compute eigenvalues and eigenvectors of derivative matrix.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine he01aa( Ngeometricaldof, position, Mvol, mn, lgeometricaldof )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  use constants, only : zero, half, one, two, ten

  use numerical, only : sqrtmachprec, small, vsmall

  use fileunits, only : ounit, lunit, hunit

  use inputlist, only : Wmacros, Whe01aa, ext, Igeometry, Nvol, tflux, pflux, helicity, mu, Energy, Lfreebound, &
                        LHevalues, LHevectors, &
                        Lperturbed, dpp, dqq, &
                        Lcheck

  use cputiming, only : The01aa

  use allglobal, only : ncpu, myid, cpus, YESstellsym, im, in, halfmm, &
                        iRbc, iZbs, iRbs, iZbc, &
                        dRbc, dZbs, dRbs, dZbc, &
                        lBBintegral, dBBdRZ, &
                        NOTstellsym, YESstellsym, &
!                       lgeometricaldof, mlocaldof, &
!                       Nt, Nz, Ntz, gij, sg, dxsg, dxgij, &
!                       cfmn, sfmn, efmn, ofmn, ijreal, ijimag, trigwk, trigm, trign, isr, &
!                       lBBintegral, dBBintegral, dAd, DenergyDrz, DanalyticDrz, &
!                       lgeometricaldof, DforcebalDrz, DspectralDrz, LforcebalDrz, &
                        dFFdRZ, dBBdmp, dmupfdx, hessian, dessian, Lhessianallocated, psifactor
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in) :: Ngeometricaldof, Mvol, mn, lgeometricaldof
  REAL                :: position(0:Ngeometricaldof) ! internal geometrical degrees of freedom;
  
  LOGICAL             :: LComputeDerivatives

  REAL                :: force(0:Ngeometricaldof), gradient(0:Ngeometricaldof)

  REAL                :: xx(0:Ngeometricaldof,-2:2), ff(0:Ngeometricaldof,-2:2), df(1:Ngeometricaldof)!, deriv

  INTEGER             :: vvol, idof, ii, mi, ni, irz, issym, isymdiff, lvol, ieval(1:1), igdof, ifd, if03aaf
  REAL                :: oldEnergy(-2:2), error, cpul

  REAL                :: oldBB(1:Mvol,-2:2), oBBdRZ(1:Mvol,0:1,1:lgeometricaldof), ohessian(1:Ngeometricaldof,1:Ngeometricaldof)

  REAL                :: oRbc(1:mn,0:Mvol), oZbs(1:mn,0:Mvol), oRbs(1:mn,0:Mvol), oZbc(1:mn,0:Mvol), determinant

  CHARACTER           :: packorunpack

  CHARACTER           :: svol*3

!  LOGICAL            :: Lderiv                                                                         ! for parallel / series construction of Hessian;
!  INTEGER            :: lvol, jvol, ivol, innout, imn, irz, jmn, jrz, tdoc, tdof, ilocaldof, jlocaldof ! for parallel / series construction of Hessian;

  INTEGER             :: tdof, tdoc, jvol, jj, jrz, jssym
  
!  REAL, allocatable  :: oldBB(:) ! original integral;
  
  INTEGER             :: Lwork, LDA, Ldvi, Ldvr, if02ebf
  REAL                :: evalr(1:Ngeometricaldof), evali(1:Ngeometricaldof)
  REAL                :: evecr(1:Ngeometricaldof,1:Ngeometricaldof), eveci(1:Ngeometricaldof,1:Ngeometricaldof)
  REAL                :: work(1:4*Ngeometricaldof) ! for construction of evalues/evectors;
  CHARACTER           :: JOB
  
  INTEGER             :: iev, jev, M1, M2, irank(1:Ngeometricaldof), im01daf ! for construction of evalues/evectors;
  CHARACTER           :: order

  REAL                :: dRZ != 1.0e-03

  REAL                :: lmu(1:Mvol), lpflux(1:Mvol), lhelicity(1:Mvol) ! original profiles; 20 Jun 14;

  INTEGER             :: IA, IAA, if04atf
  REAL                :: perturbation(1:lgeometricaldof)!, maxsoln
  REAL                :: rhs(1:Ngeometricaldof), soln(0:Ngeometricaldof)
  REAL                :: rworka(1:Ngeometricaldof), rworkb(1:Ngeometricaldof), AA(1:Ngeometricaldof,1:Ngeometricaldof)
  
  BEGIN(he01aa)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  lmu(1:Nvol) = mu(1:Nvol) ; lpflux(1:Nvol) = pflux(1:Nvol) ; lhelicity(1:Nvol) = helicity(1:Nvol) ! save original profile information; 20 Jun 14;
   
  oldEnergy(0) = Energy ! Energy was calculated in fc02aa; 26 Feb 13;

  xx(0,-2:2)= zero

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef MINIMIZE
  
  oldBB(1:Mvol,0) = lBBintegral(1:Mvol)
  
  FATALMESS(he01aa, .not.allocated(dBBdRZ), need to revise logic in al00aa where dBBdRZ is allocated )
  
  oBBdRZ(1:Mvol,0:1,1:lgeometricaldof) = dBBdRZ(1:Mvol,0:1,1:lgeometricaldof)
  
  oRbc(1:mn,0:Mvol) = iRbc(1:mn,0:Mvol)
  oZbs(1:mn,0:Mvol) = iZbs(1:mn,0:Mvol)
  oRbs(1:mn,0:Mvol) = iRbs(1:mn,0:Mvol)
  oZbc(1:mn,0:Mvol) = iZbc(1:mn,0:Mvol)
  
  FATALMESS(he01aa, Lfreebound.gt.0, this routine needs attention )
  
  do vvol = 1, Mvol-1 ! loop over volumes; 26 Feb 13;

   write(ounit,'("he01aa : ", 10x ," : vvol=",i3," ;")') vvol
   
   idof = 0
   
   do ii = 1, mn ; mi = im(ii) ; ni = in(ii) ! loop over Fourier harmonics; 26 Feb 13;
    
    do irz = 0, 1 ! loop over coordinate functions; 26 Feb 13;
     
     if( irz.eq.1 .and. Igeometry.lt.3 ) cycle ! no dependence on Z; 26 Feb 13;
     
     do issym = 0, 1 ! loop over stellarator and non-stellarator symmetric terms; 26 Feb 13;
      
      if( issym.eq.1 .and. YESstellsym ) cycle ! no dependence on non-stellarator symmetric harmonics; 26 Feb 13;
      
      if( ii.eq.1 .and. irz.eq.1 .and. issym.eq.0 ) cycle ! no dependence on Zbs_{m=0,n=0}; 26 Feb 13;
      if( ii.eq.1 .and. irz.eq.0 .and. issym.eq.1 ) cycle ! no dependence on Rbs_{m=0,n=0}; 26 Feb 13;
      
      iRbc(1:mn,0:Mvol) = oRbc(1:mn,0:Mvol)
      iZbs(1:mn,0:Mvol) = oZbs(1:mn,0:Mvol)
      iRbs(1:mn,0:Mvol) = oRbs(1:mn,0:Mvol)
      iZbc(1:mn,0:Mvol) = oZbc(1:mn,0:Mvol)

      idof = idof + 1 ! labels degree of freedom; 26 Feb 13;
      
      do isymdiff = -2, 2 ! symmetric fourth-order, finite-difference used to approximate derivatives; 26 Feb 13;

       if( isymdiff.eq.0 ) cycle
       
       dRZ = 1.0E-03

       if( issym.eq.0 ) then ! consider     stellarator symmetric harmonics; 26 Feb 13;
        if( irz.eq.0 ) iRbc(ii,vvol) = oRbc(ii,vvol) + dRZ * isymdiff
        if( irz.eq.1 ) iZbs(ii,vvol) = oZbs(ii,vvol) + dRZ * isymdiff
       else                  ! consider non-stellarator symmetric harmonics; 26 Feb 13;
        if( irz.eq.0 ) iRbs(ii,vvol) = oRbs(ii,vvol) + dRZ * isymdiff
        if( irz.eq.1 ) iZbc(ii,vvol) = oZbc(ii,vvol) + dRZ * isymdiff
       endif
       
       packorunpack = 'P' !; position(0) = zero ! this is not used; 11 Aug 14;
       WCALL(he01aa,gf00aa,( Ngeometricaldof, position(0:Ngeometricaldof), Mvol, mn, &
    iRbc(1:mn,0:Mvol), iZbs(1:mn,0:Mvol), iRbs(1:mn,0:Mvol), iZbc(1:mn,0:Mvol), packorunpack ))
       
       LComputeDerivatives = .false. !; position(0) = zero ! this is not used; 11 Aug 14;
       WCALL(he01aa,fc02aa,( Ngeometricaldof, position(0:Ngeometricaldof), gradient(0:Ngeometricaldof), LComputeDerivatives )) ! re-calculate Beltrami fields;
       
       oldBB(1:Mvol,isymdiff) = lBBintegral(1:Mvol)

       oldEnergy(isymdiff) = Energy
       
      enddo ! end of do isymdiff; 26 Feb 13;

      oldBB(1:Mvol,0) = ( - 1 * oldBB(1:Mvol,2) + 8 * oldBB(1:Mvol,1) - 8 * oldBB(1:Mvol,-1) + 1 * oldBB(1:Mvol,-2) ) / ( 12 * dRZ ) ! 4th order estimate;

      oldEnergy(0)    = ( - 1 * oldEnergy(2)    + 8 * oldEnergy(1)    - 8 * oldEnergy(-1)    + 1 * oldEnergy(-2)    ) / ( 12 * dRZ )

      cput = GETTIME
      write(ounit,1000) cput-cpus                           !12345678901234567
      write(ounit,1000) cput-cpus, myid, vvol, irz, mi, ni, "finite-difference", oldBB(vvol:vvol+1,0)
      write(ounit,1000) cput-cpus, myid, vvol, irz, mi, ni, "analytic         ", (/ oBBdRZ(vvol,1,idof), oBBdRZ(vvol+1,0,idof) /) / psifactor(ii,vvol)

      write(ounit,1001) cput-cpus, myid, vvol, irz, mi, ni, oldEnergy(0)
      write(ounit,1001) cput-cpus, myid, vvol, irz, mi, ni, ( oBBdRZ(vvol,1,idof) + oBBdRZ(vvol+1,0,idof) ) / psifactor(ii,vvol) ! ENERGY GRADIENT;
      FATALMESS(he01aa, Igeometry.eq.1, Cartesian geometry does not need regularization factor)

1000  format("he01aa : ",f10.2," : ":"myid=",i3," ; ":"vvol=",i3," ; ":"irz="i2" ; (",i3," ,",i3," ) ; "a17" ["es15.7","es15.7" ]")
1001  format("he01aa : ",f10.2," : ":"myid=",i3," ; ":"vvol=",i3," ; ":"irz="i2" ; (",i3," ,",i3," ) ;   "es15.7" ; ")
      
     enddo ! end of do issym; 26 Feb 13;
     
    enddo ! end of do irz; 26 Feb 13;
     
   enddo ! end of do ii; 26 Feb 13;
   
  enddo ! end of do vvol; 26 Feb 13;
  
  iRbc(1:mn,0:Mvol) = oRbc(1:mn,0:Mvol)
  iZbs(1:mn,0:Mvol) = oZbs(1:mn,0:Mvol)
  iRbs(1:mn,0:Mvol) = oRbs(1:mn,0:Mvol)
  iZbc(1:mn,0:Mvol) = oZbc(1:mn,0:Mvol)
  
  packorunpack = 'P' !; position(0) = zero ! this is not used; 11 Aug 14;
  WCALL(he01aa,gf00aa,( Ngeometricaldof, position(0:Ngeometricaldof), Mvol, mn, &
iRbc(1:mn,0:Mvol), iZbs(1:mn,0:Mvol), iRbs(1:mn,0:Mvol), iZbc(1:mn,0:Mvol), packorunpack ))
  
#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RALLOCATE(dFFdRZ,(1:lgeometricaldof,1:Mvol,0:1,1:lgeometricaldof,0:1))
  RALLOCATE(dBBdmp,(1:lgeometricaldof,1:Mvol,0:1,1:2                  ))
  RALLOCATE(dmupfdx,(1:Mvol,1:2,1:lgeometricaldof,0:1))

  RALLOCATE(hessian,(1:Ngeometricaldof,1:Ngeometricaldof))
  RALLOCATE(dessian,(1:Ngeometricaldof,1:lgeometricaldof)) ! part of hessian that depends on boundary variations; 18 Dec 14;
  Lhessianallocated = .true.
  
  LComputeDerivatives = .true. !; position(0) = zero ! this is not used; 11 Aug 14;
  WCALL(he01aa,fc02aa,( Ngeometricaldof, position(0:Ngeometricaldof), force(0:Ngeometricaldof), LComputeDerivatives )) ! calculate force-imbalance & hessian;
  
  ohessian(1:Ngeometricaldof,1:Ngeometricaldof) = hessian(1:Ngeometricaldof,1:Ngeometricaldof) ! internal copy; 22 Apr 15;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \end{enumerate} \subsubsection{construction of Hessian matrix} \begin{enumerate}

!latex \item The routine \verb+fc02aa+ is used to compute the derivatives, with respect to interface geometry,
!latex       of the force imbalance harmonics, $[[p+B^2/2]]_{j}$, which may be considered to be the ``physical'' constraints,
!latex       and if \verb+Igeometry.eq.3+ then also the derivatives of the ``artificial'' spectral constraints, $I_j \equiv (R_\t X + Z_\t Y)_j$.

!latex \item The input variable \verb+Lconstraint+ determines whether the derivatives are taken at 
!latex       \begin{enumerate}
!latex       \item \verb+Lconstraint.eq.0+ : fixed $\mu$ and $\delta \psi_p$, 
!latex       \item \verb+Lconstraint.eq.1+ : fixed interface transform, or 
!latex       \item \verb+Lconstraint.eq.2+ : fixed helicity, ${\cal K}$, and $\delta \psi_p$.
!latex       \end{enumerate}
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( Lcheck.eq.5 ) then ! check construction of Hessian; 01 Jul 14;
   
   xx(0,-2:2)= zero ; dRZ = 1.0E-04

   write(svol,'(i3.3)')myid
   open(lunit+myid,file="."//trim(ext)//".hessian."//svol,status="unknown")
   
!  lmu(1:Nvol) = mu(1:Nvol) ; lpflux(1:Nvol) = pflux(1:Nvol) ; lhelicity(1:Nvol) = helicity(1:Nvol) ! save original profile information; 20 Jun 14;
   
   tdof = 0 ! geometrical degree-of-freedom counter; labels derivative;
   
   do vvol = 1, Mvol-1 ! loop over internal interfaces;
    
    cput = GETTIME
    write(lunit+myid,'("he01aa : ", 10x ," : ")')
    write(lunit+myid,'("he01aa : ",f10.2," : myid=",i3," ; vvol=",i3," ; dRZ=",es9.1," ;")') cput-cpus, myid, vvol, dRZ
    write(lunit+myid,'("he01aa : ", 10x ," : ")')
    
    do ii = 1, mn ! loop over Fourier harmonics degrees of freedom;
     
     do irz = 0, 1 ! loop over R,Z degrees of freedom;
      
      if( irz.eq.1 .and. Igeometry.lt.3 ) cycle ! no dependence on Z; 26 Feb 13;
      
      do issym = 0, 1
       
       if( issym.eq.1 .and. YESstellsym ) cycle ! no dependence on non-stellarator symmetric harmonics; 26 Feb 13;
       
       if( ii.eq.1 .and. irz.eq.1 .and. issym.eq.0 ) cycle ! this is not a degree of freedom; no dependence on Zbs_{m,n} for m=0, n=0;
       if( ii.eq.1 .and. irz.eq.0 .and. issym.eq.1 ) cycle ! this is not a degree of freedom; no dependence on Rbs_{m,n} for m=0, n=0;
       
       tdof = tdof + 1
       
       do isymdiff = -2, 2
        
        if( isymdiff.eq.0 ) cycle

        xx(1:Ngeometricaldof,isymdiff) = position(1:Ngeometricaldof) ! reset geometry to original;
        
        xx(tdof,isymdiff) = position(tdof) + isymdiff * dRZ ! perturb appropriate geometric harmonic;
        
        LComputeDerivatives = .false.
        WCALL(he01aa,fc02aa,( Ngeometricaldof, xx(0:Ngeometricaldof,isymdiff), ff(0:Ngeometricaldof,isymdiff), LComputeDerivatives )) ! force-imbalance;
        
       enddo ! end of do isymdiff; 20 Jun 14;
       
       df(1:Ngeometricaldof) = ( - ff(1:Ngeometricaldof,+2) + 8 * ff(1:Ngeometricaldof,+1) - 8 * ff(1:Ngeometricaldof,-1) + ff(1:Ngeometricaldof,-2) ) / &
                               ( 12 * dRZ )
       
       tdoc = 0
       
       do jvol = 1, Mvol-1
        
        do jj = 1, mn
         
         do jrz = 0, 1 ! loop over R,Z degrees of freedom;
          
          if( jrz.eq.1 .and. Igeometry.lt.3 ) cycle ! no dependence on Z; 26 Feb 13;
          
          do jssym = 0, 1
           
           if( jssym.eq.1 .and. YESstellsym ) cycle ! no dependence on non-stellarator symmetric harmonics; 26 Feb 13;
           
           if( jj.eq.1 .and. jrz.eq.1 .and. jssym.eq.0 ) cycle ! this is not a constraint; I_{m,n} for m=0, n=0 is irrelevant;
           if( jj.eq.1 .and. jrz.eq.0 .and. jssym.eq.1 ) cycle ! this is not a constraint; I_{m,n} for m=0, n=0 is irrelevant;
           
           tdoc = tdoc + 1
           
           FATALMESS(he01aa, tdoc.lt.              1, needs attention)
           FATALMESS(he01aa, tdoc.gt.Ngeometricaldof, needs attention)
           
           cput = GETTIME
           
           error = abs( df(tdoc)-hessian(tdoc,tdof) )

           if( abs(hessian(tdoc,tdof)).gt.1.0e-05 .or. abs(df(tdoc)).gt.1.0e-05 .or. error.gt.dRZ .or. .true. ) then ! write to screen; 20 Jan 15;
            write(ounit     ,1001) cput-cpus, myid, &
          vvol, im(ii), in(ii), irz, issym, tdof, &
          jvol, im(jj), in(jj), jrz, jssym, tdoc, &
          df(tdoc), tdoc, tdof, hessian(tdoc,tdof), error
           endif

           if( abs(hessian(tdoc,tdof)).gt.1.0e-05 .or. abs(df(tdoc)).gt.1.0e-05 .or. error.gt.dRZ             ) then ! write to file; 20 Jan 15;
            write(lunit+myid,1001) cput-cpus, myid, &
          vvol, im(ii), in(ii), irz, issym, tdof, &
          jvol, im(jj), in(jj), jrz, jssym, tdoc, &
          df(tdoc), tdoc, tdof, hessian(tdoc,tdof), error
           endif
           
          enddo ! end of do jssym; 19 Sep 13;
          
         enddo ! end of do irz;
         
        enddo ! end of do jj;
        
       enddo ! end of do jvol;
       
1001   format("he01aa : ",f10.2," : myid=",i3, &
     " ; vvol=",i3," ; ("i4" ,"i4" ); irz="i2" ; issym="i2" ; tdof="i6 &
     " ; jvol="i4" ; ("i4" ,"i4" ); jrz="i2" ; jssym="i2" ; tdoc="i6 &
     " ; fd-estimate="es13.5" & hessian("2i6")="es13.5" ;":" err="es13.5" ;":,f12.4" ;")
       
      enddo ! end of do issym;
      
     enddo ! end of do irz;
     
    enddo ! end of do ii;
    
   enddo ! end of do vvol;
   
   packorunpack = 'U' !; position(0) = zero ! this is not used; 11 Aug 14;
   WCALL(he01aa,gf00aa,( Ngeometricaldof, position(0:Ngeometricaldof), Mvol, mn, &
iRbc(1:mn,0:Mvol), iZbs(1:mn,0:Mvol), iRbs(1:mn,0:Mvol), iZbc(1:mn,0:Mvol), packorunpack ))
   
   mu(1:Nvol) = lmu(1:Nvol) ; pflux(1:Nvol) = lpflux(1:Nvol) ; helicity(1:Nvol) = lhelicity(1:Nvol)
   
!  xx(1:Ngeometricaldof,0) = position(1:Ngeometricaldof) ! reset geometry to original;
!  
!  LComputeDerivatives = .false.
!  WCALL(he01aa,fc02aa,( Ngeometricaldof, xx(1:Ngeometricaldof,0), ff(0:Ngeometricaldof,0), LComputeDerivatives )) ! calculate the force-imbalance;

   close(lunit+myid)
   
  endif ! end of if( Lcheck.eq.5 ) ; 01 Jul 14;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \end{enumerate} \subsubsection{construction of eigenvalues and eigenvectors} \begin{enumerate}
  
!latex \item If \verb+LHevalues.eq.T+ then the eigenvalues of the Hessian are computed using the NAG routine \verb+F02EBF+.
!latex \item If \verb+LHevectors.eq.T+ then the eigenvalues {\em and} the eigenvectors of the Hessian are computed.
  
!latex \item Note that if \verb+Igeometry.ge.4+, then the derivative-matrix also contains information regarding how the ``artificial'' spectral constraints
!latex       vary with geometry; so, the eigenvalues and eigenvectors are not purely ``physical''.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! if( myid.eq.0 .and. ( LHevalues .or. LHevectors ) ) then ! the call to fc02aa below requires all cpus; 04 Dec 14;
  if(                 ( LHevalues .or. LHevectors ) ) then
   
   if( myid.eq.0 ) then ; cput = GETTIME ; write(ounit,'("he01aa : ",f10.2," : LHevalues="L2" , LHevectors="L2" ;")')cput-cpus, LHevalues, LHevectors
   endif
   
   evalr(1:Ngeometricaldof) = zero ; evecr(1:Ngeometricaldof,1:Ngeometricaldof) = zero
   evali(1:Ngeometricaldof) = zero ; eveci(1:Ngeometricaldof,1:Ngeometricaldof) = zero
   
   if( LHevectors ) then ; JOB='V' ; Ldvr = Ngeometricaldof ; Ldvi = Ngeometricaldof
   else                  ; JOB='N' ; Ldvr =  1              ; Ldvi =  1              ! provide dummy values when eigenvectors are not required; 04 Dec 14;
   endif
   
   cpul = GETTIME
   
   if02ebf = 1 ; LDA = Ngeometricaldof ; Lwork = 4*Ngeometricaldof
   
   hessian(1:Ngeometricaldof,1:Ngeometricaldof) = ohessian(1:Ngeometricaldof,1:Ngeometricaldof)
   
#ifdef NAG18
   call F02EBF( JOB, Ngeometricaldof, hessian(1:LDA,1:Ngeometricaldof), LDA, evalr(1:Ngeometricaldof), evali(1:Ngeometricaldof), &
                evecr(1:Ldvr,1:Ngeometricaldof), Ldvr, eveci(1:Ldvi,1:Ngeometricaldof), Ldvi, work(1:Lwork), Lwork, if02ebf )
#else
   FATALMESS(global, .true., eigenvalue solver needs updating to F08NAF)
#endif
   
   if( myid.eq.0 ) then 
   cput = GETTIME
   select case( if02ebf )
   case( 0 )    ; write(ounit,'("he01aa : ",f10.2," : computed evalues ; if02ebf="i2" ; success ;       time="f10.2"s ;")') cput-cpus, if02ebf, cput-cpul
   case( 1 )    ; write(ounit,'("he01aa : ",f10.2," : computed evalues ; if02ebf="i2" ; input error ;   time="f10.2"s ;")') cput-cpus, if02ebf, cput-cpul
   case( 2 )    ; write(ounit,'("he01aa : ",f10.2," : computed evalues ; if02ebf="i2" ; QR failed ;     time="f10.2"s ;")') cput-cpus, if02ebf, cput-cpul
   case default ; write(ounit,'("he01aa : ",f10.2," : computed evalues ; if02ebf="i2" ; illegal ifail ; time="f10.2"s ;")') cput-cpus, if02ebf, cput-cpul
   end select
  endif

! DEALLOCATE(DforcebalDrz) ! probably redundant; 04 Dec 14;
! DEALLOCATE(LforcebalDrz)
! if( Igeometry.ge.4 ) then
! DEALLOCATE(DspectralDrz)
! endif
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
!   ieval(1:1) = minloc( evalr(1:Ngeometricaldof) ) ! 04 Dec 14;
!   if( myid.eq.0) then ! screen output; 04 Dec 14;
!   write(ounit,'("he01aa : " 10x " : evalr("i4") ="es13.5" ="es13.5" ;")') ieval(1), evalr(ieval(1)), minval(evalr(1:Ngeometricaldof))
!   endif
!   
!   do ifd = -1, -6, -1 ; dRZ = ten**ifd
!    
!    xx(1:Ngeometricaldof,0) = position(1:Ngeometricaldof) + evecr(1:Ngeometricaldof,ieval(1)) * dRZ ! perturb in direction of eigenvector; 04 Dec 14;
!    
!    LComputeDerivatives = .false.
!    WCALL(he01aa,fc02aa,( Ngeometricaldof, xx(0:Ngeometricaldof,0), ff(0:Ngeometricaldof,0), LComputeDerivatives )) ! calculate the force-imbalance;
!    
!    if( myid.eq.0 ) then ! screen output; 04 Dec 14;
!     write(ounit,'("he01aa : " 10x " : Energy(old)=", es23.15," ;")') oldEnergy(0)
!     write(ounit,'("he01aa : " 10x " : Energy(new)=", es23.15," ;")')    Energy
!     write(ounit,'("he01aa : " 10x " :             ",2es23.15," ;")')    Energy - oldEnergy(0), ( Energy-oldEnergy(0) ) / dRZ**2
!     write(ounit,'("he01aa : " 10x " :             ",2es23.15," ;")')    evalr(ieval(1)) * dRZ**2, evalr(ieval(1))
!    !write(ounit,'("he01aa : " 10x " : |evector|  =",es23.15," ;")')    sum(evecr(1:Ngeometricaldof,ieval(1))*evecr(1:Ngeometricaldof,ieval(1)))
!    !write(ounit,'("he01aa : ", 10x  " : "999(" ("i3","i3")":))') (/ ( im(ii), in(ii), ii = 1, mn ) /)
!    !do lvol = 1, Mvol-1
!    ! write(ounit,'("he01aa : ",i10   " : "999es10.2         )') lvol, evecr(1+mn*(lvol-1):mn+mn*(lvol-1),ieval(1))
!    !enddo
!    !write(ounit,'("he01aa : " 10x " : ")')
!    !do lvol = 1, Mvol-1
!    ! write(ounit,'("he01aa : ",i10   " : "999es10.2         )') lvol, ff(1+mn*(lvol-1):mn+mn*(lvol-1),0) / dRZ / evalr(ieval(1))
!    !enddo
!    endif ! end of if( myid.eq.0 ) ; 04 Dec 14;
!
!   enddo ! end of do ifd; 04 Dec 14;
!
!  !do igdof = 1, Ngeometricaldof ; write(ounit,'("he01aa : " 10x " : "f15.11" ="f15.11" ;")') evecr(igdof,ieval(1)), ff(igdof,0) / dRZ / evalr(ieval(1))
!  !enddo
!   
!   packorunpack = 'U' !; position(0) = zero ! this is not used; 11 Aug 14; ! reset geometry (Rbc, Zbs, etc.) to original values; 04 Dec 14;
!   WCALL(he01aa,gf00aa,( Ngeometricaldof, position(0:Ngeometricaldof), Mvol, mn, &
!   iRbc(1:mn,0:Mvol), iZbs(1:mn,0:Mvol), iRbs(1:mn,0:Mvol), iZbc(1:mn,0:Mvol), packorunpack ))
!   
!   mu(1:Nvol) = lmu(1:Nvol) ; pflux(1:Nvol) = lpflux(1:Nvol) ; helicity(1:Nvol) = lhelicity(1:Nvol) ! reset profiles to original values; 04 Dec 14;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   do lvol = 1, Mvol-1
    do ii = 1, mn ; evecr(ii+(lvol-1)*mn,1:Ngeometricaldof) = evecr(ii+(lvol-1)*mn,1:Ngeometricaldof) * psifactor(ii,lvol) ! geometrical regularization;
    enddo
   enddo
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   if( myid.eq.0 ) then ! screen output; 04 Dec 14;
    
    if( LHevectors ) then
     
     do iev = 1, Ngeometricaldof ! loop over all eigenvalues; 04 Dec 14;
      if( evalr(iev).lt.zero ) then ! only show unstable eigenvalues; 04 Dec 14;
       write(ounit,'("he01aa : ",f10.2," : evalr="es13.5" ; ")') cput-cpus, evalr(iev)
       write(ounit,'("he01aa : ",es10.3," : "999(" ("i3","i3")":))') evalr(iev), (/ ( im(ii), in(ii), ii = 1, mn ) /)
       do lvol = 1, Mvol-1 ; write(ounit,'("he01aa : ",i10   " : "999es10.2)') lvol, evecr(1+mn*(lvol-1):mn+mn*(lvol-1),iev)
       enddo
      endif
     enddo
     
    else ! matches if( LHvectors) ; 04 Dec 14;
     
     do iev = 1, Ngeometricaldof
     !if( evalr(iev).lt.zero ) then ! 22 Apr 15;
       write(ounit,'("he01aa : ",f10.2," : ",i6," : evalue=(",es23.15," ,",es23.15," ) ;")') cput-cpus, iev, evalr(iev), evali(iev) ! 04 Dec 14;
     !endif ! 22 Apr 15;
     enddo
      
    endif ! end of if( LHevectors) ; 04 Dec 14;
    
   endif ! end of if( myid.eq.0 ) ; 04 Dec 14;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
!  !work(1:Ngeometricaldof) = sqrt( evalr(1:Ngeometricaldof)**2 + evali(1:Ngeometricaldof)**2 ) ; M1 = 1 ; M2 = Ngeometricaldof ; order = 'D'
!   work(1:Ngeometricaldof) =       evalr(1:Ngeometricaldof)                                    ; M1 = 1 ; M2 = Ngeometricaldof ; order = 'D'
!   
!   im01daf = 0
!   call M01DAF( work(1:Ngeometricaldof), M1, M2, order, irank(1:Ngeometricaldof), im01daf ) ! rank by magnitude of eigenvalue
!   
!   ;   write(ounit,'("he01aa : ",f10.2," :                           ":"    "      )')cput-cpus
!   
!   do iev = 1, Ngeometricaldof
!    do jev = 1, Ngeometricaldof
!     if( irank(jev).eq.iev ) then
!      if( LHevectors ) then
!       write(ounit,'("he01aa : ",f10.2," : ",i6," : evalr="es10.2" ; ":"evecr="999f8.3)') cput-cpus, iev, evalr(jev), evecr(1:Ngeometricaldof,jev)
!      !if( evalr(jev).lt.zero ) write(ounit,'(13es13.05)') evecr(1:Ngeometricaldof,jev)
!       write(ounit,'("he01aa : ",f10.2," : ",i6," : evali="es10.2" ; ":"eveci="999f8.3)') cput-cpus, iev, evali(jev), eveci(1:Ngeometricaldof,jev)
!       write(ounit,'("he01aa : ",f10.2," :                           ":"      "       )') cput-cpus
!      else
!       write(ounit,'("he01aa : ",f10.2," : ",i6," : evalue=(",es13.5," ,",es13.5," ) ;")') cput-cpus, iev, evalr(jev), evali(jev) ! 04 Dec 14;
!      endif
!      exit
!     endif
!    enddo
!   enddo
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
!latex \item The eigenvalues and eigenvectors (if required) are written to the file \verb+.ext.GF.ev+ as follows:
!latex \begin{verbatim}
!latex open(hunit,file="."//trim(ext)//".GF.ev",status="unknown",form="unformatted")
!latex write(hunit)Ngeometricaldof,Ldvr,Ldvi        ! integers; if only the eigenvalues were computed then Ldvr=Ldvi=1;
!latex write(hunit)evalr(1:Ngeometricaldof)         ! reals   ; real      part of eigenvalues;
!latex write(hunit)evali(1:Ngeometricaldof)         ! reals   ; imaginary part of eigenvalues; 
!latex write(hunit)evecr(1:Ngeometricaldof,1:Ngeometricaldof) ! reals   ; real      part of eigenvalues; only if Ldvr=Ngeometricaldof;
!latex write(hunit)eveci(1:Ngeometricaldof,1:Ngeometricaldof) ! reals   ; imaginary part of eigenvalues; only if Ldvi=Ngeometricaldof;
!latex close(hunit)
!latex \end{verbatim}
!latex \item The eigenvectors are saved in columns of \verb+evecr, eveci+, as described by the NAG documentation for \verb+F02EBF+.

   if( myid.eq.0 ) then ! write to file; 04 Dec 14;
    open(hunit, file="."//trim(ext)//".GF.ev", status="unknown", form="unformatted")
    write(hunit) Ngeometricaldof, Ldvr, Ldvi
    write(hunit) evalr
    write(hunit) evali
    write(hunit) evecr
    write(hunit) eveci
    close(hunit)
   endif ! end of if( myid.eq.0 ) ; 04 Dec 14;

  endif ! end of if(                 ( LHevalues .or. LHevectors ) )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( myid.eq.0 .and. Lperturbed.eq.1 ) then
   
! the following will be replaced when perturbation is supplied as input; 18 Dec 14;
   
   perturbation(1:lgeometricaldof) = zero
   
   FATALMESS(he01aa, Igeometry.gt.2 .or. NOTstellsym, only for stellarator-symmetric cylindrical)
   
   do ii = 1, mn
    if( im(ii).eq.dqq .and. in(ii).eq.dpp ) perturbation(ii) = one ! impose arbitrary perturbation; 18 Dec 14;
   enddo
   
! the   above   will be replaced when perturbation is supplied as input; 18 Dec 14;
   
   if( sqrt(sum( perturbation(1:lgeometricaldof)**2 ) ).lt.vsmall ) then
    
    write(ounit,'("he01aa : " 10x " : magnitude of perturbation is trivial ;")')
    
   else
    
    rhs(1:Ngeometricaldof) = - matmul( dessian(1:Ngeometricaldof,1:lgeometricaldof), perturbation(1:lgeometricaldof) )
    
    IA = Ngeometricaldof ; IAA = Ngeometricaldof
    
    hessian(1:Ngeometricaldof,1:Ngeometricaldof) = ohessian(1:Ngeometricaldof,1:Ngeometricaldof)

    if04atf = 1
    call F04ATF( hessian(1:IA,1:Ngeometricaldof), IA, rhs(1:Ngeometricaldof), Ngeometricaldof, soln(1:Ngeometricaldof), & ! Linear solver; 02 Jan 15;
                 AA(1:IAA,1:Ngeometricaldof), IAA, rworka(1:Ngeometricaldof), rworkb(1:Ngeometricaldof), if04atf)
    
    select case( if04atf )
    case( 0 )    ; write(ounit,'("he01aa : " 10x " : myid="i3" ; linear perturbation ; if04atf="i3" ;")') myid, if04atf
    case( 1 )    ; write(ounit,'("he01aa : " 10x " : myid="i3" ; singular matrix     ; if04atf="i3" ;")') myid, if04atf
    case( 2 )    ; write(ounit,'("he01aa : " 10x " : myid="i3" ; ill-conditioned     ; if04atf="i3" ;")') myid, if04atf
    case( 3 )    ; write(ounit,'("he01aa : " 10x " : myid="i3" ; input error         ; if04atf="i3" ;")') myid, if04atf
    case default ; FATALMESS(he01aa, .true., illegal ifail returned from F04ATF)
    end select
    
    packorunpack = 'U' ! unpack geometrical degrees-of-freedom; 13 Sep 13;
    WCALL(he01aa,gf00aa,( Ngeometricaldof,     soln(0:Ngeometricaldof), Mvol, mn, &
 dRbc(1:mn,0:Mvol), dZbs(1:mn,0:Mvol), dRbs(1:mn,0:Mvol), dZbc(1:mn,0:Mvol), packorunpack ))
    
    dRbc(1:mn,Mvol) = perturbation(1:lgeometricaldof)
    
    if( Igeometry.gt.1 ) then ! include regularization factor; 18 Dec 14;
     do lvol = 1, Mvol-1
      do ii = 1, mn ; soln(ii+(lvol-1)*lgeometricaldof) = soln(ii+(lvol-1)*lgeometricaldof) * psifactor(ii,lvol) ! unpack; 29 Apr 15;
      enddo
     enddo
    endif
    
   !maxsoln = maxval( abs( soln(1:Ngeometricaldof) ) ) ; soln(1:Ngeometricaldof) = soln(1:Ngeometricaldof) / maxsoln ! normalize screen output; 18 Dec 14;
    
    if( Whe01aa ) then ! screen output; 18 Dec 14;
     ;                   ; write(ounit,'("he01aa : " 10x " : "3x" m="999( i09   ))') im(1:mn)
     ;                   ; write(ounit,'("he01aa : " 10x " : "3x" n="999( i09   ))') in(1:mn)
     do lvol = 1, Mvol-1 
      ;                  ; write(ounit,'("he01aa : " 10x " : "i3" d="999( f09.05))') lvol, ( soln((lvol-1)*lgeometricaldof+ii), ii = 1, mn )
     enddo;              ; write(ounit,'("he01aa : " 10x " : "3x" d="999( f09.05))')         perturbation(1:lgeometricaldof) !/ maxsoln
    endif ! end of if( Whe01aa ) then; 18 Dec 14;
    
   !open(hunit, file="."//trim(ext)//".perturbed", status="unknown", form="unformatted") ! this is now placed in the .h5 file; 22 Apr 15;
   !write(hunit) Ngeometricaldof, mn, Nvol, lgeometricaldof, dpp, dqq
   !write(hunit) im(1:mn)
   !write(hunit) in(1:mn)
   !write(hunit) tflux(1:Nvol)
   !write(hunit) soln(1:Ngeometricaldof)
   !write(hunit) perturbation(1:lgeometricaldof) !/ maxsoln
   !close(hunit)
    
   endif ! end of if( sqrt(sum( perturbation(1:lgeometricaldof)**2 ) ).lt.vsmall ) then; 18 Dec 14;
   
  endif ! end of if( myid.eq.0 .and. Lperturbed.eq.1 ) then; 18 Dec 14;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( myid.eq.0 ) then
   
   hessian(1:Ngeometricaldof,1:Ngeometricaldof) = ohessian(1:Ngeometricaldof,1:Ngeometricaldof)
   
   if03aaf = 1 ; IA = Ngeometricaldof
   call F03AAF( hessian(1:IA,1:Ngeometricaldof), IA, Ngeometricaldof, determinant, evalr(1:Ngeometricaldof), if03aaf) ! evalr is used as workspace; 22 Apr 15;
   
   select case( if03aaf )
   case( 0 )    ; write(ounit,'("he01aa : " 10x " : myid="i3" ; if03aaf="i3" ;             ; determinant="es13.5" ;")') myid, if03aaf, determinant
   case( 1 )    ; write(ounit,'("he01aa : " 10x " : myid="i3" ; if03aaf="i3" ; singular    ; determinant="es13.5" ;")') myid, if03aaf, determinant
   case( 2 )    ; write(ounit,'("he01aa : " 10x " : myid="i3" ; if03aaf="i3" ; overflow    ; determinant="es13.5" ;")') myid, if03aaf, determinant
   case( 3 )    ; write(ounit,'("he01aa : " 10x " : myid="i3" ; if03aaf="i3" ; underflow   ; determinant="es13.5" ;")') myid, if03aaf, determinant
   case( 4 )    ; write(ounit,'("he01aa : " 10x " : myid="i3" ; if03aaf="i3" ; input error ; determinant="es13.5" ;")') myid, if03aaf, determinant
   case default ; FATALMESS(he01aa, .true., illegal ifail returned from F03AAF)
   end select
   
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  DEALLOCATE(dFFdRZ)
  DEALLOCATE(dBBdmp)
  DEALLOCATE(dmupfdx)

  DEALLOCATE(hessian)
  DEALLOCATE(dessian)
  Lhessianallocated = .false.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(he01aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine he01aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
