!> \file
!> \brief Computes eigenvalues and eigenvectors of derivative matrix, \f$\nabla_{\bf xi}{\bf F}\f$.

!> \brief Computes eigenvalues and eigenvectors of derivative matrix, \f$\nabla_{\bf xi}{\bf F}\f$.
!> \ingroup grp_diagnostics
!>
!> @param[in] NGdof number of global degrees of freedom
!> @param[inout] position internal geometrical degrees of freedom
!> @param[in] Mvol total number of volumes in computation
!> @param[in] mn number of Fourier harmonics
!> @param[in] LGdof what is this?
subroutine hesian( NGdof, position, Mvol, mn, LGdof )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  use constants, only : zero, half, one, two, ten

  use numerical, only : sqrtmachprec, small, vsmall

  use fileunits, only : ounit, hunit, munit

  use inputlist, only : Wmacros, Whesian, Igeometry, Nvol, pflux, helicity, mu, Lfreebound, &
                        LHevalues, LHevectors, LHmatrix, &
                        Lperturbed, dpp, dqq, &
                        Lcheck, Lfindzero

  use cputiming, only : Thesian

  use allglobal, only : ncpu, myid, cpus, MPI_COMM_SPEC, ext, hiddenext, &
                        im, in, &
                        iRbc, iZbs, iRbs, iZbc, &
                        dRbc, dZbs, dRbs, dZbc, &
                        lBBintegral, dBBdRZ, &
                        NOTstellsym, YESstellsym, Energy, &
                        dFFdRZ,HdFFdRZ, dBBdmp, dmupfdx, hessian, dessian, Lhessianallocated, psifactor, &
                        hessian2D,dessian2D,Lhessian2Dallocated, &
                        Lhessian3Dallocated,denergydrr, denergydrz,denergydzr,denergydzz, &
					            	LocalConstraint
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS

  INTEGER, intent(in) :: NGdof, Mvol, mn, LGdof
  REAL                :: position(0:NGdof)

  LOGICAL             :: LComputeDerivatives, LComputeAxis

  REAL                :: force(0:NGdof), gradient(0:NGdof)

  REAL                :: xx(0:NGdof,-2:2), ff(0:NGdof,-2:2), df(1:NGdof)!, deriv

  INTEGER             :: vvol, idof, ii, mi, ni, irz, issym, isymdiff, lvol, ieval(1:1), igdof, ifd
  REAL                :: oldEnergy(-2:2), error, cpul

  REAL                :: oldBB(1:Mvol,-2:2), oBBdRZ(1:Mvol,0:1,1:LGdof), ohessian(1:NGdof,1:NGdof)

  REAL                :: oRbc(1:mn,0:Mvol), oZbs(1:mn,0:Mvol), oRbs(1:mn,0:Mvol), oZbc(1:mn,0:Mvol), determinant

  CHARACTER           :: pack

  CHARACTER           :: svol*3

!  LOGICAL            :: Lderiv                                                                         ! for parallel / series construction of Hessian;
!  INTEGER            :: lvol, jvol, ivol, innout, imn, irz, jmn, jrz, tdoc, tdof, ilocaldof, jlocaldof ! for parallel / series construction of Hessian;

  INTEGER             :: tdof, tdoc, jvol, jj, jrz, jssym

  INTEGER             :: Lwork, LDA, Ldvi, Ldvr, if02ebf
  REAL                :: evalr(1:NGdof), evali(1:NGdof)
!  REAL                :: evecr(1:NGdof,1:NGdof), eveci(1:NGdof,1:NGdof)
  REAL                :: evecr(1:NGdof,1:NGdof), eveci(1:NGdof,1:NGdof),revecr(1:NGdof,1:2*NGdof), evecl(1:NGdof,1:NGdof)
  REAL                :: work(1:4*NGdof) ! for construction of evalues/evectors;
  CHARACTER           :: JOB

  INTEGER             :: iev, jev, M1, M2, irank(1:NGdof), im01daf ! for construction of evalues/evectors;
  CHARACTER           :: order

  REAL                :: dRZ != 1.0e-03

  REAL                :: lmu(1:Mvol), lpflux(1:Mvol), lhelicity(1:Mvol) ! original profiles; 20 Jun 14;

  INTEGER             :: IA
  INTEGER             :: idgesvx, idgetrf, ipiv(1:Ngdof), iwork4(1:NGdof)
  CHARACTER           :: equed
  REAL                :: perturbation(1:LGdof)
  REAL                :: rhs(1:NGdof), solution(0:NGdof)
  REAL                :: rworka(1:NGdof), rworkb(1:NGdof), AA(1:NGdof,1:NGdof)
  REAL                :: Rdgesvx(1:NGdof), Cdgesvx(1:NGdof), AF(1:NGdof,1:NGdof), work4(1:4*NGdof), rcond, ferr, berr, sgn

  BEGIN(hesian)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  lmu(1:Nvol) = mu(1:Nvol) ; lpflux(1:Nvol) = pflux(1:Nvol) ; lhelicity(1:Nvol) = helicity(1:Nvol) ! save original profile information; 20 Jun 14;

  oldEnergy(0) = Energy ! Energy was calculated in dforce; 26 Feb 13;

  xx(0,-2:2)= zero


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef MINIMIZE

  oldBB(1:Mvol,0) = lBBintegral(1:Mvol)

  FATAL( hesian, .not.allocated(dBBdRZ), need to revise logic in preset where dBBdRZ is allocated )

  oBBdRZ(1:Mvol,0:1,1:LGdof) = dBBdRZ(1:Mvol,0:1,1:LGdof)

  oRbc(1:mn,0:Mvol) = iRbc(1:mn,0:Mvol)
  oZbs(1:mn,0:Mvol) = iZbs(1:mn,0:Mvol)
  oRbs(1:mn,0:Mvol) = iRbs(1:mn,0:Mvol)
  oZbc(1:mn,0:Mvol) = iZbc(1:mn,0:Mvol)

  FATAL( hesian, Lfreebound.eq.1, this routine needs attention )

  do vvol = 1, Mvol-1 ! loop over volumes; 26 Feb 13;

   write(ounit,'("hesian : ", 10x ," : vvol=",i3," ;")') vvol

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

       pack = 'P' !; position(0) = zero ! this is not used; 11 Aug 14;
       LComputeAxis = .true.
       LComputeDerivatives = .false. !; position(0) = zero ! this is not used; 11 Aug 14;
       WCALL( hesian, packxi, ( NGdof, position(0:NGdof), Mvol, mn, iRbc(1:mn,0:Mvol), iZbs(1:mn,0:Mvol), &
                                iRbs(1:mn,0:Mvol), iZbc(1:mn,0:Mvol), pack, LComputeDerivatives, LComputeAxis ) )


       WCALL( hesian, dforce, ( NGdof, position(0:NGdof), gradient(0:NGdof), LComputeDerivatives, LComputeAxis ) ) ! re-calculate Beltrami fields;

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
      FATAL( hesian, Igeometry.eq.1, Cartesian geometry does not need regularization factor )

1000  format("hesian : ",f10.2," : ":"myid=",i3," ; ":"vvol=",i3," ; ":"irz="i2" ; (",i3," ,",i3," ) ; "a17" ["es15.7","es15.7" ]")
1001  format("hesian : ",f10.2," : ":"myid=",i3," ; ":"vvol=",i3," ; ":"irz="i2" ; (",i3," ,",i3," ) ;   "es15.7" ; ")

     enddo ! end of do issym; 26 Feb 13;

    enddo ! end of do irz; 26 Feb 13;

   enddo ! end of do ii; 26 Feb 13;

  enddo ! end of do vvol; 26 Feb 13;

  iRbc(1:mn,0:Mvol) = oRbc(1:mn,0:Mvol)
  iZbs(1:mn,0:Mvol) = oZbs(1:mn,0:Mvol)
  iRbs(1:mn,0:Mvol) = oRbs(1:mn,0:Mvol)
  iZbc(1:mn,0:Mvol) = oZbc(1:mn,0:Mvol)

  pack = 'P' !; position(0) = zero ! this is not used; 11 Aug 14;
  WCALL( hesian, packxi,( NGdof, position(0:NGdof), Mvol, mn, iRbc(1:mn,0:Mvol), iZbs(1:mn,0:Mvol), &
                          iRbs(1:mn,0:Mvol), iZbc(1:mn,0:Mvol), pack, LComputeDerivatives, LComputeAxis ) )

#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  if(LHmatrix .and. Igeometry.eq.2) then 
      SALLOCATE( HdFFdRZ , (1:LGdof,0:1,1:LGdof,0:1,1:Mvol), zero )
  endif
  SALLOCATE( dBBdmp , (1:LGdof,1:Mvol,0:1,        1:2), zero )
  
  SALLOCATE( denergydrr, (1:LGdof,1:Mvol,0:1,1:LGdof,0:1), zero)
  !SALLOCATE( denergydrz, (1:LGdof,1:Mvol,0:1,1:LGdof,0:1), zero)
  SALLOCATE( denergydzr, (1:LGdof,1:Mvol,0:1,1:LGdof,0:1), zero)
  !SALLOCATE( denergydzz, (1:LGdof,1:Mvol,0:1,1:LGdof,0:1), zero)

if( LocalConstraint ) then
  SALLOCATE( dmupfdx, (1:Mvol,    1:1, 1:2, 1:LGdof, 0:1), zero )
else
  SALLOCATE( dmupfdx, (1:Mvol, 1:Mvol-1, 1:2, 1:LGdof, 0:1), zero)
endif

  SALLOCATE( hessian2D, (1:NGdof,1:NGdof), zero )
  SALLOCATE( dessian2D, (1:NGdof,1:LGdof), zero ) ! part of hessian that depends on boundary variations; 18 Dec 14;

  !if (LHmatrix) then
    Lhessian3Dallocated = .true.
  !else
   ! Lhessianallocated = .true.
  !endif
   !This step cleared.

  LComputeDerivatives = .true. !; position(0) = zero ! this is not used; 11 Aug 14;
  LComputeAxis = .false.
  WCALL( hesian, dforce, ( NGdof, position(0:NGdof), force(0:NGdof), LComputeDerivatives, LComputeAxis) ) ! calculate force-imbalance & hessian;
  

  ohessian(1:NGdof,1:NGdof) = hessian2D(1:NGdof,1:NGdof) ! internal copy; 22 Apr 15;


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **construction of Hessian matrix**
!> <ul>
!> <li> The routine dforce() is used to compute the derivatives, with respect to interface geometry,
!>      of the force imbalance harmonics, \f$[[p+B^2/2]]_{j}\f$, which may be considered to be the "physical" constraints,
!>      and if \c Igeometry==3 then also the derivatives of the "artificial" spectral constraints, \f$I_j \equiv (R_\theta X + Z_\theta Y)_j\f$. </li>
!> <li> The input variable \c Lconstraint determines how the enclosed fluxes, \f$\Delta \psi_t\f$ and \f$\Delta \psi_p\f$,
!>      and the helicity multiplier, \f$\mu\f$, vary as the geometry is varied;
!>      see global.f90 and mp00ac() for more details. </li>
!> </ul>

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( Lcheck.eq.5 ) then ! check construction of Hessian; 01 Jul 14;

   xx(0,-2:2)= zero ; dRZ = 1.0E-04

   write(svol,'(i3.3)')myid
!  open(lunit+myid,file=hiddenext//".hessian."//svol,status="unknown")

!  lmu(1:Nvol) = mu(1:Nvol) ; lpflux(1:Nvol) = pflux(1:Nvol) ; lhelicity(1:Nvol) = helicity(1:Nvol) ! save original profile information; 20 Jun 14;

   tdof = 0 ! geometrical degree-of-freedom counter; labels derivative;

   do vvol = 1, Mvol-1 ! loop over internal interfaces;

    cput = GETTIME
!   write(lunit+myid,'("hesian : ", 10x ," : ")')
!   write(lunit+myid,'("hesian : ",f10.2," : myid=",i3," ; vvol=",i3," ; dRZ=",es9.1," ;")') cput-cpus, myid, vvol, dRZ
!   write(lunit+myid,'("hesian : ", 10x ," : ")')

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

        xx(1:NGdof,isymdiff) = position(1:NGdof) ! reset geometry to original;

        xx(tdof,isymdiff) = position(tdof) + isymdiff * dRZ ! perturb appropriate geometric harmonic;

        LComputeDerivatives = .false.
        LComputeAxis = .true.
        WCALL( hesian, dforce, ( NGdof, xx(0:NGdof,isymdiff), ff(0:NGdof,isymdiff), LComputeDerivatives, LComputeAxis) ) ! force-imbalance;

       enddo ! end of do isymdiff; 20 Jun 14;

       df(1:NGdof) = ( - ff(1:NGdof,+2) + 8 * ff(1:NGdof,+1) - 8 * ff(1:NGdof,-1) + ff(1:NGdof,-2) ) / ( 12 * dRZ )

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

           FATAL( hesian, tdoc.lt.              1, needs attention )
           FATAL( hesian, tdoc.gt.NGdof, needs attention )

           cput = GETTIME

           error = abs( df(tdoc)-hessian(tdoc,tdof) )

           if( abs(hessian(tdoc,tdof)).gt.1.0e-05 .or. abs(df(tdoc)).gt.1.0e-05 .or. error.gt.dRZ ) then ! write to screen; 20 Jan 15;
            write(ounit     ,1001) cput-cpus, myid, &
          vvol, im(ii), in(ii), irz, issym, tdof, &
          jvol, im(jj), in(jj), jrz, jssym, tdoc, &
          df(tdoc), tdoc, tdof, hessian(tdoc,tdof), error
           endif

!          if( abs(hessian(tdoc,tdof)).gt.1.0e-05 .or. abs(df(tdoc)).gt.1.0e-05 .or. error.gt.dRZ             ) then ! write to file; 20 Jan 15;
!           write(lunit+myid,1001) cput-cpus, myid, &
!         vvol, im(ii), in(ii), irz, issym, tdof, &
!         jvol, im(jj), in(jj), jrz, jssym, tdoc, &
!         df(tdoc), tdoc, tdof, hessian(tdoc,tdof), error
!          endif

          enddo ! end of do jssym; 19 Sep 13;

         enddo ! end of do irz;

#ifdef DEBUG
!      pause
#endif

        enddo ! end of do jj;

#ifdef DEBUG
!      pause
#endif

       enddo ! end of do jvol;

#ifdef DEBUG
!      pause
#endif

1001   format("hesian : ",f10.2," : myid=",i3, &
     " ; vvol=",i3," ; ("i4" ,"i4" ); irz="i2" ; issym="i2" ; tdof="i6 &
     " ; jvol="i4" ; ("i4" ,"i4" ); jrz="i2" ; jssym="i2" ; tdoc="i6 &
     " ; fd-estimate="es13.5" & hessian("i6","i6" )="es13.5" ;":" err="es13.5" ;":,f12.4" ;")

      enddo ! end of do issym;

     enddo ! end of do irz;

    enddo ! end of do ii;

   enddo ! end of do vvol;

   pack = 'U' !; position(0) = zero ! this is not used; 11 Aug 14;
   WCALL( hesian, packxi, ( NGdof, position(0:NGdof), Mvol, mn, iRbc(1:mn,0:Mvol), iZbs(1:mn,0:Mvol), iRbs(1:mn,0:Mvol), iZbc(1:mn,0:Mvol), pack, .false., LComputeAxis ) )

   mu(1:Nvol) = lmu(1:Nvol) ; pflux(1:Nvol) = lpflux(1:Nvol) ; helicity(1:Nvol) = lhelicity(1:Nvol)

!  xx(1:NGdof,0) = position(1:NGdof) ! reset geometry to original;
!
!  LComputeDerivatives = .false.
!  WCALL(hesian,dforce,( NGdof, xx(1:NGdof,0), ff(0:NGdof,0), LComputeDerivatives )) ! calculate the force-imbalance;

   !close(lunit+myid)

  endif ! end of if( Lcheck.eq.5 ) ; 01 Jul 14;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **construction of eigenvalues and eigenvectors**
!> <ul>
!> <li> If \c LHevalues==T then the eigenvalues of the Hessian are computed using the NAG routine \c F02EBF. </li>
!> <li> If \c LHevectors==T then the eigenvalues *and* the eigenvectors of the Hessian are computed. </li>
!> <li> Note that if \c Igeometry==3, then the derivative-matrix also contains information regarding how the "artificial" spectral constraints
!>      vary with geometry; so, the eigenvalues and eigenvectors are not purely "physical". </li>
!> <li> The eigenvalues and eigenvectors (if required) are written to the file \c .ext.GF.ev as follows:
!>
!> ```
!> open(hunit,file=hiddenext//".GF.ev",status="unknown",form="unformatted")
!> write(hunit)NGdof,Ldvr,Ldvi        ! integers; if only the eigenvalues were computed then Ldvr=Ldvi=1;
!> write(hunit)evalr(1:NGdof)         ! reals   ; real      part of eigenvalues;
!> write(hunit)evali(1:NGdof)         ! reals   ; imaginary part of eigenvalues;
!> write(hunit)evecr(1:NGdof,1:NGdof) ! reals   ; real      part of eigenvalues; only if Ldvr=NGdof;
!> write(hunit)eveci(1:NGdof,1:NGdof) ! reals   ; imaginary part of eigenvalues; only if Ldvi=NGdof;
!> close(hunit)
!> ```
!> </li>
!> <li> The eigenvectors are saved in columns of \c evecr, \c eveci, as described by the NAG documentation for \c F02EBF. </li>
!> </ul>

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( LHmatrix ) then

   if( myid.eq.0 ) then ; cput = GETTIME ; write(ounit,'("hesian : ",f10.2," : LHmatrix="L2" ;")')cput-cpus, LHmatrix ;
    open(munit, file=hiddenext//".GF.ma", status="unknown", form="unformatted")
    write(munit) NGdof
    write(munit) ohessian(1:NGdof,1:NGdof)
    close(munit)
   endif

  endif


! if( myid.eq.0 .and. ( LHevalues .or. LHevectors ) ) then ! the call to dforce below requires all cpus; 04 Dec 14;
  if(                 ( LHevalues .or. LHevectors ) ) then

   if( myid.eq.0 ) then ; cput = GETTIME ; write(ounit,'("hesian : ",f10.2," : LHevalues="L2" , LHevectors="L2" ;")')cput-cpus, LHevalues, LHevectors
   endif

   evalr(1:NGdof) = zero ; evecr(1:NGdof,1:NGdof) = zero
   evali(1:NGdof) = zero ; eveci(1:NGdof,1:NGdof) = zero

   if( LHevectors ) then ; JOB='V' ; Ldvr = NGdof ; Ldvi = NGdof
   else                  ; JOB='N' ; Ldvr =  1              ; Ldvi =  1              ! provide dummy values when eigenvectors are not required; 04 Dec 14;
   endif

   cpul = GETTIME

   if02ebf = 1 ; LDA = NGdof ; Lwork = 4*NGdof
   
   hessian2D(1:NGdof,1:NGdof) = ohessian(1:NGdof,1:NGdof)

!#ifdef NAG18
!   call F02EBF( JOB, NGdof, hessian(1:LDA,1:NGdof), LDA, evalr(1:NGdof), evali(1:NGdof), &
!                evecr(1:Ldvr,1:NGdof), Ldvr, eveci(1:Ldvi,1:NGdof), Ldvi, work(1:Lwork), Lwork, if02ebf )
!#else
!   FATAL( global, .true., eigenvalue solver needs updating to F08NAF )
!#endif
   call dgeev('N', JOB, NGdof, hessian2D(1:LDA,1:NGdof), LDA, evalr(1:NGdof), evali(1:NGdof), &
                evecl(1:Ldvr,1:NGdof), Ldvr, revecr(1:Ldvr,1:2*NGdof), Ldvr, work(1:Lwork), Lwork, if02ebf )
    evecr(1:Ldvr,1:NGdof) = revecr(1:Ldvr,1:NGdof)
    eveci(1:Ldvr,1:NGdof) = revecr(1:Ldvr,NGdof+1:2*NGdof)

   if( myid.eq.0 ) then
   cput = GETTIME
     if (if02ebf < 0) then
        write(ounit,'("hesian : ",f10.2," : DGEEV error the "i2" th argument had illegal value  ;   time="f10.2"s ;")') cput-cpus, -if02ebf, cput-cpul
     else if (if02ebf > 0) then
        write(ounit,'("hesian : ",f10.2," : DGEEV error, factorization failed  ;   time="f10.2"s ;")') cput-cpus, cput-cpul
     else
        write(ounit,'("hesian : ",f10.2," : computed evalues ; if02ebf="i2" ; success ;       time="f10.2"s ;")') cput-cpus, if02ebf, cput-cpul
     endif
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!   ieval(1:1) = minloc( evalr(1:NGdof) ) ! 04 Dec 14;
!   if( myid.eq.0) then ! screen output; 04 Dec 14;
!   write(ounit,'("hesian : " 10x " : evalr("i4") ="es13.5" ="es13.5" ;")') ieval(1), evalr(ieval(1)), minval(evalr(1:NGdof))
!   endif
!
!   do ifd = -1, -6, -1 ; dRZ = ten**ifd
!
!    xx(1:NGdof,0) = position(1:NGdof) + evecr(1:NGdof,ieval(1)) * dRZ ! perturb in direction of eigenvector; 04 Dec 14;
!
!    LComputeDerivatives = .false.
!    WCALL(hesian,dforce,( NGdof, xx(0:NGdof,0), ff(0:NGdof,0), LComputeDerivatives )) ! calculate the force-imbalance;
!
!    if( myid.eq.0 ) then ! screen output; 04 Dec 14;
!     write(ounit,'("hesian : " 10x " : Energy(old)=", es23.15," ;")') oldEnergy(0)
!     write(ounit,'("hesian : " 10x " : Energy(new)=", es23.15," ;")')    Energy
!     write(ounit,'("hesian : " 10x " :             ",2es23.15," ;")')    Energy - oldEnergy(0), ( Energy-oldEnergy(0) ) / dRZ**2
!     write(ounit,'("hesian : " 10x " :             ",2es23.15," ;")')    evalr(ieval(1)) * dRZ**2, evalr(ieval(1))
!    !write(ounit,'("hesian : " 10x " : |evector|  =",es23.15," ;")')    sum(evecr(1:NGdof,ieval(1))*evecr(1:NGdof,ieval(1)))
!    !write(ounit,'("hesian : ", 10x  " : "999(" ("i3","i3")":))') (/ ( im(ii), in(ii), ii = 1, mn ) /)
!    !do lvol = 1, Mvol-1
!    ! write(ounit,'("hesian : ",i10   " : "999es10.2         )') lvol, evecr(1+mn*(lvol-1):mn+mn*(lvol-1),ieval(1))
!    !enddo
!    !write(ounit,'("hesian : " 10x " : ")')
!    !do lvol = 1, Mvol-1
!    ! write(ounit,'("hesian : ",i10   " : "999es10.2         )') lvol, ff(1+mn*(lvol-1):mn+mn*(lvol-1),0) / dRZ / evalr(ieval(1))
!    !enddo
!    endif ! end of if( myid.eq.0 ) ; 04 Dec 14;
!
!   enddo ! end of do ifd; 04 Dec 14;
!
!  !do igdof = 1, NGdof ; write(ounit,'("hesian : " 10x " : "f15.11" ="f15.11" ;")') evecr(igdof,ieval(1)), ff(igdof,0) / dRZ / evalr(ieval(1))
!  !enddo
!
!   pack = 'U' !; position(0) = zero ! this is not used; 11 Aug 14; ! reset geometry (Rbc, Zbs, etc.) to original values; 04 Dec 14;
!   WCALL(hesian,packxi,( NGdof, position(0:NGdof), Mvol, mn, &
!   iRbc(1:mn,0:Mvol), iZbs(1:mn,0:Mvol), iRbs(1:mn,0:Mvol), iZbc(1:mn,0:Mvol), pack ))
!
!   mu(1:Nvol) = lmu(1:Nvol) ; pflux(1:Nvol) = lpflux(1:Nvol) ; helicity(1:Nvol) = lhelicity(1:Nvol) ! reset profiles to original values; 04 Dec 14;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   do lvol = 1, Mvol-1
    do ii = 1, mn ; evecr(ii+(lvol-1)*mn,1:NGdof) = evecr(ii+(lvol-1)*mn,1:NGdof) * psifactor(ii,lvol) ! geometrical regularization;
    enddo
   enddo

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   if( myid.eq.0 ) then ! screen output; 04 Dec 14;

    if( LHevectors ) then

     do iev = 1, NGdof ! loop over all eigenvalues; 04 Dec 14;
      if( evalr(iev).lt.zero ) then ! only show unstable eigenvalues; 04 Dec 14;
       write(ounit,'("hesian : ",f10.2," : evalr="es13.5" ; ")') cput-cpus, evalr(iev)
       write(ounit,'("hesian : ",es10.3," : "999(" ("i3","i3")":))') evalr(iev), (/ ( im(ii), in(ii), ii = 1, mn ) /)
       do lvol = 1, Mvol-1 ; write(ounit,'("hesian : ",i10   " : "999es10.2)') lvol, evecr(1+mn*(lvol-1):mn+mn*(lvol-1),iev)
       enddo
      endif
     enddo

    else ! matches if( LHvectors) ; 04 Dec 14;

     do iev = 1, NGdof
     !if( evalr(iev).lt.zero ) then ! 22 Apr 15;
       write(ounit,'("hesian : ",f10.2," : ",i6," : evalue=(",es23.15," ,",es23.15," ) ;")') cput-cpus, iev, evalr(iev), evali(iev) ! 04 Dec 14;
     !endif ! 22 Apr 15;
     enddo

    endif ! end of if( LHevectors) ; 04 Dec 14;

   endif ! end of if( myid.eq.0 ) ; 04 Dec 14;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!  !work(1:NGdof) = sqrt( evalr(1:NGdof)**2 + evali(1:NGdof)**2 ) ; M1 = 1 ; M2 = NGdof ; order = 'D'
!   work(1:NGdof) =       evalr(1:NGdof)                                    ; M1 = 1 ; M2 = NGdof ; order = 'D'
!
!   im01daf = 0
!   call M01DAF( work(1:NGdof), M1, M2, order, irank(1:NGdof), im01daf ) ! rank by magnitude of eigenvalue
!
!   ;   write(ounit,'("hesian : ",f10.2," :                           ":"    "      )')cput-cpus
!
!   do iev = 1, NGdof
!    do jev = 1, NGdof
!     if( irank(jev).eq.iev ) then
!      if( LHevectors ) then
!       write(ounit,'("hesian : ",f10.2," : ",i6," : evalr="es10.2" ; ":"evecr="999f8.3)') cput-cpus, iev, evalr(jev), evecr(1:NGdof,jev)
!      !if( evalr(jev).lt.zero ) write(ounit,'(13es13.05)') evecr(1:NGdof,jev)
!       write(ounit,'("hesian : ",f10.2," : ",i6," : evali="es10.2" ; ":"eveci="999f8.3)') cput-cpus, iev, evali(jev), eveci(1:NGdof,jev)
!       write(ounit,'("hesian : ",f10.2," :                           ":"      "       )') cput-cpus
!      else
!       write(ounit,'("hesian : ",f10.2," : ",i6," : evalue=(",es13.5," ,",es13.5," ) ;")') cput-cpus, iev, evalr(jev), evali(jev) ! 04 Dec 14;
!      endif
!      exit
!     endif
!    enddo
!   enddo

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


   if( myid.eq.0 ) then ! write to file; 04 Dec 14;
    open(hunit, file=hiddenext//".GF.ev", status="unknown", form="unformatted")
    write(hunit) NGdof, Ldvr, Ldvi
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

   perturbation(1:LGdof) = zero

   FATAL( hesian, Igeometry.gt.2 .or. NOTstellsym, only for stellarator-symmetric cylindrical )

   do ii = 1, mn
    if( im(ii).eq.dqq .and. in(ii).eq.dpp ) perturbation(ii) = one ! impose arbitrary perturbation; 18 Dec 14;
   enddo

! the   above   will be replaced when perturbation is supplied as input; 18 Dec 14;

   if( sqrt(sum( perturbation(1:LGdof)**2 ) ).lt.vsmall ) then

    write(ounit,'("hesian : " 10x " : magnitude of perturbation is trivial ;")')

   else

    rhs(1:NGdof) = - matmul( dessian(1:NGdof,1:LGdof), perturbation(1:LGdof) )
        
    hessian2D(1:NGdof,1:NGdof) = ohessian(1:NGdof,1:NGdof)

    call dgesvx( 'N', 'N', NGdof, 1, hessian2D(1:NGdof,1:NGdof), NGdof, AF(1:NGdof,1:NGdof),   & ! Linear solver; 09 Nov 17;
                 NGdof, ipiv(1:NGdof), equed, Rdgesvx(1:NGdof), Cdgesvx(1:NGdof),            & 
         rhs(1:NGdof), NGdof, solution(1:NGdof), NGdof, rcond, ferr, berr,           &
         work4(1:4*NGdof), iwork4(1:NGdof), idgesvx )

    select case( idgesvx )
    case( 0   )    ; write(ounit,'("hesian : " 10x " : myid="i3" ; linear perturbation ; idgesvx="i3" ;")') myid, idgesvx
    case( 1:  )    ; write(ounit,'("hesian : " 10x " : myid="i3" ; singular matrix     ; idgesvx="i3" ;")') myid, idgesvx
    case( :-1 )    ; write(ounit,'("hesian : " 10x " : myid="i3" ; input error         ; idgesvx="i3" ;")') myid, idgesvx
    case default ; FATAL( hesian, .true., illegal ifail returned from dgesvx )
    end select

    pack = 'U' ! unpack geometrical degrees-of-freedom; 13 Sep 13;
    WCALL( hesian, packxi, ( NGdof,     solution(0:NGdof), Mvol, mn, dRbc(1:mn,0:Mvol), dZbs(1:mn,0:Mvol), dRbs(1:mn,0:Mvol), dZbc(1:mn,0:Mvol), pack, .false., LComputeAxis ) )

    dRbc(1:mn,Mvol) = perturbation(1:LGdof)

    if( Igeometry.gt.1 ) then ! include regularization factor; 18 Dec 14;
     do lvol = 1, Mvol-1
      do ii = 1, mn ; solution(ii+(lvol-1)*LGdof) = solution(ii+(lvol-1)*LGdof) * psifactor(ii,lvol) ! unpack; 29 Apr 15;
      enddo
     enddo
    endif

    if( Whesian ) then ! screen output; 18 Dec 14;
     ;                   ; write(ounit,'("hesian : " 10x " : "3x" m="999( i09   ))') im(1:mn)
     ;                   ; write(ounit,'("hesian : " 10x " : "3x" n="999( i09   ))') in(1:mn)
     do lvol = 1, Mvol-1
      ;                  ; write(ounit,'("hesian : " 10x " : "i3" d="999( f09.05))') lvol, ( solution((lvol-1)*LGdof+ii), ii = 1, mn )
     enddo;              ; write(ounit,'("hesian : " 10x " : "3x" d="999( f09.05))')         perturbation(1:LGdof)
    endif ! end of if( Whesian ) then; 18 Dec 14;

   endif ! end of if( sqrt(sum( perturbation(1:LGdof)**2 ) ).lt.vsmall ) then; 18 Dec 14;

  endif ! end of if( myid.eq.0 .and. Lperturbed.eq.1 ) then; 18 Dec 14;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( myid.eq.0 ) then
   
   hessian2D(1:NGdof,1:NGdof) = ohessian(1:NGdof,1:NGdof)
   
   call dgetrf( NGdof, NGdof, hessian2D(1:NGdof,1:NGdof), NGdof, ipiv(1:NGdof), idgetrf )
     
   determinant = one

   do iev = 1,NGdof
    determinant = determinant*hessian2D(iev,iev)   !calculate determinant from factorized form of hessian; 09 Nov 17
   enddo

   sgn = one

   do iev = 1,NGdof
    if(ipiv(iev).ne. iev) then
    sgn = -sgn
    endif
   enddo

   determinant = sgn*determinant                 !correct for the sign of the determinant; 09 Nov 17

   select case( idgetrf )
   case( 0   )    ; write(ounit,'("hesian : " 10x " : myid="i3" ; idgetrf="i3" ;             ; determinant="es13.5" ;")') myid, idgetrf, determinant
   case( 1:  )    ; write(ounit,'("hesian : " 10x " : myid="i3" ; idgetrf="i3" ; singular    ; determinant="es13.5" ;")') myid, idgetrf, determinant
   case( :-1 )    ; write(ounit,'("hesian : " 10x " : myid="i3" ; idgetrf="i3" ; input error ; determinant="es13.5" ;")') myid, idgetrf, determinant
   case default ; FATAL( hesian, .true., illegal ifail returned from dgetrf )
   end select

  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  if(LHmatrix .and. Igeometry.eq.2) then 
     DALLOCATE(HdFFdRZ)
  endif
  
  DALLOCATE(dBBdmp)
  DALLOCATE(dmupfdx)
  DALLOCATE(denergydrr)
  !DALLOCATE(denergydrz)
  DALLOCATE(denergydzr)
  !DALLOCATE(denergydzz)
  Lhessian3Dallocated=.false.


  DALLOCATE(hessian2D)
  write(ounit,*) 5656

  DALLOCATE(dessian2D)


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(hesian)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine hesian

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
