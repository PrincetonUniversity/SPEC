subroutine fndiff( NGdof )

use constants, only: zero, one, half, two

use fileunits, only: ounit

use inputlist, only: Wmacros, Wfndiff, ext, &
                     Igeometry, &
                     dRZ, Lcheck

use cputiming, only: Tfndiff

use allglobal, only: ncpu, myid, cpus, &
                     Mvol, mn, im, in, &
                     iRbc, iZbs, iRbs, iZbc, &
                     LGdof, psifactor, dBdX, &
                     YESstellsym, NOTstellsym, &
                     hessian

LOCALS

INTEGER, intent(in) :: NGdof

INTEGER             :: vvol, idof, ii, irz, issym, isymdiff ! loop indices
INTEGER             :: tdof ! hessian index

REAL                :: lfactor
CHARACTER           :: packorunpack 
LOGICAL             :: LComputeAxis

REAL, allocatable   :: oRbc(:,:), oZbs(:,:), oRbs(:,:), oZbc(:,:) ! used to store original geometry;
REAL, allocatable   :: iposition(:,:), iforce(:,:) ! perturbed interfaces position and force
REAL, allocatable   :: finitediff_estimate(:,:)  ! store finite differences

BEGIN(fndiff)

if( Lcheck.eq.6 ) then
  SALLOCATE( finitediff_estimate, (1:NGdof, 1:NGdof), zero )
endif

select case( Lcheck )
  
case( 4 ) ! Finite differences of dmupfdx ( see dfp200.f90 )
  continue
case( 6 ) ! Finite differences of hessian ( see dforce.f90 )
      
  write(ounit, '("fndiff : Starting finite difference evaluation of hessian ...")')

  if( ncpu.eq.1) then

  dBdX%L = .false.
  SALLOCATE( oRbc, (1:mn,0:Mvol), iRbc(1:mn,0:Mvol) ) !save unperturbed geometry
  SALLOCATE( oZbs, (1:mn,0:Mvol), iZbs(1:mn,0:Mvol) )
  SALLOCATE( oRbs, (1:mn,0:Mvol), iRbs(1:mn,0:Mvol) )
  SALLOCATE( oZbc, (1:mn,0:Mvol), iZbc(1:mn,0:Mvol) ) 
  SALLOCATE( iforce,  (-2:2, 0:NGdof), zero)
  SALLOCATE( iposition, (-2:2, 0:NGdof), zero)
  
  
  do vvol = 1, Mvol-1 ! loop over interior surfaces;
    idof = 0
    
    do ii = 1, mn ! Loop over Fourier modes

      lfactor = psifactor(ii,vvol)   ! this "pre-conditions" the geometrical degrees-of-freedom;
    
      do irz = 0, 1 ! loop over R or Z coordinate
            
        if( irz.eq.1 .and. Igeometry.lt.3 ) cycle

        do issym = 0, 1 ! stellarator symmetry;

          if( issym.eq.1 .and. YESstellsym ) cycle ! no dependence on the non-stellarator symmetric harmonics;

          if( ii.eq.1 .and. irz.eq.1 .and. issym.eq.0 ) cycle ! no dependence on Zbs_{m=0,n=0};
          if( ii.eq.1 .and. irz.eq.0 .and. issym.eq.1 ) cycle ! no dependence on Rbs_{m=0,n=0};
      
          idof = idof + 1 ! labels degree-of-freedom;

          do isymdiff = -2, 2 ! symmetric fourth-order, finite-difference used to approximate derivatives;

            if( isymdiff.eq.0 ) cycle

            ! Reset initial geometry
            iRbc(1:mn,0:Mvol) = oRbc(1:mn,0:Mvol)
            iZbs(1:mn,0:Mvol) = oZbs(1:mn,0:Mvol)
            iRbs(1:mn,0:Mvol) = oRbs(1:mn,0:Mvol)
            iZbc(1:mn,0:Mvol) = oZbc(1:mn,0:Mvol)

            ! Perturb geometry
            if( issym.eq.0 .and. irz.eq.0 ) then
              iRbc(ii,vvol) = iRbc(ii,vvol) + dRZ * isymdiff ! perturb geometry;
            else if( issym.eq.0 .and. irz.eq.1 ) then
              iZbs(ii,vvol) = iZbs(ii,vvol) + dRZ * isymdiff ! perturb geometry;
            else if( issym.eq.1 .and. irz.eq.0 ) then
              iRbs(ii,vvol) = iRbs(ii,vvol) + dRZ * isymdiff ! perturb geometry;
            else if( issym.eq.1 .and. irz.eq.1 ) then
              iZbc(ii,vvol) = iZbc(ii,vvol) + dRZ * isymdiff ! perturb geometry;
            endif

            packorunpack = 'P' ! pack geometrical degrees-of-freedom;
            LComputeAxis = .false. ! keep axis fixed

            WCALL(fndiff, packxi,( NGdof, iposition(isymdiff,0:NGdof), Mvol, mn,iRbc(1:mn,0:Mvol),iZbs(1:mn,0:Mvol),iRbs(1:mn,0:Mvol),&
                                   iZbc(1:mn,0:Mvol),packorunpack, .false., LComputeAxis ) )
            WCALL(fndiff, dforce,( NGdof, iposition(isymdiff,0:NGdof), iforce(isymdiff,0:NGdof), .false., LComputeAxis) )
          enddo

          ! Fourth order centered finite difference scheme
          iforce(0, 0:NGdof)  = ( - 1 * iforce( 2,0:NGdof) &
                                  + 8 * iforce( 1,0:NGdof) &
                                  - 8 * iforce(-1,0:NGdof) &
                                  + 1 * iforce(-2,0:NGdof))  / ( 12 * dRZ )
            
          tdof = (vvol-1) * LGdof + idof
          finitediff_estimate(1:NGdof, tdof) = iforce(0, 1:NGdof)* lfactor

        enddo !issym
      enddo !irz
    enddo !ii
  enddo !vvol
        

  DALLOCATE(iforce)
  DALLOCATE(iposition)
  DALLOCATE(oZbc)
  DALLOCATE(oRbs)
  DALLOCATE(oZbs)
  DALLOCATE(oRbc)

  ! Print in file for diagnostics
  if(myid.eq.0) then
    ! Print hessian
    open(10, file=trim(ext)//'.Lcheck6_output.txt', status='unknown')
    write(ounit,'(A)') NEW_LINE('A')
  
    do ii=1, NGdof
      write(ounit,1345) myid, im(ii), in(ii), hessian(ii,:)
      write(10   ,1347) hessian(ii,:)
    enddo
    close(10)
        
    write(ounit,'(A)') NEW_LINE('A')

    ! Print finite differences
    open(10, file=trim(ext)//'.Lcheck6_output.FiniteDiff.txt', status='unknown')
    do ii=1, NGdof
      write(ounit,1347) myid, im(ii), in(ii), finitediff_estimate(ii,:)
      write(10   ,1347) finitediff_estimate(ii,:)
    enddo        
    write(ounit,'(A)') NEW_LINE('A')
    close(10)

    1345 format("dforce: myid=",i3," ; (",i4,",",i4," ; Hessian            = ",512f16.10 "   ;")
    1346 format("dforce: myid=",i3," ; (",i4,",",i4," ; Finite differences = ",512f16.10 "   ;")
    1347 format(512F22.16, " ")

  endif

  endif

case( 8 ) ! Finite differences of axis  ( see rzaxis.f90 )

  continue
!   threshold = 1e-8 ! print with difference between FD and analytical more than this threshold
!   dx = 1e-8 * jRbc(1,ivol)

!   newRbc = jRbc
!   newRbs = jRbs
!   newZbc = jZbc
!   newZbs = jZbs

!   if (irz .eq. 0 .and. issym .eq. 0) then
!     newRbc(imn, ivol) = jRbc(imn, ivol) + dx
!   else if (irz .eq. 1 .and. issym .eq. 0) then
!     newZbs(imn, ivol) = jZbs(imn, ivol) + dx
!   else if (irz .eq. 0 .and. issym .eq. 1) then
!     newRbs(imn, ivol) = jRbs(imn, ivol) + dx
!   else if (irz .eq. 1 .and. issym .eq. 1) then
!     newZbc(imn, ivol) = jZbc(imn, ivol) + dx
!   end if
    
!   ! call the same subroutine recursively, but do not compute derivatives
!   call rzaxis( Mvol, mn, newRbc, newZbs, newRbs, newZbc, ivol, .false. )

!   ! compare the derivatives
!   do ii = 1, Ntoraxis+1
!     if (irz.eq.0) then
!       if (abs((newRbc(ii,0) - jRbc(ii,jvol))/dx -  dRadR(ii,0,issym,imn))/jRbc(1,ivol) .ge. threshold) then
!         write(ounit, *) 'dRc/dR: ii,m,n,issym', ii, im(imn), in(imn),issym, (newRbc(ii,0) - jRbc(ii,jvol))/dx, dRadR(ii,0,issym,imn), dx, newRbc(ii,0), jRbc(ii,jvol)
!       endif
!       if (abs((newZbs(ii,0) - jZbs(ii,jvol))/dx -  dZadR(ii,1,issym,imn))/jRbc(1,ivol) .ge. threshold) then
!         write(ounit, *) 'dZs/dR: ii,m,n,issym', ii, im(imn), in(imn),issym, (newZbs(ii,0) - jZbs(ii,jvol))/dx, dZadR(ii,1,issym,imn), dx
!       endif
!       if (NOTstellsym) then
!         if (abs((newRbs(ii,0) - jRbs(ii,jvol))/dx -  dRadR(ii,1,issym,imn))/jRbc(1,ivol) .ge. threshold) then
!           write(ounit, *) 'dRs/dR: ii,m,n,issym', ii, im(imn), in(imn),issym, (newRbs(ii,0) - jRbs(ii,jvol))/dx, dRadR(ii,1,issym,imn), dx
!         endif
!         if (abs((newZbc(ii,0) - jZbc(ii,jvol))/dx -  dZadR(ii,0,issym,imn))/jRbc(1,ivol) .ge. threshold) then
!           write(ounit, *) 'dZc/dR: ii,m,n,issym', ii, im(imn), in(imn),issym, (newZbc(ii,0) - jZbc(ii,jvol))/dx, dZadR(ii,0,issym,imn), dx
!         endif
!       endif
!     else if (irz.eq.1) then
!       if (abs((newRbc(ii,0) - jRbc(ii,jvol))/dx -  dRadZ(ii,0,1-issym,imn))/jRbc(1,ivol) .ge. threshold) then
!         write(ounit, *) 'dRc/dZ: ii,m,n,issym', ii, im(imn), in(imn),issym, (newRbc(ii,0) - jRbc(ii,jvol))/dx, dRadZ(ii,0,1-issym,imn), dx
!       endif
!       if (abs((newZbs(ii,0) - jZbs(ii,jvol))/dx -  dZadZ(ii,1,1-issym,imn))/jRbc(1,ivol) .ge. threshold) then
!         write(ounit, *) 'dZs/dZ: ii,m,n,issym', ii, im(imn), in(imn),issym, (newZbs(ii,0) - jZbs(ii,jvol))/dx, dZadZ(ii,1,1-issym,imn), dx
!       endif
!       if (abs((newRbs(ii,0) - jRbs(ii,jvol))/dx -  dRadZ(ii,1,1-issym,imn))/jRbc(1,ivol) .ge. threshold) then
!         write(ounit, *) 'dRs/dZ: ii,m,n,issym', ii, im(imn), in(imn),issym, (newRbs(ii,0) - jRbs(ii,jvol))/dx, dRadZ(ii,1,1-issym,imn), dx
!       endif
!       if (abs((newZbc(ii,0) - jZbc(ii,jvol))/dx -  dZadZ(ii,0,1-issym,imn))/jRbc(1,ivol) .ge. threshold) then
!         write(ounit, *) 'dZc/dZ: ii,m,n,issym', ii, im(imn), in(imn),issym, (newZbc(ii,0) - jZbc(ii,jvol))/dx, dZadZ(ii,0,1-issym,imn), dx
!       endif
!     endif
!   enddo



case default
  
  FATAL( fndiff, .true., Invalid Lcheck )

end select


DALLOCATE(finitediff_estimate)
FATAL(fndiff, .true., Finite differences have been evaluated. )

RETURN(fndiff)

end subroutine fndiff