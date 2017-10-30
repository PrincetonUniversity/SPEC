!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Returns the magnetic field field line equations, and the tangent map if required.

!latex \item Sometimes it is required to have the full contravariant components of the field, so the quantity \verb+Bzeta +$\equiv B^\zeta$ is returned through memory.

!latex \end{enumerate} \subsubsection{equations of field line flow} \begin{enumerate}

!latex \item The field line equations are $\dot\theta=B^\theta/B^\zeta$, $\dot s=B^s/B^\zeta$.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! IT IS REQUIRED TO SET IVOL AND LTANGENT THROUGH GLOBAL MEMORY BEFORE CALLING BF00AA;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine bf00aa( zeta, st, Bst ) ! the format of this subroutine is constrained by the NAG ode integration routines. 
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, one, half, two, pi2, pi
  
  use numerical, only : vsmall, small
  
  use fileunits, only : ounit
  
  use inputlist, only : Wmacros, Wbf00aa, Nvol, Lrad, curtor, curpol
  
  use cputiming, only : Tbf00aa
  
  use allglobal, only : myid, ncpu, cpus, mn, im, in, halfmm, &
                        ivol, Ltangent, Bzeta, Ate, Aze, Ato, Azo, &
                        NOTstellsym, lchebyshev, &
                        Lcoordinatesingularity, Lvacuumregion, Lplasmaregion, Mvol
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  REAL, intent(in)  :: zeta, st(6)
  REAL, intent(out) :: Bst(6)
  
  INTEGER           :: lvol, imn, ii, jj, kk, ll, mi, ni, ideriv
  REAL              :: teta, lss, sbar, sbarhim(0:2), arg, carg, sarg, dBu(1:3,0:3), dphi(1:3,0:3), TM(2,2), Bsumn, Btumn, Bzumn, Bumn, phimn
  REAL              :: stz(1:3), RpZ(1:3), dR(0:3), dZ(0:3), jacobian, glow(1:3,1:3), gupp(1:3,1:3)
  
  REAL, allocatable :: ssTT(:,:)

  BEGIN(bf00aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  lvol = ivol ; ideriv = 0 ! the argument list of bf00aa is fixed by NAG requirements, but volume index is required by below;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef DEBUG
  FATALMESS(bf00aa, lvol.lt.1 .or. lvol.gt.Mvol, invalid ivol )
#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( abs(st(1)).gt.one ) then ; Bst(1:6) = (/ zero , zero , zero , zero , zero , zero /) ; goto 9999 ! out of domain; should speed things up . . . ;
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  lss = st(1) ; teta = st(2) ; sbar = ( lss + one ) * half ! shorthand; 20 Apr 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  lchebyshev(0,0:2) = (/ one, zero, zero /) ! Chebyshev initialization; 16 Jan 13;
  lchebyshev(1,0:2) = (/ lss,  one, zero /) ! Chebyshev initialization; 16 Jan 13;
  do ll = 2, Lrad(lvol)
   lchebyshev(ll,0:2) = (/ two * lss * lchebyshev(ll-1,0)                                  - lchebyshev(ll-2,0) , & ! Chebyshev recurrence;            ; 16 Jan 13;
                           two       * lchebyshev(ll-1,0) + two * lss * lchebyshev(ll-1,1) - lchebyshev(ll-2,1) , &
                           two * two * lchebyshev(ll-1,1) + two * lss * lchebyshev(ll-1,2) - lchebyshev(ll-2,2) /)  ! Chebyshev recurrence; derivatives; 16 Jan 13;
  enddo

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Lvacuumregion ) then ! will need metrics to convert covariant field to contravariant field; 22 Apr 13;
   
   stz(1:3) = (/ lss, teta, zeta /)
   
   CALL(bf00aa, co00aa,( lvol, stz(1:3), RpZ(1:3), dR(0:3), dZ(0:3), jacobian, glow(1:3,1:3) )) ! map to cylindrical;

   do ii = 1, 3
    
    select case( ii )
    case( 1 ) ; jj = 2 ; kk = 3
    case( 2 ) ; jj = 3 ; kk = 1
    case( 3 ) ; jj = 1 ; kk = 2
    end select
    
    gupp(ii,1) = glow(jj,2)*glow(kk,3) - glow(jj,3)*glow(kk,2) ! have ignored the jacobian factors, which will cancel; 22 Apr 13;
    gupp(ii,2) = glow(jj,3)*glow(kk,1) - glow(jj,1)*glow(kk,3)
    gupp(ii,3) = glow(jj,1)*glow(kk,2) - glow(jj,2)*glow(kk,1)
    
!    gupp(2,1) = glow(3,2)*glow(1,3) - glow(3,3)*glow(1,2)
!    gupp(2,2) = glow(3,3)*glow(1,1) - glow(3,1)*glow(1,3)
!    gupp(2,3) = glow(3,1)*glow(1,2) - glow(3,2)*glow(1,1)
!    
!    gupp(3,1) = glow(1,2)*glow(2,3) - glow(1,3)*glow(2,2)
!    gupp(3,2) = glow(1,3)*glow(2,1) - glow(1,1)*glow(2,3)
!    gupp(3,3) = glow(1,1)*glow(2,2) - glow(1,2)*glow(2,1)
   
   enddo ! end of do ii; 17 Apr 13;

  endif ! end of if( Lvacuumregion ) ; 17 Apr 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RALLOCATE(ssTT,(0:Lrad(lvol),0:2))

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  dBu(1:3,0:3) = zero ! initialize summation;

  dphi(1:3,0:3) = zero ! initialize summation; ! only required in vacuum region; 17 Apr 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  do ii = 1, mn ; mi = im(ii) ; ni = in(ii) ; arg = mi*teta - ni*zeta ; carg = cos(arg) ; sarg = sin(arg) ! shorthand; 20 Apr 13;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   if( Lcoordinatesingularity ) then
    
   !FATALMESS(bf00aa, .true., there should be a m-1 factor) ! 16 Jan 15;

   !sbarhim(0) = sbar**halfmm(ii) ; sbarhim(1) = half * halfmm(ii) * sbarhim(0) / sbar ; sbarhim(2) = half *   halfmm(ii)         * sbarhim(1) / sbar ! 16 Jan 15;
    sbarhim(0) = sbar**halfmm(ii) ; sbarhim(1) = half * halfmm(ii) * sbarhim(0) / sbar ; sbarhim(2) = half * ( halfmm(ii) - one ) * sbarhim(1) / sbar ! 16 Jan 15;

    do ll = 0, Lrad(lvol) ! shorthand for regularization factor times Chebyshev polynomials; 20 Apr 13;
     ssTT(ll,0:2) = (/ sbarhim(0) * lchebyshev(ll,0), &
                       sbarhim(0) * lchebyshev(ll,1) + sbarhim(1) * lchebyshev(ll,0), &
                       sbarhim(0) * lchebyshev(ll,2) + sbarhim(1) * lchebyshev(ll,1) * two + sbarhim(2) * lchebyshev(ll,0) /)
    enddo

   else

   !sbarhim(0) = one                 ; sbarhim(1) = zero                                     ; sbarhim(2) = zero

    do ll = 0, Lrad(lvol) ! shorthand for regularization factor times Chebyshev polynomials; 20 Apr 13;
     ssTT(ll,0:2) = (/              lchebyshev(ll,0), &
                                    lchebyshev(ll,1)                                , &
                                    lchebyshev(ll,2)                                                                       /) ! 16 Jan 15;
    enddo
    
   endif ! end of if( Lcoordinatesingularity ) ; 16 Jan 15;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
!  do ll = 0, Lrad(lvol) ! shorthand for regularization factor times Chebyshev polynomials; 20 Apr 13;
!   ssTT(ll,0:2) = (/ sbarhim(0) * lchebyshev(ll,0), &
!                     sbarhim(0) * lchebyshev(ll,1) + sbarhim(1) * lchebyshev(ll,0), &
!                     sbarhim(0) * lchebyshev(ll,2) + sbarhim(1) * lchebyshev(ll,1) * two + sbarhim(2) * lchebyshev(ll,0) /)
!  enddo

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   if( Lplasmaregion ) then
    
    do ll = 0, Lrad(lvol) ! loop over Chebyshev summation; 20 Feb 13;
     dBu(1,0) = dBu(1,0) + ( - mi * Aze(lvol,ideriv,ii)%s(ll) - ni * Ate(lvol,ideriv,ii)%s(ll) ) * ssTT(ll,0) * sarg
     dBu(2,0) = dBu(2,0) + (                                  -      Aze(lvol,ideriv,ii)%s(ll) ) * ssTT(ll,1) * carg
     dBu(3,0) = dBu(3,0) + (        Ate(lvol,ideriv,ii)%s(ll)                                  ) * ssTT(ll,1) * carg
     if( NOTstellsym ) then ! include non-symmetric harmonics; 28 Jan 13;
     dBu(1,0) = dBu(1,0) + ( + mi * Azo(lvol,ideriv,ii)%s(ll) + ni * Ato(lvol,ideriv,ii)%s(ll) ) * ssTT(ll,0) * carg
     dBu(2,0) = dBu(2,0) + (                                  -      Azo(lvol,ideriv,ii)%s(ll) ) * ssTT(ll,1) * sarg
     dBu(3,0) = dBu(3,0) + (        Ato(lvol,ideriv,ii)%s(ll)                                  ) * ssTT(ll,1) * sarg
     endif
    enddo

   endif ! end of if( Lplasmaregion ) ; 17 Apr 13;
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
!latex \item The magnetic field in the vacuum region is 
!latex       \be {\bf B} = I \nabla \t + G \nabla \z + \phi_s \nabla \s + \phi_\t \nabla \t + \phi_z \nabla \z.
!latex       \ee

   if( Lvacuumregion ) then ! recall that this is inside do ii = 1, mn loop; 01 May 13;
    
    do ll = 0, Lrad(lvol) ! loop over Chebyshev summation; 20 Feb 13;

!     phimn =        Ate(lvol,ideriv,ii)%s(ll) ; dphi(1,0) = dphi(1,0) + phimn * ssTT(ll,1) * sarg
!     phimn =   mi * Ate(lvol,ideriv,ii)%s(ll) ; dphi(2,0) = dphi(2,0) + phimn * ssTT(ll,0) * carg
!     phimn = - ni * Ate(lvol,ideriv,ii)%s(ll) ; dphi(3,0) = dphi(3,0) + phimn * ssTT(ll,0) * carg
!     if( NOTstellsym ) then ! include non-symmetric harmonics; 28 Jan 13;
!     phimn =        Ato(lvol,ideriv,ii)%s(ll) ; dphi(1,0) = dphi(1,0) + phimn * ssTT(ll,1) * carg
!     phimn = - mi * Ato(lvol,ideriv,ii)%s(ll) ; dphi(2,0) = dphi(2,0) + phimn * ssTT(ll,0) * sarg
!     phimn =   ni * Ato(lvol,ideriv,ii)%s(ll) ; dphi(3,0) = dphi(3,0) + phimn * ssTT(ll,0) * sarg
!     endif

     dphi(1,0) = dphi(1,0) +      Ate(lvol,ideriv,ii)%s(ll) * ssTT(ll,1) * sarg
     dphi(2,0) = dphi(2,0) + mi * Ate(lvol,ideriv,ii)%s(ll) * ssTT(ll,0) * carg
     dphi(3,0) = dphi(3,0) - ni * Ate(lvol,ideriv,ii)%s(ll) * ssTT(ll,0) * carg
     if( NOTstellsym ) then ! include non-symmetric harmonics; 28 Jan 13;
     dphi(1,0) = dphi(1,0) +      Ato(lvol,ideriv,ii)%s(ll) * ssTT(ll,1) * carg
     dphi(2,0) = dphi(2,0) - mi * Ato(lvol,ideriv,ii)%s(ll) * ssTT(ll,0) * sarg
     dphi(3,0) = dphi(3,0) + ni * Ato(lvol,ideriv,ii)%s(ll) * ssTT(ll,0) * sarg
     endif

    !write(ounit,'("bf00aa : ", 10x ," : ii=",i3,", mi=",i3,", ni=",i3,", ll=",i3,", Ate="es13.5" ;")') ii, mi, ni, ll, Ate(lvol,ideriv,ii)%s(ll)

    enddo

   endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   if( Ltangent.eq.1 ) then
    
    FATALMESS(bf00aa, Lvacuumregion, tangent vacuum field under construction)
    
    do ll = 0, Lrad(lvol)
     
     Bsumn = ( - mi    * Aze(lvol,ideriv,ii)%s(ll) - ni    * Ate(lvol,ideriv,ii)%s(ll) ) ; dBu(1,1) = dBu(1,1) + Bsumn * ssTT(ll,1) * sarg
     Btumn = (                                     -         Aze(lvol,ideriv,ii)%s(ll) ) ; dBu(2,1) = dBu(2,1) + Btumn * ssTT(ll,2) * carg
     Bzumn = ( +         Ate(lvol,ideriv,ii)%s(ll)                                     ) ; dBu(3,1) = dBu(3,1) + Bzumn * ssTT(ll,2) * carg 

     Bsumn = ( - mi*mi * Aze(lvol,ideriv,ii)%s(ll) - ni*mi * Ate(lvol,ideriv,ii)%s(ll) ) ; dBu(1,2) = dBu(1,2) + Bsumn * ssTT(ll,0) * carg
     Btumn = (                                     +    mi * Aze(lvol,ideriv,ii)%s(ll) ) ; dBu(2,2) = dBu(2,2) + Btumn * ssTT(ll,1) * sarg
     Bzumn = ( -    mi * Ate(lvol,ideriv,ii)%s(ll)                                     ) ; dBu(3,2) = dBu(3,2) + Bzumn * ssTT(ll,1) * sarg 

     Bsumn = ( + mi*ni * Aze(lvol,ideriv,ii)%s(ll) + ni*ni * Ate(lvol,ideriv,ii)%s(ll) ) ; dBu(1,3) = dBu(1,3) + Bsumn * ssTT(ll,0) * carg
     Btumn = (                                     -    ni * Aze(lvol,ideriv,ii)%s(ll) ) ; dBu(2,3) = dBu(2,3) + Btumn * ssTT(ll,1) * sarg
     Bzumn = ( +    ni * Ate(lvol,ideriv,ii)%s(ll)                                     ) ; dBu(3,3) = dBu(3,3) + Bzumn * ssTT(ll,1) * sarg 

    enddo
    
    if( NOTstellsym ) then
     
     do ll = 0, Lrad(lvol)
      
      Bsumn = ( + mi    * Azo(lvol,ideriv,ii)%s(ll) + ni    * Ato(lvol,ideriv,ii)%s(ll) ) ; dBu(1,1) = dBu(1,1) + Bsumn * ssTT(ll,1) * carg
      Btumn = (                                     -         Azo(lvol,ideriv,ii)%s(ll) ) ; dBu(2,1) = dBu(2,1) + Btumn * ssTT(ll,2) * sarg
      Bzumn = ( +         Ato(lvol,ideriv,ii)%s(ll)                                     ) ; dBu(3,1) = dBu(3,1) + Bzumn * ssTT(ll,2) * sarg 
      
      Bsumn = ( - mi*mi * Azo(lvol,ideriv,ii)%s(ll) - ni*mi * Ato(lvol,ideriv,ii)%s(ll) ) ; dBu(1,2) = dBu(1,2) + Bsumn * ssTT(ll,0) * sarg
      Btumn = (                                     -    mi * Azo(lvol,ideriv,ii)%s(ll) ) ; dBu(2,2) = dBu(2,2) + Btumn * ssTT(ll,1) * carg
      Bzumn = ( +    mi * Ato(lvol,ideriv,ii)%s(ll)                                     ) ; dBu(3,2) = dBu(3,2) + Bzumn * ssTT(ll,1) * carg 
      
      Bsumn = ( - mi*ni * Azo(lvol,ideriv,ii)%s(ll) - ni*ni * Ato(lvol,ideriv,ii)%s(ll) ) ; dBu(1,3) = dBu(1,3) + Bsumn * ssTT(ll,0) * sarg
      Btumn = (                                     -    ni * Azo(lvol,ideriv,ii)%s(ll) ) ; dBu(2,3) = dBu(2,3) + Btumn * ssTT(ll,1) * carg
      Bzumn = ( +    ni * Ato(lvol,ideriv,ii)%s(ll)                                     ) ; dBu(3,3) = dBu(3,3) + Bzumn * ssTT(ll,1) * carg 
      
     enddo
     
    endif
    
   endif ! end of if( Ltangent.eq.1 ) then;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
  enddo ! end of do ii = 1, mn;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Lvacuumregion ) then ! need to raise covariant components; 22 Apr 13;

   do ii = 1, 3 ; dBu(ii,0) = curtor*gupp(2,ii) + curpol*gupp(3,ii) + dphi(1,0)*gupp(1,ii) + dphi(2,0)*gupp(2,ii) + dphi(3,0)*gupp(3,ii) ! 01 May 13;
   enddo

!  dBu(1,0) = curtor*gupp(2,1)  +  curpol*gupp(3,1)  +  dphi(1,0)*gupp(1,1)  +  dphi(2,0)*gupp(2,1)  +  dphi(3,0)*gupp(3,1) ! to be deleted; 17 Apr 13;
!  dBu(2,0) = curtor*gupp(2,2)  +  curpol*gupp(3,2)  +  dphi(1,0)*gupp(1,2)  +  dphi(2,0)*gupp(2,2)  +  dphi(3,0)*gupp(3,2)
!  dBu(3,0) = curtor*gupp(2,3)  +  curpol*gupp(3,3)  +  dphi(1,0)*gupp(1,3)  +  dphi(2,0)*gupp(2,3)  +  dphi(3,0)*gupp(3,3)

  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( abs(dBu(3,0)).lt.vsmall ) then
   cput = GETTIME
   write(ounit,'("bf00aa : ",f10.2," : lvol=",i3," ; zeta="es23.15" ; (s,t)=("es23.15" ,"es23.15" ) ; B^z="es23.15" ;")') cput-cpus, lvol, zeta, st(1:2), dBu(3,0)
   FATALMESS(bf00aa, abs(dBu(3,0)).lt.vsmall, field is not toroidal)
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  Bzeta = dBu(3,0)                            ! Bzeta is returned through global; 20 Apr 13;
  
  Bst(1:2) = (/ dBu(1,0), dBu(2,0) /) / Bzeta ! normalize field line equations to toroidal field; 20 Apr 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Ltangent.eq.0 ) then ; Bst(3:6) = zero ; goto 9999 ! tangent map is not required;
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \end{enumerate} \subsubsection{definition of tangent map} \begin{enumerate}

!latex \item The tangent map satisfies 
!latex       \be \frac{d}{d \z} M = 
!latex       \left(\begin{array}{cc} \partial_\s \dot \s & \partial_\t \dot \s \\
!latex                               \partial_\s \dot \t & \partial_\t \dot \t 
!latex               \end{array}
!latex       \right)
!latex       M. \ee
!latex \item Note that the tangent map is only constructed if \verb+Ltangent.eq.1+.
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  TM(1,1) = ( dBu(1,1) - Bst(1) * dBu(3,1) ) / dBu(3,0) !  radial  derivative of  psi  dot;
  TM(1,2) = ( dBu(1,2) - Bst(1) * dBu(3,2) ) / dBu(3,0) ! poloidal derivative of  psi  dot;
  TM(2,1) = ( dBu(2,1) - Bst(2) * dBu(3,1) ) / dBu(3,0) !  radial  derivative of theta dot;
  TM(2,2) = ( dBu(2,2) - Bst(2) * dBu(3,2) ) / dBu(3,0) ! poloidal derivative of theta dot;
  
  Bst(3) = TM(1,1) * st(3) + TM(1,2) * st(5) ! tangent map obtained by matrix multiplication;
  Bst(4) = TM(1,1) * st(4) + TM(1,2) * st(6)
  Bst(5) = TM(2,1) * st(3) + TM(2,2) * st(5)
  Bst(6) = TM(2,1) * st(4) + TM(2,2) * st(6)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  DEALLOCATE(ssTT)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(bf00aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine bf00aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!subroutine bf00aa_out( zeta, st )
!
!  use constants, only : pi2
!  use inputlist, only : Lpoincare
!  use fileunits, only : ounit
!  use allglobal, only : pi2nfp, Lfieldlinedirection!, iqfmin, jgd, kgd
!
!  LOCALS
!
!  REAL, intent(inout) :: zeta
!  REAL, intent(in)    :: st(8)
!  
!  if( Lpoincare.eq.0 ) zeta = zeta + pi2nfp ! next intermediate output location;
!  if( Lpoincare.eq.1 ) zeta = zeta + pi2    ! next intermediate output location;
!  if( Lpoincare.eq.2 ) zeta = zeta + pi2nfp ! next intermediate output location;
!
!  return
!
!end subroutine bf00aa_out

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!REAL function bf00aa_end( zeta, st )
!
!  use constants, only : pi2
!  use inputlist, only : Lpoincare
!  use allglobal, only : pi2nfp, Nz
!
!  LOCALS
!
!  REAL, intent(inout) :: zeta
!  REAL, intent(in)    :: st(8)
!  
!  if( Lpoincare.eq.0 ) bf00aa_end = zeta   - pi2nfp    ! integration termination;
!  if( Lpoincare.eq.1 ) bf00aa_end = zeta   - pi2       ! integration termination;
!  if( Lpoincare.eq.2 ) bf00aa_end = abs(st(8)) - pi2nfp/Nz ! integration termination;
!
!  return
!
!end FUNCTION bf00aa_end

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
