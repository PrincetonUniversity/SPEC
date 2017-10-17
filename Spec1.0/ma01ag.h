!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Constructs matrix, and the perturbed matrix, that represents the Beltrami linear system.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \end{enumerate} \subsection{matrix elements} \begin{enumerate}

!latex \item The energy, $W \equiv \int \! dv {\; \bf B}\cdot{\bf B}$, and helicity, $K\equiv \int \! dv \; {\bf A}\cdot{\bf B}$, functionals may be written
!latex \be W & = & \frac{1}{2} \; a_i \; A_{i,j} \; a_j + a_i \; B_{i,j} \; \psi_j + \frac{1}{2} \; \psi_i \; C_{i,j} \; \psi_j \label{eq:energymatrix} \\
!latex     K & = & \frac{1}{2} \; a_i \; D_{i,j} \; a_j + a_i \; E_{i,j} \; \psi_j + \frac{1}{2} \; \psi_i \; F_{i,j} \; \psi_j \label{eq:helicitymatrix}
!latex \ee
!latex       where ${\bf a} \equiv \{ A_{\t,e,i,l}, A_{\z,e,i,l}, A_{\t,o,i,l}, A_{\z,o,i,l}, f_{e,i}, f_{o,i} \}$ contains the independent degrees of freedom
!latex       and $\boldpsi \equiv \{\Delta \psi_t,\Delta \psi_p\}$.

!latex \item The matrix elements are computed via
!latex       \be \verb+MA(i,j)+ & \equiv & A_{i,j} = \frac{\partial^2 W}{\partial    a_i \partial    a_j} \\
!latex           \verb+MB(i,j)+ & \equiv & B_{i,j} = \frac{\partial^2 W}{\partial    a_i \partial \psi_j} \\
!latex           \verb+MC(i,j)+ & \equiv & C_{i,j} = \frac{\partial^2 W}{\partial \psi_i \partial \psi_j}
!latex       \ee

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine ma01ag( lvol, mn, lrad )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half, quart, one

  use numerical, only : sqrtmachprec

  use fileunits, only : ounit, wunit

  use inputlist, only : Wmacros, Wma01ag, Nvol

  use cputiming, only : Tma01ag

  use allglobal, only : myid, cpus, &
                        NOTstellsym, &
                        im, in, &
                        Nmagneticdof, &
                        dMA, dMB, dMC, dMD, dME, dMF, &
                        Ate, Ato, Aze, Azo, Fso, Fse, &
                        TTee, TTeo, TToe, TToo, &
                        Lcoordinatesingularity, &
                        Llatex
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER, intent(in)  :: lvol, mn, lrad

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER              :: NN, ii, jj, ll, pp, mi, ni, mj, nj, mimj, minj, nimj, ninj, mjmi, mjni, njmi, njni, id, jd, Lul, Lup, Lzl, Lzp

  REAL                 :: Wtete(1:mn,0:lrad,1:mn,0:lrad), Wteto(1:mn,0:lrad,1:mn,0:lrad), Wtote(1:mn,0:lrad,1:mn,0:lrad), Wtoto(1:mn,0:lrad,1:mn,0:lrad)
  REAL                 :: Wteze(1:mn,0:lrad,1:mn,0:lrad), Wtezo(1:mn,0:lrad,1:mn,0:lrad), Wtoze(1:mn,0:lrad,1:mn,0:lrad), Wtozo(1:mn,0:lrad,1:mn,0:lrad)
  REAL                 :: Wzete(1:mn,0:lrad,1:mn,0:lrad), Wzeto(1:mn,0:lrad,1:mn,0:lrad), Wzote(1:mn,0:lrad,1:mn,0:lrad), Wzoto(1:mn,0:lrad,1:mn,0:lrad)
  REAL                 :: Wzeze(1:mn,0:lrad,1:mn,0:lrad), Wzezo(1:mn,0:lrad,1:mn,0:lrad), Wzoze(1:mn,0:lrad,1:mn,0:lrad), Wzozo(1:mn,0:lrad,1:mn,0:lrad)

  REAL                 :: Htete(1:mn,0:lrad,1:mn,0:lrad), Hteto(1:mn,0:lrad,1:mn,0:lrad), Htote(1:mn,0:lrad,1:mn,0:lrad), Htoto(1:mn,0:lrad,1:mn,0:lrad)
  REAL                 :: Hteze(1:mn,0:lrad,1:mn,0:lrad), Htezo(1:mn,0:lrad,1:mn,0:lrad), Htoze(1:mn,0:lrad,1:mn,0:lrad), Htozo(1:mn,0:lrad,1:mn,0:lrad)
  REAL                 :: Hzete(1:mn,0:lrad,1:mn,0:lrad), Hzeto(1:mn,0:lrad,1:mn,0:lrad), Hzote(1:mn,0:lrad,1:mn,0:lrad), Hzoto(1:mn,0:lrad,1:mn,0:lrad)
  REAL                 :: Hzeze(1:mn,0:lrad,1:mn,0:lrad), Hzezo(1:mn,0:lrad,1:mn,0:lrad), Hzoze(1:mn,0:lrad,1:mn,0:lrad), Hzozo(1:mn,0:lrad,1:mn,0:lrad)

  LOGICAL              :: Lzi, Lzj

 !write(0,'("ma01ag : " 10x  " : myid="i3" :                : lvol=",i3," ; mn=",i3," ; lrad      =",i3," ; ")') myid, lvol, mn, lrad

  BEGIN(ma01ag)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef DEBUG
  FATALMESS(ma01ag, lvol.lt.1 .or. lvol.gt.Nvol, invalid volume label)
  FATALMESS(ma01ag, .not.allocated(dMA), error)
  FATALMESS(ma01ag, .not.allocated(dMB), error)
  FATALMESS(ma01ag, .not.allocated(dMC), error)
  FATALMESS(ma01ag, .not.allocated(dMD), error)
  FATALMESS(ma01ag, .not.allocated(dME), error)
  FATALMESS(ma01ag, .not.allocated(dMF), error)
#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  NN = Nmagneticdof(lvol) ! shorthand; Nmagneticdof was computed in al00aa; 16 Jan 13;

  dMA(0:NN,0:NN) = zero ! initialize summation; 24 Jan 13;
  dMB(0:NN,1: 2) = zero
  dMC(1: 2,1: 2) = zero

  dMD(0:NN,0:NN) = zero
  dME(0:NN,1: 2) = zero
  dMF(1: 2,1: 2) = zero
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  Wtete(1:mn,0:lrad,1:mn,0:lrad) = zero ! these "zeros" are probably not required; 10 Mar 13;
  Wteto(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Wtote(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Wtoto(1:mn,0:lrad,1:mn,0:lrad) = zero
  
  Wteze(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Wtezo(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Wtoze(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Wtozo(1:mn,0:lrad,1:mn,0:lrad) = zero
  
  Wzete(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Wzeto(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Wzote(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Wzoto(1:mn,0:lrad,1:mn,0:lrad) = zero
  
  Wzeze(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Wzezo(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Wzoze(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Wzozo(1:mn,0:lrad,1:mn,0:lrad) = zero
  
  Htete(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Hteto(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Htote(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Htoto(1:mn,0:lrad,1:mn,0:lrad) = zero
  
  Hteze(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Htezo(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Htoze(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Htozo(1:mn,0:lrad,1:mn,0:lrad) = zero
  
  Hzete(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Hzeto(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Hzote(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Hzoto(1:mn,0:lrad,1:mn,0:lrad) = zero
  
  Hzeze(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Hzezo(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Hzoze(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Hzozo(1:mn,0:lrad,1:mn,0:lrad) = zero

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG

  if( Wma01ag ) then
   write(ounit,'("ma01ag : ", 10x ," : myid=",i3," ; lvol=",i3," ; computing second derivatives ;")') myid, lvol
  endif

!  if( Llatex ) then
!   write(wunit,'("\newpage \subsection{{second derivatives of energy and helicity integrals : volume="i3"; [ma01ag]}}")') lvol
!   write(wunit,'("\begin{{enumerate}}")')
!  endif

#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
   
   do jj = 1, mn ; mj = im(jj) ; nj = in(jj)
    
    mimj = mi * mj ; minj = mi * nj ; nimj = ni * mj ; ninj = ni * nj
    
    do ll = 0, lrad
     
     do pp = 0, lrad
      
      Wtete(jj,pp,ii,ll) = + 2 * ninj * TToo(1,ll,pp,ii,jj) - 2 * ni      * TToe(3,ll,pp,ii,jj) - 2      * nj * TToe(3,pp,ll,jj,ii) + 2 * TTee(6,ll,pp,ii,jj)
      Wzete(jj,pp,ii,ll) = + 2 * nimj * TToo(1,ll,pp,ii,jj) + 2 * ni      * TToe(2,ll,pp,ii,jj) - 2      * mj * TToe(3,pp,ll,jj,ii) - 2 * TTee(5,pp,ll,jj,ii)

      Htete(jj,pp,ii,ll) =   zero
      Hzete(jj,pp,ii,ll) = - TTee(0,pp,ll,jj,ii) + TTee(0,ll,pp,ii,jj)

      Wteze(jj,pp,ii,ll) = + 2 * minj * TToo(1,ll,pp,ii,jj) + 2      * nj * TToe(2,pp,ll,jj,ii) - 2 * mi      * TToe(3,ll,pp,ii,jj) - 2 * TTee(5,ll,pp,ii,jj)
      Wzeze(jj,pp,ii,ll) = + 2 * mimj * TToo(1,ll,pp,ii,jj) + 2 * mi      * TToe(2,ll,pp,ii,jj) + 2      * mj * TToe(2,pp,ll,jj,ii) + 2 * TTee(4,ll,pp,ii,jj)

      Hteze(jj,pp,ii,ll) = - TTee(0,ll,pp,ii,jj) + TTee(0,pp,ll,jj,ii)  
      Hzeze(jj,pp,ii,ll) =   zero  

!#ifdef DEBUG
!      if( Llatex ) then
!       if( abs(Wtete(jj,pp,ii,ll)).gt.sqrtmachprec ) then
!       write(wunit,'("\item $\ds\frac{{\partial^2 W}}{{\partial A_{{\t,"i2","i2"}}\partial A_{{\t,"i2","i2"}}}}="f9.3"$")') jj, pp, ii, ll, Wtete(jj,pp,ii,ll)
!       endif
!       if( abs(Wzete(jj,pp,ii,ll)).gt.sqrtmachprec ) then
!       write(wunit,'("\item $\ds\frac{{\partial^2 W}}{{\partial A_{{\z,"i2","i2"}}\partial A_{{\t,"i2","i2"}}}}="f9.3"$")') jj, pp, ii, ll, Wzete(jj,pp,ii,ll)
!       endif
!       if( abs(Wteze(jj,pp,ii,ll)).gt.sqrtmachprec ) then
!       write(wunit,'("\item $\ds\frac{{\partial^2 W}}{{\partial A_{{\t,"i2","i2"}}\partial A_{{\z,"i2","i2"}}}}="f9.3"$")') jj, pp, ii, ll, Wteze(jj,pp,ii,ll)
!       endif
!       if( abs(Wzeze(jj,pp,ii,ll)).gt.sqrtmachprec ) then
!       write(wunit,'("\item $\ds\frac{{\partial^2 W}}{{\partial A_{{\z,"i2","i2"}}\partial A_{{\z,"i2","i2"}}}}="f9.3"$")') jj, pp, ii, ll, Wzeze(jj,pp,ii,ll)
!       endif
!       if( abs(Hzete(jj,pp,ii,ll)).gt.sqrtmachprec ) then
!       write(wunit,'("\item $\ds\frac{{\partial^2 H}}{{\partial A_{{\z,"i2","i2"}}\partial A_{{\t,"i2","i2"}}}}="f9.3"$")') jj, pp, ii, ll, Hzete(jj,pp,ii,ll)
!       endif
!       if( abs(Hteze(jj,pp,ii,ll)).gt.sqrtmachprec ) then
!       write(wunit,'("\item $\ds\frac{{\partial^2 H}}{{\partial A_{{\t,"i2","i2"}}\partial A_{{\z,"i2","i2"}}}}="f9.3"$")') jj, pp, ii, ll, Hteze(jj,pp,ii,ll)
!       endif
!      endif
!#endif
      
      if( NOTstellsym ) then

!     Wtote(jj,pp,ii,ll) = - 2 * ninj * TToe(1,pp,ll,ii,jj) + 2      * nj * TTee(3,pp,ll,jj,ii) - 2 * ni      * TToo(3,ll,pp,ii,jj) + 2 * TToe(6,ll,pp,jj,ii)
      Wtote(jj,pp,ii,ll) = - 2 * ninj * TTeo(1,pp,ll,jj,ii) + 2      * nj * TTee(3,pp,ll,jj,ii) - 2 * ni      * TToo(3,ll,pp,ii,jj) + 2 * TTeo(6,ll,pp,ii,jj)
!     Wzote(jj,pp,ii,ll) = - 2 * nimj * TToe(1,pp,ll,ii,jj) + 2 * ni      * TToo(2,ll,pp,ii,jj) + 2      * mj * TTee(3,pp,ll,jj,ii) - 2 * TToe(5,pp,ll,jj,ii)
      Wzote(jj,pp,ii,ll) = - 2 * nimj * TTeo(1,pp,ll,jj,ii) + 2 * ni      * TToo(2,ll,pp,ii,jj) + 2      * mj * TTee(3,pp,ll,jj,ii) - 2 * TToe(5,pp,ll,jj,ii)

      Htote(jj,pp,ii,ll) =   zero
!     Hzote(jj,pp,ii,ll) = - TToe(0,pp,ll,jj,ii) + TToe(0,ll,pp,jj,ii)
      Hzote(jj,pp,ii,ll) = - TToe(0,pp,ll,jj,ii) + TTeo(0,ll,pp,ii,jj)

!     Wteto(jj,pp,ii,ll) = - 2 * ninj * TToe(1,ll,pp,jj,ii) + 2 * ni      * TTee(3,ll,pp,ii,jj) - 2      * nj * TToo(3,pp,ll,jj,ii) + 2 * TToe(6,pp,ll,ii,jj)
      Wteto(jj,pp,ii,ll) = - 2 * ninj * TTeo(1,ll,pp,ii,jj) + 2 * ni      * TTee(3,ll,pp,ii,jj) - 2      * nj * TToo(3,pp,ll,jj,ii) + 2 * TTeo(6,pp,ll,jj,ii)
!     Wtoto(jj,pp,ii,ll) = + 2 * ninj * TTee(1,ll,pp,ii,jj) + 2 * ni      * TToe(3,ll,pp,jj,ii) + 2      * nj * TToe(3,pp,ll,ii,jj) + 2 * TToo(6,ll,pp,ii,jj)
      Wtoto(jj,pp,ii,ll) = + 2 * ninj * TTee(1,ll,pp,ii,jj) + 2 * ni      * TTeo(3,ll,pp,ii,jj) + 2      * nj * TTeo(3,pp,ll,jj,ii) + 2 * TToo(6,ll,pp,ii,jj)
!     Wzeto(jj,pp,ii,ll) = - 2 * nimj * TToe(1,ll,pp,jj,ii) - 2 * ni      * TTee(2,ll,pp,ii,jj) - 2      * mj * TToo(3,pp,ll,jj,ii) - 2 * TToe(5,pp,ll,ii,jj)
      Wzeto(jj,pp,ii,ll) = - 2 * nimj * TTeo(1,ll,pp,ii,jj) - 2 * ni      * TTee(2,ll,pp,ii,jj) - 2      * mj * TToo(3,pp,ll,jj,ii) - 2 * TTeo(5,pp,ll,jj,ii)
!     Wzoto(jj,pp,ii,ll) = + 2 * nimj * TTee(1,ll,pp,ii,jj) - 2 * ni      * TToe(2,ll,pp,jj,ii) + 2      * mj * TToe(3,pp,ll,ii,jj) - 2 * TToo(5,pp,ll,jj,ii)
      Wzoto(jj,pp,ii,ll) = + 2 * nimj * TTee(1,ll,pp,ii,jj) - 2 * ni      * TTeo(2,ll,pp,ii,jj) + 2      * mj * TTeo(3,pp,ll,jj,ii) - 2 * TToo(5,pp,ll,jj,ii)

      Hteto(jj,pp,ii,ll) =   zero
      Htoto(jj,pp,ii,ll) =   zero
!     Hzeto(jj,pp,ii,ll) = - TToe(0,pp,ll,ii,jj) + TToe(0,ll,pp,ii,jj)
      Hzeto(jj,pp,ii,ll) = - TTeo(0,pp,ll,jj,ii) + TToe(0,ll,pp,ii,jj)
      Hzoto(jj,pp,ii,ll) = - TToo(0,pp,ll,jj,ii) + TToo(0,ll,pp,ii,jj)  

!     Wtoze(jj,pp,ii,ll) = - 2 * minj * TToe(1,pp,ll,ii,jj) - 2      * nj * TTee(2,pp,ll,jj,ii) - 2 * mi      * TToo(3,ll,pp,ii,jj) - 2 * TToe(5,ll,pp,jj,ii)
      Wtoze(jj,pp,ii,ll) = - 2 * minj * TTeo(1,pp,ll,jj,ii) - 2      * nj * TTee(2,pp,ll,jj,ii) - 2 * mi      * TToo(3,ll,pp,ii,jj) - 2 * TTeo(5,ll,pp,ii,jj)
!     Wzoze(jj,pp,ii,ll) = - 2 * mimj * TToe(1,pp,ll,ii,jj) - 2      * mj * TTee(2,pp,ll,jj,ii) + 2 * mi      * TToo(2,ll,pp,ii,jj) + 2 * TToe(4,ll,pp,jj,ii)
      Wzoze(jj,pp,ii,ll) = - 2 * mimj * TTeo(1,pp,ll,jj,ii) - 2      * mj * TTee(2,pp,ll,jj,ii) + 2 * mi      * TToo(2,ll,pp,ii,jj) + 2 * TTeo(4,ll,pp,ii,jj)

!     Htoze(jj,pp,ii,ll) = - TToe(0,ll,pp,jj,ii) + TToe(0,pp,ll,jj,ii)
      Htoze(jj,pp,ii,ll) = - TTeo(0,ll,pp,ii,jj) + TToe(0,pp,ll,jj,ii)
      Hzoze(jj,pp,ii,ll) =   zero

!     Wtezo(jj,pp,ii,ll) = - 2 * minj * TToe(1,ll,pp,jj,ii) + 2      * nj * TToo(2,pp,ll,jj,ii) + 2 * mi      * TTee(3,ll,pp,ii,jj) - 2 * TToe(5,ll,pp,ii,jj)
      Wtezo(jj,pp,ii,ll) = - 2 * minj * TTeo(1,ll,pp,ii,jj) + 2      * nj * TToo(2,pp,ll,jj,ii) + 2 * mi      * TTee(3,ll,pp,ii,jj) - 2 * TToe(5,ll,pp,ii,jj)
!     Wtozo(jj,pp,ii,ll) = + 2 * minj * TTee(1,ll,pp,ii,jj) - 2      * nj * TToe(2,pp,ll,ii,jj) + 2 * mi      * TToe(3,ll,pp,jj,ii) - 2 * TToo(5,ll,pp,ii,jj)
      Wtozo(jj,pp,ii,ll) = + 2 * minj * TTee(1,ll,pp,ii,jj) - 2      * nj * TTeo(2,pp,ll,jj,ii) + 2 * mi      * TTeo(3,ll,pp,ii,jj) - 2 * TToo(5,ll,pp,ii,jj)
!     Wzezo(jj,pp,ii,ll) = - 2 * mimj * TToe(1,ll,pp,jj,ii) - 2 * mi      * TTee(2,ll,pp,ii,jj) + 2      * mj * TToo(2,pp,ll,jj,ii) + 2 * TToe(4,pp,ll,ii,jj)
      Wzezo(jj,pp,ii,ll) = - 2 * mimj * TTeo(1,ll,pp,ii,jj) - 2 * mi      * TTee(2,ll,pp,ii,jj) + 2      * mj * TToo(2,pp,ll,jj,ii) + 2 * TTeo(4,pp,ll,jj,ii)
!     Wzozo(jj,pp,ii,ll) = + 2 * mimj * TTee(1,ll,pp,ii,jj) - 2 * mi      * TToe(2,ll,pp,jj,ii) - 2      * mj * TToe(2,pp,ll,ii,jj) + 2 * TToo(4,ll,pp,ii,jj)
      Wzozo(jj,pp,ii,ll) = + 2 * mimj * TTee(1,ll,pp,ii,jj) - 2 * mi      * TTeo(2,ll,pp,ii,jj) - 2      * mj * TTeo(2,pp,ll,jj,ii) + 2 * TToo(4,ll,pp,ii,jj)

!     Htezo(jj,pp,ii,ll) = - TToe(0,ll,pp,ii,jj) + TToe(0,pp,ll,ii,jj)
      Htezo(jj,pp,ii,ll) = - TToe(0,ll,pp,ii,jj) + TTeo(0,pp,ll,jj,ii)
      Htozo(jj,pp,ii,ll) = - TToo(0,ll,pp,ii,jj) + TToo(0,pp,ll,jj,ii)  
      Hzezo(jj,pp,ii,ll) =   zero
      Hzozo(jj,pp,ii,ll) =   zero

      endif ! end of if( NOTstellsym ) 26 Feb 13;
      
     enddo ! end of do pp;
     
    enddo ! end of do jj;
    
   enddo ! end of do ll;
   
  enddo ! end of do ii; 

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!#ifdef DEBUG
!
!  if( Llatex ) then
!   write(wunit,'("\end{{enumerate}}")')
!  endif
!
!#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! construct matrix elements;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!#ifdef DEBUG
!  if( Wma01ag ) write(ounit,'("ma01ag : ", 10x ," : myid=",i3," ; lvol=",i3," ; 0.0 ;")') myid, lvol
!#endif

  do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
   
  !if( Lcoordinatesingularity .and. ( mi.ne.0 .or. ni.eq.0 ) ) then ; Lzi = .true. ! 15 Jan 15; ! additional freedom; 15 Jan 15;
   if( Lcoordinatesingularity .and. ( mi.ne.0 .or. ii.eq.1 ) ) then ; Lzi = .true. ! 15 Jan 15; ! additional freedom; 15 Jan 15;
   else                                                             ; Lzi = .false.
   endif
   
   do ll = 0, lrad-1 ! this is a derivative loop;
    
    if( (ll/2)*2.eq.ll ) then ; Lul = lrad   ! ll is even; 18 Jan 13;
    else                      ; Lul = lrad-1 ! ll is odd ; 18 Jan 13;
    endif

    if( Lzi ) then ; Lzl = lrad
    else           ; Lzl = Lul
    endif
    
    do jj = 1, mn ; mj = im(jj) ; nj = in(jj)
     
     if( Lcoordinatesingularity .and. ( mj.ne.0 .or. nj.eq.0 ) ) then ; Lzj = .true.
     else                                                             ; Lzj = .false.
     endif
     
     do pp = 0, lrad-1 ! this is a derivative loop;
      
      if( (pp/2)*2.eq.pp ) then ; Lup = lrad   ! pp is even; 18 Jan 13;
      else                      ; Lup = lrad-1 ! pp is odd ; 18 Jan 13;
      endif
      
      if( Lzj ) then ; Lzp = lrad
      else           ; Lzp = Lup
      endif
      
#ifdef DEBUG
      FATALMESS( ma01ag, Ate(lvol,0,ii)%i(ll).lt.0 .or. Ate(lvol,0,ii)%i(ll).gt.NN, invalid subscript )
      FATALMESS( ma01ag, Ate(lvol,0,jj)%i(pp).lt.0 .or. Ate(lvol,0,jj)%i(pp).gt.NN, invalid subscript )
      FATALMESS( ma01ag, Aze(lvol,0,ii)%i(ll).lt.0 .or. Aze(lvol,0,ii)%i(ll).gt.NN, invalid subscript )
      FATALMESS( ma01ag, Aze(lvol,0,jj)%i(pp).lt.0 .or. Aze(lvol,0,jj)%i(pp).gt.NN, invalid subscript )
      if( NOTstellsym ) then
      FATALMESS( ma01ag, Ato(lvol,0,ii)%i(ll).lt.0 .or. Ato(lvol,0,ii)%i(ll).gt.NN, invalid subscript )
      FATALMESS( ma01ag, Ato(lvol,0,jj)%i(pp).lt.0 .or. Ato(lvol,0,jj)%i(pp).gt.NN, invalid subscript )
      FATALMESS( ma01ag, Azo(lvol,0,ii)%i(ll).lt.0 .or. Azo(lvol,0,ii)%i(ll).gt.NN, invalid subscript )
      FATALMESS( ma01ag, Azo(lvol,0,jj)%i(pp).lt.0 .or. Azo(lvol,0,jj)%i(pp).gt.NN, invalid subscript )
      endif
#endif

      id=Ate(lvol,0,ii)%i(ll) ; jd=Ate(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + Wtete(jj,pp,ii,ll) - Wtete(jj,pp,ii,Lul) - Wtete(jj,Lup,ii,ll) + Wtete(jj,Lup,ii,Lul)

      id=Ate(lvol,0,ii)%i(ll) ; jd=Aze(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + Wzete(jj,pp,ii,ll) - Wzete(jj,pp,ii,Lul) - Wzete(jj,Lzp,ii,ll) + Wzete(jj,Lzp,ii,Lul)

      id=Aze(lvol,0,ii)%i(ll) ; jd=Ate(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + Wteze(jj,pp,ii,ll) - Wteze(jj,pp,ii,Lzl) - Wteze(jj,Lup,ii,ll) + Wteze(jj,Lup,ii,Lzl)

      id=Aze(lvol,0,ii)%i(ll) ; jd=Aze(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + Wzeze(jj,pp,ii,ll) - Wzeze(jj,pp,ii,Lzl) - Wzeze(jj,Lzp,ii,ll) + Wzeze(jj,Lzp,ii,Lzl)

      id=Ate(lvol,0,ii)%i(ll) ; jd=Ate(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + Htete(jj,pp,ii,ll) - Htete(jj,pp,ii,Lul) - Htete(jj,Lup,ii,ll) + Htete(jj,Lup,ii,Lul)

      id=Ate(lvol,0,ii)%i(ll) ; jd=Aze(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + Hzete(jj,pp,ii,ll) - Hzete(jj,pp,ii,Lul) - Hzete(jj,Lzp,ii,ll) + Hzete(jj,Lzp,ii,Lul)

      id=Aze(lvol,0,ii)%i(ll) ; jd=Ate(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + Hteze(jj,pp,ii,ll) - Hteze(jj,pp,ii,Lzl) - Hteze(jj,Lup,ii,ll) + Hteze(jj,Lup,ii,Lzl)

      id=Aze(lvol,0,ii)%i(ll) ; jd=Aze(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + Hzeze(jj,pp,ii,ll) - Hzeze(jj,pp,ii,Lzl) - Hzeze(jj,Lzp,ii,ll) + Hzeze(jj,Lzp,ii,Lzl)

      if( NOTstellsym ) then

      id=Ate(lvol,0,ii)%i(ll) ; jd=Ato(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + Wtote(jj,pp,ii,ll) - Wtote(jj,pp,ii,Lul) - Wtote(jj,Lup,ii,ll) + Wtote(jj,Lup,ii,Lul)

      id=Ate(lvol,0,ii)%i(ll) ; jd=Azo(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + Wzote(jj,pp,ii,ll) - Wzote(jj,pp,ii,Lul) - Wzote(jj,Lzp,ii,ll) + Wzote(jj,Lzp,ii,Lul)

      id=Ato(lvol,0,ii)%i(ll) ; jd=Ate(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + Wteto(jj,pp,ii,ll) - Wteto(jj,pp,ii,Lul) - Wteto(jj,Lup,ii,ll) + Wteto(jj,Lup,ii,Lul)

      id=Ato(lvol,0,ii)%i(ll) ; jd=Ato(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + Wtoto(jj,pp,ii,ll) - Wtoto(jj,pp,ii,Lul) - Wtoto(jj,Lup,ii,ll) + Wtoto(jj,Lup,ii,Lul)

      id=Ato(lvol,0,ii)%i(ll) ; jd=Aze(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + Wzeto(jj,pp,ii,ll) - Wzeto(jj,pp,ii,Lul) - Wzeto(jj,Lzp,ii,ll) + Wzeto(jj,Lzp,ii,Lul)

      id=Ato(lvol,0,ii)%i(ll) ; jd=Azo(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + Wzoto(jj,pp,ii,ll) - Wzoto(jj,pp,ii,Lul) - Wzoto(jj,Lzp,ii,ll) + Wzoto(jj,Lzp,ii,Lul)

      id=Aze(lvol,0,ii)%i(ll) ; jd=Ato(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + Wtoze(jj,pp,ii,ll) - Wtoze(jj,pp,ii,Lzl) - Wtoze(jj,Lup,ii,ll) + Wtoze(jj,Lup,ii,Lzl)

      id=Aze(lvol,0,ii)%i(ll) ; jd=Azo(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + Wzoze(jj,pp,ii,ll) - Wzoze(jj,pp,ii,Lzl) - Wzoze(jj,Lzp,ii,ll) + Wzoze(jj,Lzp,ii,Lzl)

      id=Azo(lvol,0,ii)%i(ll) ; jd=Ate(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + Wtezo(jj,pp,ii,ll) - Wtezo(jj,pp,ii,Lzl) - Wtezo(jj,Lup,ii,ll) + Wtezo(jj,Lup,ii,Lzl)

      id=Azo(lvol,0,ii)%i(ll) ; jd=Ato(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + Wtozo(jj,pp,ii,ll) - Wtozo(jj,pp,ii,Lzl) - Wtozo(jj,Lup,ii,ll) + Wtozo(jj,Lup,ii,Lzl)

      id=Azo(lvol,0,ii)%i(ll) ; jd=Aze(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + Wzezo(jj,pp,ii,ll) - Wzezo(jj,pp,ii,Lzl) - Wzezo(jj,Lzp,ii,ll) + Wzezo(jj,Lzp,ii,Lzl)

      id=Azo(lvol,0,ii)%i(ll) ; jd=Azo(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + Wzozo(jj,pp,ii,ll) - Wzozo(jj,pp,ii,Lzl) - Wzozo(jj,Lzp,ii,ll) + Wzozo(jj,Lzp,ii,Lzl)

      id=Ate(lvol,0,ii)%i(ll) ; jd=Ato(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + Htote(jj,pp,ii,ll) - Htote(jj,pp,ii,Lul) - Htote(jj,Lup,ii,ll) + Htote(jj,Lup,ii,Lul)

      id=Ate(lvol,0,ii)%i(ll) ; jd=Azo(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + Hzote(jj,pp,ii,ll) - Hzote(jj,pp,ii,Lul) - Hzote(jj,Lzp,ii,ll) + Hzote(jj,Lzp,ii,Lul)

      id=Ato(lvol,0,ii)%i(ll) ; jd=Ate(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + Hteto(jj,pp,ii,ll) - Hteto(jj,pp,ii,Lul) - Hteto(jj,Lup,ii,ll) + Hteto(jj,Lup,ii,Lul)

      id=Ato(lvol,0,ii)%i(ll) ; jd=Ato(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + Htoto(jj,pp,ii,ll) - Htoto(jj,pp,ii,Lul) - Htoto(jj,Lup,ii,ll) + Htoto(jj,Lup,ii,Lul)

      id=Ato(lvol,0,ii)%i(ll) ; jd=Aze(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + Hzeto(jj,pp,ii,ll) - Hzeto(jj,pp,ii,Lul) - Hzeto(jj,Lzp,ii,ll) + Hzeto(jj,Lzp,ii,Lul)

      id=Ato(lvol,0,ii)%i(ll) ; jd=Azo(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + Hzoto(jj,pp,ii,ll) - Hzoto(jj,pp,ii,Lul) - Hzoto(jj,Lzp,ii,ll) + Hzoto(jj,Lzp,ii,Lul)

      id=Aze(lvol,0,ii)%i(ll) ; jd=Ato(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + Htoze(jj,pp,ii,ll) - Htoze(jj,pp,ii,Lzl) - Htoze(jj,Lup,ii,ll) + Htoze(jj,Lup,ii,Lzl)

      id=Aze(lvol,0,ii)%i(ll) ; jd=Azo(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + Hzoze(jj,pp,ii,ll) - Hzoze(jj,pp,ii,Lzl) - Hzoze(jj,Lzp,ii,ll) + Hzoze(jj,Lzp,ii,Lzl)

      id=Azo(lvol,0,ii)%i(ll) ; jd=Ate(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + Htezo(jj,pp,ii,ll) - Htezo(jj,pp,ii,Lzl) - Htezo(jj,Lup,ii,ll) + Htezo(jj,Lup,ii,Lzl)

      id=Azo(lvol,0,ii)%i(ll) ; jd=Ato(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + Htozo(jj,pp,ii,ll) - Htozo(jj,pp,ii,Lzl) - Htozo(jj,Lup,ii,ll) + Htozo(jj,Lup,ii,Lzl)

      id=Azo(lvol,0,ii)%i(ll) ; jd=Aze(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + Hzezo(jj,pp,ii,ll) - Hzezo(jj,pp,ii,Lzl) - Hzezo(jj,Lzp,ii,ll) + Hzezo(jj,Lzp,ii,Lzl)

      id=Azo(lvol,0,ii)%i(ll) ; jd=Azo(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + Hzozo(jj,pp,ii,ll) - Hzozo(jj,pp,ii,Lzl) - Hzozo(jj,Lzp,ii,ll) + Hzozo(jj,Lzp,ii,Lzl)

      endif

     enddo ! end of do pp; 24 Jan 13;

    enddo ! end of do jj; 24 Jan 13;
    
   enddo ! end of do ll; 24 Jan 13;

  enddo ! end of do ii; 24 Jan 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!#ifdef DEBUG
!  if( Wma01ag ) write(ounit,'("ma01ag : ", 10x ," : myid=",i3," ; lvol=",i3," ; 1.0 ;")') myid, lvol
!#endif

  do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
   
  !if( Lcoordinatesingularity .and. ( mi.ne.0 .or. ni.eq.0 ) ) then ; Lzi = .true. ! 15 Jan 15; ! additional freedom; 15 Jan 15;
   if( Lcoordinatesingularity .and. ( mi.ne.0 .or. ii.eq.1 ) ) then ; Lzi = .true. ! 15 Jan 15; ! additional freedom; 15 Jan 15;
   else                                                             ; Lzi = .false.
   endif
   
   do ll = 0, lrad-1 ! this is a derivative loop;
    
    if( (ll/2)*2.eq.ll ) then ; Lul = lrad   ! ll is even; 18 Jan 13;
    else                      ; Lul = lrad-1 ! ll is odd ; 18 Jan 13;
    endif
    
    if( Lzi ) then ; Lzl = lrad
    else           ; Lzl = Lul
    endif
    
    do jj = 1, mn ; mj = im(jj) ; nj = in(jj)
   
     if( Lcoordinatesingularity .and. ( mj.ne.0 .or. nj.eq.0 ) ) then ; Lzj = .true.
     else                                                             ; Lzj = .false.
     endif
     
     do pp = lrad-1, lrad ! this is a summation loop; 18 Feb 13;

      if( Lzj ) then ; Lzp = lrad
      else           ; Lzp = pp
      endif
      
!#ifdef DEBUG
!      FATALMESS( ma01ag, Ate(lvol,0,ii)%i(ll).lt.0 .or. Ate(lvol,0,ii)%i(ll).gt.NN, invalid subscript )
!      FATALMESS( ma01ag, Aze(lvol,0,ii)%i(ll).lt.0 .or. Aze(lvol,0,ii)%i(ll).gt.NN, invalid subscript )
!      FATALMESS( ma01ag, Fso(lvol,jj)        .lt.0 .or. Fso(lvol,jj)        .gt.NN, invalid subscript )
!      if( NOTstellsym ) then
!      FATALMESS( ma01ag, Ato(lvol,0,ii)%i(ll).lt.0 .or. Ato(lvol,0,ii)%i(ll).gt.NN, invalid subscript )
!      FATALMESS( ma01ag, Azo(lvol,0,ii)%i(ll).lt.0 .or. Azo(lvol,0,ii)%i(ll).gt.NN, invalid subscript )
!      FATALMESS( ma01ag, Fse(lvol,jj)        .lt.0 .or. Fse(lvol,jj)        .gt.NN, invalid subscript )
!      endif
!#endif
      
      id=Ate(lvol,0,ii)%i(ll) ; jd=Fso(lvol,jj)
      dMA(id,jd) = dMA(id,jd) + half * ( + mj*(Wtete(jj,pp,ii,ll)-Wtete(jj,pp,ii,Lul)) - nj*(Wzete(jj,Lzp,ii,ll)-Wzete(jj,Lzp,ii,Lul)) )

      id=Aze(lvol,0,ii)%i(ll) ; jd=Fso(lvol,jj)
      dMA(id,jd) = dMA(id,jd) + half * ( + mj*(Wteze(jj,pp,ii,ll)-Wteze(jj,pp,ii,Lzl)) - nj*(Wzeze(jj,Lzp,ii,ll)-Wzeze(jj,Lzp,ii,Lzl)) )

      id=Ate(lvol,0,ii)%i(ll) ; jd=Fso(lvol,jj)
      dMD(id,jd) = dMD(id,jd) + half * ( + mj*(Htete(jj,pp,ii,ll)-Htete(jj,pp,ii,Lul)) - nj*(Hzete(jj,Lzp,ii,ll)-Hzete(jj,Lzp,ii,Lul)) )

      id=Aze(lvol,0,ii)%i(ll) ; jd=Fso(lvol,jj)
      dMD(id,jd) = dMD(id,jd) + half * ( + mj*(Hteze(jj,pp,ii,ll)-Hteze(jj,pp,ii,Lzl)) - nj*(Hzeze(jj,Lzp,ii,ll)-Hzeze(jj,Lzp,ii,Lzl)) )

      if( NOTstellsym ) then
      id=Ate(lvol,0,ii)%i(ll) ; jd=Fse(lvol,jj)
      dMA(id,jd) = dMA(id,jd) + half * ( - mj*(Wtote(jj,pp,ii,ll)-Wtote(jj,pp,ii,Lul)) + nj*(Wzote(jj,Lzp,ii,ll)-Wzote(jj,Lzp,ii,Lul)) )

      id=Ato(lvol,0,ii)%i(ll) ; jd=Fso(lvol,jj)
      dMA(id,jd) = dMA(id,jd) + half * ( + mj*(Wteto(jj,pp,ii,ll)-Wteto(jj,pp,ii,Lul)) - nj*(Wzeto(jj,Lzp,ii,ll)-Wzeto(jj,Lzp,ii,Lul)) )

      id=Ato(lvol,0,ii)%i(ll) ; jd=Fse(lvol,jj)
      dMA(id,jd) = dMA(id,jd) + half * ( - mj*(Wtoto(jj,pp,ii,ll)-Wtoto(jj,pp,ii,Lul)) + nj*(Wzoto(jj,Lzp,ii,ll)-Wzoto(jj,Lzp,ii,Lul)) )

      id=Aze(lvol,0,ii)%i(ll) ; jd=Fse(lvol,jj)
      dMA(id,jd) = dMA(id,jd) + half * ( - mj*(Wtoze(jj,pp,ii,ll)-Wtoze(jj,pp,ii,Lzl)) + nj*(Wzoze(jj,Lzp,ii,ll)-Wzoze(jj,Lzp,ii,Lzl)) )

      id=Azo(lvol,0,ii)%i(ll) ; jd=Fso(lvol,jj)
      dMA(id,jd) = dMA(id,jd) + half * ( + mj*(Wtezo(jj,pp,ii,ll)-Wtezo(jj,pp,ii,Lzl)) - nj*(Wzezo(jj,Lzp,ii,ll)-Wzezo(jj,Lzp,ii,Lzl)) )

      id=Azo(lvol,0,ii)%i(ll) ; jd=Fse(lvol,jj)
      dMA(id,jd) = dMA(id,jd) + half * ( - mj*(Wtozo(jj,pp,ii,ll)-Wtozo(jj,pp,ii,Lzl)) + nj*(Wzozo(jj,Lzp,ii,ll)-Wzozo(jj,Lzp,ii,Lzl)) )

      id=Ate(lvol,0,ii)%i(ll) ; jd=Fse(lvol,jj)
      dMD(id,jd) = dMD(id,jd) + half * ( - mj*(Htote(jj,pp,ii,ll)-Htote(jj,pp,ii,Lul)) + nj*(Hzote(jj,Lzp,ii,ll)-Hzote(jj,Lzp,ii,Lul)) )

      id=Ato(lvol,0,ii)%i(ll) ; jd=Fso(lvol,jj)
      dMD(id,jd) = dMD(id,jd) + half * ( + mj*(Hteto(jj,pp,ii,ll)-Hteto(jj,pp,ii,Lul)) - nj*(Hzeto(jj,Lzp,ii,ll)-Hzeto(jj,Lzp,ii,Lul)) )

      id=Ato(lvol,0,ii)%i(ll) ; jd=Fse(lvol,jj)
      dMD(id,jd) = dMD(id,jd) + half * ( - mj*(Htoto(jj,pp,ii,ll)-Htoto(jj,pp,ii,Lul)) + nj*(Hzoto(jj,Lzp,ii,ll)-Hzoto(jj,Lzp,ii,Lul)) )

      id=Aze(lvol,0,ii)%i(ll) ; jd=Fse(lvol,jj)
      dMD(id,jd) = dMD(id,jd) + half * ( - mj*(Htoze(jj,pp,ii,ll)-Htoze(jj,pp,ii,Lzl)) + nj*(Hzoze(jj,Lzp,ii,ll)-Hzoze(jj,Lzp,ii,Lzl)) )

      id=Azo(lvol,0,ii)%i(ll) ; jd=Fso(lvol,jj)
      dMD(id,jd) = dMD(id,jd) + half * ( + mj*(Htezo(jj,pp,ii,ll)-Htezo(jj,pp,ii,Lzl)) - nj*(Hzezo(jj,Lzp,ii,ll)-Hzezo(jj,Lzp,ii,Lzl)) )

      id=Azo(lvol,0,ii)%i(ll) ; jd=Fse(lvol,jj)
      dMD(id,jd) = dMD(id,jd) + half * ( - mj*(Htozo(jj,pp,ii,ll)-Htozo(jj,pp,ii,Lzl)) + nj*(Hzozo(jj,Lzp,ii,ll)-Hzozo(jj,Lzp,ii,Lzl)) )

      endif

     enddo ! end of do pp; 24 Jan 13;
     
    enddo ! end of do jj; 24 Jan 13;
    
   enddo ! end of do ll; 24 Jan 13;

  enddo ! end of do ii; 24 Jan 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
!#ifdef DEBUG
!  if( Wma01ag ) write(ounit,'("ma01ag : ", 10x ," : myid=",i3," ; lvol=",i3," ; 2.0 ;")') myid, lvol
!#endif

  do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
   
  !if( Lcoordinatesingularity .and. ( mi.ne.0 .or. ni.eq.0 ) ) then ; Lzi = .true. ! 15 Jan 15; ! additional freedom; 15 Jan 15;
   if( Lcoordinatesingularity .and. ( mi.ne.0 .or. ii.eq.1 ) ) then ; Lzi = .true. ! 15 Jan 15; ! additional freedom; 15 Jan 15;
   else                                                             ; Lzi = .false.
   endif
   
   do ll = lrad-1, lrad ! this is a summation loop; 18 Feb 13;
    
    if( Lzi ) then ; Lzl = lrad
    else           ; Lzl = ll
    endif
    
    do jj = 1, mn ; mj = im(jj) ; nj = in(jj)
   
     if( Lcoordinatesingularity .and. ( mj.ne.0 .or. nj.eq.0 ) ) then ; Lzj = .true.
     else                                                             ; Lzj = .false.
     endif
     
     do pp = 0, lrad-1 ! this is a derivative loop;
      
      if( (pp/2)*2.eq.pp ) then ; Lup = lrad   ! ll is even; 18 Jan 13;
      else                      ; Lup = lrad-1 ! ll is odd ; 18 Jan 13;
      endif

      if( Lzj ) then ; Lzp = lrad
      else           ; Lzp = Lup
      endif
      
!#ifdef DEBUG
!      FATALMESS( ma01ag, Ate(lvol,0,jj)%i(pp).lt.0 .or. Ate(lvol,0,jj)%i(pp).gt.NN, invalid subscript )
!      FATALMESS( ma01ag, Aze(lvol,0,jj)%i(pp).lt.0 .or. Aze(lvol,0,jj)%i(pp).gt.NN, invalid subscript )
!      FATALMESS( ma01ag, Fso(lvol,ii)        .lt.0 .or. Fso(lvol,ii)        .gt.NN, invalid subscript )
!      if( NOTstellsym ) then
!      FATALMESS( ma01ag, Ato(lvol,0,jj)%i(pp).lt.0 .or. Ato(lvol,0,jj)%i(pp).gt.NN, invalid subscript )
!      FATALMESS( ma01ag, Azo(lvol,0,jj)%i(pp).lt.0 .or. Azo(lvol,0,jj)%i(pp).gt.NN, invalid subscript )
!      FATALMESS( ma01ag, Fse(lvol,ii)        .lt.0 .or. Fse(lvol,ii)        .gt.NN, invalid subscript )
!      endif
!#endif

      id=Fso(lvol,ii) ; jd=Ate(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + half * ( + mi*(Wtete(jj,pp,ii,ll)-Wtete(jj,Lup,ii,ll)) - ni*(Wteze(jj,pp,ii,Lzl)-Wteze(jj,Lup,ii,Lzl)) )

      id=Fso(lvol,ii) ; jd=Aze(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + half * ( + mi*(Wzete(jj,pp,ii,ll)-Wzete(jj,Lzp,ii,ll)) - ni*(Wzeze(jj,pp,ii,Lzl)-Wzeze(jj,Lzp,ii,Lzl)) )

      id=Fso(lvol,ii) ; jd=Ate(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + half * ( + mi*(Htete(jj,pp,ii,ll)-Htete(jj,Lup,ii,ll)) - ni*(Hteze(jj,pp,ii,Lzl)-Hteze(jj,Lup,ii,Lzl)) )

      id=Fso(lvol,ii) ; jd=Aze(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + half * ( + mi*(Hzete(jj,pp,ii,ll)-Hzete(jj,Lzp,ii,ll)) - ni*(Hzeze(jj,pp,ii,Lzl)-Hzeze(jj,Lzp,ii,Lzl)) )

      if( NOTstellsym ) then

      id=Fso(lvol,ii) ; jd=Ato(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + half * ( + mi*(Wtote(jj,pp,ii,ll)-Wtote(jj,Lup,ii,ll)) - ni*(Wtoze(jj,pp,ii,Lzl)-Wtoze(jj,Lup,ii,Lzl)) )

      id=Fso(lvol,ii) ; jd=Azo(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + half * ( + mi*(Wzote(jj,pp,ii,ll)-Wzote(jj,Lzp,ii,ll)) - ni*(Wzoze(jj,pp,ii,Lzl)-Wzoze(jj,Lzp,ii,Lzl)) )

      id=Fse(lvol,ii) ; jd=Ate(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + half * ( - mi*(Wteto(jj,pp,ii,ll)-Wteto(jj,Lup,ii,ll)) + ni*(Wtezo(jj,pp,ii,Lzl)-Wtezo(jj,Lup,ii,Lzl)) )

      id=Fse(lvol,ii) ; jd=Ato(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + half * ( - mi*(Wtoto(jj,pp,ii,ll)-Wtoto(jj,Lup,ii,ll)) + ni*(Wtozo(jj,pp,ii,Lzl)-Wtozo(jj,Lup,ii,Lzl)) )

      id=Fse(lvol,ii) ; jd=Aze(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + half * ( - mi*(Wzeto(jj,pp,ii,ll)-Wzeto(jj,Lzp,ii,ll)) + ni*(Wzezo(jj,pp,ii,Lzl)-Wzezo(jj,Lzp,ii,Lzl)) )

      id=Fse(lvol,ii) ; jd=Azo(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + half * ( - mi*(Wzoto(jj,pp,ii,ll)-Wzoto(jj,Lzp,ii,ll)) + ni*(Wzozo(jj,pp,ii,Lzl)-Wzozo(jj,Lzp,ii,Lzl)) )

      id=Fso(lvol,ii) ; jd=Ato(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + half * ( + mi*(Htote(jj,pp,ii,ll)-Htote(jj,Lup,ii,ll)) - ni*(Htoze(jj,pp,ii,Lzl)-Htoze(jj,Lup,ii,Lzl)) )

      id=Fso(lvol,ii) ; jd=Azo(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + half * ( + mi*(Hzote(jj,pp,ii,ll)-Hzote(jj,Lzp,ii,ll)) - ni*(Hzoze(jj,pp,ii,Lzl)-Hzoze(jj,Lzp,ii,Lzl)) )

      id=Fse(lvol,ii) ; jd=Ate(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + half * ( - mi*(Hteto(jj,pp,ii,ll)-Hteto(jj,Lup,ii,ll)) + ni*(Htezo(jj,pp,ii,Lzl)-Htezo(jj,Lup,ii,Lzl)) )

      id=Fse(lvol,ii) ; jd=Ato(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + half * ( - mi*(Htoto(jj,pp,ii,ll)-Htoto(jj,Lup,ii,ll)) + ni*(Htozo(jj,pp,ii,Lzl)-Htozo(jj,Lup,ii,Lzl)) )

      id=Fse(lvol,ii) ; jd=Aze(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + half * ( - mi*(Hzeto(jj,pp,ii,ll)-Hzeto(jj,Lzp,ii,ll)) + ni*(Hzezo(jj,pp,ii,Lzl)-Hzezo(jj,Lzp,ii,Lzl)) )

      id=Fse(lvol,ii) ; jd=Azo(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + half * ( - mi*(Hzoto(jj,pp,ii,ll)-Hzoto(jj,Lzp,ii,ll)) + ni*(Hzozo(jj,pp,ii,Lzl)-Hzozo(jj,Lzp,ii,Lzl)) )

      endif

     enddo ! end of do pp; 24 Jan 13;

    enddo ! end of do jj; 24 Jan 13;
    
   enddo ! end of do ll; 24 Jan 13;

  enddo ! end of do ii; 24 Jan 13;
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!#ifdef DEBUG
!  if( Wma01ag ) write(ounit,'("ma01ag : ", 10x ," : myid=",i3," ; lvol=",i3," ; 3.0 ;")') myid, lvol
!#endif

  do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
   
  !if( Lcoordinatesingularity .and. ( mi.ne.0 .or. ni.eq.0 ) ) then ; Lzi = .true. ! 15 Jan 15; ! additional freedom; 15 Jan 15;
   if( Lcoordinatesingularity .and. ( mi.ne.0 .or. ii.eq.1 ) ) then ; Lzi = .true. ! 15 Jan 15; ! additional freedom; 15 Jan 15;
   else                                                             ; Lzi = .false.
   endif
   
   do ll = lrad-1, lrad ! this is a summation loop; 18 Feb 13;
    
    if( Lzi ) then ; Lzl = lrad
    else           ; Lzl = ll
    endif
    
    do jj = 1, mn ; mj = im(jj) ; nj = in(jj)
   
     if( Lcoordinatesingularity .and. ( mj.ne.0 .or. nj.eq.0 ) ) then ; Lzj = .true.
     else                                                             ; Lzj = .false.
     endif
     
     mjmi = mj * mi ; mjni = mj * ni ; njmi = nj * mi ; njni = nj * ni
     
     do pp = lrad-1, lrad ! this is a summation loop; 18 Feb 13;
    
      if( Lzj ) then ; Lzp = lrad
      else           ; Lzp = pp
      endif

!#ifdef DEBUG
!      FATALMESS( ma01ag, Fso(lvol,ii)        .lt.0 .or. Fso(lvol,ii)        .gt.NN, invalid subscript )
!      FATALMESS( ma01ag, Fso(lvol,jj)        .lt.0 .or. Fso(lvol,jj)        .gt.NN, invalid subscript )
!      if( NOTstellsym ) then
!      FATALMESS( ma01ag, Fse(lvol,ii)        .lt.0 .or. Fse(lvol,ii)        .gt.NN, invalid subscript )
!      FATALMESS( ma01ag, Fse(lvol,jj)        .lt.0 .or. Fse(lvol,jj)        .gt.NN, invalid subscript )
!      endif
!#endif

      id=Fso(lvol,ii) ; jd=Fso(lvol,jj)
      dMA(id,jd) = dMA(id,jd) + quart*( +mjmi*Wtete(jj,pp,ii,ll)-mjni*Wteze(jj,pp,ii,Lzl)-njmi*Wzete(jj,Lzp,ii,ll)+njni*Wzeze(jj,Lzp,ii,Lzl))

      id=Fso(lvol,ii) ; jd=Fso(lvol,jj)
      dMD(id,jd) = dMD(id,jd) + quart*( +mjmi*Htete(jj,pp,ii,ll)-mjni*Hteze(jj,pp,ii,Lzl)-njmi*Hzete(jj,Lzp,ii,ll)+njni*Hzeze(jj,Lzp,ii,Lzl))

      if( NOTstellsym ) then

      id=Fso(lvol,ii) ; jd=Fse(lvol,jj)
      dMA(id,jd) = dMA(id,jd) + quart*( -mjmi*Wtote(jj,pp,ii,ll)+mjni*Wtoze(jj,pp,ii,Lzl)+njmi*Wzote(jj,Lzp,ii,ll)-njni*Wzoze(jj,Lzp,ii,Lzl))

      id=Fse(lvol,ii) ; jd=Fso(lvol,jj)
      dMA(id,jd) = dMA(id,jd) + quart*( -mjmi*Wteto(jj,pp,ii,ll)+mjni*Wtezo(jj,pp,ii,Lzl)+njmi*Wzeto(jj,Lzp,ii,ll)-njni*Wzezo(jj,Lzp,ii,Lzl))

      id=Fse(lvol,ii) ; jd=Fse(lvol,jj)
      dMA(id,jd) = dMA(id,jd) + quart*( +mjmi*Wtoto(jj,pp,ii,ll)-mjni*Wtozo(jj,pp,ii,Lzl)-njmi*Wzoto(jj,Lzp,ii,ll)+njni*Wzozo(jj,Lzp,ii,Lzl))

      id=Fso(lvol,ii) ; jd=Fse(lvol,jj)
      dMD(id,jd) = dMD(id,jd) + quart*( -mjmi*Htote(jj,pp,ii,ll)+mjni*Htoze(jj,pp,ii,Lzl)+njmi*Hzote(jj,Lzp,ii,ll)-njni*Hzoze(jj,Lzp,ii,Lzl))

      id=Fse(lvol,ii) ; jd=Fso(lvol,jj)
      dMD(id,jd) = dMD(id,jd) + quart*( -mjmi*Hteto(jj,pp,ii,ll)+mjni*Htezo(jj,pp,ii,Lzl)+njmi*Hzeto(jj,Lzp,ii,ll)-njni*Hzezo(jj,Lzp,ii,Lzl))

      id=Fse(lvol,ii) ; jd=Fse(lvol,jj)
      dMD(id,jd) = dMD(id,jd) + quart*( +mjmi*Htoto(jj,pp,ii,ll)-mjni*Htozo(jj,pp,ii,Lzl)-njmi*Hzoto(jj,Lzp,ii,ll)+njni*Hzozo(jj,Lzp,ii,Lzl))

      endif

     enddo ! end of do pp; 24 Jan 13;
     
    enddo ! end of do jj; 24 Jan 13;
    
   enddo ! end of do ll; 24 Jan 13;
   
  enddo ! end of do ii; 24 Jan 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!#ifdef DEBUG
!  if( Wma01ag ) write(ounit,'("ma01ag : ", 10x ," : myid=",i3," ; lvol=",i3," ; 4.0 ;")') myid, lvol
!#endif

  do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
   
  !if( Lcoordinatesingularity .and. ( mi.ne.0 .or. ni.eq.0 ) ) then ; Lzi = .true. ! 15 Jan 15; ! additional freedom; 15 Jan 15;
   if( Lcoordinatesingularity .and. ( mi.ne.0 .or. ii.eq.1 ) ) then ; Lzi = .true. ! 15 Jan 15; ! additional freedom; 15 Jan 15;
   else                                                             ; Lzi = .false.
   endif
   
   do ll = 0, lrad-1 ! this is a derivative loop;

    if( (ll/2)*2.eq.ll ) then ; Lul = lrad   ! ll is even; 18 Jan 13;
    else                      ; Lul = lrad-1 ! ll is odd ; 18 Jan 13;
    endif
    
    if( Lzi ) then ; Lzl = lrad
    else           ; Lzl = Lul
    endif
    
    jj = 1 ; mj = im(jj) ; nj = in(jj)
     
     if( Lcoordinatesingularity .and. ( mj.ne.0 .or. nj.eq.0 ) ) then ; Lzj = .true.
     else                                                             ; Lzj = .false.
     endif
     
     do pp = lrad-1, lrad ! this is a summation loop; 18 Feb 13;
      
      if( Lzj ) then ; Lzp = lrad
      else           ; Lzp = pp
      endif
      
!#ifdef DEBUG
!      FATALMESS( ma01ag, Ate(lvol,0,ii)%i(ll).lt.0 .or. Ate(lvol,0,ii)%i(ll).gt.NN, invalid subscript )
!      FATALMESS( ma01ag, Aze(lvol,0,ii)%i(ll).lt.0 .or. Aze(lvol,0,ii)%i(ll).gt.NN, invalid subscript )
!      if( NOTstellsym ) then
!      FATALMESS( ma01ag, Ato(lvol,0,ii)%i(ll).lt.0 .or. Ato(lvol,0,ii)%i(ll).gt.NN, invalid subscript )
!      FATALMESS( ma01ag, Azo(lvol,0,ii)%i(ll).lt.0 .or. Azo(lvol,0,ii)%i(ll).gt.NN, invalid subscript )
!      endif
!#endif

      id=Ate(lvol,0,ii)%i(ll) ; jd=1 ; dMB(id,jd) = dMB(id,jd) + half * ( Wtete(jj, pp,ii,ll) - Wtete(jj, pp,ii,Lul) ) !
      id=Aze(lvol,0,ii)%i(ll) ; jd=1 ; dMB(id,jd) = dMB(id,jd) + half * ( Wteze(jj, pp,ii,ll) - Wteze(jj, pp,ii,Lzl) ) ! Lzi;
      id=Ate(lvol,0,ii)%i(ll) ; jd=2 ; dMB(id,jd) = dMB(id,jd) + half * ( Wzete(jj,Lzp,ii,ll) - Wzete(jj,Lzp,ii,Lul) ) !
      id=Aze(lvol,0,ii)%i(ll) ; jd=2 ; dMB(id,jd) = dMB(id,jd) + half * ( Wzeze(jj,Lzp,ii,ll) - Wzeze(jj,Lzp,ii,Lzl) ) ! Lzi;
      id=Ate(lvol,0,ii)%i(ll) ; jd=1 ; dME(id,jd) = dME(id,jd) + half * ( Htete(jj, pp,ii,ll) - Htete(jj, pp,ii,Lul) ) !
      id=Aze(lvol,0,ii)%i(ll) ; jd=1 ; dME(id,jd) = dME(id,jd) + half * ( Hteze(jj, pp,ii,ll) - Hteze(jj, pp,ii,Lzl) ) ! Lzi;
      id=Ate(lvol,0,ii)%i(ll) ; jd=2 ; dME(id,jd) = dME(id,jd) + half * ( Hzete(jj,Lzp,ii,ll) - Hzete(jj,Lzp,ii,Lul) ) !
      id=Aze(lvol,0,ii)%i(ll) ; jd=2 ; dME(id,jd) = dME(id,jd) + half * ( Hzeze(jj,Lzp,ii,ll) - Hzeze(jj,Lzp,ii,Lzl) ) ! Lzi;
      if( NOTstellsym ) then
      id=Ato(lvol,0,ii)%i(ll) ; jd=1 ; dMB(id,jd) = dMB(id,jd) + half * ( Wteto(jj, pp,ii,ll) - Wteto(jj, pp,ii,Lul) ) !
      id=Azo(lvol,0,ii)%i(ll) ; jd=1 ; dMB(id,jd) = dMB(id,jd) + half * ( Wtezo(jj, pp,ii,ll) - Wtezo(jj, pp,ii,Lzl) ) ! Lzi;
      id=Ato(lvol,0,ii)%i(ll) ; jd=2 ; dMB(id,jd) = dMB(id,jd) + half * ( Wzeto(jj,Lzp,ii,ll) - Wzeto(jj,Lzp,ii,Lul) ) !
      id=Azo(lvol,0,ii)%i(ll) ; jd=2 ; dMB(id,jd) = dMB(id,jd) + half * ( Wzezo(jj,Lzp,ii,ll) - Wzezo(jj,Lzp,ii,Lzl) ) ! Lzi;
      id=Ato(lvol,0,ii)%i(ll) ; jd=1 ; dME(id,jd) = dME(id,jd) + half * ( Hteto(jj, pp,ii,ll) - Hteto(jj, pp,ii,Lul) ) !
      id=Azo(lvol,0,ii)%i(ll) ; jd=1 ; dME(id,jd) = dME(id,jd) + half * ( Htezo(jj, pp,ii,ll) - Htezo(jj, pp,ii,Lzl) ) ! Lzi;
      id=Ato(lvol,0,ii)%i(ll) ; jd=2 ; dME(id,jd) = dME(id,jd) + half * ( Hzeto(jj,Lzp,ii,ll) - Hzeto(jj,Lzp,ii,Lul) ) !
      id=Azo(lvol,0,ii)%i(ll) ; jd=2 ; dME(id,jd) = dME(id,jd) + half * ( Hzezo(jj,Lzp,ii,ll) - Hzezo(jj,Lzp,ii,Lzl) ) ! Lzi;
      endif

     enddo ! end of do pp; 24 Jan 13;
   
   !enddo ! end of do jj; 26 Feb 13;

   enddo ! end of do ll; 24 Jan 13;

  enddo ! end of do ii; 24 Jan 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!#ifdef DEBUG
!  if( Wma01ag ) write(ounit,'("ma01ag : ", 10x ," : myid=",i3," ; lvol=",i3," ; 5.0 ;")') myid, lvol
!#endif

  do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
   
  !if( Lcoordinatesingularity .and. ( mi.ne.0 .or. ni.eq.0 ) ) then ; Lzi = .true. ! 15 Jan 15; ! additional freedom; 15 Jan 15;
   if( Lcoordinatesingularity .and. ( mi.ne.0 .or. ii.eq.1 ) ) then ; Lzi = .true. ! 15 Jan 15; ! additional freedom; 15 Jan 15;
   else                                                             ; Lzi = .false.
   endif
   
   do ll = lrad-1, lrad ! this is a summation loop; 18 Feb 13;
    
    if( Lzi ) then ; Lzl = lrad
    else           ; Lzl = ll
    endif
    
    jj = 1 ; mj = im(jj) ; nj = in(jj)
   
     if( Lcoordinatesingularity .and. ( mj.ne.0 .or. nj.eq.0 ) ) then ; Lzj = .true.
     else                                                             ; Lzj = .false.
     endif
    
     do pp = lrad-1, lrad ! this is a summation loop; 18 Feb 13;
    
      if( Lzj ) then ; Lzp = lrad
      else           ; Lzp = pp
      endif

!#ifdef DEBUG
!      FATALMESS( ma01ag, Fso(lvol,ii)        .lt.0 .or. Fso(lvol,ii)        .gt.NN, invalid subscript )
!      if( NOTstellsym ) then
!      FATALMESS( ma01ag, Fse(lvol,ii)        .lt.0 .or. Fse(lvol,ii)        .gt.NN, invalid subscript )
!      endif
!#endif

      id=Fso(lvol,ii) ; jd=1 ; dMB(id,jd) = dMB(id,jd) + quart * ( mi * Wtete(jj, pp,ii,ll) - ni * Wteze(jj, pp,ii,Lzl) )
      id=Fso(lvol,ii) ; jd=2 ; dMB(id,jd) = dMB(id,jd) + quart * ( mi * Wzete(jj,Lzp,ii,ll) - ni * Wzeze(jj,Lzp,ii,Lzl) )
      id=Fso(lvol,ii) ; jd=1 ; dME(id,jd) = dME(id,jd) + quart * ( mi * Htete(jj, pp,ii,ll) - ni * Hteze(jj, pp,ii,Lzl) )
      id=Fso(lvol,ii) ; jd=2 ; dME(id,jd) = dME(id,jd) + quart * ( mi * Hzete(jj,Lzp,ii,ll) - ni * Hzeze(jj,Lzp,ii,Lzl) )
      if( NOTstellsym ) then
      id=Fse(lvol,ii) ; jd=1 ; dMB(id,jd) = dMB(id,jd) - quart * ( mi * Wteto(jj, pp,ii,ll) - ni * Wtezo(jj, pp,ii,Lzl) )
      id=Fse(lvol,ii) ; jd=2 ; dMB(id,jd) = dMB(id,jd) - quart * ( mi * Wzeto(jj,Lzp,ii,ll) - ni * Wzezo(jj,Lzp,ii,Lzl) )
      id=Fse(lvol,ii) ; jd=1 ; dME(id,jd) = dME(id,jd) - quart * ( mi * Hteto(jj, pp,ii,ll) - ni * Htezo(jj, pp,ii,Lzl) )
      id=Fse(lvol,ii) ; jd=2 ; dME(id,jd) = dME(id,jd) - quart * ( mi * Hzeto(jj,Lzp,ii,ll) - ni * Hzezo(jj,Lzp,ii,Lzl) )
      endif

     enddo ! end of do pp; 24 Jan 13;
   
!   end of jj;

   enddo ! end of do ll; 24 Jan 13;

  enddo ! end of do ii; 24 Jan 13;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!#ifdef DEBUG
!  if( Wma01ag ) write(ounit,'("ma01ag : ", 10x ," : myid=",i3," ; lvol=",i3," ; 6.0 ;")') myid, lvol
!#endif

  ii = 1 ; mi = im(ii) ; ni = in(ii)
  
  !if( Lcoordinatesingularity .and. ( mi.ne.0 .or. ni.eq.0 ) ) then ; Lzi = .true. ! 15 Jan 15; ! additional freedom; 15 Jan 15;
   if( Lcoordinatesingularity .and. ( mi.ne.0 .or. ii.eq.1 ) ) then ; Lzi = .true. ! 15 Jan 15; ! additional freedom; 15 Jan 15;
   else                                                             ; Lzi = .false.
   endif
   
   do ll = lrad-1, lrad ! this is a summation loop; 18 Feb 13;
    
    if( Lzi ) then ; Lzl = lrad
    else           ; Lzl = ll
    endif

    jj = 1 ; mj = im(jj) ; nj = in(jj)
   
     if( Lcoordinatesingularity .and. ( mj.ne.0 .or. nj.eq.0 ) ) then ; Lzj = .true.
     else                                                             ; Lzj = .false.
     endif
    
     do pp = lrad-1, lrad ! this is a summation loop; 18 Feb 13;
    
      if( Lzj ) then ; Lzp = lrad
      else           ; Lzp = pp
      endif
      
      id=1 ; jd=1 ; dMC(id,jd) = dMC(id,jd) + quart * Wtete(jj, pp,ii, ll)
      id=1 ; jd=2 ; dMC(id,jd) = dMC(id,jd) + quart * Wzete(jj,Lzp,ii, ll)
      id=2 ; jd=1 ; dMC(id,jd) = dMC(id,jd) + quart * Wteze(jj, pp,ii,Lzl)
      id=2 ; jd=2 ; dMC(id,jd) = dMC(id,jd) + quart * Wzeze(jj,Lzp,ii,Lzl) 
      id=1 ; jd=1 ; dMF(id,jd) = dMF(id,jd) + quart * Htete(jj, pp,ii, ll)
      id=1 ; jd=2 ; dMF(id,jd) = dMF(id,jd) + quart * Hzete(jj,Lzp,ii, ll)
      id=2 ; jd=1 ; dMF(id,jd) = dMF(id,jd) + quart * Hteze(jj, pp,ii,Lzl)
      id=2 ; jd=2 ; dMF(id,jd) = dMF(id,jd) + quart * Hzeze(jj,Lzp,ii,Lzl)
      
     enddo ! end of do pp; 24 Jan 13;
    
!   end of jj;

   enddo ! end of do ll; 24 Jan 13;

! end of ii

  RETURN(ma01ag)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine ma01ag

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

