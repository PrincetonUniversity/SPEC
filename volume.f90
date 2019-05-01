!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title ! Calculates volume of each region; ${\cal V}_i \equiv \int dv$.

!latex \briefly{Computes volume of each region; and, if required, the derivatives of the volume with respect to the interface geometry.}

!latex \calledby{\link{dforce} and \link{xspech}}
!latex \calls{\link{}}

!latex \tableofcontents

!latex \subsection{volume integral} \begin{enumerate}

!latex \item The volume enclosed by the $v$-th interface is given by the integral
!latex       \be V = \int_{{\cal V}} \; dv 
!latex             = \frac{1}{3}\int_{{\cal V}} \; \nabla \cdot {\bf x} \; dv 
!latex             = \frac{1}{3}\int_{{\cal S}} \; {\bf x} \cdot d{\bf s} 
!latex             = \frac{1}{3}        \int_{0}^{2\pi} \!\! d\t \int_{0}^{2\pi/N} \!\! d\z \;\;\; 
!latex                                        \left. {\bf x   } \cdot {\bf x_\t} \times {\bf x_\z} \right|^s
!latex       \ee
!latex       where we have used $\nabla \cdot {\bf x} = 3$, and have assumed that the domain is periodic in the angles.

!latex \end{enumerate} \subsection{representation of surfaces} \begin{enumerate}

!latex \item The coordinate functions are
!latex       \be            R(\t,\z) & = & \sum_i R_{e,i} \; \cos\a_i + \sum_i R_{o,i} \; \sin\a_i \\
!latex                      Z(\t,\z) & = & \sum_i Z_{e,i} \; \cos\a_i + \sum_i Z_{o,i} \; \sin\a_i,
!latex       \ee
!latex       where $\a_i \equiv m_i \t - n_i \z$.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine volume( lvol, vflag )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  use constants, only : zero, half, one, two, third, quart, pi2 

  use numerical, only : vsmall, small

  use fileunits, only : ounit

  use inputlist, only : Wvolume, Igeometry, Nvol, pscale

  use cputiming

  use allglobal, only : myid, cpus, &
                        YESstellsym, Mvol, &
                        mn, im, in, iRbc, iZbs, iRbs, iZbc, &
                        djkp, djkm, &
                        vvolume, dvolume, &
                        dBdX, &
                        pi2nfp, pi2pi2nfp, pi2pi2nfpquart

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS

  INTEGER, intent(in) :: lvol

  INTEGER             :: vflag

  INTEGER             :: jvol, ii, jj, kk, mi, ni, mj, nj, mk, nk, innout

  REAL                :: vol(0:1)

  REAL                :: Rei, Roi, Zei, Zoi, Rej, Roj, Zej, Zoj, Rek, Rok, Zek, Zok

  REAL                :: AA, BB, CC, DD
  
  BEGIN(volume)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef DEBUG
  FATAL( volume, lvol.lt.1 .or. lvol.gt.Mvol, invalid volume ) ! 15 Jan 13;
#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( lvol.gt.Nvol ) then ; vvolume(lvol) = one ; dvolume = zero ; goto 9999 ! this can only be the vacuum region; provide default value; 13 Sep 13;
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  ;  vol(0:1) = zero ! initialize;
  ; dvolume   = zero ! derivatives of volume wrt interface inner/outer, R/Z harmonic;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  do innout = 0, 1 ! will subtract inner volume from outer volume to obtain enclosed volume; 26 Feb 13; ! this does seem a little wasteful; 13 Sep 13;
   
   jvol = lvol - 1 + innout ! labels inner or outer interface; 13 Sep 13;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
!latex \end{enumerate} \subsection{geometry} \begin{enumerate}
   
!latex \item The geometry is controlled by the input parameter \verb+Igeometry+ as follows: 
   
   select case( Igeometry )
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
   case( 1 ) !latex \item \verb+Igeometry.eq.1+ : Cartesian : $\sqrt g = R_s$
    
!latex       \be V & = & \int_{0}^{2\pi}\!\!\!d\t \int_{0}^{2\pi/N}\!\!\!\!\! d\z \; R \nonumber \\
!latex             & = & 2\pi \; \frac{2\pi}{N} \; R_{e,1}
!latex       \ee
    
    vol(innout) = iRbc(1,jvol) ! 20 Jun 14;
    
#ifdef DEBUG
    FATAL( volume, dBdX%L .and. dBdX%irz.eq.1, volume does not depend on Z )
#endif
    
    if( dBdX%L .and. dBdX%innout.eq.innout .and. dBdX%ii.eq.1 ) then ! compute derivative of volume;
     if( dBdX%issym.eq.0 ) dvolume = one ! note that the sign factor for the lower interface is included below; 20 Jun 14;
    endif
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
   case( 2 ) !latex \item \verb+Igeometry.eq.2+ : cylindrical : $\sqrt g = R R_s = \frac{1}{2}\partial_s (R^2)$
    
!latex       \be V & = & \frac{1}{2}\int_{0}^{2\pi}\!\!\!d\t \int_{0}^{2\pi/N}\!\!\!\!\! d\z \; R^2 \nonumber \\
!latex             & = & \frac{1}{2} \; 2\pi \; \frac{2\pi}{N} \; \frac{1}{2} \; 
!latex                   \sum_i\sum_j R_{e,i}R_{e,j}\left[\cos(\a_i-\a_j)+\cos(\a_i+\a_j)\right] \nonumber \\
!latex             & + & \frac{1}{2} \; 2\pi \; \frac{2\pi}{N} \; \frac{1}{2} \;
!latex                   \sum_i\sum_j R_{o,i}R_{o,j}\left[\cos(\a_i-\a_j)-\cos(\a_i+\a_j)\right]
!latex       \ee
    
#ifdef DEBUG
    FATAL( volume, dBdX%L .and. dBdX%irz.eq.1, volume does not depend on Z for cylindrical geometry )
#endif
    
    if( YESstellsym ) then
     
     do ii = 1, mn ; mi = im(ii) ; ni = in(ii)

      do jj = 1, mn ; mj = im(jj) ; nj = in(jj)
     
       vol(innout) = vol(innout) + iRbc(ii,jvol) * iRbc(jj,jvol) * ( djkp(ii,jj) + djkm(ii,jj) )
  
      !if( mi-mj.eq.0 .and. ni-nj.eq.0 ) vol(innout) = vol(innout) + iRbc(ii,jvol) * iRbc(jj,jvol)
      !if( mi+mj.eq.0 .and. ni+nj.eq.0 ) vol(innout) = vol(innout) + iRbc(ii,jvol) * iRbc(jj,jvol)

       if( dBdX%L .and. dBdX%innout.eq.innout .and. dBdX%ii.eq.ii ) then ! compute derivative of volume;
        dvolume = dvolume + iRbc(jj,jvol) * ( djkp(jj,ii) + djkm(jj,ii) + djkp(ii,jj) + djkm(ii,jj) )
       endif

      enddo ! end of do jj; 02 Sep 14;
      
     enddo ! end of do ii; 02 Sep 14;

    else ! NOTstellsym;

     do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
      do jj = 1, mn ; mj = im(jj) ; nj = in(jj)
       
       vol(innout) = vol(innout) + iRbc(ii,jvol) * iRbc(jj,jvol) * ( djkp(ii,jj) + djkm(ii,jj) ) &
                                 + iRbs(ii,jvol) * iRbs(jj,jvol) * ( djkp(ii,jj) - djkm(ii,jj) )

      !if( mi-mj.eq.0 .and. ni-nj.eq.0 ) vol(innout) = vol(innout) + iRbc(ii,jvol) * iRbc(jj,jvol) + iRbs(ii,jvol) * iRbs(jj,jvol)
      !if( mi+mj.eq.0 .and. ni+nj.eq.0 ) vol(innout) = vol(innout) + iRbc(ii,jvol) * iRbc(jj,jvol) - iRbs(ii,jvol) * iRbs(jj,jvol)
       
       if( dBdX%L .and. dBdX%innout.eq.innout .and. dBdX%ii.eq.ii ) then ! compute derivative of volume;
        if( dBdX%issym.eq.0 ) then !     stellarator-symmetric harmonic; dV/dRei ; 13 Sep 13;
        dvolume = dvolume + iRbc(jj,jvol) * ( djkp(jj,ii) + djkm(jj,ii) + djkp(ii,jj) + djkm(ii,jj) )
        else
        FATAL( volume, .true., derivatives of volume under construction )
        dvolume = dvolume + iRbs(jj,jvol) * ( djkp(jj,ii) - djkm(jj,ii) + djkp(ii,jj) - djkm(ii,jj) ) ! needs to be checked; 02 Sep 14;
        endif
       endif
       
      enddo
     enddo

    endif ! end of if( YESstellsym ) ; 11 Aug 14;

!   FATAL( volume, dBdX%L, have not yet computed derivatives of volume wrt interface harmonics for cylindrical geometry )
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
   case( 3 ) !latex \item \verb+Igeometry.eq.3+ : toroidal : ${\bf x}\cdot {\bf e}_\theta \times {\bf e}_\zeta  = R ( Z R_\t - R Z_\t )$
    
!latex       \be V & = & \frac{1}{3} \; \int_{0}^{2\pi}\!\!\!d\t \int_{0}^{2\pi/N}\!\!\!\!\! d\z \; R \left( Z R_{\t} - R Z_{\t}  \right)              
!latex             \nonumber \\
!latex             & = & \frac{1}{3}  \; \sum_i \sum_j \sum_k R_{e,i} \left(Z_{e,j} R_{o,k} - R_{e,j} Z_{o,k} \right) (+m_k) 
!latex             \int \!\!\!\! \int \!\! d\t d\z \; \cos\a_i \cos\a_j \cos\a_k \nonumber \\
!latex             & + & \frac{1}{3}  \; \sum_i \sum_j \sum_k R_{e,i} \left(Z_{o,j} R_{e,k} - R_{o,j} Z_{e,k} \right) (-m_k) 
!latex             \int \!\!\!\! \int \!\! d\t d\z \; \cos\a_i \sin\a_j \sin\a_k \nonumber \\
!latex             & + & \frac{1}{3}  \; \sum_i \sum_j \sum_k R_{o,i} \left(Z_{e,j} R_{e,k} - R_{e,j} Z_{e,k} \right) (-m_k) 
!latex             \int \!\!\!\! \int \!\! d\t d\z \; \sin\a_i \cos\a_j \sin\a_k \nonumber \\
!latex             & + & \frac{1}{3}  \; \sum_i \sum_j \sum_k R_{o,i} \left(Z_{o,j} R_{o,k} - R_{o,j} Z_{o,k} \right) (+m_k) 
!latex             \int \!\!\!\! \int \!\! d\t d\z \; \sin\a_i \sin\a_j \cos\a_k
!latex       \ee

!latex \item (Recall that the integral over an odd function is zero, so various terms in the above expansion have been ignored.)

!latex \item The trigonometric terms are
!latex       \be \begin{array}{ccccccccccccccccccccccccccccccccccccccccccc}
!latex           4 \; \cos\a_i \cos\a_j \cos\a_k & \!\!\! = & \!\!\!
!latex       + \!\!\! & \cos(\a_i+\a_j+\a_k) & \!\!\! + \!\!\! & \cos(\a_i+\a_j-\a_k) & \!\!\! + \!\!\! & \cos(\a_i-\a_j+\a_k) & \!\!\! + \!\!\! & \cos(\a_i-\a_j-\a_k) \\
!latex           4 \; \cos\a_i \sin\a_j \sin\a_k & \!\!\! = & \!\!\!
!latex       - \!\!\! & \cos(\a_i+\a_j+\a_k) & \!\!\! + \!\!\! & \cos(\a_i+\a_j-\a_k) & \!\!\! + \!\!\! & \cos(\a_i-\a_j+\a_k) & \!\!\! - \!\!\! & \cos(\a_i-\a_j-\a_k) \\
!latex           4 \; \sin\a_i \cos\a_j \sin\a_k & \!\!\! = & \!\!\!
!latex       - \!\!\! & \cos(\a_i+\a_j+\a_k) & \!\!\! + \!\!\! & \cos(\a_i+\a_j-\a_k) & \!\!\! - \!\!\! & \cos(\a_i-\a_j+\a_k) & \!\!\! + \!\!\! & \cos(\a_i-\a_j-\a_k) \\
!latex           4 \; \sin\a_i \sin\a_j \cos\a_k & \!\!\! = & \!\!\!
!latex       - \!\!\! & \cos(\a_i+\a_j+\a_k) & \!\!\! - \!\!\! & \cos(\a_i+\a_j-\a_k) & \!\!\! + \!\!\! & \cos(\a_i-\a_j+\a_k) & \!\!\! + \!\!\! & \cos(\a_i-\a_j-\a_k)
!latex       \end{array}
!latex       \ee

!latex \item The required derivatives are
!latex       \be \begin{array}{cclccccccccccccccccccccccccccccccccccccccccccccccc}
!latex           3 \ds \frac{\partial V}{\partial R_{e,i}} & = & \left( + Z_{e,j} R_{o,k} m_k - R_{e,j} Z_{o,k} m_k - R_{e,j} Z_{o.k} m_k \right) & 
!latex                                                           \ds \int \!\!\!\! \int \!\! d\t d\z \; \cos\a_i \cos\a_j \cos\a_k \\
!latex                                                     & + & \left( - Z_{o,j} R_{e,k} m_k + R_{o,j} Z_{e,k} m_k + R_{o,j} Z_{e,k} m_k \right) &
!latex                                                           \ds \int \!\!\!\! \int \!\! d\t d\z \; \cos\a_i \sin\a_j \sin\a_k \\
!latex                                                     & + & \left( - R_{o,k} Z_{e,j} m_i                                             \right) &
!latex                                                           \ds \int \!\!\!\! \int \!\! d\t d\z \; \sin\a_i \cos\a_j \sin\a_k \\
!latex                                                     & + & \left( - R_{e,k} Z_{o,j} m_i                                             \right) &
!latex                                                           \ds \int \!\!\!\! \int \!\! d\t d\z \; \sin\a_i \sin\a_j \cos\a_k
!latex       \end{array} \ee
!latex       \be \begin{array}{cclccccccccccccccccccccccccccccccccccccccccccccccc}
!latex           3 \ds \frac{\partial V}{\partial Z_{o,i}} & = & \left( - R_{e,k} R_{e,j} m_i                                             \right) & 
!latex                                                           \ds \int \!\!\!\! \int \!\! d\t d\z \; \cos\a_i \cos\a_j \cos\a_k \\
!latex                                                     & + & \left( - R_{o,k} R_{o,j} m_i                                             \right) &
!latex                                                           \ds \int \!\!\!\! \int \!\! d\t d\z \; \cos\a_i \sin\a_j \sin\a_k \\
!latex                                                     & + & \left( - R_{e,j} R_{e,k} m_k                                             \right) &
!latex                                                           \ds \int \!\!\!\! \int \!\! d\t d\z \; \sin\a_i \cos\a_j \sin\a_k \\
!latex                                                     & + & \left( + R_{o,j} R_{o,k} m_k                                             \right) &
!latex                                                           \ds \int \!\!\!\! \int \!\! d\t d\z \; \sin\a_i \sin\a_j \cos\a_k
!latex       \end{array} \ee
    
!latex \end{enumerate}

    do ii = 1, mn ; mi = im(ii) ; ni = in(ii) ; Rei = iRbc(ii,jvol) ; Roi = iRbs(ii,jvol) ; Zei = iZbc(ii,jvol) ; Zoi = iZbs(ii,jvol)
     do jj = 1, mn ; mj = im(jj) ; nj = in(jj) ; Rej = iRbc(jj,jvol) ; Roj = iRbs(jj,jvol) ; Zej = iZbc(jj,jvol) ; Zoj = iZbs(jj,jvol)
      do kk = 1, mn ; mk = im(kk) ; nk = in(kk) ; Rek = iRbc(kk,jvol) ; Rok = iRbs(kk,jvol) ; Zek = iZbc(kk,jvol) ; Zok = iZbs(kk,jvol)
       
       
       AA = Rei * ( Zej * Rok - Rej * Zok ) * (+mk)
       BB = Rei * ( Zoj * Rek - Roj * Zek ) * (-mk)
       CC = Roi * ( Zej * Rek - Rej * Zek ) * (-mk)
       DD = Roi * ( Zoj * Rok - Roj * Zok ) * (+mk)
       
       if( mi+mj+mk.eq.0 .and. ni+nj+nk.eq.0 ) vol(innout) = vol(innout) + AA - BB - CC - DD
       if( mi+mj-mk.eq.0 .and. ni+nj-nk.eq.0 ) vol(innout) = vol(innout) + AA + BB + CC - DD
       if( mi-mj+mk.eq.0 .and. ni-nj+nk.eq.0 ) vol(innout) = vol(innout) + AA + BB - CC + DD
       if( mi-mj-mk.eq.0 .and. ni-nj-nk.eq.0 ) vol(innout) = vol(innout) + AA - BB + CC + DD
       

       if( dBdX%L .and. dBdX%innout.eq.innout .and. dBdX%ii.eq.ii ) then ! compute derivative of volume;
        
        if( dBdX%irz.eq.0 ) then ! compute derivatives wrt R; 20 Jun 14;
         
         if( dBdX%issym.eq.0 ) then !     stellarator-symmetric harmonic; dV/dRei ; 13 Sep 13;
          ;                                          ; AA = + Zej * Rok * mk - Rej * Zok * mk - Rej * Zok * mk 
          ;                                          ; BB = - Zoj * Rek * mk + Roj * Zek * mk + Roj * Zek * mk
          ;                                          ; CC = - Rok * Zej * mi
          ;                                          ; DD = - Rek * Zoj * mi
         else                                        ! non-stellarator-symmetric harmonic; dV/dRoi ; 13 Sep 13;
          ;                                          ; AA = + Rek * Zej * mi                                   
          ;                                          ; BB = + Rok * Zoj * mi
          ;                                          ; CC = + Rej * Zek * mk - Zej * Rek * mk + Rej * Zek * mk
          ;                                          ; DD = + Zoj * Rok * mk - Roj * Zok * mk - Roj * Zok * mk
         endif                                       ! end of if( dBdX%issym.eq.0 ) ; 13 Sep 13;
         
        else ! matches if( dBdX%irz.eq.0 ) then; compute derivative wrt Z; 19 Sep 13;
         
         if( dBdX%issym.eq.0 ) then ! stellarator-symmetric harmonic; dV/dZoi ; 13 Sep 13;
          ;                                          ; AA = - Rek * Rej * mi 
          ;                                          ; BB = - Rok * Roj * mi
          ;                                          ; CC = - Rej * Rek * mk
          ;                                          ; DD = + Roj * Rok * mk
         else                                        !                               ; dV/dZei ; 13 Sep 13;
          ;                                          ; AA = + Rej * Rok * mk 
          ;                                          ; BB = - Roj * Rek * mk
          ;                                          ; CC = + Rok * Rej * mi
          ;                                          ; DD = + Rek * Roj * mi
         endif                                       ! end of if( dBdX%issym.eq.0 ) ; 13 Sep 13;
         
        endif ! end of if( dBdX%irz.eq.0 ) ; 13 Sep 13;
        
        if( mi+mj+mk.eq.0 .and. ni+nj+nk.eq.0 ) dvolume = dvolume + AA - BB - CC - DD
        if( mi+mj-mk.eq.0 .and. ni+nj-nk.eq.0 ) dvolume = dvolume + AA + BB + CC - DD
        if( mi-mj+mk.eq.0 .and. ni-nj+nk.eq.0 ) dvolume = dvolume + AA + BB - CC + DD
        if( mi-mj-mk.eq.0 .and. ni-nj-nk.eq.0 ) dvolume = dvolume + AA - BB + CC + DD
        
       endif ! end of if( dBdX%L .and. dBdX%ii.eq.ii .and. dBdX%innout.eq.innout ) then; 20 Jun 14;


      enddo ! end of do kk; 13 Sep 13;
     enddo ! end of do jj; 13 Sep 13;
    enddo ! end of do ii; 13 Sep 13;
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   end select
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
#ifdef DEBUG
   if( Wvolume ) then
    if( Igeometry.eq.3 .and. lvol.eq.1 .and. innout.eq.0 ) write(ounit,1000) myid, vol(innout) * pi2pi2nfpquart * third
    if( Igeometry.eq.4 .and. lvol.eq.1 .and. innout.eq.0 ) write(ounit,1000) myid, vol(innout) * pi2pi2nfpquart * third
   endif
1000 format("volume : ", 10x ," : myid="i3" ; axis contribution to volume = ",es13.5," ;")
#endif
   
  enddo ! end of innout loop; 26 Feb 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  if( Wvolume ) then
   cput = GETTIME
   write(ounit,'("volume : ",f10.2," : myid=",i3," ; lvol=",i3," ; vol=",2f15.10," ;")') cput-cpus, myid, lvol, vol(0:1)
  endif
#endif

  select case( Igeometry)
  case( 1 ) ; vvolume(lvol) = ( vol(1) - vol(0) ) * pi2pi2nfp              ; dvolume = dvolume * pi2pi2nfp                ! 20 Jun 14;
  case( 2 ) ; vvolume(lvol) = ( vol(1) - vol(0) ) * pi2pi2nfpquart         ; dvolume = dvolume * pi2pi2nfpquart
  case( 3 ) ; vvolume(lvol) = ( vol(1) - vol(0) ) * pi2pi2nfpquart * third ; dvolume = dvolume * pi2pi2nfpquart * third 
  case( 4 ) ; vvolume(lvol) = one                                          ; dvolume = zero ! this is under construction; 04 Dec 14;
   FATAL( volume, abs(pscale).gt.vsmall,need to compute volume )
  case default
   FATAL( volume, .true., invalid Igeometry )
  end select
  
  if( dBdX%innout.eq.0 ) dvolume = - dvolume
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( Wvolume ) then
   cput = GETTIME
   write(ounit,'("volume : ",f10.2," : myid=",i3," ; Igeometry=",i2," ; vvolume(",i3," ) =",es23.15" ;")') cput-cpus, myid, Igeometry, lvol, vvolume(lvol)
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  FATAL( volume, vflag.eq.0 .and. vvolume(lvol).lt.small, volume cannot be zero or negative ) ! 15 Jan 13;
  
  if( vvolume(lvol).lt.small ) then
   write(ounit,'("volume : ", 10x ," : myid=",i3," ; lvol=",i3," ; vvolume=",es13.5," ; volume cannot be zero or negative ;")') myid, lvol, vvolume(lvol)
   vvolume(lvol) = +9.9E+09
   vflag = 1
  else
   vflag = 0
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(volume)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine volume

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

