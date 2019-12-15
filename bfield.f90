!> \defgroup diagnostics Diagnostics to check the code.

!> \file bfield.f90
!! \brief Returns \f$\dot s \equiv B^s / B^\zeta\f$ and \f$\dot \theta \equiv B^\theta / B^\zeta\f$.
!! \ingroup diagnostics
!!
!! Returns the magnetic field field line equations, \f$ d{\bf x}/d\phi = {\bf B}/B^\phi \f$ .
!!
!! **Equations of field line flow**
!! <ul>
!! <li> The equations for the fieldlines are normalized to the toroidal field, i.e. 
!!       \f{equation}{ \dot s      \equiv \frac{B^s     }{B^\zeta}, \qquad
!!                     \dot \theta \equiv \frac{B^\theta}{B^\zeta}. \label{eq:stdot} \f} </li>
!! </ul>
!! **Representation of magnetic field**
!! <ul>
!! <li> The components of the vector potential, \f${\bf A}=A_\theta \nabla + A_\zeta \nabla \zeta\f$, are
!! \f{eqnarray}{
!! A_\theta(s,\theta,\zeta) &=& \sum_{i,l} {\color{red}{A_{\theta,e,i,l}}} \; {\overline T}_{l,i}(s) \cos\alpha_i + \sum_{i,l} {\color{ora}{A_{\theta,o,i,l}}} \; {\overline T}_{l,i}(s) \sin\alpha_i, \label{eq:At} \\
!! A_\zeta( s,\theta,\zeta) &=& \sum_{i,l} {\color{blu}{A_{\zeta, e,i,l}}} \; {\overline T}_{l,i}(s) \cos\alpha_i + \sum_{i,l} {\color{gra}{A_{\zeta ,o,i,l}}} \; {\overline T}_{l,i}(s) \sin\alpha_i, \label{eq:Az}
!! \f}
!! where \f${\overline T}_{l,i}(s) \equiv \bar s^{m_i/2} \, T_l(s)\f$, \f$T_l(s)\f$ is the Chebyshev polynomial, and \f$\alpha_j \equiv m_j\theta-n_j\zeta\f$.
!! The regularity factor, \f$\bar s^{m_i/2}\f$, where \f$\bar s \equiv (1+s)/2\f$, is only included if there is a coordinate singularity in the domain
!! (i.e. only in the innermost volume, and only in cylindrical and toroidal geometry.) </li>
!! <li> The magnetic field, \f$\sqrt g \, {\bf B} = \sqrt g B^s {\bf e}_s + \sqrt g B^\theta {\bf e}_\theta + \sqrt g B^\zeta {\bf e}_\zeta\f$, is
!! \f{eqnarray}{
!! \begin{array}{ccccrcrcrcrcccccccccccccccccccccccccccccccccccccccccccccccccccc}
!! \sqrt g \, {\bf B} & = & {\bf e}_s      & \sum_{i,l} [ ( & - m_i {\color{blu}{A_{\zeta, e,i,l}}} & - & n_i {\color{red}{A_{\theta,e,i,l}}} & ) {\overline T}_{l,i}        \sin\alpha_i + ( & + m_i {\color{gra}{A_{\zeta ,o,i,l}}} & + & n_i {\color{ora}{A_{\theta,o,i,l}}} & ) {\overline T}_{l,i}        \cos\alpha_i ] \\
!!                    & + & {\bf e}_\theta & \sum_{i,l} [ ( &                                       & - &     {\color{blu}{A_{\zeta, e,i,l}}} & ) {\overline T}_{l,i}^\prime \cos\alpha_i + ( &                                       & - &     {\color{gra}{A_{\zeta ,o,i,l}}} & ) {\overline T}_{l,i}^\prime \sin\alpha_i ] \\
!!                    & + & {\bf e}_\zeta  & \sum_{i,l} [ ( &       {\color{red}{A_{\theta,e,i,l}}} &   &                                     & ) {\overline T}_{l,i}^\prime \cos\alpha_i + ( &       {\color{ora}{A_{\theta,o,i,l}}} &   &                                     & ) {\overline T}_{l,i}^\prime \sin\alpha_i ]
!! \end{array}
!! \f}
!! </li>
!! <li> In Eqn.\f$\,\eqref{eq:stdot}\f$ , the coordinate Jacobian, \f$\sqrt g\f$, cancels.
!!   No coordinate metric information is required to construct the fieldline equations from the magnetic vector potential. </li>
!! </ul>
!! IT IS REQUIRED TO SET IVOL THROUGH GLOBAL MEMORY BEFORE CALLING BFIELD.
!!
!! The format of this subroutine is constrained by the NAG ode integration routines.
!!
!! @param zeta \f$ \zeta \f$
!! @param st I have no idea.
!! @param Bst Me neither.
subroutine bfield( zeta, st, Bst )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, one, half, two
  
  use numerical, only : vsmall, small
  
  use fileunits, only : ounit
  
  use inputlist, only : Wmacros, Wbfield, Lrad
  
  use cputiming, only : Tbfield
  
  use allglobal, only : myid, ncpu, cpus, mn, im, in, halfmm, regumm, &
                        ivol, gBzeta, Ate, Aze, Ato, Azo, &
                        NOTstellsym, cheby, &
                        Lcoordinatesingularity, Mvol, &
                        Node ! 17 Dec 15;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  REAL, intent(in)   :: zeta,  st(1:Node)
  REAL, intent(out)  ::       Bst(1:Node)
  
  INTEGER            :: lvol, ii, ll, mi, ni, ideriv
  REAL               :: teta, lss, sbar, sbarhm(0:1), arg, carg, sarg, dBu(1:3)
  
  REAL, allocatable  :: TT(:,:) ! this is almost identical to cheby; 17 Dec 15;

  BEGIN(bfield)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  FATAL( bfield, ivol.lt.1 .or. ivol.gt.Mvol, invalid ivol )
#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  lvol = ivol ; ideriv = 0 ! the argument list of bfield is fixed by NAG requirements, but volume index is required below;

  Bst(1:Node) = (/ zero , zero /) ! set default intent out; this should cause a compilation error if Node.ne.2;
 
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  lss = st(1) ; teta = st(2) 
  
  if( abs(lss).gt.one ) goto 9999 ! out of domain;
    
  if( Lcoordinatesingularity ) sbar = max( ( lss + one ) * half, small )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  cheby(0,0:1) = (/ one, zero /) ! Chebyshev initialization; 16 Jan 13;
  cheby(1,0:1) = (/ lss,  one /)
  do ll = 2, Lrad(lvol) ; cheby(ll,0:1) = (/ two * lss * cheby(ll-1,0)                             - cheby(ll-2,0) , & ! Chebyshev recurrence; 17 Dec 15;
                                             two       * cheby(ll-1,0) + two * lss * cheby(ll-1,1) - cheby(ll-2,1)   /)
  enddo
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  SALLOCATE( TT, (0:Lrad(lvol),0:1), zero ) ! unless a regularization factor is required, TT = chebyshev polynomial; 17 Dec 15;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  dBu(1:3) = zero ! initialize summation;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  do ii = 1, mn ; mi = im(ii) ; ni = in(ii) ; arg = mi * teta - ni * zeta ; carg = cos(arg) ; sarg = sin(arg) ! shorthand; 20 Apr 13;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   if( Lcoordinatesingularity ) then ! regularization factor depends on mi; 17 Dec 15;

    FATAL( bfield, abs(sbar).lt.vsmall, need to avoid divide-by-zero )

    if( mi.eq.0 ) then ; sbarhm(0) = one              ; sbarhm(1) = zero
    else               ; sbarhm(0) = sbar**regumm(ii) ; sbarhm(1) = half * regumm(ii) * sbarhm(0) / sbar 
    endif   

    do ll = 0, Lrad(lvol) ; TT(ll,0:1) = (/ sbarhm(0) * cheby(ll,0), sbarhm(0) * cheby(ll,1) + sbarhm(1) * cheby(ll,0) /)
    enddo

   else

    do ll = 0, Lrad(lvol) ; TT(ll,0:1) = (/             cheby(ll,0),             cheby(ll,1)                           /)
    enddo

   endif ! end of if( Lcoordinatesingularity ) ; 16 Jan 15;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   do ll = 0, Lrad(lvol) ! loop over Chebyshev summation; 20 Feb 13;
    ;dBu(1) = dBu(1) + ( - mi * Aze(lvol,ideriv,ii)%s(ll) - ni * Ate(lvol,ideriv,ii)%s(ll) ) * TT(ll,0) * sarg
    ;dBu(2) = dBu(2) + (                                  -      Aze(lvol,ideriv,ii)%s(ll) ) * TT(ll,1) * carg
    ;dBu(3) = dBu(3) + (        Ate(lvol,ideriv,ii)%s(ll)                                  ) * TT(ll,1) * carg
    if( NOTstellsym ) then ! include non-symmetric harmonics; 28 Jan 13;
     dBu(1) = dBu(1) + ( + mi * Azo(lvol,ideriv,ii)%s(ll) + ni * Ato(lvol,ideriv,ii)%s(ll) ) * TT(ll,0) * carg
     dBu(2) = dBu(2) + (                                  -      Azo(lvol,ideriv,ii)%s(ll) ) * TT(ll,1) * sarg
     dBu(3) = dBu(3) + (        Ato(lvol,ideriv,ii)%s(ll)                                  ) * TT(ll,1) * sarg
    endif
   enddo ! end of do ll; 10 Dec 15;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
  enddo ! end of do ii = 1, mn;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  gBzeta = dBu(3) ! gBzeta is returned through global; 20 Apr 13;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( abs(gBzeta).lt.vsmall ) then

   cput = GETTIME

   write(ounit,'("bfield : ",f10.2," : lvol=",i3," ; zeta="es23.15" ; (s,t)=("es23.15" ,"es23.15" ) ; B^z="es23.15" ;")') &
                             cput-cpus, lvol,        zeta,             st(1:2),                       dBu(3)

   FATAL( bfield, abs(dBu(3)).lt.vsmall, field is not toroidal )

  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  Bst(1:2) = dBu(1:2) / gBzeta ! normalize field line equations to toroidal field; 20 Apr 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  DALLOCATE(TT)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(bfield)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine bfield

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
