!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (initialization) ! Allocates and initializes internal arrays.

!latex \briefly{Allocates and initializes internal arrays.}

!latex \calledby{\link{xspech}}
!latex \calls{\link{ra00aa}}

!latex \tableofcontents

!latex \subsection{definition of internal variables}

subroutine preset
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, one
  
  use numerical, only : sqrtmachprec, vsmall, small
  
  use fileunits, only : ounit
  
  use inputlist! only :
  
  use cputiming, only : Tpreset
  
  use allglobal! only :
  
  use fftw_interface

  use zernik

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS
  
  INTEGER   :: innout, idof, jk, ll, ii, ifail, ideriv, vvol, mi, ni, mj, nj, mk, nk, mimj, ninj, mkmj, nknj, jj, kk, lvol, mm, nn, imn
  INTEGER   :: lquad, igauleg, maxIquad, Mrad, jquad, Lcurvature, zerdof
  REAL      :: teta, zeta, arg, lss, cszeta(0:1), error
  
  BEGIN(preset)
  
  call random_seed()
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
 !FATAL( preset, Nfp.eq.0, illegal division ) ! this was checked in global: readin; SRH: 27 Feb 18;

  pi2nfp         = pi2 / Nfp
  
  pi2pi2nfp      = pi2 * pi2nfp
  pi2pi2nfphalf  = pi2 * pi2nfp * half
  pi2pi2nfpquart = pi2 * pi2nfp * quart

  Mrad  = maxval( Lrad(1:Mvol) )

  if( myid.eq.0 ) write(ounit,'("preset : ",10x," : myid=",i3," ; Mrad=",i3," : Lrad=",257(i3,",",:))') myid, Mrad, Lrad(1:Mvol)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \subsubsection{\type{LGdof} and \type{NGdof} : number of geometrical degrees-of-freedom;}

!latex \begin{enumerate}
!latex \item \type{LGdof} $\equiv$ the number of degrees-of-freedom in the geometry (i.e. Fourier harmonics) of each interface;
!latex \item \type{NGdof} $\equiv$ total number of degrees-of-freedom in geometry, i.e. of all interfaces;
!latex \end{enumerate}

!                            Rbc  Zbs    Rbs    Zbc
  select case( Igeometry )
  case( 1:2)
   if( YESstellsym ) LGdof = mn
   if( NOTstellsym ) LGdof = mn        + mn-1
  case(   3)
   if( YESstellsym ) LGdof = mn + mn-1
   if( NOTstellsym ) LGdof = mn + mn-1 + mn-1 + mn
  end select
  
  NGdof = ( Mvol-1 ) * LGdof
  
  if( Wpreset ) then ; cput = GETTIME ; write(ounit,'("preset : ",f10.2," : myid=",i3," ; NGdof=",i9," ;")') cput-cpus, myid, NGdof
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{\inputvar{iota} and \inputvar{oita} : rotational transform on interfaces;}

!latex \begin{enumerate}
!latex \item The input variables \inputvar{iota} and \inputvar{oita} are the rotational transform 
!latex       on ``inner-side'' and on the ``outer-side'' of each interface;
!latex \item These quantities are formally input
!latex \item Note that if $q_l+\gamma q_r \ne 0$, then \inputvar{iota} is given by
!latex       \be \iotabar \equiv \frac{p_l + \gamma p_r}{q_l + \gamma q_r},
!latex       \ee
!latex       where $p_l \equiv $\inputvar{pl}, $q_l \equiv $\inputvar{ql}, etc.;
!latex       and similarly for \inputvar{oita}.
!latex \end{enumerate}
  
  do vvol = 0, Nvol
   
   if( ql(vvol).eq.0 .and. qr(vvol).eq.0 ) then ; iota(vvol) = iota(vvol)
   else                                         ; iota(vvol) = ( pl(vvol) + goldenmean * pr(vvol) ) / ( ql(vvol) + goldenmean * qr(vvol) )
   endif
   
   if( lq(vvol).eq.0 .and. rq(vvol).eq.0 ) then ; oita(vvol) = oita(vvol)
   else                                         ; oita(vvol) = ( lp(vvol) + goldenmean * rp(vvol) ) / ( lq(vvol) + goldenmean * rq(vvol) )
   endif
   
   if( Wpreset .and. myid.eq.0 ) write(ounit,1002) vvol, pl(vvol), ql(vvol), pr(vvol), qr(vvol), iota(vvol), lp(vvol), lq(vvol), rp(vvol), rq(vvol), oita(vvol)
   
1002 format("preset : ",10x," :      ",3x," ; transform : ",i3," : (",i3," /",i3," ) * (",i3," /",i3," ) = ",f18.15," ; ",&
                                                                  "(",i3," /",i3," ) * (",i3," /",i3," ) = ",f18.15," ; ")
   
  enddo
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \subsubsection{\type{dtflux(1:Mvol)} and \type{dpflux(1:Mvol)} : enclosed fluxes;}

!latex \begin{enumerate}
!latex \item \type{dtflux} $\equiv \Delta \psi_{tor} / 2\pi$ and 
!latex       \type{dpflux} $\equiv \Delta \psi_{pol} / 2\pi$ in each volume.
!latex \item (Note that the total toroidal flux enclosed by the plasma boundary is $\Phi_{edge} \equiv$ \inputvar{phiedge}.)
!latex \item $\psi_{tor} \equiv$ \inputvar{tflux} and $\psi_{pol} \equiv$ \inputvar{pflux} are immediately normalized (in \link{global}) according to
!latex       $\psi_{tor,i} \rightarrow \psi_{tor,i} / \psi_{0}$ and
!latex       $\psi_{pol,i} \rightarrow \psi_{pol,i} / \psi_{0}$, where $\psi_{0} \equiv \psi_{tor,N}$ on input.
!latex \end{enumerate}

  SALLOCATE( dtflux, (1:Mvol), zero )
  SALLOCATE( dpflux, (1:Mvol), zero )

  select case( Igeometry )
  case( 1   ) ; dtflux(1) = tflux(1) ; dpflux(1) = pflux(1) ! Cartesian              ; this is the "inverse" operation defined in xspech; 09 Mar 17;
  case( 2:3 ) ; dtflux(1) = tflux(1) ; dpflux(1) = zero     ! cylindrical or toroidal; 
  end select

  dtflux(2:Mvol) = tflux(2:Mvol) - tflux(1:Mvol-1)
  dpflux(2:Mvol) = pflux(2:Mvol) - pflux(1:Mvol-1)
  
  dtflux(1:Mvol) = dtflux(1:Mvol) * phiedge / pi2
  dpflux(1:Mvol) = dpflux(1:Mvol) * phiedge / pi2
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \subsubsection{\type{sweight(1:Mvol)} : star-like angle constraint weight;}

!latex \begin{enumerate}
!latex \item the ``star-like'' poloidal angle constraint weights (only required for toroidal geometry, i.e. \type{Igeometry=3}) are given by
!latex       \be \type{sweight}_v \equiv \inputvar{upsilon} \times (l_v / N_{vol})^w,
!latex       \ee
!latex       where $l_v$ is the volume number,
!latex       and $w \equiv $ \inputvar{wpoloidal}.
!latex \end{enumerate}
  
  SALLOCATE( sweight, (1:Mvol), zero )
 !sweight(1:Mvol) = upsilon * tflux(1:Mvol)**wpoloidal ! toroidal flux in vacuum region is not constant; 11 July 18;
  do vvol = 1, Mvol ; sweight(vvol) = upsilon * (vvol*one/Nvol)**wpoloidal ! 11 July 18;
  enddo
  
#ifdef DEBUG
  write(ounit,'("preset : ",10x," : sweight=",99(es12.5,",",:))') sweight(1:Mvol)
#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \subsubsection{\type{TT(0:Mrad,0:1,0:1)} : Chebyshev polynomials at inner/outer interface;}

!latex \begin{enumerate}
!latex \item \type{TT(0:Lrad,0:1,0:1)} gives the Chebyshev polynomials, and their first derivative, evaluated at $s=-1$ and $s=+1$.
!latex \item Precisely, \type{TT(l,i,d)} $\equiv T_l^{(d)}(s_i)$ for $s_0=-1$ and $s_1=+1$.
!latex \item Note that $T_l^{(0)}(s)=s^l$ and $T_l^{(1)}(s)=s^{l+1} l^2$ for $s=\pm 1$.
!latex \item Note that 
!latex       \be T_l(-1)        = \left\{ \begin{array}{ccccccccccccccc}+1,& \mbox{\rm if $l$ is even,} \\ 
!latex                                                                  -1,& \mbox{\rm if $l$ is odd;}
!latex                                    \end{array} \right. & \; \;&
!latex           T_l(+1)        = \left\{ \begin{array}{ccccccccccccccc}+1,& \mbox{\rm if $l$ is even,} \\ 
!latex                                                                  +1,& \mbox{\rm if $l$ is odd;}
!latex                                    \end{array} \right. \\
!latex           T_l^\prime(-1) = \left\{ \begin{array}{ccccccccccccccc}-l^2,& \mbox{\rm if $l$ is even,} \\ 
!latex                                                                  +l^2,& \mbox{\rm if $l$ is odd;}
!latex                                    \end{array} \right. &\; \;&
!latex           T_l^\prime(+1) = \left\{ \begin{array}{ccccccccccccccc}+l^2,& \mbox{\rm if $l$ is even,} \\ 
!latex                                                                  +l^2,& \mbox{\rm if $l$ is odd.}
!latex                                     \end{array} \right.
!latex       \ee
!latex \item \type{TT(0:Mrad,0:1,0:1)} is used in routines that explicity require interface information, such as
!latex       \begin{enumerate}
!latex         \item the interface force-balance routine,                                \link{lforce}; 
!latex         \item the virtual casing routine,                                         \link{casing}; 
!latex         \item computing the rotational-transform on the interfaces,               \link{tr00ab}; 
!latex         \item computing the covariant components of the interface magnetic field, \link{sc00aa};
!latex         \item enforcing the constraints on the Beltrami fields,                   \link{matrix};
!latex     and \item computing the enclosed currents of the vacuum field,                \link{curent}.
!latex       \end{enumerate}
!latex \end{enumerate}

  SALLOCATE( TT, (0:Mrad,0:1,0:1), zero )
  SALLOCATE(RTT, (0:Lrad(1),0:Mpol,0:1,0:1), zero )
  SALLOCATE(RTM, (0:Lrad(1),0:Mpol), zero )
  
  do innout = 0, 1 ; lss = two * innout - one
   
   do ll = 0, Mrad ; TT(ll,innout,0) = lss**(ll  )        
    ;              ; TT(ll,innout,1) = lss**(ll+1) * ll**2 ! derivative; 26 Jan 16;
   enddo
   
  enddo ! end of do innout = 0, 1 ;

  call get_zernike( zero, Lrad(1), Mpol, RTT(:,:,0,:))
  call get_zernike( one, Lrad(1), Mpol, RTT(:,:,1,:))
  call get_zernike_rm(zero, Lrad(1), Mpol, RTM(:,:))

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \subsubsection{\type{ImagneticOK(1:Mvol)} : Beltrami/vacuum error flag;}
  
!latex \begin{enumerate}
!latex \item error flags that indicate if the magnetic field in each volume has been successfully constructed;
!latex \item \type{ImagneticOK} is initialized to \type{.false.} in \link{dforce} before the Beltrami solver routines are called.
!latex       If the construction of the Beltrami field is successful
!latex       (in either \link{ma02aa} or \link{mp00ac})
!latex       then \type{ImagneticOK} is set to \type{.true.}.
!latex \end{enumerate}

  SALLOCATE( ImagneticOK, (1:Mvol), .false. )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{\type{Lhessianallocated}}

!latex \begin{enumerate}
!latex \item The internal logical variable, \type{Lhessianallocated}, indicates whether the ``Hessian'' matrix of second-partial derivatives
!latex       (really, the first derivatives of the force-vector) has been allocated, or not!
!latex \end{enumerate}
  
  Lhessianallocated = .false.
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \subsubsection{\type{ki(1:mn,0:1)} : Fourier identification;}

!latex \begin{enumerate}
!latex \item Consider the ``abbreviated'' representation for a double Fourier series,
!latex       \be \sum_i f_i \cos(m_i \t - n_i \z) \equiv                         \sum_{n=      0  }^{     N_0} f_{0,n} \cos(    -n\z)
!latex                                                   + \sum_{m=1}^{     M_0} \sum_{n=-     N_0}^{     N_0} f_{m,n} \cos( m\t-n\z),
!latex       \ee
!latex       and the same representation but with enhanced resolution,
!latex       \be \sum_k \bar f_k \cos(\bar m_k \t - \bar n_k \z) \equiv                         \sum_{n=      0  }^{     N_1} f_{0,n} \cos(    -n\z)
!latex                                                   + \sum_{m=1}^{     M_1} \sum_{n=-     N_1}^{     N_1} f_{m,n} \cos(m\t-n\z),
!latex       \label{eq:enhancedFourierrepresentation}
!latex       \ee
!latex       with $M_1 \ge M_0$ and $N_1 \ge N_0$; 
!latex       \newline then $k_i\equiv$ \type{ki(i,0)} is defined such that $\bar m_{k_i} = m_i$ and $\bar n_{k_i} = n_i$.
!latex \end{enumerate}

  SALLOCATE( ki, (1:mn,0:1), 0 )
  
  do ii = 1, mn  ; mi = im(ii)  ; ni = in(ii)
   
   do kk = 1, mne ; mk = ime(kk) ; nk = ine(kk)
    
    if( mk.eq. mi .and. nk.eq. ni ) then
     if( mk.eq.0 .and. nk.eq.0 ) then ; ki(ii,0:1) = (/ kk, 1 /)
     else                             ; ki(ii,0:1) = (/ kk, 2 /)
     endif
    endif
    
   enddo ! end of do kk = 1, mne ;
   
  enddo ! end of do ii ;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{\type{kija(1:mn,1:mn,0:1)}, \type{kijs(1:mn,1:mn,0:1)} : Fourier identification;}

!latex \begin{enumerate}
!latex \item Consider the following quantities, which are computed in \link{ma00aa}, 
!latex       where $\bar g^{\mu\nu} = \sum_k \bar g^{\mu\nu}_k \cos \a_k$ for $\a_k \equiv m_k \t - n_k \z$, 
!latex       \be
!latex       \ooint \bar g^{\mu\nu} \cos\a_i \; \cos\a_j & = & \frac{1}{2} \ooint \bar g^{\mu\nu} ( + \cos\a_{k_{ij+}} + \cos\a_{k_{ij-}} ), \\
!latex       \ooint \bar g^{\mu\nu} \cos\a_i \; \sin\a_j & = & \frac{1}{2} \ooint \bar g^{\mu\nu} ( + \sin\a_{k_{ij+}} - \sin\a_{k_{ij-}} ), \\
!latex       \ooint \bar g^{\mu\nu} \sin\a_i \; \cos\a_j & = & \frac{1}{2} \ooint \bar g^{\mu\nu} ( + \sin\a_{k_{ij+}} + \sin\a_{k_{ij-}} ), \\
!latex       \ooint \bar g^{\mu\nu} \sin\a_i \; \sin\a_j & = & \frac{1}{2} \ooint \bar g^{\mu\nu} ( - \cos\a_{k_{ij+}} + \cos\a_{k_{ij-}} ),
!latex       \ee
!latex       where $(m_{k_{ij+}}, n_{k_{ij+}}) = (m_i + m_j, n_i + n_j)$ and $(m_{k_{ij-}}, n_{k_{ij-}}) = (m_i - m_j, n_i - n_j)$;
!latex       \newline then \mbox{\type{kija(i,j,0)}$\equiv k_{ij+}$} and \mbox{\type{kijs(i,j,0)}$\equiv k_{ij-}$}.
!latex \item Note that \Eqn{enhancedFourierrepresentation} does not include $m<0$; so,
!latex       if $m_i - m_j < 0$ then $k_{ij-}$ is re-defined such that \mbox{$(m_{k_{ij-}}, n_{k_{ij-}})$} $=$ \mbox{$ (m_j - m_i, n_j - n_i)$}; and
!latex       similarly for the case $m=0$ and $n<0$.
!latex       Also, take care that the sign of the sine harmonics in the above expressions will change for these cases.
!latex \end{enumerate}

  SALLOCATE( kija, (1:mn,1:mn,0:1), 0 )
  
  do ii = 1, mn  ; mi =  im(ii) ; ni =  in(ii)
   
   do jj = 1, mn  ; mj =  im(jj) ; nj =  in(jj) ; mimj = mi + mj ; ninj = ni + nj !   adding   ; 17 Dec 15;
    
    do kk = 1, mne ; mk = ime(kk) ; nk = ine(kk)
     
     if( mk.eq. mimj .and. nk.eq. ninj ) then
      if( mk.eq.0 .and. nk.eq.0 ) then ; kija(ii,jj,0:1) = (/ kk ,   1 /)
      else                             ; kija(ii,jj,0:1) = (/ kk ,   2 /)
      endif
     endif
     
    enddo ! end of do kk; 29 Jan 13;
    
   enddo ! end of do jj; 29 Jan 13;
   
  enddo ! end of do ii; 29 Jan 13;


  SALLOCATE( kijs, (1:mn,1:mn,0:1), 0 )
  
  do ii = 1, mn  ; mi =  im(ii) ; ni =  in(ii)
   
   do jj = 1, mn  ; mj =  im(jj) ; nj =  in(jj) ; mimj = mi - mj ; ninj = ni - nj ! subtracting; 17 Dec 15;
    
    do kk = 1, mne ; mk = ime(kk) ; nk = ine(kk)
     
     if( mimj.gt.0 .or. ( mimj.eq.0 .and. ninj.ge.0 ) ) then ! no re-definition required; 17 Dec 15;
      
      if( mk.eq. mimj .and. nk.eq. ninj ) then
       if( mk.eq.0 .and. nk.eq.0 ) then ; kijs(ii,jj,0:1) = (/ kk ,   1 /)
       else                             ; kijs(ii,jj,0:1) = (/ kk ,   2 /)
       endif
      endif
      
     else
      
      if( mk.eq.-mimj .and. nk.eq.-ninj ) then
       ;                                ; kijs(ii,jj,0:1) = (/ kk , - 2 /) ! only the sine modes need the sign factor; 17 Dec 15;
      endif
      
     endif
     
    enddo ! end of do kk; 29 Jan 13;
    
   enddo ! end of do jj; 29 Jan 13;
   
  enddo ! end of do ii; 29 Jan 13;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{\type{djkp}}

  if( Igeometry.eq.2 ) then ! standard cylindrical; 04 Dec 14;
   
   SALLOCATE( djkp, (1:mn,1:mn), 0 ) ! only used in volume; trigonometric identities; 04 Dec 14;
   SALLOCATE( djkm, (1:mn,1:mn), 0 ) ! only used in volume; trigonometric identities; 04 Dec 14;
   
   do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
    do jj = 1, mn ; mj = im(jj) ; nj = in(jj)
     if( mi-mj.eq.0 .and. ni-nj.eq.0 ) djkp(ii,jj) = 1
     if( mi+mj.eq.0 .and. ni+nj.eq.0 ) djkm(ii,jj) = 1
    enddo
   enddo
   
  endif ! end of if( Igeometry.eq.2 ) ; 04 Dec 14;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \subsubsection{\type{iotakki}}

  SALLOCATE( iotakkii, (1:mn      ), 0 ) ! used to identify matrix elements in straight-field-line angle transformation;
  
  SALLOCATE( iotaksub, (1:mn,1:mns), 0 )
  SALLOCATE( iotaksgn, (1:mn,1:mns), 0 )
  SALLOCATE( iotakadd, (1:mn,1:mns), 0 )
  
  do kk = 1, mn ; mk = im(kk) ; nk = in(kk)
   
   
   do ii = 1, mns ; mi = ims(ii) ; ni = ins(ii)
    
    if( mk.eq.mi .and. nk.eq.ni ) iotakkii(kk) = ii
    
   enddo
   
   
   do jj = 1, mns ; mj = ims(jj) ; nj = ins(jj)
    
    
    mkmj = mk - mj ; nknj = nk - nj
    
    do ii = 1, mns ; mi = ims(ii) ; ni = ins(ii)
     
     if( mkmj.gt.0 .or. ( mkmj.eq.0 .and. nknj.ge.0 ) ) then
      
      if( mi.eq. mkmj .and. ni.eq. nknj ) then ; iotaksub(kk,jj) = ii ; iotaksgn(kk,jj) =  1
      endif
      
     else
      
      if( mi.eq.-mkmj .and. ni.eq.-nknj ) then ; iotaksub(kk,jj) = ii ; iotaksgn(kk,jj) = -1
      endif
      
     endif
     
    enddo ! end of do ii; 30 Jan 13;
    
    
    mkmj = mk + mj ; nknj = nk + nj
    
    do ii = 1, mns ; mi = ims(ii) ; ni = ins(ii)
     
     if( mi.eq. mkmj .and. ni.eq. nknj ) then ; iotakadd(kk,jj) = ii
     endif
     
    enddo ! end of do ii; 29 Jan 13;
    
    
   enddo ! end of do jj; 29 Jan 13;
   
   
  enddo ! end of do kk; 29 Jan 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{\type{cheby(0:Lrad,0:2)} : Chebyshev polynomial workspace;}

!latex \begin{enumerate}
!latex \item \type{cheby(0:Lrad,0:2)} is global workspace for computing the Chebyshev polynomials, and their derivatives,
!latex       using the recurrence relations $T_0(s) = 1$, $T_1(s) = s$ and  $T_l(s) = 2 \, s \,T_{l-1}(s) - T_{l-2}(s)$.
!latex \item These are computed as required, i.e. for arbitrary $s$, in \link{bfield}, \link{jo00aa} and \link{ma00aa}.
!latex \item (Note that the quantities required for \link{ma00aa} are for fixed $s$, and so these quantities should be precomputed.)
!latex \end{enumerate}

  SALLOCATE( cheby, (0:Mrad,0:2), zero )
  SALLOCATE( zernike, (0:Lrad(1), 0:Mpol, 0:2), zero )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{\type{Iquad}, \type{gaussianweight}, \type{gaussianabscissae} : Gaussian quadrature;}

!latex \begin{enumerate}
!latex \item The volume integrals are computed using a ``Fourier'' integration over the angles and by Gaussian quadrature over the radial,
!latex       i.e. $\ds \int \!\! f(s) ds = \sum_k \omega_k f(s_k)$.
!l tex \item The numerical resolution of the Gaussian quadrature is determined primarily by \inputvar{Nquad}, 
!l tex       but \inputvar{Lrad$_v$} is also important 
!l tex       (as is \inputvar{Mpol} if the regularization factors are included in the vector potential -- only for the coordinate singularity).
!latex \item The quadrature resolution in each volume is give by \internal{Iquad(1:Mvol)} which is determined as follows:
!latex \bi
!latex \item if \inputvar{Nquad.gt.0},                                 then \internal{Iquad(vvol) =              Nquad};
!latex \item if \inputvar{Nquad.le.0 and .not.Lcoordinatesingularity}, then \internal{Iquad(vvol) = 2*Lrad(vvol)-Nquad};
!latex \item if \inputvar{Nquad.le.0 and      Lcoordinatesingularity}, then \internal{Iquad(vvol) = 2*Lrad(vvol)-Nquad+Mpol};
!latex \ei
!latex \item The Gaussian weights and abscissae are given by \internal{gaussianweight(1:maxIquad,1:Mvol)} and \internal{gaussianabscissae(1:maxIquad,1:Mvol)},
!latex       which are computed using modified Numerical Recipes routine gauleg.
!latex \item \internal{Iquad$_v$} is passed through to \link{ma00aa} to compute the volume integrals of the metric elements;
!latex       also see \link{jo00aa}, where \internal{Iquad$_v$} is used to compute the volume integrals of $||\nabla\times{\bf B} - \mu {\bf B}||$;
!latex \end{enumerate}

  SALLOCATE( Iquad, (1:Mvol), 0 ) ! 16 Jan 13;
  
  do vvol = 1, Mvol

   LREGION(vvol)
   
   if( Nquad.gt.0 ) then ;            Iquad(vvol) =                         Nquad
   else                 
    if(      Lcoordinatesingularity ) Iquad(vvol) = Mpol + 2 * Lrad(vvol) - Nquad ! NEED TO REVISE REGULARIZATION FACTORS; 26 Feb 13;
    if( .not.Lcoordinatesingularity ) Iquad(vvol) =        2 * Lrad(vvol) - Nquad
   endif
   
  enddo ! end of do vvol; 18 Feb 13;
  
  maxIquad = maxval(Iquad(1:Mvol))
  
  SALLOCATE( gaussianweight   , (1:maxIquad,1:Mvol), zero ) ! perhaps it would be neater to make this a structure; 26 Jan 16;
  SALLOCATE( gaussianabscissae, (1:maxIquad,1:Mvol), zero )

  do vvol = 1, Mvol
   
   lquad = Iquad(vvol)
   
   call gauleg( lquad, gaussianweight(1:lquad,vvol), gaussianabscissae(1:lquad,vvol), igauleg ) ! JAB; 28 Jul 17
   
   if( myid.eq.0 ) then
    cput= GETTIME
    select case( igauleg ) !                                                  123456789012345
    case( 0 )    ; if( Wpreset ) write(ounit,1000) cput-cpus, vvol, igauleg, "success        ", gaussianabscissae(1:lquad,vvol)
    case( 1 )    ;               write(ounit,1000) cput-cpus, vvol, igauleg, "failed         ", gaussianabscissae(1:lquad,vvol)
    case( 2 )    ;               write(ounit,1000) cput-cpus, vvol, igauleg, "input error    ", gaussianabscissae(1:lquad,vvol)
    case default ;               write(ounit,1000) cput-cpus, vvol, igauleg, "weird          ", gaussianabscissae(1:lquad,vvol)
     FATAL( preset, .true., weird ifail returned by gauleg )
    end select
    ;            ; if( Wpreset ) write(ounit,1001)                                              gaussianweight(1:lquad,vvol)
   endif
   
1000 format("preset : ",f10.2," : lvol=",i3," ; igauleg=",i5," ; ",a15," ; abscissae ="99f09.05)
1001 format("preset : ", 10x ," :      ",3x,"           ",5x,"   ",15x," ; weights   ="99f09.05)
   
  enddo ! end of do vvol;  7 Mar 13; 
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
 
!latex \subsubsection{\type{LBsequad}, \type{LBnewton} and \type{LBlinear}}
  
!latex \begin{enumerate}
!latex \item \type{LBsequad}, \type{LBnewton} and \type{LBlinear} depend simply on \type{LBeltrami}, which is described in \link{global}.
!latex \end{enumerate}
  
  LBsequad = .false.
  LBnewton = .false.
  LBlinear = .false.
  
  if( LBeltrami.eq.1 .or. LBeltrami.eq.3 .or. LBeltrami.eq.5 .or. LBeltrami.eq.7 ) LBsequad = .true.
  if( LBeltrami.eq.2 .or. LBeltrami.eq.3 .or. LBeltrami.eq.6 .or. LBeltrami.eq.7 ) LBnewton = .true.
  if( LBeltrami.eq.4 .or. LBeltrami.eq.5 .or. LBeltrami.eq.6 .or. LBeltrami.eq.7 ) LBlinear = .true.
  
  if (LBnewton .or. LBsequad) Lconstraint = 2

  if (Lconstraint .eq. 2) then
    FATAL( preset, Lfreebound.eq.1, The combination of helicity constraint and free boundary is under construction )
    if (Igeometry .eq. 3) then
      write(ounit, *) 'WARNING: The Hessian matrix needs further review for Igeometry = 3'
      write(ounit, *) '         However, it can still serve the purpose of Lfindzero = 2'
    endif
  endif

  if( myid.eq.0 ) then
   cput = GETTIME
   write(ounit,'("preset : ",f10.2," : LBsequad="L2" , LBnewton="L2" , LBlinear="L2" ;")')cput-cpus, LBsequad, LBnewton, LBlinear
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \subsubsection{\type{BBweight(1:mn)} : weighting of force-imbalance harmonics}

!latex \begin{enumerate}
!latex \item weight on force-imbalance harmonics;
!latex       \be \type{BBweight}_i \equiv \inputvar{opsilon} \times \exp\left[ - \inputvar{escale} \times (m_i^2 + n_i^2) \right]
!latex       \ee
!latex \item this is only used in \link{dforce} in constructing the force-imbalance vector;
!latex \end{enumerate}

  SALLOCATE( BBweight, (1:mn), opsilon * exp( - escale * ( im(1:mn)**2 + (in(1:mn)/Nfp)**2 ) ) )
  
  if( myid.eq.0 .and. escale.gt.small ) then
   do ii = 1, mn ; write(ounit,'("preset : " 10x " : myid="i3" ; ("i3","i3") : BBweight="es13.5" ;")') myid, im(ii), in(ii)/Nfp, BBweight(ii)
   enddo
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{\type{mmpp(1:mn)} : spectral condensation weight factors}

!latex \begin{enumerate}
!latex \item spectral condensation weight factors;
!latex       \be \type{mmpp(i)} \equiv m_i^p,
!latex       \ee
!latex       where $p \equiv $ \inputvar{pcondense}.
!latex \end{enumerate}
  
  SALLOCATE( mmpp, (1:mn), zero )
  
  do ii = 1, mn ; mi = im(ii)
   
   if( mi.eq.0 ) then ; mmpp(ii) = zero
   else               ; mmpp(ii) = mi**pcondense
   endif ! end of if( mi.eq.0 ) ; 11 Aug 14;
   
  enddo ! end of do ii; 08 Nov 13;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \subsubsection{\type{NAdof}, \type{Ate}, \type{Aze}, \type{Ato} and \type{Azo} : degrees-of-freedom in magnetic vector potential}
  
!latex \begin{enumerate}
!latex \item \type{NAdof(1:Mvol)} $\equiv$ total number of degrees-of-freedom in magnetic vector potential, including Lagrange multipliers, in each volume.
!latex       This can de deduced from \link{matrix}.   
!latex \item \Ais
!latex \item The Chebyshev-Fourier harmonics of the covariant components of the magnetic vector potential are kept in
!latex       \\ $\Ate{v,j,l} \equiv $ \verb!Ate(v,0,j)%s(l)!, 
!latex       \\ $\Aze{v,j,l} \equiv $ \verb!Aze(v,0,j)%s(l)!,
!latex       \\ $\Ato{v,j,l} \equiv $ \verb!Ato(v,0,j)%s(l)!, and 
!latex       \\ $\Azo{v,j,l} \equiv $ \verb!Azo(v,0,j)%s(l)!; 
!latex       \\ where $v=1,\type{Mvol}$ labels volume, $j=1,\type{mn}$ labels Fourier harmonic, and $l=0,$ \inputvar{Lrad}$(v)$ labels Chebyshev polynomial.
!latex       (These arrays also contains derivative information.)
!latex \item If \inputvar{Linitguess=1}, a guess for the initial state for the Beltrami fields is constructed. 
!latex       An initial state is required for iterative solvers of the Beltrami fields, see \inputvar{LBeltrami}.
!latex \item If \inputvar{Linitguess=2}, the initial state for the Beltrami fields is read from file (see \link{ra00aa}).
!latex       An initial state is required for iterative solvers of the Beltrami fields, see \inputvar{LBeltrami}.
!latex \end{enumerate}
  
  SALLOCATE( NAdof, (1:Mvol          ), 0 ) ! Beltrami degrees-of-freedom in each annulus;
  
  NALLOCATE( Ate  , (1:Mvol,-1:2,1:mn)    ) ! recall that this is type:sub-grid; 31 Jan 13;
  NALLOCATE( Aze  , (1:Mvol,-1:2,1:mn)    )
  NALLOCATE( Ato  , (1:Mvol,-1:2,1:mn)    )
  NALLOCATE( Azo  , (1:Mvol,-1:2,1:mn)    )
  
  SALLOCATE( Fso  , (1:Mvol,     1:mn), 0 ) ! these will become redundant if/when Lagrange multipliers are used to enforce bounday constraints; 26 Jan 16;
  SALLOCATE( Fse  , (1:Mvol,     1:mn), 0 )
  
  SALLOCATE( Lma  , (1:Mvol,     1:mn), 0 ) ! degree of freedom index; for Lagrange multiplier; 08 Feb 16;
  SALLOCATE( Lmb  , (1:Mvol,     1:mn), 0 )
  SALLOCATE( Lmc  , (1:Mvol,     1:mn), 0 ) ! only need Lmc(2:mn) ; only for NOTstellsym; 08 Feb 16;
  SALLOCATE( Lmd  , (1:Mvol,     1:mn), 0 ) ! only need Lmd(2:mn) ; only for NOTstellsym; 08 Feb 16;
  SALLOCATE( Lme  , (1:Mvol,     1:mn), 0 ) ! only need Lme(2:mn) ;
  SALLOCATE( Lmf  , (1:Mvol,     1:mn), 0 ) ! only need Lmf(2:mn) ; only for NOTstellsym; 08 Feb 16;
  SALLOCATE( Lmg  , (1:Mvol,     1:mn), 0 ) ! only need Lmg(1   ) ;
  SALLOCATE( Lmh  , (1:Mvol,     1:mn), 0 ) ! only need Lmh(1   ) ;
  
  do vvol = 1, Mvol
   
   LREGION(vvol)
   
   if( Lcoordinatesingularity ) then 
    zerdof = 0                                       ! count Zernike degree of freedom 30 Jun 19
    do ii = 1, Mpol                                  ! for m>0
     do jj = ii, Lrad(vvol), 2
      zerdof = zerdof + 2 * ntor + 1                 ! plus and minus sign for n>1, unique for n==0
      if( NOTstellsym ) zerdof = zerdof + 2*ntor + 1 ! plus and minus sign for n
     enddo
    enddo
    do jj = 0, Lrad(vvol), 2                         ! for m==0
     zerdof = zerdof + ntor + 1                      ! minus sign for n
    enddo
                                     !                                     a    c      b        d      e      f      g   h
    if( YESstellsym ) NAdof(vvol) = 2 * zerdof                           + mn        + Ntor+1        + mn-1        + 1 + 0
    if( NOTstellsym ) NAdof(vvol) = 2 * zerdof                           + mn + mn-1 + Ntor+1 + Ntor + mn-1 + mn-1 + 1 + 0
   else ! .not.Lcoordinatesingularity;                                     a    c      b        d      e      f      g   h
    if( YESstellsym ) NAdof(vvol) = 2 * ( mn        ) * ( Lrad(vvol)+1 ) + mn        + mn            + mn-1        + 1 + 1  
    if( NOTstellsym ) NAdof(vvol) = 2 * ( mn + mn-1 ) * ( Lrad(vvol)+1 ) + mn + mn-1 + mn     + mn-1 + mn-1 + mn-1 + 1 + 1
   endif ! end of if( Lcoordinatesingularity );
   
   do ii = 1, mn ! loop over Fourier harmonics;
    
    do ideriv = -1, 2 ! loop over derivatives; 14 Jan 13;
     
     SALLOCATE( Ate(vvol,ideriv,ii)%s, (0:Lrad(vvol)), zero )
     SALLOCATE( Aze(vvol,ideriv,ii)%s, (0:Lrad(vvol)), zero )
     SALLOCATE( Ato(vvol,ideriv,ii)%s, (0:Lrad(vvol)), zero )
     SALLOCATE( Azo(vvol,ideriv,ii)%s, (0:Lrad(vvol)), zero )
     
    enddo ! end of do ideriv;
    
    ;  ideriv =  0
    
     SALLOCATE( Ate(vvol,ideriv,ii)%i, (0:Lrad(vvol)), 0 ) ! degree of freedom index; 17 Jan 13;
     SALLOCATE( Aze(vvol,ideriv,ii)%i, (0:Lrad(vvol)), 0 )
     SALLOCATE( Ato(vvol,ideriv,ii)%i, (0:Lrad(vvol)), 0 )
     SALLOCATE( Azo(vvol,ideriv,ii)%i, (0:Lrad(vvol)), 0 )
    
   enddo ! end of do ii;


   select case( Linitgues ) ! for iterative solver of the Beltrami fields, an initial guess is required; 11 Mar 16;
   case( 0 )    ; 
   case( 1 )    ; Ate(vvol,0,1)%s(0:1) = dtflux(vvol) * half ! this is an integrable approximation; NEEDS CHECKING; 26 Feb 13;
    ;           ; Aze(vvol,0,1)%s(0:1) = dpflux(vvol) * half ! this is an integrable approximation; NEEDS CHECKING; 26 Feb 13;
   case( 2 )    ;                                            ! will call ra00aa below to read initial vector potential from file;
   case( 3 )    ;                                            ! the initial guess will be randomized, maximum is maxrndgues; 5 Mar 19;
    do ii = 1, mn ! loop over Fourier harmonics;
    
     do ideriv = -1, 2 ! loop over derivatives; 14 Jan 13;

      call random_number(Ate(vvol,ideriv,ii)%s)
      call random_number(Aze(vvol,ideriv,ii)%s)
      Ate(vvol,ideriv,ii)%s = Ate(vvol,ideriv,ii)%s * maxrndgues
      Aze(vvol,ideriv,ii)%s = Aze(vvol,ideriv,ii)%s * maxrndgues
      if (.not. YESstellsym) then
       call random_number(Ato(vvol,ideriv,ii)%s)
       call random_number(Azo(vvol,ideriv,ii)%s)
       Ato(vvol,ideriv,ii)%s = Ato(vvol,ideriv,ii)%s * maxrndgues
       Azo(vvol,ideriv,ii)%s = Azo(vvol,ideriv,ii)%s * maxrndgues
      endif
     
     enddo ! end of do ideriv;
    
    enddo ! end of do ii;

   end select
   
   idof = 0 ! degree of freedom index; reset to 0 in each volume;
   
   if( Lcoordinatesingularity ) then
    
    do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
     
     do ll = 0, Lrad(vvol)
      if (ll>=mi .and. mod(mi+ll,2)==0)then ; idof = idof + 1 ; Ate(vvol,0,ii)%i(ll) = idof ! Zernike 30 Jun 19
      ;                                     ; idof = idof + 1 ; Aze(vvol,0,ii)%i(ll) = idof
      if( NOTstellsym .and. ii.gt.1 ) then  ; idof = idof + 1 ; Ato(vvol,0,ii)%i(ll) = idof
       ;                                    ; idof = idof + 1 ; Azo(vvol,0,ii)%i(ll) = idof
      endif ! NOTstellsym
      endif ! Zernike 
     enddo ! end of do ll; 17 Jan 13;
     ;                                     ; idof = idof + 1 ; Lma(vvol,  ii)       = idof
     if(  mi.eq.0                   ) then ; idof = idof + 1 ; Lmb(vvol,  ii)       = idof ! 18 May 16;
     endif
     if(  ii.gt.1                   ) then ; idof = idof + 1 ; Lme(vvol,  ii)       = idof
     endif
     if(  ii.eq.1                   ) then ; idof = idof + 1 ; Lmg(vvol,  ii)       = idof
!   ! ;                                    ; idof = idof + 1 ; Lmh(vvol,  ii)       = idof ! no constraint on poloidal flux in innermost volume; 11 Mar 16;
     endif
     if( NOTstellsym ) then
     if(  ii.gt.1                   ) then ; idof = idof + 1 ; Lmc(vvol,  ii)       = idof ! 18 May 16;
      ;                                    ; idof = idof + 1 ; Lmf(vvol,  ii)       = idof ! 18 May 16;
     endif
     if(  ii.gt.1 .and. mi.eq.0     ) then ; idof = idof + 1 ; Lmd(vvol,  ii)       = idof ! 18 May 16;
     endif
     endif ! end of if( NOTstellsym ) ; 19 Jul 16;

    enddo ! end of do ii; 25 Jan 13;
    
    FATAL( preset, idof.ne.NAdof(vvol), need to count Beltrami degrees-of-freedom more carefully  for coordinate singularity )
    
   else ! .not.Lcoordinatesingularity;
        
    do ii = 1, mn
     do ll = 0, Lrad(vvol)                 ; idof = idof + 1 ; Ate(vvol,0,ii)%i(ll) = idof
      ;                                    ; idof = idof + 1 ; Aze(vvol,0,ii)%i(ll) = idof
      if( ii.gt.1 .and. NOTstellsym ) then ; idof = idof + 1 ; Ato(vvol,0,ii)%i(ll) = idof
       ;                                   ; idof = idof + 1 ; Azo(vvol,0,ii)%i(ll) = idof
      endif
     enddo ! end of do ll; 08 Feb 16;
     ;                                     ; idof = idof + 1 ; Lma(vvol,  ii)       = idof
     ;                                     ; idof = idof + 1 ; Lmb(vvol,  ii)       = idof
     if(  ii.gt.1 .and. NOTstellsym ) then ; idof = idof + 1 ; Lmc(vvol,  ii)       = idof
      ;                                    ; idof = idof + 1 ; Lmd(vvol,  ii)       = idof
     endif
     if(  ii.gt.1                   ) then ; idof = idof + 1 ; Lme(vvol,  ii)       = idof
     endif
     if(  ii.gt.1 .and. NOTstellsym ) then ; idof = idof + 1 ; Lmf(vvol,  ii)       = idof
     endif
     if(  ii.eq.1                   ) then ; idof = idof + 1 ; Lmg(vvol,  ii)       = idof
      ;                                    ; idof = idof + 1 ; Lmh(vvol,  ii)       = idof
     endif
    enddo ! end of do ii; 25 Jan 13;
    
   !if( Wpreset ) then
   ! do ii = 1, mn
   !  do ll = 0, Lrad(vvol)
   !   write(ounit,'("preset : " 10x " : myid="i3" ; ii="i4" ; ll="i4" : Ate = "i7" ;")') myid, ii, ll, Ate(vvol,0,ii)%i(ll)
   !   write(ounit,'("preset : " 10x " : myid="i3" ; ii="i4" ; ll="i4" : Aze = "i7" ;")') myid, ii, ll, Aze(vvol,0,ii)%i(ll)
   !   write(ounit,'("preset : " 10x " : myid="i3" ; ii="i4" ; ll="i4" : Ato = "i7" ;")') myid, ii, ll, Ato(vvol,0,ii)%i(ll)
   !   write(ounit,'("preset : " 10x " : myid="i3" ; ii="i4" ; ll="i4" : Azo = "i7" ;")') myid, ii, ll, Azo(vvol,0,ii)%i(ll)
   !  enddo
   !   write(ounit,'("preset : " 10x " : myid="i3" ; ii="i4" ;    "4x" : Lma = "i7" ;")') myid, ii,     Lma(vvol,  ii)
   !   write(ounit,'("preset : " 10x " : myid="i3" ; ii="i4" ;    "4x" : Lmb = "i7" ;")') myid, ii,     Lmb(vvol,  ii)
   !   write(ounit,'("preset : " 10x " : myid="i3" ; ii="i4" ;    "4x" : Lmc = "i7" ;")') myid, ii,     Lmc(vvol,  ii)
   !   write(ounit,'("preset : " 10x " : myid="i3" ; ii="i4" ;    "4x" : Lmd = "i7" ;")') myid, ii,     Lmd(vvol,  ii)
   !   write(ounit,'("preset : " 10x " : myid="i3" ; ii="i4" ;    "4x" : Lme = "i7" ;")') myid, ii,     Lme(vvol,  ii)
   !   write(ounit,'("preset : " 10x " : myid="i3" ; ii="i4" ;    "4x" : Lmf = "i7" ;")') myid, ii,     Lmf(vvol,  ii)
   !   write(ounit,'("preset : " 10x " : myid="i3" ; ii="i4" ;    "4x" : Lmg = "i7" ;")') myid, ii,     Lmg(vvol,  ii)
   !   write(ounit,'("preset : " 10x " : myid="i3" ; ii="i4" ;    "4x" : Lmh = "i7" ;")') myid, ii,     Lmh(vvol,  ii)
   ! enddo
   !endif
    
    FATAL( preset, idof.ne.NAdof(vvol), need to count degrees-of-freedom more carefully for new matrix )
    
   endif ! end of if( Lcoordinatesingularity ) ; 
   
   FATAL( preset, idof.ne.NAdof(vvol), impossible logic )
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
  enddo ! end of do vvol = 1, Nvol loop;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Linitgues.eq.2 ) then ; WCALL( preset, ra00aa, ('R') )  ! read initial guess for Beltrami field from file; 02 Jan 15;
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( myid.eq.0 ) then ! 17 Oct 12;
   cput = GETTIME
   write(ounit,'("preset : ", 10x ," : ")')         
   write(ounit,'("preset : ",f10.2," : Nquad="i4" ; mn="i5" ; NGdof="i6" ; NAdof="16(i6",")" ...")') cput-cpus, Nquad, mn, NGdof, NAdof(1:min(Mvol,16))
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \subsubsection{\type{workspace}}

! Fourier transforms;

  Nt = max( Ndiscrete*4*Mpol, 1 ) ; Nz = max( Ndiscrete*4*Ntor, 1 ) ; Ntz = Nt*Nz ; soNtz = one / sqrt( one*Ntz ) ! exaggerated discrete resolution;

  ;                  ; hNt = Nt / 2
  if( Nz.gt.1 ) then ; hNz = Nz / 2
  else               ; hNz = 0
  endif
  
  if( myid.eq.0 ) then ! 17 Oct 12;
   cput = GETTIME
   write(ounit,'("preset : ", 10x ," : ")')         
   write(ounit,'("preset : ",f10.2," : Nt="i6" ; Nz="i6" ; Ntz="i9" ;")') cput-cpus, Nt, Nz, Ntz
  endif

  SALLOCATE( iRij, (1:Ntz,0:Mvol), zero ) ! interface geometry in real space; ! 18 Jul 14;
  SALLOCATE( iZij, (1:Ntz,0:Mvol), zero ) ! 
  SALLOCATE( dRij, (1:Ntz,1:Mvol), zero ) ! interface geometry in real space; poloidal derivative; ! 18 Jul 14;
  SALLOCATE( dZij, (1:Ntz,1:Mvol), zero )
  SALLOCATE( tRij, (1:Ntz,0:Mvol), zero ) ! interface geometry in real space; poloidal derivative; ! 18 Jul 14;
  SALLOCATE( tZij, (1:Ntz,0:Mvol), zero )

  SALLOCATE(   Rij, (1:Ntz,0:3,0:3    ), zero ) ! these are used for inverse fft to reconstruct real space geometry from interpolated Fourier harmonics;
  SALLOCATE(   Zij, (1:Ntz,0:3,0:3    ), zero )
  SALLOCATE(   sg , (1:Ntz,0:3        ), zero )
  SALLOCATE( guvij, (1:Ntz,0:3,0:3,0:3), zero ) ! need this on higher resolution grid for accurate Fourier decomposition;
  SALLOCATE( gvuij, (1:Ntz,0:3,0:3    ), zero ) ! need this on higher resolution grid for accurate Fourier decomposition; 10 Dec 15;
  
  SALLOCATE( dRadR, (1:mn,0:1,0:1,1:mn), zero ) ! calculated in rzaxis; 19 Sep 16;
  SALLOCATE( dRadZ, (1:mn,0:1,0:1,1:mn), zero )
  SALLOCATE( dZadR, (1:mn,0:1,0:1,1:mn), zero )
  SALLOCATE( dZadZ, (1:mn,0:1,0:1,1:mn), zero )

  SALLOCATE( dRodR, (1:Ntz,0:1,1:mn), zero ) ! calculated in rzaxis; 19 Sep 16;
  SALLOCATE( dRodZ, (1:Ntz,0:1,1:mn), zero )
  SALLOCATE( dZodR, (1:Ntz,0:1,1:mn), zero )
  SALLOCATE( dZodZ, (1:Ntz,0:1,1:mn), zero )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \subsubsection{\type{goomne, goomno} : metric information}
!latex \subsubsection{\type{gssmne, gssmno} : metric information}
!latex \subsubsection{\type{gstmne, gstmno} : metric information}
!latex \subsubsection{\type{gszmne, gszmno} : metric information}
!latex \subsubsection{\type{gttmne, gttmno} : metric information}
!latex \subsubsection{\type{gtzmne, gtzmno} : metric information}
!latex \subsubsection{\type{gzzmne, gzzmno} : metric information}
  
!latex \begin{enumerate}
!latex \item The metric information are:
!latex \bi
!latex \item[ ] \type{goomne(0:mne)}, \type{goomno(0:mne)}
!latex \item[ ] \type{gssmne(0:mne)}, \type{gssmno(0:mne)}
!latex \item[ ] \type{gstmne(0:mne)}, \type{gstmno(0:mne)}
!latex \item[ ] \type{gszmne(0:mne)}, \type{gszmno(0:mne)}
!latex \item[ ] \type{gttmne(0:mne)}, \type{gttmno(0:mne)}
!latex \item[ ] \type{gtzmne(0:mne)}, \type{gtzmno(0:mne)}
!latex \item[ ] \type{gzzmne(0:mne)}, \type{gzzmno(0:mne)}
!latex \ei
!latex \item These are defined in \link{metrix}, and used in \link{ma00aa}.
!latex \end{enumerate}
  
  SALLOCATE( goomne, (0:mne), zero ) ! workspace for Fourier decomposition of metric terms;
  SALLOCATE( goomno, (0:mne), zero )
  SALLOCATE( gssmne, (0:mne), zero ) ! workspace for Fourier decomposition of metric terms;
  SALLOCATE( gssmno, (0:mne), zero )
  SALLOCATE( gstmne, (0:mne), zero ) ! workspace for Fourier decomposition of metric terms;
  SALLOCATE( gstmno, (0:mne), zero )
  SALLOCATE( gszmne, (0:mne), zero ) ! workspace for Fourier decomposition of metric terms;
  SALLOCATE( gszmno, (0:mne), zero )
  SALLOCATE( gttmne, (0:mne), zero ) ! workspace for Fourier decomposition of metric terms;
  SALLOCATE( gttmno, (0:mne), zero )
  SALLOCATE( gtzmne, (0:mne), zero ) ! workspace for Fourier decomposition of metric terms;
  SALLOCATE( gtzmno, (0:mne), zero )
  SALLOCATE( gzzmne, (0:mne), zero ) ! workspace for Fourier decomposition of metric terms;
  SALLOCATE( gzzmno, (0:mne), zero )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  SALLOCATE( ijreal, (1:Ntz), zero ) ! real space grid;
  SALLOCATE( ijimag, (1:Ntz), zero )
  SALLOCATE( jireal, (1:Ntz), zero )
  SALLOCATE( jiimag, (1:Ntz), zero )

  SALLOCATE( jkreal, (1:Ntz), zero )
  SALLOCATE( jkimag, (1:Ntz), zero )
  SALLOCATE( kjreal, (1:Ntz), zero )
  SALLOCATE( kjimag, (1:Ntz), zero )

  SALLOCATE( cplxin,  (1:Nt,1:Nz), zero )
  SALLOCATE( cplxout, (1:Nt,1:Nz), zero )

  ! Create and save optimal plans for forward and inverse 2D fast Fourier transforms with FFTW. -JAB; 25 Jul 2017
  planf = fftw_plan_dft_2d( Nz, Nt, cplxin, cplxout, FFTW_FORWARD,  FFTW_MEASURE + FFTW_DESTROY_INPUT )
  planb = fftw_plan_dft_2d( Nz, Nt, cplxin, cplxout, FFTW_BACKWARD, FFTW_MEASURE + FFTW_DESTROY_INPUT )

  SALLOCATE( efmn, (1:mne), zero ) ! Fourier harmonics workspace; 24 Apr 13;
  SALLOCATE( ofmn, (1:mne), zero )
  SALLOCATE( cfmn, (1:mne), zero )
  SALLOCATE( sfmn, (1:mne), zero )
  SALLOCATE( evmn, (1:mne), zero )
  SALLOCATE( odmn, (1:mne), zero )
  SALLOCATE( comn, (1:mne), zero )
  SALLOCATE( simn, (1:mne), zero )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{\type{cosi(1:Ntz,1:mn)} and \type{sini(1:Ntz,1:mn)}}

!latex \begin{enumerate}
!latex \item Trigonometric factors used in various Fast Fourier transforms, where
!latex       \be \mbox{\type{cosi}}_{j,i} & = & \cos( m_i \t_j - n_i \z_j ), \\
!latex           \mbox{\type{sini}}_{j,i} & = & \sin( m_i \t_j - n_i \z_j ).
!latex       \ee
!latex \end{enumerate}

  SALLOCATE( gteta, (1:Ntz), zero )
  SALLOCATE( gzeta, (1:Ntz), zero )
  
  SALLOCATE( cosi, (1:Ntz,1:mn), zero )
  SALLOCATE( sini, (1:Ntz,1:mn), zero )

  FATAL( preset, Nz.eq.0, illegal division )
  FATAL( preset, Nt.eq.0, illegal division )

  do ii = 1, mn ; mi = im(ii) ; ni = in(ii) ! loop over Fourier harmonics;
   
   do kk = 0, Nz-1 ; zeta = kk * pi2nfp / Nz
    do jj = 0, Nt-1 ; teta = jj * pi2    / Nt ; jk = 1 + jj + kk*Nt ; arg = mi * teta - ni * zeta 
     gteta(jk) = teta
     gzeta(jk) = zeta
     cosi(jk,ii) = cos(arg)
     sini(jk,ii) = sin(arg)
    enddo
   enddo
   
  enddo ! end of do ii; 13 May 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  
  if( Wpreset .and. myid.eq.0 ) then
   
   write(ounit,'("preset : ",10x," : checking FFT and inverse FFT ;")') 
   
   do imn = 1, mn ; mm = im(imn) ; nn = in(imn) ! in should include the Nfp factor; SRH: 27 Feb 18;
    
    ijreal(1:Ntz) = zero ; ijimag(1:Ntz) = zero
    
    do kk = 0, Nz-1 ; zeta = kk * pi2nfp / Nz
     
     do jj = 0, Nt-1 ; teta = jj * pi2    / Nt ; jk = 1 + jj + kk*Nt
      
      ijreal(jk) = cos( mm * teta - nn * zeta ) ; ijimag(jk) = sin( mm * teta - nn * zeta )
      
     enddo ! end of do jj; SRH: 27 Feb 18;
     
    enddo ! end of do kk; SRH: 27 Feb 18;
    
    jkreal = ijreal ; jkimag = ijimag

    ifail = 0 !                                                              even        odd         cos         sin
    call tfft( Nt, Nz, ijreal(1:Ntz), ijimag(1:Ntz), mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), ifail )
    
    do ii = 1, mn
     
     if( abs(efmn(ii))+abs(ofmn(ii))+abs(cfmn(ii))+abs(sfmn(ii)).gt.small ) write(ounit,2000) mm, nn, im(ii), in(ii), efmn(ii), ofmn(ii), cfmn(ii), sfmn(ii)
     
2000 format("preset : ",10x," : (",i3,",",i3," ) = (",i3,",",i3," ) : "2f15.5" ; "2f15.5" ;")
     
    enddo ! end of do ii; SRH: 27 Feb 18;
    
    call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, jireal(1:Ntz), jiimag(1:Ntz) )
    
    error = ( sum((jkreal(1:Ntz)-jireal(1:Ntz))**2) + sum((jkimag(1:Ntz)-jiimag(1:Ntz))**2) ) / Ntz

    write(ounit,'("preset : ",10x," : (",i3,",",i3," ) : error = ",es13.5," ;")') mm, nn, error

   enddo ! end of do imn; SRH: 27 Feb 18;

  endif ! end of if( myid.eq.0 ) ; SRH: 27 Feb 18;

#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( Igeometry.eq.3 .and. iRbc(1,0).lt.small ) then ! have not yet assigned coordinate axis; see global;readin for user-supplied Rac, Zas, etc. ; 19 Jul 16;
   
   select case( Linitialize )
   case( :-1 ) ; vvol = Nvol + Linitialize
   case(   0 ) ; vvol =    1 ! this is really a dummy; no interpolation of interface geometry is required; packxi calls rzaxis with lvol=1; 19 Jul 16;
   case(   1 ) ; vvol = Nvol
   case(   2 ) ; vvol = Mvol
   end select

   WCALL( preset, rzaxis, ( Mvol, mn, iRbc(1:mn,0:Mvol), iZbs(1:mn,0:Mvol), iRbs(1:mn,0:Mvol), iZbc(1:mn,0:Mvol), vvol ) ) ! set coordinate axis; 19 Jul 16;

  endif ! end of if( Igeometry.eq.3 ) then ; 19 Jul 16;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{\type{psifactor(1:mn,1:Mvol)} : coordinate ``pre-conditioning'' factor}
  
!latex \begin{enumerate}

!latex \item In toroidal geometry, the coordinate ``pre-conditioning'' factor is 
!latex       \be f_{j,v} \equiv \left\{ 
!latex       \begin{array}{lcccccc}\psi_{t,v}^{0    }&,&\mbox{for $m_j=0$}, \\ 
!latex                             \psi_{t,v}^{m_j/2}&,&\mbox{otherwise}.
!latex       \end{array}\right.
!latex       \ee
!latex       where $\psi_{t,v} \equiv $ \type{tflux} is the (normalized?) toroidal flux enclosed by the $v$-th interface.

!latex \item \type{psifactor} is used in \link{packxi}, \link{dforce} and \link{hesian}.

!latex \item \type{inifactor} is similarly constructed, with
!latex       \be f_{j,v} \equiv \left\{
!latex       \begin{array}{lcccccc}\psi_{t,v}^{ 1 /2}&,&\mbox{for $m_j=0$}, \\
!latex                             \psi_{t,v}^{m_j/2}&,&\mbox{otherwise}.
!latex       \end{array}\right.
!latex       \ee
!latex       and used only for the initialization of the surfaces taking into account axis information if provided.

!latex \end{enumerate}

  SALLOCATE( psifactor, (1:mn,1:Mvol), zero )
  SALLOCATE( inifactor, (1:mn,1:Mvol), zero )

  psifactor(1:mn,1:Mvol) = one
  inifactor(1:mn,1:Mvol) = one
  
  select case( Igeometry )
   
  case( 1 )
   
   psifactor(1:mn,1:Nvol) = one
   
  case( 2 )
   
   do vvol = 1, Nvol
    do ii = 1, mn
     if( im(ii).eq.0 ) then ; psifactor(ii,vvol) = tflux(vvol)**(          +half) ! 28 Jan 15;
     else                   ; psifactor(ii,vvol) = tflux(vvol)**(halfmm(ii)-half) ! 28 Jan 15;
     endif
    enddo
   enddo
   
  case( 3 ) 
   
   do vvol = 1, Nvol
    do ii = 1, mn
     if( im(ii).eq.0 ) then ; psifactor(ii,vvol) = Rscale * tflux(vvol)**zero       ! 08 Feb 16;
                            ; inifactor(ii,vvol) = Rscale * tflux(vvol)**half       ! 17 Dec 18;
     else                   ; psifactor(ii,vvol) = Rscale * tflux(vvol)**halfmm(ii) ! 29 Apr 15;
                            ; inifactor(ii,vvol) = Rscale * tflux(vvol)**halfmm(ii) ! 17 Dec 18
     endif
    enddo
   enddo
   
  case default
   
   FATAL( readin, .true., invalid Igeometry for construction of psifactor )
   
  end select
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Linitialize.ne.0 ) then ! interpolate / extrapolate interior interface geometry; 19 Jul 16;
   
   select case( Igeometry )
    
   case( 1 ) ! Cartesian; 29 Apr 14;
    
   !FATAL( preset, Linitialize.ne.1, geometrical initialization under construction for Cartesian ) ! 14 Apr 17;
    
    do vvol = 1, Nvol
     ;iRbc(1:mn,vvol) = iRbc(1:mn,Mvol) * tflux(vvol) / tflux(Mvol) ! 14 Apr 17;
     if( NOTstellsym ) then
      iRbs(2:mn,vvol) = iRbs(2:mn,Mvol) * tflux(vvol) / tflux(Mvol) ! 14 Apr 17;
     endif
    enddo
    
   case( 2 ) ! cylindrical - standard; 20 Apr 13;
    
    FATAL( preset, Linitialize.ne.1, geometrical initialization under construction for cylindrical )
    
    do vvol = 1, Nvol-1
     ;iRbc(1:mn,vvol) = iRbc(1:mn,Nvol) * psifactor(1:mn,vvol)
     if( NOTstellsym ) then
      iRbs(2:mn,vvol) = iRbs(2:mn,Nvol) * psifactor(2:mn,vvol)
     endif
    enddo
    
   case( 3 ) ! toroidal; 20 Apr 13;
    
    FATAL( preset, Linitialize.lt.0, geometrical initialization under construction for toroidal ) ! see commented-out source below; 19 Jul 16;
    
    lvol = Nvol-1 + Linitialize

    FATAL( preset, lvol.gt.Mvol, perhaps illegal combination of Linitialize and Lfreebound )
    
!    do vvol = 1, Nvol-1       ! 19 Jul 16;
!     ;iRbc(1:mn,vvol) = iRbc(1:mn,0) + ( iRbc(1:mn,Nvol) - iRbc(1:mn,0) ) * psifactor(1:mn,vvol) ! 19 Jul 16;
!     ;iZbs(2:mn,vvol) = iZbs(2:mn,0) + ( iZbs(2:mn,Nvol) - iZbs(2:mn,0) ) * psifactor(2:mn,vvol) ! 19 Jul 16;
!     if( NOTstellsym ) then ! 19 Jul 16;
!      iRbs(2:mn,vvol) = iRbs(2:mn,0) + ( iRbs(2:mn,Nvol) - iRbs(2:mn,0) ) * psifactor(2:mn,vvol) ! 19 Jul 16;
!      iZbc(1:mn,vvol) = iZbc(1:mn,0) + ( iZbc(1:mn,Nvol) - iZbc(1:mn,0) ) * psifactor(1:mn,vvol) ! 19 Jul 16;
!     endif ! 19 Jul 16;
!    enddo
!    
    do vvol = 1, lvol-1
     ;iRbc(1:mn,vvol) = iRbc(1:mn,0) + ( iRbc(1:mn,lvol) - iRbc(1:mn,0) ) * ( inifactor(1:mn,vvol) / Rscale ) / tflux(lvol)**halfmm(1:mn)
     ;iZbs(2:mn,vvol) = iZbs(2:mn,0) + ( iZbs(2:mn,lvol) - iZbs(2:mn,0) ) * ( inifactor(2:mn,vvol) / Rscale ) / tflux(lvol)**halfmm(2:mn)
     if( NOTstellsym ) then
      iRbs(2:mn,vvol) = iRbs(2:mn,0) + ( iRbs(2:mn,lvol) - iRbs(2:mn,0) ) * ( inifactor(2:mn,vvol) / Rscale ) / tflux(lvol)**halfmm(2:mn)
      iZbc(1:mn,vvol) = iZbc(1:mn,0) + ( iZbc(1:mn,lvol) - iZbc(1:mn,0) ) * ( inifactor(1:mn,vvol) / Rscale ) / tflux(lvol)**halfmm(1:mn)
     endif
    enddo
    
!   do vvol = 1, Nvol+Linitialize-1
!    ;iRbc(1:mn,vvol) = iRbc(1:mn,0) + ( iRbc(1:mn,Nvol+Linitialize) - iRbc(1:mn,0) ) * psifactor(1:mn,vvol) / tflux(Nvol+Linitialize)**halfmm(1:mn)
!    ;iZbs(2:mn,vvol) = iZbs(2:mn,0) + ( iZbs(2:mn,Nvol+Linitialize) - iZbs(2:mn,0) ) * psifactor(2:mn,vvol) / tflux(Nvol+Linitialize)**halfmm(2:mn)
!    if( NOTstellsym ) then
!     iRbs(2:mn,vvol) = iRbs(2:mn,0) + ( iRbs(2:mn,Nvol+Linitialize) - iRbs(2:mn,0) ) * psifactor(2:mn,vvol) / tflux(Nvol+Linitialize)**halfmm(2:mn)
!     iZbc(1:mn,vvol) = iZbc(1:mn,0) + ( iZbc(1:mn,Nvol+Linitialize) - iZbc(1:mn,0) ) * psifactor(1:mn,vvol) / tflux(Nvol+Linitialize)**halfmm(1:mn)
!    endif
!   enddo
    
   end select ! matches select case( Igeometry ); 19 Jul 16;
    
  endif ! matches if( Linitialize.ne.0 ) then; 19 Jul 16;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{\type{Bsupumn} and \type{Bsupvmn}}

  SALLOCATE( Bsupumn, (1:Nvol,0:1,1:mn), zero ) ! Fourier components of {\bf B}\cdot\nabla \theta on boundary; required for virtual casing;
  SALLOCATE( Bsupvmn, (1:Nvol,0:1,1:mn), zero ) ! Fourier components of {\bf B}\cdot\nabla \zeta  on boundary;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{\type{diotadxup} and \type{glambda} : transformation to straight fieldline angle}

!latex \begin{enumerate}
!latex \item Given the Beltrami fields in any volume, the rotational-transform on the adjacent interfaces 
!latex       may be determined (in \link{tr00ab}) by constructing the straight fieldline angle on the interfaces.
!latex \item The rotational transform on the inner or outer interface of a given volume depends on the magnetic field in that volume,
!latex       i.e. $\iotabar_\pm = \iotabar({\bf B}_\pm)$,
!latex       so that
!latex       \be \delta \iotabar_\pm = \frac{\partial \iotabar_\pm}{\partial {\bf B}_\pm} \cdot \delta {\bf B_\pm}.
!latex       \ee
!latex \item The magnetic field depends on the Fourier harmonics of both the inner and outer interface geometry (represented here as $x_j$),
!latex       the helicity multiplier, and the enclosed poloidal flux, i.e. ${\bf B_\pm} = {\bf B_\pm}(x_j, \mu, \Delta \psi_p)$, so that 
!latex       \be \delta {\bf B_\pm} = \frac{\partial {\bf B}_\pm}{\partial x_j          } \delta x_j
!latex                              + \frac{\partial {\bf B}_\pm}{\partial \mu          } \delta \mu
!latex                              + \frac{\partial {\bf B}_\pm}{\partial \Delta \psi_p} \delta \Delta \psi_p.
!latex       \ee
!latex \item The rotational-transforms, thus, can be considered to be functions of the geometry, the helicity-multiplier and the enclosed poloidal flux,
!latex       $\iotabar_{\pm} = \iotabar_{\pm}(x_j,\mu,\Delta\psi_p)$.
!latex \item The rotational-transform, and its derivatives, on the inner and outer interfaces of each volume is stored in 
!latex       \\ \type{diotadxup(0:1,-1:2,1:Mvol)}. 
!latex       The arguments label:
!latex       \begin{enumerate}
!latex       \item[i.] the first argument labels the inner or outer interface,
!latex       \item[ii.] the the second labels derivative, with 
!latex       \begin{enumerate} \item[-1 :] indicating the derivative with respect to the interface geometry,
!latex                                 i.e. $\ds \frac{\partial \iotabar_{\pm}}{\partial x_j}$,
!latex                         \item[0 :] the rotational-transform itself,
!latex                         \item[1,2 :] the derivatives with respec to $\mu$ and $\Delta \psi_p$,
!latex                                 i.e. $\ds \frac{\partial \iotabar_{\pm}}{\partial \mu}$ and
!latex                                      $\ds \frac{\partial \iotabar_{\pm}}{\partial \Delta \psi_p}$;
!latex       \end{enumerate}
!latex       \item[iii.] the third argument labels volume.
!latex       \end{enumerate}
!latex \item The values of \type{diotadxup} are assigned in \link{mp00aa} after calling \link{tr00ab}.
!latex \end{enumerate}

  SALLOCATE( diotadxup, (0:1,-1:2,1:Mvol), zero ) ! measured rotational transform on inner/outer interfaces in each annulus;
  SALLOCATE( dItGpdxtp, (0:1,-1:2,1:Mvol), zero ) ! measured plasma and linking currents                                   ;

  SALLOCATE( glambda, (1:Ntz+1,0:2,0:1,1:Mvol), zero ) ! save initial guesses for iterative calculation of rotational-transform; 21 Apr 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! Construction of ``force'';

  SALLOCATE( Bemn, (1:mn,1:Mvol,0:1), zero )
  SALLOCATE( Bomn, (1:mn,1:Mvol,0:1), zero )
  SALLOCATE( Iomn, (1:mn,1:Mvol    ), zero )
  SALLOCATE( Iemn, (1:mn,1:Mvol    ), zero )
  SALLOCATE( Somn, (1:mn,1:Mvol,0:1), zero )
  SALLOCATE( Semn, (1:mn,1:Mvol,0:1), zero )
  SALLOCATE( Pomn, (1:mn,1:Mvol,0:2), zero )
  SALLOCATE( Pemn, (1:mn,1:Mvol,0:2), zero )

  SALLOCATE( BBe , (1:Mvol-1), zero )
  SALLOCATE( IIo , (1:Mvol-1), zero )
  SALLOCATE( BBo , (1:Mvol-1), zero )
  SALLOCATE( IIe , (1:Mvol-1), zero )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  SALLOCATE( Btemn, (1:mn,0:1,1:Mvol), zero ) ! these are declared in global, calculated in sc00aa, broadcast in xspech, and written to file in hdfint;
  SALLOCATE( Bzemn, (1:mn,0:1,1:Mvol), zero )
  SALLOCATE( Btomn, (1:mn,0:1,1:Mvol), zero )
  SALLOCATE( Bzomn, (1:mn,0:1,1:Mvol), zero )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{\type{vvolume}, \type{lBBintegral} and \type{lABintegral}}

!latex \begin{enumerate}
!latex \item volume integrals
!latex       \be \type{vvolume(i)}     &=& \int_{{\cal V}_i}                       \, dv\\
!latex           \type{lBBintegral(i)} &=& \int_{{\cal V}_i} {\bf B} \cdot {\bf B} \, dv\\
!latex           \type{lABintegral(i)} &=& \int_{{\cal V}_i} {\bf A} \cdot {\bf B} \, dv
!latex       \ee
!latex \end{enumerate}
  
  SALLOCATE( vvolume    , (1:Mvol), zero ) ! volume integral of \sqrt g;
  SALLOCATE( lBBintegral, (1:Mvol), zero ) ! volume integral of B.B    ;
  SALLOCATE( lABintegral, (1:Mvol), zero ) ! volume integral of A.B    ;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( YESstellsym ) lmns = 1 + (mns-1)           ! number of independent degrees-of-freedom in angle transformation; 30 Jan 13; 
  if( NOTstellsym ) lmns = 1 + (mns-1) + (mns-1) ! number of independent degrees-of-freedom in angle transformation; 30 Jan 13; 
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if (Lfreebound > 0) then ! Only do for free-boundary; 7 Nov 18;
  
    SALLOCATE( Dxyz, (1:3,1:Ntz), zero ) ! Cartesian components of computational boundary; position; 14 Apr 17;
    SALLOCATE( Nxyz, (1:3,1:Ntz), zero ) ! Cartesian components of computational boundary; normal  ; 14 Apr 17;

    SALLOCATE( Jxyz, (1:Ntz,1:3), zero ) ! Cartesian components of virtual casing surface current; needs to be recalculated at each iteration;
  
    lvol = Mvol ; lss = one ; Lcurvature = 1 ; Lcoordinatesingularity = .false. ! will only require normal field on outer interface = computational boundary; 
  
    WCALL( preset, coords,( lvol, lss, Lcurvature, Ntz, mn ) ) ! will need Rij, Zij; THE COMPUTATIONAL BOUNDARY DOES NOT CHANGE;
  
    do kk = 0, Nz-1 ; zeta = kk * pi2nfp / Nz
   
     if( Igeometry.eq.3 ) then ; cszeta(0:1) = (/ cos(zeta), sin(zeta) /)
     endif
   
     do jj = 0, Nt-1 ; teta = jj * pi2    / Nt ; jk = 1 + jj + kk*Nt
    
      select case( Igeometry )
      case( 1 ) ! Igeometry = 1 ;
       Dxyz(1:3,jk) = (/   teta       ,  zeta       ,   Rij(jk,0,0) /)
       Nxyz(1:3,jk) = (/ - Rij(jk,2,0), -Rij(jk,3,0),   one         /)
      case( 2 ) ! Igeometry = 2 ;
       FATAL( bnorml, .true., free-boundary calculations not yet implemented in cylindrical geometry )
      case( 3 ) ! Igeometry = 3 ;
       Dxyz(1:3,jk) = (/   Rij(jk,0,0) * cszeta(0), Rij(jk,0,0) * cszeta(1), Zij(jk,0,0) /)
       Nxyz(1:3,jk) = (/   Rij(jk,2,0) * cszeta(1) * Zij(jk,3,0) - Zij(jk,2,0) * ( Rij(jk,3,0) * cszeta(1) + Rij(jk,0,0) * cszeta(0) ), &
                         - Rij(jk,2,0) * cszeta(0) * Zij(jk,3,0) + Zij(jk,2,0) * ( Rij(jk,3,0) * cszeta(0) - Rij(jk,0,0) * cszeta(1) ), &
                           Rij(jk,0,0)             * Rij(jk,2,0) /)         
      end select ! end of select case( Igeometry ) ; 09 Mar 17;
    
     enddo ! end of do jj; 14 Apr 17;
   
    enddo ! end of do kk; 14 Apr 17;
  
  endif ! Lfreebound > 1; 7 Nov 18;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(preset)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine preset

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
