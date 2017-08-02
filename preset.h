!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (initialization) ! Allocates and initializes internal arrays.

!latex \briefly{Allocates and initializes internal arrays.}

!latex \calledby{\link{xspech}}
!latex \calls{\link{ra00aa}}

!latex \tableofcontents

!latex \subsection{definition of internal variables}

subroutine preset
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero
  
  use numerical, only : sqrtmachprec, vsmall, small
  
  use fileunits, only : ounit
  
  use inputlist! only :
  
  use cputiming, only : Tpreset
  
  use allglobal! only :
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS
  
  INTEGER   :: innout, idof, jk, ll, pp, ii, jj, ij, ifail, ideriv, vvol, mi, ni, mj, nj, mk, nk, mimj, ninj, mkmj, nknj, kk, lvol
  INTEGER   :: itype, lquad, id01bcf, maxIquad, Mrad, jquad, Lcurvature
  REAL      :: teta, zeta, arg, lss, aa, bb, cc, dd, cszeta(0:1)
  
  BEGIN(preset)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  FATAL( preset, Nfp.eq.0, illegal division )

  pi2nfp         = pi2 / Nfp
  
  pi2pi2nfp      = pi2 * pi2nfp
  pi2pi2nfphalf  = pi2 * pi2nfp * half  ; opi2pi2nfphalf = one / pi2pi2nfphalf ! SRH; 01 Aug 17;
  pi2pi2nfpquart = pi2 * pi2nfp * quart

  Mrad  = maxval( Lrad(1:Mvol) )

  if( myid.eq.0 ) write(ounit,'("preset : " 10x " : myid="i3" ; Mrad="i3" : Lrad="257i3)') myid, Mrad, Lrad(1:Mvol)
  
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
  
  if( Wpreset ) then ; cput = GETTIME ; write(ounit,'("preset : ",f10.2," : myid=",i3," ; NGdof="i9" ;")') cput-cpus, myid, NGdof
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
   
1002 format("preset : " 10x " :      "3x" ; transform : "i3" : ("i3" /"i3" ) * ("i3" /"i3" ) = "f18.15" ; ("i3" /"i3" ) * ("i3" /"i3" ) = "f18.15" ;")
   
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
!latex       \be \type{sweight}_v \equiv \inputvar{upsilon} \times \psi_{v}^w,
!latex       \ee
!latex       where $\psi_v \equiv $ \type{tflux(v)} is the normalized toroidal flux enclosed by the $v$-th interface,
!latex       and $w \equiv $ \inputvar{wpoloidal}.
!latex \end{enumerate}
  
  SALLOCATE( sweight, (1:Mvol), zero )
  sweight(1:Mvol) = upsilon * tflux(1:Mvol)**wpoloidal ! found an error (2:Mvol) --> (1:Mvol) ; 19 Jul 16;

!#ifdef DEBUG
!  write(ounit,'("preset : " 10x " : sweight="99(es13.5,","))') sweight(1:Mvol)
!#endif
  
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
  
  do innout = 0, 1 ; lss = two * innout - one
   
   do ll = 0, Mrad ; TT(ll,innout,0) = lss**(ll  )        
    ;              ; TT(ll,innout,1) = lss**(ll+1) * ll**2 ! derivative; 26 Jan 16;
   enddo
   
  enddo ! end of do innout = 0, 1 ; 20 Jun 14;
  
 !write(ounit,'("preset : " 10x " : myid="i3" : TT(0:M,0,0)="999f7.2)') myid, TT(0:Mrad,0,0)
 !write(ounit,'("preset : " 10x " : myid="i3" : TT(0:M,1,0)="999f7.2)') myid, TT(0:Mrad,1,0)
 !write(ounit,'("preset : " 10x " : myid="i3" : TT(0:M,0,1)="999f7.2)') myid, TT(0:Mrad,0,1)
 !write(ounit,'("preset : " 10x " : myid="i3" : TT(0:M,1,1)="999f7.2)') myid, TT(0:Mrad,1,1)

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
!latex \item Consider the `abbreviated' representation for a double Fourier series,
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
    
   enddo ! end of do kk = 1, mne; 03 Apr 13;
   
  enddo ! end of do ii; 17 Dec 15;

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
   
   SALLOCATE( djkp, (1:mn,1:mn), 0 ) ! only used in volume; trignometric identities; 04 Dec 14;
   SALLOCATE( djkm, (1:mn,1:mn), 0 ) ! only used in volume; trignometric identities; 04 Dec 14;
   
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

  SALLOCATE( cheby    , (0:Mrad,0:2         ), zero ) ! 27 Jul 17;

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
!latex       which are computed using \nag{www.nag.co.uk/numeric/FL/manual19/pdf/D01/d01bcf_fl19.pdf}{D01BCF}.
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

  SALLOCATE( TD, (0:Mrad,0:1,1:maxIquad,1:Mvol), zero ) ! only used in ma00aa; 27 Jul 17;

  do vvol = 1, Mvol
   
   itype = 0 ; aa = -one ; bb = +one ; cc = zero ; dd = zero

   lquad = Iquad(vvol)
   
   id01bcf = 1
   call D01BCF( itype, aa, bb, cc, dd, lquad, gaussianweight(1:lquad,vvol), gaussianabscissae(1:lquad,vvol), id01bcf ) ! sets gaussian weights & abscissae;
   
   if( myid.eq.0 ) then
   cput= GETTIME
   select case( id01bcf ) !                                                        123456789012345
   case( 0 )    ; if( Wpreset ) write(ounit,1000) cput-cpus, myid, vvol, id01bcf, "success        ", gaussianabscissae(1:lquad,vvol)
   case( 1 )    ;               write(ounit,1000) cput-cpus, myid, vvol, id01bcf, "failed         ", gaussianabscissae(1:lquad,vvol)
   case( 2 )    ;               write(ounit,1000) cput-cpus, myid, vvol, id01bcf, "input error    ", gaussianabscissae(1:lquad,vvol)
   case( 3 )    ;               write(ounit,1000) cput-cpus, myid, vvol, id01bcf, "input error    ", gaussianabscissae(1:lquad,vvol)
   case( 4 )    ;               write(ounit,1000) cput-cpus, myid, vvol, id01bcf, "weight overflow", gaussianabscissae(1:lquad,vvol)
   case( 5 )    ;               write(ounit,1000) cput-cpus, myid, vvol, id01bcf, "weight zero    ", gaussianabscissae(1:lquad,vvol)
   case( 6 )    ;               write(ounit,1000) cput-cpus, myid, vvol, id01bcf, "failed         ", gaussianabscissae(1:lquad,vvol)
   case default ;               write(ounit,1000) cput-cpus, myid, vvol, id01bcf, "weird          ", gaussianabscissae(1:lquad,vvol)
    FATAL( preset, .true., weird ifail returned by D01BCF )
   end select
   ;              if( Wpreset ) write(ounit,1001)                                                    gaussianweight(1:lquad,vvol)
   endif
   
1000 format("preset : ",f10.2," : myid=",i3," ; lvol=",i3," ; id01bcf=",i5," ; "a15" ; abscissae ="99f09.05)
1001 format("preset : ", 10x ," :      "3x"        "3x"           "3x"   "15x" ; weights   ="99f09.05)
      
   do jquad = 1, lquad ! loop over radial sub-sub-grid (numerical quadrature); ! 27 Jul 17;
    
    lss = gaussianabscissae(jquad,vvol) ! 27 Jul 17;
    
    ;                     ; TD( 0,0:1,jquad,vvol) = (/ one, zero /) ! Chebyshev initialization; 16 Jan 13;   
    ;                     ; TD( 1,0:1,jquad,vvol) = (/ lss,  one /) ! 
    do ll = 2, Lrad(vvol) ; TD(ll,0:1,jquad,vvol) = (/ two * lss * TD(ll-1,0,jquad,vvol)                                     - TD(ll-2,0,jquad,vvol) , &
                                                       two       * TD(ll-1,0,jquad,vvol) + two * lss * TD(ll-1,1,jquad,vvol) - TD(ll-2,1,jquad,vvol) /)
    enddo
    
   enddo ! end of do jquad = 1, lquad;

  enddo ! end of do vvol;  7 Mar 13; 

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  SALLOCATE( llabel, (1:(Mrad+1)*(Mrad+1),1:Mvol), 0 ) ! SRH; 27 Jul 17;
  SALLOCATE( plabel, (1:(Mrad+1)*(Mrad+1),1:Mvol), 0 ) ! SRH; 27 Jul 17;
  
  do vvol = 1, Mvol 
   
   lp = 0
   do ll = 0, Lrad(vvol)
    do pp = 0, Lrad(vvol) ; lp = lp + 1 ; llabel(lp,vvol) = ll ; plabel(lp,vvol) = pp ! used in ma00aa; SRH; 27 Jul 17;
    enddo
   enddo
   
  enddo ! end of do vvol;  7 Mar 13; 
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  SALLOCATE( ilabel, (1:mnsqd), 0 ) ! SRH; 27 Jul 17;
  SALLOCATE( jlabel, (1:mnsqd), 0 ) ! SRH; 27 Jul 17;

  ij = 0
  do ii = 1, mn
   do jj = 1, mn ; ij = ij + 1 ; ilabel(ij) = ii ; jlabel(ij) = jj ! used in ma00aa; SRH; 27 Jul 17;
   enddo
  enddo

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
   
   if( Lcoordinatesingularity ) then !                                     a    c      b        d      e      f      g   h
    if( YESstellsym ) NAdof(vvol) = 2 * ( mn        ) * ( Lrad(vvol)+1 ) + mn        + Ntor+1        + mn-1        + 1 + 0
    if( NOTstellsym ) NAdof(vvol) = 2 * ( mn + mn-1 ) * ( Lrad(vvol)+1 ) + mn + mn-1 + Ntor+1 + Ntor + mn-1 + mn-1 + 1 + 0
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
   end select
   
   idof = 0 ! degree of freedom index; reset to 0 in each volume;
   
   if( Lcoordinatesingularity ) then
    
    do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
     
     do ll = 0, Lrad(vvol)
      ;                                     ; idof = idof + 1 ; Ate(vvol,0,ii)%i(ll) = idof
      ;                                     ; idof = idof + 1 ; Aze(vvol,0,ii)%i(ll) = idof
      if( NOTstellsym .and. ii.gt.1 ) then  ; idof = idof + 1 ; Ato(vvol,0,ii)%i(ll) = idof
       ;                                    ; idof = idof + 1 ; Azo(vvol,0,ii)%i(ll) = idof
      endif
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

  SALLOCATE( trigm , (1:2*Nt) , zero ) ! trignometric factors required for fast Fourier transform;
  SALLOCATE( trign , (1:2*Nz) , zero )
  SALLOCATE( trigwk, (1:2*Ntz), zero )

  SALLOCATE( ijreal, (1:Ntz), zero ) ! real space grid;
  SALLOCATE( ijimag, (1:Ntz), zero )
  SALLOCATE( jireal, (1:Ntz), zero )
  SALLOCATE( jiimag, (1:Ntz), zero )

  SALLOCATE( jkreal, (1:Ntz), zero )
  SALLOCATE( jkimag, (1:Ntz), zero )
  SALLOCATE( kjreal, (1:Ntz), zero )
  SALLOCATE( kjimag, (1:Ntz), zero )

  isr = 'I' ; ifail = 0

  WCALL( preset, C06FUF, ( Nt, Nz, ijreal(1:Ntz), ijimag(1:Ntz), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz), ifail ) )

  isr = 'S' ! prepare FFTs; 23 Jan 13;

  FATAL( preset, ifail.ne.0, error constructing Fourier transform )

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

! Rscale = one ! debugging; 03 Nov 16;

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

!latex \end{enumerate}

  SALLOCATE( psifactor, (1:mn,1:Mvol), zero )

  psifactor(1:mn,1:Mvol) = one
  
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
    !if( im(ii).eq.0 ) then ; psifactor(ii,vvol) = Rscale * tflux(vvol)**half       ! 29 Apr 15;
     if( im(ii).eq.0 ) then ; psifactor(ii,vvol) = Rscale * tflux(vvol)**zero       ! 08 Feb 16;
     else                   ; psifactor(ii,vvol) = Rscale * tflux(vvol)**halfmm(ii) ! 29 Apr 15;
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
     ;iRbc(1:mn,vvol) = iRbc(1:mn,0) + ( iRbc(1:mn,lvol) - iRbc(1:mn,0) ) * ( psifactor(1:mn,vvol) / Rscale ) / tflux(lvol)**halfmm(1:mn)
     ;iZbs(2:mn,vvol) = iZbs(2:mn,0) + ( iZbs(2:mn,lvol) - iZbs(2:mn,0) ) * ( psifactor(2:mn,vvol) / Rscale ) / tflux(lvol)**halfmm(2:mn)
     if( NOTstellsym ) then
      iRbs(2:mn,vvol) = iRbs(2:mn,0) + ( iRbs(2:mn,lvol) - iRbs(2:mn,0) ) * ( psifactor(2:mn,vvol) / Rscale ) / tflux(lvol)**halfmm(2:mn)
      iZbc(1:mn,vvol) = iZbc(1:mn,0) + ( iZbc(1:mn,lvol) - iZbc(1:mn,0) ) * ( psifactor(1:mn,vvol) / Rscale ) / tflux(lvol)**halfmm(1:mn)
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

! Construction of `force';

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
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef PRECALCULATE

  SALLOCATE( DSoocc, (0:Mrad,0:Mrad,1:mn,1:mn), zero )
  SALLOCATE( DSoocs, (0:Mrad,0:Mrad,1:mn,1:mn), zero )
  SALLOCATE( DSoosc, (0:Mrad,0:Mrad,1:mn,1:mn), zero )
  SALLOCATE( DSooss, (0:Mrad,0:Mrad,1:mn,1:mn), zero )

  SALLOCATE( DToocc, (0:Mrad,0:Mrad,1:mn,1:mn), zero )
  SALLOCATE( DToocs, (0:Mrad,0:Mrad,1:mn,1:mn), zero )
  SALLOCATE( DToosc, (0:Mrad,0:Mrad,1:mn,1:mn), zero )
  SALLOCATE( DTooss, (0:Mrad,0:Mrad,1:mn,1:mn), zero )
  
  do vvol = 1, Mvol ; lquad = Iquad(vvol)
   
   LREGION(vvol) ! assigns Lcoordinatesingularity, Lplasmaregion, etc. ; 03 Nov 16;
   
   if( Lcoordinatesingularity ) then ! additional radial factors, such as r^m, are included to "regularize" the magnetic field near the origin; SRH; 27 Jul 17;
    
    sbar(1:lquad) = ( gaussianabscissae(1:lquad,vvol) + one ) * half
    
    halfoversbar(1:lquad) = half / sbar(1:lquad)
    
    do jquad = 1, lquad ; sbarhim(jquad,1:mn) = sbar(jquad)**regumm(1:mn) ! pre-calculation of regularization factor; 12 Sep 13;
    enddo
    
    do jquad = 1, lquad ! Gaussian quadrature loop;
     
     lss = gaussianabscissae(jquad,vvol) ; jthweight = gaussianweight(jquad,vvol)
     
     FATAL( preset, .true., only goomne and goomno are required from metrix )

     WCALL( ma00aa, metrix,( vvol, lss ) ) ! compute metric elements; 16 Jan 13; ! THIS NEEDS TO BE REPLACED; SRH; 01 Aug 17;
     
     do ii = 1, mn
      
      do jj = 1, mn
       
       kks = kijs(ii,jj,0) ; kds = jthweight / kijs(ii,jj,1) ! SRH; 27 Jul 17;
       kka = kija(ii,jj,0) ; kda = jthweight / kija(ii,jj,1) ! SRH; 27 Jul 17;
       
       foocc = + goomne(kks) * abs(kds) + goomne(kka) * abs(kda) ! stell-sym; SRH; 27 Jul 17;
       foocs = - goomno(kks) *     kds  + goomno(kka) *     kda 
       foosc = + goomno(kks) *     kds  + goomno(kka) *     kda 
       fooss = + goomne(kks) * abs(kds) - goomne(kka) * abs(kda)
       
       do ll = 0, Lrad(vvol)
        
        do pp = 0, Lrad(vvol)
         
         Tl = sbarhim(jquad,ii) *                                      TD(ll,0,jquad,vvol)
         Dl = sbarhim(jquad,ii) * ( regumm(ii) * halfoversbar(jquad) * TD(ll,0,jquad,vvol) + TD(ll,1,jquad,vvol) )
         Tp = sbarhim(jquad,jj) *                                      TD(pp,0,jquad,vvol)
         Dp = sbarhim(jquad,jj) * ( regumm(jj) * halfoversbar(jquad) * TD(pp,0,jquad,vvol) + TD(pp,1,jquad,vvol) )
         
         TlTp = Tl * Tp
         TlDp = Tl * Dp
         DlTp = Dl * Tp
         DlDp = Dl * Dp
         
         DSoocc( ll, pp, ii, jj ) = DSoocc( ll, pp, ii, jj ) + DlTp * foocc ! stell-sym; SRH; 27 Jul 17;
         DSoocs( ll, pp, ii, jj ) = DSoocs( ll, pp, ii, jj ) + DlTp * foocs
         DSoosc( ll, pp, ii, jj ) = DSoosc( ll, pp, ii, jj ) + DlTp * foosc
         DSooss( ll, pp, ii, jj ) = DSooss( ll, pp, ii, jj ) + DlTp * fooss
         
        enddo ! end of do pp; SRH; 01 Aug 17;
        
       enddo ! end of do ll; SRH; 01 Aug 17;

      enddo ! end of do jj; SRH; 01 Aug 17;

     enddo ! end of do ii; SRH; 01 Aug 17;
     
    enddo ! end of do jquad; ! 16 Jan 13;
    
   else ! .not.Lcoordinatesingularity;
    
    if( Lrad(vvol).lt.Mrad ) cycle ! compute the integrals at the highest resolution; SRH; 01 Aug 17;

    do jquad = 1, lquad ! Gaussian quadrature loop;
     
     lss = gaussianabscissae(jquad,vvol) ; jthweight = gaussianweight(jquad,vvol)
     
     FATAL( preset, .true., only goomne and goomno are required from metrix )
     
     WCALL( ma00aa, metrix,( vvol, lss ) ) ! compute metric elements; 16 Jan 13; ! THIS NEEDS TO BE REPLACED; SRH; 01 Aug 17;
     
     do ii = 1, mn
      
      do jj = 1, mn
       
       kks = kijs(ii,jj,0) ; kds = jthweight / kijs(ii,jj,1) ! SRH; 27 Jul 17;
       kka = kija(ii,jj,0) ; kda = jthweight / kija(ii,jj,1) ! SRH; 27 Jul 17;
       
       foocc = + goomne(kks) * abs(kds) + goomne(kka) * abs(kda) ! stell-sym; SRH; 27 Jul 17;
       foocs = - goomno(kks) *     kds  + goomno(kka) *     kda 
       foosc = + goomno(kks) *     kds  + goomno(kka) *     kda 
       fooss = + goomne(kks) * abs(kds) - goomne(kka) * abs(kda)
       
       do ll = 0, Lrad(vvol)
        
        do pp = 0, Lrad(vvol)
         
         Tl = TD(ll,0,jquad,vvol)
         Dl = TD(ll,1,jquad,vvol)
         Tp = TD(pp,0,jquad,vvol)
         Dp = TD(pp,1,jquad,vvol)
         
         TlTp = Tl * Tp
         TlDp = Tl * Dp
         DlTp = Dl * Tp
         DlDp = Dl * Dp
         
         DToocc( ll, pp, ii, jj ) = DToocc( ll, pp, ii, jj ) + DlTp * foocc ! stell-sym; SRH; 27 Jul 17;
         DToocs( ll, pp, ii, jj ) = DToocs( ll, pp, ii, jj ) + DlTp * foocs
         DToosc( ll, pp, ii, jj ) = DToosc( ll, pp, ii, jj ) + DlTp * foosc
         DTooss( ll, pp, ii, jj ) = DTooss( ll, pp, ii, jj ) + DlTp * fooss
         
        enddo ! end of do pp; SRH; 01 Aug 17;
        
       enddo ! end of do ll; SRH; 01 Aug 17;
       
      enddo ! end of do jj; SRH; 01 Aug 17;
      
     enddo ! end of do ii; SRH; 01 Aug 17;
     
    enddo ! end of do jquad; SRH; 01 Aug 17;
    
   endif ! end of if( Lcoordinatesingularity ) ; SRH; 01 Aug 17;

  enddo ! end of do vvol; SRH; 01 Aug 17;

#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(preset)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine preset

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
