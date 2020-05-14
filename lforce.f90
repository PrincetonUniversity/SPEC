!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (&ldquo;local&rdquo; force) ! Computes $B^2$, and the spectral condensation constraints if required, on the interfaces, ${\cal I}_i$.

!latex \briefly{Computes $B^2$ and the spectral condensation constraints on the interfaces.}

!latex \calledby{\link{dforce}}
!latex \calls{\link{coords}}

!latex \tableofcontents
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
!latex \subsection{field strength} 

!latex \begin{enumerate}
!latex \item The field strength is given by $B^2 = B^\s B_\s + B^\t B_\t + B^\z B_\z$, and on the interfaces $B^\s=0$ by construction.
!latex \item The magnetic field is
!latex       $\sqrt g \; {\bf B} = (\partial_\t A_\z - \partial_\z A_\t ) {\bf e}_\s - \partial_\s A_\z {\bf e}_\t + \partial_\s A_\t {\bf e}_\z$.
!latex \item The covariant components of the field are computed via $B_\t = B^\t g_{\t\t} + B^\z g_{\t\z}$ and $B_\z = B^\t g_{\t\z} + B^\z g_{\z\z}$.
!latex \item The expression for $B^2$ is
!latex       \be
!latex       (\sqrt g)^2 B^2 = A_\z^\prime \; A_\z^\prime \; g_{\t\t} - 2 \; A_\z^\prime \; A_\t^\prime \; g_{\t\z} + A_\t^\prime \; A_\t^\prime \; g_{\z\z},
!latex       \ee
!latex       where the ``$\prime$'' denotes derivative with respect to $\s$.
!latex \item The quantity returned is
!latex       \be F \equiv \inputvar{pscale} \times \frac{P}{V^\gamma} + \frac{B^2}{2},
!latex       \ee
!latex       where $P\equiv$ \type{adiabatic} and $V\equiv$ volume.
!latex \end{enumerate} 

!latex \subsection{spectral constraints} 

!latex \begin{enumerate}
!latex \item In addition to the physical-force-balance constraints, namely that $[[p+B^2/2]]=0$ across the interfaces,
!latex       additional angle constraints are required to obtain a unique Fourier representation of the interface geometry.
!latex \item Introducing the angle functional: a weighted combination of the ``polar'' constraint;
!latex       the normalized, poloidal, spectral width 
!latex       [\paper{Hirshman \& Meier}{S.P. Hirshman \& H.K. Meier}{10.1063/1.864972}{Phys. Fluids}{28}{1387}{1985},
!latex       \paper{Hirshman \& Breslau}{S.P. Hirshman \& J. Breslau}{10.1063/1.872954}{Phys. Plasmas}{5}{2664}{1998}]
!latex       the poloidal-angle origin constraint;
!latex       and the ``length'' of the angle curves
!latex       \be F \!\equiv\! \sum_{i=1}^{N-1} \alpha_i \underbrace{\ooint \!\! \frac{1}{\Theta_{i,\t}}}_{polar-angle}
!latex          +     \sum_{i=1}^{N-1} \beta_i \, \underbrace{\;\; M_i \;\;}_{spectral-width}
!latex          +     \sum_{i=1}^{N-1} \gamma_i \int_{0}^{2\pi} \! \frac{1}{2}\left[Z_i(0,\z)-Z_{i,0}\right]^2 d\z
!latex          +     \red{\ooint \!\! \sum_{i=1}^{N} \delta_i \, \underbrace{\;\;\;\; L_i \;\;\;\;}_{poloidal-length}}
!latex       \ee
!latex       where $i$ labels the interfaces, and
!latex       \be \Theta_{i,\t} & \equiv & \frac{x \, y_\t - x_\t \, y}{x^2+y^2}, \\
!latex       M_i          & \equiv & \frac{\sum_j m_j^p ( R_{j,i}^2 + Z_{j,i}^2 )}{\sum_j       ( R_{j,i}^2 + Z_{j,i}^2 )} ,\\
!latex       \red{L_i}          & \equiv & \red{\sqrt{      [R_{i}(\t,\z)-R_{i-1}(\t,\z)]^2+[Z_{i}(\t,\z)-Z_{i-1}(\t,\z)]^2},}
!latex       \ee
!latex       and where $j$ labels the Fourier harmonics.
!latex       The $\alpha_i$, $\beta_i$, $\gamma_i$ and $\delta_i \equiv $ \internal{sweight} are user-supplied weight factors.
!latex \item The polar constraint is derived from defining $\tan \Theta \equiv y/x$, where
!latex       \be x(\t,\z) &  \equiv &  R_{i}(\t,\z)-R_{i,0}(\z), \\
!latex           y(\t,\z) &  \equiv &  Z_{i}(\t,\z)-Z_{i,0}(\z),
!latex       \ee
!latex       and where the geometric center of each interface is given by the arc-length weighted integrals, see \link{rzaxis},
!latex       \be
!latex           R_{i,0} & \equiv & \int_0^{2\pi} \!\!\!\! d\t \; R_i(\t,\z) \sqrt{R_{i,\t}(\t,\z)^2+Z_{i,\t}(\t,\z)^2}, \\
!latex           Z_{i,0} & \equiv & \int_0^{2\pi} \!\!\!\! d\t \; Z_i(\t,\z) \sqrt{R_{i,\t}(\t,\z)^2+Z_{i,\t}(\t,\z)^2},
!latex       \ee
!latex       and $\cos \Theta = x / \sqrt{x^2+y^2}$ has been used to simplify the expressions and to avoid divide-by-zero.
!latex \item Only ``poloidal tangential'' variations will be allowed to find the extremum of $F$, which are described by
!latex       \be 
!latex           \delta R_i(\t,\z) & \equiv & R_{i,\t}(\t,\z) \, \delta u_i(\t,\z),\\
!latex           \delta Z_i(\t,\z) & \equiv & Z_{i,\t}(\t,\z) \, \delta u_i(\t,\z),
!latex       \ee
!latex       from which it follows that the variation in each Fourier harmonic is
!latex       \be
!latex           \delta R_{j,i} = \ds \ooint R_{i,\t}(\t,\z) \, \delta u_i(\t,\z) \, \cos(m_j\t-n_j\z), \\
!latex           \delta Z_{j,i} = \ds \ooint Z_{i,\t}(\t,\z) \, \delta u_i(\t,\z) \, \sin(m_j\t-n_j\z),
!latex       \ee
!latex       and
!latex       \be 
!latex           \delta R_{i,\t}(\t,\z) & \equiv & R_{i,\t\t}(\t,\z) \, \delta u_i(\t,\z) + R_{i,\t}(\t,\z) \, \delta u_{i,\t} (\t,\z) \\
!latex           \delta Z_{i,\t}(\t,\z) & \equiv & Z_{i,\t\t}(\t,\z) \, \delta u_i(\t,\z) + Z_{i,\t}(\t,\z) \, \delta u_{i,\t} (\t,\z)
!latex       \ee
!latex \item The variation in $F$ is 
!latex       \be \delta F 
!latex       & = & \sum_{i=1}^{N-1} \alpha_i \;\;\;\ooint \!\! \left( \frac{-2\Theta_{i,\t\t}}{\Theta_{i,\t}^2} \right) \delta u_i \nonumber \\
!latex       & + & \sum_{i=1}^{N-1} \beta_i  \;\;\;\ooint \!\! \left( R_{i,\t} X_i + Z_{i,\t} Y_i \right) \delta u_i \nonumber \\
!latex       & + & \sum_{i=1}^{N-1} \gamma_i  \;\;\; \int \! d\z \;   \left( Z_{i}(0,\z)-Z_{i,0} \right) Z_{i,\t} \; \delta u_i \nonumber \\
!latex       & + & \red{\sum_{i=1}^{N-1} \delta_i \;\;\;\ooint \!\! \left( \frac{\Delta R_{i  } R_{i,\t} + \Delta Z_{i  } Z_{i,\t}}{L_{i  }} \right) 
!latex             \delta u_i}
!latex       \nonumber \\
!latex       & - & \red{\sum_{i=1}^{N-1} \delta_{i+1}   \ooint \!\! \left( \frac{\Delta R_{i+1} R_{i,\t} + \Delta Z_{i+1} Z_{i,\t}}{L_{i+1}} \right) 
!latex             \delta u_i}
!latex       \label{eq:firstvariation}
!latex       \ee
!latex       where, for the stellarator symmetric case,
!latex       \be
!latex           X_i & \equiv & \ds \sum\nolimits_j ( m_j^p - M_i ) \, R_{j,i} \cos(m_j\t-n_j\z), \\
!latex           Y_i & \equiv & \ds \sum\nolimits_j ( m_j^p - M_i ) \, Z_{j,i} \sin(m_j\t-n_j\z),
!latex       \ee
!latex       and 
!latex       \be \red{\Delta R_{i}} & \equiv & \red{R_{i}(\t,\z)-R_{i-1}(\t,\z)},\\
!latex           \red{\Delta Z_{i}} & \equiv & \red{Z_{i}(\t,\z)-Z_{i-1}(\t,\z)},
!latex       \ee
!latex \item The spectral constraints derived from \Eqn{firstvariation} are
!latex       \be I_i(\t,\z) & \equiv & - 2 \alpha_i \frac{\Theta_{i,\t\t}}{\Theta_{i,\t}^{2}}
!latex               + \beta_i \left( R_{i,\t} X_i + Z_{i,\t} Y_i \right)
!latex               + \gamma_i \left( Z_{i}(0,\z) - Z_{i,0} \right) Z_{i,\t}(0,\z) \nonumber \\
!latex             & + & \red{\delta_{i  } \frac{\Delta R_{i  } R_{i,\t} + \Delta Z_{i  } Z_{i,\t}}{L_i}}
!latex               -   \red{\delta_{i+1} \frac{\Delta R_{i+1} R_{i,\t} + \Delta Z_{i+1} Z_{i,\t}}{L_{i+1}}}
!latex       \label{eq:spectralconstraints}
!latex       \ee
!l!tex \item Note that choosing $p=2$ gives $X=-R_{\t\t}$ and $Y=-Z_{\t\t}$, and the spectrally condensed angle constraint, $R_\t X + Z_\t Y=0$, 
!l!tex       becomes $\partial_\t (R_\t^2+Z_\t^2)=0$, 
!l!tex       which defines the equal arc length angle.
!latex \item The poloidal-angle origin term, namely $\gamma_i \left( Z_{i}(0,\z) - Z_{i,0} \right) Z_{i,\t}(0,\z)$
!latex       is only used to constrain the $m_j=0$ harmonics.
!latex \item The construction of the angle functional was influenced by the following considerations:
!latex (i)   The minimal spectral width constraint is very desirable as it reduces the required Fourier resolution, 
!latex       but it does not constrain the $m=0$ harmonics
!latex       and the minimizing spectral-width poloidal-angle may not be consistent with the poloidal angle used on adjacent interfaces.
!latex (ii)  The regularization of the vector potential and the coordinate interpolation near the coordinate origin (see elsewhere)
!latex       assumes that the poloidal angle is the polar angle.
!latex (iii) The user will provide the Fourier harmonics of the boundary, and thus the user will implicitly define the poloidal angle used on the boundary.
!latex (iv)  Minimizing the length term will ensure that the poloidal angle used on each interface
!latex       is smoothly connected to the poloidal angle used on adjacent interfaces.
!latex \item A suitable choice of the weight factors, $\alpha_i$, $\beta_i$, $\gamma_i$ and $\delta_i$, will ensure
!latex       that the polar constraint dominates for the innermost surfaces and that this constraint rapidly becomes insignificant away from the origin;
!latex       that the minimal spectral constraint dominates in the ``middle'';
!latex       and that the minimizing length constraint will be significant near the origin and dominant near the edge,
!latex       so that the minimizing spectral width angle will be continuously connected to the polar angle on the innermost surfaces
!latex       and the user-implied angle at the plasma boundary.
!latex       The length constraint should not be insignificant where the spectral constraint is dominant (so that the $m=0$ harmonics are constrained).
!latex \item The polar constraint does not need normalization.
!latex       The spectral width constraint has already been normalized.
!latex       The length constraint is not yet normalized, but perhaps it should be.
!latex \item The spectral constraints given in \Eqn{spectralconstraints} need to be differentiated
!latex       with respect to the interface Fourier harmonics, $R_{j,i}$ and $Z_{j,i}$.
!latex       The first and second terms lead to a block diagonal hessian, and the length term leads to a block tri-diagonal hessian.
!latex \item Including the poloidal-angle origin constraint means that the polar angle constraint can probably be ignored, i.e. $\alpha_i=0$.
!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine lforce( lvol, iocons, ideriv, Ntz, dAt, dAz, XX, YY, length, DDl, MMl, iflag )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half, one, two
  
  use fileunits, only : ounit
  
  use inputlist, only : Wlforce, Igeometry, Nvol, Lrad, gamma, pscale, adiabatic, Lcheck
  
  use cputiming, only : Tlforce
  
  use allglobal, only : ncpu, myid, cpus, &
                        Lcoordinatesingularity, Mvol, &
                        iRbc, iZbs, iRbs, iZbc, &
                        YESstellsym, NOTstellsym, &
                        mn, im, in, regumm, &
                        ijreal, ijimag, jireal, jiimag, &
                        efmn, ofmn, cfmn, sfmn, evmn, odmn, comn, simn, &
                        Nt, Nz, &
                        Ate, Aze, Ato, Azo, &
                        TT, &
                        sg, guvij, iRij, iZij, dRij, dZij, tRij, tZij, &
                        mmpp, &
                        Bemn, Bomn, Iomn, Iemn, Somn, Semn, &
                        Pomn, Pemn, &
                        vvolume, & 
                        build_vector_potential
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in)  :: lvol, iocons, ideriv, Ntz, iflag
  REAL                 :: dAt(1:Ntz), dAz(1:Ntz), XX(1:Ntz), YY(1:Ntz), dRR(1:Ntz,-1:1), dZZ(1:Ntz,-1:1), DDl, MMl

  REAL                 :: IIl(1:Ntz), length(1:Ntz), dLL(1:Ntz)
  
  INTEGER              :: Lcurvature, ii, jj, kk, ll, ifail, ivol, lnn!, oicons
  REAL                 :: dBB(1:Ntz), lss, mfactor
  
  REAL                 :: dAs(1:Ntz)!, dRdt(-1:1,0:1), dZdt(-1:1,0:1)
  REAL                 :: lgvuij(1:Ntz,1:3,1:3) ! local workspace; 13 Sep 13;
  
  BEGIN(lforce)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  FATAL( lforce, lvol.lt.1 .or. lvol.gt.Mvol, illegal lvol )
  FATAL( lforce, lvol.eq.1 .and. iocons.eq.0, illegal combination )
  FATAL( lforce, lvol.eq.Mvol .and. iocons.eq.1, illegal combination )
  FATAL( lforce, iflag.lt.0 .or. iflag.gt.1, illegal iflag )
#endif
 
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  dAt(1:Ntz) = zero ! initialize intent out; 01 Jul 14;
  dAz(1:Ntz) = zero ! initialize intent out; 01 Jul 14;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  lss = two * iocons - one ! recall that iocons is effective local radial coordinate; 24 Apr 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  Lcurvature = 1

  WCALL( lforce, coords, ( lvol, lss, Lcurvature, Ntz, mn ) ) ! get coordinates and derivatives wrt Rj, Zj, at specific radial location;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! if( lvol.eq.Mvol .and. ideriv.eq.0 ) then
!  ii = 1
!  write(ounit,'("lforce : ", 10x ," : sum(Ate(",i3,",",i2,",",i2,")%s) =",99es23.15)') lvol, ideriv, ii, sum(Ate(lvol,ideriv,ii)%s(0:Lrad(lvol)))
! endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! compute B^2 on interface;
  
  call build_vector_potential(lvol, iocons, ideriv, 1)
  call invfft( mn, im, in, efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, dAt(1:Ntz), dAz(1:Ntz) ) ! map to real space;

  dBB(1:Ntz) = half * (         dAz(1:Ntz   )*dAz(1:Ntz   )*guvij(1:Ntz,2,2,0) &
                        - two * dAz(1:Ntz   )*dAt(1:Ntz   )*guvij(1:Ntz,2,3,0) &
                        +       dAt(1:Ntz   )*dAt(1:Ntz   )*guvij(1:Ntz,3,3,0) ) / sg(1:Ntz,0)**2
   
  ijreal(1:Ntz) = adiabatic(lvol) * pscale / vvolume(lvol)**gamma + dBB(1:Ntz) ! p + B^2/2; 13 Sep 13;

#ifdef DEBUG
	if( Wlforce ) then
		write(ounit, 8375) lvol, iocons, ideriv, dAz(1:Ntz), dAt(1:Ntz)
		write(ounit, 8376) lvol, iocons, ideriv, guvij(1:Ntz,2,2,0), guvij(1:Ntz,2,3,0), guvij(1:Ntz,3,3,0), sg(1:Ntz,0)
	 
8375 format("lforce : lvol=",i7,", iocons=", i7, ", ideriv=", i7 ,"; dAz=",f10.6,", dAt=", f10.6)
8376 format("lforce : lvol=",i7,", iocons=", i7, ", ideriv=", i7 ,"; g22=",f10.6,", g23=", f10.6,", g33=", f10.6,", sg=", f10.6)
	endif
#endif


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
  if( iflag .eq. 1 ) goto 9999 ! iflag = 1 indicates the derivatives of the force are to be calculated; derivatives of magnetic field calculated above;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! compute spectral constraints;
  
  select case( Igeometry )
   
  case( 1:2 ) ; dLL(1:Ntz) = zero ! placeholder; 08 Feb 16;
   ;          ; IIl(1:Ntz) = zero ! placeholder; 11 Aug 14;   
  case(   3 )
   
   do ivol = 0, 1
    
    call invfft( mn, im(1:mn), in(1:mn),            iRbc(1:mn ,lvol-1+ivol),              iRbs(1:mn ,lvol-1+ivol), &
                                                    iZbc(1:mn ,lvol-1+ivol),              iZbs(1:mn ,lvol-1+ivol), & 
                                         Nt, Nz,    iRij(1:Ntz,lvol-1+ivol),              iZij(1:Ntz,lvol-1+ivol)   )

    call invfft( mn, im(1:mn), in(1:mn), im(1:mn) * iRbs(1:mn ,lvol-1+ivol), - im(1:mn) * iRbc(1:mn ,lvol-1+ivol), &
                                         im(1:mn) * iZbs(1:mn ,lvol-1+ivol), - im(1:mn) * iZbc(1:mn ,lvol-1+ivol), &
                                         Nt, Nz,    tRij(1:Ntz,lvol-1+ivol),              tZij(1:Ntz,lvol-1+ivol) )
   enddo ! end of do ivol = 0, 1 ; 18 Jul 14;
   
   dRij(1:Ntz,lvol) = iRij(1:Ntz,lvol) - iRij(1:Ntz,lvol-1)
   dZij(1:Ntz,lvol) = iZij(1:Ntz,lvol) - iZij(1:Ntz,lvol-1)
   
   length(1:Ntz) = sqrt( dRij(1:Ntz,lvol)**2 + dZij(1:Ntz,lvol)**2 )
   
   dLL(1:Ntz) = ( dRij(1:Ntz,lvol) * tRij(1:Ntz,lvol-1+iocons) + dZij(1:Ntz,lvol) * tZij(1:Ntz,lvol-1+iocons) ) / length(1:Ntz)
   
   if( iocons.eq.1 ) then ! include spectral condensation constraints; local to interface, i.e. no tri-diagonal structure;
    ;                      ; efmn(1:mn) = ( mmpp(1:mn)            ) * iRbc(1:mn,lvol)
    ;                      ; sfmn(1:mn) = ( mmpp(1:mn)            ) * iZbs(1:mn,lvol)
    if( NOTstellsym ) then ; ofmn(1:mn) = ( mmpp(1:mn)            ) * iRbs(1:mn,lvol)
     ;                     ; cfmn(1:mn) = ( mmpp(1:mn)            ) * iZbc(1:mn,lvol)
    else                   ; ofmn(1:mn) = zero
     ;                     ; cfmn(1:mn) = zero
    endif ! end of if( NOTstellsym ) ; 20 Feb 13;
    
    call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), &
                 Nt, Nz, XX(1:Ntz), YY(1:Ntz) )
    
    if( YESstellsym ) then ; DDl = sum(              ( iRbc(1:mn,lvol)**2 + iZbs(1:mn,lvol)**2                                           ) )
     ;                     ; MMl = sum( mmpp(1:mn) * ( iRbc(1:mn,lvol)**2 + iZbs(1:mn,lvol)**2                                           ) ) / DDl
    else                   ; DDl = sum(              ( iRbc(1:mn,lvol)**2 + iZbs(1:mn,lvol)**2 + iRbs(1:mn,lvol)**2 + iZbc(1:mn,lvol)**2 ) )
     ;                     ; MMl = sum( mmpp(1:mn) * ( iRbc(1:mn,lvol)**2 + iZbs(1:mn,lvol)**2 + iRbs(1:mn,lvol)**2 + iZbc(1:mn,lvol)**2 ) ) / DDl
    endif
    
    IIl(1:Ntz) = tRij(1:Ntz,lvol) * ( XX(1:Ntz) - MMl * iRij(1:Ntz,lvol) ) &
               + tZij(1:Ntz,lvol) * ( YY(1:Ntz) - MMl * iZij(1:Ntz,lvol) ) 
    
   else ! matches if( iocons.eq.1 ) ; 11 Aug 14;

    IIl(1:Ntz) = zero ! placeholder; 11 Aug 14;

   endif ! end of if( iocons.eq.1 ) ; 20 Feb 13;
   
  end select ! end of select case( Igeometry ) ; 08 Feb 16;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!  ijreal(1:Ntz) contains the pressure + magnetic energy term;
  
#ifdef DEBUG
  FATAL( lforce, iocons.lt.0 .or. iocons.gt.2, error )
#endif
  
  ;ifail = 0
  ;call tfft( Nt, Nz, ijreal(1:Ntz), IIl(1:Ntz  ), & ! compute force-imbalance and spectral constraints;
              mn, im(1:mn), in(1:mn), Bemn(1:mn,lvol,iocons), Bomn(1:mn,lvol,iocons), Iemn(1:mn,lvol       ), Iomn(1:mn,lvol       ), ifail )
  
  if( Igeometry.ge.3 ) then ! add minimal length constraint; 18 Jul 14;
   
   ifail = 0 ; ijimag(1:Ntz) = zero

   call tfft( Nt, Nz, dLL(1:Ntz), ijimag(1:Ntz), &
              mn, im(1:mn), in(1:mn), Semn(1:mn,lvol,iocons), Somn(1:mn,lvol,iocons), Pemn(1:mn,lvol,iocons), Pomn(1:mn,lvol,iocons), ifail )

#ifdef DEBUG
   if( Wlforce ) then
    write(ounit,'("lforce : ", 10x ," : lvol=",i3," ; iocons="i2" ; Somn="999es13.5)') lvol, iocons, Somn(1:mn,lvol,iocons)
    write(ounit,'("lforce : ", 10x ," : lvol=",i3," ; iocons="i2" ; Semn="999es13.5)') lvol, iocons, Semn(1:mn,lvol,iocons)
   endif
#endif
   
  endif ! end of if( Igeometry.eq.3 ) ; 01 Jul 14;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(lforce)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
 end subroutine lforce
 
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
