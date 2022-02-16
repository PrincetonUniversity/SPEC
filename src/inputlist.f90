!> \file
!> \brief Input namelists

!> \brief Input namelists
!> \addtogroup grp_global
!> @{
module inputlist

  implicit none

  ! The following parameters set the maximum allowed resolution:
  INTEGER, parameter :: MNvol     = 256 !< The maximum value of \c Nvol is \c MNvol=256.
  INTEGER, parameter :: MMpol     =  64 !< The maximum value of \c Mpol is \c MNpol=64.
  INTEGER, parameter :: MNtor     =  64 !< The maximum value of \c Ntor is \c MNtor=64.

!> \addtogroup grp_global_physicslist physicslist
!> \brief The namelist \c physicslist controls the geometry, profiles, and numerical resolution.
!> @{
  ! note that all variables in namelist need to be broadcasted in readin
  INTEGER      :: Igeometry                  =  3        !< selects Cartesian, cylindrical or toroidal geometry;
                                                         !< <ul>
                                                         !< <li> \c Igeometry=1 : Cartesian; geometry determined by \f$R\f$; </li>
                                                         !< <li> \c Igeometry=2 : cylindrical; geometry determined by \f$R\f$; </li>
                                                         !< <li> \c Igeometry=3 : toroidal; geometry determined by \f$R\f$ *and* \f$Z\f$; </li>
                                                         !< </ul>
  INTEGER      :: Istellsym                  =  1        !< stellarator symmetry is enforced if \c Istellsym==1
  INTEGER      :: Lfreebound                 =  0        !< compute vacuum field surrounding plasma
  REAL         :: phiedge                    =  1.0      !< total enclosed toroidal magnetic flux;
  REAL         :: curtor                     =  0.0      !< total enclosed (toroidal) plasma current;
  REAL         :: curpol                     =  0.0      !< total enclosed (poloidal) linking current;
  REAL         :: gamma                      =  0.0      !< adiabatic index; cannot set \f$|\gamma| = 1\f$
  INTEGER      :: Nfp                        =  1        !< field periodicity
                                                         !< <ul>
                                                         !< <li> all Fourier representations are of the form \f$\cos(m\theta-n N \zeta)\f$, \f$\sin(m\theta-n N \zeta)\f$, where \f$N\equiv\f$\c Nfp </li>
                                                         !< <li> constraint: \c Nfp >= 1 </li>
                                                         !< </ul>
  INTEGER      :: Nvol                       =  1        !< number of volumes
                                                         !< <ul>
                                                         !< <li> each volume \f${\cal V}_l\f$ is bounded by the \f${\cal I}_{l-1}\f$ and \f${\cal I}_{l}\f$ interfaces </li>
                                                         !< <li> note that in cylindrical or toroidal geometry, \f${\cal I}_{0}\f$ is the degenerate coordinate axis </li>
                                                         !< <li> constraint: \c Nvol<=MNvol </li>
                                                         !< </ul>
  INTEGER      :: Mpol                       =  0        !< number of poloidal Fourier harmonics
                                                         !< <ul>
                                                         !< <li> all Fourier representations of doubly-periodic functions are of the form
                                                         !< \f{eqnarray}{ f(\theta,\zeta) & = & \sum_{n=0}^{\texttt{Ntor}} f_{0,n}\cos(-n \, \texttt{Nfp} \, \zeta)
                                                         !<                                   + \sum_{m=1}^{\texttt{Mpol}}\sum_{n=\texttt{-Ntor}}^{\texttt{Ntor}} f_{m,n}\cos(m\theta-n \, \texttt{Nfp} \, \zeta),
                                                         !< \f}
                                                         !< Internally these "double" summations are written as a "single" summation,
                                                         !< e.g. \f$f(\theta,\zeta) = \sum_j f_j \cos(m_j\theta-n_j\zeta)\f$. </li>
                                                         !< </ul>
  INTEGER      :: Ntor                       =  0        !< number of toroidal Fourier harmonics
                                                         !< <ul>
                                                         !< <li> all Fourier representations of doubly-periodic functions are of the form
                                                         !< \f{eqnarray}{ f(\theta,\zeta) & = & \sum_{n=0}^{\texttt{Ntor}} f_{0,n}\cos(-n \, \texttt{Nfp} \, \zeta)
                                                         !<                                   + \sum_{m=1}^{\texttt{Mpol}}\sum_{n=\texttt{-Ntor}}^{\texttt{Ntor}} f_{m,n}\cos(m\theta-n \, \texttt{Nfp} \, \zeta),
                                                         !< \f}
                                                         !< Internally these "double" summations are written as a "single" summation,
                                                         !< e.g. \f$f(\theta,\zeta) = \sum_j f_j \cos(m_j\theta-n_j\zeta)\f$. </li>
                                                         !< </ul>
  INTEGER      :: Lrad(1:MNvol+1)            =  4        !< Chebyshev resolution in each volume
                                                         !< <ul>
                                                         !< <li> constraint : \c Lrad(1:Mvol) >= 2 </li>
                                                         !< </ul>
  INTEGER      :: Lconstraint                = -1        !< selects constraints; primarily used in ma02aa() and mp00ac().
                                                         !< <ul>
                                                         !< <li> if \c Lconstraint==-1, then in the plasma regions \f$\Delta\psi_t\f$, \f$\mu\f$ and \f$\Delta \psi_p\f$ are *not* varied
                                                         !<      and in the vacuum region (only for free-boundary) \f$\Delta\psi_t\f$ and \f$\Delta \psi_p\f$ are *not* varied, and \f$\mu = 0\f$. </li>
                                                         !< <li> if \c Lconstraint==0, then in the plasma regions \f$\Delta\psi_t\f$, \f$\mu\f$ and \f$\Delta \psi_p\f$ are *not* varied
                                                         !<      and in the vacuum region (only for free-boundary) \f$\Delta\psi_t\f$ and \f$\Delta \psi_p\f$ are varied to match the
                                                         !<      prescribed plasma current, \c curtor, and the "linking" current, \c curpol, and \f$\mu = 0\f$ </li>
                                                         !< <li> if \c Lconstraint==1, then in the plasma regions \f$\mu\f$ and \f$\Delta\psi_p\f$ are adjusted
                                                         !<      in order to satisfy the inner and outer interface transform constraints
                                                         !<      (except in the simple torus, where the enclosed poloidal flux is irrelevant,
                                                         !<      and only \f$\mu\f$ is varied to satisfy the outer interface transform constraint);
                                                         !<      and in the vacuum region \f$\Delta\psi_t\f$ and \f$\Delta \psi_p\f$ are varied to match the transform constraint on the boundary
                                                         !<      and to obtain the prescribed linking current, \c curpol, and \f$\mu = 0\f$. </li>
                                                         !< <li> \todo if \c Lconstraint==2, under reconstruction.
                                                         !<
                                                         !< </li>
                                                         !< <li> if \c Lconstraint.eq.3 , then the \f$\mu\f$ and \f$\psi_p\f$ variables are adjusted
                                                         !<      in order to satisfy the volume and surface toroidal current computed with lbpol()
                                                         !<      (excepted in the inner most volume, where the volume current is irrelevant).
                                                         !<      Not implemented yet in free boundary.</li>
                                                         !< </ul>
  REAL         ::     tflux(1:MNvol+1)       =  0.0      !< toroidal flux, \f$\psi_t\f$, enclosed by each interface
                                                         !< <ul>
                                                         !< <li> For each of the plasma volumes, this is a constraint: \c tflux is *not* varied </li>
                                                         !< <li> For the vacuum region (only if \c Lfreebound==1), \c tflux  may be allowed to vary to match constraints </li>
                                                         !< <li> Note that \c tflux  will be normalized so that \c tflux(Nvol) = 1.0,
                                                         !<      so that \c tflux  is arbitrary up to a scale factor </li>
                                                         !< <li> \sa phiedge </li>
                                                         !< </ul>
  REAL         ::     pflux(1:MNvol+1)       =  0.0      !< poloidal flux, \f$\psi_p\f$, enclosed by each interface
  REAL         ::  helicity(1:MNvol)         =  0.0      !< helicity, \f${\cal K}\f$, in each volume, \f${\cal V}_i\f$
                                                         !< <ul>
                                                         !< <li> on exit, \c helicity  is set to the computed values of \f${\cal K} \equiv \int {\bf A}\cdot{\bf B}\;dv\f$ </li>
                                                         !< </ul>
  REAL         :: pscale                     =  0.0      !< pressure scale factor
                                                         !< <ul>
                                                         !< <li> the initial pressure profile is given by \c pscale  \f$*\f$ \c pressure </li>
                                                         !< </ul>
  REAL         ::  pressure(1:MNvol+1)       =  0.0      !< pressure in each volume
                                                         !< <ul>
                                                         !< <li> The pressure is *not* held constant, but \f$p_l V_l^\gamma = P_l\f$ *is* held constant,
                                                         !<      where \f$P_l\f$ is determined by the initial pressures and the initial volumes, \f$V_l\f$. </li>
                                                         !< <li> Note that if \c gamma==0.0, then \f$p_l \equiv P_l\f$. </li>
                                                         !< <li> On output, the pressure is given by \f$p_l = P_l/V_l^\gamma\f$, where \f$V_l\f$ is the final volume. </li>
                                                         !< <li> \c pressure is only used in calculation of interface force-balance. </li>
                                                         !< </ul>
  INTEGER      :: Ladiabatic                 =  0        !< logical flag
                                                         !< <ul>
                                                         !< <li> If \c Ladiabatic==0, the adiabatic constants are determined by the initial pressure and volume. </li>
                                                         !< <li> If \c Ladiabatic==1, the adiabatic constants are determined by the given input \c adiabatic. </li>
                                                         !< </ul>
  REAL         :: adiabatic(1:MNvol+1)       =  0.0      !< adiabatic constants in each volume
                                                         !< <ul>
                                                         !< <li> The pressure is *not* held constant, but \f$p_l V_l^\gamma = P_l \equiv\f$\c adiabatic is constant. </li>
                                                         !< <li> Note that if \c gamma==0.0, then \c pressure==adiabatic. </li>
                                                         !< <li> \c pressure is only used in calculation of interface force-balance. </li>
                                                         !< </ul>
  REAL         ::        mu(1:MNvol+1)       =  0.0      !< helicity-multiplier, \f$\mu\f$, in each volume
  REAL         :: Ivolume(1:MNvol+1)         =  0.0      !< Toroidal current constraint normalized by \f$\mu_0\f$ (\f$I_{volume} = \mu_0\cdot [A]\f$), in each volume.
                                                         !< This is a cumulative quantity: \f$I_{\mathcal{V},i} = \int_0^{\psi_{t,i}} \mathbf{J}\cdot\mathbf{dS}\f$.
                                                         !< Physically, it represents the sum of all non-pressure driven currents.
  REAL         :: Isurf(1:MNvol)             =  0.0      !< Toroidal current normalized by \f$\mu_0\f$ at each interface (cumulative). This is the sum of all pressure driven currents.
  INTEGER      ::        pl(0:MNvol)         =  0        !< "inside" interface rotational-transform is \f$\mbox{$\,\iota\!\!$-} = (p_l+\gamma p_r)/(q_l+\gamma q_r)\f$,
                                                         !< where \f$\gamma\f$ is the golden mean, \f$\gamma = (1 + \sqrt 5 ) / 2 \f$.
                                                         !<
                                                         !< If both \f$q_l = 0\f$ *and* \f$q_r = 0\f$, then the (inside) interface rotational-transform is defined by \c iota .
  INTEGER      ::        ql(0:MNvol)         =  0        !< "inside" interface rotational-transform is \f$\mbox{$\,\iota\!\!$-} = (p_l+\gamma p_r)/(q_l+\gamma q_r)\f$,
                                                         !< where \f$\gamma\f$ is the golden mean, \f$\gamma = (1 + \sqrt 5 ) / 2 \f$.
                                                         !<
                                                         !< If both \f$q_l = 0\f$ *and* \f$q_r = 0\f$, then the (inside) interface rotational-transform is defined by \c iota .
  INTEGER      ::        pr(0:MNvol)         =  0        !< "inside" interface rotational-transform is \f$\mbox{$\,\iota\!\!$-} = (p_l+\gamma p_r)/(q_l+\gamma q_r)\f$,
                                                         !< where \f$\gamma\f$ is the golden mean, \f$\gamma = (1 + \sqrt 5 ) / 2 \f$.
                                                         !<
                                                         !< If both \f$q_l = 0\f$ *and* \f$q_r = 0\f$, then the (inside) interface rotational-transform is defined by \c iota .
  INTEGER      ::        qr(0:MNvol)         =  0        !< "inside" interface rotational-transform is \f$\mbox{$\,\iota\!\!$-} = (p_l+\gamma p_r)/(q_l+\gamma q_r)\f$,
                                                         !< where \f$\gamma\f$ is the golden mean, \f$\gamma = (1 + \sqrt 5 ) / 2 \f$.
                                                         !<
                                                         !< If both \f$q_l = 0\f$ *and* \f$q_r = 0\f$, then the (inside) interface rotational-transform is defined by \c iota .
  REAL         ::      iota(0:MNvol)         =  0.0      !< rotational-transform, \f$\mbox{$\,\iota\!\!$-}\f$, on inner side of each interface
                                                         !< <ul>
                                                         !< <li> only relevant if illogical input for \c ql and \c qr are provided </li>
                                                         !< </ul>
  INTEGER      ::        lp(0:MNvol)         =  0        !< "outer" interface rotational-transform is \f$\mbox{$\,\iota\!\!$-} = (p_l+\gamma p_r)/(q_l+\gamma q_r)\f$,
                                                         !<  where \f$\gamma\f$ is the golden mean, \f$\gamma = (1 + \sqrt 5 ) / 2\f$.
                                                         !<
                                                         !< If both \f$q_l = 0\f$ *and* \f$q_r = 0\f$, then the (outer) interface rotational-transform is defined by \c oita .
  INTEGER      ::        lq(0:MNvol)         =  0        !< "outer" interface rotational-transform is \f$\mbox{$\,\iota\!\!$-} = (p_l+\gamma p_r)/(q_l+\gamma q_r)\f$,
                                                         !<  where \f$\gamma\f$ is the golden mean, \f$\gamma = (1 + \sqrt 5 ) / 2\f$.
                                                         !<
                                                         !< If both \f$q_l = 0\f$ *and* \f$q_r = 0\f$, then the (outer) interface rotational-transform is defined by \c oita .
  INTEGER      ::        rp(0:MNvol)         =  0        !< "outer" interface rotational-transform is \f$\mbox{$\,\iota\!\!$-} = (p_l+\gamma p_r)/(q_l+\gamma q_r)\f$,
                                                         !<  where \f$\gamma\f$ is the golden mean, \f$\gamma = (1 + \sqrt 5 ) / 2\f$.
                                                         !<
                                                         !< If both \f$q_l = 0\f$ *and* \f$q_r = 0\f$, then the (outer) interface rotational-transform is defined by \c oita .
  INTEGER      ::        rq(0:MNvol)         =  0        !< "outer" interface rotational-transform is \f$\mbox{$\,\iota\!\!$-} = (p_l+\gamma p_r)/(q_l+\gamma q_r)\f$,
                                                         !<  where \f$\gamma\f$ is the golden mean, \f$\gamma = (1 + \sqrt 5 ) / 2\f$.
                                                         !<
                                                         !< If both \f$q_l = 0\f$ *and* \f$q_r = 0\f$, then the (outer) interface rotational-transform is defined by \c oita .
  REAL         ::      oita(0:MNvol)         =  0.0      !< rotational-transform, \f$\mbox{$\,\iota\!\!$-}\f$, on outer side of each interface
                                                         !< <ul>
                                                         !< <li> only relevant if illogical input for \c ql and \c qr are provided </li>
                                                         !< </ul>
  REAL         :: mupftol                    =  1.0e-14  !< accuracy to which \f$\mu\f$ and \f$\Delta\psi_p\f$ are required
                                                         !< <ul>
                                                         !< <li> only relevant if constraints on transform, enclosed currents etc. are to be satisfied iteratively, see \c Lconstraint </li>
                                                         !< </ul>
  INTEGER      :: mupfits                    =  8        !< an upper limit on the transform/helicity constraint iterations;
                                                         !< <ul>
                                                         !< <li> only relevant if constraints on transform, enclosed currents etc. are to be satisfied iteratively, see \c Lconstraint </li>
                                                         !< <li> constraint: \c mupfits > 0 </li>
                                                         !< </ul>
  REAL         :: rpol                       =  1.0      !< poloidal extent of slab (effective radius)
                                                         !< <ul>
                                                         !< <li> only relevant if \c Igeometry==1 </li>
                                                         !< <li> poloidal size is \f$L = 2\pi*\f$\c rpol </li>
                                                         !< </ul>
  REAL         :: rtor                       =  1.0      !< toroidal extent of slab (effective radius)
                                                         !< <ul>
                                                         !< <li> only relevant if \c Igeometry==1 </li>
                                                         !< <li> toroidal size is \f$L = 2\pi*\f$\c rtor </li>
                                                         !< </ul>
  INTEGER      :: Lreflect                   =  0        !< =1 reflect the upper and lower bound in slab, =0 do not reflect

  REAL         :: Rac(     0:MNtor        )  =  0.0      !<     stellarator symmetric coordinate axis;
  REAL         :: Zas(     0:MNtor        )  =  0.0      !<     stellarator symmetric coordinate axis;
  REAL         :: Ras(     0:MNtor        )  =  0.0      !< non-stellarator symmetric coordinate axis;
  REAL         :: Zac(     0:MNtor        )  =  0.0      !< non-stellarator symmetric coordinate axis;

  REAL         :: Rbc(-MNtor:MNtor,-MMpol:MMpol)  =  0.0 !<     stellarator symmetric boundary components;
  REAL         :: Zbs(-MNtor:MNtor,-MMpol:MMpol)  =  0.0 !<     stellarator symmetric boundary components;
  REAL         :: Rbs(-MNtor:MNtor,-MMpol:MMpol)  =  0.0 !< non-stellarator symmetric boundary components;
  REAL         :: Zbc(-MNtor:MNtor,-MMpol:MMpol)  =  0.0 !< non-stellarator symmetric boundary components;

  REAL         :: Rwc(-MNtor:MNtor,-MMpol:MMpol)  =  0.0 !<     stellarator symmetric boundary components of wall;
  REAL         :: Zws(-MNtor:MNtor,-MMpol:MMpol)  =  0.0 !<     stellarator symmetric boundary components of wall;
  REAL         :: Rws(-MNtor:MNtor,-MMpol:MMpol)  =  0.0 !< non-stellarator symmetric boundary components of wall;
  REAL         :: Zwc(-MNtor:MNtor,-MMpol:MMpol)  =  0.0 !< non-stellarator symmetric boundary components of wall;

  REAL         :: Vns(-MNtor:MNtor,-MMpol:MMpol)  =  0.0 !<     stellarator symmetric normal field at boundary; vacuum component;
  REAL         :: Bns(-MNtor:MNtor,-MMpol:MMpol)  =  0.0 !<     stellarator symmetric normal field at boundary; plasma component;
  REAL         :: Vnc(-MNtor:MNtor,-MMpol:MMpol)  =  0.0 !< non-stellarator symmetric normal field at boundary; vacuum component;
  REAL         :: Bnc(-MNtor:MNtor,-MMpol:MMpol)  =  0.0 !< non-stellarator symmetric normal field at boundary; plasma component;
!> @}

!> \addtogroup grp_global_numerics numericlist
!> \brief The namelist \c numericlist controls internal resolution parameters that the user rarely needs to consider.
!> @{
  INTEGER      :: Linitialize =  0   !< Used to initialize geometry using a regularization / extrapolation method
                                     !< <ul>
                                     !< <li> if \c Linitialize = \f$-I\f$ , where \f$I\f$ is a positive integer,
                                     !<      the geometry of the \f$i=1,N_V-I\f$ surfaces constructed by an extrapolation </li>
                                     !< <li> if \c Linitialize = 0, the geometry of the interior surfaces is provided after the namelists in the input file </li>
                                     !< <li> if \c Linitialize = 1, the interior surfaces will be intialized as \f$R_{l,m,n} = R_{N,m,n} \psi_{t,l}^{m/2}\f$,
                                     !<      where \f$R_{N,m,n}\f$ is the plasma boundary and \f$\psi_{t,l}\f$ is the given toroidal flux enclosed by the
                                     !<      \f$l\f$-th interface, normalized to the total enclosed toroidal flux;
                                     !<      a similar extrapolation is used for \f$Z_{l,m,n}\f$ </li>
                                     !< <li> Note that the Fourier harmonics of the boundary is *always* given by the \c Rbc and \c Zbs
                                     !<      given in \c physicslist. </li>
                                     !< <li> if \c Linitialize = 2, the interior surfaces *and the plasma boundary* will be intialized
                                     !<       as \f$R_{l,m,n} = R_{W,m,n} \psi_{t,l}^{m/2}\f$, where \f$R_{W,m,n}\f$ is the computational boundary
                                     !<       and \f$\psi_{t,l}\f$ is the given toroidal flux enclosed by the \f$l\f$-th interface, normalized to the total enclosed toroidal flux;
                                     !<       a similar extrapolation is used for \f$Z_{l,m,n}\f$ </li>
                                     !< <li> Note that, for free-boundary calculations, the Fourier harmonics of the computational boundary
                                     !<      are *always* given by the \c Rwc and \c Zws given in \c physicslist. </li>
                                     !< <li> if \c Linitialize = 1, 2, it is not required to provide the geometry of the interfaces after the namelists </li>
                                     !< </ul>
  INTEGER      :: LautoinitBn =  1   !< Used to initialize \f$B_{ns}\f$ using an initial fixed-boundary calculation
                                     !< <ul>
                                     !< <li> only relevant if \c Lfreebound = 1 </li>
                                     !< <li> user-supplied \c Bns will only be considered if \c LautoinitBn = 0 </li>
                                     !< </ul>
  INTEGER      :: Lzerovac    =  0   !< Used to adjust vacuum field to cancel plasma field on computational boundary
                                     !< <ul>
                                     !< <li> only relevant if \c Lfreebound = 1 </li>
                                     !< </ul>
  INTEGER      :: Ndiscrete   =  2   !< resolution of the real space grid on which fast Fourier transforms are performed is given by \c Ndiscrete*Mpol*4
                                     !< <ul>
                                     !< <li> constraint \c Ndiscrete>0 </li>
                                     !< </ul>
  INTEGER      :: Nquad       = -1   !< Resolution of the Gaussian quadrature
                                     !< <ul>
                                     !< <li> The resolution of the Gaussian quadrature, \f$\displaystyle \int \!\! f(s) ds = \sum_k \omega_k f(s_k)\f$,
                                     !<      in each volume is given by \c Iquad\f$_v\f$,  </li>
                                     !< <li> \c Iquad\f$_v\f$ is set in preset() </li>
                                     !< </ul>
  INTEGER      :: iMpol       = -4   !< Fourier resolution of straight-fieldline angle on interfaces
                                     !< <ul>
                                     !< <li> the rotational-transform on the interfaces is determined by a transformation to the straight-fieldline angle,
                                     !<      with poloidal resolution given by \c iMpol </li>
                                     !< <li> if \c iMpol<=0, then \c iMpol = Mpol - iMpol </li>
                                     !< </ul>
  INTEGER      :: iNtor       = -4   !< Fourier resolution of straight-fieldline angle on interfaces;
                                     !< <ul>
                                     !< <li> the rotational-transform on the interfaces is determined by a transformation to the straight-fieldline angle,
                                     !<      with toroidal resolution given by \c iNtor </li>
                                     !< <li> if \c iNtor<=0 then \c iNtor = Ntor - iNtor </li>
                                     !< <li> if \c Ntor==0, then the toroidal resolution of the angle transformation is set \c lNtor = 0 </li>
                                     !< </ul>
  INTEGER      :: Lsparse     =  0   !< controls method used to solve for rotational-transform on interfaces
                                     !< <ul>
                                     !< <li> if \c Lsparse = 0, the transformation to the straight-fieldline angle is computed in Fourier space
                                     !<      using a dense matrix solver,  \c F04AAF </li>
                                     !< <li> if \c Lsparse = 1, the transformation to the straight-fieldline angle is computed in real space
                                     !<      using a dense matrix solver,  \c F04ATF </li>
                                     !< <li> if \c Lsparse = 2, the transformation to the straight-fieldline angle is computed in real space
                                     !<      using a sparse matrix solver, \c F11DEF </li>
                                     !< <li> if \c Lsparse = 3, the different methods for constructing the straight-fieldline angle are compared </li>
                                     !< </ul>
  INTEGER      :: Lsvdiota    =  0   !< controls method used to solve for rotational-transform on interfaces;
                                     !< only relevant if \c Lsparse = 0
                                     !< <ul>
                                     !< <li> if \c Lsvdiota = 0, use standard linear solver to construct straight fieldline angle transformation </li>
                                     !< <li> if \c Lsvdiota = 1, use SVD method to compute rotational-transform </li>
                                     !< </ul>
  INTEGER      :: imethod     =  3   !< controls iterative solution to sparse matrix
                                     !< arising in real-space transformation to the straight-fieldline angle;
                                     !< only relevant if \c Lsparse.eq.2; \see tr00ab() for details
                                     !< <ul>
                                     !< <li> if \c imethod = 1, the method is \c RGMRES   </li>
                                     !< <li> if \c imethod = 2, the method is \c CGS      </li>
                                     !< <li> if \c imethod = 3, the method is \c BICGSTAB </li>
                                     !< </ul>
  INTEGER      :: iorder      =  2   !< controls real-space grid resolution for constructing the straight-fieldline angle;
                                     !< only relevant if \c Lsparse>0
                                     !<
                                     !< determines order of finite-difference approximation to the derivatives
                                     !< <ul>
                                     !< <li> if \c iorder = 2,  </li>
                                     !< <li> if \c iorder = 4,  </li>
                                     !< <li> if \c iorder = 6,  </li>
                                     !< </ul>
  INTEGER      :: iprecon     =  0   !< controls iterative solution to sparse matrix arising in real-space transformation
                                     !< to the straight-fieldline angle;
                                     !< only relevant if \c Lsparse.eq.2; \see tr00ab() for details
                                     !< <ul>
                                     !< <li> if \c iprecon = 0, the preconditioner is `N' </li>
                                     !< <li> if \c iprecon = 1, the preconditioner is `J' </li>
                                     !< <li> if \c iprecon = 2, the preconditioner is `S' </li>
                                     !< </ul>
  REAL         :: iotatol     = -1.0 !< tolerance required for iterative construction of straight-fieldline angle;
                                     !< only relevant if \c Lsparse.ge.2
  INTEGER      :: Lextrap     =  0   !< geometry of innermost interface is defined by extrapolation
  INTEGER      :: Mregular    = -1   !< maximum regularization factor
                                     !< <ul>
                                     !< <li> if \c Mregular.ge.2, then \c regumm \f$_i\f$ = \c Mregular \f$/ 2 \f$ where \c m \f$_i > \f$ \c Mregular </li>
                                     !< </ul>
  INTEGER      :: Lrzaxis     =  1   !< controls the guess of geometry axis in the innermost volume or initialization of interfaces
                                     !< <ul>
                                     !< <li> if \c iprecon = 1, the centroid is used </li>
                                     !< <li> if \c iprecon = 2, the Jacobian \f$m=1\f$ harmonic elimination method is used </li>
                                     !< </ul>
  INTEGER      :: Ntoraxis    =  3   !< the number of \f$n\f$ harmonics used in the Jacobian \f$m=1\f$ harmonic elimination method;
                                     !< only relevant if \c Lrzaxis.ge.1 .
!> @}

!> \addtogroup grp_global_local locallist
!> \brief The namelist \c locallist controls the construction of the Beltrami fields in each volume.
!>
!> The transformation to straight-fieldline coordinates is singular when the rotational-transform of the interfaces is rational;
!> however, the rotational-transform is still well defined.
!> @{
  INTEGER      :: LBeltrami  =  4   !< Control flag for solution of Beltrami equation
                                    !< <ul>
                                    !< <li> if \c LBeltrami = 1,3,5 or 7, (SQP) then the Beltrami field in each volume is constructed
                                    !<      by minimizing the magnetic energy with the constraint of fixed helicity;
                                    !<      this is achieved by using sequential quadratic programming as provided by \c E04UFF .
                                    !<      This approach has the benefit (in theory) of robustly constructing minimum energy solutions
                                    !<      when multiple, i.e. bifurcated, solutions exist.
                                    !< <li> if \c LBeltrami = 2,3,6 or 7, (Newton) then the Beltrami fields are constructed by employing a standard Newton method
                                    !<      for locating an extremum of
                                    !<      \f$F\equiv \int B^2 dv - \mu (\int {\bf A}\cdot{\bf B}dv-{\cal K})\f$,
                                    !<      where \f$\mu\f$ is treated as an independent degree of freedom similar to the parameters describing the vector potential
                                    !<      and \f${\cal K}\f$ is the required value of the helicity;
                                    !<      this is the standard Lagrange multipler approach for locating the constrained minimum;
                                    !<      this method cannot distinguish saddle-type extrema from minima, and which solution that will be obtained depends on the initial guess;
                                    !< <li> if \c LBeltrami = 4,5,6 or 7, (linear) it is assumed that the Beltrami fields are parameterized by \f$\mu\f$;
                                    !<      in this case, it is only required to solve \f$\nabla \times {\bf B} = \mu {\bf B}\f$ which reduces to a system of linear equations;
                                    !<      \f$\mu\f$ may or may not be adjusted iteratively, depending on \c Lconstraint,
                                    !<      to satisfy either rotational-transform or helicity constraints;
                                    !< <li> for flexibility and comparison, each of the above methods can be employed; for example:
                                    !<      <ul>
                                    !<      <li> if \c LBeltrami = 1, only the SQP    method will be employed;
                                    !<      <li> if \c LBeltrami = 2, only the Newton method will be employed;
                                    !<      <li> if \c LBeltrami = 4, only the linear method will be employed;
                                    !<      <li> if \c LBeltrami = 3, the SQP and the Newton method are used;
                                    !<      <li> if \c LBeltrami = 5, the SQP and the linear method are used;
                                    !<      <li> if \c LBeltrami = 6, the Newton and the linear method are used;
                                    !<      <li> if \c LBeltrami = 7, all three methods will be employed;
                                    !<      </ul>
                                    !< </ul>
  INTEGER      :: Linitgues  =  1   !< controls how initial guess for Beltrami field is constructed
                                    !< <ul>
                                    !< <li> only relevant for routines that require an initial guess for the Beltrami fields, such as the SQP and Newton methods,
                                    !<      or the sparse linear solver;
                                    !< <li> if \c Linitgues = 0, the initial guess for the Beltrami field is trivial
                                    !< <li> if \c Linitgues = 1, the initial guess for the Beltrami field is an integrable approximation
                                    !< <li> if \c Linitgues = 2, the initial guess for the Beltrami field is read from file
                                    !< <li> if \c Linitgues = 3, the initial guess for the Beltrami field will be randomized with the maximum \c maxrndgues
                                    !< </ul>
  INTEGER      :: Lposdef    =  0   !< redundant;
  REAL         :: maxrndgues =  1.0 !< the maximum random number of the Beltrami field if \c Linitgues = 3
  INTEGER      :: Lmatsolver =  3     !< 1 for LU factorization, 2 for GMRES, 3 for GMRES matrix-free
  INTEGER      :: NiterGMRES =  200   !< number of max iteration for GMRES
  REAL         :: epsGMRES   =  1e-14 !< the precision of GMRES
  INTEGER      :: LGMRESprec =  1     !< type of preconditioner for GMRES, 1 for ILU sparse matrix
  REAL         :: epsILU     =  1e-12 !< the precision of incomplete LU factorization for preconditioning
!> @}

!> \addtogroup grp_global_global globallist
!> \brief The namelist \c globallist controls the search for global force-balance.
!>
!> Comments:
!> <ul>
!> <li> The "force" vector, \f${\bf F}\f$, which is constructed in dforce(), is a combination of pressure-imbalance Fourier harmonics,
!>       \f{eqnarray}{ F_{i,v} \equiv [[ p+B^2/2 ]]_{i,v} \times \exp\left[-\texttt{escale}(m_i^2+n_i^2) \right] \times \texttt{opsilon},
!>       \label{eq:forcebalancemn_global} \f}
!>       and spectral-condensation constraints, \f$I_{i,v}\f$, and the "star-like" angle constraints, \f$S_{i,v,}\f$, (see lforce() for details)
!>       \f{eqnarray}{ F_{i,v} \equiv \texttt{epsilon} \times I_{i,v}
!>           + \texttt{upsilon} \times \left( \psi_v^\omega S_{i,v,1} - \psi_{v+1}^\omega S_{i,v+1,0} \right),
!>       \label{eq:spectralbalancemn_global} \f}
!>       where \f$\psi_v\equiv\f$ normalized toroidal flux, \c tflux, and \f$\omega\equiv\f$ \c wpoloidal. </li>
!> </ul>
!> @{
  INTEGER      :: Lfindzero  =   0       !< use Newton methods to find zero of force-balance, which is computed by dforce()
                                         !< <ul>
                                         !< <li> if \c Lfindzero = 0, then dforce() is called once
                                         !<      to compute the Beltrami fields consistent with the given geometry and constraints </li>
                                         !< <li> if \c Lfindzero = 1, then call \c C05NDF (uses   function values only), which iteratively calls dforce() </li>
                                         !< <li> if \c Lfindzero = 2, then call \c C05PDF (uses derivative information), which iteratively calls dforce() </li>
                                         !< </ul>
  REAL         :: escale     =   0.0     !< controls the weight factor, \c BBweight, in the force-imbalance harmonics
                                         !< <ul>
                                         !< <li> \c BBweight(i) \f$\displaystyle \equiv \texttt{opsilon} \times \exp\left[-\texttt{escale} \times (m_i^2+n_i^2) \right]\f$ </li>
                                         !< <li> defined in preset() ; used in dforce() </li>
                                         !< <li> also see Eqn.\f$(\ref{eq:forcebalancemn_global})\f$ </li>
                                         !< </ul>
  REAL         :: opsilon    =   1.0     !< weighting of force-imbalance
                                         !< <ul>
                                         !< <li> used in dforce(); also see Eqn.\f$(\ref{eq:forcebalancemn_global})\f$ </li>
                                         !< </ul>
  REAL         :: pcondense  =   2.0     !< spectral condensation parameter
                                         !< <ul>
                                         !< <li> used in preset() to define \c mmpp(i) \f$\equiv m_i^p\f$, where \f$p\equiv \f$ \c pcondense </li>
                                         !< <li> the angle freedom is exploited to minimize \f$\displaystyle \texttt{epsilon} \sum_{i} m_i^p (R_{i}^2+Z_{i}^2)\f$
                                         !<      with respect to tangential variations in the interface geometry </li>
                                         !< <li> also see Eqn.\f$(\ref{eq:spectralbalancemn_global})\f$ </li>
                                         !< </ul>
  REAL         :: epsilon    =   0.0     !< weighting of spectral-width constraint
                                         !< <ul>
                                         !< <li> used in dforce(); also see Eqn.\f$(\ref{eq:spectralbalancemn_global})\f$ </li>
                                         !< </ul>
  REAL         :: wpoloidal  =   1.0     !< "star-like" poloidal angle constraint radial exponential factor
                                         !< used in preset() to construct \c sweight
  REAL         :: upsilon    =   1.0     !< weighting of "star-like" poloidal angle constraint
                                         !< used in preset() to construct \c sweight
  REAL         :: forcetol   =   1.0e-10 !< required tolerance in force-balance error; only used as an initial check
                                         !< <ul>
                                         !< <li> if the initially supplied interfaces are consistent with force-balance to within \c forcetol
                                         !<      then the geometry of the interfaces is not altered </li>
                                         !< <li> if not, then the geometry of the interfaces is changed in order to bring the configuration into force balance
                                         !<      so that the geometry of interfaces is within \c c05xtol, defined below, of the true solution </li>
                                         !< <li> to force execution of either \c C05NDF or \c C05PDF, regardless of the initial force imbalance,
                                         !<      set \c forcetol < 0 </li>
                                         !< </ul>
  REAL         :: c05xmax    =   1.0e-06 !< required tolerance in position, \f${\bf x} \equiv \{ R_{i,v}, Z_{i,v}\}\f$
  REAL         :: c05xtol    =   1.0e-12 !< required tolerance in position, \f${\bf x} \equiv \{ R_{i,v}, Z_{i,v}\}\f$
                                         !< <ul>
                                         !< <li> used by both \c C05NDF and \c C05PDF; see the NAG documents for further details on how the error is defined </li>
                                         !< <li> constraint \c c05xtol > 0.0 </li>
                                         !< </ul>
  REAL         :: c05factor  =   1.0e-02 !< used to control initial step size in
                                         !<       \c C05NDF and \c C05PDF
                                         !< <ul>
                                         !< <li> constraint \c c05factor > 0.0 </li>
                                         !< <li> only relevant if \c Lfindzero > 0 </li>
                                         !< </ul>
  LOGICAL      :: LreadGF    =  .true.   !< read \f$\nabla_{\bf x} {\bf F}\f$ from file \c ext.GF
                                         !< <ul>
                                         !< <li> only used if \c Lfindzero = 2 </li>
                                         !< <li> only used in newton() </li>
                                         !< </ul>
  INTEGER      :: mfreeits   =   0       !< maximum allowed free-boundary iterations
                                         !< <ul>
                                         !< <li> only used if \c Lfreebound = 1 </li>
                                         !< <li> only used in xspech() </li>
                                         !< </ul>
  REAL         :: bnstol     =   1.0e-06 !< redundant;
  REAL         :: bnsblend   =   0.666   !< redundant;
  REAL         :: gBntol     =   1.0e-06 !< required tolerance in free-boundary iterations
                                         !< <ul>
                                         !< <li> only used if \c Lfreebound = 1 </li>
                                         !< <li> only used in xspech() </li>
                                         !< </ul>
  REAL         :: gBnbld     =   0.666   !< normal blend
                                         !< <ul>
                                         !< <li> The "new" magnetic field at the computational boundary produced by the plasma currents is updated using a Picard scheme:
                                         !<      \f{eqnarray}{ ({\bf B}\cdot{\bf n})^{j+1} =    \texttt{gBnbld}  \times ({\bf B}\cdot{\bf n})^{j}
                                         !<                                      + (1-\texttt{gBnbld}) \times ({\bf B}\cdot{\bf n})^{*},
                                         !<      \f}
                                         !<      where \f$j\f$ labels free-boundary iterations, and \f$({\bf B}\cdot{\bf n})^{*}\f$ is computed by virtual casing. </li>
                                         !< <li> only used if \c Lfreebound = 1 </li>
                                         !< <li> only used in xspech() </li>
                                         !< </ul>
  REAL         :: vcasingeps =   1.e-12  !< regularization of Biot-Savart; see bnorml(), casing()
  REAL         :: vcasingtol =   1.e-08  !< accuracy on virtual casing integral; see bnorml(), casing()
  INTEGER      :: vcasingits =   8       !< minimum number of calls to adaptive virtual casing routine; see casing()
  INTEGER      :: vcasingper =   1       !< periods of integragion  in adaptive virtual casing routine; see casing()
  INTEGER      :: mcasingcal =   8       !< minimum number of calls to adaptive virtual casing routine; see casing(); redundant;
!> @}


!> \addtogroup grp_global_diagnostics diagnosticslist
!> \brief The namelist \c diagnosticslist controls post-processor diagnostics, such as Poincaré  plot resolution, etc.
!> @{
  REAL         :: odetol           =     1.0e-07 !< o.d.e. integration tolerance for all field line tracing routines
  REAL         :: absreq           =     1.0e-08 !< redundant
  REAL         :: relreq           =     1.0e-08 !< redundant
  REAL         :: absacc           =     1.0e-04 !< redundant
  REAL         :: epsr             =     1.0e-08 !< redundant
  INTEGER      :: nPpts            =     0       !< number of toroidal transits used (per trajectory) in following field lines
                                                 !< for constructing Poincaré plots;
                                                 !< if \c nPpts<1, no Poincaré plot is constructed;
  REAL         :: Ppts             =     0.0     !< stands for Poincare plot theta start. Chose at which angle (normalized over \f$\pi\f$) the Poincare field-line tracing start.
  INTEGER      :: nPtrj(1:MNvol+1) =    -1       !< number of trajectories in each annulus to be followed in constructing Poincaré plot
                                                 !< <ul>
                                                 !< <li> if \c nPtrj(l)<0, then \c nPtrj(l) = Ni(l),
                                                 !<       where \c Ni(l) is the grid resolution used to construct the Beltrami field in volume \f$l\f$ </li>
                                                 !< </ul>
  LOGICAL      :: LHevalues        =  .false.    !< to compute eigenvalues of \f$\nabla {\bf F}\f$
  LOGICAL      :: LHevectors       =  .false.    !< to compute eigenvectors (and also eigenvalues) of \f$\nabla {\bf F}\f$
  LOGICAL      :: LHmatrix         =  .false.    !< to compute and write to file the elements of \f$\nabla {\bf F}\f$
  INTEGER      :: Lperturbed       =     0       !< to compute linear, perturbed equilibrium
  INTEGER      :: dpp              =    -1       !< perturbed harmonic
  INTEGER      :: dqq              =    -1       !< perturbed harmonic
  INTEGER      :: Lerrortype       =     0       !< the type of error output for Lcheck=1
  INTEGER      :: Ngrid            =    -1       !< the number of points to output in the grid, -1 for Lrad(vvol)
  REAL         :: dRZ              =     1E-5    !< difference in geometry for finite difference estimate (debug only)
  INTEGER      :: Lcheck           =     0       !< implement various checks
                                                 !< <ul>
                                                 !< <li> if \c Lcheck = 0, no additional check on the calculation is performed </li>
                                                 !< <li> if \c Lcheck = 1, the error in the current, i.e. \f$\nabla\times{\bf B}-\mu{\bf B}\f$ is computed as a post-diagnostic </li>
                                                 !< <li> if \c Lcheck = 2, the analytic derivatives of the interface transform w.r.t.
                                                 !<      the helicity multiplier, \f$\mu\f$, and the enclosed poloidal flux, \f$\Delta\psi_p\f$, are compared to a finite-difference estimate
                                                 !<      <ul>
                                                 !<      <li> only if \c Lconstraint==1 </li>
                                                 !<      <li> only for \c dspec executable, i.e. must compile with \c DFLAGS = "-D DEBUG" </li>
                                                 !<      </ul> </li>
                                                 !< <li> if \c Lcheck = 3, the analytic derivatives of the volume w.r.t. interface Fourier harmonic
                                                 !<      is compared to a finite-difference estimate
                                                 !<      <ul>
                                                 !<      <li> must set \c Lfindzero\f$ = 2\f$ </li>
                                                 !<      <li> set \c forcetol sufficiently small and set \c LreadGF = F,
                                                 !<           so that the matrix of second derivatives is calculated </li>
                                                 !<      <li> only for \c dspec executable, i.e. must compile with \c DFLAGS = "-D DEBUG" </li>
                                                 !<      </ul> </li>
                                                 !< <li> if \c Lcheck = 4, the analytic calculation of the derivatives of the magnetic field, \f$B^2\f$, at the interfaces
                                                 !<      is compared to a finite-difference estimate
                                                 !<      <ul>
                                                 !<      <li> must set \c Lfindzero\f$ = 2\f$ </li>
                                                 !<      <li> set \c forcetol sufficiently small </li>
                                                 !<      <li> set \c LreadGF=F </li>
                                                 !<      <li> only for \c dspec executable, i.e. must compile with \c DFLAGS = "-D DEBUG" </li>
                                                 !<      </ul> </li>
                                                 !< <li> if \c Lcheck = 5, the analytic calculation of the matrix of the derivatives of the force imbalance
                                                 !<      is compared to a finite-difference estimate </li>
                                                 !< <li> if \c Lcheck = 6, the virtual casing calculation is compared to \c xdiagno (Lazerson 2013 \cite y2013_lazerson)
                                                 !<      <ul>
                                                 !<      <li> the input file for \c xdiagno is written by bnorml() </li>
                                                 !<      <li> this provides the Cartesian coordinates on the computational boundary where the virtual casing routine casing()
                                                 !<           computes the magnetic field, with the values of the magnetic field being written to the screen for comparison </li>
                                                 !<      <li> must set \c Freebound=1, \c Lfindzero>0, \c mfreeits!=0 </li>
                                                 !<      <li> \c xdiagno must be executed manually </li>
                                                 !<      </ul> </li>
                                                 !< </ul>
  LOGICAL      :: Ltiming          =  .false.    !< to check timing
  REAL         :: fudge            =     1.0e-00 !< redundant
  REAL         :: scaling          =     1.0e-00 !< redundant
!> @}


!> \addtogroup grp_global_screenlist screenlist
!> \brief The namelist \c screenlist controls screen output.
!> Every subroutine, e.g. \c xy00aa.h, has its own write flag, \c Wxy00aa.
!> @{
  LOGICAL :: Wmanual  = .false.
  LOGICAL :: Wrzaxis  = .false.
  LOGICAL :: Wpackxi  = .false.
  LOGICAL :: Wvolume  = .false.
  LOGICAL :: Wcoords  = .false.
  LOGICAL :: Wbasefn  = .false.
  LOGICAL :: Wmemory  = .false.
  LOGICAL :: Wmetrix  = .false.
  LOGICAL :: Wma00aa  = .false.
  LOGICAL :: Wmatrix  = .false.
  LOGICAL :: Wspsmat  = .false.
  LOGICAL :: Wspsint  = .false.
  LOGICAL :: Wmp00ac  = .false.
  LOGICAL :: Wma02aa  = .false.
  LOGICAL :: Wpackab  = .false.
  LOGICAL :: Wtr00ab  = .false.
  LOGICAL :: Wcurent  = .false.
  LOGICAL :: Wdf00ab  = .false.
  LOGICAL :: Wlforce  = .false.
  LOGICAL :: Wintghs  = .false.
  LOGICAL :: Wmtrxhs  = .false.
  LOGICAL :: Wlbpol   = .false.
  LOGICAL :: Wbrcast  = .false.
  LOGICAL :: Wdfp100  = .false.
  LOGICAL :: Wdfp200  = .false.
  LOGICAL :: Wdforce  = .false.
  LOGICAL :: Wnewton  = .false.
  LOGICAL :: Wcasing  = .false.
  LOGICAL :: Wbnorml  = .false.
  LOGICAL :: Wjo00aa  = .false.
  LOGICAL :: Wpp00aa  = .false.
  LOGICAL :: Wpp00ab  = .false.
  LOGICAL :: Wbfield  = .false.
  LOGICAL :: Wstzxyz  = .false.
  LOGICAL :: Whesian  = .false.
  LOGICAL :: Wra00aa  = .false.
  LOGICAL :: Wnumrec  = .false.
  LOGICAL :: Wdcuhre  = .false.
  LOGICAL :: Wminpack = .false.
  LOGICAL :: Wiqpack  = .false.
  LOGICAL :: Wrksuite = .false.
  LOGICAL :: Wi1mach  = .false.
  LOGICAL :: Wd1mach  = .false.
  LOGICAL :: Wilut    = .false.
  LOGICAL :: Witers   = .false.
  LOGICAL :: Wsphdf5  = .false.
  LOGICAL :: Wpreset  = .false.
  LOGICAL :: Wglobal  = .false.
  LOGICAL :: Wxspech  = .false.
  LOGICAL :: Wbuild_vector_potential = .false. !< \todo: what is this?
  LOGICAL :: Wreadin  = .false. !< write screen output of readin()
  LOGICAL :: Wwrtend  = .false. !< write screen output of wrtend()
  LOGICAL :: Wmacros  = .false. !< write screen output from expanded macros
!> @}

  namelist/physicslist/&
 Igeometry   ,&
 Istellsym   ,&
 Lfreebound  ,&
 phiedge     ,&
 curtor      ,&
 curpol      ,&
 gamma       ,&
 Nfp         ,&
 Nvol        ,&
 Mpol        ,&
 Ntor        ,&
 Lrad        ,&
 Lconstraint ,&
 tflux       ,&
 pflux       ,&
 helicity    ,&
 pscale      ,&
 pressure    ,&
 Ladiabatic  ,&
 adiabatic   ,&
 mu          ,&
 Ivolume     ,&
 Isurf       ,&
 pl          ,&
 ql          ,&
 pr          ,&
 qr          ,&
 iota        ,&
 lp          ,&
 lq          ,&
 rp          ,&
 rq          ,&
 oita        ,&
 mupftol     ,&
 mupfits     ,&
 rpol        ,&
 rtor        ,&
 Lreflect    ,&
 Rac         ,&
 Zas         ,&
 Ras         ,&
 Zac         ,&
 Rbc         ,&
 Zbs         ,&
 Rbs         ,&
 Zbc         ,&
 Rwc         ,&
 Zws         ,&
 Rws         ,&
 Zwc         ,&
 Vns         ,&
 Bns         ,&
 Vnc         ,&
 Bnc

  namelist/numericlist/&
 Linitialize ,&
 LautoinitBn ,&
 Lzerovac    ,&
 Ndiscrete   ,&
 Nquad       ,&
 iMpol       ,&
 iNtor       ,&
 Lsparse     ,&
 Lsvdiota    ,&
 imethod     ,&
 iorder      ,&
 iprecon     ,&
 iotatol     ,&
 Lextrap     ,&
 Mregular    ,&
 Lrzaxis     ,&
 Ntoraxis

  namelist/locallist/&
 LBeltrami   ,&
 Linitgues   ,&
 maxrndgues  ,&
 maxrndgues  ,&
 Lmatsolver  ,&
 NiterGMRES  ,&
 epsGMRES    ,&
 LGMRESprec  ,&
 epsILU      ,&
 Lposdef

  namelist/globallist/&
 Lfindzero   ,&
 escale      ,&
 opsilon     ,&
 pcondense   ,&
 epsilon     ,&
 wpoloidal   ,&
 upsilon     ,&
 forcetol    ,&
 c05xmax     ,&
 c05xtol     ,&
 c05factor   ,&
 LreadGF     ,&
 mfreeits    ,&
 bnstol      ,&
 bnsblend    ,&
 gBntol      ,&
 gBnbld      ,&
 vcasingeps  ,&
 vcasingtol  ,&
 vcasingits  ,&
 vcasingper  ,&
 mcasingcal

  namelist/diagnosticslist/&
 odetol     ,&
 absreq     ,&
 relreq     ,&
 absacc     ,&
 epsr       ,&
 nPpts      ,&
 Ppts       ,&
 nPtrj      ,&
 LHevalues  ,&
 LHevectors ,&
 LHmatrix   ,&
 Lperturbed ,&
 dpp        ,&
 dqq        ,&
 Lerrortype ,&
 Ngrid      ,&
 Lcheck     ,&
 dRZ        ,&
 Ltiming    ,&
 fudge      ,&
 scaling

  namelist/screenlist/&
 Wmanual , &
 Wrzaxis , &
 Wpackxi , &
 Wvolume , &
 Wcoords , &
 Wbasefn , &
 Wmemory , &
 Wmetrix , &
 Wma00aa , &
 Wmatrix , &
 Wspsmat , &
 Wspsint , &
 Wmp00ac , &
 Wma02aa , &
 Wpackab , &
 Wtr00ab , &
 Wcurent , &
 Wdf00ab , &
 Wlforce , &
 Wintghs , &
 Wmtrxhs , &
 Wlbpol  , &
 Wbrcast , &
 Wdfp100 , &
 Wdfp200 , &
 Wdforce , &
 Wnewton , &
 Wcasing , &
 Wbnorml , &
 Wjo00aa , &
 Wpp00aa , &
 Wpp00ab , &
 Wbfield , &
 Wstzxyz , &
 Whesian , &
 Wra00aa , &
 Wnumrec , &
 Wdcuhre , &
 Wminpack, &
 Wiqpack , &
 Wrksuite, &
 Wi1mach , &
 Wd1mach , &
 Wilut   , &
 Witers  , &
 Wsphdf5 , &
 Wpreset , &
 Wglobal , &
 Wxspech , &
 Wbuild_vector_potential , &
 Wreadin , &
 Wwrtend , &
 Wmacros

  contains

! The task of this subroutine is to initialize *all* input variables
! to a default known state. A subsequent call to readin() can be used
! to fill in the values from a text input file containing the
! appropriate input namelists. Alternatively, after calling this routine,
! the input variables can be modified from their default values in-memory,
! e.g., via the Python wrapper.
!
! Note that since the variables in screenlist are auto-generated,
! they do not appear here and also are not touched by this routine.
subroutine initialize_inputs

  implicit none

! physicslist

  Igeometry                  =  3
  Istellsym                  =  1
  Lfreebound                 =  0
  phiedge                    =  1.0
  curtor                     =  0.0
  curpol                     =  0.0
  gamma                      =  0.0
  Nfp                        =  1
  Nvol                       =  1
  Mpol                       =  0
  Ntor                       =  0
  Lrad(1:MNvol+1)            =  4
  Lconstraint                = -1
  tflux(1:MNvol+1)           =  0.0
      pflux(1:MNvol+1)       =  0.0
   helicity(1:MNvol)         =  0.0
  pscale                     =  0.0
   pressure(1:MNvol+1)       =  0.0
  Ladiabatic                 =  0
  adiabatic(1:MNvol+1)       =  0.0
         mu(1:MNvol+1)       =  0.0
  Ivolume(1:MNvol+1)         =  0.0
  Isurf(1:MNvol)             =  0.0
         pl(0:MNvol)         =  0
         ql(0:MNvol)         =  0
         pr(0:MNvol)         =  0
         qr(0:MNvol)         =  0
       iota(0:MNvol)         =  0.0
         lp(0:MNvol)         =  0
         lq(0:MNvol)         =  0
         rp(0:MNvol)         =  0
         rq(0:MNvol)         =  0
       oita(0:MNvol)         =  0.0
  rpol                       =  1.0
  rtor                       =  1.0

  Rac(     0:MNtor        )  =  0.0
  Zas(     0:MNtor        )  =  0.0
  Ras(     0:MNtor        )  =  0.0
  Zac(     0:MNtor        )  =  0.0

  Rbc(-MNtor:MNtor,-MMpol:MMpol)  =  0.0
  Zbs(-MNtor:MNtor,-MMpol:MMpol)  =  0.0
  Rbs(-MNtor:MNtor,-MMpol:MMpol)  =  0.0
  Zbc(-MNtor:MNtor,-MMpol:MMpol)  =  0.0

  Rwc(-MNtor:MNtor,-MMpol:MMpol)  =  0.0
  Zws(-MNtor:MNtor,-MMpol:MMpol)  =  0.0
  Rws(-MNtor:MNtor,-MMpol:MMpol)  =  0.0
  Zwc(-MNtor:MNtor,-MMpol:MMpol)  =  0.0

  Vns(-MNtor:MNtor,-MMpol:MMpol)  =  0.0
  Bns(-MNtor:MNtor,-MMpol:MMpol)  =  0.0
  Vnc(-MNtor:MNtor,-MMpol:MMpol)  =  0.0
  Bnc(-MNtor:MNtor,-MMpol:MMpol)  =  0.0

  mupftol                    =  1.0e-14
  mupfits                    =  8

  Lreflect                   =  0

! numericlist

  Linitialize =  0
  LautoinitBn =  1
  Lzerovac    =  0
  Ndiscrete   =  2
  Nquad       = -1
  iMpol       = -4
  iNtor       = -4
  Lsparse     =  0
  Lsvdiota    =  0
  imethod     =  3
  iorder      =  2
  iprecon     =  0
  iotatol     = -1.0
  Lextrap     =  0
  Mregular    = -1
  Lrzaxis     = 1
  Ntoraxis    = 3

! locallist

  LBeltrami  =  4
  Linitgues  =  1
  Lposdef    =  0
  maxrndgues =  1.0
  Lmatsolver = 3
  NiterGMRES = 200
  epsGMRES   = 1e-14
  LGMRESprec = 1
  epsILU     = 1e-12

! globallist

  Lfindzero  =   0
  escale     =   0.0
  opsilon    =   1.0
  pcondense  =   2.0
  epsilon    =   0.0
  wpoloidal  =   1.0
  upsilon    =   1.0
  forcetol   =   1.0e-10
  c05xmax    =   1.0e-06
  c05xtol    =   1.0e-12
  c05factor  =   1.0e-02
  LreadGF    =  .true.
  mfreeits   =   0
  bnstol     =   1.0e-06
  bnsblend   =   0.666
  gBntol     =   1.0e-06
  gBnbld     =   0.666
  vcasingeps =   1.e-12
  vcasingtol =   1.e-08
  vcasingits =   8
  vcasingper =   1
  mcasingcal =   8

! diagnosticslist

  odetol           =     1.0e-07
  absreq           =     1.0e-08
  relreq           =     1.0e-08
  absacc           =     1.0e-04
  epsr             =     1.0e-08
  nPpts            =     0
  Ppts             =     0.0
  nPtrj(1:MNvol+1) =    -1
  LHevalues        =  .false.
  LHevectors       =  .false.
  LHmatrix         =  .false.
  Lperturbed       =     0
  dpp              =    -1
  dqq              =    -1
  Lerrortype       =     0
  Ngrid            =    -1
  dRZ              =     1E-5
  Lcheck           =     0
  Ltiming          =  .false.
  fudge            =     1.0e-00
  scaling          =     1.0e-00

end subroutine initialize_inputs

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end module inputlist
!> @}
