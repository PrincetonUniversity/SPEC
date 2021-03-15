! This is a pure data-container holding *all* input variables
! required for a SPEC calculation.
module inputlist

!latex \subsection{input namelists}

!latex \begin{enumerate}

  implicit none

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item The input file, \verb+ext.sp+, where \verb+ext*100+ or \verb+ext.sp*100+ is given as command line input, contains the following namelists and interface geometry.

!SET MAXIMUM RESOLUTION;

  INTEGER, parameter :: MNvol     = 256 !latex \item The maximum value of \inputvar{Nvol} is \verb+MNvol=256+.
  INTEGER, parameter :: MMpol     =  64 !latex \item The maximum value of \inputvar{Mpol} is \verb+MNpol= 32+.
  INTEGER, parameter :: MNtor     =  64 !latex \item The maximum value of \inputvar{Ntor} is \verb+MNtor= 32+.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! the following variables constitute the namelist/physicslist/; note that all variables in namelist need to be broadcasted in readin;

  INTEGER      :: Igeometry
  INTEGER      :: Istellsym
  INTEGER      :: Lfreebound
  REAL         :: phiedge
  REAL         :: curtor
  REAL         :: curpol
  REAL         :: gamma
  INTEGER      :: Nfp
  INTEGER      :: Nvol
  INTEGER      :: Mpol
  INTEGER      :: Ntor
  INTEGER      :: Lrad(1:MNvol+1)
  INTEGER      :: Lconstraint
  REAL         :: tflux(1:MNvol+1)
  REAL         ::     pflux(1:MNvol+1)
  REAL         ::  helicity(1:MNvol)
  REAL         :: pscale
  REAL         ::  pressure(1:MNvol+1)
  INTEGER      :: Ladiabatic
  REAL         :: adiabatic(1:MNvol+1)
  REAL         ::        mu(1:MNvol+1)
  REAL         :: Ivolume(1:MNvol+1)
  REAL         :: Isurf(1:MNvol)
  INTEGER      ::        pl(0:MNvol)
  INTEGER      ::        ql(0:MNvol)
  INTEGER      ::        pr(0:MNvol)
  INTEGER      ::        qr(0:MNvol)
  REAL         ::      iota(0:MNvol)
  INTEGER      ::        lp(0:MNvol)
  INTEGER      ::        lq(0:MNvol)
  INTEGER      ::        rp(0:MNvol)
  INTEGER      ::        rq(0:MNvol)
  REAL         ::      oita(0:MNvol)
  REAL         :: rpol
  REAL         :: rtor

  REAL         :: Rac(     0:MNtor        ) !     stellarator symmetric coordinate axis;
  REAL         :: Zas(     0:MNtor        )
  REAL         :: Ras(     0:MNtor        ) ! non-stellarator symmetric coordinate axis;
  REAL         :: Zac(     0:MNtor        )

  REAL         :: Rbc(-MNtor:MNtor,-MMpol:MMpol) !     stellarator symmetric boundary components;
  REAL         :: Zbs(-MNtor:MNtor,-MMpol:MMpol) !     stellarator symmetric boundary components;
  REAL         :: Rbs(-MNtor:MNtor,-MMpol:MMpol) ! non-stellarator symmetric boundary components;
  REAL         :: Zbc(-MNtor:MNtor,-MMpol:MMpol) ! non-stellarator symmetric boundary components;

  REAL         :: Rwc(-MNtor:MNtor,-MMpol:MMpol) !     stellarator symmetric boundary components of wall;
  REAL         :: Zws(-MNtor:MNtor,-MMpol:MMpol) !     stellarator symmetric boundary components of wall;
  REAL         :: Rws(-MNtor:MNtor,-MMpol:MMpol) ! non-stellarator symmetric boundary components of wall;
  REAL         :: Zwc(-MNtor:MNtor,-MMpol:MMpol) ! non-stellarator symmetric boundary components of wall;

  REAL         :: Vns(-MNtor:MNtor,-MMpol:MMpol) !     stellarator symmetric normal field at boundary; vacuum component;
  REAL         :: Bns(-MNtor:MNtor,-MMpol:MMpol) !     stellarator symmetric normal field at boundary; plasma component;
  REAL         :: Vnc(-MNtor:MNtor,-MMpol:MMpol) ! non-stellarator symmetric normal field at boundary; vacuum component;
  REAL         :: Bnc(-MNtor:MNtor,-MMpol:MMpol) ! non-stellarator symmetric normal field at boundary; plasma component;

  REAL         :: mupftol
  INTEGER      :: mupfits

  INTEGER      :: Lreflect

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! the following variables constitute the namelist/numericlist/; note that all variables in namelist need to be broadcasted in readin;

  INTEGER      :: Linitialize
  INTEGER      :: LautoinitBn
  INTEGER      :: Lzerovac
  INTEGER      :: Ndiscrete
  INTEGER      :: Nquad
  INTEGER      :: iMpol
  INTEGER      :: iNtor
  INTEGER      :: Lsparse
  INTEGER      :: Lsvdiota
  INTEGER      :: imethod
  INTEGER      :: iorder
  INTEGER      :: iprecon
  REAL         :: iotatol
  INTEGER      :: Lextrap
  INTEGER      :: Mregular
  INTEGER      :: Lrzaxis
  INTEGER      :: Ntoraxis

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! the following variables constitute the namelist/locallist/; note that all variables in namelist need to be broadcasted in readin;

  INTEGER      :: LBeltrami
  INTEGER      :: Linitgues
  INTEGER      :: Lposdef
  REAL         :: maxrndgues
  INTEGER      :: Lmatsolver
  INTEGER      :: NiterGMRES
  REAL         :: epsGMRES
  INTEGER      :: LGMRESprec
  REAL         :: epsILU

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! the following variables constitute the namelist/globallist/; note that all variables in namelist need to be broadcasted in readin;

  INTEGER      :: Lfindzero
  REAL         :: escale
  REAL         :: opsilon
  REAL         :: pcondense
  REAL         :: epsilon
  REAL         :: wpoloidal
  REAL         :: upsilon
  REAL         :: forcetol
  REAL         :: c05xmax
  REAL         :: c05xtol
  REAL         :: c05factor
  LOGICAL      :: LreadGF
  INTEGER      :: mfreeits
  REAL         :: bnstol     ! redundant;
  REAL         :: bnsblend   ! redundant;
  REAL         :: gBntol
  REAL         :: gBnbld
  REAL         :: vcasingeps
  REAL         :: vcasingtol
  INTEGER      :: vcasingits
  INTEGER      :: vcasingper
  INTEGER      :: mcasingcal ! redundant;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! the following variables constitute the namelist/diagnosticslist/; note that all variables in namelist need to be broadcasted in readin;

  REAL         :: odetol
  REAL         :: absreq           ! redundant;
  REAL         :: relreq           ! redundant;
  REAL         :: absacc           ! redundant;
  REAL         :: epsr             ! redundant;
  INTEGER      :: nPpts
  REAL         :: Ppts
  INTEGER      :: nPtrj(1:MNvol+1)
  LOGICAL      :: LHevalues
  LOGICAL      :: LHevectors
  LOGICAL      :: LHmatrix
  INTEGER      :: Lperturbed
  INTEGER      :: dpp
  INTEGER      :: dqq
  INTEGER      :: Lerrortype
  INTEGER      :: Ngrid
  REAL         :: dRZ              ! For finite difference estimate
  INTEGER      :: Lcheck
  LOGICAL      :: Ltiming
  REAL         :: fudge            ! redundant;
  REAL         :: scaling          ! redundant;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! the following variables constitute the namelist/screenlist/; note that all variables in namelist need to be broadcasted in readin;

! DSCREENLIST ! define screenlist; this is expanded by Makefile; DO NOT REMOVE; each file compiled by Makefile has its own write flag;
  LOGICAL      :: Wbuild_vector_potential = .false.
  LOGICAL      :: Wreadin = .false.
  LOGICAL      :: Wwritin = .false. ! redundant;
  LOGICAL      :: Wwrtend = .false.
  LOGICAL      :: Wmacros = .false.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item In the following, all default settings are shown.
!latex \end{enumerate}

!latex \subsubsection{\type{physicslist} : }

!latex \begin{enumerate}
!latex \item The namelist \verb+physicslist+ controls the geometry, profiles, and numerical resolution.
!latex       \bi
!latex       \verb+namelist/physicslist/+

  namelist /physicslist/ &
 Igeometry   ,& !latex \item \inputvar{Igeometry = 3} : \verb!integer! : selects Cartesian, cylindrical or toroidal geometry;
                !latex \bi
                !latex \item[i.] \inputvar{Igeometry = 1} : Cartesian; geometry determined by $R$;
                !latex \item[i.] \inputvar{Igeometry = 2} : cylindrical; geometry determined by $R$;
                !latex \item[i.] \inputvar{Igeometry = 3} : toroidal; geometry determined by $R$ {\em and} $Z$;
                !latex \ei
 Istellsym   ,& !latex \item \inputvar{Istellsym = 1} : \verb!integer! : stellarator symmetry is enforced if \inputvar{Istellsym.eq.1};
 Lfreebound  ,& !latex \item \inputvar{Lfreebound = 0} : \verb!integer! : compute vacuum field surrounding plasma;
 phiedge     ,& !latex \item \inputvar{phiedge = 1.0} : \verb!real! : total enclosed toroidal magnetic flux;
 curtor      ,& !latex \item \inputvar{curtor = 0.0} : \verb!real! : total enclosed (toroidal) plasma current;
 curpol      ,& !latex \item \inputvar{curpol = 0.0} : \verb!real! : total enclosed (poloidal) linking current;
 gamma       ,& !latex \item \inputvar{gamma = 0.0} : \verb!real! : adiabatic index; cannot set $|\gamma| = 1$;
 Nfp         ,& !latex \item \inputvar{Nfp = 1} : \verb!integer! : field periodicity;
                !latex \bi
                !latex \item[i.] all Fourier representations are of the form $\cos(m\t-nN\z)$, $\sin(m\t-nN\z)$,where $N\equiv$\inputvar{Nfp};
                !latex \item[i.] constraint : \inputvar{Nfp.ge.1};
                !latex \ei
 Nvol        ,& !latex \item \inputvar{Nvol = 1} : \verb!integer! : number of volumes;
                !latex \bi
                !latex \item[i.] each volume ${\cal V}_l$ is bounded by the ${\cal I}_{l-1}$ and ${\cal I}_{l}$ interfaces;
                !latex \item[i.] note that in cylindrical or toroidal geometry, ${\cal I}_{0}$ is the degenerate coordinate axis;
                !latex \item[i.] constraint : \inputvar{Nvol.le.MNvol};
                !latex \ei
 Mpol        ,& !latex \item \inputvar{Mpol = 1} : \verb!integer! : poloidal resolution;
 Ntor        ,& !latex \item \inputvar{Ntor = 0} : \verb!integer! : toroidal resolution;
                !latex \bi
                !l tex \item[i.] all Fourier representations of doubly-periodic functions are of the form
                !l tex \be f(\t,\z)& = &\sum_{n=0}^{\type{Ntor}} f_{0,n}\cos(-n \, \type{Nfp} \, \z)
                !l tex + \sum_{m=1}^{\type{Mpol}}\sum_{n=\type{-Ntor}}^{\type{Ntor}} f_{m,n}\cos(m\t-n \, \inputvar{Nfp} \, \z),
                !l tex \ee
                !latex Internally these ``double'' summations are written as a ``single'' summation,
                !latex e.g. $f = \sum_j f_j \cos(m_j\t-n_j\z)$.
                !latex \ei
 Lrad        ,& !latex \item \inputvar{Lrad = 4} : \verb!integer(MNvol+1)! : Chebyshev resolution in each volume;
                !latex \bi
                !latex \item[i.] constraint : \inputvar{Lrad(1:Mvol)}.ge.2;
                !latex \ei
 Lconstraint ,& !latex \item \inputvar{Lconstraint = -1} : \verb!integer! : selects constraints; primarily used in \link{ma02aa} and \link{mp00ac}.
                !latex \bi
                !latex \item[i.] if \inputvar{Lconstraint}.eq.-1, then in the plasma regions $\Delta\psi_t$, $\mu$ and $\Delta \psi_p$ are {\em not} varied;
                !latex           and in the vacuum region (only for free-boundary) $\Delta\psi_t$ and $\Delta \psi_p$ are {\em not} varied, and $\mu = 0$.
                !latex \item[ii.] if \inputvar{Lconstraint}.eq.0, then in the plasma regions $\Delta\psi_t$, $\mu$ and $\Delta \psi_p$ are {\em not} varied;
                !latex            and in the vacuum region (only for free-boundary) $\Delta\psi_t$ and $\Delta \psi_p$ are varied to match the
                !latex            prescribed plasma current, \inputvar{curtor}, and the ``linking'' current, \inputvar{curpol}, and $\mu = 0$;
                !latex \item[iii.] if \inputvar{Lconstraint}.eq.1, then in the plasma regions $\mu$ and $\Delta\psi_p$ are adjusted
                !latex             in order to satisfy the inner and outer interface transform constraints
                !latex             (except in the simple torus, where the enclosed poloidal flux is irrelevant,
                !latex             and only $\mu$ is varied to satisfy the outer interface transform constraint);
                !latex             and in the vacuum region $\Delta\psi_t$ and $\Delta \psi_p$ are varied to match the transform constraint on the boundary
                !latex             and to obtain the prescribed linking current, \inputvar{curpol}, and $\mu = 0$.
                !latex \item[iv.]  if \inputvar{Lconstraint}.eq.2, under reconstruction.
        !latex \item[v.] if \inputvar{Lconstraint} eq.3, then the $\mu$ and $\psi_p$ variables are adjusted in order to satisfy the volume and surface toroidal current
        !latex        computed with \link{lbpol} (excepted in the inner most volume, where the volume current is irrelevant). Not implemented yet in free
        !latex        boundary.
                !latex \ei
 tflux       ,& !latex \item \inputvar{tflux} : \verb!real(1:MNvol+1)! : toroidal flux, $\psi_t$, enclosed by each interface;
                !latex \bi
                !latex \item[i.] For each of the plasma volumes, this is a constraint: \inputvar{tflux} is \underline{not} varied;
                !latex \item[i.] For the vacuum region (only if \inputvar{Lfreebound = 1}), \inputvar{tflux} may be allowed to vary to match constraints;
                !latex \item[i.] Note that \inputvar{tflux} will be normalized so that \inputvar{tflux(Nvol) = 1.0},
                !latex so that \inputvar{tflux} is arbitrary up to a scale factor;
                !latex \item[i.] see also \inputvar{phiedge};
                !latex \ei
 pflux       ,& !latex \item \inputvar{pflux} : \verb!real(1:MNvol+1)! : poloidal flux, $\psi_p$, enclosed by each interface;
 helicity    ,& !latex \item \inputvar{helicity} : \verb!real(1:MNvol)! : helicity, ${\cal K}$, in each volume, ${\cal V}_i$;
                !latex \bi
                !latex \item[i.] on exit, \inputvar{helicity} is set to the computed values of ${\cal K} \equiv \int {\bf A}\cdot{\bf B}\;dv$;
                !latex \ei
 pscale      ,& !latex \item \inputvar{pscale = 0.0} : \verb!real! : pressure scale factor;
                !latex \bi
                !latex \item[i.] the initial pressure profile is given by \inputvar{pscale} $*$ \inputvar{press};
                !latex \ei
 pressure    ,& !latex \item \inputvar{pressure} : \verb!real(1:MNvol+1)! : pressure in each volume;
                !latex \bi
                !latex \item[i.] the pressure is {\em not} held constant, but $p_l V_l^\gamma = P_l$ {\em is} held constant,
                !latex where $P_l$ is determined by the initial pressures and the initial volumes, $V_l$;
                !latex \item[i.] (Note that if \inputvar{gamma = 0.0}, then $p_l \equiv P_l$.)
                !latex \item[i.] on output, the pressure is given by $p_l = P_l/V_l^\gamma$, where $V_l$ is the final volume;
                !latex \item[i.] \inputvar{pressure} is only used in calculation of interface force-balance;
                !latex \ei
 Ladiabatic  ,& !latex \item \inputvar{Ladiabatic = 0} : \verb!integer! : logical flag;
                !latex \bi
                !latex \item[i.] if \inputvar{Ladiabatic = 0}, the adiabatic constants are determined by the initial pressure and volume;
                !latex \item[i.] if \inputvar{Ladiabatic = 1}, the adiabatic constants are determined by the given input \inputvar{adiabatic};
                !latex \ei
 adiabatic   ,& !latex \item \inputvar{adiabatic} : \verb!real(1:MNvol+1)! : adiabatic constants in each volume;
                !latex \bi
                !latex \item[i.] the pressure is {\em not} held constant, but $p_l V_l^\gamma = P_l \equiv$\inputvar{adiabatic} is constant,
                !latex \item[i.] note that if \inputvar{gamma = 0.0}, then \inputvar{pressure = adiabatic};
                !latex \item[i.] \inputvar{pressure} is only used in calculation of interface force-balance;
                !latex \ei
 mu          ,& !latex \item \inputvar{mu} : \verb!real(1:MNvol+1)! : helicity-multiplier, $\mu$, in each volume;
 Ivolume     ,& !latex \item \inputvar{Ivolume} : \verb!real(1:MNvol+1)! : Toroidal current constraint normalized by $\mu_0$ ($I_{volume} = \mu_0\cdot [A]$), in each volume. This is a
                !latex          cumulative quantity: $I_{\mathcal{V},i} = \int_0^{\psi_{t,i}} \mathbf{J}\cdot\mathbf{dS}$. Physically, it represents the sum of all non-pressure driven currents;
 Isurf       ,& !latex \item \inputvar{Isurf} : \verb!real(1:MNvol)! : Toroidal current normalized by $\mu_0$ at each interface (cumulative). This is the sum of all pressure driven currents.;
 pl          ,& !latex \item \inputvar{pl = 0} : \verb!integer(0:MNvol)! :
 ql          ,& !latex \item \inputvar{ql = 0} : \verb!integer(0:MNvol)! :
 pr          ,& !latex \item \inputvar{pr = 0} : \verb!integer(0:MNvol)! :
 qr          ,& !latex \item \inputvar{qr = 0} : \verb!integer(0:MNvol)! :
                !latex \bi
                !latex \item[i.] ``inside'' interface rotational-transform is $\iotabar = (p_l+\gamma p_r)/(q_l+\gamma q_r)$,
                !latex        where $\gamma$ is the golden mean, $\gamma = (1 + \sqrt 5 ) / 2 $;
                !latex \item[i.] if both $q_l = 0$ {\em and} $q_r = 0$, then the (inside) interface rotational-transform is defined by \inputvar{iota};
                !latex \ei
 iota        ,& !latex \item \inputvar{iota} : \verb!real(0:MNvol)! : rotational-transform, $\iotabar$, on inner side of each interface;
                !latex \bi
                !latex \item[i.] only relevant if illogical input for \inputvar{ql} and \inputvar{qr} are provided;
                !latex \ei
 lp          ,& !latex \item \inputvar{lp = 0} : \verb!integer(0:MNvol)! :
 lq          ,& !latex \item \inputvar{lq = 0} : \verb!integer(0:MNvol)! :
 rp          ,& !latex \item \inputvar{rp = 0} : \verb!integer(0:MNvol)! :
 rq          ,& !latex \item \inputvar{rq = 0} : \verb!integer(0:MNvol)! :
                !latex \bi
                !latex \item "outer" interface rotational-transform is $\iotabar = (p_l+\gamma p_r)/(q_l+\gamma q_r)$,
                !latex       where $\gamma$ is the golden mean, $\gamma = (1 + \sqrt 5 ) / 2 $;
                !latex \item if both $q_l = 0$ {\em and} $q_r = 0$, then the (outer) interface rotational-transform is defined by \inputvar{oita};
                !latex \ei
 oita        ,& !latex \item \inputvar{oita} : \verb!real(0:MNvol)! : rotational-transform, $\iotabar$, on outer side of each interface;
                !latex \bi
                !latex \item only relevant if illogical input for \inputvar{ql} and \inputvar{qr} are provided;
                !latex \ei
 mupftol     ,& !latex \item \inputvar{mupftol = 1.0e-14} : \verb!real! : accuracy to which $\mu$ and $\Delta\psi_p$ are required;
                !latex \bi
                !latex \item only relevant if constraints on transform, enclosed currents etc. are to be satisfied iteratively, see \inputvar{Lconstraint};
                !latex \ei
 mupfits     ,& !latex \item \inputvar{mupfits = 8} : \verb!integer! : an upper limit on the transform/helicity constraint iterations;
                !latex \item only relevant if constraints on transform, enclosed currents etc. are to be satisfied iteratively, see \inputvar{Lconstraint};
                !latex \bi
                !latex \item constraint: \inputvar{mupfits > 0};
                !latex \ei
 rpol        ,& !latex \item \inputvar{rpol = 1.0} : \verb!real! : poloidal extent of slab (effective radius);
                !latex \bi
                !latex \item[i.] only relevant if \inputvar{Igeometry} $=1$;
                !latex \item[i.] poloidal size is $L = 2\pi*$\inputvar{rpol};
                !latex \ei
 rtor        ,& !latex \item \inputvar{rtor = 1.0} : \verb!real! : toroidal extent of slab (effective radius);
                !latex \bi
                !latex \item[i.] only relevant if \inputvar{Igeometry} $=1$;
                !latex \item[i.] toroidal size is $L = 2\pi*$\inputvar{rtor};
                !latex \ei
 Rac         ,& !latex \item \inputvar{Rac} : \verb!real(     0:MNtor             )! : Fourier harmonics of axis    ;     stellarator symmetric;
 Zas         ,& !latex \item \inputvar{Zas} : \verb!real(     0:MNtor             )! : Fourier harmonics of axis    ;     stellarator symmetric;
 Ras         ,& !latex \item \inputvar{Ras} : \verb!real(     0:MNtor             )! : Fourier harmonics of axis    ; non-stellarator symmetric;
 Zac         ,& !latex \item \inputvar{Zac} : \verb!real(     0:MNtor             )! : Fourier harmonics of axis    ; non-stellarator symmetric;
 Rbc         ,& !latex \item \inputvar{Rbc} : \verb!real(-MNtor:MNtor,-MMpol:MMpol)! : Fourier harmonics of boundary;     stellarator symmetric;
 Zbs         ,& !latex \item \inputvar{Zbs} : \verb!real(-MNtor:MNtor,-MMpol:MMpol)! : Fourier harmonics of boundary;     stellarator symmetric;
 Rbs         ,& !latex \item \inputvar{Rbs} : \verb!real(-MNtor:MNtor,-MMpol:MMpol)! : Fourier harmonics of boundary; non-stellarator symmetric;
 Zbc         ,& !latex \item \inputvar{Zbc} : \verb!real(-MNtor:MNtor,-MMpol:MMpol)! : Fourier harmonics of boundary; non-stellarator symmetric;
 Rwc         ,& !latex \item \inputvar{Rwc} : \verb!real(-MNtor:MNtor,-MMpol:MMpol)! : Fourier harmonics of wall    ;     stellarator symmetric;
 Zws         ,& !latex \item \inputvar{Zws} : \verb!real(-MNtor:MNtor,-MMpol:MMpol)! : Fourier harmonics of wall    ;     stellarator symmetric;
 Rws         ,& !latex \item \inputvar{Rws} : \verb!real(-MNtor:MNtor,-MMpol:MMpol)! : Fourier harmonics of wall    ; non-stellarator symmetric;
 Zwc         ,& !latex \item \inputvar{Zwc} : \verb!real(-MNtor:MNtor,-MMpol:MMpol)! : Fourier harmonics of wall    ; non-stellarator symmetric;
 Vns         ,& !latex \item \inputvar{Vns} : \verb!real(-MNtor:MNtor,-MMpol:MMpol)! : Fourier harmonics of vacuum normal field at boundary;
 Bns         ,& !latex \item \inputvar{Bns} : \verb!real(-MNtor:MNtor,-MMpol:MMpol)! : Fourier harmonics of plasma normal field at boundary;
 Vnc         ,& !latex \item \inputvar{Vnc} : \verb!real(-MNtor:MNtor,-MMpol:MMpol)! : Fourier harmonics of vacuum normal field at boundary;
 Bnc         ,& !latex \item \inputvar{Bnc} : \verb!real(-MNtor:MNtor,-MMpol:MMpol)! : Fourier harmonics of plasma normal field at boundary;
 Lreflect       !latex \item \inputvar{Lreflect = 0} : \verb!integer! : =1 reflect the upper and lower bound in slab, =0 do not reflect
!latex \ei

!latex \end{enumerate}

!latex \subsubsection{\type{numericlist} :}

!latex \begin{enumerate}
!latex \item The namelist \type{numericlist} controls internal resolution parameters that the user rarely needs to consider.
!latex \bi
!latex \type{namelist/numericlist/}

  namelist /numericlist/ &
 Linitialize ,& !latex \item \inputvar{Linitialize = 0 : integer} : to initialize geometry using a regularization / extrapolation method;
                !latex \bi
                !latex \item if \inputvar{Linitialize = -I}, where $I$ is a positive integer,
                !latex       the geometry of the $i=1,N_V-I$ surfaces constructed by an extrapolation;
                !latex \item if \inputvar{Linitialize = 0}, the geometry of the interior surfaces is provided after the namelists in the input file;
                !latex \item if \inputvar{Linitialize = 1}, the interior surfaces will be intialized as $R_{l,m,n} = R_{N,m,n} \psi_{t,l}^{m/2}$,
                !latex where $R_{N,m,n}$ is the plasma boundary
                !latex       and $\psi_{t,l}$ is the given toroidal flux enclosed by the $l$-th interface, normalized to the total enclosed toroidal flux;
                !latex       a similar extrapolation is used for $Z_{l,m,n}$;
                !latex \item note that the Fourier harmonics of the boundary is {\em always} given by the \inputvar{Rbc} and \inputvar{Zbs}
                !latex given in \type{physicslist};
                !latex \item if \inputvar{Linitialize = 2}, the interior surfaces {\em and the plasma boundary} will be intialized
                !latex       as $R_{l,m,n} = R_{W,m,n} \psi_{t,l}^{m/2}$, where $R_{W,m,n}$ is the computational boundary
                !latex       and $\psi_{t,l}$ is the given toroidal flux enclosed by the $l$-th interface, normalized to the total enclosed toroidal flux;
                !latex       a similar extrapolation is used for $Z_{l,m,n}$;
                !latex \item note that, for free-boundary calculations, the Fourier harmonics of the computational boundary
                !latex       is {\em always} given by the \inputvar{Rwc} and \inputvar{Zws}
                !latex given in \type{physicslist};
                !latex \item if \inputvar{Linitialize = 1, 2}, it is not required to provide the geometry of the interfaces after the namelists;
                !latex \ei
 LautoinitBn ,& !latex \item \inputvar{LautoinitBn = 1 : integer} : to initialize $B_{ns}$ using an initial fixed-boundary calculation;
                !latex \bi
                !latex \item only relevant if \inputvar{Lfreebound = 1},
                !latex \item user-supplied \inputvar{Bns} will only be considered if \inputvar{LautoinitBn = 0}
                !latex \ei
 Lzerovac    ,& !latex \item \inputvar{Lzerovac = 0 : integer} : to adjust vacuum field to cancel plasma field on computational boundary;
                !latex \bi
                !latex \item only relevant if \inputvar{Lfreebound = 1},
                !latex \ei
 Ndiscrete   ,& !latex \item \inputvar{Ndiscrete = 2 : integer} :
                !latex \bi
                !latex \item resolution of the real space grid on which fast Fourier transforms are performed is given by \inputvar{Ndiscrete*Mpol*4};
                !latex \item constraint \inputvar{Ndiscrete>0};
                !latex \ei
 Nquad       ,& !latex \item \inputvar{Nquad = -1 : integer} : the resolution of the Gaussian quadrature;
                !latex \bi
                !latex \item the resolution of the Gaussian quadrature, $\ds \int \!\! f(s) ds = \sum_k \omega_k f(s_k)$,
                !latex       in each volume is given by \internal{Iquad$_v$},
                !latex \item \internal{Iquad$_v$} is set in \link{preset}.
                !l tex       and depends on \inputvar{Nquad}, \inputvar{Lrad$_v$} and \inputvar{Mpol}.
                !l tex \bi
                !l tex \item if \inputvar{Nquad.gt.0},                                 then \internal{Iquad(vvol) =              Nquad};
                !l tex \item if \inputvar{Nquad.le.0 and .not.Lcoordinatesingularity}, then \internal{Iquad(vvol) = 2*Lrad(vvol)-Nquad};
                !l tex \item if \inputvar{Nquad.le.0 and      Lcoordinatesingularity}, then \internal{Iquad(vvol) = 2*Lrad(vvol)-Nquad+Mpol};
                !l tex \ei
                !l tex \item \internal{Iquad$_v$} is passed through to \link{ma00aa} to compute various volume integrals;
                !l tex       also see \link{jo00aa}, where \internal{Iquad$_v$}
                !l tex       is also used in computing the volume integrals of $||\nabla\times{\bf B} - \mu {\bf B}||$;
                !latex \ei
 iMpol       ,& !latex \item \inputvar{iMpol = -4 : integer} : Fourier resolution of straight-fieldline angle on interfaces;
                !latex \bi
                !latex \item the rotational-transform on the interfaces is determined by a transformation to the straight-fieldline angle,
                !latex with poloidal resolution given by \inputvar{iMpol};
                !latex \item if \inputvar{iMpol.le.0}, then \inputvar{iMpol = Mpol - iMpol};
                !latex \ei
 iNtor       ,& !latex \item \inputvar{iNtor = -4 : integer} : Fourier resolution of straight-fieldline angle on interfaces;
                !latex \bi
                !latex \item the rotational-transform on the interfaces is determined by a transformation to the straight-fieldline angle,
                !latex with toroidal resolution given by \inputvar{iNtor};
                !latex \item if \inputvar{iNtor.le.0}, then \inputvar{iNtor = Ntor - iNtor};
                !latex \item if \inputvar{Ntor.eq.0}, then the toroidal resolution of the angle transformation is set \inputvar{lNtor = 0}.
                !latex \ei
 Lsparse     ,& !latex \item \inputvar{Lsparse = 0 : integer} : controls method used to solve for rotational-transform on interfaces;
                !latex \bi
                !latex \item if \inputvar{Lsparse = 0}, the transformation to the straight-fieldline angle is computed in Fourier space
                !latex using a dense matrix solver, \nag{}{F04AAF};
                !latex \item if \inputvar{Lsparse = 1}, the transformation to the straight-fieldline angle is computed in real space
                !latex using a dense matrix solver, \nag{}{F04ATF};
                !latex \item if \inputvar{Lsparse = 2}, the transformation to the straight-fieldline angle is computed in real space
                !latex using a sparse matrix solver, \nag{}{F11DEF};
                !latex \item if \inputvar{Lsparse = 3}, the different methods for constructing the straight-fieldline angle are compared;
                !latex \ei
 Lsvdiota    ,& !latex \item \inputvar{Lsvdiota = 0 : integer} : controls method used to solve for rotational-transform on interfaces;
                !latex only relevant if \inputvar{Lsparse = 0};
                !latex \bi
                !latex \item if \inputvar{Lsvdiota = 0}, use standard linear solver to construct straight fieldline angle transformation;
                !latex \item if \inputvar{Lsvdiota = 1}, use SVD method to compute rotational-transform;
                !latex \ei
 imethod     ,& !latex \item \inputvar{Imethod = 3 : integer} : controls iterative solution to sparse matrix
                !latex arising in real-space transformation to the straight-fieldline angle;
                !latex only relevant if \inputvar{Lsparse.eq.2}; see \link{tr00ab} for details;
                !latex \bi
                !latex \item if \inputvar{imethod = 1}, the method is \type{RGMRES};
                !latex \item if \inputvar{imethod = 2}, the method is \type{CGS};
                !latex \item if \inputvar{imethod = 3}, the method is \type{BICGSTAB};
                !latex \ei
 iorder      ,& !latex \item \inputvar{iorder = 2 : integer} : controls real-space grid resolution for constructing the straight-fieldline angle;
                !latex only relevant if \inputvar{Lsparse>0};
                !latex determines order of finite-difference approximation to the derivatives;
                !latex \bi
                !latex \item if \inputvar{iorder = 2},
                !latex \item if \inputvar{iorder = 4},
                !latex \item if \inputvar{iorder = 6},
                !latex \ei
 iprecon     ,& !latex \item \inputvar{Iprecon = 0 : integer} : controls iterative solution to sparse matrix arising in real-space transformation
                !latex to the straight-fieldline angle;
                !latex only relevant if \inputvar{Lsparse.eq.2}; see \link{tr00ab} for details;
                !latex \bi
                !latex \item if \inputvar{iprecon = 0}, the preconditioner is `N';
                !latex \item if \inputvar{iprecon = 1}, the preconditioner is `J';
                !latex \item if \inputvar{iprecon = 2}, the preconditioner is `S';
                !latex \ei
 iotatol     ,& !latex \item \inputvar{iotatol = -1.0 : real} : tolerance required for iterative construction of straight-fieldline angle;
                !latex only relevant if \inputvar{Lsparse.ge.2};
 Lextrap     ,& !latex \item \inputvar{Lextrap = 0 : integer} : geometry of innermost interface is defined by extrapolation;
 Mregular    ,& !latex \item \inputvar{Mregular = -1 : integer} : maximum regularization factor;                !latex \bi
                !latex \bi
                !latex \item if \inputvar{Mregular.ge.2}, then \internal{regumm}$_i$ = \inputvar{Mregular} $/ 2 $ where \internal{m}$_i > $ \inputvar{Mregular}
                !latex \ei
 Lrzaxis     ,& !latex \item \inputvar{Lrzaxis = 1 : integer} : controls the guess of geometry axis in the innermost volume or initialization of interfaces
                !latex \bi
                !latex \item if \inputvar{iprecon = 1}, the centroid is used;
                !latex \item if \inputvar{iprecon = 2}, the Jacobian $m=1$ harmonic elimination method is used;
                !latex \ei
 Ntoraxis       !latex \item \inputvar{Ntoraxis = 3 : integer} : the number of $n$ harmonics used in the Jacobian $m=1$ harmonic elimination method;
                !latex only relevant if \inputvar{Lrzaxis.ge.1};

!latex \ei

!latex \end{enumerate}

!latex \subsubsection{\type{locallist} : }

!latex \begin{enumerate}
!latex \item The namelist \type{locallist} controls the construction of the Beltrami fields in each volume.
!latex \bi
!latex \type{namelist/locallist/}

  namelist /locallist/ &
 LBeltrami,&    !latex \item\inputvar{LBeltrami = 4  integer}
                !latex \bi
                !latex \item if \inputvar{LBeltrami = 1,3,5 or 7}, (SQP) then the Beltrami field in each volume is constructed
                !latex by minimizing the magnetic energy with the constraint of fixed helicity;
                !latex this is achieved by using sequential quadratic programming as provided by \nag{}{E04UFF};
                !latex this approach has the benefit (in theory) of robustly constructing minimum energy solutions
                !latex when multiple, i.e. bifurcated, solutions exist.
                !latex \item if \inputvar{LBeltrami = 2,3,6 or 7}, (Newton) then the Beltrami fields are constructed by employing a standard Newton method
                !latex for locating an extremum of
                !latex $F\equiv \int B^2 dv - \mu (\int {\bf A}\cdot{\bf B}dv-{\cal K})$,
                !latex where $\mu$ is treated as an independent degree of freedom similar to the parameters describing the vector potential
                !latex and ${\cal K}$ is the required value of the helicity;
                !latex this is the standard Lagrange multipler approach for locating the constrained minimum;
                !latex this method cannot distinguish saddle-type extrema from minima, and which solution that will be obtained depends on the initial guess;
                !latex \item if \inputvar{LBeltrami = 4,5,6 or 7}, (linear) it is assumed that the Beltrami fields are parameterized by $\mu$;
                !latex in this case, it is only required to solve $\nabla \times {\bf B} = \mu {\bf B}$ which reduces to a system of linear equations;
                !latex $\mu$ may or may not be adjusted iteratively, depending on \inputvar{Lconstraint},
                !latex to satisfy either rotational-transform or helicity constraints;
                !latex \item for flexibility and comparison, each of the above methods can be employed; for example:
                !latex \bi
                !latex \item if \inputvar{LBeltrami = 1}, only the SQP    method will be employed;
                !latex \item if \inputvar{LBeltrami = 2}, only the Newton method will be employed;
                !latex \item if \inputvar{LBeltrami = 4}, only the linear method will be employed;
                !latex \item if \inputvar{LBeltrami = 3}, the SQP and the Newton method are used;
                !latex \item if \inputvar{LBeltrami = 5}, the SQP and the linear method are used;
                !latex \item if \inputvar{LBeltrami = 6}, the Newton and the linear method are used;
                !latex \item if \inputvar{LBeltrami = 7}, all three methods will be employed;
                !latex \ei
                !latex \ei
 Linitgues,&    !latex \item\inputvar{Linitgues = 1  integer} controls how initial guess for Beltrami field is constructed;
                !latex \bi
                !latex \item only relevant for routines that require an initial guess for the Beltrami fields, such as the SQP and Newton methods,
                !latex or the sparse linear solver;
                !latex \item if \inputvar{Linitgues = 0}, the initial guess for the Beltrami field is trivial;
                !latex \item if \inputvar{Linitgues = 1}, the initial guess for the Beltrami field is an integrable approximation;
                !latex \item if \inputvar{Linitgues = 2}, the initial guess for the Beltrami field is read from file;
                !latex \item if \inputvar{Linitgues = 3}, the initial guess for the Beltrami field will be randomized with the maximum \inputvar{maxrndgues};
                !latex \ei
 maxrndgues,&   !latex \item \inputvar{maxrndgues = 1.0} : real : the maximum random number of the Beltrami field if \inputvar{Linitgues = 3};
 Lmatsolver,&!latex \item \inputvar{Lmatsolver = 3} : integer : 1 for LU factorization, 2 for GMRES, 3 for GMRES matrix-free
 NiterGMRES,&   !latex \item \inputvar{LGMRESniter = 200} : integer : number of max iteration for GMRES
 epsGMRES,&     !latex \item \inputvar{epsGMRES = 1e-14} : real : the precision of GMRES
 LGMRESprec,&!latex \item \inputvar{LGMRESprec = 1} : integer : type of preconditioner for GMRES, 1 for ILU sparse matrix
 epsILU,&       !latex \item \inputvar{epsILU = 1e-12} : real : the precision of incomplete LU factorization for preconditioning
 Lposdef        !latex \item\inputvar{Lposdef = 0 : integer} : redundant;
!Nmaxexp        !l tex \item \inputvar{Nmaxexp = 32 : integer} : indicates maximum exponent used to precondition Beltrami linear system near singularity;
                !l tex \bi
                !l tex \item a factor of $s^{\bar m_j/2}$ is extracted, where $\bar m \equiv min(\inputvar{Nmaxexp},m_j)$;
                !l tex \ei
!latex \ei

!latex \item Comments:
!latex \begin{enumerate}
!latex \item The transformation to straight-fieldline coordinates is singular when the rotational-transform of the interfaces is rational;
!latex       however, the rotational-transform is still well defined.
!latex \end{enumerate}

!latex \end{enumerate}

!latex \hrule

!latex \subsubsection{\type{globallist} : }

!latex \begin{enumerate}
!latex \item The namelist \type{globallist} controls the search for global force-balance:
!latex \bi
!latex \type{namelist/globallist/}

  namelist /globallist/ &
 Lfindzero   ,& !latex \item \inputvar{Lfindzero = 0} : integer : use Newton methods to find zero of force-balance, which is computed by \link{dforce};
                !latex \bi
                !latex \item[o.] if \inputvar{Lfindzero = 0}, then \link{dforce} is called once
                !latex           to compute the Beltrami fields consistent with the given geometry and constraints;
                !latex \item[i.] if \inputvar{Lfindzero = 1}, then call
                !latex           \nag{www.nag.co.uk/numeric/FL/manual19/pdf/C05/c05ndf_fl19.pdf}{C05NDF} (uses function values only),
                !latex           which iteratively calls \link{dforce};
                !latex \item[ii.] if \inputvar{Lfindzero = 2}, then call
                !latex            \nag{www.nag.co.uk/numeric/FL/manual19/pdf/C05/c05pdf_fl19.pdf}{C05PDF} (uses derivative information),
                !latex            which iteratively calls \link{dforce};
                !latex \ei
 escale      ,& !latex \item \inputvar{escale = 0.0} : real : controls the weight factor, \type{BBweight}, in the force-imbalance harmonics;
                !latex \bi
                !latex \item[i.] \type{BBweight(i)} $\ds \equiv \inputvar{opsilon} \times \exp\left[-\inputvar{escale} \times (m_i^2+n_i^2) \right]$
                !latex \item[ii.] defined in \link{preset}; used in \link{dforce};
                !latex \item[iii.] also see \Eqn{forcebalancemn} below;
                !latex \ei
 opsilon     ,& !latex \item \inputvar{opsilon = 1.0} : real : weighting of force-imbalance;
                !latex \bi
                !latex \item[i.] used in \link{dforce}; also see \Eqn{forcebalancemn} below;
                !latex \ei
 pcondense   ,& !latex \item \inputvar{pcondense = 2.0} : real : spectral condensation parameter;
                !latex \bi
                !latex \item[i.] used in \link{preset} to define \type{mmpp(i)} $\equiv m_i^p$, where $p\equiv $ \inputvar{pcondense};
                !latex \item[ii.] the angle freedom is exploited to minimize $\ds \inputvar{epsilon} \sum_{i} m_i^p (R_{i}^2+Z_{i}^2)$
                !latex       with respect to tangential variations in the interface geometry;
                !latex \item[ii.] also see \Eqn{spectralbalancemn} below;
                !latex \ei
 epsilon     ,& !latex \item \inputvar{epsilon = 0.0} : real : weighting of spectral-width constraint;
                !latex \bi
                !latex \item[i.] used in \link{dforce}; also see \Eqn{spectralbalancemn} below;
                !latex \ei
 wpoloidal   ,& !latex \item \inputvar{wpoloidal = 1.0} : real : ``star-like'' poloidal angle constraint radial exponential factor;
                !latex       used in \link{preset} to construct \type{sweight}
 upsilon     ,& !latex \item \inputvar{upsilon = 1.0} : real : weighting of ``star-like'' poloidal angle constraint;
                !latex       used in \link{preset} to construct \type{sweight};
 forcetol    ,& !latex \item \inputvar{forcetol = 1.0e-10} : real : required tolerance in force-balance error; only used as an initial check;
                !latex \bi
                !latex \item[i.] if the initially supplied interfaces are consistent with force-balance to within \inputvar{forcetol},
                !latex       then the geometry of the interfaces is not altered;
                !latex \item[ii.] if not, then the geometry of the interfaces is changed in order to bring the configuration into forcebalance
                !latex       so that the geometry of interfaces is within \inputvar{c05xtol}, defined below, of the true solution;
                !latex \item[iii.] to force execution of either \nag{www.nag.co.uk/numeric/FL/manual19/pdf/C05/c05ndf_fl19.pdf}{C05NDF}
                !latex       or \nag{www.nag.co.uk/numeric/FL/manual19/pdf/C05/c05pdf_fl19.pdf}{C05PDF}, regardless of the initial force imbalance,
                !latex       set \inputvar{forcetol < 0};
                !latex \ei
 c05xmax     ,& !latex \item \inputvar{c05xmax = 1.0e-06} : real : required tolerance in position, ${\bf x} \equiv \{ R_{i,v}, Z_{i,v}\}$;
 c05xtol     ,& !latex \item \inputvar{c05xtol = 1.0e-12} : real : required tolerance in position, ${\bf x} \equiv \{ R_{i,v}, Z_{i,v}\}$;
                !latex \bi
                !latex \item[i.] used by both \nag{www.nag.co.uk/numeric/FL/manual19/pdf/C05/c05ndf_fl19.pdf}{C05NDF} and
                !latex           \nag{www.nag.co.uk/numeric/FL/manual19/pdf/C05/c05pdf_fl19.pdf}{C05PDF};
                !latex           see the NAG documents for further details on how the error is defined;
                !latex \item[ii.] constraint \inputvar{c05xtol.gt.0.0};
                !latex \ei
 c05factor   ,& !latex \item \inputvar{c05factor = 1.0e-02} : real : used to control initial step size in
                !latex       \nag{www.nag.co.uk/numeric/FL/manual19/pdf/C05/c05ndf_fl19.pdf}{C05NDF} and
                !latex       \nag{www.nag.co.uk/numeric/FL/manual19/pdf/C05/c05pdf_fl19.pdf}{C05PDF};
                !latex \bi
                !latex \item[i.] constraint \inputvar{c05factor.gt.0.0};
                !latex \item[ii.] only relevant if \inputvar{Lfindzero.gt.0};
                !latex \ei
 LreadGF     ,& !latex \item \inputvar{LreadGF = T} : logical : read $\nabla_{\bf x} {\bf F}$ from file \type{.GF};
                !latex \bi
                !latex \item[i.] only used if \inputvar{Lfindzero = 2};
                !latex \item[ii.] only used in \link{newton};
                !latex \ei
 mfreeits    ,& !latex \item \inputvar{mfreeits = 0} : integer : maximum allowed free-boundary iterations;
                !latex \bi
                !latex \item[i.] only used if \inputvar{Lfreebound = 1};
                !latex \item[ii.] only used in \link{xspech};
                !latex \ei
 bnstol      ,& !latex \item \inputvar{bnstol = 1.0e-06} : redundant;
 bnsblend    ,& !latex \item \inputvar{bnsblend = 0.666} : redundant;
 gBntol      ,& !latex \item \inputvar{gBntol = 1.0e-06} : real : required tolerance in free-boundary iterations;
                !latex \bi
                !latex \item[i.] only used if \inputvar{Lfreebound = 1};
                !latex \item[ii.] only used in \link{xspech}; see \link{xspech} for more documentation;
                !latex \ei
 gBnbld      ,& !latex \item \inputvar{gBnbld = 0.666} : real : normal blend;
                !latex \bi
                !latex \item[i.] The ``new'' magnetic field at the computational boundary produced by the plasma currents is updated using a Picard scheme:
                !latex           \be ({\bf B}\cdot{\bf n})^{j+1} =    \inputvar{gBnbld}  \times ({\bf B}\cdot{\bf n})^{j}
                !latex                                           + (1-\inputvar{gBnbld}) \times ({\bf B}\cdot{\bf n})^{*},
                !latex           \ee
                !latex           where $j$ labels free-boundary iterations, and $({\bf B}\cdot{\bf n})^{*}$ is computed by virtual casing.
                !latex \item[ii.] only used if \inputvar{Lfreebound = 1};
                !latex \item[ii.] only used in \link{xspech};
                !latex \ei
 vcasingeps  ,& !latex \item \inputvar{vcasingeps = 1.0e-12} : real : regularization of Biot-Savart; see \link{bnorml}, \link{casing};
 vcasingtol  ,& !latex \item \inputvar{vcasingtol = 1.0e-08} : real : accuracy on virtual casing integral; see \link{bnorml}, \link{casing};
 vcasingits  ,& !latex \item \inputvar{vcasingits = 8      } : integer : minimum number of calls to adaptive virtual casing routine; see \link{casing};
 vcasingper  ,& !latex \item \inputvar{vcasingper = 1      } : integer : periods of integragion  in adaptive virtual casing routine; see \link{casing};
 mcasingcal     !latex \item \inputvar{mcasingcal = 8      } : integer : minimum number of calls to adaptive virtual casing routine; see \link{casing};
!latex \ei

!latex \item Comments:
!latex \begin{enumerate}
!latex \item The ``force'' vector, ${\bf F}$, which is constructed in \link{dforce}, is a combination of pressure-imbalance Fourier harmonics,
!latex       \be F_{i,v} \equiv [[ p+B^2/2 ]]_{i,v} \times \exp\left[-\inputvar{escale}(m_i^2+n_i^2) \right] \times \inputvar{opsilon},
!latex       \label{eq:forcebalancemn} \ee
!latex       and spectral-condensation constraints, $I_{i,v}$, and the ``star-like'' angle constraints, $S_{i,v,}$, (see \link{lforce} for details)
!latex       \be F_{i,v} \equiv \inputvar{epsilon} \times I_{i,v}
!latex           + \inputvar{upsilon} \times \left( \psi_v^\omega S_{i,v,1} - \psi_{v+1}^\omega S_{i,v+1,0} \right),
!latex       \label{eq:spectralbalancemn} \ee
!latex       where $\psi_v\equiv$ normalized toroidal flux, \inputvar{tflux}, and $\omega\equiv $ \inputvar{wpoloidal}.
!latex \end{enumerate}

!latex \end{enumerate}

!latex \hrule

!latex \subsubsection{\type{diagnosticslist} : }

!latex \begin{enumerate}
!latex \item The namelist \type{diagnosticslist} controls post-processor diagnostics, such as \Poincare plot resolution, $\dots$,
!latex \bi
!latex \type{namelist/diagnosticslist/}

  namelist /diagnosticslist/ &
 odetol     ,&  !latex \item \inputvar{odetol = 1.0e-07} : real : o.d.e. integration tolerance for all field line tracing routines;
 absreq     ,&  !latex \item \inputvar{absreq = 1.0e-08} : real : redundant;
 relreq     ,&  !latex \item \inputvar{relreq = 1.0e-08} : real : redundant;
 absacc     ,&  !latex \item \inputvar{absacc = 1.0e-04} : real : redundant;
 epsr       ,&  !latex \item \inputvar{epsr = 1.0e-06} : real : redundant;
 nPpts      ,&  !latex \item \inputvar{nPpts = 0} : integer : number of toroidal transits used (per trajectory) in following field lines
                !latex for constructing \Poincare plots;
                !latex if \inputvar{nPpts<1}, no \Poincare plot is constructed;
 Ppts       ,&  !latex \item \inputvar{Ppts = 0} : stands for Poincare plot theta start. Chose at which angle (normalized over $\pi$) the Poincare field-line
                !latex tracing start.
 nPtrj      ,&  !latex \item \inputvar{nPtrj = -1 : integer(1:MNvol+1)} : number of trajectories in each annulus to be followed in constructing \Poincare plot;
                !latex \bi
                !latex \item if \inputvar{nPtrj(l)<0}, then \inputvar{nPtrj(l) = Ni(l)},
                !latex where \type{Ni(l)} is the grid resolution used to construct the Beltrami field in volume $l$;
                !latex \ei
 LHevalues  ,&  !latex \item \inputvar{LHevalues = F} : logical : to compute eigenvalues of $\nabla {\bf F}$;
 LHevectors ,&  !latex \item \inputvar{LHevectors = F} : logical : to compute eigenvectors (and also eigenvalues) of $\nabla {\bf F}$;
 LHmatrix   ,&  !latex \item \inputvar{LHmatrix = F} : logical : to compute and write to file the elements of $\nabla {\bf F}$;
 Lperturbed ,&  !latex \item \inputvar{Lperturbed = 0} : integer : to compute linear, perturbed equilibrium;
 dpp        ,&  !latex \item \inputvar{dpp = 1} : integer : perturbed harmonic;
 dqq        ,&  !latex \item \inputvar{dqq = 1} : integer : perturbed harmonic;
 Lerrortype ,&  !latex \item \inputvar{Lerrortype = 0} : integer : the type of error output for Lcheck=1
 Ngrid      ,&  !latex \item \inputvar{Ngrid=-1} : integer : the number of points to output in the grid, -1 for Lrad(vvol)
 Lcheck     ,&  !latex \item \inputvar{Lcheck = 0} : integer : implement various checks;
                !latex \bi
                !latex \item if \inputvar{Lcheck = 0}, no additional check on the calculation is performed;
                !latex \item if \inputvar{Lcheck = 1}, the error in the current, i.e. $\nabla\times{\bf B}-\mu{\bf B}$ is computed as a post-diagnostic;
                !latex \item if \inputvar{Lcheck = 2}, the analytic derivatives of the interface transform w.r.t.
                !latex       the helicity multiplier, $\mu$, and the enclosed poloidal flux, $\Delta\psi_p$, are compared to a finite-difference estimate;
                !latex \bi
                !latex \item[i.] only if \inputvar{Lconstraint.eq.1};
                !latex \item[ii.] only for \type{dspec} executable, i.e. must compile with \type{DFLAGS = "-D DEBUG"};
                !latex \ei
                !latex \item if \inputvar{Lcheck = 3}, the analytic derivatives of the volume w.r.t. interface Fourier harmonic
                !latex       is compared to a finite-difference estimate;
                !latex \bi
                !latex \item[i.] must set \inputvar{Lfindzero}$ = 2$,
                !latex \item[ii.] set \inputvar{forcetol} sufficiently small and set \inputvar{LreadGF = F},
                !latex       so that the matrix of second derivatives is calculated,
                !latex \item[iii.] only for \type{dspec} executable, i.e. must compile with \type{DFLAGS = "-D DEBUG"};
                !latex \ei
                !latex \item if \inputvar{Lcheck = 4}, the analytic calculation of the derivatives of the magnetic field, $B^2$, at the interfaces
                !latex       is compared to a finite-difference estimate;
                !latex \bi
                !latex \item[i.] must set \inputvar{Lfindzero}$ = 2$,
                !latex \item[ii.] set \inputvar{forcetol} sufficiently small,
                !latex \item[iii.] set \inputvar{LreadGF=F},
                !latex \item[iv.] only for \type{dspec} executable, i.e. must compile with \type{DFLAGS = "-D DEBUG"};
                !latex \ei
                !latex \item if \inputvar{Lcheck = 5}, the analytic calculation of the matrix of the derivatives of the force imbalance
                !latex       is compared to a finite-difference estimate;
                !latex \item if \inputvar{Lcheck = 6}, the virtual casing calculation is compared to \verb+xdiagno+;
                !latex \bi
                !latex \item[i.] the input file for \verb+xdiagno+ is written by \link{bnorml};
                !latex \item[ii.] this provides the Cartesian coordinates on the computational boundary where the virtual casing routine \link{casing}
                !latex           computes the magnetic field, with the values of the magnetic field being written to the screen for comparison;
                !latex \item[iii.] must set \inputvar{Freebound=1}, \inputvar{Lfindzero.gt.0}, \inputvar{mfreeits.ne.0};
                !latex \item[iii.] \verb+xdiagno+; must be executed manually;
                !latex \ei
                !latex \ei
 dRZ        ,&  !latex \item \inputvar{dRZ = 1E-5} : real; difference in geometry for finite difference estimate (debug only)
 Ltiming    ,&  !latex \item \inputvar{Ltiming = T} : logical : to check timing;
 fudge      ,&  !latex \item \inputvar{fudge = 1.0} : real : redundant;
 scaling        !latex \item \inputvar{scaling = 1.0} : real : redundant;

!latex \ei

!latex \end{enumerate}

!latex \subsubsection{\type{screenlist} : }

!latex \begin{enumerate}
!latex \item The namelist \type{screenlist} controls screen output.
!latex \bi
!latex \type{namelist/screenlist/}

  namelist /screenlist/ &
! NSCREENLIST ! namelist screenlist; this is expanded by Makefile; DO NOT REMOVE;
 Wbuild_vector_potential , &
 Wreadin , &  !latex \item Every subroutine, e.g. \type{xy00aa.h}, has its own write flag, \type{Wxy00aa}.
 Wwritin , & ! redundant;
 Wwrtend , &
 Wmacros

!latex \ei

!latex \end{enumerate}

!latex \subsection{input geometry}

!latex \begin{enumerate}
!latex \item The geometry of the $l$-th interface, for $l=0,N$ where $N\equiv$ \inputvar{Nvol}, is described by a set of Fourier harmonics,
!latex using an arbitrary poloidal angle,
!latex \be R_l(\t,\z)&=&\sum_{j}R_{j,l}\cos(m_j\t-n_j\z), \\ Z_l(\t,\z)&=&\sum_{j}Z_{j,l}\sin(m_j\t-n_j\z). \ee
!latex \item These harmonics are read from the \type{ext.sp} file and come directly after the namelists described above.
!latex The required format is as follows:
!latex \be \begin{array}{ccccccccc}
!latex m_1 & n_1 & R_{1,0} & Z_{1,0} & R_{1,1} & Z_{1,1} & \dots & R_{1,N} & Z_{1,N} \\
!latex m_2 & n_2 & R_{2,0} & Z_{2,0} & R_{2,1} & Z_{2,1} & \dots & R_{2,N} & Z_{2,N} \\
!latex \dots \\
!latex m_j & n_j & R_{j,0} & Z_{j,0} & R_{j,1} & Z_{j,1} & \dots & R_{j,N} & Z_{j,N} \\
!latex \dots
!latex \end{array}
!latex \ee
!latex
!latex \item The coordinate axis corresponds to $j=0$ and the outermost boundary corresponds to $j=$\inputvar{Nvol}.
!latex \item An arbitrary selection of harmonics may be inluded in any order, but only those within the range specified by \inputvar{Mpol} and \inputvar{Ntor}
!latex will be used.
!latex \item The geometry of {\em all} the interfaces, i.e. $l=0,N$, including the degenerate `coordinate-axis' interface, must be given.
!latex \end{enumerate}

! note that all variables in namelist need to be broadcasted in readin;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

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

subroutine read_inputlists_from_file()

   use fileunits
   use allglobal, only: ext, cpus, myid, MPI_COMM_SPEC

   LOCALS

   LOGICAL              :: Lspexist

   inquire( file=trim(ext)//".sp", exist=Lspexist ) ! check if file exists;
   FATAL( readin, .not.Lspexist, the input file does not exist ) ! if not, abort;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   open( iunit, file=trim(ext)//".sp", status="old")

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! read namelists one after another
   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : reading physicslist     from ext.sp ;")') cput-cpus
   endif

   read(iunit,physicslist)

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : read    physicslist     from ext.sp ;")') cput-cpus
   endif

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : reading numericlist     from ext.sp ;")') cput-cpus
   endif

   read(iunit,numericlist)

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : read    numericlist     from ext.sp ;")') cput-cpus
   endif

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : reading locallist      from ext.sp ;")') cput-cpus
   endif

   read(iunit,locallist)

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : read    locallist      from ext.sp ;")') cput-cpus
   endif

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : reading globallist   from ext.sp ;")') cput-cpus
   endif

   read(iunit,globallist)

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : read    globallist   from ext.sp ;")') cput-cpus
   endif

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : reading diagnosticslist from ext.sp ;")') cput-cpus
   endif

   read(iunit,diagnosticslist)

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : read    diagnosticslist from ext.sp ;")') cput-cpus
   endif

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : reading screenlist      from ext.sp ;")') cput-cpus
   endif

   read(iunit,screenlist)

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : read    screenlist      from ext.sp ;")') cput-cpus
   endif

end subroutine ! read_inputlists_from_file

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine check_inputs()

   use allglobal, only: cpus, MPI_COMM_SPEC, myid, Mvol
   use numerical
   use constants
   use fileunits

   LOCALS

   INTEGER              :: vvol
   REAL                 :: xx, toroidalflux, toroidalcurrent

!latex \subsubsection{reading \type{physicslist}}

!latex \begin{enumerate}
!latex \item The internal variable, \type{Mvol = Nvol + Lfreebound}, gives the number of computational domains.
!latex \item The input value for the fluxes enclosed within each interface, \inputvar{tflux(1:Mvol)} and \inputvar{tflux(1:Mvol)}, are immediately normalized:
!latex
!latex       \inputvar{tflux(1:Mvol)} $\rightarrow$ \inputvar{tflux(1:Mvol)/tflux(Nvol)}.
!latex
!latex       \inputvar{pflux(1:Mvol)} $\rightarrow$ \inputvar{pflux(1:Mvol)/tflux(Nvol)}.
!latex
!latex       (The input $\Phi_{edge} \equiv $ \inputvar{phiedge} will provide the total toroidal flux; see \link{preset}.)
!latex \item The input value for the toroidal current constraint (\inputvar{Isurf(1:Mvol)} and \inputvar{Ivolume(1:Mvol)}) are also immediately normalized, using \inputvar{curtor}.
!latex
!latex        $Ivolume \rightarrow Ivolume\cdot \frac{curtor}{\sum_i Isurf_i + Ivolume_i}$
!latex
!latex        $Isurf \rightarrow Isurf\cdot \frac{curtor}{\sum_i Isurf_i + Ivolume_i}$
!latex \end{enumerate}

   ! was: reading physicslist

   Mvol = Nvol + Lfreebound ! this is just for screen output and initial check; true assignment of Mvol appears outside if( myid.eq.0 ) then ;

   write(ounit,'("readin : ", 10x ," : ")')

   write(ounit,1010) cput-cpus, Igeometry, Istellsym, Lreflect
   write(ounit,1011)            Lfreebound, phiedge, curtor, curpol
   write(ounit,1012)            gamma
   write(ounit,1013)            Nfp, Nvol, Mvol, Mpol, Ntor
   write(ounit,1014)            pscale, Ladiabatic, Lconstraint, mupftol, mupfits
   write(ounit,1015)            Lrad(1:min(Mvol,32))

1010 format("readin : ",f10.2," : Igeometry=",i3," ; Istellsym=",i3," ; Lreflect="i3" ;")
1011 format("readin : ", 10x ," : Lfreebound=",i3," ; phiedge="es23.15" ; curtor="es23.15" ; curpol="es23.15" ;")
1012 format("readin : ", 10x ," : gamma="es23.15" ;")
1013 format("readin : ", 10x ," : Nfp=",i3," ; Nvol=",i3," ; Mvol=",i3," ; Mpol=",i3," ; Ntor=",i3," ;")
1014 format("readin : ", 10x ," : pscale="es13.5" ; Ladiabatic="i2" ; Lconstraint="i3" ; mupf: tol,its="es10.2" ,"i4" ;")
1015 format("readin : ", 10x ," : Lrad = "257(i2,",",:))

#ifdef DEBUG
   if( Wreadin ) then
    write(ounit,'("readin : ",f10.2," : tflux    ="257(es11.3" ,":))') cput-cpus, (    tflux(vvol), vvol = 1, Mvol )
    write(ounit,'("readin : ",f10.2," : pflux    ="257(es11.3" ,":))') cput-cpus, (    pflux(vvol), vvol = 1, Mvol )
    write(ounit,'("readin : ",f10.2," : helicity ="256(es11.3" ,":))') cput-cpus, ( helicity(vvol), vvol = 1, Nvol )
    write(ounit,'("readin : ",f10.2," : pressure ="257(es11.3" ,":))') cput-cpus, ( pressure(vvol), vvol = 1, Mvol )
    write(ounit,'("readin : ",f10.2," : mu       ="257(es11.3" ,":))') cput-cpus, (       mu(vvol), vvol = 1, Mvol )
    write(ounit,'("readin : ",f10.2," : Ivolume  ="257(es11.3" ,":))') cput-cpus, (  Ivolume(vvol), vvol = 1, Mvol )
    write(ounit,'("readin : ",f10.2," : Isurf    ="257(es11.3" ,":))') cput-cpus, (    Isurf(vvol), vvol = 1, Mvol )
   endif
#endif

   FATAL( readin, Igeometry.lt.1 .or. Igeometry.gt.3,      invalid geometry )
   FATAL( readin, Nfp.le.0,                                invalid Nfp )
   FATAL( readin, Mpol.lt.0 .or. Mpol.gt.MMpol,            invalid poloidal resolution: may need to recompile with higher MMpol )
   FATAL( readin, Ntor.lt.0 .or. Ntor.gt.MNtor,            invalid toroidal resolution: may need to recompile with higher MNtor )
   FATAL( readin, Nvol.lt.1 .or. Nvol.gt.MNvol,            invalid Nvol: may need to recompile with higher MNvol )
   FATAL( readin, mupftol.le.zero,                         mupftol is too small )
   FATAL( readin, abs(one+gamma).lt.vsmall,                1+gamma appears in denominator in dforce ) ! Please check this; SRH: 27 Feb 18;
   FATAL( readin, abs(one-gamma).lt.vsmall,                1-gamma appears in denominator in fu00aa ) ! Please check this; SRH: 27 Feb 18;
   FATAL( readin, Lconstraint.lt.-1 .or. Lconstraint.gt.3, illegal Lconstraint )
   FATAL( readin, Igeometry.eq.1 .and. rpol.lt.vsmall,     poloidal extent of slab too small or negative )
   FATAL( readin, Igeometry.eq.1 .and. rtor.lt.vsmall,     toroidal extent of slab too small or negative )

   if( Istellsym.eq.1 ) then
    Rbs(-MNtor:MNtor,-MMpol:MMpol) = zero
    Zbc(-MNtor:MNtor,-MMpol:MMpol) = zero
    Rws(-MNtor:MNtor,-MMpol:MMpol) = zero
    Zwc(-MNtor:MNtor,-MMpol:MMpol) = zero
    Vnc(-MNtor:MNtor,-MMpol:MMpol) = zero
    Bnc(-MNtor:MNtor,-MMpol:MMpol) = zero
   endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   FATAL( readin, abs(tflux(Nvol)).lt. vsmall, enclosed toroidal flux cannot be zero )

   toroidalflux = tflux(Nvol) ! toroidal flux is a local variable; SRH: 27 Feb 18;

   tflux(1:Mvol) = tflux(1:Mvol) / toroidalflux ! normalize toroidal flux;
   pflux(1:Mvol) = pflux(1:Mvol) / toroidalflux ! normalize poloidal flux;

   FATAL( readin, tflux(1).lt.zero, enclosed toroidal flux cannot be zero )
   do vvol = 2, Mvol
    !FATAL( readin, tflux(vvol)-tflux(vvol-1).lt.small, toroidal flux is not monotonic )
   enddo

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!latex \subsubsection{Current profiles normalization}
!latex
!latex In case of a free boundary calculation (\inputvar{Lfreebound}=1) and using a current constraint (\inputvar{Lconstraint}=3), the current profiles are
!latex renormalized in order to match the linking current \inputvar{curtor}. More specifically,
!latex \begin{align}
!latex Isurf_i\ \rightarrow\ Isurf_i\cdot\frac{curtor}{\sum_{i=1}^{Mvol-1} Isurf_i+Ivol_i}
!latex Ivol_i\ \rightarrow\ Ivol_i\cdot\frac{curtor}{\sum_{i=1}^{Mvol-1} Isurf_i+Ivol_i}.
!latex \end{align}
!latex Finally, the volume current in the vacuum region is set to $0$.

    ! Current constraint normalization

    if ((Lfreebound.EQ.1) .and. (Lconstraint.EQ.3)) then

        Ivolume(Mvol) = Ivolume(Mvol-1) ! Ensure vacuum in vacuum region

        toroidalcurrent = Ivolume(Mvol) + sum(Isurf(1:Mvol-1))

        if( curtor.NE.0 ) then
            FATAL( readin, toroidalcurrent.EQ.0 , Incompatible current profiles and toroidal linking current)

            Ivolume(1:Mvol) = Ivolume(1:Mvol) * curtor / toroidalcurrent
            Isurf(1:Mvol-1) = Isurf(1:Mvol-1) * curtor / toroidalcurrent

        else
            FATAL( readin, toroidalcurrent.NE.0, Incompatible current profiles and toroidal linking current)

            ! No rescaling if profiles have an overall zero toroidal current
        endif
    endif


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


   do vvol = 1, Mvol
    FATAL( readin, Lrad(vvol ).lt.2, require Chebyshev resolution Lrad > 2 so that Lagrange constraints can be satisfied )
   enddo

   if (Igeometry.ge.2 .and. Lrad(1).lt.Mpol) then
     write(ounit,'("readin : ",f10.2," : Minimum Lrad(1) is Mpol, automatically adjusted it to Mpol+4")') cput-cpus
     Lrad(1) = Mpol + 4
   endif
   FATAL( readin, mupfits.le.0, must give ma01aa:hybrj a postive integer value for the maximum iterations = mupfits given on input )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{reading \type{numericlist}}

   ! was: reading numericlist

   write(ounit,'("readin : ", 10x ," : ")')

   write(ounit,1020) cput-cpus, Linitialize, LautoinitBn, Lzerovac, Ndiscrete
   write(ounit,1021)            Nquad, iMpol, iNtor
   write(ounit,1022)            Lsparse, Lsvdiota, imethod, iorder, iprecon, iotatol
   write(ounit,1023)            Lextrap, Mregular, Lrzaxis, Ntoraxis

1020 format("readin : ",f10.2," : Linitialize=",i3," ;LautoinitBn=",i3," ; Lzerovac=",i2," ; Ndiscrete="i2" ;")
1021 format("readin : ", 10x ," : Nquad="i4" ; iMpol="i4" ; iNtor="i4" ;")
1022 format("readin : ", 10x ," : Lsparse="i2" ; Lsvdiota="i2" ; imethod="i2" ; iorder="i2" ; iprecon="i2" ; iotatol="es13.5" ;")
1023 format("readin : ", 10x ," : Lextrap="i2" ; Mregular="i3" ; Lrzaxis="i2" ; Ntoraxis="i2" ;")

   FATAL( readin, Ndiscrete.le.0, error )

  !FATAL(readin, Lfreebound.eq.1 .and. Lconstraint.gt.0 .and. Lsparse.eq.0, have not implemented dense Fourier angle transformation in vacuum region )

   FATAL( readin, iotatol.gt.one, illegal value for sparse tolerance ) ! I think that the sparse iota solver is no longer implemented; SRH: 27 Feb 18;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{reading \type{locallist}}

   ! was: reading locallist

   write(ounit,'("readin : ", 10x ," : ")')

   if (LBeltrami .ne. 4 .and. Lmatsolver .ne.1) then
    write(ounit,'("readin : ", 10x ," : ***Lmatsolver set to 1 for nonlinear solver***")')
    Lmatsolver = 1
   endif

   write(ounit,1030) cput-cpus, LBeltrami, Linitgues, Lmatsolver, LGMRESprec, NiterGMRES, epsGMRES, epsILU

1030 format("readin : ",f10.2," : LBeltrami="i2" ; Linitgues="i2" ; Lmatsolver="i2" ; LGMRESprec="i2" ; NiterGMRES="i4" ; epsGMRES="es13.5" ; epsILU="es13.5" ;" )

   FATAL( readin, LBeltrami.lt.0 .or. LBeltrami.gt.7, error )
   FATAL( readin, Lmatsolver.lt.0 .or. Lmatsolver.gt.3, error )
   FATAL( readin, LGMRESprec.lt.0 .or. LGMRESprec.gt.1, error )
   FATAL( readin, NiterGMRES.lt.0, error )
   FATAL( readin, abs(epsGMRES).le.machprec , error )
   FATAL( readin, abs(epsILU).le.machprec , error )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{reading \type{globallist}}

   ! was: reading globallist

   write(ounit,'("readin : ", 10x ," : ")')

   write(ounit,1040) cput-cpus, Lfindzero
   write(ounit,1041)            escale, opsilon, pcondense, epsilon, wpoloidal, upsilon
   write(ounit,1042)            forcetol, c05xmax, c05xtol, c05factor, LreadGF
   write(ounit,1043)            mfreeits, gBntol, gBnbld
   write(ounit,1044)            vcasingeps, vcasingtol, vcasingits, vcasingper

1040 format("readin : ",f10.2," : Lfindzero="i2" ;")
1041 format("readin : ", 10x ," : escale="es13.5" ; opsilon="es13.5" ; pcondense="f7.3" ; epsilon="es13.5" ; wpoloidal="f7.4" ; upsilon="es13.5" ;")
1042 format("readin : ", 10x ," : forcetol="es13.5" ; c05xmax="es13.5" ; c05xtol="es13.5" ; c05factor="es13.5" ; LreadGF="L2" ; ")
1043 format("readin : ", 10x ," : mfreeits="i4" ; gBntol="es13.5" ; gBnbld="es13.5" ;")
1044 format("readin : ", 10x ," : vcasingeps="es13.5" ; vcasingtol="es13.5" ; vcasingits="i6" ; vcasingper="i6" ;")

   FATAL( readin, escale      .lt.zero     , error )
   FATAL( readin, pcondense   .lt.one      , error )
   FATAL( readin, abs(c05xtol).le.machprec , error )
   FATAL( readin, c05factor   .le.zero     , error )
  !FATAL( readin, mfreeits    .lt.zero     , error )

   FATAL( readin, Igeometry.eq.3 .and. pcondense.le.zero, pcondense must be positive )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{reading \type{diagnosticslist}}

   ! was: reading diagnosticslist

   write(ounit,'("readin : ", 10x ," : ")')

   write(ounit,1050) cput-cpus, odetol, nPpts
   write(ounit,1051)            LHevalues, LHevectors, LHmatrix, Lperturbed, dpp, dqq, dRZ, Lcheck, Ltiming

1050 format("readin : ",f10.2," : odetol="es10.2" ; nPpts="i6" ;")
1051 format("readin : ", 10x ," : LHevalues="L2" ; LHevectors="L2" ; LHmatrix="L2" ; Lperturbed="i2" ; dpp="i3" ; dqq="i3" ; dRZ="es16.8" ; Lcheck="i3" ; Ltiming="L2" ;")

   FATAL( readin, odetol.le.zero, input error )
  !FATAL( readin, absreq.le.zero, input error )
  !FATAL( readin, relreq.le.zero, input error )
  !FATAL( readin, absacc.le.zero, input error )
  !FATAL( readin, nPpts .lt.0   , input error )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{reading \type{screenlist}}

   ! was: reading screenlist

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   write(ounit,'("readin : ", 10x ," : ")')

end subroutine ! check_inputs

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine broadcast_inputs

  use allglobal, only: ext, cpus, myid, MPI_COMM_SPEC
  use fileunits

  LOCALS


  ClBCAST( ext        ,     100, 0 )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! broadcast namelist/physicslist/

  if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : broadcasting physicslist     from ext.sp ;")') cput-cpus
  endif

  IlBCAST( Igeometry  ,       1, 0 )
  IlBCAST( Istellsym  ,       1, 0 )
  IlBCAST( Lfreebound ,       1, 0 )
  RlBCAST( phiedge    ,       1, 0 )
  RlBCAST( curtor     ,       1, 0 )
  RlBCAST( curpol     ,       1, 0 )
  RlBCAST( gamma      ,       1, 0 )
  IlBCAST( Nfp        ,       1, 0 )
  IlBCAST( Nvol       ,       1, 0 )
  IlBCAST( Mpol       ,       1, 0 )
  IlBCAST( Ntor       ,       1, 0 )
  IlBCAST( Lrad       , MNvol+1, 0 )
  RlBCAST( tflux      , MNvol+1, 0 )
  RlBCAST( pflux      , MNvol+1, 0 )
  RlBCAST( helicity   , MNvol  , 0 )
  RlBCAST( pscale     ,       1, 0 )
  RlBCAST( pressure   , MNvol+1, 0 )
  IlBCAST( Ladiabatic ,       1, 0 )
  RlBCAST( adiabatic  , MNvol+1, 0 )
  RlBCAST( mu         , MNvol+1, 0 )
  RlBCAST( Ivolume    , MNvol+1, 0 )
  RlBCAST( Isurf      , MNvol+1, 0 )
  IlBCAST( Lconstraint,       1, 0 )
  IlBCAST( pl         , MNvol  , 0 )
  IlBCAST( ql         , MNvol  , 0 )
  IlBCAST( pr         , MNvol  , 0 )
  IlBCAST( qr         , MNvol  , 0 )
  RlBCAST( iota       , MNvol  , 0 )
  IlBCAST( lp         , MNvol  , 0 )
  IlBCAST( lq         , MNvol  , 0 )
  IlBCAST( rp         , MNvol  , 0 )
  IlBCAST( rq         , MNvol  , 0 )
  RlBCAST( oita       , MNvol  , 0 )
  RlBCAST( mupftol    ,       1, 0 )
  IlBCAST( mupfits    ,       1, 0 )
  IlBCAST( Lreflect   ,       1, 0 )
  RlBCAST( rpol       ,       1, 0 )
  RlBCAST( rtor       ,       1, 0 )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! broadcast namelist/numericlist/

  if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : broadcasting numericlist     from ext.sp ;")') cput-cpus
  endif

  IlBCAST( Linitialize, 1, 0 )
  IlBCAST( LautoinitBn, 1, 0 )
  IlBCAST( Lzerovac   , 1, 0 )
  IlBCAST( Ndiscrete  , 1, 0 )
  IlBCAST( Nquad      , 1, 0 )
  IlBCAST( iMpol      , 1, 0 )
  IlBCAST( iNtor      , 1, 0 )
  IlBCAST( Lsparse    , 1, 0 )
  IlBCAST( Lsvdiota   , 1, 0 )
  IlBCAST( imethod    , 1, 0 )
  IlBCAST( iorder     , 1, 0 )
  IlBCAST( iprecon    , 1, 0 )
  RlBCAST( iotatol    , 1, 0 )
  IlBCAST( Lextrap    , 1, 0 )
  IlBCAST( Mregular   , 1, 0 )
  IlBCAST( Lrzaxis    , 1, 0 )
  IlBCAST( Ntoraxis   , 1, 0 )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! broadcast namelist/globallist/

  if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : broadcasting globallist      from ext.sp ;")') cput-cpus
  endif

  IlBCAST( Lfindzero , 1 , 0 )
  RlBCAST( escale    , 1 , 0 )
  RlBCAST( opsilon   , 1 , 0 )
  RlBCAST( pcondense , 1 , 0 )
  RlBCAST( epsilon   , 1 , 0 )
  RlBCAST( wpoloidal , 1 , 0 )
  RlBCAST( upsilon   , 1 , 0 )
  RlBCAST( forcetol  , 1 , 0 )
  RlBCAST( c05xmax   , 1 , 0 )
  RlBCAST( c05xtol   , 1 , 0 )
  RlBCAST( c05factor , 1 , 0 )
  LlBCAST( LreadGF   , 1 , 0 )
  IlBCAST( mfreeits  , 1 , 0 )
  RlBCAST( gBntol    , 1 , 0 )
  RlBCAST( gBnbld    , 1 , 0 )
  RlBCAST( vcasingeps, 1 , 0 )
  RlBCAST( vcasingtol, 1 , 0 )
  IlBCAST( vcasingits, 1 , 0 )
  IlBCAST( vcasingper, 1 , 0 )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! broadcast namelist/locallist/

  if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : broadcasting locallist       from ext.sp ;")') cput-cpus
  endif

  IlBCAST( LBeltrami    , 1, 0 )
  IlBCAST( Linitgues    , 1, 0 )
  RlBCAST( maxrndgues   , 1, 0)
  IlBCAST( Lmatsolver   , 1, 0 )
  IlBCAST( NiterGMRES   , 1, 0 )
  RlBCAST( epsGMRES     , 1, 0 )
  IlBCAST( LGMRESprec   , 1, 0 )
  RlBCAST( epsILU       , 1, 0 )
! IlBCAST( Lposdef  , 1, 0 ) ! redundant;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! broadcast namelist/diagnosticslist/

  if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : broadcasting diagnosticslist from ext.sp ;")') cput-cpus
  endif

  RlBCAST( odetol    , 1      , 0 )
 !RlBCAST( absreq    , 1      , 0 )
 !RlBCAST( relreq    , 1      , 0 )
 !RlBCAST( absacc    , 1      , 0 )
 !RlBCAST( epsr      , 1      , 0 )
  IlBCAST( nPpts     , 1      , 0 )
  RlBCAST( Ppts      , 1      , 0 )
  IlBCAST( nPtrj     , MNvol+1, 0 )
  LlBCAST( LHevalues , 1      , 0 )
  LlBCAST( LHevectors, 1      , 0 )
  LlBCAST( LHmatrix  , 1      , 0 )
  IlBCAST( Lperturbed, 1      , 0 )
  IlBCAST( dpp       , 1      , 0 )
  IlBCAST( dqq       , 1      , 0 )
  IlBCAST( Lerrortype, 1      , 0 )
  IlBCAST( Ngrid     , 1      , 0 )
  RlBCAST( dRZ       , 1      , 0 )
  IlBCAST( Lcheck    , 1      , 0 )
  LlBCAST( Ltiming   , 1      , 0 )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! broadcast namelist/screenlist/

  if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : broadcasting screenlist      from ext.sp ;")') cput-cpus
  endif

! BSCREENLIST ! broadcast screenlist; this is expanded by Makefile; do not remove;
  LlBCAST( Wreadin, 1, 0 )
  LlBCAST( Wwritin, 1, 0 ) ! redundant;
  LlBCAST( Wwrtend, 1, 0 )
  LlBCAST( Wmacros, 1, 0 )


end subroutine ! broadcast_inputs

end module inputlist
