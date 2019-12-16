!> \defgroup grp_input Input namelists and global variables

!> \file global.f90
!! \brief Defines input namelists and global variables, and opens some output files.
!! \ingroup grp_input
!! 
!! Note that all variables in namelist need to be broadcasted in readin.

!latex **input geometry**

!latex \begin{enumerate}
!latex \item The geometry of the $l$-th interface, for $l=0,N$ where $N\equiv$ \c Nvol}, is described by a set of Fourier harmonics,
!latex using an arbitrary poloidal angle,
!latex \be R_l(\t,\z)&=&\sum_{j}R_{j,l}\cos(m_j\t-n_j\z), \\ Z_l(\t,\z)&=&\sum_{j}Z_{j,l}\sin(m_j\t-n_j\z). \ee
!latex \item These harmonics are read from the \texttt{ext.sp} file and come directly after the namelists described above.
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
!latex \item The coordinate axis corresponds to $j=0$ and the outermost boundary corresponds to $j=$\c Nvol}.
!latex \item An arbitrary selection of harmonics may be inluded in any order, but only those within the range specified by \c Mpol} and \c Ntor}
!latex will be used.
!latex \item The geometry of {\em all} the interfaces, i.e. $l=0,N$, including the degenerate `coordinate-axis' interface, must be given.
!latex \end{enumerate}


! to set keyboard shortcut in emacs                                        

! (1) define macro         , e.g. \C-x \C-( . . . \C-x \C-)              
! (2) name macro           , e.g. Esc-x name-last-kbd-macro arbitraryname ! 11 Oct 12; 
! (3) set keyboard shortcut, e.g. Esc-x global-set-key F12 arbitraryname 

module constants

  implicit none

  REAL, parameter :: zero       =    0.0
  REAL, parameter :: one        =    1.0
  REAL, parameter :: two        =    2.0
  REAL, parameter :: three      =    3.0
  REAL, parameter :: four       =    4.0
  REAL, parameter :: five       =    5.0
  REAL, parameter :: six        =    6.0
  REAL, parameter :: seven      =    7.0
  REAL, parameter :: eight      =    8.0
  REAL, parameter :: nine       =    9.0
  REAL, parameter :: ten        =   10.0

  REAL, parameter :: eleven     =   11.0
  REAL, parameter :: twelve     =   12.0

  REAL, parameter :: hundred    =  100.0
  REAL, parameter :: thousand   = 1000.0

  REAL, parameter :: half       =   one / two
  REAL, parameter :: third      =   one / three 
  REAL, parameter :: quart      =   one / four
  REAL, parameter :: fifth      =   one / five
  REAL, parameter :: sixth      =   one / six

  REAL, parameter :: pi2        =   6.28318530717958623
  REAL, parameter :: pi         =   pi2 / two
  REAL, parameter :: mu0        =   2.0E-07 * pi2
  REAL, parameter :: goldenmean =   1.618033988749895 !< golden mean = \f$( 1 + \sqrt 5 ) / 2\f$ ;

  REAL, parameter :: version    =   2.00  

end module constants

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

module numerical

  implicit none

  REAL            :: machprec, vsmall, small, sqrtmachprec !< these are assigned below in readin via a call to NAG routine;
  REAL, parameter :: logtolerance = 1.0e-32 !< this is used to avoid taking alog10(zero); see e.g. dforce; 

contains

  !> \brief Duplicates NAG routine \c X02AJF (machine precision)
  !> 
  !> JAB; 27 Jul 17
  !> I suggest that this be removed; SRH: 27 Feb 18;
  REAL FUNCTION myprec()
    implicit none
    intrinsic EPSILON
    myprec = 0.5*EPSILON(small)
  END FUNCTION myprec
end module numerical

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

module fileunits

  implicit none

  INTEGER :: iunit = 10 !< input; used in global/readin:ext.sp, global/wrtend:ext.sp.end
  INTEGER :: ounit =  6 !< screen output;
  INTEGER :: gunit = 13 !< wall geometry; used in wa00aa

  INTEGER :: aunit = 11 !< vector potential; used in ra00aa:.ext.AtAzmn; 
  INTEGER :: dunit = 12 !< derivative matrix; used in newton:.ext.GF; 
  INTEGER :: hunit = 14 !< eigenvalues of Hessian; under re-construction; 
  INTEGER :: munit = 14 !< matrix elements of Hessian; 
  INTEGER :: lunit = 20 !< local unit; used in lunit+myid: pp00aa:.ext.poincare,.ext.transform; 
  INTEGER :: vunit = 15 !< for examination of adaptive quadrature; used in casing:.ext.vcint; 
 !INTEGER :: funit = 16 !< force iterations;

end module fileunits

module cputiming

! CPUVARIABLE !< this is expanded by Makefile; do not remove;

  REAL :: Treadin = 0.0
  REAL :: Twritin = 0.0 !< redundant; 
  REAL :: Twrtend = 0.0

end module cputiming

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

module typedefns
  
  type subgrid
     REAL,    allocatable :: s(:)
     INTEGER, allocatable :: i(:)
  end type subgrid

end module typedefns

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

module inputlist

  implicit none


  
  CHARACTER          :: ext*100 !< The input file is , \c ext.sp , where \c ext*100 or \c ext.sp*100 is given as command line input.

  !SET MAXIMUM RESOLUTION;
  INTEGER, parameter :: MNvol     = 256 !< The maximum value of \c Nvol is \c MNvol=256.
  INTEGER, parameter :: MMpol     =  32 !< The maximum value of \c Mpol is \c MNpol=32 .
  INTEGER, parameter :: MNtor     =  16 !< The maximum value of \c Ntor is \c MNtor=16 .

!> \addtogroup grp_input
!> @{


!> \addtogroup grp_input_physicslist physicslist
!>
!> The namelist \c physicslist controls the geometry, profiles, and numerical resolution.
!! @{
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
                                                         !< <li> \todo if \c Lconstraint==2, under reconstruction. </li>
                                                         !< 
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
  REAL         :: mupftol                    =  1.0e-16  !< accuracy to which \f$\mu\f$ and \f$\Delta\psi_p\f$ are required
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

!> \addtogroup grp_input_numerics numericlist
!>
!> The namelist \c numericlist controls internal resolution parameters that the user rarely needs to consider.
!! @{
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
                                     !< only relevant if \c Lsparse.eq.2; \see tr00ab.f90 for details
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
                                     !< only relevant if \c Lsparse.eq.2; \see tr00ab.f90 for details
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
!> @}

!> \addtogroup grp_input_local locallist
!> 
!> The namelist \c locallist controls the construction of the Beltrami fields in each volume.
!> 
!> The transformation to straight-fieldline coordinates is singular when the rotational-transform of the interfaces is rational;
!> however, the rotational-transform is still well defined.
!! @{
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
!> @}
  
!> \addtogroup grp_input_global globallist
!> 
!> The namelist \c globallist controls the search for global force-balance.
!> 
!> Comments:
!! <ul>
!! <li> The "force" vector, \f${\bf F}\f$, which is constructed in dforce(), is a combination of pressure-imbalance Fourier harmonics, 
!!       \f{eqnarray}{ F_{i,v} \equiv [[ p+B^2/2 ]]_{i,v} \times \exp\left[-\texttt{escale}(m_i^2+n_i^2) \right] \times \texttt{opsilon},
!!       \label{eq:forcebalancemn} \f}
!!       and spectral-condensation constraints, \f$I_{i,v}\f$, and the "star-like" angle constraints, \f$S_{i,v,}\f$, (see lforce.f90 for details)
!!       \f{eqnarray}{ F_{i,v} \equiv \texttt{epsilon} \times I_{i,v} 
!!           + \texttt{upsilon} \times \left( \psi_v^\omega S_{i,v,1} - \psi_{v+1}^\omega S_{i,v+1,0} \right),
!!       \label{eq:spectralbalancemn} \f}
!!       where \f$\psi_v\equiv\f$ normalized toroidal flux, \c tflux, and \f$\omega\equiv\f$ \c wpoloidal. </li>
!! </ul>
!! @{
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
                                         !< <li> defined in preset.f90; used in dforce() </li>
                                         !< <li> also see Eqn.\f$(\ref{eq:forcebalancemn})\f$ </li>
                                         !< </ul> 
  REAL         :: opsilon    =   1.0     !< weighting of force-imbalance
                                         !< <ul>
                                         !< <li> used in dforce(); also see Eqn.\f$(\ref{eq:forcebalancemn})\f$ </li>
                                         !< </ul>
  REAL         :: pcondense  =   2.0     !< spectral condensation parameter
                                         !< <ul>
                                         !< <li> used in preset() to define \c mmpp(i) \f$\equiv m_i^p\f$, where \f$p\equiv \f$ \c pcondense </li>
                                         !< <li> the angle freedom is exploited to minimize \f$\displaystyle \texttt{epsilon} \sum_{i} m_i^p (R_{i}^2+Z_{i}^2)\f$
                                         !<      with respect to tangential variations in the interface geometry </li>
                                         !< <li> also see Eqn.\f$(\ref{eq:spectralbalancemn})\f$ </li>vi 
                                         !< </ul>
  REAL         :: epsilon    =   0.0     !< weighting of spectral-width constraint
                                         !< <ul>
                                         !< <li> used in dforce(); also see Eqn.\f$(\ref{eq:spectralbalancemn})\f$ </li>
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
                                         !< <li> only used in xspech(); see xspech.f90 for more documentation </li>
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
  INTEGER      :: mcasingcal =   8       !< minimum number of calls to adaptive virtual casing routine; see casing() redundant;
!> @}


!> \addtogroup grp_input_diagnostics diagnosticslist
!> 
!> The namelist \c diagnosticslist controls post-processor diagnostics, such as Poincaré  plot resolution, etc.
!> @{
  REAL         :: odetol           =     1.0e-07 !< o.d.e. integration tolerance for all field line tracing routines
  REAL         :: absreq           =     1.0e-08 !< redundant
  REAL         :: relreq           =     1.0e-08 !< redundant
  REAL         :: absacc           =     1.0e-04 !< redundant
  REAL         :: epsr             =     1.0e-08 !< redundant
  INTEGER      :: nPpts            =     0       !< number of toroidal transits used (per trajectory) in following field lines
                                                 !< for constructing Poincaré plots;
                                                 !< if \c nPpts<1, no Poincaré plot is constructed;
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


!> \addtogroup grp_input_screenlist screenlist
!> 
!> The namelist \c screenlist controls screen output.
!> @{  
! DSCREENLIST !< define screenlist; this is expanded by Makefile; DO NOT REMOVE; each file compiled by Makefile has its own write flag;
  LOGICAL      :: Wreadin = .false. !< Every subroutine, e.g. \c xy00aa.h, has its own write flag, \c Wxy00aa.
  LOGICAL      :: Wwritin = .false. !< redundant
  LOGICAL      :: Wwrtend = .false.
  LOGICAL      :: Wmacros = .false.
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
 Mregular      

  namelist/locallist/&
 LBeltrami,&    
 Linitgues,&    
 maxrndgues,&   
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
 nPtrj      ,&  
 LHevalues  ,&  
 LHevectors ,&  
 LHmatrix   ,&  
 Lperturbed ,&  
 dpp        ,&  
 dqq        ,&  
 Lcheck     ,&  
 Ltiming    ,& 
 fudge      ,&
 scaling      

  namelist/screenlist/&
! NSCREENLIST ! namelist screenlist; this is expanded by Makefile; DO NOT REMOVE;
 Wreadin , &  
 Wwritin , & 
 Wwrtend , &
 Wmacros

end module inputlist
!> @}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! >  \brief global variables used throughout the whole SPEC code
module allglobal

  use constants
  use typedefns

  implicit none

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER              :: myid, ncpu       !< mpi variables;
  REAL                 :: cpus             !< initial time;

  REAL                 :: pi2nfp           !< pi2/nfp     ; assigned in readin;
  REAL                 :: pi2pi2nfp
  REAL                 :: pi2pi2nfphalf
  REAL                 :: pi2pi2nfpquart

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  REAL                 :: ForceErr, Energy

  INTEGER              :: Mvol

  LOGICAL              :: YESstellsym, NOTstellsym !< internal shorthand copies of Istellsym, which is an integer input; 

  REAL   , allocatable :: cheby(:,:) !< local workspace;
  
  REAL   , allocatable :: TT(:,:,:) !< derivatives of Chebyshev polynomials at the inner and outer interfaces;

  LOGICAL, allocatable :: ImagneticOK(:)   !< used to indicate if Beltrami fields have been correctly constructed;

  REAL   , allocatable :: beltramierror(:,:)  !< to store the integral of |curlB-mu*B| computed by jo00aa;
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{internal variables}

!latex \subsubsection{Fourier representation}

!latex \begin{enumerate} 

!> \ingroup grp_input
!> \addtogroup grp_internal_vars Internal Variables
!> @{
!> \addtogroup grp_fourier_repr Fourier representation
!> @{
  INTEGER              :: mn  !< total number of Fourier harmonics for coordinates/fields; calculated from Mpol,Ntor in readin;
  INTEGER, allocatable :: im(:), in(:) !< Fourier modes; set in readin;

  REAL,    allocatable :: halfmm(:), regumm(:)

  REAL                 :: Rscale
  REAL,    allocatable :: psifactor(:,:), inifactor(:,:)

  REAL,    allocatable :: BBweight(:) !< weight on force-imbalance harmonics; used in dforce; 
  
  REAL,    allocatable :: mmpp(:) !< spectral condensation factors; 
 
! INTEGER, allocatable :: dnjn(:,:)

!latex \item Enhanced resolution is required for the metric elements, $g_{ij}/\sqrt g$, which is given by \texttt{mne}, \texttt{ime}., and \texttt{ine}.
!latex The Fourier resolution here is determined by \texttt{lMpol=2*Mpol} and \texttt{lNtor=2*Ntor}.
  INTEGER              :: mne !< enhanced resolution for metric elements;
  INTEGER, allocatable :: ime(:), ine(:)

!latex \item Enhanced resolution is required for the transformation to straight-field line angle on the interfaces,
!latex which is given by \texttt{mns}, \texttt{ims}., and \texttt{ins}.
!latex The Fourier resolution here is determined by \texttt{iMpol} and \texttt{iNtor}.
!latex \end{enumerate}

  INTEGER              :: mns !< enhanced resolution for straight field line transformation;
  INTEGER, allocatable :: ims(:), ins(:)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  REAL                 :: xoffset = 1.0 !< used to normalize NAG routines;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!> @}
!> \addtogroup grp_iface_geom Interface geometry: iRbc, iZbs etc.
!> @{
!latex \subsubsection{\texttt{iRbc, iZbs, etc.} : interface geometry}

!latex \begin{enumerate}
!latex \item The Fourier harmonics of the interfaces are contained in \texttt{iRbc(1:mn,0:Mvol)} and \texttt{iZbs(1:mn,0:Mvol)}, where
!latex \texttt{iRbc(l,j)}, \texttt{iZbs(l,j)} contains the Fourier harmonics, $R_j$, $Z_j$, of the $l$-th interface.
!latex \end{enumerate}

  REAL,    allocatable :: iRbc(:,:) , iZbs(:,:)   !< interface surface geometry;     stellarator symmetric;
  REAL,    allocatable :: iRbs(:,:) , iZbc(:,:)   !< interface surface geometry; non-stellarator symmetric;

  REAL,    allocatable :: dRbc(:,:) , dZbs(:,:)   !< interface surface geometry;     stellarator symmetric; linear deformation;
  REAL,    allocatable :: dRbs(:,:) , dZbc(:,:)   !< interface surface geometry; non-stellarator symmetric;

  REAL,    allocatable :: iRij(:,:) , iZij(:,:)   !< interface surface geometry; real space;
  REAL,    allocatable :: dRij(:,:) , dZij(:,:)   !< interface surface geometry; real space;
  REAL,    allocatable :: tRij(:,:) , tZij(:,:)   !< interface surface geometry; real space;

  REAL,    allocatable :: iVns(:)                 ! 
  REAL,    allocatable :: iBns(:)                 ! 
  REAL,    allocatable :: iVnc(:)                 !
  REAL,    allocatable :: iBnc(:)                 !

  REAL,    allocatable :: lRbc(:)   , lZbs(:)     !< local workspace;
  REAL,    allocatable :: lRbs(:)   , lZbc(:)     !< local workspace;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!> @}
!> \addtogroup grp_fourier_transform Fourier Transforms
!> @{
!latex \subsubsection{Fourier Transforms}
!latex \begin{enumerate}

!latex \item The coordinate geometry and fields are mapped to/from Fourier space and real space using FFTW3.
!latex \item The resolution of the real space grid is given by \texttt{Nt=Ndiscrete*4*Mpol} and \texttt{Nz=Ndiscrete*4*Ntor}.

  INTEGER              :: Nt, Nz, Ntz, hNt, hNz !< discrete resolution; Ntz=Nt*Nz shorthand;
  REAL                 :: soNtz !< one / sqrt (one*Ntz); shorthand;


!latex \item Various workspace arrays are allocated. 
!l tex These include \texttt{Rij(1:Ntz,0:3,0:3)} and \texttt{Zij(1:Ntz,0:3,0:3)}, which contain the coordinates in real space and their derivatives;
!latex \texttt{sg(0:3,Ntz)}, which contains the Jacobian and its derivatives;
!latex and \texttt{guv(0:6,0:3,1:Ntz)}, which contains the metric elements and their derivatives.
!latex \end{enumerate}

  REAL   , allocatable :: Rij(:,:,:), Zij(:,:,:), Xij(:,:,:), Yij(:,:,:), sg(:,:), guvij(:,:,:,:), gvuij(:,:,:) !< real-space; 10 Dec 15;
  
  INTEGER, allocatable :: ki(:,:), kijs(:,:,:), kija(:,:,:) !< identification of Fourier modes; 

  INTEGER, allocatable :: iotakkii(:), iotaksub(:,:), iotakadd(:,:), iotaksgn(:,:) !< identification of Fourier modes; 

  REAL   , allocatable :: efmn(:), ofmn(:), cfmn(:), sfmn(:) !< Fourier harmonics; dummy workspace;
  REAL   , allocatable :: evmn(:), odmn(:), comn(:), simn(:) !< Fourier harmonics; dummy workspace;

  REAL   , allocatable :: ijreal(:), ijimag(:), jireal(:), jiimag(:)
  REAL   , allocatable :: jkreal(:), jkimag(:), kjreal(:), kjimag(:)

  REAL   , allocatable :: Bsupumn(:,:,:), Bsupvmn(:,:,:) !< 11 Oct 12; tangential field on interfaces; required for virtual casing construction of field;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  REAL   , allocatable :: goomne(:), goomno(:) !< described in preset.f90
  REAL   , allocatable :: gssmne(:), gssmno(:) !< described in preset.f90
  REAL   , allocatable :: gstmne(:), gstmno(:) !< described in preset.f90
  REAL   , allocatable :: gszmne(:), gszmno(:) !< described in preset.f90
  REAL   , allocatable :: gttmne(:), gttmno(:) !< described in preset.f90
  REAL   , allocatable :: gtzmne(:), gtzmno(:) !< described in preset.f90
  REAL   , allocatable :: gzzmne(:), gzzmno(:) !< described in preset.f90

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!> @}

!> \addtogroup grp_chebychev_metric Volume-integrated Chebyshev-metrics
!> @{
!latex \subsubsection{\texttt{DToocc, DToocs, DToosc, DTooss} : volume-integrated Chebyshev-metrics}
!latex \subsubsection{\texttt{TTsscc, TTsscs, TTsssc, TTssss} : volume-integrated Chebyshev-metrics}
!latex \subsubsection{\texttt{TDstcc, TDstcs, TDstsc, TDstss} : volume-integrated Chebyshev-metrics}
!latex \subsubsection{\texttt{TDszcc, TDszcs, TDszsc, TDszss} : volume-integrated Chebyshev-metrics}
!latex \subsubsection{\texttt{DDttcc, DDttcs, DDttsc, DDttss} : volume-integrated Chebyshev-metrics}
!latex \subsubsection{\texttt{DDtzcc, DDtzcs, DDtzsc, DDtzss} : volume-integrated Chebyshev-metrics}
!latex \subsubsection{\texttt{DDzzcc, DDzzcs, DDzzsc, DDzzss} : volume-integrated Chebyshev-metrics}

!latex \begin{enumerate}
!latex \item These are allocated in \link{dforce}, defined in \link{ma00aa}, and are used in \link{matrix} to construct the matrices.
!latex \end{enumerate}
  
  REAL,    allocatable :: DToocc(:,:,:,:), DToocs(:,:,:,:), DToosc(:,:,:,:), DTooss(:,:,:,:)
  REAL,    allocatable :: TTsscc(:,:,:,:), TTsscs(:,:,:,:), TTsssc(:,:,:,:), TTssss(:,:,:,:)
  REAL,    allocatable :: TDstcc(:,:,:,:), TDstcs(:,:,:,:), TDstsc(:,:,:,:), TDstss(:,:,:,:)
  REAL,    allocatable :: TDszcc(:,:,:,:), TDszcs(:,:,:,:), TDszsc(:,:,:,:), TDszss(:,:,:,:)
  REAL,    allocatable :: DDttcc(:,:,:,:), DDttcs(:,:,:,:), DDttsc(:,:,:,:), DDttss(:,:,:,:)
  REAL,    allocatable :: DDtzcc(:,:,:,:), DDtzcs(:,:,:,:), DDtzsc(:,:,:,:), DDtzss(:,:,:,:)
  REAL,    allocatable :: DDzzcc(:,:,:,:), DDzzcs(:,:,:,:), DDzzsc(:,:,:,:), DDzzss(:,:,:,:)  

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  REAL,    allocatable :: dtflux(:), dpflux(:) !< \delta \psi_{toroidal} and \delta \psi_{poloidal} in each annulus;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  REAL,    allocatable :: sweight(:) !< minimum poloidal length constraint weight; 

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!> @}

!> \addtogroup grp_vecpot Vector potential and the Beltrami linear system
!> @{
!latex \subsubsection{vector potential and the Beltrami linear system}

!latex \begin{enumerate}
!latex \item In each volume, the total degrees of freedom in the Beltrami linear system is \texttt{NAdof(1:Nvol)}.
!latex This depends on \c Mpol}, \c Ntor} and \c Lrad(vvol)}.
   
  INTEGER, allocatable :: NAdof(:) !< degrees of freedom in Beltrami fields in each annulus;

!latex \item The covariant components of the vector potential are written as
!latex       \be            A_\t & = & \sum_i \sum_{l=0}^L \Ate{i,l} \; T_{l}(s) \cos\a_i + \sum_i \sum_{l=0}^L \Ato{i,l} \; T_{l}(s) \sin\a_i \\
!latex                      A_\z & = & \sum_i \sum_{l=0}^L \Aze{i,l} \; T_{l}(s) \cos\a_i + \sum_i \sum_{l=0}^L \Azo{i,l} \; T_{l}(s) \sin\a_i ,
!latex       \ee
!latex       where $T_l(s)$ are the Chebyshev polynomials and $\a_i \equiv m_i \t - n_i \z$.

!latex \item The following internal arrays are declared in \link{preset}
!latex
!latex       \verb{dAte(0,i)%s(l){$\equiv \Ate{i,l}$
!latex
!latex       \verb{dAze(0,i)%s(l){$\equiv \Aze{i,l}$
!latex
!latex       \verb{dAto(0,i)%s(l){$\equiv \Ato{i,l}$
!latex
!latex       \verb{dAzo(0,i)%s(l){$\equiv \Azo{i,l}$
!latex \end{enumerate}

  type(subgrid), allocatable :: Ate(:,:,:), Aze(:,:,:)
  type(subgrid), allocatable :: Ato(:,:,:), Azo(:,:,:)

  INTEGER      , allocatable :: Lma(:,:), Lmb(:,:), Lmc(:,:), Lmd(:,:), Lme(:,:), Lmf(:,:), Lmg(:,:), Lmh(:,:)
! REAL         , allocatable :: Lto(    :), Lzo(    :)

! INTEGER      , allocatable :: Lmo(:,:), Lme(:,:) ! Lagrange multipliers for enforcing the boundary condition that B.n=0; 17 Dec 15;
  INTEGER      , allocatable :: Fso(:,:), Fse(:,:)

  LOGICAL                    :: Lcoordinatesingularity, Lplasmaregion, Lvacuumregion
!>@}
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!> \addtogroup grp_field_matrices Field matrices: dMA, dMB, dMC, dMD, dME, dMF
!> @{
!latex \subsubsection{\texttt{dMA, dMB, dMC, dMD, dME, dMF} : field matrices}

!latex \begin{enumerate}
!latex \item The energy, $W \equiv \int \! dv {\; \bf B}\cdot{\bf B}$, and helicity, $K\equiv \int \! dv \; {\bf A}\cdot{\bf B}$, functionals may be written
!latex      \be W & = & \frac{1}{2} \; a_i \; A_{i,j} \; a_j + a_i \; B_{i,j} \; \psi_j + \frac{1}{2} \; \psi_i \; C_{i,j} \; \psi_j \label{eq:energymatrix} \\
!latex          K & = & \frac{1}{2} \; a_i \; D_{i,j} \; a_j + a_i \; E_{i,j} \; \psi_j + \frac{1}{2} \; \psi_i \; F_{i,j} \; \psi_j \label{eq:helicitymatrix}
!latex      \ee
!latex       where ${\bf a} \equiv \{ A_{\t,e,i,l}, A_{\z,e,i,l}, A_{\t,o,i,l}, A_{\z,o,i,l}, f_{e,i}, f_{o,i} \}$ contains the independent degrees of freedom
!latex       and $\boldpsi \equiv \{\Delta \psi_t,\Delta \psi_p\}$.
!latex \item These are allocated and deallocated in \link{dforce}, assigned in \link{matrix}, and used in \link{mp00ac} and ? \link{df00aa}.
!latex \end{enumerate}

   REAL,   allocatable :: dMA(:,:), dMB(:,:)! dMC(:,:) !< energy and helicity matrices; quadratic forms; 
   REAL,   allocatable :: dMD(:,:)! dME(:,:)! dMF(:,:) !< energy and helicity matrices; quadratic forms; 

   REAL,   allocatable :: dMG(:  )

   REAL,   allocatable :: solution(:,:) !< this is allocated in dforce; used in mp00ac and ma02aa; and is passed to packab; 

!  REAL,   allocatable :: MBpsi(:), MEpsi(:) ! matrix vector products; 
   REAL,   allocatable :: MBpsi(:)           !< matrix vector products; 
!  REAL                :: psiMCpsi, psiMFpsi
!  REAL                ::           psiMFpsi

   REAL,   allocatable :: BeltramiInverse(:,:)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  REAL   , allocatable :: diotadxup(:,:,:) !< measured rotational transform on inner/outer interfaces for each volume;          d(transform)/dx; (see dforce);
  REAL   , allocatable :: dItGpdxtp(:,:,:) !< measured toroidal and poloidal current on inner/outer interfaces for each volume; d(Itor,Gpol)/dx; (see dforce);

  REAL   , allocatable :: glambda(:,:,:,:) !< save initial guesses for iterative calculation of rotational-transform; 

  INTEGER              :: lmns
!> @}
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!> \addtogroup grp_force_constr Construction of "force"
!> @{
!latex \subsubsection{construction of ``force''}
!latex \begin{enumerate}
!latex \item The force vector is comprised of \texttt{Bomn} and \texttt{Iomn}.
!latex \end{enumerate}

  REAL,    allocatable ::  Bemn(:,:,:),  Iomn(:,:), Somn(:,:,:), Pomn(:,:,:)
  REAL,    allocatable ::  Bomn(:,:,:),  Iemn(:,:), Semn(:,:,:), Pemn(:,:,:)
  REAL,    allocatable ::  BBe(:), IIo(:), BBo(:), IIe(:)
!> @}
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!> \addtogroup grp_covar_field_ifaces Covariant field on interfaces: Btemn, Bzemn, Btomn, Bzomn
!> @{
!latex \subsubsection{\texttt{Btemn, Bzemn, Btomn, Bzomn} : covariant field on interfaces}
!latex \begin{enumerate}
!latex \item The covariant field:
!latex \end{enumerate}

  REAL,    allocatable ::  Btemn(:,:,:), Bzemn(:,:,:), Btomn(:,:,:), Bzomn(:,:,:) !< covariant components of the tangential field on interfaces;
!> @}
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> \addtogroup grp_geomdof Geometrical degrees-of-freedom: LGdof, NGdof
!> @{
!latex \subsubsection{\texttt{LGdof, NGdof} : geometrical degrees-of-freedom;}    
!latex \begin{enumerate}
!latex \item The geometrical degrees-of-freedom: 
!latex \end{enumerate}

  INTEGER              :: LGdof !<       geometrical degrees of freedom associated with each interface;                   ; 
  INTEGER              :: NGdof !< total geometrical degrees of freedom                               ;                   ; 
!> @}
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!> \addtogroup grp_par_deriv_mat Parallel construction of derivative matrix
!> @{
!latex \subsubsection{parallel construction of derivative matrix}

!latex \begin{enumerate}
!latex \item The derivatives of force-balance, $[[p+B^2/2]]$, and the spectral constraints (see \link{sw03aa}), with respect to the interface geometry
!latex is constructed in parallel by \link{dforce}.
!latex \item force-balance across the $l$-th interface depends on the fields in the adjacent interfaces.
!latex \end{enumerate}

  REAL,    allocatable :: dBBdRZ(:,:,:)
  REAL,    allocatable :: dIIdRZ(:  ,:)

  REAL,    allocatable :: dFFdRZ(:,:,:,:,:) !< derivatives of B^2 at the interfaces wrt geometry     ; 
  REAL,    allocatable :: dBBdmp(:,:,:,:  ) !< derivatives of B^2 at the interfaces wrt mu and dpflux; 
!> @}
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!> \addtogroup grp_deriv_mul_polflux Derivatives of multiplier and poloidal flux with respect to geometry:  dmupfdx
!> @{
!latex \subsubsection{\texttt{dmupfdx} : derivatives of multiplier and poloidal flux with respect to geometry}

!latex \begin{enumerate}
!latex \item The information in \texttt{dmupfdx} describes how the helicity multiplier, $\mu$, and the enclosed poloidal flux, $\Delta \psi_p$, 
!latex       must vary as the geometry is varied in order to satisfy the interface transform constraint.
!latex \item The internal variable \texttt{dmupfdx(1:Mvol,1:2,1:LGdof,0:1)} is allocated/deallocated in \link{newton}, and \link{hesian} if selected.
!l tex \item The rotational transform on the inner or outer interface of a given volume depends on the magnetic field in that volume, i.e.
!l tex       \be \iotabar_\pm = \iotabar({\bf B}_\pm),
!l tex       \ee
!l tex       so that
!l tex       \be \delta \iotabar_\pm = \frac{\partial \iotabar_\pm}{\partial {\bf B}_\pm} \cdot \delta {\bf B_\pm}.
!l tex       \ee
!latex \item The magnetic field depends on the Fourier harmonics of both the inner and outer interface geometry (represented here as $x_j$),
!latex       the helicity multiplier, and the enclosed poloidal flux, i.e. ${\bf B_\pm} = {\bf B_\pm}(x_j, \mu, \Delta \psi_p)$, so that 
!latex       \be \delta {\bf B_\pm} = \frac{\partial {\bf B}_\pm}{\partial x_j          } \delta x_j
!latex                              + \frac{\partial {\bf B}_\pm}{\partial \mu          } \delta \mu
!latex                              + \frac{\partial {\bf B}_\pm}{\partial \Delta \psi_p} \delta \Delta \psi_p.
!latex       \ee
!latex \item This information is used to adjust the calculation of how force-balance, i.e. $B^2$ at the interfaces, 
!latex       varies with geometry at fixed interface rotational transform. Given 
!latex       \be B_\pm^2 = B_\pm^2 (x_j, \mu, \Delta \psi_p),
!latex       \ee
!latex       we may derive
!latex       \be \frac{\partial B_\pm^2}{\partial x_j} = \frac{\partial B_\pm^2}{\partial x_j          }                     
!latex                                                 + \frac{\partial B_\pm^2}{\partial \mu          } \frac{\partial \mu          }{\partial x_j}
!latex                                                 + \frac{\partial B_\pm^2}{\partial \Delta \psi_p} \frac{\partial \Delta \psi_p}{\partial x_j}
!latex       \ee
!latex \item The constraint to be enforced is that $\mu$ and $\Delta \psi_p$ must generally vary as the geometry is varied 
!latex       if the value of the rotational-transform constraint on the inner/outer interface is to be preserved,
!latex       \i.e. 
!latex       \be \left(\begin{array}{ccc} \ds \frac{\partial \iotabar_-}{\partial {\bf B}_-} \cdot \frac{\partial {\bf B}_-}{\partial \mu          } & , & 
!latex                                    \ds \frac{\partial \iotabar_-}{\partial {\bf B}_-} \cdot \frac{\partial {\bf B}_-}{\partial \Delta \psi_p} \\
!latex                                    \ds \frac{\partial \iotabar_+}{\partial {\bf B}_+} \cdot \frac{\partial {\bf B}_+}{\partial \mu          } & , & 
!latex                                    \ds \frac{\partial \iotabar_+}{\partial {\bf B}_+} \cdot \frac{\partial {\bf B}_+}{\partial \Delta \psi_p}
!latex                   \end{array} \right)
!latex           \left(\begin{array}{c} \ds \frac{\partial \mu}{\partial x_j} \\ \ds \frac{\partial \Delta \psi_p}{\partial x_j} \end{array} \right) = 
!latex         - \left(\begin{array}{c} \ds \frac{\partial \iotabar_-}{\partial {\bf B}_-} \cdot \frac{\partial {\bf B}_-}{\partial x_j} \\
!latex                                  \ds \frac{\partial \iotabar_+}{\partial {\bf B}_+} \cdot \frac{\partial {\bf B}_+}{\partial x_j} \end{array} \right).
!latex       \ee
!latex \item This $2\times 2$ linear equation is solved in \link{dforce};
!latex       and the derivatives of the rotational-transform are given in \internal{diotadxup}, see \link{preset}.
!latex \item A finite-difference estimate is computed if \c Lcheck.eq.4}.
!latex \end{enumerate}

  REAL,    allocatable :: dmupfdx(:,:,:,:)  !< derivatives of mu and dpflux wrt geometry at constant interface transform; 

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! REAL,    allocatable ::   DenergyDrz(:,:,:    ) ! derivatives of \int B^2 dv            wrt interface geometry;
! REAL,    allocatable :: DanalyticDrz(:,:,:    ) ! derivatives of \int B^2 dv            wrt interface geometry;
! REAL,    allocatable :: DanalyDrzDrz(:,:,:,:,:) ! derivatives of \int B^2 dv            wrt interface geometry; wrt interface geometry

! REAL,    allocatable :: DforcebalDrz(:,:,:,:,:) ! derivatives of [[p+B^2]]              wrt interface geometry;
! REAL,    allocatable :: DspectralDrz(:  ,:  ,:) ! derivatives of I \equiv R_t X + Z_t Y wrt interface geometry;
! LOGICAL, allocatable :: LforcebalDrz(:  ,:    ) ! indicates if derivative information has been previously constructed and saved in file;

  LOGICAL              :: Lhessianallocated
  REAL,    allocatable :: hessian(:,:)
  REAL,    allocatable :: dessian(:,:)
!> @}
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!> \addtogroup grp_trig Trigonometric factors
!> @{
!latex \subsubsection{trigonometric factors}

!latex \begin{enumerate}
!latex \item To facilitate construction of the metric integrals, various trigonometric identities are exploited.
!latex \item The required information is saved in
!latex

  REAL   , allocatable :: cosi(:,:), sini(:,:), gteta(:), gzeta(:)

  REAL   , allocatable :: ajk(:)   !< definition of coordinate axis; 

  REAL   , allocatable :: dRadR(:,:,:,:), dRadZ(:,:,:,:), dZadR(:,:,:,:), dZadZ(:,:,:,:) !< derivatives of coordinate axis; 
  REAL   , allocatable :: dRodR(:,  :,:), dRodZ(:,  :,:), dZodR(:,  :,:), dZodZ(:,  :,:) !< derivatives of coordinate axis; 

  INTEGER, allocatable :: djkp(:,:), djkm(:,:) !< for calculating cylindrical volume; 

!latex \item The following are used for volume integrals (see \link{volume})
!latex \be a_{i,j,k} &=& 4 \; m_k \ooint \cos(\alpha_i)\cos(\alpha_j)\cos(\alpha_k) /(2\pi)^2 , \\
!latex     b_{i,j,k} &=& 4 \; m_j \ooint \cos(\alpha_i)\sin(\alpha_j)\sin(\alpha_k) /(2\pi)^2 ,
!latex \ee 
!latex \end{enumerate}
!> @}
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> \addtogroup grp_volints Volume integrals: lBBintegral, lABintegral
!> @{
!latex \subsubsection{\texttt{lBBintegral, lABintegral} : volume integrals}
!latex \begin{enumerate}

!latex \item The energy functional, $F \equiv \sum_l F_l$, where
!latex \be F_l \equiv \left( \int_{{\cal V}_l} \frac{p_l}{\gamma-1} + \frac{B_l^2}{2} dv \right)
!latex         = \frac{P_l}{\gamma-1}V_l^{1-\gamma}+\int_{{\cal V}_l} \frac{B_l^2}{2} dv, \label{eq:energy}
!latex \ee
!latex where the second expression is derived using $p_l V_l^\gamma=P_l$, where $P_l$ is the adiabatic-constant.
!latex In \Eqn{energy}, it is implicit that ${\bf B}$ satisfies (i) the toroidal and poloidal flux constraints; 
!latex (ii) the interface constraint, ${\bf B}\cdot\nabla s=0$; and (iii) the helicity constraint (or the transform constraint)
!latex \item The derivatives of $F_l$ with respect to the inner and outer adjacent interface geometry are stored in
!latex
!latex \texttt{dFF(1:Nvol,0:1,0:mn+mn-1)}, where
!latex
!latex $         F_l                      \equiv$ \texttt{dFF(l,0,    0)}
!latex
!latex $\partial F_l / \partial R_{l-1,j} \equiv$ \verb&dFF(ll,0,   j)& 
!latex
!latex $\partial F_l / \partial Z_{l-1,j} \equiv$ \verb&dFF(ll,0,mn}j)& 
!latex
!latex $\partial F_l / \partial R_{l  ,j} \equiv$ \texttt{dFF(ll,1,   j)} 
!latex
!latex $\partial F_l / \partial Z_{l  ,j} \equiv$ \verb&dFF(ll,1,mn}j)&
!latex

!latex \item The volume integrals $\int dv$, $\int B^2 \; dv$ and $\int {\bf A}\cdot{\bf B} \; dv$ in each volume
!latex       are computed and saved in \texttt{volume(0:2,1:Nvol)}.
!latex \end{enumerate}

  REAL   , allocatable :: lBBintegral(:) !< B.B      integral;
  REAL   , allocatable :: lABintegral(:) !< A.B      integral;

! REAL                 :: dBBintegral    ! B.B      integral; derivative wrt R, Z;

! REAL   , allocatable :: oBBintegral(:) ! B.B      integral; original; used to normalize; perhaps irrelevant;


! REAL                 :: dABintegral    ! A.B      integral; derivative wrt R, Z;

  REAL   , allocatable :: vvolume(:) !< volume integral of \sqrt g; computed in volume;
  REAL                 :: dvolume    !< derivative of volume wrt interface geometry;

!> @}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> \addtogroup grp_int_global Internal global variables
!> @{
! internal global variables; internal logical variables; default values are provided here; these may be changed according to input values;
  
  INTEGER              :: ivol !< labels volume; some subroutines (called by NAG) are fixed argument list but require the volume label;
  
  REAL                 :: gBzeta !< toroidal (contravariant) field; calculated in bfield; required to convert \f$\dot \theta\f$ to \f$B^\theta\f$, \f$\dot s\f$ to \f$B^s\f$
  
  INTEGER, allocatable :: Iquad(:) !< internal copy of Nquad; 
  
  REAL   , allocatable :: gaussianweight(:,:), gaussianabscissae(:,:)
  
  LOGICAL              :: LBlinear, LBnewton, LBsequad !< controls selection of Beltrami field solver; depends on LBeltrami;
  
  REAL                 :: oRZp(1:3) !< used in mg00aa to determine (\s,\t,\z) given (R,Z,p);

!> @}
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  type derivative
     LOGICAL :: L
     INTEGER :: innout
     INTEGER :: ii
     INTEGER :: irz
     INTEGER :: issym
  end type derivative
  
  type(derivative)  :: dBdX

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! INTEGER, allocatable :: NNZ(:) ! used to count non-zero elements of Beltrami linear system; required for sparse calculations?;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! the following are miscellaneous flags required for the virtual casing field, external (vacuum) field integration,  . . .
  
! REAL                 :: gxyz(1:3) ! point at which external field is required; NAG routines employed are fixed argument list, but require position;
  
  INTEGER              :: globaljk  ! labels position       ;                ; 
  REAL, allocatable    :: Dxyz(:,:) ! computational boundary; position       ; 
  REAL, allocatable    :: Nxyz(:,:) ! computational boundary; normal         ; 
  REAL, allocatable    :: Jxyz(:,:) ! plasma        boundary; surface current; 

  REAL                 :: tetazeta(1:2)

! REAL                 :: virtualcasingfactor = one / ( four*pi * pi2 ) ! this is old factor (before toroidal flux was corrected?) ; 
  REAL                 :: virtualcasingfactor = -one / ( four*pi       ) ! this agrees with diagno; 
  
  INTEGER              :: IBerror ! for computing error in magnetic field; 

  INTEGER              :: nfreeboundaryiterations

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  INTEGER, parameter   :: Node = 2 ! best to make this global for consistency between calling and called routines; 

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOGICAL              :: first_free_bound = .false.

contains

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine readin

!latex \subsection{subroutine readin}

!latex \begin{enumerate}
!latex \item The master node reads the input namelist and sets various internal variables. The relevant quantities are then broadcast.
!latex \end{enumerate}

  use constants

  use numerical

  use fileunits

  use inputlist

  use cputiming

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS

  LOGICAL              :: Lspexist, Lchangeangle
  INTEGER              :: vvol, mm, nn, nb, imn, ix, ii, jj, ij, kk, mj, nj, mk, nk, ip, lMpol, lNtor, X02BBF, iargc, iarg, numargs, mi, ni, lvol, extlen, sppos
  REAL                 :: xx, toroidalflux
  REAL,    allocatable :: RZRZ(:,:) ! local array used for reading interface Fourier harmonics from file;
  
  CHARACTER            :: ldate*8, ltime*10, arg*100
  
  BEGIN(readin)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \subsubsection{\texttt{machprec}, \texttt{vsmall}, \texttt{small}, \texttt{sqrtmachprec} : machine precision}

!latex \begin{enumerate}
!latex \item The machine precision is determined using the Fortran 90 intrinsic function EPSILON.
!latex \end{enumerate}

  cput = GETTIME
  
  machprec = myprec() ! is this required? Why not just set real, parameter :: machprec = 1.0E-16 ? ; let's simplify the source; SRH: 27 Feb 18;

  vsmall = 100*machprec ; small = 100*vsmall ; sqrtmachprec = sqrt(machprec) ! returns machine precision;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( myid.eq.0 ) then ! only the master node reads input file and sets secondary variables;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
#ifdef CHECKNAG
   call A00AAF() ! check NAG version;
#endif
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   call date_and_time( ldate, ltime )
   
   write(ounit,'("readin : ", 10x ," : ")')
   write(ounit,1000) cput-cpus, ldate(1:4), ldate(5:6), ldate(7:8), ltime(1:2), ltime(3:4), ltime(5:6), machprec, vsmall, small
   
1000 format("readin : ",f10.2," : date="a4"/"a2"/"a2" , "a2":"a2":"a2" ; machine precision="es9.2" ; vsmall="es9.2" ; small="es9.2" ;")
   
   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : reading ext from command line ;")') cput-cpus
   endif
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{input file extension $\equiv$ command line argument}

!latex \begin{enumerate}
!latex \item The input file name, \texttt{ext}, is given as the first command line input, and the input file itself is \verb!ext.sp!
!latex \item Additional command line inputs recognized are:
!latex \begin{enumerate} 
!latex \item \texttt{-help, -h} ; will give help information to user; under construction;
!latex \item \texttt{-readin} ; will immediately set \texttt{Wreadin=T}; this may be over-ruled when \texttt{namelist/screenlist/} is read;
!latex \end{enumerate}

!latex \end{enumerate}

   call getarg( 1, ext )
   extlen = len_trim(ext)
   sppos = index(ext, ".sp", .true.) ! search for ".sp" from the back of ext
   if (sppos.eq.extlen-2) then       ! check if ext ends with ".sp";
     ext = ext(1:extlen-3)           ! if this is the case, remove ".sp" from end of ext
   endif

   if( ext .eq. "" .or. ext.eq. " " .or. ext .eq. "-h" .or. ext .eq. "-help" ) then
    ;write(ounit,'("readin : ", 10x ," : ")')
    ;write(ounit,'("readin : ", 10x ," : file extension must be given as first command line argument ; extra command line options = -help -readin ;")')
    if( ext .eq. "-h" .or. ext .eq. "-help" ) then
     write(ounit,'("readin : ", 10x ," : ")')
     write(ounit,'("readin : ", 10x ," : the input file ext.sp must contain the input namelists; see global.pdf for description ;")')
    endif
    FATAL( readin, .true., the input file does not exist) ! if not, abort;
   endif
   
   write(ounit,'("readin : ", 10x ," : ")')
   write(ounit,'("readin : ",f10.2," : ext = ",a100)') cput-cpus, ext
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   write(ounit,'("readin : ", 10x ," : ")')
   
   numargs = iargc()
   
   if( numargs.gt.1 ) then
    iarg = 1
    do while ( iarg < numargs )
     iarg = iarg + 1 ; call getarg( iarg, arg)
     select case( arg )
     case("-help","-h") ; write(ounit,'("readin : ",f10.2," : myid=",i3," : command line options = -readin ;")') cput-cpus, myid
     case("-readin"   ) ; Wreadin = .true.
     case("-p4pg"     ) ; iarg = iarg + 1 ; call getarg( iarg, arg)
     case("-p4wd"     ) ; iarg = iarg + 1 ; call getarg( iarg, arg)
     case default       ; write(ounit,'("readin : ",f10.2," : myid=",i3," : argument not recognized ; arg = ",a100)') cput-cpus, myid, arg
     end select
    enddo
   endif
      
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   inquire( file=trim(ext)//".sp", exist=Lspexist ) ! check if file exists;
   FATAL( readin, .not.Lspexist, the input file does not exist ) ! if not, abort;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   open( iunit, file=trim(ext)//".sp", status="old")
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
!latex \subsubsection{reading \texttt{physicslist}}

!latex \begin{enumerate}   
!latex \item The internal variable, \texttt{Mvol = Nvol + Lfreebound}, gives the number of computational domains.
!latex \item The input value for the fluxes enclosed within each interface, \c tflux(1:Mvol)} and \c tflux(1:Mvol)}, are immediately normalized:
!latex 
!latex       \c tflux(1:Mvol)} $\rightarrow$ \c tflux(1:Mvol)/tflux(Nvol)}.
!latex
!latex       \c pflux(1:Mvol)} $\rightarrow$ \c pflux(1:Mvol)/tflux(Nvol)}.
!latex
!latex       (The input $\Phi_{edge} \equiv $ \c phiedge} will provide the total toroidal flux; see \link{preset}.)
!latex \end{enumerate}
      
   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : reading physicslist     from ext.sp ;")') cput-cpus
   endif
   
   read(iunit,physicslist)
   
   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : read    physicslist     from ext.sp ;")') cput-cpus
   endif

   Mvol = Nvol + Lfreebound ! this is just for screen output and initial check; true assignment of Mvol appears outside if( myid.eq.0 ) then ; 
   
   write(ounit,'("readin : ", 10x ," : ")')
   
   write(ounit,1010) cput-cpus, Igeometry, Istellsym
   write(ounit,1011)            Lfreebound, phiedge, curtor, curpol
   write(ounit,1012)            gamma
   write(ounit,1013)            Nfp, Nvol, Mvol, Mpol, Ntor
   write(ounit,1014)            pscale, Ladiabatic, Lconstraint, mupftol, mupfits
   write(ounit,1015)            Lrad(1:min(Mvol,32))
   
1010 format("readin : ",f10.2," : Igeometry=",i3," ; Istellsym=",i3," ;")
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
   endif
#endif
   
   FATAL( readin, Igeometry.lt.1 .or. Igeometry.gt.3, invalid geometry )
   FATAL( readin, Nfp.le.0, invalid Nfp )
   FATAL( readin, Mpol.lt.0 .or. Mpol.gt.MMpol, invalid poloidal resolution: may need to recompile with higher MMpol )
   FATAL( readin, Ntor.lt.0 .or. Ntor.gt.MNtor, invalid toroidal resolution: may need to recompile with higher MNtor )
   FATAL( readin, Nvol.lt.1 .or. Nvol.gt.MNvol, invalid Nvol: may need to recompile with higher MNvol )
   FATAL( readin, mupftol.le.zero, mupftol is too small )
   FATAL( readin, abs(one+gamma).lt.vsmall, 1+gamma appears in denominator in dforce ) ! Please check this; SRH: 27 Feb 18;
   FATAL( readin, abs(one-gamma).lt.vsmall, 1-gamma appears in denominator in fu00aa ) ! Please check this; SRH: 27 Feb 18;
   FATAL( readin, Lconstraint.lt.-1 .or. Lconstraint.gt.2, illegal Lconstraint )
   FATAL( readin, Igeometry.eq.1 .and. rpol.lt.vsmall, poloidal extent of slab too small or negative )   
   FATAL( readin, Igeometry.eq.1 .and. rtor.lt.vsmall, toroidal extent of slab too small or negative )

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
    FATAL( readin, tflux(vvol)-tflux(vvol-1).lt.small, toroidal flux is not monotonic )
   enddo
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   do vvol = 1, Mvol
    FATAL( readin, Lrad(vvol ).lt.2, require Chebyshev resolution Lrad > 2 so that Lagrange constraints can be satisfied )
   enddo
   
   FATAL( readin, mupfits.le.0, must give ma01aa:hybrj a postive integer value for the maximum iterations = mupfits given on input )
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
!latex \subsubsection{reading \texttt{numericlist}}
   
   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : reading numericlist     from ext.sp ;")') cput-cpus
   endif
   
   read(iunit,numericlist)!,iostat=ios)
   
   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : read    numericlist     from ext.sp ;")') cput-cpus
   endif
   
   write(ounit,'("readin : ", 10x ," : ")')
   
   write(ounit,1020) cput-cpus, Linitialize, LautoinitBn, Lzerovac, Ndiscrete
   write(ounit,1021)            Nquad, iMpol, iNtor
   write(ounit,1022)            Lsparse, Lsvdiota, imethod, iorder, iprecon, iotatol
   write(ounit,1023)            Lextrap, Mregular
   
1020 format("readin : ",f10.2," : Linitialize=",i3," ;LautoinitBn=",i3," ; Lzerovac=",i2," ; Ndiscrete="i2" ;")
1021 format("readin : ", 10x ," : Nquad="i4" ; iMpol="i4" ; iNtor="i4" ;")
1022 format("readin : ", 10x ," : Lsparse="i2" ; Lsvdiota="i2" ; imethod="i2" ; iorder="i2" ; iprecon="i2" ; iotatol="es13.5" ;")
1023 format("readin : ", 10x ," : Lextrap="i2" ; Mregular="i3" ;")
   
   FATAL( readin, Ndiscrete.le.0, error )
   
  !FATAL(readin, Lfreebound.eq.1 .and. Lconstraint.gt.0 .and. Lsparse.eq.0, have not implemented dense Fourier angle transformation in vacuum region )
   
   FATAL( readin, iotatol.gt.one, illegal value for sparse tolerance ) ! I think that the sparse iota solver is no longer implemented; SRH: 27 Feb 18;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
!latex \subsubsection{reading \texttt{locallist}}
   
   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : reading locallist      from ext.sp ;")') cput-cpus
   endif
   
   read(iunit,locallist)!,iostat=ios)
   
   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : read    locallist      from ext.sp ;")') cput-cpus
   endif
   
   write(ounit,'("readin : ", 10x ," : ")')
   
   write(ounit,1030) cput-cpus, LBeltrami, Linitgues
   
1030 format("readin : ",f10.2," : LBeltrami="i2" ; Linitgues="i2" ;")
   
   FATAL( readin, LBeltrami.lt.0 .or. LBeltrami.gt.7, error )
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
!latex \subsubsection{reading \texttt{globallist}}

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : reading globallist   from ext.sp ;")') cput-cpus
   endif
   
   read(iunit,globallist)!,iostat=ios)
   
   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : read    globallist   from ext.sp ;")') cput-cpus
   endif
   
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
1044 format("readin : ", 10x ," : vcasingeps="es13.5" ; vcasingtol="es13.5" ; vcasingits="i3" ; vcasingper="i3" ;")
   
   FATAL( readin, escale      .lt.zero     , error )
   FATAL( readin, pcondense   .lt.one      , error )
   FATAL( readin, abs(c05xtol).le.machprec , error )
   FATAL( readin, c05factor   .le.zero     , error )
  !FATAL( readin, mfreeits    .lt.zero     , error )

   FATAL( readin, Igeometry.eq.3 .and. pcondense.le.zero, pcondense must be positive )
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
!latex \subsubsection{reading \texttt{diagnosticslist}}

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : reading diagnosticslist from ext.sp ;")') cput-cpus
   endif
   
   read(iunit,diagnosticslist)!,iostat=ios)
   
   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : read    diagnosticslist from ext.sp ;")') cput-cpus
   endif
   
   write(ounit,'("readin : ", 10x ," : ")')
   
   write(ounit,1050) cput-cpus, odetol, nPpts
   write(ounit,1051)            LHevalues, LHevectors, LHmatrix, Lperturbed, dpp, dqq, Lcheck, Ltiming
   
1050 format("readin : ",f10.2," : odetol="es10.2" ; nPpts="i6" ;")
1051 format("readin : ", 10x ," : LHevalues="L2" ; LHevectors="L2" ; LHmatrix="L2" ; Lperturbed="i2" ; dpp="i3" ; dqq="i3" ; Lcheck="i3" ; Ltiming="L2" ;")
   
   FATAL( readin, odetol.le.zero, input error )
  !FATAL( readin, absreq.le.zero, input error )
  !FATAL( readin, relreq.le.zero, input error )
  !FATAL( readin, absacc.le.zero, input error )
  !FATAL( readin, nPpts .lt.0   , input error )
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
!latex \subsubsection{reading \texttt{screenlist}}
   
   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : reading screenlist      from ext.sp ;")') cput-cpus
   endif
   
   read(iunit,screenlist)
   
   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : read    screenlist      from ext.sp ;")') cput-cpus
   endif
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   write(ounit,'("readin : ", 10x ," : ")')
   
  endif ! end of if myid eq 0 loop;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! broadcast command line input
  
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
  
  IlBCAST( LBeltrami, 1, 0 )
  IlBCAST( Linitgues, 1, 0 )
  RlBCAST( maxrndgues, 1, 0)
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
  IlBCAST( nPtrj     , MNvol+1, 0 )
  LlBCAST( LHevalues , 1      , 0 )
  LlBCAST( LHevectors, 1      , 0 )
  LlBCAST( LHmatrix  , 1      , 0 )
  IlBCAST( Lperturbed, 1      , 0 )
  IlBCAST( dpp       , 1      , 0 )
  IlBCAST( dqq       , 1      , 0 )
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
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! set internal parameters that depend on physicslist;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  select case( Istellsym )
  case( 0 )    ; YESstellsym = .false. ; NOTstellsym = .true.
  case( 1 )    ; YESstellsym = .true.  ; NOTstellsym = .false.
  case default ;
   FATAL( readin, .true., illegal Istellsym )
  end select

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \subsubsection{\texttt{Mvol} : total number of volumes}
!latex \begin{enumerate}
!latex \item The number of plasma volumes is \internal{Mvol}=\c Nvol}+\c Lfreebound};
!latex \end{enumerate}

  FATAL( readin, Lfreebound.lt.0 .or. Lfreebound.gt.1, illegal Lfreebound )

  Mvol = Nvol + Lfreebound

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  SALLOCATE( beltramierror,(1:Mvol,1:3), zero)  

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \subsubsection{\texttt{mn}, \texttt{im(1:mn)} and \texttt{in(1:mn)} : Fourier mode identification}
!latex \begin{enumerate}
!latex \item The Fourier description of even periodic functions is 
!latex       \be f(\t,\z) = \sum_{n=0}^{N} f_{0,n} \cos(-n\z) + \sum_{m=1}^{M}\sum_{n=-N}^{N} f_{m,n} \cos(m\t-n\z),
!latex       \ee
!latex       where the resolution is given on input, $M\equiv $ \c  Mpol} and $N\equiv $ \c  Ntor}.
!latex \item For convenience, the Fourier summations are written as
!latex       \be f(\s,\t,\z) &=& \sum_j f_j(s) \cos( m_j \t - n_j \z ),
!latex       \ee
!latex       for $j=1,$ \texttt{mn}, where \texttt{mn}$ = N + 1 +  M  ( 2 N + 1 )$.
!latex \item The integer arrays \texttt{im(1:mn)} and \texttt{in(1:mn)} contain the $m_j$ and $n_j$.
!latex \item The array \texttt{in} includes the \texttt{Nfp} factor.
!latex \end{enumerate}
  
  mn = 1 + Ntor +  Mpol * ( 2 *  Ntor + 1 ) ! Fourier resolution of interface geometry & vector potential;
  
  SALLOCATE( im, (1:mn), 0 )
  SALLOCATE( in, (1:mn), 0 )
  
  call gi00ab(  Mpol,  Ntor, Nfp, mn, im(1:mn), in(1:mn) ) ! this sets the im and in mode identification arrays;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \subsubsection{\texttt{halfmm(1:mn}, regumm(1:mn) : regularization factor}
!latex \begin{enumerate}
!latex \item The ``regularization'' factor, \texttt{halfmm(1:mn)} = \texttt{im(1:mn)} * \texttt{half}, is real.
!latex \item This is used in \link{lforce}, \link{bfield}, \link{stzxyz}, \link{coords}, \link{jo00aa}, \link{ma00aa}, \link{sc00aa} and \link{tr00ab}.
!latex \end{enumerate}
  
  SALLOCATE( halfmm, (1:mn), im(1:mn) * half )
  SALLOCATE( regumm, (1:mn), im(1:mn) * half )
  
  if( Mregular.ge.2 ) then

   where( im.gt.Mregular ) regumm = Mregular * half

  endif
  
! if( myid.eq.0 ) write(ounit,'("global : " 10x " : "i3") im ="i3" , halfmm ="f5.1" , regum ="f5.1" ;")') ( ii, im(ii), halfmm(ii), regumm(ii), ii = 1, mn )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \subsubsection{\texttt{ime} and \texttt{ine} : extended resolution Fourier mode identification}
!latex \begin{enumerate}
!latex \item The ``extended'' Fourier resolution is defined by \internal{lMpol} $ = 4 $ \c Mpol}, \internal{lNtor} $ = 4 $\c Ntor}.
!latex \end{enumerate}

! lMpol =   Mpol ; lNtor =   Ntor ! no    enhanced resolution for metrics; 
! lMpol = 2*Mpol ; lNtor = 2*Ntor !       enhanced resolution for metrics;
  lMpol = 4*Mpol ; lNtor = 4*Ntor ! extra-enhanced resolution for metrics; 
  
  mne = 1 + lNtor + lMpol * ( 2 * lNtor + 1 ) ! resolution of metrics; enhanced resolution; see metrix;

  SALLOCATE( ime, (1:mne), 0 )
  SALLOCATE( ine, (1:mne), 0 )

  call gi00ab( lMpol, lNtor, Nfp, mne, ime(1:mne), ine(1:mne) )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \subsubsection{\texttt{mns}, \texttt{ims} and \texttt{ins} : Fourier mode identification for straight-fieldline angle}

  lMpol = iMpol ; lNtor = iNtor

  if( iMpol.le.0 ) lMpol = Mpol - iMpol
  if( iNtor.le.0 ) lNtor = Ntor - iNtor
  if(  Ntor.eq.0 ) lNtor = 0
  
  mns = 1 + lNtor + lMpol * ( 2 * lNtor + 1 ) ! resolution of straight-field line transformation on interfaces; see tr00ab; soon to be redundant; 

  SALLOCATE( ims, (1:mns), 0 )
  SALLOCATE( ins, (1:mns), 0 )

  call gi00ab( lMpol, lNtor, Nfp, mns, ims(1:mns), ins(1:mns) ) ! note that the field periodicity factor is included in ins;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! set internal parameters that depend on numericlist;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! set internal parameters that depend on locallist;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! set internal parameters that depend on globallist;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! set internal parameters that depend on diagnosticslist;
  
  if( Lcheck.eq.5 ) then ; forcetol = 1.0e+12 ; nPpts = 0 ! will check Hessian using finite-differences; 
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \subsubsection{\texttt{iRbc(1:mn,0:Mvol}, \texttt{iZbs(1:mn,0:Mvol}, \texttt{iRbs(1:mn,0:Mvol} and \texttt{iZbc(1:mn,0:Mvol} : geometry}
  
!latex \begin{enumerate}
!latex \item \texttt{iRbc}, \texttt{iZbs}, \texttt{iRbs} and \texttt{iZbc} : Fourier harmonics of interface geometry;
!latex \item \texttt{iVns}, \texttt{iVnc}, \texttt{iBns} and \texttt{iBns} : Fourier harmonics of normal field at computational boundary;
!latex \end{enumerate}

  SALLOCATE( iRbc, (1:mn,0:Mvol), zero ) ! interface Fourier harmonics;
  SALLOCATE( iZbs, (1:mn,0:Mvol), zero )
  SALLOCATE( iRbs, (1:mn,0:Mvol), zero )
  SALLOCATE( iZbc, (1:mn,0:Mvol), zero )
  
  if( Lperturbed.eq.1 ) then
  SALLOCATE( dRbc, (1:mn,0:Mvol), zero ) ! interface Fourier harmonics;
  SALLOCATE( dZbs, (1:mn,0:Mvol), zero )
  SALLOCATE( dRbs, (1:mn,0:Mvol), zero )
  SALLOCATE( dZbc, (1:mn,0:Mvol), zero )
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  SALLOCATE( iVns, (1:mn), zero )
  SALLOCATE( iBns, (1:mn), zero )
  SALLOCATE( iVnc, (1:mn), zero )
  SALLOCATE( iBnc, (1:mn), zero )
  
 !SALLOCATE( lRbc, (1:mn), zero ) ! not used; SRH: 27 Feb 18;
 !SALLOCATE( lZbs, (1:mn), zero )
 !SALLOCATE( lRbs, (1:mn), zero )
 !SALLOCATE( lZbc, (1:mn), zero )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \subsubsection{\texttt{ajk} : construction of coordinate axis}

!latex \begin{enumerate}
!latex \item This is only used in \link{rzaxis} to perform the poloidal integration and is defined quite simply: \newline
!latex       \internal{ajk[i]} $\equiv 2\pi$ if $m_i =   0$, and \newline
!latex       \internal{ajk[i]} $\equiv 0   $ if $m_i \ne 0$.
!latex \end{enumerate}

  SALLOCATE( ajk, (1:mn), zero ) ! this must be allocated & assigned now, as it is used in readin; primarily used in packxi; 02 Jan 15;

  do kk = 1, mn ; mk = im(kk) ; nk = in(kk)
   
   if( mk.eq.0 ) ajk(kk) = pi2

  enddo ! end of do kk; 
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( myid.eq.0 ) then ! read plasma boundary & computational boundary; initialize interface geometry;

   if( Igeometry.eq.3 .and. Rbc(0,+1)+Rbc(0,-1).gt.zero .and. Zbs(0,+1)-Zbs(0,-1).gt.zero ) then ; Lchangeangle = .true.
   else                                                                                          ; Lchangeangle = .false.
   endif

   if( Lchangeangle ) write(ounit,'("readin : " 10x " : CHANGING ANGLE ;")')

   do ii = 1, mn ; mm = im(ii) ; nn = in(ii) / Nfp ! set plasma boundary, computational boundary; 29 Apr 15;
    
    if( Lchangeangle ) then ; jj = -1 ; kk = -nn ! change sign of poloidal angle; 
    else                    ; jj = +1 ; kk = +nn
    endif
    
    if( mm.eq.0 .and. nn.eq.0 ) then
     
     ;iRbc(ii,Nvol) = Rbc( nn, mm)                         ! plasma        boundary is ALWAYS given by namelist Rbc & Zbs;
     ;iZbs(ii,Nvol) = zero
      if( NOTstellsym ) then
     ;iRbs(ii,Nvol) = zero
     ;iZbc(ii,Nvol) = Zbc( nn, mm)
      else
     ;iRbs(ii,Nvol) = zero
     ;iZbc(ii,Nvol) = zero
      endif
     
     if( Lfreebound.eq.1 ) then

      iRbc(ii,Mvol) = Rwc( nn, mm)                         ! computational boundary is ALWAYS given by namelist Rbc & Zbs;
      iZbs(ii,Mvol) = zero
      if( NOTstellsym ) then
      iRbs(ii,Mvol) = zero
      iZbc(ii,Mvol) = Zwc( nn, mm)      
      else 
      iRbs(ii,Mvol) = zero
      iZbc(ii,Mvol) = zero
      endif

      iVns(ii     ) = zero
      iBns(ii     ) = zero
      if( NOTstellsym ) then
      iVnc(ii     ) = Vnc( nn, mm)                         ! I guess that this must be zero, because \div B = 0 ; 
      iBnc(ii     ) = Bnc( nn, mm)                         ! I guess that this must be zero, because \div B = 0 ; 
      else
      iVnc(ii     ) = zero
      iBnc(ii     ) = zero
      endif

     endif ! end of if( Lfreebound.eq.1 ) ; 
     
    else ! if( mm.eq.0 .and. nn.eq.0 ) then ; matches 
     
     ;iRbc(ii,Nvol) =   Rbc( kk, mm) + Rbc(-kk,-mm)        ! plasma        boundary is ALWAYS given by namelist Rbc & Zbs;
     ;iZbs(ii,Nvol) = ( Zbs( kk, mm) - Zbs(-kk,-mm) ) * jj 
      if( NOTstellsym ) then
     ;iRbs(ii,Nvol) = ( Rbs( kk, mm) - Rbs(-kk,-mm) ) * jj
     ;iZbc(ii,Nvol) =   Zbc( kk, mm) + Zbc(-kk,-mm)
      else
     ;iRbs(ii,Nvol) =   zero
     ;iZbc(ii,Nvol) =   zero
      endif
     
     if( Lfreebound.eq.1 ) then

      iRbc(ii,Mvol) =   Rwc( kk, mm) + Rwc(-kk,-mm)        ! computational boundary is ALWAYS given by namelist Rbc & Zbs;
      iZbs(ii,Mvol) = ( Zws( kk, mm) - Zws(-kk,-mm) ) * jj
      if( NOTstellsym ) then
      iRbs(ii,Mvol) = ( Rws( kk, mm) - Rws(-kk,-mm) ) * jj
      iZbc(ii,Mvol) =   Zwc( kk, mm) + Zwc(-kk,-mm)
      else
      iRbs(ii,Mvol) =   zero
      iZbc(ii,Mvol) =   zero
      endif

      iVns(ii     ) = ( Vns( kk, mm) - Vns(-kk,-mm) ) * jj
      iBns(ii     ) = ( Bns( kk, mm) - Bns(-kk,-mm) ) * jj
      if( NOTstellsym ) then
      iVnc(ii     ) =   Vnc( kk, mm) + Vnc(-kk,-mm)
      iBnc(ii     ) =   Bnc( kk, mm) + Bnc(-kk,-mm)
      else
      iVnc(ii     ) =   zero
      iBnc(ii     ) =   zero
      endif

     endif ! matches if( Lfreebound.eq.1 ) ; 
     
    endif ! end of if( mm.eq.0 .and. nn.eq.0 ) ; 
 
   enddo ! end of do ii = 1, mn;

     
   select case( Linitialize ) ! 24 Oct 12;
    
   case( :0 ) ! Linitialize=0 ; initial guess for geometry of the interior surfaces is given in the input file;
    
    SALLOCATE( RZRZ, (1:4,1:Nvol), zero ) ! temp array for reading input;

    if( Lchangeangle ) then ; jj = -1  ! change sign of poloidal angle; Loizu Nov 18;
    else                    ; jj = +1
    endif
    
    do ! will read in Fourier harmonics until the end of file is reached;
     
     read(iunit,*,iostat=ios) mm, nn, RZRZ(1:4,1:Nvol)   !if change of angle applies, transformation assumes m>=0 and for m=0 only n>=0;
     if( ios.ne.0 ) exit
     
     do ii = 1, mn ; mi = im(ii) ; ni = in(ii) ! loop over harmonics within range;
      if( mm.eq.0 .and. mi.eq.0 .and. nn*Nfp.eq.ni ) then
       iRbc(ii,1:Nvol-1) = RZRZ(1,1:Nvol-1) ! select relevant harmonics;
       iZbs(ii,1:Nvol-1) = RZRZ(2,1:Nvol-1) ! select relevant harmonics;
       if( NOTstellsym ) then
        iRbs(ii,1:Nvol-1) = RZRZ(3,1:Nvol-1) ! select relevant harmonics;
        iZbc(ii,1:Nvol-1) = RZRZ(4,1:Nvol-1) ! select relevant harmonics;
       else
        iRbs(ii,1:Nvol-1) = zero             ! select relevant harmonics;
        iZbc(ii,1:Nvol-1) = zero             ! select relevant harmonics;
       endif
      elseif( mm.eq.mi .and. nn*Nfp.eq.jj*ni ) then
       iRbc(ii,1:Nvol-1) = RZRZ(1,1:Nvol-1) ! select relevant harmonics;
       iZbs(ii,1:Nvol-1) = jj*RZRZ(2,1:Nvol-1) ! select relevant harmonics;
       if( NOTstellsym ) then
        iRbs(ii,1:Nvol-1) = jj*RZRZ(3,1:Nvol-1) ! select relevant harmonics;
        iZbc(ii,1:Nvol-1) = RZRZ(4,1:Nvol-1) ! select relevant harmonics;
       else
        iRbs(ii,1:Nvol-1) = zero             ! select relevant harmonics;
        iZbc(ii,1:Nvol-1) = zero             ! select relevant harmonics;
       endif
      endif
     enddo ! end of do ii;
     
    enddo ! end of do;
    
    DALLOCATE(RZRZ)
    
   end select ! end select case( Linitialize ); 
   
   if( Igeometry.eq.3 ) then
    if( Rac(0).gt.zero ) then ! user has supplied logically possible coordinate axis;
     iRbc(1:Ntor+1,0) = Rac(0:Ntor)
     iZbs(1:Ntor+1,0) = Zas(0:Ntor)
     iRbs(1:Ntor+1,0) = Ras(0:Ntor)
     iZbc(1:Ntor+1,0) = Zac(0:Ntor)
    else ! see preset for poloidal-average specification of coordinate axis and geometrical initialization; 
    endif ! end of if( Igeometry.eq.3 ) then ; 
   endif

  endif ! end of if myid.eq.0 loop; only the master will read the input file; all variables need to be broadcast;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  ; RlBCAST( iRbc(1:mn,0:Mvol), (Mvol+1)*mn, 0 )
  if( Igeometry.eq.3 ) then
   ;RlBCAST( iZbs(1:mn,0:Mvol), (Mvol+1)*mn, 0 ) ! only required for ii > 1 ;
  endif
  if( NOTstellsym ) then
   ;RlBCAST( iRbs(1:mn,0:Mvol), (Mvol+1)*mn, 0 ) ! only required for ii > 1 ;
   if( Igeometry.eq.3 ) then
    RlBCAST( iZbc(1:mn,0:Mvol), (Mvol+1)*mn, 0 )
   endif
  endif
  
  if( Lfreebound.eq.1 ) then
   ;RlBCAST( iVns(1:mn), mn, 0 ) ! only required for ii > 1 ;
   ;RlBCAST( iBns(1:mn), mn, 0 ) ! only required for ii > 1 ;
   if( NOTstellsym ) then
    RlBCAST( iVnc(1:mn), mn, 0 )
    RlBCAST( iBnc(1:mn), mn, 0 )
   endif
  endif
  
  if( Igeometry.eq.1 .or. Igeometry.eq.2 ) then
   ;iRbc(1:mn,0) = zero ! innermost volume must be trivial; this is used in volume; innermost interface is coordinate axis; 
   if( NOTstellsym ) then
    iRbs(1:mn,0) = zero ! innermost volume must be trivial; this is used in volume; 
   endif
  endif
  
  if( Igeometry.eq.3 ) then
   iZbs(1,0:Mvol) = zero ! Zbs_{m=0,n=0} is irrelevant;
  endif
  if( NOTstellsym) then
   iRbs(1,0:Mvol) = zero ! Rbs_{m=0,n=0} is irrelevant;
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  Rscale = iRbc(1,Mvol) ! this will be used to normalize the geometrical degrees-of-freedom; 

  if( myid.eq.0 ) write(ounit,'("readin : ", 10x ," : myid=",i3," ; Rscale=",es22.15," ;")') myid, Rscale

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(readin)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine readin

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine wrtend

!latex \subsection{subroutine wrtend}
!latex \begin{enumerate}
!latex \item The restart file is written.
!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  use constants, only :

  use numerical, only : machprec

  use fileunits, only : ounit, iunit

  use cputiming, only : Twrtend

  use inputlist

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS
  
  INTEGER              :: vvol, imn, ii, mm, nn
  
  BEGIN(wrtend)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( myid.ne.0 ) goto 9999

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
 
#ifdef DEBUG
  if( Wwrtend ) then ; cput = GETTIME ; write(ounit,'("wrtend : ",f10.2," : myid=",i3," ; opening/writing ext.sp.end ;")') cput-cpus, myid
  endif
#endif

  open(iunit,file=trim(ext)//".sp.end",status="unknown") ! restart input file;

#ifdef DEBUG
  if( Wwrtend ) then ; cput = GETTIME ; write(ounit,'("wrtend : ",f10.2," : myid=",i3," ; writing physicslist ;")') cput-cpus, myid
  endif
#endif
  
  write(iunit,'("&physicslist")')
  write(iunit,'(" Igeometry   = ",i9        )') Igeometry
  write(iunit,'(" Istellsym   = ",i9        )') Istellsym
  write(iunit,'(" Lfreebound  = ",i9        )') Lfreebound
  write(iunit,'(" phiedge     = ",es23.15   )') phiedge
  write(iunit,'(" curtor      = ",es23.15   )') curtor
  write(iunit,'(" curpol      = ",es23.15   )') curpol
  write(iunit,'(" gamma       = ",es23.15   )') gamma
  write(iunit,'(" Nfp         = ",i9        )') Nfp
  write(iunit,'(" Nvol        = ",i9        )') Nvol
  write(iunit,'(" Mpol        = ",i9        )') Mpol
  write(iunit,'(" Ntor        = ",i9        )') Ntor
  write(iunit,'(" Lrad        = ",257i23    )') Lrad(1:Mvol)
  write(iunit,'(" tflux       = ",257es23.15)') tflux(1:Mvol)
  write(iunit,'(" pflux       = ",257es23.15)') pflux(1:Mvol)
  write(iunit,'(" helicity    = ",256es23.15)') helicity(1:Mvol)
  write(iunit,'(" pscale      = ",es23.15   )') pscale
  write(iunit,'(" Ladiabatic  = ",i9        )') Ladiabatic
  write(iunit,'(" pressure    = ",257es23.15)') pressure(1:Mvol)
  write(iunit,'(" adiabatic   = ",257es23.15)') adiabatic(1:Mvol)
  write(iunit,'(" mu          = ",257es23.15)') mu(1:Mvol)
  write(iunit,'(" Lconstraint = ",i9        )') Lconstraint
  write(iunit,'(" pl          = ",257i23    )') pl(0:Mvol)
  write(iunit,'(" ql          = ",257i23    )') ql(0:Mvol)
  write(iunit,'(" pr          = ",257i23    )') pr(0:Mvol)
  write(iunit,'(" qr          = ",257i23    )') qr(0:Mvol)
  write(iunit,'(" iota        = ",257es23.15)') iota(0:Mvol)
  write(iunit,'(" lp          = ",257i23    )') lp(0:Mvol)
  write(iunit,'(" lq          = ",257i23    )') lq(0:Mvol)
  write(iunit,'(" rp          = ",257i23    )') rp(0:Mvol)
  write(iunit,'(" rq          = ",257i23    )') rq(0:Mvol)
  write(iunit,'(" oita        = ",257es23.15)') oita(0:Mvol)
  write(iunit,'(" mupftol     = ",es23.15   )') mupftol
  write(iunit,'(" mupfits     = ",i9        )') mupfits
  write(iunit,'(" rpol        = ",es23.15   )') rpol
  write(iunit,'(" rtor        = ",es23.15   )') rtor

  if( Lfreebound.eq.1 .or. Zbs(0,1).gt.zero ) then
   do ii = 1, mn ; mm = im(ii) ; nn = in(ii) / Nfp ; Rbc(nn,mm) = iRbc(ii,Nvol) ; Zbs(nn,mm) = iZbs(ii,Nvol) ; Vns(nn,mm) = iVns(ii) ; Bns(nn,mm) = iBns(ii)
                                                   ; Rbs(nn,mm) = iRbs(ii,Nvol) ; Zbc(nn,mm) = iZbc(ii,Nvol) ; Vnc(nn,mm) = iVnc(ii) ; Bnc(nn,mm) = iBnc(ii)
                                                   ; Rwc(nn,mm) = iRbc(ii,Mvol) ; Zws(nn,mm) = iZbs(ii,Mvol)
                                                   ; Rws(nn,mm) = iRbs(ii,Mvol) ; Zwc(nn,mm) = iZbc(ii,Mvol)
   enddo ! end of do ii = 1, mn;
  endif ! end of if( Lfreebound.eq.1 .or. . . . ) ; 

  !write(iunit,'(" Rac         = ",99es23.15)') Rac(0:Ntor)
  !write(iunit,'(" Zas         = ",99es23.15)') Zas(0:Ntor)
  !write(iunit,'(" Ras         = ",99es23.15)') Ras(0:Ntor)
  !write(iunit,'(" Zac         = ",99es23.15)') Zac(0:Ntor)
  
 write(iunit,'(" Rac         = ",99es23.15)') iRbc(1:Ntor+1,0)
 write(iunit,'(" Zas         = ",99es23.15)') iZbs(1:Ntor+1,0) 
 write(iunit,'(" Ras         = ",99es23.15)') iRbs(1:Ntor+1,0) 
 write(iunit,'(" Zac         = ",99es23.15)') iZbc(1:Ntor+1,0) 

  do mm = 0, Mpol ! will write out the plasma boundary harmonics; 
   do nn = -Ntor, Ntor
    
    if( mm.eq.0 .and. nn.lt.0 ) cycle ! these modes are always excluded; 13 Oct 12;
    
    select case( mm )
    case(   0:  9 )
     if( nn.lt.- 9 .and. nn.gt.-99 )     write(iunit,1000) nn, mm, Rbc(nn,mm), nn, mm, Zbs(nn,mm), nn, mm, Rbs(nn,mm), nn, mm, Zbc(nn,mm)
     if( nn.lt.  0 .and. nn.gt.- 9 )     write(iunit,1001) nn, mm, Rbc(nn,mm), nn, mm, Zbs(nn,mm), nn, mm, Rbs(nn,mm), nn, mm, Zbc(nn,mm)
     if( nn.ge.  0 .and. nn.le.  9 )     write(iunit,1002) nn, mm, Rbc(nn,mm), nn, mm, Zbs(nn,mm), nn, mm, Rbs(nn,mm), nn, mm, Zbc(nn,mm)
     if( nn.gt.  9 .and. nn.le. 99 )     write(iunit,1001) nn, mm, Rbc(nn,mm), nn, mm, Zbs(nn,mm), nn, mm, Rbs(nn,mm), nn, mm, Zbc(nn,mm)
    case(  10: 99 )
     if( nn.lt.- 9 .and. nn.gt.-99 )     write(iunit,1003) nn, mm, Rbc(nn,mm), nn, mm, Zbs(nn,mm), nn, mm, Rbs(nn,mm), nn, mm, Zbc(nn,mm)
     if( nn.lt.  0 .and. nn.gt.- 9 )     write(iunit,1004) nn, mm, Rbc(nn,mm), nn, mm, Zbs(nn,mm), nn, mm, Rbs(nn,mm), nn, mm, Zbc(nn,mm)
     if( nn.ge.  0 .and. nn.le.  9 )     write(iunit,1005) nn, mm, Rbc(nn,mm), nn, mm, Zbs(nn,mm), nn, mm, Rbs(nn,mm), nn, mm, Zbc(nn,mm)
     if( nn.gt.  9 .and. nn.le. 99 )     write(iunit,1004) nn, mm, Rbc(nn,mm), nn, mm, Zbs(nn,mm), nn, mm, Rbs(nn,mm), nn, mm, Zbc(nn,mm)
    end select ! end of select case( mm );

   enddo ! end of do nn;
  enddo ! end of do mm;
  
  do mm = 0, Mpol ! will write out the computation domain harmonics; (only relevant in free-boundary case); 
   do nn = -Ntor, Ntor
    
    if( mm.eq.0 .and. nn.lt.0 ) cycle ! these modes are always excluded; 13 Oct 12;
    
    select case( mm )
    case(   0:  9 )
     if( nn.lt.- 9 .and. nn.gt.-99 ) write(iunit,1010) nn, mm, Rwc(nn,mm), nn, mm, Zws(nn,mm), nn, mm, Rws(nn,mm), nn, mm, Zwc(nn,mm)
     if( nn.lt.  0 .and. nn.gt.- 9 ) write(iunit,1011) nn, mm, Rwc(nn,mm), nn, mm, Zws(nn,mm), nn, mm, Rws(nn,mm), nn, mm, Zwc(nn,mm)
     if( nn.ge.  0 .and. nn.le.  9 ) write(iunit,1012) nn, mm, Rwc(nn,mm), nn, mm, Zws(nn,mm), nn, mm, Rws(nn,mm), nn, mm, Zwc(nn,mm)
     if( nn.gt.  9 .and. nn.le. 99 ) write(iunit,1011) nn, mm, Rwc(nn,mm), nn, mm, Zws(nn,mm), nn, mm, Rws(nn,mm), nn, mm, Zwc(nn,mm)
    case(  10: 99 )
     if( nn.lt.- 9 .and. nn.gt.-99 ) write(iunit,1013) nn, mm, Rwc(nn,mm), nn, mm, Zws(nn,mm), nn, mm, Rws(nn,mm), nn, mm, Zwc(nn,mm)
     if( nn.lt.  0 .and. nn.gt.- 9 ) write(iunit,1014) nn, mm, Rwc(nn,mm), nn, mm, Zws(nn,mm), nn, mm, Rws(nn,mm), nn, mm, Zwc(nn,mm)
     if( nn.ge.  0 .and. nn.le.  9 ) write(iunit,1015) nn, mm, Rwc(nn,mm), nn, mm, Zws(nn,mm), nn, mm, Rws(nn,mm), nn, mm, Zwc(nn,mm)
     if( nn.gt.  9 .and. nn.le. 99 ) write(iunit,1014) nn, mm, Rwc(nn,mm), nn, mm, Zws(nn,mm), nn, mm, Rws(nn,mm), nn, mm, Zwc(nn,mm)
    end select ! end of select case( mm );
    
   enddo ! end of do nn;
  enddo ! end of do mm;
  
  do mm = 0, Mpol ! will write out the computation domain harmonics; (only relevant in free-boundary case); 
   do nn = -Ntor, Ntor
    
    if( mm.eq.0 .and. nn.lt.0 ) cycle ! these modes are always excluded; 13 Oct 12;
    
    select case( mm )
    case(   0:  9 )
     if( nn.lt.- 9 .and. nn.gt.-99 ) write(iunit,1020) nn, mm, Vns(nn,mm), nn, mm, Bns(nn,mm), nn, mm, Vnc(nn,mm), nn, mm, Bnc(nn,mm)
     if( nn.lt.  0 .and. nn.gt.- 9 ) write(iunit,1021) nn, mm, Vns(nn,mm), nn, mm, Bns(nn,mm), nn, mm, Vnc(nn,mm), nn, mm, Bnc(nn,mm)
     if( nn.ge.  0 .and. nn.le.  9 ) write(iunit,1022) nn, mm, Vns(nn,mm), nn, mm, Bns(nn,mm), nn, mm, Vnc(nn,mm), nn, mm, Bnc(nn,mm)
     if( nn.gt.  9 .and. nn.le. 99 ) write(iunit,1021) nn, mm, Vns(nn,mm), nn, mm, Bns(nn,mm), nn, mm, Vnc(nn,mm), nn, mm, Bnc(nn,mm)
    case(  10: 99 )
     if( nn.lt.- 9 .and. nn.gt.-99 ) write(iunit,1023) nn, mm, Vns(nn,mm), nn, mm, Bns(nn,mm), nn, mm, Vnc(nn,mm), nn, mm, Bnc(nn,mm)
     if( nn.lt.  0 .and. nn.gt.- 9 ) write(iunit,1024) nn, mm, Vns(nn,mm), nn, mm, Bns(nn,mm), nn, mm, Vnc(nn,mm), nn, mm, Bnc(nn,mm)
     if( nn.ge.  0 .and. nn.le.  9 ) write(iunit,1025) nn, mm, Vns(nn,mm), nn, mm, Bns(nn,mm), nn, mm, Vnc(nn,mm), nn, mm, Bnc(nn,mm)
     if( nn.gt.  9 .and. nn.le. 99 ) write(iunit,1024) nn, mm, Vns(nn,mm), nn, mm, Bns(nn,mm), nn, mm, Vnc(nn,mm), nn, mm, Bnc(nn,mm)
    end select ! end of select case( mm );
    
   enddo ! end of do nn;
  enddo ! end of do mm;

1000 format("Rbc(",i3,",",i1,")",2x,"=",es23.15," Zbs(",i3,",",i1,")",2x,"=",es23.15," Rbs(",i3,",",i1,")",2x,"=",es23.15," Zbc(",i3,",",i1,")",2x,"=",es23.15)
1001 format("Rbc(",i2,",",i1,")",3x,"=",es23.15," Zbs(",i2,",",i1,")",3x,"=",es23.15," Rbs(",i2,",",i1,")",3x,"=",es23.15," Zbc(",i2,",",i1,")",3x,"=",es23.15)
1002 format("Rbc(",i1,",",i1,")",4x,"=",es23.15," Zbs(",i1,",",i1,")",4x,"=",es23.15," Rbs(",i1,",",i1,")",4x,"=",es23.15," Zbc(",i1,",",i1,")",4x,"=",es23.15)
1003 format("Rbc(",i3,",",i2,")",1x,"=",es23.15," Zbs(",i3,",",i2,")",1x,"=",es23.15," Rbs(",i3,",",i2,")",1x,"=",es23.15," Zbc(",i3,",",i2,")",1x,"=",es23.15)
1004 format("Rbc(",i2,",",i2,")",2x,"=",es23.15," Zbs(",i2,",",i2,")",2x,"=",es23.15," Rbs(",i2,",",i2,")",2x,"=",es23.15," Zbc(",i2,",",i2,")",2x,"=",es23.15)
1005 format("Rbc(",i1,",",i2,")",3x,"=",es23.15," Zbs(",i1,",",i2,")",3x,"=",es23.15," Rbs(",i1,",",i2,")",3x,"=",es23.15," Zbc(",i1,",",i2,")",3x,"=",es23.15)
  
1010 format("Rwc(",i3,",",i1,")",2x,"=",es23.15," Zws(",i3,",",i1,")",2x,"=",es23.15," Rws(",i3,",",i1,")",2x,"=",es23.15," Zwc(",i3,",",i1,")",2x,"=",es23.15)
1011 format("Rwc(",i2,",",i1,")",3x,"=",es23.15," Zws(",i2,",",i1,")",3x,"=",es23.15," Rws(",i2,",",i1,")",3x,"=",es23.15," Zwc(",i2,",",i1,")",3x,"=",es23.15)
1012 format("Rwc(",i1,",",i1,")",4x,"=",es23.15," Zws(",i1,",",i1,")",4x,"=",es23.15," Rws(",i1,",",i1,")",4x,"=",es23.15," Zwc(",i1,",",i1,")",4x,"=",es23.15)
1013 format("Rwc(",i3,",",i2,")",1x,"=",es23.15," Zws(",i3,",",i2,")",1x,"=",es23.15," Rws(",i3,",",i2,")",1x,"=",es23.15," Zwc(",i3,",",i2,")",1x,"=",es23.15)
1014 format("Rwc(",i2,",",i2,")",2x,"=",es23.15," Zws(",i2,",",i2,")",2x,"=",es23.15," Rws(",i2,",",i2,")",2x,"=",es23.15," Zwc(",i2,",",i2,")",2x,"=",es23.15)
1015 format("Rwc(",i1,",",i2,")",3x,"=",es23.15," Zws(",i1,",",i2,")",3x,"=",es23.15," Rws(",i1,",",i2,")",3x,"=",es23.15," Zwc(",i1,",",i2,")",3x,"=",es23.15)
  
1020 format("Vns(",i3,",",i1,")",2x,"=",es23.15," Bns(",i3,",",i1,")",2x,"=",es23.15," Vnc(",i3,",",i1,")",2x,"=",es23.15," Bnc(",i3,",",i1,")",2x,"=",es23.15)
1021 format("Vns(",i2,",",i1,")",3x,"=",es23.15," Bns(",i2,",",i1,")",3x,"=",es23.15," Vnc(",i2,",",i1,")",3x,"=",es23.15," Bnc(",i2,",",i1,")",3x,"=",es23.15)
1022 format("Vns(",i1,",",i1,")",4x,"=",es23.15," Bns(",i1,",",i1,")",4x,"=",es23.15," Vnc(",i1,",",i1,")",4x,"=",es23.15," Bnc(",i1,",",i1,")",4x,"=",es23.15)
1023 format("Vns(",i3,",",i2,")",1x,"=",es23.15," Bns(",i3,",",i2,")",1x,"=",es23.15," Vnc(",i3,",",i2,")",1x,"=",es23.15," Bnc(",i3,",",i2,")",1x,"=",es23.15)
1024 format("Vns(",i2,",",i2,")",2x,"=",es23.15," Bns(",i2,",",i2,")",2x,"=",es23.15," Vnc(",i2,",",i2,")",2x,"=",es23.15," Bnc(",i2,",",i2,")",2x,"=",es23.15)
1025 format("Vns(",i1,",",i2,")",3x,"=",es23.15," Bns(",i1,",",i2,")",3x,"=",es23.15," Vnc(",i1,",",i2,")",3x,"=",es23.15," Bnc(",i1,",",i2,")",3x,"=",es23.15)
  
  write(iunit,'("/")')

  if( Wwrtend ) then ; cput = GETTIME ; write(ounit,'("wrtend : ",f10.2," : myid=",i3," ; writing numericlist ;")') cput-cpus, myid
  endif
  
  write(iunit,'("&numericlist")')
  write(iunit,'(" Linitialize = ",i9            )') Linitialize
  write(iunit,'(" LautoinitBn = ",i9            )') LautoinitBn
  write(iunit,'(" Lzerovac    = ",i9            )') Lzerovac       
  write(iunit,'(" Ndiscrete   = ",i9            )') Ndiscrete
  write(iunit,'(" Nquad       = ",i9            )') Nquad
  write(iunit,'(" iMpol       = ",i9            )') iMpol
  write(iunit,'(" iNtor       = ",i9            )') iNtor
  write(iunit,'(" Lsparse     = ",i9            )') Lsparse
  write(iunit,'(" Lsvdiota    = ",i9            )') Lsvdiota
  write(iunit,'(" imethod     = ",i9            )') imethod
  write(iunit,'(" iorder      = ",i9            )') iorder
  write(iunit,'(" iprecon     = ",i9            )') iprecon
  write(iunit,'(" iotatol     = ",es23.15       )') iotatol
  write(iunit,'(" Lextrap     = ",i9            )') Lextrap
  write(iunit,'(" Mregular    = ",i9            )') Mregular
  write(iunit,'("/")')

  if( Wwrtend ) then ; cput = GETTIME ; write(ounit,'("wrtend : ",f10.2," : myid=",i3," ; writing locallist ;")') cput-cpus, myid
  endif
  
  write(iunit,'("&locallist")')
  write(iunit,'(" LBeltrami   = ",i9            )') LBeltrami
  write(iunit,'(" Linitgues   = ",i9            )') Linitgues
 !write(iunit,'(" Lposdef     = ",i9            )') Lposdef ! redundant;
 !write(iunit,'(" Nmaxexp     = ",i9            )') Nmaxexp
  write(iunit,'("/")')

  if( Wwrtend ) then ; cput = GETTIME ; write(ounit,'("wrtend : ",f10.2," : myid=",i3," ; writing globallist ;")') cput-cpus, myid
  endif
  
  write(iunit,'("&globallist")')
  write(iunit,'(" Lfindzero   = ",i9            )') Lfindzero
  write(iunit,'(" escale      = ",es23.15       )') escale
  write(iunit,'(" opsilon     = ",es23.15       )') opsilon
  write(iunit,'(" pcondense   = ",es23.15       )') pcondense
  write(iunit,'(" epsilon     = ",es23.15       )') epsilon
  write(iunit,'(" wpoloidal   = ",es23.15       )') wpoloidal
  write(iunit,'(" upsilon     = ",es23.15       )') upsilon
  write(iunit,'(" forcetol    = ",es23.15       )') forcetol
  write(iunit,'(" c05xmax     = ",es23.15       )') c05xmax
  write(iunit,'(" c05xtol     = ",es23.15       )') c05xtol
  write(iunit,'(" c05factor   = ",es23.15       )') c05factor
  write(iunit,'(" LreadGF     = ",L9            )') LreadGF
  write(iunit,'(" mfreeits    = ",i9            )') mfreeits
  write(iunit,'(" gBntol      = ",es23.15       )') gBntol
  write(iunit,'(" gBnbld      = ",es23.15       )') gBnbld
  write(iunit,'(" vcasingeps  = ",es23.15       )') vcasingeps   
  write(iunit,'(" vcasingtol  = ",es23.15       )') vcasingtol   
  write(iunit,'(" vcasingits  = ",i9            )') vcasingits   
  write(iunit,'(" vcasingper  = ",i9            )') vcasingper   
  write(iunit,'("/")')

  if( Wwrtend ) then ; cput = GETTIME ; write(ounit,'("wrtend : ",f10.2," : myid=",i3," ; writing diagnosticslist ;")') cput-cpus, myid
  endif
  
  write(iunit,'("&diagnosticslist")')
  write(iunit,'(" odetol      = ",es23.15       )') odetol
 !write(iunit,'(" absreq      = ",es23.15       )') absreq
 !write(iunit,'(" relreq      = ",es23.15       )') relreq
 !write(iunit,'(" absacc      = ",es23.15       )') absacc
 !write(iunit,'(" epsr        = ",es23.15       )') epsr
  write(iunit,'(" nPpts       = ",i9            )') nPpts
  write(iunit,'(" nPtrj       = ",256i6         )') nPtrj(1:Mvol)
  write(iunit,'(" LHevalues   = ",L9            )') LHevalues
  write(iunit,'(" LHevectors  = ",L9            )') LHevectors
  write(iunit,'(" LHmatrix    = ",L9            )') LHmatrix
  write(iunit,'(" Lperturbed  = ",i9            )') Lperturbed
  write(iunit,'(" dpp         = ",i9            )') dpp
  write(iunit,'(" dqq         = ",i9            )') dqq
  write(iunit,'(" Lcheck      = ",i9            )') Lcheck
  write(iunit,'(" Ltiming     = ",L9            )') Ltiming
  write(iunit,'("/")')

  if( Wwrtend ) then ; cput = GETTIME ; write(ounit,'("wrtend : ",f10.2," : myid=",i3," ; writing screenlist ;")') cput-cpus, myid
  endif
  
  write(iunit,'("&screenlist")')
! WSCREENLIST ! write screenlist; this is expanded by Makefile ; do not remove;
  if( Wreadin           ) write(iunit,'(" Wreadin = ",L1                )') Wreadin
  if( Wwrtend           ) write(iunit,'(" Wwrtend = ",L1                )') Wwrtend
  if( Wmacros           ) write(iunit,'(" Wmacros = ",L1                )') Wmacros
  write(iunit,'("/")')
  
#ifdef DEBUG
  FATAL( wrtend, .not.allocated(iRbc), error )
  FATAL( wrtend, .not.allocated(iZbs), error )
  FATAL( wrtend, .not.allocated(iRbs), error )
  FATAL( wrtend, .not.allocated(iZbc), error )
#endif

  do imn = 1, mn ; write(iunit,'(2i6,1024es23.15)') im(imn), in(imn)/Nfp, ( iRbc(imn,vvol), iZbs(imn,vvol), iRbs(imn,vvol), iZbc(imn,vvol), vvol = 1, Nvol )
  enddo
  
  close(iunit)
  
#ifdef DEBUG
  if( Wwrtend ) then ; cput = GETTIME ; write(ounit,'("wrtend : ",f10.2," : myid=",i3," ; wrote ext.sp.end ;")') cput-cpus, myid
  endif
#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(wrtend)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine wrtend
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end module allglobal

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

module fftw_interface ! JAB; 25 Jul 17

  use, intrinsic :: iso_c_binding

  implicit none

  include 'fftw3.f03'

  TYPE(C_PTR)                            :: planf, planb
  COMPLEX(C_DOUBLE_COMPLEX), allocatable :: cplxin(:,:), cplxout(:,:)

end module fftw_interface

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
