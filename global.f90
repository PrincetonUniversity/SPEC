!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
!        1         2         3         4         5         6         7         8         9        10        11        12        13        14        15
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (input) ! Defines input namelists and global variables, and opens some output files.

!latex \briefly{briefly}

!latex \calledby{\link{}}
!latex \calls{\link{}}

!latex \tableofcontents

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! to set keyboard shortcut in emacs

! (1) define macro         , e.g. \C-x \C-( . . . \C-x \C-)
! (2) name macro           , e.g. Esc-x name-last-kbd-macro arbitraryname ! 11 Oct 12;
! (3) set keyboard shortcut, e.g. Esc-x global-set-key F12 arbitraryname

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

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
  REAL, parameter :: goldenmean =   1.618033988749895 ! golden mean = ( one + sqrt(five) ) / two ;

  ! version of SPEC
  REAL, parameter :: version    =   3.10

end module constants

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

module numerical

  implicit none

  REAL            :: machprec, vsmall, small, sqrtmachprec ! these are assigned below in readin via a call to NAG routine;
  REAL, parameter :: logtolerance = 1.0e-32 ! this is used to avoid taking alog10(zero); see e.g. dforce;

contains
  REAL FUNCTION myprec() !Duplicates NAG routine X02AJF (machine precision) ! JAB; 27 Jul 17 ! I suggest that this be removed; SRH: 27 Feb 18;
    implicit none
    intrinsic EPSILON
    myprec = 0.5*EPSILON(small)
  END FUNCTION myprec
end module numerical

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

module fileunits

  implicit none

  INTEGER :: iunit = 10 ! input; used in global/readin:ext.sp, global/wrtend:ext.sp.end
  INTEGER :: ounit =  6 ! screen output;
  INTEGER :: gunit = 13 ! wall geometry; used in wa00aa

  INTEGER :: aunit = 11 ! vector potential; used in ra00aa:.ext.AtAzmn;
  INTEGER :: dunit = 12 ! derivative matrix; used in newton:.ext.GF;
  INTEGER :: hunit = 14 ! eigenvalues of Hessian; under re-construction;
  INTEGER :: munit = 14 ! matrix elements of Hessian;
  INTEGER :: lunit = 20 ! local unit; used in lunit+myid: pp00aa:.ext.poincare,.ext.transform;
  INTEGER :: vunit = 15 ! for examination of adaptive quadrature; used in casing:.ext.vcint;
 !INTEGER :: funit = 16 ! force iterations;

  contains
    subroutine mute(action)
      implicit none

      INTEGER,intent(in) :: action
      INTEGER, parameter :: iopen = 1, iclose = 0, null = 37
      INTEGER            :: ios
      character(len=*), parameter :: nullfile="/dev/null"

      ! open a tmp file for screen output
      if (action == iopen) then
        ounit = null
        open(ounit, file=nullfile, status="unknown", action="write", iostat=ios) ! create a scratch file?
        if (ios.ne.0) print *, "something wrong with open a tmp file in focuspy.mute. IOSTAT=", ios
      else
        close(ounit)
        ounit = 6 ! recover to screen output
      endif
      return
    end subroutine mute

end module fileunits

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

module cputiming

! CPUVARIABLE ! this is expanded by Makefile; do not remove;

  REAL :: Treadin = 0.0
  REAL :: Twritin = 0.0 ! redundant;
  REAL :: Twrtend = 0.0

end module cputiming

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

module typedefns

  type subgrid
    REAL,    allocatable :: s(:)
    INTEGER, allocatable :: i(:)
  end type subgrid

  type MatrixLU
    REAL, allocatable :: mat(:,:)
    INTEGER, allocatable :: ipivot(:)
  end type MatrixLU

end module typedefns

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

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

end module inputlist

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

module allglobal

  use constants
  use typedefns

  implicit none

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!``-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER              :: myid, ncpu, MPI_COMM_SPEC       ! mpi variables;
  INTEGER              :: IsMyVolumeValue
  REAL                 :: cpus             ! initial time;

  REAL                 :: pi2nfp           !       pi2/nfp     ; assigned in readin;
  REAL                 :: pi2pi2nfp
  REAL                 :: pi2pi2nfphalf
  REAL                 :: pi2pi2nfpquart

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  CHARACTER(LEN=100)   :: ext ! extension of input filename, i.e., "G3V01L1Fi.001" for an input file G3V01L1Fi.001.sp

  REAL                 :: ForceErr, Energy

  REAL   , allocatable :: IPDt(:), IPDtDpf(:,:)    ! Toroidal pressure-driven current

  INTEGER              :: Mvol

  LOGICAL              :: YESstellsym, NOTstellsym ! internal shorthand copies of Istellsym, which is an integer input;

  LOGICAL              :: YESMatrixFree, NOTMatrixFree ! to use matrix-free method or not

  REAL   , allocatable :: cheby(:,:), zernike(:,:,:) ! local workspace;

  REAL   , allocatable :: TT(:,:,:), RTT(:,:,:,:) ! derivatives of Chebyshev and Zernike polynomials at the inner and outer interfaces;
  REAL   , allocatable :: RTM(:,:) ! r^m term of Zernike polynomials at the origin
  REAL   , allocatable :: ZernikeDof(:) ! Zernike degree of freedom for each m

  LOGICAL, allocatable :: ImagneticOK(:)   ! used to indicate if Beltrami fields have been correctly constructed;

  LOGICAL              :: IconstraintOK       ! Used to break iteration loops of slaves in the global constraint minimization.

  REAL   , allocatable :: beltramierror(:,:)  ! to store the integral of |curlB-mu*B| computed by jo00aa;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{internal variables}

!latex \subsubsection{Fourier representation}

!latex \begin{enumerate}

  INTEGER              :: mn  ! total number of Fourier harmonics for coordinates/fields; calculated from Mpol,Ntor in readin;
  INTEGER, allocatable :: im(:), in(:) ! Fourier modes; set in readin;

  REAL,    allocatable :: halfmm(:), regumm(:)

  REAL                 :: Rscale
  REAL,    allocatable :: psifactor(:,:), inifactor(:,:)

  REAL,    allocatable :: BBweight(:) ! weight on force-imbalance harmonics; used in dforce;

  REAL,    allocatable :: mmpp(:) ! spectral condensation factors;

! INTEGER, allocatable :: dnjn(:,:)

!latex \item Enhanced resolution is required for the metric elements, $g_{ij}/\sqrt g$, which is given by \type{mne}, \type{ime}., and \type{ine}.
!latex The Fourier resolution here is determined by \type{lMpol=2*Mpol} and \type{lNtor=2*Ntor}.
  INTEGER              :: mne ! enhanced resolution for metric elements;
  INTEGER, allocatable :: ime(:), ine(:)

!latex \item Enhanced resolution is required for the transformation to straight-field line angle on the interfaces,
!latex which is given by \type{mns}, \type{ims}., and \type{ins}.
!latex The Fourier resolution here is determined by \type{iMpol} and \type{iNtor}.
!latex \end{enumerate}

  INTEGER              :: mns ! enhanced resolution for straight field line transformation;
  INTEGER, allocatable :: ims(:), ins(:)

  INTEGER              :: lMpol, lNtor, sMpol, sNtor

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  REAL                 :: xoffset = 1.0 ! used to normalize NAG routines;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{\type{iRbc, iZbs, etc.} : interface geometry}

!latex \begin{enumerate}
!latex \item The Fourier harmonics of the interfaces are contained in \type{iRbc(1:mn,0:Mvol)} and \type{iZbs(1:mn,0:Mvol)}, where
!latex \type{iRbc(l,j)}, \type{iZbs(l,j)} contains the Fourier harmonics, $R_j$, $Z_j$, of the $l$-th interface.
!latex \end{enumerate}

  REAL,    allocatable :: iRbc(:,:) , iZbs(:,:)   ! interface surface geometry;     stellarator symmetric;
  REAL,    allocatable :: iRbs(:,:) , iZbc(:,:)   ! interface surface geometry; non-stellarator symmetric;

  REAL,    allocatable :: dRbc(:,:) , dZbs(:,:)   ! interface surface geometry;     stellarator symmetric; linear deformation;
  REAL,    allocatable :: dRbs(:,:) , dZbc(:,:)   ! interface surface geometry; non-stellarator symmetric;

  REAL,    allocatable :: iRij(:,:) , iZij(:,:)   ! interface surface geometry; real space;
  REAL,    allocatable :: dRij(:,:) , dZij(:,:)   ! interface surface geometry; real space;
  REAL,    allocatable :: tRij(:,:) , tZij(:,:)   ! interface surface geometry; real space;

  REAL,    allocatable :: iVns(:)                 !
  REAL,    allocatable :: iBns(:)                 !
  REAL,    allocatable :: iVnc(:)                 !
  REAL,    allocatable :: iBnc(:)                 !

  REAL,    allocatable :: lRbc(:)   , lZbs(:)     ! local workspace;
  REAL,    allocatable :: lRbs(:)   , lZbc(:)     ! local workspace;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{Fourier Transforms}
!latex \begin{enumerate}

!latex \item The coordinate geometry and fields are mapped to/from Fourier space and real space using FFTW3.
!latex \item The resolution of the real space grid is given by \type{Nt=Ndiscrete*4*Mpol} and \type{Nz=Ndiscrete*4*Ntor}.

  INTEGER              :: Nt, Nz, Ntz, hNt, hNz ! discrete resolution; Ntz=Nt*Nz shorthand;
  REAL                 :: soNtz ! one / sqrt (one*Ntz); shorthand;


!latex \item Various workspace arrays are allocated.
!l tex These include \type{Rij(1:Ntz,0:3,0:3)} and \type{Zij(1:Ntz,0:3,0:3)}, which contain the coordinates in real space and their derivatives;
!latex \type{sg(0:3,Ntz)}, which contains the Jacobian and its derivatives;
!latex and \type{guv(0:6,0:3,1:Ntz)}, which contains the metric elements and their derivatives.
!latex \end{enumerate}

  REAL   , allocatable :: Rij(:,:,:), Zij(:,:,:), Xij(:,:,:), Yij(:,:,:), sg(:,:), guvij(:,:,:,:), gvuij(:,:,:) ! real-space; 10 Dec 15;
  REAL   , allocatable :: guvijsave(:,:,:,:)

  INTEGER, allocatable :: ki(:,:), kijs(:,:,:), kija(:,:,:) ! identification of Fourier modes;

  INTEGER, allocatable :: iotakkii(:), iotaksub(:,:), iotakadd(:,:), iotaksgn(:,:) ! identification of Fourier modes;

  REAL   , allocatable :: efmn(:), ofmn(:), cfmn(:), sfmn(:) ! Fourier harmonics; dummy workspace;
  REAL   , allocatable :: evmn(:), odmn(:), comn(:), simn(:) ! Fourier harmonics; dummy workspace;

  REAL   , allocatable :: ijreal(:), ijimag(:), jireal(:), jiimag(:)
  REAL   , allocatable :: jkreal(:), jkimag(:), kjreal(:), kjimag(:)

  REAL   , allocatable :: Bsupumn(:,:,:), Bsupvmn(:,:,:) ! 11 Oct 12; tangential field on interfaces; required for virtual casing construction of field;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  REAL   , allocatable :: goomne(:,:), goomno(:,:) ! described in preset;
  REAL   , allocatable :: gssmne(:,:), gssmno(:,:) ! described in preset;
  REAL   , allocatable :: gstmne(:,:), gstmno(:,:) ! described in preset;
  REAL   , allocatable :: gszmne(:,:), gszmno(:,:) ! described in preset;
  REAL   , allocatable :: gttmne(:,:), gttmno(:,:) ! described in preset;
  REAL   , allocatable :: gtzmne(:,:), gtzmno(:,:) ! described in preset;
  REAL   , allocatable :: gzzmne(:,:), gzzmno(:,:) ! described in preset;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{\type{DToocc, DToocs, DToosc, DTooss} : volume-integrated Chebyshev-metrics}
!latex \subsubsection{\type{TTsscc, TTsscs, TTsssc, TTssss} : volume-integrated Chebyshev-metrics}
!latex \subsubsection{\type{TDstcc, TDstcs, TDstsc, TDstss} : volume-integrated Chebyshev-metrics}
!latex \subsubsection{\type{TDszcc, TDszcs, TDszsc, TDszss} : volume-integrated Chebyshev-metrics}
!latex \subsubsection{\type{DDttcc, DDttcs, DDttsc, DDttss} : volume-integrated Chebyshev-metrics}
!latex \subsubsection{\type{DDtzcc, DDtzcs, DDtzsc, DDtzss} : volume-integrated Chebyshev-metrics}
!latex \subsubsection{\type{DDzzcc, DDzzcs, DDzzsc, DDzzss} : volume-integrated Chebyshev-metrics}

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

  REAL,    allocatable :: Tsc(:,:), Tss(:,:), Dtc(:,:), Dts(:,:), Dzc(:,:), Dzs(:,:)
  REAL,    allocatable :: Ttc(:,:), Tzc(:,:) ,Tts(:,:), Tzs(:,:)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  REAL,    allocatable :: dtflux(:), dpflux(:) ! \delta \psi_{toroidal} and \delta \psi_{poloidal} in each annulus;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  REAL,    allocatable :: sweight(:) ! minimum poloidal length constraint weight;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{vector potential and the Beltrami linear system}

!latex \begin{enumerate}
!latex \item In each volume, the total degrees of freedom in the Beltrami linear system is \type{NAdof(1:Nvol)}.
!latex This depends on \inputvar{Mpol}, \inputvar{Ntor} and \inputvar{Lrad(vvol)}.

  INTEGER, allocatable :: NAdof(:) ! degrees of freedom in Beltrami fields in each annulus;
  INTEGER, allocatable :: Nfielddof(:) ! degrees of freedom in Beltrami fields in each annulus, field only, no Lagrange multipliers;

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
  REAL         , allocatable :: Lmavalue(:,:), Lmbvalue(:,:), Lmcvalue(:,:), Lmdvalue(:,:), Lmevalue(:,:), Lmfvalue(:,:)
  REAL         , allocatable :: Lmgvalue(:,:), Lmhvalue(:,:)
! REAL         , allocatable :: Lto(    :), Lzo(    :)

! INTEGER      , allocatable :: Lmo(:,:), Lme(:,:) ! Lagrange multipliers for enforcing the boundary condition that B.n=0; 17 Dec 15;
  INTEGER      , allocatable :: Fso(:,:), Fse(:,:)

  LOGICAL                    :: Lcoordinatesingularity, Lplasmaregion, Lvacuumregion
  LOGICAL                    :: Lsavedguvij        ! flag used in matrix free

  LOGICAL             :: Localconstraint

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{\type{dMA, dMB, dMC, dMD, dME, dMF} : field matrices}

!latex \begin{enumerate}
!latex \item The energy, $W \equiv \int \! dv {\; \bf B}\cdot{\bf B}$, and helicity, $K\equiv \int \! dv \; {\bf A}\cdot{\bf B}$, functionals may be written
!latex      \be W & = & \frac{1}{2} \; a_i \; A_{i,j} \; a_j + a_i \; B_{i,j} \; \psi_j + \frac{1}{2} \; \psi_i \; C_{i,j} \; \psi_j \label{eq:energymatrix} \\
!latex          K & = & \frac{1}{2} \; a_i \; D_{i,j} \; a_j + a_i \; E_{i,j} \; \psi_j + \frac{1}{2} \; \psi_i \; F_{i,j} \; \psi_j \label{eq:helicitymatrix}
!latex      \ee
!latex       where ${\bf a} \equiv \{ A_{\t,e,i,l}, A_{\z,e,i,l}, A_{\t,o,i,l}, A_{\z,o,i,l}, f_{e,i}, f_{o,i} \}$ contains the independent degrees of freedom
!latex       and $\boldpsi \equiv \{\Delta \psi_t,\Delta \psi_p\}$.
!latex \item These are allocated and deallocated in \link{dforce}, assigned in \link{matrix}, and used in \link{mp00ac} and ? \link{df00aa}.
!latex \end{enumerate}


   REAL,   allocatable :: dMA(:,:), dMB(:,:)! dMC(:,:) ! energy and helicity matrices; quadratic forms;
   REAL,   allocatable :: dMD(:,:)! dME(:,:)! dMF(:,:) ! energy and helicity matrices; quadratic forms;

   REAL,   allocatable :: dMAS(:), dMDS(:) ! sparse version of dMA and dMD, data
   INTEGER,allocatable :: idMAS(:), jdMAS(:) ! sparse version of dMA and dMD, indices
   INTEGER,allocatable :: NdMASmax(:), NdMAS(:) ! number of elements for sparse matrices

   REAL,   allocatable :: dMG(:  )

   REAL,   allocatable :: solution(:,:) ! this is allocated in dforce; used in mp00ac and ma02aa; and is passed to packab;

   REAL,   allocatable :: GMRESlastsolution(:,:,:) ! used to store the last solution for restarting GMRES

!  REAL,   allocatable :: MBpsi(:), MEpsi(:) ! matrix vector products;
   REAL,   allocatable :: MBpsi(:)           ! matrix vector products;
!  REAL                :: psiMCpsi, psiMFpsi
!  REAL                ::           psiMFpsi

   LOGICAL             :: LILUprecond        ! whether to use ILU preconditioner for GMRES
   REAL,   allocatable :: BeltramiInverse(:,:)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  REAL   , allocatable :: diotadxup(:,:,:) ! measured rotational transform on inner/outer interfaces for each volume;          d(transform)/dx; (see dforce);
  REAL   , allocatable :: dItGpdxtp(:,:,:) ! measured toroidal and poloidal current on inner/outer interfaces for each volume; d(Itor,Gpol)/dx; (see dforce);

  REAL   , allocatable :: glambda(:,:,:,:) ! save initial guesses for iterative calculation of rotational-transform;

  INTEGER              :: lmns

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{construction of ``force''}
!latex \begin{enumerate}
!latex \item The force vector is comprised of \type{Bomn} and \type{Iomn}.
!latex \end{enumerate}

  REAL,    allocatable ::  Bemn(:,:,:),  Iomn(:,:), Somn(:,:,:), Pomn(:,:,:)
  REAL,    allocatable ::  Bomn(:,:,:),  Iemn(:,:), Semn(:,:,:), Pemn(:,:,:)
  REAL,    allocatable ::  BBe(:), IIo(:), BBo(:), IIe(:)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{\type{Btemn, Bzemn, Btomn, Bzomn} : covariant field on interfaces}
!latex \begin{enumerate}
!latex \item The covariant field:
!latex \end{enumerate}

  REAL,    allocatable ::  Btemn(:,:,:), Bzemn(:,:,:), Btomn(:,:,:), Bzomn(:,:,:) ! covariant components of the tangential field on interfaces;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{\type{Bloweremn, Bloweromn} : covariant field for Hessian computation}
!latex \begin{enumerate}
!latex \item The covariant field:
!latex \end{enumerate}

  REAL,    allocatable ::  Bloweremn(:,:), Bloweromn(:,:)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{\type{LGdof, NGdof} : geometrical degrees-of-freedom;}
!latex \begin{enumerate}
!latex \item The geometrical degrees-of-freedom:
!latex \end{enumerate}

  INTEGER              :: LGdof !       geometrical degrees of freedom associated with each interface;                   ;
  INTEGER              :: NGdof ! total geometrical degrees of freedom                               ;                   ;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{parallel construction of derivative matrix}

!latex \begin{enumerate}
!latex \item The derivatives of force-balance, $[[p+B^2/2]]$, and the spectral constraints (see \link{sw03aa}), with respect to the interface geometry
!latex is constructed in parallel by \link{dforce}.
!latex \item force-balance across the $l$-th interface depends on the fields in the adjacent interfaces.
!latex \end{enumerate}

  REAL,    allocatable :: dBBdRZ(:,:,:)
  REAL,    allocatable :: dIIdRZ(:  ,:)

  REAL,    allocatable :: dFFdRZ(:,:,:,:,:) ! derivatives of B^2 at the interfaces wrt geometry     ;
  REAL,    allocatable :: dBBdmp(:,:,:,:  ) ! derivatives of B^2 at the interfaces wrt mu and dpflux;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{\type{dmupfdx} : derivatives of multiplier and poloidal flux with respect to geometry}

!latex \begin{enumerate}
!latex \item The information in \type{dmupfdx} describes how the helicity multiplier, $\mu$, and the enclosed poloidal flux, $\Delta \psi_p$,
!latex       must vary as the geometry is varied in order to satisfy the interface transform constraint.
!latex \item The internal variable \type{dmupfdx(1:Mvol,1:2,1:LGdof,0:1)} is allocated/deallocated in \link{newton}, and \link{hesian} if selected.
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
!latex \item A finite-difference estimate is computed if \inputvar{Lcheck.eq.4}.
!latex \end{enumerate}

  REAL,    allocatable :: dmupfdx(:,:,:,:,:)  ! derivatives of mu and dpflux wrt geometry at constant interface transform;

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

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{trigonometric factors}

!latex \begin{enumerate}
!latex \item To facilitate construction of the metric integrals, various trigonometric identities are exploited.
!latex \item The required information is saved in
!latex

  REAL   , allocatable :: cosi(:,:), sini(:,:), gteta(:), gzeta(:)

  REAL   , allocatable :: ajk(:)   ! definition of coordinate axis;

  REAL   , allocatable :: dRadR(:,:,:,:), dRadZ(:,:,:,:), dZadR(:,:,:,:), dZadZ(:,:,:,:) ! derivatives of coordinate axis;
  REAL   , allocatable :: dRodR(:,  :,:), dRodZ(:,  :,:), dZodR(:,  :,:), dZodZ(:,  :,:) ! derivatives of coordinate axis;

  INTEGER, allocatable :: djkp(:,:), djkm(:,:) ! for calculating cylindrical volume;

!latex \item The following are used for volume integrals (see \link{volume})
!latex \be a_{i,j,k} &=& 4 \; m_k \ooint \cos(\alpha_i)\cos(\alpha_j)\cos(\alpha_k) /(2\pi)^2 , \\
!latex     b_{i,j,k} &=& 4 \; m_j \ooint \cos(\alpha_i)\sin(\alpha_j)\sin(\alpha_k) /(2\pi)^2 ,
!latex \ee
!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{\type{lBBintegral, lABintegral} : volume integrals}
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
!latex \type{dFF(1:Nvol,0:1,0:mn+mn-1)}, where
!latex
!latex $         F_l                      \equiv$ \type{dFF(l,0,    0)}
!latex
!latex $\partial F_l / \partial R_{l-1,j} \equiv$ \verb&dFF(ll,0,   j)&
!latex
!latex $\partial F_l / \partial Z_{l-1,j} \equiv$ \verb&dFF(ll,0,mn}j)&
!latex
!latex $\partial F_l / \partial R_{l  ,j} \equiv$ \type{dFF(ll,1,   j)}
!latex
!latex $\partial F_l / \partial Z_{l  ,j} \equiv$ \verb&dFF(ll,1,mn}j)&
!latex

!latex \item The volume integrals $\int dv$, $\int B^2 \; dv$ and $\int {\bf A}\cdot{\bf B} \; dv$ in each volume
!latex       are computed and saved in \type{volume(0:2,1:Nvol)}.
!latex \end{enumerate}

  REAL   , allocatable :: lBBintegral(:) ! B.B      integral;
  REAL   , allocatable :: lABintegral(:) ! A.B      integral;

! REAL                 :: dBBintegral    ! B.B      integral; derivative wrt R, Z;

! REAL   , allocatable :: oBBintegral(:) ! B.B      integral; original; used to normalize; perhaps irrelevant;


! REAL                 :: dABintegral    ! A.B      integral; derivative wrt R, Z;

  REAL   , allocatable :: vvolume(:) ! volume integral of \sqrt g; computed in volume;
  REAL                 :: dvolume    ! derivative of volume wrt interface geometry;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! internal global variables; internal logical variables; default values are provided here; these may be changed according to input values;

  INTEGER              :: ivol ! labels volume; some subroutines (called by NAG) are fixed argument list but require the volume label;

  REAL                 :: gBzeta ! toroidal (contravariant) field; calculated in bfield; required to convert \dot \t to B^\t, \dot s to B^s;

  INTEGER, allocatable :: Iquad(:) ! internal copy of Nquad;

  REAL   , allocatable :: gaussianweight(:,:), gaussianabscissae(:,:)

  LOGICAL              :: LBlinear, LBnewton, LBsequad ! controls selection of Beltrami field solver; depends on LBeltrami;

  REAL                 :: oRZp(1:3) ! used in mg00aa to determine (\s,\t,\z) given (R,Z,p);

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  type derivative
     LOGICAL :: L
     INTEGER :: vol      ! Used in coords.f90; required for global constraint force gradient evaluation
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

subroutine build_vector_potential(lvol, iocons, aderiv, tderiv)

! Builds the covariant component of the vector potential and store them in efmn, ofmn, sfmn, cfmn.

  use constants, only: zero, half

  use fileunits, only: ounit

  use inputlist, only: Lrad, Wbuild_vector_potential, Wmacros

  use cputiming

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS

  INTEGER              :: aderiv    ! Derivative of A. -1: w.r.t geometrical degree of freedom
                                    !                   0: no derivatives
                                    !                   1: w.r.t mu
                                    !                   2: w.r.t pflux
  INTEGER              :: tderiv    ! Derivative of Chebyshev polynomialc. 0: no derivatives
                                    !                                      1: w.r.t radial coordinate s
  INTEGER              :: ii,  &    ! Loop index on Fourier harmonics
                          ll,  &    ! Loop index on radial resolution
                          mi,  &    ! Poloidal mode number
                          lvol,&    ! Volume number
                          iocons    ! inner (0) or outer (1) side of the volume
  REAL                 :: mfactor   ! Regularization factor when LcoordinateSingularity

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  BEGIN(build_vector_potential)

  efmn(1:mn) = zero ; sfmn(1:mn) = zero ; cfmn(1:mn) = zero ; ofmn(1:mn) = zero

  do ii = 1, mn ! loop over Fourier harmonics; 13 Sep 13;

   if( Lcoordinatesingularity ) then
    mi = im(ii)
    do ll = mi, Lrad(lvol),2 ! loop over Zernike polynomials; Lrad is the radial resolution; 01 Jul 19;
      ;                      ; efmn(ii) = efmn(ii) +          Ate(lvol,aderiv,ii)%s(ll) * RTT(ll,mi,iocons,1) * half
      ;                      ; cfmn(ii) = cfmn(ii) +          Aze(lvol,aderiv,ii)%s(ll) * RTT(ll,mi,iocons,1) * half
      if( NOTstellsym ) then ; ofmn(ii) = ofmn(ii) +          Ato(lvol,aderiv,ii)%s(ll) * RTT(ll,mi,iocons,1) * half
      ;                      ; sfmn(ii) = sfmn(ii) +          Azo(lvol,aderiv,ii)%s(ll) * RTT(ll,mi,iocons,1) * half
      endif
    enddo ! end of do ll; 20 Feb 13;
   else
    do ll = 0, Lrad(lvol) ! loop over Chebyshev polynomials; Lrad is the radial resolution;
      ;                      ; efmn(ii) = efmn(ii) +          Ate(lvol,aderiv,ii)%s(ll) * TT(ll,iocons,1) ! aderiv labels deriv. wrt mu, pflux;
      ;                      ; cfmn(ii) = cfmn(ii) +          Aze(lvol,aderiv,ii)%s(ll) * TT(ll,iocons,1)
      if( NOTstellsym ) then ; ofmn(ii) = ofmn(ii) +          Ato(lvol,aderiv,ii)%s(ll) * TT(ll,iocons,1)
      ;                     ; sfmn(ii) = sfmn(ii) +          Azo(lvol,aderiv,ii)%s(ll) * TT(ll,iocons,1)
      endif
    enddo ! end of do ll; 20 Feb 13;
   end if ! Lcoordinatesingularity; 01 Jul 19;
  enddo ! end of do ii; 20 Feb 13;

end subroutine build_vector_potential

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
  INTEGER              :: vvol, mm, nn, nb, imn, ix, ii, jj, ij, kk, mj, nj, mk, nk, ip, X02BBF, iargc, iarg, numargs, mi, ni, lvol, extlen, sppos
  REAL                 :: xx, toroidalflux, toroidalcurrent
  REAL,    allocatable :: RZRZ(:,:) ! local array used for reading interface Fourier harmonics from file;

  CHARACTER            :: ldate*8, ltime*10, arg*100

  BEGIN(readin)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{\type{machprec}, \type{vsmall}, \type{small}, \type{sqrtmachprec} : machine precision}

!latex \begin{enumerate}
!latex \item The machine precision is determined using the Fortran 90 intrinsic function EPSILON.
!latex \end{enumerate}

  cput = GETTIME

  machprec = myprec() ! is this required? Why not just set real, parameter :: machprec = 1.0E-16 ? ; let's simplify the source; SRH: 27 Feb 18;

  vsmall = 100*machprec ; small = 100*vsmall ; sqrtmachprec = sqrt(machprec) ! returns machine precision;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( myid.eq.0 ) then ! only the master node reads input file and sets secondary variables;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   call date_and_time( ldate, ltime )

   write(ounit,'("readin : ", 10x ," : ")')
   write(ounit,1000) cput-cpus, ldate(1:4), ldate(5:6), ldate(7:8), ltime(1:2), ltime(3:4), ltime(5:6), machprec, vsmall, small

1000 format("readin : ",f10.2," : date="a4"/"a2"/"a2" , "a2":"a2":"a2" ; machine precision="es9.2" ; vsmall="es9.2" ; small="es9.2" ;")

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : reading ext from command line ;")') cput-cpus
   endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   inquire( file=trim(ext)//".sp", exist=Lspexist ) ! check if file exists;
   FATAL( readin, .not.Lspexist, the input file does not exist ) ! if not, abort;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   open( iunit, file=trim(ext)//".sp", status="old")

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

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

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : reading physicslist     from ext.sp ;")') cput-cpus
   endif

   read(iunit,physicslist)

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : read    physicslist     from ext.sp ;")') cput-cpus
   endif

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

        Ivolume(Mvol) = Ivolume(Mvol-1) !Ensure vacuum in vacuum region

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

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : reading numericlist     from ext.sp ;")') cput-cpus
   endif

   read(iunit,numericlist)!,iostat=ios)

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : read    numericlist     from ext.sp ;")') cput-cpus
   endif

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

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : reading locallist      from ext.sp ;")') cput-cpus
   endif

   read(iunit,locallist)!,iostat=ios)

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : read    locallist      from ext.sp ;")') cput-cpus
   endif

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
1044 format("readin : ", 10x ," : vcasingeps="es13.5" ; vcasingtol="es13.5" ; vcasingits="i6" ; vcasingper="i6" ;")

   FATAL( readin, escale      .lt.zero     , error )
   FATAL( readin, pcondense   .lt.one      , error )
   FATAL( readin, abs(c05xtol).le.machprec , error )
   FATAL( readin, c05factor   .le.zero     , error )
  !FATAL( readin, mfreeits    .lt.zero     , error )

   FATAL( readin, Igeometry.eq.3 .and. pcondense.le.zero, pcondense must be positive )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{reading \type{diagnosticslist}}

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : reading diagnosticslist from ext.sp ;")') cput-cpus
   endif

   read(iunit,diagnosticslist)!,iostat=ios)

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : read    diagnosticslist from ext.sp ;")') cput-cpus
   endif

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

!latex \subsubsection{\type{Mvol} : total number of volumes}
!latex \begin{enumerate}
!latex \item The number of plasma volumes is \internal{Mvol}=\inputvar{Nvol}+\inputvar{Lfreebound};
!latex \end{enumerate}

  FATAL( readin, Lfreebound.lt.0 .or. Lfreebound.gt.1, illegal Lfreebound )

  Mvol = Nvol + Lfreebound

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  SALLOCATE( beltramierror,(1:Mvol,1:9), zero)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{\type{mn}, \type{im(1:mn)} and \type{in(1:mn)} : Fourier mode identification}
!latex \begin{enumerate}
!latex \item The Fourier description of even periodic functions is
!latex       \be f(\t,\z) = \sum_{n=0}^{N} f_{0,n} \cos(-n\z) + \sum_{m=1}^{M}\sum_{n=-N}^{N} f_{m,n} \cos(m\t-n\z),
!latex       \ee
!latex       where the resolution is given on input, $M\equiv $ \inputvar{ Mpol} and $N\equiv $ \inputvar{ Ntor}.
!latex \item For convenience, the Fourier summations are written as
!latex       \be f(\s,\t,\z) &=& \sum_j f_j(s) \cos( m_j \t - n_j \z ),
!latex       \ee
!latex       for $j=1,$ \type{mn}, where \type{mn}$ = N + 1 +  M  ( 2 N + 1 )$.
!latex \item The integer arrays \type{im(1:mn)} and \type{in(1:mn)} contain the $m_j$ and $n_j$.
!latex \item The array \type{in} includes the \type{Nfp} factor.
!latex \end{enumerate}

  mn = 1 + Ntor +  Mpol * ( 2 *  Ntor + 1 ) ! Fourier resolution of interface geometry & vector potential;

  SALLOCATE( im, (1:mn), 0 )
  SALLOCATE( in, (1:mn), 0 )

  call gi00ab(  Mpol,  Ntor, Nfp, mn, im(1:mn), in(1:mn) ) ! this sets the im and in mode identification arrays;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{\type{halfmm(1:mn}, regumm(1:mn) : regularization factor}
!latex \begin{enumerate}
!latex \item The ``regularization'' factor, \type{halfmm(1:mn)} = \type{im(1:mn)} * \type{half}, is real.
!latex \item This is used in \link{lforce}, \link{bfield}, \link{stzxyz}, \link{coords}, \link{jo00aa}, \link{ma00aa}, \link{sc00aa} and \link{tr00ab}.
!latex \end{enumerate}

  SALLOCATE( halfmm, (1:mn), im(1:mn) * half )
  SALLOCATE( regumm, (1:mn), im(1:mn) * half )

  if( Mregular.ge.2 ) then

   where( im.gt.Mregular ) regumm = Mregular * half

  endif

! if( myid.eq.0 ) write(ounit,'("global : " 10x " : "i3") im ="i3" , halfmm ="f5.1" , regum ="f5.1" ;")') ( ii, im(ii), halfmm(ii), regumm(ii), ii = 1, mn )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{\type{ime} and \type{ine} : extended resolution Fourier mode identification}
!latex \begin{enumerate}
!latex \item The ``extended'' Fourier resolution is defined by \internal{lMpol} $ = 4 $ \inputvar{Mpol}, \internal{lNtor} $ = 4 $\inputvar{Ntor}.
!latex \end{enumerate}

! lMpol =   Mpol ; lNtor =   Ntor ! no    enhanced resolution for metrics;
! lMpol = 2*Mpol ; lNtor = 2*Ntor !       enhanced resolution for metrics;
  lMpol = 4*Mpol ; lNtor = 4*Ntor ! extra-enhanced resolution for metrics;

  mne = 1 + lNtor + lMpol * ( 2 * lNtor + 1 ) ! resolution of metrics; enhanced resolution; see metrix;

  SALLOCATE( ime, (1:mne), 0 )
  SALLOCATE( ine, (1:mne), 0 )

  call gi00ab( lMpol, lNtor, Nfp, mne, ime(1:mne), ine(1:mne) )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{\type{mns}, \type{ims} and \type{ins} : Fourier mode identification for straight-fieldline angle}

  sMpol = iMpol ; sNtor = iNtor

  if( iMpol.le.0 ) sMpol = Mpol - iMpol
  if( iNtor.le.0 ) sNtor = Ntor - iNtor
  if(  Ntor.eq.0 ) sNtor = 0

  mns = 1 + sNtor + sMpol * ( 2 * sNtor + 1 ) ! resolution of straight-field line transformation on interfaces; see tr00ab; soon to be redundant;

  SALLOCATE( ims, (1:mns), 0 )
  SALLOCATE( ins, (1:mns), 0 )

  call gi00ab( sMpol, sNtor, Nfp, mns, ims(1:mns), ins(1:mns) ) ! note that the field periodicity factor is included in ins;

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

!latex \subsubsection{\type{iRbc(1:mn,0:Mvol}, \type{iZbs(1:mn,0:Mvol}, \type{iRbs(1:mn,0:Mvol} and \type{iZbc(1:mn,0:Mvol} : geometry}

!latex \begin{enumerate}
!latex \item \type{iRbc}, \type{iZbs}, \type{iRbs} and \type{iZbc} : Fourier harmonics of interface geometry;
!latex \item \type{iVns}, \type{iVnc}, \type{iBns} and \type{iBns} : Fourier harmonics of normal field at computational boundary;
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

!latex \subsubsection{\type{ajk} : construction of coordinate axis}

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

  if ( Igeometry.eq.1 .and. Lreflect.eq.1) then ! reflect upper and lower bound in slab, each take half the amplitude
    iRbc(2:mn,Mvol) = iRbc(2:mn,Mvol) * half
    iRbc(2:mn,0) = -iRbc(2:mn,Mvol)
   if( NOTstellsym ) then
    iRbs(2:mn,Mvol) = iRbs(2:mn,Mvol) * half
    iRbs(2:mn,0) = -iRbs(2:mn,Mvol)
   endif
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  Rscale = iRbc(1,Mvol) ! this will be used to normalize the geometrical degrees-of-freedom;

  if( myid.eq.0 ) write(ounit,'("readin : ", 10x ," : myid=",i3," ; Rscale=",es22.15," ;")') myid, Rscale

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(readin)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine readin

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

  INTEGER              :: vvol, imn, ii, jj, kk, jk, Lcurvature, mm, nn
  REAL                 :: lss, teta, zeta, st(1:Node), Bst(1:Node), BR, BZ, BP

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
  write(iunit,'(" Ivolume     = ",257es23.15)') Ivolume(1:Mvol)
  write(iunit,'(" Isurf       = ",257es23.15)') Isurf(1:Mvol-1)
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
  write(iunit,'(" Lreflect    = ",i9        )') Lreflect
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
     if( nn.lt.  0 .and. nn.ge.- 9 )     write(iunit,1001) nn, mm, Rbc(nn,mm), nn, mm, Zbs(nn,mm), nn, mm, Rbs(nn,mm), nn, mm, Zbc(nn,mm)
     if( nn.ge.  0 .and. nn.le.  9 )     write(iunit,1002) nn, mm, Rbc(nn,mm), nn, mm, Zbs(nn,mm), nn, mm, Rbs(nn,mm), nn, mm, Zbc(nn,mm)
     if( nn.gt.  9 .and. nn.le. 99 )     write(iunit,1001) nn, mm, Rbc(nn,mm), nn, mm, Zbs(nn,mm), nn, mm, Rbs(nn,mm), nn, mm, Zbc(nn,mm)
    case(  10: 99 )
     if( nn.lt.- 9 .and. nn.gt.-99 )     write(iunit,1003) nn, mm, Rbc(nn,mm), nn, mm, Zbs(nn,mm), nn, mm, Rbs(nn,mm), nn, mm, Zbc(nn,mm)
     if( nn.lt.  0 .and. nn.ge.- 9 )     write(iunit,1004) nn, mm, Rbc(nn,mm), nn, mm, Zbs(nn,mm), nn, mm, Rbs(nn,mm), nn, mm, Zbc(nn,mm)
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
     if( nn.lt.  0 .and. nn.ge.- 9 ) write(iunit,1011) nn, mm, Rwc(nn,mm), nn, mm, Zws(nn,mm), nn, mm, Rws(nn,mm), nn, mm, Zwc(nn,mm)
     if( nn.ge.  0 .and. nn.le.  9 ) write(iunit,1012) nn, mm, Rwc(nn,mm), nn, mm, Zws(nn,mm), nn, mm, Rws(nn,mm), nn, mm, Zwc(nn,mm)
     if( nn.gt.  9 .and. nn.le. 99 ) write(iunit,1011) nn, mm, Rwc(nn,mm), nn, mm, Zws(nn,mm), nn, mm, Rws(nn,mm), nn, mm, Zwc(nn,mm)
    case(  10: 99 )
     if( nn.lt.- 9 .and. nn.gt.-99 ) write(iunit,1013) nn, mm, Rwc(nn,mm), nn, mm, Zws(nn,mm), nn, mm, Rws(nn,mm), nn, mm, Zwc(nn,mm)
     if( nn.lt.  0 .and. nn.ge.- 9 ) write(iunit,1014) nn, mm, Rwc(nn,mm), nn, mm, Zws(nn,mm), nn, mm, Rws(nn,mm), nn, mm, Zwc(nn,mm)
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
     if( nn.lt.  0 .and. nn.ge.- 9 ) write(iunit,1021) nn, mm, Vns(nn,mm), nn, mm, Bns(nn,mm), nn, mm, Vnc(nn,mm), nn, mm, Bnc(nn,mm)
     if( nn.ge.  0 .and. nn.le.  9 ) write(iunit,1022) nn, mm, Vns(nn,mm), nn, mm, Bns(nn,mm), nn, mm, Vnc(nn,mm), nn, mm, Bnc(nn,mm)
     if( nn.gt.  9 .and. nn.le. 99 ) write(iunit,1021) nn, mm, Vns(nn,mm), nn, mm, Bns(nn,mm), nn, mm, Vnc(nn,mm), nn, mm, Bnc(nn,mm)
    case(  10: 99 )
     if( nn.lt.- 9 .and. nn.gt.-99 ) write(iunit,1023) nn, mm, Vns(nn,mm), nn, mm, Bns(nn,mm), nn, mm, Vnc(nn,mm), nn, mm, Bnc(nn,mm)
     if( nn.lt.  0 .and. nn.ge.- 9 ) write(iunit,1024) nn, mm, Vns(nn,mm), nn, mm, Bns(nn,mm), nn, mm, Vnc(nn,mm), nn, mm, Bnc(nn,mm)
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
  write(iunit,'(" Lrzaxis     = ",i9            )') Lrzaxis
  write(iunit,'(" Ntoraxis    = ",i9            )') Ntoraxis
  write(iunit,'("/")')

  if( Wwrtend ) then ; cput = GETTIME ; write(ounit,'("wrtend : ",f10.2," : myid=",i3," ; writing locallist ;")') cput-cpus, myid
  endif

  write(iunit,'("&locallist")')
  write(iunit,'(" LBeltrami   = ",i9            )') LBeltrami
  write(iunit,'(" Linitgues   = ",i9            )') Linitgues
  write(iunit,'(" Lmatsolver  = ",i9            )') Lmatsolver
  write(iunit,'(" NiterGMRES  = ",i9            )') NiterGMRES
  write(iunit,'(" LGMRESprec  = ",i9            )') LGMRESprec
  write(iunit,'(" epsGMRES    = ",es23.15       )') epsGMRES
  write(iunit,'(" epsILU      = ",es23.15       )') epsILU

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
  write(iunit,'(" Ppts        = ",es23.15       )') Ppts
  write(iunit,'(" nPtrj       = ",256i6         )') nPtrj(1:Mvol)
  write(iunit,'(" LHevalues   = ",L9            )') LHevalues
  write(iunit,'(" LHevectors  = ",L9            )') LHevectors
  write(iunit,'(" LHmatrix    = ",L9            )') LHmatrix
  write(iunit,'(" Lperturbed  = ",i9            )') Lperturbed
  write(iunit,'(" dpp         = ",i9            )') dpp
  write(iunit,'(" dqq         = ",i9            )') dqq
  write(iunit,'(" dRZ         = ",es23.15       )') dRZ
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

  ! write initial guess of interface geometry
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

subroutine IsMyVolume(vvol)

!latex \subsection{subroutine IsMyVolume}
!latex Check if volume vvol is associated to the corresponding MPI node.

LOCALS

INTEGER :: vvol
INTEGER :: lwbnd, upbnd

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

IsMyVolumeValue = -1 ! Error value - Problem with vvol / id
if( myid.ne.modulo(vvol-1,ncpu) ) then
  IsMyVolumeValue = 0
else
  IsMyVolumeValue = 1
endif

end subroutine IsMyVolume

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine WhichCpuID(vvol, cpu_id)
!latex \subsection{subroutine WhichCpuID}
!latex Returns which MPI node is associated to a given volume.

LOCALS

INTEGER            :: vvol, cpu_id

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

cpu_id = modulo(vvol-1,ncpu)

end subroutine WhichCpuID

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end module allglobal

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

module fftw_interface ! JAB; 25 Jul 17

  use, intrinsic :: iso_c_binding

  implicit none

  include 'fftw3.f03'

  TYPE(C_PTR)                            :: planf, planb
  COMPLEX(C_DOUBLE_COMPLEX), allocatable :: cplxin(:,:,:), cplxout(:,:,:)

end module fftw_interface

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
