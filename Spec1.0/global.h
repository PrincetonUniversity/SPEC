!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
!        1         2         3         4         5         6         7         8         9        10        11        12        13        14        15        
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Defines input, global variables and output files, and provides compilation instructions.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! to set keyboard shortcut in emacs                                        

! (1) define macro         , e.g. \C-x \C-( . . . \C-x \C-)              
! (2) name macro           , e.g. Esc-x name-last-kbd-macro arbitraryname ! 11 Oct 12; 
! (3) set keyboard shortcut, e.g. Esc-x global-set-key F12 arbitraryname 

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

module constants
  implicit none

  REAL, parameter :: zero       =   0.0
  REAL, parameter :: one        =   1.0
  REAL, parameter :: two        =   2.0
  REAL, parameter :: three      =   3.0
  REAL, parameter :: four       =   4.0
  REAL, parameter :: five       =   5.0
  REAL, parameter :: six        =   6.0
  REAL, parameter :: seven      =   7.0
  REAL, parameter :: eight      =   8.0
  REAL, parameter :: nine       =   9.0
  REAL, parameter :: ten        =  10.0
  REAL, parameter :: eleven     =  11.0
  REAL, parameter :: twelve     =  12.0

  REAL, parameter :: half       =  one / two
  REAL, parameter :: third      =  one / three 
  REAL, parameter :: quart      =  one / four
  REAL, parameter :: fifth      =  one / five
  REAL, parameter :: sixth      =  one / six

  REAL, parameter :: pi2        =  6.28318530717958623
  REAL, parameter :: pi         =  pi2 / two
  REAL, parameter :: mu0        =  2.0E-07 * pi2
  REAL, parameter :: goldenmean =  1.618033988749895 ! golden mean = ( one + sqrt(five) ) / two ;

end module constants

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

module numerical
  implicit none
  INTEGER         :: largestint
  REAL            :: machprec, vsmall, small, sqrtmachprec ! these are assigned below in readin via a call to NAG routine;
  REAL, parameter :: logtolerance = 1.0e-32 
end module numerical

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

module fileunits
  implicit none
  INTEGER :: aunit = 11 ! vector potential; used in ra00aa:.ext.AtAzmn; 11 Aug 14;
  INTEGER :: dunit = 12 ! derivative matrix; used in jk03aa:.ext.GF; 11 Aug 14;
  INTEGER :: gunit = 13 ! mgrid file unit; used in mg00aa:mgrid, wa00aa:wall.dat; 11 Aug 14;
  INTEGER :: hunit = 14 ! eigenvalues of Hessian; under re-construction; 11 Aug 14;
  INTEGER :: iunit = 10 ! input; used in global/readin:ext.spec, global/writin:ext.end, global/writin:.ext.grid; 11 Aug 14;
  INTEGER :: lunit = 20 ! local unit; used in lunit+myid: pp00aa:.ext.poincare,.ext.transform; 11 Aug 14;
  INTEGER :: ounit =  0 ! screen output;
  INTEGER :: vunit = 15 ! for examination of adaptive quadrature; used in vc00aa:.ext.vcint; 11 Aug 14;
  INTEGER :: wunit = 16 ! for LaTeX formatted output;
  INTEGER :: zunit = 17 ! for convergence; this file is opened in xspech:.ext.iterations, and written to in globals/writin; 11 Aug 14;
 !INTEGER :: wunit = 15 ! irrational surfaces; under re-construction; 11 Aug 14;
 !INTEGER :: qunit = 13 ! quadratic-flux minimizing surfaces; used in pq02aa; 11 Aug 14;
 !INTEGER :: munit = 14 ! periodic orbits; used in hr00aa; 11 Aug 14;
 !INTEGER :: funit = 16 ! force iterations;
 !INTEGER :: eunit = 21 ! for comparison with extender; used in 11 Aug 14;
 !INTEGER :: bunit = 23 ! for boundary harmonics;
end module fileunits

module cputiming
! CPUVARIABLE ! this is expanded by Makefile; do not remove;
  REAL :: Treadin=0.0
  REAL :: Twritin=0.0
end module cputiming

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

module typedefns
  
  type subgrid
     REAL,    allocatable :: s(:)
     INTEGER, allocatable :: i(:)
  end type subgrid

end module typedefns

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

module inputlist !latex \end{enumerate} \subsection{Input namelists} \begin{enumerate}
  implicit none

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  REAL               :: version = 2.00

!latex \item The input file, \verb+ext.spec+, where \verb+ext*100+ is given as command line input, contains the following namelists and interface geometry.

!SET RESOLUTION;

  CHARACTER          :: ext*100

  INTEGER, parameter :: MNvol     = 256 !latex \item The maximum value of \verb+Nvol+ is \verb+MNvol=256+.
  INTEGER, parameter :: MMpol     =  32 !latex \item The maximum value of \verb+Mpol+ is \verb+MNpol= 32+.
  INTEGER, parameter :: MNtor     =  32 !latex \item The maximum value of \verb+Ntor+ is \verb+MNtor= 32+.
! INTEGER, parameter :: Mmn       = MNtor+1 + MMpol*(2*MNtor+1)
  INTEGER, parameter :: maxgroups = 100 !latex \item The maximum number of independent coil currents is given by \verb+maxgroups=100+.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! the following variables constitute the namelist/physicslist/; note that all variables in namelist need to be broadcasted in readin;

  INTEGER      :: Igeometry                  =  3
  INTEGER      :: pknot                      =  1
  INTEGER      :: Istellsym                  =  1
  INTEGER      :: Lfreebound                 =  0
  REAL         :: phiedge                    =  1.0
  REAL         :: curtor                     =  0.0
  REAL         :: curpol                     =  0.0
  REAL         :: extcur(1:maxgroups)        =  0.0
  REAL         :: gamma                      =  0.0
  INTEGER      :: Nfp                        =  1
  INTEGER      :: Nvol                       =  1
  INTEGER      :: Mpol                       =  0
  INTEGER      :: Ntor                       =  0
  INTEGER      :: Lrad(1:MNvol+1)            =  2
  INTEGER      :: Lconstraint                =  2
  REAL         ::     tflux(1:MNvol)         =  0.0
  REAL         ::     pflux(1:MNvol)         =  0.0
  REAL         ::  helicity(1:MNvol)         =  0.0
  REAL         :: pscale                     =  0.0
  REAL         ::  pressure(1:MNvol)         =  0.0
  INTEGER      :: Ladiabatic                 =  0
  REAL         :: adiabatic(1:MNvol)         =  0.0
  REAL         ::        mu(1:MNvol)         =  0.0
  INTEGER      ::        pl(0:MNvol)         =  0
  INTEGER      ::        ql(0:MNvol)         =  0
  INTEGER      ::        pr(0:MNvol)         =  0
  INTEGER      ::        qr(0:MNvol)         =  0
  REAL         ::      iota(0:MNvol)         =  0.0
  INTEGER      ::        lp(0:MNvol)         =  0
  INTEGER      ::        lq(0:MNvol)         =  0
  INTEGER      ::        rp(0:MNvol)         =  0
  INTEGER      ::        rq(0:MNvol)         =  0
  REAL         ::      oita(0:MNvol)         =  0.0

  REAL         :: Rac(     0:MNtor        )  =  0.0 !     stellarator symmetric coordinate axis; 13 Sep 13;
  REAL         :: Zas(     0:MNtor        )  =  0.0
  REAL         :: Ras(     0:MNtor        )  =  0.0 ! non-stellarator symmetric coordinate axis; 13 Sep 13;
  REAL         :: Zac(     0:MNtor        )  =  0.0

  REAL         :: Rbc(-MNtor:MNtor,-MMpol:MMpol)  =  0.0 !     stellarator symmetric boundary components;  20 Jan 2016; 
  REAL         :: Zbs(-MNtor:MNtor,-MMpol:MMpol)  =  0.0 !     stellarator symmetric boundary components;
  REAL         :: Rbs(-MNtor:MNtor,-MMpol:MMpol)  =  0.0 ! non-stellarator symmetric boundary components;
  REAL         :: Zbc(-MNtor:MNtor,-MMpol:MMpol)  =  0.0 ! non-stellarator symmetric boundary components;

  REAL         :: Rwc(-MNtor:MNtor,-MMpol:MMpol)  =  0.0 !     stellarator symmetric boundary components of wall;
  REAL         :: Zws(-MNtor:MNtor,-MMpol:MMpol)  =  0.0 !     stellarator symmetric boundary components of wall;
  REAL         :: Rws(-MNtor:MNtor,-MMpol:MMpol)  =  0.0 ! non-stellarator symmetric boundary components of wall;
  REAL         :: Zwc(-MNtor:MNtor,-MMpol:MMpol)  =  0.0 ! non-stellarator symmetric boundary components of wall;

  REAL         :: Bns(-MNtor:MNtor,-MMpol:MMpol)  =  0.0 !     stellarator symmetric normal field at boundary;
  REAL         :: Bnc(-MNtor:MNtor,-MMpol:MMpol)  =  0.0 ! non-stellarator symmetric normal field at boundary;
 !REAL         :: Vns(-MNtor:MNtor,-MMpol:MMpol)  =  0.0 !     stellarator symmetric normal field at boundary; vacuum component;
 !REAL         :: Vnc(-MNtor:MNtor,-MMpol:MMpol)  =  0.0 ! non-stellarator symmetric normal field at boundary; vacuum component;

  REAL         :: mupftol                    =  1.0e-08
  INTEGER      :: mupfits                    =  8

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! the following variables constitute the namelist/numericlist/; note that all variables in namelist need to be broadcasted in readin;

  INTEGER :: Linitialize =  0 
  INTEGER :: Lwall       =  0 
  REAL    :: phiwall     =  0.9
  INTEGER :: Ndiscrete   =  2
  INTEGER :: Nquad       = -1
  INTEGER :: iMpol       = -4
  INTEGER :: iNtor       = -4
  INTEGER :: Lsparse     =  0
  INTEGER :: Lsvdiota    =  0
  INTEGER :: imethod     =  3
  INTEGER :: iorder      =  2
  INTEGER :: iprecon     =  0
  REAL    :: iotatol     = -1.0
  INTEGER :: Iswmin      =  0
  INTEGER :: Lperturb    =  0 
  REAL    :: dperturb    =  1.0e-04
  INTEGER :: Lextrap     =  0

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! the following variables constitute the namelist/locallist/; note that all variables in namelist need to be broadcasted in readin;

  INTEGER :: LBeltrami  =  4
  INTEGER :: Linitgues  =  1
  LOGICAL :: Lposdef    = .true.
  INTEGER :: Nmaxexp    = 32
! INTEGER :: Lprecon    =  1
! INTEGER :: sparseits  = -1
! REAL    :: ssoromega  =  1.0
! REAL    :: sparsetol  = -1.0
! INTEGER :: Liotasolv  =  1

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! the following variables constitute the namelist/globallist/; note that all variables in namelist need to be broadcasted in readin;

! LOGICAL   :: Lvacuum    = .false.
  LOGICAL   :: Le04dgf    = .false. ! redundant; 10 Oct 12;
  INTEGER   :: Lminimize  =  0
  INTEGER   :: Lfindzero  =  0
! INTEGER   :: pwidth     =  4 ! redundant; 26 Feb 13;
! INTEGER   :: qwidth     =  4 ! redundant; 26 Feb 13;
  INTEGER   :: Lcondense  =  1 ! redundant; 18 Jul 14;
  REAL      :: escale     =  0.0
  REAL      :: pcondense  =  2.0
  REAL      :: qcondense  = -1.0 ! 04 Dec 14; redundant; 
  REAL      :: wpoloidal  =  1.0 ! exponential weight factor on poloidal length constraint; 04 Dec 14;
  REAL      :: wcondense  =  1.0 ! redundant;
  REAL      :: forcetol   =  1.0e-10
  REAL      :: normalerr  =  1.0e-06
  REAL      :: norblend   =  0.0
  INTEGER   :: maxfbits   = 20
  REAL      :: ForceErr   = -1.0 ! this is really an output; redundant; 04 Dec 14;
  REAL      :: Energy     =  0.0 ! redundant; 02 Nov 12;
  LOGICAL   :: Le04lyf    = .false. ! redundant; 10 Oct 12;
  REAL      :: c05xtol    =  1.0e-12
  REAL      :: factor     = -9.0
  REAL      :: c05factor  =  1.0e-02
  LOGICAL   :: LreadGF    = .true.
  INTEGER   :: verify     = -1
  REAL      :: maxstep    =  1.0e-03
  REAL      :: opsilon    =  1.0e-00 ! weight factor for pressure imbalance; 04 Dec 14;
  REAL      :: epsilon    =  1.0e-00 ! weight factor for spectral constraints; 04 Dec 14;
  REAL      :: upsilon    =  1.0e-00 ! weight factor for poloidal length constraint; 04 Dec 14;
  REAL      :: apsilon    =  1.0e-00 ! redundant; 18 Jul 14;
  INTEGER   :: maxiter    = -1

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! the following variables constitute the namelist/diagnosticslist/; note that all variables in namelist need to be broadcasted in readin;

  REAL    :: odetol           =     1.0e-07
  REAL    :: absreq           =     1.0e-08
  REAL    :: relreq           =     1.0e-08
  REAL    :: absacc           =     1.0e-04
  REAL    :: epsr             =     1.0e-08
! REAL    :: divertorR        =    -1.0 ! redundant; 10 Oct 12;
! REAL    :: divertorZ        =     1.0 ! redundant; 10 Oct 12;
  INTEGER :: nPpts            =     0
  INTEGER :: nPtrj(1:MNvol+1) =    -1
  INTEGER :: Mpqits           =    10
  INTEGER :: Lpqsym           =     1
  INTEGER :: p1(1:MNvol+1)    =     0
  INTEGER :: q1(1:MNvol+1)    =     0
  INTEGER :: p2(1:MNvol+1)    =     0
  INTEGER :: q2(1:MNvol+1)    =     0
  REAL    :: pqs(1:MNvol)     =    -2.0
  REAL    :: pqt(1:MNvol)     =     0.0
! REAL    :: pqR(1:MNvol)     =    -1.0
! REAL    :: pqZ(1:MNvol)     =    -1.0
  INTEGER :: Npl              =     1
! INTEGER :: Munstable        =    10  
! REAL    :: dunstable        =     1.0e-05
! INTEGER :: Nunstable        =  1000
  INTEGER :: npq(1:MNvol+1)   =     0
! INTEGER :: Mirrits          =     0
! INTEGER :: irrMpol          =    50
! INTEGER :: irrNtor          =    25
! REAL    :: irrsvdcut        =     1.0e-12
! REAL    :: irrtol           =     1.0e-08
  LOGICAL :: Lwrpj            =  .false. ! redundant; 02 Sep 14;
! INTEGER :: Nghd             =     0
  LOGICAL :: LHessian         =  .false. ! redundant; 20 Jun 14;
  LOGICAL :: LHevalues        =  .false.
  LOGICAL :: LHevectors       =  .false.
  INTEGER :: Lperturbed       =     0   
  INTEGER :: dpp              =    -1
  INTEGER :: dqq              =    -1
  LOGICAL :: Lcurlerr         =  .false. ! redundant; 01 Jul 14;
  INTEGER :: Lcheck           =     0
  LOGICAL :: Ltiming          =  .false.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! the following variables constitute the namelist/screenlist/; note that all variables in namelist need to be broadcasted in readin;
  
! DSCREENLIST ! define screenlist; this is expanded by Makefile; DO NOT REMOVE; each file compiled by Makefile has its own write flag;
  LOGICAL :: Wreadin = .false.
  LOGICAL :: Wwritin = .false.
  LOGICAL :: Wmacros = .false.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item The (very loose) convention in naming input variables is
!latex       \begin{itemize}
!latex       \item Variables that begin with \verb+L+ serve as logical flags. These are usually integers, as integers allow more choices that simply \verb+T+ or \verb+F+.
!latex       \item Variables that begin with \verb+N+ or \verb+M+ serve as integer resolution flags.
!latex       \end{itemize}

!latex \item In the following, all default settings are shown.

!latex \end{enumerate} \subsubsection{physicslist} \begin{enumerate}

!latex \item The namelist \verb+physicslist+ controls the geometry, profiles, and numerical resolution.
!latex       \begin{itemize}
!latex       \verb+namelist/physicslist/+

  namelist/physicslist/&
 Igeometry   ,& !latex \item \verb+Igeometry=3 : integer+ : selects Cartesian, cylindrical or toroidal geometry;
                !latex \begin{itemize}
                !latex \item \verb+Igeometry=1+ : Cartesian; geometry determined by $R$;
                !latex \item \verb+Igeometry=2+ : cylindrical - standard; geometry determined by $R$;
                !latex \item \verb+Igeometry=3+ : toroidal; geometry determined by $R$ {\em and} $Z$;
                !latex \end{itemize}
 pknot       ,& !latex \item \verb+pknot=1 : integer+ knottedness;
 Istellsym   ,& !latex \item \verb+Istellsym = 1 : integer+ : stellarator symmetry is enforced if \verb+Istellsym.eq.1+;
 Lfreebound  ,& !latex \item \verb+Lfreebound = 0 : integer+ : compute vacuum field surrounding plasma; 
 phiedge     ,& !latex \item \verb+phiedge = 1.0 : real+ : total enclosed toroidal magnetic flux;
 curtor      ,& !latex \item \verb+curtor = 0.0 : real+ : total enclosed (toroidal) plasma current;
 curpol      ,& !latex \item \verb+curpol = 0.0 : real+ : total enclosed (poloidal) linking current;
 extcur      ,& !latex \item \verb+extcur : real(1:maxgroups)+ : coil currents;
 gamma       ,& !latex \item \verb+gamma = 0.0 : real+ : adiabatic index; cannot set $|\gamma|=1$;
 Nfp         ,& !latex \item \verb+Nfp = 1 : integer+ : field periodicity;
                !latex \begin{itemize}
                !latex \item all Fourier representations are of the form $\cos(m\t-nN\z)$, where $N\equiv$\verb+Nfp+;
                !latex \item constraint : \verb+Nfp.ge.1+;
                !latex \end{itemize}
 Nvol        ,& !latex \item \verb+Nvol = 1 : integer+ : number of volumes;
                !latex \begin{itemize}
                !latex \item each volume ${\cal V}_l$ is bounded by the ${\cal I}_{l-1}$ and ${\cal I}_{l}$ interfaces;
                !latex \item note that in cylindrical or toroidal geometry, ${\cal I}_{0}$ is the degenerate coordinate axis;
                !latex \item constraint : \verb+Nvol.le.MNvol+;
                !latex \end{itemize}
 Mpol        ,& !latex \item \verb+Mpol = 1 : integer+ : poloidal resolution;
 Ntor        ,& !latex \item \verb+Ntor = 0 : integer+ : toroidal resolution;
                !latex \begin{itemize}
                !latex \item all Fourier representations of doubly-periodic functions are of the form 
                !latex \be f(\t,\z)&=&\sum_{n=0}^{Ntor} f_{0,n}\cos(-nN\z) + \sum_{m=1}^{Mpol}\sum_{n=-Ntor}^{Ntor} f_{m,n}\cos(m\t-nN\z),\\
                !latex     g(\t,\z)&=&\sum_{n=1}^{Ntor} g_{0,n}\sin(-nN\z) + \sum_{m=1}^{Mpol}\sum_{n=-Ntor}^{Ntor} g_{m,n}\sin(m\t-nN\z),\ee
                !latex where $N\equiv$\verb+Nfp+; and internally these ``double'' summations are written as a ``single'' summation,
                !latex e.g. $f=\sum_j f_j \cos(m_j\t-n_j\z)$.
                !latex \end{itemize}
 Lrad        ,& !latex \item \verb+Lrad=2 : integer(MNvol)+ : Chebyshev resolution in each volume;
                !latex \begin{itemize}
                !latex \item constraint : \verb+Lrad(1:Nvol).ge.2+;
                !latex \end{itemize}
 Lconstraint ,& !latex \item \verb+Lconstraint=2 : integer+ : selects constraints; only relevant if \verb+LBeltrami.eq.4, 5, 6, or 7+;
                !latex \begin{itemize}
                !latex \item if \verb+Lconstraint.eq.0+, then $\mu$ and $\psi_p$ are {\em not} varied;
                !latex in this case \verb+helicity+ and the transform (as described by either \verb+pl+, \verb+ql+, \verb+pr+ and \verb+qr+, or \verb+iota+)
                !latex are irrelevant;
                !latex \item if \verb+Lconstraint.eq.1+, then $\mu$ and $\psi_p$ are adjusted
                !latex in order to satisfy the inner and outer interface transform constraints;
                !latex (except in the simple torus, where the enclosed poloidal flux is irrelevant,
                !latex and only $\mu$ is varied to satisfy the outer interface transform constraint);
                !latex in this case \verb+helicity+ is irrelevant;
                !latex \item if \verb+Lconstraint.eq.2+, then $\mu$ is adjusted in order to satisfy the helicity constraint, and $\psi_p$ is held constant;
                !latex in this case the  transform (as described by either \verb+pl+, \verb+ql+, \verb+pr+ and \verb+qr+, or \verb+iota+) is irrelevant;
                !latex \item an exeption: if SQP is used to construct the Beltrami fields, then \verb+Lconstraint+ is irrelevant
                !latex and $\mu$ will be adjusted to constrain the helicity;
                !latex \end{itemize}
 tflux       ,& !latex \item \verb+tflux : real(1:MNvol)+ : toroidal flux, $\psi_t$, enclosed by each interface;
                !latex \begin{itemize}
                !latex \item this is a constraint: \verb+tflux+ is \underline{not} varied
                !latex {\em except} that \verb+tflux+ will be normalized so that \verb+tflux(Nvol)=1.0+,
                !latex so that \verb+tflux+ is arbitrary up to a scale factor;
                !latex \item see also \verb+phiedge+;
                !latex \end{itemize}
 pflux       ,& !latex \item \verb+pflux : real(1:MNvol)+ : poloidal flux, $\psi_p$, enclosed by each interface;
                !latex \begin{itemize}
                !latex \item if \verb+Lconstraint.eq.0+, then \verb+pflux+ is not varied;
                !latex \item if \verb+Lconstraint.eq.1+, then \verb+pflux+ is adjusted to satisfy interface rotational-transform constraint
                !latex (except in the innermost volume if the geometry is cylindrical or toroidal);
                !latex on entry \verb+pflux+ serves as an initial guess and is changed on output;
                !latex \item if \verb+Lconstraint.eq.2+, then \verb+pflux+ is not varied;
                !latex \end{itemize}
 helicity    ,& !latex \item \verb+helicity : real(1:MNvol)+ : helicity, ${\cal K}$, in each volume, ${\cal V}_i$;
                !latex \begin{itemize}
                !latex \item on exit, \verb+helicity+ is set to the computed values of ${\cal K} \equiv \int {\bf A}\cdot{\bf B}\;dv$;
                !latex \end{itemize}
 pscale      ,& !latex \item \verb+pscale=0.0 : real+ : pressure scale factor;
                !latex \begin{itemize}
                !latex \item the initial pressure profile is given by \verb+pscale * press+;
                !latex \end{itemize}
 pressure    ,& !latex \item \verb+pressure : real(1:MNvol)+ : pressure in each volume;
                !latex \begin{itemize}
                !latex \item the pressure is {\em not} held constant, but $p_l V_l^\gamma = P_l$ {\em is} held constant,
                !latex where $P_l$ is determined by the initial pressures and the initial volumes, $V_l$;
                !latex \item (Note that if \verb+gamma=0.0+, then $p_l \equiv P_l$.)
                !latex \item on output, the pressure is given by $p_l= P_l/V_l^\gamma$, where $V_l$ is the final volume;
                !latex \item \verb+pressure+ is only used in calculation of interface force-balance;
                !latex \end{itemize}
 Ladiabatic  ,& !latex \item \verb+Ladiabatic=0 : integer+ : logical flag;
                !latex \begin{itemize}
                !latex \item if \verb+Ladiabatic=0+, the adiabatic constants are determined by the initial pressure and volume;
                !latex \item if \verb+Ladiabatic=1+, the adiabatic constants are determined by the given input \verb+adiabatic+;
                !latex \end{itemize}
 adiabatic   ,& !latex \item \verb+adiabatic : real(1:MNvol)+ : adiabatic constants in each volume;
                !latex \begin{itemize}
                !latex \item the pressure is {\em not} held constant, but $p_l V_l^\gamma = P_l \equiv$\verb+adiabatic+ is constant,
                !latex \item note that if \verb+gamma=0.0+, then \verb+pressure=adiabatic+;
                !latex \item \verb+pressure+ is only used in calculation of interface force-balance;
                !latex \end{itemize}
 mu          ,& !latex \item \verb+mu : real(1:MNvol)+ : helicity-multiplier, $\mu$, in each volume;
                !latex \begin{itemize}
                !latex \item if \verb+Lconstraint.eq.0+, then \verb+mu+ is not varied;
                !latex \item if \verb+Lconstraint.eq.1+, then \verb+mu+ is adjusted to satisfy interface rotational-transform constraint;
                !latex in this case \verb+mu+ serves as an initial guess;
                !latex \item if \verb+Lconstraint.eq.2+, then \verb+mu+ is adjusted to satisfy helicity constraint;
                !latex in this case \verb+mu+ serves as an initial guess;
                !latex \end{itemize}
 pl          ,& !latex \item \verb+pl=0 : integer(0:MNvol)+ :
 ql          ,& !latex \item \verb+ql=0 : integer(0:MNvol)+ :
 pr          ,& !latex \item \verb+pr=0 : integer(0:MNvol)+ :
 qr          ,& !latex \item \verb+qr=0 : integer(0:MNvol)+ :
                !latex \begin{itemize}
                !latex \item "inside" interface rotational-transform is $\iotabar=(p_l+\gamma p_r)/(q_l+\gamma q_r)$,
                !latex        where $\gamma$ is the golden mean, $\gamma=(1 + \sqrt 5 ) / 2 $;
                !latex \item if both $q_l=0$ {\em and} $q_r=0$, then the (inside) interface rotational-transform is defined by \verb+iota+;
                !latex \item if \verb+Lconstraint.ne.1+ then \verb+pl+, \verb+ql+, \verb+pr+ and \verb+qr+, and \verb+iota+, are irrelevant;
                !latex \end{itemize}
 iota        ,& !latex \item \verb+iota : real(0:MNvol)+ : rotational-transform, $\iotabar$, on inner side of each interface;
                !latex \begin{itemize}
                !latex \item only relevant if illogical input for \verb+ql+ and \verb+qr+ are provided;
                !latex \end{itemize}
 lp          ,& !latex \item \verb+lp=0 : integer(0:MNvol)+ :
 lq          ,& !latex \item \verb+lq=0 : integer(0:MNvol)+ :
 rp          ,& !latex \item \verb+rp=0 : integer(0:MNvol)+ :
 rq          ,& !latex \item \verb+rq=0 : integer(0:MNvol)+ :
                !latex \begin{itemize}
                !latex \item "outer" interface rotational-transform is $\iotabar=(p_l+\gamma p_r)/(q_l+\gamma q_r)$,
                !latex       where $\gamma$ is the golden mean, $\gamma=(1 + \sqrt 5 ) / 2 $;
                !latex \item if both $q_l=0$ {\em and} $q_r=0$, then the interface rotational-transform is defined by \verb+oita+;
                !latex \item if \verb+Lconstraint.ne.1+ then \verb+pl+, \verb+ql+, \verb+pr+ and \verb+qr+, and \verb+oita+, are irrelevant;
                !latex \end{itemize}
 oita        ,& !latex \item \verb+oita : real(0:MNvol)+ : rotational-transform, $\iotabar$, on outer side of each interface;
                !latex \begin{itemize}
                !latex \item only relevant if illogical input for \verb+ql+ and \verb+qr+ are provided;
                !latex \end{itemize}
 mupftol     ,& !latex \item \verb+mupftol=1.0e-06 : real+ : accuracy to which $\mu$ and $\psi_p$ are required;
                !latex \begin{itemize}
                !latex \item only relevant if \verb+Lconstraint.gt.0+;
                !latex \end{itemize}
 mupfits     ,& !latex \item \verb+mupfits=8 : integer+ : an upper limit on the transform/helicity constraint iterations;
                !latex \begin{itemize}
                !latex \item constraint: \verb+mupfits > 0+;
                !latex \end{itemize}
 Rac         ,& !latex \item \verb+Rac : real(     0:MNtor        )+ : Fourier harmonics of axis    ;     stellarator symmetric;
 Zas         ,& !latex \item \verb+Zas : real(     0:MNtor        )+ : Fourier harmonics of axis    ;     stellarator symmetric;
 Ras         ,& !latex \item \verb+Ras : real(     0:MNtor        )+ : Fourier harmonics of axis    ; non-stellarator symmetric;
 Zac         ,& !latex \item \verb+Zac : real(     0:MNtor        )+ : Fourier harmonics of axis    ; non-stellarator symmetric;
 Rbc         ,& !latex \item \verb+Rbc : real(-MNtor:MNtor,0:MMpol)+ : Fourier harmonics of boundary;     stellarator symmetric;
 Zbs         ,& !latex \item \verb+Zbs : real(-MNtor:MNtor,0:MMpol)+ : Fourier harmonics of boundary;     stellarator symmetric;
 Rbs         ,& !latex \item \verb+Rbs : real(-MNtor:MNtor,0:MMpol)+ : Fourier harmonics of boundary; non-stellarator symmetric;
 Zbc         ,& !latex \item \verb+Zbc : real(-MNtor:MNtor,0:MMpol)+ : Fourier harmonics of boundary; non-stellarator symmetric;
 Rwc         ,& !latex \item \verb+Rwc : real(-MNtor:MNtor,0:MMpol)+ : Fourier harmonics of wall    ;     stellarator symmetric;
 Zws         ,& !latex \item \verb+Zws : real(-MNtor:MNtor,0:MMpol)+ : Fourier harmonics of wall    ;     stellarator symmetric;
 Rws         ,& !latex \item \verb+Rws : real(-MNtor:MNtor,0:MMpol)+ : Fourier harmonics of wall    ; non-stellarator symmetric;
 Zwc         ,& !latex \item \verb+Zwc : real(-MNtor:MNtor,0:MMpol)+ : Fourier harmonics of wall    ; non-stellarator symmetric;
                !latex \begin{itemize}
                !latex \item The boundary, ${\bf x}=R(\t,\z) \hat {\bf R} + Z(\t,\z) \hat {\bf Z}$, is described by
                !latex \be R & = & \sum_{n=0}^{Ntor}R_{n,0}\cos(-n\z)+\sum_{m=0}^{Mpol}\sum_{n=-Ntor}^{Ntor}R_{n,m}\cos(m\t-nN\z), \\
                !latex     Z & = & \sum_{n=0}^{Ntor}Z_{n,0}\sin(-n\z)+\sum_{m=0}^{Mpol}\sum_{n=-Ntor}^{Ntor}Z_{n,m}\sin(m\t-nN\z), \ee
                !latex where $N\equiv$\verb+Nfp+ and $R_{n,m}\equiv$\verb+Rbc(n,m)+, $Z_{n,m}\equiv$\verb+Zbs(n,m)+.
                !latex \end{itemize}
 Bns         ,& !latex \item \verb+Bns : real(-MNtor:MNtor,0:MMpol)+ : Fourier harmonics of (total)  normal field at boundary;     stellarator symmetric;
 Bnc            !latex \item \verb+Bnc : real(-MNtor:MNtor,0:MMpol)+ : Fourier harmonics of (total)  normal field at boundary; non-stellarator symmetric;

!latex \end{itemize}

!latex \item Comments:
!latex \begin{enumerate}
!latex \item The helicity multipler, $\mu$, is related to the parallel current.
!latex \item The choice \verb+Lconstraint=1+ will suitably adjust $\mu$ to enforce the rotational-transform constraint.
!latex       This is analogous to setting \verb+NCURR=0+ in \verb+VMEC+.
!latex \item The choice \verb+Lconstraint=0+ will not adjust $\mu$. 
!latex       This is analogous to setting \verb+NCURR=1+ in \verb+VMEC+.
!latex \end{enumerate}

!latex \end{enumerate} \subsubsection{numericlist} \begin{enumerate}
!latex \item The namelist \verb+numericlist+ controls internal resolution parameters that the user rarely needs to consider.
!latex \begin{itemize}
!latex \verb+namelist/numericlist/+

  namelist/numericlist/&
 Linitialize ,& !latex \item \verb+Linitialize=0 : integer+ : to initialize geometry using a regularization / extrapolation method;
                !latex \begin{itemize}
                !latex \item if \verb+Linitialize=0+, the geometry of the interior surfaces is provided after the namelists in the input file;
                !latex \item if \verb+Linitialize=1+, the interior surfaces will be intialized as $R_{l,m,n}=R_{N,m,n} \psi_{t,l}^{m/2}$,
                !latex where $R_{N,m,n}$ is the boundary
                !latex       and $\psi_{t,l}$ is the given toroidal flux enclosed by the $l$-th interface, normalized to the total enclosed toroidal flux;
                !latex       a similar extrapolation is used for $Z_{l,m,n}$;
                !latex \item note that the Fourier harmonics of the boundary is {\em always} given by the \verb+Rbc+ and \verb+Zbs+ 
                !latex given in \verb+physicslist+;
                !latex \item if \verb+Linitialize=1+, it is not required to provide the geometry of the interfaces after the namelists;
                !latex \end{itemize}
 Lwall       ,& !latex \item \verb+Lwall=0 : integer+ wall is given as a set of data points; smooth approximation will be constructed;
 phiwall     ,& !latex \item \verb+phiwall=0.9 : real+ value of vacuum potential at wall;
 Ndiscrete   ,& !latex \item \verb+Ndiscrete=2 : integer+ :
                !latex \begin{itemize}
                !latex \item resolution of the real space grid on which fast Fourier transforms are performed is given by \verb+Ndiscrete*Mpol*4+;
                !latex \item constraint \verb+Ndiscrete>0+;
                !latex \end{itemize}
 Nquad       ,& !latex \item \verb+Nquad=-1 : integer+ : the resolution of the Gaussian quadrature, i.e. the radial sub-sub-grid;
                !latex \begin{itemize}
                !latex \item if \verb+Nquad.le.0+, then \verb!Iquad(vvol)=2*Lrad(vvol)-Nquad!;
                !latex \item In the innermost volume, we should use if \verb!Nquad=Mpol/2+2*Nofe+3-Nquad!;
                !latex \end{itemize}
                !latex 
 iMpol       ,& !latex \item \verb+iMpol=-4 : integer+ : Fourier resolution of straight-field-line angle on interfaces;
                !latex \begin{itemize}
                !latex \item the rotational-transform on the interfaces is determined by a transformation to the straight-field-line angle,
                !latex with poloidal resolution given by \verb+iMpol+;
                !latex \item if \verb+iMpol.le.0+, then \verb+iMpol = Mpol - iMpol+;
                !latex \end{itemize}
 iNtor       ,& !latex \item \verb+iNtor=-4 : integer+ : Fourier resolution of straight-field-line angle on interfaces;
                !latex \begin{itemize}
                !latex \item the rotational-transform on the interfaces is determined by a transformation to the straight-field-line angle,
                !latex with toroidal resolution given by \verb+iNtor+;
                !latex \item if \verb+iNtor.le.0+, then \verb+iNtor = Ntor - iNtor+;
                !latex \item if \verb+Ntor.eq.0+, then the toroidal resolution of the angle transformation is set \verb+lNtor=0+.
                !latex \end{itemize}
 Lsparse     ,& !latex \item \verb+Lsparse= 0 : integer+ : controls method used to solve for rotational-transform on interfaces;
                !latex \begin{itemize}
                !latex \item if \verb+Lsparse=0+, the transformation to the straight-field-line angle is computed in Fourier space
                !latex using a dense matrix solver, \verb+F04AAF+;
                !latex \item if \verb+Lsparse=1+, the transformation to the straight-field-line angle is computed in real space
                !latex using a dense matrix solver, \verb+F04ATF+;
                !latex \item if \verb+Lsparse=2+, the transformation to the straight-field-line angle is computed in real space
                !latex using a sparse matrix solver, \verb+F11DEF+;
                !latex \item if \verb+Lsparse=3+, the different methods for constructing the straight-field-line angle are compared;
                !latex \end{itemize}
 Lsvdiota    ,& !latex \item \verb+Lsvdiota= 0 : integer+ : controls method used to solve for rotational-transform on interfaces;
                !latex only relevant if \verb+Lsparse=0+;
                !latex \begin{itemize}
                !latex \item if \verb+Lsvdiota=0+, use standard linear solver to construct straight fieldline angle transformation;
                !latex \item if \verb+Lsvdiota=1+, use SVD method to compute rotational-transform;
                !latex \end{itemize}
 imethod     ,& !latex \item \verb+Imethod= 3 : integer+ : controls iterative solution to sparse matrix
                !latex arising in real-space transformation to the straight-field-line angle;
                !latex only relevant if \verb+Lsparse.eq.2+; see \verb+tr00ab+ for details;
                !latex \begin{itemize}
                !latex \item if \verb+imethod=1+, the method is \verb+RGMRES+;
                !latex \item if \verb+imethod=2+, the method is \verb+CGS+;
                !latex \item if \verb+imethod=3+, the method is \verb+BICGSTAB+;
                !latex \end{itemize}
 iorder      ,& !latex \item \verb+iorder = 2 : integer+ : controls real-space grid resolution for constructing the straight-field-line angle;
                !latex only relevant if \verb+Lsparse>0+;
                !latex determines order of finite-difference approximation to the derivatives;
                !latex \begin{itemize}
                !latex \item if \verb+iorder=2+, 
                !latex \item if \verb+iorder=4+, 
                !latex \item if \verb+iorder=6+, 
                !latex \end{itemize}
 iprecon     ,& !latex \item \verb+Iprecon= 0 : integer+ : controls iterative solution to sparse matrix arising in real-space transformation
                !latex to the straight-field-line angle;
                !latex only relevant if \verb+Lsparse.eq.2+; see \verb+tr00ab+ for details;
                !latex \begin{itemize}
                !latex \item if \verb+iprecon=0+, the preconditioner is `N';
                !latex \item if \verb+iprecon=1+, the preconditioner is `J'; 
                !latex \item if \verb+iprecon=2+, the preconditioner is `S';
                !latex \end{itemize}
 iotatol     ,& !latex \item \verb+iotatol=-1.0 : real+ : tolerance required for iterative construction of straight-field-line angle;
                !latex only relevant if \verb+Lsparse.ge.2+;
 Iswmin      ,& !latex \item \verb+Iswmin=0 : integer+ : construct spectrally condensed angle as a pre-processing step;
                !latex \begin{itemize}
                !latex \item if \verb+Iswmin=0+, the Fourier representation of all interfaces {\em is not} converted to spectrally condensed representation;
                !latex \item if \verb+Iswmin=1+, the Fourier representation of all interfaces, including computational boundary,
                !latex {\em is} converted to spectrally condensed form,
                !latex       but the {\em original boundary} is written to the restart file, \verb+ext.end+, on completion;
                !latex \item if \verb+Iswmin=2+, the Fourier representation of all interfaces is converted to spectrally condensed representation,
                !latex       {\em and} the spectrally condensed representation of the outer boundary is written to the restart file,
                !latex \verb+ext.end+, on completion;
                !latex \end{itemize}
 Lperturb    ,& !latex \item \verb+Lperturb=0 : integer+ : add random perturbation to input geometry;
                !latex \begin{itemize}
                !latex \item if \verb+Lperturb.eq.1+ a random deformation is added to all Fourier harmonics of all internal interfaces;
                !latex \item used for debugging, exploring domain of convergence;
                !latex \end{itemize}
 dperturb    ,& !latex \item \verb+dperturb=1.0e-04 : real+ : magnitude of random perturbation added to each Fourier harmonic;
                !latex \begin{itemize}
                !latex \item only relevant if \verb+Lperturb.eq.1+;
                !latex \end{itemize}
 Lextrap        !latex \item \verb+Lextrap=0 : integer+ : geometry of innermost interface is defined by extrapolation;

!latex \end{itemize}

!latex \item Comments:
!latex \begin{enumerate}
!latex \item The metric elements are constructed at enhanced resolution, \verb+lMpol=2*Mpol+; so \verb+Ndiscrete+ should be set carefully;
!latex \end{enumerate}

!latex \end{enumerate} \subsubsection{locallist} \begin{enumerate}
!latex \item The namelist \verb+locallist+ controls the construction of the Beltrami fields in each volume.
!latex \begin{itemize}
!latex \verb+namelist/locallist/+

  namelist/locallist/&
 LBeltrami,&    !latex \item\verb+LBeltrami=4  integer+
                !latex \begin{itemize}
                !latex \item if \verb+LBeltrami=1,3,5 or 7+, (SQP) then the Beltrami field in each volume is constructed
                !latex by minimizing the magnetic energy with the constraint of fixed helicity; 
                !latex this is achieved by using sequential quadratic programming as provided by the NAG routine \verb+E04UFF+; 
                !latex this approach has the benefit (in theory) of robustly constructing minimum energy solutions
                !latex when multiple, i.e. bifurcated, solutions exist.
                !latex \item if \verb+LBeltrami=2,3,6 or 7+, (Newton) then the Beltrami fields are constructed by employing a standard Newton method
                !latex for locating an extremum of
                !latex $F\equiv \int B^2 dv - \mu (\int {\bf A}\cdot{\bf B}dv-{\cal K})$,
                !latex where $\mu$ is treated as an independent degree of freedom similar to the parameters describing the vector potential
                !latex and ${\cal K}$ is the required value of the helicity; 
                !latex this is the standard Lagrange multipler approach for locating the constrained minimum; 
                !latex this method cannot distinguish saddle-type extrema from minima, and which solution that will be obtained depends on the initial guess;
                !latex \item if \verb+LBeltrami=4,5,6 or 7+, (linear) it is assumed that the Beltrami fields are parameterized by $\mu$;
                !latex in this case, it is only required to solve $\nabla \times {\bf B}=\mu {\bf B}$ which reduces to a system of linear equations;
                !latex $\mu$ may or may not be adjusted iteratively, depending on \verb+Lconstraint+,
                !latex to satisfy either rotational-transform or helicity constraints;
                !latex \item for flexibility and comparison, each of the above methods can be employed; for example:
                !latex \begin{itemize}
                !latex \item if \verb+LBeltrami=1+, only the SQP    method will be employed;
                !latex \item if \verb+LBeltrami=2+, only the Newton method will be employed;
                !latex \item if \verb+LBeltrami=4+, only the linear method will be employed; 
                !latex \item if \verb+LBeltrami=3+, the SQP and the Newton method are used;
                !latex \item if \verb+LBeltrami=5+, the SQP and the linear method are used;
                !latex \item if \verb+LBeltrami=6+, the Newton and the linear method are used;
                !latex \item if \verb+LBeltrami=7+, all three methods will be employed;
                !latex \end{itemize}
                !latex \end{itemize}
 Linitgues,&    !latex \item\verb+Linitgues=1  integer+ controls how initial guess for Beltrami field is constructed;
                !latex \begin{itemize}
                !latex \item only relevant for routines that require an initial guess for the Beltrami fields, such as the SQP and Newton methods,
                !latex or the sparse linear solver;
                !latex \item if \verb+Linitgues=0+, the initial guess for the Beltrami field is trivial;
                !latex \item if \verb+Linitgues=1+, the initial guess for the Beltrami field is an integrable approximation;
                !latex \item if \verb+Linitgues=2+, the initial guess for the Beltrami field is read from file; 
                !latex \end{itemize}
 Lposdef  ,&    !latex \item\verb+Lposdef=T : logical+ : indicates whether the Beltrami linear system is positive definite;
                !latex \begin{itemize}
                !latex \item if \verb+Lposdef=T+, the Beltrami linear system is assumed to be positive definite;
                !latex \item if \verb+Lposdef=F+, the Beltrami linear system is assumed not   positive definite;
                !latex \item only relevant for the linear method of constructing the Beltrami fields;
                !latex \end{itemize}
 Nmaxexp        !latex \item\verb+Nmaxexp=32 : integer+ : indicates maximum exponent used to precondition Beltrami linear system near coordinate singularity;
                !latex \begin{itemize}
                !latex \item a factor of $s^{\bar m_j/2}$ is extracted, where $\bar m \equiv $\verb+min(Nmaxexp,m_j)+;
                !latex \end{itemize}
!Lprecon ,&     !latex \item\verb+Lprecon=1 : integer+ : choice of preconditioner used in sparse linear solver, \verb+F11JEF+; REDUNDANT; 13 Sep 13;
                !latex \begin{itemize}
                !latex \item if \verb+Lprecon=0+, no preconditioner is used;
                !latex \item if \verb+Lprecon=1+, the Jacobi preconditioner is used;
                !latex \item if \verb+Lprecon=2+, successive over relaxation preconditioner is used;
                !latex \item only used if \verb+Ldense.eq.0+ (or \verb+Ldense.eq.2+);
                !latex \end{itemize}
!sparseits,&    !latex \item\verb+sparseits=-1 : integer+ : REDUNDANT; 13 Sep 13;
                !latex \begin{itemize}
                !latex \item if \verb+sparseits.le.0+, maximum allowed iterations defaults to \verb+abs(sparseits)*+$N$, where $N$ is the size of the matrix;
                !latex \item only used if \verb+Ldense.eq.0+ (or \verb+Ldense.eq.2+);
                !latex \end{itemize}
!ssoromega,&    !latex \item\verb+ssoromega=1.0 : real+ : successive over-relaxation parameter, $\omega$; REDUNDANT; 13 Sep 13;
                !latex \begin{itemize}
                !latex \item constraint $0.0 \le \omega \le 2.0$;
                !latex \item only required if \verb+Ldense.eq.0+ (or \verb+Ldense.eq.2+) and \verb+Lprecon.eq.2+;
                !latex \end{itemize}
!sparsetol      !latex \item\verb+sparsetol=-1.0 : real+ : accuracy requested in sparse solution to Beltrami linear system; REDUNDANT; 13 Sep 13;
                !latex \begin{itemize}
                !latex \item if \verb+sparsetol.le.zero+, required accuracy defaults to machine precision;
                !latex \item only used if \verb+Ldense.eq.0+ (or \verb+Ldense.eq.2+);
                !latex \end{itemize}
!Liotasolv      !latex \item\verb+Liotasolv=1 : integer+ : choice of solving for straight-field line coordinates on the interfaces; REDUNDANT; 13 Sep 13;
                !latex \begin{itemize}
                !latex \item if \verb+Liotasolv=1,3,5,7+, the linear transformation is solved using sparse LU decomposition;
                !latex \item if \verb+Liotasolv=2,3,6,7+, the linear transformation is solved using standard, iterative linear solver;
                !latex \item if \verb+Liotasolv=4,5,6,7+, the linear transformation is solved using singular value decomposition (SVD) methods;
                !latex \end{itemize}

!latex \end{itemize}

!latex \item Comments:
!latex \begin{enumerate}
!latex \item The transformation to straight-field-line coordinates is singular when the rotational-transform of the interfaces is rational;
!latex       however, the rotational-transform is still well defined.
!latex       In this case, the SVD method (\verb+Liotasolv=4,5,6,7+) must be used.
!latex \end{enumerate}

!latex \end{enumerate} \subsubsection{globallist} \begin{enumerate}
!latex \item The namelist \verb+globallist+ controls the search for global-force balance, i.e. the search for the interface geometry so that $[[p+B^2/2]]=0$.
!latex \begin{itemize}
!latex \verb+ namelist/globallist/+

  namelist/globallist/& 
 Le04dgf     ,& !latex \item \verb+Le04dgf=F : logical+ : redundant; 10 Oct 12;
 Lminimize   ,& !latex \item \verb+Lminimize=0 : integer+ : use minimization methods to construct equilibrium solutions;
                !latex \begin{itemize}
                !latex \item if \verb+Lminimize=1+, then use the preconditioned conjugate gradient method, provided by \verb+E04DGF+;
                !latex \item if \verb+Lminimize=2+, then use the homespun non-linear conjugate gradient method, described in \verb+pc01aa+;
                !latex \item if \verb+Lminimize=4+, then use the modified Newton method, employing first and second derivatives,
                !latex provided by \verb+E04LYF+; under construction;
                !latex \end{itemize}
 Lfindzero   ,& !latex \item \verb+Lfindzero=0 : integer+ : use Newton methods to find zero of force balance;
                !latex \begin{itemize}
                !latex \item if \verb+Lfindzero=1+, then use NAG routine \verb+C05NDF+, which uses function values only;
                !latex \item if \verb+Lfindzero=2+, then use NAG routine \verb+C05PDF+, which uses derivative information;
                !latex \end{itemize}
!pwidth      ,& !latex \item \verb+pwidth=4 : integer+ : ! REDUNDANT; 13 Sep 13;
!qwidth      ,& !latex \item \verb+qwidth=4 : integer+ : spectral condensation parameter; REDUNDANT; 13 Sep 13;
 Lcondense   ,& !latex \item \verb+Lcondense=   1 : integer+ : selection of spectral condensation; REDUNDANT; 18 Jul 14;
 escale      ,& !latex \item \verb+escale = 0.0 : real+ : exponential weight factor used in force-imbalance harmonics;
 pcondense   ,& !latex \item \verb+pcondense= 2.0 : real+ : spectral condensation parameter;
 qcondense   ,& !latex \item \verb+qcondense=-1.0 : real+ : spectral condensation parameter; redundant; ! 04 Dec 14;
 wcondense   ,& !latex \item \verb+wcondense= 1.0 : real+ : redundant;
 wpoloidal   ,& !latex \item \verb+wpoloidal= 1.0 : real+ : minimum poloidal length constraint radial exponential factor;
                !latex \begin{itemize}
                !latex \item the angle freedom is exploited to minimize $\sum_{m,n}(m^p+n^q)(R_{m,n}^2+Z_{m,n}^2)$ with respect to tangential variations;
                !latex \end{itemize}
 forcetol    ,& !latex \item \verb+forcetol = 1.0e-10 : real+ : required tolerance in force balance error; only used as an initial check;
                !latex \begin{itemize}
                !latex \item if the initially supplied interfaces are consistent with global pressure balance to within \verb+forcetol+,
                !latex       then the geometry of the interfaces is not altered;
                !latex \item if not, then the geometry of the interfaces is changed in order to bring the configuration into forcebalance
                !latex       so that the geometry of interfaces is within \verb+c05xtol+, defined below, of the true solution
                !latex \item to force execution of either \verb+C05NDF+ and \verb+C05PDF+ regardless of the initial force imbalance, set \verb+forcetol+ $< 0$;
                !latex \end{itemize}
 normalerr   ,& !latex \item \verb+normalerr = 1.0e-06 : real+ : required tolerance in free-boundary iterations;
 norblend    ,& !latex \item \verb+norblend  = 0.0     : real+ : normal blend;
 maxfbits    ,& !latex \item \verb+maxfbits = 20 : integer+ : maximum allowed free-boundary iterations;
 ForceErr    ,& !latex \item \verb+ForceErr       : real+ : on exit, force balance error; redundant; 02 Nov 12;
 Energy      ,& !latex \item \verb+Energy   = 0.0 : real+ : on exit, total energy; redundant; 02 Nov 12;
 Le04lyf     ,& !latex \item \verb+Le04lyf=F : logical+ : redundant; 10 Oct 12;
 c05xtol     ,& !latex \item \verb+c05xtol=1.0e-12 : real+ : required tolerance in position, {\bf x};
                !latex \begin{itemize}
                !latex \item used by both \verb+C05NDF+ and \verb+C05PDF+; see the NAG documents for further details on how the error is defined;
                !latex \item constraint \verb+c05xtol.gt.0.0+;
                !latex \end{itemize}
 factor      ,& !latex \item \verb+factor : real+ : redundant; this has been replaced by \verb+c05factor+; 10 Oct 12;
 c05factor   ,& !latex \item \verb+c05factor=1.0e-02 : real+ : used to control initial step size; provided to NAG routines \verb+C05NDF+, \verb+C05PDF+
                !latex \begin{itemize}
                !latex \item constraint \verb+c05factor.gt.0.0+;
                !latex \item only relevant if \verb+Lfindzero.gt.0+;
                !latex \item used by both \verb+C05NDF+ and \verb+C05PDF+;
                !latex \end{itemize}
 LreadGF     ,& !latex \item \verb+LreadGF=T : logical+ : read $\nabla {\bf F}$ from file \verb+.GF+;
 verify      ,& !latex \item \verb+verify=-1 : integer+ : instruct \verb+E04DGF+ to ``verify'' user supplied gradients are correct;
                !latex only relevant if \verb+Le04gdf=T+;
                !latex \begin{itemize}
                !latex \item if \verb+verify=-1+, no checks on user supplied gradients are requested;
                !latex \item if \verb+verify= 0+, simple checks on user supplied gradients are requested;
                !latex \item if \verb+verify= 1+, extensive checks on user supplied gradients are requested;
                !latex \end{itemize}
 maxstep     ,& !latex \item \verb+maxstep=1.0e-03 : real+ : provided to \verb+E04DGF+; only relevant if \verb+Le04gdf=T+;
 opsilon     ,& !latex \item \verb+opsilon=1.0e-00 : real+ : weighting of pressure imbalance;
 epsilon     ,& !latex \item \verb+epsilon=1.0e-00 : real+ : weighting of spectral-width constraints;
 upsilon     ,& !latex \item \verb+upsilon=1.0e-00 : real+ : weighting of poloidal-length constraint;
 apsilon     ,& !latex \item \verb+apsilon=1.0e-00 : real+ : weighting of minimal length spectral constraint; REDUNDANT; 18 Jul 14;
 maxiter        !latex \item \verb+maxiter=-1 : integer+ : maximum iterations allowed in \verb++E04DGF;
!latex \end{itemize}

!latex \item Comments:
!latex \begin{enumerate}
!latex \item An equilibrium is obtained when ${\bf F}({\bf x})=0$.
!latex \item The ``force'' vector, ${\bf F}$, is a combination of pressure-imbalance Fourier harmonics, i.e. $[[p+B^2/2]]_{m,n}$
!latex and spectral-condensation constraints.
!latex \item The vector ${\bf x}$ represents the internal interface geometry harmonics.
!latex \end{enumerate}

!latex \end{enumerate} \subsubsection{diagnosticslist} \begin{enumerate}
!latex \item The namelist \verb+diagnosticslist+ controls post-processor diagnostics, such as \Poincare plot resolution, $\dots$,
!latex \begin{itemize}
!latex \verb+namelist/diagnosticslist/+

  namelist/diagnosticslist/&
 odetol     ,&  !latex \item \verb+odetol=1.0e-07 : real+ : o.d.e. integration tolerance for all field line tracing routines;
 absreq     ,&  !latex \item \verb+absreq=1.0e-06 : real+ : absolute accuracy on virtual casing integral for external field;
 relreq     ,&  !latex \item \verb+relreq=1.0e-06 : real+ : relative accuracy on virtual casing integral for external field;
 absacc     ,&  !latex \item \verb+absacc=1.0e-04 : real+ : absolute accuracy required for integration of toroidal flux;
 epsr       ,&  !latex \item \verb+epsr  =1.0e-06 : real+ : absolute accuracy required for integration of Beltrami error;
!divertorR  ,&  !latex \item \verb+divertorR=-1.0 : real+ : redundant; 10 Oct 12;
!divertorZ  ,&  !latex \item \verb+divertorZ= 1.0 : real+ : redundant; 10 Oct 12;
 nPpts      ,&  !latex \item \verb+nPpts=0 : integer+ : number of toroidal transits used (per trajectory) in following field lines
                !latex for constructing \Poincare plots;
                !latex if \verb+nPpts<1+, no \Poincare plot is constructed;
 nPtrj      ,&  !latex \item \verb+nPtrj=-1 : integer(1:MNvol+1)+ : number of trajectories in each annulus following in constructing \Poincare plot;
                !latex \begin{itemize}
                !latex \item if \verb+nPtrj(l)<0+, then \verb+nPtrj(l)=Ni(l)+,
                !latex where \verb+Ni(l)+ is the grid resolution used to construct the Beltrami field in volume $l$;
                !latex \end{itemize}
 Mpqits     ,&  !latex \item \verb+Mpqits=10 : integer+ : maximum allowed iterations in finding periodic orbits; see \verb+pq00aa,pq02aa+;
 Lpqsym     ,&  !latex \item \verb+Lpqsym=1 : integer+ : indicates whether periodic orbits lie on symmetry line; see \verb+pq00aa,pq02aa+;
                !latex \begin{itemize}
                !latex \item if \verb+Lpqsym=1+, the periodic orbits are assumed to lie on the `symmetry line', $\theta=0$.
                !latex \item constraint; \verb+Lpqsym=0+ or \verb+Lpqsym=1+;
                !latex \end{itemize}
 p1         ,&  !latex \item \verb+p1=0 : integer(1:MNvol)+ : selection of periodic orbits and noble irrational surface; see \verb+pq01aa+;
 q1         ,&  !latex \item \verb+q1=0 : integer(1:MNvol)+ : selection of periodic orbits and noble irrational surface; see \verb+pq01aa+;
 p2         ,&  !latex \item \verb+p2=0 : integer(1:MNvol)+ : selection of periodic orbits and noble irrational surface; see \verb+pq01aa+;
 q2         ,&  !latex \item \verb+q2=0 : integer(1:MNvol)+ : selection of periodic orbits and noble irrational surface; see \verb+pq01aa+;
 pqs        ,&  !latex \item \verb+pqs= -2.0 : real(1:MNvol)+ : initial radial guess for location of periodic orbit;
 pqt        ,&  !latex \item \verb+pqt=  0.0 : real(1:MNvol)+ : initial angle  guess for location of periodic orbit;
!pqR        ,&  !latex \item \verb+pqR=  0.0 : real(1:MNvol)+ : output cylindrical; redundant; 10 Oct 12;
!pqZ        ,&  !latex \item \verb+pqZ=  0.0 : real(1:MNvol)+ : output cylindrical; redundant; 10 Oct 12;
!Npl        ,&  !latex \item \verb+Npl= 1 : integer+ piecewise linear resolution for Lagrangian integration; ! TO BE DELETED; 11 Aug 14;
!Nunstable  ,&  !latex \item \verb+Nunstable=1000    : integer+ : unstable manifold; resolution of line segment emanating from unstable periodic orbit;
!dunstable  ,&  !latex \item \verb+dunstable=1.0e-05 : real   + : unstable manifold; length     of line segment emanating from unstable periodic orbit;
!Munstable  ,&  !latex \item \verb+Nunstable=20      : integer+ : unstable manifold; iterations of \Poincare map;
 npq        ,&  !latex \item \verb+npq=0 : integer(1:MNvol)+ : number of required additional periodic orbits as defined by Fibonacci series
                !latex beginning from \verb+p1/q1+ and \verb+p2/q2+;
!Mirrits    ,&  !latex \item \verb+Mirrits=0 : integer+ : maximum allowed iterations in finding irrational surface; if \verb+Mirrits<1+,
                !latex irrational surfaces are not located;
!irrMpol    ,&  !latex \item \verb+irrMpol=50 : integer+ : poloidal resolution of irrational surface; only required if \verb+Mirrits>0+;
                !latex see \verb+pq01aa,ir00aa+;
!irrNtor    ,&  !latex \item \verb+irrNtor=25 : integer+ : toroidal resolution of irrational surface; only required if \verb+Mirrits>0+;
                !latex see \verb+pq01aa,ir00aa+;
!irrsvdcut  ,&  !latex \item \verb+irrsvdcut=1.0e-12 : real+ : SVD inversion singular cutoff used to locate irrational surfaces;
                !latex only required if \verb+Mirrits>0+;
                !latex see \verb+ir00aa+;
!irrtol     ,&  !latex \item \verb+irrtol=1.0e-08 : real+ : irrational surface tolerance; only required if \verb+Mirrits>0+;see \verb+ir00aa+;
 Lwrpj      ,&  !latex \item \verb+Lwrpj=F : logical+ : the pressure-jump Hamiltonian input files are written; under re-construction;
!Nghd       ,&  !latex \item \verb+Nghd=0 : integer+ : poloidal resolution of ghost-surfaces $\equiv$ quadratic-flux minimizing surfaces;
 LHessian   ,&  !latex \item \verb+LHessian=F : logical+ : to compare $\nabla {\bf F}$ constructed in parallel to series construction; redundant;
 LHevalues  ,&  !latex \item \verb+LHevalues=F : logical+ : to compute eigenvalues of $\nabla {\bf F}$;
 LHevectors ,&  !latex \item \verb+LHevectors=F : logical+ : to compute eigenvectors (and also eigenvalues) of $\nabla {\bf F}$;
 Lperturbed ,&  !latex \item \verb+Lperturbed=0 : integer+ : to compute linear, perturbed equilibrium;
 dpp        ,&  !latex \item \verb+dpp=1+ : integer+ : perturbed harmonic;
 dqq        ,&  !latex \item \verb+dqq=1+ : integer+ : perturbed harmonic;
 Lcurlerr   ,&  !latex \item \verb+Lcurlerr=F : logical+ : redundant;
 Lcheck     ,&  !latex \item \verb+Lcheck=0 : integer+ : implement various checks;
                !latex \begin{itemize}
                !latex \item if \verb+Lcheck=-1+, a LaTeX document is produced on runtime;
                !latex \item if \verb+Lcheck=0+, no additional check on the calculation is performed;
                !latex \item if \verb+Lcheck=1+, the error in the current, i.e. $\nabla\times{\bf B}-\mu{\bf B}$ is computed as a post-diagnostic;
                !latex \item if \verb+Lcheck=2+, the analytic calculation of the derivatives of the interface transform w.r.t.
                !latex       the helicity multiplier, $\mu$, and the enclosed poloidal flux, $\Delta\psi_p$, are compared to a finite-difference estimate;
                !latex       only if \verb+Lconstraint.eq.1+;
                !latex       only for \verb+dspec+ executable, i.e. must compile with \verb+DFLAGS="-D DEBUG"+;
                !latex \item if \verb+Lcheck=3+, the analytic calculation of the derivatives of the volume w.r.t. interface Fourier harmonic
                !latex       is compared to a finite-difference estimate;
                !latex       must set \verb+Lfindzero+$=2$, set \verb+forcetol+ sufficiently small and set \verb+LreadGF=F+,
                !latex       so that the matrix of second derivatives is calculated,
                !latex       only for \verb+dspec+ executable, i.e. must compile with \verb+DFLAGS="-D DEBUG"+;
                !latex \item if \verb+Lcheck=4+, the analytic calculation of the derivatives of the magnetic field, $B^2$, at the interfaces
                !latex       is compared to a finite-difference estimate;
                !latex       must set \verb+Lfindzero+$=2$, and ensure that \verb+forcetol+ is sufficiently small;
                !latex       only for \verb+dspec+ executable, i.e. must compile with \verb+DFLAGS="-D DEBUG"+;
                !latex \item if \verb+Lcheck=5+, the analytic calculation of the matrix of the derivatives of the force imbalance
                !latex       is compared to a finite-difference estimate;
                !latex \end{itemize}
 Ltiming        !latex \item \verb+Ltiming=T : logical+ : to check timing;

!latex \end{itemize}

!latex Comments:
!latex \begin{enumerate}
!latex \item No comments yet.
!latex \end{enumerate}

!latex \end{enumerate} \subsubsection{screenlist} \begin{enumerate}
!latex \item The namelist \verb+screenlist+ controls screen output.
!latex \begin{itemize}
!latex \verb+namelist/screenlist/+

  namelist/screenlist/&
! NSCREENLIST ! namelist screenlist; this is expanded by Makefile; do not remove;
 Wreadin , &  !latex \item Every subroutine, e.g. \verb+xy00aa.h+, has its own write flag, \verb+Wxy00aa+.
 Wwritin , &
 Wmacros

!latex \end{itemize}

!latex \end{enumerate} \subsection{Input geometry} \begin{enumerate}
!latex \item The geometry of the $l$-th interface, for $l=0,N$ where $N\equiv$ \verb+Nvol+, is described by a set of Fourier harmonics,
!latex using an arbitrary poloidal angle,
!latex \be R_l(\t,\z)&=&\sum_{j}R_{j,l}\cos(m_j\t-n_j\z), \\ Z_l(\t,\z)&=&\sum_{j}Z_{j,l}\sin(m_j\t-n_j\z). \ee
!latex \item These harmonics are read from the \verb+ext.spec+ file and come directly after the namelists described above.
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
!latex \item The coordinate axis corresponds to $j=0$ and the outermost boundary corresponds to $j=$\verb+Nvol+.
!latex \item An arbitrary selection of harmonics may be inluded in any order, but only those within the range specified by \verb+Mpol+ and \verb+Ntor+
!latex will be used.
!latex \item The geometry of {\em all} the interfaces, i.e. $l=0,N$, including the degenerate `coordinate-axis' interface, must be given.

! note that all variables in namelist need to be broadcasted in readin;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end module inputlist

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!module ezmgrid ! only required for free-boundary; interpolation of mgrid file; 28 Nov 12;
!
!  use ezspline_obj
!  use ezspline
!
!  implicit none
!
!  INTEGER            :: mgridnfp
!  REAL               :: pi2mgridnfp
!
!  type(ezspline3_r8) :: sTBr, sTBz, sTBp ! spline objects for total field in mgrid file;
!
!end module ezmgrid

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

module allglobal

  use constants
  use typedefns

  implicit none

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER              :: myid, ncpu       ! mpi variables;
  REAL                 :: cpus             ! initial time;

  REAL                 :: pi2nfp           !       pi2/nfp     ; assigned in readin;
  REAL                 :: pi2pi2nfp
  REAL                 :: pi2pi2nfphalf
  REAL                 :: pi2pi2nfpquart

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER              :: Mvol

  LOGICAL              :: YESstellsym, NOTstellsym ! internal shorthand copies of Istellsym, which is an integer input; 16 Jan 13;

  REAL   , allocatable :: gchebyshev(:,:,:,:) ! local workspace;

  REAL   , allocatable :: lchebyshev(:,:) ! local workspace;
  
  REAL   , allocatable :: TTll(:,:,:) ! derivatives of Chebyshev polynomials at the inner and outer interfaces;

  LOGICAL, allocatable :: ImagneticOK(:)   ! used to indicate if Beltrami fields have been correctly constructed;
  
  REAL   , allocatable :: beltramierror(:,:)  ! to store the integral of |curlB-mu*B| computed by jo00aa; 11 Apr 16;
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \end{enumerate}
!latex \subsection{Internal variables}
!latex \subsubsection{Fourier representation}
!latex \begin{enumerate} 

!latex \item The Fourier description of stellarator-symmetric (even) periodic functions is 
!latex \be f(\t,\z) = \sum_{n=0}^{N} f_{0,n} \cos(-n\z) + \sum_{m=1}^{M}\sum_{n=-N}^{N} f_{m,n} \cos(m\t-n\z),
!latex \ee
!latex where the resolution is given on input, $M\equiv$\verb+ Mpol+ and $N\equiv$\verb+ Ntor+.
!latex \item For convenience, the Fourier summations are written as
!latex \be f(\s,\t,\z) &=& \sum_j f_j(s) \cos( m_j \t - n_j \z ),\\
!latex     f(\s,\t,\z) &=& \sum_j f_j(s) \sin( m_j \t - n_j \z ),
!latex \ee
!latex for $j=1,$ \verb+mn+, where 
!latex \verb+mn+$ = N + 1 +  M  ( 2 N + 1 )$.
!latex \item The arrays \verb+im+ and \verb+in+ contain the $m_j$ and $n_j$.

  INTEGER              :: mn  ! total number of Fourier harmonics for coordinates/fields; calculated from Mpol,Ntor in readin;
  INTEGER, allocatable :: im(:), in(:) ! Fourier modes; set in readin;

  REAL,    allocatable :: halfmm(:)
! REAL,    allocatable :: regularmm(:)

  REAL,    allocatable :: psifactor(:,:)

  REAL,    allocatable :: expmmnn(:) ! exponential weight on force-imbalance harmonics; used in fc02aa; 04 Dec 14;
  
  REAL,    allocatable :: mpnq(:) ! spectral condensation factors; 18 Jul 14;
 
! INTEGER, allocatable :: dnjn(:,:)

!latex \item Enhanced resolution is required for the metric elements, $g_{ij}/\sqrt g$, which is given by \verb+mne+, \verb+ime+., and \verb+ine+.
!latex The Fourier resolution here is determined by \verb+lMpol=2*Mpol+ and \verb+lNtor=2*Ntor+.
  INTEGER              :: mne ! enhanced resolution for metric elements;
  INTEGER, allocatable :: ime(:), ine(:)

!latex \item Enhanced resolution is required for the transformation to straight-field line angle on the interfaces,
!latex which is given by \verb+mns+, \verb+ims+., and \verb+ins+.
!latex The Fourier resolution here is determined by \verb+iMpol+ and \verb+iNtor+.

  INTEGER              :: mns ! enhanced resolution for straight field line transformation;
  INTEGER, allocatable :: ims(:), ins(:)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  REAL                 :: xoffset = 1.0 ! used to normalize NAG routines;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \end{enumerate} \subsubsection{Interface geometry} \begin{enumerate}

!latex \item The Fourier harmonics of the interfaces are contained in \verb+iRbc(1:mn,0:Mvol)+ and \verb+iZbs(1:mn,0:Mvol)+, where
!latex \verb+iRbc(l,j)+, \verb+iZbs(l,j)+ contains the Fourier harmonics, $R_j$, $Z_j$, of the $l$-th interface.

  REAL,    allocatable :: iRbc(:,:) , iZbs(:,:)   ! interface surface geometry;     stellarator symmetric;
  REAL,    allocatable :: iRbs(:,:) , iZbc(:,:)   ! interface surface geometry; non-stellarator symmetric;

  REAL,    allocatable :: dRbc(:,:) , dZbs(:,:)   ! interface surface geometry;     stellarator symmetric; linear deformation;
  REAL,    allocatable :: dRbs(:,:) , dZbc(:,:)   ! interface surface geometry; non-stellarator symmetric;

  REAL,    allocatable :: iRij(:,:) , iZij(:,:)   ! interface surface geometry; real space;
  REAL,    allocatable :: dRij(:,:) , dZij(:,:)   ! interface surface geometry; real space;
  REAL,    allocatable :: tRij(:,:) , tZij(:,:)   ! interface surface geometry; real space;

  REAL,    allocatable :: iBns(:)                 ! non-zero normal field;
  REAL,    allocatable :: iBnc(:)
  
  REAL,    allocatable :: iCns(:)                 ! coil   normal field at computational boundary;
  REAL,    allocatable :: iCnc(:)
  
  REAL,    allocatable :: iPns(:)                 ! plasma normal field at computational boundary;
  REAL,    allocatable :: iPnc(:)
  
  REAL                 :: Bnserror                ! used to terminate free-boundary iterations; computed in bn00aa;
! REAL                 :: looperr                 ! used to measure free-boundary consistency; 21 Oct 12;

  REAL,    allocatable :: lRbc(:)   , lZbs(:)     ! local workspace;
  REAL,    allocatable :: lRbs(:)   , lZbc(:)     ! local workspace;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \end{enumerate} \subsubsection{Fourier Transforms} \begin{enumerate}

!latex \item The coordinate geometry and fields are mapped to/from Fourier space and real space using the NAG routine \verb+C06FUF+.
!latex \item The resolution of the real space grid is given by \verb+Nt=Ndiscrete*4*Mpol+ and \verb+Nz=Ndiscrete*4*Ntor+.
!latex \item Trigonometric information required for the fast Fourier transform's is saved in \verb+trigm(1:2*Nt)+, \verb+trign(1:2*Nz)+,
!latex       and \verb+trigwk(1:2*Ntz)+, where \verb+Ntz=Nt*Nz+.

  INTEGER              :: Nt, Nz, Ntz, hNt, hNz ! discrete resolution; Ntz=Nt*Nz shorthand;
  REAL                 :: soNtz ! one / sqrt (one*Ntz); shorthand;

  CHARACTER            :: isr ! required for C06FUF;
  REAL   , allocatable :: trigm(:), trign(:), trigwk(:) ! these are set in readin and contain trigonometric factors; cannot be changed;

!latex \item Various workspace arrays are allocated. 
!latex These include \verb+Rij(1:Ntz,0:3,0:3)+ and \verb+Zij(1:Ntz,0:3,0:3)+, which contain the coordinates in real space and their derivatives;
!latex \verb+sg(0:3,Ntz)+, which contains the Jacobian and its derivatives;
!latex and \verb+guv(0:6,0:3,1:Ntz)+, which contains the metric elements and their derivatives.

  REAL   , allocatable :: Rij(:,:,:), Zij(:,:,:), Xij(:,:,:), Yij(:,:,:), sg(:,:), guvij(:,:,:,:), gvuij(:,:,:,:) ! real-space; 11 Feb 13;

 !INTEGER              :: mnNtor
 !REAL   , allocatable :: Rj(:,:), Zj(:,:) ! 18 Jul 14;

  REAL   , allocatable :: guvmne(:,:), guvmno(:,:) ! Fourier space; 11 Feb 13;
  
  INTEGER, allocatable :: guvmnks(:,:,:), guvmnka(:,:,:) ! identification of Fourier modes; 16 Jan 13;

  INTEGER, allocatable :: guvmnk(:,:)

  INTEGER, allocatable :: iotakkii(:), iotaksub(:,:), iotakadd(:,:), iotaksgn(:,:) ! identification of Fourier modes; 29 Jan 13;

  REAL,    allocatable :: efmn(:), ofmn(:), cfmn(:), sfmn(:) ! Fourier harmonics; dummy workspace;
  REAL,    allocatable :: evmn(:), odmn(:), comn(:), simn(:) ! Fourier harmonics; dummy workspace;

  REAL,    allocatable :: ijreal(:), ijimag(:), jireal(:), jiimag(:) ! real grid; dummy workspace;

  REAL,    allocatable :: Bsupumn(:,:,:), Bsupvmn(:,:,:) ! 11 Oct 12; tangential field on interfaces; required for virtual casing construction of field;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \end{enumerate} \subsubsection{Metric information} \begin{enumerate}

!latex \item The metric information is 

  REAL,    allocatable :: TTee(:,:,:,:,:)
  REAL,    allocatable :: TTeo(:,:,:,:,:)
  REAL,    allocatable :: TToe(:,:,:,:,:)
  REAL,    allocatable :: TToo(:,:,:,:,:)

  REAL,    allocatable :: Te(:,:,:)
  REAL,    allocatable :: To(:,:,:)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \end{enumerate} \subsubsection{Enclosed fluxes} \begin{enumerate}

!latex \item The input variables are the toroidal and poloidal fluxes, given in \verb+tflux(1:Nvol)+ and \verb+pflux(1:Nvol)+.
!latex \item However, the Beltrami fields in each volume depend only on the enclosed fluxes, \verb+dtflux(1:Nvol)+ and \verb+dpflux(1:Nvol)+,
!latex which are the independent degrees of freedom required to satisfy the rotational transform constraints, see \verb+ma02aa+.

  REAL,    allocatable :: dtflux(:), dpflux(:) ! \delta \psi_{toroidal} and \delta \psi_{poloidal} in each annulus;

  REAL,    allocatable :: sweight(:) ! minimum poloidal length constraint weight; 04 Dec 14;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \end{enumerate} \subsubsection{The vector potential and the Beltrami linear system} \begin{enumerate}
  
!latex \item The covariant components of the vector potential are written as
!latex       \be            A_\t & = & \sum_i \sum_{l=0}^L A_{\t,e,i,l} \; T_{l}(s) \cos\a_i + \sum_i \sum_{l=0}^L A_{\t,o,i,l} \; T_{l}(s) \sin\a_i \\
!latex                      A_\z & = & \sum_i \sum_{l=0}^L A_{\z,e,i,l} \; T_{l}(s) \cos\a_i + \sum_i \sum_{l=0}^L A_{\z,o,i,l} \; T_{l}(s) \sin\a_i ,
!latex       \ee
!latex       where $T_l(s)$ are the Chebyshev polynomials and $\a_i \equiv m_i \t - n_i \z$.

!latex \item The following internal arrays are declared in \verb+al00aa+
!latex
!latex       \verb+dAte(0,i)%s(l)+$\equiv A_{\t,e,i,l}$
!latex
!latex       \verb+dAze(0,i)%s(l)+$\equiv A_{\z,e,i,l}$
!latex
!latex       \verb+dAto(0,i)%s(l)+$\equiv A_{\t,o,i,l}$
!latex
!latex       \verb+dAzo(0,i)%s(l)+$\equiv A_{\z,o,i,l}$

  type(subgrid), allocatable :: Ate(:,:,:), Aze(:,:,:)
  type(subgrid), allocatable :: Ato(:,:,:), Azo(:,:,:)

  INTEGER      , allocatable :: Fso(:,:), Fse(:,:)
 
! REAL                       :: IG(1:2) ! toroidal, plasma current (Itor) and "poloidal", linking current (Gpol); renamed to curtor and curpol;

!latex \item In each volume, the total degrees of freedom in the Beltrami linear system is \verb+Nmagneticdof(1:Nvol)+.
!latex This depends on \verb+Mpol+, \verb+Ntor+ and \verb+Lrad(vvol)+.
   
  INTEGER, allocatable :: Nmagneticdof(:) ! degrees of freedom in Beltrami fields in each annulus;

  LOGICAL              :: Lcoordinatesingularity
  LOGICAL              :: Lplasmaregion, Lvacuumregion
  
!latex \item The energy, $W \equiv \int \! dv {\; \bf B}\cdot{\bf B}$, and helicity, $K\equiv \int \! dv \; {\bf A}\cdot{\bf B}$, functionals may be written
!latex      \be W & = & \frac{1}{2} \; a_i \; A_{i,j} \; a_j + a_i \; B_{i,j} \; \psi_j + \frac{1}{2} \; \psi_i \; C_{i,j} \; \psi_j \label{eq:energymatrix} \\
!latex          K & = & \frac{1}{2} \; a_i \; D_{i,j} \; a_j + a_i \; E_{i,j} \; \psi_j + \frac{1}{2} \; \psi_i \; F_{i,j} \; \psi_j \label{eq:helicitymatrix}
!latex      \ee
!latex       where ${\bf a} \equiv \{ A_{\t,e,i,l}, A_{\z,e,i,l}, A_{\t,o,i,l}, A_{\z,o,i,l}, f_{e,i}, f_{o,i} \}$ contains the independent degrees of freedom
!latex       and $\boldpsi \equiv \{\Delta \psi_t,\Delta \psi_p\}$.

   REAL,    allocatable :: dMA(:,:), dMB(:,:), dMC(:,:) ! energy and helicity matrices; quadratic forms; 17 Jan 13;
   REAL,    allocatable :: dMD(:,:), dME(:,:), dMF(:,:) ! energy and helicity matrices; quadratic forms; 17 Jan 13;

   REAL,    allocatable :: solution(:,:)

   REAL,    allocatable :: MBpsi(:), MEpsi(:) ! matrix vector products; 26 Feb 13;
   REAL                 :: psiMCpsi, psiMFpsi

   REAL,    allocatable :: BeltramiInverse(:,:)

!latex \item The arrays, \verb+dAij(0,i,j)+$\equiv A_{i,j}$, etc. will be allocated and deallocated in each volume as required.
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \end{enumerate} \subsubsection{Transformation to straight-field-line angle} \begin{enumerate}

!latex \item Given the Beltrami fields in any volume, the rotational-transform on the adjacent interfaces 
!latex may be determined by constructing the straight-field line angle on the interfaces.
!latex \item This is performed in \verb+tr00ab+.
!latex \item The rotational-transform on the inner and outer interfaces of each volume is stored in \verb+diota(1:Mvol,-1:2,0:1)+. 
!latex       The first argument labels the inner or outer interface,
!latex           the second labels derivative, with -1 indicating the derivative wrt geometry,
!latex                                               0 the rotational-transform itself,
!latex                                           and 1 and 2 being the derivatives wrt $\mu$ and $\Delta \psi_p$;
!latex       and the third argument labels volume.
!latex This is assigned values in \verb+mp00aa+ after calling \verb+tr00ab+.

  REAL   , allocatable :: diota(:,:,:) ! measured rotational transform on inner/outer interfaces for each volume;

  REAL   , allocatable :: glambda(:,:,:,:) ! save initial guesses for iterative calculation of rotational-transform; 21 Apr 13;

  INTEGER              :: lmns
! REAL   , allocatable :: dlambda(:,:,:,:) ! transformation to straight fieldline angle; constructed in tr00ab; 06 Mar 15; NOWHERE USED ; 15 Sep 15;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \end{enumerate} \subsubsection{Construction of `force'} \begin{enumerate}

!latex \item The force vector is comprised of \verb+Bomn+ and \verb+Iomn+.

  REAL,    allocatable ::  Bemn(:,:,:),  Iomn(:,:), Somn(:,:,:), Pomn(:,:,:)
  REAL,    allocatable ::  Bomn(:,:,:),  Iemn(:,:), Semn(:,:,:), Pemn(:,:,:)
  REAL,    allocatable ::  BBe(:), IIo(:), BBo(:), IIe(:)

  REAL,    allocatable ::  Bsubtemn(:,:,:), Bsubzemn(:,:,:) ! covariant components of the tangential field on interfaces; 01 Jul 14;
  REAL,    allocatable ::  Bsubtomn(:,:,:), Bsubzomn(:,:,:) ! covariant components of the tangential field on interfaces; 01 Jul 14;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \end{enumerate} \subsubsection{Parallel construction of derivative matrix} \begin{enumerate}

!latex \item The derivatives of force balance, $[[p+B^2/2]]$, and the spectral constraints (see \verb+sw03aa+), with respect to the interface geometry
!latex is constructed in parallel by \verb+fc02aa+.

!latex \item Force balance across the $l$-th interface depends on the fields in the adjacent interfaces.

  INTEGER              :: lgeometricaldof !       geometrical degrees of freedom associated with each interface;                   ; 19 Apr 13;
  INTEGER              :: Ngeometricaldof ! total geometrical degrees of freedom                               ; assigned in al00aa; 19 Apr 13;

  REAL,    allocatable :: dBBdRZ(:,:,:)
  REAL,    allocatable :: dIIdRZ(:  ,:)

  REAL,    allocatable :: dFFdRZ(:,:,:,:,:) ! derivatives of B^2 at the interfaces wrt geometry     ; 20 Jun 14;
  REAL,    allocatable :: dBBdmp(:,:,:,:  ) ! derivatives of B^2 at the interfaces wrt mu and dpflux; 20 Jun 14;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \end{enumerate} \subsubsection{internal variable :: dmupfdx} \begin{enumerate}

!latex \item The internal variable \verb+dmupfdx(1:Mvol,1:2,1:lgeometricaldof,0:1)+ is allocated/deallocated in \verb+jk03aa+, and \verb+he01aa+ if selected.
!latex \item The information in \verb+dmupfdx+ describes how the helicity multiplier, $\mu$, and the enclosed poloidal flux, $\Delta \psi_p$, 
!latex       must vary as the geometry is varied in order to satisfy the interface transform constraint.
!latex \item The rotational transform on the inner or outer interface of a given volume depends on the magnetic field in that volume, i.e.
!latex       \be \iotabar_\pm = \iotabar({\bf B}_\pm),
!latex       \ee
!latex       so that
!latex       \be \delta \iotabar_\pm = \frac{\partial \iotabar_\pm}{\partial {\bf B}_\pm} \cdot \delta {\bf B_\pm}.
!latex       \ee
!latex       \item The magnetic field depends on the inner and outer interface geometry (represented here as $x_j$), the helicity multiplier, 
!latex       and the enclosed poloidal flux, i.e. \be {\bf B_\pm} = {\bf B_\pm}(x_j, \mu, \Delta \psi_p),
!latex       \ee
!latex       so that 
!latex       \be \delta {\bf B_\pm} = \frac{\partial {\bf B}_\pm}{\partial x_j          } \delta x_j
!latex                              + \frac{\partial {\bf B}_\pm}{\partial \mu          } \delta \mu
!latex                              + \frac{\partial {\bf B}_\pm}{\partial \Delta \psi_p} \delta \Delta \psi_p.
!latex       \ee
!latex       \item The constraint to be enforced is that $\mu$ and $\Delta \psi_p$ will vary as the geometry is varied 
!latex             so that the rotational-transform constraint on the inner and outer interface is enforced,
!latex       \i.e. 
!latex       \be \left(\begin{array}{ccc} \ds \frac{\partial \iotabar_-}{\partial {\bf B}_-} \cdot \frac{\partial {\bf B}_-}{\partial \mu          } & , & 
!latex                                    \ds \frac{\partial \iotabar_-}{\partial {\bf B}_-} \cdot \frac{\partial {\bf B}_-}{\partial \Delta \psi_p} \\ 
!latex                                    \ds \frac{\partial \iotabar_+}{\partial {\bf B}_+} \cdot \frac{\partial {\bf B}_+}{\partial \mu          } & , & 
!latex                                    \ds \frac{\partial \iotabar_+}{\partial {\bf B}_+} \cdot \frac{\partial {\bf B}_+}{\partial \Delta \psi_p}
!latex                   \end{array} \right)
!latex           \left(\begin{array}{c} \ds \frac{\partial \mu}{\partial x_j} \\ \ds \frac{\partial \Delta \psi_p}{\partial x_j} \end{array} \right) = 
!latex         - \left(\begin{array}{c} \ds \frac{\partial \iotabar_-}{\partial {\bf B}_-} \cdot \frac{\partial {\bf B}_-}{\partial x_j} \\
!latex                                  \ds \frac{\partial \iotabar_+}{\partial {\bf B}_+} \cdot \frac{\partial {\bf B}_+}{\partial x_j} \end{array} \right)
!latex       \ee
!latex \item This solution is constructed in \verb+fc02aa+.
!latex \item A finite-difference estimate is computed if \verb+Lcheck.eq.4+.
!latex \item This information is used to adjust the calculation of how force balance, i.e. $B^2$ at the interfaces, 
!latex       varies with geometry at fixed interface rotational transform. Given 
!latex       \be B_\pm^2 = B_\pm^2 (x_j, \mu, \Delta \psi_p),
!latex       \ee
!latex       we may derive
!latex       \be \frac{\partial B_\pm^2}{\partial x_j} = \frac{\partial B_\pm^2}{\partial x_j          }                     
!latex                                                 + \frac{\partial B_\pm^2}{\partial \mu          } \frac{\partial \mu          }{\partial x_j}
!latex                                                 + \frac{\partial B_\pm^2}{\partial \Delta \psi_p} \frac{\partial \Delta \psi_p}{\partial x_j}
!latex       \ee

  REAL,    allocatable :: dmupfdx(:,:,:,:)  ! derivatives of mu and dpflux wrt geometry at constant interface transform; 20 Jun 14;

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

!latex \end{enumerate} \subsubsection{Trigonometric factors} \begin{enumerate}

!latex \item To facilitate construction of the metric integrals, various trigonometric identities are exploited.
!latex \item The required information is saved in
!latex
!latex \verb+ijk(1:mn,1:mn)+

! REAL   , allocatable :: coszeta(:), sinzeta(:), costeta(:), sinteta(:)

  REAL   , allocatable :: cosi(:,:), sini(:,:)

! REAL   , allocatable :: dxglmn(:,:) ! these are only used in me00ab; seems wasteful on memory;

 !INTEGER, allocatable :: ijk(:,:), ikj(:,:), bmask(:)

  INTEGER, allocatable :: bjk(:,:) ! definition of coordinate axis; 11 Aug 14;

  INTEGER, allocatable :: djkp(:,:), djkm(:,:) ! for calculating cylindrical volume; 02 Sep 14;

! LOGICAL, allocatable :: doubleex(:,:)
! INTEGER, allocatable :: doubleie(:,:), doubleip(:,:), doubleim(:,:), doublein(:,:)

! INTEGER, allocatable :: tripleipp(:,:), tripleipm(:,:), tripleimp(:,:), tripleimm(:,:)

!latex \item The following are used for volume integrals (see \verb+vo00aa+)
!latex \be a_{i,j,k} &=& 4 \; m_k \ooint \cos(\alpha_i)\cos(\alpha_j)\cos(\alpha_k) /(2\pi)^2 , \\
!latex     b_{i,j,k} &=& 4 \; m_j \ooint \cos(\alpha_i)\sin(\alpha_j)\sin(\alpha_k) /(2\pi)^2 ,
!latex \ee 

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \end{enumerate} \subsubsection{Volume integrals} \begin{enumerate}

!latex \item The energy functional, $F \equiv \sum_l F_l$, where
!latex \be F_l \equiv \left( \int_{{\cal V}_l} \frac{p_l}{\gamma-1} + \frac{B_l^2}{2} dv \right)
!latex            = \frac{P_l}{\gamma-1}V_l^{1-\gamma}+\int_{{\cal V}_l} \frac{B_l^2}{2} dv, \label{eq:energy}
!latex \ee
!latex where the second expression is derived using $p_l V_l^\gamma=P_l$, where $P_l$ is the adiabatic-constant.
!latex In \Eqn{energy}, it is implicit that ${\bf B}$ satisfies (i) the toroidal and poloidal flux constraints; 
!latex (ii) the interface constraint, ${\bf B}\cdot\nabla s=0$; and (iii) the helicity constraint (or the transform constraint)
!latex \item The derivatives of $F_l$ with respect to the inner and outer adjacent interface geometry are stored in
!latex
!latex \verb&dFF(1:Nvol,0:1,0:mn+mn-1)&, where
!latex
!latex $         F_l                      \equiv$ \verb&dFF(l,0,    0)&
!latex
!latex $\partial F_l / \partial R_{l-1,j} \equiv$ \verb&dFF(ll,0,   j)& 
!latex
!latex $\partial F_l / \partial Z_{l-1,j} \equiv$ \verb&dFF(ll,0,mn+j)& 
!latex
!latex $\partial F_l / \partial R_{l  ,j} \equiv$ \verb&dFF(ll,1,   j)& 
!latex
!latex $\partial F_l / \partial Z_{l  ,j} \equiv$ \verb&dFF(ll,1,mn+j)&
!latex

!latex \item The volume integrals $\int dv$, $\int B^2 \; dv$ and $\int {\bf A}\cdot{\bf B} \; dv$ in each volume
!latex       are computed and saved in \verb+volume(0:2,1:Nvol)+.

  REAL   , allocatable :: lBBintegral(:) ! B.B      integral;
  REAL   , allocatable :: lABintegral(:) ! A.B      integral;

! REAL                 :: dBBintegral    ! B.B      integral; derivative wrt R, Z;

! REAL   , allocatable :: oBBintegral(:) ! B.B      integral; original; used to normalize; perhaps irrelevant;


! REAL                 :: dABintegral    ! A.B      integral; derivative wrt R, Z;

  REAL   , allocatable :: vvolume(:) ! volume integral of \sqrt g; computed in vo00aa;
  REAL                 :: dvolume    ! derivative of volume wrt interface geometry;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
 LOGICAL               :: Lmgridexist = .false.
 LOGICAL               :: Lmgridhasbeensplined = .false.
 INTEGER               :: nextcur = 1            ! default given in readin; re-set in mg00aa; 11 Oct 12; 
 REAL                  :: Rmin, Rmax, Zmin, Zmax ! set in mg00aa;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

 INTEGER               :: nFreeIterations ! number of free-boundary iterations; only relevant if Lvacuum.ge.3; iterations terminated by maxfbits; 18 Oct 12;
 LOGICAL               :: Lcontinueiterations ! a logical flag used for convenience; indicates whether free-boundary iterations should continue; 18 Oct 12;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \end{enumerate} \subsubsection{Diagnostics : periodic orbits} \begin{enumerate}

!latex \item Periodic orbits may be constructed. This information is contained in a \verb+periodicorbit+ structure called 
!latex
!latex \verb+pqorbit(1:Nvol,maxval(npq))+.

  type periodicorbit
     INTEGER              :: pq(2)            ! integers defining periodcity; \t(\z+2pi q) = \t(\z) + 2pi p;
     INTEGER              :: ok               ! flag indicating whether periodic orbit has been successfully located;
     REAL                 :: to, so           ! location of periodic orbit on \z=0; poloidal angle, \t, and radial coordinate, \s;
     REAL                 :: Ro, Zo           ! location of periodic orbit on \z=0; cylindrical   ,  R, and cylindrical      ,  Z;
     REAL                 :: residue          ! residue;
     REAL                 :: error            ! error;
     INTEGER              :: its              ! required iterations
     REAL                 :: wr(2)  , wi(2)   ! eigenvalues  of tangent map at periodic orbit;
     REAL                 :: vr(2,2), vi(2,2) ! eigenvectors of tangent map at periodic orbit;
     REAL,    allocatable :: s(:), t(:)       ! location of periodic orbit on all \z=const. planes;
     REAL,    allocatable :: a(:)             ! 
     REAL,    allocatable :: tt(:), ss(:)
     REAL,    allocatable :: dt(:)
  end type periodicorbit
  
  type(periodicorbit), allocatable :: pqorbit(:,:)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! REAL                 :: Dzeta ! will be used for Lagrangian integration; 14 Nov 12; ! TO BE DELETED; 11 Aug 14;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \end{enumerate} \subsubsection{Diagnostics : quadratic-flux minimizing surfaces} \begin{enumerate}

!latex \item Quadratic-flux minimizing surfaces may be constructed.
!latex All the information is contained in a \verb+QFmin+ structure called \verb+qfms(1:Nvol,mnpq)+.

  type QFmin
     INTEGER              :: i ! flag;
     INTEGER              :: pq(2) ! periodicity;
     REAL,    allocatable :: trvs(:,:,:) ! theta, radial, action-gradient;
  end type QFmin

  type(QFmin),allocatable :: qfms(:,:)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \end{enumerate} \subsubsection{Diagnostics : irrational surfaces} \begin{enumerate}

!latex \item Irrational, invariant flux surfaces may be constructed.
!latex The Fourier harmonics are described by \verb+irrmn+, \verb+irrim+ and \verb+irrin+.

  INTEGER              :: irrmn ! total number of Fourier harmonics for irrational surfaces; calculated from MNirr in readin;
  INTEGER, allocatable :: irrim(:), irrin(:) ! Fourier modes for irrational surfaces; set in readin;

  type irrational
     REAL              :: iota ! prescribed irrational transform
     REAL              :: Itor, Gpol ! secular parts;
     INTEGER           :: its ! iterations required to locate surface; screen diagnostics only I think;
     REAL              :: err ! accuracy;
     INTEGER           :: Npts ! resolution of discrete curve;
     REAL, allocatable :: sm(:), tm(:) ! radial and poloidal angle Fourier harmonics of irrational, invariant curve in toroidal coordinates;
     REAL              :: tflux ! enclosed toroidal flux;
     REAL, allocatable :: si(:), ti(:), di(:) ! invariant curve in discrete; interpolation coefficients;
     REAL, allocatable :: Rbc(:), Zbs(:) ! cylindrical Fourier harmonics of irrational/invariant surface in toroidal coordinates;
     REAL, allocatable :: Rbs(:), Zbc(:) ! cylindrical Fourier harmonics of irrational/invariant surface in toroidal coordinates;
     REAL, allocatable :: fmn(:) ! covariant field derived from surface potential;
  end type irrational

!latex The irrational surface Fourier harmonics, etc. are contained in a \verb+irrational+ structure called \verb+irrsurf+.

  type(irrational) :: irrsurf ! Fourier harmonics etc. of irrational surfaces ! shall construct one irrational surface per annulus;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! internal global variables; internal logical variables; default values are provided here; these may be changed according to input values;

  INTEGER              :: ivol ! labels volume; some subroutines (called by NAG) are fixed argument list but require the volume label;

  INTEGER              :: Ltangent = 1 ! indicates whether tangent field needs to be calculated by bf00aa, vc00aa;

  INTEGER              :: Lfieldlinedirection = 1 ! indicates whether field lines are followed or backwards;

! INTEGER              :: oiMpol, oiNtor ! save original iMpol and iNtor, which may be changed internally; 

  REAL                 :: Bzeta ! toroidal (contravariant) field; calculated in bf00aa; required to convert \dot \t to B^\t, \dot s to B^s;

! INTEGER              :: iqfmin, jgd, kgd ! I don't know if all these are required to be global; see bf00aa_out and bf00aa_end;

  REAL                 :: actiongradient ! used to construct QFmin surfaces; probably not required to be global;

! LOGICAL              :: Iwrpj ! internal copy of Lwrpj;

  INTEGER, allocatable :: Iquad(:) ! internal copy of Nquad; 16 Jan 13;

  REAL   , allocatable :: gaussianweight(:,:), gaussianabscissae(:,:)

  LOGICAL              :: LBeltramiLinear, LBeltramiNewton, LBeltramiSeQuad ! controls selection of Beltrami field solver; depends on LBeltrami;

  REAL                 :: oRZp(1:3) ! used in mg00aa to determine (\s,\t,\z) given (R,Z,p);

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  type derivative
     LOGICAL :: L
     INTEGER :: innout
     INTEGER :: ii
     INTEGER :: irz
     INTEGER :: issym
  end type derivative
  
 type(derivative)  :: DifferentiateGeometry
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! the following are required for the sparse solution to the straight-field-line coordinate on the interfaces;

! LOGICAL              :: Liotasparse, Liotastandard, LiotaSVD

! LOGICAL, allocatable :: AngleLogical(:,:)
! INTEGER              :: AngleNZ, AngleLIRN, AngleLICN, AngleIDISP(1:10)
! INTEGER, allocatable :: AngleIRN(:), AngleICN(:), AngleIKEEP(:), AngleIVECT(:), AngleJVECT(:)
! REAL   , allocatable :: AngleSparseMatrix(:)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER, allocatable :: NNZ(:) ! used to count non-zero elements of Beltrami linear system;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! the following are miscellaneous flags required for the virtual casing field, external (vacuum) field integration,  . . .
  
  REAL                 :: xyz(1:3) ! point at which external field is required; NAG routines employed are fixed argument list, but require position;
  
  REAL                 :: virtualcasingfactor = one / ( four*pi * pi2 )
  
! INTEGER              :: iexvol ! this must be set early; 24 Oct 12;
  
  INTEGER              :: IBerror ! for computing error in magnetic field; 29 Apr 13;

  LOGICAL              :: Llatex

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
contains

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine readin !latex \end{enumerate} \subsection{Subroutine readin} \begin{enumerate}

!latex \item The master node reads the input namelist and sets various secondary variables; variables are then broadcast;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item \verb+readin+ : Reads (and broadcasts) the input file.
!latex Various checks on input consistency are performed.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  use constants

  use numerical

  use fileunits

  use inputlist

  use cputiming

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS

  LOGICAL              :: inputfileexist
  INTEGER              :: vvol, mm, nn, nb, imn, ix, ii, jj, ij, kk, mj, nj, mk, nk, ip, lMpol, lNtor, X02BBF, iargc, iarg, numargs, mi, ni
  REAL                 :: X02AJF, X01AAF, machpi, xx, G05CAF
  REAL,    allocatable :: RZRZ(:,:) ! local array used for reading interface Fourier harmonics from file;
  
  CHARACTER            :: ldate*8, ltime*10, arg*100
  
  external             :: X02BBF
  
  BEGIN(readin)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  cput = GETTIME
  
  machprec = X02AJF() ; vsmall = 100*machprec ; small = 100*vsmall ; sqrtmachprec = sqrt(machprec) ! returns machine precision;
  
  largestint = X02BBF(0.0) ! this is probably nowhere used; 04 Dec 14;
  
  machpi = X01AAF(0.0) ! use NAG routine to compute pi; 04 Dec 14;
  
  FATALMESS( globals, abs(machpi-pi).gt.vsmall,incorrect value for pi)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( myid.eq.0 ) then ! only the master node reads input file and sets secondary variables;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
#ifdef CHECKNAG
   call A00AAF() ! check NAG version;
#endif
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   call date_and_time(ldate,ltime)
   
   write(ounit,'("readin : ", 10x ," : ")')
   write(ounit,1000)cput-cpus, ldate(1:4), ldate(5:6), ldate(7:8), ltime(1:2), ltime(3:4), ltime(5:6), machprec, vsmall, small, largestint
   
1000 format("readin : ",f10.2," : date="a4"/"a2"/"a2" , "a2":"a2":"a2" ; machine precision="es9.2" ; vsmall="es9.2" ; small="es9.2" ; largestint="i12" ;")
   
   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : reading ext from command line ;")')cput-cpus
   endif
   
!latex \end{enumerate}\subsubsection{command line inputs} \begin{enumerate}
   
   call getarg( 1, ext) !latex \item The input file name, \verb!ext!, is given as the first command line input, and the input file itself is \verb!ext.spec!.
   
   if( ext .eq. "" .or. ext.eq. " " .or. ext .eq. "-h" .or. ext .eq. "-help" ) then
    ;write(ounit,'("readin : ", 10x ," : ")')
    ;write(ounit,'("readin : ", 10x ," : file extension must be given as first command line argument ; extra command line options = -help -readin ;")')
    if( ext .eq. "-h" .or. ext .eq. "-help" ) then
     write(ounit,'("readin : ", 10x ," : ")')
     write(ounit,'("readin : ", 10x ," : the input file ext.spec must contain the input namelists; see global.pdf for description ;")')
    endif
    FATALMESS( readin, .true., the input file does not exist) ! if not, abort;
   endif
   
   write(ounit,'("readin : ", 10x ," : ")')
   write(ounit,'("readin : ",f10.2," : ext = ",a100)') cput-cpus, ext
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
!latex \item Additional command line inputs recognized are:
!latex \begin{enumerate}
   
   write(ounit,'("readin : ", 10x ," : ")')
   
   numargs = iargc()
   
   if( numargs.gt.1 ) then
    iarg = 1
    do while ( iarg < numargs )
     iarg = iarg + 1 ; call getarg( iarg, arg)
     select case( arg )
!latex \item \verb+-help, -h+ ; will give help information to user; under construction;
     case("-help","-h") ; write(ounit,'("readin : ",f10.2," : myid=",i3," : command line options = -readin ;")')cput-cpus, myid
!latex \item \verb+-readin+ ; will immediately set \verb+Wreadin=T+; this may be over-ruled when \verb+namelist/screenlist/+ is read;
     case("-readin"   ) ; Wreadin = .true.
     case("-p4pg"     ) ; iarg = iarg + 1 ; call getarg( iarg, arg)
     case("-p4wd"     ) ; iarg = iarg + 1 ; call getarg( iarg, arg)
     case default       ; write(ounit,'("readin : ",f10.2," : myid=",i3," : argument not recognized ; arg = ",a100)') cput-cpus, myid, arg
     end select
    enddo
   endif
   
!latex \end{enumerate}
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   inquire( file=trim(ext)//".spec", exist=inputfileexist ) ! check if file exists;
   FATALMESS( readin, .not.inputfileexist, the input file does not exist) ! if not, abort;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   open( iunit, file=trim(ext)//".spec", status="old")
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
! read and check namelist/physicslist/ ; recall this is inside the if( myid.eq.0 );
   
   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : reading physicslist     from ext.spec ;")')cput-cpus
   endif
   
   read(iunit,physicslist)
   
   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : read    physicslist     from ext.spec ;")')cput-cpus
   endif
   
   if( Lfreebound.eq.0 ) Mvol = Nvol
   if( Lfreebound.gt.0 ) Mvol = Nvol + 1 ! include vacuum region; 04 April 13;
   
   FATALMESS( readin, Igeometry.lt.1 .or. Igeometry.gt.3, invalid geometry)
   FATALMESS( readin, pknot.lt.1, invalid pknot)
   FATALMESS( readin, Nfp.le.0, invalid Nfp)
   FATALMESS( readin, Mpol.lt.0 .or. Mpol.gt.MMpol, invalid poloidal resolution: may need to recompile with higher MMpol)
   FATALMESS( readin, Ntor.lt.0 .or. Ntor.gt.MNtor, invalid toroidal resolution: may need to recompile with higher MNtor)
   FATALMESS( readin, Nvol.lt.1 .or. Nvol.gt.MNvol, invalid Nvol: may need to recompile with higher MNvol)
   FATALMESS( readin, mupftol.le.zero, mupftol is too small)
   FATALMESS( readin, abs(one+gamma).lt.vsmall, 1+gamma appears in denominator in fc02aa)
   FATALMESS( readin, abs(one-gamma).lt.vsmall, 1-gamma appears in denominator in fu00aa)
   
!  FATALMESS( readin, ql(Nvol)*qr(Nvol).eq.0 .and. abs(iota(Nvol)).lt.machprec, definition of transform has changed) ! 29 Apr 15;
   
   FATALMESS( readin, abs(tflux(Nvol)).lt. vsmall, enclosed toroidal flux cannot be zero)
   
   tflux(1:Nvol) = tflux(1:Nvol) / tflux(Nvol) ! normalize toroidal flux;  4 Feb 13;
   
   FATALMESS( readin, tflux(1).lt.zero, toroidal flux is not monotonic)
   do vvol = 2, Nvol
    FATALMESS( readin, tflux(vvol)-tflux(vvol-1).lt.small, toroidal flux is not monotonic)
   enddo
   
! Commented out since boundary input handling was generalized ; 25 Jan 16
!   do nn = -Ntor, Ntor
!    do mm = 0, Mpol
!     ; FATALMESS( readin, mm.eq.0 .and. nn.lt.0 .and. abs(Rbc(nn,mm)).gt.vsmall, illegal choice of mode)
!     if( Igeometry.eq.3 ) then
!      ;FATALMESS( readin, mm.eq.0 .and. nn.lt.0 .and. abs(Zbs(nn,mm)).gt.vsmall, illegal choice of mode)
!     endif
!     if( NOTstellsym ) then
!      ;FATALMESS( readin, mm.eq.0 .and. nn.lt.0 .and. abs(Rbs(nn,mm)).gt.vsmall, illegal choice of mode)
!      if( Igeometry.eq.3 ) then
!       FATALMESS( readin, mm.eq.0 .and. nn.lt.0 .and. abs(Zbc(nn,mm)).gt.vsmall, illegal choice of mode)
!      endif
!     endif
!    enddo
!   enddo
   
!else ! could hardwire knot representation here; 22 Oct 13;
   
   do nextcur = maxgroups, 1, -1
    if( abs(extcur(nextcur)).gt.vsmall ) exit ! this will provide default value for nextcur; this may be over-written in mg00aa;
   enddo
   
   do vvol = 1, Mvol
    FATALMESS( readin, (Lrad(vvol)/2)*2.ne.Lrad(vvol), have assumed that Chebyshev resolution is even            )
    FATALMESS( readin,  Lrad(vvol)     .lt.         2, have assumed that Chebyshev resolution is greater than two)
   enddo
   
   FATALMESS( readin, mupfits.le.0, must give C05PCF a postive integer value for the maximum iterations = mupfits given on input)
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
! read and check namelist/numericlist/ ; recall this is inside the if( myid.eq.0 );
   
   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : reading numericlist     from ext.spec ;")')cput-cpus
   endif
   
   read(iunit,numericlist)!,iostat=ios)
   
   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : read    numericlist     from ext.spec ;")')cput-cpus
   endif
   
   FATALMESS(readin, Ndiscrete.le.0, error)
   FATALMESS(readin, Iswmin.lt.0 .or. Iswmin.gt.3, error)
   
   FATALMESS(readin, Lfreebound.gt.0 .and. Lconstraint.gt.0 .and. Lsparse.eq.0, have not implemented dense Fourier angle transformation in vacuum region)
   
   FATALMESS(readin, iotatol.gt.one, illegal value for sparse tolerance)
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
! read and check namelist/locallist/
   
   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : reading locallist      from ext.spec ;")')cput-cpus
   endif
   
   read(iunit,locallist)!,iostat=ios)
   
   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : read    locallist      from ext.spec ;")')cput-cpus
   endif
   
   FATALMESS(readin,LBeltrami.lt.0 .or. LBeltrami.gt.7, error)
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
! read and check namelist/globallist/
   
   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : reading globallist   from ext.spec ;")')cput-cpus
   endif
   
   read(iunit,globallist)!,iostat=ios)
   
   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : read    globallist   from ext.spec ;")')cput-cpus
   endif
   
   FATALMESS(readin, normalerr.le.zero, error)
   FATALMESS(readin, maxfbits.lt.zero, error)
   
   FATALMESS(readin, c05xtol.le.zero, error)
   FATALMESS(readin, c05factor.le.zero, error)
   
   FATALMESS(readin, Igeometry.eq.3 .and. pcondense.le.zero, pcondense must be positive)
!  FATALMESS(readin, Igeometry.eq.3 .and. qcondense.lt.zero, qcondense must be non-negative) ! 04 Dec 14;
   
   ForceErr = one  ! 02 Nov 12;
   Energy   = zero ! 02 Nov 12;
   
   if( Le04lyf           ) write(ounit,'("readin : ", 10x ," : LE04LYF IS REDUNDANT  ; please remove from namelist/globallist/;")')
   if( factor   .ge. zero) write(ounit,'("readin : ", 10x ," : FACTOR   IS REDUNDANT ; please remove from namelist/globallist/;")')
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
! read and check namelist/diagnosticslist/
   
   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : reading diagnosticslist from ext.spec ;")')cput-cpus
   endif
   
   read(iunit,diagnosticslist)!,iostat=ios)
   
   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : read    diagnosticslist from ext.spec ;")')cput-cpus
   endif
   
   FATALMESS(readin, odetol.le.zero, input error)
   FATALMESS(readin, absreq.le.zero, input error)
   FATALMESS(readin, relreq.le.zero, input error)
   FATALMESS(readin, absacc.le.zero, input error)
   FATALMESS(readin, epsr  .le.zero, input error)
   FATALMESS(readin, nPpts .lt.0   , input error)
   FATALMESS(readin, Mpqits.lt.0   , input error)
   FATALMESS(readin, Lpqsym.lt.0   , input error)
   FATALMESS(readin, Lpqsym.gt.1   , input error)
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   Llatex = .false.
   
   select case( Lcheck )
   case( -1 ) ; Llatex = .true.
   case(  0 ) ; ! no checks            ; 01 Jul 14;
   case(  1 ) ; ! curl B - mu B        ; 01 Jul 14;
   case(  2 ) ; ! transform derivatives; 01 Jul 14;
   case(  3 ) ; ! volume derivatives   ; 01 Jul 14;
   case(  4 ) ; ! field derivatives    ; 01 Jul 14;
   case(  5 ) ; ! hessian              ; 01 Jul 14;
   case default
    FATALMESS(global, .true., invalid Lcheck)
   end select
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
#ifdef DEBUG
   
   if( Llatex ) then
    
    open( wunit, file = trim(ext)//".tex", status = "unknown" ) ! LaTeX formatted output; 14 Jan 15;
    
    ;write(wunit,'("\input{{head}} \code{{runtime diagnostics}} \item{{Diagnostic information created at runtime.}} \end{{enumerate}}")')
    
    ;write(wunit,'("\subsection{{numerical resolution [global]}}")')
    ;write(wunit,'("\begin{{enumerate}}")')
    ;write(wunit,'("\item [] $M="i4"$")') Mpol
    ;write(wunit,'("\item [] $N="i4"$")') Ntor
    do vvol = 1, Nvol
     write(wunit,'("\item [] $L("i4")="i3"$")') vvol, Lrad(vvol)
    enddo
    ;write(wunit,'("\end{{enumerate}}")')
    
   endif
   
#endif
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
! read and check namelist/screenlist
   
!  write(ounit,*) "Wreadin", Wreadin
   
   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : reading screenlist      from ext.spec ;")') cput-cpus
   endif
   
   read(iunit,screenlist)
   
   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : read    screenlist      from ext.spec ;")') cput-cpus
   endif
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
! additional checks on input logic;
   
!  FATALMESS(readin,Lvacuum.ge.3 .and. ( Lminimize.eq.0 .and. Lfindzero.eq.0 ),free-boundary iterations require force-balance) ! 26 Oct 12;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
! screen output of namelist/physicslist/
   
   if( Lfreebound.eq.0 ) Mvol = Nvol
   if( Lfreebound.gt.0 ) Mvol = Nvol + 1 ! include vacuum region; 04 April 13;
   
   write(ounit,'("readin : ", 10x ," : ")')
   
   write(ounit,1010) cput-cpus, Igeometry, pknot, Istellsym
   write(ounit,1011)            Lfreebound, nextcur, phiedge, curtor, curpol
   write(ounit,1012)            gamma
   write(ounit,1013)            Nfp, Nvol, Mpol, Ntor
   write(ounit,1014)            pscale, Ladiabatic, Lconstraint, mupftol, mupfits
   write(ounit,1015)            Lrad(1:min(Mvol,32))
   
1010 format("readin : ",f10.2," : Igeometry=",i3," ; pknot=",i3," ; Istellsym=",i3," ;")
1011 format("readin : ", 10x ," : Lfreebound=",i3," ; nextcur="i4" ; phiedge="es23.15" ; curtor="es23.15" ; curpol="es23.15" ;")
1012 format("readin : ", 10x ," : gamma="es23.15" ;")
1013 format("readin : ", 10x ," : Nfp=",i3," ; Nvol=",i3," ; Mpol=",i3," ; Ntor=",i3," ;")
1014 format("readin : ", 10x ," : pscale="es13.5" ; Ladiabatic="i2" ; Lconstraint="i2" ; mupf: tol,its="es13.5" ,"i4" ;")
1015 format("readin : ", 10x ," : Lrad = "32(i3,",")," ...")
   
#ifdef DEBUG
   if( Wreadin ) then
    write(ounit,'("readin : ",f10.2," : tflux    ="256(es11.3" ,":))') cput-cpus, (    tflux(vvol), vvol = 1, Nvol )
    write(ounit,'("readin : ",f10.2," : pflux    ="256(es11.3" ,":))') cput-cpus, (    pflux(vvol), vvol = 1, Nvol )
    write(ounit,'("readin : ",f10.2," : helicity ="256(es11.3" ,":))') cput-cpus, ( helicity(vvol), vvol = 1, Nvol )
    write(ounit,'("readin : ",f10.2," : pressure ="256(es11.3" ,":))') cput-cpus, ( pressure(vvol), vvol = 1, Nvol )
    write(ounit,'("readin : ",f10.2," : mu       ="256(es11.3" ,":))') cput-cpus, (       mu(vvol), vvol = 1, Nvol )
   endif
#endif
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
! screen output of namelist/numericlist/
   
   write(ounit,'("readin : ", 10x ," : ")')
   
   write(ounit,1020) cput-cpus, Linitialize, Lwall, phiwall, Ndiscrete
   write(ounit,1021)            Nquad, iMpol, iNtor
   write(ounit,1022)            Lsparse, Lsvdiota, imethod, iorder, iprecon, iotatol
   write(ounit,1023)            Iswmin, Lperturb, dperturb, Lextrap
   
1020 format("readin : ",f10.2," : Linitialize=",i3," ; Lwall=",i3," ; phiwall="es13.5" ; Ndiscrete="i2" ;")
1021 format("readin : ", 10x ," : Nquad="i4" ; iMpol="i4" ; iNtor="i4" ;")
1022 format("readin : ", 10x ," : Lsparse="i2" ; Lsvdiota="i2" ; imethod="i2" ; iorder="i2" ; iprecon="i2" ; iotatol="es13.5" ;")
1023 format("readin : ", 10x ," : Iswmin="i2" ; Lperturb="i2" ; dperturb="es13.5" ; Lextrap="i2" ;")
   
   if( phiwall.le.zero .or. phiwall.ge.one ) write(ounit,'("readin : ", 10x ," : phiwall="es13.5" is invalid; ")') phiwall
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
! screen output of namelist/locallist/
   
   write(ounit,'("readin : ", 10x ," : ")')
   
   write(ounit,1030) cput-cpus, LBeltrami, Linitgues, Lposdef, Nmaxexp
!  write(ounit,1031)            Lprecon, sparseits, ssoromega, sparsetol
!  write(ounit,1032)            Liotasolv
   
1030 format("readin : ",f10.2," : LBeltrami="i2" ; Linitgues="i2" ; Lposdef="L2" ; Nmaxexp="i3" ;")
!031 format("readin : ", 10x ," : Lprecon="i2" ; sparseits="i4" ; ssoromega="f6.3" ; sparsetol="es13.5" ;")
!032 format("readin : ", 10x ," : Liotasolv="i2" ;")
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
! screen output of namelist/globallist/
   
   write(ounit,'("readin : ", 10x ," : ")')
   
   write(ounit,1040) cput-cpus, Lminimize, Lfindzero
   write(ounit,1041)            escale, pcondense, wpoloidal
   write(ounit,1042)            forcetol, normalerr, norblend, maxfbits
   write(ounit,1043)            c05xtol, c05factor
   write(ounit,1044)            LreadGF, verify, maxstep, epsilon, upsilon, opsilon, maxiter
   
1040 format("readin : ",f10.2," : Lminimize="i2" ; Lfindzero="i2" ;")
1041 format("readin : ", 10x ," : escale="es13.5" ; pcondense="f7.3" ; wpoloidal="f7.4" ;")
1042 format("readin : ", 10x ," : forcetol="es13.5" ; normalerr="es13.5" ; norblend="es13.5" ; maxfbits="i4" ;")
1043 format("readin : ", 10x ," : c05xtol="es13.5" ; c05factor="es13.5" ;")
1044 format("readin : ", 10x ," : LreadGF="L2" ; verify=",i3," ; maxstep="es13.5" ; epsilon="es13.5" ; upsilon="es13.5" ; opsilon="es13.5" ; maxiter="i9" ;")
   
   FATALMESS( readin, pcondense.lt.one, invalid)
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
! screen output of namelist/diagnosticslist/
   
   write(ounit,'("readin : ", 10x ," : ")')
   
   write(ounit,1050) cput-cpus, odetol, absreq, relreq, absacc, epsr, nPpts
!  write(ounit,1051)            irrMpol, irrNtor, irrsvdcut, irrtol, Lwrpj
   write(ounit,1052)            LHevalues, LHevectors, Lperturbed, dpp, dqq, Lcheck, Ltiming
   
1050 format("readin : ",f10.2," : odetol="es10.2" ; absreq="es10.2" ; relreq="es10.2" ; absacc="es10.2" ; epsr="es10.2" ; nPpts="i6" ;")
!051 format("readin : ", 10x ," : irrMpol="i4" ; irrNtor="i4" ; irrsvdcut="es13.5" ; irrtol="es13.5" ; Lwrpj="L2" ;")
1052 format("readin : ", 10x ," : LHevalues="L2" ; LHevectors="L2" ; Lperturbed="i2" ; dpp="i2" ; dqq="i3" ; Lcheck="i3" ; Ltiming="L2" ;")
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
! screen output of namelist/screenlist/; no screen output;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   write(ounit,'("readin : ", 10x ," : ")')
   
  endif ! end of if myid eq 0 loop;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! broadcast command line input
  
  ClBCAST(ext,100,0)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! broadcast namelist/physicslist/
  
  if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : broadcasting physicslist     from ext.spec ;")') cput-cpus
  endif
  
  IlBCAST(Igeometry,1,0)
  IlBCAST(pknot,1,0)
  IlBCAST(Istellsym,1,0)
  IlBCAST(Lfreebound,1,0) ! 10 Oct 12; 
  RlBCAST(phiedge,1,0)
  RlBCAST(curtor,1,0)
  RlBCAST(curpol,1,0)
  RlBCAST(extcur,maxgroups,0)
  RlBCAST(gamma,1,0)
  IlBCAST(Nfp,1,0)
  IlBCAST(Nvol,1,0)
  IlBCAST(Mpol,1,0)
  IlBCAST(Ntor,1,0)
  IlBCAST(Lrad,MNvol,0)
  RlBCAST(tflux,MNvol,0)
  RlBCAST(pflux,MNvol,0)
  RlBCAST(helicity,MNvol,0)
  RlBCAST(pscale,1,0)
  RlBCAST(pressure,MNvol,0)
  IlBCAST(Ladiabatic,1,0)
  RlBCAST(adiabatic,MNvol,0)
  RlBCAST(mu,MNvol,0)
  IlBCAST(Lconstraint,1,0)
  IlBCAST(pl,MNvol,0)
  IlBCAST(ql,MNvol,0)
  IlBCAST(pr,MNvol,0)
  IlBCAST(qr,MNvol,0)
  RlBCAST(iota,MNvol,0)
  IlBCAST(lp,MNvol,0)
  IlBCAST(lq,MNvol,0)
  IlBCAST(rp,MNvol,0)
  IlBCAST(rq,MNvol,0)
  RlBCAST(oita,MNvol,0)
  RlBCAST(mupftol,1,0)
  IlBCAST(mupfits,1,0)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! broadcast namelist/numericlist/
  
  if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : broadcasting numericlist     from ext.spec ;")') cput-cpus
  endif
  
  IlBCAST(Linitialize,1,0)
  IlBCAST(Lwall,1,0)
  RlBCAST(phiwall,1,0)
  IlBCAST(Ndiscrete,1,0)
  IlBCAST(Nquad,1,0)
  IlBCAST(iMpol,1,0)
  IlBCAST(iNtor,1,0)
  IlBCAST(Lsparse,1,0)
  IlBCAST(Lsvdiota,1,0)
  IlBCAST(imethod,1,0)
  IlBCAST(iorder,1,0)
  IlBCAST(iprecon,1,0)
  RlBCAST(iotatol,1,0)
  IlBCAST(Iswmin,1,0)
  IlBCAST(Lperturb,1,0)
  RlBCAST(dperturb,1,0)
  IlBCAST(Lextrap,1,0)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! broadcast namelist/globallist/
  
  if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : broadcasting globallist      from ext.spec ;")') cput-cpus
  endif
  
  IlBCAST(Lminimize,1,0)
  IlBCAST(Lfindzero,1,0)
 !IlBCAST(Lcondense,1,0)
  RlBCAST(escale,1,0)
  RlBCAST(pcondense,1,0)
 !RlBCAST(qcondense,1,0) ! 04 Dec 14;
  RlBCAST(wpoloidal,1,0)
  RlBCAST(forcetol,1,0)
  RlBCAST(normalerr,1,0)
  RlBCAST(norblend,1,0)
  IlBCAST(maxfbits,1,0)
 !RlBCAST(ForceErr,1,0) ! this is an ``output'' quantity;
 !RlBCAST(Energy,1,0)   ! this is an ``output'' quantity;
  RlBCAST(c05xtol,1,0)
  RlBCAST(c05factor,1,0)
  LlBCAST(LreadGF,1,0)
  IlBCAST(verify,1,0)
  RlBCAST(maxstep,1,0)
  RlBCAST(epsilon,1,0) ! 04 Dec 14;
  RlBCAST(upsilon,1,0)
  RlBCAST(opsilon,1,0)
 !RlBCAST(apsilon,1,0) ! 18 Jul 14;
  IlBCAST(maxiter,1,0)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! broadcast namelist/locallist/
  
  if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : broadcasting locallist       from ext.spec ;")') cput-cpus
  endif
  
  IlBCAST(LBeltrami,1,0)
  IlBCAST(Linitgues,1,0)
  LlBCAST(Lposdef,1,0)
  IlBCAST(Nmaxexp,1,0)
! IlBCAST(Lprecon,1,0)
! IlBCAST(sparseits,1,0)
! RlBCAST(ssoromega,1,0)
! RlBCAST(sparsetol,1,0)
! IlBCAST(Liotasolv,1,0)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! broadcast namelist/diagnosticslist/
  
  if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : broadcasting diagnosticslist from ext.spec ;")') cput-cpus
  endif
  
  RlBCAST(odetol,1,0)
  RlBCAST(absreq,1,0)
  RlBCAST(relreq,1,0)
  RlBCAST(absacc,1,0)
  RlBCAST(epsr,1,0)
 !RlBCAST(divertorR,1,0)
 !RlBCAST(divertorZ,1,0)
  IlBCAST(nPpts,1,0)
  IlBCAST(nPtrj,(MNvol+1),0)
  IlBCAST(Mpqits,1,0)
  IlBCAST(Lpqsym,1,0)
  IlBCAST(p1,(MNvol),0)
  IlBCAST(q1,(MNvol),0)
  IlBCAST(p2,(MNvol),0)
  IlBCAST(q2,(MNvol),0)
  RlBCAST(pqs,(MNvol),0)
  RlBCAST(pqt,(MNvol),0)
! IlBCAST(Npl,1,0)
! IlBCAST(Nunstable,(MNvol),0)!
! RlBCAST(dunstable,(MNvol),0)
! IlBCAST(Munstable,(MNvol),0)
  IlBCAST(npq,(MNvol),0)
! IlBCAST(Mirrits,1,0)
! IlBCAST(irrMpol,1,0)
! IlBCAST(irrNtor,1,0)
! RlBCAST(irrsvdcut,1,0)
! RlBCAST(irrtol,1,0)
! LlBCAST(Lwrpj,1,0)
! IlBCAST(Nghd,1,0)
! LlBCAST(LHessian,1,0)
  LlBCAST(LHevalues,1,0)
  LlBCAST(LHevectors,1,0)
  IlBCAST(Lperturbed,1,0)
  IlBCAST(dpp,1,0)
  IlBCAST(dqq,1,0)
! LlBCAST(Lcurlerr,1,0)
  IlBCAST(Lcheck,1,0)
  LlBCAST(Ltiming,1,0)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! broadcast namelist/screenlist/
  
  if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : broadcasting screenlist      from ext.spec ;")') cput-cpus
  endif
  
! BSCREENLIST ! broadcast screenlist; this is expanded by Makefile; do not remove;
  LlBCAST(Wreadin,1,0)
  LlBCAST(Wwritin,1,0)
  LlBCAST(Wmacros,1,0)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! set internal parameters that depend on physicslist;

  RALLOCATE( beltramierror,(1:Nvol,1:3) ) ! 11 Apr 16;  ! replaced Mvol with Nvol; 11 Jul 17;

  select case( Istellsym )
  case( 0 ) ; NOTstellsym = .true. ; YESstellsym = .false.
  case( 1 ) ; YESstellsym = .true. ; NOTstellsym = .false.
  case default
   FATALMESS(readin,.true.,invalid Istellsym)
  end select
  
  if( Lfreebound.eq.0 ) Mvol = Nvol
  if( Lfreebound.gt.0 ) Mvol = Nvol + 1 ! include vacuum region; 04 April 13;
  
  mn  =  1 + Ntor +  Mpol * ( 2 *  Ntor + 1 ) ! Fourier resolution of interface geometry & vector potential;
  
#ifdef DEBUG
  if( Llatex ) then
   write(wunit,'("\begin{{enumerate}}")')
   write(wunit,'("\item [] \verb+mn+ $ = "i4"$")') mn
   write(wunit,'("\end{{enumerate}}")')
  endif
#endif
  
  IALLOCATE(im,(1:mn))
  IALLOCATE(in,(1:mn))
  
  call gi00ab(  Mpol,  Ntor, Nfp, mn, im, in  ) ! this sets the im and in mode identification arrays;
  
  select case( Igeometry )
  case( 1   ) 
  case( 2:3 )  ; RALLOCATE(halfmm,(1:mn))
   ;           ! RALLOCATE(regularmm,(1:mn))
   ;           ; do ii = 1, mn
   ;           ;     halfmm(ii) =               im(ii)   * half
   ;           !  regularmm(ii) = min( Nmaxexp, im(ii) ) * half
   ;           ; enddo
  case default ; FATALMESS(readin, .true., invalid Igeometry)
  end select
  
#ifdef DEBUG
  if( Wreadin .and. myid.eq.0 ) then
   write(ounit,'("readin : ", 10x ," : myid=",i3," ;    halfmm ="99(f5.1","))') myid,    halfmm(1:mn)
!  write(ounit,'("readin : ", 10x ," : myid=",i3," ; regularmm ="99(f5.1","))') myid, regularmm(1:mn)
  endif
#endif
  
! lMpol =   Mpol ; lNtor =   Ntor ! no    enhanced resolution for metrics; 14 Apr 13;
! lMpol = 2*Mpol ; lNtor = 2*Ntor !           enhanced resolution for metrics;
  lMpol = 4*Mpol ; lNtor = 4*Ntor ! extra-enhanced resolution for metrics; 14 Apr 13;
  
  mne = 1 + lNtor + lMpol * ( 2 * lNtor + 1 ) ! resolution of metrics; enhanced resolution; see me00ab;

  IALLOCATE(ime,(1:mne))
  IALLOCATE(ine,(1:mne))

  call gi00ab( lMpol, lNtor, Nfp, mne, ime, ine )
  
#ifdef DEBUG
  if( Wreadin .and. myid.eq.0 ) then
   write(ounit,'("readin : ", 10x ," : myid=",i3," ; ime="99i4)') myid, ime(1:mne)
   write(ounit,'("readin : ", 10x ," :      "3x" ; ine="99i4)')       ine(1:mne)
  endif
#endif
  
  lMpol = iMpol ; lNtor = iNtor

  if( iMpol.le.0 ) lMpol = Mpol - iMpol
  if( iNtor.le.0 ) lNtor = Ntor - iNtor
  if(  Ntor.eq.0 ) lNtor = 0
  
  mns = 1 + lNtor + lMpol * ( 2 * lNtor + 1 ) ! resolution of straight-field line transformation on interfaces; see tr00ab; soon to be redundant; 21 Apr 13;

  IALLOCATE(ims,(1:mns))
  IALLOCATE(ins,(1:mns))

  call gi00ab( lMpol, lNtor, Nfp, mns, ims, ins ) ! note that the field periodicity factor is included in ins;
  
#ifdef DEBUG
  if( Wreadin .and. myid.eq.0 ) then
   write(ounit,'("readin : ", 10x ," : ims = ",999i3)') ims(1:mns)
   write(ounit,'("readin : ", 10x ," : ins = ",999i3)') ins(1:mns)
   write(ounit,'("readin : ", 10x ," : myid=",i3," : tflux(1:Nvol)=",256es13.5)') myid, tflux(1:Nvol)
   write(ounit,'("readin : ", 10x ," : myid=",i3," : tflux(1)**(m/2)=",256es13.5)') myid, tflux(1)**(/ ( mm*half, mm = 0, Mpol ) /)
  endif
#endif
  
  
  RALLOCATE(psifactor,(1:mn,1:Nvol))
  
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
     if( im(ii).eq.0 ) then ; psifactor(ii,vvol) = tflux(vvol)**half       ! 29 Apr 15;
     else                   ; psifactor(ii,vvol) = tflux(vvol)**halfmm(ii) ! 29 Apr 15;
     endif
    enddo
   enddo
   
  case default
   
   FATALMESS(al00aa, .true., invalide Igeometry for construction of psifactor)
   
  end select
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! set internal parameters that depend on numericlist;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! set internal parameters that depend on locallist;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! set internal parameters that depend on globallist;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! set internal parameters that depend on diagnosticslist;
  
  if( Lcheck.eq.5 ) then ! will check Hessian using finite-differences; 18 Jul 14;
   forcetol = 1.0e+12
   nPpts    = 0
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RALLOCATE( iRbc,(1:mn,0:Mvol) ) ! interface Fourier harmonics;
  RALLOCATE( iZbs,(1:mn,0:Mvol) )
! if( NOTstellsym ) then ! these are always allocated, but not always used; for file consistency, they are written to output; 02 Sep 14;
  RALLOCATE( iRbs,(1:mn,0:Mvol) )
  RALLOCATE( iZbc,(1:mn,0:Mvol) )
! endif
  
  if( Lperturbed.eq.1 ) then
  RALLOCATE( dRbc,(1:mn,0:Mvol) ) ! interface Fourier harmonics;
  RALLOCATE( dZbs,(1:mn,0:Mvol) )
! if( NOTstellsym ) then ! these are always allocated, but not always used; for file consistency, they are written to output; 02 Sep 14;
  RALLOCATE( dRbs,(1:mn,0:Mvol) )
  RALLOCATE( dZbc,(1:mn,0:Mvol) )
! endif
  endif
  
  RALLOCATE( iBns,(1:mn) ) ! non-zero normal field at outermost interface;
  RALLOCATE( iBnc,(1:mn) )
  
  RALLOCATE( iCns,(1:mn) ) ! coil   normal field at computational boundary;
  RALLOCATE( iCnc,(1:mn) )
  
  RALLOCATE( iPns,(1:mn) ) ! plasma normal field at computational boundary;
  RALLOCATE( iPnc,(1:mn) )
  
  RALLOCATE( lRbc,(1:mn) ) ! interface Fourier harmonics; 18 Apr 13; local workspace; perhaps redundant; 14 Apr 13;
  RALLOCATE( lZbs,(1:mn) )
  RALLOCATE( lRbs,(1:mn) ) ! interface Fourier harmonics;
  RALLOCATE( lZbc,(1:mn) )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  IALLOCATE(bjk,(1:mn,1:mn)) ! this must be allocated & assigned now, as it is used in readin; primarily used in gf00aa; 02 Jan 15;
  
  do kk = 1, mn ; mk = im(kk) ; nk = in(kk)
   
   if( mk.ne.0 ) cycle
   
   do jj = 1, mn ; mj = im(jj) ; nj = in(jj)
    
    if( (mj/2)*2.eq.mj ) then ! mj is even; 11 Aug 14;
     
     if( nj.eq.-nk) bjk(jj,kk) = -1
     if( nj.eq. nk) bjk(jj,kk) = +1 ! this overrides 0.eq.0 case; 11 Aug 14;
     
    endif ! end of if( mj.is.even ) ; 11 Aug 14;
    
   enddo ! end of do jj; 11 Aug 14;
   
  enddo ! end of do kk; 11 Aug 14;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( myid.eq.0 ) then ! initialize interface geometry;
   
   do ii = 1, mn ; mm = im(ii) ; nn = in(ii) / Nfp ! set plasma boundary, computational boundary; Modified 20 Jan 2016;
    
    if( Igeometry.eq.3 .and. Rbc(0,+1)+Rbc(0,-1).gt.zero .and. Zbs(0,+1)-Zbs(0,-1).gt.zero ) then
     jj = -1; kk = -nn
    else
     jj = +1; kk = +nn
    endif
    
    if(mm.eq.0 .and. nn.eq.0) then
     ;iRbc(ii,Nvol) =   Rbc(nn,mm)    
     ;iZbs(ii,Nvol) =   zero
     ;iRbs(ii,Nvol) =   zero
     ;iZbc(ii,Nvol) =   Zbc(nn,mm)
    else
     ;iRbc(ii,Nvol) =   Rbc(kk,mm) + Rbc(-kk,-mm)    
     ;iZbs(ii,Nvol) =   (Zbs(kk,mm) - Zbs(-kk,-mm)) * jj
     ;iRbs(ii,Nvol) =   (Rbs(kk,mm) - Rbs(-kk,-mm)) * jj
     ;iZbc(ii,Nvol) =   Zbc(kk,mm) + Zbc(-kk,-mm)    
    endif 
  
   enddo ! end of do ii = 1, mn;
   
   
! --- Commented on 20 Jan 2016 ---
!   if( Igeometry.eq.3 .and. Zbs(0,1).gt.zero ) then ! change sign of poloidal angle; 11 Aug 14;
!    
!    write(ounit,'("readin : ", 10x ," : myid=",i3," ; CHANGING POLOIDAL ANGLE ;")') myid
!    write(ounit,'("readin : ", 10x ," : ")')
!    
!    do ii = 1, mn ; mm = im(ii) ; nn = in(ii) / Nfp
!     
!     ;if( mm.eq.0 .or.  nn.eq.0 ) then ; iRbc(ii,Nvol) =   Rbc( nn,mm)
!     ;else                             ; iRbc(ii,Nvol) =   Rbc(-nn,mm)
!     ;endif
!     ;if( mm.eq.0 .and. nn.eq.0 ) then ; iZbs(ii,Nvol) =   Zbs( nn,mm)
!     ;else                             ; iZbs(ii,Nvol) = - Zbs(-nn,mm)
!     ;endif
!     
!     if( NOTstellsym ) then
!      if( mm.eq.0 .or.  nn.eq.0 ) then ; iZbc(ii,Nvol) =   Zbc( nn,mm)
!      else                             ; iZbc(ii,Nvol) =   Zbc(-nn,mm)
!      endif
!      if( mm.eq.0 .and. nn.eq.0 ) then ; iRbs(ii,Nvol) =   Rbs( nn,mm)
!      else                             ; iRbs(ii,Nvol) = - Rbs(-nn,mm)
!      endif
!     endif
!     
!    enddo ! end of do ii; 11 Aug 14;
!   
!   endif ! end of if( Igeometry.eq.3 .and. Zbs(0,1).lt.zero ) then ; 11 Aug 14;
! --------------------------
   

#ifdef DEBUG
   if( Wreadin ) then
    ; write(ounit,'("readin : ", 10x ," : myid=",i3," : Rbc(1:mn,",i3,")=",99es13.5)') myid, Nvol, iRbc(1:mn,Nvol)
    ; write(ounit,'("readin : ", 10x ," : myid=",i3," : Zbs(1:mn,",i3,")=",99es13.5)') myid, Nvol, iZbs(1:mn,Nvol)
    if( NOTstellsym ) then
     ;write(ounit,'("readin : ", 10x ," : myid=",i3," : Rbs(1:mn,",i3,")=",99es13.5)') myid, Nvol, iRbs(1:mn,Nvol)
     ;write(ounit,'("readin : ", 10x ," : myid=",i3," : Zbc(1:mn,",i3,")=",99es13.5)') myid, Nvol, iZbc(1:mn,Nvol)
    endif
    ; write(ounit,'("readin : ", 10x ," : myid=",i3," : Rbc(1:mn,",i3,")=",99es13.5)') myid, Mvol, iRbc(1:mn,Mvol)
    ; write(ounit,'("readin : ", 10x ," : myid=",i3," : Zbs(1:mn,",i3,")=",99es13.5)') myid, Mvol, iZbs(1:mn,Mvol)
    if( NOTstellsym ) then
     ;write(ounit,'("readin : ", 10x ," : myid=",i3," : Rbs(1:mn,",i3,")=",99es13.5)') myid, Mvol, iRbs(1:mn,Mvol)
     ;write(ounit,'("readin : ", 10x ," : myid=",i3," : Zbc(1:mn,",i3,")=",99es13.5)') myid, Mvol, iZbc(1:mn,Mvol)
    endif
    if( Lfreebound.gt.0 ) then
     ;write(ounit,'("readin : ", 10x ," : myid=",i3," : Bns(1:mn ",3x,")=",99es13.5)') myid      , iBns(1:mn     )
     if( NOTstellsym ) then
      write(ounit,'("readin : ", 10x ," : myid=",i3," : Bnc(1:mn ",3x,")=",99es13.5)') myid      , iBnc(1:mn     )
     endif
    endif
   endif
#endif
   
   if( Rac(0).lt.small ) then ! inside myid.eq.0; supplied coordinate axis is invalid; assume user wants reasonable axis to be provided; 13 Sep 13;
    
    vvol = 0
    
    do nb = 0, Ntor
     
     do mm = 0, Mpol, 2 ! ensure mm is even; 11 Aug 14;
      do nn = -Ntor, Ntor
       if(                   nn.eq. nb ) then
        ; iRbc(nb+1,vvol) = iRbc(nb+1,vvol) + Rbc(nn,mm)
        if( Igeometry.eq.3 ) then
         ;iZbs(nb+1,vvol) = iZbs(nb+1,vvol) + Zbs(nn,mm)
        endif
        if( NOTstellsym ) then
         ;iRbs(nb+1,vvol) = iRbs(nb+1,vvol) + Rbs(nn,mm)
         if( Igeometry.eq.3 ) then
          iZbc(nb+1,vvol) = iZbc(nb+1,vvol) + Zbc(nn,mm)
         endif
        endif ! end of if( NOTstellsym ) ; 11 Aug 14;
       elseif( nb.ne.0 .and. nn.eq.-nb ) then
        ; iRbc(nb+1,vvol) = iRbc(nb+1,vvol) + Rbc(nn,mm)
        if( Igeometry.eq.3 ) then
         ;iZbs(nb+1,vvol) = iZbs(nb+1,vvol) - Zbs(nn,mm)
        endif
        if( NOTstellsym ) then
         ;iRbs(nb+1,vvol) = iRbs(nb+1,vvol) - Rbs(nn,mm)
         if( Igeometry.eq.3 ) then
          iZbc(nb+1,vvol) = iZbc(nb+1,vvol) + Zbc(nn,mm)
         endif
        endif ! end of if( NOTstellsym ) ; 11 Aug 14;
       endif
      enddo ! end of do nn; 11 Aug 14;
     enddo ! end of do mm; 11 Aug 14;
     
    enddo ! end of do nb; 11 Aug 14;
    
   else ! matches if( Rac(0).lt.small) ; 29 Apr 15;
    
    vvol = 0
    
    ; iRbc(1:Ntor+1,vvol) = Rac(0:Ntor) ! this is inside if( myid.eq.0 ) ; set location of degenerate interface; initialize interface geometry; 11 Aug 14;
    if( Igeometry.eq.3 ) then
     ;iZbs(2:Ntor+1,vvol) = Zas(1:Ntor)
    endif
    if( NOTstellsym ) then
     ;iRbs(2:Ntor+1,vvol) = Ras(1:Ntor)
     if( Igeometry.eq.3 ) then
      iZbc(1:Ntor+1,vvol) = Zac(0:Ntor)
     endif
    endif
    
   endif ! end of if( Rac(0).lt.small ) ; 13 Sep 13;
   

!  if( Wreadin ) then 
!   cput = GETTIME
!   ; write(ounit,'("readin : ",f10.2," : myid=",i3," : Rac =",99es13.5)') cput-cpus, myid, Rac(0:Ntor)
!   if( Igeometry.eq.3 ) then
!    ;write(ounit,'("readin : ",f10.2," : myid=",i3," : Zas =",99es13.5)') cput-cpus, myid, Zas(0:Ntor)
!   endif
!   if( NOTstellsym ) then    
!    ;write(ounit,'("readin : ",f10.2," : myid=",i3," : Ras =",99es13.5)') cput-cpus, myid, Ras(0:Ntor)
!    if( Igeometry.eq.3 ) then
!     write(ounit,'("readin : ",f10.2," : myid=",i3," : Zac =",99es13.5)') cput-cpus, myid, Zac(0:Ntor)
!    endif
!   endif
!  endif

   
   select case( Linitialize ) ! 24 Oct 12;
    
   case( 0 ) ! initial guess for geometry of the interior surfaces is given in the input file;
    
    RALLOCATE(RZRZ,(1:4,1:Nvol)) ! temp array for reading input;
    
    do ! will read in Fourier harmonics until the end of file is reached;
     
     read(iunit,*,iostat=ios) mm, nn, RZRZ(1:4,1:Nvol)
     if( ios.ne.0 ) exit
     
     do ii = 1, mn ; mi = im(ii) ; ni = in(ii) ! loop over harmonics within range;
      if( mm.eq.mi .and. nn*Nfp.eq.ni ) then
       ; iRbc(ii,1:Nvol-1) = RZRZ(1,1:Nvol-1) ! select relevant harmonics;
       if( Igeometry.eq.3 ) then
        ;iZbs(ii,1:Nvol-1) = RZRZ(2,1:Nvol-1) ! select relevant harmonics;
       endif
       if( NOTstellsym ) then
        ;iRbs(ii,1:Nvol-1) = RZRZ(3,1:Nvol-1) ! select relevant harmonics;
        if( Igeometry.eq.3 ) then
         iZbc(ii,1:Nvol-1) = RZRZ(4,1:Nvol-1) ! select relevant harmonics;
        endif
       endif
      endif
     enddo ! end of do ii;
     
    enddo ! end of do;
    
    DEALLOCATE(RZRZ)
    
   case( 1 ) ! Linitialize=1; initialize internal geometry using regular extrapolation;
    
    select case( Igeometry ) ! initialization extrapolation depends on geometry; 20 Apr 13;
     
    case( 1 ) ! Cartesian; 29 Apr 14;
     
     do vvol = 1, Nvol-1
      ;iRbc(1:mn,vvol) = iRbc(1:mn,Nvol) * tflux(vvol)
      if( NOTstellsym ) then
       iRbs(2:mn,vvol) = iRbs(2:mn,Nvol) * tflux(vvol)
      endif
     enddo
     
    case( 2 ) ! cylindrical - standard; 20 Apr 13;
     
     do vvol = 1, Nvol-1
      
      ;iRbc(1:mn,vvol) = iRbc(1:mn,Nvol) * psifactor(1:mn,vvol)
      if( NOTstellsym ) then
       iRbs(2:mn,vvol) = iRbs(2:mn,Nvol) * psifactor(2:mn,vvol)
      endif
     enddo
     
    case( 3 ) ! toroidal; 20 Apr 13;
     
     do vvol = 1, Nvol-1
      
      ; iRbc(1:mn,vvol) = iRbc(1:mn,0) + ( iRbc(1:mn,Nvol) - iRbc(1:mn,0) ) * psifactor(1:mn,vvol)
      if( Igeometry.eq.3 ) then
       ;iZbs(2:mn,vvol) = iZbs(2:mn,0) + ( iZbs(2:mn,Nvol) - iZbs(2:mn,0) ) * psifactor(2:mn,vvol)
      endif
      if( NOTstellsym ) then
       ;iRbs(2:mn,vvol) = iRbs(2:mn,0) + ( iRbs(2:mn,Nvol) - iRbs(2:mn,0) ) * psifactor(2:mn,vvol)
       if( Igeometry.eq.3 ) then
        iZbc(1:mn,vvol) = iZbc(1:mn,0) + ( iZbc(1:mn,Nvol) - iZbc(1:mn,0) ) * psifactor(1:mn,vvol)
       endif
      endif
      
     enddo
          
    case default ! matches select case( Igeometry ) ; 01 May 13;
     
     FATALMESS(readin, .true., Linitialize is not supported for this Igeometry)
     
    end select ! matches select case( Igeometry ) ; 01 May 13;
    
   case default ! matches select case( Linitialize ) ; 01 May 13;
    
    FATALMESS(readin, .true., invalid Linitialize)
    
   end select ! matches select case( Linitialize ) ; 01 May 13;
   
   
   close(iunit) ! finished reading input file; 13 Sep 13;

   
   if( Rac(0).lt.small .and. Linitialize.eq.0 ) then ! presumably the coordinate axis has not been given; calculate default; 02 Jan 15;
    
    do kk = 1, mn ! see also global:writin; 11 Aug 14;
     ; iRbc(kk,0) = sum( iRbc(1:mn,1) * abs( bjk(1:mn,kk) ) )
     if( Igeometry.eq.3 ) then
      ;iZbs(kk,0) = sum( iZbs(1:mn,1) *      bjk(1:mn,kk)   )
     endif
     if( NOTstellsym ) then
      ;iRbs(kk,0) = sum( iRbs(1:mn,1) *      bjk(1:mn,kk)   )
      if( Igeometry.eq.3 ) then
       iZbc(kk,0) = sum( iZbc(1:mn,1) * abs( bjk(1:mn,kk) ) )
      endif
     endif
    enddo
    
   endif ! end of if( Rac(0).lt.small .and. Linitialize.eq.0 ) then ; 11 Aug 14;
   
#ifdef OLDNAG

   if( Lperturb.eq.1 ) then ! add random perturbation to interior interface Fourier harmonics;
    
    do vvol = 1, Mvol-1
     
     do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
      
      ; xx = G05CAF(xx) ; iRbc(ii,vvol) = iRbc(ii,vvol) + ( xx - half ) * dperturb / max(1,mi)
      if( Igeometry.eq.3 ) then
       ;xx = G05CAF(xx) ; iZbs(ii,vvol) = iZbs(ii,vvol) + ( xx - half ) * dperturb / max(1,mi) ! only required for ii > 1;
      endif
      if( NOTstellsym ) then
       ;xx = G05CAF(xx) ; iRbs(ii,vvol) = iRbs(ii,vvol) + ( xx - half ) * dperturb / max(1,mi) ! only required for ii > 1;
       if( Igeometry.eq.3 ) then
        xx = G05CAF(xx) ; iZbc(ii,vvol) = iZbc(ii,vvol) + ( xx - half ) * dperturb / max(1,mi)
       endif
      endif
      
     enddo ! end of do ii; 26 Feb 13;
     
    enddo ! end of do vvol; 26 Feb 13;
    
   endif ! end of if( Lperturb) ; 26 Feb 13;
   
#endif

  endif ! end of if myid.eq.0 loop; only the master will read the input file; all variables need to be broadcast;
  
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
   ;RlBCAST( iBns(1:mn), mn, 0 ) ! only required for ii > 1 ;
   if( NOTstellsym ) then
    RlBCAST( iBnc(1:mn), mn, 0 )
   endif
  endif
  
  if( Igeometry.eq.1 .or. Igeometry.eq.2 ) then
   ;iRbc(1:mn,0) = zero ! innermost volume must be trivial; this is used in vo00aa; 26 Feb 13; innermost interface is coordinate axis; 13 Sep 13;
   if( NOTstellsym ) then
    iRbs(1:mn,0) = zero ! innermost volume must be trivial; this is used in vo00aa; 26 Feb 13;
   endif
  endif
  
  if( Igeometry.eq.3 ) then
   iZbs(1,0:Mvol) = zero ! Zbs_{m=0,n=0} is irrelevant;
  endif
  if( NOTstellsym) then
   iRbs(1,0:Mvol) = zero ! Rbs_{m=0,n=0} is irrelevant;
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(readin)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine readin

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine writin( wflag, iflag, rflag ) !latex \end{enumerate} \subsubsection{writin} \begin{enumerate}

!latex \item The restart file is written.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  use constants, only :

  use numerical, only : machprec

  use fileunits, only : ounit, iunit, zunit

  use cputiming, only : Twritin

  use inputlist

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS

  INTEGER, intent(in)  :: wflag, iflag
  REAL   , intent(in)  :: rflag
  
  LOGICAL              :: Lverbose = .true. ! will eventually have this as user input; set Lverbose = .false. to have concise restart file;

  INTEGER              :: vvol, imn, ii, jj, kk, Lcurvature, mm, nn
  REAL                 :: lss
  
  BEGIN(writin)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( myid.ne.0 ) goto 9999

  FATALMESS(writin, myid.ne.0, error)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
 !WCALL(writin,hdfint) ! 18 Jul 14;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Wwritin ) then ; cput = GETTIME ; write(ounit,'("writin : ",f10.2," : myid=",i3," ; opening/writing ext.end ;")') cput-cpus, myid
  endif

  open(iunit,file=trim(ext)//".end",status="unknown") ! restart input file;

  if( Wwritin ) then ; cput = GETTIME ; write(ounit,'("writin : ",f10.2," : myid=",i3," ; writing physicslist ;")') cput-cpus, myid
  endif
  
  write(iunit,'("&physicslist")')
  write(iunit,'(" Igeometry   = ",i9        )') Igeometry
  write(iunit,'(" pknot       = ",i9        )') pknot
  write(iunit,'(" Istellsym   = ",i9        )') Istellsym
  write(iunit,'(" Lfreebound  = ",i9        )') Lfreebound
  write(iunit,'(" phiedge     = ",es23.15   )') phiedge
  write(iunit,'(" curtor      = ",es23.15   )') curtor
  write(iunit,'(" curpol      = ",es23.15   )') curpol
  write(iunit,'(" extcur      = ",999es23.15)') extcur(1:max(nextcur,1))
  write(iunit,'(" gamma       = ",es23.15   )') gamma
  write(iunit,'(" Nfp         = ",i9        )') Nfp
  write(iunit,'(" Nvol        = ",i9        )') Nvol
  write(iunit,'(" Mpol        = ",i9        )') Mpol
  write(iunit,'(" Ntor        = ",i9        )') Ntor
  write(iunit,'(" Lrad        = ",256i23    )') Lrad(1:Mvol)
  write(iunit,'(" tflux       = ",256es23.15)') tflux(1:Nvol)
  write(iunit,'(" pflux       = ",256es23.15)') pflux(1:Nvol)
  write(iunit,'(" helicity    = ",256es23.15)') helicity(1:Nvol)
  write(iunit,'(" pscale      = ",es23.15   )') pscale
  write(iunit,'(" Ladiabatic  = ",i9        )') Ladiabatic
  write(iunit,'(" pressure    = ",256es23.15)') pressure(1:Nvol)
  write(iunit,'(" adiabatic   = ",256es23.15)') adiabatic(1:Nvol)
  write(iunit,'(" mu          = ",256es23.15)') mu(1:Nvol)
  write(iunit,'(" Lconstraint = ",i9        )') Lconstraint
  write(iunit,'(" pl          = ",257i23    )') pl(0:Nvol)
  write(iunit,'(" ql          = ",257i23    )') ql(0:Nvol)
  write(iunit,'(" pr          = ",257i23    )') pr(0:Nvol)
  write(iunit,'(" qr          = ",257i23    )') qr(0:Nvol)
  write(iunit,'(" iota        = ",257es23.15)') iota(0:Nvol)
  write(iunit,'(" lp          = ",257i23    )') lp(0:Nvol)
  write(iunit,'(" lq          = ",257i23    )') lq(0:Nvol)
  write(iunit,'(" rp          = ",257i23    )') rp(0:Nvol)
  write(iunit,'(" rq          = ",257i23    )') rq(0:Nvol)
  write(iunit,'(" oita        = ",257es23.15)') oita(0:Nvol)
  write(iunit,'(" mupftol     = ",es23.15   )') mupftol
  write(iunit,'(" mupfits     = ",i9        )') mupfits

  if( Lfreebound.gt.0 .or. Zbs(0,1).gt.zero ) then
   do ii = 1, mn ; mm = im(ii) ; nn = in(ii) / Nfp ; Rbc(nn,mm) = iRbc(ii,Nvol) ; Zbs(nn,mm) = iZbs(ii,Nvol) ; Bns(nn,mm) = iBns(ii)
    if( NOTstellsym ) then                         ; Rbs(nn,mm) = iRbs(ii,Nvol) ; Zbc(nn,mm) = iZbc(ii,Nvol) ; Bnc(nn,mm) = iBnc(ii)
    endif
   enddo ! end of do ii = 1, mn;
  endif ! end of if( Lfreebound ) ; 15 May 13;

  write(iunit,'(" Rac         = ",99es23.15)') Rac(0:Ntor)
  write(iunit,'(" Zas         = ",99es23.15)') Zas(0:Ntor)
  write(iunit,'(" Ras         = ",99es23.15)') Ras(0:Ntor)
  write(iunit,'(" Zac         = ",99es23.15)') Zac(0:Ntor)
  
 !write(iunit,'(" Rac         = ",99es23.15)') iRbc(1:Ntor+1,0)
 !write(iunit,'(" Zas         = ",99es23.15)') iZbs(1:Ntor+1,0) 
 !write(iunit,'(" Ras         = ",99es23.15)') iRbs(1:Ntor+1,0) 
 !write(iunit,'(" Zac         = ",99es23.15)') iZbc(1:Ntor+1,0) 

  do mm = 0, Mpol ! will write out the plasma boundary harmonics; 01 May 13;
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


  do mm = 0, Mpol ! will write out the computation domain harmonics; (only relevant in free-boundary case); 01 May 13;
   do nn = -Ntor, Ntor
    
    if( mm.eq.0 .and. nn.lt.0 ) cycle ! these modes are always excluded; 13 Oct 12;
    
    select case( mm )
    case(   0:  9 )
     if( nn.lt.- 9 .and. nn.gt.-99 ) &
  write(iunit,1010) nn, mm, Rwc(nn,mm), nn, mm, Zws(nn,mm), nn, mm, Rws(nn,mm), nn, mm, Zwc(nn,mm), nn, mm, Bns(nn,mm), nn, mm, Bnc(nn,mm)
     if( nn.lt.  0 .and. nn.gt.- 9 ) &
  write(iunit,1011) nn, mm, Rwc(nn,mm), nn, mm, Zws(nn,mm), nn, mm, Rws(nn,mm), nn, mm, Zwc(nn,mm), nn, mm, Bns(nn,mm), nn, mm, Bnc(nn,mm)
     if( nn.ge.  0 .and. nn.le.  9 ) &
  write(iunit,1012) nn, mm, Rwc(nn,mm), nn, mm, Zws(nn,mm), nn, mm, Rws(nn,mm), nn, mm, Zwc(nn,mm), nn, mm, Bns(nn,mm), nn, mm, Bnc(nn,mm)
     if( nn.gt.  9 .and. nn.le. 99 ) &
  write(iunit,1011) nn, mm, Rwc(nn,mm), nn, mm, Zws(nn,mm), nn, mm, Rws(nn,mm), nn, mm, Zwc(nn,mm), nn, mm, Bns(nn,mm), nn, mm, Bnc(nn,mm)
    case(  10: 99 )
     if( nn.lt.- 9 .and. nn.gt.-99 ) &
  write(iunit,1013) nn, mm, Rwc(nn,mm), nn, mm, Zws(nn,mm), nn, mm, Rws(nn,mm), nn, mm, Zwc(nn,mm), nn, mm, Bns(nn,mm), nn, mm, Bnc(nn,mm)
     if( nn.lt.  0 .and. nn.gt.- 9 ) &
  write(iunit,1014) nn, mm, Rwc(nn,mm), nn, mm, Zws(nn,mm), nn, mm, Rws(nn,mm), nn, mm, Zwc(nn,mm), nn, mm, Bns(nn,mm), nn, mm, Bnc(nn,mm)
     if( nn.ge.  0 .and. nn.le.  9 ) &
  write(iunit,1015) nn, mm, Rwc(nn,mm), nn, mm, Zws(nn,mm), nn, mm, Rws(nn,mm), nn, mm, Zwc(nn,mm), nn, mm, Bns(nn,mm), nn, mm, Bnc(nn,mm)
     if( nn.gt.  9 .and. nn.le. 99 ) &
  write(iunit,1014) nn, mm, Rwc(nn,mm), nn, mm, Zws(nn,mm), nn, mm, Rws(nn,mm), nn, mm, Zwc(nn,mm), nn, mm, Bns(nn,mm), nn, mm, Bnc(nn,mm)
    end select ! end of select case( mm );

   enddo ! end of do nn;
  enddo ! end of do mm;


1000 format("Rbc(",i3,",",i1,")",2x,"=",es23.15," Zbs(",i3,",",i1,")",2x,"=",es23.15," Rbs(",i3,",",i1,")",2x,"=",es23.15," Zbc(",i3,",",i1,")",2x,"=",es23.15)
1001 format("Rbc(",i2,",",i1,")",3x,"=",es23.15," Zbs(",i2,",",i1,")",3x,"=",es23.15," Rbs(",i2,",",i1,")",3x,"=",es23.15," Zbc(",i2,",",i1,")",3x,"=",es23.15)
1002 format("Rbc(",i1,",",i1,")",4x,"=",es23.15," Zbs(",i1,",",i1,")",4x,"=",es23.15," Rbs(",i1,",",i1,")",4x,"=",es23.15," Zbc(",i1,",",i1,")",4x,"=",es23.15)
1003 format("Rbc(",i3,",",i2,")",1x,"=",es23.15," Zbs(",i3,",",i2,")",1x,"=",es23.15," Rbs(",i3,",",i2,")",1x,"=",es23.15," Zbc(",i3,",",i2,")",1x,"=",es23.15)
1004 format("Rbc(",i2,",",i2,")",2x,"=",es23.15," Zbs(",i2,",",i2,")",2x,"=",es23.15," Rbs(",i2,",",i2,")",2x,"=",es23.15," Zbc(",i2,",",i2,")",2x,"=",es23.15)
1005 format("Rbc(",i1,",",i2,")",3x,"=",es23.15," Zbs(",i1,",",i2,")",3x,"=",es23.15," Rbs(",i1,",",i2,")",3x,"=",es23.15," Zbc(",i1,",",i2,")",3x,"=",es23.15)
  
1010 format("Rwc(",i3,",",i1,")",2x,"=",es23.15," Zws(",i3,",",i1,")",2x,"=",es23.15," Rws(",i3,",",i1,")",2x,"=",es23.15," Zwc(",i3,",",i1,")",2x,"=",es23.15&
  ,         " Bns(",i3,",",i1,")",2x,"=",es23.15," Bnc(",i3,",",i1,")",2x,"=",es23.15)
1011 format("Rwc(",i2,",",i1,")",3x,"=",es23.15," Zws(",i2,",",i1,")",3x,"=",es23.15," Rws(",i2,",",i1,")",3x,"=",es23.15," Zwc(",i2,",",i1,")",3x,"=",es23.15&
  ,         " Bns(",i2,",",i1,")",3x,"=",es23.15," Bnc(",i2,",",i1,")",3x,"=",es23.15)
1012 format("Rwc(",i1,",",i1,")",4x,"=",es23.15," Zws(",i1,",",i1,")",4x,"=",es23.15," Rws(",i1,",",i1,")",4x,"=",es23.15," Zwc(",i1,",",i1,")",4x,"=",es23.15&
  ,         " Bns(",i1,",",i1,")",4x,"=",es23.15," Bnc(",i1,",",i1,")",4x,"=",es23.15)
1013 format("Rwc(",i3,",",i2,")",1x,"=",es23.15," Zws(",i3,",",i2,")",1x,"=",es23.15," Rws(",i3,",",i2,")",1x,"=",es23.15," Zwc(",i3,",",i2,")",1x,"=",es23.15&
  ,         " Bns(",i3,",",i2,")",1x,"=",es23.15," Bnc(",i3,",",i2,")",1x,"=",es23.15)
1014 format("Rwc(",i2,",",i2,")",2x,"=",es23.15," Zws(",i2,",",i2,")",2x,"=",es23.15," Rws(",i2,",",i2,")",2x,"=",es23.15," Zwc(",i2,",",i2,")",2x,"=",es23.15&
  ,         " Bns(",i2,",",i2,")",2x,"=",es23.15," Bnc(",i2,",",i2,")",2x,"=",es23.15)
1015 format("Rwc(",i1,",",i2,")",3x,"=",es23.15," Zws(",i1,",",i2,")",3x,"=",es23.15," Rws(",i1,",",i2,")",3x,"=",es23.15," Zwc(",i1,",",i2,")",3x,"=",es23.15&
  ,         " Bns(",i1,",",i2,")",3x,"=",es23.15," Bnc(",i1,",",i2,")",3x,"=",es23.15)
  
  write(iunit,'("/")')

  if( Wwritin ) then ; cput = GETTIME ; write(ounit,'("writin : ",f10.2," : myid=",i3," ; writing numericlist ;")') cput-cpus, myid
  endif
  
  write(iunit,'("&numericlist")')
  write(iunit,'(" Linitialize = ",i9            )') Linitialize
  write(iunit,'(" Lwall       = ",i9            )') Lwall
  write(iunit,'(" phiwall     = ",es23.15       )') phiwall
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
  write(iunit,'(" Iswmin      = ",i9            )') Iswmin
  write(iunit,'(" Lperturb    = ",i9            )') Lperturb
  write(iunit,'(" dperturb    = ",es23.15       )') dperturb
  write(iunit,'(" Lextrap     = ",i9            )') Lextrap
  write(iunit,'("/")')

  if( Wwritin ) then ; cput = GETTIME ; write(ounit,'("writin : ",f10.2," : myid=",i3," ; writing locallist ;")') cput-cpus, myid
  endif
  
  write(iunit,'("&locallist")')
  write(iunit,'(" LBeltrami   = ",i9            )') LBeltrami
  write(iunit,'(" Linitgues   = ",i9            )') Linitgues
  write(iunit,'(" Lposdef     = ",L9            )') Lposdef
  write(iunit,'(" Nmaxexp     = ",i9            )') Nmaxexp
! write(iunit,'(" Liotasolv   = ",i9            )') Liotasolv
  write(iunit,'("/")')

  if( Wwritin ) then ; cput = GETTIME ; write(ounit,'("writin : ",f10.2," : myid=",i3," ; writing globallist ;")') cput-cpus, myid
  endif
  
  write(iunit,'("&globallist")')
  write(iunit,'(" Lminimize   = ",i9            )') Lminimize
  write(iunit,'(" Lfindzero   = ",i9            )') Lfindzero
 !write(iunit,'(" Lcondense   = ",i9            )') Lcondense ! 18 Jul 14;
  write(iunit,'(" escale      = ",es23.15       )') escale
  write(iunit,'(" pcondense   = ",es23.15       )') pcondense
 !write(iunit,'(" qcondense   = ",es23.15       )') qcondense ! 18 Jul 14; ! 04 Dec 14;
  write(iunit,'(" wpoloidal   = ",es23.15       )') wpoloidal ! minimum poloidal length constraint exponential factor; 04 Dec 14;
  write(iunit,'(" forcetol    = ",es23.15       )') forcetol
  write(iunit,'(" normalerr   = ",es23.15       )') normalerr
  write(iunit,'(" norblend    = ",es23.15       )') norblend
  write(iunit,'(" maxfbits    = ",i9            )') maxfbits
  write(iunit,'(" ForceErr    = ",es23.15       )') ForceErr ! 02 Nov 12;
 !write(iunit,'(" Energy      = ",es23.15       )') Energy ! 02 Nov 12;
  write(iunit,'(" c05xtol     = ",es23.15       )') c05xtol
  write(iunit,'(" c05factor   = ",es23.15       )') c05factor
  write(iunit,'(" LreadGF     = ",L9            )') LreadGF
  write(iunit,'(" verify      = ",i9            )') verify
  write(iunit,'(" maxstep     = ",es23.15       )') maxstep
  write(iunit,'(" opsilon     = ",es23.15       )') opsilon ! weight factor for pressure imbalance constraint; 04 Dec 14;
  write(iunit,'(" epsilon     = ",es23.15       )') epsilon ! weight factor for spectral condensation constraint; 04 Dec 14;
  write(iunit,'(" upsilon     = ",es23.15       )') upsilon ! weight factor for poloidal lenght constraint; 04 Dec 14;
 !write(iunit,'(" apsilon     = ",es23.15       )') apsilon ! 18 Jul 14;
  write(iunit,'(" maxiter     = ",i9            )') maxiter
  write(iunit,'("/")')

  if( Wwritin ) then ; cput = GETTIME ; write(ounit,'("writin : ",f10.2," : myid=",i3," ; writing diagnosticslist ;")') cput-cpus, myid
  endif
  
  write(iunit,'("&diagnosticslist")')
  write(iunit,'(" odetol      = ",es23.15       )') odetol
  write(iunit,'(" absreq      = ",es23.15       )') absreq
  write(iunit,'(" relreq      = ",es23.15       )') relreq
  write(iunit,'(" absacc      = ",es23.15       )') absacc
  write(iunit,'(" epsr        = ",es23.15       )') epsr
  write(iunit,'(" nPpts       = ",i9            )') nPpts
  write(iunit,'(" Mpqits      = ",i9            )') Mpqits
  write(iunit,'(" Lpqsym      = ",i9            )') Lpqsym
  write(iunit,'(" p1          = ",256i6         )') p1(1:Mvol)
  write(iunit,'(" q1          = ",256i6         )') q1(1:Mvol)
  write(iunit,'(" p2          = ",256i6         )') p2(1:Mvol)
  write(iunit,'(" q2          = ",256i6         )') q2(1:Mvol)
  write(iunit,'(" pqs         = ",256es23.15    )') pqs(1:Mvol)
  write(iunit,'(" pqt         = ",256es23.15    )') pqt(1:Mvol)
! write(iunit,'(" Npl         = ",i9            )') Npl
! write(iunit,'(" Nunstable   = ",i9            )') Nunstable
! write(iunit,'(" dunstable   = ",es23.15       )') dunstable
! write(iunit,'(" Munstable   = ",i9            )') Munstable
  write(iunit,'(" npq         = ",256i6         )') npq(1:Mvol)
  write(iunit,'(" nPtrj       = ",256i6         )') nPtrj(1:Mvol)
! write(iunit,'(" Mirrits     = ",i9            )') Mirrits
! write(iunit,'(" irrMpol     = ",i9            )') irrMpol
! write(iunit,'(" irrNtor     = ",i9            )') irrNtor
! write(iunit,'(" irrsvdcut   = ",es23.15       )') irrsvdcut
! write(iunit,'(" irrtol      = ",es23.15       )') irrtol
! write(iunit,'(" Lwrpj       = ",L9            )') Lwrpj
! write(iunit,'(" Nghd        = ",i9            )') Nghd
! write(iunit,'(" LHessian    = ",L9            )') LHessian
  write(iunit,'(" LHevalues   = ",L9            )') LHevalues
  write(iunit,'(" LHevectors  = ",L9            )') LHevectors
  write(iunit,'(" Lperturbed  = ",i9            )') Lperturbed
  write(iunit,'(" dpp         = ",i9            )') dpp
  write(iunit,'(" dqq         = ",i9            )') dqq
! write(iunit,'(" Lcurlerr    = ",L9            )') Lcurlerr
  write(iunit,'(" Lcheck      = ",i9            )') Lcheck
  write(iunit,'(" Ltiming     = ",L9            )') Ltiming
  write(iunit,'("/")')

  if( Wwritin ) then ; cput = GETTIME ; write(ounit,'("writin : ",f10.2," : myid=",i3," ; writing screenlist ;")') cput-cpus, myid
  endif
  
  write(iunit,'("&screenlist")')
! WSCREENLIST ! write screenlist; this is expanded by Makefile ; do not remove;
  if( Wreadin           ) write(iunit,'(" Wreadin = ",L1                )') Wreadin
  if( Wwritin           ) write(iunit,'(" Wwritin = ",L1                )') Wwritin
  if( Wmacros           ) write(iunit,'(" Wmacros = ",L1                )') Wmacros
  write(iunit,'("/")')
  
!latex \item Note that the boundary harmonics are changed by the spectral condensation routine \verb+sw00ab+, which is called by \verb+spec+.
!latex At finite resolution, as the boundary harmonics are changed, the boundary surface itself may change. 
!latex For restarts at potentially higher resolution, the original boundary harmonics (saved in \verb+Rbm+ and \verb+Zbm+, see \verb+global/readin+)
!latex are written to file.
  
#ifdef DEBUG
  FATALMESS(writin, .not.allocated(iRbc), error)
  FATALMESS(writin, .not.allocated(iZbs), error)
  FATALMESS(writin, .not.allocated(iRbs), error)
  FATALMESS(writin, .not.allocated(iZbc), error)
#endif

  do imn = 1, mn ; write(iunit,'(2i6,1024es23.15)') im(imn), in(imn)/Nfp, ( iRbc(imn,vvol), iZbs(imn,vvol), iRbs(imn,vvol), iZbc(imn,vvol), vvol = 1, Nvol )
  enddo
  
  close(iunit)
  
  if( Wwritin ) then ; cput = GETTIME ; write(ounit,'("writin : ",f10.2," : myid=",i3," ; wrote ext.end ;")') cput-cpus, myid
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \end{enumerate} \subsubsection{iterations towards convergence} \begin{enumerate}
!latex \item The file \verb+.ext.convergence+ is created.
!latex This file contains the evolution of the geometry towards convergence: i.e., the geometry of the interfaces at every iteration is saved.
!latex The format of this file is: needs revision;
!latex \begin{verbatim}
!latex open(zunit,file="."//trim(ext)//".convergence",status="unknown")!,form="unformatted")
!latex write(zunit,'(2i9)')mn, Nvol
!latex write(zunit,'(999i4)')im(1:mn)
!latex write(zunit,'(999i4)')in(1:mn)
!latex ! begin do loop over iterations;
!latex write(zunit)wflag,iflag,Energy,rflag
!latex do imn = 1, mn
!latex  write(zunit,'(999es23.15)') ( iRbc(imn,vvol), iZbs(imn,vvol), vvol = 0, Nvol )
!latex enddo
!latex !  end  do loop over iterations;
!latex \end{verbatim}
!latex \item The \verb+wflag+ integer indicates algorithm: \verb+wflag=1+ indicates \verb+E04DGF+; \verb+wflag=2+ indicates \verb+C05PDF+.
!latex \item The \verb+iflag+ integer indicates number of derivative evaluations required by \verb+C05PDF+. 

  if( wflag.gt.0 ) then ! this writes to convergence file;

   if( Wwritin ) then ; cput = GETTIME ; write(ounit,'("writin : ",f10.2," : myid=",i3," ; writing to zunit ;")') cput-cpus, myid
   endif

   write(zunit) wflag, iflag, Energy, rflag ! this file is opened in xspech; 11 Aug 14;

   if( Igeometry.eq.3 ) then

!   FATALMESS(global, .true., what is iRbc etc.)

!#ifdef OLDAXIS

!    iRbc(1:Ntor+1,0) = iRbc(1:Ntor+1,1) ; iRbc(Ntor+2:mn,0) = zero ! regularization extrapolation; 19 Sep 13;
!    iZbs(1:Ntor+1,0) = iZbs(1:Ntor+1,1) ; iZbs(Ntor+2:mn,0) = zero
!    iRbs(1:Ntor+1,0) = iRbs(1:Ntor+1,1) ; iRbs(Ntor+2:mn,0) = zero
!    iZbc(1:Ntor+1,0) = iZbc(1:Ntor+1,1) ; iZbc(Ntor+2:mn,0) = zero

!#else

    FATALMESS(writin, .not.allocated(bjk), error)

    do kk = 1, mn ! see also gf00aa; 11 Aug 14;
     iRbc(kk,0) = sum( iRbc(1:mn,1) * abs( bjk(1:mn,kk) ) )
     iZbs(kk,0) = sum( iZbs(1:mn,1) *      bjk(1:mn,kk)   )
     iRbs(kk,0) = sum( iRbs(1:mn,1) *      bjk(1:mn,kk)   )
     iZbc(kk,0) = sum( iZbc(1:mn,1) * abs( bjk(1:mn,kk) ) )
    enddo

!#endif

   endif ! end of if( Igeometry.eq.3 ) ; 04 Dec 14;
    
   write(zunit) iRbc(1:mn,0:Mvol)
   write(zunit) iZbs(1:mn,0:Mvol)
   write(zunit) iRbs(1:mn,0:Mvol)
   write(zunit) iZbc(1:mn,0:Mvol)

   call flush(zunit) ! this file is opened in xspech; 11 Aug 14;

  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Wwritin ) then ; cput = GETTIME ; write(ounit,'("writin : ",f10.2," : myid=",i3," ; writing .ext.grid ;")') cput-cpus, myid
  endif

  open(iunit, file="."//trim(ext)//".grid", status="unknown", form="unformatted" ) ! coordinate grid;
  
  write(iunit) Nt, Nz, Ntz, Mvol, Igeometry
  
  do vvol = 1, Mvol ! volumes;
   
   if( Igeometry.eq.1 .or. vvol.gt.1 ) then ; Lcoordinatesingularity = .false. ! 22 Apr 13;
   else                                     ; Lcoordinatesingularity = .true.  ! 22 Apr 13;
   endif
   
   if( vvol.le.Nvol ) then ; Lplasmaregion = .true.
   else                    ; Lplasmaregion = .false.
   endif                   ; Lvacuumregion = .not. Lplasmaregion
   
   write(iunit) Lrad(vvol) ! sub-grid radial resolution; not really sub-grid resolution, but really the Chebyshev resolution; 29 Jan 13; 18 Dec 14;
   
   do ii = 0, Lrad(vvol) ! sub-grid;
    
    lss = - one + ii * two / Lrad(vvol)
    
    Lcurvature = 0 ; WCALL( writin, co01aa, ( vvol, lss, Lcurvature, Ntz, mn ) ) ! only Rij(0,:) and Zij(0,:) are required; Rmn & Zmn are available;
    
#ifdef DEBUG
    FATALMESS(writin, .not.allocated(Rij), error)
    FATALMESS(writin, .not.allocated(Zij), error)
#endif

    write(iunit) Rij(1:Ntz,0,0)
    write(iunit) Zij(1:Ntz,0,0)
    
   enddo ! end of do ii; 29 Jan 13;
   
  enddo ! end of do vvol; 29 Jan 13;
  
  close(iunit)
  
  if( Wwritin ) then ; cput = GETTIME ; write(ounit,'("writin : ",f10.2," : myid=",i3," ; opened /wrote   .ext.grid ;")') cput-cpus, myid
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(writin)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine writin

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine option( intvariable, fltvariable ) ! this is an example of optional arguments; 03 Apr 13;
  
  use constants, only :
  use numerical, only :
  use fileunits, only : ounit
  use inputlist, only : ext
  use cputiming, only : 
  
  LOCALS
  
  INTEGER, intent(in), optional :: intvariable
  REAL   , intent(in), optional :: fltvariable

  cput = GETTIME

  if( present(intvariable) ) write(ounit,'("option : ",f10.2," : myid=",i3," ; integer = "i9" ;")') cput-cpus, myid, intvariable
  if( present(fltvariable) ) write(ounit,'("option : ",f10.2," : myid=",i3," ; real    = "es23.15" ;")') cput-cpus, myid, fltvariable

  return
  
end subroutine option
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end module allglobal

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

module pjhamilto

  INTEGER              :: io

  REAL,    allocatable :: Itor(:,:), Gpol(:,:), spmn(:,:,:)

end module pjhamilto

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

module laplaces

  LOGICAL              :: stage1, exterior, dorm
  INTEGER              :: Nintervals, Nsegments, IC, NP4, NP1
  INTEGER, allocatable :: icint(:)
  REAL                 :: originalalpha
  REAL, allocatable    :: xpoly(:), ypoly(:), phi(:), phid(:), CC(:,:)

  INTEGER              :: ilength
  REAL                 :: totallength

  INTEGER              :: niterations ! counter; eventually redundant; 24 Oct 12;
  
  INTEGER              :: iangle ! angle ! eventually redundant; 24 Oct 12;

  REAL                 :: Rmid, Zmid ! used to define local polar coordinate; eventually redundant; 24 Oct 12;

  REAL                 :: alpha ! eventually redundant; 24 Oct 12;

end module laplaces

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \end{enumerate} \subsection{Compilation instructions} \begin{enumerate}

!latex \item The source is kept under CVS: \verb+>cvs -d /u/shudson/cvs_Spec/ checkout Spec+

!latex \item Compilation is provided by a Makefile: \verb+>make xspec+. Try \verb+>make help+ for compilation options.

!latex \begin{enumerate}

!latex \item The compilation flags are given by \verb+FLAGS+. These may be over-ruled by command line arguments.

!latex \item Compilation flags must be set that convert single precision to double precision, e.g. \verb+make FLAGS="--dbl"+.

!latex \item The NAG library is used and must be correctly linked.

!latex \item To apply a patch: \verb+>patch -p1 < file+

!latex \item The free-boundary capability employs \verb+mgrid+. For this, \verb+stellinstall_sanssiesta.zip+ is included.
!latex Simply open a new directory, e.g. \verb+>mkdir Stellopt+, unzip, e.g. \verb+>unzip stellinstall_sanssiesta.zip+, and execute \verb+setup+.
!latex This will place the required files into \verb+$(HOME)/bin_847/+.

!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
