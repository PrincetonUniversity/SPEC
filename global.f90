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

! queries the numerical precision of the machine you are building SPEC on
module numerical

!latex \subsubsection{\type{machprec}, \type{vsmall}, \type{small}, \type{sqrtmachprec} : machine precision}

!latex \begin{enumerate}
!latex \item The machine precision was determined using the Fortran 90 intrinsic function EPSILON and the return value was put into the code statically.
!latex \end{enumerate}

  implicit none

  REAL, parameter :: machprec = 1.11e-16 ! 0.5*epsilon(one) for 64 bit double precision
  REAL, parameter :: vsmall = 100*machprec
  REAL, parameter :: small = 10000*machprec
  REAL, parameter :: sqrtmachprec = sqrt(machprec) ! these were previously assigned below in readin via a call to NAG routine;
  REAL, parameter :: logtolerance = 1.0e-32 ! this is used to avoid taking alog10(zero); see e.g. dforce;

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

   type derivative
     LOGICAL :: L
     INTEGER :: vol      ! Used in coords.f90; required for global constraint force gradient evaluation
     INTEGER :: innout
     INTEGER :: ii
     INTEGER :: irz
     INTEGER :: issym
  end type derivative

end module typedefns

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

  LOGICAL                    :: Localconstraint

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

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

contains

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

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

! Use this routine to run SPEC on a different MPI communicator.
! It takes care to re-assign the rank and size information
! that is used internally in SPEC.
subroutine set_mpi_comm(comm)

  implicit none

  integer, intent(in) :: comm
  integer             :: ierr

  ! MPI_COMM_SPEC is the global variable for the SPEC communicator.
  MPI_COMM_SPEC = comm

  myid = 0 ; ncpu = 1

  ierr=0
  call MPI_COMM_RANK( MPI_COMM_SPEC, myid, ierr )
  if (ierr.ne.0) write(*,*) "error in call to MPI_COMM_RANK"

  ierr=0
  call MPI_COMM_SIZE( MPI_COMM_SPEC, ncpu, ierr )
  if (ierr.ne.0) write(*,*) "error in call to MPI_COMM_SIZE"

end subroutine

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

  LOGICAL              :: Lchangeangle
  INTEGER              :: mm, nn, nb, imn, ix, ii, jj, ij, kk, mj, nj, mk, nk, ip, X02BBF, iargc, iarg, numargs, mi, ni, lvol
  REAL                 :: xx
  REAL,    allocatable :: RZRZ(:,:) ! local array used for reading interface Fourier harmonics from file;

  BEGIN(readin)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  cput = GETTIME

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( myid.eq.0 ) then ! only the master node reads input file and sets secondary variables;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   call read_inputlists_from_file

   call check_inputs

  endif ! end of if myid eq 0 loop;






!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
! broadcast command line input
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  call broadcast_inputs

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

      iRbc(ii,Mvol) = Rwc( nn, mm)                         ! computational boundary is ALWAYS given by namelist Rwc & Zws;
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

      iRbc(ii,Mvol) =   Rwc( kk, mm) + Rwc(-kk,-mm)        ! computational boundary is ALWAYS given by namelist Rwc & Zws;
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
