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




















































subroutine read_inputlists_from_file(num_modes, mmRZRZ, nnRZRZ, allRZRZ)

   use constants
   use fileunits
   use inputlist

   LOCALS

   integer, intent(out) :: num_modes
   integer, allocatable, intent(out) :: mmRZRZ(:), nnRZRZ(:)
   REAL,    allocatable, intent(out) :: allRZRZ(:,:,:) ! local array used for reading interface Fourier harmonics from file;

   LOGICAL              :: Lspexist
   integer :: filepos, seek_status, cpfile, instat, idx_mode

   INTEGER              :: mm, nn
   REAL,    allocatable :: RZRZ(:,:) ! local array used for reading interface Fourier harmonics from file;

   inquire( file=trim(ext)//".sp", exist=Lspexist ) ! check if file exists;
   FATAL( readin, .not.Lspexist, the input file does not exist ) ! if not, abort;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   open( iunit, file=trim(ext)//".sp", status="old")

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   instat = 0 ! initially, no error

! read namelists one after another
   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : reading physicslist     from ext.sp ;")') cput-cpus
   endif

   read(iunit,physicslist, iostat=instat)

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : read    physicslist     from ext.sp ;")') cput-cpus
   endif

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : reading numericlist     from ext.sp ;")') cput-cpus
   endif

   read(iunit,numericlist, iostat=instat)

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : read    numericlist     from ext.sp ;")') cput-cpus
   endif

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : reading locallist      from ext.sp ;")') cput-cpus
   endif

   read(iunit,locallist, iostat=instat)

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : read    locallist      from ext.sp ;")') cput-cpus
   endif

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : reading globallist   from ext.sp ;")') cput-cpus
   endif

   read(iunit,globallist, iostat=instat)

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : read    globallist   from ext.sp ;")') cput-cpus
   endif

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : reading diagnosticslist from ext.sp ;")') cput-cpus
   endif

   read(iunit,diagnosticslist, iostat=instat)

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : read    diagnosticslist from ext.sp ;")') cput-cpus
   endif

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : reading screenlist      from ext.sp ;")') cput-cpus
   endif

   read(iunit,screenlist, iostat=instat)

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : read    screenlist      from ext.sp ;")') cput-cpus
   endif

   ! At this point, the input namelists are read.
   ! It remains to read the initial guess for the interface geometry,
   ! which follows after the namelists in the input file.

   num_modes = 0
   if (Linitialize .le. 0) then

     SALLOCATE( RZRZ, (1:4,1:Nvol), zero ) ! temp array for reading input;

     ! determine how many modes are specified by reading them one
     call ftell(iunit, filepos)


     do ! will read in Fourier harmonics until the end of file is reached;
       read(iunit,*,iostat=instat) mm, nn, RZRZ(1:4,1:Nvol)   !if change of angle applies, transformation assumes m>=0 and for m=0 only n>=0;
       if( instat.ne.0 ) exit

       num_modes = num_modes + 1
     enddo ! end of do;

     ! rewind file to reset EOF flag
     ! and seek back to (start of modes) == (end of input namelists)
     rewind(iunit)
     call fseek(iunit, filepos, 0, seek_status)
     FATAL(inplst, seek_status.ne.0, failed to seek back to end of input namelists )

     ! now allocate arrays and read...
     allocate(mmRZRZ(1:num_modes), nnRZRZ(1:num_modes), allRZRZ(1:4,1:Nvol,1:num_modes))

     do idx_mode = 1, num_modes
       read(iunit,*,iostat=instat) mmRZRZ(idx_mode), nnRZRZ(idx_mode), allRZRZ(1:4,1:Nvol, idx_mode)
     enddo ! end of do;

     ! rewind file to reset EOF flag
     ! and seek back to (start of modes) == (end of input namelists)
     rewind(iunit)
     call fseek(iunit, filepos, 0, seek_status)
     FATAL(inplst, seek_status.ne.0, failed to seek back to end of input namelists )

     ! no need for temporary RZRZ anymore
     DALLOCATE(RZRZ)

    end if ! Linitialize .le. 0

    close(iunit)

end subroutine ! read_inputlists_from_file

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine check_inputs()

   use numerical
   use constants
   use fileunits
   use inputlist

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

  use fileunits
  use inputlist

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
