!> \defgroup grp_global Input namelists and global variables
!>
!> \latexonly
!> \definecolor{Orange}{rgb}{1.0,0.5,0.0}
!> \definecolor{Cerulean}{rgb}{0.0,0.5,1.0}
!> \endlatexonly
!>
!> \file
!> \brief Defines input namelists and global variables, and opens some output files.
!>
!> Note that all variables in namelist need to be broadcasted in readin.
!>
!> **Input geometry**
!> <ul>
!> <li> The geometry of the \f$l\f$-th interface, for \f$l=0,N\f$ where \f$N\equiv\f$ Nvol, is described by a set of Fourier harmonics,
!>      using an arbitrary poloidal angle,
!>      \f{eqnarray}{ R_l(\theta,\zeta)&=&\sum_{j}R_{j,l}\cos(m_j\theta-n_j\zeta), \\
!>                    Z_l(\theta,\zeta)&=&\sum_{j}Z_{j,l}\sin(m_j\theta-n_j\zeta). \f}
!> <li> These harmonics are read from the \c ext.sp file and come directly after the namelists described above.
!>      The required format is as follows:
!>      \f{eqnarray}{ \begin{array}{ccccccccc}
!>      m_1 & n_1 & R_{1,0} & Z_{1,0} & R_{1,1} & Z_{1,1} & ... & R_{1,N} & Z_{1,N} \\
!>      m_2 & n_2 & R_{2,0} & Z_{2,0} & R_{2,1} & Z_{2,1} & ... & R_{2,N} & Z_{2,N} \\
!>      ... \\
!>      m_j & n_j & R_{j,0} & Z_{j,0} & R_{j,1} & Z_{j,1} & ... & R_{j,N} & Z_{j,N} \\
!>      ...
!>      \end{array}
!>      \f}
!> <li> The coordinate axis corresponds to \f$j=0\f$ and the outermost boundary corresponds to \f$j=\f$ Nvol.
!> <li> An arbitrary selection of harmonics may be inluded in any order, but only those within the range specified by Mpol and Ntor will be used.
!> <li> The geometry of *all* the interfaces, i.e. \f$l=0,N\f$, including the degenerate "coordinate-axis" interface, must be given.
!> </ul>

!> \ingroup grp_global
!> \brief some constants used throughout the code
module constants
    use mod_kinds, only: wp => dp
    implicit none

    real(wp), parameter :: zero = 0.0 !< 0
    real(wp), parameter :: one = 1.0 !< 1
    real(wp), parameter :: two = 2.0 !< 2
    real(wp), parameter :: three = 3.0 !< 3
    real(wp), parameter :: four = 4.0 !< 4
    real(wp), parameter :: five = 5.0 !< 5
    real(wp), parameter :: six = 6.0 !< 6
    real(wp), parameter :: seven = 7.0 !< 7
    real(wp), parameter :: eight = 8.0 !< 8
    real(wp), parameter :: nine = 9.0 !< 9
    real(wp), parameter :: ten = 10.0 !< 10

    real(wp), parameter :: eleven = 11.0 !< 11
    real(wp), parameter :: twelve = 12.0 !< 12

    real(wp), parameter :: hundred = 100.0 !< 100
    real(wp), parameter :: thousand = 1000.0 !< 1000

    real(wp), parameter :: half = one/two   !< 1/2
    real(wp), parameter :: third = one/three !< 1/3
    real(wp), parameter :: quart = one/four  !< 1/4
    real(wp), parameter :: fifth = one/five  !< 1/5
    real(wp), parameter :: sixth = one/six   !< 1/6

    real(wp), parameter :: pi2 = 6.28318530717958623 !< \f$2\pi\f$
    real(wp), parameter :: pi = pi2/two           !< \f$\pi\f$
    real(wp), parameter :: mu0 = 2.0E-07*pi2       !< \f$4\pi\cdot10^{-7}\f$
    real(wp), parameter :: goldenmean = 1.618033988749895   !< golden mean = \f$( 1 + \sqrt 5 ) / 2\f$ ;

    real(wp), parameter :: version = 3.20  !< version of SPEC

end module constants

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> \brief platform-dependant numerical resolution
!> \ingroup grp_global
module numerical
    use mod_kinds, only: wp => dp
    implicit none

    real(wp), parameter :: machprec = 1.11e-16           !< machine precision: 0.5*epsilon(one) for 64 bit double precision
    real(wp), parameter :: vsmall = 100*machprec         !< very small number
    real(wp), parameter :: small = 10000*machprec        !< small number
    real(wp), parameter :: sqrtmachprec = sqrt(machprec) !< square root of machine precision
    real(wp), parameter :: logtolerance = 1.0e-32        !< this is used to avoid taking alog10(zero); see e.g. dforce;

end module numerical

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> \brief central definition of file units to avoid conflicts
!> \ingroup grp_global
module fileunits
    use mod_kinds, only: wp => dp
    implicit none

    integer :: iunit = 10 !< input; used in global/readin:ext.sp, global/wrtend:ext.sp.end
    integer :: ounit = 6 !< screen output;
    integer :: gunit = 13 !< wall geometry; used in wa00aa

    integer :: aunit = 11 !< vector potential; used in ra00aa:.ext.AtAzmn;
    integer :: dunit = 12 !< derivative matrix; used in newton:.ext.GF;
    integer :: hunit = 14 !< eigenvalues of Hessian; under re-construction;
    integer :: munit = 14 !< matrix elements of Hessian;
    integer :: lunit = 20 !< local unit; used in lunit+myid: pp00aa:.ext.poincare,.ext.transform;
    integer :: vunit = 15 !< for examination of adaptive quadrature; used in casing:.ext.vcint;

contains
    subroutine mute(action)
        implicit none

        integer, intent(in) :: action
        integer, parameter :: iopen = 1, iclose = 0, null = 37
        integer :: ios
        character(len=*), parameter :: nullfile = "/dev/null"

        ! open a tmp file for screen output
        if (action == iopen) then
            ounit = null
            open (ounit, file=nullfile, status="unknown", action="write", iostat=ios) ! create a scratch file?
            if (ios .ne. 0) print *, "something wrong with open a tmp file in focuspy.mute. IOSTAT=", ios
        else
            close (ounit)
            ounit = 6 ! recover to screen output
        end if
        return
    end subroutine mute

end module fileunits

!> \brief timing variables
!> \ingroup grp_global
module cputiming
    use mod_kinds, only: wp => dp
    real(wp) :: Tmanual = 0.0, manualT = 0.0
    real(wp) :: Trzaxis = 0.0, rzaxisT = 0.0
    real(wp) :: Tpackxi = 0.0, packxiT = 0.0
    real(wp) :: Tvolume = 0.0, volumeT = 0.0
    real(wp) :: Tcoords = 0.0, coordsT = 0.0
    real(wp) :: Tbasefn = 0.0, basefnT = 0.0
    real(wp) :: Tmemory = 0.0, memoryT = 0.0
    real(wp) :: Tmetrix = 0.0, metrixT = 0.0
    real(wp) :: Tma00aa = 0.0, ma00aaT = 0.0
    real(wp) :: Tmatrix = 0.0, matrixT = 0.0
    real(wp) :: Tspsmat = 0.0, spsmatT = 0.0
    real(wp) :: Tspsint = 0.0, spsintT = 0.0
    real(wp) :: Tmp00ac = 0.0, mp00acT = 0.0
    real(wp) :: Tma02aa = 0.0, ma02aaT = 0.0
    real(wp) :: Tpackab = 0.0, packabT = 0.0
    real(wp) :: Ttr00ab = 0.0, tr00abT = 0.0
    real(wp) :: Tcurent = 0.0, curentT = 0.0
    real(wp) :: Tdf00ab = 0.0, df00abT = 0.0
    real(wp) :: Tlforce = 0.0, lforceT = 0.0
    real(wp) :: Tintghs = 0.0, intghsT = 0.0
    real(wp) :: Tmtrxhs = 0.0, mtrxhsT = 0.0
    real(wp) :: Tlbpol = 0.0, lbpolT = 0.0
    real(wp) :: Tbrcast = 0.0, brcastT = 0.0
    real(wp) :: Tdfp100 = 0.0, dfp100T = 0.0
    real(wp) :: Tdfp200 = 0.0, dfp200T = 0.0
    real(wp) :: Tdforce = 0.0, dforceT = 0.0
    real(wp) :: Tnewton = 0.0, newtonT = 0.0
    real(wp) :: Tcasing = 0.0, casingT = 0.0
    real(wp) :: Tbnorml = 0.0, bnormlT = 0.0
    real(wp) :: Tjo00aa = 0.0, jo00aaT = 0.0
    real(wp) :: Tpp00aa = 0.0, pp00aaT = 0.0
    real(wp) :: Tpp00ab = 0.0, pp00abT = 0.0
    real(wp) :: Tbfield = 0.0, bfieldT = 0.0
    real(wp) :: Tstzxyz = 0.0, stzxyzT = 0.0
    real(wp) :: Thesian = 0.0, hesianT = 0.0
    real(wp) :: Tra00aa = 0.0, ra00aaT = 0.0
    real(wp) :: Tnumrec = 0.0, numrecT = 0.0
    real(wp) :: Tdcuhre = 0.0, dcuhreT = 0.0
    real(wp) :: Tminpack = 0.0, minpackT = 0.0
    real(wp) :: Tiqpack = 0.0, iqpackT = 0.0
    real(wp) :: Trksuite = 0.0, rksuiteT = 0.0
    real(wp) :: Ti1mach = 0.0, i1machT = 0.0
    real(wp) :: Td1mach = 0.0, d1machT = 0.0
    real(wp) :: Tilut = 0.0, ilutT = 0.0
    real(wp) :: Titers = 0.0, itersT = 0.0
    real(wp) :: Tsphdf5 = 0.0, sphdf5T = 0.0
    real(wp) :: Tpreset = 0.0, presetT = 0.0
    real(wp) :: Tglobal = 0.0, globalT = 0.0
    real(wp) :: Txspech = 0.0, xspechT = 0.0
    real(wp) :: Tinputlist = 0.0, inputlistT = 0.0
    real(wp) :: Twa00aa = 0.0, wa00aaT = 0.0
    real(wp) :: Tpc00ab = 0.0, pc00abT = 0.0

    real(wp) :: Treadin = 0.0
!  REAL :: Twritin = 0.0 ! redundant;
    real(wp) :: Twrtend = 0.0

end module cputiming

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> \brief type definitions for custom datatypes
!> \ingroup grp_global
module typedefns
    use mod_kinds, only: wp => dp
    !> \brief used for quantities which have different resolutions in different volumes, e.g. the vector potential
    type subgrid
        real(wp), allocatable :: s(:) !< coefficients
        integer, allocatable :: i(:) !< indices
    end type subgrid

    type MatrixLU
        real(wp), allocatable :: mat(:, :)
        integer, allocatable :: ipivot(:)
    end type MatrixLU

    !> \brief \f${\rm d}\mathbf{B}/{\rm d}\mathbf{X}\f$ (?)
    type derivative
        LOGICAL :: L      !< what is this?
        integer :: vol    !< Used in coords(); required for global constraint force gradient evaluation
        integer :: innout !< what is this?
        integer :: ii     !< what is this?
        integer :: irz    !< what is this?
        integer :: issym  !< what is this?
    end type derivative

end module typedefns

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> \brief global variable storage used as "workspace" throughout the code
module allglobal
    use mod_kinds, only: wp => dp
    use constants
    use typedefns

    implicit none

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!``-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    integer :: myid !< MPI rank of current CPU
    integer :: ncpu !< number of MPI tasks
    integer :: IsMyVolumeValue !< flag to indicate if a CPU is operating on its assigned volume
    real(wp) :: cpus !< initial time
    integer :: MPI_COMM_SPEC !< SPEC MPI communicator

    LOGICAL :: skip_write = .false. ! flag to disable any HDF5-related calls

    real(wp) :: pi2nfp           !       pi2/nfp     ; assigned in readin;
    real(wp) :: pi2pi2nfp
    real(wp) :: pi2pi2nfphalf
    real(wp) :: pi2pi2nfpquart

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    character(len=255) :: ext ! extension of input filename, i.e., "G3V01L1Fi.001" for an input file G3V01L1Fi.001.sp

    real(wp) :: ForceErr !< total force-imbalance
    real(wp) :: Energy   !< MHD energy

    real(wp), allocatable :: IPDt(:), IPDtDpf(:, :)  !< Toroidal pressure-driven current

    integer :: Mvol

    LOGICAL :: YESstellsym !< internal shorthand copies of Istellsym, which is an integer input;
    LOGICAL :: NOTstellsym !< internal shorthand copies of Istellsym, which is an integer input;

    LOGICAL :: YESMatrixFree, NOTMatrixFree !< to use matrix-free method or not

    real(wp), allocatable :: cheby(:, :) !< local workspace for evaluation of Chebychev polynomials
    real(wp), allocatable :: zernike(:, :, :) !< local workspace for evaluation of Zernike polynomials

    real(wp), allocatable :: TT(:, :, :)    !< derivatives of Chebyshev polynomials at the inner and outer interfaces;
    real(wp), allocatable :: RTT(:, :, :, :) !< derivatives of Zernike   polynomials at the inner and outer interfaces;

    real(wp), allocatable :: RTM(:, :) !< \f$r^m\f$ term of Zernike polynomials at the origin
    real(wp), allocatable :: ZernikeDof(:) !< Zernike degree of freedom for each \f$m\f$

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> \addtogroup grp_enh_res_metr Enhanced resolution for metric elements
!> Enhanced resolution is required for the metric elements, \f$g_{ij}/\sqrt g\f$, which is given by mne, ime, and ine.
!> The Fourier resolution here is determined by \c lMpol=2*Mpol  and \c lNtor=2*Ntor.
!> @{
    integer :: mne    !< enhanced resolution for metric elements
    integer, allocatable :: ime(:) !< enhanced poloidal mode numbers for metric elements
    integer, allocatable :: ine(:) !< enhanced toroidal mode numbers for metric elements
!> @}

!> \addtogroup grp_enh_res_sfl Enhanced resolution for transformation to straight-field line angle
!> Enhanced resolution is required for the transformation to straight-field line angle on the interfaces,
!> which is given by mns, ims  and ins.
!> The Fourier resolution here is determined by \c iMpol  and \c iNtor.
!> @{
    integer :: mns    !< enhanced resolution for straight field line transformation
    integer, allocatable :: ims(:) !< enhanced poloidal mode numbers for straight field line transformation
    integer, allocatable :: ins(:) !< enhanced toroidal mode numbers for straight field line transformation
!> @}

    integer :: lMpol !< what is this?
    integer :: lNtor !< what is this?
    integer :: sMpol !< what is this?
    integer :: sNtor !< what is this?

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    real(wp) :: xoffset = 1.0 !< used to normalize NAG routines (which ones exacly where?)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    LOGICAL, allocatable :: ImagneticOK(:) !< used to indicate if Beltrami fields have been correctly constructed;

    LOGICAL :: IconstraintOK !< Used to break iteration loops of slaves in the global constraint minimization.

    real(wp), allocatable :: beltramierror(:, :)  !< to store the integral of |curlB-mu*B| computed by jo00aa;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> \addtogroup grp_internal_vars Internal Variables
!> @{
!>
!> \addtogroup grp_fourier_repr Fourier representation
!> @{
    integer :: mn    !< total number of Fourier harmonics for coordinates/fields; calculated from Mpol, Ntor in readin()
    integer, allocatable :: im(:) !< poloidal mode numbers for Fourier representation
    integer, allocatable :: in(:) !< toroidal mode numbers for Fourier representation

    real(wp), allocatable :: halfmm(:) !< I saw this already somewhere...
    real(wp), allocatable :: regumm(:) !< I saw this already somewhere...

    real(wp) :: Rscale    !< no idea
    real(wp), allocatable :: psifactor(:, :) !< no idea
    real(wp), allocatable :: inifactor(:, :) !< no idea

    real(wp), allocatable :: BBweight(:) !< weight on force-imbalance harmonics; used in dforce()

    real(wp), allocatable :: mmpp(:) !< spectral condensation factors

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!> @}

!> \addtogroup grp_iface_geom Interface geometry: iRbc, iZbs etc.
!> The Fourier harmonics of the interfaces are contained in \c iRbc(1:mn,0:Mvol) and \c iZbs(1:mn,0:Mvol), where
!> \c iRbc(l,j), \c iZbs(l,j) contains the Fourier harmonics, \f$R_j\f$, \f$Z_j\f$, of the \f$l\f$-th interface.
!> @{
    real(wp), allocatable :: iRbc(:, :) !< cosine R harmonics of interface surface geometry;     stellarator symmetric
    real(wp), allocatable :: iZbs(:, :) !<   sine Z harmonics of interface surface geometry;     stellarator symmetric
    real(wp), allocatable :: iRbs(:, :) !<   sine R harmonics of interface surface geometry; non-stellarator symmetric
    real(wp), allocatable :: iZbc(:, :) !< cosine Z harmonics of interface surface geometry; non-stellarator symmetric

    real(wp), allocatable :: dRbc(:, :) !< cosine R harmonics of interface surface geometry;     stellarator symmetric; linear deformation
    real(wp), allocatable :: dZbs(:, :) !<   sine Z harmonics of interface surface geometry;     stellarator symmetric; linear deformation
    real(wp), allocatable :: dRbs(:, :) !<   sine R harmonics of interface surface geometry; non-stellarator symmetric; linear deformation
    real(wp), allocatable :: dZbc(:, :) !< cosine Z harmonics of interface surface geometry; non-stellarator symmetric; linear deformation

    real(wp), allocatable :: iRij(:, :) !< interface surface geometry; real space
    real(wp), allocatable :: iZij(:, :) !< interface surface geometry; real space
    real(wp), allocatable :: dRij(:, :) !< interface surface geometry; real space
    real(wp), allocatable :: dZij(:, :) !< interface surface geometry; real space
    real(wp), allocatable :: tRij(:, :) !< interface surface geometry; real space
    real(wp), allocatable :: tZij(:, :) !< interface surface geometry; real space

    real(wp), allocatable :: iVns(:)   !<   sine harmonics of vacuum normal magnetic field on interfaces;     stellarator symmetric
    real(wp), allocatable :: iBns(:)   !<   sine harmonics of plasma normal magnetic field on interfaces;     stellarator symmetric
    real(wp), allocatable :: iVnc(:)   !< cosine harmonics of vacuum normal magnetic field on interfaces; non-stellarator symmetric
    real(wp), allocatable :: iBnc(:)   !< cosine harmonics of plasma normal magnetic field on interfaces; non-stellarator symmetric

    real(wp), allocatable :: lRbc(:)   !< local workspace
    real(wp), allocatable :: lZbs(:)   !< local workspace
    real(wp), allocatable :: lRbs(:)   !< local workspace
    real(wp), allocatable :: lZbc(:)   !< local workspace

    ! local array used for reading interface Fourier harmonics from file;
    integer :: num_modes
    integer, allocatable :: mmRZRZ(:), nnRZRZ(:)
    real(wp), allocatable :: allRZRZ(:, :, :)
!> @}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> \addtogroup grp_fourier_transform Fourier Transforms
!> The coordinate geometry and fields are mapped to/from Fourier space and real space using FFTW3.
!> The resolution of the real space grid is given by \c Nt=Ndiscrete*4*Mpol and \c Nz=Ndiscrete*4*Ntor.
!>
!> Various workspace arrays are allocated.
!> These include \c Rij(1:Ntz,0:3,0:3) and \c Zij(1:Ntz,0:3,0:3), which contain the coordinates in real space and their derivatives;
!> \c sg(0:3,Ntz), which contains the Jacobian and its derivatives;
!> and \c guv(0:6,0:3,1:Ntz), which contains the metric elements and their derivatives.
!> @{
    integer :: Nt  !< discrete resolution along \f$\theta\f$ of grid in real space
    integer :: Nz  !< discrete resolution along \f$\zeta\f$  of grid in real space
    integer :: Ntz !< discrete resolution; Ntz=Nt*Nz shorthand
    integer :: hNt !< discrete resolution; Ntz=Nt*Nz shorthand
    integer :: hNz !< discrete resolution; Ntz=Nt*Nz shorthand
    real(wp) :: soNtz !< one / sqrt (one*Ntz); shorthand

    real(wp), allocatable :: Rij(:, :, :) !< real-space grid; R
    real(wp), allocatable :: Zij(:, :, :) !< real-space grid; Z
    real(wp), allocatable :: Xij(:, :, :) !< what is this?
    real(wp), allocatable :: Yij(:, :, :) !< what is this?
    real(wp), allocatable :: sg(:, :)    !< real-space grid; jacobian and its derivatives
    real(wp), allocatable :: guvij(:, :, :, :) !< real-space grid; metric elements
    real(wp), allocatable :: gvuij(:, :, :)   !< real-space grid; metric elements (?); 10 Dec 15;
    real(wp), allocatable :: guvijsave(:, :, :, :) !< what is this?

    integer, allocatable :: ki(:, :)     !< identification of Fourier modes
    integer, allocatable :: kijs(:, :, :) !< identification of Fourier modes
    integer, allocatable :: kija(:, :, :) !< identification of Fourier modes

    integer, allocatable :: iotakkii(:)   !< identification of Fourier modes
    integer, allocatable :: iotaksub(:, :) !< identification of Fourier modes
    integer, allocatable :: iotakadd(:, :) !< identification of Fourier modes
    integer, allocatable :: iotaksgn(:, :) !< identification of Fourier modes

    real(wp), allocatable :: efmn(:) !< Fourier harmonics; dummy workspace
    real(wp), allocatable :: ofmn(:) !< Fourier harmonics; dummy workspace
    real(wp), allocatable :: cfmn(:) !< Fourier harmonics; dummy workspace
    real(wp), allocatable :: sfmn(:) !< Fourier harmonics; dummy workspace
    real(wp), allocatable :: evmn(:) !< Fourier harmonics; dummy workspace
    real(wp), allocatable :: odmn(:) !< Fourier harmonics; dummy workspace
    real(wp), allocatable :: comn(:) !< Fourier harmonics; dummy workspace
    real(wp), allocatable :: simn(:) !< Fourier harmonics; dummy workspace

    real(wp), allocatable :: ijreal(:) !< what is this ?
    real(wp), allocatable :: ijimag(:) !< what is this ?
    real(wp), allocatable :: jireal(:) !< what is this ?
    real(wp), allocatable :: jiimag(:) !< what is this ?

    real(wp), allocatable :: jkreal(:) !< what is this ?
    real(wp), allocatable :: jkimag(:) !< what is this ?
    real(wp), allocatable :: kjreal(:) !< what is this ?
    real(wp), allocatable :: kjimag(:) !< what is this ?

    real(wp), allocatable :: Bsupumn(:, :, :) !< tangential field on interfaces; \f$\theta\f$-component; required for virtual casing construction of field; 11 Oct 12
    real(wp), allocatable :: Bsupvmn(:, :, :) !< tangential field on interfaces; \f$\zeta\f$ -component; required for virtual casing construction of field; 11 Oct 12

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    real(wp), allocatable :: goomne(:, :) !< described in preset()
    real(wp), allocatable :: goomno(:, :) !< described in preset()
    real(wp), allocatable :: gssmne(:, :) !< described in preset()
    real(wp), allocatable :: gssmno(:, :) !< described in preset()
    real(wp), allocatable :: gstmne(:, :) !< described in preset()
    real(wp), allocatable :: gstmno(:, :) !< described in preset()
    real(wp), allocatable :: gszmne(:, :) !< described in preset()
    real(wp), allocatable :: gszmno(:, :) !< described in preset()
    real(wp), allocatable :: gttmne(:, :) !< described in preset()
    real(wp), allocatable :: gttmno(:, :) !< described in preset()
    real(wp), allocatable :: gtzmne(:, :) !< described in preset()
    real(wp), allocatable :: gtzmno(:, :) !< described in preset()
    real(wp), allocatable :: gzzmne(:, :) !< described in preset()
    real(wp), allocatable :: gzzmno(:, :) !< described in preset()
!> @}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> \addtogroup grp_chebychev_metric Volume-integrated Chebyshev-metrics
!> These are allocated in dforce(), defined in ma00aa(), and are used in matrix() to construct the matrices.
!> @{
    real(wp), allocatable :: DToocc(:, :, :, :) !< volume-integrated Chebychev-metrics; see matrix()
    real(wp), allocatable :: DToocs(:, :, :, :) !< volume-integrated Chebychev-metrics; see matrix()
    real(wp), allocatable :: DToosc(:, :, :, :) !< volume-integrated Chebychev-metrics; see matrix()
    real(wp), allocatable :: DTooss(:, :, :, :) !< volume-integrated Chebychev-metrics; see matrix()
    real(wp), allocatable :: TTsscc(:, :, :, :) !< volume-integrated Chebychev-metrics; see matrix()
    real(wp), allocatable :: TTsscs(:, :, :, :) !< volume-integrated Chebychev-metrics; see matrix()
    real(wp), allocatable :: TTsssc(:, :, :, :) !< volume-integrated Chebychev-metrics; see matrix()
    real(wp), allocatable :: TTssss(:, :, :, :) !< volume-integrated Chebychev-metrics; see matrix()
    real(wp), allocatable :: TDstcc(:, :, :, :) !< volume-integrated Chebychev-metrics; see matrix()
    real(wp), allocatable :: TDstcs(:, :, :, :) !< volume-integrated Chebychev-metrics; see matrix()
    real(wp), allocatable :: TDstsc(:, :, :, :) !< volume-integrated Chebychev-metrics; see matrix()
    real(wp), allocatable :: TDstss(:, :, :, :) !< volume-integrated Chebychev-metrics; see matrix()
    real(wp), allocatable :: TDszcc(:, :, :, :) !< volume-integrated Chebychev-metrics; see matrix()
    real(wp), allocatable :: TDszcs(:, :, :, :) !< volume-integrated Chebychev-metrics; see matrix()
    real(wp), allocatable :: TDszsc(:, :, :, :) !< volume-integrated Chebychev-metrics; see matrix()
    real(wp), allocatable :: TDszss(:, :, :, :) !< volume-integrated Chebychev-metrics; see matrix()
    real(wp), allocatable :: DDttcc(:, :, :, :) !< volume-integrated Chebychev-metrics; see matrix()
    real(wp), allocatable :: DDttcs(:, :, :, :) !< volume-integrated Chebychev-metrics; see matrix()
    real(wp), allocatable :: DDttsc(:, :, :, :) !< volume-integrated Chebychev-metrics; see matrix()
    real(wp), allocatable :: DDttss(:, :, :, :) !< volume-integrated Chebychev-metrics; see matrix()
    real(wp), allocatable :: DDtzcc(:, :, :, :) !< volume-integrated Chebychev-metrics; see matrix()
    real(wp), allocatable :: DDtzcs(:, :, :, :) !< volume-integrated Chebychev-metrics; see matrix()
    real(wp), allocatable :: DDtzsc(:, :, :, :) !< volume-integrated Chebychev-metrics; see matrix()
    real(wp), allocatable :: DDtzss(:, :, :, :) !< volume-integrated Chebychev-metrics; see matrix()
    real(wp), allocatable :: DDzzcc(:, :, :, :) !< volume-integrated Chebychev-metrics; see matrix()
    real(wp), allocatable :: DDzzcs(:, :, :, :) !< volume-integrated Chebychev-metrics; see matrix()
    real(wp), allocatable :: DDzzsc(:, :, :, :) !< volume-integrated Chebychev-metrics; see matrix()
    real(wp), allocatable :: DDzzss(:, :, :, :) !< volume-integrated Chebychev-metrics; see matrix()

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    real(wp), allocatable :: Tsc(:, :) !< what is this?
    real(wp), allocatable :: Tss(:, :) !< what is this?
    real(wp), allocatable :: Dtc(:, :) !< what is this?
    real(wp), allocatable :: Dts(:, :) !< what is this?
    real(wp), allocatable :: Dzc(:, :) !< what is this?
    real(wp), allocatable :: Dzs(:, :) !< what is this?
    real(wp), allocatable :: Ttc(:, :) !< what is this?
    real(wp), allocatable :: Tzc(:, :) !< what is this?
    real(wp), allocatable :: Tts(:, :) !< what is this?
    real(wp), allocatable :: Tzs(:, :) !< what is this?

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    real(wp), allocatable :: dtflux(:) !< \f$\delta \psi_{toroidal}\f$ in each annulus
    real(wp), allocatable :: dpflux(:) !< \f$\delta \psi_{poloidal}\f$ in each annulus

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    real(wp), allocatable :: sweight(:) !< minimum poloidal length constraint weight
!> @}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> \addtogroup grp_vecpot Vector potential and the Beltrami linear system
!> <ul>
!> <li> In each volume, the total degrees of freedom in the Beltrami linear system is \c NAdof(1:Nvol).
!>      This depends on \c Mpol, \c Ntor and \c Lrad(vvol). </li>
!> <li> The covariant components of the vector potential are written as
!>       \f{eqnarray}{ A_\theta & = & \sum_i \sum_{l=0}^L {\color{red}  A_{\theta,e,i,l}} \; T_{l}(s) \cos\alpha_i + \sum_i \sum_{l=0}^L {\color{Orange}  A_{\theta,o,i,l}} \; T_{l}(s) \sin\alpha_i \\
!>                     A_\zeta  & = & \sum_i \sum_{l=0}^L {\color{blue} A_{\zeta, e,i,l}} \; T_{l}(s) \cos\alpha_i + \sum_i \sum_{l=0}^L {\color{Cerulean}A_{\zeta ,o,i,l}} \; T_{l}(s) \sin\alpha_i ,
!>       \f}
!>       where \f$T_l(s)\f$ are the Chebyshev polynomials and \f$\alpha_i \equiv m_i \theta - n_i \zeta\f$. </li>
!> <li> The following internal arrays are declared in preset() :
!>
!>       \c dAte(0,i)%%s(l) \f$\equiv {\color{red}     A_{\theta,e,i,l}}\f$
!>
!>       \c dAze(0,i)%%s(l) \f$\equiv {\color{blue}    A_{\zeta, e,i,l}}\f$
!>
!>       \c dAto(0,i)%%s(l) \f$\equiv {\color{Orange}  A_{\theta,o,i,l}}\f$
!>
!>       \c dAzo(0,i)%%s(l) \f$\equiv {\color{Cerulean}A_{\zeta ,o,i,l}}\f$ </li>
!> </ul>
!> @{
    integer, allocatable :: NAdof(:) !< degrees of freedom in Beltrami fields in each annulus
    integer, allocatable :: Nfielddof(:) !< degrees of freedom in Beltrami fields in each annulus, field only, no Lagrange multipliers

    type(subgrid), allocatable :: Ate(:, :, :) !< magnetic vector potential cosine Fourier harmonics;     stellarator-symmetric
    type(subgrid), allocatable :: Aze(:, :, :) !< magnetic vector potential cosine Fourier harmonics;     stellarator-symmetric
    type(subgrid), allocatable :: Ato(:, :, :) !< magnetic vector potential   sine Fourier harmonics; non-stellarator-symmetric
    type(subgrid), allocatable :: Azo(:, :, :) !< magnetic vector potential   sine Fourier harmonics; non-stellarator-symmetric

    integer, allocatable :: Lma(:, :) !< Lagrange multipliers (?)
    integer, allocatable :: Lmb(:, :) !< Lagrange multipliers (?)
    integer, allocatable :: Lmc(:, :) !< Lagrange multipliers (?)
    integer, allocatable :: Lmd(:, :) !< Lagrange multipliers (?)
    integer, allocatable :: Lme(:, :) !< Lagrange multipliers (?)
    integer, allocatable :: Lmf(:, :) !< Lagrange multipliers (?)
    integer, allocatable :: Lmg(:, :) !< Lagrange multipliers (?)
    integer, allocatable :: Lmh(:, :) !< Lagrange multipliers (?)

    real(wp), allocatable :: Lmavalue(:, :) !< what is this?
    real(wp), allocatable :: Lmbvalue(:, :) !< what is this?
    real(wp), allocatable :: Lmcvalue(:, :) !< what is this?
    real(wp), allocatable :: Lmdvalue(:, :) !< what is this?
    real(wp), allocatable :: Lmevalue(:, :) !< what is this?
    real(wp), allocatable :: Lmfvalue(:, :) !< what is this?
    real(wp), allocatable :: Lmgvalue(:, :) !< what is this?
    real(wp), allocatable :: Lmhvalue(:, :) !< what is this?

    integer, allocatable :: Fso(:, :) !< what is this?
    integer, allocatable :: Fse(:, :) !< what is this?

    LOGICAL :: Lcoordinatesingularity !< set by \c LREGION macro; true if inside the innermost volume
    LOGICAL :: Lplasmaregion          !< set by \c LREGION macro; true if inside the plasma region
    LOGICAL :: Lvacuumregion          !< set by \c LREGION macro; true if inside the vacuum region
    LOGICAL :: Lsavedguvij            !< flag used in matrix free
    LOGICAL :: Localconstraint        !< what is this?
!> @}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> \addtogroup grp_field_matrices Field matrices: dMA, dMB, dMC, dMD, dME, dMF
!> <ul>
!> <li> The energy, \f$W \equiv \int \! dv {\; \bf B}\cdot{\bf B}\f$, and helicity, \f$K\equiv \int \! dv \; {\bf A}\cdot{\bf B}\f$, functionals may be written
!>      \f{eqnarray}{ W & = & \frac{1}{2} \; a_i \; A_{i,j} \; a_j + a_i \; B_{i,j} \; \psi_j + \frac{1}{2} \; \psi_i \; C_{i,j} \; \psi_j \label{eq:energy_globalmatrix_global} \\
!>                    K & = & \frac{1}{2} \; a_i \; D_{i,j} \; a_j + a_i \; E_{i,j} \; \psi_j + \frac{1}{2} \; \psi_i \; F_{i,j} \; \psi_j \label{eq:helicitymatrix_global}
!>      \f}
!>       where \f${\bf a} \equiv \{ {\color{red} A_{\theta,e,i,l}}, {\color{blue} A_{\zeta, e,i,l}}, {\color{Orange}  A_{\theta,o,i,l}}, {\color{Cerulean}A_{\zeta ,o,i,l}}, f_{e,i}, f_{o,i} \}\f$
!>       contains the independent degrees of freedom and \f$\boldsymbol{\psi} \equiv \{\Delta \psi_t,\Delta \psi_p\}\f$. </li>
!> <li> These are allocated and deallocated in dforce(), assigned in matrix(), and used in mp00ac() and (?) df00aa(). </li>
!> </ul>
!> @{
    real(wp), allocatable :: dMA(:, :) !< energy and helicity matrices; quadratic forms
    real(wp), allocatable :: dMB(:, :) !< energy and helicity matrices; quadratic forms
!  REAL,   allocatable :: dMC(:,:) !< energy and helicity matrices; quadratic forms
    real(wp), allocatable :: dMD(:, :) !< energy and helicity matrices; quadratic forms
!  REAL,   allocatable :: dME(:,:) !< energy and helicity matrices; quadratic forms
!  REAL,   allocatable :: dMF(:,:) !< energy and helicity matrices; quadratic forms

    real(wp), allocatable :: dMAS(:)     !< sparse version of dMA, data
    real(wp), allocatable :: dMDS(:)     !< sparse version of dMD, data
    integer, allocatable :: idMAS(:)    !< sparse version of dMA and dMD, indices
    integer, allocatable :: jdMAS(:)    !< sparse version of dMA and dMD, indices
    integer, allocatable :: NdMASmax(:) !< number of elements for sparse matrices
    integer, allocatable :: NdMAS(:)    !< number of elements for sparse matrices

    real(wp), allocatable :: dMG(:) !< what is this?

    real(wp), allocatable :: AdotX(:) !< the matrix-vector product
    real(wp), allocatable :: DdotX(:) !< the matrix-vector product

    real(wp), allocatable :: solution(:, :) !< this is allocated in dforce; used in mp00ac and ma02aa; and is passed to packab

    real(wp), allocatable :: GMRESlastsolution(:, :, :) !< used to store the last solution for restarting GMRES

    real(wp), allocatable :: MBpsi(:)      !< matrix vector products

    LOGICAL :: LILUprecond        !< whether to use ILU preconditioner for GMRES

    real(wp), allocatable :: BeltramiInverse(:, :) !< Beltrami inverse matrix

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    real(wp), allocatable :: diotadxup(:, :, :) !< measured rotational transform on inner/outer interfaces for each volume;          d(transform)/dx; (see dforce)
    real(wp), allocatable :: dItGpdxtp(:, :, :) !< measured toroidal and poloidal current on inner/outer interfaces for each volume; d(Itor,Gpol)/dx; (see dforce)

    real(wp), allocatable :: glambda(:, :, :, :) !< save initial guesses for iterative calculation of rotational-transform

    integer :: lmns !< number of independent degrees of freedom in angle transformation;

    real(wp), allocatable :: dlambdaout(:, :, :)
!> @}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> \addtogroup grp_force_constr Construction of "force"
!> The force vector is comprised of \c Bomn and \c Iomn.
!> @{
    real(wp), allocatable :: Bemn(:, :, :) !< force vector;     stellarator-symmetric (?)
    real(wp), allocatable :: Iomn(:, :)   !< force vector;     stellarator-symmetric (?)
    real(wp), allocatable :: Somn(:, :, :) !< force vector; non-stellarator-symmetric (?)
    real(wp), allocatable :: Pomn(:, :, :) !< force vector; non-stellarator-symmetric (?)

    real(wp), allocatable :: Bomn(:, :, :) !< force vector;     stellarator-symmetric (?)
    real(wp), allocatable :: Iemn(:, :)   !< force vector;     stellarator-symmetric (?)
    real(wp), allocatable :: Semn(:, :, :) !< force vector; non-stellarator-symmetric (?)
    real(wp), allocatable :: Pemn(:, :, :) !< force vector; non-stellarator-symmetric (?)

    real(wp), allocatable :: BBe(:) !< force vector (?);     stellarator-symmetric (?)
    real(wp), allocatable :: IIo(:) !< force vector (?);     stellarator-symmetric (?)
    real(wp), allocatable :: BBo(:) !< force vector (?); non-stellarator-symmetric (?)
    real(wp), allocatable :: IIe(:) !< force vector (?); non-stellarator-symmetric (?)
!> @}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!> \addtogroup grp_covar_field_ifaces Covariant field on interfaces: Btemn, Bzemn, Btomn, Bzomn
!> The covariant field
!> @{
    real(wp), allocatable :: Btemn(:, :, :) !< covariant \f$\theta\f$ cosine component of the tangential field on interfaces;     stellarator-symmetric
    real(wp), allocatable :: Bzemn(:, :, :) !< covariant \f$\zeta\f$  cosine component of the tangential field on interfaces;     stellarator-symmetric
    real(wp), allocatable :: Btomn(:, :, :) !< covariant \f$\theta\f$   sine component of the tangential field on interfaces; non-stellarator-symmetric
    real(wp), allocatable :: Bzomn(:, :, :) !< covariant \f$\zeta\f$    sine component of the tangential field on interfaces; non-stellarator-symmetric
!> @}
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> \addtogroup grp_covar_field_hessian covariant field for Hessian computation: Bloweremn, Bloweromn
!> @{
    real(wp), allocatable :: Bloweremn(:, :) !< covariant field for Hessian computation
    real(wp), allocatable :: Bloweromn(:, :) !< covariant field for Hessian computation
!> @}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!> \addtogroup grp_geomdof Geometrical degrees-of-freedom: LGdof, NGdof
!> The geometrical degrees-of-freedom
!> @{
    integer :: LGdof !<       geometrical degrees of freedom associated with each interface
    integer :: NGdof !< total geometrical degrees of freedom
!> @}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!> \addtogroup grp_par_deriv_mat Parallel construction of derivative matrix
!> <ul>
!> <li> The derivatives of force-balance, \f$[[p+B^2/2]]\f$, and the spectral constraints (see sw03aa()), with respect to the interface geometry
!>      is constructed in parallel by dforce(). </li>
!> <li> force-balance across the \f$l\f$-th interface depends on the fields in the adjacent interfaces. </li>
!> </ul>
!> @{
    real(wp), allocatable :: dBBdRZ(:, :, :) !< derivative of magnetic field w.r.t. geometry (?)
    real(wp), allocatable :: dIIdRZ(:, :) !< derivative of spectral constraints w.r.t. geometry (?)

    real(wp), allocatable :: dFFdRZ(:, :, :, :, :) !< derivatives of B^2 at the interfaces wrt geometry
    real(wp), allocatable :: dBBdmp(:, :, :, :) !< derivatives of B^2 at the interfaces wrt mu and dpflux

    real(wp), allocatable :: HdFFdRZ(:, :, :, :, :) !< derivatives of B^2 at the interfaces wrt geometry 2D Hessian;

    real(wp), allocatable :: denergydrr(:, :, :, :, :) !< derivatives of energy at the interfaces wrt geometry 3D Hessian;
    real(wp), allocatable :: denergydrz(:, :, :, :, :) !< derivatives of energy at the interfaces wrt geometry 3D Hessian;
    real(wp), allocatable :: denergydzr(:, :, :, :, :) !< derivatives of energy at the interfaces wrt geometry 3D Hessian;
    real(wp), allocatable :: denergydzz(:, :, :, :, :) !< derivatives of energy at the interfaces wrt geometry 3D Hessian;
!> @}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> \addtogroup grp_deriv_mul_polflux Derivatives of multiplier and poloidal flux with respect to geometry: dmupfdx
!> <ul>
!> <li> The information in \c dmupfdx describes how the helicity multiplier, \f$\mu\f$, and the enclosed poloidal flux, \f$\Delta \psi_p\f$,
!>      must vary as the geometry is varied in order to satisfy the interface transform constraint. </li>
!> <li> The internal variable \c dmupfdx(1:Mvol,1:2,1:LGdof,0:1) is allocated/deallocated in newton(), and hesian() if selected. </li>
!> <li> The magnetic field depends on the Fourier harmonics of both the inner and outer interface geometry (represented here as \f$x_j\f$),
!>      the helicity multiplier, and the enclosed poloidal flux, i.e. \f${\bf B_\pm} = {\bf B_\pm}(x_j, \mu, \Delta \psi_p)\f$, so that
!>      \f{eqnarray}{ \delta {\bf B_\pm} = \frac{\partial {\bf B}_\pm}{\partial x_j          } \delta x_j
!>                                       + \frac{\partial {\bf B}_\pm}{\partial \mu          } \delta \mu
!>                                       + \frac{\partial {\bf B}_\pm}{\partial \Delta \psi_p} \delta \Delta \psi_p.
!>      \f} </li>
!> <li> This information is used to adjust the calculation of how force-balance, i.e. \f$B^2\f$ at the interfaces,
!>      varies with geometry at fixed interface rotational transform. Given
!>      \f{eqnarray}{ B_\pm^2 = B_\pm^2 (x_j, \mu, \Delta \psi_p),
!>      \f}
!>      we may derive
!>      \f{eqnarray}{ \frac{\partial B_\pm^2}{\partial x_j} = \frac{\partial B_\pm^2}{\partial x_j          }
!>                                                          + \frac{\partial B_\pm^2}{\partial \mu          } \frac{\partial \mu          }{\partial x_j}
!>                                                          + \frac{\partial B_\pm^2}{\partial \Delta \psi_p} \frac{\partial \Delta \psi_p}{\partial x_j}
!>      \f} </li>
!> <li> The constraint to be enforced is that \f$\mu\f$ and \f$\Delta \psi_p\f$ must generally vary as the geometry is varied
!>      if the value of the rotational-transform constraint on the inner/outer interface is to be preserved,
!>      i.e.
!>      \f{eqnarray}{ \left(\begin{array}{ccc} \displaystyle \frac{\partial {{\,\iota\!\!\!}-}_-}{\partial {\bf B}_-} \cdot \frac{\partial {\bf B}_-}{\partial \mu          } & , &
!>                                             \displaystyle \frac{\partial {{\,\iota\!\!\!}-}_-}{\partial {\bf B}_-} \cdot \frac{\partial {\bf B}_-}{\partial \Delta \psi_p} \\
!>                                             \displaystyle \frac{\partial {{\,\iota\!\!\!}-}_+}{\partial {\bf B}_+} \cdot \frac{\partial {\bf B}_+}{\partial \mu          } & , &
!>                                             \displaystyle \frac{\partial {{\,\iota\!\!\!}-}_+}{\partial {\bf B}_+} \cdot \frac{\partial {\bf B}_+}{\partial \Delta \psi_p}
!>                    \end{array} \right)
!>                      \left(\begin{array}{c} \displaystyle \frac{\partial \mu}{\partial x_j} \\
!>                                             \displaystyle \frac{\partial \Delta \psi_p}{\partial x_j} \end{array} \right) =
!>                    - \left(\begin{array}{c} \displaystyle \frac{\partial {{\,\iota\!\!\!}-}_-}{\partial {\bf B}_-} \cdot \frac{\partial {\bf B}_-}{\partial x_j} \\
!>                                             \displaystyle \frac{\partial {{\,\iota\!\!\!}-}_+}{\partial {\bf B}_+} \cdot \frac{\partial {\bf B}_+}{\partial x_j} \end{array} \right).
!>      \f} </li>
!> <li> This \f$2\times 2\f$ linear equation is solved in dforce()
!>      and the derivatives of the rotational-transform are given in \c diotadxup, see preset.f90 . </li>
!> <li> A finite-difference estimate is computed if \c Lcheck==4. </li>
!> </ul>
!> @{
    real(wp), allocatable :: dmupfdx(:, :, :, :, :)  !< derivatives of mu and dpflux wrt geometry at constant interface transform

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    LOGICAL :: Lhessianallocated !< flag to indicate that force gradient matrix is allocated (?)
    real(wp), allocatable :: hessian(:, :)      !<               force gradient matrix (?)
    real(wp), allocatable :: dessian(:, :)      !< derivative of force gradient matrix (?)

    LOGICAL :: Lhessian2Dallocated !< flag to indicate that 2D Hessian matrix is allocated (?)
    real(wp), allocatable :: hessian2D(:, :) !< Hessian 2D
    real(wp), allocatable :: dessian2D(:, :) !< derivative Hessian 2D

    LOGICAL :: Lhessian3Dallocated !< flag to indicate that 2D Hessian matrix is allocated (?)
    real(wp), allocatable :: hessian3D(:, :) !< Hessian 3D
    real(wp), allocatable :: dessian3D(:, :) !< derivative Hessian 3D
!> @}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!> \addtogroup grp_trig Trigonometric factors
!> <ul>
!> <li> To facilitate construction of the metric integrals, various trigonometric identities are exploited. </li>
!> <li> The following are used for volume integrals (see volume() ):
!>      \f{eqnarray}{ a_{i,j,k} &=& 4 \; m_k \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \cos(\alpha_i)\cos(\alpha_j)\cos(\alpha_k) /(2\pi)^2 , \\
!>                    b_{i,j,k} &=& 4 \; m_j \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \cos(\alpha_i)\sin(\alpha_j)\sin(\alpha_k) /(2\pi)^2 ,
!>      \f} </li>
!> </ul>
!> @{
    real(wp), allocatable :: cosi(:, :) !< some precomputed cosines
    real(wp), allocatable :: sini(:, :) !< some precomputed sines
    real(wp), allocatable :: gteta(:)  !< something related to \f$\sqrt g\f$ and \f$\theta\f$ ?
    real(wp), allocatable :: gzeta(:)  !< something related to \f$\sqrt g\f$ and \f$\zeta\f$ ?

    real(wp), allocatable :: ajk(:)    !< definition of coordinate axis

    real(wp), allocatable :: dRadR(:, :, :, :) !< derivatives of coordinate axis
    real(wp), allocatable :: dRadZ(:, :, :, :) !< derivatives of coordinate axis
    real(wp), allocatable :: dZadR(:, :, :, :) !< derivatives of coordinate axis
    real(wp), allocatable :: dZadZ(:, :, :, :) !< derivatives of coordinate axis

    real(wp), allocatable :: dRodR(:, :, :) !< derivatives of coordinate axis
    real(wp), allocatable :: dRodZ(:, :, :) !< derivatives of coordinate axis
    real(wp), allocatable :: dZodR(:, :, :) !< derivatives of coordinate axis
    real(wp), allocatable :: dZodZ(:, :, :) !< derivatives of coordinate axis

    integer, allocatable :: djkp(:, :) !< for calculating cylindrical volume
    integer, allocatable :: djkm(:, :) !< for calculating cylindrical volume
!> @}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> \addtogroup grp_volints Volume integrals: lBBintegral, lABintegral
!> <ul>
!> <li> The energy functional, \f$F \equiv \sum_l F_l\f$, where
!>      \f{eqnarray}{ F_l \equiv \left( \int_{{\cal V}_l} \frac{p_l}{\gamma-1} + \frac{B_l^2}{2} dv \right)
!>                    = \frac{P_l}{\gamma-1}V_l^{1-\gamma}+\int_{{\cal V}_l} \frac{B_l^2}{2} dv, \label{eq:energy_global}
!>      \f}
!>      where the second expression is derived using \f$p_l V_l^\gamma=P_l\f$, where \f$P_l\f$ is the adiabatic-constant.
!>      In Eqn.\f$(\ref{eq:energy_global})\f$, it is implicit that \f${\bf B}\f$ satisfies (i) the toroidal and poloidal flux constraints;
!>      (ii) the interface constraint, \f${\bf B}\cdot\nabla s=0\f$; and (iii) the helicity constraint (or the transform constraint). </li>
!> <li> The derivatives of \f$F_l\f$ with respect to the inner and outer adjacent interface geometry are stored in
!>      \c dFF(1:Nvol,0:1,0:mn+mn-1), where
!>
!>      \f$         F_l                      \equiv\f$ \c dFF(l,0,    0)
!>
!>      \f$\partial F_l / \partial R_{l-1,j} \equiv\f$ \c dFF(ll,0,   j)
!>
!>      \f$\partial F_l / \partial Z_{l-1,j} \equiv\f$ \c dFF(ll,0,mn j)
!>
!>      \f$\partial F_l / \partial R_{l  ,j} \equiv\f$ \c dFF(ll,1,   j)
!>
!>      \f$\partial F_l / \partial Z_{l  ,j} \equiv\f$ \c dFF(ll,1,mn j)
!>      </li>
!> <li> The volume integrals \f$\int dv\f$, \f$\int B^2 \; dv\f$ and \f$\int {\bf A}\cdot{\bf B} \; dv\f$ in each volume
!>      are computed and saved in \c volume(0:2,1:Nvol). </li>
!> </ul>
!> @{
    real(wp), allocatable :: lBBintegral(:) !< B.B integral
    real(wp), allocatable :: lABintegral(:) !< A.B integral

    real(wp), allocatable :: vvolume(:) !< volume integral of \f$\sqrt g\f$; computed in volume
    real(wp) :: dvolume    !< derivative of volume w.r.t. interface geometry
!> @}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> \addtogroup grp_int_global Internal global variables
!> internal global variables; internal logical variables; default values are provided here; these may be changed according to input values
!> @{
    integer :: ivol !< labels volume; some subroutines (called by NAG) are fixed argument list but require the volume label

    real(wp) :: gBzeta !< toroidal (contravariant) field; calculated in bfield; required to convert \f$\dot \theta\f$ to \f$B^\theta\f$, \f$\dot s\f$ to \f$B^s\f$

    integer, allocatable :: Iquad(:) !< internal copy of Nquad

    real(wp), allocatable :: gaussianweight(:, :)    !<   weights for Gaussian quadrature
    real(wp), allocatable :: gaussianabscissae(:, :) !< abscissae for Gaussian quadrature

    LOGICAL :: LBlinear !< controls selection of Beltrami field solver; depends on LBeltrami
    LOGICAL :: LBnewton !< controls selection of Beltrami field solver; depends on LBeltrami
    LOGICAL :: LBsequad !< controls selection of Beltrami field solver; depends on LBeltrami

    real(wp) :: oRZp(1:3) !< used in mg00aa() to determine \f$(s,\theta,\zeta)\f$ given \f$(R,Z,\varphi)\f$

!> @}

    type(derivative) :: dBdX !< \f${\rm d}\mathbf{B}/{\rm d}\mathbf{X}\f$ (?)

!> \addtogroup grp_misc Miscellaneous
!> The following are miscellaneous flags required for the virtual casing field, external (vacuum) field integration, ...
!> @{
    integer :: globaljk  !< labels position
    real(wp), allocatable :: Dxyz(:, :) !< computational boundary; position
    real(wp), allocatable :: Nxyz(:, :) !< computational boundary; normal
    real(wp), allocatable :: Jxyz(:, :) !< plasma        boundary; surface current

    real(wp) :: tetazeta(1:2) !< what is this?

    real(wp) :: virtualcasingfactor = -one/(four*pi) !< this agrees with diagno

    integer :: IBerror !< for computing error in magnetic field

    integer :: nfreeboundaryiterations !< number of free-boundary iterations already performed

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    integer, parameter :: Node = 2 !< best to make this global for consistency between calling and called routines

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    LOGICAL :: first_free_bound = .false. !< flag to indicate that this is the first free-boundary iteration
!> @}
!>
!> @}

contains

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    subroutine build_vector_potential(lvol, iocons, aderiv, tderiv)
! Builds the covariant component of the vector potential and store them in efmn, ofmn, sfmn, cfmn.

        use constants, only: zero, half

        use fileunits, only: ounit

        use inputlist, only: Lrad, Wbuild_vector_potential, Wmacros

        use cputiming

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef OPENMP
        USE OMP_LIB
#endif
        use mpi
        implicit none
        integer :: ierr, astat, ios, nthreads, ithread
        real(wp) :: cput, cpui, cpuo = 0 ! cpu time; cpu initial; cpu old; 31 Jan 13;

        integer :: aderiv    ! Derivative of A. -1: w.r.t geometrical degree of freedom
        !                   0: no derivatives
        !                   1: w.r.t mu
        !                   2: w.r.t pflux
        integer :: tderiv    ! Derivative of Chebyshev polynomialc. 0: no derivatives
        !                                      1: w.r.t radial coordinate s
        integer :: ii, &    ! Loop index on Fourier harmonics
                   ll, &    ! Loop index on radial resolution
                   mi, &    ! Poloidal mode number
                   lvol, &    ! Volume number
                   iocons    ! inner (0) or outer (1) side of the volume
        real(wp) :: mfactor   ! Regularization factor when LcoordinateSingularity

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

        cpui = MPI_WTIME()
        cpuo = cpui
#ifdef OPENMP
        nthreads = omp_get_max_threads()
#else
        nthreads = 1
#endif

        efmn(1:mn) = zero
        sfmn(1:mn) = zero
        cfmn(1:mn) = zero
        ofmn(1:mn) = zero

        do ii = 1, mn ! loop over Fourier harmonics; 13 Sep 13;

            if (Lcoordinatesingularity) then
                mi = im(ii)
                do ll = mi, Lrad(lvol), 2 ! loop over Zernike polynomials; Lrad is the radial resolution; 01 Jul 19;
                    efmn(ii) = efmn(ii) + Ate(lvol, aderiv, ii)%s(ll)*RTT(ll, mi, iocons, 1)*half
                    cfmn(ii) = cfmn(ii) + Aze(lvol, aderiv, ii)%s(ll)*RTT(ll, mi, iocons, 1)*half
                    if (NOTstellsym) then
                        ofmn(ii) = ofmn(ii) + Ato(lvol, aderiv, ii)%s(ll)*RTT(ll, mi, iocons, 1)*half
                        sfmn(ii) = sfmn(ii) + Azo(lvol, aderiv, ii)%s(ll)*RTT(ll, mi, iocons, 1)*half
                    end if
                end do ! end of do ll; 20 Feb 13;
            else
                do ll = 0, Lrad(lvol) ! loop over Chebyshev polynomials; Lrad is the radial resolution;
                    efmn(ii) = efmn(ii) + Ate(lvol, aderiv, ii)%s(ll)*TT(ll, iocons, 1) ! aderiv labels deriv. wrt mu, pflux;
                    cfmn(ii) = cfmn(ii) + Aze(lvol, aderiv, ii)%s(ll)*TT(ll, iocons, 1)
                    if (NOTstellsym) then
                        ofmn(ii) = ofmn(ii) + Ato(lvol, aderiv, ii)%s(ll)*TT(ll, iocons, 1)
                        sfmn(ii) = sfmn(ii) + Azo(lvol, aderiv, ii)%s(ll)*TT(ll, iocons, 1)
                    end if
                end do ! end of do ll; 20 Feb 13;
            end if ! Lcoordinatesingularity; 01 Jul 19;
        end do ! end of do ii; 20 Feb 13;

    end subroutine build_vector_potential

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! Use this routine to run SPEC on a different MPI communicator.
! It takes care to re-assign the rank and size information
! that is used internally in SPEC.
    subroutine set_mpi_comm(comm)

        implicit none

        integer, intent(in) :: comm
        integer :: ierr

        ! MPI_COMM_SPEC is the global variable for the SPEC communicator.
        MPI_COMM_SPEC = comm

        myid = 0
        ncpu = 1

        ierr = 0
        call MPI_COMM_RANK(MPI_COMM_SPEC, myid, ierr)
        if (ierr .ne. 0) write (*, *) "error in call to MPI_COMM_RANK"

        ierr = 0
        call MPI_COMM_SIZE(MPI_COMM_SPEC, ncpu, ierr)
        if (ierr .ne. 0) write (*, *) "error in call to MPI_COMM_SIZE"

    end subroutine

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    subroutine read_inputlists_from_file()

        use constants
        use fileunits
        use inputlist

#ifdef IFORT
        use ifport ! for fseek, ftell with Intel compiler
#endif

#ifdef OPENMP
        USE OMP_LIB
#endif
        use mpi
        implicit none
        integer :: ierr, astat, ios, nthreads, ithread
        real(wp) :: cput, cpui, cpuo = 0 ! cpu time; cpu initial; cpu old; 31 Jan 13;

        LOGICAL :: Lspexist
        integer :: filepos, seek_status, cpfile, instat, idx_mode

        character(len=1000) :: line

        integer :: mm, nn
        real(wp), allocatable :: RZRZ(:, :) ! local array used for reading interface Fourier harmonics from file;

        inquire (file=trim(ext)//".sp", exist=Lspexist) ! check if file exists;

        if (.not. Lspexist) then
            write (6, '("readin :      fatal : myid=",i3," ; .not.Lspexist ; the input file does not exist ;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : .not.Lspexist : the input file does not exist  ;"
        end if
        ! if not, abort;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

        open (iunit, file=trim(ext)//".sp", status="old")

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

        instat = 0 ! initially, no error

! read namelists one after another
        if (Wreadin) then
            cput = MPI_WTIME()
            write (ounit, '("readin : ",f10.2," : reading physicslist     from ext.sp ;")') cput - cpus
        end if

        read (iunit, physicslist, iostat=instat)
        if (instat .ne. 0) then
            ! help to debug invalid inputs:
            ! re-read last line that lead to error and print it to screen
            backspace (iunit)
            read (iunit, fmt='(A)') line
            write (*, '(A)') 'Invalid line in physicslist: '//trim(line)
        end if

        if (Wreadin) then
            cput = MPI_WTIME()
            write (ounit, '("readin : ",f10.2," : read    physicslist     from ext.sp ;")') cput - cpus
        end if

        if (Wreadin) then
            cput = MPI_WTIME()
            write (ounit, '("readin : ",f10.2," : reading numericlist     from ext.sp ;")') cput - cpus
        end if

        read (iunit, numericlist, iostat=instat)
        if (instat .ne. 0) then
            ! help to debug invalid inputs:
            ! re-read last line that lead to error and print it to screen
            backspace (iunit)
            read (iunit, fmt='(A)') line
            write (*, '(A)') 'Invalid line in numericlist: '//trim(line)
        end if

        if (Wreadin) then
            cput = MPI_WTIME()
            write (ounit, '("readin : ",f10.2," : read    numericlist     from ext.sp ;")') cput - cpus
        end if

        if (Wreadin) then
            cput = MPI_WTIME()
            write (ounit, '("readin : ",f10.2," : reading locallist      from ext.sp ;")') cput - cpus
        end if

        read (iunit, locallist, iostat=instat)
        if (instat .ne. 0) then
            ! help to debug invalid inputs:
            ! re-read last line that lead to error and print it to screen
            backspace (iunit)
            read (iunit, fmt='(A)') line
            write (*, '(A)') 'Invalid line in locallist: '//trim(line)
        end if

        if (Wreadin) then
            cput = MPI_WTIME()
            write (ounit, '("readin : ",f10.2," : read    locallist      from ext.sp ;")') cput - cpus
        end if

        if (Wreadin) then
            cput = MPI_WTIME()
            write (ounit, '("readin : ",f10.2," : reading globallist   from ext.sp ;")') cput - cpus
        end if

        read (iunit, globallist, iostat=instat)
        if (instat .ne. 0) then
            ! help to debug invalid inputs:
            ! re-read last line that lead to error and print it to screen
            backspace (iunit)
            read (iunit, fmt='(A)') line
            write (*, '(A)') 'Invalid line in globallist: '//trim(line)
        end if

        if (Wreadin) then
            cput = MPI_WTIME()
            write (ounit, '("readin : ",f10.2," : read    globallist   from ext.sp ;")') cput - cpus
        end if

        if (Wreadin) then
            cput = MPI_WTIME()
            write (ounit, '("readin : ",f10.2," : reading diagnosticslist from ext.sp ;")') cput - cpus
        end if

        read (iunit, diagnosticslist, iostat=instat)
        if (instat .ne. 0) then
            ! help to debug invalid inputs:
            ! re-read last line that lead to error and print it to screen
            backspace (iunit)
            read (iunit, fmt='(A)') line
            write (*, '(A)') 'Invalid line in diagnosticslist: '//trim(line)
        end if

        if (Wreadin) then
            cput = MPI_WTIME()
            write (ounit, '("readin : ",f10.2," : read    diagnosticslist from ext.sp ;")') cput - cpus
        end if

        if (Wreadin) then
            cput = MPI_WTIME()
            write (ounit, '("readin : ",f10.2," : reading screenlist      from ext.sp ;")') cput - cpus
        end if

        read (iunit, screenlist, iostat=instat)
        if (instat .ne. 0) then
            ! help to debug invalid inputs:
            ! re-read last line that lead to error and print it to screen
            backspace (iunit)
            read (iunit, fmt='(A)') line
            write (*, '(A)') 'Invalid line in screenlist: '//trim(line)
        end if

        if (Wreadin) then
            cput = MPI_WTIME()
            write (ounit, '("readin : ",f10.2," : read    screenlist      from ext.sp ;")') cput - cpus
        end if

        ! At this point, the input namelists are read.
        ! It remains to read the initial guess for the interface geometry,
        ! which follows after the namelists in the input file.

        ! need to reset status flag for below logic to work
        instat = 0

        num_modes = 0
        if (Linitialize .le. 0) then

            ! duplicate of checks required for below code

            if (Nvol .lt. 1 .or. Nvol .gt. MNvol) then
                write (6, '("readin :      fatal : myid=",i3," ; Nvol.lt.1 .or. Nvol.gt.MNvol ; invalid Nvol: may need to recompile with higher MNvol ;")') myid
                call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                stop "readin : Nvol.lt.1 .or. Nvol.gt.MNvol : invalid Nvol: may need to recompile with higher MNvol  ;"
            end if

            allocate (RZRZ(1:4, 1:Nvol), stat=astat)
            RZRZ(1:4, 1:Nvol) = zero
            ! temp array for reading input;

            ! determine how many modes are specified by reading them once
#ifdef IFORT
            filepos = ftell(iunit) + 1
#else
            call ftell(iunit, filepos)
#endif
            do ! will read in Fourier harmonics until the end of file is reached;
                read (iunit, *, iostat=instat) mm, nn, RZRZ(1:4, 1:Nvol)   !if change of angle applies, transformation assumes m>=0 and for m=0 only n>=0;
                if (instat .ne. 0) exit

                num_modes = num_modes + 1
            end do

            ! rewind file to reset EOF flag
            rewind (iunit)

            ! seek back to (start of modes) == (end of input namelists)
#ifdef IFORT
            seek_status = fseek(iunit, filepos, 0)
#else
            call fseek(iunit, filepos, 0, seek_status)
#endif

            if (seek_status .ne. 0) then
                write (6, '("inplst :      fatal : myid=",i3," ; seek_status.ne.0 ; failed to seek to end of input namelists ;")') myid
                call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                stop "inplst : seek_status.ne.0 : failed to seek to end of input namelists  ;"
            end if

            ! now allocate arrays and read...
            ! Need to free memory, in case preset() called multiple times via python wrappers
            if (allocated(mmRZRZ)) deallocate (mmRZRZ, nnRZRZ, allRZRZ)
            allocate (mmRZRZ(1:num_modes), nnRZRZ(1:num_modes), allRZRZ(1:4, 1:Nvol, 1:num_modes))

            do idx_mode = 1, num_modes
                read (iunit, *, iostat=instat) mmRZRZ(idx_mode), nnRZRZ(idx_mode), allRZRZ(1:4, 1:Nvol, idx_mode)
            end do

            ! no need for temporary RZRZ anymore

            deallocate (RZRZ, stat=astat)

        end if ! Linitialize .le. 0

        close (iunit)

    end subroutine ! read_inputlists_from_file

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    subroutine check_inputs()

        use numerical
        use constants
        use fileunits
        use inputlist
        use cputiming, only: Treadin

#ifdef OPENMP
        USE OMP_LIB
#endif
        use mpi
        implicit none
        integer :: ierr, astat, ios, nthreads, ithread
        real(wp) :: cput, cpui, cpuo = 0 ! cpu time; cpu initial; cpu old; 31 Jan 13;

        integer :: vvol
        real(wp) :: xx, toroidalflux, toroidalcurrent

        cpui = MPI_WTIME()
        cpuo = cpui
#ifdef OPENMP
        nthreads = omp_get_max_threads()
#else
        nthreads = 1
#endif

        Mvol = Nvol + Lfreebound ! this is just for screen output and initial check; true assignment of Mvol appears outside if( myid.eq.0 ) then ;

        write (ounit, '("readin : ", 10x ," : ")')

        cput = MPI_WTIME()

        write (ounit, 1010) cput - cpus, Igeometry, Istellsym, Lreflect
        write (ounit, 1011) Lfreebound, phiedge, curtor, curpol
        write (ounit, 1012) gamma
        write (ounit, 1013) Nfp, Nvol, Mvol, Mpol, Ntor
        write (ounit, 1014) pscale, Ladiabatic, Lconstraint, mupftol, mupfits
        write (ounit, 1015) Lrad(1:min(Mvol, 32))

1010    format("readin : ", f10.2, " : Igeometry=", i3, " ; Istellsym=", i3, " ; Lreflect="i3" ;")
1011    format("readin : ", 10x, " : Lfreebound=", i3, " ; phiedge="es23.15" ; curtor="es23.15" ; curpol="es23.15" ;")
1012    format("readin : ", 10x, " : gamma="es23.15" ;")
1013    format("readin : ", 10x, " : Nfp=", i3, " ; Nvol=", i3, " ; Mvol=", i3, " ; Mpol=", i3, " ; Ntor=", i3, " ;")
1014    format("readin : ", 10x, " : pscale="es13.5" ; Ladiabatic="i2" ; Lconstraint="i3" ; mupf: tol,its="es10.2" ,"i4" ;")
1015    format("readin : ", 10x, " : Lrad = "257(i2, ",", :))

#ifdef DEBUG
        if (Wreadin) then
            write (ounit, '("readin : ",f10.2," : tflux    ="257(es11.3" ,":))') cput - cpus, (tflux(vvol), vvol=1, Mvol)
            write (ounit, '("readin : ",f10.2," : pflux    ="257(es11.3" ,":))') cput - cpus, (pflux(vvol), vvol=1, Mvol)
            write (ounit, '("readin : ",f10.2," : helicity ="256(es11.3" ,":))') cput - cpus, (helicity(vvol), vvol=1, Nvol)
            write (ounit, '("readin : ",f10.2," : pressure ="257(es11.3" ,":))') cput - cpus, (pressure(vvol), vvol=1, Mvol)
            write (ounit, '("readin : ",f10.2," : mu       ="257(es11.3" ,":))') cput - cpus, (mu(vvol), vvol=1, Mvol)
            write (ounit, '("readin : ",f10.2," : Ivolume  ="257(es11.3" ,":))') cput - cpus, (Ivolume(vvol), vvol=1, Mvol)
            write (ounit, '("readin : ",f10.2," : Isurf    ="257(es11.3" ,":))') cput - cpus, (Isurf(vvol), vvol=1, Mvol)
        end if
#endif

        if (Igeometry .lt. 1 .or. Igeometry .gt. 3) then
            write (6, '("readin :      fatal : myid=",i3," ; Igeometry.lt.1 .or. Igeometry.gt.3 ; invalid geometry ;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : Igeometry.lt.1 .or. Igeometry.gt.3 : invalid geometry  ;"
        end if

        if (Nfp .le. 0) then
            write (6, '("readin :      fatal : myid=",i3," ; Nfp.le.0 ; invalid Nfp ;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : Nfp.le.0 : invalid Nfp  ;"
        end if

        if (Mpol .lt. 0 .or. Mpol .gt. MMpol) then
            write (6, '("readin :      fatal : myid=",i3," ; Mpol.lt.0 .or. Mpol.gt.MMpol ; invalid poloidal resolution: may need to recompile with higher MMpol ;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : Mpol.lt.0 .or. Mpol.gt.MMpol : invalid poloidal resolution: may need to recompile with higher MMpol  ;"
        end if

        if (Ntor .lt. 0 .or. Ntor .gt. MNtor) then
            write (6, '("readin :      fatal : myid=",i3," ; Ntor.lt.0 .or. Ntor.gt.MNtor ; invalid toroidal resolution: may need to recompile with higher MNtor ;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : Ntor.lt.0 .or. Ntor.gt.MNtor : invalid toroidal resolution: may need to recompile with higher MNtor  ;"
        end if

        if (Nvol .lt. 1 .or. Nvol .gt. MNvol) then
            write (6, '("readin :      fatal : myid=",i3," ; Nvol.lt.1 .or. Nvol.gt.MNvol ; invalid Nvol: may need to recompile with higher MNvol ;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : Nvol.lt.1 .or. Nvol.gt.MNvol : invalid Nvol: may need to recompile with higher MNvol  ;"
        end if

        if (mupftol .le. zero) then
            write (6, '("readin :      fatal : myid=",i3," ; mupftol.le.zero ; mupftol is too small ;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : mupftol.le.zero : mupftol is too small  ;"
        end if

        if (abs(one + gamma) .lt. vsmall) then
            write (6, '("readin :      fatal : myid=",i3," ; abs(one+gamma).lt.vsmall ; 1+gamma appears in denominator in dforce ;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : abs(one+gamma).lt.vsmall : 1+gamma appears in denominator in dforce  ;"
        end if
        !< \todo Please check this; SRH: 27 Feb 18;

        if (abs(one - gamma) .lt. vsmall) then
            write (6, '("readin :      fatal : myid=",i3," ; abs(one-gamma).lt.vsmall ; 1-gamma appears in denominator in fu00aa ;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : abs(one-gamma).lt.vsmall : 1-gamma appears in denominator in fu00aa  ;"
        end if
        !< \todo Please check this; SRH: 27 Feb 18;

        if (Lconstraint .lt. -1 .or. Lconstraint .gt. 3) then
            write (6, '("readin :      fatal : myid=",i3," ; Lconstraint.lt.-1 .or. Lconstraint.gt.3 ; illegal Lconstraint ;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : Lconstraint.lt.-1 .or. Lconstraint.gt.3 : illegal Lconstraint  ;"
        end if

        if (Igeometry .eq. 1 .and. rpol .lt. vsmall) then
            write (6, '("readin :      fatal : myid=",i3," ; Igeometry.eq.1 .and. rpol.lt.vsmall ; poloidal extent of slab too small or negative ;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : Igeometry.eq.1 .and. rpol.lt.vsmall : poloidal extent of slab too small or negative  ;"
        end if

        if (Igeometry .eq. 1 .and. rtor .lt. vsmall) then
            write (6, '("readin :      fatal : myid=",i3," ; Igeometry.eq.1 .and. rtor.lt.vsmall ; toroidal extent of slab too small or negative ;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : Igeometry.eq.1 .and. rtor.lt.vsmall : toroidal extent of slab too small or negative  ;"
        end if

        if (Istellsym .eq. 1) then
            Rbs(-MNtor:MNtor, -MMpol:MMpol) = zero
            Zbc(-MNtor:MNtor, -MMpol:MMpol) = zero
            Rws(-MNtor:MNtor, -MMpol:MMpol) = zero
            Zwc(-MNtor:MNtor, -MMpol:MMpol) = zero
            Vnc(-MNtor:MNtor, -MMpol:MMpol) = zero
            Bnc(-MNtor:MNtor, -MMpol:MMpol) = zero
        end if

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **reading of physicslist**
!> <ul>
!> <li> The internal variable, \c Mvol=Nvol+Lfreebound , gives the number of computational domains. </li>
!> <li> The input value for the fluxes enclosed within each interface, \c tflux(1:Mvol) and \c tflux(1:Mvol), are immediately normalized:
!>
!>       \c tflux(1:Mvol) \f$\rightarrow\f$ \c tflux(1:Mvol)/tflux(Nvol).
!>
!>       \c pflux(1:Mvol) \f$\rightarrow\f$ \c pflux(1:Mvol)/tflux(Nvol).
!>
!>       The input \f$\Phi_{edge} \equiv \f$ \c phiedge will provide the total toroidal flux; see preset(). </li>
!> <li> The input value for the toroidal current constraint (\c Isurf(1:Mvol) and \c Ivolume(1:Mvol) ) are also immediately normalized, using \c curtor .
!>       \f$Ivolume \rightarrow Ivolume \cdot \frac{curtor}{\sum_i Isurf_i + Ivolume_i}\f$
!>       \f$Isurf   \rightarrow Isurf   \cdot \frac{curtor}{\sum_i Isurf_i + Ivolume_i}\f$
!> </ul>

        if (abs(tflux(Nvol)) .lt. vsmall) then
            write (6, '("readin :      fatal : myid=",i3," ; abs(tflux(Nvol)).lt. vsmall ; enclosed toroidal flux cannot be zero ;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : abs(tflux(Nvol)).lt. vsmall : enclosed toroidal flux cannot be zero  ;"
        end if

        toroidalflux = tflux(Nvol) ! toroidal flux is a local variable; SRH: 27 Feb 18

        tflux(1:Mvol) = tflux(1:Mvol)/toroidalflux ! normalize toroidal flux
        pflux(1:Mvol) = pflux(1:Mvol)/toroidalflux ! normalize poloidal flux

        if (tflux(1) .lt. zero) then
            write (6, '("readin :      fatal : myid=",i3," ; tflux(1).lt.zero ; enclosed toroidal flux cannot be zero ;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : tflux(1).lt.zero : enclosed toroidal flux cannot be zero  ;"
        end if

        do vvol = 2, Mvol
            !FATAL( readin, tflux(vvol)-tflux(vvol-1).lt.small, toroidal flux is not monotonic )
        end do

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **Current profiles normalization**
!>
!> In case of a free boundary calculation (\c Lfreebound=1) and using a current constraint (\c Lconstraint=3),
!> the current profiles are renormalized in order to match the linking current \c curtor.
!> More specifically,
!> \f{eqnarray}{
!> Isurf_i & \rightarrow\ Isurf_i \cdot \frac{curtor}{\sum_{i=1}^{Mvol-1} Isurf_i+Ivol_i}
!> Ivol_i  & \rightarrow\ Ivol_i  \cdot \frac{curtor}{\sum_{i=1}^{Mvol-1} Isurf_i+Ivol_i}
!> \f}
!> Finally, the volume current in the vacuum region is set to \f$0\f$.

        if ((Lfreebound .EQ. 1) .and. (Lconstraint .EQ. 3)) then

            Ivolume(Mvol) = Ivolume(Mvol - 1) ! Ensure vacuum in vacuum region

            toroidalcurrent = Ivolume(Mvol) + sum(Isurf(1:Mvol - 1))

            if (curtor .NE. 0) then

                if (toroidalcurrent .EQ. 0) then
                    write (6, '("readin :      fatal : myid=",i3," ; toroidalcurrent.EQ.0  ; Incompatible current profiles and toroidal linking current;")') myid
                    call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                    stop "readin : toroidalcurrent.EQ.0  : Incompatible current profiles and toroidal linking current ;"
                end if

                Ivolume(1:Mvol) = Ivolume(1:Mvol)*curtor/toroidalcurrent
                Isurf(1:Mvol - 1) = Isurf(1:Mvol - 1)*curtor/toroidalcurrent

            else

                if (toroidalcurrent .NE. 0) then
                    write (6, '("readin :      fatal : myid=",i3," ; toroidalcurrent.NE.0 ; Incompatible current profiles and toroidal linking current;")') myid
                    call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                    stop "readin : toroidalcurrent.NE.0 : Incompatible current profiles and toroidal linking current ;"
                end if

                ! No rescaling if profiles have an overall zero toroidal current
            end if
        end if

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

        do vvol = 1, Mvol

            if (Lrad(vvol) .lt. 2) then
                write (6, '("readin :      fatal : myid=",i3," ; Lrad(vvol ).lt.2 ; require Chebyshev resolution Lrad > 2 so that Lagrange constraints can be satisfied ;")') myid
                call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                stop "readin : Lrad(vvol ).lt.2 : require Chebyshev resolution Lrad > 2 so that Lagrange constraints can be satisfied  ;"
            end if

        end do

        if (Igeometry .ge. 2 .and. Lrad(1) .lt. Mpol) then
            write (ounit, '("readin : ",f10.2," : Minimum Lrad(1) is Mpol, automatically adjusted it to Mpol+4")') cput - cpus
            Lrad(1) = Mpol + 4
        end if

        if (mupfits .le. 0) then
            write (6, '("readin :      fatal : myid=",i3," ; mupfits.le.0 ; must give ma01aa:hybrj a postive integer value for the maximum iterations = mupfits given on input ;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : mupfits.le.0 : must give ma01aa:hybrj a postive integer value for the maximum iterations = mupfits given on input  ;"
        end if

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **reading of numericlist**

        write (ounit, '("readin : ", 10x ," : ")')

        write (ounit, 1020) cput - cpus, Linitialize, LautoinitBn, Lzerovac, Ndiscrete
        write (ounit, 1021) Nquad, iMpol, iNtor
        write (ounit, 1022) Lsparse, Lsvdiota, imethod, iorder, iprecon, iotatol
        write (ounit, 1023) Lextrap, Mregular, Lrzaxis, Ntoraxis

1020    format("readin : ", f10.2, " : Linitialize=", i3, " ;LautoinitBn=", i3, " ; Lzerovac=", i2, " ; Ndiscrete="i2" ;")
1021    format("readin : ", 10x, " : Nquad="i4" ; iMpol="i4" ; iNtor="i4" ;")
1022    format("readin : ", 10x, " : Lsparse="i2" ; Lsvdiota="i2" ; imethod="i2" ; iorder="i2" ; iprecon="i2" ; iotatol="es13.5" ;")
1023    format("readin : ", 10x, " : Lextrap="i2" ; Mregular="i3" ; Lrzaxis="i2" ; Ntoraxis="i2" ;")

        if (Ndiscrete .le. 0) then
            write (6, '("readin :      fatal : myid=",i3," ; Ndiscrete.le.0 ; error ;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : Ndiscrete.le.0 : error  ;"
        end if

        !FATAL(readin, Lfreebound.eq.1 .and. Lconstraint.gt.0 .and. Lsparse.eq.0, have not implemented dense Fourier angle transformation in vacuum region )

        if (iotatol .gt. one) then
            write (6, '("readin :      fatal : myid=",i3," ; iotatol.gt.one ; illegal value for sparse tolerance ;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : iotatol.gt.one : illegal value for sparse tolerance  ;"
        end if
        ! I think that the sparse iota solver is no longer implemented; SRH: 27 Feb 18;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **reading of locallist**

        write (ounit, '("readin : ", 10x ," : ")')

        if (LBeltrami .ne. 4 .and. Lmatsolver .ne. 1) then
            write (ounit, '("readin : ", 10x ," : ***Lmatsolver set to 1 for nonlinear solver***")')
            Lmatsolver = 1
        end if

        write (ounit, 1030) cput - cpus, LBeltrami, Linitgues, Lmatsolver, LGMRESprec, NiterGMRES, epsGMRES, epsILU

1030    format("readin : ", f10.2, " : LBeltrami="i2" ; Linitgues="i2" ; Lmatsolver="i2" ; LGMRESprec="i2" ; NiterGMRES="i4" ; epsGMRES="es13.5" ; epsILU="es13.5" ;")

        if (LBeltrami .lt. 0 .or. LBeltrami .gt. 7) then
            write (6, '("readin :      fatal : myid=",i3," ; LBeltrami.lt.0 .or. LBeltrami.gt.7 ; error ;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : LBeltrami.lt.0 .or. LBeltrami.gt.7 : error  ;"
        end if

        if (Lmatsolver .lt. 0 .or. Lmatsolver .gt. 3) then
            write (6, '("readin :      fatal : myid=",i3," ; Lmatsolver.lt.0 .or. Lmatsolver.gt.3 ; error ;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : Lmatsolver.lt.0 .or. Lmatsolver.gt.3 : error  ;"
        end if

        if (LGMRESprec .lt. 0 .or. LGMRESprec .gt. 1) then
            write (6, '("readin :      fatal : myid=",i3," ; LGMRESprec.lt.0 .or. LGMRESprec.gt.1 ; error ;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : LGMRESprec.lt.0 .or. LGMRESprec.gt.1 : error  ;"
        end if

        if (NiterGMRES .lt. 0) then
            write (6, '("readin :      fatal : myid=",i3," ; NiterGMRES.lt.0 ; error ;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : NiterGMRES.lt.0 : error  ;"
        end if

        if (abs(epsGMRES) .le. machprec) then
            write (6, '("readin :      fatal : myid=",i3," ; abs(epsGMRES).le.machprec  ; error ;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : abs(epsGMRES).le.machprec  : error  ;"
        end if

        if (abs(epsILU) .le. machprec) then
            write (6, '("readin :      fatal : myid=",i3," ; abs(epsILU).le.machprec  ; error ;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : abs(epsILU).le.machprec  : error  ;"
        end if

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **reading of globallist**

        write (ounit, '("readin : ", 10x ," : ")')

        write (ounit, 1040) cput - cpus, Lfindzero
        write (ounit, 1041) escale, opsilon, pcondense, epsilon, wpoloidal, upsilon
        write (ounit, 1042) forcetol, c05xmax, c05xtol, c05factor, LreadGF
        write (ounit, 1043) mfreeits, gBntol, gBnbld
        write (ounit, 1044) vcasingeps, vcasingtol, vcasingits, vcasingper

1040    format("readin : ", f10.2, " : Lfindzero="i2" ;")
1041    format("readin : ", 10x, " : escale="es13.5" ; opsilon="es13.5" ; pcondense="f7.3" ; epsilon="es13.5" ; wpoloidal="f7.4" ; upsilon="es13.5" ;")
1042    format("readin : ", 10x, " : forcetol="es13.5" ; c05xmax="es13.5" ; c05xtol="es13.5" ; c05factor="es13.5" ; LreadGF="L2" ; ")
1043    format("readin : ", 10x, " : mfreeits="i4" ; gBntol="es13.5" ; gBnbld="es13.5" ;")
1044    format("readin : ", 10x, " : vcasingeps="es13.5" ; vcasingtol="es13.5" ; vcasingits="i6" ; vcasingper="i6" ;")

        if (escale .lt. zero) then
            write (6, '("readin :      fatal : myid=",i3," ; escale      .lt.zero      ; error ;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : escale      .lt.zero      : error  ;"
        end if

        if (pcondense .lt. one) then
            write (6, '("readin :      fatal : myid=",i3," ; pcondense   .lt.one       ; error ;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : pcondense   .lt.one       : error  ;"
        end if

        if (abs(c05xtol) .le. machprec) then
            write (6, '("readin :      fatal : myid=",i3," ; abs(c05xtol).le.machprec  ; error ;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : abs(c05xtol).le.machprec  : error  ;"
        end if

        if (c05factor .le. zero) then
            write (6, '("readin :      fatal : myid=",i3," ; c05factor   .le.zero      ; error ;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : c05factor   .le.zero      : error  ;"
        end if

        !FATAL( readin, mfreeits    .lt.zero     , error )

        if (Igeometry .eq. 3 .and. pcondense .le. zero) then
            write (6, '("readin :      fatal : myid=",i3," ; Igeometry.eq.3 .and. pcondense.le.zero ; pcondense must be positive ;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : Igeometry.eq.3 .and. pcondense.le.zero : pcondense must be positive  ;"
        end if

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **reading of diagnosticslist**

        write (ounit, '("readin : ", 10x ," : ")')

        write (ounit, 1050) cput - cpus, odetol, nPpts
        write (ounit, 1051) LHevalues, LHevectors, LHmatrix, Lperturbed, dpp, dqq, dRZ, Lcheck, Ltiming

1050    format("readin : ", f10.2, " : odetol="es10.2" ; nPpts="i6" ;")
1051    format("readin : ", 10x, " : LHevalues="L2" ; LHevectors="L2" ; LHmatrix="L2" ; Lperturbed="i2" ; dpp="i3" ; dqq="i3" ; dRZ="es16.8" ; Lcheck="i3" ; Ltiming="L2" ;")

        if (odetol .le. zero) then
            write (6, '("readin :      fatal : myid=",i3," ; odetol.le.zero ; input error ;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : odetol.le.zero : input error  ;"
        end if

        !FATAL( readin, absreq.le.zero, input error )
        !FATAL( readin, relreq.le.zero, input error )
        !FATAL( readin, absacc.le.zero, input error )
        !FATAL( readin, nPpts .lt.0   , input error )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **reading of screenlist**

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

        write (ounit, '("readin : ", 10x ," : ")')

9999    continue
        cput = MPI_WTIME()
        Treadin = Treadin + (cput - cpuo)
        return

    end subroutine ! check_inputs

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    subroutine broadcast_inputs

        use fileunits
        use inputlist

#ifdef OPENMP
        USE OMP_LIB
#endif
        use mpi
        implicit none
        integer :: ierr, astat, ios, nthreads, ithread
        real(wp) :: cput, cpui, cpuo = 0 ! cpu time; cpu initial; cpu old; 31 Jan 13;

        call MPI_BCAST(ext, 100, MPI_CHARACTER, 0, MPI_COMM_SPEC, ierr)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **broadcast physicslist**

        if (Wreadin) then; cput = MPI_WTIME(); write (ounit, '("readin : ",f10.2," : broadcasting physicslist     from ext.sp ;")') cput - cpus
        end if

        call MPI_BCAST(Igeometry, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(Istellsym, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(Lfreebound, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(phiedge, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(curtor, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(curpol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(gamma, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(Nfp, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(Nvol, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(Mpol, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(Ntor, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(Lrad, MNvol + 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(tflux, MNvol + 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(pflux, MNvol + 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(helicity, MNvol, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(pscale, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(pressure, MNvol + 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(Ladiabatic, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(adiabatic, MNvol + 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(mu, MNvol + 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(Ivolume, MNvol + 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(Isurf, MNvol + 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(Lconstraint, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(pl, MNvol, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(ql, MNvol, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(pr, MNvol, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(qr, MNvol, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(iota, MNvol, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(lp, MNvol, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(lq, MNvol, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(rp, MNvol, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(rq, MNvol, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(oita, MNvol, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(mupftol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(mupfits, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(Lreflect, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(rpol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(rtor, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **broadcast numericlist**

        if (Wreadin) then; cput = MPI_WTIME(); write (ounit, '("readin : ",f10.2," : broadcasting numericlist     from ext.sp ;")') cput - cpus
        end if

        call MPI_BCAST(Linitialize, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(LautoinitBn, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(Lzerovac, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(Ndiscrete, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(Nquad, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(iMpol, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(iNtor, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(Lsparse, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(Lsvdiota, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(imethod, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(iorder, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(iprecon, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(iotatol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(Lextrap, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(Mregular, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(Lrzaxis, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(Ntoraxis, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **broadcast globallist**

        if (Wreadin) then; cput = MPI_WTIME(); write (ounit, '("readin : ",f10.2," : broadcasting globallist      from ext.sp ;")') cput - cpus
        end if

        call MPI_BCAST(Lfindzero, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(escale, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(opsilon, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(pcondense, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(epsilon, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(wpoloidal, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(upsilon, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(forcetol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(c05xmax, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(c05xtol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(c05factor, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(LreadGF, 1, MPI_LOGICAL, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(mfreeits, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(gBntol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(gBnbld, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(vcasingeps, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(vcasingtol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(vcasingits, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(vcasingper, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **broadcast locallist**

        if (Wreadin) then
            cput = MPI_WTIME()
            write (ounit, '("readin : ",f10.2," : broadcasting locallist       from ext.sp ;")') cput - cpus
        end if

        call MPI_BCAST(LBeltrami, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(Linitgues, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(maxrndgues, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(Lmatsolver, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(NiterGMRES, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(epsGMRES, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(LGMRESprec, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(epsILU, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

! IlBCAST( Lposdef  , 1, 0 ) ! redundant;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **broadcast diagnosticslist**

        if (Wreadin) then; cput = MPI_WTIME(); write (ounit, '("readin : ",f10.2," : broadcasting diagnosticslist from ext.sp ;")') cput - cpus
        end if

        call MPI_BCAST(odetol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

        !RlBCAST( absreq    , 1      , 0 )
        !RlBCAST( relreq    , 1      , 0 )
        !RlBCAST( absacc    , 1      , 0 )
        !RlBCAST( epsr      , 1      , 0 )

        call MPI_BCAST(nPpts, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(Ppts, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(nPtrj, MNvol + 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(LHevalues, 1, MPI_LOGICAL, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(LHevectors, 1, MPI_LOGICAL, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(Ltransform, 1, MPI_LOGICAL, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(LHmatrix, 1, MPI_LOGICAL, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(Lperturbed, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(dpp, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(dqq, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(Lerrortype, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(Ngrid, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(dRZ, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(Lcheck, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(Ltiming, 1, MPI_LOGICAL, 0, MPI_COMM_SPEC, ierr)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **broadcast screenlist**

        if (Wreadin) then; cput = MPI_WTIME(); write (ounit, '("readin : ",f10.2," : broadcasting screenlist      from ext.sp ;")') cput - cpus
        end if

! BSCREENLIST ! broadcast screenlist; this is expanded by Makefile; do not remove;

        call MPI_BCAST(Wreadin, 1, MPI_LOGICAL, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(Wwrtend, 1, MPI_LOGICAL, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(Wmacros, 1, MPI_LOGICAL, 0, MPI_COMM_SPEC, ierr)

    end subroutine ! broadcast_inputs

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> \brief The restart file is written.
    subroutine wrtend

        use constants, only:

        use numerical, only: machprec

        use fileunits, only: ounit, iunit

        use cputiming, only: Twrtend

        use inputlist

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef OPENMP
        USE OMP_LIB
#endif
        use mpi
        implicit none
        integer :: ierr, astat, ios, nthreads, ithread
        real(wp) :: cput, cpui, cpuo = 0 ! cpu time; cpu initial; cpu old; 31 Jan 13;

        integer :: vvol !< iteration variable over all nested volumes
        integer :: imn  !< iteration variable for all Fourier harmonics
        integer :: ii   !< iteration variable for all Fourier harmonics
        integer :: jj   !< iteration variable
        integer :: kk   !< iteration variable
        integer :: jk   !< iteration variable
        integer :: Lcurvature !< curvature flag (?)
        integer :: mm   !< current poloidal mode number
        integer :: nn   !< current toroidal mode number

        real(wp) :: lss !< (?)
        real(wp) :: teta !< (?)
        real(wp) :: zeta !< (?)
        real(wp) :: st(1:Node) !< (?)
        real(wp) :: Bst(1:Node) !< (?)
        real(wp) :: BR !< (?)
        real(wp) :: BZ !< (?)
        real(wp) :: BP !< (?)

        cpui = MPI_WTIME()
        cpuo = cpui
#ifdef OPENMP
        nthreads = omp_get_max_threads()
#else
        nthreads = 1
#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

        if (myid .ne. 0) goto 9999

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef DEBUG
        if (Wwrtend) then; cput = MPI_WTIME(); write (ounit, '("wrtend : ",f10.2," : myid=",i3," ; opening/writing ext.sp.end ;")') cput - cpus, myid
        end if
#endif

        open (iunit, file=trim(ext)//".sp.end", status="unknown") ! restart input file;

#ifdef DEBUG
        if (Wwrtend) then; cput = MPI_WTIME(); write (ounit, '("wrtend : ",f10.2," : myid=",i3," ; writing physicslist ;")') cput - cpus, myid
        end if
#endif

        write (iunit, '("&physicslist")')
        write (iunit, '(" Igeometry   = ",i9        )') Igeometry
        write (iunit, '(" Istellsym   = ",i9        )') Istellsym
        write (iunit, '(" Lfreebound  = ",i9        )') Lfreebound
        write (iunit, '(" phiedge     = ",es23.15   )') phiedge
        write (iunit, '(" curtor      = ",es23.15   )') curtor
        write (iunit, '(" curpol      = ",es23.15   )') curpol
        write (iunit, '(" gamma       = ",es23.15   )') gamma
        write (iunit, '(" Nfp         = ",i9        )') Nfp
        write (iunit, '(" Nvol        = ",i9        )') Nvol
        write (iunit, '(" Mpol        = ",i9        )') Mpol
        write (iunit, '(" Ntor        = ",i9        )') Ntor
        write (iunit, '(" Lrad        = ",257i23    )') Lrad(1:Mvol)
        write (iunit, '(" tflux       = ",257es23.15)') tflux(1:Mvol)
        write (iunit, '(" pflux       = ",257es23.15)') pflux(1:Mvol)
        write (iunit, '(" helicity    = ",256es23.15)') helicity(1:Mvol)
        write (iunit, '(" pscale      = ",es23.15   )') pscale
        write (iunit, '(" Ladiabatic  = ",i9        )') Ladiabatic
        write (iunit, '(" pressure    = ",257es23.15)') pressure(1:Mvol)
        write (iunit, '(" adiabatic   = ",257es23.15)') adiabatic(1:Mvol)
        write (iunit, '(" mu          = ",257es23.15)') mu(1:Mvol)
        write (iunit, '(" Ivolume     = ",257es23.15)') Ivolume(1:Mvol)
        write (iunit, '(" Isurf       = ",257es23.15)') Isurf(1:Mvol - 1)
        write (iunit, '(" Lconstraint = ",i9        )') Lconstraint
        write (iunit, '(" pl          = ",257i23    )') pl(0:Mvol)
        write (iunit, '(" ql          = ",257i23    )') ql(0:Mvol)
        write (iunit, '(" pr          = ",257i23    )') pr(0:Mvol)
        write (iunit, '(" qr          = ",257i23    )') qr(0:Mvol)
        write (iunit, '(" iota        = ",257es23.15)') iota(0:Mvol)
        write (iunit, '(" lp          = ",257i23    )') lp(0:Mvol)
        write (iunit, '(" lq          = ",257i23    )') lq(0:Mvol)
        write (iunit, '(" rp          = ",257i23    )') rp(0:Mvol)
        write (iunit, '(" rq          = ",257i23    )') rq(0:Mvol)
        write (iunit, '(" oita        = ",257es23.15)') oita(0:Mvol)
        write (iunit, '(" mupftol     = ",es23.15   )') mupftol
        write (iunit, '(" mupfits     = ",i9        )') mupfits
        write (iunit, '(" Lreflect    = ",i9        )') Lreflect
        write (iunit, '(" rpol        = ",es23.15   )') rpol
        write (iunit, '(" rtor        = ",es23.15   )') rtor

        if (Lfreebound .eq. 1 .or. Zbs(0, 1) .gt. zero) then
            do ii = 1, mn
                mm = im(ii)
                nn = in(ii)/Nfp

                ! geometry of boundary
                Rbc(nn, mm) = iRbc(ii, Nvol)
                Zbs(nn, mm) = iZbs(ii, Nvol)
                Rbs(nn, mm) = iRbs(ii, Nvol)
                Zbc(nn, mm) = iZbc(ii, Nvol)

                ! vacuum and plasma fields
                Vns(nn, mm) = iVns(ii)
                Bns(nn, mm) = iBns(ii)
                Vnc(nn, mm) = iVnc(ii)
                Bnc(nn, mm) = iBnc(ii)

                ! geometry of wall
                Rwc(nn, mm) = iRbc(ii, Mvol)
                Zws(nn, mm) = iZbs(ii, Mvol)
                Rws(nn, mm) = iRbs(ii, Mvol)
                Zwc(nn, mm) = iZbc(ii, Mvol)
            end do ! end of do ii = 1, mn;
        end if ! end of if( Lfreebound.eq.1 .or. . . . ) ;

        !write(iunit,'(" Rac         = ",99es23.15)') Rac(0:Ntor)
        !write(iunit,'(" Zas         = ",99es23.15)') Zas(0:Ntor)
        !write(iunit,'(" Ras         = ",99es23.15)') Ras(0:Ntor)
        !write(iunit,'(" Zac         = ",99es23.15)') Zac(0:Ntor)

        ! geometry of axis
        write (iunit, '(" Rac         = ",99es23.15)') iRbc(1:Ntor + 1, 0)
        write (iunit, '(" Zas         = ",99es23.15)') iZbs(1:Ntor + 1, 0)
        write (iunit, '(" Ras         = ",99es23.15)') iRbs(1:Ntor + 1, 0)
        write (iunit, '(" Zac         = ",99es23.15)') iZbc(1:Ntor + 1, 0)

        do mm = 0, Mpol ! will write out the plasma boundary harmonics;
            do nn = -Ntor, Ntor

                if (mm .eq. 0 .and. nn .lt. 0) cycle ! these modes are always excluded; 13 Oct 12;

                select case (mm)
                case (0:9)
                    if (nn .lt. -9 .and. nn .gt. -99) write (iunit, 1000) nn, mm, Rbc(nn, mm), nn, mm, Zbs(nn, mm), nn, mm, Rbs(nn, mm), nn, mm, Zbc(nn, mm)
                    if (nn .lt. 0 .and. nn .ge. -9) write (iunit, 1001) nn, mm, Rbc(nn, mm), nn, mm, Zbs(nn, mm), nn, mm, Rbs(nn, mm), nn, mm, Zbc(nn, mm)
                    if (nn .ge. 0 .and. nn .le. 9) write (iunit, 1002) nn, mm, Rbc(nn, mm), nn, mm, Zbs(nn, mm), nn, mm, Rbs(nn, mm), nn, mm, Zbc(nn, mm)
                    if (nn .gt. 9 .and. nn .le. 99) write (iunit, 1001) nn, mm, Rbc(nn, mm), nn, mm, Zbs(nn, mm), nn, mm, Rbs(nn, mm), nn, mm, Zbc(nn, mm)
                case (10:99)
                    if (nn .lt. -9 .and. nn .gt. -99) write (iunit, 1003) nn, mm, Rbc(nn, mm), nn, mm, Zbs(nn, mm), nn, mm, Rbs(nn, mm), nn, mm, Zbc(nn, mm)
                    if (nn .lt. 0 .and. nn .ge. -9) write (iunit, 1004) nn, mm, Rbc(nn, mm), nn, mm, Zbs(nn, mm), nn, mm, Rbs(nn, mm), nn, mm, Zbc(nn, mm)
                    if (nn .ge. 0 .and. nn .le. 9) write (iunit, 1005) nn, mm, Rbc(nn, mm), nn, mm, Zbs(nn, mm), nn, mm, Rbs(nn, mm), nn, mm, Zbc(nn, mm)
                    if (nn .gt. 9 .and. nn .le. 99) write (iunit, 1004) nn, mm, Rbc(nn, mm), nn, mm, Zbs(nn, mm), nn, mm, Rbs(nn, mm), nn, mm, Zbc(nn, mm)
                end select ! end of select case( mm );

            end do ! end of do nn;
        end do ! end of do mm;

        do mm = 0, Mpol ! will write out the computation domain harmonics; (only relevant in free-boundary case);
            do nn = -Ntor, Ntor

                if (mm .eq. 0 .and. nn .lt. 0) cycle ! these modes are always excluded; 13 Oct 12;

                select case (mm)
                case (0:9)
                    if (nn .lt. -9 .and. nn .gt. -99) write (iunit, 1010) nn, mm, Rwc(nn, mm), nn, mm, Zws(nn, mm), nn, mm, Rws(nn, mm), nn, mm, Zwc(nn, mm)
                    if (nn .lt. 0 .and. nn .ge. -9) write (iunit, 1011) nn, mm, Rwc(nn, mm), nn, mm, Zws(nn, mm), nn, mm, Rws(nn, mm), nn, mm, Zwc(nn, mm)
                    if (nn .ge. 0 .and. nn .le. 9) write (iunit, 1012) nn, mm, Rwc(nn, mm), nn, mm, Zws(nn, mm), nn, mm, Rws(nn, mm), nn, mm, Zwc(nn, mm)
                    if (nn .gt. 9 .and. nn .le. 99) write (iunit, 1011) nn, mm, Rwc(nn, mm), nn, mm, Zws(nn, mm), nn, mm, Rws(nn, mm), nn, mm, Zwc(nn, mm)
                case (10:99)
                    if (nn .lt. -9 .and. nn .gt. -99) write (iunit, 1013) nn, mm, Rwc(nn, mm), nn, mm, Zws(nn, mm), nn, mm, Rws(nn, mm), nn, mm, Zwc(nn, mm)
                    if (nn .lt. 0 .and. nn .ge. -9) write (iunit, 1014) nn, mm, Rwc(nn, mm), nn, mm, Zws(nn, mm), nn, mm, Rws(nn, mm), nn, mm, Zwc(nn, mm)
                    if (nn .ge. 0 .and. nn .le. 9) write (iunit, 1015) nn, mm, Rwc(nn, mm), nn, mm, Zws(nn, mm), nn, mm, Rws(nn, mm), nn, mm, Zwc(nn, mm)
                    if (nn .gt. 9 .and. nn .le. 99) write (iunit, 1014) nn, mm, Rwc(nn, mm), nn, mm, Zws(nn, mm), nn, mm, Rws(nn, mm), nn, mm, Zwc(nn, mm)
                end select ! end of select case( mm );

            end do ! end of do nn;
        end do ! end of do mm;

        do mm = 0, Mpol ! will write out the computation domain harmonics; (only relevant in free-boundary case);
            do nn = -Ntor, Ntor

                if (mm .eq. 0 .and. nn .lt. 0) cycle ! these modes are always excluded; 13 Oct 12;

                select case (mm)
                case (0:9)
                    if (nn .lt. -9 .and. nn .gt. -99) write (iunit, 1020) nn, mm, Vns(nn, mm), nn, mm, Bns(nn, mm), nn, mm, Vnc(nn, mm), nn, mm, Bnc(nn, mm)
                    if (nn .lt. 0 .and. nn .ge. -9) write (iunit, 1021) nn, mm, Vns(nn, mm), nn, mm, Bns(nn, mm), nn, mm, Vnc(nn, mm), nn, mm, Bnc(nn, mm)
                    if (nn .ge. 0 .and. nn .le. 9) write (iunit, 1022) nn, mm, Vns(nn, mm), nn, mm, Bns(nn, mm), nn, mm, Vnc(nn, mm), nn, mm, Bnc(nn, mm)
                    if (nn .gt. 9 .and. nn .le. 99) write (iunit, 1021) nn, mm, Vns(nn, mm), nn, mm, Bns(nn, mm), nn, mm, Vnc(nn, mm), nn, mm, Bnc(nn, mm)
                case (10:99)
                    if (nn .lt. -9 .and. nn .gt. -99) write (iunit, 1023) nn, mm, Vns(nn, mm), nn, mm, Bns(nn, mm), nn, mm, Vnc(nn, mm), nn, mm, Bnc(nn, mm)
                    if (nn .lt. 0 .and. nn .ge. -9) write (iunit, 1024) nn, mm, Vns(nn, mm), nn, mm, Bns(nn, mm), nn, mm, Vnc(nn, mm), nn, mm, Bnc(nn, mm)
                    if (nn .ge. 0 .and. nn .le. 9) write (iunit, 1025) nn, mm, Vns(nn, mm), nn, mm, Bns(nn, mm), nn, mm, Vnc(nn, mm), nn, mm, Bnc(nn, mm)
                    if (nn .gt. 9 .and. nn .le. 99) write (iunit, 1024) nn, mm, Vns(nn, mm), nn, mm, Bns(nn, mm), nn, mm, Vnc(nn, mm), nn, mm, Bnc(nn, mm)
                end select ! end of select case( mm );

            end do ! end of do nn;
        end do ! end of do mm;

1000    format("Rbc(", i3, ",", i1, ")", 2x, "=", es23.15, " Zbs(", i3, ",", i1, ")", 2x, "=", es23.15, " Rbs(", i3, ",", i1, ")", 2x, "=", es23.15, " Zbc(", i3, ",", i1, ")", 2x, "=", es23.15)
1001    format("Rbc(", i2, ",", i1, ")", 3x, "=", es23.15, " Zbs(", i2, ",", i1, ")", 3x, "=", es23.15, " Rbs(", i2, ",", i1, ")", 3x, "=", es23.15, " Zbc(", i2, ",", i1, ")", 3x, "=", es23.15)
1002    format("Rbc(", i1, ",", i1, ")", 4x, "=", es23.15, " Zbs(", i1, ",", i1, ")", 4x, "=", es23.15, " Rbs(", i1, ",", i1, ")", 4x, "=", es23.15, " Zbc(", i1, ",", i1, ")", 4x, "=", es23.15)
1003    format("Rbc(", i3, ",", i2, ")", 1x, "=", es23.15, " Zbs(", i3, ",", i2, ")", 1x, "=", es23.15, " Rbs(", i3, ",", i2, ")", 1x, "=", es23.15, " Zbc(", i3, ",", i2, ")", 1x, "=", es23.15)
1004    format("Rbc(", i2, ",", i2, ")", 2x, "=", es23.15, " Zbs(", i2, ",", i2, ")", 2x, "=", es23.15, " Rbs(", i2, ",", i2, ")", 2x, "=", es23.15, " Zbc(", i2, ",", i2, ")", 2x, "=", es23.15)
1005    format("Rbc(", i1, ",", i2, ")", 3x, "=", es23.15, " Zbs(", i1, ",", i2, ")", 3x, "=", es23.15, " Rbs(", i1, ",", i2, ")", 3x, "=", es23.15, " Zbc(", i1, ",", i2, ")", 3x, "=", es23.15)

1010    format("Rwc(", i3, ",", i1, ")", 2x, "=", es23.15, " Zws(", i3, ",", i1, ")", 2x, "=", es23.15, " Rws(", i3, ",", i1, ")", 2x, "=", es23.15, " Zwc(", i3, ",", i1, ")", 2x, "=", es23.15)
1011    format("Rwc(", i2, ",", i1, ")", 3x, "=", es23.15, " Zws(", i2, ",", i1, ")", 3x, "=", es23.15, " Rws(", i2, ",", i1, ")", 3x, "=", es23.15, " Zwc(", i2, ",", i1, ")", 3x, "=", es23.15)
1012    format("Rwc(", i1, ",", i1, ")", 4x, "=", es23.15, " Zws(", i1, ",", i1, ")", 4x, "=", es23.15, " Rws(", i1, ",", i1, ")", 4x, "=", es23.15, " Zwc(", i1, ",", i1, ")", 4x, "=", es23.15)
1013    format("Rwc(", i3, ",", i2, ")", 1x, "=", es23.15, " Zws(", i3, ",", i2, ")", 1x, "=", es23.15, " Rws(", i3, ",", i2, ")", 1x, "=", es23.15, " Zwc(", i3, ",", i2, ")", 1x, "=", es23.15)
1014    format("Rwc(", i2, ",", i2, ")", 2x, "=", es23.15, " Zws(", i2, ",", i2, ")", 2x, "=", es23.15, " Rws(", i2, ",", i2, ")", 2x, "=", es23.15, " Zwc(", i2, ",", i2, ")", 2x, "=", es23.15)
1015    format("Rwc(", i1, ",", i2, ")", 3x, "=", es23.15, " Zws(", i1, ",", i2, ")", 3x, "=", es23.15, " Rws(", i1, ",", i2, ")", 3x, "=", es23.15, " Zwc(", i1, ",", i2, ")", 3x, "=", es23.15)

1020    format("Vns(", i3, ",", i1, ")", 2x, "=", es23.15, " Bns(", i3, ",", i1, ")", 2x, "=", es23.15, " Vnc(", i3, ",", i1, ")", 2x, "=", es23.15, " Bnc(", i3, ",", i1, ")", 2x, "=", es23.15)
1021    format("Vns(", i2, ",", i1, ")", 3x, "=", es23.15, " Bns(", i2, ",", i1, ")", 3x, "=", es23.15, " Vnc(", i2, ",", i1, ")", 3x, "=", es23.15, " Bnc(", i2, ",", i1, ")", 3x, "=", es23.15)
1022    format("Vns(", i1, ",", i1, ")", 4x, "=", es23.15, " Bns(", i1, ",", i1, ")", 4x, "=", es23.15, " Vnc(", i1, ",", i1, ")", 4x, "=", es23.15, " Bnc(", i1, ",", i1, ")", 4x, "=", es23.15)
1023    format("Vns(", i3, ",", i2, ")", 1x, "=", es23.15, " Bns(", i3, ",", i2, ")", 1x, "=", es23.15, " Vnc(", i3, ",", i2, ")", 1x, "=", es23.15, " Bnc(", i3, ",", i2, ")", 1x, "=", es23.15)
1024    format("Vns(", i2, ",", i2, ")", 2x, "=", es23.15, " Bns(", i2, ",", i2, ")", 2x, "=", es23.15, " Vnc(", i2, ",", i2, ")", 2x, "=", es23.15, " Bnc(", i2, ",", i2, ")", 2x, "=", es23.15)
1025    format("Vns(", i1, ",", i2, ")", 3x, "=", es23.15, " Bns(", i1, ",", i2, ")", 3x, "=", es23.15, " Vnc(", i1, ",", i2, ")", 3x, "=", es23.15, " Bnc(", i1, ",", i2, ")", 3x, "=", es23.15)

        write (iunit, '("/")')

        if (Wwrtend) then; cput = MPI_WTIME(); write (ounit, '("wrtend : ",f10.2," : myid=",i3," ; writing numericlist ;")') cput - cpus, myid
        end if

        write (iunit, '("&numericlist")')
        write (iunit, '(" Linitialize = ",i9            )') Linitialize
        write (iunit, '(" LautoinitBn = ",i9            )') LautoinitBn
        write (iunit, '(" Lzerovac    = ",i9            )') Lzerovac
        write (iunit, '(" Ndiscrete   = ",i9            )') Ndiscrete
        write (iunit, '(" Nquad       = ",i9            )') Nquad
        write (iunit, '(" iMpol       = ",i9            )') iMpol
        write (iunit, '(" iNtor       = ",i9            )') iNtor
        write (iunit, '(" Lsparse     = ",i9            )') Lsparse
        write (iunit, '(" Lsvdiota    = ",i9            )') Lsvdiota
        write (iunit, '(" imethod     = ",i9            )') imethod
        write (iunit, '(" iorder      = ",i9            )') iorder
        write (iunit, '(" iprecon     = ",i9            )') iprecon
        write (iunit, '(" iotatol     = ",es23.15       )') iotatol
        write (iunit, '(" Lextrap     = ",i9            )') Lextrap
        write (iunit, '(" Mregular    = ",i9            )') Mregular
        write (iunit, '(" Lrzaxis     = ",i9            )') Lrzaxis
        write (iunit, '(" Ntoraxis    = ",i9            )') Ntoraxis
        write (iunit, '("/")')

        if (Wwrtend) then; cput = MPI_WTIME(); write (ounit, '("wrtend : ",f10.2," : myid=",i3," ; writing locallist ;")') cput - cpus, myid
        end if

        write (iunit, '("&locallist")')
        write (iunit, '(" LBeltrami   = ",i9            )') LBeltrami
        write (iunit, '(" Linitgues   = ",i9            )') Linitgues
        write (iunit, '(" Lmatsolver  = ",i9            )') Lmatsolver
        write (iunit, '(" NiterGMRES  = ",i9            )') NiterGMRES
        write (iunit, '(" LGMRESprec  = ",i9            )') LGMRESprec
        write (iunit, '(" epsGMRES    = ",es23.15       )') epsGMRES
        write (iunit, '(" epsILU      = ",es23.15       )') epsILU

        !write(iunit,'(" Lposdef     = ",i9            )') Lposdef ! redundant;
        !write(iunit,'(" Nmaxexp     = ",i9            )') Nmaxexp
        write (iunit, '("/")')

        if (Wwrtend) then; cput = MPI_WTIME(); write (ounit, '("wrtend : ",f10.2," : myid=",i3," ; writing globallist ;")') cput - cpus, myid
        end if

        write (iunit, '("&globallist")')
        write (iunit, '(" Lfindzero   = ",i9            )') Lfindzero
        write (iunit, '(" escale      = ",es23.15       )') escale
        write (iunit, '(" opsilon     = ",es23.15       )') opsilon
        write (iunit, '(" pcondense   = ",es23.15       )') pcondense
        write (iunit, '(" epsilon     = ",es23.15       )') epsilon
        write (iunit, '(" wpoloidal   = ",es23.15       )') wpoloidal
        write (iunit, '(" upsilon     = ",es23.15       )') upsilon
        write (iunit, '(" forcetol    = ",es23.15       )') forcetol
        write (iunit, '(" c05xmax     = ",es23.15       )') c05xmax
        write (iunit, '(" c05xtol     = ",es23.15       )') c05xtol
        write (iunit, '(" c05factor   = ",es23.15       )') c05factor
        write (iunit, '(" LreadGF     = ",L9            )') LreadGF
        write (iunit, '(" mfreeits    = ",i9            )') mfreeits
        write (iunit, '(" gBntol      = ",es23.15       )') gBntol
        write (iunit, '(" gBnbld      = ",es23.15       )') gBnbld
        write (iunit, '(" vcasingeps  = ",es23.15       )') vcasingeps
        write (iunit, '(" vcasingtol  = ",es23.15       )') vcasingtol
        write (iunit, '(" vcasingits  = ",i9            )') vcasingits
        write (iunit, '(" vcasingper  = ",i9            )') vcasingper
        write (iunit, '("/")')

        if (Wwrtend) then; cput = MPI_WTIME(); write (ounit, '("wrtend : ",f10.2," : myid=",i3," ; writing diagnosticslist ;")') cput - cpus, myid
        end if

        write (iunit, '("&diagnosticslist")')
        write (iunit, '(" odetol      = ",es23.15       )') odetol
        !write(iunit,'(" absreq      = ",es23.15       )') absreq
        !write(iunit,'(" relreq      = ",es23.15       )') relreq
        !write(iunit,'(" absacc      = ",es23.15       )') absacc
        !write(iunit,'(" epsr        = ",es23.15       )') epsr
        write (iunit, '(" nPpts       = ",i9            )') nPpts
        write (iunit, '(" Ppts        = ",es23.15       )') Ppts
        write (iunit, '(" nPtrj       = ",256i6         )') nPtrj(1:Mvol)
        write (iunit, '(" LHevalues   = ",L9            )') LHevalues
        write (iunit, '(" LHevectors  = ",L9            )') LHevectors
        write (iunit, '(" LHmatrix    = ",L9            )') LHmatrix
        write (iunit, '(" Lperturbed  = ",i9            )') Lperturbed
        write (iunit, '(" dpp         = ",i9            )') dpp
        write (iunit, '(" dqq         = ",i9            )') dqq
        write (iunit, '(" dRZ         = ",es23.15       )') dRZ
        write (iunit, '(" Lcheck      = ",i9            )') Lcheck
        write (iunit, '(" Ltiming     = ",L9            )') Ltiming
        write (iunit, '("/")')

        if (Wwrtend) then; cput = MPI_WTIME(); write (ounit, '("wrtend : ",f10.2," : myid=",i3," ; writing screenlist ;")') cput - cpus, myid
        end if

        write (iunit, '("&screenlist")')
! WSCREENLIST ! write screenlist; this is expanded by Makefile ; do not remove;
        if (Wreadin) write (iunit, '(" Wreadin = ",L1                )') Wreadin
        if (Wwrtend) write (iunit, '(" Wwrtend = ",L1                )') Wwrtend
        if (Wmacros) write (iunit, '(" Wmacros = ",L1                )') Wmacros
        write (iunit, '("/")')

#ifdef DEBUG

        if (.not. allocated(iRbc)) then
            write (6, '("wrtend :      fatal : myid=",i3," ; .not.allocated(iRbc) ; error ;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "wrtend : .not.allocated(iRbc) : error  ;"
        end if

        if (.not. allocated(iZbs)) then
            write (6, '("wrtend :      fatal : myid=",i3," ; .not.allocated(iZbs) ; error ;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "wrtend : .not.allocated(iZbs) : error  ;"
        end if

        if (.not. allocated(iRbs)) then
            write (6, '("wrtend :      fatal : myid=",i3," ; .not.allocated(iRbs) ; error ;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "wrtend : .not.allocated(iRbs) : error  ;"
        end if

        if (.not. allocated(iZbc)) then
            write (6, '("wrtend :      fatal : myid=",i3," ; .not.allocated(iZbc) ; error ;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "wrtend : .not.allocated(iZbc) : error  ;"
        end if

#endif

        ! write initial guess of interface geometry
        do imn = 1, mn; write (iunit, '(2i6,1024es23.15)') im(imn), in(imn)/Nfp, (iRbc(imn, vvol), iZbs(imn, vvol), iRbs(imn, vvol), iZbc(imn, vvol), vvol=1, Nvol)
        end do

        close (iunit)

#ifdef DEBUG
        if (Wwrtend) then; cput = MPI_WTIME(); write (ounit, '("wrtend : ",f10.2," : myid=",i3," ; wrote ext.sp.end ;")') cput - cpus, myid
        end if
#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

9999    continue
        cput = MPI_WTIME()
        Twrtend = Twrtend + (cput - cpuo)
        return

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    end subroutine wrtend

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> \brief Check if volume vvol is associated to the corresponding MPI node.
!>
!> The global variable \c IsMyVolumeValue is updated to 0 or 1,
!> depending on \c vvol . A value of -1 is set if an error occured.
!>
!> @param vvol volume to check
    subroutine IsMyVolume(vvol)

#ifdef OPENMP
        USE OMP_LIB
#endif
        use mpi
        implicit none
        integer :: ierr, astat, ios, nthreads, ithread
        real(wp) :: cput, cpui, cpuo = 0 ! cpu time; cpu initial; cpu old; 31 Jan 13;

        integer, intent(in) :: vvol

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

        IsMyVolumeValue = -1 ! Error value - Problem with vvol / id
        if (myid .ne. modulo(vvol - 1, ncpu)) then
            IsMyVolumeValue = 0
        else
            IsMyVolumeValue = 1
        end if

    end subroutine IsMyVolume

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> \brief Returns which MPI node is associated to a given volume.
    subroutine WhichCpuID(vvol, cpu_id)

#ifdef OPENMP
        USE OMP_LIB
#endif
        use mpi
        implicit none
        integer :: ierr, astat, ios, nthreads, ithread
        real(wp) :: cput, cpui, cpuo = 0 ! cpu time; cpu initial; cpu old; 31 Jan 13;

        integer :: vvol, cpu_id

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

        cpu_id = modulo(vvol - 1, ncpu)

    end subroutine WhichCpuID

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end module allglobal

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> \brief Interface to FFTW library
module fftw_interface ! JAB; 25 Jul 17
    use mod_kinds, only: wp => dp
    use, intrinsic :: iso_c_binding

    implicit none

    include 'fftw3.f03'

    TYPE(C_PTR) :: planf        !< FFTW-related (?)
    TYPE(C_PTR) :: planb        !< FFTW-related (?)
    COMPLEX(C_DOUBLE_COMPLEX), allocatable :: cplxin(:, :, :)  !< FFTW-related (?)
    COMPLEX(C_DOUBLE_COMPLEX), allocatable :: cplxout(:, :, :) !< FFTW-related (?)

end module fftw_interface

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
