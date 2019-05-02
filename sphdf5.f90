!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (output) ! Writes all the output information to ext.h5.

!latex \briefly{All the input and output information is contained in \type{ext.h5}.}
!latex \calledby{\link{xspech}}
!l tex \calls{\link{}}

!latex \tableofcontents

!latex \newcommand{\pb}[1]{\parbox[t]{13cm}{#1}}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

module sphdf5

  use inputlist , only : ext
  use fileunits , only : ounit
  use allglobal , only : myid
  use hdf5

  implicit none

  INTEGER                        :: hdfier, rank  ! error flag for HDF5 library
  integer(hid_t)                 :: file_id, space_id, dset_id
  integer(hsize_t)               :: onedims(1:1), twodims(1:2), threedims(1:3)

contains

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
subroutine init_outfile

  LOCALS
  integer(HID_T) :: plist_id      ! Property list identifier

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  ! initialize Fortran interface to the HDF5 library;
  call h5open_f( hdfier )
  FATAL( sphdf5, hdfier.ne.0, error calling h5open_f )

  ! Create file access property list to be able to tell HDF5 about MPI I/O
  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdfier)
  FATAL( sphdf5, hdfier.ne.0, h5pcreate_f returned an error )

  ! enable MPI I/O for parallel I/O access in file access property list
  call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, hdfier)
  FATAL( sphdf5, hdfier.ne.0, h5pset_fapl_mpio_f returned an error )

  ! Create the file (collectively, since called by all processes)
  call h5fcreate_f(trim(ext)//".h5", H5F_ACC_TRUNC_F, file_id, hdfier, access_prp = plist_id)
  FATAL( sphdf5, hdfier.ne.0, h5fcreate_f returned an error )

  ! file access property list is not needed after been used to specify MPI I/O during opening of file
  call h5pclose_f(plist_id, hdfier)
  FATAL( sphdf5, hdfier.ne.0, h5pclose_f returned an error )

end subroutine init_outfile

! mirror input variables into output file
subroutine mirror_input_to_outfile

  use inputlist
  use allglobal , only : Mvol

  LOCALS

  integer(hid_t) :: grpInput
  integer(hid_t) :: grpInputPhysics, grpInputNumerics, grpInputLocal, grpInputGlobal, grpInputDiagnostics

  HDEFGRP( file_id, input, grpInput )

! the following variables constitute the namelist/physicslist/; note that all variables in namelist need to be broadcasted in readin;
! they go into ext.h5/input/physics

  HDEFGRP( grpInput, physics, grpInputPhysics )

  HWRITEIV( grpInputPhysics,           1, Igeometry         , (/ Igeometry      /))
  HWRITEIV( grpInputPhysics,           1, Istellsym         , (/ Istellsym      /))
  HWRITEIV( grpInputPhysics,           1, Lfreebound        , (/ Lfreebound     /))
  HWRITERV( grpInputPhysics,           1, phiedge           , (/ phiedge        /))
  HWRITERV( grpInputPhysics,           1, curtor            , (/ curtor         /))
  HWRITERV( grpInputPhysics,           1, curpol            , (/ curpol         /))
  HWRITERV( grpInputPhysics,           1, gamma             , (/ gamma          /))
  HWRITEIV( grpInputPhysics,           1, Nfp               , (/ Nfp            /))
  HWRITEIV( grpInputPhysics,           1, Nvol              , (/ Nvol           /))
  HWRITEIV( grpInputPhysics,           1, Mpol              , (/ Mpol           /))
  HWRITEIV( grpInputPhysics,           1, Ntor              , (/ Ntor           /))
  HWRITEIV( grpInputPhysics,        Mvol, Lrad              ,      Lrad(1:Mvol)   )
  HWRITEIV( grpInputPhysics,           1, Lconstraint       , (/ Lconstraint    /))
  HWRITERV( grpInputPhysics,        Mvol, tflux             ,     tflux(1:Mvol)   )
  HWRITERV( grpInputPhysics,        Mvol, pflux             ,     pflux(1:Mvol)   )
  HWRITERV( grpInputPhysics,        Nvol, helicity          ,  helicity(1:Nvol)   )
  HWRITERV( grpInputPhysics,           1, pscale            , (/ pscale         /))
  HWRITERV( grpInputPhysics,        Nvol, pressure          ,  pressure(1:Nvol)   )
  HWRITEIV( grpInputPhysics,           1, Ladiabatic        , (/ Ladiabatic     /))
  HWRITERV( grpInputPhysics,        Mvol, adiabatic         , adiabatic(1:Nvol)   )
  HWRITERV( grpInputPhysics,    (1+Nvol), mu                ,        mu(1:Mvol)   )
  HWRITEIV( grpInputPhysics,    (1+Mvol), pl                ,        pl(0:Nvol)   )
  HWRITEIV( grpInputPhysics,    (1+Mvol), ql                ,        ql(0:Nvol)   )
  HWRITEIV( grpInputPhysics,    (1+Mvol), pr                ,        pr(0:Nvol)   )
  HWRITEIV( grpInputPhysics,    (1+Mvol), qr                ,        qr(0:Nvol)   )
  HWRITERV( grpInputPhysics,    (1+Nvol), iota              ,      iota(0:Nvol)   )
  HWRITEIV( grpInputPhysics,    (1+Mvol), lp                ,        lp(0:Nvol)   )
  HWRITEIV( grpInputPhysics,    (1+Mvol), lq                ,        lq(0:Nvol)   )
  HWRITEIV( grpInputPhysics,    (1+Mvol), rp                ,        rp(0:Nvol)   )
  HWRITEIV( grpInputPhysics,    (1+Mvol), rq                ,        rq(0:Nvol)   )
  HWRITERV( grpInputPhysics,    (1+Nvol), oita              ,      oita(0:Nvol)   )

  HWRITERV( grpInputPhysics,    (1+Ntor), Rac               ,       Rac(0:Ntor)   ) !     stellarator symmetric coordinate axis;
  HWRITERV( grpInputPhysics,    (1+Ntor), Zas               ,       Zas(0:Ntor)   )
  HWRITERV( grpInputPhysics,    (1+Ntor), Ras               ,       Ras(0:Ntor)   ) ! non-stellarator symmetric coordinate axis;
  HWRITERV( grpInputPhysics,    (1+Ntor), Zac               ,       Zac(0:Ntor)   )

  HWRITERA( grpInputPhysics,   (2*Ntor+1), (2*Mpol+1),  Rbc, Rbc(-Ntor:Ntor,-Mpol:Mpol) ) !     stellarator symmetric boundary components;
  HWRITERA( grpInputPhysics,   (2*Ntor+1), (2*Mpol+1),  Zbs, Zbs(-Ntor:Ntor,-Mpol:Mpol) ) !     stellarator symmetric boundary components;
  HWRITERA( grpInputPhysics,   (2*Ntor+1), (2*Mpol+1),  Rbs, Rbs(-Ntor:Ntor,-Mpol:Mpol) ) ! non-stellarator symmetric boundary components;
  HWRITERA( grpInputPhysics,   (2*Ntor+1), (2*Mpol+1),  Zbc, Zbc(-Ntor:Ntor,-Mpol:Mpol) ) ! non-stellarator symmetric boundary components;

  HWRITERA( grpInputPhysics,   (2*Ntor+1), (2*Mpol+1),  Rwc, Rwc(-Ntor:Ntor,-Mpol:Mpol) ) !     stellarator symmetric boundary components of wall;
  HWRITERA( grpInputPhysics,   (2*Ntor+1), (2*Mpol+1),  Zws, Zws(-Ntor:Ntor,-Mpol:Mpol) ) !     stellarator symmetric boundary components of wall;
  HWRITERA( grpInputPhysics,   (2*Ntor+1), (2*Mpol+1),  Rws, Rws(-Ntor:Ntor,-Mpol:Mpol) ) ! non-stellarator symmetric boundary components of wall;
  HWRITERA( grpInputPhysics,   (2*Ntor+1), (2*Mpol+1),  Zwc, Zwc(-Ntor:Ntor,-Mpol:Mpol) ) ! non-stellarator symmetric boundary components of wall;

  HWRITERA( grpInputPhysics,   (2*Ntor+1), (2*Mpol+1),  Vns, Vns(-Ntor:Ntor,-Mpol:Mpol) ) !     stellarator symmetric normal field at boundary; vacuum component;
  HWRITERA( grpInputPhysics,   (2*Ntor+1), (2*Mpol+1),  Bns, Bns(-Ntor:Ntor,-Mpol:Mpol) ) !     stellarator symmetric normal field at boundary; plasma component;
  HWRITERA( grpInputPhysics,   (2*Ntor+1), (2*Mpol+1),  Vnc, Vnc(-Ntor:Ntor,-Mpol:Mpol) ) ! non-stellarator symmetric normal field at boundary; vacuum component;
  HWRITERA( grpInputPhysics,   (2*Ntor+1), (2*Mpol+1),  Bnc, Bnc(-Ntor:Ntor,-Mpol:Mpol) ) ! non-stellarator symmetric normal field at boundary; plasma component;

  HWRITERV( grpInputPhysics,           1, mupftol           , (/ mupftol        /))
  HWRITEIV( grpInputPhysics,           1, mupfits           , (/ mupfits        /))

  HCLOSEGRP( grpInputPhysics )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! the following variables constitute the namelist/numericlist/; note that all variables in namelist need to be broadcasted in readin;
! they go into ext.h5/input/numerics

  HDEFGRP( grpInput, numerics, grpInputNumerics )

  HWRITEIV( grpInputNumerics,          1, Linitialize        , (/ Linitialize /))
  HWRITEIV( grpInputNumerics,          1, Lzerovac           , (/ Lzerovac    /))
  HWRITEIV( grpInputNumerics,          1, Ndiscrete          , (/ Ndiscrete   /))
  HWRITEIV( grpInputNumerics,          1, Nquad              , (/ Nquad       /))
  HWRITEIV( grpInputNumerics,          1, iMpol              , (/ iMpol       /))
  HWRITEIV( grpInputNumerics,          1, iNtor              , (/ iNtor       /))
  HWRITEIV( grpInputNumerics,          1, Lsparse            , (/ Lsparse     /))
  HWRITEIV( grpInputNumerics,          1, Lsvdiota           , (/ Lsvdiota    /))
  HWRITEIV( grpInputNumerics,          1, imethod            , (/ imethod     /))
  HWRITEIV( grpInputNumerics,          1, iorder             , (/ iorder      /))
  HWRITEIV( grpInputNumerics,          1, iprecon            , (/ iprecon     /))
  HWRITERV( grpInputNumerics,          1, iotatol            , (/ iotatol     /))
  HWRITEIV( grpInputNumerics,          1, Lextrap            , (/ Lextrap     /))
  HWRITEIV( grpInputNumerics,          1, Mregular           , (/ Mregular    /))

  HCLOSEGRP( grpInputNumerics )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! the following variables constitute the namelist/locallist/; note that all variables in namelist need to be broadcasted in readin;
! they go into ext.h5/input/local

  HDEFGRP( grpInput, local, grpInputLocal )

  HWRITEIV( grpInputLocal,             1, LBeltrami          , (/ LBeltrami   /))
  HWRITEIV( grpInputLocal,             1, Linitgues          , (/ Linitgues   /))
  HWRITEIV( grpInputLocal,             1, Lposdef            , (/ Lposdef     /)) ! redundant;
  HWRITERV( grpInputLocal,             1, maxrndgues         , (/ maxrndgues  /))

  HCLOSEGRP( grpInputLocal )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! the following variables constitute the namelist/globallist/; note that all variables in namelist need to be broadcasted in readin;
! they go into ext.h5/input/global

  HDEFGRP( grpInput, global, grpInputGlobal )

  HWRITEIV( grpInputGlobal,            1,  Lfindzero         , (/ Lfindzero   /))
  HWRITERV( grpInputGlobal,            1,  escale            , (/ escale      /))
  HWRITERV( grpInputGlobal,            1,  opsilon           , (/ opsilon     /))
  HWRITERV( grpInputGlobal,            1,  pcondense         , (/ pcondense   /))
  HWRITERV( grpInputGlobal,            1,  epsilon           , (/ epsilon     /))
  HWRITERV( grpInputGlobal,            1,  wpoloidal         , (/ wpoloidal   /))
  HWRITERV( grpInputGlobal,            1,  upsilon           , (/ upsilon     /))
  HWRITERV( grpInputGlobal,            1,  forcetol          , (/ forcetol    /))
  HWRITERV( grpInputGlobal,            1,  c05xmax           , (/ c05xmax     /))
  HWRITERV( grpInputGlobal,            1,  c05xtol           , (/ c05xtol     /))
  HWRITERV( grpInputGlobal,            1,  c05factor         , (/ c05factor   /))
  HWRITELV( grpInputGlobal,            1,  LreadGF           , (/ LreadGF     /))
  HWRITEIV( grpInputGlobal,            1,  mfreeits          , (/ mfreeits    /))
  HWRITERV( grpInputGlobal,            1,  bnstol            , (/ bnstol      /))  ! redundant;
  HWRITERV( grpInputGlobal,            1,  bnsblend          , (/ bnsblend    /))  ! redundant;
  HWRITERV( grpInputGlobal,            1,  gBntol            , (/ gBntol      /))
  HWRITERV( grpInputGlobal,            1,  gBnbld            , (/ gBnbld      /))
  HWRITERV( grpInputGlobal,            1,  vcasingeps        , (/ vcasingeps  /))
  HWRITERV( grpInputGlobal,            1,  vcasingtol        , (/ vcasingtol  /))
  HWRITEIV( grpInputGlobal,            1,  vcasingits        , (/ vcasingits  /))
  HWRITEIV( grpInputGlobal,            1,  vcasingper        , (/ vcasingper  /))
  HWRITEIV( grpInputGlobal,            1,  mcasingcal        , (/ mcasingcal  /))  ! redundant;

  HCLOSEGRP( grpInputGlobal )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! the following variables constitute the namelist/diagnosticslist/; note that all variables in namelist need to be broadcasted in readin;
! they go into ext.h5/input/diagnostics

  HDEFGRP( grpInput, diagnostics, grpInputDiagnostics )

  HWRITERV( grpInputDiagnostics,       1,  odetol            , (/ odetol         /))
  HWRITERV( grpInputDiagnostics,       1,  absreq            , (/ absreq         /))           ! redundant;
  HWRITERV( grpInputDiagnostics,       1,  relreq            , (/ relreq         /))           ! redundant;
  HWRITERV( grpInputDiagnostics,       1,  absacc            , (/ absacc         /))           ! redundant;
  HWRITERV( grpInputDiagnostics,       1,  epsr              , (/ epsr           /))           ! redundant;
  HWRITEIV( grpInputDiagnostics,       1,  nPpts             , (/ nPpts          /))
  HWRITEIV( grpInputDiagnostics,       1,  nPtrj             ,    nPtrj(1:Nvol+1)  )
  HWRITELV( grpInputDiagnostics,       1,  LHevalues         , (/ LHevalues      /))
  HWRITELV( grpInputDiagnostics,       1,  LHevectors        , (/ LHevectors     /))
  HWRITELV( grpInputDiagnostics,       1,  LHmatrix          , (/ LHmatrix       /))
  HWRITEIV( grpInputDiagnostics,       1,  Lperturbed        , (/ Lperturbed     /))
  HWRITEIV( grpInputDiagnostics,       1,  dpp               , (/ dpp            /))
  HWRITEIV( grpInputDiagnostics,       1,  dqq               , (/ dqq            /))
  HWRITEIV( grpInputDiagnostics,       1,  Lcheck            , (/ Lcheck         /))
  HWRITELV( grpInputDiagnostics,       1,  Ltiming           , (/ Ltiming        /))
  HWRITERV( grpInputDiagnostics,       1,  fudge             , (/ fudge          /))         ! redundant;
  HWRITERV( grpInputDiagnostics,       1,  scaling           , (/ scaling        /))          ! redundant;

  HCLOSEGRP( grpInputDiagnostics )

  HCLOSEGRP( grpInput )

end subroutine mirror_input_to_outfile



!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
subroutine hdfint

  use fileunits, only : ounit
  use inputlist
  use allglobal, only : ncpu, cpus, &
                        Mvol, ForceErr, &
                        mn, im, in, iRbc, iZbs, iRbs, iZbc, &
                        dRbc, dZbs, dRbs, dZbc, &
                        vvolume, dvolume, &
                        Bsupumn, Bsupvmn, &
                        Btemn, Bzemn, Btomn, Bzomn, &
                        iVns, iBns, iVnc, iBnc, &
                        lmns, &
                        TT, &
                        beltramierror

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS

  INTEGER                        :: Mrad
  REAL                           :: tvolume

  integer(hid_t) :: grpOutput

  HDEFGRP( file_id, output, grpOutput )

  HWRITERV( grpOutput,           mn,               Vns,       iVns(1:mn)   ) !     stellarator symmetric normal field at boundary; vacuum component;
  HWRITERV( grpOutput,           mn,               Bns,       iBns(1:mn)   ) !     stellarator symmetric normal field at boundary; plasma component;
  HWRITERV( grpOutput,           mn,               Vnc,       iVnc(1:mn)   ) ! non-stellarator symmetric normal field at boundary; vacuum component;
  HWRITERV( grpOutput,           mn,               Bnc,       iBnc(1:mn)   ) ! non-stellarator symmetric normal field at boundary; plasma component;

!latex \begin{enumerate}

!latex \item In addition to the input variables, which are described in \link{global}, the following quantities are written to \type{ext.h5} :
!latex
!latex \begin{tabular}{|l|l|l|} \hline

!latex \type{variable}               & type    & \pb{description} \\ \hline

!latex \type{mn}                     & integer & \pb{number of Fourier modes} \\
  HWRITEIV( grpOutput, 1, mn, (/ mn /)  )
!latex \type{im(1:mn)}               & integer & \pb{poloidal mode numbers} \\
  HWRITEIV( grpOutput, mn, im, im(1:mn) )
!latex \type{in(1:mn)}               & integer & \pb{toroidal mode numbers} \\
  HWRITEIV( grpOutput, mn, in, in(1:mn) )
!latex \type{Mvol}                   & integer & \pb{number of interfaces = number of volumes} \\
  HWRITEIV( grpOutput, 1, Mvol, (/ Mvol /))
!latex \type{iRbc(1:mn,0:Mvol)}      & real    & \pb{Fourier harmonics, $R_{m,n}$, of interfaces} \\
  HWRITERA( grpOutput, mn, (Mvol+1), Rbc, iRbc(1:mn,0:Mvol) )
!latex \type{iZbs(1:mn,0:Mvol)}      & real    & \pb{Fourier harmonics, $Z_{m,n}$, of interfaces} \\
  HWRITERA( grpOutput, mn, (Mvol+1), Zbs, iZbs(1:mn,0:Mvol) )
!latex \type{iRbs(1:mn,0:Mvol)}      & real    & \pb{Fourier harmonics, $R_{m,n}$, of interfaces} \\
  HWRITERA( grpOutput, mn, (Mvol+1), Rbs, iRbs(1:mn,0:Mvol) )
!latex \type{iZbc(1:mn,0:Mvol)}      & real    & \pb{Fourier harmonics, $Z_{m,n}$, of interfaces} \\
  HWRITERA( grpOutput, mn, (Mvol+1), Zbc, iZbc(1:mn,0:Mvol) )
!latex \type{forcetol}               & real    & \pb{force-balance error across interfaces} \\
  HWRITERV( grpOutput, 1, forcetol, (/ forcetol /))
!latex \type{ForceErr}               & real    & \pb{force-balance error across interfaces} \\
  HWRITERV( grpOutput, 1, ForceErr, (/ ForceErr /))

  if( Lcheck.eq.1 ) then
!latex \type{beltramierror}          & real    & \pb{error in beltrami field (volume integral)} \\
   HWRITERA( grpOutput, Mvol, 3, beltramierror, beltramierror(1:Mvol,1:3) )
  endif

  if( allocated(vvolume) ) then ! why is it required to confirm that vvolume has been allocated ; 24 Nov 16;

   tvolume = sum(vvolume(1:Nvol) )
!latex \type{volume}                 & real    & \pb{total volume = $\sum V_v$} \\
   HWRITERV( grpOutput, 1, volume, (/ tvolume /))

  else

   if( Wsphdf5 ) write(ounit,'("hdfint : ", 10x ," : myid=",i3," ; vvolume is not allocated ;")') myid

  endif ! end of if( allocated(vvolume) ) ; 11 Aug 14;

  Mrad  = maxval( Lrad(1:Mvol) )
!latex \type{Mrad}                   & integer & \pb{the maximum radial (Chebyshev) resolution} \\
  HWRITEIV( grpOutput, 1, Mrad, (/ Mrad /))
!latex \type{TT(0:Mrad,0:1,0:1)}     & real    & \pb{the Chebyshev polynomials, $T_l$, and their derivatives, evaluated at $s=\pm 1$} \\
  HWRITERC( grpOutput, (Mrad+1), 2, 2, TT, TT(0:Mrad,0:1,0:1) )
!latex \type{Btemn(1:mn,0:1,1:Mvol)} & real    & \pb{the cosine harmonics of the covariant poloidal field, \\
!latex                                           i.e. $[[B_{\t,j}]]$ evaluated on the inner and outer interface in each volume} \\
  HWRITERC( grpOutput, mn, 2, Mvol, Btemn, Btemn(1:mn,0:1,1:Mvol) )
!latex \type{Bzemn(1:mn,0:1,1:Mvol)} & real    & \pb{the cosine harmonics of the covariant toroidal field, \\
!latex                                           i.e. $[[B_{\z,j}]]$ evaluated on the inner and outer interface in each volume} \\
  HWRITERC( grpOutput, mn, 2, Mvol, Bzemn, Bzemn(1:mn,0:1,1:Mvol) )
!latex \type{Btomn(1:mn,0:1,1:Mvol)} & real    & \pb{the sine harmonics of the covariant poloidal field, \\
!latex                                           i.e. $[[B_{\t,j}]]$ evaluated on the inner and outer interface in each volume} \\
  HWRITERC( grpOutput, mn, 2, Mvol, Btomn, Btomn(1:mn,0:1,1:Mvol) )
!latex \type{Bzomn(1:mn,0:1,1:Mvol)} & real    & \pb{the sine harmonics of the covariant toroidal field, \\
!latex                                           i.e. $[[B_{\z,j}]]$ evaluated on the inner and outer interface in each volume} \\
  HWRITERC( grpOutput, mn, 2, Mvol, Bzomn, Bzomn(1:mn,0:1,1:Mvol) )

  if( Lperturbed.eq.1 ) then

!latex \type{dRbc(1:mn,0:Nvol)}      & real    & \pb{Fourier harmonics, $R_{j}$, of interfaces; linearly perturbed solution} \\
  HWRITERA( grpOutput, mn, (Nvol+1), dRbc, dRbc(1:mn,0:Nvol) )
!latex \type{dZbs(1:mn,0:Nvol)}      & real    & \pb{Fourier harmonics, $Z_{j}$, of interfaces; linearly perturbed solution} \\
  HWRITERA( grpOutput, mn, (Nvol+1), dZbs, dZbs(1:mn,0:Nvol) )
!latex \type{dRbs(1:mn,0:Nvol)}      & real    & \pb{Fourier harmonics, $R_{j}$, of interfaces; linearly perturbed solution} \\
  HWRITERA( grpOutput, mn, (Nvol+1), dRbs, dRbs(1:mn,0:Nvol) )
!latex \type{dZbc(1:mn,0:Nvol)}      & real    & \pb{Fourier harmonics, $Z_{j}$, of interfaces; linearly perturbed solution} \\
  HWRITERA( grpOutput, mn, (Nvol+1), dZbc, dZbc(1:mn,0:Nvol) )

  endif

!latex \type{lmns}                   & integer & \pb{resolution of straight fieldline transformation} \\
  HWRITEIV( grpOutput, 1, lmns, (/ lmns /))

!latex \hline \end{tabular}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item All quantities marked as real should be treated as double precision.

!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  !RETURN(hdfint)
  HCLOSEGRP( grpOutput )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine hdfint

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
subroutine finish_outfile

  LOCALS

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  call h5fclose_f( file_id, hdfier ) ! terminate access on output file;
  FATAL( sphdf5, hdfier.ne.0, error calling h5fclose_f )

  call h5close_f( hdfier ) ! close Fortran interface to the HDF5 library;
  FATAL( sphdf5, hdfier.ne.0, error calling h5close_f )
end subroutine finish_outfile

end module sphdf5
