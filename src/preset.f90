!> \defgroup grp_initialization Initialization of the code
!>
!> \latexonly
!> \definecolor{Orange}{rgb}{1.0,0.5,0.0}
!> \definecolor{Cerulean}{rgb}{0.0,0.5,1.0}
!> \endlatexonly
!>
!> \file
!> \brief Allocates and initializes internal arrays.

!> \brief Allocates and initializes internal arrays.
!> \ingroup grp_initialization
!>
subroutine preset

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  use constants, only : zero, one, mu0

  use numerical, only : sqrtmachprec, vsmall, small

  use fileunits, only : ounit

  use inputlist

  use cputiming, only : Tpreset

  use allglobal

  use fftw_interface

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS

  INTEGER   :: innout, idof, jk, ll, ii, ifail, ideriv, vvol, mi, ni, mj, nj, mk, nk, mimj, ninj, mkmj, nknj, jj, kk, lvol, mm, nn, imn
  INTEGER   :: lquad, igauleg, maxIquad, Mrad, jquad, Lcurvature, zerdof, iret, work1, work2
  REAL      :: teta, zeta, arg, lss, cszeta(0:1), error
  LOGICAL   :: LComputeAxis

  LOGICAL              :: Lchangeangle
  INTEGER              :: nb, ix, ij, ip, idx_mode
  REAL                 :: xx


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  BEGIN(preset)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! the following was in global:readin previously
! set internal parameters that depend on physicslist;


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

    if( Lchangeangle ) then ; jj = -1  ! change sign of poloidal angle; Loizu Nov 18;
    else                    ; jj = +1
    endif

    do idx_mode=1, num_modes! will read in Fourier harmonics until the end of file is reached;
     mm = mmRZRZ(idx_mode)
     nn = nnRZRZ(idx_mode)

     do ii = 1, mn ; mi = im(ii) ; ni = in(ii) ! loop over harmonics within range;
      if( mm.eq.0 .and. mi.eq.0 .and. nn*Nfp.eq.ni ) then
       iRbc(ii,1:Nvol-1) = allRZRZ(1,1:Nvol-1, idx_mode) ! select relevant harmonics;
       iZbs(ii,1:Nvol-1) = allRZRZ(2,1:Nvol-1, idx_mode) ! select relevant harmonics;
       if( NOTstellsym ) then
        iRbs(ii,1:Nvol-1) = allRZRZ(3,1:Nvol-1, idx_mode) ! select relevant harmonics;
        iZbc(ii,1:Nvol-1) = allRZRZ(4,1:Nvol-1, idx_mode) ! select relevant harmonics;
       else
        iRbs(ii,1:Nvol-1) = zero             ! select relevant harmonics;
        iZbc(ii,1:Nvol-1) = zero             ! select relevant harmonics;
       endif
      elseif( mm.eq.mi .and. nn*Nfp.eq.jj*ni ) then
       iRbc(ii,1:Nvol-1) = allRZRZ(1,1:Nvol-1, idx_mode) ! select relevant harmonics;
       iZbs(ii,1:Nvol-1) = jj*allRZRZ(2,1:Nvol-1, idx_mode) ! select relevant harmonics;
       if( NOTstellsym ) then
        iRbs(ii,1:Nvol-1) = jj*allRZRZ(3,1:Nvol-1, idx_mode) ! select relevant harmonics;
        iZbc(ii,1:Nvol-1) = allRZRZ(4,1:Nvol-1, idx_mode) ! select relevant harmonics;
       else
        iRbs(ii,1:Nvol-1) = zero             ! select relevant harmonics;
        iZbc(ii,1:Nvol-1) = zero             ! select relevant harmonics;
       endif
      endif
     enddo ! end of do ii;

    enddo ! end of do;

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




  call random_seed()
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

 !FATAL( preset, Nfp.eq.0, illegal division ) ! this was checked in global: readin; SRH: 27 Feb 18;

  pi2nfp         = pi2 / Nfp

  pi2pi2nfp      = pi2 * pi2nfp
  pi2pi2nfphalf  = pi2 * pi2nfp * half
  pi2pi2nfpquart = pi2 * pi2nfp * quart

  Mrad  = maxval( Lrad(1:Mvol) )

  if( myid.eq.0 ) write(ounit,'("preset : ",10x," : myid=",i3," ; Mrad=",i3," : Lrad=",257(i3,",",:))') myid, Mrad, Lrad(1:Mvol)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **LGdof and NGdof : number of geometrical degrees-of-freedom**
!>
!> <ul>
!> <li> \c LGdof \f$\equiv\f$ the number of degrees-of-freedom in the geometry (i.e. Fourier harmonics) of each interface </li>
!> <li> \c NGdof \f$\equiv\f$ total number of degrees-of-freedom in geometry, i.e. of all interfaces </li>
!> </ul>

!                            Rbc  Zbs    Rbs    Zbc
  select case( Igeometry )
  case( 1:2)
   if( YESstellsym ) LGdof = mn
   if( NOTstellsym ) LGdof = mn        + mn-1
  case(   3)
   if( YESstellsym ) LGdof = mn + mn-1
   if( NOTstellsym ) LGdof = mn + mn-1 + mn-1 + mn
  end select

  NGdof = ( Mvol-1 ) * LGdof

  if( Wpreset ) then ; cput = GETTIME ; write(ounit,'("preset : ",f10.2," : myid=",i3," ; NGdof=",i9," ;")') cput-cpus, myid, NGdof
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **iota and oita: rotational transform on interfaces**
!>
!> <ul>
!> <li> The input variables \c iota and \c oita are the rotational transform
!>       on "inner-side" and on the "outer-side" of each interface. </li>
!> <li> These quantities are formally inputs. </li>
!> <li> Note that if \f$q_l+\gamma q_r \ne 0\f$, then \c iota is given by
!>       \f{eqnarray}{ {{\,\iota\!\!\!}-} \equiv \frac{p_l + \gamma p_r}{q_l + \gamma q_r},
!>       \f}
!>       where \f$p_l \equiv\,\f$\c pl, \f$q_l \equiv\,\f$\c ql , etc.;
!>       and similarly for \c oita . </li>
!> </ul>

  do vvol = 0, Nvol

   if( ql(vvol).eq.0 .and. qr(vvol).eq.0 ) then ; iota(vvol) = iota(vvol)
   else                                         ; iota(vvol) = ( pl(vvol) + goldenmean * pr(vvol) ) / ( ql(vvol) + goldenmean * qr(vvol) )
   endif

   if( lq(vvol).eq.0 .and. rq(vvol).eq.0 ) then ; oita(vvol) = oita(vvol)
   else                                         ; oita(vvol) = ( lp(vvol) + goldenmean * rp(vvol) ) / ( lq(vvol) + goldenmean * rq(vvol) )
   endif

   if( Wpreset .and. myid.eq.0 ) write(ounit,1002) vvol, pl(vvol), ql(vvol), pr(vvol), qr(vvol), iota(vvol), lp(vvol), lq(vvol), rp(vvol), rq(vvol), oita(vvol)

1002 format("preset : ",10x," :      ",3x," ; transform : ",i3," : (",i3," /",i3," ) * (",i3," /",i3," ) = ",f18.15," ; ",&
                                                                  "(",i3," /",i3," ) * (",i3," /",i3," ) = ",f18.15," ; ")

  enddo

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **dtflux(1:Mvol) and dpflux(1:Mvol): enclosed fluxes**
!>
!> <ul>
!> <li> \c dtflux \f$\equiv \Delta \psi_{tor} / 2\pi\f$ and
!>      \c dpflux \f$\equiv \Delta \psi_{pol} / 2\pi\f$ in each volume. </li>
!> <li> Note that the total toroidal flux enclosed by the plasma boundary is \f$\Phi_{edge} \equiv\,\f$\c phiedge . </li>
!> <li> \f$\psi_{tor} \equiv\,\f$\c tflux and \f$\psi_{pol} \equiv\,\f$\c pflux are immediately normalized (in readin() ) according to
!>      \f$\psi_{tor,i} \rightarrow \psi_{tor,i} / \psi_{0}\f$ and
!>      \f$\psi_{pol,i} \rightarrow \psi_{pol,i} / \psi_{0}\f$, where \f$\psi_{0} \equiv \psi_{tor,N}\f$ on input. </li>
!> </ul>

  SALLOCATE( dtflux, (1:Mvol), zero )
  SALLOCATE( dpflux, (1:Mvol), zero )

  if(Lconstraint.eq.3 .and. Igeometry.eq.1) then
      total_pflux = pflux(Mvol) * phiedge / pi2
      pflux(Mvol) = 0
  endif

  select case( Igeometry )
  case( 1   ) ; dtflux(1) = tflux(1) ; dpflux(1) = pflux(1) ! Edit by Erol, to keep total pflux 0
  case( 2:3 ) ; dtflux(1) = tflux(1) ; dpflux(1) = zero     ! cylindrical or toroidal;
  end select

  dtflux(2:Mvol) = tflux(2:Mvol) - tflux(1:Mvol-1)
  dpflux(2:Mvol) = pflux(2:Mvol) - pflux(1:Mvol-1)

  dtflux(1:Mvol) = dtflux(1:Mvol) * phiedge / pi2
  dpflux(1:Mvol) = dpflux(1:Mvol) * phiedge / pi2


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{\type{mu(1:Mvol)} evaluation of \inputvar{mu} from \inputvar{Ivolume};}
!latex Only used when $Lconstraint = 3$. The coefficients \inputvar{mu} are evaluated as
!latex \begin{equation}
!latex \mu_1 = \mu_0 \frac{I_{volume}(1)}{\psi_{t,1}}
!latex \end{equation}
!latex \begin{equation}
!latex \mu_n = \mu_0 \frac{I_{volume}(n) - I_{volume}(n-1)}{\psi_{t,n}-\psi_{t,n-1}},\qquad \forall\ n>1;
!latex \end{equation}

if (Lconstraint.EQ.3) then

  mu(1) = Ivolume(1) / (tflux(1) * phiedge)

  do vvol = 2, Mvol
    mu(vvol) = (Ivolume(vvol) - Ivolume(vvol-1)) / ((tflux(vvol) - tflux(vvol-1)) * phiedge)
  end do

#ifdef DEBUG
  if (myid.eq.0) then
    write(*,*) " "
    write(ounit,'("preset : ", 10x ," : Ivolume = "257(es11.3",",:))') (Ivolume(vvol), vvol=1, Mvol)
    write(ounit,'("preset : ", 10x ," : tflux   = "257(es11.3",",:))') (  tflux(vvol), vvol=1, Mvol)
    write(ounit,'("preset : ", 10x ," : phiedge = "257(es11.3",",:))') phiedge
    write(ounit,'("preset : ", 10x ," : mu      = "257(es11.3",",:))') (     mu(vvol), vvol=1, Mvol)
  end if
#endif
endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **sweight(1:Mvol): star-like angle constraint weight**
!>
!> <ul>
!> <li> the "star-like" poloidal angle constraint weights (only required for toroidal geometry, i.e. \c Igeometry=3) are given by
!>       \f{eqnarray}{ \texttt{sweight}_v \equiv \texttt{upsilon} \times (l_v / N_{vol})^w,
!>       \f}
!>       where \f$l_v\f$ is the volume number,
!>       and \f$w \equiv\,\f$\c wpoloidal. </li>
!> </ul>

  SALLOCATE( sweight, (1:Mvol), zero )
 !sweight(1:Mvol) = upsilon * tflux(1:Mvol)**wpoloidal ! toroidal flux in vacuum region is not constant; 11 July 18;
  do vvol = 1, Mvol ; sweight(vvol) = upsilon * (vvol*one/Nvol)**wpoloidal ! 11 July 18;
  enddo

#ifdef DEBUG
  if (myid.eq.0) then
    write(ounit,'("preset : ",10x," : sweight =",99(es12.5,",",:))') sweight(1:Mvol)
  end if
#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **TT(0:Mrad,0:1,0:1): Chebyshev polynomials at inner/outer interface**
!>
!> <ul>
!> <li> \c TT(0:Lrad,0:1,0:1) gives the Chebyshev polynomials, and their first derivative, evaluated at \f$s=-1\f$ and \f$s=+1\f$. </li>
!> <li> Precisely, \c TT(l,i,d) \f$\equiv T_l^{(d)}(s_i)\f$ for \f$s_0=-1\f$ and \f$s_1=+1\f$. </li>
!> <li> Note that \f$T_l^{(0)}(s)=s^l\f$ and \f$T_l^{(1)}(s)=s^{l+1} l^2\f$ for \f$s=\pm 1\f$. </li>
!> <li> Note that
!>       \f{eqnarray}{ T_l(-1)        = \left\{ \begin{array}{ccccccccccccccc}+1,& \textrm{ if $l$ is even,} \\
!>                                                                            -1,& \textrm{ if $l$ is odd;}
!>                                              \end{array} \right. & \; \;&
!>                     T_l(+1)        = \left\{ \begin{array}{ccccccccccccccc}+1,& \textrm{ if $l$ is even,} \\
!>                                                                            +1,& \textrm{ if $l$ is odd;}
!>                                              \end{array} \right. \\
!>                     T_l^\prime(-1) = \left\{ \begin{array}{ccccccccccccccc}-l^2,& \textrm{ if $l$ is even,} \\
!>                                                                            +l^2,& \textrm{ if $l$ is odd;}
!>                                              \end{array} \right. &\; \;&
!>                     T_l^\prime(+1) = \left\{ \begin{array}{ccccccccccccccc}+l^2,& \textrm{ if $l$ is even,} \\
!>                                                                            +l^2,& \textrm{ if $l$ is odd.}
!>                                               \end{array} \right.
!>       \f} </li>
!> <li> \c TT(0:Mrad,0:1,0:1) is used in routines that explicity require interface information, such as
!>       <ul>
!>         <li> the interface force-balance routine,                                lforce() </li>
!>         <li> the virtual casing routine,                                         casing() </li>
!>         <li> computing the rotational-transform on the interfaces,               tr00ab() </li>
!>         <li> computing the covariant components of the interface magnetic field, sc00aa() </li>
!>         <li> enforcing the constraints on the Beltrami fields,                   matrix() </li>
!>     and <li> computing the enclosed currents of the vacuum field,                curent(). </li>
!>       </ul> </li>
!> </ul>

  SALLOCATE( TT, (0:Mrad,0:1,0:1), zero )
  SALLOCATE(RTT, (0:Lrad(1),0:Mpol,0:1,0:1), zero )
  SALLOCATE(RTM, (0:Lrad(1),0:Mpol), zero )

  call get_cheby( -one, Mrad, TT(:,0,:))
  call get_cheby( one , Mrad, TT(:,1,:))

  call get_zernike( zero, Lrad(1), Mpol, RTT(:,:,0,:))
  call get_zernike( one, Lrad(1), Mpol, RTT(:,:,1,:))
  call get_zernike_rm(zero, Lrad(1), Mpol, RTM(:,:))

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **ImagneticOK(1:Mvol): Beltrami/vacuum error flag**
!>
!> <ul>
!> <li> error flags that indicate if the magnetic field in each volume has been successfully constructed </li>
!> <li> \c ImagneticOK is initialized to \c .false. in dforce() before the Beltrami solver routines are called.
!>       If the construction of the Beltrami field is successful
!>       (in either ma02aa() or mp00ac() )
!>       then \c ImagneticOK is set to \c .true. . </li>
!> </ul>

  SALLOCATE( ImagneticOK, (1:Mvol), .false. )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsubsection{\type{IconstraintOK} : Global constraint flag;}

!latex \begin{enumerate}
!latex \item error flags that indicate if solution found is consistent with the global constraint
!latex \item \type{ImagneticOK} is initialized to \type{.false.}
!latex       If the construction the difference between the evaluated global constraint and in the input is small enough,
!latex       then \type{ImagneticOK} is set to \type{.true.}.
!latex \end{enumerate}

!> **Lhessianallocated**
!>
!> <ul>
!> <li> The internal logical variable, \c Lhessianallocated, indicates whether the ``Hessian'' matrix of second-partial derivatives
!>       (really, the first derivatives of the force-vector) has been allocated, or not! </li>
!> </ul>
  Lhessianallocated = .false.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **ki(1:mn,0:1): Fourier identification**
!>
!> <ul>
!> <li> Consider the "abbreviated" representation for a double Fourier series,
!>       \f{eqnarray}{ \sum_i      f_i \cos(     m_i \theta -      n_i \zeta) \equiv                         \sum_{n=      0  }^{     N_0} f_{0,n} \cos(       -n\zeta)
!>                                                                                   + \sum_{m=1}^{     M_0} \sum_{n=-     N_0}^{     N_0} f_{m,n} \cos(m\theta-n\zeta),
!>       \f}
!>       and the same representation but with enhanced resolution,
!>       \f{eqnarray}{ \sum_k \bar f_k \cos(\bar m_k \theta - \bar n_k \zeta) \equiv                         \sum_{n=      0  }^{     N_1} f_{0,n} \cos(       -n\zeta)
!>                                                                                   + \sum_{m=1}^{     M_1} \sum_{n=-     N_1}^{     N_1} f_{m,n} \cos(m\theta-n\zeta),
!>       \label{eq:enhancedFourierrepresentation_preset}
!>       \f}
!>       with \f$M_1 \ge M_0\f$ and \f$N_1 \ge N_0\f$;
!>       then \f$k_i\equiv\,\f$\c ki(i,0) is defined such that \f$\bar m_{k_i} = m_i\f$ and \f$\bar n_{k_i} = n_i\f$. </li>
!> </ul>

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **kija(1:mn,1:mn,0:1), kijs(1:mn,1:mn,0:1): Fourier identification**
!>
!> <ul>
!> <li> Consider the following quantities, which are computed in ma00aa(),
!>       where \f$\bar g^{\mu\nu} = \sum_k \bar g^{\mu\nu}_k \cos \alpha_k\f$ for \f$\alpha_k \equiv m_k \theta - n_k \zeta\f$,
!>       \f{eqnarray}{
!>       \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \bar g^{\mu\nu} \cos\alpha_i \; \cos\alpha_j & = & \frac{1}{2} \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \bar g^{\mu\nu} ( + \cos\alpha_{k_{ij+}} + \cos\alpha_{k_{ij-}} ), \\
!>       \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \bar g^{\mu\nu} \cos\alpha_i \; \sin\alpha_j & = & \frac{1}{2} \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \bar g^{\mu\nu} ( + \sin\alpha_{k_{ij+}} - \sin\alpha_{k_{ij-}} ), \\
!>       \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \bar g^{\mu\nu} \sin\alpha_i \; \cos\alpha_j & = & \frac{1}{2} \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \bar g^{\mu\nu} ( + \sin\alpha_{k_{ij+}} + \sin\alpha_{k_{ij-}} ), \\
!>       \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \bar g^{\mu\nu} \sin\alpha_i \; \sin\alpha_j & = & \frac{1}{2} \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \bar g^{\mu\nu} ( - \cos\alpha_{k_{ij+}} + \cos\alpha_{k_{ij-}} ),
!>       \f}
!>       where \f$(m_{k_{ij+}}, n_{k_{ij+}}) = (m_i + m_j, n_i + n_j)\f$ and \f$(m_{k_{ij-}}, n_{k_{ij-}}) = (m_i - m_j, n_i - n_j)\f$;
!>       then \c kija(i,j,0) \f$\equiv k_{ij+}\f$ and \c kijs(i,j,0) \f$\equiv k_{ij-}\f$. </li>
!> <li> Note that Eqn.\f$(\ref{eq:enhancedFourierrepresentation_preset})\f$ does not include \f$m<0\f$; so,
!>       if \f$m_i - m_j < 0\f$ then \f$k_{ij-}\f$ is re-defined such that \f$(m_{k_{ij-}}, n_{k_{ij-}}) = (m_j - m_i, n_j - n_i)\f$; and
!>       similarly for the case \f$m=0\f$ and \f$n<0\f$.
!>       Also, take care that the sign of the sine harmonics in the above expressions will change for these cases. </li>
!> </ul>

  SALLOCATE( ki, (1:mn,0:1), 0 )
  SALLOCATE( kija, (1:mn,1:mn,0:1), 0 )
  SALLOCATE( kijs, (1:mn,1:mn,0:1), 0 )

  do ii = 1, mn  ; mi =  im(ii) ; ni =  in(ii)

    call getimn(lMpol, lNtor, Nfp, mi, ni, kk)
    if (kk.gt.0) then
      if( mi.eq.0 .and. ni.eq.0 ) then ; ki(ii,0:1) = (/ kk, 1 /)
      else                             ; ki(ii,0:1) = (/ kk, 2 /)
      endif
    endif

    do jj = 1, mn  ; mj =  im(jj) ; nj =  in(jj) ; mimj = mi + mj ; ninj = ni + nj !   adding   ; 17 Dec 15;

      call getimn(lMpol, lNtor, Nfp, mimj, ninj, kk)
      if (kk.gt.0) then
        if( mimj.eq.0 .and. ninj.eq.0 ) then ; kija(ii,jj,0:1) = (/ kk, 1 /)
        else                                 ; kija(ii,jj,0:1) = (/ kk, 2 /)
        endif
      endif
      ;                                           ; mimj = mi - mj ; ninj = ni - nj ! subtracting; 17 Dec 15;

      if( mimj.gt.0 .or. ( mimj.eq.0 .and. ninj.ge.0 ) ) then
        call getimn(lMpol, lNtor, Nfp, mimj, ninj, kk)
        if (kk.gt.0) then
          if( mimj.eq.0 .and. ninj.eq.0 ) then ; kijs(ii,jj,0:1) = (/ kk, 1 /)
          else                                 ; kijs(ii,jj,0:1) = (/ kk, 2 /)
          endif
        endif
      else
        call getimn(lMpol, lNtor, Nfp, -mimj, -ninj, kk)
        if (kk.gt.0) then
          ;                                    ; kijs(ii,jj,0:1) = (/ kk , - 2 /) ! only the sine modes need the sign factor; 17 Dec 15;
        endif
      endif

    enddo ! end of do jj; 29 Jan 13;

  enddo ! end of do ii; 29 Jan 13;
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **djkp**

  if( Igeometry.eq.2 ) then ! standard cylindrical; 04 Dec 14;

   SALLOCATE( djkp, (1:mn,1:mn), 0 ) ! only used in volume; trigonometric identities; 04 Dec 14;
   SALLOCATE( djkm, (1:mn,1:mn), 0 ) ! only used in volume; trigonometric identities; 04 Dec 14;

   do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
    do jj = 1, mn ; mj = im(jj) ; nj = in(jj)
     if( mi-mj.eq.0 .and. ni-nj.eq.0 ) djkp(ii,jj) = 1
     if( mi+mj.eq.0 .and. ni+nj.eq.0 ) djkm(ii,jj) = 1
    enddo
   enddo

  endif ! end of if( Igeometry.eq.2 ) ; 04 Dec 14;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **iotakki**

  SALLOCATE( iotakkii, (1:mn      ), 0 ) ! used to identify matrix elements in straight-field-line angle transformation;

  SALLOCATE( iotaksub, (1:mn,1:mns), 0 )
  SALLOCATE( iotaksgn, (1:mn,1:mns), 0 )
  SALLOCATE( iotakadd, (1:mn,1:mns), 0 )

  do kk = 1, mn ; mk = im(kk) ; nk = in(kk)

    call getimn(sMpol, sNtor, Nfp, mk, nk, ii)
    if (ii.gt.0) iotakkii(kk) = ii

    do jj = 1, mns ; mj = ims(jj) ; nj = ins(jj)

      mkmj = mk - mj ; nknj = nk - nj

      if( mkmj.gt.0 .or. ( mkmj.eq.0 .and. nknj.ge.0 ) ) then

        call getimn(sMpol, sNtor, Nfp, mkmj, nknj, ii)
        if (ii.gt.0) then ; iotaksub(kk,jj) = ii ; iotaksgn(kk,jj) =  1
        endif

      else

        call getimn(sMpol, sNtor, Nfp, -mkmj, -nknj, ii)
        if (ii.gt.0) then ; iotaksub(kk,jj) = ii ; iotaksgn(kk,jj) =  -1
        endif

      endif

      mkmj = mk + mj ; nknj = nk + nj

      call getimn(sMpol, sNtor, Nfp, mkmj, nknj, ii)
      if (ii.gt.0) then ; iotakadd(kk,jj) = ii
      endif

    enddo ! end of do jj; 29 Jan 13;

  enddo ! end of do kk; 29 Jan 13;


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **cheby(0:Lrad,0:2): Chebyshev polynomial workspace**
!>
!> <ul>
!> <li> \c cheby(0:Lrad,0:2) is global workspace for computing the Chebyshev polynomials, and their derivatives,
!>       using the recurrence relations \f$T_0(s) = 1\f$, \f$T_1(s) = s\f$ and  \f$T_l(s) = 2 \, s \,T_{l-1}(s) - T_{l-2}(s)\f$. </li>
!> <li> These are computed as required, i.e. for arbitrary \f$s\f$, in bfield(), jo00aa() and ma00aa(). </li>
!> <li> Note that the quantities required for ma00aa() are for fixed \f$s\f$, and so these quantities should be precomputed. </li>
!> </ul>

! Allocate space for the toroidal current array in each interface

  SALLOCATE( IPDt, (1:Mvol), zero)
  if( Lfreebound.eq.1 ) then
    SALLOCATE( IPDtDpf, (1:Mvol  , 1:Mvol  ), zero)
  else
    if(Igeometry.eq.1) then
      ! add an additional constraint to make the total pflux = 0 -- Edit Erol
      SALLOCATE( IPDtDpf, (1:Mvol, 1:Mvol), zero) 
    else
      SALLOCATE( IPDtDpf, (1:Mvol-1, 1:Mvol-1), zero)
    endif
  endif

  SALLOCATE( cheby, (0:Mrad,0:2), zero )
  SALLOCATE( zernike, (0:Lrad(1), 0:Mpol, 0:2), zero )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **Iquad, gaussianweight, gaussianabscissae: Gauss-Legendre quadrature**
!>
!> <ul>
!> <li> The volume integrals are computed using a "Fourier" integration over the angles and by Gauss-Legendre quadrature over the radial,
!>       i.e. \f$\displaystyle \int \!\! f(s) ds = \sum_k \omega_k f(s_k)\f$. </li>
!> <li> The quadrature resolution in each volume is give by \c Iquad(1:Mvol) which is determined as follows:
!> <ul>
!> <li> if \c Nquad.gt.0                                    , then \c Iquad(vvol)=Nquad                   </li>
!> <li> if \c Nquad.le.0 and \c .not.Lcoordinatesingularity , then \c Iquad(vvol)=2*Lrad(vvol)-Nquad      </li>
!> <li> if \c Nquad.le.0 and      \c Lcoordinatesingularity , then \c Iquad(vvol)=2*Lrad(vvol)-Nquad+Mpol </li>
!> </ul> </li>
!> <li> The Gaussian weights and abscissae are given by \c gaussianweight(1:maxIquad,1:Mvol) and \c gaussianabscissae(1:maxIquad,1:Mvol),
!>       which are computed using modified Numerical Recipes routine gauleg() . </li>
!> <li> \c Iquad\f$_v\f$ is passed through to ma00aa() to compute the volume integrals of the metric elements;
!>       also see jo00aa(), where \c Iquad\f$_v\f$ is used to compute the volume integrals of \f$||\nabla\times{\bf B} - \mu {\bf B}||\f$. </li>
!> </ul>

  SALLOCATE( Iquad, (1:Mvol), 0 ) ! 16 Jan 13;

  do vvol = 1, Mvol

   LREGION(vvol)

   if( Nquad.gt.0 ) then ;            Iquad(vvol) =                         Nquad
   else
    if(      Lcoordinatesingularity ) Iquad(vvol) = Mpol + 2 * Lrad(vvol) - Nquad ! NEED TO REVISE REGULARIZATION FACTORS; 26 Feb 13;
    if( .not.Lcoordinatesingularity ) Iquad(vvol) =        2 * Lrad(vvol) - Nquad
   endif

  enddo ! end of do vvol; 18 Feb 13;

  maxIquad = maxval(Iquad(1:Mvol))

  SALLOCATE( gaussianweight   , (1:maxIquad,1:Mvol), zero ) ! perhaps it would be neater to make this a structure; 26 Jan 16;
  SALLOCATE( gaussianabscissae, (1:maxIquad,1:Mvol), zero )

  do vvol = 1, Mvol

   lquad = Iquad(vvol)

   call gauleg( lquad, gaussianweight(1:lquad,vvol), gaussianabscissae(1:lquad,vvol), igauleg ) ! JAB; 28 Jul 17

   if( myid.eq.0 ) then
    cput= GETTIME
    select case( igauleg ) !                                                  123456789012345
    case( 0 )    ; if( Wpreset ) write(ounit,1000) cput-cpus, vvol, igauleg, "success        ", gaussianabscissae(1:lquad,vvol)
    case( 1 )    ;               write(ounit,1000) cput-cpus, vvol, igauleg, "failed         ", gaussianabscissae(1:lquad,vvol)
    case( 2 )    ;               write(ounit,1000) cput-cpus, vvol, igauleg, "input error    ", gaussianabscissae(1:lquad,vvol)
    case default ;               write(ounit,1000) cput-cpus, vvol, igauleg, "weird          ", gaussianabscissae(1:lquad,vvol)
     FATAL( preset, .true., weird ifail returned by gauleg )
    end select
    ;            ; if( Wpreset ) write(ounit,1001)                                              gaussianweight(1:lquad,vvol)
   endif

1000 format("preset : ",f10.2," : lvol=",i3," ; igauleg=",i5," ; ",a15," ; abscissae ="99f09.05)
1001 format("preset : ", 10x ," :      ",3x,"           ",5x,"   ",15x," ; weights   ="99f09.05)

  enddo ! end of do vvol;  7 Mar 13;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **LBsequad, LBnewton and LBlinear**
!>
!> <ul>
!> <li> \c LBsequad, \c LBnewton and \c LBlinear depend simply on \c LBeltrami , which is described in global.f90 . </li>
!> </ul>

  LBsequad = .false.
  LBnewton = .false.
  LBlinear = .false.

  if( LBeltrami.eq.1 .or. LBeltrami.eq.3 .or. LBeltrami.eq.5 .or. LBeltrami.eq.7 ) LBsequad = .true.
  if( LBeltrami.eq.2 .or. LBeltrami.eq.3 .or. LBeltrami.eq.6 .or. LBeltrami.eq.7 ) LBnewton = .true.
  if( LBeltrami.eq.4 .or. LBeltrami.eq.5 .or. LBeltrami.eq.6 .or. LBeltrami.eq.7 ) LBlinear = .true.

  if (LBnewton .or. LBsequad) Lconstraint = 2

  if (Lconstraint .eq. 2) then
    FATAL( preset, Lfreebound.eq.1, The combination of helicity constraint and free boundary is under construction )
    if (Igeometry .eq. 3 .and. myid.eq.0) then
      write(ounit, *) 'WARNING: The Hessian matrix needs further review for Igeometry = 3'
      write(ounit, *) '         However, it can still serve the purpose of Lfindzero = 2'
    endif
  endif

  if( myid.eq.0 ) then
   cput = GETTIME
   write(ounit,'("preset : ",f10.2," : LBsequad="L2" , LBnewton="L2" , LBlinear="L2" ;")')cput-cpus, LBsequad, LBnewton, LBlinear
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **BBweight(1:mn): weighting of force-imbalance harmonics**
!>
!> <ul>
!> <li> weight on force-imbalance harmonics;
!>       \f{eqnarray}{ \texttt{BBweight}_i \equiv \texttt{opsilon} \times \exp\left[ - \texttt{escale} \times (m_i^2 + n_i^2) \right]
!>       \f} </li>
!> <li> this is only used in dforce() in constructing the force-imbalance vector </li>
!> </ul>

  SALLOCATE( BBweight, (1:mn), opsilon * exp( - escale * ( im(1:mn)**2 + (in(1:mn)/Nfp)**2 ) ) )

  if( myid.eq.0 .and. escale.gt.small ) then
   do ii = 1, mn ; write(ounit,'("preset : " 10x " : myid="i3" ; ("i3","i3") : BBweight="es13.5" ;")') myid, im(ii), in(ii)/Nfp, BBweight(ii)
   enddo
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **mmpp(1:mn): spectral condensation weight factors**
!>
!> <ul>
!> <li> spectral condensation weight factors;
!>       \f{eqnarray}{ \texttt{mmpp(i)} \equiv m_i^p,
!>       \f}
!>       where \f$p \equiv\,\f$\c pcondense . </li>
!> </ul>

  SALLOCATE( mmpp, (1:mn), zero )

  do ii = 1, mn ; mi = im(ii)

   if( mi.eq.0 ) then ; mmpp(ii) = zero
   else               ; mmpp(ii) = mi**pcondense
   endif ! end of if( mi.eq.0 ) ; 11 Aug 14;

  enddo ! end of do ii; 08 Nov 13;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **NAdof, Ate, Aze, Ato and Azo: degrees-of-freedom in magnetic vector potential**
!>
!> <ul>
!> <li> \c NAdof(1:Mvol) \f$\equiv\f$ total number of degrees-of-freedom in magnetic vector potential, including Lagrange multipliers, in each volume.
!>       This can de deduced from matrix(). </li>
!> <li> The components of the vector potential, \f${\bf A}=A_\theta \nabla + A_\zeta \nabla \zeta\f$, are
!>      \f{eqnarray}{
!>        A_\theta(s,\theta,\zeta) &=& \sum_{i,l} {\color{red}  A_{\theta,e,i,l}} \; {\overline T}_{l,i}(s) \cos\alpha_i + \sum_{i,l} {\color{Orange}  A_{\theta,o,i,l}} \; {\overline T}_{l,i}(s) \sin\alpha_i, \label{eq:At_preset} \\
!>        A_\zeta( s,\theta,\zeta) &=& \sum_{i,l} {\color{blue} A_{\zeta, e,i,l}} \; {\overline T}_{l,i}(s) \cos\alpha_i + \sum_{i,l} {\color{Cerulean}A_{\zeta ,o,i,l}} \; {\overline T}_{l,i}(s) \sin\alpha_i, \label{eq:Az_preset}
!>      \f}
!>      where \f${\overline T}_{l,i}(s) \equiv \bar s^{m_i/2} \, T_l(s)\f$, \f$T_l(s)\f$ is the Chebyshev polynomial, and \f$\alpha_j \equiv m_j\theta-n_j\zeta\f$.
!>      The regularity factor, \f$\bar s^{m_i/2}\f$, where \f$\bar s \equiv (1+s)/2\f$, is only included if there is a coordinate singularity in the domain
!>      (i.e. only in the innermost volume, and only in cylindrical and toroidal geometry.) </li>
!> <li> The Chebyshev-Fourier harmonics of the covariant components of the magnetic vector potential are kept in
!>      \f{eqnarray}{
!>          {\color{red}     A_{\theta,e,i,l}} &\equiv& \texttt{Ate(v,0,j)}\%\texttt{s(l)} , \\
!>          {\color{blue}    A_{\zeta, e,i,l}} &\equiv& \texttt{Aze(v,0,j)}\%\texttt{s(l)} , \\
!>          {\color{Orange}  A_{\theta,o,i,l}} &\equiv& \texttt{Ato(v,0,j)}\%\texttt{s(l)} , \mathrm{and} \\
!>          {\color{Cerulean}A_{\zeta ,o,i,l}} &\equiv& \texttt{Azo(v,0,j)}\%\texttt{s(l)} ;
!>      \f}
!>      where \f$v=1,\texttt{Mvol}\f$ labels volume, \f$j=1,\texttt{mn}\f$ labels Fourier harmonic, and \f$l=0,\,\f$\c Lrad \f$(v)\f$ labels Chebyshev polynomial.
!>      (These arrays also contains derivative information.) </li>
!> <li> If \c Linitguess=1 , a guess for the initial state for the Beltrami fields is constructed.
!>      An initial state is required for iterative solvers of the Beltrami fields, see \c LBeltrami . </li>
!> <li> If \c Linitguess=2 , the initial state for the Beltrami fields is read from file (see ra00aa() ).
!>      An initial state is required for iterative solvers of the Beltrami fields, see \c LBeltrami . </li>
!> </ul>

  SALLOCATE( NAdof, (1:Mvol          ), 0 ) ! Beltrami degrees-of-freedom in each annulus;
  SALLOCATE( Nfielddof,(1:Mvol       ), 0 ) ! Beltrami degrees-of-freedom in each annulus, field only;
  SALLOCATE( NdMASmax, (1:Mvol       ), 0 ) ! The maximum size of sparse matrix for GMRES preconditioning;
  SALLOCATE( NdMAS   , (1:Mvol       ), 0 ) ! The actual size of sparse matrix for GMRES preconditioning;

  NALLOCATE( Ate  , (1:Mvol,-2:2,1:mn)    ) ! recall that this is type:sub-grid; 31 Jan 13;
  NALLOCATE( Aze  , (1:Mvol,-2:2,1:mn)    ) ! -2 : for use of matrix-free solver ; -1 : for use of force gradient
  NALLOCATE( Ato  , (1:Mvol,-2:2,1:mn)    ) !  0 : normal data
  NALLOCATE( Azo  , (1:Mvol,-2:2,1:mn)    ) ! 1:2: use to compute derivative w.r.t. fluxes

  SALLOCATE( Fso  , (1:Mvol,     1:mn), 0 ) ! these will become redundant if/when Lagrange multipliers are used to enforce bounday constraints; 26 Jan 16;
  SALLOCATE( Fse  , (1:Mvol,     1:mn), 0 )

  SALLOCATE( Lma  , (1:Mvol,     1:mn), 0 ) ! degree of freedom index; for Lagrange multiplier; 08 Feb 16;
  SALLOCATE( Lmb  , (1:Mvol,     1:mn), 0 )
  SALLOCATE( Lmc  , (1:Mvol,     1:mn), 0 ) ! only need Lmc(2:mn) ; only for NOTstellsym; 08 Feb 16;
  SALLOCATE( Lmd  , (1:Mvol,     1:mn), 0 ) ! only need Lmd(2:mn) ; only for NOTstellsym; 08 Feb 16;
  SALLOCATE( Lme  , (1:Mvol,     1:mn), 0 ) ! only need Lme(2:mn) ;
  SALLOCATE( Lmf  , (1:Mvol,     1:mn), 0 ) ! only need Lmf(2:mn) ; only for NOTstellsym; 08 Feb 16;
  SALLOCATE( Lmg  , (1:Mvol,     1:mn), 0 ) ! only need Lmg(1   ) ;
  SALLOCATE( Lmh  , (1:Mvol,     1:mn), 0 ) ! only need Lmh(1   ) ;

  SALLOCATE( Lmavalue, (1:Mvol,     1:mn), zero )
  SALLOCATE( Lmbvalue, (1:Mvol,     1:mn), zero )
  SALLOCATE( Lmcvalue, (1:Mvol,     1:mn), zero )
  SALLOCATE( Lmdvalue, (1:Mvol,     1:mn), zero )
  SALLOCATE( Lmevalue, (1:Mvol,     1:mn), zero )
  SALLOCATE( Lmfvalue, (1:Mvol,     1:mn), zero )
  SALLOCATE( Lmgvalue, (1:Mvol,     1:mn), zero )
  SALLOCATE( Lmhvalue, (1:Mvol,     1:mn), zero )

  do vvol = 1, Mvol

   LREGION(vvol)

   if( Lcoordinatesingularity ) then
    zerdof = 0                                       ! count Zernike degree of freedom 30 Jun 19
    do ii = 2, Mpol                                  ! for m>1
     do jj = ii, Lrad(vvol), 2
      zerdof = zerdof + 2 * ntor + 1                 ! plus and minus sign for n>1, unique for n==0
      if( NOTstellsym ) zerdof = zerdof + 2*ntor + 1 ! plus and minus sign for n
     enddo
    enddo
    zerdof = zerdof * 2                              ! we have one for At and one for Az

    do jj = 0, Lrad(vvol), 2                         ! for m==0
     zerdof = zerdof + ntor + 1                      ! minus sign for n, Aze
     if (jj .ge. 2) zerdof = zerdof + ntor + 1       ! minus sign for n, Ate, without l=0 due to recombination

     if( NOTstellsym ) then
      zerdof = zerdof + ntor                         ! sin component minus sign for n, Azo
      if (jj .ge. 2) zerdof = zerdof + ntor          ! minus sign for n, Ato, without l=0 due to recombination
     endif
    enddo

    if (Mpol .ge. 1) then ! for m==1
      do jj = 1, Lrad(vvol), 2
        zerdof = zerdof + 2 * ntor + 1                  ! minus and plus sign for n, Aze
        if (jj .ge. 2) zerdof = zerdof + 2 * ntor + 1   ! minus sign for n, Ate, without l=0 due to recombination

        if( NOTstellsym ) then
          zerdof = zerdof + 2 * ntor + 1                 ! sin component minus and plus sign for n, Azo
          if (jj .ge. 2) zerdof = zerdof + 2 * ntor + 1  ! minus and plus sign for n, Ato, without l=0 due to recombination
        endif
      enddo
    endif

    ! the degree of freedom in the Beltrami field without Lagrange multipliers
    Nfielddof(vvol) = zerdof
                                     !                                     a    c      b        d      e      f      g   h
    if( YESstellsym ) NAdof(vvol) = zerdof                               + mn        + Ntor+1        + mn-1        + 1 + 0
    if( NOTstellsym ) NAdof(vvol) = zerdof                               + mn + mn-1 + Ntor+1 + Ntor + mn-1 + mn-1 + 1 + 0 ! this is broken at the moment

    ! due to basis recombination, Lma will not have the m=0 and m=1 harmonics. We substract them now
    ! m = 0
    NAdof(vvol) = NAdof(vvol) - (ntor + 1)
    if (NOTstellsym) NAdof(vvol) = NAdof(vvol) - ntor

    ! m = 1
    if (Mpol .ge. 1) then
      NAdof(vvol) = NAdof(vvol) - (2 * ntor + 1)
      if (NOTstellsym) NAdof(vvol) = NAdof(vvol) - (2 * ntor + 1)
    endif

    ! Guess the size of the sparse matrix ! 28 Jan 20
    ! If an iterative method is used and requires an preconditioner, we need to construct it as a sparse matrix
    if (Lmatsolver.ge.2 .and. LGMRESprec.gt.0) then
      if( YESstellsym ) NdMASmax(vvol) = (2 * (Lrad(vvol)/2 + 1))**2 * mn + 2 * 2 * 5 * Lrad(vvol) * mn ! Ate, Aze
      if( NOTstellsym ) NdMASmax(vvol) = (4 * (Lrad(vvol)/2 + 1))**2 * mn + 2 * 4 * 8 * Lrad(vvol) * mn ! Ate, Aze, Ato, Azo
    end if
   else ! .not.Lcoordinatesingularity;                                     a    c      b        d      e      f      g   h
    if( YESstellsym ) NAdof(vvol) = 2 * ( mn        ) * ( Lrad(vvol)    )                            + mn-1        + 1 + 1
    if( NOTstellsym ) NAdof(vvol) = 2 * ( mn + mn-1 ) * ( Lrad(vvol)    )                            + mn-1 + mn-1 + 1 + 1

    ! dof for field variables only
    if( YESstellsym ) Nfielddof(vvol) = 2 * ( mn        ) * ( Lrad(vvol)    )
    if( NOTstellsym ) Nfielddof(vvol) = 2 * ( mn + mn-1 ) * ( Lrad(vvol)    )

    ! Guess the size of the sparse matrix ! 28 Jan 20
    ! If an iterative method is used and requires an preconditioner, we need to construct it as a sparse matrix
    if (Lmatsolver.ge.2 .and. LGMRESprec.gt.0) then
      if( YESstellsym ) NdMASmax(vvol) = (2 * (Lrad(vvol) + 1))**2 * mn + 2 * 2 * 5 * Lrad(vvol) * mn        ! Ate, Aze
      if( NOTstellsym ) NdMASmax(vvol) = (4 * (Lrad(vvol) + 1))**2 * mn + 2 * 4 * 8 * Lrad(vvol) * mn        ! Ate, Aze, Ato, Azo
    end if
   endif ! end of if( Lcoordinatesingularity );

   do ii = 1, mn ! loop over Fourier harmonics;

    do ideriv = -2, 2 ! loop over derivatives; 14 Jan 13;

     SALLOCATE( Ate(vvol,ideriv,ii)%s, (0:Lrad(vvol)), zero )
     SALLOCATE( Aze(vvol,ideriv,ii)%s, (0:Lrad(vvol)), zero )
     SALLOCATE( Ato(vvol,ideriv,ii)%s, (0:Lrad(vvol)), zero )
     SALLOCATE( Azo(vvol,ideriv,ii)%s, (0:Lrad(vvol)), zero )

    enddo ! end of do ideriv;

    ;  ideriv =  0

     SALLOCATE( Ate(vvol,ideriv,ii)%i, (0:Lrad(vvol)), 0 ) ! degree of freedom index; 17 Jan 13;
     SALLOCATE( Aze(vvol,ideriv,ii)%i, (0:Lrad(vvol)), 0 )
     SALLOCATE( Ato(vvol,ideriv,ii)%i, (0:Lrad(vvol)), 0 )
     SALLOCATE( Azo(vvol,ideriv,ii)%i, (0:Lrad(vvol)), 0 )

   enddo ! end of do ii;
   
   select case( Linitgues ) ! for iterative solver of the Beltrami fields, an initial guess is required; 11 Mar 16;
   case( 0 )    ;
   case( 1 )    ; Ate(vvol,0,1)%s(0:1) = dtflux(vvol) * half ! this is an integrable approximation; NEEDS CHECKING; 26 Feb 13;
    ;           ; Aze(vvol,0,1)%s(0:1) = dpflux(vvol) * half ! this is an integrable approximation; NEEDS CHECKING; 26 Feb 13;
    if (Lcoordinatesingularity) then
    ;           ; Ate(vvol,0,1)%s(2) = dtflux(vvol) * half * half
    endif
   case( 2 )    ;                                            ! will call ra00aa below to read initial vector potential from file;
   case( 3 )    ;                                            ! the initial guess will be randomized, maximum is maxrndgues; 5 Mar 19;
    do ii = 1, mn ! loop over Fourier harmonics;

     do ideriv = -2, 2 ! loop over derivatives; 14 Jan 13;

      call random_number(Ate(vvol,ideriv,ii)%s)
      call random_number(Aze(vvol,ideriv,ii)%s)
      Ate(vvol,ideriv,ii)%s = Ate(vvol,ideriv,ii)%s * maxrndgues
      Aze(vvol,ideriv,ii)%s = Aze(vvol,ideriv,ii)%s * maxrndgues
      if (.not. YESstellsym) then
       call random_number(Ato(vvol,ideriv,ii)%s)
       call random_number(Azo(vvol,ideriv,ii)%s)
       Ato(vvol,ideriv,ii)%s = Ato(vvol,ideriv,ii)%s * maxrndgues
       Azo(vvol,ideriv,ii)%s = Azo(vvol,ideriv,ii)%s * maxrndgues
      endif

     enddo ! end of do ideriv;

    enddo ! end of do ii;

   end select

   idof = 0 ! degree of freedom index; reset to 0 in each volume;

   if( Lcoordinatesingularity ) then

    do ii = 1, mn ; mi = im(ii) ; ni = in(ii)

     do ll = 0, Lrad(vvol)
      ! Zernike is non zero only if ll>=mi and when they have the same parity
      if (ll>=mi .and. mod(mi+ll,2)==0)then
      ! We use the basis combination for m=0 and 1. They don't have ll=0 component.
      if (.not.((ll==0.and.mi==0).or.(ll==1.and.mi==1))) then
                                            ; idof = idof + 1 ; Ate(vvol,0,ii)%i(ll) = idof ! Zernike 30 Jun 19
      endif
      ;                                     ; idof = idof + 1 ; Aze(vvol,0,ii)%i(ll) = idof
      if( NOTstellsym .and. ii.gt.1 ) then
        if (.not.((ll==0.and.mi==0).or.(ll==1.and.mi==1))) then
                                            ; idof = idof + 1 ; Ato(vvol,0,ii)%i(ll) = idof ! Zernike 30 Jun 19
        endif
       ;                                    ; idof = idof + 1 ; Azo(vvol,0,ii)%i(ll) = idof
      endif ! NOTstellsym
      endif ! Zernike
     enddo ! end of do ll; 17 Jan 13;

    enddo ! end of do ii

    do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
     ! Lma is for Ate boundary condition on axis. For m=0 and 1, the boundary condition has been satisfied by basis recombination, so they are excluded.
     if ( mi.ne.0 .and. mi.ne.1     )  then ; idof = idof + 1 ; Lma(vvol,  ii)       = idof
     endif
     ! Lmb is for Aze boundary condition on axis. We only have that for m=0.
     if(  mi.eq.0                   ) then ; idof = idof + 1 ; Lmb(vvol,  ii)       = idof ! 18 May 16;
     endif
     ! Lme is for B.n at the outer boundary cos component. We don't have it for m=n=0.
     if(  ii.gt.1                   ) then ; idof = idof + 1 ; Lme(vvol,  ii)       = idof
     endif
     ! Lmg is for dtflux. We only have it for m=n=0.
     if(  ii.eq.1                   ) then ; idof = idof + 1 ; Lmg(vvol,  ii)       = idof
!   ! ;                                    ; idof = idof + 1 ; Lmh(vvol,  ii)       = idof ! no constraint on poloidal flux in innermost volume; 11 Mar 16;
     endif
     if( NOTstellsym ) then
      ! Lmc is for Ato boundary condition on axis. Same as Lma.
      if(  mi.ne.0 .and. mi.ne.1    ) then ; idof = idof + 1 ; Lmc(vvol,  ii)       = idof ! 18 May 16;
      endif
      ! Lmf is for B.n at the outer boundary sin component. Same as Lme.
      if(  ii.gt.1                  ) then ; idof = idof + 1 ; Lmf(vvol,  ii)       = idof ! 18 May 16;
     endif
     ! Lmd is for Azo on axis. We only have it for m=0, but not m=n=0.
     if(  ii.gt.1 .and. mi.eq.0     ) then ; idof = idof + 1 ; Lmd(vvol,  ii)       = idof ! 18 May 16;
     endif
     endif ! end of if( NOTstellsym ) ; 19 Jul 16;

    enddo ! end of do ii; 25 Jan 13;

    FATAL( preset, idof.ne.NAdof(vvol), need to count Beltrami degrees-of-freedom more carefully  for coordinate singularity )
    FATAL( preset, (idof+1)**2.ge.HUGE(idof)), NAdof too big, should be smaller than maximum of int32 type )

   else ! .not.Lcoordinatesingularity;

    do ii = 1, mn
     ! We use basis recombination method to ensure the inner boundary has At=Az=0. Therefore they don't have ll=0 component.
     do ll = 1, Lrad(vvol)                 ; idof = idof + 1 ; Ate(vvol,0,ii)%i(ll) = idof
      ;                                    ; idof = idof + 1 ; Aze(vvol,0,ii)%i(ll) = idof
      if( ii.gt.1 .and. NOTstellsym ) then ; idof = idof + 1 ; Ato(vvol,0,ii)%i(ll) = idof
       ;                                   ; idof = idof + 1 ; Azo(vvol,0,ii)%i(ll) = idof
      endif
     enddo ! end of do ll; 08 Feb 16;
    enddo

    do ii = 1, mn
     !;                                     ; idof = idof + 1 ; Lma(vvol,  ii)       = idof
     !;                                     ; idof = idof + 1 ; Lmb(vvol,  ii)       = idof
     !if(  ii.gt.1 .and. NOTstellsym ) then ; idof = idof + 1 ; Lmc(vvol,  ii)       = idof
     !;                                     ; idof = idof + 1 ; Lmd(vvol,  ii)       = idof
     !endif
     ! Lme is for B.n at the outer boundary cos component. We don't have it for m=n=0.
     if(  ii.gt.1                   ) then ; idof = idof + 1 ; Lme(vvol,  ii)       = idof
     endif
     ! Lmf is for B.n at the outer boundary sin component. Same as Lme
     if(  ii.gt.1 .and. NOTstellsym ) then ; idof = idof + 1 ; Lmf(vvol,  ii)       = idof
     endif
     ! Lmg and Lmh are the dtflux and dpflux constraint. Only present for m=n=0
     if(  ii.eq.1                   ) then ; idof = idof + 1 ; Lmg(vvol,  ii)       = idof
      ;                                    ; idof = idof + 1 ; Lmh(vvol,  ii)       = idof
     endif
    enddo ! end of do ii; 25 Jan 13;

   !if( Wpreset ) then
   ! do ii = 1, mn
   !  do ll = 0, Lrad(vvol)
   !   write(ounit,'("preset : " 10x " : myid="i3" ; ii="i4" ; ll="i4" : Ate = "i7" ;")') myid, ii, ll, Ate(vvol,0,ii)%i(ll)
   !   write(ounit,'("preset : " 10x " : myid="i3" ; ii="i4" ; ll="i4" : Aze = "i7" ;")') myid, ii, ll, Aze(vvol,0,ii)%i(ll)
   !   write(ounit,'("preset : " 10x " : myid="i3" ; ii="i4" ; ll="i4" : Ato = "i7" ;")') myid, ii, ll, Ato(vvol,0,ii)%i(ll)
   !   write(ounit,'("preset : " 10x " : myid="i3" ; ii="i4" ; ll="i4" : Azo = "i7" ;")') myid, ii, ll, Azo(vvol,0,ii)%i(ll)
   !  enddo
   !   write(ounit,'("preset : " 10x " : myid="i3" ; ii="i4" ;    "4x" : Lma = "i7" ;")') myid, ii,     Lma(vvol,  ii)
   !   write(ounit,'("preset : " 10x " : myid="i3" ; ii="i4" ;    "4x" : Lmb = "i7" ;")') myid, ii,     Lmb(vvol,  ii)
   !   write(ounit,'("preset : " 10x " : myid="i3" ; ii="i4" ;    "4x" : Lmc = "i7" ;")') myid, ii,     Lmc(vvol,  ii)
   !   write(ounit,'("preset : " 10x " : myid="i3" ; ii="i4" ;    "4x" : Lmd = "i7" ;")') myid, ii,     Lmd(vvol,  ii)
   !   write(ounit,'("preset : " 10x " : myid="i3" ; ii="i4" ;    "4x" : Lme = "i7" ;")') myid, ii,     Lme(vvol,  ii)
   !   write(ounit,'("preset : " 10x " : myid="i3" ; ii="i4" ;    "4x" : Lmf = "i7" ;")') myid, ii,     Lmf(vvol,  ii)
   !   write(ounit,'("preset : " 10x " : myid="i3" ; ii="i4" ;    "4x" : Lmg = "i7" ;")') myid, ii,     Lmg(vvol,  ii)
   !   write(ounit,'("preset : " 10x " : myid="i3" ; ii="i4" ;    "4x" : Lmh = "i7" ;")') myid, ii,     Lmh(vvol,  ii)
   ! enddo
   !endif

    FATAL( preset, idof.ne.NAdof(vvol), need to count degrees-of-freedom more carefully for new matrix )
    FATAL( preset, (idof+1)**2.ge.HUGE(idof)), NAdof too big, should be smaller than maximum of int32 type )

   endif ! end of if( Lcoordinatesingularity ) ;

   FATAL( preset, idof.ne.NAdof(vvol), impossible logic )

   do ii = 1, mn
      do jj = 0, Lrad(vvol)
        if (Ate(vvol,0,ii)%i(jj) == 0) Ate(vvol,0,ii)%s(jj) = zero
        if (Aze(vvol,0,ii)%i(jj) == 0) Aze(vvol,0,ii)%s(jj) = zero
        if (.not. YESstellsym) then
          if (Ato(vvol,0,ii)%i(jj) == 0) Azo(vvol,0,ii)%s(jj) = zero
          if (Azo(vvol,0,ii)%i(jj) == 0) Azo(vvol,0,ii)%s(jj) = zero
        end if
      end do !jj
   end do !ii
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  enddo ! end of do vvol = 1, Nvol loop;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( Linitgues.eq.2 ) then ; WCALL( preset, ra00aa, ('R') )  ! read initial guess for Beltrami field from file; 02 Jan 15;
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( myid.eq.0 ) then ! 17 Oct 12;
   cput = GETTIME
   write(ounit,'("preset : ", 10x ," : ")')
   write(ounit,'("preset : ",f10.2," : Nquad="i4" ; mn="i5" ; NGdof="i6" ; NAdof="16(i6",")" ...")') cput-cpus, Nquad, mn, NGdof, NAdof(1:min(Mvol,16))
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **workspace**

! Fourier transforms;

  Nt = max( Ndiscrete*4*Mpol, 1 ) ; Nz = max( Ndiscrete*4*Ntor, 1 ) ; Ntz = Nt*Nz ; soNtz = one / sqrt( one*Ntz ) ! exaggerated discrete resolution;

  ;                  ; hNt = Nt / 2
  if( Nz.gt.1 ) then ; hNz = Nz / 2
  else               ; hNz = 0
  endif

  if( myid.eq.0 ) then ! 17 Oct 12;
   cput = GETTIME
   write(ounit,'("preset : ", 10x ," : ")')
   write(ounit,'("preset : ",f10.2," : Nt="i6" ; Nz="i6" ; Ntz="i9" ;")') cput-cpus, Nt, Nz, Ntz
  endif

  SALLOCATE( iRij, (1:Ntz,0:Mvol), zero ) ! interface geometry in real space; ! 18 Jul 14;
  SALLOCATE( iZij, (1:Ntz,0:Mvol), zero ) !
  SALLOCATE( dRij, (1:Ntz,1:Mvol), zero ) ! interface geometry in real space; poloidal derivative; ! 18 Jul 14;
  SALLOCATE( dZij, (1:Ntz,1:Mvol), zero )
  SALLOCATE( tRij, (1:Ntz,0:Mvol), zero ) ! interface geometry in real space; poloidal derivative; ! 18 Jul 14;
  SALLOCATE( tZij, (1:Ntz,0:Mvol), zero )

  SALLOCATE(   Rij, (1:Ntz,0:3,0:3    ), zero ) ! these are used for inverse fft to reconstruct real space geometry from interpolated Fourier harmonics;
  SALLOCATE(   Zij, (1:Ntz,0:3,0:3    ), zero )
  SALLOCATE(   sg , (1:Ntz,0:3        ), zero )
  SALLOCATE( guvij, (1:Ntz,0:3,0:3,-1:3), zero ) ! need this on higher resolution grid for accurate Fourier decomposition;
  SALLOCATE( gvuij, (1:Ntz,0:3,0:3    ), zero ) ! need this on higher resolution grid for accurate Fourier decomposition; 10 Dec 15;

  if ((Lfindzero .eq. 2) .or. (Lcheck.eq.5 .or. LHevalues .or. LHevectors .or. LHmatrix .or. Lperturbed.eq.1)) then
    SALLOCATE( dRadR, (1:mn,0:1,0:1,1:mn), zero ) ! calculated in rzaxis; 19 Sep 16;
    SALLOCATE( dRadZ, (1:mn,0:1,0:1,1:mn), zero )
    SALLOCATE( dZadR, (1:mn,0:1,0:1,1:mn), zero )
    SALLOCATE( dZadZ, (1:mn,0:1,0:1,1:mn), zero )

    SALLOCATE( dRodR, (1:Ntz,0:3,1:mn), zero ) ! calculated in rzaxis; 19 Sep 16;
    SALLOCATE( dRodZ, (1:Ntz,0:3,1:mn), zero )
    SALLOCATE( dZodR, (1:Ntz,0:3,1:mn), zero )
    SALLOCATE( dZodZ, (1:Ntz,0:3,1:mn), zero )
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **goomne, goomno: metric information**
!> These are defined in metrix() , and used in ma00aa().
!>
!> **gssmne, gssmno: metric information**
!> These are defined in metrix() , and used in ma00aa().
!>
!> **gstmne, gstmno: metric information**
!> These are defined in metrix() , and used in ma00aa().
!>
!> **gszmne, gszmno: metric information**
!> These are defined in metrix() , and used in ma00aa().
!>
!> **gttmne, gttmno: metric information**
!> These are defined in metrix() , and used in ma00aa().
!>
!> **gtzmne, gtzmno: metric information**
!> These are defined in metrix() , and used in ma00aa().
!>
!> **gzzmne, gzzmno: metric information**
!> These are defined in metrix() , and used in ma00aa().
!>
  SALLOCATE( goomne, (0:mne, maxIquad), zero ) ! workspace for Fourier decomposition of metric terms;
  SALLOCATE( goomno, (0:mne, maxIquad), zero )
  SALLOCATE( gssmne, (0:mne, maxIquad), zero ) ! workspace for Fourier decomposition of metric terms;
  SALLOCATE( gssmno, (0:mne, maxIquad), zero )
  SALLOCATE( gstmne, (0:mne, maxIquad), zero ) ! workspace for Fourier decomposition of metric terms;
  SALLOCATE( gstmno, (0:mne, maxIquad), zero )
  SALLOCATE( gszmne, (0:mne, maxIquad), zero ) ! workspace for Fourier decomposition of metric terms;
  SALLOCATE( gszmno, (0:mne, maxIquad), zero )
  SALLOCATE( gttmne, (0:mne, maxIquad), zero ) ! workspace for Fourier decomposition of metric terms;
  SALLOCATE( gttmno, (0:mne, maxIquad), zero )
  SALLOCATE( gtzmne, (0:mne, maxIquad), zero ) ! workspace for Fourier decomposition of metric terms;
  SALLOCATE( gtzmno, (0:mne, maxIquad), zero )
  SALLOCATE( gzzmne, (0:mne, maxIquad), zero ) ! workspace for Fourier decomposition of metric terms;
  SALLOCATE( gzzmno, (0:mne, maxIquad), zero )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  SALLOCATE( ijreal, (1:Ntz), zero ) ! real space grid;
  SALLOCATE( ijimag, (1:Ntz), zero )
  SALLOCATE( jireal, (1:Ntz), zero )
  SALLOCATE( jiimag, (1:Ntz), zero )

  SALLOCATE( jkreal, (1:Ntz), zero )
  SALLOCATE( jkimag, (1:Ntz), zero )
  SALLOCATE( kjreal, (1:Ntz), zero )
  SALLOCATE( kjimag, (1:Ntz), zero )

  SALLOCATE( cplxin,  (1:Nt,1:Nz,nthreads), zero )
  SALLOCATE( cplxout, (1:Nt,1:Nz,nthreads), zero )

  ! Create and save optimal plans for forward and inverse 2D fast Fourier transforms with FFTW. -JAB; 25 Jul 2017
  planf = fftw_plan_dft_2d( Nz, Nt, cplxin(:,:,1), cplxout(:,:,1), FFTW_FORWARD,  FFTW_MEASURE + FFTW_DESTROY_INPUT )
  planb = fftw_plan_dft_2d( Nz, Nt, cplxin(:,:,1), cplxout(:,:,1), FFTW_BACKWARD, FFTW_MEASURE + FFTW_DESTROY_INPUT )

  SALLOCATE( efmn, (1:mne), zero ) ! Fourier harmonics workspace; 24 Apr 13;
  SALLOCATE( ofmn, (1:mne), zero )
  SALLOCATE( cfmn, (1:mne), zero )
  SALLOCATE( sfmn, (1:mne), zero )
  SALLOCATE( evmn, (1:mne), zero )
  SALLOCATE( odmn, (1:mne), zero )
  SALLOCATE( comn, (1:mne), zero )
  SALLOCATE( simn, (1:mne), zero )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **cosi(1:Ntz,1:mn) and sini(1:Ntz,1:mn)**
!>
!> <ul>
!> <li> Trigonometric factors used in various Fast Fourier transforms, where
!>       \f{eqnarray}{ \texttt{cosi}_{j,i} & = & \cos( m_i \theta_j - n_i \zeta_j ), \\
!>                     \texttt{sini}_{j,i} & = & \sin( m_i \theta_j - n_i \zeta_j ).
!>       \f} </li>
!> </ul>

  SALLOCATE( gteta, (1:Ntz), zero )
  SALLOCATE( gzeta, (1:Ntz), zero )

  SALLOCATE( cosi, (1:Ntz,1:mn), zero )
  SALLOCATE( sini, (1:Ntz,1:mn), zero )

  FATAL( preset, Nz.eq.0, illegal division )
  FATAL( preset, Nt.eq.0, illegal division )

  do ii = 1, mn ; mi = im(ii) ; ni = in(ii) ! loop over Fourier harmonics;

  do kk = 0, Nz-1 ; zeta = kk * pi2nfp / Nz
    do jj = 0, Nt-1 ; teta = jj * pi2    / Nt ; jk = 1 + jj + kk*Nt ; arg = mi * teta - ni * zeta
    gteta(jk) = teta
    gzeta(jk) = zeta
    cosi(jk,ii) = cos(arg)
    sini(jk,ii) = sin(arg)
    enddo
  enddo

  enddo ! end of do ii; 13 May 13;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef DEBUG

  if( Wpreset .and. myid.eq.0 ) then

   write(ounit,'("preset : ",10x," : checking FFT and inverse FFT ;")')

   do imn = 1, mn ; mm = im(imn) ; nn = in(imn) ! in should include the Nfp factor; SRH: 27 Feb 18;

    ijreal(1:Ntz) = zero ; ijimag(1:Ntz) = zero

    do kk = 0, Nz-1 ; zeta = kk * pi2nfp / Nz

     do jj = 0, Nt-1 ; teta = jj * pi2    / Nt ; jk = 1 + jj + kk*Nt

      ijreal(jk) = cos( mm * teta - nn * zeta ) ; ijimag(jk) = sin( mm * teta - nn * zeta )

     enddo ! end of do jj; SRH: 27 Feb 18;

    enddo ! end of do kk; SRH: 27 Feb 18;

    jkreal = ijreal ; jkimag = ijimag

    ifail = 0 !                                                              even        odd         cos         sin
    call tfft( Nt, Nz, ijreal(1:Ntz), ijimag(1:Ntz), mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), ifail )

    do ii = 1, mn

     if( abs(efmn(ii))+abs(ofmn(ii))+abs(cfmn(ii))+abs(sfmn(ii)).gt.small ) write(ounit,2000) mm, nn, im(ii), in(ii), efmn(ii), ofmn(ii), cfmn(ii), sfmn(ii)

2000 format("preset : ",10x," : (",i3,",",i3," ) = (",i3,",",i3," ) : "2f15.5" ; "2f15.5" ;")

    enddo ! end of do ii; SRH: 27 Feb 18;

    call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, jireal(1:Ntz), jiimag(1:Ntz) )

    error = ( sum((jkreal(1:Ntz)-jireal(1:Ntz))**2) + sum((jkimag(1:Ntz)-jiimag(1:Ntz))**2) ) / Ntz

    write(ounit,'("preset : ",10x," : (",i3,",",i3," ) : error = ",es13.5," ;")') mm, nn, error

   enddo ! end of do imn; SRH: 27 Feb 18;

  endif ! end of if( myid.eq.0 ) ; SRH: 27 Feb 18;

#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( Igeometry.eq.3 .and. iRbc(1,0).lt.small ) then ! have not yet assigned coordinate axis; see global;readin for user-supplied Rac, Zas, etc. ; 19 Jul 16;

   select case( Linitialize )
   case( :-1 ) ; vvol = Nvol + Linitialize
   case(   0 ) ; vvol =    1 ! this is really a dummy; no interpolation of interface geometry is required; packxi calls rzaxis with lvol=1; 19 Jul 16;
   case(   1 ) ; vvol = Nvol
   case(   2 ) ; vvol = Mvol
   end select

   WCALL( preset, rzaxis, ( Mvol, mn, iRbc(1:mn,0:Mvol), iZbs(1:mn,0:Mvol), iRbs(1:mn,0:Mvol), iZbc(1:mn,0:Mvol), vvol, .false. ) ) ! set coordinate axis; 19 Jul 16;

  endif ! end of if( Igeometry.eq.3 ) then ; 19 Jul 16;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **psifactor(1:mn,1:Mvol): coordinate "pre-conditioning" factor**
!>
!> <ul>
!> <li> In toroidal geometry, the coordinate "pre-conditioning" factor is
!>       \f{eqnarray}{ f_{j,v} \equiv \left\{
!>       \begin{array}{lcccccc}\psi_{t,v}^{0    }&,&\mbox{for $m_j=0$}, \\
!>                             \psi_{t,v}^{m_j/2}&,&\mbox{otherwise}.
!>       \end{array}\right.
!>       \f}
!>       where \f$\psi_{t,v} \equiv\,\f$\c tflux is the (normalized?) toroidal flux enclosed by the \f$v\f$-th interface. </li>
!> <li> \c psifactor is used in packxi(), dforce() and hesian(). </li>
!> <li> \c inifactor is similarly constructed, with
!>       \f{eqnarray}{ f_{j,v} \equiv \left\{
!>       \begin{array}{lcccccc}\psi_{t,v}^{ 1 /2}&,&\mbox{for $m_j=0$}, \\
!>                             \psi_{t,v}^{m_j/2}&,&\mbox{otherwise}.
!>       \end{array}\right.
!>       \f}
!>       and used only for the initialization of the surfaces taking into account axis information if provided. </li>
!> </ul>

  SALLOCATE( psifactor, (1:mn,1:Mvol), zero )
  SALLOCATE( inifactor, (1:mn,1:Mvol), zero )

  psifactor(1:mn,1:Mvol) = one
  inifactor(1:mn,1:Mvol) = one

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
     if( im(ii).eq.0 ) then ; psifactor(ii,vvol) = Rscale * tflux(vvol)**zero       ! 08 Feb 16;
                            ; inifactor(ii,vvol) = Rscale * tflux(vvol)**half       ! 17 Dec 18;
     else                   ; psifactor(ii,vvol) = Rscale * tflux(vvol)**halfmm(ii) ! 29 Apr 15;
                            ; inifactor(ii,vvol) = Rscale * tflux(vvol)**halfmm(ii) ! 17 Dec 18
     endif
    enddo
   enddo

  case default

   FATAL( readin, .true., invalid Igeometry for construction of psifactor )

  end select

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( Linitialize.ne.0 ) then ! interpolate / extrapolate interior interface geometry; 19 Jul 16;

   select case( Igeometry )

   case( 1 ) ! Cartesian; 29 Apr 14;

   !FATAL( preset, Linitialize.ne.1, geometrical initialization under construction for Cartesian ) ! 14 Apr 17;

    do vvol = 1, Nvol
     ;iRbc(1:mn,vvol) = iRbc(1:mn,Mvol) * tflux(vvol) / tflux(Mvol) ! 14 Apr 17;
     if( NOTstellsym ) then
      iRbs(2:mn,vvol) = iRbs(2:mn,Mvol) * tflux(vvol) / tflux(Mvol) ! 14 Apr 17;
     endif
    enddo

   case( 2 ) ! cylindrical - standard; 20 Apr 13;

    FATAL( preset, Linitialize.ne.1, geometrical initialization under construction for cylindrical )

    do vvol = 1, Nvol-1
     ;iRbc(1:mn,vvol) = iRbc(1:mn,Nvol) * psifactor(1:mn,vvol)
     if( NOTstellsym ) then
      iRbs(2:mn,vvol) = iRbs(2:mn,Nvol) * psifactor(2:mn,vvol)
     endif
    enddo

   case( 3 ) ! toroidal; 20 Apr 13;

    FATAL( preset, Linitialize.lt.0, geometrical initialization under construction for toroidal ) ! see commented-out source below; 19 Jul 16;

    lvol = Nvol-1 + Linitialize

    FATAL( preset, lvol.gt.Mvol, perhaps illegal combination of Linitialize and Lfreebound )

!    do vvol = 1, Nvol-1       ! 19 Jul 16;
!     ;iRbc(1:mn,vvol) = iRbc(1:mn,0) + ( iRbc(1:mn,Nvol) - iRbc(1:mn,0) ) * psifactor(1:mn,vvol) ! 19 Jul 16;
!     ;iZbs(2:mn,vvol) = iZbs(2:mn,0) + ( iZbs(2:mn,Nvol) - iZbs(2:mn,0) ) * psifactor(2:mn,vvol) ! 19 Jul 16;
!     if( NOTstellsym ) then ! 19 Jul 16;
!      iRbs(2:mn,vvol) = iRbs(2:mn,0) + ( iRbs(2:mn,Nvol) - iRbs(2:mn,0) ) * psifactor(2:mn,vvol) ! 19 Jul 16;
!      iZbc(1:mn,vvol) = iZbc(1:mn,0) + ( iZbc(1:mn,Nvol) - iZbc(1:mn,0) ) * psifactor(1:mn,vvol) ! 19 Jul 16;
!     endif ! 19 Jul 16;
!    enddo
!
    do vvol = 1, lvol-1
     ;iRbc(1:mn,vvol) = iRbc(1:mn,0) + ( iRbc(1:mn,lvol) - iRbc(1:mn,0) ) * ( inifactor(1:mn,vvol) / Rscale ) / tflux(lvol)**halfmm(1:mn)
     ;iZbs(2:mn,vvol) = iZbs(2:mn,0) + ( iZbs(2:mn,lvol) - iZbs(2:mn,0) ) * ( inifactor(2:mn,vvol) / Rscale ) / tflux(lvol)**halfmm(2:mn)
     if( NOTstellsym ) then
      iRbs(2:mn,vvol) = iRbs(2:mn,0) + ( iRbs(2:mn,lvol) - iRbs(2:mn,0) ) * ( inifactor(2:mn,vvol) / Rscale ) / tflux(lvol)**halfmm(2:mn)
      iZbc(1:mn,vvol) = iZbc(1:mn,0) + ( iZbc(1:mn,lvol) - iZbc(1:mn,0) ) * ( inifactor(1:mn,vvol) / Rscale ) / tflux(lvol)**halfmm(1:mn)
     endif
    enddo

!   do vvol = 1, Nvol+Linitialize-1
!    ;iRbc(1:mn,vvol) = iRbc(1:mn,0) + ( iRbc(1:mn,Nvol+Linitialize) - iRbc(1:mn,0) ) * psifactor(1:mn,vvol) / tflux(Nvol+Linitialize)**halfmm(1:mn)
!    ;iZbs(2:mn,vvol) = iZbs(2:mn,0) + ( iZbs(2:mn,Nvol+Linitialize) - iZbs(2:mn,0) ) * psifactor(2:mn,vvol) / tflux(Nvol+Linitialize)**halfmm(2:mn)
!    if( NOTstellsym ) then
!     iRbs(2:mn,vvol) = iRbs(2:mn,0) + ( iRbs(2:mn,Nvol+Linitialize) - iRbs(2:mn,0) ) * psifactor(2:mn,vvol) / tflux(Nvol+Linitialize)**halfmm(2:mn)
!     iZbc(1:mn,vvol) = iZbc(1:mn,0) + ( iZbc(1:mn,Nvol+Linitialize) - iZbc(1:mn,0) ) * psifactor(1:mn,vvol) / tflux(Nvol+Linitialize)**halfmm(1:mn)
!    endif
!   enddo

   end select ! matches select case( Igeometry ); 19 Jul 16;

  endif ! matches if( Linitialize.ne.0 ) then; 19 Jul 16;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **Bsupumn and Bsupvmn**

  SALLOCATE( Bsupumn, (1:Nvol,0:1,1:mn), zero ) ! Fourier components of {\bf B}\cdot\nabla \theta on boundary; required for virtual casing;
  SALLOCATE( Bsupvmn, (1:Nvol,0:1,1:mn), zero ) ! Fourier components of {\bf B}\cdot\nabla \zeta  on boundary;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **diotadxup and glambda: transformation to straight fieldline angle**
!>
!> <ul>
!> <li> Given the Beltrami fields in any volume, the rotational-transform on the adjacent interfaces
!>       may be determined (in tr00ab()) by constructing the straight fieldline angle on the interfaces. </li>
!> <li> The rotational transform on the inner or outer interface of a given volume depends on the magnetic field in that volume,
!>       i.e. \f${{\,\iota\!\!\!}-}_\pm = {{\,\iota\!\!\!}-}({\bf B}_\pm)\f$,
!>       so that
!>       \f{eqnarray}{ \delta {{\,\iota\!\!\!}-}_\pm = \frac{\partial {{\,\iota\!\!\!}-}_\pm}{\partial {\bf B}_\pm} \cdot \delta {\bf B_\pm}.
!>       \f} </li>
!> <li> The magnetic field depends on the Fourier harmonics of both the inner and outer interface geometry (represented here as \f$x_j\f$),
!>       the helicity multiplier, and the enclosed poloidal flux, i.e. \f${\bf B_\pm} = {\bf B_\pm}(x_j, \mu, \Delta \psi_p)\f$, so that
!>       \f{eqnarray}{ \delta {\bf B_\pm} = \frac{\partial {\bf B}_\pm}{\partial x_j          } \delta x_j
!>                                        + \frac{\partial {\bf B}_\pm}{\partial \mu          } \delta \mu
!>                                        + \frac{\partial {\bf B}_\pm}{\partial \Delta \psi_p} \delta \Delta \psi_p.
!>       \f} </li>
!> <li> The rotational-transforms, thus, can be considered to be functions of the geometry, the helicity-multiplier and the enclosed poloidal flux,
!>       \f${{\,\iota\!\!\!}-}_{\pm} = {{\,\iota\!\!\!}-}_{\pm}(x_j,\mu,\Delta\psi_p)\f$. </li>
!> <li> The rotational-transform, and its derivatives, on the inner and outer interfaces of each volume is stored in
!>       \c diotadxup(0:1,-1:2,1:Mvol) .
!>       The indices label:
!>       <ul>
!>       <li> the first index labels the inner or outer interface, </li>
!>       <li> the the second one labels derivative, with  </li>
!>       <ul> <li>\c -1 : indicating the derivative with respect to the interface geometry,
!>                                 i.e. \f$\displaystyle \frac{\partial {{\,\iota\!\!\!}-}_{\pm}}{\partial x_j}\f$, </li>
!>            <li>\c 0 : the rotational-transform itself, </li>
!>            <li>\c 1,2 : the derivatives with respec to \f$\mu\f$ and \f$\Delta \psi_p\f$,
!>                    i.e. \f$\displaystyle \frac{\partial {{\,\iota\!\!\!}-}_{\pm}}{\partial \mu}\f$ and
!>                         \f$\displaystyle \frac{\partial {{\,\iota\!\!\!}-}_{\pm}}{\partial \Delta \psi_p}\f$; </li>
!>       </ul> </li>
!>       <li>The third index labels volume. </li>
!>       </ul> </li>
!> <li> The values of \c diotadxup are assigned in mp00aa() after calling tr00ab(). </li>
!> </ul>

  SALLOCATE( diotadxup, (0:1,-1:2,1:Mvol), zero ) ! measured rotational transform on inner/outer interfaces in each annulus;
  SALLOCATE( dItGpdxtp, (0:1,-1:2,1:Mvol), zero ) ! measured plasma and linking currents                                   ;

  SALLOCATE( glambda, (1:Ntz+1,0:2,0:1,1:Mvol), zero ) ! save initial guesses for iterative calculation of rotational-transform; 21 Apr 13;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! Construction of ``force'';

  SALLOCATE( Bemn, (1:mn,1:Mvol,0:1), zero )
  SALLOCATE( Bomn, (1:mn,1:Mvol,0:1), zero )
  SALLOCATE( Iomn, (1:mn,1:Mvol    ), zero )
  SALLOCATE( Iemn, (1:mn,1:Mvol    ), zero )
  SALLOCATE( Somn, (1:mn,1:Mvol,0:1), zero )
  SALLOCATE( Semn, (1:mn,1:Mvol,0:1), zero )
  SALLOCATE( Pomn, (1:mn,1:Mvol,0:2), zero )
  SALLOCATE( Pemn, (1:mn,1:Mvol,0:2), zero )

  SALLOCATE( BBe , (1:Mvol-1), zero )
  SALLOCATE( IIo , (1:Mvol-1), zero )
  SALLOCATE( BBo , (1:Mvol-1), zero )
  SALLOCATE( IIe , (1:Mvol-1), zero )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  SALLOCATE( Btemn, (1:mn,0:1,1:Mvol), zero ) ! these are declared in global, calculated in sc00aa, broadcast in xspech, and written to file in hdfint;
  SALLOCATE( Bzemn, (1:mn,0:1,1:Mvol), zero )
  SALLOCATE( Btomn, (1:mn,0:1,1:Mvol), zero )
  SALLOCATE( Bzomn, (1:mn,0:1,1:Mvol), zero )

  SALLOCATE( Bloweremn, (1:mn, 3), zero) ! these are declared in global, calculated in getbco, used in mtrxhs
  SALLOCATE( Bloweromn, (1:mn, 3), zero)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **vvolume, lBBintegral and lABintegral**
!>
!> <ul>
!> <li> volume integrals
!>       \f{eqnarray}{ \texttt{vvolume(i)}     &=& \int_{{\cal V}_i}                       \, dv\\
!>                     \texttt{lBBintegral(i)} &=& \int_{{\cal V}_i} {\bf B} \cdot {\bf B} \, dv\\
!>                     \texttt{lABintegral(i)} &=& \int_{{\cal V}_i} {\bf A} \cdot {\bf B} \, dv
!>       \f} </li>
!> </ul>

! Allocate matrix to store the last solution of GMRES as initialization
  LILUprecond = .false.
  if (Lmatsolver.eq.2 .or. Lmatsolver.eq.3) then ! use GMRES
    SALLOCATE(GMRESlastsolution, (MAXVAL(NAdof),0:2,1:Mvol), zero )
    GMRESlastsolution = zero
    if (LGMRESprec .eq. 1) LILUprecond = .true.
  endif

  if (Lmatsolver.eq.3) then
    YESMatrixFree = .true.
    NOTMatrixFree = .false.
  else
    YESMatrixFree = .false.
    NOTMatrixFree = .true.
  endif

  SALLOCATE( vvolume    , (1:Mvol), zero ) ! volume integral of \sqrt g;
  SALLOCATE( lBBintegral, (1:Mvol), zero ) ! volume integral of B.B    ;
  SALLOCATE( lABintegral, (1:Mvol), zero ) ! volume integral of A.B    ;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( YESstellsym ) lmns = 1 + (mns-1)           ! number of independent degrees-of-freedom in angle transformation; 30 Jan 13;
  if( NOTstellsym ) lmns = 1 + (mns-1) + (mns-1) ! number of independent degrees-of-freedom in angle transformation; 30 Jan 13;

  SALLOCATE( dlambdaout, (1:lmns,1:Mvol,0:1), zero )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if (Lfreebound > 0) then ! Only do for free-boundary; 7 Nov 18;

    SALLOCATE( Dxyz, (1:3,1:Ntz), zero ) ! Cartesian components of computational boundary; position; 14 Apr 17;
    SALLOCATE( Nxyz, (1:3,1:Ntz), zero ) ! Cartesian components of computational boundary; normal  ; 14 Apr 17;

    SALLOCATE( Jxyz, (1:Ntz,1:3), zero ) ! Cartesian components of virtual casing surface current; needs to be recalculated at each iteration;

    lvol = Mvol ; lss = one ; Lcurvature = 1 ; Lcoordinatesingularity = .false. ! will only require normal field on outer interface = computational boundary;

    WCALL( preset, coords,( lvol, lss, Lcurvature, Ntz, mn ) ) ! will need Rij, Zij; THE COMPUTATIONAL BOUNDARY DOES NOT CHANGE;

    do kk = 0, Nz-1 ; zeta = kk * pi2nfp / Nz

     if( Igeometry.eq.3 ) then ; cszeta(0:1) = (/ cos(zeta), sin(zeta) /)
     endif

     do jj = 0, Nt-1 ; teta = jj * pi2    / Nt ; jk = 1 + jj + kk*Nt

      select case( Igeometry )
      case( 1 ) ! Igeometry = 1 ;
       Dxyz(1:3,jk) = (/   teta       ,  zeta       ,   Rij(jk,0,0) /)
       Nxyz(1:3,jk) = (/ - Rij(jk,2,0), -Rij(jk,3,0),   one         /)
      case( 2 ) ! Igeometry = 2 ;
       FATAL( bnorml, .true., free-boundary calculations not yet implemented in cylindrical geometry )
      case( 3 ) ! Igeometry = 3 ;
       Dxyz(1:3,jk) = (/   Rij(jk,0,0) * cszeta(0), Rij(jk,0,0) * cszeta(1), Zij(jk,0,0) /)
       Nxyz(1:3,jk) = (/   Rij(jk,2,0) * cszeta(1) * Zij(jk,3,0) - Zij(jk,2,0) * ( Rij(jk,3,0) * cszeta(1) + Rij(jk,0,0) * cszeta(0) ), &
                         - Rij(jk,2,0) * cszeta(0) * Zij(jk,3,0) + Zij(jk,2,0) * ( Rij(jk,3,0) * cszeta(0) - Rij(jk,0,0) * cszeta(1) ), &
                           Rij(jk,0,0)             * Rij(jk,2,0) /)
      end select ! end of select case( Igeometry ) ; 09 Mar 17;

     enddo ! end of do jj; 14 Apr 17;

    enddo ! end of do kk; 14 Apr 17;

  endif ! Lfreebound > 1; 7 Nov 18;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( Lconstraint .EQ. 3) then
    Localconstraint = .false.
  else
    Localconstraint = .true.
  endif
  
  RETURN(preset)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine preset

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
