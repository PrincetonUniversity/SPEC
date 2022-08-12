!> \file
!> \brief Calculates volume integrals of Chebyshev-polynomials and covariant field for Hessian computation.

!latex \briefly{Calculates volume integrals of Chebyshev polynomials and covariant field products.}

!latex \calledby{\link{dforce}}
!latex \calls{\link{getbco}}

!latex \tableofcontents

!l tex \newcommand{\bT}[1]{{\overline T}_{#1}}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{chebyshev-metric information}

!latex \begin{enumerate}

!latex \item The following quantities are calculated:

!latex       \be \verb+DToocc(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}  \; \ooint \cos\a_i \cos\a_j                  \\
!latex           \verb+DToocs(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}  \; \ooint \cos\a_i \sin\a_j                  \\
!latex           \verb+DToosc(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}  \; \ooint \sin\a_i \cos\a_j                  \\
!latex           \verb+DTooss(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}  \; \ooint \sin\a_i \sin\a_j
!latex       \ee

!latex       \be \verb+TTsscc(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}  \; \bT{p,j}  \; \ooint \cos\a_i \cos\a_j \; \bar g_{\s\s} \\
!latex           \verb+TTsscs(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}  \; \bT{p,j}  \; \ooint \cos\a_i \sin\a_j \; \bar g_{\s\s} \\
!latex           \verb+TTsssc(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}  \; \bT{p,j}  \; \ooint \sin\a_i \cos\a_j \; \bar g_{\s\s} \\
!latex           \verb+TTssss(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}  \; \bT{p,j}  \; \ooint \sin\a_i \sin\a_j \; \bar g_{\s\s}
!latex       \ee

!latex       \be \verb+TDstcc(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}  \; \bT{p,j}' \; \ooint \cos\a_i \cos\a_j \; \bar g_{\s\t} \\
!latex           \verb+TDstcs(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}  \; \bT{p,j}' \; \ooint \cos\a_i \sin\a_j \; \bar g_{\s\t} \\
!latex           \verb+TDsTsc(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}  \; \bT{p,j}' \; \ooint \sin\a_i \cos\a_j \; \bar g_{\s\t} \\
!latex           \verb+TDsTss(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}  \; \bT{p,j}' \; \ooint \sin\a_i \sin\a_j \; \bar g_{\s\t}
!latex       \ee

!latex       \be \verb+TDstcc(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}  \; \bT{p,j}' \; \ooint \cos\a_i \cos\a_j \; \bar g_{\s\z} \\
!latex           \verb+TDstcs(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}  \; \bT{p,j}' \; \ooint \cos\a_i \sin\a_j \; \bar g_{\s\z} \\
!latex           \verb+TDsTsc(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}  \; \bT{p,j}' \; \ooint \sin\a_i \cos\a_j \; \bar g_{\s\z} \\
!latex           \verb+TDsTss(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}  \; \bT{p,j}' \; \ooint \sin\a_i \sin\a_j \; \bar g_{\s\z}
!latex       \ee

!latex       \be \verb+DDstcc(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}' \; \ooint \cos\a_i \cos\a_j \; \bar g_{\t\t} \\
!latex           \verb+DDstcs(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}' \; \ooint \cos\a_i \sin\a_j \; \bar g_{\t\t} \\
!latex           \verb+DDsTsc(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}' \; \ooint \sin\a_i \cos\a_j \; \bar g_{\t\t} \\
!latex           \verb+DDsTss(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}' \; \ooint \sin\a_i \sin\a_j \; \bar g_{\t\t}
!latex       \ee

!latex       \be \verb+DDstcc(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}' \; \ooint \cos\a_i \cos\a_j \; \bar g_{\t\z} \\
!latex           \verb+DDstcs(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}' \; \ooint \cos\a_i \sin\a_j \; \bar g_{\t\z} \\
!latex           \verb+DDsTsc(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}' \; \ooint \sin\a_i \cos\a_j \; \bar g_{\t\z} \\
!latex           \verb+DDsTss(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}' \; \ooint \sin\a_i \sin\a_j \; \bar g_{\t\z}
!latex       \ee

!latex       \be \verb+DDstcc(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}' \; \ooint \cos\a_i \cos\a_j \; \bar g_{\z\z} \\
!latex           \verb+DDstcs(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}' \; \ooint \cos\a_i \sin\a_j \; \bar g_{\z\z} \\
!latex           \verb+DDsTsc(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}' \; \ooint \sin\a_i \cos\a_j \; \bar g_{\z\z} \\
!latex           \verb+DDsTss(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}' \; \ooint \sin\a_i \sin\a_j \; \bar g_{\z\z}
!latex       \ee

!latex       where $\bT{l,i}\equiv T_l \, \bar s^{m_i/2}$ if the domain includes the coordinate singularity, and $\bT{l,i}\equiv T_l$ if not;
!latex       and $\bar g_{\mu\nu} \equiv g_{\mu\nu} / \sqrt g$.

!latex \item The double-angle formulae are used to reduce the above expressions to the Fourier harmonics of $\bar g_{\mu\nu}$:
!latex       see \internal{kija} and \internal{kijs}, which are defined in \link{preset}.

!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

module intghs_module
  use mod_kinds, only: wp => dp
  !> \brief This calculates the integral of something related to matrix-vector-multiplication.
  !> \todo Zhisong might need to update the documentation of this type.
  type intghs_workspace
    real(wp), allocatable   :: efmn(:,:)        !< This is efmn.
    real(wp), allocatable   :: ofmn(:,:)        !< This is ofmn.
    real(wp), allocatable   :: cfmn(:,:)        !<
    real(wp), allocatable   :: sfmn(:,:)        !<
    real(wp), allocatable   :: evmn(:,:)        !<
    real(wp), allocatable   :: odmn(:,:)        !<
    real(wp), allocatable   :: ijreal(:,:)      !<
    real(wp), allocatable   :: jireal(:,:)      !<
    real(wp), allocatable   :: jkreal(:,:)      !<
    real(wp), allocatable   :: kjreal(:,:)      !<
    real(wp), allocatable   :: Bloweremn(:,:,:) !<
    real(wp), allocatable   :: Bloweromn(:,:,:) !<
    real(wp), allocatable   :: gBupper(:,:,:)   !<
    real(wp), allocatable   :: Blower(:,:,:)    !<
    real(wp), allocatable   :: basis(:,:,:,:)   !<
  end type

  TYPE(intghs_workspace) :: wk !< This is an instance of the intghs_workspace type.

end module intghs_module

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> \brief Calculates volume integrals of Chebyshev polynomials and covariant field products.
!>
!> @param lquad
!> @param mn
!> @param lvol
!> @param lrad
!> @param idx
subroutine intghs( lquad, mn, lvol, lrad, idx )
  use mod_kinds, only: wp => dp
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  use constants, only : zero, half, one, two, pi, pi2

  use numerical, only : vsmall, small, sqrtmachprec

  use fileunits, only : ounit

  use inputlist, only : mpol, Wintghs, Wmacros

  use cputiming, only : Tintghs

  use allglobal, only : myid, ncpu, cpus, MPI_COMM_SPEC, &
                        Mvol, im, in, mne, Ntz, Nt, Nz, &
                        YESstellsym, NOTstellsym, &
                        gaussianweight, gaussianabscissae, &
                        Tsc, Tss, Dtc, Dts, Dzc, Dzs, &
                        Ttc, Tts, Tzc, Tzs, &
                        Lcoordinatesingularity, &
                        pi2pi2nfp, pi2pi2nfphalf, dBdX, &
                        Ntz, NOTstellsym, dBdX, Lsavedguvij, &
                        Ate, Ato, Aze, Azo, &
                        sg, guvij, guvijsave

  use intghs_module
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


#ifdef OPENMP
  USE OMP_LIB
#endif
  use mpi
  implicit none
  integer   :: ierr, astat, ios, nthreads, ithread
  real(wp)      :: cput, cpui, cpuo=0 ! cpu time; cpu initial; cpu old; 31 Jan 13;


  integer, intent(in) :: lquad, mn, lvol, lrad, idx

  integer             :: jquad, ll, pp, ll1, pp1, uv, ii, jj, io, mn2, lp2, mn2_max, lp2_max, nele, ideriv, ifail, Lcurvature
  integer             :: mi, ni, bid

  real(wp)                :: lss, jthweight, Tl, Dl, sbar, dfactor, ik, w(1:lquad)


  cpui = MPI_WTIME()
  cpuo = cpui
#ifdef OPENMP
  nthreads = omp_get_max_threads()
#else
  nthreads = 1
#endif


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  Tss = zero
  Dtc = zero
  Dzc = zero

  if (.not.dBdX%L) then
    Ttc = zero
    Tzc = zero
  endif

  if (NOTstellsym) then
    Tsc = zero
    Dts = zero
    Dzs = zero
    if (.not.dBdX%L) then
      Tts = zero
      Tzs = zero
    endif

  endif !NOTstellsym

  if( dBdX%L ) then ; Lcurvature = 3 ; ideriv = 1
  else              ; Lcurvature = 1 ; ideriv = 0
  endif

  if (.not. Lsavedguvij) then

   cput = MPI_WTIME()
   Tintghs = Tintghs + ( cput-cpuo )
   call compute_guvijsave(lquad, lvol, ideriv, Lcurvature)
   cpuo = MPI_WTIME()

  endif

  do jquad = 1, lquad
    lss = gaussianabscissae(jquad,lvol) ; jthweight = gaussianweight(jquad,lvol)
    sbar = (lss + one) * half
    if (Lcoordinatesingularity) then
      call get_zernike(sbar, lrad, mpol, wk%basis(:,:,0:1,jquad)) ! use Zernike polynomials 29 Jun 19;
    else
      call get_cheby(lss, lrad, wk%basis(:,0,0:1,jquad))
    endif
  enddo

  wk%efmn = zero ; wk%ofmn = zero ; wk%cfmn = zero ; wk%sfmn = zero
  wk%evmn = zero ; wk%odmn = zero ;

  wk%ijreal = zero ; wk%jireal = zero ; wk%jkreal = zero ; wk%kjreal = zero

  wk%gBupper = zero; wk%Blower = zero;
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!$OMP PARALLEL DO PRIVATE(jquad,ll1,lss,sbar,Tl,Dl,ii,jj,ll,mi,ni,ik) SHARED(lquad,mn,lvol,lrad,idx)
  do jquad = 1, lquad

    ! A_t  : wk%ijreal, wk%jireal
    ! A_z  : wk%jkreal, wk%kjreal
    ! gB^s : wk%evmn, wk%odmn
    ! gB^t : wk%cfmn, wk%sfmn
    ! gB^z : wk%efmn, wk%ofmn

    do ii = 1, mn ; mi = im(ii) ; ni = in(ii) ! loop over Fourier harmonics;

      if (Lcoordinatesingularity) then
        do ll = mi, lrad, 2 ! loop over Zernike polynomials; Lrad is the radial resolution;
        ;                      ; wk%efmn(ii,jquad) = wk%efmn(ii,jquad) + Ate(lvol,idx,ii)%s(ll) * ( wk%basis(ll,mi,1,jquad)*half)
        ;                      ; wk%cfmn(ii,jquad) = wk%cfmn(ii,jquad) - Aze(lvol,idx,ii)%s(ll) * ( wk%basis(ll,mi,1,jquad)*half)
        ;                      ; wk%odmn(ii,jquad) = wk%odmn(ii,jquad) - Ate(lvol,idx,ii)%s(ll) * ( wk%basis(ll,mi,0,jquad)) * ni &
                                                    - Aze(lvol,idx,ii)%s(ll) * ( wk%basis(ll,mi,0,jquad)) * mi
        ;                      ; wk%ijreal(ii,jquad) = wk%ijreal(ii,jquad) + Ate(lvol,idx,ii)%s(ll) * wk%basis(ll,mi,0,jquad)
        ;                      ; wk%jkreal(ii,jquad) = wk%jkreal(ii,jquad) + Aze(lvol,idx,ii)%s(ll) * wk%basis(ll,mi,0,jquad)
        if( NOTstellsym ) then ; wk%ofmn(ii,jquad) = wk%ofmn(ii,jquad) + Ato(lvol,idx,ii)%s(ll) * ( wk%basis(ll,mi,1,jquad)*half)
          ;                    ; wk%sfmn(ii,jquad) = wk%sfmn(ii,jquad) - Azo(lvol,idx,ii)%s(ll) * ( wk%basis(ll,mi,1,jquad)*half)
          ;                    ; wk%evmn(ii,jquad) = wk%evmn(ii,jquad) + Ato(lvol,idx,ii)%s(ll) * ( wk%basis(ll,mi,0,jquad)) * ni &
                                                    + Azo(lvol,idx,ii)%s(ll) * ( wk%basis(ll,mi,0,jquad)) * mi
          ;                    ; wk%jireal(ii,jquad) = wk%jireal(ii,jquad) + Ato(lvol,idx,ii)%s(ll) * wk%basis(ll,mi,0,jquad)
          ;                    ; wk%kjreal(ii,jquad) = wk%kjreal(ii,jquad) + Azo(lvol,idx,ii)%s(ll) * wk%basis(ll,mi,0,jquad)
        endif
        enddo ! end of do ll; 20 Feb 13;
      else
        do ll = 0, lrad ! loop over Chebyshev polynomials; Lrad is the radial resolution;
        ;                      ; wk%efmn(ii,jquad) = wk%efmn(ii,jquad) + Ate(lvol,idx,ii)%s(ll) * ( wk%basis(ll,0,1,jquad))
        ;                      ; wk%cfmn(ii,jquad) = wk%cfmn(ii,jquad) - Aze(lvol,idx,ii)%s(ll) * ( wk%basis(ll,0,1,jquad))
        ;                      ; wk%odmn(ii,jquad) = wk%odmn(ii,jquad) - Ate(lvol,idx,ii)%s(ll) * ( wk%basis(ll,0,0,jquad)) * ni &
                                                    - Aze(lvol,idx,ii)%s(ll) * ( wk%basis(ll,0,0,jquad)) * mi
        ;                      ; wk%ijreal(ii,jquad) = wk%ijreal(ii,jquad) + Ate(lvol,idx,ii)%s(ll) * wk%basis(ll,0,0,jquad)
        ;                      ; wk%jkreal(ii,jquad) = wk%jkreal(ii,jquad) + Aze(lvol,idx,ii)%s(ll) * wk%basis(ll,0,0,jquad)
        if( NOTstellsym ) then ; wk%ofmn(ii,jquad) = wk%ofmn(ii,jquad) + Ato(lvol,idx,ii)%s(ll) * ( wk%basis(ll,0,1,jquad))
          ;                    ; wk%sfmn(ii,jquad) = wk%sfmn(ii,jquad) - Azo(lvol,idx,ii)%s(ll) * ( wk%basis(ll,0,1,jquad))
          ;                    ; wk%evmn(ii,jquad) = wk%evmn(ii,jquad) + Ato(lvol,idx,ii)%s(ll) * ( wk%basis(ll,0,0,jquad)) * ni &
                                                    + Azo(lvol,idx,ii)%s(ll) * ( wk%basis(ll,0,0,jquad)) * mi
          ;                    ; wk%jireal(ii,jquad) = wk%jireal(ii,jquad) + Ato(lvol,idx,ii)%s(ll) * wk%basis(ll,0,0,jquad)
          ;                    ; wk%kjreal(ii,jquad) = wk%kjreal(ii,jquad) + Azo(lvol,idx,ii)%s(ll) * wk%basis(ll,0,0,jquad)
        endif
        enddo ! end of do ll; 20 Feb 13;
      end if ! Lcoordinatesingularity; 01 Jul 19
    enddo ! end of do ii; 20 Feb 13;

    call invfft( mn, im, in, wk%efmn(1:mn,jquad), wk%ofmn(1:mn,jquad), wk%cfmn(1:mn,jquad), wk%sfmn(1:mn,jquad), Nt, Nz, wk%gBupper(1:Ntz,3,jquad), wk%gBupper(1:Ntz,2,jquad) )
    call invfft( mn, im, in, wk%evmn(1:mn,jquad), wk%odmn(1:mn,jquad), wk%cfmn(1:mn,jquad), wk%sfmn(1:mn,jquad), Nt, Nz, wk%gBupper(1:Ntz,1,jquad), wk%gBupper(1:Ntz,2,jquad) )

    do ii = 1, 3
      do jj = 1, 3
        wk%Blower(:,ii,jquad) = wk%Blower(:,ii,jquad) + wk%gBupper(:,jj,jquad) * guvijsave(1:Ntz,jj,ii,jquad)
      enddo
    enddo

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    ifail = 0

    call tfft( Nt, Nz, wk%Blower(1:Ntz,2,jquad), wk%Blower(1:Ntz,3,jquad), &
               mn, im(1:mn), in(1:mn), wk%Bloweremn(1:mn,2,jquad), wk%Bloweromn(1:mn,2,jquad), wk%Bloweremn(1:mn,3,jquad), wk%Bloweromn(1:mn,3,jquad), ifail )
    call tfft( Nt, Nz, wk%Blower(1:Ntz,1,jquad), wk%Blower(1:Ntz,3,jquad), &
               mn, im(1:mn), in(1:mn), wk%Bloweromn(1:mn,1,jquad), wk%Bloweremn(1:mn,1,jquad), wk%Bloweremn(1:mn,3,jquad), wk%Bloweromn(1:mn,3,jquad), ifail )
    wk%Bloweremn(1,1  ,jquad) = zero
    wk%Bloweromn(1,2:3,jquad) = zero

  end do
!$OMP END PARALLEL DO

w(1:lquad) = gaussianweight(1:lquad,lvol)

!$OMP PARALLEL DO SHARED(mn,lrad,w) PRIVATE(ii,ik,ll,ll1,bid,dfactor)
do ii = 1, mn

  if (ii==1) then ;ik = two
  else            ;ik = one
  endif

    do ll = 0, lrad

      if (Lcoordinatesingularity) then
        if (ll.lt.im(ii) .or. mod(ll+im(ii),2).ne.0) cycle
        ll1 = (ll - mod(ll,2))/2 ! shrinked dof for Zernike; 02 Jul 19
        bid = im(ii)
        dfactor = half
      else
        ll1 = ll
        bid = 0
        dfactor = one
      endif

      Tss(ll1,ii) = sum(wk%basis(ll, bid, 0, :) * wk%Bloweremn(ii,1,:) * w) * ik
      Dtc(ll1,ii) = sum(wk%basis(ll, bid, 1, :) * wk%Bloweremn(ii,2,:) * w) * dfactor * ik
      Dzc(ll1,ii) = sum(wk%basis(ll, bid, 1, :) * wk%Bloweremn(ii,3,:) * w) * dfactor * ik

      if (NOTstellsym) then
        Tsc(ll1,ii) =  sum(wk%basis(ll, bid, 0, :) * wk%Bloweromn(ii,1,:) * w) * ik
        Dts(ll1,ii) =  sum(wk%basis(ll, bid, 1, :) * wk%Bloweromn(ii,2,:) * w) * dfactor * ik
        Dzs(ll1,ii) =  sum(wk%basis(ll, bid, 1, :) * wk%Bloweromn(ii,3,:) * w) * dfactor * ik
      endif

      if (dBdX%L) cycle ! dMD matrix does not depend on geometry
      Ttc(ll1,ii) = (sum(wk%basis(ll, bid, 0, :) * wk%cfmn(ii,:) * w) + sum(wk%basis(ll, bid, 1, :) * wk%jkreal(ii,:) * w) * dfactor) * ik
      Tzc(ll1,ii) = (sum(wk%basis(ll, bid, 0, :) * wk%efmn(ii,:) * w) - sum(wk%basis(ll, bid, 1, :) * wk%ijreal(ii,:) * w) * dfactor) * ik

      if (NOTstellsym) then
        Tts(ll1,ii) = (sum(wk%basis(ll, bid, 0, :) * wk%sfmn(ii,:) * w) + sum(wk%basis(ll, bid, 1, :) * wk%kjreal(ii,:) * w) * dfactor) * ik
        Tzs(ll1,ii) = (sum(wk%basis(ll, bid, 0, :) * wk%ofmn(ii,:) * w) - sum(wk%basis(ll, bid, 1, :) * wk%jireal(ii,:) * w) * dfactor) * ik
      endif

    enddo !ll
  enddo !ii
!$OMP END PARALLEL DO

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  Tss = Tss * pi2pi2nfphalf
  Dtc = Dtc * pi2pi2nfphalf
  Dzc = Dzc * pi2pi2nfphalf

  if (.not.dBdX%L) then
    Ttc = Ttc * pi2pi2nfphalf
    Tzc = Tzc * pi2pi2nfphalf
  endif

  if (NOTstellsym) then

    Tsc = Tsc * pi2pi2nfphalf
    Dts = Dts * pi2pi2nfphalf
    Dzs = Dzs * pi2pi2nfphalf

    if (.not.dBdX%L) then
      Tts = Tts * pi2pi2nfphalf
      Tzs = Tzs * pi2pi2nfphalf
    endif

  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


9999 continue
  cput = MPI_WTIME()
  Tintghs  = Tintghs  + ( cput-cpuo )
  return


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine intghs

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> \brief init workspace
!>
!> @param lvol
subroutine intghs_workspace_init(lvol)
  use mod_kinds, only: wp => dp
  use constants, only : zero
  use inputlist, only : Mpol, Lrad, Wmacros, Wintghs
  use fileunits, only : ounit
  use cputiming, only : Tintghs
  use allglobal, only : Ntz, mn, Iquad, myid, ncpu, cpus, MPI_COMM_SPEC
  use intghs_module


#ifdef OPENMP
  USE OMP_LIB
#endif
  use mpi
  implicit none
  integer   :: ierr, astat, ios, nthreads, ithread
  real(wp)      :: cput, cpui, cpuo=0 ! cpu time; cpu initial; cpu old; 31 Jan 13;


  integer, INTENT(IN) :: lvol
  integer             :: lquad


  cpui = MPI_WTIME()
  cpuo = cpui
#ifdef OPENMP
  nthreads = omp_get_max_threads()
#else
  nthreads = 1
#endif


  lquad = Iquad(lvol)


   allocate( wk%gBupper(1:Ntz,3,lquad), stat=astat )
   wk%gBupper(1:Ntz,3,lquad) = zero


   allocate( wk%Blower(1:Ntz,3,lquad), stat=astat )
   wk%Blower(1:Ntz,3,lquad) = zero


   allocate( wk%Bloweremn(1:mn,3,lquad), stat=astat )
   wk%Bloweremn(1:mn,3,lquad) = zero


   allocate( wk%Bloweromn(1:mn,3,lquad), stat=astat )
   wk%Bloweromn(1:mn,3,lquad) = zero


   allocate( wk%efmn(1:mn,lquad), stat=astat )
   wk%efmn(1:mn,lquad) = zero


   allocate( wk%ofmn(1:mn,lquad), stat=astat )
   wk%ofmn(1:mn,lquad) = zero


   allocate( wk%evmn(1:mn,lquad), stat=astat )
   wk%evmn(1:mn,lquad) = zero


   allocate( wk%odmn(1:mn,lquad), stat=astat )
   wk%odmn(1:mn,lquad) = zero


   allocate( wk%cfmn(1:mn,lquad), stat=astat )
   wk%cfmn(1:mn,lquad) = zero


   allocate( wk%sfmn(1:mn,lquad), stat=astat )
   wk%sfmn(1:mn,lquad) = zero


   allocate( wk%ijreal(1:mn,lquad), stat=astat )
   wk%ijreal(1:mn,lquad) = zero


   allocate( wk%jkreal(1:mn,lquad), stat=astat )
   wk%jkreal(1:mn,lquad) = zero


   allocate( wk%jireal(1:mn,lquad), stat=astat )
   wk%jireal(1:mn,lquad) = zero


   allocate( wk%kjreal(1:mn,lquad), stat=astat )
   wk%kjreal(1:mn,lquad) = zero


   allocate( wk%basis(0:Lrad(lvol),0:mpol,0:1, lquad), stat=astat )
   wk%basis(0:Lrad(lvol),0:mpol,0:1, lquad) = zero



9999 continue
  cput = MPI_WTIME()
  Tintghs = Tintghs + ( cput-cpuo )
  return


end subroutine intghs_workspace_init

!> \brief free workspace
!>
!> @param lvol
subroutine intghs_workspace_destroy()
  use mod_kinds, only: wp => dp
  use inputlist, only : Wmacros, Wintghs
  use fileunits, only : ounit
  use cputiming, only : Tintghs
  use allglobal, only : myid, ncpu, cpus, MPI_COMM_SPEC
  use intghs_module


#ifdef OPENMP
  USE OMP_LIB
#endif
  use mpi
  implicit none
  integer   :: ierr, astat, ios, nthreads, ithread
  real(wp)      :: cput, cpui, cpuo=0 ! cpu time; cpu initial; cpu old; 31 Jan 13;



  cpui = MPI_WTIME()
  cpuo = cpui
#ifdef OPENMP
  nthreads = omp_get_max_threads()
#else
  nthreads = 1
#endif



   deallocate(wk%gBupper,stat=astat)


   deallocate(wk%Blower,stat=astat)


   deallocate(wk%Bloweremn,stat=astat)


   deallocate(wk%Bloweromn,stat=astat)


   deallocate(wk%efmn,stat=astat)


   deallocate(wk%ofmn,stat=astat)


   deallocate(wk%evmn,stat=astat)


   deallocate(wk%odmn,stat=astat)


   deallocate(wk%cfmn,stat=astat)


   deallocate(wk%sfmn,stat=astat)


   deallocate(wk%ijreal,stat=astat)


   deallocate(wk%jkreal,stat=astat)


   deallocate(wk%jireal,stat=astat)


   deallocate(wk%kjreal,stat=astat)


   deallocate(wk%basis,stat=astat)



9999 continue
  cput = MPI_WTIME()
  Tintghs = Tintghs + ( cput-cpuo )
  return


end subroutine intghs_workspace_destroy
