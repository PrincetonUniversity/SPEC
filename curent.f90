!> \defgroup grp_currents Plasma Currents
!>
!> \file
!> \brief Computes the plasma current, \f$I \equiv \int B_\theta \, d\theta\f$, and the "linking" current, \f$G \equiv \int B_\zeta \, d\zeta\f$.

!> \brief Computes the plasma current, \f$I \equiv \int B_\theta \, d\theta\f$, and the "linking" current, \f$G \equiv \int B_\zeta \, d\zeta\f$.
!> \ingroup grp_currents
!>
!> **enclosed currents**
!>
!> <ul>
!> <li>  In the vacuum region, the enclosed currents are given by either surface integrals of the current density or line integrals of the magnetic field,
!>       \f{eqnarray}{ \int_{\cal S} {\bf j}\cdot d{\bf s} = \int_{\partial \cal S} {\bf B}\cdot d{\bf l},
!>       \f}
!>       and line integrals are usually easier to compute than surface integrals. </li>
!> <li>  The magnetic field is given by the curl of the magnetic vector potential, as described in e.g. bfield() . </li>
!> <li>  The toroidal, plasma current is obtained by taking a "poloidal" loop, \f$d{\bf l}={\bf e}_\theta \, d\theta\f$,
!>       on the plasma boundary, where \f$B^s=0\f$, to obtain
!>       \f{eqnarray}{ I   \equiv   \int_{0}^{2\pi} \! {\bf B} \cdot {\bf e}_\theta \, d\theta
!>                            =     \int_{0}^{2\pi} \! \left( - \partial_s A_\zeta \; \bar g_{\theta\theta} + \partial_s A_\theta \; \bar g_{\theta\zeta} \right) \, d\theta,
!>       \f}
!>       where \f$\bar g_{\mu\nu} \equiv g_{\mu\nu} / \sqrt g\f$. </li>
!> <li>  The poloidal, "linking" current through the torus is obtained by taking a "toroidal" loop, \f$d{\bf l}={\bf e}_\zeta \, d\zeta\f$,
!>       on the plasma boundary to obtain
!>       \f{eqnarray}{ G \equiv \int_{0}^{2\pi} \! {\bf B}\cdot {\bf e}_\zeta \, d\zeta
!>                        =     \int_{0}^{2\pi} \! \left( - \partial_s A_\zeta \; \bar g_{\theta\zeta} + \partial_s A_\theta \; \bar g_{\zeta\zeta} \right) \, d\zeta.
!>       \f} </li>
!> </ul>
!>
!> **Fourier integration**
!>
!> <ul>
!> <li>  Using \f$f \equiv - \partial_s A_\zeta \; \bar g_{\theta\theta} + \partial_s A_\theta \; \bar g_{\theta\zeta}\f$,
!>       the integral for the plasma current is
!>       \f{eqnarray}{ I = {\sum_i}^\prime f_i \cos(n_i \zeta) 2\pi, \label{eq:plasmacurrent_curent}
!>       \f}
!>       where \f$\sum^\prime\f$ includes only the \f$m_i=0\f$ harmonics. </li>
!> <li>  Using \f$g \equiv - \partial_s A_\zeta \; \bar g_{\theta\zeta} + \partial_s A_\theta \; \bar g_{\zeta\zeta}\f$,
!>       the integral for the linking current is
!>       \f{eqnarray}{ G = {\sum_i}^\prime g_i \cos(m_i \zeta) 2\pi, \label{eq:linkingcurrent_curent}
!>       \f}
!>       where \f$\sum^\prime\f$ includes only the \f$n_i=0\f$ harmonics. </li>
!> <li>  The plasma  current, Eqn.\f$(\ref{eq:plasmacurrent_curent})\f$, should be independent of \f$\zeta\f$,
!>       and the linking current, Eqn.\f$(\ref{eq:linkingcurrent_curent})\f$, should be independent of \f$\theta\f$.
!>       \todo Perhaps this can be proved analytically; in any case it should be confirmed numerically.
!>
!> </li>
!> </ul>
!>
!> @param[in] lvol index of volume
!> @param[in] mn number of Fourier harmonics
!> @param[in] Nt number of grid points along \f$\theta\f$
!> @param[in] Nz number of grid points along \f$\zeta\f$
!> @param[in] iflag some integer flag
!> @param[out] ldItGp plasma and linking current
subroutine curent( lvol, mn, Nt, Nz, iflag, ldItGp )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  use constants, only : zero, one, two, pi2

  use numerical, only :

  use fileunits, only : ounit

  use inputlist, only : Wmacros, Wcurent, Lrad, Lconstraint

  use cputiming, only : Tcurent

  use allglobal, only : ncpu, cpus, myid, MPI_COMM_SPEC, &
                        Mvol, im, in, mne, ime, ine, &
                        YESstellsym, NOTstellsym, &
                        sg, guvij, &
                        Ntz, ijreal, ijimag, jireal, jiimag, &
                        efmn, ofmn, cfmn, sfmn, evmn, odmn, comn, simn, &
                        Ate, Aze, Ato, Azo, TT, &
                        build_vector_potential

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS

  INTEGER, intent(in)  :: lvol, mn, Nt, Nz, iflag
  REAL   , intent(out) :: ldItGp(0:1,-1:2)

  INTEGER              :: innout, ideriv, ii, ll, Lcurvature, ifail
  REAL                 :: lss
  REAL                 :: Bsupt(1:Nt*Nz,-1:2), Bsupz(1:Nt*Nz,-1:2), Bsups(1:Nt*Nz,-1:2), Bsups_2(1:Nt*Nz,-1:2)
  
  BEGIN(curent)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef DEBUG
  FATAL( curent, lvol.ne.Mvol, this is only defined in the vacuum region )
#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

 !ldItGp(0:1,-1:2) = zero ! initialize intent out; 12 Sep 16;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  innout = 0 ; lss = two * innout - one ! lvol = Mvol ! internal shorthand/variables; 08 Feb 16; note that lvol = Mvol is hardwired; 12 Sep 16;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if (lconstraint .eq. -2) then
    innout = 1.0
    lss = 1.0
  end if

  do ideriv = -1, 2 ! labels derivative of magnetic field wrt enclosed fluxes; 20 Apr 13;

   if( iflag.eq. 1 .and. ideriv.ne.0 ) cycle ! derivatives of currents                                                   are not required; 20 Jun 14;
   if( iflag.eq. 2 .and. ideriv.lt.0 ) cycle ! derivatives of currents  wrt geometry                                     is  not required; 20 Jun 14;
   if( iflag.eq.-1 .and. ideriv.gt.0 ) cycle ! derivatives of currents  wrt enclosed toroidal and enclosed poloidal flux are not required; 20 Jun 14;

   call build_vector_potential(lvol, innout, ideriv, 1)

   call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), &
                Nt, Nz, Bsupz(1:Ntz,ideriv), Bsupt(1:Ntz,ideriv) ) ! map to real space;

   if (Lconstraint .eq. -2) then
        call build_vector_potential(lvol, innout, ideriv, 0)

        do ii = 1, mn ! loop over Fourier harmonics; 17 May 21;
            efmn(ii) = -im(ii)*efmn(ii) ! theta derivative of Aze
            cfmn(ii) = -in(ii)*cfmn(ii) ! zeta  derivative of Ate
        enddo ! end of do ii; 

        call invfft( mn, im(1:mn), in(1:mn), ofmn(1:mn), efmn(1:mn), -sfmn(1:mn), cfmn(1:mn), &
                Nt, Nz, Bsups(1:Ntz,ideriv), Bsups_2(1:Ntz,ideriv))
        Bsups(1:Ntz,ideriv) = Bsups(1:Ntz,ideriv) + Bsups_2(1:ntz,ideriv)
   else
        Bsups = zero
   end if

  enddo ! end of do ideriv; 31 Jan 13;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  select case( iflag )
  case( -1 ) ; Lcurvature = 3 ! will need derivatives of metrics w.r.t. interface deformation; 15 Sep 16;
  case(  1 ) ; Lcurvature = 1
  case(  2 ) ; Lcurvature = 1
  end select

  WCALL( curent, coords,( lvol, lss, Lcurvature, Ntz, mn ) ) ! get "lower" metric elements evaluated on innout interface;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  do ideriv = -1, 2 ! labels derivative of magnetic field wrt enclosed fluxes; 20 Apr 13;

   if( iflag.eq. 1 .and. ideriv.ne.0 ) cycle ! derivatives of currents                                                   are not required; 20 Jun 14;
   if( iflag.eq. 2 .and. ideriv.lt.0 ) cycle ! derivatives of currents  wrt geometry                                     is  not required; 20 Jun 14;
   if( iflag.eq.-1 .and. ideriv.gt.0 ) cycle ! derivatives of currents  wrt enclosed toroidal and enclosed poloidal flux are not required; 20 Jun 14;

   if (Lconstraint .eq. -2) then
   ijreal(1:Ntz) =  (Bsups(1:Ntz,ideriv) * guvij(1:Ntz,2,1,0) - Bsupt(1:Ntz,ideriv) * guvij(1:Ntz,2,2,0) + Bsupz(1:Ntz,ideriv) * guvij(1:Ntz,2,3,0) ) / sg(1:Ntz,0)
   ijimag(1:Ntz) =  (Bsups(1:Ntz,ideriv) * guvij(1:Ntz,1,3,0) - Bsupt(1:Ntz,ideriv) * guvij(1:Ntz,2,3,0) + Bsupz(1:Ntz,ideriv) * guvij(1:Ntz,3,3,0) ) / sg(1:Ntz,0)
   else
   ijreal(1:Ntz) =  ( - Bsupt(1:Ntz,ideriv) * guvij(1:Ntz,2,2,0) + Bsupz(1:Ntz,ideriv) * guvij(1:Ntz,2,3,0) ) / sg(1:Ntz,0)
   ijimag(1:Ntz) =  ( - Bsupt(1:Ntz,ideriv) * guvij(1:Ntz,2,3,0) + Bsupz(1:Ntz,ideriv) * guvij(1:Ntz,3,3,0) ) / sg(1:Ntz,0)
   end if 

   if( ideriv.eq.-1 ) then ! add derivatives of metrics with respect to interface geometry; 15 Sep 16;
   ijreal(1:Ntz) = ijreal(1:Ntz) + ( - Bsupt(1:Ntz,     0) * guvij(1:Ntz,2,2,1) + Bsupz(1:Ntz,     0) * guvij(1:Ntz,2,3,1) ) / sg(1:Ntz,0)
   ijimag(1:Ntz) = ijimag(1:Ntz) + ( - Bsupt(1:Ntz,     0) * guvij(1:Ntz,2,3,1) + Bsupz(1:Ntz,     0) * guvij(1:Ntz,3,3,1) ) / sg(1:Ntz,0)
   endif

   ifail = 0
   call tfft( Nt, Nz, ijreal(1:Ntz), ijimag(1:Ntz), &
              mne, ime(1:mne), ine(1:mne), efmn(1:mne), ofmn(1:mne), cfmn(1:mne), sfmn(1:mne), ifail )

   ldItGp(0,ideriv) = efmn(1) * pi2 ! "toroidal", plasma  current; 12 Sep 16;
   ldItGp(1,ideriv) = cfmn(1) * pi2 ! "poloidal", linking current; 12 Sep 16;

#ifdef DEBUG
   if( Wcurent ) then
    write(ounit,'("curent : ", 10x ," : myid=",i3," ; dItGp(0:1,",i2,")=",es23.15,",",es23.15," ;")') myid, ideriv, ldItGp(0:1,ideriv)
   endif
#endif

  enddo ! end of do ideriv; 31 Jan 13;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!#ifdef DEBUG
!
!  if( Wcurent ) then
!
!   select case( ideriv )
!   case( 0 )
!    do ii = 1, mne
!     if( ime(ii).eq.0 ) write(ounit,1000) myid, ideriv, ime(ii), ine(ii), efmn(ii)
!    enddo
!    do ii = 1, mne
!     if( ine(ii).eq.0 ) write(ounit,1010) myid, ideriv, ime(ii), ine(ii), cfmn(ii)
!    enddo
!   end select
!
!  endif ! end of if( Wcurent ) ; 05 Feb 16;
!
!1000 format("curent : " 10x " : myid="i3" ; ideriv ="i2" ; ("i3","i3" ) : I_n ="es13.5" ;")
!1010 format("curent : " 10x " : myid="i3" ; ideriv ="i2" ; ("i3","i3" ) : G_m ="es13.5" ;")
!
!#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(curent)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine curent

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
