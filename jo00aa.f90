!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (diagnostic) ! Measures error in Beltrami equation, $\nabla \times {\bf B} - \mu {\bf B}$.

!latex \briefly{Measures error in Beltrami field, $||\nabla\times{\bf B}-\mu {\bf B}||$.}

!latex \calledby{\link{xspech}}
!latex \calls{\link{coords}}

!latex \tableofcontents

!latex \subsection{overview}

!latex This routine is called by \link{xspech} as a post diagnostic and only if \inputvar{Lcheck} $= 1$.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{construction of current, ${\bf j} \equiv \nabla \times \nabla \times {\bf A}$}

!latex \begin{enumerate}
!latex \item \Ais
!latex \item \Bis
!latex \item The current is
!latex       \be \sqrt g \, {\bf j} = ( \partial_\t B_\z - \partial_\z B_\t) \; {\bf e}_\s
!latex                           + ( \partial_\z B_\s - \partial_\s B_\z) \; {\bf e}_\t
!latex                           + ( \partial_\s B_\t - \partial_\t B_\s) \; {\bf e}_\z,
!latex       \ee
!latex       where (for computational convenience) the covariant components of ${\bf B}$ are computed as
!latex       \be
!latex       B_\s & = & (\sqrt g B^\s) \, g_{\s\s} / \sqrt g + (\sqrt g B^\t) \, g_{\s\t} / \sqrt g + (\sqrt g B^\z) \, g_{\s\z} / \sqrt g, \\
!latex       B_\t & = & (\sqrt g B^\s) \, g_{\s\t} / \sqrt g + (\sqrt g B^\t) \, g_{\t\t} / \sqrt g + (\sqrt g B^\z) \, g_{\t\z} / \sqrt g, \\
!latex       B_\z & = & (\sqrt g B^\s) \, g_{\s\z} / \sqrt g + (\sqrt g B^\t) \, g_{\t\z} / \sqrt g + (\sqrt g B^\z) \, g_{\z\z} / \sqrt g.
!latex       \ee
!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{quantification of the error}

!latex \begin{enumerate}
!latex \item The measures of the error are
!latex       \be ||\left( {\bf j}-\mu {\bf B}\right)\cdot\nabla \s || &  \equiv  &
!latex         \int \!\! ds \ooint \left| \sqrt g \, {\bf j}\cdot\nabla \s-\mu \; \sqrt g \, {\bf B} \cdot\nabla \s \right|, \label{eq:Es} \\
!latex           ||\left( {\bf j}-\mu {\bf B}\right)\cdot\nabla \t || &  \equiv  &
!latex         \int \!\! ds \ooint \left| \sqrt g \, {\bf j}\cdot\nabla \t-\mu \; \sqrt g \, {\bf B} \cdot\nabla \t \right|, \label{eq:Et} \\
!latex           ||\left( {\bf j}-\mu {\bf B}\right)\cdot\nabla \z || &  \equiv  &
!latex         \int \!\! ds \ooint \left| \sqrt g \, {\bf j}\cdot\nabla \z-\mu \; \sqrt g \, {\bf B} \cdot\nabla \z \right|. \label{eq:Ez}
!latex       \ee
!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{details of the numerics}

!latex \begin{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine jo00aa( lvol, Ntz, lquad, mn )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  use constants, only : zero, half, one, two, pi2

  use fileunits, only : ounit

  use inputlist, only : Wmacros, Wjo00aa, Nvol, Lrad, mu, mpol, Igeometry, Nfp, Lerrortype

  use cputiming, only : Tjo00aa

  use allglobal, only : myid, cpus, MPI_COMM_SPEC, ext, ivol, &
                        im, in, regumm, &
                        Mvol, &
                        cheby, zernike, &
                        Ate, Aze, Ato, Azo, &
                        sg, guvij, Rij, Zij, &
                        Nt, Nz, efmn, ofmn, cfmn, sfmn, &
                        NOTstellsym, &
                        Lcoordinatesingularity, &
                        beltramierror, Rij, Zij, gBzeta, Node, &
                        pi2nfp, ivol, RTT, TT, dtflux, dpflux

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS

  INTEGER, intent(in) :: lvol, Ntz, lquad, mn ! these are really global, but are included in argument list to remove allocations;

  INTEGER             :: jquad, Lcurvature, ll, ii, jj, kk, uu, ideriv, twolquad, mm, jk

  REAL                :: lss, sbar, sbarhim(0:2), gBu(1:Ntz,1:3,0:3), gJu(1:Ntz,1:3), jerror(1:3), jerrormax(1:3), intvol
  REAL                :: B_cartesian(1:Ntz,1:3), J_cartesian(1:Ntz,1:3)

  REAL                :: Atemn(1:mn,0:2), Azemn(1:mn,0:2), Atomn(1:mn,0:2), Azomn(1:mn,0:2)

  INTEGER             :: itype, icdgqf
  REAL                :: aa, bb, cc, dd, weight(1:lquad+1), abscis(1:lquad), workfield(1:2*lquad)

  REAL                :: zeta, teta, st(2), Bst(2)

  BEGIN(jo00aa)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef DEBUG
  FATAL( jo00aa, lquad.lt.1, invalid lquad supplied to jo00aa )
  FATAL( jo00aa, lvol.lt.1 .or. lvol.gt.Mvol, invalid interface label )
#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item The integration over $s$ is performed using Gaussian integration, e.g., $\ds \int \!\! f(s) ds \approx \sum_k \omega_k f(s_k)$;
!latex       with the abscissae, $s_k$, and the weights, $\omega_k$, for $k=1,$ \internal{Iquad$_v$},
!latex       determined by \verb+CDGQF+.
!latex       The resolution, \type{N} $ \equiv $ \internal{Iquad$_v$}, is determined by \inputvar{Nquad} (see \link{global} and \link{preset}).
!latex       A fatal error is enforced by \link{jo00aa}
!latex       if \verb+CDGQF+ returns an \type{ifail} $\ne 0$.

  itype = 1 ; aa = -one ; bb = +one ; cc = zero ; dd = zero ; twolquad = 2 * lquad

  call CDGQF( lquad, abscis(1:lquad), weight(1:lquad), itype, aa, bb, twolquad, workfield(1:twolquad), icdgqf ) ! prepare Gaussian quadrature;
  weight(lquad+1) = zero

! write(ounit,'("jo00aa :  WARNING ! : THE ERROR FLAGS RETURNED BY CDGQF SEEM DIFFERENT TO NAG:D01BCF (may be trivial, but please revise); 2018/01/10;")')

  cput= GETTIME
  select case( icdgqf ) !                                                         123456789012345
  case( 0 )    ;  if( Wjo00aa ) write(ounit,1000) cput-cpus, myid, lvol, icdgqf, "success        ", abscis(1:lquad)
  case( 1 )    ;                write(ounit,1000) cput-cpus, myid, lvol, icdgqf, "failed         ", abscis(1:lquad)
  case( 2 )    ;                write(ounit,1000) cput-cpus, myid, lvol, icdgqf, "input error    ", abscis(1:lquad)
  case( 3 )    ;                write(ounit,1000) cput-cpus, myid, lvol, icdgqf, "input error    ", abscis(1:lquad)
  case( 4 )    ;                write(ounit,1000) cput-cpus, myid, lvol, icdgqf, "weight overflow", abscis(1:lquad)
  case( 5 )    ;                write(ounit,1000) cput-cpus, myid, lvol, icdgqf, "weight zero    ", abscis(1:lquad)
  case( 6 )    ;                write(ounit,1000) cput-cpus, myid, lvol, icdgqf, "failed         ", abscis(1:lquad)
  case default ;                write(ounit,1000) cput-cpus, myid, lvol, icdgqf, "weird          ", abscis(1:lquad)
  end select

  if( Wjo00aa )                 write(ounit,1001)                                                   weight(1:lquad)

1000 format("jo00aa : ",f10.2," : myid=",i3," ; lvol=",i3," ; icdgqf=",i3," ; "a15" ;":" abscissae ="99f10.6)
1001 format("jo00aa : ", 10x ," :       "3x"          "3x"            "3x"    "15x" ;":" weights   ="99f10.6)

  FATAL( jo00aa, icdgqf.ne.0, failed to construct Gaussian integration abscisae and weights )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  jerror(1:3) = zero ; ideriv = 0 ! three components of the error in \curl B - mu B; initialize summation;
  jerrormax(1:3) = zero; intvol = zero

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Inside the Gaussian quadrature loop, i.e. for each $s_k$,

!latex \begin{enumerate}

  do jquad = 1, lquad+1 ! loop over radial sub-sub-grid (numerical quadrature);

   if (jquad.eq.lquad+1) then
    lss = one           ; sbar = one
   else
    lss = abscis(jquad) ; sbar = ( lss + one ) * half
   endif

!latex \item The metric elements, $g_{\mu,\nu} \equiv $ \type{gij(1:6,0,1:Ntz)}, and the Jacobian, $\sqrt g \equiv $ \type{sg(0,1:Ntz)},
!latex       are calculated on a regular angular grid, $(\t_i,\z_j)$, in \link{coords}.
!latex       The derivatives $\partial_i g_{\mu,\nu} \equiv$ \type{gij(1:6,i,1:Ntz)} and $\partial_i \sqrt g \equiv $ \type{sg(i,1:Ntz)},
!latex       with respect to $ i \in \{ \s,\t,\z \}$ are also returned.

   Lcurvature = 2

   WCALL( jo00aa, coords, ( lvol, lss, Lcurvature, Ntz, mn ) ) ! returns coordinates, metrics, . . .

   if (Lcoordinatesingularity) then ! Zernike 1 Jul 2019
     call get_zernike_d2(sbar, Lrad(lvol), mpol, zernike)
   else
     call get_cheby_d2(lss, Lrad(lvol), cheby(0:Lrad(lvol),0:2))
   endif

   Atemn(1:mn,0:2) = zero ! initialize summation over Chebyshev/Zernike polynomials;
   Azemn(1:mn,0:2) = zero
   if( NOTstellsym ) then
   Atomn(1:mn,0:2) = zero
   Azomn(1:mn,0:2) = zero
   else
   Atomn(1:mn,0:2) = zero ! these are used below;
   Azomn(1:mn,0:2) = zero
   endif

  !latex \item The Fourier components of the vector potential given in \Eqn{At} and \Eqn{Az}, and their first and second radial derivatives, are summed.

   if (Lcoordinatesingularity) then
    do ll = 0, Lrad(lvol) ! radial (Chebyshev) resolution of magnetic vector potential;

      do ii = 1, mn  ! Fourier resolution of magnetic vector potential;
      mm = im(ii)
      if (ll < mm) cycle
      if (mod(ll+mm, 2) > 0) cycle

      ;Atemn(ii,0) = Atemn(ii,0) + Ate(lvol,ideriv,ii)%s(ll) * zernike(ll,mm,0)
      ;Atemn(ii,1) = Atemn(ii,1) + Ate(lvol,ideriv,ii)%s(ll) * zernike(ll,mm,1) * half
      ;Atemn(ii,2) = Atemn(ii,2) + Ate(lvol,ideriv,ii)%s(ll) * zernike(ll,mm,2) * half * half

      ;Azemn(ii,0) = Azemn(ii,0) + Aze(lvol,ideriv,ii)%s(ll) * zernike(ll,mm,0)
      ;Azemn(ii,1) = Azemn(ii,1) + Aze(lvol,ideriv,ii)%s(ll) * zernike(ll,mm,1) * half
      ;Azemn(ii,2) = Azemn(ii,2) + Aze(lvol,ideriv,ii)%s(ll) * zernike(ll,mm,2) * half * half

      if( NOTstellsym ) then

        Atomn(ii,0) = Atomn(ii,0) + Ato(lvol,ideriv,ii)%s(ll) * zernike(ll,mm,0)
        Atomn(ii,1) = Atomn(ii,1) + Ato(lvol,ideriv,ii)%s(ll) * zernike(ll,mm,1) * half
        Atomn(ii,2) = Atomn(ii,2) + Ato(lvol,ideriv,ii)%s(ll) * zernike(ll,mm,2) * half * half

        Azomn(ii,0) = Azomn(ii,0) + Azo(lvol,ideriv,ii)%s(ll) * zernike(ll,mm,0)
        Azomn(ii,1) = Azomn(ii,1) + Azo(lvol,ideriv,ii)%s(ll) * zernike(ll,mm,1) * half
        Azomn(ii,2) = Azomn(ii,2) + Azo(lvol,ideriv,ii)%s(ll) * zernike(ll,mm,2) * half * half

      endif ! end of if( NOTstellsym) ; 20 Jun 14;

      enddo ! end of do ii;

    enddo ! end of do ll;

   else

    do ll = 0, Lrad(lvol) ! radial (Chebyshev) resolution of magnetic vector potential;

      do ii = 1, mn  ! Fourier resolution of magnetic vector potential;

      ;Atemn(ii,0) = Atemn(ii,0) + Ate(lvol,ideriv,ii)%s(ll) * cheby(ll,0)
      ;Atemn(ii,1) = Atemn(ii,1) + Ate(lvol,ideriv,ii)%s(ll) * cheby(ll,1)
      ;Atemn(ii,2) = Atemn(ii,2) + Ate(lvol,ideriv,ii)%s(ll) * cheby(ll,2)

      ;Azemn(ii,0) = Azemn(ii,0) + Aze(lvol,ideriv,ii)%s(ll) * cheby(ll,0)
      ;Azemn(ii,1) = Azemn(ii,1) + Aze(lvol,ideriv,ii)%s(ll) * cheby(ll,1)
      ;Azemn(ii,2) = Azemn(ii,2) + Aze(lvol,ideriv,ii)%s(ll) * cheby(ll,2)

      if( NOTstellsym ) then

        Atomn(ii,0) = Atomn(ii,0) + Ato(lvol,ideriv,ii)%s(ll) * cheby(ll,0)
        Atomn(ii,1) = Atomn(ii,1) + Ato(lvol,ideriv,ii)%s(ll) * cheby(ll,1)
        Atomn(ii,2) = Atomn(ii,2) + Ato(lvol,ideriv,ii)%s(ll) * cheby(ll,2)

        Azomn(ii,0) = Azomn(ii,0) + Azo(lvol,ideriv,ii)%s(ll) * cheby(ll,0)
        Azomn(ii,1) = Azomn(ii,1) + Azo(lvol,ideriv,ii)%s(ll) * cheby(ll,1)
        Azomn(ii,2) = Azomn(ii,2) + Azo(lvol,ideriv,ii)%s(ll) * cheby(ll,2)

      endif ! end of if( NOTstellsym) ; 20 Jun 14;

      enddo ! end of do ii;

    enddo ! end of do ll;
   end if
!latex \item The quantities $\sqrt g B^\s$, $\sqrt g B^\t$ and $\sqrt g B^\z$, and their first and second derivatives with respect to $(\s,\t,\z)$,
!latex       are computed on the regular angular grid (using FFTs).

   ofmn(1:mn) = -          im(1:mn)*Azemn(1:mn,0) -          in(1:mn)*Atemn(1:mn,0)
   efmn(1:mn) = +          im(1:mn)*Azomn(1:mn,0) +          in(1:mn)*Atomn(1:mn,0)
   sfmn(1:mn) = -          im(1:mn)*Azemn(1:mn,1) -          in(1:mn)*Atemn(1:mn,1)
   cfmn(1:mn) = +          im(1:mn)*Azomn(1:mn,1) +          in(1:mn)*Atomn(1:mn,1)

   call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, gBu(1:Ntz,1,0), gBu(1:Ntz,1,1) ) !  (gB^s)   , d(gB^s)/ds;

   efmn(1:mn) = - im(1:mn)*im(1:mn)*Azemn(1:mn,0) - im(1:mn)*in(1:mn)*Atemn(1:mn,0)
   ofmn(1:mn) = - im(1:mn)*im(1:mn)*Azomn(1:mn,0) - im(1:mn)*in(1:mn)*Atomn(1:mn,0)
   cfmn(1:mn) = + in(1:mn)*im(1:mn)*Azemn(1:mn,0) + in(1:mn)*in(1:mn)*Atemn(1:mn,0)
   sfmn(1:mn) = + in(1:mn)*im(1:mn)*Azomn(1:mn,0) + in(1:mn)*in(1:mn)*Atomn(1:mn,0)

   call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, gBu(1:Ntz,1,2), gBu(1:Ntz,1,3) ) ! d(gB^s)/dt, d(gB^s)/dz;

   efmn(1:mn) =                                   -                   Azemn(1:mn,1)
   ofmn(1:mn) =                                   -                   Azomn(1:mn,1)
   cfmn(1:mn) =                                   -                   Azemn(1:mn,2)
   sfmn(1:mn) =                                   -                   Azomn(1:mn,2)

   call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, gBu(1:Ntz,2,0), gBu(1:Ntz,2,1) ) !  (gB^t)   , d(gB^t)/ds;

   ofmn(1:mn) =                                   +          im(1:mn)*Azemn(1:mn,1)
   efmn(1:mn) =                                   -          im(1:mn)*Azomn(1:mn,1)
   sfmn(1:mn) =                                   -          in(1:mn)*Azemn(1:mn,1)
   cfmn(1:mn) =                                   +          in(1:mn)*Azomn(1:mn,1)

   call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, gBu(1:Ntz,2,2), gBu(1:Ntz,2,3) ) ! d(gB^t)/dt, d(gB^t)/dz;

   efmn(1:mn) = +                   Atemn(1:mn,1)
   ofmn(1:mn) = +                   Atomn(1:mn,1)
   cfmn(1:mn) = +                   Atemn(1:mn,2)
   sfmn(1:mn) = +                   Atomn(1:mn,2)

   call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, gBu(1:Ntz,3,0), gBu(1:Ntz,3,1) ) !  (gB^z)   , d(gB^z)/ds;

   ofmn(1:mn) = -          im(1:mn)*Atemn(1:mn,1)
   efmn(1:mn) = +          im(1:mn)*Atomn(1:mn,1)
   sfmn(1:mn) = +          in(1:mn)*Atemn(1:mn,1)
   cfmn(1:mn) = -          in(1:mn)*Atomn(1:mn,1)

   call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, gBu(1:Ntz,3,2), gBu(1:Ntz,3,3) ) ! d(gB^z)/dt, d(gB^z)/dz;

!latex \item The following quantities are then computed on the regular angular grid
!latex       \be \sqrt g j^\s & = & \sum_u \left[
!latex           \partial_\t (\sqrt g B^u) \; g_{u,\z} + (\sqrt g B^u) \; \partial_\t g_{u,\z} - (\sqrt g B^u) g_{u,\z} \; \partial_\t \sqrt g / \sqrt g
!latex                                    \right] / \sqrt g \nonumber \\
!latex                       & - & \sum_u \left[
!latex           \partial_\z (\sqrt g B^u) \; g_{u,\t} + (\sqrt g B^u) \; \partial_\z g_{u,\t} - (\sqrt g B^u) g_{u,\t} \; \partial_\z \sqrt g / \sqrt g
!latex                                    \right] / \sqrt g, \\
!latex           \sqrt g j^\t & = & \sum_u \left[
!latex           \partial_\z (\sqrt g B^u) \; g_{u,\s} + (\sqrt g B^u) \; \partial_\z g_{u,\s} - (\sqrt g B^u) g_{u,\s} \; \partial_\z \sqrt g / \sqrt g
!latex                                    \right] / \sqrt g \nonumber \\
!latex                       & - & \sum_u \left[
!latex           \partial_\s (\sqrt g B^u) \; g_{u,\z} + (\sqrt g B^u) \; \partial_\s g_{u,\z} - (\sqrt g B^u) g_{u,\z} \; \partial_\s \sqrt g / \sqrt g
!latex                                    \right] / \sqrt g, \\
!latex           \sqrt g j^\z & = & \sum_u \left[
!latex           \partial_\s (\sqrt g B^u) \; g_{u,\t} + (\sqrt g B^u) \; \partial_\s g_{u,\t} - (\sqrt g B^u) g_{u,\t} \; \partial_\s \sqrt g / \sqrt g
!latex                                    \right] / \sqrt g \nonumber \\
!latex                       & - & \sum_u \left[
!latex           \partial_\t (\sqrt g B^u) \; g_{u,\s} + (\sqrt g B^u) \; \partial_\t g_{u,\s} - (\sqrt g B^u) g_{u,\s} \; \partial_\t \sqrt g / \sqrt g
!latex                                    \right] / \sqrt g.
!latex       \ee

   do ii = 1, 3

    select case( ii )
    case( 1 ) ; jj = 2 ; kk = 3
    case( 2 ) ; jj = 3 ; kk = 1
    case( 3 ) ; jj = 1 ; kk = 2
    end select

    gJu(1:Ntz,ii) = zero

    do uu = 1, 3 ! summation over uu;

     gJu(1:Ntz,ii) = gJu(1:Ntz,ii) &
                   + ( gBu(1:Ntz,uu,jj) * guvij(1:Ntz,uu,kk, 0) &
                     + gBu(1:Ntz,uu, 0) * guvij(1:Ntz,uu,kk,jj) &
                     - gBu(1:Ntz,uu, 0) * guvij(1:Ntz,uu,kk, 0) * sg(1:Ntz,jj) / sg(1:Ntz,0) &
                     ) &
                   - ( gBu(1:Ntz,uu,kk) * guvij(1:Ntz,uu,jj, 0) &
                     + gBu(1:Ntz,uu, 0) * guvij(1:Ntz,uu,jj,kk) &
                     - gBu(1:Ntz,uu, 0) * guvij(1:Ntz,uu,jj, 0) * sg(1:Ntz,kk) / sg(1:Ntz,0) &
                     )
    enddo ! end of do uu;

    gJu(1:Ntz,ii) = gJu(1:Ntz,ii) / sg(1:Ntz,0)

   enddo ! end of do ii;

   select case (Igeometry)
   case (1)
   case (2)
   case (3)
    B_cartesian(:,1) = (gBu(1:Ntz,1,0)*Rij(1:Ntz,1,0) + gBu(1:Ntz,2,0)*Rij(1:Ntz,2,0) + gBu(1:Ntz,3,0)*Rij(1:Ntz,3,0))/sg(1:Ntz,0)
    B_cartesian(:,2) = (gBu(1:Ntz,1,0)*Zij(1:Ntz,1,0) + gBu(1:Ntz,2,0)*Zij(1:Ntz,2,0) + gBu(1:Ntz,3,0)*Zij(1:Ntz,3,0))/sg(1:Ntz,0)
    B_cartesian(:,3) = (                                                                gBu(1:Ntz,3,0)*Rij(1:Ntz,0,0))/sg(1:Ntz,0)

    J_cartesian(:,1) = (gJu(1:Ntz,1)*Rij(1:Ntz,1,0) + gJu(1:Ntz,2)*Rij(1:Ntz,2,0) + gJu(1:Ntz,3)*Rij(1:Ntz,3,0))/sg(1:Ntz,0)
    J_cartesian(:,2) = (gJu(1:Ntz,1)*Zij(1:Ntz,1,0) + gJu(1:Ntz,2)*Zij(1:Ntz,2,0) + gJu(1:Ntz,3)*Zij(1:Ntz,3,0))/sg(1:Ntz,0)
    J_cartesian(:,3) = (                                                            gJu(1:Ntz,3)*Rij(1:Ntz,0,0))/sg(1:Ntz,0)
   end select

!latex \end{enumerate}

!latex \item The final calculation of the error, which is written to screen, is a sum over the angular grid:
!latex       \be E^\s & \equiv & \frac{1}{N} \sum_k \omega_k \sum_{i,j} | \sqrt g j^\s - \mu \sqrt g B^\s |, \\
!latex           E^\t & \equiv & \frac{1}{N} \sum_k \omega_k \sum_{i,j} | \sqrt g j^\t - \mu \sqrt g B^\t |, \\
!latex           E^\z & \equiv & \frac{1}{N} \sum_k \omega_k \sum_{i,j} | \sqrt g j^\z - \mu \sqrt g B^\z |,
!latex       \ee
!latex       where $N\equiv \sum_{i,j}1$.
   if (Lerrortype.eq.1 .and. Igeometry .eq. 3) then
    do ii = 1, 3 ; jerror(ii) = jerror(ii) + weight(jquad) *sum(sg(1:Ntz,0)*abs(  J_cartesian(1:Ntz,ii) - mu(lvol) * B_cartesian(1:Ntz,ii)  ) )
                  !if (maxval(abs(  J_cartesian(1:Ntz,ii) - mu(lvol) * B_cartesian(1:Ntz,ii)  )) > jerrormax(ii)) write(ounit,*) ii, lss,maxval(abs(  J_cartesian(1:Ntz,ii) - mu(lvol) * B_cartesian(1:Ntz,ii)  ))
    ;            ; jerrormax(ii) = max(jerrormax(ii), maxval(abs(  J_cartesian(1:Ntz,ii) - mu(lvol) * B_cartesian(1:Ntz,ii)  )))
    enddo
   else
    do ii = 1, 3 ; jerror(ii) = jerror(ii) + weight(jquad) * sum( abs(  gJu(1:Ntz,ii) - mu(lvol) * gBu(1:Ntz,ii,0)  ) )
    ;            ; jerrormax(ii) = max(jerrormax(ii), maxval(abs(  gJu(1:Ntz,ii) - mu(lvol) * gBu(1:Ntz,ii,0)  )/ sg(1:Ntz,0)))
    enddo
   endif
   intvol = intvol + weight(jquad) * sum(sg(1:Ntz,0))
  enddo ! end of do jquad;

  beltramierror(lvol,1:3) = jerror(1:3) / Ntz  ! the 'tranditional' SPEC error
  jerror(1:3) = jerror(1:3) / intvol
  beltramierror(lvol,4:6) = jerror(1:3)        ! the volume average error
  beltramierror(lvol,7:9) = jerrormax(1:3)     ! the max error

  if (Lerrortype.eq.1 .and. Igeometry .eq. 3) then
    cput = GETTIME ; write(ounit,1002) cput-cpus, myid, lvol, Lrad(lvol), jerror(1:3), cput-cpui ! write error to screen;
    ;              ; write(ounit,1003) cput-cpus, myid, lvol, Lrad(lvol), jerrormax(1:3), cput-cpui ! write error to screen;
  else
    cput = GETTIME ; write(ounit,1004) cput-cpus, myid, lvol, Lrad(lvol), jerror(1:3), cput-cpui ! write error to screen;
    ;              ; write(ounit,1005) cput-cpus, myid, lvol, Lrad(lvol), jerrormax(1:3), cput-cpui ! write error to screen;
  endif

!latex \item The error is stored into an array called \type{beltramierror} which is then written to the HDF5 file in \link{hdfint}.


1002 format("jo00aa : ",f10.2," : myid=",i3," ; lvol =",i3," ; lrad =",i3," ; AVG E^\R="es23.15" , E^\Z="es23.15" , E^\phi="es23.15" ; time="f8.2"s ;")
1003 format("jo00aa : ",f10.2," : myid=",i3," ; lvol =",i3," ; lrad =",i3," ; MAX E^\R="es23.15" , E^\Z="es23.15" , E^\phi="es23.15" ; time="f8.2"s ;")
1004 format("jo00aa : ",f10.2," : myid=",i3," ; lvol =",i3," ; lrad =",i3," ; AVG E^\s="es23.15" , E^\t="es23.15" , E^\z="es23.15" ; time="f8.2"s ;")
1005 format("jo00aa : ",f10.2," : myid=",i3," ; lvol =",i3," ; lrad =",i3," ; MAX E^\s="es23.15" , E^\t="es23.15" , E^\z="es23.15" ; time="f8.2"s ;")

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! check boundary condition
  jerror = zero
  ivol = lvol

  do kk = 0, Nz-1 ; zeta = kk * pi2nfp / Nz
    do jj = 0, Nt-1 ; teta = jj * pi2    / Nt ; jk = 1 + jj + kk*Nt ; st(1:2) = (/ one, teta /)

    WCALL( jo00aa, bfield, ( zeta, st(1:Node), Bst(1:Node) ) )

    jerror(2) = max(jerror(2), abs(Bst(1) * gBzeta))

    enddo
  enddo

  if (.not.Lcoordinatesingularity) then
    do kk = 0, Nz-1 ; zeta = kk * pi2nfp / Nz
      do jj = 0, Nt-1 ; teta = jj * pi2    / Nt ; jk = 1 + jj + kk*Nt ; st(1:2) = (/ -one, teta /)

      WCALL( jo00aa, bfield, ( zeta, st(1:Node), Bst(1:Node) ) )

      jerror(1) = max(jerror(1), abs(Bst(1) * gBzeta))

      enddo
    enddo
  endif
  cput = GETTIME ; write(ounit,1006) cput-cpus, myid, lvol, Lrad(lvol), jerror(1:2), cput-cpui ! write error to screen;

  ! check fluxes
  Bst = zero

  if (Lcoordinatesingularity) then
    do ll = 0, Lrad(lvol)
        Bst(1) = Bst(1) + Ate(lvol,0,1)%s(ll) * RTT(ll,0,1,0)
      enddo
    Bst(1) = abs(Bst(1) - dtflux(lvol))
  else
    do ll = 0, Lrad(lvol)
        Bst(1) = Bst(1) + Ate(lvol,0,1)%s(ll) * TT(ll,1,0)
    enddo
    Bst(1) = abs(Bst(1) - dtflux(lvol))
    do ll = 0, Lrad(lvol)
        Bst(2) = Bst(2) - Aze(lvol,0,1)%s(ll) * TT(ll,1,0)
    enddo
    Bst(2) = abs(Bst(2) - dpflux(lvol))
  endif

  cput = GETTIME ; write(ounit,1007) cput-cpus, myid, lvol, Lrad(lvol), Bst(1:2), cput-cpui ! write error to screen;

1006 format("jo00aa : ",f10.2," : myid=",i3," ; lvol =",i3," ; lrad =",i3," ; MAX gB^s(-1)="es23.15" , gB^s(+1) ="es23.15" ; time="f8.2"s ;")
1007 format("jo00aa : ",f10.2," : myid=",i3," ; lvol =",i3," ; lrad =",i3," ; dtfluxERR   ="es23.15" , dpfluxERR="es23.15" ; time="f8.2"s ;")

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(jo00aa)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine jo00aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \end{enumerate}

!latex \subsection{comments}

!latex \begin{enumerate}
!latex \item Is there a better definition and quantification of the error? For example, should we employ an error measure that is dimensionless?
!latex \item If the coordinate singularity is in the domain, then $|\nabla\t|\rightarrow\infty$ at the coordinate origin.
!latex       What then happens to $||\left( {\bf j}-\mu {\bf B}\right)\cdot\nabla \t ||$ as defined in \Eqn{Et}?
!latex \item What is the predicted scaling of the error in the Chebyshev-Fourier representation scale with numerical resolution?
!latex       Note that the predicted error scaling for $E^\s$, $E^\t$ and $E^\z$ may not be standard,
!latex       as various radial derivatives are taken to compute the components of ${\bf j}$.
!latex       (See for example the discussion in Sec.IV.C in
!latex       [\doi{10.1063/1.4765691}{Hudson, Dewar {\em et al.}, Phys. Plasmas {\bf 19}, 112502 (2012)}],
!latex       where the expected scaling of the error for a finite-element implementation is confirmed numerically.)
!latex \item Instead of using Gaussian integration to compute the integral over $s$, an adaptive quadrature algorithm may be preferable.
!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
