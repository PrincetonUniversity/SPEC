!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (diagnostic) ! Measures error in Beltrami equation, $\nabla \times {\bf B} - \mu {\bf B}$.

!latex \briefly{Measures error in magnetic field, $||\nabla\times{\bf B}-\mu {\bf B}||$.}

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
!latex \item Define the following measures of the error:
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

  use constants, only : zero, half, one, two

  use fileunits, only : ounit, lunit

  use inputlist, only : Wmacros, Wjo00aa, ext, Nvol, Lrad, mu

  use cputiming, only : Tjo00aa

  use allglobal, only : myid, cpus, ivol, &
                        im, in, regumm, &
                        Mvol, &
                        cheby, &
                        Ate, Aze, Ato, Azo, &
                        sg, guvij, Rij, Zij, &
                        Nt, Nz, efmn, ofmn, cfmn, sfmn, &
                        NOTstellsym, &
                        Lcoordinatesingularity, &
                        IBerror

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS
  
  INTEGER, intent(in) :: lvol, Ntz, lquad, mn ! these are really global, but are included in argument list to remove allocations;
  
  INTEGER             :: jquad, Lcurvature, ll, ii, jj, kk, uu, ideriv
  
  REAL                :: lss, gBu(1:Ntz,1:3,0:3), gJu(1:Ntz,1:3), jerror(1:3)
  
  REAL                :: Atemn(1:mn,0:2), Azemn(1:mn,0:2), Atomn(1:mn,0:2), Azomn(1:mn,0:2), sbar, sbarhim(0:2)
  
  INTEGER             :: itype, id01bcf
  REAL                :: aa, bb, cc, dd, weight(1:lquad), abscis(1:lquad), workfield(1:2*lquad)
  
  BEGIN(jo00aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  FATAL( jo00aa, lquad.lt.1, invalid lquad supplied to jo00aa )
  FATAL( jo00aa, lvol.lt.1 .or. lvol.gt.Mvol, invalid interface label )
#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \item The integration over $s$i s performed using Gaussian integration, e.g. $\ds \int \!\! f(s) ds \approx \sum_k \omega_k f(s_k)$;
!latex       with the abscissae, $s_k$, and the weights, $\omega_k$, for $k=1,$ \internal{Iquad$_v$},
!latex       determined by the NAG routine \nag{www.nag.co.uk/numeric/FL/manual19/pdf/D01/d01bcf_fl19.pdf}{D01BCF}
!latex       (which is called with the input parameters \type{itype = 0}, \type{A = -1.}, \type{B = +1.}, \type{C = 0.} and \type{D = 0.}).
!latex       The resolution, \type{N} $ \equiv $ \internal{Iquad$_v$}, is determined by \inputvar{Nquad} 
!latex       (see \link{global} and \link{preset}).
!latex       A fatal error is enforced by \link{jo00aa} 
!latex       if \nag{www.nag.co.uk/numeric/FL/manual19/pdf/D01/d01bcf_fl19.pdf}{D01BCF} returns an \type{ifail} $\ne 0$.

  itype = 0 ; aa = -one ; bb = +one ; cc = zero ; dd = zero ! prepare Gaussian quadrature;  6 Feb 13;
  
!  id01bcf = 1
!  call D01BCF( itype, aa, bb, cc, dd, lquad, weight(1:lquad), abscis(1:lquad), id01bcf ) ! sets gaussian weights & abscissae;
!  do ii = 1, lquad
!     print *, ii, weight(ii),abscis(ii)
!  end do
!      SUBROUTINE CDGQF(NT,T,WTS,KIND,ALPHA,BETA,NWF,WF,IER)
  call CDGQF(lquad, abscis(1:lquad), weight(1:lquad), itype+1, aa, bb, 2*lquad, workfield, id01bcf)
!  do ii = 1, lquad
!     print *, ii, weight(ii),abscis(ii)
!  end do
  
  cput= GETTIME
  select case( id01bcf ) !                                                         123456789012345
  case( 0 )    ;  if( Wjo00aa ) write(ounit,1000) cput-cpus, myid, lvol, id01bcf, "success        ", abscis(1:lquad)
  case( 1 )    ;                write(ounit,1000) cput-cpus, myid, lvol, id01bcf, "failed         ", abscis(1:lquad)
  case( 2 )    ;                write(ounit,1000) cput-cpus, myid, lvol, id01bcf, "input error    ", abscis(1:lquad)
  case( 3 )    ;                write(ounit,1000) cput-cpus, myid, lvol, id01bcf, "input error    ", abscis(1:lquad)
  case( 4 )    ;                write(ounit,1000) cput-cpus, myid, lvol, id01bcf, "weight overflow", abscis(1:lquad)
  case( 5 )    ;                write(ounit,1000) cput-cpus, myid, lvol, id01bcf, "weight zero    ", abscis(1:lquad)
  case( 6 )    ;                write(ounit,1000) cput-cpus, myid, lvol, id01bcf, "failed         ", abscis(1:lquad)
  case default ;                write(ounit,1000) cput-cpus, myid, lvol, id01bcf, "weird          ", abscis(1:lquad)
  end select
  
  if( Wjo00aa )                 write(ounit,1001)                                                    weight(1:lquad)
  
1000 format("jo00aa : ",f10.2," : myid=",i3," ; lvol=",i3," ; id01bcf=",i3," ; "a15" ;":" abscissae ="99f10.6)
1001 format("jo00aa : ", 10x ," :      "3x"        "3x"           "3x"   "15x" ;":" weights   ="99f10.6)
  
  FATAL( jo00aa, id01bcf.ne.0, failed to construct Gaussian integration abscisae and weights )
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  jerror(1:3) = zero ; ideriv = 0 ! three components of the error in \curl B - mu B; initialize summation; 6 Feb 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \item Inside the Gaussian quadrature loop, i.e. for each $s_k$, 

!latex \begin{enumerate}

  do jquad = 1, lquad ! loop over radial sub-sub-grid (numerical quadrature);
   
   lss = abscis(jquad) ; sbar = ( lss + one ) * half
   
!latex \item The metric elements, $g_{\mu,\nu} \equiv $ \type{gij(1:6,0,1:Ntz)}, and the Jacobian, $\sqrt g \equiv $ \type{sg(0,1:Ntz)}, 
!latex       are calculated on a regular angular grid, $(\t_i,\z_j)$, in \link{coords}.
!latex       The derivatives $\partial_i g_{\mu,\nu} \equiv$ \type{gij(1:6,i,1:Ntz)} and $\partial_i \sqrt g \equiv $ \type{sg(i,1:Ntz)}, 
!latex       with respect to $ i \in \{ \s,\t,\z \}$ are also returned.

   Lcurvature = 2
   
   WCALL( jo00aa, coords, ( lvol, lss, Lcurvature, Ntz, mn ) ) ! returns coordinates, metrics, . . .
   
   cheby(0,0:2) = (/ one, zero, zero /) ! Chebyshev initialization; 16 Jan 13;
   cheby(1,0:2) = (/ lss,  one, zero /) ! Chebyshev initialization; 16 Jan 13;
   do ll = 2, Lrad(lvol) ; cheby(ll,0:2) = (/ two * lss * cheby(ll-1,0)                                                         - cheby(ll-2,0) , &
                                              two       * cheby(ll-1,0) + two * lss * cheby(ll-1,1)                             - cheby(ll-2,1) , &
                                              two       * cheby(ll-1,1) + two       * cheby(ll-1,1) + two * lss * cheby(ll-1,2) - cheby(ll-2,2) /)
   enddo ! end of do ll; 20 Jun 14;
    
   Atemn(1:mn,0:2) = zero ! initialize summation over Chebyshev polynomials;  6 Feb 13;
   Azemn(1:mn,0:2) = zero
   Atomn(1:mn,0:2) = zero
   Azomn(1:mn,0:2) = zero
   
   do ll = 0, Lrad(lvol)
     
    do ii = 1, mn 
      
!latex \item The Fourier components of the vector potential given in \Eqn{At} and \Eqn{Az}, and their first and second radial derivatives, are summed.

     if( Lcoordinatesingularity ) then
      sbarhim(0) = sbar**regumm(ii) ; sbarhim(1) = half * regumm(ii) * sbarhim(0) / sbar ; sbarhim(2) = half * ( regumm(ii)-one ) * sbarhim(1) / sbar
     else                              
      sbarhim(0) = one              ; sbarhim(1) = zero                                  ; sbarhim(2) = zero
     endif
      
     ;Atemn(ii,0) = Atemn(ii,0) + Ate(lvol,ideriv,ii)%s(ll) * ( cheby(ll,0) * sbarhim(0)                                                             )
     ;Atemn(ii,1) = Atemn(ii,1) + Ate(lvol,ideriv,ii)%s(ll) * ( cheby(ll,1) * sbarhim(0)                                  + cheby(ll,0) * sbarhim(1) )
     ;Atemn(ii,2) = Atemn(ii,2) + Ate(lvol,ideriv,ii)%s(ll) * ( cheby(ll,2) * sbarhim(0) + two * cheby(ll,1) * sbarhim(1) + cheby(ll,0) * sbarhim(2) )

     ;Azemn(ii,0) = Azemn(ii,0) + Aze(lvol,ideriv,ii)%s(ll) * ( cheby(ll,0) * sbarhim(0)                                                             )
     ;Azemn(ii,1) = Azemn(ii,1) + Aze(lvol,ideriv,ii)%s(ll) * ( cheby(ll,1) * sbarhim(0)                                  + cheby(ll,0) * sbarhim(1) )
     ;Azemn(ii,2) = Azemn(ii,2) + Aze(lvol,ideriv,ii)%s(ll) * ( cheby(ll,2) * sbarhim(0) + two * cheby(ll,1) * sbarhim(1) + cheby(ll,0) * sbarhim(2) )

     if( NOTstellsym ) then

      Atomn(ii,0) = Atomn(ii,0) + Ato(lvol,ideriv,ii)%s(ll) * ( cheby(ll,0) * sbarhim(0)                                                             )
      Atomn(ii,1) = Atomn(ii,1) + Ato(lvol,ideriv,ii)%s(ll) * ( cheby(ll,1) * sbarhim(0)                                  + cheby(ll,0) * sbarhim(1) )
      Atomn(ii,2) = Atomn(ii,2) + Ato(lvol,ideriv,ii)%s(ll) * ( cheby(ll,2) * sbarhim(0) + two * cheby(ll,1) * sbarhim(1) + cheby(ll,0) * sbarhim(2) )

      Azomn(ii,0) = Azomn(ii,0) + Azo(lvol,ideriv,ii)%s(ll) * ( cheby(ll,0) * sbarhim(0)                                                             )
      Azomn(ii,1) = Azomn(ii,1) + Azo(lvol,ideriv,ii)%s(ll) * ( cheby(ll,1) * sbarhim(0)                                  + cheby(ll,0) * sbarhim(1) )
      Azomn(ii,2) = Azomn(ii,2) + Azo(lvol,ideriv,ii)%s(ll) * ( cheby(ll,2) * sbarhim(0) + two * cheby(ll,1) * sbarhim(1) + cheby(ll,0) * sbarhim(2) )

     endif ! end of if( NOTstellsym) ; 20 Jun 14;
      
    enddo ! end of do ii;  6 Feb 13;
    
   enddo ! end of do ll;  6 Feb 13;
    
!latex \item The quantities $\sqrt g B^\s$, $\sqrt g B^\t$ and $\sqrt g B^\z$, and their first and second derivatives with respect to $(\s,\t,\z)$, 
!latex       are computed on the regular angular grid (using FFTs).

   ofmn(1:mn) = -          im(1:mn)*Azemn(1:mn,0) -          in(1:mn)*Atemn(1:mn,0)
   efmn(1:mn) = +          im(1:mn)*Azomn(1:mn,0) +          in(1:mn)*Atomn(1:mn,0)
   sfmn(1:mn) = -          im(1:mn)*Azemn(1:mn,1) -          in(1:mn)*Atemn(1:mn,1)
   cfmn(1:mn) = +          im(1:mn)*Azomn(1:mn,1) +          in(1:mn)*Atomn(1:mn,1)
   
   call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), &
                Nt, Nz, gBu(1:Ntz,1,0), gBu(1:Ntz,1,1) ) !  (gB^s)   , d(gB^s)/ds; 11 Mar 16;
   
   efmn(1:mn) = - im(1:mn)*im(1:mn)*Azemn(1:mn,0) - im(1:mn)*in(1:mn)*Atemn(1:mn,0)
   ofmn(1:mn) = - im(1:mn)*im(1:mn)*Azomn(1:mn,0) - im(1:mn)*in(1:mn)*Atomn(1:mn,0)
   cfmn(1:mn) = + in(1:mn)*im(1:mn)*Azemn(1:mn,0) + in(1:mn)*in(1:mn)*Atemn(1:mn,0)
   sfmn(1:mn) = + in(1:mn)*im(1:mn)*Azomn(1:mn,0) + in(1:mn)*in(1:mn)*Atomn(1:mn,0)
   
   call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), &
                Nt, Nz, gBu(1:Ntz,1,2), gBu(1:Ntz,1,3) ) ! d(gB^s)/dt, d(gB^s)/dz; 11 Mar 16;
   
   efmn(1:mn) =                                   -                   Azemn(1:mn,1)
   ofmn(1:mn) =                                   -                   Azomn(1:mn,1)
   cfmn(1:mn) =                                   -                   Azemn(1:mn,2)
   sfmn(1:mn) =                                   -                   Azomn(1:mn,2)
   
   call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), &
                Nt, Nz, gBu(1:Ntz,2,0), gBu(1:Ntz,2,1) ) !  (gB^t)   , d(gB^t)/ds; 11 Mar 16;

   ofmn(1:mn) =                                   +          im(1:mn)*Azemn(1:mn,1)
   efmn(1:mn) =                                   -          im(1:mn)*Azomn(1:mn,1)
   sfmn(1:mn) =                                   -          in(1:mn)*Azemn(1:mn,1)
   cfmn(1:mn) =                                   +          in(1:mn)*Azomn(1:mn,1)
   
   call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), &
                Nt, Nz, gBu(1:Ntz,2,2), gBu(1:Ntz,2,3) ) ! d(gB^t)/dt, d(gB^t)/dz; 11 Mar 16;
   
   efmn(1:mn) = +                   Atemn(1:mn,1)
   ofmn(1:mn) = +                   Atomn(1:mn,1)
   cfmn(1:mn) = +                   Atemn(1:mn,2)
   sfmn(1:mn) = +                   Atomn(1:mn,2)
   
   call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), &
                Nt, Nz, gBu(1:Ntz,3,0), gBu(1:Ntz,3,1) ) !  (gB^z)   , d(gB^z)/ds; 11 Mar 16;

   ofmn(1:mn) = -          im(1:mn)*Atemn(1:mn,1)
   efmn(1:mn) = +          im(1:mn)*Atomn(1:mn,1)
   sfmn(1:mn) = +          in(1:mn)*Atemn(1:mn,1)
   cfmn(1:mn) = -          in(1:mn)*Atomn(1:mn,1)
   
   call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), &
                Nt, Nz, gBu(1:Ntz,3,2), gBu(1:Ntz,3,3) ) ! d(gB^z)/dt, d(gB^z)/dz; 11 Mar 16;
     
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
    
    do uu = 1, 3 ! summation over uu;  6 Feb 13;
     
     gJu(1:Ntz,ii) = gJu(1:Ntz,ii) &
                   + ( gBu(1:Ntz,uu,jj) * guvij(1:Ntz,uu,kk, 0) &
                     + gBu(1:Ntz,uu, 0) * guvij(1:Ntz,uu,kk,jj) &
                     - gBu(1:Ntz,uu, 0) * guvij(1:Ntz,uu,kk, 0) * sg(1:Ntz,jj) / sg(1:Ntz,0) &
                     ) &
                   - ( gBu(1:Ntz,uu,kk) * guvij(1:Ntz,uu,jj, 0) &
                     + gBu(1:Ntz,uu, 0) * guvij(1:Ntz,uu,jj,kk) &
                     - gBu(1:Ntz,uu, 0) * guvij(1:Ntz,uu,jj, 0) * sg(1:Ntz,kk) / sg(1:Ntz,0) &
                     )
    enddo ! end of do uu;  6 Feb 13;
     
    gJu(1:Ntz,ii) = gJu(1:Ntz,ii) / sg(1:Ntz,0)
    
   enddo ! end of do ii;  5 Feb 13;

!latex \end{enumerate}
    
!latex \item The final calculation of the error, which is written to screen, is a sum over the angular grid:
!latex       \be E^\s & \equiv & \frac{1}{N} \sum_k \omega_k \sum_{i,j} | \sqrt g j^\s - \mu \sqrt g B^\s |, \\
!latex           E^\t & \equiv & \frac{1}{N} \sum_k \omega_k \sum_{i,j} | \sqrt g j^\t - \mu \sqrt g B^\t |, \\
!latex           E^\z & \equiv & \frac{1}{N} \sum_k \omega_k \sum_{i,j} | \sqrt g j^\z - \mu \sqrt g B^\z |,
!latex       \ee
!latex       where $N\equiv \sum_{i,j}$.

   do ii = 1, 3 ; jerror(ii) = jerror(ii) + weight(jquad) * sum( abs(  gJu(1:Ntz,ii) - mu(lvol) * gBu(1:Ntz,ii,0)  ) )
   enddo
   
  enddo ! end of do jquad;  5 Feb 13;
   

  jerror(1:3) = jerror(1:3) / Ntz
   
  cput = GETTIME ; write(ounit,1002) cput-cpus, myid, lvol, Lrad(lvol), jerror(1:3), cput-cpui ! write error to screen;
   
1002 format("jo00aa : ",f10.2," : myid=",i3," ; lvol =",i3," ; lrad =",i3," ; E^\s="es12.5", E^\t="es12.5", E^\z="es12.5" ; time="f8.2"s ;")
  
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
!latex       [\paper{Hudson, Dewar {\em et al.}}{S.R. Hudson, R.L. Dewar {\em et al}.}{10.1063/1.4765691}{Phys. Plasmas}{19}{112502}{2012}],
!latex       where the expected scaling of the error for a finite-element implementation is confirmed numerically.)
!latex \item Instead of using Gaussian integration to compute the integral over $s$, an adaptive quadrature algorithm may be preferable.
!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
