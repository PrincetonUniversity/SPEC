!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Calculates volume integrals of Chebyshev polynomials and metric element products.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \end{enumerate} \subsection{Chebyshev-metric information} \begin{enumerate}

!latex \item When the {\em vector} potential is used to construct the field the following quantities are calculated:

!latex       \be \verb+TTee(0,l,p,i,j)+ & \equiv & \int ds \; T_l^\prime \; T_p        \; \ooint               \cos\a_i \cos\a_j \\
!latex           \verb+TTeo(0,l,p,i,j)+ & \equiv & \int ds \; T_l^\prime \; T_p        \; \ooint               \cos\a_i \sin\a_j \\
!latex           \verb+TToe(0,l,p,i,j)+ & \equiv & \int ds \; T_l^\prime \; T_p        \; \ooint               \sin\a_i \cos\a_j \\
!latex           \verb+TToo(0,l,p,i,j)+ & \equiv & \int ds \; T_l^\prime \; T_p        \; \ooint               \sin\a_i \sin\a_j
!latex       \ee

!latex       \be \verb+TTee(1,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p        \; \ooint \bar g_{\s\s} \cos\a_i \cos\a_j \\
!latex           \verb+TTeo(1,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p        \; \ooint \bar g_{\s\s} \cos\a_i \sin\a_j \\
!latex           \verb+TToe(1,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p        \; \ooint \bar g_{\s\s} \sin\a_i \cos\a_j \\
!latex           \verb+TToo(1,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p        \; \ooint \bar g_{\s\s} \sin\a_i \sin\a_j
!latex       \ee

!latex       \be \verb+TTee(2,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p^\prime \; \ooint \bar g_{\s\t} \cos\a_i \cos\a_j \\
!latex           \verb+TTeo(2,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p^\prime \; \ooint \bar g_{\s\t} \cos\a_i \sin\a_j \\
!latex           \verb+TToe(2,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p^\prime \; \ooint \bar g_{\s\t} \sin\a_i \cos\a_j \\
!latex           \verb+TToo(2,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p^\prime \; \ooint \bar g_{\s\t} \sin\a_i \sin\a_j
!latex       \ee

!latex       \be \verb+TTee(3,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p^\prime \; \ooint \bar g_{\s\z} \cos\a_i \cos\a_j \\
!latex           \verb+TTeo(3,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p^\prime \; \ooint \bar g_{\s\z} \cos\a_i \sin\a_j \\
!latex           \verb+TToe(3,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p^\prime \; \ooint \bar g_{\s\z} \sin\a_i \cos\a_j \\
!latex           \verb+TToo(3,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p^\prime \; \ooint \bar g_{\s\z} \sin\a_i \sin\a_j
!latex       \ee

!latex       \be \verb+TTee(4,l,p,i,j)+ & \equiv & \int ds \; T_l^\prime \; T_p^\prime \; \ooint \bar g_{\t\t} \cos\a_i \cos\a_j \\
!latex           \verb+TTeo(4,l,p,i,j)+ & \equiv & \int ds \; T_l^\prime \; T_p^\prime \; \ooint \bar g_{\t\t} \cos\a_i \sin\a_j \\
!latex           \verb+TToe(4,l,p,i,j)+ & \equiv & \int ds \; T_l^\prime \; T_p^\prime \; \ooint \bar g_{\t\t} \sin\a_i \cos\a_j \\
!latex           \verb+TToo(4,l,p,i,j)+ & \equiv & \int ds \; T_l^\prime \; T_p^\prime \; \ooint \bar g_{\t\t} \sin\a_i \sin\a_j
!latex       \ee

!latex       \be \verb+TTee(5,l,p,i,j)+ & \equiv & \int ds \; T_l^\prime \; T_p^\prime \; \ooint \bar g_{\t\z} \cos\a_i \cos\a_j \\
!latex           \verb+TTeo(5,l,p,i,j)+ & \equiv & \int ds \; T_l^\prime \; T_p^\prime \; \ooint \bar g_{\t\z} \cos\a_i \sin\a_j \\
!latex           \verb+TToe(5,l,p,i,j)+ & \equiv & \int ds \; T_l^\prime \; T_p^\prime \; \ooint \bar g_{\t\z} \sin\a_i \cos\a_j \\
!latex           \verb+TToo(5,l,p,i,j)+ & \equiv & \int ds \; T_l^\prime \; T_p^\prime \; \ooint \bar g_{\t\z} \sin\a_i \sin\a_j
!latex       \ee

!latex       \be \verb+TTee(6,l,p,i,j)+ & \equiv & \int ds \; T_l^\prime \; T_p^\prime \; \ooint \bar g_{\z\z} \cos\a_i \cos\a_j \\
!latex           \verb+TTeo(6,l,p,i,j)+ & \equiv & \int ds \; T_l^\prime \; T_p^\prime \; \ooint \bar g_{\z\z} \cos\a_i \sin\a_j \\
!latex           \verb+TToe(6,l,p,i,j)+ & \equiv & \int ds \; T_l^\prime \; T_p^\prime \; \ooint \bar g_{\z\z} \sin\a_i \cos\a_j \\
!latex           \verb+TToo(6,l,p,i,j)+ & \equiv & \int ds \; T_l^\prime \; T_p^\prime \; \ooint \bar g_{\z\z} \sin\a_i \sin\a_j
!latex       \ee

!latex       where $\bar g_{uv} \equiv g_{uv} / \sqrt g$.

!latex \item When the {\em scalar} potential is used to construct the field the following quantities are calculated:

!!atex       \be \verb+TTee(0,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p        \; \ooint               \cos\a_i \cos\a_j \\
!!atex           \verb+TTeo(0,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p        \; \ooint               \cos\a_i \sin\a_j \\
!!atex           \verb+TToe(0,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p        \; \ooint               \sin\a_i \cos\a_j \\
!!atex           \verb+TToo(0,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p        \; \ooint               \sin\a_i \sin\a_j
!!atex       \ee

!latex       \be \verb+TTee(1,l,p,i,j)+ & \equiv & \int ds \; T_l^\prime \; T_p^\prime \; \ooint \bar g^{\s\s} \cos\a_i \cos\a_j \\
!latex           \verb+TTeo(1,l,p,i,j)+ & \equiv & \int ds \; T_l^\prime \; T_p^\prime \; \ooint \bar g^{\s\s} \cos\a_i \sin\a_j \\
!latex           \verb+TToe(1,l,p,i,j)+ & \equiv & \int ds \; T_l^\prime \; T_p^\prime \; \ooint \bar g^{\s\s} \sin\a_i \cos\a_j \\
!latex           \verb+TToo(1,l,p,i,j)+ & \equiv & \int ds \; T_l^\prime \; T_p^\prime \; \ooint \bar g^{\s\s} \sin\a_i \sin\a_j
!latex       \ee

!latex       \be \verb+TTee(2,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p^\prime \; \ooint \bar g^{\s\t} \cos\a_i \cos\a_j \\
!latex           \verb+TTeo(2,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p^\prime \; \ooint \bar g^{\s\t} \cos\a_i \sin\a_j \\
!latex           \verb+TToe(2,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p^\prime \; \ooint \bar g^{\s\t} \sin\a_i \cos\a_j \\
!latex           \verb+TToo(2,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p^\prime \; \ooint \bar g^{\s\t} \sin\a_i \sin\a_j
!latex       \ee

!latex       \be \verb+TTee(3,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p^\prime \; \ooint \bar g^{\s\z} \cos\a_i \cos\a_j \\
!latex           \verb+TTeo(3,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p^\prime \; \ooint \bar g^{\s\z} \cos\a_i \sin\a_j \\
!latex           \verb+TToe(3,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p^\prime \; \ooint \bar g^{\s\z} \sin\a_i \cos\a_j \\
!latex           \verb+TToo(3,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p^\prime \; \ooint \bar g^{\s\z} \sin\a_i \sin\a_j
!latex       \ee

!latex       \be \verb+TTee(4,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p        \; \ooint \bar g^{\t\t} \cos\a_i \cos\a_j \\
!latex           \verb+TTeo(4,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p        \; \ooint \bar g^{\t\t} \cos\a_i \sin\a_j \\
!latex           \verb+TToe(4,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p        \; \ooint \bar g^{\t\t} \sin\a_i \cos\a_j \\
!latex           \verb+TToo(4,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p        \; \ooint \bar g^{\t\t} \sin\a_i \sin\a_j
!latex       \ee

!latex       \be \verb+TTee(5,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p        \; \ooint \bar g^{\t\z} \cos\a_i \cos\a_j \\
!latex           \verb+TTeo(5,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p        \; \ooint \bar g^{\t\z} \cos\a_i \sin\a_j \\
!latex           \verb+TToe(5,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p        \; \ooint \bar g^{\t\z} \sin\a_i \cos\a_j \\
!latex           \verb+TToo(5,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p        \; \ooint \bar g^{\t\z} \sin\a_i \sin\a_j
!latex       \ee

!latex       \be \verb+TTee(6,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p        \; \ooint \bar g^{\z\z} \cos\a_i \cos\a_j \\
!latex           \verb+TTeo(6,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p        \; \ooint \bar g^{\z\z} \cos\a_i \sin\a_j \\
!latex           \verb+TToe(6,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p        \; \ooint \bar g^{\z\z} \sin\a_i \cos\a_j \\
!latex           \verb+TToo(6,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p        \; \ooint \bar g^{\z\z} \sin\a_i \sin\a_j
!latex       \ee

!latex       where $\bar g^{uv} \equiv g^{uv} \sqrt g$.

!latex \end{enumerate} \subsubsection{angle integrations} \begin{enumerate}

!latex \item The integrals over the angles can be computed using double-angle formulae:
!latex       \be
!latex       \ooint \bar g_{\mu\nu} \; \cos\a_i \; \cos\a_j & = & \frac{1}{2} \ooint \bar g_{\mu\nu} ( \cos\a_{k_{ij}^-} + \cos\a_{k_{ij}^+} ) \\
!latex       \ooint \bar g_{\mu\nu} \; \sin\a_i \; \cos\a_j & = & \frac{1}{2} \ooint \bar g_{\mu\nu} ( \sin\a_{k_{ij}^-} + \sin\a_{k_{ij}^+} ) \\
!latex       \ooint \bar g_{\mu\nu} \; \sin\a_i \; \sin\a_j & = & \frac{1}{2} \ooint \bar g_{\mu\nu} ( \cos\a_{k_{ij}^-} - \cos\a_{k_{ij}^+} )
!latex       \ee
!latex       where $(m_{k_{ij}^-},n_{k_{ij}^-})=(m_i-m_j,n_i-n_j)$ and $(m_{k_{ij}^+},n_{k_{ij}^+})=(m_i+m_j,n_i+n_j)$.

!latex \end{enumerate} \subsubsection{regularization factors} \begin{enumerate}

!latex \item The regularization factors are easy to include.
!latex       Introduce $\bar s \equiv (s+1)/2$, and then 
!latex       \be \bar T_{i,l}        & \equiv & \bar s^{m_i/2} T_l \\
!latex           \bar T_{i,l}^\prime & \equiv & \frac{m_i}{2} \bar s^{m_i/2} \frac{1}{\bar s} \frac{1}{2} T_l + \bar s^{m_i/2} T_l^\prime.
!latex       \ee

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine ma00aa( lquad, mn, lvol, lrad )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half, one, two, pi, pi2

  use numerical, only : vsmall, small, sqrtmachprec

  use fileunits, only : ounit, wunit

  use inputlist, only : Wma00aa, Igeometry

  use cputiming, only : Tma00aa

  use allglobal, only : myid, ncpu, cpus, &
                        Mvol, im, in, mne, &
                        gaussianweight, gaussianabscissae, &
                        TTee, TTeo, TToe, TToo, Te, To, &
                        guvmnks, guvmnka, guvmnk, &
                        guvmne, guvmno, &
                        lchebyshev, &
!                       gchebyshev, &
                        Lcoordinatesingularity, Lvacuumregion, Lplasmaregion, halfmm, & !regularmm, &
                        pi2pi2nfp, pi2pi2nfphalf, &
                        Llatex
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in) :: lquad, mn, lvol, lrad
  
! required for Gaussian integration routine;
  INTEGER             :: itype, id01bcf
  REAL                :: aa, bb, cc, dd

!loop control;
  INTEGER             :: jquad, ll, pp, uv, ii, jj, llp(0:6), ppp(0:6), kksk, kkak, kksd, kkad

  INTEGER             :: mi, ni, mj, nj ! shorthand; 17 Apr 13;

  REAL                :: lss, fee, feo, foe, foo, fe, fo ! 17 Apr 13;

  REAL                :: Tl, Tp, TlTp ! Chebyshev polynomials; 17 Apr 13;
  REAL                :: lcpu, sbar(1:lquad), halfoversbar(1:lquad), sbarhim(1:lquad,1:mn), factor(1:lquad,1:mn,0:lrad)
  REAl                :: jthweight ! Gaussian quadrature; 17 Apr 13;
  
  BEGIN(ma00aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  FATALMESS(ma00aa, lvol.lt.1 .or. lvol.gt.Mvol, illegal volume label)
#endif
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  TTee( 0:6, 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
  TTeo( 0:6, 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
  TToe( 0:6, 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
  TToo( 0:6, 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero

  if( Lvacuumregion ) then ; Te( 2:6, 0:lrad, 1:mn ) = zero
   ;                       ; To( 2:6, 0:lrad, 1:mn ) = zero
  endif
  
  if( Lplasmaregion ) then ; llp(0:6) = (/ 1, 0, 0, 0, 1, 1, 1 /) ! order of derivative of Tl;  2 Feb 13;  first argument; 20 Jan 15;
   ;                       ; ppp(0:6) = (/ 0, 0, 1, 1, 1, 1, 1 /) ! order of derivative of Tp;  2 Feb 13; second argument; 20 Jan 15;
  else                     ; llp(0:6) = (/ 0, 1, 0, 0, 0, 0, 0 /) ! order of derivative of Tl;  2 Feb 13;
   ;                       ; ppp(0:6) = (/ 0, 1, 1, 1, 0, 0, 0 /) ! order of derivative of Tp;  2 Feb 13;
  endif
  
  if( Lcoordinatesingularity ) then ! prepare miscellaneous factors;  2 Feb 13;
   
   sbar(1:lquad) = ( gaussianabscissae(1:lquad,lvol) + one ) * half

#ifdef DEBUG
   FATALMESS(ma00aa,minval(abs(sbar(1:lquad))).lt.small, divide by zero in ma00aa)
#endif

   halfoversbar(1:lquad) = half / sbar(1:lquad)
   
   do jquad = 1, lquad ; sbarhim(jquad,1:mn) = sbar(jquad)**halfmm(1:mn) ! regularization factor; 12 Sep 13;
   enddo

  endif ! end of if( Lcoordinatesingularity) ;  2 Feb 13;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! BEGIN Gaussian quadrature loop;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  do jquad = 1, lquad ! loop over radial sub-sub-grid (numerical quadrature);
   
   lss = gaussianabscissae(jquad,lvol) ; jthweight = gaussianweight(jquad,lvol) ! radial sub-sub-grid = Gassian integration; shorthand;  2 Feb 13;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   ;lchebyshev( 0,0:1) = (/ one, zero /) ! Chebyshev initialization; 16 Jan 13;
   ;lchebyshev( 1,0:1) = (/ lss,  one /) ! Chebyshev initialization; 16 Jan 13;   
   do ll = 2, lrad
    lchebyshev(ll,0:1) = (/ two * lss * lchebyshev(ll-1,0)                                  - lchebyshev(ll-2,0) , & ! Chebyshev recurrence;            ; 16 Jan 13;
                            two       * lchebyshev(ll-1,0) + two * lss * lchebyshev(ll-1,1) - lchebyshev(ll-2,1) /)  ! Chebyshev recurrence; derivatives; 16 Jan 13;
   enddo

!WHAT IS THE DIFFERENCE BETWEEN LCHEBYSHEV AS CALCULATED HERE AND GCHEBYSHEV WHICH IS CALCULATED IN AL00AA; 02 Apr 13;

   WCALL(ma00aa,me00ab,( lvol, lss )) ! compute metric elements; 16 Jan 13;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
   do uv = 0, 6 ! include "trivial" metric; 22 Jan 13;

    if( Lvacuumregion .and. uv.eq.0 ) cycle ! these quantities are only required for constructing the helicity; 17 Apr 13;

    do ii = 1, mn

     do jj = 1, mn

      kksk = guvmnks(ii,jj,0) ; kksd = guvmnks(ii,jj,1) 
      kkak = guvmnka(ii,jj,0) ; kkad = guvmnka(ii,jj,1) 
      
      fee = jthweight * (   guvmne(kksk,uv) / abs(kksd) + guvmne(kkak,uv) / abs(kkad) ) ! derived from \cos\a_i . \cos\a_j ; even-even; 17 Jan 13;
      feo = jthweight * ( - guvmno(kksk,uv) /     kksd  + guvmno(kkak,uv) /     kkad  ) ! derived from \sin\a_i . \cos\a_j ; odd -even; 17 Jan 13;
      foe = jthweight * (   guvmno(kksk,uv) /     kksd  + guvmno(kkak,uv) /     kkad  ) ! derived from \sin\a_i . \cos\a_j ; odd -even; 17 Jan 13;
      foo = jthweight * (   guvmne(kksk,uv) / abs(kksd) - guvmne(kkak,uv) / abs(kkad) ) ! derived from \sin\a_i . \sin\a_j ; odd -odd ; 17 Jan 13;

      do ll = 0, lrad
       
       if( Lcoordinatesingularity ) then ! must include additional radial factor; 15 Feb 13;
        select case( llp(uv) )
        case( 0 )    ; Tl = sbarhim(jquad,ii) * (                                    lchebyshev(ll,0)                    )
        case( 1 )    ; Tl = sbarhim(jquad,ii) * ( halfmm(ii) * halfoversbar(jquad) * lchebyshev(ll,0) + lchebyshev(ll,1) )
        case default ; FATALMESS(ma00aa, .true., illegal llp)
        end select
       else
        select case( llp(uv) ) ! 20 Jan 15;
        case( 0 )    ; Tl =                                                             lchebyshev(ll,0)                      ! 20 Jan 15;
        case( 1 )    ; Tl =                                                                                lchebyshev(ll,1)   ! 20 Jan 15;
        case default ; FATALMESS(ma00aa, .true., illegal llp) ! 20 Jan 15;
        end select
       !;            ; Tl = lchebyshev(ll,llp(uv)) ! 20 Jan 15;
       endif ! end of if( Lcoordinatesingularity );  1 Feb 13;
       
       do pp = 0, lrad

        if( Lcoordinatesingularity ) then ! then must include additional radial factor; 20 Jan 15;
         select case( ppp(uv) )
         case( 0 )    ; Tp = sbarhim(jquad,jj) * (                                    lchebyshev(pp,0)                    )
         case( 1 )    ; Tp = sbarhim(jquad,jj) * ( halfmm(jj) * halfoversbar(jquad) * lchebyshev(pp,0) + lchebyshev(pp,1) )
         case default ; FATALMESS(ma00aa, .true., illegal ppp)
         end select
        else
         select case( ppp(uv) ) ! 20 Jan 15;
         case( 0 )    ; Tp =                                                             lchebyshev(pp,0)                      ! 20 Jan 15;
         case( 1 )    ; Tp =                                                                                lchebyshev(pp,1)   ! 20 Jan 15;
         case default ; FATALMESS(ma00aa, .true., illegal ppp) ! 20 Jan 15;
         end select ! 20 Jan 15;
        !;            ; Tp = lchebyshev(pp,ppp(uv)) ! 20 Jan 15;
        endif ! end of if( Lcoordinatesingularity );  1 Feb 13;

        TlTp = Tl * Tp

        TTee( uv, ll, pp, ii, jj ) = TTee( uv, ll, pp, ii, jj ) + TlTp * fee
        TTeo( uv, ll, pp, ii, jj ) = TTeo( uv, ll, pp, ii, jj ) + TlTp * feo
        TToe( uv, ll, pp, ii, jj ) = TToe( uv, ll, pp, ii, jj ) + TlTp * foe
        TToo( uv, ll, pp, ii, jj ) = TToo( uv, ll, pp, ii, jj ) + TlTp * foo

       enddo ! end of do pp;  1 Feb 13;
      enddo ! end of do ll;  1 Feb 13;
      
     enddo ! end of do jj;  1 Feb 13;
    enddo ! end of do ii;  1 Feb 13;
    
   enddo ! end of do uv;  1 Feb 13;
   
   if( Lvacuumregion ) then ! ! still inside do jquad loop; 22 Apr 13;
    
   do uv = 2, 6
     
    do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
     
     kksk = guvmnk(ii,0) ; kksd = guvmnk(ii,1) 

     fe = jthweight * guvmne(kksk,uv) / kksd
     fo = jthweight * guvmno(kksk,uv) / kksd

     do ll = 0, lrad
      
      Tl = lchebyshev(ll,ppp(uv)) ! NEEDS CHECKING; 17 Apr 13;
      
      Te( uv, ll, ii ) = Te( uv, ll, ii ) + Tl * fe ! THERE ARE SYMMETRIES THAT CAN BE EXPLOITED ;  2 Feb 13;
      To( uv, ll, ii ) = To( uv, ll, ii ) + Tl * fo ! THERE ARE SYMMETRIES THAT CAN BE EXPLOITED ;  2 Feb 13;
      
     enddo ! end of do ll;  1 Feb 13;
     
    enddo ! end of do ii;  1 Feb 13;
    
   enddo ! end of do uv;  1 Feb 13;
   
  endif ! end of if( Lvacuumregion ) ; 17 Apr 13;
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  enddo ! end of do jquad; ! 16 Jan 13;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!END Gaussian quadrature loop;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  TTee( 0:6, 0:lrad, 0:lrad, 1:mn, 1:mn ) = TTee( 0:6, 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
  TTeo( 0:6, 0:lrad, 0:lrad, 1:mn, 1:mn ) = TTeo( 0:6, 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
  TToe( 0:6, 0:lrad, 0:lrad, 1:mn, 1:mn ) = TToe( 0:6, 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
  TToo( 0:6, 0:lrad, 0:lrad, 1:mn, 1:mn ) = TToo( 0:6, 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf

  if( Lvacuumregion ) then ; Te( 2:6, 0:lrad, 1:mn ) = Te( 2:6, 0:lrad, 1:mn ) * pi2pi2nfp
   ;                       ; To( 2:6, 0:lrad, 1:mn ) = To( 2:6, 0:lrad, 1:mn ) * pi2pi2nfp
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG

  if( Wma00aa ) then
   if( Lvacuumregion ) then
    do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
     if( ni.eq.0 ) cycle !for debugging non-axisymmetric harmonics; 03 Apr 13;
     do ll = 0, lrad
      write(ounit,'("ma00aa : ", 10x ," : mi=",i3," ; ni=",i3," ; Te="5f21.16" ; To="5f21.16" ;")') mi, ni, Te(2:6,ll,ii), To(2:6,ll,ii)
     enddo
    enddo
   endif
  endif

  if( Llatex ) then ! 20 Jan 15;
   
   write(wunit,'("\newpage \subsection{{volume integrated metric elements : volume="i3"; [ma00aa]}}")') lvol

   write(wunit,'("\begin{{enumerate}}")')
   
  !if( lvol.eq.1 ) then ! include regularization factors; 20 Jan 15;
    
    do ll = 0, lrad
     do pp = 0, lrad
      write(wunit,'("\item $l="i2"$ ; $p="i2"$ ;")') ll, pp ! 20 Jan 15;
      do ii = 1, mn
       do jj = 1, mn
        if( abs(TTee(0,ll,pp,ii,jj)).gt.sqrtmachprec ) then
        write(wunit,1000) ll, pp, ii, jj, ll, im(ii), pp, im(jj), im(ii), in(ii), im(jj), in(jj), TTee(0,ll,pp,ii,jj)
        endif
       enddo
      enddo
     enddo
    enddo
    
  !endif ! end of if( lvol.eq.1 ) ; 20 Jan 15;
  
   write(wunit,'("\end{{enumerate}}")') 

1000 format("\item [] \verb+TTee(0,"3(i1","),i1")+$\ds \equiv \int_{{-1}}^{{+1}} \!\!\!\!\! ds \;\; [T_{{"i2"}} \; \bar s^{{"i2"/2}}]^\prime \;\;",&
                                                                                                   "[T_{{"i2"}} \; \bar s^{{"i2"/2}}]            ",&
                                                          "\ooint \!\!\! \cos("i2"\;\t-"i2"\;\z) \cos("i2"\;\t-"i2"\;\z)="f15.3"$") ! 20 Jan 15;

  endif ! end of if( Llatex ) ; 20 Jan 15;
  
#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(ma00aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine ma00aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
