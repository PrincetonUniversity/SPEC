!> \defgroup grp_integrals Integrals

!> \file ma00aa.f90
!> \brief Calculates volume integrals of Chebyshev polynomials and metric element products.
!> \ingroup grp_integrals

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> **Chebyshev-metric information**
!> <ul>
!> <li> The following quantities are calculated:
!>
!>       \f{eqnarray}{ \verb+DToocc(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}^\prime \; {\overline T}_{p,j}        \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \cos\alpha_i \cos\alpha_j                  \\
!>                     \verb+DToocs(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}^\prime \; {\overline T}_{p,j}        \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \cos\alpha_i \sin\alpha_j                  \\
!>                     \verb+DToosc(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}^\prime \; {\overline T}_{p,j}        \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \sin\alpha_i \cos\alpha_j                  \\
!>                     \verb+DTooss(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}^\prime \; {\overline T}_{p,j}        \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \sin\alpha_i \sin\alpha_j                  
!>       \f}
!>
!>       \f{eqnarray}{ \verb+TTsscc(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}        \; {\overline T}_{p,j}        \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \cos\alpha_i \cos\alpha_j \; \bar g_{ss} \\
!>                     \verb+TTsscs(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}        \; {\overline T}_{p,j}        \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \cos\alpha_i \sin\alpha_j \; \bar g_{ss} \\
!>                     \verb+TTsssc(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}        \; {\overline T}_{p,j}        \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \sin\alpha_i \cos\alpha_j \; \bar g_{ss} \\
!>                     \verb+TTssss(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}        \; {\overline T}_{p,j}        \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \sin\alpha_i \sin\alpha_j \; \bar g_{ss}
!>       \f}
!>
!>       \f{eqnarray}{ \verb+TDstcc(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}        \; {\overline T}_{p,j}^\prime \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \cos\alpha_i \cos\alpha_j \; \bar g_{s\theta} \\
!>                     \verb+TDstcs(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}        \; {\overline T}_{p,j}^\prime \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \cos\alpha_i \sin\alpha_j \; \bar g_{s\theta} \\
!>                     \verb+TDstsc(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}        \; {\overline T}_{p,j}^\prime \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \sin\alpha_i \cos\alpha_j \; \bar g_{s\theta} \\
!>                     \verb+TDstss(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}        \; {\overline T}_{p,j}^\prime \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \sin\alpha_i \sin\alpha_j \; \bar g_{s\theta}
!>       \f}
!>
!>       \f{eqnarray}{ \verb+TDstcc(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}        \; {\overline T}_{p,j}^\prime \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \cos\alpha_i \cos\alpha_j \; \bar g_{s\zeta} \\
!>                     \verb+TDstcs(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}        \; {\overline T}_{p,j}^\prime \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \cos\alpha_i \sin\alpha_j \; \bar g_{s\zeta} \\
!>                     \verb+TDstsc(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}        \; {\overline T}_{p,j}^\prime \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \sin\alpha_i \cos\alpha_j \; \bar g_{s\zeta} \\
!>                     \verb+TDstss(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}        \; {\overline T}_{p,j}^\prime \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \sin\alpha_i \sin\alpha_j \; \bar g_{s\zeta}
!>       \f}
!>
!>       \f{eqnarray}{ \verb+DDstcc(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}^\prime \; {\overline T}_{p,j}^\prime \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \cos\alpha_i \cos\alpha_j \; \bar g_{\theta\theta} \\
!>                     \verb+DDstcs(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}^\prime \; {\overline T}_{p,j}^\prime \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \cos\alpha_i \sin\alpha_j \; \bar g_{\theta\theta} \\
!>                     \verb+DDstsc(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}^\prime \; {\overline T}_{p,j}^\prime \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \sin\alpha_i \cos\alpha_j \; \bar g_{\theta\theta} \\
!>                     \verb+DDstss(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}^\prime \; {\overline T}_{p,j}^\prime \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \sin\alpha_i \sin\alpha_j \; \bar g_{\theta\theta}
!>       \f}
!>
!>       \f{eqnarray}{ \verb+DDstcc(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}^\prime \; {\overline T}_{p,j}^\prime \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \cos\alpha_i \cos\alpha_j \; \bar g_{\theta\zeta} \\
!>                     \verb+DDstcs(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}^\prime \; {\overline T}_{p,j}^\prime \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \cos\alpha_i \sin\alpha_j \; \bar g_{\theta\zeta} \\
!>                     \verb+DDstsc(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}^\prime \; {\overline T}_{p,j}^\prime \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \sin\alpha_i \cos\alpha_j \; \bar g_{\theta\zeta} \\
!>                     \verb+DDstss(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}^\prime \; {\overline T}_{p,j}^\prime \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \sin\alpha_i \sin\alpha_j \; \bar g_{\theta\zeta}
!>       \f}
!>
!>       \f{eqnarray}{ \verb+DDstcc(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}^\prime \; {\overline T}_{p,j}^\prime \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \cos\alpha_i \cos\alpha_j \; \bar g_{\zeta\zeta} \\
!>                     \verb+DDstcs(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}^\prime \; {\overline T}_{p,j}^\prime \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \cos\alpha_i \sin\alpha_j \; \bar g_{\zeta\zeta} \\
!>                     \verb+DDstsc(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}^\prime \; {\overline T}_{p,j}^\prime \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \sin\alpha_i \cos\alpha_j \; \bar g_{\zeta\zeta} \\
!>                     \verb+DDstss(l,p,i,j)+ & \equiv & \int ds \; {\overline T}_{l,i}^\prime \; {\overline T}_{p,j}^\prime \; \oint\!\!\!\oint \!d\theta d\zeta \,\,\, \sin\alpha_i \sin\alpha_j \; \bar g_{\zeta\zeta}
!>       \f}
!>
!>       where \f${\overline T}_{l,i}\equiv T_l \, \bar s^{m_i/2}\f$ if the domain includes the coordinate singularity, and \f${\overline T}_{l,i}\equiv T_l\f$ if not;
!>       and \f$\bar g_{\mu\nu} \equiv g_{\mu\nu} / \sqrt g\f$. </li>
!>
!> <li> The double-angle formulae are used to reduce the above expressions to the Fourier harmonics of \f$\bar g_{\mu\nu}\f$:
!>       see \c kija and \c kijs, which are defined in preset.f90 . </li>
!>
!> </ul>

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> \brief Calculates volume integrals of Chebyshev polynomials and metric element products.
!> 
!> @param[in] lquad degree of quadrature
!> @param[in] mn    number of Fourier harmonics
!> @param[in] lvol  index of nested volume
!> @param[in] lrad  order of Chebychev polynomials
subroutine ma00aa( lquad, mn, lvol, lrad )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half, one, two, pi, pi2
  
  use numerical, only : vsmall, small, sqrtmachprec
  
  use fileunits, only : ounit
  
  use inputlist, only : Wma00aa
  
  use cputiming, only : Tma00aa
  
  use allglobal, only : myid, ncpu, cpus, &
                        Mvol, im, in, mne, &
                        gaussianweight, gaussianabscissae, &
                        DToocc, DToocs, DToosc, DTooss, &
                        TTsscc, TTsscs, TTsssc, TTssss, &
                        TDstcc, TDstcs, TDstsc, TDstss, &
                        TDszcc, TDszcs, TDszsc, TDszss, &
                        DDttcc, DDttcs, DDttsc, DDttss, &
                        DDtzcc, DDtzcs, DDtzsc, DDtzss, &
                        DDzzcc, DDzzcs, DDzzsc, DDzzss, &
                        ki, kija, kijs, &
                        goomne, goomno, &
                        gssmne, gssmno, &
                        gstmne, gstmno, &
                        gszmne, gszmno, &
                        gttmne, gttmno, &
                        gtzmne, gtzmno, &
                        gzzmne, gzzmno, &
                        cheby, &
                        Lcoordinatesingularity, regumm, &
                        pi2pi2nfp, pi2pi2nfphalf
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in) :: lquad, mn, lvol, lrad
  
  INTEGER             :: jquad, ll, pp, uv, ii, jj, io, mn2, lp2, mn2_max, lp2_max
  
  INTEGER             :: kk, kd, kka, kks, kda, kds
  
  REAL                :: lss, jthweight, fee, feo, foe, foo, Tl, Dl, Tp, Dp, TlTp, TlDp, DlTp, DlDp, ikda, ikds, imn2, ilrad

  REAL                :: foocc, foocs, foosc, fooss
  REAL                :: fsscc, fsscs, fsssc, fssss
  REAL                :: fstcc, fstcs, fstsc, fstss
  REAL                :: fszcc, fszcs, fszsc, fszss
  REAL                :: fttcc, fttcs, fttsc, fttss
  REAL                :: ftzcc, ftzcs, ftzsc, ftzss
  REAL                :: fzzcc, fzzcs, fzzsc, fzzss
  
  REAL                :: sbar(1:lquad), halfoversbar(1:lquad), sbarhim(1:lquad,1:mn) ! regularization factors; 10 Dec 15;
  
  BEGIN( ma00aa )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  FATAL( ma00aa, lvol.lt.1 .or. lvol.gt.Mvol, illegal volume label )
#endif
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  mn2_max = mn*mn
  lp2_max = (lrad+1)*(lrad+1)
  imn2    =  one/real(mn)
  ilrad = one/real(lrad+1)

  DToocc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
  DToocs( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
  DToosc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
  DTooss( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero

  TTsscc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
  TTsscs( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
  TTsssc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
  TTssss( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero

  TDstcc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
  TDstcs( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
  TDstsc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
  TDstss( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero

  TDszcc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
  TDszcs( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
  TDszsc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
  TDszss( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero

  DDttcc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
  DDttcs( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
  DDttsc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
  DDttss( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero

  DDtzcc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
  DDtzcs( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
  DDtzsc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
  DDtzss( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero

  DDzzcc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
  DDzzcs( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
  DDzzsc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
  DDzzss( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Lcoordinatesingularity ) then
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   sbar(1:lquad) = ( gaussianabscissae(1:lquad,lvol) + one ) * half
   
   halfoversbar(1:lquad) = half / sbar(1:lquad)
   
   do jquad = 1, lquad ; sbarhim(jquad,1:mn) = sbar(jquad)**regumm(1:mn) ! pre-calculation of regularization factor; 12 Sep 13;
   enddo
   
   do jquad = 1, lquad ! Gaussian quadrature loop;
    
    lss = gaussianabscissae(jquad,lvol) ; jthweight = gaussianweight(jquad,lvol)
    
    ;                 cheby( 0,0:1) = (/ one                                       , zero                                                            /)
    ;                 cheby( 1,0:1) = (/ lss                                       , one                                                             /)
    do ll = 2, lrad ; cheby(ll,0:1) = (/ two * lss * cheby(ll-1,0) - cheby(ll-2,0) , two * cheby(ll-1,0) + two * lss * cheby(ll-1,1) - cheby(ll-2,1) /)
    enddo
    
    WCALL( ma00aa, metrix,( lvol, lss ) ) ! compute metric elements; 16 Jan 13;

    do mn2 = 1, mn2_max
      ii = mod(mn2-1,mn)+1
      jj = floor(real(mn2-1) * imn2)+1
      
      kks = kijs(ii,jj,0) !; kds = kijs(ii,jj,1) 
      kka = kija(ii,jj,0) !; kda = kija(ii,jj,1) 
      ikds = jthweight * sbarhim(jquad,ii)* sbarhim(jquad,jj) / kijs(ii,jj,1)
      ikda = jthweight * sbarhim(jquad,ii)* sbarhim(jquad,jj) / kija(ii,jj,1)

      foocc = ( + goomne(kks) * abs(ikds) + goomne(kka) * abs(ikda) )
      foocs = ( - goomno(kks) *     ikds  + goomno(kka) *     ikda  )
      foosc = ( + goomno(kks) *     ikds  + goomno(kka) *     ikda  )
      fooss = ( + goomne(kks) * abs(ikds) - goomne(kka) * abs(ikda) )
      
      fsscc = ( + gssmne(kks) * abs(ikds) + gssmne(kka) * abs(ikda) )
      fsscs = ( - gssmno(kks) *     ikds  + gssmno(kka) *     ikda  )
      fsssc = ( + gssmno(kks) *     ikds  + gssmno(kka) *     ikda  )
      fssss = ( + gssmne(kks) * abs(ikds) - gssmne(kka) * abs(ikda) )
      
      fstcc = ( + gstmne(kks) * abs(ikds) + gstmne(kka) * abs(ikda) )
      fstcs = ( - gstmno(kks) *     ikds  + gstmno(kka) *     ikda  )
      fstsc = ( + gstmno(kks) *     ikds  + gstmno(kka) *     ikda  )
      fstss = ( + gstmne(kks) * abs(ikds) - gstmne(kka) * abs(ikda) )
      
      fszcc = ( + gszmne(kks) * abs(ikds) + gszmne(kka) * abs(ikda) )
      fszcs = ( - gszmno(kks) *     ikds  + gszmno(kka) *     ikda  )
      fszsc = ( + gszmno(kks) *     ikds  + gszmno(kka) *     ikda  )
      fszss = ( + gszmne(kks) * abs(ikds) - gszmne(kka) * abs(ikda) )
      
      fttcc = ( + gttmne(kks) * abs(ikds) + gttmne(kka) * abs(ikda) )
      fttcs = ( - gttmno(kks) *     ikds  + gttmno(kka) *     ikda  )
      fttsc = ( + gttmno(kks) *     ikds  + gttmno(kka) *     ikda  )
      fttss = ( + gttmne(kks) * abs(ikds) - gttmne(kka) * abs(ikda) )
      
      ftzcc = ( + gtzmne(kks) * abs(ikds) + gtzmne(kka) * abs(ikda) )
      ftzcs = ( - gtzmno(kks) *     ikds  + gtzmno(kka) *     ikda  )
      ftzsc = ( + gtzmno(kks) *     ikds  + gtzmno(kka) *     ikda  )
      ftzss = ( + gtzmne(kks) * abs(ikds) - gtzmne(kka) * abs(ikda) )
      
      fzzcc = ( + gzzmne(kks) * abs(ikds) + gzzmne(kka) * abs(ikda) )
      fzzcs = ( - gzzmno(kks) *     ikds  + gzzmno(kka) *     ikda  )
      fzzsc = ( + gzzmno(kks) *     ikds  + gzzmno(kka) *     ikda  )
      fzzss = ( + gzzmne(kks) * abs(ikds) - gzzmne(kka) * abs(ikda) )

      do lp2 = 1, lp2_max 
        ll = mod(lp2-1,lrad+1)
        pp = floor(real(lp2-1) * ilrad) 
       
        Tl =                                      cheby(ll,0)                 ! this is the only difference for Lcoordinatesingularity;
        Dl = ( regumm(ii) * halfoversbar(jquad) * cheby(ll,0) + cheby(ll,1) ) ! this is the only difference for Lcoordinatesingularity;
          
        Tp =                                       cheby(pp,0)                 ! this is the only difference for Lcoordinatesingularity;
        Dp =  ( regumm(jj) * halfoversbar(jquad) * cheby(pp,0) + cheby(pp,1) ) ! this is the only difference for Lcoordinatesingularity;
        
        TlTp = Tl * Tp
        TlDp = Tl * Dp
        DlTp = Dl * Tp
        DlDp = Dl * Dp 

        DToocc( ll, pp, ii, jj ) = DToocc( ll, pp, ii, jj ) + DlTp * foocc
        DToocs( ll, pp, ii, jj ) = DToocs( ll, pp, ii, jj ) + DlTp * foocs
        DToosc( ll, pp, ii, jj ) = DToosc( ll, pp, ii, jj ) + DlTp * foosc
        DTooss( ll, pp, ii, jj ) = DTooss( ll, pp, ii, jj ) + DlTp * fooss
        
        TTsscc( ll, pp, ii, jj ) = TTsscc( ll, pp, ii, jj ) + TlTp * fsscc
        TTsscs( ll, pp, ii, jj ) = TTsscs( ll, pp, ii, jj ) + TlTp * fsscs
        TTsssc( ll, pp, ii, jj ) = TTsssc( ll, pp, ii, jj ) + TlTp * fsssc
        TTssss( ll, pp, ii, jj ) = TTssss( ll, pp, ii, jj ) + TlTp * fssss
        
        TDstcc( ll, pp, ii, jj ) = TDstcc( ll, pp, ii, jj ) + TlDp * fstcc
        TDstcs( ll, pp, ii, jj ) = TDstcs( ll, pp, ii, jj ) + TlDp * fstcs
        TDstsc( ll, pp, ii, jj ) = TDstsc( ll, pp, ii, jj ) + TlDp * fstsc
        TDstss( ll, pp, ii, jj ) = TDstss( ll, pp, ii, jj ) + TlDp * fstss
        
        TDszcc( ll, pp, ii, jj ) = TDszcc( ll, pp, ii, jj ) + TlDp * fszcc
        TDszcs( ll, pp, ii, jj ) = TDszcs( ll, pp, ii, jj ) + TlDp * fszcs
        TDszsc( ll, pp, ii, jj ) = TDszsc( ll, pp, ii, jj ) + TlDp * fszsc
        TDszss( ll, pp, ii, jj ) = TDszss( ll, pp, ii, jj ) + TlDp * fszss
        
        DDttcc( ll, pp, ii, jj ) = DDttcc( ll, pp, ii, jj ) + DlDp * fttcc
        DDttcs( ll, pp, ii, jj ) = DDttcs( ll, pp, ii, jj ) + DlDp * fttcs
        DDttsc( ll, pp, ii, jj ) = DDttsc( ll, pp, ii, jj ) + DlDp * fttsc
        DDttss( ll, pp, ii, jj ) = DDttss( ll, pp, ii, jj ) + DlDp * fttss
        
        DDtzcc( ll, pp, ii, jj ) = DDtzcc( ll, pp, ii, jj ) + DlDp * ftzcc
        DDtzcs( ll, pp, ii, jj ) = DDtzcs( ll, pp, ii, jj ) + DlDp * ftzcs
        DDtzsc( ll, pp, ii, jj ) = DDtzsc( ll, pp, ii, jj ) + DlDp * ftzsc
        DDtzss( ll, pp, ii, jj ) = DDtzss( ll, pp, ii, jj ) + DlDp * ftzss
        
        DDzzcc( ll, pp, ii, jj ) = DDzzcc( ll, pp, ii, jj ) + DlDp * fzzcc
        DDzzcs( ll, pp, ii, jj ) = DDzzcs( ll, pp, ii, jj ) + DlDp * fzzcs
        DDzzsc( ll, pp, ii, jj ) = DDzzsc( ll, pp, ii, jj ) + DlDp * fzzsc
        DDzzss( ll, pp, ii, jj ) = DDzzss( ll, pp, ii, jj ) + DlDp * fzzss
       
      enddo ! end of do lp2; 08 Feb 16;
     
    enddo ! end of do mn2; 08 Feb 16;
    
   enddo ! end of do jquad; ! 16 Jan 13;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
  else ! .not.Lcoordinatesingularity; 17 Dec 15;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   do jquad = 1, lquad ! Gaussian quadrature loop;
    
    lss = gaussianabscissae(jquad,lvol) ; jthweight = gaussianweight(jquad,lvol)
    
    ;                 cheby( 0,0:1) = (/ one                                       , zero                                                            /)
    ;                 cheby( 1,0:1) = (/ lss                                       , one                                                             /)
    do ll = 2, lrad ; cheby(ll,0:1) = (/ two * lss * cheby(ll-1,0) - cheby(ll-2,0) , two * cheby(ll-1,0) + two * lss * cheby(ll-1,1) - cheby(ll-2,1) /)
    enddo
    
    WCALL( ma00aa, metrix,( lvol, lss ) ) ! compute metric elements; 16 Jan 13;

    do mn2 = 1, mn2_max

      ii = mod(mn2-1,mn)+1

      jj = floor(real(mn2-1) * imn2)+1
      
      kks = kijs(ii,jj,0) ! ; kds = kijs(ii,jj,1) 
      kka = kija(ii,jj,0) !; kda = kija(ii,jj,1) 
      ikds = jthweight / kijs(ii,jj,1)
      ikda = jthweight / kija(ii,jj,1)
      
      foocc = ( + goomne(kks) * abs(ikds) + goomne(kka) * abs(ikda) )
      foocs = ( - goomno(kks) *     ikds  + goomno(kka) *     ikda  )
      foosc = ( + goomno(kks) *     ikds  + goomno(kka) *     ikda  )
      fooss = ( + goomne(kks) * abs(ikds) - goomne(kka) * abs(ikda) )
      
      fsscc = ( + gssmne(kks) * abs(ikds) + gssmne(kka) * abs(ikda) )
      fsscs = ( - gssmno(kks) *     ikds  + gssmno(kka) *     ikda  )
      fsssc = ( + gssmno(kks) *     ikds  + gssmno(kka) *     ikda  )
      fssss = ( + gssmne(kks) * abs(ikds) - gssmne(kka) * abs(ikda) )
      
      fstcc = ( + gstmne(kks) * abs(ikds) + gstmne(kka) * abs(ikda) )
      fstcs = ( - gstmno(kks) *     ikds  + gstmno(kka) *     ikda  )
      fstsc = ( + gstmno(kks) *     ikds  + gstmno(kka) *     ikda  )
      fstss = ( + gstmne(kks) * abs(ikds) - gstmne(kka) * abs(ikda) )
      
      fszcc = ( + gszmne(kks) * abs(ikds) + gszmne(kka) * abs(ikda) )
      fszcs = ( - gszmno(kks) *     ikds  + gszmno(kka) *     ikda  )
      fszsc = ( + gszmno(kks) *     ikds  + gszmno(kka) *     ikda  )
      fszss = ( + gszmne(kks) * abs(ikds) - gszmne(kka) * abs(ikda) )
      
      fttcc = ( + gttmne(kks) * abs(ikds) + gttmne(kka) * abs(ikda) )
      fttcs = ( - gttmno(kks) *     ikds  + gttmno(kka) *     ikda  )
      fttsc = ( + gttmno(kks) *     ikds  + gttmno(kka) *     ikda  )
      fttss = ( + gttmne(kks) * abs(ikds) - gttmne(kka) * abs(ikda) )
      
      ftzcc = ( + gtzmne(kks) * abs(ikds) + gtzmne(kka) * abs(ikda) )
      ftzcs = ( - gtzmno(kks) *     ikds  + gtzmno(kka) *     ikda  )
      ftzsc = ( + gtzmno(kks) *     ikds  + gtzmno(kka) *     ikda  )
      ftzss = ( + gtzmne(kks) * abs(ikds) - gtzmne(kka) * abs(ikda) )
      
      fzzcc = ( + gzzmne(kks) * abs(ikds) + gzzmne(kka) * abs(ikda) )
      fzzcs = ( - gzzmno(kks) *     ikds  + gzzmno(kka) *     ikda  )
      fzzsc = ( + gzzmno(kks) *     ikds  + gzzmno(kka) *     ikda  )
      fzzss = ( + gzzmne(kks) * abs(ikds) - gzzmne(kka) * abs(ikda) )
      
      do lp2 = 1, lp2_max 
        ll = mod(lp2-1,lrad+1)
        pp = floor(real(lp2-1) * ilrad) 
        
        Tl = cheby(ll,0)
        Dl = cheby(ll,1)
        Tp = cheby(pp,0)
        Dp = cheby(pp,1)
        
        TlTp = Tl * Tp
        TlDp = Tl * Dp
        DlTp = Dl * Tp
        DlDp = Dl * Dp

        DToocc( ll, pp, ii, jj ) = DToocc( ll, pp, ii, jj ) + DlTp * foocc
        DToocs( ll, pp, ii, jj ) = DToocs( ll, pp, ii, jj ) + DlTp * foocs
        DToosc( ll, pp, ii, jj ) = DToosc( ll, pp, ii, jj ) + DlTp * foosc
        DTooss( ll, pp, ii, jj ) = DTooss( ll, pp, ii, jj ) + DlTp * fooss
        
        TTsscc( ll, pp, ii, jj ) = TTsscc( ll, pp, ii, jj ) + TlTp * fsscc
        TTsscs( ll, pp, ii, jj ) = TTsscs( ll, pp, ii, jj ) + TlTp * fsscs
        TTsssc( ll, pp, ii, jj ) = TTsssc( ll, pp, ii, jj ) + TlTp * fsssc
        TTssss( ll, pp, ii, jj ) = TTssss( ll, pp, ii, jj ) + TlTp * fssss
        
        TDstcc( ll, pp, ii, jj ) = TDstcc( ll, pp, ii, jj ) + TlDp * fstcc
        TDstcs( ll, pp, ii, jj ) = TDstcs( ll, pp, ii, jj ) + TlDp * fstcs
        TDstsc( ll, pp, ii, jj ) = TDstsc( ll, pp, ii, jj ) + TlDp * fstsc
        TDstss( ll, pp, ii, jj ) = TDstss( ll, pp, ii, jj ) + TlDp * fstss
        
        TDszcc( ll, pp, ii, jj ) = TDszcc( ll, pp, ii, jj ) + TlDp * fszcc
        TDszcs( ll, pp, ii, jj ) = TDszcs( ll, pp, ii, jj ) + TlDp * fszcs
        TDszsc( ll, pp, ii, jj ) = TDszsc( ll, pp, ii, jj ) + TlDp * fszsc
        TDszss( ll, pp, ii, jj ) = TDszss( ll, pp, ii, jj ) + TlDp * fszss
        
        DDttcc( ll, pp, ii, jj ) = DDttcc( ll, pp, ii, jj ) + DlDp * fttcc
        DDttcs( ll, pp, ii, jj ) = DDttcs( ll, pp, ii, jj ) + DlDp * fttcs
        DDttsc( ll, pp, ii, jj ) = DDttsc( ll, pp, ii, jj ) + DlDp * fttsc
        DDttss( ll, pp, ii, jj ) = DDttss( ll, pp, ii, jj ) + DlDp * fttss
        
        DDtzcc( ll, pp, ii, jj ) = DDtzcc( ll, pp, ii, jj ) + DlDp * ftzcc
        DDtzcs( ll, pp, ii, jj ) = DDtzcs( ll, pp, ii, jj ) + DlDp * ftzcs
        DDtzsc( ll, pp, ii, jj ) = DDtzsc( ll, pp, ii, jj ) + DlDp * ftzsc
        DDtzss( ll, pp, ii, jj ) = DDtzss( ll, pp, ii, jj ) + DlDp * ftzss
        
        DDzzcc( ll, pp, ii, jj ) = DDzzcc( ll, pp, ii, jj ) + DlDp * fzzcc
        DDzzcs( ll, pp, ii, jj ) = DDzzcs( ll, pp, ii, jj ) + DlDp * fzzcs
        DDzzsc( ll, pp, ii, jj ) = DDzzsc( ll, pp, ii, jj ) + DlDp * fzzsc
        DDzzss( ll, pp, ii, jj ) = DDzzss( ll, pp, ii, jj ) + DlDp * fzzss
       
      enddo ! end of do lp2;  1 Feb 13;
      
    enddo ! end of do mn2;  1 Feb 13;
    
   enddo ! end of do jquad; ! 16 Jan 13;
    
  endif ! end of if( Lcoordinatesingularity ) ; 17 Dec 15;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  DToocc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = DToocc( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
  DToocs( 0:lrad, 0:lrad, 1:mn, 1:mn ) = DToocs( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
  DToosc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = DToosc( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
  DTooss( 0:lrad, 0:lrad, 1:mn, 1:mn ) = DTooss( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf

  TTsscc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = TTsscc( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
  TTsscs( 0:lrad, 0:lrad, 1:mn, 1:mn ) = TTsscs( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
  TTsssc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = TTsssc( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
  TTssss( 0:lrad, 0:lrad, 1:mn, 1:mn ) = TTssss( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf

  TDstcc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = TDstcc( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
  TDstcs( 0:lrad, 0:lrad, 1:mn, 1:mn ) = TDstcs( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
  TDstsc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = TDstsc( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
  TDstss( 0:lrad, 0:lrad, 1:mn, 1:mn ) = TDstss( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf

  TDszcc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = TDszcc( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
  TDszcs( 0:lrad, 0:lrad, 1:mn, 1:mn ) = TDszcs( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
  TDszsc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = TDszsc( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
  TDszss( 0:lrad, 0:lrad, 1:mn, 1:mn ) = TDszss( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf

  DDttcc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = DDttcc( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
  DDttcs( 0:lrad, 0:lrad, 1:mn, 1:mn ) = DDttcs( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
  DDttsc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = DDttsc( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
  DDttss( 0:lrad, 0:lrad, 1:mn, 1:mn ) = DDttss( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf

  DDtzcc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = DDtzcc( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
  DDtzcs( 0:lrad, 0:lrad, 1:mn, 1:mn ) = DDtzcs( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
  DDtzsc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = DDtzsc( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
  DDtzss( 0:lrad, 0:lrad, 1:mn, 1:mn ) = DDtzss( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf

  DDzzcc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = DDzzcc( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
  DDzzcs( 0:lrad, 0:lrad, 1:mn, 1:mn ) = DDzzcs( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
  DDzzsc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = DDzzsc( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
  DDzzss( 0:lrad, 0:lrad, 1:mn, 1:mn ) = DDzzss( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  
  if( Wma00aa ) then
   do ll = 0, lrad
    do pp = 0, lrad
     do ii = 1, mn
      do jj = 1, mn
       write(ounit,1000) ll, pp, ii, jj, (/ DToocc( ll, pp, ii, jj ), DToocs( ll, pp, ii, jj ), DToosc( ll, pp, ii, jj ), DTooss( ll, pp, ii, jj ) /) / pi/pi
      enddo
     enddo
    enddo
   enddo
  endif
  
1000 format("ma00aa : ll="i3" ; pp="i3" ; ii="i3" ; jj="i3" ; DToocc="f15.10" ; DToocs="f15.10" ; DToosc="f15.10" ; DTooss="f15.10" ;")
  
#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN( ma00aa )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine ma00aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
