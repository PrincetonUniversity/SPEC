!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (integrals) ! Calculates volume integrals of Chebyshev-polynomials and metric elements.

!latex \briefly{Calculates volume integrals of Chebyshev polynomials and metric element products.}

!latex \calledby{\link{dforce}}
!latex \calls{\link{metrix}}

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
!latex           \verb+TDstsc(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}  \; \bT{p,j}' \; \ooint \sin\a_i \cos\a_j \; \bar g_{\s\t} \\
!latex           \verb+TDstss(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}  \; \bT{p,j}' \; \ooint \sin\a_i \sin\a_j \; \bar g_{\s\t}
!latex       \ee

!latex       \be \verb+TDstcc(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}  \; \bT{p,j}' \; \ooint \cos\a_i \cos\a_j \; \bar g_{\s\z} \\
!latex           \verb+TDstcs(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}  \; \bT{p,j}' \; \ooint \cos\a_i \sin\a_j \; \bar g_{\s\z} \\
!latex           \verb+TDstsc(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}  \; \bT{p,j}' \; \ooint \sin\a_i \cos\a_j \; \bar g_{\s\z} \\
!latex           \verb+TDstss(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}  \; \bT{p,j}' \; \ooint \sin\a_i \sin\a_j \; \bar g_{\s\z}
!latex       \ee

!latex       \be \verb+DDstcc(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}' \; \ooint \cos\a_i \cos\a_j \; \bar g_{\t\t} \\
!latex           \verb+DDstcs(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}' \; \ooint \cos\a_i \sin\a_j \; \bar g_{\t\t} \\
!latex           \verb+DDstsc(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}' \; \ooint \sin\a_i \cos\a_j \; \bar g_{\t\t} \\
!latex           \verb+DDstss(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}' \; \ooint \sin\a_i \sin\a_j \; \bar g_{\t\t}
!latex       \ee

!latex       \be \verb+DDstcc(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}' \; \ooint \cos\a_i \cos\a_j \; \bar g_{\t\z} \\
!latex           \verb+DDstcs(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}' \; \ooint \cos\a_i \sin\a_j \; \bar g_{\t\z} \\
!latex           \verb+DDstsc(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}' \; \ooint \sin\a_i \cos\a_j \; \bar g_{\t\z} \\
!latex           \verb+DDstss(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}' \; \ooint \sin\a_i \sin\a_j \; \bar g_{\t\z}
!latex       \ee

!latex       \be \verb+DDstcc(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}' \; \ooint \cos\a_i \cos\a_j \; \bar g_{\z\z} \\
!latex           \verb+DDstcs(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}' \; \ooint \cos\a_i \sin\a_j \; \bar g_{\z\z} \\
!latex           \verb+DDstsc(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}' \; \ooint \sin\a_i \cos\a_j \; \bar g_{\z\z} \\
!latex           \verb+DDstss(l,p,i,j)+ & \equiv & \int ds \; \bT{l,i}' \; \bT{p,j}' \; \ooint \sin\a_i \sin\a_j \; \bar g_{\z\z}
!latex       \ee

!latex       where $\bT{l,i}\equiv T_l \, \bar s^{m_i/2}$ if the domain includes the coordinate singularity, and $\bT{l,i}\equiv T_l$ if not;
!latex       and $\bar g_{\mu\nu} \equiv g_{\mu\nu} / \sqrt g$.

!latex \item The double-angle formulae are used to reduce the above expressions to the Fourier harmonics of $\bar g_{\mu\nu}$:
!latex       see \internal{kija} and \internal{kijs}, which are defined in \link{preset}.

!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

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
  
  INTEGER             :: jquad, ll, pp, uv, ii, jj, io
  
  INTEGER             :: kk, kd, kka, kks, kda, kds
  
  REAL                :: lss, jthweight, fee, feo, foe, foo, Tl, Dl, Tp, Dp, TlTp, TlDp, DlTp, DlDp

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
    
    do ii = 1, mn
     
     do jj = 1, mn
      
      kks = kijs(ii,jj,0) ; kds = kijs(ii,jj,1) 
      kka = kija(ii,jj,0) ; kda = kija(ii,jj,1) 
      
      foocc = jthweight * ( + goomne(kks) / abs(kds) + goomne(kka) / abs(kda) )
      foocs = jthweight * ( - goomno(kks) /     kds  + goomno(kka) /     kda  )
      foosc = jthweight * ( + goomno(kks) /     kds  + goomno(kka) /     kda  )
      fooss = jthweight * ( + goomne(kks) / abs(kds) - goomne(kka) / abs(kda) )
      
      fsscc = jthweight * ( + gssmne(kks) / abs(kds) + gssmne(kka) / abs(kda) )
      fsscs = jthweight * ( - gssmno(kks) /     kds  + gssmno(kka) /     kda  )
      fsssc = jthweight * ( + gssmno(kks) /     kds  + gssmno(kka) /     kda  )
      fssss = jthweight * ( + gssmne(kks) / abs(kds) - gssmne(kka) / abs(kda) )
      
      fstcc = jthweight * ( + gstmne(kks) / abs(kds) + gstmne(kka) / abs(kda) )
      fstcs = jthweight * ( - gstmno(kks) /     kds  + gstmno(kka) /     kda  )
      fstsc = jthweight * ( + gstmno(kks) /     kds  + gstmno(kka) /     kda  )
      fstss = jthweight * ( + gstmne(kks) / abs(kds) - gstmne(kka) / abs(kda) )
      
      fszcc = jthweight * ( + gszmne(kks) / abs(kds) + gszmne(kka) / abs(kda) )
      fszcs = jthweight * ( - gszmno(kks) /     kds  + gszmno(kka) /     kda  )
      fszsc = jthweight * ( + gszmno(kks) /     kds  + gszmno(kka) /     kda  )
      fszss = jthweight * ( + gszmne(kks) / abs(kds) - gszmne(kka) / abs(kda) )
      
      fttcc = jthweight * ( + gttmne(kks) / abs(kds) + gttmne(kka) / abs(kda) )
      fttcs = jthweight * ( - gttmno(kks) /     kds  + gttmno(kka) /     kda  )
      fttsc = jthweight * ( + gttmno(kks) /     kds  + gttmno(kka) /     kda  )
      fttss = jthweight * ( + gttmne(kks) / abs(kds) - gttmne(kka) / abs(kda) )
      
      ftzcc = jthweight * ( + gtzmne(kks) / abs(kds) + gtzmne(kka) / abs(kda) )
      ftzcs = jthweight * ( - gtzmno(kks) /     kds  + gtzmno(kka) /     kda  )
      ftzsc = jthweight * ( + gtzmno(kks) /     kds  + gtzmno(kka) /     kda  )
      ftzss = jthweight * ( + gtzmne(kks) / abs(kds) - gtzmne(kka) / abs(kda) )
      
      fzzcc = jthweight * ( + gzzmne(kks) / abs(kds) + gzzmne(kka) / abs(kda) )
      fzzcs = jthweight * ( - gzzmno(kks) /     kds  + gzzmno(kka) /     kda  )
      fzzsc = jthweight * ( + gzzmno(kks) /     kds  + gzzmno(kka) /     kda  )
      fzzss = jthweight * ( + gzzmne(kks) / abs(kds) - gzzmne(kka) / abs(kda) )
      
      do ll = 0, lrad
       
       Tl = sbarhim(jquad,ii) *                                      cheby(ll,0)                 ! this is the only difference for Lcoordinatesingularity;
       Dl = sbarhim(jquad,ii) * ( regumm(ii) * halfoversbar(jquad) * cheby(ll,0) + cheby(ll,1) ) ! this is the only difference for Lcoordinatesingularity;
       
       do pp = 0, lrad
          
        Tp = sbarhim(jquad,jj) *                                      cheby(pp,0)                 ! this is the only difference for Lcoordinatesingularity;
        Dp = sbarhim(jquad,jj) * ( regumm(jj) * halfoversbar(jquad) * cheby(pp,0) + cheby(pp,1) ) ! this is the only difference for Lcoordinatesingularity;
        
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

       enddo ! end of do pp ; 08 Feb 16;
       
      enddo ! end of do ll; 08 Feb 16;
      
     enddo ! end of do jj ; 08 Feb 16;
     
    enddo ! end of do ii; 08 Feb 16;
    
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
   
    do ii = 1, mn

     do jj = 1, mn
      
      kks = kijs(ii,jj,0) ; kds = kijs(ii,jj,1) 
      kka = kija(ii,jj,0) ; kda = kija(ii,jj,1) 
      
      foocc = jthweight * ( + goomne(kks) / abs(kds) + goomne(kka) / abs(kda) )
      foocs = jthweight * ( - goomno(kks) /     kds  + goomno(kka) /     kda  )
      foosc = jthweight * ( + goomno(kks) /     kds  + goomno(kka) /     kda  )
      fooss = jthweight * ( + goomne(kks) / abs(kds) - goomne(kka) / abs(kda) )
      
      fsscc = jthweight * ( + gssmne(kks) / abs(kds) + gssmne(kka) / abs(kda) )
      fsscs = jthweight * ( - gssmno(kks) /     kds  + gssmno(kka) /     kda  )
      fsssc = jthweight * ( + gssmno(kks) /     kds  + gssmno(kka) /     kda  )
      fssss = jthweight * ( + gssmne(kks) / abs(kds) - gssmne(kka) / abs(kda) )
      
      fstcc = jthweight * ( + gstmne(kks) / abs(kds) + gstmne(kka) / abs(kda) )
      fstcs = jthweight * ( - gstmno(kks) /     kds  + gstmno(kka) /     kda  )
      fstsc = jthweight * ( + gstmno(kks) /     kds  + gstmno(kka) /     kda  )
      fstss = jthweight * ( + gstmne(kks) / abs(kds) - gstmne(kka) / abs(kda) )
      
      fszcc = jthweight * ( + gszmne(kks) / abs(kds) + gszmne(kka) / abs(kda) )
      fszcs = jthweight * ( - gszmno(kks) /     kds  + gszmno(kka) /     kda  )
      fszsc = jthweight * ( + gszmno(kks) /     kds  + gszmno(kka) /     kda  )
      fszss = jthweight * ( + gszmne(kks) / abs(kds) - gszmne(kka) / abs(kda) )
      
      fttcc = jthweight * ( + gttmne(kks) / abs(kds) + gttmne(kka) / abs(kda) )
      fttcs = jthweight * ( - gttmno(kks) /     kds  + gttmno(kka) /     kda  )
      fttsc = jthweight * ( + gttmno(kks) /     kds  + gttmno(kka) /     kda  )
      fttss = jthweight * ( + gttmne(kks) / abs(kds) - gttmne(kka) / abs(kda) )
      
      ftzcc = jthweight * ( + gtzmne(kks) / abs(kds) + gtzmne(kka) / abs(kda) )
      ftzcs = jthweight * ( - gtzmno(kks) /     kds  + gtzmno(kka) /     kda  )
      ftzsc = jthweight * ( + gtzmno(kks) /     kds  + gtzmno(kka) /     kda  )
      ftzss = jthweight * ( + gtzmne(kks) / abs(kds) - gtzmne(kka) / abs(kda) )
      
      fzzcc = jthweight * ( + gzzmne(kks) / abs(kds) + gzzmne(kka) / abs(kda) )
      fzzcs = jthweight * ( - gzzmno(kks) /     kds  + gzzmno(kka) /     kda  )
      fzzsc = jthweight * ( + gzzmno(kks) /     kds  + gzzmno(kka) /     kda  )
      fzzss = jthweight * ( + gzzmne(kks) / abs(kds) - gzzmne(kka) / abs(kda) )
      
      do ll = 0, lrad
       
       Tl = cheby(ll,0)
       Dl = cheby(ll,1)
       
       do pp = 0, lrad
        
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
        
       enddo ! end of do pp;  1 Feb 13;
       
      enddo ! end of do ll;  1 Feb 13;
      
     enddo ! end of do jj;  1 Feb 13;
    enddo ! end of do ii;  1 Feb 13;
    
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
