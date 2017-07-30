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

!latex \subsection{chebyshev-metric information}

!latex \begin{enumerate}

!latex \item There are various symmetries that can be exploited.

!latex \item Most simply, \verb+DToocc+, \verb+DToocs+, \verb+DToosc+ and \verb+DTooss+ do not depend on geometry and need only be computed once.

!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine ma00aa( lquad, mn, lvol, lrad )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half, one, pi
  
  use numerical, only : vsmall
  
  use fileunits, only : ounit
  
  use inputlist, only : Wma00aa
  
  use cputiming, only : Tma00aa
  
  use allglobal, only : myid, ncpu, cpus, &
                        Mvol, mnsqd, ilabel, jlabel, llabel, plabel, &
                        gaussianweight, gaussianabscissae, &
                        DToocc, DToocs, DToosc, DTooss, &
                        TTsscc, TTsscs, TTsssc, TTssss, &
                        TDstcc, TDstcs, TDstsc, TDstss, &
                        TDszcc, TDszcs, TDszsc, TDszss, &
                        DDttcc, DDttcs, DDttsc, DDttss, &
                        DDtzcc, DDtzcs, DDtzsc, DDtzss, &
                        DDzzcc, DDzzcs, DDzzsc, DDzzss, &
                        kija, kijs, &
                        goomne, goomno, &
                        gssmne, gssmno, &
                        gstmne, gstmno, &
                        gszmne, gszmno, &
                        gttmne, gttmno, &
                        gtzmne, gtzmno, &
                        gzzmne, gzzmno, &
                        TD, &
                        Lcoordinatesingularity, regumm, &
                       !pi2pi2nfphalf, & ! SRH; 27 Jul 17;
                        YESstellsym
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in) :: lquad, mn, lvol, lrad

  LOGICAL             :: Lpause ! for debugging only; SRH; 30 Jul 17;
  
  INTEGER             :: jquad, ii, jj, ij, ll, pp, lp, kka, kks ! SRH; 27 Jul 17;
  REAL                :: lss, jthweight, fee, feo, foe, foo, Tl, Dl, Tp, Dp, TlTp, TlDp, DlTp, DlDp, kda, kds ! SRH; 27 Jul 17;

  REAL                :: foocc, foocs, foosc, fooss
  REAL                :: fsscc, fsscs, fsssc, fssss
  REAL                :: fstcc, fstcs, fstsc, fstss
  REAL                :: fszcc, fszcs, fszsc, fszss
  REAL                :: fttcc, fttcs, fttsc, fttss
  REAL                :: ftzcc, ftzcs, ftzsc, ftzss
  REAL                :: fzzcc, fzzcs, fzzsc, fzzss
  
  REAL                :: sbar(1:lquad), halfoversbar(1:lquad), sbarhim(1:lquad,1:mn) ! regularization factors required for coordinate singularity; 10 Dec 15;
  
  BEGIN( ma00aa )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  FATAL( ma00aa, lvol.lt.1 .or. lvol.gt.Mvol, illegal volume label ) ! its good to keep internal checks like this in the DEBUG version; SRH; 27 Jul 17;
#endif
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( YESstellsym ) then
   
   DToocc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero ! initialize summation of Gaussian quadrature (loop over jquad); SRH; 27 Jul 17;
  !DToocs( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero ! non-stellarator-symmetric terms are commented; SRH; 27 Jul 17;
  !DToosc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero ! I am happy to delete this commented-out source eventually; "a short code is a good code"; SRH; 27 Jul 17;
  !DTooss( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
   
  !TTsscc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
  !TTsscs( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
  !TTsssc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
   TTssss( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
   
  !TDstcc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
  !TDstcs( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
   TDstsc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
  !TDstss( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
   
  !TDszcc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
  !TDszcs( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
   TDszsc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
  !TDszss( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
   
   DDttcc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
  !DDttcs( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
  !DDttsc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
  !DDttss( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
   
   DDtzcc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
  !DDtzcs( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
  !DDtzsc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
  !DDtzss( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
   
   DDzzcc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
  !DDzzcs( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
  !DDzzsc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
  !DDzzss( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero
   
  else ! NOTstellsym ; SRH; 27 Jul 17;
   
   DToocc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = zero ! initialize summation of Gaussian quadrature (loop over jquad); SRH; 27 Jul 17;
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
   
  endif ! end of if( YESstellsym ) ; SRH; 27 Jul 17;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Lcoordinatesingularity ) then ! additional radial factors, such as r^m, are included to "regularize" the magnetic field near the origin; SRH; 27 Jul 17;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   sbar(1:lquad) = ( gaussianabscissae(1:lquad,lvol) + one ) * half
   
   halfoversbar(1:lquad) = half / sbar(1:lquad)
   
   do jquad = 1, lquad ; sbarhim(jquad,1:mn) = sbar(jquad)**regumm(1:mn) ! pre-calculation of regularization factor; 12 Sep 13;
   enddo
   
   if( YESstellsym ) then
    
    do jquad = 1, lquad ! Gaussian quadrature loop;
     
     lss = gaussianabscissae(jquad,lvol) ; jthweight = gaussianweight(jquad,lvol)
     
     WCALL( ma00aa, metrix,( lvol, lss ) ) ! compute metric elements; 16 Jan 13;
     
     do ij = 1, mnsqd ; ii = ilabel(ij) ; jj = jlabel(ij) ! SRH; 27 Jul 17;
      
      kks = kijs(ii,jj,0) ; kds = jthweight / kijs(ii,jj,1) ! SRH; 27 Jul 17;
      kka = kija(ii,jj,0) ; kda = jthweight / kija(ii,jj,1) ! SRH; 27 Jul 17;
      
      foocc = + goomne(kks) * abs(kds) + goomne(kka) * abs(kda) ! stell-sym; SRH; 27 Jul 17;
     !foocs = - goomno(kks) *     kds  + goomno(kka) *     kda 
     !foosc = + goomno(kks) *     kds  + goomno(kka) *     kda 
     !fooss = + goomne(kks) * abs(kds) - goomne(kka) * abs(kda)
      
     !fsscc = + gssmne(kks) * abs(kds) + gssmne(kka) * abs(kda)
     !fsscs = - gssmno(kks) *     kds  + gssmno(kka) *     kda 
     !fsssc = + gssmno(kks) *     kds  + gssmno(kka) *     kda 
      fssss = + gssmne(kks) * abs(kds) - gssmne(kka) * abs(kda) ! stell-sym; 
      
     !fstcc = + gstmne(kks) * abs(kds) + gstmne(kka) * abs(kda)
     !fstcs = - gstmno(kks) *     kds  + gstmno(kka) *     kda 
      fstsc = + gstmno(kks) *     kds  + gstmno(kka) *     kda  ! stell-sym; 
     !fstss = + gstmne(kks) * abs(kds) - gstmne(kka) * abs(kda)
      
     !fszcc = + gszmne(kks) * abs(kds) + gszmne(kka) * abs(kda)
     !fszcs = - gszmno(kks) *     kds  + gszmno(kka) *     kda 
      fszsc = + gszmno(kks) *     kds  + gszmno(kka) *     kda  ! stell-sym; 
     !fszss = + gszmne(kks) * abs(kds) - gszmne(kka) * abs(kda)
      
      fttcc = + gttmne(kks) * abs(kds) + gttmne(kka) * abs(kda) ! stell-sym; 
     !fttcs = - gttmno(kks) *     kds  + gttmno(kka) *     kda 
     !fttsc = + gttmno(kks) *     kds  + gttmno(kka) *     kda 
     !fttss = + gttmne(kks) * abs(kds) - gttmne(kka) * abs(kda)
      
      ftzcc = + gtzmne(kks) * abs(kds) + gtzmne(kka) * abs(kda) ! stell-sym; 
     !ftzcs = - gtzmno(kks) *     kds  + gtzmno(kka) *     kda 
     !ftzsc = + gtzmno(kks) *     kds  + gtzmno(kka) *     kda 
     !ftzss = + gtzmne(kks) * abs(kds) - gtzmne(kka) * abs(kda)
      
      fzzcc = + gzzmne(kks) * abs(kds) + gzzmne(kka) * abs(kda) ! stell-sym; 
     !fzzcs = - gzzmno(kks) *     kds  + gzzmno(kka) *     kda 
     !fzzsc = + gzzmno(kks) *     kds  + gzzmno(kka) *     kda 
     !fzzss = + gzzmne(kks) * abs(kds) - gzzmne(kka) * abs(kda)
      
      do lp = 1, (1+lrad)*(1+lrad) ; ll = llabel(lp,lvol) ; pp = plabel(lp,lvol) ! SRH; 27 Jul 17;
       
       Tl = sbarhim(jquad,ii) *                                      TD(ll,0,jquad,lvol)
       Dl = sbarhim(jquad,ii) * ( regumm(ii) * halfoversbar(jquad) * TD(ll,0,jquad,lvol) + TD(ll,1,jquad,lvol) )
       Tp = sbarhim(jquad,jj) *                                      TD(pp,0,jquad,lvol)
       Dp = sbarhim(jquad,jj) * ( regumm(jj) * halfoversbar(jquad) * TD(pp,0,jquad,lvol) + TD(pp,1,jquad,lvol) )
       
       TlTp = Tl * Tp
       TlDp = Tl * Dp
       DlTp = Dl * Tp
       DlDp = Dl * Dp
       
       DToocc( ll, pp, ii, jj ) = DToocc( ll, pp, ii, jj ) + DlTp * foocc ! stell-sym; SRH; 27 Jul 17;
      !DToocs( ll, pp, ii, jj ) = DToocs( ll, pp, ii, jj ) + DlTp * foocs
      !DToosc( ll, pp, ii, jj ) = DToosc( ll, pp, ii, jj ) + DlTp * foosc
      !DTooss( ll, pp, ii, jj ) = DTooss( ll, pp, ii, jj ) + DlTp * fooss
       
      !TTsscc( ll, pp, ii, jj ) = TTsscc( ll, pp, ii, jj ) + TlTp * fsscc
      !TTsscs( ll, pp, ii, jj ) = TTsscs( ll, pp, ii, jj ) + TlTp * fsscs
      !TTsssc( ll, pp, ii, jj ) = TTsssc( ll, pp, ii, jj ) + TlTp * fsssc
       TTssss( ll, pp, ii, jj ) = TTssss( ll, pp, ii, jj ) + TlTp * fssss ! stell-sym; SRH; 27 Jul 17;
       
      !TDstcc( ll, pp, ii, jj ) = TDstcc( ll, pp, ii, jj ) + TlDp * fstcc
      !TDstcs( ll, pp, ii, jj ) = TDstcs( ll, pp, ii, jj ) + TlDp * fstcs
       TDstsc( ll, pp, ii, jj ) = TDstsc( ll, pp, ii, jj ) + TlDp * fstsc ! stell-sym; SRH; 27 Jul 17;
      !TDstss( ll, pp, ii, jj ) = TDstss( ll, pp, ii, jj ) + TlDp * fstss
       
      !TDszcc( ll, pp, ii, jj ) = TDszcc( ll, pp, ii, jj ) + TlDp * fszcc
      !TDszcs( ll, pp, ii, jj ) = TDszcs( ll, pp, ii, jj ) + TlDp * fszcs
       TDszsc( ll, pp, ii, jj ) = TDszsc( ll, pp, ii, jj ) + TlDp * fszsc ! stell-sym; SRH; 27 Jul 17;
      !TDszss( ll, pp, ii, jj ) = TDszss( ll, pp, ii, jj ) + TlDp * fszss
       
       DDttcc( ll, pp, ii, jj ) = DDttcc( ll, pp, ii, jj ) + DlDp * fttcc ! stell-sym; SRH; 27 Jul 17;
      !DDttcs( ll, pp, ii, jj ) = DDttcs( ll, pp, ii, jj ) + DlDp * fttcs
      !DDttsc( ll, pp, ii, jj ) = DDttsc( ll, pp, ii, jj ) + DlDp * fttsc
      !DDttss( ll, pp, ii, jj ) = DDttss( ll, pp, ii, jj ) + DlDp * fttss
       
       DDtzcc( ll, pp, ii, jj ) = DDtzcc( ll, pp, ii, jj ) + DlDp * ftzcc ! stell-sym; SRH; 27 Jul 17;
      !DDtzcs( ll, pp, ii, jj ) = DDtzcs( ll, pp, ii, jj ) + DlDp * ftzcs
      !DDtzsc( ll, pp, ii, jj ) = DDtzsc( ll, pp, ii, jj ) + DlDp * ftzsc
      !DDtzss( ll, pp, ii, jj ) = DDtzss( ll, pp, ii, jj ) + DlDp * ftzss
       
       DDzzcc( ll, pp, ii, jj ) = DDzzcc( ll, pp, ii, jj ) + DlDp * fzzcc ! stell-sym; SRH; 27 Jul 17;
      !DDzzcs( ll, pp, ii, jj ) = DDzzcs( ll, pp, ii, jj ) + DlDp * fzzcs
      !DDzzsc( ll, pp, ii, jj ) = DDzzsc( ll, pp, ii, jj ) + DlDp * fzzsc
      !DDzzss( ll, pp, ii, jj ) = DDzzss( ll, pp, ii, jj ) + DlDp * fzzss
       
      enddo ! end of do lp ; 08 Feb 16;
      
     enddo ! end of do ij ; 08 Feb 16;
     
    enddo ! end of do jquad; ! 16 Jan 13;
    
   else ! NOTstellsym ; SRH; 27 Jul 17;
    
    do jquad = 1, lquad ! Gaussian quadrature loop;
     
     lss = gaussianabscissae(jquad,lvol) ; jthweight = gaussianweight(jquad,lvol)
     
     WCALL( ma00aa, metrix,( lvol, lss ) ) ! compute metric elements; 16 Jan 13;
     
     do ij = 1, mnsqd ; ii = ilabel(ij) ; jj = jlabel(ij) ! SRH; 27 Jul 17;
      
      kks = kijs(ii,jj,0) ; kds = jthweight / kijs(ii,jj,1) ! SRH; 27 Jul 17;
      kka = kija(ii,jj,0) ; kda = jthweight / kija(ii,jj,1) ! SRH; 27 Jul 17;
      
      foocc = + goomne(kks) * abs(kds) + goomne(kka) * abs(kda) ! SRH; 27 Jul 17;
      foocs = - goomno(kks) *     kds  + goomno(kka) *     kda 
      foosc = + goomno(kks) *     kds  + goomno(kka) *     kda 
      fooss = + goomne(kks) * abs(kds) - goomne(kka) * abs(kda)
      
      fsscc = + gssmne(kks) * abs(kds) + gssmne(kka) * abs(kda)
      fsscs = - gssmno(kks) *     kds  + gssmno(kka) *     kda 
      fsssc = + gssmno(kks) *     kds  + gssmno(kka) *     kda 
      fssss = + gssmne(kks) * abs(kds) - gssmne(kka) * abs(kda)
      
      fstcc = + gstmne(kks) * abs(kds) + gstmne(kka) * abs(kda)
      fstcs = - gstmno(kks) *     kds  + gstmno(kka) *     kda 
      fstsc = + gstmno(kks) *     kds  + gstmno(kka) *     kda 
      fstss = + gstmne(kks) * abs(kds) - gstmne(kka) * abs(kda)
      
      fszcc = + gszmne(kks) * abs(kds) + gszmne(kka) * abs(kda)
      fszcs = - gszmno(kks) *     kds  + gszmno(kka) *     kda 
      fszsc = + gszmno(kks) *     kds  + gszmno(kka) *     kda 
      fszss = + gszmne(kks) * abs(kds) - gszmne(kka) * abs(kda)
      
      fttcc = + gttmne(kks) * abs(kds) + gttmne(kka) * abs(kda)
      fttcs = - gttmno(kks) *     kds  + gttmno(kka) *     kda 
      fttsc = + gttmno(kks) *     kds  + gttmno(kka) *     kda 
      fttss = + gttmne(kks) * abs(kds) - gttmne(kka) * abs(kda)
      
      ftzcc = + gtzmne(kks) * abs(kds) + gtzmne(kka) * abs(kda)
      ftzcs = - gtzmno(kks) *     kds  + gtzmno(kka) *     kda 
      ftzsc = + gtzmno(kks) *     kds  + gtzmno(kka) *     kda 
      ftzss = + gtzmne(kks) * abs(kds) - gtzmne(kka) * abs(kda)
      
      fzzcc = + gzzmne(kks) * abs(kds) + gzzmne(kka) * abs(kda)
      fzzcs = - gzzmno(kks) *     kds  + gzzmno(kka) *     kda 
      fzzsc = + gzzmno(kks) *     kds  + gzzmno(kka) *     kda 
      fzzss = + gzzmne(kks) * abs(kds) - gzzmne(kka) * abs(kda)
      
      do lp = 1, (1+lrad)*(1+lrad) ; ll = llabel(lp,lvol) ; pp = plabel(lp,lvol) ! SRH; 27 Jul 17;
       
       Tl = sbarhim(jquad,ii) *                                      TD(ll,0,jquad,lvol)
       Dl = sbarhim(jquad,ii) * ( regumm(ii) * halfoversbar(jquad) * TD(ll,0,jquad,lvol) + TD(ll,1,jquad,lvol) )
       Tp = sbarhim(jquad,jj) *                                      TD(pp,0,jquad,lvol)
       Dp = sbarhim(jquad,jj) * ( regumm(jj) * halfoversbar(jquad) * TD(pp,0,jquad,lvol) + TD(pp,1,jquad,lvol) )
       
       TlTp = Tl * Tp
       TlDp = Tl * Dp
       DlTp = Dl * Tp
       DlDp = Dl * Dp
       
       DToocc( ll, pp, ii, jj ) = DToocc( ll, pp, ii, jj ) + DlTp * foocc ! stell-sym; SRH; 27 Jul 17;
       DToocs( ll, pp, ii, jj ) = DToocs( ll, pp, ii, jj ) + DlTp * foocs
       DToosc( ll, pp, ii, jj ) = DToosc( ll, pp, ii, jj ) + DlTp * foosc
       DTooss( ll, pp, ii, jj ) = DTooss( ll, pp, ii, jj ) + DlTp * fooss
       
       TTsscc( ll, pp, ii, jj ) = TTsscc( ll, pp, ii, jj ) + TlTp * fsscc
       TTsscs( ll, pp, ii, jj ) = TTsscs( ll, pp, ii, jj ) + TlTp * fsscs
       TTsssc( ll, pp, ii, jj ) = TTsssc( ll, pp, ii, jj ) + TlTp * fsssc
       TTssss( ll, pp, ii, jj ) = TTssss( ll, pp, ii, jj ) + TlTp * fssss ! stell-sym; SRH; 27 Jul 17;
       
       TDstcc( ll, pp, ii, jj ) = TDstcc( ll, pp, ii, jj ) + TlDp * fstcc
       TDstcs( ll, pp, ii, jj ) = TDstcs( ll, pp, ii, jj ) + TlDp * fstcs
       TDstsc( ll, pp, ii, jj ) = TDstsc( ll, pp, ii, jj ) + TlDp * fstsc ! stell-sym; SRH; 27 Jul 17;
       TDstss( ll, pp, ii, jj ) = TDstss( ll, pp, ii, jj ) + TlDp * fstss
       
       TDszcc( ll, pp, ii, jj ) = TDszcc( ll, pp, ii, jj ) + TlDp * fszcc
       TDszcs( ll, pp, ii, jj ) = TDszcs( ll, pp, ii, jj ) + TlDp * fszcs
       TDszsc( ll, pp, ii, jj ) = TDszsc( ll, pp, ii, jj ) + TlDp * fszsc ! stell-sym; SRH; 27 Jul 17;
       TDszss( ll, pp, ii, jj ) = TDszss( ll, pp, ii, jj ) + TlDp * fszss
       
       DDttcc( ll, pp, ii, jj ) = DDttcc( ll, pp, ii, jj ) + DlDp * fttcc ! stell-sym; SRH; 27 Jul 17;
       DDttcs( ll, pp, ii, jj ) = DDttcs( ll, pp, ii, jj ) + DlDp * fttcs
       DDttsc( ll, pp, ii, jj ) = DDttsc( ll, pp, ii, jj ) + DlDp * fttsc
       DDttss( ll, pp, ii, jj ) = DDttss( ll, pp, ii, jj ) + DlDp * fttss
       
       DDtzcc( ll, pp, ii, jj ) = DDtzcc( ll, pp, ii, jj ) + DlDp * ftzcc ! stell-sym; SRH; 27 Jul 17;
       DDtzcs( ll, pp, ii, jj ) = DDtzcs( ll, pp, ii, jj ) + DlDp * ftzcs
       DDtzsc( ll, pp, ii, jj ) = DDtzsc( ll, pp, ii, jj ) + DlDp * ftzsc
       DDtzss( ll, pp, ii, jj ) = DDtzss( ll, pp, ii, jj ) + DlDp * ftzss
       
       DDzzcc( ll, pp, ii, jj ) = DDzzcc( ll, pp, ii, jj ) + DlDp * fzzcc ! stell-sym; SRH; 27 Jul 17;
       DDzzcs( ll, pp, ii, jj ) = DDzzcs( ll, pp, ii, jj ) + DlDp * fzzcs
       DDzzsc( ll, pp, ii, jj ) = DDzzsc( ll, pp, ii, jj ) + DlDp * fzzsc
       DDzzss( ll, pp, ii, jj ) = DDzzss( ll, pp, ii, jj ) + DlDp * fzzss
       
      enddo ! end of do lp ; 08 Feb 16;
      
     enddo ! end of do ij ; 08 Feb 16;
     
    enddo ! end of do jquad; ! 16 Jan 13;
    
   endif ! end of if( YESstellsym) ; 
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
  else ! .not.Lcoordinatesingularity; 17 Dec 15;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   if( YESstellsym ) then
    
    do jquad = 1, lquad ! Gaussian quadrature loop;
     
     lss = gaussianabscissae(jquad,lvol) ; jthweight = gaussianweight(jquad,lvol)
     
     WCALL( ma00aa, metrix,( lvol, lss ) ) ! compute metric elements; 16 Jan 13;
     
     do ij = 1, mnsqd ; ii = ilabel(ij) ; jj = jlabel(ij) ! SRH; 27 Jul 17;
      
      kks = kijs(ii,jj,0) ; kds = jthweight / kijs(ii,jj,1) ! SRH; 27 Jul 17;
      kka = kija(ii,jj,0) ; kda = jthweight / kija(ii,jj,1) ! SRH; 27 Jul 17;
      
      foocc = + goomne(kks) * abs(kds) + goomne(kka) * abs(kda) ! stell-sym; SRH; 27 Jul 17;
     !foocs = - goomno(kks) *     kds  + goomno(kka) *     kda 
     !foosc = + goomno(kks) *     kds  + goomno(kka) *     kda 
     !fooss = + goomne(kks) * abs(kds) - goomne(kka) * abs(kda)
      
     !fsscc = + gssmne(kks) * abs(kds) + gssmne(kka) * abs(kda)
     !fsscs = - gssmno(kks) *     kds  + gssmno(kka) *     kda 
     !fsssc = + gssmno(kks) *     kds  + gssmno(kka) *     kda 
      fssss = + gssmne(kks) * abs(kds) - gssmne(kka) * abs(kda) ! stell-sym; 
      
     !fstcc = + gstmne(kks) * abs(kds) + gstmne(kka) * abs(kda)
     !fstcs = - gstmno(kks) *     kds  + gstmno(kka) *     kda 
      fstsc = + gstmno(kks) *     kds  + gstmno(kka) *     kda  ! stell-sym; 
     !fstss = + gstmne(kks) * abs(kds) - gstmne(kka) * abs(kda)
      
     !fszcc = + gszmne(kks) * abs(kds) + gszmne(kka) * abs(kda)
     !fszcs = - gszmno(kks) *     kds  + gszmno(kka) *     kda 
      fszsc = + gszmno(kks) *     kds  + gszmno(kka) *     kda  ! stell-sym; 
     !fszss = + gszmne(kks) * abs(kds) - gszmne(kka) * abs(kda)
      
      fttcc = + gttmne(kks) * abs(kds) + gttmne(kka) * abs(kda) ! stell-sym; 
     !fttcs = - gttmno(kks) *     kds  + gttmno(kka) *     kda 
     !fttsc = + gttmno(kks) *     kds  + gttmno(kka) *     kda 
     !fttss = + gttmne(kks) * abs(kds) - gttmne(kka) * abs(kda)
      
      ftzcc = + gtzmne(kks) * abs(kds) + gtzmne(kka) * abs(kda) ! stell-sym; 
     !ftzcs = - gtzmno(kks) *     kds  + gtzmno(kka) *     kda 
     !ftzsc = + gtzmno(kks) *     kds  + gtzmno(kka) *     kda 
     !ftzss = + gtzmne(kks) * abs(kds) - gtzmne(kka) * abs(kda)
      
      fzzcc = + gzzmne(kks) * abs(kds) + gzzmne(kka) * abs(kda) ! stell-sym; 
     !fzzcs = - gzzmno(kks) *     kds  + gzzmno(kka) *     kda 
     !fzzsc = + gzzmno(kks) *     kds  + gzzmno(kka) *     kda 
     !fzzss = + gzzmne(kks) * abs(kds) - gzzmne(kka) * abs(kda)
      
      do lp = 1, (1+lrad)*(1+lrad) ; ll = llabel(lp,lvol) ; pp = plabel(lp,lvol) ! SRH; 27 Jul 17;
       
       Tl = TD(ll,0,jquad,lvol)
       Dl = TD(ll,1,jquad,lvol)
       Tp = TD(pp,0,jquad,lvol)
       Dp = TD(pp,1,jquad,lvol)
       
       TlTp = Tl * Tp
       TlDp = Tl * Dp
       DlTp = Dl * Tp
       DlDp = Dl * Dp
       
       DToocc( ll, pp, ii, jj ) = DToocc( ll, pp, ii, jj ) + DlTp * foocc ! stell-sym; SRH; 27 Jul 17;
      !DToocs( ll, pp, ii, jj ) = DToocs( ll, pp, ii, jj ) + DlTp * foocs
      !DToosc( ll, pp, ii, jj ) = DToosc( ll, pp, ii, jj ) + DlTp * foosc
      !DTooss( ll, pp, ii, jj ) = DTooss( ll, pp, ii, jj ) + DlTp * fooss
       
      !TTsscc( ll, pp, ii, jj ) = TTsscc( ll, pp, ii, jj ) + TlTp * fsscc
      !TTsscs( ll, pp, ii, jj ) = TTsscs( ll, pp, ii, jj ) + TlTp * fsscs
      !TTsssc( ll, pp, ii, jj ) = TTsssc( ll, pp, ii, jj ) + TlTp * fsssc
       TTssss( ll, pp, ii, jj ) = TTssss( ll, pp, ii, jj ) + TlTp * fssss ! stell-sym; SRH; 27 Jul 17;
       
      !TDstcc( ll, pp, ii, jj ) = TDstcc( ll, pp, ii, jj ) + TlDp * fstcc
      !TDstcs( ll, pp, ii, jj ) = TDstcs( ll, pp, ii, jj ) + TlDp * fstcs
       TDstsc( ll, pp, ii, jj ) = TDstsc( ll, pp, ii, jj ) + TlDp * fstsc ! stell-sym; SRH; 27 Jul 17;
      !TDstss( ll, pp, ii, jj ) = TDstss( ll, pp, ii, jj ) + TlDp * fstss
       
      !TDszcc( ll, pp, ii, jj ) = TDszcc( ll, pp, ii, jj ) + TlDp * fszcc
      !TDszcs( ll, pp, ii, jj ) = TDszcs( ll, pp, ii, jj ) + TlDp * fszcs
       TDszsc( ll, pp, ii, jj ) = TDszsc( ll, pp, ii, jj ) + TlDp * fszsc ! stell-sym; SRH; 27 Jul 17;
      !TDszss( ll, pp, ii, jj ) = TDszss( ll, pp, ii, jj ) + TlDp * fszss
       
       DDttcc( ll, pp, ii, jj ) = DDttcc( ll, pp, ii, jj ) + DlDp * fttcc ! stell-sym; SRH; 27 Jul 17;
      !DDttcs( ll, pp, ii, jj ) = DDttcs( ll, pp, ii, jj ) + DlDp * fttcs
      !DDttsc( ll, pp, ii, jj ) = DDttsc( ll, pp, ii, jj ) + DlDp * fttsc
      !DDttss( ll, pp, ii, jj ) = DDttss( ll, pp, ii, jj ) + DlDp * fttss
       
       DDtzcc( ll, pp, ii, jj ) = DDtzcc( ll, pp, ii, jj ) + DlDp * ftzcc ! stell-sym; SRH; 27 Jul 17;
      !DDtzcs( ll, pp, ii, jj ) = DDtzcs( ll, pp, ii, jj ) + DlDp * ftzcs
      !DDtzsc( ll, pp, ii, jj ) = DDtzsc( ll, pp, ii, jj ) + DlDp * ftzsc
      !DDtzss( ll, pp, ii, jj ) = DDtzss( ll, pp, ii, jj ) + DlDp * ftzss
       
       DDzzcc( ll, pp, ii, jj ) = DDzzcc( ll, pp, ii, jj ) + DlDp * fzzcc ! stell-sym; SRH; 27 Jul 17;
      !DDzzcs( ll, pp, ii, jj ) = DDzzcs( ll, pp, ii, jj ) + DlDp * fzzcs
      !DDzzsc( ll, pp, ii, jj ) = DDzzsc( ll, pp, ii, jj ) + DlDp * fzzsc
      !DDzzss( ll, pp, ii, jj ) = DDzzss( ll, pp, ii, jj ) + DlDp * fzzss
       
      enddo ! end of do lp; 27 Jul 17;
      
     enddo ! end of do ij; 27 Jul 17;
     
    enddo ! end of do jquad; ! 16 Jan 13;
    
   else ! NOTstellsym
    
    do jquad = 1, lquad ! Gaussian quadrature loop;
     
     lss = gaussianabscissae(jquad,lvol) ; jthweight = gaussianweight(jquad,lvol)
     
     WCALL( ma00aa, metrix,( lvol, lss ) ) ! compute metric elements; 16 Jan 13;
     
     do ij = 1, mnsqd ; ii = ilabel(ij) ; jj = jlabel(ij) ! SRH; 27 Jul 17;
      
      kks = kijs(ii,jj,0) ; kds = jthweight / kijs(ii,jj,1) ! SRH; 27 Jul 17;
      kka = kija(ii,jj,0) ; kda = jthweight / kija(ii,jj,1) ! SRH; 27 Jul 17;
      
      foocc = + goomne(kks) * abs(kds) + goomne(kka) * abs(kda) ! SRH; 27 Jul 17;
      foocs = - goomno(kks) *     kds  + goomno(kka) *     kda 
      foosc = + goomno(kks) *     kds  + goomno(kka) *     kda 
      fooss = + goomne(kks) * abs(kds) - goomne(kka) * abs(kda)
      
      fsscc = + gssmne(kks) * abs(kds) + gssmne(kka) * abs(kda)
      fsscs = - gssmno(kks) *     kds  + gssmno(kka) *     kda 
      fsssc = + gssmno(kks) *     kds  + gssmno(kka) *     kda 
      fssss = + gssmne(kks) * abs(kds) - gssmne(kka) * abs(kda)
      
      fstcc = + gstmne(kks) * abs(kds) + gstmne(kka) * abs(kda)
      fstcs = - gstmno(kks) *     kds  + gstmno(kka) *     kda 
      fstsc = + gstmno(kks) *     kds  + gstmno(kka) *     kda 
      fstss = + gstmne(kks) * abs(kds) - gstmne(kka) * abs(kda)
      
      fszcc = + gszmne(kks) * abs(kds) + gszmne(kka) * abs(kda)
      fszcs = - gszmno(kks) *     kds  + gszmno(kka) *     kda 
      fszsc = + gszmno(kks) *     kds  + gszmno(kka) *     kda 
      fszss = + gszmne(kks) * abs(kds) - gszmne(kka) * abs(kda)
      
      fttcc = + gttmne(kks) * abs(kds) + gttmne(kka) * abs(kda)
      fttcs = - gttmno(kks) *     kds  + gttmno(kka) *     kda 
      fttsc = + gttmno(kks) *     kds  + gttmno(kka) *     kda 
      fttss = + gttmne(kks) * abs(kds) - gttmne(kka) * abs(kda)
      
      ftzcc = + gtzmne(kks) * abs(kds) + gtzmne(kka) * abs(kda)
      ftzcs = - gtzmno(kks) *     kds  + gtzmno(kka) *     kda 
      ftzsc = + gtzmno(kks) *     kds  + gtzmno(kka) *     kda 
      ftzss = + gtzmne(kks) * abs(kds) - gtzmne(kka) * abs(kda)
      
      fzzcc = + gzzmne(kks) * abs(kds) + gzzmne(kka) * abs(kda)
      fzzcs = - gzzmno(kks) *     kds  + gzzmno(kka) *     kda 
      fzzsc = + gzzmno(kks) *     kds  + gzzmno(kka) *     kda 
      fzzss = + gzzmne(kks) * abs(kds) - gzzmne(kka) * abs(kda)
      
      do lp = 1, (1+lrad)*(1+lrad) ; ll = llabel(lp,lvol) ; pp = plabel(lp,lvol) ! SRH; 27 Jul 17;
       
       Tl = TD(ll,0,jquad,lvol)
       Dl = TD(ll,1,jquad,lvol)
       Tp = TD(pp,0,jquad,lvol)
       Dp = TD(pp,1,jquad,lvol)
       
       TlTp = Tl * Tp
       TlDp = Tl * Dp
       DlTp = Dl * Tp
       DlDp = Dl * Dp
       
       DToocc( ll, pp, ii, jj ) = DToocc( ll, pp, ii, jj ) + DlTp * foocc ! stell-sym; SRH; 27 Jul 17;
       DToocs( ll, pp, ii, jj ) = DToocs( ll, pp, ii, jj ) + DlTp * foocs
       DToosc( ll, pp, ii, jj ) = DToosc( ll, pp, ii, jj ) + DlTp * foosc
       DTooss( ll, pp, ii, jj ) = DTooss( ll, pp, ii, jj ) + DlTp * fooss
       
       TTsscc( ll, pp, ii, jj ) = TTsscc( ll, pp, ii, jj ) + TlTp * fsscc
       TTsscs( ll, pp, ii, jj ) = TTsscs( ll, pp, ii, jj ) + TlTp * fsscs
       TTsssc( ll, pp, ii, jj ) = TTsssc( ll, pp, ii, jj ) + TlTp * fsssc
       TTssss( ll, pp, ii, jj ) = TTssss( ll, pp, ii, jj ) + TlTp * fssss ! stell-sym; SRH; 27 Jul 17;
       
       TDstcc( ll, pp, ii, jj ) = TDstcc( ll, pp, ii, jj ) + TlDp * fstcc
       TDstcs( ll, pp, ii, jj ) = TDstcs( ll, pp, ii, jj ) + TlDp * fstcs
       TDstsc( ll, pp, ii, jj ) = TDstsc( ll, pp, ii, jj ) + TlDp * fstsc ! stell-sym; SRH; 27 Jul 17;
       TDstss( ll, pp, ii, jj ) = TDstss( ll, pp, ii, jj ) + TlDp * fstss
       
       TDszcc( ll, pp, ii, jj ) = TDszcc( ll, pp, ii, jj ) + TlDp * fszcc
       TDszcs( ll, pp, ii, jj ) = TDszcs( ll, pp, ii, jj ) + TlDp * fszcs
       TDszsc( ll, pp, ii, jj ) = TDszsc( ll, pp, ii, jj ) + TlDp * fszsc ! stell-sym; SRH; 27 Jul 17;
       TDszss( ll, pp, ii, jj ) = TDszss( ll, pp, ii, jj ) + TlDp * fszss
       
       DDttcc( ll, pp, ii, jj ) = DDttcc( ll, pp, ii, jj ) + DlDp * fttcc ! stell-sym; SRH; 27 Jul 17;
       DDttcs( ll, pp, ii, jj ) = DDttcs( ll, pp, ii, jj ) + DlDp * fttcs
       DDttsc( ll, pp, ii, jj ) = DDttsc( ll, pp, ii, jj ) + DlDp * fttsc
       DDttss( ll, pp, ii, jj ) = DDttss( ll, pp, ii, jj ) + DlDp * fttss
       
       DDtzcc( ll, pp, ii, jj ) = DDtzcc( ll, pp, ii, jj ) + DlDp * ftzcc ! stell-sym; SRH; 27 Jul 17;
       DDtzcs( ll, pp, ii, jj ) = DDtzcs( ll, pp, ii, jj ) + DlDp * ftzcs
       DDtzsc( ll, pp, ii, jj ) = DDtzsc( ll, pp, ii, jj ) + DlDp * ftzsc
       DDtzss( ll, pp, ii, jj ) = DDtzss( ll, pp, ii, jj ) + DlDp * ftzss
       
       DDzzcc( ll, pp, ii, jj ) = DDzzcc( ll, pp, ii, jj ) + DlDp * fzzcc ! stell-sym; SRH; 27 Jul 17;
       DDzzcs( ll, pp, ii, jj ) = DDzzcs( ll, pp, ii, jj ) + DlDp * fzzcs
       DDzzsc( ll, pp, ii, jj ) = DDzzsc( ll, pp, ii, jj ) + DlDp * fzzsc
       DDzzss( ll, pp, ii, jj ) = DDzzss( ll, pp, ii, jj ) + DlDp * fzzss
       
      enddo ! end of do lp; 27 Jul 17;
      
     enddo ! end of do ij; 27 Jul 17;
     
    enddo ! end of do jquad; ! 16 Jan 13;
    
   endif ! end of if( YESstellsym ) ; 
   
  endif ! end of if( Lcoordinatesingularity ) ; 17 Dec 15;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!  DToocc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = DToocc( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf ! this factor is now included in matrix; SRH; 27 Jul 17;
!  DToocs( 0:lrad, 0:lrad, 1:mn, 1:mn ) = DToocs( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf ! I guessed that this would be faster; SRH; 27 Jul 17;
!  DToosc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = DToosc( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf ! but I did not confirm; SRH; 27 Jul 17;
!  DTooss( 0:lrad, 0:lrad, 1:mn, 1:mn ) = DTooss( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
!
!  TTsscc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = TTsscc( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
!  TTsscs( 0:lrad, 0:lrad, 1:mn, 1:mn ) = TTsscs( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
!  TTsssc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = TTsssc( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
!  TTssss( 0:lrad, 0:lrad, 1:mn, 1:mn ) = TTssss( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
!
!  TDstcc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = TDstcc( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
!  TDstcs( 0:lrad, 0:lrad, 1:mn, 1:mn ) = TDstcs( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
!  TDstsc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = TDstsc( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
!  TDstss( 0:lrad, 0:lrad, 1:mn, 1:mn ) = TDstss( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
!
!  TDszcc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = TDszcc( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
!  TDszcs( 0:lrad, 0:lrad, 1:mn, 1:mn ) = TDszcs( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
!  TDszsc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = TDszsc( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
!  TDszss( 0:lrad, 0:lrad, 1:mn, 1:mn ) = TDszss( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
!
!  DDttcc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = DDttcc( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
!  DDttcs( 0:lrad, 0:lrad, 1:mn, 1:mn ) = DDttcs( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
!  DDttsc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = DDttsc( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
!  DDttss( 0:lrad, 0:lrad, 1:mn, 1:mn ) = DDttss( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
!
!  DDtzcc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = DDtzcc( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
!  DDtzcs( 0:lrad, 0:lrad, 1:mn, 1:mn ) = DDtzcs( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
!  DDtzsc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = DDtzsc( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
!  DDtzss( 0:lrad, 0:lrad, 1:mn, 1:mn ) = DDtzss( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
!
!  DDzzcc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = DDzzcc( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
!  DDzzcs( 0:lrad, 0:lrad, 1:mn, 1:mn ) = DDzzcs( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
!  DDzzsc( 0:lrad, 0:lrad, 1:mn, 1:mn ) = DDzzsc( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
!  DDzzss( 0:lrad, 0:lrad, 1:mn, 1:mn ) = DDzzss( 0:lrad, 0:lrad, 1:mn, 1:mn ) * pi2pi2nfphalf
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  
  if( Wma00aa ) then ! check symmetries; SRH; 30 Jul 17;
   
   if( YESstellsym ) then
    
    do ll = 0, lrad
     do pp = 0, lrad
      Lpause = .false.
      do ii = 1, mn
       do jj = 1, mn
        if( abs(TTssss(ll,pp,ii,jj)) .gt. vsmall ) then
         write(ounit,1000) myid, lvol, ll, pp, ii, jj, TTssss(ll,pp,ii,jj), TTssss(pp,ll,jj,ii), TTssss(ll,pp,ii,jj)-TTssss(pp,ll,jj,ii)
         FATAL( ma00aa, abs(TTssss(ll,pp,ii,jj)-TTssss(pp,ll,jj,ii)).gt.vsmall ), symmetry error )
        !Lpause = .true.
        endif
       enddo
      enddo
     !if( Lpause ) pause
     enddo
    enddo
    
   else ! NOTstellsym; SRH; 30 Jul 17;
    
    do ll = 0, lrad
     do pp = 0, lrad
      do ii = 1, mn
       do jj = 1, mn
       !write(ounit,1000) myid, lvol, ll, pp, ii, jj, DToocc(ll,pp,ii,jj), DToocc(pp,ll,jj,ii), DTooss(ll,pp,ii,jj), DTooss(pp,ll,jj,ii)
       enddo
      enddo
     enddo
    enddo
    
   endif ! end of if( YESstelsym ) ; SRH; 30 Jul 17;

  endif ! end of if( Wma00aa ) ; SRH; 30 Jul 17;

1000 format("ma00aa : " 10x " : myid ="i3" ; lvol ="i3" ; ll ="i3" ; pp ="i3" ; ii ="i3" ; jj ="i3" ; TTssss ="2es23.15" ; error ="es10.2" ;")

#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN( ma00aa )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine ma00aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
