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
  
  use inputlist, only : mpol, Wma00aa
  
  use cputiming, only : Tma00aa
  
  use allglobal, only : myid, ncpu, cpus, &
                        Mvol, im, in, mne, &
                        YESstellsym, NOTstellsym, &
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
                        Lcoordinatesingularity, regumm, &
                        pi2pi2nfp, pi2pi2nfphalf
                        
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in) :: lquad, mn, lvol, lrad
  
  INTEGER             :: jquad, ll, pp, ll1, pp1, uv, ii, jj, io, mn2, lp2, mn2_max, lp2_max, nele 
  
  INTEGER             :: kk, kd, kka, kks, kda, kds
  
  REAL                :: lss, jthweight, fee, feo, foe, foo, Tl, Dl, Tp, Dp, TlTp, TlDp, DlTp, DlDp, ikda, ikds, imn2, ilrad, lssm

  REAL                :: foocc, foocs, foosc, fooss
  REAL                :: fsscc, fsscs, fsssc, fssss
  REAL                :: fstcc, fstcs, fstsc, fstss
  REAL                :: fszcc, fszcs, fszsc, fszss
  REAL                :: fttcc, fttcs, fttsc, fttss
  REAL                :: ftzcc, ftzcs, ftzsc, ftzss
  REAL                :: fzzcc, fzzcs, fzzsc, fzzss
  
  REAL                :: sbar
  REAL, allocatable   :: basis(:,:,:)
  
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

  DToocc = zero
  TTssss = zero
  TDstsc = zero
  TDszsc = zero
  DDttcc = zero
  DDtzcc = zero
  DDzzcc = zero

  if (NOTstellsym) then
    DToocs = zero
    DToosc = zero
    DTooss = zero

    TTsscc = zero
    TTsscs = zero
    TTsssc = zero

    TDstcc = zero
    TDstcs = zero
    TDstss = zero

    TDszcc = zero
    TDszcs = zero
    TDszss = zero

    DDttcs = zero
    DDttsc = zero
    DDttss = zero

    DDtzcs = zero
    DDtzsc = zero
    DDtzss = zero

    DDzzcs = zero
    DDzzsc = zero
    DDzzss = zero
  endif !NOTstellsym
  
  SALLOCATE(basis, (0:lrad,0:mpol,0:1), zero)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
  do jquad = 1, lquad ! Gaussian quadrature loop;
    
    lss = gaussianabscissae(jquad,lvol) ; jthweight = gaussianweight(jquad,lvol)
    sbar = (lss + one) * half
    
    if (Lcoordinatesingularity) then
      call get_zernike(sbar, lrad, mpol, basis(:,:,0:1)) ! use Zernike polynomials 29 Jun 19;
    else
      call get_cheby(lss, lrad, basis(:,0,0:1))
    endif

    WCALL( ma00aa, metrix,( lvol, lss ) ) ! compute metric elements; 16 Jan 13;

    do mn2 = 1, mn2_max
      ii = mod(mn2-1,mn)+1
      jj = (mn2-ii) / mn + 1
      
      kks = kijs(ii,jj,0) !; kds = kijs(ii,jj,1) 
      kka = kija(ii,jj,0) !; kda = kija(ii,jj,1) 
      ikds = jthweight / kijs(ii,jj,1)
      ikda = jthweight / kija(ii,jj,1)

      foocc = ( + goomne(kks) * abs(ikds) + goomne(kka) * abs(ikda) )
      fssss = ( + gssmne(kks) * abs(ikds) - gssmne(kka) * abs(ikda) )
      fstsc = ( + gstmno(kks) *     ikds  + gstmno(kka) *     ikda  )
      fszsc = ( + gszmno(kks) *     ikds  + gszmno(kka) *     ikda  )
      fttcc = ( + gttmne(kks) * abs(ikds) + gttmne(kka) * abs(ikda) )
      ftzcc = ( + gtzmne(kks) * abs(ikds) + gtzmne(kka) * abs(ikda) )
      fzzcc = ( + gzzmne(kks) * abs(ikds) + gzzmne(kka) * abs(ikda) )

      if (NOTstellsym) then
        foocs = ( - goomno(kks) *     ikds  + goomno(kka) *     ikda  )
        foosc = ( + goomno(kks) *     ikds  + goomno(kka) *     ikda  )
        fooss = ( + goomne(kks) * abs(ikds) - goomne(kka) * abs(ikda) )

        fsscc = ( + gssmne(kks) * abs(ikds) + gssmne(kka) * abs(ikda) )
        fsscs = ( - gssmno(kks) *     ikds  + gssmno(kka) *     ikda  )
        fsssc = ( + gssmno(kks) *     ikds  + gssmno(kka) *     ikda  )

        fstcc = ( + gstmne(kks) * abs(ikds) + gstmne(kka) * abs(ikda) )
        fstcs = ( - gstmno(kks) *     ikds  + gstmno(kka) *     ikda  )
        fstss = ( + gstmne(kks) * abs(ikds) - gstmne(kka) * abs(ikda) )

        fszcc = ( + gszmne(kks) * abs(ikds) + gszmne(kka) * abs(ikda) )
        fszcs = ( - gszmno(kks) *     ikds  + gszmno(kka) *     ikda  )
        fszss = ( + gszmne(kks) * abs(ikds) - gszmne(kka) * abs(ikda) )

        fttcs = ( - gttmno(kks) *     ikds  + gttmno(kka) *     ikda  )
        fttsc = ( + gttmno(kks) *     ikds  + gttmno(kka) *     ikda  )
        fttss = ( + gttmne(kks) * abs(ikds) - gttmne(kka) * abs(ikda) )

        ftzcs = ( - gtzmno(kks) *     ikds  + gtzmno(kka) *     ikda  )
        ftzsc = ( + gtzmno(kks) *     ikds  + gtzmno(kka) *     ikda  )
        ftzss = ( + gtzmne(kks) * abs(ikds) - gtzmne(kka) * abs(ikda) )

        fzzcs = ( - gzzmno(kks) *     ikds  + gzzmno(kka) *     ikda  )
        fzzsc = ( + gzzmno(kks) *     ikds  + gzzmno(kka) *     ikda  )
        fzzss = ( + gzzmne(kks) * abs(ikds) - gzzmne(kka) * abs(ikda) )
      end if !NOTstellsym

      do lp2 = 1, lp2_max 
        ll = mod(lp2-1,lrad+1)
        pp = (lp2-ll-1)/(lrad+1)

        if (Lcoordinatesingularity) then

          ll1 = (ll - mod(ll,2))/2 ! shrinked dof for Zernike; 02 Jul 19
          pp1 = (pp - mod(pp,2))/2 ! shrinked dof for Zernike; 02 Jul 19

          if (ll < im(ii)) cycle ! zernike only non-zero for ll>=ii
          if (pp < im(jj)) cycle ! zernike only non-zero for pp>=jj
          if (mod(ll+im(ii),2)/=0) cycle ! zernike only non-zero if ll and ii have the same parity
          if (mod(pp+im(jj),2)/=0) cycle ! zernike only non-zero if pp and jj have the same parity

          Tl = basis(ll, im(ii), 0)         ! use Zernike polynomials 29 Jun 19;
          Dl = basis(ll, im(ii), 1) * half  ! use Zernike polynomials 29 Jun 19;
            
          Tp = basis(pp, im(jj), 0)         ! use Zernike polynomials 29 Jun 19;
          Dp = basis(pp, im(jj), 1) * half  ! use Zernike polynomials 29 Jun 19;
        
        else

          ll1 = ll
          pp1 = pp

          Tl = basis(ll, 0, 0)
          Dl = basis(ll, 0, 1)
            
          Tp = basis(pp, 0, 0)
          Dp = basis(pp, 0, 1)

        endif ! Lcoordinatesingularity

        TlTp = Tl * Tp
        TlDp = Tl * Dp
        DlTp = Dl * Tp
        DlDp = Dl * Dp 

        DToocc( ll1, pp1, ii, jj ) = DToocc( ll1, pp1, ii, jj ) + DlTp * foocc
        TTssss( ll1, pp1, ii, jj ) = TTssss( ll1, pp1, ii, jj ) + TlTp * fssss
        TDstsc( ll1, pp1, ii, jj ) = TDstsc( ll1, pp1, ii, jj ) + TlDp * fstsc
        TDszsc( ll1, pp1, ii, jj ) = TDszsc( ll1, pp1, ii, jj ) + TlDp * fszsc
        DDttcc( ll1, pp1, ii, jj ) = DDttcc( ll1, pp1, ii, jj ) + DlDp * fttcc
        DDtzcc( ll1, pp1, ii, jj ) = DDtzcc( ll1, pp1, ii, jj ) + DlDp * ftzcc
        DDzzcc( ll1, pp1, ii, jj ) = DDzzcc( ll1, pp1, ii, jj ) + DlDp * fzzcc

        if (NOTstellsym) then
          DToocs( ll1, pp1, ii, jj ) = DToocs( ll1, pp1, ii, jj ) + DlTp * foocs
          DToosc( ll1, pp1, ii, jj ) = DToosc( ll1, pp1, ii, jj ) + DlTp * foosc
          DTooss( ll1, pp1, ii, jj ) = DTooss( ll1, pp1, ii, jj ) + DlTp * fooss

          TTsscc( ll1, pp1, ii, jj ) = TTsscc( ll1, pp1, ii, jj ) + TlTp * fsscc
          TTsscs( ll1, pp1, ii, jj ) = TTsscs( ll1, pp1, ii, jj ) + TlTp * fsscs
          TTsssc( ll1, pp1, ii, jj ) = TTsssc( ll1, pp1, ii, jj ) + TlTp * fsssc

          TDstcc( ll1, pp1, ii, jj ) = TDstcc( ll1, pp1, ii, jj ) + TlDp * fstcc
          TDstcs( ll1, pp1, ii, jj ) = TDstcs( ll1, pp1, ii, jj ) + TlDp * fstcs
          TDstss( ll1, pp1, ii, jj ) = TDstss( ll1, pp1, ii, jj ) + TlDp * fstss

          TDszcc( ll1, pp1, ii, jj ) = TDszcc( ll1, pp1, ii, jj ) + TlDp * fszcc
          TDszcs( ll1, pp1, ii, jj ) = TDszcs( ll1, pp1, ii, jj ) + TlDp * fszcs
          TDszss( ll1, pp1, ii, jj ) = TDszss( ll1, pp1, ii, jj ) + TlDp * fszss

          DDttcs( ll1, pp1, ii, jj ) = DDttcs( ll1, pp1, ii, jj ) + DlDp * fttcs
          DDttsc( ll1, pp1, ii, jj ) = DDttsc( ll1, pp1, ii, jj ) + DlDp * fttsc
          DDttss( ll1, pp1, ii, jj ) = DDttss( ll1, pp1, ii, jj ) + DlDp * fttss

          DDtzcs( ll1, pp1, ii, jj ) = DDtzcs( ll1, pp1, ii, jj ) + DlDp * ftzcs
          DDtzsc( ll1, pp1, ii, jj ) = DDtzsc( ll1, pp1, ii, jj ) + DlDp * ftzsc
          DDtzss( ll1, pp1, ii, jj ) = DDtzss( ll1, pp1, ii, jj ) + DlDp * ftzss

          DDzzcs( ll1, pp1, ii, jj ) = DDzzcs( ll1, pp1, ii, jj ) + DlDp * fzzcs
          DDzzsc( ll1, pp1, ii, jj ) = DDzzsc( ll1, pp1, ii, jj ) + DlDp * fzzsc
          DDzzss( ll1, pp1, ii, jj ) = DDzzss( ll1, pp1, ii, jj ) + DlDp * fzzss
        end if !NOTstellsym
       
      enddo ! end of do lp2; 08 Feb 16;
     
    enddo ! end of do mn2; 08 Feb 16;
!!$OMP END DO
!!$OMP END PARALLEL    
  enddo ! end of do jquad; ! 16 Jan 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  DALLOCATE(basis)
  
  nele = SIZE(TTssss)
 
  call DSCAL(nele, pi2pi2nfphalf, DToocc, 1)
  call DSCAL(nele, pi2pi2nfphalf, TTssss, 1)
  call DSCAL(nele, pi2pi2nfphalf, TDstsc, 1)
  call DSCAL(nele, pi2pi2nfphalf, TDszsc, 1)
  call DSCAL(nele, pi2pi2nfphalf, DDttcc, 1)
  call DSCAL(nele, pi2pi2nfphalf, DDtzcc, 1)
  call DSCAL(nele, pi2pi2nfphalf, DDzzcc, 1)

  if (NOTstellsym) then

    call DSCAL(nele, pi2pi2nfphalf, DToocs, 1)
    call DSCAL(nele, pi2pi2nfphalf, DToosc, 1)
    call DSCAL(nele, pi2pi2nfphalf, DTooss, 1)

    call DSCAL(nele, pi2pi2nfphalf, TTsscc, 1)
    call DSCAL(nele, pi2pi2nfphalf, TTsscs, 1)
    call DSCAL(nele, pi2pi2nfphalf, TTsssc, 1)

    call DSCAL(nele, pi2pi2nfphalf, TDstcc, 1)
    call DSCAL(nele, pi2pi2nfphalf, TDstcs, 1)
    call DSCAL(nele, pi2pi2nfphalf, TDstss, 1)

    call DSCAL(nele, pi2pi2nfphalf, TDszcc, 1)
    call DSCAL(nele, pi2pi2nfphalf, TDszcs, 1)
    call DSCAL(nele, pi2pi2nfphalf, TDszss, 1)

    call DSCAL(nele, pi2pi2nfphalf, DDttsc, 1)
    call DSCAL(nele, pi2pi2nfphalf, DDttcs, 1)
    call DSCAL(nele, pi2pi2nfphalf, DDttss, 1)

    call DSCAL(nele, pi2pi2nfphalf, DDtzsc, 1)
    call DSCAL(nele, pi2pi2nfphalf, DDtzcs, 1)
    call DSCAL(nele, pi2pi2nfphalf, DDtzss, 1)

    call DSCAL(nele, pi2pi2nfphalf, DDzzsc, 1)
    call DSCAL(nele, pi2pi2nfphalf, DDzzcs, 1)
    call DSCAL(nele, pi2pi2nfphalf, DDzzss, 1)

  end if !NOTstellsym
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN( ma00aa )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine ma00aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
