!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (integrals) ! Calculates volume integrals of Chebyshev-polynomials and metric elements for preconditioner.


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine spsint( lquad, mn, lvol, lrad )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half, one, two, pi, pi2
  
  use numerical, only : vsmall, small, sqrtmachprec
  
  use fileunits, only : ounit
  
  use inputlist, only : mpol, Wspsint
  
  use cputiming, only : Tspsint
  
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
                        cheby, zernike, &
                        Lcoordinatesingularity, regumm, &
                        pi2pi2nfp, pi2pi2nfphalf
                        
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in) :: lquad, mn, lvol, lrad
  
  INTEGER             :: jquad, ll, pp, ll1, pp1, uv, ii, jj, io, mn2, lp2, mn2_max, lp2_max, nele, mi
  
  INTEGER             :: kk, kd, kka, kks, kda, kds
  
  REAL                :: lss, jthweight, fee, feo, foe, foo, Tl, Dl, Tp, Dp, TlTp, TlDp, DlTp, DlDp, ikda, ikds, imn2, ilrad, lssm

  REAL                :: foocc, foocs, foosc, fooss
  REAL                :: fsscc, fsscs, fsssc, fssss
  REAL                :: fstcc, fstcs, fstsc, fstss
  REAL                :: fszcc, fszcs, fszsc, fszss
  REAL                :: fttcc, fttcs, fttsc, fttss
  REAL                :: ftzcc, ftzcs, ftzsc, ftzss
  REAL                :: fzzcc, fzzcs, fzzsc, fzzss
  
  REAL                :: sbar(1:lquad), halfoversbar(1:lquad), sbarhim(1:lquad,1:mn) ! regularization factors; 10 Dec 15;
  
  BEGIN( spsint )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  mn2_max = mn*mn
  lp2_max = (lrad+1)*(lrad+1)
  imn2    =  one/real(mn)
  ilrad = one/real(lrad+1)

  DToocc = zero
  TTssss = zero
  !TDstsc = zero
  !TDszsc = zero
  DDttcc = zero
  DDtzcc = zero
  DDzzcc = zero

  if (NOTstellsym) then
    !DToocs = zero
    !DToosc = zero
    DTooss = zero

    TTsscc = zero
    !TTsscs = zero
    !TTsssc = zero

    TDstcc = zero
    !TDstcs = zero
    TDstss = zero

    TDszcc = zero
    !TDszcs = zero
    TDszss = zero

    !DDttcs = zero
    !DDttsc = zero
    DDttss = zero

    !DDtzcs = zero
    !DDtzsc = zero
    DDtzss = zero

    !DDzzcs = zero
    !DDzzsc = zero
    DDzzss = zero
  endif !NOTstellsym

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Lcoordinatesingularity ) then
   ! switch to sbar=r; 29 Jun 19
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   sbar(1:lquad) = ( gaussianabscissae(1:lquad,lvol) + one ) * half
   
   halfoversbar(1:lquad) = half / sbar(1:lquad)
  
   do jquad = 1, lquad ! Gaussian quadrature loop;
    
    lss = gaussianabscissae(jquad,lvol) ; jthweight = gaussianweight(jquad,lvol)

    call get_zernike(sbar(jquad), lrad, mpol, zernike(:,:,0:1)) ! use Zernike polynomials 29 Jun 19;

    WCALL( spsint, spsmetrix,( lvol, lss, jquad ) ) ! compute metric elements; 16 Jan 13;
 
    foocc = goomne(0) * jthweight
    fssss = gssmne(0) * jthweight
    !fstsc = zero * jthweight
    !fszsc = zero * jthweight
    fttcc = gttmne(0) * jthweight
    ftzcc = gtzmne(0) * jthweight
    fzzcc = gzzmne(0) * jthweight


    if (NOTstellsym) then
      !foocs = zero * jthweight
      !foosc = zero * jthweight
      fooss = goomne(0) * jthweight

      fsscc = gssmne(0) * jthweight
      !fsscs = zero * jthweight
      !fsssc = zero * jthweight

      fstcc = gstmne(0) * jthweight
      !fstcs = zero * jthweight
      fstss = gstmne(0) * jthweight

      fszcc = gszmne(0) * jthweight
      !fszcs = zero * jthweight
      fszss = gszmne(0) * jthweight

      !fttcs = zero * jthweight
      !fttsc = zero * jthweight
      fttss = gttmne(0) * jthweight

      !ftzcs = zero * jthweight
      !ftzsc = zero * jthweight
      ftzss = gtzmne(0) * jthweight

      !fzzcs = zero * jthweight
      !fzzsc = zero * jthweight
      fzzss = gzzmne(0) * jthweight
    end if !NOTstellsym

    do mi = 0, mpol

      do lp2 = 1, lp2_max 
        ll = mod(lp2-1,lrad+1)
        pp = (lp2-ll-1)/(lrad+1)

        ll1 = (ll - mod(ll,2))/2 ! shrinked dof for Zernike; 02 Jul 19
        pp1 = (pp - mod(pp,2))/2 ! shrinked dof for Zernike; 02 Jul 19

        if (ll < mi) cycle ! zernike only non-zero for ll>=ii
        if (pp < mi) cycle ! zernike only non-zero for pp>=jj
        if (mod(ll+mi,2)/=0) cycle ! zernike only non-zero if ll and ii have the same parity
        if (mod(pp+mi,2)/=0) cycle ! zernike only non-zero if pp and jj have the same parity

        Tl = zernike(ll, mi, 0)         ! use Zernike polynomials 29 Jun 19;
        Dl = zernike(ll, mi, 1) * half  ! use Zernike polynomials 29 Jun 19;

        Tp = zernike(pp, mi, 0)         ! use Zernike polynomials 29 Jun 19;
        Dp = zernike(pp, mi, 1) * half  ! use Zernike polynomials 29 Jun 19;

        TlTp = Tl * Tp
        TlDp = Tl * Dp
        DlTp = Dl * Tp
        DlDp = Dl * Dp 

        DToocc( ll1, pp1, mi, 1 ) = DToocc( ll1, pp1, mi, 1 ) + DlTp * foocc
        TTssss( ll1, pp1, mi, 1 ) = TTssss( ll1, pp1, mi, 1 ) + TlTp * fssss
        !TDstsc( ll1, pp1, mi, 1 ) = TDstsc( ll1, pp1, mi, 1 ) + TlDp * fstsc
        !TDszsc( ll1, pp1, mi, 1 ) = TDszsc( ll1, pp1, mi, 1 ) + TlDp * fszsc
        DDttcc( ll1, pp1, mi, 1 ) = DDttcc( ll1, pp1, mi, 1 ) + DlDp * fttcc
        DDtzcc( ll1, pp1, mi, 1 ) = DDtzcc( ll1, pp1, mi, 1 ) + DlDp * ftzcc
        DDzzcc( ll1, pp1, mi, 1 ) = DDzzcc( ll1, pp1, mi, 1 ) + DlDp * fzzcc

        if (NOTstellsym) then
          !DToocs( ll1, pp1, mi, 1 ) = DToocs( ll1, pp1, mi, 1 ) + DlTp * foocs
          !DToosc( ll1, pp1, mi, 1 ) = DToosc( ll1, pp1, mi, 1 ) + DlTp * foosc
          DTooss( ll1, pp1, mi, 1 ) = DTooss( ll1, pp1, mi, 1 ) + DlTp * fooss

          TTsscc( ll1, pp1, mi, 1 ) = TTsscc( ll1, pp1, mi, 1 ) + TlTp * fsscc
          !TTsscs( ll1, pp1, mi, 1 ) = TTsscs( ll1, pp1, mi, 1 ) + TlTp * fsscs
          !TTsssc( ll1, pp1, mi, 1 ) = TTsssc( ll1, pp1, mi, 1 ) + TlTp * fsssc

          TDstcc( ll1, pp1, mi, 1 ) = TDstcc( ll1, pp1, mi, 1 ) + TlDp * fstcc
          !TDstcs( ll1, pp1, mi, 1 ) = TDstcs( ll1, pp1, mi, 1 ) + TlDp * fstcs
          !TDstss( ll1, pp1, mi, 1 ) = TDstss( ll1, pp1, mi, 1 ) + TlDp * fstss

          TDszcc( ll1, pp1, mi, 1 ) = TDszcc( ll1, pp1, mi, 1 ) + TlDp * fszcc
          !TDszcs( ll1, pp1, mi, 1 ) = TDszcs( ll1, pp1, mi, 1 ) + TlDp * fszcs
          TDszss( ll1, pp1, mi, 1 ) = TDszss( ll1, pp1, mi, 1 ) + TlDp * fszss

          DDttcs( ll1, pp1, mi, 1 ) = DDttcs( ll1, pp1, mi, 1 ) + DlDp * fttcs
          !DDttsc( ll1, pp1, mi, 1 ) = DDttsc( ll1, pp1, mi, 1 ) + DlDp * fttsc
          DDttss( ll1, pp1, mi, 1 ) = DDttss( ll1, pp1, mi, 1 ) + DlDp * fttss

          !DDtzcs( ll1, pp1, mi, 1 ) = DDtzcs( ll1, pp1, mi, 1 ) + DlDp * ftzcs
          !DDtzsc( ll1, pp1, mi, 1 ) = DDtzsc( ll1, pp1, mi, 1 ) + DlDp * ftzsc
          DDtzss( ll1, pp1, mi, 1 ) = DDtzss( ll1, pp1, mi, 1 ) + DlDp * ftzss

          !DDzzcs( ll1, pp1, mi, 1 ) = DDzzcs( ll1, pp1, mi, 1 ) + DlDp * fzzcs
          !DDzzsc( ll1, pp1, mi, 1 ) = DDzzsc( ll1, pp1, mi, 1 ) + DlDp * fzzsc
          DDzzss( ll1, pp1, mi, 1 ) = DDzzss( ll1, pp1, mi, 1 ) + DlDp * fzzss
        end if !NOTstellsym
       
      enddo ! end of do lp2; 08 Feb 16;
     
    enddo ! end of do mi
 
   enddo ! end of do jquad; ! 16 Jan 13;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
  else ! .not.Lcoordinatesingularity; 17 Dec 15;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   do jquad = 1, lquad ! Gaussian quadrature loop;

    lss = gaussianabscissae(jquad,lvol) ; jthweight = gaussianweight(jquad,lvol)

    call get_cheby(lss, lrad, cheby(0:lrad,0:1))
    
    WCALL( spsint, spsmetrix, ( lvol, lss, jquad ) ) ! compute metric elements; 16 Jan 13;
    
    mi = 0
   
    foocc = goomne(0)
    fssss = gssmne(0)
    !fstsc = zero
    !fszsc = zero
    fttcc = gttmne(0)
    ftzcc = gtzmne(0)
    fzzcc = gtzmne(0)

    if (NOTstellsym) then
      !foocs = zero
      !foosc = zero
      fooss = goomne(0)

      fsscc = gssmne(0)
      !fsscs = zero
      !fsssc = zero

      fstcc = gstmne(0)
      !fstcs = zero 
      fstss = gstmne(0)

      fszcc = gszmne(0)
      !fszcs = zero 
      fszss = gszmne(0)

      !fttcs = zero
      !fttsc = zero
      fttss = gttmne(0)

      !ftzcs = zero
      !ftzsc = zero
      ftzss = gtzmne(0)

      !fzzcs = zero
      !fzzsc = zero
      fzzss = gzzmne(0)
    end if !NOTstellsym

    do lp2 = 1, lp2_max 
      ll = mod(lp2-1,lrad+1)
      pp = (lp2-ll-1)/(lrad+1)
      
      Tl = cheby(ll,0)
      Dl = cheby(ll,1)
      Tp = cheby(pp,0)
      Dp = cheby(pp,1)
      
      TlTp = Tl * Tp
      TlDp = Tl * Dp
      DlTp = Dl * Tp
      DlDp = Dl * Dp

      DToocc( ll1, pp1, mi, 1 ) = DToocc( ll1, pp1, mi, 1 ) + DlTp * foocc
      TTssss( ll1, pp1, mi, 1 ) = TTssss( ll1, pp1, mi, 1 ) + TlTp * fssss
      !TDstsc( ll1, pp1, mi, 1 ) = TDstsc( ll1, pp1, mi, 1 ) + TlDp * fstsc
      !TDszsc( ll1, pp1, mi, 1 ) = TDszsc( ll1, pp1, mi, 1 ) + TlDp * fszsc
      DDttcc( ll1, pp1, mi, 1 ) = DDttcc( ll1, pp1, mi, 1 ) + DlDp * fttcc
      DDtzcc( ll1, pp1, mi, 1 ) = DDtzcc( ll1, pp1, mi, 1 ) + DlDp * ftzcc
      DDzzcc( ll1, pp1, mi, 1 ) = DDzzcc( ll1, pp1, mi, 1 ) + DlDp * fzzcc

      if (NOTstellsym) then
        !DToocs( ll1, pp1, mi, 1 ) = DToocs( ll1, pp1, mi, 1 ) + DlTp * foocs
        !DToosc( ll1, pp1, mi, 1 ) = DToosc( ll1, pp1, mi, 1 ) + DlTp * foosc
        DTooss( ll1, pp1, mi, 1 ) = DTooss( ll1, pp1, mi, 1 ) + DlTp * fooss

        TTsscc( ll1, pp1, mi, 1 ) = TTsscc( ll1, pp1, mi, 1 ) + TlTp * fsscc
        !TTsscs( ll1, pp1, mi, 1 ) = TTsscs( ll1, pp1, mi, 1 ) + TlTp * fsscs
        !TTsssc( ll1, pp1, mi, 1 ) = TTsssc( ll1, pp1, mi, 1 ) + TlTp * fsssc

        TDstcc( ll1, pp1, mi, 1 ) = TDstcc( ll1, pp1, mi, 1 ) + TlDp * fstcc
        !TDstcs( ll1, pp1, mi, 1 ) = TDstcs( ll1, pp1, mi, 1 ) + TlDp * fstcs
        !TDstss( ll1, pp1, mi, 1 ) = TDstss( ll1, pp1, mi, 1 ) + TlDp * fstss

        TDszcc( ll1, pp1, mi, 1 ) = TDszcc( ll1, pp1, mi, 1 ) + TlDp * fszcc
        !TDszcs( ll1, pp1, mi, 1 ) = TDszcs( ll1, pp1, mi, 1 ) + TlDp * fszcs
        TDszss( ll1, pp1, mi, 1 ) = TDszss( ll1, pp1, mi, 1 ) + TlDp * fszss

        DDttcs( ll1, pp1, mi, 1 ) = DDttcs( ll1, pp1, mi, 1 ) + DlDp * fttcs
        !DDttsc( ll1, pp1, mi, 1 ) = DDttsc( ll1, pp1, mi, 1 ) + DlDp * fttsc
        DDttss( ll1, pp1, mi, 1 ) = DDttss( ll1, pp1, mi, 1 ) + DlDp * fttss

        !DDtzcs( ll1, pp1, mi, 1 ) = DDtzcs( ll1, pp1, mi, 1 ) + DlDp * ftzcs
        !DDtzsc( ll1, pp1, mi, 1 ) = DDtzsc( ll1, pp1, mi, 1 ) + DlDp * ftzsc
        DDtzss( ll1, pp1, mi, 1 ) = DDtzss( ll1, pp1, mi, 1 ) + DlDp * ftzss

        !DDzzcs( ll1, pp1, mi, 1 ) = DDzzcs( ll1, pp1, mi, 1 ) + DlDp * fzzcs
        !DDzzsc( ll1, pp1, mi, 1 ) = DDzzsc( ll1, pp1, mi, 1 ) + DlDp * fzzsc
        DDzzss( ll1, pp1, mi, 1 ) = DDzzss( ll1, pp1, mi, 1 ) + DlDp * fzzss
      end if !NOTstellsym
       
     enddo ! end of do lp2;  1 Feb 13;

    enddo ! end of do jquad; ! 16 Jan 13;


  endif ! end of if( Lcoordinatesingularity ) ; 17 Dec 15;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  nele = SIZE(TTssss)
 
  call DSCAL(nele, pi2pi2nfphalf, DToocc, 1)
  call DSCAL(nele, pi2pi2nfphalf, TTssss, 1)
  !call DSCAL(nele, pi2pi2nfphalf, TDstsc, 1)
  !call DSCAL(nele, pi2pi2nfphalf, TDszsc, 1)
  call DSCAL(nele, pi2pi2nfphalf, DDttcc, 1)
  call DSCAL(nele, pi2pi2nfphalf, DDtzcc, 1)
  call DSCAL(nele, pi2pi2nfphalf, DDzzcc, 1)

  if (NOTstellsym) then

    !call DSCAL(nele, pi2pi2nfphalf, DToocs, 1)
    !call DSCAL(nele, pi2pi2nfphalf, DToosc, 1)
    call DSCAL(nele, pi2pi2nfphalf, DTooss, 1)

    call DSCAL(nele, pi2pi2nfphalf, TTsscc, 1)
    !call DSCAL(nele, pi2pi2nfphalf, TTsscs, 1)
    !call DSCAL(nele, pi2pi2nfphalf, TTsssc, 1)

    call DSCAL(nele, pi2pi2nfphalf, TDstcc, 1)
    !call DSCAL(nele, pi2pi2nfphalf, TDstcs, 1)
    call DSCAL(nele, pi2pi2nfphalf, TDstss, 1)

    call DSCAL(nele, pi2pi2nfphalf, TDszcc, 1)
    !call DSCAL(nele, pi2pi2nfphalf, TDszcs, 1)
    call DSCAL(nele, pi2pi2nfphalf, TDszss, 1)

    !call DSCAL(nele, pi2pi2nfphalf, DDttsc, 1)
    !call DSCAL(nele, pi2pi2nfphalf, DDttcs, 1)
    call DSCAL(nele, pi2pi2nfphalf, DDttss, 1)

    !call DSCAL(nele, pi2pi2nfphalf, DDtzsc, 1)
    !call DSCAL(nele, pi2pi2nfphalf, DDtzcs, 1)
    call DSCAL(nele, pi2pi2nfphalf, DDtzss, 1)

    !call DSCAL(nele, pi2pi2nfphalf, DDzzsc, 1)
    !call DSCAL(nele, pi2pi2nfphalf, DDzzcs, 1)
    call DSCAL(nele, pi2pi2nfphalf, DDzzss, 1)

  end if !NOTstellsym

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN( spsint )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine spsint

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
