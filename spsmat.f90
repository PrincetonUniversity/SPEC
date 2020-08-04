!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (build matrices) ! Constructs matrices for the precondtioner.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine spsmat( lvol, mn, lrad )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, one, two
  
  use numerical, only : small
  
  use fileunits, only : ounit
  
  use inputlist, only : Wmacros, Wspsmat, mpol
  
  use cputiming, only : Tspsmat
  
  use allglobal, only : ncpu, myid, cpus, &
                        YESstellsym, NOTstellsym, &
                        im, in, &
                        NAdof, &
                        dMA, dMD, dMB, dMG, &
                        LILUprecond, NdMASmax, NdMAS, dMAS, dMDS, idMAS, jdMAS, & ! preconditioning matrix
                        Ate, Ato, Aze, Azo, &
                        iVns, iBns, iVnc, iBnc, &
                        Lma, Lmb, Lmc, Lmd, Lme, Lmf, Lmg, Lmh, &
                        Lcoordinatesingularity, TT, RTT, RTM, &
                        DToocc, DToocs, DToosc, DTooss, &
                        TTsscc, TTsscs, TTsssc, TTssss, &
                        TDstcc, TDstcs, TDstsc, TDstss, &
                        TDszcc, TDszcs, TDszsc, TDszss, &
                        DDttcc, DDttcs, DDttsc, DDttss, &
                        DDtzcc, DDtzcs, DDtzsc, DDtzss, &
                        DDzzcc, DDzzcs, DDzzsc, DDzzss
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  INTEGER, intent(in)  :: lvol, mn, lrad
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  INTEGER              :: NN, ii, jj, ll, kk, pp, ll1, pp1, mi, ni, mj, nj, mimj, minj, nimj, ninj, mjmi, mjni, njmi, njni, id, jd, idx
  
  REAL                 :: Wtete, Wteto, Wtote, Wtoto
  REAL                 :: Wteze, Wtezo, Wtoze, Wtozo
  REAL                 :: Wzete, Wzeto, Wzote, Wzoto
  REAL                 :: Wzeze, Wzezo, Wzoze, Wzozo
  
  REAL                 :: Htete, Hteto, Htote, Htoto
  REAL                 :: Hteze, Htezo, Htoze, Htozo
  REAL                 :: Hzete, Hzeto, Hzote, Hzoto
  REAL                 :: Hzeze, Hzezo, Hzoze, Hzozo
  REAL                 :: adata, ddata, factorcc, factorss
  
  REAL,allocatable     :: dMASqueue(:,:), dMDSqueue(:,:), TTdata(:,:,:), TTMdata(:,:) ! queues to construct sparse matrices
  INTEGER,allocatable  :: jdMASqueue(:,:) ! indices
  INTEGER              :: nqueue(4), nrow, ns, nmaxqueue

  BEGIN(spsmat)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  NN = NAdof(lvol) ! shorthand;

  NdMAS(lvol) = 0
  nrow = 0
  dMAS = zero
  dMDS = zero
  idMAS = 0
  ns = 0

  nmaxqueue = 4 * lrad + 10 ! estimate the size of the queue
  
  SALLOCATE( dMASqueue, (1:nmaxqueue, 4), zero)
  SALLOCATE( dMDSqueue, (1:nmaxqueue, 4), zero)
  SALLOCATE( jdMASqueue, (1:nmaxqueue, 4), zero)

  SALLOCATE( TTdata, (0:lrad, 0:mpol, 0:1), zero)
  SALLOCATE( TTMdata, (0:lrad, 0:mpol), zero)

  ! fill in Zernike/Chebyshev polynomials depending on Lcooridnatesingularity
  if (Lcoordinatesingularity) then
    TTdata(0:lrad,0:mpol,0:1) = RTT(0:lrad,0:mpol,0:1,0)
    TTMdata(0:lrad,0:mpol) = RTM(0:lrad,0:mpol)
  else
    do ii = 0, mpol
      TTdata(0:lrad,ii,0:1) = TT(0:lrad,0:1,0)
      TTMdata(0:lrad,ii) = TT(0:lrad,0,0)
    enddo
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( YESstellsym ) then
 
    do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
    
      jj = ii ; mj = im(jj) ; nj = in(jj) ; mimj = mi * mj ; minj = mi * nj ; nimj = ni * mj ; ninj = ni * nj

      if (Lcoordinatesingularity) then
        idx = mi + 1
      else
        idx = 1
      endif

      if (ii.eq.1) then
        factorcc = two
        factorss = zero
      else
        factorcc = one
        factorss = one
      endif
      
      do ll = 0, lrad
        ! clean up the queue
        call clean_queue(nqueue, nmaxqueue, dMASqueue, dMDSqueue, jdMASqueue)
        
        if (Lcoordinatesingularity) then
          if (ll.lt.mi .or. mod(ll+mi,2).gt.0) cycle
        endif

        do pp = 0, lrad
        
        if (Lcoordinatesingularity) then
          if (pp.lt.mj .or. mod(pp+mj,2).gt.0) cycle ! rule out zero components of Zernike; 02 Jul 19
          ll1 = (ll - mod(ll, 2)) / 2 ! shrinked dof for Zernike; 02 Jul 19
          pp1 = (pp - mod(pp, 2)) / 2 ! shrinked dof for Zernike; 02 Jul 19
        else
          ll1 = ll
          pp1 = pp
        endif

        Wtete = + 2 * ninj * TTssss(ll1,pp1,idx,1) * factorss + 2 * DDzzcc(ll1,pp1,idx,1) * factorcc
        Wzete = + 2 * nimj * TTssss(ll1,pp1,idx,1) * factorss - 2 * DDtzcc(pp1,ll1,idx,1) * factorcc
        Wteze = + 2 * minj * TTssss(ll1,pp1,idx,1) * factorss - 2 * DDtzcc(ll1,pp1,idx,1) * factorcc
        Wzeze = + 2 * mimj * TTssss(ll1,pp1,idx,1) * factorss + 2 * DDttcc(ll1,pp1,idx,1) * factorcc
        
        Htete =   zero
        Hzete = (- DToocc(pp1,ll1,idx,1) + DToocc(ll1,pp1,idx,1)) * factorcc
        Hteze = (+ DToocc(pp1,ll1,idx,1) - DToocc(ll1,pp1,idx,1)) * factorcc
        Hzeze =   zero
        
        id = Ate(lvol,0,ii)%i(ll) ; jd = Ate(lvol,0,jj)%i(pp)
        if (id.ne.0 .and. jd.ne.0) call push_back(1,nqueue,nmaxqueue,Wtete,Htete,jd,dMASqueue,dMDSqueue,jdMASqueue)
        ;                         ; jd = Aze(lvol,0,jj)%i(pp)
        if (id.ne.0 .and. jd.ne.0) call push_back(1,nqueue,nmaxqueue,Wzete,Hzete,jd,dMASqueue,dMDSqueue,jdMASqueue)
        id = Aze(lvol,0,ii)%i(ll) ; jd = Ate(lvol,0,jj)%i(pp)
        if (id.ne.0 .and. jd.ne.0) call push_back(2,nqueue,nmaxqueue,Wteze,Hteze,jd,dMASqueue,dMDSqueue,jdMASqueue)
        ;                         ; jd = Aze(lvol,0,jj)%i(pp)
        if (id.ne.0 .and. jd.ne.0) call push_back(2,nqueue,nmaxqueue,Wzeze,Hzeze,jd,dMASqueue,dMDSqueue,jdMASqueue)
        
        enddo ! end of do pp ;
        
        if( Lcoordinatesingularity .and. ii.eq.1 ) then ; kk = 1
        else                                            ; kk = 0
        endif

        ddata = zero
        ;                  ; id = Ate(lvol,0,ii)%i(ll) ; jd = Lma(lvol,  ii) ; adata = +      TTMdata(ll, mi)
        if (id.ne.0 .and. jd.ne.0)   call push_back(1,nqueue,nmaxqueue,adata,ddata,jd,dMASqueue,dMDSqueue,jdMASqueue)
        ;                  ; id = Aze(lvol,0,ii)%i(ll) ; jd = Lmb(lvol,  ii) ; adata = +      TTdata(ll, mi,kk)
        if (id.ne.0 .and. jd.ne.0)   call push_back(2,nqueue,nmaxqueue,adata,ddata,jd,dMASqueue,dMDSqueue,jdMASqueue)
        if( ii.gt.1 ) then ; id = Ate(lvol,0,ii)%i(ll) ; jd = Lme(lvol,  ii) ; adata = - ni * TTdata(ll, mi, 1)
          if (id.ne.0 .and. jd.ne.0) call push_back(1,nqueue,nmaxqueue,adata,ddata,jd,dMASqueue,dMDSqueue,jdMASqueue)
        ;                  ; id = Aze(lvol,0,ii)%i(ll) ; jd = Lme(lvol,  ii) ; adata = - mi * TTdata(ll, mi, 1)
          if (id.ne.0 .and. jd.ne.0) call push_back(2,nqueue,nmaxqueue,adata,ddata,jd,dMASqueue,dMDSqueue,jdMASqueue)
        else               ; id = Ate(lvol,0,ii)%i(ll) ; jd = Lmg(lvol,  ii) ; adata = +      TTdata(ll, mi, 1)
          if (id.ne.0 .and. jd.ne.0) call push_back(1,nqueue,nmaxqueue,adata,ddata,jd,dMASqueue,dMDSqueue,jdMASqueue)
        ;                  ; id = Aze(lvol,0,ii)%i(ll) ; jd = Lmh(lvol,  ii) ; adata = +      TTdata(ll, mi, 1)
          if (id.ne.0 .and. jd.ne.0) call push_back(2,nqueue,nmaxqueue,adata,ddata,jd,dMASqueue,dMDSqueue,jdMASqueue)
        endif
        
        ! putting things in the sparse matrix
        call addline(nqueue, nmaxqueue, dMASqueue, dMDSqueue, jdMASqueue, ns, nrow, dMAS, dMDS, jdMAS, idMAS)
        
      enddo ! end of do ll ;
 
    enddo ! end of do ii ;

  else ! NOTstellsym

    do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
    
      jj = ii ; mj = im(jj) ; nj = in(jj) ; mimj = mi * mj ; minj = mi * nj ; nimj = ni * mj ; ninj = ni * nj

      if (Lcoordinatesingularity) then
        idx = mi + 1
      else
        idx = 1
      endif

      if (ii.eq.1) then
        factorcc = two
        factorss = zero
      else
        factorcc = one
        factorss = one
      endif
      
      do ll = 0, lrad
        ! clean up the queue
        call clean_queue(nqueue, nmaxqueue, dMASqueue, dMDSqueue, jdMASqueue)
        
        if (Lcoordinatesingularity) then
          if (ll.lt.mi .or. mod(ll+mi,2).gt.0) cycle
        endif

        do pp = 0, lrad
        
        if (Lcoordinatesingularity) then
          if (pp.lt.mj .or. mod(pp+mj,2).gt.0) cycle ! rule out zero components of Zernike; 02 Jul 19
          ll1 = (ll - mod(ll, 2)) / 2 ! shrinked dof for Zernike; 02 Jul 19
          pp1 = (pp - mod(pp, 2)) / 2 ! shrinked dof for Zernike; 02 Jul 19
        else
          ll1 = ll
          pp1 = pp
        endif

        Wtete = + 2 * ninj * TTssss(ll1,pp1,idx,1) * factorss + 2 * DDzzcc(ll1,pp1,idx,1) * factorcc
        Wzete = + 2 * nimj * TTssss(ll1,pp1,idx,1) * factorss - 2 * DDtzcc(pp1,ll1,idx,1) * factorcc
        Wteze = + 2 * minj * TTssss(ll1,pp1,idx,1) * factorss - 2 * DDtzcc(ll1,pp1,idx,1) * factorcc
        Wzeze = + 2 * mimj * TTssss(ll1,pp1,idx,1) * factorss + 2 * DDttcc(ll1,pp1,idx,1) * factorcc
              
        Wtote = 2 * (   + nj * TDszcc(pp1,ll1,idx,1) * factorcc - ni * TDszss(ll1,pp1,idx,1) * factorss )
        Wzote = 2 * (   + mj * TDszcc(pp1,ll1,idx,1) * factorcc + ni * TDstss(ll1,pp1,idx,1) * factorss )
       
        Wteto = 2 * (   - nj * TDszss(pp1,ll1,idx,1) * factorss + ni * TDszcc(ll1,pp1,idx,1) * factorcc )
        Wtoto = 2 * ( + njni * TTsscc(pp1,ll1,idx,1) * factorcc +      DDzzss(pp1,ll1,idx,1) * factorss )
        Wzeto = 2 * (   - mj * TDszss(pp1,ll1,idx,1) * factorss - ni * TDstcc(ll1,pp1,idx,1) * factorcc )
        Wzoto = 2 * ( + mjni * TTsscc(pp1,ll1,idx,1) * factorcc -      DDtzss(pp1,ll1,idx,1) * factorss )
       
        Wtoze = 2 * (   - nj * TDstcc(pp1,ll1,idx,1) * factorcc - mi * TDszss(ll1,pp1,idx,1) * factorss )
        Wzoze = 2 * (   - mj * TDstcc(pp1,ll1,idx,1) * factorcc + mi * TDstss(ll1,pp1,idx,1) * factorss )
       
        Wtezo = 2 * (   + nj * TDstss(pp1,ll1,idx,1) * factorss + mi * TDszcc(ll1,pp1,idx,1) * factorcc )
        Wtozo = 2 * ( + njmi * TTsscc(pp1,ll1,idx,1) * factorcc -      DDtzss(pp1,ll1,idx,1) * factorss )
        Wzezo = 2 * (   + mj * TDstss(pp1,ll1,idx,1) * factorss - mi * TDstcc(ll1,pp1,idx,1) * factorcc )
        Wzozo = 2 * ( + mjmi * TTsscc(pp1,ll1,idx,1) * factorcc +      DDttss(pp1,ll1,idx,1) * factorss )
        
        Htete =   zero
        Hzete = (- DToocc(pp1,ll1,idx,1) + DToocc(ll1,pp1,idx,1)) * factorcc
        Hteze = (+ DToocc(pp1,ll1,idx,1) - DToocc(ll1,pp1,idx,1)) * factorcc
        Hzeze =   zero


        Htote =   zero
        Hzote =   zero
       
        Hteto =   zero
        Htoto =   zero
        Hzeto =   zero
        Hzoto = (- DTooss(pp1,ll1,idx,1) + DTooss(ll1,pp1,idx,1)) * factorss
       
        Htoze =   zero
        Hzeze =   zero 
        Hzoze =   zero
       
        Htezo =   zero
        Htozo = (+ DTooss(pp1,ll1,idx,1) - DTooss(ll1,pp1,idx,1)) * factorss
        Hzezo =   zero
        Hzozo =   zero

        
        id = Ate(lvol,0,ii)%i(ll) ; jd = Ate(lvol,0,jj)%i(pp)
        if (id.ne.0 .and. jd.ne.0) call push_back(1,nqueue,nmaxqueue,Wtete,Htete,jd,dMASqueue,dMDSqueue,jdMASqueue)
        ;                         ; jd = Aze(lvol,0,jj)%i(pp)
        if (id.ne.0 .and. jd.ne.0) call push_back(1,nqueue,nmaxqueue,Wzete,Hzete,jd,dMASqueue,dMDSqueue,jdMASqueue)
        ;                         ; jd = Ato(lvol,0,jj)%i(pp)
        if (id.ne.0 .and. jd.ne.0) call push_back(1,nqueue,nmaxqueue,Wtote,Htote,jd,dMASqueue,dMDSqueue,jdMASqueue)
        ;                         ; jd = Azo(lvol,0,jj)%i(pp)
        if (id.ne.0 .and. jd.ne.0) call push_back(1,nqueue,nmaxqueue,Wzote,Hzote,jd,dMASqueue,dMDSqueue,jdMASqueue)

        id = Aze(lvol,0,ii)%i(ll) ; jd = Ate(lvol,0,jj)%i(pp)
        if (id.ne.0 .and. jd.ne.0) call push_back(2,nqueue,nmaxqueue,Wteze,Hteze,jd,dMASqueue,dMDSqueue,jdMASqueue)
        ;                         ; jd = Aze(lvol,0,jj)%i(pp)
        if (id.ne.0 .and. jd.ne.0) call push_back(2,nqueue,nmaxqueue,Wzeze,Hzeze,jd,dMASqueue,dMDSqueue,jdMASqueue)
        ;                         ; jd = Ato(lvol,0,jj)%i(pp)
        if (id.ne.0 .and. jd.ne.0) call push_back(2,nqueue,nmaxqueue,Wtoze,Htoze,jd,dMASqueue,dMDSqueue,jdMASqueue)
        ;                         ; jd = Azo(lvol,0,jj)%i(pp)
        if (id.ne.0 .and. jd.ne.0) call push_back(2,nqueue,nmaxqueue,Wzoze,Hzoze,jd,dMASqueue,dMDSqueue,jdMASqueue)

        id = Ato(lvol,0,ii)%i(ll) ; jd = Ate(lvol,0,jj)%i(pp)
        if (id.ne.0 .and. jd.ne.0) call push_back(3,nqueue,nmaxqueue,Wteto,Hteto,jd,dMASqueue,dMDSqueue,jdMASqueue)
        ;                         ; jd = Aze(lvol,0,jj)%i(pp)
        if (id.ne.0 .and. jd.ne.0) call push_back(3,nqueue,nmaxqueue,Wzeto,Hzeto,jd,dMASqueue,dMDSqueue,jdMASqueue)
        ;                         ; jd = Ato(lvol,0,jj)%i(pp)
        if (id.ne.0 .and. jd.ne.0) call push_back(3,nqueue,nmaxqueue,Wtoto,Htoto,jd,dMASqueue,dMDSqueue,jdMASqueue)
        ;                         ; jd = Azo(lvol,0,jj)%i(pp)
        if (id.ne.0 .and. jd.ne.0) call push_back(3,nqueue,nmaxqueue,Wzoto,Hzoto,jd,dMASqueue,dMDSqueue,jdMASqueue)

        id = Azo(lvol,0,ii)%i(ll) ; jd = Ate(lvol,0,jj)%i(pp)
        if (id.ne.0 .and. jd.ne.0) call push_back(4,nqueue,nmaxqueue,Wtezo,Htezo,jd,dMASqueue,dMDSqueue,jdMASqueue)
        ;                         ; jd = Aze(lvol,0,jj)%i(pp)
        if (id.ne.0 .and. jd.ne.0) call push_back(4,nqueue,nmaxqueue,Wzezo,Hzezo,jd,dMASqueue,dMDSqueue,jdMASqueue)
        ;                         ; jd = Ato(lvol,0,jj)%i(pp)
        if (id.ne.0 .and. jd.ne.0) call push_back(4,nqueue,nmaxqueue,Wtozo,Htozo,jd,dMASqueue,dMDSqueue,jdMASqueue)
        ;                         ; jd = Azo(lvol,0,jj)%i(pp)
        if (id.ne.0 .and. jd.ne.0) call push_back(4,nqueue,nmaxqueue,Wzozo,Hzozo,jd,dMASqueue,dMDSqueue,jdMASqueue)
        
        enddo ! end of do pp ;
        
        if( Lcoordinatesingularity .and. ii.eq.1 ) then ; kk = 1
        else                                            ; kk = 0
        endif

        ddata = zero
        ! inner boundary Lagrange multiplier terms
        ;                  ; id = Ate(lvol,0,ii)%i(ll) ; jd = Lma(lvol,  ii) ; adata = +      TTMdata(ll, mi)
        if (id.ne.0 .and. jd.ne.0)   call push_back(1,nqueue,nmaxqueue,adata,ddata,jd,dMASqueue,dMDSqueue,jdMASqueue)
        ;                  ; id = Aze(lvol,0,ii)%i(ll) ; jd = Lmb(lvol,  ii) ; adata = +      TTdata(ll, mi,kk)
        if (id.ne.0 .and. jd.ne.0)   call push_back(2,nqueue,nmaxqueue,adata,ddata,jd,dMASqueue,dMDSqueue,jdMASqueue)

        ;                  ; id = Ato(lvol,0,ii)%i(ll) ; jd = Lmc(lvol,  ii) ; adata = +      TTMdata(ll, mi)
        if (id.ne.0 .and. jd.ne.0)   call push_back(3,nqueue,nmaxqueue,adata,ddata,jd,dMASqueue,dMDSqueue,jdMASqueue)
        ;                  ; id = Azo(lvol,0,ii)%i(ll) ; jd = Lmd(lvol,  ii) ; adata = +      TTdata(ll, mi,0)
        if (id.ne.0 .and. jd.ne.0)   call push_back(4,nqueue,nmaxqueue,adata,ddata,jd,dMASqueue,dMDSqueue,jdMASqueue)

        ! outer boundary Lagrange multiplier terms & fluxes
        if( ii.gt.1 ) then ; id = Ate(lvol,0,ii)%i(ll) ; jd = Lme(lvol,  ii) ; adata = - ni * TTdata(ll, mi, 1)
          if (id.ne.0 .and. jd.ne.0) call push_back(1,nqueue,nmaxqueue,adata,ddata,jd,dMASqueue,dMDSqueue,jdMASqueue)
        ;                  ; id = Aze(lvol,0,ii)%i(ll) ; jd = Lme(lvol,  ii) ; adata = - mi * TTdata(ll, mi, 1)
          if (id.ne.0 .and. jd.ne.0) call push_back(2,nqueue,nmaxqueue,adata,ddata,jd,dMASqueue,dMDSqueue,jdMASqueue)
        ;                  ; id = Ato(lvol,0,ii)%i(ll) ; jd = Lmf(lvol,  ii) ; adata = + ni * TTdata(ll, mi, 1)
          if (id.ne.0 .and. jd.ne.0) call push_back(3,nqueue,nmaxqueue,adata,ddata,jd,dMASqueue,dMDSqueue,jdMASqueue)
        ;                  ; id = Azo(lvol,0,ii)%i(ll) ; jd = Lmf(lvol,  ii) ; adata = + mi * TTdata(ll, mi, 1)
          if (id.ne.0 .and. jd.ne.0) call push_back(4,nqueue,nmaxqueue,adata,ddata,jd,dMASqueue,dMDSqueue,jdMASqueue)

        else               ; id = Ate(lvol,0,ii)%i(ll) ; jd = Lmg(lvol,  ii) ; adata = +      TTdata(ll, mi, 1)
          if (id.ne.0 .and. jd.ne.0) call push_back(1,nqueue,nmaxqueue,adata,ddata,jd,dMASqueue,dMDSqueue,jdMASqueue)
        ;                  ; id = Aze(lvol,0,ii)%i(ll) ; jd = Lmh(lvol,  ii) ; adata = +      TTdata(ll, mi, 1)
          if (id.ne.0 .and. jd.ne.0) call push_back(2,nqueue,nmaxqueue,adata,ddata,jd,dMASqueue,dMDSqueue,jdMASqueue)
        endif
        
        ! putting things in the sparse matrix
        call addline(nqueue, nmaxqueue, dMASqueue, dMDSqueue, jdMASqueue, ns, nrow, dMAS, dMDS, jdMAS, idMAS)
        
      enddo ! end of do ll ;
 
    enddo ! end of do ii ;

  endif ! if (YESstellsym)

  ! deal with rest of the columes

  do ii = 1, mn; mi = im(ii) ; ni = in(ii)

    if( Lcoordinatesingularity .and. ii.eq.1 ) then ; kk = 1
    else                                            ; kk = 0
    endif

    call clean_queue(nqueue, nmaxqueue, dMASqueue, dMDSqueue, jdMASqueue)
    ! the inner boundary condition terms, YESstellsym
    do ll = 0, lrad ; jd = Ate(lvol,0,ii)%i(ll) ; id = Lma(lvol,  ii) ; adata = +      TTMdata(ll, mi)
    if (id.ne.0 .and. jd.ne.0)   call push_back(1,nqueue,nmaxqueue,adata,zero,jd,dMASqueue,dMDSqueue,jdMASqueue)
    enddo
    do ll = 0, lrad ; jd = Aze(lvol,0,ii)%i(ll) ; id = Lmb(lvol,  ii) ; adata = +      TTdata(ll, mi,kk)
    if (id.ne.0 .and. jd.ne.0)   call push_back(2,nqueue,nmaxqueue,adata,zero,jd,dMASqueue,dMDSqueue,jdMASqueue)
    enddo
    call addline(nqueue, nmaxqueue, dMASqueue, dMDSqueue, jdMASqueue, ns, nrow, dMAS, dMDS, jdMAS, idMAS)

    ! the outerboundary and flux terms, YESstellsym
    call clean_queue(nqueue, nmaxqueue, dMASqueue, dMDSqueue, jdMASqueue)
    if (ii.gt.1) then
      do ll = 0, lrad ; jd = Ate(lvol,0,ii)%i(ll) ; id = Lme(lvol,  ii) ; adata = - ni * TTdata(ll, mi, 1)
      if (id.ne.0 .and. jd.ne.0)   call push_back(1,nqueue,nmaxqueue,adata,zero,jd,dMASqueue,dMDSqueue,jdMASqueue)
      ;               ; jd = Aze(lvol,0,ii)%i(ll) ; id = Lme(lvol,  ii) ; adata = - mi * TTdata(ll, mi, 1)
      if (id.ne.0 .and. jd.ne.0)   call push_back(1,nqueue,nmaxqueue,adata,zero,jd,dMASqueue,dMDSqueue,jdMASqueue)
      enddo
    else
      do ll = 0, lrad ; jd = Ate(lvol,0,ii)%i(ll) ; id = Lmg(lvol,  ii) ; adata = +      TTdata(ll, mi, 1)
      if (id.ne.0 .and. jd.ne.0)   call push_back(1,nqueue,nmaxqueue,adata,zero,jd,dMASqueue,dMDSqueue,jdMASqueue)
      ;               ; jd = Aze(lvol,0,ii)%i(ll) ; id = Lmh(lvol,  ii) ; adata = +      TTdata(ll, mi, 1)
      if (id.ne.0 .and. jd.ne.0)   call push_back(2,nqueue,nmaxqueue,adata,zero,jd,dMASqueue,dMDSqueue,jdMASqueue)
      enddo
    endif
    call addline(nqueue, nmaxqueue, dMASqueue, dMDSqueue, jdMASqueue, ns, nrow, dMAS, dMDS, jdMAS, idMAS)

    if (NOTstellsym) then
      call clean_queue(nqueue, nmaxqueue, dMASqueue, dMDSqueue, jdMASqueue)
      ! the inner boundary condition terms, NOTstellsym
      do ll = 0, lrad ; jd = Ato(lvol,0,ii)%i(ll) ; id = Lmc(lvol,  ii) ; adata = +      TTMdata(ll, mi)
      if (id.ne.0 .and. jd.ne.0)   call push_back(1,nqueue,nmaxqueue,adata,zero,jd,dMASqueue,dMDSqueue,jdMASqueue)
      enddo
      do ll = 0, lrad ; jd = Azo(lvol,0,ii)%i(ll) ; id = Lmd(lvol,  ii) ; adata = +      TTdata(ll, mi,0)
      if (id.ne.0 .and. jd.ne.0)   call push_back(2,nqueue,nmaxqueue,adata,zero,jd,dMASqueue,dMDSqueue,jdMASqueue)
      enddo
      call addline(nqueue, nmaxqueue, dMASqueue, dMDSqueue, jdMASqueue, ns, nrow, dMAS, dMDS, jdMAS, idMAS)

      ! the outer boundary and flux terms, NOTstellsym
      call clean_queue(nqueue, nmaxqueue, dMASqueue, dMDSqueue, jdMASqueue)
      if (ii.gt.1) then
        do ll = 0, lrad ; jd = Ato(lvol,0,ii)%i(ll) ; id = Lmf(lvol,  ii) ; adata = + ni * TTdata(ll, mi, 1)
        if (id.ne.0 .and. jd.ne.0)   call push_back(1,nqueue,nmaxqueue,adata,zero,jd,dMASqueue,dMDSqueue,jdMASqueue)
        ;               ; jd = Azo(lvol,0,ii)%i(ll) ; id = Lmf(lvol,  ii) ; adata = + mi * TTdata(ll, mi, 1)
        if (id.ne.0 .and. jd.ne.0)   call push_back(1,nqueue,nmaxqueue,adata,zero,jd,dMASqueue,dMDSqueue,jdMASqueue)
        enddo
      endif
      call addline(nqueue, nmaxqueue, dMASqueue, dMDSqueue, jdMASqueue, ns, nrow, dMAS, dMDS, jdMAS, idMAS)
    endif ! if (NOTstellsym)
  enddo

  NdMAS(lvol) = ns
  idMAS(NN+1) = ns+1

  ! write(ounit,*) nrow, NN
  ! write(ounit,*) dMAS(1:ns)
  ! write(ounit,*) jdMAS(1:ns)
  ! write(ounit,*) idMAS(1:NN+1)
  ! stop
  
  ! dMB and dMG are constructed elsewhere

  DALLOCATE( dMASqueue )
  DALLOCATE( dMDSqueue )
  DALLOCATE( jdMASqueue )

  DALLOCATE( TTdata )
  DALLOCATE( TTMdata )
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(spsmat)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine spsmat

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! some subroutines to help construct the sparse matrix

subroutine push_back(iq, nq, NN, vA, vD, vjA, qA, qD, qjA)
  ! push a new element at the back of the queue
  ! INPUTS:
  ! iq - INTEGER, which queue (1-4)
  ! nq - INTEGER(4), number of element in the queue
  ! NN - INTEGER, max length of matrix
  ! vA, vD - REAL, values to put in
  ! vjA - INTEGER, column index to put in
  ! qA, qD, qjA - the queues
  use constants, only : zero
  implicit none

  REAL, INTENT(IN)       :: vA, vD
  INTEGER, INTENT(IN)    :: vjA, iq, NN
  REAL, INTENT(INOUT)    :: qA(NN,4), qD(NN,4)
  INTEGER, INTENT(INOUT) :: qjA(NN,4), nq(4)
  
  if (abs(vA).gt.zero .or. abs(vD).gt.zero) then

    nq(iq) = nq(iq) + 1
    qA(nq(iq),iq) = vA
    qD(nq(iq),iq) = vD
    qjA(nq(iq),iq) = vjA

  end if
  return
end subroutine push_back

subroutine clean_queue(nq, NN, qA, qD, qjA)
  ! clean the queue
  use constants, only : zero
  implicit none

  INTEGER, INTENT(IN)    :: NN
  REAL, INTENT(INOUT)    :: qA(NN,4), qD(NN,4)
  INTEGER, INTENT(INOUT) :: qjA(NN,4), nq(4)

  nq = 0
  qA = zero
  qD = zero
  qjA = 0

  return
end subroutine clean_queue

subroutine addline(nq, NN, qA, qD, qjA, ns, nrow, dMAS, dMDS, jdMAS, idMAS)
  ! add the content from the queue to the real matrices
  implicit none

  INTEGER, INTENT(INOUT)    :: NN, ns, nrow
  REAL, INTENT(INOUT)    :: qA(NN,4), qD(NN,4)
  INTEGER, INTENT(INOUT) :: qjA(NN,4), nq(4)
  REAL :: dMAS(*), dMDS(*)
  INTEGER :: jdMAS(*), idMAS(*)

  INTEGER :: pp

  do pp = 1, 4
    if (nq(pp) .eq. 0) cycle
    nrow = nrow + 1
    dMAS(ns+1:ns+nq(pp)) = qA(1:nq(pp),pp)
    dMDS(ns+1:ns+nq(pp)) = qD(1:nq(pp),pp)
    jdMAS(ns+1:ns+nq(pp)) = qjA(1:nq(pp),pp)
    idMAS(nrow) = ns + 1
    ns = ns + nq(pp)
  end do ! pp = 1:4

end subroutine addline