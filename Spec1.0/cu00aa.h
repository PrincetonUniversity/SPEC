!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Constructs loop current and plasma current.

!latex \end{enumerate} \subsection{loop integrals} \begin{enumerate}

!latex \item The currents are given by loop integrals: $\int_{{\cal S}} {\bf j}\cdot d{\bf s} = \int_{\partial {\cal S}} {\bf B}\cdot d{\bf l}$.

!latex \item Choosing a loop $d{\bf l}\equiv {\bf e}_\z \; d\z$ we find the ``poloidal'' current linking the torus:
!latex       \be \verb+Jt+ \equiv \int B_\z \; d\z = \int d\z \left( -A_\z^\prime \; g_{\t\z} + A_\t^\prime \; g_{\z\z} \right) / \sqrt g
!latex       \ee

!latex \item Choosing a loop $d{\bf l}\equiv {\bf e}_\t \; d\t$ we find the ``toroidal'' plasma current:
!latex       \be \verb+Jz+ \equiv \int B_\z \; d\t = \int d\t \left( -A_\z^\prime \; g_{\t\t} + A_\t^\prime \; g_{\z\t} \right) / \sqrt g
!latex       \ee

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine cu00aa( lvol, Ntz, mn )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half, one, two, pi2
  
  use numerical, only :
  
  use fileunits, only : ounit
  
  use inputlist, only : Wcu00aa, Lrad, Nvol
  
  use cputiming, only : Tcu00aa
  
  use allglobal, only : myid, ncpu, cpus, Mvol, pi2nfp, im, in, halfmm, &
                        sg, guvij, Rij, Zij, ijreal, ijimag, efmn, ofmn, cfmn, sfmn, &
                        Ate, Aze, Ato, Azo, TTll, &
                        isr, trigm, trign, trigwk, Nt, Nz, &
                        Lplasmaregion, Lvacuumregion, Lcoordinatesingularity, NOTstellsym
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in)  :: lvol, Ntz, mn

  
  INTEGER              :: Lcurvature, ii, ll, innout, ifail, ideriv
  REAL                 :: lss, mfactor, dAt(1:Ntz), dAz(1:Ntz), Jzemn(1:mn), Jzomn(1:mn), Jtemn(1:mn), Jtomn(1:mn), lcurtor, lcurpol

  BEGIN(cu00aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  FATALMESS(cu00aa, lvol.lt.1 .or. lvol.gt.Mvol, invalid volume )
#endif
  
  ideriv = 0
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Lplasmaregion ) then
   
   do innout = 0, 1
    
    if( Lcoordinatesingularity .and. innout.eq.0 ) cycle

    lss = two * innout - one
    
    Lcurvature = 1
    WCALL(cu00aa, co01aa,( lvol, lss, Lcurvature, Ntz, mn )) ! get coordinates and derivatives wrt Rj, Zj, at specific radial location;
    
    efmn(1:mn) = zero ; ofmn(1:mn) = zero ; cfmn(1:mn) = zero ; sfmn(1:mn) = zero ! ! initialize summation; 24 Apr 13;
    
    do ii = 1, mn ! compute radial derivatives of vector potential; 20 Feb 13;
     
     if( Lcoordinatesingularity ) then ; mfactor = halfmm(ii) * half
     else                              ; mfactor = zero
     endif
     
     do ll = 0, Lrad(lvol)
      ;                      ; efmn(ii) = efmn(ii) + Ate(lvol,ideriv,ii)%s(ll) * ( TTll(ll,innout,1) + mfactor )
      ;                      ; cfmn(ii) = cfmn(ii) + Aze(lvol,ideriv,ii)%s(ll) * ( TTll(ll,innout,1) + mfactor )
      if( NOTstellsym ) then ; ofmn(ii) = ofmn(ii) + Ato(lvol,ideriv,ii)%s(ll) * ( TTll(ll,innout,1) + mfactor )
       ;                     ; sfmn(ii) = sfmn(ii) + Azo(lvol,ideriv,ii)%s(ll) * ( TTll(ll,innout,1) + mfactor )
!     else                   ; ofmn(ii) = zero
!      ;                     ; sfmn(ii) = zero
      endif
     enddo ! end of do ll; 20 Feb 13;
    
    enddo ! end of do ii; 20 Feb 13;
   
    call invfft( mn, im, in, efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, dAt(1:Ntz), dAz(1:Ntz), isr, trigm, trign, trigwk ) ! map to real space; 20 Feb 13;
    
    ijreal(1:Ntz) = ( - dAz(1:Ntz) * guvij(1:Ntz,2,2,0) + dAt(1:Ntz) * guvij(1:Ntz,2,3,0) ) / sg(1:Ntz,0)
    ijimag(1:Ntz) = ( - dAz(1:Ntz) * guvij(1:Ntz,2,3,0) + dAt(1:Ntz) * guvij(1:Ntz,3,3,0) ) / sg(1:Ntz,0)
    
    ifail = 0
    call tfft( Nt, Nz, ijreal(1:Ntz), ijimag(1:Ntz), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz), & ! compute force-imbalance and spectral constraints; 20 Feb 13;
               mn, im(1:mn), in(1:mn), Jzemn(1:mn), Jzomn(1:mn), Jtemn(1:mn), Jtomn(1:mn), ifail )

#ifdef DEBUG
    if( Wcu00aa ) then
     cput = GETTIME
     write(ounit,'("cu00aa : ",f10.2," : myid=",i3," ; lvol=",i3," ; innout="i2" ; Jtemn="99es13.5)') cput-cpus, myid, lvol, innout, Jtemn(1:mn)
     write(ounit,'("cu00aa : ",f10.2," : myid=",i3," ; lvol=",i3," ; innout="i2" ; Jzemn="99es13.5)') cput-cpus, myid, lvol, innout, Jzemn(1:mn)
    endif
#endif
    
    if( lvol.eq.Nvol .and. innout.eq.1 ) then
     lcurpol = Jtemn(1) ; lcurtor = Jzemn(1)
     write(ounit,'("cu00aa : ", 10x ," : ")')
     write(ounit,'("cu00aa : ", 10x ," : myid=",i3," ; curpol="es23.15" ; curtor="es23.15" ;")') myid, (/ lcurpol, lcurtor /)
    endif

   enddo ! ! end of do innout = 0, 1 ; 24 Apr 13;
   
  else ! ! matches if( Lplasmregion ) ; 24 Apr 13;
   
   cput = GETTIME
   write(ounit,'("cu00aa : ",f10.2," : myid=",i3," ; lvol=",i3," ; vacuum region under construction ;")') cput-cpus, myid, lvol

  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(cu00aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine cu00aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
