!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Packs, and unpacks, geometrical degrees of freedom.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \end{enumerate} \subsection{Geometrical degrees of freedom} \begin{enumerate}

!latex \item The geometrical degrees-of-freedom are represented by the vector $\boldxi$.

!latex \item This represents the $R_{e,i,j}$ and $Z_{o,i,j}$ harmonics,
!latex       where $i$ labels the interface and $j$ labels the Fourier harmonic.

!latex \item A factor of $\psi_{t,i}^{m_j/2}$ is included, e.g.
!latex       \be \boldxi_k \equiv \frac{R_{e,i,j}}{\psi_{t,i}^{m_j/2}},
!latex       \ee
!latex       where $\psi_{t,i}$ is the toroidal flux enclosed by the $i$-th interface, and similarly for the other harmonics.

!latex \item The coordinate axis, if required, is set according to 
!latex       \be R_{j,0} & = & \sum_{i} R_{i,1} |\beta_{i,j}| \\
!latex           Z_{j,0} & = & \sum_{i} Z_{i,1}  \beta_{i,j}  \\
!latex       \ee
!latex       where the $\beta_{i,j}$ are set in \verb+global+.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine gf00aa( Ngeometricaldof, position, Mvol, mn, iRbc, iZbs, iRbs, iZbc, packorunpack )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, one
  
  use numerical, only :
  
  use fileunits, only : ounit
  
  use inputlist, only : Wgf00aa, Igeometry, Ntor, tflux
  
  use cputiming, only : Tgf00aa
  
  use allglobal, only : ncpu, myid, cpus, im, in, halfmm, YESstellsym, NOTstellsym, &
                        bjk, &
                        trigm, trign, trigwk, isr, Nt, Nz, Ntz, iRij, iZij, tRij, tZij, &
                        ijreal, ijimag, jireal, jiimag, efmn, ofmn, cfmn, sfmn, evmn, odmn, comn, simn, psifactor
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in)    :: Ngeometricaldof, Mvol, mn
  REAL                   :: position(0:Ngeometricaldof), iRbc(1:mn,0:Mvol), iZbs(1:mn,0:Mvol), iRbs(1:mn,0:Mvol), iZbc(1:mn,0:Mvol)
  CHARACTER              :: packorunpack
  
  INTEGER                :: lvol, ii, kk, irz, issym, igeometricaldegreeoffreedom, ifail
  
  BEGIN(gf00aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  igeometricaldegreeoffreedom = 0 ! initialize counter; 14 Jan 13;
  
  
  do lvol = 1, Mvol-1 ! loop over internal interfaces;
   
   do ii = 1, mn ! loop over Fourier harmonics;
    
    do irz = 0, 1 ! loop over R & Z;
     
     if( Igeometry.lt.3 .and. irz.eq.1 ) cycle ! no dependence on Z; 14 Jan 13;
     
     do issym = 0, 1 ! loop over even & odd;
      
      if( YESstellsym .and. issym.eq.1 ) cycle
      
      if( issym.eq.0 .and. irz.eq.1 .and. ii.eq.1 ) cycle ! no dependence on Zbs_{0,0}; 14 Jan 13;
      if( issym.eq.1 .and. irz.eq.0 .and. ii.eq.1 ) cycle ! no dependence on Rbs_{0,0}; 14 Jan 13;
      
      igeometricaldegreeoffreedom = igeometricaldegreeoffreedom + 1
      
#ifdef DEBUG
      FATALMESS(gf00aa, igeometricaldegreeoffreedom.le.0 .or. igeometricaldegreeoffreedom.gt.Ngeometricaldof, out of bounds)
#endif
      
      select case( packorunpack )
       
      case( 'P' ) !   pack vector of unknowns;
       
       if( irz.eq.0 .and. issym.eq.0 ) position(igeometricaldegreeoffreedom) = iRbc(ii,lvol) / psifactor(ii,lvol)
       if( irz.eq.1 .and. issym.eq.0 ) position(igeometricaldegreeoffreedom) = iZbs(ii,lvol) / psifactor(ii,lvol)       
       if( irz.eq.0 .and. issym.eq.1 ) position(igeometricaldegreeoffreedom) = iRbs(ii,lvol) / psifactor(ii,lvol)
       if( irz.eq.1 .and. issym.eq.1 ) position(igeometricaldegreeoffreedom) = iZbc(ii,lvol) / psifactor(ii,lvol)
       
      case( 'U' ) ! unpack vector of unknowns;
       
       if( irz.eq.0 .and. issym.eq.0 ) iRbc(ii,lvol) = position(igeometricaldegreeoffreedom) * psifactor(ii,lvol)
       if( irz.eq.1 .and. issym.eq.0 ) iZbs(ii,lvol) = position(igeometricaldegreeoffreedom) * psifactor(ii,lvol)
       if( irz.eq.0 .and. issym.eq.1 ) iRbs(ii,lvol) = position(igeometricaldegreeoffreedom) * psifactor(ii,lvol)
       if( irz.eq.1 .and. issym.eq.1 ) iZbc(ii,lvol) = position(igeometricaldegreeoffreedom) * psifactor(ii,lvol)
       
      case default
       
       FATALMESS(gf00aa, .true., given value of packorunpack is not supported)
       
      end select
      
     enddo ! end of do issym;
     
    enddo ! end of do irz;
    
   enddo ! end of do ii;
   
  enddo ! end of do lvol;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! if( packorunpack.eq.'U' ) then ! set degenerate surface = coordinate axis; 18 Jul 14;
  
  select case( Igeometry ) 
   
  case( 1 )
   
   iRbc(1:mn,0) = zero
  !iZbs(1:mn,0) = zero
   if( NOTstellsym ) then
   iRbs(1:mn,0) = zero
   !iZbc(1:mn,0) = zero
   endif
    
  case( 2 )
   
   iRbc(1:mn,0) = zero
  !iZbs(1:mn,0) = zero
   if( NOTstellsym ) then
   iRbs(1:mn,0) = zero
   !iZbc(1:mn,0) = zero
   endif
    
  case( 3 )
    
!#ifdef NEWAXIS
!   
!   call invfft( mn, im(1:mn), in(1:mn),            iRbc(1:mn,1),              iRbs(1:mn,1), &
!                                                   iZbc(1:mn,1),              iZbs(1:mn,1), &
!                Nt, Nz, iRij(1:Ntz,1), iZij(1:Ntz,1), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz) )
!
!   call invfft( mn, im(1:mn), in(1:mn), im(1:mn) * iRbs(1:mn,1), - im(1:mn) * iRbc(1:mn,1), &
!                                        im(1:mn) * iZbs(1:mn,1), - im(1:mn) * iZbc(1:mn,1), &
!                Nt, Nz, tRij(1:Ntz,1), tZij(1:Ntz,1), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz) )
!
!   ijreal(1:Ntz) = sqrt( tRij(1:Ntz,1)**2 + tZij(1:Ntz,1)**2 ) ! differential poloidal length; 11 Aug 14;
!   ijimag(1:Ntz) = ijreal(1:Ntz)
!   
!   jireal(1:Ntz) = ijreal(1:Ntz)
!
!   ifail = 0
!   call tfft( Nt, Nz, ijreal(1:Ntz), ijimag(1:Ntz), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz), & ! compute force-imbalance and spectral constraints; 20 Feb 13;
!              mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), ifail )
!
!   ijreal(1:Ntz) = jireal(1:Ntz) * iRij(1:Ntz,1) ! weighted differential poloidal length; 11 Aug 14;
!   ijimag(1:Ntz) = jireal(1:Ntz) * iZij(1:Ntz,1) ! weighted differential poloidal length; 11 Aug 14;
!
!   ifail = 0
!   call tfft( Nt, Nz, ijreal(1:Ntz), ijimag(1:Ntz), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz), & ! compute force-imbalance and spectral constraints; 20 Feb 13;
!              mn, im(1:mn), in(1:mn), evmn(1:mn), odmn(1:mn), comn(1:mn), simn(1:mn), ifail )
!
!   efmn(1:mn) = efmn(1:mn) * ajk(1:mn) ! poloidal integration; only take m=0 harmonics; 11 Aug 14;
!   ofmn(1:mn) = zero
!   cfmn(1:mn) = cfmn(1:mn) * ajk(1:mn) ! poloidal integration; only take m=0 harmonics; 11 Aug 14;
!   sfmn(1:mn) = zero
!   evmn(1:mn) = evmn(1:mn) * ajk(1:mn) ! poloidal integration; only take m=0 harmonics; 11 Aug 14;
!   odmn(1:mn) = zero
!   comn(1:mn) = zero
!   simn(1:mn) = simn(1:mn) * ajk(1:mn) ! poloidal integration; only take m=0 harmonics; 11 Aug 14;
!   
!   call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), &
!                Nt, Nz, ijreal(1:Ntz), ijimag(1:Ntz), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz) )
!   
!   call invfft( mn, im(1:mn), in(1:mn), evmn(1:mn), odmn(1:mn), comn(1:mn), simn(1:mn), &
!                Nt, Nz, jireal(1:Ntz), jiimag(1:Ntz), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz) )
!
!   ijreal(1:Ntz) = jireal(1:Ntz) / ijimag(1:Ntz)
!   ijimag(1:Ntz) = jiimag(1:Ntz) / ijimag(1:Ntz)
!
!   ifail = 0
!   call tfft( Nt, Nz, ijreal(1:Ntz), ijimag(1:Ntz), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz), & ! compute force-imbalance and spectral constraints; 20 Feb 13;
!              mn, im(1:mn), in(1:mn), iRbc(1:mn,0), iRbs(1:mn,0), iZbc(1:mn,0), iZbs(1:mn,0), ifail )
!
!   write(ounit,'("gf00aa : ",f10.2," : myid=",i3," ; axis :      ",999(i12   ,","))') cput-cpus, myid,   im(1:mn)
!   write(ounit,'("gf00aa : ",f10.2," : myid=",i3," ; axis :      ",999(i12   ,","))') cput-cpus, myid,   in(1:mn)
!   write(ounit,'("gf00aa : ",f10.2," : myid=",i3," ; axis : Rbc=[",999(f12.08,","))') cput-cpus, myid, iRbc(1:mn,0)
!   write(ounit,'("gf00aa : ",f10.2," : myid=",i3," ; axis : Zbs=[",999(f12.08,","))') cput-cpus, myid, iZbs(1:mn,0)
!   
!#endif
   
   do kk = 1, mn ! see also global:writin; 11 Aug 14;
    iRbc(kk,0) = sum( iRbc(1:mn,1) * abs( bjk(1:mn,kk) ) ) ! sets coordinate axis as midpoint of innermost interface; 04 Dec 14;
    iZbs(kk,0) = sum( iZbs(1:mn,1) *      bjk(1:mn,kk)   ) ! sets coordinate axis as midpoint of innermost interface; 04 Dec 14;
    iRbs(kk,0) = sum( iRbs(1:mn,1) *      bjk(1:mn,kk)   ) ! sets coordinate axis as midpoint of innermost interface; 04 Dec 14;
    iZbc(kk,0) = sum( iZbc(1:mn,1) * abs( bjk(1:mn,kk) ) ) ! sets coordinate axis as midpoint of innermost interface; 04 Dec 14;
   enddo
   
!#ifdef NEWAXIS
!   cput = GETTIME
!   write(ounit,'("gf00aa : ", 10x ," : ")')
!   write(ounit,'("gf00aa : ",f10.2," : myid=",i3," ; orig : Rbc=[",999(f12.08,","))') cput-cpus, myid, iRbc(1:mn,1)
!   write(ounit,'("gf00aa : ",f10.2," : myid=",i3," ; orig : Zbs=[",999(f12.08,","))') cput-cpus, myid, iZbs(1:mn,1)
!   write(ounit,'("gf00aa : ",f10.2," : myid=",i3," ; axis : Rbc=[",999(f12.08,","))') cput-cpus, myid, iRbc(1:mn,0)
!   write(ounit,'("gf00aa : ",f10.2," : myid=",i3," ; axis : Zbs=[",999(f12.08,","))') cput-cpus, myid, iZbs(1:mn,0)
!#endif
   
  case( 4 )
   
   do kk = 1, mn ! see also global:writin; 11 Aug 14;
    iRbc(kk,0) = sum( iRbc(1:mn,1) * abs( bjk(1:mn,kk) ) ) ! sets coordinate axis as midpoint of innermost interface; 04 Dec 14;
    iZbs(kk,0) = sum( iZbs(1:mn,1) *      bjk(1:mn,kk)   ) ! sets coordinate axis as midpoint of innermost interface; 04 Dec 14;
    iRbs(kk,0) = sum( iRbs(1:mn,1) *      bjk(1:mn,kk)   ) ! sets coordinate axis as midpoint of innermost interface; 04 Dec 14;
    iZbc(kk,0) = sum( iZbc(1:mn,1) * abs( bjk(1:mn,kk) ) ) ! sets coordinate axis as midpoint of innermost interface; 04 Dec 14;
   enddo
   
  case default
   
   FATALMESS(gf00aa, .true., Igeometry not supported)
   
  end select
  
!endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  FATALMESS(gf00aa, igeometricaldegreeoffreedom.ne.Ngeometricaldof, counting error)
#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(gf00aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine gf00aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


