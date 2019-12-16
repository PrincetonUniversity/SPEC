!> \defgroup free-boundary Free-Boundary Computation

!> \file bnorml.f90
!! \brief Computes \f${\bf B}_{Plasma} \cdot {\bf e}_\theta \times {\bf e}_\zeta \;\f$ on the computational boundary, \f$\partial {\cal D}\f$.
!! \ingroup free-boundary
!!
!! **free-boundary constraint**
!! <ul>
!! <li> The normal field at the computational boundary, \f$\partial {\cal D}\f$, should be equal to 
!!      \f$\left({\bf B}_P + {\bf B}_C\right)\cdot {\bf e}_\theta \times {\bf e}_\zeta\f$,
!!      where \f${\bf B}_P\f$ is the "plasma" field (produced by internal plasma currents) and is computed using virtual casing, 
!!      and \f${\bf B}_C\f$ is the "vacuum" field (produced by the external coils) and is given on input. </li>
!! <li> The plasma field, \f${\bf B}_P\f$, can only be computed after the equilibrium is determined, 
!!      but this information is required to compute the equilibrium to begin with; and so there is an iteration involved. </li>
!! <li> Suggested values of the vacuum field can be self generated; see xspech.f90 for more documentation on this. </li>
!! </ul>
!!
!! **compute the normal field on a regular grid on the computational boundary**
!! <ul>
!! <li> For each point on the compuational boundary, \ref casing() is called to compute the normal field produced by the plasma currents. </li>
!! <li> (There is a very clumsy attempt to parallelize this which could be greatly improved.) </li>
!! <li> An FFT gives the required Fourier harmonics. </li>
!! </ul>
!! @see casing.f90


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!l tex \subsection{construction of normal field} 

!l tex \begin{enumerate}
!l tex \item The area-weighted normal vector to the computational domain is given as follows:
!l tex       \bi
!l tex       \item[] \inputvar{Igeometry.eq.1} : Cartesian \\
!l tex       $\bx = \t \; \hat i + \z \; \hat j + R(\t,\z) \; \hat k$ \\
!l tex       ${\bf e}_\t \times {\bf e}_\z = - R_\t \; \hat i - R_\z \; \hat j + \hat k$
!l tex       \item[] \inputvar{Igeometry.eq.2} : Cylindrical \\
!l tex       \item[] \inputvar{Igeometry.eq.3} : Toroidal \\
!l tex       $\bx = R(\t,\z) \cos \z \; \hat i + R(\t,\z) \sin \z \; \hat j + Z(\t,\z) \; \hat k$ \\
!l tex       ${\bf e}_\t \times {\bf e}_\z = - R \, Z_\theta \, \hat r + (Z_\theta \,R_\zeta - R_\theta \,Z_\zeta) \hat \phi + R \,R_\theta \,\hat z$
!l tex       \ei
!l tex \item NOTE: it is ${\bf e}_\t \times {\bf e}_\z$ that is required, not the unit normal vector, 
!l tex       ${\bf n} \equiv ( {\bf e}_\t \times {\bf e}_\z ) / | {\bf e}_\t \times {\bf e}_\z |$.
!l tex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!l tex \subsection{outline} 

!l tex \begin{enumerate}
!l tex \item The computational boundary is obtained using \link{coords}.
!l tex       (Note that the computational boundary does not change, so this needs only to be determined once.)
!l tex \item At each point on the computational boundary (i.e., on the discrete grid), 
!l tex       \link{casing} is used to compute the plasma field using the virtual casing principle.
!l tex       I think that \link{casing} returns the field in Cartesian coordinates, i.e., ${\bf B} = B_x {\bf i} + B_y {\bf j} + B_z {\bf k}$.
!l tex \item In toroidal geometry, the vector transformation from Cartesian to cylindrical is given by
!l tex       \be \begin{array}{cccccccccccccccccccccc}
!l tex           B^R    & = &   & + B_x \cos \z & + & B_y \sin \z &       & \\
!l tex           B^\phi & = & ( & - B_x \sin \z & + & B_y \cos \z & ) / R & \\
!l tex           B^Z    & = &   &               &   &             &       & B_z
!l tex           \end{array}
!l tex       \ee
!l tex       so that ${\bf B} = B^R {\bf e}_R + B^\phi {\bf e}_\phi + B^Z {\bf e}_Z$.
!l tex \end{enumerate}
!!
!!  -!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> @brief Computes \f${\bf B}_{Plasma} \cdot {\bf e}_\theta \times {\bf e}_\zeta \;\f$ on the computational boundary, \f$\partial {\cal D}\f$.
!!
!! Computes \f${\bf B}_{Plasma} \cdot {\bf e}_\theta \times {\bf e}_\zeta \;\f$ on the computational boundary, \f$\partial {\cal D}\f$.
!!
!! @param[in] mn total number of Fourier harmonics
!! @param[in] Ntz total number of grid points in \f$\theta\f$ and \f$zeta\f$
!! @param[out] efmn even Fourier coefficients
!! @param[out] ofmn odd Fouier coefficients
subroutine bnorml( mn, Ntz, efmn, ofmn )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half, one, two, pi, pi2, ten
  
  use numerical, only : small
  
  use fileunits, only : ounit, lunit
  
  use inputlist, only : Wmacros, Wbnorml, Igeometry, Lcheck, vcasingtol, vcasingper, Lrad

  use cputiming, only : Tbnorml
  
  use allglobal, only : ncpu, myid, cpus, pi2nfp, Mvol, &
                        Nt, Nz, &
                        Rij, Zij, guvij, sg, TT, &
                        NOTstellsym, Lcoordinatesingularity, &
                        im, in, Ate, Aze, Ato, Azo, &
                        Nt, Nz, cfmn, sfmn, &
                        ijreal, ijimag, jireal, jiimag, &
                        globaljk, tetazeta, virtualcasingfactor, gteta, gzeta, Dxyz, Nxyz
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in)  :: mn, Ntz
  REAL   , intent(out) :: efmn(1:mn), ofmn(1:mn)
  
  INTEGER              :: lvol, Lcurvature, Lparallel, ii, jj, kk, jk, ll, kkmodnp, jkmodnp, ifail, id01daf, nvccalls, icasing, ideriv
  REAL                 :: lss, zeta, teta, cszeta(0:1), tetalow, tetaupp, absacc, gBn
  REAL                 :: Jxyz(1:Ntz,1:3), Bxyz(1:Ntz,1:3), dAt(1:Ntz), dAz(1:Ntz), distance(1:Ntz)
  
 !REAL                 :: vcintegrand, zetalow, zetaupp
! external             :: vcintegrand, zetalow, zetaupp
  
  BEGIN(bnorml)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  Lparallel = 1 ! controls choice of parallelization; see below;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!#ifdef CASING
!  
!  lvol = Mvol ; lss = zero ; Lcurvature = 1 ; Lcoordinatesingularity = .false.
!  
!  WCALL( bnorml, coords,( lvol, lss, Lcurvature, Ntz, mn ) ) ! re-compute plasma boundary and metrics; this may not be required; 14 Apr 17;
!  
!  jireal(1:Ntz) = Rij(1:Ntz,0,0) ! save plasma boundary; 04 May 17;
!  jiimag(1:Ntz) = Zij(1:Ntz,0,0)
!
!  efmn(1:mn) = zero ; sfmn(1:mn) = zero ; cfmn(1:mn) = zero ; ofmn(1:mn) = zero
!  
!  ideriv = 1
!
!  do ii = 1, mn ! loop over Fourier harmonics; 13 Sep 13;
!   
!   do ll = 0, Lrad(Mvol) ! loop over Chebyshev polynomials; Lrad is the radial resolution;
!    ;                      ; efmn(ii) = efmn(ii) + Ate(Mvol,ideriv,ii)%s(ll) * TT(ll,0,1) ! ideriv labels deriv. wrt mu, pflux; 
!    ;                      ; cfmn(ii) = cfmn(ii) + Aze(Mvol,ideriv,ii)%s(ll) * TT(ll,0,1) 
!    if( NOTstellsym ) then ; ofmn(ii) = ofmn(ii) + Ato(Mvol,ideriv,ii)%s(ll) * TT(ll,0,1) 
!     ;                     ; sfmn(ii) = sfmn(ii) + Azo(Mvol,ideriv,ii)%s(ll) * TT(ll,0,1) 
!    endif
!   enddo ! end of do ll; 20 Feb 13;
!   
!  enddo ! end of do ii; 20 Feb 13;
!  
!  call invfft( mn, im, in, efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, dAt(1:Ntz), dAz(1:Ntz) ) ! map to real space;
!
!  ijreal(1:Ntz) = ( dAt(1:Ntz) * guvij(1:Ntz,2,3,0) + dAz(1:Ntz) * guvij(1:Ntz,3,3,0) ) / sg(1:Ntz,0) ! \alpha; 14 Apr 17;
!  ijimag(1:Ntz) = ( dAt(1:Ntz) * guvij(1:Ntz,2,2,0) + dAz(1:Ntz) * guvij(1:Ntz,2,3,0) ) / sg(1:Ntz,0) ! \beta ; 14 Apr 17;
!  
!  select case( Igeometry )
!   
!  case( 1 ) ! Igeometry = 1; 14 Apr 17;
!   
!   Bxyz(1:Ntz,1) = gteta(1:Ntz) ! plasma boundary; Cartesian; 04 May 17;
!   Bxyz(1:Ntz,2) = gzeta(1:Ntz)
!   Bxyz(1:Ntz,3) = Rij(1:Ntz,0,0)
!   
!   Jxyz(1:Ntz,1) = ijreal(1:Ntz)                                               ! surface current; Cartesian; 04 May 17;
!   Jxyz(1:Ntz,2) =                              - ijimag(1:Ntz)                ! surface current; Cartesian; 04 May 17;
!   Jxyz(1:Ntz,3) = ijreal(1:Ntz)*Rij(1:Ntz,2,0) - ijimag(1:Ntz)*Rij(1:Ntz,3,0) ! surface current; Cartesian; 04 May 17;
!   
!  case( 2 ) ! Igeometry = 2 ; 08 May 17;
!   
!   FATAL( bnorml, .true., new evaluation of virtual casing for Igeometry = 2 is under construction )
!   
!  case( 3 ) ! Igeometry = 3; 14 Apr 17;
!   
!   Bxyz(1:Ntz,1) = Rij(1:Ntz,0,0)*cos(gzeta(1:Ntz)) ! plasma boundary; Cartesian; 08 May 17;
!   Bxyz(1:Ntz,2) = Rij(1:Ntz,0,0)*sin(gzeta(1:Ntz))
!   Bxyz(1:Ntz,3) = Zij(1:Ntz,0,0)                  
!   
!   Jxyz(1:Ntz,1) = ijreal(1:Ntz)*Rij(1:Ntz,2,0)*cos(gzeta(1:Ntz)) - ijimag(1:Ntz)*(Rij(1:Ntz,3,0)*cos(gzeta(1:Ntz))-Rij(1:Ntz,0,0)*sin(gzeta(1:Ntz)))
!   Jxyz(1:Ntz,2) = ijreal(1:Ntz)*Rij(1:Ntz,2,0)*sin(gzeta(1:Ntz)) - ijimag(1:Ntz)*(Rij(1:Ntz,3,0)*sin(gzeta(1:Ntz))+Rij(1:Ntz,0,0)*cos(gzeta(1:Ntz)))
!   Jxyz(1:Ntz,3) = ijreal(1:Ntz)*Zij(1:Ntz,2,0)                   - ijimag(1:Ntz)*(Zij(1:Ntz,3,0)                                                   )
!   
!  case default
!   
!   FATAL( bnorml, .true., illegal Igeometry )
!   
!  end select ! end of select case( Igeometry ) ; 08 May 17;
!  
!  do jk = 1, Ntz
!   
!   distance(1:Ntz) = (Dxyz(jk,1)-Bxyz(1:Ntz,1))**2 + (Dxyz(jk,2)-Bxyz(1:Ntz,2))**2 + (Dxyz(jk,3)-Bxyz(1:Ntz,3))**2
!   
!   distance(1:Ntz) = sqrt(distance(1:Ntz)) * distance(1:Ntz)
!
!   ijimag(jk) = sum( ( Jxyz(1:Ntz,2)*(Dxyz(jk,3)-Bxyz(1:Ntz,3)) - Jxyz(1:Ntz,3)*(Dxyz(jk,2)-Bxyz(1:Ntz,2)) ) / distance(1:Ntz) ) * Nxyz(jk,1) &
!              + sum( ( Jxyz(1:Ntz,3)*(Dxyz(jk,1)-Bxyz(1:Ntz,1)) - Jxyz(1:Ntz,1)*(Dxyz(jk,3)-Bxyz(1:Ntz,3)) ) / distance(1:Ntz) ) * Nxyz(jk,2) &
!              + sum( ( Jxyz(1:Ntz,1)*(Dxyz(jk,2)-Bxyz(1:Ntz,2)) - Jxyz(1:Ntz,2)*(Dxyz(jk,1)-Bxyz(1:Ntz,1)) ) / distance(1:Ntz) ) * Nxyz(jk,3)
!   
!  enddo ! end of do jk = 1, Ntz; 08 May 17;
!
!  ijimag(1:Ntz) = ijimag(1:Ntz) * pi2 * pi2nfp / Ntz ! WHAT ABOUT FIELD PERIODICITY; 08 May 17;
!    
!#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  ijreal(1:Ntz) = zero ! normal plasma field; 15 Oct 12;
 !ijimag(1:Ntz) = zero
  
 !jireal(1:Ntz) = zero
 !jiimag(1:Ntz) = zero
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  if( myid.eq.0 .and. Lcheck.eq.6 ) then
   write(ounit,'("bnorml : " 10x " : writing input for xdiagno ; screen comparison ; Ntz =",i7," ;")') Ntz
   open(lunit, file="btest.diagno", status="unknown" )
   write(lunit,'(i9)') Ntz
  endif
#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  do kk = 0, Nz-1 ; zeta = kk * pi2nfp / Nz
   
   if( Igeometry.eq.3 ) then ; cszeta(0:1) = (/ cos(zeta), sin(zeta) /)
   endif
   
   do jj = 0, Nt-1 ; teta = jj * pi2    / Nt ; jk = 1 + jj + kk*Nt

    globaljk = jk ! this is global; passed through to vcintegrand & casing;

    select case( Lparallel ) ! perform in parallel;
    case( 0 ) ! Lparallel = 0 ; 09 Mar 17;
     if( myid.ne.modulo(kk,ncpu) ) cycle
    case( 1 ) ! Lparallel = 1 ; 09 Mar 17;
     if( myid.ne.modulo(jk-1,ncpu) ) cycle ! 11 Oct 12; this is a weird parallelization, but perhaps better exploits all available cpus;
    case default ! Lparallel; 09 Mar 17;
     FATAL( bnorml, .true., invalid Lparallel in parallelization loop )
    end select ! end of select case( Lparallel ) ; 09 Mar 17;

    tetazeta(1:2) = (/ teta, zeta /) ! this is global; passed through to zetalow & zetaupp; 14 Apr 17;
    
!#ifdef COMPARECASING
!    
!    tetalow = tetazeta(1) - vcasingper * pi ; tetaupp = tetazeta(1) + vcasingper * pi ; absacc = vcasingtol 
!    
!    id01daf = 1 ; call D01DAF( tetalow, tetaupp, zetalow, zetaupp, vcintegrand, absacc, gBn, nvccalls, id01daf ) ! 04 May 17;
!    
!    ijimag(jk) = gBn
!    
!#endif
    
    WCALL( bnorml, casing, ( teta, zeta, gBn, icasing ) ) ! tetazeta is global; 26 Apr 17;
    
    ijreal(jk) = gBn
    
!#ifdef COMPARECASING
!    write(ounit,1000) myid, zeta, teta, ijreal(jk), ijimag(jk), ijreal(jk)-ijimag(jk)
!#endif
!    
!#ifdef CASING
!    write(ounit,1000) myid, zeta, teta, ijreal(jk), ijimag(jk), ijreal(jk)-ijimag(jk)
!#endif
    
1000 format("bnorml : ", 10x ," : myid=",i3," : \z =",f6.3," ; \t =",f6.3," ; B . x_t x x_z =",2f22.15," ; ":"err =",es13.5," ;")
    
   enddo ! end of do jj;
  enddo ! end of do kk;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

1001 format("bnorml : ", 10x ," : "a1" : (t,z) = ("f8.4","f8.4" ) ; gBn=",f23.15," ; ":" error =",f23.15" ;")

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( myid.eq.0 .and. Lcheck.eq.6 ) then ! THIS WAS CORRUPTED; see before 14 Apr 17 for complete source; 
   close(lunit)
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  do kk = 0, Nz-1
   
   kkmodnp = modulo(kk,ncpu)
   
   select case( Lparallel )
    
   case( 0 ) ! Lparallel = 0 ; 09 Mar 17;
    
    RlBCAST(ijreal(1+kk*Nt:Nt+kk*Nt),Nt,kkmodnp) ! plasma; 03 Apr 13;
   !RlBCAST(ijimag(1+kk*Nt:Nt+kk*Nt),Nt,kkmodnp)
    
   !RlBCAST(jireal(1+kk*Nt:Nt+kk*Nt),Nt,kkmodnp)
   !RlBCAST(jiimag(1+kk*Nt:Nt+kk*Nt),Nt,kkmodnp)

   case( 1 ) ! Lparallel = 1 ; 09 Mar 17;
    
    do jj = 0, Nt-1
     
     jk = 1 + jj + kk*Nt
     
     jkmodnp = modulo(jk-1,ncpu)
     
     RlBCAST(ijreal(jk),1,jkmodnp) ! plasma; 03 Apr 13;
    !RlBCAST(ijimag(jk),1,jkmodnp)
 
    !RlBCAST(jireal(jk),1,jkmodnp)
    !RlBCAST(jiimag(jk),1,jkmodnp)
     
    enddo
    
   case default ! Lparallel; 09 Mar 17;
    
    FATAL( bnorml, .true., invalid Lparallel for broadcasting )
    
   end select ! end of select case( Lparallel ) ; 09 Mar 17;
   
  enddo ! 11 Oct 12;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  ijreal(1:Ntz) = ijreal(1:Ntz) * virtualcasingfactor
  ijimag(1:Ntz) = zero

  call tfft( Nt, Nz, ijreal(1:Ntz), ijimag(1:Ntz), &
             mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), ifail ) ! Fourier decompose normal field;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(bnorml)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine bnorml
