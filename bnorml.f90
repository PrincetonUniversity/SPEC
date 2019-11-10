!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (free-boundary) ! Computes ${\bf B}_{Plasma} \cdot {\bf e}_\theta \times {\bf e}_\zeta \;$ on the computational boundary, $\partial {\cal D}$.

!latex \briefly{Computes ${\bf B}_P \cdot {\bf e}_\t \times {\bf e}_\z$ on computational boundary, $\partial {\cal D}$.}

!latex \calledby{\link{xspech}}

!latex \calls{\link{coords} and 
!latex        \link{casing}}

!latex \tableofcontents

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{free-boundary constraint} 

!latex \begin{enumerate}
!latex \item The normal field at the computational boundary, $\partial {\cal D}$, should be equal to 
!latex       $\left({\bf B}_P + {\bf B}_C\right)\cdot {\bf e}_\theta \times {\bf e}_\zeta$,
!latex       where ${\bf B}_P$ is the ``plasma'' field (produced by internal plasma currents) and is computed using virtual casing, 
!latex       and ${\bf B}_C$ is the ``vacuum'' field (produced by the external coils) and is given on input.
!latex \item The plasma field, ${\bf B}_P$, can only be computed after the equilibrium is determined, 
!latex       but this information is required to compute the equilibrium to begin with; and so there is an iteration involved.
!latex \item Suggested values of the vacuum field can be self generated; see \link{xspech} for more documentation on this.
!latex \end{enumerate}

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

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{compute the normal field on a regular grid on the computational boundary}
!latex \begin{enumerate}
!latex \item For each point on the compuational boundary, \link{casing} is called to compute the normal field produced by the plasma currents.
!latex \item (There is a very clumsy attempt to parallelize this which could be greatly improved.)
!latex \item An FFT gives the required Fourier harmonics.
!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

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

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!l tex \subsection{theory and numerics}

!l tex \begin{enumerate}
!l tex \item Required inputs to this subroutine are the geometry of the plasma boundary, 
!l tex       \be {\bf x}(\t,\z) \equiv x(\t,\z) {\bf i} + y(\t,\z) {\bf j} + z(\t,\z) {\bf k},
!l tex       \ee
!l tex       and the tangential field on this boundary, 
!l tex       \be {\bf B}_s=B^\t {\bf e}_\t + B^\z {\bf e}_\z,
!l tex       \ee
!l tex       where $\t$ and $\z$ are arbitrary poloidal and toroidal angles,
!l tex       and ${\bf e}_\t \equiv \partial {\bf x}/\partial \t$, ${\bf e}_\z \equiv \partial {\bf x}/\partial \z$.
!l tex       This routine assumes that the plasma boundary is a flux surface, i.e. ${\bf B} \cdot {\bf e}_\t \times {\bf e}_\z = 0$.
!l tex \item The virtual casing principle 
!l tex       [\paper{Shafranov \& Zakharov}{V.D. Shafranov \& L.E. Zakharov}{10.1088/0029-5515/12/5/009}{Nucl. Fusion}{12}{599}{1972},
!l tex        \paper{Lazerson}{S.A. Lazerson}{10.1088/0741-3335/54/12/122002}{Plasma Phys. Control. Fusion}{54}{122002}{2012},
!l tex        \paper{Hanson}{J.D. Hanson}{10.1088/0741-3335/57/11/115006}{Plasma Phys. Control. Fusion}{57}{115006}{2015}]
!l tex       shows that the field outside/inside the plasma arising from plasma currents inside/outside the boundary
!l tex       is equivalent to the field generated by a surface current,
!l tex       \be {\bf j} = {\bf B}_s \times {\bf n},
!l tex       \ee
!l tex       where ${\bf n}$ is normal to the surface.
!l tex \item The field at some arbitrary point, $\bar {\bf x}$, created by this surface current is given by
!l tex       \be {\bf B}({\bf \bar x}) = \int_{\cal S} \frac{\left( {\bf B}_s \times d{\bf s} \right) \times {\bf \hat r}}{r^2},
!l tex       \ee
!l tex       where $d{\bf s} \equiv {\bf e}_\t \times {\bf e}_\z \; d\t d\z$.
!l tex \item For ease of notation introduce
!l tex       \be {\bf J} & \equiv & {\bf B}_s \times d{\bf s} = \alpha \; {\bf e}_\t - \beta \; {\bf e}_\z,
!l tex       \ee
!l tex       where $\alpha \equiv B_\z =  B^\t g_{\t\z} + B^\z g_{\z\z}$ and $\beta \equiv B_\t = B^\t g_{\t\t} + B^\z g_{\t\z}$,
!l tex \item We may write in Cartesian coordinates ${\bf J} = j_x \; {\bf i} + j_y \; {\bf j} + j_z \; {\bf k}$, where
!l tex       \be j_x     &    =   & \alpha \; x_\t - \beta \; x_\z \\
!l tex           j_y     &    =   & \alpha \; y_\t - \beta \; y_\z \\
!l tex           j_z     &    =   & \alpha \; z_\t - \beta \; z_\z .
!l tex       \ee
!l tex \item Requiring that the current,
!l tex       \be {\bf j} & \equiv & \nabla \times {\bf B}
!l tex                       =      {\sqrt g}^{-1}(\partial_\t B_\z-\partial_\z B_\t) \; {\bf e}_\s
!l tex                       +      {\sqrt g}^{-1}(\partial_\z B_\s-\partial_\s B_\z) \; {\bf e}_\t
!l tex                       +      {\sqrt g}^{-1}(\partial_\s B_\t-\partial_\t B_\s) \; {\bf e}_\z,
!l tex       \ee
!l tex       has no normal component to the surface, i.e. ${\bf j}\cdot \nabla s=0$,
!l tex       we obtain the condition $\partial_\t B_\z = \partial_\z B_\t$, or $\partial_\t \alpha = \partial_\z \beta$.
!l tex       In axisymmetric configurations, where $\partial_\z \beta=0$, we must have $\partial_\t \alpha=0$.
!l tex \item The displacement from an arbitrary point, $(X,Y,Z)$, to a point, $(x,y,z)$, that lies on the surface is given 
!l tex       \be {\bf r} \equiv r_x \; {\bf i} + r_y \; {\bf j} + r_z \; {\bf k} = (X-x) \; {\bf i} + (Y-y) \; {\bf j} + (Z-z) \; {\bf k}.
!l tex       \ee
!l tex \item The components of the magnetic field produced by the surface current are then
!l tex       \be B^x &=&  \ooint (j_y r_z - j_z r_y)/r^3,\\
!l tex           B^y &=&  \ooint (j_z r_x - j_x r_z)/r^3,\\
!l tex           B^z &=&  \ooint (j_x r_y - j_y r_x)/r^3
!l tex       \ee
!l tex \item The surface integral is performed using \nag{www.nag.co.uk/numeric/FL/manual19/pdf/D01/d01eaf_fl19.pdf}{D01EAF}, which 
!l tex       uses an adaptive subdivision strategy and also computes absolute error estimates.
!l tex       The absolute and relative accuracy required are provided by the input \inputvar{vcasingtol}.
!l tex       The minimum number of function evaluations is provided by the input \inputvar{vcasingits}.
!l tex \item It may be convenient to have the derivatives:
!l tex       \be \frac{\partial B^x}{\partial x} & = & \ooint \left[ -3 (j_y r_z - j_z r_y) (X-x) / r^5 \;\;\;\;\;\;\;\;\;\;\;\;     \right], \\
!l tex           \frac{\partial B^x}{\partial y} & = & \ooint \left[ -3 (j_y r_z - j_z r_y) (Y-y) / r^5 - j_z/r^3                    \right], \\
!l tex           \frac{\partial B^x}{\partial z} & = & \ooint \left[ -3 (j_y r_z - j_z r_y) (Z-z) / r^5 + j_y/r^3                    \right], \\
!l tex           \frac{\partial B^y}{\partial x} & = & \ooint \left[ -3 (j_z r_x - j_x r_z) (X-x) / r^5 + j_z/r^3                    \right], \\
!l tex           \frac{\partial B^y}{\partial y} & = & \ooint \left[ -3 (j_z r_x - j_x r_z) (Y-y) / r^5 \;\;\;\;\;\;\;\;\;\;\;\;\;   \right], \\
!l tex           \frac{\partial B^y}{\partial z} & = & \ooint \left[ -3 (j_z r_x - j_x r_z) (Z-z) / r^5 - j_x/r^3                    \right], \\
!l tex           \frac{\partial B^z}{\partial x} & = & \ooint \left[ -3 (j_x r_y - j_y r_x) (X-x) / r^5 - j_y/r^3                    \right], \\
!l tex           \frac{\partial B^z}{\partial y} & = & \ooint \left[ -3 (j_x r_y - j_y r_x) (Y-y) / r^5 + j_x/r^3                    \right], \\
!l tex           \frac{\partial B^z}{\partial z} & = & \ooint \left[ -3 (j_x r_y - j_y r_x) (Z-z) / r^5 \;\;\;\;\;\;\;\;\;\;\;\;\;   \right].
!l tex       \ee
!l tex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!l tex \subsection{calculation of integrand}

!l tex \begin{enumerate}
!l tex \item Presently, the NAG routine D01DAF is used to compute the virtual casing integral, and this uses an adaptive integration. 
!l tex       Consequently, the magnetic field tangential to the plasma boundary is required at an arbitrary point.
!l tex       This is computed, as always, from ${\bf B} = \nabla \times {\bf A}$, and this provides ${\bf B} = B^\t {\bf e}_\t + B^\z {\bf e}_\z$.
!l tex       (Recall that $B^s=0$ by construction on the plasma boundary.)
!l tex       (It would be MUCH faster to only require the tangential field on a regular grid!!!)
!l tex \item Then, the metric elements $g_{\t\t}$, $g_{\t\z}$ and $g_{\z\z}$ are computed.
!l tex       These are used to ``lower'' the components of the magnetic field, ${\bf B} = B_\t \nabla \t + B_\z \nabla \z$.
!l tex       (Please check why $B_s$ is not computed. Is it because $B_s \nabla s \times {\bf n} = 0$ ?)
!l tex \item The distance between the ``evaluate'' point, $(X,Y,Z)$, and the given point on the surface, $(x,y,z)$ is computed.
!l tex \item If the computational boundary becomes too close to the plasma boundary, the distance is small and this causes problems for the numerics.
!l tex       I have tried to regularize this problem by introducing $\epsilon \equiv $\inputvar{vcasingeps}.
!l tex       Let the ``distance'' be
!l tex       \be D \equiv \sqrt{(X-x)^2 + (Y-y)^2 + (Z-Z)^2} + \epsilon^2.
!l tex       \ee
!l tex \item On taking the limit that $\epsilon \rightarrow 0$, the virtual casing integrand is 
!l tex        \be \verb+vcintegrand+ \equiv ( B_x n_x + B_y n_y + B_z n_z ) ( 1 + 3 \epsilon^2 / D^2 ) / D^3,
!l tex        \ee
!l tex       where the normal vector is ${\bf n} \equiv n_x {\bf i} + n_y {\bf j} + n_z {\bf k}$.
!l tex       This needs to be revised.
!l tex \end{enumerate}

REAL function vcintegrand( lteta, lzeta ) ! THIS ROUTINE IS NOT USED; differential virtual-casing field; format is fixed by NAG requirements;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half, one, three, pi
  
  use numerical, only :
  
  use fileunits, only : ounit, vunit
  
  use inputlist, only : Nvol, Igeometry, Lrad, vcasingeps
  
  use cputiming, only : 

  use allglobal, only : myid, ncpu, cpus, &
                        Mvol, &
                        mn, im, in, &
                        iRbc, iZbs, iRbs, iZbc, &
                        Ate, Aze, Ato, Azo, &
                        TT, &
                        YESstellsym, NOTstellsym, &
                        globaljk, Dxyz, Nxyz, &
                        tetazeta
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  REAL                :: lteta, lzeta, teta, zeta
  INTEGER             :: ii, mi, ni, ll, ideriv, jk
  REAL                :: dR(0:3), dZ(0:3), gBut, gBuz, gtt, gtz, gzz, sqrtg, Blt, Blz, Bxyz(1:3)
  REAL                :: arg, carg, sarg, czeta, szeta, XX, YY, ZZ, XXt, XXz, YYt, YYz, ZZt, ZZz, jj(1:3), rr(1:3), distance(1:3), firstorderfactor

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  teta = tetazeta(1) + lteta ; zeta = tetazeta(2) + lzeta
 !teta =               lteta ; zeta =               lzeta
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  dR(0:3) = zero ; ideriv = 0 ; gBut = zero ; gBuz = zero ! initialize summation of coordinates and tangential field;
  dZ(0:3) = zero
  
  jk = globaljk

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  select case( Igeometry )
   
  case( 1 ) ! Igeometry = 1 ; 09 Mar 17;
   
   if( YESstellsym ) then
    
    do ii = 1, mn ; mi = im(ii) ; ni = in(ii) ! loop over Fourier modes; construct surface current; slow transform required as position is arbitrary;
     
     arg = mi * teta - ni * zeta ; carg = cos(arg) ; sarg = sin(arg)
     
     dR(0) = dR(0) +          (                   iRbc(ii,Nvol) ) * carg                                            
     dR(1) = dR(1) +      (   (   iRbc(ii,Mvol) - iRbc(ii,Nvol) ) * carg                                            ) * half
     dR(2) = dR(2) + mi * ( - (                   iRbc(ii,Nvol) ) * sarg                                            )
     dR(3) = dR(3) - ni * ( - (                   iRbc(ii,Nvol) ) * sarg                                            )

    !dZ(0) = dZ(0) +                                                       (                 iZbs(ii,Nvol) ) * sarg 
    !dZ(1) = dZ(1) +      (                                                ( iZbs(ii,Mvol) - iZbs(ii,Nvol) ) * sarg ) * half
    !dZ(2) = dZ(2) + mi * (                                                (                 iZbs(ii,Nvol) ) * carg )
    !dZ(3) = dZ(3) - ni * (                                                (                 iZbs(ii,Nvol) ) * carg )
    
     do ll = 0, Lrad(Mvol)
      gBut = gBut - ( Aze(Mvol,ideriv,ii)%s(ll) * carg                                    ) * TT(ll,0,1) ! contravariant; Jacobian comes later; 
      gBuz = gBuz + ( Ate(Mvol,ideriv,ii)%s(ll) * carg                                    ) * TT(ll,0,1)
     enddo
     
    enddo ! end of do ii = 1, mn ;
    
   else ! NOTstellsym ; 08 Feb 16;
    
    do ii = 1, mn ; mi = im(ii) ; ni = in(ii) ! loop over Fourier modes; construct surface current; slow transform required as position is arbitrary;
     
     arg = mi * teta - ni * zeta ; carg = cos(arg) ; sarg = sin(arg)
     
     dR(0) = dR(0) +                              iRbc(ii,Nvol)   * carg +                   iRbs(ii,Nvol)   * sarg 
     dR(1) = dR(1) +      (   (   iRbc(ii,Mvol) - iRbc(ii,Nvol) ) * carg + ( iRbs(ii,Mvol) - iRbs(ii,Nvol) ) * sarg ) * half
     dR(2) = dR(2) + mi * ( -                     iRbc(ii,Nvol)   * sarg +                   iRbs(ii,Nvol)   * carg )
     dR(3) = dR(3) - ni * ( -                     iRbc(ii,Nvol)   * sarg +                   iRbs(ii,Nvol)   * carg )

    !dZ(0) = dZ(0) +                              iZbc(ii,Nvol)   * carg +                   iZbs(ii,Nvol)   * sarg 
    !dZ(1) = dZ(1) +      (   (   iZbc(ii,Mvol) - iZbc(ii,Nvol) ) * carg + ( iZbs(ii,Mvol) - iZbs(ii,Nvol) ) * sarg ) * half
    !dZ(2) = dZ(2) + mi * ( -                     iZbc(ii,Nvol)   * sarg +                   iZbs(ii,Nvol)   * carg )
    !dZ(3) = dZ(3) - ni * ( -                     iZbc(ii,Nvol)   * sarg +                   iZbs(ii,Nvol)   * carg )
    
     do ll = 0, Lrad(Mvol)
      gBut = gBut - ( Aze(Mvol,ideriv,ii)%s(ll) * carg + Azo(Mvol,ideriv,ii)%s(ll) * sarg ) * TT(ll,0,1) ! contravariant; Jacobian comes later; 
      gBuz = gBuz + ( Ate(Mvol,ideriv,ii)%s(ll) * carg + Ato(Mvol,ideriv,ii)%s(ll) * sarg ) * TT(ll,0,1)
     enddo
     
    enddo ! end of do ii = 1, mn ;
    
   endif ! end of if( YESstellsym ) ; 08 Feb 16;
   
  case( 2 ) ! Igeometry = 2 ; 09 Mar 17;
   
   FATAL( casing, .true., virtual casing under construction for cylindrical geometry )
   
  case( 3 ) ! Igeometry = 3 ; 09 Mar 17;
   
   if( YESstellsym ) then
    
    do ii = 1, mn ; mi = im(ii) ; ni = in(ii) ! loop over Fourier modes; construct surface current; slow transform required as position is arbitrary;
     
     arg = mi * teta - ni * zeta ; carg = cos(arg) ; sarg = sin(arg)
     
     dR(0) = dR(0) +          (                   iRbc(ii,Nvol) ) * carg                                            
     dR(1) = dR(1) +      (   (   iRbc(ii,Mvol) - iRbc(ii,Nvol) ) * carg                                            ) * half
     dR(2) = dR(2) + mi * ( - (                   iRbc(ii,Nvol) ) * sarg                                            )
     dR(3) = dR(3) - ni * ( - (                   iRbc(ii,Nvol) ) * sarg                                            )
     
     dZ(0) = dZ(0) +                                                       (                 iZbs(ii,Nvol) ) * sarg 
     dZ(1) = dZ(1) +      (                                                ( iZbs(ii,Mvol) - iZbs(ii,Nvol) ) * sarg ) * half
     dZ(2) = dZ(2) + mi * (                                                (                 iZbs(ii,Nvol) ) * carg )
     dZ(3) = dZ(3) - ni * (                                                (                 iZbs(ii,Nvol) ) * carg )
     
     do ll = 0, Lrad(Mvol)
      gBut = gBut - ( Aze(Mvol,ideriv,ii)%s(ll) * carg                                    ) * TT(ll,0,1) ! contravariant; Jacobian comes later; 
      gBuz = gBuz + ( Ate(Mvol,ideriv,ii)%s(ll) * carg                                    ) * TT(ll,0,1)
     enddo
     
    enddo ! end of do ii = 1, mn ;
    
   else ! NOTstellsym ; 08 Feb 16;
    
    do ii = 1, mn ; mi = im(ii) ; ni = in(ii) ! loop over Fourier modes; construct surface current; slow transform required as position is arbitrary;
     
     arg = mi * teta - ni * zeta ; carg = cos(arg) ; sarg = sin(arg)
     
     dR(0) = dR(0) +                              iRbc(ii,Nvol)   * carg +                   iRbs(ii,Nvol)   * sarg 
     dR(1) = dR(1) +      (   (   iRbc(ii,Mvol) - iRbc(ii,Nvol) ) * carg + ( iRbs(ii,Mvol) - iRbs(ii,Nvol) ) * sarg ) * half
     dR(2) = dR(2) + mi * ( -                     iRbc(ii,Nvol)   * sarg +                   iRbs(ii,Nvol)   * carg )
     dR(3) = dR(3) - ni * ( -                     iRbc(ii,Nvol)   * sarg +                   iRbs(ii,Nvol)   * carg )
     
     dZ(0) = dZ(0) +                              iZbc(ii,Nvol)   * carg +                   iZbs(ii,Nvol)   * sarg 
     dZ(1) = dZ(1) +      (   (   iZbc(ii,Mvol) - iZbc(ii,Nvol) ) * carg + ( iZbs(ii,Mvol) - iZbs(ii,Nvol) ) * sarg ) * half
     dZ(2) = dZ(2) + mi * ( -                     iZbc(ii,Nvol)   * sarg +                   iZbs(ii,Nvol)   * carg )
     dZ(3) = dZ(3) - ni * ( -                     iZbc(ii,Nvol)   * sarg +                   iZbs(ii,Nvol)   * carg )
     
     do ll = 0, Lrad(Mvol)
      gBut = gBut - ( Aze(Mvol,ideriv,ii)%s(ll) * carg + Azo(Mvol,ideriv,ii)%s(ll) * sarg ) * TT(ll,0,1) ! contravariant; Jacobian comes later; 
      gBuz = gBuz + ( Ate(Mvol,ideriv,ii)%s(ll) * carg + Ato(Mvol,ideriv,ii)%s(ll) * sarg ) * TT(ll,0,1)
     enddo
     
    enddo ! end of do ii = 1, mn ;
    
   endif ! end of if( YESstellsym ) ; 08 Feb 16;
  
  end select ! end of select case( Igeometry ) ; 09 Mar 17;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  select case( Igeometry )
   
  case( 1 ) ! Igeometry = 1 ; 09 Mar 17;
   
   gtt = one + dR(2)*dR(2)
   gtz =       dR(2)*dR(3)
   gzz = one + dR(3)*dR(3)
   
   sqrtg = dR(1)
   
  case( 2 ) ! Igeometry = 2 ; 09 Mar 17;
   
   FATAL( casing, .true., virtual casing under construction for cylindrical geometry )
   
  case( 3 ) ! Igeometry = 3 ; 09 Mar 17;
   
   gtt = dR(2)*dR(2) + dZ(2)*dZ(2)
   gtz = dR(2)*dR(3) + dZ(2)*dZ(3)
   gzz = dR(3)*dR(3) + dZ(3)*dZ(3) + dR(0)*dR(0)
   
   sqrtg = dR(0) * ( dZ(1) * dR(2) - dR(1) * dZ(2) )
   
  end select
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  Blt = ( gBut * gtt + gBuz * gtz ) / sqrtg
  Blz = ( gBut * gtz + gBuz * gzz ) / sqrtg
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  select case( Igeometry )
   
  case( 1 ) ! Igeometry = 1 ; 09 Mar 17;
   
   XX =          teta ; XXt =           one ; XXz =          zero
   YY =          zeta ; YYt =          zero ; YYz =           one
   ZZ = dR(0)         ; ZZt = dR(2)         ; ZZz = dR(3)
   
  case( 2 ) ! Igeometry = 2 ; 09 Mar 17;
   
   FATAL( casing, .true., virtual casing under construction for cylindrical geometry )
   
  case( 3 ) ! Igeometry = 3 ; toroidal geometry;
   
   czeta = cos( zeta ) ; szeta = sin( zeta )

   XX = dR(0) * czeta ; XXt = dR(2) * czeta ; XXz = dR(3) * czeta - dR(0) * szeta ! 10 Apr 13;
   YY = dR(0) * szeta ; YYt = dR(2) * szeta ; YYz = dR(3) * szeta + dR(0) * czeta   
   ZZ = dZ(0)         ; ZZt = dZ(2)         ; ZZz = dZ(3)
   
  end select

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  rr(1:3) = (/ Dxyz(1,jk) - XX, &
               Dxyz(2,jk) - YY, &
               Dxyz(3,jk) - ZZ /)

  jj(1:3) = (/ Blz * XXt - Blt * XXz, &
               Blz * YYt - Blt * YYz, &
               Blz * ZZt - Blt * ZZz /)

  distance(2) = sum( rr(1:3) * rr(1:3) ) + vcasingeps**2 ! 04 May 17;

  distance(1) = sqrt( distance(2) ) ; distance(3) = distance(1) * distance(2) ! powers of distance; 24 Nov 16;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  firstorderfactor = ( one + three * vcasingeps**2 / distance(2) ) / distance(3) ! 04 May 17;

  Bxyz(1:3) = (/ jj(2) * rr(3) - jj(3) * rr(2), &
                 jj(3) * rr(1) - jj(1) * rr(3), &
                 jj(1) * rr(2) - jj(2) * rr(1)  /) ! 04 May 17;

  vcintegrand = sum( Bxyz(1:3) * Nxyz(1:3,jk) ) * firstorderfactor ! 04 May 17;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!   FF( 4) = - three * FF(1) * rr(1) / distance(2)         ! dBxdx; ! NOTE THAT DERIVATIVES NEED TO BE DIVIDED BY DISTANCE(3) ; 14 Apr 17;
!   FF( 5) = - three * FF(1) * rr(2) / distance(2) - jj(3) ! dBxdy; ! NOTE THAT DERIVATIVES NEED TO BE DIVIDED BY DISTANCE(3) ; 14 Apr 17;
!   FF( 6) = - three * FF(1) * rr(3) / distance(2) + jj(2) ! dBxdz; ! NOTE THAT DERIVATIVES NEED TO BE DIVIDED BY DISTANCE(3) ; 14 Apr 17;
!   
!   FF( 7) = - three * FF(2) * rr(1) / distance(2) + jj(3) ! dBydx; ! NOTE THAT DERIVATIVES NEED TO BE DIVIDED BY DISTANCE(3) ; 14 Apr 17;
!   FF( 8) = - three * FF(2) * rr(2) / distance(2)         ! dBydy; ! NOTE THAT DERIVATIVES NEED TO BE DIVIDED BY DISTANCE(3) ; 14 Apr 17;
!   FF( 9) = - three * FF(2) * rr(3) / distance(2) - jj(1) ! dBydz; ! NOTE THAT DERIVATIVES NEED TO BE DIVIDED BY DISTANCE(3) ; 14 Apr 17;
!   
!   FF(10) = - three * FF(3) * rr(1) / distance(2) - jj(2) ! dBzdx; ! NOTE THAT DERIVATIVES NEED TO BE DIVIDED BY DISTANCE(3) ; 14 Apr 17;
!   FF(11) = - three * FF(3) * rr(2) / distance(2) + jj(1) ! dBzdy; ! NOTE THAT DERIVATIVES NEED TO BE DIVIDED BY DISTANCE(3) ; 14 Apr 17;
!   FF(12) = - three * FF(3) * rr(3) / distance(2)         ! dBzdz; ! NOTE THAT DERIVATIVES NEED TO BE DIVIDED BY DISTANCE(3) ; 14 Apr 17;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  return
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end function vcintegrand

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

REAL function zetalow( teta ) ! THIS ROUTINE IS NOT USED;
  
  use constants, only : pi
  
  use allglobal, only : tetazeta
  
  LOCALS
  
  REAL :: teta
  
  zetalow = tetazeta(2) - pi
  
  return
  
end function zetalow

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

REAL function zetaupp( teta ) ! THIS ROUTINE IS NOT USED;
  
  use constants, only : pi
  
  use allglobal, only : tetazeta
  
  LOCALS
  
  REAL :: teta
  
  zetaupp = tetazeta(2) + pi
  
  return
  
end function zetaupp

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
