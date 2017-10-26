!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Given vector position, returns force and derivative matrix.

!latex \item The force vector, ${\bf F}({\bf x})$, is a combination of the pressure-imbalance Fourier harmonics, $[[p+B^2/2]]_{m,n}$,
!latex       and the spectral condensation constraints, $\partial M/\partial \lambda_j$.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine fc02aa( Ngeometricaldof, position, force, LComputeDerivatives )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half, one, two, pi2, pi
  
  use numerical, only : vsmall, small, logtolerance
  
  use fileunits, only : ounit, wunit
  
  use inputlist, only : Wmacros, Wfc02aa, ext, Nvol, Mpol, Ntor, Lrad, tflux, Igeometry, &
                        gamma, adiabatic, pscale, mu, &
                        Energy, ForceErr, &
                        epsilon, opsilon, &
                        Lminimize, Lfindzero, &
                        Lconstraint, Lcheck, &
                        Lextrap
  
  use cputiming, only : Tfc02aa
  
  use allglobal, only : ncpu, myid, cpus, pi2nfp, &
                        Lcoordinatesingularity, Lvacuumregion, Lplasmaregion, Mvol, &
                        Iquad, &                  ! convenience; provided to ma00aa as argument to avoid allocations;
                        iRbc, iZbs, iRbs, iZbc, & ! Fourier harmonics of geometry; vector of independent variables, position, is "unpacked" into iRbc,iZbs;
                        Nmagneticdof, &
                        ImagneticOK, &
                        YESstellsym, NOTstellsym, &
                        mn, im, in, mns, halfmm, Ntz, &
                        ijreal, ijimag, jireal, jiimag, &
                        efmn, ofmn, cfmn, sfmn, &
                        evmn, odmn, comn, simn, &
                        trigm, trign, trigwk, isr, Nt, Nz, &
                        cosi, sini, & ! FFT workspace;
                        DifferentiateGeometry, &
                        dMA, dMB, dMC, dMD, dME, dMF, solution, &
                        MBpsi, MEpsi, &
                        dtflux, dpflux, sweight, &
                        mpnq, &
                        Bemn, Bomn, Iomn, Iemn, Somn, Semn, &
!                       Pomn, Pemn, & ! computed in bb00aa; 20 Jun 14;
                        BBe, IIo, BBo, IIe, &
                        lgeometricaldof, &
                        dBBdRZ, dIIdRZ, &
                        vvolume, dvolume, lBBintegral, lABintegral, &
                        Rij, Zij, sg, guvij, iRij, iZij, dRij, dZij, tRij, tZij, & ! Jacobian and metrics; computed in co01aa;
                        diota, &
                        dFFdRZ, dBBdmp, dmupfdx, hessian, dessian,Lhessianallocated, &
                        expmmnn, & ! exponential weight on force-imbalance harmonics; 04 Dec 14;
                        psifactor, &
                        Llatex, &
                        lmns
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in)  :: Ngeometricaldof               ! dimensions;
  REAL,    intent(in)  :: position(0:Ngeometricaldof)   ! degrees-of-freedom = internal geometry;
  REAL,    intent(out) :: force(0:Ngeometricaldof)      ! force;
  LOGICAL, intent(in)  :: LComputeDerivatives           ! indicates whether derivatives are to be calculated;
  
  INTEGER              :: Ndof, IA, ifail, if01adf, vflag
  
  INTEGER              :: vvol, innout, ii, jj, lnn, irz, issym, iocons, tdoc, idoc, idof, tdof, jdof, ivol, imn

  INTEGER              :: Lcurvature, ideriv, id
  
  REAL                 :: lastcpu, lss, lfactor, dRdt(-1:1,0:1), dZdt(-1:1,0:1), DDl, MMl

  INTEGER              :: iflag
  REAL                 :: Itor, Gpol, det
  
  REAL                 :: dpsi(1:2)
  REAL   , allocatable :: oBI(:,:), rhs(:)

  REAL   , allocatable :: dAt(:,:), dAz(:,:), XX(:), YY(:), dBB(:,:), dII(:), dLL(:), dPP(:), length(:), dRR(:,:), dZZ(:,:), constraint(:)

  CHARACTER            :: packorunpack 

#ifdef DEBUG
  INTEGER              :: isymdiff
  REAL                 :: dRZ = 1.0e-05, dvol(-1:+1), evolume, imupf(1:2,-2:2)
  REAL,    allocatable :: oRbc(:,:), oZbs(:,:), oRbs(:,:), oZbc(:,:) ! original geometry; 20 Jun 14;
  REAL,    allocatable :: isolution(:,:)
#endif

  BEGIN(fc02aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  packorunpack = 'U' ! unpack geometrical degrees-of-freedom; 13 Sep 13;
  WCALL(fc02aa,gf00aa,( Ngeometricaldof, position(0:Ngeometricaldof), Mvol, mn, &
iRbc(1:mn,0:Mvol), iZbs(1:mn,0:Mvol), iRbs(1:mn,0:Mvol), iZbc(1:mn,0:Mvol), packorunpack ))
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef DEBUG
  
  if( Lcheck.eq.3 .or. Lcheck.eq.4 ) then ! will check volume derivatives; 01 Jul 14;
   
   RALLOCATE(oRbc,(1:mn,0:Mvol))
   RALLOCATE(oZbs,(1:mn,0:Mvol))
   RALLOCATE(oRbs,(1:mn,0:Mvol))
   RALLOCATE(oZbc,(1:mn,0:Mvol))  
   
   oRbc(1:mn,0:Mvol) = iRbc(1:mn,0:Mvol)
   oZbs(1:mn,0:Mvol) = iZbs(1:mn,0:Mvol)
   oRbs(1:mn,0:Mvol) = iRbs(1:mn,0:Mvol)
   oZbc(1:mn,0:Mvol) = iZbc(1:mn,0:Mvol)
   
  endif ! end of if( Lcheck.eq.3 ) ; 01 Jul 14;

#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RALLOCATE(dAt,(1:Ntz,-1:2))
  RALLOCATE(dAz,(1:Ntz,-1:2))

  RALLOCATE(dBB,(1:Ntz,-1:2)) ! magnetic field strength (on interfaces) in real space and derivatives; 01 Jul 14;

  RALLOCATE( XX,(1:Ntz))
  RALLOCATE( YY,(1:Ntz))

  RALLOCATE( length,(1:Ntz)) ! this is calculated in bb00aa; 04 Dec 14;

  RALLOCATE( dRR,(1:Ntz,-1:1))
  RALLOCATE( dZZ,(1:Ntz,-1:1))

  RALLOCATE(dII,(1:Ntz)) ! spectral constraint; 08 Nov 13;
  RALLOCATE(dLL,(1:Ntz)) ! length   constraint; 08 Nov 13;
  RALLOCATE(dPP,(1:Ntz)) ! poloidal constraint; 08 Nov 13;

  RALLOCATE(constraint,(1:Ntz))
 
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef DEBUG  
  FATALMESS(fc02aa, Mvol.lt.1, illegal allocation)
#endif

  ImagneticOK(1:Mvol) = .false. ! default is to assume that the construction of the magnetic field in each volume has failed; 22 Apr 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( LcomputeDerivatives ) then
#ifdef DEBUG
   FATALMESS( fc02aa, .not.allocated(dBBdmp), do not pass go)
#endif
   dBBdmp(1:lgeometricaldof,1:Mvol,0:1,1:2) = zero
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  do vvol = 1, Mvol ! loop over volumes;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   if( myid.ne.modulo(vvol-1,ncpu) ) cycle ! construct Beltrami fields in parallel; 16 Jan 13;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
  !if( vvol.eq.1 .and. Igeometry.gt.1 ) then ; Lcoordinatesingularity = .true.  ! 04 Dec 14;
  !else                                      ; Lcoordinatesingularity = .false. ! 04 Dec 14;
  !endif

  !if( Igeometry.eq.1 .or.  vvol.gt.1 ) Lcoordinatesingularity = .false. ! 04 Dec 14;
  !if( Igeometry.gt.1 .and. vvol.eq.1 ) Lcoordinatesingularity = .true.  ! 04 Dec 14;
   
   if( Igeometry.eq.1 .or. vvol.gt.1 ) then ; Lcoordinatesingularity = .false. ! 14 Jan 15;
   else                                     ; Lcoordinatesingularity = .true.  ! 14 Jan 15;
   endif
   
   if( vvol.le.Nvol ) then ; Lplasmaregion = .true.  ; Lvacuumregion = .not.Lplasmaregion
   else                    ; Lplasmaregion = .false. ; Lvacuumregion = .not.Lplasmaregion
   endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   Ndof = Nmagneticdof(vvol) ! shorthand; ! Nmagneticdof was computed in al00aa; 16 Jan 13;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
! allocate matrices, etc. ! 16 Jan 13;
   
   RALLOCATE(dMA,(0:Ndof,0:Ndof)) ! required for both plasma region and vacuum region; 18 Apr 13;
   RALLOCATE(dMB,(0:Ndof,0:   2))
   
   RALLOCATE(solution,(1:Ndof,-1:2)) ! this will contain the vector/scalar potential from the linear solver and it's derivatives; 01 May 13;

   
   if( Lplasmaregion ) then ! required only for plasma region; 18 Apr 13;

    RALLOCATE(dMD,(0:Ndof,0:Ndof))
    RALLOCATE(dME,(0:Ndof,1:   2))
    RALLOCATE(dMC,(1:   2,1:   2))
    RALLOCATE(dMF,(1:   2,1:   2))
    
    RALLOCATE(MBpsi,(1:Ndof))
    RALLOCATE(MEpsi,(1:Ndof))

   endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   if( LcomputeDerivatives ) then ! allocate some additional memory; 18 Dec 14;
    
    RALLOCATE(oBI,(1:Ndof,1:Ndof)) ! inverse of ``original'', i.e. unperturbed, Beltrami matrix; 19 Sep 13;
    RALLOCATE(rhs,(1:Ndof))
    
#ifdef DEBUG
    if( Lcheck.eq.4 ) then
     RALLOCATE(isolution,(1:Ndof,-2:2))
    endif
#endif

   endif ! end of if( LcomputeDerivatives ) 19 Sep 13;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   WCALL(fc02aa,ma00ab,( 'A', vvol )) ! allocate volume integrated metric elements;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   DifferentiateGeometry%L = .false. ! first, compute Beltrami fields; 02 Sep 14;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   lastcpu = GETTIME ! 30 Jan 13;
   
   WCALL(fc02aa,ma00aa,( Iquad(vvol), mn, vvol, Lrad(vvol) )) ! compute volume integrals of metric elements;
   
   cput = GETTIME
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   lastcpu = GETTIME ! 30 Jan 13;
   
   if( Lplasmaregion ) then
    WCALL(fc02aa,ma01ag,( vvol, mn, Lrad(vvol) )) ! construct generalized Beltrami matrices; dMA, dMB & dMC; dMD, dME & dMF;
   else ! vacuum region; 13 Sep 13;
    WCALL(fc02aa,va00aa,(       mn, Lrad(vvol) )) ! construct             Laplace  matrices;
   endif
  
   cput = GETTIME
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
#ifdef DEBUG

   if( Llatex ) then

    write(wunit,'("\newpage \subsection{{energy matrices and helicity matrices : volume=$"i3"$; [fc02aa]}}")') vvol

    write(wunit,'("\be A = + \frac{{ 8 \pi^2 }}{{3}} \left( \begin{{array}}{{rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr}}")')
    do ii = 1, Ndof
     write(wunit,'(f12.3,999(:",&"f12.3))') dMA(ii,1:Ndof) / ( 8 * pi**2 / 3 )
     write(wunit,'("\\")') 
    enddo
    write(wunit,'("\end{{array}} \right) \ee")')

    write(wunit,'("\be B = + \frac{{ 8 \pi^2 }}{{3}} \left( \begin{{array}}{{rr}}")')
    do ii = 1, Ndof
     write(wunit,'(f12.3,999(:",&"f12.3))') dMB(ii,1:   2) / ( 8 * pi**2 / 3 )
     write(wunit,'("\\")') 
    enddo
    write(wunit,'("\end{{array}} \right) \ee")')

    write(wunit,'("\be C = + \frac{{ 8 \pi^2 }}{{3}} \left( \begin{{array}}{{rr}}")')
    do ii = 1,    2
     write(wunit,'(f12.3,999(:",&"f12.3))') dMC(ii,1:   2) / ( 8 * pi**2 / 3 )
     write(wunit,'("\\")') 
    enddo
    write(wunit,'("\end{{array}} \right) \ee")')

    write(wunit,'("\be D = + \frac{{ 8 \pi^2 }}{{3}} \left( \begin{{array}}{{rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr}}")')
    do ii = 1, Ndof
     write(wunit,'(f12.3,999(:",&"f12.3))') dMD(ii,1:Ndof) / ( 8 * pi**2 / 3 )
     write(wunit,'("\\")') 
    enddo
    write(wunit,'("\end{{array}} \right) \ee")')

    write(wunit,'("\be E = + \frac{{ 8 \pi^2 }}{{3}} \left( \begin{{array}}{{rr}}")')
    do ii = 1, Ndof
     write(wunit,'(f12.3,999(:",&"f12.3))') dME(ii,1:   2) / ( 8 * pi**2 / 3 )
     write(wunit,'("\\")') 
    enddo
    write(wunit,'("\end{{array}} \right) \ee")')

    write(wunit,'("\be F = - \frac{{20 \pi^2 }}{{3}} \left( \begin{{array}}{{rr}}")')
    do ii = 1,    2
     write(wunit,'(f12.3,999(:",&"f12.3))') dMF(ii,1:   2) / ( 8 * pi**2 / 3 )
     write(wunit,'("\\")') 
    enddo
    write(wunit,'("\end{{array}} \right) \ee")')

   endif ! 20 Jan 15;

#endif
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   lastcpu = GETTIME ! 30 Jan 13;
   
   WCALL(fc02aa,ma02aa,( vvol, Ndof )) ! calls mp00ac to compute the magnetic fields (perhaps iteratively to satisfy transform/helicity constraints);

   cput = GETTIME

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef DEBUG
   if( Lcheck.eq.2 ) then
    goto 2000 ! will take no other action except a finite-difference comparison on the derivatives of the rotational-transform wrt mu and dpflux; 20 Jun 14;
   endif
#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   if( LcomputeDerivatives ) then ! compute inverse of matrix; 13 Sep 13; ! compute inverse of Beltrami matrices; 29 Apr 14;
    
    lastcpu = GETTIME ! 30 Jan 13;

    dMA(0:Ndof-1,1:Ndof) = dMA(1:Ndof,1:Ndof) - mu(vvol) * dMD(1:Ndof,1:Ndof) ! this corrupts dMA, but dMA is no longer used;  7 Mar 13; 
    dMA(  Ndof  ,1:Ndof) = zero

    dMD(1:Ndof  ,1:Ndof) = dMA(0:Ndof-1,1:Ndof) ! copy of original matrix; this is used below; 13 Sep 13;
    
    IA = Ndof + 1
    
    if01adf = 1
    call F01ADF( Ndof, dMA(0:Ndof,1:Ndof), IA, if01adf ) ! assumes symmetric, positive definite matrix; 29 Jan 13; ! dMA is corrupted;  7 Mar 13; 
    
    cput = GETTIME
    
    select case( if01adf ) !                                                                    0123456789012345678
    case( 0 )    ; if( Wfc02aa ) write(ounit,1010) cput-cpus, myid, vvol,cput-lastcpu, if01adf, "success ;         "
    case( 1 )    ;               write(ounit,1010) cput-cpus, myid, vvol,cput-lastcpu, if01adf, "not +ve definite ;"
    case( 2 )    ;               write(ounit,1010) cput-cpus, myid, vvol,cput-lastcpu, if01adf, "input error ;     "
    case default ;               FATALMESS(ma01ga, .true., illegal ifail returned from F01ADF)
    end select
    
1010 format("fc02aa : ",f10.2," : myid=",i3," ; vvol=",i3," ; called F01ADF ; time=",f10.2,"s ; computed inverse of Beltrami matrix; if01adf=",i2," ; ",a18)
    
    do ii = 1, Ndof
     do jj = 1, ii
      oBI(ii,jj) = dMA(ii,jj) ! inverse matrix is returned in lower triangle; see NAG documentation for details; 29 Jan 13;
      oBI(jj,ii) = dMA(ii,jj) ! I assume the inverse is symmetric; 29 Jan 13;
     enddo
    enddo
    
    if( Lconstraint.eq.1 ) then ! first, determine how B^2 varies with mu and dpflux; 20 Jun 14;
     
     do iocons = 0, 1 ! labels constraint; 20 Jun 14;
      
      if( vvol.eq.   1 .and. iocons.eq.0 ) cycle ! fixed inner boundary (or coordinate axis); no force-balance constraints; 19 Sep 13;
      if( vvol.eq.Mvol .and. iocons.eq.1 ) cycle ! fixed outer boundary                     ; no force-balance constraints; 19 Sep 13;
      
      do ideriv = 0, 2 ; id = ideriv ! derivatives wrt helicity multiplier and differential poloidal flux; 02 Sep 14;

       iflag = 1 ! bb00aa will only return dAt(1:Ntz,id) and dAz(1:Ntz,id);
       
       WCALL(fc02aa, bb00aa,( vvol, iocons, ideriv, Ntz, dAt(1:Ntz,id), dAz(1:Ntz,id), XX(1:Ntz), YY(1:Ntz), length(1:Ntz), DDl, MMl, iflag ) )
       
       if( ideriv.eq.0 ) then
        dBB(1:Ntz,id) = half * ( &
     dAz(1:Ntz, 0)*dAz(1:Ntz, 0)*guvij(1:Ntz,2,2,0) - two * dAz(1:Ntz, 0)*dAt(1:Ntz, 0)*guvij(1:Ntz,2,3,0) + dAt(1:Ntz, 0)*dAt(1:Ntz, 0)*guvij(1:Ntz,3,3,0) &
     ) / sg(1:Ntz,0)**2
       else
        dBB(1:Ntz,id) = half * ( &
     dAz(1:Ntz,id)*dAz(1:Ntz, 0)*guvij(1:Ntz,2,2,0) - two * dAz(1:Ntz,id)*dAt(1:Ntz, 0)*guvij(1:Ntz,2,3,0) + dAt(1:Ntz,id)*dAt(1:Ntz, 0)*guvij(1:Ntz,3,3,0) &
   + dAz(1:Ntz, 0)*dAz(1:Ntz,id)*guvij(1:Ntz,2,2,0) - two * dAz(1:Ntz, 0)*dAt(1:Ntz,id)*guvij(1:Ntz,2,3,0) + dAt(1:Ntz, 0)*dAt(1:Ntz,id)*guvij(1:Ntz,3,3,0) &
     ) / sg(1:Ntz,0)**2
       endif ! end of if( ideriv.gt.0 ) ; 20 Jun 14;
      
      enddo ! end of do ideriv = 0, 2; 20 Jun 14;

      call tfft( Nt, Nz, dBB(1:Ntz,1), dBB(1:Ntz,2), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz), & ! derivatives of B^2 wrt mu and dpflux; 02 Sep 14;
                 mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), ifail )
         
      ; idoc = 0
      ; dBBdmp(idoc+1:idoc+mn  ,vvol,iocons,1) = efmn(1:mn) * expmmnn(1:mn) * opsilon ! pressure; 04 Dec 14;
      ; dBBdmp(idoc+1:idoc+mn  ,vvol,iocons,2) = cfmn(1:mn) * expmmnn(1:mn) * opsilon ! pressure; 04 Dec 14;
      ; idoc = idoc + mn   ! even; 19 Sep 13;
      if( Igeometry.ge.3 ) then ! add spectral constraints; 19 Sep 13; spectral constraints do not depend on mu or dpflux; 01 Jul 14;
       ;idoc = idoc + mn-1 ! oddd; 19 Sep 13;
      endif ! end of if( Igeometry.ge.3) ; 20 Jun 14;
      if( NOTstellsym ) then
       ;dBBdmp(idoc+1:idoc+mn-1,vvol,iocons,1) = ofmn(2:mn) * expmmnn(1:mn) * opsilon ! pressure; 04 Dec 14;
       ;dBBdmp(idoc+1:idoc+mn-1,vvol,iocons,2) = sfmn(2:mn) * expmmnn(1:mn) * opsilon ! pressure; 04 Dec 14;
       ;idoc = idoc + mn-1 ! oddd; 19 Sep 13;
       if( Igeometry.ge.3 ) then ! add spectral constraints; 19 Sep 13;
        idoc = idoc + mn   ! even; 19 Sep 13;
       endif ! end of if( Igeometry.ge.3) ; 20 Jun 14;
      endif ! end of if( NOTstellsym) ; 20 Jun 14;
      
     enddo ! end of do iocons; 20 Jun 14;
     
    endif ! end of if( Lconstraint.eq.1 ) ; 20 Jun 14;
    
   endif ! end of if( LcomputeDerivatives ) ;  7 Mar 13; 
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   vflag = 1 ! this flag instructs vo00aa to continue even if the volume is invalid; 04 Dec 14;
   WCALL(fc02aa, vo00aa,( vvol, vflag ) ) ! compute volume; 13 Sep 13;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   do iocons = 0, 1 ! construct field magnitude on inner and outer interfaces; inside do vvol;
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
    if( vvol.eq.1    .and. iocons.eq.0 ) cycle ! fixed inner boundary (or coordinate axis);
    if( vvol.eq.Mvol .and. iocons.eq.1 ) cycle ! fixed outer boundary                     ; there are no constraints at outer boundary;
    
    ideriv = 0 ; id = ideriv

    iflag = 0 ! dAt, dAz, XX & YY are returned by bb00aa; Bemn(1:mn,vvol,iocons), Iomn(1:mn,vvol) etc. are returned through global; 01 Jul 14;
    
    WCALL(fc02aa, bb00aa,( vvol, iocons, ideriv, Ntz, dAt(1:Ntz,id), dAz(1:Ntz,id), XX(1:Ntz), YY(1:Ntz), length(1:Ntz), DDl, MMl, iflag ) )
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
!latex \end{enumerate} \subsection{energy functional gradient} \begin{enumerate} 
    
!latex \item Consider 
!latex       \be \delta F = - \int_{\cal I} [[p+B^2/2]] \; \boldxi \cdot {\bf dS} = 
!latex       - \int \!\!\!\!\! \int [[p+B^2/2]] \; \boldxi \cdot {\bf e}_\t \times {\bf e}_\z \; d\t d\z.
!latex       \ee
    
!    if( Lminimize.gt.0 ) then 
!
!     select case( Igeometry )
!     
!    !case( 1 )
!     
!    !case( 2 )
!     
!     case( 3 )
!      
!      ijreal(1:Ntz) = ( adiabatic(vvol) * pscale / vvolume(vvol)**gamma + half * dBB(1:Ntz,0) ) * Rij(1:Ntz,0,0) * Zij(1:Ntz,2,0) ! sign factor below;
!      ijimag(1:Ntz) = ( adiabatic(vvol) * pscale / vvolume(vvol)**gamma + half * dBB(1:Ntz,0) ) * Rij(1:Ntz,0,0) * Rij(1:Ntz,2,0) ! sign factor below;
!      
!      call tfft( Nt, Nz, ijreal(1:Ntz), ijimag(1:Ntz), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz), & ! Fourier harmonics of force-imbalance;
!                 mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), ifail )
!
!
!      if( iocons.eq.1 ) then ! compute spectral constraint;  7 Mar 13; 
!       
!       jireal(1:Ntz) = Rij(1:Ntz,2,0) * dII(1:Ntz) ! weighted spectral constraints; 26 Feb 13; ! this is  odd * odd is even (for stell-sym); 26 Feb 13;
!       jiimag(1:Ntz) = Zij(1:Ntz,2,0) * dII(1:Ntz) ! weighted spectral constraints; 26 Feb 13; ! this is even * odd is odd  (for stell-sym); 26 Feb 13;
!       
!       call tfft( Nt, Nz, jireal(1:Ntz), jiimag(1:Ntz), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz), &
!                  mn, im(1:mn), in(1:mn), evmn(1:mn), odmn(1:mn), comn(1:mn), simn(1:mn), ifail )
!        
!      endif ! end of if( iocons.eq.1 ) ; 19 Sep 13;
!
!
!      idof = 0 ! local to interface;
!      
!      do ii = 1, mn
!       
!       lfactor = pi * pi2nfp * tflux(vvol-1+iocons)**regularmm(ii) * lss ! DEFINITION OF TFLUX MAY HAVE CHANGED; include factor of pi2; 01 Jul 14;
!       
!       do irz = 0, 1
!        
!        if( irz.eq.1 .and. Igeometry.ne.3 ) cycle ! no dependence on Z; 26 Feb 13;
!        
!        do issym = 0, 1 ! loop over stellarator and non-stellarator symmetric terms; 26 Feb 13;
!         
!         if( issym.eq.1 .and. YESstellsym ) cycle ! no dependence on non-stellarator symmetric harmonics; 26 Feb 13;
!         
!         FATALMESS(fc02aa, NOTstellsym, under construction)
!         
!         if( ii.eq.1 .and. irz.eq.1 .and. issym.eq.0 ) cycle ! no dependence on Z_{m,n} \sin( m\t-n\z) for m=0, n=0;
!         if( ii.eq.1 .and. irz.eq.0 .and. issym.eq.1 ) cycle ! no dependence on Z_{m,n} \sin( m\t-n\z) for m=0, n=0;
!         
!         idof = idof + 1
!         
!         if( irz.eq.0 .and. ii.eq.1 ) dBBdRZ(vvol,iocons,idof) =   efmn(ii) * two * lfactor ! R;
!         if( irz.eq.0 .and. ii.gt.1 ) dBBdRZ(vvol,iocons,idof) =   efmn(ii)       * lfactor ! R;
!        !if( irz.eq.1 .and. ii.eq.1 ) dBBdRZ(vvol,iocons,idof) = - sfmn(ii) * two * lfactor ! Z; this mode is always irrelevant;
!         if( irz.eq.1 .and. ii.gt.1 ) dBBdRZ(vvol,iocons,idof) = - sfmn(ii)       * lfactor ! Z;
!        
!         if( iocons.eq.1 ) then ! compute spectral constraint; 19 Sep 13;
!         if( irz.eq.0 .and. ii.eq.1 ) dIIdRZ(vvol       ,idof) =   evmn(ii) * two * lfactor ! R;
!         if( irz.eq.0 .and. ii.gt.1 ) dIIdRZ(vvol       ,idof) =   evmn(ii)       * lfactor ! R;
!        !if( irz.eq.1 .and. ii.eq.1 ) dIIdRZ(vvol       ,idof) = - simn(ii) * two * lfactor ! Z; this mode is always irrelevant;
!         if( irz.eq.1 .and. ii.gt.1 ) dIIdRZ(vvol       ,idof) = - simn(ii)       * lfactor ! Z;
!         endif
!
!        enddo ! endo of do issym; 26 Feb 13;
!        
!       enddo ! end of do irz;
!       
!      enddo ! end of do ii;
!      
!     case default ! 10 Mar 13;
!      
!      FATALMESS(fc02aa, .true., supplied Igeometry not yet supported)
!      
!     end select !end of select case( Igeometry ) ; 26 Feb 13;
!     
!    endif ! end of if( Lminimize ) ; 

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   enddo ! end of do iocons = 0, 1;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! still inside do vvol = 1, Mvol loop;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   if( LComputeDerivatives ) then

    DifferentiateGeometry%L = .true. ! will need derivatives;
    
    idof = 0 ! labels degree of freedom; local to interface;

    
    do ii = 1, mn ! loop over deformations in Fourier harmonics; inside do vvol;
     
     DifferentiateGeometry%ii = ii ! controls construction of derivatives in subroutines called below;
     

     do irz = 0, 1 ! loop over deformations in R and Z; inside do vvol; inside do ii;
      
      if( irz.eq.1 .and. Igeometry.lt.3 ) cycle ! no dependence on Z; 26 Feb 13;
      
      DifferentiateGeometry%irz = irz ! controls construction of derivatives;


      do issym = 0, 1 ! loop over stellarator and non-stellarator symmetric terms; 26 Feb 13;
       
       if( issym.eq.1 .and. YESstellsym ) cycle ! no dependence on non-stellarator symmetric harmonics; 26 Feb 13;
       
       if( ii.eq.1 .and. irz.eq.1 .and. issym.eq.0 ) cycle ! no dependence on Zbs_{m=0,n=0}; 26 Feb 13;
       if( ii.eq.1 .and. irz.eq.0 .and. issym.eq.1 ) cycle ! no dependence on Rbs_{m=0,n=0}; 26 Feb 13;
       
       DifferentiateGeometry%issym = issym ! controls construction of derivatives;

       
       idof = idof + 1 ! this labels the degree-of-freedom that the derivative is taken with respect to; 19 Sep 13; ! this is outside do innout; 19 Sep 13;
       
       
       do innout = 0, 1 ! loop over deformations to inner and outer interface; inside do vvol; inside do ii; inside do irz;
        
        if( vvol.eq.1    .and. innout.eq.0 ) cycle ! no degrees of freedom at coordinate axis / fixed inner boundary;
        
!      !if( vvol.eq.Mvol .and. innout.eq.1 ) cycle ! no degress of freedom                      fixed outer boundary; for linearized displacement; 18 Dec 14;
        
        DifferentiateGeometry%innout = innout

        WCALL(fc02aa,ma00aa,( Iquad(vvol), mn, vvol, Lrad(vvol) )) ! compute volume integrals of metric elements;

        if( Lplasmaregion ) then ! plasma region; 19 Sep 13;
         WCALL(fc02aa,ma01ag,( vvol, mn, Lrad(vvol) )) ! construct generalized          matrices;
        else                     ! vacuum region; 19 Sep 13;
         WCALL(fc02aa,va00aa,(       mn, Lrad(vvol) )) ! construct             Laplace  matrices;
        endif

        dpsi(1:2) = (/ dtflux(vvol), dpflux(vvol) /) ! local enclosed toroidal and poloidal fluxes; 18 Dec 14;
        
        rhs(1:Ndof) = - matmul( dMB(1:Ndof,1:2   ) - mu(vvol) * dME(1:Ndof,1:2   ), dpsi(1:2)          ) & ! 20 Jun 14;
                      - matmul( dMA(1:Ndof,1:Ndof) - mu(vvol) * dMD(1:Ndof,1:Ndof), solution(1:Ndof,0) )
        
        solution(1:Ndof,-1) = matmul( oBI(1:Ndof,1:Ndof), rhs(1:Ndof) ) ! this is the perturbed, packed solution; 01 Jul 14;

        ideriv = -1 ; dpsi(1:2) = (/ dtflux(vvol), dpflux(vvol) /) ! these are also used below; 26 Feb 13;
        
        packorunpack = 'U'
        
        WCALL(fc02aa,up00aa,( packorunpack, vvol, Ndof,  solution(1:Ndof,-1), dpsi(1:2), ideriv )) ! derivatives placed in Ate(vvol,ideriv,1:mn)%s(0:Lrad),

        lfactor = psifactor(ii,vvol-1+innout) ! this "pre-conditions" the geometrical degrees-of-freedom; 02 Sep 14;
        
        if( Lconstraint.eq.1 ) then ! need to determine how mu and dplux vary with geometry at fixed interface transform; 20 Jun 14;
         
         iflag = -1 ! instructs tr00ab to compute derivative wrt geometry; 20 Jun 14;

        !if( YESstellsym ) lmns = 1 + (mns-1)           ! number of independent degrees of freedom in angle transformation; 30 Jan 13; 
        !if( NOTstellsym ) lmns = 1 + (mns-1) + (mns-1) ! number of independent degrees of freedom in angle transformation; 30 Jan 13; 

         Itor = zero ; Gpol = zero
         
         FATALMESS(fc02aa, Lvacuumregion, have not set Itor and Gpol)
         
! given d(vector potential)/dx, compute d(transform)/dx; diota(ideriv,innout);
         WCALL(fc02aa,tr00ab,( vvol, mn, lmns, Nt, Nz, Itor, Gpol, iflag, diota(0:1,-1:2,vvol) ) ) 
         
         if( Lcoordinatesingularity ) then ! solution does not depend on dpflux, and only outer transform is a constraint; 01 Jul 14;
          det = diota(1,1,vvol)
          FATALMESS(fc02aa, abs(det).lt.small, error computing derivatives of mu          wrt geometry at fixed transform)
          dmupfdx(vvol,1,idof,innout) = - lfactor * (                                                          diota(1,-1,vvol) ) / det
          dmupfdx(vvol,2,idof,innout) =   zero
         else
          det = diota(0,1,vvol) * diota(1,2,vvol) - diota(0,2,vvol) * diota(1,1,vvol)
          FATALMESS(fc02aa, abs(det).lt.small, error computing derivatives of mu & dpflux wrt geometry at fixed transform)
          dmupfdx(vvol,1,idof,innout) = - lfactor * ( + diota(1, 2,vvol) * diota(0,-1,vvol) - diota(0, 2,vvol) * diota(1,-1,vvol) ) / det
          dmupfdx(vvol,2,idof,innout) = - lfactor * ( - diota(1, 1,vvol) * diota(0,-1,vvol) + diota(0, 1,vvol) * diota(1,-1,vvol) ) / det
         endif
         
        endif ! end of if( Lconstraint.eq.1 ) then; 20 Jun 14;
        
#ifdef DEBUG

        if( Lcheck.eq.4 ) then ! check derivatives of field; 01 Jul 14;
         
         DifferentiateGeometry%L = .false.
         
         do isymdiff = -2, 2 ! symmetric fourth-order, finite-difference used to approximate derivatives; 26 Feb 13;
          
          if( isymdiff.eq.0 ) cycle
          
          iRbc(1:mn,0:Mvol) = oRbc(1:mn,0:Mvol)
          iZbs(1:mn,0:Mvol) = oZbs(1:mn,0:Mvol)
          iRbs(1:mn,0:Mvol) = oRbs(1:mn,0:Mvol)
          iZbc(1:mn,0:Mvol) = oZbc(1:mn,0:Mvol)
          
          if( issym.eq.0 .and. irz.eq.0 ) iRbc(ii,vvol-1+innout) = iRbc(ii,vvol-1+innout) + dRZ * isymdiff ! perturb geometry; 10 Mar 13;
          if( issym.eq.0 .and. irz.eq.1 ) iZbs(ii,vvol-1+innout) = iZbs(ii,vvol-1+innout) + dRZ * isymdiff ! perturb geometry; 10 Mar 13;
          if( issym.eq.1 .and. irz.eq.0 ) iRbs(ii,vvol-1+innout) = iRbs(ii,vvol-1+innout) + dRZ * isymdiff ! perturb geometry; 10 Mar 13;
          if( issym.eq.1 .and. irz.eq.1 ) iZbc(ii,vvol-1+innout) = iZbc(ii,vvol-1+innout) + dRZ * isymdiff ! perturb geometry; 10 Mar 13;
          
          WCALL(fc02aa,ma00aa,( Iquad(vvol), mn, vvol, Lrad(vvol) )) ! calls me00ab, which calls co01aa, which computes geometrical quantities; 10 Mar 13;
          WCALL(fc02aa,ma01ag,( vvol, mn, Lrad(vvol) )) ! construct Beltrami matrices; 13 Sep 13;
          WCALL(fc02aa,ma02aa,( vvol, Ndof )) ! this may or may not iterate on mu and dpflux to enforce transform constraints; 20 Jun 14;
          
          imupf(1:2,isymdiff) = (/ mu(vvol), dpflux(vvol) /) ! mu and dpflux are computed for the perturbed geometry by ma02aa/mp00ac if Lconstraint=1;
          
          isolution(1:Ndof,isymdiff) = solution(1:Ndof,0) ! solution is computed in mp00ac, which is called by ma02aa; 01 Jul 14;
          
         enddo ! end of do isymdiff; 10 Mar 13;
         
         isolution(1:Ndof,0) = ( - 1 * isolution(1:Ndof, 2) + 8 * isolution(1:Ndof, 1) - 8 * isolution(1:Ndof,-1) + 1 * isolution(1:Ndof,-2) ) / ( 12 * dRZ )
         imupf(1:2,0)        = ( - 1 * imupf(1:2, 2)        + 8 * imupf(1:2, 1)        - 8 * imupf(1:2,-1)        + 1 * imupf(1:2,-2)        ) / ( 12 * dRZ )
         
         ;solution(1:Ndof,-1) = abs(  solution(1:Ndof,-1) ) ! 01 Jul 14;
         isolution(1:Ndof, 0) = abs( isolution(1:Ndof, 0) )
        
         ifail = 0 ; call M01CAF(  solution(1:Ndof,-1), 1, Ndof, 'D', ifail ) ! sorting screen output; this corrupts; 01 Jul 14;
         ifail = 0 ; call M01CAF( isolution(1:Ndof, 0), 1, Ndof, 'D', ifail ) ! sorting screen output; this corrupts;
         
         cput = GETTIME

         select case( Lconstraint )
         case( 0 )
          write(ounit,3002) cput-cpus, myid, vvol, im(ii), in(ii), irz, issym, innout,  solution(1:min(Ndof,16),-1) ! 01 Jul 14;
          write(ounit,3002) cput-cpus, myid, vvol, im(ii), in(ii), irz, issym, innout, isolution(1:min(Ndof,16), 0)
         case( 1 )
          write(ounit,3003)
          write(ounit,3003) cput-cpus, myid, vvol, im(ii), in(ii), irz, issym, innout, "finite-diff", imupf(1:2,0)
          write(ounit,3003) cput-cpus, myid, vvol, im(ii), in(ii), irz, issym, innout, "analytic   ", dmupfdx(vvol,1:2,idof,innout) / lfactor
         end select ! end of select case( Lconstraint ) ; 02 Sep 14;
         
3002     format("fc02aa : ",f10.2," : ",:,"myid=",i3," ; vvol=",i3," ; (",i3," ,",i3," ) ; irz=",i2," ; issym=",i2," ; innout=",i2," ; dxi=",99es13.5)
3003     format("fc02aa : ",f10.2," : ",:,"myid=",i3," ; vvol=",i3," ; (",i3," ,",i3," ) ; irz=",i2," ; issym=",i2," ; innout=",i2," ; ",a11,&
      " : dmupf=",2f23.16" ;")
         
         DifferentiateGeometry%L = .true.
         
         iRbc(1:mn,0:Mvol) = oRbc(1:mn,0:Mvol)
         iZbs(1:mn,0:Mvol) = oZbs(1:mn,0:Mvol)
         iRbs(1:mn,0:Mvol) = oRbs(1:mn,0:Mvol)
         iZbc(1:mn,0:Mvol) = oZbc(1:mn,0:Mvol)
         
        endif ! end of if( Lcheck.eq.4 ) ; 01 Jul 14;

#endif
        
        vflag = 1 ! this flag instructs vo00aa to continue even if the volume is invalid; 04 Dec 14;
        WCALL(fc02aa,vo00aa,( vvol, vflag )) ! compute derivative of volume; wrt to harmonic described by DifferentiateGeometry structure; 20 Jun 14;
        
#ifdef DEBUG

        if( Lcheck.eq.3 ) then
         
         dvol(0) = dvolume
         
         cput = GETTIME
         write(ounit,1001) cput-cpus, myid, vvol, im(ii), in(ii), irz, issym, innout, "analytic", dvolume
         
1001     format("fc02aa : ",f10.2," : myid=",i3," ; vvol=",i3," ; (",i3," ,",i3,") ; irz=",i2," ; issym=",i2," ; innout=",i2,&
       " : ",a8," ; dvolume=",f23.15," ;",:," error=",es13.5," ;")
         
         DifferentiateGeometry%L = .false. ! instruct vo00aa to not calculate derivatives; 02 Sep 14;
         
         do isymdiff = -1, 1, 2 ! symmetric finite-difference estimate of derivative of volume wrt geometrical degree-of-freedom; 02 Sep 14;
          
          if( DifferentiateGeometry%issym.eq.0 ) then !     stellarator symmetric harmonics; 20 Jun 14;
           if( DifferentiateGeometry%irz.eq.0 ) iRbc(DifferentiateGeometry%ii,vvol-1+innout) = oRbc(DifferentiateGeometry%ii,vvol-1+innout) + isymdiff * dRZ * half
           if( DifferentiateGeometry%irz.eq.1 ) iZbs(DifferentiateGeometry%ii,vvol-1+innout) = oZbs(DifferentiateGeometry%ii,vvol-1+innout) + isymdiff * dRZ * half
          else                                        ! non-stellarator symmetric harmonics; 20 Jun 14;
           if( DifferentiateGeometry%irz.eq.0 ) iRbs(DifferentiateGeometry%ii,vvol-1+innout) = oRbs(DifferentiateGeometry%ii,vvol-1+innout) + isymdiff * dRZ * half
           if( DifferentiateGeometry%irz.eq.1 ) iZbc(DifferentiateGeometry%ii,vvol-1+innout) = oZbc(DifferentiateGeometry%ii,vvol-1+innout) + isymdiff * dRZ * half
          endif
          
          vflag = 1 ! this flag instructs vo00aa to continue even if the volume is invalid; 04 Dec 14;
          CALL(fc02aa,vo00aa,( vvol, vflag )) ! compute volume; 20 Jun 14; ! this corrupts calculation of dvolume; 20 Jun 14;
          
          dvol(isymdiff) = vvolume(vvol)
          
         enddo ! end of do isymdiff; 20 Jun 14;
         
         evolume = abs( ( dvol(+1)-dvol(-1) ) / dRZ - dvol(0) ) ! error in finite-difference calculation and analytic derivative; 20 Jun 14;
         
         cput = GETTIME
         write(ounit,1001) cput-cpus, myid, vvol, im(ii), in(ii), irz, issym, innout, "finite-d", ( dvol(+1)-dvol(-1) ) / dRZ, evolume
         
         FATALMESS(fc02aa, evolume.gt.dRZ, unacceptable error in volume derivative)

         iRbc(1:mn,0:Mvol) = oRbc(1:mn,0:Mvol)
         iZbs(1:mn,0:Mvol) = oZbs(1:mn,0:Mvol)
         iRbs(1:mn,0:Mvol) = oRbs(1:mn,0:Mvol)
         iZbc(1:mn,0:Mvol) = oZbc(1:mn,0:Mvol)
         
         DifferentiateGeometry%L = .true.
         
         dvolume = dvol(0)
        
        endif ! end of if( Lcheck.eq.3 ) ; 01 Jul 14;
        
#endif
        
        do iocons = 0, 1
         
         if( vvol.eq.   1 .and. iocons.eq.0 ) cycle ! fixed inner boundary (or coordinate axis); no constraints; 19 Sep 13;
         if( vvol.eq.Mvol .and. iocons.eq.1 ) cycle ! fixed outer boundary                     ; no constraints; 19 Sep 13;
    
         ideriv = 0 ; id = ideriv ; iflag = 0 ! why does bb00aa need to be re-called; need to determine which quantity (if any) has been corrupted; 18 Jul 14;

         WCALL(fc02aa, bb00aa,( vvol, iocons, ideriv, Ntz, dAt(1:Ntz,id), dAz(1:Ntz,id), XX(1:Ntz), YY(1:Ntz), length(1:Ntz), DDl, MMl, iflag ) )

         dBB(1:Ntz,0) = half * ( &
      dAz(1:Ntz, 0)*dAz(1:Ntz, 0)*guvij(1:Ntz,2,2,0) - two * dAz(1:Ntz, 0)*dAt(1:Ntz, 0)*guvij(1:Ntz,2,3,0) + dAt(1:Ntz, 0)*dAt(1:Ntz, 0)*guvij(1:Ntz,3,3,0) &
      ) / sg(1:Ntz,0)**2

         ideriv = -1 ; id = ideriv ; iflag = 1 ! compute derivatives of magnetic field; 18 Jul 14;
         
         WCALL(fc02aa, bb00aa,( vvol, iocons, ideriv, Ntz, dAt(1:Ntz,id), dAz(1:Ntz,id), XX(1:Ntz), YY(1:Ntz), length(1:Ntz), DDl, MMl, iflag ) )
         
         lss = two * iocons - one ; Lcurvature = 4
         WCALL(fc02aa,co01aa,( vvol, lss, Lcurvature, Ntz, mn )) ! get coordinate metrics and their derivatives wrt Rj, Zj on interface;
         
         dBB(1:Ntz,id) = half * ( &
      dAz(1:Ntz,id)*dAz(1:Ntz, 0)*guvij(1:Ntz,2,2,0) - two * dAz(1:Ntz,id)*dAt(1:Ntz, 0)*guvij(1:Ntz,2,3,0) + dAt(1:Ntz,id)*dAt(1:Ntz, 0)*guvij(1:Ntz,3,3,0) &
    + dAz(1:Ntz, 0)*dAz(1:Ntz,id)*guvij(1:Ntz,2,2,0) - two * dAz(1:Ntz, 0)*dAt(1:Ntz,id)*guvij(1:Ntz,2,3,0) + dAt(1:Ntz, 0)*dAt(1:Ntz,id)*guvij(1:Ntz,3,3,0) &
    + dAz(1:Ntz, 0)*dAz(1:Ntz, 0)*guvij(1:Ntz,2,2,1) - two * dAz(1:Ntz, 0)*dAt(1:Ntz, 0)*guvij(1:Ntz,2,3,1) + dAt(1:Ntz, 0)*dAt(1:Ntz, 0)*guvij(1:Ntz,3,3,1) &
    ) / sg(1:Ntz,0)**2                                                                                                                                       &
    - dBB(1:Ntz,0) * two * sg(1:Ntz,1) / sg(1:Ntz,0)

         FATALMESS(fc02aa, vvolume(vvol).lt.small, shall divide by vvolume(vvol)**(gamma+one))

         ijreal(1:Ntz) = - adiabatic(vvol) * pscale * gamma * dvolume / vvolume(vvol)**(gamma+one) + dBB(1:Ntz,-1) ! derivatives of force wrt geometry;

         dLL(1:Ntz) = zero ! either no spectral constraint, or not the appropriate interface; 19 Sep 13;
         dPP(1:Ntz) = zero ! either no spectral constraint, or not the appropriate interface; 19 Sep 13;
         
         if( Igeometry.ge.3 ) then ! spectral constraints are only required in toroidal or extended-cylindrical geometry;
          
          if( innout.eq.1 .and. iocons.eq.1 ) then ! include derivatives of spectral constraints; 19 Sep 13;
#ifdef DEBUG
             FATALMESS(fc02aa, abs(DDl).lt.small, divide by zero on spectral constraint)
#endif
           if( issym.eq.0 ) then ! take derivatives wrt Rbc and Zbs, 30 Oct 13;
            if( irz.eq.0 ) then ! take derivative wrt Rbc; 11 Aug 14;
             dII(1:Ntz) = - im(ii) * sini(1:Ntz,ii) * ( XX(1:Ntz) - MMl * iRij(1:Ntz,vvol) ) &
                          - two * ( mpnq(ii) - MMl ) * iRbc(ii,vvol) * ( Rij(1:Ntz,2,0) * iRij(1:Ntz,vvol) + Zij(1:Ntz,2,0) * iZij(1:Ntz,vvol) ) / DDl &
                          + Rij(1:Ntz,2,0) * ( mpnq(ii) - MMl ) * cosi(1:Ntz,ii)
            else ! take derivative wrt Zbs; 11 Aug 14;
             dII(1:Ntz) = + im(ii) * cosi(1:Ntz,ii) * ( YY(1:Ntz) - MMl * iZij(1:Ntz,vvol) ) &
                          - two * ( mpnq(ii) - MMl ) * iZbs(ii,vvol) * ( Rij(1:Ntz,2,0) * iRij(1:Ntz,vvol) + Zij(1:Ntz,2,0) * iZij(1:Ntz,vvol) ) / DDl &
                          + Zij(1:Ntz,2,0) * ( mpnq(ii) - MMl ) * sini(1:Ntz,ii)
            endif ! end of if( irz.eq.0 ) ; 11 Aug 14;
           else                  ! take derivatives wrt Rbs and Zbc; 30 Oct 13;
            if( irz.eq.0 ) then 
             dII(1:Ntz) = + im(ii) * cosi(1:Ntz,ii) * ( XX(1:Ntz) - MMl * iRij(1:Ntz,vvol) ) &
                          - two * ( mpnq(ii) - MMl ) * iRbs(ii,vvol) * ( Rij(1:Ntz,2,0) * iRij(1:Ntz,vvol) + Zij(1:Ntz,2,0) * iZij(1:Ntz,vvol) ) / DDl &
                          + Rij(1:Ntz,2,0) * ( mpnq(ii) - MMl ) * sini(1:Ntz,ii)
            else                
             dII(1:Ntz) = - im(ii) * sini(1:Ntz,ii) * ( YY(1:Ntz) - MMl * iZij(1:Ntz,vvol) ) &
                          - two * ( mpnq(ii) - MMl ) * iZbc(ii,vvol) * ( Rij(1:Ntz,2,0) * iRij(1:Ntz,vvol) + Zij(1:Ntz,2,0) * iZij(1:Ntz,vvol) ) / DDl &
                          + Zij(1:Ntz,2,0) * ( mpnq(ii) - MMl ) * cosi(1:Ntz,ii)
            endif
           endif

          else

           dII(1:Ntz) = zero ! either no spectral constraint, or not the appropriate interface; 19 Sep 13;

          endif ! end of if( innout.eq.1 .and. iocons.eq.1 ) ; 08 Nov 13;
         
          if( iocons.eq.0 ) then ! take derivatives of constraints at inner boundary; 08 Nov 13;

           constraint(1:Ntz) = + ( dRij(1:Ntz,vvol) * tRij(1:Ntz,vvol-1) + dZij(1:Ntz,vvol) * tZij(1:Ntz,vvol-1) ) / length(1:Ntz)
           
           if( innout.eq.0 ) then ! derivative wrt inner boundary coefficient; 08 Nov 13;
            if( issym.eq.0 ) then ! take derivatives wrt Rbc and Zbs, 30 Oct 13;
             if( irz.eq.0 ) then ; dLL(1:Ntz) = + ( - cosi(1:Ntz,ii) * tRij(1:Ntz,vvol-1) - dRij(1:Ntz,vvol) * im(ii) * sini(1:Ntz,ii) ) / length(1:Ntz) & 
                                                + constraint(1:Ntz) * dRij(1:Ntz,vvol) * cosi(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
             else                ; dLL(1:Ntz) = + ( - sini(1:Ntz,ii) * tZij(1:Ntz,vvol-1) + dZij(1:Ntz,vvol) * im(ii) * cosi(1:Ntz,ii) ) / length(1:Ntz) & 
                                                + constraint(1:Ntz) * dZij(1:Ntz,vvol) * sini(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
             endif
            else ! if issym.eq.1 ; take derivatives wrt Rbs and Zbc; 30 Oct 13;
             if( irz.eq.0 ) then ; dLL(1:Ntz) = + ( - sini(1:Ntz,ii) * tRij(1:Ntz,vvol-1) + dRij(1:Ntz,vvol) * im(ii) * cosi(1:Ntz,ii) ) / length(1:Ntz) & 
                                                + constraint(1:Ntz) * dRij(1:Ntz,vvol) * sini(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
             else                ; dLL(1:Ntz) = + ( - cosi(1:Ntz,ii) * tZij(1:Ntz,vvol-1) - dZij(1:Ntz,vvol) * im(ii) * sini(1:Ntz,ii) ) / length(1:Ntz) & 
                                                + constraint(1:Ntz) * dZij(1:Ntz,vvol) * cosi(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
             endif
            endif
           else ! if innout.eq.1 ; derivative wrt outer boundary coefficient; 08 Nov 13;
            if( issym.eq.0 ) then ! take derivatives wrt Rbc and Zbs, 30 Oct 13;
             if( irz.eq.0 ) then ; dLL(1:Ntz) = + ( + cosi(1:Ntz,ii) * tRij(1:Ntz,vvol-1)                                              ) / length(1:Ntz) & 
                                                - constraint(1:Ntz) * dRij(1:Ntz,vvol) * cosi(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
             else                ; dLL(1:Ntz) = + ( + sini(1:Ntz,ii) * tZij(1:Ntz,vvol-1)                                              ) / length(1:Ntz) & 
                                                - constraint(1:Ntz) * dZij(1:Ntz,vvol) * sini(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
             endif
            else ! if issym.eq.1 ; take derivatives wrt Rbs and Zbc; 30 Oct 13;
             if( irz.eq.0 ) then ; dLL(1:Ntz) = + ( + sini(1:Ntz,ii) * tRij(1:Ntz,vvol-1)                                              ) / length(1:Ntz) & 
                                                - constraint(1:Ntz) * dRij(1:Ntz,vvol) * sini(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
             else                ; dLL(1:Ntz) = + ( + cosi(1:Ntz,ii) * tZij(1:Ntz,vvol-1)                                              ) / length(1:Ntz) & 
                                                - constraint(1:Ntz) * dZij(1:Ntz,vvol) * cosi(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
             endif
            endif
           endif

          else ! if iocons.eq.1 ; take derivatives of constraints at outer boundary; 08 Nov 13;

           constraint(1:Ntz) = + ( dRij(1:Ntz,vvol) * tRij(1:Ntz,vvol  ) + dZij(1:Ntz,vvol) * tZij(1:Ntz,vvol  ) ) / length(1:Ntz)
           
           if( innout.eq.0 ) then ! derivative wrt inner boundary coefficient; 08 Nov 13;
            if( issym.eq.0 ) then ! take derivatives wrt Rbc and Zbs, 30 Oct 13;
             if( irz.eq.0 ) then ; dLL(1:Ntz) = + ( - cosi(1:Ntz,ii) * tRij(1:Ntz,vvol  )                                              ) / length(1:Ntz) & 
                                                + constraint(1:Ntz) * dRij(1:Ntz,vvol) * cosi(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
             else                ; dLL(1:Ntz) = + ( - sini(1:Ntz,ii) * tZij(1:Ntz,vvol  )                                              ) / length(1:Ntz) & 
                                                + constraint(1:Ntz) * dZij(1:Ntz,vvol) * sini(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
             endif
            else                  ! take derivatives wrt Rbs and Zbc; 30 Oct 13;
             if( irz.eq.0 ) then ; dLL(1:Ntz) = + ( - sini(1:Ntz,ii) * tRij(1:Ntz,vvol  )                                              ) / length(1:Ntz) & 
                                                + constraint(1:Ntz) * dRij(1:Ntz,vvol) * sini(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
             else                ; dLL(1:Ntz) = + ( - cosi(1:Ntz,ii) * tZij(1:Ntz,vvol  )                                              ) / length(1:Ntz) & 
                                                + constraint(1:Ntz) * dZij(1:Ntz,vvol) * cosi(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
             endif
            endif
           else ! if innout.eq.1 ; derivative wrt outer boundary coefficient; 08 Nov 13;
            if( issym.eq.0 ) then ! take derivatives wrt Rbc and Zbs, 30 Oct 13;
             if( irz.eq.0 ) then ; dLL(1:Ntz) = + ( + cosi(1:Ntz,ii) * tRij(1:Ntz,vvol  ) - dRij(1:Ntz,vvol) * im(ii) * sini(1:Ntz,ii) ) / length(1:Ntz) & 
                                                - constraint(1:Ntz) * dRij(1:Ntz,vvol) * cosi(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
             else                ; dLL(1:Ntz) = + ( + sini(1:Ntz,ii) * tZij(1:Ntz,vvol  ) + dZij(1:Ntz,vvol) * im(ii) * cosi(1:Ntz,ii) ) / length(1:Ntz) & 
                                                - constraint(1:Ntz) * dZij(1:Ntz,vvol) * sini(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
             endif
            else                  ! take derivatives wrt Rbs and Zbc; 30 Oct 13;
             if( irz.eq.0 ) then ; dLL(1:Ntz) = + ( + sini(1:Ntz,ii) * tRij(1:Ntz,vvol  ) + dRij(1:Ntz,vvol) * im(ii) * cosi(1:Ntz,ii) ) / length(1:Ntz) & 
                                                - constraint(1:Ntz) * dRij(1:Ntz,vvol) * sini(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
             else                ; dLL(1:Ntz) = + ( + cosi(1:Ntz,ii) * tZij(1:Ntz,vvol  ) - dZij(1:Ntz,vvol) * im(ii) * sini(1:Ntz,ii) ) / length(1:Ntz) & 
                                                - constraint(1:Ntz) * dZij(1:Ntz,vvol) * cosi(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
             endif
            endif
           endif
          endif ! end of if( iocons.eq.0 ) ; 01 Jul 14;
          
         endif ! end of if( Igeometry.ge.3 ) ; 08 Nov 13;

         call tfft( Nt, Nz, ijreal(1:Ntz), dII(1:Ntz), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz), & ! recall that ijreal contains pressure term;
                    mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), ifail )
         

         call tfft( Nt, Nz, dPP(1:Ntz)   , dLL(1:Ntz), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz), & ! recall that ijreal is probably just a dummy;
                    mn, im(1:mn), in(1:mn), evmn(1:mn), odmn(1:mn), comn(1:mn), simn(1:mn), ifail )          ! evmn and odmn are available as workspace;


         ; idoc = 0
         
         ;  dFFdRZ(idoc+1:idoc+mn    ,vvol,iocons,idof,innout) = + efmn(1:mn    ) * psifactor(ii,vvol-1+innout) * expmmnn(1:mn) * opsilon ! pressure;
         ; idoc = idoc + mn   ! even; 19 Sep 13;
         ;if( Igeometry.ge.3 ) then ! add spectral constraints; 19 Sep 13;
         ;  dFFdRZ(idoc+1:idoc+mn-1  ,vvol,iocons,idof,innout) = - sfmn(2:mn    ) * psifactor(ii,vvol-1+innout) * epsilon       & ! spectral condensation;
                                                                 - simn(2:mn    ) * psifactor(ii,vvol-1+innout) * sweight(vvol)   ! poloidal length constraint;
         ! if( Ntor.gt.0 ) then
         !  dFFdRZ(idoc+1:idoc+Ntor  ,vvol,iocons,idof,innout) = + odmn(2:Ntor+1) * psifactor(ii,vvol-1+innout) * apsilon
         ! endif
         ; ;idoc = idoc + mn-1 ! oddd; 19 Sep 13;
         ;endif ! end of if( Igeometry.ge.3) ; 20 Jun 14;
         if( NOTstellsym ) then
          ; dFFdRZ(idoc+1:idoc+mn-1  ,vvol,iocons,idof,innout) = + ofmn(2:mn    ) * psifactor(ii,vvol-1+innout) * expmmnn(1:mn) * opsilon ! pressure;
          ;idoc = idoc + mn-1 ! oddd; 19 Sep 13;
          if( Igeometry.ge.3 ) then ! add spectral constraints; 19 Sep 13;
           ;dFFdRZ(idoc+1:idoc+mn    ,vvol,iocons,idof,innout) = - cfmn(1:mn    ) * psifactor(ii,vvol-1+innout) * epsilon       & ! spectral condensation;
                                                                 - comn(1:mn    ) * psifactor(ii,vvol-1+innout) * sweight(vvol)   ! poloidal length constraint;
          !if( Ntor.ge.0 ) then
          ! dFFdRZ(idoc+1:idoc+Ntor+1,vvol,iocons,idof,innout) = + evmn(1:Ntor+1) * psifactor(ii,vvol-1+innout) * apsilon ! poloidal origin      ; 18 Jul 14;
          !endif
           idoc = idoc + mn   ! even; 19 Sep 13;
          endif ! end of if( Igeometry.ge.3) ; 20 Jun 14;
         endif ! end of if( NOTstellsym) ; 20 Jun 14;
         
#ifdef DEBUG
         FATALMESS( fc02aa, idoc.ne.lgeometricaldof, counting error)
#endif
         
        enddo ! end of do iocons, 19 Sep 13;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
        
!#ifdef SECONDDERIVATIVES
!
! inside do vvol; inside do ii; inside do irz; inside do innout; inside do iocons;
!        
!        select case( irz )
!        case( 0 ) ! compute derivative wrt R;
!         if( innout.eq.iocons ) then 
!          ijreal(1:Ntz) = ( - adiabatic(vvol) * pscale * gamma * dvolume / vvolume(vvol)**(one+gamma) + half * ddBB(1:Ntz) ) * Rij(0,:)    * Zij(2,:) &
!                        + (   adiabatic(vvol) * pscale                   / vvolume(vvol)**(    gamma) + half *  dBB(1:Ntz) ) * cosi(ii,:)  * Zij(2,:) 
!
!          ijimag(1:Ntz) = ( - adiabatic(vvol) * pscale * gamma * dvolume / vvolume(vvol)**(one+gamma) + half * ddBB(1:Ntz) ) * Rij(0,:)    * Rij(2,:) &
!                        + (   adiabatic(vvol) * pscale                   / vvolume(vvol)**(    gamma) + half *  dBB(1:Ntz) ) * cosi(ii,:)  * Rij(2,:) &
!                        + (   adiabatic(vvol) * pscale                   / vvolume(vvol)**(    gamma) + half *  dBB(1:Ntz) ) * Rij(0,:)    * sini(ii,:) * ( - im(ii) )
!         else
!          ijreal(1:Ntz) = ( - adiabatic(vvol) * pscale * gamma * dvolume / vvolume(vvol)**(one+gamma) + half * ddBB(1:Ntz) ) * Rij(0,:)    * Zij(2,:)
!
!          ijimag(1:Ntz) = ( - adiabatic(vvol) * pscale * gamma * dvolume / vvolume(vvol)**(one+gamma) + half * ddBB(1:Ntz) ) * Rij(0,:)    * Rij(2,:)
!         endif
!        case( 1 ) ! compute derivative wrt Z;
!         if( innout.eq.iocons ) then
!          ijreal(1:Ntz) = ( - adiabatic(vvol) * pscale * gamma * dvolume / vvolume(vvol)**(one+gamma) + half * ddBB(1:Ntz) ) * Rij(0,:)    * Zij(2,:) &
!                        + (   adiabatic(vvol) * pscale                   / vvolume(vvol)**(    gamma) + half *  dBB(1:Ntz) ) * Rij(0,:)    * cosi(ii,:) * (   im(ii) )
!
!          ijimag(1:Ntz) = ( - adiabatic(vvol) * pscale * gamma * dvolume / vvolume(vvol)**(one+gamma) + half * ddBB(1:Ntz) ) * Rij(0,:)    * Rij(2,:)
!         else
!          ijreal(1:Ntz) = ( - adiabatic(vvol) * pscale * gamma * dvolume / vvolume(vvol)**(one+gamma) + half * ddBB(1:Ntz) ) * Rij(0,:)    * Zij(2,:)
!
!          ijimag(1:Ntz) = ( - adiabatic(vvol) * pscale * gamma * dvolume / vvolume(vvol)**(one+gamma) + half * ddBB(1:Ntz) ) * Rij(0,:)    * Rij(2,:)
!         endif
!        case default
!         FATALMESS(fc02aa, .true., must be R or Z)
!        end select
!        
!        call tfft( Nt,Nz,ijreal(1:Ntz),ijimag(1:Ntz),isr,trigm,trign,trigwk,mn,im(1:mn),in(1:mn),efmn(1:mn),ofmn(1:mn),cfmn(1:mn),sfmn(1:mn),ifail )       
!        
!        jdof = 0
!        do jmn = 1, mn
!         lfactor = pi * pi2nfp * tflux(vvol-1+innout)**(im(jmn)*half) * (-1)**(innout-1) * ( -Itoroidal ) ! ! DEFINITION OF TFLUX MAY HAVE CHANGED; included pi2; 01 Jul 14;
!         do jrz = 0, 1
!          if( Igeometry.le.3 .and. jrz.eq.1 ) cycle ! no dependence on Z;
!          if( jmn.eq.1 .and. jrz.eq.1 ) cycle
!          jdof = jdof + 1
!          if( jrz.eq.0 .and. jmn.eq.1 ) DanalyDrzDrz(vvol,iocons,jdof,innout,idof) = - efmn(jmn) * two * lfactor ! R;
!          if( jrz.eq.0 .and. jmn.gt.1 ) DanalyDrzDrz(vvol,iocons,jdof,innout,idof) = - efmn(jmn)       * lfactor ! R;
!         !if( jrz.eq.1 .and. jmn.eq.1 ) DanalyDrzDrz(vvol,iocons,jdof,innout,idof) =   sfmn(jmn) * two * lfactor ! Z; this mode is always irrelevant;
!          if( jrz.eq.1 .and. jmn.gt.1 ) DanalyDrzDrz(vvol,iocons,jdof,innout,idof) =   sfmn(jmn)       * lfactor ! Z;
!         enddo ! end of do jrz;
!        enddo ! end of do jmn;
!
!#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
        
       enddo ! matches do innout;  7 Mar 13;
       
      enddo ! matches do issym;  7 Mar 13;

     enddo ! matches do irz;  7 Mar 13;

    enddo ! matches do ii;  7 Mar 13;
    
    DifferentiateGeometry%L = .false. ! probably not needed, but included anyway; 13 Mar 13;
    
   endif ! end of if( LComputeDerivatives) ;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

2000 continue
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   WCALL(fc02aa,ma00ab,( 'D', vvol )) ! deallocate volume integrated metric elements; 16 Jan 13;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   DEALLOCATE(dMA) ! 17 Jan 13;
   DEALLOCATE(dMB) ! 17 Jan 13;
   
   DEALLOCATE(solution) ! 21 Apr 13;

   
   if( Lplasmaregion ) then
    
    DEALLOCATE(dMC) ! 17 Jan 13;
    DEALLOCATE(dMD) ! 17 Jan 13;
    DEALLOCATE(dME) ! 17 Jan 13;
    DEALLOCATE(dMF) ! 17 Jan 13;
    
    DEALLOCATE(MBpsi) ! 26 Feb 13;
    DEALLOCATE(MEpsi) ! 26 Feb 13;
    
   endif
   
   if( LcomputeDerivatives) then
    
    DEALLOCATE(oBI) ! 17 Jan 13;
    DEALLOCATE(rhs)
    
#ifdef DEBUG
    if( Lcheck.eq.4 ) then
     DEALLOCATE(isolution)
    endif
#endif
    
   endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  enddo ! end of do vvol = 1, Mvol ; 16 Jan 13;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  DEALLOCATE(dAt)
  DEALLOCATE(dAz)
  
  DEALLOCATE( XX) ! spectral constraints; 20 Feb 13; ! not used; 19 Sep 13;
  DEALLOCATE( YY)
  
  DEALLOCATE(length)

  DEALLOCATE(dRR)
  DEALLOCATE(dZZ)

  DEALLOCATE(dII)
  DEALLOCATE(dLL)
  DEALLOCATE(dPP)

  DEALLOCATE(constraint)

  DEALLOCATE(dBB)

#ifdef DEBUG
  if( Lcheck.eq.3 .or. Lcheck.eq.4 ) then
   DEALLOCATE(oRbc)
   DEALLOCATE(oZbs)
   DEALLOCATE(oRbs)
   DEALLOCATE(oZbc)
  endif
#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  if( Lcheck.eq.2 ) then
   write(ounit,'("fc02aa : ", 10x ," : myid=",i3," ; finished computing derivatives of rotational-transform wrt mu and dpflux ")') myid
   stop "fc02aa :            : myid=    ; finished computing derivatives of rotational-transform wrt mu and dpflux " ! this will allow other cpus to finish;
  endif
  FATALMESS(fc02aa, Lcheck.eq.2, finished computing derivatives of rotational-transform wrt mu and dpflux)
#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  do vvol = 1, Mvol ! broadcast;
   
   if( vvol.le.Nvol ) Lplasmaregion = .true.
   if( vvol.gt.Nvol ) Lplasmaregion = .false.

   Lvacuumregion = .not. Lplasmaregion

   WCALL(fc02aa,bc00aa,( vvol ))
   
  enddo
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  lBBintegral(1:Nvol) = lBBintegral(1:Nvol) * half
  
  Energy = sum( lBBintegral(1:Nvol) ) ! should also compute beta;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! CONSTRUCT FORCE; 15 Jan 15;

 !packorunpack = 'U' ! unpack geometrical degrees-of-freedom; 13 Sep 13;
 !WCALL(fc02aa,gf00aa,( Ngeometricaldof, position(0:Ngeometricaldof), Mvol, mn, iRbc(1:mn,0:Mvol), iZbs(1:mn,0:Mvol), iRbs(1:mn,0:Mvol), iZbc(1:mn,0:Mvol), packorunpack ))

 !write(ounit,'("fc02aa : " 10x " : position="999f15.10)') position(1:lgeometricaldof)
 !write(ounit,'("fc02aa : " 10x " : iRbc    ="999f15.10)') iRbc(1:mn,1)
 !pause

  ;   force(0:Ngeometricaldof) = zero
  
  do vvol = 1, Mvol-1 ! loop over interior surfaces
   
   if( Igeometry.eq.1 .or. vvol.gt.1 ) then ; Lcoordinatesingularity = .false. ! 14 Jan 15;
   else                                     ; Lcoordinatesingularity = .true.  ! 14 Jan 15;
   endif
   
   tdoc = (vvol-1) * lgeometricaldof 
   
   if( ImagneticOK(vvol) .and. ImagneticOK(vvol+1) ) then ! the magnetic fields in the volumes adjacent to this interface are valid; 19 Sep 13;
    
    ;  idoc = 0
    
!latex \item The extrapolation constraint is $R_{j,1} = R_{j,2} \, \psi_1^{m/2} / \psi_2^{m/2}$.
!latex       Combining this with the regularization factor for the geometry, i.e. $R_{j,i}=\psi_{i}^{m/2} \xi_{j,i}$, we obtain
!latex       \be \xi_{j,1} = R_{j,2} / \psi_2^{m/2}.
!latex       \ee

    if( Lextrap.eq.1 .and. vvol.eq.1 ) then ! 20 Jan 15;
     force(tdoc+idoc+1:tdoc+idoc+mn) = position(1:mn) - ( iRbc(1:mn,2) / psifactor(1:mn,2) )
    else
     force(tdoc+idoc+1:tdoc+idoc+mn    ) = ( Bemn(1:mn    ,vvol+1,0) - Bemn(1:mn    ,vvol+0,1) ) * expmmnn(1:mn) * opsilon ! pressure imbalance; 29 Apr 14;
    endif
    
    ;  BBe(vvol) = max( sum( abs( force(tdoc+idoc+1:tdoc+idoc+mn  ) ) ) / (mn  ), logtolerance ) ! screen diagnostics; 29 Apr 14;
    
    ;  idoc = idoc + mn
    
    if( Igeometry.ge.3 ) then ! add spectral constraints; 29 Apr 14;
     
     ; force(tdoc+idoc+1:tdoc+idoc+mn-1  ) = (                           Iomn(2:mn    ,vvol+0  ) ) * epsilon         & ! spectral constraints; 29 Apr 14;
                                           + (                         + Somn(2:mn    ,vvol+0,1) ) * sweight(vvol+0) & ! poloidal length constraint; sweight contains upsilon;
                                           - ( Somn(2:mn    ,vvol+1,0)                           ) * sweight(vvol+1)

!     if( Ntor.gt.0 ) then ! poloidal angle origin is not otherwise constrained ; 18 Jul 14;
!      force(tdoc+idoc+1:tdoc+idoc+Ntor  ) = ( Pomn(2:Ntor+1,vvol+1,0) - Pomn(2:Ntor+1,vvol+0,1) ) * apsilon ! choice of spectral constraint can be enforced; 08 Nov 13;
!     endif
     
     ;IIo(vvol) = max( sum( abs( force(tdoc+idoc+1:tdoc+idoc+mn-1) ) ) / (mn-1), logtolerance ) ! screen diagnostics; 29 Apr 14;

     ; idoc = idoc + mn-1
     
    endif ! end of if( Igeometry.ge.3 ) ; 15 Jan 15;
    
    if( NOTstellsym ) then
     
     if( Lcoordinatesingularity ) then
      FATALMESS(fc02aa, .true., coordinate singularity under re-construction)
     else
      ; force(tdoc+idoc+1:tdoc+idoc+mn-1  ) = ( Bomn(2:mn    ,vvol+1,0) - Bomn(2:mn    ,vvol+0,1) ) * expmmnn(1:mn) * opsilon ! pressure imbalance; 29 Apr 14;
     endif
     
     ; BBo(vvol) = max( sum( abs( force(tdoc+idoc+1:tdoc+idoc+mn-1) ) ) / (mn-1), logtolerance ) ! screen diagnostics; 29 Apr 14;
     
     ; idoc = idoc + mn-1
     
     if( Igeometry.ge.3 ) then ! add spectral constraints; 29 Apr 14;
     
      ;force(tdoc+idoc+1:tdoc+idoc+mn    ) = (                           Iemn(1:mn    ,vvol+0  ) ) * epsilon         & ! spectral constraints; 29 Apr 14;
                                           + (                         + Semn(1:mn    ,vvol+0,1) ) * sweight(vvol+0) & ! poloidal length constraint; sweight contains upsilon;
                                           - ( Semn(1:mn    ,vvol+1,0)                           ) * sweight(vvol+1)
      
!     ;force(tdoc+idoc+1:tdoc+idoc+mn    ) = (                           Iemn(1:mn    ,vvol+0  ) ) * epsilon &
!                                          - ( Semn(1:mn    ,vvol+1,0) - Semn(1:mn    ,vvol+0,1) ) * sweight(vvol-1+innout)

!     if( Ntor.ge.0 ) then
!      force(tdoc+idoc+1:tdoc+idoc+Ntor+1) = ( Pemn(1:Ntor+1,vvol+1,0) - Pemn(1:Ntor+1,vvol+0,1) ) * apsilon ! choice of spectral constraint can be enforced; 08 Nov 13;
!     endif

      ;IIe(vvol) = max( sum( abs( force(tdoc+idoc+1:tdoc+idoc+mn  ) ) ) / (mn  ), logtolerance ) ! screen diagnostics; 29 Apr 14;

      ;idoc = idoc + mn

     endif ! end of if( Igeometry.ge.3 ) ; 15 Jan 15;

    endif ! end of if( NOTstellsym ) ; 15 Jan 15;
    
#ifdef DEBUG
    FATALMESS( fc02aa, idoc.ne.lgeometricaldof, counting error) ! this has caught bugs; 08 Nov 13;
#endif
    
   else ! magnetic field calculation in an adjacent volume failed; 20 Feb 13;
    
    ;  BBe(vvol) = 9.9E+09
    ;  IIo(vvol) = 9.9E+09
    if ( NOTstellsym ) then
     ; BBo(vvol) = 9.9E+09
     ; IIe(vvol) = 9.9E+09
    endif

    ; force(tdoc+1:tdoc+lgeometricaldof) = 9.9E+09
    
   endif ! end of if( ImagneticOK(vvol) .and. ImagneticOK(vvol+1) ) ; 20 Feb 13;
   
  enddo ! end of do vvol; 20 Feb 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Ngeometricaldof.ne.0 ) then ; ForceErr = sqrt( sum( force(1:Ngeometricaldof)*force(1:Ngeometricaldof) ) / Ngeometricaldof ) ! this includes spectral constraints; 26 Feb 13;
  else                            ; ForceErr = zero
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG

  if( Wfc02aa ) then
   
   if( myid.eq.0 ) then
    cput = GETTIME
    ; write(ounit,4000) cput-cpus,                                   ForceErr,  cput-lastcpu, "|BB|e",        BBe(1:min(Mvol-1,28))
    if( Igeometry.ge.3 ) then
     ;write(ounit,4001)                                                                       "|II|o",        IIo(1:min(Mvol-1,28))
    endif
    if( NOTstellsym ) then
     ;write(ounit,4001)                                                                       "|BB|o",        BBo(1:min(Mvol-1,28))
     if( Igeometry.ge.3 ) then
      write(ounit,4001)                                                                       "|II|e",        IIe(1:min(Mvol-1,28))
     endif
    endif
   endif

  endif ! end of if( Wfc02aa ) ; 11 Aug 14;
  
#endif

4000 format("fc02aa : ",f10.2," : ",6x,3x," ",:,"|f|=",es12.5," ; ",:,"time=",f10.2,"s ;",:,"    ",a5,"=",28es13.5," ...")
4001 format("fc02aa : ", 10x ," : ",6x,3x," ",:,"    ",  12x ,"   ",:,"     ", 10x ,"  ;",:,"    ",a5,"=",28es13.5," ...")

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( LcomputeDerivatives ) then ! construct Hessian; 19 Sep 13;
   
#ifdef DEBUG
   FATALMESS(fc02aa, .not.Lhessianallocated, need to allocate hessian)
#endif
   
   hessian(1:Ngeometricaldof,1:Ngeometricaldof) = zero 
   
   do vvol = 1, Mvol-1 ! loop over interior surfaces;
    
   !if( Igeometry.eq.1 .or. vvol.gt.1 ) then ; Lcoordinatesingularity = .false. ! 14 Jan 15;
   !else                                     ; Lcoordinatesingularity = .true.  ! 14 Jan 15;
   !endif
    
    if( ImagneticOK(vvol) .and. ImagneticOK(vvol+1) ) then ! the magnetic fields in the volumes adjacent to this interface are valid; 19 Sep 13;
     
     idof = 0 ! labels degree-of-freedom = Fourier harmonic of surface geometry; 18 Dec 14;
     
     do ii = 1, mn ! loop over degrees-of-freedom; 19 Sep 13;
      
      do irz = 0, 1 ! Fourier harmonic of R, Fourier harmonic of Z; 18 Dec 14;
       
       if( irz.eq.1 .and. Igeometry.lt.3 ) cycle ! no dependence on Z; 26 Feb 13;
       
       do issym = 0, 1 ! stellarator symmetry; 18 Dec 14;
        
        if( issym.eq.1 .and. YESstellsym ) cycle ! no dependence on the non-stellarator symmetric harmonics; 26 Feb 13;
        
        if( ii.eq.1 .and. irz.eq.1 .and. issym.eq.0 ) cycle ! no dependence on Zbs_{m=0,n=0}; 26 Feb 13;
        if( ii.eq.1 .and. irz.eq.0 .and. issym.eq.1 ) cycle ! no dependence on Rbs_{m=0,n=0}; 26 Feb 13;
        
        idof = idof + 1 ! labels degree-of-freedom; 18 Dec 14;
        
        if( vvol.gt.1 ) then
         
         tdof = (vvol-2) * lgeometricaldof + idof ! labels degree-of-freedom in internal interface geometry   ; 18 Dec 14;
         tdoc = (vvol-1) * lgeometricaldof        ! labels force-balance constraint across internal interfaces; 19 Sep 13;
         idoc = 0                                 ! local  force-balance constraint across internal interface ; 18 Dec 14;
         hessian(tdoc+idoc+1:tdoc+idoc+lgeometricaldof,tdof) =                                                     - dFFdRZ(idoc+1:idoc+lgeometricaldof,vvol+0,1,idof,0)
         if( Lconstraint.eq.1 ) then
         hessian(tdoc+idoc+1:tdoc+idoc+lgeometricaldof,tdof) = hessian(tdoc+idoc+1:tdoc+idoc+lgeometricaldof,tdof)                     &
                                                             - dBBdmp(idoc+1:idoc+lgeometricaldof,vvol+0,1,1) * dmupfdx(vvol,1,idof,0) &
                                                             - dBBdmp(idoc+1:idoc+lgeometricaldof,vvol+0,1,2) * dmupfdx(vvol,2,idof,0)

         endif ! end of if( Lconstraint.eq.1 ) ; 20 Jun 14;
      
        endif ! end of if( vvol.gt.1 ) ; 19 Sep 13;
        

        ;tdof = (vvol-1) * lgeometricaldof + idof
        ;tdoc = (vvol-1) * lgeometricaldof ! shorthand; 19 Sep 13;         
        ;idoc = 0
        if( Lextrap.eq.1 .and. vvol.eq.1 ) then ! 20 Jan 15;
        ;hessian(tdoc+idof                            ,tdof) = one ! diagonal elements; 20 Jan 15;
        else
        ;hessian(tdoc+idoc+1:tdoc+idoc+lgeometricaldof,tdof) = dFFdRZ(idoc+1:idoc+lgeometricaldof,vvol+1,0,idof,0) - dFFdRZ(idoc+1:idoc+lgeometricaldof,vvol+0,1,idof,1)
         if( Lconstraint.eq.1 ) then
         hessian(tdoc+idoc+1:tdoc+idoc+lgeometricaldof,tdof) = hessian(tdoc+idoc+1:tdoc+idoc+lgeometricaldof,tdof)                       &
                                                             + dBBdmp(idoc+1:idoc+lgeometricaldof,vvol+1,0,1) * dmupfdx(vvol+1,1,idof,0) &
                                                             + dBBdmp(idoc+1:idoc+lgeometricaldof,vvol+1,0,2) * dmupfdx(vvol+1,2,idof,0) &
                                                             - dBBdmp(idoc+1:idoc+lgeometricaldof,vvol+0,1,1) * dmupfdx(vvol+0,1,idof,1) &
                                                             - dBBdmp(idoc+1:idoc+lgeometricaldof,vvol+0,1,2) * dmupfdx(vvol+0,2,idof,1)
         endif ! end of if( Lconstraint.eq.1 ); 20 Jun 14;
         endif
        
         if( vvol.lt.Mvol-1 ) then

         tdof = (vvol+0) * lgeometricaldof + idof
         tdoc = (vvol-1) * lgeometricaldof ! shorthand; 19 Sep 13;         
         idoc = 0
         if( Lextrap.eq.1 .and. vvol.eq.1 ) then ! 20 Jan 15;
         if    ( im(idof).le.0                     ) then ; hessian(tdoc+idof,tdof) = - one
         else                                             ; hessian(tdoc+idof,tdof) = - one
         endif
         else
         hessian(tdoc+idoc+1:tdoc+idoc+lgeometricaldof,tdof) = dFFdRZ(idoc+1:idoc+lgeometricaldof,vvol+1,0,idof,1)
         if( Lconstraint.eq.1 ) then ! 01 Jul 14;
         hessian(tdoc+idoc+1:tdoc+idoc+lgeometricaldof,tdof) = hessian(tdoc+idoc+1:tdoc+idoc+lgeometricaldof,tdof)                       & 
                                                             + dBBdmp(idoc+1:idoc+lgeometricaldof,vvol+1,0,1) * dmupfdx(vvol+1,1,idof,1) &
                                                             + dBBdmp(idoc+1:idoc+lgeometricaldof,vvol+1,0,2) * dmupfdx(vvol+1,2,idof,1)
         endif ! end of if( Lconstraint.eq.1 ) then; 20 Jun 14;
         endif

        endif ! end of if( vvol.lt.Mvol-1 ) ; 19 Sep 13;

        if( vvol.eq.Mvol-1 ) then
        !tdof = (vvol+0) * lgeometricaldof + idof
         tdoc = (vvol-1) * lgeometricaldof ! shorthand; 19 Sep 13;         
         idoc = 0
         dessian(tdoc+idoc+1:tdoc+idoc+lgeometricaldof,idof) = dFFdRZ(idoc+1:idoc+lgeometricaldof,vvol+1,0,idof,1)
         if( Lconstraint.eq.1 ) then ! 01 Jul 14;
         dessian(tdoc+idoc+1:tdoc+idoc+lgeometricaldof,idof) = dessian(tdoc+idoc+1:tdoc+idoc+lgeometricaldof,idof)                       & 
                                                             + dBBdmp(idoc+1:idoc+lgeometricaldof,vvol+1,0,1) * dmupfdx(vvol+1,1,idof,1) &
                                                             + dBBdmp(idoc+1:idoc+lgeometricaldof,vvol+1,0,2) * dmupfdx(vvol+1,2,idof,1)
         endif ! end of if( Lconstraint.eq.1 ) then; 20 Jun 14;

         
        endif ! end of if( vvol.lt.Mvol-1 ) ; 19 Sep 13;
        
       enddo ! matches do issym ; 19 Sep 13;
       
      enddo ! matches do irz; 19 Sep 13;
      
     enddo ! matches do ii; 19 Sep 13;
     
    else ! magnetic field calculation in an adjacent volume failed; 20 Feb 13;
     
     FATALMESS(fc02aa, .true., need to provide suitable values for hessian in case of field failure)
     
    endif ! end of if( ImagneticOK(vvol) .and. ImagneticOK(vvol+1) ) ; 20 Feb 13;
    
   enddo ! end of do vvol; 20 Feb 13;
   
  endif ! end of if( LcomputeDerivatives ) ; 19 Sep 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(fc02aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine fc02aa
 
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
