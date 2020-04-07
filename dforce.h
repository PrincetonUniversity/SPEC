!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (&ldquo;global&rdquo; force) ! Given &ldquo;position&rdquo;, ${\bf \xi}$, computes ${\bf F}({\bf \xi})$ and $\nabla_{\bf \xi}{\bf F}$.

!latex \briefly{Calculates ${\bf F}({\bf x})$, where ${\bf x} \equiv \{\mbox{\rm geometry}\} \equiv \{ R_{i,v}, Z_{i,v}\}$ 
!latex          and ${\bf F}\equiv[[p+B^2/2]] + \{\mbox{\rm spectral constraints}\} $, and $\nabla {\bf F}$.}

!latex \calledby{\link{hesian}, 
!latex           \link{newton}, 
!latex           \link{pc00aa}, 
!latex           \link{pc00ab} and 
!latex           \link{xspech}} \\

!latex \calls{\link{packxi}, 
!latex        \link{ma00aa}, 
!latex        \link{matrix}, 
!latex        \link{ma02aa}, 
!latex        \link{lforce}, 
!latex        \link{volume}, 
!latex        \link{packab}, 
!latex        \link{tr00ab}, 
!latex        \link{coords} and 
!latex        \link{brcast}}

!latex \tableofcontents

!latex \subsection{unpacking}

!latex \begin{enumerate}

!latex \item The geometrical degrees of freedom are represented as a vector, ${\bf x} \equiv \{ R_{i,v}, Z_{i,v}\}$, 
!latex       where $i=1,$ \internal{mn} labels the Fourier harmonic and $v=1,$ \internal{Mvol}$-1$ is the interface label.
!latex       This vector is ``unpacked'' using \link{packxi}.
!latex       (Note that \link{packxi} also sets the coordinate axis, i.e. the $R_{i,0}$ and $Z_{i,0}$.)

!latex \end{enumerate}

!latex \subsection{parallelization over volumes}

!latex \begin{enumerate}
!latex \item In each volume, \internal{vvol = 1, Mvol}, 
!latex       \begin{enumerate}
!latex       \item the logical array \internal{ImagneticOK(vvol)} is set to \internal{.false.}
!latex       \item the energy and helicity matrices, \internal{dMA(0:NN,0:NN)}, \internal{dMB(0:NN,0:2)}, etc. are allocated;
!latex       \item the volume-integrated metric arrays, \internal{DToocc}, etc. are allocated;
!latex       \item calls \link{ma00aa} to compute the volume-integrated metric arrays;
!latex       \item calls \link{matrix} to construct the energy and helicity matrices;
!latex       \item calls \link{ma02aa} to solve for the magnetic fields consistent with the appropriate constraints, perhaps by iterating on \link{mp00ac};
!latex       \item calls \link{volume} to compute the volume of the $v$-th region;
!latex       \item calls \link{lforce} to compute $p+B^2/2$ (and the spectral constraints if required) on the inner and outer interfaces;
!latex       \item the derivatives of the force-balance will also be computed if \internal{LComputeDerivatives = 1};
!latex       \end{enumerate}
!latex \item After the parallelization loop over the volumes, \link{brcast} is called to broadcast the required information.
!latex \end{enumerate}

!latex \subsection{broadcasting}

!latex \begin{enumerate}
!latex \item The required quantities are broadcast by \link{brcast}.
!latex \end{enumerate}

!latex \subsection{construction of force}

!latex \begin{enumerate}
!latex \item The force vector, ${\bf F}({\bf x})$, is a combination of the pressure-imbalance Fourier harmonics, $[[p+B^2/2]]_{i,v}$,
!latex       where $i$ labels Fourier harmonic and $v$ is the interface label:
!latex       \be F_{i,v} \equiv \left[ ( p_{v+1}+B^2_{i,v+1}/2 ) - ( p_v + B^2_{i,v}/2 ) \right] \times \internal{BBweight}_i,
!latex       \ee
!latex       where \internal{BBweight(i)} is defined in \link{preset};
!latex       and the spectral condensation constraints, 
!latex       \be F_{i,v} \equiv I_{i,v} \times \inputvar{epsilon} + S_{i,v,1} \times \internal{sweight}_v - S_{i,v+1,0} \times \internal{sweight}_{v+1},
!latex       \ee
!latex       where the spectral condensation constraints, $I_{i,v}$, and the ``star-like'' poloidal angle constraints, $S_{i,v,\pm 1}$,
!latex       are calculated and defined in \link{lforce};
!latex       and the \internal{sweight}$_v$ are defined in \link{preset}.
!latex \end{enumerate}

!latex \subsection{construct derivatives of matrix equation}

!latex \begin{enumerate}
!latex \item Matrix perturbation theory is used to compute the derivatives of the solution, i.e. the Beltrami fields, as the geometry of the 
!latex       interfaces changes:
!latex \end{enumerate}

!latex \subsection{extrapolation: planned redundant}

!latex \begin{enumerate}
!latex \item The extrapolation constraint is $R_{j,1} = R_{j,2} \, \psi_1^{m/2} / \psi_2^{m/2}$.
!latex       Combining this with the regularization factor for the geometry, i.e. $R_{j,i}=\psi_{i}^{m/2} \xi_{j,i}$, we obtain
!latex       \be \xi_{j,1} = R_{j,2} / \psi_2^{m/2}.
!latex       \ee
!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine dforce( NGdof, position, force, LComputeDerivatives)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half, one, two, pi2, pi
  
  use numerical, only : vsmall, small, logtolerance
  
  use fileunits, only : ounit
  
  use inputlist, only : Wmacros, Wdforce, ext, Nvol, Mpol, Ntor, Lrad, tflux, Igeometry, &
                        gamma, adiabatic, pscale, mu, &
                        epsilon, &
                        Lfindzero, &
                        Lconstraint, Lcheck, &
                        Lextrap, Mpol
  
  use cputiming, only : Tdforce
  
  use allglobal, only : ncpu, myid, cpus, pi2nfp, &
                        Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, &
                        Mvol, &
                        Iquad, &                  ! convenience; provided to ma00aa as argument to avoid allocations;
                        iRbc, iZbs, iRbs, iZbc, & ! Fourier harmonics of geometry; vector of independent variables, position, is "unpacked" into iRbc,iZbs;
                        NAdof, &
                        ImagneticOK, &
                        Energy, ForceErr, &
                        YESstellsym, NOTstellsym, &
                        mn, im, in, mns, Ntz, &
                        Ate, Aze, Ato, Azo, & ! only required for debugging;
                        ijreal, ijimag, jireal, jiimag, &
                        efmn, ofmn, cfmn, sfmn, &
                        evmn, odmn, comn, simn, &
                        Nt, Nz, &
                        cosi, sini, & ! FFT workspace;
                        dBdX, &
                       !dMA, dMB, dMC, dMD, dME, dMF, dMG, solution, &
                        dMA, dMB,      dMD,           dMG, solution, &
                        NdMASmax, NdMAS, dMAS, dMDS, idMAS, jdMAS, & ! preconditioning matrix
                        MBpsi,        &
                        dtflux, dpflux, sweight, &
                        mmpp, &
                        Bemn, Bomn, Iomn, Iemn, Somn, Semn, &
                        BBe, IIo, BBo, IIe, & ! these are just used for screen diagnostics;
                        LGdof, &
                        dBBdRZ, dIIdRZ, &
                        vvolume, dvolume, lBBintegral, lABintegral, &
                        Rij, Zij, sg, guvij, iRij, iZij, dRij, dZij, tRij, tZij, & ! Jacobian and metrics; computed in coords;
                        diotadxup, dItGpdxtp, &
                        dFFdRZ, dBBdmp, dmupfdx, hessian, dessian, Lhessianallocated, &
                        BBweight, & ! exponential weight on force-imbalance harmonics;
                        psifactor, Rscale, &
                        lmns, &
                        mn, mne, &
                        DToocc, DToocs, DToosc, DTooss, &
                        TTsscc, TTsscs, TTsssc, TTssss, &
                        TDstcc, TDstcs, TDstsc, TDstss, &
                        TDszcc, TDszcs, TDszsc, TDszss, &
                        DDttcc, DDttcs, DDttsc, DDttss, &
                        DDtzcc, DDtzcs, DDtzsc, DDtzss, &
                        DDzzcc, DDzzcs, DDzzsc, DDzzss, &
                        Tss, Tsc, Dts, Dtc, Dzs, Dzc, &
                        Ttc, Tzc, Tts, Tzs, &
                        dRodR, dRodZ, dZodR, dZodZ, &
                        LILUprecond, NOTMatrixFree, YESMatrixFree, &
                        gaussianabscissae, guvijsave, Lsavedguvij
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in)  :: NGdof               ! dimensions;
  REAL,    intent(in)  :: position(0:NGdof)   ! degrees-of-freedom = internal geometry;
  REAL,    intent(out) :: force(0:NGdof)      ! force;
  LOGICAL, intent(in)  :: LComputeDerivatives ! indicates whether derivatives are to be calculated;
  
  INTEGER              :: NN, IA, ifail, if01adf, vflag, MM, LDA, idgetrf, idgetrs, Lwork
  INTEGER, allocatable :: ipivot(:)
  REAL   , allocatable :: work(:)
  
  INTEGER              :: vvol, innout, ii, jj, irz, issym, iocons, tdoc, idoc, idof, tdof, jdof, ivol, imn, ll, lldof, iidof, jjdof

  INTEGER              :: Lcurvature, ideriv, id, jquad
  
  REAL                 :: lastcpu, lss, lfactor, DDl, MMl

  INTEGER              :: iflag
  REAL                 :: det
  
  REAL                 :: dpsi(1:2)
  REAL   , allocatable :: oBI(:,:), rhs(:) ! original Beltrami-matrix inverse; used to compute derivatives of matrix equation;

  REAL   , allocatable :: dAt(:,:), dAz(:,:), XX(:), YY(:), dBB(:,:), dII(:), dLL(:), dPP(:), length(:), dRR(:,:), dZZ(:,:), constraint(:)

  CHARACTER            :: packorunpack 

#ifdef DEBUG
  INTEGER              :: isymdiff
  REAL                 :: dRZ = 1.0e-05, dvol(-1:+1), evolume, imupf(1:2,-2:2)
  REAL,    allocatable :: oRbc(:,:), oZbs(:,:), oRbs(:,:), oZbc(:,:) ! original geometry;
  REAL,    allocatable :: isolution(:,:)
#endif

  BEGIN(dforce)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  packorunpack = 'U' ! unpack geometrical degrees-of-freedom;
  
  WCALL( dforce, packxi,( NGdof, position(0:NGdof), Mvol, mn, iRbc(1:mn,0:Mvol), iZbs(1:mn,0:Mvol), iRbs(1:mn,0:Mvol), iZbc(1:mn,0:Mvol), packorunpack, LcomputeDerivatives ) )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  
  if( Lcheck.eq.3 .or. Lcheck.eq.4 ) then ! will check volume derivatives;
   
   SALLOCATE( oRbc, (1:mn,0:Mvol), iRbc(1:mn,0:Mvol) )
   SALLOCATE( oZbs, (1:mn,0:Mvol), iZbs(1:mn,0:Mvol) )
   SALLOCATE( oRbs, (1:mn,0:Mvol), iRbs(1:mn,0:Mvol) )
   SALLOCATE( oZbc, (1:mn,0:Mvol), iZbc(1:mn,0:Mvol) )  
   
  endif ! end of if( Lcheck.eq.3 .or. Lcheck.eq.4 ) ;
  
#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  SALLOCATE( dAt       , (1:Ntz,-1:2), zero )
  SALLOCATE( dAz       , (1:Ntz,-1:2), zero )
  SALLOCATE( dBB       , (1:Ntz,-1:2), zero ) ! magnetic field strength (on interfaces) in real space and derivatives;
  SALLOCATE(  XX       , (1:Ntz     ), zero )
  SALLOCATE(  YY       , (1:Ntz     ), zero )
  SALLOCATE( length    , (1:Ntz     ), zero ) ! this is calculated in lforce;
  SALLOCATE( dRR       , (1:Ntz,-1:1), zero )
  SALLOCATE( dZZ       , (1:Ntz,-1:1), zero )
  SALLOCATE( dII       , (1:Ntz     ), zero ) ! spectral constraint;
  SALLOCATE( dLL       , (1:Ntz     ), zero ) ! length   constraint;
  SALLOCATE( dPP       , (1:Ntz     ), zero ) ! poloidal constraint;
  SALLOCATE( constraint, (1:Ntz     ), zero )
 
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( LcomputeDerivatives ) then
#ifdef DEBUG
   FATAL( dforce, .not.allocated(dBBdmp), do not pass go )
#endif
   dBBdmp(1:LGdof,1:Mvol,0:1,1:2) = zero
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  do vvol = 1, Mvol

   LREGION(vvol) ! assigns Lcoordinatesingularity, Lplasmaregion, etc. ;

   ImagneticOK(vvol) = .false.
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   if( myid.ne.modulo(vvol-1,ncpu) ) cycle ! construct Beltrami fields in parallel;
   
   NN = NAdof(vvol) ! shorthand;
   
   if (NOTMatrixFree) then
     SALLOCATE( dMA, (0:NN,0:NN), zero ) ! required for both plasma region and vacuum region;
     SALLOCATE( dMD, (0:NN,0:NN), zero )
   endif

   ! we will the rest even with or without matrix-free
   SALLOCATE( dMB, (0:NN,0: 2), zero )
   SALLOCATE( dMG, (0:NN     ), zero )  
   
   SALLOCATE( solution, (1:NN,-1:2), zero ) ! this will contain the vector potential from the linear solver and its derivatives;
   
   SALLOCATE( MBpsi, (1:NN), zero )
!  SALLOCATE( MEpsi, (1:NN), zero )

   if (LILUprecond) then
     SALLOCATE( dMAS, (1:NdMASmax(vvol)), zero)
     SALLOCATE( dMDS, (1:NdMASmax(vvol)), zero)
     SALLOCATE( idMAS, (1:NN+1), 0)
     SALLOCATE( jdMAS, (1:NdMASmax(vvol)), 0)
   endif ! if we use GMRES and ILU preconditioner
   
   if( LcomputeDerivatives ) then ! allocate some additional memory;

    FATAL(dforce, YESMatrixFree, we do not support Lfindzero=2 with matrix free yet)

    SALLOCATE( oBI, (0:NN,0:NN), zero ) ! inverse of ``original'', i.e. unperturbed, Beltrami matrix;
    SALLOCATE( rhs, (0:NN     ), zero )
    
#ifdef DEBUG
    if( Lcheck.eq.4 ) then
     SALLOCATE( isolution, (1:NN,-2:2), zero )
    endif
#endif

   endif ! end of if( LcomputeDerivatives ) ;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   ll = Lrad(vvol)

   if (Lcoordinatesingularity) then ! different radial dof for Zernike; 02 Jul 19
     lldof = (Lrad(vvol) - mod(Lrad(vvol),2)) / 2
     if (YESMatrixFree) then
       iidof = Mpol + 1
       jjdof = 1
     else
       iidof = mn
       jjdof = mn
     endif
   else
     lldof = Lrad(vvol)
     if (YESMatrixFree) then
       iidof = 1
       jjdof = 1
     else
       iidof = mn
       jjdof = mn
     endif
   end if

   SALLOCATE( DToocc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
   SALLOCATE( TTssss, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
   SALLOCATE( TDstsc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
   SALLOCATE( TDszsc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
   SALLOCATE( DDttcc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
   SALLOCATE( DDtzcc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
   SALLOCATE( DDzzcc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )

   SALLOCATE( guvijsave, (1:Ntz,1:3,1:3,1:Iquad(vvol)), zero)

   SALLOCATE( Tss, (0:lldof,1:mn), zero )
   SALLOCATE( Dtc, (0:lldof,1:mn), zero )
   SALLOCATE( Dzc, (0:lldof,1:mn), zero )
   SALLOCATE( Ttc, (0:lldof,1:mn), zero )
   SALLOCATE( Tzc, (0:lldof,1:mn), zero )

   if (NOTstellsym) then

    SALLOCATE( DToocs, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( DToosc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( DTooss, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )

    SALLOCATE( TTsscc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( TTsscs, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( TTsssc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )

    SALLOCATE( TDstcc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( TDstcs, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( TDstss, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )

    SALLOCATE( TDszcc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( TDszcs, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( TDszss, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )

    SALLOCATE( DDttcs, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( DDttsc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( DDttss, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )

    SALLOCATE( DDtzcs, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( DDtzsc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( DDtzss, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )

    SALLOCATE( DDzzcs, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( DDzzsc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( DDzzss, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )

    SALLOCATE( Tsc, (0:lldof,1:mn), zero )
    SALLOCATE( Dts, (0:lldof,1:mn), zero )
    SALLOCATE( Dzs, (0:lldof,1:mn), zero )
    SALLOCATE( Tts, (0:lldof,1:mn), zero )
    SALLOCATE( Tzs, (0:lldof,1:mn), zero )
   end if !NOTstellsym

   call intghs_workspace_init(vvol)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   dBdX%L = .false. ! first, compute Beltrami fields;
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   ! pre-compute guvij
   ideriv = 0 ; Lcurvature = 1
   WCALL( dforce, compute_guvijsave, (Iquad(vvol), vvol, ideriv, Lcurvature) )
   Lsavedguvij = .true.

  ! we need to construct the preconditioner if needed
   if (LILUprecond) then
      WCALL( dforce, spsint, ( Iquad(vvol), mn, vvol, ll ) )
      WCALL( dforce, spsmat, ( vvol, mn, ll) )
   endif

   if (NOTMatrixFree) then ! construct Beltrami matrix
     WCALL( dforce, ma00aa, ( Iquad(vvol), mn, vvol, ll ) ) ! compute volume integrals of metric elements;
     WCALL( dforce, matrix, ( vvol, mn, ll ) )
   else ! matrix free, so construct something else
     ! we will still need to construct the dMB and dMG matrix
     WCALL( dforce, matrixBG, ( vvol, mn, ll ) )
   endif

   WCALL( dforce, ma02aa, ( vvol, NN ) )
   
   Lsavedguvij = .false.
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef DEBUG
   if( Lcheck.eq.2 ) then
    goto 2000 ! will take no other action except a finite-difference comparison on the derivatives of the rotational-transform wrt mu and dpflux;
   endif
#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   if( LcomputeDerivatives ) then ! compute inverse of Beltrami matrices;
    
    lastcpu = GETTIME
    
    if(Lconstraint .eq. 2) then
      !dMA(1:NN,1:NN) = dMA(1:NN,1:NN) - mu(vvol) * dMD(1:NN,1:NN) ! this corrupts dMA, but dMA is no longer used;
      call DAXPY((NN+1)*(NN+1), -mu(vvol), dMD, 1, dMA, 1) ! BLAS version; 24 Jul 2019
      dMA(0,0)       = zero
      !dMA(1:NN,0)    = zero
      dMA(1:NN,0)    = -matmul(dMD(1:NN,1:NN),solution(1:NN,0))
      ! There is something wrong with the BLAS version for more than one thread, I don't know why, but let's not use it for the moment
      !call DGEMV('N', NN, NN, -one, dMD(1,1), NN+1, solution(1,0), 1, zero, dMA(1,0), 1) ! BLAS version; 24 Jul 2019
      dMA(0,1:NN)    = dMA(1:NN,0) !-matmul(solution(1:NN,0),dMD(1:NN,1:NN))
      IA = NN + 1 + 1
      MM = NN ; LDA = NN+1 ; Lwork = NN+1
      SALLOCATE( ipivot, (0:NN), 0 )
      idgetrf = 1 ; call DGETRF( MM+1, NN+1, dMA(0:LDA-1,0:NN), LDA, ipivot(0:NN), idgetrf )

      cput = GETTIME
      select case( idgetrf ) !                                                                     0123456789012345678
      case(  :-1 ) ;               write(ounit,1010) cput-cpus, myid, vvol, cput-lastcpu, idgetrf, "input error;      "
      case(  0   ) ; if( Wdforce ) write(ounit,1010) cput-cpus, myid, vvol, cput-lastcpu, idgetrf, "success;          "
      case( 1:   ) ;               write(ounit,1010) cput-cpus, myid, vvol, cput-lastcpu, idgetrf, "singular;         "
      case default ;               FATAL( dforce, .true., illegal ifail returned from F07ADF )
      end select    

      ! Using DGETRS instead of inverting matrix ; 22 Jul 09
      !SALLOCATE( work, (1:Lwork), zero )
      !idgetri = 1 ; call DGETRI( NN+1, dMA(0:LDA-1,0:NN), LDA, ipivot(0:NN), work(1:Lwork), Lwork, idgetri )
      !DALLOCATE(work)
      !DALLOCATE(ipivot)

      !cput = GETTIME
      !select case( idgetri ) !                                                                     0123456789012345678
      !case(  :-1 ) ;               write(ounit,1011) cput-cpus, myid, vvol, cput-lastcpu, idgetri, "input error;      "
      !case(  0   ) ; if( Wdforce ) write(ounit,1011) cput-cpus, myid, vvol, cput-lastcpu, idgetri, "success;          "
      !case( 1:   ) ;               write(ounit,1011) cput-cpus, myid, vvol, cput-lastcpu, idgetri, "singular;         "
      !case default ;               FATAL( dforce, .true., illegal ifail returned from F07AJF )
      !end select

      !oBI(0:NN,0:NN) = dMA(0:NN,0:NN)
      call DCOPY((1+NN)*(1+NN), dMA, 1, oBI, 1)            ! BLAS version; 24 Jul 2019

    else ! for LBlinear = T (Linear-solver)

      !dMA(0:NN-1,1:NN) = dMA(1:NN,1:NN) - mu(vvol) * dMD(1:NN,1:NN) ! this corrupts dMA, but dMA is no longer used;
      !dMA(  NN  ,1:NN) = zero
      call DAXPY((NN+1)*(NN+1), -mu(vvol), dMD, 1, dMA, 1) ! BLAS version; 24 Jul 2019

      !dMD(1:NN  ,1:NN) = dMA(0:NN-1,1:NN) ! copy of original matrix; this is used below;
      call DCOPY((NN+1)*(NN+1), dMA, 1, dMD, 1)            ! BLAS version; 24 Jul 2019

      IA = NN + 1
    
      MM = NN ; LDA = NN + 1 ; Lwork = NN
    
      SALLOCATE( ipivot, (1:NN), 0 )
    
      idgetrf = 1 ; call DGETRF( MM, NN, dMA(1,1), LDA, ipivot(1), idgetrf )
    
      cput = GETTIME
      select case( idgetrf ) !                                                                     0123456789012345678
      case(  :-1 ) ;               write(ounit,1010) cput-cpus, myid, vvol, cput-lastcpu, idgetrf, "input error;      "
      case(  0   ) ; if( Wdforce ) write(ounit,1010) cput-cpus, myid, vvol, cput-lastcpu, idgetrf, "success;          "
      case( 1:   ) ;               write(ounit,1010) cput-cpus, myid, vvol, cput-lastcpu, idgetrf, "singular;         "
      case default ;               FATAL( dforce, .true., illegal ifail returned from F07ADF )
      end select
    
      ! Using DGETRS instead of inverting matrix ; 22 Jul 09
      !SALLOCATE( work, (1:Lwork), zero )
    
      !idgetri = 1 ; call DGETRI( NN, dMA(0:LDA-1,1:NN), LDA, ipivot(1:NN), work(1:Lwork), Lwork, idgetri )
      !DALLOCATE(work)
      !DALLOCATE(ipivot)
      
      !cput = GETTIME
      !select case( idgetri ) !                                                                     0123456789012345678
      !case(  :-1 ) ;               write(ounit,1011) cput-cpus, myid, vvol, cput-lastcpu, idgetri, "input error;      "
      !case(  0   ) ; if( Wdforce ) write(ounit,1011) cput-cpus, myid, vvol, cput-lastcpu, idgetri, "success;          "
      !case( 1:   ) ;               write(ounit,1011) cput-cpus, myid, vvol, cput-lastcpu, idgetri, "singular;         "
      !case default ;               FATAL( dforce, .true., illegal ifail returned from F07AJF )
      !end select
    
      oBI(1:NN,1:NN) = dMA(1:NN,1:NN)
    
    endif

1011 format("dforce : ",f10.2," : myid=",i3," ; vvol=",i3," ; called DGETRI ; time=",f10.2,"s ; inverse of Beltrami matrix; idgetrf=",i2," ; ",a18)
1010 format("dforce : ",f10.2," : myid=",i3," ; vvol=",i3," ; called DGETRF ; time=",f10.2,"s ; LU factorization of matrix; idgetrf=",i2," ; ",a18)

    if( Lconstraint.eq.1 ) then ! first, determine how B^2 varies with mu and dpflux;
     
     do iocons = 0, 1 ! labels constraint;
      
      if( vvol.eq.   1 .and. iocons.eq.0 ) cycle ! fixed inner boundary (or coordinate axis); no force-balance constraints;
      if( vvol.eq.Mvol .and. iocons.eq.1 ) cycle ! fixed outer boundary                     ; no force-balance constraints;
      
      do ideriv = 0, 2 ; id = ideriv ! derivatives wrt helicity multiplier and differential poloidal flux;
       
       iflag = 1 ! lforce will only return dAt(1:Ntz,id) and dAz(1:Ntz,id);
       
       WCALL( dforce, lforce, ( vvol, iocons, ideriv, Ntz, dAt(1:Ntz,id), dAz(1:Ntz,id), XX(1:Ntz), YY(1:Ntz), length(1:Ntz), DDl, MMl, iflag ) )
       
       if( ideriv.eq.0 ) then
        dBB(1:Ntz,id) = half * ( &
     dAz(1:Ntz, 0)*dAz(1:Ntz, 0)*guvij(1:Ntz,2,2,0) - two * dAz(1:Ntz, 0)*dAt(1:Ntz, 0)*guvij(1:Ntz,2,3,0) + dAt(1:Ntz, 0)*dAt(1:Ntz, 0)*guvij(1:Ntz,3,3,0) &
     ) / sg(1:Ntz,0)**2
       else
        dBB(1:Ntz,id) = half * ( &
     dAz(1:Ntz,id)*dAz(1:Ntz, 0)*guvij(1:Ntz,2,2,0) - two * dAz(1:Ntz,id)*dAt(1:Ntz, 0)*guvij(1:Ntz,2,3,0) + dAt(1:Ntz,id)*dAt(1:Ntz, 0)*guvij(1:Ntz,3,3,0) &
   + dAz(1:Ntz, 0)*dAz(1:Ntz,id)*guvij(1:Ntz,2,2,0) - two * dAz(1:Ntz, 0)*dAt(1:Ntz,id)*guvij(1:Ntz,2,3,0) + dAt(1:Ntz, 0)*dAt(1:Ntz,id)*guvij(1:Ntz,3,3,0) &
     ) / sg(1:Ntz,0)**2
       endif ! end of if( ideriv.gt.0 ) ;
      
      enddo ! end of do ideriv = 0, 2;

      call tfft( Nt, Nz, dBB(1:Ntz,1), dBB(1:Ntz,2), & ! derivatives of B^2 wrt mu and dpflux;
                 mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), ifail )
         
      ; idoc = 0
      ; dBBdmp(idoc+1:idoc+mn  ,vvol,iocons,1) = efmn(1:mn) * BBweight(1:mn) ! pressure;
      ; dBBdmp(idoc+1:idoc+mn  ,vvol,iocons,2) = cfmn(1:mn) * BBweight(1:mn) ! pressure;
      ; idoc = idoc + mn   ! even;
      if( Igeometry.ge.3 ) then ! add spectral constraints; spectral constraints do not depend on mu or dpflux;
       ;idoc = idoc + mn-1 ! oddd;
      endif ! end of if( Igeometry.ge.3) ;
      if( NOTstellsym ) then
       ;dBBdmp(idoc+1:idoc+mn-1,vvol,iocons,1) = ofmn(2:mn) * BBweight(2:mn) ! pressure;
       ;dBBdmp(idoc+1:idoc+mn-1,vvol,iocons,2) = sfmn(2:mn) * BBweight(2:mn) ! pressure;
       ;idoc = idoc + mn-1 ! oddd;
       if( Igeometry.ge.3 ) then ! add spectral constraints;
        idoc = idoc + mn   ! even;
       endif ! end of if( Igeometry.ge.3) ;
      endif ! end of if( NOTstellsym) ;
      
     enddo ! end of do iocons;
     
    endif ! end of if( Lconstraint.eq.1 ) ;
    
   endif ! end of if( LcomputeDerivatives ) ;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   vflag = 1
   WCALL( dforce, volume, ( vvol, vflag ) ) ! compute volume;
      
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   do iocons = 0, 1 ! construct field magnitude on inner and outer interfaces; inside do vvol;
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
    if( vvol.eq.1    .and. iocons.eq.0 ) cycle ! fixed inner boundary (or coordinate axis);
    if( vvol.eq.Mvol .and. iocons.eq.1 ) cycle ! fixed outer boundary                     ; there are no constraints at outer boundary;
    
    ideriv = 0 ; id = ideriv

    iflag = 0 ! dAt, dAz, XX & YY are returned by lforce; Bemn(1:mn,vvol,iocons), Iomn(1:mn,vvol) etc. are returned through global;
    
    WCALL( dforce, lforce, ( vvol, iocons, ideriv, Ntz, dAt(1:Ntz,id), dAz(1:Ntz,id), XX(1:Ntz), YY(1:Ntz), length(1:Ntz), DDl, MMl, iflag ) )
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
   enddo ! end of do iocons = 0, 1;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! still inside do vvol = 1, Mvol loop;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   if( LComputeDerivatives ) then

    dBdX%L = .true. ! will need derivatives;
    
    idof = 0 ! labels degree of freedom; local to interface;
    
    do ii = 1, mn ! loop over deformations in Fourier harmonics; inside do vvol;
     
     dBdX%ii = ii ! controls construction of derivatives in subroutines called below;     

     do irz = 0, 1 ! loop over deformations in R and Z; inside do vvol; inside do ii;
      
      if( irz.eq.1 .and. Igeometry.lt.3 ) cycle ! no dependence on Z;
      
      dBdX%irz = irz ! controls construction of derivatives;

      do issym = 0, 1 ! loop over stellarator and non-stellarator symmetric terms;
       
       if( issym.eq.1 .and. YESstellsym ) cycle ! no dependence on non-stellarator symmetric harmonics;
       
       if( ii.eq.1 .and. irz.eq.1 .and. issym.eq.0 ) cycle ! no dependence on Zbs_{m=0,n=0};
       if( ii.eq.1 .and. irz.eq.0 .and. issym.eq.1 ) cycle ! no dependence on Rbs_{m=0,n=0};
       
       dBdX%issym = issym ! controls construction of derivatives;
       
       idof = idof + 1 ! this labels the degree-of-freedom that the derivative is taken with respect to; this is outside do innout;
       
#ifdef DEBUG
       FATAL( dforce, idof.gt.LGdof, illegal degree-of-freedom index constructing derivatives ) ! this can be deleted;
#endif
       
       do innout = 0, 1 ! loop over deformations to inner and outer interface; inside do vvol; inside do ii; inside do irz;
        
        if( vvol.eq.1    .and. innout.eq.0 ) cycle ! no degrees of freedom at coordinate axis / fixed inner boundary;
        if( vvol.eq.Mvol .and. innout.eq.1 ) cycle ! no degress of freedom                      fixed outer boundary; for linearized displacement;
        
        dBdX%innout = innout

        !WCALL( dforce, ma00aa,( Iquad(vvol), mn, vvol, ll ) ) ! compute volume integrals of metric elements;
        
        !WCALL( dforce, matrix,( vvol, mn, ll ) ) ! construct Beltrami matrices;
        
        !dpsi(1:2) = (/ dtflux(vvol), dpflux(vvol) /) ! local enclosed toroidal and poloidal fluxes;
        
        ! this is the original without BLAS
        !rhs(0)    =   half*sum(solution(1:NN,0)*matmul(dMD(1:NN,1:NN),solution(1:NN,0)))
        !rhs(1:NN) = - matmul( dMB(1:NN,1:2 )                            , dpsi(1:2)        ) &
        !            - matmul( dMA(1:NN,1:NN) - mu(vvol) * dMD(1:NN,1:NN), solution(1:NN,0) )

        ! BLAS version 17 Jul 2019
        !call DGEMV('N',NN,NN, one,dMD(1,1),NN+1,solution(1,0),1,zero,rhs(1),1)
        !rhs(0) = half * sum(solution(1:NN,0) * rhs(1:NN))
        !call DGEMV('N',NN,NN,-one,dMA(1,1),NN+1,solution(1,0),1,mu(vvol),rhs(1),1)
        !rhs(1:NN) = rhs(1:NN) - matmul( dMB(1:NN,1:2 ), dpsi(1:2) )

        WCALL( dforce, intghs, ( Iquad(vvol), mn, vvol, ll, 0 ) )

        WCALL( dforce, mtrxhs, ( vvol, mn, ll, dMA(0:NN,0), dMD(0:NN,0), 0) )

        rhs(1:NN) = -dMA(1:NN,0)

        if (Lconstraint .eq. 2) then
          SALLOCATE( work, (1:NN+1), zero )
          work(1:NN+1) = rhs(0:NN)
          !work(1:NN+1)  =  matmul( oBI(0:NN,0:NN), rhs(0:NN)) ! original version
          call DGETRS('N',NN+1,1,oBI(0,0),NN+1,ipivot(0:NN),work(1),NN+1,idgetrs) ! Change to DGETRS; 22 Jul 19
          solution(1:NN,-1) = work(2:NN+1)
          DALLOCATE(work)
        else
          !solution(1:NN,-1) = matmul( oBI(1:NN,1:NN), rhs(1:NN) ) ! original version
          solution(1:NN,-1) = rhs(1:NN)
          call DGETRS('N',NN,1,oBI(1,1),NN+1,ipivot(1:NN),solution(1,-1),NN,idgetrs) ! Change to DGETRS; 22 Jul 19
        endif

        ideriv = -1 ; dpsi(1:2) = (/ dtflux(vvol), dpflux(vvol) /) ! these are also used below;
        
        packorunpack = 'U'
        
        WCALL( dforce, packab,( packorunpack, vvol, NN,  solution(1:NN,-1), ideriv ) ) ! derivatives placed in Ate(vvol,ideriv,1:mn)%s(0:Lrad),

#ifdef DEBUG
        FATAL( dforce, vvol-1+innout.gt.Mvol, psifactor needs attention )
#endif

        lfactor = psifactor(ii,vvol-1+innout) ! this "pre-conditions" the geometrical degrees-of-freedom;
        

        dmupfdx(vvol,1,idof,innout) = zero
        dmupfdx(vvol,2,idof,innout) = zero
        
        if( Lconstraint.eq.1 .or. ( Lvacuumregion .and. Lconstraint.ge.0 ) ) then ! will need to accommodate constraints;
         
         if(                     Lconstraint.eq.1 ) then
          iflag = -1 ; WCALL( dforce, tr00ab, ( vvol, mn, lmns, Nt, Nz, iflag, diotadxup(0:1,-1:2,vvol) ) ) ! compute d(transform)/dx;
         endif
         
         if( Lvacuumregion .and. Lconstraint.ge.0 ) then
          iflag = -1 ; WCALL( dforce, curent, ( vvol, mn,       Nt, Nz, iflag, dItGpdxtp(0:1,-1:2,vvol) ) ) ! compute d(Itor,Gpol)/dx;
         endif
         
         if( Lplasmaregion .and. Lconstraint.eq.1 ) then
          
          if( Lcoordinatesingularity ) then ! solution does not depend on dpflux, and only outer transform is a constraint;
           
           det = diotadxup(1,1,vvol)
           FATAL( dforce, abs(det).lt.small, error computing derivatives of mu          wrt geometry at fixed transform )
           
           dmupfdx(vvol,1,idof,innout) = - lfactor * (                                                                        diotadxup(1,-1,vvol) ) / det
           dmupfdx(vvol,2,idof,innout) =   zero
           
          else ! if( .not.Lcoordinatesingularity ) ;
           
           det = diotadxup(0,1,vvol) * diotadxup(1,2,vvol) - diotadxup(0,2,vvol) * diotadxup(1,1,vvol)
           FATAL( dforce, abs(det).lt.small, error computing derivatives of mu & dpflux wrt geometry at fixed transform )
           
           dmupfdx(vvol,1,idof,innout) = - lfactor * ( + diotadxup(1, 2,vvol) * diotadxup(0,-1,vvol) - diotadxup(0, 2,vvol) * diotadxup(1,-1,vvol) ) / det
           dmupfdx(vvol,2,idof,innout) = - lfactor * ( - diotadxup(1, 1,vvol) * diotadxup(0,-1,vvol) + diotadxup(0, 1,vvol) * diotadxup(1,-1,vvol) ) / det
           
          endif ! end of if( Lcoordinatesingularity ) ;
          
         endif ! end of if( Lplasmaregion .and. Lconstraint.eq.1 ) ;
         
         if( Lvacuumregion ) then
          
          if    ( Lconstraint.eq.0 ) then ! THIS NEEDS ATTENTION;

           det = dItGpdxtp(0,1,vvol) * dItGpdxtp(1,2,vvol) - dItGpdxtp(0,2,vvol) * dItGpdxtp(1,1,vvol)
           FATAL( dforce, abs(det).lt.small, error computing derivatives of dtflux & dpflux wrt geometry at fixed Itor and Gpol )
           
           dmupfdx(vvol,1,idof,innout) = - lfactor * ( + dItGpdxtp(1, 2,vvol) * dItGpdxtp(0,-1,vvol) - dItGpdxtp(0, 2,vvol) * dItGpdxtp(1,-1,vvol) ) / det
           dmupfdx(vvol,2,idof,innout) = - lfactor * ( - dItGpdxtp(1, 1,vvol) * dItGpdxtp(0,-1,vvol) + dItGpdxtp(0, 1,vvol) * dItGpdxtp(1,-1,vvol) ) / det
           
          elseif( Lconstraint.eq.1 ) then

           det = diotadxup(0,1,vvol) * dItGpdxtp(1,2,vvol) - diotadxup(0,2,vvol) * dItGpdxtp(1,1,vvol)
           FATAL( dforce, abs(det).lt.small, error computing derivatives of dtflux & dpflux wrt geometry at fixed Itor and Gpol )
           
           dmupfdx(vvol,1,idof,innout) = - lfactor * ( + dItGpdxtp(1, 2,vvol) * diotadxup(0,-1,vvol) - diotadxup(0, 2,vvol) * dItGpdxtp(1,-1,vvol) ) / det
           dmupfdx(vvol,2,idof,innout) = - lfactor * ( - dItGpdxtp(1, 1,vvol) * diotadxup(0,-1,vvol) + diotadxup(0, 1,vvol) * dItGpdxtp(1,-1,vvol) ) / det

          endif
          
         endif ! end of if( Lvacuumregion ) ;
         
        endif ! end of if( Lconstraint.eq.1 .or. ( Lvacuumregion .and. Lconstraint.ge.0 ) then;
        

#ifdef DEBUG

        if( Lcheck.eq.4 ) then ! check derivatives of field;
         
         dBdX%L = .false.
         
         do isymdiff = -2, 2 ! symmetric fourth-order, finite-difference used to approximate derivatives;
          
          if( isymdiff.eq.0 ) cycle
          
          iRbc(1:mn,0:Mvol) = oRbc(1:mn,0:Mvol)
          iZbs(1:mn,0:Mvol) = oZbs(1:mn,0:Mvol)
          iRbs(1:mn,0:Mvol) = oRbs(1:mn,0:Mvol)
          iZbc(1:mn,0:Mvol) = oZbc(1:mn,0:Mvol)
          
          if( issym.eq.0 .and. irz.eq.0 ) iRbc(ii,vvol-1+innout) = iRbc(ii,vvol-1+innout) + dRZ * isymdiff ! perturb geometry;
          if( issym.eq.0 .and. irz.eq.1 ) iZbs(ii,vvol-1+innout) = iZbs(ii,vvol-1+innout) + dRZ * isymdiff ! perturb geometry;
          if( issym.eq.1 .and. irz.eq.0 ) iRbs(ii,vvol-1+innout) = iRbs(ii,vvol-1+innout) + dRZ * isymdiff ! perturb geometry;
          if( issym.eq.1 .and. irz.eq.1 ) iZbc(ii,vvol-1+innout) = iZbc(ii,vvol-1+innout) + dRZ * isymdiff ! perturb geometry;
          
          WCALL( dforce, ma00aa, ( Iquad(vvol), mn, vvol, ll ) )
          
          WCALL( dforce, matrix, ( vvol, mn, ll ) )

          WCALL( dforce, ma02aa, ( vvol, NN ) ) ! this may or may not iterate on mu and dpflux to enforce transform constraints;
          
          if( Lplasmaregion ) then
          imupf(1:2,isymdiff) = (/     mu(vvol), dpflux(vvol) /) ! mu     and dpflux are computed for the perturbed geometry by ma02aa/mp00ac if Lconstraint=1;
          else ! if( Lvacuumregion ) ;
          imupf(1:2,isymdiff) = (/ dtflux(vvol), dpflux(vvol) /) ! dtflux and dpflux are computed for the perturbed geometry by ma02aa/mp00ac if Lconstraint=1;
          endif
          
          isolution(1:NN,isymdiff) = solution(1:NN,0) ! solution is computed in mp00ac, which is called by ma02aa;
          
         enddo ! end of do isymdiff;
         
         isolution(1:NN,0) = ( - 1 * isolution(1:NN, 2) + 8 * isolution(1:NN, 1) - 8 * isolution(1:NN,-1) + 1 * isolution(1:NN,-2) ) / ( 12 * dRZ )
         imupf(1:2,0)      = ( - 1 * imupf(1:2, 2)      + 8 * imupf(1:2, 1)      - 8 * imupf(1:2,-1)      + 1 * imupf(1:2,-2)      ) / ( 12 * dRZ )
         
          solution(1:NN,-1) = abs(  solution(1:NN,-1) )
         isolution(1:NN, 0) = abs( isolution(1:NN, 0) )
        
!        ifail = 0 ; call M01CAF(  solution(1:NN,-1), 1, NN, 'D', ifail ) ! sorting screen output; this corrupts;
!        ifail = 0 ; call M01CAF( isolution(1:NN, 0), 1, NN, 'D', ifail ) ! sorting screen output; this corrupts;
         
         ifail = 0 ; call dlasrt( 'D', NN,  solution(1:NN,-1), ifail ) ! sorting screen output; this corrupts;
         ifail = 0 ; call dlasrt( 'D', NN, isolution(1:NN, 0), ifail ) ! sorting screen output; this corrupts;        
         
         cput = GETTIME

        !select case( Lconstraint )
        !case( 0 )
        ! write(ounit,3002) cput-cpus, myid, vvol, im(ii), in(ii), irz, issym, innout,  solution(1:min(NN,16),-1)
        ! write(ounit,3002) cput-cpus, myid, vvol, im(ii), in(ii), irz, issym, innout, isolution(1:min(NN,16), 0)
        !case( 1 )
          write(ounit,3003)
          write(ounit,3003) cput-cpus, myid, vvol, im(ii), in(ii), irz, issym, innout, "finite-diff", imupf(1:2,0)
          write(ounit,3003) cput-cpus, myid, vvol, im(ii), in(ii), irz, issym, innout, "analytic   ", dmupfdx(vvol,1:2,idof,innout) / lfactor
        !case default
        ! FATAL( dforce, .true., Lconstraint not supported for Lcheck=4 )
        !end select ! end of select case( Lconstraint ) ;
         
!3002     format("dforce : ",f10.2," : ",:,"myid=",i3," ; vvol=",i2," ; (",i2,",",i3," ) ; irz=",i1," ; issym=",i1," ; innout=",i1," ; dxi=",99es10.02)
3003     format("dforce : ",f10.2," : ",:,"myid=",i3," ; vvol=",i2," ; (",i2,",",i3," ) ; irz=",i1," ; issym=",i1," ; innout=",i1," ; ",a11,&
      " : dmupf=",2f11.05" ;")
         
         dBdX%L = .true.
         
         iRbc(1:mn,0:Mvol) = oRbc(1:mn,0:Mvol)
         iZbs(1:mn,0:Mvol) = oZbs(1:mn,0:Mvol)
         iRbs(1:mn,0:Mvol) = oRbs(1:mn,0:Mvol)
         iZbc(1:mn,0:Mvol) = oZbc(1:mn,0:Mvol)
         
        endif ! end of if( Lcheck.eq.4 ) ;

#endif
        
        vflag = 1 ! this flag instructs volume to continue even if the volume is invalid;
        WCALL( dforce, volume, ( vvol, vflag ) ) ! compute derivative of volume; wrt to harmonic described by dBdX structure;
        
#ifdef DEBUG

        if( Lcheck.eq.3 ) then
         
         dvol(0) = dvolume
         
         cput = GETTIME
         write(ounit,1001) cput-cpus, myid, vvol, im(ii), in(ii), irz, issym, innout, "analytic", dvolume
         
1001     format("dforce : ",f10.2," : myid=",i3," ; vvol=",i3," ; (",i3," ,",i3,") ; irz=",i2," ; issym=",i2," ; innout=",i2,&
       " : ",a8," ; dvolume=",f23.15," ;",:," error=",es13.5," ;")
         
         dBdX%L = .false. ! instruct volume to not calculate derivatives;
         
         do isymdiff = -1, 1, 2 ! symmetric finite-difference estimate of derivative of volume wrt geometrical degree-of-freedom;
          
          if( dBdX%issym.eq.0 ) then !     stellarator symmetric harmonics;
           if( dBdX%irz.eq.0 ) iRbc(dBdX%ii,vvol-1+innout) = oRbc(dBdX%ii,vvol-1+innout) + isymdiff * dRZ * half
           if( dBdX%irz.eq.1 ) iZbs(dBdX%ii,vvol-1+innout) = oZbs(dBdX%ii,vvol-1+innout) + isymdiff * dRZ * half
          else                                        ! non-stellarator symmetric harmonics;
           if( dBdX%irz.eq.0 ) iRbs(dBdX%ii,vvol-1+innout) = oRbs(dBdX%ii,vvol-1+innout) + isymdiff * dRZ * half
           if( dBdX%irz.eq.1 ) iZbc(dBdX%ii,vvol-1+innout) = oZbc(dBdX%ii,vvol-1+innout) + isymdiff * dRZ * half
          endif
          
          vflag = 1 ! this flag instructs volume to continue even if the volume is invalid;
          WCALL( dforce, volume, ( vvol, vflag ) ) ! compute volume; this corrupts calculation of dvolume;
          
          dvol(isymdiff) = vvolume(vvol)
          
         enddo ! end of do isymdiff;
         
         evolume = abs( ( dvol(+1)-dvol(-1) ) / dRZ - dvol(0) ) ! error in finite-difference calculation and analytic derivative;
         
         cput = GETTIME
         write(ounit,1001) cput-cpus, myid, vvol, im(ii), in(ii), irz, issym, innout, "finite-d", ( dvol(+1)-dvol(-1) ) / dRZ, evolume
         
         FATAL( dforce, evolume.gt.dRZ, unacceptable error in volume derivative )

         iRbc(1:mn,0:Mvol) = oRbc(1:mn,0:Mvol)
         iZbs(1:mn,0:Mvol) = oZbs(1:mn,0:Mvol)
         iRbs(1:mn,0:Mvol) = oRbs(1:mn,0:Mvol)
         iZbc(1:mn,0:Mvol) = oZbc(1:mn,0:Mvol)
         
         dBdX%L = .true.
         
         dvolume = dvol(0)
        
        endif ! end of if( Lcheck.eq.3 ) ;
        
#endif
        
        do iocons = 0, 1
         
         if( vvol.eq.   1 .and. iocons.eq.0 ) cycle ! fixed inner boundary (or coordinate axis); no constraints;
         if( vvol.eq.Mvol .and. iocons.eq.1 ) cycle ! fixed outer boundary                     ; no constraints;
    
         ideriv = 0 ; id = ideriv ; iflag = 0 ! why does lforce need to be re-called; need to determine which quantity (if any) has been corrupted;

         WCALL( dforce, lforce, ( vvol, iocons, ideriv, Ntz, dAt(1:Ntz,id), dAz(1:Ntz,id), XX(1:Ntz), YY(1:Ntz), length(1:Ntz), DDl, MMl, iflag ) )

         dBB(1:Ntz,0) = half * ( &
      dAz(1:Ntz, 0)*dAz(1:Ntz, 0)*guvij(1:Ntz,2,2,0) - two * dAz(1:Ntz, 0)*dAt(1:Ntz, 0)*guvij(1:Ntz,2,3,0) + dAt(1:Ntz, 0)*dAt(1:Ntz, 0)*guvij(1:Ntz,3,3,0) &
      ) / sg(1:Ntz,0)**2

         ideriv = -1 ; id = ideriv ; iflag = 1 ! compute derivatives of magnetic field;
         
         WCALL( dforce, lforce, ( vvol, iocons, ideriv, Ntz, dAt(1:Ntz,id), dAz(1:Ntz,id), XX(1:Ntz), YY(1:Ntz), length(1:Ntz), DDl, MMl, iflag ) )
         
         lss = two * iocons - one ; Lcurvature = 4
         WCALL( dforce, coords, ( vvol, lss, Lcurvature, Ntz, mn ) ) ! get coordinate metrics and their derivatives wrt Rj, Zj on interface;
         
         dBB(1:Ntz,id) = half * ( &
      dAz(1:Ntz,id)*dAz(1:Ntz, 0)*guvij(1:Ntz,2,2,0) - two * dAz(1:Ntz,id)*dAt(1:Ntz, 0)*guvij(1:Ntz,2,3,0) + dAt(1:Ntz,id)*dAt(1:Ntz, 0)*guvij(1:Ntz,3,3,0) &
    + dAz(1:Ntz, 0)*dAz(1:Ntz,id)*guvij(1:Ntz,2,2,0) - two * dAz(1:Ntz, 0)*dAt(1:Ntz,id)*guvij(1:Ntz,2,3,0) + dAt(1:Ntz, 0)*dAt(1:Ntz,id)*guvij(1:Ntz,3,3,0) &
    + dAz(1:Ntz, 0)*dAz(1:Ntz, 0)*guvij(1:Ntz,2,2,1) - two * dAz(1:Ntz, 0)*dAt(1:Ntz, 0)*guvij(1:Ntz,2,3,1) + dAt(1:Ntz, 0)*dAt(1:Ntz, 0)*guvij(1:Ntz,3,3,1) &
    ) / sg(1:Ntz,0)**2                                                                                                                                       &
    - dBB(1:Ntz,0) * two * sg(1:Ntz,1) / sg(1:Ntz,0)

         FATAL( dforce, vvolume(vvol).lt.small, shall divide by vvolume(vvol)**(gamma+one) )

         ijreal(1:Ntz) = - adiabatic(vvol) * pscale * gamma * dvolume / vvolume(vvol)**(gamma+one) + dBB(1:Ntz,-1) ! derivatives of force wrt geometry;

         dLL(1:Ntz) = zero ! either no spectral constraint, or not the appropriate interface;
         dPP(1:Ntz) = zero ! either no spectral constraint, or not the appropriate interface;
         
         if( Igeometry.ge.3 ) then ! spectral constraints are only required in toroidal or extended-cylindrical geometry;
          
          if( innout.eq.1 .and. iocons.eq.1 ) then ! include derivatives of spectral constraints;
#ifdef DEBUG
             FATAL( dforce, abs(DDl).lt.small, divide by zero on spectral constraint )
#endif
           if( issym.eq.0 ) then ! take derivatives wrt Rbc and Zbs;
            if( irz.eq.0 ) then ! take derivative wrt Rbc;
             dII(1:Ntz) = - im(ii) * sini(1:Ntz,ii) * ( XX(1:Ntz) - MMl * iRij(1:Ntz,vvol) ) &
                          - two * ( mmpp(ii) - MMl ) * iRbc(ii,vvol) * ( Rij(1:Ntz,2,0) * iRij(1:Ntz,vvol) + Zij(1:Ntz,2,0) * iZij(1:Ntz,vvol) ) / DDl &
                          + Rij(1:Ntz,2,0) * ( mmpp(ii) - MMl ) * cosi(1:Ntz,ii)
            else ! take derivative wrt Zbs;
             dII(1:Ntz) = + im(ii) * cosi(1:Ntz,ii) * ( YY(1:Ntz) - MMl * iZij(1:Ntz,vvol) ) &
                          - two * ( mmpp(ii) - MMl ) * iZbs(ii,vvol) * ( Rij(1:Ntz,2,0) * iRij(1:Ntz,vvol) + Zij(1:Ntz,2,0) * iZij(1:Ntz,vvol) ) / DDl &
                          + Zij(1:Ntz,2,0) * ( mmpp(ii) - MMl ) * sini(1:Ntz,ii)
            endif ! end of if( irz.eq.0 ) ;
           else                  ! take derivatives wrt Rbs and Zbc;
            if( irz.eq.0 ) then 
             dII(1:Ntz) = + im(ii) * cosi(1:Ntz,ii) * ( XX(1:Ntz) - MMl * iRij(1:Ntz,vvol) ) &
                          - two * ( mmpp(ii) - MMl ) * iRbs(ii,vvol) * ( Rij(1:Ntz,2,0) * iRij(1:Ntz,vvol) + Zij(1:Ntz,2,0) * iZij(1:Ntz,vvol) ) / DDl &
                          + Rij(1:Ntz,2,0) * ( mmpp(ii) - MMl ) * sini(1:Ntz,ii)
            else                
             dII(1:Ntz) = - im(ii) * sini(1:Ntz,ii) * ( YY(1:Ntz) - MMl * iZij(1:Ntz,vvol) ) &
                          - two * ( mmpp(ii) - MMl ) * iZbc(ii,vvol) * ( Rij(1:Ntz,2,0) * iRij(1:Ntz,vvol) + Zij(1:Ntz,2,0) * iZij(1:Ntz,vvol) ) / DDl &
                          + Zij(1:Ntz,2,0) * ( mmpp(ii) - MMl ) * cosi(1:Ntz,ii)
            endif
           endif

          else

           dII(1:Ntz) = zero ! either no spectral constraint, or not the appropriate interface;

          endif ! end of if( innout.eq.1 .and. iocons.eq.1 ) ;

          constraint(1:Ntz) = + ( dRij(1:Ntz,vvol) * tRij(1:Ntz,vvol-1+iocons) + dZij(1:Ntz,vvol) * tZij(1:Ntz,vvol-1+iocons) ) / length(1:Ntz)
         
          if( iocons.eq.0 ) then ! take derivatives of constraints at inner boundary;
           
           if( innout.eq.0 ) then ! derivative wrt inner boundary coefficient;
           !write(ounit,'("dforce : " 10x " : A ; vvol="i3" ; iocons="i2" ; innout="i2" ;")') vvol, iocons, innout
            if( issym.eq.0 ) then ! take derivatives wrt Rbc and Zbs;
             if( irz.eq.0 ) then ; dLL(1:Ntz) = + ( - cosi(1:Ntz,ii) * tRij(1:Ntz,vvol-1) - dRij(1:Ntz,vvol) * im(ii) * sini(1:Ntz,ii) ) / length(1:Ntz) & 
                                                + constraint(1:Ntz) * dRij(1:Ntz,vvol) * cosi(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
             else                ; dLL(1:Ntz) = + ( - sini(1:Ntz,ii) * tZij(1:Ntz,vvol-1) + dZij(1:Ntz,vvol) * im(ii) * cosi(1:Ntz,ii) ) / length(1:Ntz) & 
                                                + constraint(1:Ntz) * dZij(1:Ntz,vvol) * sini(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
             endif
            else ! if issym.eq.1 ; take derivatives wrt Rbs and Zbc;
             if( irz.eq.0 ) then ; dLL(1:Ntz) = + ( - sini(1:Ntz,ii) * tRij(1:Ntz,vvol-1) + dRij(1:Ntz,vvol) * im(ii) * cosi(1:Ntz,ii) ) / length(1:Ntz) & 
                                                + constraint(1:Ntz) * dRij(1:Ntz,vvol) * sini(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
             else                ; dLL(1:Ntz) = + ( - cosi(1:Ntz,ii) * tZij(1:Ntz,vvol-1) - dZij(1:Ntz,vvol) * im(ii) * sini(1:Ntz,ii) ) / length(1:Ntz) & 
                                                + constraint(1:Ntz) * dZij(1:Ntz,vvol) * cosi(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
             endif
            endif
           else ! if innout.eq.1 ; derivative wrt outer boundary coefficient;
           !write(ounit,'("dforce : " 10x " : B ; vvol="i3" ; iocons="i2" ; innout="i2" ;")') vvol, iocons, innout
            if( issym.eq.0 ) then ! take derivatives wrt Rbc and Zbs;
             if( irz.eq.0 ) then ; dLL(1:Ntz) = + ( + cosi(1:Ntz,ii) * tRij(1:Ntz,vvol-1)                                              ) / length(1:Ntz) & 
                                                - constraint(1:Ntz) * dRij(1:Ntz,vvol) * cosi(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
             else                ; dLL(1:Ntz) = + ( + sini(1:Ntz,ii) * tZij(1:Ntz,vvol-1)                                              ) / length(1:Ntz) & 
                                                - constraint(1:Ntz) * dZij(1:Ntz,vvol) * sini(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
             endif
            else ! if issym.eq.1 ; take derivatives wrt Rbs and Zbc;
             if( irz.eq.0 ) then ; dLL(1:Ntz) = + ( + sini(1:Ntz,ii) * tRij(1:Ntz,vvol-1)                                              ) / length(1:Ntz) & 
                                                - constraint(1:Ntz) * dRij(1:Ntz,vvol) * sini(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
             else                ; dLL(1:Ntz) = + ( + cosi(1:Ntz,ii) * tZij(1:Ntz,vvol-1)                                              ) / length(1:Ntz) & 
                                                - constraint(1:Ntz) * dZij(1:Ntz,vvol) * cosi(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
             endif
            endif
           endif

          else ! if iocons.eq.1 ; take derivatives of constraints at outer boundary;

           if( innout.eq.0 ) then ! derivative wrt inner boundary coefficient;
           !write(ounit,'("dforce : " 10x " : C ; vvol="i3" ; iocons="i2" ; innout="i2" ;")') vvol, iocons, innout
            if( issym.eq.0 ) then ! take derivatives wrt Rbc and Zbs;
             if( irz.eq.0 ) then ; dLL(1:Ntz) = + ( - cosi(1:Ntz,ii) * tRij(1:Ntz,vvol  )                                              ) / length(1:Ntz) & 
                                                + constraint(1:Ntz) * dRij(1:Ntz,vvol) * cosi(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
             else                ; dLL(1:Ntz) = + ( - sini(1:Ntz,ii) * tZij(1:Ntz,vvol  )                                              ) / length(1:Ntz) & 
                                                + constraint(1:Ntz) * dZij(1:Ntz,vvol) * sini(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
             endif
            else                  ! take derivatives wrt Rbs and Zbc;
             if( irz.eq.0 ) then ; dLL(1:Ntz) = + ( - sini(1:Ntz,ii) * tRij(1:Ntz,vvol  )                                              ) / length(1:Ntz) & 
                                                + constraint(1:Ntz) * dRij(1:Ntz,vvol) * sini(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
             else                ; dLL(1:Ntz) = + ( - cosi(1:Ntz,ii) * tZij(1:Ntz,vvol  )                                              ) / length(1:Ntz) & 
                                                + constraint(1:Ntz) * dZij(1:Ntz,vvol) * cosi(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
             endif
            endif
           else ! if innout.eq.1 ; derivative wrt outer boundary coefficient;
!#ifdef AXIS
            if( Igeometry.eq.3 .and. vvol.eq.1 ) then ! need to accomodate derivatives of coordinate axis;
!#else
!            if( Igeometry.eq.3 .and. vvol.lt.1 ) then ! need to accomodate derivatives of coordinate axis;
!#endif
            !write(ounit,'("dforce : " 10x " : dRodR(1: ,0,"i2")=",99es11.3)') ii, dRodR(1:20,0,ii)
            !write(ounit,'("dforce : " 10x " : dRodR(1: ,1,"i2")=",99es11.3)') ii, dRodR(1:20,1,ii)
            !write(ounit,'("dforce : " 10x " : dRodZ(1: ,0,"i2")=",99es11.3)') ii, dRodZ(1:20,0,ii)
            !write(ounit,'("dforce : " 10x " : dRodZ(1: ,1,"i2")=",99es11.3)') ii, dRodZ(1:20,1,ii)
            !write(ounit,'("dforce : " 10x " : dZodR(1: ,0,"i2")=",99es11.3)') ii, dZodR(1:20,0,ii)
            !write(ounit,'("dforce : " 10x " : dZodR(1: ,1,"i2")=",99es11.3)') ii, dZodR(1:20,1,ii)
            !write(ounit,'("dforce : " 10x " : dZodZ(1: ,0,"i2")=",99es11.3)') ii, dZodZ(1:20,0,ii)
            !write(ounit,'("dforce : " 10x " : dZodZ(1: ,1,"i2")=",99es11.3)') ii, dZodZ(1:20,1,ii)
             if( issym.eq.0 ) then ! take derivatives wrt Rbc and Zbs;
              if( irz.eq.0 ) then ; dLL(1:Ntz) = ( &   ! d/dRbc ;
                                                   + ( cosi(1:Ntz,ii) - dRodR(1:Ntz,0,ii) ) * tRij(1:Ntz,vvol) - dRij(1:Ntz,vvol) * im(ii) * sini(1:Ntz,ii) &
                                                   + (                - dZodR(1:Ntz,0,ii) ) * tZij(1:Ntz,vvol) &
                                                   - constraint(1:Ntz) &
                                                   * ( dRij(1:Ntz,vvol) * ( cosi(1:Ntz,ii) - dRodR(1:Ntz,0,ii) )   &
                                                     + dZij(1:Ntz,vvol) * (                - dZodR(1:Ntz,0,ii) ) ) / length(1:Ntz) ) / length(1:Ntz)
              else                ; dLL(1:Ntz) = ( &   ! d/dZbs ;
                                                   + (                - dRodZ(1:Ntz,1,ii) ) * tRij(1:Ntz,vvol) + dZij(1:Ntz,vvol) * im(ii) * cosi(1:Ntz,ii) &
                                                   + ( sini(1:Ntz,ii) - dZodZ(1:Ntz,1,ii) ) * tZij(1:Ntz,vvol) &
                                                   - constraint(1:Ntz) &
                                                   * ( dRij(1:Ntz,vvol) * (                - dRodZ(1:Ntz,1,ii) )   &
                                                     + dZij(1:Ntz,vvol) * ( sini(1:Ntz,ii) - dZodZ(1:Ntz,1,ii) ) ) / length(1:Ntz) ) / length(1:Ntz)
              endif ! end of if( irz.eq.0 ) ;
             else
              if( irz.eq.0 ) then ; dLL(1:Ntz) = ( &   ! d/dRbs ;
                                                   + ( sini(1:Ntz,ii) - dRodR(1:Ntz,1,ii) ) * tRij(1:Ntz,vvol) + dRij(1:Ntz,vvol) * im(ii) * cosi(1:Ntz,ii) &
                                                   + (                - dZodR(1:Ntz,1,ii) ) * tZij(1:Ntz,vvol) &
                                                   - constraint(1:Ntz) &
                                                   * ( dRij(1:Ntz,vvol) * ( sini(1:Ntz,ii) - dRodR(1:Ntz,1,ii) )   &
                                                     + dZij(1:Ntz,vvol) * (                - dZodR(1:Ntz,1,ii) ) ) / length(1:Ntz) ) / length(1:Ntz)
              else                ; dLL(1:Ntz) = ( &   ! d/dZbs ;
                                                   + (                - dRodZ(1:Ntz,0,ii) ) * tRij(1:Ntz,vvol) - dZij(1:Ntz,vvol) * im(ii) * sini(1:Ntz,ii) &
                                                   + ( cosi(1:Ntz,ii) - dZodZ(1:Ntz,0,ii) ) * tZij(1:Ntz,vvol) &
                                                   - constraint(1:Ntz) &
                                                   * ( dRij(1:Ntz,vvol) * (                - dRodZ(1:Ntz,0,ii) )   &
                                                     + dZij(1:Ntz,vvol) * ( cosi(1:Ntz,ii) - dZodZ(1:Ntz,0,ii) ) ) / length(1:Ntz) ) / length(1:Ntz)
              endif ! end of if( irz.eq.0 ) ;
             endif
            else
             if( issym.eq.0 ) then ! take derivatives wrt Rbc and Zbs;
              if( irz.eq.0 ) then ; dLL(1:Ntz) = + ( + cosi(1:Ntz,ii) * tRij(1:Ntz,vvol  ) - dRij(1:Ntz,vvol) * im(ii) * sini(1:Ntz,ii) ) / length(1:Ntz) & 
                                                - constraint(1:Ntz) * dRij(1:Ntz,vvol) * cosi(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
              else                ; dLL(1:Ntz) = + ( + sini(1:Ntz,ii) * tZij(1:Ntz,vvol  ) + dZij(1:Ntz,vvol) * im(ii) * cosi(1:Ntz,ii) ) / length(1:Ntz) & 
                                                 - constraint(1:Ntz) * dZij(1:Ntz,vvol) * sini(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
              endif
             else                  ! take derivatives wrt Rbs and Zbc;
              if( irz.eq.0 ) then ; dLL(1:Ntz) = + ( + sini(1:Ntz,ii) * tRij(1:Ntz,vvol  ) + dRij(1:Ntz,vvol) * im(ii) * cosi(1:Ntz,ii) ) / length(1:Ntz) & 
                                                 - constraint(1:Ntz) * dRij(1:Ntz,vvol) * sini(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
              else                ; dLL(1:Ntz) = + ( + cosi(1:Ntz,ii) * tZij(1:Ntz,vvol  ) - dZij(1:Ntz,vvol) * im(ii) * sini(1:Ntz,ii) ) / length(1:Ntz) & 
                                                 - constraint(1:Ntz) * dZij(1:Ntz,vvol) * cosi(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
              endif ! end of if( irz.eq.0 ) ;
             endif ! end of if( issym.eq.0 ) ;
            endif  ! end of if( Igeometry.eq.3 .and. vvol.eq.1 ) ;
           endif ! end of if( innout.eq.0 ) ;

          endif ! end of if( iocons.eq.0 ) ;

         endif ! end of if( Igeometry.ge.3 ) ;

         call tfft( Nt, Nz, ijreal(1:Ntz), dII(1:Ntz), & ! recall that ijreal contains pressure term;
                    mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), ifail )
         

         call tfft( Nt, Nz, dPP(1:Ntz)   , dLL(1:Ntz), & ! recall that ijreal is probably just a dummy;
                    mn, im(1:mn), in(1:mn), evmn(1:mn), odmn(1:mn), comn(1:mn), simn(1:mn), ifail )          ! evmn and odmn are available as workspace;


         FATAL( dforce, vvol-1+innout.gt.Mvol, psifactor needs attention )

         ; idoc = 0
         
         ;  dFFdRZ(idoc+1:idoc+mn    ,iocons,idof,innout,vvol) = + efmn(1:mn    ) * psifactor(ii,vvol-1+innout) * BBweight(1:mn) ! pressure;
         ; idoc = idoc + mn   ! even;
         ;if( Igeometry.ge.3 ) then ! add spectral constraints;
         ;  dFFdRZ(idoc+1:idoc+mn-1  ,iocons,idof,innout,vvol) = - sfmn(2:mn    ) * psifactor(ii,vvol-1+innout) * epsilon       & ! spectral condensation;
                                                                 - simn(2:mn    ) * psifactor(ii,vvol-1+innout) * sweight(vvol)   ! poloidal length constraint;
         ! if( Ntor.gt.0 ) then
         !  dFFdRZ(idoc+1:idoc+Ntor  ,iocons,idof,innout,vvol) = + odmn(2:Ntor+1) * psifactor(ii,vvol-1+innout) * apsilon
         ! endif
         ; ;idoc = idoc + mn-1 ! oddd;
         ;endif ! end of if( Igeometry.ge.3) ;
         if( NOTstellsym ) then
          ; dFFdRZ(idoc+1:idoc+mn-1  ,iocons,idof,innout,vvol) = + ofmn(2:mn    ) * psifactor(ii,vvol-1+innout) * BBweight(2:mn) ! pressure;
          ;idoc = idoc + mn-1 ! oddd;
          if( Igeometry.ge.3 ) then ! add spectral constraints;
           ;dFFdRZ(idoc+1:idoc+mn    ,iocons,idof,innout,vvol) = - cfmn(1:mn    ) * psifactor(ii,vvol-1+innout) * epsilon       & ! spectral condensation;
                                                                 - comn(1:mn    ) * psifactor(ii,vvol-1+innout) * sweight(vvol)   ! poloidal length constraint;
          !if( Ntor.ge.0 ) then
          ! dFFdRZ(idoc+1:idoc+Ntor+1,iocons,idof,innout,vvol) = + evmn(1:Ntor+1) * psifactor(ii,vvol-1+innout) * apsilon ! poloidal origin      ;
          !endif
           idoc = idoc + mn   ! even;
          endif ! end of if( Igeometry.ge.3) ;
         endif ! end of if( NOTstellsym) ;
         
#ifdef DEBUG
         FATAL( dforce, idoc.ne.LGdof, counting error )
#endif
         
        enddo ! end of do iocons;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
        
       enddo ! matches do innout;
       
      enddo ! matches do issym;

     enddo ! matches do irz;

    enddo ! matches do ii;
    
    dBdX%L = .false. ! probably not needed, but included anyway;
    DALLOCATE( ipivot )
    
   endif ! end of if( LComputeDerivatives ) ;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

2000 continue

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   call intghs_workspace_destroy()

   if (NOTMatrixFree .or. LILUprecond) then
    DALLOCATE(DToocc)
    DALLOCATE(TTssss)
    DALLOCATE(TDstsc)
    DALLOCATE(TDszsc)
    DALLOCATE(DDttcc)
    DALLOCATE(DDtzcc)
    DALLOCATE(DDzzcc)
   endif

   Lsavedguvij = .false.
   DALLOCATE(guvijsave)

   DALLOCATE(Tss)
   DALLOCATE(Dtc)
   DALLOCATE(Dzc)
   DALLOCATE(Ttc)
   DALLOCATE(Tzc)

   if (NOTstellsym) then
    if (NOTMatrixFree .or. LILUprecond) then
      DALLOCATE(DToocs)
      DALLOCATE(DToosc)
      DALLOCATE(DTooss)

      DALLOCATE(TTsscc)
      DALLOCATE(TTsscs)
      DALLOCATE(TTsssc)

      DALLOCATE(TDstcc)
      DALLOCATE(TDstcs)
      DALLOCATE(TDstss)

      DALLOCATE(TDszcc)
      DALLOCATE(TDszcs)
      DALLOCATE(TDszss)

      DALLOCATE(DDttcs)
      DALLOCATE(DDttsc)
      DALLOCATE(DDttss)

      DALLOCATE(DDtzcs)
      DALLOCATE(DDtzsc)
      DALLOCATE(DDtzss)

      DALLOCATE(DDzzcs)
      DALLOCATE(DDzzsc)
      DALLOCATE(DDzzss)
     endif
       
    DALLOCATE(Tsc)
    DALLOCATE(Dts)
    DALLOCATE(Dzs)
    DALLOCATE(Tts)
    DALLOCATE(Tzs)

   end if !NOTstellsym
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   if (NOTMatrixFree) then
     DALLOCATE(dMA)
     DALLOCATE(dMD)
   endif

   DALLOCATE(dMB)
! !DALLOCATE(dMC)
! !DALLOCATE(dME)
! !DALLOCATE(dMF)

   DALLOCATE(dMG)
   
   DALLOCATE(solution)
   
   DALLOCATE(MBpsi)
!  DALLOCATE(MEpsi)
   
   if (LILUprecond) then
     DALLOCATE(dMAS)
     DALLOCATE(dMDS)
     DALLOCATE(idMAS)
     DALLOCATE(jdMAS)
   endif ! if we use GMRES and ILU preconditioner

   if( LcomputeDerivatives) then
    
    DALLOCATE(oBI)
    DALLOCATE(rhs)
    
#ifdef DEBUG
    if( Lcheck.eq.4 ) then
     DALLOCATE(isolution)
    endif
#endif
    
   endif ! end of if( LcomputeDerivatives ) ;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  enddo ! end of do vvol = 1, Mvol (this is the parallelization loop);

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  DALLOCATE(dAt)
  DALLOCATE(dAz)
  
  DALLOCATE( XX) ! spectral constraints; not used;
  DALLOCATE( YY)
  
  DALLOCATE(length)

  DALLOCATE(dRR)
  DALLOCATE(dZZ)

  DALLOCATE(dII)
  DALLOCATE(dLL)
  DALLOCATE(dPP)

  DALLOCATE(constraint)

  DALLOCATE(dBB)

#ifdef DEBUG
  if( Lcheck.eq.3 .or. Lcheck.eq.4 ) then
   DALLOCATE(oRbc)
   DALLOCATE(oZbs)
   DALLOCATE(oRbs)
   DALLOCATE(oZbc)
  endif
#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  if( Lcheck.eq.2 ) then
   write(ounit,'("dforce : ", 10x ," : myid=",i3," ; finished computing derivatives of rotational-transform wrt mu and dpflux ;")') myid
   stop "dforce :            : myid=    ; finished computing derivatives of rotational-transform wrt mu and dpflux ;" ! this will allow other cpus to finish;
  endif
  FATAL( dforce, Lcheck.eq.2, finished computing derivatives of rotational-transform wrt mu and dpflux )
#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef DEBUG
  if( Wdforce ) write(ounit,'("dforce : " 10x " : myid="i3" ; LComputeDerivatives="L2" ; ImagneticOK="999L2)') myid, LComputeDerivatives, ImagneticOK(1:Mvol)
#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  do vvol = 1, Mvol

   LREGION( vvol )

   WCALL( dforce, brcast, ( vvol ) )

  enddo
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef DEBUG
  if( Wdforce ) write(ounit,'("dforce : " 10x " : myid="i3" ; LComputeDerivatives="L2" ; ImagneticOK="999L2)') myid, LComputeDerivatives, ImagneticOK(1:Mvol)
#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  lBBintegral(1:Nvol) = lBBintegral(1:Nvol) * half
  
  Energy = sum( lBBintegral(1:Nvol) ) ! should also compute beta;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! construct force;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  ;   force(0:NGdof) = zero
  
  do vvol = 1, Mvol-1

   LREGION(vvol)
   
   tdoc = (vvol-1) * LGdof 
   
   if( ImagneticOK(vvol) .and. ImagneticOK(vvol+1) ) then ! the magnetic fields in the volumes adjacent to this interface are valid;
    
    ;  idoc = 0           ! degree-of-constraint counter; set;
    
    if( Lextrap.eq.1 .and. vvol.eq.1 ) then ! to be made redundant;
     FATAL( dforce, 2.gt.Mvol, psifactor needs attention )
     ;force(tdoc+idoc+1:tdoc+idoc+mn) = position(1:mn) - ( iRbc(1:mn,2) / psifactor(1:mn,2) )
    else
     ;force(tdoc+idoc+1:tdoc+idoc+mn    ) = ( Bemn(1:mn    ,vvol+1,0) - Bemn(1:mn    ,vvol+0,1) ) * BBweight(1:mn) ! pressure imbalance;
    endif
    
    ;  BBe(vvol) = max( sum( abs( force(tdoc+idoc+1:tdoc+idoc+mn  ) ) ) / (mn  ), logtolerance ) ! screen diagnostics;
    
    ;  idoc = idoc + mn   ! degree-of-constraint counter; increment;
    
    if( Igeometry.ge.3 ) then ! add spectral constraints;
    
     ;force(tdoc+idoc+1:tdoc+idoc+mn-1  ) = (                           Iomn(2:mn    ,vvol+0  ) ) * epsilon         & ! spectral constraints;
                                          + (                         + Somn(2:mn    ,vvol+0,1) ) * sweight(vvol+0) & ! poloidal length constraint;
                                          - ( Somn(2:mn    ,vvol+1,0)                           ) * sweight(vvol+1)
          
!     if( Ntor.gt.0 ) then ! poloidal angle origin is not otherwise constrained ;
!      force(tdoc+idoc+1:tdoc+idoc+Ntor  ) = ( Pomn(2:Ntor+1,vvol+1,0) - Pomn(2:Ntor+1,vvol+0,1) ) * apsilon ! choice of spectral constraint can be enforced;
!     endif
     
     ;IIo(vvol) = max( sum( abs( force(tdoc+idoc+1:tdoc+idoc+mn-1) ) ) / (mn-1), logtolerance ) ! screen diagnostics;
     
     ; idoc = idoc + mn-1
     
    endif ! end of if( Igeometry.ge.3 ) ;
    
    if( NOTstellsym ) then
     
     ;force(tdoc+idoc+1:tdoc+idoc+mn-1  ) = ( Bomn(2:mn    ,vvol+1,0) - Bomn(2:mn    ,vvol+0,1) ) * BBweight(2:mn) ! pressure imbalance;
     
     ; BBo(vvol) = max( sum( abs( force(tdoc+idoc+1:tdoc+idoc+mn-1) ) ) / (mn-1), logtolerance ) ! screen diagnostics;
     
     ; idoc = idoc + mn-1 ! degree-of-constraint counter; increment;
     
     if( Igeometry.ge.3 ) then ! add spectral constraints;
      
      force(tdoc+idoc+1:tdoc+idoc+mn    ) = (                           Iemn(1:mn    ,vvol+0  ) ) * epsilon         & ! spectral constraints;
                                          + (                         + Semn(1:mn    ,vvol+0,1) ) * sweight(vvol+0) & ! poloidal length constraint;
                                          - ( Semn(1:mn    ,vvol+1,0)                           ) * sweight(vvol+1)
      
!     if( Ntor.ge.0 ) then
!      force(tdoc+idoc+1:tdoc+idoc+Ntor+1) = ( Pemn(1:Ntor+1,vvol+1,0) - Pemn(1:Ntor+1,vvol+0,1) ) * apsilon ! choice of spectral constraint can be enforced;
!     endif
      
      ;IIe(vvol) = max( sum( abs( force(tdoc+idoc+1:tdoc+idoc+mn  ) ) ) / (mn  ), logtolerance ) ! screen diagnostics;
      
      ;idoc = idoc + mn   ! degree-of-constraint counter; increment;
      
     endif ! end of if( Igeometry.ge.3 ) ;
     
    endif ! end of if( NOTstellsym ) ;
    
#ifdef DEBUG
    FATAL( dforce, idoc.ne.LGdof, counting error ) ! this has caught bugs;
#endif
    
   else ! matches if( ImagneticOK(vvol) .and. ImagneticOK(vvol+1) );
    
    ;                       ; BBe(vvol) = 9.9E+09
    ;                       ; IIo(vvol) = 9.9E+09
    if ( NOTstellsym ) then ; BBo(vvol) = 9.9E+09
     ;                      ; IIe(vvol) = 9.9E+09
    endif
    
    ; force(tdoc+1:tdoc+LGdof) = 9.9E+09
    
   endif ! end of if( ImagneticOK(vvol) .and. ImagneticOK(vvol+1) ) ;
   
  enddo ! end of do vvol;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( NGdof.ne.0 ) then ; ForceErr = sqrt( sum( force(1:NGdof)*force(1:NGdof) ) / NGdof ) ! this includes spectral constraints;
  else                  ; ForceErr = zero
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  
  if( Wdforce .and. myid.eq.0 ) then
   
   cput = GETTIME
   ;                   ; write(ounit,4000) cput-cpus, ForceErr, cput-cpuo, "|BB|e", alog10(BBe(1:min(Mvol-1,28)))
   if( Igeometry.ge.3 )  write(ounit,4001)                                 "|II|o", alog10(IIo(1:min(Mvol-1,28)))
   if( NOTstellsym ) then
    ;                  ; write(ounit,4001)                                 "|BB|o", alog10(BBo(1:min(Mvol-1,28)))
    if( Igeometry.ge.3 ) write(ounit,4001)                                 "|II|e", alog10(IIe(1:min(Mvol-1,28)))
   endif
   
  endif ! end of if( Wdforce .and. myid.eq.0 ) ;
  
#endif

4000 format("dforce : ",f10.2," : ",6x,3x,"; ",:,"|f|=",es12.5," ; ",:,"time=",f10.2,"s ;",:," log",a5,"=",28f6.2  ," ...")
4001 format("dforce : ", 10x ," : ",6x,3x,"; ",:,"    ",  12x ,"   ",:,"     ", 10x ,"  ;",:," log",a5,"=",28f6.2  ," ...")

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( LcomputeDerivatives ) then ! construct Hessian;
   
#ifdef DEBUG
   FATAL( dforce, .not.Lhessianallocated, need to allocate hessian )
#endif
   
   hessian(1:NGdof,1:NGdof) = zero 
   
   do vvol = 1, Mvol-1 ! loop over interior surfaces;
    
    if( ImagneticOK(vvol) .and. ImagneticOK(vvol+1) ) then ! the magnetic fields in the volumes adjacent to this interface are valid;
     
     idof = 0 ! labels degree-of-freedom = Fourier harmonic of surface geometry;
     
#ifdef DEBUG
     if( idof.gt.LGdof ) write(ounit,1000) myid, vvol, -1, -1, -1, idof, LGdof ! can be deleted;
#endif

     do ii = 1, mn ! loop over degrees-of-freedom;
      
#ifdef DEBUG
      if( idof.gt.LGdof ) write(ounit,1000) myid, vvol, ii, -1, -1, idof, LGdof ! can be deleted;
#endif

      do irz = 0, 1 ! Fourier harmonic of R, Fourier harmonic of Z;
       
#ifdef DEBUG
       if( idof.gt.LGdof ) write(ounit,1000) myid, vvol, ii, irz, -1, idof, LGdof ! can be deleted;
#endif

       if( irz.eq.1 .and. Igeometry.lt.3 ) cycle ! no dependence on Z;
       
       do issym = 0, 1 ! stellarator symmetry;

#ifdef DEBUG
        if( idof.gt.LGdof ) write(ounit,1000) myid, vvol, ii, irz, issym, idof, LGdof ! can be deleted;
#endif
        
        if( issym.eq.1 .and. YESstellsym ) cycle ! no dependence on the non-stellarator symmetric harmonics;
        
        if( ii.eq.1 .and. irz.eq.1 .and. issym.eq.0 ) cycle ! no dependence on Zbs_{m=0,n=0};
        if( ii.eq.1 .and. irz.eq.0 .and. issym.eq.1 ) cycle ! no dependence on Rbs_{m=0,n=0};
        
        idof = idof + 1 ! labels degree-of-freedom;

#ifdef DEBUG
        if( idof.gt.LGdof ) write(ounit,1000) myid, vvol, ii, irz, issym, idof, LGdof ! can be deleted;

1000 format("hforce : " 10x " : myid=",i3," ; vvol=",i3," ; ii= ",i3," ; irz="i3" ; issym="i3" ; idof="i3" ; LGdof="i3" ;")
        
        FATAL( hforce, idof.gt.LGdof, illegal degree-of-freedom index constructing hessian ) ! can be deleted;
#endif
        
        if( vvol.gt.1 ) then
         
         tdof = (vvol-2) * LGdof + idof ! labels degree-of-freedom in internal interface geometry   ;
         tdoc = (vvol-1) * LGdof        ! labels force-balance constraint across internal interfaces;
         idoc = 0                       ! local  force-balance constraint across internal interface ;
         hessian(tdoc+idoc+1:tdoc+idoc+LGdof,tdof) =                                           - dFFdRZ(idoc+1:idoc+LGdof,1,idof,0,vvol+0)
         if( Lconstraint.eq.1 ) then ! this is a little clumsy; could include Lfreebound or something . . . ;
         hessian(tdoc+idoc+1:tdoc+idoc+LGdof,tdof) = hessian(tdoc+idoc+1:tdoc+idoc+LGdof,tdof)                     &
                                                   - dBBdmp(idoc+1:idoc+LGdof,vvol+0,1,1) * dmupfdx(vvol,1,idof,0) &
                                                   - dBBdmp(idoc+1:idoc+LGdof,vvol+0,1,2) * dmupfdx(vvol,2,idof,0)
         endif ! end of if( Lconstraint.eq.1 ) ; 
      
        endif ! end of if( vvol.gt.1 ) ;
        

        ;tdof = (vvol-1) * LGdof + idof
        ;tdoc = (vvol-1) * LGdof ! shorthand;
        ;idoc = 0
        if( Lextrap.eq.1 .and. vvol.eq.1 ) then
        ;hessian(tdoc+idof                  ,tdof) = one ! diagonal elements;
        else
        ;hessian(tdoc+idoc+1:tdoc+idoc+LGdof,tdof) = dFFdRZ(idoc+1:idoc+LGdof,0,idof,0,vvol+1) - dFFdRZ(idoc+1:idoc+LGdof,1,idof,1,vvol+0)
         if( Lconstraint.eq.1 ) then ! this is a little clumsy;
         hessian(tdoc+idoc+1:tdoc+idoc+LGdof,tdof) = hessian(tdoc+idoc+1:tdoc+idoc+LGdof,tdof)                       &
                                                   + dBBdmp(idoc+1:idoc+LGdof,vvol+1,0,1) * dmupfdx(vvol+1,1,idof,0) &
                                                   + dBBdmp(idoc+1:idoc+LGdof,vvol+1,0,2) * dmupfdx(vvol+1,2,idof,0) &
                                                   - dBBdmp(idoc+1:idoc+LGdof,vvol+0,1,1) * dmupfdx(vvol+0,1,idof,1) &
                                                   - dBBdmp(idoc+1:idoc+LGdof,vvol+0,1,2) * dmupfdx(vvol+0,2,idof,1)
         endif ! end of if( Lconstraint.eq.1 );
         endif
        
         if( vvol.lt.Mvol-1 ) then

         tdof = (vvol+0) * LGdof + idof
         tdoc = (vvol-1) * LGdof ! shorthand;
         idoc = 0
         if( Lextrap.eq.1 .and. vvol.eq.1 ) then
         if    ( im(idof).le.0                     ) then ; hessian(tdoc+idof,tdof) = - one
         else                                             ; hessian(tdoc+idof,tdof) = - one
         endif
         else
         hessian(tdoc+idoc+1:tdoc+idoc+LGdof,tdof) = dFFdRZ(idoc+1:idoc+LGdof,0,idof,1,vvol+1)
         if( Lconstraint.eq.1 ) then ! this is a little clumsy;
         hessian(tdoc+idoc+1:tdoc+idoc+LGdof,tdof) = hessian(tdoc+idoc+1:tdoc+idoc+LGdof,tdof)                       &
                                                   + dBBdmp(idoc+1:idoc+LGdof,vvol+1,0,1) * dmupfdx(vvol+1,1,idof,1) &
                                                   + dBBdmp(idoc+1:idoc+LGdof,vvol+1,0,2) * dmupfdx(vvol+1,2,idof,1)
         endif ! end of if( Lconstraint.eq.1 ) then;
         endif

        endif ! end of if( vvol.lt.Mvol-1 ) ;

        if( vvol.eq.Mvol-1 ) then
        !tdof = (vvol+0) * LGdof + idof
         tdoc = (vvol-1) * LGdof ! shorthand ;
         idoc = 0
         dessian(tdoc+idoc+1:tdoc+idoc+LGdof,idof) = dFFdRZ(idoc+1:idoc+LGdof,0,idof,1,vvol+1)
         if( Lconstraint.eq.1 ) then ! this is a little clumsy;
         dessian(tdoc+idoc+1:tdoc+idoc+LGdof,idof) = dessian(tdoc+idoc+1:tdoc+idoc+LGdof,idof)                       &
                                                   + dBBdmp(idoc+1:idoc+LGdof,vvol+1,0,1) * dmupfdx(vvol+1,1,idof,1) &
                                                   + dBBdmp(idoc+1:idoc+LGdof,vvol+1,0,2) * dmupfdx(vvol+1,2,idof,1)
         endif ! end of if( Lconstraint.eq.1 ) then;

         
        endif ! end of if( vvol.lt.Mvol-1 ) ;
        
       enddo ! matches do issym ;
       
      enddo ! matches do irz ;
      
     enddo ! matches do ii ;
     
    else ! matches if( ImagneticOK(vvol) .and. ImagneticOK(vvol+1) ) ; 
     
     FATAL( dforce, .true., need to provide suitable values for hessian in case of field failure )
     
    endif ! end of if( ImagneticOK(vvol) .and. ImagneticOK(vvol+1) ) ;
    
   enddo ! end of do vvol;
   
  endif ! end of if( LcomputeDerivatives ) ;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(dforce)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine dforce
 
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
