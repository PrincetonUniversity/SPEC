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
!latex        \link{coords},
!latex        \link{curent} and 
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

subroutine dforce( NGdof, position, force, LComputeDerivatives )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half, one, two, pi2, pi
  
  use numerical, only : vsmall, small, logtolerance
  
  use fileunits, only : ounit
  
  use inputlist, only : Wmacros, Wdforce, ext, Nvol, Mpol, Ntor, Lrad, tflux, Igeometry, &
                        gamma, adiabatic, pscale, mu, &
                        epsilon, &
                        Lfindzero, &
                        Lconstraint, Lcheck, &
                        Lextrap, &
			mupftol
  
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
                        dMA, dMB,      dMD,           dMG, solution,  &
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
                        dRodR, dRodZ, dZodR, dZodZ, &
                        LocalConstraint, xoffset
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in)  :: NGdof               ! dimensions;
  REAL,    intent(in)  :: position(0:NGdof)   ! degrees-of-freedom = internal geometry;
  REAL,    intent(out) :: force(0:NGdof)      ! force;
  LOGICAL, intent(in)  :: LComputeDerivatives ! indicates whether derivatives are to be calculated;
  
  INTEGER              :: NN, IA, ifail, if01adf, vflag, MM, LDA, idgetrf, idgetri, Lwork
  INTEGER, allocatable :: ipivot(:)
  REAL   , allocatable :: work(:)
  
  INTEGER              :: vvol, innout, ii, jj, irz, issym, iocons, tdoc, idoc, idof, tdof, jdof, ivol, imn, ll, ihybrd1, lwa, Ndofgl
  INTEGER              :: maxfev, ml, muhybr, mode, nprint, nfev, ldfjac, lr
  DOUBLE PRECISION     :: epsfcn, factor
  DOUBLE PRECISION     :: Fdof(1:Mvol-1), Xdof(1:Mvol-1), Fvec(1:Mvol-1)
  DOUBLE PRECISION     :: diag(1:Mvol-1), qtf(1:Mvol-1), wa1(1:Mvol-1), wa2(1:Mvol-1), wa3(1:Mvol-1), wa4(1:mvol-1)
  DOUBLE PRECISION, allocatable :: fjac(:, :), r(:) 

  INTEGER              :: Lcurvature, ideriv, id
  
  REAL                 :: lastcpu, lss, lfactor, DDl, MMl

  INTEGER              :: iflag
  REAL                 :: det
  
  REAL                 :: dpsi(1:2)
  REAL   , allocatable :: oBI(:,:), rhs(:) ! original Beltrami-matrix inverse; used to compute derivatives of matrix equation;

  REAL   , allocatable :: dAt(:,:), dAz(:,:), XX(:), YY(:), dBB(:,:), dII(:), dLL(:), dPP(:), length(:), dRR(:,:), dZZ(:,:), constraint(:)

  CHARACTER            :: packorunpack 

  EXTERNAL 	       :: dfp100

#ifdef DEBUG
  INTEGER              :: isymdiff
  REAL                 :: dRZ = 1.0e-05, dvol(-1:+1), evolume, imupf(1:2,-2:2)
  REAL,    allocatable :: oRbc(:,:), oZbs(:,:), oRbs(:,:), oZbc(:,:) ! original geometry;
  REAL,    allocatable :: isolution(:,:)
#endif

  BEGIN(dforce)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! Unpack position to generate arrays iRbc, iZbs, IRbs, iZbc.

  packorunpack = 'U' ! unpack geometrical degrees-of-freedom;
  
  WCALL( dforce, packxi,( NGdof, position(0:NGdof), Mvol, mn, iRbc(1:mn,0:Mvol), iZbs(1:mn,0:Mvol), iRbs(1:mn,0:Mvol), iZbc(1:mn,0:Mvol), packorunpack ) )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG

! Store initial arrays for debug purposes.  

  if( Lcheck.eq.3 .or. Lcheck.eq.4 ) then ! will check volume derivatives;
   
   SALLOCATE( oRbc, (1:mn,0:Mvol), iRbc(1:mn,0:Mvol) )
   SALLOCATE( oZbs, (1:mn,0:Mvol), iZbs(1:mn,0:Mvol) )
   SALLOCATE( oRbs, (1:mn,0:Mvol), iRbs(1:mn,0:Mvol) )
   SALLOCATE( oZbc, (1:mn,0:Mvol), iZbc(1:mn,0:Mvol) )  
   
  endif ! end of if( Lcheck.eq.3 .or. Lcheck.eq.4 ) ;
  
#endif
  
 
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

   NN = NAdof(vvol) ! shorthand;

   ll = Lrad(vvol)

   SALLOCATE( DToocc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( DToocs, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( DToosc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( DTooss, (0:ll,0:ll,1:mn,1:mn), zero )

   SALLOCATE( TTsscc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( TTsscs, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( TTsssc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( TTssss, (0:ll,0:ll,1:mn,1:mn), zero )

   SALLOCATE( TDstcc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( TDstcs, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( TDstsc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( TDstss, (0:ll,0:ll,1:mn,1:mn), zero )

   SALLOCATE( TDszcc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( TDszcs, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( TDszsc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( TDszss, (0:ll,0:ll,1:mn,1:mn), zero )

   SALLOCATE( DDttcc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( DDttcs, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( DDttsc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( DDttss, (0:ll,0:ll,1:mn,1:mn), zero )

   SALLOCATE( DDtzcc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( DDtzcs, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( DDtzsc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( DDtzss, (0:ll,0:ll,1:mn,1:mn), zero )

   SALLOCATE( DDzzcc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( DDzzcs, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( DDzzsc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( DDzzss, (0:ll,0:ll,1:mn,1:mn), zero )

   WCALL( dforce, ma00aa, ( Iquad(vvol), mn, vvol, ll ) ) ! compute volume integrals of metric elements - evaluate TD, DT, DD, ...;
   WCALL( dforce, matrix, ( vvol, mn, ll ) )

   DALLOCATE(DToocc)
   DALLOCATE(DToocs)
   DALLOCATE(DToosc)
   DALLOCATE(DTooss)

   DALLOCATE(TTsscc)
   DALLOCATE(TTsscs)
   DALLOCATE(TTsssc)
   DALLOCATE(TTssss)

   DALLOCATE(TDstcc)
   DALLOCATE(TDstcs)
   DALLOCATE(TDstsc)
   DALLOCATE(TDstss)

   DALLOCATE(TDszcc)
   DALLOCATE(TDszcs)
   DALLOCATE(TDszsc)
   DALLOCATE(TDszss)

   DALLOCATE(DDttcc)
   DALLOCATE(DDttcs)
   DALLOCATE(DDttsc)
   DALLOCATE(DDttss)

   DALLOCATE(DDtzcc)
   DALLOCATE(DDtzcs)
   DALLOCATE(DDtzsc)
   DALLOCATE(DDtzss)

   DALLOCATE(DDzzcc)
   DALLOCATE(DDzzcs)
   DALLOCATE(DDzzsc)
   DALLOCATE(DDzzss)


  enddo

  if( LocalConstraint ) then

	Ndofgl = 0; Fvec(1:Mvol-1) = 0; iflag = 0;
	Xdof(1:Mvol-1) = dpflux(2:Mvol) + xoffset

	WCALL(dforce, dfp100, (Ndofgl, Xdof, Fvec, iflag) )

	do vvol = 1, Mvol

		WCALL(dforce, dfp200, ( NGdof, position, LcomputeDerivatives, vvol) )
   
	enddo ! end of do vvol = 1, Mvol (this is the parallelization loop);

  else
	
!   If global constraint, start the minimization of the constraint using hybrd1
    Ndofgl = Mvol-1; lwa = 8 * Ndofgl * Ndofgl; maxfev = 1000; nfev=0; lr=Mvol*(Mvol-1); ldfjac=Mvol-1
    ml = Mvol-2; muhybr = Mvol-2; epsfcn=1E-16; diag=0.0; mode=1; factor=0.01; nprint=-1;

    Xdof(1:Mvol-1)   = dpflux(2:Mvol) + xoffset

    SALLOCATE(fjac, (1:ldfjac,1:Mvol-1), 0)
    SALLOCATE(r, (1:lr), 0)

    WCALL( dforce,  hybrd, (dfp100, Ndofgl, Xdof(1:Ndofgl), Fvec(1:Ndofgl), mupftol, maxfev, ml, muhybr, epsfcn, diag(1:Ndofgl), mode, &
			    factor, nprint, ihybrd1, nfev, fjac(1:Ndofgl,1:Ndofgl), ldfjac, r(1:lr), lr, qtf(1:Ndofgl), wa1(1:Ndofgl), &
			    wa2(1:Ndofgl), wa3(1:Ndofgl), wa4(1:Ndofgl)) ) 
 
    dpflux(2:Mvol) = Xdof(1:Ndofgl) - xoffset
 
    do vvol = 1, Mvol
	WCALL(dforce, dfp200, ( NGdof, position, LcomputeDerivatives, vvol) )
    enddo
  
  endif

#ifdef DEBUG
  select case( ihybrd1 )
    case( 1   )  ; write(ounit,'("dforce : ",f10.2," : finished ; success        ; dpflux = ", es12.5, ", its="i7";")') cput-cpus, dpflux, nfev
    case( 0   )  ; write(ounit,'("dforce : ",f10.2," : finished ; input error    ; dpflux = ", es12.5, ", its="i7";")') cput-cpus, dpflux, nfev
    case( 2   )  ; write(ounit,'("dforce : ",f10.2," : finished ; max. iter      ; dpflux = ", es12.5, ", its="i7";")') cput-cpus, dpflux, nfev
    case( 3   )  ; write(ounit,'("dforce : ",f10.2," : finished ; xtol too small ; dpflux = ", es12.5, ", its="i7";")') cput-cpus, dpflux, nfev
    case( 4:5 )  ; write(ounit,'("dforce : ",f10.2," : finished ; bad progress   ; dpflux = ", es12.5, ", its="i7";")') cput-cpus, dpflux, nfev
    case default ; write(ounit,'("dforce : ",f10.2," : finished ; illegal ifail  ; dpflux = ", es12.5, ", its="i7";")') cput-cpus, dpflux, nfev
  end select
#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
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
         hessian(tdoc+idoc+1:tdoc+idoc+LGdof,tdof) =                                           - dFFdRZ(idoc+1:idoc+LGdof,vvol+0,1,idof,0)
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
        ;hessian(tdoc+idoc+1:tdoc+idoc+LGdof,tdof) = dFFdRZ(idoc+1:idoc+LGdof,vvol+1,0,idof,0) - dFFdRZ(idoc+1:idoc+LGdof,vvol+0,1,idof,1)
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
         hessian(tdoc+idoc+1:tdoc+idoc+LGdof,tdof) = dFFdRZ(idoc+1:idoc+LGdof,vvol+1,0,idof,1)
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
         dessian(tdoc+idoc+1:tdoc+idoc+LGdof,idof) = dFFdRZ(idoc+1:idoc+LGdof,vvol+1,0,idof,1)
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
