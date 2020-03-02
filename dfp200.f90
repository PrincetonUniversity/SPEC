!title (&ldquo;global&rdquo; dfp200) ! Given the field consistent with the constraints and the geometry, computes local quantites related to the force evaluation.

!latex \briefly{Calculates ${ F_i}({\bf x})$, where ${\bf x} \equiv \{\mbox{\rm geometry}\} \equiv \{ R_{i,v}, Z_{i,v}\}$ 
!latex          and ${ F_i}\equiv p_i+B_i^2/2 + \{\mbox{\rm spectral constraints}\} $, and $\nabla {\bf F_i}$.}

!latex \calledby{\link{dforce}} \\

!latex \calls{\link{ma00aa}, 
!latex        \link{ma02aa}, 
!latex        \link{coords}, 
!latex        \link{dlasrt}, 
!latex        \link{lforce}, 
!latex        \link{volume}, 
!latex        \link{packab}, 
!latex        \link{tr00ab}, 
!latex        \link{curent}, 
!latex        \link{matrix},}


!latex \tableofcontents


!latex \subsection{Construction of local force}
!latex See \link{dforce} documentation for more details

!latex \subsection{Construction of matrix equation derivatives}
!latex \begin{enumerate}
!latex \item Matrix perturbation theory is used to compute the derivatives of the solution, i.e. the Beltrami fields, as the geometry of the 
!latex       interfaces changes:
!latex \end{enumerate}

!latex \subsection{Extrapolation: planned redundant}
!latex \begin{enumerate}
!latex \item The extrapolation constraint is $R_{j,1} = R_{j,2} \, \psi_1^{m/2} / \psi_2^{m/2}$.
!latex       Combining this with the regularization factor for the geometry, i.e. $R_{j,i}=\psi_{i}^{m/2} \xi_{j,i}$, we obtain
!latex       \be \xi_{j,1} = R_{j,2} / \psi_2^{m/2}.
!latex       \ee
!latex \end{enumerate}

subroutine dfp200( LcomputeDerivatives, vvol)

  use constants, only : zero, half, one, two
  
  use numerical, only : small
  
  use fileunits, only : ounit
  
  use inputlist, only : Wmacros, Wdfp200, ext, Nvol, Mpol, Ntor, Lrad, tflux, Igeometry, &
                        gamma, adiabatic, pscale, mu, &
                        epsilon, &
                        Lfindzero, &
                        Lconstraint, Lcheck, &
                        Lextrap
  
  use cputiming, only : Tdfp200
  
  use allglobal, only : ncpu, myid, cpus, &
                        Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, &
                        Mvol, &
                        Iquad, &                  ! convenience; provided to ma00aa as argument to avoid allocations;
                        iRbc, iZbs, iRbs, iZbc, & ! Fourier harmonics of geometry; vector of independent variables, position, is "unpacked" into iRbc,iZbs;
                        NAdof, &
                        YESstellsym, NOTstellsym, &
                        mn, im, in, mns, Ntz, &
                        Ate, Aze, Ato, Azo, & ! only required for debugging;
                        ijreal, &
                        efmn, ofmn, cfmn, sfmn, &
                        evmn, odmn, comn, simn, &
                        Nt, Nz, &
                        cosi, sini, & ! FFT workspace;
                        dBdX, &
                        dMA, dMB, dMD, dMG, solution, &
                        dtflux, dpflux, sweight, &
                        mmpp, &
                        Bemn, Bomn, Iomn, Iemn, Somn, Semn, &
                        LGdof, &
                        vvolume, dvolume, &
                        Rij, Zij, sg, guvij, iRij, iZij, dRij, dZij, tRij, tZij, & ! Jacobian and metrics; computed in coords;
                        diotadxup, dItGpdxtp, &
                        dFFdRZ, dBBdmp, dmupfdx, hessian, dessian, Lhessianallocated, &
                        BBweight, & ! exponential weight on force-imbalance harmonics;
                        psifactor, &
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
                        LocalConstraint, &
						IsMyVolume, IsMyVolumeValue, IndMatrixArray, Btemn
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  LOGICAL, intent(in)  :: LComputeDerivatives ! indicates whether derivatives are to be calculated;
  
  INTEGER              :: NN, IA, ifail, if01adf, vflag, MM, idgetrf, idgetri, Lwork, lvol
  INTEGER              :: vvol, innout, ii, jj, irz, issym, iocons, idoc, idof, imn, ll
  INTEGER              :: Lcurvature, ideriv, id, ind_matrix
  INTEGER              :: iflag
  INTEGER, allocatable :: ipivot(:)

  REAL                 :: lastcpu, lss, lfactor, DDl, MMl
  REAL                 :: det
  REAL   , allocatable :: oBI(:,:), rhs(:) ! original Beltrami-matrix inverse; used to compute derivatives of matrix equation;
  REAL   , allocatable :: dAt(:,:), dAz(:,:), XX(:), YY(:), dBB(:,:), dII(:), dLL(:), dPP(:), length(:), dRR(:,:), dZZ(:,:), constraint(:)
  REAL   , allocatable :: work(:)

  CHARACTER            :: packorunpack 




BEGIN(dfp200)


#ifdef DEBUG
	if( Lcheck.eq.2 ) then
		goto 2000 ! will take no other action except a finite-difference comparison on the derivatives of the rotational-transform wrt mu and dpflux; (in dforce)
	endif
#endif

if( LocalConstraint ) then
	do vvol = 1, Mvol

		WCALL(dfp200, IsMyVolume, (vvol))

		if( IsMyVolumeValue .EQ. 0 ) then
			cycle
		else if( IsMyVolumeValue .EQ. -1) then
			FATAL(dfp200, .true., Unassociated volume)
		endif
					
			
		LREGION(vvol) ! assigns Lcoordinatesingularity, Lplasmaregion, etc. ;
		
		dBdX%vol = vvol  ! Useless when local constraint
		ll = Lrad(vvol)  ! Shorthand
	   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

	   	NN = NAdof(vvol) ! shorthand;
		SALLOCATE( dAt       , (1:Ntz,-1:2), zero )
		SALLOCATE( dAz       , (1:Ntz,-1:2), zero )
		SALLOCATE(  XX       , (1:Ntz     ), zero )
		SALLOCATE(  YY       , (1:Ntz     ), zero )
		SALLOCATE( length    , (1:Ntz     ), zero ) ! this is calculated in lforce;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

		vflag = 1
		WCALL( dfp200, volume, ( vvol, vflag ) ) ! compute volume;


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
		do iocons = 0, 1 ! construct field magnitude on inner and outer interfaces; inside do vvol;
		
			if( vvol.eq.1    .and. iocons.eq.0 ) cycle ! fixed inner boundary (or coordinate axis);
			if( vvol.eq.Mvol .and. iocons.eq.1 ) cycle ! fixed outer boundary                     ; there are no constraints at outer boundary;
			
			ideriv = 0 ; id = ideriv
			iflag = 0 ! dAt, dAz, XX & YY are returned by lforce; Bemn(1:mn,vvol,iocons), Iomn(1:mn,vvol) etc. are returned through global;
			WCALL( dfp200, lforce, ( vvol, iocons, ideriv, Ntz, dAt(1:Ntz,id), dAz(1:Ntz,id), XX(1:Ntz), YY(1:Ntz), length(1:Ntz), DDl, MMl, iflag ) )
		   
		enddo ! end of do iocons = 0, 1;


	   	if( LcomputeDerivatives ) then ! compute inverse of Beltrami matrices;


			SALLOCATE( dBB       , (1:Ntz,-1:2), zero ) ! magnetic field strength (on interfaces) in real space and derivatives;
			SALLOCATE( dRR       , (1:Ntz,-1:1), zero )
			SALLOCATE( dZZ       , (1:Ntz,-1:1), zero )
			SALLOCATE( dII       , (1:Ntz     ), zero ) ! spectral constraint;
			SALLOCATE( dLL       , (1:Ntz     ), zero ) ! length   constraint;
			SALLOCATE( dPP       , (1:Ntz     ), zero ) ! poloidal constraint;
			SALLOCATE( constraint, (1:Ntz     ), zero )
			SALLOCATE( oBI		 , (1:NN,1:NN ), zero ) ! inverse of ``original'', i.e. unperturbed, Beltrami matrix;

			SALLOCATE( rhs, (1:NN     ), zero )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

	   		call get_inverse_Beltrami_matrices(vvol, oBI, NN)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
			dBdX%L = .true. ! will need derivatives;
			idof = 0 ! labels degree of freedom; local to interface;
			   
				
			do ii = 1, mn ! loop over deformations in Fourier harmonics; inside do vvol;
				
				dBdX%ii = ii ! controls construction of derivatives in subroutines called below;     
				do irz = 0, 1 ! loop over deformations in R and Z; inside do vvol; inside do ii;
					
					if( irz.eq.1 .and. Igeometry.lt.3 ) cycle ! no dependence on Z;
					dBdX%irz = irz ! controls construction of derivatives;
					
					do issym = 0, 1 ! loop over stellarator and non-stellarator symmetric terms;

						if( issym.eq.1 .and. YESstellsym               ) cycle ! no dependence on non-stellarator symmetric harmonics;
						if( ii.eq.1    .and. irz.eq.1 .and. issym.eq.0 ) cycle ! no dependence on Zbs_{m=0,n=0};
						if( ii.eq.1    .and. irz.eq.0 .and. issym.eq.1 ) cycle ! no dependence on Rbs_{m=0,n=0};

						dBdX%issym = issym ! controls construction of derivatives;
						idof = idof + 1 ! this labels the degree-of-freedom that the derivative is taken with respect to; this is outside do innout;

#ifdef DEBUG
						FATAL( dfp200, idof.gt.LGdof, illegal degree-of-freedom index constructing derivatives ) ! this can be deleted;
#endif

						do innout = 0, 1 ! loop over deformations to inner and outer interface; inside do vvol; inside do ii; inside do irz;

							if( vvol.eq.1    .and. innout.eq.0 ) cycle ! no degrees of freedom at coordinate axis / fixed inner boundary;
							if( vvol.eq.Mvol .and. innout.eq.1 ) cycle ! no degress of freedom                      fixed outer boundary; for linearized displacement;

							dBdX%innout = innout

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

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
							! Perturbed solution

							call get_perturbed_solution(vvol, rhs, oBI, NN)

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

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
							! Helicity multiplier and poloidal flux derivatives

							call evaluate_dmupfdx(vvol, innout, idof, ii, issym, irz)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
							! EVALUATE dBB

							call evaluate_dBB(vvol, idof, innout, issym, irz, ii, dAt, dAz, dBB, XX, YY, length, dRR, dZZ, dII, dLL, dPP, Ntz)



						enddo ! matches do innout;
					enddo ! matches do issym;
				enddo ! matches do irz;
			enddo ! matches do ii;
				
			
			dBdX%L = .false. ! probably not needed, but included anyway;



			DALLOCATE(rhs)

			DALLOCATE(oBI)
			DALLOCATE(dRR)
			DALLOCATE(dZZ)
			DALLOCATE(dII)
			DALLOCATE(dLL)
			DALLOCATE(dPP)
			DALLOCATE(constraint)
			DALLOCATE(dBB)
	   
		endif ! end of if( LComputeDerivatives ) ;

		DALLOCATE(dAt)
		DALLOCATE(dAz)
		DALLOCATE( XX) ! spectral constraints; not used;
		DALLOCATE( YY)
		DALLOCATE(length)

	enddo ! matches do vvol = 1, Mvol


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

else ! CASE SEMI GLOBAL CONSTRAINT


	SALLOCATE( dAt       , (1:Ntz,-1:2), zero )
	SALLOCATE( dAz       , (1:Ntz,-1:2), zero )
	SALLOCATE(  XX       , (1:Ntz     ), zero )
	SALLOCATE(  YY       , (1:Ntz     ), zero )
	SALLOCATE( length    , (1:Ntz     ), zero ) ! this is calculated in lforce;

	do vvol = 1, Mvol
			WCALL(dfp200, IsMyVolume, (vvol))

			if( IsMyVolumeValue .EQ. 0 ) then
				cycle
			else if( IsMyVolumeValue .EQ. -1) then
				FATAL(dfp200, .true., Unassociated volume)
			endif
						
				
			LREGION(vvol) ! assigns Lcoordinatesingularity, Lplasmaregion, etc. ;
			ll = Lrad(vvol)  ! Shorthand
		   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

		   	NN = NAdof(vvol) ! shorthand;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

			vflag = 1
			WCALL( dfp200, volume, ( vvol, vflag ) ) ! compute volume;


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
			do iocons = 0, 1 ! construct field magnitude on inner and outer interfaces; inside do vvol;
			
				if( vvol.eq.1    .and. iocons.eq.0 ) cycle ! fixed inner boundary (or coordinate axis);
				if( vvol.eq.Mvol .and. iocons.eq.1 ) cycle ! fixed outer boundary                     ; there are no constraints at outer boundary;
				
				ideriv = 0 ; id = ideriv
				iflag = 0 ! dAt, dAz, XX & YY are returned by lforce; Bemn(1:mn,vvol,iocons), Iomn(1:mn,vvol) etc. are returned through global;
				WCALL( dfp200, lforce, ( vvol, iocons, ideriv, Ntz, dAt(1:Ntz,id), dAz(1:Ntz,id), XX(1:Ntz), YY(1:Ntz), length(1:Ntz), DDl, MMl, iflag ) )
			   
			enddo ! end of do iocons = 0, 1;
	enddo



	if( LComputeDerivatives ) then
		dBdX%L = .true. ! will need derivatives;
						

! Evaluation of froce gradient in case of semi-global constraint (for now, this is specific to Lconstraint=3)
! Loop over all geometrical degrees of freedom. Thi means we loop over
! - each interface (vvol)
! - each Fourier mode (ii)
! - on R or Z (irz)
! - on symmetric and non stellarator symmetric terms (issym)
! 
! Then the perturbed solution is evaluated on volumes sharing the interface vvol (volume vvol and vvol+1).
! i.e. relevant routines are run with lvol=vvol, innout=1 and lvol=vvol+1, innout=0

		! Allocate memory
		SALLOCATE( dBB       , (1:Ntz,-1:2), zero ) ! magnetic field strength (on interfaces) in real space and derivatives;
		SALLOCATE( dRR       , (1:Ntz,-1:1), zero )
		SALLOCATE( dZZ       , (1:Ntz,-1:1), zero )
		SALLOCATE( dII       , (1:Ntz     ), zero ) ! spectral constraint;
		SALLOCATE( dLL       , (1:Ntz     ), zero ) ! length   constraint;
		SALLOCATE( dPP       , (1:Ntz     ), zero ) ! poloidal constraint;
		SALLOCATE( constraint, (1:Ntz     ), zero )

		do vvol = 1, Mvol-1 !labels which interface is perturbed
		dBdX%vol = vvol
		idof = 0 ! labels degree of freedom of interface vvol

		do ii = 1, mn ! loop over deformations in Fourier harmonics; inside do vvol;
			dBdX%ii = ii ! controls construction of derivatives in subroutines called below;     

			do irz = 0, 1 ! loop over deformations in R and Z; inside do vvol; inside do ii;
				if( irz.eq.1 .and. Igeometry.lt.3 ) cycle ! no dependence on Z;
				dBdX%irz = irz ! controls construction of derivatives;	

				do issym = 0, 1 ! loop over stellarator and non-stellarator symmetric terms;
					if( issym.eq.1 .and. YESstellsym               ) cycle ! no dependence on non-stellarator symmetric harmonics;
					if( ii.eq.1    .and. irz.eq.1 .and. issym.eq.0 ) cycle ! no dependence on Zbs_{m=0,n=0};
					if( ii.eq.1    .and. irz.eq.0 .and. issym.eq.1 ) cycle ! no dependence on Rbs_{m=0,n=0};
					dBdX%issym = issym ! controls construction of derivatives;

					idof = idof + 1 ! this labels the degree-of-freedom that the derivative is taken with respect to; this is outside do innout;
#ifdef DEBUG
					FATAL( dfp200, idof.gt.LGdof, illegal degree-of-freedom index constructing derivatives ) ! this can be deleted;
#endif

					do lvol = vvol, vvol+1
					  do innout = 0, 1 ! inside do vvol; inside do ii; inside do irz; 
							if( lvol.eq.1    .and. innout.eq.0 ) cycle ! no degrees of freedom at coordinate axis / fixed inner boundary;
							if( lvol.eq.Mvol .and. innout.eq.1 ) cycle ! no degress of freedom                      fixed outer boundary; for linearized displacement;

                            if( lvol.eq.vvol   .and. innout.eq.0 ) cycle ! do not perturb inner interface of volume vvol
                            if( lvol.eq.vvol+1 .and. innout.eq.1 ) cycle ! do not perturn outer interface of volume vvol+1


							dBdX%innout = innout

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
! Set up volume information
							WCALL(dfp200, IsMyVolume, (lvol))

							if( IsMyVolumeValue .EQ. 0 ) then
								cycle
							else if( IsMyVolumeValue .EQ. -1) then
								FATAL(dfp200, .true., Unassociated volume)
							endif
										
								
							LREGION(lvol) ! assigns Lcoordinatesingularity, Lplasmaregion, etc. ;
							ll = Lrad(lvol)  ! Shorthand
						   	NN = NAdof(lvol) ! shorthand;
						   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

							vflag = 1
							WCALL( dfp200, volume, ( lvol, vflag ) ) ! compute volume;				! TODO: IS THIS NECESSARY?

							! Allocate memory
							SALLOCATE( oBI,    (1:NN,1:NN          ), zero ) ! inverse of ``original'', i.e. unperturbed, Beltrami matrix;
							SALLOCATE( rhs,    (1:NN               ), zero )

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

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

					   		call get_inverse_Beltrami_matrices(lvol, oBI, NN)! Can't be moved outside of loops - need OBI after.

							! Perturbed solution
							call get_perturbed_solution(lvol, rhs, oBI, NN)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

							! Free memory
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

							DALLOCATE(rhs)
							DALLOCATE(oBI)

						enddo ! end of do innout=0,1
					enddo ! end of do lvol = vvol, vvol+1



! At this point, we have the inverse Beltrami matrices dMA, ..., the origin Beltrami matrix oBI 
! and the perturbed solution for each volume neighboring the perturbed interface.
! We now loop again on everything to compute the derivatives of mu and psip w.r.t the position and dBB.

						! Helicity multiplier and poloidal flux derivatives
						call evaluate_dmupfdx(vvol, 1, idof, ii, issym, irz)


						do lvol = vvol, vvol+1
							if( lvol.eq.vvol ) then
                                innout      = 1
								dBdx%innout = innout
							else
                                innout      = 0
								dBdX%innout = innout
							endif

							WCALL(dfp200, IsMyVolume, (lvol))

							if( IsMyVolumeValue .EQ. 0 ) then
								cycle
							else if( IsMyVolumeValue .EQ. -1) then
								FATAL(dfp200, .true., Unassociated volume)
							endif
										
								
							LREGION(lvol) ! assigns Lcoordinatesingularity, Lplasmaregion, etc. ;
						   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


							! EVALUATE dBB
							call evaluate_dBB(lvol, idof, innout, issym, irz, ii, dAt, dAz, dBB, XX, YY, length, dRR, dZZ, dII, dLL, dPP, Ntz)
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
									

					enddo ! matches do lvol = vvol, vvol+1
				enddo ! matches do issym;
			enddo ! matches do irz;
		enddo ! matches do ii;
    enddo ! matches do vvol;
								
	! Free memory
	DALLOCATE(dRR)
	DALLOCATE(dZZ)
	DALLOCATE(dII)
	DALLOCATE(dLL)
	DALLOCATE(dPP)
	DALLOCATE(constraint)
	DALLOCATE(dBB)

	dBdX%L = .false. ! probably not needed, but included anyway;

	endif ! End of if( LComputeDerivatives ) 

DALLOCATE(dAt)
DALLOCATE(dAz)
DALLOCATE( XX) ! spectral constraints; not used;
DALLOCATE( YY)
DALLOCATE(length)


endif ! End of if( LocalConstraint )



2000 continue

RETURN(dfp200)

end subroutine










!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
!																LOCAL SUBROUTINES
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine get_inverse_beltrami_matrices(vvol, oBI, NN)

! Evaluate the inverse of Beltrami matrices and store the original one in oBI.

! INPUT
! -----
! 	vvol: 	Volume number
!   NN:		is equal to NAdof(vvol), is a shorthand.


! MODULES
! -------

	use constants, only :	zero, half, one, two

	use fileunits, only : 	ounit

	use cputiming, only :   Tdfp200

	use inputlist, only :	Wmacros, Wdfp200, Lrad, mu

	use allglobal, only :   ncpu, myid, cpus, &
		                    Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, &
		                    Nt, Nz, &
		                    dBdX, &
		                    dMA, dMB, dMD, dMG, solution, &
		                    mn, mne, IndMatrixArray


  LOCALS
! ------

INTEGER 			:: lastcpu, NN, ind_matrix, vvol, IA, MM, LDA, Lwork
INTEGER 			:: idgetrf, idgetri
REAL, allocatable   :: ipivot(:), work(:)
REAL	 			:: oBI(1:NN, 1:NN)


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

		lastcpu = GETTIME

		ind_matrix = IndMatrixArray(vvol, 2)

		dMA(ind_matrix)%mat(0:NN-1,1:NN) = dMA(ind_matrix)%mat(1:NN,1:NN) - mu(vvol) * dMD(ind_matrix)%mat(1:NN,1:NN) ! this corrupts dMA, but dMA is no longer used;
		dMA(ind_matrix)%mat(  NN  ,1:NN) = zero
		dMD(ind_matrix)%mat(1:NN  ,1:NN) = dMA(ind_matrix)%mat(0:NN-1,1:NN) ! copy of original matrix; this is used below;

    
        IA = NN + 1
        MM = NN ; LDA = NN ; Lwork = NN
        SALLOCATE( ipivot, (1:NN), 0 )
    
        idgetrf = 1 ; call DGETRF( MM, NN, dMA(ind_matrix)%mat(0:LDA-1,1:NN), LDA, ipivot(1:NN), idgetrf ) !LU decomposition

		cput = GETTIME
		select case( idgetrf ) !                                                                     0123456789012345678
			case(  :-1 ) ;               write(ounit,1010) cput-cpus, myid, vvol, cput-lastcpu, idgetrf, "input error;      "
			case(  0   ) ; if( Wdfp200 ) write(ounit,1010) cput-cpus, myid, vvol, cput-lastcpu, idgetrf, "success;          "
			case( 1:   ) ;               write(ounit,1010) cput-cpus, myid, vvol, cput-lastcpu, idgetrf, "singular;         "
			case default ;               FATAL( dfp200, .true., illegal ifail returned from F07ADF )
		end select

1010 format("dfp200 : ",f10.2," : myid=",i3," ; vvol=",i3," ; called DGETRF ; time=",f10.2,"s ; LU factorization of matrix; idgetrf=",i2," ; ",a18)

		SALLOCATE( work, (1:NN), zero )

		! inverse of MA using the LU decomposition of DGETRF
		idgetri = 1 ; call DGETRI( NN, dMA(ind_matrix)%mat(0:LDA-1,1:NN), LDA, ipivot(1:NN), work(1:Lwork), Lwork, idgetri ) 

		DALLOCATE(work)
		DALLOCATE(ipivot)

		cput = GETTIME
		select case( idgetri ) !                                                                     0123456789012345678
			case(  :-1 ) ;               write(ounit,1011) cput-cpus, myid, vvol, cput-lastcpu, idgetri, "input error;      "
			case(  0   ) ; if( Wdfp200 ) write(ounit,1011) cput-cpus, myid, vvol, cput-lastcpu, idgetri, "success;          "
			case( 1:   ) ;               write(ounit,1011) cput-cpus, myid, vvol, cput-lastcpu, idgetri, "singular;         "
			case default ;               FATAL( dfp200, .true., illegal ifail returned from F07AJF )
		end select

1011 format("dfp200 : ",f10.2," : myid=",i3," ; vvol=",i3," ; called DGETRI ; time=",f10.2,"s ; inverse of Beltrami matrix; idgetrf=",i2," ; ",a18)

		oBI(1:NN,1:NN) = dMA(ind_matrix)%mat(0:LDA-1,1:NN)

end subroutine get_inverse_beltrami_matrices






!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-


subroutine get_perturbed_solution(vvol, rhs, oBI, NN)

! This routine evaluates the value of the magnetic field once the interface is perturbed using matrix perturbation theory.
! Separated from the main dfp200 core to allow local and semi-global constraints to be prescribed
!
! Attention: Solution is perturbed for a given degree of freedom (information stored in dBdX)

	use constants, only :	zero, half, one, two

	use fileunits, only : 	ounit

	use cputiming, only :   Tdfp200

	use inputlist, only :	Wmacros, Wdfp200, Lrad, mu

	use allglobal, only :   ncpu, myid, cpus, &
							mn, Iquad, NAdof, &
	                    	dMA, dMB, dMD, dMG, solution, &
							DToocc, DToocs, DToosc, DTooss, &
		                    TTsscc, TTsscs, TTsssc, TTssss, &
	                    	TDstcc, TDstcs, TDstsc, TDstss, &
		                    TDszcc, TDszcs, TDszsc, TDszss, &
		                    DDttcc, DDttcs, DDttsc, DDttss, &
		                    DDtzcc, DDtzcs, DDtzsc, DDtzss, &
		                    DDzzcc, DDzzcs, DDzzsc, DDzzss, &
		                    dtflux, dpflux, IndMatrixArray


 LOCALS
!------

INTEGER					:: ideriv, vvol, ind_matrix, ll, NN
REAL					:: dpsi(1:2)
REAL					:: rhs(1:NN)
REAL 					:: oBI(1:NN,1:NN)
CHARACTER           	:: packorunpack 

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
	
	ll = Lrad(vvol)  ! Shorthand
	ind_matrix = IndMatrixArray(vvol, 2)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
	! Evaluate perturbed geometry dependent matrices

	WCALL( dfp200, ma00aa,( Iquad(vvol), mn, vvol, ll ) ) ! compute volume integrals of metric elements;
	WCALL( dfp200, matrix,( vvol, mn, ll ) ) ! construct Beltrami matrices;

	dpsi(1:2) = (/ dtflux(vvol), dpflux(vvol) /) ! local enclosed toroidal and poloidal fluxes;
	rhs(1:NN) = - matmul( dMB(ind_matrix)%mat(1:NN,1:2 )                            , dpsi(1:2)        ) &
			    - matmul( dMA(ind_matrix)%mat(1:NN,1:NN) - mu(vvol) * dMD(ind_matrix)%mat(1:NN,1:NN), solution(vvol)%mat(1:NN,0) )

	! Evaluate perturbed solution
	solution(vvol)%mat(1:NN,-1) = matmul( oBI(1:NN,1:NN), rhs(1:NN) ) ! this is the perturbed, packxi solution;

	ideriv = -1 ;

	packorunpack = 'U'

	WCALL( dfp200, packab,( packorunpack, vvol, NN,  solution(vvol)%mat(1:NN,-1), ideriv ) ) ! derivatives placed in Ate(vvol,ideriv,1:mn)%s(0:Lrad),

end subroutine get_perturbed_solution





!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine evaluate_dmupfdx(vvol, innout, idof, ii, issym, irz)

! Evaluate mu and psip derivatives and store them in dmupfdx.
! If debug and Lcheck=4, compare with finite difference approximation
!
! INPUT:
! ------
! 	vvol: 	Perturbed interface number
! 	innout: inner or outer boundary. Not used in case of Lconstraint=3 since inside loop
! 	idof:  	???
! 	ii: 	index of Fourier harmonics, ii=1:mn
! 	issym: 	loop over stellerator symmetric and non-symmetric terms
! 	irz: 	loop over R or Z

! MODULES:
! -------
	use constants, only :   zero, half, one, two

	use fileunits, only : 	ounit

	use numerical, only : 	small

	use cputiming, only :   Tdfp200

	use inputlist, only :	Wmacros, Wdfp200, Lrad, mu, Lconstraint, Lcheck, Nvol, Igeometry, mupftol, Lfreebound

  	use allglobal, only : 	ncpu, myid, cpus, &
                        	Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, &
	                        Mvol, Iquad, NGdof, &
                           	iRbc, iZbs, iRbs, iZbc, & ! Fourier harmonics of geometry; vector of independent variables, position, is "unpacked" into iRbc,iZbs;
	                    	NAdof, &
    	                   	mn, im, in, mns, &
    	                   	Ate, Aze, Ato, Azo, & 													! only required for debugging;
    	                   	Nt, Nz, &
    	                   	dBdX, &
    	                   	solution, &
    	                   	dtflux, dpflux, sweight, &
    	                   	Rij, Zij, & 															
    	                   	diotadxup, dItGpdxtp, dmupfdx, &										
    	                   	psifactor, &
    	                   	lmns, &
    	                   	mn, mne, &
    	                   	LocalConstraint, &
							vvolume, dvolume, &
							IsMyVolume, IsMyVolumeValue, IndMatrixArray, &
							DToocc, DToocs, DToosc, DTooss, &
							TTsscc, TTsscs, TTsssc, TTssss, &
					       	TDstcc, TDstcs, TDstsc, TDstss, &
							TDszcc, TDszcs, TDszsc, TDszss, &
							DDttcc, DDttcs, DDttsc, DDttss, &
							DDtzcc, DDtzcs, DDtzsc, DDtzss, &
							DDzzcc, DDzzcs, DDzzsc, DDzzss, &
							Btemn, xoffset


  LOCALS:
! -------

	INTEGER				:: 	vvol, innout, idof, iflag, ii, issym, irz, ll, NN, ifail, vflag, N, iwork(1:Nvol-1), idgesvx, pvol
	REAL				::  det, lfactor, dBdmpf(1:Mvol-1,1:Mvol-1), dBdx2(1:Mvol-1), AF(1:Nvol-1,1:Nvol-1), IPIV(1:Nvol-1)
	REAL 				::  R(1:Nvol-1), C(1:Nvol-1), work(1:4*Nvol-4), ferr, berr, rcond, tmp(2:Nvol)
	LOGICAL 			::  Lonlysolution, LcomputeDerivatives


#ifdef DEBUG
	INTEGER             :: isymdiff, maxfev, nfev, lr, ldfjac, ml, muhybr, epsfcn, mode, nprint
	INTEGER             :: jj, tdoc, idoc, tdof, jdof, imn, ihybrd1, lwa, Ndofgl, llmodnp
	REAL                :: dRZ = 1.0e-5, dvol(-1:+1), evolume, imupf_global(1:Mvol,1:2,-2:2), imupf_local(1:2,-2:2), factor, Btemn_debug(1:mn, 0:1, 1:Mvol, -1:2)
    DOUBLE PRECISION     :: diag(1:Mvol-1), qtf(1:Mvol-1), wa1(1:Mvol-1), wa2(1:Mvol-1), wa3(1:Mvol-1), wa4(1:mvol-1)
	REAL,   allocatable :: oRbc(:,:), oZbs(:,:), oRbs(:,:), oZbc(:,:) ! original geometry;
	REAL,   allocatable :: isolution(:,:)
	REAL                :: position(0:NGdof), force(0:NGdof)
	CHARACTER			:: packorunpack
    DOUBLE PRECISION    ::  Fdof(1:Mvol-1), Xdof(1:Mvol-1), Fvec(1:Mvol-1)
    DOUBLE PRECISION, allocatable :: fjac(:, :), r_deb(:)

	EXTERNAL 			:: dfp100
#endif


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

	ll = Lrad(vvol)		! shorthand
	NN = NAdof(vvol) 	! shorthand;

#ifdef DEBUG
	FATAL( dfp200, vvol-1+innout.gt.Mvol, psifactor needs attention )
#endif


#ifdef DEBUG

! Store initial arrays for debug purposes.  

	if( Lcheck.eq.3 .or. Lcheck.eq.4 ) then ! will check volume derivatives;
		SALLOCATE( oRbc, (1:mn,0:Mvol), iRbc(1:mn,0:Mvol) )
		SALLOCATE( oZbs, (1:mn,0:Mvol), iZbs(1:mn,0:Mvol) )
		SALLOCATE( oRbs, (1:mn,0:Mvol), iRbs(1:mn,0:Mvol) )
		SALLOCATE( oZbc, (1:mn,0:Mvol), iZbc(1:mn,0:Mvol) )  
	endif ! end of if( Lcheck.eq.3 .or. Lcheck.eq.4 ) ;
  
#endif



	lfactor = psifactor(ii,vvol-1+innout) 	! this "pre-conditions" the geometrical degrees-of-freedom;

	if( Lconstraint.eq.1 .or. Lconstraint.eq.3 .or. ( Lvacuumregion .and. Lconstraint.ge.0 ) ) then ! will need to accommodate constraints;

		if(                     Lconstraint.eq.1 ) then
			iflag = -1 ; WCALL( dfp200, tr00ab, ( vvol, mn, lmns, Nt, Nz, iflag, diotadxup(0:1,-1:2,vvol) ) ) ! compute d(transform)/dx;
		endif

		if( Lvacuumregion .and. Lconstraint.ge.0 ) then
			iflag = -1 ; WCALL( dfp200, curent, ( vvol, mn,       Nt, Nz, iflag, dItGpdxtp(0:1,-1:2,vvol) ) ) ! compute d(Itor,Gpol)/dx;
		endif

		if (Lconstraint.eq.3) then		  
	
			! Derivatives of helicity multipliers
			dmupfdx(1:Mvol,vvol,1,idof,1) = zero ! The helicity multiplier is not varied (constrained)


			! Derivatives of poloidal flux. Need to solve a linear system of 5 equations - by hand for now
			! Matrix coefficients evaluation
			do pvol = 1, Mvol
				LREGION(pvol)
				WCALL(dfp200, lbpol, (pvol, 2)) ! Stores derivative in global variable Btemn
			enddo

#ifdef DEBUG
			Btemn_debug(1:mn, 0:1, 1:Mvol, 2) = zero
			Btemn_debug(1:mn, 0:1, 1:Mvol, 2) = Btemn(1:mn, 0:1, 1:Mvol)
#endif
			
			dBdmpf(1:Mvol-1,1:Mvol-1) = zero ! Initialize
			do pvol=1, Mvol-2
				dBdmpf(pvol,   pvol) =  Btemn(1, 0, pvol+1)
				dBdmpf(pvol+1, pvol) = -Btemn(1, 1, pvol+1)
			enddo
			dBdmpf(Mvol-1,Mvol-1) = Btemn(1, 0, Mvol)


#ifdef DEBUG
			do pvol=1,Mvol
                LREGION(pvol)
				WCALL(dfp200, lbpol, (pvol, 0))
			enddo
			Btemn_debug(1:mn, 0:1, 1:Mvol,  0) = Btemn(1:mn, 0:1, 1:Mvol)
#endif

			! RHS coefficients evaluation
			do pvol = vvol, vvol+1
				LREGION(pvol)
                if( pvol.eq.vvol ) then
                    dBdX%innout = 1 ! take derivative w.r.t outer interface
                else !pvol.eq.vvol+1
                    dBdX%innout = 0 ! w.r.t inner interface
                endif
				WCALL(dfp200, lbpol, (pvol, -1)) ! derivate w.r.t geometry
			enddo

#ifdef DEBUG
			Btemn_debug(1:mn, 0:1, 1:Mvol, -1) = Btemn(1:mn, 0:1, 1:Mvol)
#endif

			dBdx2(1:Mvol-1) = zero
            if( vvol.gt.1 ) then
			    dBdx2(vvol-1)   =                     - Btemn(1, 0, vvol  )
			endif
			;   dBdx2(vvol  )   = Btemn(1, 1, vvol  ) - Btemn(1, 0, vvol+1)
		    ;   dBdx2(vvol+1)   = Btemn(1, 1, vvol+1)


			! TODO: use standard library to solve linear system
			dmupfdx(1:Mvol, vvol, 2, idof, 1) = zero ! Set all components to zero excepted below
			dmupfdx(2     , vvol, 2, idof, 1) = dBdx2(1) / dBdmpf(1, 1)

			do pvol = 3, Mvol
				dmupfdx(pvol, vvol, 2, idof, 1) = ( dBdx2(pvol-1) - dmupfdx(pvol-1, vvol, 2, idof, 1) * dBdmpf(pvol-1, pvol-2) ) / dBdmpf(pvol-1, pvol-1) 
			enddo

            dmupfdx(1:Mvol, vvol, 2, idof, 1) = lfactor * dmupfdx(1:Mvol, vvol, 2, idof, 1) ! multiply by coordinate pre-conditioning factor

	else ! LocalConstraint

        dmupfdx(vvol,1,1,idof,innout) = zero    	! Prepare array
        dmupfdx(vvol,1,2,idof,innout) = zero		! Prepare array

		if( Lplasmaregion ) then

			if( Lconstraint.eq.1) then

				if( Lcoordinatesingularity ) then ! solution does not depend on dpflux, and only outer transform is a constraint;

					det = diotadxup(1,1,vvol)
					FATAL( dfp200, abs(det).lt.small, error computing derivatives of mu          wrt geometry at fixed transform )

					dmupfdx(vvol,1,1,idof,innout) = - lfactor * (                                                                        diotadxup(1,-1,vvol) ) / det
					dmupfdx(vvol,1,2,idof,innout) =   zero

				else ! if( .not.Lcoordinatesingularity ) ;

					det = diotadxup(0,1,vvol) * diotadxup(1,2,vvol) - diotadxup(0,2,vvol) * diotadxup(1,1,vvol)
					FATAL( dfp200, abs(det).lt.small, error computing derivatives of mu & dpflux wrt geometry at fixed transform )

					dmupfdx(vvol,1,1,idof,innout) = - lfactor * ( + diotadxup(1, 2,vvol) * diotadxup(0,-1,vvol) - diotadxup(0, 2,vvol) * diotadxup(1,-1,vvol) ) / det
					dmupfdx(vvol,1,2,idof,innout) = - lfactor * ( - diotadxup(1, 1,vvol) * diotadxup(0,-1,vvol) + diotadxup(0, 1,vvol) * diotadxup(1,-1,vvol) ) / det

				endif ! end of if( Lcoordinatesingularity ) ;

			endif ! end of if( Lconstraint.eq.1 ) ;
		
		else ! Vacuum region

			if    ( Lconstraint.eq.0 ) then ! THIS NEEDS ATTENTION;

				det = dItGpdxtp(0,1,vvol) * dItGpdxtp(1,2,vvol) - dItGpdxtp(0,2,vvol) * dItGpdxtp(1,1,vvol)
				FATAL( dfp200, abs(det).lt.small, error computing derivatives of dtflux & dpflux wrt geometry at fixed Itor and Gpol )

				dmupfdx(vvol,1,1,idof,innout) = - lfactor * ( + dItGpdxtp(1, 2,vvol) * dItGpdxtp(0,-1,vvol) - dItGpdxtp(0, 2,vvol) * dItGpdxtp(1,-1,vvol) ) / det
				dmupfdx(vvol,1,2,idof,innout) = - lfactor * ( - dItGpdxtp(1, 1,vvol) * dItGpdxtp(0,-1,vvol) + dItGpdxtp(0, 1,vvol) * dItGpdxtp(1,-1,vvol) ) / det

			else if( Lconstraint.eq.1 ) then

				det = diotadxup(0,1,vvol) * dItGpdxtp(1,2,vvol) - diotadxup(0,2,vvol) * dItGpdxtp(1,1,vvol)
				FATAL( dfp200, abs(det).lt.small, error computing derivatives of dtflux & dpflux wrt geometry at fixed Itor and Gpol )

				dmupfdx(vvol,1,1,idof,innout) = - lfactor * ( + dItGpdxtp(1, 2,vvol) * diotadxup(0,-1,vvol) - diotadxup(0, 2,vvol) * dItGpdxtp(1,-1,vvol) ) / det
				dmupfdx(vvol,1,2,idof,innout) = - lfactor * ( - dItGpdxtp(1, 1,vvol) * diotadxup(0,-1,vvol) + diotadxup(0, 1,vvol) * dItGpdxtp(1,-1,vvol) ) / det

			endif

		endif ! end of if( Lplasmaregion ) ;

	endif ! end of if Lconstraint.eq.3

	endif ! end of if( Lconstraint.eq.1 .or. ( Lvacuumregion .and. Lconstraint.ge.0 ) then;
        

#ifdef DEBUG

	if( Lcheck.eq.4 ) then ! check derivatives of field;

		SALLOCATE( isolution, (1:NN,-2:2), zero )
		dBdX%L = .false.

		do isymdiff = -2, 2 ! symmetric fourth-order, finite-difference used to approximate derivatives;

            if( isymdiff.eq.0 ) cycle

			iRbc(1:mn,0:Mvol) = oRbc(1:mn,0:Mvol)
			iZbs(1:mn,0:Mvol) = oZbs(1:mn,0:Mvol)
			iRbs(1:mn,0:Mvol) = oRbs(1:mn,0:Mvol)
			iZbc(1:mn,0:Mvol) = oZbc(1:mn,0:Mvol)

			if( LocalConstraint ) then
				if( issym.eq.0 .and. irz.eq.0 ) iRbc(ii,vvol-1+innout) = iRbc(ii,vvol-1+innout) + dRZ * isymdiff ! perturb geometry;
				if( issym.eq.0 .and. irz.eq.1 ) iZbs(ii,vvol-1+innout) = iZbs(ii,vvol-1+innout) + dRZ * isymdiff ! perturb geometry;
				if( issym.eq.1 .and. irz.eq.0 ) iRbs(ii,vvol-1+innout) = iRbs(ii,vvol-1+innout) + dRZ * isymdiff ! perturb geometry;
				if( issym.eq.1 .and. irz.eq.1 ) iZbc(ii,vvol-1+innout) = iZbc(ii,vvol-1+innout) + dRZ * isymdiff ! perturb geometry;
			else
				if( issym.eq.0 .and. irz.eq.0 ) iRbc(ii,vvol) = iRbc(ii,vvol) + dRZ * isymdiff ! perturb geometry;
				if( issym.eq.0 .and. irz.eq.1 ) iZbs(ii,vvol) = iZbs(ii,vvol) + dRZ * isymdiff ! perturb geometry;
				if( issym.eq.1 .and. irz.eq.0 ) iRbs(ii,vvol) = iRbs(ii,vvol) + dRZ * isymdiff ! perturb geometry;
				if( issym.eq.1 .and. irz.eq.1 ) iZbc(ii,vvol) = iZbc(ii,vvol) + dRZ * isymdiff ! perturb geometry;
			endif


			do pvol = 1,Mvol

			    LREGION(pvol)
				ll = Lrad(vvol)		! shorthand
				NN = NAdof(vvol) 	! shorthand;

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

					  WCALL( dfp200, ma00aa, ( Iquad(pvol), mn, pvol, ll ) )
					  
					  WCALL( dfp200, matrix, ( pvol, mn, ll ) )

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
				WCALL(dfp200, dfp100, (Ndofgl, Xdof, Fvec, iflag) )
			else
				Ndofgl = Mvol-1; 
				lwa = 8 * Ndofgl * Ndofgl; maxfev = 1000; nfev=0; lr=Mvol*(Mvol-1); ldfjac=Mvol-1
				ml = Mvol-2; muhybr = Mvol-2; epsfcn=1E-16; diag=0.0; mode=1; factor=0.01; nprint=1e5;	!nprint=1e5 to force last call - used for MPI communications

				Xdof(1:Mvol-1)   = dpflux(2:Mvol) + xoffset

				SALLOCATE(fjac, (1:ldfjac,1:Mvol-1), 0)
				SALLOCATE(r_deb, (1:lr), 0)

                ! TODO PROBLEM WITH MPI - not important, this is only for debug
				! Hybrid-Powell method, iterates on all poloidal fluxes to match the global constraint
				WCALL( dfp200,  hybrd, (dfp100, Ndofgl, Xdof(1:Ndofgl), Fvec(1:Ndofgl), mupftol, maxfev, ml, muhybr, epsfcn, diag(1:Ndofgl), mode, &
						  factor, nprint, ihybrd1, nfev, fjac(1:Ndofgl,1:Ndofgl), ldfjac, r_deb(1:lr), lr, qtf(1:Ndofgl), wa1(1:Ndofgl), &
						  wa2(1:Ndofgl), wa3(1:Ndofgl), wa4(1:Ndofgl)) ) 
				dpflux(2:Mvol) = Xdof(1:Ndofgl) - xoffset
		        Xdof(1:Mvol-1) = dpflux(2:Mvol) + xoffset
                WCALL( dfp200, dfp100, (Ndofgl, Xdof(1:Ndofgl), Fvec(1:Ndofgl), 6) )

				DALLOCATE(fjac)
				DALLOCATE(r_deb)

			endif




			if( LocalConstraint ) then
                LREGION(vvol)
                ll = Lrad(vvol)		! shorthand
			    NN = NAdof(vvol) 	! shorthand;

				if( Lplasmaregion ) then
					imupf_local(1:2,isymdiff) = (/     mu(vvol), dpflux(vvol) /) ! mu     and dpflux are computed for the perturbed geometry by ma02aa/mp00ac if Lconstraint=1;
				else ! if( Lvacuumregion ) ;
					imupf_local(1:2,isymdiff) = (/ dtflux(vvol), dpflux(vvol) /) ! dtflux and dpflux are computed for the perturbed geometry by ma02aa/mp00ac if Lconstraint=1;
				endif
			
                isolution(1:NN,isymdiff) = solution(vvol)%mat(1:NN,0) ! solution is computed in mp00ac, which is called by ma02aa;

			else
				if( Lplasmaregion ) then
					imupf_global(1:Mvol,1,isymdiff) = (/  mu(1:Mvol)     /) ! mu     is computed for the perturbed geometry by ma02aa/mp00ac
					imupf_global(1:Mvol,2,isymdiff) = (/  dpflux(1:Mvol) /) ! dpflux is computed for the perturbed geometry by ma02aa/mp00ac
				else ! if( Lvacuumregion ) ;
					imupf_global(1:Mvol,1,isymdiff) = (/ dtflux(1:Mvol) /) ! dtflux is computed for the perturbed geometry by ma02aa/mp00ac 
					imupf_global(1:Mvol,2,isymdiff) = (/ dpflux(1:Mvol) /) ! dpflux is computed for the perturbed geometry by ma02aa/mp00ac 
				endif

			endif
		enddo ! end of do isymdiff;


        ! Evaluate derivatives using finite differences
		if( LocalConstraint ) then
		    isolution(1:NN,0) = ( - 1 * isolution(1:NN, 2) + 8 * isolution(1:NN, 1)&
                                  - 8 * isolution(1:NN,-1) + 1 * isolution(1:NN,-2)   ) / ( 12 * dRZ )
			imupf_local(1:2,0)      = ( - 1 * imupf_local(1:2, 2)  + 8 * imupf_local(1:2, 1)&
                                        - 8 * imupf_local(1:2,-1)  + 1 * imupf_local(1:2,-2)   ) / ( 12 * dRZ )

			! TODO: see if necesseray to take absolute value - this is not used anyway...
		    solution(vvol)%mat(1:NN,-1) = abs(  solution(vvol)%mat(1:NN,-1) ) !WHY?
		    isolution(1:NN, 0) = abs( isolution(1:NN, 0) )                    !WHY?

		else
			imupf_global(1:Mvol,1:2,0)      = ( - 1 * imupf_global(1:Mvol,1:2, 2)&
                                                + 8 * imupf_global(1:Mvol,1:2, 1)&
                                                - 8 * imupf_global(1:Mvol,1:2,-1)&
                                                + 1 * imupf_global(1:Mvol,1:2,-2)   ) / ( 12 * dRZ )
		endif

		cput = GETTIME

		write(ounit,3003)
		if( LocalConstraint ) then
			!ifail = 0 ; call dlasrt( 'D', NN,  solution(vvol)%mat(1:NN,-1), ifail ) ! sorting screen output; this corrupts;
			!ifail = 0 ; call dlasrt( 'D', NN, isolution(1:NN, 0), ifail ) ! sorting screen output; this corrupts;   
			write(ounit,3003) cput-cpus, myid, vvol, im(ii), in(ii), irz, issym, innout, "finite-diff", imupf_local(1:2,0)
			write(ounit,3003) cput-cpus, myid, vvol, im(ii), in(ii), irz, issym, innout, "analytic   ", dmupfdx(vvol,1,1:2,idof,innout) / lfactor
		else
			write(ounit,3004) cput-cpus, myid, vvol, im(ii), in(ii), irz, issym, 1, "dmu finite-diff", imupf_global(1:Mvol,1,0)
			write(ounit,3004) cput-cpus, myid, vvol, im(ii), in(ii), irz, issym, 1, "dmu analytic   ", dmupfdx(1:Mvol,vvol,1,idof,1) / lfactor
			write(ounit,3004) cput-cpus, myid, vvol, im(ii), in(ii), irz, issym, 1, "dpflux finite-diff", imupf_global(1:Mvol,2,0)
			write(ounit,3004) cput-cpus, myid, vvol, im(ii), in(ii), irz, issym, 1, "dpflux analytic   ", dmupfdx(1:Mvol,vvol,2,idof,1) / lfactor
		endif

3003    format("dfp200 : ",f10.2," : ",:,"myid=",i3," ; vvol=",i2," ; (",i2,",",i3," ) ; irz=",i1," ; issym=",i1," ; innout=",i1," ; ",a18," : dmupf=",8f11.05" ;")
3004    format("dfp200 : ",f10.2," : ",:,"myid=",i3," ; vvol=",i2," ; (",i2,",",i3," ) ; irz=",i1," ; issym=",i1," ; innout=",i1," ; ",a18," : dmupf=",8f11.05" ;")


		! Re-evaluate unperturbed solution
		dBdX%L = .true.

		iRbc(1:mn,0:Mvol) = oRbc(1:mn,0:Mvol)
		iZbs(1:mn,0:Mvol) = oZbs(1:mn,0:Mvol)
		iRbs(1:mn,0:Mvol) = oRbs(1:mn,0:Mvol)
		iZbc(1:mn,0:Mvol) = oZbc(1:mn,0:Mvol)


         DALLOCATE(isolution)

!		LcomputeDerivatives = .false.
!		Lonlysolution = .true.
!		WCALL( dfp200, dforce, ( NGdof, position, force, LComputeDerivatives ))


	endif ! end of if( Lcheck.eq.4 ) ;

#endif

	vflag = 1 ! this flag instructs volume to continue even if the volume is invalid;
	WCALL( dfp200, volume, ( vvol, vflag ) ) ! compute derivative of volume; wrt to harmonic described by dBdX structure;
        
#ifdef DEBUG

	if( Lcheck.eq.3 ) then

		dvol(0) = dvolume

		cput = GETTIME
		write(ounit,1001) cput-cpus, myid, vvol, im(ii), in(ii), irz, issym, innout, "analytic", dvolume

		1001     format("dfp200 : ",f10.2," : myid=",i3," ; vvol=",i3," ; (",i3," ,",i3,") ; irz=",i2," ; issym=",i2," ; innout=",i2,&
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
			WCALL( dfp200, volume, ( vvol, vflag ) ) ! compute volume; this corrupts calculation of dvolume;

			dvol(isymdiff) = vvolume(vvol)

		enddo ! end of do isymdiff;

		evolume = abs( ( dvol(+1)-dvol(-1) ) / dRZ - dvol(0) ) ! error in finite-difference calculation and analytic derivative;

		cput = GETTIME
		write(ounit,1001) cput-cpus, myid, vvol, im(ii), in(ii), irz, issym, innout, "finite-d", ( dvol(+1)-dvol(-1) ) / dRZ, evolume

		FATAL( dfp200, evolume.gt.dRZ, unacceptable error in volume derivative )

		iRbc(1:mn,0:Mvol) = oRbc(1:mn,0:Mvol)
		iZbs(1:mn,0:Mvol) = oZbs(1:mn,0:Mvol)
		iRbs(1:mn,0:Mvol) = oRbs(1:mn,0:Mvol)
		iZbc(1:mn,0:Mvol) = oZbc(1:mn,0:Mvol)

		dBdX%L = .true.

		dvolume = dvol(0)

	endif ! end of if( Lcheck.eq.3 ) ;
        
#endif

  
#ifdef DEBUG
	if( Lcheck.eq.3 .or. Lcheck.eq.4 ) then
		DALLOCATE(oRbc)
		DALLOCATE(oZbs)
		DALLOCATE(oRbs)
		DALLOCATE(oZbc)
	endif
#endif

end subroutine evaluate_dmupfdx






!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


subroutine evaluate_dBB(vvol, idof, innout, issym, irz, ii, dAt, dAz, dBB, XX, YY, length, dRR, dZZ, dII, dLL, dPP, Ntz)

! Evaluate the derivative of the square of the magnetic field modulus. Add spectral constraint derivatives if
! required. 

! INPUT
! -----
!	vvol: 	Volume number
!   idof: 	Labels degree of freedom
!   innout: Index for loop on inner / outer interface (which interface is perturbed)
!   issym:  Index for loop on stellarator non-symmetric terms
!   irz: 	Index for loop on R or Z modes
!   ii: 	Index for loop on Fourier mode

! MODULES
! -------
  use constants, only : zero, half, one, two
  
  use numerical, only : small
  
  use fileunits, only : ounit
  
  use inputlist, only : Wmacros, Wdfp200, ext, Ntor, Igeometry, Lconstraint, &
                        gamma, adiabatic, pscale, &
                        epsilon
  
  use cputiming, only : Tdfp200
  
  use allglobal, only : ncpu, myid, cpus, &
                        Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, &
                        Mvol, &
                        iRbc, iZbs, iRbs, iZbc, & ! Fourier harmonics of geometry; vector of independent variables, position, is "unpacked" into iRbc,iZbs;
                        NOTstellsym, &
                        mn, im, in, mns, &
                        ijreal, &
                        efmn, ofmn, cfmn, sfmn, &
                        evmn, odmn, comn, simn, &
                        Nt, Nz, &
                        cosi, sini, & ! FFT workspace
                        sweight, &
                        mmpp, &
                        LGdof, &
                        vvolume, dvolume, &
                        Rij, Zij, sg, guvij, iRij, iZij, dRij, dZij, tRij, tZij, & ! Jacobian and metrics; computed in coords;
                        dFFdRZ, dBBdmp, &
                        BBweight, & ! exponential weight on force-imbalance harmonics;
                        psifactor, &
                        mn, &
                        dRodR, dRodZ, dZodR, dZodZ, dBdX

 LOCALS
!------

INTEGER		 	 		:: iocons, vvol, ideriv, id, iflag, Lcurvature, innout, issym, irz, ii, ifail, idoc, idof, Ntz
REAL         	        :: lss, DDl, MMl
REAL 					:: dAt(1:Ntz,-1:2), dAz(1:Ntz,-1:2), XX(1:Ntz), YY(1:Ntz), dBB(1:Ntz,-1:2), dII(1:Ntz), dLL(1:Ntz)
REAL					:: dPP(1:Ntz), length(1:Ntz), dRR(1:Ntz,-1:2), dZZ(1:Ntz,-1:2), constraint(1:Ntz)


do iocons = 0, 1

	if( Lconstraint.eq.1 .OR. Lconstraint.eq.3 ) then ! first, determine how B^2 varies with mu and dpflux;

		if( vvol.eq.   1 .and. iocons.eq.0 ) cycle ! fixed inner boundary (or coordinate axis); no force-balance constraints;
		if( vvol.eq.Mvol .and. iocons.eq.1 ) cycle ! fixed outer boundary                     ; no force-balance constraints;

		
		call evaluate_derivative_BB(iocons, vvol, dBB)! In a subroutine; called somewhere else when semi global constraint


		call tfft(	Nt, Nz, dBB(1:Ntz,1), dBB(1:Ntz,2), & ! derivatives of B^2 wrt mu and dpflux;
					mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), ifail )
		idoc = 0
		dBBdmp(idoc+1:idoc+mn  ,vvol,iocons,1) = efmn(1:mn) * BBweight(1:mn) ! pressure;
		dBBdmp(idoc+1:idoc+mn  ,vvol,iocons,2) = cfmn(1:mn) * BBweight(1:mn) ! pressure;
		idoc = idoc + mn   ! even;

		if( Igeometry.ge.3 ) then ! add spectral constraints; spectral constraints do not depend on mu or dpflux;
			idoc = idoc + mn-1 ! oddd;
		endif ! end of if( Igeometry.ge.3) ;

		if( NOTstellsym ) then
			dBBdmp(idoc+1:idoc+mn-1,vvol,iocons,1) = ofmn(2:mn) * BBweight(2:mn) ! pressure;
			dBBdmp(idoc+1:idoc+mn-1,vvol,iocons,2) = sfmn(2:mn) * BBweight(2:mn) ! pressure;
			idoc = idoc + mn-1 ! oddd;

			if( Igeometry.ge.3 ) then ! add spectral constraints;
				idoc = idoc + mn   ! even;
			endif ! end of if( Igeometry.ge.3) ;

		endif ! end of if( NOTstellsym) ;

	endif ! end of if( Lconstraint.eq.1 .OR. Lconstraint.eq.3 ) ;


	if( vvol.eq.   1 .and. iocons.eq.0 ) cycle ! fixed inner boundary (or coordinate axis); no constraints;
	if( vvol.eq.Mvol .and. iocons.eq.1 ) cycle ! fixed outer boundary                     ; no constraints;

	ideriv = 0 ; id = ideriv ; iflag = 0 ! why does lforce need to be re-called; need to determine which quantity (if any) has been corrupted;

	WCALL( dfp200, lforce, ( vvol, iocons, ideriv, Ntz, dAt(1:Ntz,id), dAz(1:Ntz,id), XX(1:Ntz), YY(1:Ntz), length(1:Ntz), DDl, MMl, iflag ) )

	dBB(1:Ntz,0) = 	half * ( &
							dAz(1:Ntz, 0)*dAz(1:Ntz, 0)*guvij(1:Ntz,2,2,0) &
					- two * dAz(1:Ntz, 0)*dAt(1:Ntz, 0)*guvij(1:Ntz,2,3,0) &
					+ 		dAt(1:Ntz, 0)*dAt(1:Ntz, 0)*guvij(1:Ntz,3,3,0)  ) / sg(1:Ntz,0)**2

	ideriv = -1 ; id = ideriv ; iflag = 1 ! compute derivatives of magnetic field;

    !                        lvol  iocons  ideriv  Ntz  dAt            dAz
	WCALL( dfp200, lforce, ( vvol, iocons, ideriv, Ntz, dAt(1:Ntz,id), dAz(1:Ntz,id), XX(1:Ntz), YY(1:Ntz), length(1:Ntz), DDl, MMl, iflag ) )

	lss = two * iocons - one ; Lcurvature = 4
	WCALL( dfp200, coords, ( vvol, lss, Lcurvature, Ntz, mn ) ) ! get coordinate metrics and their derivatives wrt Rj, Zj on interface;

	dBB(1:Ntz,id) = half * ( &
							dAz(1:Ntz,id)*dAz(1:Ntz, 0)*guvij(1:Ntz,2,2,0) &
					- two * dAz(1:Ntz,id)*dAt(1:Ntz, 0)*guvij(1:Ntz,2,3,0) &
					+ 		dAt(1:Ntz,id)*dAt(1:Ntz, 0)*guvij(1:Ntz,3,3,0) &
					+ 		dAz(1:Ntz, 0)*dAz(1:Ntz,id)*guvij(1:Ntz,2,2,0) &
					- two * dAz(1:Ntz, 0)*dAt(1:Ntz,id)*guvij(1:Ntz,2,3,0) &
					+ 		dAt(1:Ntz, 0)*dAt(1:Ntz,id)*guvij(1:Ntz,3,3,0) &
					+ 		dAz(1:Ntz, 0)*dAz(1:Ntz, 0)*guvij(1:Ntz,2,2,1) &
					- two * dAz(1:Ntz, 0)*dAt(1:Ntz, 0)*guvij(1:Ntz,2,3,1) &
					+ 		dAt(1:Ntz, 0)*dAt(1:Ntz, 0)*guvij(1:Ntz,3,3,1)) / sg(1:Ntz,0)**2 -	dBB(1:Ntz,0) * two * sg(1:Ntz,1) / sg(1:Ntz,0)

	FATAL( dfp200, vvolume(vvol).lt.small, shall divide by vvolume(vvol)**(gamma+one) )

	ijreal(1:Ntz) = - adiabatic(vvol) * pscale * gamma * dvolume / vvolume(vvol)**(gamma+one) + dBB(1:Ntz,-1) ! derivatives of force wrt geometry;

	dLL(1:Ntz) = zero ! either no spectral constraint, or not the appropriate interface;
	dPP(1:Ntz) = zero ! either no spectral constraint, or not the appropriate interface;

	if( Igeometry.ge.3 ) then ! spectral constraints are only required in toroidal or extended-cylindrical geometry;

		if( innout.eq.1 .and. iocons.eq.1 ) then ! include derivatives of spectral constraints;
#ifdef DEBUG
		FATAL( dfp200, abs(DDl).lt.small, divide by zero on spectral constraint )
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
				endif !matches ( irz.eq.0 )
			endif! matches ( issym.eq.0 )

		else

		dII(1:Ntz) = zero ! either no spectral constraint, or not the appropriate interface;

		endif ! end of if( innout.eq.1 .and. iocons.eq.1 ) ;

		constraint(1:Ntz) = + ( dRij(1:Ntz,vvol) * tRij(1:Ntz,vvol-1+iocons) + dZij(1:Ntz,vvol) * tZij(1:Ntz,vvol-1+iocons) ) / length(1:Ntz)

		if( iocons.eq.0 ) then ! take derivatives of constraints at inner boundary;

			if( innout.eq.0 ) then ! derivative wrt inner boundary coefficient;
				!write(ounit,'("dfp200 : " 10x " : A ; vvol="i3" ; iocons="i2" ; innout="i2" ;")') vvol, iocons, innout
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
				!write(ounit,'("dfp200 : " 10x " : B ; vvol="i3" ; iocons="i2" ; innout="i2" ;")') vvol, iocons, innout
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
				!write(ounit,'("dfp200 : " 10x " : C ; vvol="i3" ; iocons="i2" ; innout="i2" ;")') vvol, iocons, innout
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
					!write(ounit,'("dfp200 : " 10x " : dRodR(1: ,0,"i2")=",99es11.3)') ii, dRodR(1:20,0,ii)
					!write(ounit,'("dfp200 : " 10x " : dRodR(1: ,1,"i2")=",99es11.3)') ii, dRodR(1:20,1,ii)
					!write(ounit,'("dfp200 : " 10x " : dRodZ(1: ,0,"i2")=",99es11.3)') ii, dRodZ(1:20,0,ii)
					!write(ounit,'("dfp200 : " 10x " : dRodZ(1: ,1,"i2")=",99es11.3)') ii, dRodZ(1:20,1,ii)
					!write(ounit,'("dfp200 : " 10x " : dZodR(1: ,0,"i2")=",99es11.3)') ii, dZodR(1:20,0,ii)
					!write(ounit,'("dfp200 : " 10x " : dZodR(1: ,1,"i2")=",99es11.3)') ii, dZodR(1:20,1,ii)
					!write(ounit,'("dfp200 : " 10x " : dZodZ(1: ,0,"i2")=",99es11.3)') ii, dZodZ(1:20,0,ii)
					!write(ounit,'("dfp200 : " 10x " : dZodZ(1: ,1,"i2")=",99es11.3)') ii, dZodZ(1:20,1,ii)
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
						if( irz.eq.0 ) then ; dLL(1:Ntz) =    ( &   ! d/dRbs ;
															+ ( sini(1:Ntz,ii) - dRodR(1:Ntz,1,ii) ) * tRij(1:Ntz,vvol) + dRij(1:Ntz,vvol) * im(ii) * cosi(1:Ntz,ii) 	&
															+ (                - dZodR(1:Ntz,1,ii) ) * tZij(1:Ntz,vvol) 												&
															- constraint(1:Ntz) 																						&
															* (   dRij(1:Ntz,vvol) * ( sini(1:Ntz,ii)   - dRodR(1:Ntz,1,ii) )  											&
																+ dZij(1:Ntz,vvol) * (                  - dZodR(1:Ntz,1,ii) ) ) / length(1:Ntz) ) / length(1:Ntz)
						else                ; dLL(1:Ntz) = 	( &   ! d/dZbs ;
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

	call tfft( 	Nt, Nz, ijreal(1:Ntz), dII(1:Ntz), & ! recall that ijreal contains derivatives of pressure term;
				mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), ifail )


	call tfft( 	Nt, Nz, dPP(1:Ntz)   , dLL(1:Ntz), & ! recall that ijreal is probably just a dummy;
				mn, im(1:mn), in(1:mn), evmn(1:mn), odmn(1:mn), comn(1:mn), simn(1:mn), ifail )          ! evmn and odmn are available as workspace;


	FATAL( dfp200, vvol-1+innout.gt.Mvol, psifactor needs attention )

	; idoc = 0

	;  dFFdRZ(idoc+1:idoc+mn    ,vvol,iocons,idof,innout) = + efmn(1:mn    ) * psifactor(ii,vvol-1+innout) * BBweight(1:mn) ! pressure;
	; idoc = idoc + mn   ! even;
	;if( Igeometry.ge.3 ) then ! add spectral constraints;
	;  dFFdRZ(idoc+1:idoc+mn-1  ,vvol,iocons,idof,innout) = - sfmn(2:mn    ) * psifactor(ii,vvol-1+innout) * epsilon       & ! spectral condensation;
	- simn(2:mn    ) * psifactor(ii,vvol-1+innout) * sweight(vvol)   ! poloidal length constraint;
	! if( Ntor.gt.0 ) then
	!  dFFdRZ(idoc+1:idoc+Ntor  ,vvol,iocons,idof,innout) = + odmn(2:Ntor+1) * psifactor(ii,vvol-1+innout) * apsilon
	! endif
	; ;idoc = idoc + mn-1 ! oddd;
	;endif ! end of if( Igeometry.ge.3) ;

	if( NOTstellsym ) then
		; dFFdRZ(idoc+1:idoc+mn-1  ,vvol,iocons,idof,innout) = + ofmn(2:mn    ) * psifactor(ii,vvol-1+innout) * BBweight(2:mn) ! pressure;
		;idoc = idoc + mn-1 ! oddd;
		if( Igeometry.ge.3 ) then ! add spectral constraints;
			;dFFdRZ(idoc+1:idoc+mn    ,vvol,iocons,idof,innout) = - cfmn(1:mn    ) * psifactor(ii,vvol-1+innout) * epsilon       & ! spectral condensation;
			- comn(1:mn    ) * psifactor(ii,vvol-1+innout) * sweight(vvol)   ! poloidal length constraint;
			!if( Ntor.ge.0 ) then
				! dFFdRZ(idoc+1:idoc+Ntor+1,vvol,iocons,idof,innout) = + evmn(1:Ntor+1) * psifactor(ii,vvol-1+innout) * apsilon ! poloidal origin      ;
			!endif
			idoc = idoc + mn   ! even;
		endif ! end of if( Igeometry.ge.3) ;
	endif ! end of if( NOTstellsym) ;

#ifdef DEBUG
	FATAL( dfp200, idoc.ne.LGdof, counting error )
#endif
     
enddo ! end of do iocons;

end subroutine evaluate_dBB



subroutine evaluate_derivative_BB(iocons, vvol, dBB)

! INPUT
! -----
! 	iocons: 	0: inner interface, 1: outer interface
! 	vvol:		Volume number

! MODULES
! -------

  use constants, only : zero, half, one, two
  
  use numerical, only : small
  
  use fileunits, only : ounit
  
  use inputlist, only : Wmacros, Wdfp200
  
  use cputiming, only : Tdfp200
  
  use allglobal, only : ncpu, myid, cpus, &
                        Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, &
                        sg, guvij, Ntz


 LOCALS
!------

INTEGER 				:: ideriv, id, iocons, vvol, iflag
REAL 					:: dAt(1:Ntz,-1:2), dAz(1:Ntz,-1:2), XX(1:Ntz), YY(1:Ntz), dBB(1:Ntz,-1:2), length(1:Ntz), DDl, MMl


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
do ideriv = 0, 2 ; id = ideriv ! derivatives wrt helicity multiplier and differential poloidal flux;

	iflag = 1 ! lforce will only return dAt(1:Ntz,id) and dAz(1:Ntz,id);
	WCALL( dfp200, lforce, ( vvol, iocons, ideriv, Ntz, dAt(1:Ntz,id), dAz(1:Ntz,id), XX(1:Ntz), YY(1:Ntz), length(1:Ntz), DDl, MMl, iflag ) )

	if( ideriv.eq.0 ) then
		dBB(1:Ntz,id) = half * ( 	dAz(1:Ntz, 0)*dAz(1:Ntz, 0)*guvij(1:Ntz,2,2,0) - &
							  two * dAz(1:Ntz, 0)*dAt(1:Ntz, 0)*guvij(1:Ntz,2,3,0) + &
									dAt(1:Ntz, 0)*dAt(1:Ntz, 0)*guvij(1:Ntz,3,3,0) 		) / sg(1:Ntz,0)**2
	else
		dBB(1:Ntz,id) = half * (    dAz(1:Ntz,id)*dAz(1:Ntz, 0)*guvij(1:Ntz,2,2,0) - &
							  two * dAz(1:Ntz,id)*dAt(1:Ntz, 0)*guvij(1:Ntz,2,3,0) + &
									dAt(1:Ntz,id)*dAt(1:Ntz, 0)*guvij(1:Ntz,3,3,0) + &
									dAz(1:Ntz, 0)*dAz(1:Ntz,id)*guvij(1:Ntz,2,2,0) - &
							  two * dAz(1:Ntz, 0)*dAt(1:Ntz,id)*guvij(1:Ntz,2,3,0) + &
									dAt(1:Ntz, 0)*dAt(1:Ntz,id)*guvij(1:Ntz,3,3,0) ) / sg(1:Ntz,0)**2
	endif ! end of if( ideriv.gt.0 ) ;
enddo ! end of do ideriv = 0, 2;

end subroutine evaluate_derivative_BB





