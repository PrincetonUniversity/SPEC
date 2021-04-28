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
                        Lconstraint, Lcheck, LHmatrix, &
                        Lextrap, &
                        Lfreebound
  
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
                        fijreal,&
                        efmn, ofmn, cfmn, sfmn, &
                        evmn, odmn, comn, simn, &
                        Nt, Nz, &
                        cosi, sini, & ! FFT workspace;
                        dBdX, &
                        dtflux, dpflux, sweight, &
                        mmpp, &
                        dMA, dMB, dMD, dMG, &
                        Bemn, Bomn, Iomn, Iemn, Somn, Semn, &
                        LGdof, &
                        vvolume, dvolume, &
                        Rij, Zij, sg, guvij, iRij, iZij, dRij, dZij, tRij, tZij, & ! Jacobian and metrics; computed in coords;
                        diotadxup, dItGpdxtp, &
                        dFFdRZ, dBBdmp, dmupfdx, hessian, dessian, Lhessianallocated,Lhessian2Dallocated,Lhessian3Dallocated, &
                        BBweight, & ! exponential weight on force-imbalance harmonics;
                        psifactor, &
                        lmns, &
                        mn, mne, &
                        dRodR, dRodZ, dZodR, dZodZ, &
                        LocalConstraint, solution, &
                        IsMyVolume, IsMyVolumeValue, Btemn, WhichCpuID

  use typedefns
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  LOGICAL, intent(in)  :: LComputeDerivatives ! indicates whether derivatives are to be calculated;
  LOGICAL              :: LInnerVolume
  
  INTEGER              :: NN, IA, ifail, if01adf, vflag, MM, idgetrf, idgetri, Lwork, lvol, pvol
  INTEGER              :: vvol, innout, ii, jj, irz, issym, iocons, idoc, idof, imn, ll
  INTEGER              :: Lcurvature, ideriv, id
  INTEGER              :: iflag, cpu_id, cpu_id1, even_or_odd, vol_parity
  INTEGER              :: stat(MPI_STATUS_SIZE), tag, tag2, req1, req2, req3, req4

  REAL                 :: lastcpu, lss, lfactor, DDl, MMl
  REAL                 :: det
  REAL   , allocatable :: XX(:), YY(:), dBB(:,:), dII(:), dLL(:), dPP(:), length(:), dRR(:,:), dZZ(:,:), constraint(:)
  REAL   , allocatable :: ddFcol1(:), ddFcol2(:), ddFcol3(:), ddFcol4(:)


  CHARACTER            :: packorunpack 

  type(MatrixLU)  :: oBI(1:Mvol)

  BEGIN(dfp200)


#ifdef DEBUG
  if( Lcheck.eq.2 ) then
    goto 2000 ! will take no other action except a finite-difference comparison on the derivatives of the rotational-transform wrt mu and dpflux; (in dforce)
  endif
#endif

  SALLOCATE( dBB       , (1:Ntz,-1:2), zero ) ! magnetic field strength (on interfaces) in real space and derivatives;
  SALLOCATE(  XX       , (1:Ntz     ), zero )
  SALLOCATE(  YY       , (1:Ntz     ), zero )
  SALLOCATE( length    , (1:Ntz     ), zero ) ! this is calculated in lforce;

  if( LComputeDerivatives ) then
    SALLOCATE( dRR       , (1:Ntz,-1:1), zero )
    SALLOCATE( dZZ       , (1:Ntz,-1:1), zero )
    SALLOCATE( dII       , (1:Ntz     ), zero ) ! spectral constraint;
    SALLOCATE( dLL       , (1:Ntz     ), zero ) ! length   constraint;
    SALLOCATE( dPP       , (1:Ntz     ), zero ) ! poloidal constraint;
    SALLOCATE( constraint, (1:Ntz     ), zero )
    SALLOCATE( ddFcol1   , (1:Ntz     ), zero )
    SALLOCATE( ddFcol2   , (1:Ntz     ), zero )
    SALLOCATE( ddFcol3   , (1:Ntz     ), zero )
    SALLOCATE( ddFcol4   , (1:Ntz     ), zero )
  endif

  if( LocalConstraint ) then
    
  do vvol = 1, Mvol

    WCALL(dfp200, IsMyVolume, (vvol))

    if( IsMyVolumeValue .EQ. 0 ) then
      cycle
    else if( IsMyVolumeValue .EQ. -1) then
      FATAL(dfp200, .true., Unassociated volume)
    endif
                    
            
    LREGION(vvol) ! assigns Lcoordinatesingularity, Lplasmaregion, etc. ;
        
    dBdX%vol = vvol  ! Label
    ll = Lrad(vvol)  ! Shorthand
    NN = NAdof(vvol) ! shorthand;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    vflag = 1
    WCALL( dfp200, volume, ( vvol, vflag ) ) ! compute volume;
     
     !!----Hessian2D cleared-----

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    do iocons = 0, 1 ! construct field magnitude on inner and outer interfaces; inside do vvol;
        
      if( vvol.eq.1    .and. iocons.eq.0 ) cycle ! fixed inner boundary (or coordinate axis);
      if( vvol.eq.Mvol .and. iocons.eq.1 ) cycle ! fixed outer boundary                     ; there are no constraints at outer boundary;
        
      ideriv = 0 ; id = ideriv
      iflag = 0 ! XX & YY are returned by lforce; Bemn(1:mn,vvol,iocons), Iomn(1:mn,vvol) etc. are returned through global;
      WCALL( dfp200, lforce, ( vvol, iocons, ideriv, Ntz, dBB(1:Ntz,id), XX(1:Ntz), YY(1:Ntz), length(1:Ntz), DDl, MMl, iflag ) )
        
    enddo ! end of do iocons = 0, 1;
    
    if( LcomputeDerivatives ) then ! compute inverse of Beltrami matrices;

      ! Allocate some memory for storing the LU matrix and ipivot
      ! we need one extra dimension for helicity constraint, so dimension is NN+1 -> 0:NN
      SALLOCATE( oBI(vvol)%mat, (0:NN,0:NN ), zero ) ! inverse of ``original'', i.e. unperturbed, Beltrami matrix;
      SALLOCATE( oBI(vvol)%ipivot, (0:NN), zero)

      ! initialize matrices for Beltrami linear system
      call allocate_geometry_matrices(vvol, LcomputeDerivatives)
      call allocate_Beltrami_matrices(vvol, LcomputeDerivatives)
      call intghs_workspace_init(vvol)

      packorunpack = 'P'
      WCALL( dfp200, packab, ( packorunpack, vvol, NN, solution(1:NN,0), 0 ) ) ! packing, put the solutions back to the solution matrix;
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

      ! Invert beltrami matrix. Required for matrix perturbation theory
      call get_LU_Beltrami_matrices(vvol, oBI(vvol), NN)
   
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

              ! Evaluate derivatives of field w.r.t geometry
              call get_perturbed_solution(vvol, oBI(vvol), NN)


              ! Helicity multiplier and poloidal flux derivatives w.r.t geometry
              call evaluate_dmupfdx(innout, idof, ii, issym, irz)

              !-----Hessian 2D cleared

            !if (LHmatrix .and. Lhessian2Dallocated) then
             !   call hessian_dFFdRZ(vvol, idof, innout, issym, irz, ii, dBB, XX, YY, length, dRR, dZZ, dII, dLL, dPP, Ntz)
            !else
                ! Evaluate derivatives of B square 
              if(Lhessianallocated) then
                call evaluate_dBB(vvol, idof, innout, issym, irz, ii, dBB, XX, YY, length, dRR, dZZ, dII, dLL, dPP, Ntz)
              endif

                if (LHmatrix .and. Lhessian3Dallocated .and. Igeometry.ge.3) then
                    call hessian3D_dFFdRZ(vvol, idof, innout, issym, irz, ii, dBB, XX, YY, length, dRR, dZZ, dII, dLL, dPP, Ntz)
                endif

                !if (LHmatrix .and. Lhessian2Dallocated .and. Igeometry.eq.2) then
                !    call hessian_dFFdRZ(vvol, idof, innout, issym, irz, ii, dBB, XX, YY, length, dRR, dZZ, dII, dLL, dPP, Ntz)
                !endif

            enddo ! matches do innout;
          enddo ! matches do issym;
        enddo ! matches do irz;
      enddo ! matches do ii;    
            
      dBdX%L = .false. ! probably not needed, but included anyway;

      call intghs_workspace_destroy()
      call deallocate_Beltrami_matrices(LcomputeDerivatives)
      call deallocate_geometry_matrices(LcomputeDerivatives)

      DALLOCATE(oBI(vvol)%mat)
      DALLOCATE(oBI(vvol)%ipivot)
       
    endif ! end of if( LComputeDerivatives ) ;

  enddo ! matches do vvol = 1, Mvol

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

else ! CASE SEMI GLOBAL CONSTRAINT

     do vvol = 1, Mvol
         WCALL(dfp200, IsMyVolume, (vvol))

         if( IsMyVolumeValue .EQ. 0 ) then
             cycle
         else if( IsMyVolumeValue .EQ. -1) then
             FATAL(dfp200, .true., Unassociated volume)
         endif
                    
            
         LREGION(vvol) ! assigns Lcoordinatesingularity, Lplasmaregion, etc. ;
         ll = Lrad(vvol)  ! Shorthand
         NN = NAdof(vvol) ! shorthand;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

         vflag = 1
         WCALL( dfp200, volume, ( vvol, vflag ) ) ! compute volume;

         do iocons = 0, 1 ! construct field magnitude on inner and outer interfaces; inside do vvol;
        
             if( vvol.eq.1    .and. iocons.eq.0 ) cycle ! fixed inner boundary (or coordinate axis);
             if( vvol.eq.Mvol .and. iocons.eq.1 ) cycle ! fixed outer boundary                     ; there are no constraints at outer boundary;
            
             ideriv = 0 ; id = ideriv
             iflag = 0 ! dAt, dAz, XX & YY are returned by lforce; Bemn(1:mn,vvol,iocons), Iomn(1:mn,vvol) etc. are returned through global;
             WCALL( dfp200, lforce, ( vvol, iocons, ideriv, Ntz, dBB(1:Ntz,id), XX(1:Ntz), YY(1:Ntz), length(1:Ntz), DDl, MMl, iflag ) )
            
         enddo ! end of do iocons = 0, 1;
     enddo



    if( LComputeDerivatives ) then
        dBdX%L = .true. ! will need derivatives;
                        
        ! Each volume will need solution and derivative information in vacuum region. Thus, broadcast!
        ! Derivative of w.r.t geometry only required if plasma interface is perturbed - this is 
        ! calculated below
        if( Lfreebound.eq.1 ) then
            do ideriv = 0, 2
                do ii = 1, mn  

                    call WhichCpuID(Mvol, cpu_id)
                    
                    RlBCAST( Ate(Mvol,ideriv,ii)%s(0:Lrad(Mvol)), Lrad(Mvol)+1, cpu_id )
                    RlBCAST( Aze(Mvol,ideriv,ii)%s(0:Lrad(Mvol)), Lrad(Mvol)+1, cpu_id )
                enddo
             enddo  
        endif



        ! Evaluation of froce gradient in case of semi-global constraint (for now, this is specific to Lconstraint=3)
        ! Loop over all geometrical degrees of freedom. Thi means we loop over
        ! - each interface (vvol)
        ! - each Fourier mode (ii)
        ! - on R or Z (irz)
        ! - on symmetric and non stellarator symmetric terms (issym)

        ! Then the perturbed solution is evaluated on volumes sharing the interface vvol (volume vvol and vvol+1).
        ! i.e. relevant routines are run with lvol=vvol, innout=1 and lvol=vvol+1, innout=0


        ! First invert Beltrami matrices and store them in OBI
        do vvol = 1, Mvol
            LREGION(vvol) ! assigns Lcoordinatesingularity, Lplasmaregion, etc. ; TODO: maybe not necessary, remove
            NN = NAdof(vvol) ! shorthand;
            ll = Lrad(vvol)

            SALLOCATE( oBI(vvol)%mat, (0:NN,0:NN ), zero ) ! inverse of ``original'', i.e. unperturbed, Beltrami matrix;
            SALLOCATE( oBI(vvol)%ipivot, (0:NN), zero)

            ! Parallelization
            WCALL(dfp200, IsMyVolume, (vvol))
            if( IsMyVolumeValue .EQ. 0 ) then
                cycle
            else if( IsMyVolumeValue .EQ. -1) then
                FATAL(dfp200, .true., Unassociated volume)
            endif


            ! Invert LHS of Beltrami system and store it in oBI. This will be used to 
            ! evaluate derivatives of solution.
            ! TODO: Storage of OBi is not optimal. Find a way around?

            call allocate_geometry_matrices(vvol, LcomputeDerivatives)
            call allocate_Beltrami_matrices(vvol, LcomputeDerivatives)

            call get_LU_beltrami_matrices(vvol, oBI(vvol), NN)

            call deallocate_Beltrami_matrices(LcomputeDerivatives)
            call deallocate_geometry_matrices(LcomputeDerivatives)
        enddo ! end of vvol = 1, Mvol
            
        ! Broadcast oBI to all CPU
        ! TODO: THIS is bad! oBI can be very large... Change parallelization?
        !       Another possible parallelization: Each CPU keep derivatives info about its solution, 
        !                                         CPU1 build linear system and solve it?
        if( ncpu .gt. 1 ) then
            do vvol = 1, Mvol
                NN = NAdof(vvol)
                call WhichCpuID(vvol, cpu_id)
                call MPI_BCAST( oBI(vvol)%mat(0:NN,0:NN), (NN+1)**2, MPI_DOUBLE_PRECISION, cpu_id , MPI_COMM_WORLD , ierr)
                call MPI_BCAST( oBI(vvol)%ipivot(0:NN)  ,  NN+1    , MPI_INTEGER         , cpu_id , MPI_COMM_WORLD , ierr)
            enddo
        endif

        ! Loop on perturbed interfaces
        do even_or_odd = 0, 1 ! First loop on even interfaces perturbation, then on odd interfaces. This allow efficient parallelization
        
            do vvol = 1, Mvol-1 !labels which interface is perturbed

            
            ! Parallelization - parallelization on perturbed interface and not on volume (outermost loop)
            vol_parity = MODULO(vvol,2)
            if( (vol_parity.eq.0 ) .and. (even_or_odd.eq.1) ) cycle ! even_or_odd=1 thus perturb only odd  interfaces
            if( (vol_parity.eq.1 ) .and. (even_or_odd.eq.0) ) cycle ! even_or_odd=0 thus perturb only even interfaces
            
            WCALL(dfp200, IsMyVolume, (vvol))

            if( IsMyVolumeValue.EQ.0 ) then ! This CPU does not deal with interface's inner volume
    
                WCALL(dfp200, IsMyVolume, (vvol+1))
    
                if( IsMyVolumeValue.eq.0 ) then ! This CPU does not deal with interface's outer volume either - cycle
                    cycle
                else if( IsMyVolumeValue.eq.-1 ) then
                    FATAL(dfp200, .true., Unassociated volume)
                else
                    LInnerVolume = .false.
                endif
    
            else if( IsMyVolumeValue.EQ.-1 ) then
                FATAL(dfp200, .true., Unassociated volume)
            else
                LinnerVolume = .true.    
            endif


            dBdX%vol = vvol     ! Perturbed interface
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

                        ! Only interested in one interface perturbation
                        if( ncpu.gt.1 ) then
                            if(      LinnerVolume .and. (lvol.eq.vvol+1) ) cycle
                            if( .not.LinnerVolume .and. (lvol.eq.vvol  ) ) cycle
                        endif

                        ! Set up perturbation information
                        if( lvol.eq.vvol   ) innout=1 ! Perturb w.r.t outer interface
                        if( lvol.eq.vvol+1 ) innout=0 ! Perturb w.r.t inner interface

                        dBdX%innout = innout
                        dBdX%L      = .true.

                        ! Set up volume information                        
                        LREGION(lvol) ! assigns Lcoordinatesingularity, Lplasmaregion, etc. ;
                        NN = NAdof(lvol)
                       
                        ! Allocate memory. This cannot be moved outside due to NN and ll dependence on volume.
                        call allocate_geometry_matrices(lvol, LcomputeDerivatives)
                        call allocate_Beltrami_matrices(lvol, LcomputeDerivatives)

                        ! Get derivative of vector potential w.r.t geometry. Matrix perturbation theory.
                        call intghs_workspace_init(lvol)
                        call get_perturbed_solution(lvol, oBI(lvol), NN)
                        call intghs_workspace_destroy()
                    
                        ! Free memory
                        call deallocate_Beltrami_matrices(LcomputeDerivatives)
                        call deallocate_geometry_matrices(LcomputeDerivatives)

                    enddo ! end of do lvol = vvol, vvol+1


                    ! Now volumes neighbouring the interface share perturbed solution
                    call WhichCpuID(vvol  , cpu_id ) 
                    call WhichCpuID(vvol+1, cpu_id1) 


                    ! TODO IMPROVE MPI COMMUNICATIONS. Could directly send dAt and dAz?
                    ! Gather everything in inner volume
                    if( ncpu.gt. 1) then
                        if( LinnerVolume ) then    
                            do jj = 1, mn  
                                tag  = 1 ! Tags for MPI communications
                                tag2 = 2

                                call MPI_RECV(Ate(vvol+1,-1,jj)%s(0:Lrad(vvol+1)), Lrad(vvol+1)+1, MPI_DOUBLE_PRECISION, cpu_id1, tag , MPI_COMM_WORLD, stat, ierr)
                                call MPI_RECV(Aze(vvol+1,-1,jj)%s(0:Lrad(vvol+1)), Lrad(vvol+1)+1, MPI_DOUBLE_PRECISION, cpu_id1, tag2, MPI_COMM_WORLD, stat, ierr)
                            enddo

                            ! Non-stellarator symmetric terms
                            if( NOTstellsym ) then
                                do jj = 1, mn  
                                    tag  = 3
                                    tag2 = 4

                                    call MPI_RECV(Ato(vvol+1,-1,jj)%s(0:Lrad(vvol+1)), Lrad(vvol+1)+1, MPI_DOUBLE_PRECISION, cpu_id1, tag , MPI_COMM_WORLD, stat, ierr)
                                    call MPI_RECV(Azo(vvol+1,-1,jj)%s(0:Lrad(vvol+1)), Lrad(vvol+1)+1, MPI_DOUBLE_PRECISION, cpu_id1, tag2, MPI_COMM_WORLD, stat, ierr)
                                enddo
                            endif

                        else
                            do jj = 1, mn  
                                tag  = 1
                                tag2 = 2

                                call MPI_SEND(Ate(vvol+1,-1,jj)%s(0:Lrad(vvol+1)), Lrad(vvol+1)+1, MPI_DOUBLE_PRECISION, cpu_id , tag , MPI_COMM_WORLD, ierr)
                                call MPI_SEND(Aze(vvol+1,-1,jj)%s(0:Lrad(vvol+1)), Lrad(vvol+1)+1, MPI_DOUBLE_PRECISION, cpu_id , tag2, MPI_COMM_WORLD, ierr)
                            enddo

                            if( NOTstellsym ) then
                                do jj = 1, mn  
                                    tag  = 3
                                    tag2 = 4

                                    call MPI_SEND(Ato(vvol+1,-1,jj)%s(0:Lrad(vvol+1)), Lrad(vvol+1)+1, MPI_DOUBLE_PRECISION, cpu_id , tag , MPI_COMM_WORLD, ierr)
                                    call MPI_SEND(Azo(vvol+1,-1,jj)%s(0:Lrad(vvol+1)), Lrad(vvol+1)+1, MPI_DOUBLE_PRECISION, cpu_id , tag2, MPI_COMM_WORLD, ierr)
                                enddo
                            endif
                        endif
                    endif

                    ! At this point, the inverted original Beltrami matrix oBI and the perturbed solution for each volume 
                    ! neighboring the perturbed interface in inner volume cpu. We now compute the derivatives of mu and 
                    ! psip w.r.t the position and dBB.
                    if( LinnerVolume) then
                        ! Helicity multiplier and poloidal flux derivatives
                        call evaluate_dmupfdx(1, idof, ii, issym, irz)
                    else
                        ! Only inner volume computes the derivatives in all volumes (it will be broadcasted later)
                        dmupfdx(1:Mvol, vvol, 1:2, idof, 1) = zero 
                    endif

                    do lvol = vvol, vvol+1
                        WCALL(dfp200, IsMyVolume, (lvol))

                        if( IsMyVolumeValue .EQ. 0 ) then
                            cycle
                        else if( IsMyVolumeValue .EQ. -1) then
                            FATAL(dfp200, .true., Unassociated volume)
                        endif

                        if( lvol.eq.vvol ) then
                            innout      = 1
                        else
                            innout      = 0
                        endif
                        dBdX%innout = innout                        
                        
                        LREGION(lvol) ! assigns Lcoordinatesingularity, Lplasmaregion, etc. ;

                        ! EVALUATE dBB
                        call evaluate_dBB(lvol, idof, innout, issym, irz, ii, dBB, XX, YY, length, dRR, dZZ, dII, dLL, dPP, Ntz)

                    enddo     ! matches do lvol = vvol, vvol+1 

                enddo ! matches do issym;
            enddo ! matches do irz;
        enddo ! matches do ii;
      enddo ! matches do vvol;
    enddo ! matches do even_or_odd;

    ! Free memory                          
    do vvol = 1, Mvol
        DALLOCATE(oBi(vvol)%mat)
        DALLOCATE(oBI(vvol)%ipivot)
    enddo

    dBdX%L = .false. ! probably not needed, but included anyway;

    ! We know need to broadcast the vectors dmupfdx and dFFdRZ and dBBdmp
    do vvol = 1, Mvol
        call WhichCpuID(vvol, cpu_id)

        if( vvol.ne.Mvol ) then
            call MPI_BCAST( dmupfdx(1:Mvol, vvol ,1:2, 1:LGdof,   1), Mvol*2*LGdof  , MPI_DOUBLE_PRECISION, cpu_id, MPI_COMM_WORLD, ierr )
        endif
    
        call MPI_BCAST( dFFdRZ(1:LGdof, 0:1, 1:LGdof, 0:1, vvol), 2*2*(LGdof**2), MPI_DOUBLE_PRECISION, cpu_id, MPI_COMM_WORLD, ierr )
        call MPI_BCAST( dBBdmp(1:LGdof, vvol ,0:1, 1:2         ), 2*2*LGdof, MPI_DOUBLE_PRECISION, cpu_id, MPI_COMM_WORLD, ierr )
    enddo

    endif ! End of if( LComputeDerivatives ) 



  endif ! End of if( LocalConstraint )

  if( LcomputeDerivatives ) then
      DALLOCATE(constraint)
      DALLOCATE(dPP)
      DALLOCATE(dLL)
      DALLOCATE(dII)
      DALLOCATE(dZZ)
      DALLOCATE(dRR)
      DALLOCATE(ddFcol1)
  DALLOCATE(ddFcol2)
  DALLOCATE(ddFcol3)
  DALLOCATE(ddFcol4)
  endif

  DALLOCATE(dBB)
  DALLOCATE( XX) ! spectral constraints; not used;
  DALLOCATE( YY)
  DALLOCATE(length)


2000 continue

  RETURN(dfp200)

end subroutine dfp200



!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
!                                                                LOCAL SUBROUTINES
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine get_LU_beltrami_matrices(vvol, oBI, NN)

! Evaluate the LU factorization of Beltrami matrices and store the original one in oBI.

! INPUT
! -----
!     vvol:     Volume number
!   NN:        is equal to NAdof(vvol), is a shorthand.


! MODULES
! -------

  use constants, only :   zero, half, one, two

  use fileunits, only :   ounit

  use cputiming, only :   Tdfp200

  use inputlist, only :   Wmacros, Wdfp200, Wdforce, Lrad, mu, Lconstraint

  use allglobal, only :   ncpu, myid, cpus, &
                          Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, &
                          Nt, Nz, &
                          dBdX, &
                          dMA, dMB, dMD, dMG, &
                          mn, mne, Iquad, solution, Lsavedguvij

  use typedefns

  LOCALS
! ------
  TYPE(MatrixLU), intent(inout) :: oBI
  INTEGER, intent(in)    :: vvol, NN

  INTEGER             :: IA, MM, LDA, Lwork, ll
  INTEGER             :: idgetrf, idgetri
  REAL                :: lastcpu
  REAL, allocatable   :: work(:)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  lastcpu = GETTIME

  ll = Lrad(vvol)

  ! Evaluate dMa, dMD
  dBdX%L = .false.

  Lsavedguvij = .false.

  WCALL( dfp200, ma00aa, (Iquad(vvol), mn, vvol, ll) )
  WCALL( dfp200, matrix, (vvol, mn, ll) )

  lastcpu = GETTIME
    
  if(Lconstraint .eq. 2) then  ! for helicity constraint
    
    !dMA(1:NN,1:NN) = dMA(1:NN,1:NN) - mu(vvol) * dMD(1:NN,1:NN) ! this corrupts dMA, but dMA is no longer used;
    call DAXPY((NN+1)*(NN+1), -mu(vvol), dMD, 1, dMA, 1) ! BLAS version; 24 Jul 2019
    dMA(0,0)       = zero
      
    dMA(1:NN,0)    = -matmul(dMD(1:NN,1:NN),solution(1:NN,0))
    dMA(0,1:NN)    = dMA(1:NN,0) ! This is the transpose of the above
      
    IA = NN + 1 + 1
    MM = NN ; LDA = NN+1 ; Lwork = NN+1
    idgetrf = 1 ; call DGETRF( MM+1, NN+1, dMA(0:LDA-1,0:NN), LDA, oBI%ipivot, idgetrf )

    cput = GETTIME
    select case( idgetrf ) !                                                                     0123456789012345678
    case(  :-1 ) ;               write(ounit,1010) cput-cpus, myid, vvol, cput-lastcpu, idgetrf, "input error;      "
    case(  0   ) ; if( Wdforce ) write(ounit,1010) cput-cpus, myid, vvol, cput-lastcpu, idgetrf, "success;          "
    case( 1:   ) ;               write(ounit,1010) cput-cpus, myid, vvol, cput-lastcpu, idgetrf, "singular;         "
    case default ;               FATAL( dforce, .true., illegal ifail returned from F07ADF )
    end select    

    call DCOPY((1+NN)*(1+NN), dMA, 1, oBI%mat, 1)            ! BLAS version; 24 Jul 2019

  else ! for other constraints

    !dMA(0:NN-1,1:NN) = dMA(1:NN,1:NN) - mu(vvol) * dMD(1:NN,1:NN) ! this corrupts dMA, but dMA is no longer used;
    !dMA(  NN  ,1:NN) = zero
    call DAXPY((NN+1)*(NN+1), -mu(vvol), dMD, 1, dMA, 1) ! BLAS version; 24 Jul 2019

    !dMD(1:NN  ,1:NN) = dMA(0:NN-1,1:NN) ! copy of original matrix; this is used below;
    call DCOPY((NN+1)*(NN+1), dMA, 1, dMD, 1)            ! BLAS version; 24 Jul 2019

    IA = NN + 1
    
    MM = NN ; LDA = NN + 1 ; Lwork = NN
    
    idgetrf = 1 ; call DGETRF( MM, NN, dMA(1,1), LDA, oBI%ipivot(1:NN), idgetrf )
    
    cput = GETTIME
    select case( idgetrf ) !                                                                     0123456789012345678
    case(  :-1 ) ;               write(ounit,1010) cput-cpus, myid, vvol, cput-lastcpu, idgetrf, "input error;      "
    case(  0   ) ; if( Wdforce ) write(ounit,1010) cput-cpus, myid, vvol, cput-lastcpu, idgetrf, "success;          "
    case( 1:   ) ;               write(ounit,1010) cput-cpus, myid, vvol, cput-lastcpu, idgetrf, "singular;         "
    case default ;               FATAL( dforce, .true., illegal ifail returned from F07ADF )
    end select
  
    oBI%mat(1:NN,1:NN) = dMA(1:NN,1:NN)
    
  endif

1011 format("dforce : ",f10.2," : myid=",i3," ; vvol=",i3," ; called DGETRI ; time=",f10.2,"s ; inverse of Beltrami matrix; idgetrf=",i2," ; ",a18)
1010 format("dforce : ",f10.2," : myid=",i3," ; vvol=",i3," ; called DGETRF ; time=",f10.2,"s ; LU factorization of matrix; idgetrf=",i2," ; ",a18)

end subroutine get_LU_beltrami_matrices

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine get_perturbed_solution(vvol, oBI, NN)

! This routine evaluates the value of the magnetic field once the interface is perturbed using matrix perturbation theory.
! Separated from the main dfp200 core to allow local and semi-global constraints to be prescribed
!
! Attention: Solution is perturbed for a given degree of freedom (information stored in dBdX)

  use constants, only :   zero, half, one, two

  use fileunits, only :   ounit

  use cputiming, only :   Tdfp200

  use inputlist, only :   Wmacros, Wdfp200, Lrad, mu, Lconstraint

  use allglobal, only :   ncpu, myid, cpus, &
                          mn, Iquad, NAdof, &
                          dMA, dMB, dMD, dMG, solution, &
                          dtflux, dpflux, dBdX

  use typedefns

 LOCALS
!------

  INTEGER, intent(in)     :: vvol, NN
  TYPE(MatrixLU),intent(inout) :: oBI

  INTEGER                 :: ideriv, ll, idgetrs
  REAL                    :: dpsi(1:2), work(1:NN+1), rhs(0:NN), dVA(0:NN), dVD(0:NN)
  CHARACTER               :: packorunpack

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
    
  ll = Lrad(vvol)  ! Shorthand
  dBdX%L = .true.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  WCALL( dfp200, intghs, ( Iquad(vvol), mn, vvol, ll, 0 ) )

  WCALL( dfp200, mtrxhs, ( vvol, mn, ll, dVA, dVD, 0) )

  rhs(0)    = zero
  rhs(1:NN) = -dVA(1:NN)

  if (Lconstraint .eq. 2) then
    work(1:NN+1) = rhs(0:NN)
    !work(1:NN+1)  =  matmul( oBI(0:NN,0:NN), rhs(0:NN)) ! original version
    call DGETRS('N',NN+1,1,oBI%mat(0,0),NN+1,oBI%ipivot(0:NN),work(1),NN+1,idgetrs) ! Change to DGETRS; 22 Jul 19
    solution(1:NN,-1) = work(2:NN+1)
  else
    !solution(1:NN,-1) = matmul( oBI(1:NN,1:NN), rhs(1:NN) ) ! original version
    solution(1:NN,-1) = rhs(1:NN)
    call DGETRS('N',NN,1,oBI%mat(1,1),NN+1,oBI%ipivot(1:NN),solution(1,-1),NN,idgetrs) ! Change to DGETRS; 22 Jul 19
  endif

  ! Unpack derivatives of solution
  packorunpack = 'U'
  WCALL( dfp200, packab,( packorunpack, vvol, NN,  solution(1:NN,-1), -1 ) ) ! derivatives placed in Ate(lvol,ideriv,1:mn)%s(0:Lrad),

end subroutine get_perturbed_solution





!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine evaluate_dmupfdx(innout, idof, ii, issym, irz)

! Evaluate mu and psip derivatives and store them in dmupfdx.
! If debug and Lcheck=4, compare with finite difference approximation
!
! INPUT:
! ------
!     innout: inner or outer boundary. Not used in case of Lconstraint=3 since inside loop
!     idof:      ???
!     ii:     index of Fourier harmonics, ii=1:mn
!     issym:     loop over stellerator symmetric and non-symmetric terms
!     irz:     loop over R or Z

! MODULES:
! -------
    use constants, only :   zero, half, one, two

    use fileunits, only :   ounit

    use numerical, only :   small

    use cputiming, only :   Tdfp200

    use inputlist, only :   Wmacros, Wdfp200, Lrad, mu, Lconstraint, Lcheck, Nvol, Igeometry, mupftol, Lfreebound, dRZ

    use allglobal, only :   ncpu, myid, cpus, &
                            Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, &
                            Mvol, Iquad, NGdof, &
                            iRbc, iZbs, iRbs, iZbc, & ! Fourier harmonics of geometry; vector of independent variables, position, is "unpacked" into iRbc,iZbs;
                            NAdof, &
                            mn, im, in, mns, &
                            Ate, Aze, Ato, Azo, &     ! only required for debugging;
                            Nt, Nz, dBdX, &
                            dtflux, dpflux, sweight, &
                            Rij, Zij, &
                            diotadxup, dItGpdxtp, dmupfdx, &
                            psifactor, &
                            lmns, &
                            mn, mne, &
                            LocalConstraint, &
                            vvolume, dvolume, &
                            IsMyVolume, IsMyVolumeValue, &
                            Btemn, xoffset, &
                            IPdtdPf, &
                            dMA, dMB, dMG, dMD


  LOCALS:
! -------

    INTEGER             ::  vvol, innout, idof, iflag, ii, issym, irz, ll, NN, ifail
    INTEGER             ::  vflag, N, iwork(1:Nvol-1), idgesvx, pvol, order, IDGESV
    INTEGER             ::  iocons
    INTEGER, allocatable::  IPIV(:)
    REAL                ::  det, lfactor, Bt00(1:Mvol, 0:1, -1:2)
    REAL                ::  R(1:Nvol-1), C(1:Nvol-1), work(1:4*Nvol-4), ferr, berr, rcond, tmp(2:Nvol)
    LOGICAL             ::  Lonlysolution, LcomputeDerivatives
    REAL, allocatable   ::  dBdmpf(:,:), dBdx2(:)

#ifdef DEBUG
    INTEGER             :: isymdiff, lr, ml, mode
    INTEGER             :: jj, tdoc, idoc, tdof, jdof, imn, Ndofgl
    REAL                :: dvol(-1:+1), evolume, imupf_global(1:Mvol,1:2,-2:2), imupf_local(1:2,-2:2), factor, Btemn_debug(1:mn, 0:1, 1:Mvol, -1:2)
    REAL                :: position(0:NGdof), force(0:NGdof)
    REAL                :: Fdof(1:Mvol-1), Xdof(1:Mvol-1)
    REAL, allocatable   :: fjac(:, :), r_deb(:), Fvec(:), dpfluxout(:)
    REAL, allocatable   :: oRbc(:,:), oZbs(:,:), oRbs(:,:), oZbc(:,:) ! original geometry;

    CHARACTER           :: packorunpack

    EXTERNAL            :: dfp100
#endif


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    vvol = dBdX%vol     ! shorthand

    ll = Lrad(vvol)        ! shorthand
    NN = NAdof(vvol)     ! shorthand;

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

    lfactor = psifactor(ii,vvol-1+innout)     ! this "pre-conditions" the geometrical degrees-of-freedom;

    if( Lconstraint.eq.1 .or. Lconstraint.eq.3 .or. ( Lvacuumregion .and. Lconstraint.ge.0 ) ) then ! will need to accommodate constraints
        
        if (Lconstraint.eq.3) then    
    
            ! In case of freeboundary, there is an additional equation in the linear system, related to the poloidal linking current.
            if( Lfreebound.eq.1 ) then
                order = Mvol
            else
                order = Mvol-1
            endif
            
            SALLOCATE(dBdmpf , ( 1:order, 1:order ), zero)
            SALLOCATE(dBdx2  , ( 1:order ), zero)
            SALLOCATE(IPIV   , ( 1:order ), zero)


            ! Derivatives of helicity multipliers
            dmupfdx(1:Nvol,vvol,1,idof,1) = zero ! The helicity multiplier is not varied (constrained). However dtflux is varied in vacuum region, if free boundary

            ! Derivatives of poloidal flux. Need to solve a linear system of 5 equations - by hand for now
            ! Matrix coefficients evaluation
            do pvol = 1, Mvol
                LREGION(pvol)

                do iocons = 0, 1
            if( ( Lcoordinatesingularity .and. iocons.eq.0 ) .or. ( Lvacuumregion .and. iocons.eq.1 ) ) cycle
                    WCALL(dfp200, lbpol, (pvol, Bt00(1:Mvol, 0:1, -1:2), 2, iocons)) ! Stores derivative in global variable Btemn
                enddo
#ifdef DEBUG
                if( .false. ) then
                    write(ounit, 8375) myid, dBdX%vol, dBdX%innout, 2, pvol, Bt00(pvol, 0:1, 2)
                endif
#endif
            enddo
            
            dBdmpf(1:order,1:order) = zero ! Initialize. TODO: useless?
            do pvol=1, Mvol-2
                dBdmpf(pvol,   pvol) =  Bt00(pvol+1, 0, 2)
                dBdmpf(pvol+1, pvol) = -Bt00(pvol+1, 1, 2)
            enddo
            dBdmpf(Mvol-1,Mvol-1) = Bt00(Mvol, 0, 2)


!             do pvol=1,Mvol
!                 LREGION(pvol)

!                 do iocons = 0, 1
!                     WCALL(dfp200, lbpol, (pvol, Bt00(1:Mvol, 0:1, -1:2), 0, iocons))
!                 enddo
! #ifdef DEBUG
!                 if( .false. ) then
!                     write(ounit, 8375) myid, dBdX%vol, dBdX%innout, 0, pvol, Bt00(pvol, 0:1, 0)
!                 endif
! #endif
!             enddo

            ! RHS coefficients evaluation
            do pvol = vvol, vvol+1
                LREGION(pvol)
                if( pvol.eq.vvol ) then
                    dBdX%innout = 1 ! take derivative w.r.t outer interface
                else !pvol.eq.vvol+1
                    dBdX%innout = 0 ! w.r.t inner interface
                endif

                do iocons = 0, 1
           if( ( Lcoordinatesingularity .and. iocons.eq.0 ) .or. ( Lvacuumregion .and. iocons.eq.1 ) ) cycle
                    WCALL(dfp200, lbpol, (pvol, Bt00(1:Mvol, 0:1, -1:2), -1, iocons)) ! derivate w.r.t geometry
                enddo

#ifdef DEBUG
                if( .false. ) then
                    write(ounit, 8375) myid, dBdX%vol, dBdX%innout, -1, pvol, Bt00(pvol, 0:1, -1)        
8375                format("dfp200  : myid=",i3, ", vvol=",i3,", innout=", i3 ,", ideriv=", i3, ", lvol=",i3,";  Bt00=",2f10.6)
                endif
#endif
            enddo
            
            dBdx2(1:Mvol-1) = zero
            if( vvol.gt.1 ) then
                dBdx2(vvol-1)   =                 - Bt00(vvol,   0, -1)
            endif
            ;   dBdx2(vvol  )   = Bt00(vvol  , 1, -1) - Bt00(vvol+1, 0, -1)
            if (vvol.lt.Mvol-1) then
            ; dBdx2(vvol+1) = Bt00(vvol+1, 1, -1)
            endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
            if( Lfreebound.eq.1 ) then ! Need to modify last two equations

                ! Set all last column to zero - maybe not necessary ?
                dBdmpf(1:Mvol, Mvol  ) = zero

                ! Get derivatives of B_theta w.r.t the toroidal flux in vacuum region
                iocons = 0
                WCALL(dfp200, lbpol, (Mvol, Bt00(1:Mvol, 0:1, -1:2), 1, iocons))

                ! compute d(Itor,Gpol)/dpsip and d(Itor,Gpol)/dpsit 
                ! TODO: this should already be evaluated in mp00ac...
                ! TODO: THIS COULD BE MOVED OUTSIDE THE LOOPS
                iflag =  2 ; WCALL( dfp200, curent, ( Mvol, mn, Nt, Nz, iflag, dItGpdxtp(0:1,-1:2,Mvol) ) ) 

                dBdmpf(Mvol-1, Mvol  ) =  Bt00(Mvol, 0, 1)         !dBdpsit
                dBdmpf(Mvol  , Mvol-1) =  dItGpdxtp( 1, 2, Mvol)    !dIpdpsip
                dBdmpf(Mvol  , Mvol  ) =  dItGpdxtp( 1, 1, Mvol)    !dIpdpsit

                ! compute d(Itor,Gpol)/dx
                if( vvol.eq.Mvol-1 ) then ! Plasma interface is perturbed
                    iflag = -1 ; WCALL( dfp200, curent, ( Mvol, mn, Nt, Nz, iflag, dItGpdxtp(0:1,-1:2,Mvol) ) ) 
                    dBdx2( Mvol ) = -dItGpdxtp( 1,-1, Mvol)    !-dIpdxj    
                else ! Inner interface is perturbed
                    dBdx2( Mvol ) = zero
                endif
#ifdef DEBUG
                if( .false. ) then
                    write(ounit, 8827) vvol, dBdmpf(Mvol-1, Mvol-1), Btemn(1, 0, Mvol), dItGpdxtp( 1, 2, Mvol), dItGpdxtp( 1, 1, Mvol), dBdx2(Mvol-1), dBdx2(Mvol)
8827                format("dfp200: vvol = ", i3, ", dBdpsip = ", f10.8, ", dBdpsit = ", f10.8, ", dIpdpsip = ", f10.8, ", dIpdpsit = ", f16.10, ", dBdx2(N-1) = ", f10.8, ", dIpdxj = ", f10.8)
                endif
#endif       
            endif


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
            ! Solve linear system dBdmpf * x = dBdx2. Solution (x) is stored in dBdx2 on exit
            call DGESV(order, 1   , dBdmpf(1:order,1:order), order, IPIV, dBdx2(1:order), order, IDGESV )

            ! Fill dmupfdx array with solution and multiply by coordinate conditioning factor
            dmupfdx(1     , vvol, 2, idof, 1) = zero ! First poloidal flux is always zero
            dmupfdx(2:Mvol, vvol, 2, idof, 1) = lfactor * dBdx2(1:Mvol-1) ! These are the derivatives of pflux

            if( Lfreebound.eq.1 ) then
                dmupfdx(Mvol, vvol, 1, idof, 1) = lfactor * dBdx2(Mvol) ! This is the derivative of tflux
            endif

            ! Free memory
            DALLOCATE( dBdmpf )
            DALLOCATE( dBdx2  )
            DALLOCATE( IPIV   )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    else ! LocalConstraint

        if( Lconstraint.eq.1 ) then
            iflag = -1 ; WCALL( dfp200, tr00ab, ( vvol, mn, lmns, Nt, Nz, iflag, diotadxup(0:1,-1:2,vvol) ) ) ! compute d(transform)/dx;
        endif

        if( Lvacuumregion .and. Lconstraint.ge.0 ) then
            iflag = -1 ; WCALL( dfp200, curent, ( vvol, mn,       Nt, Nz, iflag, dItGpdxtp(0:1,-1:2,vvol) ) ) ! compute d(Itor,Gpol)/dx;
        endif

        dmupfdx(vvol,1,1,idof,innout) = zero        ! Prepare array
        dmupfdx(vvol,1,2,idof,innout) = zero        ! Prepare array

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
        
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
#ifdef DEBUG

    if( Lcheck.eq.4 ) then ! check derivatives of field;

        if( (ncpu.gt.1) .and. (Lconstraint.eq.3) ) then
            goto 8294
        endif

        dBdX%L = .false.

        if( LocalConstraint ) then
            WCALL(dfp200, deallocate_geometry_matrices, (LcomputeDerivatives))
            WCALL(dfp200, deallocate_Beltrami_matrices, (LcomputeDerivatives))
            WCALL(dfp200, intghs_workspace_destroy, ()))
        endif

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

            ! Solve Beltrami equation consistently with constraints
            Xdof(1:Mvol-1) = zero;
            if( LocalConstraint ) then

                SALLOCATE( Fvec, (1:Mvol-1), zero)

                Ndofgl = 0; Fvec(1:Mvol-1) = 0; iflag = 0;
                Xdof(1:Mvol-1) = dpflux(2:Mvol) + xoffset
                
                ! Solve for field
                dBdX%L = .false. ! No need for derivatives in this context
                WCALL(dfp200, dfp100, (Ndofgl, Xdof, Fvec, iflag) )

                DALLOCATE( Fvec )

            ! --------------------------------------------------------------------------------------------------
            ! Global constraint - call the master thread calls hybrd1 on dfp100, others call dfp100_loop.
            else
        
                IPDtdPf = zero
                Xdof(1:Mvol-1)   = dpflux(2:Mvol) + xoffset

                if( Lfreebound.eq.1 ) then
                    ! Mvol-1 surface current and 1 poloidal linking current constraints
                    Ndofgl = Mvol
                else
                    ! Mvol-1 surface current                                constraints
                    Ndofgl = Mvol-1
                endif 

                SALLOCATE(dpfluxout, (1:Ndofgl), zero )
                SALLOCATE(     Fvec, (1:Ndofgl), zero )
                SALLOCATE(     IPIV, (1:Mvol-1), zero )

                WCALL(dfp200, dfp100, (Ndofgl, Xdof(1:Mvol-1), Fvec(1:Ndofgl), 1))

                ! Only one cpu with this test - thus no need for broadcast
                dpfluxout = Fvec
                call DGESV( Ndofgl, 1, IPdtdPf, Ndofgl, IPIV, dpfluxout, Ndofgl, idgesv )

                ! one step Newton's method
                dpflux(2:Mvol)   = dpflux(2:Mvol) - dpfluxout(1:Mvol-1)
                if( Lfreebound.eq.1 ) then
                    dtflux(Mvol) = dtflux(Mvol  ) - dpfluxout(Mvol    )
                endif

                DALLOCATE(IPIV)
                DALLOCATE(Fvec)
                DALLOCATE(dpfluxout)
            endif



            if( LocalConstraint ) then
                LREGION(vvol)
                ll = Lrad(vvol)        ! shorthand
                NN = NAdof(vvol)     ! shorthand;

                if( Lplasmaregion ) then
                    imupf_local(1:2,isymdiff) = (/     mu(vvol), dpflux(vvol) /) ! mu     and dpflux are computed for the perturbed geometry by ma02aa/mp00ac if Lconstraint=1;
                else ! if( Lvacuumregion ) ;
                    imupf_local(1:2,isymdiff) = (/ dtflux(vvol), dpflux(vvol) /) ! dtflux and dpflux are computed for the perturbed geometry by ma02aa/mp00ac if Lconstraint=1;
                endif


            else ! global constraint
                if( Lplasmaregion ) then
                    imupf_global(1:Mvol,1,isymdiff) = (/  mu(1:Mvol)     /) ! mu     is computed for the perturbed geometry by ma02aa/mp00ac
                    imupf_global(1:Mvol,2,isymdiff) = (/  dpflux(1:Mvol) /) ! dpflux is computed for the perturbed geometry by ma02aa/mp00ac
                else ! if( Lvacuumregion ) ;
                    imupf_global(1:Mvol,1,isymdiff) = (/ dtflux(1:Mvol) /) ! dtflux is computed for the perturbed geometry by ma02aa/mp00ac 
                    imupf_global(1:Mvol,2,isymdiff) = (/ dpflux(1:Mvol) /) ! dpflux is computed for the perturbed geometry by ma02aa/mp00ac 
                endif

            endif
        enddo ! end of do isymdiff;

        if( LocalConstraint ) then
            ! reallocate matrices for next iteration
            WCALL(dfp200, intghs_workspace_init, (vvol))
            WCALL(dfp200, allocate_Beltrami_matrices, (vvol,LcomputeDerivatives))
            WCALL(dfp200, allocate_geometry_matrices, (vvol,LcomputeDerivatives))
        endif

8294        continue
        ! Evaluate derivatives using finite differences
        if( LocalConstraint ) then
            imupf_local(1:2,0)      = ( - 1 * imupf_local(1:2, 2)  + 8 * imupf_local(1:2, 1)&
                                        - 8 * imupf_local(1:2,-1)  + 1 * imupf_local(1:2,-2)   ) / ( 12 * dRZ )
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
            if( Lplasmaregion ) then
                write(ounit,3004) cput-cpus, myid, vvol, im(ii), in(ii), irz, issym, 1, "dmu finite-diff", imupf_global(1:Mvol,1,0)
                write(ounit,3004) cput-cpus, myid, vvol, im(ii), in(ii), irz, issym, 1, "dmu analytic   ", dmupfdx(1:Mvol,vvol,1,idof,1) / lfactor
            else ! vacuum
                write(ounit,3004) cput-cpus, myid, vvol, im(ii), in(ii), irz, issym, 1, "dtflux finite-diff", imupf_global(1:Mvol,1,0)
                write(ounit,3004) cput-cpus, myid, vvol, im(ii), in(ii), irz, issym, 1, "dtflux analytic   ", dmupfdx(1:Mvol,vvol,1,idof,1) / lfactor
            endif

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


    endif ! end of if( Lcheck.eq.4 ) ;

#endif
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

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


subroutine evaluate_dBB(lvol, idof, innout, issym, irz, ii, dBB, XX, YY, length, dRR, dZZ, dII, dLL, dPP, Ntz)

! Evaluate the derivative of the square of the magnetic field modulus. Add spectral constraint derivatives if
! required. 

! INPUT
! -----
!    lvol:     Volume number
!   idof:     Labels degree of freedom
!   innout: Index for loop on inner / outer interface (which interface is perturbed)
!   issym:  Index for loop on stellarator non-symmetric terms
!   irz:     Index for loop on R or Z modes
!   ii:     Index for loop on Fourier mode

! MODULES
! -------
  use constants, only : zero, half, one, two
  
  use numerical, only : small
  
  use fileunits, only : ounit
  
  use inputlist, only : Wmacros, Wdfp200, ext, Ntor, Igeometry, Lconstraint, &
                        gamma, adiabatic, pscale, Lcheck, dRZ, Nvol, &
                        epsilon, Lrad
  
  use cputiming, only : Tdfp200
  
  use allglobal, only : ncpu, myid, cpus, &
                        Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, LocalConstraint, &
                        Mvol, dpflux, &
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
                        LGdof, NGdof, NAdof, &
                        vvolume, dvolume, &
                        Rij, Zij, sg, guvij, iRij, iZij, dRij, dZij, tRij, tZij, & ! Jacobian and metrics; computed in coords;
                        dFFdRZ,HdFFdRZ, dBBdmp, &
                        BBweight, & ! exponential weight on force-imbalance harmonics;
                        psifactor, &
                        mn, Iquad, &
                        dRodR, dRodZ, dZodR, dZodZ, dBdX, &
                        xoffset

 LOCALS
!------

INTEGER                 :: iocons, lvol, ideriv, id, iflag, Lcurvature, innout, issym, irz, ii, ifail, idoc, idof, Ntz
REAL                    :: lss, DDl, MMl
REAL                    :: dAt(1:Ntz,-1:2), dAz(1:Ntz,-1:2), XX(1:Ntz), YY(1:Ntz), dBB(1:Ntz,-1:2), dII(1:Ntz), dLL(1:Ntz)
REAL                    :: dPP(1:Ntz), length(1:Ntz), dRR(1:Ntz,-1:2), dZZ(1:Ntz,-1:2), constraint(1:Ntz)
#ifdef DEBUG
INTEGER					:: isymdiff, ll, NN, ndofgl, pvol
REAL                    :: Fvec(1:Mvol-1), Xdof(1:Mvol-1)
REAL,   allocatable 	:: oRbc(:,:), oZbs(:,:), oRbs(:,:), oZbc(:,:) ! original geometry;
REAL,   allocatable 	:: idBB(:,:), iforce(:,:), iposition(:,:)
CHARACTER               :: packorunpack
#endif





do iocons = 0, 1

    if( lvol.eq.   1 .and. iocons.eq.0 ) cycle ! fixed inner boundary (or coordinate axis); no constraints;
    if( lvol.eq.Mvol .and. iocons.eq.1 ) cycle ! fixed outer boundary                     ; no constraints;


! dBBdmp CONSTRUCTION
! ===================

! Evaluate derivatives of B^2 w.r.t mu and pflux
! ----------------------------------------------
    if( Lconstraint.eq.1 .OR. Lconstraint.eq.3 ) then ! first, determine how B^2 varies with mu and dpflux;

        iflag = 1
        do ideriv=1, 2
            !call evaluate_Bsquare(iocons, lvol, dBB, dAt, dAz, XX, YY, length, DDl, MMl, ideriv)! In a subroutine; called somewhere else when semi global constraint
            WCALL(dfp200, lforce, (lvol, iocons, ideriv, Ntz, dBB, XX, YY, length, DDl, MMl, iflag) )
        enddo

        call tfft(    Nt, Nz, dBB(1:Ntz,1), dBB(1:Ntz,2), & ! derivatives of B^2 wrt mu and dpflux;
                    mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), ifail )
        idoc = 0
        dBBdmp(idoc+1:idoc+mn  ,lvol,iocons,1) = efmn(1:mn) * BBweight(1:mn) ! pressure;
        dBBdmp(idoc+1:idoc+mn  ,lvol,iocons,2) = cfmn(1:mn) * BBweight(1:mn) ! pressure;
        idoc = idoc + mn   ! even;

        ! Add spectral constraints; spectral constraints do not depend on mu or dpflux; thus add nothing
        ! Commented - kept for understanding
        !if( Igeometry.ge.3 ) then 
        !    idoc = idoc + mn-1 ! oddd;
        !endif ! end of if( Igeometry.ge.3) ;

        if( NOTstellsym ) then
            dBBdmp(idoc+1:idoc+mn-1,lvol,iocons,1) = ofmn(2:mn) * BBweight(2:mn) ! pressure;
            dBBdmp(idoc+1:idoc+mn-1,lvol,iocons,2) = sfmn(2:mn) * BBweight(2:mn) ! pressure;
            idoc = idoc + mn-1 ! oddd;

            ! Add spectral constraints; spectral constraints do not depend on mu or dpflux; thus add nothing
            ! Commented - kept for understanding
            !if( Igeometry.ge.3 ) then 
            !    idoc = idoc + mn   ! even;
            !endif ! end of if( Igeometry.ge.3) ;

        endif ! end of if( NOTstellsym) ;

    endif ! end of if( Lconstraint.eq.1 .OR. Lconstraint.eq.3 ) ;



! Evaluate B^2 (no derivatives)
! -----------------------------------
    ideriv = 0; iflag=1

    WCALL( dfp200, lforce, (lvol, iocons, ideriv, Ntz, dBB, XX, YY, length, DDl, MMl, iflag) )


! dFFdRZ CONSTRUCTION
! ===================

! B square contribution
! ---------------------
    ideriv = -1; iflag=0

    WCALL( dfp200, lforce, (lvol, iocons, ideriv, Ntz, dBB, XX, YY, length, DDl, MMl, iflag) )

    ! Add derivatives of pressure as well
    FATAL( dfp200, vvolume(lvol).lt.small, shall divide by vvolume(lvol)**(gamma+one) )

    ! Derivatives of force wrt geometry; In real space.
    ijreal(1:Ntz) = - adiabatic(lvol) * pscale * gamma * dvolume / vvolume(lvol)**(gamma+one) + dBB(1:Ntz,-1)



! Spectral condensation contribution
! ----------------------------------

    dLL(1:Ntz) = zero ! either no spectral constraint, or not the appropriate interface;
    dPP(1:Ntz) = zero ! either no spectral constraint, or not the appropriate interface;

    if( Igeometry.ge.3 ) then ! spectral constraints are only required in toroidal or extended-cylindrical geometry;

        if( innout.eq.1 .and. iocons.eq.1 ) then ! include derivatives of spectral constraints;

#ifdef DEBUG
            FATAL( dfp200, abs(DDl).lt.small, divide by zero on spectral constraint )
#endif

            if( issym.eq.0 ) then ! take derivatives wrt Rbc and Zbs;
                if( irz.eq.0 ) then ! take derivative wrt Rbc;
                    dII(1:Ntz) = - im(ii) * sini(1:Ntz,ii) * ( XX(1:Ntz) - MMl * iRij(1:Ntz,lvol) ) &
                    - two * ( mmpp(ii) - MMl ) * iRbc(ii,lvol) * ( Rij(1:Ntz,2,0) * iRij(1:Ntz,lvol) + Zij(1:Ntz,2,0) * iZij(1:Ntz,lvol) ) / DDl &
                    + Rij(1:Ntz,2,0) * ( mmpp(ii) - MMl ) * cosi(1:Ntz,ii)
                else ! take derivative wrt Zbs;
                    dII(1:Ntz) = + im(ii) * cosi(1:Ntz,ii) * ( YY(1:Ntz) - MMl * iZij(1:Ntz,lvol) ) &
                    - two * ( mmpp(ii) - MMl ) * iZbs(ii,lvol) * ( Rij(1:Ntz,2,0) * iRij(1:Ntz,lvol) + Zij(1:Ntz,2,0) * iZij(1:Ntz,lvol) ) / DDl &
                    + Zij(1:Ntz,2,0) * ( mmpp(ii) - MMl ) * sini(1:Ntz,ii)
                endif ! end of if( irz.eq.0 ) ;
            else                  ! take derivatives wrt Rbs and Zbc;
                if( irz.eq.0 ) then 
                    dII(1:Ntz) = + im(ii) * cosi(1:Ntz,ii) * ( XX(1:Ntz) - MMl * iRij(1:Ntz,lvol) ) &
                    - two * ( mmpp(ii) - MMl ) * iRbs(ii,lvol) * ( Rij(1:Ntz,2,0) * iRij(1:Ntz,lvol) + Zij(1:Ntz,2,0) * iZij(1:Ntz,lvol) ) / DDl &
                    + Rij(1:Ntz,2,0) * ( mmpp(ii) - MMl ) * sini(1:Ntz,ii)
                else                
                    dII(1:Ntz) = - im(ii) * sini(1:Ntz,ii) * ( YY(1:Ntz) - MMl * iZij(1:Ntz,lvol) ) &
                    - two * ( mmpp(ii) - MMl ) * iZbc(ii,lvol) * ( Rij(1:Ntz,2,0) * iRij(1:Ntz,lvol) + Zij(1:Ntz,2,0) * iZij(1:Ntz,lvol) ) / DDl &
                    + Zij(1:Ntz,2,0) * ( mmpp(ii) - MMl ) * cosi(1:Ntz,ii)
                endif !matches ( irz.eq.0 )
            endif! matches ( issym.eq.0 )

        else

        dII(1:Ntz) = zero ! either no spectral constraint, or not the appropriate interface;

        endif ! end of if( innout.eq.1 .and. iocons.eq.1 ) ;

        constraint(1:Ntz) = + ( dRij(1:Ntz,lvol) * tRij(1:Ntz,lvol-1+iocons) + dZij(1:Ntz,lvol) * tZij(1:Ntz,lvol-1+iocons) ) / length(1:Ntz)

        if( iocons.eq.0 ) then ! take derivatives of constraints at inner boundary;

            if( innout.eq.0 ) then ! derivative wrt inner boundary coefficient;
                !write(ounit,'("dfp200 : " 10x " : A ; lvol="i3" ; iocons="i2" ; innout="i2" ;")') lvol, iocons, innout
                if( issym.eq.0 ) then ! take derivatives wrt Rbc and Zbs;
                    if( irz.eq.0 ) then ; dLL(1:Ntz) = + ( - cosi(1:Ntz,ii) * tRij(1:Ntz,lvol-1) - dRij(1:Ntz,lvol) * im(ii) * sini(1:Ntz,ii) ) / length(1:Ntz) & 
                                                                        + constraint(1:Ntz) * dRij(1:Ntz,lvol) * cosi(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
                    else                ; dLL(1:Ntz) = + ( - sini(1:Ntz,ii) * tZij(1:Ntz,lvol-1) + dZij(1:Ntz,lvol) * im(ii) * cosi(1:Ntz,ii) ) / length(1:Ntz) & 
                                                                        + constraint(1:Ntz) * dZij(1:Ntz,lvol) * sini(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
                    endif
                else ! if issym.eq.1 ; take derivatives wrt Rbs and Zbc;
                    if( irz.eq.0 ) then ; dLL(1:Ntz) = + ( - sini(1:Ntz,ii) * tRij(1:Ntz,lvol-1) + dRij(1:Ntz,lvol) * im(ii) * cosi(1:Ntz,ii) ) / length(1:Ntz) & 
                                                                        + constraint(1:Ntz) * dRij(1:Ntz,lvol) * sini(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
                    else                ; dLL(1:Ntz) = + ( - cosi(1:Ntz,ii) * tZij(1:Ntz,lvol-1) - dZij(1:Ntz,lvol) * im(ii) * sini(1:Ntz,ii) ) / length(1:Ntz) & 
                                                                        + constraint(1:Ntz) * dZij(1:Ntz,lvol) * cosi(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
                    endif
                endif
            else ! if innout.eq.1 ; derivative wrt outer boundary coefficient;
                !write(ounit,'("dfp200 : " 10x " : B ; lvol="i3" ; iocons="i2" ; innout="i2" ;")') lvol, iocons, innout
                if( issym.eq.0 ) then ! take derivatives wrt Rbc and Zbs;
                    if( irz.eq.0 ) then ; dLL(1:Ntz) = + ( + cosi(1:Ntz,ii) * tRij(1:Ntz,lvol-1)                                              ) / length(1:Ntz) & 
                                                                        - constraint(1:Ntz) * dRij(1:Ntz,lvol) * cosi(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
                    else                ; dLL(1:Ntz) = + ( + sini(1:Ntz,ii) * tZij(1:Ntz,lvol-1)                                              ) / length(1:Ntz) & 
                                                                        - constraint(1:Ntz) * dZij(1:Ntz,lvol) * sini(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
                    endif
                else ! if issym.eq.1 ; take derivatives wrt Rbs and Zbc;
                    if( irz.eq.0 ) then ; dLL(1:Ntz) = + ( + sini(1:Ntz,ii) * tRij(1:Ntz,lvol-1)                                              ) / length(1:Ntz) & 
                                                                        - constraint(1:Ntz) * dRij(1:Ntz,lvol) * sini(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
                    else                ; dLL(1:Ntz) = + ( + cosi(1:Ntz,ii) * tZij(1:Ntz,lvol-1)                                              ) / length(1:Ntz) & 
                                                                        - constraint(1:Ntz) * dZij(1:Ntz,lvol) * cosi(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
                    endif
                endif
            endif

        else ! if iocons.eq.1 ; take derivatives of constraints at outer boundary;

            if( innout.eq.0 ) then ! derivative wrt inner boundary coefficient;
                !write(ounit,'("dfp200 : " 10x " : C ; lvol="i3" ; iocons="i2" ; innout="i2" ;")') lvol, iocons, innout
                if( issym.eq.0 ) then ! take derivatives wrt Rbc and Zbs;
                    if( irz.eq.0 ) then ; dLL(1:Ntz) = + ( - cosi(1:Ntz,ii) * tRij(1:Ntz,lvol  )                                              ) / length(1:Ntz) & 
                                                                        + constraint(1:Ntz) * dRij(1:Ntz,lvol) * cosi(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
                    else                ; dLL(1:Ntz) = + ( - sini(1:Ntz,ii) * tZij(1:Ntz,lvol  )                                              ) / length(1:Ntz) & 
                                                                        + constraint(1:Ntz) * dZij(1:Ntz,lvol) * sini(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
                    endif
                else                  ! take derivatives wrt Rbs and Zbc;
                    if( irz.eq.0 ) then ; dLL(1:Ntz) = + ( - sini(1:Ntz,ii) * tRij(1:Ntz,lvol  )                                              ) / length(1:Ntz) & 
                                                                        + constraint(1:Ntz) * dRij(1:Ntz,lvol) * sini(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
                    else                ; dLL(1:Ntz) = + ( - cosi(1:Ntz,ii) * tZij(1:Ntz,lvol  )                                              ) / length(1:Ntz) & 
                                                                        + constraint(1:Ntz) * dZij(1:Ntz,lvol) * cosi(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
                    endif
                endif
            else ! if innout.eq.1 ; derivative wrt outer boundary coefficient;
                !#ifdef AXIS
                if( Igeometry.eq.3 .and. lvol.eq.1 ) then ! need to accomodate derivatives of coordinate axis;
                    !#else
                    !            if( Igeometry.eq.3 .and. lvol.lt.1 ) then ! need to accomodate derivatives of coordinate axis;
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
                                                            + ( cosi(1:Ntz,ii) - dRodR(1:Ntz,0,ii) ) * tRij(1:Ntz,lvol) - dRij(1:Ntz,lvol) * im(ii) * sini(1:Ntz,ii) &
                                                            + (                - dZodR(1:Ntz,0,ii) ) * tZij(1:Ntz,lvol) &
                                                            - constraint(1:Ntz) &
                                                            * ( dRij(1:Ntz,lvol) * ( cosi(1:Ntz,ii) - dRodR(1:Ntz,0,ii) )   &
                                                            + dZij(1:Ntz,lvol) * (                - dZodR(1:Ntz,0,ii) ) ) / length(1:Ntz) ) / length(1:Ntz)
                        else                ; dLL(1:Ntz) = ( &   ! d/dZbs ;
                                                            + (                - dRodZ(1:Ntz,1,ii) ) * tRij(1:Ntz,lvol) + dZij(1:Ntz,lvol) * im(ii) * cosi(1:Ntz,ii) &
                                                            + ( sini(1:Ntz,ii) - dZodZ(1:Ntz,1,ii) ) * tZij(1:Ntz,lvol) &
                                                            - constraint(1:Ntz) &
                                                            * ( dRij(1:Ntz,lvol) * (                - dRodZ(1:Ntz,1,ii) )   &
                                                            + dZij(1:Ntz,lvol) * ( sini(1:Ntz,ii) - dZodZ(1:Ntz,1,ii) ) ) / length(1:Ntz) ) / length(1:Ntz)
                        endif ! end of if( irz.eq.0 ) ;
                    else
                        if( irz.eq.0 ) then ; dLL(1:Ntz) =    ( &   ! d/dRbs ;
                                                            + ( sini(1:Ntz,ii) - dRodR(1:Ntz,1,ii) ) * tRij(1:Ntz,lvol) + dRij(1:Ntz,lvol) * im(ii) * cosi(1:Ntz,ii)     &
                                                            + (                - dZodR(1:Ntz,1,ii) ) * tZij(1:Ntz,lvol)                                                 &
                                                            - constraint(1:Ntz)                                                                                         &
                                                            * (   dRij(1:Ntz,lvol) * ( sini(1:Ntz,ii)   - dRodR(1:Ntz,1,ii) )                                              &
                                                                + dZij(1:Ntz,lvol) * (                  - dZodR(1:Ntz,1,ii) ) ) / length(1:Ntz) ) / length(1:Ntz)
                        else                ; dLL(1:Ntz) =     ( &   ! d/dZbs ;
                                                            + (                - dRodZ(1:Ntz,0,ii) ) * tRij(1:Ntz,lvol) - dZij(1:Ntz,lvol) * im(ii) * sini(1:Ntz,ii) &
                                                            + ( cosi(1:Ntz,ii) - dZodZ(1:Ntz,0,ii) ) * tZij(1:Ntz,lvol) &
                                                            - constraint(1:Ntz) &
                                                            * ( dRij(1:Ntz,lvol) * (                - dRodZ(1:Ntz,0,ii) )   &
                                                            + dZij(1:Ntz,lvol) * ( cosi(1:Ntz,ii) - dZodZ(1:Ntz,0,ii) ) ) / length(1:Ntz) ) / length(1:Ntz)
                        endif ! end of if( irz.eq.0 ) ;
                    endif
                else
                    if( issym.eq.0 ) then ! take derivatives wrt Rbc and Zbs;
                        if( irz.eq.0 ) then ; dLL(1:Ntz) = + ( + cosi(1:Ntz,ii) * tRij(1:Ntz,lvol  ) - dRij(1:Ntz,lvol) * im(ii) * sini(1:Ntz,ii) ) / length(1:Ntz) & 
                                                                            - constraint(1:Ntz) * dRij(1:Ntz,lvol) * cosi(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
                        else                ; dLL(1:Ntz) = + ( + sini(1:Ntz,ii) * tZij(1:Ntz,lvol  ) + dZij(1:Ntz,lvol) * im(ii) * cosi(1:Ntz,ii) ) / length(1:Ntz) & 
                                                                            - constraint(1:Ntz) * dZij(1:Ntz,lvol) * sini(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
                        endif
                    else                  ! take derivatives wrt Rbs and Zbc;
                        if( irz.eq.0 ) then ; dLL(1:Ntz) = + ( + sini(1:Ntz,ii) * tRij(1:Ntz,lvol  ) + dRij(1:Ntz,lvol) * im(ii) * cosi(1:Ntz,ii) ) / length(1:Ntz) & 
                                                                            - constraint(1:Ntz) * dRij(1:Ntz,lvol) * sini(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
                        else                ; dLL(1:Ntz) = + ( + cosi(1:Ntz,ii) * tZij(1:Ntz,lvol  ) - dZij(1:Ntz,lvol) * im(ii) * sini(1:Ntz,ii) ) / length(1:Ntz) & 
                                                                            - constraint(1:Ntz) * dZij(1:Ntz,lvol) * cosi(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
                        endif ! end of if( irz.eq.0 ) ;
                    endif ! end of if( issym.eq.0 ) ;
                endif  ! end of if( Igeometry.eq.3 .and. lvol.eq.1 ) ;
            endif ! end of if( innout.eq.0 ) ;

        endif ! end of if( iocons.eq.0 ) ;

    endif ! end of if( Igeometry.ge.3 ) ;

    ! Map to Fourier space
    call tfft(  Nt, Nz, ijreal(1:Ntz), dII(1:Ntz), & ! recall that ijreal contains derivatives of pressure term;
                mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), ifail )


    call tfft(  Nt, Nz, dPP(1:Ntz)   , dLL(1:Ntz), &
                mn, im(1:mn), in(1:mn), evmn(1:mn), odmn(1:mn), comn(1:mn), simn(1:mn), ifail ) ! evmn and odmn are available as workspace;


    FATAL( dfp200, lvol-1+innout.gt.Mvol, psifactor needs attention )

    
    idoc = 0
  
    ! Plasma and magnetic pressure;
    ;   dFFdRZ(idoc+1:idoc+mn    ,iocons,idof,innout,lvol) = + efmn(1:mn) * psifactor(ii,lvol-1+innout) * BBweight(1:mn)
    
    idoc = idoc + mn   ! even;
    if( Igeometry.ge.3 ) then ! Add spectral constraints;
        dFFdRZ(idoc+1:idoc+mn-1  ,iocons,idof,innout,lvol) = - sfmn(2:mn) * psifactor(ii,lvol-1+innout) * epsilon       & ! spectral condensation;
                                                             - simn(2:mn) * psifactor(ii,lvol-1+innout) * sweight(lvol)   ! poloidal length constraint;
    ! if( Ntor.gt.0 ) then
    !  dFFdRZ(idoc+1:idoc+Ntor  ,iocons,idof,innout,lvol) = + odmn(2:Ntor+1) * psifactor(ii,lvol-1+innout) * apsilon
    ! endif
      idoc = idoc + mn-1 ! odd;
    endif ! end of if( Igeometry.ge.3) ;

    if( NOTstellsym ) then ! Construct non-stellarator symmetric terms

    ! Plasma and magnetic pressure;
    ;       dFFdRZ(idoc+1:idoc+mn-1  ,iocons,idof,innout,lvol) = + ofmn(2:mn) * psifactor(ii,lvol-1+innout) * BBweight(2:mn) 
      
        idoc = idoc + mn-1 ! odd;
        if( Igeometry.ge.3 ) then ! Add spectral constraints;
            dFFdRZ(idoc+1:idoc+mn    ,iocons,idof,innout,lvol) = - cfmn(1:mn) * psifactor(ii,lvol-1+innout) * epsilon       & ! spectral condensation;
                                                                 - comn(1:mn) * psifactor(ii,lvol-1+innout) * sweight(lvol)   ! poloidal length constraint;
            !if( Ntor.ge.0 ) then
                ! dFFdRZ(idoc+1:idoc+Ntor+1,iocons,idof,innout,lvol) = + evmn(1:Ntor+1) * psifactor(ii,lvol-1+innout) * apsilon ! poloidal origin      ;
            !endif
            idoc = idoc + mn   ! even;
        endif ! end of if( Igeometry.ge.3) ;

    endif ! end of if( NOTstellsym) ;

#ifdef DEBUG
    FATAL( dfp200, idoc.ne.LGdof, counting error )
#endif
     
enddo ! end of do iocons;

end subroutine evaluate_dBB


subroutine hessian_dFFdRZ(lvol, idof, innout, issym, irz, ii, dBB, XX, YY, length, dRR, dZZ, dII, dLL, dPP, Ntz)

! Evaluate the derivative of the square of the magnetic field modulus. Add spectral constraint derivatives if
! required. 

! INPUT
! -----
!    lvol:     Volume number
!   idof:     Labels degree of freedom
!   innout: Index for loop on inner / outer interface (which interface is perturbed)
!   issym:  Index for loop on stellarator non-symmetric terms
!   irz:     Index for loop on R or Z modes
!   ii:     Index for loop on Fourier mode

! MODULES
! -------
  use constants, only : zero, half, one, two
  
  use numerical, only : small
  
  use fileunits, only : ounit
  
  use inputlist, only : Wmacros, Wdfp200, ext, Ntor, Igeometry, Lconstraint, &
                        gamma, adiabatic, pscale, Lcheck, dRZ, Nvol, &
                        epsilon, Lrad
  
  use cputiming, only : Tdfp200
  
  use allglobal, only : ncpu, myid, cpus, &
                        Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, LocalConstraint, &
                        Mvol, dpflux, &
                        iRbc, iZbs, iRbs, iZbc, & ! Fourier harmonics of geometry; vector of independent variables, position, is "unpacked" into iRbc,iZbs;
                        NOTstellsym, &
                        mn, im, in, mns, &
                        ijreal, &
                        !hijreal, hijmn,&
                        efmn, ofmn, cfmn, sfmn, &
                        evmn, odmn, comn, simn, &
                        Nt, Nz, &
                        cosi, sini, & ! FFT workspace
                        sweight, &
                        mmpp, &
                        LGdof, NGdof, NAdof, &
                        vvolume, dvolume, &
                        Rij, Zij, sg, guvij, iRij, iZij, dRij, dZij, tRij, tZij, & ! Jacobian and metrics; computed in coords;
                        dFFdRZ, dBBdmp, HdFFdRZ, &
                        BBweight, & ! exponential weight on force-imbalance harmonics;
                        psifactor, &
                        mn, Iquad, &
                        dRodR, dRodZ, dZodR, dZodZ, dBdX, &
                        xoffset

 LOCALS
!------

INTEGER                 :: iocons, lvol, ideriv, id, iflag, Lcurvature, innout, issym, irz, ii, ifail, idoc, idof, Ntz
REAL                    :: lss, DDl, MMl
REAL                    :: dAt(1:Ntz,-1:2), dAz(1:Ntz,-1:2), XX(1:Ntz), YY(1:Ntz), dBB(1:Ntz,-1:2), dII(1:Ntz), dLL(1:Ntz)
REAL                    :: dPP(1:Ntz), length(1:Ntz), dRR(1:Ntz,-1:2), dZZ(1:Ntz,-1:2), constraint(1:Ntz)
#ifdef DEBUG
INTEGER					:: isymdiff, ll, NN, ndofgl, pvol
REAL                    :: Fvec(1:Mvol-1), Xdof(1:Mvol-1)
REAL,   allocatable 	:: oRbc(:,:), oZbs(:,:), oRbs(:,:), oZbc(:,:) ! original geometry;
REAL,   allocatable 	:: idBB(:,:), iforce(:,:), iposition(:,:)
CHARACTER               :: packorunpack
#endif



do iocons = 0, 1

    if( lvol.eq.   1 .and. iocons.eq.0 ) cycle ! fixed inner boundary (or coordinate axis); no constraints;
    if( lvol.eq.Mvol .and. iocons.eq.1 ) cycle ! fixed outer boundary                     ; no constraints;


! Evaluate B^2 (no derivatives)
! -----------------------------------
    ideriv = 0; iflag=1

    WCALL( dfp200, lforce, (lvol, iocons, ideriv, Ntz, dBB, XX, YY, length, DDl, MMl, iflag) )

! hessian_dFFdRZ CONSTRUCTION
! ===================

! B square contribution
! ---------------------
    ideriv = -1; iflag=0

    WCALL( dfp200, lforce, (lvol, iocons, ideriv, Ntz, dBB, XX, YY, length, DDl, MMl, iflag) )

    ! Add derivatives of pressure as well
    FATAL( dfp200, vvolume(lvol).lt.small, shall divide by vvolume(lvol)**(gamma+one) )

    ! Derivatives of force wrt geometry; In real space.
    ijreal(1:Ntz) = - adiabatic(lvol) * pscale * gamma * dvolume / vvolume(lvol)**(gamma+one) + dBB(1:Ntz,-1)*Rij(1:Ntz,0,0) 

! Spectral condensation contribution
! ----------------------------------

    dLL(1:Ntz) = zero ! either no spectral constraint, or not the appropriate interface;
    dPP(1:Ntz) = zero ! either no spectral constraint, or not the appropriate interface;
    dII(1:Ntz) = zero

   
    ! Map to Fourier space
    call tfft(  Nt, Nz, ijreal(1:Ntz), dII(1:Ntz), & ! recall that ijreal contains derivatives of pressure term;
                mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), ifail )


    FATAL( dfp200, lvol-1+innout.gt.Mvol, psifactor needs attention )
    
    idoc = 0
  
    ! Plasma and magnetic pressure;
    ;   HdFFdRZ(idoc+1:idoc+mn    ,iocons,idof,innout,lvol) = + efmn(1:mn) !* psifactor(ii,lvol-1+innout) * BBweight(1:mn)
    idoc = idoc + mn   ! even;
    !if( Igeometry.ge.3 ) then ! Add spectral constraints;
     !   HdFFdRZ(idoc+1:idoc+mn-1  ,iocons,idof,innout,lvol) = - sfmn(2:mn) * psifactor(ii,lvol-1+innout) * epsilon       & ! spectral condensation;
      !                                                       - simn(2:mn) * psifactor(ii,lvol-1+innout) * sweight(lvol)   ! poloidal length constraint;
    ! if( Ntor.gt.0 ) then
    !  dFFdRZ(idoc+1:idoc+Ntor  ,iocons,idof,innout,lvol) = + odmn(2:Ntor+1) * psifactor(ii,lvol-1+innout) * apsilon
    ! endif
     ! idoc = idoc + mn-1 ! odd;
    !endif ! end of if( Igeometry.ge.3) ;

    if( NOTstellsym ) then ! Construct non-stellarator symmetric terms

    ! Plasma and magnetic pressure;
    ;       HdFFdRZ(idoc+1:idoc+mn-1  ,iocons,idof,innout,lvol) = + ofmn(2:mn) !* psifactor(ii,lvol-1+innout) * BBweight(2:mn) 
      
        idoc = idoc + mn-1 ! odd;
        ! if( Igeometry.ge.3 ) then ! Add spectral constraints;
        !     HdFFdRZ(idoc+1:idoc+mn    ,iocons,idof,innout,lvol) = - cfmn(1:mn) * psifactor(ii,lvol-1+innout) * epsilon       & ! spectral condensation;
        !                                                          - comn(1:mn) * psifactor(ii,lvol-1+innout) * sweight(lvol)   ! poloidal length constraint;
        !     !if( Ntor.ge.0 ) then
        !         ! dFFdRZ(idoc+1:idoc+Ntor+1,iocons,idof,innout,lvol) = + evmn(1:Ntor+1) * psifactor(ii,lvol-1+innout) * apsilon ! poloidal origin      ;
        !     !endif
        !     idoc = idoc + mn   ! even;
        ! endif ! end of if( Igeometry.ge.3) ;

    endif ! end of if( NOTstellsym) ;

#ifdef DEBUG
    FATAL( dfp200, idoc.ne.LGdof, counting error )
#endif
     
enddo ! end of do iocons;

end subroutine hessian_dFFdRZ

subroutine hessian3D_dFFdRZ(lvol, idof, innout, issym, irz, ii, dBB, XX, YY, length, dRR, dZZ, dII, dLL, dPP, Ntz)

! Evaluate the derivative of the square of the magnetic field modulus. Add spectral constraint derivatives if
! required. 

! INPUT
! -----
!    lvol:     Volume number
!   idof:     Labels degree of freedom
!   innout: Index for loop on inner / outer interface (which interface is perturbed)
!   issym:  Index for loop on stellarator non-symmetric terms
!   irz:     Index for loop on R or Z modes
!   ii:     Index for loop on Fourier mode

! MODULES
! -------
  use constants, only : zero, half, one, two
  
  use numerical, only : small
  
  use fileunits, only : ounit
  
  use inputlist, only : Wmacros, Wdfp200, ext, Ntor, Igeometry, Lconstraint, &
                        gamma, adiabatic, pscale, Lcheck, dRZ, Nvol, &
                        epsilon, Lrad
  
  use cputiming, only : Tdfp200
  
  use allglobal, only : ncpu, myid, cpus, YESstellsym, NOTstellsym,&
                        Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, LocalConstraint, &
                        Mvol, dpflux, &
                        iRbc, iZbs, iRbs, iZbc, & ! Fourier harmonics of geometry; vector of independent variables, position, is "unpacked" into iRbc,iZbs;
                        NOTstellsym, &
                        mn, im, in, mns, &
                        ijreal, &
                        fijreal, &
                        efmn, ofmn, cfmn, sfmn, &
                        evmn, odmn, comn, simn, &
                        Nt, Nz, &
                        cosi, sini, & ! FFT workspace
                        sweight, Lhessian3Dallocated, &
                        mmpp, &
                        LGdof, NGdof, NAdof, &
                        vvolume, dvolume, &
                        Rij, Zij, sg, guvij, iRij, iZij, dRij, dZij, tRij, tZij, & ! Jacobian and metrics; computed in coords;
                        dFFdRZ, dBBdmp, HdFFdRZ, &
                        denergydrr, denergydrz,denergydzr,denergydzz,&
!                        efcol1mn,efcol2mn,efcol3mn,efcol4mn,ofcol1mn,ofcol2mn,ofcol3mn,ofcol4mn, &
                        BBweight, & ! exponential weight on force-imbalance harmonics;
                        psifactor, &
                        mn, Iquad, &
                        dRodR, dRodZ, dZodR, dZodZ, dBdX, &
                        xoffset

 LOCALS
!------

INTEGER                 :: iocons, lvol, ideriv, id, iflag, Lcurvature, innout, issym, irz, ii, ifail, idoc, idof, Ntz
REAL                    :: lss, DDl, MMl
REAL                    :: dAt(1:Ntz,-1:2), dAz(1:Ntz,-1:2), XX(1:Ntz), YY(1:Ntz), dBB(1:Ntz,-1:2), dII(1:Ntz), dLL(1:Ntz)
REAL                    :: dPP(1:Ntz), length(1:Ntz), dRR(1:Ntz,-1:2), dZZ(1:Ntz,-1:2), constraint(1:Ntz)
REAL                    :: ddFcol1(1:Ntz), ddFcol2(1:Ntz),ddFcol3(1:Ntz), ddFcol4(1:Ntz)

#ifdef DEBUG
INTEGER					:: isymdiff, ll, NN, ndofgl, pvol
REAL                    :: Fvec(1:Mvol-1), Xdof(1:Mvol-1)
REAL,   allocatable 	:: oRbc(:,:), oZbs(:,:), oRbs(:,:), oZbc(:,:) ! original geometry;
REAL,   allocatable 	:: idBB(:,:), iforce(:,:), iposition(:,:)
CHARACTER               :: packorunpack
#endif



do iocons = 0, 1

    if( lvol.eq.   1 .and. iocons.eq.0 ) cycle ! fixed inner boundary (or coordinate axis); no constraints;
    if( lvol.eq.Mvol .and. iocons.eq.1 ) cycle ! fixed outer boundary                     ; no constraints;


! Evaluate B^2 (no derivatives)
! -----------------------------------
    ideriv = 0; iflag=0

    WCALL( dfp200, lforce, (lvol, iocons, ideriv, Ntz, dBB, XX, YY, length, DDl, MMl,iflag) )


! hessian_dFFdRZ CONSTRUCTION
! ===================

! B square contribution
! ---------------------
    ideriv = -1; iflag=1

    WCALL( dfp200, lforce, (lvol, iocons, ideriv, Ntz, dBB, XX, YY, length, DDl, MMl, iflag) )

    ! Add derivatives of pressure as well
    FATAL( dfp200, vvolume(lvol).lt.small, shall divide by vvolume(lvol)**(gamma+one) )

    ! Derivatives of force wrt geometry; In real space.
    ijreal(1:Ntz) = - adiabatic(lvol) * pscale * gamma * dvolume / vvolume(lvol)**(gamma+one) + dBB(1:Ntz,-1)

        if(YESstellsym ) then
       !if( issym.eq.0 ) then

            if( irz.eq.0 ) then !derivative wrt Rbs of d/dr(dF)

                 if(innout.eq.iocons) then
                    ddFcol1(1:Ntz) = ijreal(1:Ntz) * Zij(1:Ntz,2,0) * Rij(1:Ntz,0,0) & !Frr       d/dr(dw/dr)
                                   +  (dBB(1:Ntz,0)+adiabatic(lvol) * pscale/vvolume(lvol)**gamma ) * Zij(1:Ntz,2,0) * cosi(1:Ntz,ii)

                    ddFcol2(1:Ntz) = ijreal(1:Ntz) * Rij(1:Ntz,0,0) * Rij(1:Ntz,2,0) &  !Fzr        d/dr(dw/dz)
                                    + (dBB(1:Ntz,0) +adiabatic(lvol) * pscale/vvolume(lvol)**gamma )* Rij(1:Ntz,2,0) * cosi(1:Ntz,ii) &
                                    + (dBB(1:Ntz,0) +adiabatic(lvol) * pscale/vvolume(lvol)**gamma )* Rij(1:Ntz,0,0) * sini(1:Ntz,ii) * (-im(ii))
                  else

                  ddFcol1(1:Ntz) = ijreal(1:Ntz) * Zij(1:Ntz,2,0) * Rij(1:Ntz,0,0) !Frr       d/dr(dw/dr)
                  ddFcol2(1:Ntz) =  ijreal(1:Ntz) * Rij(1:Ntz,0,0) * Rij(1:Ntz,2,0) !Fzr        d/dr(dw/dz)
                 endif
            else !derivative wrt to Zbs of d/dz(dF)
                
                if(innout.eq.iocons) then
                   ddFcol3(1:Ntz) = ijreal(1:Ntz) * Rij(1:Ntz,0,0) * Zij(1:Ntz,2,0) &   !Frz      d/dz(dw/dr)
                                  + (dBB(1:Ntz,0) +adiabatic(lvol) * pscale/vvolume(lvol)**gamma )* Rij(1:Ntz,0,0) * im(ii) *cosi(1:Ntz,ii)
                   ddFcol4(1:Ntz) =  ijreal(1:Ntz) * Rij(1:Ntz,0,0) * Rij(1:Ntz,2,0)  !Fzz    d/dz(dw/dz)
                else
                   ddFcol3(1:Ntz) = ijreal(1:Ntz) * Rij(1:Ntz,0,0) * Zij(1:Ntz,2,0)  !Frz      d/dz(dw/dr)
                    ddFcol4(1:Ntz) =  ijreal(1:Ntz) * Rij(1:Ntz,0,0) * Rij(1:Ntz,2,0)  !Fzz    d/dz(dw/dz)
                endif
           endif 
        else
              !   if( irz.eq.0 ) then
              !       ddFcol1(1:Ntz) = -ijreal(1:Ntz) * Zij(1:Ntz,2,0) * Rij(1:Ntz,0,0) & !Frr       d/dr(dw/dr)
             !                       -(adiabatic(vvol) * pscale/vvolume(vvol)**gamma + dBB(1:Ntz,0)) * Zij(1:Ntz,2,0) * cosi(1:Ntz,ii)

              !       ddFcol2(1:Ntz) =ijreal(1:Ntz) * Rij(1:Ntz,0,0) * Rij(1:Ntz,2,0) &  !Fzr        d/dr(dw/dz)
              !                    + (adiabatic(vvol) * pscale/vvolume(vvol)**gamma + dBB(1:Ntz,0)) * Rij(1:Ntz,2,0) * cosi(1:Ntz,ii) &
             !                    + (adiabatic(vvol) * pscale/vvolume(vvol)**gamma + dBB(1:Ntz,0)) * Rij(1:Ntz,0,0) * sini(1:Ntz,ii) * (-im(ii))
             !   else
             !       ddFcol3(1:Ntz) = -ijreal(1:Ntz) * Rij(1:Ntz,0,0) * Zij(1:Ntz,2,0) &   !Frz      d/dz(dw/dr)
              !                        -(adiabatic(vvol) * pscale/vvolume(vvol)**gamma + dBB(1:Ntz,0)) * Rij(1:Ntz,0,0) * im(ii) *cosi(1:Ntz,ii)

              !       ddFcol4(1:Ntz) = ijreal(1:Ntz) * Rij(1:Ntz,0,0) * Rij(1:Ntz,2,0)      !Fzz    d/dz(dw/dz)
             !   endif
                 ddFcol1(1:Ntz) = zero
                 ddFcol2(1:Ntz) = zero
                 ddFcol3(1:Ntz) = zero
                 ddFcol4(1:Ntz) = zero
              !FATAL(dfp200, .true. work progress for hessian axisymmetric )
        endif

        !Map to Fourier space

                   call tfft(  Nt, Nz, ddFcol1(1:Ntz),ddFcol2(1:Ntz), & ! recall that ijreal contains derivatives of pressure term;
                               mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), ifail )

                   call tfft(  Nt, Nz, ddFcol3(1:Ntz),ddFcol4(1:Ntz), &
                               mn, im(1:mn), in(1:mn), evmn(1:mn), odmn(1:mn), comn(1:mn), simn(1:mn), ifail ) ! evmn and odmn are available as workspace;
        

          !call tfft(Nt , Nz, ddFcol1(1:Ntz),ddFcol2(1:Ntz), &
          !                                 mn, im(1:mn), in(1:mn),  efcol1mn(1:mn), ofcol1mn(1:mn), ofcol2mn(1:mn), efcol2mn(1:mn), ifail )
          !call tfft(Nt , Nz, ddFcol3(1:Ntz),ddFcol4(1:Ntz), &
          !                                  mn, im(1:mn), in(1:mn),  efcol3mn(1:mn), ofcol3mn(1:mn), ofcol4mn(1:mn), efcol4mn(1:mn), ifail )
             !                          hessian3d= [1 2;3 4]

           if (irz.eq.0 .and. ii.eq.1) then
               ;idoc=0
                      ;denergydrr(idoc+1:idoc+mn ,lvol,iocons,idof,innout) = -efmn(1:mn)*half !1 wrr
               ;idoc=idoc+mn
                      ;denergydrr(idoc+1:idoc+mn-1,lvol,iocons,idof,innout) = -sfmn(2:mn)*half ! 2 wrz
               ;idoc=idoc+mn-1
           endif
          if (irz.eq.0 .and. ii .gt.1) then
              ;idoc=0
                     ;denergydrr(idoc+1:idoc+mn ,lvol,iocons,idof,innout) = -efmn(1:mn) !1 wrr
              ;idoc=idoc+mn
                     ;denergydrr(idoc+1:idoc+mn-1,lvol,iocons,idof,innout) = -sfmn(2:mn) ! 2 wrz
              ;idoc=idoc+mn-1
            endif

          if(irz.eq.1 .and. ii .gt.1) then
              ;idoc=0
                     ;denergydzr(idoc+1:idoc+mn-1 ,lvol,iocons,idof,innout) =  evmn(1:mn) !wzr
                        !write(ounit,*) im(ii),in(ii), evmn(1:mn) !vvol, im(ii), in(ii), irz, issym, tdofr, tdofz 

              ;idoc=idoc+mn
                      ;denergydzr(idoc+1:idoc+mn-1 ,lvol,iocons,idof,innout) =  simn(2:mn) !wzz
              ;idoc=idoc+mn-1
          end if


          !if (NOTstellsym) then
           !  FATAL(dfp200, .true. work progress for hessian axisymmetric )
          !endif

            !write(ounit,'("hesian3D : ",f10.2," : efcol1="f10.2")') cput-cpus, efcol1mn(1:mn)
            !print * , Mvol, efcol1mn(1:mn)
            !write(ounit,1000) 'values are:' Mvol, efcol1mn(1:mn)
             !write(90,1000) efcol1mn(1:mn)
             !1000 format(" "10x" "es23.15" ")
             !open(nm1unit, file="."//trim(ext)//".GF.hcol1", status="unknown", form="unformatted")
             !write(nm1unit) NGdof, Mvol
             !write(nm1unit) efcol1mn(1:Ntz)
             !close(nm1unit)

             !open(nm2unit, file="."//trim(ext)//".GF.hcol2", status="unknown", form="unformatted")
             !!write(nm2unit) NGdof, Mvol
             !write(nm2unit) efcol2mn(1:Ntz)
             !close(nm2unit)
             !write(15,1020) denergydrr(idoc+1:idoc+mn ,vvol,iocons,idof,innout)
             !1020 format(""10x" " es23.15" ")

! Spectral condensation contribution
! ----------------------------------
    dLL(1:Ntz) = zero !  no spectral constraint, or not the appropriate interface;
    dPP(1:Ntz) = zero !  no spectral constraint, or not the appropriate interface;
    dII(1:Ntz) = zero !  no angle/spectral width constraint

#ifdef DEBUG
    FATAL( dfp200, idoc.ne.LGdof, counting error )
#endif
     
enddo ! end of do iocons;

end subroutine hessian3D_dFFdRZ 
