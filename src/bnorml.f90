!> \defgroup grp_free-boundary Free-Boundary Computation
!>
!> \file
!> \brief Computes \f${\bf B}_{Plasma} \cdot {\bf e}_\theta \times {\bf e}_\zeta \;\f$ on the computational boundary, \f$\partial {\cal D}\f$.

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
!> \ingroup grp_free-boundary
!>
!> **free-boundary constraint**
!> <ul>
!> <li> The normal field at the computational boundary, \f$\partial {\cal D}\f$, should be equal to
!>      \f$\left({\bf B}_P + {\bf B}_C\right)\cdot {\bf e}_\theta \times {\bf e}_\zeta\f$,
!>      where \f${\bf B}_P\f$ is the "plasma" field (produced by internal plasma currents) and is computed using virtual casing,
!>      and \f${\bf B}_C\f$ is the "vacuum" field (produced by the external coils) and is given on input. </li>
!> <li> The plasma field, \f${\bf B}_P\f$, can only be computed after the equilibrium is determined,
!>      but this information is required to compute the equilibrium to begin with; and so there is an iteration involved. </li>
!> <li> Suggested values of the vacuum field can be self generated; see xspech() for more documentation on this. </li>
!> </ul>
!>
!> **compute the normal field on a regular grid on the computational boundary**
!> <ul>
!> <li> For each point on the compuational boundary, casing() is called to compute the normal field produced by the plasma currents. </li>
!> <li> \todo There is a very clumsy attempt to parallelize this which could be greatly improved.
!>
!> </li>
!> <li> An FFT gives the required Fourier harmonics. </li>
!> </ul>
!> \see casing.f90
!>
!> @param[in] mn total number of Fourier harmonics
!> @param[in] Ntz total number of grid points in \f$\theta\f$ and \f$zeta\f$
!> @param[out] efmn even Fourier coefficients
!> @param[out] ofmn odd Fouier coefficients
subroutine bnorml( mn, Ntz, efmn, ofmn )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  use constants, only : zero, half, one, two, pi, pi2, ten

  use numerical, only : small, vsmall

  use fileunits, only : ounit, lunit

  use inputlist, only : Wmacros, Wbnorml, Igeometry, Lcheck, vcasingtol, vcasingits, vcasingper, Lrad, Lvcgrid  

  use cputiming, only : Tbnorml

  use inputlist, only : vcNz, vcNt

  use allglobal, only : ncpu, myid, cpus, MPI_COMM_SPEC, pi2nfp, Mvol, &
                        Nt, Nz, &
                        Rij, Zij, guvij, sg, TT, &
                        NOTstellsym, Lcoordinatesingularity, &
                        im, in, Ate, Aze, Ato, Azo, &
                        Nt, Nz, cfmn, sfmn, &
                        ijreal, ijimag, jireal, jiimag, &
                        globaljk, virtualcasingfactor, gteta, gzeta, Dxyz, Nxyz, &
                        Jxyz, Pbxyz, YESstellsym

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS

  INTEGER, intent(in)  :: mn, Ntz
  REAL   , intent(out) :: efmn(1:mn), ofmn(1:mn)

  INTEGER              :: lvol, Lparallel, ii, jj, kk, jk, ll, kkmodnp, jkmodnp, ifail, id01daf, nvccalls, icasing
  REAL                 :: zeta, teta, gBn
  REAL                 :: Bxyz(1:Ntz,1:3), distance(1:Ntz)
  REAL                 :: absvcerr, relvcerr, resulth, resulth2, resulth4, deltah4h2, deltah2h
  INTEGER              :: vcstride, Nzwithsym 


  BEGIN(bnorml)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  Lparallel = 1 ! controls choice of parallelization; see below;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  ijreal(1:Ntz) = zero ! normal plasma field; 15 Oct 12;
  ijimag(1:Ntz) = zero 
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef DEBUG
  if( myid.eq.0 .and. Lcheck.eq.6 ) then
   write(ounit,'("bnorml : " 10x " : writing input for xdiagno ; screen comparison ; Ntz =",i7," ;")') Ntz
   open(lunit, file="btest.diagno", status="unknown" )
   write(lunit,'(i9)') Ntz
  endif
#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  ! In stellerator symmetric geometry, the virtual casing field also exhibits symmetriy which we can exploit
  if (YESstellsym) then
    ! we have to compute approximately half points due to symmetry, but including the middle and the edge points! (Closed interval [0, pi])
    Nzwithsym = min(Nz, (Nz + 3) / 2 ) 
  else
    Nzwithsym = Nz
  endif

#ifdef COMPARECASING
! When comparing results, both methods should always run
if (.true.) then
#else
if ( Lvcgrid.eq.1 ) then
#endif
  ! Precompute Jxyz(:,1:3) and the corresponding positions on the high resolution plasma boundary

  !$OMP PARALLEL DO SHARED(Pbxyz, Jxyz) PRIVATE(jk, jj, kk, teta, zeta) COLLAPSE(2)
  do kk = 0, vcNz-1 ; 
    do jj = 0, vcNt-1 ; 
      zeta = kk * pi2 / vcNz
      teta = jj * pi2  / vcNt ; 
      jk = 1 + jj + kk*vcNt

      ! Each MPI rank only computes every a 1/ncpu surfacecurrent() calls 
      select case( Lparallel ) 
      case( 0 ) ! Lparallel = 0 
      if( myid.ne.modulo(kk,ncpu) ) cycle
      case( 1 ) ! Lparallel = 1 
      if( myid.ne.modulo(jk-1,ncpu) ) cycle
      case default ! Lparallel
        FATAL( bnorml, .true., invalid Lparallel in parallelization loop )
      end select 

      ! pbxyz and jxyz are  both [out] parameters
      call surfacecurrent( teta, zeta, Pbxyz(jk,1:3), Jxyz(jk,1:3)  )
    enddo
  enddo

  ! MPI reductions for positions and currents to accumulate them on all ranks (valid because initialized to zero) 
  ! and Broadcast the total currents and evaluation points back to all ranks
  call MPI_Allreduce(MPI_IN_PLACE, Pbxyz(:,1:3), 3*vcNt*vcNz, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_SPEC, ierr )
  if (ierr.ne.MPI_SUCCESS) then
    FATAL( bnorml, .true., error in MPI_Allreduce for Pbxyz )
  endif
  call MPI_Allreduce(MPI_IN_PLACE, Jxyz(:,1:3),  3*vcNt*vcNz, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_SPEC, ierr )
  if (ierr.ne.MPI_SUCCESS) then
    FATAL( bnorml, .true., error in MPI_Allreduce for Jxyz )
  endif

  ! iterate over resolutions of the virtual casing grid to get an estimate of the accuracy. Write the result into ijimag
  deltah2h = 1e12
  !> The surface currents get precomputed at all points, but maybe we don't have to integrate over the whole grid. 
  !> Start with a large stride over the plasmaboundary (= how many points to skip at each step), then get progressively finer.
  !> Do at least vcasingits iterations (sqrt(vcasingits) resolution per dimension) and halve the stride until the accuracy is reached. 
  !> At least one step in the resolution cascade is required to determine an estimate of the current accuracy.

  ! The starting stride is chosen such that it is the largest possible stride log2(sqrt(vcNt*vcNz)), but at least considers vcasingits points
  do vcstride = INT(0.5*log(real(vcNt*vcNz/vcasingits))/log(2.0)), 0, -1
    !$OMP PARALLEL DO SHARED(Dxyz, Nxyz, Pbxyz, Jxyz, ijreal, ijimag) FIRSTPRIVATE(jk, gBn) COLLAPSE(2)
    do kk = 0, Nzwithsym-1 ; 
      do jj = 0, Nt-1 ; 
        jk = 1 + jj + kk*Nt

        ! kk = 0 terms are also symmetryic (relevant e.g. for tokamaks) 
        if ((Nz.eq.1) .and. (jj.gt.((Nt+1)/2))) cycle 

        ! Each MPI rank only computes every a 1/ncpu surfacecurrent() calls 
        ! Identical MPI parallelization scheme as for Lvcgrid=0
        select case( Lparallel ) 
        case( 0 ) ! Lparallel = 0 
        if( myid.ne.modulo(kk,ncpu) ) cycle
        case( 1 ) ! Lparallel = 1 
        if( myid.ne.modulo(jk-1,ncpu) ) cycle 
        case default ! Lparallel; 
        FATAL( bnorml, .true., invalid Lparallel in parallelization loop )
        end select ! end of select case( Lparallel ) 

        gBn = zero
        globaljk = jk ! only to check against the precomputed values against dvcfield
        call casinggrid( Dxyz(1:3,jk), Nxyz(1:3,jk), Pbxyz, Jxyz, 2**vcstride,  gBn)
        
        ijimag(jk) = ijreal(jk) ! previous solution (lower resolution)
        ijreal(jk) = gBn ! current solution (higher resolution)
      enddo ! end of do jj
    enddo ! end of do kk

    deltah4h2 = deltah2h
    deltah2h =  maxval(abs(ijimag - ijreal)) ! mean delta between the h and h/2 solutions

    ! Order of the integration method: log(deltah4h2/deltah2h)/log(2.0) = 1
    absvcerr = deltah2h 
    relvcerr = 2 * deltah2h / maxval(abs(ijreal)) ! overestimate the relative error by a factor of two
    if (myid.eq.0) then
      write(ounit, '("bnorml : ", 10x ," : relvcerr = ",es13.5," ; absvcerr = ",es13.5," ; resolution reduction = ",i8)') relvcerr, absvcerr, 2**vcstride
      ! print *, "Convergence order: ",  log(deltah4h2/deltah2h)/log(2.0)
    endif
    ! Tolerance was already reached, exit the loop. (Different threads may exit at different times)
    if (absvcerr.lt.vcasingtol .and. relvcerr.lt.vcasingtol) then
      exit
    endif
  enddo ! end of do vcstride
#ifdef COMPARECASING
  ! To compare with the casing() implementation, copy the results into ijimag
  if(Lvcgrid.eq.1) then
    ijimag(1:Ntz) = ijreal(1:Ntz)
  endif
endif
! When comparing results, both methods should always run
if (.true.) then
#else
else ! if not Lvcgrid  
#endif

  do kk = 0, Nzwithsym-1 ; 
   zeta = kk * pi2 / Nz

   do jj = 0, Nt-1 ; 
    teta = jj * pi2 / Nt ; jk = 1 + jj + kk*Nt

    ! kk = 0 terms are also symmetryic (relevant e.g. for tokamaks) 
    if ((Nz.eq.1) .and. (jj.gt.((Nt+1)/2))) cycle 

    globaljk = jk ! this is global; passed through to vcintegrand & casing;

    select case( Lparallel ) ! perform in parallel;
    case( 0 ) ! Lparallel = 0 ; 09 Mar 17;
     if( myid.ne.modulo(kk,ncpu) ) cycle
    case( 1 ) ! Lparallel = 1 ; 09 Mar 17;
     if( myid.ne.modulo(jk-1,ncpu) ) cycle ! 11 Oct 12; this is a weird parallelization, but perhaps better exploits all available cpus;
    case default ! Lparallel; 09 Mar 17;
     FATAL( bnorml, .true., invalid Lparallel in parallelization loop )
    end select ! end of select case( Lparallel ) ; 09 Mar 17;

    WCALL( bnorml, casing, ( teta, zeta, gBn, icasing ) )
    ijreal(jk) = gBn

#ifdef COMPARECASING
  ! write(ounit,1000) myid, zeta, teta, ijreal(jk), ijimag(jk), ijreal(jk)-ijimag(jk)
  write(ounit,1000) myid, zeta, teta, ijreal(jk), ijimag(jk), ijreal(jk)-ijimag(jk)
  1000 format("bnorml : ", 10x ," : myid=",i3," : \z =",f6.3," ; \t =",f6.3," ; B . x_t x x_z =",2f22.15," ; ":"err =",es13.5," ;")
#endif

   enddo ! end of do jj;
  enddo ! end of do kk;
endif ! end of if (Lvcgrid  )
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
   case( 1 ) ! Lparallel = 1 ; 09 Mar 17;

    do jj = 0, Nt-1

     jk = 1 + jj + kk*Nt

     jkmodnp = modulo(jk-1,ncpu)

     RlBCAST(ijreal(jk),1,jkmodnp) ! plasma; 03 Apr 13;

    enddo

   case default ! Lparallel; 09 Mar 17;

    FATAL( bnorml, .true., invalid Lparallel for broadcasting )

   end select ! end of select case( Lparallel ) ; 09 Mar 17;

  enddo ! 11 Oct 12;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  ! In stellerator symmetric geometry, the virtual casing field also exhibits symmetriy which we can exploit
  if (YESstellsym) then
    ! We are dealing with half open intervals in theta and phi [0, 2pi[ so we need to skip the first row and column when mirroring 
    if (Nz.eq.1) then 
      do jj = 0, Nt/2; 
        jk = 1 + jj
        ijreal(Ntz + 1 - jk) = -ijreal(jk+1)
      enddo    
    endif

    do kk = 0, Nzwithsym-2; 
      do jj = 0, Nt-1;  
        jk = 1 + jj + kk*Nt

        ijreal(Ntz + 1 - jk) = -ijreal(jk+Nt+1)
      enddo
    enddo
  endif
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  ijreal(1:Ntz) = ijreal(1:Ntz) * virtualcasingfactor
  ijimag(1:Ntz) = zero

  call tfft( Nt, Nz, ijreal(1:Ntz), ijimag(1:Ntz), &
             mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), ifail ) ! Fourier decompose normal field;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(bnorml)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine bnorml
