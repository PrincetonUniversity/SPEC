!> \defgroup grp_fourier_mapping

module bndRep
    implicit none  
    PUBLIC ! everything is public, excepted stated otherwise.
  
    INTEGER              :: Mpol_field, Ntor_field, Mpol_force, Ntor_force, Mpol_max, Ntor_max
    INTEGER              :: nel_m1_i, nel_m1_j, nel_m2_i, nel_m2_j ! dimension of matrices
    LOGICAL              :: Lchangeangle
    REAL   , allocatable :: Mat1(:,:), Mat2(:,:)                !< MAPPING MATRIX, TIMES FOUR
    REAL   , allocatable :: RHS1(:), RHS2(:), LHS1(:), LHS2(:)  !< Right hand side and Left and side of linear systems
    REAL   , allocatable :: dRZdhenn(:,:)                       !< Derivatives of R,Z_mn w.r.t r,z,b_n and rho_mn
    REAL   , allocatable :: precond_rho(:,:), precond_b(:)      !< Similar to psifactor for Rmn, Zmn for rhomn

    !------- public / private statement ----------
  
    PRIVATE :: Mat1, Mat2
    PRIVATE :: nel_m1_i, nel_m1_j, nel_m2_i, nel_m2_j
    PRIVATE :: build_mapping_matrices ! only used in forwardMap and backwardMap subroutines.
    PRIVATE :: pack_rhomn_bn, unpack_rhomn_bn, pack_rmn_zmn, unpack_rmn_zmn ! specific for this module
    PRIVATE :: RHS1, RHS2, LHS1, LHS2
    PRIVATE :: Lchangeangle
  
    contains
  
  
  
    ! ------------------------------------------------------------------
    !                     PUBLIC SUBROUTINES
    
    !> \brief Initialize mapping matrices, and set truncation values
    !> \ingroup grp_fourier_mapping
    !>
    !> @param[in] Langle
    subroutine initialize_mapping( Langle )
      ! In this subroutine we compute the mapping matrix, and allocate necessary memory
      ! This should only be called once at the beginning of preset.

      use constants, only: zero, half, one
      use inputlist, only: Wmacros, Mpol, Ntor, twoalpha, Nfp, Lboundary, tflux
      use fileunits, only: ounit, lunit
      use allglobal, only: myid, mn_field, im_field, in_field, &
                           mn_rho, im_rho, in_rho, &
                           mn_force, im_force, in_force, &
                           mn_max, im_max, in_max, &
                           MPI_COMM_SPEC, cpus, &
                           Rscale, Mvol

      LOCALS

      LOGICAL, INTENT(IN):: Langle !< Input; used to determine if angle has to be reversed
      INTEGER            :: vvol   !< Loop index on volumes
      INTEGER            :: ii     !< Loop index on Fourier modes \f$\rho_{mn}\f$


      Lchangeangle = Langle

      !> <ul> Define Fourier series truncation
      !> <li> If \c Lboundary=0 (Hudson representation), all physics quantities
      !>      are truncated at \c Mpol and \c Ntor
      !> <li> If \c Lboundary=1 (Henneberg representation), the truncation has
      !>      to be selected carefully so that the number of boundary dofs is 
      !>      equal to the number of force harmonics. We set truncation as follows:
      !>      <ul>
      !>      <li> \f$\rho_{mn},\ b_n,\ r_n,\ z_n\f$ at \c Mpol and \c Ntor (input)
      !>      <li> The \f$(R_{mn},Z_{mn})\f$ harmonics, as well as all \c Ate, \c Ato, \c Aze and \c Azo at
      !>           \f$Mpol_{field} = Mpol\f$ and  \f$Ntor_{field} = Ntor+2|\alpha|\f$
      !>      <li> The force at \f$Mpol_{force}=Mpol+1\f$, \f$Ntor_{force}=Ntor\f$
      !>      </ul>
      !> </ul>
      if( Lboundary.eq.0 ) then  
        Mpol_field = Mpol
        Ntor_field = Ntor
  
        Mpol_force = Mpol
        Ntor_force = Ntor

        Mpol_max   = Mpol
        Ntor_max   = Ntor

      elseif( Lboundary.eq.1 ) then
        !Mpol_field = Mpol + 1
        Mpol_field = Mpol
        Ntor_field = Ntor + abs(twoalpha)
  
        Mpol_force = Mpol + 1
        !Ntor_force = Ntor + abs(twoalpha)
        Ntor_force = Ntor

        Mpol_max = Mpol+1
        Ntor_max = Ntor+abs(twoalpha)
        
      else 
        FATAL( bndRep, .true., Invalid Lboundary )

      endif


      mn_field = 1 + Ntor_field +  Mpol_field  * ( 2 *  Ntor_field  + 1 ) ! Fourier resolution of interface geometry & vector potential;
      mn_rho   =                   Mpol        * ( 2 *  Ntor        + 1 )
      mn_force = 1 + Ntor_force +  Mpol_force  * ( 2 *  Ntor_force  + 1 ) ! Fourier resolution of interface geometry & vector potential;
      mn_max   = 1 + Ntor_max   +  Mpol_max    * ( 2 *  Ntor_max    + 1 ) ! Fourier resolution of interface geometry & vector potential;

      !latex \subsubsection{\type{mn}, \type{im(1:mn)} and \type{in(1:mn)} : Fourier mode identification}
      !latex \begin{enumerate}
      !latex \item The Fourier description of even periodic functions is
      !latex       \be f(\t,\z) = \sum_{n=0}^{N} f_{0,n} \cos(-n\z) + \sum_{m=1}^{M}\sum_{n=-N}^{N} f_{m,n} \cos(m\t-n\z),
      !latex       \ee
      !latex       where the resolution is given on input, $M\equiv $ \inputvar{ Mpol} and $N\equiv $ \inputvar{ Ntor}.
      !latex \item For convenience, the Fourier summations are written as
      !latex       \be f(\s,\t,\z) &=& \sum_j f_j(s) \cos( m_j \t - n_j \z ),
      !latex       \ee
      !latex       for $j=1,$ \type{mn}, where \type{mn}$ = N + 1 +  M  ( 2 N + 1 )$.
      !latex \item The integer arrays \type{im(1:mn)} and \type{in(1:mn)} contain the $m_j$ and $n_j$.
      !latex \item The array \type{in} includes the \type{Nfp} factor.
      !latex \end{enumerate}
      
      SALLOCATE( im_field  , (1:mn_field  ), zero )
      SALLOCATE( in_field  , (1:mn_field  ), zero )
      SALLOCATE( im_rho    , (1:mn_rho    ), zero )
      SALLOCATE( in_rho    , (1:mn_rho    ), zero )
      SALLOCATE( im_force  , (1:mn_force  ), zero )
      SALLOCATE( in_force  , (1:mn_force  ), zero )
      SALLOCATE( im_max    , (1:mn_max    ), zero )
      SALLOCATE( in_max    , (1:mn_max    ), zero )
    
      call gi00ab(  Mpol_field,  Ntor_field, Nfp, mn_field, im_field(1:mn_field), in_field(1:mn_field), .true.  ) ! this sets the im and in mode identification arrays;
      call gi00ab(  Mpol_force,  Ntor_force, Nfp, mn_force, im_force(1:mn_force), in_force(1:mn_force), .true.  ) ! this sets the im and in mode identification arrays;
      call gi00ab(  Mpol      ,  Ntor      , Nfp, mn_rho  , im_rho(  1:mn_rho  ), in_rho(  1:mn_rho  ), .false. ) ! this sets the im and in mode identification arrays;
      call gi00ab(  Mpol_max  ,  Ntor_max  , Nfp, mn_max  , im_max(  1:mn_max  ), in_max(  1:mn_max  ), .false. ) ! this sets the im and in mode identification arrays;

      nel_m1_i = 2*(2*Ntor_field+1)
      nel_m1_j = 3*Ntor+2
      SALLOCATE( Mat1, (1:nel_m1_i,1:nel_m1_j), zero )
      SALLOCATE( RHS1, (1:nel_m1_j           ), zero )
      SALLOCATE( LHS1, (1:nel_m1_i           ), zero )

      nel_m2_i = 2*(Mpol_field-1)*(2*Ntor_field+1)
      nel_m2_j = (Mpol-1)*(2*Ntor+1)
      SALLOCATE( Mat2, (1:nel_m2_i,1:nel_m2_j), zero )
      SALLOCATE( RHS2, (1:nel_m2_j           ), zero )
      SALLOCATE( LHS2, (1:nel_m2_i           ), zero )

      call build_mapping_matrices()

      SALLOCATE( precond_rho, (1:mn_rho, 1:Mvol), one )
      do vvol = 1, Mvol
        do ii = 1, mn_rho
          precond_rho( ii, vvol ) = Rscale * tflux(vvol)**(im_rho(ii) * half) 
        enddo
      enddo

      SALLOCATE( precond_b, (1:Mvol), one )
      do vvol = 1, Mvol
        precond_b( vvol ) = Rscale * tflux(vvol)**half
      enddo

    end subroutine initialize_mapping
  

      subroutine change_mapping_angle( Langle )

        use fileunits, only: ounit, lunit

        LOCALS

        LOGICAL, INTENT(IN):: Langle

        Lchangeangle = Langle


      end subroutine change_mapping_angle
  
  
      subroutine forwardMap( rhomn, bn, R0c, Z0s, Rmn, Zmn )
  
        use constants, only: zero
        use inputlist, only: Wmacros, Mpol, Ntor, twoalpha, Nfp
        use fileunits, only: ounit, lunit
        use allglobal, only: myid, im_rho, in_rho, im_field, in_field, mn_rho, mn_field, &
                             MPI_COMM_SPEC, cpus

        LOCALS

        ! INPUTS
        REAL, intent(in) :: R0c(0:Ntor), Z0s(1:Ntor)
        REAL, intent(in) :: rhomn(1:mn_rho), bn(0:Ntor)
        REAL             :: Rmn(1:mn_field), Zmn(1:mn_field)
        INTEGER          :: ii, mm, nn

        ! m=0 modes, from R0c and Z0c
        do ii=1,mn_field
          mm=im_field(ii)
          nn=in_field(ii) / Nfp

          if( mm.ne.0    .or. nn.lt.0    ) cycle

          if( nn.gt.Ntor ) then
            Rmn(ii) = zero
          else
            Rmn(ii) = R0c(nn)
          endif

          if( (nn.eq.0) .or. (nn.gt.Ntor) ) then
            Zmn(ii) = zero
          else
            Zmn(ii) =-Z0s(nn)
          endif
        enddo
  
        ! m=1 modes
        call pack_rhomn_bn( rhomn(1:mn_rho), bn(0:Ntor) )
        LHS1 = MATMUL( Mat1, RHS1 )
  
        ! m>1 modes
        if( Mpol_field.gt.1 ) then
          LHS2 = MATMUL( Mat2, RHS2 )
        endif        
  
        call unpack_rmn_zmn( Rmn(1:mn_field), Zmn(1:mn_field) )


        ! Change angle
        if( Lchangeangle ) then
          call change_angle( Rmn(1:mn_field), Zmn(1:mn_field) )
        endif
  
      end subroutine forwardMap
  
  
      subroutine backwardMap( Rmn, Zmn, rhomn, bn, R0c, Z0s )
  
        use constants, only: zero
        use inputlist, only: Wmacros, Mpol, Ntor, twoalpha, Nfp
        use fileunits, only: ounit, lunit
        use allglobal, only: myid, im_rho, in_rho, im_field, in_field, mn_rho, mn_field, &
                             MPI_COMM_SPEC, cpus

        LOCALS

        ! INPUTS
        REAL, intent(in)    :: Rmn(1:mn_field), Zmn(1:mn_field)
        REAL, intent(inout) :: rhomn(1:mn_rho), bn(0:Ntor)
        REAL, intent(inout) :: R0c(0:Ntor), Z0s(1:Ntor)


        ! LOCAL VARIABLES
        CHARACTER           :: TRANS
        REAL                :: Rwork(1:mn_field), Zwork(1:mn_field)
        INTEGER             :: ii, mm, nn
        REAL, allocatable   :: A(:,:), WORK(:), B(:)
        INTEGER             :: NRHS, LDA, LDB, LWORK, INFO

        ! First, change angle if necessary.
        Rwork = Rmn
        Zwork = Zmn

        if( Lchangeangle ) then
          call change_angle( Rwork(1:mn_field), Zwork(1:mn_field) )
        endif
  
        ! m=0 modes - set R0c and Z0s
        do ii=1,mn_field
          mm = im_field(ii)
          nn = in_field(ii) / Nfp


          if( mm.ne.0    .or. nn.lt.0    ) cycle ! Only interested in m=0, n>=0 modes
          if( mm.gt.Mpol .or. nn.gt.Ntor ) cycle ! other elements are zeros by construction. otherwise seg fault.

          R0c(nn) = Rwork(ii)
          if( nn.ne.0 ) then
            Z0s(nn) =-Zwork(ii)
          endif
        enddo

        ! m>0 modes
        SALLOCATE( A, (1:nel_m1_i,1:nel_m1_j), zero )  
        SALLOCATE( B, (1:nel_m1_i           ), zero )    
  
        call DCOPY(nel_m1_i*nel_m1_j, Mat1, 1, A, 1)
        call pack_rhomn_bn( rhomn(1:mn_rho), bn(0:Ntor) )
        TRANS = 'N'
        NRHS = 1
        LDA = nel_m1_i
        LDB = nel_m1_i
  
        ! Work query to get optimal block size
        LWORK = -1
        SALLOCATE( WORK, (1:1), zero )
        call DGELS( TRANS, nel_m1_i, nel_m1_j, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )
  
        select case( INFO )
        case( 0 )
          LWORK = WORK(1)
  
        case( :-1)
          FATAL( bndRep, .true., Illegal value in DGELS )
  
        case(1: )
          FATAL( bndRep, .true., Rank zero for backward mapping )
  
        end select
        DALLOCATE( WORK )
  
  
        ! Now solve linear system
        call pack_rmn_zmn( Rwork, Zwork )

        SALLOCATE( WORK, (1:LWORK), zero )
        LDA = nel_m1_i
        LDB = nel_m1_i
        call DCOPY( nel_m1_i*nel_m1_j, Mat1, 1, A, 1)
        B(1:nel_m1_i) = LHS1(1:nel_m1_i)
        call DGELS( TRANS, nel_m1_i, nel_m1_j, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )
        
        select case( INFO )
        case( 0 )
          RHS1(1:nel_m1_j) = B(1:nel_m1_j)
  
        case( :-1)
          FATAL( bndRep, .true., Illegal value in DGELS )
  
        case(1: )
          FATAL( bndRep, .true., Rank zero for backward mapping )
  
        end select
        
        DALLOCATE( A )
        DALLOCATE( B )
        DALLOCATE( WORK )

        
        ! Now solve second system
        if( Mpol>1 ) then
          ! First figure out optimal block size - see dgels documentation 
          LWORK = -1
          LDA = nel_m2_i
          LDB = nel_m2_i
          SALLOCATE( A, (1:nel_m2_i,1:nel_m2_j), zero )
          SALLOCATE( B, (1:nel_m2_i), zero )
          SALLOCATE( WORK, (1:1), zero )
          call DGELS( TRANS, nel_m2_i, nel_m2_j, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )

  
          select case( INFO )
          case( 0 )
            LWORK = WORK(1)
    
          case( :-1)
            FATAL( bndRep, .true., Illegal value in DGELS )
    
          case(1: )
            FATAL( bndRep, .true., Rank zero for backward mapping )
    
          end select
          DALLOCATE( WORK )

          ! Then solve system
          SALLOCATE( WORK, (1:LWORK), zero )
          call DCOPY( nel_m2_i*nel_m2_j, Mat2, 1, A, 1 )
          B(1:nel_m2_i) = LHS2(1:nel_m2_i)
          call DGELS( TRANS, nel_m2_i, nel_m2_j, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )
          select case( INFO )
          case( 0 )
            RHS2(1:nel_m2_j) = B(1:nel_m2_j)
    
          case( :-1)
            FATAL( bndRep, .true., Illegal value in DGELS )
    
          case(1: )
            FATAL( bndRep, .true., Rank zero for backward mapping )
    
          end select
          DALLOCATE( A )
          DALLOCATE( B )
          DALLOCATE( WORK )
        endif


        ! Unpack RHS1, RHS2 and store in rhomn, bn.
        call unpack_rhomn_bn( rhomn, bn )

 
      end subroutine backwardMap

      subroutine initialize_force_gradient_transformation( iisize, jjsize )

        use constants, only: zero, one, quart, half
        use fileunits, only: ounit, lunit
        use inputlist, only: Igeometry, Nfp, twoalpha, Wmacros, Ntor, Lcheck
        use allglobal, only: YESstellsym, Mvol, LGdof_field, LGdof_bnd, &
                             im_rho, in_rho, mn_rho, &
                             im_field, in_field, mn_field, &
                             MPI_COMM_SPEC, cpus, myid, &
                             Rscale, &
                             ext, NGdof_bnd, NGdof_field !Only for debug purposes.

        LOCALS

        INTEGER, INTENT(IN) :: iisize, jjsize
        INTEGER             :: vvol                 ! volume index
        INTEGER             :: ii                   ! Fourier harmonics of field and hudson representation index
        INTEGER             :: irz                  ! 0: Rmn harmonics, 1: Zmn harmonics
        INTEGER             :: issym                ! 0: stellarator symmetric terms, 1: non-stellarator symmetric terms
        INTEGER             :: idof                 ! degree of freedom on the interface, hudson representation
        INTEGER             :: jdof                 ! degree of freedom on the interface, henneberg representation
        INTEGER             :: tdof                 ! global geometrical degree of freedom, hudson representation
        INTEGER             :: udof                 ! global geometrical degree of freedom, henneberg representation
        INTEGER             :: mm, nn               ! poloidal, toroidal mode number
        INTEGER             :: sgn_delta_twoalpha   ! Sign dependence of some terms.. + or -1
        INTEGER             :: kk, ll               ! Set to -1 if Lchangeangle is set to true

#ifdef DEBUG
        REAL                :: orhoc(1:mn_rho, 1:Mvol), obc(0:Ntor, 1:Mvol), orn(0:Ntor, 1:Mvol), ozn(0:Ntor, 1:Mvol)
        REAL                :: oRbc(1:mn_field, 0:Mvol), oZbs(1:mn_field, 0:Mvol), oRbs(1:mn_field, 0:Mvol), oZbc(1:mn_field, 0:Mvol)
        REAL                :: iposition(0:NGdof_field)
        REAL                :: dRZdhenn_debug( 1:NGdof_field, 1:NGdof_bnd )
        INTEGER             :: jj
#endif

        ! --------------------------------------------------------------------------------------

        if( Lcheck.eq.62 ) then
          precond_rho( 1:mn_rho, 1:Mvol ) = one
          precond_b( 1:Mvol ) = one
          Rscale = 1
        endif

        if( Lchangeangle ) then 
          kk = -1
        else
          kk = +1
        endif


        ! Allocate memory
        SALLOCATE( dRZdhenn, (1:iisize, 1:jjsize), zero )

        ! Build matrix
        do vvol = 1, Mvol-1 ! loop over interior surfaces;
          idof = 0
      
          do ii = 1, mn_field ! Loop over Fourier modes
      
            do irz = 0, 1 ! loop over R or Z coordinate
      
              if( irz.eq.1 .and. Igeometry.lt.3 ) cycle

              mm = im_field( ii )
              nn = in_field( ii ) / Nfp

              ll = 1
              if( irz.eq.0 ) then
                sgn_delta_twoalpha = -1
              else
                sgn_delta_twoalpha =  1
                if( Lchangeangle ) ll = -1
              endif
      
              do issym = 0, 1 ! stellarator symmetry;
      
                if( issym.eq.1 .and. YESstellsym ) cycle ! no dependence on the non-stellarator symmetric harmonics;
      
                if( ii.eq.1 .and. irz.eq.1 .and. issym.eq.0 ) cycle ! no dependence on Zbs_{m=0,n=0};
                if( ii.eq.1 .and. irz.eq.0 .and. issym.eq.1 ) cycle ! no dependence on Rbs_{m=0,n=0};
      

                idof = idof + 1 ! labels degree-of-freedom;
                tdof = (vvol-1) * LGdof_field + idof

                ! dR / drn
                if( irz.eq.0 .and. mm.eq.0 .and. nn.ge.0 .and. nn.le.Ntor ) then
                  if( nn.eq.0 ) then
                    jdof = mn_rho + 2
                  else
                    jdof = mn_rho + 3*nn + 1
                  endif
#ifdef DEBUG
                  FATAL( bndRep, jdof.gt.LGdof_bnd .or. jdof.lt.1, Index is out of bound 001 )
#endif

                  udof = (vvol-1) * LGdof_bnd + jdof

                  dRZdhenn( tdof, udof ) = one * Rscale
                endif

                ! dZ / dn
                if( irz.eq.1 .and. mm.eq.0 .and. nn.ge.0 .and. nn.le.Ntor  ) then
                  if( nn.ne.0 ) then
                    jdof = mn_rho + 3*nn + 2
                  endif
#ifdef DEBUG
                  FATAL( bndRep, jdof.gt.LGdof_bnd .or. jdof.lt.1, Index is out of bound 002 )
#endif      
                  udof = (vvol-1) * LGdof_bnd + jdof

                  dRZdhenn( tdof, udof ) = -one * Rscale
                endif

                ! dR / dbn    and     dZ / dbn
                if( mm.eq.1 ) then                     
                  if( kk*nn.ge.-Ntor .and. kk*nn.le.-1 ) then
                    jdof = mn_rho - 3*nn*kk
                    udof = (vvol-1) * LGdof_bnd + jdof                  
                    dRZdhenn( tdof, udof ) = dRZdhenn( tdof, udof ) + quart * precond_b( vvol ) * ll
                  elseif( kk*nn.ge.1 .and. kk*nn.le.Ntor) then
                    jdof = mn_rho + 3*nn*kk
                    udof = (vvol-1) * LGdof_bnd + jdof                  
                    dRZdhenn( tdof, udof ) = dRZdhenn( tdof, udof ) + quart * precond_b( vvol ) * ll
                  elseif( kk*nn.eq. 0) then
                    jdof = mn_rho + 1
                    udof = (vvol-1) * LGdof_bnd + jdof                  
                    dRZdhenn( tdof, udof ) = dRZdhenn( tdof, udof ) + half * precond_b( vvol ) * ll
                  endif

                  if( kk*nn-twoalpha.gt.0 .and. kk*nn-twoalpha.le.Ntor ) then
                    jdof = mn_rho + 3*( kk*nn-twoalpha )
                    udof = (vvol-1) * LGdof_bnd + jdof                  
                    dRZdhenn( tdof, udof ) = dRZdhenn( tdof, udof ) + sgn_delta_twoalpha * quart * precond_b( vvol ) * ll
                  elseif( kk*nn-twoalpha.eq.0 ) then
                    jdof = mn_rho + 1
                    udof = (vvol-1) * LGdof_bnd + jdof                  
                    dRZdhenn( tdof, udof ) = dRZdhenn( tdof, udof ) + sgn_delta_twoalpha * quart * precond_b( vvol ) * ll
                  endif

                  if(-kk*nn+twoalpha.gt.0 .and. -kk*nn+twoalpha.le.Ntor ) then
                    jdof = mn_rho + 3*(-kk*nn+twoalpha)
                    udof = (vvol-1) * LGdof_bnd + jdof                  
                    dRZdhenn( tdof, udof ) = dRZdhenn( tdof, udof ) + sgn_delta_twoalpha * quart * precond_b( vvol ) * ll
                  elseif(-kk*nn+twoalpha.eq.0 ) then
                    jdof = mn_rho + 1
                    udof = (vvol-1) * LGdof_bnd + jdof                  
                    dRZdhenn( tdof, udof ) = dRZdhenn( tdof, udof ) + sgn_delta_twoalpha * quart * precond_b( vvol ) * ll
                  endif
                endif


                ! dR / drhomn     and     dZ / drhomn
                if( mm.ne.0 ) then
                  call find_index(  mn_rho, im_rho, in_rho / Nfp, mm, -kk*nn, jdof )
                  if( jdof.gt.0 ) then
                    udof = (vvol-1) * LGdof_bnd + jdof                  
                    dRZdhenn( tdof, udof ) = half * precond_rho( jdof, vvol ) * ll
                  endif


                  call find_index(  mn_rho, im_rho, in_rho / Nfp, mm, -kk*nn+twoalpha, jdof )
                  if( jdof.gt.0 ) then
                    udof = (vvol-1) * LGdof_bnd + jdof                  
                    dRZdhenn( tdof, udof ) = dRZdhenn( tdof, udof ) - sgn_delta_twoalpha * half * precond_rho( jdof, vvol ) * ll
                  endif
                endif
              enddo !issym
            enddo !irz
          enddo !ii
        enddo !vvol

#ifdef DEBUG

        if( Lcheck.eq.9 ) then
          idof = 0
          do vvol=1,Mvol-1
            do jj=1,mn_rho
              idof = idof+1

              orhoc(1:mn_rho, 1:Mvol)  = zero
              obc(0:Ntor, 1:Mvol)      = zero
              orn(0:Ntor, 1:Mvol)      = zero
              ozn(0:Ntor, 1:Mvol)      = zero
              oRbc(1:mn_field, 0:Mvol) = zero
              oZbs(1:mn_field, 0:Mvol) = zero
              oRbs(1:mn_field, 0:Mvol) = zero
              oZbc(1:mn_field, 0:Mvol) = zero
  
              ! Take derivative of mapping matrix equation
              orhoc(jj,vvol)           = precond_rho( jj, vvol )
              call forwardMap( orhoc(1:mn_rho, vvol), obc(0:Ntor, vvol), orn(0:Ntor, vvol), &
                               ozn(1:Ntor, vvol), oRbc(1:mn_field, vvol), oZbs(1:mn_field, vvol))

              ! Pack in dofs
              call packxi( NGdof_field, iposition(0:NGdof_field), Mvol, mn_field, &
                           oRbc(1:mn_field, 0:Mvol), oZbs(1:mn_field, 0:Mvol), &
                           oRbs(1:mn_field, 0:Mvol), oZbc(1:mn_field, 0:Mvol), &
                           'P', .false., .false. )
              
              ! Store derivatives
              dRZdhenn_debug( 1:NGdof_field, idof ) = iposition( 1:NGdof_field )

            enddo
  
            do jj=0,Ntor
              ! --- derivatives w.r.t b_n
              idof = idof+1
              orhoc(1:mn_rho, 1:Mvol)  = zero
              obc(0:Ntor, 1:Mvol)      = zero
              orn(0:Ntor, 1:Mvol)      = zero
              ozn(0:Ntor, 1:Mvol)      = zero
              oRbc(1:mn_field, 0:Mvol) = zero
              oZbs(1:mn_field, 0:Mvol) = zero
              oRbs(1:mn_field, 0:Mvol) = zero
              oZbc(1:mn_field, 0:Mvol) = zero

              ! Take derivative of mapping matrix equation
              obc(jj, vvol) = precond_b( vvol )
              call forwardMap( orhoc(1:mn_rho, vvol), obc(0:Ntor, vvol), orn(0:Ntor, vvol), &
                               ozn(1:Ntor, vvol), oRbc(1:mn_field, vvol), oZbs(1:mn_field, vvol))

              ! Pack in dofs
              call packxi(  NGdof_field, iposition(0:NGdof_field), Mvol, mn_field, &
                            oRbc(1:mn_field, 0:Mvol), oZbs(1:mn_field, 0:Mvol), &
                            oRbs(1:mn_field, 0:Mvol), oZbc(1:mn_field, 0:Mvol), &
                            'P', .false., .false. )
              ! Store derivatives
              dRZdhenn_debug( 1:NGdof_field, idof ) = iposition( 1:NGdof_field )


              ! --- derivatives w.r.t r_n
              idof = idof+1
              orhoc(1:mn_rho, 1:Mvol)  = zero
              obc(0:Ntor, 1:Mvol)      = zero
              orn(0:Ntor, 1:Mvol)      = zero
              ozn(0:Ntor, 1:Mvol)      = zero
              oRbc(1:mn_field, 0:Mvol) = zero
              oZbs(1:mn_field, 0:Mvol) = zero
              oRbs(1:mn_field, 0:Mvol) = zero
              oZbc(1:mn_field, 0:Mvol) = zero
  
              ! Take derivative of mapping matrix equation
              orn(jj, vvol) = Rscale
              call forwardMap( orhoc(1:mn_rho, vvol), obc(0:Ntor, vvol), orn(0:Ntor, vvol), &
                               ozn(1:Ntor, vvol), oRbc(1:mn_field, vvol), oZbs(1:mn_field, vvol))

              ! Pack in dofs
              call packxi(  NGdof_field, iposition(0:NGdof_field), Mvol, mn_field, &
                            oRbc(1:mn_field, 0:Mvol), oZbs(1:mn_field, 0:Mvol), &
                            oRbs(1:mn_field, 0:Mvol), oZbc(1:mn_field, 0:Mvol), &
                            'P', .false., .false. )

              ! Store derivatives
              dRZdhenn_debug( 1:NGdof_field, idof ) = iposition( 1:NGdof_field )


              ! --- derivatives w.r.t r_n
              if( jj.gt.0 ) then 
                idof = idof+1
                orhoc(1:mn_rho, 1:Mvol)  = zero
                obc(0:Ntor, 1:Mvol)      = zero
                orn(0:Ntor, 1:Mvol)      = zero
                ozn(0:Ntor, 1:Mvol)      = zero
                oRbc(1:mn_field, 0:Mvol) = zero
                oZbs(1:mn_field, 0:Mvol) = zero
                oRbs(1:mn_field, 0:Mvol) = zero
                oZbc(1:mn_field, 0:Mvol) = zero

                ! Take derivative of mapping matrix equation
                ozn(jj, vvol) = Rscale
                call forwardMap( orhoc(1:mn_rho, vvol), obc(0:Ntor, vvol), orn(0:Ntor, vvol), &
                                 ozn(1:Ntor, vvol), oRbc(1:mn_field, vvol), oZbs(1:mn_field, vvol))
  
                ! Pack in dofs
                call packxi(  NGdof_field, iposition(0:NGdof_field), Mvol, mn_field, &
                              oRbc(1:mn_field, 0:Mvol), oZbs(1:mn_field, 0:Mvol), &
                              oRbs(1:mn_field, 0:Mvol), oZbc(1:mn_field, 0:Mvol), &
                              'P', .false., .false. )
                              
                ! Store derivatives
                dRZdhenn_debug( 1:NGdof_field, idof ) = iposition( 1:NGdof_field )  
              endif
            enddo
          enddo !vvol

          open(10, file=trim(ext)//'.dRZdhenn.txt', status='unknown')
          do ii=1,NGdof_field
            write(10, 2801) dRZdhenn( ii, : )
          enddo
          close(10)

          open(10, file=trim(ext)//'.dRZdhenn_debug.txt', status='unknown')
          do ii=1,NGdof_field
            write(10, 2801) dRZdhenn_debug( ii, : )
          enddo
          close(10)

2801  format(512F22.16, " ")

          FATAL( bndRep, .true., End of Lcheck.eq.9 )

        endif
#endif


      end subroutine initialize_force_gradient_transformation

      subroutine pack_henneberg_to_hudson( position, bndDofs )

        use constants, only: zero
        use inputlist, only: Wmacros, Mpol, Ntor, twoalpha, Nfp
        use fileunits, only: ounit, lunit
        use allglobal, only: myid, im_rho, in_rho, im_field, in_field, mn_rho, mn_field, &
                             MPI_COMM_SPEC, cpus, Mvol, &
                             NGdof_field, NGdof_bnd, &
                             IsMyVolume, IsMyVolumeValue, WhichCpuID, &
                             Rscale, NOTstellsym

        LOCALS

        ! --------------------------------------------------------------
        ! Local Variables

        REAL, INTENT(OUT)      :: position( 0:NGdof_field )
        REAL, INTENT(IN)       :: bndDofs( 0:NGdof_bnd )
        CHARACTER              :: packorunpack ! Either 'RZ_to_H' or 'H_to_RZ'
        LOGICAL                :: LComputeDerivatives, LComputeAxis, work
        INTEGER                :: lvol, idof, jj, cpu_id
        REAL                   :: iRbc(1:mn_field,0:Mvol), iRbs(1:mn_field,0:Mvol), iZbc(1:mn_field,0:Mvol), iZbs(1:mn_field,0:Mvol)
        REAL                   :: irhoc(1:mn_rho,1:Mvol), ibc(0:Ntor,1:Mvol), iR0c(0:Ntor,1:Mvol), iZ0s(1:Ntor,1:Mvol)

        ! First unpack the array bndDofs
        idof = 0

        do lvol=1,Mvol-1

          ! Each CPU does the mapping for its own interface
          call IsMyVolume(lvol)

          if( IsMyVolumeValue.EQ.0 ) then
            cycle
          elseif( IsMyVolumeValue.EQ.-1 ) then
            FATAL( bndRep, .true., Unassociated volume)
          endif

          do jj=1,mn_rho

            idof = idof+1
#ifdef DEBUG
            FATAL( bndRep, idof.le.0 .or. idof.gt.NGdof_bnd, out of bounds )
#endif

            irhoc(jj,lvol) = bndDofs(idof) * precond_rho( jj, lvol )
          enddo

          do jj=0,Ntor
            idof = idof+1
            ibc(jj, lvol) = bndDofs(idof) * precond_b( lvol )
            
            idof = idof+1
            iR0c(jj, lvol) = bndDofs(idof) * Rscale

            if( jj.gt.0 ) then 
              idof = idof+1
#ifdef DEBUG
              FATAL( bndRep, idof.le.0 .or. idof.gt.NGdof_bnd, out of bounds )
#endif
              iZ0s(jj, lvol) = bndDofs(idof) * Rscale
            endif
          enddo

        
          call forwardMap( irhoc(1:mn_rho, lvol), ibc(0:Ntor, lvol), &
                          iR0c(0:Ntor, lvol), iZ0s(1:Ntor, lvol), &
                          iRbc(1:mn_field, lvol), iZbs(1:mn_field, lvol) )
  
          FATAL( bndRep, idof.ne.NGdof_bnd, incorrect number of dofs. )
        enddo

        ! Broadcast the interface Fourier harmonics
        do lvol=1,Mvol-1
          call WhichCpuID(lvol, cpu_id)
          RlBCAST( iRbc(1:mn_field, lvol), mn_field, cpu_id )
          RlBCAST( iZbs(1:mn_field, lvol), mn_field, cpu_id )

          if( NOTstellsym ) then
            RlBCAST( iRbc(1:mn_field, lvol), mn_field, cpu_id )
            RlBCAST( iRbc(1:mn_field, lvol), mn_field, cpu_id )
          endif
        enddo
      

        ! Finally, build position array
        packorunpack = 'P'
        LcomputeDerivatives = .FALSE.
        LComputeAxis = .FALSE.
        call packxi( NGdof_field, position( 0:NGdof_field ), Mvol, mn_field, &
                     iRbc(1:mn_field,0:Mvol), iZbs(1:mn_field,0:Mvol), iRbs(1:mn_field,0:Mvol), &
                     iZbc(1:mn_field,0:Mvol), packorunpack, LComputeDerivatives, LComputeAxis )

      end subroutine pack_henneberg_to_hudson


      

      subroutine pack_hudson_to_henneberg( position, bndDofs )

        use constants, only: zero
        use inputlist, only: Wmacros, Mpol, Ntor, twoalpha, Nfp
        use fileunits, only: ounit, lunit
        use allglobal, only: myid, im_rho, in_rho, im_field, in_field, mn_rho, mn_field, &
                             MPI_COMM_SPEC, cpus, Mvol, &
                             NGdof_field, NGdof_bnd, &
                             Rscale, &
                             IsMyVolume, IsMyVolumeValue, WhichCpuID

        LOCALS

        ! --------------------------------------------------------------
        ! Local Variables

        REAL, INTENT(IN)       :: position( 0:NGdof_field )
        REAL, INTENT(OUT)      :: bndDofs( 0:NGdof_bnd )
        CHARACTER              :: packorunpack ! Either 'RZ_to_H' or 'H_to_RZ'
        LOGICAL                :: LComputeDerivatives, LComputeAxis, work
        INTEGER                :: lvol, idof, jj, cpu_id
        REAL                   :: iRbc(1:mn_field,0:Mvol), iRbs(1:mn_field,0:Mvol), iZbc(1:mn_field,0:Mvol), iZbs(1:mn_field,0:Mvol)
        REAL                   :: irhoc(1:mn_rho,1:Mvol), ibc(0:Ntor,1:Mvol), iR0c(0:Ntor,1:Mvol), iZ0s(1:Ntor,1:Mvol)
          
        
        
        ! First unpack position in Rmn, Zmn
        packorunpack = 'U'
        LComputeDerivatives = .FALSE.
        LComputeAxis = .FALSE.
        call packxi( NGdof_field, position( 0:NGdof_field ), Mvol, mn_field, & 
                     iRbc(1:mn_field,0:Mvol), iZbs(1:mn_field,0:Mvol), iRbs(1:mn_field,0:Mvol), &
                     iZbc(1:mn_field,0:Mvol), packorunpack, LComputeDerivatives, LComputeAxis )

        ! Then map to rhomn, bn, R0n, Z0n
        do lvol=1,Mvol-1
          ! Each CPU does the backward mapping for its own interface
          call IsMyVolume(lvol)
          
          if( IsMyVolumeValue.EQ.0 ) then
            cycle
          elseif( IsMyVolumeValue.EQ.-1 ) then
            FATAL( bndRep, .true., Unassociated volume)
          endif

          call backwardMap( iRbc(1:mn_field,lvol), iZbs(1:mn_field,lvol), &
                            irhoc(1:mn_rho,lvol), ibc(0:Ntor, lvol), &
                            iR0c(0:Ntor, lvol), iZ0s(1:Ntor, lvol) )      
        enddo

        do lvol=1,Mvol-1
          call WhichCpuID(lvol, cpu_id)
          RlBCAST( irhoc(1:mn_rho, lvol), mn_rho, cpu_id )
          RlBCAST( ibc( 0:Ntor, lvol), Ntor+1, cpu_id )
          RlBCAST( iR0c(0:Ntor, lvol), Ntor+1, cpu_id )
          RlBCAST( iZ0s(1:Ntor, lvol), Ntor  , cpu_id )
        enddo

        ! Finally construct the array bndDofs
        idof = 0
        do lvol=1,Mvol-1
          do jj=1,mn_rho
            idof = idof + 1
#ifdef DEBUG
            FATAL( bndRep, idof.le.0 .or. idof.gt.NGdof_bnd, out of bounds )
#endif

            bndDofs(idof) = irhoc(jj,lvol) / precond_rho( jj, lvol )

          enddo

          do jj=0,Ntor
            idof = idof+1
            bndDofs(idof) = ibc(jj, lvol) / precond_b( lvol )

            idof = idof+1
            bndDofs(idof) = iR0c(jj, lvol) / Rscale

            if( jj.ne.0 ) then
              idof = idof+1
#ifdef DEBUG
              FATAL( bndRep, idof.le.0 .or. idof.gt.NGdof_bnd, out of bounds )
#endif
              bndDofs(idof) = iZ0s(jj, lvol) / Rscale
            endif
          enddo
        enddo !lvol
      end subroutine pack_hudson_to_henneberg
  
  













    ! ------------------------------------------------------------------
    !                     PRIVATE SUBROUTINES
  

      ! This routine pack rhomn, bn into RHS1, RHS2
      subroutine pack_rhomn_bn( rhomn, bn )
        use constants, only: zero
        use inputlist, only: Wmacros, Mpol, Ntor, twoalpha, Nfp
        use fileunits, only: ounit, lunit
        use allglobal, only: myid, im_rho, in_rho, im_field, in_field, mn_rho, mn_field, &
                             MPI_COMM_SPEC, cpus

        LOCALS

        REAL, intent(in)  :: rhomn( 1:mn_rho )
        REAL, intent(in)  :: bn( 0:Ntor )

        REAL  :: rho_work(1:Mpol,-Ntor:Ntor)

        INTEGER :: ii, mm, nn

        ! ---------------------------------------------------------------
  
        ! Ensure working variables are set to zero
        rho_work(1:Mpol,-Ntor:Ntor) = zero

        ! Build rho_work in format (m,n) from inpu
        do ii = 1, mn_rho
          mm = im_rho(ii)
          nn = in_rho(ii) / Nfp

#ifdef DEBUG
          FATAL( bndRep, (mm.lt.1        ).or.(mm.gt.Mpol    ), Error in resolution )
          FATAL( bndRep, (nn/Nfp.lt.-Ntor).or.(nn/Nfp.gt.Ntor), Error in resolution )
#endif

          rho_work(mm,nn) = rhomn( ii )
        enddo

        ! Build RHS
        RHS1(        1:2*Ntor+1 ) = rho_work(1,-Ntor:Ntor)
        RHS1( 2*Ntor+2:3*Ntor+2 ) = bn(0:Ntor)

        ii=0
        do mm=2,Mpol
          RHS2( ii+1:ii+2*Ntor+1 ) = rho_work( mm, -Ntor:Ntor )
          ii = ii+2*Ntor+1
        enddo
  
      end subroutine pack_rhomn_bn
  

      ! This routine unpack RHS1, RHS2 into rhomn, bn
      subroutine unpack_rhomn_bn( rhomn, bn )
        use constants, only: zero
        use inputlist, only: Wmacros, Mpol, Ntor, twoalpha, Nfp
        use fileunits, only: ounit, lunit
        use allglobal, only: myid, im_rho, in_rho, im_field, in_field, mn_rho, mn_field, &
                             MPI_COMM_SPEC, cpus
  
        LOCALS

        REAL, INTENT(OUT) :: rhomn( 1:mn_rho )
        REAL, INTENT(OUT) :: bn(0:Ntor)

        REAL              :: rho_work(1:Mpol,-Ntor:Ntor)

        INTEGER           :: ii, jj, mm, nn
  
        ! Ensure working variables are set to zero
        rho_work(1:Mpol,-Ntor:Ntor) = zero
  
        ! Unpack RHS into a format (m,n)
        rho_work( 1, -Ntor:Ntor ) = RHS1(        1:2*Ntor+1 )
        bn( 0:Ntor )              = RHS1( 2*Ntor+2:3*Ntor+2 )
  
        jj=0
        do mm=2,Mpol
          rho_work( mm, -Ntor:Ntor ) = RHS2( jj+1:jj+2*Ntor+1 )
          jj = jj+2*Ntor+1
        enddo

        ! Build output as a single array
        do ii = 1, mn_rho
          mm = im_rho(ii)
          nn = in_rho(ii) / Nfp

          FATAL( bndRep, (mm.lt.1        ).or.(mm.gt.Mpol    ), Error in resolution )
          FATAL( bndRep, (nn/Nfp.lt.-Ntor).or.(nn/Nfp.gt.Ntor), Error in resolution )

          rhomn( ii ) = rho_work( mm, nn )
        enddo
  
      end subroutine unpack_rhomn_bn
  

      ! This routine pack rmn, zmn into LHS1, LHS2
      subroutine pack_rmn_zmn( rmn, zmn )
        use constants, only: zero
        use inputlist, only: Wmacros, Mpol, Ntor, twoalpha, Nfp
        use fileunits, only: ounit, lunit
        use allglobal, only: myid, im_rho, in_rho, im_field, in_field, mn_rho, mn_field, &
                             MPI_COMM_SPEC, cpus
  
        LOCALS

        REAL, INTENT(IN)  :: rmn(1:mn_field)
        REAL, INTENT(IN)  :: zmn(1:mn_field)

        REAL              :: rmn_work(0:Mpol_field,-Ntor_field:Ntor_field)
        REAL              :: zmn_work(0:Mpol_field,-Ntor_field:Ntor_field)
        INTEGER           :: ii, jj, irz, mm, nn

        ! Ensure work variables are set to zero
        rmn_work(1:Mpol_field,-Ntor_field:Ntor_field) = zero
        zmn_work(1:Mpol_field,-Ntor_field:Ntor_field) = zero

        ! Build array in format (m,n)
        do ii=1,mn_field
          mm=im_field(ii)
          nn=in_field(ii) / Nfp

          if( mm.eq.0 ) cycle ! These are determined by R0c, Z0s.

          rmn_work(mm,nn) = rmn( ii )
          zmn_work(mm,nn) = zmn( ii )
        enddo

        ! Build LHS from input
        LHS1(1:2*Ntor_field+1) = rmn_work(1, -Ntor_field:Ntor_field)
        LHS1(2*Ntor_field+2:2*(2*Ntor_field+1)) = zmn_work(1, -Ntor_field:Ntor_field)
  
        jj=0
        do irz=0,1
          do mm=2,Mpol_field
  
            if( irz==0 ) then
              LHS2( jj+1: jj+2*Ntor_field+1 ) = rmn_work(mm,-Ntor_field:Ntor_field)
            else
              LHS2( jj+1: jj+2*Ntor_field+1 ) = zmn_work(mm,-Ntor_field:Ntor_field)
            endif
  
            jj = jj+2*Ntor_field+1
  
          enddo
        enddo
  
  
      end subroutine pack_rmn_zmn
  
      ! This routine pack LHS1, LHS2 into rmn, zmn.
      subroutine unpack_rmn_zmn( rmn, zmn )
        use constants, only: zero
        use inputlist, only: Wmacros, Mpol, Ntor, twoalpha, Nfp
        use fileunits, only: ounit, lunit
        use allglobal, only: myid, im_rho, in_rho, im_field, in_field, mn_rho, mn_field, &
                             MPI_COMM_SPEC, cpus
  
        LOCALS

        REAL, INTENT(OUT)   :: rmn(1:mn_field)
        REAL, INTENT(OUT)   :: zmn(1:mn_field)

        REAL                :: rmn_work(1:Mpol_field, -Ntor_field:Ntor_field)
        REAL                :: zmn_work(1:Mpol_field, -Ntor_field:Ntor_field)

        INTEGER             :: ii, jj, irz, mm, nn

        LOGICAL             :: Lchangeangle

        ! Ensure work variable are set to zero
        rmn_work(1:Mpol_field,-Ntor_field:Ntor_field) = zero
        zmn_work(1:Mpol_field,-Ntor_field:Ntor_field) = zero

        ! Unpack LHS into work variable in format (m,n)
        rmn_work(1, -Ntor_field:Ntor_field) = LHS1(1:2*Ntor_field+1)
        zmn_work(1, -Ntor_field:Ntor_field) = LHS1(2*Ntor_field+2:2*(2*Ntor_field+1))
  
        jj=0
        do irz=0,1
          do mm=2,Mpol_field
  
            if( irz==0 ) then
              rmn_work(mm,-Ntor_field:Ntor_field) = LHS2( jj+1: jj+2*Ntor_field+1 )
            else
              zmn_work(mm,-Ntor_field:Ntor_field) = LHS2( jj+1: jj+2*Ntor_field+1 )
            endif
  
            jj = jj+2*Ntor_field+1
  
          enddo
        enddo

        ! Build output one-dimensional array
        do ii=1,mn_field
          mm=im_field(ii)
          nn=in_field(ii) / Nfp

          if( mm.eq.0 ) cycle ! These modes already filled by R0c, Z0s in forwardMap

          rmn(ii) = rmn_work(mm,nn)
          zmn(ii) = zmn_work(mm,nn)
        enddo
  
      end subroutine unpack_rmn_zmn
  
      subroutine build_mapping_matrices()
        use constants, only: four
        use inputlist, only: Wmacros, Mpol, Ntor, twoalpha, Nfp
        use fileunits, only: ounit, lunit
        use allglobal, only: myid, im_rho, in_rho, im_field, in_field, mn_rho, mn_rho, &
                             MPI_COMM_SPEC, cpus
  
        LOCALS

        INTEGER :: jjmin, jjmax, nn, jj, ii, mm

        ! Mat1 - mapping matrix for m=1 modes
        do nn = -Ntor_field,Ntor_field
          ! n: toroidal mode number
          ! ii: line index in M1
          ii = Ntor_field + nn + 1
  
          ! Rhomn elements
          ! jj: corresponding column index
          jj = Ntor-nn+1
          if( (jj.le.2*Ntor+1) .and. (jj.gt.0) ) then !Check that indices don't overflow in over elements
            Mat1(                  ii, jj ) = Mat1(                  ii, jj ) + 2 ! Rmn elements
            Mat1( 2*Ntor_field+1 + ii, jj ) = Mat1( 2*Ntor_field+1 + ii, jj ) + 2 ! Zmn  elements
          endif
  
          ! jj: corresponding column index
          jj = Ntor-nn+twoalpha+1
          if( (jj.le.2*Ntor+1) .and. (jj.gt.0) ) then !Check that indices don't overflow in over elements
            Mat1(                  ii, jj ) = Mat1(                  ii, jj ) + 2 ! Rmn elements
            Mat1( 2*Ntor_field+1 + ii, jj ) = Mat1( 2*Ntor_field+1 + ii, jj ) - 2 ! Zmn  elements
          endif
  
          !bn elements
          ! b n
          jj =  nn + 1;
          if( (jj.le.Ntor+1) .and. (jj.gt.0) ) then
              Mat1(               ii, 2*Ntor+1+jj) = Mat1(               ii, 2*Ntor+1+jj) + 1;
              Mat1(2*Ntor_field+1+ii, 2*Ntor+1+jj) = Mat1(2*Ntor_field+1+ii, 2*Ntor+1+jj) + 1;
          endif
  
          ! b -n
          jj = -nn + 1;
          if( (jj.le.Ntor+1) .and. (jj.gt.0) ) then
              Mat1(               ii, 2*Ntor+1+jj) = Mat1(               ii, 2*Ntor+1+jj) + 1;
              Mat1(2*Ntor_field+1+ii, 2*Ntor+1+jj) = Mat1(2*Ntor_field+1+ii, 2*Ntor+1+jj) + 1;
          endif
  
          ! b n- 2 alpha
          jj =  nn - twoalpha + 1;
          if( (jj.le.Ntor+1) .and. (jj.gt.0) ) then
              Mat1(               ii, 2*Ntor+1+jj) = Mat1(               ii, 2*Ntor+1+jj) - 1;
              Mat1(2*Ntor_field+1+ii, 2*Ntor+1+jj) = Mat1(2*Ntor_field+1+ii, 2*Ntor+1+jj) + 1;
          endif
  
          ! b -n + 2 alpha
          jj = -nn + twoalpha + 1;
          if( (jj.le.Ntor+1) .and. (jj.gt.0) ) then
              Mat1(               ii, 2*Ntor+1+jj) = Mat1(               ii, 2*Ntor+1+jj) - 1;
              Mat1(2*Ntor_field+1+ii, 2*Ntor+1+jj) = Mat1(2*Ntor_field+1+ii, 2*Ntor+1+jj) + 1;
          endif
        end do
  
  
  
        ! Mat2 - mapping matrice for modes m>1
        if( Mpol>1 ) then
  
          do mm=2,Mpol_field
            do nn=-Ntor_field,Ntor_field
  
              if( mm.gt.Mpol ) cycle

              ii =    (mm-2)*(2*Ntor_field +1) + Ntor_field + nn + 1;
              jjmin = (mm-2)*(2*Ntor   +1)                       + 1;
              jjmax = (mm-2)*(2*Ntor   +1)    + 2*Ntor           + 1;
  
              jj =    (mm-2)*(2*Ntor   +1)    +   Ntor      - nn + 1;
              ! rho m -n
              if ((jj.ge.jjmin) .and. (jj.le.jjmax)) then
                  Mat2(                             ii, jj ) =  Mat2(                             ii, jj ) + 2; ! For Rmn equation
                  Mat2( (Mpol-1)*(2*Ntor_field+1) + ii, jj ) =  Mat2( (Mpol-1)*(2*Ntor_field+1) + ii, jj ) + 2; ! For Zmn equation
              endif
  
  
              ! rho m -n+2alpha
              jj =    (mm-2)*(2*Ntor   +1)     +  Ntor      - nn + twoalpha + 1;
              if ((jj.ge.jjmin) .and. (jj.le.jjmax)) then
  
                  Mat2(                             ii, jj ) = Mat2(                             ii, jj ) + 2; ! For Rmn equation
                  Mat2( (Mpol-1)*(2*Ntor_field+1) + ii, jj ) = Mat2( (Mpol-1)*(2*Ntor_field+1) + ii, jj ) - 2; ! For Zmn equation
              endif
  
            enddo
          enddo
        endif

        Mat1 = Mat1 / four
        Mat2 = Mat2 / four
  
  
      end subroutine build_mapping_matrices

      subroutine change_angle( Rmn, Zmn )
        use constants, only: zero
        use inputlist, only: Wmacros, Mpol, Ntor, twoalpha, Nfp
        use fileunits, only: ounit, lunit
        use allglobal, only: myid, im_rho, in_rho, im_field, in_field, mn_field, mn_rho, mn_rho, &
                             MPI_COMM_SPEC, cpus
  
        LOCALS

        REAL               :: Rmn(1:mn_field), Zmn(1:mn_field)
        REAL               :: Rmn_tmp(1:mn_field), Zmn_tmp(1:mn_field), mode_r, mode_z
        INTEGER            :: ind, ii, mm, nn

        Rmn_tmp = zero
        Zmn_tmp = zero

        do ii = 1, mn_field
          mm = im_field(ii)
          nn = in_field(ii)

          if( mm.eq.0 ) then
            Rmn_tmp(ii) = Rmn(ii)
            Zmn_tmp(ii) = Zmn(ii)
          else ! ( mm.eq.0 )

            call find_index( mn_field, im_field, in_field, mm, -nn, ind )
            if( ind.eq.-1 ) then ! mode not found, fill with zero
              Rmn_tmp(ii) = zero
              Zmn_tmp(ii) = zero
            else
              Rmn_tmp(ii) =  Rmn( ind )
              Zmn_tmp(ii) = -Zmn( ind ) 
            endif

          endif

        enddo ! ii=1, mn_field

        Rmn(1:mn_field) = Rmn_tmp(1:mn_field)
        Zmn(1:mn_field) = Zmn_tmp(1:mn_field)

      end subroutine change_angle

      subroutine find_index( mn, im_array, in_array, m_search, n_search, ind )
        use inputlist, only: Wmacros
        use fileunits, only: ounit, lunit
        use allglobal, only: mn_field

        LOCALS

        INTEGER, INTENT(IN) :: m_search, n_search, mn, im_array(1:mn), in_array(1:mn)
        INTEGER, INTENT(OUT):: ind

        do ind = 1, mn
          if( im_array(ind).eq.m_search .and. in_array(ind).eq.n_search ) then
            return
          endif
        enddo

        ind = -1
      end subroutine find_index
  
  
  end module bndRep
