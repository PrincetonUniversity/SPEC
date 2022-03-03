module bndRep
    implicit none  
    PUBLIC ! everything is public, excepted stated otherwise.
  
    INTEGER              :: Mpol_field, Ntor_field, Mpol_force, Ntor_force 
    INTEGER              :: ii, jj, mm, nn, jjmin, jjmax, irz ! Loop indices used in subroutines
    INTEGER              :: nel_m1_i, nel_m1_j, nel_m2_i, nel_m2_j ! dimension of matrices
    LOGICAL              :: Lchangeangle
    REAL   , allocatable :: Mat1(:,:), Mat2(:,:), RHS1(:), RHS2(:), LHS1(:), LHS2(:)   !< MAPPING MATRIX, TIMES FOUR
  
    !------- public / private statement ----------
  
    PRIVATE :: Mat1, Mat2
    PRIVATE :: nel_m1_i, nel_m1_j, nel_m2_i, nel_m2_j
    PRIVATE :: ii, jj, mm, nn, jjmin, jjmax, irz
    PRIVATE :: build_mapping_matrices ! only used in forwardMap and backwardMap subroutines.
    PRIVATE :: pack_rhomn_bn, unpack_rhomn_bn, pack_rmn_zmn, unpack_rmn_zmn ! specific for this module
    PRIVATE :: RHS1, RHS2, LHS1, LHS2
    PRIVATE :: Lchangeangle
  
    contains
  
  
  
    ! ------------------------------------------------------------------
    !                     PUBLIC SUBROUTINES
      subroutine initialize_mapping( Langle )
        ! In this subroutine we compute the mapping matrix, and allocate necessary memory
        ! This should only be called once at the beginning of preset.
  
        use inputlist, only: Wmacros, Mpol, Ntor, twoalpha, Nfp, Lboundary
        use fileunits, only: ounit, lunit
        use allglobal, only: myid, mn_field, im_field, in_field, &
                             mn_rho, im_rho, in_rho, &
                             mn_force, im_force, in_force, &
                             MPI_COMM_SPEC, cpus

        LOCALS

        LOGICAL, INTENT(IN):: Langle


        Lchangeangle = Langle

        ! If RZ representation, truncation is the same for all quantities.
        if( Lboundary.eq.0 ) then  
          Mpol_field = Mpol
          Ntor_field = Ntor
    
          Mpol_force  = Mpol
          Ntor_force  = Ntor

        elseif( Lboundary.eq.1 ) then
          Mpol_field = Mpol
          Ntor_field = Ntor + abs(twoalpha)
    
          Mpol_force  = Mpol + 1
          Ntor_force  = Ntor

        endif

  
        mn_field = 1 + Ntor_field +  Mpol_field  * ( 2 *  Ntor_field  + 1 ) ! Fourier resolution of interface geometry & vector potential;
        mn_rho   =                   Mpol        * ( 2 *  Ntor        + 1 )
        mn_force = 1 + Ntor_force +  Mpol_force  * ( 2 *  Ntor_force  + 1 ) ! Fourier resolution of interface geometry & vector potential;
  
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
        
        SALLOCATE( im_field  , (1:mn_field  ), 0 )
        SALLOCATE( in_field  , (1:mn_field  ), 0 )
        SALLOCATE( im_rho    , (1:mn_rho    ), 0 )
        SALLOCATE( in_rho    , (1:mn_rho    ), 0 )
        SALLOCATE( im_force  , (1:mn_force  ), 0 )
        SALLOCATE( in_force  , (1:mn_force  ), 0 )
      
        call gi00ab(  Mpol_field,  Ntor_field, Nfp, mn_field, im_field(1:mn_field), in_field(1:mn_field), .true.  ) ! this sets the im and in mode identification arrays;
        call gi00ab(  Mpol_force,  Ntor_force, Nfp, mn_force, im_force(1:mn_force), in_force(1:mn_force), .true.  ) ! this sets the im and in mode identification arrays;
        call gi00ab(  Mpol      ,  Ntor      , Nfp, mn_rho  , im_rho(  1:mn_rho  ), in_rho(  1:mn_rho  ), .false. ) ! this sets the im and in mode identification arrays;
  
        nel_m1_i = 2*(2*Ntor_field+1)
        nel_m1_j = 3*Ntor+2
        SALLOCATE( Mat1, (1:nel_m1_i,1:nel_m1_j), 0 )
        SALLOCATE( RHS1, (1:nel_m1_j), 0 )
        SALLOCATE( LHS1, (1:nel_m1_i), 0 )
  
        nel_m2_i = 2*(Mpol_field-1)*(2*Ntor_field+1)
        nel_m2_j = (Mpol-1)*(2*Ntor+1)
        SALLOCATE( Mat2, (1:nel_m2_i,1:nel_m2_j), 0 )
        SALLOCATE( RHS2, (1:nel_m2_j), 0 )
        SALLOCATE( LHS2, (1:nel_m2_i), 0 )

        call build_mapping_matrices()
  
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
        INTEGER          :: ind

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
            Zmn(ii) = 0.0
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
        call find_index_field(  1, 0, ind )

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

        REAL                :: Rwork(1:mn_field), Zwork(1:mn_field)

        ! LOCAL VARIABLES
        CHARACTER :: TRANS
        REAL, allocatable :: A(:,:), WORK(:), B(:)
        INTEGER :: NRHS, LDA, LDB, LWORK, INFO

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


      subroutine pack_henneberg( pack, position, bndDofs )

        use constants, only: zero
        use inputlist, only: Wmacros, Mpol, Ntor, twoalpha, Nfp
        use fileunits, only: ounit, lunit
        use allglobal, only: myid, im_rho, in_rho, im_field, in_field, mn_rho, mn_field, &
                             MPI_COMM_SPEC, cpus, Mvol, &
                             NGdof_field, NGdof_bnd

        LOCALS

        ! --------------------------------------------------------------
        ! Local Variables

        CHARACTER, INTENT(IN)  :: pack
        CHARACTER              :: packorunpack ! Either 'RZ_to_H' or 'H_to_RZ'
        LOGICAL                :: LComputeDerivatives, LComputeAxis, work
        INTEGER                :: lvol, idof
        REAL                   :: position( 0:NGdof_field )
        REAL                   :: bndDofs( 0:NGdof_bnd )
        REAL                   :: iRbc(1:mn_field,0:Mvol), iRbs(1:mn_field,0:Mvol), iZbc(1:mn_field,0:Mvol), iZbs(1:mn_field,0:Mvol)
        REAL                   :: irhoc(1:mn_rho,1:Mvol), ibc(0:Ntor,1:Mvol), iR0c(0:Ntor,1:Mvol), iZ0s(1:Ntor,1:Mvol)

        ! --------------------------------------------------------------
        ! Subroutine core

        select case( pack )
        case( 'H' )

          ! First unpack position in Rmn, Zmn
          packorunpack = 'U'
          LComputeDerivatives = .FALSE.
          LComputeAxis = .FALSE.
          call packxi( NGdof_field, position( 0:NGdof_field ), Mvol, mn_field, & 
                       iRbc(1:mn_field,0:Mvol), iZbs(1:mn_field,0:Mvol), iRbs(1:mn_field,0:Mvol), iZbc(1:mn_field,0:Mvol), &
                       packorunpack, LComputeDerivatives, LComputeAxis )

          ! Then map to rhomn, bn, R0n, Z0n
          do lvol=1,Mvol-1
            call backwardMap( iRbc(1:mn_field,lvol), iZbs(1:mn_field,lvol), &
                              irhoc(1:mn_rho,lvol), ibc(0:Ntor, lvol), &
                              iR0c(0:Ntor, lvol), iZ0s(1:Ntor, lvol) )      
          enddo

          ! Finally construct the array bndDofs
          idof = 0
          do lvol=1,Mvol-1
            do jj=1,mn_rho
              idof = idof + 1
#ifdef DEBUG
              FATAL( bndRep, idof.le.0 .or. idof.gt.NGdof_bnd, out of bounds )
#endif

              bndDofs(idof) = irhoc(jj,lvol)

            enddo

            do jj=0,Ntor
              idof = idof+1
              bndDofs(idof) = ibc(jj, lvol)

              idof = idof+1
              bndDofs(idof) = iR0c(jj, lvol)

              if( jj.ne.0 ) then
                idof = idof+1
#ifdef DEBUG
                FATAL( bndRep, idof.le.0 .or. idof.gt.NGdof_bnd, out of bounds )
#endif
                bndDofs(idof) = iZ0s(jj, lvol)
              endif
            enddo
          enddo !lvol


        case( 'R' )

          ! First unpack the array bndDofs
          idof = 0
          do lvol=1,Mvol-1
            do jj=1,mn_rho

              idof = idof+1
#ifdef DEBUG
              FATAL( bndRep, idof.le.0 .or. idof.gt.NGdof_bnd, out of bounds )
#endif

              irhoc(jj,lvol) = bndDofs(idof)
            enddo

            do jj=0,Ntor
              idof = idof+1
              ibc(jj, lvol) = bndDofs(idof)
              
              idof = idof+1
              iR0c(jj, lvol) = bndDofs(idof)

              if( jj.gt.0 ) then 
                idof = idof+1
#ifdef DEBUG
                FATAL( bndRep, idof.le.0 .or. idof.gt.NGdof_bnd, out of bounds )
#endif
                iZ0s(jj, lvol) = bndDofs(idof)
              endif
            enddo
          enddo !lvol

          FATAL( bndRep, idof.ne.NGdof_bnd, incorrect number of dofs. )

          
        ! Then map rhomn, bn, R0n, Z0n to Rmn, Zmn
          do lvol=1,Mvol-1
            call forwardMap( irhoc(1:mn_rho, lvol), ibc(0:Ntor, lvol), &
                            iR0c(0:Ntor, lvol), iZ0s(1:Ntor, lvol), &
                            iRbc(1:mn_field, lvol), iZbs(1:mn_field, lvol) )

#ifdef DEBUG
            if( work ) then
              write(ounit,'("Changing angle in pack_henneberg for lvol=", i3 ,)')
            endif
#endif

          enddo
        

          ! Finally, build position array
          packorunpack = 'P'
          LcomputeDerivatives = .FALSE.
          LComputeAxis = .FALSE.
          call packxi( NGdof_field, position( 0:NGdof_field ), Mvol, mn_field, &
                       iRbc(1:mn_field,0:Mvol), iZbs(1:mn_field,0:Mvol), iRbs(1:mn_field,0:Mvol), iZbc(1:mn_field,0:Mvol), &
                       packorunpack, LComputeDerivatives, LComputeAxis )


        end select 

      end subroutine pack_henneberg
  
  
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
  
        ! Ensure working variables are set to zero
        rho_work(1:Mpol,-Ntor:Ntor) = zero
  
        ! Unpack RHS into a format (m,n)
        rho_work( 1, -Ntor:Ntor ) = RHS1(        1:2*Ntor+1 )
        bn( 0:Ntor )              = RHS1( 2*Ntor+2:3*Ntor+2 )
  
        ii=0
        do mm=2,Mpol
          rho_work( mm, -Ntor:Ntor ) = RHS2( ii+1:ii+2*Ntor+1 )
          ii = ii+2*Ntor+1
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
  
        ii=0
        do irz=0,1
          do mm=2,Mpol_field
  
            if( irz==0 ) then
              LHS2( ii+1: ii+2*Ntor_field+1 ) = rmn_work(mm,-Ntor_field:Ntor_field)
            else
              LHS2( ii+1: ii+2*Ntor_field+1 ) = zmn_work(mm,-Ntor_field:Ntor_field)
            endif
  
            ii = ii+2*Ntor_field+1
  
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

        LOGICAL             :: Lchangeangle

        ! Ensure work variable are set to zero
        rmn_work(1:Mpol_field,-Ntor_field:Ntor_field) = zero
        zmn_work(1:Mpol_field,-Ntor_field:Ntor_field) = zero

        ! Unpack LHS into work variable in format (m,n)
        rmn_work(1, -Ntor_field:Ntor_field) = LHS1(1:2*Ntor_field+1)
        zmn_work(1, -Ntor_field:Ntor_field) = LHS1(2*Ntor_field+2:2*(2*Ntor_field+1))
  
        ii=0
        do irz=0,1
          do mm=2,Mpol_field
  
            if( irz==0 ) then
              rmn_work(mm,-Ntor_field:Ntor_field) = LHS2( ii+1: ii+2*Ntor_field+1 )
            else
              zmn_work(mm,-Ntor_field:Ntor_field) = LHS2( ii+1: ii+2*Ntor_field+1 )
            endif
  
            ii = ii+2*Ntor_field+1
  
          enddo
        enddo

        ! ! Check for changing angle
        ! if( rmn_work(1,0).gt.0 .and. zmn_work(1,0).gt.0 ) then ; Lchangeangle=.true.
        ! else                                                   ; Lchangeangle=.false.
        ! endif

        ! if( Lchangeangle ) then
        !   do ii=1,mn_field
        !     mm=im_field(ii)
        !     nn=in_field(ii) / Nfp

        !     if( mm.eq.0 ) cycle ! These modes already filled in forwardMap

        !     rmn(ii) =  rmn_work(mm,-nn)
        !     zmn(ii) = -zmn_work(mm,-nn)
        !   enddo

        ! else

          ! Build output one-dimensional array
          do ii=1,mn_field
            mm=im_field(ii)
            nn=in_field(ii) / Nfp

            if( mm.eq.0 ) cycle ! These modes already filled by R0c, Z0s in forwardMap

            rmn(ii) = rmn_work(mm,nn)
            zmn(ii) = zmn_work(mm,nn)
          enddo
        ! endif
  
      end subroutine unpack_rmn_zmn
  
      subroutine build_mapping_matrices()
        use inputlist, only: Wmacros, Mpol, Ntor, twoalpha, Nfp
        use fileunits, only: ounit, lunit
        use allglobal, only: myid, im_rho, in_rho, im_field, in_field, mn_rho, mn_rho, &
                             MPI_COMM_SPEC, cpus
  
        LOCALS

        ! Mat1 - mapping matrix for m=1 modes
        do nn = -Ntor_field,Ntor_field
          ! n: toroidal mode number
          ! ii: line index in M1
          ii = Ntor_field + nn + 1
  
          ! Rhomn elements
          ! jj: corresponding column index
          jj = Ntor-nn+1
          if( (jj.le.2*Ntor+1) .and. (jj.gt.0) ) then !Check that indices don't overflow in over elements
            Mat1(              ii, jj ) = Mat1(              ii, jj ) + 2 ! Rmn elements
            Mat1( 2*Ntor_field+1 + ii, jj ) = Mat1( 2*Ntor_field+1 + ii, jj ) + 2 ! Zmn  elements
          endif
  
          ! jj: corresponding column index
          jj = Ntor-nn+twoalpha+1
          if( (jj.le.2*Ntor+1) .and. (jj.gt.0) ) then !Check that indices don't overflow in over elements
            Mat1(              ii, jj ) = Mat1(              ii, jj ) + 2 ! Rmn elements
            Mat1( 2*Ntor_field+1 + ii, jj ) = Mat1( 2*Ntor_field+1 + ii, jj ) - 2 ! Zmn  elements
          endif
  
          !bn elements
          ! b n
          jj =  nn + 1;
          if( (jj.le.Ntor+1) .and. (jj.gt.0) ) then
              Mat1(           ii, 2*Ntor+1+jj) = Mat1(           ii, 2*Ntor+1+jj) + 1;
              Mat1(2*Ntor_field+1+ii, 2*Ntor+1+jj) = Mat1(2*Ntor_field+1+ii, 2*Ntor+1+jj) + 1;
          endif
  
          ! b -n
          jj = -nn + 1;
          if( (jj.le.Ntor+1) .and. (jj.gt.0) ) then
              Mat1(           ii, 2*Ntor+1+jj) = Mat1(           ii, 2*Ntor+1+jj) + 1;
              Mat1(2*Ntor_field+1+ii, 2*Ntor+1+jj) = Mat1(2*Ntor_field+1+ii, 2*Ntor+1+jj) + 1;
          endif
  
          ! b n- 2 alpha
          jj =  nn - twoalpha + 1;
          if( (jj.le.Ntor+1) .and. (jj.gt.0) ) then
              Mat1(           ii, 2*Ntor+1+jj) = Mat1(           ii, 2*Ntor+1+jj) - 1;
              Mat1(2*Ntor_field+1+ii, 2*Ntor+1+jj) = Mat1(2*Ntor_field+1+ii, 2*Ntor+1+jj) + 1;
          endif
  
          ! b -n + 2 alpha
          jj = -nn + twoalpha + 1;
          if( (jj.le.Ntor+1) .and. (jj.gt.0) ) then
              Mat1(           ii, 2*Ntor+1+jj) = Mat1(           ii, 2*Ntor+1+jj) - 1;
              Mat1(2*Ntor_field+1+ii, 2*Ntor+1+jj) = Mat1(2*Ntor_field+1+ii, 2*Ntor+1+jj) + 1;
          endif
        end do
  
  
  
        ! Mat2 - mapping matrice for modes m>1
        if( Mpol>1 ) then
  
          do mm=2,Mpol_field
            do nn=-Ntor_field,Ntor_field
  
              ii =    (mm-2)*(2*Ntor_field +1) + Ntor_field + nn + 1;
              jjmin = (mm-2)*(2*Ntor   +1)               + 1;
              jjmax = (mm-2)*(2*Ntor   +1) + 2*Ntor      + 1;
  
              jj =    (mm-2)*(2*Ntor   +1) + Ntor   - nn + 1;
              ! rho m -n
              if ((jj.ge.jjmin) .and. (jj.le.jjmax)) then
                  Mat2(                         ii, jj ) =  Mat2(                         ii, jj ) + 2; ! For Rmn equation
                  Mat2( (Mpol-1)*(2*Ntor_field+1) + ii, jj ) =  Mat2( (Mpol-1)*(2*Ntor_field+1) + ii, jj ) + 2; ! For Zmn equation
              endif
  
  
              ! rho m -n+2alpha
              jj =    (mm-2)*(2*Ntor   +1) + Ntor    - nn + twoalpha + 1;
              if ((jj.ge.jjmin) .and. (jj.le.jjmax)) then
  
                  Mat2(                         ii, jj ) = Mat2(                         ii, jj ) + 2; ! For Rmn equation
                  Mat2( (Mpol-1)*(2*Ntor_field+1) + ii, jj ) = Mat2( (Mpol-1)*(2*Ntor_field+1) + ii, jj ) - 2; ! For Zmn equation
              endif
  
            enddo
          enddo
        endif

        Mat1 = Mat1 / 4.0
        Mat2 = Mat2 / 4.0
  
  
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
        INTEGER            :: ind

        Rmn_tmp = zero
        Zmn_tmp = zero

        do ii = 1, mn_field
          mm = im_field(ii)
          nn = in_field(ii)

          if( mm.eq.0 ) then
            Rmn_tmp(ii) = Rmn(ii)
            Zmn_tmp(ii) = Zmn(ii)
          else ! ( mm.eq.0 )

            call find_index_field( mm, -nn, ind )
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

      subroutine find_index_field( m_search, n_search, ind )
        use inputlist, only: Wmacros
        use fileunits, only: ounit, lunit
        use allglobal, only: mn_field, im_field, in_field

        LOCALS

        INTEGER, INTENT(IN) :: m_search, n_search
        INTEGER, INTENT(OUT):: ind

        do ind = 1, mn_field
          if( im_field(ind).eq.m_search .and. in_field(ind).eq.n_search ) then
            return
          endif
        enddo

        ind = -1
      end subroutine find_index_field
  
  
  end module bndRep
