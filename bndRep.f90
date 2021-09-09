module bndRep
    implicit none  
    PUBLIC ! everything is public, excepted stated otherwise.
  
    INTEGER              :: MpolRZ, NtorRZ, MpolF, NtorF 
    INTEGER              :: ii, jj, mm, nn, jjmin, jjmax, irz ! Loop indices used in subroutines
    INTEGER              :: nel_m1_i, nel_m1_j, nel_m2_i, nel_m2_j ! dimension of matrices
    INTEGER, allocatable :: Mat1(:,:), Mat2(:,:), RHS1(:), RHS2(:), LHS1(:), LHS2(:)   !< MAPPING MATRIX, TIMES FOUR
  
    !------- public / private statement ----------
  
    PRIVATE :: Mat1, Mat2
    PRIVATE :: nel_m1_i, nel_m1_j, nel_m2_i, nel_m2_j
    PRIVATE :: ii, jj, mm, nn, jjmin, jjmax, irz
    PRIVATE :: build_mapping_matrices ! only used in forwardMap and backwardMap subroutines.
    PRIVATE :: pack_rhomn_bn, unpack_rhomn_bn, pack_rmn_zmn, unpack_rmn_zmn ! specific for this module
    PRIVATE :: RHS1, RHS2, LHS1, LHS2
  
    contains
  
  
  
    ! ------------------------------------------------------------------
    !                     PUBLIC SUBROUTINES
      subroutine initialize_mapping()
        ! In this subroutine we compute the mapping matrix, and allocate necessary memory
        ! This should only be called once at the beginning of preset.
  
        use inputlist, only: Mpol, Ntor, twoalpha, Nfp
        use allglobal, only: myid, im, in, imRZ, inRZ, imf, inf, mn, mnRZ, mnf, &
                             MPI_COMM_SPEC

        LOCALS


        MpolRZ = Mpol
        NtorRZ = Ntor + twoalpha
  
        MpolF  = Mpol + 1
        NtorF  = Ntor
  
        mn   = 1 + Ntor   +  Mpol   * ( 2 *  Ntor   + 1 ) ! Fourier resolution of interface geometry & vector potential;
        mnRZ = 1 + NtorRZ +  MpolRZ * ( 2 *  NtorRZ + 1 ) ! Fourier resolution of interface geometry & vector potential;
        mnf  = 1 + NtorF  +  MpolF  * ( 2 *  NtorF  + 1 ) ! Fourier resolution of interface geometry & vector potential;
  
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
        
        SALLOCATE( im  , (1:mn  ), 0 )
        SALLOCATE( in  , (1:mn  ), 0 )
        SALLOCATE( imRZ, (1:mnRZ), 0 )
        SALLOCATE( inRZ, (1:mnRZ), 0 )
        SALLOCATE( imf , (1:mnf ), 0 )
        SALLOCATE( inf , (1:mnf ), 0 )
      
        call gi00ab(  Mpol,  Ntor, Nfp, mn  , im(1:mn  ), in(1:mn  ) ) ! this sets the im and in mode identification arrays;
        call gi00ab(  Mpol,  Ntor, Nfp, mnRZ, imRZ(1:mnRZ), inRZ(1:mnRZ) ) ! this sets the im and in mode identification arrays;
        call gi00ab(  Mpol,  Ntor, Nfp, mnf , imf(1:mnf ), inf(1:mnf ) ) ! this sets the im and in mode identification arrays;
  
        nel_m1_i = 2*(2*NtorRZ+1)
        nel_m1_j = 3*Ntor+2
        SALLOCATE( Mat1, (1:nel_m1_i,1:nel_m1_j), 0 )
        SALLOCATE( RHS1, (1:nel_m1_j), 0 )
        SALLOCATE( LHS1, (1:nel_m1_i), 0 )
  
        nel_m2_i = 2*(MpolRZ-1)*(2*NtorRZ+1)
        nel_m2_j = (Mpol-1)*(2*Ntor+1)
        SALLOCATE( Mat2, (1:nel_m2_i,1:nel_m2_j), 0 )
        SALLOCATE( RHS2, (1:nel_m2_j), 0 )
        SALLOCATE( LHS2, (1:nel_m2_i), 0 )
  
      end subroutine initialize_mapping
  
  
  
      subroutine forwardMap( rhomn, bn, R0c, Z0s, Rmn, Zmn )
  
        use inputlist, only: Mpol, Ntor, twoalpha
        use allglobal, only: myid, im, in, imRZ, inRZ, mn, mnRZ, &
                             MPI_COMM_SPEC

        LOCALS

        ! INPUTS
        REAL, intent(in) :: R0c(0:Ntor), Z0s(1:Ntor)
        REAL, intent(out):: Rmn(0:MpolRZ, -NtorRZ:NtorRZ), Zmn(0:MpolRZ, -NtorRZ:NtorRZ)
        REAL, intent(in) :: rhomn(0:Mpol,-Ntor:Ntor), bn(0:Ntor)

  
        ! ---------------------------------------------------------------
  
        ! m=0 modes
        Rmn(0,0:Ntor) = R0c(0:Ntor)
        Zmn(0,0)      = 0.0
        Zmn(0,1:Ntor) =-Z0s(1:Ntor)
  
        ! m=1 modes
        call pack_rhomn_bn( rhomn, bn )
        LHS1 = MATMUL( Mat1, RHS1 )
  
        ! m>1 modes
        if( MpolRZ.gt.1 ) then
          LHS2 = MATMUL( Mat2, RHS2 )
        endif
  
        call unpack_rmn_zmn( Rmn, Zmn )
  
      end subroutine forwardMap
  
  
  
      subroutine backwardMap( Rmn, Zmn, rhomn, bn, R0c, Z0s )
  
        use constants, only: zero
        use inputlist, only: Mpol, Ntor, twoalpha
        use allglobal, only: myid, im, in, imRZ, inRZ, mn, mnRZ, &
                             MPI_COMM_SPEC

        LOCALS

        ! INPUTS
        REAL, intent(in)  :: Rmn(0:MpolRZ, -NtorRZ:NtorRZ), Zmn(0:MpolRZ, -NtorRZ:NtorRZ)
        REAL, intent(out) :: rhomn(0:Mpol,-Ntor:Ntor), bn(0:Ntor)
        REAL, intent(out) :: R0c(0:Ntor), Z0s(1:Ntor)
  
        ! LOCAL VARIABLES
        CHARACTER :: TRANS
        REAL, allocatable :: A(:,:), WORK(:), B(:)
        INTEGER :: NRHS, LDA, LDB, LWORK, INFO

  
        ! m=0 modes
        R0c(0:Ntor) = Rmn(0,0:Ntor)
        Z0s(1:Ntor) = Zmn(0,1:Ntor)
  
        ! m>0 modes
        SALLOCATE( A, (1:nel_m1_i,1:nel_m1_j), zero )  
        SALLOCATE( B, (1:nel_m1_i), zero )    
  
        A = Mat1
        call pack_rhomn_bn( rhomn, bn )
        TRANS = 'N'
        NRHS = 1
        LDA = nel_m1_i
        LDB = nel_m1_i
  
        ! Work query to get optimal block size
        LWORK = -1
        SALLOCATE( WORK, (1:1), zero )
        call DGELS( TRANS, nel_m1_i, nel_m1_j, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )
  
        select case( LWORK )
        case( 0 )
          LWORK = WORK(1)
          DALLOCATE( WORK )
  
        case( :-1)
          FATAL( global, .true., Illegal value in DGELS )
  
        case(1: )
          FATAL( global, .true., Rank zero for backward mapping )
  
        end select
  
  
        ! Now solve linear system
        call pack_rmn_zmn( rmn, zmn )

        SALLOCATE( WORK, (1:LWORK), zero )
        LDA = nel_m1_i
        LDB = nel_m1_i
        A = Mat1
        B(1:nel_m1_i) = LHS1(1:nel_m1_i)
        call DGELS( TRANS, nel_m1_i, nel_m1_j, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )
        
        select case( LWORK )
        case( 0 )
          RHS1(1:nel_m1_j) = B(1:nel_m1_j)
  
        case( :-1)
          FATAL( global, .true., Illegal value in DGELS )
  
        case(1: )
          FATAL( global, .true., Rank zero for backward mapping )
  
        end select
        

        
        ! Now solve second system
        DALLOCATE( A )
        DALLOCATE( B )
        SALLOCATE( A, (1:nel_m2_i,1:nel_m2_j), zero )
        SALLOCATE( B, (1:nel_m2_i), zero )
        A = Mat2
        LDA = nel_m2_i
        LDB = nel_m2_i
        B(1:nel_m2_i) = LHS2(1:nel_m2_i)
        call DGELS( TRANS, nel_m2_i, nel_m2_j, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )
        select case( LWORK )
        case( 0 )
          RHS2(1:nel_m2_j) = B(1:nel_m2_j)
  
        case( :-1)
          FATAL( global, .true., Illegal value in DGELS )
  
        case(1: )
          FATAL( global, .true., Rank zero for backward mapping )
  
        end select

        ! Unpack RHS1, RHS2 and store in rhomn, bn.
        call unpack_rhomn_bn( rhomn, bn )


        DALLOCATE( WORK )
        DALLOCATE( A )
        DALLOCATE( B )
  
  
      end subroutine backwardMap
  
  
    ! ------------------------------------------------------------------
    !                     PRIVATE SUBROUTINES
  
      subroutine pack_rhomn_bn( rhomn, bn )
        use inputlist, only: Mpol, Ntor, twoalpha
        use allglobal, only: myid, im, in, imRZ, inRZ, mn, mnRZ, &
                             MPI_COMM_SPEC
  
        LOCALS

        REAL, INTENT(IN)  :: rhomn(0:Mpol, -Ntor:Ntor)
        REAL, INTENT(IN)  :: bn(0:Ntor)

  
        RHS1(        1:2*Ntor+1 ) = rhomn(1,-Ntor:Ntor)
        RHS1( 2*Ntor+2:3*Ntor+3 ) = bn(0:Ntor)
  
        ii=0
        do mm=2,Mpol
          RHS2( ii+1:ii+2*Ntor+1 ) = rhomn( mm, -Ntor:Ntor )
          ii = ii+2*Ntor+1
        enddo
  
      end subroutine pack_rhomn_bn
  
      subroutine unpack_rhomn_bn( rhomn, bn )
        use inputlist, only: Mpol, Ntor, twoalpha
        use allglobal, only: myid, im, in, imRZ, inRZ, mn, mnRZ, &
                             MPI_COMM_SPEC
  
        LOCALS

        REAL, INTENT(OUT) :: rhomn(0:Mpol, -Ntor:Ntor)
        REAL, INTENT(OUT) :: bn(0:Ntor)
  
  
        rhomn( 1, -Ntor:Ntor ) = RHS1(        1:2*Ntor+1 )
        bn( 0:Ntor )           = RHS1( 2*Ntor+2:3*Ntor+3 )
  
        ii=0
        do mm=2,Mpol
          rhomn( mm, -Ntor:Ntor ) = RHS2( ii+1:ii+2*Ntor+1 )
          ii = ii+2*Ntor+1
        enddo
  
      end subroutine unpack_rhomn_bn
  
      subroutine pack_rmn_zmn( rmn, zmn )
        use inputlist, only: Mpol, Ntor, twoalpha
        use allglobal, only: myid, im, in, imRZ, inRZ, mn, mnRZ, &
                             MPI_COMM_SPEC
  
        LOCALS

        REAL, INTENT(IN)  :: rmn(0:MpolRZ,-NtorRZ:NtorRZ)
        REAL, INTENT(IN)  :: zmn(0:MpolRZ,-NtorRZ:NtorRZ)
  

        LHS1(1:2*NtorRZ+1) = rmn(1, -NtorRZ:NtorRZ)
        LHS1(2*NtorRZ+2:2*(2*NtorRZ+1)) = zmn(1, -NtorRZ:NtorRZ)
  
        ii=0
        do irz=0,1
          do mm=2,MpolRZ
  
            if( irz==1 ) then
              LHS2( ii+1: ii+2*NtorRZ+1 ) = rmn(mm,-NtorRZ:NtorRZ)
            else
              LHS2( ii+1: ii+2*NtorRZ+1 ) = zmn(mm,-NtorRZ:NtorRZ)
            endif
  
            ii = ii+2*NtorRZ+1
  
          enddo
        enddo
  
  
      end subroutine pack_rmn_zmn
  
      subroutine unpack_rmn_zmn( rmn, zmn )
        use inputlist, only: Mpol, Ntor, twoalpha
        use allglobal, only: myid, im, in, imRZ, inRZ, mn, mnRZ, &
                             MPI_COMM_SPEC
  
        LOCALS

        REAL, INTENT(OUT) :: rmn(0:MpolRZ,-NtorRZ:NtorRZ)
        REAL, INTENT(OUT) :: zmn(0:MpolRZ,-NtorRZ:NtorRZ)
  

        rmn(1, -NtorRZ:NtorRZ) = LHS1(1:2*NtorRZ+1)
        zmn(1, -NtorRZ:NtorRZ) = LHS1(2*NtorRZ+2:2*(2*NtorRZ+1))
  
        ii=0
        do irz=0,1
          do mm=2,MpolRZ
  
            if( irz==1 ) then
              rmn(mm,-NtorRZ:NtorRZ) = LHS2( ii+1: ii+2*NtorRZ+1 )
            else
              zmn(mm,-NtorRZ:NtorRZ) = LHS2( ii+1: ii+2*NtorRZ+1 )
            endif
  
            ii = ii+2*NtorRZ+1
  
          enddo
        enddo
  
  
  
      end subroutine unpack_rmn_zmn
  
      subroutine build_mapping_matrices()
        use inputlist, only: Mpol, Ntor, twoalpha
        use allglobal, only: myid, im, in, imRZ, inRZ, mn, mnRZ, &
                             MPI_COMM_SPEC
  
        LOCALS

        ! Mat1 - mapping matrix for m=1 modes
        do nn = -NtorRZ,NtorRZ
          ! n: toroidal mode number
          ! ii: line index in M1
          ii = NtorRZ + nn + 1
  
          ! Rhomn elements
          ! jj: corresponding column index
          jj = Ntor-nn+1
          if( (jj.le.2*Ntor+1) .and. (jj.gt.0) ) then !Check that indices don't overflow in over elements
            Mat1(              ii, jj ) = Mat1(              ii, jj ) + 2 ! Rmn elements
            Mat1( 2*NtorRZ+1 + ii, jj ) = Mat1( 2*NtorRZ+1 + ii, jj ) + 2 ! Zmn  elements
          endif
  
          ! jj: corresponding column index
          jj = Ntor-nn+twoalpha+1
          if( (jj.le.2*Ntor+1) .and. (jj.gt.0) ) then !Check that indices don't overflow in over elements
            Mat1(              ii, jj ) = Mat1(              ii, jj ) + 2 ! Rmn elements
            Mat1( 2*NtorRZ+1 + ii, jj ) = Mat1( 2*NtorRZ+1 + ii, jj ) - 2 ! Zmn  elements
          endif
  
          !bn elements
          ! b n
          jj =  nn + 1;
          if( (jj.le.Ntor+1) .and. (jj.gt.0) ) then
              Mat1(           ii, 2*Ntor+1+jj) = Mat1(           ii, 2*Ntor+1+jj) + 1;
              Mat1(2*NtorRZ+1+ii, 2*Ntor+1+jj) = Mat1(2*NtorRZ+1+ii, 2*Ntor+1+jj) + 1;
          endif
  
          ! b -n
          jj = -nn + 1;
          if( (jj.le.Ntor+1) .and. (jj.gt.0) ) then
              Mat1(           ii, 2*Ntor+1+jj) = Mat1(           ii, 2*Ntor+1+jj) + 1;
              Mat1(2*NtorRZ+1+ii, 2*Ntor+1+jj) = Mat1(2*NtorRZ+1+ii, 2*Ntor+1+jj) + 1;
          endif
  
          ! b n- 2 alpha
          jj =  nn - twoalpha + 1;
          if( (jj.le.Ntor+1) .and. (jj.gt.0) ) then
              Mat1(           ii, 2*Ntor+1+jj) = Mat1(           ii, 2*Ntor+1+jj) - 1;
              Mat1(2*NtorRZ+1+ii, 2*Ntor+1+jj) = Mat1(2*NtorRZ+1+ii, 2*Ntor+1+jj) + 1;
          endif
  
          ! b -n + 2 alpha
          jj = -nn + twoalpha + 1;
          if( (jj.le.Ntor+1) .and. (jj.gt.0) ) then
              Mat1(           ii, 2*Ntor+1+jj) = Mat1(           ii, 2*Ntor+1+jj) - 1;
              Mat1(2*NtorRZ+1+ii, 2*Ntor+1+jj) = Mat1(2*NtorRZ+1+ii, 2*Ntor+1+jj) + 1;
          endif
        end do
  
  
  
        ! Mat2 - mapping matrice for modes m>1
        if( Mpol>1 ) then
  
          do mm=2,MpolRZ
            do nn=-NtorRZ,NtorRZ
  
              ii =    (mm-2)*(2*NtorRZ +1) + NtorRZ + nn + 1;
              jjmin = (mm-2)*(2*Ntor   +1)               + 1;
              jjmax = (mm-2)*(2*Ntor   +1) + 2*Ntor      + 1;
  
              jj =    (mm-2)*(2*Ntor   +1) + Ntor   - nn + 1;
              ! rho m -n
              if ((jj.ge.jjmin) .and. (jj.le.jjmax)) then
                  Mat2(                         ii, jj ) =  Mat2(                         ii, jj ) + 2; ! For Rmn equation
                  Mat2( (Mpol-1)*(2*NtorRZ+1) + ii, jj ) =  Mat2( (Mpol-1)*(2*NtorRZ+1) + ii, jj ) + 2; ! For Zmn equation
              endif
  
  
              ! rho m -n+2alpha
              jj =    (mm-2)*(2*Ntor   +1) + Ntor    - nn + twoalpha + 1;
              if ((jj.ge.jjmin) .and. (jj.le.jjmax)) then
  
                  Mat2(                         ii, jj ) = Mat2(                         ii, jj ) + 2; ! For Rmn equation
                  Mat2( (Mpol-1)*(2*NtorRZ+1) + ii, jj ) = Mat2( (Mpol-1)*(2*NtorRZ+1) + ii, jj ) - 2; ! For Zmn equation
              endif
  
            enddo
          enddo
        endif
  
  
      end subroutine build_mapping_matrices
  
  
  end module bndRep