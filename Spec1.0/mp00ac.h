!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Solves Beltrami linear system for given helicity multiplier and poloidal flux and returns an error function.

!latex \subsection{Plasma region}

!latex \item If \verb+Lplasmaregion+, then the ``Beltrami'' matrix (see 1.18 of manual.pdf) is constructed from \verb+dMA+ and \verb+dMD+.
!latex \item If the user has selected \verb+Lposdef+, then the user has assumed that the Beltrami matrix is positive-definite, 
!latex       and so the NAG routine \verb+F04ABF+ is used to solve the linear system, i.e. to solve for the vector-potential in ``packed'' format.
!latex \item If the user has selected \verb+.not.Lposdef+, then the user has assumed that the Beltrami matrix is {\em not} positive-definite, 
!latex       and so the NAG routine \verb+F04AEF+ is used.
!latex \item If it is known that the Beltrami matrix is positive definite, then it is usually more computationally efficient to exploit this,
!latex       so the {\em fastest} option is to choose \verb+Lposdef=T+; 
!latex       however, there is no a-priori reason why the Beltrami matrix is positive definite, so the {\em safest} option is to choose \verb+Lposdef=F+.

!latex \subsection{Vacuum region}

!latex \item If \verb+Lplasmaregion+, then the ``Beltrami'' matrix is constructed from \verb+dMA+.
!latex \item The Beltrami matrix for the vacuum region is identically the Laplacian, and is always positive-definite, 
!latex       so the NAG routine \verb+F04ABF+ is always used to solve the linear system, i.e. to solve for the scalar-potential in ``packed-format''.
!latex \item The positive-definitess is guaranteed because the functional $\int B^2 dv$ is guaranteed to have a minimum. [Is this correct?]

!latex \subsection{Error messages}

!latex \item The \verb+F04ABF+ routine may fail if the provided matrix is not positive definite or ill-conditioned, 
!latex       and such error messages are given as screen output.
!latex \item The \verb+F04AEF+ routine may fail if the provided matrix is singular              or ill-conditioned. 

!latex \subsection{Energy and helicity integrals}

!latex \item It is most convenient to calculate the energy and helicity integrals with the vector potential in the packed format, 
!latex       as these integrals reduce to vector-matrix-vector products.

!latex \subsection{``Unpacking''}

!latex \item The routine \verb+up00aa+ is then called to ``unpack'' the linear solution into the more easily identifiable arrays.

!latex \subsection{Returned values}

!latex \item This routine may be called iteratively to find the particular $(\Delta \psi_{p},\mu)$ that satisfies the required constraints, 
!latex       e.g. for \verb+Lconstraint=1+, the user wishes to enforce the rotational-transform constraint on the adjacent interfaces, 
!latex       which is calculated by calling \verb+tr00ab+.

!latex \item Because this routine is called by NAG the input/ouput arguments are constrained, and \verb+mp00ac+ returns either a function,
!latex       which is equal to zero when the appropriate constraints are satisfied, or the derivative of the same function 
!latex       with respect to $\Delta \psi_p$ and/or $\mu$, as the case may be.
!latex \item Note that the derivatives of the function are determined by matrix perturbation methods.
!latex       If the derivatives are required, then the derivatives of the Beltrami matrix must be provided in the \verb+dMA+ and \verb+dMD+ arrays.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine mp00ac( Nmudp, Xmudp, Fmudp, Dmudp, ldfjac, iflag ) ! argument list is fixed by NAG; ma02aa calls mp00ac through C05PCF;

! if iflag.eq.0 : Xmudp and Fmudp are available for PRINTING ; 28 Jan 13; Fmudp MUST NOT BE CHANGED; Dmudp MUST NOT BE CHANGED;
! if iflag.eq.1 :           Fmudp is to be          UPDATED  ; 28 Jan 13;                          ; Dmudp MUST NOT BE CHANGED;
! if iflag.eq.2 :           Dmudp is to be          UPDATED  ; 28 Jan 13; Fmudp MUST NOT BE CHANGED;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half, one, two, goldenmean
  
  use numerical, only : machprec, sqrtmachprec, small
  
  use fileunits, only : ounit
  
  use inputlist, only : Wmacros, Wmp00ac, Wtr00ab, Wma02aa, Nvol, &
                        mu, helicity, iota, oita, curtor, curpol, &
                        Lposdef, &
                        Lconstraint, mupftol
  
  use cputiming, only : Tmp00ac
  
  use allglobal, only : myid, ncpu, cpus, ivol, &
                        YESstellsym, NOTstellsym, &
                        Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, &
                        Mvol, mn, im, in, mns, &
                        Nt, Nz, & ! only required to pass through as arguments to tr00ab; 23 Apr 13;
                        Nmagneticdof, &
                        dMA, dMB, dMC, dMD, dME, dMF, &
                        solution, &
                        dtflux, dpflux, &
                        diota, & ! this is global so that rotational-transform estimate provided by tr00ab can be compared to field-line tracing in pp00aa;
                        lBBintegral, lABintegral, &
                        xoffset, &
                        ImagneticOK
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in)  :: Nmudp, ldfjac
  REAL   , intent(in)  :: Xmudp(1:Nmudp)
  REAL                 :: Fmudp(1:Nmudp), Dmudp(1:ldfjac,1:Nmudp)
  INTEGER              :: iflag ! indicates whether (i) iflag=1: ``function'' values are required; or (ii) iflag=2: ``derivative'' values are required;
  

  INTEGER              :: lvol, NN, MM, ideriv, IA, IB, IC, IBB, IAA, lmns, if04abf(0:1), if04aef(0:1), ii
  
  REAL                 :: lmu, dpf, dpsi(1:2), lpsi(1:2), Itor, Gpol, lcpu
  
  CHARACTER            :: packorunpack
  
  REAL   , allocatable :: matrix(:,:), rhs(:,:)

  REAL   , allocatable :: rwork(:), residual(:,:), LU(:,:)
  
  BEGIN(mp00ac)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  lvol = ivol ! RECALL THAT vvol is global; 24 Jan 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  FATALMESS(mp00ac, iflag.ne.1 .and. iflag.ne.2, invalid iflag) ! 25 Jan 13;
  FATALMESS(mp00ac, lvol.lt.1 .or. lvol.gt.Mvol, invalid lvol) ! 25 Jan 13;  
#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Lplasmaregion ) then ! inside plasma; 19 Apr 13;
   
   ;                     ; lmu  = Xmudp(1) - xoffset ! note construction of Xmudp in ma02aa; included so NAG normalizes error to one + mu;
   if( Nmudp.eq.2 ) then ; dpf  = Xmudp(2) - xoffset ! enclosed poloidal flux provided as argument; 28 Jan 13;
   else                  ; dpf  = dpflux(lvol)       ! value provided through global variable; 28 Jan 13;
   endif
   
   dpsi(1:2) = (/ dtflux(lvol), dpf /) ; lpsi(1:2) = (/ zero, one /) ! enclosed toroidal fluxes; 28 Jan 13;
   
  else ! inside vacuum region; 19 Apr 13;
   
   ;                     ; Itor = Xmudp(1) - xoffset ! shall iterate on plasma current in order to match transform; 22 Apr 13;
   ;                     ; Gpol = curpol

  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  diota(0:1,-1:2,lvol) = zero ! rotational-transform, and its derivatives with respect to lmu and dpf, or toroidal current, on the inner and outer interface;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  NN = Nmagneticdof(lvol) ! shorthand;
  
  RALLOCATE(matrix,(1:NN,1:NN))

  RALLOCATE(rhs,(1:NN,0:2))

  solution(1:NN,-1:2) = zero ! this is a global array allocated in fc02aa; 20 Jun 14;
  
  RALLOCATE(rwork,(1:NN))
  RALLOCATE(residual,(1:NN,0:2))

  if( .not. Lposdef ) then
  RALLOCATE(LU,(1:NN,1:NN))
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if04abf(0:1) = 0 ; if04aef(0:1) = 0 ! error flags;  4 Feb 13;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! ideriv labels derivative as follows:
! ideriv = 0 : compute Beltrami field; ideriv = 1 : compute d Beltrami field / d \mu ; ideriv = 2 : compute d Beltrami field / d \delta \psi_p ;
  

  do ideriv = 0, 1 ! loop over derivatives with respect to helicity multipler and poloidal flux, or the toroidal plasma current;
      

   if( iflag.eq.1 .and. ideriv.gt.0 ) cycle ! only need to return function; recall the derivative estimate requires function evaluation;
   

   if( Lplasmaregion ) then

!#ifdef DEBUG
!     if( Wmp00ac ) then
!      do ii = 1, NN
!      write(ounit,'("mp00ac : ", 10x ," : dMA   =",999es08.01)') dMA(ii,1:NN)
!      enddo
!     !pause
!      do ii = 1, NN
!      write(ounit,'("mp00ac : ", 10x ," : dMD   =",999es08.01)') dMD(ii,1:NN)
!      enddo
!     !pause
!      write(ounit,'("mp00ac : ", 10x ," : lmu   =",999es08.01)') lmu
!     !pause
!     endif
!#endif
    
    matrix(1:NN,1:NN) = dMA(1:NN,1:NN) - lmu * dMD(1:NN,1:NN) ! ma01ag must have been called prior to mp00ac; 28 Jan 13;
    
    select case( ideriv )
    case( 0 )    ; rhs(1:NN,0) = - matmul(  dMB(1:NN,1:2) - lmu  * dME(1:NN,1:2), dpsi(1:2) )
    case( 1 )    ; rhs(1:NN,1) = - matmul(                - one  * dME(1:NN,1:2), dpsi(1:2) ) - matmul( - one  * dMD(1:NN,1:NN), solution(1:NN,0) )
     ;           ; rhs(1:NN,2) = - matmul(  dMB(1:NN,1:2) - lmu  * dME(1:NN,1:2), lpsi(1:2) )
    case default ; FATALMESS(mp00ac, .true., invalid ideriv)
    end select

   else ! inside vacuum region; 19 Apr 13;

    matrix(1:NN,1:NN) = dMA(1:NN,1:NN) ! this was constructed in va00aa, which is called by fc02aa; 21 Apr 13;
    
    select case( ideriv )
    case( 0 )    ; rhs(1:NN,0) = dMB(1:NN,0) - Itor * dMB(1:NN,1) - Gpol * dMB(1:NN,2)
    case( 1 )    ; rhs(1:NN,1) =             -        dMB(1:NN,1)                      ! perturbed matrix is trivial;
     ;           ; rhs(1:NN,2) =                                  -        dMB(1:NN,2) ! not actually required; 21 Apr 13;
    case default ; FATALMESS(mp00ac, .true., invalid ideriv)
    end select
    
   endif ! end of if( Lplasmaregion ) ; 19 Apr 13;
   

   IA = NN ; IB = NN ; IC = NN ; IBB = NN ; IAA = NN
   
   
   lcpu = GETTIME

!#ifdef DEBUG
!    if( Wmp00ac ) then
!     do ii = 1, NN
!      write(ounit,'("mp00ac : matrix=",999es08.01)') matrix(ii,1:NN)
!     enddo
!    endif
!#endif
   
   if( Lposdef .or. Lvacuumregion ) then ! vacuum calculation is always positive definite; 17 Apr 13;
    
!#ifdef DEBUG
!    if( Wmp00ac ) then
!     write(ounit,'("mp00ac : ", 10x ," : calling F04ABF ;")')
!    !pause
!    endif
!#endif
    
    select case( ideriv )
     
    case( 0 ) ! matches select case( ideriv ) ; 22 Apr 13;
     
     if04abf(ideriv) = 1 ; MM = 1
     
    !call F04ABF( A                , IA, B            , IB, N , M , C                 , IC, WKSPCE     , BB                 , IBB, IFAIL           )
     call F04ABF( matrix(1:IA,1:NN), IA, rhs(1:IB,0:0), IB, NN, MM, solution(1:IC,0:0), IC, rwork(1:NN), residual(1:IBB,0:0), IBB, if04abf(ideriv) )
     
    case( 1 ) ! matches select case( ideriv ) ; 22 Apr 13;
     
     if04abf(ideriv) = 1 ; MM = 2

     call F04ABF( matrix(1:IA,1:NN), IA, rhs(1:NN,1:2), IB, NN, MM, solution(1:NN,1:2), IC, rwork(1:NN), residual(1:IBB,1:2), IBB, if04abf(ideriv) )
     
    end select ! matches select case( ideriv ) ; 22 Apr 13;
    
    cput = GETTIME

    select case( if04abf(ideriv) ) !                                                                           123456789012345678
    case(  0  )  ; if( Wmp00ac ) write(ounit,1010) cput-cpus, myid, lvol, ideriv, "if04abf", if04abf(ideriv), "success ;         ", cput-lcpu
    case(  1  )  ;               write(ounit,1010) cput-cpus, myid, lvol, ideriv, "if04abf", if04abf(ideriv), "not +ve definite ;"            ! pause
    case(  2  )  ;               write(ounit,1010) cput-cpus, myid, lvol, ideriv, "if04abf", if04abf(ideriv), "ill conditioned ; "            ! pause
    case(  3  )  ;               write(ounit,1010) cput-cpus, myid, lvol, ideriv, "if04abf", if04abf(ideriv), "input error ;     "            ! pause
    case default ;               write(ounit,1010) cput-cpus, myid, lvol, ideriv, "if04abf", if04abf(ideriv), "invalid if04abf ; "            ! pause
    end select
    
   !FATALMESS(mp00ac, .true., de-bugging: under-construction)
    
   else ! matches if( Lposdef .or. Lvacuumregion ) ; 

!#ifdef DEBUG
!    if( Wmp00ac ) then
!     write(ounit,'("mp00ac : ", 10x ," : calling F04AEF ;")')
!    !pause
!    endif
!#endif

    select case( ideriv ) ! still inside do ideriv = 0, 1; 21 Apr 13;

    case( 0 ) ! matches select case( ideriv ) ; 22 Apr 13;

     if04aef(ideriv) = 1 ; MM = 1
     call F04AEF( matrix(1:NN,1:NN), IA, rhs(1:NN,  0), IB, NN, MM, solution(1:NN,  0), &
  IC, rwork(1:NN), LU(1:NN,1:NN), IAA, residual(1:NN,  0), IBB, if04aef(ideriv) )

    case( 1 ) ! matches select case( ideriv ) ; 22 Apr 13;

     if04aef(ideriv) = 1 ; MM = 2
     call F04AEF( matrix(1:NN,1:NN), IA, rhs(1:NN,1:2), IB, NN, MM, solution(1:NN,1:2), &
  IC, rwork(1:NN), LU(1:NN,1:NN), IAA, residual(1:NN,1:2), IBB, if04aef(ideriv) )

    end select ! matches select case( ideriv ) ; 22 Apr 13;
    
    cput = GETTIME

    select case( if04aef(ideriv) )                                                                            !123456789012345678
    case(  0  )  ; if( Wmp00ac ) write(ounit,1010) cput-cpus, myid, lvol, ideriv, "if04aef", if04aef(ideriv), "success ;         ", cput-lcpu
    case(  1  )  ;               write(ounit,1010) cput-cpus, myid, lvol, ideriv, "if04aef", if04aef(ideriv), "singular ;        "            ! pause
    case(  2  )  ;               write(ounit,1010) cput-cpus, myid, lvol, ideriv, "if04aef", if04aef(ideriv), "ill conditioned ; "            ! pause
    case(  3  )  ;               write(ounit,1010) cput-cpus, myid, lvol, ideriv, "if04aef", if04aef(ideriv), "input error ;     "            ! pause
    case default ;               write(ounit,1010) cput-cpus, myid, lvol, ideriv, "if04aef", if04aef(ideriv), "invalid if04aef ; "            ! pause
    end select


   endif ! end of if( Lposdef .or. Lvacuumregion ) ; 31 Jan 13;
   
1010 format("mp00ac : ",f10.2," : myid=",i3," ; lvol=",i3," ; ideriv="i2" ; "a7"=",i3," ; "a18,:" time=",f10.2," ;")
   
  enddo ! end of do ideriv; 25 Jan 13;
  

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! can compute the energy and helicity integrals; easiest to do this with solution in packed format; 20 Feb 13;
  
  if( Lplasmaregion ) then
   
   lBBintegral(lvol) = half * sum( solution(1:NN,0) * matmul( dMA(1:NN,1:NN), solution(1:NN,0) ) ) & 
                     +        sum( solution(1:NN,0) * matmul( dMB(1:NN,1: 2),     dpsi(1: 2  ) ) ) &
                     + half * sum(     dpsi(1: 2  ) * matmul( dMC(1: 2,1: 2),     dpsi(1: 2  ) ) )
  
   lABintegral(lvol) = half * sum( solution(1:NN,0) * matmul( dMD(1:NN,1:NN), solution(1:NN,0) ) ) & 
                     +        sum( solution(1:NN,0) * matmul( dME(1:NN,1: 2),     dpsi(1: 2  ) ) ) &
                     + half * sum(     dpsi(1: 2  ) * matmul( dMF(1: 2,1: 2),     dpsi(1: 2  ) ) )

#ifdef DEBUG
   if( Wmp00ac ) then
    cput = GETTIME
    write(ounit,'("mp00ac : ",f10.2," : myid=",i3," ; lvol=",i3," ; BB="es23.15" ; AB="es23.15" ; ")') &
 cput-cpus, myid, lvol, lBBintegral(lvol), lABintegral(lvol)
   endif
#endif
   
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  do ideriv = 0, 2
   
   if( iflag.eq.1 .and. ideriv.gt.0 ) cycle
   
   packorunpack = 'U'
   WCALL(mp00ac, up00aa,( packorunpack, lvol, NN, solution(1:NN,ideriv), dpsi(1:2), ideriv )) ! unpacking; this assigns oAt, oAz through common;
   
  enddo

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  DEALLOCATE(matrix)
  DEALLOCATE(rhs)
  
  DEALLOCATE(rwork)
  DEALLOCATE(residual)

  if( .not. Lposdef ) then
  DEALLOCATE(LU)
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( if04abf(0).ne.0 .or. if04aef(0).ne.0 .or. if04abf(1).ne.0 .or. if04aef(1).ne.0 ) then ! failed to construct Beltrami field and/or derivatives;
   
   ImagneticOK(lvol) = .false. ! set error flag;
   
   if( iflag.eq.1 ) Fmudp(1:Nmudp        ) = zero ! provide dummy intent out;
   if( iflag.eq.2 ) Dmudp(1:Nmudp,1:Nmudp) = zero ! provide dummy intent out;

   iflag = -1 ! this value will be returned by C05PCF to ma02aa;
   
   goto 9999
   
  endif
  
  ImagneticOK(lvol) = .true. ! set error flag; used in fc02aa;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( YESstellsym ) lmns = 1 + (mns-1)           ! number of independent degrees of freedom in angle transformation; 30 Jan 13; 

#ifdef NEWIOTA
! perhaps a rigid shift in the angle does not change the rotational-transform; 02 Sep 14;
  if( NOTstellsym ) lmns = 1 + (mns-1) + (mns  ) ! included non-stellarator symmetric angle transformation; 02 Sep 14;
#else
  if( NOTstellsym ) lmns = 1 + (mns-1) + (mns-1) ! only required for dense, Fourier angle transformation; 21 Apr 13;
#endif
 
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
     
  select case( Lconstraint ) 
   
  case( 0 ) ! Lconstraint=0; 30 Jan 13;
   
   if( Wtr00ab ) then ! construct rotational transform only for diagnostic purposes; 21 Apr 13;
    WCALL(mp00ac,tr00ab,( lvol, mn, lmns, Nt, Nz, Itor, Gpol, iflag, diota(0:1,-1:2,lvol) ) ) ! only required if user requests screen output;
   endif

   Fmudp(1:Nmudp        ) = zero ! provide dummy intent out; Lconstraint.eq.0 indicates no iterations over mu, dpflux are required;
   Dmudp(1:Nmudp,1:Nmudp) = zero ! provide dummy intent out;   
   
  case( 1 ) ! Lconstraint=1; 30 Jan 13;
   
   lcpu = GETTIME
   
   WCALL(mp00ac,tr00ab,( lvol, mn, lmns, Nt, Nz, Itor, Gpol, iflag, diota(0:1,-1:2,lvol) ) ) ! given vector potential, computes transform on interfaces;
   
   cput = GETTIME
   
   if( Lcoordinatesingularity ) then ! Nmudp=1;
    if( iflag.eq.1 ) Fmudp(1    ) = diota(1,  0,lvol) - iota(        lvol )
    if( iflag.eq.2 ) Dmudp(1  ,1) = diota(1,  1,lvol)                        ! derivative of outer rotational-transform wrt helicity multiplier
   endif
   
   if( Nmudp.eq.2 ) then ! this only happens in plasmaregion; 21 Apr 13;
    if( iflag.eq.1 ) Fmudp(1:2  ) = diota(0:1,0,lvol) - (/ oita(lvol-1), iota(lvol) /)
    if( iflag.eq.2 ) Dmudp(1:2,1) = diota(0:1,1,lvol)
    if( iflag.eq.2 ) Dmudp(1:2,2) = diota(0:1,2,lvol)
   endif
   
   if( Lvacuumregion ) then
    if( iflag.eq.1 ) Fmudp(1    ) = diota(0,  0,lvol) -    oita(Nvol  )
    if( iflag.eq.2 ) Dmudp(1  ,1) = diota(0,  1,lvol)                       ! derivative of inner rotational-transform wrt toroidal plasma current;
   endif

! case( 2 ) ! Lconstraint=2; 30 Jan 13;
   
!   if( Wtr00ab ) then
!    lflag = 1
!    WCALL(mp00ac,tr00ab,( lvol, mn, lmns, Nt, Nz, Itor, Gpol, iflag, diota(0:1,-1:2,lvol) ) ) ! only required if user requests screen output;
!   endif
!
!   if( iflag.eq.1 ) then
!    Ldo = 0
!    WCALL(mp00ac,fu00aa,( lvol, Ldo )) ! calculate volume integrals of B.B and A.B;
!    Fmudp(1  ) = lABintegral(lvol) - helicity(lvol)
!   endif
!   
!   if( iflag.eq.2 ) then
!    Ldo = 1
!    WCALL(mp00ac,fu00aa,( lvol, Ldo )) ! calculate volume integrals of B.B and A.B;
!    Dmudp(1,1) = dABintegral
!   endif
   
  case default
   
   FATALMESS(mp00ac, .true., selected Lconstraint not supported)
   
  end select
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Wmp00ac .or. Wma02aa ) then ! the following is screen output; 03 Apr 13; ! 04 Dec 14;
 !if( Wmp00ac              ) then ! the following is screen output; 03 Apr 13; 20 Jun 14;
   
   cput = GETTIME

   if( Lplasmaregion ) then
    select case( iflag )
    case( 0 )    ; write(ounit,3000) cput-cpus, myid, lvol, lmu, dpf, iflag                         ! this is impossible by above logic;
    case( 1 )    ; write(ounit,3000) cput-cpus, myid, lvol, lmu, dpf, iflag, Fmudp(1:Nmudp)
    case( 2 )    ; write(ounit,3010) cput-cpus, myid, lvol, lmu, dpf, iflag, Dmudp(1:Nmudp,1:Nmudp)
    case default ; FATALMESS(mp00ac, .true., illegal iflag on entry)
    end select
   endif
  
3000 format("mp00ac : ",f10.2," : myid=",i3," ; lvol=",i3," ; (mu,dp)=("es23.15" ,"es23.15" ) ; iflag="i2" ;":" F="2es13.05" ;")
3010 format("mp00ac : ",f10.2," : myid=",i3," ; lvol=",i3," ; (mu,dp)=("es23.15" ,"es23.15" ) ; iflag="i2" ;":" D="4es13.05" ;")
  
   if( Lvacuumregion ) then
    select case( iflag )
    case( 0 )    ; write(ounit,3001) cput-cpus, myid, lvol, Itor, Gpol, iflag                         ! this is impossible by above logic;
    case( 1 )    ; write(ounit,3001) cput-cpus, myid, lvol, Itor, Gpol, iflag, Fmudp(1:Nmudp)
    case( 2 )    ; write(ounit,3011) cput-cpus, myid, lvol, Itor, Gpol, iflag, Dmudp(1:Nmudp,1:Nmudp)
    case default ; FATALMESS(mp00ac, .true., illegal iflag on entry)
    end select
   endif
  
3001 format("mp00ac : ",f10.2," : myid=",i3," ; lvol=",i3," ; (I ,G )=("es23.15" ,"es23.15" ) ; iflag="i2" ;":" F="2es13.05" ;")
3011 format("mp00ac : ",f10.2," : myid=",i3," ; lvol=",i3," ; (I ,G )=("es23.15" ,"es23.15" ) ; iflag="i2" ;":" D="4es13.05" ;")

  endif ! end of if( Wmp00ac .or. Wma02aa ) ; 21 Apr 13;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( iflag.eq.1 ) then ! only in this case is Fmudp defined; 17 Apr 13;
   
   if( sum( abs( Fmudp(1:Nmudp) ) ) / Nmudp .lt. mupftol ) then ! satisfactory;  1 Feb 13;
    
    if( Lplasmaregion ) then ; mu(lvol) = lmu  ; dpflux(lvol) = dpf
    else                     ;   curtor = Itor ;!    curpol   = Gpol
    endif
    
    iflag = -2 ! return "acceptance" flag through to ma02aa via ifail;  1 Feb 13; early termination; 
    
   endif ! end of if( sum(Fmudp) ) ; 17 Apr 13;
   
  endif ! end of if( iflag.eq.1 ) ; 17 Apr 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(mp00ac)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine mp00ac

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
