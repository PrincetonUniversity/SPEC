!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (solver) ! Solves Beltrami/vacuum (linear) system, given matrices.

!latex \briefly{Solves for magnetic vector potential given $\boldmu\equiv(\Delta\psi_t,\Delta\psi_p,\mu)^T$.}

!latex \calledby{\link{ma02aa}}
!latex \calls{\link{packab}, \link{curent} and \link{tr00ab}}

!latex \tableofcontents

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{unpacking fluxes, helicity multiplier}

!latex \begin{enumerate}
!latex \item The vector of ``parameters'', $\boldmu$, is unpacked. (Recall that $\boldmu$ was ``packed'' in \link{ma02aa}.)
!latex       In the following, $\boldpsi \equiv (\Delta \psi_t, \Delta \psi_p)^T$.
!latex \end{enumerate}

!latex \subsection{construction of linear system}

!latex \begin{enumerate}
!latex \item The equation $\nabla \times {\bf B} = \mu {\bf B}$ is cast as a matrix equation, \be {\cal M} \cdot {\bf a} = {\cal R},\ee
!latex       where ${\bf a}$ represents the degrees-of-freedom in the magnetic vector potential, ${\bf a} \equiv \{ \Ate{i,l}, \Aze{i,l}, \dots \}$.
!latex \item The matrix ${\cal M}$ is constructed from ${\cal A}\equiv$\internal{dMA} and ${\cal D}\equiv$\internal{dMD},
!latex       which were constructed in \link{matrix}, according to 
!latex       \be {\cal M} \equiv {\cal A} - \mu {\cal D}.
!latex       \ee
!latex       Note that in the vacuum region, $\mu=0$, so ${\cal M} $ reduces to ${\cal M} \equiv {\cal A}$.
!latex \item The construction of the vector ${\cal R}$ is as follows:
!latex       \bi \item [i.] if \internal{Lcoordinatesingularity=T}, then 
!latex           \be {\cal R} \equiv - \left( {\cal B} - \mu {\cal E} \right) \cdot \boldpsi 
!latex           \ee
!latex           \item [ii.] if \internal{Lcoordinatesingularity=F} and \internal{Lplasmaregion=T}, then 
!latex           \be {\cal R} \equiv - {\cal B} \cdot \boldpsi 
!latex           \ee
!latex           \item [iii.] if \internal{Lcoordinatesingularity=F} and \internal{Lvacuumregion=T}, then 
!latex           \be {\cal R} \equiv - {\cal G} - {\cal B} \cdot \boldpsi
!latex           \ee
!latex       \ei 
!latex       The quantities ${\cal B}\equiv$\internal{dMB}, ${\cal E}\equiv$\internal{dME} and ${\cal G}\equiv$\internal{dMG} are constructed in \link{matrix}.
!latex \end{enumerate}

!latex \subsection{solving linear system}

!latex It is {\em not} assumed that the linear system is positive definite.
!latex The \verb+LAPACK+ routine \verb+DSYSVX+ is used to solve the linear system.

!latex \subsection{unpacking, . . .}

!latex \begin{enumerate}
!latex \item The magnetic degrees-of-freedom are unpacked by \link{packab}.
!latex \item The error flag, \internal{ImagneticOK}, is set that indicates if the Beltrami fields were successfully constructed.
!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{construction of ``constraint'' function}

!latex \begin{enumerate}

!latex \item The construction of the function ${\bf f}(\boldmu)$ is required so that iterative methods can be used to construct the Beltrami field
!latex       consistent with the required constraints (e.g. on the enclosed fluxes, helicity, rotational-transform,. . .).
!latex       See \link{ma02aa} for additional details.

!latex \subsubsection{plasma region}

!latex \begin{enumerate}
!latex \item For \type{Lcoordinatesingularity = T}, the returned function is:
!latex \be {\bf f}(\mu,\Delta\psi_p) \equiv 
!latex \left\{ \begin{array}{cccccccr}
!latex (&                                   0&,&                                 0&)^T, & \mbox{\rm if \inputvar{Lconstraint}} &=& -1 \\
!latex (&                                   0&,&                                 0&)^T, & \mbox{\rm if \inputvar{Lconstraint}} &=&  0 \\
!latex (&\iotabar(+1)-\inputvar{iota(lvol  )}&,&                                 0&)^T, & \mbox{\rm if \inputvar{Lconstraint}} &=&  1 \\
!latex (&                                   ?&,&                                 ?&)^T, & \mbox{\rm if \inputvar{Lconstraint}} &=&  2
!latex           \end{array}\right.
!latex \ee
!latex \item For \type{Lcoordinatesingularity = F}, the returned function is:
!latex \be {\bf f}(\mu,\Delta\psi_p) \equiv 
!latex \left\{ \begin{array}{cccccccr}
!latex (&                                   0&,&                                 0&)^T, & \mbox{\rm if \inputvar{Lconstraint}} &=& -1 \\
!latex (&                                   0&,&                                 0&)^T, & \mbox{\rm if \inputvar{Lconstraint}} &=&  0 \\
!latex (&\iotabar(-1)-\inputvar{oita(lvol-1)}&,&\iotabar(+1)-\inputvar{iota(lvol)}&)^T, & \mbox{\rm if \inputvar{Lconstraint}} &=&  1 \\
!latex (&                                   ?&,&                                 ?&)^T, & \mbox{\rm if \inputvar{Lconstraint}} &=&  2
!latex           \end{array}\right.
!latex \ee
!latex \end{enumerate}

!latex \subsubsection{vacuum region}

!latex \begin{enumerate}
!latex \item For the vacuum region, the returned function is:
!latex \be {\bf f}(\Delta\psi_t,\Delta\psi_p) \equiv 
!latex \left\{ \begin{array}{cccccccr}
!latex (&                                   0&,&                                 0&)^T, & \mbox{\rm if \inputvar{Lconstraint}} &=& -1 \\
!latex (&           I-\inputvar{curtor}      &,&           G-\inputvar{curpol}    &)^T, & \mbox{\rm if \inputvar{Lconstraint}} &=&  0 \\
!latex (&\iotabar(-1)-\inputvar{oita(lvol-1)}&,&           G-\inputvar{curpol}    &)^T, & \mbox{\rm if \inputvar{Lconstraint}} &=&  1 \\
!latex (&                                   ?&,&                                 ?&)^T, & \mbox{\rm if \inputvar{Lconstraint}} &=&  2
!latex           \end{array}\right.
!latex \ee
!latex \end{enumerate}

!latex \item The rotational-transform, $\iotabar$, is computed by \link{tr00ab}; and the enclosed currents, $I$ and $G$, are computed by \link{curent}.

!latex \end{enumerate}

!latex \subsection{early termination}

!latex \begin{enumerate}
!latex \item If $|{\bf f}| < $ \inputvar{mupftol}, then early termination is enforced (i.e., \internal{iflag} is set to negative integer).
!latex       (See \link{ma02aa} for details of how \link{mp00ac} is called iteratively.)
!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine mp00ac( Ndof, Xdof, Fdof, Ddof, Ldfjac, iflag ) ! argument list is fixed by NAG; ma02aa calls mp00ac through C05PCF;

! if iflag.eq.0 : Xdof and Fdof are available for PRINTING ; Fdof MUST NOT BE CHANGED; Ddof MUST NOT BE CHANGED;
! if iflag.eq.1 :          Fdof is to be          UPDATED  ;                         ; Ddof MUST NOT BE CHANGED;
! if iflag.eq.2 :          Ddof is to be          UPDATED  ; Fdof MUST NOT BE CHANGED;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half, one

  use numerical, only : small, machprec
  
  use fileunits, only : ounit
  
  use inputlist, only : Wmacros, Wmp00ac, Wtr00ab, Wcurent, Wma02aa, &
                        mu, helicity, iota, oita, curtor, curpol, Lrad, Ntor,&
                       !Lposdef, &
                        Lconstraint, mupftol, &
                        Lmatsolver, NiterGMRES, epsGMRES, LGMRESprec, epsILU
  
  use cputiming, only : Tmp00ac
  
  use allglobal, only : myid, ncpu, cpus, ivol, &
                        YESstellsym, NOTstellsym, &
                        Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, &
                        mn, im, in, mns, &
                        Nt, Nz, & ! only required to pass through as arguments to tr00ab;
                        NAdof, &
!                       dMA, dMB, dMC, dMD, dME, dMF, dMG, &
                        dMA, dMB,      dMD,           dMG, &
                        NdMASmax, NdMAS, dMAS, dMDS, idMAS, jdMAS, & ! preconditioning matrix
                        solution, GMRESlastsolution, &
                        dtflux, dpflux, &
                        diotadxup, dItGpdxtp, &
                        lBBintegral, lABintegral, &
                        xoffset, &
                        ImagneticOK, &
                        Ate, Aze, Ato, Azo, Mvol, Iquad, &
                        LILUprecond, GMRESlastsolution, NOTMatrixFree
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS

  INTEGER, intent(in)  :: Ndof, Ldfjac
  REAL   , intent(in)  :: Xdof(1:Ndof)
  REAL                 :: Fdof(1:Ndof), Ddof(1:Ldfjac,1:Ndof)
  INTEGER              :: iflag ! indicates whether (i) iflag=1: ``function'' values are required; or (ii) iflag=2: ``derivative'' values are required;
  
  
  INTEGER, parameter   :: NB = 4 ! optimal workspace block size for LAPACK:DGECON;
  
  INTEGER              :: lvol, NN, MM, ideriv, lmns, ii, jj, nnz, Lwork

  INTEGER              :: idgetrf(0:1), idgetrs(0:1), idgerfs(0:1), idgecon(0:1)
  
  REAL                 :: lmu, dpf, dtf, dpsi(1:2), tpsi(1:2), ppsi(1:2), lcpu, test(2,2)
  
  REAL                 :: anorm, rcond, ferr(2), berr(2), signfactor

  CHARACTER            :: packorunpack
  
  ! For direct LU decompose
  INTEGER, allocatable :: ipiv(:), Iwork(:)

  REAL   , allocatable :: matrix(:,:), rhs(:,:), LU(:,:)

  REAL   , allocatable :: RW(:), RD(:,:)

  REAL   , allocatable :: matrixC(:,:)

! For GMRES + ILU
  INTEGER, parameter   :: nrestart = 5 ! do GMRES restart after nrestart iterations 
  INTEGER              :: maxfil   ! bandwidth for ILU subroutines, will be estimated

  INTEGER              :: NS, itercount, Nbilut

  REAL   , allocatable :: matrixS(:), bilut(:)
  INTEGER, allocatable :: ibilut(:),  jbilut(:)
  
  INTEGER, parameter   :: ipar_SIZE = 128
  INTEGER              :: ipar(ipar_SIZE), iluierr, RCI_REQUEST, nw, t1, t2, t3
  REAL                 :: fpar(ipar_SIZE), v1
  REAL, allocatable    :: wk(:)
  INTEGER,allocatable  :: jw(:), iperm(:)

  BEGIN(mp00ac)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  lvol = ivol ! recall that ivol is global;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  FATAL( mp00ac, iflag.ne.1 .and. iflag.ne.2, invalid iflag ) ! see nprint=0 in ma02aa and C05PCF; perhaps NAG:C05PCF is no longer used;
#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Lplasmaregion ) then
   
   ;                    ; lmu  = Xdof(1) - xoffset
   ;                    ; dtf  = dtflux(lvol)
   if( Ndof.eq.2 ) then ; dpf  = Xdof(2) - xoffset
   else                 ; dpf  = dpflux(lvol)
   endif
   
  else ! Lvacuumregion;
   
#ifdef FORCEFREEVACUUM
   ;                    ; lmu  = mu(lvol)           ! generalize for arbitrary force-free field;
#else
   ;                    ; lmu  = zero               ! restrict attention to strict vacuum field;
#endif
   
   ;                    ; dtf  = Xdof(1) - xoffset
   if( Ndof.eq.2 ) then ; dpf  = Xdof(2) - xoffset
   else                 ; dpf  = dpflux(lvol)
   endif
   
  endif ! end of if( Lplasmaregion ) ;
  
  dpsi(1:2) = (/  dtf,  dpf /) ! enclosed poloidal fluxes and their derivatives;
  tpsi(1:2) = (/  one, zero /) ! enclosed toroidal fluxes and their derivatives;
  ppsi(1:2) = (/ zero,  one /) ! enclosed toroidal fluxes and their derivatives;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  diotadxup(0:1,-1:2,lvol) = zero ! rotational-transform, and its derivatives with respect to lmu and dpf, or toroidal current, on the inner/outer interface;
  dItGpdxtp(0:1,-1:2,lvol) = zero ! plasma and linking currents;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  NN = NAdof(lvol) ! shorthand;
  
  SALLOCATE( rhs   , (1:NN,0:2 ), zero )
  if (NOTMatrixFree) then ! create the full size matrix
    SALLOCATE( matrix, (1:NN,1:NN), zero )
  else                    ! create a dummy variable
    SALLOCATE( matrix, (1:1,1:1), zero )
  endif

  solution(1:NN,-1:2) = zero ! this is a global array allocated in dforce;

  ! allocate work space
  select case (Lmatsolver)
  case (1) ! direct matrix solver
    Lwork = NB*NN

    SALLOCATE( RW,    (1:Lwork ),  zero )
    SALLOCATE( RD,    (1:NN,0:2),  zero )
    SALLOCATE( LU,    (1:NN,1:NN), zero )
    SALLOCATE( ipiv,  (1:NN),         0 )
    SALLOCATE( Iwork, (1:NN),         0 )
  case (2:3) ! GMRES
    if (LILUprecond) then
      NS = NdMAS(lvol) ! shorthand
      SALLOCATE( matrixS, (1:NS), zero )
      
      ! estimate bandwidth
      if (Lcoordinatesingularity) then
        maxfil = Lrad(lvol) + 10
        if (NOTstellsym) maxfil = maxfil + Lrad(lvol) + 10 
      else
        maxfil = 2 * Lrad(lvol) + 10
        if (NOTstellsym) maxfil = maxfil + 2 * Lrad(lvol) + 10
      end if

      Nbilut = (2*maxfil+2)*NN
      SALLOCATE( bilut, (1:Nbilut), zero)
      SALLOCATE( jbilut, (1:Nbilut), 0)
      SALLOCATE( ibilut, (1:NN+1), 0)

    endif
    nw = (NN+3)*(nrestart+2) + (nrestart+1)*nrestart
    SALLOCATE( wk, (1:nw), zero)
    SALLOCATE( jw, (1:2*NN), 0)
    SALLOCATE( iperm, (1:2*NN), 0)
  end select
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  idgetrf(0:1) = 0 ! error flags;
  idgetrs(0:1) = 0 ! error flags;
  idgerfs(0:1) = 0 ! error flags;
  idgecon(0:1) = 0 ! error flags;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! ideriv labels derivative as follows:
! ideriv = 0 : compute Beltrami field; ideriv = 1 : compute d Beltrami field / d \mu           ; ideriv = 2 : compute d Beltrami field / d \Delta \psi_p ;
! ideriv = 0 : compute Vacuum   field; ideriv = 1 : compute d Vacuum   field / d \Delta \psi_t ; ideriv = 2 : compute d Vacuum   field / d \Delta \psi_p ;
    
  do ideriv = 0, 1 ! loop over derivatives;
   
   if( iflag.eq.1 .and. ideriv.eq.1 ) cycle ! only need to return function; recall the derivative estimate requires function evaluation;
   
   if( Lcoordinatesingularity ) then
    
    if (NOTMatrixFree) then
      ;matrix(1:NN,1:NN) = dMA(1:NN,1:NN) - lmu * dMD(1:NN,1:NN)
    endif
    
!  !;select case( ideriv )
!  !;case( 0 )    ; rhs(1:NN,0) = - matmul(  dMB(1:NN,1:2) - lmu  * dME(1:NN,1:2), dpsi(1:2) )
!  !;case( 1 )    ; rhs(1:NN,1) = - matmul(                - one  * dME(1:NN,1:2), dpsi(1:2) ) - matmul( - one  * dMD(1:NN,1:NN), solution(1:NN,0) )
!  !; ;           ; rhs(1:NN,2) = - matmul(  dMB(1:NN,1:2) - lmu  * dME(1:NN,1:2), ppsi(1:2) )
!  !;end select
    
    ;select case( ideriv )
    ;case( 0 )    ; rhs(1:NN,0) = - matmul(  dMB(1:NN,1:2)                       , dpsi(1:2) )
    ;case( 1 )    ! construct dMD*solution
    if (NOTMatrixFree) then
    ; ;           ; rhs(1:NN,1) =                                                              - matmul( - one  * dMD(1:NN,1:NN), solution(1:NN,0) )
    else ! Matrix free version
    ; ;call intghs(Iquad(lvol), mn, lvol, Lrad(lvol), 0) ! compute the integrals of B_lower
    ; ;call mtrxhs(lvol, mn, Lrad(lvol), wk(1:NN+1), wk(NN+2:2*NN+2), 0) ! construct a.x from the integral
    ; ;           ; rhs(1:NN,1) = wk(NN+3:2*NN+2)
    endif ! NOTMatrixFree
    ; ;           ; rhs(1:NN,2) = - matmul(  dMB(1:NN,1:2)                       , ppsi(1:2) )
    ;end select
    
   else ! .not.Lcoordinatesingularity; 
    
    if( Lplasmaregion ) then

     if (NOTMatrixFree) then
       matrix(1:NN,1:NN) = dMA(1:NN,1:NN) - lmu * dMD(1:NN,1:NN)
     endif

    ;select case( ideriv )
    ;case( 0 )    ; rhs(1:NN,0) = - matmul(  dMB(1:NN,1:2)                       , dpsi(1:2) )
    ;case( 1 )    ! construct dMD*solution
    if (NOTMatrixFree) then
    ; ;           ; rhs(1:NN,1) =                                                              - matmul( - one  * dMD(1:NN,1:NN), solution(1:NN,0) )
    else ! Matrix free version
    ; ;call intghs(Iquad(lvol), mn, lvol, Lrad(lvol), 0) ! compute the integrals of B_lower
    ; ;call mtrxhs(lvol, mn, Lrad(lvol), wk(1:NN+1), wk(NN+2:2*NN+2), 0) ! construct a.x from the integral
    ; ;           ; rhs(1:NN,1) = wk(NN+3:2*NN+2)
    endif ! NOTMatrixFree
    ; ;           ; rhs(1:NN,2) = - matmul(  dMB(1:NN,1:2)                       , ppsi(1:2) )
    ;end select
     
    else ! Lvacuumregion ;
     
#ifdef FORCEFREEVACUUM
     FATAL( mp00ac, .true., need to revise Beltrami matrices in vacuum region for arbitrary force-free field )
#else
     
     if (NOTMatrixFree) then
       matrix(1:NN,1:NN) = dMA(1:NN,1:NN) ! - lmu * dMD(1:NN,1:NN) ;
     endif

     select case( ideriv )
     case( 0 )    ; rhs(1:NN,0) = - dMG(1:NN) - matmul( dMB(1:NN,1:2), dpsi(1:2) ) ! perhaps there is an lmu term missing here;
     case( 1 )    ; rhs(1:NN,1) =             - matmul( dMB(1:NN,1:2), tpsi(1:2) ) ! perhaps there is an lmu term missing here;
      ;           ; rhs(1:NN,2) =             - matmul( dMB(1:NN,1:2), ppsi(1:2) ) ! perhaps there is an lmu term missing here;
     end select
#endif

    endif ! end of if( Lplasmaregion ) ;

   endif ! end of if( Lcoordinatesingularity ) ;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   select case( Lmatsolver )

   case(1) ! Using direct matrix solver (LU factorization), must not be matrix free
    
    select case( ideriv )
      
    case( 0 ) ! ideriv=0;
      
      MM = 1
      !LU = matrix
      call DCOPY(NN*NN, matrix, 1, LU, 1) ! BLAS version
      solution(1:NN,0   ) = rhs(:,0   )
      call DGETRF(NN, NN, LU, NN, ipiv, idgetrf(ideriv) ) ! LU factorization

      anorm=maxval(sum(abs(matrix),1))
      call DGECON('I', NN, LU, NN, anorm, rcond, RW, Iwork, idgecon(ideriv)) ! estimate the condition number

      call DGETRS('N', NN, MM, LU, NN, ipiv, solution(1:NN,0   ), NN, idgetrs(ideriv) ) ! sovle linear equation
      call DGERFS('N', NN, MM, matrix, NN, LU, NN, ipiv, rhs(1,0), NN, solution(1,0), NN, FERR, BERR, RW, Iwork, idgerfs(ideriv)) ! refine the solution
    case( 1 ) ! ideriv=1;
      
      MM = 2
      solution(1:NN,1:MM) = rhs(:,1:MM)
      call DGETRS( 'N', NN, MM, LU, NN, ipiv, solution(1:NN,1:MM), NN, idgetrs(ideriv) )
      call DGERFS('N', NN, MM, matrix, NN, LU, NN, ipiv, rhs(1,1), NN, solution(1,1), NN, FERR, BERR, RW, Iwork, idgerfs(ideriv))

    end select ! ideriv;
    
    cput = GETTIME

    if(     idgetrf(ideriv) .eq. 0 .and. idgetrs(ideriv) .eq. 0 .and. idgerfs(ideriv) .eq. 0 .and. rcond .ge. machprec) then
      if( Wmp00ac ) write(ounit,1010) cput-cpus, myid, lvol, ideriv, "idgetrf idgetrs idgerfs", idgetrf(ideriv), idgetrs(ideriv), idgetrf(ideriv), "success ;         ", cput-lcpu	   
    elseif( idgetrf(ideriv) .lt. 0 .or. idgetrs(ideriv) .lt. 0 .or. idgerfs(ideriv) .lt. 0   ) then
      ;             write(ounit,1010) cput-cpus, myid, lvol, ideriv, "idgetrf idgetrs idgerfs", idgetrf(ideriv), idgetrs(ideriv), idgetrf(ideriv), "input error ;     "
    elseif( idgetrf(ideriv) .gt. 0 ) then
      ;             write(ounit,1010) cput-cpus, myid, lvol, ideriv, "idgetrf idgetrs idgerfs", idgetrf(ideriv), idgetrs(ideriv), idgetrf(ideriv), "singular ;        "
    elseif( rcond .le. machprec) then
      ;             write(ounit,1010) cput-cpus, myid, lvol, ideriv, "idgetrf idgetrs idgerfs", idgetrf(ideriv), idgetrs(ideriv), idgetrf(ideriv), "ill conditioned ; "
    else
      ;             write(ounit,1010) cput-cpus, myid, lvol, ideriv, "idgetrf idgetrs idgerfs", idgetrf(ideriv), idgetrs(ideriv), idgetrf(ideriv), "invalid error ; "
    endif
    
  1010 format("mp00ac : ",f10.2," : myid=",i3," ; lvol=",i3," ; ideriv="i2" ; "a23"=",i3,' ',i3,' ',i3," ; "a34,:" time=",f10.2," ;")
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   case(2:3) ! Using GMRES, can be matrix free

    select case (ideriv)

    case (0) ! ideriv = 0

      if (LILUprecond) then ! Using ILU preconditioner
        matrixS(1:NS) = dMAS(1:NS) - lmu * dMDS(1:NS)  ! construct the sparse precondtioner matrix
      end if ! if (LILUprecond)
      
      if (MAXVAL(abs(GMRESlastsolution(1:NN,0,lvol))) .le. small) then
        GMRESlastsolution(1:NN,0,lvol) = zero
      endif

      if (LILUprecond) then
      ! ILU factorization
        call ilutp(NN,matrixS,jdMAS,idMAS,maxfil,epsILU,0.1,NN,bilut,jbilut,ibilut,Nbilut,wk,jw,iperm,iluierr)
        FATAL(mp00ac, iluierr.ne.0, construction of preconditioner failed)
      endif

      call rungmres(NN,nrestart,lmu,lvol,rhs(1:NN,0),solution(1:NN,0),ipar,fpar,wk,nw,GMRESlastsolution(1:NN,0,lvol),matrix,bilut,jbilut,ibilut,iperm,ierr)
      
      if (ierr .gt. 0) then
        ImagneticOK(lvol) = .true.
        GMRESlastsolution(1:NN,0,lvol) = solution(1:NN,0)
      elseif (ierr .eq. 0) then
        solution(1:NN,0) = GMRESlastsolution(1:NN,0,lvol)
      else
        ImagneticOK(lvol) = .false.
      endif 

    case (1) ! ideriv = 1

      do ii = 1, 2
        if (MAXVAL(abs(GMRESlastsolution(1:NN,ii,lvol))) .le. small) then
          GMRESlastsolution(1:NN,ii,lvol) = zero
        endif

        call rungmres(NN,nrestart,lmu,lvol,rhs(1:NN,ii),solution(1:NN,ii),ipar,fpar,wk,nw,GMRESlastsolution(1:NN,ii,lvol),matrix,bilut,jbilut,ibilut,iperm,ierr)

        if (ierr .gt. 0) then
          GMRESlastsolution(1:NN,ii,lvol) = solution(1:NN,ii)
        elseif (ierr .eq. 0) then
          solution(1:NN,0) = GMRESlastsolution(1:NN,ii,lvol)
        else
          exit
        endif
        
      end do ! ii
    end select ! ideriv

    cput = GETTIME

    if (ierr.ge.0) then
      if( Wmp00ac ) write(ounit,1011) cput-cpus, myid, lvol, ideriv, ierr, " successful ; "
    elseif (ierr.eq.-1) then
        ;           write(ounit,1011) cput-cpus, myid, lvol, ideriv, ierr, " max niter, epsGMRES not reached ; "
    elseif (ierr.eq.-2) then
        ;           write(ounit,1011) cput-cpus, myid, lvol, ideriv, ierr, " more workspace needed ; "
    elseif (ierr.eq.-3) then
        ;           write(ounit,1011) cput-cpus, myid, lvol, ideriv, ierr, " solver internal break down ; "  
    else
        ;           write(ounit,1011) cput-cpus, myid, lvol, ideriv, ierr, " check iters.f for error code ; "
    end if

1011 format("mp00ac : ",f10.2," : myid=",i3," ; lvol=",i3," ; ideriv=",i2," ; GMRES ierr=",i4, " ; "a34" ")

   end select ! Lmatsolver 

  enddo ! end of do ideriv;
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  do ideriv = 0, 2
   
   if( iflag.eq.1 .and. ideriv.gt.0 ) cycle
   
   packorunpack = 'U'
   WCALL( mp00ac, packab, ( packorunpack, lvol, NN, solution(1:NN,ideriv), ideriv ) ) ! unpacking; this assigns oAt, oAz through common;
   
  enddo ! do ideriv = 0, 2;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
! can compute the energy and helicity integrals; easiest to do this with solution in packed format;
  
  if (NOTMatrixFree) then

   lBBintegral(lvol) = half * sum( solution(1:NN,0) * matmul( dMA(1:NN,1:NN), solution(1:NN,0) ) ) & 
                     +        sum( solution(1:NN,0) * matmul( dMB(1:NN,1: 2),     dpsi(1: 2  ) ) ) !
!                    + half * sum(     dpsi(1: 2  ) * matmul( dMC(1: 2,1: 2),     dpsi(1: 2  ) ) )
  
   lABintegral(lvol) = half * sum( solution(1:NN,0) * matmul( dMD(1:NN,1:NN), solution(1:NN,0) ) ) ! 
!                    +        sum( solution(1:NN,0) * matmul( dME(1:NN,1: 2),     dpsi(1: 2  ) ) ) !
!                    + half * sum(     dpsi(1: 2  ) * matmul( dMF(1: 2,1: 2),     dpsi(1: 2  ) ) )
  else
    call intghs(Iquad(lvol), mn, lvol, Lrad(lvol), 0) ! compute the integrals of B_lower
    call mtrxhs(lvol, mn, Lrad(lvol), wk(1:NN+1), wk(NN+2:2*NN+2), 0) ! construct a.x from the integral

    lBBintegral(lvol) = half * sum( solution(1:NN,0) * wk(2:NN+1) ) & 
                      +        sum( solution(1:NN,0) * matmul( dMB(1:NN,1: 2),     dpsi(1: 2  ) ) ) !

    lABintegral(lvol) = half * sum( solution(1:NN,0) * wk(NN+3:2*NN+2) )
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!  
  DALLOCATE( matrix )
  DALLOCATE( rhs    )

  select case (Lmatsolver)
  case (1) ! LU
      DALLOCATE( RW )
      DALLOCATE( RD )
      DALLOCATE( LU )
      DALLOCATE( ipiv )
      DALLOCATE( Iwork )
  case (2:3) ! GMRES
    if (LILUprecond) then
      DALLOCATE( matrixS )
      DALLOCATE( bilut )
      DALLOCATE( jbilut )
      DALLOCATE( ibilut )
    endif
    DALLOCATE( wk )
    DALLOCATE( jw )
    DALLOCATE( iperm )
  end select    

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  idgetrf(0:1) = abs(idgetrf(0:1)) + abs(idgetrs(0:1)) + abs(idgerfs(0:1)) + abs(idgecon(0:1))
  if( idgetrf(0).ne.0 .or. idgetrf(1).ne.0 ) then ! failed to construct Beltrami/vacuum field and/or derivatives;
   
   ImagneticOK(lvol) = .false. ! set error flag;
   
   if( iflag.eq.1 ) Fdof(1:Ndof       ) = zero ! provide dummy intent out;
   if( iflag.eq.2 ) Ddof(1:Ndof,1:Ndof) = zero ! provide dummy intent out;
   
   iflag = -1 ! this value will be returned by C05PCF to ma02aa;
   
   goto 9999
   
  else
   
   ImagneticOK(lvol) = .true. ! set error flag; used in dforce;
  
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( YESstellsym ) then ; lmns = 1 + (mns-1)           ! number of independent degrees of freedom in angle transformation;
  else                   ; lmns = 1 + (mns-1) + (mns-1) ! only required for dense, Fourier angle transformation;
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!#ifdef NEWIOTA
!! perhaps a rigid shift in the angle does not change the rotational-transform;
!  if( NOTstellsym ) lmns = 1 + (mns-1) + (mns  ) ! included non-stellarator symmetric angle transformation;
!#else
!  if( NOTstellsym ) lmns = 1 + (mns-1) + (mns-1) ! only required for dense, Fourier angle transformation;
!#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  select case( Lconstraint ) 
   
  case( -1 ) ! Lconstraint=-1;
   
   if( Lplasmaregion ) then
    
    if( Wtr00ab ) then ! compute rotational transform only for diagnostic purposes;
     WCALL( mp00ac, tr00ab, ( lvol, mn, lmns, Nt, Nz, iflag, diotadxup(0:1,-1:2,lvol) ) )
    endif
    
    Fdof(1:Ndof       ) = zero ! provide dummy intent out; Lconstraint=-1 indicates no iterations over mu   , dpflux are required;
    Ddof(1:Ndof,1:Ndof) = zero ! provide dummy intent out;   
    
   else ! Lvacuumregion
    
    if( Wtr00ab ) then ! compute rotational transform only for diagnostic purposes;
     WCALL( mp00ac, tr00ab, ( lvol, mn, lmns, Nt, Nz, iflag, diotadxup(0:1,-1:2,lvol) ) )
    endif
    
    if( Wcurent ) then ! compute enclosed currents    only for diagnostic purposes;
     WCALL( mp00ac, curent,( lvol, mn, Nt, Nz, iflag, dItGpdxtp(0:1,-1:2,lvol) ) )
     curtor = dItGpdxtp(0,0,lvol) ! icurrent(0) ! update input variables;
     curpol = dItGpdxtp(1,0,lvol) ! gcurrent(0)
    endif
    
    Fdof(1:Ndof       ) = zero ! provide dummy intent out;Lconstraint=-1 indicates no iterations over dtflux, dpflux are required;
    Ddof(1:Ndof,1:Ndof) = zero ! provide dummy intent out;   
    
   endif ! end of if( Lplasmaregion) ;

  case(  0 ) ! Lconstraint= 0;
   
   if( Lplasmaregion ) then
    
    if( Wtr00ab ) then ! compute rotational transform only for diagnostic purposes;
     WCALL( mp00ac, tr00ab, ( lvol, mn, lmns, Nt, Nz, iflag, diotadxup(0:1,-1:2,lvol) ) )
    endif
    
    Fdof(1:Ndof       ) = zero ! provide dummy intent out; Lconstraint= 0 indicates no iterations over mu, dpflux are required;
    Ddof(1:Ndof,1:Ndof) = zero ! provide dummy intent out;   
    
   else ! Lvacuumregion
    
    WCALL( mp00ac, curent,( lvol, mn, Nt, Nz, iflag, dItGpdxtp(0:1,-1:2,lvol) ) )
    
    if( iflag.eq.1 ) Fdof(1:2  ) = (/ dItGpdxtp(0,0,lvol) - curtor, dItGpdxtp(1,0,lvol) - curpol /)
    if( iflag.eq.2 ) Ddof(1:2,1) = (/ dItGpdxtp(0,1,lvol)         , dItGpdxtp(1,1,lvol)          /)
   !if( iflag.eq.2 ) Ddof(1:2,2) = (/ dItGpdxtp(0,0,lvol)         , dItGpdxtp(1,2,lvol)          /)
    if( iflag.eq.2 ) Ddof(1:2,2) = (/ dItGpdxtp(0,2,lvol)         , dItGpdxtp(1,2,lvol)          /)
    
   endif ! end of if( Lplasmaregion) ;
   
  case(  1 ) ! Lconstraint= 1;

   WCALL( mp00ac, tr00ab,( lvol, mn, lmns, Nt, Nz, iflag, diotadxup(0:1,-1:2,lvol) ) ) ! required for both plasma and vacuum region;
   
   if( Lplasmaregion ) then
    
    if( Lcoordinatesingularity ) then ! Ndof = 1;
     if( iflag.eq.1 ) Fdof(1    ) = diotadxup(1,  0,lvol) - iota(        lvol )
     if( iflag.eq.2 ) Ddof(1  ,1) = diotadxup(1,  1,lvol)                        ! derivative of outer rotational-transform wrt helicity multiplier
    endif
    
    if( Ndof.eq.2 ) then
     if( iflag.eq.1 ) Fdof(1:2  ) = diotadxup(0:1,0,lvol) - (/ oita(lvol-1), iota(lvol) /)
     if( iflag.eq.2 ) Ddof(1:2,1) = diotadxup(0:1,1,lvol)
     if( iflag.eq.2 ) Ddof(1:2,2) = diotadxup(0:1,2,lvol)
    endif
    
   else ! Lvacuumregion

    WCALL( mp00ac, curent, ( lvol, mn,     Nt, Nz, iflag, dItGpdxtp(0:1,-1:2,lvol) ) )

    curtor = dItGpdxtp(0,0,lvol) ! update input variables; 08 Jun 16; 
   !curpol = dItGpdxtp(1,0,lvol)

    if( iflag.eq.1 ) Fdof(1:2  ) = (/ diotadxup(0,0,lvol) - oita(lvol-1), dItGpdxtp(1,0,lvol) - curpol /)
    if( iflag.eq.2 ) Ddof(1:2,1) = (/ diotadxup(0,1,lvol)               , dItGpdxtp(1,1,lvol)          /)
    if( iflag.eq.2 ) Ddof(1:2,2) = (/ diotadxup(0,2,lvol)               , dItGpdxtp(1,2,lvol)          /)

   endif ! end of if( Lplasmaregion) ;

  case(  2 )

   if ( iflag.eq.1 ) Fdof(1     ) = lABintegral(lvol) - helicity(lvol)
   if ( iflag.eq.2 ) Ddof(1   ,1) = half * sum( solution(1:NN,1) * matmul( dMD(1:NN,1:NN), solution(1:NN,0) ) ) &
                                  + half * sum( solution(1:NN,0) * matmul( dMD(1:NN,1:NN), solution(1:NN,1) ) )

  case(  3 )

   if( Lplasmaregion ) then
        
    if( Wtr00ab ) then ! compute rotational transform only for diagnostic purposes;
     WCALL( mp00ac, tr00ab, ( lvol, mn, lmns, Nt, Nz, iflag, diotadxup(0:1,-1:2,lvol) ) )
    endif
    
    Fdof(1:Ndof       ) = zero ! provide dummy intent out; no iteration other mu and psip locally
    Ddof(1:Ndof,1:Ndof) = zero ! provide dummy intent out;   
    
   else ! Lvacuumregion
    
    !WCALL( mp00ac, curent,( lvol, mn, Nt, Nz, iflag, dItGpdxtp(0:1,-1:2,lvol) ) )
    
    ! Iteration only on toroidal flux.
    ! if( iflag.eq.1 ) Fdof(1:Ndof  ) = dItGpdxtp(1,0,lvol) - curpol
    ! if( iflag.eq.2 ) Ddof(1:Ndof,1:Ndof) = dItGpdxtp(1,1,lvol) 
    
    Fdof(1:Ndof       ) = zero ! provide dummy intent out; no iteration other mu and psip locally
    Ddof(1:Ndof,1:Ndof) = zero ! provide dummy intent out;
    
   endif ! end of if( Lplasmaregion) ;

  end select ! end of select case( Lconstraint ) ;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Wmp00ac .or. Wma02aa ) then ! the following is screen output;
   
   cput = GETTIME
   
   if( Lplasmaregion ) then
    select case( iflag )
    case( 0 )    ; write(ounit,3000) cput-cpus, myid, lvol, lmu, dpf, iflag                         ! this is impossible by above logic;
    case( 1 )    ; write(ounit,3000) cput-cpus, myid, lvol, lmu, dpf, iflag, Fdof(1:Ndof)
    case( 2 )    ; write(ounit,3010) cput-cpus, myid, lvol, lmu, dpf, iflag, Ddof(1:Ndof,1:Ndof)
    case default ; FATAL( mp00ac, .true., illegal iflag on entry )
    end select
   else ! Lvacuumregion
    select case( iflag )
    case( 0 )    ; write(ounit,3001) cput-cpus, myid, lvol, dtf, dpf, iflag                         ! this is impossible by above logic;
    case( 1 )    ; write(ounit,3001) cput-cpus, myid, lvol, dtf, dpf, iflag, Fdof(1:Ndof)
    case( 2 )    ; write(ounit,3011) cput-cpus, myid, lvol, dtf, dpf, iflag, Ddof(1:Ndof,1:Ndof)
    case default ; FATAL( mp00ac, .true., illegal iflag on entry )
    end select
   endif
  
3000 format("mp00ac : ",f10.2," : myid=",i3," ; lvol=",i3," ; (mu,dp)=("es23.15" ,"es23.15" ) ; iflag="i2" ;":" F="2es13.05" ;")
3010 format("mp00ac : ",f10.2," : myid=",i3," ; lvol=",i3," ; (mu,dp)=("es23.15" ,"es23.15" ) ; iflag="i2" ;":" D="4es13.05" ;")
3001 format("mp00ac : ",f10.2," : myid=",i3," ; lvol=",i3," ; (dt,dp)=("es23.15" ,"es23.15" ) ; iflag="i2" ;":" F="2es13.05" ;")
3011 format("mp00ac : ",f10.2," : myid=",i3," ; lvol=",i3," ; (dt,dp)=("es23.15" ,"es23.15" ) ; iflag="i2" ;":" D="4es13.05" ;")
  
  endif ! end of if( Wmp00ac .or. Wma02aa ) ;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( iflag.eq.1 ) then ! only in this case is Fdof defined;
   
   if( sum( abs( Fdof(1:Ndof) ) ) / Ndof .lt. mupftol ) then ! satisfactory;
    
    if ( Lplasmaregion ) then ; mu(lvol) = lmu  ;                    ; dpflux(lvol) = dpf
#ifdef FORCEFREEVACUUM
    else                      ; mu(lvol) = lmu  ; dtflux(lvol) = dtf ; dpflux(lvol) = dpf
#else
    else                      ; mu(lvol) = zero ; dtflux(lvol) = dtf ; dpflux(lvol) = dpf
#endif
    endif
    
    iflag = -2 ! return "acceptance" flag through to ma02aa via ifail; early termination;
    
   endif ! end of if( sum(Fdof) ) ;
   
  endif ! end of if( iflag.eq.1 ) ;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(mp00ac)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine mp00ac

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine rungmres(n,nrestart,mu,vvol,rhs,sol,ipar,fpar,wk,nw,guess,a,au,jau,ju,iperm,ierr)
  ! Driver subroutine for GMRES
  ! modified from riters.f from SPARSKIT v2.0 
  ! by ZSQ 02 Feb 2020
  use constants, only : zero, one
  use inputlist, only : NiterGMRES, epsGMRES
  use allglobal, only : LILUprecond
  use fileunits
  implicit none
  INTEGER  :: n, nrestart, nw, vvol, ju(*), jau(*), iperm(*)
  INTEGER  :: ipar(16)
  INTEGER  :: ierr
  REAL     :: guess(n), au(*), mu
  REAL     :: fpar(16), rhs(1:n), sol(1:n), wk(1:nw), a(*)

  INTEGER :: i, its
  REAL :: res, tmprhs(1:n)

  its = 0
  res = zero
  wk = zero

  ! setup solver parameters
  ipar = 0
  fpar = zero
  wk = zero
  if (LILUprecond) then
    ipar(2) = 1 ! tell GMRES to use preconditioner
  else
    ipar(2) = 0 ! do not use preconditioner, not recommended
  end if
  ipar(3) = 2           ! stopping test type, see iters.f
  ipar(4) = nw          ! length of work array
  ipar(5) = nrestart    ! restart parameter, size of Krylov subspace
  ipar(6) = NiterGMRES  ! maximum number of iteration
  fpar(1) = epsGMRES    ! stop criterion, relative error
  fpar(2) = epsGMRES    ! stop criterion, absolute error

  sol(1:n) = guess(1:n)
  tmprhs(1:n) = rhs(1:n) ! copy to a temp vector because it will be destroyed by gmres

  ipar(1) = 999
  do while (ipar(1) .gt. 0) ! main reversed communication loop

    call gmres(n,tmprhs,sol,ipar,fpar,wk)
    res = fpar(5) ! shorthand for resolution
    its = ipar(7) ! shorthand for iteration number

    if (ipar(1).eq.1) then ! compute A.x
      ! we should compute wk(ipar(9) = matmul(matrix, wk(ipar(8)))
      call matvec(n, wk(ipar(8)), wk(ipar(9)), a, mu, vvol)

    else if (ipar(1).eq.3) then
      ! apply the preconditioner
      call prec_solve(n,wk(ipar(8)),wk(ipar(9)),au,jau,ju,iperm)

    else if (ipar(1).lt.0) then 
      ! error or max iter reached, exit
      exit
    endif
  
  end do ! end main loop

  ierr = ipar(1)
  if (ierr.eq.0) ierr = its

  return
end subroutine rungmres

subroutine matvec(n, x, ax, a, mu, vvol)
  ! compute a.x by either by coumputing it directly, 
  ! or using a matrix free method
  use constants, only : zero, one
  use fileunits, only : ounit
  use inputlist, only : Lrad, Wmp00ac
  use cputiming, only : Tmp00ac
  use allglobal, only : NOTMatrixFree, Iquad, mn, dmd, cpus, myid
! implicit none
  LOCALS
  INTEGER, intent(in) :: n, vvol
  REAL                :: ax(1:n), x(1:n), a(*), mu
  INTEGER             :: ideriv
  REAL                :: dax(0:n), ddx(0:n), lastcpu
  CHARACTER           :: packorunpack

  BEGIN(mp00ac) 

  if (NOTMatrixFree) then ! if we have the matrix, then just multiply it to x
    call DGEMV('N', n, n, one, dMD(1,1), n+1, x, 1, zero, ddx(1), 1)
    call DGEMV('N', n, n, one, a, n, x, 1, zero, ax, 1)
  else ! if we are matrix-free, then we construct a.x directly
    ideriv = -2        ! this is used for matrix-free only
    packorunpack = 'U'
    WCALL(mp00ac, packab, (packorunpack, vvol, n, x, ideriv) )         ! unpack solution to construct a.x
!   call packab(packorunpack, vvol, n, x, ideriv)          ! unpack solution to construct a.x
    call intghs(Iquad(vvol), mn, vvol, Lrad(vvol), ideriv) ! compute the integrals of B_lower
    call mtrxhs(vvol, mn, Lrad(vvol), dax, ddx, ideriv)    ! construct a.x from the integral
    ax = dax(1:n) - mu * ddx(1:n)                          ! put in the mu factor
  endif

  return

end subroutine matvec

subroutine prec_solve(n,vecin,vecout,au,jau,ju,iperm)
  ! apply the preconditioner
  implicit none
  INTEGER :: n, iperm(*), jau(*), ju(*)
  REAL    :: vecin(*), au(*)
  REAL    :: vecout(*)

  INTEGER :: ii
  REAL :: tempv(n)

  call lusol(n,vecin,tempv,au,jau,ju) ! sparse LU solve
  !  apply permutation
  do ii = 1, n
    vecout(ii) = tempv(iperm(ii+n))
  enddo

  return
end subroutine prec_solve
