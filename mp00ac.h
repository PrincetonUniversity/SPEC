!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (solver) ! Solves Beltrami/vacuum (linear) system, given matrices.

!latex \briefly{Solves for magnetic vector potential given $\boldmu\equiv(\Delta\psi_t,\Delta\psi_p,\mu)^T$.}

!latex \calledby{\link{ma02aa}}
!latex \calls{\link{packab}, \link{curent} and \link{tr00ab}}

!latex \tableofcontents

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{unpacking fluxes, helicity multiplier}

!latex \begin{enumerate}
!latex \item The vector of ``parameters'', $\boldmu$, is unpackxi. (Recall that $\boldmu$ was ``packxi'' in \link{ma02aa}.)
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

!latex \begin{enumerate}
!latex \item If \inputvar{Lposdef=0}, then it is {\em not} assumed that the linear system is positive definite and
!latex       \nag{www.nag.co.uk/numeric/FL/manual19/pdf/F04/f04aef_fl19.pdf}{F04AEF}
!latex       is used to solve the linear system.
!latex \item If \inputvar{Lposdef=1}, then it {\em is}     assumed that the linear system is positive definite and
!latex       \nag{www.nag.co.uk/numeric/FL/manual19/pdf/F04/f04abf_fl19.pdf}{F04ABF}
!latex       is used to solve the linear system.
!latex \end{enumerate}

!latex \subsection{unpacking, . . .}

!latex \begin{enumerate}
!latex \item The magnetic degrees-of-freedom are ``unpackxi'' by \link{packab}.
!latex \item The error flag, \internal{ImagneticOK}, is set that indicates if the Beltrami fields were successfully constructed.
!latex \end{enumerate}

!l tex \subsection{linear system}

!l tex \begin{enumerate}
!l tex \item The energy, $W \equiv \int \! dv {\; \bf B}\cdot{\bf B}$, and helicity, $K\equiv \int \! dv \; {\bf A}\cdot{\bf B}$, functionals may be written
!l tex \be W & = & \frac{1}{2} \; a_i \; {\cal A}_{i,j} \; a_j + a_i \; {\cal B}_{i,j} \; \psi_j + \frac{1}{2} \; \psi_i \; {\cal C}_{i,j} \; \psi_j \\
!l tex     K & = & \frac{1}{2} \; a_i \; {\cal D}_{i,j} \; a_j + a_i \; {\cal E}_{i,j} \; \psi_j + \frac{1}{2} \; \psi_i \; {\cal F}_{i,j} \; \psi_j 
!l tex \ee
!l tex       where ${\bf a} \equiv \{ \Ate{i,l}, \Aze{i,l}, \Ato{i,l}, \Azo{i,l}, f_{e,i}, f_{o,i} \}$ contains the independent degrees of freedom
!l tex       and $\boldpsi \equiv \{\Delta \psi_t,\Delta \psi_p\}$.
!l tex \item The matrices ${\cal A}_{i,j}$, etc. are constructed in 
!l tex \item The linear system that defines the magnetic field is
!l tex       \be ( {\cal A}_{i,j} - \mu {\cal D}_{i,j} ) \; a_j = - ( {\cal B}_{i,j} - \mu {\cal E}_{i,j} ) \; \psi_j
!l tex       \ee
!l tex \item This routine is written as a function of \mbox{${\bf x}\equiv(\mu,\Delta \psi_p)^T$} for the plasma regions,
!l tex       where $\Delta\psi_t$ is assumed to be given,\\
!l tex       or \mbox{${\bf x}\equiv(\Delta \psi_t,\Delta\psi_p)^T$} for the vacuum regions,
!l tex       where $\mu=0$ is assumed.
!l tex \end{enumerate}

!l tex \subsection{plasma region}

!l tex \begin{enumerate}
!l tex \item If \verb+Lvectorpotential+, then the ``Beltrami'' matrix (see 1.18 of manual.pdf) is constructed from \verb+dMA+ and \verb+dMD+.
!l tex \item If the user has selected \inputvar{Lposdef=1}, then the user has assumed that the Beltrami matrix is positive-definite, 
!l tex       and the NAG routine \nag{www.nag.co.uk/numeric/FL/manual19/pdf/F04/f04abf_fl19.pdf}{F04ABF} is used to solve the linear system,
!l tex       i.e. to solve for the vector-potential in ``packxi'' format.
!l tex \item If the user has selected \inputvar{Lposdef=0}, then the user has assumed that the Beltrami matrix is {\em not} positive-definite, 
!l tex       and so the NAG routine \nag{www.nag.co.uk/numeric/FL/manual19/pdf/F04/f04aef_fl19.pdf}{F04AEF} is used.
!l tex \item If it is known that the Beltrami matrix is positive definite, then it is usually more computationally efficient to exploit this,
!l tex       so the {\em fastest} option is to choose \inputvar{Lposdef=1}; 
!l tex       however, there is no a-priori reason why the Beltrami matrix is positive definite, so the {\em safest} option is to choose \inputvar{Lposdef=0}.
!l tex \end{enumerate}

!l tex \subsection{vacuum region}

!l tex \begin{enumerate}
!l tex \item If \verb+Lvectorpotential+, then the ``Beltrami'' matrix is constructed from \verb+dMA+.
!l tex \item The Beltrami matrix for the vacuum region is identically the Laplacian, and is always positive-definite, 
!l tex       and the NAG routine \nag{www.nag.co.uk/numeric/FL/manual19/pdf/F04/f04abf_fl19.pdf}{F04ABF}
!l tex       is always used to solve the linear system, i.e. to solve for the scalar-potential in ``packxi-format''.
!l tex \item The positive-definitess is guaranteed because the functional $\int B^2 dv$ is guaranteed to have a minimum. [Is this correct?]
!l tex \end{enumerate}

!l tex \subsection{error messages}

!l tex \begin{enumerate}
!l tex \item The \nag{www.nag.co.uk/numeric/FL/manual19/pdf/F04/f04abf_fl19.pdf}{F04ABF}
!l tex       routine may fail if the provided matrix is not positive definite or ill-conditioned, 
!l tex       and such error messages are given as screen output.
!l tex \item The \nag{www.nag.co.uk/numeric/FL/manual19/pdf/F04/f04aef_fl19.pdf}{F04AEF}
!l tex       routine may fail if the provided matrix is singular              or ill-conditioned. 
!l tex \end{enumerate}

!l tex \subsection{energy and helicity integrals}

!l tex \begin{enumerate}
!l tex \item It is most convenient to calculate the energy and helicity integrals with the vector potential in the packxi format, 
!l tex       as these integrals reduce to vector-matrix-vector products.
!l tex \end{enumerate}

!l tex \subsection{``unpacking''}

!l tex \begin{enumerate}
!l tex \item The routine \link{packab} is then called to ``unpack'' the linear solution into the more easily identifiable arrays.
!l tex \end{enumerate}

!l tex \subsection{returned values}

!l tex \begin{enumerate}
!l tex \item This routine may be called iteratively to find the particular $\boldmu$ that satisfies the required constraints, 
!l tex       e.g. for \inputvar{Lconstraint=1}, the user wishes to enforce the rotational-transform constraint on the adjacent interfaces, 
!l tex       which is calculated by calling \link{tr00ab}.
!l tex \item Because this routine is called by NAG the input/ouput arguments are constrained, and \link{mp00ac} returns either a function,
!l tex       which is equal to zero when the appropriate constraints are satisfied, or the derivative of the same function 
!l tex       with respect to $\Delta \psi_p$ and/or $\mu$, as the case may be.
!l tex \item Note that the derivatives of the function are determined by matrix perturbation methods.
!l tex       If the derivatives are required, then the derivatives of the Beltrami matrix must be provided in the \verb+dMA+ and \verb+dMD+ arrays.
!l tex \end{enumerate}
 
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
!latex \item If $|{\bf f}| < $ \inputvar{mupftol}, then early termination is enforced (i.e. \internal{iflag} is set to negative integer).
!latex       (See \link{ma02aa} for details of how \link{mp00ac} is called iteratively.)
!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine mp00ac( Ndof, Xdof, Fdof, Ddof, Ldfjac, iflag ) ! argument list is fixed by NAG; ma02aa calls mp00ac through C05PCF;

! if iflag.eq.0 : Xdof and Fdof are available for PRINTING ; 28 Jan 13; Fdof MUST NOT BE CHANGED; Ddof MUST NOT BE CHANGED;
! if iflag.eq.1 :          Fdof is to be          UPDATED  ; 28 Jan 13;                         ; Ddof MUST NOT BE CHANGED;
! if iflag.eq.2 :          Ddof is to be          UPDATED  ; 28 Jan 13; Fdof MUST NOT BE CHANGED;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half, one, two, goldenmean
  
  use numerical, only : machprec, sqrtmachprec, small
  
  use fileunits, only : ounit
  
  use inputlist, only : Wmacros, Wmp00ac, Wtr00ab, Wcurent, Wma02aa, &
                        mu, helicity, iota, oita, curtor, curpol, Lrad, &
                        Lposdef, &
                        Lconstraint, mupftol
  
  use cputiming, only : Tmp00ac
  
  use allglobal, only : myid, ncpu, cpus, ivol, &
                        YESstellsym, NOTstellsym, &
                        Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, &
                        mn, im, in, mns, &
                        Nt, Nz, & ! only required to pass through as arguments to tr00ab; 23 Apr 13;
                        NAdof, &
                        dMA, dMB, dMC, dMD, dME, dMF, dMG, &
                        solution, &
                        dtflux, dpflux, &
                        diotadxup, dItGpdxtp, &
                        lBBintegral, lABintegral, &
                        xoffset, &
                        ImagneticOK, &
                        Ate, Aze, Ato, Azo, Mvol
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in)  :: Ndof, Ldfjac
  REAL   , intent(in)  :: Xdof(1:Ndof)
  REAL                 :: Fdof(1:Ndof), Ddof(1:Ldfjac,1:Ndof)
  INTEGER              :: iflag ! indicates whether (i) iflag=1: ``function'' values are required; or (ii) iflag=2: ``derivative'' values are required;


  INTEGER              :: lvol, NN, MM, ideriv, IA, IB, IC, IBB, IAA, lmns, if04abf(0:1), if04aef(0:1), ii, jj, nnz
  
  REAL                 :: lmu, dpf, dtf, dpsi(1:2), tpsi(1:2), ppsi(1:2), lcpu !, icurrent(0:2), gcurrent(0:2) ! 12 Sep 16;
  
  CHARACTER            :: packorunpack
  
  REAL   , allocatable :: matrix(:,:), rhs(:,:)

  REAL   , allocatable :: RW(:), RD(:,:), LU(:,:)
  
  BEGIN(mp00ac)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  lvol = ivol ! recall that ivol is global; 24 Jan 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  FATAL( mp00ac, iflag.ne.1 .and. iflag.ne.2, invalid iflag ) ! see nprint=0 in ma02aa and C05PCF; 12 Sep 16;
#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Lplasmaregion ) then
   
   ;                    ; lmu  = Xdof(1) - xoffset
   ;                    ; dtf  = dtflux(lvol)
   if( Ndof.eq.2 ) then ; dpf  = Xdof(2) - xoffset
   else                 ; dpf  = dpflux(lvol)
   endif
   
  else ! Lvacuumregion; 26 Jan 16;
   
#ifdef FORCEFREEVACUUM
   ;                    ; lmu  = mu(lvol)           ! generalize for arbitrary force-free field; 04 May 17;
#else
   ;                    ; lmu  = zero               ! restrict attention to strict vacuum field; 04 May 17;
#endif

   ;                    ; dtf  = Xdof(1) - xoffset
   if( Ndof.eq.2 ) then ; dpf  = Xdof(2) - xoffset
   else                 ; dpf  = dpflux(lvol)
   endif
      
  endif ! end of if( Lplasmaregion ) ; 11 Mar 16;

  dpsi(1:2) = (/  dtf,  dpf /) ! enclosed toroidal fluxes and their derivatives;
  tpsi(1:2) = (/  one, zero /) ! enclosed toroidal fluxes and their derivatives;
  ppsi(1:2) = (/ zero,  one /) ! enclosed toroidal fluxes and their derivatives;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  diotadxup(0:1,-1:2,lvol) = zero ! rotational-transform, and its derivatives with respect to lmu and dpf, or toroidal current, on the inner/outer interface;
  dItGpdxtp(0:1,-1:2,lvol) = zero ! plasma and linking currents;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  NN = NAdof(lvol) ! shorthand;
  
  SALLOCATE( matrix, (1:NN,1:NN), zero )
  SALLOCATE( rhs   , (1:NN,0:2 ), zero )

  solution(1:NN,-1:2) = zero ! this is a global array allocated in dforce; 20 Jun 14;
  
  SALLOCATE( RW, (1:NN    ), zero )
  SALLOCATE( RD, (1:NN,0:2), zero )

  if( Lposdef.eq.0 ) then
  SALLOCATE( LU, (1:NN,1:NN), zero )
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if04abf(0:1) = 0 ; if04aef(0:1) = 0 ! error flags;  4 Feb 13;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! ideriv labels derivative as follows:
! ideriv = 0 : compute Beltrami field; ideriv = 1 : compute d Beltrami field / d \mu           ; ideriv = 2 : compute d Beltrami field / d \Delta \psi_p ;
! ideriv = 0 : compute Vacuum   field; ideriv = 1 : compute d Vacuum   field / d \Delta \psi_t ; ideriv = 2 : compute d Vacuum   field / d \Delta \psi_p ;
  
  do ideriv = 0, 1 ! loop over derivatives;
   
   if( iflag.eq.1 .and. ideriv.eq.1 ) cycle ! only need to return function; recall the derivative estimate requires function evaluation;
   
   if( Lcoordinatesingularity ) then
    
    matrix(1:NN,1:NN) = dMA(1:NN,1:NN) - lmu * dMD(1:NN,1:NN)
    
    select case( ideriv )
    case( 0 )    ; rhs(1:NN,0) = - matmul(  dMB(1:NN,1:2) - lmu  * dME(1:NN,1:2), dpsi(1:2) )
    case( 1 )    ; rhs(1:NN,1) = - matmul(                - one  * dME(1:NN,1:2), dpsi(1:2) ) - matmul( - one  * dMD(1:NN,1:NN), solution(1:NN,0) )
     ;           ; rhs(1:NN,2) = - matmul(  dMB(1:NN,1:2) - lmu  * dME(1:NN,1:2), ppsi(1:2) )
    end select
    
   else ! .not.Lcoordinatesingularity; 
    
    if( Lplasmaregion ) then
     
     matrix(1:NN,1:NN) = dMA(1:NN,1:NN) - lmu * dMD(1:NN,1:NN)
     
     select case( ideriv )
     case( 0 )    ; rhs(1:NN,0) = - matmul( dMB(1:NN,1:2 ), dpsi(1:2) )
     case( 1 )    ; rhs(1:NN,1) =                                       - matmul( - one * dMD(1:NN,1:NN), solution(1:NN,0) )
      ;           ; rhs(1:NN,2) = - matmul( dMB(1:NN,1:2 ), ppsi(1:2) )
     end select
     
    else ! Lvacuumregion ; 08 Feb 16;
     
#ifdef FORCEFREEVACUUM
     FATAL( mp00ac, .true., need to revise Beltrami matrices in vacuum region for arbitrary force-free field )
#else
     matrix(1:NN,1:NN) = dMA(1:NN,1:NN) ! - lmu * dMD(1:NN,1:NN) ; 04 May 17;

     select case( ideriv )
     case( 0 )    ; rhs(1:NN,0) = - dMG(1:NN) - matmul( dMB(1:NN,1:2), dpsi(1:2) ) ! perhaps there is an lmu term missing here; 04 May 17;
     case( 1 )    ; rhs(1:NN,1) =             - matmul( dMB(1:NN,1:2), tpsi(1:2) ) ! perhaps there is an lmu term missing here; 04 May 17;
      ;           ; rhs(1:NN,2) =             - matmul( dMB(1:NN,1:2), ppsi(1:2) ) ! perhaps there is an lmu term missing here; 04 May 17;
     end select
#endif

    endif ! end of if( Lplasmaregion ) ; 08 Feb 16;

   endif ! end of if( Lcoordinatesingularity ) ; 08 Feb 16;
   

   IA = NN ; IB = NN ; IC = NN ; IBB = NN ; IAA = NN
      
   lcpu = GETTIME
   
   
   select case( Lposdef )
    
   case( 0 ) ! Lposdef=0;
    
    select case( ideriv )
     
    case( 0 ) ! Lposdef=0; ideriv=0;

     if04aef(ideriv) = 1 ; MM = 1
     call F04AEF( matrix(1:NN,1:NN), IA, rhs(1:NN,  0), IB, NN, MM, solution(1:NN,  0), IC, RW(1:NN), LU(1:NN,1:NN), IAA, RD(1:NN,  0), IBB, if04aef(ideriv) )

    case( 1 ) ! Lposdef=0; ideriv=1;

     if04aef(ideriv) = 1 ; MM = 2
     call F04AEF( matrix(1:NN,1:NN), IA, rhs(1:NN,1:2), IB, NN, MM, solution(1:NN,1:2), IC, RW(1:NN), LU(1:NN,1:NN), IAA, RD(1:NN,1:2), IBB, if04aef(ideriv) )

    end select ! ideriv;
    
    cput = GETTIME

    select case( if04aef(ideriv) )                                                                            !123456789012345678
    case(  0  )  ; if( Wmp00ac ) write(ounit,1010) cput-cpus, myid, lvol, ideriv, "if04aef", if04aef(ideriv), "success ;         ", cput-lcpu
    case(  1  )  ;               write(ounit,1010) cput-cpus, myid, lvol, ideriv, "if04aef", if04aef(ideriv), "singular ;        "
    case(  2  )  ;               write(ounit,1010) cput-cpus, myid, lvol, ideriv, "if04aef", if04aef(ideriv), "ill conditioned ; "
    case(  3  )  ;               write(ounit,1010) cput-cpus, myid, lvol, ideriv, "if04aef", if04aef(ideriv), "input error ;     "
    case default ;               write(ounit,1010) cput-cpus, myid, lvol, ideriv, "if04aef", if04aef(ideriv), "invalid if04aef ; "
    end select
    
   case( 1 ) ! Lposdef=1;
    
    select case( ideriv )
     
    case( 0 ) ! Lposdef=1; ideriv=0;
     
     if04abf(ideriv) = 1 ; MM = 1
     call F04ABF( matrix(1:IA,1:NN), IA, rhs(1:IB,0:0), IB, NN, MM, solution(1:IC,0:0), IC, RW(1:NN), RD(1:IBB,0:0), IBB, if04abf(ideriv) )
     
    case( 1 ) ! Lposdef=1; ideriv=1;
     
     if04abf(ideriv) = 1 ; MM = 2
     call F04ABF( matrix(1:IA,1:NN), IA, rhs(1:NN,1:2), IB, NN, MM, solution(1:NN,1:2), IC, RW(1:NN), RD(1:IBB,1:2), IBB, if04abf(ideriv) )
     
    end select ! ideriv;
    
    cput = GETTIME
    
    select case( if04abf(ideriv) ) !                                                                           1234567890123456789012345678901234
    case(  0  )  ; if( Wmp00ac ) write(ounit,1010) cput-cpus, myid, lvol, ideriv, "if04abf", if04abf(ideriv), "success ;                         ", cput-lcpu
    case(  1  )  ;               write(ounit,1010) cput-cpus, myid, lvol, ideriv, "if04abf", if04abf(ideriv), "not +ve definite ; try Lposdef=0 ;"
    case(  2  )  ;               write(ounit,1010) cput-cpus, myid, lvol, ideriv, "if04abf", if04abf(ideriv), "ill conditioned ;                 "
    case(  3  )  ;               write(ounit,1010) cput-cpus, myid, lvol, ideriv, "if04abf", if04abf(ideriv), "input error ;                     "
    case default ;               write(ounit,1010) cput-cpus, myid, lvol, ideriv, "if04abf", if04abf(ideriv), "invalid if04abf ;                 "
    end select

   end select ! Lposdef; 17 Dec 15;
   
1010 format("mp00ac : ",f10.2," : myid=",i3," ; lvol=",i3," ; ideriv="i2" ; "a7"=",i3," ; "a34,:" time=",f10.2," ;")
   
  enddo ! end of do ideriv; 25 Jan 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! can compute the energy and helicity integrals; easiest to do this with solution in packxi format; 20 Feb 13;
  
   lBBintegral(lvol) = half * sum( solution(1:NN,0) * matmul( dMA(1:NN,1:NN), solution(1:NN,0) ) ) & 
                     +        sum( solution(1:NN,0) * matmul( dMB(1:NN,1: 2),     dpsi(1: 2  ) ) ) &
                     + half * sum(     dpsi(1: 2  ) * matmul( dMC(1: 2,1: 2),     dpsi(1: 2  ) ) )
  
   lABintegral(lvol) = half * sum( solution(1:NN,0) * matmul( dMD(1:NN,1:NN), solution(1:NN,0) ) ) & 
                     +        sum( solution(1:NN,0) * matmul( dME(1:NN,1: 2),     dpsi(1: 2  ) ) ) &
                     + half * sum(     dpsi(1: 2  ) * matmul( dMF(1: 2,1: 2),     dpsi(1: 2  ) ) )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  do ideriv = 0, 2
   
   if( iflag.eq.1 .and. ideriv.gt.0 ) cycle
   
   packorunpack = 'U'
   WCALL( mp00ac, packab, ( packorunpack, lvol, NN, solution(1:NN,ideriv), ideriv ) ) ! unpacking; this assigns oAt, oAz through common;
   
  enddo ! do ideriv = 0, 2; 12 Sep 16;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  DALLOCATE( matrix )
  DALLOCATE( rhs    )
  
  DALLOCATE( RW )
  DALLOCATE( RD )

  if( Lposdef.eq.0 ) then
  DALLOCATE( LU )
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( if04abf(0).ne.0 .or. if04aef(0).ne.0 .or. if04abf(1).ne.0 .or. if04aef(1).ne.0 ) then ! failed to construct Beltrami/vacuum field and/or derivatives;
   
   ImagneticOK(lvol) = .false. ! set error flag;
   
   if( iflag.eq.1 ) Fdof(1:Ndof       ) = zero ! provide dummy intent out;
   if( iflag.eq.2 ) Ddof(1:Ndof,1:Ndof) = zero ! provide dummy intent out;
   
   iflag = -1 ! this value will be returned by C05PCF to ma02aa;
   
   goto 9999
   
  else
   
   ImagneticOK(lvol) = .true. ! set error flag; used in dforce;
  
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( YESstellsym ) then ; lmns = 1 + (mns-1)           ! number of independent degrees of freedom in angle transformation; 30 Jan 13; 
  else                   ; lmns = 1 + (mns-1) + (mns-1) ! only required for dense, Fourier angle transformation; 21 Apr 13;
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!#ifdef NEWIOTA ! 19 Sep 16;
!! perhaps a rigid shift in the angle does not change the rotational-transform; 02 Sep 14; ! 19 Sep 16;
!  if( NOTstellsym ) lmns = 1 + (mns-1) + (mns  ) ! included non-stellarator symmetric angle transformation; 02 Sep 14; ! 19 Sep 16;
!#else ! 19 Sep 16;
!  if( NOTstellsym ) lmns = 1 + (mns-1) + (mns-1) ! only required for dense, Fourier angle transformation; 21 Apr 13;
!#endif ! 19 Sep 16;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  select case( Lconstraint ) 
   
  case( -1 ) ! Lconstraint=-1; 30 Jan 13;
   
   if( Lplasmaregion ) then
    
    if( Wtr00ab ) then ! compute rotational transform only for diagnostic purposes; 21 Apr 13;
     WCALL( mp00ac, tr00ab, ( lvol, mn, lmns, Nt, Nz, iflag, diotadxup(0:1,-1:2,lvol) ) )
    endif
    
    Fdof(1:Ndof       ) = zero ! provide dummy intent out; Lconstraint=-1 indicates no iterations over mu   , dpflux are required;
    Ddof(1:Ndof,1:Ndof) = zero ! provide dummy intent out;   
    
   else ! Lvacuumregion
    
    if( Wtr00ab ) then ! compute rotational transform only for diagnostic purposes; 21 Apr 13;
     WCALL( mp00ac, tr00ab, ( lvol, mn, lmns, Nt, Nz, iflag, diotadxup(0:1,-1:2,lvol) ) )
    endif
    
    if( Wcurent ) then ! compute enclosed currents    only for diagnostic purposes; 21 Apr 13;
     WCALL( mp00ac, curent,( lvol, mn, Nt, Nz, iflag, dItGpdxtp(0:1,-1:2,lvol) ) )     ! 12 Sep 16;
     curtor = dItGpdxtp(0,0,lvol) ! icurrent(0) ! update input variables; 08 Jun 16;  ! 12 Sep 16;
     curpol = dItGpdxtp(1,0,lvol) ! gcurrent(0)
    endif
    
    Fdof(1:Ndof       ) = zero ! provide dummy intent out;Lconstraint=-1 indicates no iterations over dtflux, dpflux are required;
    Ddof(1:Ndof,1:Ndof) = zero ! provide dummy intent out;   
    
   endif ! end of if( Lplasmaregion) ; 26 Jan 16;
   
  case(  0 ) ! Lconstraint= 0; 30 Jan 13;
   
   if( Lplasmaregion ) then
    
    if( Wtr00ab ) then ! compute rotational transform only for diagnostic purposes; 21 Apr 13;
     WCALL( mp00ac, tr00ab, ( lvol, mn, lmns, Nt, Nz, iflag, diotadxup(0:1,-1:2,lvol) ) )
    endif
    
    Fdof(1:Ndof       ) = zero ! provide dummy intent out; Lconstraint= 0 indicates no iterations over mu, dpflux are required;
    Ddof(1:Ndof,1:Ndof) = zero ! provide dummy intent out;   
    
   else ! Lvacuumregion
    
    WCALL( mp00ac, curent,( lvol, mn, Nt, Nz, iflag, dItGpdxtp(0:1,-1:2,lvol) ) ) ! 12 Sep 16;
    
    if( iflag.eq.1 ) Fdof(1:2  ) = (/ dItGpdxtp(0,0,lvol) - curtor, dItGpdxtp(1,0,lvol) - curpol /)
    if( iflag.eq.2 ) Ddof(1:2,1) = (/ dItGpdxtp(0,1,lvol)         , dItGpdxtp(1,1,lvol)          /)
   !if( iflag.eq.2 ) Ddof(1:2,2) = (/ dItGpdxtp(0,0,lvol)         , dItGpdxtp(1,2,lvol)          /) ! 19 Sep 16;
    if( iflag.eq.2 ) Ddof(1:2,2) = (/ dItGpdxtp(0,2,lvol)         , dItGpdxtp(1,2,lvol)          /) ! 19 Sep 16;
    
   endif ! end of if( Lplasmaregion) ; 26 Jan 16;
   
  case(  1 ) ! Lconstraint= 1; 30 Jan 13;
   
   WCALL( mp00ac, tr00ab,( lvol, mn, lmns, Nt, Nz, iflag, diotadxup(0:1,-1:2,lvol) ) ) ! required for both plasma and vacuum region; 12 Sep 16;
   
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

    WCALL( mp00ac, curent, ( lvol, mn,     Nt, Nz, iflag, dItGpdxtp(0:1,-1:2,lvol) ) ) ! 12 Sep 16;

    curtor = dItGpdxtp(0,0,lvol) ! update input variables; 08 Jun 16; 
   !curpol = dItGpdxtp(1,0,lvol)

    if( iflag.eq.1 ) Fdof(1:2  ) = (/ diotadxup(0,0,lvol) - oita(lvol-1), dItGpdxtp(1,0,lvol) - curpol /)
    if( iflag.eq.2 ) Ddof(1:2,1) = (/ diotadxup(0,1,lvol)               , dItGpdxtp(1,1,lvol)          /)
    if( iflag.eq.2 ) Ddof(1:2,2) = (/ diotadxup(0,2,lvol)               , dItGpdxtp(1,2,lvol)          /)

   endif ! end of if( Lplasmaregion) ; 26 Jan 16;

  case(  2 )

   FATAL( mp00ac, .true., where is helicity calculated )

  end select ! end of select case( Lconstraint ) ; 08 Feb 16;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Wmp00ac .or. Wma02aa ) then ! the following is screen output; 03 Apr 13; ! 04 Dec 14;
   
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
  
  endif ! end of if( Wmp00ac .or. Wma02aa ) ; 21 Apr 13;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( iflag.eq.1 ) then ! only in this case is Fdof defined; 17 Apr 13;
   
   if( sum( abs( Fdof(1:Ndof) ) ) / Ndof .lt. mupftol ) then ! satisfactory;  1 Feb 13;
    
    if ( Lplasmaregion ) then ; mu(lvol) = lmu  ;                    ; dpflux(lvol) = dpf
#ifdef FORCEFREEVACUUM
    else                      ; mu(lvol) = lmu  ; dtflux(lvol) = dtf ; dpflux(lvol) = dpf ! 04 May 17;
#else
    else                      ; mu(lvol) = zero ; dtflux(lvol) = dtf ; dpflux(lvol) = dpf
#endif
    endif
    
    iflag = -2 ! return "acceptance" flag through to ma02aa via ifail;  1 Feb 13; early termination; 
    
   endif ! end of if( sum(Fdof) ) ; 17 Apr 13;
   
  endif ! end of if( iflag.eq.1 ) ; 17 Apr 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(mp00ac)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine mp00ac

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
