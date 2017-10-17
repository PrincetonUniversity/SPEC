!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Calls a NAG routine (either C05NDF or C05PDF) to find a zero of the force-balance vector.

!latex \end{enumerate} \subsubsection{construction of force-balance vector} \begin{enumerate}

!latex \item The force-balance vector is a combination of the Fourier harmonics of the pressure imbalance,
!latex \be F_j \equiv [[p+B^2/2]]_j
!latex \ee
!latex and the spectral condensation constraints,
!latex \be F_j \equiv (R_\t X + Z_\t Y)_j.
!latex \ee

!latex \end{enumerate} \subsubsection{relevant input parameters} \begin{enumerate}

!latex \item \verb+Lfindzero+ : \begin{itemize} \item If \verb+Lfindzero.eq.1+, then \verb+C05NDF+ is used.
!latex                                                This routine uses function values only.
!latex                                          \item If \verb+Lfindzero.eq.2+, then \verb+C05PDF+ is used.
!latex                                                This routine uses function values and user supplied derivatives.
!latex                          \end{itemize}
!latex \item \verb+forcetol+  : \begin{itemize} \item The accuracy in force-balance to which the solution is required.
!latex                          \end{itemize}
!latex \item \verb+c05xtol+   : \begin{itemize} \item The accuracy in position to which the solution is required. See NAG document for more detail.
!latex                          \end{itemize}
!latex \item \verb+c05factor+ : \begin{itemize} \item A quantity to be used in determining the initial step bound. See NAG document for more detail.
!latex                          \end{itemize}
!latex \item \verb+LreadGF+   : \begin{itemize} \item If \verb+LreadGF.eq.T+, then the derivative matrix will be read from the file \verb+.GF+.
!latex                          \end{itemize}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine jk03aa( Ngeometricaldof, position, ic05pxf )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, one, ten

  use numerical, only : sqrtmachprec

  use fileunits, only : ounit

  use inputlist, only : Wmacros, Wjk03aa, ext, &
                        Igeometry, & ! only for screen output; 04 Dec 14;
                        Nvol,                    &
                        Lfindzero, forcetol, c05xtol, c05factor, LreadGF, &
                        ForceErr, Lcheck, Energy

  use cputiming, only : Tjk03aa

  use allglobal, only : writin, myid, ncpu, cpus, &
                        YESstellsym, NOTstellsym, &
                        mn, im, in, iRbc, iZbs, iRbs, iZbc, Mvol, &
                        BBe, IIo, BBo, IIe, &
                        lgeometricaldof, dFFdRZ, dBBdmp, dmupfdx, hessian, dessian, Lhessianallocated 
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in)    :: Ngeometricaldof
  REAL   , intent(inout) :: position(0:Ngeometricaldof)
  INTEGER, intent(out)   :: ic05pxf
  
  LOGICAL                :: LComputeDerivatives
  INTEGER                :: nFcalls, nDcalls, wflag, iflag, idof, jdof, ijdof, ireadhessian, igdof, lvol, ii, imn
  REAL                   :: rflag, lastcpu
  CHARACTER              :: packorunpack

  INTEGER                :: irevcm, mode, Ldfjac, LR
  REAL                   :: xtol, epsfcn, factor
  REAL                   :: diag(1:Ngeometricaldof), RR(1:Ngeometricaldof*(Ngeometricaldof+1)/2), QTF(1:Ngeometricaldof), workspace(1:Ngeometricaldof,1:4)

  REAL                   :: force(0:Ngeometricaldof), fjac(1:Ngeometricaldof,1:Ngeometricaldof)
  
  INTEGER                :: ML, MU ! required for only Lc05ndf;
  
  LOGICAL                :: Lexit = .true. ! perhaps this could be made user input;
  
  BEGIN(jk03aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! screen output; 20 Jun 14;

  if( myid.eq.0 ) then
   cput = GETTIME
   write(ounit,'("jk03aa : ", 10x ," : ")')
   write(ounit,'("jk03aa : ",f10.2," : Lfindzero="i2" ; forcetol="es13.5" ; c05xtol="es13.5" ; c05factor="es13.5" ; LreadGF="L2" ; Ngeometricaldof="i6" ;")')&
                           cput-cpus,  Lfindzero,       forcetol,           c05xtol,           c05factor,           LreadGF,       Ngeometricaldof
   write(ounit,'("jk03aa : ", 10x ," : ")')
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! set some miscellaneous; 20 Jun 14;

  xtol = c05xtol ! tolerance in position supplied to NAG;
  
  Ldfjac = Ngeometricaldof ; LR = Ngeometricaldof * (Ngeometricaldof+1) / 2 ! supplied to NAG;
  
  mode = 0 ; diag(1:Ngeometricaldof) = one ! if mode=2, multiplicative scale factors need to be provided in diag; if mode=0, factors computed internally;
  
  factor = c05factor ! used to determine initial step bound; ! supplied to NAG;
  
  select case( Lfindzero )
  case( 1 )    ; ML = Ngeometricaldof-1 ; MU = Ngeometricaldof-1 ; epsfcn = sqrtmachprec ! only required for C05NDF; ! supplied to NAG;
  case( 2 )    ;
  case default ; FATALMESS(jk03aa, .true., supplied Lfindzero not supported)
  end select
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! initialize some local counters; 20 Jun 14;

  nFcalls = 0 ; nDcalls= 0 ! counters;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  lastcpu = GETTIME
  
  if( Lexit ) then ! will call initial force, and if ForceErr.lt.forcetol will immediately exit; 04 Dec 14;

   LComputeDerivatives= .false.
   WCALL(jk03aa,fc02aa,( Ngeometricaldof, position(0:Ngeometricaldof), force(0:Ngeometricaldof), LComputeDerivatives )) ! calculate the force-imbalance;
   
   if( myid.eq.0 ) then ! screen output; 20 Jun 14;
    cput = GETTIME
    ; write(ounit,1000) cput-cpus, nFcalls, nDcalls, ForceErr,  cput-lastcpu, "|BB|e", alog10(BBe(1:min(Mvol-1,28)))
    if( Igeometry.ge.3 ) then ! include spectral constraints; 04 Dec 14;
     ;write(ounit,1001)                                                                       "|II|o", alog10(IIo(1:min(Mvol-1,28)))
    endif
    if( NOTstellsym ) then
     ;write(ounit,1001)                                                                       "|BB|o", alog10(BBo(1:min(Mvol-1,28)))
     if( Igeometry.ge.3 ) then ! include spectral constraints; 04 Dec 14;
      write(ounit,1001)                                                                       "|II|e", alog10(IIe(1:min(Mvol-1,28)))
     endif
    endif
   endif

   if( ForceErr.lt.forcetol ) then ; ic05pxf = 0 ; goto 9999 ! force-balance is satisfied;  
   endif
   
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

1000 format("jk03aa : ",f10.2," : "i5,i3," ; ":"|f|="es12.5" ; ":"time=",f10.2,"s ;":" log"a5"="28f6.2" ...")
1001 format("jk03aa : ", 10x ," : "5x,3x" ; ":"    "  12x "   ":"     ", 10x ,"  ;":" log"a5"="28f6.2" ...")
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
  irevcm = 0 ; ic05pxf = 1 ! required for initial entry; herefater unchanged by user;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( Lfindzero.eq.2 ) then
   RALLOCATE( dFFdRZ,(1:lgeometricaldof,1:Mvol,0:1,1:lgeometricaldof,0:1))
   RALLOCATE( dBBdmp,(1:lgeometricaldof,1:Mvol,0:1,1:2))
   RALLOCATE( dmupfdx,(1:Mvol,1:2,1:lgeometricaldof,0:1))
   RALLOCATE( hessian,(1:Ngeometricaldof,1:Ngeometricaldof))
   RALLOCATE( dessian,(1:Ngeometricaldof,1:lgeometricaldof))
   Lhessianallocated = .true.
  else
   Lhessianallocated = .false.
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  do ! reverse communication loop; controlled by irevcm;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   select case( Lfindzero )
   
 
   case( 1 ) ! use function values                               to find x st f(x)=0, where x is the geometry of the interfaces, and f is the force;
    
    call C05NDF( irevcm, Ngeometricaldof, position(1:Ngeometricaldof), force(1:Ngeometricaldof), &
                 xtol, ML, MU, epsfcn, diag(1:Ngeometricaldof), mode, factor, fjac(1:Ldfjac,1:Ngeometricaldof), Ldfjac, &
                 RR(1:LR), LR, QTF(1:Ngeometricaldof), workspace(1:Ngeometricaldof,1:4), ic05pxf )
     

   case( 2 ) ! use function values and user-supplied derivatives to find x st f(x)=0, where x is the geometry of the interfaces, and f is the force;
     
    call C05PDF( irevcm, Ngeometricaldof, position(1:Ngeometricaldof), force(1:Ngeometricaldof), &
                                                                              fjac(1:Ldfjac,1:Ngeometricaldof), Ldfjac, &
                 xtol                , diag(1:Ngeometricaldof), mode, factor, &
                 RR(1:LR), LR, QTF(1:Ngeometricaldof), workspace(1:Ngeometricaldof,1:4), ic05pxf )

   case default
    
    FATALMESS( jk03aa, .true., value of Lfindzero not supported)
    
   end select
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   select case( irevcm )
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
   case( 0 ) ! final exit;
    
    if( myid.eq.0 ) then
     cput = GETTIME
     ;              write(ounit,'("jk03aa : ", 10x ," :")')
     select case( ic05pxf )
     case( 0   )  ; write(ounit,'("jk03aa : ",f10.2," : finished ; success        ; ic05p*f="i2" ; its="i7" ,"i4" ;")') cput-cpus, ic05pxf, nFcalls, nDcalls
     case( 1   )  ; write(ounit,'("jk03aa : ",f10.2," : finished ; input error    ; ic05p*f="i2" ; its="i7" ,"i4" ;")') cput-cpus, ic05pxf, nFcalls, nDcalls
     case( 2   )  ; write(ounit,'("jk03aa : ",f10.2," : finished ; irevcm error   ; ic05p*f="i2" ; its="i7" ,"i4" ;")') cput-cpus, ic05pxf, nFcalls, nDcalls
     case( 3   )  ; write(ounit,'("jk03aa : ",f10.2," : finished ; xtol too small ; ic05p*f="i2" ; its="i7" ,"i4" ;")') cput-cpus, ic05pxf, nFcalls, nDcalls
     case( 4:5 )  ; write(ounit,'("jk03aa : ",f10.2," : finished ; bad progress   ; ic05p*f="i2" ; its="i7" ,"i4" ;")') cput-cpus, ic05pxf, nFcalls, nDcalls
     case default ; write(ounit,'("jk03aa : ",f10.2," : finished ; illegal ifail  ; ic05p*f="i2" ; its="i7" ,"i4" ;")') cput-cpus, ic05pxf, nFcalls, nDcalls
     end select
    endif ! end of if( myid.eq.0 ) then;
    
    exit ! escape from infinite do loop;
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
   case( 1 ) ! indicates start of new iteration; no action is required; position and force available for printing; force must not be changed;
    
    packorunpack = 'U' !! unpack geometrical degrees of freedom; 13 Sep 13; ! I guess this is just for writin; 11 Aug 14;
    WCALL(jk03aa,gf00aa,( Ngeometricaldof, position(0:Ngeometricaldof), Mvol, mn, &
 iRbc(1:mn,0:Mvol), iZbs(1:mn,0:Mvol), iRbs(1:mn,0:Mvol), iZbc(1:mn,0:Mvol), packorunpack ))
    
    if( myid.eq.0 ) then
     
     cput = GETTIME
     
     ; write(ounit,1000) cput-cpus, nFcalls, nDcalls, ForceErr, cput-lastcpu, "|BB|e", alog10(BBe(1:min(Mvol-1,28)))
     if( Igeometry.ge.3 ) then ! include spectral constraints; 04 Dec 14;
      ;write(ounit,1001)                                                                      "|II|o", alog10(IIo(1:min(Mvol-1,28)))
     endif
     if( NOTstellsym ) then
      ;write(ounit,1001)                                                                      "|BB|o", alog10(BBo(1:min(Mvol-1,28)))
      if( Igeometry.ge.3 ) then ! include spectral constraints; 04 Dec 14;
       write(ounit,1001)                                                                      "|II|e", alog10(IIe(1:min(Mvol-1,28)))
      endif
     endif
     lastcpu = GETTIME
     
     wflag = 2 ; iflag = nDcalls ; rflag = ForceErr
     WCALL(jk03aa,writin,( wflag, iflag, rflag )) ! write restart file; save geometry to ext.end; iRbc, iZbs consistent with position;
     
    endif ! end of if( myid.eq.0 );
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
   case( 2 ) ! before re-entry to C05NDF / C05PDF, force must contain the function values;
    
    nFcalls = nFcalls + 1
    
    LComputeDerivatives = .false.
    WCALL(jk03aa,fc02aa,( Ngeometricaldof, position(0:Ngeometricaldof), force(0:Ngeometricaldof), LComputeDerivatives )) ! calculate the force-imbalance;
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
   case( 3 ) ! before re-entry to          C05PDF, fjac must contain the derivatives;
    
#ifdef DEBUG
    FATALMESS( jk03aa, .not.Lhessianallocated, need to allocate hessian)
#endif
    
    nDcalls = nDcalls + 1
    
    if( LreadGF .and. nDcalls.eq.1 ) then ! this is the first iteration; will check to see if derivative matrix already exists in file .GF;

     if( myid.eq.0 ) call writereadgf( 'R', Ngeometricaldof, ireadhessian ) ! reads derivatives matrix = hessian from file;
     
     IlBCAST( ireadhessian, 1, 0 )
     
     if( ireadhessian.eq.1 ) then ! hessian has been read from file;
      RlBCAST( hessian(1:Ngeometricaldof,1:Ngeometricaldof), Ngeometricaldof*Ngeometricaldof, 0 )
     endif
     
    else ! matches if( LreadGF .and. nDcalls.eq.1 ) then;
     
     ireadhessian = 0 ! hessian has not been read from file;
     
    endif ! end of if( LreadGF .and. nDcalls.eq.1 ) then;
    
    if( ireadhessian.eq.0 ) then
     
     LComputeDerivatives = .true.
     WCALL(jk03aa,fc02aa,( Ngeometricaldof, position(0:Ngeometricaldof), force(0:Ngeometricaldof), LComputeDerivatives )) ! calculate the force-imbalance;

#ifdef DEBUG
     FATALMESS(jk03aa, Lcheck.eq.4, derivatives of Beltrami field have been computed)
#endif
     
    endif
    
    fjac(1:Ngeometricaldof,1:Ngeometricaldof) = hessian(1:Ngeometricaldof,1:Ngeometricaldof) ! hessian is passed through global; CAN SAVE MEMORY;
    
    if( myid.eq.0 ) call writereadgf( 'W', Ngeometricaldof, ireadhessian ) ! will always save derivative matrix;
    
#ifdef DEBUG

    if( Lcheck.eq.3 ) then
     write(ounit,'("jk03aa : ", 10x ," : myid=",i3," ; volume derivatives have been compared ;")') myid
     stop "jk03aa :            : myid=    ; volume derivatives have been compared ;"
    endif

    FATALMESS(jk03aa, Lcheck.eq.3, volume derivatives have been compared) ! the first process will terminate all processes; 02 Sep 14;

    if( Lcheck.eq.4 ) then
     write(ounit,'("jk03aa : ", 10x ," : myid=",i3," ; field derivatives have been compared ;")') myid
     stop "jk03aa :            : myid=    ; field derivatives have been compared ;"
    endif

    FATALMESS(jk03aa, Lcheck.eq.4, field derivatives have been compared) ! the first process will terminate all processes; 02 Sep 14;

#endif
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
   case default
    
    FATALMESS(jk03aa, .true., illegal irevcm : C05P*F error)
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
   end select ! end of select case(irevcm);
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
  enddo ! end of infinite reverse-communication do loop;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Lfindzero.eq.2 .and. myid.eq.0 .and. irevcm.eq.0 ) then ! will save derivative matrix for future use;

   if( Wjk03aa ) write(ounit,'("jk03aa : ", 10x ," : saving hessian to file ;")')

#ifdef DEBUG
   FATALMESS( jk03aa, .not.Lhessianallocated, error)
#endif

   hessian(1:Ngeometricaldof,1:Ngeometricaldof) = zero
   ijdof = 0
   do idof = 1, Ngeometricaldof
    do jdof = idof, Ngeometricaldof ; ijdof = ijdof + 1 ; hessian(idof,jdof) = RR(ijdof) ! un-pack R matrix;
    enddo
   enddo

!  derivative matrix = Q R;
   hessian(1:Ngeometricaldof,1:Ngeometricaldof) = matmul( fjac(1:Ngeometricaldof,1:Ngeometricaldof), hessian(1:Ngeometricaldof,1:Ngeometricaldof) ) 
   
   call writereadgf( 'W', Ngeometricaldof, ireadhessian ) ! write derivative matrix = hessian to file;

   if( Wjk03aa ) write(ounit,'("jk03aa : ", 10x ," : saved  hessian to file ;")')

  endif ! end of if( myid.eq.0 ) then;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Lfindzero.eq.2 ) then 
   DEALLOCATE( dFFdRZ )
   DEALLOCATE( dBBdmp )
   DEALLOCATE( dmupfdx )
   DEALLOCATE( hessian )
   DEALLOCATE( dessian )
   Lhessianallocated = .false.
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(jk03aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine jk03aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine writereadgf( readorwrite, Ngeometricaldof , ireadhessian )
  
  use constants, only : zero
  
  use numerical, only :
  
  use fileunits, only : ounit, dunit
  
  use inputlist, only : Wjk03aa, ext, Nvol, Mpol, Ntor
  
  use cputiming, only : Tjk03aa
  
  use allglobal, only : myid, cpus, mn, im, in, hessian, Lhessianallocated
  
  LOCALS
  
  CHARACTER, intent(in) :: readorwrite
  LOGICAL               :: exist
  INTEGER, intent(in)   :: Ngeometricaldof
  INTEGER, intent(out)  :: ireadhessian
  
  INTEGER               :: lNvol, lMpol, lNtor, lNgeometricaldof
  
  FATALMESS( jk03aa, .not.Lhessianallocated, error)
  
  ireadhessian = 0 ! set default intent out;
  
  select case( readorwrite ) 
   
  case( 'W' ) ! will write derivative matrix to file;
   
   open( dunit, file="."//trim(ext)//".GF", status="replace", form="unformatted", iostat=ios ) ! save derivative matrix to file;
   FATALMESS( jk03aa, ios.ne.0, error opening derivative matrix file)
   
   write( dunit, iostat=ios ) Nvol, Mpol, Ntor, Ngeometricaldof ! enable resolution consistency check;
   FATALMESS( jk03aa, ios.ne.0, error writing Nvol, Mpol, Ntor, Ngeometricaldof)
   
   write( dunit, iostat=ios ) hessian(1:Ngeometricaldof,1:Ngeometricaldof)
   FATALMESS( jk03aa, ios.ne.0, error writing hessian to file)
   
   close( dunit, iostat=ios )   
   FATALMESS( jk03aa, ios.ne.0, error closing derivative matrix file)
   
  case( 'R' )
   
   cput = GETTIME
   
   inquire( file=".GF", exist=exist ) ! the derivative matrix;

   if( .not.exist ) then ; write(ounit,1000) cput-cpus, myid," does not exist ;                                  " ; goto 9999
   endif
   
   open( dunit, file=".GF", status="old", form="unformatted", iostat=ios)
   
   if( ios  .ne.   0 ) then ; write(ounit,1000) cput-cpus, myid, "error opening file ;                           " ; goto 9999
   endif
   
   read( dunit, iostat=ios ) lNvol, lMpol, lNtor, lNgeometricaldof ! resolution consistency check;
   
   if( ios  .ne.   0 ) then ; write(ounit,1000) cput-cpus, myid, "error reading Nvol, Mpol, Ntor, Ngeometricaldof" ; goto 9998
   endif
   if( lNvol.ne.Nvol ) then ; write(ounit,1000) cput-cpus, myid, "inconsistent resolution ; Nvol ;               ", lNvol, Nvol ; goto 9998
   endif
   if( lMpol.ne.Mpol ) then ; write(ounit,1000) cput-cpus, myid, "inconsistent resolution ; Mpol ;               ", lMpol, Mpol ; goto 9998
   endif
   if( lNtor.ne.Ntor ) then ; write(ounit,1000) cput-cpus, myid, "inconsistent resolution ; Ntor ;               ", lNtor, Ntor ; goto 9998
   endif
   
   read( dunit, iostat=ios ) hessian(1:Ngeometricaldof,1:Ngeometricaldof)
   
   if( ios  .ne.   0 ) then ; write(ounit,1000) cput-cpus, myid, "error reading hessian ;                        " ; goto 9998
   endif
   
   ireadhessian = 1
   
   ;                          write(ounit,1000) cput-cpus, myid, "ok ;                                           "
   
9998 close( dunit, iostat=ios )
   
  case default
   
   FATALMESS( jk03aa, .true., invalid readorwrite)
   
  end select
  
9999 return
  
1000 format("jk03aa : ",f10.2," : myid=",i3," ; .GF file (contains derivative matrix) ; "a37,:"; ":"old="i4" ; new="i4" ;")
  
end subroutine writereadgf

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
