!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Locates periodic orbits and their residues

!latex \item This routine uses a simple trajectory following method to locate fieldlines that satisfy
!latex       \be \s( \zeta + 2\pi q) &=& \s(\zeta) \\
!latex           \t( \zeta + 2\pi q) &=& \t(\zeta)+2\pi p,
!latex       \ee
!latex       where the integers $p$ and $q$ are given on input.

!latex \item The tangent map, $M$, is exploited.
!latex       This is defined by
!latex       \be \bar \s( \s+\delta \s, \t+\delta \t) =\bar \s( \s, \t) + M \cdot \left(\delta \s, \delta \t\right)^T,
!latex       \ee
!latex       where $\bar \s \equiv \s( \zeta + 2\pi q)$ and $\bar \t \equiv \t( \zeta + 2\pi q)$.

!latex \item The field and the tangent map are provided by \verb!bf00aa!.

!latex \item The residue is related to the tangent map, and provides information regarding the stability of the orbit \cite{Greene_79}.

!latex \item More details of the algorithm are supplied in \cite{Hudson_04}.
!latex        To locate high-order periodic orbits in strongly chaotic fields, Lagrangian integration will usually be preferred \cite{Hudson_06a}.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \end{enumerate} \subsubsection{input parameters} \begin{enumerate}

!latex \item \verb+odetol+ : o.d.e. integration accuracy provided to NAG routine \verb+D02BJF+;
!latex \item \verb!Mpqits! : maximum iterations allowed in search;
!latex \item \verb!Lpqsym! : controls whether search is constrained along symmetry line, $\t=0$;
!latex \item \verb!pqs, pqt! : initial guess for the $(\s,\t)$; only relevant if \verb+Lpqsym=0+;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

module periodicitylength

  REAL :: zstart
  
end module periodicitylength

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine pq00aa( lvol, nc, low, upp, ipqfail )

  use constants, only : zero, one, two, four, pi, pi2

  use numerical, only : small, sqrtmachprec

  use fileunits, only : ounit

  use inputlist, only : Wpq00aa, Nfp, odetol, Mpqits, Lpqsym

  use cputiming, only : Tpq00aa

  use allglobal, only : myid, cpus, pi2nfp, ivol, pqorbit, Ltangent

  use periodicitylength

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS

  INTEGER, intent(in)    :: lvol, nc
  REAL   , intent(in)    :: low, upp                         ! lower/upper bounds for search;
  INTEGER, intent(inout) :: ipqfail                          ! error flag;

  INTEGER                :: Node, ii, ifail, ioddeven, iev, its, pp, qq, kk
  REAL                   :: zst, zend, tol, Tsti(1:6), Tst(1:6), TWk(20*6), serr, terr, snew, tnew, determinant, harvest, TangentMap(1:2,1:2)
  REAL                   :: stz(1:3), RpZ(1:3), dR(0:3), dZ(0:3) , jacobian, guv(1:3,1:3) ! coordinate tranformation;
  CHARACTER              :: RR, failed*6

  INTEGER                :: Nmat, LDA, LDVR, LDVI, Lwk
  REAL                   :: wr(1:2), wi(1:2), vr(1:2,1:2), vi(1:2,1:2), wk(1:8)
  CHARACTER*1            :: job
  
  EXTERNAL               :: bf00aa, Tbfout, Tbfend
  
  BEGIN(pq00aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  ivol = lvol ; Ltangent = 1 ! this needs to be passed through to bf00aa;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  its = 0
  
  pp = pqorbit(lvol,nc)%pq(1)
  qq = pqorbit(lvol,nc)%pq(2)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  Tsti(1:6) = (/ pqorbit(lvol,nc)%so, pqorbit(lvol,nc)%to, one, zero, zero, one /) ! INITIAL GUESS; 26 Feb 13;
  
  pqorbit(lvol,nc)%ok = 1 ! CAREFUL: default intent out is to indicate that periodic orbit has been correctly located;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \end{enumerate} \subsubsection{symmetry line} \begin{enumerate}

!latex \item The default is to assume stellarator symmetry, and only periodic orbits along $\theta=0$ line are located. 
!latex \item In this case, there is only 1 degree-of-freedom in the numerical search, which is the radial location, $s$.
!latex \item The condition that the trajectory is periodic is
!latex       \be \theta( 2 \pi q)=\theta(0)+2 \pi p,
!latex       \ee
!latex       which must be satisfied to within \verb!q*odetol!.
!latex \item Note that stellarator symmetric orbits need only be followed {\em half} of the full periodicity distance, so that the periodicity condition reduces to
!latex       \be \theta(   \pi q)=\theta(0)+  \pi p.
!latex       \ee
!latex       This is under construction.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  Node = 6 ; tol = odetol ; RR = 'D' ; TWk(1:120) = zero ; terr = one ; serr = one

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
  if( Tsti(1).le.low .or. Tsti(1).ge.upp ) then

   failed='failed'

   cput = GETTIME
   write(ounit,1002)cput-cpus, myid, lvol, pp, qq, failed, its, terr, serr, Tsti(2), low, Tsti(1), upp ; goto 9999

  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  do ; its = its + 1 ! iterate until periodic orbit is located;
   
   select case( Lpqsym ) ! controls whether search is general or whether search is along symmetry line, \t = zero;
   case( 0 ) ! Tsti(1:6) = (/ Tsti(1), Tsti(2), one, zero, zero, one /) ! search is general        ; initial guess provided elsewhere;
   case( 1 ) ; Tsti(1:6) = (/ Tsti(1),  zero  , one, zero, zero, one /) ! search is along \t = zero; initial guess for angle is provided in Tsti(2) is irrelevant;
   case default
    FATALMESS(pq00aa, .true., invalid Lpqsym )
   end select
   
   Tst(1:6) = Tsti(1:6)
   
   pqorbit(lvol,nc)%s(0) = Tst(1) ! this information is returned (elsewhere, a smooth curve is fit through the periodic orbits);
   pqorbit(lvol,nc)%t(0) = Tst(2)
   
   do ii = 1, qq ! integrate along magnetic fieldlines;
   
    zstart = (ii-1) * pi2nfp ! this is common; used in Tbfend to terminate o.d.e. integration;

    zst    = (ii-1) * pi2nfp
    zend   =  ii    * pi2nfp + sqrtmachprec

    ifail = 1
    call D02BJF( zst, zend, Node, Tst(1:Node), bf00aa, tol, RR, Tbfout, Tbfend, TWk(1:120), ifail ) ! NAG ode integration;

    cput = GETTIME
    select case( ifail )
    case( 0 ) ; 
    case( 1 ) ; write(ounit,'("pq00aa : ",f10.2," : myid=",i3," : ifail=",i3," : input error ;                 ")')cput-cpus,myid,ifail
    case( 2 ) ; write(ounit,'("pq00aa : ",f10.2," : myid=",i3," : ifail=",i3," : no further progress possible ;")')cput-cpus,myid,ifail
    case( 3 ) ; write(ounit,'("pq00aa : ",f10.2," : myid=",i3," : ifail=",i3," : tol too small ;               ")')cput-cpus,myid,ifail
    case( 4 ) ; write(ounit,'("pq00aa : ",f10.2," : myid=",i3," : ifail=",i3," : xsol not reset ;              ")')cput-cpus,myid,ifail
    case( 5 ) ; write(ounit,'("pq00aa : ",f10.2," : myid=",i3," : ifail=",i3," : xsol not reset ;              ")')cput-cpus,myid,ifail
    case( 6 ) ; write(ounit,'("pq00aa : ",f10.2," : myid=",i3," : ifail=",i3," : function did not change sign ;")')cput-cpus,myid,ifail
    case( 7 ) ; write(ounit,'("pq00aa : ",f10.2," : myid=",i3," : ifail=",i3," : serious error ;               ")')cput-cpus,myid,ifail
    case default
     FATALMESS(pq00aa,.true.,invalid ifail returned from D03BJF : error integrating field to locate periodic orbit)
    end select

    pqorbit(lvol,nc)%s(ii) = Tst(1) ! this information is returned (elsewhere, a smooth curve is fit through the periodic orbits);
    pqorbit(lvol,nc)%t(ii) = Tst(2)
    
   enddo ! end of do i = 1, qq;
   
   serr = ( Tsti(1)          ) - Tst(1)
   terr = ( Tsti(2) + pp*pi2 ) - Tst(2)

!  !half-period --> full-period tangent map;
  
!   tm(1,1)=Tst(6) ; tm(1,2)=Tst(5) ; tm(2,1)=Tst(4) ; tm(2,2)=Tst(3) ; it=tm ; it(1,1)=tm(2,2) ; it(2,2)=tm(1,1)
!   tm=matmul(it,tm) 
!   determinant=tm(1,1)*tm(2,2)-tm(1,2)*tm(2,1) ; residue=(two-tm(1,1)-tm(2,2))/four 
   
   determinant = Tst(3)*Tst(6) - Tst(4)*Tst(5) ! serves as check on calculation; determinant of tangent map at periodic orbit should be equal to one;
   
   pqorbit(lvol,nc)%residue = ( two - (Tst(3)+Tst(6)) ) / four ! IS THIS CORRECT WHEN NFP > 1 ; 04 Dec 12;
   
   if( Wpq00aa ) then
    cput = GETTIME ; failed='      '
    write(ounit,1002)cput-cpus, myid, lvol, pp, qq, failed, its, terr, serr, Tsti(2), low, Tsti(1), upp, abs(determinant-one), pqorbit(lvol,nc)%residue
   endif
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   if( its.ge.Mpqits ) then ; pqorbit(lvol,nc)%ok = 0 ; exit ! failed;
   endif

   select case( Lpqsym )
   case( 0 ) ; pqorbit(lvol,nc)%error = abs(terr)+abs(serr) !  two   degrees of freedom and  two   constraints; this should really be pqtol?
   case( 1 ) ; pqorbit(lvol,nc)%error = abs(terr)           ! single degree  of freedom and single constraint ; this should really be pqtol?
   case default
    FATALMESS(pq00aa,.true.,invalid Lpqsym)
   end select
   
   if( pqorbit(lvol,nc)%error.lt.qq*odetol ) exit ! finished searching for periodic orbit; recall that pqorbit(lvol,nc)%ok = 1 was set earlier;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
!latex \end{enumerate} \subsubsection{algorithm details} \begin{enumerate}
   
!latex \item A standard Newton search is performed to locate periodic orbits. 
!latex       \be \left( \begin{array}{c} \bar \s(\s+ \delta s, \t+\delta \t ) \\ \bar \t( \s+ \delta s, \t+\delta \t )
!latex                  \end{array}
!latex           \right)
!latex         = \left( \begin{array}{c} \bar \s(\s          , \t           ) \\ \bar \t( \s          , \t           )
!latex                  \end{array}
!latex           \right)
!latex         + \left( \begin{array}{ccc} \partial_\s \bar \s &,& \partial_\t \bar \s \\
!latex                                     \partial_\s \bar \t &,& \partial_\t \bar \t
!latex                  \end{array}
!latex           \right)
!latex           \left( \begin{array}{c}             \delta \s                \\                           \delta \t
!latex                  \end{array}
!latex           \right)
!latex         = \left( \begin{array}{l}         \s+ \delta s                 \\                        \t+\delta \t + 2 \pi p
!latex                  \end{array}
!latex           \right)
!latex       \ee
   
   TangentMap(1,1) = Tst(3) - one ! \partial_\s \bar \s - 1 ; note that this is a `modified' tangent map required for Newton iterations;
   TangentMap(1,2) = Tst(4)       ! \partial_\t \bar \s
   TangentMap(2,1) = Tst(5)       ! \partial_\s \bar \t
   TangentMap(2,2) = Tst(6) - one ! \partial_\t \bar \t - 1
   
   determinant = ( TangentMap(1,1) * TangentMap(2,2) - TangentMap(1,2) * TangentMap(2,1) )

   FATALMESS(pq00aa, abs(determinant).lt.small, error constructing tangent map )

   select case( Lpqsym )

   case( 0 ) ; snew = Tsti(1) + (   TangentMap(2,2) * serr - TangentMap(1,2) * terr ) / determinant
    ;        ; tnew = Tsti(2) + ( - TangentMap(2,1) * serr + TangentMap(1,1) * terr ) / determinant
    ;        ; tnew = modulo( tnew, pi2 )

   case( 1 ) ; snew = Tsti(1) +                                                terr   / TangentMap(2,1)
    ;        ; tnew = zero    

   case default

    FATALMESS(pq00aa, .true., invalid Lpqsym )

   end select
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
!latex \item If the Newton correction takes the radial guess outside [\verb!low,upp!], which define the computational boundary, the new radial guess is random.
   
   if( snew.ge.low .and. snew.le.upp ) then ;                               Tsti(1:2) = (/ snew                    , tnew /) ! accept Newton correction;
   else                                     ; call random_number(harvest) ; Tsti(1:2) = (/ low + (upp-low)*harvest , tnew /) ! randomly restart;
   endif
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
  enddo ! end of do ; its=its+1 loop;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
 
  if( pqorbit(lvol,nc)%ok.eq.0 ) then ; ipqfail = 1 ; failed = 'failed' ; ! failed; no further action required;
  else                                ; ipqfail = 0 ; failed = '  ok  '
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  cput = GETTIME
  write(ounit,1002)cput-cpus, myid, lvol, pp, qq, failed, its, terr, serr, Tsti(2),  low, Tsti(1), upp, abs(determinant-one), pqorbit(lvol,nc)%residue
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( pqorbit(lvol,nc)%ok.eq.0 ) goto 9999   ! failed; no further action required;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  stz(1:3) = (/ Tsti(1) , Tsti(2) , zero /)
                                                         
#ifdef DEBUG
   FATALMESS(pq00aa, abs(stz(1)).gt.one, illegal value)
#endif

  call co00aa( lvol, stz(1:3), RpZ(1:3), dR(0:3), dZ(0:3), jacobian, guv(1:3,1:3) ) ! map to cylindrical coordinates; for Poincare diagnostic;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef OLDNAG

!latex \end{enumerate} \subsubsection{constructing unstable manifold} \begin{enumerate}

!latex \item Given the tangent map at the unstable periodic orbit, all that is required to identify the unstable manifold is to compute the eigenvectors.
!latex       The NAG routine \verb+F02EBF+ is employed for this.

  job = 'V' ; Nmat = 2 ; LDA = Nmat ; LDVR = Nmat ; LDVI = Nmat ; Lwk = 4 * Nmat
  
  TangentMap(1,1) = Tst(3) ! \partial_\s \bar \s
  TangentMap(1,2) = Tst(4) ! \partial_\t \bar \s
  TangentMap(2,1) = Tst(5) ! \partial_\s \bar \t
  TangentMap(2,2) = Tst(6) ! \partial_\t \bar \t
  
  ifail = 1 ; call F02EBF( job, Nmat, TangentMap(1:2,1:2), LDA, wr(1:2), wi(1:2), vr(1:2,1:2), LDVR, vi(1:2,1:2), LDVI, wk(1:Lwk), Lwk, ifail )
  
  cput = GETTIME
  select case( ifail ) !                                                12345678901
  case( 0 ) ; if( Wpq00aa ) write(ounit,1001)cput-cpus,myid,lvol,ifail,"success    ", pqorbit(lvol,nc)%residue, ( wr(iev), wi(iev), iev = 1, 2 )
  case( 1 ) ;               write(ounit,1001)cput-cpus,myid,lvol,ifail,"input error"
  case( 2 ) ;               write(ounit,1001)cput-cpus,myid,lvol,ifail,"QR failed  "
  case default
   FATALMESS(pq00aa,.true.,illegal ifail returned by F02EBF)
  end select

1001 format("pq00aa : ",f10.2," : myid=",i3," ; lvol=",i3," ; ifail="i2" ; constructing eigenvalues/vectors ; "a11" ; residue="f9.2" ; eval="2(f9.5" +"f9.5" i ,"))

#else
  
  wr(1:2) = zero ; wi(1:2) = zero ; vr(1:2,1:2) = zero ; vi(1:2,1:2) = zero

#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! save all information relevant to periodic orbits in pqorbit structure array;
  
  pqorbit(lvol,nc)%so          = Tsti(1)                                  ! location of periodic orbit on \z=0 plane; radial  ;
  pqorbit(lvol,nc)%to          = Tsti(2)                                  ! location of periodic orbit on \z=0 plane; poloidal;
! pqorbit(lvol,nc)%residue     = pqorbit(lvol,nc)%residue                 ! residue; ! this was assigned above; 11 Aug 13;
! pqorbit(lvol,nc)%error       = pqorbit(lvol,nc)%error                   ! error locating periodic orbit;
  pqorbit(lvol,nc)%its         = its                                      ! iterations required to locate periodic orbit;
  pqorbit(lvol,nc)%a(0:qq)     = (/ ( kk, kk = 0, qq ) /)      * pi2 / qq ! straight fieldline angle; this may be incorrect if Nfp > 1; 27 Nov 12;
  pqorbit(lvol,nc)%Ro          = RpZ(1)                                   ! cylindrical R; 
  pqorbit(lvol,nc)%Zo          = RpZ(3)                                   ! cylindrical Z; 
  pqorbit(lvol,nc)%wr(1:2)     = wr(1:2)                                  ! eigenvalues : real;
  pqorbit(lvol,nc)%wi(1:2)     = wi(1:2)                                  ! eigenvalues : imag;
  pqorbit(lvol,nc)%vr(1:2,1:2) = vr(1:2,1:2)                              ! eigenvector : real;
  pqorbit(lvol,nc)%vi(1:2,1:2) = vi(1:2,1:2)                              ! eigenvector : imag;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(pq00aa)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

1002 format("pq00aa : ":,f10.2," : myid=",i3," ; lvol=",i3," ; ",i3," /"i4" ; "a6" its=",i3," ; err="2es8.0" ; t="f8.4," ; [l,s,u]=["3f17.13" ] ;":" |det-1|="es7.0" ; R="f15.9" ;")
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine pq00aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine Tbfout( integrationparameter, Tst ) ! intermediate output along fieldline;

  use constants , only :
  use fileunits , only :
  use inputlist , only : 
  use allglobal , only : pi2nfp

  implicit none

  REAL, intent(inout) :: integrationparameter
  REAL, intent(in)    :: Tst(6)

  integrationparameter = integrationparameter + pi2nfp

  return

end subroutine Tbfout

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

REAL function Tbfend( integrationparameter, Tst ) ! termination of fieldline integration;

  use constants , only :
  use fileunits , only :
  use inputlist , only :
  use allglobal , only : pi2nfp

  use periodicitylength

  implicit none

  REAL, intent(in) :: integrationparameter, Tst(6)

  Tbfend = integrationparameter - ( zstart + pi2nfp )

  return

end function Tbfend

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
