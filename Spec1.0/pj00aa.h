!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Checks pressure balance across arbitrary surface. 

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine pj00aa( mn , Rmn , Zmn , Itori , Gpoli , fmni , presi , Itoro , Gpolo , fmno , preso , sumerr )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  use constants, only : zero, half, one, two, pi2
  use fileunits, only : ounit
  use numerical, only : small
  use inputlist, only : Wpj00aa, Igeometry
  use cputiming, only : Tpj00aa
  use allglobal, only : myid, cpus, im, in, pi2nfp, Nt, Nz, Ntz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS

  INTEGER , intent(in)  :: mn
  REAL    , intent(in)  :: Rmn(1:mn), Zmn(1:mn), Itori, Gpoli, fmni(1:mn), Itoro, Gpolo, fmno(1:mn), presi, preso
  REAL    , intent(out) :: sumerr(2)

  INTEGER               :: jj, kk, imn
  REAL                  :: stz(1:3), dR(0:3), dZ(0:3), fi(0:3), fo(0:3), gl(3, 3), arg, carg, sarg, Bsqdi, Bsqdo, denom

  BEGIN(pj00aa)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  sumerr(1:2) = (/ zero , zero /)

  do kk = 0 , Nz-1
   do jj = 0 , Nt-1
    
    stz(2:3) = (/ jj , kk /) * (/ pi2 , pi2nfp /) / (/ Nt , Nz /)

    dR(0:3) = zero
    dZ(0:3) = zero

    fi(0:3) = (/ zero , zero , Itori , Gpoli /)
    fo(0:3) = (/ zero , zero , Itoro , Gpolo /)

    do imn = 1 , mn
     
     arg = im(imn) * stz(2) - in(imn) * stz(3) ; carg = cos(arg) ; sarg = sin(arg)
     
!latex \item The geometry of the interface, ${\bf x}(\t,\z)$, is given, where $\t$ and $\z$ are arbitrary angle parameters.

     dR(0) = dR(0)+ Rmn(imn) * carg
!    dR(1) ! not required; not even defined in this subroutine;
     dR(2) = dR(2)+ Rmn(imn) * sarg * ( -im(imn) )
     dR(3) = dR(3)+ Rmn(imn) * sarg * ( +in(imn) )

     dZ(0) = dZ(0)+ Zmn(imn) * sarg
!    dZ(1) ! not required; not even defined in this subroutine;
     dZ(2) = dZ(2)+ Zmn(imn) * carg * ( +im(imn) )
     dZ(3) = dZ(3)+ Zmn(imn) * carg * ( -in(imn) )

!latex \item The ``surface potential'', $f(\t,\z)=I\t+G\z+\tilde f(\t,\z)$, is given, where $I$ and $G$ are constants and $\tilde f$ is periodic.

!latex \item The covariant components of the field are determined by
!latex \be B_\t & = & \partial_\t f,\\
!latex     B_\z & = & \partial_\z f,\\
!latex     B_\s & = & ( - g^{\s\t} B_\t - g^{\s\z} B_\z ) / g^{\s\s},
!latex \ee
!latex where we have assumed that $B^\s=0$.

!    fi(0) = fi(0) + fmni(imn) * sarg                                
!    fi(1) ! not required or defined;
     fi(2) = fi(2) + fmni(imn) * carg * ( +im(imn) )
     fi(3) = fi(3) + fmni(imn) * carg * ( -in(imn) )

!    fo(0) = fo(0) + fmno(imn) * sarg                                
!    fo(1) ! not required or defined;
     fo(2) = fo(2) + fmno(imn) * carg * ( +im(imn) )
     fo(3) = fo(3) + fmno(imn) * carg * ( -in(imn) )

    enddo ! end of do imn=1,mn ;


    select case( Igeometry ) !latex \item The choice of geometry determines the metrics, and we may have: \begin{itemize}
     
    case(   1) !latex \item \verb+Igeometry.eq.1+ : Cartesian : ${\bf x} = \t \; \hat {\bf i} + \z \; \hat {\bf j} + R(\s,\t,\z) \; \hat {\bf k}$;
     
     gl(2,2) = dR(2) * dR(2) + one
     gl(2,3) = dR(2) * dR(3)
     gl(3,2) = dR(3) * dR(2)
     gl(3,3) = dR(3) * dR(3) + one

!   case(   2) ! cylindrical;

!   case(   3) ! toroidal;

!   case(   4) ! Cartesian;

!   case(   5) ! cylindrical;

!   case(-6,6) !latex \item \verb+Igeometry.eq.6+ : toroidal : ${\bf x} = R(\s,\t,\z) \; \hat {\bf R} + Z(\s,\t,\z) \; \hat {\bf k}$;

!    gl(2,2) = dR(2) * dR(2) + dZ(2) * dZ(2)                 ! g_t,t;
!    gl(2,3) = dR(2) * dR(3) + dZ(2) * dZ(3)                 ! g_t,z;
!    gl(3,2) = dR(3) * dR(2) + dZ(3) * dZ(2)                 ! g_z,t;
!    gl(3,3) = dR(3) * dR(3) + dZ(3) * dZ(3) + dR(0) * dR(0) ! g_z,z;

    case default

     FATALMESS(pj00aa,.true.,Igeometry is not supported)

    end select !latex \end{itemize}
    
!latex \item The expresion for $B^2$ can be written
!latex \be B^2 = \frac{ g_{\z\z} B_\t B_\t - 2 g_{\t\z} B_\t B_\z + g_{\t\t} B_\z B_\z}{ g_{\t\t} g_{\z\z} - g_{\t\z} g_{\t\z} }
!latex \ee

    denom = gl(2,2) * gl(3,3) - gl(2,3) * gl(3,2)

    FATALMESS(pj00aa,abs(denom).lt.small,zero denominator)

    Bsqdi = ( gl(3,3) * fi(2) * fi(2) - two * gl(2,3) * fi(2) * fi(3) + gl(2,2) * fi(3) * fi(3) ) / denom
    Bsqdo = ( gl(3,3) * fo(2) * fo(2) - two * gl(2,3) * fo(2) * fo(3) + gl(2,2) * fo(3) * fo(3) ) / denom

    sumerr(1) = sumerr(1) +   two * presi + Bsqdi + two * preso + Bsqdo 
    sumerr(2) = sumerr(2) + ( two * presi + Bsqdi - two * preso - Bsqdo )**2
    
   enddo ! end of poloidal grid loop
  enddo ! end of toroidal grid loop
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  sumerr(1) = half * sumerr(1) / Ntz
  sumerr(2) = sqrt(  sumerr(2) / Ntz )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(pj00aa)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine pj00aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
