!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Measures error in magnetic field.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \end{enumerate} \subsubsection{construction of current} \begin{enumerate}

!latex \item The magnetic field is given as the curl of the vector potential,
!latex \be {\bf B} = \frac{\partial_\t A_\z - \partial_\z A_\t}{\sqrt g}{\bf e}_\s
!latex             - \frac{                   \partial_\s A_\z}{\sqrt g}{\bf e}_\t
!latex             + \frac{\partial_\s A_\t                   }{\sqrt g}{\bf e}_\z
!latex             = B^\s {\bf e}_\s + B^\t {\bf e}_\t + B^\z {\bf e}_\z.
!latex \ee

!latex \item To take the curl of this, it is convenient to `lower' the metric elements
!latex \be B_\s&=&B^\s g_{\s\s}+B^\t g_{\s\t}+B^\z g_{\s\z},\\
!latex     B_\t&=&B^\s g_{\s\t}+B^\t g_{\s\t}+B^\z g_{\s\t},\\
!latex     B_\z&=&B^\s g_{\s\z}+B^\t g_{\s\z}+B^\z g_{\s\z}.
!latex \ee

!latex \item The curl is then 
!latex \be \nabla \times {\bf B} = \frac{\partial_\t B_\z - \partial_\z B_\t}{\sqrt g} \; {\bf e}_\s
!latex                           + \frac{\partial_\z B_\s - \partial_\s B_\z}{\sqrt g} \; {\bf e}_\t
!latex                           + \frac{\partial_\s B_\t - \partial_\t B_\s}{\sqrt g} \; {\bf e}_\z.
!latex \ee

!latex \item The derivatives of the metric elements, \verb+gij(1:6,0:3,1:Ntz)+, and the Jacobian, \verb+sg(0:3,1:Ntz)+, are calculated in \verb+co01aa+.

!latex \end{enumerate} \subsubsection{volume integrals} \begin{enumerate}
!latex \item Define 
!latex \be ||\left( \nabla \times {\bf B}-\mu {\bf B}\right)\cdot\nabla \s || & \!\!\!\!\!\! \equiv \!\!\!\!\!\! & 
!latex   \int \!\! ds \ooint \sqrt g \left| \left( \nabla \times {\bf B}-\mu {\bf B}\right)\cdot\nabla \s \right|
!latex = \int \!\! ds \ooint \left[(\partial_\t B_\z - \partial_\z B_\t) - \mu (\partial_\t A_\z - \partial_\z A_\t)\right], \nonumber \\
!latex     ||\left( \nabla \times {\bf B}-\mu {\bf B}\right)\cdot\nabla \t || & \!\!\!\!\!\! \equiv \!\!\!\!\!\! & 
!latex   \int \!\! ds \ooint \sqrt g \left|\left( \nabla \times {\bf B}-\mu {\bf B}\right)\cdot\nabla \t \right|
!latex =\int \!\! ds \ooint \left[(\partial_\z B_\s - \partial_\s B_\z) - \mu (\;\;\;\;\;\;\;\; - \partial_\s A_\z)\right], \nonumber \\
!latex     ||\left( \nabla \times {\bf B}-\mu {\bf B}\right)\cdot\nabla \z || & \!\!\!\!\!\! \equiv \!\!\!\!\!\! & 
!latex   \int \!\! ds \ooint \sqrt g \left|\left( \nabla \times {\bf B}-\mu {\bf B}\right)\cdot\nabla \z \right|
!latex =\int \!\! ds \ooint \left[(\partial_\s B_\t - \partial_\t B_\s) - \mu (\partial_\s A_\t   \;\;\;\;\;\;\;\;)\right], \nonumber
!latex \ee

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine jo00aa( lvol, Ntz, lrad, lquad, mn )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half, one, two

  use fileunits, only : ounit, lunit

  use inputlist, only : Wmacros, Wjo00aa, ext, Nvol, mu, curtor, curpol, epsr

  use cputiming, only : Tjo00aa

  use allglobal, only : myid, cpus, ivol, &
                        im, in, halfmm, &
                        Mvol, Lplasmaregion, Lvacuumregion, &
                        Ate, Aze, Ato, Azo, &
                        sg, guvij, Rij, Zij, &
                        trigwk, trigm, trign, isr, Nt, Nz, efmn, ofmn, cfmn, sfmn, &
                        NOTstellsym, &
                        Lcoordinatesingularity, & ! 11 Oct 12; 
                        IBerror, beltramierror  ! 11 Apr 16;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS
  
  INTEGER, intent(in) :: lvol, Ntz, lrad, lquad, mn ! these are really global, but are included in argument list to remove allocations;
  
  INTEGER             :: jquad, Lcurvature, ll, ii, jj, kk, uu, ideriv
  
  REAL                :: lss, gBu(1:Ntz,1:3,0:3), gJu(1:Ntz,1:3), jerror(1:3)
  
  REAL                :: chebyshev(0:lrad,0:2)
  
  REAL                :: Atemn(1:mn,0:2), Azemn(1:mn,0:2), Atomn(1:mn,0:2), Azomn(1:mn,0:2), sbar, sbarhim(0:2)
  
! required for Gaussian integration routine;
  INTEGER             :: itype, id01bcf
  REAL                :: aa, bb, cc, dd, jthweight, weight(1:lquad), abscis(1:lquad)
  
!required for adaptive integration routine;
  INTEGER             :: nfunctioncalls, nfunctioncallslimit=100000, ifail
  REAL                :: lowerlimit, upperlimit, lepsr, relerr, error(0:1), lcpu
  
  REAL                :: D01AHF, VacuumError
  external            :: D01AHF, VacuumError

  CHARACTER           :: svol*3
  
  BEGIN(jo00aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  FATALMESS(jo00aa, lquad.lt.1, invalid lquad supplied to jo00aa )
  FATALMESS(jo00aa, lvol.lt.1 .or. lvol.gt.Mvol, invalid interface label )
#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
 !write(svol,'(,i3,.3)')lvol
 !open(lunit+myid,file="."//trim(ext)//".curlerror."//svol,status="unknown")
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Lplasmaregion ) then ! this is the old integration method; new adaptive integration under construction in vacuum region; 10 Apr 13;
   
   itype = 0 ; aa = -one ; bb = +one ; cc = zero ; dd = zero ! prepare Gaussian quadrature;  6 Feb 13;
   
   id01bcf = 1
   call D01BCF( itype, aa, bb, cc, dd, lquad, weight(1:lquad), abscis(1:lquad), id01bcf ) ! sets gaussian weights & abscissae;
   
   cput= GETTIME
   select case( id01bcf ) !                                                         123456789012345
   case( 0 )    ;  if( Wjo00aa ) write(ounit,1000) cput-cpus, myid, lvol, id01bcf, "success        ", abscis(1:lquad)
   case( 1 )    ;                write(ounit,1000) cput-cpus, myid, lvol, id01bcf, "failed         ", abscis(1:lquad)
   case( 2 )    ;                write(ounit,1000) cput-cpus, myid, lvol, id01bcf, "input error    ", abscis(1:lquad)
   case( 3 )    ;                write(ounit,1000) cput-cpus, myid, lvol, id01bcf, "input error    ", abscis(1:lquad)
   case( 4 )    ;                write(ounit,1000) cput-cpus, myid, lvol, id01bcf, "weight overflow", abscis(1:lquad)
   case( 5 )    ;                write(ounit,1000) cput-cpus, myid, lvol, id01bcf, "weight zero    ", abscis(1:lquad)
   case( 6 )    ;                write(ounit,1000) cput-cpus, myid, lvol, id01bcf, "failed         ", abscis(1:lquad)
   case default ;                write(ounit,1000) cput-cpus, myid, lvol, id01bcf, "weird          ", abscis(1:lquad)
   end select
   
   ;               if( Wjo00aa ) write(ounit,1001)                                                    weight(1:lquad)
   
1000 format("jo00aa : ",f10.2," : myid=",i3," ; lvol=",i3," ; id01bcf=",i3," ; "a15" ;":" abscissae ="99f10.6)
1001 format("jo00aa : ", 10x ," :      "3x"        "3x"           "3x"   "15x" ;":" weights   ="99f10.6)
   
   FATALMESS(jo00aa,id01bcf.ne.0,failed to construct Gaussian integration abscisae and weights)
   
  endif ! end of if( Lplasmaregion) ; 20 Jun 14;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  jerror(1:3) = zero ; ideriv = 0 ! three components of the error in \curl B - mu B; initialize summation; 6 Feb 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Lplasmaregion ) then
   
   do jquad = 1, lquad ! loop over radial sub-sub-grid (numerical quadrature);
    
    lss = abscis(jquad) ; jthweight = weight(jquad) ; sbar = ( lss + one ) * half
    
    Lcurvature = 2
    
    WCALL(jo00aa,co01aa( lvol, lss, Lcurvature, Ntz, mn )) ! returns coordinates, metrics, . . .
    
    chebyshev(0,0:2) = (/ one, zero, zero /) ! Chebyshev initialization; 16 Jan 13;
    chebyshev(1,0:2) = (/ lss,  one, zero /) ! Chebyshev initialization; 16 Jan 13;
    
    do ll = 2, lrad
     chebyshev(ll,0:2) = (/ two * lss * chebyshev(ll-1,0)                                                                 - chebyshev(ll-2,0) , &
                            two       * chebyshev(ll-1,0) + two * lss * chebyshev(ll-1,1)                                 - chebyshev(ll-2,1) , &
                            two       * chebyshev(ll-1,1) + two       * chebyshev(ll-1,1) + two * lss * chebyshev(ll-1,2) - chebyshev(ll-2,2) /)
    enddo ! end of do ll; 20 Jun 14;
    
    Atemn(1:mn,0:2) = zero ! initialize summation over Chebyshev polynomials;  6 Feb 13;
    Azemn(1:mn,0:2) = zero
    Atomn(1:mn,0:2) = zero
    Azomn(1:mn,0:2) = zero
    
    do ll = 0, lrad
     
     do ii = 1, mn 
      
      if( Lcoordinatesingularity ) then
      !FATALMESS(jo00aa, .true., need to revise construction of current) ! 04 Dec 14; ! 16 Jan 15;
       sbarhim(0) = sbar**halfmm(ii) ; sbarhim(1) = half * halfmm(ii) * sbarhim(0) / sbar ; sbarhim(2) = half * ( halfmm(ii)-one ) * sbarhim(1) / sbar
      else                              
       sbarhim(0) = one              ; sbarhim(1) = zero                                  ; sbarhim(2) = zero
      endif
      
      ;Atemn(ii,0) = Atemn(ii,0) + Ate(lvol,ideriv,ii)%s(ll) * ( chebyshev(ll,0) * sbarhim(0) )
      ;Atemn(ii,1) = Atemn(ii,1) + Ate(lvol,ideriv,ii)%s(ll) * ( chebyshev(ll,1) * sbarhim(0) + chebyshev(ll,0) * sbarhim(1) )
      ;Atemn(ii,2) = Atemn(ii,2) + Ate(lvol,ideriv,ii)%s(ll) * ( chebyshev(ll,2) * sbarhim(0) + two * chebyshev(ll,1) * sbarhim(1) + chebyshev(ll,0) * sbarhim(2) )
      
      ;Azemn(ii,0) = Azemn(ii,0) + Aze(lvol,ideriv,ii)%s(ll) * ( chebyshev(ll,0) * sbarhim(0) )
      ;Azemn(ii,1) = Azemn(ii,1) + Aze(lvol,ideriv,ii)%s(ll) * ( chebyshev(ll,1) * sbarhim(0) + chebyshev(ll,0) * sbarhim(1) )
      ;Azemn(ii,2) = Azemn(ii,2) + Aze(lvol,ideriv,ii)%s(ll) * ( chebyshev(ll,2) * sbarhim(0) + two * chebyshev(ll,1) * sbarhim(1) + chebyshev(ll,0) * sbarhim(2) )
      if( NOTstellsym ) then
       Atomn(ii,0) = Atomn(ii,0) + Ato(lvol,ideriv,ii)%s(ll) * ( chebyshev(ll,0) * sbarhim(0) )
       Atomn(ii,1) = Atomn(ii,1) + Ato(lvol,ideriv,ii)%s(ll) * ( chebyshev(ll,1) * sbarhim(0) + chebyshev(ll,0) * sbarhim(1) )
       Atomn(ii,2) = Atomn(ii,2) + Ato(lvol,ideriv,ii)%s(ll) * ( chebyshev(ll,2) * sbarhim(0) + two * chebyshev(ll,1) * sbarhim(1) + chebyshev(ll,0) * sbarhim(2) )
       
       Azomn(ii,0) = Azomn(ii,0) + Azo(lvol,ideriv,ii)%s(ll) * ( chebyshev(ll,0) * sbarhim(0) )
       Azomn(ii,1) = Azomn(ii,1) + Azo(lvol,ideriv,ii)%s(ll) * ( chebyshev(ll,1) * sbarhim(0) + chebyshev(ll,0) * sbarhim(1) )
       Azomn(ii,2) = Azomn(ii,2) + Azo(lvol,ideriv,ii)%s(ll) * ( chebyshev(ll,2) * sbarhim(0) + two * chebyshev(ll,1) * sbarhim(1) + chebyshev(ll,0) * sbarhim(2) )
      endif ! end of if( NOTstellsym) ; 20 Jun 14;
      
     !;Atemn(ii,0:2) = Atemn(ii,0:2) + Ate(lvol,ideriv,ii)%s(ll) * chebyshev(ll,0:2)
     !;Azemn(ii,0:2) = Azemn(ii,0:2) + Aze(lvol,ideriv,ii)%s(ll) * chebyshev(ll,0:2)
     !;if( NOTstellsym ) then
     !;Atomn(ii,0:2) = Atomn(ii,0:2) + Ato(lvol,ideriv,ii)%s(ll) * chebyshev(ll,0:2)
     !;Azomn(ii,0:2) = Azomn(ii,0:2) + Azo(lvol,ideriv,ii)%s(ll) * chebyshev(ll,0:2)
     !;endif 
      
     enddo ! end of do ii;  6 Feb 13;
    
    enddo ! end of do ll;  6 Feb 13;
    
    ofmn(1:mn) = -          im(1:mn)*Azemn(1:mn,0) -          in(1:mn)*Atemn(1:mn,0)
    efmn(1:mn) = +          im(1:mn)*Azomn(1:mn,0) +          in(1:mn)*Atomn(1:mn,0)
    sfmn(1:mn) = -          im(1:mn)*Azemn(1:mn,1) -          in(1:mn)*Atemn(1:mn,1)
    cfmn(1:mn) = +          im(1:mn)*Azomn(1:mn,1) +          in(1:mn)*Atomn(1:mn,1)
    
    call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), &
                 Nt, Nz, gBu(1:Ntz,1,0), gBu(1:Ntz,1,1), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz) )

    efmn(1:mn) = - im(1:mn)*im(1:mn)*Azemn(1:mn,0) - im(1:mn)*in(1:mn)*Atemn(1:mn,0)
    ofmn(1:mn) = - im(1:mn)*im(1:mn)*Azomn(1:mn,0) - im(1:mn)*in(1:mn)*Atomn(1:mn,0)
    cfmn(1:mn) = + in(1:mn)*im(1:mn)*Azemn(1:mn,0) + in(1:mn)*in(1:mn)*Atemn(1:mn,0)
    sfmn(1:mn) = + in(1:mn)*im(1:mn)*Azomn(1:mn,0) + in(1:mn)*in(1:mn)*Atomn(1:mn,0)
    
    call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), &
                 Nt, Nz, gBu(1:Ntz,1,2), gBu(1:Ntz,1,3), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz) )
    
    efmn(1:mn) =                                   -                   Azemn(1:mn,1)
    ofmn(1:mn) =                                   -                   Azomn(1:mn,1)
    cfmn(1:mn) =                                   -                   Azemn(1:mn,2)
    sfmn(1:mn) =                                   -                   Azomn(1:mn,2)
    
    call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), &
                 Nt, Nz, gBu(1:Ntz,2,0), gBu(1:Ntz,2,1), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz) )

    ofmn(1:mn) =                                   +          im(1:mn)*Azemn(1:mn,1)
    efmn(1:mn) =                                   -          im(1:mn)*Azomn(1:mn,1)
    sfmn(1:mn) =                                   -          in(1:mn)*Azemn(1:mn,1)
    cfmn(1:mn) =                                   +          in(1:mn)*Azomn(1:mn,1)
    
    call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), &
                 Nt, Nz, gBu(1:Ntz,2,2), gBu(1:Ntz,2,3), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz) )
   
    efmn(1:mn) = +                   Atemn(1:mn,1)
    ofmn(1:mn) = +                   Atomn(1:mn,1)
    cfmn(1:mn) = +                   Atemn(1:mn,2)
    sfmn(1:mn) = +                   Atomn(1:mn,2)
    
    call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), &
                 Nt, Nz, gBu(1:Ntz,3,0), gBu(1:Ntz,3,1), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz) )

    ofmn(1:mn) = -          im(1:mn)*Atemn(1:mn,1)
    efmn(1:mn) = +          im(1:mn)*Atomn(1:mn,1)
    sfmn(1:mn) = +          in(1:mn)*Atemn(1:mn,1)
    cfmn(1:mn) = -          in(1:mn)*Atomn(1:mn,1)
    
    call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), &
                 Nt, Nz, gBu(1:Ntz,3,2), gBu(1:Ntz,3,3), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz) )
     
    do ii = 1, 3
     
     select case( ii )
     case( 1 ) ; jj = 2 ; kk = 3
     case( 2 ) ; jj = 3 ; kk = 1
     case( 3 ) ; jj = 1 ; kk = 2
     end select
     
     gJu(1:Ntz,ii) = zero
     
     do uu = 1, 3 ! summation over uu;  6 Feb 13;
      
      gJu(1:Ntz,ii) = gJu(1:Ntz,ii) &
                    + ( gBu(1:Ntz,uu,jj) * guvij(1:Ntz,uu,kk,0) + gBu(1:Ntz,uu,0) * guvij(1:Ntz,uu,kk,jj) - gBu(1:Ntz,uu,0) * guvij(1:Ntz,uu,kk,0) * sg(1:Ntz,jj) / sg(1:Ntz,0) ) &
                    - ( gBu(1:Ntz,uu,kk) * guvij(1:Ntz,uu,jj,0) + gBu(1:Ntz,uu,0) * guvij(1:Ntz,uu,jj,kk) - gBu(1:Ntz,uu,0) * guvij(1:Ntz,uu,jj,0) * sg(1:Ntz,kk) / sg(1:Ntz,0) )
     enddo ! end of do uu;  6 Feb 13;
     
     gJu(1:Ntz,ii) = gJu(1:Ntz,ii) / sg(1:Ntz,0)
     
    enddo ! end of do ii;  5 Feb 13;
    
    do ii = 1, 3 ; jerror(ii) = jerror(ii) + jthweight * sum( abs( gJu(1:Ntz,ii) - mu(lvol) * gBu(1:Ntz,ii,0) ) )
    enddo
    
   enddo ! end of do jquad;  5 Feb 13;
   
   jerror(1:3) = jerror(1:3) / Ntz
   
   cput = GETTIME ; write(ounit     ,1002) cput-cpus, myid, lvol, lrad, jerror(1:3), cput-cpui ! write error to screen;
   
   beltramierror(lvol,1:3) = jerror(1:3)   ! store error to write in h5 file; 11 Apr 16;
   
  endif ! end of if( Lplasmaregion ) ; 10 Apr 13;
  
1002 format("jo00aa : ",f10.2," : myid=",i3," ; lvol=",i3," ; lrad=",i3," ; error "3es13.5 " ; time="f8.2"s ;")
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Lvacuumregion ) then
   
   do IBerror = 0, 1 ! calculate numerator and denominator; 29 Apr 13;
    
    ivol = lvol ! volume identification label passed through global to VacuumError; 10 Apr 13;
    
    lcpu = GETTIME ! record how long the adaptive integration takes; 10 Apr 13;
    
    ifail = 1 ; lowerlimit = -one ; upperlimit = one ; lepsr = epsr ; relerr = -one ! epsr is given on input; 29 Apr 13;
    error(IBerror) = D01AHF( lowerlimit, upperlimit, lepsr, nfunctioncalls, relerr, VacuumError, nfunctioncallslimit, ifail )
    
    cput = GETTIME
    
    select case( ifail )                                                                                                                        !123456789012345678
    case( 0 ) ;               write(ounit,1003) cput-cpus, myid, lvol, lrad, error(IBerror), IBerror, ifail, nfunctioncalls, relerr, cput-lcpu, "success ;         "
    case( 1 ) ;               write(ounit,1003) cput-cpus, myid, lvol, lrad, error(IBerror), IBerror, ifail, nfunctioncalls, relerr, cput-lcpu, "NLIMIT too small ;"
    case( 2 ) ;               write(ounit,1003) cput-cpus, myid, lvol, lrad, error(IBerror), IBerror, ifail, nfunctioncalls, relerr, cput-lcpu, "error ;           "
    case( 3 ) ;               write(ounit,1003) cput-cpus, myid, lvol, lrad, error(IBerror), IBerror, ifail, nfunctioncalls, relerr, cput-lcpu, "epsr < 0.0 ;      "
    case default
     FATALMESS(jo00aa, .true., invalid ifail )
    end select
    
1003 format("jo00aa : ",f10.2," : myid=",i3," ; ivol=",i3," ; lrad=",i3," ; error=" es13.5" ; IB="i2" ; ifail="i2" ; #calls="i6" ; relerr="es8.1" ; time="f7.1"s ; ",a18)
    
   enddo ! end of do( IBerror ) ; 29 Apr 13;
   
   write(ounit,1004) cput-cpus, myid, lvol, lrad, error(0)*(/1,1,1/)/error(1)
1004 format("jo00aa : ",f10.2," : myid=",i3," ; lvol=",i3," ; lrad=",i3," ; error="3es13.5" ; ")
   
  endif ! end of if( Lvacuumregion ) ; 10 Apr 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
 !open(lunit+myid)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(jo00aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine jo00aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

REAL function BeltramiError( lss )

  use constants, only : zero
  use allglobal, only : ncpu, myid, cpus, ivol, pi2nfp
  
  LOCALS
  
  REAL                :: lss

  FATALMESS(jo00aa, .true., under construction )
  
  BeltramiError = zero

  return
  
end function BeltramiError

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \end{enumerate} \subsubsection{construction of ``error function" for scalar potential} \begin{enumerate}

!latex \item We may quantify the error in the scalar potential by computing the following integral:
!latex       \be     \int_{\cal V} \left| \nabla \cdot \nabla \Phi \right| dv 
!latex            =  \int \!\!\! \int \!\!\! \int \!\! d\s \, d\t \, d\z \left[\partial_\s ( \sqrt g B^\s ) + \partial_\t ( \sqrt g B^\t ) + \partial_\z ( \sqrt g B^\z ) \right]
!latex       ,\ee
!latex       where $\sqrt g B^\mu = I \sqrt g \, g^{\t\mu} + G \sqrt g \, g^{\z\mu} + \sum_{\nu=\s,\t,\z} \phi_\nu \sqrt g \, g^{\nu\mu}$ for $\mu=\s,\t,\z$.

!latex \item This quantity is normalized by
!latex       \be     \int_{\cal V} B^2 dv.
!latex       \ee

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

REAL function VacuumError( lss )
  
  use constants, only : zero, half, one, two, pi2
  
  use fileunits, only : ounit
  
  use inputlist, only : Wmacros, Wjo00aa, ext, Nvol, mu, curtor, curpol, Lrad, epsr
  
  use cputiming, only : Tjo00aa

  use allglobal, only : ncpu, myid, cpus, ivol, pi2nfp, &
                        lchebyshev, &
                        mn, im, in, &
                        Mvol, Lplasmaregion, Lvacuumregion, &
                        Ate, Aze, Ato, Azo, &
                        Ntz, sg, guvij, Rij, Zij, gvuij, cosi, sini, &
                        ijreal, mne, ime, ine, &
                        trigwk, trigm, trign, isr, Nt, Nz, efmn, ofmn, cfmn, sfmn, &
                        NOTstellsym, &
                        Lcoordinatesingularity, & ! 11 Oct 12; 
                        IBerror
  
  LOCALS
  
  REAL                :: lss, lcpu

  INTEGER             :: lvol, ll, Lcurvature, lNtz, lmn, ii, jj, kk, dd, mi, ni, ifail, ideriv

  REAL, allocatable   :: dPHIij(:,:,:), errij(:)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  BEGIN(jo00aa)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! VacuumError = zero

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  lvol = ivol ; lmn = mn ; lNtz = Ntz ; ideriv = 0 ! shorthand; 10 Apr 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! construct Chebyshev polynomials; 10 Apr 13;
  
  lchebyshev(0,0:2) = (/ one, zero, zero /) ! Chebyshev initialization; 16 Jan 13;
  lchebyshev(1,0:2) = (/ lss,  one, zero /) ! Chebyshev initialization; 16 Jan 13;
  do ll = 2, Lrad(lvol)
   lchebyshev(ll,0:2) = (/ two * lss * lchebyshev(ll-1,0)                                  - lchebyshev(ll-2,0) , & ! Chebyshev recurrence;            ; 16 Jan 13;
                           two       * lchebyshev(ll-1,0) + two * lss * lchebyshev(ll-1,1) - lchebyshev(ll-2,1) , &
                           two * two * lchebyshev(ll-1,1) + two * lss * lchebyshev(ll-1,2) - lchebyshev(ll-2,2) /)  ! Chebyshev recurrence; derivatives; 16 Jan 13;
  enddo
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! determine "lower" metric elements and their derivatives; 10 Apr 13;
  
  Lcurvature = 2
  CALL(jo00aa, co01aa, ( lvol, lss, Lcurvature, lNtz, lmn )) ! provides metric elements on grid; 10 Apr 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! construct "upper" metric elements; note the extra Jacobian factor; 10 Apr 13;
  
  do ii = 1, 3
   
   select case( ii )
   case( 1 ) ; jj = 2 ; kk = 3
   case( 2 ) ; jj = 3 ; kk = 1
   case( 3 ) ; jj = 1 ; kk = 2
   end select
   
   gvuij(1:Ntz,ii, 1, 0) = ( guvij(1:Ntz,jj, 2, 0)*guvij(1:Ntz,kk, 3, 0) - guvij(1:Ntz,jj, 3, 0)*guvij(1:Ntz,kk, 2, 0) ) / sg(1:Ntz, 0)
   gvuij(1:Ntz,ii, 2, 0) = ( guvij(1:Ntz,jj, 3, 0)*guvij(1:Ntz,kk, 1, 0) - guvij(1:Ntz,jj, 1, 0)*guvij(1:Ntz,kk, 3, 0) ) / sg(1:Ntz, 0)
   gvuij(1:Ntz,ii, 3, 0) = ( guvij(1:Ntz,jj, 1, 0)*guvij(1:Ntz,kk, 2, 0) - guvij(1:Ntz,jj, 2, 0)*guvij(1:Ntz,kk, 1, 0) ) / sg(1:Ntz, 0)
   
  enddo
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! differentiate upper metric elements; note the extra Jacobian factor; 10 Apr 13;
  
  do dd = 1, 3 ! labels derivative; 10 Apr 13;
   
   do ii = 1, 3
    
    select case( ii )
    case( 1 ) ; jj = 2 ; kk = 3
    case( 2 ) ; jj = 3 ; kk = 1
    case( 3 ) ; jj = 1 ; kk = 2
    end select
    
    gvuij(1:Ntz,ii, 1,dd) = ( guvij(1:Ntz,jj, 2,dd)*guvij(1:Ntz,kk, 3, 0) - guvij(1:Ntz,jj,3,dd)*guvij(1:Ntz,kk,2, 0) &
                            + guvij(1:Ntz,jj, 2, 0)*guvij(1:Ntz,kk, 3,dd) - guvij(1:Ntz,jj,3, 0)*guvij(1:Ntz,kk,2,dd) - gvuij(1:Ntz,ii, 1, 0)*sg(1:Ntz,dd) ) / sg(1:Ntz, 0)
    
    gvuij(1:Ntz,ii, 2,dd) = ( guvij(1:Ntz,jj, 3,dd)*guvij(1:Ntz,kk, 1, 0) - guvij(1:Ntz,jj,1,dd)*guvij(1:Ntz,kk,3, 0) &
                            + guvij(1:Ntz,jj, 3, 0)*guvij(1:Ntz,kk, 1,dd) - guvij(1:Ntz,jj,1, 0)*guvij(1:Ntz,kk,3,dd) - gvuij(1:Ntz,ii, 2, 0)*sg(1:Ntz,dd) ) / sg(1:Ntz, 0)
    
    gvuij(1:Ntz,ii, 3,dd) = ( guvij(1:Ntz,jj, 1,dd)*guvij(1:Ntz,kk, 2, 0) - guvij(1:Ntz,jj,2,dd)*guvij(1:Ntz,kk,1, 0) &
                            + guvij(1:Ntz,jj, 1, 0)*guvij(1:Ntz,kk, 2,dd) - guvij(1:Ntz,jj,2, 0)*guvij(1:Ntz,kk,1,dd) - gvuij(1:Ntz,ii, 3, 0)*sg(1:Ntz,dd) ) / sg(1:Ntz, 0)

   enddo

  enddo

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!  construct contra-variant components of (the non-secular components of) the field on regular angle grid; 10 Apr 13;
  
  RALLOCATE(dPHIij,(1:Ntz,1:3,0:3)) ! recall that this initializes to zero; 10 Apr 13;
  
  do ii = 2, mn ; mi = im(ii) ; ni = in(ii) ! loop over Fourier resolution; 10 Apr 13;
   
   do ll = 0, Lrad(lvol) ! loop over Chebyshev summation; 20 Feb 13;
    
    dPHIij(1:Ntz,1,0) = dPHIij(1:Ntz,1,0) + Ate(lvol,ideriv,ii)%s(ll) * lchebyshev(ll,1) * sini(1:Ntz,ii)
    
    dPHIij(1:Ntz,1,1) = dPHIij(1:Ntz,1,1) + Ate(lvol,ideriv,ii)%s(ll) * lchebyshev(ll,2) * sini(1:Ntz,ii)
    dPHIij(1:Ntz,1,2) = dPHIij(1:Ntz,1,2) + Ate(lvol,ideriv,ii)%s(ll) * lchebyshev(ll,1) * cosi(1:Ntz,ii) * ( + mi )
    dPHIij(1:Ntz,1,3) = dPHIij(1:Ntz,1,3) + Ate(lvol,ideriv,ii)%s(ll) * lchebyshev(ll,1) * cosi(1:Ntz,ii) * ( - ni )

    dPHIij(1:Ntz,2,0) = dPHIij(1:Ntz,2,0) + Ate(lvol,ideriv,ii)%s(ll) * lchebyshev(ll,0) * cosi(1:Ntz,ii) * ( + mi )

    dPHIij(1:Ntz,2,1) = dPHIij(1:Ntz,2,1) + Ate(lvol,ideriv,ii)%s(ll) * lchebyshev(ll,1) * cosi(1:Ntz,ii) * ( + mi )
    dPHIij(1:Ntz,2,2) = dPHIij(1:Ntz,2,2) + Ate(lvol,ideriv,ii)%s(ll) * lchebyshev(ll,0) * sini(1:Ntz,ii) * ( + mi ) * ( - mi )
    dPHIij(1:Ntz,2,3) = dPHIij(1:Ntz,2,3) + Ate(lvol,ideriv,ii)%s(ll) * lchebyshev(ll,0) * sini(1:Ntz,ii) * ( + mi ) * ( + ni )
    
    dPHIij(1:Ntz,3,0) = dPHIij(1:Ntz,3,0) + Ate(lvol,ideriv,ii)%s(ll) * lchebyshev(ll,0) * cosi(1:Ntz,ii) * ( - ni )

    dPHIij(1:Ntz,3,1) = dPHIij(1:Ntz,3,1) + Ate(lvol,ideriv,ii)%s(ll) * lchebyshev(ll,1) * cosi(1:Ntz,ii) * ( - ni )
    dPHIij(1:Ntz,3,2) = dPHIij(1:Ntz,3,2) + Ate(lvol,ideriv,ii)%s(ll) * lchebyshev(ll,0) * sini(1:Ntz,ii) * ( - ni ) * ( - mi )
    dPHIij(1:Ntz,3,3) = dPHIij(1:Ntz,3,3) + Ate(lvol,ideriv,ii)%s(ll) * lchebyshev(ll,0) * sini(1:Ntz,ii) * ( - ni ) * ( + ni )
    
   enddo
   
  enddo
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! compute \nabla \cdot \nabla \Phi, which is a measure of the numerical accuracy; 10 Apr 13;
  
  RALLOCATE(errij,(1:Ntz))
  
  if( IBerror.eq.0 ) then ! compute numerator; 29 Apr 13;
   
   errij(1:Ntz) = curtor * gvuij(1:Ntz,2,1,1) + curpol * gvuij(1:Ntz,3,1,1) & 
                + dPHIij(1:Ntz,1,1) * gvuij(1:Ntz,1,1,0) + dPHIij(1:Ntz,2,1) * gvuij(1:Ntz,2,1,0) + dPHIij(1:Ntz,3,1) * gvuij(1:Ntz,3,1,0) & 
                + dPHIij(1:Ntz,1,0) * gvuij(1:Ntz,1,1,1) + dPHIij(1:Ntz,2,0) * gvuij(1:Ntz,2,1,1) + dPHIij(1:Ntz,3,0) * gvuij(1:Ntz,3,1,1) & 
                + curtor * gvuij(1:Ntz,2,2,2) + curpol * gvuij(1:Ntz,3,2,2) & 
                + dPHIij(1:Ntz,1,2) * gvuij(1:Ntz,1,2,0) + dPHIij(1:Ntz,2,2) * gvuij(1:Ntz,2,2,0) + dPHIij(1:Ntz,3,2) * gvuij(1:Ntz,3,2,0) & 
                + dPHIij(1:Ntz,1,0) * gvuij(1:Ntz,1,2,2) + dPHIij(1:Ntz,2,0) * gvuij(1:Ntz,2,2,2) + dPHIij(1:Ntz,3,0) * gvuij(1:Ntz,3,2,2) & 
                + curtor * gvuij(1:Ntz,2,3,3) + curpol * gvuij(1:Ntz,3,3,3) & 
                + dPHIij(1:Ntz,1,3) * gvuij(1:Ntz,1,3,0) + dPHIij(1:Ntz,2,3) * gvuij(1:Ntz,2,3,0) + dPHIij(1:Ntz,3,3) * gvuij(1:Ntz,3,3,0) & 
                + dPHIij(1:Ntz,1,0) * gvuij(1:Ntz,1,3,3) + dPHIij(1:Ntz,2,0) * gvuij(1:Ntz,2,3,3) + dPHIij(1:Ntz,3,0) * gvuij(1:Ntz,3,3,3)
   
  endif ! end of if( IBerror.eq.0 ) ; 29 Apr 13;
  
  if( IBerror.eq.1 ) then ! compute denominator; 29 Apr 13;
   
   errij(1:Ntz) = curtor            * curtor            * gvuij(1:Ntz,2,2,0) &
                + curtor            * curpol            * gvuij(1:Ntz,3,2,0) &
                + curtor            * dPHIij(1:Ntz,1,0) * gvuij(1:Ntz,1,2,0) &
                + curtor            * dPHIij(1:Ntz,2,0) * gvuij(1:Ntz,2,2,0) &
                + curtor            * dPHIij(1:Ntz,3,0) * gvuij(1:Ntz,3,2,0) &
                + curpol            * curtor            * gvuij(1:Ntz,2,3,0) &
                + curpol            * curpol            * gvuij(1:Ntz,3,3,0) &
                + curpol            * dPHIij(1:Ntz,1,0) * gvuij(1:Ntz,1,3,0) &
                + curpol            * dPHIij(1:Ntz,2,0) * gvuij(1:Ntz,2,3,0) &
                + curpol            * dPHIij(1:Ntz,3,0) * gvuij(1:Ntz,3,3,0) &
                + dPHIij(1:Ntz,1,0) * curtor            * gvuij(1:Ntz,2,1,0) &
                + dPHIij(1:Ntz,1,0) * curpol            * gvuij(1:Ntz,3,1,0) &
                + dPHIij(1:Ntz,1,0) * dPHIij(1:Ntz,1,0) * gvuij(1:Ntz,1,1,0) &
                + dPHIij(1:Ntz,1,0) * dPHIij(1:Ntz,2,0) * gvuij(1:Ntz,2,1,0) &
                + dPHIij(1:Ntz,1,0) * dPHIij(1:Ntz,3,0) * gvuij(1:Ntz,3,1,0) &
                + dPHIij(1:Ntz,2,0) * curtor            * gvuij(1:Ntz,2,2,0) &
                + dPHIij(1:Ntz,2,0) * curpol            * gvuij(1:Ntz,3,2,0) &
                + dPHIij(1:Ntz,2,0) * dPHIij(1:Ntz,1,0) * gvuij(1:Ntz,1,2,0) &
                + dPHIij(1:Ntz,2,0) * dPHIij(1:Ntz,2,0) * gvuij(1:Ntz,2,2,0) &
                + dPHIij(1:Ntz,2,0) * dPHIij(1:Ntz,3,0) * gvuij(1:Ntz,3,2,0) &
                + dPHIij(1:Ntz,3,0) * curtor            * gvuij(1:Ntz,2,3,0) &
                + dPHIij(1:Ntz,3,0) * curpol            * gvuij(1:Ntz,3,3,0) &
                + dPHIij(1:Ntz,3,0) * dPHIij(1:Ntz,1,0) * gvuij(1:Ntz,1,3,0) &
                + dPHIij(1:Ntz,3,0) * dPHIij(1:Ntz,2,0) * gvuij(1:Ntz,2,3,0) &
                + dPHIij(1:Ntz,3,0) * dPHIij(1:Ntz,3,0) * gvuij(1:Ntz,3,3,0)
   
  endif ! end of if( IBerror.eq.1 ) ; 29 Apr 13;


  errij(1:Ntz) = abs( errij(1:Ntz) )
  
  VacuumError = sum( errij(1:Ntz) ) * pi2*pi2nfp / Ntz ! normalization of error function; 10 Apr 13;


  DEALLOCATE(errij)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  DEALLOCATE(dPHIij)

  return
  
end function VacuumError

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
