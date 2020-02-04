!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (diagnostic) ! Returns the magnetic field field line equations, $d{\bf x}/d\phi = {\bf B}/B^\phi$.

!latex \briefly{Returns $\dot \s \equiv B^\s / B^\z$ and $\dot \t \equiv B^\t / B^\z$.}

!latex \calledby{\link{pp00ab}}
!      \calls{}

!latex \tableofcontents

!latex \subsection{equations of field line flow} \begin{enumerate}

!latex \item The equations for the fieldlines are normalized to the toroidal field, i.e. 
!latex       \be \dot \s \equiv \frac{B^\s}{B^\z}, \qquad
!latex           \dot \t \equiv \frac{B^\t}{B^\z}. \label{eq:stdot}
!latex       \ee
     
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \end{enumerate} \subsection{representation of magnetic field} \begin{enumerate}

!latex \item \Ais
!latex \item \Bis
!latex \item In \Eqn{stdot}, the coordinate Jacobian, $\sqrt g$, cancels.
!latex       No coordinate metric information is required to construct the fieldline equations from the magnetic vector potential.

!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! IT IS REQUIRED TO SET IVOL THROUGH GLOBAL MEMORY BEFORE CALLING BF00AA;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine bfield( zeta, st, Bst ) ! the format of this subroutine is constrained by the NAG ode integration routines. 
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, one, half, two
  
  use numerical, only : vsmall, small
  
  use fileunits, only : ounit
  
  use inputlist, only : Wmacros, Wbfield, Lrad, Mpol
  
  use cputiming, only : Tbfield
  
  use allglobal, only : myid, ncpu, cpus, mn, im, in, halfmm, regumm, &
                        ivol, gBzeta, Ate, Aze, Ato, Azo, &
                        NOTstellsym, &
                        Lcoordinatesingularity, Mvol, &
                        Node ! 17 Dec 15;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  REAL, intent(in)   :: zeta,  st(1:Node)
  REAL, intent(out)  ::       Bst(1:Node)
  
  INTEGER            :: lvol, ii, ll, mi, ni, ideriv
  REAL               :: teta, lss, sbar, sbarhm(0:1), arg, carg, sarg, dBu(1:3)
  REAL               :: cheby(0:Lrad(ivol),0:1), zernike(0:Lrad(1),0:Mpol,0:1)
  
  REAL               :: TT(0:Lrad(ivol),0:1) ! this is almost identical to cheby; 17 Dec 15;

  BEGIN(bfield)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  FATAL( bfield, ivol.lt.1 .or. ivol.gt.Mvol, invalid ivol )
#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  lvol = ivol ; ideriv = 0 ! the argument list of bfield is fixed by NAG requirements, but volume index is required below;

  Bst(1:Node) = (/ zero , zero /) ! set default intent out; this should cause a compilation error if Node.ne.2;
 
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  lss = st(1) ; teta = st(2) 
  
  if( abs(lss).gt.one ) goto 9999 ! out of domain;
    
  if( Lcoordinatesingularity ) sbar = max( ( lss + one ) * half, small )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if (Lcoordinatesingularity) then
    call get_zernike(sbar, Lrad(lvol), Mpol, zernike(:,:,0:1))
  else
    call get_cheby(lss, Lrad(lvol), cheby(0:Lrad(lvol),0:1))
  end if

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  dBu(1:3) = zero ! initialize summation;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  do ii = 1, mn ; mi = im(ii) ; ni = in(ii) ; arg = mi * teta - ni * zeta ; carg = cos(arg) ; sarg = sin(arg) ! shorthand; 20 Apr 13;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   if( Lcoordinatesingularity ) then ! regularization factor depends on mi; 17 Dec 15;

    FATAL( bfield, abs(sbar).lt.vsmall, need to avoid divide-by-zero )

    do ll = 0, Lrad(lvol) ; TT(ll,0:1) = (/        zernike(ll,mi,0),        zernike(ll,mi,1)*half                      /)
    enddo

   else

    do ll = 0, Lrad(lvol) ; TT(ll,0:1) = (/             cheby(ll,0),             cheby(ll,1)                           /)
    enddo

   endif ! end of if( Lcoordinatesingularity ) ; 16 Jan 15;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   do ll = 0, Lrad(lvol) ! loop over Chebyshev summation; 20 Feb 13;
    ;dBu(1) = dBu(1) + ( - mi * Aze(lvol,ideriv,ii)%s(ll) - ni * Ate(lvol,ideriv,ii)%s(ll) ) * TT(ll,0) * sarg
    ;dBu(2) = dBu(2) + (                                  -      Aze(lvol,ideriv,ii)%s(ll) ) * TT(ll,1) * carg
    ;dBu(3) = dBu(3) + (        Ate(lvol,ideriv,ii)%s(ll)                                  ) * TT(ll,1) * carg
    if( NOTstellsym ) then ! include non-symmetric harmonics; 28 Jan 13;
     dBu(1) = dBu(1) + ( + mi * Azo(lvol,ideriv,ii)%s(ll) + ni * Ato(lvol,ideriv,ii)%s(ll) ) * TT(ll,0) * carg
     dBu(2) = dBu(2) + (                                  -      Azo(lvol,ideriv,ii)%s(ll) ) * TT(ll,1) * sarg
     dBu(3) = dBu(3) + (        Ato(lvol,ideriv,ii)%s(ll)                                  ) * TT(ll,1) * sarg
    endif
   enddo ! end of do ll; 10 Dec 15;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
  enddo ! end of do ii = 1, mn;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  gBzeta = dBu(3) ! gBzeta is returned through global; 20 Apr 13;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( abs(gBzeta).lt.vsmall ) then

   cput = GETTIME

   write(ounit,'("bfield : ",f10.2," : lvol=",i3," ; zeta="es23.15" ; (s,t)=("es23.15" ,"es23.15" ) ; B^z="es23.15" ;")') &
                             cput-cpus, lvol,        zeta,             st(1:2),                       dBu(3)

   FATAL( bfield, abs(dBu(3)).lt.vsmall, field is not toroidal )

  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  Bst(1:2) = dBu(1:2) / gBzeta ! normalize field line equations to toroidal field; 20 Apr 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(bfield)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine bfield

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
