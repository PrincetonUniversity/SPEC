!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Constructs spectrally-condensed Fourier representation of interfaces using a stream function.

!latex \subsubsection{angle transformation}

!latex \item The geometry of each interface is given (on input) as
!latex \be \begin{array}{ccccccccccccccccccccc}
!latex     R(\theta,\zeta) & = &\sum_i R_i \cos(m_i \theta - n_i \zeta), \\
!latex     Z(\theta,\zeta) & = &\sum_i Z_i \sin(m_i \theta - n_i \zeta).
!latex \end{array} \label{eq:initialgeometry}\ee

!latex \item A new angle, $\bar \t$, shall be introduced via a stream function, $\lambda(\bar \t,\z)$, according to 
!latex \be \theta &=&\bar \theta + \sum_j \lambda_j \sin(m_j \bar \theta - n_j \zeta), \label{eq:angletransformation} \ee
!latex where the $\lambda_j$ are, as yet, unknown degrees of freedom.

!latex \item The Fourier harmonics in the new angle are
!latex \be \bar R_k &=& \oint \!\!\! \oint \! d\bar \theta d \zeta \; R \cos(m_k \bar \theta - n_k \zeta),\\
!latex     \bar Z_k &=& \oint \!\!\! \oint \! d\bar \theta d \zeta \; Z \sin(m_k \bar \theta - n_k \zeta),
!latex \ee
!latex where, by combining \Eqn{initialgeometry} and \Eqn{angletransformation}, it is understood that $R\equiv R(\bar \t,\z)$and $Z\equiv Z(\bar \t,\z)$.

!latex \item The spectral-width (in the new angle) is defined
!latex \be M = \frac{1}{2}\sum_k (m_k^p+n_k^q) \left(\bar R_k^2 + \bar Z_k^2\right),\ee
!latex where $m_k^p = 0$ for $m_k=0$, and $n_k^q = 0$ for $n_k=0$, and where $p\equiv$\verb+pwidth+ and $q\equiv$\verb+qwidth+ are given on input.

!latex \item The variation in spectral-width due to variations, $\delta \lambda_j$ is
!latex \be \frac{\partial M}{\partial \lambda_j} = 
!latex     \sum_k (m_k^p+n_k^q)\left( \bar R_k\frac{\partial \bar R_k}{\partial \lambda_j} + \bar Z_k\frac{\partial \bar Z_k}{\partial \lambda_j} \right). 
!latex \ee

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

module spectralm

  REAL, allocatable :: dRmndl(:,:), dZmndl(:,:), dRjkdl(:,:), dZjkdl(:,:)

end module spectralm

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine sw00ac( lvol, mn, Ntz, Rbc, Zbs, Rbs, Zbc )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, one

  use numerical, only : machprec, sqrtmachprec

  use fileunits, only : ounit

  use inputlist, only : Wmacros, Wsw00ac

  use cputiming, only : Tsw00ac

  use allglobal, only : myid, cpus, im, in, lRbc, lZbs, lRbs, lZbc, xoffset, &
                        NOTstellsym

  use spectralm
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS

  INTEGER, intent(in) :: lvol, mn, Ntz
  REAL                :: Rbc(1:mn), Zbs(1:mn), Rbs(1:mn), Zbc(1:mn)
  
  INTEGER             :: Ndof, irevcm, ML, MU, mode, LDfjac, LR, ic05ndf, nc05ndf

  REAL                :: lambda(2:mn), specwidth, dMdl(2:mn), xtol, dR, dZ
  REAL                :: epsfcn, diag(2:mn), factor, fjac(2:mn,2:mn), RR( 1:mn*(mn+1)/2 ), QTF(2:mn), work(2:mn,1:4)
  
! INTEGER             :: LWA
! REAL                :: FVEC(2:mn), WA((mn-1)*(3*(mn-1)+13)/2)
! external            :: FCN

! INTEGER             :: ITER, IWORK(N+1), IUSER(1:1), IFAIL
! real                :: OBJF, OBJGRD(N), X(N), WORK(13*N), USER(1:1)
! external            :: OBJFUN
  
  BEGIN(sw00ac)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  FATALMESS(sw00ac, NOTstellsym, under construction )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  Ndof = mn-1 ! #degrees-of-freedom; \lambda_1 = \lambda_{m=0,n=0} is irrelevant;

  lRbc(1:mn) = Rbc(1:mn) ! this is passed through global;
  lZbs(1:mn) = Zbs(1:mn)
  lRbs(1:mn) = Rbs(1:mn) ! this is passed through global;
  lZbc(1:mn) = Zbc(1:mn)
  
  RALLOCATE(dRjkdl,(1:Ntz,1:mn))
  RALLOCATE(dZjkdl,(1:Ntz,1:mn))
  
  RALLOCATE(dRmndl,(1:mn,1:mn))
  RALLOCATE(dZmndl,(1:mn,1:mn))

!latex \subsubsection{numerical implementation}
  
!latex \item This routine seeks a zero of a vector function, ${\bf F}(\boldlambda)$, where $F_j \equiv \partial M / \partial \lambda_j$.
!latex \item The NAG routine \verb+c05nbf+ is employed
!latex (This routine uses function values only: perhaps the derivatives could be calculated and more efficient routines enabled.)

!latex \item It is probably preferable to use \verb+E04LYF+.
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!  lambda(2:mn) = offset ; xtol = machprec ; LWA = (mn-1)*(3*(mn-1)+13)/2 ; ic05nbf = 1
!  call C05NBF( FCN, Ndof, lambda(2:mn), FVEC(2:mn), XTOL, WA(1:LWA), LWA, ic05nbf ) ! NEEDS FCN TO BE DEFINED; FCN NEEDS TO BE CONSISTENT WITH NAG;
!
!  cput = GETTIME
!  select case( ic05nbf )
!  case(  0 )   ; write(ounit,'("sw00ac : ",f10.2," : myid=",i3," ; lvol=",i3," : ic05nbf=",i3," ; success ;         ")')cput-cpus,myid,lvol,ic05nbf
!  case( -1 )   ; write(ounit,'("sw00ac : ",f10.2," : myid=",i3," ; lvol=",i3," : ic05nbf=",i3," ; terminated ;      ")')cput-cpus,myid,lvol,ic05nbf
!  case(  1 )   ; write(ounit,'("sw00ac : ",f10.2," : myid=",i3," ; lvol=",i3," : ic05nbf=",i3," ; input error ;     ")')cput-cpus,myid,lvol,ic05nbf
!  case(  2 )   ; write(ounit,'("sw00ac : ",f10.2," : myid=",i3," ; lvol=",i3," : ic05nbf=",i3," ; consider restart ;")')cput-cpus,myid,lvol,ic05nbf
!  case(  3 )   ; write(ounit,'("sw00ac : ",f10.2," : myid=",i3," ; lvol=",i3," : ic05nbf=",i3," ; xtol too small ;  ")')cput-cpus,myid,lvol,ic05nbf
!  case(  4 )   ; write(ounit,'("sw00ac : ",f10.2," : myid=",i3," ; lvol=",i3," : ic05nbf=",i3," ; bad progress ;    ")')cput-cpus,myid,lvol,ic05nbf
!  case default ; write(ounit,'("sw00ac : ",f10.2," : myid=",i3," ; lvol=",i3," : ic05nbf=",i3," ; illegal ifail ;   ")')cput-cpus,myid,lvol,ic05nbf
!  end select

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!  ie04dgf = 1
!  call E04DGF( OBJFUN, ITER, OBJF, OBJGRD, lambda(2:mn), IWORK(2:mn+1), WORK(1:, IUSER(1:1), USER(1:1), ie04dgf ) ! THIS IS NOT THREAD SAFE;
!
!  cput = GETTIME
!  select case( ie04dgf )
!  case(  0 )   ; write(ounit,'("sw00ac : ",f10.2," : myid=",i3," ; lvol=",i3," : ie04dgf=",i3," ; success ;               ")')cput-cpus,myid,lvol,ie04dgf
!  case( -1 )   ; write(ounit,'("sw00ac : ",f10.2," : myid=",i3," ; lvol=",i3," : ie04dgf=",i3," ; terminated ;            ")')cput-cpus,myid,lvol,ie04dgf
!  case(  3 )   ; write(ounit,'("sw00ac : ",f10.2," : myid=",i3," ; lvol=",i3," : ie04dgf=",i3," ; iteration limit ;       ")')cput-cpus,myid,lvol,ie04dgf
!  case(  4 )   ; write(ounit,'("sw00ac : ",f10.2," : myid=",i3," ; lvol=",i3," : ie04dgf=",i3," ; small step length ;     ")')cput-cpus,myid,lvol,ie04dgf
!  case(  6 )   ; write(ounit,'("sw00ac : ",f10.2," : myid=",i3," ; lvol=",i3," : ie04dgf=",i3," ; almost satisfactory ;   ")')cput-cpus,myid,lvol,ie04dgf
!  case(  7 )   ; write(ounit,'("sw00ac : ",f10.2," : myid=",i3," ; lvol=",i3," : ie04dgf=",i3," ; gradient error ;        ")')cput-cpus,myid,lvol,ie04dgf
!  case(  8 )   ; write(ounit,'("sw00ac : ",f10.2," : myid=",i3," ; lvol=",i3," : ie04dgf=",i3," ; initial gradient error ;")')cput-cpus,myid,lvol,ie04dgf
!  case(  9 )   ; write(ounit,'("sw00ac : ",f10.2," : myid=",i3," ; lvol=",i3," : ie04dgf=",i3," ; input error ;           ")')cput-cpus,myid,lvol,ie04dgf
!  case default ; write(ounit,'("sw00ac : ",f10.2," : myid=",i3," ; lvol=",i3," : ie04dgf=",i3," ; illegal ifail ;         ")')cput-cpus,myid,lvol,ie04dgf
!  end select

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  lambda(2:mn) = xoffset ; xtol = machprec ; ML = Ndof-1 ; MU = Ndof-1 ; epsfcn = zero
  
  mode = 0 ; diag(2:mn) = zero ; factor = one

  LDfjac = Ndof ; fjac(2:mn,2:mn) = zero ; RR(1:mn*(mn+1)/2) = zero ; LR = mn*(mn+1)/2 ; QTF(2:mn) = zero ; work(2:mn,1:4) = zero
  
  irevcm = 0 ; ic05ndf = 1 ; nc05ndf = 0
  
  do ! reverse communication loop;
   
   call C05NDF( irevcm, Ndof, lambda(2:mn), dMdl(2:mn), xtol, ML, MU, epsfcn, diag(2:mn), mode, factor, fjac(2:mn,2:mn), LDfjac, &
                RR(1:mn*(mn+1)/2), LR, QTF(2:mn), work(2:mn,1:4), ic05ndf )
   
   select case(irevcm)
    
   case( 0 ) ! final termination;
    
    cput = GETTIME
    
    select case( ic05ndf )                                                                                    !1234567890123456
    case( 0 ) ; if( Wsw00ac) write(ounit,1002) cput-cpus, myid, lvol, ic05ndf, sqrt(sum(dMdl(2:mn)**2)/Ndof), "success ;       "
    case( 1 ) ;              write(ounit,1002) cput-cpus, myid, lvol, ic05ndf, sqrt(sum(dMdl(2:mn)**2)/Ndof), "input error ;   "
    case( 2 ) ;              write(ounit,1002) cput-cpus, myid, lvol, ic05ndf, sqrt(sum(dMdl(2:mn)**2)/Ndof), "irevcm error ;  "
    case( 3 ) ;              write(ounit,1002) cput-cpus, myid, lvol, ic05ndf, sqrt(sum(dMdl(2:mn)**2)/Ndof), "xtol too small ;"
    case( 4 ) ;              write(ounit,1002) cput-cpus, myid, lvol, ic05ndf, sqrt(sum(dMdl(2:mn)**2)/Ndof), "bad progress ;  "
    case( 5 ) ;              write(ounit,1002) cput-cpus, myid, lvol, ic05ndf, sqrt(sum(dMdl(2:mn)**2)/Ndof), "bad progress ;  "
    case default
     FATALMESS(sw00ac,.true.,invalid ic05ndf)
    end select
    
1002 format("sw00ac : ",f10.2," : myid=",i3," ; lvol=",i3," ; ic05ndf="i2" ; |dMdl|="es13.5" ; "a16)
    
    nc05ndf = nc05ndf + 0 ; call specwidthderivatives( mn, lambda(2:mn), specwidth, dMdl(2:mn) ) ! to make sure best guess is dRmndl(1:mn,1), dZmndl(1:mn,1)

    if( Wsw00ac ) then
     cput = GETTIME
     if( Ndof.gt. 14 ) write(ounit,1003) cput-cpus, myid, lvol, Ndof, nc05ndf, ic05ndf, specwidth, sqrt(sum( (lambda(2:mn)-xoffset)**2 )/Ndof)
     if( Ndof.le. 14 ) write(ounit,1004) cput-cpus, myid, lvol, Ndof, nc05ndf, ic05ndf, specwidth,            lambda(2:mn)-xoffset
    endif
    
    exit !! exit reverse communication iterations; finished; 15 May 13;
    
   case( 1 ) ! no action required;
    
   case( 2 ) ! return function values;
    
    nc05ndf = nc05ndf + 1 ; call specwidthderivatives( mn, lambda(2:mn), specwidth, dMdl(2:mn) )
    
   case default
    
    FATALMESS(sw00ac, .true., invalid irevcm returned by C05NDF )
    
   end select
   
  enddo ! end of reverse communication loop;
  
1003 format("sw00ac : ",f10.2," : myid=",i3," ; lvol=",i3," ; Ndof="i4" ; nc05ndf="i6" ; ic05ndf="i2" ; M="es23.15" ;":" |lambda|="es10.2" ;")
1004 format("sw00ac : ",f10.2," : myid=",i3," ; lvol=",i3," ; Ndof="i4" ; nc05ndf="i6" ; ic05ndf="i2" ; M="es23.15" ;":"  lambda ="999es10.2" ;")
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  DEALLOCATE(dRjkdl)
  DEALLOCATE(dZjkdl)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Wsw00ac ) then
   cput = GETTIME
   dR = sqrt( sum(( Rbc(1:mn)-dRmndl(1:mn,1) )**2) /  mn    ) ! this measures how much the Fourier representation has changed from the initial;
   dZ = sqrt( sum(( Zbs(2:mn)-dZmndl(2:mn,1) )**2) / (mn-1) )
   write(ounit,1005) cput-cpus, myid, lvol, ic05ndf, mn, sqrt(sum(dMdl**2)/Ndof), dR, dZ, cput-cpui
  endif
  
1005 format("sw00ac : ",f10.2," : myid=",i3," ; lvol=",i3," ; ic05ndf="i2" ; mn="i4" ; |dMdl|="es9.2" ; dR,dZ="2es9.2" ; time=",f10.2,"s ;")
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  cput = GETTIME
  if( maxval(abs(dMdl(2:mn))) .lt. sqrtmachprec ) then !latex \item Condensed representation only accepted if $|\partial M / \partial \lambda_j|<$\verb+small+.
   if( Wsw00ac ) write(ounit,'("sw00ac : ",f10.2," : myid=",i3," ; lvol=",i3," ; updating Fourier harmonics ;        ")') cput-cpus, myid, lvol
   Rbc(1:mn) = dRmndl(1:mn,1) ! update Fourier representation; note that the Fourier representations are calculated in specwidthderivatives;
   Zbs(2:mn) = dZmndl(2:mn,1)
  else
   ;             write(ounit,'("sw00ac : ",f10.2," : myid=",i3," ; lvol=",i3," ; updating Fourier harmonics ; FAIL ; |dMdl|="es10.2" ;")') &
cput-cpus, myid, lvol, maxval(abs(dMdl(2:mn)))
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  DEALLOCATE(dRmndl)
  DEALLOCATE(dZmndl)
  
  if( Wsw00ac ) then
   cput = GETTIME
   write(ounit,'("sw00ac : ",f10.2," : Rbc="999f15.11)')cput-cpus,Rbc(1:mn)
   write(ounit,'("sw00ac : ",f10.2," : Zbs="999f15.11)')cput-cpus,Zbs(1:mn)
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(sw00ac)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine sw00ac

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine specwidthderivatives( mn, lambda, specwidth, dMdl )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, pi2

  use numerical, only :

  use fileunits, only : ounit

  use inputlist, only : Wmacros

  use allglobal, only : myid, cpus, pi2nfp, &
                        im, in, mpnq, & ! nq, & ! 18 Jul 14;
                        lRbc, lZbs, &
                        Nt, Nz, Ntz, isr, trigm, trign, trigwk, efmn, ofmn, cfmn, sfmn, ijreal, ijimag, xoffset

  use spectralm
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in)    :: mn
  REAL   , intent(in)    :: lambda(2:mn)
  REAL   , intent(out)   :: specwidth, dMdl(2:mn)
  
  INTEGER                :: imn, jmn, kmn, jj, kk, jk, ifail
  REAL                   :: tb, tt, zz, arg, ca, sa, barg, sba
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! construct lambda on regular \bar\theta grid;
  efmn(1:mn) = zero ; ofmn(1:mn) = (/ zero, lambda(2:mn) - xoffset /) ; cfmn(1:mn) = zero ; sfmn(1:mn) = zero 
  
  call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, ijreal(1:Ntz), ijimag(1:Ntz), isr, trigm, trign, trigwk )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! construct R, Z (and their derivatives) on regular \bar \theta grid, which is an irregular \theta grid and so must be solved using a slow transform;
  
  dRjkdl(1:Ntz,1:mn) = zero
  dZjkdl(1:Ntz,1:mn) = zero
  
  do kk = 0,Nz-1 ; zz = kk * pi2nfp / Nz
   do jj = 0,Nt-1 ; jk = 1 + jj + kk*Nt ; tb = jj * pi2 / Nt ; tt = tb + ijreal(jk)
    
    do imn = 1,mn ; arg = im(imn)*tt - in(imn)*zz ; ca = cos(arg) ; sa = sin(arg) ! slow transform as tt is not regular;

! need to include mask;
     
     dRjkdl(jk, 1 ) = dRjkdl(jk, 1 ) + lRbc(imn) * ca ! R;
     dZjkdl(jk, 1 ) = dZjkdl(jk, 1 ) + lZbs(imn) * sa ! Z;
     
     do jmn = 2,mn ; barg = im(jmn)*tb - in(jmn)*zz ; sba = sin(barg) ! NOTE: lambda MAY be higher resolution, 
      
!latex \item Differentiating the Fourier harmonic $\bar R_k$ with respect to $\lambda_j$ is equivalent to Fourier decomposing the derivative:
!latex \be \frac{\partial \bar R_k}{\partial \lambda_j} & \equiv & \left(\frac{\partial \bar R}{\partial \lambda_j}\right)_k, \\
!latex     \frac{\partial \bar Z_k}{\partial \lambda_j} & \equiv & \left(\frac{\partial \bar Z}{\partial \lambda_j}\right)_k. \ee
      
      dRjkdl(jk,jmn) = dRjkdl(jk,jmn) - lRbc(imn) * sa * im(imn) * sba ! d R / d l_j
      dZjkdl(jk,jmn) = dZjkdl(jk,jmn) + lZbs(imn) * ca * im(imn) * sba ! d Z / d l_j
      
     enddo ! end of do jmn = 2,mn
     
    enddo ! end of do imn = 1,mn
    
   enddo ! end of do jj = 0,Nt-1
  enddo ! end of do kk = 0 Nz-1
  
! now need to combine this information to get derivatives of spectral-width wrt lambda harmonics; these Fourier decompositions are reduced resolution;
  
  do jmn = 1, mn
   call tfft( Nt, Nz, dRjkdl(1:Ntz,jmn), dZjkdl(1:Ntz,jmn), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz), &
              mn, im(1:mn), in(1:mn), dRmndl(1:mn,jmn), ofmn(1:mn), cfmn(1:mn), dZmndl(1:mn,jmn), ifail )
  enddo ! end of do jmn = 1,mn
  

  specwidth = zero ; dMdl(2:mn) = zero
  
  do jmn = 2, mn ! loop over non-trivial lambda harmonics;
   
   specwidth = specwidth + ( mpnq(jmn)         ) * ( dRmndl(jmn,1)**2 + dZmndl(jmn,1)**2 )

   do kmn = 2, mn

! derivatives of spectral-width wrt m.ne.0 lambda harmonics;    
    dMdl(jmn) = dMdl(jmn) + ( mpnq(kmn)         ) * ( dRmndl(kmn,1) * dRmndl(kmn,jmn) + dZmndl(kmn,1) * dZmndl(kmn,jmn) )

   enddo
   
  enddo
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  return
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine specwidthderivatives

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

