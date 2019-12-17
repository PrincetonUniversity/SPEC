!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (coordinate axis) ! Specifies position of coordinate axis; ${\bf x}_a(\zeta) \equiv \int {\bf x}_1(\theta,\zeta) dl \, / \int dl$.

!latex \briefly{The coordinate axis is assigned via a poloidal average over an arbitrary surface.}

!latex \calledby{\link{preset}, \link{packxi}}

!latex \tableofcontents

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{coordinate axis}

!latex \begin{enumerate}
!latex \item The coordinate axis is {\em not} an independent degree-of-freedom of the geometry.
!latex       It is constructed by extrapolating the geometry of a given interface, as determined by $i \equiv $ \internal{ivol} which is given on input,
!latex       down to a line.
!latex \item If the coordinate axis depends only on the {\em geometry} of the interface and not the angle parameterization,
!latex       then the block tri-diagonal structure of the the force-derivative matrix is preserved.
!latex \item Define the arc-length-weighted averages,
!latex       \be R_0(\z) \equiv \frac{\ds \int_{0}^{2\pi} R_i(\t,\z) \, dl}{\ds \int_{0}^{2\pi} \!\!\!\! dl}, \qquad
!latex           Z_0(\z) \equiv \frac{\ds \int_{0}^{2\pi} Z_i(\t,\z) \, dl}{\ds \int_{0}^{2\pi} \!\!\!\! dl},
!latex       \ee
!latex       where $dl \equiv \dot l \, d\t = \sqrt{ \partial_\t R_i(\t,\z)^2 + \partial_\t Z_i(\t,\z)^2 } \, d\t$.
!latex \item (Note that if $\dot l$ does not depend on $\t$, i.e. if $\t$ is the equal arc-length angle, then the expressions simplify.
!latex        This constraint is not enforced.)
!latex \item The geometry of the coordinate axis thus constructed only depends on the geometry of the interface, i.e. 
!latex       the angular parameterization of the interface is irrelevant.
!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{coordinate axis: derivatives}

!latex \begin{enumerate}
!latex \item The derivatives of the coordinate axis with respect to the Fourier harmonics of the given interface are given by
!latex       \be
!latex       \ds \frac{\partial R_0}{\partial R_{i,j}^c} & = & \ds \int \left( \cos\a_j \; \dot l
!latex                                                     -       \Delta R_i R_{i,\t} \, m_j \sin\a_j / \; \dot l \right) d\t / L \\
!latex       \ds \frac{\partial R_0}{\partial R_{i,j}^s} & = & \ds \int \left( \sin\a_j \; \dot l
!latex                                                     +       \Delta R_i R_{i,\t} \, m_j \cos\a_j / \; \dot l \right) d\t / L \\
!latex       \ds \frac{\partial R_0}{\partial Z_{i,j}^c} & = & \ds \int \left( \;\;\;\;\;\;\;\;\;\;\;\;\,
!latex                                                     -       \Delta R_i Z_{i,\t} \, m_j \sin\a_j / \; \dot l \right) d\t / L \\
!latex       \ds \frac{\partial R_0}{\partial Z_{i,j}^s} & = & \ds \int \left( \;\;\;\;\;\;\;\;\;\;\;\;\,                             
!latex                                                     +       \Delta R_i Z_{i,\t} \, m_j \cos\a_j / \; \dot l \right) d\t / L \\ \nonumber \\
!latex       \ds \frac{\partial Z_0}{\partial R_{i,j}^c} & = & \ds \int \left( \;\;\;\;\;\;\;\;\;\;\;\;\,                             
!latex                                                     -       \Delta Z_i R_{i,\t} \, m_j \sin\a_j / \; \dot l \right) d\t / L \\
!latex       \ds \frac{\partial Z_0}{\partial R_{i,j}^s} & = & \ds \int \left( \;\;\;\;\;\;\;\;\;\;\;\;                            
!latex                                                     +       \Delta Z_i R_{i,\t} \, m_j \cos\a_j / \; \dot l \right) d\t / L \\
!latex       \ds \frac{\partial Z_0}{\partial Z_{i,j}^c} & = & \ds \int \left( \cos\a_j \; \dot l
!latex                                                     -       \Delta Z_i Z_{i,\t} \, m_j \sin\a_j / \; \dot l \right) d\t / L \\
!latex       \ds \frac{\partial Z_0}{\partial Z_{i,j}^s} & = & \ds \int \left( \sin\a_j \; \dot l
!latex                                                     +       \Delta Z_i Z_{i,\t} \, m_j \cos\a_j / \; \dot l \right) d\t / L
!latex       \ee
!latex       where $\ds L(\z) \equiv \int_{0}^{2\pi} \!\!\!\! dl$.
!latex \end{enumerate}

!latex \subsection{some numerical comments}

!latex \begin{enumerate}
!latex \item First, the differential poloidal length, $\dot l \equiv \sqrt{ R_\t^2 + Z_\t^2 }$, is computed in real space using 
!latex       an inverse FFT from the Fourier harmonics of $R$ and $Z$.
!latex \item Second, the Fourier harmonics of $dl$ are computed using an FFT.
!latex       The integration over $\t$ to construct $L\equiv \int dl$ is now trivial: just multiply the $m=0$ harmonics of $dl$ by $2\pi$.
!latex       The \internal{ajk(1:mn)} variable is used, and this is assigned in \link{global}.
!latex \item Next, the weighted $R \, dl$ and $Z \, dl$ are computed in real space, and the poloidal integral is similarly taken.
!latex \item Last, the Fourier harmonics are constructed using an FFT after dividing in real space.
!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine rzaxis( Mvol, mn, iRbc, iZbs, iRbs, iZbc, ivol )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero
  
  use numerical, only :
  
  use fileunits, only : ounit
  
  use inputlist, only : Wrzaxis, Igeometry, Ntor
  
  use cputiming, only : Trzaxis
  
  use allglobal, only : ncpu, myid, cpus, im, in, &
                        ajk, Nt, Nz, Ntz, &
                        ijreal, ijimag, jireal, jiimag, jkreal, jkimag, kjreal, kjimag, &
                        efmn, ofmn, cfmn, sfmn, evmn, odmn, comn, simn, cosi, sini, &
                        NOTstellsym, &
                        dRodR, dRodZ, dZodR, dZodZ, &
                        dRadR, dRadZ, dZadR, dZadZ
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in)    :: Mvol, mn, ivol
  REAL                   :: iRbc(1:mn,0:Mvol), iZbs(1:mn,0:Mvol), iRbs(1:mn,0:Mvol), iZbc(1:mn,0:Mvol)
  
  INTEGER                :: jvol, ii, ifail
  
  BEGIN(rzaxis)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef DEBUG
  FATAL( rzaxis, ivol.gt.Mvol, perhaps illegal combination Linitialize=2 and Lfreebound=0 )
#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  jvol = 0 ! this identifies the "surface" in which the poloidal averaged harmonics will be placed; 19 Jul 16; 

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  select case( Igeometry ) 
   
  case( 1:2 )
   
   iRbc(1:mn,jvol) = zero
   iRbs(1:mn,jvol) = zero
   
  case(   3 )
   
   call invfft( mn, im(1:mn), in(1:mn), im(1:mn) * iRbs(1:mn,ivol), - im(1:mn) * iRbc(1:mn,ivol), &
                                        im(1:mn) * iZbs(1:mn,ivol), - im(1:mn) * iZbc(1:mn,ivol), &
                Nt, Nz, jkreal(1:Ntz), jkimag(1:Ntz) ) ! R_\t, Z_\t; 03 Nov 16;
 
   ijreal(1:Ntz) = sqrt( jkreal(1:Ntz)**2 + jkimag(1:Ntz)**2 ) ! dl ; 11 Aug 14;
   ijimag(1:Ntz) = zero

   jireal(1:Ntz) = ijreal(1:Ntz) ! dl ; 19 Sep 16;
   
   ifail = 0
   call tfft( Nt, Nz, ijreal(1:Ntz), ijimag(1:Ntz), &
              mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), ifail ) ! Fourier harmonics of differential poloidal length; 11 Mar 16;

   efmn(1:mn) = efmn(1:mn) * ajk(1:mn) ! poloidal integration of length; only take m=0 harmonics; 11 Aug 14;
   ofmn(1:mn) = ofmn(1:mn) * ajk(1:mn)
   cfmn(1:mn) = zero
   sfmn(1:mn) = zero
      
   call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), & ! map length = "integrated dl" back to real space; 19 Sep 16;
                Nt, Nz, ijreal(1:Ntz), ijimag(1:Ntz) )
    
   jiimag(1:Ntz) = ijreal(1:Ntz) !  L ; 19 Sep 16;


   call invfft( mn, im(1:mn), in(1:mn),            iRbc(1:mn,ivol),              iRbs(1:mn,ivol), &
                                                   iZbc(1:mn,ivol),              iZbs(1:mn,ivol), &
                Nt, Nz, kjreal(1:Ntz), kjimag(1:Ntz) ) ! R, Z; 03 Nov 16;
   
   ijreal(1:Ntz) = kjreal(1:Ntz) * jireal(1:Ntz) ! R dl;
   ijimag(1:Ntz) = kjimag(1:Ntz) * jireal(1:Ntz) ! Z dl;
   
   ifail = 0
   call tfft( Nt, Nz, ijreal(1:Ntz), ijimag(1:Ntz), &
              mn, im(1:mn), in(1:mn), evmn(1:mn), odmn(1:mn), comn(1:mn), simn(1:mn), ifail ) ! Fourier harmonics of weighted R & Z; 11 Mar 16;

   evmn(1:mn) = evmn(1:mn) * ajk(1:mn) ! poloidal integration of R dl; 19 Sep 16;
   odmn(1:mn) = odmn(1:mn) * ajk(1:mn)
   comn(1:mn) = comn(1:mn) * ajk(1:mn) ! poloidal integration of Z dl; 19 Sep 16;
   simn(1:mn) = simn(1:mn) * ajk(1:mn)
   
   call invfft( mn, im(1:mn), in(1:mn), evmn(1:mn), odmn(1:mn), comn(1:mn), simn(1:mn), &
                Nt, Nz, ijreal(1:Ntz), ijimag(1:Ntz) )
 
   ijreal(1:Ntz) = ijreal(1:Ntz) / jiimag(1:Ntz) ! Ro; 19 Sep 16;
   ijimag(1:Ntz) = ijimag(1:Ntz) / jiimag(1:Ntz) ! Zo; 19 Sep 16;
   
   kjreal(1:Ntz) = kjreal(1:Ntz) - ijreal(1:Ntz) ! \Delta R = R_1 - R_0 ; 03 Nov 16;
   kjimag(1:Ntz) = kjimag(1:Ntz) - ijimag(1:Ntz) ! \Delta R = Z_1 - Z_0 ; 03 Nov 16;

   ifail = 0
   call tfft( Nt, Nz, ijreal(1:Ntz), ijimag(1:Ntz), &
              mn, im(1:mn), in(1:mn), iRbc(1:mn,jvol), iRbs(1:mn,jvol), iZbc(1:mn,jvol), iZbs(1:mn,jvol), ifail )
   
#ifdef DEBUG
   if( Wrzaxis ) then
    cput = GETTIME
    write(ounit,'("rzaxis : ", 10x ," : ")')
    write(ounit,'("rzaxis : ",f10.2," : myid=",i3," ; inner : Rbc=[",    999(es23.15," ,"))') cput-cpus, myid, iRbc(1:Ntor+1,ivol)
    write(ounit,'("rzaxis : ",f10.2," : myid=",i3," ; axis  : Rbc=[",    999(es23.15," ,"))') cput-cpus, myid, iRbc(1:Ntor+1,jvol)
    if( Ntor.gt.0 ) then
    write(ounit,'("rzaxis : ",f10.2," : myid=",i3," ; inner : Zbs=[",25x,998(es23.15," ,"))') cput-cpus, myid, iZbs(2:Ntor+1,ivol)
    write(ounit,'("rzaxis : ",f10.2," : myid=",i3," ; axis  : Zbs=[",25x,998(es23.15," ,"))') cput-cpus, myid, iZbs(2:Ntor+1,jvol)
    endif
    if( NOTstellsym ) then
    if( Ntor.gt.0 ) then
    write(ounit,'("rzaxis : ",f10.2," : myid=",i3," ; inner : Rbs=[",25x,998(es23.15," ,"))') cput-cpus, myid, iRbs(2:Ntor+1,ivol)
    write(ounit,'("rzaxis : ",f10.2," : myid=",i3," ; axis  : Rbs=[",25x,998(es23.15," ,"))') cput-cpus, myid, iRbs(2:Ntor+1,jvol)
    endif
    write(ounit,'("rzaxis : ",f10.2," : myid=",i3," ; inner : Zbc=[",    999(es23.15," ,"))') cput-cpus, myid, iZbc(1:Ntor+1,ivol)
    write(ounit,'("rzaxis : ",f10.2," : myid=",i3," ; axis  : Zbc=[",    999(es23.15," ,"))') cput-cpus, myid, iZbc(1:Ntor+1,jvol)
    endif
   endif
#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef DEBUG
   FATAL( rzaxis, .not.allocated(cosi), fatal )
   FATAL( rzaxis, .not.allocated(sini), fatal )
#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! compute derivatives of axis; 03 Nov 16;

   do ii = 1, mn

    
    dRodR(1:Ntz,0,ii) = cosi(1:Ntz,ii) * jireal(1:Ntz) - kjreal(1:Ntz) * jkreal(1:Ntz) * im(ii) * sini(1:Ntz,ii) / jireal(1:Ntz) ! dRodRjc;
    dRodR(1:Ntz,1,ii) = sini(1:Ntz,ii) * jireal(1:Ntz) + kjreal(1:Ntz) * jkreal(1:Ntz) * im(ii) * cosi(1:Ntz,ii) / jireal(1:Ntz) ! dRodRjs;
    
    ifail = 0
    call tfft( Nt, Nz, dRodR(1:Ntz,0,ii), dRodR(1:Ntz,1,ii), &
              mn, im(1:mn), in(1:mn), dRadR(1:mn,0,0,ii), dRadR(1:mn,1,0,ii), dRadR(1:mn,0,1,ii), dRadR(1:mn,1,1,ii), ifail )

    dRadR(1:mn,0,0,ii) = dRadR(1:mn,0,0,ii) * ajk(1:mn) ! poloidal integration; 03 Nov 16;
    dRadR(1:mn,1,0,ii) = dRadR(1:mn,1,0,ii) * ajk(1:mn)
    dRadR(1:mn,0,1,ii) = dRadR(1:mn,0,1,ii) * ajk(1:mn)
    dRadR(1:mn,1,1,ii) = dRadR(1:mn,1,1,ii) * ajk(1:mn)

    call invfft( mn, im(1:mn), in(1:mn), dRadR(1:mn,0,0,ii), dRadR(1:mn,1,0,ii), dRadR(1:mn,0,1,ii), dRadR(1:mn,1,1,ii), &
                 Nt, Nz, dRodR(1:Ntz,0,ii), dRodR(1:Ntz,1,ii) ) ! R, Z; 03 Nov 16;

    dRodR(1:Ntz,0,ii) = dRodR(1:Ntz,0,ii) / jiimag(1:Ntz) ! divide by length; 03 Nov 16;
    dRodR(1:Ntz,1,ii) = dRodR(1:Ntz,1,ii) / jiimag(1:Ntz)


    dRodZ(1:Ntz,0,ii) =                                - kjreal(1:Ntz) * jkimag(1:Ntz) * im(ii) * sini(1:Ntz,ii) / jireal(1:Ntz) ! dRodZjc;
    dRodZ(1:Ntz,1,ii) =                                + kjreal(1:Ntz) * jkimag(1:Ntz) * im(ii) * cosi(1:Ntz,ii) / jireal(1:Ntz) ! dRodZjs;

    ifail = 0
    call tfft( Nt, Nz, dRodZ(1:Ntz,0,ii), dRodZ(1:Ntz,1,ii), &
              mn, im(1:mn), in(1:mn), dRadZ(1:mn,0,0,ii), dRadZ(1:mn,1,0,ii), dRadZ(1:mn,0,1,ii), dRadZ(1:mn,1,1,ii), ifail )

    dRadZ(1:mn,0,0,ii) = dRadZ(1:mn,0,0,ii) * ajk(1:mn) ! poloidal integration; 03 Nov 16;
    dRadZ(1:mn,1,0,ii) = dRadZ(1:mn,1,0,ii) * ajk(1:mn)
    dRadZ(1:mn,0,1,ii) = dRadZ(1:mn,0,1,ii) * ajk(1:mn)
    dRadZ(1:mn,1,1,ii) = dRadZ(1:mn,1,1,ii) * ajk(1:mn)

    call invfft( mn, im(1:mn), in(1:mn), dRadZ(1:mn,0,0,ii), dRadZ(1:mn,1,0,ii), dRadZ(1:mn,0,1,ii), dRadZ(1:mn,1,1,ii), &
                 Nt, Nz, dRodZ(1:Ntz,0,ii), dRodZ(1:Ntz,1,ii) ) ! R, Z; 03 Nov 16;

    dRodZ(1:Ntz,0,ii) = dRodZ(1:Ntz,0,ii) / jiimag(1:Ntz) ! divide by length; 03 Nov 16;
    dRodZ(1:Ntz,1,ii) = dRodZ(1:Ntz,1,ii) / jiimag(1:Ntz)



    dZodR(1:Ntz,0,ii) =                                - kjimag(1:Ntz) * jkreal(1:Ntz) * im(ii) * sini(1:Ntz,ii) / jireal(1:Ntz) ! dZodRjc;
    dZodR(1:Ntz,1,ii) =                                + kjimag(1:Ntz) * jkreal(1:Ntz) * im(ii) * cosi(1:Ntz,ii) / jireal(1:Ntz) ! dZodRjs;

    ifail = 0
    call tfft( Nt, Nz, dZodR(1:Ntz,0,ii), dZodR(1:Ntz,1,ii), &
              mn, im(1:mn), in(1:mn), dZadR(1:mn,0,0,ii), dZadR(1:mn,1,0,ii), dZadR(1:mn,0,1,ii), dZadR(1:mn,1,1,ii), ifail )

    dZadR(1:mn,0,0,ii) = dZadR(1:mn,0,0,ii) * ajk(1:mn) ! poloidal integration; 03 Nov 16;
    dZadR(1:mn,1,0,ii) = dZadR(1:mn,1,0,ii) * ajk(1:mn)
    dZadR(1:mn,0,1,ii) = dZadR(1:mn,0,1,ii) * ajk(1:mn)
    dZadR(1:mn,1,1,ii) = dZadR(1:mn,1,1,ii) * ajk(1:mn)

    call invfft( mn, im(1:mn), in(1:mn), dZadR(1:mn,0,0,ii), dZadR(1:mn,1,0,ii), dZadR(1:mn,0,1,ii), dZadR(1:mn,1,1,ii), &
                 Nt, Nz, dZodR(1:Ntz,0,ii), dZodR(1:Ntz,1,ii) ) ! R, Z; 03 Nov 16;

    dZodR(1:Ntz,0,ii) = dZodR(1:Ntz,0,ii) / jiimag(1:Ntz) ! divide by length; 03 Nov 16;
    dZodR(1:Ntz,1,ii) = dZodR(1:Ntz,1,ii) / jiimag(1:Ntz)


    dZodZ(1:Ntz,0,ii) = cosi(1:Ntz,ii) * jireal(1:Ntz) - kjimag(1:Ntz) * jkimag(1:Ntz) * im(ii) * sini(1:Ntz,ii) / jireal(1:Ntz) ! dZodZjc;
    dZodZ(1:Ntz,1,ii) = sini(1:Ntz,ii) * jireal(1:Ntz) + kjimag(1:Ntz) * jkimag(1:Ntz) * im(ii) * cosi(1:Ntz,ii) / jireal(1:Ntz) ! dZodZjs;

    ifail = 0
    call tfft( Nt, Nz, dZodZ(1:Ntz,0,ii), dZodZ(1:Ntz,1,ii), &
              mn, im(1:mn), in(1:mn), dZadZ(1:mn,0,0,ii), dZadZ(1:mn,1,0,ii), dZadZ(1:mn,0,1,ii), dZadZ(1:mn,1,1,ii), ifail )

    dZadZ(1:mn,0,0,ii) = dZadZ(1:mn,0,0,ii) * ajk(1:mn) ! poloidal integration; 03 Nov 16;
    dZadZ(1:mn,1,0,ii) = dZadZ(1:mn,1,0,ii) * ajk(1:mn)
    dZadZ(1:mn,0,1,ii) = dZadZ(1:mn,0,1,ii) * ajk(1:mn)
    dZadZ(1:mn,1,1,ii) = dZadZ(1:mn,1,1,ii) * ajk(1:mn)

    call invfft( mn, im(1:mn), in(1:mn), dZadZ(1:mn,0,0,ii), dZadZ(1:mn,1,0,ii), dZadZ(1:mn,0,1,ii), dZadZ(1:mn,1,1,ii), &
                 Nt, Nz, dZodZ(1:Ntz,0,ii), dZodZ(1:Ntz,1,ii) ) ! R, Z; 03 Nov 16;

    dZodZ(1:Ntz,0,ii) = dZodZ(1:Ntz,0,ii) / jiimag(1:Ntz) ! divide by length; 03 Nov 16;
    dZodZ(1:Ntz,1,ii) = dZodZ(1:Ntz,1,ii) / jiimag(1:Ntz)


   enddo ! end of do ii; 03 Nov 16;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  end select ! end of select case( Igeometry ) ; 08 Feb 16;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(rzaxis)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine rzaxis

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
