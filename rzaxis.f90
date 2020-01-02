!> \defgroup grp_coord_axis Coordinate axis

!> \file rzaxis.f90
!> \brief The coordinate axis is assigned via a poloidal average over an arbitrary surface.

!> \brief The coordinate axis is assigned via a poloidal average over an arbitrary surface.
!> \ingroup grp_coord_axis
!> 
!> Specifies position of coordinate axis; \f${\bf x}_a(\zeta) \equiv \int {\bf x}_1(\theta,\zeta) dl \, / \int dl\f$.
!> 
!> **coordinate axis**
!>
!> <ul>
!> <li> The coordinate axis is _not_ an independent degree-of-freedom of the geometry.
!>       It is constructed by extrapolating the geometry of a given interface, as determined by \f$i \equiv\,\f$\c ivol which is given on input,
!>       down to a line.
!> <li> If the coordinate axis depends only on the _geometry_ of the interface and not the angle parameterization,
!>       then the block tri-diagonal structure of the the force-derivative matrix is preserved.
!> <li> Define the arc-length-weighted averages,
!>       \f{eqnarray}{ R_0(\zeta) \equiv \frac{\displaystyle \int_{0}^{2\pi} R_i(\theta,\zeta) \, dl}{\displaystyle \int_{0}^{2\pi} \!\!\!\! dl}, \qquad
!>                     Z_0(\zeta) \equiv \frac{\displaystyle \int_{0}^{2\pi} Z_i(\theta,\zeta) \, dl}{\displaystyle \int_{0}^{2\pi} \!\!\!\! dl},
!>       \f}
!>       where \f$dl \equiv \dot l \, d\theta = \sqrt{ \partial_\theta R_i(\theta,\zeta)^2 + \partial_\theta Z_i(\theta,\zeta)^2 } \, d\theta\f$.
!> <li> (Note that if \f$\dot l\f$ does not depend on \f$\theta\f$, i.e. if \f$\theta\f$ is the equal arc-length angle, then the expressions simplify.
!>        This constraint is not enforced.)
!> <li> The geometry of the coordinate axis thus constructed only depends on the geometry of the interface, i.e. 
!>       the angular parameterization of the interface is irrelevant.
!> </ul>
!>
!> **coordinate axis: derivatives**
!>
!> <ul>
!> <li> The derivatives of the coordinate axis with respect to the Fourier harmonics of the given interface are given by
!>       \f{eqnarray}{
!>       \displaystyle \frac{\partial R_0}{\partial R_{i,j}^c} & = & \displaystyle \int \left( \cos\alpha_j \; \dot l
!>                                                     -       \Delta R_i R_{i,\theta} \, m_j \sin\alpha_j / \; \dot l \right) d\theta / L \\
!>       \displaystyle \frac{\partial R_0}{\partial R_{i,j}^s} & = & \displaystyle \int \left( \sin\alpha_j \; \dot l
!>                                                     +       \Delta R_i R_{i,\theta} \, m_j \cos\alpha_j / \; \dot l \right) d\theta / L \\
!>       \displaystyle \frac{\partial R_0}{\partial Z_{i,j}^c} & = & \displaystyle \int \left( \;\;\;\;\;\;\;\;\;\;\;\;\,
!>                                                     -       \Delta R_i Z_{i,\theta} \, m_j \sin\alpha_j / \; \dot l \right) d\theta / L \\
!>       \displaystyle \frac{\partial R_0}{\partial Z_{i,j}^s} & = & \displaystyle \int \left( \;\;\;\;\;\;\;\;\;\;\;\;\,                             
!>                                                     +       \Delta R_i Z_{i,\theta} \, m_j \cos\alpha_j / \; \dot l \right) d\theta / L \\ \nonumber \\
!>       \displaystyle \frac{\partial Z_0}{\partial R_{i,j}^c} & = & \displaystyle \int \left( \;\;\;\;\;\;\;\;\;\;\;\;\,                             
!>                                                     -       \Delta Z_i R_{i,\theta} \, m_j \sin\alpha_j / \; \dot l \right) d\theta / L \\
!>       \displaystyle \frac{\partial Z_0}{\partial R_{i,j}^s} & = & \displaystyle \int \left( \;\;\;\;\;\;\;\;\;\;\;\;                            
!>                                                     +       \Delta Z_i R_{i,\theta} \, m_j \cos\alpha_j / \; \dot l \right) d\theta / L \\
!>       \displaystyle \frac{\partial Z_0}{\partial Z_{i,j}^c} & = & \displaystyle \int \left( \cos\alpha_j \; \dot l
!>                                                     -       \Delta Z_i Z_{i,\theta} \, m_j \sin\alpha_j / \; \dot l \right) d\theta / L \\
!>       \displaystyle \frac{\partial Z_0}{\partial Z_{i,j}^s} & = & \displaystyle \int \left( \sin\alpha_j \; \dot l
!>                                                     +       \Delta Z_i Z_{i,\theta} \, m_j \cos\alpha_j / \; \dot l \right) d\theta / L
!>       \f}
!>       where \f$\displaystyle L(\zeta) \equiv \int_{0}^{2\pi} \!\!\!\! dl\f$.
!> </ul>
!>
!> **some numerical comments**
!>
!> <ul>
!> <li> First, the differential poloidal length, \f$\dot l \equiv \sqrt{ R_\theta^2 + Z_\theta^2 }\f$, is computed in real space using 
!>       an inverse FFT from the Fourier harmonics of \f$R\f$ and \f$Z\f$.
!> <li> Second, the Fourier harmonics of \f$dl\f$ are computed using an FFT.
!>       The integration over \f$\theta\f$ to construct \f$L\equiv \int dl\f$ is now trivial: just multiply the \f$m=0\f$ harmonics of \f$dl\f$ by \f$2\pi\f$.
!>       The \c ajk(1:mn) variable is used, and this is assigned in readin() .
!> <li> Next, the weighted \f$R \, dl\f$ and \f$Z \, dl\f$ are computed in real space, and the poloidal integral is similarly taken.
!> <li> Last, the Fourier harmonics are constructed using an FFT after dividing in real space.
!> </ul>
!>
!> @param[in]  Mvol
!> @param[in]  mn
!> @param      iRbc
!> @param      iZbs
!> @param      iRbs
!> @param      iZbc
!> @param[in]  ivol
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
