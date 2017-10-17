!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Computes normal field, and loop integrals to determine enclosed currents, on computational boundary.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \end{enumerate} \subsubsection{free-boundary constraint} \begin{enumerate}
    
!latex \item The normal field at the computational boundary should be equal to $\left({\bf B}_P + {\bf B}_C\right)\cdot {\bf n}$,
!latex       where ${\bf B}_P$ is the plasma field and is computed using virtual casing, 
!latex       and ${\bf B}_C$ is the magnetic field produced by the coils and is computed using Biot-Savart
!latex       (using the \verb+mgrid+--\verb+ezspline+ interpolation).

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \end{enumerate} \subsubsection{construction of normal field} \begin{enumerate}

!latex \item The normal field depends on geometry:

!latex \begin{enumerate}

!latex \item \verb+Igeometry.eq.1+ : Cartesian ;
!latex \item \verb+Igeometry.eq.2+ : Cylindrical ;
!latex \item \verb+Igeometry.eq.3+ : Toroidal ; 
!latex ${\bf e}_\t \times {\bf e}_\z = - R \, Z_\theta \, \hat r + (Z_\theta \,R_\zeta - R_\theta \,Z_\zeta) \hat \phi + R \,R_\theta \,\hat z$.
!latex
!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine bn00aa( mn, Ntz )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half, one, pi2
  
  use numerical, only : small
  
  use fileunits, only : ounit
  
  use inputlist, only : Wmacros, Wbn00aa, Igeometry, Nvol, Ntor, curtor, curpol, &
                        norblend, normalerr, maxfbits
  
  use cputiming, only : Tbn00aa
  
  use allglobal, only : ncpu, myid, cpus, pi2pi2nfp, pi2nfp, Mvol, &
                        Nt, Nz, &
                        Rij, Zij, sg, guvij, &
                        Lcoordinatesingularity, &
                        im, in, &
                        trigwk, trigm, trign, isr, Nt, Nz, efmn, ofmn, cfmn, sfmn, &
                        ijreal, ijimag, jireal, jiimag ,&
                        iBns, iBnc, iCns, iCnc, iPns, iPnc, &
                        Ltangent, Bnserror,  &
                        nFreeIterations, Lcontinueiterations, NOTstellsym, &
                        Rmin, Zmin, Rmax, Zmax
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in) :: mn, Ntz
  
  INTEGER             :: lvol, Lcurvature, Lparallel, jj, kk, jk, kkmodnp, jkmodnp, imgridfail, ivirtualfail, ifail
  
  REAL                :: lss, zeta, teta, coszeta, sinzeta
  
  REAL                :: DRmin, DRmax, DZmin, DZmax ! 11 Oct 12; computational domain; shall compare to mgrid domain;
   
  REAL                :: rzp(1:3) ! cylindrical coordinates; 12 Oct 12;
  REAL                :: xyz(1:3) ! Cartesian   coordinates; 12 Oct 12;
  
  REAL                :: Brzp(1:3), Bxyz(1:3), dBxyzdxyz(1:3,1:3), dBrzpdrzp(1:3,1:3)
  
  REAL                :: vacBrzp(1:3) ! vacuum field as computed by interpolation of mgrid; 12 Oct 12;
  REAL                :: virBrzp(1:3) ! plasma field as computed by virtual casing        ; 12 Oct 12;
  
  BEGIN(bn00aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  FATALMESS(bn00aa, Igeometry.ne.3, need to update geometry )
#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  Lparallel = 1 ! controls choice of parallelization; see below;
  
  Ltangent = 0 ! used in virtual field calculation; the derivatives of the field wrt position are not computed if Ltangent=0;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  lvol = Mvol ; lss = one ; Lcurvature = 1 ; Lcoordinatesingularity = .false. ! will only require normal field on outer interface = computational boundary; 
  
  WCALL(bn00aa, co01aa,( lvol, lss, Lcurvature, Ntz, mn )) ! will need Rij, Zij;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( myid.eq.0 .and. nFreeIterations.eq.0 ) then ! 11 Oct 12; just for convenience; to show computational boundary; ! this can be deleted after debugging;
   
   DRmin = minval(Rij(:,0,0)) ; DRmax = maxval(Rij(:,0,0)) ! 11 Oct 12; 
   DZmin = minval(Zij(:,0,0)) ; DZmax = maxval(Zij(:,0,0)) ! 11 Oct 12; 
   
   cput = GETTIME
   write(ounit,'("bn00aa : ", 10x ," : ")') ! 11 Oct 12; 
   write(ounit,'("bn00aa : ", 10x ," : myid=",i3," ; [ Rmin, Rmax]=["es13.5" ,"es13.5" ] ; [ Zmin, Zmax]=["es13.5" ,"es13.5" ] ;")') &
myid,  Rmin,  Rmax,  Zmin,  Zmax ! 23 Oct 12;
   write(ounit,'("bn00aa : ", 10x ," : myid=",i3," ; [DRmin,DRmax]=["es13.5" ,"es13.5" ] ; [DZmin,DZmax]=["es13.5" ,"es13.5" ] ;")') &
myid, DRmin, DRmax, DZmin, DZmax
   
   if( DRmin.lt.Rmin .or. DRmax.gt.Rmax .or. DZmin.lt.Zmin .or. DZmax.gt.Zmax ) then ! 11 Oct 12; 
    write(ounit,'("bn00aa : ", 10x ," : ")') ! 23 Oct 12;
    write(ounit,'("bn00aa : ", 10x ," : myid=",i3," ; computational boundary exceeds mgrid domain ; this will certainly cause an error below ;")') myid
   endif
   
  endif ! 11 Oct 12; end of if( myid.eq.0 .and. nFreeIterations.eq.0 ) then
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  ijreal(1:Ntz) = zero ! normal plasma field; 15 Oct 12;
  ijimag(1:Ntz) = zero ! normal coils  field; 15 Oct 12;

  jireal(1:Ntz) = zero ! vacuum  covariant poloidal field; 15 Oct 12;
  jiimag(1:Ntz) = zero ! virtual covariant poloidal field; 15 Oct 12;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  do kk = 0, Nz-1 ; zeta = kk * pi2nfp / Nz ! loop over real space grid; toroidal;
   
   do jj = 0, Nt-1 ; teta = jj * pi2    / Nt ! loop over real space grid; poloidal; note that this is sequential in jj;
    
    jk = 1 + jj + kk*Nt
    
    
    select case( Lparallel ) ! perform in parallel;
    case( 0 )
     if( myid.ne.modulo(kk,ncpu) ) cycle
    case( 1 ) 
     if( myid.ne.modulo(jk-1,ncpu) ) cycle ! 11 Oct 12; this is a weird parallelization, but perhaps better exploits all available cpus;
    case default
     FATALMESS(bn00aa, .true., invalid Lparallel in parallelization loop )
    end select
    
    
    rzp(1:3) = (/ Rij(jk,0,0), Zij(jk,0,0), zeta /) ! shorthand; 17 Apr 13;
    
#ifdef DEBUG
    FATALMESS(bn00aa, rzp(1).lt.small, divide by zero in coordinate transformation )
#endif
    
    FATALMESS(bn00aa, .true., need to revise mgridfield)
    vacBrzp(1:3) = zero ; imgridfail = 2

   !WCALL( bn00aa, mgridfield, ( rzp(1:3), vacBrzp(1:3), dBrzpdrzp(1:3,1:3), imgridfail )) ! returns coil field in cylindrical coordinates;
    
#ifdef DEBUG
    FATALMESS(bn00aa, imgridfail.eq.1, outside mgrid domain       ) ! 11 Oct 12; 
    FATALMESS(bn00aa, imgridfail.eq.2, mgrid has not been splined ) ! 11 Oct 12; 
#endif
    
!latex \end{enumerate} \subsubsection{coordinate and vector transformation} \begin{enumerate}

!latex \item The transformation from cylindrical to Cartesian is 
!latex       \be x = R \cos \z, \quad
!latex           y = R \sin \z, \quad
!latex           z = Z.
!latex       \ee
!latex \item This induces the vector transformation
!latex       \be B^R    = B_x \cos \z + B_y \sin \z, \quad
!latex           B^\phi = ( - B_x \sin \z + B_y \cos \z ) / R, \quad
!latex           B^Z    = B_z.
!latex       \ee

    coszeta = cos( zeta ) ; sinzeta = sin( zeta ) ! 10 Apr 13; required for coordinate transformation, and vector transformation;
    
    xyz(1:3) = (/ rzp(1)*coszeta, rzp(1)*sinzeta, rzp(2) /) ! cylindrical to Cartesian coordinate transformation;
    
!#ifdef DEBUG
!    if( Wbn00aa ) then
!     cput = GETTIME
!     write(ounit,'("bn00aa : ",f10.2," : myid=",i3," ; calling vc00aa : xyz="3es13.5" ;")') cput-cpus, myid, xyz(1:3)
!    endif
!#endif
    
! virtual casing: given tangential field compute field; returns field in Cartesian;
    WCALL( bn00aa, vc00aa, ( xyz(1:3), Bxyz(1:3), dBxyzdxyz(1:3,1:3), ivirtualfail )) 
    
#ifdef DEBUG
    FATALMESS(bn00aa, ivirtualfail.ne.0, an error has occurred in vc00aa )
#endif
    
    virBrzp(1) =     coszeta*Bxyz(1) + sinzeta*Bxyz(2)            ! 10 Oct 12; Cartesian to cylindrical vector;
    virBrzp(2) =             Bxyz(3)                              ! 10 Oct 12; Cartesian to cylindrical vector;
    virBrzp(3) = ( - sinzeta*Bxyz(1) + coszeta*Bxyz(2) ) / rzp(1) ! 10 Oct 12; Cartesian to cylindrical vector;
    
    ijreal(jk) = -   Zij(jk,2,0)                                       * virBrzp(1) * Rij(jk,0,0) &
                 + ( Zij(jk,2,0)*Rij(jk,3,0)-Rij(jk,2,0)*Zij(jk,3,0) ) * virBrzp(3) * Rij(jk,0,0) &
                 +   Rij(jk,2,0)                                       * virBrzp(2) * Rij(jk,0,0) !plasma-normal;

    ijimag(jk) = -   Zij(jk,2,0)                                       * vacBrzp(1) * Rij(jk,0,0) &
                 + ( Zij(jk,2,0)*Rij(jk,3,0)-Rij(jk,2,0)*Zij(jk,3,0) ) * vacBrzp(3) * Rij(jk,0,0) &
                 +   Rij(jk,2,0)                                       * vacBrzp(2) * Rij(jk,0,0) !coils -normal;

!latex \end{enumerate} \subsubsection{enclosed currents} \begin{enumerate}
    
!latex \item The enclosed currents are given by line integrals,
!latex       \be \int_{{\cal S}} {\bf j} \cdot d{\bf S}  =  \int_{\partial {\cal S}} {\bf B} \cdot d{\bf l}.
!latex       \ee
!latex       \be               {\bf e}_\t &=&     R_\t {\bf \hat r} +     Z_\t {\bf \hat z}                           ,\\
!latex                         {\bf e}_\z &=&     R_\z {\bf \hat r} +     Z_\z {\bf \hat z} +        R \, \hatboldphi ,\\
!latex           {\bf B} \cdot {\bf e}_\t &=& B^R R_\t              + B^Z Z_\t                                         \\
!latex           {\bf B} \cdot {\bf e}_\z &=& B^R R_\z              + B^Z Z_\z              + B^\phi R^2
!latex       \ee
    
    jireal(jk) = vacBrzp(1) * Rij(jk,2,0) + vacBrzp(2) * Zij(jk,2,0)                          ! B \cdot e_\t from coils ;
    jiimag(jk) = vacBrzp(1) * Rij(jk,3,0) + vacBrzp(2) * Zij(jk,3,0) + vacBrzp(3) * rzp(1)**2 ! B \cdot e_\z from coils ;
    
   enddo ! end of do jj;
   
  enddo ! end of do kk;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  do kk = 0, Nz-1
   
   kkmodnp = modulo(kk,ncpu)
   
   select case( Lparallel )
    
   case( 0 )
    
    RlBCAST(ijreal(1+kk*Nt:Nt+kk*Nt),Nt,kkmodnp) ! plasma; 03 Apr 13;
    RlBCAST(ijimag(1+kk*Nt:Nt+kk*Nt),Nt,kkmodnp) ! coils ; 03 Apr 13;

    RlBCAST(jireal(1+kk*Nt:Nt+kk*Nt),Nt,kkmodnp)
    RlBCAST(jiimag(1+kk*Nt:Nt+kk*Nt),Nt,kkmodnp)

   case( 1 )
    
    do jj = 0, Nt-1
     
     jk = 1 + jj + kk*Nt
     
     jkmodnp = modulo(jk-1,ncpu)
     
     RlBCAST(ijreal(jk),1,jkmodnp) ! plasma; 03 Apr 13;
     RlBCAST(ijimag(jk),1,jkmodnp) ! coils ; 03 Apr 13;
 
     RlBCAST(jireal(jk),1,jkmodnp)
     RlBCAST(jiimag(jk),1,jkmodnp)
     
    enddo
    
   case default 
    
    FATALMESS(bn00aa, .true., invalid Lparallel for broadcasting )
    
   end select
   
  enddo ! 11 Oct 12;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  FATALMESS(bn00aa, NOTstellsym, under construction )
#endif
  
  call tfft( Nt, Nz, ijreal(1:Ntz), ijimag(1:Ntz), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz), &
             mn, im(1:mn), in(1:mn), iPnc(1:mn), iPns(1:mn), iCnc(1:mn), iCns(1:mn), ifail ) ! Fourier decompose normal field;

  ofmn(2:mn) = + iPns(2:mn) + iCns(2:mn) * half / pi2 ! add plamsa field and coil field at computational boundary; 28 Apr 13;

  Bnserror = sum( abs( iBns(2:mn) - ofmn(2:mn) ) ) / (mn-1) ! quantify error in assumed normal field at computational boundary; ! 17 Oct 12;
  
  if( myid.eq.0 ) then ! 28 Apr 13;
   cput = GETTIME
   write(ounit,'("bn00aa : ", 10x ," : ")')
   write(ounit,'("bn00aa : ",f10.2," : ",i3," iBns=["99(es8.1","))') cput-cpus, nFreeIterations, iBns(2:mn) ! 18 Oct 12;
   write(ounit,'("bn00aa : ",f10.2," : ",i3," iPns=["99(es8.1","))') cput-cpus, nFreeIterations, iPns(2:mn)
   write(ounit,'("bn00aa : ",f10.2," : ",i3," iCns=["99(es8.1","))') cput-cpus, nFreeIterations, iCns(2:mn)
  !write(ounit,'("bn00aa : ",f10.2," : ",i3," iCns=["99(es8.1","))') cput-cpus, nFreeIterations, iBns(2:mn)-iPns(2:mn) ! the coil field;
   write(ounit,'("bn00aa : ", 10x ," : ")')
   write(ounit,'("bn00aa : ",f10.2," : ",i3," |dBns|="es10.2" ; norblend="f6.3" ;")') cput-cpus, nFreeIterations, Bnserror, norblend
  endif
  
  if( maxfbits.gt.0 ) iBns(2:mn) = norblend * iBns(2:mn) + ( one - norblend ) * ofmn(2:mn) ! Picard `update' of normal field at computational boundary;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( myid.eq.0 ) then ! construct loop integrals to identify enclosed coil currents; 03 Apr 13;
   
   lvol = Mvol ; lss = -one ; Lcurvature = 1 ; Lcoordinatesingularity = .false.
   
   WCALL(bn00aa, co01aa,( lvol, lss, Lcurvature, Ntz, mn )) ! will need Rij, Zij;
   
   call tfft( Nt, Nz, jireal(1:Ntz), jiimag(1:Ntz), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz), &
              mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), ifail ) ! Fourier decompose normal field;
   
! write(ounit,'("bn00aa : ", 10x ," : ")')
! write(ounit,'("bn00aa : ", 10x ," : B_t mn="99es13.5)') efmn(1:mn)
! write(ounit,'("bn00aa : ", 10x ," : B_t mn="99es13.5)') ofmn(1:mn)
! write(ounit,'("bn00aa : ", 10x ," : B_z mn="99es13.5)') cfmn(1:mn)
! write(ounit,'("bn00aa : ", 10x ," : B_z mn="99es13.5)') sfmn(1:mn)
   
   write(ounit,'("bn00aa : ", 10x ," : ")')
   write(ounit,'("bn00aa : ", 10x ," : G ="es23.15" ;":" curtor="es13.5" ;")') cfmn(1)

   curpol = cfmn(1) ! 17 Apr 13;
   
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RlBCAST(curpol,1,0)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Bnserror.gt.normalerr .and. nFreeIterations.lt.maxfbits ) then ; Lcontinueiterations = .true.
  else                                                               ; Lcontinueiterations = .false.
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(bn00aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine bn00aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
