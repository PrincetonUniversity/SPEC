!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Returns the magnetic field field line equations in external domain (free-boundary).

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine bf10aa( arblabel, rzp, Brzp )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  use constants, only : zero, one

  use numerical, only : small

  use fileunits, only : ounit, eunit

  use inputlist, only : Wbf10aa

  use cputiming, only : Tbf10aa

  use allglobal, only : myid, cpus, Lmgridexist, Ltangent
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  REAL, intent(in)  :: arblabel, rzp(1:7)
  REAL, intent(out) :: Brzp(1:7)
  
! LOGICAL           :: opened
  LOGICAL           :: Lcomputevirtualfield, Lincludevirtualfield, Lcomputevacuumfield, Lincludevacuumfield
  INTEGER           :: ncalls = 0, ifail, imgridfail, ivirtualfail
  REAL              :: xyz(1:3), Bxyz(1:3), dBxyzdxyz(1:3,1:3), vacuumBrzp(1:3), dBrzpdrzp(1:3,1:3), virtualBrzp(1:3), cosphi, sinphi, cpup, cpuc, TM(1:2,1:2)
  
  BEGIN(bf10aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  Lcomputevirtualfield = .true.  ; Lincludevirtualfield = .true.  ! local logical variables;
  Lcomputevacuumfield  = .true.  ; Lincludevacuumfield  = .true.  ! local logical variables;
  
  if( Lincludevirtualfield ) Lcomputevirtualfield = .true.        ! local logical variables;
  if( Lincludevacuumfield )  Lcomputevacuumfield = .true.         ! local logical variables;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  FATALMESS(bf10aa,rzp(1).lt.small,divide by zero in coordinate transformation)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  cpuc = zero ! records time taken to compute (i.e. interpolate mgrid) magnetic field produced by external coils;
  
  vacuumBrzp(1:3) = zero ; dBrzpdrzp(1:3,1:3) = zero

  if( Lcomputevacuumfield .and. Lmgridexist ) then ! include field due to external coil currents, as provided by mgrid and interpolated using ezsplines; see mg00aa;

   cpuo = GETTIME
   call mgridfield( rzp(1:3), vacuumBrzp(1:3), dBrzpdrzp(1:3,1:3), imgridfail ) ! ESSENTIAL THAT MGRIDFIELD DOES NOT CHANGE RZP;
   cput = GETTIME ; cpuc = cput-cpuo

   select case( imgridfail ) 
   case( 0 ) ;                                                                       ! mgrid field has been successfully constructed;
   case( 1 ) ; Brzp(1:7) = (/ zero, zero, one, zero, zero, zero, zero /) ; goto 9999 ! the point [R,Z,p] lies outside the mgrid domain;
   case default
    FATALMESS(bf10aa,.true.,illegal value of imgridfail returned from mgridfield)    ! perhaps need to update mg00aa and/or bf10aa;
   end select
   
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \item The transformation from cylindrical $(R,\phi,Z)$ to Cartesian $(x,y,z)$ is described by $x = R \cos \phi$, $y = R \sin\phi$ and $z=z$.
  
  cpup = zero ! records time taken to compute field from virtual casing given tangential field on boundary;
  
  if( Lcomputevirtualfield ) then
   
   cpuo = GETTIME
   
   cosphi = cos(rzp(3)) ; sinphi = sin(rzp(3))
   
   xyz(1:3) = (/ rzp(1)*cosphi, rzp(1)*sinphi, rzp(2) /) ! cylindrical to Cartesian coordinate transformation;
   
   call vc00aa( xyz(1:3), Bxyz(1:3), dBxyzdxyz(1:3,1:3), ivirtualfail ) ! virtual casing: given tangential field on surface, compute field;
   
   FATALMESS(bf10aa,ivirtualfail.ne.0,an error has occurred in vc00aa)

!latex \item The vector transformation from Cartesian to cylindrical (as required by the field line integration routine \verb+pp10aa+) is given:
!latex       \be B^x \nabla x + B^y \nabla y + B^z \nabla z = 
!latex           B^x (\cos \phi \nabla R - \sin \phi R \nabla \phi) + B^y (\sin \phi \nabla R + \cos \phi R \nabla \phi) + B^z \nabla z.
!latex       \ee
   
   virtualBrzp(1) = (   cosphi * Bxyz(1) + sinphi * Bxyz(2) )          ! Cartesian to cylindrical vector transformation;
   virtualBrzp(2) =              Bxyz(3)                               ! Cartesian to cylindrical vector transformation;
   virtualBrzp(3) = ( - sinphi * Bxyz(1) + cosphi * Bxyz(2) ) / rzp(1) ! Cartesian to cylindrical vector transformation;
   
! this is included for comparision with extender; only required if rphiz.in pointsfile is present (see pp10aa);
!   inquire( unit=eunit+myid, opened=opened )
!   if( opened ) then
!    write(eunit+myid,'(i8,7es16.8,i2)')&
!     ncalls,rzp(1),rzp(3),rzp(2),virtualBrzp(1),virtualBrzp(3)*rzp(1),virtualBrzp(2),sqrt(virtualBrzp(1)**2+virtualBrzp(3)**2+virtualBrzp(2)**2),0
!    ncalls=ncalls+1
!   endif

   cput = GETTIME ; cpup = cput-cpuo
   
  endif ! end of if( Lcomputevirtualfield );
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! residual(1:3) = vacuumBrzp(1:3)/vacuumBrzp(3) - virtualBrzp(1:3)/virtualBrzp(3) ! this is the difference between the vacuum field and the virtual casing field;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  ;                          Brzp(1:3) = zero                         ! initialize summation ;
  if( Lincludevacuumfield  ) Brzp(1:3) = Brzp(1:3) +  vacuumBrzp(1:3) ! vacuum  field        ;
  if( Lincludevirtualfield ) Brzp(1:3) = Brzp(1:3) + virtualBrzp(1:3) ! virtual field        ;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! FATALMESS(bf10aa,abs(Brzp(3)).lt.small,zero toroidal field)
  
  if( abs(Brzp(3)).gt.small ) then ; Brzp(1:3) = Brzp(1:3) / Brzp(3)     ! normalize to toroidal field; 
  else                             ; Brzp(1:3) = (/ zero, zero, one /) ! provide dummy values;
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  select case( Ltangent ) ! determines whether derivatives are required;
   
  case( 0 )
   
   Brzp(4:7) = zero ! provide dummy values that should be simple to integrate;
   
  case( 1 )
   
   FATALMESS(bf10aa,.true.,need to implement tangent field)
   
   TM(1:2,1:2) = zero
   
!   TM(1,1) = ( dBu(1,1) - Bst(1) * dBu(3,1) ) / dBu(3,0) !  d \dot R / d R ;
!   TM(1,2) = ( dBu(1,2) - Bst(1) * dBu(3,2) ) / dBu(3,0) !  d \dot R / d Z ;
!   TM(2,1) = ( dBu(2,1) - Bst(2) * dBu(3,1) ) / dBu(3,0) !  d \dot Z / d R ;
!   TM(2,2) = ( dBu(2,2) - Bst(2) * dBu(3,2) ) / dBu(3,0) !  d \dot Z / d Z ;
   
   Brzp(3) = TM(1,1) * rzp(4) + TM(1,2) * rzp(6) ! tangent map;
   Brzp(4) = TM(1,1) * rzp(5) + TM(1,2) * rzp(7) ! tangent map;
   Brzp(5) = TM(2,1) * rzp(4) + TM(2,2) * rzp(6) ! tangent map;
   Brzp(6) = TM(2,1) * rzp(5) + TM(2,2) * rzp(7) ! tangent map;
   
   Brzp(4:7) = zero
   
  case default
   
   FATALMESS(bf10aa,.true.,invalid Ltangent)
   
  end select
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Wbf10aa ) then ; cput = GETTIME ; write(ounit,1001)cput-cpus,rzp(1:3),cpuc,cpup,vacuumBrzp(1:2)/vacuumBrzp(3),virtualBrzp(1:2)/virtualBrzp(3)
  endif
  
1001 format("bf10aa : ",f10.2," : [R,Z,phi] = ["3es13.05" ] ; time="2f8.4" ; Bc= ["es23.15" ,"es23.15" ] & Bp= ["es23.15" ,"es23.15" ] ;":,3es23.15)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(bf10aa)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine bf10aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine bf10aa_out( arblabel, rzp )

  use allglobal, only : pi2nfp

  LOCALS

  REAL, intent(inout) :: arblabel
  REAL, intent(in)    :: rzp(1:3)
  
  arblabel = arblabel + pi2nfp ! next intermediate output location;

  return

end subroutine bf10aa_out

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

REAL function bf10aa_end( arblabel, rzp )

  use allglobal, only : pi2nfp

  LOCALS

  REAL, intent(inout) :: arblabel
  REAL, intent(in)    :: rzp(1:3)
  
  bf10aa_end = arblabel - pi2nfp ! integration termination;

  return

end FUNCTION bf10aa_end

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
