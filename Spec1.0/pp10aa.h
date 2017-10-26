!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Constructs \Poincare plot in external domain (free-boundary).

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine pp10aa
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half, one, pi2, goldenmean

  use numerical, only : small

  use fileunits, only : ounit, lunit, eunit

  use inputlist, only : Wpp10aa, Wbf10aa, ext, Nvol, Ni, nPtrj, nPpts, odetol

  use cputiming, only : Tpp10aa

  use allglobal, only : myid, ncpu, cpus, mn, im, in, Nz, Rmin, Rmax, Zmin, Zmax, Ntz, Rij, Zij, pi2nfp, Ltangent, interfacelabel

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS

  LOGICAL             :: exist
  INTEGER             :: lvol, Lcurv, lnPtrj, ii, jj, kk, Node, id02bjf, Lstart, iioff
  REAL                :: lss, phistart, phiend, lrzp(1:7), rzps(1:3), rzpe(1:3), rzpo(1:3), workarray(1:20*7), tol
  REAL, allocatable   :: poincaredata(:,:,:) ! this perhaps should be changed to be consistent with pp00aa, but this routine will be redundant at some point anyway . . .
  CHARACTER           :: RA, suf*4

  external            :: bf10aa, bf10aa_out, bf10aa_end

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  BEGIN(pp10aa)

  FATALMESS(pp10aa,.true.,definition of Rij Zij has changed)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!  if( myid.eq.0 ) then ! for comparison with extender, diagno, . . .
!   
!   inquire(file="rphiz.in",exist=exist)
!   if( exist ) then 
!    cput = GETTIME ; write(ounit,'("pp10aa : ",f10.2," : myid=",i3," ; computing plasma field at [R,p,Z] points given in rphiz.in ;")')cput-cpus,myid
!    open(eunit+myid,file="rphiz.out."//trim(ext)//".spec",status="unknown") ! this file will be written to in bf10aa;
!    open(lunit+myid,file="rphiz.in",status="old")                           ! arbitrary selection of points at which field will be computed (to be compared to e.g. extender);
!    do ! will read file until error (end of file) is encountered;
!     read(lunit+myid,*,iostat=ios)lrzp(1:3) ; lrzp(1:3) = (/ lrzp(1), lrzp(3), lrzp(2) /)
!     if( ios.ne.0 ) exit
!     arblabel = zero ; call bf10aa( arblabel, lrzp(1:3), Brzp(1:3) )
!    enddo
!    close(lunit+myid)
!    close(eunit+myid)
!    cput = GETTIME ; write(ounit,'("pp10aa : ",f10.2," : myid=",i3," ; computed  plasma field at [R,p,Z] points given in rphiz.in ; data saved in rphiz.ext ;")')cput-cpus,myid
!   endif
!   
!  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  Ltangent = 0 ; Node = 7 ; RA = 'D' ; tol = odetol ! construct outer boundary; do not require metrics; do not require tangential field;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( nPtrj(Nvol+1).gt.0 ) then ; lnPtrj = nPtrj(Nvol+1) ! set number of trajectories to be followed; required number of trajectories explicitly given;
  else                          ; lnPtrj =    Ni(Nvol  ) ! set number of trajectories to be followed; defaults to radial resolution in outermost volume;
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  Lstart = 0 ! used to select starting locations . . . .

  select case( Lstart )

  case( 0 )

   lvol = Nvol ; lss = interfacelabel(lvol       ) ; Lcurv = 0 ; WCALL(pp10aa,co01aa,( lvol, lss, Lcurv, Ntz, mn )) ! construct geometry of interface in R,Z;

   rzps(1:3) = (/       Rij(0,1)  , zero, zero /) ! from edge   of plasma;

  case( 1 )

   rzps(1:3) = (/ half*(Rmin+Rmax), zero, zero /) ! from center of mgrid ;

  case default

   FATALMESS( pp10aa, .true., invalid Lstart )

  end select
  
  rzpe(1:3) = (/ Rmax , zero, zero /) ! to mgrid boundary;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  write(suf,'(i4.4)') myid

  open( lunit+myid, file="."//trim(ext)//".external."//suf, status="unknown", form="unformatted" )

  write(lunit+myid) Rmin, Rmax, Zmin, Zmax ! this line is additional to the other Poincare data files;
  
  call flush( lunit+myid )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
  RALLOCATE(poincaredata,(1:4,0:Nz-1,1:nPpts)) ! for block writing to file (allows faster reading for post-processing plotting routines);
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  iioff = 2 ! don't start too close to plasma edge; virtual casing calculation becomes too slow;
  
  FATALMESS(pp10aa,iioff.le.0,cannot start inside virtual casing boundary)

  do ii = iioff, lnPtrj-1
   
   if( myid.ne.modulo(ii-iioff,ncpu) ) cycle ! parallel;
   
   rzpo(1:3) = rzps(1:3) + ii * ( rzpe(1:3) - rzps(1:3) ) / lnPtrj ! original starting location;
  
   lrzp(1:7) = (/ rzpo(1), rzpo(2), rzpo(3), one, zero, zero, one /)
   
   poincaredata(1:4,0:Nz-1,1:nPpts) = zero ! reset output data at start of following each field line;
   
   if( Wpp10aa ) then ; cput = GETTIME ; write(ounit,1000) cput-cpus, myid, ii, rzpo(1:3)
   endif

   do jj = 1, nPpts ! loop over iterations;
    
    phistart = zero
    
    do kk = 0, Nz-1 ! loop over toroidal Poincare cross sections;
     
     poincaredata(1:4,kk,jj) = (/ zero, zero, lrzp(1), lrzp(2) /)         ! save Poincare information to array;
     
     phiend = phistart + pi2nfp / Nz
     
     id02bjf = 1
     CALL(D02BJF,( phistart, phiend, Node, lrzp(1:7), bf10aa, tol, RA, bf10aa_out, bf10aa_end, workarray, id02bjf ))
     
     cput = GETTIME
     select case( id02bjf ) !                                                          123456789012345678901
     case(   0 ) ;!if( Wpp10aa ) write(ounit,1000)cput-cpus,myid,ii,rzpo(1:3),id02bjf,"success              "
     case(   1 ) ;               write(ounit,1000)cput-cpus,myid,ii,rzpo(1:3),id02bjf,"input error          "
     case(   2 ) ;!if( Wpp10aa ) write(ounit,1000)cput-cpus,myid,ii,rzpo(1:3),id02bjf,"no further progress  "
     case(   3 ) ;               write(ounit,1000)cput-cpus,myid,ii,rzpo(1:3),id02bjf,"tol too small        "
     case(   4 ) ;               write(ounit,1000)cput-cpus,myid,ii,rzpo(1:3),id02bjf,"xsol not reset       "
     case(   5 ) ;               write(ounit,1000)cput-cpus,myid,ii,rzpo(1:3),id02bjf,"invalid xsol         "
     case(   6 ) ;!if( Wpp10aa ) write(ounit,1000)cput-cpus,myid,ii,rzpo(1:3),id02bjf,"g did not change sign"
     case(   7 ) ;               write(ounit,1000)cput-cpus,myid,ii,rzpo(1:3),id02bjf,"serious error        "
     case default
      FATALMESS(pp10aa,.true.,invalid ifail returned from D02BJF)
     end select
     
     if( id02bjf.ne.0 .and. id02bjf.ne.6 ) exit ! an error has occured integrating this field line; skip to next field line;
     
    enddo ! end of do kk = 0, Nz-1 ! loop over toroidal Poincare cross sections;
    
    if( id02bjf.ne.0 .and. id02bjf.ne.6 ) exit ! an error has occured integrating this field line; skip to next field line;
    
   enddo ! end of do jj = 1, nPpts ! loop over iterations;
   
   if( Wpp10aa ) then ; cput = GETTIME ; write(ounit,1000)cput-cpus,myid,ii,rzpo(1:3),id02bjf,"success              "
   endif
   
   if( id02bjf.eq.0 ) then ; write(lunit+myid) Nz, nPpts                        ! write field line trajectory to file;
    ;                      ; write(lunit+myid) poincaredata(1:4,0:Nz-1,1:nPpts) ! write field line trajectory to file;
    ;                      ; call flush(lunit+myid)                             ! after debugging, can leave this line out;
   endif
   
  enddo ! end of do ii = iioff, lnPtrj-1;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  close( lunit+myid )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  DEALLOCATE(poincaredata)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(pp10aa)
  
1000 format("pp10aa : ",f10.2," : myid=",i3," : ii="i4" ; [R,Z,p] = ["es13.5" ,"es13.5" ,"es13.5" ] ;":" id02bjf=",i3," ; "a21" ;")
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine pp10aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
