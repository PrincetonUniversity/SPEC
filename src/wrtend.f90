!> \brief The restart file is written.
subroutine wrtend

    use constants, only :
  
    use numerical, only : machprec
  
    use fileunits, only : ounit, iunit
  
    use cputiming, only : Twrtend

    use allglobal
  
    use bndRep
  
    use inputlist
  
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
    LOCALS
  
    INTEGER              :: vvol !< iteration variable over all nested volumes
    INTEGER              :: imn  !< iteration variable for all Fourier harmonics
    INTEGER              :: ii   !< iteration variable for all Fourier harmonics
    INTEGER              :: jj   !< iteration variable
    INTEGER              :: kk   !< iteration variable
    INTEGER              :: jk   !< iteration variable
    INTEGER              :: Lcurvature !< curvature flag (?)
    INTEGER              :: mm   !< current poloidal mode number
    INTEGER              :: nn   !< current toroidal mode number
  
    REAL                 :: lss !< (?)
    REAL                 :: teta !< (?)
    REAL                 :: zeta !< (?)
    REAL                 :: st(1:Node) !< (?)
    REAL                 :: Bst(1:Node) !< (?)
    REAL                 :: BR !< (?)
    REAL                 :: BZ !< (?)
    REAL                 :: BP !< (?)
  
    BEGIN(wrtend)
  
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
    if( myid.ne.0 ) goto 9999
  
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
    if( Wwrtend ) then ; cput = GETTIME ; write(ounit,'("wrtend : ",f10.2," : myid=",i3," ; opening/writing ext.sp.end ;")') cput-cpus, myid
    endif
#endif
  
    open(iunit,file=trim(ext)//".sp.end",status="unknown") ! restart input file;
  
#ifdef DEBUG
    if( Wwrtend ) then ; cput = GETTIME ; write(ounit,'("wrtend : ",f10.2," : myid=",i3," ; writing physicslist ;")') cput-cpus, myid
    endif
#endif
 
    write(iunit,'("&physicslist")')
    write(iunit,'(" Igeometry   = ",i9        )') Igeometry
    write(iunit,'(" Istellsym   = ",i9        )') Istellsym
    write(iunit,'(" Lfreebound  = ",i9        )') Lfreebound
    write(iunit,'(" Lboundary   = ",i9        )') Lboundary
    write(iunit,'(" twoalpha    = ",i9        )') twoalpha
    write(iunit,'(" phiedge     = ",es23.15   )') phiedge
    write(iunit,'(" curtor      = ",es23.15   )') curtor
    write(iunit,'(" curpol      = ",es23.15   )') curpol
    write(iunit,'(" gamma       = ",es23.15   )') gamma
    write(iunit,'(" Nfp         = ",i9        )') Nfp
    write(iunit,'(" Nvol        = ",i9        )') Nvol
    write(iunit,'(" Mpol        = ",i9        )') Mpol
    write(iunit,'(" Ntor        = ",i9        )') Ntor
    write(iunit,'(" Lrad        = ",257i23    )') Lrad(1:Mvol)
    write(iunit,'(" tflux       = ",257es23.15)') tflux(1:Mvol)
    write(iunit,'(" pflux       = ",257es23.15)') pflux(1:Mvol)
    write(iunit,'(" helicity    = ",256es23.15)') helicity(1:Mvol)
    write(iunit,'(" pscale      = ",es23.15   )') pscale
    write(iunit,'(" Ladiabatic  = ",i9        )') Ladiabatic
    write(iunit,'(" pressure    = ",257es23.15)') pressure(1:Mvol)
    write(iunit,'(" adiabatic   = ",257es23.15)') adiabatic(1:Mvol)
    write(iunit,'(" mu          = ",257es23.15)') mu(1:Mvol)
    write(iunit,'(" Ivolume     = ",257es23.15)') Ivolume(1:Mvol)
    write(iunit,'(" Isurf       = ",257es23.15)') Isurf(1:Mvol-1)
    write(iunit,'(" Lconstraint = ",i9        )') Lconstraint
    write(iunit,'(" pl          = ",257i23    )') pl(0:Mvol)
    write(iunit,'(" ql          = ",257i23    )') ql(0:Mvol)
    write(iunit,'(" pr          = ",257i23    )') pr(0:Mvol)
    write(iunit,'(" qr          = ",257i23    )') qr(0:Mvol)
    write(iunit,'(" iota        = ",257es23.15)') iota(0:Mvol)
    write(iunit,'(" lp          = ",257i23    )') lp(0:Mvol)
    write(iunit,'(" lq          = ",257i23    )') lq(0:Mvol)
    write(iunit,'(" rp          = ",257i23    )') rp(0:Mvol)
    write(iunit,'(" rq          = ",257i23    )') rq(0:Mvol)
    write(iunit,'(" oita        = ",257es23.15)') oita(0:Mvol)
    write(iunit,'(" mupftol     = ",es23.15   )') mupftol
    write(iunit,'(" mupfits     = ",i9        )') mupfits
    write(iunit,'(" Lreflect    = ",i9        )') Lreflect
    write(iunit,'(" rpol        = ",es23.15   )') rpol
    write(iunit,'(" rtor        = ",es23.15   )') rtor
  
    if( Lfreebound.eq.1 .or. Zbs(0,1).gt.zero ) then
     do ii = 1, mn_field ; mm = im_field(ii) ; nn = in_field(ii) / Nfp 
        Rbc(nn,mm) = iRbc(ii,Nvol) ; Zbs(nn,mm) = iZbs(ii,Nvol)
        Rbs(nn,mm) = iRbs(ii,Nvol) ; Zbc(nn,mm) = iZbc(ii,Nvol)
        Rwc(nn,mm) = iRbc(ii,Mvol) ; Zws(nn,mm) = iZbs(ii,Mvol)
        Rws(nn,mm) = iRbs(ii,Mvol) ; Zwc(nn,mm) = iZbc(ii,Mvol)
     enddo ! end of do ii = 1, mn_field;

     if( Lfreebound.eq.1 ) then
      do ii = 1, mn_field ; 
        mm = im_field(ii) ; nn = in_field(ii) / Nfp 
        Vns(nn,mm) = iVns(ii)
        Bns(nn,mm) = iBns(ii)
        Vnc(nn,mm) = iVnc(ii) 
        Bnc(nn,mm) = iBnc(ii)
      enddo ! end of do ii = 1, mn_field;
     endif
    endif ! end of if( Lfreebound.eq.1 .or. . . . ) ;
  
    !write(iunit,'(" Rac         = ",99es23.15)') Rac(0:Ntor)
    !write(iunit,'(" Zas         = ",99es23.15)') Zas(0:Ntor)
    !write(iunit,'(" Ras         = ",99es23.15)') Ras(0:Ntor)
    !write(iunit,'(" Zac         = ",99es23.15)') Zac(0:Ntor)
  
   write(iunit,'(" Rac         = ",99es23.15)') iRbc(1:Ntor+1,0)
   write(iunit,'(" Zas         = ",99es23.15)') iZbs(1:Ntor+1,0)
   write(iunit,'(" Ras         = ",99es23.15)') iRbs(1:Ntor+1,0)
   write(iunit,'(" Zac         = ",99es23.15)') iZbc(1:Ntor+1,0)
  
  if( Lboundary.eq.0 ) then
    do mm = 0, Mpol_field ! will write out the plasma boundary harmonics;
     do nn = -Ntor_field, Ntor_field
  
      if( mm.eq.0 .and. nn.lt.0 ) cycle ! these modes are always excluded; 13 Oct 12;
  
      select case( mm )
      case(   0:  9 )
       if( nn.lt.- 9 .and. nn.gt.-99 )     write(iunit,1000) nn, mm, Rbc(nn,mm), nn, mm, Zbs(nn,mm), nn, mm, Rbs(nn,mm), nn, mm, Zbc(nn,mm)
       if( nn.lt.  0 .and. nn.ge.- 9 )     write(iunit,1001) nn, mm, Rbc(nn,mm), nn, mm, Zbs(nn,mm), nn, mm, Rbs(nn,mm), nn, mm, Zbc(nn,mm)
       if( nn.ge.  0 .and. nn.le.  9 )     write(iunit,1002) nn, mm, Rbc(nn,mm), nn, mm, Zbs(nn,mm), nn, mm, Rbs(nn,mm), nn, mm, Zbc(nn,mm)
       if( nn.gt.  9 .and. nn.le. 99 )     write(iunit,1001) nn, mm, Rbc(nn,mm), nn, mm, Zbs(nn,mm), nn, mm, Rbs(nn,mm), nn, mm, Zbc(nn,mm)
      case(  10: 99 )
       if( nn.lt.- 9 .and. nn.gt.-99 )     write(iunit,1003) nn, mm, Rbc(nn,mm), nn, mm, Zbs(nn,mm), nn, mm, Rbs(nn,mm), nn, mm, Zbc(nn,mm)
       if( nn.lt.  0 .and. nn.ge.- 9 )     write(iunit,1004) nn, mm, Rbc(nn,mm), nn, mm, Zbs(nn,mm), nn, mm, Rbs(nn,mm), nn, mm, Zbc(nn,mm)
       if( nn.ge.  0 .and. nn.le.  9 )     write(iunit,1005) nn, mm, Rbc(nn,mm), nn, mm, Zbs(nn,mm), nn, mm, Rbs(nn,mm), nn, mm, Zbc(nn,mm)
       if( nn.gt.  9 .and. nn.le. 99 )     write(iunit,1004) nn, mm, Rbc(nn,mm), nn, mm, Zbs(nn,mm), nn, mm, Rbs(nn,mm), nn, mm, Zbc(nn,mm)
      end select ! end of select case( mm );
  
     enddo ! end of do nn;
    enddo ! end of do mm;
  else ! Lboundary.eq.1
   
    write(iunit,'(" bn(0:",i3,")         = ",99es23.15)') Ntor, ibc(0:Ntor,Nvol)
    write(iunit,'(" R0c(0:",i3,")        = ",99es23.15)') Ntor, iR0c(0:Ntor,Nvol)
    write(iunit,'(" Z0s(0:",i3,")        = ",99es23.15)') Ntor, iZ0s(0:Ntor,Nvol)

    do ii=1,mn_rho
      if( im_rho(ii).le.9 ) then
        if( in_rho(ii)/Nfp.lt.-9 )                           write(iunit,8925) in_rho(ii)/Nfp, im_rho(ii), irhoc(ii,Nvol)
        if( in_rho(ii)/Nfp.ge.-9 .and. in_rho(ii)/Nfp.lt.0 ) write(iunit,8924) in_rho(ii)/Nfp, im_rho(ii), irhoc(ii,Nvol)
        if( in_rho(ii)/Nfp.ge. 0 .and. in_rho(ii)/Nfp.lt.9 ) write(iunit,8923) in_rho(ii)/Nfp, im_rho(ii), irhoc(ii,Nvol)
        if( in_rho(ii)/Nfp.gt. 9 )                           write(iunit,8924) in_rho(ii)/Nfp, im_rho(ii), irhoc(ii,Nvol)
      else
        if( in_rho(ii)/Nfp.lt.-9 )                           write(iunit,8935) in_rho(ii)/Nfp, im_rho(ii), irhoc(ii,Nvol)
        if( in_rho(ii)/Nfp.ge.-9 .and. in_rho(ii)/Nfp.lt.0 ) write(iunit,8934) in_rho(ii)/Nfp, im_rho(ii), irhoc(ii,Nvol)
        if( in_rho(ii)/Nfp.ge. 0 .and. in_rho(ii)/Nfp.lt.9 ) write(iunit,8933) in_rho(ii)/Nfp, im_rho(ii), irhoc(ii,Nvol)
        if( in_rho(ii)/Nfp.gt. 9 )                           write(iunit,8934) in_rho(ii)/Nfp, im_rho(ii), irhoc(ii,Nvol)
      endif
    enddo

  endif ! Lboundary.eq.0
  
  8923 format(" rhomn(",i1,","i1") = ",es23.15)
  8924 format(" rhomn(",i2,","i1") = ",es23.15)
  8925 format(" rhomn(",i3,","i1") = ",es23.15)
  8933 format(" rhomn(",i1,","i2") = ",es23.15)
  8934 format(" rhomn(",i2,","i2") = ",es23.15)
  8935 format(" rhomn(",i3,","i2") = ",es23.15)


  if( Lfreebound.eq.1 ) then
    do mm = 0, Mpol ! will write out the computation domain harmonics; (only relevant in free-boundary case);
     do nn = -Ntor, Ntor
  
      if( mm.eq.0 .and. nn.lt.0 ) cycle ! these modes are always excluded; 13 Oct 12;
  
      select case( mm )
      case(   0:  9 )
       if( nn.lt.- 9 .and. nn.gt.-99 ) write(iunit,1010) nn, mm, Rwc(nn,mm), nn, mm, Zws(nn,mm), nn, mm, Rws(nn,mm), nn, mm, Zwc(nn,mm)
       if( nn.lt.  0 .and. nn.ge.- 9 ) write(iunit,1011) nn, mm, Rwc(nn,mm), nn, mm, Zws(nn,mm), nn, mm, Rws(nn,mm), nn, mm, Zwc(nn,mm)
       if( nn.ge.  0 .and. nn.le.  9 ) write(iunit,1012) nn, mm, Rwc(nn,mm), nn, mm, Zws(nn,mm), nn, mm, Rws(nn,mm), nn, mm, Zwc(nn,mm)
       if( nn.gt.  9 .and. nn.le. 99 ) write(iunit,1011) nn, mm, Rwc(nn,mm), nn, mm, Zws(nn,mm), nn, mm, Rws(nn,mm), nn, mm, Zwc(nn,mm)
      case(  10: 99 )
       if( nn.lt.- 9 .and. nn.gt.-99 ) write(iunit,1013) nn, mm, Rwc(nn,mm), nn, mm, Zws(nn,mm), nn, mm, Rws(nn,mm), nn, mm, Zwc(nn,mm)
       if( nn.lt.  0 .and. nn.ge.- 9 ) write(iunit,1014) nn, mm, Rwc(nn,mm), nn, mm, Zws(nn,mm), nn, mm, Rws(nn,mm), nn, mm, Zwc(nn,mm)
       if( nn.ge.  0 .and. nn.le.  9 ) write(iunit,1015) nn, mm, Rwc(nn,mm), nn, mm, Zws(nn,mm), nn, mm, Rws(nn,mm), nn, mm, Zwc(nn,mm)
       if( nn.gt.  9 .and. nn.le. 99 ) write(iunit,1014) nn, mm, Rwc(nn,mm), nn, mm, Zws(nn,mm), nn, mm, Rws(nn,mm), nn, mm, Zwc(nn,mm)
      end select ! end of select case( mm );
  
     enddo ! end of do nn;
    enddo ! end of do mm;
  
    do mm = 0, Mpol ! will write out the computation domain harmonics; (only relevant in free-boundary case);
      do nn = -Ntor, Ntor
    
        if( mm.eq.0 .and. nn.lt.0 ) cycle ! these modes are always excluded; 13 Oct 12;
    
        select case( mm )
        case(   0:  9 )
        if( nn.lt.- 9 .and. nn.gt.-99 ) write(iunit,1020) nn, mm, Vns(nn,mm), nn, mm, Bns(nn,mm), nn, mm, Vnc(nn,mm), nn, mm, Bnc(nn,mm)
        if( nn.lt.  0 .and. nn.ge.- 9 ) write(iunit,1021) nn, mm, Vns(nn,mm), nn, mm, Bns(nn,mm), nn, mm, Vnc(nn,mm), nn, mm, Bnc(nn,mm)
        if( nn.ge.  0 .and. nn.le.  9 ) write(iunit,1022) nn, mm, Vns(nn,mm), nn, mm, Bns(nn,mm), nn, mm, Vnc(nn,mm), nn, mm, Bnc(nn,mm)
        if( nn.gt.  9 .and. nn.le. 99 ) write(iunit,1021) nn, mm, Vns(nn,mm), nn, mm, Bns(nn,mm), nn, mm, Vnc(nn,mm), nn, mm, Bnc(nn,mm)
        case(  10: 99 )
        if( nn.lt.- 9 .and. nn.gt.-99 ) write(iunit,1023) nn, mm, Vns(nn,mm), nn, mm, Bns(nn,mm), nn, mm, Vnc(nn,mm), nn, mm, Bnc(nn,mm)
        if( nn.lt.  0 .and. nn.ge.- 9 ) write(iunit,1024) nn, mm, Vns(nn,mm), nn, mm, Bns(nn,mm), nn, mm, Vnc(nn,mm), nn, mm, Bnc(nn,mm)
        if( nn.ge.  0 .and. nn.le.  9 ) write(iunit,1025) nn, mm, Vns(nn,mm), nn, mm, Bns(nn,mm), nn, mm, Vnc(nn,mm), nn, mm, Bnc(nn,mm)
        if( nn.gt.  9 .and. nn.le. 99 ) write(iunit,1024) nn, mm, Vns(nn,mm), nn, mm, Bns(nn,mm), nn, mm, Vnc(nn,mm), nn, mm, Bnc(nn,mm)
        end select ! end of select case( mm );
    
      enddo ! end of do nn;
    enddo ! end of do mm;
  endif
  
  1000 format("Rbc(",i3,",",i1,")",2x,"=",es23.15," Zbs(",i3,",",i1,")",2x,"=",es23.15," Rbs(",i3,",",i1,")",2x,"=",es23.15," Zbc(",i3,",",i1,")",2x,"=",es23.15)
  1001 format("Rbc(",i2,",",i1,")",3x,"=",es23.15," Zbs(",i2,",",i1,")",3x,"=",es23.15," Rbs(",i2,",",i1,")",3x,"=",es23.15," Zbc(",i2,",",i1,")",3x,"=",es23.15)
  1002 format("Rbc(",i1,",",i1,")",4x,"=",es23.15," Zbs(",i1,",",i1,")",4x,"=",es23.15," Rbs(",i1,",",i1,")",4x,"=",es23.15," Zbc(",i1,",",i1,")",4x,"=",es23.15)
  1003 format("Rbc(",i3,",",i2,")",1x,"=",es23.15," Zbs(",i3,",",i2,")",1x,"=",es23.15," Rbs(",i3,",",i2,")",1x,"=",es23.15," Zbc(",i3,",",i2,")",1x,"=",es23.15)
  1004 format("Rbc(",i2,",",i2,")",2x,"=",es23.15," Zbs(",i2,",",i2,")",2x,"=",es23.15," Rbs(",i2,",",i2,")",2x,"=",es23.15," Zbc(",i2,",",i2,")",2x,"=",es23.15)
  1005 format("Rbc(",i1,",",i2,")",3x,"=",es23.15," Zbs(",i1,",",i2,")",3x,"=",es23.15," Rbs(",i1,",",i2,")",3x,"=",es23.15," Zbc(",i1,",",i2,")",3x,"=",es23.15)
  
  1010 format("Rwc(",i3,",",i1,")",2x,"=",es23.15," Zws(",i3,",",i1,")",2x,"=",es23.15," Rws(",i3,",",i1,")",2x,"=",es23.15," Zwc(",i3,",",i1,")",2x,"=",es23.15)
  1011 format("Rwc(",i2,",",i1,")",3x,"=",es23.15," Zws(",i2,",",i1,")",3x,"=",es23.15," Rws(",i2,",",i1,")",3x,"=",es23.15," Zwc(",i2,",",i1,")",3x,"=",es23.15)
  1012 format("Rwc(",i1,",",i1,")",4x,"=",es23.15," Zws(",i1,",",i1,")",4x,"=",es23.15," Rws(",i1,",",i1,")",4x,"=",es23.15," Zwc(",i1,",",i1,")",4x,"=",es23.15)
  1013 format("Rwc(",i3,",",i2,")",1x,"=",es23.15," Zws(",i3,",",i2,")",1x,"=",es23.15," Rws(",i3,",",i2,")",1x,"=",es23.15," Zwc(",i3,",",i2,")",1x,"=",es23.15)
  1014 format("Rwc(",i2,",",i2,")",2x,"=",es23.15," Zws(",i2,",",i2,")",2x,"=",es23.15," Rws(",i2,",",i2,")",2x,"=",es23.15," Zwc(",i2,",",i2,")",2x,"=",es23.15)
  1015 format("Rwc(",i1,",",i2,")",3x,"=",es23.15," Zws(",i1,",",i2,")",3x,"=",es23.15," Rws(",i1,",",i2,")",3x,"=",es23.15," Zwc(",i1,",",i2,")",3x,"=",es23.15)
  
  1020 format("Vns(",i3,",",i1,")",2x,"=",es23.15," Bns(",i3,",",i1,")",2x,"=",es23.15," Vnc(",i3,",",i1,")",2x,"=",es23.15," Bnc(",i3,",",i1,")",2x,"=",es23.15)
  1021 format("Vns(",i2,",",i1,")",3x,"=",es23.15," Bns(",i2,",",i1,")",3x,"=",es23.15," Vnc(",i2,",",i1,")",3x,"=",es23.15," Bnc(",i2,",",i1,")",3x,"=",es23.15)
  1022 format("Vns(",i1,",",i1,")",4x,"=",es23.15," Bns(",i1,",",i1,")",4x,"=",es23.15," Vnc(",i1,",",i1,")",4x,"=",es23.15," Bnc(",i1,",",i1,")",4x,"=",es23.15)
  1023 format("Vns(",i3,",",i2,")",1x,"=",es23.15," Bns(",i3,",",i2,")",1x,"=",es23.15," Vnc(",i3,",",i2,")",1x,"=",es23.15," Bnc(",i3,",",i2,")",1x,"=",es23.15)
  1024 format("Vns(",i2,",",i2,")",2x,"=",es23.15," Bns(",i2,",",i2,")",2x,"=",es23.15," Vnc(",i2,",",i2,")",2x,"=",es23.15," Bnc(",i2,",",i2,")",2x,"=",es23.15)
  1025 format("Vns(",i1,",",i2,")",3x,"=",es23.15," Bns(",i1,",",i2,")",3x,"=",es23.15," Vnc(",i1,",",i2,")",3x,"=",es23.15," Bnc(",i1,",",i2,")",3x,"=",es23.15)
  
    write(iunit,'("/")')
  
    if( Wwrtend ) then ; cput = GETTIME ; write(ounit,'("wrtend : ",f10.2," : myid=",i3," ; writing numericlist ;")') cput-cpus, myid
    endif
  
    write(iunit,'("&numericlist")')
    write(iunit,'(" Linitialize = ",i9            )') Linitialize
    write(iunit,'(" LautoinitBn = ",i9            )') LautoinitBn
    write(iunit,'(" Lzerovac    = ",i9            )') Lzerovac
    write(iunit,'(" Ndiscrete   = ",i9            )') Ndiscrete
    write(iunit,'(" Nquad       = ",i9            )') Nquad
    write(iunit,'(" iMpol       = ",i9            )') iMpol
    write(iunit,'(" iNtor       = ",i9            )') iNtor
    write(iunit,'(" Lsparse     = ",i9            )') Lsparse
    write(iunit,'(" Lsvdiota    = ",i9            )') Lsvdiota
    write(iunit,'(" imethod     = ",i9            )') imethod
    write(iunit,'(" iorder      = ",i9            )') iorder
    write(iunit,'(" iprecon     = ",i9            )') iprecon
    write(iunit,'(" iotatol     = ",es23.15       )') iotatol
    write(iunit,'(" Lextrap     = ",i9            )') Lextrap
    write(iunit,'(" Mregular    = ",i9            )') Mregular
    write(iunit,'(" Lrzaxis     = ",i9            )') Lrzaxis
    write(iunit,'(" Ntoraxis    = ",i9            )') Ntoraxis
    write(iunit,'("/")')
  
    if( Wwrtend ) then ; cput = GETTIME ; write(ounit,'("wrtend : ",f10.2," : myid=",i3," ; writing locallist ;")') cput-cpus, myid
    endif
  
    write(iunit,'("&locallist")')
    write(iunit,'(" LBeltrami   = ",i9            )') LBeltrami
    write(iunit,'(" Linitgues   = ",i9            )') Linitgues
    write(iunit,'(" Lmatsolver  = ",i9            )') Lmatsolver
    write(iunit,'(" NiterGMRES  = ",i9            )') NiterGMRES
    write(iunit,'(" LGMRESprec  = ",i9            )') LGMRESprec
    write(iunit,'(" epsGMRES    = ",es23.15       )') epsGMRES
    write(iunit,'(" epsILU      = ",es23.15       )') epsILU
  
   !write(iunit,'(" Lposdef     = ",i9            )') Lposdef ! redundant;
   !write(iunit,'(" Nmaxexp     = ",i9            )') Nmaxexp
    write(iunit,'("/")')
  
    if( Wwrtend ) then ; cput = GETTIME ; write(ounit,'("wrtend : ",f10.2," : myid=",i3," ; writing globallist ;")') cput-cpus, myid
    endif
  
    write(iunit,'("&globallist")')
    write(iunit,'(" Lfindzero   = ",i9            )') Lfindzero
    write(iunit,'(" escale      = ",es23.15       )') escale
    write(iunit,'(" opsilon     = ",es23.15       )') opsilon
    write(iunit,'(" pcondense   = ",es23.15       )') pcondense
    write(iunit,'(" epsilon     = ",es23.15       )') epsilon
    write(iunit,'(" wpoloidal   = ",es23.15       )') wpoloidal
    write(iunit,'(" upsilon     = ",es23.15       )') upsilon
    write(iunit,'(" forcetol    = ",es23.15       )') forcetol
    write(iunit,'(" c05xmax     = ",es23.15       )') c05xmax
    write(iunit,'(" c05xtol     = ",es23.15       )') c05xtol
    write(iunit,'(" c05factor   = ",es23.15       )') c05factor
    write(iunit,'(" LreadGF     = ",L9            )') LreadGF
    write(iunit,'(" mfreeits    = ",i9            )') mfreeits
    write(iunit,'(" gBntol      = ",es23.15       )') gBntol
    write(iunit,'(" gBnbld      = ",es23.15       )') gBnbld
    write(iunit,'(" vcasingeps  = ",es23.15       )') vcasingeps
    write(iunit,'(" vcasingtol  = ",es23.15       )') vcasingtol
    write(iunit,'(" vcasingits  = ",i9            )') vcasingits
    write(iunit,'(" vcasingper  = ",i9            )') vcasingper
    write(iunit,'("/")')
  
    if( Wwrtend ) then ; cput = GETTIME ; write(ounit,'("wrtend : ",f10.2," : myid=",i3," ; writing diagnosticslist ;")') cput-cpus, myid
    endif
  
    write(iunit,'("&diagnosticslist")')
    write(iunit,'(" odetol      = ",es23.15       )') odetol
   !write(iunit,'(" absreq      = ",es23.15       )') absreq
   !write(iunit,'(" relreq      = ",es23.15       )') relreq
   !write(iunit,'(" absacc      = ",es23.15       )') absacc
   !write(iunit,'(" epsr        = ",es23.15       )') epsr
    write(iunit,'(" nPpts       = ",i9            )') nPpts
    write(iunit,'(" Ppts        = ",es23.15       )') Ppts
    write(iunit,'(" nPtrj       = ",256i6         )') nPtrj(1:Mvol)
    write(iunit,'(" LHevalues   = ",L9            )') LHevalues
    write(iunit,'(" LHevectors  = ",L9            )') LHevectors
    write(iunit,'(" LHmatrix    = ",L9            )') LHmatrix
    write(iunit,'(" Lperturbed  = ",i9            )') Lperturbed
    write(iunit,'(" dpp         = ",i9            )') dpp
    write(iunit,'(" dqq         = ",i9            )') dqq
    write(iunit,'(" dRZ         = ",es23.15       )') dRZ
    write(iunit,'(" Lcheck      = ",i9            )') Lcheck
    write(iunit,'(" Ltiming     = ",L9            )') Ltiming
    write(iunit,'("/")')
  
    if( Wwrtend ) then ; cput = GETTIME ; write(ounit,'("wrtend : ",f10.2," : myid=",i3," ; writing screenlist ;")') cput-cpus, myid
    endif
  
    write(iunit,'("&screenlist")')
  ! WSCREENLIST ! write screenlist; this is expanded by Makefile ; do not remove;
    if( Wreadin           ) write(iunit,'(" Wreadin = ",L1                )') Wreadin
    if( Wwrtend           ) write(iunit,'(" Wwrtend = ",L1                )') Wwrtend
    if( Wmacros           ) write(iunit,'(" Wmacros = ",L1                )') Wmacros
    write(iunit,'("/")')
  
#ifdef DEBUG
    FATAL( wrtend, .not.allocated(iRbc), error )
    FATAL( wrtend, .not.allocated(iZbs), error )
    FATAL( wrtend, .not.allocated(iRbs), error )
    FATAL( wrtend, .not.allocated(iZbc), error )
#endif
  
    ! write initial guess of interface geometry
    if( Lboundary.eq.0 ) then
      do imn = 1, mn_field ; write(iunit,'(2i6,1024es23.15)') im_field(imn), in_field(imn)/Nfp, ( iRbc(imn,vvol), iZbs(imn,vvol), iRbs(imn,vvol), iZbc(imn,vvol), vvol = 1, Nvol )
      enddo
    else    
      ! Backward map last geometry to henneberg's representation
      do vvol=1, Nvol-1  
        call backwardMap( iRbc(1:mn_field,vvol), iZbs(1:mn_field,vvol), irhoc(1:mn_rho,vvol), ibc(0:Ntor,vvol), iR0c(0:Ntor,vvol), iZ0s(1:Ntor,vvol) )
      enddo !vvol=1:Nvol-1
  
      do nn = 0, Ntor
        write(iunit,'(2i6, 1024es23.15)') 0, nn, ( ibc(nn,vvol), iR0c(nn,vvol), iZ0s(nn,vvol), zero, vvol=1,Nvol )
      enddo! nn=0,Ntor
  
      do imn = 1, mn_rho
        write(iunit,'(2i6, 1024es23.15)') im_rho(imn), in_rho(imn)/Nfp, ( zero, zero, zero, irhoc( imn, vvol ), vvol=1,Nvol )
      enddo!imn=1,mn_rho
  
    endif! Lboundary.eq.0
  
    close(iunit)
  
#ifdef DEBUG
    if( Wwrtend ) then ; cput = GETTIME ; write(ounit,'("wrtend : ",f10.2," : myid=",i3," ; wrote ext.sp.end ;")') cput-cpus, myid
    endif
#endif
  
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
    RETURN(wrtend)
  
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  end subroutine wrtend
