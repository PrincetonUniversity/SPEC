!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (output) ! Writes all the output information to ext.sp.h5.

!latex \briefly{All the output information is contained in \type{ext.sp.h5}.}
!latex \calledby{\link{xspech}}
!l tex \calls{\link{}}

!latex \tableofcontents

!latex \newcommand{\pb}[1]{\parbox[t]{13cm}{#1}}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine hdfint
  
  use constants, only : 
  use numerical, only : 
  use fileunits, only : ounit
  use inputlist
  use cputiming, only : Thdfint
  use allglobal, only : myid, ncpu, cpus, &
                        Mvol, ForceErr, &
                        mn, im, in, iRbc, iZbs, iRbs, iZbc, &
                        dRbc, dZbs, dRbs, dZbc, &
                        vvolume, dvolume, &
                        Bsupumn, Bsupvmn, &
                        Btemn, Bzemn, Btomn, Bzomn, &
                        iVns, iBns, iVnc, iBnc, &
                        lmns, &
                        TT, &
                        beltramierror, &
                        IPDt
  
#ifdef NOHDF5
#else
  use hdf5
#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER                        :: Mrad
  REAL                           :: tvolume
  
#ifdef NOHDF5
#else
  INTEGER                        :: hdfier, rank
  integer(hid_t)                 :: file_id, space_id, dset_id
  integer(hsize_t)               :: onedims(1:1), twodims(1:2), threedims(1:3)
#endif
  
  BEGIN(hdfint)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef NOHDF5

#else

  call h5open_f( hdfier ) ! initialize Fortran interface;

  FATAL( hdfint, hdfier.ne.0, error calling h5open_f )

  call h5fcreate_f( trim(ext)//".sp.h5", H5F_ACC_TRUNC_F, file_id, hdfier ) ! create new file;

  FATAL( hdfint, hdfier.ne.0, error calling h5fcreate_f )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   HWRITEIV(           1, Igeometry         , Igeometry           )
   HWRITEIV(           1, Istellsym         , Istellsym           )
   HWRITEIV(           1, Lfreebound        , Lfreebound          )
   HWRITERV(           1, phiedge           , phiedge             )
   HWRITERV(           1, curtor            , curtor              )
   HWRITERV(           1, curpol            , curpol              )
   HWRITERV(           1, gamma             , gamma               )
   HWRITEIV(           1, Nfp               , Nfp                 )
   HWRITEIV(           1, Nvol              , Nvol                )
   HWRITEIV(           1, Mpol              , Mpol                )
   HWRITEIV(           1, Ntor              , Ntor                )
   HWRITEIV(        Mvol, Lrad              ,      Lrad(1:Mvol)   )
   HWRITEIV(           1, Lconstraint       , Lconstraint         )
   HWRITERV(        Mvol, tflux             ,     tflux(1:Mvol)   )
   HWRITERV(        Mvol, pflux             ,     pflux(1:Mvol)   )
   HWRITERV(        Nvol, helicity          ,  helicity(1:Nvol)   )
   HWRITERV(           1, pscale            , pscale              )
   HWRITERV(        Nvol, pressure          ,  pressure(1:Nvol)   )
   HWRITEIV(           1, Ladiabatic        , Ladiabatic          )
   HWRITERV(        Mvol, adiabatic         , adiabatic(1:Nvol)   )
   HWRITERV(      1+Nvol, mu                ,        mu(1:Mvol)   )
   HWRITEIV(      1+Mvol, pl                ,        pl(0:Nvol)   )
   HWRITEIV(      1+Mvol, ql                ,        ql(0:Nvol)   )
   HWRITEIV(      1+Mvol, pr                ,        pr(0:Nvol)   )
   HWRITEIV(      1+Mvol, qr                ,        qr(0:Nvol)   )
   HWRITERV(      1+Nvol, iota              ,      iota(0:Nvol)   )
   HWRITEIV(      1+Mvol, lp                ,        lp(0:Nvol)   )
   HWRITEIV(      1+Mvol, lq                ,        lq(0:Nvol)   )
   HWRITEIV(      1+Mvol, rp                ,        rp(0:Nvol)   )
   HWRITEIV(      1+Mvol, rq                ,        rq(0:Nvol)   )
   HWRITERV(      1+Nvol, oita              ,      oita(0:Nvol)   )
   HWRITERV(           1, mupftol           , mupftol             )
   HWRITEIV(           1, mupfits           , mupfits             )
   HWRITERV(        Mvol, IPDt              , IPDt(1:Mvol-1)        )
!  HWRITERV(     MNtor+1, Rac               , Rac(0:MNtor)        )
!  HWRITERV(     MNtor+1, Zas               , Zas(0:MNtor)        )
!  HWRITERV(     MNtor+1, Ras               , Ras(0:MNtor)        )
!  HWRITERV(     MNtor+1, Zac               , Zac(0:MNtor)        )
!  HWRITERA(   2*MNtor+1, 2*MMpol+1,  Rbc, Rbc(-MNtor:MNtor,-MMpol:MMpol) )
!  HWRITERA(   2*MNtor+1, 2*MMpol+1,  Zbs, Zbs(-MNtor:MNtor,-MMpol:MMpol) )
!  HWRITERA(   2*MNtor+1, 2*MMpol+1,  Rbs, Rbs(-MNtor:MNtor,-MMpol:MMpol) )
!  HWRITERA(   2*MNtor+1, 2*MMpol+1,  Zbc, Zbc(-MNtor:MNtor,-MMpol:MMpol) )
!  HWRITERA(   2*MNtor+1, 2*MMpol+1,  Rwc, Rwc(-MNtor:MNtor,-MMpol:MMpol) )
!  HWRITERA(   2*MNtor+1, 2*MMpol+1,  Zws, Zws(-MNtor:MNtor,-MMpol:MMpol) )
!  HWRITERA(   2*MNtor+1, 2*MMpol+1,  Rws, Rws(-MNtor:MNtor,-MMpol:MMpol) )
!  HWRITERA(   2*MNtor+1, 2*MMpol+1,  Zwc, Zwc(-MNtor:MNtor,-MMpol:MMpol) )
!  HWRITERA(   2*MNtor+1, 2*MMpol+1,  Vns, Vns(-MNtor:MNtor,-MMpol:MMpol) )
!  HWRITERA(   2*MNtor+1, 2*MMpol+1,  Bns, Bns(-MNtor:MNtor,-MMpol:MMpol) )
!  HWRITERA(   2*MNtor+1, 2*MMpol+1,  Vnc, Vnc(-MNtor:MNtor,-MMpol:MMpol) )
!  HWRITERA(   2*MNtor+1, 2*MMpol+1,  Bnc, Bnc(-MNtor:MNtor,-MMpol:MMpol) )

  HWRITEIV( 1, Lperturbed, Lperturbed )
  HWRITEIV( 1, dpp       , dpp        )
  HWRITEIV( 1, dqq       , dqq        )

  HWRITERV( mn, Vns, iVns(1:mn) )
  HWRITERV( mn, Bns, iBns(1:mn) )
  HWRITERV( mn, Vnc, iVnc(1:mn) )
  HWRITERV( mn, Bnc, iBnc(1:mn) )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \begin{enumerate}

!latex \item In addition to the input variables, which are described in \link{global}, the following quantities are written to \type{ext.sp.h5} :
!latex 
!latex \begin{tabular}{|l|l|l|} \hline

!latex \type{variable}               & type    & \pb{description} \\ \hline

!latex \type{mn}                     & integer & \pb{number of Fourier modes} \\
  HWRITEIV( 1, mn, mn )                          
!latex \type{im(1:mn)}               & integer & \pb{poloidal mode numbers} \\
  HWRITEIV( mn, im, im(1:mn) )                   
!latex \type{in(1:mn)}               & integer & \pb{toroidal mode numbers} \\
  HWRITEIV( mn, in, in(1:mn) )                   
!latex \type{Mvol}                   & integer & \pb{number of interfaces = number of volumes} \\
  HWRITEIV( 1, Mvol, Mvol )                      
!latex \type{iRbc(1:mn,0:Mvol)}      & real    & \pb{Fourier harmonics, $R_{m,n}$, of interfaces} \\
  HWRITERA( mn, Mvol+1, Rbc, iRbc(1:mn,0:Mvol) ) 
!latex \type{iZbs(1:mn,0:Mvol)}      & real    & \pb{Fourier harmonics, $Z_{m,n}$, of interfaces} \\
  HWRITERA( mn, Mvol+1, Zbs, iZbs(1:mn,0:Mvol) ) 
!latex \type{iRbs(1:mn,0:Mvol)}      & real    & \pb{Fourier harmonics, $R_{m,n}$, of interfaces} \\
  HWRITERA( mn, Mvol+1, Rbs, iRbs(1:mn,0:Mvol) ) 
!latex \type{iZbc(1:mn,0:Mvol)}      & real    & \pb{Fourier harmonics, $Z_{m,n}$, of interfaces} \\
  HWRITERA( mn, Mvol+1, Zbc, iZbc(1:mn,0:Mvol) ) 
!latex \type{forcetol}               & real    & \pb{force-balance error across interfaces} \\
  HWRITERV( 1, forcetol, forcetol )               
!latex \type{ForceErr}               & real    & \pb{force-balance error across interfaces} \\
  HWRITERV( 1, ForceErr, ForceErr )              

  if( Lcheck.eq.1 ) then
!latex \type{beltramierror}          & real    & \pb{error in beltrami field (volume integral)} \\ 
   HWRITERA( Mvol, 3, beltramierror, beltramierror(1:Mvol,1:3) )  
  endif
  
  if( allocated(vvolume) ) then ! why is it required to confirm that vvolume has been allocated ; 24 Nov 16;
   
   tvolume = sum(vvolume(1:Nvol) )
!latex \type{volume}                 & real    & \pb{total volume = $\sum V_v$} \\
   HWRITERV( 1, volume, tvolume) 
   
  else
   
   if( Whdfint ) write(ounit,'("hdfint : ", 10x ," : myid=",i3," ; vvolume is not allocated ;")') myid
   
  endif ! end of if( allocated(vvolume) ) ; 11 Aug 14;

  Mrad  = maxval( Lrad(1:Mvol) )
!latex \type{Mrad}                   & integer & \pb{the maximum radial (Chebyshev) resolution}\\
  HWRITEIV( 1, Mrad, Mrad )
!latex \type{TT(0:Mrad,0:1,0:1)}     & real    & \pb{the Chebyshev polynomials, $T_l$, and their derivatives, evaluated at $s=\pm 1$}\\
  HWRITERC( Mrad+1, 2, 2, TT, TT, (0:Mrad,0:1,0:1) )
!latex \type{Btemn(1:mn,0:1,1:Mvol)} & real    & \pb{the cosine harmonics of the covariant poloidal field, \\
!latex                                           i.e. $[[B_{\t,j}]]$ evaluated on the inner and outer interface in each volume} \\
  HWRITERC( mn, 2, Mvol, Btemn, Btemn(1:mn,0:1,1:Mvol) )
!latex \type{Bzemn(1:mn,0:1,1:Mvol)} & real    & \pb{the cosine harmonics of the covariant toroidal field, \\
!latex                                           i.e. $[[B_{\z,j}]]$ evaluated on the inner and outer interface in each volume} \\
  HWRITERC( mn, 2, Mvol, Bzemn, Bzemn(1:mn,0:1,1:Mvol) )
!latex \type{Btomn(1:mn,0:1,1:Mvol)} & real    & \pb{the sine harmonics of the covariant poloidal field, \\
!latex                                           i.e. $[[B_{\t,j}]]$ evaluated on the inner and outer interface in each volume} \\
  HWRITERC( mn, 2, Mvol, Btomn, Btomn(1:mn,0:1,1:Mvol) )
!latex \type{Bzomn(1:mn,0:1,1:Mvol)} & real    & \pb{the sine harmonics of the covariant toroidal field, \\
!latex                                           i.e. $[[B_{\z,j}]]$ evaluated on the inner and outer interface in each volume} \\
  HWRITERC( mn, 2, Mvol, Bzomn, Bzomn(1:mn,0:1,1:Mvol) )

  if( Lperturbed.eq.1 ) then

!latex \type{dRbc(1:mn,0:Nvol)}      & real    & \pb{Fourier harmonics, $R_{j}$, of interfaces; linearly perturbed solution} \\
  HWRITERA( mn, Nvol+1, dRbc, dRbc(1:mn,0:Nvol) )
!latex \type{dZbs(1:mn,0:Nvol)}      & real    & \pb{Fourier harmonics, $Z_{j}$, of interfaces; linearly perturbed solution} \\
  HWRITERA( mn, Nvol+1, dZbs, dZbs(1:mn,0:Nvol) )
!latex \type{dRbs(1:mn,0:Nvol)}      & real    & \pb{Fourier harmonics, $R_{j}$, of interfaces; linearly perturbed solution} \\
  HWRITERA( mn, Nvol+1, dRbs, dRbs(1:mn,0:Nvol) )
!latex \type{dZbc(1:mn,0:Nvol)}      & real    & \pb{Fourier harmonics, $Z_{j}$, of interfaces; linearly perturbed solution} \\
  HWRITERA( mn, Nvol+1, dZbc, dZbc(1:mn,0:Nvol) )

  endif

!latex \type{lmns}                   & integer & \pb{resolution of straight fieldline transformation} \\
  HWRITEIV( 1, lmns, lmns ) 
  
!latex \hline \end{tabular}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item All quantities marked as real should be treated as double precision.

!latex \end{enumerate}
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  call h5fclose_f( file_id, hdfier ) ! terminate access;

  FATAL( hdfint, hdfier.ne.0, error calling h5fclose_f )

  call h5close_f( hdfier ) ! close Fortran interface;

  FATAL( hdfint, hdfier.ne.0, error calling h5close_f )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#endif
  
  RETURN(hdfint)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine hdfint

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
