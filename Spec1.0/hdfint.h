!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Writes ext.h5 file.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!subroutine h5fcreate_f(name, access_flags, file_id, hdferr)              ! creates an HDF5 file.;
!  CHARACTER(LEN=*), INTENT(IN)  :: name                                  ! name of the file;
!  INTEGER         , INTENT(IN)  :: access_flag ! file access flags: H5F_ACC_RDWR_F, H5F_ACC_RDONLY_F, H5F_ACC_TRUNC_F, H5F_ACC_EXCL_F, H5F_ACC_DEBUG_F;
!  INTEGER(HID_T)  , INTENT(OUT) :: file_id                               ! file identifier;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!subroutine h5screate_simple_f(rank, dims, space_id, hdferr)              ! creates a new simple dataspace and opens it for access;
!  INTEGER         , INTENT(IN)  :: rank                                  ! number of dataspace dimensions;
!  INTEGER(HSIZE_T), INTENT(IN)  :: dims(*)                               ! array with current dimension sizes;
!  INTEGER(HID_T)  , INTENT(OUT) :: space_id                              ! dataspace identifier;

!subroutine h5dcreate_f(loc_id, name, type_id, space_id, dset_id, hdferr) ! creates a new dataset and links it to a location in the file;
!  INTEGER(HID_T)  , INTENT(IN)  :: loc_id                                ! file or group identifier 
!  CHARACTER(LEN=*), INTENT(IN)  :: name                                  ! name of the dataset 
!  INTEGER(HID_T)  , INTENT(IN)  :: type_id                               ! datatype identifier 
!  INTEGER(HID_T)  , INTENT(IN)  :: space_id                              ! dataspace identifier 
!  INTEGER(HID_T)  , INTENT(OUT) :: dset_id                               ! dataset identifier 

!subroutine h5dwrite_f(dset_id, mem_type_id, buf, dims, hdferr )          ! writes raw data from a buffer to a dataset;
!  INTEGER(HID_T)  , INTENT(IN)  :: dset_id                               ! dataset identifier
!  INTEGER(HID_T)  , INTENT(IN)  :: mem_type_id                           ! memory datatype identifier
!  TYPE            , INTENT(IN)  :: buf                                   ! data buffer; may be a scalar or an array
!  INTEGER(HSIZE_T), INTENT(IN)  :: dims(*) ! array to hold corresponding dimension sizes of data buffer buf;
!                                             dim(k) has value of the k-th dimension of buffer buf; values are ignored if buf is a scalar

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine hdfint
  
  use constants, only : zero, goldenmean
  use numerical, only : 
  use fileunits, only : ounit
  use inputlist, only : Wmacros, Whdfint, version, ext, Nvol, tflux, pflux, pressure, mu, helicity, iota, oita, forcetol, ForceErr, dpp, dqq, &
                        Lperturbed, pscale
  use cputiming, only : Thdfint
  use allglobal, only : myid, ncpu, cpus, &
                        Mvol, &
                        mn, im, in, iRbc, iZbs, iRbs, iZbc, & ! 10 Oct 12; 
                        dRbc, dZbs, dRbs, dZbc, &
                        vvolume, dvolume, &
                        Bsupumn, Bsupvmn, Bsubtemn, Bsubzemn, Bsubtomn, Bsubzomn, &
                        lmns, &
			beltramierror

#ifdef NOHDF5

#else

  use hdf5

#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS

  INTEGER                        :: lvol
  
#ifdef NOHDF5

#else

! the following are used by the macros HWRITE below; do not alter/remove;
  INTEGER                        :: hdfier, rank
  integer(hid_t)                 :: file_id, space_id, dset_id                 ! warning: the string "INTEGER" is a macro;
  integer(hsize_t)               :: onedims(1:1), twodims(1:2), threedims(1:3) ! warning: the string "INTEGER" is a macro;

  INTEGER                        :: imn
  REAL                           :: tvolume

#endif
  
  BEGIN(hdfint)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef NOHDF5

#else

  call h5open_f( hdfier ) ! initialize Fortran interface;

  FATALMESS(hdfint, hdfier.ne.0, error calling h5open_f )

  call h5fcreate_f( trim(ext)//".h5", H5F_ACC_TRUNC_F, file_id, hdfier ) ! create new file;

  FATALMESS(hdfint, hdfier.ne.0, error calling h5fcreate_f )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item The following quantities are written to \verb+ext.h5+ :

!latex \begin{itemize}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item \verb+version+; real; identifies executable \verb+xspec+;

  HWRITEIV( 1, version, version )

!latex \item \verb+mn+; integer;

  HWRITEIV( 1, mn, mn )

!latex \item \verb+im(1:mn)+; integer; poloidal mode numbers;

  HWRITEIV( mn, im, im(1:mn) )

!latex \item \verb+in(1:mn)+; integer; toroidal mode numbers;

  HWRITEIV( mn, in, in(1:mn) )

!latex \item \verb+Nvol+; integer; number of interfaces = number of volumes;

  HWRITEIV( 1, Nvol, Nvol )
  
!latex \item \verb+beltramierror+; real; error in beltrami field (volume integral);  11 Apr 16;

  HWRITERA( Nvol, 3, beltramierror, beltramierror(1:Nvol,1:3) )  ! replaced Mvol with Nvol; 11 Jul 17;

!latex \item \verb+iRbc(1:mn,0:Nvol)+; real; Fourier harmonics, $R_{m,n}$, of interfaces;

  HWRITERA( mn, Nvol+1, Rbc, iRbc(1:mn,0:Nvol) )

!latex \item \verb+iZbs(1:mn,0:Nvol)+; real; Fourier harmonics, $Z_{m,n}$, of interfaces;

  HWRITERA( mn, Nvol+1, Zbs, iZbs(1:mn,0:Nvol) )

!latex \item \verb+iRbs(1:mn,0:Nvol)+; real; Fourier harmonics, $R_{m,n}$, of interfaces;

  HWRITERA( mn, Nvol+1, Rbs, iRbs(1:mn,0:Nvol) )

!latex \item \verb+iZbc(1:mn,0:Nvol)+; real; Fourier harmonics, $Z_{m,n}$, of interfaces;

  HWRITERA( mn, Nvol+1, Zbc, iZbc(1:mn,0:Nvol) )

!latex \item \verb+volume+; real; volume;

  if( allocated(vvolume) ) then
   
   tvolume = sum(vvolume(1:Nvol) ) ! may or may not be neccessary; 11 Aug 14;
   HWRITERV( 1, volume, tvolume)
   
  else
   
#ifdef DEBUG
   if( Whdfint ) write(ounit,'("hdfint : ", 10x ," : myid=",i3," ; vvolume is not allocated ;")') myid
#endif
   
  endif ! end of if( allocated(vvolume) ) ; 11 Aug 14;

!latex \item \verb+tflux(1:Nvol)+; real; toroidal flux enclosed by each interface; this is (probably) normalized so that \verb+tflux(Nvol)=1+;

  HWRITERV( Nvol, tflux, tflux(1:Nvol) )
  
!latex \item \verb+pflux(1:Nvol)+; real; poloidal flux enclosed by each interface; 18 Apr 17;

  HWRITERV( Nvol, pflux, pflux(1:Nvol) )
  
!latex \item \verb+pressure(1:Nvol)+; real; pressure in each volume;

  HWRITERV( Nvol, pressure, pressure(1:Nvol) )

!latex \item \verb+iota(1:Nvol)+; real; rotational-transform on each interface (inner);

  HWRITERV( Nvol, iota, iota(1:Nvol) )

!latex \item \verb+oita(1:Nvol)+; real; rotational-transform on each interface (outer);

  HWRITERV( Nvol, oita, oita(1:Nvol) )

!latex \item \verb+mu(1:Nvol)+; real; Lagrange multiplier, parallel current in each volume;

  HWRITERV( Nvol, mu, mu(1:Nvol) )

!latex \item \verb+helicity(1:Nvol)+; real; helicity in each volume;

  HWRITERV( Nvol, helicity, helicity(1:Nvol) )

!latex \item \verb+forcetol+; real; force-balance error across interfaces;

  HWRITERV( 1, forcetol, forcetol )

!latex \item \verb+ForceErr+; real; force-balance error across interfaces;

  HWRITERV( 1, ForceErr, ForceErr )

!latex \item \verb+pscale+; real; pressure scale factor;

  HWRITERV( 1, pscale, pscale )

!latex \item \verb+Bsubtemn(1:mn,0:1,1:Mvol)+; real; the even Fourier harmonics, \verb+j=1,mn+, of the covariant poloidal field, i.e. $[[B_{\t,j}]]$
!latex       evaluated on the inner and outer, \verb+i=0,1+, interface in each volume, \verb+l=1,Mvol+;
  HWRITERC( mn, 2, Mvol, Bsubtemn, Bsubtemn(1:mn,0:1,1:Mvol) )
  
!latex \item \verb+Bsubzemn(1:mn,0:1,1:Mvol)+; real; the even Fourier harmonics, \verb+j=1,mn+, of the covariant toroidal field, i.e. $[[B_{\z,j}]]$
!latex       evaluated on the inner and outer, \verb+i=0,1+, interface in each volume, \verb+l=1,Mvol+;
  HWRITERC( mn, 2, Mvol, Bsubzemn, Bsubzemn(1:mn,0:1,1:Mvol) )
  
!latex \item \verb+Bsubtomn(1:mn,0:1,1:Mvol)+; real; the odd  Fourier harmonics, \verb+j=1,mn+, of the covariant poloidal field, i.e. $[[B_{\t,j}]]$
!latex       evaluated on the inner and outer, \verb+i=0,1+, interface in each volume, \verb+l=1,Mvol+;
  HWRITERC( mn, 2, Mvol, Bsubtomn, Bsubtomn(1:mn,0:1,1:Mvol) )
  
!latex \item \verb+Bsubzomn(1:mn,0:1,1:Mvol)+; real; the odd  Fourier harmonics, \verb+j=1,mn+, of the covariant toroidal field, i.e. $[[B_{\z,j}]]$
!latex       evaluated on the inner and outer, \verb+i=0,1+, interface in each volume, \verb+l=1,Mvol+;
  HWRITERC( mn, 2, Mvol, Bsubzomn, Bsubzomn(1:mn,0:1,1:Mvol) )

!latex \item \verb+Lperturbed+; integer;

  HWRITEIV( 1, Lperturbed, Lperturbed )

  if( Lperturbed.eq.1 ) then

!latex \item \verb+dRbc(1:mn,0:Nvol)+; real; Fourier harmonics, $R_{m,n}$, of interfaces; linearly perturbed solution;

  HWRITERA( mn, Nvol+1, dRbc, dRbc(1:mn,0:Nvol) )

!latex \item \verb+dZbs(1:mn,0:Nvol)+; real; Fourier harmonics, $Z_{m,n}$, of interfaces; linearly perturbed solution;

  HWRITERA( mn, Nvol+1, dZbs, dZbs(1:mn,0:Nvol) )

!latex \item \verb+dRbs(1:mn,0:Nvol)+; real; Fourier harmonics, $R_{m,n}$, of interfaces; linearly perturbed solution;

  HWRITERA( mn, Nvol+1, dRbs, dRbs(1:mn,0:Nvol) )

!latex \item \verb+dZbc(1:mn,0:Nvol)+; real; Fourier harmonics, $Z_{m,n}$, of interfaces; linearly perturbed solution;

  HWRITERA( mn, Nvol+1, dZbc, dZbc(1:mn,0:Nvol) )

  endif

!latex \item \verb+dpp+; integer;

  HWRITEIV( 1, dpp, dpp )

!latex \item \verb+dpp+; integer;

  HWRITEIV( 1, dqq, dqq )

!latex \item \verb+lmns+; resolution of straight fieldline transformation;

  HWRITEIV( 1, lmns, lmns )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \end{itemize}

!latex \item Note that all quantities marked as real should be treated as double precision.
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  call h5fclose_f( file_id, hdfier ) ! terminate access;

  FATALMESS(hdfint, hdfier.ne.0, error calling h5fclose_f )

  call h5close_f( hdfier ) ! close Fortran interface;

  FATALMESS(hdfint, hdfier.ne.0, error calling h5close_f )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#endif
  
  RETURN(hdfint)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine hdfint

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
