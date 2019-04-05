subroutine dfp100(vvol, Fdof, Ddof)

use constants, only : zero, half, one, two, pi2, pi

use fileunits, only : ounit

use inputlist, only : Wmacros, Wdfp100, Igeometry, Nvol, Lrad, Isurf, &
					  Lconstraint

use cputiming, only : Tdfp100

use allglobal, only : ncpu, myid, cpus, &
                      ImagneticOK, NAdof, mn, &
                      Mvol, &
                      dBdX, Iquad, &
                      Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, Localconstraint, &
                      DToocc, DToocs, DToosc, DTooss, &
                      TTsscc, TTsscs, TTsssc, TTssss, &
                      TDstcc, TDstcs, TDstsc, TDstss, &
                      TDszcc, TDszcs, TDszsc, TDszss, &
                      DDttcc, DDttcs, DDttsc, DDttss, &
                      DDtzcc, DDtzcs, DDtzsc, DDtzss, &
                      DDzzcc, DDzzcs, DDzzsc, DDzzss, &
                      IPDt


 LOCALS
!------

INTEGER              :: vvol, NN, ll, ii

REAL                 :: Fdof(1:Mvol-1), Ddof(1:Mvol-1)



  BEGIN(dfp100)


   LREGION(vvol) ! assigns Lcoordinatesingularity, Lplasmaregion, etc. ;

   ImagneticOK(vvol) = .false.
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   if( myid.ne.modulo(vvol-1,ncpu) ) goto 5000 ! construct Beltrami fields in parallel;
   
   NN = NAdof(vvol) ! shorthand;

   ll = Lrad(vvol)

   SALLOCATE( DToocc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( DToocs, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( DToosc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( DTooss, (0:ll,0:ll,1:mn,1:mn), zero )

   SALLOCATE( TTsscc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( TTsscs, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( TTsssc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( TTssss, (0:ll,0:ll,1:mn,1:mn), zero )

   SALLOCATE( TDstcc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( TDstcs, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( TDstsc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( TDstss, (0:ll,0:ll,1:mn,1:mn), zero )

   SALLOCATE( TDszcc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( TDszcs, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( TDszsc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( TDszss, (0:ll,0:ll,1:mn,1:mn), zero )

   SALLOCATE( DDttcc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( DDttcs, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( DDttsc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( DDttss, (0:ll,0:ll,1:mn,1:mn), zero )

   SALLOCATE( DDtzcc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( DDtzcs, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( DDtzsc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( DDtzss, (0:ll,0:ll,1:mn,1:mn), zero )

   SALLOCATE( DDzzcc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( DDzzcs, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( DDzzsc, (0:ll,0:ll,1:mn,1:mn), zero )
   SALLOCATE( DDzzss, (0:ll,0:ll,1:mn,1:mn), zero )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   dBdX%L = .false. ! first, compute Beltrami fields;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   write( ounit, '("dforce, vvol = ", i3, ". starting ma00aa ...")') vvol
   WCALL( dfp100, ma00aa, ( Iquad(vvol), mn, vvol, ll ) ) ! compute volume integrals of metric elements - evaluate TD, DT, DD, ...;
   write( ounit, '("dforce, ma00aa done. starting matrix ... ")')
   WCALL( dfp100, matrix, ( vvol, mn, ll ) )
   write( ounit, '("dforce, matrix done. starting ma02aa ... ")')
   WCALL( dfp100, ma02aa, ( vvol, NN ) )
   write( ounit, '("dforce, ma02aa done.")')

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   DALLOCATE(DToocc)
   DALLOCATE(DToocs)
   DALLOCATE(DToosc)
   DALLOCATE(DTooss)

   DALLOCATE(TTsscc)
   DALLOCATE(TTsscs)
   DALLOCATE(TTsssc)
   DALLOCATE(TTssss)

   DALLOCATE(TDstcc)
   DALLOCATE(TDstcs)
   DALLOCATE(TDstsc)
   DALLOCATE(TDstss)

   DALLOCATE(TDszcc)
   DALLOCATE(TDszcs)
   DALLOCATE(TDszsc)
   DALLOCATE(TDszss)

   DALLOCATE(DDttcc)
   DALLOCATE(DDttcs)
   DALLOCATE(DDttsc)
   DALLOCATE(DDttss)

   DALLOCATE(DDtzcc)
   DALLOCATE(DDtzcs)
   DALLOCATE(DDtzsc)
   DALLOCATE(DDtzss)

   DALLOCATE(DDzzcc)
   DALLOCATE(DDzzcs)
   DALLOCATE(DDzzsc)
   DALLOCATE(DDzzss)

   if( .not.LocalConstraint ) then

     select case (Lconstraint)

	 case( 3 )
		! Compute IPDt on each interface. Eventually need to put the do
		! loop inside the surfcurent subroutine... TODO
		do ii = 1, Mvol-1
			WCALL(dfp100, surfcurent, (ii, mn))
		enddo
		
		Fdof = IPDt - Isurf
		!Ddof = ??? TODO: SEE IF AN ANALYTICAL FORMULATION EXISTS...

	 case default
		FATAL(dfp100, .true., Unaccepted value for Lconstraint)

	 end select
   endif


5000 continue


RETURN(dfp100)

end subroutine dfp100
