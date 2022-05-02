diff --git a/bndRep.f90 b/bndRep.f90
index 0850c02b..3fe7c657 100644
--- a/bndRep.f90
+++ b/bndRep.f90
@@ -36,7 +36,7 @@ module bndRep
       ! In this subroutine we compute the mapping matrix, and allocate necessary memory
       ! This should only be called once at the beginning of preset.
 
-      use constants, only: zero, half
+      use constants, only: zero, half, one
       use inputlist, only: Wmacros, Mpol, Ntor, twoalpha, Nfp, Lboundary, tflux
       use fileunits, only: ounit, lunit
       use allglobal, only: myid, mn_field, im_field, in_field, &
@@ -130,16 +130,16 @@ module bndRep
 
       call build_mapping_matrices()
 
-      SALLOCATE( precond_rho, (1:mn_rho, 1:Mvol), zero )
+      SALLOCATE( precond_rho, (1:mn_rho, 1:Mvol), one )
       do vvol = 1, Mvol
         do ii = 1, mn_rho
-          precond_rho( ii, vvol ) = Rscale * tflux(vvol)**(im_rho(ii) * half) 
+          !precond_rho( ii, vvol ) = Rscale * tflux(vvol)**(im_rho(ii) * half) 
         enddo
       enddo
 
-      SALLOCATE( precond_b, (1:Mvol), zero )
+      SALLOCATE( precond_b, (1:Mvol), one )
       do vvol = 1, Mvol
-        precond_b( vvol ) = Rscale * tflux(vvol)**half
+        !precond_b( vvol ) = Rscale * tflux(vvol)**half
       enddo
 
     end subroutine initialize_mapping
diff --git a/dforce.f90 b/dforce.f90
index d2c36969..0f1cefb4 100644
--- a/dforce.f90
+++ b/dforce.f90
@@ -671,8 +671,30 @@ endif
 
     ! Evaluate force gradient
 #ifdef DEBUG
-    if( Lcheck.eq.6 ) then
-       WCALL(dforce, fndiff_dforce, ( NGdof_force, NGdof_field ) )
+    if( Lcheck.eq.6 .or. Lcheck.eq.7 ) then
+       WCALL(dforce, fndiff_dforce, ( NGdof_force, NGdof_field, dforcedRZ ) )
+    endif
+
+    if( Lcheck.eq.8 ) then
+      open(10, file=trim(ext)//'.dforcedRZ.txt', status='unknown')
+      write(ounit,'(A)') NEW_LINE('A')
+
+      do ii=1, NGdof_force
+        write(10   ,5438) dforcedRZ(ii,:)
+      enddo
+      close(10)
+
+      open(10, file=trim(ext)//'.dRZdhenn.txt', status='unknown')
+      write(ounit,'(A)') NEW_LINE('A')
+
+      do ii=1, NGdof_field
+        write(10   ,5438) dRZdhenn(ii,:)
+      enddo
+      close(10)
+
+      FATAL( dforce, .true., Lcheck.eq.8 completed. )
+
+5438 format(512F22.16, " ")
     endif
 #endif
 
@@ -696,7 +718,7 @@ end subroutine dforce
 !>
 !> @param[in] NGdof_force
 !> @param[in] NGdof_field
-subroutine fndiff_dforce( NGdof_force, NGdof_field )
+subroutine fndiff_dforce( NGdof_force, NGdof_field, dforcedRZ )
 
 use constants, only: zero, one, half, two
 
@@ -714,7 +736,7 @@ use allglobal, only: ncpu, myid, cpus, MPI_COMM_SPEC, &
                      mn_field, im_field, in_field, &
                      mn_rho, im_rho, in_rho, &
                      iRbc, iZbs, iRbs, iZbc, &
-                     LGdof_force, LGdof_field, psifactor, dBdX, &
+                     LGdof_force, LGdof_field, LGdof_bnd, psifactor, dBdX, &
                      YESstellsym, NOTstellsym, &
                      hessian, ext, &
                      irhoc, iR0c, iZ0s, ibc, Rscale
@@ -726,6 +748,7 @@ LOCALS
 
 INTEGER, intent(in) :: NGdof_force !< first dimension of the force gradient
 INTEGER, intent(in) :: NGdof_field !< second dimension of the force gradient
+REAL,    intent(in) :: dforcedRZ(1:NGdof_force, 1:NGdof_field)
 
 INTEGER             :: vvol !< Loop index on perturbed interfaces
 INTEGER             :: lvol !< Loop index on interfaces for forward mapping 
@@ -759,10 +782,9 @@ BEGIN(dforce)
 
   if( ncpu.eq.1) then
 
-  SALLOCATE( finitediff_estimate, (1:NGdof_force, 1:NGdof_force), zero )
 
   dBdX%L = .false.
-  if( Lboundary.eq.0 ) then
+  if( Lcheck.eq.6 ) then
     SALLOCATE( oRbc, (1:mn_field,0:Mvol), iRbc(1:mn_field,0:Mvol) ) !save unperturbed geometry
     SALLOCATE( oZbs, (1:mn_field,0:Mvol), iZbs(1:mn_field,0:Mvol) )
     SALLOCATE( oRbs, (1:mn_field,0:Mvol), iRbs(1:mn_field,0:Mvol) )
@@ -770,6 +792,7 @@ BEGIN(dforce)
     SALLOCATE( iforce,  (-2:2, 0:NGdof_force), zero)
     SALLOCATE( iposition, (-2:2, 0:NGdof_field), zero)
 
+    SALLOCATE( finitediff_estimate, (1:NGdof_force, 1:NGdof_field), zero )
 
     do vvol = 1, Mvol-1 ! loop over interior surfaces;
       idof = 0
@@ -812,8 +835,8 @@ BEGIN(dforce)
               endif
 
               packorunpack = 'P' ! pack geometrical degrees-of-freedom;
-              !LComputeAxis = .false. ! keep axis fixed
-              LComputeAxis = .true.
+              LComputeAxis = .false. ! keep axis fixed
+              !LComputeAxis = .true.
 
               WCALL(dforce, packxi,( NGdof_field, iposition(isymdiff,0:NGdof_field), Mvol, mn_field, iRbc(1:mn_field,0:Mvol),&
                                     iZbs(1:mn_field,0:Mvol), iRbs(1:mn_field,0:Mvol), &
@@ -842,7 +865,44 @@ BEGIN(dforce)
     DALLOCATE(oZbs)
     DALLOCATE(oRbc)
 
-  else ! Lboundary.eq.1
+    ! Print in file for diagnostics
+    if(myid.eq.0) then
+      ! Print hessian
+      open(10, file=trim(ext)//'.Lcheck6_output.txt', status='unknown')
+      write(ounit,'(A)') NEW_LINE('A')
+
+      do ii=1, NGdof_force
+        write(ounit,1345) myid, dforcedRZ(ii,:)
+        write(10   ,1347) dforcedRZ(ii,:)
+      enddo
+      close(10)
+
+      write(ounit,'(A)') NEW_LINE('A')
+
+      ! Print finite differences
+      open(10, file=trim(ext)//'.Lcheck6_output.FiniteDiff.txt', status='unknown')
+      do ii=1, NGdof_force
+        write(ounit,1346) myid, finitediff_estimate(ii,:)
+        write(10   ,1347) finitediff_estimate(ii,:)
+      enddo
+      write(ounit,'(A)') NEW_LINE('A')
+      close(10)
+
+1345 format("dforce: myid=",i3,"  ; Hessian            = ",512f16.10 "   ;")
+1346 format("dforce: myid=",i3,"  ; Finite differences = ",512f16.10 "   ;")
+1347 format(512F22.16, " ")
+
+    endif ! myid.eq.0
+
+
+
+  elseif( Lcheck.eq.7 ) then
+
+    if( Lboundary.eq.0 ) then
+      FATAL( dforce, .true., Lcheck.eq.7 only compatible with Lboundary.eq.1 )
+    endif
+
+    SALLOCATE( finitediff_estimate, (1:NGdof_force, 1:NGdof_force), zero )
 
     ! First save the initial geometry
     SALLOCATE( orhoc, (1:mn_rho, 1:Mvol), irhoc(1:mn_rho, 1:Mvol) )
@@ -895,7 +955,7 @@ BEGIN(dforce)
                                       - 8 * iforce(-1,0:NGdof_force) &
                                       + 1 * iforce(-2,0:NGdof_force))  / ( 12 * dRZ )
 
-        tdof = (vvol-1) * LGdof_field + idof
+        tdof = (vvol-1) * LGdof_bnd + idof
 
         finitediff_estimate(1:NGdof_force, tdof) = iforce(0, 1:NGdof_force) * precond_rho( ii, vvol )
       enddo ! ii
@@ -932,8 +992,14 @@ BEGIN(dforce)
   
             ! Pack to dofs
             packorunpack = 'P' ! pack geometrical degrees-of-freedom;
-            LComputeAxis = .false. ! keep axis fixed
-  
+
+            ! if( irz.ne.0 .and. vvol.eq.1 ) then
+            !   LComputeAxis = .true.
+            ! else
+            !   LComputeAxis = .false. ! keep axis fixed
+            ! endif
+            LComputeAxis = .false.
+
             WCALL(dforce, packxi,( NGdof_field, iposition(isymdiff,0:NGdof_field), Mvol, mn_field, iRbc(1:mn_field,0:Mvol),&
                                   iZbs(1:mn_field,0:Mvol), iRbs(1:mn_field,0:Mvol), iZbc(1:mn_field,0:Mvol), packorunpack, .false., LComputeAxis ) )
             WCALL(dforce, dforce,( NGdof_field, iposition(isymdiff,0:NGdof_field), iforce(isymdiff,0:NGdof_force), .false., LComputeAxis) )
@@ -946,7 +1012,7 @@ BEGIN(dforce)
                                         - 8 * iforce(-1,0:NGdof_force) &
                                         + 1 * iforce(-2,0:NGdof_force))  / ( 12 * dRZ )
   
-          tdof = (vvol-1) * LGdof_field + idof
+          tdof = (vvol-1) * LGdof_bnd + idof
 
           select case( irz )
           case( 0 ) ! bn modes
@@ -969,40 +1035,39 @@ BEGIN(dforce)
     DALLOCATE( iposition )
     DALLOCATE( iforce )
 
-  endif
-
+    ! Print in file for diagnostics
+    if(myid.eq.0) then
+      ! Print hessian
+      open(10, file=trim(ext)//'.Lcheck6_output.txt', status='unknown')
+      write(ounit,'(A)') NEW_LINE('A')
 
-  ! Print in file for diagnostics
-  if(myid.eq.0) then
-    ! Print hessian
-    open(10, file=trim(ext)//'.Lcheck6_output.txt', status='unknown')
-    write(ounit,'(A)') NEW_LINE('A')
+      do ii=1, mn_force
+        write(ounit,2345) myid, im_force(ii), in_force(ii), hessian(ii,:)
+        write(10   ,2347) hessian(ii,:)
+      enddo
+      close(10)
 
-    do ii=1, mn_force
-      write(ounit,1345) myid, im_force(ii), in_force(ii), hessian(ii,:)
-      write(10   ,1347) hessian(ii,:)
-    enddo
-    close(10)
+      write(ounit,'(A)') NEW_LINE('A')
 
-    write(ounit,'(A)') NEW_LINE('A')
+      ! Print finite differences
+      open(10, file=trim(ext)//'.Lcheck6_output.FiniteDiff.txt', status='unknown')
+      do ii=1, mn_force
+        write(ounit,2346) myid, im_force(ii), in_force(ii), finitediff_estimate(ii,:)
+        write(10   ,2347) finitediff_estimate(ii,:)
+      enddo
+      write(ounit,'(A)') NEW_LINE('A')
+      close(10)
 
-    ! Print finite differences
-    open(10, file=trim(ext)//'.Lcheck6_output.FiniteDiff.txt', status='unknown')
-    do ii=1, mn_force
-      write(ounit,1346) myid, im_force(ii), in_force(ii), finitediff_estimate(ii,:)
-      write(10   ,1347) finitediff_estimate(ii,:)
-    enddo
-    write(ounit,'(A)') NEW_LINE('A')
-    close(10)
+2345 format("dforce: myid=",i3," ; (",i4,",",i4," ; Hessian            = ",512f16.10 "   ;")
+2346 format("dforce: myid=",i3," ; (",i4,",",i4," ; Finite differences = ",512f16.10 "   ;")
+2347 format(512F22.16, " ")
 
-    1345 format("dforce: myid=",i3," ; (",i4,",",i4," ; Hessian            = ",512f16.10 "   ;")
-    1346 format("dforce: myid=",i3," ; (",i4,",",i4," ; Finite differences = ",512f16.10 "   ;")
-    1347 format(512F22.16, " ")
+    endif ! myid.eq.0
 
-  endif
+  endif ! Lcheck.eq.6
 
   DALLOCATE(finitediff_estimate)
-  endif
+  endif !ncpu.eq.1
 
 
 FATAL(fndiff, .true., Finite differences have been evaluated. )
diff --git a/dfp200.f90 b/dfp200.f90
index 4a417697..35be77e5 100644
--- a/dfp200.f90
+++ b/dfp200.f90
@@ -1114,6 +1114,7 @@ subroutine evaluate_dmupfdx(innout, idof, ii, issym, irz)
         dBdX%L = .false.
 
         if( LocalConstraint ) then
+            Lcomputederivatives = .False.
             WCALL(dfp200, deallocate_geometry_matrices, (LcomputeDerivatives))
             WCALL(dfp200, deallocate_Beltrami_matrices, (LcomputeDerivatives))
             WCALL(dfp200, intghs_workspace_destroy, ()))
@@ -1127,6 +1128,7 @@ subroutine evaluate_dmupfdx(innout, idof, ii, issym, irz)
             iZbs(1:mn_field,0:Mvol) = oZbs(1:mn_field,0:Mvol)
             iRbs(1:mn_field,0:Mvol) = oRbs(1:mn_field,0:Mvol)
             iZbc(1:mn_field,0:Mvol) = oZbc(1:mn_field,0:Mvol)
+            
 
             if( LocalConstraint ) then
                 if( issym.eq.0 .and. irz.eq.0 ) iRbc(ii,vvol-1+innout) = iRbc(ii,vvol-1+innout) + dRZ * isymdiff ! perturb geometry;
@@ -1139,6 +1141,10 @@ subroutine evaluate_dmupfdx(innout, idof, ii, issym, irz)
                 if( issym.eq.1 .and. irz.eq.0 ) iRbs(ii,vvol) = iRbs(ii,vvol) + dRZ * isymdiff ! perturb geometry;
                 if( issym.eq.1 .and. irz.eq.1 ) iZbc(ii,vvol) = iZbc(ii,vvol) + dRZ * isymdiff ! perturb geometry;
             endif
+            
+            if( im_field(ii).eq.0 ) then
+                WCALL( dfp200, rzaxis, ( Mvol, mn_field, iRbc(1:mn_field,0:Mvol), iZbs(1:mn_field,0:Mvol), iRbs(1:mn_field,0:Mvol), iZbc(1:mn_field,0:Mvol), vvol, LComputeDerivatives ) ) ! set coordinate axis; 19 Jul 16;            endif
+            endif
 
             ! Solve Beltrami equation consistently with constraints
             Xdof(1:Mvol-1) = zero;
@@ -1221,6 +1227,7 @@ subroutine evaluate_dmupfdx(innout, idof, ii, issym, irz)
 
         if( LocalConstraint ) then
             ! reallocate matrices for next iteration
+            LcomputeDerivatives = .FALSE.
             WCALL(dfp200, intghs_workspace_init, (vvol))
             WCALL(dfp200, allocate_Beltrami_matrices, (vvol,LcomputeDerivatives))
             WCALL(dfp200, allocate_geometry_matrices, (vvol,LcomputeDerivatives))
diff --git a/preset.f90 b/preset.f90
index bd0406be..c253e556 100644
--- a/preset.f90
+++ b/preset.f90
@@ -94,6 +94,8 @@ subroutine set_global_variables()
     Rscale = Rwc(0,0)
   endif
 
+  Rscale = one
+
   ! set up Henneberg's mapping
   if( Lboundary.eq.0 ) then
     twoalpha = 0.0
@@ -1415,7 +1417,8 @@ subroutine set_global_variables()
 !>       \end{array}\right.
 !>       \f}
 !>       where \f$\psi_{t,v} \equiv\,\f$\c tflux is the (normalized?) toroidal flux enclosed by the \f$v\f$-th interface. </li>
-!> <li> \c psifactor is used in packxi(), dforce() and hesian(). </li>
+!> <li> If \c Lboundary equals 1 (Henneberg's representation), then psifactor=1.
+!> <li> \c psifactor is used in packxi(), dforce(), dfp200() and hesian(). </li>
 !> <li> \c inifactor is similarly constructed, with
 !>       \f{eqnarray}{ f_{j,v} \equiv \left\{
 !>       \begin{array}{lcccccc}\psi_{t,v}^{ 1 /2}&,&\mbox{for $m_j=0$}, \\
@@ -1445,21 +1448,35 @@ subroutine set_global_variables()
   enddo
 
   case( 3 )
-  do vvol = 1, Nvol
-    do ii = 1, mn_field
-    if( im_field(ii).eq.0 ) then ; psifactor(ii,vvol) = Rscale * tflux(vvol)**zero       ! 08 Feb 16;
-                                 ; inifactor(ii,vvol) = Rscale * tflux(vvol)**half       ! 17 Dec 18;
-    else                         ; psifactor(ii,vvol) = Rscale * tflux(vvol)**halfmm(ii) ! 29 Apr 15;
-                                 ; inifactor(ii,vvol) = Rscale * tflux(vvol)**halfmm(ii) ! 17 Dec 18
+    if( Lboundary.eq.0 ) then ! Otherwise leave inifactor = 1
+      do vvol = 1, Nvol
+        do ii = 1, mn_field
+        if( im_field(ii).eq.0 ) then
+          psifactor(ii,vvol) = Rscale * tflux(vvol)**zero       ! 08 Feb 16;
+        else
+          psifactor(ii,vvol) = Rscale * tflux(vvol)**halfmm(ii) ! 29 Apr 15;
+        endif
+        enddo
+      enddo
     endif
-    enddo
-  enddo
 
   case default
   FATAL( readin, .true., invalid Igeometry for construction of psifactor )
 
   end select
 
+  if( Igeometry.eq.3 ) then 
+    do vvol = 1, Nvol
+      do ii = 1, mn_field
+      if( im_field(ii).eq.0 ) then
+        inifactor(ii,vvol) = Rscale * tflux(vvol)**half       ! 17 Dec 18;
+      else
+        inifactor(ii,vvol) = Rscale * tflux(vvol)**halfmm(ii) ! 17 Dec 18
+      endif
+      enddo
+    enddo
+  endif
+
   !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
 
   if( Lconstraint .EQ. 3) then
