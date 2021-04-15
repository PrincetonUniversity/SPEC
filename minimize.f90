!=======================================================================
subroutine steepest_descent(ndim,x,f,xtol,gtol,ftol,maxiter &
     ,iprint,iflag,func,grad)
  implicit none
  integer,intent(in):: ndim,maxiter,iprint
  integer,intent(inout):: iflag
  real(8),intent(in):: xtol,gtol,ftol
  real(8),intent(inout):: f,x(ndim)
!!$    real(8):: func,grad
  interface
    function func(n,x)
      integer,intent(in):: n
      real(8),intent(in):: x(n)
      real(8):: func
    end function func
    function grad(n,x)
      integer,intent(in):: n
      real(8),intent(in):: x(n)
      real(8):: grad(n)
    end function grad
  end interface

  integer:: iter
  real(8):: alpha,fp,gnorm
  real(8),external:: sprod
  real(8),save,allocatable:: g(:)

  if( .not. allocated(g) ) allocate(g(ndim))

  iter= 0
  f= func(ndim,x)
  g= grad(ndim,x)
  gnorm= sprod(ndim,g,g)
  g(1:ndim)= -g(1:ndim)/sqrt(gnorm)
  gnorm= gnorm/ndim
  if( iprint.eq.1 ) then
    write(6,'(a,i8,10es15.7)') ' iter,f,gnorm=' &
         ,iter,f,gnorm
  else if( iprint.ge.2 ) then
    write(6,'(a,i8,10es15.7)') ' iter,x,f,gnorm=' &
         ,iter,x(1:ndim),f,gnorm
  endif

  do iter=1,maxiter
    fp= f
!.....line minimization
    call quad_interpolate(ndim,x,g,f,xtol,gtol,ftol,alpha,iprint &
         ,iflag,func)
!.....if quad interpolation failed, perform golden section
    if( iflag/100.ne.0 ) then
      iflag= iflag -(iflag/100)*100
      call golden_section(ndim,x,g,f,xtol,gtol,ftol,alpha &
           ,iprint,iflag,func)
    endif
    if( iflag/100.ne.0 ) return
    x(1:ndim)= x(1:ndim) +alpha*g(1:ndim)
    f= func(ndim,x)
    g= grad(ndim,x)
    gnorm= sprod(ndim,g,g)
    g(1:ndim)= -g(1:ndim)/sqrt(gnorm)
    gnorm= gnorm/ndim
    if( iprint.eq.1 ) then
      write(6,'(a,i8,10es15.7)') ' iter,f,gnorm=' &
           ,iter,f,gnorm
    else if( iprint.ge.2 ) then
      write(6,'(a,i8,10es15.7)') ' iter,x,f,gnorm=' &
           ,iter,x(1:ndim),f,gnorm
    endif
!.....check convergence 
!!$    if( abs(alpha).lt.xtol ) then
!!$      print *,'>>> SD converged wrt xtol'
!!$      write(6,'(a,2es15.7)') '   alpha,xtol=',alpha,xtol
!!$      iflag= iflag +1
!!$      return
    if( gnorm.lt.gtol ) then
      print *,'>>> SD converged wrt gtol'
      write(6,'(a,2es15.7)') '   gnorm,gtol=',gnorm,gtol
      iflag= iflag +2
      return
    else if( abs(f-fp)/abs(fp).lt.ftol ) then
      print *,'>>> Sd converged wrt ftol'
      write(6,'(a,2es15.7)') '   f-fp,ftol=',abs(f-fp)/abs(fp),ftol
      iflag= iflag +3
      return
    endif
  enddo

  print *,'maxiter exceeded in steepest_descent'
  print *,'steepest_descent done.'
  iflag= iflag +10
  return
end subroutine steepest_descent
!=======================================================================
subroutine frcg(ndim,x,f,xtol,gtol,ftol,maxiter,iprint,iflag &
     ,func,grad)
!
!  Conjugate gradient minimization
!
  implicit none
  integer,intent(in):: ndim,maxiter,iprint
  integer,intent(inout):: iflag
  real(8),intent(in):: xtol,gtol,ftol
  real(8),intent(inout):: f,x(ndim)
  interface
    function func(n,x)
      integer,intent(in):: n
      real(8),intent(in):: x(n)
      real(8):: func
    end function func
    function grad(n,x)
      integer,intent(in):: n
      real(8),intent(in):: x(n)
      real(8):: grad(n)
    end function grad
  end interface

  integer:: iter
  real(8):: alpha,fp,gnorm,gnormp
  real(8),external:: sprod
  real(8),save,allocatable:: g(:),u(:)

  if( .not. allocated(g) ) allocate(g(ndim),u(ndim))

  iter= 0
  f= func(ndim,x)
  g= grad(ndim,x)
  gnorm= sprod(ndim,g,g)
!!$    g(1:ndim)= g(1:ndim)/sqrt(gnorm)
  gnorm= gnorm/ndim
  if( iprint.eq.1 ) then
    write(6,'(a,i8,10es15.7)') ' iter,f,gnorm=' &
         ,iter,f,gnorm
  else if( iprint.ge.2 ) then
    write(6,'(a,i8,10es15.7)') ' iter,x,f,gnorm=' &
         ,iter,x(1:ndim),f,gnorm
  endif
  u(1:ndim)= -g(1:ndim)

  do iter=1,maxiter
    fp= f
!.....line minimization
    call quad_interpolate(ndim,x,u,f,xtol,gtol,ftol,alpha,iprint &
         ,iflag,func)
!.....if quad interpolation failed, perform golden section
    if( iflag/100.ne.0 ) then
      iflag= iflag -(iflag/100)*100
      call golden_section(ndim,x,u,f,xtol,gtol,ftol,alpha &
           ,iprint,iflag,func)
    endif
    if( iflag/100.ne.0 ) return
    x(1:ndim)= x(1:ndim) +alpha*u(1:ndim)
    f= func(ndim,x)
    g= grad(ndim,x)
!.....store previous gnorm
    gnormp= gnorm
    gnorm= sprod(ndim,g,g)
!!$      g(1:ndim)= g(1:ndim)/sqrt(gnorm)
    gnorm= gnorm/ndim
    u(1:ndim)= -g(1:ndim) +gnorm/gnormp *u(1:ndim)
    if( iprint.eq.1 ) then
      write(6,'(a,i8,10es15.7)') ' iter,f,gnorm=' &
           ,iter,f,gnorm
    else if( iprint.ge.2 ) then
      write(6,'(a,i8,10es15.7)') ' iter,x,f,gnorm=' &
           ,iter,x(1:ndim),f,gnorm
    endif
!.....check convergence 
!!$    if( abs(alpha).lt.xtol ) then
!!$      print *,'>>> CG converged wrt xtol'
!!$      write(6,'(a,2es15.7)') '   alpha,xtol=',alpha,xtol
!!$      iflag= iflag +1
!!$      return
    if( gnorm.lt.gtol ) then
      print *,'>>> CG converged wrt gtol'
      write(6,'(a,2es15.7)') '   gnorm,gtol=',gnorm,gtol
      iflag= iflag +2
      return
    else if( abs(f-fp)/abs(fp).lt.ftol ) then
      print *,'>>> CG converged wrt ftol'
      write(6,'(a,2es15.7)') '   f-fp,ftol=',abs(f-fp)/abs(fp),ftol
      iflag= iflag +3
      return
    endif
  enddo

  print *,'*** maxiter exceeded ***'
  print *,'steepest_descent done.'
  iflag= iflag +10
  return
end subroutine frcg
!=======================================================================
subroutine bfgs(ndim,x0,f,xtol,gtol,ftol,maxiter &
     ,iprint,iflag,func,grad)
!
!  Broyden-Fletcher-Goldfarb-Shanno type of Quasi-Newton method.
!
  implicit none
  integer,intent(in):: ndim,maxiter,iprint
  integer,intent(inout):: iflag
  real(8),intent(in):: xtol,gtol,ftol
  real(8),intent(inout):: f,x0(ndim)
!!$    real(8):: func,grad
  interface
    function func(n,x)
      integer,intent(in):: n
      real(8),intent(in):: x(n)
      real(8):: func
    end function func
    function grad(n,x)
      integer,intent(in):: n
      real(8),intent(in):: x(n)
      real(8):: grad(n)
    end function grad
  end interface

  real(8),external:: sprod
  real(8),save,allocatable:: gg(:,:),x(:),u(:),v(:),y(:),gp(:) &
       ,ggy(:),ygg(:),aa(:,:),cc(:,:),g(:)
  real(8):: tmp1,tmp2,b,svy,svyi,fp,alpha,gnorm
  integer:: i,j,iter

  if( .not.allocated(gg) ) allocate(gg(ndim,ndim),x(ndim),u(ndim)&
       ,v(ndim),y(ndim),g(ndim),gp(ndim),ggy(ndim),ygg(ndim) &
       ,aa(ndim,ndim),cc(ndim,ndim))

!.....initial G = I
  gg(1:ndim,1:ndim)= 0d0
  do i=1,ndim
    gg(i,i)= 1d0
  enddo
  f= func(ndim,x0)
  g= grad(ndim,x0)
  gnorm= sprod(ndim,g,g)
!!$    g(1:ndim)= g(1:ndim)/sqrt(gnorm)
  x(1:ndim)= x0(1:ndim)

  iter= 0
  gnorm= gnorm/ndim
  if( iprint.eq.1 ) then
    write(6,'(a,i8,100es15.7)') ' iter,f,gnorm=',iter,f,gnorm
  else if( iprint.ge.2 ) then
    write(6,'(a,i8,100es15.7)') ' iter,x,f,gnorm=' &
         ,iter,x(1:ndim),f,gnorm
  endif

  do iter=1,maxiter
    u(1:ndim)= 0d0
    do i=1,ndim
      u(1:ndim)= u(1:ndim) -gg(1:ndim,i)*g(i)
    enddo
!!$      print *,' u =',u(1:2)
!.....store previous func and grad values
    fp= f
    gp(1:ndim)= g(1:ndim)
!.....line minimization
    call quad_interpolate(ndim,x,u,f,xtol,gtol,ftol,alpha &
         ,iprint,iflag,func)
!.....if quad interpolation failed, perform golden section
    if( iflag/100.ne.0 ) then
      iflag= iflag -(iflag/100)*100
      call golden_section(ndim,x,u,f,xtol,gtol,ftol,alpha &
           ,iprint,iflag,func)
    endif
    if( iflag/100.ne.0 ) then
      x0(1:ndim)= x(1:ndim)
      return
    endif
    x(1:ndim)= x(1:ndim) +alpha*u(1:ndim)
    g= grad(ndim,x)
    gnorm= sprod(ndim,g,g)
!!$      g(1:ndim)= g(1:ndim)/sqrt(gnorm)
    gnorm= gnorm/ndim
    if( iprint.eq.1 ) then
      write(6,'(a,i8,100es15.7)') ' iter,f,gnorm=',iter,f,gnorm
    else if( iprint.ge.2 ) then
      write(6,'(a,i8,100es15.7)') ' iter,x,f,gnorm=' &
           ,iter,x(1:ndim),f,gnorm
    endif
!.....check convergence 
!!$    if( abs(alpha).lt.xtol ) then
!!$      print *,'>>> BFGS converged wrt xtol'
!!$      write(6,'(a,2es15.7)') '   alpha,xtol=',alpha,xtol
!!$      x0(1:ndim)= x(1:ndim)
!!$      iflag= iflag +1
!!$      return
    if( gnorm.lt.gtol ) then
      print *,'>>> BFGS converged wrt gtol'
      write(6,'(a,2es15.7)') '   gnorm,gtol=',gnorm,gtol
      x0(1:ndim)= x(1:ndim)
      iflag= iflag +2
      return
    else if( abs(f-fp)/abs(fp).lt.ftol ) then
      print *,'>>> BFGS converged wrt ftol'
      write(6,'(a,2es15.7)') '   f-fp,ftol=',abs(f-fp)/abs(fp),ftol
      x0(1:ndim)= x(1:ndim)
      iflag= iflag +3
      return
    endif

    v(1:ndim)= alpha *u(1:ndim)
    y(1:ndim)= g(1:ndim) -gp(1:ndim)

!.....update G matrix, gg, according to BFGS
    svy= sprod(ndim,v,y)
    svyi= 1d0/svy
    do i=1,ndim
      tmp1= 0d0
      tmp2= 0d0
      do j=1,ndim
        aa(j,i)= v(j)*v(i) *svyi
        tmp1= tmp1 +gg(i,j)*y(j)
        tmp2= tmp2 +y(j)*gg(j,i)
      enddo
      ggy(i)= tmp1
      ygg(i)= tmp2
    enddo
    b= 1d0
    do i=1,ndim
      b=b +y(i)*ggy(i) *svyi
    enddo
    aa(1:ndim,1:ndim)= aa(1:ndim,1:ndim) *b
    cc(1:ndim,1:ndim)= 0d0
    do j=1,ndim
      do i=1,ndim
        cc(i,j)=cc(i,j) +(v(i)*ygg(j) +ggy(i)*v(j)) *svyi
      enddo
    enddo
    gg(1:ndim,1:ndim)=gg(1:ndim,1:ndim) +aa(1:ndim,1:ndim) &
         -cc(1:ndim,1:ndim)
  enddo

  print *,'maxiter exceeded in bfgs'
  iflag= iflag +10
  x0(1:ndim)= x(1:ndim)
  return
end subroutine bfgs
!=======================================================================
subroutine get_bracket(ndim,x0,d,a,b,c,fa,fb,fc,iprint,iflag,func)
  implicit none
  integer,intent(in):: ndim,iprint
  integer,intent(inout):: iflag
  real(8),intent(in):: x0(ndim),d(ndim)
  real(8),intent(inout):: a,b,fa,fb
  real(8),intent(out):: c,fc
  real(8):: func

  real(8),parameter:: RATIO = 1.61803398875d0
  real(8),parameter:: RATIOI= 1d0/RATIO
  real(8),parameter:: TINY= 1d-12
  real(8),parameter:: GLIMIT= 100d0
  real(8),parameter:: MAXITER= 50
  real(8):: dum,r,q,u,ulim,fu
  integer:: iter

  fa= func(ndim,x0+a*d)
  fb= func(ndim,x0+b*d)
  iter= 0
10 continue
  iter= iter +1
  if( iter.gt.MAXITER ) then
    print *,'[Error] iter.gt.MAXITER in get_bracket'
    print *,'  Search direction may not be a descent direction.'
    iflag= iflag +1000
    return
!!$      stop
  endif
  if( abs(b-a).lt.1d-12) then
    print *,'[Error] a and b is too close in get_bracket'
    print *,'  Search direction may not be a descent direction.'
    iflag= iflag +2000
    return
!!$      stop
  endif
  if( fa.lt.fb ) then
    c= a +RATIOI*(b-a)
    fc= func(ndim,x0+c*d)
    call exchange(c,b)
    call exchange(fc,fb)
    if( iprint.eq.3 ) then
      write(6,'(a,2(1x,3es12.4))') ' a,b,c,fa,fb,fc=',a,b,c,fa,fb,fc
    endif
    goto 10
  else
    c= a +RATIO*(b-a)
    fc= func(ndim,x0+c*d)
    if( fb.gt.fc ) then
      b= a +RATIO*(c-a)
      fb= func(ndim,x0+b*d)
      call exchange(b,c)
      call exchange(fb,fc)
      if( iprint.eq.3 ) then
        write(6,'(a,2(1x,3es12.4))') ' a,b,c,fa,fb,fc=',a,b,c,fa,fb,fc
      endif
      goto 10
    endif
!!$      write(6,'(a,2(1x,3es12.4))') ' a,b,c,fa,fb,fc=',a,b,c,fa,fb,fc
  endif

  return
end subroutine get_bracket
!=======================================================================
subroutine exchange(a,b)
  implicit none
  real(8),intent(inout):: a,b
  real(8):: tmp
  tmp= a
  a= b
  b= tmp
  return
end subroutine exchange
!=======================================================================
subroutine quad_interpolate(ndim,x0,g,f,xtol,gtol,ftol,a,iprint &
     ,iflag,func)
  implicit none
  integer,intent(in):: ndim,iprint
  integer,intent(inout):: iflag
  real(8),intent(in):: x0(ndim),xtol,gtol,ftol
  real(8),intent(out):: f,g(ndim),a

  real(8),parameter:: STP0    = 1d-4
  real(8),parameter:: STPMAX  = 1d-1
  real(8),parameter:: TINY    = 1d-15
  integer,parameter:: MAXITER = 100

  interface
    function func(n,x)
      integer,intent(in):: n
      real(8),intent(in):: x(n)
      real(8):: func
    end function func
  end interface

  integer:: iter,imin,imax,ix
  real(8):: r,q,fmin,fmax,dmin,dmax,d,xmin
  real(8),save,allocatable:: xi(:),fi(:)

  if( .not. allocated(xi) ) allocate(xi(4),fi(4))

  xi(1)= 0d0
  xi(2)= xi(1) +STP0
  call get_bracket(ndim,x0,g,xi(1),xi(2),xi(3),fi(1),fi(2),fi(3)&
       ,iprint,iflag,func)
  if( iflag/1000.ne.0 ) return
!!$    fi(1)= func(ndim,x0+xi(1)*g)
!!$    fi(2)= func(ndim,x0+xi(2)*g)
!!$    if( fi(1).gt.fi(2) ) then
!!$      xi(3)= xi(1) +2*STP0
!!$      fi(3)= func(ndim,x0+xi(3)*g)
!!$    else
!!$      xi(3)= xi(1) -STP0
!!$      fi(3)= func(ndim,x0+xi(3)*g)
!!$    endif

  iter= 0
10 continue
  iter= iter +1
  if( iter.gt.MAXITER ) then
    print *,' [Error] iter.gt.MAXITER in quad_interpolate !!!'
    print *,'   iter,MAXITER= ',iter,MAXITER
    iflag= iflag +100
    return
!!$      stop
  endif
!.....step 3; compute turning point
  r= (xi(2)-xi(1))*(fi(2)-fi(3))
  q= (xi(2)-xi(3))*(fi(2)-fi(1))
  xi(4)= xi(2) -((xi(2)-xi(3))*q -(xi(2)-xi(1))*r) &
       /(2d0*sign(max(abs(q-r),TINY),q-r))
  fi(4)= func(ndim,x0+xi(4)*g)
!!$    write(6,'(a,2(2x,4f11.2))') ' xi,fi=',xi(1:4),fi(1:4)

!.....step4
  fmin= min(fi(1),fi(2),fi(3))
  fmax= max(fi(1),fi(2),fi(3))
  dmin= min(abs(xi(4)-xi(1)),abs(xi(4)-xi(2)),abs(xi(4)-xi(3)))
  if( fi(4).lt.fmin .and. dmin.gt.STPMAX ) then
!!$      print *,' 01'
    imax= 0
    dmax= 0d0
    do ix=1,3
      d= abs(xi(4)-xi(ix))
      if( dmax.lt.d ) then
        dmax= d
        imax= ix
      endif
    enddo
!.....eliminate max function point and set point at STPMAX
    imax= 0
    fmax= -1d+30
    do ix=1,3
      if( fmax.lt.fi(ix) ) then
        fmax= fi(ix)
        imax= ix
      endif
    enddo
    do ix=imax+1,3
      xi(ix-1)= xi(ix)
      fi(ix-1)= fi(ix)
    enddo
    if( fi(2).gt.fi(1) ) then
      xi(3)= xi(1) +STPMAX
    else
      xi(3)= xi(2) +STPMAX
    endif
    fi(3)= func(ndim,x0+xi(3)*g)
    goto 10
  else if( fi(4).gt.fmax ) then ! fi(4) is maximum
!!$      print *,' 02'
    imin= 0
    dmin= 1d+30
    do ix=1,3
      d= abs(xi(4)-xi(ix))
      if( dmin.gt.d ) then
        dmin= d
        xmin= xi(ix)
        imin= ix
      endif
    enddo
!.....eliminate nearest point to xi(4) and add a new point
    do ix=imin+1,3
      xi(ix-1)= xi(ix)
      fi(ix-1)= fi(ix)
    enddo
!!$      if( fi(2).gt.fi(1) ) then
!!$        xi(3)= xi(1) +STPMAX
!!$      else
!!$        xi(3)= xi(2) +STPMAX
!!$      endif
    xi(3)= (xmin +xi(4))*0.5
    fi(3)= func(ndim,x0+xi(3)*g)
    goto 10
  endif

!.....step 5: check convergence
  if( dmin.lt.xtol ) then
    imin= 0
    fmin= 1d+30
    do ix=1,4
      if( fmin.gt.fi(ix) ) then
        imin= ix
        fmin= fi(ix)
      endif
    enddo
    f= fi(imin)
    a= xi(imin)
    return
  endif

!.....step 6: discard point of highest f value and replace it by xi(4)
!!$    print *,' 03'
  imax= 0
  fmax= -1d+30
  do ix=1,3
    if( fmax.lt.fi(ix) ) then
      imax= ix
      fmax= fi(ix)
    endif
  enddo
  xi(imax)= xi(4)
  fi(imax)= fi(4)
  goto 10

end subroutine quad_interpolate
!=======================================================================
subroutine golden_section(ndim,x0,g,f,xtol,gtol,ftol,alpha,iprint &
     ,iflag,func)
  implicit none
  integer,intent(in):: ndim,iprint
  integer,intent(inout):: iflag
  real(8),intent(in):: xtol,gtol,ftol,x0(ndim),g(ndim)
  real(8),intent(inout):: f,alpha
  real(8):: func

  real(8),parameter:: STP0 = 1d-1
  real(8),parameter:: GR   = 0.61803398875d0
  real(8),parameter:: GR2  = 1d0 -GR
  integer,parameter:: MAXITER= 100

  integer:: iter
  real(8):: a,b1,b2,c,fa,fb1,fb2,fc,xl

  a= 0d0
  b1= STP0
  call get_bracket(ndim,x0,g,a,b1,c,fa,fb1,fc,iprint,iflag,func)
  if( iflag/1000.ne.0 ) return
  xl= (c-a)
  b1= a +GR2*xl
  b2= a +GR *xl
  fb1= func(ndim,x0+b1*g)
  fb2= func(ndim,x0+b2*g)

  iter= 0
10 continue
  iter= iter +1
  if( iter.gt.MAXITER ) then
    print *,'[Error] iter.gt.MAXITER in golden_section.'
    print *,'  iter,MAXITER = ',iter,MAXITER
    iflag= iflag +100
    return
!!$      stop
  endif
!!$    write(6,'(a,2(2x,4es11.3))') ' a,b1,b2,c,fa,fb1,fb2,fc=' &
!!$         ,a,b1,b2,c,fa,fb1,fb2,fc
  if( fb1.gt.fb2 ) then
    a= b1
    fa= fb1
    b1= b2
    fb1= fb2
    xl= (c-a)
    b2= a +GR*xl
    fb2= func(ndim,x0+b2*g)
  else
    c= b2
    fc= fb2
    b2= b1
    fb2= fb1
    xl= (c-a)
    b1= a +GR2*xl
    fb1= func(ndim,x0+b1*g)
  endif
!!$    print *,' xl,c,a,xtol=',xl,c,a,xtol
  if( xl.lt.xtol ) then
    if( fb1.gt.fb2 ) then
      alpha= b2
      f= fb2
    else
      alpha= b1
      f= fb1
    endif
    return
  endif
  goto 10

end subroutine golden_section
!=======================================================================
function sprod(n,a,b)
  implicit none
  integer,intent(in):: n
  real(8),intent(in):: a(n),b(n)
  real(8):: sprod

  integer:: i
  sprod= 0d0
  do i=1,n
    sprod= sprod +a(i)*b(i)
  enddo
  return
end function sprod
