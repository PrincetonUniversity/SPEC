!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (numerics) ! Some miscellaneous numerical routines.

!latex \briefly{miscellaneous ``numerical'' routines}

!l tex \calledby{\link{}}
!l tex \calls{\link{}}

!latex \tableofcontents

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{Outline}

!latex This file contains various miscellaneous ``numerical'' routines as described below.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!l tex \begin{itemize}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!l tex \item \type{gi00aa}

!subroutine gi00aa( ii, jj, ig ) ! not used; SRH: 27 Feb 18;
!  
!  implicit none
!  
!  INTEGER, intent(in)  :: ii,jj
!  INTEGER, intent(out) :: ig
!  
!  if( ( ii.eq.1 .and. jj.eq.1 )                                ) ig = 1
!  if( ( ii.eq.1 .and. jj.eq.2 ) .or. ( ii.eq.2 .and. jj.eq.1 ) ) ig = 2
!  if( ( ii.eq.1 .and. jj.eq.3 ) .or. ( ii.eq.3 .and. jj.eq.1 ) ) ig = 3
!  if( ( ii.eq.2 .and. jj.eq.2 )                                ) ig = 4
!  if( ( ii.eq.2 .and. jj.eq.3 ) .or. ( ii.eq.3 .and. jj.eq.2 ) ) ig = 5
!  if( ( ii.eq.3 .and. jj.eq.3 )                                ) ig = 6
!  
!  return
!  
!end subroutine gi00aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{\type{gi00ab}}

!latex \begin{enumerate}

!latex \item This routine assigns the Fourier mode labels that converts a double-sum into a single sum; i.e., the $m_j$ and $n_j$ are assigned where
!latex \be f(\t,\z) & = & \sum_{n=0}^{N} f_{0,n}\cos(-n \, N_P \, \z) 
!latex + \sum_{m=1}^{M} \sum_{n=-N}^{N} f_{m,n}\cos(m\t-n \, N_P \, \z) \\
!latex              & = & \sum_j f_j \cos(m_j\t-n_j\z),
!latex \ee
!latex where $N\equiv $ \type{Ntor} and $M\equiv $ \type{Mpol} are given on input, and $N_P \equiv $ \type{Nfp} is the field periodicity.

!latex \end{enumerate}

subroutine gi00ab( Mpol, Ntor, Nfp, mn, im, in )
  
  implicit none
  
  INTEGER, intent(in)  :: Mpol, Ntor, Nfp, mn
  INTEGER, intent(out) :: im(mn), in(mn)
  
  INTEGER              :: imn, mm, nn
  
  imn = 0
  
  ;  mm = 0  
  ;do nn = 0, Ntor
  ; imn = imn+1 ; im(imn) = mm ; in(imn) = nn*Nfp
  ;enddo
  ;

  do mm = 1, Mpol
   do nn = -Ntor, Ntor
    imn = imn+1 ; im(imn) = mm ; in(imn) = nn*Nfp
   enddo
  enddo
  
  return
  
end subroutine gi00ab

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{\type{tfft}}

!latex \begin{enumerate}

!latex \item This constructs the ``forward'' Fourier transform. Given a set of data $(f_{i,j},g_{i,j})$ on a regular grid angle grid, 
!latex the Fourier harmonics are constructed.

!latex \end{enumerate}

subroutine tfft( Nt, Nz, ijreal, ijimag, mn, im, in, efmn, ofmn, cfmn, sfmn, ifail )

  use constants

  use inputlist, only : Nfp
  use allglobal, only : pi2nfp

  use fftw_interface

  implicit none

  intrinsic aimag
  
  INTEGER :: Nt, Nz, mn, im(mn), in(mn), Ntz, imn, ifail, mm, nn
  REAL    :: ijreal(1:Nt*Nz), ijimag(1:Nt*Nz), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn)
  
  LOGICAL :: check=.false.
  INTEGER :: jj, kk
  REAL    :: jireal(1:Nt*Nz), jiimag(1:Nt*Nz), arg, ca, sa

  if( check ) then ; jireal=ijreal ; jiimag=ijimag
  endif

  Ntz = Nt * Nz

  !Copy real arrays to complex
  do jj=1,Nz
    cplxin(:,jj) = CMPLX(ijreal((jj-1)*Nt+1:jj*Nt), ijimag((jj-1)*Nt+1:jj*Nt),KIND=C_DOUBLE_COMPLEX)
  enddo

  call fftw_execute_dft(planf, cplxin, cplxout) !Forward transform

  !Copy complex result back to real arrays, normalize
  do jj=1,Nz
    ijreal((jj-1)*Nt+1:jj*Nt) = real(cplxout(:,jj),KIND=C_DOUBLE_COMPLEX)/Ntz
    ijimag((jj-1)*Nt+1:jj*Nt) = aimag(cplxout(:,jj))/Ntz
  enddo
  
  cfmn=zero ; sfmn=zero ; efmn=zero ; ofmn=zero
  
  do imn = 1,mn ; mm = im(imn) ; nn = in(imn) / Nfp
   
   if    ( mm.gt.0 .and. nn.gt.0 ) then
    
    efmn(imn) =   ijreal(1+(Nt-mm)+(   nn)*Nt) + ijreal(1+(   mm)+(Nz-nn)*Nt)
    ofmn(imn) =   ijimag(1+(Nt-mm)+(   nn)*Nt) - ijimag(1+(   mm)+(Nz-nn)*Nt)
    cfmn(imn) =   ijimag(1+(Nt-mm)+(   nn)*Nt) + ijimag(1+(   mm)+(Nz-nn)*Nt)
    sfmn(imn) = - ijreal(1+(Nt-mm)+(   nn)*Nt) + ijreal(1+(   mm)+(Nz-nn)*Nt) 
    
   elseif( mm.gt.0 .and. nn.eq.0 ) then
    
    efmn(imn) =   ijreal(1+(Nt-mm)           ) + ijreal(1+(   mm)           )
    ofmn(imn) =   ijimag(1+(Nt-mm)           ) - ijimag(1+(   mm)           )
    cfmn(imn) =   ijimag(1+(Nt-mm)           ) + ijimag(1+(   mm)           )
    sfmn(imn) = - ijreal(1+(Nt-mm)           ) + ijreal(1+(   mm)           )
    
   elseif( mm.gt.0 .and. nn.lt.0 ) then
    
    efmn(imn) =   ijreal(1+(Nt-mm)+(Nz+nn)*Nt) + ijreal(1+(   mm)+(  -nn)*Nt)
    ofmn(imn) =   ijimag(1+(Nt-mm)+(Nz+nn)*Nt) - ijimag(1+(   mm)+(  -nn)*Nt)
    cfmn(imn) =   ijimag(1+(Nt-mm)+(Nz+nn)*Nt) + ijimag(1+(   mm)+(  -nn)*Nt)
    sfmn(imn) = - ijreal(1+(Nt-mm)+(Nz+nn)*Nt) + ijreal(1+(   mm)+(  -nn)*Nt) 
    
   elseif( mm.eq.0 .and. nn.gt.0 ) then
    
    efmn(imn) =   ijreal(1+        (Nz-nn)*Nt) + ijreal(1+        (   nn)*Nt)
    ofmn(imn) = - ijimag(1+        (Nz-nn)*Nt) + ijimag(1+        (   nn)*Nt)
    cfmn(imn) =   ijimag(1+        (Nz-nn)*Nt) + ijimag(1+        (   nn)*Nt)
    sfmn(imn) =   ijreal(1+        (Nz-nn)*Nt) - ijreal(1+        (   nn)*Nt)
    
   elseif( mm.eq.0 .and. nn.eq.0 ) then
    
    efmn(imn) =   ijreal(1)
    ofmn(imn) =     zero
    cfmn(imn) =   ijimag(1)
    sfmn(imn) =     zero
    
   endif
   
  enddo

  if( .not.check ) return
  
  ijreal(1:Ntz)=zero ; ijimag(1:Ntz)=zero
  
  do jj=0,Nt-1
   do kk=0,Nz-1
    do imn=1,mn ; arg=im(imn)*jj*pi2/Nt-in(imn)*kk*pi2nfp/Nz ; ca=cos(arg) ; sa=sin(arg)
     
     ijreal(1+jj+kk*Nt) = ijreal(1+jj+kk*Nt) + efmn(imn)*ca + ofmn(imn)*sa
     ijimag(1+jj+kk*Nt) = ijimag(1+jj+kk*Nt) + cfmn(imn)*ca + sfmn(imn)*sa
     
    enddo
   enddo
  enddo
  
  write(*,'("tfft : Fourier reconstruction error = "2es15.5)')sqrt(sum((ijreal-jireal)**2)/Ntz),sqrt(sum((ijimag-jiimag)**2)/Ntz)

  return

end subroutine tfft

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{\type{invfft}}

!latex \begin{enumerate}

!latex \item Given the Fourier harmonics, the data on a regular angular grid are constructed.

!latex \item This is the inverse routine to \type{tfft}.

!latex \end{enumerate}

subroutine invfft( mn, im, in, efmn, ofmn, cfmn, sfmn, Nt, Nz, ijreal, ijimag )

  use constants, only : zero, half, one
  use inputlist, only : Nfp
  use fftw_interface
  implicit none
  
  INTEGER, intent(in)  :: mn, im(mn), in(mn)
  REAL   , intent(in)  :: efmn(mn), ofmn(mn), cfmn(mn), sfmn(mn)
  INTEGER, intent(in)  :: Nt, Nz
  REAL   , intent(out) :: ijreal(Nt*Nz), ijimag(Nt*Nz) ! output real space;
  
  INTEGER              :: Ntz, imn, jj, kk, mm, nn
  
  Ntz = Nt*Nz ; ijreal(1:Ntz) = zero ; ijimag(1:Ntz) = zero

  do imn = 1,mn ; mm = im(imn) ; nn = in(imn) / Nfp
   
   if    ( mm.gt.0 .and. nn.gt.0 ) then
    
    ijreal(1+(Nt-mm)+(   nn)*Nt) = (   efmn(imn) - sfmn(imn) ) * half
    ijreal(1+(   mm)+(Nz-nn)*Nt) = (   efmn(imn) + sfmn(imn) ) * half

    ijimag(1+(Nt-mm)+(   nn)*Nt) = (   ofmn(imn) + cfmn(imn) ) * half
    ijimag(1+(   mm)+(Nz-nn)*Nt) = ( - ofmn(imn) + cfmn(imn) ) * half

   elseif( mm.gt.0 .and. nn.eq.0 ) then
    
    ijreal(1+(Nt-mm)           ) = (   efmn(imn) - sfmn(imn) ) * half
    ijreal(1+(   mm)           ) = (   efmn(imn) + sfmn(imn) ) * half

    ijimag(1+(Nt-mm)           ) = ( + ofmn(imn) + cfmn(imn) ) * half
    ijimag(1+(   mm)           ) = ( - ofmn(imn) + cfmn(imn) ) * half

   elseif( mm.gt.0 .and. nn.lt.0 ) then
    
    ijreal(1+(Nt-mm)+(Nz+nn)*Nt) = (   efmn(imn) - sfmn(imn) ) * half
    ijreal(1+(   mm)+(  -nn)*Nt) = (   efmn(imn) + sfmn(imn) ) * half

    ijimag(1+(Nt-mm)+(Nz+nn)*Nt) = ( + ofmn(imn) + cfmn(imn) ) * half
    ijimag(1+(   mm)+(  -nn)*Nt) = ( - ofmn(imn) + cfmn(imn) ) * half

   elseif( mm.eq.0 .and. nn.gt.0 ) then
    
    ijreal(1+        (Nz-nn)*Nt) = ( + efmn(imn) + sfmn(imn) ) * half
    ijreal(1+        (   nn)*Nt) = ( + efmn(imn) - sfmn(imn) ) * half
    
    ijimag(1+        (Nz-nn)*Nt) = ( - ofmn(imn) + cfmn(imn) ) * half
    ijimag(1+        (   nn)*Nt) = ( + ofmn(imn) + cfmn(imn) ) * half

   elseif( mm.eq.0 .and. nn.eq.0 ) then
    
    ijreal(1) = efmn(imn)
    ijimag(1) = cfmn(imn)
    
   endif
   
  enddo

  !Copy real arrays to complex
  do jj=1,Nz
    cplxin(:,jj) = CMPLX(ijreal((jj-1)*Nt+1:jj*Nt), ijimag((jj-1)*Nt+1:jj*Nt),KIND=C_DOUBLE_COMPLEX)
  enddo

  call fftw_execute_dft(planb, cplxin, cplxout) !Inverse transform

  !Copy complex result back to real arrays
  do jj=1,Nz
    ijreal((jj-1)*Nt+1:jj*Nt) = real(cplxout(:,jj),KIND=C_DOUBLE_COMPLEX)
    ijimag((jj-1)*Nt+1:jj*Nt) = aimag(cplxout(:,jj))
  enddo

  return

end subroutine invfft

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef DELETETHIS

!l tex \subsection{\type{svdcmp}} ! not used; SRH: 27 Feb 18;

subroutine svdcmp(a,m,n,mp,np,w,v)
  use constants,only:zero,one

      implicit none
      integer nMaX,M,n,MP,nP,i,j,jj,k,l,its,nM
      parameter (nMaX=500)
      REAL ::  a(MP,nP),w(nP),v(nP,nP),rv1(nMaX)
      REAL ::  c,F,g,H,s,X,Y,Z,scale,anorM,pythag,oone

     !stop "svdcmp : to be deleted?"

      g=zero
      scale=zero
      anorm=zero
      do 25 i=1,n
         l=i+1
         rv1(i)=scale*g
         g=zero
         s=zero
         scale=zero
         if(i.le.m)then
            do 11 k=i,m
               scale=scale+abs(a(k,i))
11          continue
            if(scale.ne.zero)then
               do 12 k=i,m
                  a(k,i)=a(k,i)/scale
                  s=s+a(k,i)*a(k,i)
12             continue
               f=a(i,i)
               g=-sign(sqrt(s),f)
               h=f*g-s
               a(i,i)=f-g
               do 15 j=l,n
                  s=zero
                  do 13 k=i,m
                     s=s+a(k,i)*a(k,j)
13                continue
                  f=s/h
                  do 14 k=i,m
                     a(k,j)=a(k,j)+f*a(k,i)
14                continue
15             continue
               do 16 k=i,m
                  a(k,i)=scale*a(k,i)
16             continue
            endif
         endif
         w(i)=scale *g
         g=zero
         s=zero
         scale=zero
         if((i.le.m).and.(i.ne.n))then
            do 17 k=l,n
               scale=scale+abs(a(i,k))
17          continue
            if(scale.ne.zero)then
               do 18 k=l,n
                  a(i,k)=a(i,k)/scale
                  s=s+a(i,k)*a(i,k)
18             continue
               f=a(i,l)
               g=-sign(sqrt(s),f)
               h=f*g-s
               a(i,l)=f-g
               do 19 k=l,n
                  rv1(k)=a(i,k)/h
19             continue
               do 23 j=l,m
                  s=zero
                  do 21 k=l,n
                     s=s+a(j,k)*a(i,k)
21                continue
                  do 22 k=l,n
                     a(j,k)=a(j,k)+s*rv1(k)
22                continue
23             continue
               do 24 k=l,n
                  a(i,k)=scale*a(i,k)
24             continue
            endif
         endif
         anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
25    continue
      do 32 i=n,1,-1
         if(i.lt.n)then
            if(g.ne.zero)then
               do 26 j=l,n
                  v(j,i)=(a(i,j)/a(i,l))/g
26             continue
               do 29 j=l,n
                  s=zero
                  do 27 k=l,n
                     s=s+a(i,k)*v(k,j)
27                continue
                  do 28 k=l,n
                     v(k,j)=v(k,j)+s*v(k,i)
28                continue
29             continue
            endif
            do 31 j=l,n
               v(i,j)=zero
               v(j,i)=zero
31          continue
         endif
         v(i,i)=one
         g=rv1(i)
         l=i
32    continue
      do 39 i=min(m,n),1,-1
         l=i+1
         g=w(i)
         do 33 j=l,n
            a(i,j)=zero
33       continue
         if(g.ne.zero)then
            g=one/g
            do 36 j=l,n
               s=zero
               do 34 k=l,m
                  s=s+a(k,i)*a(k,j)
34             continue
               f=(s/a(i,i))*g
               do 35 k=i,m
                  a(k,j)=a(k,j)+f*a(k,i)
35             continue
36          continue
            do 37 j=i,m
               a(j,i)=a(j,i)*g
37          continue        
         else
            do 38 j= i,m
               a(j,i)=zero
38          continue
         endif
         a(i,i)=a(i,i)+one
39    continue
      do 49 k=n,1,-1
         do 48 its=1,30
            do 41 l=k,1,-1
               nm=l-1
               if((abs(rv1(l))+anorm).eq.anorm)  goto 2
               if((abs(w(nm))+anorm).eq.anorm)  goto 1
41          continue
1           c=zero
            s=one
            do 43 i=l,k
               f=s*rv1(i)
               rv1(i)=c*rv1(i)
               if((abs(f)+anorm).eq.anorm) goto 2
               g=w(i)
               h=pythag(f,g)
               w(i)=h
               h=one/h
               c= (g*h)
               s=-(f*h)
               do 42 j=1,m
                  y=a(j,nm)
                  z=a(j,i)
                  a(j,nm)=(y*c)+(z*s)
                  a(j,i)=-(y*s)+(z*c)
42             continue
43          continue
2           z=w(k)
            if(l.eq.k)then
               if(z.lt.zero)then
                  w(k)=-z
                  do 44 j=1,n
                     v(j,k)=-v(j,k)
44                continue
               endif
               goto 3
            endif
            if(its.eq.30) stop "svdcmp : no convergence"
            x=w(l)
            nm=k-1
            y=w(nm)
            g=rv1(nm)
            h=rv1(k)
            f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
            oone=one
            g=pythag(f,oone)
            f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
            c=one
            s=one
            do 47 j=l,nm
               i=j+1
               g=rv1(i)
               y=w(i)
               h=s*g
               g=c*g
               z=pythag(f,h)
               rv1(j)=z
               c=f/z
               s=h/z
               f= (x*c)+(g*s)
               g=-(x*s)+(g*c)
               h=y*s
               y=y*c
               do 45 jj=1,n
                  x=v(jj,j)
                  z=v(jj,i)
                  v(jj,j)= (x*c)+(z*s)
                  v(jj,i)=-(x*s)+(z*c)
45             continue
               z=pythag(f,h)
               w(j)=z
               if(z.ne.zero)then
                  z=one/z
                  c=f*z
                  s=h*z
               endif
               f= (c*g)+(s*y)
               x=-(s*g)+(c*y)
               do 46 jj=1,m
                  y=a(jj,j)
                  z=a(jj,i)
                  a(jj,j)= (y*c)+(z*s)
                  a(jj,i)=-(y*s)+(z*c)
46             continue
47          continue
            rv1(l)=zero
            rv1(k)=f
            w(k)=x
48       continue
3        continue
49    continue
      return
      end

#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
#ifdef DELETETHIS
   
!l tex \subsection{\type{pythag}} ! not used; SRH: 27 Feb 18;

REAL function pythag(a,b)  
  implicit none
  REAL ::  a,b
  REAL ::  absa,absb

 !stop "pythag : to be deleted?"

  absa=abs(a)
  absb=abs(b)
  if(absa.gt.absb) then
   pythag=absa*sqrt(1.+(absb/absa)**2)
  else
   if(absb.eq.0.) then
    pythag=0.
   else
    pythag=absb*sqrt(1.+(absa/absb)**2)
   endif
  endif
  return
end function pythag

#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
#ifdef DELETETHIS

!l tex \subsection{\type{svbksb}} ! not used; SRH: 27 Feb 18;

subroutine svbksb(u,w,v,M,n,MP,nP,b,X)
  
      implicit none
      integer nMaX,M,n,MP,nP,i,j,jj
      parameter (nMaX=10000)
      REAL, intent(in) :: b(MP)
      REAL ::  u(MP,nP),w(nP),v(nP,nP),X(nP),tMP(nMaX)
      REAL ::  s

     !stop "svbksb : to be deleted?"
  
      do 12 j=1,n
         s=0.
         if(w(j).ne.0.)then
            do 11 i=1,M
               s=s+u(i,j)*b(i)
11          continue
            s=s/w(j)
         endif
         tMP(j)=s
12    continue
      do 14 j=1,n
         s=0.
         do 13 jj=1,n
            s=s+v(j,jj)*tMP(jj)
13       continue
         X(j)=s
14    continue
      return
      end

#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
      
#ifdef DELETETHIS

!l tex \subsection{\type{sort}} ! not used; SRH: 27 Feb 18;

      subroutine sort(n,ra)

      implicit none
      integer n,l,ir,i,j
      REAL ::  ra(n),rra

     !stop "sort : to be deleted?"

      if(n.eq.1) return
      l=n/2+1
      ir=n
10    continue
      if(l.gt.1)then
         l=l-1
         rra=ra(l)
      else
         rra=ra(ir)
         ra(ir)=ra(1)
         ir=ir-1
         if(ir.eq.1)then
            ra(1)=rra
            return
         endif
      endif
      i=l
      j=l+l
20    if(j.le.ir)then
         if(j.lt.ir)then
            if(ra(j).lt.ra(j+1))j=j+1
         endif
         if(rra.lt.ra(j))then
            ra(i)=ra(j)
            i=j
            j=j+j
         else
            j=ir+1
         endif
         goto 20
      endif
      ra(i)=rra
      goto 10
      end

#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef DELETETHIS
     
!l tex \subsection{\type{singvalues}} ! not used; SRH: 27 Feb 18;

      subroutine singvalues(nrow,ncol,Mat,b,sx,cutoff,wsvd) ! nrow = nconstraints ; ncol = nfreedom
      implicit none
      INTEGER, intent(in) :: nrow,ncol
      REAL,intent(in) ::  Mat(nrow,ncol),b(nrow)
      integer i,nev

      REAL :: Mato(nrow,ncol)
      REAL ::  sx(ncol)
      REAL ::  vsvd(ncol,ncol),wsvd(ncol),wsvdc(ncol)
      REAL ::  cutoff,wmax,wmin

      sx=0.0;wsvd=0.0;vsvd=0.0;wsvdc=0.0;wmax=0.0;Mato=Mat

      call svdcmp(Mato,nrow,ncol,nrow,ncol,wsvd,vsvd)
      wsvdc=wsvd
      call sort(ncol,wsvdc)
      wmax=wsvdc(ncol)
      wmin=abs(wmax)*cutoff
      wsvdc=0.0;nev=0
      do i=1,ncol
       if(abs(wsvd(i)).ge.wmin) then;wsvdc(i)=wsvd(i);nev=nev+1
       endif
      enddo
      call svbksb(Mato,wsvdc,vsvd,nrow,ncol,nrow,ncol,b,sx)
     !Mat=Mato
      call sort(ncol,wsvd)
      return
      end subroutine singvalues

#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!l tex \end{itemize}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
