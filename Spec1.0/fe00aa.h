!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Calculates finite element functions, and their derivatives, on radial sub-sub-grid.
!latex The quantities calculated are the $\varphi_{l,p}$ for $l=0,1$ and $p=0,$\verb+Nofe+
!latex which are stored in array \verb+felr(0:1,0:Nofe,0:1,1:Nj)+ (which is allocated in \verb+readin+).

!latex \item This routine is only called once, from \verb+xspech+.

!latex \item Should check that \mbox{$\varphi_{1,p}(x)=(-1)^p \varphi_{0,p}(1-x)$} is satisfied.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine fe00aa( Iquad )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, one, two, three, four, six

  use numerical, only : 

  use fileunits, only : ounit

  use inputlist, only : Wfe00aa, Nofe

  use cputiming, only : Tfe00aa

  use allglobal, only : myid, cpus, felr
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS

  INTEGER, intent(in) :: Iquad

  INTEGER             :: id01bcf, itype, pp
  REAL                :: aa, bb, cc, dd ! supplemental parameters for Gaussian integration;
  REAL                :: gwt(1:Iquad), xx1(1:Iquad), xx2(1:Iquad), xx3(1:Iquad), xx4(1:Iquad), xx5(1:Iquad), xx6(1:Iquad), xx7(1:Iquad)

  BEGIN(fe00aa)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Radial quadrature is performed using Gaussian quadrature.
!latex The resolution of the approximation is controlled by \verb+Iquad+, and the weights and abscissae are returned by NAG routine \verb+D01BCF+.
!latex At present, only \verb+itype=0+ (Gauss-Legrendre) and \verb+a=0+, \verb+b=1+, \verb+c=0+ and \verb+d=0+ are supplied as input to \verb+D01BCF+. 

 !ALLOCATE(gwt,(1:Iquad))
 !ALLOCATE(xx1,(1:Iquad))
  
  itype = 0 ; aa = zero ; bb = one ; cc = zero ; dd = zero
  id01bcf=0 ; call D01BCF( itype , aa , bb , cc , dd , Iquad , gwt(1:Iquad) , xx1(1:Iquad) , id01bcf ) ! everything here local to sub-grid; x in [0,1];
  
  FATALMESS(fe00aa,id01bcf.ne.0,an error has occured with the Gaussian quadrature routine D01CBF)

  if( Wfe00aa .and. myid.eq.0 ) then
   cput = GETTIME
   write(ounit,'("fe00aa : ",f10.2," : [ a , labscis , b ]   ="9999f15.11)')cput-cpus,aa,xx1,bb
   write(ounit,'("fe00aa : ",f10.2," :       weight          ="15x,9998f15.11)')cput-cpus,gwt
  endif

 !DEALLOCATE(gwt) ! only the felr are assigned here; the weights are not required;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item The array \verb+felr(0:1,0:Nofe,0:1,1:Iquad)+ is constructed as follows:
!latex \begin{enumerate} 
!latex \item the first  argument is the left-right index:            \verb+felr(0,p,0,j)+$=\varphi_{0,p}(x_j)$ and \verb+felr(1,p,0,j)+$=\varphi_{1,p}(x_j)$; 
!latex \item the second argument is the order of the basis function:
!latex \item the third  argument is the derivative:                  \verb+felr(0,p,1,j)+$=\varphi_{L,p}^{\prime}(x_j)$;
!latex \item and the fourth argument indicates the radial sub-sub-grid index.
!latex \end{enumerate}

!latex \item Note that the derivatives {\bf with respect to x} are returned, where $x=(s-s_{l,i-1})/\Delta s_l$ and $\Delta s_l=(s_{l,i}-s_{l,i-1})$.
!latex The radial scaling factor $\Delta s_l$ needs to be included elsewhere.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item The radial basis functions used are: \begin{itemize}

  select case( Nofe )
   
  case(0) !latex \item \verb+Nofe=0+: linear  basis functions;
   
   felr(0,0,0,1:Iquad) = - xx1 + one
   felr(0,0,1,1:Iquad) = - one           
   
   felr(1,0,0,1:Iquad) =   xx1       
   felr(1,0,1,1:Iquad) =   one           
   
  case(1) !latex \item \verb+Nofe=1+: cubic basis functions;
   
  !ALLOCATE(xx2,(1:Iquad))
  !ALLOCATE(xx3,(1:Iquad))
   
   xx2 = xx1 * xx1 ; xx3 = xx2 * xx1
   
!latex \be \begin{array}{cccccccccccccccccccccc}
!latex \varphi_{0,0} &=&  2 & x^3 & - & 3 & x^2 &   &   & + & 1 \\
!latex \varphi_{0,1} &=&    & x^3 & - & 2 & x^2 & + & x &   &   \\
!latex \varphi_{1,0} &=& -2 & x^3 & + & 3 & x^2 &   &   &   &   \\
!latex \varphi_{1,1} &=&    & x^3 & - &   & x^2 &   &   &   &   \end{array} \ee

   felr(0,0,0,1:Iquad) =     two * xx3 - three * xx2            + one 
   felr(0,0,1,1:Iquad) =     six * xx2 -   six * xx1                 
   
   felr(0,1,0,1:Iquad) =           xx3 -   two * xx2 +       xx1      
   felr(0,1,1,1:Iquad) =   three * xx2 -  four * xx1 + one           
   
   felr(1,0,0,1:Iquad) = -   two * xx3 + three * xx2                  
   felr(1,0,1,1:Iquad) = -   six * xx2 +   six * xx1
   
   felr(1,1,0,1:Iquad) =           xx3 -         xx2                  
   felr(1,1,1,1:Iquad) =   three * xx2 -   two * xx1

  !DEALLOCATE(xx3)
  !DEALLOCATE(xx2)
   
  case(2) !latex \item \verb+Nofe=2+: quintic basis function;
   
  !ALLOCATE(xx2,(1:Iquad))
  !ALLOCATE(xx3,(1:Iquad))
  !ALLOCATE(xx4,(1:Iquad))
  !ALLOCATE(xx5,(1:Iquad))
   
   xx2 = xx1 * xx1 ; xx3 = xx2 * xx1 ; xx4 = xx3 * xx1 ; xx5 = xx4 * xx1

!latex \be \begin{array}{cccccccccccccccccccccccccccccccccc}
!latex \varphi_{L,0} & = & - 6    & x^5   & + & 15   & x^4   & - & 10   & x^3 &   &     &     &   &   & + & 1 \\
!latex \varphi_{L,1} & = & - 3    & x^5   & + &  8   & x^4   & - &  6   & x^3 &   &     &     & + & x &   &   \\
!latex \varphi_{L,2} & = & - 1/2  & x^5   & + &  3/2 & x^4   & - &  3/2 & x^3 & + & 1/2 & x^2 &   &   &   &   \\
!latex \varphi_{R,0} & = &   6    & x^5   & - & 15   & x^4   & + & 10   & x^3 &   &     &     &   &   &   &   \\
!latex \varphi_{R,0} & = & - 3    & x^5   & + &  7   & x^4   & - &  4   & x^3 &   &     &     &   &   &   &   \\
!latex \varphi_{R,0} & = &   1/2  & x^5   & - &      & x^4   & + & 1/2  & x^3 &   &     &     &   &   &   &   \end{array} \ee   

   felr(0,0,0,1:Iquad) =  -  6   * xx5 + 15   * xx4 - 10   * xx3                           +  1  
   felr(0,0,1,1:Iquad) =  - 30   * xx4 + 60   * xx3 - 30   * xx2                            

   felr(0,1,0,1:Iquad) =  -  3   * xx5 +  8   * xx4 -  6   * xx3               +        xx1       
   felr(0,1,1,1:Iquad) =  - 15   * xx4 + 32   * xx3 - 18   * xx2               +  1              

   felr(0,2,0,1:Iquad) =  -  0.5 * xx5 +  1.5 * xx4 -  1.5 * xx3 +  0.5 * xx2                   
   felr(0,2,1,1:Iquad) =  -  2.5 * xx4 +  6   * xx3 -  4.5 * xx2 +        xx1                   

   felr(1,0,0,1:Iquad) =     6   * xx5 - 15   * xx4 + 10   * xx3                                  
   felr(1,0,1,1:Iquad) =    30   * xx4 - 60   * xx3 + 30   * xx2                                  

   felr(1,1,0,1:Iquad) =  -  3   * xx5 +  7   * xx4 -  4   * xx3                                  
   felr(1,1,1,1:Iquad) =  - 15   * xx4 + 28   * xx3 - 12   * xx2                                  

   felr(1,2,0,1:Iquad) =     0.5 * xx5 -        xx4 +  0.5 * xx3                                  
   felr(1,2,1,1:Iquad) =     2.5 * xx4 -  4   * xx3 +  1.5 * xx2                                  

  !DEALLOCATE(xx5)
  !DEALLOCATE(xx4)
  !DEALLOCATE(xx3)
  !DEALLOCATE(xx2)
   
  case(3) !latex \item \verb+Nofe=3+: septemtic basis function;
   
  !ALLOCATE(xx2,(1:Iquad))
  !ALLOCATE(xx3,(1:Iquad))
  !ALLOCATE(xx4,(1:Iquad))
  !ALLOCATE(xx5,(1:Iquad))
  !ALLOCATE(xx6,(1:Iquad))
  !ALLOCATE(xx7,(1:Iquad))
   
   xx2 = xx1 * xx1 ; xx3 = xx2 * xx1 ; xx4 = xx3 * xx1 ; xx5 = xx4 * xx1 ; xx6=xx5 * xx1 ; xx7=xx6 * xx1
   
!latex \be \begin{array}{cccccccccccccccccccccccccccccccccc}
!latex \varphi_{0,0} &=&  20   & x^7 & - & 70   & x^6 & + & 84   & x^5 & - & 35   & x^4 &   &     &     &   &     &     &   &   &   & + & 1 \\
!latex \varphi_{0,1} &=&  10   & x^7 & - & 36   & x^6 & + & 45   & x^5 & - & 20   & x^4 &   &     &     &   &     &     & + & 1 & x &   &   \\
!latex \varphi_{0,2} &=&   2   & x^7 & - & 15/2 & x^6 & + & 10   & x^5 & - &  5   & x^4 &   &     &     & + & 1/2 & x^2 &   &   &   &   &   \\
!latex \varphi_{0,3} &=&   1/6 & x^7 & - &  2/3 & x^6 & + &  1   & x^5 & - &  2/3 & x^4 & + & 1/6 & x^3 &   &     &     &   &   &   &   &   \\
!latex \varphi_{1,0} &=& -20   & x^7 & + & 70   & x^6 & - & 84   & x^5 & + & 35   & x^4 &   &     &     &   &     &     &   &   &   &   &   \\
!latex \varphi_{1,1} &=&  10   & x^7 & - & 34   & x^6 & + & 39   & x^5 & - & 15   & x^4 &   &     &     &   &     &     &   &   &   &   &   \\
!latex \varphi_{1,2} &=& - 2   & x^7 & + & 13/2 & x^6 & - &  7   & x^5 & + &  5/2 & x^4 &   &     &     &   &     &     &   &   &   &   &   \\
!latex \varphi_{1,3} &=&   1/6 & x^7 & - &  1/2 & x^6 & + &  1/2 & x^5 & - &  1/6 & x^4 &   &     &     &   &     &     &   &   &   &   &   \end{array} \ee

   felr(0,0,0,1:Iquad) = +  20     * xx7 -  70     * xx6 +  84   * xx5 -  35     * xx4                                           + 1    
   felr(0,0,1,1:Iquad) = + 140     * xx6 - 420     * xx5 + 420   * xx4 - 140     * xx3
   
   felr(0,1,0,1:Iquad) = +  10     * xx7 -  36     * xx6 +  45   * xx5 -  20     * xx4                          + 1 * xx1
   felr(0,1,1,1:Iquad) = +  70     * xx6 - 216     * xx5 + 225   * xx4 -  80     * xx3                          + 1    
   
   felr(0,2,0,1:Iquad) = +   2     * xx7 -   7.5   * xx6 +  10   * xx5 -  5     * xx4               + 0.5 * xx2                       
   felr(0,2,1,1:Iquad) = +  14     * xx6 -  45     * xx5 +  50   * xx4 - 20     * xx3               +       xx1
   
   felr(0,3,0,1:Iquad) = +   1.0/6 * xx7 -   2.0/3 * xx6 +         xx5 -   2.0/3 * xx4 + 1.0/6 * xx3                                    
   felr(0,3,1,1:Iquad) = +   7.0/6 * xx6 -   4     * xx5 +   5   * xx4 -   8.0/3 * xx3 + 3.0/6 * xx2
   
   felr(1,0,0,1:Iquad) = -  20     * xx7 +  70     * xx6 -  84   * xx5 +  35     * xx4
   felr(1,0,1,1:Iquad) = - 140     * xx6 + 420     * xx5 - 420   * xx4 + 140     * xx3
   
   felr(1,1,0,1:Iquad) = +  10     * xx7 -  34     * xx6 +  39   * xx5 -  15     * xx4
   felr(1,1,1,1:Iquad) = +  70     * xx6 - 204     * xx5 + 195   * xx4 -  60     * xx3

   felr(1,2,0,1:Iquad) = -   2     * xx7 +   6.5   * xx6 -   7   * xx5 +   2.5   * xx4
   felr(1,2,1,1:Iquad) = -  14     * xx6 +  39     * xx5 -  35   * xx4 +  10     * xx3
   
   felr(1,3,0,1:Iquad) = +   1.0/6 * xx7 -   0.5   * xx6 +   0.5 * xx5 -   1.0/6 * xx4
   felr(1,3,1,1:Iquad) = +   7.0/6 * xx6 -   3     * xx5 +   2.5 * xx4 -   4.0/6 * xx3
   
  !DEALLOCATE(xx7)
  !DEALLOCATE(xx6)
  !DEALLOCATE(xx5)
  !DEALLOCATE(xx4)
  !DEALLOCATE(xx3)
  !DEALLOCATE(xx2)
   
  case default  !latex \item otherwise there is an error; \end{itemize}
   
   FATAL(fe00aa : illegal Nofe,.true.)
   
  end select
  
  if( Wfe00aa .and. myid.eq.0 ) then
   cput = GETTIME
   do pp = 0,Nofe
    write(ounit,'("fe00aa : ",f10.2," : (-1)^"i1" varphi_{L,"i1"}(1-x)="15x,99f15.11)')cput-cpus,pp,pp,felr(0,pp,0,Iquad:1:-1) * (-1)**pp
    write(ounit,'("fe00aa : ",f10.2," :        varphi_{R,"i1"}(  x)="15x,99f15.11)')cput-cpus,   pp,felr(1,pp,0,1:Iquad:+1)
   enddo
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
 !DEALLOCATE(xx1)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(fe00aa)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine fe00aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
