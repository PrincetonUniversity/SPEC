!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Calculates finite element functions, and their derivatives, at arbitrary radius.
!latex \item This routine is a direct copy of \verb+fe00aa+,
!latex whereas \verb+fe00aa+ calculates the radial basis functions on the fixed radial sub-sub-grid,
!latex and assigns values to the \verb+felr+ arrays,
!latex \verb+fe00ab+ calculates the radial basis functions at an arbitrary radial location.
!latex This information is is required for field line integrations, for example.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine fe00ab( lvol , xx , Nofe , vphi )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  use constants, only : zero, one
  use numerical, only :
  use fileunits, only : ounit
  use inputlist, only : Wfe00ab
  use cputiming, only : Tfe00ab
  use allglobal, only : myid, cpus

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS

  INTEGER, intent(in)  :: lvol, Nofe
  REAL   , intent(in)  :: xx ! local interpolation coordinate;
  REAL   , intent(out) :: vphi(0:1,0:3,0:2)

  REAL                 :: xx2, xx3, xx4, xx5, xx6, xx7

  BEGIN(fe00ab)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item The output quantity \verb+vphi(0:1,0:Nofe,0:2)+ contains the following: \newline
!latex \verb+vphi(0,p,0)+$=\varphi_{0,p}$                                          \newline
!latex \verb+vphi(0,p,1)+$=\varphi_{0,p}^{\prime}$                                 \newline
!latex \verb+vphi(0,p,2)+$=\varphi_{0,p}^{\prime\prime}$                           \newline
!latex \verb+vphi(1,p,0)+$=\varphi_{1,p}$                                          \newline
!latex \verb+vphi(1,p,1)+$=\varphi_{1,p}^{\prime}$                                 \newline
!latex \verb+vphi(1,p,2)+$=\varphi_{1,p}^{\prime\prime}$

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Note that the derivatives {\bf with respect to x} are returned.
!latex Therefore, the scaling factor needs to be included elsewhere.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( xx.lt.zero .or. xx.gt.one ) then ; cput = GETTIME ; write(ounit,1000)cput-cpus,myid,lvol,xx
  endif

1000 format("fe00ab : ",f10.2,"s : WF(?) ; myid=",i3," ; lvol=",i3," ; xx="es24.16" ; ")

  FATALMESS(fe00ab,xx.lt.zero.or.xx.gt.one,local interpolation variable is out of range)

  FATALMESS(fe0baa,Nofe.gt.3,need to change dimensions of vphi,Nofe.gt.3)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  vphi(0:1,0:3,0:2) = zero ! ensure intent out quantity is assigned;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  select case(Nofe)
   
  case(0) ! linear basis functions;
   
   vphi(0,0,0) = - xx + 1
   vphi(0,0,1) = - 1             
   
   vphi(1,0,0) =   xx 
   vphi(1,0,1) =   1        
   
  case(1)
   
   xx2 = xx * xx ; xx3 = xx2 * xx 
   
   vphi(0,0,0) =    2 * xx3 - 3 * xx2          + 1   
   vphi(0,0,1) =    6 * xx2 - 6 * xx                  
   vphi(0,0,2) =   12 * xx  - 6                    
   
   vphi(0,1,0) =        xx3 - 2 * xx2 + 1 * xx       
   vphi(0,1,1) =    3 * xx2 - 4 * xx  + 1            
   vphi(0,1,2) =    6 * xx  - 4   
   
   vphi(1,0,0) = -  2 * xx3 + 3 * xx2                
   vphi(1,0,1) = -  6 * xx2 + 6 * xx 
   vphi(1,0,2) = - 12 * xx  + 6                     
   
   vphi(1,1,0) =        xx3 -     xx2                
   vphi(1,1,1) =    3 * xx2 - 2 * xx 
   vphi(1,1,2) =    6 * xx  - 2   
   
  case(2) ! quintic basis functions;
   
   xx2 = xx * xx ; xx3 = xx2 * xx ; xx4 = xx3 * xx ; xx5 = xx4 * xx 
   
   vphi(0,0,0) = -   6   * xx5 + 15   * xx4 - 10   * xx3                       + 1  
   vphi(0,0,1) = -  30   * xx4 + 60   * xx3 - 30   * xx2                            
   vphi(0,0,2) = - 120   * xx3 +180   * xx2 - 60   * xx                               
   
   vphi(0,1,0) = -   3   * xx5 +  8   * xx4 -  6   * xx3                  + xx       
   vphi(0,1,1) = -  15   * xx4 + 32   * xx3 - 18   * xx2              + 1              
   vphi(0,1,2) = -  60   * xx3 + 96   * xx2 - 36   * xx                                     
   
   vphi(0,2,0) = -   0.5 * xx5 +  1.5 * xx4 -  1.5 * xx3 +  0.5 * xx2                   
   vphi(0,2,1) = -   2.5 * xx4 +  6   * xx3 -  4.5 * xx2 +        xx                      
   vphi(0,2,2) = -  10   * xx3 + 18   * xx2 -  9   * xx  +  1                             
   
   vphi(1,0,0) =     6   * xx5 - 15   * xx4 + 10   * xx3                                  
   vphi(1,0,1) =    30   * xx4 - 60   * xx3 + 30   * xx2                                  
   vphi(1,0,2) =   120   * xx3 -180   * xx2 + 60   * xx                                     
   
   vphi(1,1,0) = -   3   * xx5 +  7   * xx4 -  4   * xx3                                  
   vphi(1,1,1) = -  15   * xx4 + 28   * xx3 - 12   * xx2                                  
   vphi(1,1,2) = -  60   * xx3 + 84   * xx2 - 24   * xx                                     
   
   vphi(1,2,0) =     0.5 * xx5 -        xx4 +  0.5 * xx3                                  
   vphi(1,2,1) =     2.5 * xx4 -  4   * xx3 +  1.5 * xx2
   vphi(1,2,2) =    10   * xx3 - 12   * xx2 +  3   * xx 
   
  case(3) ! septemtic basis functions;
   
   xx2 = xx * xx  ; xx3 = xx2 * xx ; xx4 = xx3 * xx ; xx5 = xx4 * xx ; xx6 = xx5 * xx ; xx7 = xx6 * xx 
   
   vphi(0,0,0) = +  20     * xx7 -   70     * xx6 +   84   * xx5 -  35     * xx4                                    + 1    
   vphi(0,0,1) = + 140     * xx6 -  420     * xx5 +  420   * xx4 - 140     * xx3
   vphi(0,0,2) = + 840     * xx5 - 2100     * xx4 + 1680   * xx3 - 420     * xx2
   
   vphi(0,1,0) = +  10     * xx7 -   36     * xx6 +   45   * xx5 -  20     * xx4                           + 1 * xx 
   vphi(0,1,1) = +  70     * xx6 -  216     * xx5 +  225   * xx4 -  80     * xx3                           + 1    
   vphi(0,1,2) = + 420     * xx5 - 1080     * xx4 +  900   * xx3 - 240     * xx2
   
   vphi(0,2,0) = +   2     * xx7 -    7.5   * xx6 +   10   * xx5 -   5     * xx4               + 0.5 * xx2                       
   vphi(0,2,1) = +  14     * xx6 -   45     * xx5 +   50   * xx4 -  20     * xx3               + 1   * xx 
   vphi(0,2,2) = +  84     * xx5 -  225     * xx4 +  200   * xx3 -  60     * xx2               + 1  
   
   vphi(0,3,0) = +   1.0/6 * xx7 -    2.0/3 * xx6 +    1   * xx5 -   2.0/3 * xx4 + 1.0/6 * xx3                                    
   vphi(0,3,1) = +   7.0/6 * xx6 -    4     * xx5 +    5   * xx4 -   8.0/3 * xx3 + 3.0/6 * xx2
   vphi(0,3,2) = +   7     * xx5 -   20     * xx4 +   20   * xx3 -   8     * xx2 + 1     * xx 
   
   vphi(1,0,0) = -  20     * xx7 +   70     * xx6 -   84   * xx5 +  35     * xx4
   vphi(1,0,1) = - 140     * xx6 +  420     * xx5 -  420   * xx4 + 140     * xx3
   vphi(1,0,2) = - 840     * xx5 + 2100     * xx4 - 1680   * xx3 + 420     * xx2
   
   vphi(1,1,0) = +  10     * xx7 -   34     * xx6 +   39   * xx5 -  15     * xx4
   vphi(1,1,1) = +  70     * xx6 -  204     * xx5 +  195   * xx4 -  60     * xx3
   vphi(1,1,2) = + 420     * xx5 - 1020     * xx4 +  780   * xx3 - 180     * xx2
   
   vphi(1,2,0) = -   2     * xx7 +    6.5   * xx6 -    7   * xx5 +   2.5   * xx4
   vphi(1,2,1) = -  14     * xx6 +   39     * xx5 -   35   * xx4 +  10     * xx3
   vphi(1,2,2) = -  84     * xx5 +  195     * xx4 -  140   * xx3 +  30     * xx2
   
   vphi(1,3,0) = +   1.0/6 * xx7 -    0.5   * xx6 +    0.5 * xx5 -   1.0/6 * xx4
   vphi(1,3,1) = +   7.0/6 * xx6 -    3     * xx5 +    2.5 * xx4 -   4.0/6 * xx3
   vphi(1,3,2) = +   7     * xx5 -   15     * xx4 +   10   * xx3 -   2     * xx2
   
  case default

   FATAL(fe00ab : invalid Nofe,.true.)
   
  end select

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(fe00ab)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine fe00ab

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

