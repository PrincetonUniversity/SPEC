!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Lagrangian integration: locate high-order periodic field lines that approximate cantori.

!latex \item Cantori form important barriers that can severely restrict field line transport and thus anisotropic heat transport \cite{Hudson_Breslau_08}.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \end{enumerate} \subsubsection{Lagrangian integration} \begin{enumerate}

!latex \item Magnetic field lines are curves that extremize the action integral \cite{Cary_Littlejohn_83},
!latex       \be S \equiv \int_{\cal C} {\bf A}\cdot d{\bf l},
!latex       \ee
!latex       where ${\bf A}=A_\t \nabla \t + A_z \nabla \z$ is the magnetic vector potential, and 
!latex       $d{\bf l}\equiv d\s \, {\bf e_\s} + d\t \, {\bf e_\t} + d\z \, {\bf e_\z}$ is a line segment along a `trial' curve, ${\cal C}$,
!latex       which is described by $\s(\z)$ and $\t(\z)$ with $\z$ used to describe position along the curve.

!latex \item In the following, it is assumed that the vector potential is given,
!latex       \be   A_\t \equiv A_\t(\s,\t,\z) = \sum_j A_{\t,j}(\s)\cos(m_j \t - n_j \z), \\
!latex             A_\t \equiv A_\z(\s,\t,\z) = \sum_j A_{\z,j}(\s)\cos(m_j \t - n_j \z).
!latex       \ee
!latex \item The computational task is to construct extremizing periodic curves.

!latex \end{enumerate} \subsubsection{discretization of trial curve} \begin{enumerate}

!latex \item A practical discretization of ${\cal C}$ is given \cite{Hudson_06a}
!latex       \be \left. \begin{array}{ccl} \s & = & \s_i                                          \\
!latex                                     \t & = & \t_{i-1} + \dot \t_i \left(\z-\z_{i-1}\right)
!latex                  \end{array}
!latex           \right\} \; {\rm for } \; \z \in [\z_{i-1},\z_i],
!latex       \ee
!latex       where $\z_i\equiv i \, \Delta \z$, $\Delta \z\equiv \pi/N$ where $N$ is a resolution parameter, and $\dot \t_i \equiv (\t_{i}-\t_{i-1})/\Delta \z$.

!latex \item The curve is now described by the $\s_i$ and the $\t_i$.

!latex \end{enumerate} \subsubsection{periodicity constraint} \begin{enumerate}

!latex \item The periodicity constraint, $\t(\z+2\pi \, q) = \t(\z) + 2\pi \, p$, is enforced by the constraint
!latex       \be \t_{2qN} = \t_{0} + 2\pi \, p
!latex       \ee

!latex \item The degrees of freedom in the curve are thus $\s_i$ for $i=1,\dots,2qN$ and $\t_i$ for $i=0,\dots,2qN-1$.

!latex \end{enumerate} \subsubsection{piecewise action integral} \begin{enumerate}

!latex \item Using this representation for the trial curve, the action integral becomes
!latex       \be S = \sum_{i=1}^{2qN} S_{i}(\t_{i-1},\t_{i},\s_{i})
!latex       \ee
!latex       where
!latex       \be S_i & \equiv &        \int_{\z_{i-1}}^{\z_i} {\bf A}\cdot d{\bf l}                 \\
!latex               &    =   &        \int_{\z_{i-1}}^{\z_i} \left( A_\t \, \dot \t_i + A_\z \right) d\z \\
!latex               &    =   & a_1 \Delta \z + \sum_j a_j \lambda_{j,i},
!latex       \ee
!latex       where $a_j\equiv A_{\t,j}(s_i) \, \dot \t_i \, + A_{\z,j}(s_i)$,
!latex       $\lambda_{j,i} \equiv [ \sin(\alpha_{j,i})-\sin(\alpha_{j,i-1}) ] / (m_j \dot \t_i - n_j)$,
!latex       and $\alpha_{j,i}\equiv(m_j \t_i - n_j \z_i)$,
!latex       and where the summation over $j$ excludes the $(m_j,n_j)=(0,0)$ component.

!latex \end{enumerate} \subsubsection{conditions for extrema} \begin{enumerate}

!latex \item The action integral is extremized when 
!latex       \be \frac{\partial S}{\partial \s_i} &=& 0 \\
!latex           \frac{\partial S}{\partial \t_i} &=& 0.
!latex       \ee

!latex \item Assume that the $\t$ curve is given and the extremizing $\s$ curve is to be constructed.
!latex       We must solve
!latex       \be \frac{\partial S}{\partial s_i} = \frac{\partial S_i}{\partial s_i} = a_1^\prime \Delta \z + \sum_j a_j^\prime \lambda_{j,i}.
!latex      %\sum_j \left( A_{\t,j}^\prime \, \dot \t_i \, + A_{\z,j}^\prime \right) \lambda_{j,i} = 0,
!latex       \ee
!latex      %where $\lambda_{j,i} \equiv [ \sin(\alpha_{j,i})-\sin(\alpha_{j,i-1}) ] / (m_j \dot \t_i - n_j)$.
!latex \item Define $f(s_i) \equiv a_1^\prime \Delta \z + \sum_j a_j^\prime \lambda_{j,i}$.
!latex       A one-dimensional Newton method can be employed to find $f(s_i + \delta s_i) \approx f(s_i) + f^\prime(s_i) \, \delta s_i = 0$,
!latex       where $f^\prime(s_i) \equiv a_1^{\prime\prime} \Delta \z + \sum_j a_j^{\prime\prime} \lambda_{j,i}$.
!latex \item Note that the solution, $s_i$, depends only on $\t_{i-1}$ and $\t_{i}$, so that $s_i=s_i(\t_{i-1},\t_{i})$.
!latex       The action integral now becomes a function only of the $\t_i$.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine li00aa( lvol, nc )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half, one, pi, pi2

  use numerical, only : vsmall, small

  use fileunits, only : ounit

  use inputlist, only : Wli00aa, Npl, Mpqits

  use cputiming, only : Tli00aa

  use allglobal, only : myid, ncpu, cpus, interfacelabel, mn, im, in, pqorbit, mn, im, in, doAtmn, doAzmn, Dzeta
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS

  INTEGER, intent(in)    :: lvol, nc

  INTEGER                :: pp, qq, ii, tqN, Ldo, imn, its
  REAL                   :: lss, oargji, nargji, lji, ff, df, low, upp, harvest, snew
  CHARACTER              :: warn

  BEGIN(li00aa)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  pp = pqorbit(lvol,nc)%pq(1)
  qq = pqorbit(lvol,nc)%pq(2)

  low   = interfacelabel(lvol-1)
  upp   = interfacelabel(lvol  )

  tqN = 2*qq*Npl ; Ldo = 0 ! shorthand;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  cput = GETTIME  
  
  write(ounit,'("li00aa : "10x" : pp=",i3," ; qq="i4" ; tt="   99f9.4)') pp, qq, ( pqorbit(lvol,nc)%tt(ii), ii = 0, 2*qq*Npl )
  write(ounit,'("li00aa : "10x" : pp=",i3," ; qq="i4" ; ss="9x,99f9.4)') pp, qq, ( pqorbit(lvol,nc)%ss(ii), ii = 1, 2*qq*Npl )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! begin iteration loop;
  
  pqorbit(lvol,nc)%dt(1:tqN) = ( pqorbit(lvol,nc)%tt(1:tqN) - pqorbit(lvol,nc)%tt(0:tqN-1) ) / Dzeta
  

  do ii = 1, tqN
   
   its = 0 ; warn=" " ! iteration counter;
   
   do ; its = its+1
    
    lss = pqorbit(lvol,nc)%ss(ii) ! shorthand;
    
    call sa00aa( lvol, lss, Ldo ) ! interpolate Fourier harmonics of vector potential; returns doAtmn(0,0:2,1:mn) and doAzmn(0,0:2,1:mn);
    
    imn = 1
    ff = ( doAtmn(0,1,imn) + pqorbit(lvol,nc)%dt(ii) * doAzmn(0,1,imn) ) * Dzeta ! initialize summation;
    df = ( doAtmn(0,2,imn) + pqorbit(lvol,nc)%dt(ii) * doAzmn(0,2,imn) ) * Dzeta ! initialize summation;
    
    do imn = 2, mn ! exclude the (m,n)=(0,0) harmonic;
     
     oargji = im(imn) * pqorbit(lvol,nc)%tt(ii-1) - in(imn) * (ii-1)*Dzeta                 ! this does not depend on ss; need only be computed once;
     nargji = im(imn) * pqorbit(lvol,nc)%tt(ii  ) - in(imn) * (ii  )*Dzeta                 ! this does not depend on ss; need only be computed once;
     
     FATALMESS(li00aa,abs(im(imn) - pqorbit(lvol,nc)%dt(ii) * in(imn)).lt.small,zero denominator in Lagrangian integration)
     
     lji = ( sin(nargji) - sin(oargji) ) / ( im(imn) - pqorbit(lvol,nc)%dt(ii) * in(imn) ) ! this does not depend on ss; need only be computed once;
     
     ff = ff + ( doAtmn(0,1,imn) + pqorbit(lvol,nc)%dt(ii) * doAzmn(0,1,imn) ) * lji
     df = df + ( doAtmn(0,2,imn) + pqorbit(lvol,nc)%dt(ii) * doAzmn(0,2,imn) ) * lji
     
    enddo ! end of do imn;

    write(ounit,1000)myid,lvol,ii,its,lss,ff,df,warn
    
    FATALMESS(li00aa,abs(df).lt.small,zero denominator in Newton iterations)
    
    snew = lss - ff / df ! trial correction;
    
    if( snew.ge.low .and. snew.le.upp ) then ;                               pqorbit(lvol,nc)%ss(ii)= snew                    ; warn=" " ! accept Newton correction;
    else ;                                   ; call random_number(harvest) ; pqorbit(lvol,nc)%ss(ii)= low + (upp-low)*harvest ; warn="!" ! randomly restart;
    endif
    
    if( abs(ff).lt.vsmall .or. its.gt.Mpqits ) exit
    
   enddo ! end of do ; its = its + 1;
   
  enddo ! end of do ii;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(li00aa)
 
1000 format("li00aa : ", 10x ," : myid=",i3," ; lvol=",i3," ; ii="i6" ; its="i5" ; ss="f15.12" ; ff="f15.10" ; df="f15.10" ; "a1)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine li00aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
