!> \defgroup grp_build_matrices Build matrices
!>
!> \latexonly
!> \definecolor{Orange}{rgb}{1.0,0.5,0.0}
!> \definecolor{Cerulean}{rgb}{0.0,0.5,1.0}
!> \endlatexonly
!>
!> \file matrix.f90
!> \brief Constructs energy and helicity matrices that represent the Beltrami linear system.

!> **gauge conditions**
!>
!> <ul>
!> <li> In the \f$v\f$-th annulus, bounded by the \f$(v-1)\f$-th and \f$v\f$-th interfaces,
!>       a general covariant representation of the magnetic vector-potential is written
!>       \f{eqnarray}{ {\bf \bar A}=\bar A_{s} \nabla s + \bar A_{\theta} \nabla \theta + \bar A_{\zeta} \nabla \zeta eta. \f} </li>
!> <li> To this add \f$\nabla g(s,\theta,\zeta)\f$, where \f$g\f$ satisfies
!>       \f{eqnarray}{ \begin{array}{cccccccccccccccccccccccccccccccccccccc}
!>           \partial_s      g( s,\theta,\zeta) & = & - & \bar A_{s     }( s,\theta,\zeta) \\
!>           \partial_\theta g(-1,\theta,\zeta) & = & - & \bar A_{\theta}(-1,\theta,\zeta) \\
!>           \partial_\zeta  g(-1,     0,\zeta) & = & - & \bar A_{\zeta }(-1,     0,\zeta).
!>       \end{array}
!>      \f} </li>
!> <li> Then \f${\bf A}={\bf \bar A}+\nabla g\f$ is given by \f${\bf A}=A_{\theta}\nabla\theta+A_{\zeta}\nabla\zeta\f$ with
!>       \f{eqnarray}{
!>        A_{\theta}(-1,\theta,\zeta) &=& 0 \label{eq:Atgauge} \\
!>        A_{\zeta }(-1,     0,\zeta) &=& 0 \label{eq:Azgauge}
!>       \f} </li>
!> <li> This specifies the gauge: to see this, notice that no gauge term can be added without violating the conditions in Eqn.\f$(\ref{eq:Atgauge})\f$ or Eqn.\f$(\ref{eq:Azgauge})\f$. </li>
!> <li> Note that the gauge employed in each volume is distinct. </li>
!> </ul> 
!>
!> **boundary conditions**
!>
!> <ul>
!> <li> The magnetic field is 
!>       \f$\sqrt g \, {\bf B} = (\partial_\theta A_\zeta - \partial_\zeta A_\theta)\;{\bf e}_s - \partial_s A_\zeta \;{\bf e}_\theta + \partial_s A_\theta \;{\bf e}_\zeta\f$. </li>
!> <li> In the annular volumes, the condition that the field is tangential to the inner interface, \f$\sqrt g {\bf B}\cdot\nabla s=0\f$ at \f$s=-1\f$,
!>       gives \f$\partial_\theta A_\zeta - \partial_\zeta A_\theta = 0\f$.
!>       With the above condition on \f$A_\theta\f$ given in Eqn.\f$(\ref{eq:Atgauge})\f$, this gives \f$\partial_\theta A_\zeta=0\f$, which with Eqn.\f$(\ref{eq:Azgauge})\f$ gives
!>       \f{eqnarray}{ A_\zeta(-1,\theta,\zeta)=0.
!>       \f} </li>
!> <li> The condition at the outer interface, \f$s=+1\f$, is that the field is \f$\sqrt g \, {\bf B}\cdot\nabla s = \partial_\theta A_\zeta - \partial_\zeta A_\theta = b\f$,
!>       where \f$b\f$ is supplied by the user.
!>       For each of the plasma regions, \f$b=0\f$.
!>       For the vacuum region, generally \f$b\ne0\f$. </li>
!> </ul> 
!>
!> **enclosed fluxes**
!>
!> <ul>
!> <li> In the plasma regions, the enclosed fluxes must be constrained. </li>
!> <li> The toroidal and poloidal fluxes enclosed in each volume are determined using
!>       \f{eqnarray}{ \int_S {\bf B}\cdot{\bf ds}=\int_{\partial S}{\bf A}\cdot {\bf dl}.
!>       \f} </li>
!> </ul> 
!>
!> **Fourier-Chebyshev representation**
!>
!> <ul>
!> <li> The components of the vector potential, \f${\bf A}=A_\theta \nabla + A_\zeta \nabla \zeta eta\f$, are
!>      \f{eqnarray}{
!>        A_\theta(s,\theta,\zeta) &=& \sum_{i,l} {\color{red} A_{\theta,e,i,l}} \; {\overline T}_{l,i}(s) \cos\alpha_i + \sum_{i,l} {\color{Orange}  A_{\theta,o,i,l}} \; {\overline T}_{l,i}(s) \sin\alpha_i, \label{eq:At} \\
!>        A_\zeta( s,\theta,\zeta) &=& \sum_{i,l} {\color{blue}A_{\zeta, e,i,l}} \; {\overline T}_{l,i}(s) \cos\alpha_i + \sum_{i,l} {\color{Cerulean}A_{\zeta ,o,i,l}} \; {\overline T}_{l,i}(s) \sin\alpha_i, \label{eq:Az}
!>      \f}
!>      where \f${\overline T}_{l,i}(s) \equiv \bar s^{m_i/2} \, T_l(s)\f$, \f$T_l(s)\f$ is the Chebyshev polynomial, and \f$\alpha_j \equiv m_j\theta-n_j\zeta\f$.
!>      The regularity factor, \f$\bar s^{m_i/2}\f$, where \f$\bar s \equiv (1+s)/2\f$, is only included if there is a coordinate singularity in the domain
!>      (i.e. only in the innermost volume, and only in cylindrical and toroidal geometry.) </li>
!> <li> The magnetic field, \f$\sqrt g \, {\bf B} = \sqrt g B^s {\bf e}_s + \sqrt g B^\theta {\bf e}_\theta + \sqrt g B^\zeta {\bf e}_\zeta\f$, is
!>      \f{eqnarray}{
!>        \begin{array}{ccccrcrcrcrcccccccccccccccccccccccccccccccccccccccccccccccccccc}
!>        \sqrt g \, {\bf B} & = & {\bf e}_s      & \sum_{i,l} [ ( & - m_i {\color{blue}A_{\zeta, e,i,l}} & - & n_i {\color{red} A_{\theta,e,i,l}} & ) {\overline T}_{l,i}        \sin\alpha_i + ( & + m_i {\color{Cerulean}A_{\zeta ,o,i,l}} & + & n_i {\color{Orange}  A_{\theta,o,i,l}} & ) {\overline T}_{l,i}        \cos\alpha_i ] \\
!>                           & + & {\bf e}_\theta & \sum_{i,l} [ ( &                                      & - &     {\color{blue}A_{\zeta, e,i,l}} & ) {\overline T}_{l,i}^\prime \cos\alpha_i + ( &                                          & - &     {\color{Cerulean}A_{\zeta ,o,i,l}} & ) {\overline T}_{l,i}^\prime \sin\alpha_i ] \\
!>                           & + & {\bf e}_\zeta  & \sum_{i,l} [ ( &       {\color{red} A_{\theta,e,i,l}} &   &                                    & ) {\overline T}_{l,i}^\prime \cos\alpha_i + ( &       {\color{Orange}  A_{\theta,o,i,l}} &   &                                        & ) {\overline T}_{l,i}^\prime \sin\alpha_i ]
!>        \end{array}
!>      \f}
!> </li>
!> <li> The components of the velocity, \f${\bf v} \equiv v_s \nabla s + v_\theta \nabla \theta + v_\zeta \nabla \zeta eta\f$, are
!>       \f{eqnarray}{ v_s     (s,\theta,\zeta) &=& \sum_{i,l} {\color{red}  v_{     s,e,i,l}} \; {\overline T}_{l,i}(s) \cos\alpha_i + \sum_{i,l} {\color{Orange}   v_{     s,o,i,l}} \; {\overline T}_{l,i}(s) \sin\alpha_i, \\
!>                     v_\theta(s,\theta,\zeta) &=& \sum_{i,l} {\color{red}  v_{\theta,e,i,l}} \; {\overline T}_{l,i}(s) \cos\alpha_i + \sum_{i,l} {\color{Orange}   v_{\theta,o,i,l}} \; {\overline T}_{l,i}(s) \sin\alpha_i, \\
!>                     v_\zeta (s,\theta,\zeta) &=& \sum_{i,l} {\color{blue} v_{\zeta ,e,i,l}} \; {\overline T}_{l,i}(s) \cos\alpha_i + \sum_{i,l} {\color{Cerulean} v_{\zeta ,o,i,l}} \; {\overline T}_{l,i}(s) \sin\alpha_i.
!>       \f} </li>
!> </ul> 
!>
!> **constrained energy functional**
!>
!> <ul>
!>
!> <li> The constrained energy functional in each volume depends on the vector potential and the Lagrange multipliers, 
!> \f{eqnarray}{ {\cal F} \equiv
!> {\cal F}[{\color{red} A_{\theta,e,i,l}},{\color{blue}A_{\zeta, e,i,l}},{\color{Orange}  A_{\theta,o,i,l}},{\color{Cerulean}A_{\zeta ,o,i,l}},
!>  {\color{red} v_{s,e,i,l}},{\color{Orange} v_{s,o,i,l}},{\color{red} v_{\theta,e,i,l}},{\color{Orange} v_{\theta,o,i,l}},{\color{blue} v_{\zeta,e,i,l}},{\color{Cerulean} v_{\zeta,o,i,l}},\mu,a_i,b_i,c_i,d_i,e_i,f_i,g_1,h_1],
!> \f}
!>       and is given by:
!> \f{eqnarray}{ \begin{array}{cclcclcccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc}
!>     {\cal F}
!>     & \equiv & \displaystyle                        \int {\bf B} \cdot {\bf B} \, dv
!>        +                                            \int {\bf v} \cdot {\bf v} \, dv
!>        -       \displaystyle       \mu       \left[ \int {\bf A} \cdot {\bf B} \, dv - K \right]                                                                                                              \\
!>     &  + &   & \displaystyle \sum_{i=1} & a_i       & \displaystyle \left[  \sum_l {\color{red}     A_{\theta,e,i,l}} T_l(-1) - 0                                                                     \right] \\
!>     &  + &   & \displaystyle \sum_{i=1} & b_i       & \displaystyle \left[  \sum_l {\color{blue}    A_{\zeta, e,i,l}} T_l(-1) - 0                                                                     \right] \\
!>     &  + &   & \displaystyle \sum_{i=2} & c_i       & \displaystyle \left[  \sum_l {\color{Orange}  A_{\theta,o,i,l}} T_l(-1) - 0                                                                     \right] \\
!>     &  + &   & \displaystyle \sum_{i=2} & d_i       & \displaystyle \left[  \sum_l {\color{Cerulean}A_{\zeta ,o,i,l}} T_l(-1) - 0                                                                     \right] \\
!>     &  + &   & \displaystyle \sum_{i=2} & e_i       & \displaystyle \left[  \sum_l \left( - m_i {\color{blue}    A_{\zeta, e,i,l}} - n_i {\color{red}     A_{\theta,e,i,l}} \right) T_l(+1) - b_{s,i} \right] \\
!>     &  + &   & \displaystyle \sum_{i=2} & f_i       & \displaystyle \left[  \sum_l \left( + m_i {\color{Cerulean}A_{\zeta ,o,i,l}} + n_i {\color{Orange}  A_{\theta,o,i,l}} \right) T_l(+1) - b_{c,i} \right] \\
!>     &  + &   & \displaystyle            & g_1       & \displaystyle \left[  \sum_l {\color{red}     A_{\theta,e,1,l}} T_l(+1) - \Delta \psi_t \right]                                \\
!>     &  + &   & \displaystyle            & h_1       & \displaystyle \left[  \sum_l {\color{blue}    A_{\zeta, e,1,l}} T_l(+1) + \Delta \psi_p \right]                             
!>     \end{array}
!> \f}   
!>       where
!>       <ul>
!>       <li>\f$a_i\f$, \f$b_i\f$, \f$c_i\f$ and \f$d_i\f$ are Lagrange multipliers used to enforce the combined gauge and interface boundary condition
!>            on the inner interface, </li>
!>       <li>\f$e_i\f$ and \f$f_i\f$                       are Lagrange multipliers used to enforce the interface boundary condition on the outer interface, 
!>            namely \f$\sqrt g \, {\bf B}\cdot\nabla s = b\f$; and </li>
!>       <li>\f$g_1\f$ and \f$h_1\f$                       are Lagrange multipliers used to enforce the constraints on the enclosed fluxes. </li>
!>       </ul>
!> <li> In each plasma volume the boundary condition on the outer interface is \f$b=0\f$.
!> <li> In the vacuum volume (only for free-boundary), we may set \f$\mu=0\f$.
!>
!> </ul> 
!>
!> **derivatives of magnetic energy integrals**
!>
!> <ul>
!>
!> <li> The first  derivatives of \f$\int \! dv \; {\bf B} \! \cdot \! {\bf B}\f$ with respect to 
!> \f${\color{red} A_{\theta,e,i,l}}\f$, \f${\color{Orange}  A_{\theta,o,i,l}}\f$, \f${\color{blue}A_{\zeta, e,i,l}}\f$ and \f${\color{Cerulean}A_{\zeta ,o,i,l}}\f$ are
!> \f{eqnarray}{
!> \frac{\partial}{\partial {\color{red} A_{\theta,e,i,l}}} \int \!\! dv \; {\bf B} \cdot {\bf B} & = & 
!> 2 \int \!\! dv \; {\bf B} \cdot \frac{\partial {\bf B}}{\partial {\color{red} A_{\theta,e,i,l}}} =
!>      2 \int \!\! dv \; {\bf B} \cdot \left[ - n_i {\overline T}_{l,i}  \sin\alpha_i \, {\bf e}_s + {\overline T}_{l,i}^\prime \cos\alpha_i \, {\bf e}_\zeta \right] / \sqrt g \\
!> \frac{\partial}{\partial {\color{Orange}  A_{\theta,o,i,l}}} \int \!\! dv \; {\bf B} \cdot {\bf B} & = & 
!> 2 \int \!\! dv \; {\bf B} \cdot \frac{\partial {\bf B}}{\partial {\color{Orange}  A_{\theta,o,i,l}}} =
!>      2 \int \!\! dv \; {\bf B} \cdot \left[ + n_i {\overline T}_{l,i}  \cos\alpha_i \, {\bf e}_s + {\overline T}_{l,i}^\prime \sin\alpha_i \, {\bf e}_\zeta \right] / \sqrt g \\
!> \frac{\partial}{\partial {\color{blue}A_{\zeta, e,i,l}}} \int \!\! dv \; {\bf B} \cdot {\bf B} & = & 
!> 2 \int \!\! dv \; {\bf B} \cdot \frac{\partial {\bf B}}{\partial {\color{blue}A_{\zeta, e,i,l}}} =
!>      2 \int \!\! dv \; {\bf B} \cdot \left[ - m_i {\overline T}_{l,i}  \sin\alpha_i \, {\bf e}_s - {\overline T}_{l,i}^\prime \cos\alpha_i \, {\bf e}_\theta \right] / \sqrt g \\
!> \frac{\partial}{\partial {\color{Cerulean}A_{\zeta ,o,i,l}}} \int \!\! dv \; {\bf B} \cdot {\bf B} & = & 
!> 2 \int \!\! dv \; {\bf B} \cdot \frac{\partial {\bf B}}{\partial {\color{Cerulean}A_{\zeta ,o,i,l}}} =
!>      2 \int \!\! dv \; {\bf B} \cdot \left[ + m_i {\overline T}_{l,i}  \cos\alpha_i \, {\bf e}_s - {\overline T}_{l,i}^\prime \sin\alpha_i \, {\bf e}_\theta \right] / \sqrt g 
!> \f} </li>
!> <li> The second derivatives of \f$\int \! dv \; {\bf B} \! \cdot \! {\bf B}\f$ with respect to 
!>  \f${\color{red} A_{\theta,e,i,l}}\f$, \f${\color{Orange}  A_{\theta,o,i,l}}\f$, \f${\color{blue}A_{\zeta, e,i,l}}\f$ and \f${\color{Cerulean}A_{\zeta ,o,i,l}}\f$ are
!> \f{eqnarray}{ 
!> \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial {\color{red}     A_{\theta,e,j,p}}} \frac{\partial}{\partial {\color{red} A_{\theta,e,i,l}}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!> \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!>  ( + n_j n_i {\overline T}_{p,j}  {\overline T}_{l,i}  s_j s_i g_{ss} - n_j {\overline T}_{p,j}  {\overline T}_{l,i}^\prime s_j c_i g_{s\zeta}
!>    -     n_i {\overline T}_{l,i}  {\overline T}_{p,j}^\prime s_i c_j g_{s\zeta} +     {\overline T}_{p,j}^\prime {\overline T}_{l,i}^\prime c_j c_i g_{\zeta\zeta}) / \sqrt g^2\nonumber \\
!> \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial {\color{Orange}  A_{\theta,o,j,p}}} \frac{\partial}{\partial {\color{red} A_{\theta,e,i,l}}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!> \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!>  ( - n_j n_i {\overline T}_{p,j}  {\overline T}_{l,i}  c_j s_i g_{ss} + n_j {\overline T}_{p,j}  {\overline T}_{l,i}^\prime c_j c_i g_{s\zeta}
!>    -     n_i {\overline T}_{l,i}  {\overline T}_{p,j}^\prime s_i s_j g_{s\zeta} +     {\overline T}_{p,j}^\prime {\overline T}_{l,i}^\prime s_j c_i g_{\zeta\zeta}) / \sqrt g^2\nonumber \\
!> \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial {\color{blue}    A_{\zeta, e,j,p}}} \frac{\partial}{\partial {\color{red} A_{\theta,e,i,l}}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!> \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!>  ( + m_j n_i {\overline T}_{p,j}  {\overline T}_{l,i}  s_j s_i g_{ss} - m_j {\overline T}_{p,j}  {\overline T}_{l,i}^\prime s_j c_i g_{s\zeta}
!>    +     n_i {\overline T}_{l,i}  {\overline T}_{p,j}^\prime s_i c_j g_{s\theta} -     {\overline T}_{p,j}^\prime {\overline T}_{l,i}^\prime c_j c_i g_{\theta\zeta}) / \sqrt g^2\nonumber \\
!> \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial {\color{Cerulean}A_{\zeta ,o,j,p}}} \frac{\partial}{\partial {\color{red} A_{\theta,e,i,l}}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!> \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!>  ( - m_j n_i {\overline T}_{p,j}  {\overline T}_{l,i}  c_j s_i g_{ss} + m_j {\overline T}_{p,j}  {\overline T}_{l,i}^\prime c_j c_i g_{s\zeta}
!>    +     n_i {\overline T}_{l,i}  {\overline T}_{p,j}^\prime s_i s_j g_{s\theta} -     {\overline T}_{p,j}^\prime {\overline T}_{l,i}^\prime s_j c_i g_{\theta\zeta}) / \sqrt g^2\nonumber
!> \f} 
!> \f{eqnarray}{
!> \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial {\color{red}     A_{\theta,e,j,p}}} \frac{\partial}{\partial {\color{Orange}  A_{\theta,o,i,l}}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!> \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!>  ( - n_j n_i {\overline T}_{p,j}  {\overline T}_{l,i}  s_j c_i g_{ss} - n_j {\overline T}_{p,j}  {\overline T}_{l,i}^\prime s_j s_i g_{s\zeta}
!>    +     n_i {\overline T}_{l,i}  {\overline T}_{p,j}^\prime c_i c_j g_{s\zeta} +     {\overline T}_{p,j}^\prime {\overline T}_{l,i}^\prime c_j s_i g_{\zeta\zeta}) / \sqrt g^2\nonumber \\
!> \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial {\color{Orange}  A_{\theta,o,j,p}}} \frac{\partial}{\partial {\color{Orange}  A_{\theta,o,i,l}}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!> \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!>  ( + n_j n_i {\overline T}_{p,j}  {\overline T}_{l,i}  c_j c_i g_{ss} + n_j {\overline T}_{p,j}  {\overline T}_{l,i}^\prime c_j s_i g_{s\zeta}
!>    +     n_i {\overline T}_{l,i}  {\overline T}_{p,j}^\prime c_i s_j g_{s\zeta} +     {\overline T}_{p,j}^\prime {\overline T}_{l,i}^\prime s_j s_i g_{\zeta\zeta}) / \sqrt g^2\nonumber \\
!> \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial {\color{blue}    A_{\zeta, e,j,p}}} \frac{\partial}{\partial {\color{Orange}  A_{\theta,o,i,l}}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!> \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!>  ( - m_j n_i {\overline T}_{p,j}  {\overline T}_{l,i}  s_j c_i g_{ss} - m_j {\overline T}_{p,j}  {\overline T}_{l,i}^\prime s_j s_i g_{s\zeta}
!>    -     n_i {\overline T}_{l,i}  {\overline T}_{p,j}^\prime c_i c_j g_{s\theta} -     {\overline T}_{p,j}^\prime {\overline T}_{l,i}^\prime c_j s_i g_{\theta\zeta}) / \sqrt g^2\nonumber \\
!> \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial {\color{Cerulean}A_{\zeta ,o,j,p}}} \frac{\partial}{\partial {\color{Orange}  A_{\theta,o,i,l}}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!> \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!>  ( + m_j n_i {\overline T}_{p,j}  {\overline T}_{l,i}  c_j c_i g_{ss} + m_j {\overline T}_{p,j}  {\overline T}_{l,i}^\prime c_j s_i g_{s\zeta}
!>    -     n_i {\overline T}_{l,i}  {\overline T}_{p,j}^\prime c_i s_j g_{s\theta} -     {\overline T}_{p,j}^\prime {\overline T}_{l,i}^\prime s_j s_i g_{\theta\zeta}) / \sqrt g^2\nonumber   
!> \f} 
!> \f{eqnarray}{
!> \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial {\color{red}     A_{\theta,e,j,p}}} \frac{\partial}{\partial {\color{blue}A_{\zeta, e,i,l}}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!> \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!>  ( + n_j m_i {\overline T}_{p,j}  {\overline T}_{l,i}  s_j s_i g_{ss} + n_j {\overline T}_{p,j}  {\overline T}_{l,i}^\prime s_j c_i g_{s\theta}
!>    -     m_i {\overline T}_{l,i}  {\overline T}_{p,j}^\prime s_i c_j g_{s\zeta} -     {\overline T}_{p,j}^\prime {\overline T}_{l,i}^\prime c_j c_i g_{\theta\zeta}) / \sqrt g^2\nonumber \\
!> \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial {\color{Orange}  A_{\theta,o,j,p}}} \frac{\partial}{\partial {\color{blue}A_{\zeta, e,i,l}}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!> \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!>  ( - n_j m_i {\overline T}_{p,j}  {\overline T}_{l,i}  c_j s_i g_{ss} - n_j {\overline T}_{p,j}  {\overline T}_{l,i}^\prime c_j c_i g_{s\theta}
!>    -     m_i {\overline T}_{l,i}  {\overline T}_{p,j}^\prime s_i s_j g_{s\zeta} -     {\overline T}_{p,j}^\prime {\overline T}_{l,i}^\prime s_j c_i g_{\theta\zeta}) / \sqrt g^2\nonumber \\
!> \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial {\color{blue}    A_{\zeta, e,j,p}}} \frac{\partial}{\partial {\color{blue}A_{\zeta, e,i,l}}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!> \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!>  ( + m_j m_i {\overline T}_{p,j}  {\overline T}_{l,i}  s_j s_i g_{ss} + m_j {\overline T}_{p,j}  {\overline T}_{l,i}^\prime s_j c_i g_{s\theta}
!>    +     m_i {\overline T}_{l,i}  {\overline T}_{p,j}^\prime s_i c_j g_{s\theta} +     {\overline T}_{p,j}^\prime {\overline T}_{l,i}^\prime c_j c_i g_{\theta\theta}) / \sqrt g^2\nonumber \\
!> \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial {\color{Cerulean}A_{\zeta ,o,j,p}}} \frac{\partial}{\partial {\color{blue}A_{\zeta, e,i,l}}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!> \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!>  ( - m_j m_i {\overline T}_{p,j}  {\overline T}_{l,i}  c_j s_i g_{ss} - m_j {\overline T}_{p,j}  {\overline T}_{l,i}^\prime c_j c_i g_{s\theta}
!>    +     m_i {\overline T}_{l,i}  {\overline T}_{p,j}^\prime s_i s_j g_{s\theta} +     {\overline T}_{p,j}^\prime {\overline T}_{l,i}^\prime s_j c_i g_{\theta\theta}) / \sqrt g^2\nonumber
!> \f} 
!> \f{eqnarray}{
!> \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial {\color{red}     A_{\theta,e,j,p}}} \frac{\partial}{\partial {\color{Cerulean}A_{\zeta ,o,i,l}}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!> \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!>  ( - n_j m_i {\overline T}_{p,j}  {\overline T}_{l,i}  s_j c_i g_{ss} + n_j {\overline T}_{p,j}  {\overline T}_{l,i}^\prime s_j s_i g_{s\theta}
!>    +     m_i {\overline T}_{l,i}  {\overline T}_{p,j}^\prime c_i c_j g_{s\zeta} -     {\overline T}_{p,j}^\prime {\overline T}_{l,i}^\prime c_j s_i g_{\theta\zeta}) / \sqrt g^2\nonumber \\
!> \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial {\color{Orange}  A_{\theta,o,j,p}}} \frac{\partial}{\partial {\color{Cerulean}A_{\zeta ,o,i,l}}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!> \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!>  ( + n_j m_i {\overline T}_{p,j}  {\overline T}_{l,i}  c_j c_i g_{ss} - n_j {\overline T}_{p,j}  {\overline T}_{l,i}^\prime c_j s_i g_{s\theta}
!>    +     m_i {\overline T}_{l,i}  {\overline T}_{p,j}^\prime c_i s_j g_{s\zeta} -     {\overline T}_{p,j}^\prime {\overline T}_{l,i}^\prime s_j s_i g_{\theta\zeta}) / \sqrt g^2\nonumber \\
!> \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial {\color{blue}    A_{\zeta, e,j,p}}} \frac{\partial}{\partial {\color{Cerulean}A_{\zeta ,o,i,l}}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!> \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!>
!>  ( - m_j m_i {\overline T}_{p,j}  {\overline T}_{l,i}  s_j c_i g_{ss} + m_j {\overline T}_{p,j}  {\overline T}_{l,i}^\prime s_j s_i g_{s\theta}
!>    -     m_i {\overline T}_{l,i}  {\overline T}_{p,j}^\prime c_i c_j g_{s\theta} +     {\overline T}_{p,j}^\prime {\overline T}_{l,i}^\prime c_j s_i g_{\theta\theta}) / \sqrt g^2\nonumber \\
!> \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial {\color{Cerulean}A_{\zeta ,o,j,p}}} \frac{\partial}{\partial {\color{Cerulean}A_{\zeta ,o,i,l}}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!> \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!>
!>  ( + m_j m_i {\overline T}_{p,j}  {\overline T}_{l,i}  c_j c_i g_{ss} - m_j {\overline T}_{p,j}  {\overline T}_{l,i}^\prime c_j s_i g_{s\theta}
!>    -     m_i {\overline T}_{l,i}  {\overline T}_{p,j}^\prime c_i s_j g_{s\theta} +     {\overline T}_{p,j}^\prime {\overline T}_{l,i}^\prime s_j s_i g_{\theta\theta}) / \sqrt g^2\nonumber
!> \f} </li>
!>
!> </ul>
!>
!> **derivatives of helicity        integrals**
!>
!> <ul>
!>
!> <li> The first  derivatives of \f$\int \! dv \; {\bf A} \cdot{\bf B}\f$ with respect to 
!> \f${\color{red} A_{\theta,e,i,l}}\f$, \f${\color{Orange}  A_{\theta,o,i,l}}\f$, \f${\color{blue}A_{\zeta, e,i,l}}\f$ and \f${\color{Cerulean}A_{\zeta ,o,i,l}}\f$ are
!> \f{eqnarray}{ 
!> \frac{\partial}{\partial {\color{red} A_{\theta,e,i,l}}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &
!> \int \!\! dv \; \left( \frac{\partial {\bf A}}{\partial {\color{red} A_{\theta,e,i,l}}} \cdot {\bf B} + {\bf A} \cdot \frac{\partial {\bf B}}{\partial {\color{red} A_{\theta,e,i,l}}} \right) =
!> \int \!\! dv \; ( {\overline T}_{l,i} \cos\alpha_i \nabla \theta \cdot {\bf B} + {\bf A} \cdot {\overline T}_{l,i}^\prime \cos\alpha_i \, {\bf e}_\zeta / \sqrt g ) \\
!> \frac{\partial}{\partial {\color{Orange}  A_{\theta,o,i,l}}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &
!> \int \!\! dv \; \left( \frac{\partial {\bf A}}{\partial {\color{Orange}  A_{\theta,o,i,l}}} \cdot {\bf B} + {\bf A} \cdot \frac{\partial {\bf B}}{\partial {\color{Orange}  A_{\theta,o,i,l}}} \right) =
!> \int \!\! dv \; ( {\overline T}_{l,i} \sin\alpha_i \nabla \theta \cdot {\bf B} + {\bf A} \cdot {\overline T}_{l,i}^\prime \sin\alpha_i \, {\bf e}_\zeta / \sqrt g ) \\
!> \frac{\partial}{\partial {\color{blue}A_{\zeta, e,i,l}}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &
!> \int \!\! dv \; \left( \frac{\partial {\bf A}}{\partial {\color{blue}A_{\zeta, e,i,l}}} \cdot {\bf B} + {\bf A} \cdot \frac{\partial {\bf B}}{\partial {\color{blue}A_{\zeta, e,i,l}}} \right) =
!> \int \!\! dv \; ( {\overline T}_{l,i} \cos\alpha_i \nabla \zeta  \cdot {\bf B} - {\bf A} \cdot {\overline T}_{l,i}^\prime \cos\alpha_i \, {\bf e}_\theta / \sqrt g ) \\
!> \frac{\partial}{\partial {\color{Cerulean}A_{\zeta ,o,i,l}}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &
!> \int \!\! dv \; \left( \frac{\partial {\bf A}}{\partial {\color{Cerulean}A_{\zeta ,o,i,l}}} \cdot {\bf B} + {\bf A} \cdot \frac{\partial {\bf B}}{\partial {\color{Cerulean}A_{\zeta ,o,i,l}}} \right) =
!> \int \!\! dv \; ( {\overline T}_{l,i} \sin\alpha_i \nabla \zeta  \cdot {\bf B} - {\bf A} \cdot {\overline T}_{l,i}^\prime \sin\alpha_i \, {\bf e}_\theta / \sqrt g )
!> \f} </li>
!> <li> Note that in the above expressions, \f${\bf A}\cdot{\bf e}_s=0\f$ has been used. </li>
!> <li> The second derivatives of \f$\int \! dv \; {\bf A} \cdot{\bf B}\f$ with respect to 
!> \f${\color{red} A_{\theta,e,i,l}}\f$, \f${\color{Orange}  A_{\theta,o,i,l}}\f$, \f${\color{blue}A_{\zeta, e,i,l}}\f$ and \f${\color{Cerulean}A_{\zeta ,o,i,l}}\f$ are
!> \f{eqnarray}{
!> \frac{\partial}{\partial {\color{red}     A_{\theta,e,j,p}}} \frac{\partial}{\partial {\color{red} A_{\theta,e,i,l}}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &   \int \!\! dv \; \left[
!> \cancel{
!> + {\overline T}_{l,i} \cos\alpha_i \nabla \theta \cdot {\overline T}_{p,j}^\prime \cos\alpha_j \, {\bf e}_\zeta}
!> \cancel{
!> + {\overline T}_{p,j} \cos\alpha_j \nabla \theta \cdot {\overline T}_{l,i}^\prime \cos\alpha_i \, {\bf e}_\zeta}
!> \right] / \sqrt g \\
!> \frac{\partial}{\partial {\color{Orange}  A_{\theta,o,j,p}}} \frac{\partial}{\partial {\color{red} A_{\theta,e,i,l}}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &   \int \!\! dv \;  \left[
!> \cancel{
!> + {\overline T}_{l,i} \cos\alpha_i \nabla \theta \cdot {\overline T}_{p,j}^\prime \sin\alpha_j \, {\bf e}_\zeta}
!> \cancel{
!> + {\overline T}_{p,j} \sin\alpha_j \nabla \theta \cdot {\overline T}_{l,i}^\prime \cos\alpha_i \, {\bf e}_\zeta}
!> \right] / \sqrt g \\
!> \frac{\partial}{\partial {\color{blue}    A_{\zeta, e,j,p}}} \frac{\partial}{\partial {\color{red} A_{\theta,e,i,l}}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &   \int \!\! dv \;  \left[
!> - {\overline T}_{l,i} \cos\alpha_i \nabla \theta \cdot {\overline T}_{p,j}^\prime \cos\alpha_j \, {\bf e}_\theta + {\overline T}_{p,j} \cos\alpha_j \nabla \zeta  \cdot {\overline T}_{l,i}^\prime \cos\alpha_i \, {\bf e}_\zeta \right] / \sqrt g \\
!> \frac{\partial}{\partial {\color{Cerulean}A_{\zeta ,o,j,p}}} \frac{\partial}{\partial {\color{red} A_{\theta,e,i,l}}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &   \int \!\! dv \;  \left[
!> - {\overline T}_{l,i} \cos\alpha_i \nabla \theta \cdot {\overline T}_{p,j}^\prime \sin\alpha_j \, {\bf e}_\theta + {\overline T}_{p,j} \sin\alpha_j \nabla \zeta  \cdot {\overline T}_{l,i}^\prime \cos\alpha_i \, {\bf e}_\zeta \right] / \sqrt g \\
!> \nonumber \\
!> \frac{\partial}{\partial {\color{red}     A_{\theta,e,j,p}}} \frac{\partial}{\partial {\color{Orange}  A_{\theta,o,i,l}}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &   \int \!\! dv \;  \left[
!> \cancel{
!> + {\overline T}_{l,i} \sin\alpha_i \nabla \theta \cdot {\overline T}_{p,j}^\prime \cos\alpha_j \, {\bf e}_\zeta}
!> \cancel{
!> + {\overline T}_{p,j} \cos\alpha_j \nabla \theta \cdot {\overline T}_{l,i}^\prime \sin\alpha_i \, {\bf e}_\zeta}
!> \right] / \sqrt g \\
!> \frac{\partial}{\partial {\color{Orange}  A_{\theta,o,j,p}}} \frac{\partial}{\partial {\color{Orange}  A_{\theta,o,i,l}}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &   \int \!\! dv \;  \left[
!> \cancel{
!> + {\overline T}_{l,i} \sin\alpha_i \nabla \theta \cdot {\overline T}_{p,j}^\prime \sin\alpha_j \, {\bf e}_\zeta}
!> \cancel{
!> + {\overline T}_{p,j} \sin\alpha_j \nabla \theta \cdot {\overline T}_{l,i}^\prime \sin\alpha_i \, {\bf e}_\zeta}
!> \right] / \sqrt g \\
!> \frac{\partial}{\partial {\color{blue}    A_{\zeta, e,j,p}}} \frac{\partial}{\partial {\color{Orange}  A_{\theta,o,i,l}}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &   \int \!\! dv \;  \left[
!> - {\overline T}_{l,i} \sin\alpha_i \nabla \theta \cdot {\overline T}_{p,j}^\prime \cos\alpha_j \, {\bf e}_\theta + {\overline T}_{p,j} \cos\alpha_j \nabla \zeta  \cdot {\overline T}_{l,i}^\prime \sin\alpha_i \, {\bf e}_\zeta \right] / \sqrt g \\
!> \frac{\partial}{\partial {\color{Cerulean}A_{\zeta ,o,j,p}}} \frac{\partial}{\partial {\color{Orange}  A_{\theta,o,i,l}}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &   \int \!\! dv \;  \left[
!> - {\overline T}_{l,i} \sin\alpha_i \nabla \theta \cdot {\overline T}_{p,j}^\prime \sin\alpha_j \, {\bf e}_\theta + {\overline T}_{p,j} \sin\alpha_j \nabla \zeta  \cdot {\overline T}_{l,i}^\prime \sin\alpha_i \, {\bf e}_\zeta \right] / \sqrt g \\
!> \nonumber \\
!> \frac{\partial}{\partial {\color{red}     A_{\theta,e,j,p}}} \frac{\partial}{\partial {\color{blue}A_{\zeta, e,i,l}}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &   \int \!\! dv \;  \left[
!> + {\overline T}_{l,i} \cos\alpha_i \nabla \zeta  \cdot {\overline T}_{p,j}^\prime \cos\alpha_j \, {\bf e}_\zeta - {\overline T}_{p,j} \cos\alpha_j \nabla \theta \cdot {\overline T}_{l,i}^\prime \cos\alpha_i \, {\bf e}_\theta \right] / \sqrt g \\
!> \frac{\partial}{\partial {\color{Orange}  A_{\theta,o,j,p}}} \frac{\partial}{\partial {\color{blue}A_{\zeta, e,i,l}}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &   \int \!\! dv \;  \left[
!> + {\overline T}_{l,i} \cos\alpha_i \nabla \zeta  \cdot {\overline T}_{p,j}^\prime \sin\alpha_j \, {\bf e}_\zeta - {\overline T}_{p,j} \sin\alpha_j \nabla \theta \cdot {\overline T}_{l,i}^\prime \cos\alpha_i \, {\bf e}_\theta \right] / \sqrt g \\
!> \frac{\partial}{\partial {\color{blue}    A_{\zeta, e,j,p}}} \frac{\partial}{\partial {\color{blue}A_{\zeta, e,i,l}}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &   \int \!\! dv \;  \left[
!> \cancel{
!> - {\overline T}_{l,i} \cos\alpha_i \nabla \zeta  \cdot {\overline T}_{p,j}^\prime \cos\alpha_j \, {\bf e}_\theta}
!> \cancel{
!> - {\overline T}_{p,j} \cos\alpha_j \nabla \zeta  \cdot {\overline T}_{l,i}^\prime \cos\alpha_i \, {\bf e}_\theta}
!> \right] / \sqrt g \\
!> \frac{\partial}{\partial {\color{Cerulean}A_{\zeta ,o,j,p}}} \frac{\partial}{\partial {\color{blue}A_{\zeta, e,i,l}}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &   \int \!\! dv \;  \left[
!> \cancel{
!> - {\overline T}_{l,i} \cos\alpha_i \nabla \zeta  \cdot {\overline T}_{p,j}^\prime \sin\alpha_j \, {\bf e}_\theta}
!> \cancel{
!> - {\overline T}_{p,j} \sin\alpha_j \nabla \zeta  \cdot {\overline T}_{l,i}^\prime \cos\alpha_i \, {\bf e}_\theta}
!> \right] / \sqrt g \\
!> \nonumber \\
!> \frac{\partial}{\partial {\color{red}     A_{\theta,e,j,p}}} \frac{\partial}{\partial {\color{Cerulean}A_{\zeta ,o,i,l}}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &   \int \!\! dv \;  \left[
!> + {\overline T}_{l,i} \sin\alpha_i \nabla \zeta  \cdot {\overline T}_{p,j}^\prime \cos\alpha_j \, {\bf e}_\zeta - {\overline T}_{p,j} \cos\alpha_j \nabla \theta \cdot {\overline T}_{l,i}^\prime \sin\alpha_i \, {\bf e}_\theta \right] / \sqrt g \\
!> \frac{\partial}{\partial {\color{Orange}  A_{\theta,o,j,p}}} \frac{\partial}{\partial {\color{Cerulean}A_{\zeta ,o,i,l}}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &   \int \!\! dv \;  \left[
!> + {\overline T}_{l,i} \sin\alpha_i \nabla \zeta  \cdot {\overline T}_{p,j}^\prime \sin\alpha_j \, {\bf e}_\zeta - {\overline T}_{p,j} \sin\alpha_j \nabla \theta \cdot {\overline T}_{l,i}^\prime \sin\alpha_i \, {\bf e}_\theta \right] / \sqrt g \\
!> \frac{\partial}{\partial {\color{blue}    A_{\zeta, e,j,p}}} \frac{\partial}{\partial {\color{Cerulean}A_{\zeta ,o,i,l}}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &   \int \!\! dv \;  \left[
!> \cancel{
!> - {\overline T}_{l,i} \sin\alpha_i \nabla \zeta  \cdot {\overline T}_{p,j}^\prime \cos\alpha_j \, {\bf e}_\theta}
!> \cancel{
!> - {\overline T}_{p,j} \cos\alpha_j \nabla \zeta  \cdot {\overline T}_{l,i}^\prime \sin\alpha_i \, {\bf e}_\theta}
!> \right] / \sqrt g \\
!> \frac{\partial}{\partial {\color{Cerulean}A_{\zeta ,o,j,p}}} \frac{\partial}{\partial {\color{Cerulean}A_{\zeta ,o,i,l}}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &   \int \!\! dv \;  \left[
!> \cancel{
!> - {\overline T}_{l,i} \sin\alpha_i \nabla \zeta  \cdot {\overline T}_{p,j}^\prime \sin\alpha_j \, {\bf e}_\theta}
!> \cancel{
!> - {\overline T}_{p,j} \sin\alpha_j \nabla \zeta  \cdot {\overline T}_{l,i}^\prime \sin\alpha_i \, {\bf e}_\theta}
!> \right] / \sqrt g
!> \f} </li>
!> <li> In these expressions the terms  \f$\nabla \theta \cdot {\bf e}_\theta =        \nabla \zeta \cdot {\bf e}_\zeta =1 \f$, 
!>                                  and \f$\cancel{\nabla \theta \cdot {\bf e}_\zeta} =\cancel{\nabla \zeta \cdot {\bf e}_\theta}=0\f$ 
!>       have been included to show the structure of the derivation. </li>
!> </ul>
!>
!> **derivatives of kinetic energy integrals**
!>
!> <ul>
!> <li> The first derivatives of \f$\int \! dv \; v^2\f$ with respect to \f${\color{red}  v_{     s,e,i,l}}\f$ etc. are
!> \f{eqnarray}{ \frac{\partial}{\partial {\color{red}      v_{     s,e,i,l}}} \int \!\! dv \; {\bf v} \cdot {\bf v} & = & 2 \int \!\! dv \; {\bf v} \cdot {\overline T}_{l,i} \cos\alpha_i \nabla s \\
!>               \frac{\partial}{\partial {\color{Orange}   v_{     s,o,i,l}}} \int \!\! dv \; {\bf v} \cdot {\bf v} & = & 2 \int \!\! dv \; {\bf v} \cdot {\overline T}_{l,i} \sin\alpha_i \nabla s \\
!>               \frac{\partial}{\partial {\color{red}      v_{\theta,e,i,l}}} \int \!\! dv \; {\bf v} \cdot {\bf v} & = & 2 \int \!\! dv \; {\bf v} \cdot {\overline T}_{l,i} \cos\alpha_i \nabla \theta \\
!>               \frac{\partial}{\partial {\color{Orange}   v_{\theta,o,i,l}}} \int \!\! dv \; {\bf v} \cdot {\bf v} & = & 2 \int \!\! dv \; {\bf v} \cdot {\overline T}_{l,i} \sin\alpha_i \nabla \theta \\
!>               \frac{\partial}{\partial {\color{blue}     v_{\zeta ,e,i,l}}} \int \!\! dv \; {\bf v} \cdot {\bf v} & = & 2 \int \!\! dv \; {\bf v} \cdot {\overline T}_{l,i} \cos\alpha_i \nabla \zeta  \\
!>               \frac{\partial}{\partial {\color{Cerulean} v_{\zeta ,o,i,l}}} \int \!\! dv \; {\bf v} \cdot {\bf v} & = & 2 \int \!\! dv \; {\bf v} \cdot {\overline T}_{l,i} \sin\alpha_i \nabla \zeta  \\
!> \f} </li>
!> </ul>
!>
!> **calculation of volume-integrated basis-function-weighted metric information**
!>
!> <ul>
!> <li> The required geometric information is calculated in ma00aa(). </li>as
!> </ul>

!> \brief Constructs energy and helicity matrices that represent the Beltrami linear system.
!> \ingroup grp_build_matrices
!> @param[in] lvol
!> @param[in] mn
!> @param[in] lrad
subroutine matrix( lvol, mn, lrad )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, one
  
  use numerical, only : small
  
  use fileunits, only : ounit
  
  use inputlist, only : Wmacros, Wmatrix
  
  use cputiming, only : Tmatrix
  
  use allglobal, only : ncpu, myid, cpus, &
                        YESstellsym, NOTstellsym, &
                        im, in, &
                        NAdof, &
                        dMA, dMD, dMB, dMG, &
                        Ate, Ato, Aze, Azo, &
                        iVns, iBns, iVnc, iBnc, &
                        Lma, Lmb, Lmc, Lmd, Lme, Lmf, Lmg, Lmh, &
                        Lcoordinatesingularity, TT, &
                        DToocc, DToocs, DToosc, DTooss, &
                        TTsscc, TTsscs, TTsssc, TTssss, &
                        TDstcc, TDstcs, TDstsc, TDstss, &
                        TDszcc, TDszcs, TDszsc, TDszss, &
                        DDttcc, DDttcs, DDttsc, DDttss, &
                        DDtzcc, DDtzcs, DDtzsc, DDtzss, &
                        DDzzcc, DDzzcs, DDzzsc, DDzzss
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  INTEGER, intent(in)  :: lvol, mn, lrad
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  INTEGER              :: NN, ii, jj, ll, kk, pp, mi, ni, mj, nj, mimj, minj, nimj, ninj, mjmi, mjni, njmi, njni, id, jd
  
  REAL                 :: Wtete, Wteto, Wtote, Wtoto
  REAL                 :: Wteze, Wtezo, Wtoze, Wtozo
  REAL                 :: Wzete, Wzeto, Wzote, Wzoto
  REAL                 :: Wzeze, Wzezo, Wzoze, Wzozo
  
  REAL                 :: Htete, Hteto, Htote, Htoto
  REAL                 :: Hteze, Htezo, Htoze, Htozo
  REAL                 :: Hzete, Hzeto, Hzote, Hzoto
  REAL                 :: Hzeze, Hzezo, Hzoze, Hzozo
  
  BEGIN(matrix)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  FATAL( matrix, .not.allocated(dMA), error )
  FATAL( matrix, .not.allocated(dMD), error )
  FATAL( matrix, .not.allocated(dMB), error )
  FATAL( matrix, .not.allocated(dMG), error )  
#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  NN = NAdof(lvol) ! shorthand;
  
  dMA(0:NN,0:NN) = zero
  dMD(0:NN,0:NN) = zero
  dMB(0:NN,1: 2) = zero
  dMG(0:NN     ) = zero
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( YESstellsym ) then
   
   do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
    
    do jj = 1, mn ; mj = im(jj) ; nj = in(jj) ; mimj = mi * mj ; minj = mi * nj ; nimj = ni * mj ; ninj = ni * nj
     
     do ll = 0, lrad
      
      do pp = 0, lrad
       
       Wtete = + 2 * ninj * TTssss(ll,pp,ii,jj) - 2 * ni      * TDszsc(ll,pp,ii,jj) - 2      * nj * TDszsc(pp,ll,jj,ii) + 2 * DDzzcc(ll,pp,ii,jj)
       Wzete = + 2 * nimj * TTssss(ll,pp,ii,jj) + 2 * ni      * TDstsc(ll,pp,ii,jj) - 2      * mj * TDszsc(pp,ll,jj,ii) - 2 * DDtzcc(pp,ll,jj,ii)
       Wteze = + 2 * minj * TTssss(ll,pp,ii,jj) + 2      * nj * TDstsc(pp,ll,jj,ii) - 2 * mi      * TDszsc(ll,pp,ii,jj) - 2 * DDtzcc(ll,pp,ii,jj)
       Wzeze = + 2 * mimj * TTssss(ll,pp,ii,jj) + 2 * mi      * TDstsc(ll,pp,ii,jj) + 2      * mj * TDstsc(pp,ll,jj,ii) + 2 * DDttcc(ll,pp,ii,jj)
       
       Htete =   zero
       Hzete = - DToocc(pp,ll,jj,ii) + DToocc(ll,pp,ii,jj)
       Hteze = + DToocc(pp,ll,jj,ii) - DToocc(ll,pp,ii,jj)
       Hzeze =   zero 
       
       id = Ate(lvol,0,ii)%i(ll) ; jd = Ate(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wtete ; dMD(id,jd) = Htete
       ;                         ; jd = Aze(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wzete ; dMD(id,jd) = Hzete
       id = Aze(lvol,0,ii)%i(ll) ; jd = Ate(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wteze ; dMD(id,jd) = Hteze
       ;                         ; jd = Aze(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wzeze ; dMD(id,jd) = Hzeze
       
      enddo ! end of do pp ;
      
     enddo ! end of do ll ;
     
    enddo ! end of do jj ;
    
    if( Lcoordinatesingularity .and. ii.eq.1 ) then ; kk = 1
    else                                            ; kk = 0
    endif
    
    do ll = 0, lrad ! Chebyshev polynomial ;
     
     ;                  ; id = Ate(lvol,0,ii)%i(ll) ; jd = Lma(lvol,  ii)       ; dMA(id,jd) = +      TT(ll, 0,0)
     ;                  ; id = Aze(lvol,0,ii)%i(ll) ; jd = Lmb(lvol,  ii)       ; dMA(id,jd) = +      TT(ll,kk,0) ! check coordinate singularity ;
     if( ii.gt.1 ) then ; id = Ate(lvol,0,ii)%i(ll) ; jd = Lme(lvol,  ii)       ; dMA(id,jd) = - ni * TT(ll, 1,0)
      ;                 ; id = Aze(lvol,0,ii)%i(ll) ; jd = Lme(lvol,  ii)       ; dMA(id,jd) = - mi * TT(ll, 1,0)
     else               ; id = Ate(lvol,0,ii)%i(ll) ; jd = Lmg(lvol,  ii)       ; dMA(id,jd) = +      TT(ll, 1,0)
      ;                 ; id = Aze(lvol,0,ii)%i(ll) ; jd = Lmh(lvol,  ii)       ; dMA(id,jd) = +      TT(ll, 1,0)
     endif
     
    enddo ! end of do ll ;
    
    do pp = 0, lrad     ; id = Lma(lvol,  ii)       ; jd = Ate(lvol,0,ii)%i(pp) ; dMA(id,jd) = +      TT(pp, 0,0)
     ;                  ; id = Lmb(lvol,  ii)       ; jd = Aze(lvol,0,ii)%i(pp) ; dMA(id,jd) = +      TT(pp,kk,0) ! check coordinate singularity ;
     if( ii.gt.1 ) then ; id = Lme(lvol,  ii)       ; jd = Ate(lvol,0,ii)%i(pp) ; dMA(id,jd) = - ni * TT(pp, 1,0)
      ;                 ; id = Lme(lvol,  ii)       ; jd = Aze(lvol,0,ii)%i(pp) ; dMA(id,jd) = - mi * TT(pp, 1,0)
     else               ; id = Lmg(lvol,  ii)       ; jd = Ate(lvol,0,ii)%i(pp) ; dMA(id,jd) = +      TT(pp, 1,0)
      ;                 ; id = Lmh(lvol,  ii)       ; jd = Aze(lvol,0,ii)%i(pp) ; dMA(id,jd) = +      TT(pp, 1,0)
     endif
    enddo ! end of do pp ;
    
    ;if( ii.gt.1 ) then ; id = Lme(lvol,  ii)       ;                           ; dMG(id   ) = - ( iVns(ii) + iBns(ii) )
    ;else               ; id = Lmg(lvol,  ii)       ;                           ; dMB(id, 1) = -       one
!   ;                   ; id = Lmh(lvol,  ii)       ;                           ; dMB(id, 2) = -       one ! to be deleted;
    ;                   ; id = Lmh(lvol,  ii)       ;                           ; dMB(id, 2) = +       one ! changed sign;
    ;endif
    
   enddo ! end of do ii ;
   
  else ! NOTstellsym ;
   
   do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
    
    do jj = 1, mn ; mj = im(jj) ; nj = in(jj) ; mjmi = mi * mj ; mjni = ni * mj ; njmi = mi * nj ; njni = ni * nj
     
     do ll = 0, lrad
      
      do pp = 0, lrad
       
       Wtete = 2 * ( + njni * TTssss(pp,ll,jj,ii) - nj * TDszsc(pp,ll,jj,ii) - ni * TDszsc(ll,pp,ii,jj) + DDzzcc(pp,ll,jj,ii) )
       Wtote = 2 * ( - njni * TTsscs(pp,ll,jj,ii) + nj * TDszcc(pp,ll,jj,ii) - ni * TDszss(ll,pp,ii,jj) + DDzzsc(pp,ll,jj,ii) )
       Wzete = 2 * ( + mjni * TTssss(pp,ll,jj,ii) - mj * TDszsc(pp,ll,jj,ii) + ni * TDstsc(ll,pp,ii,jj) - DDtzcc(pp,ll,jj,ii) )
       Wzote = 2 * ( - mjni * TTsscs(pp,ll,jj,ii) + mj * TDszcc(pp,ll,jj,ii) + ni * TDstss(ll,pp,ii,jj) - DDtzsc(pp,ll,jj,ii) )
       
       Wteto = 2 * ( - njni * TTsssc(pp,ll,jj,ii) - nj * TDszss(pp,ll,jj,ii) + ni * TDszcc(ll,pp,ii,jj) + DDzzcs(pp,ll,jj,ii) )
       Wtoto = 2 * ( + njni * TTsscc(pp,ll,jj,ii) + nj * TDszcs(pp,ll,jj,ii) + ni * TDszcs(ll,pp,ii,jj) + DDzzss(pp,ll,jj,ii) )
       Wzeto = 2 * ( - mjni * TTsssc(pp,ll,jj,ii) - mj * TDszss(pp,ll,jj,ii) - ni * TDstcc(ll,pp,ii,jj) - DDtzcs(pp,ll,jj,ii) )
       Wzoto = 2 * ( + mjni * TTsscc(pp,ll,jj,ii) + mj * TDszcs(pp,ll,jj,ii) - ni * TDstcs(ll,pp,ii,jj) - DDtzss(pp,ll,jj,ii) )
       
       Wteze = 2 * ( + njmi * TTssss(pp,ll,jj,ii) + nj * TDstsc(pp,ll,jj,ii) - mi * TDszsc(ll,pp,ii,jj) - DDtzcc(pp,ll,jj,ii) )
       Wtoze = 2 * ( - njmi * TTsscs(pp,ll,jj,ii) - nj * TDstcc(pp,ll,jj,ii) - mi * TDszss(ll,pp,ii,jj) - DDtzsc(pp,ll,jj,ii) )
       Wzeze = 2 * ( + mjmi * TTssss(pp,ll,jj,ii) + mj * TDstsc(pp,ll,jj,ii) + mi * TDstsc(ll,pp,ii,jj) + DDttcc(pp,ll,jj,ii) )
       Wzoze = 2 * ( - mjmi * TTsscs(pp,ll,jj,ii) - mj * TDstcc(pp,ll,jj,ii) + mi * TDstss(ll,pp,ii,jj) + DDttsc(pp,ll,jj,ii) )
       
       Wtezo = 2 * ( - njmi * TTsssc(pp,ll,jj,ii) + nj * TDstss(pp,ll,jj,ii) + mi * TDszcc(ll,pp,ii,jj) - DDtzcs(pp,ll,jj,ii) )
       Wtozo = 2 * ( + njmi * TTsscc(pp,ll,jj,ii) - nj * TDstcs(pp,ll,jj,ii) + mi * TDszcs(ll,pp,ii,jj) - DDtzss(pp,ll,jj,ii) )
       Wzezo = 2 * ( - mjmi * TTsssc(pp,ll,jj,ii) + mj * TDstss(pp,ll,jj,ii) - mi * TDstcc(ll,pp,ii,jj) + DDttcs(pp,ll,jj,ii) )
       Wzozo = 2 * ( + mjmi * TTsscc(pp,ll,jj,ii) - mj * TDstcs(pp,ll,jj,ii) - mi * TDstcs(ll,pp,ii,jj) + DDttss(pp,ll,jj,ii) )
       
       Htete =   zero
       Htote =   zero
       Hzete = - DToocc(pp,ll,jj,ii) + DToocc(ll,pp,ii,jj)
       Hzote = - DToosc(pp,ll,jj,ii) + DToocs(ll,pp,ii,jj)
       
       Hteto =   zero
       Htoto =   zero
       Hzeto = - DToocs(pp,ll,jj,ii) + DToosc(ll,pp,ii,jj)
       Hzoto = - DTooss(pp,ll,jj,ii) + DTooss(ll,pp,ii,jj)  
       
       Hteze = + DToocc(pp,ll,jj,ii) - DToocc(ll,pp,ii,jj)
       Htoze = + DToosc(pp,ll,jj,ii) - DToocs(ll,pp,ii,jj)
       Hzeze =   zero 
       Hzoze =   zero
       
       Htezo = + DToocs(pp,ll,jj,ii) - DToosc(ll,pp,ii,jj)
       Htozo = + DTooss(pp,ll,jj,ii) - DTooss(ll,pp,ii,jj)
       Hzezo =   zero
       Hzozo =   zero
     
       id = Ate(lvol,0,ii)%i(ll) ; jd = Ate(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wtete ; dMD(id,jd) = Htete
       ;                         ; jd = Ato(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wtote ; dMD(id,jd) = Htote
       ;                         ; jd = Aze(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wzete ; dMD(id,jd) = Hzete
       ;                         ; jd = Azo(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wzote ; dMD(id,jd) = Hzote
       id = Ato(lvol,0,ii)%i(ll) ; jd = Ate(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wteto ; dMD(id,jd) = Hteto
       ;                         ; jd = Ato(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wtoto ; dMD(id,jd) = Htoto
       ;                         ; jd = Aze(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wzeto ; dMD(id,jd) = Hzeto
       ;                         ; jd = Azo(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wzoto ; dMD(id,jd) = Hzoto
       id = Aze(lvol,0,ii)%i(ll) ; jd = Ate(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wteze ; dMD(id,jd) = Hteze
       ;                         ; jd = Ato(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wtoze ; dMD(id,jd) = Htoze
       ;                         ; jd = Aze(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wzeze ; dMD(id,jd) = Hzeze
       ;                         ; jd = Azo(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wzoze ; dMD(id,jd) = Hzoze
       id = Azo(lvol,0,ii)%i(ll) ; jd = Ate(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wtezo ; dMD(id,jd) = Htezo
       ;                         ; jd = Ato(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wtozo ; dMD(id,jd) = Htozo
       ;                         ; jd = Aze(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wzezo ; dMD(id,jd) = Hzezo
       ;                         ; jd = Azo(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wzozo ; dMD(id,jd) = Hzozo
       
      enddo ! end of do pp ;
      
     enddo ! end of do jj ;
     
    enddo ! end of do ll ;
    
    if( Lcoordinatesingularity .and. ii.eq.1 ) then ; kk = 1
    else                                            ; kk = 0
    endif
    
    do ll = 0, lrad ! Chebyshev polynomial ;

     ;                  ; id = Ate(lvol,0,ii)%i(ll) ; jd = Lma(lvol,ii)         ; dMA(id,jd) = +      TT(ll, 0,0)
     ;                  ; id = Aze(lvol,0,ii)%i(ll) ; jd = Lmb(lvol,ii)         ; dMA(id,jd) = +      TT(ll,kk,0)
     if( ii.gt.1 ) then ; id = Ato(lvol,0,ii)%i(ll) ; jd = Lmc(lvol,ii)         ; dMA(id,jd) = +      TT(ll, 0,0)
      ;                 ; id = Azo(lvol,0,ii)%i(ll) ; jd = Lmd(lvol,ii)         ; dMA(id,jd) = +      TT(ll, 0,0)
      ;                 ; id = Ate(lvol,0,ii)%i(ll) ; jd = Lme(lvol,ii)         ; dMA(id,jd) = - ni * TT(ll, 1,0)
      ;                 ; id = Aze(lvol,0,ii)%i(ll) ; jd = Lme(lvol,ii)         ; dMA(id,jd) = - mi * TT(ll, 1,0)
      ;                 ; id = Ato(lvol,0,ii)%i(ll) ; jd = Lmf(lvol,ii)         ; dMA(id,jd) = + ni * TT(ll, 1,0)
      ;                 ; id = Azo(lvol,0,ii)%i(ll) ; jd = Lmf(lvol,ii)         ; dMA(id,jd) = + mi * TT(ll, 1,0)
     else               ; id = Ate(lvol,0,ii)%i(ll) ; jd = Lmg(lvol,ii)         ; dMA(id,jd) = +      TT(ll, 1,0)
      ;                 ; id = Aze(lvol,0,ii)%i(ll) ; jd = Lmh(lvol,ii)         ; dMA(id,jd) = +      TT(ll, 1,0)
     endif
     
    enddo ! end of do ll;
    
    do pp = 0, lrad
     ;                  ; id = Lma(lvol,ii)         ; jd = Ate(lvol,0,ii)%i(pp) ; dMA(id,jd) = +      TT(pp, 0,0)
     ;                  ; id = Lmb(lvol,ii)         ; jd = Aze(lvol,0,ii)%i(pp) ; dMA(id,jd) = +      TT(pp,kk,0)
     if( ii.gt.1 ) then ; id = Lmc(lvol,ii)         ; jd = Ato(lvol,0,ii)%i(pp) ; dMA(id,jd) = +      TT(pp, 0,0)
      ;                 ; id = Lmd(lvol,ii)         ; jd = Azo(lvol,0,ii)%i(pp) ; dMA(id,jd) = +      TT(pp, 0,0)
      ;                 ; id = Lme(lvol,ii)         ; jd = Ate(lvol,0,ii)%i(pp) ; dMA(id,jd) = - ni * TT(pp, 1,0)
      ;                 ; id = Lme(lvol,ii)         ; jd = Aze(lvol,0,ii)%i(pp) ; dMA(id,jd) = - mi * TT(pp, 1,0)
      ;                 ; id = Lmf(lvol,ii)         ; jd = Ato(lvol,0,ii)%i(pp) ; dMA(id,jd) = + ni * TT(pp, 1,0)
      ;                 ; id = Lmf(lvol,ii)         ; jd = Azo(lvol,0,ii)%i(pp) ; dMA(id,jd) = + mi * TT(pp, 1,0)
     else               ; id = Lmg(lvol,ii)         ; jd = Ate(lvol,0,ii)%i(pp) ; dMA(id,jd) = +      TT(pp, 1,0)
      ;                 ; id = Lmh(lvol,ii)         ; jd = Aze(lvol,0,ii)%i(pp) ; dMA(id,jd) = +      TT(pp, 1,0)
     endif
    enddo ! end of do pp ;
    
    ;if( ii.gt.1 ) then ; id = Lme(lvol,ii)         ;                           ; dMG(id   ) = - ( iVns(ii) + iBns(ii) )
    ;                   ; id = Lmf(lvol,ii)         ;                           ; dMG(id   ) = - ( iVnc(ii) + iBnc(ii) )
    ;else               ; id = Lmg(lvol,ii)         ;                           ; dMB(id, 1) = -       one
!   ;                   ; id = Lmh(lvol,ii)         ;                           ; dMB(id, 2) = -       one ! to be deleted;
    ;                   ; id = Lmh(lvol,ii)         ;                           ; dMB(id, 2) = +       one ! changed sign;
    ;endif
    
   enddo ! end of do ii ;
   
  endif ! end of if( YESstellsym ) ;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  
  if( Wmatrix ) then ! check symmetry of matrices ;
   
   do ii = 1, NN
    
    do jj = 1, NN
     if( abs(dMA(ii,jj)-dMA(jj,ii)) .gt. small*abs(dMA(ii,jj)+dMA(jj,ii)) ) then
      write(ounit,1000) myid, dMA(ii,jj), dMA(jj,ii), dMA(ii,jj)-dMA(jj,ii)
     endif
     
     if( abs(dMD(ii,jj)-dMD(jj,ii)) .gt. small*abs(dMD(ii,jj)+dMD(jj,ii)) ) then
      write(ounit,1001) myid, dMD(ii,jj), dMD(jj,ii), dMD(ii,jj)-dMD(jj,ii)
     endif
     
    enddo ! end of do jj;
    
   enddo ! end of do ii;
   
  endif ! end of if( Wmatrix ) ;
  
1000 format("matrix : " 10x " : myid="i3" : dMA(ii,jj)="es23.15", dMA(jj,ii)="es23.15", dMA(ii,jj)-dMA(jj,ii)="es13.5" ;")
1001 format("matrix : " 10x " : myid="i3" : dMD(ii,jj)="es23.15", dMD(jj,ii)="es23.15", dMD(ii,jj)-dMD(jj,ii)="es13.5" ;")
  
#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(matrix)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine matrix

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
