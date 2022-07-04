!> \defgroup grp_build_matrices Build matrices
!>
!> \latexonly
!> \definecolor{Orange}{rgb}{1.0,0.5,0.0}
!> \definecolor{Cerulean}{rgb}{0.0,0.5,1.0}
!> \endlatexonly
!>
!> \file
!> \brief Constructs energy and helicity matrices that represent the Beltrami linear system.

!> \brief Constructs energy and helicity matrices that represent the Beltrami linear system.
!> \ingroup grp_build_matrices
!> **gauge conditions**
!>
!> <ul>
!> <li> In the \f$v\f$-th annulus, bounded by the \f$(v-1)\f$-th and \f$v\f$-th interfaces,
!>       a general covariant representation of the magnetic vector-potential is written
!>       \f{eqnarray}{ {\bf \bar A}=\bar A_{s} \nabla s + \bar A_{\theta} \nabla \theta + \bar A_{\zeta} \nabla \zeta. \f} </li>
!> <li> To this add \f$\nabla g(s,\theta,\zeta)\f$, where \f$g\f$ satisfies
!>       \f{eqnarray}{ \begin{array}{cccccccccccccccccccccccccccccccccccccc}
!>           \partial_s      g( s,\theta,\zeta) & = & - & \bar A_{s     }( s,\theta,\zeta) \\
!>           \partial_\theta g(-1,\theta,\zeta) & = & - & \bar A_{\theta}(-1,\theta,\zeta) \\
!>           \partial_\zeta  g(-1,     0,\zeta) & = & - & \bar A_{\zeta }(-1,     0,\zeta).
!>       \end{array}
!>      \f} </li>
!> <li> Then \f${\bf A}={\bf \bar A}+\nabla g\f$ is given by \f${\bf A}=A_{\theta}\nabla\theta+A_{\zeta}\nabla\zeta\f$ with
!>       \f{eqnarray}{
!>        A_{\theta}(-1,\theta,\zeta) &=& 0 \label{eq:At_matrixgauge_matrix} \\
!>        A_{\zeta }(-1,     0,\zeta) &=& 0 \label{eq:Az_matrixgauge_matrix}
!>       \f} </li>
!> <li> This specifies the gauge: to see this, notice that no gauge term can be added without violating the conditions in Eqn.\f$(\ref{eq:At_matrixgauge_matrix})\f$ or Eqn.\f$(\ref{eq:Az_matrixgauge_matrix})\f$. </li>
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
!>       With the above condition on \f$A_\theta\f$ given in Eqn.\f$(\ref{eq:At_matrixgauge_matrix})\f$, this gives \f$\partial_\theta A_\zeta=0\f$, which with Eqn.\f$(\ref{eq:Az_matrixgauge_matrix})\f$ gives
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
!> <li> The components of the vector potential, \f${\bf A}=A_\theta \nabla + A_\zeta \nabla \zeta \f$, are
!>      \f{eqnarray}{
!>        A_\theta(s,\theta,\zeta) &=& \sum_{i,l} {\color{red} A_{\theta,e,i,l}} \; {\overline T}_{l,i}(s) \cos\alpha_i + \sum_{i,l} {\color{Orange}  A_{\theta,o,i,l}} \; {\overline T}_{l,i}(s) \sin\alpha_i, \label{eq:At_matrix} \\
!>        A_\zeta( s,\theta,\zeta) &=& \sum_{i,l} {\color{blue}A_{\zeta, e,i,l}} \; {\overline T}_{l,i}(s) \cos\alpha_i + \sum_{i,l} {\color{Cerulean}A_{\zeta ,o,i,l}} \; {\overline T}_{l,i}(s) \sin\alpha_i, \label{eq:Az_matrix}
!>      \f}
!>      where \f${\overline T}_{l,i}(s) \f$ is the __recombined__ Chebyshev polynomial in a volume without an axis, or __modified__ Zernike polynomial in a volume with an axis
!>      (i.e. only in the innermost volume, and only in cylindrical and toroidal geometry.)
!>      , and \f$\alpha_j \equiv m_j\theta-n_j\zeta\f$.</li>
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
!> <li> __Note:__ in SPEC version >3.00, the basis recombination method is used to ensure the boundary condition on the inner side of an interface.
!>      The lagrange multipliers \f$a_i, b_i, c_i, d_i\f$ are no longer used in volumes without a coordinate singularity.
!>      In a volume with a coordinate singularity, they are used only\f$a_i, c_i\f$ with $m=0,1$ are excluded also due to Zernike basis recombination.
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
!> <li> The required geometric information is calculated in ma00aa(). </li>
!> </ul>
!>
!> @param[in] lvol
!> @param[in] mn
!> @param[in] lrad
subroutine matrix( lvol, mn, lrad )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  use constants, only : zero, one, two

  use numerical, only : small

  use fileunits, only : ounit

  use inputlist, only : Wmacros, Wmatrix, mpol

  use cputiming, only : Tmatrix

  use allglobal, only : ncpu, myid, cpus, MPI_COMM_SPEC, &
                        YESstellsym, NOTstellsym, &
                        im, in, &
                        NAdof, &
                        dMA, dMD, dMB, dMG, &
                        Ate, Ato, Aze, Azo, &
                        iVns, iBns, iVnc, iBnc, &
                        Lma, Lmb, Lmc, Lmd, Lme, Lmf, Lmg, Lmh, &
                        Lcoordinatesingularity, TT, RTT, RTM, &
                        DToocc, DToocs, DToosc, DTooss, &
                        TTsscc, TTsscs, TTsssc, TTssss, &
                        TDstcc, TDstcs, TDstsc, TDstss, &
                        TDszcc, TDszcs, TDszsc, TDszss, &
                        DDttcc, DDttcs, DDttsc, DDttss, &
                        DDtzcc, DDtzcs, DDtzsc, DDtzss, &
                        DDzzcc, DDzzcs, DDzzsc, DDzzss, &
                        dBdX

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER, intent(in)  :: lvol, mn, lrad

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER              :: NN, ii, jj, ll, kk, pp, ll1, pp1, mi, ni, mj, nj, mimj, minj, nimj, ninj, mjmi, mjni, njmi, njni, id, jd

  REAL                 :: Wtete, Wteto, Wtote, Wtoto
  REAL                 :: Wteze, Wtezo, Wtoze, Wtozo
  REAL                 :: Wzete, Wzeto, Wzote, Wzoto
  REAL                 :: Wzeze, Wzezo, Wzoze, Wzozo

  REAL                 :: Htete, Hteto, Htote, Htoto
  REAL                 :: Hteze, Htezo, Htoze, Htozo
  REAL                 :: Hzete, Hzeto, Hzote, Hzoto
  REAL                 :: Hzeze, Hzezo, Hzoze, Hzozo

  REAL,allocatable     :: TTdata(:,:,:), TTMdata(:,:) ! queues to construct sparse matrices

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

  SALLOCATE( TTdata, (0:lrad, 0:mpol, 0:1), zero)
  SALLOCATE( TTMdata, (0:lrad, 0:mpol), zero)

  ! fill in Zernike/Chebyshev polynomials depending on Lcooridnatesingularity
  if (Lcoordinatesingularity) then
    TTdata(0:lrad,0:mpol,0:1) = RTT(0:lrad,0:mpol,0:1,0)
    TTMdata(0:lrad,0:mpol) = RTM(0:lrad,0:mpol)
  else
    do ii = 0, mpol
      TTdata(0:lrad,ii,0:1) = TT(0:lrad,0:1,0)
      TTMdata(0:lrad,ii) = TT(0:lrad,0,0)
    enddo
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( YESstellsym ) then
!$OMP PARALLEL DO PRIVATE(ii,jj,ll,pp,mi,ni,mj,nj,mimj,minj,nimj,ninj,ll1,pp1,Wtete,Wzete,Wteze,Wzeze,Htete,Hzete,Hteze,Hzeze,id,jd,kk) SHARED(lvol,lrad)
   do ii = 1, mn ; mi = im(ii) ; ni = in(ii)

    do jj = 1, mn ; mj = im(jj) ; nj = in(jj) ; mimj = mi * mj ; minj = mi * nj ; nimj = ni * mj ; ninj = ni * nj

     do ll = 0, lrad

      do pp = 0, lrad

       if (Lcoordinatesingularity) then
        if (ll < mi .or. pp < mj) cycle ! rule out zero components of Zernike; 02 Jul 19
        if (mod(ll+mi,2)+mod(pp+mj,2)>0) cycle ! rule out zero components of Zernike; 02 Jul 19
        ll1 = (ll - mod(ll, 2)) / 2 ! shrinked dof for Zernike; 02 Jul 19
        pp1 = (pp - mod(pp, 2)) / 2 ! shrinked dof for Zernike; 02 Jul 19
       else
        ll1 = ll
        pp1 = pp
       end if

       Wtete = + 2 * ninj * TTssss(ll1,pp1,ii,jj) - 2 * ni      * TDszsc(ll1,pp1,ii,jj) - 2      * nj * TDszsc(pp1,ll1,jj,ii) + 2 * DDzzcc(ll1,pp1,ii,jj)
       Wzete = + 2 * nimj * TTssss(ll1,pp1,ii,jj) + 2 * ni      * TDstsc(ll1,pp1,ii,jj) - 2      * mj * TDszsc(pp1,ll1,jj,ii) - 2 * DDtzcc(pp1,ll1,jj,ii)
       Wteze = + 2 * minj * TTssss(ll1,pp1,ii,jj) + 2      * nj * TDstsc(pp1,ll1,jj,ii) - 2 * mi      * TDszsc(ll1,pp1,ii,jj) - 2 * DDtzcc(ll1,pp1,ii,jj)
       Wzeze = + 2 * mimj * TTssss(ll1,pp1,ii,jj) + 2 * mi      * TDstsc(ll1,pp1,ii,jj) + 2      * mj * TDstsc(pp1,ll1,jj,ii) + 2 * DDttcc(ll1,pp1,ii,jj)

       Htete =   zero
       Hzete = - DToocc(pp1,ll1,jj,ii) + DToocc(ll1,pp1,ii,jj)
       Hteze = + DToocc(pp1,ll1,jj,ii) - DToocc(ll1,pp1,ii,jj)
       Hzeze =   zero

       id = Ate(lvol,0,ii)%i(ll) ; jd = Ate(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wtete ; dMD(id,jd) = Htete
       ;                         ; jd = Aze(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wzete ; dMD(id,jd) = Hzete
       id = Aze(lvol,0,ii)%i(ll) ; jd = Ate(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wteze ; dMD(id,jd) = Hteze
       ;                         ; jd = Aze(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wzeze ; dMD(id,jd) = Hzeze

      enddo ! end of do pp ;

     enddo ! end of do ll ;

    enddo ! end of do jj ;

    if (dBdX%L) cycle

    if( Lcoordinatesingularity .and. ii.eq.1 ) then ; kk = 1
    else                                            ; kk = 0
    endif

    do ll = 0, lrad
    ;                  ; id = Ate(lvol,0,ii)%i(ll) ; jd = Lma(lvol,  ii)       ; dMA(id,jd) = +      TTMdata(ll, mi)
    ;                  ; id = Aze(lvol,0,ii)%i(ll) ; jd = Lmb(lvol,  ii)       ; dMA(id,jd) = +      TTdata(ll, mi,kk)
    if( ii.gt.1 ) then ; id = Ate(lvol,0,ii)%i(ll) ; jd = Lme(lvol,  ii)       ; dMA(id,jd) = - ni * TTdata(ll, mi, 1)
      ;                 ; id = Aze(lvol,0,ii)%i(ll) ; jd = Lme(lvol,  ii)       ; dMA(id,jd) = - mi * TTdata(ll, mi, 1)
    else               ; id = Ate(lvol,0,ii)%i(ll) ; jd = Lmg(lvol,  ii)       ; dMA(id,jd) = +      TTdata(ll, mi, 1)
      ;                 ; id = Aze(lvol,0,ii)%i(ll) ; jd = Lmh(lvol,  ii)       ; dMA(id,jd) = +      TTdata(ll, mi, 1)
    endif

    enddo ! end of do ll ;

    do pp = 0, lrad     ; id = Lma(lvol,  ii)       ; jd = Ate(lvol,0,ii)%i(pp) ; dMA(id,jd) = +      TTMdata(pp, mi)
    ;                  ; id = Lmb(lvol,  ii)       ; jd = Aze(lvol,0,ii)%i(pp) ; dMA(id,jd) = +      TTdata(pp, mi,kk)
    if( ii.gt.1 ) then ; id = Lme(lvol,  ii)       ; jd = Ate(lvol,0,ii)%i(pp) ; dMA(id,jd) = - ni * TTdata(pp, mi, 1)
      ;                 ; id = Lme(lvol,  ii)       ; jd = Aze(lvol,0,ii)%i(pp) ; dMA(id,jd) = - mi * TTdata(pp, mi, 1)
    else               ; id = Lmg(lvol,  ii)       ; jd = Ate(lvol,0,ii)%i(pp) ; dMA(id,jd) = +      TTdata(pp, mi, 1)
      ;                 ; id = Lmh(lvol,  ii)       ; jd = Aze(lvol,0,ii)%i(pp) ; dMA(id,jd) = +      TTdata(pp, mi, 1)
    endif
    enddo ! end of do pp ;

   enddo ! end of do ii ;
!$OMP END PARALLEL DO

  else ! NOTstellsym ;

!$OMP PARALLEL DO PRIVATE(ii,jj,ll,pp,mi,ni,mj,nj,mjmi,mjni,njmi,njni,ll1,pp1,Wtete,Wzete,Wteze,Wzeze,Htete,Hzete,Hteze,Hzeze,Wteto,Wzeto,Wtezo,Wzezo,Hteto,Hzeto,Htezo,Hzezo,Wtote,Wzote,Wtoze,Wzoze,Htote,Hzote,Htoze,Hzoze,Wtoto,Wzoto,Wtozo,Wzozo,Htoto,Hzoto,Htozo,Hzozo,id,jd,kk) SHARED(lvol,lrad)
   do ii = 1, mn ; mi = im(ii) ; ni = in(ii)

    do jj = 1, mn ; mj = im(jj) ; nj = in(jj) ; mjmi = mi * mj ; mjni = ni * mj ; njmi = mi * nj ; njni = ni * nj

     do ll = 0, lrad

      do pp = 0, lrad

       if (Lcoordinatesingularity) then
        if (ll < mi .or. pp < mj) cycle ! rule out zero components of Zernike; 02 Jul 19
        if (mod(ll+mi,2)+mod(pp+mj,2)>0) cycle ! rule out zero components of Zernike; 02 Jul 19
        ll1 = (ll - mod(ll, 2)) / 2! shrinked dof for Zernike; 02 Jul 19
        pp1 = (pp - mod(pp, 2)) / 2 ! shrinked dof for Zernike; 02 Jul 19
       else
        ll1 = ll
        pp1 = pp
       end if

       Wtete = 2 * ( + njni * TTssss(pp1,ll1,jj,ii) - nj * TDszsc(pp1,ll1,jj,ii) - ni * TDszsc(ll1,pp1,ii,jj) + DDzzcc(pp1,ll1,jj,ii) )
       Wtote = 2 * ( - njni * TTsscs(pp1,ll1,jj,ii) + nj * TDszcc(pp1,ll1,jj,ii) - ni * TDszss(ll1,pp1,ii,jj) + DDzzsc(pp1,ll1,jj,ii) )
       Wzete = 2 * ( + mjni * TTssss(pp1,ll1,jj,ii) - mj * TDszsc(pp1,ll1,jj,ii) + ni * TDstsc(ll1,pp1,ii,jj) - DDtzcc(pp1,ll1,jj,ii) )
       Wzote = 2 * ( - mjni * TTsscs(pp1,ll1,jj,ii) + mj * TDszcc(pp1,ll1,jj,ii) + ni * TDstss(ll1,pp1,ii,jj) - DDtzsc(pp1,ll1,jj,ii) )

       Wteto = 2 * ( - njni * TTsssc(pp1,ll1,jj,ii) - nj * TDszss(pp1,ll1,jj,ii) + ni * TDszcc(ll1,pp1,ii,jj) + DDzzcs(pp1,ll1,jj,ii) )
       Wtoto = 2 * ( + njni * TTsscc(pp1,ll1,jj,ii) + nj * TDszcs(pp1,ll1,jj,ii) + ni * TDszcs(ll1,pp1,ii,jj) + DDzzss(pp1,ll1,jj,ii) )
       Wzeto = 2 * ( - mjni * TTsssc(pp1,ll1,jj,ii) - mj * TDszss(pp1,ll1,jj,ii) - ni * TDstcc(ll1,pp1,ii,jj) - DDtzcs(pp1,ll1,jj,ii) )
       Wzoto = 2 * ( + mjni * TTsscc(pp1,ll1,jj,ii) + mj * TDszcs(pp1,ll1,jj,ii) - ni * TDstcs(ll1,pp1,ii,jj) - DDtzss(pp1,ll1,jj,ii) )

       Wteze = 2 * ( + njmi * TTssss(pp1,ll1,jj,ii) + nj * TDstsc(pp1,ll1,jj,ii) - mi * TDszsc(ll1,pp1,ii,jj) - DDtzcc(pp1,ll1,jj,ii) )
       Wtoze = 2 * ( - njmi * TTsscs(pp1,ll1,jj,ii) - nj * TDstcc(pp1,ll1,jj,ii) - mi * TDszss(ll1,pp1,ii,jj) - DDtzsc(pp1,ll1,jj,ii) )
       Wzeze = 2 * ( + mjmi * TTssss(pp1,ll1,jj,ii) + mj * TDstsc(pp1,ll1,jj,ii) + mi * TDstsc(ll1,pp1,ii,jj) + DDttcc(pp1,ll1,jj,ii) )
       Wzoze = 2 * ( - mjmi * TTsscs(pp1,ll1,jj,ii) - mj * TDstcc(pp1,ll1,jj,ii) + mi * TDstss(ll1,pp1,ii,jj) + DDttsc(pp1,ll1,jj,ii) )

       Wtezo = 2 * ( - njmi * TTsssc(pp1,ll1,jj,ii) + nj * TDstss(pp1,ll1,jj,ii) + mi * TDszcc(ll1,pp1,ii,jj) - DDtzcs(pp1,ll1,jj,ii) )
       Wtozo = 2 * ( + njmi * TTsscc(pp1,ll1,jj,ii) - nj * TDstcs(pp1,ll1,jj,ii) + mi * TDszcs(ll1,pp1,ii,jj) - DDtzss(pp1,ll1,jj,ii) )
       Wzezo = 2 * ( - mjmi * TTsssc(pp1,ll1,jj,ii) + mj * TDstss(pp1,ll1,jj,ii) - mi * TDstcc(ll1,pp1,ii,jj) + DDttcs(pp1,ll1,jj,ii) )
       Wzozo = 2 * ( + mjmi * TTsscc(pp1,ll1,jj,ii) - mj * TDstcs(pp1,ll1,jj,ii) - mi * TDstcs(ll1,pp1,ii,jj) + DDttss(pp1,ll1,jj,ii) )

       Htete =   zero
       Htote =   zero
       Hzete = - DToocc(pp1,ll1,jj,ii) + DToocc(ll1,pp1,ii,jj)
       Hzote = - DToosc(pp1,ll1,jj,ii) + DToocs(ll1,pp1,ii,jj)

       Hteto =   zero
       Htoto =   zero
       Hzeto = - DToocs(pp1,ll1,jj,ii) + DToosc(ll1,pp1,ii,jj)
       Hzoto = - DTooss(pp1,ll1,jj,ii) + DTooss(ll1,pp1,ii,jj)

       Hteze = + DToocc(pp1,ll1,jj,ii) - DToocc(ll1,pp1,ii,jj)
       Htoze = + DToosc(pp1,ll1,jj,ii) - DToocs(ll1,pp1,ii,jj)
       Hzeze =   zero
       Hzoze =   zero

       Htezo = + DToocs(pp1,ll1,jj,ii) - DToosc(ll1,pp1,ii,jj)
       Htozo = + DTooss(pp1,ll1,jj,ii) - DTooss(ll1,pp1,ii,jj)
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

    if (dBdX%L) cycle

    if( Lcoordinatesingularity .and. ii.eq.1 ) then ; kk = 1
    else                                            ; kk = 0
    endif

    do ll = 0, lrad ! Chebyshev polynomial ;

    ;                  ; id = Ate(lvol,0,ii)%i(ll) ; jd = Lma(lvol,ii)         ; dMA(id,jd) = TTMdata(ll, mi)
    ;                  ; id = Aze(lvol,0,ii)%i(ll) ; jd = Lmb(lvol,ii)         ; dMA(id,jd) = TTdata(ll, mi,kk)
    if( ii.gt.1 ) then ; id = Ato(lvol,0,ii)%i(ll) ; jd = Lmc(lvol,ii)         ; dMA(id,jd) = TTMdata(ll, mi)
      ;                 ; id = Azo(lvol,0,ii)%i(ll) ; jd = Lmd(lvol,ii)         ; dMA(id,jd) = TTdata(ll, mi,0)
      ;                 ; id = Ate(lvol,0,ii)%i(ll) ; jd = Lme(lvol,ii)         ; dMA(id,jd) = - ni * TTdata(ll, mi, 1)
      ;                 ; id = Aze(lvol,0,ii)%i(ll) ; jd = Lme(lvol,ii)         ; dMA(id,jd) = - mi * TTdata(ll, mi, 1)
      ;                 ; id = Ato(lvol,0,ii)%i(ll) ; jd = Lmf(lvol,ii)         ; dMA(id,jd) = + ni * TTdata(ll, mi, 1)
      ;                 ; id = Azo(lvol,0,ii)%i(ll) ; jd = Lmf(lvol,ii)         ; dMA(id,jd) = + mi * TTdata(ll, mi, 1)
    else               ; id = Ate(lvol,0,ii)%i(ll) ; jd = Lmg(lvol,ii)         ; dMA(id,jd) = +      TTdata(ll, mi, 1)
      ;                 ; id = Aze(lvol,0,ii)%i(ll) ; jd = Lmh(lvol,ii)         ; dMA(id,jd) = +      TTdata(ll, mi, 1)
    endif

    enddo ! end of do ll;

    do pp = 0, lrad
    ;                  ; id = Lma(lvol,ii)         ; jd = Ate(lvol,0,ii)%i(pp) ; dMA(id,jd) = TTMdata(pp, mi)
    ;                  ; id = Lmb(lvol,ii)         ; jd = Aze(lvol,0,ii)%i(pp) ; dMA(id,jd) = TTdata(pp, mi,kk)
    if( ii.gt.1 ) then ; id = Lmc(lvol,ii)         ; jd = Ato(lvol,0,ii)%i(pp) ; dMA(id,jd) = TTMdata(pp, mi)
      ;                 ; id = Lmd(lvol,ii)         ; jd = Azo(lvol,0,ii)%i(pp) ; dMA(id,jd) = TTdata(pp, mi,0)
      ;                 ; id = Lme(lvol,ii)         ; jd = Ate(lvol,0,ii)%i(pp) ; dMA(id,jd) = - ni * TTdata(pp, mi, 1)
      ;                 ; id = Lme(lvol,ii)         ; jd = Aze(lvol,0,ii)%i(pp) ; dMA(id,jd) = - mi * TTdata(pp, mi, 1)
      ;                 ; id = Lmf(lvol,ii)         ; jd = Ato(lvol,0,ii)%i(pp) ; dMA(id,jd) = + ni * TTdata(pp, mi, 1)
      ;                 ; id = Lmf(lvol,ii)         ; jd = Azo(lvol,0,ii)%i(pp) ; dMA(id,jd) = + mi * TTdata(pp, mi, 1)
    else               ; id = Lmg(lvol,ii)         ; jd = Ate(lvol,0,ii)%i(pp) ; dMA(id,jd) = +      TTdata(pp, mi, 1)
      ;                 ; id = Lmh(lvol,ii)         ; jd = Aze(lvol,0,ii)%i(pp) ; dMA(id,jd) = +      TTdata(pp, mi, 1)
    endif
    enddo ! end of do pp ;

   enddo ! end of do ii ;
!$OMP END PARALLEL DO
  endif ! end of if( YESstellsym ) ;

  ! call subroutine matrixBG to construct dMB and dMG
  WCALL( matrix, matrixBG, ( lvol, mn, lrad ) )

  DALLOCATE( TTdata )
  DALLOCATE( TTMdata )

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

subroutine matrixBG( lvol, mn, lrad )
  ! only compute the dMB and dMG matrix for matrix-free mode
  use constants, only : zero, one
  use allglobal, only : NAdof, im, in,&
                        dMG, dMB, YESstellsym, &
                        iVnc, iVns, iBnc, iBns, &
                        Lme, Lmf, Lmg, Lmh
  implicit none
  INTEGER, intent(in)  :: lvol, mn, lrad

  INTEGER :: NN, ii, id, mi, ni

  NN = NAdof(lvol) ! shorthand;

  dMB(0:NN,1: 2) = zero
  dMG(0:NN     ) = zero

  if( YESstellsym ) then

   do ii = 1, mn ; mi = im(ii) ; ni = in(ii)

    ;if( ii.gt.1 ) then ; id = Lme(lvol,  ii)       ;                           ; dMG(id   ) = - ( iVns(ii) + iBns(ii) )
    ;else               ; id = Lmg(lvol,  ii)       ;                           ; dMB(id, 1) = -       one
!   ;                   ; id = Lmh(lvol,  ii)       ;                           ; dMB(id, 2) = -       one ! to be deleted;
    ;                   ; id = Lmh(lvol,  ii)       ;                           ; dMB(id, 2) = +       one ! changed sign;
    ;endif

   enddo ! end of do ii ;

  else ! NOTstellsym ;

   do ii = 1, mn ; mi = im(ii) ; ni = in(ii)

    ;if( ii.gt.1 ) then ; id = Lme(lvol,ii)         ;                           ; dMG(id   ) = - ( iVns(ii) + iBns(ii) )
    ;                   ; id = Lmf(lvol,ii)         ;                           ; dMG(id   ) = - ( iVnc(ii) + iBnc(ii) )
    ;else               ; id = Lmg(lvol,ii)         ;                           ; dMB(id, 1) = -       one
!   ;                   ; id = Lmh(lvol,ii)         ;                           ; dMB(id, 2) = -       one ! to be deleted;
    ;                   ; id = Lmh(lvol,ii)         ;                           ; dMB(id, 2) = +       one ! changed sign;
    ;endif

   enddo ! end of do ii ;

  endif ! end of if( YESstellsym ) ;

  return

end subroutine matrixBG

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
