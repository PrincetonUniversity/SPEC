!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (build matrices) ! Constructs matrices that represent the Beltrami linear system.

!latex \briefly{Constructs energy and helicity matrices that represent the Beltrami linear system.}

!latex \calledby{\link{dforce}}
!latex \calls{\link{}}

!latex \tableofcontents

!latex \subsection{gauge conditions} 

!latex \begin{enumerate}

!latex \item In the \mbox{$v$}-th annulus, bounded by the \mbox{$(v-1)$}-th and \mbox{$v$}-th interfaces,
!latex       a general covariant representation of the magnetic vector-potential is written
!latex       \be {\bf \bar A}=\bar A_{s} \nabla s + \bar A_{\t} \nabla \t + \bar A_{\z} \nabla \z. \ee
!latex \item To this add \mbox{$\nabla g(s,\t,\z)$}, where $g$ satisfies
!latex       \be \begin{array}{cccccccccccccccccccccccccccccccccccccc}
!latex           \partial_\s g(\s,\t,\z) & = & - & \bar A_{\s}(\s,\t,\z) \\
!latex           \partial_\t g(-1,\t,\z) & = & - & \bar A_{\t}(-1,\t,\z) \\
!latex           \partial_\z g(-1, 0,\z) & = & - & \bar A_{\z}(-1, 0,\z).
!latex       \end{array}\ee
!latex \item Then \mbox{${\bf A}={\bf \bar A}+\nabla g$} is given by \mbox{${\bf A}=A_{\t}\nabla\t+A_{\z}\nabla\z$} with
!latex       \be \
!latex        A_{\t}(-1,\t,\z) &=& 0 \label{eq:Atgauge} \\
!latex        A_{\z}(-1, 0,\z) &=& 0 \label{eq:Azgauge}
!latex       \ee
!latex \item This specifies the gauge: to see this, notice that no gauge term can be added without violating the conditions in \Eqn{Atgauge} or \Eqn{Azgauge}.
!latex \item Note that the gauge employed in each volume is distinct.

!latex \end{enumerate} 

!latex \subsection{boundary conditions} 

!latex \begin{enumerate}

!latex \item The magnetic field is 
!latex       $\sqrt g \, {\bf B} = (\partial_\t A_\z - \partial_\z A_\t)\;{\bf e}_\s - \partial_\s A_\z \;{\bf e}_\t + \partial_\s A_\t \;{\bf e}_\z$.
!latex \item In the annular volumes, the condition that the field is tangential to the inner interface, $\sqrt g {\bf B}\cdot\nabla s=0$ at $s=-1$,
!latex       gives $\partial_\t A_\z - \partial_\z A_\t = 0$.
!latex       With the above condition on $A_\t$ given in \Eqn{Atgauge}, this gives $\partial_\t A_\z=0$, which with \Eqn{Azgauge} gives
!latex       \be A_\z(-1,\t,\z)=0.
!latex       \ee
!latex \item The condition at the outer interface, $s=+1$, is that the field is $\sqrt g \, {\bf B}\cdot\nabla s = \partial_\t A_\z - \partial_\z A_\t = b$,
!latex       where $b$ is supplied by the user.
!latex       For each of the plasma regions, $b=0$.
!latex       For the vacuum region, generally $b\ne0$.

!latex \end{enumerate} 

!latex \subsection{enclosed fluxes} 

!latex \begin{enumerate}

!latex \item In the plasma regions, the enclosed fluxes must be constrained.
!latex \item The toroidal and poloidal fluxes enclosed in each volume are determined using
!latex       \be \int_S {\bf B}\cdot{\bf ds}=\int_{\partial S}{\bf A}\cdot {\bf dl}.
!latex       \ee

!latex \end{enumerate} 

!latex \subsection{Fourier-Chebyshev representation} 

!latex \begin{enumerate}

!latex \item \Ais
!latex \item \Bis
!latex \item The components of the velocity, ${\bf v} \equiv v_\s \nabla \s + v_\t \nabla \t + v_\z \nabla \z$, are
!latex       \be v_\s(\s,\t,\z) &=& \sum_{i,l} \vse{i,l} \; \bT{l,i}(s) \ci + \sum_{i,l} \vso{i,l} \; \bT{l,i}(s) \si, \\
!latex           v_\t(\s,\t,\z) &=& \sum_{i,l} \vte{i,l} \; \bT{l,i}(s) \ci + \sum_{i,l} \vto{i,l} \; \bT{l,i}(s) \si, \\
!latex           v_\z(\s,\t,\z) &=& \sum_{i,l} \vze{i,l} \; \bT{l,i}(s) \ci + \sum_{i,l} \vzo{i,l} \; \bT{l,i}(s) \si.
!latex       \ee

!latex \end{enumerate} 

!latex \subsection{constrained energy functional} 

!latex \begin{enumerate}

!latex \item The constrained energy functional in each volume depends on the vector potential and the Lagrange multipliers, 
!latex \be {\cal F} \equiv
!latex {\cal F}[\Ate{i,l},\Aze{i,l},\Ato{i,l},\Azo{i,l},\vse{i,l},\vso{i,l},\vte{i,l},\vto{i,l},\vze{i,l},\vzo{i,l},\mu,a_i,b_i,c_i,d_i,e_i,f_i,g_1,h_1],
!latex \ee
!latex       and is given by:
!latex \be \begin{array}{cclcclcccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc}
!latex     {\cal F}
!latex     & \equiv & \ds                        \int {\bf B} \cdot {\bf B} \, dv
!latex       +                                   \int {\bf v} \cdot {\bf v} \, dv
!latex        -       \ds       \mu       \left[ \int {\bf A} \cdot {\bf B} \, dv - K \right]                                                  \\
!latex     &  + &   & \ds \sum_{i=1} & a_i       & \ds \left[  \sum_l \Ate{i,l} T_l(-1) - 0             \right]                                \\
!latex     &  + &   & \ds \sum_{i=1} & b_i       & \ds \left[  \sum_l \Aze{i,l} T_l(-1) - 0             \right]                                \\
!latex     &  + &   & \ds \sum_{i=2} & c_i       & \ds \left[  \sum_l \Ato{i,l} T_l(-1) - 0             \right]                                \\
!latex     &  + &   & \ds \sum_{i=2} & d_i       & \ds \left[  \sum_l \Azo{i,l} T_l(-1) - 0             \right]                                \\
!latex     &  + &   & \ds \sum_{i=2} & e_i       & \ds \left[  \sum_l \left( - m_i \Aze{i,l} - n_i \Ate{i,l} \right) T_l(+1) - b_{s,i} \right] \\
!latex     &  + &   & \ds \sum_{i=2} & f_i       & \ds \left[  \sum_l \left( + m_i \Azo{i,l} + n_i \Ato{i,l} \right) T_l(+1) - b_{c,i} \right] \\
!latex     &  + &   & \ds            & g_1       & \ds \left[  \sum_l \Ate{1,l} T_l(+1) - \Delta \psi_t \right]                                \\
!latex     &  + &   & \ds            & h_1       & \ds \left[  \sum_l \Aze{1,l} T_l(+1) + \Delta \psi_p \right]                             
!latex     \end{array}
!latex \ee   
!latex       where
!latex       \bi
!latex       \item[  i.] $a_i$, $b_i$, $c_i$ and $d_i$ are Lagrange multipliers used to enforce the combined gauge and interface boundary condition
!latex                   on the inner interface,
!latex       \item[ ii.] $e_i$ and $f_i$               are Lagrange multipliers used to enforce the interface boundary condition on the outer interface, 
!latex                   namely \mbox{$\sqrt g \, {\bf B}\cdot\nabla s = b$}; and
!latex       \item[iii.] $g_1$    and $h_1$            are Lagrange multipliers used to enforce the constraints on the enclosed fluxes.

!latex       \ei
!latex \item In each plasma volume the boundary condition on the outer interface is $b=0$.
!latex \item In the vacuum volume (only for free-boundary), we may set $\mu=0$.

!latex \end{enumerate} 

!latex \subsection{derivatives of magnetic energy integrals}

!latex \begin{enumerate}

!latex \item The first  derivatives of $\int \! dv \; {\bf B} \! \cdot \! {\bf B}$ with respect to $\Ate{i,l}$, $\Ato{i,l}$, $\Aze{i,l}$ and $\Azo{i,l}$ are
!latex \be 
!latex \frac{\partial}{\partial \Ate{i,l}} \int \!\! dv \; {\bf B} \cdot {\bf B} & = & 
!latex 2 \int \!\! dv \; {\bf B} \cdot \frac{\partial {\bf B}}{\partial \Ate{i,l}} =
!latex      2 \int \!\! dv \; {\bf B} \cdot \left[ - n_i \bT{l,i}  \si \, {\bf e}_\s + \bT{l,i}' \ci \, {\bf e}_\z \right] / \sqrt g \\
!latex \frac{\partial}{\partial \Ato{i,l}} \int \!\! dv \; {\bf B} \cdot {\bf B} & = & 
!latex 2 \int \!\! dv \; {\bf B} \cdot \frac{\partial {\bf B}}{\partial \Ato{i,l}} =
!latex      2 \int \!\! dv \; {\bf B} \cdot \left[ + n_i \bT{l,i}  \ci \, {\bf e}_\s + \bT{l,i}' \si \, {\bf e}_\z \right] / \sqrt g \\
!latex \frac{\partial}{\partial \Aze{i,l}} \int \!\! dv \; {\bf B} \cdot {\bf B} & = & 
!latex 2 \int \!\! dv \; {\bf B} \cdot \frac{\partial {\bf B}}{\partial \Aze{i,l}} =
!latex      2 \int \!\! dv \; {\bf B} \cdot \left[ - m_i \bT{l,i}  \si \, {\bf e}_\s - \bT{l,i}' \ci \, {\bf e}_\t \right] / \sqrt g \\
!latex \frac{\partial}{\partial \Azo{i,l}} \int \!\! dv \; {\bf B} \cdot {\bf B} & = & 
!latex 2 \int \!\! dv \; {\bf B} \cdot \frac{\partial {\bf B}}{\partial \Azo{i,l}} =
!latex      2 \int \!\! dv \; {\bf B} \cdot \left[ + m_i \bT{l,i}  \ci \, {\bf e}_\s - \bT{l,i}' \si \, {\bf e}_\t \right] / \sqrt g 
!latex \ee
!latex \item The second derivatives of $\int \! dv \; {\bf B} \! \cdot \! {\bf B}$ with respect to $\Ate{i,l}$, $\Ato{i,l}$, $\Aze{i,l}$ and $\Azo{i,l}$ are
!latex \be 
!latex \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial \Ate{j,p}} \frac{\partial}{\partial \Ate{i,l}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!latex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!latex  ( + n_j n_i \bT{p,j}  \bT{l,i}  s_j s_i g_{\s\s} - n_j \bT{p,j}  \bT{l,i}' s_j c_i g_{\s\z}
!latex    -     n_i \bT{l,i}  \bT{p,j}' s_i c_j g_{\s\z} +     \bT{p,j}' \bT{l,i}' c_j c_i g_{\z\z}) / \sqrt g^2\nonumber \\
!latex \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial \Ato{j,p}} \frac{\partial}{\partial \Ate{i,l}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!latex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!latex  ( - n_j n_i \bT{p,j}  \bT{l,i}  c_j s_i g_{\s\s} + n_j \bT{p,j}  \bT{l,i}' c_j c_i g_{\s\z}
!latex    -     n_i \bT{l,i}  \bT{p,j}' s_i s_j g_{\s\z} +     \bT{p,j}' \bT{l,i}' s_j c_i g_{\z\z}) / \sqrt g^2\nonumber \\
!latex \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial \Aze{j,p}} \frac{\partial}{\partial \Ate{i,l}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!latex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!latex  ( + m_j n_i \bT{p,j}  \bT{l,i}  s_j s_i g_{\s\s} - m_j \bT{p,j}  \bT{l,i}' s_j c_i g_{\s\z}
!latex    +     n_i \bT{l,i}  \bT{p,j}' s_i c_j g_{\s\t} -     \bT{p,j}' \bT{l,i}' c_j c_i g_{\t\z}) / \sqrt g^2\nonumber \\
!latex \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial \Azo{j,p}} \frac{\partial}{\partial \Ate{i,l}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!latex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!latex  ( - m_j n_i \bT{p,j}  \bT{l,i}  c_j s_i g_{\s\s} + m_j \bT{p,j}  \bT{l,i}' c_j c_i g_{\s\z}
!latex    +     n_i \bT{l,i}  \bT{p,j}' s_i s_j g_{\s\t} -     \bT{p,j}' \bT{l,i}' s_j c_i g_{\t\z}) / \sqrt g^2\nonumber
!latex \ee 
!latex \be
!latex \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial \Ate{j,p}} \frac{\partial}{\partial \Ato{i,l}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!latex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!latex  ( - n_j n_i \bT{p,j}  \bT{l,i}  s_j c_i g_{\s\s} - n_j \bT{p,j}  \bT{l,i}' s_j s_i g_{\s\z}
!latex    +     n_i \bT{l,i}  \bT{p,j}' c_i c_j g_{\s\z} +     \bT{p,j}' \bT{l,i}' c_j s_i g_{\z\z}) / \sqrt g^2\nonumber \\
!latex \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial \Ato{j,p}} \frac{\partial}{\partial \Ato{i,l}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!latex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!latex  ( + n_j n_i \bT{p,j}  \bT{l,i}  c_j c_i g_{\s\s} + n_j \bT{p,j}  \bT{l,i}' c_j s_i g_{\s\z}
!latex    +     n_i \bT{l,i}  \bT{p,j}' c_i s_j g_{\s\z} +     \bT{p,j}' \bT{l,i}' s_j s_i g_{\z\z}) / \sqrt g^2\nonumber \\
!latex \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial \Aze{j,p}} \frac{\partial}{\partial \Ato{i,l}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!latex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!latex  ( - m_j n_i \bT{p,j}  \bT{l,i}  s_j c_i g_{\s\s} - m_j \bT{p,j}  \bT{l,i}' s_j s_i g_{\s\z}
!latex    -     n_i \bT{l,i}  \bT{p,j}' c_i c_j g_{\s\t} -     \bT{p,j}' \bT{l,i}' c_j s_i g_{\t\z}) / \sqrt g^2\nonumber \\
!latex \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial \Azo{j,p}} \frac{\partial}{\partial \Ato{i,l}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!latex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!latex  ( + m_j n_i \bT{p,j}  \bT{l,i}  c_j c_i g_{\s\s} + m_j \bT{p,j}  \bT{l,i}' c_j s_i g_{\s\z}
!latex    -     n_i \bT{l,i}  \bT{p,j}' c_i s_j g_{\s\t} -     \bT{p,j}' \bT{l,i}' s_j s_i g_{\t\z}) / \sqrt g^2\nonumber   
!latex \ee 
!latex \be
!latex \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial \Ate{j,p}} \frac{\partial}{\partial \Aze{i,l}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!latex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!latex  ( + n_j m_i \bT{p,j}  \bT{l,i}  s_j s_i g_{\s\s} + n_j \bT{p,j}  \bT{l,i}' s_j c_i g_{\s\t}
!latex    -     m_i \bT{l,i}  \bT{p,j}' s_i c_j g_{\s\z} -     \bT{p,j}' \bT{l,i}' c_j c_i g_{\t\z}) / \sqrt g^2\nonumber \\
!latex \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial \Ato{j,p}} \frac{\partial}{\partial \Aze{i,l}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!latex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!latex  ( - n_j m_i \bT{p,j}  \bT{l,i}  c_j s_i g_{\s\s} - n_j \bT{p,j}  \bT{l,i}' c_j c_i g_{\s\t}
!latex    -     m_i \bT{l,i}  \bT{p,j}' s_i s_j g_{\s\z} -     \bT{p,j}' \bT{l,i}' s_j c_i g_{\t\z}) / \sqrt g^2\nonumber \\
!latex \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial \Aze{j,p}} \frac{\partial}{\partial \Aze{i,l}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!latex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!latex  ( + m_j m_i \bT{p,j}  \bT{l,i}  s_j s_i g_{\s\s} + m_j \bT{p,j}  \bT{l,i}' s_j c_i g_{\s\t}
!latex    +     m_i \bT{l,i}  \bT{p,j}' s_i c_j g_{\s\t} +     \bT{p,j}' \bT{l,i}' c_j c_i g_{\t\t}) / \sqrt g^2\nonumber \\
!latex \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial \Azo{j,p}} \frac{\partial}{\partial \Aze{i,l}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!latex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!latex  ( - m_j m_i \bT{p,j}  \bT{l,i}  c_j s_i g_{\s\s} - m_j \bT{p,j}  \bT{l,i}' c_j c_i g_{\s\t}
!latex    +     m_i \bT{l,i}  \bT{p,j}' s_i s_j g_{\s\t} +     \bT{p,j}' \bT{l,i}' s_j c_i g_{\t\t}) / \sqrt g^2\nonumber
!latex \ee 
!latex \be
!latex \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial \Ate{j,p}} \frac{\partial}{\partial \Azo{i,l}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!latex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!latex  ( - n_j m_i \bT{p,j}  \bT{l,i}  s_j c_i g_{\s\s} + n_j \bT{p,j}  \bT{l,i}' s_j s_i g_{\s\t}
!latex    +     m_i \bT{l,i}  \bT{p,j}' c_i c_j g_{\s\z} -     \bT{p,j}' \bT{l,i}' c_j s_i g_{\t\z}) / \sqrt g^2\nonumber \\
!latex \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial \Ato{j,p}} \frac{\partial}{\partial \Azo{i,l}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!latex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!latex  ( + n_j m_i \bT{p,j}  \bT{l,i}  c_j c_i g_{\s\s} - n_j \bT{p,j}  \bT{l,i}' c_j s_i g_{\s\t}
!latex    +     m_i \bT{l,i}  \bT{p,j}' c_i s_j g_{\s\z} -     \bT{p,j}' \bT{l,i}' s_j s_i g_{\t\z}) / \sqrt g^2\nonumber \\
!latex \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial \Aze{j,p}} \frac{\partial}{\partial \Azo{i,l}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!latex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 

!latex  ( - m_j m_i \bT{p,j}  \bT{l,i}  s_j c_i g_{\s\s} + m_j \bT{p,j}  \bT{l,i}' s_j s_i g_{\s\t}
!latex    -     m_i \bT{l,i}  \bT{p,j}' c_i c_j g_{\s\t} +     \bT{p,j}' \bT{l,i}' c_j s_i g_{\t\t}) / \sqrt g^2\nonumber \\
!latex \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial \Azo{j,p}} \frac{\partial}{\partial \Azo{i,l}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!latex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 

!latex  ( + m_j m_i \bT{p,j}  \bT{l,i}  c_j c_i g_{\s\s} - m_j \bT{p,j}  \bT{l,i}' c_j s_i g_{\s\t}
!latex    -     m_i \bT{l,i}  \bT{p,j}' c_i s_j g_{\s\t} +     \bT{p,j}' \bT{l,i}' s_j s_i g_{\t\t}) / \sqrt g^2\nonumber
!latex \ee

!latex \end{enumerate}

!latex \subsection{derivatives of helicity        integrals}

!latex \begin{enumerate}

!latex \item The first  derivatives of $\int \! dv \; {\bf A} \cdot{\bf B}$ with respect to $\Ate{i,l}$, $\Ato{i,l}$, $\Aze{i,l}$ and $\Azo{i,l}$ are
!latex \be 
!latex \frac{\partial}{\partial \Ate{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &
!latex \int \!\! dv \; \left( \frac{\partial {\bf A}}{\partial \Ate{i,l}} \cdot {\bf B} + {\bf A} \cdot \frac{\partial {\bf B}}{\partial \Ate{i,l}} \right) =
!latex \int \!\! dv \; ( \bT{l,i} \ci \nabla \t \cdot {\bf B} + {\bf A} \cdot \bT{l,i}' \ci \, {\bf e}_\z / \sqrt g ) \\
!latex \frac{\partial}{\partial \Ato{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &
!latex \int \!\! dv \; \left( \frac{\partial {\bf A}}{\partial \Ato{i,l}} \cdot {\bf B} + {\bf A} \cdot \frac{\partial {\bf B}}{\partial \Ato{i,l}} \right) =
!latex \int \!\! dv \; ( \bT{l,i} \si \nabla \t \cdot {\bf B} + {\bf A} \cdot \bT{l,i}' \si \, {\bf e}_\z / \sqrt g ) \\
!latex \frac{\partial}{\partial \Aze{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &
!latex \int \!\! dv \; \left( \frac{\partial {\bf A}}{\partial \Aze{i,l}} \cdot {\bf B} + {\bf A} \cdot \frac{\partial {\bf B}}{\partial \Aze{i,l}} \right) =
!latex \int \!\! dv \; ( \bT{l,i} \ci \nabla \z \cdot {\bf B} - {\bf A} \cdot \bT{l,i}' \ci \, {\bf e}_\t / \sqrt g ) \\
!latex \frac{\partial}{\partial \Azo{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &
!latex \int \!\! dv \; \left( \frac{\partial {\bf A}}{\partial \Azo{i,l}} \cdot {\bf B} + {\bf A} \cdot \frac{\partial {\bf B}}{\partial \Azo{i,l}} \right) =
!latex \int \!\! dv \; ( \bT{l,i} \si \nabla \z \cdot {\bf B} - {\bf A} \cdot \bT{l,i}' \si \, {\bf e}_\t / \sqrt g )
!latex \ee
!latex \item Note that in the above expressions, ${\bf A}\cdot{\bf e}_s=0$ has been used.
!latex \item The second derivatives of $\int \! dv \; {\bf A} \cdot{\bf B}$ with respect to $\Ate{i,l}$, $\Ato{i,l}$, $\Aze{i,l}$ and $\Azo{i,l}$ are
!latex \be 
!latex \frac{\partial}{\partial \Ate{j,p}} \frac{\partial}{\partial \Ate{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &   \int \!\! dv \; \left[
!latex \cancel{
!latex + \bT{l,i} \ci \nabla \t \cdot \bT{p,j}' \cj \, {\bf e}_\z}
!latex \cancel{
!latex + \bT{p,j} \cj \nabla \t \cdot \bT{l,i}' \ci \, {\bf e}_\z}
!latex \right] / \sqrt g \\
!latex \frac{\partial}{\partial \Ato{j,p}} \frac{\partial}{\partial \Ate{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &   \int \!\! dv \;  \left[
!latex \cancel{
!latex + \bT{l,i} \ci \nabla \t \cdot \bT{p,j}' \sj \, {\bf e}_\z}
!latex \cancel{
!latex + \bT{p,j} \sj \nabla \t \cdot \bT{l,i}' \ci \, {\bf e}_\z}
!latex \right] / \sqrt g \\
!latex \frac{\partial}{\partial \Aze{j,p}} \frac{\partial}{\partial \Ate{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &   \int \!\! dv \;  \left[
!latex - \bT{l,i} \ci \nabla \t \cdot \bT{p,j}' \cj \, {\bf e}_\t + \bT{p,j} \cj \nabla \z \cdot \bT{l,i}' \ci \, {\bf e}_\z \right] / \sqrt g \\
!latex \frac{\partial}{\partial \Azo{j,p}} \frac{\partial}{\partial \Ate{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &   \int \!\! dv \;  \left[
!latex - \bT{l,i} \ci \nabla \t \cdot \bT{p,j}' \sj \, {\bf e}_\t + \bT{p,j} \sj \nabla \z \cdot \bT{l,i}' \ci \, {\bf e}_\z \right] / \sqrt g \\
!latex \nonumber \\
!latex \frac{\partial}{\partial \Ate{j,p}} \frac{\partial}{\partial \Ato{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &   \int \!\! dv \;  \left[
!latex \cancel{
!latex + \bT{l,i} \si \nabla \t \cdot \bT{p,j}' \cj \, {\bf e}_\z}
!latex \cancel{
!latex + \bT{p,j} \cj \nabla \t \cdot \bT{l,i}' \si \, {\bf e}_\z}
!latex \right] / \sqrt g \\
!latex \frac{\partial}{\partial \Ato{j,p}} \frac{\partial}{\partial \Ato{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &   \int \!\! dv \;  \left[
!latex \cancel{
!latex + \bT{l,i} \si \nabla \t \cdot \bT{p,j}' \sj \, {\bf e}_\z}
!latex \cancel{
!latex + \bT{p,j} \sj \nabla \t \cdot \bT{l,i}' \si \, {\bf e}_\z}
!latex \right] / \sqrt g \\
!latex \frac{\partial}{\partial \Aze{j,p}} \frac{\partial}{\partial \Ato{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &   \int \!\! dv \;  \left[
!latex - \bT{l,i} \si \nabla \t \cdot \bT{p,j}' \cj \, {\bf e}_\t + \bT{p,j} \cj \nabla \z \cdot \bT{l,i}' \si \, {\bf e}_\z \right] / \sqrt g \\
!latex \frac{\partial}{\partial \Azo{j,p}} \frac{\partial}{\partial \Ato{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &   \int \!\! dv \;  \left[
!latex - \bT{l,i} \si \nabla \t \cdot \bT{p,j}' \sj \, {\bf e}_\t + \bT{p,j} \sj \nabla \z \cdot \bT{l,i}' \si \, {\bf e}_\z \right] / \sqrt g \\
!latex \nonumber \\
!latex \frac{\partial}{\partial \Ate{j,p}} \frac{\partial}{\partial \Aze{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &   \int \!\! dv \;  \left[
!latex + \bT{l,i} \ci \nabla \z \cdot \bT{p,j}' \cj \, {\bf e}_\z - \bT{p,j} \cj \nabla \t \cdot \bT{l,i}' \ci \, {\bf e}_\t \right] / \sqrt g \\
!latex \frac{\partial}{\partial \Ato{j,p}} \frac{\partial}{\partial \Aze{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &   \int \!\! dv \;  \left[
!latex + \bT{l,i} \ci \nabla \z \cdot \bT{p,j}' \sj \, {\bf e}_\z - \bT{p,j} \sj \nabla \t \cdot \bT{l,i}' \ci \, {\bf e}_\t \right] / \sqrt g \\
!latex \frac{\partial}{\partial \Aze{j,p}} \frac{\partial}{\partial \Aze{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &   \int \!\! dv \;  \left[
!latex \cancel{
!latex - \bT{l,i} \ci \nabla \z \cdot \bT{p,j}' \cj \, {\bf e}_\t}
!latex \cancel{
!latex - \bT{p,j} \cj \nabla \z \cdot \bT{l,i}' \ci \, {\bf e}_\t}
!latex \right] / \sqrt g \\
!latex \frac{\partial}{\partial \Azo{j,p}} \frac{\partial}{\partial \Aze{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &   \int \!\! dv \;  \left[
!latex \cancel{
!latex - \bT{l,i} \ci \nabla \z \cdot \bT{p,j}' \sj \, {\bf e}_\t}
!latex \cancel{ma00aa
!latex - \bT{l,i} \si \nabla \z \cdot \bT{p,j}' \sj \, {\bf e}_\t}
!latex \cancel{
!latex - \bT{p,j} \sj \nabla \z \cdot \bT{l,i}' \si \, {\bf e}_\t}
!latex \right] / \sqrt g
!latex \ee
!latex \item In these expressions the terms $        \nabla \t \cdot {\bf e}_\t =        \nabla \z \cdot {\bf e}_\z =1$, 
!latex                                  and $\cancel{\nabla \t \cdot {\bf e}_\z}=\cancel{\nabla \z \cdot {\bf e}_\t}=0$ 
!latex       have been included to show the structure of the derivation.

!latex \end{enumerate}

!latex \subsection{derivatives of kinetic energy integrals}

!latex \begin{enumerate}

!latex \item The first derivatives of $\int \! dv \; v^2$ with respect to $\vse{i,l}$ etc. are
!latex \be \frac{\partial}{\partial \vse{i,l}} \int \!\! dv \; {\bf v} \cdot {\bf v} & = & 2 \int \!\! dv \; {\bf v} \cdot \bT{l,i} \ci \nabla \s \\
!latex     \frac{\partial}{\partial \vso{i,l}} \int \!\! dv \; {\bf v} \cdot {\bf v} & = & 2 \int \!\! dv \; {\bf v} \cdot \bT{l,i} \si \nabla \s \\
!latex     \frac{\partial}{\partial \vte{i,l}} \int \!\! dv \; {\bf v} \cdot {\bf v} & = & 2 \int \!\! dv \; {\bf v} \cdot \bT{l,i} \ci \nabla \t \\
!latex     \frac{\partial}{\partial \vto{i,l}} \int \!\! dv \; {\bf v} \cdot {\bf v} & = & 2 \int \!\! dv \; {\bf v} \cdot \bT{l,i} \si \nabla \t \\
!latex     \frac{\partial}{\partial \vze{i,l}} \int \!\! dv \; {\bf v} \cdot {\bf v} & = & 2 \int \!\! dv \; {\bf v} \cdot \bT{l,i} \ci \nabla \z \\
!latex     \frac{\partial}{\partial \vzo{i,l}} \int \!\! dv \; {\bf v} \cdot {\bf v} & = & 2 \int \!\! dv \; {\bf v} \cdot \bT{l,i} \si \nabla \z \\
!latex \ee

!latex \end{enumerate}

!latex \subsection{calculation of volume-integrated basis-function-weighted metric information}

!latex \begin{enumerate}

!latex \item The required geometric information is calculated in \link{ma00aa}.

!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine matrix( lvol, mn, lrad )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, one, two, half, pi2
  
  use numerical, only : small
  
  use fileunits, only : ounit
  
  use inputlist, only : Wmacros, Wmatrix, mpol, Nvol, Lvcvacuum, Cteta, Czeta
  
  use cputiming, only : Tmatrix
  
  use allglobal, only : ncpu, myid, cpus, &
                        YESstellsym, NOTstellsym, &
                        im, in, &
                        NAdof, &
                        dMA, dMD, dMB, dMG, &
                        Lvacuumregion, iRbc, iZbs, Dxyz, iZbs, dVC, &
!                       VRij, VZij, Vsg, Vguvij, &
                        Ctz, kjzeta, kjicos, kjisin, & ! 03/03/21 ;
                        Dxyz, Nxyz, virtualcasingfactor, &
                        Ctz, Ntz, Nt, Nz, Mvol, &
!                       Ateij, Azeij, & ! 03/03/21 ;
                        efmn, ofmn, cfmn, sfmn, & ! 03/03/21 ;
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
  
  INTEGER              :: jk, kj, ifail ! 03/03/21 ;

  REAL                 :: Wtete, Wteto, Wtote, Wtoto
  REAL                 :: Wteze, Wtezo, Wtoze, Wtozo
  REAL                 :: Wzete, Wzeto, Wzote, Wzoto
  REAL                 :: Wzeze, Wzezo, Wzoze, Wzozo
  
  REAL                 :: Htete, Hteto, Htote, Htoto
  REAL                 :: Hteze, Htezo, Htoze, Htozo
  REAL                 :: Hzete, Hzeto, Hzote, Hzoto
  REAL                 :: Hzeze, Hzezo, Hzoze, Hzozo
  
  REAL,allocatable     :: TTdata(:,:,:), TTMdata(:,:) ! queues to construct sparse matrices

  REAL                 :: VRij(1:Ctz,0:3), VZij(1:Ctz,0:3), Vsg(1:Ctz), Vguvij(1:Ctz,1:3,1:3) ! 03/03/21 ;
  REAL                 :: rx(1:Ctz), ry(1:Ctz), rz(1:Ctz), rrr(1:Ctz) ! 03/03/21 ;
  REAL                 :: xtrr(1:Ctz,1:3), xzrr(1:Ctz,1:3), xtrrdotds(1:Ctz), xzrrdotds(1:Ctz) ! 03/03/21 ;
  REAL                 :: Atejkj(1:Ntz,1:mn), Azejkj(1:Ntz,1:mn) ! 03/03/21 ;
  REAL                 :: Ateij(1:mn,1:mn), Azeij(1:mn,1:mn)

  BEGIN(matrix)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  FATAL( matrix, .not.allocated(dMA), error )
  FATAL( matrix, .not.allocated(dMD), error )
  FATAL( matrix, .not.allocated(dMB), error )
  FATAL( matrix, .not.allocated(dMG), error )  
  if( Lvacuumregion .and. Lvcvacuum.eq.1 ) then
  FATAL( matrix, .not.allocated(dVC), error )  
  endif
#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  NN = NAdof(lvol) ! shorthand;
  
  dMA(0:NN,0:NN) = zero
  dMD(0:NN,0:NN) = zero

  if( Lvacuumregion .and. Lvcvacuum.eq.1 ) dVC(0:NN,0:NN) = zero

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Lvacuumregion .and. Lvcvacuum.eq.1 ) then ! 03/03/21 ;
   
   write(ounit,'("matrix : " 10x " : dBdX: L =",L2," ; vol =",i3," ; innout =",i2," ; ii =",i3," ; irz =",i2," ; issym =",i2," ;")') &
                                     dBdX%L, dBdX%vol, dBdX%innout, dBdX%ii, dBdX%irz, dBdX%issym

!  pause

   VRij(1:Ctz,0:3) = zero ! initialize summation ; 03/03/21 ;
   VZij(1:Ctz,0:3) = zero
   
   do kk = 0, Czeta-1
    do jj = 0, Cteta-1 ; kj = 1 + jj + kk * Cteta ! construct high-resolution plasma boundary (for virtual-casing integral) ; 03/03/21 ;
     
     do ii = 1, mn
      VRij(kj,0) = VRij(kj,0) +   iRbc(ii,Nvol)                   * kjicos(kj,ii) 
      VRij(kj,1) = VRij(kj,1) + ( iRbc(ii,Mvol) - iRbc(ii,Nvol) ) * kjicos(kj,ii) * half
      VRij(kj,2) = VRij(kj,2) -   iRbc(ii,Nvol)                   * kjisin(kj,ii) * ( + im(ii) )
      VRij(kj,3) = VRij(kj,3) -   iRbc(ii,Nvol)                   * kjisin(kj,ii) * ( - in(ii) )
      VZij(kj,0) = VZij(kj,0) +   iZbs(ii,Nvol)                   * kjisin(kj,ii)
      VZij(kj,1) = VZij(kj,1) + ( iZbs(ii,Mvol) - iZbs(ii,Nvol) ) * kjisin(kj,ii) * half
      VZij(kj,2) = VZij(kj,2) +   iZbs(ii,Nvol)                   * kjicos(kj,ii) * ( + im(ii) )
      VZij(kj,3) = VZij(kj,3) +   iZbs(ii,Nvol)                   * kjicos(kj,ii) * ( - in(ii) )
     enddo
     
    enddo ! end of do kk ; 07/29/20;
   enddo ! end of do jj ; 07/29/20;
   
   Vsg(1:Ctz) = VRij(1:Ctz,0) * ( VZij(1:Ctz,1)*VRij(1:Ctz,2) - VRij(1:Ctz,1)*VZij(1:Ctz,2) ) ; Vsg = one / Vsg
   do ii = 1, 3
    do jj = ii, 3 ; Vguvij(1:Ctz,ii,jj) = VRij(1:Ctz,ii) * VRij(1:Ctz,jj) + VZij(1:Ctz,ii) * VZij(1:Ctz,jj)
    enddo
   enddo
   Vguvij(1:Ctz,3,3) = Vguvij(1:Ctz,3,3) + VRij(1:Ctz,0) * VRij(1:Ctz,0)
   
   do jk = 1, Ntz ! loop over plasma boundary (to construct Fourier harmonics of B_{plasma} \cdot {\bf n} ; 03/03/21 ;
    
    rx(1:Ctz) = Dxyz(1,jk) - VRij(1:Ctz,0) * cos(kjzeta(1:Ctz))
    ry(1:Ctz) = Dxyz(2,jk) - VRij(1:Ctz,0) * sin(kjzeta(1:Ctz))
    rz(1:Ctz) = Dxyz(3,jk) - VZij(1:Ctz,0)
    
    rrr(1:Ctz) = ( sqrt( rx(1:Ctz) * rx(1:Ctz) + ry(1:Ctz) * ry(1:Ctz) + rz(1:Ctz) * rz(1:Ctz) ) )**3 ; rrr(1:Ctz) = one / rrr(1:Ctz)
    
    xtrr(1:Ctz,1) = ( VRij(1:Ctz,2)*sin(kjzeta(1:Ctz))                                    ) * rz(1:Ctz) & ! cross product \bx_\t \times \br ;
                  - ( VZij(1:Ctz,2)                                                       ) * ry(1:Ctz)
    xtrr(1:Ctz,2) = ( VZij(1:Ctz,2)                                                       ) * rx(1:Ctz) &
                  - ( VRij(1:Ctz,2)*cos(kjzeta(1:Ctz))                                    ) * rz(1:Ctz)  
    xtrr(1:Ctz,3) = ( VRij(1:Ctz,2)*cos(kjzeta(1:Ctz))                                    ) * ry(1:Ctz) &
                  - ( VRij(1:Ctz,2)*sin(kjzeta(1:Ctz))                                    ) * rx(1:Ctz)
  
    xzrr(1:Ctz,1) = ( VRij(1:Ctz,3)*sin(kjzeta(1:Ctz)) + VRij(1:Ctz,0)*cos(kjzeta(1:Ctz)) ) * rz(1:Ctz) & ! cross product \bx_\z \times \br ;
                  - ( VZij(1:Ctz,3)                                                       ) * ry(1:Ctz)
    xzrr(1:Ctz,2) = ( VZij(1:Ctz,3)                                                       ) * rx(1:Ctz) &
                  - ( VRij(1:Ctz,3)*cos(kjzeta(1:Ctz)) - VRij(1:Ctz,0)*sin(kjzeta(1:Ctz)) ) * rz(1:Ctz)
    xzrr(1:Ctz,3) = ( VRij(1:Ctz,3)*cos(kjzeta(1:Ctz)) - VRij(1:Ctz,0)*sin(kjzeta(1:Ctz)) ) * ry(1:Ctz) &
                  - ( VRij(1:Ctz,3)*sin(kjzeta(1:Ctz)) + VRij(1:Ctz,0)*cos(kjzeta(1:Ctz)) ) * rx(1:Ctz)
     
    xtrrdotds(1:Ctz) = xtrr(1:Ctz,1) * Nxyz(1,jk) + xtrr(1:Ctz,2) * Nxyz(2,jk) + xtrr(1:Ctz,3) * Nxyz(3,jk)
    xzrrdotds(1:Ctz) = xzrr(1:Ctz,1) * Nxyz(1,jk) + xzrr(1:Ctz,2) * Nxyz(2,jk) + xzrr(1:Ctz,3) * Nxyz(3,jk)
    
    do jj = 1, mn
     Atejkj(jk,jj) = + sum( kjicos(1:Ctz,jj) * ( Vguvij(1:Ctz,3,3) * xtrrdotds(1:Ctz) - Vguvij(1:Ctz,2,3) * xzrrdotds(1:Ctz) ) * rrr(1:Ctz) * Vsg(1:Ctz) )
     Azejkj(jk,jj) = - sum( kjicos(1:Ctz,jj) * ( Vguvij(1:Ctz,2,3) * xtrrdotds(1:Ctz) - Vguvij(1:Ctz,2,2) * xzrrdotds(1:Ctz) ) * rrr(1:Ctz) * Vsg(1:Ctz) )  
    enddo ! end of do jj ; 07/29/20;
    
   enddo ! end of do jk ; 07/29/20;
   
   Atejkj(1:Ntz,1:mn) = Atejkj(1:Ntz,1:mn) * pi2 * pi2 * virtualcasingfactor / Ctz
   Azejkj(1:Ntz,1:mn) = Azejkj(1:Ntz,1:mn) * pi2 * pi2 * virtualcasingfactor / Ctz
   
   do jj = 1, mn
    call tfft( Nt, Nz, Atejkj(1:Ntz,jj), Azejkj(1:Ntz,jj), mn, im(1:mn), in(1:mn), efmn(1:mn), Ateij(1:mn,jj), cfmn(1:mn), Azeij(1:mn,jj), ifail )     
   enddo
   
   do ii = 1, mn
    
    do ll = 0, lrad ! Chebyshev polynomial ;
     
     if( ii.gt.1 ) then ; id = Lme(lvol,  ii)
      do jj = 1, mn     ; jd = Ate(lvol,0,jj)%i(ll) ; dVC(id,jd) = Ateij(ii,jj) * TT(ll,0,1)
       ;                ; jd = Aze(lvol,0,jj)%i(ll) ; dVC(id,jd) = Azeij(ii,jj) * TT(ll,0,1)
      enddo
     endif
     
    enddo ! end of do ll ;
    
   enddo ! end of do ii ;
    
  endif ! end of if( Lvacuumregion .and. Lvcvacuum.eq.1 ) ! 03/03/21 ;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

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

  if( Lvacuumregion .and. Lvcvacuum.eq.1 ) dMA(0:NN,0:NN) = dMA(0:NN,0:NN) - dVC(0:NN,0:NN) ! 03/03/21 ;

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


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine matrixBG( lvol, mn, lrad )

  ! only compute the dMB and dMG matrix for matrix-free mode

  use constants, only : zero, one

  use fileunits, only : ounit

  use inputlist, only : Wmatrix, Lvcvacuum

  use cputiming, only : Tmatrix

  use allglobal, only : myid, cpus, NAdof, im, in,&
                        dMG, dMB, YESstellsym, &
                        iVnc, iVns, iBnc, iBns, &
                        Lme, Lmf, Lmg, Lmh, &
                        Lvacuumregion

  LOCALS

  INTEGER, intent(in)  :: lvol, mn, lrad

  INTEGER :: NN, ii, id, mi, ni

  BEGIN( matrix )

  NN = NAdof(lvol) ! shorthand;
  
  dMB(0:NN,1: 2) = zero
  dMG(0:NN     ) = zero
  
  if( YESstellsym ) then
   
   select case( Lvcvacuum ) 
    
   case( 1 ) ! Lvcvacuum.eq.1 ; 03/03/21 ;
    
    do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
     
     ;if( ii.gt.1 ) then ; id = Lme(lvol,  ii)       ;                           ; dMG(id   ) = - ( iVns(ii)            )
     ;else               ; id = Lmg(lvol,  ii)       ;                           ; dMB(id, 1) = -       one
!    ;                   ; id = Lmh(lvol,  ii)       ;                           ; dMB(id, 2) = -       one ! to be deleted;
     ;                   ; id = Lmh(lvol,  ii)       ;                           ; dMB(id, 2) = +       one ! changed sign;
     ;endif
     
    enddo ! end of do ii ;
    
   case default ! Lvcvacuum ; 03/03/21 ;
    
    do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
     
     ;if( ii.gt.1 ) then ; id = Lme(lvol,  ii)       ;                           ; dMG(id   ) = - ( iVns(ii) + iBns(ii) )
     ;else               ; id = Lmg(lvol,  ii)       ;                           ; dMB(id, 1) = -       one
!    ;                   ; id = Lmh(lvol,  ii)       ;                           ; dMB(id, 2) = -       one ! to be deleted;
     ;                   ; id = Lmh(lvol,  ii)       ;                           ; dMB(id, 2) = +       one ! changed sign;
     ;endif
     
    enddo ! end of do ii ;
    
   end select ! end of select case( Lvcvacuum ) ; 03/03/21 ;
   
  else ! NOTstellsym ;
   
   select case( Lvcvacuum ) 
    
   case( 1 ) ! Lvcvacuum.eq.1 ; 03/03/21 ;
    
    FATAL( matrix, .true., Lvcvacuum not implemented for non-stellarator-symmetry )
    
   case default ! Lvcvacuum ; 03/03/21 ;
    
    do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
     
     ;if( ii.gt.1 ) then ; id = Lme(lvol,ii)         ;                           ; dMG(id   ) = - ( iVns(ii) + iBns(ii) )
     ;                   ; id = Lmf(lvol,ii)         ;                           ; dMG(id   ) = - ( iVnc(ii) + iBnc(ii) )
     ;else               ; id = Lmg(lvol,ii)         ;                           ; dMB(id, 1) = -       one
!    ;                   ; id = Lmh(lvol,ii)         ;                           ; dMB(id, 2) = -       one ! to be deleted;
     ;                   ; id = Lmh(lvol,ii)         ;                           ; dMB(id, 2) = +       one ! changed sign;
     ;endif
    
    enddo ! end of do ii ;
    
   end select ! end of select case( Lvcvacuum ) ; 03/03/21 ;

   endif ! end of if( YESstellsym ) ;

  RETURN( matrix )

end subroutine matrixBG

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
