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

!l tex \subsection{enclosed currents} 

!l tex \begin{enumerate}
!l tex \item In the vacuum region, the enclosed currents must be constrained.
!l tex \item The plasma current is
!l tex       \be I \equiv \int_{\cal S} {\bf j}\cdot d{\bf s} = \int_{\partial \cal S} {\bf B}\cdot d{\bf l} = \int_{0}^{2\pi} {\bf B}\cdot {\bf e}_\t \, d\t
!l tex       \ee
!l tex \item Choosing to take the line integral to lie on the inner surface (the plasma boundary), where $B^s=0$, this is
!l tex       \be I = \int_{0}^{2\pi} \left( - \partial_s A_\z \; \bar g_{\t\t} + \partial_s A_\t \; \bar g_{\t\z} \right) \, d\t,
!l tex       \ee
!l tex       where $\bar g_{\mu\nu} \equiv g_{\mu\nu} / \sqrt g$.
!l tex \item The ``linking'' current through the torus is
!l tex       \be G \equiv \int_{\cal S} {\bf j}\cdot d{\bf s} = \int_{\partial \cal S} {\bf B}\cdot d{\bf l} = \int_{0}^{2\pi} {\bf B}\cdot {\bf e}_\z \, d\z
!l tex       \ee
!l tex \item Choosing to take the line integral to lie on the inner surface (the plasma boundary), where $B^s=0$, this is
!l tex       \be G = \int_{0}^{2\pi} \left( -\partial_s A_\z \; \bar g_{\t\z} + \partial_s A_\t \; \bar g_{\z\z} \right) \, d\z.
!l tex       \ee
!l tex \end{enumerate} 

!latex \subsection{Fourier-Chebyshev representation} 

!latex \begin{enumerate}
!latex \item \Ais %mponents of the vector potential are
!l tex       \be A_\t(\s,\t,\z) &=& \sum_{i,l} \Ate{i,l} \; \bT{l,i}(s) \ci + \sum_{i,l} \Ato{i,l} \; \bT{l,i}(s) \si, \\
!l tex           A_\z(\s,\t,\z) &=& \sum_{i,l} \Aze{i,l} \; \bT{l,i}(s) \ci + \sum_{i,l} \Azo{i,l} \; \bT{l,i}(s) \si,
!l tex       \ee
!l tex       where $\bT{l,i}(s) \equiv s^{m_i/2} \, T_l(s)$ and $\a_j \equiv m_j\t-n_j\z$.
!l tex \item (The regularity factor, $s^{m_i/2}$, is only required if the coordinate singularity is in the domain.)
!latex \item \Bis % netic field is then
!l tex       \be \begin{array}{ccccrcrcrcrcccccccccccccccccccccccccccccccccccccccccccccccccccc}
!l tex       \!\!\!\!\!\!\!\!\!\! \sqrt g \, {\bf B} & = & {\bf e}_\s & \sum_{i,l} [ 
!l tex       ( & - m_i \Aze{i,l} & - & n_i \Ate{i,l} & ) \bT{l,i}  \si + ( & + m_i \Azo{i,l} & + & n_i \Ato{i,l} & ) \bT{l,i} \ci
!l tex       ] \\
!l tex                          & + & {\bf e}_\t & \sum_{i,l} [ 
!l tex       ( &                 & - &     \Aze{i,l} & ) \bT{l,i}' \ci + ( &                 & - &     \Azo{i,l} & ) \bT{l,i}' \si
!l tex       ] \\
!l tex                          & + & {\bf e}_\z & \sum_{i,l} [ 
!l tex       ( &       \Ate{i,l} &   &               & ) \bT{l,i}' \ci + ( &       \Ato{i,l} &   &               & ) \bT{l,i}' \si
!l tex       ] 
!l tex       \end{array} \ee
!latex \item The components of the velocity, ${\bf v} \equiv v_\s \nabla \s + v_\t \nabla \t + v_\z \nabla \z$,  are
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
!latex and is given by:
!latex \be \begin{array}{cclcclcccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc}                                                   
!latex     {\cal F}
!latex                 & \equiv & \ds                        \int {\bf B} \cdot {\bf B} \, dv                                                              
!latex                   +                                   \int {\bf v} \cdot {\bf v} \, dv                                                                
!latex                    -       \ds       \mu       \left[ \int {\bf A} \cdot {\bf B} \, dv - K \right]                                                  \\
!latex                 &  + &   & \ds \sum_{i=1} & a_i       & \ds \left[  \sum_l \Ate{i,l} T_l(-1) - 0             \right]                                \\
!latex                 &  + &   & \ds \sum_{i=1} & b_i       & \ds \left[  \sum_l \Aze{i,l} T_l(-1) - 0             \right]                                \\
!latex                 &  + &   & \ds \sum_{i=2} & c_i       & \ds \left[  \sum_l \Ato{i,l} T_l(-1) - 0             \right]                                \\
!latex                 &  + &   & \ds \sum_{i=2} & d_i       & \ds \left[  \sum_l \Azo{i,l} T_l(-1) - 0             \right]                                \\
!latex                 &  + &   & \ds \sum_{i=2} & e_i       & \ds \left[  \sum_l \left( - m_i \Aze{i,l} - n_i \Ate{i,l} \right) T_l(+1) - b_{s,i} \right] \\
!latex                 &  + &   & \ds \sum_{i=2} & f_i       & \ds \left[  \sum_l \left( + m_i \Azo{i,l} + n_i \Ato{i,l} \right) T_l(+1) - b_{c,i} \right] \\
!latex                 &  + &   & \ds            & g_1       & \ds \left[  \sum_l \Ate{1,l} T_l(+1) - \Delta \psi_t \right]                                \\
!latex                 &  + &   & \ds            & h_1       & \ds \left[  \sum_l \Aze{1,l} T_l(+1) - \Delta \psi_p \right]                             
!l tex                 &  + &   & \ds            & \gamma    & \ds \left[ \int_{0}^{2\pi} {\bf B} \cdot {\bf e}_\t \, d\t - I \right]                      \\
!l tex                 &  + &   & \ds            & \delta    & \ds \left[ \int_{0}^{2\pi} {\bf B} \cdot {\bf e}_\z \, d\z - G \right],
!latex     \end{array}
!latex \ee   
!latex       where
!latex       \bi
!latex       \item[  i.] $a_i$, $b_i$, $c_i$ and $d_i$ are Lagrange multipliers used to enforce the combined gauge and interface boundary condition
!latex                   on the inner interface,
!latex       \item[ ii.] $e_i$ and $f_i$               are Lagrange multipliers used to enforce the interface boundary condition on the outer interface, 
!latex                   namely \mbox{$\sqrt g \, {\bf B}\cdot\nabla s = b$}; and
!latex       \item[iii.] $g_1$    and $h_1$            are Lagrange multipliers used to enforce the constraints on the enclosed fluxes.
!l tex       \item[ iv.] $\gamma$ and $\delta$         are Lagrange multipliers used to enforce the condition on the enclosed currents.
!latex       \ei
!latex \item In each plasma volume
!latex           the boundary condition on the outer interface is $b=0$.
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
!l tex ( - n_j \bT{p,j}  \sj \, {\bf e}_\s + \bT{p,j}' \cj \, {\bf e}_\z ) \cdot
!l tex ( - n_i \bT{l,i}  \si \, {\bf e}_\s + \bT{l,i}' \ci \, {\bf e}_\z ) / \sqrt g^2 \nonumber \\
!l tex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \;
!latex  ( + n_j n_i \bT{p,j}  \bT{l,i}  s_j s_i g_{\s\s} - n_j \bT{p,j}  \bT{l,i}' s_j c_i g_{\s\z}
!latex    -     n_i \bT{l,i}  \bT{p,j}' s_i c_j g_{\s\z} +     \bT{p,j}' \bT{l,i}' c_j c_i g_{\z\z}) / \sqrt g^2\nonumber \\
!latex \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial \Ato{j,p}} \frac{\partial}{\partial \Ate{i,l}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!latex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!l tex ( + n_j \bT{p,j}  \cj \, {\bf e}_\s + \bT{p,j}' \sj \, {\bf e}_\z ) \cdot
!l tex ( - n_i \bT{l,i}  \si \, {\bf e}_\s + \bT{l,i}' \ci \, {\bf e}_\z ) / \sqrt g^2 \nonumber \\
!l tex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \;
!latex  ( - n_j n_i \bT{p,j}  \bT{l,i}  c_j s_i g_{\s\s} + n_j \bT{p,j}  \bT{l,i}' c_j c_i g_{\s\z}
!latex    -     n_i \bT{l,i}  \bT{p,j}' s_i s_j g_{\s\z} +     \bT{p,j}' \bT{l,i}' s_j c_i g_{\z\z}) / \sqrt g^2\nonumber \\
!latex \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial \Aze{j,p}} \frac{\partial}{\partial \Ate{i,l}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!latex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!l tex ( - m_j \bT{p,j}  \sj \, {\bf e}_\s - \bT{p,j}' \cj \, {\bf e}_\t ) \cdot
!l tex ( - n_i \bT{l,i}  \si \, {\bf e}_\s + \bT{l,i}' \ci \, {\bf e}_\z ) / \sqrt g^2 \nonumber \\
!l tex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \;
!latex  ( + m_j n_i \bT{p,j}  \bT{l,i}  s_j s_i g_{\s\s} - m_j \bT{p,j}  \bT{l,i}' s_j c_i g_{\s\z}
!latex    +     n_i \bT{l,i}  \bT{p,j}' s_i c_j g_{\s\t} -     \bT{p,j}' \bT{l,i}' c_j c_i g_{\t\z}) / \sqrt g^2\nonumber \\
!latex \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial \Azo{j,p}} \frac{\partial}{\partial \Ate{i,l}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!latex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!l tex ( + m_j \bT{p,j}  \cj \, {\bf e}_\s - \bT{p,j}' \sj \, {\bf e}_\t ) \cdot
!l tex ( - n_i \bT{l,i}  \si \, {\bf e}_\s + \bT{l,i}' \ci \, {\bf e}_\z ) / \sqrt g^2 \nonumber \\
!l tex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \;
!latex  ( - m_j n_i \bT{p,j}  \bT{l,i}  c_j s_i g_{\s\s} + m_j \bT{p,j}  \bT{l,i}' c_j c_i g_{\s\z}
!latex    +     n_i \bT{l,i}  \bT{p,j}' s_i s_j g_{\s\t} -     \bT{p,j}' \bT{l,i}' s_j c_i g_{\t\z}) / \sqrt g^2\nonumber   
!latex \ee 
!latex \be
!latex \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial \Ate{j,p}} \frac{\partial}{\partial \Ato{i,l}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!latex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!l tex ( - n_j \bT{p,j}  \sj \, {\bf e}_\s + \bT{p,j}' \cj \, {\bf e}_\z ) \cdot
!l tex ( + n_i \bT{l,i}  \ci \, {\bf e}_\s + \bT{l,i}' \si \, {\bf e}_\z ) / \sqrt g^2 \nonumber \\
!l tex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \;
!latex  ( - n_j n_i \bT{p,j}  \bT{l,i}  s_j c_i g_{\s\s} - n_j \bT{p,j}  \bT{l,i}' s_j s_i g_{\s\z}
!latex    +     n_i \bT{l,i}  \bT{p,j}' c_i c_j g_{\s\z} +     \bT{p,j}' \bT{l,i}' c_j s_i g_{\z\z}) / \sqrt g^2\nonumber \\
!latex \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial \Ato{j,p}} \frac{\partial}{\partial \Ato{i,l}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!latex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!l tex ( + n_j \bT{p,j}  \cj \, {\bf e}_\s + \bT{p,j}' \sj \, {\bf e}_\z ) \cdot
!l tex ( + n_i \bT{l,i}  \ci \, {\bf e}_\s + \bT{l,i}' \si \, {\bf e}_\z ) / \sqrt g^2 \nonumber \\
!l tex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \;
!latex  ( + n_j n_i \bT{p,j}  \bT{l,i}  c_j c_i g_{\s\s} + n_j \bT{p,j}  \bT{l,i}' c_j s_i g_{\s\z}
!latex    +     n_i \bT{l,i}  \bT{p,j}' c_i s_j g_{\s\z} +     \bT{p,j}' \bT{l,i}' s_j s_i g_{\z\z}) / \sqrt g^2\nonumber \\
!latex \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial \Aze{j,p}} \frac{\partial}{\partial \Ato{i,l}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!latex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!l tex ( - m_j \bT{p,j}  \sj \, {\bf e}_\s - \bT{p,j}' \cj \, {\bf e}_\t ) \cdot
!l tex ( + n_i \bT{l,i}  \ci \, {\bf e}_\s + \bT{l,i}' \si \, {\bf e}_\z ) / \sqrt g^2 \nonumber \\
!l tex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \;
!latex  ( - m_j n_i \bT{p,j}  \bT{l,i}  s_j c_i g_{\s\s} - m_j \bT{p,j}  \bT{l,i}' s_j s_i g_{\s\z}
!latex    -     n_i \bT{l,i}  \bT{p,j}' c_i c_j g_{\s\t} -     \bT{p,j}' \bT{l,i}' c_j s_i g_{\t\z}) / \sqrt g^2\nonumber \\
!latex \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial \Azo{j,p}} \frac{\partial}{\partial \Ato{i,l}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!latex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!l tex ( + m_j \bT{p,j}  \cj \, {\bf e}_\s - \bT{p,j}' \sj \, {\bf e}_\t ) \cdot
!l tex ( + n_i \bT{l,i}  \ci \, {\bf e}_\s + \bT{l,i}' \si \, {\bf e}_\z ) / \sqrt g^2 \nonumber \\
!l tex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \;
!latex  ( + m_j n_i \bT{p,j}  \bT{l,i}  c_j c_i g_{\s\s} + m_j \bT{p,j}  \bT{l,i}' c_j s_i g_{\s\z}
!latex    -     n_i \bT{l,i}  \bT{p,j}' c_i s_j g_{\s\t} -     \bT{p,j}' \bT{l,i}' s_j s_i g_{\t\z}) / \sqrt g^2\nonumber   
!latex \ee 
!latex \be
!latex \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial \Ate{j,p}} \frac{\partial}{\partial \Aze{i,l}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!latex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!l tex ( - n_j \bT{p,j}  \sj \, {\bf e}_\s + \bT{p,j}' \cj \, {\bf e}_\z ) \cdot
!l tex ( - m_i \bT{l,i}  \si \, {\bf e}_\s - \bT{l,i}' \ci \, {\bf e}_\t ) / \sqrt g^2 \nonumber \\
!l tex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \;
!latex  ( + n_j m_i \bT{p,j}  \bT{l,i}  s_j s_i g_{\s\s} + n_j \bT{p,j}  \bT{l,i}' s_j c_i g_{\s\t}
!latex    -     m_i \bT{l,i}  \bT{p,j}' s_i c_j g_{\s\z} -     \bT{p,j}' \bT{l,i}' c_j c_i g_{\t\z}) / \sqrt g^2\nonumber \\
!latex \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial \Ato{j,p}} \frac{\partial}{\partial \Aze{i,l}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!latex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!l tex ( + n_j \bT{p,j}  \cj \, {\bf e}_\s + \bT{p,j}' \sj \, {\bf e}_\z ) \cdot
!l tex ( - m_i \bT{l,i}  \si \, {\bf e}_\s - \bT{l,i}' \ci \, {\bf e}_\t ) / \sqrt g^2 \nonumber \\
!l tex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \;
!latex  ( - n_j m_i \bT{p,j}  \bT{l,i}  c_j s_i g_{\s\s} - n_j \bT{p,j}  \bT{l,i}' c_j c_i g_{\s\t}
!latex    -     m_i \bT{l,i}  \bT{p,j}' s_i s_j g_{\s\z} -     \bT{p,j}' \bT{l,i}' s_j c_i g_{\t\z}) / \sqrt g^2\nonumber \\
!latex \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial \Aze{j,p}} \frac{\partial}{\partial \Aze{i,l}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!latex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!l tex ( - m_j \bT{p,j}  \sj \, {\bf e}_\s - \bT{p,j}' \cj \, {\bf e}_\t ) \cdot
!l tex ( - m_i \bT{l,i}  \si \, {\bf e}_\s - \bT{l,i}' \ci \, {\bf e}_\t ) / \sqrt g^2 \nonumber \\
!l tex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \;
!latex  ( + m_j m_i \bT{p,j}  \bT{l,i}  s_j s_i g_{\s\s} + m_j \bT{p,j}  \bT{l,i}' s_j c_i g_{\s\t}
!latex    +     m_i \bT{l,i}  \bT{p,j}' s_i c_j g_{\s\t} +     \bT{p,j}' \bT{l,i}' c_j c_i g_{\t\t}) / \sqrt g^2\nonumber \\
!latex \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial \Azo{j,p}} \frac{\partial}{\partial \Aze{i,l}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!latex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!l tex ( + m_j \bT{p,j}  \cj \, {\bf e}_\s - \bT{p,j}' \sj \, {\bf e}_\t ) \cdot
!l tex ( - m_i \bT{l,i}  \si \, {\bf e}_\s - \bT{l,i}' \ci \, {\bf e}_\t ) / \sqrt g^2 \nonumber \\
!l tex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \;
!latex  ( - m_j m_i \bT{p,j}  \bT{l,i}  c_j s_i g_{\s\s} - m_j \bT{p,j}  \bT{l,i}' c_j c_i g_{\s\t}
!latex    +     m_i \bT{l,i}  \bT{p,j}' s_i s_j g_{\s\t} +     \bT{p,j}' \bT{l,i}' s_j c_i g_{\t\t}) / \sqrt g^2\nonumber   
!latex \ee 
!latex \be
!latex \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial \Ate{j,p}} \frac{\partial}{\partial \Azo{i,l}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!latex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!l tex ( - n_j \bT{p,j}  \sj \, {\bf e}_\s + \bT{p,j}' \cj \, {\bf e}_\z ) \cdot
!l tex ( + m_i \bT{l,i}  \ci \, {\bf e}_\s - \bT{l,i}' \si \, {\bf e}_\t ) / \sqrt g^2 \nonumber \\
!l tex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \;
!latex  ( - n_j m_i \bT{p,j}  \bT{l,i}  s_j c_i g_{\s\s} + n_j \bT{p,j}  \bT{l,i}' s_j s_i g_{\s\t}
!latex    +     m_i \bT{l,i}  \bT{p,j}' c_i c_j g_{\s\z} -     \bT{p,j}' \bT{l,i}' c_j s_i g_{\t\z}) / \sqrt g^2\nonumber \\
!latex \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial \Ato{j,p}} \frac{\partial}{\partial \Azo{i,l}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!latex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!l tex ( + n_j \bT{p,j}  \cj \, {\bf e}_\s + \bT{p,j}' \sj \, {\bf e}_\z ) \cdot
!l tex ( + m_i \bT{l,i}  \ci \, {\bf e}_\s - \bT{l,i}' \si \, {\bf e}_\t ) / \sqrt g^2 \nonumber \\
!l tex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \;
!latex  ( + n_j m_i \bT{p,j}  \bT{l,i}  c_j c_i g_{\s\s} - n_j \bT{p,j}  \bT{l,i}' c_j s_i g_{\s\t}
!latex    +     m_i \bT{l,i}  \bT{p,j}' c_i s_j g_{\s\z} -     \bT{p,j}' \bT{l,i}' s_j s_i g_{\t\z}) / \sqrt g^2\nonumber \\
!latex \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial \Aze{j,p}} \frac{\partial}{\partial \Azo{i,l}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!latex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!l tex ( - m_j \bT{p,j}  \sj \, {\bf e}_\s - \bT{p,j}' \cj \, {\bf e}_\t ) \cdot
!l tex ( + m_i \bT{l,i}  \ci \, {\bf e}_\s - \bT{l,i}' \si \, {\bf e}_\t ) / \sqrt g^2 \nonumber \\
!l tex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \;
!latex  ( - m_j m_i \bT{p,j}  \bT{l,i}  s_j c_i g_{\s\s} + m_j \bT{p,j}  \bT{l,i}' s_j s_i g_{\s\t}
!latex    -     m_i \bT{l,i}  \bT{p,j}' c_i c_j g_{\s\t} +     \bT{p,j}' \bT{l,i}' c_j s_i g_{\t\t}) / \sqrt g^2\nonumber \\
!latex \!\!\!\!\!\!\!\!\! \frac{\partial}{\partial \Azo{j,p}} \frac{\partial}{\partial \Azo{i,l}} \! \int \!\! dv \; {\bf B} \! \cdot \! {\bf B}
!latex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \; 
!l tex ( + m_j \bT{p,j}  \cj \, {\bf e}_\s - \bT{p,j}' \sj \, {\bf e}_\t ) \cdot
!l tex ( + m_i \bT{l,i}  \ci \, {\bf e}_\s - \bT{l,i}' \si \, {\bf e}_\t ) / \sqrt g^2 \nonumber \\
!l tex \!\!&\!\!=\!\!&\!\!2 \int \!\! dv \;
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
!latex + \bT{l,i} \ci \nabla \t \cdot \bT{p,j}' \cj \, {\bf e}_\z + \bT{p,j} \cj \nabla \t \cdot \bT{l,i}' \ci \, {\bf e}_\z \right] / \sqrt g \\
!latex \frac{\partial}{\partial \Ato{j,p}} \frac{\partial}{\partial \Ate{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &   \int \!\! dv \;  \left[
!latex + \bT{l,i} \ci \nabla \t \cdot \bT{p,j}' \sj \, {\bf e}_\z + \bT{p,j} \sj \nabla \t \cdot \bT{l,i}' \ci \, {\bf e}_\z \right] / \sqrt g \\
!latex \frac{\partial}{\partial \Aze{j,p}} \frac{\partial}{\partial \Ate{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &   \int \!\! dv \;  \left[
!latex - \bT{l,i} \ci \nabla \t \cdot \bT{p,j}' \cj \, {\bf e}_\t + \bT{p,j} \cj \nabla \z \cdot \bT{l,i}' \ci \, {\bf e}_\z \right] / \sqrt g \\
!latex \frac{\partial}{\partial \Azo{j,p}} \frac{\partial}{\partial \Ate{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &   \int \!\! dv \;  \left[
!latex - \bT{l,i} \ci \nabla \t \cdot \bT{p,j}' \sj \, {\bf e}_\t + \bT{p,j} \sj \nabla \z \cdot \bT{l,i}' \ci \, {\bf e}_\z \right] / \sqrt g \\
!latex \nonumber \\
!latex \frac{\partial}{\partial \Ate{j,p}} \frac{\partial}{\partial \Ato{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &   \int \!\! dv \;  \left[
!latex + \bT{l,i} \si \nabla \t \cdot \bT{p,j}' \cj \, {\bf e}_\z + \bT{p,j} \cj \nabla \t \cdot \bT{l,i}' \si \, {\bf e}_\z \right] / \sqrt g \\
!latex \frac{\partial}{\partial \Ato{j,p}} \frac{\partial}{\partial \Ato{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &   \int \!\! dv \;  \left[
!latex + \bT{l,i} \si \nabla \t \cdot \bT{p,j}' \sj \, {\bf e}_\z + \bT{p,j} \sj \nabla \t \cdot \bT{l,i}' \si \, {\bf e}_\z \right] / \sqrt g \\
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
!latex - \bT{l,i} \ci \nabla \z \cdot \bT{p,j}' \cj \, {\bf e}_\t - \bT{p,j} \cj \nabla \z \cdot \bT{l,i}' \ci \, {\bf e}_\t \right] / \sqrt g \\
!latex \frac{\partial}{\partial \Azo{j,p}} \frac{\partial}{\partial \Aze{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &   \int \!\! dv \;  \left[
!latex - \bT{l,i} \ci \nabla \z \cdot \bT{p,j}' \sj \, {\bf e}_\t - \bT{p,j} \sj \nabla \z \cdot \bT{l,i}' \ci \, {\bf e}_\t \right] / \sqrt g \\
!latex \nonumber \\
!latex \frac{\partial}{\partial \Ate{j,p}} \frac{\partial}{\partial \Azo{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &   \int \!\! dv \;  \left[
!latex + \bT{l,i} \si \nabla \z \cdot \bT{p,j}' \cj \, {\bf e}_\z - \bT{p,j} \cj \nabla \t \cdot \bT{l,i}' \si \, {\bf e}_\t \right] / \sqrt g \\
!latex \frac{\partial}{\partial \Ato{j,p}} \frac{\partial}{\partial \Azo{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &   \int \!\! dv \;  \left[
!latex + \bT{l,i} \si \nabla \z \cdot \bT{p,j}' \sj \, {\bf e}_\z - \bT{p,j} \sj \nabla \t \cdot \bT{l,i}' \si \, {\bf e}_\t \right] / \sqrt g \\
!latex \frac{\partial}{\partial \Aze{j,p}} \frac{\partial}{\partial \Azo{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &   \int \!\! dv \;  \left[
!latex - \bT{l,i} \si \nabla \z \cdot \bT{p,j}' \cj \, {\bf e}_\t - \bT{p,j} \cj \nabla \z \cdot \bT{l,i}' \si \, {\bf e}_\t \right] / \sqrt g \\
!latex \frac{\partial}{\partial \Azo{j,p}} \frac{\partial}{\partial \Azo{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = &   \int \!\! dv \;  \left[
!latex - \bT{l,i} \si \nabla \z \cdot \bT{p,j}' \sj \, {\bf e}_\t - \bT{p,j} \sj \nabla \z \cdot \bT{l,i}' \si \, {\bf e}_\t \right] / \sqrt g
!latex \ee
!latex \end{enumerate}

!latex \subsection{derivatives of kinetic  energy integrals}

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

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine matrix( lvol, mn, lrad )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, one

  use numerical, only : sqrtmachprec, small

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

  LOGICAL              :: Lpause

  INTEGER              :: NN, ii, jj, ll, kk, pp, mi, ni, mj, nj, mimj, minj, nimj, ninj, mjmi, mjni, njmi, njni, id, jd

  REAL                 :: Wtete(1:mn,0:lrad,1:mn,0:lrad), Wteto(1:mn,0:lrad,1:mn,0:lrad), Wtote(1:mn,0:lrad,1:mn,0:lrad), Wtoto(1:mn,0:lrad,1:mn,0:lrad)
  REAL                 :: Wteze(1:mn,0:lrad,1:mn,0:lrad), Wtezo(1:mn,0:lrad,1:mn,0:lrad), Wtoze(1:mn,0:lrad,1:mn,0:lrad), Wtozo(1:mn,0:lrad,1:mn,0:lrad)
  REAL                 :: Wzete(1:mn,0:lrad,1:mn,0:lrad), Wzeto(1:mn,0:lrad,1:mn,0:lrad), Wzote(1:mn,0:lrad,1:mn,0:lrad), Wzoto(1:mn,0:lrad,1:mn,0:lrad)
  REAL                 :: Wzeze(1:mn,0:lrad,1:mn,0:lrad), Wzezo(1:mn,0:lrad,1:mn,0:lrad), Wzoze(1:mn,0:lrad,1:mn,0:lrad), Wzozo(1:mn,0:lrad,1:mn,0:lrad)

  REAL                 :: Htete(1:mn,0:lrad,1:mn,0:lrad), Hteto(1:mn,0:lrad,1:mn,0:lrad), Htote(1:mn,0:lrad,1:mn,0:lrad), Htoto(1:mn,0:lrad,1:mn,0:lrad)
  REAL                 :: Hteze(1:mn,0:lrad,1:mn,0:lrad), Htezo(1:mn,0:lrad,1:mn,0:lrad), Htoze(1:mn,0:lrad,1:mn,0:lrad), Htozo(1:mn,0:lrad,1:mn,0:lrad)
  REAL                 :: Hzete(1:mn,0:lrad,1:mn,0:lrad), Hzeto(1:mn,0:lrad,1:mn,0:lrad), Hzote(1:mn,0:lrad,1:mn,0:lrad), Hzoto(1:mn,0:lrad,1:mn,0:lrad)
  REAL                 :: Hzeze(1:mn,0:lrad,1:mn,0:lrad), Hzezo(1:mn,0:lrad,1:mn,0:lrad), Hzoze(1:mn,0:lrad,1:mn,0:lrad), Hzozo(1:mn,0:lrad,1:mn,0:lrad)
  
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
  
  dMA(0:NN,0:NN) = zero ! initialize summation; 24 Jan 13;
  dMD(0:NN,0:NN) = zero
  dMB(0:NN,1: 2) = zero
  dMG(0:NN     ) = zero
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  Wtete(1:mn,0:lrad,1:mn,0:lrad) = zero ! these "zeros" are probably not required; 10 Mar 13;
  Wtote(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Wzete(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Wzote(1:mn,0:lrad,1:mn,0:lrad) = zero
  
  Wteto(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Wtoto(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Wzeto(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Wzoto(1:mn,0:lrad,1:mn,0:lrad) = zero
  
  Wteze(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Wtoze(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Wzeze(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Wzoze(1:mn,0:lrad,1:mn,0:lrad) = zero
  
  Wtezo(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Wtozo(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Wzezo(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Wzozo(1:mn,0:lrad,1:mn,0:lrad) = zero
  
  Htete(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Htote(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Hzete(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Hzote(1:mn,0:lrad,1:mn,0:lrad) = zero
  
  Hteto(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Htoto(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Hzeto(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Hzoto(1:mn,0:lrad,1:mn,0:lrad) = zero
  
  Hteze(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Htoze(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Hzeze(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Hzoze(1:mn,0:lrad,1:mn,0:lrad) = zero
  
  Htezo(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Htozo(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Hzezo(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Hzozo(1:mn,0:lrad,1:mn,0:lrad) = zero
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! construct derivatives of energy and helicity; 08 Feb 16;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( YESstellsym ) then
   
   do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
    
    do jj = 1, mn ; mj = im(jj) ; nj = in(jj) ; mimj = mi * mj ; minj = mi * nj ; nimj = ni * mj ; ninj = ni * nj
     
     do ll = 0, lrad
      
      do pp = 0, lrad
       
       Wtete(jj,pp,ii,ll) = + 2 * ninj * TTssss(ll,pp,ii,jj) - 2 * ni      * TDszsc(ll,pp,ii,jj) - 2      * nj * TDszsc(pp,ll,jj,ii) + 2 * DDzzcc(ll,pp,ii,jj)
       Wzete(jj,pp,ii,ll) = + 2 * nimj * TTssss(ll,pp,ii,jj) + 2 * ni      * TDstsc(ll,pp,ii,jj) - 2      * mj * TDszsc(pp,ll,jj,ii) - 2 * DDtzcc(pp,ll,jj,ii)
       Wteze(jj,pp,ii,ll) = + 2 * minj * TTssss(ll,pp,ii,jj) + 2      * nj * TDstsc(pp,ll,jj,ii) - 2 * mi      * TDszsc(ll,pp,ii,jj) - 2 * DDtzcc(ll,pp,ii,jj)
       Wzeze(jj,pp,ii,ll) = + 2 * mimj * TTssss(ll,pp,ii,jj) + 2 * mi      * TDstsc(ll,pp,ii,jj) + 2      * mj * TDstsc(pp,ll,jj,ii) + 2 * DDttcc(ll,pp,ii,jj)
       
!     !Htete(jj,pp,ii,ll) =   zero
       Hzete(jj,pp,ii,ll) = - DToocc(pp,ll,jj,ii) + DToocc(ll,pp,ii,jj)
       Hteze(jj,pp,ii,ll) = + DToocc(pp,ll,jj,ii) - DToocc(ll,pp,ii,jj)
!     !Hzeze(jj,pp,ii,ll) =   zero  
      
      enddo ! end of do pp; 08 Feb 16;

     enddo ! end of do ll; 08 Feb 16;

    enddo ! end of do jj; 08 Feb 16;

   enddo ! end of do ii; 08 Feb 16;
   
  else ! NOTstellsym ; 08 Feb 16;
   
   do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
    
    do jj = 1, mn ; mj = im(jj) ; nj = in(jj) ; mjmi = mi * mj ; mjni = ni * mj ; njmi = mi * nj ; njni = ni * nj
     
     do ll = 0, lrad
      
      do pp = 0, lrad
     
      Wtete(jj,pp,ii,ll) = 2 * ( + njni * TTssss(pp,ll,jj,ii) - nj * TDszsc(pp,ll,jj,ii) - ni * TDszsc(ll,pp,ii,jj) + DDzzcc(pp,ll,jj,ii) )
      Wtote(jj,pp,ii,ll) = 2 * ( - njni * TTsscs(pp,ll,jj,ii) + nj * TDszcc(pp,ll,jj,ii) - ni * TDszss(ll,pp,ii,jj) + DDzzsc(pp,ll,jj,ii) )
      Wzete(jj,pp,ii,ll) = 2 * ( + mjni * TTssss(pp,ll,jj,ii) - mj * TDszsc(pp,ll,jj,ii) + ni * TDstsc(ll,pp,ii,jj) - DDtzcc(pp,ll,jj,ii) )
      Wzote(jj,pp,ii,ll) = 2 * ( - mjni * TTsscs(pp,ll,jj,ii) + mj * TDszcc(pp,ll,jj,ii) + ni * TDstss(ll,pp,ii,jj) - DDtzsc(pp,ll,jj,ii) )
     
      Wteto(jj,pp,ii,ll) = 2 * ( - njni * TTsssc(pp,ll,jj,ii) - nj * TDszss(pp,ll,jj,ii) + ni * TDszcc(ll,pp,ii,jj) + DDzzcs(pp,ll,jj,ii) )
      Wtoto(jj,pp,ii,ll) = 2 * ( + njni * TTsscc(pp,ll,jj,ii) + nj * TDszcs(pp,ll,jj,ii) + ni * TDszcs(ll,pp,ii,jj) + DDzzss(pp,ll,jj,ii) )
!     Wzeto(jj,pp,ii,ll) = 2 * ( - mjni * TTsssc(pp,ll,jj,ii) - mj * TDszcc(pp,ll,jj,ii) - ni * TDstcc(ll,pp,ii,jj) - DDtzcs(pp,ll,jj,ii) ) ! 19 Jul 16;
      Wzeto(jj,pp,ii,ll) = 2 * ( - mjni * TTsssc(pp,ll,jj,ii) - mj * TDszss(pp,ll,jj,ii) - ni * TDstcc(ll,pp,ii,jj) - DDtzcs(pp,ll,jj,ii) ) ! 19 Jul 16;
      Wzoto(jj,pp,ii,ll) = 2 * ( + mjni * TTsscc(pp,ll,jj,ii) + mj * TDszcs(pp,ll,jj,ii) - ni * TDstcs(ll,pp,ii,jj) - DDtzss(pp,ll,jj,ii) )

      Wteze(jj,pp,ii,ll) = 2 * ( + njmi * TTssss(pp,ll,jj,ii) + nj * TDstsc(pp,ll,jj,ii) - mi * TDszsc(ll,pp,ii,jj) - DDtzcc(pp,ll,jj,ii) )
      Wtoze(jj,pp,ii,ll) = 2 * ( - njmi * TTsscs(pp,ll,jj,ii) - nj * TDstcc(pp,ll,jj,ii) - mi * TDszss(ll,pp,ii,jj) - DDtzsc(pp,ll,jj,ii) )
      Wzeze(jj,pp,ii,ll) = 2 * ( + mjmi * TTssss(pp,ll,jj,ii) + mj * TDstsc(pp,ll,jj,ii) + mi * TDstsc(ll,pp,ii,jj) + DDttcc(pp,ll,jj,ii) )
      Wzoze(jj,pp,ii,ll) = 2 * ( - mjmi * TTsscs(pp,ll,jj,ii) - mj * TDstcc(pp,ll,jj,ii) + mi * TDstss(ll,pp,ii,jj) + DDttsc(pp,ll,jj,ii) )

      Wtezo(jj,pp,ii,ll) = 2 * ( - njmi * TTsssc(pp,ll,jj,ii) + nj * TDstss(pp,ll,jj,ii) + mi * TDszcc(ll,pp,ii,jj) - DDtzcs(pp,ll,jj,ii) )
      Wtozo(jj,pp,ii,ll) = 2 * ( + njmi * TTsscc(pp,ll,jj,ii) - nj * TDstcs(pp,ll,jj,ii) + mi * TDszcs(ll,pp,ii,jj) - DDtzss(pp,ll,jj,ii) )
      Wzezo(jj,pp,ii,ll) = 2 * ( - mjmi * TTsssc(pp,ll,jj,ii) + mj * TDstss(pp,ll,jj,ii) - mi * TDstcc(ll,pp,ii,jj) + DDttcs(pp,ll,jj,ii) )
      Wzozo(jj,pp,ii,ll) = 2 * ( + mjmi * TTsscc(pp,ll,jj,ii) - mj * TDstcs(pp,ll,jj,ii) - mi * TDstcs(ll,pp,ii,jj) + DDttss(pp,ll,jj,ii) )

      Htete(jj,pp,ii,ll) =   zero
      Htote(jj,pp,ii,ll) =   zero
      Hzete(jj,pp,ii,ll) = - DToocc(pp,ll,jj,ii) + DToocc(ll,pp,ii,jj)
      Hzote(jj,pp,ii,ll) = - DToosc(pp,ll,jj,ii) + DToocs(ll,pp,ii,jj)

      Hteto(jj,pp,ii,ll) =   zero
      Htoto(jj,pp,ii,ll) =   zero
      Hzeto(jj,pp,ii,ll) = - DToocs(pp,ll,jj,ii) + DToosc(ll,pp,ii,jj)
      Hzoto(jj,pp,ii,ll) = - DTooss(pp,ll,jj,ii) + DTooss(ll,pp,ii,jj)  

      Hteze(jj,pp,ii,ll) = + DToocc(pp,ll,jj,ii) - DToocc(ll,pp,ii,jj)
      Htoze(jj,pp,ii,ll) = + DToosc(pp,ll,jj,ii) - DToocs(ll,pp,ii,jj)
      Hzeze(jj,pp,ii,ll) =   zero  
      Hzoze(jj,pp,ii,ll) =   zero

      Htezo(jj,pp,ii,ll) = + DToocs(pp,ll,jj,ii) - DToosc(ll,pp,ii,jj)
      Htozo(jj,pp,ii,ll) = + DTooss(pp,ll,jj,ii) - DTooss(ll,pp,ii,jj)
      Hzezo(jj,pp,ii,ll) =   zero
      Hzozo(jj,pp,ii,ll) =   zero
      
      enddo ! end of do pp;
      
     enddo ! end of do jj;
     
    enddo ! end of do ll;
    
   enddo ! end of do ii;
   
  endif ! end of if( YESstellsym ) ; 08 Feb 16;
      
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! construct matrix elements;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( YESstellsym ) then
   
   do ii = 1, mn ; mi = im(ii) ; ni = in(ii) ! Fourier harmonics; 17 Jan 13;
    
    if( Lcoordinatesingularity .and. ii.eq.1 ) then ; kk = 1
    else                                            ; kk = 0
    endif
    
    do ll = 0, lrad ! Chebyshev polynomial; 08 Feb 16;
     
     do jj = 1, mn
      do pp = 0, lrad   ; id = Ate(lvol,0,ii)%i(ll) ; jd = Ate(lvol,0,jj)%i(pp) ; dMA(id,jd) = + Wtete(jj,pp,ii,ll)
       ;                ;                           ;                           ; dMD(id,jd) = - Htete(jj,pp,ii,ll)
       ;                ; id = Ate(lvol,0,ii)%i(ll) ; jd = Aze(lvol,0,jj)%i(pp) ; dMA(id,jd) = + Wzete(jj,pp,ii,ll)
       ;                ;                                                       ; dMD(id,jd) = - Hzete(jj,pp,ii,ll)
       ;                ; id = Aze(lvol,0,ii)%i(ll) ; jd = Ate(lvol,0,jj)%i(pp) ; dMA(id,jd) = + Wteze(jj,pp,ii,ll)
       ;                ;                                                       ; dMD(id,jd) = - Hteze(jj,pp,ii,ll)
       ;                ; id = Aze(lvol,0,ii)%i(ll) ; jd = Aze(lvol,0,jj)%i(pp) ; dMA(id,jd) = + Wzeze(jj,pp,ii,ll)
       ;                ;                                                       ; dMD(id,jd) = - Hzeze(jj,pp,ii,ll)
      enddo ! end of do pp; 08 Feb 16;
     enddo ! end of do jj; 08 Feb 16;
     
     ;                  ; id = Ate(lvol,0,ii)%i(ll) ; jd = Lma(lvol,  ii)       ; dMA(id,jd) = +      TT(ll, 0,0)
     ;                  ; id = Aze(lvol,0,ii)%i(ll) ; jd = Lmb(lvol,  ii)       ; dMA(id,jd) = +      TT(ll,kk,0) ! check coordinate singularity; 18 May 16;
     if( ii.gt.1 ) then ; id = Ate(lvol,0,ii)%i(ll) ; jd = Lme(lvol,  ii)       ; dMA(id,jd) = - ni * TT(ll, 1,0)
      ;                 ; id = Aze(lvol,0,ii)%i(ll) ; jd = Lme(lvol,  ii)       ; dMA(id,jd) = - mi * TT(ll, 1,0)
     else               ; id = Ate(lvol,0,ii)%i(ll) ; jd = Lmg(lvol,  ii)       ; dMA(id,jd) = +      TT(ll, 1,0)
      ;                 ; id = Aze(lvol,0,ii)%i(ll) ; jd = Lmh(lvol,  ii)       ; dMA(id,jd) = +      TT(ll, 1,0)
     endif
     
    enddo ! end of do ll; 08 Feb 16;
    
    do pp = 0, lrad     ; id = Lma(lvol,  ii)       ; jd = Ate(lvol,0,ii)%i(pp) ; dMA(id,jd) = +      TT(pp, 0,0)
     ;                  ; id = Lmb(lvol,  ii)       ; jd = Aze(lvol,0,ii)%i(pp) ; dMA(id,jd) = +      TT(pp,kk,0) ! check coordinate singularity; 18 May 16;
     if( ii.gt.1 ) then ; id = Lme(lvol,  ii)       ; jd = Ate(lvol,0,ii)%i(pp) ; dMA(id,jd) = - ni * TT(pp, 1,0)
      ;                 ; id = Lme(lvol,  ii)       ; jd = Aze(lvol,0,ii)%i(pp) ; dMA(id,jd) = - mi * TT(pp, 1,0)
     else               ; id = Lmg(lvol,  ii)       ; jd = Ate(lvol,0,ii)%i(pp) ; dMA(id,jd) = +      TT(pp, 1,0)
      ;                 ; id = Lmh(lvol,  ii)       ; jd = Aze(lvol,0,ii)%i(pp) ; dMA(id,jd) = +      TT(pp, 1,0)
     endif
    enddo ! end of do pp; 08 Feb 16;
    
!   ;if( ii.gt.1 ) then ; id = Lme(lvol,  ii)       ;                           ; dMG(id   ) = - ( iVns(ii) + iBns(ii) )
    ;if( ii.gt.1 ) then ; id = Lme(lvol,  ii)       ;                           ; dMG(id   ) = - ( iVns(ii) - iBns(ii) ) ! 14 Apr 17;
    ;else               ; id = Lmg(lvol,  ii)       ;                           ; dMB(id, 1) = -       one
    ;                   ; id = Lmh(lvol,  ii)       ;                           ; dMB(id, 2) = -       one
    ;endif
    
   enddo ! end of do ii; 08 Feb 16;
   
  else ! NOTstellsym;
   
   do ii = 1, mn ; mi = im(ii) ; ni = in(ii) ! Fourier harmonics; 17 Jan 13;
    
    if( Lcoordinatesingularity .and. ii.eq.1 ) then ; kk = 1
    else                                            ; kk = 0
    endif
    
    do ll = 0, lrad ! Chebyshev polynomial; 08 Feb 16;
     
     do jj = 1, mn
      do pp = 0, lrad   ; id = Ate(lvol,0,ii)%i(ll) ; jd = Ate(lvol,0,jj)%i(pp) ; dMA(id,jd) = + Wtete(jj,pp,ii,ll)
       ;                ;                           ;                           ; dMD(id,jd) = - Htete(jj,pp,ii,ll)
       ;                ; id = Ate(lvol,0,ii)%i(ll) ; jd = Ato(lvol,0,jj)%i(pp) ; dMA(id,jd) = + Wtote(jj,pp,ii,ll)
       ;                ;                           ;                           ; dMD(id,jd) = - Htote(jj,pp,ii,ll)
       ;                ; id = Ate(lvol,0,ii)%i(ll) ; jd = Aze(lvol,0,jj)%i(pp) ; dMA(id,jd) = + Wzete(jj,pp,ii,ll)
       ;                ;                           ;                           ; dMD(id,jd) = - Hzete(jj,pp,ii,ll)
       ;                ; id = Ate(lvol,0,ii)%i(ll) ; jd = Azo(lvol,0,jj)%i(pp) ; dMA(id,jd) = + Wzote(jj,pp,ii,ll)
       ;                ;                           ;                           ; dMD(id,jd) = - Hzote(jj,pp,ii,ll)

       ;                ; id = Ato(lvol,0,ii)%i(ll) ; jd = Ate(lvol,0,jj)%i(pp) ; dMA(id,jd) = + Wteto(jj,pp,ii,ll)
       ;                ;                           ;                           ; dMD(id,jd) = - Hteto(jj,pp,ii,ll)
       ;                ; id = Ato(lvol,0,ii)%i(ll) ; jd = Ato(lvol,0,jj)%i(pp) ; dMA(id,jd) = + Wtoto(jj,pp,ii,ll)
       ;                ;                           ;                           ; dMD(id,jd) = - Htoto(jj,pp,ii,ll)
       ;                ; id = Ato(lvol,0,ii)%i(ll) ; jd = Aze(lvol,0,jj)%i(pp) ; dMA(id,jd) = + Wzeto(jj,pp,ii,ll)
       ;                ;                           ;                           ; dMD(id,jd) = - Hzeto(jj,pp,ii,ll)
       ;                ; id = Ato(lvol,0,ii)%i(ll) ; jd = Azo(lvol,0,jj)%i(pp) ; dMA(id,jd) = + Wzoto(jj,pp,ii,ll)
       ;                ;                           ;                           ; dMD(id,jd) = - Hzoto(jj,pp,ii,ll)

       ;                ; id = Aze(lvol,0,ii)%i(ll) ; jd = Ate(lvol,0,jj)%i(pp) ; dMA(id,jd) = + Wteze(jj,pp,ii,ll)
       ;                ;                           ;                           ; dMD(id,jd) = - Hteze(jj,pp,ii,ll)
       ;                ; id = Aze(lvol,0,ii)%i(ll) ; jd = Ato(lvol,0,jj)%i(pp) ; dMA(id,jd) = + Wtoze(jj,pp,ii,ll)
       ;                ;                           ;                           ; dMD(id,jd) = - Htoze(jj,pp,ii,ll)
       ;                ; id = Aze(lvol,0,ii)%i(ll) ; jd = Aze(lvol,0,jj)%i(pp) ; dMA(id,jd) = + Wzeze(jj,pp,ii,ll)
       ;                ;                           ;                           ; dMD(id,jd) = - Hzeze(jj,pp,ii,ll)
       ;                ; id = Aze(lvol,0,ii)%i(ll) ; jd = Azo(lvol,0,jj)%i(pp) ; dMA(id,jd) = + Wzoze(jj,pp,ii,ll)
       ;                ;                           ;                           ; dMD(id,jd) = - Hzoze(jj,pp,ii,ll)

       ;                ; id = Azo(lvol,0,ii)%i(ll) ; jd = Ate(lvol,0,jj)%i(pp) ; dMA(id,jd) = + Wtezo(jj,pp,ii,ll)
       ;                ;                           ;                           ; dMD(id,jd) = - Htezo(jj,pp,ii,ll)
       ;                ; id = Azo(lvol,0,ii)%i(ll) ; jd = Ato(lvol,0,jj)%i(pp) ; dMA(id,jd) = + Wtozo(jj,pp,ii,ll)
       ;                ;                           ;                           ; dMD(id,jd) = - Htozo(jj,pp,ii,ll)
       ;                ; id = Azo(lvol,0,ii)%i(ll) ; jd = Aze(lvol,0,jj)%i(pp) ; dMA(id,jd) = + Wzezo(jj,pp,ii,ll)
       ;                ;                           ;                           ; dMD(id,jd) = - Hzezo(jj,pp,ii,ll)
       ;                ; id = Azo(lvol,0,ii)%i(ll) ; jd = Azo(lvol,0,jj)%i(pp) ; dMA(id,jd) = + Wzozo(jj,pp,ii,ll)
       ;                ;                           ;                           ; dMD(id,jd) = - Hzozo(jj,pp,ii,ll)
      enddo ! end of do pp; 08 Feb 16;
     enddo ! end of do jj; 08 Feb 16;
     
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
     
    enddo ! end of do ll; 08 Feb 16;
    
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
    enddo ! end of do pp; 08 Feb 16;
    
!   ;if( ii.gt.1 ) then ; id = Lme(lvol,ii)         ;                           ; dMG(id   ) = - ( iVns(ii) + iBns(ii) )
!   ;                   ; id = Lmf(lvol,ii)         ;                           ; dMG(id   ) = - ( iVnc(ii) + iBnc(ii) )
    ;if( ii.gt.1 ) then ; id = Lme(lvol,ii)         ;                           ; dMG(id   ) = - ( iVns(ii) - iBns(ii) ) ! 14 Apr 17;
    ;                   ; id = Lmf(lvol,ii)         ;                           ; dMG(id   ) = - ( iVnc(ii) - iBnc(ii) ) ! 14 Apr 17;
    ;else               ; id = Lmg(lvol,ii)         ;                           ; dMB(id, 1) = -       one
    ;                   ; id = Lmh(lvol,ii)         ;                           ; dMB(id, 2) = -       one
    ;endif
    
   enddo ! end of do ii; 08 Feb 16;
   
  endif ! matches if( YESstellsym ) ; 08 Feb 16;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  
  if( Wmatrix ) then

   do ii = 1, NN

    Lpause = .false.

    do jj = 1, NN
     if( abs(dMA(ii,jj)-dMA(jj,ii)) .gt. small*abs(dMA(ii,jj)+dMA(jj,ii)) ) then
      write(ounit,1000) myid, dMA(ii,jj), dMA(jj,ii), dMA(ii,jj)-dMA(jj,ii)
      Lpause = .true.
     endif

     if( abs(dMD(ii,jj)-dMD(jj,ii)) .gt. small*abs(dMD(ii,jj)+dMD(jj,ii)) ) then
      write(ounit,1001) myid, dMD(ii,jj), dMD(jj,ii), dMD(ii,jj)-dMD(jj,ii)
      Lpause = .true.
     endif

    enddo ! end of do jj; 18 May 16;

    if( ncpu.eq.1 .and. Lpause ) pause

   enddo ! end of do ii; 18 May 16;

  endif ! end of if( Wmatrix ) ; 18 May 16;
   
1000 format("matrix : " 10x " : myid="i3" : dMA(ii,jj)="es23.15", dMA(jj,ii)="es23.15", dMA(ii,jj)-dMA(jj,ii)="es13.5" ;")
1001 format("matrix : " 10x " : myid="i3" : dMD(ii,jj)="es23.15", dMD(jj,ii)="es23.15", dMD(ii,jj)-dMD(jj,ii)="es13.5" ;")

#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(matrix)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine matrix

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

