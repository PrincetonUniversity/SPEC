!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (build matrices) ! Constructs matrices that represent the Beltrami linear system.

!latex \briefly{Constructs matrices that represent the Beltrami linear system.}

!latex \calledby{\link{dforce}}
!      \calls{\link{}}

!latex \tableofcontents

!latex \subsection{gauge conditions} \begin{enumerate}

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

!latex \end{enumerate} \subsection{boundary conditions} \begin{enumerate}

!latex \item The magnetic field is 
!latex       $\sqrt g \, {\bf B} = (\partial_\t A_\z - \partial_\z A_\t)\;{\bf e}_\s - \partial_\s A_\z \;{\bf e}_\t + \partial_\s A_\t \;{\bf e}_\z$.

!latex \item In the annular volumes, the condition that the field is tangential to the inner interface gives $\partial_\t A_\z - \partial_\z A_\t = 0$.
!latex       With the above gauge condition on $A_\t$ given in \Eqn{Atgauge}, this gives $\partial_\t A_\z=0$, which with \Eqn{Azgauge} gives
!latex       \be A_\z(-1,\t,\z)=0.
!latex       \ee

!latex \item The condition at the outer interface is that the field is $\sqrt g \, {\bf B}\cdot\nabla s = \partial_\t A_\z - \partial_\z A_\t = b$,
!latex       where $b$ is supplied by the user.
!latex       For each of the plasma regions, $b=0$.
!latex       For the vacuum region, generally $b\ne0$.

!latex \end{enumerate} \subsection{enclosed fluxes} \begin{enumerate}

!latex \item In the plasma regions, the enclosed fluxes must be constrained.

!latex \item The toroidal and poloidal fluxes enclosed in each volume are determined using
!latex       \be \int_S {\bf B}\cdot{\bf ds}=\int_{\partial S}{\bf A}\cdot {\bf dl}.
!latex       \ee

!l tex \end{enumerate} \subsection{enclosed currents} \begin{enumerate}
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

!latex \end{enumerate} \subsection{Fourier-Chebyshev representation} \begin{enumerate}

!latex \item The components of the vector potential are
!latex       \be A_\t(\s,\t,\z) &=& \sum_{i,l} \Ate{i,l} \; T_l(s) \cos\alpha_i + \sum_{i,l} \Ato{i,l} \; T_l(s) \sin\alpha_i, \\
!latex           A_\z(\s,\t,\z) &=& \sum_{i,l} \Aze{i,l} \; T_l(s) \cos\alpha_i + \sum_{i,l} \Azo{i,l} \; T_l(s) \sin\alpha_i,
!latex       \ee
!latex       where $\alpha_j \equiv m_j\t-n_j\z$.

!latex \end{enumerate} \subsection{constrained energy functional} \begin{enumerate}

!latex \item The constrained energy functional in each volume is
!latex \be \begin{array}{cclcclcccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc}                                                   
!latex     {\cal F}_l & \equiv & \ds                         \int {\bf B} \cdot {\bf B} \, dv                                                      
!latex                   -           \frac{\mu}{2}    \left[ \int {\bf A} \cdot {\bf B} \, dv - K_{l} \right]                                    \\
!latex                 &  + &   & \ds \sum_i & a_i       & \ds \left[  \sum_l \Ate{i,l} T_l(-1) - 0             \right]                          \\
!latex                 &  + &   & \ds \sum_i & b_i       & \ds \left[  \sum_l \Aze{i,l} T_l(-1) - 0             \right]                          \\
!latex                 &  + &   & \ds \sum_i & c_i       & \ds \left[  \sum_l \Ato{i,l} T_l(-1) - 0             \right]                          \\
!latex                 &  + &   & \ds \sum_i & d_i       & \ds \left[  \sum_l \Azo{i,l} T_l(-1) - 0             \right]                          \\
!latex                 &  + &   & \ds \sum_i & e_i       & \ds \left[  \sum_l \left( - m_i \Aze{i,l} - n_i \Ate{i,l} \right) T_l(+1) - b_{s,i} \right] \\
!latex                 &  + &   & \ds \sum_i & f_i       & \ds \left[  \sum_l \left( + m_i \Azo{i,l} + n_i \Ato{i,l} \right) T_l(+1) - b_{c,i} \right] \\
!latex                 &  + &   & \ds        & \alpha    & \ds \left[  \sum_l \Ate{1,l} T_l(+1) - \Delta \psi_t \right]                          \\
!latex                 &  + &   & \ds        & \beta     & \ds \left[  \sum_l \Aze{1,l} T_l(+1) - \Delta \psi_p \right]                             
!l tex                 &  + &   & \ds        & \gamma    & \ds \left[ \int_{0}^{2\pi} {\bf B} \cdot {\bf e}_\t \, d\t - I \right]                \\
!l tex                 &  + &   & \ds        & \delta    & \ds \left[ \int_{0}^{2\pi} {\bf B} \cdot {\bf e}_\z \, d\z - G \right],
!latex     \end{array}
!latex \ee   
!latex       where
!latex       \bi
!latex       \item[  i.] $a_i$, $b_i$, $c_i$ and $d_i$ are Lagrange multipliers used to enforce the combined gauge and interface boundary condition
!latex                 on the inner interface,
!latex       \item[ ii.] $e_i$ and $f_i$               are Lagrange multipliers used to enforce the interface boundary condition on the outer interface, 
!latex       namely \mbox{$\sqrt g {\bf B}\cdot\nabla s = b$}; and
!latex       \item[iii.] $\alpha$ and $\beta$          are Lagrange multipliers used to enforce the constraints on the enclosed fluxes.
!l tex       \item[ iv.] $\gamma$ and $\delta$         are Lagrange multipliers used to enforce the condition on the enclosed currents.
!latex       \ei

!latex \item In each plasma volume
!latex           the boundary condition on the outer interface is $b=0$.
!latex \item In the vacuum volume (only for free-boundary), we may set $\mu=0$.

!latex \end{enumerate} \subsection{energy and helicity integrands} \begin{enumerate}

!latex \item The integrands are

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \sqrt g \; {\bf B} \cdot {\bf B} 

!latex & = & + &   & ( m_i \Azo{i,l} + n_i \Ato{i,l} ) & ( m_j \Azo{j,p} + n_j \Ato{j,p} ) & T_{l}        \; T_{p}        \; \bgss \; \cos\a_i \; \cos\a_j \\
!latex &   & - & 2 & ( m_i \Azo{i,l} + n_i \Ato{i,l} ) & ( m_j \Aze{j,p} + n_j \Ate{j,p} ) & T_{l}        \; T_{p}        \; \bgss \; \cos\a_i \; \sin\a_j \\
!latex &   & + &   & ( m_i \Aze{i,l} + n_i \Ate{i,l} ) & ( m_j \Aze{j,p} + n_j \Ate{j,p} ) & T_{l}        \; T_{p}        \; \bgss \; \sin\a_i \; \sin\a_j \\ 
!latex \\

!latex &   & - & 2 & ( m_i \Azo{i,l} + n_i \Ato{i,l} ) &       \Aze{j,p}                   & T_{l}        \; T_{p}^\prime \; \bgst \; \cos\a_i \; \cos\a_j \\
!latex &   & - & 2 & ( m_i \Azo{i,l} + n_i \Ato{i,l} ) &       \Azo{j,p}                   & T_{l}        \; T_{p}^\prime \; \bgst \; \cos\a_i \; \sin\a_j \\
!latex &   & + & 2 & ( m_i \Aze{i,l} + n_i \Ate{i,l} ) &       \Aze{j,p}                   & T_{l}        \; T_{p}^\prime \; \bgst \; \sin\a_i \; \cos\a_j \\
!latex &   & + & 2 & ( m_i \Aze{i,l} + n_i \Ate{i,l} ) &       \Azo{j,p}                   & T_{l}        \; T_{p}^\prime \; \bgst \; \sin\a_i \; \sin\a_j \\  
!latex \\

!latex &   & + & 2 & ( m_i \Azo{i,l} + n_i \Ato{i,l} ) &       \Ate{j,p}                   & T_{l}        \; T_{p}^\prime \; \bgsz \; \cos\a_i \; \cos\a_j \\
!latex &   & + & 2 & ( m_i \Azo{i,l} + n_i \Ato{i,l} ) &       \Ato{j,p}                   & T_{l}        \; T_{p}^\prime \; \bgsz \; \cos\a_i \; \sin\a_j \\
!latex &   & - & 2 & ( m_i \Aze{i,l} + n_i \Ate{i,l} ) &       \Ate{j,p}                   & T_{l}        \; T_{p}^\prime \; \bgsz \; \sin\a_i \; \cos\a_j \\
!latex &   & - & 2 & ( m_i \Aze{i,l} + n_i \Ate{i,l} ) &       \Ato{j,p}                   & T_{l}        \; T_{p}^\prime \; \bgsz \; \sin\a_i \; \sin\a_j \\  
!latex \\

!latex &   & + &   &       \Aze{i,l}                   &       \Aze{j,p}                   & T_{l}^\prime \; T_{p}^\prime \; \bgtt \; \cos\a_i \; \cos\a_j \\
!latex &   & + & 2 &       \Aze{i,l}                   &       \Azo{j,p}                   & T_{l}^\prime \; T_{p}^\prime \; \bgtt \; \cos\a_i \; \sin\a_j \\
!latex &   & + &   &       \Azo{i,l}                   &       \Azo{j,p}                   & T_{l}^\prime \; T_{p}^\prime \; \bgtt \; \sin\a_i \; \sin\a_j \\  
!latex \\

!latex &   & - & 2 &       \Aze{i,l}                   &       \Ate{j,p}                   & T_{l}^\prime \; T_{p}^\prime \; \bgtz \; \cos\a_i \; \cos\a_j \\
!latex &   & - & 2 &       \Aze{i,l}                   &       \Ato{j,p}                   & T_{l}^\prime \; T_{p}^\prime \; \bgtz \; \cos\a_i \; \sin\a_j \\
!latex &   & - & 2 &       \Azo{i,l}                   &       \Ate{j,p}                   & T_{l}^\prime \; T_{p}^\prime \; \bgtz \; \sin\a_i \; \cos\a_j \\
!latex &   & - & 2 &       \Azo{i,l}                   &       \Ato{j,p}                   & T_{l}^\prime \; T_{p}^\prime \; \bgtz \; \sin\a_i \; \sin\a_j \\  
!latex \\

!latex &   & + &   &       \Ate{i,l}                   &       \Ate{j,p}                   & T_{l}^\prime \; T_{p}^\prime \; \bgzz \; \cos\a_i \; \cos\a_j \\
!latex &   & + & 2 &       \Ate{i,l}                   &       \Ato{j,p}                   & T_{l}^\prime \; T_{p}^\prime \; \bgzz \; \cos\a_i \; \sin\a_j \\
!latex &   & + &   &       \Ato{i,l}                   &       \Ato{j,p}                   & T_{l}^\prime \; T_{p}^\prime \; \bgzz \; \sin\a_i \; \sin\a_j 

!latex \end{array}        
!latex \ee

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \sqrt g \; {\bf A} \cdot {\bf B} 
!latex & = & - &         \Aze{i,l} \; \Ate{j,p}    &         T_{l}^\prime \; T_{p} \; \cos\a_i \; \cos\a_j   \\
!latex &   & - &         \Aze{i,l} \; \Ato{j,p}    &         T_{l}^\prime \; T_{p} \; \cos\a_i \; \sin\a_j   \\
!latex &   & - &         \Azo{i,l} \; \Ate{j,p}    &         T_{l}^\prime \; T_{p} \; \sin\a_i \; \cos\a_j   \\
!latex &   & - &         \Azo{i,l} \; \Ato{j,p}    &         T_{l}^\prime \; T_{p} \; \sin\a_i \; \sin\a_j   \\ \\
!latex &   & + &         \Ate{i,l} \; \Aze{j,p}    &         T_{l}^\prime \; T_{p} \; \cos\a_i \; \cos\a_j   \\
!latex &   & + &         \Ate{i,l} \; \Azo{j,p}    &         T_{l}^\prime \; T_{p} \; \cos\a_i \; \sin\a_j   \\
!latex &   & + &         \Ato{i,l} \; \Aze{j,p}    &         T_{l}^\prime \; T_{p} \; \sin\a_i \; \cos\a_j   \\
!latex &   & + &         \Ato{i,l} \; \Azo{j,p}    &         T_{l}^\prime \; T_{p} \; \sin\a_i \; \sin\a_j             
!latex \end{array}        
!latex \ee

!latex \newpage

!latex \end{enumerate} \subsection{first derivatives of energy and helicity integrands with respect to $\Ate{i,l}$ and $\Ato{i,l}$} \begin{enumerate}

!latex \item The first derivatives with respect to $\Ate{i,l}$ and $\Ato{i,l}$ are

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \ds \frac{\partial}{\partial \Ate{i,l}} \int \!\! dv \; {\bf B} \cdot {\bf B} & = & \\
!latex & - & 2 & n_i & ( m_j \Azo{j,p} + n_j \Ato{j,p} ) & \int ds & T_{p}        \; T_{l}        \; \ooint \bgss \; \cos\a_j \; \sin\a_i \\
!latex & + & 2 & n_i & ( m_j \Aze{j,p} + n_j \Ate{j,p} ) & \int ds & T_{l}        \; T_{p}        \; \ooint \bgss \; \sin\a_i \; \sin\a_j \\
!latex & + & 2 & n_i &       \Aze{j,p}                   & \int ds & T_{l}        \; T_{p}^\prime \; \ooint \bgst \; \sin\a_i \; \cos\a_j \\
!latex & + & 2 & n_i &       \Azo{j,p}                   & \int ds & T_{l}        \; T_{p}^\prime \; \ooint \bgst \; \sin\a_i \; \sin\a_j \\
!latex & + & 2 &     & ( m_j \Azo{j,p} + n_j \Ato{j,p} ) & \int ds & T_{p}        \; T_{l}^\prime \; \ooint \bgsz \; \cos\a_j \; \cos\a_i \\
!latex & - & 2 & n_i &       \Ate{j,p}                   & \int ds & T_{l}        \; T_{p}^\prime \; \ooint \bgsz \; \sin\a_i \; \cos\a_j \\
!latex & - & 2 &     & ( m_j \Aze{j,p} + n_j \Ate{j,p} ) & \int ds & T_{p}        \; T_{l}^\prime \; \ooint \bgsz \; \sin\a_j \; \cos\a_i \\
!latex & - & 2 & n_i &       \Ato{j,p}                   & \int ds & T_{l}        \; T_{p}^\prime \; \ooint \bgsz \; \sin\a_i \; \sin\a_j \\
!latex & - & 2 &     &       \Aze{j,p}                   & \int ds & T_{p}^\prime \; T_{l}^\prime \; \ooint \bgtz \; \cos\a_j \; \cos\a_i \\
!latex & - & 2 &     &       \Azo{j,p}                   & \int ds & T_{p}^\prime \; T_{l}^\prime \; \ooint \bgtz \; \sin\a_j \; \cos\a_i \\
!latex & + & 2 &     &       \Ate{j,p}                   & \int ds & T_{l}^\prime \; T_{p}^\prime \; \ooint \bgzz \; \cos\a_i \; \cos\a_j \\
!latex & + & 2 &     &       \Ato{j,p}                   & \int ds & T_{l}^\prime \; T_{p}^\prime \; \ooint \bgzz \; \cos\a_i \; \sin\a_j 
!latex \end{array}        
!latex \ee

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \ds \frac{\partial}{\partial \Ate{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = & \\
!latex & - & \Aze{j,p} & \int ds & T_{p}^\prime \; T_{l}    \; \ooint \cos\a_j \; \cos\a_i \\
!latex & - & \Azo{j,p} & \int ds & T_{p}^\prime \; T_{l}    \; \ooint \sin\a_j \; \cos\a_i \\
!latex & + & \Aze{j,p} & \int ds & T_{l}^\prime \; T_{p}    \; \ooint \cos\a_i \; \cos\a_j \\
!latex & + & \Azo{j,p} & \int ds & T_{l}^\prime \; T_{p}    \; \ooint \cos\a_i \; \sin\a_j       
!latex \end{array}         
!latex \ee

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \ds \frac{\partial}{\partial \Ato{i,l}} \int \!\! dv \; {\bf B} \cdot {\bf B} & = & \\
!latex & + & 2 & n_i & ( m_j \Azo{j,p} + n_j \Ato{j,p} ) & \int ds & T_{l}        \; T_{p}        \; \ooint \bgss \; \cos\a_i \; \cos\a_j \\
!latex & - & 2 & n_i & ( m_j \Aze{j,p} + n_j \Ate{j,p} ) & \int ds & T_{l}        \; T_{p}        \; \ooint \bgss \; \cos\a_i \; \sin\a_j \\
!latex & - & 2 & n_i &       \Aze{j,p}                   & \int ds & T_{l}        \; T_{p}^\prime \; \ooint \bgst \; \cos\a_i \; \cos\a_j \\
!latex & - & 2 & n_i &       \Azo{j,p}                   & \int ds & T_{l}        \; T_{p}^\prime \; \ooint \bgst \; \cos\a_i \; \sin\a_j \\
!latex & + & 2 & n_i &       \Ate{j,p}                   & \int ds & T_{l}        \; T_{p}^\prime \; \ooint \bgsz \; \cos\a_i \; \cos\a_j \\
!latex & + & 2 & n_i &       \Ato{j,p}                   & \int ds & T_{l}        \; T_{p}^\prime \; \ooint \bgsz \; \cos\a_i \; \sin\a_j \\
!latex & + & 2 &     & ( m_j \Azo{j,p} + n_j \Ato{j,p} ) & \int ds & T_{p}        \; T_{l}^\prime \; \ooint \bgsz \; \cos\a_j \; \sin\a_i \\
!latex & - & 2 &     & ( m_j \Aze{j,p} + n_j \Ate{j,p} ) & \int ds & T_{p}        \; T_{l}^\prime \; \ooint \bgsz \; \sin\a_j \; \sin\a_i \\
!latex & - & 2 &     &       \Aze{j,p}                   & \int ds & T_{p}^\prime \; T_{l}^\prime \; \ooint \bgtz \; \cos\a_j \; \sin\a_i \\
!latex & - & 2 &     &       \Azo{j,p}                   & \int ds & T_{p}^\prime \; T_{l}^\prime \; \ooint \bgtz \; \sin\a_j \; \sin\a_i \\
!latex & + & 2 &     &       \Ate{j,p}                   & \int ds & T_{p}^\prime \; T_{l}^\prime \; \ooint \bgzz \; \cos\a_j \; \sin\a_i \\
!latex & + & 2 &     &       \Ato{j,p}                   & \int ds & T_{l}^\prime \; T_{p}^\prime \; \ooint \bgzz \; \sin\a_i \; \sin\a_j 
!latex \end{array}
!latex \ee

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \ds \frac{\partial}{\partial \Ato{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = & \\
!latex & - & \Aze{j,p} & \int ds & T_{p}^\prime \; T_{l}        \; \ooint \cos\a_j \; \sin\a_i \\
!latex & - & \Azo{j,p} & \int ds & T_{p}^\prime \; T_{l}        \; \ooint \sin\a_j \; \sin\a_i \\
!latex & + & \Aze{j,p} & \int ds & T_{l}^\prime \; T_{p}        \; \ooint \sin\a_i \; \cos\a_j \\
!latex & + & \Azo{j,p} & \int ds & T_{l}^\prime \; T_{p}        \; \ooint \sin\a_i \; \sin\a_j       
!latex \end{array}         
!latex \ee

!latex \newpage

!latex \end{enumerate} \subsection{first derivatives of energy and helicity integrands with respect to $\Aze{i,l}$ and $\Azo{i,l}$} \begin{enumerate}

!latex \item The first derivatives with respect to $\Aze{i,l}$ and $\Azo{i,l}$ are

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \ds \frac{\partial}{\partial \Aze{i,l}} \int \!\! dv \; {\bf B} \cdot {\bf B} & = & \\
!latex & - & 2 & m_i & ( m_j \Azo{j,p} + n_j \Ato{j,p} ) & \int ds & T_{p}        \; T_{l}        \; \ooint \bgss \; \cos\a_j \; \sin\a_i \\
!latex & + & 2 & m_i & ( m_j \Aze{j,p} + n_j \Ate{j,p} ) & \int ds & T_{l}        \; T_{p}        \; \ooint \bgss \; \sin\a_i \; \sin\a_j \\
!latex & - & 2 &     & ( m_j \Azo{j,p} + n_j \Ato{j,p} ) & \int ds & T_{p}        \; T_{l}^\prime \; \ooint \bgst \; \cos\a_j \; \cos\a_i \\
!latex & + & 2 & m_i &       \Aze{j,p}                   & \int ds & T_{l}        \; T_{p}^\prime \; \ooint \bgst \; \sin\a_i \; \cos\a_j \\
!latex & + & 2 &     & ( m_j \Aze{j,p} + n_j \Ate{j,p} ) & \int ds & T_{p}        \; T_{l}^\prime \; \ooint \bgst \; \sin\a_j \; \cos\a_i \\
!latex & + & 2 & m_i &       \Azo{j,p}                   & \int ds & T_{l}        \; T_{p}^\prime \; \ooint \bgst \; \sin\a_i \; \sin\a_j \\
!latex & - & 2 & m_i &       \Ate{j,p}                   & \int ds & T_{l}        \; T_{p}^\prime \; \ooint \bgsz \; \sin\a_i \; \cos\a_j \\
!latex & - & 2 & m_i &       \Ato{j,p}                   & \int ds & T_{l}        \; T_{p}^\prime \; \ooint \bgsz \; \sin\a_i \; \sin\a_j \\
!latex & + & 2 &     &       \Aze{j,p}                   & \int ds & T_{l}^\prime \; T_{p}^\prime \; \ooint \bgtt \; \cos\a_i \; \cos\a_j \\
!latex & + & 2 &     &       \Azo{j,p}                   & \int ds & T_{l}^\prime \; T_{p}^\prime \; \ooint \bgtt \; \cos\a_i \; \sin\a_j \\
!latex & - & 2 &     &       \Ate{j,p}                   & \int ds & T_{l}^\prime \; T_{p}^\prime \; \ooint \bgtz \; \cos\a_i \; \cos\a_j \\
!latex & - & 2 &     &       \Ato{j,p}                   & \int ds & T_{l}^\prime \; T_{p}^\prime \; \ooint \bgtz \; \cos\a_i \; \sin\a_j 
!latex \end{array}        
!latex \ee

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \ds \frac{\partial}{\partial \Aze{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = & \\
!latex & - & \Ate{j,p} & \int ds & T_{l}^\prime \; T_{p}        \; \ooint \cos\a_i \; \cos\a_j \\
!latex & - & \Ato{j,p} & \int ds & T_{l}^\prime \; T_{p}        \; \ooint \cos\a_i \; \sin\a_j \\
!latex & + & \Ate{j,p} & \int ds & T_{p}^\prime \; T_{l}        \; \ooint \cos\a_j \; \cos\a_i \\
!latex & + & \Ato{j,p} & \int ds & T_{p}^\prime \; T_{l}        \; \ooint \sin\a_j \; \cos\a_i       
!latex \end{array}         
!latex \ee

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \ds \frac{\partial}{\partial \Azo{i,l}} \int \!\! dv \; {\bf B} \cdot {\bf B} & = & \\
!latex & + & 2 & m_i & ( m_j \Azo{j,p} + n_j \Ato{j,p} ) & \int ds & T_{l}        \; T_{p}        \; \ooint \bgss \; \cos\a_i \; \cos\a_j \\
!latex & - & 2 & m_i & ( m_j \Aze{j,p} + n_j \Ate{j,p} ) & \int ds & T_{l}        \; T_{p}        \; \ooint \bgss \; \cos\a_i \; \sin\a_j \\
!latex & - & 2 & m_i &       \Aze{j,p}                   & \int ds & T_{l}        \; T_{p}^\prime \; \ooint \bgst \; \cos\a_i \; \cos\a_j \\
!latex & - & 2 & m_i &       \Azo{j,p}                   & \int ds & T_{l}        \; T_{p}^\prime \; \ooint \bgst \; \cos\a_i \; \sin\a_j \\
!latex & - & 2 &     & ( m_j \Azo{j,p} + n_j \Ato{j,p} ) & \int ds & T_{p}        \; T_{l}^\prime \; \ooint \bgst \; \cos\a_j \; \sin\a_i \\
!latex & + & 2 &     & ( m_j \Aze{j,p} + n_j \Ate{j,p}   & \int ds & T_{p}        \; T_{l}^\prime \; \ooint \bgst \; \sin\a_j \; \sin\a_i \\
!latex & + & 2 & m_i &       \Ate{j,p}                   & \int ds & T_{l}        \; T_{p}^\prime \; \ooint \bgsz \; \cos\a_i \; \cos\a_j \\
!latex & + & 2 & m_i &       \Ato{j,p}                   & \int ds & T_{l}        \; T_{p}^\prime \; \ooint \bgsz \; \cos\a_i \; \sin\a_j \\
!latex & + & 2 &     &       \Aze{j,p}                   & \int ds & T_{p}^\prime \; T_{l}^\prime \; \ooint \bgtt \; \cos\a_j \; \sin\a_i \\
!latex & + & 2 &     &       \Azo{j,p}                   & \int ds & T_{l}^\prime \; T_{p}^\prime \; \ooint \bgtt \; \sin\a_i \; \sin\a_j \\
!latex & - & 2 &     &       \Ate{j,p}                   & \int ds & T_{l}^\prime \; T_{p}^\prime \; \ooint \bgtz \; \sin\a_i \; \cos\a_j \\
!latex & - & 2 &     &       \Ato{j,p}                   & \int ds & T_{l}^\prime \; T_{p}^\prime \; \ooint \bgtz \; \sin\a_i \; \sin\a_j 
!latex \end{array}        
!latex \ee

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \ds \frac{\partial}{\partial \Azo{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B} & = & \\
!latex & - & \Ate{j,p} & \int ds & T_{l}^\prime \; T_{p}        \; \ooint \sin\a_i \; \cos\a_j \\
!latex & - & \Ato{j,p} & \int ds & T_{l}^\prime \; T_{p}        \; \ooint \sin\a_i \; \sin\a_j \\
!latex & + & \Ate{j,p} & \int ds & T_{p}^\prime \; T_{l}        \; \ooint \cos\a_j \; \sin\a_i \\
!latex & + & \Ato{j,p} & \int ds & T_{p}^\prime \; T_{l}        \; \ooint \sin\a_j \; \sin\a_i       
!latex \end{array}         
!latex \ee

!latex \newpage

!latex \end{enumerate} \subsection{second derivatives of energy and helicity integrands} \begin{enumerate}

!latex \item The second derivatives wrt $\Ate{i,l}$ (stellarator symmetric) are

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \ds \frac{\partial}{\partial \Ate{j,p}} \frac{\partial}{\partial \Ate{i,l}} \int \!\! dv \; {\bf B} \cdot {\bf B}  & = & \\
!latex & + & 2 & n_i & n_j & \int ds & T_{l}        \; T_{p}        \; \ooint \bgss \; \sin\a_i \; \sin\a_j \\
!latex & - & 2 & n_i &     & \int ds & T_{l}        \; T_{p}^\prime \; \ooint \bgsz \; \sin\a_i \; \cos\a_j \\
!latex & - & 2 &     & n_j & \int ds & T_{p}        \; T_{l}^\prime \; \ooint \bgsz \; \sin\a_j \; \cos\a_i \\
!latex & + & 2 &     &     & \int ds & T_{l}^\prime \; T_{p}^\prime \; \ooint \bgzz \; \cos\a_i \; \cos\a_j
!latex \end{array}
!latex \ee

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \ds \frac{\partial}{\partial \Ato{j,p}} \frac{\partial}{\partial \Ate{i,l}} \int \!\! dv \; {\bf B} \cdot {\bf B}  & = & \\
!latex & - & 2 & n_i & n_j & \int ds & T_{p}        \; T_{l}        \; \ooint \bgss \; \cos\a_j \; \sin\a_i \; ! \\
!latex & + & 2 &     & n_j & \int ds & T_{p}        \; T_{l}^\prime \; \ooint \bgsz \; \cos\a_j \; \cos\a_i \\
!latex & - & 2 & n_i &     & \int ds & T_{l}        \; T_{p}^\prime \; \ooint \bgsz \; \sin\a_i \; \sin\a_j \\
!latex & + & 2 &     &     & \int ds & T_{l}^\prime \; T_{p}^\prime \; \ooint \bgzz \; \cos\a_i \; \sin\a_j \; !
!latex \end{array}        
!latex \ee

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \ds \frac{\partial}{\partial \Aze{j,p}} \frac{\partial}{\partial \Ate{i,l}} \int \!\! dv \; {\bf B} \cdot {\bf B}  & = & \\
!latex & + & 2 & n_i & m_j & \int ds & T_{l}        \; T_{p}        \; \ooint \bgss \; \sin\a_i \; \sin\a_j \\
!latex & + & 2 & n_i &     & \int ds & T_{l}        \; T_{p}^\prime \; \ooint \bgst \; \sin\a_i \; \cos\a_j \\
!latex & - & 2 &     & m_j & \int ds & T_{p}        \; T_{l}^\prime \; \ooint \bgsz \; \sin\a_j \; \cos\a_i \\
!latex & - & 2 &     &     & \int ds & T_{p}^\prime \; T_{l}^\prime \; \ooint \bgtz \; \cos\a_j \; \cos\a_i
!latex \end{array}        
!latex \ee

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \ds \frac{\partial}{\partial \Azo{j,p}} \frac{\partial}{\partial \Ate{i,l}} \int \!\! dv \; {\bf B} \cdot {\bf B}  & = & \\
!latex & - & 2 & n_i & m_j & \int ds & T_{p}        \; T_{l}        \; \ooint \bgss \; \cos\a_j \; \sin\a_i \; ! \\
!latex & + & 2 & n_i &     & \int ds & T_{l}        \; T_{p}^\prime \; \ooint \bgst \; \sin\a_i \; \sin\a_j \\
!latex & + & 2 &     & m_j & \int ds & T_{p}        \; T_{l}^\prime \; \ooint \bgsz \; \cos\a_j \; \cos\a_i \\
!latex & - & 2 &     &     & \int ds & T_{p}^\prime \; T_{l}^\prime \; \ooint \bgtz \; \sin\a_j \; \cos\a_i
!latex \end{array}        
!latex \ee

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \ds \frac{\partial}{\partial \Ate{j,p}} \frac{\partial}{\partial \Ate{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B}  & = & 0
!latex \end{array}        
!latex \ee

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \ds \frac{\partial}{\partial \Ato{j,p}} \frac{\partial}{\partial \Ate{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B}  & = & 0
!latex \end{array}        
!latex \ee

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \ds \frac{\partial}{\partial \Aze{j,p}} \frac{\partial}{\partial \Ate{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B}  & = & \\
!latex & - & \int ds & T_{p}^\prime \; T_{l}        \; \ooint \cos\a_j \; \cos\a_i \\
!latex & + & \int ds & T_{l}^\prime \; T_{p}        \; \ooint \cos\a_i \; \cos\a_j   
!latex \end{array}        
!latex \ee

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \ds \frac{\partial}{\partial \Azo{j,p}} \frac{\partial}{\partial \Ate{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B}  & = & \\
!latex & - & \int ds & T_{p}^\prime \; T_{l}        \; \ooint \sin\a_j \; \cos\a_i \\
!latex & + & \int ds & T_{l}^\prime \; T_{p}        \; \ooint \cos\a_i \; \sin\a_j \; !
!latex \end{array}        
!latex \ee

!latex \newpage

!latex \end{enumerate} \subsection{second derivatives of energy and helicity integrands} \begin{enumerate}

!latex \item The second derivatives wrt $\Ato{i,l}$ (non-stellarator symmetric) are

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \ds \frac{\partial}{\partial \Ate{j,p}} \frac{\partial}{\partial \Ato{i,l}} \int \!\! dv \; {\bf B} \cdot {\bf B}  & = & \\
!latex & - & 2 & n_i & n_j & \int ds & T_{l}        \; T_{p}        \; \ooint \bgss \; \cos\a_i \; \sin\a_j \; ! \\
!latex & + & 2 & n_i &     & \int ds & T_{l}        \; T_{p}^\prime \; \ooint \bgsz \; \cos\a_i \; \cos\a_j \\
!latex & - & 2 &     & n_j & \int ds & T_{p}        \; T_{l}^\prime \; \ooint \bgsz \; \sin\a_j \; \sin\a_i \\
!latex & + & 2 &     &     & \int ds & T_{p}^\prime \; T_{l}^\prime \; \ooint \bgzz \; \cos\a_j \; \sin\a_i \; !
!latex \end{array}        
!latex \ee

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \ds \frac{\partial}{\partial \Ato{j,p}} \frac{\partial}{\partial \Ato{i,l}} \int \!\! dv \; {\bf B} \cdot {\bf B}  & = & \\
!latex & + & 2 & n_i & n_j & \int ds & T_{l}        \; T_{p}        \; \ooint \bgss \; \cos\a_i \; \cos\a_j \\
!latex & + & 2 & n_i &     & \int ds & T_{l}        \; T_{p}^\prime \; \ooint \bgsz \; \cos\a_i \; \sin\a_j \; ! \\
!latex & + & 2 &     & n_j & \int ds & T_{p}        \; T_{l}^\prime \; \ooint \bgsz \; \cos\a_j \; \sin\a_i \; ! \\
!latex & + & 2 &     &     & \int ds & T_{l}^\prime \; T_{p}^\prime \; \ooint \bgzz \; \sin\a_i \; \sin\a_j
!latex \end{array}        
!latex \ee

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \ds \frac{\partial}{\partial \Aze{j,p}} \frac{\partial}{\partial \Ato{i,l}} \int \!\! dv \; {\bf B} \cdot {\bf B}  & = & \\
!latex & - & 2 & n_i & m_j & \int ds & T_{l}        \; T_{p}        \; \ooint \bgss \; \cos\a_i \; \sin\a_j \; ! \\
!latex & - & 2 & n_i &     & \int ds & T_{l}        \; T_{p}^\prime \; \ooint \bgst \; \cos\a_i \; \cos\a_j \\
!latex & - & 2 &     & m_j & \int ds & T_{p}        \; T_{l}^\prime \; \ooint \bgsz \; \sin\a_j \; \sin\a_i \\
!latex & - & 2 &     &     & \int ds & T_{p}^\prime \; T_{l}^\prime \; \ooint \bgtz \; \cos\a_j \; \sin\a_i \; !
!latex \end{array}        
!latex \ee

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \ds \frac{\partial}{\partial \Azo{j,p}} \frac{\partial}{\partial \Ato{i,l}} \int \!\! dv \; {\bf B} \cdot {\bf B}  & = & \\
!latex & + & 2 & n_i & m_j & \int ds & T_{l}        \; T_{p}        \; \ooint \bgss \; \cos\a_i \; \cos\a_j \\
!latex & - & 2 & n_i &     & \int ds & T_{l}        \; T_{p}^\prime \; \ooint \bgst \; \cos\a_i \; \sin\a_j \; ! \\
!latex & + & 2 &     & m_j & \int ds & T_{p}        \; T_{l}^\prime \; \ooint \bgsz \; \cos\a_j \; \sin\a_i \; ! \\
!latex & - & 2 &     &     & \int ds & T_{p}^\prime \; T_{l}^\prime \; \ooint \bgtz \; \sin\a_j \; \sin\a_i
!latex \end{array}        
!latex \ee

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \ds \frac{\partial}{\partial \Ate{j,p}} \frac{\partial}{\partial \Ato{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B}  & = & 0
!latex \end{array}        
!latex \ee

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \ds \frac{\partial}{\partial \Ato{j,p}} \frac{\partial}{\partial \Ato{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B}  & = & 0
!latex \end{array}        
!latex \ee

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \ds \frac{\partial}{\partial \Aze{j,p}} \frac{\partial}{\partial \Ato{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B}  & = & \\
!latex & - & \int ds & T_{p}^\prime \; T_{l}        \; \ooint \cos\a_j \; \sin\a_i \; ! \\
!latex & + & \int ds & T_{l}^\prime \; T_{p}        \; \ooint \sin\a_i \; \cos\a_j   
!latex \end{array}        
!latex \ee

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \ds \frac{\partial}{\partial \Azo{j,p}} \frac{\partial}{\partial \Ato{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B}  & = & \\
!latex & - & \int ds & T_{p}^\prime \; T_{l}        \; \ooint \sin\a_j \; \sin\a_i \\
!latex & + & \int ds & T_{l}^\prime \; T_{p}        \; \ooint \sin\a_i \; \sin\a_j   
!latex \end{array}        
!latex \ee

!latex \newpage

!latex \end{enumerate} \subsection{second derivatives of energy and helicity integrands} \begin{enumerate}

!latex \item The second derivatives wrt $\Aze{i,l}$ (stellarator symmetric) are

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \ds \frac{\partial}{\partial \Ate{j,p}} \frac{\partial}{\partial \Aze{i,l}} \int \!\! dv \; {\bf B} \cdot {\bf B}  & = & \\
!latex & + & 2 & m_i & n_j & \int ds & T_{l}        \; T_{p}        \; \ooint \bgss \; \sin\a_i \; \sin\a_j \\
!latex & + & 2 &     & n_j & \int ds & T_{p}        \; T_{l}^\prime \; \ooint \bgst \; \sin\a_j \; \cos\a_i \\
!latex & - & 2 & m_i &     & \int ds & T_{l}        \; T_{p}^\prime \; \ooint \bgsz \; \sin\a_i \; \cos\a_j \\
!latex & - & 2 &     &     & \int ds & T_{l}^\prime \; T_{p}^\prime \; \ooint \bgtz \; \cos\a_i \; \cos\a_j
!latex \end{array}        
!latex \ee

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \ds \frac{\partial}{\partial \Ato{j,p}} \frac{\partial}{\partial \Aze{i,l}} \int \!\! dv \; {\bf B} \cdot {\bf B}  & = & \\
!latex & - & 2 & m_i & n_j & \int ds & T_{p}        \; T_{l}        \; \ooint \bgss \; \cos\a_j \; \sin\a_i \; ! \\
!latex & - & 2 &     & n_j & \int ds & T_{p}        \; T_{l}^\prime \; \ooint \bgst \; \cos\a_j \; \cos\a_i \\
!latex & - & 2 & m_i &     & \int ds & T_{l}        \; T_{p}^\prime \; \ooint \bgsz \; \sin\a_i \; \sin\a_j \\
!latex & - & 2 &     &     & \int ds & T_{l}^\prime \; T_{p}^\prime \; \ooint \bgtz \; \cos\a_i \; \sin\a_j \; !
!latex \end{array}        
!latex \ee

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \ds \frac{\partial}{\partial \Aze{j,p}} \frac{\partial}{\partial \Aze{i,l}} \int \!\! dv \; {\bf B} \cdot {\bf B}  & = & \\
!latex & + & 2 & m_i & m_j & \int ds & T_{l}        \; T_{p}        \; \ooint \bgss \; \sin\a_i \; \sin\a_j \\
!latex & + & 2 & m_i &     & \int ds & T_{l}        \; T_{p}^\prime \; \ooint \bgst \; \sin\a_i \; \cos\a_j \\
!latex & + & 2 &     & m_j & \int ds & T_{p}        \; T_{l}^\prime \; \ooint \bgst \; \sin\a_j \; \cos\a_i \\
!latex & + & 2 &     &     & \int ds & T_{l}^\prime \; T_{p}^\prime \; \ooint \bgtt \; \cos\a_i \; \cos\a_j
!latex \end{array}        
!latex \ee

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \ds \frac{\partial}{\partial \Azo{j,p}} \frac{\partial}{\partial \Aze{i,l}} \int \!\! dv \; {\bf B} \cdot {\bf B}  & = & \\
!latex & - & 2 & m_i & m_j & \int ds & T_{p}        \; T_{l}        \; \ooint \bgss \; \cos\a_j \; \sin\a_i \; ! \\
!latex & - & 2 &     & m_j & \int ds & T_{p}        \; T_{l}^\prime \; \ooint \bgst \; \cos\a_j \; \cos\a_i \\
!latex & + & 2 & m_i &     & \int ds & T_{l}        \; T_{p}^\prime \; \ooint \bgst \; \sin\a_i \; \sin\a_j \\
!latex & + & 2 &     &     & \int ds & T_{l}^\prime \; T_{p}^\prime \; \ooint \bgtt \; \cos\a_i \; \sin\a_j \; !
!latex \end{array}        
!latex \ee

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \ds \frac{\partial}{\partial \Ate{j,p}} \frac{\partial}{\partial \Aze{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B}  & = & \\
!latex & - & \int ds & T_{l}^\prime \; T_{p}        \; \ooint \cos\a_i \; \cos\a_j \\
!latex & + & \int ds & T_{p}^\prime \; T_{l}        \; \ooint \cos\a_j \; \cos\a_i   
!latex \end{array}        
!latex \ee

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \ds \frac{\partial}{\partial \Ato{j,p}} \frac{\partial}{\partial \Aze{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B}  & = & \\
!latex & - & \int ds & T_{l}^\prime \; T_{p}        \; \ooint \cos\a_i \; \sin\a_j \; ! \\
!latex & + & \int ds & T_{p}^\prime \; T_{l}        \; \ooint \sin\a_j \; \cos\a_i   
!latex \end{array}        
!latex \ee

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \ds \frac{\partial}{\partial \Aze{j,p}} \frac{\partial}{\partial \Aze{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B}  & = & 0
!latex \end{array}        
!latex \ee

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \ds \frac{\partial}{\partial \Azo{j,p}} \frac{\partial}{\partial \Aze{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B}  & = & 0
!latex \end{array}        
!latex \ee

!latex \newpage

!latex \end{enumerate} \subsection{second derivatives of energy and helicity integrands} \begin{enumerate}

!latex \item The second derivatives wrt $\Azo{i,l}$ (non-stellarator symmetric) are

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \ds \frac{\partial}{\partial \Ate{j,p}} \frac{\partial}{\partial \Azo{i,l}} \int \!\! dv \; {\bf B} \cdot {\bf B}  & = & \\
!latex & - & 2 & m_i & n_j & \int ds & T_{l}        \; T_{p}        \; \ooint \bgss \; \cos\a_i \; \sin\a_j \; ! \\
!latex & + & 2 &     & n_j & \int ds & T_{p}        \; T_{l}^\prime \; \ooint \bgst \; \sin\a_j \; \sin\a_i \\
!latex & + & 2 & m_i &     & \int ds & T_{l}        \; T_{p}^\prime \; \ooint \bgsz \; \cos\a_i \; \cos\a_j \\
!latex & - & 2 &     &     & \int ds & T_{l}^\prime \; T_{p}^\prime \; \ooint \bgtz \; \sin\a_i \; \cos\a_j
!latex \end{array}        
!latex \ee

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \ds \frac{\partial}{\partial \Ato{j,p}} \frac{\partial}{\partial \Azo{i,l}} \int \!\! dv \; {\bf B} \cdot {\bf B}  & = & \\
!latex & + & 2 & m_i & n_j & \int ds & T_{l}        \; T_{p}        \; \ooint \bgss \; \cos\a_i \; \cos\a_j \\
!latex & - & 2 &     & n_j & \int ds & T_{p}        \; T_{l}^\prime \; \ooint \bgst \; \cos\a_j \; \sin\a_i \; ! \\
!latex & + & 2 & m_i &     & \int ds & T_{l}        \; T_{p}^\prime \; \ooint \bgsz \; \cos\a_i \; \sin\a_j \; ! \\
!latex & - & 2 &     &     & \int ds & T_{l}^\prime \; T_{p}^\prime \; \ooint \bgtz \; \sin\a_i \; \sin\a_j
!latex \end{array}        
!latex \ee

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \ds \frac{\partial}{\partial \Aze{j,p}} \frac{\partial}{\partial \Azo{i,l}} \int \!\! dv \; {\bf B} \cdot {\bf B}  & = & \\
!latex & - & 2 & m_i & m_j & \int ds & T_{l}        \; T_{p}        \; \ooint \bgss \; \cos\a_i \; \sin\a_j \; ! \\
!latex & - & 2 & m_i &     & \int ds & T_{l}        \; T_{p}^\prime \; \ooint \bgst \; \cos\a_i \; \cos\a_j \\
!latex & + & 2 &     & m_j & \int ds & T_{p}        \; T_{l}^\prime \; \ooint \bgst \; \sin\a_j \; \sin\a_i \\
!latex & + & 2 &     &     & \int ds & T_{p}^\prime \; T_{l}^\prime \; \ooint \bgtt \; \cos\a_j \; \sin\a_i \; !
!latex \end{array}        
!latex \ee

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \ds \frac{\partial}{\partial \Azo{j,p}} \frac{\partial}{\partial \Azo{i,l}} \int \!\! dv \; {\bf B} \cdot {\bf B}  & = & \\
!latex & + & 2 & m_i & m_j & \int ds & T_{l}        \; T_{p}        \; \ooint \bgss \; \cos\a_i \; \cos\a_j \\
!latex & - & 2 & m_i &     & \int ds & T_{l}        \; T_{p}^\prime \; \ooint \bgst \; \cos\a_i \; \sin\a_j \; ! \\
!latex & - & 2 &     & m_j & \int ds & T_{p}        \; T_{l}^\prime \; \ooint \bgst \; \cos\a_j \; \sin\a_i \; ! \\
!latex & + & 2 &     &     & \int ds & T_{l}^\prime \; T_{p}^\prime \; \ooint \bgtt \; \sin\a_i \; \sin\a_j
!latex \end{array}        
!latex \ee

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \ds \frac{\partial}{\partial \Ate{j,p}} \frac{\partial}{\partial \Azo{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B}  & = & \\
!latex & - & \int ds & T_{l}^\prime \; T_{p}        \; \ooint \sin\a_i \; \cos\a_j \\
!latex & + & \int ds & T_{p}^\prime \; T_{l}        \; \ooint \cos\a_j \; \sin\a_i \; !
!latex \end{array}        
!latex \ee

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \ds \frac{\partial}{\partial \Ato{j,p}} \frac{\partial}{\partial \Azo{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B}  & = & \\
!latex & - & \int ds & T_{l}^\prime \; T_{p}        \; \ooint \sin\a_i \; \sin\a_j \\
!latex & + & \int ds & T_{p}^\prime \; T_{l}        \; \ooint \sin\a_j \; \sin\a_i   
!latex \end{array}        
!latex \ee

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \ds \frac{\partial}{\partial \Aze{j,p}} \frac{\partial}{\partial \Azo{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B}  & = & 0
!latex \end{array}        
!latex \ee

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \ds \frac{\partial}{\partial \Azo{j,p}} \frac{\partial}{\partial \Azo{i,l}} \int \!\! dv \; {\bf A} \cdot {\bf B}  & = & 0
!latex \end{array}        
!latex \ee

!latex \end{enumerate} \subsection{matrix elements} \begin{enumerate}

!latex \item The energy, $W \equiv \int \! dv {\; \bf B}\cdot{\bf B}$, and helicity, $K\equiv \int \! dv \; {\bf A}\cdot{\bf B}$, functionals may be written
!latex \be W & = & \frac{1}{2} \; a_i \; A_{i,j} \; a_j + a_i \; B_{i,j} \; \psi_j + \frac{1}{2} \; \psi_i \; C_{i,j} \; \psi_j \label{eq:energymatrix} \\
!latex     K & = & \frac{1}{2} \; a_i \; D_{i,j} \; a_j + a_i \; E_{i,j} \; \psi_j + \frac{1}{2} \; \psi_i \; F_{i,j} \; \psi_j \label{eq:helicitymatrix}
!latex \ee
!latex       where ${\bf a} \equiv \{ A_{\t,e,i,l}, A_{\z,e,i,l}, A_{\t,o,i,l}, A_{\z,o,i,l}, f_{e,i}, f_{o,i} \}$ contains the independent degrees of freedom
!latex       and $\boldpsi \equiv \{\Delta \psi_t,\Delta \psi_p\}$.

!latex \item The matrix elements are computed via
!latex       \be \verb+MA(i,j)+ & \equiv & A_{i,j} = \frac{\partial^2 W}{\partial    a_i \partial    a_j} \\
!latex           \verb+MB(i,j)+ & \equiv & B_{i,j} = \frac{\partial^2 W}{\partial    a_i \partial \psi_j} \\
!latex           \verb+MC(i,j)+ & \equiv & C_{i,j} = \frac{\partial^2 W}{\partial \psi_i \partial \psi_j}
!latex       \ee

!latex \item The energy functionals can also be represented as
!latex       \be \begin{array}{ccccccccccccccccccccccccccccc}
!latex           W & = & \frac{1}{2} & {\bf a}^T & \cdot & A[{\bf x}] & \cdot & {\bf a} & + & B[{\bf x}, \boldpsi] & \cdot & {\bf a} & + & C[{\bf x}, \boldpsi], \\
!latex           K & = & \frac{1}{2} & {\bf a}^T & \cdot & D          & \cdot & {\bf a} & + & E[         \boldpsi] & \cdot & {\bf a} & + & F[         \boldpsi].
!latex           \end{array}
!latex       \ee

!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine ma01ag( lvol, mn, lrad )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half, quart, one

  use numerical, only : sqrtmachprec

  use fileunits, only : ounit

  use inputlist, only : Wmacros, Wma01ag

  use cputiming, only : Tma01ag

  use allglobal, only : myid, cpus, &
                        NOTstellsym, &
                        im, in, &
                        NAdof, &
                        dMA, dMB, dMC, dMD, dME, dMF, &
                        Ate, Ato, Aze, Azo, Fso, Fse, &
                        DToocc, DToocs, DToosc, DTooss, &
                        TTsscc, TTsscs, TTsssc, TTssss, &
                        TDstcc, TDstcs, TDstsc, TDstss, &
                        TDszcc, TDszcs, TDszsc, TDszss, &
                        DDttcc, DDttcs, DDttsc, DDttss, &
                        DDtzcc, DDtzcs, DDtzsc, DDtzss, &
                        DDzzcc, DDzzcs, DDzzsc, DDzzss, &
                        Lcoordinatesingularity
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER, intent(in)  :: lvol, mn, lrad

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER              :: NN, ii, jj, ll, pp, mi, ni, mj, nj, mimj, minj, nimj, ninj, mjmi, mjni, njmi, njni, id, jd, Lul, Lup, Lzl, Lzp

  REAL                 :: Wtete(1:mn,0:lrad,1:mn,0:lrad), Wteto(1:mn,0:lrad,1:mn,0:lrad), Wtote(1:mn,0:lrad,1:mn,0:lrad), Wtoto(1:mn,0:lrad,1:mn,0:lrad)
  REAL                 :: Wteze(1:mn,0:lrad,1:mn,0:lrad), Wtezo(1:mn,0:lrad,1:mn,0:lrad), Wtoze(1:mn,0:lrad,1:mn,0:lrad), Wtozo(1:mn,0:lrad,1:mn,0:lrad)
  REAL                 :: Wzete(1:mn,0:lrad,1:mn,0:lrad), Wzeto(1:mn,0:lrad,1:mn,0:lrad), Wzote(1:mn,0:lrad,1:mn,0:lrad), Wzoto(1:mn,0:lrad,1:mn,0:lrad)
  REAL                 :: Wzeze(1:mn,0:lrad,1:mn,0:lrad), Wzezo(1:mn,0:lrad,1:mn,0:lrad), Wzoze(1:mn,0:lrad,1:mn,0:lrad), Wzozo(1:mn,0:lrad,1:mn,0:lrad)

  REAL                 :: Htete(1:mn,0:lrad,1:mn,0:lrad), Hteto(1:mn,0:lrad,1:mn,0:lrad), Htote(1:mn,0:lrad,1:mn,0:lrad), Htoto(1:mn,0:lrad,1:mn,0:lrad)
  REAL                 :: Hteze(1:mn,0:lrad,1:mn,0:lrad), Htezo(1:mn,0:lrad,1:mn,0:lrad), Htoze(1:mn,0:lrad,1:mn,0:lrad), Htozo(1:mn,0:lrad,1:mn,0:lrad)
  REAL                 :: Hzete(1:mn,0:lrad,1:mn,0:lrad), Hzeto(1:mn,0:lrad,1:mn,0:lrad), Hzote(1:mn,0:lrad,1:mn,0:lrad), Hzoto(1:mn,0:lrad,1:mn,0:lrad)
  REAL                 :: Hzeze(1:mn,0:lrad,1:mn,0:lrad), Hzezo(1:mn,0:lrad,1:mn,0:lrad), Hzoze(1:mn,0:lrad,1:mn,0:lrad), Hzozo(1:mn,0:lrad,1:mn,0:lrad)

  LOGICAL              :: Lzi, Lzj

  BEGIN(ma01ag)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef DEBUG
  FATAL(ma01ag, .not.allocated(dMA), error)
  FATAL(ma01ag, .not.allocated(dMB), error)
  FATAL(ma01ag, .not.allocated(dMC), error)
  FATAL(ma01ag, .not.allocated(dMD), error)
  FATAL(ma01ag, .not.allocated(dME), error)
  FATAL(ma01ag, .not.allocated(dMF), error)
#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  NN = NAdof(lvol) ! shorthand;

  dMA(0:NN,0:NN) = zero ! initialize summation; 24 Jan 13;
  dMB(0:NN,1: 2) = zero
  dMC(1: 2,1: 2) = zero

  dMD(0:NN,0:NN) = zero
  dME(0:NN,1: 2) = zero
  dMF(1: 2,1: 2) = zero
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  Wtete(1:mn,0:lrad,1:mn,0:lrad) = zero ! these "zeros" are probably not required; 10 Mar 13;
  Wteto(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Wtote(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Wtoto(1:mn,0:lrad,1:mn,0:lrad) = zero
  
  Wteze(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Wtezo(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Wtoze(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Wtozo(1:mn,0:lrad,1:mn,0:lrad) = zero
  
  Wzete(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Wzeto(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Wzote(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Wzoto(1:mn,0:lrad,1:mn,0:lrad) = zero
  
  Wzeze(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Wzezo(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Wzoze(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Wzozo(1:mn,0:lrad,1:mn,0:lrad) = zero
  
  Htete(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Hteto(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Htote(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Htoto(1:mn,0:lrad,1:mn,0:lrad) = zero
  
  Hteze(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Htezo(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Htoze(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Htozo(1:mn,0:lrad,1:mn,0:lrad) = zero
  
  Hzete(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Hzeto(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Hzote(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Hzoto(1:mn,0:lrad,1:mn,0:lrad) = zero
  
  Hzeze(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Hzezo(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Hzoze(1:mn,0:lrad,1:mn,0:lrad) = zero 
  Hzozo(1:mn,0:lrad,1:mn,0:lrad) = zero

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
   
   do jj = 1, mn ; mj = im(jj) ; nj = in(jj)
    
    mimj = mi * mj ; minj = mi * nj ; nimj = ni * mj ; ninj = ni * nj
    
    do ll = 0, lrad
     
     do pp = 0, lrad
     
      Wtete(jj,pp,ii,ll) = + 2 * ninj * TTssss(ll,pp,ii,jj) - 2 * ni      * TDszsc(ll,pp,ii,jj) - 2      * nj * TDszsc(pp,ll,jj,ii) + 2 * DDzzcc(ll,pp,ii,jj)
      Wzete(jj,pp,ii,ll) = + 2 * nimj * TTssss(ll,pp,ii,jj) + 2 * ni      * TDstsc(ll,pp,ii,jj) - 2      * mj * TDszsc(pp,ll,jj,ii) - 2 * DDtzcc(pp,ll,jj,ii)

      Htete(jj,pp,ii,ll) =   zero
      Hzete(jj,pp,ii,ll) = - DToocc(pp,ll,jj,ii) + DToocc(ll,pp,ii,jj)

      Wteze(jj,pp,ii,ll) = + 2 * minj * TTssss(ll,pp,ii,jj) + 2      * nj * TDstsc(pp,ll,jj,ii) - 2 * mi      * TDszsc(ll,pp,ii,jj) - 2 * DDtzcc(ll,pp,ii,jj)
      Wzeze(jj,pp,ii,ll) = + 2 * mimj * TTssss(ll,pp,ii,jj) + 2 * mi      * TDstsc(ll,pp,ii,jj) + 2      * mj * TDstsc(pp,ll,jj,ii) + 2 * DDttcc(ll,pp,ii,jj)

      Hteze(jj,pp,ii,ll) = - DToocc(ll,pp,ii,jj) + DToocc(pp,ll,jj,ii)  
      Hzeze(jj,pp,ii,ll) =   zero  

      if( NOTstellsym ) then

!     Wtote(jj,pp,ii,ll) = - 2 * ninj * TTsssc(pp,ll,ii,jj) + 2      * nj * TDszcc(pp,ll,jj,ii) - 2 * ni      * TDszss(ll,pp,ii,jj) + 2 * DDzzsc(ll,pp,jj,ii)
      Wtote(jj,pp,ii,ll) = - 2 * ninj * TTsscs(pp,ll,jj,ii) + 2      * nj * TDszcc(pp,ll,jj,ii) - 2 * ni      * TDszss(ll,pp,ii,jj) + 2 * DDzzcs(ll,pp,ii,jj)
!     Wzote(jj,pp,ii,ll) = - 2 * nimj * TTsssc(pp,ll,ii,jj) + 2 * ni      * TDstss(ll,pp,ii,jj) + 2      * mj * TDszcc(pp,ll,jj,ii) - 2 * DDtzsc(pp,ll,jj,ii)
      Wzote(jj,pp,ii,ll) = - 2 * nimj * TTsscs(pp,ll,jj,ii) + 2 * ni      * TDstss(ll,pp,ii,jj) + 2      * mj * TDszcc(pp,ll,jj,ii) - 2 * DDtzsc(pp,ll,jj,ii)

      Htote(jj,pp,ii,ll) =   zero
!     Hzote(jj,pp,ii,ll) = - DToosc(pp,ll,jj,ii) + DToosc(ll,pp,jj,ii)
      Hzote(jj,pp,ii,ll) = - DToosc(pp,ll,jj,ii) + DToocs(ll,pp,ii,jj)

!     Wteto(jj,pp,ii,ll) = - 2 * ninj * TTsssc(ll,pp,jj,ii) + 2 * ni      * TDszcc(ll,pp,ii,jj) - 2      * nj * TDszss(pp,ll,jj,ii) + 2 * DDzzsc(pp,ll,ii,jj)
      Wteto(jj,pp,ii,ll) = - 2 * ninj * TTsscs(ll,pp,ii,jj) + 2 * ni      * TDszcc(ll,pp,ii,jj) - 2      * nj * TDszss(pp,ll,jj,ii) + 2 * DDzzcs(pp,ll,jj,ii)
!     Wtoto(jj,pp,ii,ll) = + 2 * ninj * TTsscc(ll,pp,ii,jj) + 2 * ni      * TDszsc(ll,pp,jj,ii) + 2      * nj * TDszsc(pp,ll,ii,jj) + 2 * DDzzss(ll,pp,ii,jj)
      Wtoto(jj,pp,ii,ll) = + 2 * ninj * TTsscc(ll,pp,ii,jj) + 2 * ni      * TDszcs(ll,pp,ii,jj) + 2      * nj * TDszcs(pp,ll,jj,ii) + 2 * DDzzss(ll,pp,ii,jj)
!     Wzeto(jj,pp,ii,ll) = - 2 * nimj * TTsssc(ll,pp,jj,ii) - 2 * ni      * TDstcc(ll,pp,ii,jj) - 2      * mj * TDszss(pp,ll,jj,ii) - 2 * DDtzsc(pp,ll,ii,jj)
      Wzeto(jj,pp,ii,ll) = - 2 * nimj * TTsscs(ll,pp,ii,jj) - 2 * ni      * TDstcc(ll,pp,ii,jj) - 2      * mj * TDszss(pp,ll,jj,ii) - 2 * DDtzcs(pp,ll,jj,ii)
!     Wzoto(jj,pp,ii,ll) = + 2 * nimj * TTsscc(ll,pp,ii,jj) - 2 * ni      * TDstsc(ll,pp,jj,ii) + 2      * mj * TDszsc(pp,ll,ii,jj) - 2 * DDtzss(pp,ll,jj,ii)
      Wzoto(jj,pp,ii,ll) = + 2 * nimj * TTsscc(ll,pp,ii,jj) - 2 * ni      * TDstcs(ll,pp,ii,jj) + 2      * mj * TDszcs(pp,ll,jj,ii) - 2 * DDtzss(pp,ll,jj,ii)

      Hteto(jj,pp,ii,ll) =   zero
      Htoto(jj,pp,ii,ll) =   zero
!     Hzeto(jj,pp,ii,ll) = - DToosc(pp,ll,ii,jj) + DToosc(ll,pp,ii,jj)
      Hzeto(jj,pp,ii,ll) = - DToocs(pp,ll,jj,ii) + DToosc(ll,pp,ii,jj)
      Hzoto(jj,pp,ii,ll) = - DTooss(pp,ll,jj,ii) + DTooss(ll,pp,ii,jj)  

!     Wtoze(jj,pp,ii,ll) = - 2 * minj * TTsssc(pp,ll,ii,jj) - 2      * nj * TDstcc(pp,ll,jj,ii) - 2 * mi      * TDszss(ll,pp,ii,jj) - 2 * DDtzsc(ll,pp,jj,ii)
      Wtoze(jj,pp,ii,ll) = - 2 * minj * TTsscs(pp,ll,jj,ii) - 2      * nj * TDstcc(pp,ll,jj,ii) - 2 * mi      * TDszss(ll,pp,ii,jj) - 2 * DDtzcs(ll,pp,ii,jj)
!     Wzoze(jj,pp,ii,ll) = - 2 * mimj * TTsssc(pp,ll,ii,jj) - 2      * mj * TDstcc(pp,ll,jj,ii) + 2 * mi      * TDstss(ll,pp,ii,jj) + 2 * DDttsc(ll,pp,jj,ii)
      Wzoze(jj,pp,ii,ll) = - 2 * mimj * TTsscs(pp,ll,jj,ii) - 2      * mj * TDstcc(pp,ll,jj,ii) + 2 * mi      * TDstss(ll,pp,ii,jj) + 2 * DDttcs(ll,pp,ii,jj)

!     Htoze(jj,pp,ii,ll) = - DToosc(ll,pp,jj,ii) + DToosc(pp,ll,jj,ii)
      Htoze(jj,pp,ii,ll) = - DToocs(ll,pp,ii,jj) + DToosc(pp,ll,jj,ii)
      Hzoze(jj,pp,ii,ll) =   zero

!     Wtezo(jj,pp,ii,ll) = - 2 * minj * TTsssc(ll,pp,jj,ii) + 2      * nj * TDstss(pp,ll,jj,ii) + 2 * mi      * TDszcc(ll,pp,ii,jj) - 2 * DDtzsc(ll,pp,ii,jj)
      Wtezo(jj,pp,ii,ll) = - 2 * minj * TTsscs(ll,pp,ii,jj) + 2      * nj * TDstss(pp,ll,jj,ii) + 2 * mi      * TDszcc(ll,pp,ii,jj) - 2 * DDtzsc(ll,pp,ii,jj)
!     Wtozo(jj,pp,ii,ll) = + 2 * minj * TTsscc(ll,pp,ii,jj) - 2      * nj * TDstsc(pp,ll,ii,jj) + 2 * mi      * TDszsc(ll,pp,jj,ii) - 2 * DDtzss(ll,pp,ii,jj)
      Wtozo(jj,pp,ii,ll) = + 2 * minj * TTsscc(ll,pp,ii,jj) - 2      * nj * TDstcs(pp,ll,jj,ii) + 2 * mi      * TDszcs(ll,pp,ii,jj) - 2 * DDtzss(ll,pp,ii,jj)
!     Wzezo(jj,pp,ii,ll) = - 2 * mimj * TTsssc(ll,pp,jj,ii) - 2 * mi      * TDstcc(ll,pp,ii,jj) + 2      * mj * TDstss(pp,ll,jj,ii) + 2 * DDttsc(pp,ll,ii,jj)
      Wzezo(jj,pp,ii,ll) = - 2 * mimj * TTsscs(ll,pp,ii,jj) - 2 * mi      * TDstcc(ll,pp,ii,jj) + 2      * mj * TDstss(pp,ll,jj,ii) + 2 * DDttcs(pp,ll,jj,ii)
!     Wzozo(jj,pp,ii,ll) = + 2 * mimj * TTsscc(ll,pp,ii,jj) - 2 * mi      * TDstsc(ll,pp,jj,ii) - 2      * mj * TDstsc(pp,ll,ii,jj) + 2 * DDttss(ll,pp,ii,jj)
      Wzozo(jj,pp,ii,ll) = + 2 * mimj * TTsscc(ll,pp,ii,jj) - 2 * mi      * TDstcs(ll,pp,ii,jj) - 2      * mj * TDstcs(pp,ll,jj,ii) + 2 * DDttss(ll,pp,ii,jj)

!     Htezo(jj,pp,ii,ll) = - DToosc(ll,pp,ii,jj) + DToosc(pp,ll,ii,jj)
      Htezo(jj,pp,ii,ll) = - DToosc(ll,pp,ii,jj) + DToocs(pp,ll,jj,ii)
      Htozo(jj,pp,ii,ll) = - DTooss(ll,pp,ii,jj) + DTooss(pp,ll,jj,ii)  
      Hzezo(jj,pp,ii,ll) =   zero
      Hzozo(jj,pp,ii,ll) =   zero

      endif ! end of if( NOTstellsym ) 26 Feb 13;

#ifdef DEBUG
      if( Wma01ag ) then
       write(ounit,'("ma01ag : "5i3" : Wtete="4es23.15)') lvol, ii,jj,ll,pp, Wtete(jj,pp,ii,ll), Wzete(jj,pp,ii,ll), Wteze(jj,pp,ii,ll), Wzeze(jj,pp,ii,ll)
       write(ounit,'("ma01ag : "5i3" : Wteto="4es23.15)') lvol, ii,jj,ll,pp, Wteto(jj,pp,ii,ll), Wzeto(jj,pp,ii,ll), Wtezo(jj,pp,ii,ll), Wzezo(jj,pp,ii,ll)
       write(ounit,'("ma01ag : "5i3" : Wtote="4es23.15)') lvol, ii,jj,ll,pp, Wtote(jj,pp,ii,ll), Wzote(jj,pp,ii,ll), Wtoze(jj,pp,ii,ll), Wzoze(jj,pp,ii,ll)
       write(ounit,'("ma01ag : "5i3" : Wtoto="4es23.15)') lvol, ii,jj,ll,pp, Wtoto(jj,pp,ii,ll), Wzoto(jj,pp,ii,ll), Wtozo(jj,pp,ii,ll), Wzozo(jj,pp,ii,ll)
       write(ounit,'("ma01ag : "5i3" : Htete="4es23.15)') lvol, ii,jj,ll,pp, Htete(jj,pp,ii,ll), Hzete(jj,pp,ii,ll), Hteze(jj,pp,ii,ll), Hzeze(jj,pp,ii,ll)
       write(ounit,'("ma01ag : "5i3" : Hteto="4es23.15)') lvol, ii,jj,ll,pp, Hteto(jj,pp,ii,ll), Hzeto(jj,pp,ii,ll), Htezo(jj,pp,ii,ll), Hzezo(jj,pp,ii,ll)
       write(ounit,'("ma01ag : "5i3" : Htote="4es23.15)') lvol, ii,jj,ll,pp, Htote(jj,pp,ii,ll), Hzote(jj,pp,ii,ll), Htoze(jj,pp,ii,ll), Hzoze(jj,pp,ii,ll)
       write(ounit,'("ma01ag : "5i3" : Htoto="4es23.15)') lvol, ii,jj,ll,pp, Htoto(jj,pp,ii,ll), Hzoto(jj,pp,ii,ll), Htozo(jj,pp,ii,ll), Hzozo(jj,pp,ii,ll)
      endif
#endif

     enddo ! end of do pp;
     
    enddo ! end of do jj;
    
   enddo ! end of do ll;
   
  enddo ! end of do ii; 

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! construct matrix elements;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
   
   if( Lcoordinatesingularity .and. ( mi.ne.0 .or. ii.eq.1 ) ) then ; Lzi = .true. ! 15 Jan 15; ! additional freedom; 15 Jan 15;
   else                                                             ; Lzi = .false.
   endif
   
   do ll = 0, lrad-1 ! this is a derivative loop;
    
    if( (ll/2)*2.eq.ll ) then ; Lul = lrad   ! ll is even; 18 Jan 13;
    else                      ; Lul = lrad-1 ! ll is odd ; 18 Jan 13;
    endif

    if( Lzi ) then ; Lzl = lrad
    else           ; Lzl = Lul
    endif
    
    do jj = 1, mn ; mj = im(jj) ; nj = in(jj)
     
     if( Lcoordinatesingularity .and. ( mj.ne.0 .or. nj.eq.0 ) ) then ; Lzj = .true.
     else                                                             ; Lzj = .false.
     endif
     
     do pp = 0, lrad-1 ! this is a derivative loop;
      
      if( (pp/2)*2.eq.pp ) then ; Lup = lrad   ! pp is even; 18 Jan 13;
      else                      ; Lup = lrad-1 ! pp is odd ; 18 Jan 13;
      endif
      
      if( Lzj ) then ; Lzp = lrad
      else           ; Lzp = Lup
      endif
      
#ifdef DEBUG
      FATAL( ma01ag, Ate(lvol,0,ii)%i(ll).lt.0 .or. Ate(lvol,0,ii)%i(ll).gt.NN, invalid subscript )
      FATAL( ma01ag, Ate(lvol,0,jj)%i(pp).lt.0 .or. Ate(lvol,0,jj)%i(pp).gt.NN, invalid subscript )
      FATAL( ma01ag, Aze(lvol,0,ii)%i(ll).lt.0 .or. Aze(lvol,0,ii)%i(ll).gt.NN, invalid subscript )
      FATAL( ma01ag, Aze(lvol,0,jj)%i(pp).lt.0 .or. Aze(lvol,0,jj)%i(pp).gt.NN, invalid subscript )
      if( NOTstellsym ) then
      FATAL( ma01ag, Ato(lvol,0,ii)%i(ll).lt.0 .or. Ato(lvol,0,ii)%i(ll).gt.NN, invalid subscript )
      FATAL( ma01ag, Ato(lvol,0,jj)%i(pp).lt.0 .or. Ato(lvol,0,jj)%i(pp).gt.NN, invalid subscript )
      FATAL( ma01ag, Azo(lvol,0,ii)%i(ll).lt.0 .or. Azo(lvol,0,ii)%i(ll).gt.NN, invalid subscript )
      FATAL( ma01ag, Azo(lvol,0,jj)%i(pp).lt.0 .or. Azo(lvol,0,jj)%i(pp).gt.NN, invalid subscript )
      endif
#endif

      id=Ate(lvol,0,ii)%i(ll) ; jd=Ate(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + Wtete(jj,pp,ii,ll) - Wtete(jj,pp,ii,Lul) - Wtete(jj,Lup,ii,ll) + Wtete(jj,Lup,ii,Lul)

      id=Ate(lvol,0,ii)%i(ll) ; jd=Aze(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + Wzete(jj,pp,ii,ll) - Wzete(jj,pp,ii,Lul) - Wzete(jj,Lzp,ii,ll) + Wzete(jj,Lzp,ii,Lul)

      id=Aze(lvol,0,ii)%i(ll) ; jd=Ate(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + Wteze(jj,pp,ii,ll) - Wteze(jj,pp,ii,Lzl) - Wteze(jj,Lup,ii,ll) + Wteze(jj,Lup,ii,Lzl)

      id=Aze(lvol,0,ii)%i(ll) ; jd=Aze(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + Wzeze(jj,pp,ii,ll) - Wzeze(jj,pp,ii,Lzl) - Wzeze(jj,Lzp,ii,ll) + Wzeze(jj,Lzp,ii,Lzl)

      id=Ate(lvol,0,ii)%i(ll) ; jd=Ate(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + Htete(jj,pp,ii,ll) - Htete(jj,pp,ii,Lul) - Htete(jj,Lup,ii,ll) + Htete(jj,Lup,ii,Lul)

      id=Ate(lvol,0,ii)%i(ll) ; jd=Aze(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + Hzete(jj,pp,ii,ll) - Hzete(jj,pp,ii,Lul) - Hzete(jj,Lzp,ii,ll) + Hzete(jj,Lzp,ii,Lul)

      id=Aze(lvol,0,ii)%i(ll) ; jd=Ate(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + Hteze(jj,pp,ii,ll) - Hteze(jj,pp,ii,Lzl) - Hteze(jj,Lup,ii,ll) + Hteze(jj,Lup,ii,Lzl)

      id=Aze(lvol,0,ii)%i(ll) ; jd=Aze(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + Hzeze(jj,pp,ii,ll) - Hzeze(jj,pp,ii,Lzl) - Hzeze(jj,Lzp,ii,ll) + Hzeze(jj,Lzp,ii,Lzl)

      if( NOTstellsym ) then

      id=Ate(lvol,0,ii)%i(ll) ; jd=Ato(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + Wtote(jj,pp,ii,ll) - Wtote(jj,pp,ii,Lul) - Wtote(jj,Lup,ii,ll) + Wtote(jj,Lup,ii,Lul)

      id=Ate(lvol,0,ii)%i(ll) ; jd=Azo(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + Wzote(jj,pp,ii,ll) - Wzote(jj,pp,ii,Lul) - Wzote(jj,Lzp,ii,ll) + Wzote(jj,Lzp,ii,Lul)

      id=Ato(lvol,0,ii)%i(ll) ; jd=Ate(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + Wteto(jj,pp,ii,ll) - Wteto(jj,pp,ii,Lul) - Wteto(jj,Lup,ii,ll) + Wteto(jj,Lup,ii,Lul)

      id=Ato(lvol,0,ii)%i(ll) ; jd=Ato(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + Wtoto(jj,pp,ii,ll) - Wtoto(jj,pp,ii,Lul) - Wtoto(jj,Lup,ii,ll) + Wtoto(jj,Lup,ii,Lul)

      id=Ato(lvol,0,ii)%i(ll) ; jd=Aze(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + Wzeto(jj,pp,ii,ll) - Wzeto(jj,pp,ii,Lul) - Wzeto(jj,Lzp,ii,ll) + Wzeto(jj,Lzp,ii,Lul)

      id=Ato(lvol,0,ii)%i(ll) ; jd=Azo(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + Wzoto(jj,pp,ii,ll) - Wzoto(jj,pp,ii,Lul) - Wzoto(jj,Lzp,ii,ll) + Wzoto(jj,Lzp,ii,Lul)

      id=Aze(lvol,0,ii)%i(ll) ; jd=Ato(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + Wtoze(jj,pp,ii,ll) - Wtoze(jj,pp,ii,Lzl) - Wtoze(jj,Lup,ii,ll) + Wtoze(jj,Lup,ii,Lzl)

      id=Aze(lvol,0,ii)%i(ll) ; jd=Azo(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + Wzoze(jj,pp,ii,ll) - Wzoze(jj,pp,ii,Lzl) - Wzoze(jj,Lzp,ii,ll) + Wzoze(jj,Lzp,ii,Lzl)

      id=Azo(lvol,0,ii)%i(ll) ; jd=Ate(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + Wtezo(jj,pp,ii,ll) - Wtezo(jj,pp,ii,Lzl) - Wtezo(jj,Lup,ii,ll) + Wtezo(jj,Lup,ii,Lzl)

      id=Azo(lvol,0,ii)%i(ll) ; jd=Ato(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + Wtozo(jj,pp,ii,ll) - Wtozo(jj,pp,ii,Lzl) - Wtozo(jj,Lup,ii,ll) + Wtozo(jj,Lup,ii,Lzl)

      id=Azo(lvol,0,ii)%i(ll) ; jd=Aze(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + Wzezo(jj,pp,ii,ll) - Wzezo(jj,pp,ii,Lzl) - Wzezo(jj,Lzp,ii,ll) + Wzezo(jj,Lzp,ii,Lzl)

      id=Azo(lvol,0,ii)%i(ll) ; jd=Azo(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + Wzozo(jj,pp,ii,ll) - Wzozo(jj,pp,ii,Lzl) - Wzozo(jj,Lzp,ii,ll) + Wzozo(jj,Lzp,ii,Lzl)

      id=Ate(lvol,0,ii)%i(ll) ; jd=Ato(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + Htote(jj,pp,ii,ll) - Htote(jj,pp,ii,Lul) - Htote(jj,Lup,ii,ll) + Htote(jj,Lup,ii,Lul)

      id=Ate(lvol,0,ii)%i(ll) ; jd=Azo(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + Hzote(jj,pp,ii,ll) - Hzote(jj,pp,ii,Lul) - Hzote(jj,Lzp,ii,ll) + Hzote(jj,Lzp,ii,Lul)

      id=Ato(lvol,0,ii)%i(ll) ; jd=Ate(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + Hteto(jj,pp,ii,ll) - Hteto(jj,pp,ii,Lul) - Hteto(jj,Lup,ii,ll) + Hteto(jj,Lup,ii,Lul)

      id=Ato(lvol,0,ii)%i(ll) ; jd=Ato(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + Htoto(jj,pp,ii,ll) - Htoto(jj,pp,ii,Lul) - Htoto(jj,Lup,ii,ll) + Htoto(jj,Lup,ii,Lul)

      id=Ato(lvol,0,ii)%i(ll) ; jd=Aze(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + Hzeto(jj,pp,ii,ll) - Hzeto(jj,pp,ii,Lul) - Hzeto(jj,Lzp,ii,ll) + Hzeto(jj,Lzp,ii,Lul)

      id=Ato(lvol,0,ii)%i(ll) ; jd=Azo(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + Hzoto(jj,pp,ii,ll) - Hzoto(jj,pp,ii,Lul) - Hzoto(jj,Lzp,ii,ll) + Hzoto(jj,Lzp,ii,Lul)

      id=Aze(lvol,0,ii)%i(ll) ; jd=Ato(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + Htoze(jj,pp,ii,ll) - Htoze(jj,pp,ii,Lzl) - Htoze(jj,Lup,ii,ll) + Htoze(jj,Lup,ii,Lzl)

      id=Aze(lvol,0,ii)%i(ll) ; jd=Azo(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + Hzoze(jj,pp,ii,ll) - Hzoze(jj,pp,ii,Lzl) - Hzoze(jj,Lzp,ii,ll) + Hzoze(jj,Lzp,ii,Lzl)

      id=Azo(lvol,0,ii)%i(ll) ; jd=Ate(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + Htezo(jj,pp,ii,ll) - Htezo(jj,pp,ii,Lzl) - Htezo(jj,Lup,ii,ll) + Htezo(jj,Lup,ii,Lzl)

      id=Azo(lvol,0,ii)%i(ll) ; jd=Ato(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + Htozo(jj,pp,ii,ll) - Htozo(jj,pp,ii,Lzl) - Htozo(jj,Lup,ii,ll) + Htozo(jj,Lup,ii,Lzl)

      id=Azo(lvol,0,ii)%i(ll) ; jd=Aze(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + Hzezo(jj,pp,ii,ll) - Hzezo(jj,pp,ii,Lzl) - Hzezo(jj,Lzp,ii,ll) + Hzezo(jj,Lzp,ii,Lzl)

      id=Azo(lvol,0,ii)%i(ll) ; jd=Azo(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + Hzozo(jj,pp,ii,ll) - Hzozo(jj,pp,ii,Lzl) - Hzozo(jj,Lzp,ii,ll) + Hzozo(jj,Lzp,ii,Lzl)

      endif

     enddo ! end of do pp; 24 Jan 13;

    enddo ! end of do jj; 24 Jan 13;
    
   enddo ! end of do ll; 24 Jan 13;

  enddo ! end of do ii; 24 Jan 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!#ifdef DEBUG
!  if( Wma01ag ) write(ounit,'("ma01ag : ", 10x ," : myid=",i3," ; lvol=",i3," ; 1.0 ;")') myid, lvol
!#endif

  do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
   
  !if( Lcoordinatesingularity .and. ( mi.ne.0 .or. ni.eq.0 ) ) then ; Lzi = .true. ! 15 Jan 15; ! additional freedom; 15 Jan 15;
   if( Lcoordinatesingularity .and. ( mi.ne.0 .or. ii.eq.1 ) ) then ; Lzi = .true. ! 15 Jan 15; ! additional freedom; 15 Jan 15;
   else                                                             ; Lzi = .false.
   endif
   
   do ll = 0, lrad-1 ! this is a derivative loop;
    
    if( (ll/2)*2.eq.ll ) then ; Lul = lrad   ! ll is even; 18 Jan 13;
    else                      ; Lul = lrad-1 ! ll is odd ; 18 Jan 13;
    endif
    
    if( Lzi ) then ; Lzl = lrad
    else           ; Lzl = Lul
    endif
    
    do jj = 1, mn ; mj = im(jj) ; nj = in(jj)
   
     if( Lcoordinatesingularity .and. ( mj.ne.0 .or. nj.eq.0 ) ) then ; Lzj = .true.
     else                                                             ; Lzj = .false.
     endif
     
     do pp = lrad-1, lrad ! this is a summation loop; 18 Feb 13;

      if( Lzj ) then ; Lzp = lrad
      else           ; Lzp = pp
      endif
      
!#ifdef DEBUG
!      FATAL( ma01ag, Ate(lvol,0,ii)%i(ll).lt.0 .or. Ate(lvol,0,ii)%i(ll).gt.NN, invalid subscript )
!      FATAL( ma01ag, Aze(lvol,0,ii)%i(ll).lt.0 .or. Aze(lvol,0,ii)%i(ll).gt.NN, invalid subscript )
!      FATAL( ma01ag, Fso(lvol,jj)        .lt.0 .or. Fso(lvol,jj)        .gt.NN, invalid subscript )
!      if( NOTstellsym ) then
!      FATAL( ma01ag, Ato(lvol,0,ii)%i(ll).lt.0 .or. Ato(lvol,0,ii)%i(ll).gt.NN, invalid subscript )
!      FATAL( ma01ag, Azo(lvol,0,ii)%i(ll).lt.0 .or. Azo(lvol,0,ii)%i(ll).gt.NN, invalid subscript )
!      FATAL( ma01ag, Fse(lvol,jj)        .lt.0 .or. Fse(lvol,jj)        .gt.NN, invalid subscript )
!      endif
!#endif
      
      id=Ate(lvol,0,ii)%i(ll) ; jd=Fso(lvol,jj)
      dMA(id,jd) = dMA(id,jd) + half * ( + mj*(Wtete(jj,pp,ii,ll)-Wtete(jj,pp,ii,Lul)) - nj*(Wzete(jj,Lzp,ii,ll)-Wzete(jj,Lzp,ii,Lul)) )

      id=Aze(lvol,0,ii)%i(ll) ; jd=Fso(lvol,jj)
      dMA(id,jd) = dMA(id,jd) + half * ( + mj*(Wteze(jj,pp,ii,ll)-Wteze(jj,pp,ii,Lzl)) - nj*(Wzeze(jj,Lzp,ii,ll)-Wzeze(jj,Lzp,ii,Lzl)) )

      id=Ate(lvol,0,ii)%i(ll) ; jd=Fso(lvol,jj)
      dMD(id,jd) = dMD(id,jd) + half * ( + mj*(Htete(jj,pp,ii,ll)-Htete(jj,pp,ii,Lul)) - nj*(Hzete(jj,Lzp,ii,ll)-Hzete(jj,Lzp,ii,Lul)) )

      id=Aze(lvol,0,ii)%i(ll) ; jd=Fso(lvol,jj)
      dMD(id,jd) = dMD(id,jd) + half * ( + mj*(Hteze(jj,pp,ii,ll)-Hteze(jj,pp,ii,Lzl)) - nj*(Hzeze(jj,Lzp,ii,ll)-Hzeze(jj,Lzp,ii,Lzl)) )

      if( NOTstellsym ) then
      id=Ate(lvol,0,ii)%i(ll) ; jd=Fse(lvol,jj)
      dMA(id,jd) = dMA(id,jd) + half * ( - mj*(Wtote(jj,pp,ii,ll)-Wtote(jj,pp,ii,Lul)) + nj*(Wzote(jj,Lzp,ii,ll)-Wzote(jj,Lzp,ii,Lul)) )

      id=Ato(lvol,0,ii)%i(ll) ; jd=Fso(lvol,jj)
      dMA(id,jd) = dMA(id,jd) + half * ( + mj*(Wteto(jj,pp,ii,ll)-Wteto(jj,pp,ii,Lul)) - nj*(Wzeto(jj,Lzp,ii,ll)-Wzeto(jj,Lzp,ii,Lul)) )

      id=Ato(lvol,0,ii)%i(ll) ; jd=Fse(lvol,jj)
      dMA(id,jd) = dMA(id,jd) + half * ( - mj*(Wtoto(jj,pp,ii,ll)-Wtoto(jj,pp,ii,Lul)) + nj*(Wzoto(jj,Lzp,ii,ll)-Wzoto(jj,Lzp,ii,Lul)) )

      id=Aze(lvol,0,ii)%i(ll) ; jd=Fse(lvol,jj)
      dMA(id,jd) = dMA(id,jd) + half * ( - mj*(Wtoze(jj,pp,ii,ll)-Wtoze(jj,pp,ii,Lzl)) + nj*(Wzoze(jj,Lzp,ii,ll)-Wzoze(jj,Lzp,ii,Lzl)) )

      id=Azo(lvol,0,ii)%i(ll) ; jd=Fso(lvol,jj)
      dMA(id,jd) = dMA(id,jd) + half * ( + mj*(Wtezo(jj,pp,ii,ll)-Wtezo(jj,pp,ii,Lzl)) - nj*(Wzezo(jj,Lzp,ii,ll)-Wzezo(jj,Lzp,ii,Lzl)) )

      id=Azo(lvol,0,ii)%i(ll) ; jd=Fse(lvol,jj)
      dMA(id,jd) = dMA(id,jd) + half * ( - mj*(Wtozo(jj,pp,ii,ll)-Wtozo(jj,pp,ii,Lzl)) + nj*(Wzozo(jj,Lzp,ii,ll)-Wzozo(jj,Lzp,ii,Lzl)) )

      id=Ate(lvol,0,ii)%i(ll) ; jd=Fse(lvol,jj)
      dMD(id,jd) = dMD(id,jd) + half * ( - mj*(Htote(jj,pp,ii,ll)-Htote(jj,pp,ii,Lul)) + nj*(Hzote(jj,Lzp,ii,ll)-Hzote(jj,Lzp,ii,Lul)) )

      id=Ato(lvol,0,ii)%i(ll) ; jd=Fso(lvol,jj)
      dMD(id,jd) = dMD(id,jd) + half * ( + mj*(Hteto(jj,pp,ii,ll)-Hteto(jj,pp,ii,Lul)) - nj*(Hzeto(jj,Lzp,ii,ll)-Hzeto(jj,Lzp,ii,Lul)) )

      id=Ato(lvol,0,ii)%i(ll) ; jd=Fse(lvol,jj)
      dMD(id,jd) = dMD(id,jd) + half * ( - mj*(Htoto(jj,pp,ii,ll)-Htoto(jj,pp,ii,Lul)) + nj*(Hzoto(jj,Lzp,ii,ll)-Hzoto(jj,Lzp,ii,Lul)) )

      id=Aze(lvol,0,ii)%i(ll) ; jd=Fse(lvol,jj)
      dMD(id,jd) = dMD(id,jd) + half * ( - mj*(Htoze(jj,pp,ii,ll)-Htoze(jj,pp,ii,Lzl)) + nj*(Hzoze(jj,Lzp,ii,ll)-Hzoze(jj,Lzp,ii,Lzl)) )

      id=Azo(lvol,0,ii)%i(ll) ; jd=Fso(lvol,jj)
      dMD(id,jd) = dMD(id,jd) + half * ( + mj*(Htezo(jj,pp,ii,ll)-Htezo(jj,pp,ii,Lzl)) - nj*(Hzezo(jj,Lzp,ii,ll)-Hzezo(jj,Lzp,ii,Lzl)) )

      id=Azo(lvol,0,ii)%i(ll) ; jd=Fse(lvol,jj)
      dMD(id,jd) = dMD(id,jd) + half * ( - mj*(Htozo(jj,pp,ii,ll)-Htozo(jj,pp,ii,Lzl)) + nj*(Hzozo(jj,Lzp,ii,ll)-Hzozo(jj,Lzp,ii,Lzl)) )

      endif

     enddo ! end of do pp; 24 Jan 13;
     
    enddo ! end of do jj; 24 Jan 13;
    
   enddo ! end of do ll; 24 Jan 13;

  enddo ! end of do ii; 24 Jan 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
!#ifdef DEBUG
!  if( Wma01ag ) write(ounit,'("ma01ag : ", 10x ," : myid=",i3," ; lvol=",i3," ; 2.0 ;")') myid, lvol
!#endif

  do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
   
  !if( Lcoordinatesingularity .and. ( mi.ne.0 .or. ni.eq.0 ) ) then ; Lzi = .true. ! 15 Jan 15; ! additional freedom; 15 Jan 15;
   if( Lcoordinatesingularity .and. ( mi.ne.0 .or. ii.eq.1 ) ) then ; Lzi = .true. ! 15 Jan 15; ! additional freedom; 15 Jan 15;
   else                                                             ; Lzi = .false.
   endif
   
   do ll = lrad-1, lrad ! this is a summation loop; 18 Feb 13;
    
    if( Lzi ) then ; Lzl = lrad
    else           ; Lzl = ll
    endif
    
    do jj = 1, mn ; mj = im(jj) ; nj = in(jj)
   
     if( Lcoordinatesingularity .and. ( mj.ne.0 .or. nj.eq.0 ) ) then ; Lzj = .true.
     else                                                             ; Lzj = .false.
     endif
     
     do pp = 0, lrad-1 ! this is a derivative loop;
      
      if( (pp/2)*2.eq.pp ) then ; Lup = lrad   ! ll is even; 18 Jan 13;
      else                      ; Lup = lrad-1 ! ll is odd ; 18 Jan 13;
      endif

      if( Lzj ) then ; Lzp = lrad
      else           ; Lzp = Lup
      endif
      
!#ifdef DEBUG
!      FATAL( ma01ag, Ate(lvol,0,jj)%i(pp).lt.0 .or. Ate(lvol,0,jj)%i(pp).gt.NN, invalid subscript )
!      FATAL( ma01ag, Aze(lvol,0,jj)%i(pp).lt.0 .or. Aze(lvol,0,jj)%i(pp).gt.NN, invalid subscript )
!      FATAL( ma01ag, Fso(lvol,ii)        .lt.0 .or. Fso(lvol,ii)        .gt.NN, invalid subscript )
!      if( NOTstellsym ) then
!      FATAL( ma01ag, Ato(lvol,0,jj)%i(pp).lt.0 .or. Ato(lvol,0,jj)%i(pp).gt.NN, invalid subscript )
!      FATAL( ma01ag, Azo(lvol,0,jj)%i(pp).lt.0 .or. Azo(lvol,0,jj)%i(pp).gt.NN, invalid subscript )
!      FATAL( ma01ag, Fse(lvol,ii)        .lt.0 .or. Fse(lvol,ii)        .gt.NN, invalid subscript )
!      endif
!#endif

      id=Fso(lvol,ii) ; jd=Ate(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + half * ( + mi*(Wtete(jj,pp,ii,ll)-Wtete(jj,Lup,ii,ll)) - ni*(Wteze(jj,pp,ii,Lzl)-Wteze(jj,Lup,ii,Lzl)) )

      id=Fso(lvol,ii) ; jd=Aze(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + half * ( + mi*(Wzete(jj,pp,ii,ll)-Wzete(jj,Lzp,ii,ll)) - ni*(Wzeze(jj,pp,ii,Lzl)-Wzeze(jj,Lzp,ii,Lzl)) )

      id=Fso(lvol,ii) ; jd=Ate(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + half * ( + mi*(Htete(jj,pp,ii,ll)-Htete(jj,Lup,ii,ll)) - ni*(Hteze(jj,pp,ii,Lzl)-Hteze(jj,Lup,ii,Lzl)) )

      id=Fso(lvol,ii) ; jd=Aze(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + half * ( + mi*(Hzete(jj,pp,ii,ll)-Hzete(jj,Lzp,ii,ll)) - ni*(Hzeze(jj,pp,ii,Lzl)-Hzeze(jj,Lzp,ii,Lzl)) )

      if( NOTstellsym ) then

      id=Fso(lvol,ii) ; jd=Ato(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + half * ( + mi*(Wtote(jj,pp,ii,ll)-Wtote(jj,Lup,ii,ll)) - ni*(Wtoze(jj,pp,ii,Lzl)-Wtoze(jj,Lup,ii,Lzl)) )

      id=Fso(lvol,ii) ; jd=Azo(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + half * ( + mi*(Wzote(jj,pp,ii,ll)-Wzote(jj,Lzp,ii,ll)) - ni*(Wzoze(jj,pp,ii,Lzl)-Wzoze(jj,Lzp,ii,Lzl)) )

      id=Fse(lvol,ii) ; jd=Ate(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + half * ( - mi*(Wteto(jj,pp,ii,ll)-Wteto(jj,Lup,ii,ll)) + ni*(Wtezo(jj,pp,ii,Lzl)-Wtezo(jj,Lup,ii,Lzl)) )

      id=Fse(lvol,ii) ; jd=Ato(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + half * ( - mi*(Wtoto(jj,pp,ii,ll)-Wtoto(jj,Lup,ii,ll)) + ni*(Wtozo(jj,pp,ii,Lzl)-Wtozo(jj,Lup,ii,Lzl)) )

      id=Fse(lvol,ii) ; jd=Aze(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + half * ( - mi*(Wzeto(jj,pp,ii,ll)-Wzeto(jj,Lzp,ii,ll)) + ni*(Wzezo(jj,pp,ii,Lzl)-Wzezo(jj,Lzp,ii,Lzl)) )

      id=Fse(lvol,ii) ; jd=Azo(lvol,0,jj)%i(pp)
      dMA(id,jd) = dMA(id,jd) + half * ( - mi*(Wzoto(jj,pp,ii,ll)-Wzoto(jj,Lzp,ii,ll)) + ni*(Wzozo(jj,pp,ii,Lzl)-Wzozo(jj,Lzp,ii,Lzl)) )

      id=Fso(lvol,ii) ; jd=Ato(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + half * ( + mi*(Htote(jj,pp,ii,ll)-Htote(jj,Lup,ii,ll)) - ni*(Htoze(jj,pp,ii,Lzl)-Htoze(jj,Lup,ii,Lzl)) )

      id=Fso(lvol,ii) ; jd=Azo(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + half * ( + mi*(Hzote(jj,pp,ii,ll)-Hzote(jj,Lzp,ii,ll)) - ni*(Hzoze(jj,pp,ii,Lzl)-Hzoze(jj,Lzp,ii,Lzl)) )

      id=Fse(lvol,ii) ; jd=Ate(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + half * ( - mi*(Hteto(jj,pp,ii,ll)-Hteto(jj,Lup,ii,ll)) + ni*(Htezo(jj,pp,ii,Lzl)-Htezo(jj,Lup,ii,Lzl)) )

      id=Fse(lvol,ii) ; jd=Ato(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + half * ( - mi*(Htoto(jj,pp,ii,ll)-Htoto(jj,Lup,ii,ll)) + ni*(Htozo(jj,pp,ii,Lzl)-Htozo(jj,Lup,ii,Lzl)) )

      id=Fse(lvol,ii) ; jd=Aze(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + half * ( - mi*(Hzeto(jj,pp,ii,ll)-Hzeto(jj,Lzp,ii,ll)) + ni*(Hzezo(jj,pp,ii,Lzl)-Hzezo(jj,Lzp,ii,Lzl)) )

      id=Fse(lvol,ii) ; jd=Azo(lvol,0,jj)%i(pp)
      dMD(id,jd) = dMD(id,jd) + half * ( - mi*(Hzoto(jj,pp,ii,ll)-Hzoto(jj,Lzp,ii,ll)) + ni*(Hzozo(jj,pp,ii,Lzl)-Hzozo(jj,Lzp,ii,Lzl)) )

      endif

     enddo ! end of do pp; 24 Jan 13;

    enddo ! end of do jj; 24 Jan 13;
    
   enddo ! end of do ll; 24 Jan 13;

  enddo ! end of do ii; 24 Jan 13;
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!#ifdef DEBUG
!  if( Wma01ag ) write(ounit,'("ma01ag : ", 10x ," : myid=",i3," ; lvol=",i3," ; 3.0 ;")') myid, lvol
!#endif

  do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
   
  !if( Lcoordinatesingularity .and. ( mi.ne.0 .or. ni.eq.0 ) ) then ; Lzi = .true. ! 15 Jan 15; ! additional freedom; 15 Jan 15;
   if( Lcoordinatesingularity .and. ( mi.ne.0 .or. ii.eq.1 ) ) then ; Lzi = .true. ! 15 Jan 15; ! additional freedom; 15 Jan 15;
   else                                                             ; Lzi = .false.
   endif
   
   do ll = lrad-1, lrad ! this is a summation loop; 18 Feb 13;
    
    if( Lzi ) then ; Lzl = lrad
    else           ; Lzl = ll
    endif
    
    do jj = 1, mn ; mj = im(jj) ; nj = in(jj)
   
     if( Lcoordinatesingularity .and. ( mj.ne.0 .or. nj.eq.0 ) ) then ; Lzj = .true.
     else                                                             ; Lzj = .false.
     endif
     
     mjmi = mj * mi ; mjni = mj * ni ; njmi = nj * mi ; njni = nj * ni
     
     do pp = lrad-1, lrad ! this is a summation loop; 18 Feb 13;
    
      if( Lzj ) then ; Lzp = lrad
      else           ; Lzp = pp
      endif

!#ifdef DEBUG
!      FATAL( ma01ag, Fso(lvol,ii)        .lt.0 .or. Fso(lvol,ii)        .gt.NN, invalid subscript )
!      FATAL( ma01ag, Fso(lvol,jj)        .lt.0 .or. Fso(lvol,jj)        .gt.NN, invalid subscript )
!      if( NOTstellsym ) then
!      FATAL( ma01ag, Fse(lvol,ii)        .lt.0 .or. Fse(lvol,ii)        .gt.NN, invalid subscript )
!      FATAL( ma01ag, Fse(lvol,jj)        .lt.0 .or. Fse(lvol,jj)        .gt.NN, invalid subscript )
!      endif
!#endif

      id=Fso(lvol,ii) ; jd=Fso(lvol,jj)
      dMA(id,jd) = dMA(id,jd) + quart*( +mjmi*Wtete(jj,pp,ii,ll)-mjni*Wteze(jj,pp,ii,Lzl)-njmi*Wzete(jj,Lzp,ii,ll)+njni*Wzeze(jj,Lzp,ii,Lzl))

      id=Fso(lvol,ii) ; jd=Fso(lvol,jj)
      dMD(id,jd) = dMD(id,jd) + quart*( +mjmi*Htete(jj,pp,ii,ll)-mjni*Hteze(jj,pp,ii,Lzl)-njmi*Hzete(jj,Lzp,ii,ll)+njni*Hzeze(jj,Lzp,ii,Lzl))

      if( NOTstellsym ) then

      id=Fso(lvol,ii) ; jd=Fse(lvol,jj)
      dMA(id,jd) = dMA(id,jd) + quart*( -mjmi*Wtote(jj,pp,ii,ll)+mjni*Wtoze(jj,pp,ii,Lzl)+njmi*Wzote(jj,Lzp,ii,ll)-njni*Wzoze(jj,Lzp,ii,Lzl))

      id=Fse(lvol,ii) ; jd=Fso(lvol,jj)
      dMA(id,jd) = dMA(id,jd) + quart*( -mjmi*Wteto(jj,pp,ii,ll)+mjni*Wtezo(jj,pp,ii,Lzl)+njmi*Wzeto(jj,Lzp,ii,ll)-njni*Wzezo(jj,Lzp,ii,Lzl))

      id=Fse(lvol,ii) ; jd=Fse(lvol,jj)
      dMA(id,jd) = dMA(id,jd) + quart*( +mjmi*Wtoto(jj,pp,ii,ll)-mjni*Wtozo(jj,pp,ii,Lzl)-njmi*Wzoto(jj,Lzp,ii,ll)+njni*Wzozo(jj,Lzp,ii,Lzl))

      id=Fso(lvol,ii) ; jd=Fse(lvol,jj)
      dMD(id,jd) = dMD(id,jd) + quart*( -mjmi*Htote(jj,pp,ii,ll)+mjni*Htoze(jj,pp,ii,Lzl)+njmi*Hzote(jj,Lzp,ii,ll)-njni*Hzoze(jj,Lzp,ii,Lzl))

      id=Fse(lvol,ii) ; jd=Fso(lvol,jj)
      dMD(id,jd) = dMD(id,jd) + quart*( -mjmi*Hteto(jj,pp,ii,ll)+mjni*Htezo(jj,pp,ii,Lzl)+njmi*Hzeto(jj,Lzp,ii,ll)-njni*Hzezo(jj,Lzp,ii,Lzl))

      id=Fse(lvol,ii) ; jd=Fse(lvol,jj)
      dMD(id,jd) = dMD(id,jd) + quart*( +mjmi*Htoto(jj,pp,ii,ll)-mjni*Htozo(jj,pp,ii,Lzl)-njmi*Hzoto(jj,Lzp,ii,ll)+njni*Hzozo(jj,Lzp,ii,Lzl))

      endif

     enddo ! end of do pp; 24 Jan 13;
     
    enddo ! end of do jj; 24 Jan 13;
    
   enddo ! end of do ll; 24 Jan 13;
   
  enddo ! end of do ii; 24 Jan 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!#ifdef DEBUG
!  if( Wma01ag ) write(ounit,'("ma01ag : ", 10x ," : myid=",i3," ; lvol=",i3," ; 4.0 ;")') myid, lvol
!#endif

  do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
   
  !if( Lcoordinatesingularity .and. ( mi.ne.0 .or. ni.eq.0 ) ) then ; Lzi = .true. ! 15 Jan 15; ! additional freedom; 15 Jan 15;
   if( Lcoordinatesingularity .and. ( mi.ne.0 .or. ii.eq.1 ) ) then ; Lzi = .true. ! 15 Jan 15; ! additional freedom; 15 Jan 15;
   else                                                             ; Lzi = .false.
   endif
   
   do ll = 0, lrad-1 ! this is a derivative loop;

    if( (ll/2)*2.eq.ll ) then ; Lul = lrad   ! ll is even; 18 Jan 13;
    else                      ; Lul = lrad-1 ! ll is odd ; 18 Jan 13;
    endif
    
    if( Lzi ) then ; Lzl = lrad
    else           ; Lzl = Lul
    endif
    
    jj = 1 ; mj = im(jj) ; nj = in(jj)
     
     if( Lcoordinatesingularity .and. ( mj.ne.0 .or. nj.eq.0 ) ) then ; Lzj = .true.
     else                                                             ; Lzj = .false.
     endif
     
     do pp = lrad-1, lrad ! this is a summation loop; 18 Feb 13;
      
      if( Lzj ) then ; Lzp = lrad
      else           ; Lzp = pp
      endif
      
!#ifdef DEBUG
!      FATAL( ma01ag, Ate(lvol,0,ii)%i(ll).lt.0 .or. Ate(lvol,0,ii)%i(ll).gt.NN, invalid subscript )
!      FATAL( ma01ag, Aze(lvol,0,ii)%i(ll).lt.0 .or. Aze(lvol,0,ii)%i(ll).gt.NN, invalid subscript )
!      if( NOTstellsym ) then
!      FATAL( ma01ag, Ato(lvol,0,ii)%i(ll).lt.0 .or. Ato(lvol,0,ii)%i(ll).gt.NN, invalid subscript )
!      FATAL( ma01ag, Azo(lvol,0,ii)%i(ll).lt.0 .or. Azo(lvol,0,ii)%i(ll).gt.NN, invalid subscript )
!      endif
!#endif

      id=Ate(lvol,0,ii)%i(ll) ; jd=1 ; dMB(id,jd) = dMB(id,jd) + half * ( Wtete(jj, pp,ii,ll) - Wtete(jj, pp,ii,Lul) ) !
      id=Aze(lvol,0,ii)%i(ll) ; jd=1 ; dMB(id,jd) = dMB(id,jd) + half * ( Wteze(jj, pp,ii,ll) - Wteze(jj, pp,ii,Lzl) ) ! Lzi;
      id=Ate(lvol,0,ii)%i(ll) ; jd=2 ; dMB(id,jd) = dMB(id,jd) + half * ( Wzete(jj,Lzp,ii,ll) - Wzete(jj,Lzp,ii,Lul) ) !
      id=Aze(lvol,0,ii)%i(ll) ; jd=2 ; dMB(id,jd) = dMB(id,jd) + half * ( Wzeze(jj,Lzp,ii,ll) - Wzeze(jj,Lzp,ii,Lzl) ) ! Lzi;
      id=Ate(lvol,0,ii)%i(ll) ; jd=1 ; dME(id,jd) = dME(id,jd) + half * ( Htete(jj, pp,ii,ll) - Htete(jj, pp,ii,Lul) ) !
      id=Aze(lvol,0,ii)%i(ll) ; jd=1 ; dME(id,jd) = dME(id,jd) + half * ( Hteze(jj, pp,ii,ll) - Hteze(jj, pp,ii,Lzl) ) ! Lzi;
      id=Ate(lvol,0,ii)%i(ll) ; jd=2 ; dME(id,jd) = dME(id,jd) + half * ( Hzete(jj,Lzp,ii,ll) - Hzete(jj,Lzp,ii,Lul) ) !
      id=Aze(lvol,0,ii)%i(ll) ; jd=2 ; dME(id,jd) = dME(id,jd) + half * ( Hzeze(jj,Lzp,ii,ll) - Hzeze(jj,Lzp,ii,Lzl) ) ! Lzi;
      if( NOTstellsym ) then
      id=Ato(lvol,0,ii)%i(ll) ; jd=1 ; dMB(id,jd) = dMB(id,jd) + half * ( Wteto(jj, pp,ii,ll) - Wteto(jj, pp,ii,Lul) ) !
      id=Azo(lvol,0,ii)%i(ll) ; jd=1 ; dMB(id,jd) = dMB(id,jd) + half * ( Wtezo(jj, pp,ii,ll) - Wtezo(jj, pp,ii,Lzl) ) ! Lzi;
      id=Ato(lvol,0,ii)%i(ll) ; jd=2 ; dMB(id,jd) = dMB(id,jd) + half * ( Wzeto(jj,Lzp,ii,ll) - Wzeto(jj,Lzp,ii,Lul) ) !
      id=Azo(lvol,0,ii)%i(ll) ; jd=2 ; dMB(id,jd) = dMB(id,jd) + half * ( Wzezo(jj,Lzp,ii,ll) - Wzezo(jj,Lzp,ii,Lzl) ) ! Lzi;
      id=Ato(lvol,0,ii)%i(ll) ; jd=1 ; dME(id,jd) = dME(id,jd) + half * ( Hteto(jj, pp,ii,ll) - Hteto(jj, pp,ii,Lul) ) !
      id=Azo(lvol,0,ii)%i(ll) ; jd=1 ; dME(id,jd) = dME(id,jd) + half * ( Htezo(jj, pp,ii,ll) - Htezo(jj, pp,ii,Lzl) ) ! Lzi;
      id=Ato(lvol,0,ii)%i(ll) ; jd=2 ; dME(id,jd) = dME(id,jd) + half * ( Hzeto(jj,Lzp,ii,ll) - Hzeto(jj,Lzp,ii,Lul) ) !
      id=Azo(lvol,0,ii)%i(ll) ; jd=2 ; dME(id,jd) = dME(id,jd) + half * ( Hzezo(jj,Lzp,ii,ll) - Hzezo(jj,Lzp,ii,Lzl) ) ! Lzi;
      endif

     enddo ! end of do pp; 24 Jan 13;
   
   !enddo ! end of do jj; 26 Feb 13;

   enddo ! end of do ll; 24 Jan 13;

  enddo ! end of do ii; 24 Jan 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!#ifdef DEBUG
!  if( Wma01ag ) write(ounit,'("ma01ag : ", 10x ," : myid=",i3," ; lvol=",i3," ; 5.0 ;")') myid, lvol
!#endif

  do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
   
  !if( Lcoordinatesingularity .and. ( mi.ne.0 .or. ni.eq.0 ) ) then ; Lzi = .true. ! 15 Jan 15; ! additional freedom; 15 Jan 15;
   if( Lcoordinatesingularity .and. ( mi.ne.0 .or. ii.eq.1 ) ) then ; Lzi = .true. ! 15 Jan 15; ! additional freedom; 15 Jan 15;
   else                                                             ; Lzi = .false.
   endif
   
   do ll = lrad-1, lrad ! this is a summation loop; 18 Feb 13;
    
    if( Lzi ) then ; Lzl = lrad
    else           ; Lzl = ll
    endif
    
    jj = 1 ; mj = im(jj) ; nj = in(jj)
   
     if( Lcoordinatesingularity .and. ( mj.ne.0 .or. nj.eq.0 ) ) then ; Lzj = .true.
     else                                                             ; Lzj = .false.
     endif
    
     do pp = lrad-1, lrad ! this is a summation loop; 18 Feb 13;
    
      if( Lzj ) then ; Lzp = lrad
      else           ; Lzp = pp
      endif

!#ifdef DEBUG
!      FATAL( ma01ag, Fso(lvol,ii)        .lt.0 .or. Fso(lvol,ii)        .gt.NN, invalid subscript )
!      if( NOTstellsym ) then
!      FATAL( ma01ag, Fse(lvol,ii)        .lt.0 .or. Fse(lvol,ii)        .gt.NN, invalid subscript )
!      endif
!#endif

      id=Fso(lvol,ii) ; jd=1 ; dMB(id,jd) = dMB(id,jd) + quart * ( mi * Wtete(jj, pp,ii,ll) - ni * Wteze(jj, pp,ii,Lzl) )
      id=Fso(lvol,ii) ; jd=2 ; dMB(id,jd) = dMB(id,jd) + quart * ( mi * Wzete(jj,Lzp,ii,ll) - ni * Wzeze(jj,Lzp,ii,Lzl) )
      id=Fso(lvol,ii) ; jd=1 ; dME(id,jd) = dME(id,jd) + quart * ( mi * Htete(jj, pp,ii,ll) - ni * Hteze(jj, pp,ii,Lzl) )
      id=Fso(lvol,ii) ; jd=2 ; dME(id,jd) = dME(id,jd) + quart * ( mi * Hzete(jj,Lzp,ii,ll) - ni * Hzeze(jj,Lzp,ii,Lzl) )
      if( NOTstellsym ) then
      id=Fse(lvol,ii) ; jd=1 ; dMB(id,jd) = dMB(id,jd) - quart * ( mi * Wteto(jj, pp,ii,ll) - ni * Wtezo(jj, pp,ii,Lzl) )
      id=Fse(lvol,ii) ; jd=2 ; dMB(id,jd) = dMB(id,jd) - quart * ( mi * Wzeto(jj,Lzp,ii,ll) - ni * Wzezo(jj,Lzp,ii,Lzl) )
      id=Fse(lvol,ii) ; jd=1 ; dME(id,jd) = dME(id,jd) - quart * ( mi * Hteto(jj, pp,ii,ll) - ni * Htezo(jj, pp,ii,Lzl) )
      id=Fse(lvol,ii) ; jd=2 ; dME(id,jd) = dME(id,jd) - quart * ( mi * Hzeto(jj,Lzp,ii,ll) - ni * Hzezo(jj,Lzp,ii,Lzl) )
      endif

     enddo ! end of do pp; 24 Jan 13;
   
!   end of jj;

   enddo ! end of do ll; 24 Jan 13;

  enddo ! end of do ii; 24 Jan 13;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!#ifdef DEBUG
!  if( Wma01ag ) write(ounit,'("ma01ag : ", 10x ," : myid=",i3," ; lvol=",i3," ; 6.0 ;")') myid, lvol
!#endif

  ii = 1 ; mi = im(ii) ; ni = in(ii)
  
  !if( Lcoordinatesingularity .and. ( mi.ne.0 .or. ni.eq.0 ) ) then ; Lzi = .true. ! 15 Jan 15; ! additional freedom; 15 Jan 15;
   if( Lcoordinatesingularity .and. ( mi.ne.0 .or. ii.eq.1 ) ) then ; Lzi = .true. ! 15 Jan 15; ! additional freedom; 15 Jan 15;
   else                                                             ; Lzi = .false.
   endif
   
   do ll = lrad-1, lrad ! this is a summation loop; 18 Feb 13;
    
    if( Lzi ) then ; Lzl = lrad
    else           ; Lzl = ll
    endif

    jj = 1 ; mj = im(jj) ; nj = in(jj)
   
     if( Lcoordinatesingularity .and. ( mj.ne.0 .or. nj.eq.0 ) ) then ; Lzj = .true.
     else                                                             ; Lzj = .false.
     endif
    
     do pp = lrad-1, lrad ! this is a summation loop; 18 Feb 13;
    
      if( Lzj ) then ; Lzp = lrad
      else           ; Lzp = pp
      endif
      
      id=1 ; jd=1 ; dMC(id,jd) = dMC(id,jd) + quart * Wtete(jj, pp,ii, ll)
      id=1 ; jd=2 ; dMC(id,jd) = dMC(id,jd) + quart * Wzete(jj,Lzp,ii, ll)
      id=2 ; jd=1 ; dMC(id,jd) = dMC(id,jd) + quart * Wteze(jj, pp,ii,Lzl)
      id=2 ; jd=2 ; dMC(id,jd) = dMC(id,jd) + quart * Wzeze(jj,Lzp,ii,Lzl) 
      id=1 ; jd=1 ; dMF(id,jd) = dMF(id,jd) + quart * Htete(jj, pp,ii, ll)
      id=1 ; jd=2 ; dMF(id,jd) = dMF(id,jd) + quart * Hzete(jj,Lzp,ii, ll)
      id=2 ; jd=1 ; dMF(id,jd) = dMF(id,jd) + quart * Hteze(jj, pp,ii,Lzl)
      id=2 ; jd=2 ; dMF(id,jd) = dMF(id,jd) + quart * Hzeze(jj,Lzp,ii,Lzl)
      
     enddo ! end of do pp; 24 Jan 13;
    
!   end of jj;

   enddo ! end of do ll; 24 Jan 13;

! end of ii

  RETURN(ma01ag)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine ma01ag

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
