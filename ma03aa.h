!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title ! Constructs matrices that define the constrained energy functional in each volume.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

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

!latex \end{enumerate} \subsection{enclosed currents} \begin{enumerate}

!latex \item In the vacuum region, the enclosed currents must be constrained.
!latex \item The plasma current is
!latex       \be I \equiv \int_{\cal S} {\bf j}\cdot d{\bf s} = \int_{\partial \cal S} {\bf B}\cdot d{\bf l} = \int_{0}^{2\pi} {\bf B}\cdot {\bf e}_\t \, d\t
!latex       \ee
!latex \item Choosing to take the line integral to lie on the inner surface (the plasma boundary), where $B^s=0$, this is
!latex       \be I = \int_{0}^{2\pi} \left( - \partial_s A_\z \; \bar g_{\t\t} + \partial_s A_\t \; \bar g_{\t\z} \right) \, d\t,
!latex       \ee
!latex       where $\bar g_{\mu\nu} \equiv g_{\mu\nu} / \sqrt g$.

!latex \item The ``linking'' current through the torus is
!latex       \be G \equiv \int_{\cal S} {\bf j}\cdot d{\bf s} = \int_{\partial \cal S} {\bf B}\cdot d{\bf l} = \int_{0}^{2\pi} {\bf B}\cdot {\bf e}_\z \, d\z
!latex       \ee
!latex \item Choosing to take the line integral to lie on the inner surface (the plasma boundary), where $B^s=0$, this is
!latex       \be G = \int_{0}^{2\pi} \left( -\partial_s A_\z \; \bar g_{\t\z} + \partial_s A_\t \; \bar g_{\z\z} \right) \, d\z.
!latex       \ee

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
!latex                 &  + &   & \ds        & \beta     & \ds \left[  \sum_l \Aze{1,l} T_l(+1) - \Delta \psi_p \right]                          \\
!latex                 &  + &   & \ds        & \gamma    & \ds \left[ \int_{0}^{2\pi} {\bf B} \cdot {\bf e}_\t \, d\t - I \right]                \\
!latex                 &  + &   & \ds        & \delta    & \ds \left[ \int_{0}^{2\pi} {\bf B} \cdot {\bf e}_\z \, d\z - G \right],
!latex     \end{array}
!latex \ee   
!latex       where
!latex       \bi
!latex       \item[  i.] $a_i$, $b_i$, $c_i$ and $d_i$ are Lagrange multipliers used to enforce the combined gauge and interface boundary condition
!latex                 on the inner interface,
!latex       \item[ ii.] $e_i$ and $f_i$               are Lagrange multipliers used to enforce the interface boundary condition on the outer interface, 
!latex       namely $\sqrt g {\bf B}\cdot\nabla s = b$;
!latex       \item[iii.] $\alpha$ and $\beta$          are Lagrange multipliers used to enforce the constraints on the enclosed fluxes; and
!latex       \item[ iv.] $\gamma$ and $\delta$         are Lagrange multipliers used to enforce the condition on the enclosed currents.
!latex       \ei

!latex \item In each plasma volume, the constraints on the enclosed currents are not required,
!latex       and the boundary condition on the outer interface is $b=0$.
!latex \item In the vacuum volume (only for free-boundary), we may set $\mu=0$ and the constraints on the enclosed fluxes are not required.

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

!latex \end{enumerate} \subsection{constructing Beltrami fields} \begin{enumerate}

!latex \item The energy, $W \equiv \int \! dv {\; \bf B}\cdot{\bf B}$, and helicity, $K\equiv \int \! dv \; {\bf A}\cdot{\bf B}$, functionals may be written
!latex       \be W & = & \frac{1}{2} \; a_i \; A_{i,j} \; a_j + a_i \; B_{i,j} \; \psi_j + \frac{1}{2} \; \psi_i \; C_{i,j} \; \psi_j \label{eq:energymatrix} \\
!latex           K & = & \frac{1}{2} \; a_i \; D_{i,j} \; a_j + a_i \; E_{i,j} \; \psi_j + \frac{1}{2} \; \psi_i \; F_{i,j} \; \psi_j \label{eq:helicitymatrix}
!latex       \ee
!latex       where ${\bf a} \equiv \{ \Ate{i,l}, \Aze{i,l}, \Ato{i,l}, \Azo{i,l}, f_{e,i}, f_{o,i} \}$ contains the independent degrees of freedom
!latex       and $\boldpsi \equiv \{\Delta \psi_t,\Delta \psi_p\}$.

!latex \item The matrix elements are computed via
!latex       \be \frac{\partial^2 W}{\partial    a_i \partial    a_j} = A_{i,j}, \qquad
!latex           \frac{\partial^2 W}{\partial    a_i \partial \psi_j} = B_{i,j}, \qquad
!latex           \frac{\partial^2 W}{\partial \psi_i \partial \psi_j} = C_{i,j}.
!latex       \ee

!latex \item The energy functionals can also be represented as
!latex       \be \begin{array}{ccccccccccccccccccccccccccccc}
!latex           W & = & \frac{1}{2} & {\bf a}^T & \cdot & A[{\bf x}] & \cdot & {\bf a} & + & B[{\bf x}, \boldpsi] & \cdot & {\bf a} & + & C[{\bf x}, \boldpsi], \\
!latex           K & = & \frac{1}{2} & {\bf a}^T & \cdot & D          & \cdot & {\bf a} & + & E[         \boldpsi] & \cdot & {\bf a} & + & F[         \boldpsi].
!latex           \end{array}
!latex       \ee

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine ma03aa( lvol, mn, lrad )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half, quart, one

  use numerical, only : sqrtmachprec

  use fileunits, only : ounit

  use inputlist, only : Wmacros, Wma03aa, Nvol

  use cputiming, only : Tma03aa

  use allglobal, only : myid, cpus, &
                        NOTstellsym, &
                        im, in, &
                        NAdof, Mvol, &
                        dMA, dMB, dMC, dMD, dME, dMF, &
                        Ate, Ato, Aze, Azo, &
                        TTee, TTeo, TToe, TToo, &
                        Lcoordinatesingularity
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER, intent(in)  :: lvol, mn, lrad

  INTEGER              :: NN

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  BEGIN(ma03aa)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef DEBUG
  FATAL(ma03aa, lvol.lt.1 .or. lvol.gt.Mvol, invalid volume label)
#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  NN = NAdof(lvol) ! shorthand;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(ma03aa)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine ma03aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
