!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Description of numerical discretization and outline of algorithm.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item For a description of the code input, see \verb+global.pdf+.

!latex \be F \equiv \sum_{i=1}^{N} \Bigg [ \; \underbrace{\int_{{\cal V}_i} \!\!\!\!dv \left(\frac{p}{\gamma-1}+\frac{B^2}{2}\right)}_{energy} - \frac{\mu_i}{2} \Bigg ( \underbrace{\int_{{\cal V}_i} \!\!\!\!dv \; {\bf A}\cdot{\bf B}}_{helicity} - K_i \Bigg ) \; \Bigg ]
!latex \ee
!latex \be \nabla \times {\bf B} = \mu {\bf B}
!latex \ee
!latex \be [[p+B^2/2]]=0
!latex \ee
!latex If $N = 1$, globally relaxed Taylor state.
!latex
!latex As $N\rightarrow \infty$, recover $\nabla p = {\bf j} \times {\bf B}$.

!latex \end{enumerate} \subsection{coordinates} \begin{enumerate}

!latex \item We shall work in coordinates, $(\s,\t,\z)$, which are be defined inversely via a transformation to Cartesian coordinates, $(x,y,z)$.
!latex \item The toroidal angle, $\z$, is identical to the cylindrical angle, $\z\equiv\phi$.
!latex \item The radial coordinate, $s$, is {\em not} a global variable: it only needs to be defined in each volume, and in each volume $\s \in [-1,1]$.
!latex \item The poloidal angle, $\t$, remains arbitrary.
!latex \item The geometry of the interfaces, ${\bf x}_v(\t,\z)$, is given by $R(\t,\z)$ and $Z(\t,\z)$ as follows:
!latex       \begin{itemize}
!latex       \item \verb+Igeometry=1+ : Cartesian
!latex       \be {\bf x} & \equiv &  \t \; {\bf \hat i} + \z \; {\bf \hat j}+ R \; {\bf \hat k}
!latex       \ee
!latex       \item \verb+Igeometry=2+ : Cylindrical
!latex       \be {\bf x} & = & R \; \cos\t \; {\bf \hat i} + R \sin\t \; {\bf \hat j} + \z \; {\bf \hat k}
!latex       \ee
!latex       \item \verb+Igeometry=3+ : Toroidal
!latex       \be {\bf x} & \equiv & R \; {\bf \hat r} + Z \; {\bf \hat k}
!latex       \ee
!latex       where ${\bf \hat r}\equiv \cos \phi \; {\bf \hat i} + \sin \phi \; {\bf \hat j}$ and
!latex             ${\bf \hat \phi} \equiv - \sin \phi \; {\bf \hat i} + \cos \phi \; {\bf \hat j}$.
!latex       \item \verb+Igeometry=4+ : Cylindrical
!latex       \be {\bf x} & \equiv & R \; {\bf \hat i} + Z \; {\bf \hat j} + \z \; {\bf \hat k}
!latex       \ee
!latex       \end{itemize}
!latex \item The "coordinate" functions, $R(\t,\z)$ and $Z(\t,\z)$, are given by a Fourier representation, as will be described below.
!latex \item The $v$-th volume is bounded by the ${\bf x}_{v-1}$ and ${\bf x}_{v}$.
!latex \item The coordinate axis $\equiv$ magnetic axis $\equiv$ degenerate interface is ${\bf x}_{0}$.
!latex \item In each {\em annular} volume, the coordinates are constructed by linear interpolation:
!latex       \be R(\s,\t,\z) &\equiv& \left[ \; (1-s) \; R_{v-1} + (1+s) \; R_{v}\; \right] \; / \; 2,\\
!latex           Z(\s,\t,\z) &\equiv& \left[ \; (1-s) \; Z_{v-1} + (1+s) \; Z_{v}\; \right] \; / \; 2.
!latex       \ee
!latex \item In the innermost, {\em simple-torus} volume, the coordinates are constructed by an interpolation that
!latex       hopefully guarantees that the interpolated coordinate surfaces do not intersect, e.g. 
!latex       \be R_{m,n}(s)-R_{m,n}(-1) = X_{m,n}(s) \; \rho^{m/2},
!latex       \ee
!latex       where $\rho \equiv (s+1)/2$ and $X_{m,n}(s)$ is an arbitrary polynomial which we choose to be constant,
!latex       \be X_{m,n}(s)=R_{m,n}(+1)-R_{m,n}(-1).
!latex       \ee
!latex       The coordinate interpolation becomes
!latex       \be R_{m,n}(s) & \equiv & (1-\rho^{m/2}) \; R_{v-1,m,n} + \rho^{m/2} \; R_{v,m,n}, \\
!latex           Z_{m,n}(s) & \equiv & (1-\rho^{m/2}) \; Z_{v-1,m,n} + \rho^{m/2} \; Z_{v,m,n}.
!latex       \ee
!latex \item The coordinate Jacobian is given by
!latex       \begin{itemize}
!latex       \item \verb+Igeometry=1+ : Cartesian
!latex       \be {\bf e}_\theta \times {\bf e}_\zeta & = & -R_\t \; \hat {\bf i} -  R_\z \; \hat {\bf j} + \hat {\bf k} \\
!latex           \boldxi \cdot {\bf e}_\theta \times {\bf e}_\zeta & =&  \delta R \\
!latex           \sqrt g                                           & =&         R_s
!latex       \ee
!latex       \item \verb+Igeometry=2+ : Cylindrical
!latex       \be {\bf e}_\theta \times {\bf e}_\zeta & = & (R_\t \sin \t + R \cos\t ) \; {\bf \hat i} + (R    \sin \t - R_\t \cos\t ) \; {\bf \hat j} - R R_\z \; {\bf \hat k} \\
!latex           \boldxi\cdot {\bf e}_\theta \times {\bf e}_\zeta & = & \delta R \; R \\
!latex           \sqrt g                                          & = & R_s \; R
!latex       \ee
!latex       \item \verb+Igeometry=3+ : Toroidal
!latex       \be {\bf e}_\theta \times {\bf e}_\zeta & = & - R \, Z_\theta \, \hat r + (Z_\theta \,R_\zeta - R_\theta \,Z_\zeta) \hat \phi + R \,R_\theta \,\hat z\\
!latex           \boldxi\cdot {\bf e}_\theta \times {\bf e}_\zeta & = & R ( \delta Z \; R_\t - \delta R \; Z_\t ) \\
!latex           \sqrt g                                         & = & R ( Z_s      \; R_\t - R_s      \; Z_\t )
!latex       \ee
!latex       \item \verb+Igeometry=4+ : Cylindrical
!latex       \be {\bf e}_\theta \times {\bf e}_\zeta & = & Z_\theta \, \hat {\bf i} - R_\t {\bf j} + ( R_\t Z_\z - Z_\t R_\z ) \,\hat {\bf k}\\
!latex           \boldxi \cdot {\bf e}_\theta \times {\bf e}_\zeta & = & \delta R \; Z_\t - \delta Z \; R_\t \\
!latex           \sqrt g                                         & = & R_s      \; Z_\t - Z_s      \; R_\t
!latex       \ee
!latex       \end{itemize}

!latex \newpage

!latex \end{enumerate} \subsection{gauge and boundary conditions} \begin{enumerate}

!latex \item In the \mbox{$v$}-th annulus, bounded by the \mbox{$(v-1)$}-th and \mbox{$v$}-th interfaces,
!latex       a general covariant representation of the magnetic vector-potential is written
!latex       \be {\bf \bar A}=\bar A_{s} \nabla s + \bar A_{\t} \nabla \t + \bar A_{\z} \nabla \z. \ee

!latex \item To this add \mbox{$\nabla g(s,\t,\z)$}, where $g$ satisfies
!latex       \be \begin{array}{cccccccccccccccccccccccccccccccccccccc}
!latex           \partial_\s g(\s,\t,\z) & = & - & \bar A_{\s}(\s,\t,\z)&   &              \\
!latex           \partial_\t g(-1,\t,\z) & = & - & \bar A_{\t}(-1,\t,\z)& + & \psi_{t,v-1} \\
!latex           \partial_\z g(-1, 0,\z) & = & - & \bar A_{\z}(-1, 0,\z)& + & \psi_{p,v-1}
!latex       \end{array}\ee
!latex       for arbitrary constants \mbox{$\psi_{t,v-1}$}, \mbox{$\psi_{p,v-1}$}. 

!latex \item Then \mbox{${\bf A}={\bf \bar A}+\nabla g$} is given by \mbox{${\bf A}=A_{\t}\nabla\t+A_{\z}\nabla\z$} with
!latex       \be \begin{array}{ccccccc}
!latex        A_{\t}(-1,\t,\z) &=& \ds \frac{\psi_{t,v-1}}{2\pi},\\
!latex        A_{\z}(-1, 0,\z) &=& \ds \frac{\psi_{p,v-1}}{2\pi}. \end{array} \label{eq:gauge}
!latex       \ee

!latex \item This specifies the gauge: to see this, notice that no gauge term can be added without violating the conditions in \Eqn{gauge}.

!latex \item Note that the gauge employed in each volume is distinct.

!latex \item The magnetic field is $\sqrt g {\bf B} = (\partial_\t A_\z - \partial_\z A_\t)\;{\bf e}_\s - \partial_\s A_\z \;{\bf e}_\t + \partial_\s A_\t \;{\bf e}_\z$.

!latex \item In the annular volumes, the condition that the field is tangential to the inner interface gives $\partial_\t A_\z - \partial_\z A_\t = 0$.
!latex       With the above gauge condition on $A_\t$ given in \Eqn{gauge}, this entails that $\partial_\t A_\z=0$, which with \Eqn{gauge} entails that
!latex       $A_\z(-1,\t,\z)=\psi_{p,v-1}/2\pi$.

!latex \item In the annular volumes, the condition that the field is tangential to the outer interface cannot be so simply enforced
!latex       (all of the gauge freedom was exploited to simplify $A_\t$ and $A_\z$ at the inner interface).
!latex       At the outer interface we must constrain the vector potential to be of the form
!latex       \be A_{\t}(1,\t,\z) &=& \partial_\t f(\t,\z), \\
!latex           A_{\z}(1,\t,\z) &=& \partial_\z f(\t,\z), \label{eq:outersurfaceconstraint} \ee
!latex       for arbitrary \mbox{$f$} of the form
!latex       \be \label{eq:surfacepotential} f = \frac{\psi_{t,v}}{2\pi} \; \t + \frac{\psi_{p,v}}{2\pi} \; \z + \tilde f(\t,\z),
!latex       \ee
!latex       where $\tilde f$ is periodic, i.e. $\tilde f(\t+2\pi,\z)=f(\t,\z)$ and $\tilde f(\t,\z+2\pi)=f(\t,\z)$.

!latex \item The toroidal and poloidal fluxes are determined using
!latex       \be \int_S {\bf B}\cdot{\bf ds}=\int_{\partial S}{\bf A}\cdot {\bf dl},
!latex       \ee
!latex       which shows that the toroidal flux is $\psi_{t,v}$ and the poloidal flux is $\psi_{p,v}$.

!latex \newpage

!latex \end{enumerate} \subsection{coordinate singularity: polar coordinates} \begin{enumerate}

!latex \item Consider the following magnetic vector potential, which is given in Cartesian coordinates,
!latex       \be {\bf A} = A_x \nabla x + A_y \nabla y +  A_z \nabla z + \nabla g \label{eq:CartesianVectorPotential}
!latex       \ee
!latex       where $A_x$, $A_y$, $A_z$, and the as-yet-arbitrary gauge function, $g$, are regular at $(x,y)=(0,0)$, i.e. they can be expanded as a Taylor series, 
!latex       \be A_x = \sum_{i,j} \alpha_{i,j} x^i y^j, \qquad
!latex           A_y = \sum_{i,j}  \beta_{i,j} x^i y^j, \qquad
!latex           A_z = \sum_{i,j} \gamma_{i,j} x^i y^j, \qquad
!latex            g  = \sum_{i,j} \delta_{i,j} x^i y^j, \label{eq:Taylorexpansion}
!latex       \ee
!latex       for small $x$ and small $y$.
!latex \item Note that we have restricted attention to the ``axisymmetric'' case, as there is no dependence on $z$.
!latex \item Consider the coordinate transformation $x = r \cos \t$, $y = r \sin \t$, $z = \z$,
!latex       which induces the vector transformation
!latex       \be \nabla x = \cos \t \; \nabla r - r \sin \t \; \nabla \t, \qquad
!latex           \nabla y = \sin \t \; \nabla r + r \cos \t \; \nabla \t, \qquad
!latex           \nabla z = \nabla \z.
!latex       \ee
!latex \item By repeated applications of the double-angle formula, the expressions for $A_x$, $A_y$ and $g$ can be cast as functions of the polar coordinates,
!latex       \be A_x & = & \sum_m r^m [ a_{m,0} + a_{m,1} \; r^2 + a_{m,2} \; r^4 + \dots ] \sin(m\t), \\
!latex           A_y & = & \sum_m r^m [ b_{m,0} + b_{m,1} \; r^2 + b_{m,2} \; r^4 + \dots ] \cos(m\t), \\
!latex           A_z & = & \sum_m r^m [ c_{m,0} + c_{m,1} \; r^2 + c_{m,2} \; r^4 + \dots ] \cos(m\t), \label{eq:regularAz} \\
!latex            g  & = & \sum_m r^m [ g_{m,0} + g_{m,1} \; r^2 + g_{m,2} \; r^4 + \dots ] \sin(m\t), 
!latex       \ee
!latex       where we have restricted attention to stellarator symmetric geometry.
!latex \item Collecting these expressions, the vector potential can be expressed in polar coordinates,
!latex       \be {\bf A} = A_r \nabla r + A_\t \nabla \t + A_\z \nabla \z + \partial_r g \; \nabla r + \partial_\t g \; \nabla \t,
!latex       \ee
!latex       where
!latex       \be \begin{array}{ccccccclcccccccccccccccccccccccccccccccccccccccccccccccccccc}
!latex       A_r & = & r^0 & [ & ( &           &   & b_{0,0}   & + &   g_{1,0} & ) & + & (\dots) r^2 + (\dots)r^4 + \dots & ] & \sin   \t \\
!latex           & + & r^1 & [ & ( & a_{1,0}/2 & + & b_{1,0}/2 & + & 2 g_{2,0} & ) & + & (\dots) r^2 + (\dots)r^4 + \dots & ] & \sin 2 \t \\
!latex           & + & r^2 & [ & ( & a_{2,0}/2 & + & b_{2,0}/2 & + & 3 g_{3,0} & ) & + & (\dots) r^2 + (\dots)r^4 + \dots & ] & \sin 3 \t \\
!latex           & + & \dots
!latex       \end{array} \ee
!latex       (Note: Mathematica was used to perform the algebraic manipulations.)
!latex \item There is precisely enough gauge freedom so that we may choose $A_r=0$.
!latex       For example, the choice
!latex       \be \begin{array}{cccccclccccccccccccccccccccccccccccccccc}
!latex           g_{1,0} & = & - &   &           &   & b_{0,0},   &   &   &    \\
!latex           g_{2,0} & = & - & ( & a_{1,0}/2 & + & b_{1,0}/2  & ) & / & 2, \\
!latex           g_{3,0} & = & - & ( & a_{2,0}/2 & + & b_{2,0}/2  & ) & / & 3, \\
!latex           \dots   & = &   & \dots
!latex       \end{array} \ee
!latex       eliminates the lowest order $r$ dependence in each harmonic.
!latex \item By working through the algebra (again, using Mathematica) the expression for $A_\t$ becomes
!latex       \be A_\t = r^2 f_0(s) + r^3 f_1(s) \cos(\t) + r^4 f_2(s) \cos(2\t) + r^5 f_3(s) \cos(2\t) + \dots
!latex       \ee
!latex       where $s\equiv r^2$ and the $f_m(s)$ are abitrary polynomials in $s$.
!latex \item The expression for $A_\z$ is unchanged from \Eqn{regularAz},
!latex       \be A_\z = r^0 g_0(s) + r^1 g_1(s) \cos(\t) + r^2 g_2(s) \cos(2\t) + r^3 g_3(s) \cos(2\t) + \dots \label{eq:nearoriginAz}
!latex       \ee

!latex \end{enumerate} \subsection{coordinate axis coinciding with the magnetic axis constraint} \begin{enumerate}

!latex \item Consider the case where the magnetic axis is forced to coincide with the coordinate axis.
!latex       Returning to Cartesian coordinates and the expression for the vector potential given in \Eqn{CartesianVectorPotential}, the magnetic field is
!latex       \be {\bf B} = \partial_y A_z \; {\bf e}_x - \partial_x A_z \; {\bf e}_y + ( \partial_x A_y - \partial_y A_x ) \; {\bf e}_z.
!latex       \ee
!latex \item The constraint that the magnetic axis coincides with the coordinate axis is
!latex       \be \partial_x A_z(0,0) = 0, \qquad \partial_y A_z(0,0) = 0,
!latex       \ee
!latex       which, using \Eqn{Taylorexpansion}, is enforced by setting $\gamma_{1,0}=0$ and $\gamma_{0,1}=0$.
!latex \item The expression for $A_\z$ becomes
!latex       \be A_\z = r^0 g_0(s) + r^3 g_1(s) \cos(\t) + r^2 g_2(s) \cos(2\t) + r^3 g_3(s) \cos(2\t) + \dots \label{eq:nearoriginAzaxisconstraint}
!latex       \ee
!latex \item Note that if this constraint is enforced, in addition to the independent degrees-of-freedom that describe the magnetic vector potential, 
!latex       the location of the coordinate axis should be considered as a degree-of-freedom in the description of the magnetic field.

!latex \newpage 

!latex \end{enumerate} \subsection{mixed Fourier, Chebyshev representation} \begin{enumerate}

!latex \item Doubly periodic, even functions are represented as
!latex       \be f(\t,\z) & = & \sum_{n=0}^{N} f_{0,n}\cos( -n \z ) + \sum_{m=1}^{M}\sum_{n=-N}^{N} f_{m,n}\cos( m \t - n \z ) \equiv \sum_i f_i \cos\a_i,
!latex       \ee
!latex       where $\a_i \equiv m_i \t - n_i \z$, and similarly for odd functions.

!latex \item The coordinate functions are
!latex       \be            R_v(\t,\z) & = & \sum_i R_{v,e,i} \; \cos\a_i + \sum_i R_{v,o,i} \; \sin\a_i \\
!latex                      Z_v(\t,\z) & = & \sum_i Z_{v,e,i} \; \cos\a_i + \sum_i Z_{v,o,i} \; \sin\a_i
!latex       \ee

!latex \item The covariant components of the vector potential are written as the following sum, 
!latex       \be            A_\t & = & \sum_i \sum_{l=0}^L A_{\t,e,i,l} \; T_{l}(s) \cos\a_i + \sum_i \sum_{l=0}^L A_{\t,o,i,l} \; T_{l}(s) \sin\a_i \\
!latex                      A_\z & = & \sum_i \sum_{l=0}^L A_{\z,e,i,l} \; T_{l}(s) \cos\a_i + \sum_i \sum_{l=0}^L A_{\z,o,i,l} \; T_{l}(s) \sin\a_i ,
!latex       \ee
!latex       where $T_l(s)$ are the Chebyshev polynomials.

!latex \item The boundary condition on the outer iterface is constrained by
!latex       \be f = \sum_i f_{e,i}\cos\alpha_i+ \sum_i f_{o,i}\sin\alpha_i
!latex       \ee

!latex \item The Chebyshev polynomials may be defined recursively:
!latex       \be T_{  0}(s)&=&1,\\
!latex           T_{  1}(s)&=&s,\\
!latex           T_{l+1}(s)&=&2 s \; T_{l}(s) - T_{l-1}(s).
!latex       \ee

!latex \item Note that $T_l(-1) = +1$ if $l$ is even, and $T_l(-1) = -1$ if $l$ is odd, and that $T_l(+1) = +1$.

!latex \item Also, note that $T_l^\prime(-1)=l^2(-1)^{(l+1)}$ and $T_l^\prime(+1)=l^2$. ! 14 Jan 15;

!latex \item The magnetic field is, with summation over $i$ and $l$ hereafter implied, 
!latex       \be \begin{array}{ccccccccccccccccccccccccccccc}
!latex       \sqrt g B^\s=&+& (m_i A_{\z,o,i,l}+n_i A_{\t,o,i,l}) & T_{l}        \; \cos\a_i &-& (m_i A_{\z,e,i,l}+n_i A_{\t,e,i,l}) & T_{l}        \; \sin\a_i \\
!latex       \sqrt g B^\t=&-&      A_{\z,e,i,l}                   & T_{l}^\prime \; \cos\a_i &-&      A_{\z,o,i,l}                   & T_{l}^\prime \; \sin\a_i \\
!latex       \sqrt g B^\z=&+&      A_{\t,e,i,l}                   & T_{l}^\prime \; \cos\a_i &+&      A_{\t,o,i,l}                   & T_{l}^\prime \; \sin\a_i  
!latex       \end{array}
!latex       \ee

!latex \newpage

!latex \end{enumerate} \subsection{interface boundary conditions: annular volumes} \begin{enumerate}

!latex \item The boundary conditions are, where $\Delta\psi_{t,v} \equiv \psi_{t,v} - \psi_{t,v-1}$ and $\Delta\psi_{p,v} \equiv \psi_{p,v} - \psi_{p,v-1}$, 

!latex \be A_{\t,e,1}(-1) = 0 \; & \; \mbox{\rm and} \; & \;\; A_{\t,e,1}(+1) = \Delta\psi_{t,v} \\
!latex     A_{\t,e,i}(-1) = 0 \; & \; \mbox{\rm and} \; & \;\; A_{\t,e,i}(+1) = +m_i f_{o,i}     \\
!latex     A_{\t,o,i}(-1) = 0 \; & \; \mbox{\rm and} \; & \;\; A_{\t,o,i}(+1) = -m_i f_{e,i}     \\
!latex     A_{\z,e,1}(-1) = 0 \; & \; \mbox{\rm and} \; & \;\; A_{\z,e,1}(+1) = \Delta\psi_{p,v} \\
!latex     A_{\z,e,i}(-1) = 0 \; & \; \mbox{\rm and} \; & \;\; A_{\z,e,i}(+1) = -n_i f_{o,i}     \\
!latex     A_{\z,o,i}(-1) = 0 \; & \; \mbox{\rm and} \; & \;\; A_{\z,o,i}(+1) = +n_i f_{e,i}       
!latex \ee

!latex \item The boundary conditions may be enforced using the Chebyshev representation, $f(s) = \sum f_l T_l(s)$, as follows:
!latex       \be \begin{array}{cccccccccccccccccccccccccccccccccccccccc}
!latex           f(-1) & = & \ds \sum_{l=0}^{L-2} & f_l & (-1)^l & + & f_{L-1} & (-1)^{L-1} & + & f_{L} & (-1)^{L} \\ \\
!latex           f(+1) & = & \ds \sum_{l=0}^{L-2} & f_l &        & + & f_{L-1} &            & + & f_{L} &         
!latex       \end{array} \ee

!latex \item We consider the $f_{L-1}$ and $f_{L}$ to depend on the $f_{0}$, $f_1$, $\dots$, $f_{L-2}$, as follows
!latex       \be f_{L-1} & = & \left\{ \begin{array}{ccc} \left[ f(+1) - f(-1) \right] / 2 - f_1 - f_3 \dots - f_{L-3},& L \;\; {\rm even} \\
!latex                                                    \left[ f(+1) + f(-1) \right] / 2 - f_0 - f_2 \dots - f_{L-3},& L \;\; {\rm odd} 
!latex                                 \end{array} \right. \\
!latex           f_{L  } & = & \left\{ \begin{array}{ccc} \left[ f(+1) + f(-1) \right] / 2 - f_0 - f_2 \dots - f_{L-2},& L \;\; {\rm even} \\
!latex                                                    \left[ f(+1) - f(-1) \right] / 2 - f_1 - f_3 \dots - f_{L-2},& L \;\; {\rm odd} 
!latex                                 \end{array} \right.
!latex       \ee

!latex \item Hereafter, we will assume that $L$ is even, so we have
!latex       \be \frac{\partial}{\partial f_l} & \equiv & \frac{\partial}{\partial f_l} - \frac{\partial}{\partial f_{L  }}, \;\;\;\; l \;\; {\rm even} \\
!latex           \frac{\partial}{\partial f_l} & \equiv & \frac{\partial}{\partial f_l} - \frac{\partial}{\partial f_{L-1}}, \;\;\;\; l \;\; {\rm odd }
!latex       \ee

!l!tex \item For $L$  odd, we have
!l!tex       \be \frac{\partial}{\partial f_l} & \equiv & \frac{\partial}{\partial f_l} - \frac{\partial}{\partial f_{L-1}}, \;\;\;\; l \;\; {\rm even} \\
!l!tex           \frac{\partial}{\partial f_l} & \equiv & \frac{\partial}{\partial f_l} - \frac{\partial}{\partial f_{L  }}, \;\;\;\; l \;\; {\rm odd }
!l!tex       \ee

!latex \item Using loose notation, the first derivatives with respect to the vector potential degress of freedom are
!latex       \be \frac{\partial}{\partial A_{\alpha,u,i,l}} = \frac{\partial}{\partial A_{\alpha,u,i,l  }} - \frac{\partial}{\partial A_{\alpha,u,i,L_l}},
!latex       \ee
!latex       where $L_l=L$ if $l$ is even and $L_l=L-1$ if $l$ is odd.

!latex \item The derivatives with respect to the surface potential are
!latex       \be \frac{\partial}{\partial f_{o,i}} & = & + \frac{1}{2} m_i \sum_{l=L-1}^{L} \frac{\partial}{\partial \Ate{i,l}} 
!latex                                                   - \frac{1}{2} n_i \sum_{l=L-1}^{L} \frac{\partial}{\partial \Aze{i,l}} \\
!latex           \frac{\partial}{\partial f_{e,i}} & = & - \frac{1}{2} m_i \sum_{l=L-1}^{L} \frac{\partial}{\partial \Ato{i,l}} 
!latex                                                   + \frac{1}{2} n_i \sum_{l=L-1}^{L} \frac{\partial}{\partial \Azo{i,l}}
!latex       \ee

!latex \item The derivatives with respect to the enclosed fluxes are
!latex       \be \frac{\partial}{\partial \Delta \psi_t} & = & \frac{1}{2}\sum_{l=L-1}^{L} \frac{\partial}{\partial \Ate{1,l}} \\
!latex           \frac{\partial}{\partial \Delta \psi_p} & = & \frac{1}{2}\sum_{l=L-1}^{L} \frac{\partial}{\partial \Aze{1,l}}
!latex       \ee

!latex \item The required second derivatives are
!latex       \be  \!\!\!\!\!\!\!\! \frac{\partial}{\partial A_{\beta,v,j,p}}\frac{\partial}{\partial A_{\alpha,u,i,l}} & = &
!latex            + \frac{\partial}{\partial A_{\beta,v,j,  p}}\frac{\partial}{\partial A_{\alpha,u,i,l  }}
!latex            - \frac{\partial}{\partial A_{\beta,v,j,  p}}\frac{\partial}{\partial A_{\alpha,u,i,L_l}}
!latex            - \frac{\partial}{\partial A_{\beta,v,j,L_p}}\frac{\partial}{\partial A_{\alpha,u,i,l  }}
!latex            + \frac{\partial}{\partial A_{\beta,v,j,L_p}}\frac{\partial}{\partial A_{\alpha,u,i,L_l}} \nonumber
!latex       \ee

!latex       \be \!\!\!\!\!\!\!\! \frac{\partial}{\partial f_{o,j}}\frac{\partial}{\partial A_{\alpha,u,i,l}} =
!latex           & + & \frac{1}{2} m_j \sum_{p=L-1}^{L}\left(\frac{\partial}{\partial \Ate{j,p}}\frac{\partial}{\partial A_{\alpha,u,i,  l}}
!latex                                                      -\frac{\partial}{\partial \Ate{j,p}}\frac{\partial}{\partial A_{\alpha,u,i,L_l}}
!latex                                                 \right) \nonumber \\
!latex           & - & \frac{1}{2} n_j \sum_{p=L-1}^{L}\left(\frac{\partial}{\partial \Aze{j,p}}\frac{\partial}{\partial A_{\alpha,u,i,  l}}
!latex                                                      -\frac{\partial}{\partial \Aze{j,p}}\frac{\partial}{\partial A_{\alpha,u,i,L_l}}
!latex                                                 \right)
!latex       \ee

!latex       \be \!\!\!\!\!\!\!\! \frac{\partial}{\partial f_{e,j}}\frac{\partial}{\partial A_{\alpha,u,i,l}} =
!latex           & - & \frac{1}{2} m_j \sum_{p=L-1}^{L}\left(\frac{\partial}{\partial \Ato{j,p}}\frac{\partial}{\partial A_{\alpha,u,i,  l}}
!latex                                                      -\frac{\partial}{\partial \Ato{j,p}}\frac{\partial}{\partial A_{\alpha,u,i,L_l}}
!latex                                                 \right) \nonumber \\
!latex           & + & \frac{1}{2} n_j \sum_{p=L-1}^{L}\left(\frac{\partial}{\partial \Azo{j,p}}\frac{\partial}{\partial A_{\alpha,u,i,  l}}
!latex                                                      -\frac{\partial}{\partial \Azo{j,p}}\frac{\partial}{\partial A_{\alpha,u,i,L_l}}
!latex                                                 \right)
!latex       \ee

!latex       \be  \!\!\!\!\!\!\!\! \frac{\partial}{\partial \Delta \psi_{t}}\frac{\partial}{\partial A_{\alpha,u,i,l}} =
!latex           & + & \frac{1}{2} \sum_{p=L-1}^{L} \left( \frac{\partial}{\partial \Ate{1,p}}\frac{\partial}{\partial A_{\alpha,u,i,l  }}
!latex                                                   - \frac{\partial}{\partial \Ate{1,p}}\frac{\partial}{\partial A_{\alpha,u,i,L_l}}
!latex                                              \right)
!latex       \ee

!latex       \be  \!\!\!\!\!\!\!\! \frac{\partial}{\partial \Delta \psi_{p}}\frac{\partial}{\partial A_{\alpha,u,i,l}} =
!latex           & + & \frac{1}{2} \sum_{p=L-1}^{L} \left( \frac{\partial}{\partial \Aze{1,p}}\frac{\partial}{\partial A_{\alpha,u,i,l  }}
!latex                                                   - \frac{\partial}{\partial \Aze{1,p}}\frac{\partial}{\partial A_{\alpha,u,i,L_l}}
!latex                                              \right)
!latex       \ee

!latex \item The required second derivatives are

!latex       \be  \!\!\!\!\!\!\!\! \frac{\partial}{\partial A_{\a,u,j,p}}\frac{\partial}{\partial f_{o,i}} =
!latex           & + & \frac{1}{2} m_i \sum_{l=L-1}^{L} \left( \frac{\partial}{\partial A_{\a,u,j,  p}}\frac{\partial}{\partial \Ate{i,l}}
!latex                                                       - \frac{\partial}{\partial A_{\a,u,j,L_p}}\frac{\partial}{\partial \Ate{i,l}}
!latex                                                  \right) \nonumber \\
!latex           & - & \frac{1}{2} n_i \sum_{l=L-1}^{L} \left( \frac{\partial}{\partial A_{\a,u,j,  p}}\frac{\partial}{\partial \Aze{i,l}}
!latex                                                       - \frac{\partial}{\partial A_{\a,u,j,L_p}}\frac{\partial}{\partial \Aze{i,l}}
!latex                                                  \right) 
!latex       \ee

!latex       \be  \!\!\!\!\!\!\!\! \frac{\partial}{\partial f_{o,j}}\frac{\partial}{\partial f_{o,i}} = 
!latex            & + & \frac{1}{4} m_j m_i \sum_{p=L-1}^{L} \sum_{l=L-1}^{L} \frac{\partial}{\partial \Ate{j,p}} \frac{\partial}{\partial \Ate{i,l}} 
!latex              -   \frac{1}{4} m_j n_i \sum_{p=L-1}^{L} \sum_{l=L-1}^{L} \frac{\partial}{\partial \Ate{j,p}} \frac{\partial}{\partial \Aze{i,l}} \nonumber \\
!latex            & - & \frac{1}{4} n_j m_i \sum_{p=L-1}^{L} \sum_{l=L-1}^{L} \frac{\partial}{\partial \Aze{j,p}} \frac{\partial}{\partial \Ate{i,l}} 
!latex              +   \frac{1}{4} n_j n_i \sum_{p=L-1}^{L} \sum_{l=L-1}^{L} \frac{\partial}{\partial \Aze{j,p}} \frac{\partial}{\partial \Aze{i,l}}
!latex       \ee

!latex       \be  \!\!\!\!\!\!\!\! \frac{\partial}{\partial f_{e,j}}\frac{\partial}{\partial f_{o,i}} = \nonumber
!latex            & - & \frac{1}{4} m_j m_i \sum_{p=L-1}^{L} \sum_{l=L-1}^{L} \frac{\partial}{\partial \Ato{j,p}} \frac{\partial}{\partial \Ate{i,l}} 
!latex              +   \frac{1}{4} m_j n_i \sum_{p=L-1}^{L} \sum_{l=L-1}^{L} \frac{\partial}{\partial \Ato{j,p}} \frac{\partial}{\partial \Aze{i,l}} \nonumber \\
!latex            & + & \frac{1}{4} n_j m_i \sum_{p=L-1}^{L} \sum_{l=L-1}^{L} \frac{\partial}{\partial \Azo{j,p}} \frac{\partial}{\partial \Ate{i,l}} 
!latex              -   \frac{1}{4} n_j n_i \sum_{p=L-1}^{L} \sum_{l=L-1}^{L} \frac{\partial}{\partial \Azo{j,p}} \frac{\partial}{\partial \Aze{i,l}}
!latex       \ee

!latex       \be  \!\!\!\!\!\!\!\! \frac{\partial}{\partial \Delta \psi_t}\frac{\partial}{\partial f_{o,i}} = \nonumber
!latex            & + & \frac{1}{4}     m_i \sum_{p=L-1}^{L} \sum_{l=L-1}^{L} \frac{\partial}{\partial \Ate{1,p}} \frac{\partial}{\partial \Ate{i,l}} 
!latex              -   \frac{1}{4}     n_i \sum_{p=L-1}^{L} \sum_{l=L-1}^{L} \frac{\partial}{\partial \Ate{1,p}} \frac{\partial}{\partial \Aze{i,l}}
!latex       \ee

!latex       \be  \!\!\!\!\!\!\!\! \frac{\partial}{\partial \Delta \psi_p}\frac{\partial}{\partial f_{o,i}} = \nonumber
!latex            & + & \frac{1}{4}     m_i \sum_{p=L-1}^{L} \sum_{l=L-1}^{L} \frac{\partial}{\partial \Aze{1,p}} \frac{\partial}{\partial \Ate{i,l}} 
!latex              -   \frac{1}{4}     n_i \sum_{p=L-1}^{L} \sum_{l=L-1}^{L} \frac{\partial}{\partial \Aze{1,p}} \frac{\partial}{\partial \Aze{i,l}}
!latex       \ee

!latex \item The required second derivatives are

!latex       \be  \!\!\!\!\!\!\!\! \frac{\partial}{\partial A_{\a,u,j,p}}\frac{\partial}{\partial f_{e,i}} =
!latex           & - & \frac{1}{2} m_i \sum_{l=L-1}^{L} \left( \frac{\partial}{\partial A_{\a,u,j,  p}}\frac{\partial}{\partial \Ato{i,l}}
!latex                                                       - \frac{\partial}{\partial A_{\a,u,j,L_p}}\frac{\partial}{\partial \Ato{i,l}}
!latex                                                  \right) \nonumber \\
!latex           & + & \frac{1}{2} n_i \sum_{l=L-1}^{L} \left( \frac{\partial}{\partial A_{\a,u,j,  p}}\frac{\partial}{\partial \Azo{i,l}}
!latex                                                       - \frac{\partial}{\partial A_{\a,u,j,L_p}}\frac{\partial}{\partial \Azo{i,l}}
!latex                                                  \right) 
!latex       \ee

!latex       \be  \!\!\!\!\!\!\!\! \frac{\partial}{\partial f_{o,j}}\frac{\partial}{\partial f_{e,i}} = 
!latex            & - & \frac{1}{4} m_j m_i \sum_{p=L-1}^{L} \sum_{l=L-1}^{L} \frac{\partial}{\partial \Ate{j,p}} \frac{\partial}{\partial \Ato{i,l}} 
!latex              +   \frac{1}{4} m_j n_i \sum_{p=L-1}^{L} \sum_{l=L-1}^{L} \frac{\partial}{\partial \Ate{j,p}} \frac{\partial}{\partial \Azo{i,l}} \nonumber \\
!latex            & + & \frac{1}{4} n_j m_i \sum_{p=L-1}^{L} \sum_{l=L-1}^{L} \frac{\partial}{\partial \Aze{j,p}} \frac{\partial}{\partial \Ato{i,l}} 
!latex              -   \frac{1}{4} n_j n_i \sum_{p=L-1}^{L} \sum_{l=L-1}^{L} \frac{\partial}{\partial \Aze{j,p}} \frac{\partial}{\partial \Azo{i,l}}
!latex       \ee

!latex       \be  \!\!\!\!\!\!\!\! \frac{\partial}{\partial f_{e,j}}\frac{\partial}{\partial f_{e,i}} = \nonumber
!latex            & + & \frac{1}{4} m_j m_i \sum_{p=L-1}^{L} \sum_{l=L-1}^{L} \frac{\partial}{\partial \Ato{j,p}} \frac{\partial}{\partial \Ato{i,l}} 
!latex              -   \frac{1}{4} m_j n_i \sum_{p=L-1}^{L} \sum_{l=L-1}^{L} \frac{\partial}{\partial \Ato{j,p}} \frac{\partial}{\partial \Azo{i,l}} \nonumber \\
!latex            & - & \frac{1}{4} n_j m_i \sum_{p=L-1}^{L} \sum_{l=L-1}^{L} \frac{\partial}{\partial \Azo{j,p}} \frac{\partial}{\partial \Ato{i,l}} 
!latex              +   \frac{1}{4} n_j n_i \sum_{p=L-1}^{L} \sum_{l=L-1}^{L} \frac{\partial}{\partial \Azo{j,p}} \frac{\partial}{\partial \Azo{i,l}}
!latex       \ee

!latex       \be  \!\!\!\!\!\!\!\! \frac{\partial}{\partial \Delta \psi_t}\frac{\partial}{\partial f_{e,i}} = \nonumber
!latex            & - & \frac{1}{4}     m_i \sum_{p=L-1}^{L} \sum_{l=L-1}^{L} \frac{\partial}{\partial \Ate{1,p}} \frac{\partial}{\partial \Ato{i,l}} 
!latex              +   \frac{1}{4}     n_i \sum_{p=L-1}^{L} \sum_{l=L-1}^{L} \frac{\partial}{\partial \Ate{1,p}} \frac{\partial}{\partial \Azo{i,l}}
!latex       \ee

!latex       \be  \!\!\!\!\!\!\!\! \frac{\partial}{\partial \Delta \psi_p}\frac{\partial}{\partial f_{e,i}} = \nonumber
!latex            & - & \frac{1}{4}     m_i \sum_{p=L-1}^{L} \sum_{l=L-1}^{L} \frac{\partial}{\partial \Aze{1,p}} \frac{\partial}{\partial \Ato{i,l}} 
!latex              +   \frac{1}{4}     n_i \sum_{p=L-1}^{L} \sum_{l=L-1}^{L} \frac{\partial}{\partial \Aze{1,p}} \frac{\partial}{\partial \Azo{i,l}}
!latex       \ee


!latex \item The required second derivatives are

!latex       \be  \!\!\!\!\!\!\!\! \frac{\partial}{\partial A_{\a,u,j,l}}\frac{\partial}{\partial \Delta \psi_t} =
!latex            & + & \frac{1}{2}                \sum_{l=L-1}^{L} \frac{\partial}{\partial A_{\alpha,u,j,  p}} \frac{\partial}{\partial \Ate{1,l}} 
!latex              -   \frac{1}{2}                \sum_{l=L-1}^{L} \frac{\partial}{\partial A_{\alpha,u,j,L_p}} \frac{\partial}{\partial \Ate{1,l}} 
!latex       \ee

!latex       \be  \!\!\!\!\!\!\!\! \frac{\partial}{\partial f_{o,i}}\frac{\partial}{\partial \Delta \psi_t} =
!latex            & + & \frac{1}{4} m_j     \sum_{p=L-1}^{L} \sum_{l=L-1}^{L} \frac{\partial}{\partial \Ate{j,p}} \frac{\partial}{\partial \Ate{1,l}} 
!latex              -   \frac{1}{4} n_j     \sum_{p=L-1}^{L} \sum_{l=L-1}^{L} \frac{\partial}{\partial \Aze{j,p}} \frac{\partial}{\partial \Ate{1,l}} 
!latex       \ee

!latex       \be  \!\!\!\!\!\!\!\! \frac{\partial}{\partial f_{e,i}}\frac{\partial}{\partial \Delta \psi_t} =
!latex            & - & \frac{1}{4} m_j     \sum_{p=L-1}^{L} \sum_{l=L-1}^{L} \frac{\partial}{\partial \Ato{j,p}} \frac{\partial}{\partial \Ate{1,l}} 
!latex              +   \frac{1}{4} n_j     \sum_{p=L-1}^{L} \sum_{l=L-1}^{L} \frac{\partial}{\partial \Azo{j,p}} \frac{\partial}{\partial \Ate{1,l}} 
!latex       \ee

!latex       \be  \!\!\!\!\!\!\!\! \frac{\partial}{\partial \Delta \psi_t}\frac{\partial}{\partial \Delta \psi_t} =
!latex            \frac{1}{4} \sum_{p=L-1}^{L} \sum_{l=L-1}^{L} \frac{\partial}{\partial \Ate{1,p}}   \frac{\partial}{\partial \Ate{1,l}}
!latex       \ee

!latex       \be  \!\!\!\!\!\!\!\! \frac{\partial}{\partial \Delta \psi_p}\frac{\partial}{\partial \Delta \psi_t} =
!latex            \frac{1}{4} \sum_{p=L-1}^{L} \sum_{l=L-1}^{L} \frac{\partial}{\partial \Aze{1,p}}   \frac{\partial}{\partial \Ate{1,l}}
!latex       \ee

!latex \item The required second derivatives are

!latex       \be  \!\!\!\!\!\!\!\! \frac{\partial}{\partial A_{\a,u,j,l}}\frac{\partial}{\partial \Delta \psi_p} =
!latex            & + & \frac{1}{2}                \sum_{l=L-1}^{L} \frac{\partial}{\partial A_{\alpha,u,j,  p}} \frac{\partial}{\partial \Aze{1,l}} 
!latex              -   \frac{1}{2}                \sum_{l=L-1}^{L} \frac{\partial}{\partial A_{\alpha,u,j,L_p}} \frac{\partial}{\partial \Aze{1,l}} 
!latex       \ee

!latex       \be  \!\!\!\!\!\!\!\! \frac{\partial}{\partial f_{o,i}}\frac{\partial}{\partial \Delta \psi_p} =
!latex            & + & \frac{1}{4} m_j     \sum_{p=L-1}^{L} \sum_{l=L-1}^{L} \frac{\partial}{\partial \Ate{j,p}} \frac{\partial}{\partial \Aze{1,l}} 
!latex              -   \frac{1}{4} n_j     \sum_{p=L-1}^{L} \sum_{l=L-1}^{L} \frac{\partial}{\partial \Aze{j,p}} \frac{\partial}{\partial \Aze{1,l}} 
!latex       \ee

!latex       \be  \!\!\!\!\!\!\!\! \frac{\partial}{\partial f_{e,i}}\frac{\partial}{\partial \Delta \psi_p} =
!latex            & - & \frac{1}{4} m_j     \sum_{p=L-1}^{L} \sum_{l=L-1}^{L} \frac{\partial}{\partial \Ato{j,p}} \frac{\partial}{\partial \Aze{1,l}} 
!latex              +   \frac{1}{4} n_j     \sum_{p=L-1}^{L} \sum_{l=L-1}^{L} \frac{\partial}{\partial \Azo{j,p}} \frac{\partial}{\partial \Aze{1,l}} 
!latex       \ee

!latex       \be  \!\!\!\!\!\!\!\! \frac{\partial}{\partial \Delta \psi_t}\frac{\partial}{\partial \Delta \psi_p} =
!latex            \frac{1}{4} \sum_{p=L-1}^{L} \sum_{l=L-1}^{L} \frac{\partial}{\partial \Ate{1,p}}   \frac{\partial}{\partial \Aze{1,l}}
!latex       \ee

!latex       \be  \!\!\!\!\!\!\!\! \frac{\partial}{\partial \Delta \psi_p}\frac{\partial}{\partial \Delta \psi_p} =
!latex            \frac{1}{4} \sum_{p=L-1}^{L} \sum_{l=L-1}^{L} \frac{\partial}{\partial \Aze{1,p}}   \frac{\partial}{\partial \Aze{1,l}}
!latex       \ee

!latex \newpage

!latex \end{enumerate} \subsection{interface boundary conditions: coordinate singularity} \begin{enumerate}

!latex \item In the innermost volume, we shall still use the interpolation coordinate $-1 \le s \le 1$,
!latex       where $s=-1$ is the coordinate singularity and $s=1$ is the innermost interface.
!latex \item For convenience, the following shall employ $\rho \equiv (s+1)/2$.
!latex \item The gauge and regularity conditions on the vector potential, which may generally be written
!latex       \be {\bf A} = \sum_{m,n}A_{\t,m,n} e^{i(m\t-n\z)} \nabla \theta + \sum_{m,n}A_{\z,m,n} e^{i(m\t-n\z)} \nabla \zeta,
!latex       \ee
!latex       require that near the coordinate origin, when the magnetic axis is not forced to coincide with the coordinate origin,
!latex       the vector potential harmonics have the following radial dependencies:
!latex       \be \begin{array}{cclcccccccccccccccccccccccccc}
!latex           A_{\t,m,n} & \sim & r^{m+2} & f_{m,n}(\rho),                    &                 &          \\
!latex           A_{\z,m,n} & \sim &         & g_{m,n}(\rho), & \qquad {\rm for} & \qquad m  =  0, & n  =  0, \\
!latex           A_{\z,m,n} & \sim & r^{m+2} & g_{m,n}(\rho), & \qquad {\rm for} & \qquad m  =  0, & n \ne 0, \\
!latex           A_{\z,m,n} & \sim & r^{m  } & g_{m,n}(\rho), & \qquad {\rm for} & \qquad m \ne 0, &
!latex       \end{array} \label{eq:originconstraints} \ee
!latex       where the $f_{m,n}(\rho)$ and $g_{m,n}(\rho)$ are arbitrary polynomials in $\rho$.
!latex \item (Note that the representation for $A_\z$ has changed slightly from that given in \Eqn{nearoriginAz} and \Eqn{nearoriginAzaxisconstraint},
!latex       because the additional gauge freedom, $\partial_\z g$, has been used.)
!latex \item To simplify the algorithmic implementation of these combined regularity-gauge conditions, we shall use
!latex       \be A_{\t,m,n} & = & \rho^{m/2} \alpha_{m,n}(\rho), \\
!latex           A_{\z,m,n} & = & \rho^{m/2} \beta_{m,n}(\rho), \label{eq:regularitygauge}
!latex       \ee
!latex       where $\alpha_{m,n}(\rho)$ and $\beta_{m,n}(\rho)$ are arbitrary polynomials in $\rho$,
!latex       except for the constraints $\alpha_{m,n}(0)=0$ for all $(m,n)$, and $\beta_{m,n}(0)=0$ for $m=0$ and $n \ne 0$.
!latex \item Note that the magnetic flux linking the torus is $2 \pi A_{\z,0,0}(1)$. This is arbitrary, so we choose $A_{\z,0,0}(1)=0$.

!latex \item When the magnetic axis is forced to coincide with the coordinate origin, \Eqn{originconstraints} becomes
!latex       \be \begin{array}{cclcccccccccccccccccccccccccc}
!latex           A_{\t,m,n} & \sim & r^{m+2} & f_{m,n}(\rho),                    &                 &          \\
!latex           A_{\z,m,n} & \sim &         & g_{m,n}(\rho), & \qquad {\rm for} & \qquad m  =  0, & n  =  0, \\
!latex           A_{\z,m,n} & \sim & r^{m+2} & g_{m,n}(\rho), & \qquad {\rm for} & \qquad m  =  0, & n \ne 0, \\
!latex           A_{\z,m,n} & \sim & r^{m+2} & g_{m,n}(\rho), & \qquad {\rm for} & \qquad m  =  1, &          \\
!latex           A_{\z,m,n} & \sim & r^{m  } & g_{m,n}(\rho), & \qquad {\rm for} & \qquad m  >  1. &
!latex       \end{array} \label{eq:originconstraintsaxisconstraint} \ee
!latex \item Note that, in the stellarator symmetric case, this eliminates $2N+1$ degrees-of-freedom in the representation for the vector potential.
!latex       This is precisely equal to the degrees-of-freedom in the representation for the coordinate axis,
!latex       $R = R_{0,n} \cos(-n\z)$ for $n=0,\dots,N$, and $Z = Z_{0,n} \sin(-n\z)$ for $n=1,\dots,N$.
!latex \item The boundary conditions are the same as in \Eqn{regularitygauge}, but now $\beta_{m,n}(0)=0$ for $m=1$.

!latex \newpage

!latex \end{enumerate} \subsection{interface boundary conditions: coordinate singularity} \begin{enumerate}

!latex \item The boundary conditions must accomodate the coordinate singularity by including the regularization factors:

!latex \be A_{\t,e,1}(-1) = 0 \; & \; \mbox{\rm and} \; & \;\; A_{\t,e,1}(+1) = \Delta\psi_{t,1} \\
!latex     A_{\t,e,i}(-1) = 0 \; & \; \mbox{\rm and} \; & \;\; A_{\t,e,i}(+1) = +m_i f_{o,i}     \\
!latex     A_{\t,o,i}(-1) = 0 \; & \; \mbox{\rm and} \; & \;\; A_{\t,o,i}(+1) = -m_i f_{e,i}     \\
!latex     A_{\z,e,1}(-1) = ? \; & \; \mbox{\rm and} \; & \;\; A_{\z,e,1}(+1) = \Delta\psi_{p,1} \\
!latex     A_{\z,e,i}(-1) = 0 \; & \; \mbox{\rm and} \; & \;\; A_{\z,e,i}(+1) = -n_i f_{o,i} \mbox{ for } m_i =   0,             n_i \ne 0 \\
!latex     A_{\z,e,i}(-1) = ? \; & \; \mbox{\rm and} \; & \;\; A_{\z,e,i}(+1) = -n_i f_{o,i} \mbox{ for } m_i \ne 0,                       \\
!latex     A_{\z,o,i}(-1) = 0 \; & \; \mbox{\rm and} \; & \;\; A_{\z,o,i}(+1) = +n_i f_{e,i} \mbox{ for } m_i =   0,             n_i \ne 0 \\
!latex     A_{\z,o,i}(-1) = ? \; & \; \mbox{\rm and} \; & \;\; A_{\z,o,i}(+1) = +n_i f_{e,i} \mbox{ for } m_i \ne 0,
!latex \ee

!latex       The `$?$' symbol indicates that there is no boundary condition, and $\Delta\psi_{p,1}$ is the magnetic flux linking the torus:
!latex       it must be specified, but it is otherwise irrelevant.

!latex \item The ``double'' boundary conditions may be enforced using as before.
!latex       Where there is only a boundary condition at the outer interface, i.e. for the $A_{\z,*,i}$ harmonics with $m_i = 0$, $n_i \ne 0$, 
!latex       the boundary conditions are enforced simply
!latex       \be f_{L  } & = & f(+1) - f_0 - f_1 - f_2 - \dots - f_{L-1}
!latex       \ee

!latex \item Using loose notation, the first derivatives with respect to the vector potential degress of freedom are
!latex       \be \frac{\partial}{\partial A_{\t,u,i,l}} & = & \frac{\partial}{\partial A_{\t,u,i,l}} - \frac{\partial}{\partial A_{\t,u,i,L_l}} \\
!latex           \frac{\partial}{\partial A_{\z,u,i,l}} & = & \frac{\partial}{\partial A_{\z,u,i,l}} - \frac{\partial}{\partial A_{\z,u,i,L_l}} 
!latex           \mbox{ for } m_i  =  0 \mbox{ and } n_i \ne 0 \\
!latex           \frac{\partial}{\partial A_{\z,u,i,l}} & = & \frac{\partial}{\partial A_{\z,u,i,l}} - \frac{\partial}{\partial A_{\z,u,i,L  }} 
!latex           \mbox{ for } m_i \ne 0 \mbox{ or  } n_i  =  0
!latex       \ee
!latex       where $L_l=L$ if $l$ is even and $L_l=L-1$ if $l$ is odd.

!latex \item The derivatives with respect to the surface potential are
!latex       \be \frac{\partial}{\partial f_{o,i}} & = & + \frac{1}{2} m_i \sum_{l=L-1}^{L} \frac{\partial}{\partial \Ate{i,l}} 
!latex                                                   - \frac{1}{2} n_i \sum_{l=L-1}^{L} \frac{\partial}{\partial \Aze{i,l}} \mbox{ for } m_i  =  0  \mbox{ and } n_i \ne 0 \\
!latex           \frac{\partial}{\partial f_{o,i}} & = & + \frac{1}{2} m_i \sum_{l=L-1}^{L} \frac{\partial}{\partial \Ate{i,l}} 
!latex                                                   - \frac{1}{2} n_i \sum_{l=L-1}^{L} \frac{\partial}{\partial \Aze{i,L}} \mbox{ for } m_i \ne 0  \mbox{  or } n_i  =  0 \\
!latex           \frac{\partial}{\partial f_{e,i}} & = & - \frac{1}{2} m_i \sum_{l=L-1}^{L} \frac{\partial}{\partial \Ato{i,l}} 
!latex                                                   + \frac{1}{2} n_i \sum_{l=L-1}^{L} \frac{\partial}{\partial \Azo{i,l}} \mbox{ for } m_i  =  0  \mbox{ and } n_i \ne 0 \\
!latex           \frac{\partial}{\partial f_{e,i}} & = & - \frac{1}{2} m_i \sum_{l=L-1}^{L} \frac{\partial}{\partial \Ato{i,l}} 
!latex                                                   + \frac{1}{2} n_i \sum_{l=L-1}^{L} \frac{\partial}{\partial \Azo{i,L}} \mbox{ for } m_i \ne 0  \mbox{  or } n_i  =  0
!latex       \ee

!latex \item The derivatives with respect to the enclosed fluxes are
!latex       \be \frac{\partial}{\partial \Delta \psi_t} & = & \frac{1}{2}\sum_{l=L-1}^{L} \frac{\partial}{\partial \Ate{1,l}} \\
!latex           \frac{\partial}{\partial \Delta \psi_p} & = & \frac{1}{2}\sum_{l=L-1}^{L} \frac{\partial}{\partial \Aze{1,L}}
!latex       \ee

!latex \newpage

!latex \end{enumerate} \subsection{integrands} \begin{enumerate}

!latex \item The integrands are

!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \sqrt g \; {\bf B} \cdot {\bf B} 

!latex & = & + &   & ( m_i \Azo{i,l} + n_i \Ato{i,l} ) & ( m_j \Azo{j,p} + n_j \Ato{j,p} ) & T_{l}        \; T_{p}        \; \bgss \; \cos\a_i \; \cos\a_j \\
!latex &   & - & 2 & ( m_i \Azo{i,l} + n_i \Ato{i,l} ) & ( m_j \Aze{j,p} + n_j \Ate{j,p} ) & T_{l}        \; T_{p}        \; \bgss \; \cos\a_i \; \sin\a_j \\
!latex &   & + &   & ( m_i \Aze{i,l} + n_i \Ate{i,l} ) & ( m_j \Aze{j,p} + n_j \Ate{j,p} ) & T_{l}        \; T_{p}        \; \bgss \; \sin\a_i \; \sin\a_j \\ \\

!latex &   & - & 2 & ( m_i \Azo{i,l} + n_i \Ato{i,l} ) &       \Aze{j,p}                   & T_{l}        \; T_{p}^\prime \; \bgst \; \cos\a_i \; \cos\a_j \\
!latex &   & - & 2 & ( m_i \Azo{i,l} + n_i \Ato{i,l} ) &       \Azo{j,p}                   & T_{l}        \; T_{p}^\prime \; \bgst \; \cos\a_i \; \sin\a_j \\
!latex &   & + & 2 & ( m_i \Aze{i,l} + n_i \Ate{i,l} ) &       \Aze{j,p}                   & T_{l}        \; T_{p}^\prime \; \bgst \; \sin\a_i \; \cos\a_j \\
!latex &   & + & 2 & ( m_i \Aze{i,l} + n_i \Ate{i,l} ) &       \Azo{j,p}                   & T_{l}        \; T_{p}^\prime \; \bgst \; \sin\a_i \; \sin\a_j \\  \\

!latex &   & + & 2 & ( m_i \Azo{i,l} + n_i \Ato{i,l} ) &       \Ate{j,p}                   & T_{l}        \; T_{p}^\prime \; \bgsz \; \cos\a_i \; \cos\a_j \\
!latex &   & + & 2 & ( m_i \Azo{i,l} + n_i \Ato{i,l} ) &       \Ato{j,p}                   & T_{l}        \; T_{p}^\prime \; \bgsz \; \cos\a_i \; \sin\a_j \\
!latex &   & - & 2 & ( m_i \Aze{i,l} + n_i \Ate{i,l} ) &       \Ate{j,p}                   & T_{l}        \; T_{p}^\prime \; \bgsz \; \sin\a_i \; \cos\a_j \\
!latex &   & - & 2 & ( m_i \Aze{i,l} + n_i \Ate{i,l} ) &       \Ato{j,p}                   & T_{l}        \; T_{p}^\prime \; \bgsz \; \sin\a_i \; \sin\a_j \\  \\

!latex &   & + &   &       \Aze{i,l}                   &       \Aze{j,p}                   & T_{l}^\prime \; T_{p}^\prime \; \bgtt \; \cos\a_i \; \cos\a_j \\
!latex &   & + & 2 &       \Aze{i,l}                   &       \Azo{j,p}                   & T_{l}^\prime \; T_{p}^\prime \; \bgtt \; \cos\a_i \; \sin\a_j \\
!latex &   & + &   &       \Azo{i,l}                   &       \Azo{j,p}                   & T_{l}^\prime \; T_{p}^\prime \; \bgtt \; \sin\a_i \; \sin\a_j \\  \\

!latex &   & - & 2 &       \Aze{i,l}                   &       \Ate{j,p}                   & T_{l}^\prime \; T_{p}^\prime \; \bgtz \; \cos\a_i \; \cos\a_j \\
!latex &   & - & 2 &       \Aze{i,l}                   &       \Ato{j,p}                   & T_{l}^\prime \; T_{p}^\prime \; \bgtz \; \cos\a_i \; \sin\a_j \\
!latex &   & - & 2 &       \Azo{i,l}                   &       \Ate{j,p}                   & T_{l}^\prime \; T_{p}^\prime \; \bgtz \; \sin\a_i \; \cos\a_j \\
!latex &   & - & 2 &       \Azo{i,l}                   &       \Ato{j,p}                   & T_{l}^\prime \; T_{p}^\prime \; \bgtz \; \sin\a_i \; \sin\a_j \\  \\

!latex &   & + &   &       \Ate{i,l}                   &       \Ate{j,p}                   & T_{l}^\prime \; T_{p}^\prime \; \bgzz \; \cos\a_i \; \cos\a_j \\
!latex &   & + & 2 &       \Ate{i,l}                   &       \Ato{j,p}                   & T_{l}^\prime \; T_{p}^\prime \; \bgzz \; \cos\a_i \; \sin\a_j \\
!latex &   & + &   &       \Ato{i,l}                   &       \Ato{j,p}                   & T_{l}^\prime \; T_{p}^\prime \; \bgzz \; \sin\a_i \; \sin\a_j           

!latex \end{array}        
!latex \ee
!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccc}
!latex \sqrt g \; {\bf A} \cdot {\bf B} 
!latex & = & - & \fbox{$ \Aze{i,l} \; \Ate{j,p} $} & \fbox{$ T_{l}^\prime \; T_{p} \; \cos\a_i \; \cos\a_j $}\\
!latex &   & - &         \Aze{i,l} \; \Ato{j,p}    &         T_{l}^\prime \; T_{p} \; \cos\a_i \; \sin\a_j   \\
!latex &   & - &         \Azo{i,l} \; \Ate{j,p}    &         T_{l}^\prime \; T_{p} \; \sin\a_i \; \cos\a_j   \\
!latex &   & - &         \Azo{i,l} \; \Ato{j,p}    &         T_{l}^\prime \; T_{p} \; \sin\a_i \; \sin\a_j   \\ \\
!latex &   & + & \fbox{$ \Ate{i,l} \; \Aze{j,p} $} & \fbox{$ T_{l}^\prime \; T_{p} \; \cos\a_i \; \cos\a_j $}\\
!latex &   & + &         \Ate{i,l} \; \Azo{j,p}    &         T_{l}^\prime \; T_{p} \; \cos\a_i \; \sin\a_j   \\
!latex &   & + &         \Ato{i,l} \; \Aze{j,p}    &         T_{l}^\prime \; T_{p} \; \sin\a_i \; \cos\a_j   \\
!latex &   & + &         \Ato{i,l} \; \Azo{j,p}    &         T_{l}^\prime \; T_{p} \; \sin\a_i \; \sin\a_j             
!latex \end{array}        
!latex \ee

!latex \newpage

!latex \end{enumerate} \subsection{first derivatives with respect to $\Ate{i,l}$ and $\Ato{i,l}$} \begin{enumerate}

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

!latex \end{enumerate} \subsection{first derivatives with respect to $\Aze{i,l}$ and $\Azo{i,l}$} \begin{enumerate}

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

!latex \end{enumerate} \subsection{second derivatives} \begin{enumerate}

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

!latex \end{enumerate} \subsection{second derivatives} \begin{enumerate}

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

!latex \end{enumerate} \subsection{second derivatives} \begin{enumerate}

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

!latex \end{enumerate} \subsection{second derivatives} \begin{enumerate}

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

!latex \newpage

!latex \end{enumerate} \subsection{reduction using FFTs and double angle formulae} \begin{enumerate}

!latex \item The metric elements may be represented as a Fourier series using fast Fourier transforms
!latex       \be \bar g_{\mu\nu}(\s,\t,\z) & = & \sum_k \bar g_{\mu\nu,e,k}(s) \; \cos(m_k\t-n_k\z) + \sum_k \bar g_{\mu\nu,o,k}(s) \; \sin(m_k\t-n_k\z) 
!latex       \ee

!latex \item The integrals over the angles can be computed using double angle formulae:
!latex       \be
!latex       \ooint g_{\mu\nu} \; \cos\a_i \; \cos\a_j & = & \frac{1}{2} \ooint g_{\mu\nu} ( \cos\a_{k_{ij}^-} + \cos\a_{k_{ij}^+} ) \\
!l!tex       \ooint g_{\mu\nu} \; \cos\a_i \; \sin\a_j & = & \ooint g_{\mu\nu} \; \sin\a_j \; \cos\a_i                               \\
!latex       \ooint g_{\mu\nu} \; \sin\a_i \; \cos\a_j & = & \frac{1}{2} \ooint g_{\mu\nu} ( \sin\a_{k_{ij}^-} + \sin\a_{k_{ij}^+} ) \\
!latex       \ooint g_{\mu\nu} \; \sin\a_i \; \sin\a_j & = & \frac{1}{2} \ooint g_{\mu\nu} ( \cos\a_{k_{ij}^-} - \cos\a_{k_{ij}^+} )
!latex       \ee
!latex       where $(m_{k_{ij}^-},n_{k_{ij}^-})=(m_i-m_j,n_i-n_j)$ and $(m_{k_{ij}^+},n_{k_{ij}^+})=(m_i+m_j,n_i+n_j)$.

!latex \newpage

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

!latex \newpage

!latex \end{enumerate} \subsection{spectral condensation} \begin{enumerate}

!latex \item The spectral width is defined
!latex       \be M = \frac{1}{2}\sum_i (m_i^p+n_i^q) \left(R_{e,i}^2+R_{o,i}^2+Z_{e,i}^2+Z_{o,i}^2\right).
!latex       \ee
!latex       where $p\equiv$ \verb+pwidth+, $q\equiv$ \verb+qwidth+ are positive integers given on input, and $m_i^p=0$ for $m_i=0$, $n_i^q=0$ for $n_i=0$.

!latex \item We seek to extremize the spectral width without changing the geometry of the interface.
!latex       We restrict attention to tangential, poloidal variations, i.e. variations of the form
!latex       \be \delta R &=& R_\t \; \delta u , \label{eq:deltaR}\\
!latex           \delta Z &=& Z_\t \; \delta u . \label{eq:deltaZ}
!latex       \ee
!latex       where $\delta u = \sum_j u_{e,j} \cos \alpha_j + \sum_j u_{o,j} \sin \alpha_j$.

!latex \item From \Eqn{deltaR} and \Eqn{deltaZ}, the variations in the Fourier harmonics of $R$ and $Z$ are given by
!latex       \be \delta R_{e,i} &=& \ooint R_\t \; \delta u \; \cos \alpha_i, \\
!latex           \delta R_{o,i} &=& \ooint R_\t \; \delta u \; \sin \alpha_i, \\
!latex           \delta Z_{e,i} &=& \ooint Z_\t \; \delta u \; \cos \alpha_i, \\
!latex           \delta Z_{o,i} &=& \ooint Z_\t \; \delta u \; \sin \alpha_i.
!latex       \ee

!latex \item The first variation in $M$ as
!latex       \be \delta M = \ooint \left(R_\t X +  Z_\t Y \right) \delta u,
!latex       \ee
!latex       where $X = \sum_i (m_i^p+n_i^q)(R_{e,i} \cos\alpha_i+R_{o,i} \sin\alpha_i)$ and $Y = \sum_i (m_i^p+n_i^q)(Z_{e,i} \cos\alpha_i+Z_{o,i} \sin\alpha_i)$.

!latex \item The condition that $\delta M = 0$ for arbitrary $\delta u$ is 
!latex       \be I \equiv R_\t \; X + Z_\t \; Y = 0. \ee

!latex \end{enumerate} \subsubsection{poloidal angle origin} \begin{enumerate}

!latex \item In the non-stellarator-symmetric case, it seems that the poloidal angle origin, $\t=0$, on the $\z=0$ plane is not constrained, and so
!latex       an additional constraint must be introduced.
!latex \item We seek to minimize the total length of the curve between $(R_0(0,0),Z_0(0,0))$ and $(R_N(0,0),Z_N(0,0))$,
!latex       where $R_l(\t,\z)$ and $Z_l(\t,\z)$ describe the geometry of the $l$-th interface.

!latex \item Define the length-squared of the poloidal angle origin curve as
!latex       \be L \equiv \sum_{l=1}^{N} \frac{1}{2}\left\{\left[R_{l}(0,0)-R_{l-1}(0,0)\right]^2+\left[Z_{l}(0,0)-Z_{l-1}(0,0)\right]^2\right\},
!latex       \ee
!latex       and the additional constraint is
!latex       \be \frac{\partial L}{\partial \t_l} \equiv \frac{\partial R_l}{\partial \t}\left(-R_{l-1}+2R_{l}-R_{l+1}\right) 
!latex                                                 + \frac{\partial Z_l}{\partial \t}\left(-Z_{l-1}+2Z_{l}-Z_{l+1}\right) = 0.
!latex       \ee

!l!tex \item The derivatives of $M$ with respect to the $u_k$ are given
!l!tex       \be \frac{\partial M}{\partial u_{e,j}} & = & \ooint (R_\t X + Z_\t Y) \cos\alpha_j \\
!l!tex           \frac{\partial M}{\partial u_{o,j}} & = & \ooint (R_\t X + Z_\t Y) \sin\alpha_j
!l!tex       \ee
!l!tex       These quantities are provided by an fast Fourier transform of $I$.

!l!tex \item The derivatives of the spectral constraints, $I =R_\t X + Z_\t Y$,  are derived using $\sin(\alpha+\beta)=\left[\sin(\alpha+\beta)+\sin(\alpha-\beta)\right]/2$
!l!tex       to give
!l!tex       \be \frac{\partial I}{\partial R_j} &=& \frac{\partial R_\t}{\partial R_j} X + R_\t \frac{\partial X}{\partial R_j} \nonumber \\
!l!tex                                           &=& \frac{1}{2} \sum_k \left[ - ( m_j \lambda_k + m_k \lambda_j ) R_k \sin (\alpha_j + \alpha_k )
!l!tex                                                                         - ( m_j \lambda_k - m_k \lambda_j ) R_k \sin (\alpha_j - \alpha_k ) \right]\\
!l!tex           \frac{\partial I}{\partial Z_j} &=& \frac{\partial Z_\t}{\partial R_j} Y + Z_\t \frac{\partial Y}{\partial R_j} \nonumber \\
!l!tex                                           &=& \frac{1}{2} \sum_k \left[ + ( m_j \lambda_k + m_k \lambda_j ) Z_k \sin (\alpha_j + \alpha_k )
!l!tex                                                                         - ( m_j \lambda_k - m_k \lambda_j ) Z_k \sin (\alpha_j - \alpha_k ) \right].
!l!tex       \ee

!latex \newpage

!latex \end{enumerate} \subsection{energy functional gradient and derivatives} \begin{enumerate} 
     
!latex \item Consider 
!latex       \be \delta F = - \int_{\cal I} [[p+B^2/2]] \; \boldxi \cdot {\bf dS} = - \int \!\!\!\!\! \int [[p+B^2/2]] \; \boldxi \cdot {\bf e}_\t \times {\bf e}_\z \; d\t d\z.
!latex       \ee
!latex       \begin{enumerate}     

!latex \item \verb+Igeometry.eq.1+ Cartesian geometry: 
!latex       \be \frac{\partial F}{\partial R_j} = [[p+B^2/2]]_j\ee

!latex \item \verb+Igeometry.eq.2+ cylindrical:
!latex       \be \frac{\partial F}{\partial R_j} &=& - \left( [[ p + B^2/2 ]] \; R \right)_j
!latex       \ee

!latex \item \verb+Igeometry.eq.3+ toroidal geometry:
!latex       \be \frac{\partial F}{\partial R_j} &=& - \left( [[ p + B^2/2 ]] \; R \; Z_\t \right)_j\\
!latex           \frac{\partial F}{\partial Z_j} &=& + \left( [[ p + B^2/2 ]] \; R \; R_\t \right)_j.
!latex       \ee  

!latex \item \verb+Igeometry.eq.4+ cylindrical geometry;
!latex       \be \frac{\partial F}{\partial R_j} &=& - \left( [[ p + B^2/2 ]] \; R \; Z_\t \right)_j\\
!latex           \frac{\partial F}{\partial Z_j} &=& + \left( [[ p + B^2/2 ]] \; R \; R_\t \right)_j.
!latex       \ee  
    
!latex       \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine manual
 
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero

  use numerical, only : 

  use fileunits, only : ounit

  use inputlist, only : Wmanual

  use cputiming, only : Tmanual

  use allglobal, only : myid, cpus
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS

  BEGIN(manual)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(manual)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine manual

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
