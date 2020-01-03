!> \file manual.f90
!> \brief Code development issues and future physics applications.
!> \see \ref grp_documentation


!> \page grp_documentation Manual / Documentation
!>
!> \section grp_polFluxRotTransform Poloidal flux and rotational transform
!>
!> Given the canonical integrable form, 
!> \f${\bf A} = \psi \nabla \theta - \chi(\psi) \nabla \zeta\f$, 
!> we can derive 
!> \f${\bf B} = \nabla \psi \times \nabla \theta +\nabla \zeta \times \nabla \psi \;\chi^\prime\f$.
!> The poloidal flux is given by 
!> \f{eqnarray}{ \Psi_p = \int \!\! \int {\bf B} \cdot {\bf e}_\zeta \times{\bf e}_\psi \; d\zeta d\psi = 2 \pi \int \chi^\prime d\psi.
!> \f}
!> The rotational-transform is 
!> \f{eqnarray}{ {{\,\iota\!\!\!}-} = \frac{ {\bf B} \cdot \nabla \theta}{ {\bf B} \cdot \nabla \zeta} = \chi^\prime.
!> \f}
!> The rotational-transform has the same sign as the poloidal flux.
!> 
!> The SPEC representation for the magnetic vector potential is
!> \f{eqnarray}{ {\bf A} = A_\theta \nabla \theta + A_\zeta \nabla \zeta,
!> \f}
!> where we can see that \f$A_\zeta = - \chi\f$.
!> The poloidal flux is 
!> \f{eqnarray}{ \int {\bf B}\cdot d{\bf s} = \oint A_\zeta d\zeta.
!> \f}
!> It would seem that the rotational-transform has opposite sign to \f$A_\zeta\f$.
!> To be honest, I am a little confused regarding the sign.
!> 
!> \section grp_outline Outline
!> 
!> This document is intended to organise the different potentially valuable improvements to the SPEC code, which could make it more robust, faster, and increase its capabilities. 
!> 
!> The document is divided in two categories:
!> 
!> \ref sec_NumericalImprovements : independent improvements that are of numerical importance but have no added physics value *per se*, although they may allow new or better physics investigations. 
!> 
!> \ref sec_PhysicsApplications : research topics that could be addressed with the code, either in its present form or after the completion of one or more topics listed in \ref sec_NumericalImprovements .
!> 
!> \section sec_NumericalImprovements Numerical Improvements
!> 
!> \subsection sec_gcc Compile code with GCC for error checking
!> 
!> Has been implemented in Makefile for most platforms. Checks against Intel version show small differences on the order of \f$10^{-15}\f$ relative deviation,
!> which are likely due so slighly different optimization strategies.
!> 
!> \subsection sec_gprof Profile code with gprof to find inefficient lines of code
!>
!> \subsection sec_valgrind Run code with Valgrind to identify memory leaks
!>
!> \subsection sec_denag De-NAG-ification
!> 
!> Compilation of SPEC does not rely on NAG anymore; some functionality (e.g. SQP in ma02aa.f90) might need replacements for the NAG routines to be re-enabled.
!>
!> \subsection sec_spectralConstraints Revision of spectral-constraints
!> 
!> This is bit of a mess.
!> All the mathematics is standard, and all that is required is for someone to calmly go through lots of algebra.
!> This task should be high priority, as SRH suspects that the spectral constraints as presently enforced result in an ill-conditioned force vector,
!> which means that the code is overly sensitive to the initial guess and does not converge robustly.
!> Potential speed improvements are tremendous.
!>
!> \subsection sec_torAngle Extension to arbitrary toroidal angle
!> 
!> This can further reduce the required Fourier resolution, and so this can reduce the computation.
!> SRH is particularly interested in this as it will allow for exotic configurations (knots, figure-8, etc.) that cannot presently be computed.
!>
!> \subsection sec_metricSymmetry Exploit symmetry of the metric
!> 
!> This is easy, but somewhat tedious.
!> Take a look at ma00aa() to see what is required.
!> Potential speed improvement is considerable.
!>
!> \subsection sec_beltrami symmetry of "local" Beltrami matrices
!> 
!> This is easy. Take a look at matrix(), which constructs the Beltrami matrices, and mp00ac(), which performs the inversion. 
!> Potential speed improvement is considerable.
!>
!> \subsection sec_tridiagnonal Exploit block tri-diagonal structure of "global" linearized force balance matrix
!> 
!> This requires an efficient subroutine.
!> SRH believes that Hirshman constructed such a routine (Hirshman et al. (2010) \cite y2010_hirshman).
!> The potential speed improvement is tremendous.
!> See newton() for where the tri-diagonal, linearized force-balance matrix is inverted.
!>
!> \subsection sec_heliciy_constraint Enforce Helicity constraint
!> 
!> This will allow investigation of different, arguably more-physical classes of equilibria. 
!> See ma02aa() .
!>
!> \subsection sec_tests Establish test-cases
!> 
!> A suite of test cases should be constructed, with different geometries etc., that run fast, and that can be benchmarked to machine precision.
!> In the InputFiles/TestCases directory, some input files for SPEC are available for this purpose.
!> One should write routines which execute these input files and compare the output data against a publicy-available set of output files
!> to check SPEC before a new release is made.
!>
!> \subsection sec_freeb Verify free-boundary
!> 
!> This is almost complete. The corresponding publication is being written.
!> The virtual casing routines need to be investigated and made more efficient.
!> The virtual casing routine in slab geometry needs revision (because of an integral over an infinite domain).
!>
!> \subsection sec_toroidal_current Enforcement of toroidal current profile
!> 
!> Adjust \f$\mu\f$'s, fluxes and/or rotational transform to obtain desired current profile (without singular currents).
!> This is implemented and needs to be merged into the master branch.
!> An additional routine is required to iterate on the helicity multipliers etc. as required *after* the local Beltrami fields have been calculated
!> and *before* the global force balance iterations proceed.
!>
!> \subsection sec_hessian Interpret eigenvectors and eigenvalues of Hessian
!> 
!> This is already completed: see hesian(). However, this actually computes the force gradient matrix.
!> For toroidal geometry there is a complication; namely that the hessian matrix includes the derivatives of the spectral constraints.
!> For Cartesian geometry, it is ready to go. 
!> SRH will begin writing a paper on the stability of slab MRxMHD equilibria.
!>
!> \section sec_PhysicsApplications Physics Applications
!>
!> \subsection sec_hires Calculate high-resolution equilibria, e.g. W7-X
!> 
!> requires: \ref sec_metricSymmetry , \ref sec_beltrami , and other improvements that can make the code faster at high Fourier resolution
!> 
!> \subsection sec_calc_consvHelicity Calculate equilibria by conserving helicity and fluxes
!> 
!> Applications to saturated island studies, sawteeth, etc.
!> requires: \ref sec_calc_consvHelicity
!> 
!> \subsection sec_calc_freeb Calculate free-boundary stellarator equilibria
!> 
!> to predict scrape-off-layer (SOL) topologies and \f$\beta\f$-limits.
!> requires: \ref sec_freeb
!> Mostly complete.
!> 
!> \subsection sec_eval_stability Evaluate stability of MRxMHD equilibria
!> 
!> perhaps starting from simplest system (slab tearing).
!> requires: \ref sec_hessian
!> 
!> \section grp_coord_singularity Revision of coordinate singularity: axisymmetric; polar coordinates
!>
!> <ul>
!> <li> Consider a general, magnetic vector potential given in Cartesian coordinates,
!>       \f{eqnarray}{{\bf A} = A_x \nabla x + A_y \nabla y +  A_z \nabla z + \nabla g \label{eq:CartesianVectorPotential_manual}
!>       \f}
!>       where \f$A_x\f$, \f$A_y\f$, \f$A_z\f$, and the as-yet-arbitrary gauge function, \f$g\f$, are regular at \f$(x,y)=(0,0)\f$,
!>       i.e. they can be expanded as a Taylor series, e.g.
!>       \f{eqnarray}{A_x = \sum_{i,j} \alpha_{i,j} x^i y^j, \qquad
!>           A_y = \sum_{i,j}  \beta_{i,j} x^i y^j, \qquad
!>           A_z = \sum_{i,j} \gamma_{i,j} x^i y^j, \qquad
!>            g  = \sum_{i,j} \delta_{i,j} x^i y^j, \label{eq:Taylorexpansion_manual}
!>       \f}
!>       for small \f$x\f$ and small \f$y\f$. </li>
!> <li> Note that we have restricted attention to the "axisymmetric" case, as there is no dependence on \f$z\f$. </li>
!> <li> The "polar" coordinate transformation,
!>       \f{eqnarray}{x &=& r \cos \theta, \nonumber \\ y &=& r \sin \theta, \\ z &=& \zeta, \nonumber
!>       \f}
!>       induces the vector transformation
!>       \f{eqnarray}{\begin{array}{ccccccccccccccccccccccccccccccc} \nabla x & = & \cos \theta \; \nabla r & - & r \sin \theta \; \nabla \theta & &, \\
!>                                                                   \nabla y & = & \sin \theta \; \nabla r & + & r \cos \theta \; \nabla \theta & &, \\
!>                                                                   \nabla z & = &                         &   &                                & \nabla \zeta&.
!>       \end{array} \f} </li>
!> <li> By repeated applications of the double-angle formula, the expressions for \f$A_x\f$, \f$A_y\f$ and \f$g\f$ can be cast as functions of \f$(r,\theta)\f$,
!>       \f{eqnarray}{A_x & = & \sum_m r^m [ a_{m,0} + a_{m,1} \; r^2 + a_{m,2} \; r^4 + ... ] \sin(m\theta), \\
!>                    A_y & = & \sum_m r^m [ b_{m,0} + b_{m,1} \; r^2 + b_{m,2} \; r^4 + ... ] \cos(m\theta), \\
!>                    A_z & = & \sum_m r^m [ c_{m,0} + c_{m,1} \; r^2 + c_{m,2} \; r^4 + ... ] \cos(m\theta), \label{eq:regularAz_manual} \\
!>                     g  & = & \sum_m r^m [ g_{m,0} + g_{m,1} \; r^2 + g_{m,2} \; r^4 + ... ] \sin(m\theta), 
!>       \f}
!>       where attention is restricted to stellarator symmetric geometry, but similar expressions hold for the non-stellarator symmetric terms. </li>
!> <li> Collecting these expressions, the vector potential can be expressed
!>       \f{eqnarray}{{\bf A} = A_r \nabla r + A_\theta \nabla \theta + A_\zeta \nabla \zeta + \partial_r g \; \nabla r + \partial_\theta g \; \nabla \theta,
!>       \f}
!>       where
!>       \f{eqnarray}{\begin{array}{ccccccclcccccccccccccccccccccccccccccccccccccccccccccccccccc}
!>       A_r & = & r^0 & [ & ( &           &   & b_{0,0}   & + &   g_{1,0} & ) & + & (...) r^2 + (...)r^4 + ... & ] & \sin   \theta \\
!>           & + & r^1 & [ & ( & a_{1,0}/2 & + & b_{1,0}/2 & + & 2 g_{2,0} & ) & + & (...) r^2 + (...)r^4 + ... & ] & \sin 2 \theta \\
!>           & + & r^2 & [ & ( & a_{2,0}/2 & + & b_{2,0}/2 & + & 3 g_{3,0} & ) & + & (...) r^2 + (...)r^4 + ... & ] & \sin 3 \theta \\
!>           & + & ...
!>       \end{array} \f}
!>       (Note: Mathematica was used to perform the algebraic manipulations,
!>       and the relevant notebook was included as part of the SPEC \c CVS repository.) </li>
!> <li> There is precisely enough gauge freedom so that we may choose \f$A_r=0\f$.
!>       For example, the choice
!>       \f{eqnarray}{\begin{array}{cccccclccccccccccccccccccccccccccccccccc}
!>           g_{1,0} & = & - &   &           &   & b_{0,0}    &   &   &  &, \\
!>           g_{2,0} & = & - & ( & a_{1,0}/2 & + & b_{1,0}/2  & ) & / & 2&, \\
!>           g_{3,0} & = & - & ( & a_{2,0}/2 & + & b_{2,0}/2  & ) & / & 3&, \\
!>           ...   & = &   & ...
!>       \end{array} \f}
!>       eliminates the lowest order \f$r\f$ dependence in each harmonic. </li>
!> <li> By working through the algebra (again, using Mathematica) the expressions for \f$A_\theta\f$ and \f$A_\zeta\f$ become
!>       \f{eqnarray}{A_\theta &=& r^2     f_0(\rho) + r^3 f_1(\rho) \cos(\theta) + r^4 f_2(\rho) \cos(2\theta) + r^5 f_3(\rho) \cos(3\theta) + ... \label{eq:nearoriginAt_manual} \\
!>                    A_\zeta  &=& \;\;\;\,g_0(\rho) + r^1 g_1(\rho) \cos(\theta) + r^2 g_2(\rho) \cos(2\theta) + r^3 g_3(\rho) \cos(3\theta) + ... \label{eq:nearoriginAz_manual}
!>       \f}
!>       where \f$\rho\equiv r^2\f$ and the \f$f_m(\rho)\f$ and \f$g_m(\rho)\f$ are abitrary polynomials in \f$\rho\f$.
!>       [The expression for \f$A_\zeta\f$ is unchanged from Eqn.\f$(\ref{eq:regularAz_manual})\f$.] </li>
!> </ul>
!>
!> \subsection sec_generally somewhat generally, ...
!>
!> <ul>
!> <li> For stellarator-symmetric configurations,
!>       \f{eqnarray}{{\bf A} = \sum_{m,n} A_{\theta,m,n} \cos(m\theta-n\zeta) \nabla \theta + \sum_{m,n} A_{\zeta,m,n} \cos(m\theta-n\zeta) \nabla \zeta,
!>       \f}
!>       where now the dependence on \f$\zeta\f$ is included, and the angles are arbitrary. </li>
!> <li> The near-origin behaviour of \f$A_\theta\f$ and \f$A_\zeta\f$ given in Eqn.\f$(\ref{eq:nearoriginAt_manual})\f$ and Eqn.\f$(\ref{eq:nearoriginAz_manual})\f$ are flippantly generalized to
!>       \f{eqnarray}{
!>           A_{\theta,m,n} & = & r^{m+2}          f_{m,n}(\rho), \label{eq:Atmn_manual} \\
!>           A_{\zeta ,m,n} & = & r^{m  } \;\;\;\; g_{m,n}(\rho), \label{eq:Azmn_manual}
!>       \f}
!>       where the \f$f_{m,n}(\rho)\f$ and \f$g_{m,n}(\rho)\f$ are arbitrary polynomials in \f$\rho\f$. </li>
!> <li> Additional gauge freedom can be exploited: including an additional gauge term \f$\nabla h\f$ where \f$h\f$ only depends on \f$\zeta\f$, e.g.
!>       \f{eqnarray}{h(\zeta) = h_{0,0} \, \zeta + \sum h_{0,n} \sin( - n\zeta),\f} 
!>       does not change the magnetic field and does not change any of the above discussion. </li>
!> <li> The representation for the \f$A_{\theta,m,n}\f$ does not change, but we must clarify that Eqn.\f$(\ref{eq:Azmn_manual})\f$ holds for only the \f$m\ne0\f$ harmonics:
!>       \f{eqnarray}{A_{\zeta,m,n} & = & r^{m} \;\;\;\;g_{m,n}(\rho), \;\;\;\mathrm{for}\;\;\; m \ne 0.
!>       \f} </li>
!> <li> For the \f$m=0\f$, \f$n\ne0\f$ harmonics of \f$A_\zeta\f$, including the additional gauge gives \f$A_{\zeta,0,n} = g_{0,n}(\rho) + n \, h_{0,n}\f$.
!>       Recall that \f$g_{0,n}(\rho) = g_{0,n,0} + g_{0,n,1} \rho + g_{0,n,2} \rho^2 + ...\f$, and we can choose \f$h_{0,n} = - g_{0,n,0} / n\f$ to obtain
!>       \f{eqnarray}{A_{\zeta,m,n} & = & r^{m  } \;\;\;g_{m,n}(\rho), \;\;\;\mathrm{for}\;\; m=0, n\ne0, \;\;\mathrm{with}\;\; g_{m,n}(0)=0.
!>       \f} </li>
!> <li> For the \f$m=0\f$, \f$n=0\f$  harmonic of \f$A_\zeta\f$, we have \f$A_{\zeta,0,0} = g_{0,0}(\rho) + h_{0,0}\f$.
!>       Similarly, choose \f$h_{0,0} = - g_{0,n,0}\f$ to obtain
!>       \f{eqnarray}{A_{\zeta,m,m} & = & r^{m  } \;\;\;g_{m,n}(\rho), \;\;\;\mathrm{for}\;\; m=0, n=0, \;\;\mathrm{with}\;\; g_{m,n}(0)=0.
!>       \f} </li>
!> <li> To simplify the algorithmic implementation of these conditions,
!>       we shall introduce a "regularization" factor, \f$\rho^{m/2} = r^m\f$. </li>
!> <li> Note that the representation for \f$A_{\theta,m,n}\f$ given in Eqn.\f$(\ref{eq:Atmn_manual})\f$,
!>       with an arbitrary polynomial \f$f_{m,n}(\rho) = f_{m,n,0} + f_{m,n,1}\rho + f_{m,n,2}\rho^2 + ...\f$, 
!>       is equivalent to \f$A_{\theta,m,n} = \rho^{m/2} \alpha_{m,n}(\rho)\f$ where \f$\alpha_{m,n}(\rho)\f$ is an arbitrary polynomial
!>       with the constraint \f$\alpha_{m,n}(0)=0\f$. </li>
!> <li> We can write the vector potential as
!>       \f{eqnarray}{A_{\theta,m,n} & = & \rho^{m/2} \alpha_{m,n}(\rho), \;\;\; \mathrm{with}\;\; \alpha_{m,n}(0)=0 \;\;\mathrm{for\;all}\;\; (m,n),\\
!>                    A_{\zeta ,m,n} & = & \rho^{m/2}  \beta_{m,n}(\rho), \;\;\; \mathrm{with}\;\;  \beta_{m,n}(0)=0 \;\;\mathrm{for}\;\; m=0.
!>       \f} </li>
!> </ul>
!>
!> \subsection sec_nonstellsym non-stellarator symmetric terms
!>
!> <ul>
!> <li> Just guessing, for the non-stellarator-symmetric configurations,
!>       \f{eqnarray}{A_{\theta,m,n} & = & \rho^{m/2} \alpha_{m,n}(\rho), \;\;\; \mathrm{with}\;\; \alpha_{m,n}(0)=0 \;\;\mathrm{for\;all}\;\; (m,n),\\
!>                    A_{\zeta ,m,n} & = & \rho^{m/2}  \beta_{m,n}(\rho), \;\;\; \mathrm{with}\;\;  \beta_{m,n}(0)=0 \;\;\mathrm{for}\;\;     m=0.
!>       \f} </li>
!> </ul>
!> 
!>

!l tex **constraints on Fourier harmonics}

!l tex \begin{enumerate}
!l tex <li> The boundary conditions must accomodate the coordinate singularity by including the regularization factors:
!l tex \f{eqnarray}{A_{\t,e,1}(-1) = 0 \; & \; \mbox{\rm and} \; & \;\; A_{\t,e,1}(+1) = \Delta\psi_{t,1} \\
!l tex     A_{\t,e,i}(-1) = 0 \; & \; \mbox{\rm and} \; & \;\; A_{\t,e,i}(+1) = +m_i f_{o,i}     \\
!l tex     A_{\t,o,i}(-1) = 0 \; & \; \mbox{\rm and} \; & \;\; A_{\t,o,i}(+1) = -m_i f_{e,i}     \\
!l tex     A_{\z,e,1}(-1) = ? \; & \; \mbox{\rm and} \; & \;\; A_{\z,e,1}(+1) = \Delta\psi_{p,1} \\
!l tex     A_{\z,e,i}(-1) = 0 \; & \; \mbox{\rm and} \; & \;\; A_{\z,e,i}(+1) = -n_i f_{o,i} \mbox{ for } m_i =   0,             n_i \ne 0 \\
!l tex     A_{\z,e,i}(-1) = ? \; & \; \mbox{\rm and} \; & \;\; A_{\z,e,i}(+1) = -n_i f_{o,i} \mbox{ for } m_i \ne 0,                       \\
!l tex     A_{\z,o,i}(-1) = 0 \; & \; \mbox{\rm and} \; & \;\; A_{\z,o,i}(+1) = +n_i f_{e,i} \mbox{ for } m_i =   0,             n_i \ne 0 \\
!l tex     A_{\z,o,i}(-1) = ? \; & \; \mbox{\rm and} \; & \;\; A_{\z,o,i}(+1) = +n_i f_{e,i} \mbox{ for } m_i \ne 0,
!l tex \f}
!l tex       The `$?$' symbol indicates that there is no boundary condition, and $\Delta\psi_{p,1}$ is the magnetic flux linking the torus:
!l tex       it must be specified, but it is otherwise irrelevant.
!l tex </ul>

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! subroutine manual
 
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!  use constants, only : zero
!
!  use numerical, only : 
!
!  use fileunits, only : ounit
!
!  use inputlist, only : Wmanual
!
!  use cputiming, only : Tmanual
!
!  use allglobal, only : myid, cpus
!  
!!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!  
!  LOCALS
!
!  BEGIN(manual)
!
!!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!
! ! this "routine" is purely for documentation; 08 Feb 16;
!  
!!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!
!  RETURN(manual)
!
!!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!
!end subroutine manual

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
