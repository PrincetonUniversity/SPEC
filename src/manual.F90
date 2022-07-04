!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (documentation) ! Code development schedule.

!latex \briefly{code development issues and future physics applications.}

!l tex \calledby{\link{}}
!l tex \calls{\link{}}

!latex \tableofcontents

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{poloidal flux and rotational transform}

!latex Given the canonical integrable form,
!latex ${\bf A} = \psi \nabla \t - \chi(\psi) \nabla \z$,
!latex we can derive
!latex ${\bf B} = \nabla \psi \times \nabla \t +\nabla \zeta \times \nabla \psi \;\chi^\prime$.
!latex The poloidal flux is given by
!latex \begin{eqnarray} \Psi_p = \int \!\! \int {\bf B} \cdot {\bf e}_\zeta \times{\bf e}_\psi \; d\zeta d\psi = 2 \pi \int \chi^\prime d\psi.
!latex \end{eqnarray}
!latex The rotational-transform is
!latex \begin{eqnarray} \iotabar = \frac{ {\bf B} \cdot \nabla \t}{ {\bf B} \cdot \nabla \zeta} = \chi^\prime.
!latex \end{eqnarray}
!latex The rotational-transform is the same sign as the poloidal flux.
!latex
!latex The SPEC representation for the magnetic vector potential is
!latex \begin{eqnarray} {\bf A} = A_\t \nabla \t + A_\z \nabla \z,
!latex \end{eqnarray}
!latex where we can see that $A_\z = - \chi$.
!latex The poloidal flux is
!latex \begin{eqnarray} \int {\bf B}\cdot d{\bf s} = \oint A_\z d\z.
!latex \end{eqnarray}
!latex It would seem that the rotational-transform has opposite sign to $A_\z$.
!latex To be honest, I am a little confused regarding the sign.

!latex \subsection{Outline}

!latex This document is intended to organise the different potentially valuable improvements to the SPEC code, which could make it more robust, faster, and increase its capabilities.
!latex The document is divided in two categories.
!latex In \Sec{NumericalImprovements}, Numerical Improvements: independent improvements that are of numerical importance but have no added physics value \emph{per se}, although they may allow new or better physics investigations.
!latex In \Sec{PhysicsApplications}, research topics that could be addressed with the code, either in its present form or after the completion of one or more topics listed in \Sec{NumericalImprovements}.

!latex \subsection{Numerical Improvements} \label{sec:NumericalImprovements}

!latex \subsubsection{Downloading \SPEC to \type{Git}}
!latex This should be straight forward and easy.
!latex SRH will work on this. Estimated date of completion 06/2017.

!latex \subsubsection{Compile code with \type{GCC} for error checking} \label{sec:recompile}

!latex \subsubsection{Profile code with \type{gprof} to find inefficient lines of code} \label{sec:gprof}

!latex \subsubsection{Run code with \type{Valgrind} to identify memory leaks} \label{sec:valgrind}

!latex \subsubsection{De-NAG-ification}
!latex This will take time . . .
!latex It is not, at present, crucial, as we all have access to NAG.

!latex \subsubsection{Revision of spectral-constraints}
!latex This is bit of a mess.
!latex All the mathematics is standard, and all that is required is for someone to calmly go through lots of algebra.
!latex This task should be high priority, as SRH suspects that the spectral constraints as presently enforced result in an ill-conditioned force vector,
!latex which means that the code is overly sensitive to the initial guess and does not converge robustly.
!latex Potential speed improvements are tremendous.

!latex \subsubsection{Extension to arbitrary toroidal angle}
!latex This can further reduce the required Fourier resolution, and so this can reduce the computation.
!latex SRH is particularly interested in this as it will allow for exotic configurations (knots, figure-8, etc.) that cannot presently be computed.

!latex \subsubsection{Exploit symmetry of the metric} \label{sec:metric}
!latex This is easy, but somewhat tedious.
!latex Take a look at \link{ma00aa} to see what is required.
!latex Potential speed improvement is considerable.

!latex \subsubsection{Exploit symmetry of ``local'' Beltrami matrices} \label{sec:beltrami}
!latex This is easy. Take a look at \link{matrix}, which constructs the Beltrami matrices, and \link{mp00ac}, which performs the inversion.
!latex Potential speed improvement is considerable.

!latex \subsubsection{Exploit block tri-diagonal structure of ``global'' linearized force balance matrix}
!latex This requires an efficient subroutine.
!latex SRH believes that Hirshman constructed such a routine
!latex [\paper{S.P. Hirshman {\em et al.}}{S.P. Hirshman, K.S. Purumalla {\em et al.}}{10.1016/j.jcp.2010.04.049}{J. Comput. Phys.}{229}{6392}{2010}].
!latex The potential speed improvement is tremendous.
!latex See \link{newton} for where the tri-diagonal, linearized force-balance matrix is inverted.

!latex \subsubsection{Enforce Helicity constraint} \label{sec:L2}
!latex This will allow investigation of different, arguably more-physical classes of equilibria.
!latex See \link{ma02aa}

!latex \subsubsection{Establish test-cases} \label{sec:testcase}
!latex A suite of test cases should be constructed, with different geometries etc., that run fast, and that can be benchmarked to machine precision.

!latex \subsubsection{Verify free-boundary} \label{sec:freeb}
!latex This is almost complete.
!latex The virtual casing routines need to be investigated and made more efficient.
!latex The virtual casing routine in slab geometry needs revision (because of an integral over an infinite domain).

!latex \subsubsection{Enforcement current ``profile''}
!latex i.e., adjust $\mu$'s, fluxes and/or rotational transform to obtain desired current profile (without singular currents).
!latex This should be a reasonably simple task.
!latex An additional routine is required to iterate on the helicity multipliers etc. as required {\em after} the local Beltrami fields have been calculated
!latex and {\em before} the global force balance iterations proceed.

!latex \subsubsection{Interpret eigenvectors and eigenvalues of Hessian} \label{sec:stability}
!latex This is already completed: see \link{hesian}.
!latex For toroidal geometry there is a complication; namely that the hessian matrix includes the derivatives of the spectral constraints.
!latex For Cartesian geometry, it is ready to go.
!latex SRH will begin writing a paper on the stability of slab MRxMHD equilibria.

!latex \subsection{Physics Applications} \label{sec:PhysicsApplications}

!latex \subsubsection{Calculate high-resolution equilibria} e.g. W7-X.
!latex requires: \Sec{metric}, \Sec{beltrami}, and other improvements that can make the code faster at high Fourier resolution
!latex \subsubsection{Calculate equilibria by conserving helicity and fluxes} Applications to saturated island studies, sawteeth, etc.
!latex requires: \Sec{L2}
!latex \subsubsection{Calculate free-boundary stellarator equilibria} to predict scrape-off-layer (SOL) topologies and $\beta$-limits.
!latex requires: \Sec{freeb}
!latex \subsubsection{Evaluate stability of MRxMHD equilibria} perhaps starting from simplest system (slab tearing).
!latex requires: \Sec{stability}

!latex \subsection{Revision of coordinate singularity: axisymmetric; polar coordinates;}

!latex \begin{enumerate}
!latex \item Consider a general, magnetic vector potential given in Cartesian coordinates,
!latex       \be {\bf A} = A_x \nabla x + A_y \nabla y +  A_z \nabla z + \nabla g \label{eq:CartesianVectorPotential}
!latex       \ee
!latex       where $A_x$, $A_y$, $A_z$, and the as-yet-arbitrary gauge function, $g$, are regular at $(x,y)=(0,0)$,
!latex       i.e. they can be expanded as a Taylor series, e.g.
!latex       \be A_x = \sum_{i,j} \alpha_{i,j} x^i y^j, \qquad
!latex           A_y = \sum_{i,j}  \beta_{i,j} x^i y^j, \qquad
!latex           A_z = \sum_{i,j} \gamma_{i,j} x^i y^j, \qquad
!latex            g  = \sum_{i,j} \delta_{i,j} x^i y^j, \label{eq:Taylorexpansion}
!latex       \ee
!latex       for small $x$ and small $y$.
!latex \item Note that we have restricted attention to the ``axisymmetric'' case, as there is no dependence on $z$.
!latex \item The `polar' coordinate transformation,
!latex       \be x &=& r \cos \t, \nonumber \\ y &=& r \sin \t, \\ z &=& \z, \nonumber
!latex       \ee
!latex       induces the vector transformation
!latex       \be \begin{array}{ccccccccccccccccccccccccccccccc} \nabla x & = & \cos \t \; \nabla r & - & r \sin \t \; \nabla \t & &, \\
!latex                                                          \nabla y & = & \sin \t \; \nabla r & + & r \cos \t \; \nabla \t & &, \\
!latex                                                          \nabla z & = &                     &   &                        & \nabla \z&.
!latex       \end{array} \ee
!latex \item By repeated applications of the double-angle formula, the expressions for $A_x$, $A_y$ and $g$ can be cast as functions of $(r,\t)$,
!latex       \be A_x & = & \sum_m r^m [ a_{m,0} + a_{m,1} \; r^2 + a_{m,2} \; r^4 + \dots ] \sin(m\t), \\
!latex           A_y & = & \sum_m r^m [ b_{m,0} + b_{m,1} \; r^2 + b_{m,2} \; r^4 + \dots ] \cos(m\t), \\
!latex           A_z & = & \sum_m r^m [ c_{m,0} + c_{m,1} \; r^2 + c_{m,2} \; r^4 + \dots ] \cos(m\t), \label{eq:regularAz} \\
!latex            g  & = & \sum_m r^m [ g_{m,0} + g_{m,1} \; r^2 + g_{m,2} \; r^4 + \dots ] \sin(m\t),
!latex       \ee
!latex       where attention is restricted to stellarator symmetric geometry, but similar expressions hold for the non-stellarator symmetric terms.
!latex \item Collecting these expressions, the vector potential can be expressed
!latex       \be {\bf A} = A_r \nabla r + A_\t \nabla \t + A_\z \nabla \z + \partial_r g \; \nabla r + \partial_\t g \; \nabla \t,
!latex       \ee
!latex       where
!latex       \be \begin{array}{ccccccclcccccccccccccccccccccccccccccccccccccccccccccccccccc}
!latex       A_r & = & r^0 & [ & ( &           &   & b_{0,0}   & + &   g_{1,0} & ) & + & (\dots) r^2 + (\dots)r^4 + \dots & ] & \sin   \t \\
!latex           & + & r^1 & [ & ( & a_{1,0}/2 & + & b_{1,0}/2 & + & 2 g_{2,0} & ) & + & (\dots) r^2 + (\dots)r^4 + \dots & ] & \sin 2 \t \\
!latex           & + & r^2 & [ & ( & a_{2,0}/2 & + & b_{2,0}/2 & + & 3 g_{3,0} & ) & + & (\dots) r^2 + (\dots)r^4 + \dots & ] & \sin 3 \t \\
!latex           & + & \dots
!latex       \end{array} \ee
!latex       (Note: Mathematica was used to perform the algebraic manipulations,
!latex       and the relevant notebook is included as part of the \SPEC {\small CVS} repository.)
!latex \item There is precisely enough gauge freedom so that we may choose $A_r=0$.
!latex       For example, the choice
!latex       \be \begin{array}{cccccclccccccccccccccccccccccccccccccccc}
!latex           g_{1,0} & = & - &   &           &   & b_{0,0}    &   &   &  &, \\
!latex           g_{2,0} & = & - & ( & a_{1,0}/2 & + & b_{1,0}/2  & ) & / & 2&, \\
!latex           g_{3,0} & = & - & ( & a_{2,0}/2 & + & b_{2,0}/2  & ) & / & 3&, \\
!latex           \dots   & = &   & \dots
!latex       \end{array} \ee
!latex       eliminates the lowest order $r$ dependence in each harmonic.
!latex \item By working through the algebra (again, using Mathematica) the expressions for $A_\t$ and $A_\z$ become
!latex       \be A_\t &=& r^2     f_0(\rho) + r^3 f_1(\rho) \cos(\t) + r^4 f_2(\rho) \cos(2\t) + r^5 f_3(\rho) \cos(3\t) + \dots \label{eq:nearoriginAt} \\
!latex           A_\z &=& \;\;\;\,g_0(\rho) + r^1 g_1(\rho) \cos(\t) + r^2 g_2(\rho) \cos(2\t) + r^3 g_3(\rho) \cos(3\t) + \dots \label{eq:nearoriginAz}
!latex       \ee
!latex       where $\rho\equiv r^2$ and the $f_m(\rho)$ and $g_m(\rho)$ are abitrary polynomials in $\rho$.
!latex       [The expression for $A_\z$ is unchanged from \Eqn{regularAz}.]
!latex \end{enumerate}

!latex \subsubsection{somewhat generally, . . . }

!latex \begin{enumerate}
!latex \item For stellarator-symmetric configurations,
!latex       \be {\bf A} = \sum_{m,n} A_{\t,m,n} \cos(m\t-n\z) \nabla \theta + \sum_{m,n} A_{\z,m,n} \cos(m\t-n\z) \nabla \zeta,
!latex       \ee
!latex       where now the dependence on $\z$ is included, and the angles are arbitrary.
!latex \item The near-origin behaviour of $A_\t$ and $A_\z$ given in \Eqn{nearoriginAt} and \Eqn{nearoriginAz} are flippantly generalized to
!latex       \be
!latex           A_{\t,m,n} & = & r^{m+2}          f_{m,n}(\rho), \label{eq:Atmn} \\
!latex           A_{\z,m,n} & = & r^{m  } \;\;\;\; g_{m,n}(\rho), \label{eq:Azmn}
!latex       \ee
!latex       where the $f_{m,n}(\rho)$ and $g_{m,n}(\rho)$ are arbitrary polynomials in $\rho$.
!latex \item Additional gauge freedom can be exploited: including an additional gauge term $\nabla h$ where $h$ only depends on $\z$, e.g.
!latex       \be h(\z) = h_{0,0} \, \z + \sum h_{0,n} \sin( - n\z),\ee
!latex       does not change the magnetic field and does not change any of the above discussion.
!latex \item The representation for the $A_{\t,m,n}$ does not change, but we must clarify that \Eqn{Azmn} holds for only the $m\ne0$ harmonics:
!latex       \be A_{\z,m,n} & = & r^{m} \;\;\;\;g_{m,n}(\rho), \;\;\;\mbox{\rm for $m \ne 0$}.
!latex       \ee
!latex \item For the $m=0$, $n\ne0$ harmonics of $A_\z$, including the additional gauge gives $A_{\z,0,n} = g_{0,n}(\rho) + n \, h_{0,n}$.
!latex       Recall that $g_{0,n}(\rho) = g_{0,n,0} + g_{0,n,1} \rho + g_{0,n,2} \rho^2 + \dots$, and we can choose $h_{0,n} = - g_{0,n,0} / n$ to obtain
!latex       \be A_{\z,m,n} & = & r^{m  } \;\;\;g_{m,n}(\rho), \;\;\;\mbox{\rm for $m=0$, $n\ne0$, with $g_{m,n}(0)=0$}.
!latex       \ee
!latex \item For the $m=0$, $n=0$  harmonic of $A_\z$, we have $A_{\z,0,0} = g_{0,0}(\rho) + h_{0,0}$.
!latex       Similarly, choose $h_{0,0} = - g_{0,n,0}$ to obtain
!latex       \be A_{\z,m,m} & = & r^{m  } \;\;\;g_{m,n}(\rho), \;\;\;\mbox{\rm for $m=0$, $n=0$, with $g_{m,n}(0)=0$}.
!latex       \ee
!latex \item To simplify the algorithmic implementation of these conditions,
!latex       we shall introduce a `regularization' factor, $\rho^{m/2} = r^m$.
!latex \item Note that the representation for $A_{\t,m,n}$ given in \Eqn{Atmn},
!latex       with an arbitrary polynomial $f_{m,n}(\rho) = f_{m,n,0} + f_{m,n,1}\rho + f_{m,n,2}\rho^2 + \dots$,
!latex       is equivalent to $A_{\t,m,n} = \rho^{m/2} \alpha_{m,n}(\rho)$ where $\alpha_{m,n}(\rho)$ is an arbitrary polynomial
!latex       with the constraint $\alpha_{m,n}(0)=0$.
!latex \item We can write the vector potential as
!latex       \be A_{\t,m,n} & = & \rho^{m/2} \alpha_{m,n}(\rho), \;\;\; \mbox{\rm with $\alpha_{m,n}(0)=0$ for all $(m,n)$},\\
!latex           A_{\z,m,n} & = & \rho^{m/2}  \beta_{m,n}(\rho), \;\;\; \mbox{\rm with $ \beta_{m,n}(0)=0$ for $m=0$.}
!latex       \ee
!latex \end{enumerate}

!latex \subsubsection{non-stellarator symmetric terms}

!latex \begin{enumerate}
!latex \item Just guessing, for the non-stellarator-symmetric configurations,
!latex       \be A_{\t,m,n} & = & \rho^{m/2} \alpha_{m,n}(\rho), \;\;\; \mbox{\rm with $\alpha_{m,n}(0)=0$ for all $(m,n)$},\\
!latex           A_{\z,m,n} & = & \rho^{m/2}  \beta_{m,n}(\rho), \;\;\; \mbox{\rm with $ \beta_{m,n}(0)=0$ for $m=0$.}
!latex       \ee
!latex \end{enumerate}

!l tex \subsubsection{constraints on Fourier harmonics}

!l tex \begin{enumerate}
!l tex \item The boundary conditions must accomodate the coordinate singularity by including the regularization factors:
!l tex \be A_{\t,e,1}(-1) = 0 \; & \; \mbox{\rm and} \; & \;\; A_{\t,e,1}(+1) = \Delta\psi_{t,1} \\
!l tex     A_{\t,e,i}(-1) = 0 \; & \; \mbox{\rm and} \; & \;\; A_{\t,e,i}(+1) = +m_i f_{o,i}     \\
!l tex     A_{\t,o,i}(-1) = 0 \; & \; \mbox{\rm and} \; & \;\; A_{\t,o,i}(+1) = -m_i f_{e,i}     \\
!l tex     A_{\z,e,1}(-1) = ? \; & \; \mbox{\rm and} \; & \;\; A_{\z,e,1}(+1) = \Delta\psi_{p,1} \\
!l tex     A_{\z,e,i}(-1) = 0 \; & \; \mbox{\rm and} \; & \;\; A_{\z,e,i}(+1) = -n_i f_{o,i} \mbox{ for } m_i =   0,             n_i \ne 0 \\
!l tex     A_{\z,e,i}(-1) = ? \; & \; \mbox{\rm and} \; & \;\; A_{\z,e,i}(+1) = -n_i f_{o,i} \mbox{ for } m_i \ne 0,                       \\
!l tex     A_{\z,o,i}(-1) = 0 \; & \; \mbox{\rm and} \; & \;\; A_{\z,o,i}(+1) = +n_i f_{e,i} \mbox{ for } m_i =   0,             n_i \ne 0 \\
!l tex     A_{\z,o,i}(-1) = ? \; & \; \mbox{\rm and} \; & \;\; A_{\z,o,i}(+1) = +n_i f_{e,i} \mbox{ for } m_i \ne 0,
!l tex \ee
!l tex       The `$?$' symbol indicates that there is no boundary condition, and $\Delta\psi_{p,1}$ is the magnetic flux linking the torus:
!l tex       it must be specified, but it is otherwise irrelevant.
!l tex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine manual
    use mod_kinds, only: wp => dp
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    use constants, only: zero

    use numerical, only:

    use fileunits, only: ounit

    use inputlist, only: Wmanual

    use cputiming, only: Tmanual

    use allglobal, only: myid, cpus

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef OPENMP
    USE OMP_LIB
#endif
    use mpi
    implicit none
    integer :: ierr, astat, ios, nthreads, ithread
    real(wp) :: cput, cpui, cpuo = 0 ! cpu time; cpu initial; cpu old; 31 Jan 13;

    cpui = MPI_WTIME()
    cpuo = cpui
#ifdef OPENMP
    nthreads = omp_get_max_threads()
#else
    nthreads = 1
#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    ! this "routine" is purely for documentation; 08 Feb 16;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

9999 continue
    cput = MPI_WTIME()
    Tmanual = Tmanual + (cput - cpuo)
    return

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine manual

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
