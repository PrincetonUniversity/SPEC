!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Supporting documentation.

!latex \end{enumerate}

!l!t!x \subsection{A strongly non-axisymmetric example \SPEC calculation.}

!l!t!x A simple, verfication calculation is performed to confirm that \SPEC is correctly treating
!l!t!x (i) the coordinate singularity and (ii) the spectral-condensation constraints in (iii) stellarator geometry.

!l!t!x \begin{enumerate}

!l!t!x \item We consider the strongly-shaped boundary:
!l!t!x \be R &=& R_0            + r_o \cos\t + \delta \cos(\t-\z) + \delta \cos(2\t), \\
!l!t!x     Z &=& \;\;\;\;\;\;   - r_o \sin\t + \delta \sin(\t-\z),
!l!t!x \ee
!l!t!x with $R_0=1$, $r_o=0.3$ and $\delta=0.06$; which is shown in \Fig{Poincare}.

!l!t!x \item The pressure, which for reasons that should become clear below, is set to zero, i.e. $p=0$.

!l!t!x \item \SPEC is executed in `fixed-transform' mode, for which the parallel current (i.e. $\mu$, where $\nabla\times{\bf B}=\mu{\bf B}$)
!l!t!x       and enclosed poloidal flux, $\Delta\psi_p$, (only for annulur volumes, i.e. not in the innermost toroidal volume) are adjusted
!l!t!x       to match the prescribed values of the rotational-transform on the ideal interfaces.

!l!t!x \item Given these inputs, the geometry of the internal interfaces, if any, are adjusted in order to satisfy force balance.

!l!t!x \item Two \SPEC calculations are performed, the first with a single relaxed volume, i.e. $N_V=1$, and the second with $N_V=2$ relaxed volumes.

!l!t!x \item The first calculation, with $N_V=1$ and thus with no internal interfaces,
!l!t!x       requires only the total enclosed toroidal flux, which is normalized to $\psi_{t,edge}=1$, and $\iotabar_{edge}$ to be supplied.

!l!t!x \item Consider the case where $\iotabar_{edge}$ has been iteratively adjusted so that $\mu$ for the $N_V=1$ calculation is equal to a presribed value, e.g. $\bar\mu=6$.
!l!t!x       A \Poincare plot of this field is shown in \Fig{Poincare}.
!l!t!x       This solution satisfies $\nabla\times{\bf B}=\bar\mu{\bf B}$ over the entire domain.
!l!t!x \item (Note: the value of $\mu=\bar\mu$ could have been fixed, and then \SPEC could be executed in `fixed-parallel-current' mode,
!l!t!x       but for some rather-trivial numerical reasons it was easier to iterate on $\iotabar$ in order to get the required $\mu$,
!l!t!x       and this method is more consistent with the approach used for the $N_V=2$ calculation that follows.)

!l!t!x \item From following the magnetic fieldlines, the approximate rotational-transform in the domain can be determined,
!l!t!x       and we observe that the $\iotabar = (3+\gamma 4)/(1+\gamma 1) = 3.382\dots$, where $\gamma$ is the golden mean, irrational flux surface is present.

!l!t!x \item The second calculation, with $N_V=2$ and thus with one internal interface, requires $\psi_{t,1}$, the toroidal flux enclosed in the innermost toroidal volume
!l!t!x       to be specified, in addition to $\psi_{t,edge}$.
!l!t!x       The transform immediately inside and outside the ideal interface, $\iotabar_{1}$, must also be supplied, and we choose $\iotabar_{1} = 3.382\dots$

!l!t!x \item Consider now the case where, for the $N_V=2$ equilibrium, that both $\iotabar_{edge}$ and $\psi_{t,1}$ have been iteratively adjusted so that 
!l!t!x       $\mu_i=\bar\mu$ in {\em both} volumes.
!l!t!x       This solution satisfies $\nabla\times{\bf B}=\bar\mu{\bf B}$ in the first volume and $\nabla\times{\bf B}=\bar\mu{\bf B}$ in the second volume.

!l!t!x \item The geometry of the internal interface is shown in \Fig{Poincare}.

!l!t!x \item The calculation is encouraging for the following reasons: 

!l!t!x \item With no pressure discontinuity at the internal interface, and with the values of the enclosed toroidal fluxes and transforms chosen accordingly,
!l!t!x       (and assuming that the $\iotabar=3.382\dots$ irrational surface exists in the $N_V=1$ calculation,)
!l!t!x       the $N_V=2$ calculation should be identical (up to Fourier resolution) with the $N_V=1$ calculation;
!l!t!x       and so the geometry of the ideal interface in the $N_V=2$ calculation should coincide with the irrational flux surface with $\iotabar = 3.382\dots$
!l!t!x       in the $N_V=1$ calculation, which it does.

!l!t!x \item The are two important numerical differences in the calculations:
!l!t!x       (i) The calculation of the Beltrami field in the toroidal volume requires careful treatment of the coordinate singularity,
!l!t!x           whereas the annular regions do not have a coordinate singularity.
!l!t!x       (ii) The $N_V=2$ calculation requires the use of spectral-condensation constraints to adequately constrain the Fourier harmonics that describe the geometry of the interior interface.
!l!t!x       So, different subroutines (that enforce regularity at the origin or not) and different constraints (the inclusion of angle constraints or not) 
!l!t!x       are used in the different calculations, but the results appear the same.

!l!t!x \item This is just a `reality-check', as I have not yet quantified the error etc.
!l!t!x       There are still several aspects of the algorithm that I wish to revise, but I am cautiously hopeful that \SPEC should be ready to perform W7-X, LHD calculations etc.

!l!t!x \insertdblfigure{Nv=01a.opt.poin.c.triple.ps}{A \Poincare plot of the $N_V=1$ calculation (grey dots), and the interior interface of the $N_V=2$ calculation (black line), on three different toroidal planes.}{Poincare}

!l!t!x \item An indication of the robustness of the \SPEC convergence properties is provided by \Fig{Initialization}, and \Fig{iterations},
!l!t!x       in which the initial geometry of the interior interface
!l!t!x       (obtained by a simple extrapolation of the Fourier harmonics of the supplied boundary) in the $N_V=2$ calculation is compared to the geometry of the final solution.
!l!t!x       \Fig{iterations} shows the `force-error', as a function of iteration.
!l!t!x       With the Fourier resolution \verb+Mpol=6+ and \verb+Ntor=6+, the $N_V=2$ calculation required about $2$ minutes.
!l!t!x       (Note this is after the `derivative' matrix, $\nabla_{\bf x}{\bf F}$, was calculated, which for this example required about $10$ minutes.
!l!t!x        The construction of the derivative matrix is an `overhead' calculation, by which I mean that the derivative matrix need {\em not} be recomputed
!l!t!x        when multiple equilibrium calculations are performed with the same Fourier resolution
!l!t!x        and, say, slowing varying pressure, rotational-transform, or boundary shaping $\dots$)

!l!t!x \item Adding additional volumes and pressure is almost trivial: it was the coordinate singularity, the spectral constraints, and the strongly shaped boundary that this calculation was intended to investigate.
!l!t!x       The routines for computing the Beltrami fields in the annulur regions has been well tested;
!l!t!x       and \SPEC is parallelized over the volumes, so adding additional volumes can actually make the calculation faster provided additional cpus are available.
!l!t!x       Adding pressure is almost as simple as adding a constant to the calculation of $B^2$ at the interfaces.

!l!t!x \insertdblfigure{Nv=02b.poin.c.triple.ps}{Initial guess for the geometry of the interior interface of the $N_V=2$ calculation (red), as obtained by extrapolating the Fourier harmonics of the given boundary; and the final geometry that is consistent with force-balance (black).}{Initialization}

!l!t!x \insertsglfigure{Nv=02c.fits.ps}{Reduction of force-imbalance, i.e. $|{\bf F}|$, against iteration.}{iterations}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!l!t!x \subsection{rotational-transform}

!latex \begin{enumerate}

!l!t!x \item This brief note investigates the relationship between a singular current density and the rotational-transform, assuming an integrable magnetic field.

!l!t!x \item An integrable magnetic field may be written
!l!t!x       \be {\bf B} = \nabla \psi \times \nabla \bar \theta + \nabla \zeta \times \nabla \chi(\psi). \label{eq:canonical}
!l!t!x       \ee
!l!t!x \item This satisfies both $\nabla\cdot{\bf B}=0$, because ${\bf B}=\nabla\times(\psi\nabla\t-\chi\nabla\z)$, and ${\bf B}\cdot\nabla\psi=0$ by construction.
!l!t!x \item By taking certain surface integrals, the functions $\psi$ and $\chi$ may be identified as the toroidal and poloidal flux enclosed by a flux surface.
!l!t!x \item The flux surfaces may be labelled with an arbitrary parameter, e.g. $s\equiv$ the volume enclosed by a flux surface, which we assume is a smooth function of postion;
!l!t!x       and the toroidal and poloidal fluxes are surface functions, i.e. $\psi\equiv\psi(s)$ and $\chi\equiv\chi(s)$.
!l!t!x \item The rotational-transform is an angle-independent quantity.
!l!t!x       It may be defined as the ratio of the differential fluxes, 
!l!t!x       \be \iotabar \equiv \frac{d\chi}{d\psi} = \frac{\chi^\prime}{\psi^\prime}.
!l!t!x       \ee
!l!t!x \item Near a surface that may have a discontinuity in $\iotabar$, one-sided limits will need to be taken.
!l!t!x \item The $\bar\t$ and $\z$ form straight-fieldline angles. We are free to choose the geometrical cylindrical angle as the toroidal angle, i.e. $\z\equiv\phi$.
!l!t!x       The straight-fieldline poloidal angle, $\bar\t$, on every surface with {\em irrational} rotational-transform is then {\em almost} unique.
!l!t!x \item The straight-fieldline angle is only almost unique, because a rigid shift of the angle, e.g. $\bar\t\rightarrow\bar\t+c$, for some constant $c$,
!l!t!x       also gives a straight-fieldline angle.
!latex \item On rational surfaces, where $\iotabar=n/m$, the straight fieldline poloidal angle is {\em not} unique.
!latex \item To prove this, let $\tilde\t =  \t + \delta_{m,n} \sin(m\t-n\z)$ where $\t$ satisfies the straight fieldline angle condition, 
!latex       ${\bf B}\cdot\nabla\t = \iotabar \, {\bf B}\cdot\nabla\z$.
!latex \item Then,
!latex       ${\bf B}\cdot\nabla\tilde\t = {\bf B}\cdot\nabla\t + ( m \, {\bf B} \cdot \nabla  \t - n \, {\bf B} \cdot \nabla \z ) \,
!latex       \delta_{m,n} \cos( m  \t - n \z )$, 
!latex       which reduces to ${\bf B}\cdot\nabla\tilde\t = \iotabar \, {\bf B}\cdot\nabla\z$ if $\iotabar=n/m$.
!latex \item If required, it is possible to define a `preferred' straight fieldline angle on the rational surfaces as
!latex       that angle which is smoothly connected to the straight fieldline angle on the infinitesimally-nearby irrational surfaces.

!l!t!x \end{enumerate}

!l!t!x \subsection{suitable coordinates}

!l!t!x \begin{enumerate}

!l!t!x \item It is convenient {\em not} to use $\bar\t$ as a coordinate, as will consider the possibility that $\bar\t$ may not be a smooth, continuous function of position.
!l!t!x       Instead, we use $(\s,\t,\z)$ as coordinates, where $\t$ is an arbitrary poloidal angle, e.g. the geometric polar angle, or the equal-arc-length angle.
!l!t!x \item Assuming that $s$, $\t$ and $\z$ are continuous and smooth functions of position,
!l!t!x       the metric elements, $g_{\s\s}$, $g_{\s\t}$, $g_{\s\z}$, $g_{\t\t}$, $g_{\t\z}$ and $g_{\z\z}$, are each continuous.
!l!t!x \item The transformation from $\t$ to $\bar\t$ is generally
!l!t!x       \be \bar\t=\t+\lambda(\s,\t,\z), \label{eq:lambda}
!l!t!x       \ee
!l!t!x       where $\lambda=\sum \lambda_{m,n}(s)\exp(im\t-in\z)$ is a unique, single-valued function of position, with the following exceptions:
!l!t!x \item As noted above, a rigid shift of one straight fieldline angle creates another, i.e. $\lambda_{0,0}$ is not constrained.
!l!t!x       This amounts to a somewhat-trivial re-definition of the where on the flux surface that $\bar\t=0$, and this degree-of-freedom will hereafter be ignored.
!l!t!x       To be precise, we can assume `stellarator-symmetry', and restrict attention to transformations of the form $\lambda=\sum \lambda_{m,n}(s)\sin(m\t-n\z)$.
!l!t!x \item Also as noted above, on the rational surfaces the straight fieldline angle is not unique.
!l!t!x       For an angle transformation from one straight fieldline angle to another,
!l!t!x       this freedom means that $\lambda_{m,n}$ is not constrained for all $(m,n)$ that `resonate', i.e. for which $n/m=\iotabar$.
!l!t!x       An angle transformation from an arbitrary poloidal angle to a straight fieldline angle is similarly, not-fully constrained;
!l!t!x       it is, however, more difficult to identify the unconstrained harmonics (i.e. singular value decomposition methods are required).
!l!t!x \item We allow for the possibility that there is a discontinuity in the tangential field at a given surface, $\s=s_0$.
!l!t!x       We define the single-valued functions 
!l!t!x       $\lambda^{-}(\t,\z)\equiv\lambda(\s_0-\epsilon,\t,\z)$ and $\lambda^{+}(\t,\z)\equiv\lambda(\s_0+\epsilon,\t,\z)$,
!l!t!x       in the limit that $\epsilon\rightarrow0$.
!l!t!x \item Assuming the angle transformation given in \Eqn{lambda}, the magnetic field, \Eqn{canonical}, may be written in $(\s,\t,\z)$ coordinates as
!l!t!x       \be {\bf B} = \left(\chi' - \psi' \lambda_\z \right) \frac{{\bf e}_\t}{\sqrt g} + \left( \psi' + \psi' \lambda_\t \right) \frac{{\bf e}_\z}{\sqrt g}.
!l!t!x       \ee
!l!t!x \item The covariant form is ${\bf B} = B_\s \nabla \s + B_\t \nabla \t + B_\z \nabla \z$,
!l!t!x       where
!l!t!x       \be
!l!t!x       B_\mu = \left(\chi' - \psi' \lambda_\z \right) \frac{g_{\t\mu}}{\sqrt g} + \left( \psi' + \psi' \lambda_\t \right) \frac{g_{\z\mu}}{\sqrt g}, 
!l!t!x           \label{eq:Bt}
!l!t!x       \ee
!l!t!x       for $\mu$ in $\{ \s,\t,\z \}$.

!l!t!x \end{enumerate}

!l!t!x \subsection{current density singularity}

!l!t!x \begin{enumerate}

!l!t!x \item Imagine that there is current density that satisfies
!l!t!x       (i) $\nabla \cdot {\bf j}=0$, and (ii) ${\bf j} \cdot \nabla s=0$, and (iii) is localized to a flux surface $s=s_0$.
!l!t!x \item Using an appropriate `gauge', such a current density may be represented
!l!t!x       \be {\bf j} = \nabla \times ( \kappa \nabla s )% = \kappa_\z \nabla \z \times \nabla \s - \kappa_\t \nabla \s \times \nabla \t,
!l!t!x       \ee
!l!t!x       where generally $\kappa(\s,\t\,z) \equiv I \, \delta(s-s_0) \, \t + G \, \delta(s-s_0) \, \z + \tilde \kappa(\s,\t,\z)$,
!l!t!x       where $\tilde \kappa$ is single-valued.
!l!t!x \item The toroidal current passing through an infinitesimal surface about $s=s_0$ is 
!l!t!x       \be \delta j_\z 
!l!t!x           & \equiv &   \int_{\cal S} {\bf j} \cdot d{\bf s} 
!l!t!x               =        \int_{s_0-\epsilon}^{s_0+\epsilon} \int_{0}^{2\pi} {\bf j} \cdot {\bf e}_\s \times {\bf e}_\t \;\; d\t \; d\s \label{eq:surfacecurrent}
!l!t!x               =        \int_{s_0-\epsilon}^{s_0+\epsilon} \int_{0}^{2\pi} {\bf j} \cdot \nabla \z \sqrt g \;\; d\t \; d\s \\
!l!t!x           &   =    & - \int_{s_0-\epsilon}^{s_0+\epsilon} \int_{0}^{2\pi} \left[ I \, \delta(s-s_0) + \tilde \kappa_\t \right]       \;\; d\t \; d\s \\
!l!t!x           &   =    & - 2\pi I.\label{eq:surfacecurrentc}
!l!t!x       \ee
!l!t!x \item Now, assume that ${\bf j} = \nabla \times {\bf B}$, where ${\bf B}$ is single valued:
!l!t!x       \be \delta j_\z 
!l!t!x           & \equiv &   \int_{\cal S} {\bf j} \cdot d{\bf s} \nonumber \\
!l!t!x           &   =    &   \int_{\partial \cal S} {\bf B}\cdot d{\bf l} \nonumber \\
!l!t!x           &   =    &   \int_{0}^{2\pi} [[B_\t]] \, d\t. \label{eq:surfacecurrentb}
!l!t!x       \ee
!l!t!x \item To obtain the result $\delta j_\z = - 2 \pi I$, it is neccessary to recognize that the magnetic field given in \Eqn{canonical},
!l!t!x       and therefore also the expressions given in \Eqn{Bt}, are not yet completely constrained.
!l!t!x       Enforcing the condition that ${\bf j}=\nabla\times{\bf B}$ is tangential to the flux surfaces, i.e. ${\bf j}\cdot\nabla s=0$,
!l!t!x       requires that $\partial_\t B_\z - \partial_\z B_\t=0$.
!l!t!x       This is satisfied by
!l!t!x       $B_\t \equiv \partial_\t f$ and $B_\z \equiv \partial_\z f$,
!l!t!x       for $f\equiv I \, \delta(s-s_0) \, \t + G \, \delta(s-s_0) \, \z + \tilde f(\s,\t,\z)$ and where $\tilde f$ is single-valued.
!l!t!x \item Then \Eqn{surfacecurrentb} yields  $\delta j_\z = 2 \pi I$, which is identical to \Eqn{surfacecurrentc}, except for perhaps a sign change!
!l!t!x \item The magnetic field given in \Eqn{canonical} simultaneously satisfies $\nabla\cdot{\bf B}=0$ and ${\bf B}\cdot\nabla\psi=0$ by construction.
!l!t!x       The condition that ${\bf j}\cdot\nabla\psi=0$ for ${\bf j}=\nabla\times{\bf B}$ reduces to
!l!t!x       \be I + \tilde f_\t = \left(\chi' - \psi' \lambda_\z \right) \frac{g_{\t\t}}{\sqrt g} + \left( \psi' + \psi' \lambda_\t \right) \frac{g_{\t\z}}{\sqrt g}, 
!l!t!x           \label{eq:jt} \\
!l!t!x           G + \tilde f_\z = \left(\chi' - \psi' \lambda_\z \right) \frac{g_{\t\z}}{\sqrt g} + \left( \psi' + \psi' \lambda_\t \right) \frac{g_{\z\z}}{\sqrt g},
!l!t!x           \label{eq:jz}
!l!t!x       \ee
!l!t!x       which would seem to place a constraint on the geometry.
!l!t!x \item Note that the constraint of force-balance, ${\bf j}\times{\bf B}=\nabla p$, has not yet been enforced.
!l!t!x \end{enumerate}

!l!t!x \subsection{continuous rotational-transform constraint}\label{sec:singular current and transform}

!l!t!x \begin{enumerate}

!l!t!x \item Using \Eqn{Bt}, \Eqn{surfacecurrentb} becomes
!l!t!x       \be 2 \pi I = \int_{0}^{2\pi} \left( [[\chi' - \psi' \lambda_\z ]] \frac{g_{\t\t}}{\sqrt g}
!l!t!x                                          + [[\psi' + \psi' \lambda_\t ]] \frac{g_{\t\z}}{\sqrt g} \right) \; d\t.
!l!t!x       \ee
!l!t!x \item If the flux surface $s=s_0$ has irrational transform, there is no ambiguity in the angle transformation.
!l!t!x \item In order to avoid the ambiguity in the angle transformation at rational sufaces,
!l!t!x let $\epsilon$ in \Eqn{surfacecurrent} be restricted so that $s=s_0\pm\epsilon$ identifies irrational surfaces.

!l!t!x \item Let us assume that both $\psi'$ and $\chi'$ are continuous.
!l!t!x       In this case, the rotational-transform is continuous.
!l!t!x \item Then, 
!l!t!x       \be 2 \pi I = \psi' \int_{0}^{2\pi} \!\! \left( [[\lambda_\t]] \frac{g_{\t\z}}{\sqrt g} - [[\lambda_\z]] \frac{g_{\t\t}}{\sqrt g}\right) \, d\t.
!l!t!x       \label{eq:angleintegral}
!l!t!x       \ee
!l!t!x \item If the discontinuity in the straight-fieldline angle resonates with the geometry, i.e. if the integral in \Eqn{angleintegral} is non-zero,
!l!t!x       then the conditions (i) $\delta j_\z\ne0$, and (ii) both $\psi'$ and $\chi'$ are continuous, can be satisfied simultaneously.
!l!t!x       This would appear to be possible in arbitrary geometry,
!l!t!x       and a singular current with a nonvanishing integral {\em can} exist {\em without} there being a discontinuity in the rotational-transform.
!l!t!x \item However, perhaps it is more complicated.
!l!t!x       It has not yet been enforced that ${\bf j} \cdot \nabla s =0$, or that ${\bf j}\times{\bf B}=\nabla p$.
!l!t!x \item Consider the stellarator symmetric case, in which $\lambda$ is an odd function, and so $\lambda_\t$ and $\lambda_\z$ are both even.
!l!t!x       Assuming that $g_{\t\t}/\sqrt g$ and $g_{\t\z}/\sqrt g$ are both even, the integrand in \Eqn{angleintegral},
!l!t!x       $f \equiv ( [[\lambda_\t]] g_{\t\z} - [[\lambda_\z]] g_{\t\t})/\sqrt g$, is even, and so it may be written $f \equiv \sum f_{m,n}\cos(m\t-n\z)$.
!l!t!x       The integral reduces to 
!l!t!x       \be 2 \pi I = \psi' \sum_{n} f_{0,n} \cos(-n\z)2\pi,
!l!t!x       \ee
!l!t!x       which, given that $2\pi I$ does not depend on $\z$, must reduce to 
!l!t!x       \be 2 \pi I = \psi' f_{0,0} 2\pi.
!l!t!x       \ee
!l!t!x \item Given the constraints in \Eqn{jt} and \Eqn{jz}, and perhaps also the constraints implied by force-balance, is there any freedom to satisfy this?

!l!t!x \end{enumerate}

!l!t!x \subsection{axisymmetric}\label{sec:axisymmetriccase}

!l!t!x \begin{enumerate}

!l!t!x \item It is interesting to consider the axisymmetric, cylindrical case, e.g. $x=s\cos\t$, $y=s\sin\t$ and $z=\z$.
!l!t!x       The Jacobian and relevant metric elements are $\sqrt g = s$, $g_{\t\z}=0$ and $g_{\t\t}=s^2$; and \Eqn{angleintegral} becomes
!l!t!x       \be \delta j_\z = - s_0 \, \psi' \int_{0}^{2\pi} \! \! [[\lambda_\z]] \, d\t,
!l!t!x       \ee
!l!t!x       The function $[[\lambda]]$ must be single valued, and so $[[\lambda_\z]]$ cannot have a non-zero, $m=0$ harmonic;
!l!t!x       and so $\delta j_\z$ must be zero if $\psi'$ and $\chi'$ are continuous.
!l!t!x       A singular current with a nonvanishing integral {\em cannot} exist in the symmetric case
!l!t!x       without there being a discontinuity in either $\psi'$, $\chi'$, or both.

!l!t!x \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!l!t!x \section{Rosenbluth analysis in Cartesian geometry}

!l!t!x \subsection{unperturbed equilibrium}

!l!t!x \begin{enumerate}

!l!t!x \item Position is generally
!l!t!x       \be {\bf x} = \t \, {\bf i} + \z \, {\bf j} + R(\s,\t,\z) \, {\bf k}.
!l!t!x       \ee
!l!t!x       The unperturbed geometry is $R=s$.

!l!t!x \item The calculation domain is $\s\in[-1,1]$.

!l!t!x \item The magnetic field is 
!l!t!x       \be {\bf B} = B^\t(s) \, {\bf i} + B^\z(s) \, {\bf j},
!l!t!x       \ee
!l!t!x       where we choose $B^\t \equiv \iotabar^\prime s B_0$ and, to satisfy ${\bf j}\times{\bf B}=0$, we must have $B^\z = \sqrt{B_0^2-(B^\t)^2}$ for some constant $B_0$.

!l!t!x \end{enumerate}

!l!t!x \subsection{linearized displacement}

!l!t!x \begin{enumerate}

!l!t!x \item Consider the boundary displacement $R(1,\t,\z) \equiv 1 + \epsilon \cos(m\t-n\z)$, where we choose $(m,n)=(1,0)$.

!l!t!x \item Solving $({\bf j}+\delta {\bf j})\times({\bf B}+\delta{\bf B})=0$ consistent with
!l!t!x \bi
!l!t!x \item[i.] the displacment of the flux surfaces and the perturbed magnetic field being ``ideally'' related via $\delta {\bf B}\equiv\nabla\times(\boldxi \times {\bf B})$;
!l!t!x \item[ii.]the boundary condition ${\bf B}\cdot\nabla s=0$ at $s=\pm 1$;
!l!t!x \ei
!l!t!x to ${\cal O}(\epsilon)$ gives the Newcomb equation:
!l!t!x \be s^2 \xi^{\prime\prime} + 2 s \, \xi^\prime - s^2 \xi = 0.
!l!t!x \ee

!l!t!x \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!l!t!x \newpage

!l!t!x    \section{Cartesian geometry - general}

!l!t!x \begin{enumerate}

!l!t!x \item Position is given by
!l!t!x       \be {\bf x} \equiv \t \; {\bf \hat i} + \z \; {\bf \hat j}+ R(\s,\t,\z) \; {\bf \hat k}.
!l!t!x       \ee
!l!t!x       The calculation domain is $\s\in[-1,+1]$, $\t\in[0,2\pi]$ and $\z\in[0,2\pi]$.

!l!t!x \item The position derivatives are
!l!t!x       \be  {\bf e}_\s  =                                    R_\s \; {\bf \hat k}, \;\;
!l!t!x            {\bf e}_\t  =  {\bf \hat i}                   +  R_\t \; {\bf \hat k}, \;\;
!l!t!x            {\bf e}_\z  =                   {\bf \hat j}  +  R_\z \; {\bf \hat k},
!l!t!x       \ee 
!l!t!x       and the Jacobian is $\sqrt g \equiv {\bf e}_\s \cdot {\bf e}_\t \times {\bf e}_\z = R_\s$.

!l!t!x \item The coordinate gradients are
!l!t!x          $\nabla \s = ( - R_\t {\bf \hat i} - R_\z {\bf \hat j} + {\bf \hat k} ) / R_\s$
!l!t!x          $\nabla \t =           {\bf \hat i}$, and
!l!t!x          $\nabla \z =           {\bf \hat j}$.

!l!t!x \item The metric elements are 
!l!t!x       \be g_{\s\s}  =  R_\s R_\s    , \;\;
!l!t!x           g_{\s\t}  =  R_\s R_\t    , \;\;
!l!t!x           g_{\s\z}  =  R_\s R_\z    , \;\;
!l!t!x           g_{\t\t}  =  R_\t R_\t + 1, \;\;
!l!t!x           g_{\t\z}  =  R_\t R_\z    , \;\;
!l!t!x           g_{\z\z}  =  R_\z R_\z + 1.  
!l!t!x       \ee

!l!t!x \item The coordinate function $R$ in the $v$-th volume is defined by a linear interpolation between the $v-1$ and $v$-th interfaces,
!l!t!x       \be R(\s,\t,\z) \equiv R_{v-1}(\t,\z) \frac{(1-s)}{2} + R_{v}(\t,\z) \frac{(1+s)}{2}.
!l!t!x       \ee

!l!t!x \item The general expression for the magnetic field is
!l!t!x       \be {\bf B} = \nabla \times \left[ A_\t(\s,\t,\z) \nabla \t + A_\z(\s,\t,\z) \nabla \z \right],
!l!t!x       \ee
!l!t!x       with the constraint $A_\t(-1,\t,\z)=0$ and$A_\z(-1,\t,\z)=0$.

!l!t!x \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!l!t!x \subsection{Cartesian geometry - symmetric solution} \label{sec:symmetric}

!l!t!x \begin{enumerate}

!l!t!x \item Choose $R_{v}(\t,\z) = R_{v,0}$.

!l!t!x \item Introduce the shorthand notation $\Delta_0 \equiv R_{v,0}-R_{v-1,0}$.

!l!t!x \item The metric elements reduce to
!l!t!x       \be
!l!t!x           g_{\s\s}  =  \Delta_0^2 / 4 , \;\;
!l!t!x           g_{\s\t}  =  0, \;\;
!l!t!x           g_{\s\z}  =  0, \;\;
!l!t!x           g_{\t\t}  =  1, \;\;
!l!t!x           g_{\t\z}  =  0, \;\;
!l!t!x           g_{\z\z}  =  1,
!l!t!x       \ee
!l!t!x       and the Jacobian is
!l!t!x       $\sqrt g = \Delta_0/2.
!l!t!x       $

!l!tex \item The inverse Jacobian is
!l!tex       \be (\sqrt g)^{-1} = \frac{2}{\Delta_0}.
!l!tex       \ee

!l!t!x \item Assuming poloidal and toroidal symmetry, so that \mbox{$A_\t(\s) \equiv \psi(s)$} and \mbox{$A_\z(\s) \equiv -\chi(s)$}, the magnetic field is
!l!t!x       \be {\bf B} = \frac{2\chi'}{\Delta_0} \, {\bf \hat i} + \frac{2\psi'}{\Delta_0} \, {\bf \hat j}.
!l!t!x       \ee

!l!t!x \item The enclosed toroidal and poloidal fluxes are
!l!t!x       \be \Delta \psi \equiv \int_{-1}^{+1} \!\!\! ds \int_{0}^{2\pi} \!\!\! d\t \;\;\; {\bf B} \cdot {\bf e}_\s \times {\bf e}_\t = 2 \pi \psi(1),\nonumber \\
!l!t!x           \Delta \chi \equiv \int_{-1}^{+1} \!\!\! ds \int_{0}^{2\pi} \!\!\! d\t \;\;\; {\bf B} \cdot {\bf e}_\z \times {\bf e}_\s = 2 \pi \chi(1).\nonumber
!l!t!x       \ee
!l!t!x \item The equation $\nabla \times {\bf B} = \mu {\bf B}$ becomes,
!l!t!x       \be  \begin{array}{ccccc} \psi^{\prime\prime} & = & - \bar \mu \, \chi^{\prime},\\
!l!t!x                                 \chi^{\prime\prime} & = & + \bar \mu \, \psi^{\prime},
!l!t!x            \end{array} \label{eq:unperturbed}
!l!t!x       \ee
!l!t!x       where $\bar \mu \equiv \mu \Delta_0 / 2$.
!l!t!x \item This has solution
!l!t!x       \be \begin{array}{ccccc} \psi(s) & = & a \cos(\bar\mu s) - b \sin(\bar\mu s) + \psi_0,\\
!l!t!x                                \chi(s) & = & a \sin(\bar\mu s) + b \cos(\bar\mu s) + \chi_0,
!l!t!x            \end{array}\label{eq:unperturbedsolution}
!l!t!x       \ee
!l!t!x       where $\psi_0$ and $\chi_0$ are constants.
!l!t!x \item The toroidal and poloidal flux constraints require that
!l!t!x       \be \begin{array}{rcccccccc}
!l!t!x           a \cos\bar\mu &+& b\sin\bar\mu &+& \psi_0 & = & 0, \\
!l!t!x           a \cos\bar\mu &-& b\sin\bar\mu &+& \psi_0 & = & \Delta\psi / 2\pi,\\
!l!t!x          -a \sin\bar\mu &+& b\cos\bar\mu &+& \chi_0 & = & 0,\\
!l!t!x           a \sin\bar\mu &+& b\cos\bar\mu &+& \chi_0 & = & \Delta\chi / 2\pi,
!l!t!x       \end{array} \ee
!l!t!x       which will be partially solved by writing $\psi_0\equiv - a \cos\bar\mu - b\sin\bar\mu$ and $\chi_0\equiv a \sin\bar\mu - b\cos\bar\mu$, leaving
!l!t!x       \be 
!l!t!x           -2 b \sin \bar\mu = \Delta\psi/2\pi, \label{eq:gaugefluxconstraints} \\
!l!t!x           +2 a \sin \bar\mu = \Delta\chi/2\pi.
!l!t!x       \ee
!l!t!x \item The rotational-transform is $\iotabar(s)\equiv \chi^\prime / \psi^\prime$, and the inner and outer interface transform constraints are
!l!t!x       \be \frac{a \cos(\bar\mu  ) + b \sin(\bar\mu  )}{+a \sin(\bar\mu  ) - b \cos(\bar\mu  )} & = & \iotabar_{v-1}, \label{eq:innertransform} \\
!l!t!x           \frac{a \cos(\bar\mu  ) - b \sin(\bar\mu  )}{-a \sin(\bar\mu  ) - b \cos(\bar\mu  )} & = & \iotabar_{v  }. \label{eq:outertransform}
!l!t!x       \ee
!l!t!x \item Given $\Delta \psi$, $\Delta \chi$ and $\mu$, the solution is determined. (This corresponds to using \verb+lconstraint=0+ in \SPEC.)
!l!t!x \item Alternatively, given $\Delta\psi$, $\iotabar_{v-1}$ and $\iotabar_v$,
!l!t!x       the equations \Eqn{gaugefluxconstraints},  \Eqn{innertransform} and \Eqn{outertransform} constrain the $3$ unknown constants, $a$, $b$ and $\bar\mu$.
!l!t!x       (This corresponds to using \verb+lconstraint=1+.)
!l!t!x       These non-linear equations require an iterative solution.

!l!t!x \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!l!t!x \subsection{Cartesian geometry - perturbed solution}

!l!t!x \begin{enumerate}

!l!t!x \item Choose $R_{v}(\t,\z) \equiv R_{v,0} + \epsilon \, R_{v,1} \cos(m\t-n\z)$,
!l!t!x       where $\epsilon$ is an ordering parameter and quantities ${\cal O}(\epsilon^2)$ will be ignored.
!l!t!x \item It is convenient to introduce rotating coordinates defined via
!l!t!x       \be \a = m \t - n \z, \;\; \bar \z = \z.
!l!t!x       \ee
!l!t!x \item Position is generally given by
!l!t!x       \be {\bf x} \equiv (\a + n\bar \z / m) \; {\bf \hat i} + \bar \z \; {\bf \hat j}+ R(\s,\a,\bar \z) \; {\bf \hat k}.
!l!t!x       \ee
!l!t!x \item The coordinate derivatives are
!l!t!x       \be  {\bf e}_\s  =                                             R_\s \; {\bf \hat k}, \;\;
!l!t!x            {\bf e}_\a  =              {\bf \hat i}                +  R_\a \; {\bf \hat k}, \;\;
!l!t!x            {\bf e}_\z  =  \frac{n}{m} {\bf \hat i} + {\bf \hat j} +  R_\z \; {\bf \hat k}.
!l!t!x       \ee
!l!t!x \item The metric elements are 
!l!t!x       \be g_{\s\s}  =  R_\s R_\s    , \;\;
!l!t!x           g_{\s\a}  =  R_\s R_\a    , \;\;
!l!t!x           g_{\s\z}  =  R_\s R_\z    , \;\;
!l!t!x           g_{\a\a}  =  R_\a R_\a + 1, \;\;
!l!t!x           g_{\a\z}  =  R_\a R_\z + n/m, \;\;
!l!t!x           g_{\z\z}  =  R_\z R_\z + 1 + n^2/m^2.  
!l!t!x       \ee
!l!t!x \item The Jacobian is
!l!t!x       \be \sqrt g = R_\s.
!l!t!x       \ee
!l!t!x \item Introduce the shorthand notation $\Delta_i \equiv R_{v,i}-R_{v-1,i}$ and $\bar\Delta_i \equiv R_{v,i}+R_{v-1,i}$.
!l!t!x \item The metric elements reduce to
!l!t!x       \be
!l!t!x           g_{\s\s}  \approx  \frac{1}{4}\Delta_0^2 + \epsilon \, \frac{1}{2}\cos\a \, \Delta_0 \, \Delta_1, \;\;
!l!t!x           g_{\s\a}  \approx  -\epsilon \, \frac{1}{4} \sin\a \, m \, \Delta_0 \; \left(\bar\Delta_1+\Delta_1 \; s\right),
!l!t!x       \ee
!l!t!x       with $g_{\s\z}  =        0$, 
!l!t!x            $g_{\a\a}  \approx  1$,
!l!t!x            $g_{\a\z}  = n/m$, and
!l!t!x            $g_{\z\z}  =        1+n^2/m^2$.
!l!t!x \item The Jacobian is
!l!t!x       \be \sqrt g = \frac{\Delta_0}{2} + \epsilon \frac{\cos\a \Delta_1}{2},
!l!t!x       \ee
!l!t!x       and the inverse Jacobian is
!l!t!x       \be (\sqrt g)^{-1} \approx \frac{2}{\Delta_0} - \epsilon \frac{2 \cos\a \Delta_1}{\Delta_0^2}.
!l!t!x       \ee
!l!t!x \item Assume the solution takes the form
!l!t!x       \be \begin{array}{cccccccccccccccccc}
!l!t!x           A_\a & = &   & \psi(\s) & + & \epsilon \, f(\s)\cos\a, \\
!l!t!x           A_\z & = & - & \chi(\s) & + & \epsilon \, g(\s)\cos\a.
!l!t!x       \end{array} \ee
!l!t!x \item The magnetic field, ${\bf B}\equiv\nabla\times{\bf A}$, is
!l!t!x       \be B^\s & = &                                                               - \epsilon \frac{ 2 g \, \sin \a}{\Delta_0}, \label{eq:radialfield} \\
!l!t!x           B^\a & = & \frac{2 \chi^\prime}{\Delta_0} - \epsilon \frac{2 \chi^\prime \cos\a \Delta_1}{\Delta_0^2} - \epsilon \frac{2 g^\prime \cos \a}{\Delta_0}, \\
!l!t!x           B^\z & = & \frac{2 \psi^\prime}{\Delta_0} - \epsilon \frac{2 \psi^\prime \cos\a \Delta_1}{\Delta_0^2} + \epsilon \frac{2 f^\prime \cos \a}{\Delta_0}. \label{eq:toroidalfield}
!l!t!x       \ee
!l!t!x       (The equation for $B^\z$ will later be shown to simplify to $B^\z = 2 \psi^\prime / \Delta_0 + \epsilon \mu g \cos \t $.)
!l!t!x \item Force balance, $\nabla \times {\bf B} = \mu \nabla \times {\bf A}$, assuming $\partial_\z\equiv0$, reduces to
!l!tex       \be \begin{array}{lllccccccccccccccc}
!l!tex           \partial_\t \; B^\z                &   &                                              & = & \mu \, \partial_\t A_\z &   & \\
!l!tex                                              & - & \partial_\s \; B^\z                          & = &                         & - & \mu \, \partial_\s A_\z \\
!l!tex           \partial_\s (B^\s g_{\s\t} + B^\t) & - & \partial_\t ( B^\s g_{\s\s} + B^\t g_{\s\t}) & = & \mu \, \partial_\s A_\t &   &  
!l!tex       \end{array} \ee
!l!tex \item This becomes
!l!t!x       \be \epsilon \sin \a \left[ \left(\frac{2 \psi^\prime \Delta_1}{\Delta_0^2} - \frac{2 f^\prime}{\Delta_0}\right) \left(\frac{m^2+n^2}{m^2}\right) + \left(\frac{2 \chi^\prime \Delta_1}{\Delta_0^2} + \frac{2 g^\prime}{\Delta_0}\right)\left(\frac{n}{m}\right)\right] & = & -\mu \epsilon g \sin \t, \\
!l!t!x         - \frac{2 \psi^{\prime\prime} }{\Delta_0} + \epsilon \frac{2 \psi^{\prime\prime} \cos \t \Delta_1}{\Delta_0^2} - \epsilon \frac{2 f^{\prime\prime} \cos \t}{\Delta_0} & = & \mu \chi^{\prime} - \mu \epsilon g^\prime \cos\t, \\
!l!t!x          \frac{2 \chi^{\prime\prime}}{\Delta_0} - \epsilon \frac{2 \chi^{\prime\prime} \cos \t \Delta_1}{\Delta_0^2} - \epsilon \frac{2 g^{\prime\prime} \cos\t}{\Delta_0} + \epsilon \frac{g \cos\t \Delta_0}{2} + \epsilon \frac{\chi^{\prime} \cos \t [(R_{v,1}+R_{v-1,1})+\Delta_1s]}{2} & = & \mu \psi^{\prime} + \mu \epsilon f^{\prime} \cos \t.
!l!t!x       \ee
!l!t!x \item The ${\cal O}(\epsilon^0)$ part is identical to the unperturbed case, as described by \Eqn{unperturbed}.
!l!t!x \item The ${\cal O}(\epsilon^1)$ part is
!l!t!x       \be \psi^{\prime}\frac{\Delta_1}{\Delta_0} - f^{\prime} & = & - \bar \mu g \label{eq:aa}\\
!l!t!x           \psi^{\prime\prime}\frac{\Delta_1}{\Delta_0}-f^{\prime\prime} & = & -\bar \mu g^\prime \label{eq:bb}\\
!l!t!x          -\chi^{\prime\prime}\frac{\Delta_1}{\Delta_0^2}-g^{\prime\prime}\frac{2}{\Delta_0}+g\frac{\Delta_0}{2}+\chi^{\prime}\frac{(R_{v,1}+R_{v-1,1})+\Delta_1 s}{2} & = & \mu f^\prime.\label{eq:cc}
!l!t!x       \ee
!l!t!x \item Using \Eqn{aa} to write $f^\prime = \bar \mu g + \psi^\prime \Delta_1 / \Delta_0$, \Eqn{cc} becomes a second-order, ordinary differential equation for $g$, 
!l!t!x       \be g^{\prime\prime} + g \left(\bar\mu^2-\frac{\Delta_0^2}{4}\right) = \chi^\prime \frac{\left[ ( R_{v,1} + R_{v-1,1} ) + \Delta_1 s \right]\Delta_0}{4} - 2 \chi^{\prime\prime} \frac{\Delta_1}{\Delta_0},
!l!t!x       \label{eq:perturbedgeneral}
!l!t!x       \ee
!l!t!x       where we have used \Eqn{unperturbed} to write $\bar\mu\psi^\prime = \chi^{\prime\prime}$.
!l!t!x \item This is to be solved subject to the boundary conditions $g(-1)=g(+1)=0$.
!l!t!x       (These boundary conditions are derived from requiring that the field normal to the interfaces, i.e. $B^s$ as given by \Eqn{radialfield}, is zero.)
!l!t!x \item With $\chi$ given by \Eqn{unperturbedsolution}, and assuming $\lambda^2 \equiv \bar \mu^2 - \Delta_0^2/2^2 > 0$, the solution is
!l!t!x       \be g(s) & = & c \cos(\lambda \s) + d \sin(\lambda \s) + \bar g(s), \label{eq:perturbedfield}
!l!t!x       \ee
!l!t!x       where $\bar g(s)$ is the particular solution, 
!l!t!x       \be \bar g(s) \equiv \bar \mu \left[b\sin(\bar\mu s)-a\cos(\bar\mu s)\right] \left[(R_{v,1}+R_{v-1,1})+\Delta_1 s \right] / \Delta_0,
!l!t!x       \ee
!l!t!x       and
!l!t!x          $ c  =  \left[-\bar g(-1) - \bar g(+1) \right]/2 \cos \lambda$ and
!l!t!x          $ d  =  \left[+\bar g(-1) - \bar g(+1) \right]/2 \sin \lambda$,
!l!t!x       where it is required that $\sin\lambda\ne 0$ and $\cos\lambda\ne 0$.
!l!t!x \item Using \Eqn{aa}, the expression for the toroidal field, $B^\z$, as given by \Eqn{toroidalfield} reduces to 
!l!t!x       \be B^\z & = & \frac{2 \psi^\prime}{\Delta_0} + \epsilon \cos \t \, \mu g.
!l!t!x       \ee
!l!t!x       This shows that there is no variation in the toroidal field on the interfaces.
!l!t!x \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!l!t!x \subsubsection{rotational-transform}

!l!t!x \begin{enumerate}

!l!t!x \item The `local' rotational-transform, $\dot \t \equiv B^\t / B^\z$, is not constant on the inner and outer interfaces,
!l!t!x       \be \dot \t = \frac{\chi^\prime}{\psi^\prime} - \epsilon \cos \t \; h \label{eq:localtransform}
!l!t!x       \ee
!l!t!x       where $h\equiv (g^\prime/\psi^\prime+\chi^\prime \Delta_1 / \psi^\prime \Delta_0 )$.
!l!t!x \item To calculate the rotational-transform, it is convenient to introduce a straight-fieldline angle, $\t_s \equiv \t + \epsilon \lambda \sin\t$,
!l!t!x       where $\lambda$ is to be constrained by requiring that
!l!t!x       \be \frac{{\bf B}\cdot\nabla\t_s}{{\bf B}\cdot\nabla\z} = \dot \t \; ( 1 + \epsilon \lambda \cos \t) = const.
!l!t!x       \ee
!l!t!x       where the constant will be identified as the rotational-transform.
!l!t!x \item With the expression for $\dot\t$ given in \Eqn{localtransform}, this becomes
!l!t!x       \be \iotabar = \frac{\chi^\prime}{\psi^\prime} + \epsilon \cos \t \frac{\chi^\prime}{\psi^\prime} \lambda - \epsilon \cos \t \,  h.
!l!t!x       \ee
!l!t!x \item The solution is $\lambda \equiv g'/\chi'+\Delta_1 / \Delta_0$, and the rotational-transform is unchanged from the unperturbed case, i.e. $\iotabar=\chi'/\psi'$.

!l!t!x \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!l!t!x \subsubsection{field strength}

!l!t!x \begin{enumerate}

!l!t!x \item The magnetic field strength on the interfaces is 
!l!t!x       \be B^2 & = & \frac{4}{\Delta_0^2} ( \chi^{\prime 2} + \psi^{\prime 2} ) - \epsilon \frac{8}{\Delta_0^2} \cos\t \chi^\prime \left( \chi^{\prime} \frac{\Delta_1}{\Delta_0}+g^\prime\right).
!l!t!x       \ee
!l!t!x \item The $m=0$ component reduces to 
!l!t!x       \be B_{0}^{2} = \frac{4\bar\mu^2(a^2+b^2)}{\Delta_0^2}.
!l!t!x       \ee

!l!t!x \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!l!t!x    \section{cylindrical geometry - general}

!l!t!x \begin{enumerate}

!l!t!x \item Position is given by
!l!t!x       \be {\bf x} \equiv R(\s,\t,\z) \cos(\t) \; {\bf \hat i} + R(\s,\t,\z) \sin(\t) \; {\bf \hat j} + \z \; {\bf \hat k},
!l!t!x       \ee
!l!t!x        where $s\in[-1,+1]$.
!l!t!x \item The (lower) metric elements and Jacobian are
!l!t!x       \be g_{\s\s} & = & R_\s R_\s      , \;\; g_{\s\t} = R_\s R_\t, \;\; g_{\s\z} = R_\s R_\z, \;\;
!l!t!x           g_{\t\t}   =   R_\t R_\t + R^2, \;\; g_{\t\z} = R_\t R_\z, \;\;
!l!t!x           g_{\z\z}   =   R_\s R_\s + 1  ; \\
!l!t!x           \sqrt g  & = & R R_\s.
!l!t!x       \ee
!l!t!x \item The (upper) metric elements are
!l!t!x       \be g^{\s\s} &=& \frac{R_\t R_\t + R^2 R_\z R_\z + R^2}{R^2 R_\s^2} \\
!l!t!x           g^{\s\t} &=& \frac{   - R_\t                      }{R   R_\s  } \\
!l!t!x           g^{\s\z} &=& \frac{   - R_\z                      }{    R_\s  } \\
!l!t!x           g^{\t\t} &=& \frac{        1                      }{R^2       } \\
!l!t!x           g^{\t\z} &=&               0                                    \\
!l!t!x           g^{\z\z} &=& \frac{        1                      }{    R_\s^2}
!l!t!x       \ee

!l!t!x \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!l!t!x \subsection{cylindrical geometry - general representation for magnetic field}

!l!t!x \begin{enumerate}

!l!t!x \item As usual, the magnetic field ${\bf B}\equiv \nabla \times \left( A_\t \nabla \t + A_\z \nabla \z \right)$ is
!l!t!x       \be B^\s & = & (\partial_\t A_\z        - \partial_\z A_\t       )/\sqrt g \\
!l!t!x           B^\t & = & (\;\;\;\;\;\;\;  \;      -  \;\;\;\;   A_\z^\prime)/\sqrt g \\
!l!t!x           B^\z & = & (  \;\;\;\;  A_\t^\prime    \;\;\;\;\;\;\;\;\;\;\;)/\sqrt g
!l!t!x       \ee
!l!t!x       \be B_\s & = & B^\s \; g_{\s\s} + B^\t \; g_{\s\t} + B^\z \; g_{\s\z} \\
!l!t!x           B_\t & = & B^\s \; g_{\s\t} + B^\t \; g_{\t\t} + B^\z \; g_{\t\z} \\
!l!t!x           B_\z & = & B^\s \; g_{\s\z} + B^\t \; g_{\t\z} + B^\z \; g_{\z\z}
!l!t!x       \ee

!l!t!x \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!l!t!x \subsection{cylindrical geometry - Beltrami equation}

!l!t!x \begin{enumerate}

!l!t!x \item The Beltrami equation reduces to the completely-general form:
!l!t!x       \be \mu \sqrt g B^\s &=& \partial_\t B_\z        - \partial_\z B_\t        \\
!l!t!x           \mu \sqrt g B^\t &=& \partial_\z B_\s        - \;\;\;\;    B_\z^\prime \\
!l!t!x           \mu \sqrt g B^\z &=& \;\;\;\;    B_\t^\prime - \partial_\t B_\s
!l!t!x       \ee

!l!t!x \item In terms of the vector potential, this becomes
!l!t!x       \be \begin{array}{ccccccccccccccccccccc}
!l!t!x        \mu \partial_\t A_\z - \mu \partial_\z A_\t &=&
!l!t!x     \partial_\t \partial_\t A_\z        \;             \bar g_{\s\z} 
!l!t!x &-& \partial_\t \partial_\z A_\t        \;             \bar g_{\s\z} 
!l!t!x &+&             \partial_\t A_\z        \; \partial_\t \bar g_{\s\z} 
!l!t!x &-&             \partial_\z A_\t        \; \partial_\t \bar g_{\s\z} \\
!l!t!x &-& \partial_\t             A_\z^\prime \;             \bar g_{\t\z}
!l!t!x &-&                         A_\z^\prime \; \partial_\t \bar g_{\t\z}
!l!t!x &+& \partial_\t             A_\t^\prime \;             \bar g_{\z\z}
!l!t!x &+&                         A_\t^\prime \; \partial_\t \bar g_{\z\z} \\
!l!t!x &-& \partial_\z \partial_\t A_\z        \;             \bar g_{\s\t} 
!l!t!x &-&             \partial_\t A_\z        \; \partial_\z \bar g_{\s\t} 
!l!t!x &+& \partial_\z \partial_\z A_\t        \;             \bar g_{\s\t} 
!l!t!x &+&             \partial_\z A_\t        \; \partial_\z \bar g_{\s\t} \\
!l!t!x &+& \partial_\z             A_\z^\prime \;             \bar g_{\t\t}
!l!t!x &+&                         A_\z^\prime \; \partial_\z \bar g_{\t\t}
!l!t!x &-& \partial_\z             A_\t^\prime \;             \bar g_{\t\z}
!l!t!x &-&                         A_\t^\prime \; \partial_\z \bar g_{\t\z},
!l!t!x \end{array}
!l!t!x       \ee
!l!t!x       \be \begin{array}{ccccccccccccccccccccc}
!l!t!x                             - \mu             A_\z^\prime &=&
!l!t!x     \partial_\z \partial_\t A_\z        \;             \bar g_{\s\s} 
!l!t!x &-& \partial_\z \partial_\z A_\t        \;             \bar g_{\s\s} 
!l!t!x &+&             \partial_\t A_\z        \; \partial_\z \bar g_{\s\s} 
!l!t!x &-&             \partial_\z A_\t        \; \partial_\z \bar g_{\s\s} \\
!l!t!x &-& \partial_\z             A_\z^\prime \;             \bar g_{\s\t}
!l!t!x &-&                         A_\z^\prime \; \partial_\z \bar g_{\s\t}
!l!t!x &+& \partial_\z             A_\t^\prime \;             \bar g_{\s\z}
!l!t!x &+&                         A_\t^\prime \; \partial_\z \bar g_{\s\z} \\
!l!t!x &-& \partial_\s \partial_\t A_\z        \;             \bar g_{\s\z} 
!l!t!x &-&             \partial_\t A_\z        \; \partial_\s \bar g_{\s\z} 
!l!t!x &+& \partial_\s \partial_\z A_\t        \;             \bar g_{\s\z} 
!l!t!x &+&             \partial_\z A_\t        \; \partial_\s \bar g_{\s\z} \\
!l!t!x &+& \partial_\s             A_\z^\prime \;             \bar g_{\t\z}
!l!t!x &+&                         A_\z^\prime \; \partial_\s \bar g_{\t\z}
!l!t!x &-& \partial_\s             A_\t^\prime \;             \bar g_{\z\z}
!l!t!x &-&                         A_\t^\prime \; \partial_\s \bar g_{\z\z},
!l!t!x \end{array}
!l!t!x       \ee
!l!t!x       \be \begin{array}{ccccccccccccccccccccc}
!l!t!x                             + \mu             A_\t^\prime &=&
!l!t!x     \partial_\s \partial_\t A_\z        \;             \bar g_{\s\t} 
!l!t!x &-& \partial_\s \partial_\z A_\t        \;             \bar g_{\s\t} 
!l!t!x &+&             \partial_\t A_\z        \; \partial_\s \bar g_{\s\t} 
!l!t!x &-&             \partial_\z A_\t        \; \partial_\s \bar g_{\s\t} \\
!l!t!x &-& \partial_\s             A_\z^\prime \;             \bar g_{\t\t}
!l!t!x &-&                         A_\z^\prime \; \partial_\s \bar g_{\t\t}
!l!t!x &+& \partial_\s             A_\t^\prime \;             \bar g_{\t\z}
!l!t!x &+&                         A_\t^\prime \; \partial_\s \bar g_{\t\z} \\
!l!t!x &-& \partial_\t \partial_\t A_\z        \;             \bar g_{\s\s} 
!l!t!x &-&             \partial_\t A_\z        \; \partial_\t \bar g_{\s\s} 
!l!t!x &+& \partial_\t \partial_\z A_\t        \;             \bar g_{\s\s} 
!l!t!x &+&             \partial_\z A_\t        \; \partial_\t \bar g_{\s\s} \\
!l!t!x &+& \partial_\t             A_\z^\prime \;             \bar g_{\s\t}
!l!t!x &+&                         A_\z^\prime \; \partial_\t \bar g_{\s\t}
!l!t!x &-& \partial_\t             A_\t^\prime \;             \bar g_{\s\z}
!l!t!x &-&                         A_\t^\prime \; \partial_\t \bar g_{\s\z},
!l!t!x \end{array}
!l!t!x       \ee
!l!t!x       where $\bar g_{\a\b}\equiv g_{\a\b}/\sqrt g $.

!l!t!x \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!l!t!x \subsection{cylindrical geometry - coordinate interpolation}

!l!t!x \begin{enumerate}

!l!t!x \item The coordinate harmonics are interpolated according to
!l!t!x       \be R_j(s) &=& R_j(1) \; \bar s^{1/2}\;, \; \mbox{\rm for } \; m_j=0,\\
!l!t!x           R_j(s) &=& R_j(1) \; \bar s^{m_j/2}, \; \mbox{\rm for } \; m_j\ne0,
!l!t!x       \ee
!l!t!x       where $\bar s \equiv (s+1)/2$.

!l!t!x \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!l!t!x \subsection{cylindrical geometry - axisymmetric}

!l!t!x \begin{enumerate}

!l!t!x \item For axisymmetry, the derivatives w.r.t. $\t$ and $\z$ vanish, and 
!l!t!x       \be - \mu A_\z^\prime &=& - A_\t^{\prime\prime} \bar g_{\z\z} - A_\t^\prime \bar g_{\z\z}^\prime \\
!l!t!x           + \mu A_\t^\prime &=& - A_\z^{\prime\prime} \bar g_{\t\t} - A_\z^\prime \bar g_{\t\t}^\prime
!l!t!x       \ee

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine verify
 
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero

  use numerical, only : 

  use fileunits, only : ounit

  use inputlist, only : Wverify

  use cputiming, only : Tverify

  use allglobal, only : myid, cpus
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS

  BEGIN(verify)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  
  
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  RETURN(verify)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine verify

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
