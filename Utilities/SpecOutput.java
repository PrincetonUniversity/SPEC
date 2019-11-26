package ipp.w7x.equi;
// AUTO-GENERATED; DO NOT COMMIT CHANGES TO THIS FILE !
// auto-created by a user called 'jonathan' on a machine called 'Nebuchadnezaar' at 26/11/2019 15:27:37 UTC

import java.io.IOException;
import java.util.Locale;
import ucar.ma2.DataType;
import ucar.nc2.NetcdfFile;

/** Output file of the SPEC MHD equilibrium code */
public class SpecOutput {

	/** The namelist \verb+physicslist+ controls the geometry, profiles, and numerical resolution. */
	public static class InputPhysics {
		/** selects Cartesian, cylindrical or toroidal geometry
		\bi
		\item[i.] \inputvar{Igeometry = 1} : Cartesian; geometry determined by $R$
		\item[i.] \inputvar{Igeometry = 2} : cylindrical; geometry determined by $R$
		\item[i.] \inputvar{Igeometry = 3} : toroidal; geometry determined by $R$ {\em and} $Z$
		\ei */
		public int Igeometry;
		/** stellarator symmetry is enforced if \inputvar{Istellsym.eq.1} */
		public int Istellsym;
		/** compute vacuum field surrounding plasma */
		public int Lfreebound;
		/** total enclosed toroidal magnetic flux */
		public double phiedge;
		/** total enclosed (toroidal) plasma current */
		public double curtor;
		/** total enclosed (poloidal) linking current */
		public double curpol;
		/** adiabatic index; cannot set $|\gamma| = 1$ */
		public double gamma;
		/** field periodicity
		\bi
		\item[i.] all Fourier representations are of the form $\cos(m\t-nN\z)$, $\sin(m\t-nN\z)$,where $N\equiv$\inputvar{Nfp}
		\item[i.] constraint : \inputvar{Nfp.ge.1}
		\ei */
		public int Nfp;
		/** number of volumes
		\bi
		\item[i.] each volume ${\cal V}_l$ is bounded by the ${\cal I}_{l-1}$ and ${\cal I}_{l}$ interfaces
		\item[i.] note that in cylindrical or toroidal geometry, ${\cal I}_{0}$ is the degenerate coordinate axis
		\item[i.] constraint : \inputvar{Nvol.le.MNvol}
		\ei */
		public int Nvol;
		/** number of poloidal Fourier harmonics */
		public int Mpol;
		/** number of toroidal Fourier harmonics */
		public int Ntor;
		/** Chebyshev resolution in each volume
		\bi
		\item[i.] constraint : \inputvar{Lrad(1:Mvol)}.ge.2
		\ei */
		public int[] Lrad;
		/** selects constraints; primarily used in \link{ma02aa} and \link{mp00ac}.
		\bi
		\item[i.]   if \inputvar{Lconstraint}.eq.-1, then in the plasma regions $\Delta\psi_t$, $\mu$ and $\Delta \psi_p$ are {\em not} varied;
		            and in the vacuum region (only for free-boundary) $\Delta\psi_t$ and $\Delta \psi_p$ are {\em not} varied, and $\mu = 0$.
		\item[ii.]  if \inputvar{Lconstraint}.eq.0, then in the plasma regions $\Delta\psi_t$, $\mu$ and $\Delta \psi_p$ are {\em not} varied;
		            and in the vacuum region (only for free-boundary) $\Delta\psi_t$ and $\Delta \psi_p$ are varied to match the 
		            prescribed plasma current, \inputvar{curtor}, and the ``linking'' current, \inputvar{curpol}, and $\mu = 0$;
		\item[iii.] if \inputvar{Lconstraint}.eq.1, then in the plasma regions $\mu$ and $\Delta\psi_p$ are adjusted
		            in order to satisfy the inner and outer interface transform constraints
		            (except in the simple torus, where the enclosed poloidal flux is irrelevant,
		            and only $\mu$ is varied to satisfy the outer interface transform constraint);
		            and in the vacuum region $\Delta\psi_t$ and $\Delta \psi_p$ are varied to match the transform constraint on the boundary
		            and to obtain the prescribed linking current, \inputvar{curpol}, and $\mu = 0$.
		\item[iv.]  if \inputvar{Lconstraint}.eq.2, under reconstruction.
		\ei */
		public int Lconstraint;
		/** toroidal flux, $\psi_t$, enclosed by each interface
		\bi
		\item[i.] For each of the plasma volumes, this is a constraint: \inputvar{tflux} is \emph{not} varied
		\item[i.] For the vacuum region (only if \inputvar{Lfreebound = 1}), \inputvar{tflux} may be allowed to vary to match constraints
		\item[i.] Note that \inputvar{tflux} will be normalized so that \inputvar{tflux(Nvol) = 1.0},
		          so that \inputvar{tflux} is arbitrary up to a scale factor
		\item[i.] see also \inputvar{phiedge}
		\ei */
		public double[] tflux;
		/** poloidal flux, $\psi_p$, enclosed by each interface */
		public double[] pflux;
		/** helicity, ${\cal K}$, in each volume, ${\cal V}_i$
		\bi
		\item[i.] on exit, \inputvar{helicity} is set to the computed values of ${\cal K} \equiv \int {\bf A}\cdot{\bf B}\;dv$
		\ei */
		public double[] helicity;
		/** pressure scale factor
		\bi
		\item[i.] the initial pressure profile is given by \inputvar{pscale} $*$ \inputvar{press}
		\ei */
		public double pscale;
		/** pressure in each volume
		\bi
		\item[i.] the pressure is {\em not} held constant, but $p_l V_l^\gamma = P_l$ {\em is} held constant,
		          where $P_l$ is determined by the initial pressures and the initial volumes, $V_l$
		\item[i.] (Note that if \inputvar{gamma = 0.0}, then $p_l \equiv P_l$.)
		\item[i.] on output, the pressure is given by $p_l = P_l/V_l^\gamma$, where $V_l$ is the final volume
		\item[i.] \inputvar{pressure} is only used in calculation of interface force-balance
		\ei */
		public double[] pressure;
		/** logical flag
		\bi
		\item[i.] if \inputvar{Ladiabatic = 0}, the adiabatic constants are determined by the initial pressure and volume
		\item[i.] if \inputvar{Ladiabatic = 1}, the adiabatic constants are determined by the given input \inputvar{adiabatic}
		\ei */
		public int Ladiabatic;
		/** adiabatic constants in each volume
		\bi
		\item[i.] the pressure is {\em not} held constant, but $p_l V_l^\gamma = P_l \equiv$\inputvar{adiabatic} is constant,
		\item[i.] note that if \inputvar{gamma = 0.0}, then \inputvar{pressure = adiabatic}
		\item[i.] \inputvar{pressure} is only used in calculation of interface force-balance
		\ei */
		public double[] adiabatic;
		/** helicity-multiplier, $\mu$, in each volume */
		public double[] mu;
		/** part of noble irrational defining iota
		\bi
		\item[i.] ``inside'' interface rotational-transform is $\iotabar = (p_l+\gamma p_r)/(q_l+\gamma q_r)$,
		          where $\gamma$ is the golden mean, $\gamma = (1 + \sqrt 5 ) / 2 $
		\item[i.] if both $q_l = 0$ {\em and} $q_r = 0$, then the (inside) interface rotational-transform is defined by \inputvar{iota}
		\ei */
		public int[] pl;
		/** part of noble irrational defining iota
		\bi
		\item[i.] ``inside'' interface rotational-transform is $\iotabar = (p_l+\gamma p_r)/(q_l+\gamma q_r)$,
		          where $\gamma$ is the golden mean, $\gamma = (1 + \sqrt 5 ) / 2 $
		\item[i.] if both $q_l = 0$ {\em and} $q_r = 0$, then the (inside) interface rotational-transform is defined by \inputvar{iota}
		\ei */
		public int[] ql;
		/** part of noble irrational defining iota
		\bi
		\item[i.] ``inside'' interface rotational-transform is $\iotabar = (p_l+\gamma p_r)/(q_l+\gamma q_r)$,
		          where $\gamma$ is the golden mean, $\gamma = (1 + \sqrt 5 ) / 2 $
		\item[i.] if both $q_l = 0$ {\em and} $q_r = 0$, then the (inside) interface rotational-transform is defined by \inputvar{iota}
		\ei */
		public int[] pr;
		/** part of noble irrational defining iota
		\bi
		\item[i.] ``inside'' interface rotational-transform is $\iotabar = (p_l+\gamma p_r)/(q_l+\gamma q_r)$,
		          where $\gamma$ is the golden mean, $\gamma = (1 + \sqrt 5 ) / 2 $
		\item[i.] if both $q_l = 0$ {\em and} $q_r = 0$, then the (inside) interface rotational-transform is defined by \inputvar{iota}
		\ei */
		public int[] qr;
		/** rotational-transform, $\iotabar$, on inner side of each interface
		\bi
		\item[i.] only relevant if \inputvar{ql}=0 and \inputvar{qr}=0
		\ei */
		public double[] iota;
		/** part of noble irrational defining oita
		\bi
		\item "outer" interface rotational-transform is $\iotabar = (p_l+\gamma p_r)/(q_l+\gamma q_r)$,
		      where $\gamma$ is the golden mean, $\gamma = (1 + \sqrt 5 ) / 2 $
		\item if both $q_l = 0$ {\em and} $q_r = 0$, then the (outer) interface rotational-transform is defined by \inputvar{oita}
		\ei */
		public int[] lp;
		/** part of noble irrational defining oita
		\bi
		\item "outer" interface rotational-transform is $\iotabar = (p_l+\gamma p_r)/(q_l+\gamma q_r)$,
		      where $\gamma$ is the golden mean, $\gamma = (1 + \sqrt 5 ) / 2 $
		\item if both $q_l = 0$ {\em and} $q_r = 0$, then the (outer) interface rotational-transform is defined by \inputvar{oita}
		\ei */
		public int[] lq;
		/** part of noble irrational defining oita
		\bi
		\item "outer" interface rotational-transform is $\iotabar = (p_l+\gamma p_r)/(q_l+\gamma q_r)$,
		      where $\gamma$ is the golden mean, $\gamma = (1 + \sqrt 5 ) / 2 $
		\item if both $q_l = 0$ {\em and} $q_r = 0$, then the (outer) interface rotational-transform is defined by \inputvar{oita}
		\ei */
		public int[] rp;
		/** part of noble irrational defining oita
		\bi
		\item "outer" interface rotational-transform is $\iotabar = (p_l+\gamma p_r)/(q_l+\gamma q_r)$,
		      where $\gamma$ is the golden mean, $\gamma = (1 + \sqrt 5 ) / 2 $
		\item if both $q_l = 0$ {\em and} $q_r = 0$, then the (outer) interface rotational-transform is defined by \inputvar{oita}
		\ei */
		public int[] rq;
		/** rotational-transform, $\iotabar$, on outer side of each interface
		\bi
		\item only relevant if illogical input for \inputvar{ql} and \inputvar{qr} are provided
		\ei */
		public double[] oita;
		/** accuracy to which $\mu$ and $\Delta\psi_p$ are required
		\bi
		\item only relevant if constraints on transform, enclosed currents etc. are to be satisfied iteratively, see \inputvar{Lconstraint}
		\ei */
		public double mupftol;
		/** an upper limit on the transform/helicity constraint iterations
		\item only relevant if constraints on transform, enclosed currents etc. are to be satisfied iteratively, see \inputvar{Lconstraint}
		\bi
		\item constraint: \inputvar{mupfits > 0}
		\ei */
		public int mupfits;
		/** stellarator symmetric coordinate axis R cosine Fourier coefficients */
		public double[] Rac;
		/** stellarator symmetric coordinate axis Z sine Fourier coefficients */
		public double[] Zas;
		/** non-stellarator symmetric coordinate axis R sine Fourier coefficients */
		public double[] Ras;
		/** non-stellarator symmetric coordinate axis Z cosine Fourier coefficients */
		public double[] Zac;
		/** stellarator symmetric boundary R cosine Fourier coefficients */
		public double[][] Rbc;
		/** stellarator symmetric boundary Z sine Fourier coefficients */
		public double[][] Zbs;
		/** non-stellarator symmetric boundary R sine Fourier coefficients */
		public double[][] Rbs;
		/** non-stellarator symmetric boundary Z cosine Fourier coefficients */
		public double[][] Zbc;
		/** stellarator symmetric boundary R cosine Fourier coefficients of wall */
		public double[][] Rwc;
		/** stellarator symmetric boundary Z sine Fourier coefficients of wall */
		public double[][] Zws;
		/** non-stellarator symmetric boundary R sine Fourier coefficients of wall */
		public double[][] Rws;
		/** non-stellarator symmetric boundary Z cosine Fourier coefficients of wall */
		public double[][] Zwc;
		/** stellarator symmetric normal field at boundary; vacuum component */
		public double[][] Vns;
		/** stellarator symmetric normal field at boundary; plasma component */
		public double[][] Bns;
		/** non-stellarator symmetric normal field at boundary; vacuum component */
		public double[][] Vnc;
		/** non-stellarator symmetric normal field at boundary; plasma component */
		public double[][] Bnc;
	} // end of InputPhysics

	/** The namelist \verb+numericlist+ controls internal resolution parameters that the user rarely needs to consider. */
	public static class InputNumerics {
		/** to initialize geometry using a regularization / extrapolation method
		\bi
		\item if \inputvar{Linitialize = -I}, where $I$ is a positive integer, 
		      the geometry of the $i=1,N_V-I$ surfaces constructed by an extrapolation
		\item if \inputvar{Linitialize = 0}, the geometry of the interior surfaces is provided after the namelists in the input file
		\item if \inputvar{Linitialize = 1}, the interior surfaces will be intialized as $R_{l,m,n} = R_{N,m,n} \psi_{t,l}^{m/2}$,
		      where $R_{N,m,n}$ is the plasma boundary
		      and $\psi_{t,l}$ is the given toroidal flux enclosed by the $l$-th interface, normalized to the total enclosed toroidal flux
		      a similar extrapolation is used for $Z_{l,m,n}$
		\item note that the Fourier harmonics of the boundary is {\em always} given by the \inputvar{Rbc} and \inputvar{Zbs} 
		      given in \type{physicslist}
		\item if \inputvar{Linitialize = 2}, the interior surfaces {\em and the plasma boundary} will be intialized
		      as $R_{l,m,n} = R_{W,m,n} \psi_{t,l}^{m/2}$, where $R_{W,m,n}$ is the computational boundary
		      and $\psi_{t,l}$ is the given toroidal flux enclosed by the $l$-th interface, normalized to the total enclosed toroidal flux
		      a similar extrapolation is used for $Z_{l,m,n}$
		\item note that, for free-boundary calculations, the Fourier harmonics of the computational boundary
		      is {\em always} given by the \inputvar{Rwc} and \inputvar{Zws}
		      given in \type{physicslist}
		\item if \inputvar{Linitialize = 1, 2}, it is not required to provide the geometry of the interfaces after the namelists
		\ei */
		public int Linitialize;
		/** to adjust vacuum field to cancel plasma field on computational boundary
		\bi
		\item only relevant if \inputvar{Lfreebound = 1},
		\ei */
		public int Lzerovac;
		/** \bi
		\item resolution of the real space grid on which fast Fourier transforms are performed is given by \inputvar{Ndiscrete*Mpol*4}
		\item constraint \inputvar{Ndiscrete>0}
		\ei */
		public int Ndiscrete;
		/** the resolution of the Gaussian quadrature
		\bi
		\item the resolution of the Gaussian quadrature, $\ds \int \!\! f(s) ds = \sum_k \omega_k f(s_k)$,
		      in each volume is given by \internal{Iquad$_v$}, 
		\item \internal{Iquad$_v$} is set in \link{preset}.
		%      and depends on \inputvar{Nquad}, \inputvar{Lrad$_v$} and \inputvar{Mpol}.
		%\bi
		%\item if \inputvar{Nquad.gt.0},                                 then \internal{Iquad(vvol) =              Nquad};
		%\item if \inputvar{Nquad.le.0 and .not.Lcoordinatesingularity}, then \internal{Iquad(vvol) = 2*Lrad(vvol)-Nquad};
		%\item if \inputvar{Nquad.le.0 and      Lcoordinatesingularity}, then \internal{Iquad(vvol) = 2*Lrad(vvol)-Nquad+Mpol};
		%\ei
		%\item \internal{Iquad$_v$} is passed through to \link{ma00aa} to compute various volume integrals; 
		%      also see \link{jo00aa}, where \internal{Iquad$_v$} 
		%      is also used in computing the volume integrals of $||\nabla\times{\bf B} - \mu {\bf B}||$;
		\ei */
		public int Nquad;
		/** Fourier resolution of straight-fieldline angle on interfaces
		\bi
		\item the rotational-transform on the interfaces is determined by a transformation to the straight-fieldline angle,
		      with poloidal resolution given by \inputvar{iMpol}
		\item if \inputvar{iMpol.le.0}, then \inputvar{iMpol = Mpol - iMpol}
		\ei */
		public int iMpol;
		/** Fourier resolution of straight-fieldline angle on interfaces
		\bi
		\item the rotational-transform on the interfaces is determined by a transformation to the straight-fieldline angle,
		      with toroidal resolution given by \inputvar{iNtor}
		\item if \inputvar{iNtor.le.0}, then \inputvar{iNtor = Ntor - iNtor}
		\item if \inputvar{Ntor.eq.0}, then the toroidal resolution of the angle transformation is set \inputvar{lNtor = 0}.
		\ei */
		public int iNtor;
		/** controls method used to solve for rotational-transform on interfaces
		\bi
		\item if \inputvar{Lsparse = 0}, the transformation to the straight-fieldline angle is computed in Fourier space
		      using a dense matrix solver, \nag{}{F04AAF}
		\item if \inputvar{Lsparse = 1}, the transformation to the straight-fieldline angle is computed in real space
		      using a dense matrix solver, \nag{}{F04ATF}
		\item if \inputvar{Lsparse = 2}, the transformation to the straight-fieldline angle is computed in real space
		      using a sparse matrix solver, \nag{}{F11DEF}
		\item if \inputvar{Lsparse = 3}, the different methods for constructing the straight-fieldline angle are compared
		\ei */
		public int Lsparse;
		/** controls method used to solve for rotational-transform on interfaces
		only relevant if \inputvar{Lsparse = 0}
		\bi
		\item if \inputvar{Lsvdiota = 0}, use standard linear solver to construct straight fieldline angle transformation
		\item if \inputvar{Lsvdiota = 1}, use SVD method to compute rotational-transform
		\ei */
		public int Lsvdiota;
		/** controls iterative solution to sparse matrix
		arising in real-space transformation to the straight-fieldline angle
		only relevant if \inputvar{Lsparse.eq.2}; see \link{tr00ab} for details
		\bi
		\item if \inputvar{imethod = 1}, the method is \type{RGMRES}
		\item if \inputvar{imethod = 2}, the method is \type{CGS}
		\item if \inputvar{imethod = 3}, the method is \type{BICGSTAB}
		\ei */
		public int imethod;
		/** controls real-space grid resolution for constructing the straight-fieldline angle
		only relevant if \inputvar{Lsparse>0}
		determines order of finite-difference approximation to the derivatives
		\bi
		\item if \inputvar{iorder = 2}, 
		\item if \inputvar{iorder = 4}, 
		\item if \inputvar{iorder = 6}, 
		\ei */
		public int iorder;
		/** controls iterative solution to sparse matrix arising in real-space transformation
		to the straight-fieldline angle
		only relevant if \inputvar{Lsparse.eq.2}; see \link{tr00ab} for details
		\bi
		\item if \inputvar{iprecon = 0}, the preconditioner is `N'
		\item if \inputvar{iprecon = 1}, the preconditioner is `J'
		\item if \inputvar{iprecon = 2}, the preconditioner is `S'
		\ei */
		public int iprecon;
		/** tolerance required for iterative construction of straight-fieldline angle
		only relevant if \inputvar{Lsparse.ge.2} */
		public double iotatol;
		/** geometry of innermost interface is defined by extrapolation */
		public int Lextrap;
		/** maximum regularization factor
		\bi
		\item if \inputvar{Mregular.ge.2}, then \internal{regumm}$_i$ = \inputvar{Mregular} $/ 2 $ where \internal{m}$_i > $ \inputvar{Mregular}
		\ei */
		public int Mregular;
	} // end of InputNumerics

	/** The namelist \verb+locallist+ controls the construction of the Beltrami fields in each volume. */
	public static class InputLocal {
		/** \bi
		\item if \inputvar{LBeltrami = 1,3,5 or 7}, (SQP) then the Beltrami field in each volume is constructed
		      by minimizing the magnetic energy with the constraint of fixed helicity;
		      this is achieved by using sequential quadratic programming as provided by \nag{}{E04UFF}; 
		      this approach has the benefit (in theory) of robustly constructing minimum energy solutions
		      when multiple, i.e. bifurcated, solutions exist.
		\item if \inputvar{LBeltrami = 2,3,6 or 7}, (Newton) then the Beltrami fields are constructed by employing a standard Newton method
		      for locating an extremum of
		      $F\equiv \int B^2 dv - \mu (\int {\bf A}\cdot{\bf B}dv-{\cal K})$,
		      where $\mu$ is treated as an independent degree of freedom similar to the parameters describing the vector potential
		      and ${\cal K}$ is the required value of the helicity; 
		      this is the standard Lagrange multipler approach for locating the constrained minimum; 
		      this method cannot distinguish saddle-type extrema from minima, and which solution that will be obtained depends on the initial guess;
		\item if \inputvar{LBeltrami = 4,5,6 or 7}, (linear) it is assumed that the Beltrami fields are parameterized by $\mu$;
		      in this case, it is only required to solve $\nabla \times {\bf B} = \mu {\bf B}$ which reduces to a system of linear equations;
		      $\mu$ may or may not be adjusted iteratively, depending on \inputvar{Lconstraint},
		      to satisfy either rotational-transform or helicity constraints;
		\item for flexibility and comparison, each of the above methods can be employed; for example:
		      \bi
		      \item if \inputvar{LBeltrami = 1}, only the SQP    method will be employed;
		      \item if \inputvar{LBeltrami = 2}, only the Newton method will be employed;
		      \item if \inputvar{LBeltrami = 4}, only the linear method will be employed; 
		      \item if \inputvar{LBeltrami = 3}, the SQP and the Newton method are used;
		      \item if \inputvar{LBeltrami = 5}, the SQP and the linear method are used;
		      \item if \inputvar{LBeltrami = 6}, the Newton and the linear method are used;
		      \item if \inputvar{LBeltrami = 7}, all three methods will be employed;
		      \ei
		\ei */
		public int LBeltrami;
		/** controls how initial guess for Beltrami field is constructed;
		\bi
		\item only relevant for routines that require an initial guess for the Beltrami fields, such as the SQP and Newton methods,
		or the sparse linear solver;
		\item if \inputvar{Linitgues = 0}, the initial guess for the Beltrami field is trivial;
		\item if \inputvar{Linitgues = 1}, the initial guess for the Beltrami field is an integrable approximation;
		\item if \inputvar{Linitgues = 2}, the initial guess for the Beltrami field is read from file; 
		\item if \inputvar{Linitgues = 3}, the initial guess for the Beltrami field will be randomized with the maximum \inputvar{maxrndgues};
		\ei */
		public int Linitgues;
		/** the maximum random number of the Beltrami field if \inputvar{Linitgues = 3} */
		public double maxrndgues;
		/** redundant */
		public int Lposdef;
	} // end of InputLocal

	/** The namelist \verb+globallist+ controls the search for global force-balance */
	public static class InputGlobal {
		/** use Newton methods to find zero of force-balance, which is computed by \link{dforce};
		\bi
		\item[o.] if \inputvar{Lfindzero = 0}, then \link{dforce} is called once 
		          to compute the Beltrami fields consistent with the given geometry and constraints;
		\item[i.] if \inputvar{Lfindzero = 1}, then call
		          \nag{www.nag.co.uk/numeric/FL/manual19/pdf/C05/c05ndf_fl19.pdf}{C05NDF} (uses function values only),
		          which iteratively calls \link{dforce};
		\item[ii.] if \inputvar{Lfindzero = 2}, then call
		           \nag{www.nag.co.uk/numeric/FL/manual19/pdf/C05/c05pdf_fl19.pdf}{C05PDF} (uses derivative information),
		           which iteratively calls \link{dforce};
		\ei */
		public int Lfindzero;
		/** controls the weight factor, \type{BBweight}, in the force-imbalance harmonics;
		\bi
		\item[i.] \type{BBweight(i)} $\ds \equiv \inputvar{opsilon} \times \exp\left[-\inputvar{escale} \times (m_i^2+n_i^2) \right]$
		\item[ii.] defined in \link{preset}; used in \link{dforce};
		\item[iii.] also see \Eqn{forcebalancemn} below;
		\ei */
		public double escale;
		/** weighting of force-imbalance; 
		\bi
		\item[i.] used in \link{dforce}; also see \Eqn{forcebalancemn} below;
		\ei */
		public double opsilon;
		/** spectral condensation parameter; 
		\bi
		\item[i.] used in \link{preset} to define \type{mmpp(i)} $\equiv m_i^p$, where $p\equiv $ \inputvar{pcondense};
		\item[ii.] the angle freedom is exploited to minimize $\ds \inputvar{epsilon} \sum_{i} m_i^p (R_{i}^2+Z_{i}^2)$
		      with respect to tangential variations in the interface geometry;
		\item[ii.] also see \Eqn{spectralbalancemn} below;
		\ei */
		public double pcondense;
		/** weighting of spectral-width constraint
		\bi
		\item[i.] used in \link{dforce}; also see \Eqn{spectralbalancemn} below
		\ei */
		public double epsilon;
		/** ``star-like'' poloidal angle constraint radial exponential factor;
		used in \link{preset} to construct \type{sweight} */
		public double wpoloidal;
		/** weighting of ``star-like'' poloidal angle constraint;
		used in \link{preset} to construct \type{sweight}; */
		public double upsilon;
		/** required tolerance in force-balance error; only used as an initial check;
		\bi
		\item[i.]   if the initially supplied interfaces are consistent with force-balance to within \inputvar{forcetol},
		            then the geometry of the interfaces is not altered;
		\item[ii.]  if not, then the geometry of the interfaces is changed in order to bring the configuration into forcebalance
		            so that the geometry of interfaces is within \inputvar{c05xtol}, defined below, of the true solution;
		\item[iii.] to force execution of either \nag{www.nag.co.uk/numeric/FL/manual19/pdf/C05/c05ndf_fl19.pdf}{C05NDF}
		            or \nag{www.nag.co.uk/numeric/FL/manual19/pdf/C05/c05pdf_fl19.pdf}{C05PDF}, regardless of the initial force imbalance, 
		            set \inputvar{forcetol < 0};
		\ei */
		public double forcetol;
		/** required tolerance in position, ${\bf x} \equiv \{ R_{i,v}, Z_{i,v}\}$; */
		public double c05xmax;
		/** required tolerance in position, ${\bf x} \equiv \{ R_{i,v}, Z_{i,v}\}$;
		\bi
		\item[i.] used by both \nag{www.nag.co.uk/numeric/FL/manual19/pdf/C05/c05ndf_fl19.pdf}{C05NDF} and 
		          \nag{www.nag.co.uk/numeric/FL/manual19/pdf/C05/c05pdf_fl19.pdf}{C05PDF};
		          see the NAG documents for further details on how the error is defined;
		\item[ii.] constraint \inputvar{c05xtol.gt.0.0};
		\ei */
		public double c05xtol;
		/** used to control initial step size in
		      \nag{www.nag.co.uk/numeric/FL/manual19/pdf/C05/c05ndf_fl19.pdf}{C05NDF} and 
		      \nag{www.nag.co.uk/numeric/FL/manual19/pdf/C05/c05pdf_fl19.pdf}{C05PDF};
		\bi
		\item[i.] constraint \inputvar{c05factor.gt.0.0};
		\item[ii.] only relevant if \inputvar{Lfindzero.gt.0};
		\ei */
		public double c05factor;
		/** read $\nabla_{\bf x} {\bf F}$ from file \type{.GF}; 
		\bi
		\item[i.] only used if \inputvar{Lfindzero = 2}; 
		\item[ii.] only used in \link{newton};
		\ei */
		public boolean LreadGF;
		/** maximum allowed free-boundary iterations;
		\bi
		\item[i.] only used if \inputvar{Lfreebound = 1}; 
		\item[ii.] only used in \link{xspech};
		\ei */
		public int mfreeits;
		/** redundant */
		public double bnstol;
		/** redundant */
		public double bnsblend;
		/** equired tolerance in free-boundary iterations;
		\bi
		\item[i.] only used if \inputvar{Lfreebound = 1}; 
		\item[ii.] only used in \link{xspech}; see \link{xspech} for more documentation;
		\ei */
		public double gBntol;
		/** normal blend;
		\bi
		\item[i.] The ``new'' magnetic field at the computational boundary produced by the plasma currents is updated using a Picard scheme:
		          \be ({\bf B}\cdot{\bf n})^{j+1} =    \inputvar{gBnbld}  \times ({\bf B}\cdot{\bf n})^{j} 
		                                          + (1-\inputvar{gBnbld}) \times ({\bf B}\cdot{\bf n})^{*},
		          \ee
		          where $j$ labels free-boundary iterations, and $({\bf B}\cdot{\bf n})^{*}$ is computed by virtual casing.
		\item[ii.] only used if \inputvar{Lfreebound = 1}; 
		\item[ii.] only used in \link{xspech};
		\ei */
		public double gBnbld;
		/** minimum number of calls to adaptive virtual casing routine; see \link{casing} */
		public double vcasingeps;
		/** accuracy on virtual casing integral; see \link{bnorml}, \link{casing} */
		public double vcasingtol;
		public int vcasingits;
		public int vcasingper;
		public int mcasingcal;
	} // end of InputGlobal

	/** The namelist \type{diagnosticslist} controls post-processor diagnostics, such as \Poincare plot resolution, $\dots$,... */
	public static class InputDiagnostics {
		/** o.d.e. integration tolerance for all field line tracing routines */
		public double odetol;
		/** redundant */
		public double absreq;
		/** redundant */
		public double relreq;
		/** redundant */
		public double absacc;
		/** redundant */
		public double epsr;
		/** number of toroidal transits used (per trajectory) in following field lines
		for constructing \Poincare plots;
		if \inputvar{nPpts<1}, no \Poincare plot is constructed; */
		public int nPpts;
		/** number of trajectories in each annulus to be followed in constructing \Poincare plot;
		\bi
		\item if \inputvar{nPtrj(l)<0}, then \inputvar{nPtrj(l) = Ni(l)},
		      where \type{Ni(l)} is the grid resolution used to construct the Beltrami field in volume $l$;
		\ei */
		public int[] nPtrj;
		/** to compute eigenvalues of $\nabla {\bf F}$ */
		public boolean LHevalues;
		/** to compute eigenvectors (and also eigenvalues) of $\nabla {\bf F}$ */
		public boolean LHevectors;
		/** to compute and write to file the elements of $\nabla {\bf F}$ */
		public boolean LHmatrix;
		/** to compute linear, perturbed equilibrium */
		public int Lperturbed;
		/** perturbed harmonic */
		public int dpp;
		/** perturbed harmonic */
		public int dqq;
		/** implement various checks;
		\bi
		\item if \inputvar{Lcheck = 0}, no additional check on the calculation is performed;
		\item if \inputvar{Lcheck = 1}, the error in the current, i.e. $\nabla\times{\bf B}-\mu{\bf B}$ is computed as a post-diagnostic;
		\item if \inputvar{Lcheck = 2}, the analytic derivatives of the interface transform w.r.t.
		      the helicity multiplier, $\mu$, and the enclosed poloidal flux, $\Delta\psi_p$, are compared to a finite-difference estimate;
		\bi
		\item[i.] only if \inputvar{Lconstraint.eq.1};
		\item[ii.] only for \type{dspec} executable, i.e. must compile with \type{DFLAGS = "-D DEBUG"};
		\ei
		\item if \inputvar{Lcheck = 3}, the analytic derivatives of the volume w.r.t. interface Fourier harmonic
		      is compared to a finite-difference estimate;
		\bi
		\item[i.] must set \inputvar{Lfindzero}$ = 2$, 
		\item[ii.] set \inputvar{forcetol} sufficiently small and set \inputvar{LreadGF = F},
		      so that the matrix of second derivatives is calculated,
		\item[iii.] only for \type{dspec} executable, i.e. must compile with \type{DFLAGS = "-D DEBUG"};
		\ei
		\item if \inputvar{Lcheck = 4}, the analytic calculation of the derivatives of the magnetic field, $B^2$, at the interfaces
		      is compared to a finite-difference estimate;
		\bi
		\item[i.] must set \inputvar{Lfindzero}$ = 2$, 
		\item[ii.] set \inputvar{forcetol} sufficiently small,
		\item[iii.] set \inputvar{LreadGF=F},
		\item[iv.] only for \type{dspec} executable, i.e. must compile with \type{DFLAGS = "-D DEBUG"};
		\ei
		\item if \inputvar{Lcheck = 5}, the analytic calculation of the matrix of the derivatives of the force imbalance
		      is compared to a finite-difference estimate;
		\item if \inputvar{Lcheck = 6}, the virtual casing calculation is compared to \verb+xdiagno+;
		\bi
		\item[i.] the input file for \verb+xdiagno+ is written by \link{bnorml};
		\item[ii.] this provides the Cartesian coordinates on the computational boundary where the virtual casing routine \link{casing} 
		          computes the magnetic field, with the values of the magnetic field being written to the screen for comparison;
		\item[iii.] must set \inputvar{Freebound=1}, \inputvar{Lfindzero.gt.0}, \inputvar{mfreeits.ne.0};
		\item[iii.] \verb+xdiagno+; must be executed manually;
		\ei
		\ei */
		public int Lcheck;
		/** to check timing */
		public boolean Ltiming;
		/** redundant */
		public double fudge;
		/** redundant */
		public double scaling;
	} // end of InputDiagnostics

	/** input data for this SPEC run */
	public static class Input {
		/** initialize complex datatypes */
		public Input() {
			physics = new InputPhysics();
			numerics = new InputNumerics();
			local = new InputLocal();
			global = new InputGlobal();
			diagnostics = new InputDiagnostics();
		}

		/** The namelist \verb+physicslist+ controls the geometry, profiles, and numerical resolution. */
		public InputPhysics physics;
		/** The namelist \verb+numericlist+ controls internal resolution parameters that the user rarely needs to consider. */
		public InputNumerics numerics;
		/** The namelist \verb+locallist+ controls the construction of the Beltrami fields in each volume. */
		public InputLocal local;
		/** The namelist \verb+globallist+ controls the search for global force-balance */
		public InputGlobal global;
		/** The namelist \type{diagnosticslist} controls post-processor diagnostics, such as \Poincare plot resolution, $\dots$,... */
		public InputDiagnostics diagnostics;
	} // end of Input

	/** output data; content of the previous HDF5 file */
	public static class Output {
		/** stellarator symmetric normal field at boundary; vacuum component */
		public double[] Vns;
		/** stellarator symmetric normal field at boundary; plasma component */
		public double[] Bns;
		/** non-stellarator symmetric normal field at boundary; vacuum component */
		public double[] Vnc;
		/** non-stellarator symmetric normal field at boundary; plasma component */
		public double[] Bnc;
		/** number of Fourier modes */
		public int mn;
		/** poloidal mode numbers */
		public int[] im;
		/** toroidal mode numbers */
		public int[] in;
		/** number of interfaces = number of volumes */
		public int Mvol;
		/** stellarator symmetric Fourier harmonics, $R_{m,n}$, of interfaces */
		public double[][] Rbc;
		/** stellarator symmetric Fourier harmonics, $Z_{m,n}$, of interfaces */
		public double[][] Zbs;
		/** non-stellarator symmetric Fourier harmonics, $R_{m,n}$, of interfaces */
		public double[][] Rbs;
		/** non-stellarator symmetric Fourier harmonics, $Z_{m,n}$, of interfaces */
		public double[][] Zbc;
		/** force-balance error across interfaces */
		public double ForceErr;
		public double[] adiabatic;
		public double[] helicity;
		public double[] mu;
		public double[] tflux;
		public double[] pflux;
		public double volume;
		/** the maximum radial (Chebyshev) resolution */
		public int Mrad;
		/** the Chebyshev polynomials, $T_l$, and their derivatives, evaluated at $s=\pm 1$ */
		public double[][][] TT;
		/** the cosine harmonics of the covariant poloidal field,
		i.e. $[[B_{\t,j}]]$ evaluated on the inner and outer interface in each volume */
		public double[][][] Btemn;
		/** the cosine harmonics of the covariant toroidal field,
		i.e. $[[B_{\z,j}]]$ evaluated on the inner and outer interface in each volume */
		public double[][][] Bzemn;
		/** the sine harmonics of the covariant poloidal field,
		i.e. $[[B_{\t,j}]]$ evaluated on the inner and outer interface in each volume */
		public double[][][] Btomn;
		/** the sine harmonics of the covariant toroidal field,
		i.e. $[[B_{\t,j}]]$ evaluated on the inner and outer interface in each volume */
		public double[][][] Bzomn;
		/** resolution of the straight fieldline transformation */
		public double lmns;
	} // end of Output

	/** The covariant components of the vector potential are written as
	\be            A_\t & = & \sum_i \sum_{l=0}^L \Ate{i,l} \; T_{l}(s) \cos\a_i + \sum_i \sum_{l=0}^L \Ato{i,l} \; T_{l}(s) \sin\a_i \\
	               A_\z & = & \sum_i \sum_{l=0}^L \Aze{i,l} \; T_{l}(s) \cos\a_i + \sum_i \sum_{l=0}^L \Azo{i,l} \; T_{l}(s) \sin\a_i ,
	\ee
	where $T_l(s)$ are the Chebyshev polynomials and $\a_i \equiv m_i \t - n_i \z$.
	The following internal arrays are declared in \link{preset}
	\verb{dAte(0,i)%s(l){$\equiv \Ate{i,l}$
	\verb{dAze(0,i)%s(l){$\equiv \Aze{i,l}$
	\verb{dAto(0,i)%s(l){$\equiv \Ato{i,l}$
	\verb{dAzo(0,i)%s(l){$\equiv \Azo{i,l}$ */
	public static class VectorPotential {
		public double[][] Ate;
		public double[][] Aze;
		public double[][] Ato;
		public double[][] Azo;
	} // end of VectorPotential

	/** Initialize complex datatypes. */
	public SpecOutput() {
		input = new Input();
		output = new Output();
		vector_potential = new VectorPotential();
	}

	/**
	 * Initalize complex datatypes and load SpecOutput contents from a HDF5 file identified by {@code filename}.
	 * @param filename path to the HDF5 file to load
	 */
	public SpecOutput(String filename) {
		this();
		try {
			NetcdfFile file = NetcdfFile.open(filename);
			loadFrom(file);
			file.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Initalize complex datatypes and load SpecOutput contents from an already-open NetCDF file identified by {@code file}.
	 * @param file open file to load the data from
	 */
	public SpecOutput(NetcdfFile file) {
		this();
		try {
			loadFrom(file);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/** SPEC version */
	public double version;
	/** input data for this SPEC run */
	public Input input;
	/** output data; content of the previous HDF5 file */
	public Output output;
	/** The covariant components of the vector potential are written as
	\be            A_\t & = & \sum_i \sum_{l=0}^L \Ate{i,l} \; T_{l}(s) \cos\a_i + \sum_i \sum_{l=0}^L \Ato{i,l} \; T_{l}(s) \sin\a_i \\
	               A_\z & = & \sum_i \sum_{l=0}^L \Aze{i,l} \; T_{l}(s) \cos\a_i + \sum_i \sum_{l=0}^L \Azo{i,l} \; T_{l}(s) \sin\a_i ,
	\ee
	where $T_l(s)$ are the Chebyshev polynomials and $\a_i \equiv m_i \t - n_i \z$.
	The following internal arrays are declared in \link{preset}
	\verb{dAte(0,i)%s(l){$\equiv \Ate{i,l}$
	\verb{dAze(0,i)%s(l){$\equiv \Aze{i,l}$
	\verb{dAto(0,i)%s(l){$\equiv \Ato{i,l}$
	\verb{dAzo(0,i)%s(l){$\equiv \Azo{i,l}$ */
	public VectorPotential vector_potential;

	/**
	 * Load SpecOutput contents from an already-open NetCDF file identified by {@code file}.
	 * @param file open file to load the data from
	 * @return initialized SpecOutput object
	 */
	public SpecOutput loadFrom(NetcdfFile file) throws IOException {
		input.physics.Igeometry = file.findVariable("/input/physics/Igeometry").readScalarInt();
		input.physics.Istellsym = file.findVariable("/input/physics/Istellsym").readScalarInt();
		input.physics.Lfreebound = file.findVariable("/input/physics/Lfreebound").readScalarInt();
		input.physics.phiedge = file.findVariable("/input/physics/phiedge").readScalarDouble();
		input.physics.curtor = file.findVariable("/input/physics/curtor").readScalarDouble();
		input.physics.curpol = file.findVariable("/input/physics/curpol").readScalarDouble();
		input.physics.gamma = file.findVariable("/input/physics/gamma").readScalarDouble();
		input.physics.Nfp = file.findVariable("/input/physics/Nfp").readScalarInt();
		input.physics.Nvol = file.findVariable("/input/physics/Nvol").readScalarInt();
		input.physics.Mpol = file.findVariable("/input/physics/Mpol").readScalarInt();
		input.physics.Ntor = file.findVariable("/input/physics/Ntor").readScalarInt();
		input.physics.Lrad = (int[])file.findVariable("/input/physics/Lrad").read().get1DJavaArray(DataType.INT);
		input.physics.Lconstraint = file.findVariable("/input/physics/Lconstraint").readScalarInt();
		input.physics.tflux = (double[])file.findVariable("/input/physics/tflux").read().get1DJavaArray(DataType.DOUBLE);
		input.physics.pflux = (double[])file.findVariable("/input/physics/pflux").read().get1DJavaArray(DataType.DOUBLE);
		input.physics.helicity = (double[])file.findVariable("/input/physics/helicity").read().get1DJavaArray(DataType.DOUBLE);
		input.physics.pscale = file.findVariable("/input/physics/pscale").readScalarDouble();
		input.physics.pressure = (double[])file.findVariable("/input/physics/pressure").read().get1DJavaArray(DataType.DOUBLE);
		input.physics.Ladiabatic = file.findVariable("/input/physics/Ladiabatic").readScalarInt();
		input.physics.adiabatic = (double[])file.findVariable("/input/physics/adiabatic").read().get1DJavaArray(DataType.DOUBLE);
		input.physics.mu = (double[])file.findVariable("/input/physics/mu").read().get1DJavaArray(DataType.DOUBLE);
		input.physics.pl = (int[])file.findVariable("/input/physics/pl").read().get1DJavaArray(DataType.INT);
		input.physics.ql = (int[])file.findVariable("/input/physics/ql").read().get1DJavaArray(DataType.INT);
		input.physics.pr = (int[])file.findVariable("/input/physics/pr").read().get1DJavaArray(DataType.INT);
		input.physics.qr = (int[])file.findVariable("/input/physics/qr").read().get1DJavaArray(DataType.INT);
		input.physics.iota = (double[])file.findVariable("/input/physics/iota").read().get1DJavaArray(DataType.DOUBLE);
		input.physics.lp = (int[])file.findVariable("/input/physics/lp").read().get1DJavaArray(DataType.INT);
		input.physics.lq = (int[])file.findVariable("/input/physics/lq").read().get1DJavaArray(DataType.INT);
		input.physics.rp = (int[])file.findVariable("/input/physics/rp").read().get1DJavaArray(DataType.INT);
		input.physics.rq = (int[])file.findVariable("/input/physics/rq").read().get1DJavaArray(DataType.INT);
		input.physics.oita = (double[])file.findVariable("/input/physics/oita").read().get1DJavaArray(DataType.DOUBLE);
		input.physics.mupftol = file.findVariable("/input/physics/mupftol").readScalarDouble();
		input.physics.mupfits = file.findVariable("/input/physics/mupfits").readScalarInt();
		input.physics.Rac = (double[])file.findVariable("/input/physics/Rac").read().get1DJavaArray(DataType.DOUBLE);
		input.physics.Zas = (double[])file.findVariable("/input/physics/Zas").read().get1DJavaArray(DataType.DOUBLE);
		input.physics.Ras = (double[])file.findVariable("/input/physics/Ras").read().get1DJavaArray(DataType.DOUBLE);
		input.physics.Zac = (double[])file.findVariable("/input/physics/Zac").read().get1DJavaArray(DataType.DOUBLE);
		input.physics.Rbc = (double[][])file.findVariable("/input/physics/Rbc").read().copyToNDJavaArray();
		input.physics.Zbs = (double[][])file.findVariable("/input/physics/Zbs").read().copyToNDJavaArray();
		input.physics.Rbs = (double[][])file.findVariable("/input/physics/Rbs").read().copyToNDJavaArray();
		input.physics.Zbc = (double[][])file.findVariable("/input/physics/Zbc").read().copyToNDJavaArray();
		input.physics.Rwc = (double[][])file.findVariable("/input/physics/Rwc").read().copyToNDJavaArray();
		input.physics.Zws = (double[][])file.findVariable("/input/physics/Zws").read().copyToNDJavaArray();
		input.physics.Rws = (double[][])file.findVariable("/input/physics/Rws").read().copyToNDJavaArray();
		input.physics.Zwc = (double[][])file.findVariable("/input/physics/Zwc").read().copyToNDJavaArray();
		input.physics.Vns = (double[][])file.findVariable("/input/physics/Vns").read().copyToNDJavaArray();
		input.physics.Bns = (double[][])file.findVariable("/input/physics/Bns").read().copyToNDJavaArray();
		input.physics.Vnc = (double[][])file.findVariable("/input/physics/Vnc").read().copyToNDJavaArray();
		input.physics.Bnc = (double[][])file.findVariable("/input/physics/Bnc").read().copyToNDJavaArray();
		input.numerics.Linitialize = file.findVariable("/input/numerics/Linitialize").readScalarInt();
		input.numerics.Lzerovac = file.findVariable("/input/numerics/Lzerovac").readScalarInt();
		input.numerics.Ndiscrete = file.findVariable("/input/numerics/Ndiscrete").readScalarInt();
		input.numerics.Nquad = file.findVariable("/input/numerics/Nquad").readScalarInt();
		input.numerics.iMpol = file.findVariable("/input/numerics/iMpol").readScalarInt();
		input.numerics.iNtor = file.findVariable("/input/numerics/iNtor").readScalarInt();
		input.numerics.Lsparse = file.findVariable("/input/numerics/Lsparse").readScalarInt();
		input.numerics.Lsvdiota = file.findVariable("/input/numerics/Lsvdiota").readScalarInt();
		input.numerics.imethod = file.findVariable("/input/numerics/imethod").readScalarInt();
		input.numerics.iorder = file.findVariable("/input/numerics/iorder").readScalarInt();
		input.numerics.iprecon = file.findVariable("/input/numerics/iprecon").readScalarInt();
		input.numerics.iotatol = file.findVariable("/input/numerics/iotatol").readScalarDouble();
		input.numerics.Lextrap = file.findVariable("/input/numerics/Lextrap").readScalarInt();
		input.numerics.Mregular = file.findVariable("/input/numerics/Mregular").readScalarInt();
		input.local.LBeltrami = file.findVariable("/input/local/LBeltrami").readScalarInt();
		input.local.Linitgues = file.findVariable("/input/local/Linitgues").readScalarInt();
		input.local.maxrndgues = file.findVariable("/input/local/maxrndgues").readScalarDouble();
		input.local.Lposdef = file.findVariable("/input/local/Lposdef").readScalarInt();
		input.global.Lfindzero = file.findVariable("/input/global/Lfindzero").readScalarInt();
		input.global.escale = file.findVariable("/input/global/escale").readScalarDouble();
		input.global.opsilon = file.findVariable("/input/global/opsilon").readScalarDouble();
		input.global.pcondense = file.findVariable("/input/global/pcondense").readScalarDouble();
		input.global.epsilon = file.findVariable("/input/global/epsilon").readScalarDouble();
		input.global.wpoloidal = file.findVariable("/input/global/wpoloidal").readScalarDouble();
		input.global.upsilon = file.findVariable("/input/global/upsilon").readScalarDouble();
		input.global.forcetol = file.findVariable("/input/global/forcetol").readScalarDouble();
		input.global.c05xmax = file.findVariable("/input/global/c05xmax").readScalarDouble();
		input.global.c05xtol = file.findVariable("/input/global/c05xtol").readScalarDouble();
		input.global.c05factor = file.findVariable("/input/global/c05factor").readScalarDouble();
		input.global.LreadGF = (file.findVariable("/input/global/LreadGF").readScalarInt() > 0 ? true : false);
		input.global.mfreeits = file.findVariable("/input/global/mfreeits").readScalarInt();
		input.global.bnstol = file.findVariable("/input/global/bnstol").readScalarDouble();
		input.global.bnsblend = file.findVariable("/input/global/bnsblend").readScalarDouble();
		input.global.gBntol = file.findVariable("/input/global/gBntol").readScalarDouble();
		input.global.gBnbld = file.findVariable("/input/global/gBnbld").readScalarDouble();
		input.global.vcasingeps = file.findVariable("/input/global/vcasingeps").readScalarDouble();
		input.global.vcasingtol = file.findVariable("/input/global/vcasingtol").readScalarDouble();
		input.global.vcasingits = file.findVariable("/input/global/vcasingits").readScalarInt();
		input.global.vcasingper = file.findVariable("/input/global/vcasingper").readScalarInt();
		input.global.mcasingcal = file.findVariable("/input/global/mcasingcal").readScalarInt();
		input.diagnostics.odetol = file.findVariable("/input/diagnostics/odetol").readScalarDouble();
		input.diagnostics.absreq = file.findVariable("/input/diagnostics/absreq").readScalarDouble();
		input.diagnostics.relreq = file.findVariable("/input/diagnostics/relreq").readScalarDouble();
		input.diagnostics.absacc = file.findVariable("/input/diagnostics/absacc").readScalarDouble();
		input.diagnostics.epsr = file.findVariable("/input/diagnostics/epsr").readScalarDouble();
		input.diagnostics.nPpts = file.findVariable("/input/diagnostics/nPpts").readScalarInt();
		input.diagnostics.nPtrj = (int[])file.findVariable("/input/diagnostics/nPtrj").read().get1DJavaArray(DataType.INT);
		input.diagnostics.LHevalues = (file.findVariable("/input/diagnostics/LHevalues").readScalarInt() > 0 ? true : false);
		input.diagnostics.LHevectors = (file.findVariable("/input/diagnostics/LHevectors").readScalarInt() > 0 ? true : false);
		input.diagnostics.LHmatrix = (file.findVariable("/input/diagnostics/LHmatrix").readScalarInt() > 0 ? true : false);
		input.diagnostics.Lperturbed = file.findVariable("/input/diagnostics/Lperturbed").readScalarInt();
		input.diagnostics.dpp = file.findVariable("/input/diagnostics/dpp").readScalarInt();
		input.diagnostics.dqq = file.findVariable("/input/diagnostics/dqq").readScalarInt();
		input.diagnostics.Lcheck = file.findVariable("/input/diagnostics/Lcheck").readScalarInt();
		input.diagnostics.Ltiming = (file.findVariable("/input/diagnostics/Ltiming").readScalarInt() > 0 ? true : false);
		input.diagnostics.fudge = file.findVariable("/input/diagnostics/fudge").readScalarDouble();
		input.diagnostics.scaling = file.findVariable("/input/diagnostics/scaling").readScalarDouble();
		output.Vns = (double[])file.findVariable("/output/Vns").read().get1DJavaArray(DataType.DOUBLE);
		output.Bns = (double[])file.findVariable("/output/Bns").read().get1DJavaArray(DataType.DOUBLE);
		output.Vnc = (double[])file.findVariable("/output/Vnc").read().get1DJavaArray(DataType.DOUBLE);
		output.Bnc = (double[])file.findVariable("/output/Bnc").read().get1DJavaArray(DataType.DOUBLE);
		output.mn = file.findVariable("/output/mn").readScalarInt();
		output.im = (int[])file.findVariable("/output/im").read().get1DJavaArray(DataType.INT);
		output.in = (int[])file.findVariable("/output/in").read().get1DJavaArray(DataType.INT);
		output.Mvol = file.findVariable("/output/Mvol").readScalarInt();
		output.Rbc = (double[][])file.findVariable("/output/Rbc").read().copyToNDJavaArray();
		output.Zbs = (double[][])file.findVariable("/output/Zbs").read().copyToNDJavaArray();
		output.Rbs = (double[][])file.findVariable("/output/Rbs").read().copyToNDJavaArray();
		output.Zbc = (double[][])file.findVariable("/output/Zbc").read().copyToNDJavaArray();
		output.ForceErr = file.findVariable("/output/ForceErr").readScalarDouble();
		output.adiabatic = (double[])file.findVariable("/output/adiabatic").read().get1DJavaArray(DataType.DOUBLE);
		output.helicity = (double[])file.findVariable("/output/helicity").read().get1DJavaArray(DataType.DOUBLE);
		output.mu = (double[])file.findVariable("/output/mu").read().get1DJavaArray(DataType.DOUBLE);
		output.tflux = (double[])file.findVariable("/output/tflux").read().get1DJavaArray(DataType.DOUBLE);
		output.pflux = (double[])file.findVariable("/output/pflux").read().get1DJavaArray(DataType.DOUBLE);
		output.volume = file.findVariable("/output/volume").readScalarDouble();
		output.Mrad = file.findVariable("/output/Mrad").readScalarInt();
		output.TT = (double[][][])file.findVariable("/output/TT").read().copyToNDJavaArray();
		output.Btemn = (double[][][])file.findVariable("/output/Btemn").read().copyToNDJavaArray();
		output.Bzemn = (double[][][])file.findVariable("/output/Bzemn").read().copyToNDJavaArray();
		output.Btomn = (double[][][])file.findVariable("/output/Btomn").read().copyToNDJavaArray();
		output.Bzomn = (double[][][])file.findVariable("/output/Bzomn").read().copyToNDJavaArray();
		output.lmns = file.findVariable("/output/lmns").readScalarDouble();
		vector_potential.Ate = (double[][])file.findVariable("/vector_potential/Ate").read().copyToNDJavaArray();
		vector_potential.Aze = (double[][])file.findVariable("/vector_potential/Aze").read().copyToNDJavaArray();
		vector_potential.Ato = (double[][])file.findVariable("/vector_potential/Ato").read().copyToNDJavaArray();
		vector_potential.Azo = (double[][])file.findVariable("/vector_potential/Azo").read().copyToNDJavaArray();
		version = file.findVariable("/version").readScalarDouble();
		return this;
	}

    public static void main(String[] args) {
        SpecOutput s = new SpecOutput("/home/jonathan/Uni/04_PhD/00_programs/SPEC/SPEC/InputFiles/TestCases/G3V02L1Fi.001.h5");
        System.out.printf(Locale.ENGLISH, "SPEC version: %.2f\n", s.version);
    }
} // end of SpecOutput
