#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 20:22:55 2019

@author: Jonathan Schilling (jonathan.schilling@ipp.mpg.de)
"""

#%% define the output file structure of the SPEC HDF5 file
from Hdf5File import Hdf5File, Group, Dataset

# Define SPEC output file format
s = Hdf5File("SpecOutput")
s.rootGroup.setDescription("Output file of the SPEC MHD equilibrium code")

d_version = Dataset(s, "version", "double")
d_version.setDescription("SPEC version")
d_version.setStorageLocation({"w": "global"})
d_version.setReadWriteRoutines({'w': "mirror_input_to_outfile"})

g_input = Group(s, "input")
g_input.setDescription("input data for this SPEC run")

g_input_physics = Group(g_input, "physics")
g_input_physics.setDescription(r"The namelist \verb+physicslist+ controls the geometry, profiles, and numerical resolution.")

d_input_physics_Igeometry = Dataset(g_input_physics, "Igeometry", "int")
d_input_physics_Igeometry.setDefaultValue(3)
d_input_physics_Igeometry.setDescription(r"""selects Cartesian, cylindrical or toroidal geometry
\bi
\item[i.] \inputvar{Igeometry = 1} : Cartesian; geometry determined by $R$
\item[i.] \inputvar{Igeometry = 2} : cylindrical; geometry determined by $R$
\item[i.] \inputvar{Igeometry = 3} : toroidal; geometry determined by $R$ {\em and} $Z$
\ei""")

d_input_physics_Istellsym = Dataset(g_input_physics, "Istellsym", "int")
d_input_physics_Istellsym.setDefaultValue(1)
d_input_physics_Istellsym.setDescription(r"stellarator symmetry is enforced if \inputvar{Istellsym.eq.1}")

d_input_physics_Lfreebound = Dataset(g_input_physics, "Lfreebound", "int")
d_input_physics_Lfreebound.setDefaultValue(0)
d_input_physics_Lfreebound.setDescription("compute vacuum field surrounding plasma")

d_input_physics_phiedge = Dataset(g_input_physics, "phiedge", "double")
d_input_physics_phiedge.setDefaultValue(1.0)
d_input_physics_phiedge.setDescription("total enclosed toroidal magnetic flux")

d_input_physics_curtor = Dataset(g_input_physics, "curtor", "double")
d_input_physics_curtor.setDefaultValue(0.0)
d_input_physics_curtor.setDescription("total enclosed (toroidal) plasma current")

d_input_physics_curpol = Dataset(g_input_physics, "curpol", "double")
d_input_physics_curpol.setDefaultValue(0.0)
d_input_physics_curpol.setDescription("total enclosed (poloidal) linking current")

d_input_physics_gamma = Dataset(g_input_physics, "gamma", "double")
d_input_physics_gamma.setDefaultValue(0.0)
d_input_physics_gamma.setDescription(r"adiabatic index; cannot set $|\gamma| = 1$")

d_input_physics_Nfp = Dataset(g_input_physics, "Nfp", "int")
d_input_physics_Nfp.setDefaultValue(1)
d_input_physics_Nfp.setDescription(r"""field periodicity
\bi
\item[i.] all Fourier representations are of the form $\cos(m\t-nN\z)$, $\sin(m\t-nN\z)$,where $N\equiv$\inputvar{Nfp}
\item[i.] constraint : \inputvar{Nfp.ge.1}
\ei""")

d_input_physics_Nvol = Dataset(g_input_physics, "Nvol", "int")
d_input_physics_Nvol.setDefaultValue(1)
d_input_physics_Nvol.setDescription(r"""number of volumes
\bi
\item[i.] each volume ${\cal V}_l$ is bounded by the ${\cal I}_{l-1}$ and ${\cal I}_{l}$ interfaces
\item[i.] note that in cylindrical or toroidal geometry, ${\cal I}_{0}$ is the degenerate coordinate axis
\item[i.] constraint : \inputvar{Nvol.le.MNvol}
\ei""")

d_input_physics_mpol = Dataset(g_input_physics, "Mpol", "int")
d_input_physics_mpol.setDescription("number of poloidal Fourier harmonics")
d_input_physics_mpol.setStorageLocation({"w": "global"})
d_input_physics_mpol.setReadWriteRoutines({'w': "mirror_input_to_outfile"})

d_input_physics_ntor = Dataset(g_input_physics, "Ntor", "int")
d_input_physics_ntor.setDescription("number of toroidal Fourier harmonics")
d_input_physics_ntor.setStorageLocation({"w": "global"})
d_input_physics_ntor.setReadWriteRoutines({'w': "mirror_input_to_outfile"})

d_input_physics_Lrad = Dataset(g_input_physics, "Lrad", "int", 1)
d_input_physics_Lrad.setDefaultValue([4])
d_input_physics_Lrad.setDescription(r"""Chebyshev resolution in each volume
\bi
\item[i.] constraint : \inputvar{Lrad(1:Mvol)}.ge.2
\ei""")
d_input_physics_Lrad.setStorageLocation({"w": "global"})
d_input_physics_Lrad.setReadWriteRoutines({'w': "mirror_input_to_outfile"})

d_input_physics_Lconstraint = Dataset(g_input_physics, "Lconstraint", "int")
d_input_physics_Lconstraint.setDefaultValue(-1)
d_input_physics_Lconstraint.setDescription(r"""selects constraints; primarily used in \link{ma02aa} and \link{mp00ac}.
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
\ei""")

d_input_physics_tflux = Dataset(g_input_physics, "tflux", "double", 1)
d_input_physics_tflux.setDescription(r"""toroidal flux, $\psi_t$, enclosed by each interface
\bi
\item[i.] For each of the plasma volumes, this is a constraint: \inputvar{tflux} is \emph{not} varied
\item[i.] For the vacuum region (only if \inputvar{Lfreebound = 1}), \inputvar{tflux} may be allowed to vary to match constraints
\item[i.] Note that \inputvar{tflux} will be normalized so that \inputvar{tflux(Nvol) = 1.0},
          so that \inputvar{tflux} is arbitrary up to a scale factor
\item[i.] see also \inputvar{phiedge}
\ei""")

d_input_physics_pflux = Dataset(g_input_physics, "pflux", "double", 1)
d_input_physics_pflux.setDescription(r"poloidal flux, $\psi_p$, enclosed by each interface")

d_input_physics_helicity = Dataset(g_input_physics, "helicity", "double", 1)
d_input_physics_helicity.setDescription(r"""helicity, ${\cal K}$, in each volume, ${\cal V}_i$
\bi
\item[i.] on exit, \inputvar{helicity} is set to the computed values of ${\cal K} \equiv \int {\bf A}\cdot{\bf B}\;dv$
\ei""")

d_input_physics_pscale = Dataset(g_input_physics, "pscale", "double")
d_input_physics_pscale.setDefaultValue(0.0)
d_input_physics_pscale.setDescription(r"""pressure scale factor
\bi
\item[i.] the initial pressure profile is given by \inputvar{pscale} $*$ \inputvar{press}
\ei""")

d_input_physics_pressure = Dataset(g_input_physics, "pressure", "double", 1)
d_input_physics_pressure.setDescription(r"""pressure in each volume
\bi
\item[i.] the pressure is {\em not} held constant, but $p_l V_l^\gamma = P_l$ {\em is} held constant,
          where $P_l$ is determined by the initial pressures and the initial volumes, $V_l$
\item[i.] (Note that if \inputvar{gamma = 0.0}, then $p_l \equiv P_l$.)
\item[i.] on output, the pressure is given by $p_l = P_l/V_l^\gamma$, where $V_l$ is the final volume
\item[i.] \inputvar{pressure} is only used in calculation of interface force-balance
\ei""")

d_input_physics_Ladiabatic = Dataset(g_input_physics, "Ladiabatic", "int")
d_input_physics_Ladiabatic.setDefaultValue(0)
d_input_physics_Ladiabatic.setDescription(r"""logical flag
\bi
\item[i.] if \inputvar{Ladiabatic = 0}, the adiabatic constants are determined by the initial pressure and volume
\item[i.] if \inputvar{Ladiabatic = 1}, the adiabatic constants are determined by the given input \inputvar{adiabatic}
\ei""")

d_input_physics_adiabatic = Dataset(g_input_physics, "adiabatic", "double", 1)
d_input_physics_adiabatic.setDescription(r"""adiabatic constants in each volume
\bi
\item[i.] the pressure is {\em not} held constant, but $p_l V_l^\gamma = P_l \equiv$\inputvar{adiabatic} is constant,
\item[i.] note that if \inputvar{gamma = 0.0}, then \inputvar{pressure = adiabatic}
\item[i.] \inputvar{pressure} is only used in calculation of interface force-balance
\ei""")

d_input_physics_mu = Dataset(g_input_physics, "mu", "double", 1)
d_input_physics_mu.setDescription(r"helicity-multiplier, $\mu$, in each volume")

d_input_physics_pl = Dataset(g_input_physics, "pl", "int", 1)
d_input_physics_pl.setDefaultValue(0)
d_input_physics_pl.setDescription(r"""part of noble irrational defining iota
\bi
\item[i.] ``inside'' interface rotational-transform is $\iotabar = (p_l+\gamma p_r)/(q_l+\gamma q_r)$,
          where $\gamma$ is the golden mean, $\gamma = (1 + \sqrt 5 ) / 2 $
\item[i.] if both $q_l = 0$ {\em and} $q_r = 0$, then the (inside) interface rotational-transform is defined by \inputvar{iota}
\ei""")

d_input_physics_ql = Dataset(g_input_physics, "ql", "int", 1)
d_input_physics_ql.setDefaultValue(0)
d_input_physics_ql.setDescription(r"""part of noble irrational defining iota
\bi
\item[i.] ``inside'' interface rotational-transform is $\iotabar = (p_l+\gamma p_r)/(q_l+\gamma q_r)$,
          where $\gamma$ is the golden mean, $\gamma = (1 + \sqrt 5 ) / 2 $
\item[i.] if both $q_l = 0$ {\em and} $q_r = 0$, then the (inside) interface rotational-transform is defined by \inputvar{iota}
\ei""")

d_input_physics_pr = Dataset(g_input_physics, "pr", "int", 1)
d_input_physics_pr.setDefaultValue(0)
d_input_physics_pr.setDescription(r"""part of noble irrational defining iota
\bi
\item[i.] ``inside'' interface rotational-transform is $\iotabar = (p_l+\gamma p_r)/(q_l+\gamma q_r)$,
          where $\gamma$ is the golden mean, $\gamma = (1 + \sqrt 5 ) / 2 $
\item[i.] if both $q_l = 0$ {\em and} $q_r = 0$, then the (inside) interface rotational-transform is defined by \inputvar{iota}
\ei""")

d_input_physics_qr = Dataset(g_input_physics, "qr", "int", 1)
d_input_physics_qr.setDefaultValue(0)
d_input_physics_qr.setDescription(r"""part of noble irrational defining iota
\bi
\item[i.] ``inside'' interface rotational-transform is $\iotabar = (p_l+\gamma p_r)/(q_l+\gamma q_r)$,
          where $\gamma$ is the golden mean, $\gamma = (1 + \sqrt 5 ) / 2 $
\item[i.] if both $q_l = 0$ {\em and} $q_r = 0$, then the (inside) interface rotational-transform is defined by \inputvar{iota}
\ei""")

d_input_physics_iota = Dataset(g_input_physics, "iota", "double", 1)
d_input_physics_iota.setDescription(r"""rotational-transform, $\iotabar$, on inner side of each interface
\bi
\item[i.] only relevant if \inputvar{ql}=0 and \inputvar{qr}=0
\ei""")

d_input_physics_lp = Dataset(g_input_physics, "lp", "int", 1)
d_input_physics_lp.setDefaultValue(0)
d_input_physics_lp.setDescription(r"""part of noble irrational defining oita
\bi
\item "outer" interface rotational-transform is $\iotabar = (p_l+\gamma p_r)/(q_l+\gamma q_r)$,
      where $\gamma$ is the golden mean, $\gamma = (1 + \sqrt 5 ) / 2 $
\item if both $q_l = 0$ {\em and} $q_r = 0$, then the (outer) interface rotational-transform is defined by \inputvar{oita}
\ei""")

d_input_physics_lq = Dataset(g_input_physics, "lq", "int", 1)
d_input_physics_lq.setDefaultValue(0)
d_input_physics_lq.setDescription(r"""part of noble irrational defining oita
\bi
\item "outer" interface rotational-transform is $\iotabar = (p_l+\gamma p_r)/(q_l+\gamma q_r)$,
      where $\gamma$ is the golden mean, $\gamma = (1 + \sqrt 5 ) / 2 $
\item if both $q_l = 0$ {\em and} $q_r = 0$, then the (outer) interface rotational-transform is defined by \inputvar{oita}
\ei""")

d_input_physics_rp = Dataset(g_input_physics, "rp", "int", 1)
d_input_physics_rp.setDefaultValue(0)
d_input_physics_rp.setDescription(r"""part of noble irrational defining oita
\bi
\item "outer" interface rotational-transform is $\iotabar = (p_l+\gamma p_r)/(q_l+\gamma q_r)$,
      where $\gamma$ is the golden mean, $\gamma = (1 + \sqrt 5 ) / 2 $
\item if both $q_l = 0$ {\em and} $q_r = 0$, then the (outer) interface rotational-transform is defined by \inputvar{oita}
\ei""")

d_input_physics_rq = Dataset(g_input_physics, "rq", "int", 1)
d_input_physics_rq.setDefaultValue(0)
d_input_physics_rq.setDescription(r"""part of noble irrational defining oita
\bi
\item "outer" interface rotational-transform is $\iotabar = (p_l+\gamma p_r)/(q_l+\gamma q_r)$,
      where $\gamma$ is the golden mean, $\gamma = (1 + \sqrt 5 ) / 2 $
\item if both $q_l = 0$ {\em and} $q_r = 0$, then the (outer) interface rotational-transform is defined by \inputvar{oita}
\ei""")

d_input_physics_oita = Dataset(g_input_physics, "oita", "double", 1)
d_input_physics_oita.setDescription(r"""rotational-transform, $\iotabar$, on outer side of each interface
\bi
\item only relevant if illogical input for \inputvar{ql} and \inputvar{qr} are provided
\ei""")

d_input_physics_mupftol = Dataset(g_input_physics, "mupftol", "double")
d_input_physics_mupftol.setDefaultValue(1.0e-16)
d_input_physics_mupftol.setDescription(r"""accuracy to which $\mu$ and $\Delta\psi_p$ are required
\bi
\item only relevant if constraints on transform, enclosed currents etc. are to be satisfied iteratively, see \inputvar{Lconstraint}
\ei""")

d_input_physics_mupfits = Dataset(g_input_physics, "mupfits", "int")
d_input_physics_mupfits.setDefaultValue(8)
d_input_physics_mupfits.setDescription(r"""an upper limit on the transform/helicity constraint iterations
\item only relevant if constraints on transform, enclosed currents etc. are to be satisfied iteratively, see \inputvar{Lconstraint}
\bi
\item constraint: \inputvar{mupfits > 0}
\ei""")

# commented out until they are (always or never) written by SPEC
#d_input_physics_rpol = Dataset(g_input_physics, "rpol", "double")
#d_input_physics_rpol.setDefaultValue(1.0)
#d_input_physics_rpol.setDescription(r"""poloidal extent of slab (effective radius)
#\bi
#\item[i.] only relevant if \inputvar{Igeometry} $=1$
#\item[i.] poloidal size is $L = 2\pi*$\inputvar{rpol}
#\ei""")
#
#d_input_physics_rtor = Dataset(g_input_physics, "rtor", "double")
#d_input_physics_rtor.setDefaultValue(1.0)
#d_input_physics_rtor.setDescription(r"""toroidal extent of slab (effective radius)
#\bi
#\item[i.] only relevant if \inputvar{Igeometry} $=1$
#\item[i.] toroidal size is $L = 2\pi*$\inputvar{rtor}
#\ei""")

d_input_physics_Rac = Dataset(g_input_physics, "Rac", "double", 1)
d_input_physics_Rac.setDescription("stellarator symmetric coordinate axis R cosine Fourier coefficients")

d_input_physics_Zas = Dataset(g_input_physics, "Zas", "double", 1)
d_input_physics_Zas.setDescription("stellarator symmetric coordinate axis Z sine Fourier coefficients")

d_input_physics_Ras = Dataset(g_input_physics, "Ras", "double", 1)
d_input_physics_Ras.setDescription("non-stellarator symmetric coordinate axis R sine Fourier coefficients")

d_input_physics_Zac = Dataset(g_input_physics, "Zac", "double", 1)
d_input_physics_Zac.setDescription("non-stellarator symmetric coordinate axis Z cosine Fourier coefficients")

d_input_physics_Rbc = Dataset(g_input_physics, "Rbc", "double", 2)
d_input_physics_Rbc.setDescription("stellarator symmetric boundary R cosine Fourier coefficients")

d_input_physics_Zbs = Dataset(g_input_physics, "Zbs", "double", 2)
d_input_physics_Zbs.setDescription("stellarator symmetric boundary Z sine Fourier coefficients")

d_input_physics_Rbs = Dataset(g_input_physics, "Rbs", "double", 2)
d_input_physics_Rbs.setDescription("non-stellarator symmetric boundary R sine Fourier coefficients")

d_input_physics_Zbc = Dataset(g_input_physics, "Zbc", "double", 2)
d_input_physics_Zbc.setDescription("non-stellarator symmetric boundary Z cosine Fourier coefficients")

d_input_physics_Rwc = Dataset(g_input_physics, "Rwc", "double", 2)
d_input_physics_Rwc.setDescription("stellarator symmetric boundary R cosine Fourier coefficients of wall")

d_input_physics_Zws = Dataset(g_input_physics, "Zws", "double", 2)
d_input_physics_Zws.setDescription("stellarator symmetric boundary Z sine Fourier coefficients of wall")

d_input_physics_Rws = Dataset(g_input_physics, "Rws", "double", 2)
d_input_physics_Rws.setDescription("non-stellarator symmetric boundary R sine Fourier coefficients of wall")

d_input_physics_Zwc = Dataset(g_input_physics, "Zwc", "double", 2)
d_input_physics_Zwc.setDescription("non-stellarator symmetric boundary Z cosine Fourier coefficients of wall")

d_input_physics_Vns = Dataset(g_input_physics, "Vns", "double", 2)
d_input_physics_Vns.setDescription("stellarator symmetric normal field at boundary; vacuum component")

d_input_physics_Bns = Dataset(g_input_physics, "Bns", "double", 2)
d_input_physics_Bns.setDescription("stellarator symmetric normal field at boundary; plasma component")

d_input_physics_Vnc = Dataset(g_input_physics, "Vnc", "double", 2)
d_input_physics_Vnc.setDescription("non-stellarator symmetric normal field at boundary; vacuum component")

d_input_physics_Bnc = Dataset(g_input_physics, "Bnc", "double", 2)
d_input_physics_Bnc.setDescription("non-stellarator symmetric normal field at boundary; plasma component")

g_input_numerics = Group(g_input, "numerics")
g_input_numerics.setDescription(r"The namelist \verb+numericlist+ controls internal resolution parameters that the user rarely needs to consider.")

d_input_numerics_Linitialize = Dataset(g_input_numerics, "Linitialize", "int")
d_input_numerics_Linitialize.setDefaultValue(0)
d_input_numerics_Linitialize.setDescription(r"""to initialize geometry using a regularization / extrapolation method
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
\ei""")

d_input_numerics_Lzerovac = Dataset(g_input_numerics, "Lzerovac", "int")
d_input_numerics_Lzerovac.setDefaultValue(0)
d_input_numerics_Lzerovac.setDescription(r"""to adjust vacuum field to cancel plasma field on computational boundary
\bi
\item only relevant if \inputvar{Lfreebound = 1},
\ei""")

d_input_numerics_Ndiscrete = Dataset(g_input_numerics, "Ndiscrete", "int")
d_input_numerics_Ndiscrete.setDefaultValue(2)
d_input_numerics_Ndiscrete.setDescription(r"""\bi
\item resolution of the real space grid on which fast Fourier transforms are performed is given by \inputvar{Ndiscrete*Mpol*4}
\item constraint \inputvar{Ndiscrete>0}
\ei""")

d_input_numerics_Nquad = Dataset(g_input_numerics, "Nquad", "int")
d_input_numerics_Nquad.setDefaultValue(-1)
d_input_numerics_Nquad.setDescription(r"""the resolution of the Gaussian quadrature
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
\ei""")
                                      
d_input_numerics_iMpol = Dataset(g_input_numerics, "iMpol", "int")
d_input_numerics_iMpol.setDefaultValue(-4)
d_input_numerics_iMpol.setDescription(r"""Fourier resolution of straight-fieldline angle on interfaces
\bi
\item the rotational-transform on the interfaces is determined by a transformation to the straight-fieldline angle,
      with poloidal resolution given by \inputvar{iMpol}
\item if \inputvar{iMpol.le.0}, then \inputvar{iMpol = Mpol - iMpol}
\ei""")

d_input_numerics_iNtor = Dataset(g_input_numerics, "iNtor", "int")
d_input_numerics_iNtor.setDefaultValue(-4)
d_input_numerics_iNtor.setDescription(r"""Fourier resolution of straight-fieldline angle on interfaces
\bi
\item the rotational-transform on the interfaces is determined by a transformation to the straight-fieldline angle,
      with toroidal resolution given by \inputvar{iNtor}
\item if \inputvar{iNtor.le.0}, then \inputvar{iNtor = Ntor - iNtor}
\item if \inputvar{Ntor.eq.0}, then the toroidal resolution of the angle transformation is set \inputvar{lNtor = 0}.
\ei""")

d_input_numerics_Lsparse = Dataset(g_input_numerics, "Lsparse", "int")
d_input_numerics_Lsparse.setDefaultValue(0)
d_input_numerics_Lsparse.setDescription(r"""controls method used to solve for rotational-transform on interfaces
\bi
\item if \inputvar{Lsparse = 0}, the transformation to the straight-fieldline angle is computed in Fourier space
      using a dense matrix solver, \nag{}{F04AAF}
\item if \inputvar{Lsparse = 1}, the transformation to the straight-fieldline angle is computed in real space
      using a dense matrix solver, \nag{}{F04ATF}
\item if \inputvar{Lsparse = 2}, the transformation to the straight-fieldline angle is computed in real space
      using a sparse matrix solver, \nag{}{F11DEF}
\item if \inputvar{Lsparse = 3}, the different methods for constructing the straight-fieldline angle are compared
\ei""")

d_input_numerics_Lsvdiota = Dataset(g_input_numerics, "Lsvdiota", "int")
d_input_numerics_Lsvdiota.setDefaultValue(0)
d_input_numerics_Lsvdiota.setDescription(r"""controls method used to solve for rotational-transform on interfaces
only relevant if \inputvar{Lsparse = 0}
\bi
\item if \inputvar{Lsvdiota = 0}, use standard linear solver to construct straight fieldline angle transformation
\item if \inputvar{Lsvdiota = 1}, use SVD method to compute rotational-transform
\ei""")

d_input_numerics_imethod = Dataset(g_input_numerics, "imethod", "int")
d_input_numerics_imethod.setDefaultValue(3)
d_input_numerics_imethod.setDescription(r"""controls iterative solution to sparse matrix
arising in real-space transformation to the straight-fieldline angle
only relevant if \inputvar{Lsparse.eq.2}; see \link{tr00ab} for details
\bi
\item if \inputvar{imethod = 1}, the method is \type{RGMRES}
\item if \inputvar{imethod = 2}, the method is \type{CGS}
\item if \inputvar{imethod = 3}, the method is \type{BICGSTAB}
\ei""")

d_input_numerics_iorder = Dataset(g_input_numerics, "iorder", "int")
d_input_numerics_iorder.setDefaultValue(2)
d_input_numerics_iorder.setDescription(r"""controls real-space grid resolution for constructing the straight-fieldline angle
only relevant if \inputvar{Lsparse>0}
determines order of finite-difference approximation to the derivatives
\bi
\item if \inputvar{iorder = 2}, 
\item if \inputvar{iorder = 4}, 
\item if \inputvar{iorder = 6}, 
\ei""")

d_input_numerics_iprecon = Dataset(g_input_numerics, "iprecon", "int")
d_input_numerics_iprecon.setDefaultValue(0)
d_input_numerics_iprecon.setDescription(r"""controls iterative solution to sparse matrix arising in real-space transformation
to the straight-fieldline angle
only relevant if \inputvar{Lsparse.eq.2}; see \link{tr00ab} for details
\bi
\item if \inputvar{iprecon = 0}, the preconditioner is `N'
\item if \inputvar{iprecon = 1}, the preconditioner is `J'
\item if \inputvar{iprecon = 2}, the preconditioner is `S'
\ei""")

d_input_numerics_iotatol = Dataset(g_input_numerics, "iotatol", "double")
d_input_numerics_iotatol.setDefaultValue(-1.0)
d_input_numerics_iotatol.setDescription(r"""tolerance required for iterative construction of straight-fieldline angle
only relevant if \inputvar{Lsparse.ge.2}""")

d_input_numerics_Lextrap = Dataset(g_input_numerics, "Lextrap", "int")
d_input_numerics_Lextrap.setDefaultValue(0)
d_input_numerics_Lextrap.setDescription("geometry of innermost interface is defined by extrapolation")

d_input_numerics_Mregular = Dataset(g_input_numerics, "Mregular", "int")
d_input_numerics_Mregular.setDefaultValue(-1)
d_input_numerics_Mregular.setDescription(r"""maximum regularization factor
\bi
\item if \inputvar{Mregular.ge.2}, then \internal{regumm}$_i$ = \inputvar{Mregular} $/ 2 $ where \internal{m}$_i > $ \inputvar{Mregular}
\ei""")

g_input_local = Group(g_input, "local")
g_input_local.setDescription(r"The namelist \verb+locallist+ controls the construction of the Beltrami fields in each volume.")

d_input_local_LBeltrami = Dataset(g_input_local, "LBeltrami", "int")
d_input_local_LBeltrami.setDefaultValue(4)
d_input_local_LBeltrami.setDescription(r"""\bi
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
\ei""")

d_input_local_Linitgues = Dataset(g_input_local, "Linitgues", "int")
d_input_local_Linitgues.setDefaultValue(1)
d_input_local_Linitgues.setDescription(r"""controls how initial guess for Beltrami field is constructed;
\bi
\item only relevant for routines that require an initial guess for the Beltrami fields, such as the SQP and Newton methods,
or the sparse linear solver;
\item if \inputvar{Linitgues = 0}, the initial guess for the Beltrami field is trivial;
\item if \inputvar{Linitgues = 1}, the initial guess for the Beltrami field is an integrable approximation;
\item if \inputvar{Linitgues = 2}, the initial guess for the Beltrami field is read from file; 
\item if \inputvar{Linitgues = 3}, the initial guess for the Beltrami field will be randomized with the maximum \inputvar{maxrndgues};
\ei""")

d_input_local_maxrndgues = Dataset(g_input_local, "maxrndgues", "double")
d_input_local_maxrndgues.setDefaultValue(1.0)
d_input_local_maxrndgues.setDescription(r"the maximum random number of the Beltrami field if \inputvar{Linitgues = 3}")

d_input_local_Lposdef = Dataset(g_input_local, "Lposdef", "int")
d_input_local_Lposdef.setDefaultValue(0)
d_input_local_Lposdef.setDescription("redundant")

g_input_global = Group(g_input, "global")
g_input_global.setDescription(r"The namelist \verb+globallist+ controls the search for global force-balance")

d_input_global_Lfindzero = Dataset(g_input_global, "Lfindzero", "int")
d_input_global_Lfindzero.setDefaultValue(0)
d_input_global_Lfindzero.setDescription(r"""use Newton methods to find zero of force-balance, which is computed by \link{dforce};
\bi
\item[o.] if \inputvar{Lfindzero = 0}, then \link{dforce} is called once 
          to compute the Beltrami fields consistent with the given geometry and constraints;
\item[i.] if \inputvar{Lfindzero = 1}, then call
          \nag{www.nag.co.uk/numeric/FL/manual19/pdf/C05/c05ndf_fl19.pdf}{C05NDF} (uses function values only),
          which iteratively calls \link{dforce};
\item[ii.] if \inputvar{Lfindzero = 2}, then call
           \nag{www.nag.co.uk/numeric/FL/manual19/pdf/C05/c05pdf_fl19.pdf}{C05PDF} (uses derivative information),
           which iteratively calls \link{dforce};
\ei""")

d_input_global_escale = Dataset(g_input_global, "escale", "double")
d_input_global_escale.setDefaultValue(0.0)
d_input_global_escale.setDescription(r"""controls the weight factor, \type{BBweight}, in the force-imbalance harmonics;
\bi
\item[i.] \type{BBweight(i)} $\ds \equiv \inputvar{opsilon} \times \exp\left[-\inputvar{escale} \times (m_i^2+n_i^2) \right]$
\item[ii.] defined in \link{preset}; used in \link{dforce};
\item[iii.] also see \Eqn{forcebalancemn} below;
\ei""")

d_input_global_opsilon = Dataset(g_input_global, "opsilon", "double")
d_input_global_opsilon.setDefaultValue(1.0)
d_input_global_opsilon.setDescription(r"""weighting of force-imbalance; 
\bi
\item[i.] used in \link{dforce}; also see \Eqn{forcebalancemn} below;
\ei""")

d_input_global_pcondense = Dataset(g_input_global, "pcondense", "double")
d_input_global_pcondense.setDefaultValue(2.0)
d_input_global_pcondense.setDescription(r"""spectral condensation parameter; 
\bi
\item[i.] used in \link{preset} to define \type{mmpp(i)} $\equiv m_i^p$, where $p\equiv $ \inputvar{pcondense};
\item[ii.] the angle freedom is exploited to minimize $\ds \inputvar{epsilon} \sum_{i} m_i^p (R_{i}^2+Z_{i}^2)$
      with respect to tangential variations in the interface geometry;
\item[ii.] also see \Eqn{spectralbalancemn} below;
\ei""")

d_input_global_epsilon = Dataset(g_input_global, "epsilon", "double")
d_input_global_epsilon.setDefaultValue(0.0)
d_input_global_epsilon.setDescription(r"""weighting of spectral-width constraint
\bi
\item[i.] used in \link{dforce}; also see \Eqn{spectralbalancemn} below
\ei""")

d_input_global_wpoloidal = Dataset(g_input_global, "wpoloidal", "double")
d_input_global_wpoloidal.setDefaultValue(1.0)
d_input_global_wpoloidal.setDescription(r"""``star-like'' poloidal angle constraint radial exponential factor;
used in \link{preset} to construct \type{sweight}""")

d_input_global_upsilon = Dataset(g_input_global, "upsilon", "double")
d_input_global_upsilon.setDefaultValue(1.0)
d_input_global_upsilon.setDescription(r"""weighting of ``star-like'' poloidal angle constraint;
used in \link{preset} to construct \type{sweight};""")

d_input_global_forcetol = Dataset(g_input_global, "forcetol", "double")
d_input_global_forcetol.setDefaultValue(1.0e-10)
d_input_global_forcetol.setDescription(r"""required tolerance in force-balance error; only used as an initial check;
\bi
\item[i.]   if the initially supplied interfaces are consistent with force-balance to within \inputvar{forcetol},
            then the geometry of the interfaces is not altered;
\item[ii.]  if not, then the geometry of the interfaces is changed in order to bring the configuration into forcebalance
            so that the geometry of interfaces is within \inputvar{c05xtol}, defined below, of the true solution;
\item[iii.] to force execution of either \nag{www.nag.co.uk/numeric/FL/manual19/pdf/C05/c05ndf_fl19.pdf}{C05NDF}
            or \nag{www.nag.co.uk/numeric/FL/manual19/pdf/C05/c05pdf_fl19.pdf}{C05PDF}, regardless of the initial force imbalance, 
            set \inputvar{forcetol < 0};
\ei""")

d_input_global_c05xmax = Dataset(g_input_global, "c05xmax", "double")
d_input_global_c05xmax.setDefaultValue(1.0e-6)
d_input_global_c05xmax.setDescription(r"required tolerance in position, ${\bf x} \equiv \{ R_{i,v}, Z_{i,v}\}$;")

d_input_global_c05xtol = Dataset(g_input_global, "c05xtol", "double")
d_input_global_c05xtol.setDefaultValue(1.0e-12)
d_input_global_c05xtol.setDescription(r"""required tolerance in position, ${\bf x} \equiv \{ R_{i,v}, Z_{i,v}\}$;
\bi
\item[i.] used by both \nag{www.nag.co.uk/numeric/FL/manual19/pdf/C05/c05ndf_fl19.pdf}{C05NDF} and 
          \nag{www.nag.co.uk/numeric/FL/manual19/pdf/C05/c05pdf_fl19.pdf}{C05PDF};
          see the NAG documents for further details on how the error is defined;
\item[ii.] constraint \inputvar{c05xtol.gt.0.0};
\ei""")

d_input_global_c05factor = Dataset(g_input_global, "c05factor", "double")
d_input_global_c05factor.setDefaultValue(1.0e-2)
d_input_global_c05factor.setDescription(r"""used to control initial step size in
      \nag{www.nag.co.uk/numeric/FL/manual19/pdf/C05/c05ndf_fl19.pdf}{C05NDF} and 
      \nag{www.nag.co.uk/numeric/FL/manual19/pdf/C05/c05pdf_fl19.pdf}{C05PDF};
\bi
\item[i.] constraint \inputvar{c05factor.gt.0.0};
\item[ii.] only relevant if \inputvar{Lfindzero.gt.0};
\ei""")

d_input_global_LreadGF = Dataset(g_input_global, "LreadGF", "boolean")
d_input_global_LreadGF.setDefaultValue(True)
d_input_global_LreadGF.setDescription(r"""read $\nabla_{\bf x} {\bf F}$ from file \type{.GF}; 
\bi
\item[i.] only used if \inputvar{Lfindzero = 2}; 
\item[ii.] only used in \link{newton};
\ei""")

d_input_global_mfreeits = Dataset(g_input_global, "mfreeits", "int")
d_input_global_mfreeits.setDefaultValue(0)
d_input_global_mfreeits.setDescription(r"""maximum allowed free-boundary iterations;
\bi
\item[i.] only used if \inputvar{Lfreebound = 1}; 
\item[ii.] only used in \link{xspech};
\ei""")

d_input_global_bnstol = Dataset(g_input_global, "bnstol", "double")
d_input_global_bnstol.setDefaultValue(1.0e-6)
d_input_global_bnstol.setDescription("redundant")

d_input_global_bnsblend = Dataset(g_input_global, "bnsblend", "double")
d_input_global_bnsblend.setDefaultValue(0.666)
d_input_global_bnsblend.setDescription("redundant")

d_input_global_gBntol = Dataset(g_input_global, "gBntol", "double")
d_input_global_gBntol.setDefaultValue(1.0e-6)
d_input_global_gBntol.setDescription(r"""equired tolerance in free-boundary iterations;
\bi
\item[i.] only used if \inputvar{Lfreebound = 1}; 
\item[ii.] only used in \link{xspech}; see \link{xspech} for more documentation;
\ei""")

d_input_global_gBnbld = Dataset(g_input_global, "gBnbld", "double")
d_input_global_gBnbld.setDefaultValue(0.666)
d_input_global_gBnbld.setDescription(r"""normal blend;
\bi
\item[i.] The ``new'' magnetic field at the computational boundary produced by the plasma currents is updated using a Picard scheme:
          \be ({\bf B}\cdot{\bf n})^{j+1} =    \inputvar{gBnbld}  \times ({\bf B}\cdot{\bf n})^{j} 
                                          + (1-\inputvar{gBnbld}) \times ({\bf B}\cdot{\bf n})^{*},
          \ee
          where $j$ labels free-boundary iterations, and $({\bf B}\cdot{\bf n})^{*}$ is computed by virtual casing.
\item[ii.] only used if \inputvar{Lfreebound = 1}; 
\item[ii.] only used in \link{xspech};
\ei""")

d_input_global_vcasingeps = Dataset(g_input_global, "vcasingeps", "double")
d_input_global_vcasingeps.setDefaultValue(1.0e-12)
d_input_global_vcasingeps.setDescription(r"regularization of Biot-Savart; see \link{bnorml}, \link{casing}")

d_input_global_vcasingtol = Dataset(g_input_global, "vcasingtol", "double")
d_input_global_vcasingtol.setDefaultValue(1.0e-8)
d_input_global_vcasingtol.setDescription(r"accuracy on virtual casing integral; see \link{bnorml}, \link{casing}")

d_input_global_vcasingits = Dataset(g_input_global, "vcasingits", "int")
d_input_global_vcasingeps.setDefaultValue(8)
d_input_global_vcasingeps.setDescription(r"minimum number of calls to adaptive virtual casing routine; see \link{casing}")

d_input_global_vcasingper = Dataset(g_input_global, "vcasingper", "int")
d_input_global_vcasingeps.setDefaultValue(1)
d_input_global_vcasingeps.setDescription(r"periods of integragion  in adaptive virtual casing routine; see \link{casing}")

d_input_global_mcasingcal = Dataset(g_input_global, "mcasingcal", "int")
d_input_global_vcasingeps.setDefaultValue(8)
d_input_global_vcasingeps.setDescription(r"minimum number of calls to adaptive virtual casing routine; see \link{casing}")

g_input_diagnostics = Group(g_input, "diagnostics")
g_input_diagnostics.setDescription(r"The namelist \type{diagnosticslist} controls post-processor diagnostics, such as \Poincare plot resolution, $\dots$,...")

d_input_diagnostics_odetol = Dataset(g_input_diagnostics, "odetol", "double")
d_input_diagnostics_odetol.setDefaultValue(1.0e-7)
d_input_diagnostics_odetol.setDescription("o.d.e. integration tolerance for all field line tracing routines")

d_input_diagnostics_absreq = Dataset(g_input_diagnostics, "absreq", "double")
d_input_diagnostics_absreq.setDefaultValue(1.0e-8)
d_input_diagnostics_absreq.setDescription("redundant")

d_input_diagnostics_relreq = Dataset(g_input_diagnostics, "relreq", "double")
d_input_diagnostics_relreq.setDefaultValue(1.0e-8)
d_input_diagnostics_relreq.setDescription("redundant")

d_input_diagnostics_absacc = Dataset(g_input_diagnostics, "absacc", "double")
d_input_diagnostics_absacc.setDefaultValue(1.0e-4)
d_input_diagnostics_absacc.setDescription("redundant")

d_input_diagnostics_epsr = Dataset(g_input_diagnostics, "epsr", "double")
d_input_diagnostics_epsr.setDefaultValue(1.0e-6)
d_input_diagnostics_epsr.setDescription("redundant")

d_input_diagnostics_nPpts = Dataset(g_input_diagnostics, "nPpts", "int")
d_input_diagnostics_nPpts.setDefaultValue(0)
d_input_diagnostics_nPpts.setDescription(r"""number of toroidal transits used (per trajectory) in following field lines
for constructing \Poincare plots;
if \inputvar{nPpts<1}, no \Poincare plot is constructed;""")

d_input_diagnostics_nPtrj = Dataset(g_input_diagnostics, "nPtrj", "int", 1)
d_input_diagnostics_nPtrj.setDefaultValue(-1)
d_input_diagnostics_nPtrj.setDescription(r"""number of trajectories in each annulus to be followed in constructing \Poincare plot;
\bi
\item if \inputvar{nPtrj(l)<0}, then \inputvar{nPtrj(l) = Ni(l)},
      where \type{Ni(l)} is the grid resolution used to construct the Beltrami field in volume $l$;
\ei""")

d_input_diagnostics_LHevalues = Dataset(g_input_diagnostics, "LHevalues", "boolean")
d_input_diagnostics_LHevalues.setDefaultValue(False)
d_input_diagnostics_LHevalues.setDescription(r"to compute eigenvalues of $\nabla {\bf F}$")

d_input_diagnostics_LHevectors = Dataset(g_input_diagnostics, "LHevectors", "boolean")
d_input_diagnostics_LHevectors.setDefaultValue(False)
d_input_diagnostics_LHevectors.setDescription(r"to compute eigenvectors (and also eigenvalues) of $\nabla {\bf F}$")

d_input_diagnostics_LHmatrix = Dataset(g_input_diagnostics, "LHmatrix", "boolean")
d_input_diagnostics_LHmatrix.setDefaultValue(False)
d_input_diagnostics_LHmatrix.setDescription(r"to compute and write to file the elements of $\nabla {\bf F}$")

d_input_diagnostics_Lperturbed = Dataset(g_input_diagnostics, "Lperturbed", "int")
d_input_diagnostics_Lperturbed.setDefaultValue(0)
d_input_diagnostics_Lperturbed.setDescription("to compute linear, perturbed equilibrium")

d_input_diagnostics_dpp = Dataset(g_input_diagnostics, "dpp", "int")
d_input_diagnostics_dpp.setDefaultValue(1)
d_input_diagnostics_dpp.setDescription("perturbed harmonic")

d_input_diagnostics_dqq = Dataset(g_input_diagnostics, "dqq", "int")
d_input_diagnostics_dqq.setDefaultValue(1)
d_input_diagnostics_dqq.setDescription("perturbed harmonic")

d_input_diagnostics_Lcheck = Dataset(g_input_diagnostics, "Lcheck", "int")
d_input_diagnostics_Lcheck.setDefaultValue(0)
d_input_diagnostics_Lcheck.setDescription(r"""implement various checks;
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
\ei""")

d_input_diagnostics_Ltiming = Dataset(g_input_diagnostics, "Ltiming", "boolean")
d_input_diagnostics_Ltiming.setDefaultValue(True)
d_input_diagnostics_Ltiming.setDescription("to check timing")

d_input_diagnostics_fudge = Dataset(g_input_diagnostics, "fudge", "double")
d_input_diagnostics_fudge.setDefaultValue(1.0)
d_input_diagnostics_fudge.setDescription("redundant")

d_input_diagnostics_scaling = Dataset(g_input_diagnostics, "scaling", "double")
d_input_diagnostics_scaling.setDefaultValue(1.0)
d_input_diagnostics_scaling.setDescription("redundant")

g_output = Group(s, "output")
g_output.setDescription("output data; content of the previous HDF5 file")

d_output_Vns = Dataset(g_output, "Vns", "double", 1)
d_output_Vns.setDescription("stellarator symmetric normal field at boundary; vacuum component")

d_output_Bns = Dataset(g_output, "Bns", "double", 1)
d_output_Bns.setDescription("stellarator symmetric normal field at boundary; plasma component")

d_output_Vnc = Dataset(g_output, "Vnc", "double", 1)
d_output_Vnc.setDescription("non-stellarator symmetric normal field at boundary; vacuum component")

d_output_Bnc = Dataset(g_output, "Bnc", "double", 1)
d_output_Bnc.setDescription("non-stellarator symmetric normal field at boundary; plasma component")

d_output_mn = Dataset(g_output, "mn", "int")
d_output_mn.setDescription("number of Fourier modes")

d_output_im = Dataset(g_output, "im", "int", 1)
d_output_im.setDescription("poloidal mode numbers")

d_output_in = Dataset(g_output, "in", "int", 1)
d_output_in.setDescription("toroidal mode numbers")

d_output_Mvol = Dataset(g_output, "Mvol", "int")
d_output_Mvol.setDescription("number of interfaces = number of volumes")

d_output_Rbc = Dataset(g_output, "Rbc", "double", 2)
d_output_Rbc.setDescription(r"stellarator symmetric Fourier harmonics, $R_{m,n}$, of interfaces")

d_output_Zbs = Dataset(g_output, "Zbs", "double", 2)
d_output_Zbs.setDescription(r"stellarator symmetric Fourier harmonics, $Z_{m,n}$, of interfaces")

d_output_Rbs = Dataset(g_output, "Rbs", "double", 2)
d_output_Rbs.setDescription(r"non-stellarator symmetric Fourier harmonics, $R_{m,n}$, of interfaces")

d_output_Zbc = Dataset(g_output, "Zbc", "double", 2)
d_output_Zbc.setDescription(r"non-stellarator symmetric Fourier harmonics, $Z_{m,n}$, of interfaces")

d_output_ForceErr = Dataset(g_output, "ForceErr", "double")
d_output_ForceErr.setDescription("force-balance error across interfaces")

d_output_adiabatic = Dataset(g_output, "adiabatic", "double", 1)
d_output_helicity = Dataset(g_output, "helicity", "double", 1)
d_output_mu = Dataset(g_output, "mu", "double", 1)
d_output_tflux = Dataset(g_output, "tflux", "double", 1)
d_output_pflux = Dataset(g_output, "pflux", "double", 1)

d_output_volume = Dataset(g_output, "volume", "double")

d_output_Mrad = Dataset(g_output, "Mrad", "int")
d_output_Mrad.setDescription("the maximum radial (Chebyshev) resolution")

d_output_TT = Dataset(g_output, "TT", "double", 3)
d_output_TT.setDescription(r"the Chebyshev polynomials, $T_l$, and their derivatives, evaluated at $s=\pm 1$")

d_output_Btemn = Dataset(g_output, "Btemn", "double", 3)
d_output_Btemn.setDescription(r"""the cosine harmonics of the covariant poloidal field,
i.e. $[[B_{\t,j}]]$ evaluated on the inner and outer interface in each volume""")

d_output_Bzemn = Dataset(g_output, "Bzemn", "double", 3)
d_output_Bzemn.setDescription(r"""the cosine harmonics of the covariant toroidal field,
i.e. $[[B_{\z,j}]]$ evaluated on the inner and outer interface in each volume""")

d_output_Btomn = Dataset(g_output, "Btomn", "double", 3)
d_output_Btomn.setDescription(r"""the sine harmonics of the covariant poloidal field,
i.e. $[[B_{\t,j}]]$ evaluated on the inner and outer interface in each volume""")

d_output_Bzomn = Dataset(g_output, "Bzomn", "double", 3)
d_output_Bzomn.setDescription(r"""the sine harmonics of the covariant toroidal field,
i.e. $[[B_{\t,j}]]$ evaluated on the inner and outer interface in each volume""")

d_output_lmns = Dataset(g_output, "lmns", "double")
d_output_lmns.setDescription("resolution of the straight fieldline transformation")

g_vectorPotential = Group(s, "vector_potential")
g_vectorPotential.setDescription(r"""The covariant components of the vector potential are written as
\be            A_\t & = & \sum_i \sum_{l=0}^L \Ate{i,l} \; T_{l}(s) \cos\a_i + \sum_i \sum_{l=0}^L \Ato{i,l} \; T_{l}(s) \sin\a_i \\
               A_\z & = & \sum_i \sum_{l=0}^L \Aze{i,l} \; T_{l}(s) \cos\a_i + \sum_i \sum_{l=0}^L \Azo{i,l} \; T_{l}(s) \sin\a_i ,
\ee
where $T_l(s)$ are the Chebyshev polynomials and $\a_i \equiv m_i \t - n_i \z$.
The following internal arrays are declared in \link{preset}
\verb{dAte(0,i)%s(l){$\equiv \Ate{i,l}$
\verb{dAze(0,i)%s(l){$\equiv \Aze{i,l}$
\verb{dAto(0,i)%s(l){$\equiv \Ato{i,l}$
\verb{dAzo(0,i)%s(l){$\equiv \Azo{i,l}$""")

d_vectorPotential_Ate = Dataset(g_vectorPotential, "Ate", "double", 2)
d_vectorPotential_Aze = Dataset(g_vectorPotential, "Aze", "double", 2)
d_vectorPotential_Ato = Dataset(g_vectorPotential, "Ato", "double", 2)
d_vectorPotential_Azo = Dataset(g_vectorPotential, "Azo", "double", 2)

#s.inventory()



#%% generate Fortran reading module
from genFortran import genFortranReader

# target Fortran module
moduleName   = "read_spec"
outdir = "/home/jonathan/Uni/04_PhD/00_programs/SPEC/SPEC/Utilities/" # read_spec.f90

genFortranReader(outdir, moduleName, s)

#%% generate Java reading class
from genJava import genJavaReader

# target Java class
packageName = "ipp.w7x.equi"
className   = "SpecOutput"
outdir = "/home/jonathan/Uni/04_PhD/00_programs/SPEC/SPEC/Utilities/" # SpecOutput.java

genJavaReader(outdir, packageName, className, s)



