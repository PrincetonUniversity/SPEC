#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is the definition script which controls the generation of the input file
and output file reading and writing routines for the SPEC MRxMHD
equilibrium code.

In the first part of this script, the input file format
and the output file format are declared in an abstract way.

Based on these declarations, reading and writing routines in various
programming languages are generated automagically.
The code generation functionality is encapsulated in the
Interface Definition Framework, which can be found at:
https://github.com/jonathanschilling/idf

If you want to have an input or output variable added to SPEC, this script
is the place to put it. Then run it, commit your changed source code to Git
and compile SPEC again to use your changes.

In order to get the correct order of the comments, you should use Python >=3.7.
From Python 3.7 on, it is guaranteed that the insertion order of dict()
items is kept and this is relied on in this code when iterating over the keys
of the documentation parts.
More on this: https://docs.python.org/3/whatsnew/3.7.html

@author: Jonathan Schilling (jonathan.schilling@ipp.mpg.de)
"""

import os # linesep

# interface definition framework for source code generation
idfPath = "/home/jonathan/work/code/idf"
import sys
if not idfPath in sys.path:
    sys.path.insert(0, idfPath)

try:
    from idf import Variable, indented, toDoc, get_creation_tag, relname
    from idf import Fortran
except ImportError:
    raise ImportError("'idf' python package missing")

# define the version of SPEC
version = Variable("version")
version.setDescription(r"version of SPEC")
version.setType("int")
version.setRank(1)
version.setMaximumIndices(["3"])
version.setIsParameter(True)
version.setDefaultValue([3,0,0])

# define the input quantities for SPEC

MNvol = Variable("MNvol")
MNvol.setDescription(r"maximum value of \c Nvol")
MNvol.setType("int")
MNvol.setDefaultValue(256)
MNvol.setIsParameter(True)

MMpol = Variable("MMpol")
MMpol.setDescription(r"maximum value of \c Mpol")
MMpol.setType("int")
MMpol.setDefaultValue(64)
MMpol.setIsParameter(True)

MNtor = Variable("MNtor")
MNtor.setDescription(r"maximum value of \c Ntor")
MNtor.setType("int")
MNtor.setDefaultValue(64)
MNtor.setIsParameter(True)

params_maxDims = [
        MNvol,
        MMpol,
        MNtor
        ]

###############################################################################
# physicslist 
###############################################################################
input_physics_Igeometry = Variable("Igeometry")
input_physics_Igeometry.setDescription({r"selects Cartesian, cylindrical or toroidal geometry":
                                        [r"\c Igeometry=1 : Cartesian; geometry determined by \f$R\f$",
                                         r"\c Igeometry=2 : cylindrical; geometry determined by \f$R\f$",
                                         r"\c Igeometry=3 : toroidal; geometry determined by \f$R\f$ *and* \f$Z\f$"]
                                        })
input_physics_Igeometry.setType("int")
input_physics_Igeometry.setDefaultValue(3)

input_physics_Istellsym = Variable("Istellsym")
input_physics_Istellsym.setDescription(r"stellarator symmetry is enforced if \c Istellsym=1")
input_physics_Istellsym.setType("int")
input_physics_Istellsym.setDefaultValue(1)

input_physics_Lfreebound = Variable("Lfreebound")
input_physics_Lfreebound.setDescription(r"compute vacuum field surrounding plasma")
input_physics_Lfreebound.setType("int")
input_physics_Lfreebound.setDefaultValue(0)

input_physics_phiedge = Variable("phiedge")
input_physics_phiedge.setDescription(r"total enclosed toroidal magnetic flux")
input_physics_phiedge.setType("double")
input_physics_phiedge.setDefaultValue(1.0)
input_physics_phiedge.setUnit("Vs")

input_physics_curtor = Variable("curtor")
input_physics_curtor.setDescription(r"total enclosed (toroidal) plasma current")
input_physics_curtor.setType("double")
input_physics_curtor.setDefaultValue(0.0)

input_physics_curpol = Variable("curpol")
input_physics_curpol.setDescription(r"total enclosed (poloidal) linking current")
input_physics_curpol.setType("double")
input_physics_curpol.setDefaultValue(0.0)

input_physics_gamma = Variable("gamma")
input_physics_gamma.setDescription(r"adiabatic index; cannot set \f$|\gamma| = 1\f$")
input_physics_gamma.setType("double")
input_physics_gamma.setDefaultValue(0.0)

input_physics_Nfp = Variable("Nfp")
input_physics_Nfp.setDescription({r"field periodicity":
                                  [r"all Fourier representations are of the form \f$\cos(m\theta-n N \zeta)\f$, \f$\sin(m\theta-n N \zeta)\f$, where \f$N\equiv\f$\c Nfp",
                                   r"constraint: \c Nfp >= 1"]
                                  })
input_physics_Nfp.setType("int")
input_physics_Nfp.setDefaultValue(1)

input_physics_Nvol = Variable("Nvol")
input_physics_Nvol.setDescription({r"number of volumes":
                                   [r"each volume \f${\cal V}_l\f$ is bounded by the \f${\cal I}_{l-1}\f$ and \f${\cal I}_{l}\f$ interfaces",
                                    r"note that in cylindrical or toroidal geometry, \f${\cal I}_{0}\f$ is the degenerate coordinate axis",
                                    r"constraint: \c Nvol<=MNvol"]
                                   })
input_physics_Nvol.setType("int")
input_physics_Nvol.setDefaultValue(1)

input_physics_Mpol = Variable("Mpol")
input_physics_Mpol.setDescription({r"number of poloidal Fourier harmonics":
                                   [ r"all Fourier representations of doubly-periodic functions are of the form"+os.linesep
                                    +r"\f{eqnarray}{ f(\theta,\zeta) & = & \sum_{n=0}^{\texttt{Ntor}} f_{0,n}\cos(-n \, \texttt{Nfp} \, \zeta)"+os.linesep
                                    +r"\sum_{m=1}^{\texttt{Mpol}}\sum_{n=\texttt{-Ntor}}^{\texttt{Ntor}} f_{m,n}\cos(m\theta-n \, \texttt{Nfp} \, \zeta),"+os.linesep
                                    +r"\f}"+os.linesep
                                    +r"Internally these \"double\" summations are written as a \"single\" summation,"+os.linesep
                                    +r"e.g. \f$f(\theta,\zeta) = \sum_j f_j \cos(m_j\theta-n_j\zeta)\f$."]
                                   })
input_physics_Mpol.setType("int")
input_physics_Mpol.setDefaultValue(0)

input_physics_Ntor = Variable("Ntor")
input_physics_Ntor.setDescription({r"number of toroidal Fourier harmonics":
                                   [ r"all Fourier representations of doubly-periodic functions are of the form"+os.linesep
                                    +r"\f{eqnarray}{ f(\theta,\zeta) & = & \sum_{n=0}^{\texttt{Ntor}} f_{0,n}\cos(-n \, \texttt{Nfp} \, \zeta)"+os.linesep
                                    +r"\sum_{m=1}^{\texttt{Mpol}}\sum_{n=\texttt{-Ntor}}^{\texttt{Ntor}} f_{m,n}\cos(m\theta-n \, \texttt{Nfp} \, \zeta),"+os.linesep
                                    +r"\f}"+os.linesep
                                    +r"Internally these \"double\" summations are written as a \"single\" summation,"+os.linesep
                                    +r"e.g. \f$f(\theta,\zeta) = \sum_j f_j \cos(m_j\theta-n_j\zeta)\f$."]
                                   })
input_physics_Ntor.setType("int")
input_physics_Ntor.setDefaultValue(0)

input_physics_Lrad = Variable("Lrad")
input_physics_Lrad.setDescription({r"Chebyshev resolution in each volume":
                                   [r"constraint: \c Lrad(1:Mvol) >= 2"]
                                   })
input_physics_Lrad.setType("int")
input_physics_Lrad.setRank(1)
input_physics_Lrad.setDefaultValue(4)
input_physics_Lrad.setMaximumIndices(["MNvol+1"])

input_physics_Lconstraint = Variable("Lconstraint")
input_physics_Lconstraint.setDescription({r"selects constraints; primarily used in ma02aa() and mp00ac()":
                                          [ r"if \c Lconstraint=-1, then in the plasma regions \f$\Delta\psi_t\f$, \f$\mu\f$ and \f$\Delta \psi_p\f$ are *not* varied"+os.linesep
                                           +r"and in the vacuum region (only for free-boundary) \f$\Delta\psi_t\f$ and \f$\Delta \psi_p\f$ are *not* varied, and \f$\mu = 0\f$",
                                            r"if \c Lconstraint=0, then in the plasma regions \f$\Delta\psi_t\f$, \f$\mu\f$ and \f$\Delta \psi_p\f$ are *not* varied"+os.linesep
                                           +r"and in the vacuum region (only for free-boundary) \f$\Delta\psi_t\f$ and \f$\Delta \psi_p\f$ are varied to match the"+os.linesep
                                           +r"prescribed plasma current, \c curtor, and the \"linking\" current, \c curpol, and \f$\mu = 0\f$",
                                            r"if \c Lconstraint=1, then in the plasma regions \f$\mu\f$ and \f$\Delta\psi_p\f$ are adjusted"+os.linesep
                                           +r"in order to satisfy the inner and outer interface transform constraints"+os.linesep
                                           +r"(except in the simple torus, where the enclosed poloidal flux is irrelevant,"+os.linesep
                                           +r"and only \f$\mu\f$ is varied to satisfy the outer interface transform constraint);"+os.linesep
                                           +r"and in the vacuum region \f$\Delta\psi_t\f$ and \f$\Delta \psi_p\f$ are varied to match the transform constraint on the boundary"+os.linesep
                                           +r"and to obtain the prescribed linking current, \c curpol, and \f$\mu = 0\f$",
                                            r"\todo if \c Lconstraint=2, under reconstruction"]
                                          })
input_physics_Lconstraint.setType("int")
input_physics_Lconstraint.setDefaultValue(-1)

input_physics_tflux = Variable("tflux")
input_physics_tflux.setDescription({r"toroidal flux, \f$\psi_t\f$, enclosed by each interface":
                                    [ r"For each of the plasma volumes, this is a constraint: \c tflux is *not* varied",
                                      r"For the vacuum region (only if \c Lfreebound==1), \c tflux  may be allowed to vary to match constraints",
                                      r"Note that \c tflux  will be normalized so that \c tflux(Nvol) = 1.0,"+os.linesep
                                     +r"so that \c tflux  is arbitrary up to a scale factor",
                                      r"\sa phiedge"]
                                    })
input_physics_tflux.setUnit("Wb")
input_physics_tflux.setType("double")
input_physics_tflux.setRank(1)
input_physics_tflux.setDefaultValue(0.0)
input_physics_tflux.setMaximumIndices(["MNvol+1"])

input_physics_pflux = Variable("pflux")
input_physics_pflux.setDescription(r"poloidal flux, \f$\psi_p\f$, enclosed by each interface")
input_physics_pflux.setUnit("Wb")
input_physics_pflux.setType("double")
input_physics_pflux.setRank(1)
input_physics_pflux.setDefaultValue(0.0)
input_physics_pflux.setMaximumIndices(["MNvol+1"])

input_physics_helicity = Variable("helicity")
input_physics_helicity.setDescription(r"helicity, \f${\cal K}\f$, in each volume, \f${\cal V}_i\f$")
input_physics_helicity.setType("double")
input_physics_helicity.setRank(1)
input_physics_helicity.setDefaultValue(0.0)
input_physics_helicity.setMaximumIndices(["MNvol+1"])

input_physics_pscale = Variable("pscale")
input_physics_pscale.setDescription({r"pressure scale factor":
                                     [r"the initial pressure profile is given by \c pscale  \f$*\f$ \c pressure"]
                                    })
input_physics_pscale.setType("double")
input_physics_pscale.setDefaultValue(0.0)

input_physics_pressure = Variable("pressure")
input_physics_pressure.setDescription({r"pressure in each volume":
                                       [ r"The pressure is *not* held constant, but \f$p_l V_l^\gamma = P_l\f$ *is* held constant,"+os.linesep
                                        +r"where \f$P_l\f$ is determined by the initial pressures and the initial volumes, \f$V_l\f$.",
                                         r"Note that if \c gamma==0.0, then \f$p_l \equiv P_l\f$.",
                                         r"On output, the pressure is given by \f$p_l = P_l/V_l^\gamma\f$, where \f$V_l\f$ is the final volume.",
                                         r"\c pressure is only used in calculation of interface force-balance."]
                                       })
input_physics_pressure.setType("double")
input_physics_pressure.setRank(1)
input_physics_pressure.setDefaultValue(0.0)
input_physics_pressure.setMaximumIndices(["MNvol+1"])

input_physics_Ladiabatic = Variable("Ladiabatic")
input_physics_Ladiabatic.setDescription({r"logical flag":
                                         [r"If \c Ladiabatic==0, the adiabatic constants are determined by the initial pressure and volume.",
                                          r"If \c Ladiabatic==1, the adiabatic constants are determined by the given input \c adiabatic."]
                                         })
input_physics_Ladiabatic.setType("int")
input_physics_Ladiabatic.setDefaultValue(0)

input_physics_adiabatic = Variable("adiabatic")
input_physics_adiabatic.setDescription({r"adiabatic constants in each volume":
                                        [r"The pressure is *not* held constant, but \f$p_l V_l^\gamma = P_l \equiv\f$\c adiabatic is constant.",
                                         r"Note that if \c gamma==0.0, then \c pressure==adiabatic.",
                                         r"\c pressure is only used in calculation of interface force-balance."]
                                        })
input_physics_adiabatic.setType("double")
input_physics_adiabatic.setRank(1)
input_physics_adiabatic.setDefaultValue(0.0)
input_physics_adiabatic.setMaximumIndices(["MNvol+1"])

input_physics_mu = Variable("mu")
input_physics_mu.setDescription(r"helicity-multiplier, \f$\mu\f$, in each volume")
input_physics_mu.setType("double")
input_physics_mu.setRank(1)
input_physics_mu.setDefaultValue(0.0)
input_physics_mu.setMaximumIndices(["MNvol+1"])

input_physics_Ivolume = Variable("Ivolume")
input_physics_Ivolume.setDescription( r"Toroidal current constraint normalized by \f$\mu_0\f$ (\f$I_{volume} = \mu_0\cdot [A]\f$), in each volume. "+os.linesep
                                     +r"This is a cumulative quantity: \f$I_{\mathcal{V},i} = \int_0^{\psi_{t,i}} \mathbf{J}\cdot\mathbf{dS}\f$. "+os.linesep
                                     +r"Physically, it represents the sum of all non-pressure driven currents.")
input_physics_Ivolume.setType("double")
input_physics_Ivolume.setRank(1)
input_physics_Ivolume.setDefaultValue(0.0)
input_physics_Ivolume.setMaximumIndices(["MNvol+1"])

input_physics_Isurf = Variable("Isurf")
input_physics_Isurf.setDescription(r"Toroidal current normalized by \f$\mu_0\f$ at each interface (cumulative). This is the sum of all pressure driven currents.")
input_physics_Isurf.setType("double")
input_physics_Isurf.setRank(1)
input_physics_Isurf.setDefaultValue(0.0)
input_physics_Isurf.setMaximumIndices(["MNvol+1"])

input_physics_pl = Variable("pl")
input_physics_pl.setDescription( r"\"inside\" interface rotational-transform is \f$\mbox{$\,\iota\!\!$-} = (p_l+\gamma p_r)/(q_l+\gamma q_r)\f$,"+os.linesep
                                +r"where \f$\gamma\f$ is the golden mean, \f$\gamma = (1 + \sqrt 5 ) / 2 \f$."+os.linesep
                                +r"If both \f$q_l = 0\f$ *and* \f$q_r = 0\f$, then the (inside) interface rotational-transform is defined by \c iota .")
input_physics_pl.setType("int")
input_physics_pl.setRank(1)
input_physics_pl.setDefaultValue(0)
input_physics_pl.setStartingIndices([r"0"])
input_physics_pl.setMaximumIndices(["MNvol"])

input_physics_ql = Variable("ql")
input_physics_ql.setDescription( r"\"inside\" interface rotational-transform is \f$\mbox{$\,\iota\!\!$-} = (p_l+\gamma p_r)/(q_l+\gamma q_r)\f$,"+os.linesep
                                +r"where \f$\gamma\f$ is the golden mean, \f$\gamma = (1 + \sqrt 5 ) / 2 \f$."+os.linesep
                                +r"If both \f$q_l = 0\f$ *and* \f$q_r = 0\f$, then the (inside) interface rotational-transform is defined by \c iota .")
input_physics_ql.setType("int")
input_physics_ql.setRank(1)
input_physics_ql.setDefaultValue(0)
input_physics_ql.setStartingIndices([r"0"])
input_physics_ql.setMaximumIndices(["MNvol"])

input_physics_pr = Variable("pr")
input_physics_pr.setDescription( r"\"inside\" interface rotational-transform is \f$\mbox{$\,\iota\!\!$-} = (p_l+\gamma p_r)/(q_l+\gamma q_r)\f$,"+os.linesep
                                +r"where \f$\gamma\f$ is the golden mean, \f$\gamma = (1 + \sqrt 5 ) / 2 \f$."+os.linesep
                                +r"If both \f$q_l = 0\f$ *and* \f$q_r = 0\f$, then the (inside) interface rotational-transform is defined by \c iota .")
input_physics_pr.setType("int")
input_physics_pr.setRank(1)
input_physics_pr.setDefaultValue(0)
input_physics_pr.setStartingIndices([r"0"])
input_physics_pr.setMaximumIndices(["MNvol"])

input_physics_qr = Variable("qr")
input_physics_qr.setDescription( r"\"inside\" interface rotational-transform is \f$\mbox{$\,\iota\!\!$-} = (p_l+\gamma p_r)/(q_l+\gamma q_r)\f$,"+os.linesep
                                +r"where \f$\gamma\f$ is the golden mean, \f$\gamma = (1 + \sqrt 5 ) / 2 \f$."+os.linesep
                                +r"If both \f$q_l = 0\f$ *and* \f$q_r = 0\f$, then the (inside) interface rotational-transform is defined by \c iota .")
input_physics_qr.setType("int")
input_physics_qr.setRank(1)
input_physics_qr.setDefaultValue(0)
input_physics_qr.setStartingIndices([r"0"])
input_physics_qr.setMaximumIndices(["MNvol"])

input_physics_iota = Variable("iota")
input_physics_iota.setDescription({r"rotational-transform, \f$\mbox{$\,\iota\!\!$-}\f$, on inner side of each interface":
                                   [r"only relevant if illogical input for \c ql and \c qr are provided"]
                                   })
input_physics_iota.setType("double")
input_physics_iota.setRank(1)
input_physics_iota.setDefaultValue(0.0)
input_physics_iota.setStartingIndices([r"0"])
input_physics_iota.setMaximumIndices(["MNvol"])

input_physics_lp = Variable("lp")
input_physics_lp.setDescription( r"\"outer\" interface rotational-transform is \f$\mbox{$\,\iota\!\!$-} = (p_l+\gamma p_r)/(q_l+\gamma q_r)\f$,"+os.linesep
                                +r"where \f$\gamma\f$ is the golden mean, \f$\gamma = (1 + \sqrt 5 ) / 2\f$."+os.linesep
                                +r"If both \f$q_l = 0\f$ *and* \f$q_r = 0\f$, then the (outer) interface rotational-transform is defined by \c oita .")
input_physics_lp.setType("int")
input_physics_lp.setRank(1)
input_physics_lp.setDefaultValue(0)
input_physics_lp.setStartingIndices([r"0"])
input_physics_lp.setMaximumIndices(["MNvol"])

input_physics_lq = Variable("lq")
input_physics_lq.setDescription( r"\"outer\" interface rotational-transform is \f$\mbox{$\,\iota\!\!$-} = (p_l+\gamma p_r)/(q_l+\gamma q_r)\f$,"+os.linesep
                                +r"where \f$\gamma\f$ is the golden mean, \f$\gamma = (1 + \sqrt 5 ) / 2\f$."+os.linesep
                                +r"If both \f$q_l = 0\f$ *and* \f$q_r = 0\f$, then the (outer) interface rotational-transform is defined by \c oita .")
input_physics_lq.setType("int")
input_physics_lq.setRank(1)
input_physics_lq.setDefaultValue(0)
input_physics_lq.setStartingIndices([r"0"])
input_physics_lq.setMaximumIndices(["MNvol"])

input_physics_rp = Variable("rp")
input_physics_rp.setDescription( r"\"outer\" interface rotational-transform is \f$\mbox{$\,\iota\!\!$-} = (p_l+\gamma p_r)/(q_l+\gamma q_r)\f$,"+os.linesep
                                +r"where \f$\gamma\f$ is the golden mean, \f$\gamma = (1 + \sqrt 5 ) / 2\f$."+os.linesep
                                +r"If both \f$q_l = 0\f$ *and* \f$q_r = 0\f$, then the (outer) interface rotational-transform is defined by \c oita .")
input_physics_rp.setType("int")
input_physics_rp.setRank(1)
input_physics_rp.setDefaultValue(0)
input_physics_rp.setStartingIndices([r"0"])
input_physics_rp.setMaximumIndices(["MNvol"])

input_physics_rq = Variable("rq")
input_physics_rq.setDescription( r"\"outer\" interface rotational-transform is \f$\mbox{$\,\iota\!\!$-} = (p_l+\gamma p_r)/(q_l+\gamma q_r)\f$,"+os.linesep
                                +r"where \f$\gamma\f$ is the golden mean, \f$\gamma = (1 + \sqrt 5 ) / 2\f$."+os.linesep
                                +r"If both \f$q_l = 0\f$ *and* \f$q_r = 0\f$, then the (outer) interface rotational-transform is defined by \c oita .")
input_physics_rq.setType("int")
input_physics_rq.setRank(1)
input_physics_rq.setDefaultValue(0)
input_physics_rq.setStartingIndices([r"0"])
input_physics_rq.setMaximumIndices(["MNvol"])

input_physics_oita = Variable("oita")
input_physics_oita.setDescription({r"rotational-transform, \f$\mbox{$\,\iota\!\!$-}\f$, on outer side of each interface":
                                   [r"only relevant if illogical input for \c ql and \c qr are provided"]
                                   })
input_physics_oita.setType("double")
input_physics_oita.setRank(1)
input_physics_oita.setDefaultValue(0.0)
input_physics_oita.setStartingIndices([r"0"])
input_physics_oita.setMaximumIndices(["MNvol"])

input_physics_mupftol = Variable("mupftol")
input_physics_mupftol.setDescription({r"accuracy to which \f$\mu\f$ and \f$\Delta\psi_p\f$ are required":
                                      [r"only relevant if constraints on transform, enclosed currents etc. are to be satisfied iteratively, see \c Lconstraint"]
                                      })
input_physics_mupftol.setType("double")
input_physics_mupftol.setDefaultValue(1.0e-16)

input_physics_mupfits = Variable("mupfits")
input_physics_mupfits.setDescription({r"an upper limit on the transform/helicity constraint iterations":
                                      [r"only relevant if constraints on transform, enclosed currents etc. are to be satisfied iteratively, see \c Lconstraint",
                                       r"constraint: \c mupfits > 0"]
                                      })
input_physics_mupfits.setType("int")
input_physics_mupfits.setDefaultValue(8)

input_physics_rpol = Variable("rpol")
input_physics_rpol.setDescription({r"poloidal extent of slab (effective radius)":
                                   [r"only relevant if \c Igeometry==1",
                                    r"poloidal size is \f$L = 2\pi*\f$\c rpol"]
                                   })
input_physics_rpol.setType("double")
input_physics_rpol.setDefaultValue(1.0)
input_physics_rpol.setUnit("m")

input_physics_rtor = Variable("rtor")
input_physics_rtor.setDescription({r"toroidal extent of slab (effective radius)":
                                   [r"only relevant if \c Igeometry==1",
                                    r"toroidal size is \f$L = 2\pi*\f$\c rtor"]
                                   })
input_physics_rtor.setType("double")
input_physics_rtor.setDefaultValue(1.0)
input_physics_rtor.setUnit("m")

input_physics_Lreflect = Variable("Lreflect")
input_physics_Lreflect.setDescription(r"=1 reflect the upper and lower bound in slab, =0 do not reflect")
input_physics_Lreflect.setType("int")
input_physics_Lreflect.setDefaultValue(0)

input_physics_Rac = Variable("Rac")
input_physics_Rac.setDescription(r"    stellarator symmetric coordinate axis; R; cosine")
input_physics_Rac.setType("double")
input_physics_Rac.setRank(1)
input_physics_Rac.setDefaultValue(0.0)
input_physics_Rac.setUnit("m")
input_physics_Rac.setStartingIndices([r"0"])
input_physics_Rac.setMaximumIndices(["MNtor"])

input_physics_Zas = Variable("Zas")
input_physics_Zas.setDescription(r"    stellarator symmetric coordinate axis; Z;   sine")
input_physics_Zas.setType("double")
input_physics_Zas.setRank(1)
input_physics_Zas.setDefaultValue(0.0)
input_physics_Zas.setUnit("m")
input_physics_Zas.setStartingIndices([r"0"])
input_physics_Zas.setMaximumIndices(["MNtor"])

input_physics_Ras = Variable("Ras")
input_physics_Ras.setDescription(r"non-stellarator symmetric coordinate axis; R;   sine")
input_physics_Ras.setType("double")
input_physics_Ras.setRank(1)
input_physics_Ras.setDefaultValue(0.0)
input_physics_Ras.setUnit("m")
input_physics_Ras.setStartingIndices([r"0"])
input_physics_Ras.setMaximumIndices(["MNtor"])

input_physics_Zac = Variable("Zac")
input_physics_Zac.setDescription(r"non-stellarator symmetric coordinate axis; Z; cosine")
input_physics_Zac.setType("double")
input_physics_Zac.setRank(1)
input_physics_Zac.setDefaultValue(0.0)
input_physics_Zac.setUnit("m")
input_physics_Zac.setStartingIndices([r"0"])
input_physics_Zac.setMaximumIndices(["MNtor"])

input_physics_Rbc = Variable("Rbc")
input_physics_Rbc.setDescription(r"    stellarator symmetric boundary components; R; cosine")
input_physics_Rbc.setType("double")
input_physics_Rbc.setRank(2)
input_physics_Rbc.setDefaultValue(0.0)
input_physics_Rbc.setUnit("m")
input_physics_Rbc.setStartingIndices(["-MNtor", "-MMpol"])
input_physics_Rbc.setMaximumIndices(["MNtor", "MMpol"])

input_physics_Zbs = Variable("Zbs")
input_physics_Zbs.setDescription(r"    stellarator symmetric boundary components; Z;   sine")
input_physics_Zbs.setType("double")
input_physics_Zbs.setRank(2)
input_physics_Zbs.setDefaultValue(0.0)
input_physics_Zbs.setUnit("m")
input_physics_Zbs.setStartingIndices(["-MNtor", "-MMpol"])
input_physics_Zbs.setMaximumIndices(["MNtor", "MMpol"])

input_physics_Rbs = Variable("Rbs")
input_physics_Rbs.setDescription(r"non-stellarator symmetric boundary components; R;   sine")
input_physics_Rbs.setType("double")
input_physics_Rbs.setRank(2)
input_physics_Rbs.setDefaultValue(0.0)
input_physics_Rbs.setUnit("m")
input_physics_Rbs.setStartingIndices(["-MNtor", "-MMpol"])
input_physics_Rbs.setMaximumIndices(["MNtor", "MMpol"])

input_physics_Zbc = Variable("Zbc")
input_physics_Zbc.setDescription(r"non-stellarator symmetric boundary components; Z; cosine")
input_physics_Zbc.setType("double")
input_physics_Zbc.setRank(2)
input_physics_Zbc.setDefaultValue(0.0)
input_physics_Zbc.setUnit("m")
input_physics_Zbc.setStartingIndices(["-MNtor", "-MMpol"])
input_physics_Zbc.setMaximumIndices(["MNtor", "MMpol"])

input_physics_Rwc = Variable("Rwc")
input_physics_Rwc.setDescription(r"    stellarator symmetric boundary components of wall; R; cosine")
input_physics_Rwc.setType("double")
input_physics_Rwc.setRank(2)
input_physics_Rwc.setDefaultValue(0.0)
input_physics_Rwc.setUnit("m")
input_physics_Rwc.setStartingIndices(["-MNtor", "-MMpol"])
input_physics_Rwc.setMaximumIndices(["MNtor", "MMpol"])

input_physics_Zws = Variable("Zws")
input_physics_Zws.setDescription(r"    stellarator symmetric boundary components of wall; Z;   sine")
input_physics_Zws.setType("double")
input_physics_Zws.setRank(2)
input_physics_Zws.setDefaultValue(0.0)
input_physics_Zws.setUnit("m")
input_physics_Zws.setStartingIndices(["-MNtor", "-MMpol"])
input_physics_Zws.setMaximumIndices(["MNtor", "MMpol"])

input_physics_Rws = Variable("Rws")
input_physics_Rws.setDescription(r"non-stellarator symmetric boundary components of wall; R;   sine")
input_physics_Rws.setType("double")
input_physics_Rws.setRank(2)
input_physics_Rws.setDefaultValue(0.0)
input_physics_Rws.setUnit("m")
input_physics_Rws.setStartingIndices(["-MNtor", "-MMpol"])
input_physics_Rws.setMaximumIndices(["MNtor", "MMpol"])

input_physics_Zwc = Variable("Zwc")
input_physics_Zwc.setDescription(r"non-stellarator symmetric boundary components of wall; Z; cosine")
input_physics_Zwc.setType("double")
input_physics_Zwc.setRank(2)
input_physics_Zwc.setDefaultValue(0.0)
input_physics_Zwc.setUnit("m")
input_physics_Zwc.setStartingIndices(["-MNtor", "-MMpol"])
input_physics_Zwc.setMaximumIndices(["MNtor", "MMpol"])

input_physics_Vns = Variable("Vns")
input_physics_Vns.setDescription(r"    stellarator symmetric normal field at boundary; vacuum component;   sine")
input_physics_Vns.setType("double")
input_physics_Vns.setRank(2)
input_physics_Vns.setDefaultValue(0.0)
input_physics_Vns.setUnit("T")
input_physics_Vns.setStartingIndices(["-MNtor", "-MMpol"])
input_physics_Vns.setMaximumIndices(["MNtor", "MMpol"])

input_physics_Bns = Variable("Bns")
input_physics_Bns.setDescription(r"    stellarator symmetric normal field at boundary; plasma component;   sine")
input_physics_Bns.setType("double")
input_physics_Bns.setRank(2)
input_physics_Bns.setDefaultValue(0.0)
input_physics_Bns.setUnit("T")
input_physics_Bns.setStartingIndices(["-MNtor", "-MMpol"])
input_physics_Bns.setMaximumIndices(["MNtor", "MMpol"])

input_physics_Vnc = Variable("Vnc")
input_physics_Vnc.setDescription(r"non-stellarator symmetric normal field at boundary; vacuum component; cosine")
input_physics_Vnc.setType("double")
input_physics_Vnc.setRank(2)
input_physics_Vnc.setDefaultValue(0.0)
input_physics_Vnc.setUnit("T")
input_physics_Vnc.setStartingIndices(["-MNtor", "-MMpol"])
input_physics_Vnc.setMaximumIndices(["MNtor", "MMpol"])

input_physics_Bnc = Variable("Bnc")
input_physics_Bnc.setDescription(r"non-stellarator symmetric normal field at boundary; plasma component; cosine")
input_physics_Bnc.setType("double")
input_physics_Bnc.setRank(2)
input_physics_Bnc.setDefaultValue(0.0)
input_physics_Bnc.setUnit("T")
input_physics_Bnc.setStartingIndices(["-MNtor", "-MMpol"])
input_physics_Bnc.setMaximumIndices(["MNtor", "MMpol"])

vars_physicslist = [
        input_physics_Igeometry,
        input_physics_Istellsym,
        input_physics_Lfreebound,
        input_physics_phiedge,
        input_physics_curtor,
        input_physics_curpol,
        input_physics_gamma,
        input_physics_Nfp,
        input_physics_Nvol,
        input_physics_Mpol,
        input_physics_Ntor,
        input_physics_Lrad,
        input_physics_Lconstraint,
        input_physics_tflux,
        input_physics_pflux,
        input_physics_helicity,
        input_physics_pscale,
        input_physics_pressure,
        input_physics_Ladiabatic,
        input_physics_adiabatic,
        input_physics_mu,
        input_physics_pl,
        input_physics_ql,
        input_physics_pr,
        input_physics_qr,
        input_physics_iota,
        input_physics_lp,
        input_physics_lq,
        input_physics_rp,
        input_physics_rq,
        input_physics_oita,
        input_physics_mupftol,
        input_physics_mupfits,
        input_physics_rpol,
        input_physics_rtor,
        input_physics_Rac,
        input_physics_Zas,
        input_physics_Ras,
        input_physics_Zac,
        input_physics_Rbc,
        input_physics_Zbs,
        input_physics_Rbs,
        input_physics_Zbc,
        input_physics_Vns,
        input_physics_Bns,
        input_physics_Vnc,
        input_physics_Bnc
        ]

physicslist = Fortran.Namelist("physicslist")
physicslist.setDescription(r"The namelist \c physicslist controls the geometry, profiles, and numerical resolution.")
physicslist.addVariables(vars_physicslist)

###############################################################################
# numericlist 
###############################################################################

input_numeric_Linitialize = Variable("Linitialize")
input_numeric_Linitialize.setDescription({r"Used to initialize geometry using a regularization / extrapolation method":
                                           [ r"if \c Linitialize = \f$-I\f$ , where \f$I\f$ is a positive integer,"+os.linesep
                                            +r"the geometry of the \f$i=1,N_V-I\f$ surfaces constructed by an extrapolation",
                                             r"if \c Linitialize=0, the geometry of the interior surfaces is provided after the namelists in the input file",
                                             r"if \c Linitialize=1, the interior surfaces will be intialized as \f$R_{l,m,n} = R_{N,m,n} \psi_{t,l}^{m/2}\f$,"+os.linesep
                                            +r"where \f$R_{N,m,n}\f$ is the plasma boundary and \f$\psi_{t,l}\f$ is the given toroidal flux enclosed by the"+os.linesep
                                            +r"\f$l\f$-th interface, normalized to the total enclosed toroidal flux;"+os.linesep
                                            +r"a similar extrapolation is used for \f$Z_{l,m,n}\f$",
                                             r"Note that the Fourier harmonics of the boundary is *always* given by the \c Rbc and \c Zbs"+os.linesep
                                            +r"given in \c physicslist.",
                                             r"if \c Linitialize=2, the interior surfaces *and the plasma boundary* will be intialized"+os.linesep
                                            +r"as \f$R_{l,m,n} = R_{W,m,n} \psi_{t,l}^{m/2}\f$, where \f$R_{W,m,n}\f$ is the computational boundary"+os.linesep
                                            +r"and \f$\psi_{t,l}\f$ is the given toroidal flux enclosed by the \f$l\f$-th interface, normalized to the total enclosed toroidal flux;"+os.linesep
                                            +r"a similar extrapolation is used for \f$Z_{l,m,n}\f$",
                                             r"Note that, for free-boundary calculations, the Fourier harmonics of the computational boundary"+os.linesep
                                            +r"are *always* given by the \c Rwc and \c Zws given in \c physicslist.",
                                             r"if \c Linitialize=1,2 , it is not required to provide the geometry of the interfaces after the namelists"]
                                          })
input_numeric_Linitialize.setType("int")
input_numeric_Linitialize.setDefaultValue(0)

input_numeric_LautoinitBn = Variable("LautoinitBn")
input_numeric_LautoinitBn.setDescription({r"Used to initialize \f$B_{ns}\f$ using an initial fixed-boundary calculation":
                                          [r"only relevant if \c Lfreebound=1",
                                           r"user-supplied \c Bns will only be considered if \c LautoinitBn=0"]
                                          })
input_numeric_LautoinitBn.setType("int")
input_numeric_LautoinitBn.setDefaultValue(1)

input_numeric_Lzerovac = Variable("Lzerovac")
input_numeric_Lzerovac.setDescription({r"Used to adjust vacuum field to cancel plasma field on computational boundary":
                                       [r"only relevant if \c Lfreebound=1"]
                                       })
input_numeric_Lzerovac.setType("int")
input_numeric_Lzerovac.setDefaultValue(0)

input_numeric_Ndiscrete = Variable("Ndiscrete")
input_numeric_Ndiscrete.setDescription({r"resolution of the real space grid on which fast Fourier transforms are performed is given by \c Ndiscrete*Mpol*4":
                                        [r"constraint \c Ndiscrete>0"]
                                        })
input_numeric_Ndiscrete.setType("int")
input_numeric_Ndiscrete.setDefaultValue(2)

input_numeric_Nquad = Variable("Nquad")
input_numeric_Nquad.setDescription({r"Resolution of the Gaussian quadrature":
                                    [ r"The resolution of the Gaussian quadrature, \f$\displaystyle \int \!\! f(s) ds = \sum_k \omega_k f(s_k)\f$,"+os.linesep
                                     +r"in each volume is given by \c Iquad\f$_v\f$",
                                      r"\c Iquad\f$_v\f$ is set in preset()"]
                                    })
input_numeric_Nquad.setType("int")
input_numeric_Nquad.setDefaultValue(-1)

input_numeric_iMpol = Variable("iMpol")
input_numeric_iMpol.setDescription({r"Fourier resolution of straight-fieldline angle on interfaces":
                                    [ r"the rotational-transform on the interfaces is determined by a transformation to the straight-fieldline angle,"+os.linesep
                                     +r"with poloidal resolution given by \c iMpol",
                                      r"if \c iMpol<=0, then \c iMpol = Mpol - iMpol"]
                                    })
input_numeric_iMpol.setType("int")
input_numeric_iMpol.setDefaultValue(-4)

input_numeric_iNtor = Variable("iNtor")
input_numeric_iNtor.setDescription({r"Fourier resolution of straight-fieldline angle on interfaces":
                                    [ r"the rotational-transform on the interfaces is determined by a transformation to the straight-fieldline angle,"+os.linesep
                                     +r"with toroidal resolution given by \c iNtor",
                                      r"if \c iNtor<=0 then \c iNtor = Ntor - iNtor",
                                      r"if \c Ntor=0, then the toroidal resolution of the angle transformation is set \c lNtor=0"]
                                    })
input_numeric_iNtor.setType("int")
input_numeric_iNtor.setDefaultValue(-4)

input_numeric_Lsparse = Variable("Lsparse")
input_numeric_Lsparse.setDescription({r"controls method used to solve for rotational-transform on interfaces":
                                      [ r"if \c Lsparse=0, the transformation to the straight-fieldline angle is computed in Fourier space"+os.linesep
                                       +r"using a dense matrix solver, \c F04AAF",
                                        r"if \c Lsparse=1, the transformation to the straight-fieldline angle is computed in real space"+os.linesep
                                       +r"using a dense matrix solver, \c F04ATF",
                                        r"if \c Lsparse=2, the transformation to the straight-fieldline angle is computed in real space"+os.linesep
                                       +r"using a sparse matrix solver, \c F11DEF",
                                        r"if \c Lsparse=3, the different methods for constructing the straight-fieldline angle are compared"]
                                      })
input_numeric_Lsparse.setType("int")
input_numeric_Lsparse.setDefaultValue(0)

input_numeric_Lsvdiota = Variable("Lsvdiota")
input_numeric_Lsvdiota.setDescription({r"controls method used to solve for rotational-transform on interfaces":
                                       [r"if \c Lsvdiota=0, use standard linear solver to construct straight fieldline angle transformation",
                                        r"if \c Lsvdiota=1, use SVD method to compute rotational-transform"],
                                       r"only relevant if \c Lsparse=0": None
                                       })
input_numeric_Lsvdiota.setType("int")
input_numeric_Lsvdiota.setDefaultValue(0)

input_numeric_imethod = Variable("imethod")
input_numeric_imethod.setDescription({ r"controls iterative solution to sparse matrix"+os.linesep
                                      +r"arising in real-space transformation to the straight-fieldline angle":
                                       [r"if \c imethod=1, the method is \c RGMRES",
                                        r"if \c imethod=2, the method is \c CGS",
                                        r"if \c imethod=3, the method is \c BICGSTAB"],
                                       r"only relevant if \c Lsparse=2; \see tr00ab() for details": None
                                      })
input_numeric_imethod.setType("int")
input_numeric_imethod.setDefaultValue(3)

input_numeric_iorder = Variable("iorder")
input_numeric_iorder.setDescription({r"determines order of finite-difference approximation to the derivatives":
                                     [r"if \c iorder=2, second-order",
                                      r"if \c iorder=4, fourth-order",
                                      r"if \c iorder=6, sixth-order"],
                                     r"controls real-space grid resolution for constructing the straight-fieldline angle": None,
                                     r"only relevant if \c Lsparse>0": None
                                    })
input_numeric_iorder.setType("int")
input_numeric_iorder.setDefaultValue(2)

input_numeric_iprecon = Variable("iprecon")
input_numeric_iprecon.setDescription({ r"controls iterative solution to sparse matrix arising in real-space transformation"+os.linesep
                                      +r"to the straight-fieldline angle":
                                       [r"if \c iprecon=0, the preconditioner is `N'",
                                        r"if \c iprecon=1, the preconditioner is `J'",
                                        r"if \c iprecon=2, the preconditioner is `S'"],
                                       r"only relevant if \c Lsparse=2; \see tr00ab() for details": None
                                      })
input_numeric_iprecon.setType("int")
input_numeric_iprecon.setDefaultValue(0)

input_numeric_iotatol = Variable("iotatol")
input_numeric_iotatol.setDescription({r"tolerance required for iterative construction of straight-fieldline angle":
                                      r"only relevant if \c Lsparse.ge.2"})
input_numeric_iotatol.setType("double")
input_numeric_iotatol.setDefaultValue(-1.0)

input_numeric_Lextrap = Variable("Lextrap")
input_numeric_Lextrap.setDescription(r"geometry of innermost interface is defined by extrapolation")
input_numeric_Lextrap.setType("int")
input_numeric_Lextrap.setDefaultValue(0)

input_numeric_Mregular = Variable("Mregular")
input_numeric_Mregular.setDescription({r"maximum regularization factor":
                                       [r"if \c Mregular>=2, then \c regumm \f$_i\f$ = \c Mregular \f$/ 2 \f$ where \c m \f$_i > \f$ \c Mregular"]
                                       })
input_numeric_Mregular.setType("int")
input_numeric_Mregular.setDefaultValue(-1)

input_numeric_Lrzaxis = Variable("Lrzaxis")
input_numeric_Lrzaxis.setDescription({r"controls the guess of geometry axis in the innermost volume or initialization of interfaces":
                                      [r"if \c iprecon = 1, the centroid is used",
                                       r"if \c iprecon = 2, the Jacobian \f$m=1\f$ harmonic elimination method is used"]
                                      })
input_numeric_Lrzaxis.setType("int")
input_numeric_Lrzaxis.setDefaultValue(1)

input_numeric_Ntoraxis = Variable("Ntoraxis")
input_numeric_Ntoraxis.setDescription( r"the number of \f$n\f$ harmonics used in the Jacobian \f$m=1\f$ harmonic elimination method;"+os.linesep
                                      +r"only relevant if \c Lrzaxis.ge.1 .")
input_numeric_Ntoraxis.setType("int")
input_numeric_Ntoraxis.setDefaultValue(3)

vars_numericlist = [
        input_numeric_Linitialize,
        input_numeric_LautoinitBn,
        input_numeric_Lzerovac,
        input_numeric_Ndiscrete,
        input_numeric_Nquad,
        input_numeric_iMpol,
        input_numeric_iNtor,
        input_numeric_Lsparse,
        input_numeric_Lsvdiota,
        input_numeric_imethod,
        input_numeric_iorder,
        input_numeric_iprecon,
        input_numeric_iotatol,
        input_numeric_Lextrap,
        input_numeric_Mregular
        ]

numericlist = Fortran.Namelist("numericlist")
numericlist.setDescription(r"The namelist \c numericlist controls internal resolution parameters that the user rarely needs to consider.")
numericlist.addVariables(vars_numericlist)

###############################################################################
# locallist 
###############################################################################

input_local_LBeltrami = Variable("LBeltrami")
input_local_LBeltrami.setDescription({r"Control flag for solution of Beltrami equation":
                                      [ r"if \c LBeltrami = 1,3,5 or 7, (SQP) then the Beltrami field in each volume is constructed"+os.linesep
                                       +r"by minimizing the magnetic energy with the constraint of fixed helicity;"+os.linesep
                                       +r"this is achieved by using sequential quadratic programming as provided by \c E04UFF ."+os.linesep
                                       +r"This approach has the benefit (in theory) of robustly constructing minimum energy solutions"+os.linesep
                                       +r"when multiple, i.e. bifurcated, solutions exist.",
                                        r"if \c LBeltrami = 2,3,6 or 7, (Newton) then the Beltrami fields are constructed by employing a standard Newton method"+os.linesep
                                       +r"for locating an extremum of"+os.linesep
                                       +r"\f$F\equiv \int B^2 dv - \mu (\int {\bf A}\cdot{\bf B}dv-{\cal K})\f$,"+os.linesep
                                       +r"where \f$\mu\f$ is treated as an independent degree of freedom similar to the parameters describing the vector potential"+os.linesep
                                       +r"and \f${\cal K}\f$ is the required value of the helicity;"+os.linesep
                                       +r"this is the standard Lagrange multipler approach for locating the constrained minimum;"+os.linesep
                                       +r"this method cannot distinguish saddle-type extrema from minima, and which solution that will be obtained depends on the initial guess",
                                        r"if \c LBeltrami = 4,5,6 or 7, (linear) it is assumed that the Beltrami fields are parameterized by \f$\mu\f$;"+os.linesep
                                       +r"in this case, it is only required to solve \f$\nabla \times {\bf B} = \mu {\bf B}\f$ which reduces to a system of linear equations;"+os.linesep
                                       +r"\f$\mu\f$ may or may not be adjusted iteratively, depending on \c Lconstraint,"+os.linesep
                                       +r"to satisfy either rotational-transform or helicity constraints",
                                       {r"for flexibility and comparison, each of the above methods can be employed; for example:":
                                        [r"if \c LBeltrami=1 , only the SQP    method will be employed",
                                         r"if \c LBeltrami=2 , only the Newton method will be employed",
                                         r"if \c LBeltrami=4 , only the linear method will be employed",
                                         r"if \c LBeltrami=3 , the SQP and the Newton method are used",
                                         r"if \c LBeltrami=5 , the SQP and the linear method are used",
                                         r"if \c LBeltrami=6 , the Newton and the linear method are used",
                                         r"if \c LBeltrami=7 , all three methods will be employed"]
                                        }]
                                      })
input_local_LBeltrami.setType("int")
input_local_LBeltrami.setDefaultValue(4)

input_local_Linitgues = Variable("Linitgues")
input_local_Linitgues.setDescription({r"controls how initial guess for Beltrami field is constructed":
                                      [ r"only relevant for routines that require an initial guess for the Beltrami fields, such as the SQP and Newton methods,"+os.linesep
                                       +r"or the sparse linear solver",
                                        r"if \c Linitgues=0, the initial guess for the Beltrami field is trivial",
                                        r"if \c Linitgues=1, the initial guess for the Beltrami field is an integrable approximation",
                                        r"if \c Linitgues=2, the initial guess for the Beltrami field is read from file",
                                        r"if \c Linitgues=3, the initial guess for the Beltrami field will be randomized with the maximum \c maxrndgues"]
                                      })
input_local_Linitgues.setType("int")
input_local_Linitgues.setDefaultValue(1)

input_local_Lposdef = Variable("Lposdef")
input_local_Lposdef.setDescription(r"redundant")
input_local_Lposdef.setType("int")
input_local_Lposdef.setDefaultValue(0)

input_local_maxrndgues = Variable("maxrndgues")
input_local_maxrndgues.setDescription(r"the maximum random number of the Beltrami field if \c Linitgues = 3")
input_local_maxrndgues.setType("double")
input_local_maxrndgues.setDefaultValue(1.0)

input_local_Lmatsolver = Variable("Lmatsolver")
input_local_Lmatsolver.setDescription(r"1 for LU factorization, 2 for GMRES, 3 for GMRES matrix-free")
input_local_Lmatsolver.setType("int")
input_local_Lmatsolver.setDefaultValue(3)

input_local_NiterGMRES = Variable("NiterGMRES")
input_local_NiterGMRES.setDescription(r"number of max iteration for GMRES")
input_local_NiterGMRES.setType("int")
input_local_NiterGMRES.setDefaultValue(200)

input_local_epsGMRES = Variable("epsGMRES")
input_local_epsGMRES.setDescription(r"the precision of GMRES")
input_local_epsGMRES.setType("double")
input_local_epsGMRES.setDefaultValue(1.0e-14)

input_local_LGMRESprec = Variable("LGMRESprec")
input_local_LGMRESprec.setDescription(r"type of preconditioner for GMRES, 1 for ILU sparse matrix")
input_local_LGMRESprec.setType("int")
input_local_LGMRESprec.setDefaultValue(1)

input_local_epsILU = Variable("epsILU")
input_local_epsILU.setDescription(r"the precision of incomplete LU factorization for preconditioning")
input_local_epsILU.setType("double")
input_local_epsILU.setDefaultValue(1.0e-12)

vars_locallist = [
        input_local_LBeltrami,
        input_local_Linitgues,
        input_local_Lposdef,
        input_local_maxrndgues
        ]

locallist = Fortran.Namelist("locallist")
locallist.setDescription({r"The namelist \c locallist controls the construction of the Beltrami fields in each volume.":
                          [ r"The transformation to straight-fieldline coordinates is singular when the rotational-transform of the interfaces is rational;"+os.linesep
                           +r"however, the rotational-transform is still well defined."]
                          })
locallist.addVariables(vars_locallist)

###############################################################################
# globallist 
###############################################################################

input_global_Lfindzero = Variable("Lfindzero")
input_global_Lfindzero.setDescription({r"use Newton methods to find zero of force-balance, which is computed by dforce()":
                                       [ r"if \c Lfindzero=0 , then dforce() is called once"+os.linesep
                                        +r"to compute the Beltrami fields consistent with the given geometry and constraints",
                                         r"if \c Lfindzero=1 , then call \c C05NDF (uses   function values only), which iteratively calls dforce()",
                                         r"if \c Lfindzero=2 , then call \c C05PDF (uses derivative information), which iteratively calls dforce()"]
                                       })
input_global_Lfindzero.setType("int")
input_global_Lfindzero.setDefaultValue(0)

input_global_escale = Variable("escale")
input_global_escale.setDescription({r"controls the weight factor, \c BBweight, in the force-imbalance harmonics":
                                    [r"\c BBweight(i) \f$\displaystyle \equiv \texttt{opsilon} \times \exp\left[-\texttt{escale} \times (m_i^2+n_i^2) \right]\f$",
                                     r"defined in preset() ; used in dforce()",
                                     r"\sa Eqn.\f$(\ref{eq:forcebalancemn_global})\f$"]
                                    })
input_global_escale.setType("double")
input_global_escale.setDefaultValue(0.0)

input_global_opsilon = Variable("opsilon")
input_global_opsilon.setDescription({r"weighting of force-imbalance":
                                     [r"used in dforce(); \sa Eqn.\f$(\ref{eq:forcebalancemn_global})\f$"]
                                     })
input_global_opsilon.setType("double")
input_global_opsilon.setDefaultValue(1.0)

input_global_pcondense = Variable("pcondense")
input_global_pcondense.setDescription({r"spectral condensation parameter":
                                       [ r"used in preset() to define \c mmpp(i) \f$\equiv m_i^p\f$, where \f$p\equiv \f$ \c pcondense",
                                         r"the angle freedom is exploited to minimize \f$\displaystyle \texttt{epsilon} \sum_{i} m_i^p (R_{i}^2+Z_{i}^2)\f$"+os.linesep
                                        +r"with respect to tangential variations in the interface geometry",
                                         r"\sa Eqn.\f$(\ref{eq:spectralbalancemn_global})\f$"]
                                       })
input_global_pcondense.setType("double")
input_global_pcondense.setDefaultValue(2.0)

input_global_epsilon = Variable("epsilon")
input_global_epsilon.setDescription({r"weighting of spectral-width constraint":
                                     [r"used in dforce(); \sa Eqn.\f$(\ref{eq:spectralbalancemn_global})\f$"]
                                     })
input_global_epsilon.setType("double")
input_global_epsilon.setDefaultValue(0.0)

input_global_wpoloidal = Variable("wpoloidal")
input_global_wpoloidal.setDescription([r"\"star-like\" poloidal angle constraint radial exponential factor",
                                       r"used in preset() to construct \c sweight"])
input_global_wpoloidal.setType("double")
input_global_wpoloidal.setDefaultValue(1.0)

input_global_upsilon = Variable("upsilon")
input_global_upsilon.setDescription([r"weighting of \"star-like\" poloidal angle constraint",
                                     r"used in preset() to construct \c sweight"])
input_global_upsilon.setType("double")
input_global_upsilon.setDefaultValue(1.0)

input_global_forcetol = Variable("forcetol")
input_global_forcetol.setDescription({r"required tolerance in force-balance error; only used as an initial check":
                                      [ r"if the initially supplied interfaces are consistent with force-balance to within \c forcetol"+os.linesep
                                       +r"then the geometry of the interfaces is not altered",
                                        r"if not, then the geometry of the interfaces is changed in order to bring the configuration into force balance"+os.linesep
                                       +r"so that the geometry of interfaces is within \c c05xtol, defined below, of the true solution",
                                        r"to force execution of either \c C05NDF or \c C05PDF, regardless of the initial force imbalance,"+os.linesep
                                       +r"set \c forcetol<0"]
                                      })
input_global_forcetol.setType("double")
input_global_forcetol.setDefaultValue(1.0e-10)

input_global_c05xmax = Variable("c05xmax")
input_global_c05xmax.setDescription(r"required tolerance in position, \f${\bf x} \equiv \{ R_{i,v}, Z_{i,v}\}\f$")
input_global_c05xmax.setType("double")
input_global_c05xmax.setDefaultValue(1.0e-6)

input_global_c05xtol = Variable("c05xtol")
input_global_c05xtol.setDescription({r"required tolerance in position, \f${\bf x} \equiv \{ R_{i,v}, Z_{i,v}\}\f$":
                                     [r"used by both \c C05NDF and \c C05PDF; see the NAG documents for further details on how the error is defined",
                                      r"constraint \c c05xtol>0.0"]
                                     })
input_global_c05xtol.setType("double")
input_global_c05xtol.setDefaultValue(1.0e-12)

input_global_c05factor = Variable("c05factor")
input_global_c05factor.setDescription({r"used to control initial step size in \c C05NDF and \c C05PDF":
                                       [r"constraint \c c05factor>0.0",
                                        r"only relevant if \c Lfindzero>0"]
                                       })
input_global_c05factor.setType("double")
input_global_c05factor.setDefaultValue(1.0e-2)

input_global_LreadGF = Variable("LreadGF")
input_global_LreadGF.setDescription({r"read \f$\nabla_{\bf x} {\bf F}\f$ from file \c ext.GF ":
                                     [r"only used if \c Lfindzero=2",
                                      r"only used in newton()"]
                                     })
input_global_LreadGF.setType("boolean")
input_global_LreadGF.setDefaultValue(True)

input_global_mfreeits = Variable("mfreeits")
input_global_mfreeits.setDescription({r"maximum allowed free-boundary iterations":
                                      [r"only used if \c Lfreebound=1",
                                       r"only used in xspech()"]
                                      })
input_global_mfreeits.setType("int")
input_global_mfreeits.setDefaultValue(0)

input_global_bnstol = Variable("bnstol")
input_global_bnstol.setDescription(r"redundant")
input_global_bnstol.setType("double")
input_global_bnstol.setDefaultValue(1.0e-6)

input_global_bnsblend = Variable("bnsblend")
input_global_bnsblend.setDescription(r"redundant")
input_global_bnsblend.setType("double")
input_global_bnsblend.setDefaultValue(0.666)

input_global_gBntol = Variable("gBntol")
input_global_gBntol.setDescription({r"required tolerance in free-boundary iterations":
                                    [r"only used if \c Lfreebound=1",
                                     r"only used in xspech()"]
                                    })
input_global_gBntol.setType("double")
input_global_gBntol.setDefaultValue(1.0e-6)

input_global_gBnbld = Variable("gBnbld")
input_global_gBnbld.setDescription({r"normal blend":
                                    [ r"The \"new\" magnetic field at the computational boundary produced by the plasma currents is updated using a Picard scheme:"+os.linesep
                                     +r"\f{eqnarray}{ ({\bf B}\cdot{\bf n})^{j+1} =    \texttt{gBnbld}  \times ({\bf B}\cdot{\bf n})^{j}"+os.linesep
                                     +r"                                          + (1-\texttt{gBnbld}) \times ({\bf B}\cdot{\bf n})^{*},"+os.linesep
                                     +r"\f}"+os.linesep
                                     +r"where \f$j\f$ labels free-boundary iterations, and \f$({\bf B}\cdot{\bf n})^{*}\f$ is computed by virtual casing.",
                                      r"only used if \c Lfreebound=1",
                                      r"only used in xspech()"]
                                    })
input_global_gBnbld.setType("double")
input_global_gBnbld.setDefaultValue(0.666)

input_global_vcasingeps = Variable("vcasingeps")
input_global_vcasingeps.setDescription(r"regularization of Biot-Savart; see bnorml(), casing()")
input_global_vcasingeps.setType("double")
input_global_vcasingeps.setDefaultValue(1.0e-12)

input_global_vcasingtol = Variable("vcasingtol")
input_global_vcasingtol.setDescription(r"accuracy on virtual casing integral; see bnorml(), casing()")
input_global_vcasingtol.setType("double")
input_global_vcasingtol.setDefaultValue(1.0e-8)

input_global_vcasingits = Variable("vcasingits")
input_global_vcasingits.setDescription(r"minimum number of calls to adaptive virtual casing routine; see casing()")
input_global_vcasingits.setType("int")
input_global_vcasingits.setDefaultValue(8)

input_global_vcasingper = Variable("vcasingper")
input_global_vcasingper.setDescription(r"periods of integragion  in adaptive virtual casing routine; see casing()")
input_global_vcasingper.setType("int")
input_global_vcasingper.setDefaultValue(1)

input_global_mcasingcal = Variable("mcasingcal")
input_global_mcasingcal.setDescription(r'minimum number of calls to adaptive virtual casing routine; see casing() redundant')
input_global_mcasingcal.setType("int")
input_global_mcasingcal.setDefaultValue(8)

vars_globallist = [
        input_global_Lfindzero,
        input_global_escale,
        input_global_opsilon,
        input_global_pcondense,
        input_global_epsilon,
        input_global_wpoloidal,
        input_global_upsilon,
        input_global_forcetol,
        input_global_c05xmax,
        input_global_c05xtol,
        input_global_c05factor,
        input_global_LreadGF,
        input_global_mfreeits,
        input_global_bnstol,
        input_global_bnsblend,
        input_global_gBntol,
        input_global_gBnbld,
        input_global_vcasingeps,
        input_global_vcasingtol,
        input_global_vcasingits,
        input_global_vcasingper,
        input_global_mcasingcal
        ]

globallist = Fortran.Namelist("globallist")
globallist.setDescription({r"The namelist \c globallist controls the search for global force-balance.": None,
                           r"Comments:":
                           [ r'The "force" vector, \f${\bf F}\f$, which is constructed in dforce(), is a combination of pressure-imbalance Fourier harmonics,'+os.linesep
                            +r"\f{eqnarray}{ F_{i,v} \equiv [[ p+B^2/2 ]]_{i,v} \times \exp\left[-\texttt{escale}(m_i^2+n_i^2) \right] \times \texttt{opsilon},"+os.linesep
                            +r"\label{eq:forcebalancemn_global} \f}"+os.linesep
                            +r'and spectral-condensation constraints, \f$I_{i,v}\f$, and the "star-like" angle constraints, \f$S_{i,v,}\f$, (see lforce() for details)'+os.linesep
                            +r"\f{eqnarray}{ F_{i,v} \equiv \texttt{epsilon} \times I_{i,v}"+os.linesep
                            +r"                           + \texttt{upsilon} \times \left( \psi_v^\omega S_{i,v,1} - \psi_{v+1}^\omega S_{i,v+1,0} \right),"+os.linesep
                            +r"\label{eq:spectralbalancemn_global} \f}"+os.linesep
                            +r"where \f$\psi_v\equiv\f$ normalized toroidal flux, \c tflux, and \f$\omega\equiv\f$ \c wpoloidal."]
                           })
globallist.addVariables(vars_globallist)

###############################################################################
# diagnosticslist 
###############################################################################

input_diagnostics_odetol = Variable("odetol")
input_diagnostics_odetol.setDescription(r"o.d.e. integration tolerance for all field line tracing routines")
input_diagnostics_odetol.setType("double")
input_diagnostics_odetol.setDefaultValue(1.0e-7)

input_diagnostics_absreq = Variable("absreq")
input_diagnostics_absreq.setDescription(r"redundant")
input_diagnostics_absreq.setType("double")
input_diagnostics_absreq.setDefaultValue(1.0e-8)

input_diagnostics_relreq = Variable("relreq")
input_diagnostics_relreq.setDescription(r"redundant")
input_diagnostics_relreq.setType("double")
input_diagnostics_relreq.setDefaultValue(1.0e-8)

input_diagnostics_absacc = Variable("absacc")
input_diagnostics_absacc.setDescription(r"redundant")
input_diagnostics_absacc.setType("double")
input_diagnostics_absacc.setDefaultValue(1.0e-4)

input_diagnostics_epsr = Variable("epsr")
input_diagnostics_epsr.setDescription(r"redundant")
input_diagnostics_epsr.setType("double")
input_diagnostics_epsr.setDefaultValue(1.0e-8)

input_diagnostics_nPpts = Variable("nPpts")
input_diagnostics_nPpts.setDescription({ r"number of toroidal transits used (per trajectory) in following field lines"+os.linesep
                                        +r"for constructing Poincaré plots":
                                         [r"if \c nPpts<1, no Poincaré plot is constructed"]
                                         })
input_diagnostics_nPpts.setType("int")
input_diagnostics_nPpts.setDefaultValue(0)

input_diagnostics_Ppts = Variable("Ppts")
input_diagnostics_Ppts.setDescription(r"stands for Poincare plot theta start. Chose at which angle (normalized over \f$\pi\f$) the Poincare field-line tracing start.")
input_diagnostics_Ppts.setType("double")
input_diagnostics_Ppts.setDefaultValue(0.0)

input_diagnostics_nPtrj = Variable("nPtrj")
input_diagnostics_nPtrj.setDescription({r"number of trajectories in each annulus to be followed in constructing Poincaré plot":
                                        [ r"if \c nPtrj(l)<0, then \c nPtrj(l) = Ni(l),"+os.linesep
                                         +r"where \c Ni(l) is the grid resolution used to construct the Beltrami field in volume \f$l\f$"]
                                        })
input_diagnostics_nPtrj.setType("int")
input_diagnostics_nPtrj.setRank(1)
input_diagnostics_nPtrj.setDefaultValue(-1)
input_diagnostics_nPtrj.setMaximumIndices(["MNvol+1"])

input_diagnostics_LHevalues = Variable("LHevalues")
input_diagnostics_LHevalues.setDescription(r"to compute eigenvalues of \f$\nabla {\bf F}\f$")
input_diagnostics_LHevalues.setType("boolean")
input_diagnostics_LHevalues.setDefaultValue(False)

input_diagnostics_LHevectors = Variable("LHevectors")
input_diagnostics_LHevectors.setDescription(r"to compute eigenvectors (and also eigenvalues) of \f$\nabla {\bf F}\f$")
input_diagnostics_LHevectors.setType("boolean")
input_diagnostics_LHevectors.setDefaultValue(False)

input_diagnostics_LHmatrix = Variable("LHmatrix")
input_diagnostics_LHmatrix.setDescription(r"to compute and write to file the elements of \f$\nabla {\bf F}\f$")
input_diagnostics_LHmatrix.setType("boolean")
input_diagnostics_LHmatrix.setDefaultValue(False)

input_diagnostics_Lperturbed = Variable("Lperturbed")
input_diagnostics_Lperturbed.setDescription(r"to compute linear, perturbed equilibrium")
input_diagnostics_Lperturbed.setType("int")
input_diagnostics_Lperturbed.setDefaultValue(0)

input_diagnostics_dpp = Variable("dpp")
input_diagnostics_dpp.setDescription(r"perturbed harmonic")
input_diagnostics_dpp.setType("int")
input_diagnostics_dpp.setDefaultValue(-1)

input_diagnostics_dqq = Variable("dqq")
input_diagnostics_dqq.setDescription(r"perturbed harmonic")
input_diagnostics_dqq.setType("int")
input_diagnostics_dqq.setDefaultValue(-1)

input_diagnostics_Lerrortype = Variable("Lerrortype")
input_diagnostics_Lerrortype.setDescription(r"the type of error output for Lcheck=1")
input_diagnostics_Lerrortype.setType("int")
input_diagnostics_Lerrortype.setDefaultValue(0)

input_diagnostics_Ngrid = Variable("Ngrid")
input_diagnostics_Ngrid.setDescription(r"the number of points to output in the grid, -1 for Lrad(vvol)")
input_diagnostics_Ngrid.setType("int")
input_diagnostics_Ngrid.setDefaultValue(-1)

input_diagnostics_dRZ = Variable("dRZ")
input_diagnostics_dRZ.setDescription(r"difference in geometry for finite difference estimate (debug only)")
input_diagnostics_dRZ.setType("double")
input_diagnostics_dRZ.setDefaultValue(1.0e-5)

input_diagnostics_Lcheck = Variable("Lcheck")
input_diagnostics_Lcheck.setDescription({r"implement various checks":
                                         [ r"if \c Lcheck = 0, no additional check on the calculation is performed",
                                           r"if \c Lcheck = 1, the error in the current, i.e. \f$\nabla\times{\bf B}-\mu{\bf B}\f$ is computed as a post-diagnostic",
                                          {r"if \c Lcheck = 2, the analytic derivatives of the interface transform w.r.t."+os.linesep
                                          +r"the helicity multiplier, \f$\mu\f$, and the enclosed poloidal flux, \f$\Delta\psi_p\f$, are compared to a finite-difference estimate":
                                          [ r"only if \c Lconstraint=1",
                                            r"only for \c dspec executable, i.e. must compile with \c DFLAGS=\"-D DEBUG\""]},
                                          {r"if \c Lcheck = 3, the analytic derivatives of the volume w.r.t. interface Fourier harmonic"+os.linesep
                                          +r"is compared to a finite-difference estimate":
                                          [ r"must set \c Lfindzero=2",
                                            r"set \c forcetol sufficiently small and set \c LreadGF=F,"+os.linesep
                                           +r"so that the matrix of second derivatives is calculated",
                                            r"only for \c dspec executable, i.e. must compile with \c DFLAGS=\"-D DEBUG\""]},
                                          {r"if \c Lcheck = 4, the analytic calculation of the derivatives of the magnetic field, \f$B^2\f$, at the interfaces"+os.linesep
                                          +r"is compared to a finite-difference estimate":
                                          [ r"must set \c Lfindzero=2",
                                            r"set \c forcetol sufficiently small",
                                            r"set \c LreadGF=F",
                                            r"only for \c dspec executable, i.e. must compile with \c DFLAGS=\"-D DEBUG\""]},
                                           r"if \c Lcheck = 5, the analytic calculation of the matrix of the derivatives of the force imbalance"+os.linesep
                                          +r"is compared to a finite-difference estimate",
                                          {r"if \c Lcheck = 6, the virtual casing calculation is compared to \c xdiagno (Lazerson 2013 \cite y2013_lazerson)":
                                           [ r"the input file for \c xdiagno is written by bnorml()",
                                             r"this provides the Cartesian coordinates on the computational boundary where the virtual casing routine casing()"+os.linesep
                                            +r"computes the magnetic field, with the values of the magnetic field being written to the screen for comparison",
                                             r"must set \c Freebound=1, \c Lfindzero>0, \c mfreeits!=0",
                                             r"\c xdiagno must be executed manually"]}
                                          ]
                                         })
input_diagnostics_Lcheck.setType("int")
input_diagnostics_Lcheck.setDefaultValue(0)

input_diagnostics_Ltiming = Variable("Ltiming")
input_diagnostics_Ltiming.setDescription(r"to check timing")
input_diagnostics_Ltiming.setType("boolean")
input_diagnostics_Ltiming.setDefaultValue(False)

input_diagnostics_fudge = Variable("fudge")
input_diagnostics_fudge.setDescription(r"redundant")
input_diagnostics_fudge.setType("double")
input_diagnostics_fudge.setDefaultValue(1.0)

input_diagnostics_scaling = Variable("scaling")
input_diagnostics_scaling.setDescription(r"redundant")
input_diagnostics_scaling.setType("double")
input_diagnostics_scaling.setDefaultValue(1.0)

vars_diagnosticslist = [
        input_diagnostics_odetol,
        input_diagnostics_absreq,
        input_diagnostics_relreq,
        input_diagnostics_absacc,
        input_diagnostics_epsr,
        input_diagnostics_nPpts,
        input_diagnostics_nPtrj,
        input_diagnostics_LHevalues,
        input_diagnostics_LHevectors,
        input_diagnostics_LHmatrix,
        input_diagnostics_Lperturbed,
        input_diagnostics_dpp,
        input_diagnostics_dqq,
        input_diagnostics_Lcheck,
        input_diagnostics_Ltiming,
        input_diagnostics_fudge,
        input_diagnostics_scaling
        ]

diagnosticslist = Fortran.Namelist("diagnosticslist")
diagnosticslist.setDescription(r"The namelist \c diagnosticslist controls post-processor diagnostics, such as Poincaré  plot resolution, etc.")
diagnosticslist.addVariables(vars_diagnosticslist)

###############################################################################
# screenlist --> not in output file
###############################################################################


###############################################################################
# initial guess for geometry of ideal interfaces
###############################################################################

# TODO











###############################################################################
# final list of input namelists for SPEC
###############################################################################
input_namelists = [physicslist,
                   numericlist,
                   locallist,
                   globallist,
                   diagnosticslist]

#-----------------------------------------------------------------------------#

###############################################################################
# output quantities
###############################################################################

iterations_nDcalls = Variable("nDcalls")
iterations_nDcalls.setDescription(r"number of calls to something (?)")
iterations_nDcalls.setType("int")

iterations_Energy = Variable("Energy")
iterations_Energy.setDescription(r"MRxMHD energy in the full plasma")
iterations_Energy.setType("double")

iterations_ForceErr = Variable("ForceErr")
iterations_ForceErr.setDescription(r"residual force on the ideal interfaces in the plasma")
iterations_ForceErr.setType("double")

iterations_iRbc = Variable("iRbc")
iterations_iRbc.setDescription(r"stellarator symmetric interface components; R; cosine")
iterations_iRbc.setType("double")
iterations_iRbc.setUnit("m")
iterations_iRbc.setRank(3)
iterations_iRbc.setStartingIndices(["-MNtor", "-MMpol", "1"])
iterations_iRbc.setMaximumIndices(["MNtor", "MMpol", "MNvol+1"])

iterations_iZbs = Variable("iZbs")
iterations_iZbs.setDescription(r"stellarator symmetric interface components; Z;   sine")
iterations_iZbs.setType("double")
iterations_iZbs.setUnit("m")
iterations_iZbs.setRank(3)
iterations_iZbs.setStartingIndices(["-MNtor", "-MMpol", "1"])
iterations_iZbs.setMaximumIndices(["MNtor", "MMpol", "MNvol+1"])

iterations_iRbs = Variable("iRbs")
iterations_iRbs.setDescription(r"non-stellarator symmetric interface components; R;   sine")
iterations_iRbs.setType("double")
iterations_iRbs.setUnit("m")
iterations_iRbs.setRank(3)
iterations_iRbs.setStartingIndices(["-MNtor", "-MMpol", "1"])
iterations_iRbs.setMaximumIndices(["MNtor", "MMpol", "MNvol+1"])

iterations_iZbc = Variable("iZbc")
iterations_iZbc.setDescription(r"non-stellarator symmetric interface components; Z; cosine")
iterations_iZbc.setType("double")
iterations_iZbc.setUnit("m")
iterations_iZbc.setRank(3)
iterations_iZbc.setStartingIndices(["-MNtor", "-MMpol", "1"])
iterations_iZbc.setMaximumIndices(["MNtor", "MMpol", "MNvol+1"])

# this one is a little bit special:
# 1-dim array of a compound datatype, unlimited length (to allow convergence until eternity)
iterations = Variable("iterations")
iterations.setDescription(r"convergence log of force, energy and interface geometry")
iterations.setType([iterations_nDcalls, iterations_Energy, iterations_ForceErr,
                    iterations_iRbc, iterations_iZbs, iterations_iRbs, iterations_iZbc])
iterations.setRank(1)
iterations.setMaximumIndices(["UNLIMITED"])


grid_Rij = Variable("Rij")
grid_Rij.setDescription(r"R positions at which the magnetic field is evaluated")
grid_Rij.setType("double")
grid_Rij.setUnit("m")
grid_Rij.setRank(3)
grid_Rij.setMaximumIndices(["Mvol", "Ngrid_local", "Ntz"])

grid_Zij = Variable("Zij")
grid_Zij.setDescription(r"Z positions at which the magnetic field is evaluated")
grid_Zij.setType("double")
grid_Zij.setUnit("m")
grid_Zij.setRank(3)
grid_Zij.setMaximumIndices(["Mvol", "Ngrid_local", "Ntz"])

grid_sg = Variable("sg")
grid_sg.setDescription(r"jacobian at positions at which the magnetic field is evaluated")
grid_sg.setType("double")
grid_sg.setRank(3)
grid_sg.setMaximumIndices(["Mvol", "Ngrid_local", "Ntz"])

grid_BR = Variable("BR")
grid_BR.setDescription(r"cylindrical R component of magnetic field")
grid_BR.setType("double")
grid_BR.setUnit("T")
grid_BR.setRank(3)
grid_BR.setMaximumIndices(["Mvol", "Ngrid_local", "Ntz"])

grid_Bp = Variable("Bp")
grid_Bp.setDescription(r"cylindrical phi component of magnetic field")
grid_Bp.setType("double")
grid_Bp.setUnit("T")
grid_Bp.setRank(3)
grid_Bp.setMaximumIndices(["Mvol", "Ngrid_local", "Ntz"])

grid_BZ = Variable("BZ")
grid_BZ.setDescription(r"Z component of magnetic field")
grid_BZ.setType("double")
grid_BZ.setUnit("T")
grid_BZ.setRank(3)
grid_BZ.setMaximumIndices(["Mvol", "Ngrid_local", "Ntz"])

# grid group
vars_grid = [grid_Rij,
        grid_Zij,
        grid_sg,
        grid_BR,
        grid_Bp,
        grid_BZ]

poincare_t = Variable("t")
poincare_t.setDescription(r"theta positions of field-line tracing result")
poincare_t.setType("double")
poincare_t.setRank(3)
poincare_t.setMaximumIndices(["Nz", "nPpts", "numTrajTotal"])

poincare_s = Variable("s")
poincare_s.setDescription(r"s positions of field-line tracing result")
poincare_s.setType("double")
poincare_s.setRank(3)
poincare_s.setMaximumIndices(["Nz", "nPpts", "numTrajTotal"])

poincare_R = Variable("R")
poincare_R.setDescription(r"R positions of field-line tracing result")
poincare_R.setType("double")
poincare_R.setUnit("m")
poincare_R.setRank(3)
poincare_R.setMaximumIndices(["Nz", "nPpts", "numTrajTotal"])

poincare_Z = Variable("Z")
poincare_Z.setDescription(r"Z positions of field-line tracing result")
poincare_Z.setType("double")
poincare_Z.setUnit("m")
poincare_Z.setRank(3)
poincare_Z.setMaximumIndices(["Nz", "nPpts", "numTrajTotal"])

poincare_success = Variable("success")
poincare_success.setDescription(r"flag to indicate if a given trajectory was successfully followed")
poincare_success.setType("boolean")
poincare_success.setRank(1)
poincare_success.setMaximumIndices(["numTrajTotal"])

poincare_diotadxup = Variable("diotadxup")
poincare_diotadxup.setDescription(r"measured rotational transform on inner/outer interfaces for each volume; d(transform)/dx")
poincare_diotadxup.setType("double")
poincare_diotadxup.setRank(2)
poincare_diotadxup.setMaximumIndices(["2", "Mvol"])

poincare_fiota = Variable("fiota")
poincare_fiota.setDescription(r"rotational transform from field-line tracing")
poincare_fiota.setType("double")
poincare_fiota.setRank(2)
poincare_fiota.setMaximumIndices(["numTrajTotal", "2"])

vars_poincare = [poincare_t,
                 poincare_s,
                 poincare_R,
                 poincare_Z,
                 poincare_success,
                 poincare_diotadxup,
                 poincare_fiota]

vector_potential_Ate = Variable("Ate")
vector_potential_Ate.setDescription(r"theta component of magnetic vector potential cosine Fourier harmonics; stellarator-symmetric")
vector_potential_Ate.setType("double")
vector_potential_Ate.setUnit("Tm")
vector_potential_Ate.setRank(3)
vector_potential_Ate.setMaximumIndices(["Mvol", "Lrad", "mn"])

vector_potential_Aze = Variable("Aze")
vector_potential_Aze.setDescription(r"zeta component of magnetic vector potential cosine Fourier harmonics; stellarator-symmetric")
vector_potential_Aze.setType("double")
vector_potential_Aze.setUnit("Tm")
vector_potential_Aze.setRank(3)
vector_potential_Aze.setMaximumIndices(["Mvol", "Lrad", "mn"])

vector_potential_Ato = Variable("Ato")
vector_potential_Ato.setDescription(r"theta component of magnetic vector potential sine Fourier harmonics; non-stellarator-symmetric")
vector_potential_Ato.setType("double")
vector_potential_Ato.setUnit("Tm")
vector_potential_Ato.setRank(3)
vector_potential_Ato.setMaximumIndices(["Mvol", "Lrad", "mn"])

vector_potential_Azo = Variable("Azo")
vector_potential_Azo.setDescription(r"zeta component of magnetic vector potential sine Fourier harmonics; non-stellarator-symmetric")
vector_potential_Azo.setType("double")
vector_potential_Azo.setUnit("Tm")
vector_potential_Azo.setRank(3)
vector_potential_Azo.setMaximumIndices(["Mvol", "Lrad", "mn"])

vars_vector_potential = [vector_potential_Ate,
                         vector_potential_Aze,
                         vector_potential_Ato,
                         vector_potential_Azo]

output_mn = Variable("mn")
output_mn.setDescription(r"number of combined polodial and toroidal Fourier harmonics")
output_mn.setType("int")

output_im = Variable("im")
output_im.setDescription(r"poloidal mode number array")
output_im.setType("int")

output_in = Variable("in")
output_in.setDescription(r"toroidal mode number array")
output_in.setType("int")

output_Mvol = Variable("Mvol")
output_Mvol.setDescription(r"number of interfaces")
output_Mvol.setType("int")

output_Rbc = Variable("Rbc")
output_Rbc.setDescription(r"stellarator symmetric boundary components; R; cosine")
output_Rbc.setType("double")
output_Rbc.setRank(2)
output_Rbc.setStartingIndices(["1", "0"])
output_Rbc.setMaximumIndices(["mn", "Mvol"])

output_Zbs = Variable("Zbs")
output_Zbs.setDescription(r"stellarator symmetric boundary components; Z; sine")
output_Zbs.setType("double")
output_Zbs.setRank(2)
output_Zbs.setStartingIndices(["1", "0"])
output_Zbs.setMaximumIndices(["mn", "Mvol"])

output_Rbs = Variable("Rbs")
output_Rbs.setDescription(r"non-stellarator symmetric boundary components; R; sine")
output_Rbs.setType("double")
output_Rbs.setRank(2)
output_Rbs.setStartingIndices(["1", "0"])
output_Rbs.setMaximumIndices(["mn", "Mvol"])

output_Zbc = Variable("Zbc")
output_Zbc.setDescription(r"non-stellarator symmetric boundary components; Z; cosine")
output_Zbc.setType("double")
output_Zbc.setRank(2)
output_Zbc.setStartingIndices(["1", "0"])
output_Zbc.setMaximumIndices(["mn", "Mvol"])

output_ForceErr = Variable("ForceErr")
output_ForceErr.setDescription(r"residual force on the ideal interfaces in the plasma")
output_ForceErr.setType("double")

output_Ivolume = Variable("Ivolume")
output_Ivolume.setDescription(r"Volume current at output (parallel, externally induced)")
output_Ivolume.setType("double")
output_Ivolume.setRank(1)
output_Ivolume.setMaximumIndices(["Mvol"])

output_IPDt = Variable("IPDt")
output_IPDt.setDescription(r"Surface current at output")
output_IPDt.setType("double")
output_IPDt.setRank(1)
output_IPDt.setMaximumIndices(["Mvol"])

output_adiabatic = Variable("adiabatic")
output_adiabatic.setDescription(r"adiabatic constants in each volume")
output_adiabatic.setType("double")
output_adiabatic.setRank(1)
output_adiabatic.setMaximumIndices(["Nvol"])

output_helicity = Variable("helicity")
output_helicity.setDescription(r"the computed values of \f${\cal K} \equiv \int {\bf A}\cdot{\bf B}\;dv\f$")
output_helicity.setType("double")
output_helicity.setRank(1)
output_helicity.setMaximumIndices(["Nvol"])

output_mu = Variable("mu")
output_mu.setDescription(r"helicity-multiplier, \f$\mu\f$, in each volume")
output_mu.setType("double")
output_mu.setRank(1)
output_mu.setMaximumIndices(["Nvol"])

output_tflux = Variable("tflux")
output_tflux.setDescription(r"toroidal flux, \f$\psi_t\f$, enclosed by each interface")
output_tflux.setUnit("Wb")
output_tflux.setType("double")
output_tflux.setRank(1)
output_tflux.setMaximumIndices(["Nvol"])

output_pflux = Variable("pflux")
output_pflux.setDescription(r"poloidal flux, \f$\psi_p\f$, enclosed by each interface")
output_pflux.setUnit("Wb")
output_pflux.setType("double")
output_pflux.setRank(1)
output_pflux.setMaximumIndices(["Nvol"])

output_volume = Variable("volume")
output_volume.setDescription(r"total volume = $\sum V_v$")
output_volume.setType("double")
output_volume.setRank(1)
output_volume.setMaximumIndices(["Nvol"])

output_Mrad = Variable("Mrad")
output_Mrad.setDescription(r"maximum radial (Chebyshev) resolution")
output_Mrad.setType("int")

output_TT = Variable("TT")
output_TT.setDescription(r"Chebyshev polynomials, $T_l$, and their derivatives, evaluated at $s=\pm 1$")
output_TT.setType("double")
output_TT.setRank(3)
output_TT.setStartingIndices(["0", "0", "0"])
output_TT.setMaximumIndices(["Mrad", "1", "1"])

output_Btemn = Variable("Btemn")
output_Btemn.setDescription( r"cosine harmonics of the covariant poloidal field, "
                            +r"i.e. $[[B_{\t,j}]]$ evaluated on the inner and outer interface in each volume")
output_Btemn.setUnit("T")
output_Btemn.setType("double")
output_Btemn.setRank(3)
output_Btemn.setStartingIndices(["1", "0", "1"])
output_Btemn.setMaximumIndices(["mn", "1", "Mvol"])

output_Bzemn = Variable("Bzemn")
output_Bzemn.setDescription( r"cosine harmonics of the covariant toroidal field, "
                            +r"i.e. $[[B_{\z,j}]]$ evaluated on the inner and outer interface in each volume")
output_Bzemn.setUnit("T")
output_Bzemn.setType("double")
output_Bzemn.setRank(3)
output_Bzemn.setStartingIndices(["1", "0", "1"])
output_Bzemn.setMaximumIndices(["mn", "1", "Mvol"])

output_Btomn = Variable("Btomn")
output_Btomn.setDescription( r"the sine harmonics of the covariant poloidal field, "
                            +r"i.e. $[[B_{\t,j}]]$ evaluated on the inner and outer interface in each volume")
output_Btomn.setUnit("T")
output_Btomn.setType("double")
output_Btomn.setRank(3)
output_Btomn.setStartingIndices(["1", "0", "1"])
output_Btomn.setMaximumIndices(["mn", "1", "Mvol"])

output_Bzomn = Variable("Bzomn")
output_Bzomn.setDescription( r"the sine harmonics of the covariant toroidal field, "
                            +r"i.e. $[[B_{\z,j}]]$ evaluated on the inner and outer interface in each volume")
output_Bzomn.setUnit("T")
output_Bzomn.setType("double")
output_Bzomn.setRank(3)
output_Bzomn.setStartingIndices(["1", "0", "1"])
output_Bzomn.setMaximumIndices(["mn", "1", "Mvol"])

output_lmns = Variable("lmns")
output_lmns.setDescription(r"resolution of the straight-fieldline transformation")
output_lmns.setType("int")

vars_output = [output_mn,
               output_im,
               output_in,
               output_Mvol,
               output_Rbc,
               output_Zbs,
               output_Rbs,
               output_Zbc,
               output_ForceErr,
               output_Ivolume,
               output_IPDt,
               output_adiabatic,
               output_helicity,
               output_mu,
               output_tflux,
               output_pflux,
               output_volume,
               output_Mrad,
               output_TT,
               output_Btemn,
               output_Bzemn,
               output_Btomn,
               output_Bzomn,
               output_lmns]


#-----------------------------------------------------------------------------#
# code generation starts here
#-----------------------------------------------------------------------------#
print("definition done, now starting code generation")

# functionality to be covered in the following languages:
# Fortran, Python, Matlab, Java, C, C++
# 1.1 input variables declaration
# 1.2  read input variables from namelist
# 1.3 write input variables   to namelist
# 1.4  read input variables from HDF5 file
# 1.5 write input variables   to HDF5 file
# 2.1 output variables declaration
# 2.2  read output variables from HDF5 file
# 2.3 write output variables   to HDF5 file

# to be used in SPEC itself (all in Fortran):
# 1.1: allglobal#inputlist
# 1.4: readin_hdf5()
# 1.5: mirror_input_to_outfile()
# 2.1: allglobal
# 2.3: sphdf5


###############################################################################
# 1.1: generate Fortran declarations of the input quantities of SPEC
###############################################################################

# actually generate Fortran module for reading SPEC output files
def genFortranDefInputlist():
    
    creation_tag = get_creation_tag()
    moduleName = "inputlist"
    
    # dry-run declaration to determine maximum declaration length for doc indentation
    maxLength = 0
    for nml in input_namelists:
        for var in nml.variables:
            declLen = len(Fortran.declareVariable(var, attachDescription=False))
            if declLen>maxLength: maxLength = declLen
    #print("maximum decl. length: "+str(maxLength))
    
    fortranFilename = os.path.join(".", moduleName+".f90")
    print("creating Fortran inputlist definition into '"+fortranFilename+"'")
    
    relative_path_to_this_file = relname(__file__, fortranFilename)
    with open(fortranFilename, "w") as f:
        f.write("! AUTO-GENERATED BY "+relative_path_to_this_file
                +"; DO NOT COMMIT CHANGES TO THIS FILE !"+"\n"
                "! "+creation_tag+"\n")

        f.write(r"!> @file "+moduleName+".f90"+r"""
!> \brief Input namelists
!> \addtogroup grp_global
!> @{
""")
    
        f.write(r"""module """+moduleName+r"""
!> \brief Input namelists
!> \addtogroup grp_global
!> @{
""")
        
        # parameters: maximum array dimensions
        for param in params_maxDims:
            f.write(Fortran.declareVariable(param, refDeclLength=maxLength)+"\n")

        brief = r"\brief "

        # input Variables, i.e. namelist contents
        for nml in input_namelists:
            
            nml_desc_indented = indented(len(brief), toDoc(nml.description), " ")
            commented_brief = Fortran.commentOut(brief+nml_desc_indented[len(brief):])
            
            f.write(r"!> \addtogroup grp_global_"+nml.name+" "+nml.name+"\n"
                    +commented_brief+"\n"
                    +"!> @{\n")
                                           
            for var in nml.variables:
                f.write(Fortran.declareVariable(var, refDeclLength=maxLength)+"\n")
            f.write("!> @}\n")
    
        # namelist declarations
        for nml in input_namelists:
            f.write(Fortran.declareNamelist(nml)+"\n\n")
        
        f.write("!> @}\n")
        f.write("end module inputlist\n")
# end of genFortranDefInputlist

###############################################################################
# 1.4: generate Fortran routine to read input data from HDF5 file
###############################################################################
def genFortranReadInputlistFromHDF5():
    
    creation_tag = get_creation_tag()
    moduleName = "readin_h5"
    
    fortranFilename = os.path.join(".", moduleName+".f90")
    print("creating Fortran inputlist reading module (HDF5) into '"+fortranFilename+"'")
    
    relative_path_to_this_file = relname(__file__, fortranFilename)
    with open(fortranFilename, "w") as f:
        f.write("! AUTO-GENERATED BY "+relative_path_to_this_file
                +"; DO NOT COMMIT CHANGES TO THIS FILE !"+"\n"
                "! "+creation_tag+"\n")
        
        f.write("subroutine "+moduleName+"(filename)\n")
        f.write("  use hdf5\n")
        f.write("  use inputlist\n")
        f.write("  use allglobal, only: version\n")
        f.write("  implicit none\n")
        f.write("  character, intent(in) :: filename*255\n")
        f.write("  int(hid_t) :: file_id\n")
        
        f.write(indented(2, Fortran.readHdf5Group("file_id", "input", input_namelists, 2, " "), " ")+"\n")
        
        f.write("end subroutine readin_h5\n")
# end genFortranReadInputlistFromHDF5












###############################################################################
# 2.2: generate Fortran routine to read output data from HDF5 file
###############################################################################
    
#     # begin code for root group (== enclosing class)
#     with open(fortranFilename, "w") as f:
    
        
        
#         # custom datatypes come first
#         for dtype in s.getDatatypes():
#             f.write(Fortran.genType(dtype.name, dtype.items)+'\n')
            
#         # we need to reverse the definition order so that types which are used inside other types
#         # are already defined when used
#         reverse_groupStack =  []
        
#         groupStack = []
#         groupStack.append(s.rootGroup)
#         while len(groupStack)>0:
#             currentGroup = groupStack[-1]
#             groupStack = groupStack[:-1]
            
#             if type(currentGroup)==Group:
#                 reverse_groupStack.append(currentGroup)
        
#             for item in currentGroup.items:
#                 if type(item)==Group:
#                     groupStack.append(item)
        
#         # iterate in reverse order over the discovered variables to generate type definitions in correct order
#         for currentGroup in reverse_groupStack[::-1]:
#             f.write(Fortran.genType(currentGroup.name, currentGroup.items)+'\n')
        
#         f.write("contains\n")
        
#         # initial code of loading routine
#         Fortran.startLoader(f)
        
#         # loop over all variables again and put the loader code for each of them one after another
#         for currentGroup in reverse_groupStack[::-1]:
#             for item in currentGroup.items:
#                 if type(item)==Dataset:
#                     Fortran.loadItem(f, item)
            
#         # finalizing code of loading routine
#         Fortran.endLoader(f)
        
#         # write the freeSpec subroutine to free the memory it occupied
#         Fortran.startFree(f)
        
#         for currentGroup in reverse_groupStack[::-1]:
#             for item in currentGroup.items:
#                 if type(item)==Dataset:
#                     Fortran.freeItem(f, item)
        
#         # finalizing code of freeing routine
#         Fortran.endFree(f)
        
#         f.write("end module "+moduleName+"\n")
    
#         # write demo code
#         #fortran_demoLoader(f)









# # from adf import indented, toDoc
# # from genFortran import declareVariable, declareNamelist






if __name__=="__main__":
    genFortranDefInputlist()
    genFortranReadInputlistFromHDF5()