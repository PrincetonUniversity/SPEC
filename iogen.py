#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is the definition script which controls the generation of the input file
and output file reading and writing routines for the SPEC MRxMHD
equilibrium code.

In the first part of this script, the input file format and the output file
format are declared in an abstract way.

Based on these declarations, reading and writing routines in various
programming languages are generated automagically.
The code generation functionality is encapsulated in the
Interface Definition Framework, which can be found at:
https://github.com/jonathanschilling/idf

If you want to have an input or output idf.Variable added to SPEC, this script
is the place to put it. Then run it, commit your changed source code to Git
and compile SPEC again to use your changes.

In order to get the correct order of the comments, you should use Python >=3.7.
From Python 3.7 on, it is guaranteed that the insertion order of dict()
items is kept and this is relied on in this code when iterating over the keys
of the documentation parts.
More on this: https://docs.python.org/3/whatsnew/3.7.html

@author: Jonathan Schilling (jonathan.schilling@ipp.mpg.de)
"""

# framework for source code generation
idfPath = "/data/jonathan/work/code/idf"
import sys
if not idfPath in sys.path:
    sys.path.insert(0, idfPath)
import idf

# define the input quantities for SPEC

MNvol = idf.Variable("MNvol")
MNvol.setDescription(r"maximum value of \c Nvol")
MNvol.setType("int")
MNvol.setDefaultValue(256)
MNvol.setIsParameter(True)


MMpol = idf.Variable("MMpol")
MMpol.setDescription(r"maximum value of \c Mpol")
MMpol.setType("int")
MMpol.setDefaultValue(64)
MMpol.setIsParameter(True)


MNtor = idf.Variable("MNtor")
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
input_physics_Igeometry = idf.Variable("Igeometry")
input_physics_Igeometry.setDescription({r"selects Cartesian, cylindrical or toroidal geometry":
                                        [r"\c Igeometry=1 : Cartesian; geometry determined by \f$R\f$",
                                         r"\c Igeometry=2 : cylindrical; geometry determined by \f$R\f$",
                                         r"\c Igeometry=3 : toroidal; geometry determined by \f$R\f$ *and* \f$Z\f$"]
                                        })
input_physics_Igeometry.setType("int")
input_physics_Igeometry.setDefaultValue(3)


input_physics_Istellsym = idf.Variable("Istellsym")
input_physics_Istellsym.setDescription(r"stellarator symmetry is enforced if \c Istellsym=1")
input_physics_Istellsym.setType("int")
input_physics_Istellsym.setDefaultValue(1)


input_physics_Lfreebound = idf.Variable("Lfreebound")
input_physics_Lfreebound.setDescription(r"compute vacuum field surrounding plasma")
input_physics_Lfreebound.setType("int")
input_physics_Lfreebound.setDefaultValue(0)


input_physics_phiedge = idf.Variable("phiedge")
input_physics_phiedge.setDescription(r"total enclosed toroidal magnetic flux")
input_physics_phiedge.setType("double")
input_physics_phiedge.setDefaultValue(1.0)
input_physics_phiedge.setUnit("Vs")


input_physics_curtor = idf.Variable("curtor")
input_physics_curtor.setDescription(r"total enclosed (toroidal) plasma current")
input_physics_curtor.setType("double")
input_physics_curtor.setDefaultValue(0.0)


input_physics_curpol = idf.Variable("curpol")
input_physics_curpol.setDescription(r"total enclosed (poloidal) linking current")
input_physics_curpol.setType("double")
input_physics_curpol.setDefaultValue(0.0)


input_physics_gamma = idf.Variable("gamma")
input_physics_gamma.setDescription(r"adiabatic index; cannot set \f$|\gamma| = 1\f$")
input_physics_gamma.setType("double")
input_physics_gamma.setDefaultValue(0.0)


input_physics_Nfp = idf.Variable("Nfp")
input_physics_Nfp.setDescription({r"field periodicity":
                                  [r"all Fourier representations are of the form \f$\cos(m\theta-n N \zeta)\f$, \f$\sin(m\theta-n N \zeta)\f$, where \f$N\equiv\f$\c Nfp",
                                   r"constraint: \c Nfp >= 1"]
                                  })
input_physics_Nfp.setType("int")
input_physics_Nfp.setDefaultValue(1)


input_physics_Nvol = idf.Variable("Nvol")
input_physics_Nvol.setDescription({r"number of volumes":
                                   [r"each volume \f${\cal V}_l\f$ is bounded by the \f${\cal I}_{l-1}\f$ and \f${\cal I}_{l}\f$ interfaces",
                                    r"note that in cylindrical or toroidal geometry, \f${\cal I}_{0}\f$ is the degenerate coordinate axis",
                                    r"constraint: \c Nvol<=MNvol"]
                                   })
input_physics_Nvol.setType("int")
input_physics_Nvol.setDefaultValue(1)


input_physics_Mpol = idf.Variable("Mpol")
input_physics_Mpol.setDescription({r"number of poloidal Fourier harmonics":
                                   [ r"all Fourier representations of doubly-periodic functions are of the form"+"\n"
                                    +r"\f{eqnarray}{ f(\theta,\zeta) & = & \sum_{n=0}^{\texttt{Ntor}} f_{0,n}\cos(-n \, \texttt{Nfp} \, \zeta)"+"\n"
                                    +r"\sum_{m=1}^{\texttt{Mpol}}\sum_{n=\texttt{-Ntor}}^{\texttt{Ntor}} f_{m,n}\cos(m\theta-n \, \texttt{Nfp} \, \zeta),"+"\n"
                                    +r"\f}"+"\n"
                                    +r"Internally these \"double\" summations are written as a \"single\" summation,"+"\n"
                                    +r"e.g. \f$f(\theta,\zeta) = \sum_j f_j \cos(m_j\theta-n_j\zeta)\f$."]
                                   })
input_physics_Mpol.setType("int")
input_physics_Mpol.setDefaultValue(0)


input_physics_Ntor = idf.Variable("Ntor")
input_physics_Ntor.setDescription({r"number of toroidal Fourier harmonics":
                                   [ r"all Fourier representations of doubly-periodic functions are of the form"+"\n"
                                    +r"\f{eqnarray}{ f(\theta,\zeta) & = & \sum_{n=0}^{\texttt{Ntor}} f_{0,n}\cos(-n \, \texttt{Nfp} \, \zeta)"+"\n"
                                    +r"\sum_{m=1}^{\texttt{Mpol}}\sum_{n=\texttt{-Ntor}}^{\texttt{Ntor}} f_{m,n}\cos(m\theta-n \, \texttt{Nfp} \, \zeta),"+"\n"
                                    +r"\f}"+"\n"
                                    +r"Internally these \"double\" summations are written as a \"single\" summation,"+"\n"
                                    +r"e.g. \f$f(\theta,\zeta) = \sum_j f_j \cos(m_j\theta-n_j\zeta)\f$."]
                                   })
input_physics_Ntor.setType("int")
input_physics_Ntor.setDefaultValue(0)


input_physics_Lrad = idf.Variable("Lrad")
input_physics_Lrad.setDescription({r"Chebyshev resolution in each volume":
                                   [r"constraint : \c Lrad(1:Mvol) >= 2"]
                                   })
input_physics_Lrad.setType("int")
input_physics_Lrad.setRank(1)
input_physics_Lrad.setDefaultValue(4)
input_physics_Lrad.setMaximumIndices([r"MNvol+1"])


input_physics_Lconstraint = idf.Variable("Lconstraint")
input_physics_Lconstraint.setDescription({r"selects constraints; primarily used in ma02aa() and mp00ac()":
                                          [ r"if \c Lconstraint=-1, then in the plasma regions \f$\Delta\psi_t\f$, \f$\mu\f$ and \f$\Delta \psi_p\f$ are *not* varied"+"\n"
                                           +r"and in the vacuum region (only for free-boundary) \f$\Delta\psi_t\f$ and \f$\Delta \psi_p\f$ are *not* varied, and \f$\mu = 0\f$",
                                            r"if \c Lconstraint=0, then in the plasma regions \f$\Delta\psi_t\f$, \f$\mu\f$ and \f$\Delta \psi_p\f$ are *not* varied"+"\n"
                                           +r"and in the vacuum region (only for free-boundary) \f$\Delta\psi_t\f$ and \f$\Delta \psi_p\f$ are varied to match the"+"\n"
                                           +r"prescribed plasma current, \c curtor, and the \"linking\" current, \c curpol, and \f$\mu = 0\f$",
                                            r"if \c Lconstraint=1, then in the plasma regions \f$\mu\f$ and \f$\Delta\psi_p\f$ are adjusted"+"\n"
                                           +r"in order to satisfy the inner and outer interface transform constraints"+"\n"
                                           +r"(except in the simple torus, where the enclosed poloidal flux is irrelevant,"+"\n"
                                           +r"and only \f$\mu\f$ is varied to satisfy the outer interface transform constraint);"+"\n"
                                           +r"and in the vacuum region \f$\Delta\psi_t\f$ and \f$\Delta \psi_p\f$ are varied to match the transform constraint on the boundary"+"\n"
                                           +r"and to obtain the prescribed linking current, \c curpol, and \f$\mu = 0\f$",
                                            r"\todo if \c Lconstraint=2, under reconstruction"]
                                          })
input_physics_Lconstraint.setType("int")
input_physics_Lconstraint.setDefaultValue(-1)


input_physics_tflux = idf.Variable("tflux")
input_physics_tflux.setDescription({r"toroidal flux, \f$\psi_t\f$, enclosed by each interface":
                                    [ r"For each of the plasma volumes, this is a constraint: \c tflux is *not* varied",
                                      r"For the vacuum region (only if \c Lfreebound==1), \c tflux  may be allowed to vary to match constraints",
                                      r"Note that \c tflux  will be normalized so that \c tflux(Nvol) = 1.0,"+"\n"
                                     +r"so that \c tflux  is arbitrary up to a scale factor",
                                      r"\sa phiedge"]
                                    })
input_physics_tflux.setType("double")
input_physics_tflux.setRank(1)
input_physics_tflux.setDefaultValue(0.0)
input_physics_tflux.setMaximumIndices([r"MNvol+1"])


input_physics_pflux = idf.Variable("pflux")
input_physics_pflux.setDescription(r"poloidal flux, \f$\psi_p\f$, enclosed by each interface")
input_physics_pflux.setType("double")
input_physics_pflux.setRank(1)
input_physics_pflux.setDefaultValue(0.0)
input_physics_pflux.setMaximumIndices([r"MNvol+1"])


input_physics_helicity = idf.Variable("helicity")
input_physics_helicity.setDescription({r"helicity, \f${\cal K}\f$, in each volume, \f${\cal V}_i\f$":
                                       [r"on exit, \c helicity  is set to the computed values of \f${\cal K} \equiv \int {\bf A}\cdot{\bf B}\;dv\f$"]
                                       })
input_physics_helicity.setType("double")
input_physics_helicity.setRank(1)
input_physics_helicity.setDefaultValue(0.0)
input_physics_helicity.setMaximumIndices([r"MNvol+1"])


input_physics_pscale = idf.Variable("pscale")
input_physics_pscale.setDescription({r"pressure scale factor":
                                     [r"the initial pressure profile is given by \c pscale  \f$*\f$ \c pressure"]
                                    })
input_physics_pscale.setType("double")
input_physics_pscale.setDefaultValue(0.0) #TODO maybe this should be 1.0?


input_physics_pressure = idf.Variable("pressure")
input_physics_pressure.setDescription({r"pressure in each volume":
                                       [ r"The pressure is *not* held constant, but \f$p_l V_l^\gamma = P_l\f$ *is* held constant,"+"\n"
                                        +r"where \f$P_l\f$ is determined by the initial pressures and the initial volumes, \f$V_l\f$.",
                                         r"Note that if \c gamma==0.0, then \f$p_l \equiv P_l\f$.",
                                         r"On output, the pressure is given by \f$p_l = P_l/V_l^\gamma\f$, where \f$V_l\f$ is the final volume.",
                                         r"\c pressure is only used in calculation of interface force-balance."]
                                       })
input_physics_pressure.setType("double")
input_physics_pressure.setRank(1)
input_physics_pressure.setDefaultValue(0.0)
input_physics_pressure.setMaximumIndices([r"MNvol+1"])


input_physics_Ladiabatic = idf.Variable("Ladiabatic")
input_physics_Ladiabatic.setDescription({r"logical flag":
                                         [r"If \c Ladiabatic==0, the adiabatic constants are determined by the initial pressure and volume.",
                                          r"If \c Ladiabatic==1, the adiabatic constants are determined by the given input \c adiabatic."]
                                         })
input_physics_Ladiabatic.setType("int")
input_physics_Ladiabatic.setDefaultValue(0)


input_physics_adiabatic = idf.Variable("adiabatic")
input_physics_adiabatic.setDescription({r"adiabatic constants in each volume":
                                        [r"The pressure is *not* held constant, but \f$p_l V_l^\gamma = P_l \equiv\f$\c adiabatic is constant.",
                                         r"Note that if \c gamma==0.0, then \c pressure==adiabatic.",
                                         r"\c pressure is only used in calculation of interface force-balance."]
                                        })
input_physics_adiabatic.setType("double")
input_physics_adiabatic.setRank(1)
input_physics_adiabatic.setDefaultValue(0.0)
input_physics_adiabatic.setMaximumIndices([r"MNvol+1"])


input_physics_mu = idf.Variable("mu")
input_physics_mu.setDescription(r"helicity-multiplier, \f$\mu\f$, in each volume")
input_physics_mu.setType("double")
input_physics_mu.setRank(1)
input_physics_mu.setDefaultValue(0.0)
input_physics_mu.setMaximumIndices([r"MNvol+1"])

# TODO: Ivolume
# TODO: Isurf



input_physics_pl = idf.Variable("pl")
input_physics_pl.setDescription( r"\"inside\" interface rotational-transform is \f$\mbox{$\,\iota\!\!$-} = (p_l+\gamma p_r)/(q_l+\gamma q_r)\f$,"+"\n"
                                +r"where \f$\gamma\f$ is the golden mean, \f$\gamma = (1 + \sqrt 5 ) / 2 \f$."+"\n"
                                +r"If both \f$q_l = 0\f$ *and* \f$q_r = 0\f$, then the (inside) interface rotational-transform is defined by \c iota .")
input_physics_pl.setType("int")
input_physics_pl.setRank(1)
input_physics_pl.setDefaultValue(0)
input_physics_pl.setStartingIndices([r"0"])
input_physics_pl.setMaximumIndices([r"MNvol"])


input_physics_ql = idf.Variable("ql")
input_physics_ql.setDescription( r"\"inside\" interface rotational-transform is \f$\mbox{$\,\iota\!\!$-} = (p_l+\gamma p_r)/(q_l+\gamma q_r)\f$,"+"\n"
                                +r"where \f$\gamma\f$ is the golden mean, \f$\gamma = (1 + \sqrt 5 ) / 2 \f$."+"\n"
                                +r"If both \f$q_l = 0\f$ *and* \f$q_r = 0\f$, then the (inside) interface rotational-transform is defined by \c iota .")
input_physics_ql.setType("int")
input_physics_ql.setRank(1)
input_physics_ql.setDefaultValue(0)
input_physics_ql.setStartingIndices([r"0"])
input_physics_ql.setMaximumIndices([r"MNvol"])


input_physics_pr = idf.Variable("pr")
input_physics_pr.setDescription( r"\"inside\" interface rotational-transform is \f$\mbox{$\,\iota\!\!$-} = (p_l+\gamma p_r)/(q_l+\gamma q_r)\f$,"+"\n"
                                +r"where \f$\gamma\f$ is the golden mean, \f$\gamma = (1 + \sqrt 5 ) / 2 \f$."+"\n"
                                +r"If both \f$q_l = 0\f$ *and* \f$q_r = 0\f$, then the (inside) interface rotational-transform is defined by \c iota .")
input_physics_pr.setType("int")
input_physics_pr.setRank(1)
input_physics_pr.setDefaultValue(0)
input_physics_pr.setStartingIndices([r"0"])
input_physics_pr.setMaximumIndices([r"MNvol"])


input_physics_qr = idf.Variable("qr")
input_physics_qr.setDescription( r"\"inside\" interface rotational-transform is \f$\mbox{$\,\iota\!\!$-} = (p_l+\gamma p_r)/(q_l+\gamma q_r)\f$,"+"\n"
                                +r"where \f$\gamma\f$ is the golden mean, \f$\gamma = (1 + \sqrt 5 ) / 2 \f$."+"\n"
                                +r"If both \f$q_l = 0\f$ *and* \f$q_r = 0\f$, then the (inside) interface rotational-transform is defined by \c iota .")
input_physics_qr.setType("int")
input_physics_qr.setRank(1)
input_physics_qr.setDefaultValue(0)
input_physics_qr.setStartingIndices([r"0"])
input_physics_qr.setMaximumIndices([r"MNvol"])

input_physics_iota = idf.Variable("iota")
input_physics_iota.setDescription({r"rotational-transform, \f$\mbox{$\,\iota\!\!$-}\f$, on inner side of each interface":
                                   [r"only relevant if illogical input for \c ql and \c qr are provided"]
                                   })
input_physics_iota.setType("double")
input_physics_iota.setRank(1)
input_physics_iota.setDefaultValue(0.0)
input_physics_iota.setStartingIndices([r"0"])
input_physics_iota.setMaximumIndices([r"MNvol"])


input_physics_lp = idf.Variable("lp")
input_physics_lp.setDescription( r"\"outer\" interface rotational-transform is \f$\mbox{$\,\iota\!\!$-} = (p_l+\gamma p_r)/(q_l+\gamma q_r)\f$,"+"\n"
                                +r"where \f$\gamma\f$ is the golden mean, \f$\gamma = (1 + \sqrt 5 ) / 2\f$."+"\n"
                                +r"If both \f$q_l = 0\f$ *and* \f$q_r = 0\f$, then the (outer) interface rotational-transform is defined by \c oita .")
input_physics_lp.setType("int")
input_physics_lp.setRank(1)
input_physics_lp.setDefaultValue(0)
input_physics_lp.setStartingIndices([r"0"])
input_physics_lp.setMaximumIndices([r"MNvol"])


input_physics_lq = idf.Variable("lq")
input_physics_lq.setDescription( r"\"outer\" interface rotational-transform is \f$\mbox{$\,\iota\!\!$-} = (p_l+\gamma p_r)/(q_l+\gamma q_r)\f$,"+"\n"
                                +r"where \f$\gamma\f$ is the golden mean, \f$\gamma = (1 + \sqrt 5 ) / 2\f$."+"\n"
                                +r"If both \f$q_l = 0\f$ *and* \f$q_r = 0\f$, then the (outer) interface rotational-transform is defined by \c oita .")
input_physics_lq.setType("int")
input_physics_lq.setRank(1)
input_physics_lq.setDefaultValue(0)
input_physics_lq.setStartingIndices([r"0"])
input_physics_lq.setMaximumIndices([r"MNvol"])


input_physics_rp = idf.Variable("rp")
input_physics_rp.setDescription( r"\"outer\" interface rotational-transform is \f$\mbox{$\,\iota\!\!$-} = (p_l+\gamma p_r)/(q_l+\gamma q_r)\f$,"+"\n"
                                +r"where \f$\gamma\f$ is the golden mean, \f$\gamma = (1 + \sqrt 5 ) / 2\f$."+"\n"
                                +r"If both \f$q_l = 0\f$ *and* \f$q_r = 0\f$, then the (outer) interface rotational-transform is defined by \c oita .")
input_physics_rp.setType("int")
input_physics_rp.setRank(1)
input_physics_rp.setDefaultValue(0)
input_physics_rp.setStartingIndices([r"0"])
input_physics_rp.setMaximumIndices([r"MNvol"])


input_physics_rq = idf.Variable("rq")
input_physics_rq.setDescription( r"\"outer\" interface rotational-transform is \f$\mbox{$\,\iota\!\!$-} = (p_l+\gamma p_r)/(q_l+\gamma q_r)\f$,"+"\n"
                                +r"where \f$\gamma\f$ is the golden mean, \f$\gamma = (1 + \sqrt 5 ) / 2\f$."+"\n"
                                +r"If both \f$q_l = 0\f$ *and* \f$q_r = 0\f$, then the (outer) interface rotational-transform is defined by \c oita .")
input_physics_rq.setType("int")
input_physics_rq.setRank(1)
input_physics_rq.setDefaultValue(0)
input_physics_rq.setStartingIndices([r"0"])
input_physics_rq.setMaximumIndices([r"MNvol"])


input_physics_oita = idf.Variable("oita")
input_physics_oita.setDescription({r"rotational-transform, \f$\mbox{$\,\iota\!\!$-}\f$, on outer side of each interface":
                                   [r"only relevant if illogical input for \c ql and \c qr are provided"]
                                   })
input_physics_oita.setType("double")
input_physics_oita.setRank(1)
input_physics_oita.setDefaultValue(0.0)
input_physics_oita.setStartingIndices([r"0"])
input_physics_oita.setMaximumIndices([r"MNvol"])


input_physics_mupftol = idf.Variable("mupftol")
input_physics_mupftol.setDescription({r"accuracy to which \f$\mu\f$ and \f$\Delta\psi_p\f$ are required":
                                      [r"only relevant if constraints on transform, enclosed currents etc. are to be satisfied iteratively, see \c Lconstraint"]
                                      })
input_physics_mupftol.setType("double")
input_physics_mupftol.setDefaultValue(1.0e-16)


input_physics_mupfits = idf.Variable("mupfits")
input_physics_mupfits.setDescription({r"an upper limit on the transform/helicity constraint iterations":
                                      [r"only relevant if constraints on transform, enclosed currents etc. are to be satisfied iteratively, see \c Lconstraint",
                                       r"constraint: \c mupfits > 0"]
                                      })
input_physics_mupfits.setType("int")
input_physics_mupfits.setDefaultValue(8)


input_physics_rpol = idf.Variable("rpol")
input_physics_rpol.setDescription({r"poloidal extent of slab (effective radius)":
                                   [r"only relevant if \c Igeometry==1",
                                    r"poloidal size is \f$L = 2\pi*\f$\c rpol"]
                                   })
input_physics_rpol.setType("double")
input_physics_rpol.setDefaultValue(1.0)
input_physics_rpol.setUnit("m")


input_physics_rtor = idf.Variable("rtor")
input_physics_rtor.setDescription({r"toroidal extent of slab (effective radius)":
                                   [r"only relevant if \c Igeometry==1",
                                    r"toroidal size is \f$L = 2\pi*\f$\c rtor"]
                                   })
input_physics_rtor.setType("double")
input_physics_rtor.setDefaultValue(1.0)
input_physics_rtor.setUnit("m")

# TODO: Lreflect

input_physics_Rac = idf.Variable("Rac")
input_physics_Rac.setDescription(r"    stellarator symmetric coordinate axis; R; cosine")
input_physics_Rac.setType("double")
input_physics_Rac.setRank(1)
input_physics_Rac.setDefaultValue(0.0)
input_physics_Rac.setUnit("m")
input_physics_Rac.setStartingIndices([r"0"])
input_physics_Rac.setMaximumIndices([r"MNtor"])

input_physics_Zas = idf.Variable("Zas")
input_physics_Zas.setDescription(r"    stellarator symmetric coordinate axis; Z;   sine")
input_physics_Zas.setType("double")
input_physics_Zas.setRank(1)
input_physics_Zas.setDefaultValue(0.0)
input_physics_Zas.setUnit("m")
input_physics_Zas.setStartingIndices([r"0"])
input_physics_Zas.setMaximumIndices([r"MNtor"])

input_physics_Ras = idf.Variable("Ras")
input_physics_Ras.setDescription(r"non-stellarator symmetric coordinate axis; R;   sine")
input_physics_Ras.setType("double")
input_physics_Ras.setRank(1)
input_physics_Ras.setDefaultValue(0.0)
input_physics_Ras.setUnit("m")
input_physics_Ras.setStartingIndices([r"0"])
input_physics_Ras.setMaximumIndices([r"MNtor"])

input_physics_Zac = idf.Variable("Zac")
input_physics_Zac.setDescription(r"non-stellarator symmetric coordinate axis; Z; cosine")
input_physics_Zac.setType("double")
input_physics_Zac.setRank(1)
input_physics_Zac.setDefaultValue(0.0)
input_physics_Zac.setUnit("m")
input_physics_Zac.setStartingIndices([r"0"])
input_physics_Zac.setMaximumIndices([r"MNtor"])



input_physics_Rbc = idf.Variable("Rbc")
input_physics_Rbc.setDescription(r"    stellarator symmetric boundary components; R; cosine")
input_physics_Rbc.setType("double")
input_physics_Rbc.setRank(2)
input_physics_Rbc.setDefaultValue(0.0)
input_physics_Rbc.setUnit("m")
input_physics_Rbc.setStartingIndices([r"-MNtor", r"-MMpol"])
input_physics_Rbc.setMaximumIndices([r"MNtor", r"MMpol"])

input_physics_Zbs = idf.Variable("Zbs")
input_physics_Zbs.setDescription(r"    stellarator symmetric boundary components; Z;   sine")
input_physics_Zbs.setType("double")
input_physics_Zbs.setRank(2)
input_physics_Zbs.setDefaultValue(0.0)
input_physics_Zbs.setUnit("m")
input_physics_Zbs.setStartingIndices([r"-MNtor", r"-MMpol"])
input_physics_Zbs.setMaximumIndices([r"MNtor", r"MMpol"])

input_physics_Rbs = idf.Variable("Rbs")
input_physics_Rbs.setDescription(r"non-stellarator symmetric boundary components; R;   sine")
input_physics_Rbs.setType("double")
input_physics_Rbs.setRank(2)
input_physics_Rbs.setDefaultValue(0.0)
input_physics_Rbs.setUnit("m")
input_physics_Rbs.setStartingIndices([r"-MNtor", r"-MMpol"])
input_physics_Rbs.setMaximumIndices([r"MNtor", r"MMpol"])

input_physics_Zbc = idf.Variable("Zbc")
input_physics_Zbc.setDescription(r"non-stellarator symmetric boundary components; Z; cosine")
input_physics_Zbc.setType("double")
input_physics_Zbc.setRank(2)
input_physics_Zbc.setDefaultValue(0.0)
input_physics_Zbc.setUnit("m")
input_physics_Zbc.setStartingIndices([r"-MNtor", r"-MMpol"])
input_physics_Zbc.setMaximumIndices([r"MNtor", r"MMpol"])



input_physics_Rwc = idf.Variable("Rwc")
input_physics_Rwc.setDescription(r"    stellarator symmetric boundary components of wall; R; cosine")
input_physics_Rwc.setType("double")
input_physics_Rwc.setRank(2)
input_physics_Rwc.setDefaultValue(0.0)
input_physics_Rwc.setUnit("m")
input_physics_Rwc.setStartingIndices([r"-MNtor", r"-MMpol"])
input_physics_Rwc.setMaximumIndices([r"MNtor", r"MMpol"])

input_physics_Zws = idf.Variable("Zws")
input_physics_Zws.setDescription(r"    stellarator symmetric boundary components of wall; Z;   sine")
input_physics_Zws.setType("double")
input_physics_Zws.setRank(2)
input_physics_Zws.setDefaultValue(0.0)
input_physics_Zws.setUnit("m")
input_physics_Zws.setStartingIndices([r"-MNtor", r"-MMpol"])
input_physics_Zws.setMaximumIndices([r"MNtor", r"MMpol"])

input_physics_Rws = idf.Variable("Rws")
input_physics_Rws.setDescription(r"non-stellarator symmetric boundary components of wall; R;   sine")
input_physics_Rws.setType("double")
input_physics_Rws.setRank(2)
input_physics_Rws.setDefaultValue(0.0)
input_physics_Rws.setUnit("m")
input_physics_Rws.setStartingIndices([r"-MNtor", r"-MMpol"])
input_physics_Rws.setMaximumIndices([r"MNtor", r"MMpol"])

input_physics_Zwc = idf.Variable("Zwc")
input_physics_Zwc.setDescription(r"non-stellarator symmetric boundary components of wall; Z; cosine")
input_physics_Zwc.setType("double")
input_physics_Zwc.setRank(2)
input_physics_Zwc.setDefaultValue(0.0)
input_physics_Zwc.setUnit("m")
input_physics_Zwc.setStartingIndices([r"-MNtor", r"-MMpol"])
input_physics_Zwc.setMaximumIndices([r"MNtor", r"MMpol"])


input_physics_Vns = idf.Variable("Vns")
input_physics_Vns.setDescription(r"    stellarator symmetric normal field at boundary; vacuum component;   sine")
input_physics_Vns.setType("double")
input_physics_Vns.setRank(2)
input_physics_Vns.setDefaultValue(0.0)
input_physics_Vns.setUnit("T")
input_physics_Vns.setStartingIndices([r"-MNtor", r"-MMpol"])
input_physics_Vns.setMaximumIndices([r"MNtor", r"MMpol"])

input_physics_Bns = idf.Variable("Bns")
input_physics_Bns.setDescription(r"    stellarator symmetric normal field at boundary; plasma component;   sine")
input_physics_Bns.setType("double")
input_physics_Bns.setRank(2)
input_physics_Bns.setDefaultValue(0.0)
input_physics_Bns.setUnit("T")
input_physics_Bns.setStartingIndices([r"-MNtor", r"-MMpol"])
input_physics_Bns.setMaximumIndices([r"MNtor", r"MMpol"])

input_physics_Vnc = idf.Variable("Vnc")
input_physics_Vnc.setDescription(r"non-stellarator symmetric normal field at boundary; vacuum component; cosine")
input_physics_Vnc.setType("double")
input_physics_Vnc.setRank(2)
input_physics_Vnc.setDefaultValue(0.0)
input_physics_Vnc.setUnit("T")
input_physics_Vnc.setStartingIndices([r"-MNtor", r"-MMpol"])
input_physics_Vnc.setMaximumIndices([r"MNtor", r"MMpol"])

input_physics_Bnc = idf.Variable("Bnc")
input_physics_Bnc.setDescription(r"non-stellarator symmetric normal field at boundary; plasma component; cosine")
input_physics_Bnc.setType("double")
input_physics_Bnc.setRank(2)
input_physics_Bnc.setDefaultValue(0.0)
input_physics_Bnc.setUnit("T")
input_physics_Bnc.setStartingIndices([r"-MNtor", r"-MMpol"])
input_physics_Bnc.setMaximumIndices([r"MNtor", r"MMpol"])


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

physicslist = idf.Namelist("physicslist")
physicslist.setDescription(r"The namelist \c physicslist controls the geometry, profiles, and numerical resolution.")
physicslist.addidf.Variables(vars_physicslist)

###############################################################################
# numericlist 
###############################################################################

input_numeric_Linitialize = idf.Variable("Linitialize")
input_numeric_Linitialize.setDescription({r"Used to initialize geometry using a regularization / extrapolation method":
                                           [ r"if \c Linitialize = \f$-I\f$ , where \f$I\f$ is a positive integer,"+"\n"
                                            +r"the geometry of the \f$i=1,N_V-I\f$ surfaces constructed by an extrapolation",
                                             r"if \c Linitialize=0, the geometry of the interior surfaces is provided after the namelists in the input file",
                                             r"if \c Linitialize=1, the interior surfaces will be intialized as \f$R_{l,m,n} = R_{N,m,n} \psi_{t,l}^{m/2}\f$,"+"\n"
                                            +r"where \f$R_{N,m,n}\f$ is the plasma boundary and \f$\psi_{t,l}\f$ is the given toroidal flux enclosed by the"+"\n"
                                            +r"\f$l\f$-th interface, normalized to the total enclosed toroidal flux;"+"\n"
                                            +r"a similar extrapolation is used for \f$Z_{l,m,n}\f$",
                                             r"Note that the Fourier harmonics of the boundary is *always* given by the \c Rbc and \c Zbs"+"\n"
                                            +r"given in \c physicslist.",
                                             r"if \c Linitialize=2, the interior surfaces *and the plasma boundary* will be intialized"+"\n"
                                            +r"as \f$R_{l,m,n} = R_{W,m,n} \psi_{t,l}^{m/2}\f$, where \f$R_{W,m,n}\f$ is the computational boundary"+"\n"
                                            +r"and \f$\psi_{t,l}\f$ is the given toroidal flux enclosed by the \f$l\f$-th interface, normalized to the total enclosed toroidal flux;"+"\n"
                                            +r"a similar extrapolation is used for \f$Z_{l,m,n}\f$",
                                             r"Note that, for free-boundary calculations, the Fourier harmonics of the computational boundary"+"\n"
                                            +r"are *always* given by the \c Rwc and \c Zws given in \c physicslist.",
                                             r"if \c Linitialize=1,2 , it is not required to provide the geometry of the interfaces after the namelists"]
                                          })
input_numeric_Linitialize.setType("int")
input_numeric_Linitialize.setDefaultValue(0)


input_numeric_LautoinitBn = idf.Variable("LautoinitBn")
input_numeric_LautoinitBn.setDescription({r"Used to initialize \f$B_{ns}\f$ using an initial fixed-boundary calculation":
                                          [r"only relevant if \c Lfreebound=1",
                                           r"user-supplied \c Bns will only be considered if \c LautoinitBn=0"]
                                          })
input_numeric_LautoinitBn.setType("int")
input_numeric_LautoinitBn.setDefaultValue(1)


input_numeric_Lzerovac = idf.Variable("Lzerovac")
input_numeric_Lzerovac.setDescription({r"Used to adjust vacuum field to cancel plasma field on computational boundary":
                                       [r"only relevant if \c Lfreebound=1"]
                                       })
input_numeric_Lzerovac.setType("int")
input_numeric_Lzerovac.setDefaultValue(0)


input_numeric_Ndiscrete = idf.Variable("Ndiscrete")
input_numeric_Ndiscrete.setDescription({r"resolution of the real space grid on which fast Fourier transforms are performed is given by \c Ndiscrete*Mpol*4":
                                        [r"constraint \c Ndiscrete>0"]
                                        })
input_numeric_Ndiscrete.setType("int")
input_numeric_Ndiscrete.setDefaultValue(2)


input_numeric_Nquad = idf.Variable("Nquad")
input_numeric_Nquad.setDescription({r"Resolution of the Gaussian quadrature":
                                    [ r"The resolution of the Gaussian quadrature, \f$\displaystyle \int \!\! f(s) ds = \sum_k \omega_k f(s_k)\f$,"+"\n"
                                     +r"in each volume is given by \c Iquad\f$_v\f$",
                                      r"\c Iquad\f$_v\f$ is set in preset()"]
                                    })
input_numeric_Nquad.setType("int")
input_numeric_Nquad.setDefaultValue(-1)


input_numeric_iMpol = idf.Variable("iMpol")
input_numeric_iMpol.setDescription({r"Fourier resolution of straight-fieldline angle on interfaces":
                                    [ r"the rotational-transform on the interfaces is determined by a transformation to the straight-fieldline angle,"+"\n"
                                     +r"with poloidal resolution given by \c iMpol",
                                      r"if \c iMpol<=0, then \c iMpol = Mpol - iMpol"]
                                    })
input_numeric_iMpol.setType("int")
input_numeric_iMpol.setDefaultValue(-4)


input_numeric_iNtor = idf.Variable("iNtor")
input_numeric_iNtor.setDescription({r"Fourier resolution of straight-fieldline angle on interfaces":
                                    [ r"the rotational-transform on the interfaces is determined by a transformation to the straight-fieldline angle,"+"\n"
                                     +r"with toroidal resolution given by \c iNtor",
                                      r"if \c iNtor<=0 then \c iNtor = Ntor - iNtor",
                                      r"if \c Ntor=0, then the toroidal resolution of the angle transformation is set \c lNtor=0"]
                                    })
input_numeric_iNtor.setType("int")
input_numeric_iNtor.setDefaultValue(-4)


input_numeric_Lsparse = idf.Variable("Lsparse")
input_numeric_Lsparse.setDescription({r"controls method used to solve for rotational-transform on interfaces":
                                      [ r"if \c Lsparse=0, the transformation to the straight-fieldline angle is computed in Fourier space"+"\n"
                                       +r"using a dense matrix solver, \c F04AAF",
                                        r"if \c Lsparse=1, the transformation to the straight-fieldline angle is computed in real space"+"\n"
                                       +r"using a dense matrix solver, \c F04ATF",
                                        r"if \c Lsparse=2, the transformation to the straight-fieldline angle is computed in real space"+"\n"
                                       +r"using a sparse matrix solver, \c F11DEF",
                                        r"if \c Lsparse=3, the different methods for constructing the straight-fieldline angle are compared"]
                                      })
input_numeric_Lsparse.setType("int")
input_numeric_Lsparse.setDefaultValue(0)


input_numeric_Lsvdiota = idf.Variable("Lsvdiota")
input_numeric_Lsvdiota.setDescription({r"controls method used to solve for rotational-transform on interfaces":
                                       [r"if \c Lsvdiota=0, use standard linear solver to construct straight fieldline angle transformation",
                                        r"if \c Lsvdiota=1, use SVD method to compute rotational-transform"],
                                       r"only relevant if \c Lsparse=0": None
                                       })
input_numeric_Lsvdiota.setType("int")
input_numeric_Lsvdiota.setDefaultValue(0)


input_numeric_imethod = idf.Variable("imethod")
input_numeric_imethod.setDescription({ r"controls iterative solution to sparse matrix"+"\n"
                                      +r"arising in real-space transformation to the straight-fieldline angle":
                                       [r"if \c imethod=1, the method is \c RGMRES",
                                        r"if \c imethod=2, the method is \c CGS",
                                        r"if \c imethod=3, the method is \c BICGSTAB"],
                                       r"only relevant if \c Lsparse=2; \see tr00ab() for details": None
                                      })
input_numeric_imethod.setType("int")
input_numeric_imethod.setDefaultValue(3)


input_numeric_iorder = idf.Variable("iorder")
input_numeric_iorder.setDescription({r"determines order of finite-difference approximation to the derivatives":
                                     [r"if \c iorder=2, second-order",
                                      r"if \c iorder=4, fourth-order",
                                      r"if \c iorder=6, sixth-order"],
                                     r"controls real-space grid resolution for constructing the straight-fieldline angle": None,
                                     r"only relevant if \c Lsparse>0": None
                                    })
input_numeric_iorder.setType("int")
input_numeric_iorder.setDefaultValue(2)


input_numeric_iprecon = idf.Variable("iprecon")
input_numeric_iprecon.setDescription({ r"controls iterative solution to sparse matrix arising in real-space transformation"+"\n"
                                      +r"to the straight-fieldline angle":
                                       [r"if \c iprecon=0, the preconditioner is `N'",
                                        r"if \c iprecon=1, the preconditioner is `J'",
                                        r"if \c iprecon=2, the preconditioner is `S'"],
                                       r"only relevant if \c Lsparse=2; \see tr00ab() for details": None
                                      })
input_numeric_iprecon.setType("int")
input_numeric_iprecon.setDefaultValue(0)


input_numeric_iotatol = idf.Variable("iotatol")
input_numeric_iotatol.setDescription({r"tolerance required for iterative construction of straight-fieldline angle":
                                      r"only relevant if \c Lsparse.ge.2"})
input_numeric_iotatol.setType("double")
input_numeric_iotatol.setDefaultValue(-1.0)


input_numeric_Lextrap = idf.Variable("Lextrap")
input_numeric_Lextrap.setDescription(r"geometry of innermost interface is defined by extrapolation")
input_numeric_Lextrap.setType("int")
input_numeric_Lextrap.setDefaultValue(0)


input_numeric_Mregular = idf.Variable("Mregular")
input_numeric_Mregular.setDescription({r"maximum regularization factor":
                                       [r"if \c Mregular>=2, then \c regumm \f$_i\f$ = \c Mregular \f$/ 2 \f$ where \c m \f$_i > \f$ \c Mregular"]
                                       })
input_numeric_Mregular.setType("int")
input_numeric_Mregular.setDefaultValue(-1)

# TODO: Lrzaxis
# TODO: Ntoraxis

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

numericlist = idf.Namelist("numericlist")
numericlist.setDescription(r"The namelist \c numericlist controls internal resolution parameters that the user rarely needs to consider.")
numericlist.addidf.Variables(vars_numericlist)

###############################################################################
# locallist 
###############################################################################

input_local_LBeltrami = idf.Variable("LBeltrami")
input_local_LBeltrami.setDescription({r"Control flag for solution of Beltrami equation":
                                      [ r"if \c LBeltrami = 1,3,5 or 7, (SQP) then the Beltrami field in each volume is constructed"+"\n"
                                       +r"by minimizing the magnetic energy with the constraint of fixed helicity;"+"\n"
                                       +r"this is achieved by using sequential quadratic programming as provided by \c E04UFF ."+"\n"
                                       +r"This approach has the benefit (in theory) of robustly constructing minimum energy solutions"+"\n"
                                       +r"when multiple, i.e. bifurcated, solutions exist.",
                                        r"if \c LBeltrami = 2,3,6 or 7, (Newton) then the Beltrami fields are constructed by employing a standard Newton method"+"\n"
                                       +r"for locating an extremum of"+"\n"
                                       +r"\f$F\equiv \int B^2 dv - \mu (\int {\bf A}\cdot{\bf B}dv-{\cal K})\f$,"+"\n"
                                       +r"where \f$\mu\f$ is treated as an independent degree of freedom similar to the parameters describing the vector potential"+"\n"
                                       +r"and \f${\cal K}\f$ is the required value of the helicity;"+"\n"
                                       +r"this is the standard Lagrange multipler approach for locating the constrained minimum;"+"\n"
                                       +r"this method cannot distinguish saddle-type extrema from minima, and which solution that will be obtained depends on the initial guess",
                                        r"if \c LBeltrami = 4,5,6 or 7, (linear) it is assumed that the Beltrami fields are parameterized by \f$\mu\f$;"+"\n"
                                       +r"in this case, it is only required to solve \f$\nabla \times {\bf B} = \mu {\bf B}\f$ which reduces to a system of linear equations;"+"\n"
                                       +r"\f$\mu\f$ may or may not be adjusted iteratively, depending on \c Lconstraint,"+"\n"
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

input_local_Linitgues = idf.Variable("Linitgues")
input_local_Linitgues.setDescription({r"controls how initial guess for Beltrami field is constructed":
                                      [ r"only relevant for routines that require an initial guess for the Beltrami fields, such as the SQP and Newton methods,"+"\n"
                                       +r"or the sparse linear solver",
                                        r"if \c Linitgues=0, the initial guess for the Beltrami field is trivial",
                                        r"if \c Linitgues=1, the initial guess for the Beltrami field is an integrable approximation",
                                        r"if \c Linitgues=2, the initial guess for the Beltrami field is read from file",
                                        r"if \c Linitgues=3, the initial guess for the Beltrami field will be randomized with the maximum \c maxrndgues"]
                                      })
input_local_Linitgues.setType("int")
input_local_Linitgues.setDefaultValue(1)

input_local_Lposdef = idf.Variable("Lposdef")
input_local_Lposdef.setDescription(r"redundant")
input_local_Lposdef.setType("int")
input_local_Lposdef.setDefaultValue(0)

input_local_maxrndgues = idf.Variable("maxrndgues")
input_local_maxrndgues.setDescription(r"the maximum random number of the Beltrami field if \c Linitgues = 3")
input_local_maxrndgues.setType("double")
input_local_maxrndgues.setDefaultValue(1.0)

# TODO: Lmatsolver
# TODO: NiterGMRES
# TODO: epsGMRES
# TODO: LGMRESprec
# TODO: epsILU


vars_locallist = [
        input_local_LBeltrami,
        input_local_Linitgues,
        input_local_Lposdef,
        input_local_maxrndgues
        ]

locallist = idf.Namelist("locallist")
locallist.setDescription({r"The namelist \c locallist controls the construction of the Beltrami fields in each volume.":
                          [ r"The transformation to straight-fieldline coordinates is singular when the rotational-transform of the interfaces is rational;"+"\n"
                           +r"however, the rotational-transform is still well defined."]
                          })
locallist.addidf.Variables(vars_locallist)

###############################################################################
# globallist 
###############################################################################

input_global_Lfindzero = idf.Variable("Lfindzero")
input_global_Lfindzero.setDescription({r"use Newton methods to find zero of force-balance, which is computed by dforce()":
                                       [ r"if \c Lfindzero=0 , then dforce() is called once"+"\n"
                                        +r"to compute the Beltrami fields consistent with the given geometry and constraints",
                                         r"if \c Lfindzero=1 , then call \c C05NDF (uses   function values only), which iteratively calls dforce()",
                                         r"if \c Lfindzero=2 , then call \c C05PDF (uses derivative information), which iteratively calls dforce()"]
                                       })
input_global_Lfindzero.setType("int")
input_global_Lfindzero.setDefaultValue(0)

input_global_escale = idf.Variable("escale")
input_global_escale.setDescription({r"controls the weight factor, \c BBweight, in the force-imbalance harmonics":
                                    [r"\c BBweight(i) \f$\displaystyle \equiv \texttt{opsilon} \times \exp\left[-\texttt{escale} \times (m_i^2+n_i^2) \right]\f$",
                                     r"defined in preset() ; used in dforce()",
                                     r"\sa Eqn.\f$(\ref{eq:forcebalancemn_global})\f$"]
                                    })
input_global_escale.setType("double")
input_global_escale.setDefaultValue(0.0)

input_global_opsilon = idf.Variable("opsilon")
input_global_opsilon.setDescription({r"weighting of force-imbalance":
                                     [r"used in dforce(); \sa Eqn.\f$(\ref{eq:forcebalancemn_global})\f$"]
                                     })
input_global_opsilon.setType("double")
input_global_opsilon.setDefaultValue(1.0)

input_global_pcondense = idf.Variable("pcondense")
input_global_pcondense.setDescription({r"spectral condensation parameter":
                                       [ r"used in preset() to define \c mmpp(i) \f$\equiv m_i^p\f$, where \f$p\equiv \f$ \c pcondense",
                                         r"the angle freedom is exploited to minimize \f$\displaystyle \texttt{epsilon} \sum_{i} m_i^p (R_{i}^2+Z_{i}^2)\f$"+"\n"
                                        +r"with respect to tangential variations in the interface geometry",
                                         r"\sa Eqn.\f$(\ref{eq:spectralbalancemn_global})\f$"]
                                       })
input_global_pcondense.setType("double")
input_global_pcondense.setDefaultValue(2.0)

input_global_epsilon = idf.Variable("epsilon")
input_global_epsilon.setDescription({r"weighting of spectral-width constraint":
                                     [r"used in dforce(); \sa Eqn.\f$(\ref{eq:spectralbalancemn_global})\f$"]
                                     })
input_global_epsilon.setType("double")
input_global_epsilon.setDefaultValue(0.0)

input_global_wpoloidal = idf.Variable("wpoloidal")
input_global_wpoloidal.setDescription([r"\"star-like\" poloidal angle constraint radial exponential factor",
                                       r"used in preset() to construct \c sweight"])
input_global_wpoloidal.setType("double")
input_global_wpoloidal.setDefaultValue(1.0)

input_global_upsilon = idf.Variable("upsilon")
input_global_upsilon.setDescription([r"weighting of \"star-like\" poloidal angle constraint",
                                     r"used in preset() to construct \c sweight"])
input_global_upsilon.setType("double")
input_global_upsilon.setDefaultValue(1.0)

input_global_forcetol = idf.Variable("forcetol")
input_global_forcetol.setDescription({r"required tolerance in force-balance error; only used as an initial check":
                                      [ r"if the initially supplied interfaces are consistent with force-balance to within \c forcetol"+"\n"
                                       +r"then the geometry of the interfaces is not altered",
                                        r"if not, then the geometry of the interfaces is changed in order to bring the configuration into force balance"+"\n"
                                       +r"so that the geometry of interfaces is within \c c05xtol, defined below, of the true solution",
                                        r"to force execution of either \c C05NDF or \c C05PDF, regardless of the initial force imbalance,"+"\n"
                                       +r"set \c forcetol<0"]
                                      })
input_global_forcetol.setType("double")
input_global_forcetol.setDefaultValue(1.0e-10)

input_global_c05xmax = idf.Variable("c05xmax")
input_global_c05xmax.setDescription(r"required tolerance in position, \f${\bf x} \equiv \{ R_{i,v}, Z_{i,v}\}\f$")
input_global_c05xmax.setType("double")
input_global_c05xmax.setDefaultValue(1.0e-6)

input_global_c05xtol = idf.Variable("c05xtol")
input_global_c05xtol.setDescription({r"required tolerance in position, \f${\bf x} \equiv \{ R_{i,v}, Z_{i,v}\}\f$":
                                     [r"used by both \c C05NDF and \c C05PDF; see the NAG documents for further details on how the error is defined",
                                      r"constraint \c c05xtol>0.0"]
                                     })
input_global_c05xtol.setType("double")
input_global_c05xtol.setDefaultValue(1.0e-12)

input_global_c05factor = idf.Variable("c05factor")
input_global_c05factor.setDescription({r"used to control initial step size in \c C05NDF and \c C05PDF":
                                       [r"constraint \c c05factor>0.0",
                                        r"only relevant if \c Lfindzero>0"]
                                       })
input_global_c05factor.setType("double")
input_global_c05factor.setDefaultValue(1.0e-2)

input_global_LreadGF = idf.Variable("LreadGF")
input_global_LreadGF.setDescription({r"read \f$\nabla_{\bf x} {\bf F}\f$ from file \c ext.GF ":
                                     [r"only used if \c Lfindzero=2",
                                      r"only used in newton()"]
                                     })
input_global_LreadGF.setType("boolean")
input_global_LreadGF.setDefaultValue(True)

input_global_mfreeits = idf.Variable("mfreeits")
input_global_mfreeits.setDescription({r"maximum allowed free-boundary iterations":
                                      [r"only used if \c Lfreebound=1",
                                       r"only used in xspech()"]
                                      })
input_global_mfreeits.setType("int")
input_global_mfreeits.setDefaultValue(0)

input_global_bnstol = idf.Variable("bnstol")
input_global_bnstol.setDescription(r"redundant")
input_global_bnstol.setType("double")
input_global_bnstol.setDefaultValue(1.0e-6)

input_global_bnsblend = idf.Variable("bnsblend")
input_global_bnsblend.setDescription(r"redundant")
input_global_bnsblend.setType("double")
input_global_bnsblend.setDefaultValue(0.666)

input_global_gBntol = idf.Variable("gBntol")
input_global_gBntol.setDescription({r"required tolerance in free-boundary iterations":
                                    [r"only used if \c Lfreebound=1",
                                     r"only used in xspech()"]
                                    })
input_global_gBntol.setType("double")
input_global_gBntol.setDefaultValue(1.0e-6)

input_global_gBnbld = idf.Variable("gBnbld")
input_global_gBnbld.setDescription({r"normal blend":
                                    [ r"The \"new\" magnetic field at the computational boundary produced by the plasma currents is updated using a Picard scheme:"+"\n"
                                     +r"\f{eqnarray}{ ({\bf B}\cdot{\bf n})^{j+1} =    \texttt{gBnbld}  \times ({\bf B}\cdot{\bf n})^{j}"+"\n"
                                     +r"                                          + (1-\texttt{gBnbld}) \times ({\bf B}\cdot{\bf n})^{*},"+"\n"
                                     +r"\f}"+"\n"
                                     +r"where \f$j\f$ labels free-boundary iterations, and \f$({\bf B}\cdot{\bf n})^{*}\f$ is computed by virtual casing.",
                                      r"only used if \c Lfreebound=1",
                                      r"only used in xspech()"]
                                    })
input_global_gBnbld.setType("double")
input_global_gBnbld.setDefaultValue(0.666)

input_global_vcasingeps = idf.Variable("vcasingeps")
input_global_vcasingeps.setDescription(r"regularization of Biot-Savart; see bnorml(), casing()")
input_global_vcasingeps.setType("double")
input_global_vcasingeps.setDefaultValue(1.0e-12)

input_global_vcasingtol = idf.Variable("vcasingtol")
input_global_vcasingtol.setDescription(r"accuracy on virtual casing integral; see bnorml(), casing()")
input_global_vcasingtol.setType("double")
input_global_vcasingtol.setDefaultValue(1.0e-8)

input_global_vcasingits = idf.Variable("vcasingits")
input_global_vcasingits.setDescription(r"minimum number of calls to adaptive virtual casing routine; see casing()")
input_global_vcasingits.setType("int")
input_global_vcasingits.setDefaultValue(8)

input_global_vcasingper = idf.Variable("vcasingper")
input_global_vcasingper.setDescription(r"periods of integragion  in adaptive virtual casing routine; see casing()")
input_global_vcasingper.setType("int")
input_global_vcasingper.setDefaultValue(1)

input_global_mcasingcal = idf.Variable("mcasingcal")
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

globallist = idf.Namelist("globallist")
globallist.setDescription({r"The namelist \c globallist controls the search for global force-balance.": None,
                           r"Comments:":
                           [ r'The "force" vector, \f${\bf F}\f$, which is constructed in dforce(), is a combination of pressure-imbalance Fourier harmonics,'+"\n"
                            +r"\f{eqnarray}{ F_{i,v} \equiv [[ p+B^2/2 ]]_{i,v} \times \exp\left[-\texttt{escale}(m_i^2+n_i^2) \right] \times \texttt{opsilon},"+"\n"
                            +r"\label{eq:forcebalancemn_global} \f}"+"\n"
                            +r'and spectral-condensation constraints, \f$I_{i,v}\f$, and the "star-like" angle constraints, \f$S_{i,v,}\f$, (see lforce() for details)'+"\n"
                            +r"\f{eqnarray}{ F_{i,v} \equiv \texttt{epsilon} \times I_{i,v}"+"\n"
                            +r"                           + \texttt{upsilon} \times \left( \psi_v^\omega S_{i,v,1} - \psi_{v+1}^\omega S_{i,v+1,0} \right),"+"\n"
                            +r"\label{eq:spectralbalancemn_global} \f}"+"\n"
                            +r"where \f$\psi_v\equiv\f$ normalized toroidal flux, \c tflux, and \f$\omega\equiv\f$ \c wpoloidal."]
                           })
globallist.addidf.Variables(vars_globallist)

###############################################################################
# diagnosticslist 
###############################################################################

input_diagnostics_odetol = idf.Variable("odetol")
input_diagnostics_odetol.setDescription(r"o.d.e. integration tolerance for all field line tracing routines")
input_diagnostics_odetol.setType("double")
input_diagnostics_odetol.setDefaultValue(1.0e-7)

input_diagnostics_absreq = idf.Variable("absreq")
input_diagnostics_absreq.setDescription(r"redundant")
input_diagnostics_absreq.setType("double")
input_diagnostics_absreq.setDefaultValue(1.0e-8)

input_diagnostics_relreq = idf.Variable("relreq")
input_diagnostics_relreq.setDescription(r"redundant")
input_diagnostics_relreq.setType("double")
input_diagnostics_relreq.setDefaultValue(1.0e-8)

input_diagnostics_absacc = idf.Variable("absacc")
input_diagnostics_absacc.setDescription(r"redundant")
input_diagnostics_absacc.setType("double")
input_diagnostics_absacc.setDefaultValue(1.0e-4)

input_diagnostics_epsr = idf.Variable("epsr")
input_diagnostics_epsr.setDescription(r"redundant")
input_diagnostics_epsr.setType("double")
input_diagnostics_epsr.setDefaultValue(1.0e-8)

input_diagnostics_nPpts = idf.Variable("nPpts")
input_diagnostics_nPpts.setDescription({ r"number of toroidal transits used (per trajectory) in following field lines"+"\n"
                                        +r"for constructing Poincaré plots":
                                         [r"if \c nPpts<1, no Poincaré plot is constructed"]
                                         })
input_diagnostics_nPpts.setType("int")
input_diagnostics_nPpts.setDefaultValue(0)

# TODO: Ppts

input_diagnostics_nPtrj = idf.Variable("nPtrj")
input_diagnostics_nPtrj.setDescription({r"number of trajectories in each annulus to be followed in constructing Poincaré plot":
                                        [ r"if \c nPtrj(l)<0, then \c nPtrj(l) = Ni(l),"+"\n"
                                         +r"where \c Ni(l) is the grid resolution used to construct the Beltrami field in volume \f$l\f$"]
                                        })
input_diagnostics_nPtrj.setType("int")
input_diagnostics_nPtrj.setRank(1)
input_diagnostics_nPtrj.setDefaultValue(-1)
input_diagnostics_nPtrj.setMaximumIndices([r"MNvol+1"])

input_diagnostics_LHevalues = idf.Variable("LHevalues")
input_diagnostics_LHevalues.setDescription(r"to compute eigenvalues of \f$\nabla {\bf F}\f$")
input_diagnostics_LHevalues.setType("boolean")
input_diagnostics_LHevalues.setDefaultValue(False)

input_diagnostics_LHevectors = idf.Variable("LHevectors")
input_diagnostics_LHevectors.setDescription(r"to compute eigenvectors (and also eigenvalues) of \f$\nabla {\bf F}\f$")
input_diagnostics_LHevectors.setType("boolean")
input_diagnostics_LHevectors.setDefaultValue(False)

input_diagnostics_LHmatrix = idf.Variable("LHmatrix")
input_diagnostics_LHmatrix.setDescription(r"to compute and write to file the elements of \f$\nabla {\bf F}\f$")
input_diagnostics_LHmatrix.setType("boolean")
input_diagnostics_LHmatrix.setDefaultValue(False)

input_diagnostics_Lperturbed = idf.Variable("Lperturbed")
input_diagnostics_Lperturbed.setDescription(r"to compute linear, perturbed equilibrium")
input_diagnostics_Lperturbed.setType("int")
input_diagnostics_Lperturbed.setDefaultValue(0)

input_diagnostics_dpp = idf.Variable("dpp")
input_diagnostics_dpp.setDescription(r"perturbed harmonic")
input_diagnostics_dpp.setType("int")
input_diagnostics_dpp.setDefaultValue(-1)

input_diagnostics_dqq = idf.Variable("dqq")
input_diagnostics_dqq.setDescription(r"perturbed harmonic")
input_diagnostics_dqq.setType("int")
input_diagnostics_dqq.setDefaultValue(-1)

# TODO: Lerrortype
# TODO: Ngrid
# TODO: dRZ

input_diagnostics_Lcheck = idf.Variable("Lcheck")
input_diagnostics_Lcheck.setDescription({r"implement various checks":
                                         [ r"if \c Lcheck = 0, no additional check on the calculation is performed",
                                           r"if \c Lcheck = 1, the error in the current, i.e. \f$\nabla\times{\bf B}-\mu{\bf B}\f$ is computed as a post-diagnostic",
                                          {r"if \c Lcheck = 2, the analytic derivatives of the interface transform w.r.t."+"\n"
                                          +r"the helicity multiplier, \f$\mu\f$, and the enclosed poloidal flux, \f$\Delta\psi_p\f$, are compared to a finite-difference estimate":
                                          [ r"only if \c Lconstraint=1",
                                            r"only for \c dspec executable, i.e. must compile with \c DFLAGS=\"-D DEBUG\""]},
                                          {r"if \c Lcheck = 3, the analytic derivatives of the volume w.r.t. interface Fourier harmonic"+"\n"
                                          +r"is compared to a finite-difference estimate":
                                          [ r"must set \c Lfindzero=2",
                                            r"set \c forcetol sufficiently small and set \c LreadGF=F,"+"\n"
                                           +r"so that the matrix of second derivatives is calculated",
                                            r"only for \c dspec executable, i.e. must compile with \c DFLAGS=\"-D DEBUG\""]},
                                          {r"if \c Lcheck = 4, the analytic calculation of the derivatives of the magnetic field, \f$B^2\f$, at the interfaces"+"\n"
                                          +r"is compared to a finite-difference estimate":
                                          [ r"must set \c Lfindzero=2",
                                            r"set \c forcetol sufficiently small",
                                            r"set \c LreadGF=F",
                                            r"only for \c dspec executable, i.e. must compile with \c DFLAGS=\"-D DEBUG\""]},
                                           r"if \c Lcheck = 5, the analytic calculation of the matrix of the derivatives of the force imbalance"+"\n"
                                          +r"is compared to a finite-difference estimate",
                                          {r"if \c Lcheck = 6, the virtual casing calculation is compared to \c xdiagno (Lazerson 2013 \cite y2013_lazerson)":
                                           [ r"the input file for \c xdiagno is written by bnorml()",
                                             r"this provides the Cartesian coordinates on the computational boundary where the virtual casing routine casing()"+"\n"
                                            +r"computes the magnetic field, with the values of the magnetic field being written to the screen for comparison",
                                             r"must set \c Freebound=1, \c Lfindzero>0, \c mfreeits!=0",
                                             r"\c xdiagno must be executed manually"]}
                                          ]
                                         })
input_diagnostics_Lcheck.setType("int")
input_diagnostics_Lcheck.setDefaultValue(0)

input_diagnostics_Ltiming = idf.Variable("Ltiming")
input_diagnostics_Ltiming.setDescription(r"to check timing")
input_diagnostics_Ltiming.setType("boolean")
input_diagnostics_Ltiming.setDefaultValue(False)

input_diagnostics_fudge = idf.Variable("fudge")
input_diagnostics_fudge.setDescription(r"redundant")
input_diagnostics_fudge.setType("double")
input_diagnostics_fudge.setDefaultValue(1.0)

input_diagnostics_scaling = idf.Variable("scaling")
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

diagnosticslist = idf.Namelist("diagnosticslist")
diagnosticslist.setDescription(r"The namelist \c diagnosticslist controls post-processor diagnostics, such as Poincaré  plot resolution, etc.")
diagnosticslist.addidf.Variables(vars_diagnosticslist)

###############################################################################
# screenlist --> not in output file
###############################################################################

###############################################################################
# TODO: output contents
###############################################################################

# TODO: version
# TODO: iterations
# TODO: grid
# TODO: poincare
# TODO: vector_potential
# TODO: output















###############################################################################
# generate Fortran declarations of the input quantities of SPEC
###############################################################################
from adf import indented, toDoc
from genFortran import declareVariable, declareNamelist

# dry-run declaration to determine maximum declaration length for doc indentation
maxLength = 0
for var in vars_physicslist:
    declLen = len(declareVariable(var, attachDescription=False))
    if declLen>maxLength: maxLength = declLen
for var in vars_numericlist:
    declLen = len(declareVariable(var, attachDescription=False))
    if declLen>maxLength: maxLength = declLen
for var in vars_locallist:
    declLen = len(declareVariable(var, attachDescription=False))
    if declLen>maxLength: maxLength = declLen
for var in vars_globallist:
    declLen = len(declareVariable(var, attachDescription=False))
    if declLen>maxLength: maxLength = declLen
for var in vars_diagnosticslist:
    declLen = len(declareVariable(var, attachDescription=False))
    if declLen>maxLength: maxLength = declLen

#print("maximum decl. length: "+str(maxLength))

module_inputlist = "implicit none\n"

# parameters: maximum array dimensions
for param in params_maxDims:
    module_inputlist += declareVariable(param, refDeclLength=maxLength)+"\n"

module_inputlist += "\n"

# input idf.Variables, i.e. namelist contents
module_inputlist += r"""!> \addtogroup grp_global_physicslist physicslist
!> \brief """+toDoc(physicslist.description)+r"""
!> @{
"""
for var in physicslist.idf.Variables:
    module_inputlist += declareVariable(var, refDeclLength=maxLength)+"\n"
module_inputlist += "!> @}\n"


for var in numericlist.idf.Variables:
    module_inputlist += declareVariable(var, refDeclLength=maxLength)+"\n"
module_inputlist += "\n"
for var in locallist.idf.Variables:
    module_inputlist += declareVariable(var, refDeclLength=maxLength)+"\n"
module_inputlist += "\n"
for var in globallist.idf.Variables:
    module_inputlist += declareVariable(var, refDeclLength=maxLength)+"\n"
module_inputlist += "\n"
for var in diagnosticslist.idf.Variables:
    module_inputlist += declareVariable(var, refDeclLength=maxLength)+"\n"

module_inputlist += "\n"

# namelist declarations
module_inputlist += declareNamelist(physicslist)+"\n"
module_inputlist += "\n"
module_inputlist += declareNamelist(numericlist)+"\n"
module_inputlist += "\n"
module_inputlist += declareNamelist(locallist)+"\n"
module_inputlist += "\n"
module_inputlist += declareNamelist(globallist)+"\n"
module_inputlist += "\n"
module_inputlist += declareNamelist(diagnosticslist)+"\n"

with open("../inplst.f90", "w") as f:
    
    f.write(r"""!> @file inplst.f90
!> \brief Input namelists
!> \addtogroup grp_global
!> @{
""")
    
    f.write(r"""module inputlist
!> \brief Input namelists
!> \addtogroup grp_global
!> @{
""")
    f.write(indented(2, module_inputlist, " "))
    f.write("end module inputlist\n")
