#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Jonathan Schilling (jonathan.schilling@ipp.mpg.de)
@author: Caoxiang Zhu (czhu@pppl.gov)
@author: Ksenia Aleynikova (ksenia.aleynikov@ipp.mpg.de)
@author: Zhisong Qu (zhisong.qu@anu.edu.au)
"""

import h5py
import numpy as np  # for isscalar
import os  # for path.abspath
import keyword  # for getting python keywords
SPEC_MAJOR_VERSION = 3.00

# reader class for Stepped Pressure Equilibrium Code output file
# S. Hudson et al., Physics of Plasmas 19, 112502 (2012); doi: 10.1063/1.4765691
class SPECout:
    """
    Class that contains the output of a SPEC calculation.
    Call signature:
        myspec = SPECout(filename) (e.g. myspec=SPECout("/path/to/GxVxxLx.sp.h5") )

    This class contains other post-processing functions in separate file.
    You can use them directly as class attributes, like myspec.plot_pressure().
    """

    # plot functions and others
    from ._plot_kam_surface import plot_kam_surface
    from ._plot_poincare import plot_poincare
    from ._plot_pressure import plot_pressure
    from ._processing import (
        get_grid_and_jacobian_and_metric,
        grid,
        jacobian,
        metric,
        get_B,
        get_modB,
        get_B_covariant
    )
    from ._plot_modB import plot_modB
    from ._plot_iota import plot_iota

    def __init__(self, *args, **kwargs):
        # args[0] should always be the name of a file or an item inside the root object
        # if args[0] is not a filename, kwargs['content'] should be the content to be added
        # as self.`args[0]`

        _content = None
        if kwargs.get("content") == None:
            # assume arg[0] is a filename
            _content = h5py.File(args[0], "r")

            # keep track of which file this object corresponds to
            self.filename = os.path.abspath(args[0])

            # check version and print warning
            try:
                if _content['version'][()][0] < SPEC_MAJOR_VERSION:
                    print("!!!Warning: this python package is used for SPEC!")
            except KeyError:
                print("!!!Warning: you might be not reading a SPEC HDF5 file!")
        elif isinstance(kwargs["content"], h5py.Group):
            _content = kwargs["content"]

        if _content != None:
            for key in _content:
                if isinstance(_content[key], h5py.Group):
                    # recurse into group
                    if key in keyword.kwlist:  # avoid assign python keywords
                        setattr(self, key + "1", SPECout(content=_content[key]))
                    else:
                        setattr(self, key, SPECout(content=_content[key]))
                elif isinstance(_content[key], h5py.Dataset):  # read dataset
                    if (
                        key in keyword.kwlist
                    ):  # add underscore avoiding assigning python keywords
                        setattr(self, key + "_", _content[key][()])
                    else:
                        if (
                            len(_content[key][()]) == 1
                        ):  # if just one element, use the value directly
                            setattr(self, key, _content[key][0])
                        else:
                            setattr(self, key, _content[key][()])

        if isinstance(_content, h5py.File):
            _content.close()

            # make sure that Lrad is always an array
            if np.isscalar(self.input.physics.Lrad):
                self.input.physics.Lrad = np.array([self.input.physics.Lrad])
            # make sure that im always an array
            if np.isscalar(self.output.im):
                self.output.im = np.array([self.output.im])
            # make sure that in_ is always an array
            if np.isscalar(self.output.in_):
                self.output.in_ = np.array([self.output.in_])

            # these define the target dimensions in the radial direction
            Nvol = self.input.physics.Nvol
            if self.input.physics.Lfreebound:
                Nvol += 1
                self.input.physics.Nvol += 1

            Lrad = self.input.physics.Lrad

            # lists for vector potential
            cAte = []
            cAto = []
            cAze = []
            cAzo = []

            # lists for grid
            cRij = []
            cZij = []
            csg = []
            cBR = []
            cBp = []
            cBZ = []

            # split up radial matrix dimension into list of matrices for each of the nested volumes
            start = 0
            for i in range(Nvol):
                # vector potential
                cAte.append(
                    np.atleast_2d(self.vector_potential.Ate)[
                        :, start : start + Lrad[i] + 1
                    ]
                )
                cAto.append(
                    np.atleast_2d(self.vector_potential.Ato)[
                        :, start : start + Lrad[i] + 1
                    ]
                )
                cAze.append(
                    np.atleast_2d(self.vector_potential.Aze)[
                        :, start : start + Lrad[i] + 1
                    ]
                )
                cAzo.append(
                    np.atleast_2d(self.vector_potential.Azo)[
                        :, start : start + Lrad[i] + 1
                    ]
                )

                # grid
                cRij.append(
                    np.atleast_2d(self.grid.Rij)[:, start : start + Lrad[i] + 1]
                )
                cZij.append(
                    np.atleast_2d(self.grid.Zij)[:, start : start + Lrad[i] + 1]
                )
                csg.append(np.atleast_2d(self.grid.sg)[:, start : start + Lrad[i] + 1])
                cBR.append(np.atleast_2d(self.grid.BR)[:, start : start + Lrad[i] + 1])
                cBp.append(np.atleast_2d(self.grid.Bp)[:, start : start + Lrad[i] + 1])
                cBZ.append(np.atleast_2d(self.grid.BZ)[:, start : start + Lrad[i] + 1])

                # move along the merged array dimension
                start = start + Lrad[i] + 1

            # replace original content in data structure
            self.vector_potential.Ate = cAte
            self.vector_potential.Ato = cAto
            self.vector_potential.Aze = cAze
            self.vector_potential.Azo = cAzo

            self.grid.Rij = cRij
            self.grid.Zij = cZij
            self.grid.sg = csg
            self.grid.BR = cBR
            self.grid.Bp = cBp
            self.grid.BZ = cBZ

            if hasattr(self, "poincare"):
                # remove unsuccessful Poincare trajectories
                self.poincare.R = self.poincare.R[self.poincare.success == 1, :, :]
                self.poincare.Z = self.poincare.Z[self.poincare.success == 1, :, :]
                self.poincare.t = self.poincare.t[self.poincare.success == 1, :, :]
                self.poincare.s = self.poincare.s[self.poincare.success == 1, :, :]
                
        return

    # needed for iterating over the contents of the file
    def __iter__(self):
        return iter(self.__dict__)

    def __next__(self):
        return next(self.__dict__)

    # needed for using SPECout with 'with' statement
    def __enter__(self):
        return self

    def __exit__(self, t, v, tb):
        return

    # print a list of items contained in this object
    def inventory(self, prefix=""):
        _prefix = ""
        if prefix != "":
            _prefix = prefix + "/"

        for a in self:
            try:
                # recurse into member
                getattr(self, a).inventory(prefix=_prefix + a)
            except:
                # print item name
                print(_prefix + a)
