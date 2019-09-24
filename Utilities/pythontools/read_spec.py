#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Jonathan Schilling (jonathan.schilling@ipp.mpg.de)
with additions from Caoxiang Zhu to avoid conflicts with Python keywords
"""

import h5py
import numpy as np  # for isscalar
import os           # for path.abspath
import keyword      # for getting python keywords


# reader class for Stepped Pressure Equilibrium Code output file
# S. Hudson et al., Physics of Plasmas 19, 112502 (2012); doi: 10.1063/1.4765691
class SPEC:
    
    # use as s = SPEC(filename), e.g. s=SPEC("ext.h5") or s=SPEC("/path/to/ext.h5")
    def __init__(self, *args, **kwargs):
        # args[0] should always be the name of a file or an item inside the root object
        # if args[0] is not a filename, kwargs['content'] should be the content to be added
        # as self.`args[0]`
        
        _content = None
        if kwargs.get('content') == None:
            # assume arg[0] is a filename
            _content = h5py.File(args[0], "r")
            
            # keep track of which file this object corresponds to
            self.filename = os.path.abspath(args[0])
        elif isinstance(kwargs['content'], h5py.Group):
            _content = kwargs['content']
        
        if (_content != None):
            for key in _content:
                if isinstance(_content[key], h5py.Group):
                    # recurse into group
                    setattr(self, key, SPEC(content=_content[key]))
                elif isinstance(_content[key], h5py.Dataset):  # read dataset
                    if key in keyword.kwlist:  # avoid assign python keywords
                        setattr(self, key + '1', _content[key][()])
                    else:
                        if len(_content[key][()]) == 1:  # if just one element, use the value directly
                            setattr(self, key, _content[key][0])
                        else: 
                            setattr(self, key, _content[key][()])
        
        if isinstance(_content, h5py.File):
            _content.close()
            
            # make sure that Lrad is always an array
            if np.isscalar(self.input.physics.Lrad):
                self.input.physics.Lrad = np.array([self.input.physics.Lrad])
            
            # these define the target dimensions in the radial direction
            Nvol = self.input.physics.Nvol
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
              cAte.append(self.vector_potential.Ate[:, start:start + Lrad[i] + 1])
              cAto.append(self.vector_potential.Ato[:, start:start + Lrad[i] + 1])
              cAze.append(self.vector_potential.Aze[:, start:start + Lrad[i] + 1])
              cAzo.append(self.vector_potential.Azo[:, start:start + Lrad[i] + 1])

              # grid
              cRij.append(self.grid.Rij[:, start:start + Lrad[i] + 1])
              cZij.append(self.grid.Zij[:, start:start + Lrad[i] + 1])
              csg.append(self.grid.sg[:, start:start + Lrad[i] + 1])
              cBR.append(self.grid.BR[:, start:start + Lrad[i] + 1])
              cBp.append(self.grid.Bp[:, start:start + Lrad[i] + 1])
              cBZ.append(self.grid.BZ[:, start:start + Lrad[i] + 1])
            
              # move along the merged array dimension
              start = start + Lrad[i] + 1;
            
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
            
            # remove unsuccessful Poincare trajectories
            self.poincare.R = self.poincare.R[self.poincare.success == 1, :, :]
            self.poincare.Z = self.poincare.Z[self.poincare.success == 1, :, :]
            self.poincare.t = self.poincare.t[self.poincare.success == 1, :, :]
            self.poincare.s = self.poincare.s[self.poincare.success == 1, :, :]
    
    # needed for iterating over the contents of the file
    def __iter__(self):
        return iter(self.__dict__)

    def __next__(self):
        return next(self.__dict__)
    
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
                
