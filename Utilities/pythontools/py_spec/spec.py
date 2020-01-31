#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Jonathan Schilling (jonathan.schilling@ipp.mpg.de)
@author: Caoxiang Zhu (czhu@pppl.gov)
@author: Ksenia Aleynikova (ksenia.aleynikov@ipp.mpg.de)
"""

import h5py
import numpy as np  # for isscalar
import os           # for path.abspath
import keyword      # for getting python keywords

# reader class for Stepped Pressure Equilibrium Code output file
# S. Hudson et al., Physics of Plasmas 19, 112502 (2012); doi: 10.1063/1.4765691
class SPEC:    
    def __init__(self, *args, **kwargs):
        '''Initialization
        use as s = SPEC(filename), e.g. s=SPEC("ext.h5") or s=SPEC("/path/to/ext.h5")

        args[0] should always be the name of a file or an item inside the root object
        if args[0] is not a filename, kwargs['content'] should be the content to be added
        as self.`args[0]`
        '''        
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
              cAte.append(np.atleast_2d(self.vector_potential.Ate)[:, start:start + Lrad[i] + 1])
              cAto.append(np.atleast_2d(self.vector_potential.Ato)[:, start:start + Lrad[i] + 1])
              cAze.append(np.atleast_2d(self.vector_potential.Aze)[:, start:start + Lrad[i] + 1])
              cAzo.append(np.atleast_2d(self.vector_potential.Azo)[:, start:start + Lrad[i] + 1])

              # grid
              cRij.append(np.atleast_2d(self.grid.Rij)[:, start:start + Lrad[i] + 1])
              cZij.append(np.atleast_2d(self.grid.Zij)[:, start:start + Lrad[i] + 1])
              csg.append(np.atleast_2d(self.grid.sg)[:, start:start + Lrad[i] + 1])
              cBR.append(np.atleast_2d(self.grid.BR)[:, start:start + Lrad[i] + 1])
              cBp.append(np.atleast_2d(self.grid.Bp)[:, start:start + Lrad[i] + 1])
              cBZ.append(np.atleast_2d(self.grid.BZ)[:, start:start + Lrad[i] + 1])
            
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
            
            if hasattr(self, 'poincare'):
                # remove unsuccessful Poincare trajectories
            	self.poincare.R = self.poincare.R[self.poincare.success == 1, :, :]
            	self.poincare.Z = self.poincare.Z[self.poincare.success == 1, :, :]
            	self.poincare.t = self.poincare.t[self.poincare.success == 1, :, :]
            	self.poincare.s = self.poincare.s[self.poincare.success == 1, :, :]
        return
    def plot_pressure(self, normalize=True, **kwargs):
        '''Plot stepped pressure profile
        Parameters:
           normalize -- Boolean, True (default). SPEC normalizes the pressure with mu_0. 
                        If False, multiply mu_0 back to pressue.
        '''
        import matplotlib.pyplot as plt
        pressure = self.input.physics.pressure * self.input.physics.pscale
        tflux = self.output.tflux
        if not normalize :
            #  remove  mu_0
            pressure /= (4*np.pi*1.0E-7)
        # get figure and ax data
        if plt.get_fignums():
            fig = plt.gcf()
            ax = plt.gca()
        else :
            fig, ax = plt.subplots()
        # set default plotting parameters
        if kwargs.get('linewidth') == None:
            kwargs.update({'linewidth': 2.0}) # prefer thicker lines
        if kwargs.get('label') == None:
            kwargs.update({'label': 'SPEC_pressure'}) # default label 
        # plots
        for ivol in range(len(pressure)):
            if ivol == 0:
                ax.plot([0, tflux[ivol]],[pressure[ivol], pressure[ivol]],**kwargs)
            else:
                ax.plot([tflux[ivol-1], tflux[ivol]],[pressure[ivol], pressure[ivol]],**kwargs)
                ax.plot([tflux[ivol-1], tflux[ivol-1]],[pressure[ivol-1], pressure[ivol]],**kwargs)
        # Figure properties
        plt.xlabel('Normalized flux',fontsize=20)
        plt.ylabel('Pressure',fontsize=20)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        return
    
    def plot_kam_surface(self, **kwargs):
        # should we use FourSurf class?
        pass
    
    def plot_poincare(self, toroidalIdx=0, prange='full', **kwargs):
        '''Poincare plots
        parameters:
            toroidalIdx -- int, default: 0. The index of toroidal cross-section to be plotted.
            prange -- str, ['full'(default), 'upper', 'lower']. Range of plotted points.
            **kwargs  -- keyword arguments. Matplotlib.pyplot scatter keyword arguments.
            
        return:
            dots -- Matplotlib.pyplot scatter class
        '''
        import matplotlib.pyplot as plt
        import matplotlib.lines as mlines
        # extract slice corresponding to the given toroidal cutplane
        rr = self.poincare.R[:,:,toroidalIdx]
        zz = self.poincare.Z[:,:,toroidalIdx]
        # get current figure or build new one;
        if plt.get_fignums():
            fig = plt.gcf()
            ax = plt.gca()
        else :
            fig, ax = plt.subplots()
        # set default plotting parameters
        # use dots
        if kwargs.get('marker') == None:
            kwargs.update({'marker': '.'}) 
        # use gray color
        if kwargs.get('c') == None:
            kwargs.update({'c': 'gray'}) 
        # make plot depending on the 'range'
        if prange == 'full':
            dots = ax.scatter(rr, zz, **kwargs)
        elif prange == 'upper':
            dots = ax.scatter(rr[zz>=0], zz[zz>=0], **kwargs)
        elif prange == 'lower':
            dots = ax.scatter(rr[zz<=0], zz[zz<=0], **kwargs)
        else :
            raise ValueError("prange should be one of ['full'(default), 'upper', 'lower'].")
        # adjust figure properties
        plt.xlabel('R [m]',fontsize=20)
        plt.ylabel('Z [m]',fontsize=20)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.axis('equal')
        return dots
        
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
                
