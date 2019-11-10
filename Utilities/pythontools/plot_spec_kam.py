#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Caoxiang Zhu (czhu@pppl.gov)
copied from surface.py at https://github.com/zhucaoxiang/CoilPy
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import xarray as ncdata # read netcdf file
from read_spec import SPEC


########################################################################

class FourSurf(object):
    '''
    toroidal surface in Fourier representation
    R = \sum RBC cos(mu-nv) + RBS sin(mu-nv)
    Z = \sum ZBC cos(mu-nv) + ZBS sin(mu-nv)
    '''
    def __init__(self, xm=[], xn=[], rbc=[], zbs=[], rbs=[], zbc=[]):
        """Initialization with Fourier harmonics.
        
        Parameters:
          xm -- list or numpy array, array of m index (default: [])
          xn -- list or numpy array, array of n index (default: [])
          rbc -- list or numpy array, array of radial cosine harmonics (default: [])
          zbs -- list or numpy array, array of z sine harmonics (default: [])   
          rbs -- list or numpy array, array of radial sine harmonics (default: [])
          zbc -- list or numpy array, array of z cosine harmonics (default: [])      
        
        """
        self.xm  = np.atleast_1d(xm)
        self.xn  = np.atleast_1d(xn)
        self.rbc = np.atleast_1d(rbc)
        self.rbs = np.atleast_1d(rbs)
        self.zbc = np.atleast_1d(zbc)
        self.zbs = np.atleast_1d(zbs)
        self.mn = len(self.xn)
        return

    @classmethod
    def read_spec_output(cls, spec_out, ns=-1):
        """initialize surface from the ns-th interface SPEC output 
        
        Parameters:
          spec_out -- SPEC class, SPEC hdf5 results
          ns -- integer, the index of SPEC interface (default: -1)
        
        Returns:
          fourier_surface class
        """
        # check if spec_out is in correct format
        if not isinstance(spec_out, SPEC):
            raise TypeError("Invalid type of input data, should be SPEC type.")
        # get required data
        xm = spec_out.output.im
        xn = spec_out.output.in1
        rbc = spec_out.output.Rbc[ns,:]
        zbs = spec_out.output.Zbs[ns,:]
        if spec_out.input.physics.Istellsym:
            # stellarator symmetry enforced
            rbs = np.zeros_like(rbc)
            zbc = np.zeros_like(rbc)
        else:
            rbs = spec_out.output.Rbs[ns,:]
            zbc = spec_out.output.Zbc[ns,:]
        return cls(xm=xm, xn=xn, rbc=rbc, rbs=rbs, zbc=zbc, zbs=zbs)

    @classmethod
    def read_vmec_output(cls, woutfile, ns=-1):
        """initialize surface from the ns-th interface SPEC output 
        
        Parameters:
          woutfile -- string, path + name to the wout file from VMEC output
          ns -- integer, the index of VMEC nested flux surfaces (default: -1)
        
        Returns:
          fourier_surface class
        """
        vmec = ncdata.open_dataset(woutfile)
        xm = vmec['xm'].values
        xn = vmec['xn'].values
        rmnc = vmec['rmnc'].values
        zmns = vmec['zmns'].values
        rbc = rmnc[ns,:]
        zbs = zmns[ns,:]

        if vmec['lasym__logical__'].values:
            # stellarator symmetry enforced
            zmnc = vmec['zmnc'].values
            rmns = vmec['rmns'].values
            rbs = rmns[ns,:]
            zbc = zmnc[ns,:]
        else :
            rbs = np.zeros_like(rbc)
            zbc = np.zeros_like(rbc)
        return cls(xm=xm, xn=xn, rbc=rbc, rbs=rbs, zbc=zbc, zbs=zbs)

    def rz(self, theta, zeta):
        """ get r,z position of list of (theta, zeta)
        
        Parameters:
          theta -- float array_like, poloidal angle
          zeta -- float array_like, toroidal angle value

        Returns:
           r, z -- float array_like
        """
        assert len(np.atleast_1d(theta)) == len(np.atleast_1d(zeta)), "theta, zeta should be equal size"
        # mt - nz (in matrix)
        _mtnz = np.matmul( np.reshape(self.xm, (-1,1)), np.reshape(theta, (1,-1)) ) \
              - np.matmul( np.reshape(self.xn, (-1,1)), np.reshape( zeta, (1,-1)) ) 
        _cos = np.cos(_mtnz)
        _sin = np.sin(_mtnz)

        r = np.matmul( np.reshape(self.rbc, (1,-1)), _cos ) \
          + np.matmul( np.reshape(self.rbs, (1,-1)), _sin )

        z = np.matmul( np.reshape(self.zbc, (1,-1)), _cos ) \
          + np.matmul( np.reshape(self.zbs, (1,-1)), _sin )
        return (r.ravel(), z.ravel())

    def xyz(self, theta, zeta):
        """ get x,y,z position of list of (theta, zeta)
        
        Parameters:
          theta -- float array_like, poloidal angle
          zeta -- float array_like, toroidal angle value

        Returns:
           x, y, z -- float array_like
        """
        r, z = self.rz(theta, zeta)
        return (r*np.cos(np.ravel(zeta)), r*np.sin(np.ravel(zeta)), z)

    def plot(self, zeta=0.0, npoints=360, **kwargs):
        """ plot the cross-section at zeta using matplotlib.pyplot
        
        Parameters:       
          zeta -- float, toroidal angle value
          npoints -- integer, number of discretization points (default: 360)
          kwargs -- optional keyword arguments for pyplot

        Returns:
           line class in matplotlib.pyplot
        """
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
            kwargs.update({'label': 'toroidal surface'}) # default label 
        # get (r,z) data
        _r, _z = self.rz( np.linspace(0, 2*np.pi, npoints), zeta*np.ones(npoints) )
        line = ax.plot(_r, _z, **kwargs)
        plt.axis('equal')
        plt.xlabel('R [m]',fontsize=20)
        plt.ylabel('Z [m]',fontsize=20)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        return line

    def plot3d(self, engine='pyplot', theta0=0.0, theta1=2*np.pi, zeta0=0.0, zeta1=2*np.pi, \
                   npol=360, ntor=360, **kwargs):
        """ plot 3D shape of the surface
        
        Parameters: 
          engine -- string, plotting engine {'pyplot' (default), 'mayavi', 'noplot'}
          theta0 -- float, starting poloidal angle (default: 0.0)
          theta1 -- float, ending poloidal angle (default: 2*np.pi)
          zeta0 -- float, starting toroidal angle (default: 0.0)
          zeta1 -- float, ending toroidal angle (default: 2*np.pi)
          npol -- integer, number of poloidal discretization points (default: 360)
          ntor -- integer, number of toroidal discretization points (default: 360)
          kwargs -- optional keyword arguments for plotting

        Returns:
           xsurf, ysurf, zsurf -- arrays of x,y,z coordinates on the surface
        """
        # get mesh data
        _theta = np.linspace(theta0, theta1, npol)
        _zeta = np.linspace(zeta0, zeta1, ntor)
        _tv, _zv = np.meshgrid(_theta, _zeta, indexing='ij')
        _x, _y, _z = self.xyz(_tv, _zv)
        xsurf = np.reshape(_x, (npol, ntor))
        ysurf = np.reshape(_y, (npol, ntor))
        zsurf = np.reshape(_z, (npol, ntor))
        if engine == 'noplot':
            # just return xyz data
            pass
        elif engine == 'pyplot':
            # plot in matplotlib.pyplot
            if plt.get_fignums():
                fig = plt.gcf()
                ax = plt.gca()
            else :
                fig = plt.figure()
                ax = fig.add_subplot(111, projection='3d')
            ax.plot_surface(xsurf, ysurf, zsurf, **kwargs)
        elif engine == 'mayavi':
            # plot 3D surface in mayavi.mlab
            mlab.mesh(xsurf, ysurf, zsurf, **kwargs)
        else:
            raise ValueError('Invalid engine option {pyplot, mayavi, noplot}')
        return (xsurf, ysurf, zsurf)

    def __del__(self):
        class_name = self.__class__.__name__
        
########################################################################

def plot_spec_kam(spec, ns='all', zeta=0.0, **kwargs):
    """ plot SPEC KAM surfaces from hdf5 file
    
    Parameters:
      spec -- SPEC class, hdf5 file read by using SPEC class
      ns -- integer or array or 'all', index of KAM surfaces, 
            if 'all', plotting all the surfaces (default: 'all')
      zeta -- float, toroidal angle value
      kwargs -- optional keyword arguments for pyplot

    Returns:
       surfs -- array of surfaces in FourSurf type
    """
    # check if spec_out is in correct format
    #if not isinstance(spec, SPEC):
    #    raise TypeError("Invalid type of input data, should be SPEC type.")
    # get surface index
    if ns == 'all':
        # get the number of all surfaces
        Nsurf = spec.input.physics.Nvol + 1 + spec.input.physics.Lfreebound 
        ns = np.arange(Nsurf)
    else:
        ns = np.atleast_1d(ns)
        Nsurf = len(ns)
    # construct surface
    surfs = []
    for i in ns:
        _surf = FourSurf.read_spec_output(spec, i)
        if i==0:
            # plot axis
            _r, _z = _surf.rz(0.0, zeta)
            plt.scatter(_r, _z, **kwargs)
        else:
            _surf.plot(zeta=zeta, **kwargs)
        surfs.append(_surf)
    return surfs
