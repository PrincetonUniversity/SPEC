#!/usr/bin/env python
# coding: utf-8
"""
SPEC postprocessing library

# Place holder for math introduction of the SPEC data: <br>
# $R_{ac} cos(m \theta - n \phi)$

@author: Ksenia Aleynikova (ksenia.aleynikov@ipp.mpg.de)
@author: Caoxiang Zhu (czhu@pppl.gov)
@author: Jonathan Schilling (jonathan.schilling@ipp.mpg.de)
"""

import numpy as _np

def get_grid_and_jacobian_and_metric(s,lvol=0,sarr=_np.linspace(1,1,1),tarr=_np.linspace(0,0,1),zarr=_np.linspace(0,0,1)):
    
    Rac,Rbc = s.output.Rbc[lvol:lvol+2]
    Zas,Zbs = s.output.Zbs[lvol:lvol+2]
    
    mn = Rac.size #s.output.mn
    im = s.output.im
    in1 = s.output.in1
    
    sbar    = (sarr+1)/2;
    fac = []
    for j in range(mn):
        if (lvol>1 or im[j]==0):
            fac.append([sbar,0.5*_np.ones(sarr.size)])
        else:
            fac.append([sbar**(im[j]/2.),(im[j]/4.)*sbar**(im[j]/2.-1.)])
    fac=_np.array(fac)
    
    nax = _np.newaxis
    ang_arg = im[:,nax,nax]*tarr[nax,:,nax]-in1[:,nax,nax]*zarr[nax,nax,:]
    cos = _np.cos(ang_arg); sin = _np.sin(ang_arg)
    dR1 = Rac[:,nax] + fac[:,0,:]*(Rbc[:,nax]-Rac[:,nax])
    dZ1 = Zas[:,nax] + fac[:,0,:]*(Zbs[:,nax]-Zas[:,nax])

    Rarr0 = _np.sum(                                           dR1[:,:,nax,nax]*cos[:,nax,:,:],axis=0)
    Zarr0 = _np.sum(                                           dZ1[:,:,nax,nax]*sin[:,nax,:,:],axis=0)

    Rarr1 = _np.sum( fac[:,1,:,nax,nax]*(Rbc[:,nax,nax,nax]-Rac[:,nax,nax,nax])*cos[:,nax,:,:],axis=0)
    Zarr1 = _np.sum( fac[:,1,:,nax,nax]*(Zbs[:,nax,nax,nax]-Zas[:,nax,nax,nax])*sin[:,nax,:,:],axis=0)

    Rarr2 = _np.sum(                        -im[:,nax,nax,nax]*dR1[:,:,nax,nax]*sin[:,nax,:,:],axis=0)
    Zarr2 = _np.sum(                         im[:,nax,nax,nax]*dZ1[:,:,nax,nax]*cos[:,nax,:,:],axis=0)

    Rarr3 = _np.sum(                        in1[:,nax,nax,nax]*dR1[:,:,nax,nax]*sin[:,nax,:,:],axis=0)
    Zarr3 = _np.sum(                       -in1[:,nax,nax,nax]*dZ1[:,:,nax,nax]*cos[:,nax,:,:],axis=0)
    
      
    jacobian = Rarr0*(Rarr2*Zarr1 - Rarr1*Zarr2) # from matlab
    
    g11 = Rarr1**2 + Zarr1**2;                #gss
    g22 = Rarr2**2 + Zarr2**2;                #gtt
    g33 = Rarr0**2 + Rarr3**2 + Zarr3**2;     #gzz
    g12 = Rarr1*Rarr2 + Zarr1*Zarr2;          #gst
    g13 = Rarr1*Rarr3 + Zarr1*Zarr3;          #gsz
    g23 = Rarr2*Rarr3 + Zarr2*Zarr3;          #gtz

    g = _np.array([[g11,g12,g13],
                  [g12,g22,g23],
                  [g13,g23,g33]])
    
    
    return Rarr0, Zarr0, jacobian, g

def grid(s,lvol=0,sarr=_np.linspace(1,1,1),tarr=_np.linspace(0,0,1),zarr=_np.linspace(0,0,1)):
  
    Rarr0, Zarr0, _, _ = get_grid_and_jacobian_and_metric(s,lvol,sarr,tarr,zarr)
    return Rarr0, Zarr0

def jacobian(s,lvol=0,sarr=_np.linspace(1,1,1),tarr=_np.linspace(0,0,1),zarr=_np.linspace(0,0,1)):

    _, _, jacobian, _ = get_grid_and_jacobian_and_metric(s,lvol,sarr,tarr,zarr)
    return jacobian

def metric(s,lvol=0,sarr=_np.linspace(1,1,1),tarr=_np.linspace(0,0,1),zarr=_np.linspace(0,0,1)):
    
    _, _, _, g = get_grid_and_jacobian_and_metric(s,lvol,sarr,tarr,zarr)
    return g

def get_B(s,lvol=0,jacobian=None,sarr=_np.linspace(0,0,1),tarr=_np.linspace(0,0,1),zarr=_np.linspace(0,0,1)):
    
    if jacobian is None:
        R, Z, jacobian, g = get_grid_and_jacobian_and_metric(s,lvol=lvol,sarr=sarr,tarr=tarr,zarr=zarr)

    #Lrad = s.i_nput.physics.Lrad[lvol]
    Ate  = s.vector_potential.Ate[lvol]
    Aze  = s.vector_potential.Aze[lvol]
    Ato  = s.vector_potential.Ato[lvol]
    Azo  = s.vector_potential.Azo[lvol]
    
    mn = Ate.shape[0]
    im = s.output.im
    in1 = s.output.in1
    
    fac = []
    sbar    = (sarr+1)/2;
    for j in range(mn):
        if (lvol>1 or im[j]==0):
            fac.append([_np.ones(sarr.size),_np.zeros(sarr.size)])
        else:
            fac.append([sbar**(im[j]/2.),(im[j]/4.)*sbar**(im[j]/2.-1.)])
    fac=_np.array(fac)
    
    import numpy.polynomial.chebyshev as Cheb
    
    nax = _np.newaxis
                #[mn,it,iz]
    ang_arg = im[:,nax,nax]*tarr[nax,:,nax]-in1[:,nax,nax]*zarr[nax,nax,:]
    cosa = _np.cos(ang_arg); sina = _np.sin(ang_arg)
        
            #Ch ,mn,t ,z                                  
    c = ((im[nax,:,nax,nax]*Azo.T[:,:,nax,nax] + in1[nax,:,nax,nax]*Ato.T[:,:,nax,nax])*cosa[nax,:,:,:] 
        -(im[nax,:,nax,nax]*Aze.T[:,:,nax,nax] + in1[nax,:,nax,nax]*Ate.T[:,:,nax,nax])*sina[nax,:,:,:])

    Bs = _np.rollaxis(_np.sum( fac[:,0,nax,nax,:]*Cheb.chebval(sarr, c) , axis=0),2)

    c1 = Aze.T[:,:,nax,nax]*cosa[nax,:,:,:] + Azo.T[:,:,nax,nax]*sina[nax,:,:,:]
    Bt = _np.rollaxis(_np.sum(-fac[:,0,nax,nax,:]*Cheb.chebval(sarr, Cheb.chebder(c1))
                            -fac[:,1,nax,nax,:]*Cheb.chebval(sarr, c1), axis=0),2)
    
    c2 = Ate.T[:,:,nax,nax]*cosa[nax,:,:,:] + Ato.T[:,:,nax,nax]*sina[nax,:,:,:]
    Bz = _np.rollaxis(_np.sum( fac[:,0,nax,nax,:]*Cheb.chebval(sarr, Cheb.chebder(c2))
                            +fac[:,1,nax,nax,:]*Cheb.chebval(sarr, c2), axis=0),2)
    
    Bcontrav = _np.array([Bs,Bt,Bz])/jacobian
    return Bcontrav

#Bcontrav = get_B(s,lvol=lvol,jacobian=jacobian,sarr=sarr,tarr=tarr,zarr=zarr)


def get_modB(Bcontrav,g):
    """Input - Bcontrav has to come from get_B function
    """
    modB = _np.sqrt(_np.einsum('iabc,jiabc,jabc->abc',Bcontrav,g,Bcontrav))
    return modB






#    #!/usr/bin/env python
#    # coding: utf-8
#
#    # # Plot of |B| in (R,Z) cross-section by K. Aleynikova, ksenia.aleynikova@ipp.mpg.de, 2019
#
#    # In[86]:
#
#
#    """
#    @author Ksenia Aleynikova ksenia.aleynikov@ipp.mpg.de
#    """
#
#    import numpy as np
#    from read_spec import SPEC
#
#    import SPEC_postprocessing_lib
#
#
#
#    s = SPEC('G3V01L0Fi.002.ksu1.h5')
#
#
#
#    sarr    = np.linspace(-0.999,1,4)
#    tarr    = np.linspace(0,2*np.pi,5)
#    zarr    = np.linspace(0,0,6)
#    lvol    = 0
#
#
#
#    R, Z, jacobian, g, jacobian_from_metric = SPEC_lib.get_grid_and_jacobian_and_metric(s,lvol=lvol,sarr=sarr,tarr=tarr,zarr=zarr)
#    Bcontrav          = SPEC_lib.get_B(s,lvol=lvol,jacobian=jacobian,sarr=sarr,tarr=tarr,zarr=zarr)
#
#
#
#    modB              = SPEC_lib.get_modB(Bcontrav,g)
#
#
#
#    import matplotlib.pyplot as plt
#
#    fig = plt.figure(figsize=(3*2,5*2))
#    ax = fig.add_subplot(111)
#    plot = ax.pcolormesh(R[:,:,0],Z[:,:,0],modB[:,:,0], vmin=2.8, vmax=3.5)#,linewidths=0.01/3,edgecolors='y')
#    ax.set_aspect('equal')
#    ax.set_xlabel(r"$R$")
#    ax.set_ylabel(r"$Z$")
#    cbar = fig.colorbar(plot, ax=ax)
#    cbar.set_label(r"$|\vec B|$", rotation=0)
#    #plt.show()
#
#
#
#
#
#    #!/usr/bin/env python3
#    # -*- coding: utf-8 -*-
#    """
#    @author: Caoxiang Zhu (czhu@pppl.gov)
#    copied from surface.py at https://github.com/zhucaoxiang/CoilPy
#    """
#
#    import numpy as np
#    import matplotlib.pyplot as plt
#    from mpl_toolkits.mplot3d import Axes3D
#    import xarray as ncdata # read netcdf file
#    from read_spec import SPEC
#
#
#    ########################################################################
#
#    class FourSurf(object):
#        '''
#        toroidal surface in Fourier representation
#        R = \sum RBC cos(mu-nv) + RBS sin(mu-nv)
#        Z = \sum ZBC cos(mu-nv) + ZBS sin(mu-nv)
#        '''
#        def __init__(self, xm=[], xn=[], rbc=[], zbs=[], rbs=[], zbc=[]):
#            """Initialization with Fourier harmonics.
#            
#            Parameters:
#              xm -- list or numpy array, array of m index (default: [])
#              xn -- list or numpy array, array of n index (default: [])
#              rbc -- list or numpy array, array of radial cosine harmonics (default: [])
#              zbs -- list or numpy array, array of z sine harmonics (default: [])   
#              rbs -- list or numpy array, array of radial sine harmonics (default: [])
#              zbc -- list or numpy array, array of z cosine harmonics (default: [])      
#            
#            """
#            self.xm  = np.atleast_1d(xm)
#            self.xn  = np.atleast_1d(xn)
#            self.rbc = np.atleast_1d(rbc)
#            self.rbs = np.atleast_1d(rbs)
#            self.zbc = np.atleast_1d(zbc)
#            self.zbs = np.atleast_1d(zbs)
#            self.mn = len(self.xn)
#            return
#
#        @classmethod
#        def read_spec_output(cls, spec_out, ns=-1):
#            """initialize surface from the ns-th interface SPEC output 
#            
#            Parameters:
#              spec_out -- SPEC class, SPEC hdf5 results
#              ns -- integer, the index of SPEC interface (default: -1)
#            
#            Returns:
#              fourier_surface class
#            """
#            # check if spec_out is in correct format
#            if not isinstance(spec_out, SPEC):
#                raise TypeError("Invalid type of input data, should be SPEC type.")
#            # get required data
#            xm = spec_out.output.im
#            xn = spec_out.output.in1
#            rbc = spec_out.output.Rbc[ns,:]
#            zbs = spec_out.output.Zbs[ns,:]
#            if spec_out.input.physics.Istellsym:
#                # stellarator symmetry enforced
#                rbs = np.zeros_like(rbc)
#                zbc = np.zeros_like(rbc)
#            else:
#                rbs = spec_out.output.Rbs[ns,:]
#                zbc = spec_out.output.Zbc[ns,:]
#            return cls(xm=xm, xn=xn, rbc=rbc, rbs=rbs, zbc=zbc, zbs=zbs)
#
#        @classmethod
#        def read_vmec_output(cls, woutfile, ns=-1):
#            """initialize surface from the ns-th interface SPEC output 
#            
#            Parameters:
#              woutfile -- string, path + name to the wout file from VMEC output
#              ns -- integer, the index of VMEC nested flux surfaces (default: -1)
#            
#            Returns:
#              fourier_surface class
#            """
#            vmec = ncdata.open_dataset(woutfile)
#            xm = vmec['xm'].values
#            xn = vmec['xn'].values
#            rmnc = vmec['rmnc'].values
#            zmns = vmec['zmns'].values
#            rbc = rmnc[ns,:]
#            zbs = zmns[ns,:]
#
#            if vmec['lasym__logical__'].values:
#                # stellarator symmetry enforced
#                zmnc = vmec['zmnc'].values
#                rmns = vmec['rmns'].values
#                rbs = rmns[ns,:]
#                zbc = zmnc[ns,:]
#            else :
#                rbs = np.zeros_like(rbc)
#                zbc = np.zeros_like(rbc)
#            return cls(xm=xm, xn=xn, rbc=rbc, rbs=rbs, zbc=zbc, zbs=zbs)
#
#        def rz(self, theta, zeta):
#            """ get r,z position of list of (theta, zeta)
#            
#            Parameters:
#              theta -- float array_like, poloidal angle
#              zeta -- float array_like, toroidal angle value
#
#            Returns:
#               r, z -- float array_like
#            """
#            assert len(np.atleast_1d(theta)) == len(np.atleast_1d(zeta)), "theta, zeta should be equal size"
#            # mt - nz (in matrix)
#            _mtnz = np.matmul( np.reshape(self.xm, (-1,1)), np.reshape(theta, (1,-1)) ) \
#                  - np.matmul( np.reshape(self.xn, (-1,1)), np.reshape( zeta, (1,-1)) ) 
#            _cos = np.cos(_mtnz)
#            _sin = np.sin(_mtnz)
#
#            r = np.matmul( np.reshape(self.rbc, (1,-1)), _cos ) \
#              + np.matmul( np.reshape(self.rbs, (1,-1)), _sin )
#
#            z = np.matmul( np.reshape(self.zbc, (1,-1)), _cos ) \
#              + np.matmul( np.reshape(self.zbs, (1,-1)), _sin )
#            return (r.ravel(), z.ravel())
#
#        def xyz(self, theta, zeta):
#            """ get x,y,z position of list of (theta, zeta)
#            
#            Parameters:
#              theta -- float array_like, poloidal angle
#              zeta -- float array_like, toroidal angle value
#
#            Returns:
#               x, y, z -- float array_like
#            """
#            r, z = self.rz(theta, zeta)
#            return (r*np.cos(np.ravel(zeta)), r*np.sin(np.ravel(zeta)), z)
#
#        def plot(self, zeta=0.0, npoints=360, **kwargs):
#            """ plot the cross-section at zeta using matplotlib.pyplot
#            
#            Parameters:       
#              zeta -- float, toroidal angle value
#              npoints -- integer, number of discretization points (default: 360)
#              kwargs -- optional keyword arguments for pyplot
#
#            Returns:
#               line class in matplotlib.pyplot
#            """
#            # get figure and ax data
#            if plt.get_fignums():
#                fig = plt.gcf()
#                ax = plt.gca()
#            else :
#                fig, ax = plt.subplots()
#            # set default plotting parameters
#            if kwargs.get('linewidth') == None:
#                kwargs.update({'linewidth': 2.0}) # prefer thicker lines
#            if kwargs.get('label') == None:
#                kwargs.update({'label': 'toroidal surface'}) # default label 
#            # get (r,z) data
#            _r, _z = self.rz( np.linspace(0, 2*np.pi, npoints), zeta*np.ones(npoints) )
#            line = ax.plot(_r, _z, **kwargs)
#            plt.axis('equal')
#            plt.xlabel('R [m]',fontsize=20)
#            plt.ylabel('Z [m]',fontsize=20)
#            plt.xticks(fontsize=16)
#            plt.yticks(fontsize=16)
#            return line
#
#        def plot3d(self, engine='pyplot', theta0=0.0, theta1=2*np.pi, zeta0=0.0, zeta1=2*np.pi, \
#                       npol=360, ntor=360, **kwargs):
#            """ plot 3D shape of the surface
#            
#            Parameters: 
#              engine -- string, plotting engine {'pyplot' (default), 'mayavi', 'noplot'}
#              theta0 -- float, starting poloidal angle (default: 0.0)
#              theta1 -- float, ending poloidal angle (default: 2*np.pi)
#              zeta0 -- float, starting toroidal angle (default: 0.0)
#              zeta1 -- float, ending toroidal angle (default: 2*np.pi)
#              npol -- integer, number of poloidal discretization points (default: 360)
#              ntor -- integer, number of toroidal discretization points (default: 360)
#              kwargs -- optional keyword arguments for plotting
#
#            Returns:
#               xsurf, ysurf, zsurf -- arrays of x,y,z coordinates on the surface
#            """
#            # get mesh data
#            _theta = np.linspace(theta0, theta1, npol)
#            _zeta = np.linspace(zeta0, zeta1, ntor)
#            _tv, _zv = np.meshgrid(_theta, _zeta, indexing='ij')
#            _x, _y, _z = self.xyz(_tv, _zv)
#            xsurf = np.reshape(_x, (npol, ntor))
#            ysurf = np.reshape(_y, (npol, ntor))
#            zsurf = np.reshape(_z, (npol, ntor))
#            if engine == 'noplot':
#                # just return xyz data
#                pass
#            elif engine == 'pyplot':
#                # plot in matplotlib.pyplot
#                if plt.get_fignums():
#                    fig = plt.gcf()
#                    ax = plt.gca()
#                else :
#                    fig = plt.figure()
#                    ax = fig.add_subplot(111, projection='3d')
#                ax.plot_surface(xsurf, ysurf, zsurf, **kwargs)
#            elif engine == 'mayavi':
#                # plot 3D surface in mayavi.mlab
#                mlab.mesh(xsurf, ysurf, zsurf, **kwargs)
#            else:
#                raise ValueError('Invalid engine option {pyplot, mayavi, noplot}')
#            return (xsurf, ysurf, zsurf)
#
#        def __del__(self):
#            class_name = self.__class__.__name__
#            
#    ########################################################################
#
#    def plot_spec_kam(spec, ns='all', zeta=0.0, **kwargs):
#        """ plot SPEC KAM surfaces from hdf5 file
#        
#        Parameters:
#          spec -- SPEC class, hdf5 file read by using SPEC class
#          ns -- integer or array or 'all', index of KAM surfaces, 
#                if 'all', plotting all the surfaces (default: 'all')
#          zeta -- float, toroidal angle value
#          kwargs -- optional keyword arguments for pyplot
#
#        Returns:
#           surfs -- array of surfaces in FourSurf type
#        """
#        # check if spec_out is in correct format
#        #if not isinstance(spec, SPEC):
#        #    raise TypeError("Invalid type of input data, should be SPEC type.")
#        # get surface index
#        if ns == 'all':
#            # get the number of all surfaces
#            Nsurf = spec.input.physics.Nvol + 1 + spec.input.physics.Lfreebound 
#            ns = np.arange(Nsurf)
#        else:
#            ns = np.atleast_1d(ns)
#            Nsurf = len(ns)
#        # construct surface
#        surfs = []
#        for i in ns:
#            _surf = FourSurf.read_spec_output(spec, i)
#            if i==0:
#                # plot axis
#                _r, _z = _surf.rz(0.0, zeta)
#                plt.scatter(_r, _z, **kwargs)
#            else:
#                _surf.plot(zeta=zeta, **kwargs)
#            surfs.append(_surf)
#        return surfs


# plot the Poincare data
def plot_poincare(s, toroidalIdx=0, plt=None):
    
    # extract slice corresponding to the given toroidal cutplane
    pltR=s.poincare.R[:,:,toroidalIdx]
    pltZ=s.poincare.Z[:,:,toroidalIdx]
    
    # do the plot
    plt.plot(pltR, pltZ, 'k.', markersize=1)

