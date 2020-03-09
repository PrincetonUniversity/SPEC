########################################
# fourier_surface.py: class that contains the fourier representation of a surface.
# coded by @zhucaoxiang (czhu@pppl.gov)
# adapted by @smiet (csmiet@pppl.gov)
#


class fourier_surface(object):
    import numpy as np
    '''
    toroidal surface in Fourier representation
    R = \sum RBC cos(mu-nv) + RBS sin(mu-nv)
    Z = \sum ZBC cos(mu-nv) + ZBS sin(mu-nv)
    '''
    def __init__(self, xm=[], xn=[], rbc=[], zbs=[], rbs=[], zbc=[]):
        import numpy as np
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
        import numpy as np
        # check if spec_out is in correct format
        #if not isinstance(spec_out, SPEC):
        #    raise TypeError("Invalid type of input data, should be SPEC type.")
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
        import numpy as np
        import xarray as ncdata # read netcdf file
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


    @classmethod
    def read_focus_input(cls, filename, Mpol=9999, Ntor=9999):
        """initialize surface from the FOCUS format input file 'plasma.boundary'

        Parameters:
          filename -- string, path + name to the FOCUS input boundary file
          Mpol -- maximum truncated poloidal mode number (default: 9999)
          Ntol -- maximum truncated toroidal mode number (default: 9999)

        Returns:
          fourier_surface class
        """
        import numpy as np
        with open(filename, 'r') as f:
            line = f.readline() #skip one line
            line = f.readline()
            num = int(line.split()[0]) #harmonics number
            nfp = int(line.split()[1]) #number of field periodicity
            nbn = int(line.split()[2]) #number of Bn harmonics
            xm = []
            xn = []
            rbc = []
            rbs = []
            zbc = []
            zbs = []
            line = f.readline() #skip one line
            line = f.readline() #skip one line
            for i in range(num):
                line = f.readline()
                line_list = line.split()
                n = int(line_list[0])
                m = int(line_list[1])
                if abs(m)>Mpol or abs(n)>Ntor:
                    continue
                xm.append(m)
                xn.append(n)
                rbc.append(float(line_list[2]))
                rbs.append(float(line_list[3]))
                zbc.append(float(line_list[4]))
                zbs.append(float(line_list[5]))
        return cls(xm=np.array(xm), xn=-np.array(xn)*nfp,
                       rbc=np.array(rbc), rbs=np.array(rbs),
                       zbc=np.array(zbc), zbs=np.array(zbs))

    @classmethod
    def read_winding_surfce(cls, filename, Mpol=9999, Ntor=9999):
        """initialize surface from the NESCOIL format input file 'nescin.xxx'

        Parameters:
          filename -- string, path + name to the NESCOIL input boundary file
          Mpol -- maximum truncated poloidal mode number (default: 9999)
          Ntol -- maximum truncated toroidal mode number (default: 9999)

        Returns:
          fourier_surface class
        """
        import numpy as np
        with open(filename, 'r') as f:
            line = ''
            while "phip_edge" not in line:
                line = f.readline()
            line = f.readline()
            nfp = int(line.split()[0])
            #print "nfp:",nfp

            line = ''
            while "Current Surface" not in line:
                line = f.readline()
            line = f.readline()
            line = f.readline()
            #print "Number of Fourier modes in coil surface from nescin file: ",line
            num = int(line)
            xm = []
            xn = []
            rbc = []
            rbs = []
            zbc = []
            zbs = []
            line = f.readline() #skip one line
            line = f.readline() #skip one line
            for i in range(num):
                line = f.readline()
                line_list = line.split()
                m = int(line_list[0])
                n = int(line_list[1])
                if abs(m)>Mpol or abs(n)>Ntor:
                    continue
                xm.append(m)
                xn.append(n)
                rbc.append(float(line_list[2]))
                zbs.append(float(line_list[3]))
                rbs.append(float(line_list[4]))
                zbc.append(float(line_list[5]))
            # NESCOIL uses mu+nv, minus sign is added
            return cls(xm=np.array(xm), xn=-np.array(xn)*nfp,
                       rbc=np.array(rbc), rbs=np.array(rbs),
                       zbc=np.array(zbc), zbs=np.array(zbs))

    def rz(self, theta, zeta, normal=False):
        """ get r,z position of list of (theta, zeta)

        Parameters:
          theta -- float array_like, poloidal angle
          zeta -- float array_like, toroidal angle value
          normal -- logical, calculate the normal vector or not (default: False)
        Returns:
           r, z -- float array_like
           r, z, [rt, zt], [rz, zz] -- if normal
        """
        import numpy as np
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

        if not normal :
            return (r.ravel(), z.ravel())
        else:
            rt = np.matmul( np.reshape(self.xm * self.rbc, (1,-1)), -_sin ) \
               + np.matmul( np.reshape(self.xm * self.rbs, (1,-1)),  _cos )
            zt = np.matmul( np.reshape(self.xm * self.zbc, (1,-1)), -_sin ) \
               + np.matmul( np.reshape(self.xm * self.zbs, (1,-1)),  _cos )

            rz = np.matmul( np.reshape(-self.xn * self.rbc, (1,-1)), -_sin ) \
               + np.matmul( np.reshape(-self.xn * self.rbs, (1,-1)),  _cos )
            zz = np.matmul( np.reshape(-self.xn * self.zbc, (1,-1)), -_sin ) \
               + np.matmul( np.reshape(-self.xn * self.zbs, (1,-1)),  _cos )
            return (r.ravel(), z.ravel(), [rt.ravel(), zt.ravel()], [rz.ravel(), zz.ravel()])

    def xyz(self, theta, zeta, normal=False):
        """ get x,y,z position of list of (theta, zeta)

        Parameters:
          theta -- float array_like, poloidal angle
          zeta -- float array_like, toroidal angle value
          normal -- logical, calculate the normal vector or not (default: False)
        Returns:
           x, y, z -- float array_like
           x, y, z, [nx, ny, nz] -- if normal
        """
        import numpy as np
        data = self.rz(theta, zeta, normal)
        r = data[0]
        z = data[1]
        _sin = np.sin(np.ravel(zeta))
        _cos = np.cos(np.ravel(zeta))
        if not normal:
            return (r*_cos, r*_sin, z)
        else:
            _xt = data[2][0]*_cos #  dx/dtheta
            _yt = data[2][0]*_sin # dy/dtheta
            _zt = data[2][1]      # dz/dtheta
            _xz = data[3][0]*_cos - r*_sin # dx/dzeta
            _yz = data[3][0]*_sin + r*_cos # dy/dzeta
            _zz = data[3][1]               # dz/dzeta
            # n = dr/dz x  dr/dt
            n = np.cross(np.transpose([_xz, _yz, _zz]), np.transpose([_xt, _yt, _zt]))
            return (r*_cos, r*_sin, z, n)

    def areaVolume(self, theta0=0.0, theta1=2*np.pi, zeta0=0.0, zeta1=2*np.pi, \
                   npol=360, ntor=360):
        """ Internel function to get surface area and volume

        Parameters:
          theta0 -- float, starting poloidal angle (default: 0.0)
          theta1 -- float, ending poloidal angle (default: 2*np.pi)
          zeta0 -- float, starting toroidal angle (default: 0.0)
          zeta1 -- float, ending toroidal angle (default: 2*np.pi)
          npol -- integer, number of poloidal discretization points (default: 360)
          ntor -- integer, number of toroidal discretization points (default: 360)
        Returns:
          area -- surface area
          volume -- surface volume
        """
        import numpy as np
        # get mesh data
        _theta = np.linspace(theta0, theta1, npol, endpoint=False)
        _zeta = np.linspace(zeta0, zeta1, ntor, endpoint=False)
        _tv, _zv = np.meshgrid(_theta, _zeta, indexing='ij')
        _x, _y, _z, _n = self.xyz(_tv, _zv, normal=True)
        # calculates the area and volume
        _dt = (theta1-theta0)/npol
        _dz = (zeta1 - zeta0)/ntor
        _nn = np.linalg.norm(_n, axis=1)
        area = np.sum(_nn)*_dt*_dz
        volume = abs(np.sum(_x*_n[:,0]))*_dt*_dz
        return area, volume

    def get_area(self):
        """ Get the surface area and saved in self.area
        More comprehensive options can be found in self._areaVolume()
        Parameters:
           None
        Returns:
           area
        """
        import numpy as np
        self.area, _volume = self._areaVolume()
        return self.area

    def get_volume(self):
        """ Get the surface volume and saved in self.volume
        More comprehensive options can be found in self._areaVolume()
        Parameters:
           None
        Returns:
           volume
        """
        _area, self.volume = self._areaVolume()
        return self.volume

    def plot(self, zeta=0.0, npoints=360, **kwargs):
        """ plot the cross-section at zeta using matplotlib.pyplot

        Parameters:
          zeta -- float, toroidal angle value
          npoints -- integer, number of discretization points (default: 360)
          kwargs -- optional keyword arguments for pyplot
        Returns:
           line class in matplotlib.pyplot
        """
        import matplotlib.pyplot as plt
        import numpy as np
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
        import numpy as np
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
            import matplotlib.pyplot as plt
            from mpl_toolkits.mplot3d import Axes3D
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
            from mayavi import mlab # to overrid plt.mlab
            mlab.mesh(xsurf, ysurf, zsurf, **kwargs)
        else:
            raise ValueError('Invalid engine option {pyplot, mayavi, noplot}')
        return (xsurf, ysurf, zsurf)

    def tovtk(self, vtkname, npol=360, ntor=360, **kwargs):
        """ save surface shape a vtk grid file

        Parameters:
          vtkname -- string, the filename you want to save, final name is 'vtkname.vts'
          npol -- integer, number of poloidal discretization points (default: 360)
          ntor -- integer, number of toroidal discretization points (default: 360)
          kwargs -- optional keyword arguments for saving as pointdata
        Returns:

        """
        import numpy as np
        from pyevtk.hl import gridToVTK # save to binary vtk
        _xx, _yy, _zz = self.plot3d('noplot', zeta0=0.0, zeta1=2*np.pi,
                                    theta0=0.0, theta1=2*np.pi, npol=npol, ntor=ntor)
        _xx = _xx.reshape((1, npol, ntor))
        _yy = _yy.reshape((1, npol, ntor))
        _zz = _zz.reshape((1, npol, ntor))

        if kwargs:
            gridToVTK(vtkname, _xx, _yy, _zz, pointData=kwargs)
        else:
            gridToVTK(vtkname, _xx, _yy, _zz)
        return

    def __del__(self):
        class_name = self.__class__.__name__
