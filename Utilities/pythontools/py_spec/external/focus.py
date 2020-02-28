###############################################################################
# Focus.py: read in the output from FOCUS runs.
#

from .hdf5 import HDF5
from .misc import get_figure, map_matrix
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter


# the FOCUS class is a wrapper over the HDF5 class defined in this same folder
class FOCUS(HDF5):
    """FOCUS output hdf5 file

    """
    # initialization, test = FOCUSHDF5('focus_test.h5')
    def __init__(self, filename, periodic=False, **kwargs):
        """ Initialization

        Keyword arguments:
            filenmae -- string, path and name to FOCUS output hdf5 file,
                        usually in the format of 'focus_*.h5'
            periodic -- logical, map all 2D surface data automatically,
                        (default: True)
        """
        # read data
        super().__init__(filename)
        # print version
        try:
            print("FOCUS version: "+self.version)
        except AttributeError:
            print(filename + " is not a valid FOCUS output. Please check.")
            raise
        # add additional colume and row for plotting
        if periodic:
            self.xsurf = map_matrix(self.xsurf)
            self.ysurf = map_matrix(self.ysurf)
            self.zsurf = map_matrix(self.zsurf)
            self.nx = map_matrix(self.nx)
            self.ny = map_matrix(self.ny)
            self.nz = map_matrix(self.nz)
            self.nn = map_matrix(self.nn)
            self.Bx = map_matrix(self.Bx)
            self.By = map_matrix(self.By)
            self.Bz = map_matrix(self.Bz)
            self.Bn = map_matrix(self.Bn)
            self.plas_Bn = map_matrix(self.plas_Bn)
    # convergence plot
    def convergence(self, term='bnorm', iteration=True, axes=None, **kwargs):
        # get figure
        fig, axes = get_figure(axes)
        # set default plotting parameters
        kwargs['linewidth'] = kwargs.get('linewidth', 2.5) # line width
        #kwargs['marker'] = kwargs.get('marker', 'o') # extent
        # get iteration data
        if iteration:
            abscissa = 1+ np.arange(self.iout)
            _xlabel = 'iteration'
        else :
            abscissa = self.evolution[0,:] # be careful; DF is not saving wall-time
            _xlabel = 'wall time [Second]'
        # plot data
        if term.lower() == 'chi':
            data = self.evolution[1,:]
            kwargs['label'] = kwargs.get('label', r'$\chi^2$') # line label
        elif term.lower() == 'gradient':
            data = self.evolution[2,:]
            kwargs['label'] = kwargs.get('label', r'$|d \chi^2 / d {\bf X}|$') # line label
        elif term.lower() == 'bnorm':
            data = self.evolution[3,:]
            kwargs['label'] = kwargs.get('label', r'$f_{B_n}$') # line label
        elif term.lower() == 'bharm':
            data = self.evolution[4,:]
            kwargs['label'] = kwargs.get('label', r'$f_{B_{mn}}$') # line label
        elif term.lower() == 'tflux':
            data = self.evolution[5,:]
            kwargs['label'] = kwargs.get('label', r'$f_{\Psi}$') # line label
        elif term.lower() == 'ttlen':
            data = self.evolution[6,:]
            kwargs['label'] = kwargs.get('label', r'$f_L$') # line label
        elif term.lower() == 'cssep':
            data = self.evolution[7,:]
            kwargs['label'] = kwargs.get('label', r'$f_{CS}$') # line label
        elif term.lower() == 'all':
            lines = []
            lines.append(self.convergence(term='chi', iteration=iteration, axes=axes,
                                    linestyle='-', color='k', **kwargs))
            lines.append(self.convergence(term='gradient', iteration=iteration, axes=axes,
                                    linestyle='--', color='k', **kwargs))
            lines.append(self.convergence(term='bnorm', iteration=iteration, axes=axes,
                                    linestyle='-.', color='r', **kwargs))
            lines.append(self.convergence(term='tflux', iteration=iteration, axes=axes,
                                    linestyle='-.', color='g', **kwargs))
            lines.append(self.convergence(term='bharm', iteration=iteration, axes=axes,
                                    linestyle='-.', color='b', **kwargs))
            lines.append(self.convergence(term='ttlen', iteration=iteration, axes=axes,
                                    linestyle='-.', color='c', **kwargs))
            lines.append(self.convergence(term='cssep', iteration=iteration, axes=axes,
                                    linestyle='-.', color='m', **kwargs))
            #fig.legend(loc='upper right', frameon=False, ncol=2, prop={'size':16})
            plt.legend()
            return lines
        else :
            raise ValueError('unsupported option for term')
        line = axes.semilogy(abscissa, data, **kwargs)
        axes.tick_params(axis='both', which='major', labelsize=15)
        axes.set_xlabel(_xlabel, fontsize=15)
        axes.set_ylabel('cost functions', fontsize=15)
        # fig.legend(loc='upper right', frameon=False, prop={'size':24, 'weight':'bold'})
        plt.legend()
        return line
    # poincare plot
    def poincare(self, prange='full', **kwargs):
        """
        Adaptation of CoilPy[old]'s poincare plot.
        prange -- str, ['full'(default), 'upper', 'lower']. Range of plotted points.
        all kwargs are passed to the scatter function
        """
        r = self.ppr
        z = self.ppz
        if kwargs.get('s') is None:
            kwargs.update({'s': .1})  # set size sensibly
        if kwargs.get('c') is None:
            kwargs.update({'c': 'k'})  # set color black
        if kwargs.get('marker') is None:
            kwargs.update({'marker': '.'})

        if plt.get_fignums():
            fig = plt.gcf()
            ax = plt.gca()
        else:
            fig, ax = plt.subplots()
        if prange == 'full':
            ax.scatter(r, z, **kwargs)
        elif prange == 'upper':
            ax.scatter(r[z>=0.0],  z[z>=0.0], **kwargs)
        elif prange == 'lower':
            ax.scatter(r[z<=0.0],  z[z<=0.0], **kwargs)
        return

    # Bnorm plot
    def Bnorm(self, plottype='2D', source='all', axes=None, flip=False, **kwargs):
        """ Plot Bn distribution.

        Keyword arguments:
            plottype -- string, '2D' (default) or '3D', determine the plottype
            source -- string, 'coil', 'plasma', 'sum' or 'all' (default), data source
            axes  -- matplotlib.pyplot or mayavi.mlab axis, axis to be
                     plotted on  (default None)
            flip -- logical, determine how to calculate Bn from coil,
                    True: coil_Bn = Bn - plas_Bn; False: coil_Bn = Bn + plas_Bn

        Returns:
           obj -- matplotlib.pyplot or mayavi.mlab plotting object
        """
        obj = []
        # check coil_Bn
        if flip:
            coil_Bn = self.Bn - self.plas_Bn
        else :
            coil_Bn = self.Bn + self.plas_Bn
        # 2D plots
        if plottype.lower() == '2d':
            # prepare coordinates
            ntheta  = self.Nteta
            nzeta = self.Nzeta
            nt, nz = self.Bn.shape
            if self.IsSymmetric >0 :
                zeta_end = 2*np.pi/self.Nfp
            else :
                zeta_end = 2*np.pi
            def theta(x, pos):
                return '{:3.2f}'.format(np.pi/ntheta + x*(2*np.pi)/ntheta)
            def zeta(x, pos):
                return '{:3.2f}'.format(0.5*zeta_end/nzeta + x*zeta_end/nzeta)
            theta_format = FuncFormatter(theta)
            zeta_format = FuncFormatter(zeta)
            # planar plotting
            if source.lower() == 'all':
                fig, axes = get_figure(axes, ncols=3, sharex=True)
                axes = np.atleast_1d(axes)
                #kwargs['aspect'] = kwargs.get('aspect', nt/(2.0*nz*self.Nfp))
                kwargs['aspect'] = kwargs.get('aspect', 'auto')
            else :
                fig, axes = get_figure(axes)
                axes = np.atleast_1d(axes)
                kwargs['aspect'] = kwargs.get('aspect', float(nt/nz))
            # prepare axes
            for ax in axes:
                ax.xaxis.set_major_formatter(zeta_format)
                ax.yaxis.set_major_formatter(theta_format)
            # set default plotting parameters
            kwargs['cmap'] = kwargs.get('cmap','RdBu_r') # colormap
            kwargs['origin'] = kwargs.get('origin', 'lower') # number of contours
            kwargs['extent'] = kwargs.get('extent', [0, nt, 0, nz]) # extent
            plt.subplots_adjust(hspace=0.05)
            # imshow
            if source.lower() == 'coil':
                obj.append(axes[0].imshow(np.transpose(coil_Bn), **kwargs))
                axes[0].set_title('Bn from coils', fontsize=15)
                axes[0].set_ylabel(r'$\theta$', fontsize=14)
                axes[0].set_xlabel(r'$\phi$', fontsize=14)
            elif source.lower() == 'plasma':
                obj.append(axes[0].imshow(np.transpose(self.plas_Bn), **kwargs))
                axes[0].set_title('Bn from plasma', fontsize=15)
                axes[0].set_ylabel(r'$\theta$', fontsize=14)
                axes[0].set_xlabel(r'$\phi$', fontsize=14)
            elif source.lower() == 'sum':
                obj.append(axes[0].imshow(np.transpose(self.Bn), **kwargs))
                axes[0].set_title('Residual Bn', fontsize=15)
                axes[0].set_ylabel(r'$\theta$', fontsize=14)
                axes[0].set_xlabel(r'$\phi$', fontsize=14)
            elif source.lower() == 'all':
                vmin = np.min([self.plas_Bn, self.Bn, coil_Bn])
                vmax = np.max([self.plas_Bn, self.Bn, coil_Bn])
                obj.append(axes[0].imshow(np.transpose(self.plas_Bn),
                                          vmin=vmin, vmax=vmax, **kwargs))
                #plt.setp(axes[0].get_xticklabels(), visible=False)
                axes[1].set_title('Bn from plasma (left), coil (mid), overall(right).', fontsize=14)
                obj.append(axes[1].imshow(np.transpose(coil_Bn),
                                          vmin=vmin, vmax=vmax, **kwargs))
                #axes[1].set_title('Bn from coils', fontsize=15)
                plt.setp(axes[1].get_yticklabels(), visible=False)
                obj.append(axes[2].imshow(np.transpose(self.Bn),
                                          vmin=vmin, vmax=vmax, **kwargs))
                plt.setp(axes[2].get_yticklabels(), visible=False)
                #axes[2].set_title('Residual Bn', fontsize=15)
                axes[0].set_ylabel(r'$\theta$', fontsize=14)
                #axes[1].set_ylabel(r'$\theta$', fontsize=14)
                #axes[2].set_ylabel(r'$\theta$', fontsize=14)
                axes[0].set_xlabel(r'$\phi$', fontsize=14)
                axes[1].set_xlabel(r'$\phi$', fontsize=14)
                axes[2].set_xlabel(r'$\phi$', fontsize=14)
            else:
                raise ValueError('unsupported option for source')
            fig.subplots_adjust(right=0.85)
            cbar_ax = fig.add_axes([0.86, 0.10, 0.04, 0.8])
            fig.colorbar(obj[0], cax=cbar_ax)
        return obj
    # Bmod plot
    def Bmod(self):
        return
    # write vtk
    def toVTK(self, name=None, full=False, **kwargs):
        """Save surface and magnetic field data into VTK file
        Arguments:
          name -- string, VTK file name. default: None, if None, using self.filename
          full -- boolean, if save the entire torus
          **kwargs -- external data will be saved in VTK file

        Return:
          VTK file name
        """
        from pyevtk.hl import gridToVTK
        # automatically get a file name, focus_*.h5 -> vtk_*
        if name is None:
            name = self.filename[:-3].replace('focus_', 'vtk_')
        xx = np.atleast_3d(map_matrix(self.xsurf, first=False))
        yy = np.atleast_3d(map_matrix(self.ysurf, first=False))
        zz = np.atleast_3d(map_matrix(self.zsurf, first=False))
        Bn = np.atleast_3d(map_matrix(self.Bn, first=False))
        plas_Bn =  np.atleast_3d(map_matrix(self.plas_Bn, first=False))
        B = (np.atleast_3d(map_matrix(self.Bx, first=False)),
             np.atleast_3d(map_matrix(self.By, first=False)),
             np.atleast_3d(map_matrix(self.Bz, first=False)))
        if full :
            pass
        data = {"Bn":Bn, "plas_Bn":plas_Bn, "B":B}
        data.update(kwargs)
        return gridToVTK(name, xx, yy, zz, pointData=data)
