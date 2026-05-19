from .spec import SPECout

# helper functions for plotting field-lines and interface boundary
class SPECoutplus(SPECout):
    def __init__(self, *args, **kwargs):
        super().__init__(*args,**kwargs)
        if hasattr(self, "poincare"):
            self.poincare.npoinc=self.input.numerics.Ndiscrete*4*self.input.physics.Ntor
            
    def plot_poincare(self,phi_indices,ncols=3,axes=None,**kwargs):
        import numpy as np
        import matplotlib.pyplot as plt
        
        R=self.poincare.R
        Z=self.poincare.Z
        
        nlines = R.shape[0]
        nplots = len(phi_indices)

        # colormap: one color per field line
        cmap = plt.cm.hsv
        colors = cmap(np.linspace(0, 1, nlines))
    
        # layout
        nrows = int(np.ceil(nplots / ncols))
    
        if axes is None:
            _, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 4*nrows))
            if ncols>1:
                axes = axes.flatten()
            else:
                axes=[axes]
        
        for k, i in enumerate(phi_indices):
            ax = axes[k]
        
            ## plot field-lines
            for j in range(nlines):
                mask = (R[j, :, i] != 0) & (Z[j, :, i] != 0)
                Ri = R[j,:,i][mask]
                Zi = Z[j,:,i][mask]
                ax.plot(Ri,Zi,'.',color=colors[j],markersize=0.5)
    
            ## styling
            ax.set_title(f"φ={i/self.grid.Nz:.2f}2π/nfp")
            ax.set_xlabel("R")
            ax.set_ylabel("Z")
            ax.set_aspect('equal')
        return axes

    def plot_grid(self,phi_indices,ncols=3,axes=None,**kwargs):
        import numpy as np
        import matplotlib.pyplot as plt

        nphi=len(phi_indices)
        nplots = len(phi_indices)
        nt=56
        phi_vals=phi_indices/self.grid.Nz*2*np.pi/self.input.physics.Nfp
        tarr=np.linspace(0,2*np.pi,nt)

        # layout
        nrows = int(np.ceil(nplots / ncols))
        if axes is None:
            _, axes = plt.subplots(nrows,ncols,figsize=(4*ncols,4*nrows))
            if ncols>1:
                axes = axes.flatten()
            else:
                axes=[axes]

        for lv in range(self.input.physics.Nvol):

            if self.input.physics.Nvol==1:
                ns=self.input.diagnostics.nPtrj+1
            else:
                ns=self.input.diagnostics.nPtrj[lv]+1
            sarr=np.linspace(-1,1,ns)
            for k,i in enumerate(phi_indices):
                ax = axes[k]
                Rb,Zb=self.get_grid(lvol=lv,sarr=sarr,tarr=tarr,zarr=phi_vals)

                if "c" not in kwargs:
                    kwargs.update({"c":"grey"})
                
                ax.plot(Rb[:,:,k],Zb[:,:,k],lw=0.7,**kwargs)
                ax.plot(Rb[:,:,k].T,Zb[:,:,k].T,lw=0.7,**kwargs)

                # cosmestics
                ax.set_title(f"φ={i/self.grid.Nz:.2f}2π/nfp")
                ax.set_xlabel("R")
                ax.set_ylabel("Z")
                ax.set_aspect('equal')
                ax.set_xlim((self.grid.Rmin,self.grid.Rmax))
                ax.set_ylim((self.grid.Zmin,self.grid.Zmax))
                ax.grid("on")
        return axes
    

    def plot_boundary_from_grid(self,phi_indices,ncols=3,axes=None,**kwargs):
        import numpy as np
        import matplotlib.pyplot as plt

        nphi=len(phi_indices)
        nplots = len(phi_indices)
        #memory efficient boundary array
        Rb=np.empty(self.grid.Nt+1,dtype=self.grid.Rij[0].dtype)
        Zb=np.empty(self.grid.Nt+1,dtype=self.grid.Rij[0].dtype)

        # layout
        nrows = int(np.ceil(nplots / ncols))
        if axes is None:
            _, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 4*nrows))
            if ncols>1:
                axes = axes.flatten()
            else:
                axes=[axes]

        for lv in range(self.input.physics.Nvol):
            for k, i in enumerate(phi_indices):
                ax = axes[k]
        
                ## plot outer boundary
                Rb[:-1] = self.grid.Rij[lv][i*self.grid.Nt:(i+1)*self.grid.Nt,-1]
                Zb[:-1] = self.grid.Zij[lv][i*self.grid.Nt:(i+1)*self.grid.Nt,-1]
                # periodize
                Rb[-1] = self.grid.Rij[lv][i*self.grid.Nt,-1]
                Zb[-1] = self.grid.Zij[lv][i*self.grid.Nt,-1]

                if "c" not in kwargs:
                    kwargs.update({"c":"grey"})
                ax.plot(Rb,Zb,**kwargs)

                ## plot inner boundary (axis)
                if lv == 0:
                    Ra= self.grid.Rij[lv][i*self.grid.Nt,0]
                    Za = self.grid.Zij[lv][i*self.grid.Nt,0]
                    ax.scatter(Ra,Za,marker='x',c="black",s=50,zorder=2)

                ax.set_title(f"φ={i/self.grid.Nz:.2f}2π/nfp")
                ax.set_xlabel("R")
                ax.set_ylabel("Z")
                ax.set_aspect('equal')
                ax.set_xlim((self.grid.Rmin,self.grid.Rmax))
                ax.set_ylim((self.grid.Zmin,self.grid.Zmax))
                ax.grid("on")
        return axes

    def extract_boundary(self):
        import numpy as np
        from ..input.boundary_diagnostics import ToroidalSurface
        m_modes = np.arange(-self.input.physics.Mpol,self.input.physics.Mpol+1)
        n_modes = np.arange(-self.input.physics.Ntor,self.input.physics.Ntor+1)*self.input.physics.Nfp
        
        # assuming stellarator symmetry
        Rmn = self.input.physics.Rbc  # major + minor radius
        Zmn = self.input.physics.Zbs

        return ToroidalSurface(m_modes, n_modes, Rmn, Zmn)
    
    def plot(self,zetastep=1,title=None,ncols=None,outfile=None,axes=None,show_grid=False,**kwargs):
        import numpy as np
        import matplotlib.pyplot as plt
        nphi = self.poincare.R.shape[2]
        phi_indices = list(range(0, nphi, zetastep))
        phi_vals=phi_indices/self.grid.Nz*2*np.pi/self.input.physics.Nfp

        if ncols is None:
            ncols=len(phi_indices)

        
        if show_grid:
            if axes is None:
                axes=self.plot_grid(phi_indices,ncols=ncols,**kwargs)
            else:
                axes=self.plot_grid(phi_indices,axes=axes,**kwargs)
        
        if axes is None:
            axes=self.plot_poincare(phi_indices,ncols=ncols,**kwargs)
        else:
            axes=self.plot_poincare(phi_indices,axes=axes,**kwargs)
        
        axes=self.plot_boundary_from_grid(phi_indices,axes=axes,**kwargs)
        if title is not None:
            plt.suptitle(title,fontsize=14)
            # Adjust layout so suptitle doesn't overlap subplots
            plt.tight_layout(rect=[0, 0, 1, 0.98])  # leave space on top for suptitle
        
        # Save figure as PNG
        if outfile is not None:
            plt.savefig(outfile, dpi=300, bbox_inches='tight')  # dpi for resolution
            
        return axes,phi_vals


    def minimum_jacobian(self,ns=52,nt=64,nz=96):
        import numpy as np
        nfp = self.input.physics.Nfp
        nvol = self.input.physics.Nvol

        tarr    = np.linspace(0,2*np.pi,nt)
        zarr    = np.linspace(0,2*np.pi/nfp,nz)

        minjac = 9999999.0
        smin = np.nan
        tmin = np.nan
        zmin = np.nan
        lvolmin = np.nan
        
        for ivol in range(0,nvol):
            if ivol==0: sarr=np.linspace(-0.999,1,ns)
            else: sarr=np.linspace(-1,1,ns)
            
            _,_, jacobian,_ = self.get_grid_and_jacobian_and_metric(lvol=ivol,sarr=sarr,tarr=tarr,zarr=zarr)

            si,ti,zi = np.unravel_index(np.argmin(jacobian), jacobian.shape)
            J_min = jacobian[si,ti,zi]

            print(jacobian.shape)
            
            if minjac>J_min:
                minjac = J_min
                smin = sarr[si]
                tmin = tarr[ti]
                zmin = zarr[zi]
                lvolmin = ivol 
        return minjac,smin,tmin,zmin, lvolmin
            
    
    def __str__(self):
        import numpy as np
        jmin,smin,tmin,zmin,lvolmin = self.minimum_jacobian()
        torfluxtot = 0.
        polfluxtot = 0.
        for lv in range(self.input.physics.Nvol):
            torfluxtot+=self.get_torflux(lv)
            polfluxtot+=self.get_polflux(lv)
        
        s=f"""\
        SPEC output: {self.filename}
        ===                     # 
        Inputs:
        * number of volumes: {self.input.physics.Nvol}
        * poloidal mode numbers: {self.input.physics.Mpol}
        * toroidal mode numbers: {self.input.physics.Ntor}
        * radial polynomials: {self.input.physics.Lrad}
        * field periodicity: {self.input.physics.Nfp}
        * toroidal flux: {self.input.physics.phiedge}
        * helicity: {self.input.physics.helicity}
        * mu Lagrange multiplier: {self.input.physics.mu}
        ---
        Outputs:
        * helicity: {self.output.helicity}
        * mu Lagrange multiplier: {self.output.mu}
        * total plasma volume : {self.output.volume*self.input.physics.Nfp}
        * total toroidal flux: {torfluxtot}
        * total poloidal flux: {polfluxtot}
        * poloidal flux ratio: {self.output.pflux}
        * toroidal flux ratio: {self.output.tflux}
        * minimum jacobian: {np.round(jmin,3)} in lvol={lvolmin} at (s={np.round(smin,3)},th={np.round(tmin,3)},z={np.round(zmin,3)})
        """
        if hasattr(self,"poincare"):
            s += f"""* poincare trajectories: {self.poincare.npoinc}\n"""
        s+=f"        ---\n"
        return s


