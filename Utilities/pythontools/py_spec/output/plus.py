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
    
        for k, i in enumerate(phi_indices):
            ax = axes[k]
        
            ## plot outer boundary
            Rb[:-1] = self.grid.Rij[0][i*self.grid.Nt:(i+1)*self.grid.Nt,-1]
            Zb[:-1] = self.grid.Zij[0][i*self.grid.Nt:(i+1)*self.grid.Nt,-1]
            # periodize
            Rb[-1] = self.grid.Rij[0][i*self.grid.Nt,-1]
            Zb[-1] = self.grid.Zij[0][i*self.grid.Nt,-1]

            if "c" not in kwargs:
                kwargs.update({"c":"grey"})
            ax.plot(Rb,Zb,**kwargs)
    
            ## plot inner boundary (axis)
            Ra= self.grid.Rij[0][i*self.grid.Nt,0]
            Za = self.grid.Zij[0][i*self.grid.Nt,0]
            ax.scatter(Ra,Za,marker='x',c="black",s=50,zorder=2)

            ax.set_title(f"φ={i/self.grid.Nz:.2f}2π/nfp")
            ax.set_xlabel("R")
            ax.set_ylabel("Z")
            ax.set_aspect('equal')
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
    
    def plot(self,zetastep=1,title=None,ncols=None,outfile=None,axes=None,**kwargs):
        import numpy as np
        import matplotlib.pyplot as plt
        nphi = self.poincare.R.shape[2]
        phi_indices = list(range(0, nphi, zetastep))

        if ncols is None:
            ncols=len(phi_indices)

        if axes is None:
            axes=self.plot_boundary_from_grid(phi_indices,ncols=ncols,**kwargs)
        else:
            axes=self.plot_boundary_from_grid(phi_indices,axes=axes,**kwargs)
            
        self.plot_poincare(phi_indices,axes=axes,**kwargs)
        if title is not None:
            plt.suptitle(title,fontsize=14)
            # Adjust layout so suptitle doesn't overlap subplots
            plt.tight_layout(rect=[0, 0, 1, 0.98])  # leave space on top for suptitle
        
        # Save figure as PNG
        if outfile is not None:
            plt.savefig(outfile, dpi=300, bbox_inches='tight')  # dpi for resolution
            
        return axes
    
    def __str__(self):
        s=f"""\
        SPEC output: {self.filename}
        ===
        Inputs:
        * number of volumes: {self.input.physics.Nvol}
        * poloidal mode numbers: {self.input.physics.Mpol}
        * toroidal mode numbers: {self.input.physics.Ntor}
        * radial Chebyshev polynomials: {self.input.physics.Lrad}
        * field periodicity: {self.input.physics.Nfp}
        * toroidal flux: {self.input.physics.phiedge}
        * helicity: {self.input.physics.helicity}
        * mu Lagrange multiplier: {self.input.physics.mu}
        ---
        Outputs:
        * helicity: {self.output.helicity}
        * mu Lagrange multiplier: {self.output.mu}
        * total plasma volume : {self.output.volume*self.input.physics.Nfp}
        * toroidal flux: {self.get_torflux(0)}
        """
        if hasattr(self,"poincare"):
            s += f"""* poincare trajectories: {self.poincare.npoinc}\n"""
        s+=f"        ---\n"
        return s


