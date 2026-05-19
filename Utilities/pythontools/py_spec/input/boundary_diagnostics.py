import numpy as np

class ToroidalSurface:
    @classmethod
    def from_grid(cls, m_vals, n_vals, Rmn, Zmn):
        """
        Construct directly from dense (m,n) coefficient grids.
        
        Parameters
        ----------
        Rmn, Zmn : array_like, shape (Nm, Nn)
            Fourier coefficient grids
        m_vals : array_like
            m mode numbers
        n_vals : array_like
            n mode numbers
        """
        return cls(
            m_vals=np.asarray(m_vals),
            n_vals=np.asarray(n_vals),
            Rmn=np.asarray(Rmn),
            Zmn=np.asarray(Zmn),
        )
    
    @classmethod
    def from_sparse(cls, Rcoeffs, Zcoeffs, im, in_, m_vals=None, n_vals=None):
        """
        Construct from sparse mode-list representation.
        
        Parameters
        ----------
        Rcoeffs, Zcoeffs : (N,)
        Fourier coefficients
        im, in_ : (N,)
        mode numbers
        m_vals, n_vals : optional
        full mode grids. If omitted, inferred from sparse modes.
        """
        
        Rcoeffs = np.asarray(Rcoeffs)
        Zcoeffs = np.asarray(Zcoeffs)
        im = np.asarray(im)
        in_ = np.asarray(in_)
        
        if m_vals is None:
            m_vals = np.arange(im.min(), im.max() + 1)

        if n_vals is None:
            n_vals = np.arange(in_.min(), in_.max() + 1)

        m_vals = np.asarray(m_vals)
        n_vals = np.asarray(n_vals)

        def modes_to_grid(coeffs):
            grid = np.zeros((len(m_vals), len(n_vals)))

            m_index = {m: i for i, m in enumerate(m_vals)}
            n_index = {n: j for j, n in enumerate(n_vals)}

            for c, m, n in zip(coeffs, im, in_):
                grid[m_index[m], n_index[n]] += c

            return grid

        Rmn = modes_to_grid(Rcoeffs)
        Zmn = modes_to_grid(Zcoeffs)
        
        return cls(
            m_vals=m_vals,
            n_vals=n_vals,
            Rmn=Rmn,
            Zmn=Zmn
        )

    def __init__(self, m_vals, n_vals, Rmn, Zmn):
        self.Rmn = Rmn
        self.Zmn = Zmn
        self.m = m_vals
        self.n = n_vals
         # sanity checks
        assert self.Rmn.shape == (len(self.m), len(self.n)), \
            f"Rmn shape {self.Rmn.shape} incompatible with m,n"
        assert self.Zmn.shape == (len(self.m), len(self.n)), \
            f"Zmn shape {self.Zmn.shape} incompatible with m,n"
        self.Ntor=(len(self.n)-1)//2  # this is not general
        self.Mpol=(len(self.m)-1)//2  # this is not general
        
    def evaluate(self, u, v):
        """
        u, v can be scalars or arrays (same shape)
        """
        u = np.asarray(u)
        v = np.asarray(v)

        # broadcast to (m, n, points)
        phase = (
            self.m[:, None, None] * u[None, None, :] -
            self.n[None, :, None] * v[None, None, :]
        )

        cos_phase = np.cos(phase)
        sin_phase = np.sin(phase)

        R = np.sum(self.Rmn[:, :, None] * cos_phase, axis=(0, 1))
        Z = np.sum(self.Zmn[:, :, None] * sin_phase, axis=(0, 1))

        return R, Z

    def derivatives(self, u, v):
        u = np.asarray(u)
        v = np.asarray(v)

        phase = (
            self.m[:, None, None] * u[None, None, :] -
            self.n[None, :, None] * v[None, None, :]
        )

        cos_phase = np.cos(phase)
        sin_phase = np.sin(phase)

        m = self.m[:, None, None]
        n = self.n[None, :, None]

        dR_du = np.sum(-self.Rmn[:, :, None] * m * sin_phase, axis=(0,1))
        dR_dv = np.sum( self.Rmn[:, :, None] * n * sin_phase, axis=(0,1))

        dZ_du = np.sum( self.Zmn[:, :, None] * m * cos_phase, axis=(0,1))
        dZ_dv = np.sum(-self.Zmn[:, :, None] * n * cos_phase, axis=(0,1))

        return dR_du, dR_dv, dZ_du, dZ_dv

    def tangent_vectors(self, u, v):
        R, Z = self.evaluate(u, v)
        dR_du, dR_dv, dZ_du, dZ_dv = self.derivatives(u, v)

        cosv = np.cos(v)
        sinv = np.sin(v)

        ru = np.vstack([
            dR_du * cosv,
            dR_du * sinv,
            dZ_du
        ])

        rv = np.vstack([
            dR_dv * cosv - R * sinv,
            dR_dv * sinv + R * cosv,
            dZ_dv
        ])

        return ru, rv

    def jacobian(self, u, v):
        ru, rv = self.tangent_vectors(u, v)
        cross = np.cross(ru.T, rv.T)
        return np.linalg.norm(cross, axis=1)


    # ---------- core geometry ----------
    def _compute_geometry(self, U, V):
        m = self.m[:, None, None, None]
        n = self.n[None, :, None, None]

        phase = m * U[None,None,:,:] - n * V[None,None,:,:]

        R = np.sum(self.Rmn[:, :, None, None] * np.cos(phase), axis=(0,1))
        Z = np.sum(self.Zmn[:, :, None, None] * np.sin(phase), axis=(0,1))

        X = R * np.cos(V)
        Y = R * np.sin(V)

        # derivatives
        dR_du = np.sum(-self.Rmn[:, :, None, None] * m * np.sin(phase), axis=(0,1))
        dR_dv = np.sum( self.Rmn[:, :, None, None] * n * np.sin(phase), axis=(0,1))

        dZ_du = np.sum( self.Zmn[:, :, None, None] * m * np.cos(phase), axis=(0,1))
        dZ_dv = np.sum(-self.Zmn[:, :, None, None] * n * np.cos(phase), axis=(0,1))

        cosv = np.cos(V)
        sinv = np.sin(V)

        ru = np.stack([
            dR_du * cosv,
            dR_du * sinv,
            dZ_du
        ], axis=-1)

        rv = np.stack([
            dR_dv * cosv - R * sinv,
            dR_dv * sinv + R * cosv,
            dZ_dv
        ], axis=-1)

        normals = np.cross(rv, ru)
        J = np.linalg.norm(normals, axis=2)

        # normalize normals
        normals_unit = normals / (J[..., None] + 1e-12)

         # enforce outward orientation automatically
        centroid = np.stack([X, Y, Z], axis=-1)
        if np.mean(np.sum(normals * centroid, axis=2)) < 0:
            print("flipping inward pointing normals")
            normals_unit = -normals_unit
        else:
            print("normals are outward pointing")

        return X, Y, Z, J, normals_unit

    def is_curve(self, tol=1e-14):
        """
        True if all m>0 modes vanish.
        """
        mpos = self.m != 0
        return np.all(np.abs(self.Rmn[mpos, :]) < tol) and np.all(np.abs(self.Zmn[mpos, :]) < tol)


    def plot_cross_section(self, v0, npts=300,ax=None,**kwargs):
        
        import matplotlib.pyplot as plt

        is_curve=self.is_curve()
        
        if is_curve:
            u = np.linspace(0,2*np.pi,1)
        else:
            u = np.linspace(0, 2*np.pi, npts)
            
        v = np.full_like(u, v0)
    
        R, Z = self.evaluate(u, v)

        if ax is None:
            _,ax = plt.subplots()

        if is_curve:
            ax.scatter(R,Z,**kwargs)
        else:
            ax.plot(R, Z,**kwargs)
        ax.set_xlabel("R")
        ax.set_ylabel("Z")
        ax.set_title(f"φ={v0:.2f}")
        ax.set_aspect('equal')
        return ax

     # ---------- PyVista export ----------
    def to_pyvista(self, Nu=100, Nv=100):
        if self.is_curve():
            return self._curve_to_pyvista(Nv=Nv)
        else:
            return self._surface_to_pyvista(Nu=Nu, Nv=Nv)

    def _curve_to_pyvista(self, Nv=400):
        import pyvista as pv
        import numpy as np

        v = np.linspace(0, 2*np.pi, Nv)
        
        n = self.n[None, :]
        V = v[:, None]
        
        phase = -n * V
        
        # only m=0 row
        i0 = np.where(self.m == 0)[0][0]
        
        R = np.sum(self.Rmn[i0, :] * np.cos(phase), axis=1)
        Z = np.sum(self.Zmn[i0, :] * np.sin(phase), axis=1)
        
        X = R * np.cos(v)
        Y = R * np.sin(v)
        
        points = np.column_stack([X, Y, Z])

        # closed polyline
        curve = pv.lines_from_points(points, close=True)
        
        return curve

    
    def _surface_to_pyvista(self,Nu=100,Nv=100):
        
        import pyvista as pv

        u = np.linspace(0, 2*np.pi, Nu)
        v = np.linspace(0, 2*np.pi, Nv)
        U, V = np.meshgrid(u, v, indexing='ij')

        X, Y, Z, J, normals = self._compute_geometry(U, V)

        grid = pv.StructuredGrid(X, Y, Z)

        # Fortran consistent ordering
        grid.point_data["Jacobian"] = J.flatten(order="F")
        grid.point_data["Normals"] = normals.reshape(-1, 3,order="F")

        grid.set_active_scalars("Jacobian")
        grid.set_active_vectors("Normals")

        return grid

    
    def export_vtk(self, filename, Nu=100, Nv=100):
        """
        Export surface using class pipeline (Fortran-consistent).
        """
        grid = self.to_pyvista(Nu=Nu, Nv=Nv)
    
        if not filename.endswith(".vts"):
            filename += ".vts"
    
        grid.save(filename)
    
        print(f"Saved VTK surface to: {filename}")
        
    # ---------- convenience plotting ----------
    def plot(self, Nu=100, Nv=100, show_normals=False):
        import pyvista as pv
        
        grid = self.to_pyvista(Nu, Nv)

        plotter = pv.Plotter()
        plotter.add_mesh(grid, scalars="Jacobian", cmap="viridis")

        if show_normals:
            glyphs = grid.glyph(orient="Normals", scale=False, factor=0.3)
            plotter.add_mesh(glyphs, color="black")

        plotter.add_axes()
        plotter.show()

    # ---------- highlight min Jacobian ----------
    def plot_min_jacobian(self, Nu=100, Nv=100, tol=1.05):
        import pyvista as pv
        
        grid = self.to_pyvista(Nu, Nv)

        J = grid["Jacobian"]
        J_min = J.min()

        mask = J < tol * J_min
        points = grid.points[mask]

        plotter = pv.Plotter()
        plotter.add_mesh(grid, scalars="Jacobian", cmap="viridis")

        plotter.add_mesh(
            pv.PolyData(points),
            color="red",
            point_size=10,
            render_points_as_spheres=True
        )

        plotter.add_axes()
        plotter.show()

        print("Min Jacobian:", J_min)
        print("Points highlighted:", len(points))
        
    def plot_fouriermap(self,m_max,n_max,threshold=1.e-5):
        import matplotlib.pyplot as plt
        
        # --- colormaps ---
        cmap_div = plt.cm.seismic.copy()
        cmap_div.set_bad(color='gray')

        cmap_amp = plt.cm.viridis.copy()
        cmap_amp.set_bad(color='gray')

        # assuming stellarator symmetry 
        R=self.Rmn
        Z=self.Zmn
         # --- combined amplitude ---
        A = np.sqrt(R**2 + Z**2)

        fig, axes = plt.subplots(nrows=1,ncols=3,figsize=(18, 5),sharey=True)

        vmax=max(np.max(np.abs(R)),np.max(np.abs(Z)),np.max(np.abs(A)))
        # aliases
        ntor=self.Ntor
        mpol=self.Mpol
        
        for ax,data,title in zip(
            axes[:2],
            [R,Z],
            [r"$R_{mn}$",r"$Z_{mn}$"]
        ):
            d_masked = np.ma.masked_where(np.abs(data)<threshold,data)
            im=ax.imshow(d_masked,
                         cmap=cmap_div,
                         origin="lower",
                         aspect="auto",
                         vmin=-vmax,vmax=vmax,
                         extent=[-ntor-0.5,ntor+0.5,-mpol-0.5,mpol+0.5]
                         )
            # --- restrict n range ---
            ax.set_xlim(-n_max, n_max)
            ax.set_ylim(-1, m_max)
             # --- labels ---
            ax.set_xlabel(r"$n$")
            ax.set_title(title)
        
        A_masked=np.ma.masked_where(np.abs(A)<threshold,A)

        im_amp=axes[2].imshow(A_masked,
                              cmap=cmap_amp,
                              origin="lower",
                              aspect="auto",
                              extent=[-ntor-0.5,ntor+0.5,-mpol-0.5,mpol+0.5]
                              )
        
        axes[2].set_xlim(-n_max, n_max)
        axes[2].set_ylim(-1, m_max)
        axes[2].set_xlabel(r"$n$")
        axes[2].set_title(r"$\sqrt{R_{mn}^2+Z_{mn}^2}$")

        # shared y-axis
        axes[0].set_ylabel(r"$m$")
        
        # shared colorbar
        cbar1 = fig.colorbar(im, ax=axes[:2],location="left")
        cbar1.set_label(r"$R_{mn}, Z_{mn}$")

        cbar2 = fig.colorbar(im_amp, ax=axes[2])
        cbar2.set_label(r"Amplitude")

        fig.suptitle("Fourier modes of boundary")

        return fig, axes

