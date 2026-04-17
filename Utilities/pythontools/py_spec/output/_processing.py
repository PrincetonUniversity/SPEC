import numpy as np
from scipy import integrate
import typing

def get_RZ_derivatives(
    self,
    lvol=0,
    sarr=np.linspace(1, 1, 1),
    tarr=np.linspace(0, 0, 1),
    zarr=np.linspace(0, 0, 1),
    input1D=False
):

    Igeometry = self.input.physics.Igeometry 

    if lvol==0 and (sarr==-1).any() and Igeometry!=0:
        raise ValueError('Cannot evaluate coordinate derivative on magnetic axis !')

    sym = self.input.physics.Istellsym == 1

    Rac, Rbc = self.output.Rbc[lvol : lvol + 2]
    Ras, Rbs = self.output.Rbs[lvol : lvol + 2]
    Zas, Zbs = self.output.Zbs[lvol : lvol + 2]
    Zac, Zbc = self.output.Zbc[lvol : lvol + 2]

    mn = Rac.size  # s.output.mn
    im = self.output.im
    in_ = self.output.in_

    #sbar = (sarr + 1) / 2
    sbar = 0.5*(sarr+1.0)
    fac = []

    rpol = self.input.physics.rpol
    rtor = self.input.physics.rtor

    if Igeometry == 1:
        for j in range(mn):
            fac.append([sbar, 0.5 * np.ones(sarr.size), np.zeros(sarr.size)])
    elif Igeometry == 2:
        for j in range(mn):
            if lvol > 0 or im[j] == 0:
                fac.append([sbar, 0.5 * np.ones(sarr.size), np.zeros(sarr.size)])
            else:
                fac.append(
                    [
                        sbar ** (im[j] + 1.0),
                        (im[j] + 1.0) / 2.0 * sbar ** (im[j]),
                        (im[j] + 1.0) * (im[j]) / 4.0 * sbar ** (im[j] - 1),
                    ]
                )
    elif Igeometry == 3:
        for j in range(mn):
            if lvol == 0 and im[j] == 0:
                fac.append([sbar ** 2, sbar, 0.5 * np.ones(sarr.size)])
            elif lvol == 0 and im[j] > 0:
                fac.append(
                    [
                        sbar ** im[j],
                        (im[j] / 2.0) * sbar ** (im[j] - 1.0),
                        (im[j] * (im[j] - 1) / 4.0) * sbar ** (im[j] - 2.0),
                    ]
                )
            else:
                fac.append([sbar, 0.5 * np.ones(sarr.size), np.zeros(sarr.size)])

    # now fac has the dimension (number of modes, number of derivatives, number of s points)
    fac = np.array(fac)
    # transpose to (number of derivatives, number of modes, number of s points)
    fac = np.moveaxis(fac, 0, 1)

    nax = np.newaxis
    if not input1D:
        im = im[:, nax, nax, nax]
        in_ = in_[:, nax, nax, nax]
        ang_arg = +im * tarr[nax, nax, :, nax] - in_ * zarr[nax, nax, nax, :]
    else:
        im = im[:, nax]
        in_ = in_[:, nax]
        ang_arg = im * tarr[nax, :] - in_ * zarr[nax, :]

    cos = np.cos(ang_arg)
    sin = np.sin(ang_arg)

    if not input1D:
        fac = fac[:, :, :, nax, nax]
        Rac = Rac[:, nax, nax, nax]
        Rbc = Rbc[:, nax, nax, nax]
        Zas = Zas[:, nax, nax, nax]
        Zbs = Zbs[:, nax, nax, nax]
        if not sym:
            Ras = Ras[:, nax, nax, nax]
            Rbs = Rbs[:, nax, nax, nax]
            Zac = Zac[:, nax, nax, nax]
            Zbc = Zbc[:, nax, nax, nax]
    else:
        Rac = Rac[:, nax]
        Rbc = Rbc[:, nax]
        Zas = Zas[:, nax]
        Zbs = Zbs[:, nax]
        if not sym:
            Ras = Ras[:, nax]
            Rbs = Rbs[:, nax]
            Zac = Zac[:, nax]
            Zbc = Zbc[:, nax]

    dR1 = Rac + fac[0] * (Rbc - Rac)
    Rarr0 = np.sum(dR1 * cos, axis=0)

    Rarr1 = np.sum(fac[1] * (Rbc - Rac) * cos, axis=0)
    Rarr2 = np.sum(-im * dR1 * sin, axis=0)
    Rarr3 = np.sum(in_ * dR1 * sin, axis=0)

    Rarr = np.array([Rarr1, Rarr2, Rarr3])

    # We only need Z for Igeometry=3
    if Igeometry == 3:
        dZ1 = Zas + fac[0] * (Zbs - Zas)
        Zarr0 = np.sum(dZ1 * sin, axis=0)

        Zarr1 = np.sum(fac[1] * (Zbs - Zas) * sin, axis=0)
        Zarr2 = np.sum(im * dZ1 * cos, axis=0)
        Zarr3 = np.sum(-in_ * dZ1 * cos, axis=0)

        Zarr = np.array([Zarr1, Zarr2, Zarr3])
    else:
        Zarr0 = None

    return Rarr0, Rarr1, Rarr2, Rarr3, Zarr0, Zarr1, Zarr2, Zarr3

def get_grid_and_jacobian_and_metric(
    self,
    lvol=0,
    sarr=np.linspace(1, 1, 1),
    tarr=np.linspace(0, 0, 1),
    zarr=np.linspace(0, 0, 1),
    input1D=False,
    derivative=False,
):
    r"""!Compute the metric and Jacobian on a given grid

    @param lvol (int, optional): The SPEC volume of interest, starting from 0. Defaults to 0.
    @param sarr (1D numpy array, optional): The s grid. Defaults to np.linspace(1,1,1).
    @param tarr (1D numpy array, optional): The $\theta$ grid. Defaults to np.linspace(0,0,1).
    @param zarr (1D numpy array, optional): The $\zeta$ grid. Defaults to np.linspace(0,0,1).
    @param input1D If sarr, tarr and zarr should be treated as a grid in 3D or just 1D input. Defaults to False
    @param derivative (bool, optional): If the derivatives of jacobian and $g_{ij}$ is needed.

    @returns Rarr0, Zarr0, jacobian, g, [djacobian, dg]: $R, Z, J, g_{ij}$. If derivative==True also return the derivative of $J, g_{ij}$ (derivative is the first dimension).
    """
    sym = self.input.physics.Istellsym == 1

    Rac, Rbc = self.output.Rbc[lvol : lvol + 2]
    Ras, Rbs = self.output.Rbs[lvol : lvol + 2]
    Zas, Zbs = self.output.Zbs[lvol : lvol + 2]
    Zac, Zbc = self.output.Zbc[lvol : lvol + 2]

    mn = Rac.size  # s.output.mn
    im = self.output.im
    in_ = self.output.in_

    #sbar = (sarr + 1) / 2
    sbar = 0.5*(sarr+1.0)

    Igeometry = self.input.physics.Igeometry
    rpol = self.input.physics.rpol
    rtor = self.input.physics.rtor

    im_arr = im[:,None] # shape (mn,1)
    sbar_arr = sbar[None, :]  # shape (1, ns)

    sbar_full = np.broadcast_to(sbar_arr, (im.size, sbar.size)) # shape (mm,ns)
    fac0 =  sbar_full.copy()
    fac1 = 0.5 * np.ones_like(sbar_full)
    fac2 = np.zeros_like(sbar_full)

    if Igeometry == 2:

        mask = (lvol == 0) & (im_arr != 0)
        mask = np.broadcast_to(mask, fac0.shape)
        
        fac0[mask] = sbar_full[mask] ** (im_arr[mask] + 1)
        fac1[mask] = (im_arr[mask] + 1)/2 * sbar_full[mask] ** im_arr[mask]
        fac2[mask] = (im_arr[mask]*(im_arr[mask]-1)/4) * sbar_full[mask] ** (im_arr[mask]-1)

    elif Igeometry == 3:
        
        if lvol == 0:
            mask0 = (im_arr == 0)
            maskp = (im_arr > 0)

            mask0 = np.broadcast_to(mask0, fac0.shape)
            maskp = np.broadcast_to(maskp, fac0.shape)

            # precompute powers safely (with broadcasting)
            sbar_pow_im   = sbar_full ** im_arr
            sbar_pow_im1  = sbar_full ** (im_arr - 1)
            sbar_pow_im2  = sbar_full ** (im_arr - 2)

            fac0[mask0] = sbar_full[mask0]**2
            fac1[mask0] = sbar_full[mask0]
            fac2[mask0] = 0.5

            fac0[maskp] = sbar_pow_im[maskp]
            fac1[maskp] = (im_arr/2.0 * sbar_pow_im1)[maskp]
            fac2[maskp] = ((im_arr*(im_arr-1)/4.0) * sbar_pow_im2)[maskp]
        
    # stack to match your original structur
    fac = np.stack([fac0, fac1, fac2], axis=0)

    ## old code
    # now fac has the dimension (number of modes, number of derivatives, number of s points)
    # fac = np.array(fac)  ### this is for old code
    # transpose to (number of derivatives, number of modes, number of s points)
    # fac = np.moveaxis(fac, 0, 1)

    nax = np.newaxis
    if not input1D:
        im = im[:, nax, nax, nax]
        in_ = in_[:, nax, nax, nax]
        ang_arg = +im * tarr[nax, nax, :, nax] - in_ * zarr[nax, nax, nax, :]
    else:
        im = im[:, nax]
        in_ = in_[:, nax]
        ang_arg = im * tarr[nax, :] - in_ * zarr[nax, :]
    
    cos = np.cos(ang_arg)
    sin = np.sin(ang_arg)

    if not input1D:
        fac = fac[:, :, :, nax, nax]
        Rac = Rac[:, nax, nax, nax]
        Rbc = Rbc[:, nax, nax, nax]
        Zas = Zas[:, nax, nax, nax]
        Zbs = Zbs[:, nax, nax, nax]
        if not sym:
            Ras = Ras[:, nax, nax, nax]
            Rbs = Rbs[:, nax, nax, nax]
            Zac = Zac[:, nax, nax, nax]
            Zbc = Zbc[:, nax, nax, nax]
    else:
        Rac = Rac[:, nax]
        Rbc = Rbc[:, nax]
        Zas = Zas[:, nax]
        Zbs = Zbs[:, nax]
        if not sym:
            Ras = Ras[:, nax]
            Rbs = Rbs[:, nax]
            Zac = Zac[:, nax]
            Zbc = Zbc[:, nax]

    dR = (Rbc - Rac)
    fac1_dR = fac[1] * dR
    dR1 = Rac + fac[0] * dR
    Rarr0 = np.einsum("m...,m...->...", dR1, cos)
    Rarr1 = np.einsum("m...,m...->...", fac1_dR, cos)
    Rarr2 = -np.einsum("m...,m...,m...->...", im, dR1, sin)
    Rarr3 =  np.einsum("m...,m...,m...->...", in_, dR1, sin)

    Rarr = np.array([Rarr1, Rarr2, Rarr3])

    # We only need Z for Igeometry=3
    if Igeometry == 3:
        dZ=(Zbs-Zas)        
        dZ1 = Zas + fac[0] *dZ
        fac1_dZ=fac[1]*dZ
        Zarr0 = np.einsum("m...,m...->...", dZ1, sin)
        Zarr1 = np.einsum("m...,m...->...", fac1_dZ, sin)
        Zarr2 = np.einsum("m...,m...,m...->...", im, dZ1, cos)
        Zarr3 =-np.einsum("m...,m...,m...->...", in_, dZ1, cos)
        
        Zarr = np.array([Zarr1, Zarr2, Zarr3])
    else:
        Zarr0 = None

    # If the derivative of g and jacobian is needed
    if derivative:

        im2 = im * im
        in2 = in_ * in_
        imin = im * in_
        fac2_dR = fac[2] * dR
        dR1_cos = dR1 * cos

        Rarr11 = np.einsum("m...,m...->...", fac2_dR, cos)
        Rarr12 = -np.einsum("m...,m...,m...->...", im, fac1_dR, sin)
        Rarr13 =  np.einsum("m...,m...,m...->...", in_, fac1_dR, sin)
        
        Rarr22 = -np.einsum("m...,m...->...", im2, dR1_cos)
        Rarr23 =  np.einsum("m...,m...->...", imin, dR1_cos)
        Rarr33 = -np.einsum("m...,m...->...", in2, dR1_cos)

        dRarr = np.array(
            [
                [Rarr11, Rarr12, Rarr13],
                [Rarr12, Rarr22, Rarr23],
                [Rarr13, Rarr23, Rarr33],
            ]
        )

        if Igeometry == 3:
            
            fac2_dZ = fac[2] * dZ
            dZ1_sin = dZ1 * sin

            Zarr11 = np.einsum("m...,m...->...", fac2_dZ, sin)
            Zarr12 = np.einsum("m...,m...,m...->...", im, fac1_dZ, cos)
            Zarr13 =-np.einsum("m...,m...,m...->...", in_, fac1_dZ, cos)
            Zarr22 =-np.einsum("m...,m...->...", im2, dZ1_sin)
            Zarr23 = np.einsum("m...,m...->...", imin, dZ1_sin)
            Zarr33 =-np.einsum("m...,m...->...", in2, dZ1_sin)

            dZarr = np.array(
                [
                    [Zarr11, Zarr12, Zarr13],
                    [Zarr12, Zarr22, Zarr23],
                    [Zarr13, Zarr23, Zarr33],
                ]
            )

    if Igeometry == 1:
        jacobian = Rarr1 * rpol * rtor

        g = Rarr[:, nax, :] * Rarr[nax, :, :]
        # g22
        g[1, 1, :] += rpol ** 2
        # g33
        g[2, 2, :] += rtor ** 2

        if derivative:
            djacobian = dRarr[0, :] * rpol * rtor
            dg = (
                dRarr[:, :, nax, :] * Rarr[nax, nax, :, :]
                + dRarr[:, nax, :, :] * Rarr[nax, :, nax, :]
            )

    if Igeometry == 2:
        jacobian = Rarr1 * Rarr0

        g = Rarr[:, nax, :] * Rarr[nax, :, :]
        # g22
        g[1, 1, :] += Rarr0 ** 2
        # g33
        g[2, 2, :] += 1.0

        if derivative:
            djacobian = dRarr[0, :] * Rarr0[nax, :] + Rarr1[nax, :] * Rarr
            dg = (
                dRarr[:, :, nax, :] * Rarr[nax, nax, :, :]
                + dRarr[:, nax, :, :] * Rarr[nax, :, nax, :]
            )
            dg[:, 1, 1, :] += 2.0 * Rarr * Rarr0[nax, :]

    elif Igeometry == 3:
        jacobian = Rarr0 * (Rarr2 * Zarr1 - Rarr1 * Zarr2)  # from matlab

        g = Rarr[:, nax, :] * Rarr[nax, :, :] + Zarr[:, nax, :] * Zarr[nax, :, :]
        g[2, 2, :] += Rarr0 ** 2

        if derivative:
            djacobian = (
                Rarr * (Rarr2 * Zarr1 - Rarr1 * Zarr2)[nax, :]
                + Rarr0[nax, :]
                * (dRarr[1, :] * Zarr1[nax, :] - dRarr[0, :] * Zarr2[nax, :])
                + Rarr0[nax, :]
                * (Rarr2[nax, :] * dZarr[0, :] - Rarr1[nax, :] * dZarr[1, :])
            )
            dg = (
                dRarr[:, :, nax, :] * Rarr[nax, nax, :, :]
                + dRarr[:, nax, :, :] * Rarr[nax, :, nax, :]
                + dZarr[:, :, nax, :] * Zarr[nax, nax, :, :]
                + dZarr[:, nax, :, :] * Zarr[nax, :, nax, :]
            )
            dg[:, 2, 2, :] += 2 * Rarr0[nax, :] * Rarr

    # moving axis - move dofs of coordinates to the front
    g = np.moveaxis(g, (0,1), (-2,-1))
    if derivative:
        djacobian = np.moveaxis(djacobian, 0, -1)
        dg = np.moveaxis(dg, (0,1,2), (-3,-2,-1))

    if derivative:
        return Rarr0, Zarr0, jacobian, g, djacobian, dg
    else:
        return Rarr0, Zarr0, jacobian, g

def get_grid(
    self,
    lvol=0,
    sarr=np.linspace(1, 1, 1),
    tarr=np.linspace(0, 0, 1),
    zarr=np.linspace(0, 0, 1),
    input1D=False,
):

    Rarr0, Zarr0, _, _ = get_grid_and_jacobian_and_metric(
        self, lvol, sarr, tarr, zarr, input1D=input1D
    )
    return Rarr0, Zarr0

def get_jacobian(
    self,
    lvol=0,
    sarr=np.linspace(1, 1, 1),
    tarr=np.linspace(0, 0, 1),
    zarr=np.linspace(0, 0, 1),
    input1D=False,
):

    _, _, jacobian, _ = get_grid_and_jacobian_and_metric(
        self, lvol, sarr, tarr, zarr, input1D=input1D
    )
    return jacobian

def get_metric(
    self,
    lvol=0,
    sarr=np.linspace(1, 1, 1),
    tarr=np.linspace(0, 0, 1),
    zarr=np.linspace(0, 0, 1),
    input1D=False,
):

    _, _, _, g = get_grid_and_jacobian_and_metric(
        self, lvol, sarr, tarr, zarr, input1D=input1D
    )
    return g

def get_B(
    self,
    lvol=0,
    jacobian=None,
    sarr=np.linspace(0, 0, 1),
    tarr=np.linspace(0, 0, 1),
    zarr=np.linspace(0, 0, 1),
    input1D=False,
    derivative=False,
    djacobian=None,
):
    r"""!Compute the contravariant components of the magnetic field $(B^s, B^\theta, B^\zeta)$

    @param lvol (int, optional): The SPEC volume of interest, starting from 0. Defaults to 0.
    @param jacobian(numpy array, optional): if jacobian is already computed, provide it here
    @param sarr (1D numpy array, optional): The s grid. Defaults to np.linspace(1,1,1).
    @param tarr (1D numpy array, optional): The $\theta$ grid. Defaults to np.linspace(0,0,1).
    @param zarr (1D numpy array, optional): The $\zeta$ grid. Defaults to np.linspace(0,0,1).
    @param input1D If sarr, tarr and zarr should be treated as a grid in 3D or just 1D input. Defaults to False
    @param derivative (bool, optional): If the derivatives is needed.

    @returns Bcontrav, [dBcontrav]: $(B^s, B^\theta, B^\zeta)$. If derivative==True also return the derivative (derivative is the first dimension).
    """

    if not derivative:
        if jacobian is None:
            R, Z, jacobian, g = get_grid_and_jacobian_and_metric(
                self,
                lvol=lvol,
                sarr=sarr,
                tarr=tarr,
                zarr=zarr,
                input1D=input1D,
                derivative=derivative,
            )
    else:
        if jacobian is None or djacobian is None:
            R, Z, jacobian, g, djacobian, _ = get_grid_and_jacobian_and_metric(
                self,
                lvol=lvol,
                sarr=sarr,
                tarr=tarr,
                zarr=zarr,
                input1D=input1D,
                derivative=derivative,
            )

    nax = np.newaxis

    from pyoculus.problems import SPECBfield

    eq = SPECBfield(self, lvol=lvol + 1)
    if not derivative:
        B = eq.B_many(sarr, tarr, zarr, input1D=input1D)
    else:
        B, dBdX = eq.dBdX_many(sarr, tarr, zarr, input1D=input1D)

    Bcontrav = B / jacobian[...,nax]

    if derivative:
        dBcontrav = dBdX / jacobian[...,nax,nax] - djacobian[..., :, nax] * Bcontrav[..., nax, :] / jacobian[...,nax,nax]
        return Bcontrav, dBcontrav
    else:
        return Bcontrav

def get_modB(self, Bcontrav, g, derivative=False, dBcontrav=None, dg=None):
    """Input - Bcontrav has to come from get_B function"""
    modB = np.sqrt(np.einsum("...i,...ji,...j->...", Bcontrav, g, Bcontrav))
    if not derivative:
        return modB
    else:
        dmodB2 = 2 * np.einsum(
            "...ki,...ji,...j->...k", dBcontrav, g, Bcontrav
        ) + np.einsum("...i,...kji,...j->...k", Bcontrav, dg, Bcontrav)
        return modB, dmodB2

def get_B_covariant(self, Bcontrav=None, g=None, derivative=False):
    """Get covariant component of B"""
    Bco = np.einsum("...i,...ji->...j", Bcontrav, g)
    return Bco

def get_volume(self, ivol=0, ns=64, nt=64, nz=64):
    """Returns volume occupied by volume ivol"""
    
    # Create coordinate grid
    nfp = self.input.physics.Nfp
    tarr = np.linspace(0, 2*np.pi, nt, endpoint=True)
    zarr = np.linspace(0, 2*np.pi / nfp, nz, endpoint=True)

    if ivol==0: sarr=np.linspace(-0.999,1,ns)
    else: sarr=np.linspace(-1,    1, ns)

    # Get jacobian
    j = self.get_jacobian(lvol=ivol, sarr=sarr, tarr=tarr, zarr=zarr)

    # Integrate
    dt = tarr[1]-tarr[0]
    dz = zarr[1]-zarr[0]
    ds = sarr[1]-sarr[0]
    return nfp * integrate.simpson( y=integrate.simpson( y=integrate.simpson( y=j, x=zarr ), x=tarr ), x=sarr )

def get_area(self, ivol=0,ns=64,nt=64,phi0=0):
    """
    Calculates cross-sectional area of a given volume at fixed phi

    INPUT
    -----
    -data    : must be produced by calling read_spec(filename)
    -lvol    : volume number
    -smax    : max s
    -ns      : is the resolution in the s-coordinate     (e.g. 64)
    -nt      : is the resolution in the theta-coordinate (e.g. 64)
    -phi0    : toroidal angle defining a toroidal plane

    OUTPUT
    ------
    -Avol    : area in m^2 if geometrical dimensions (R,Z) are interpreted in meters."""
    tarr = np.linspace(0, 2*np.pi, nt, endpoint=True)

    if ivol==0: sarr=np.linspace(-0.999,1,ns)
    else: sarr=np.linspace(-1,    1, ns)

    # Get jacobian
    j = self.get_jacobian(lvol=ivol, sarr=sarr, tarr=tarr, zarr=np.linspace(phi0,phi0,1))

    # Integrate
    return integrate.simpson( y=integrate.simpson( y=j, x=tarr ), x=sarr )

def get_average_beta(self, ns=64, nt=64, nz=64):
    """Get beta averaged in plasma volume"""

    # Read pressure
    press = np.atleast_1d(self.input.physics.pressure) * self.input.physics.pscale

    # Create coordinate grid
    nfp = self.input.physics.Nfp
    tarr = np.linspace(0, 2*np.pi, nt)
    zarr = np.linspace(0, 2*np.pi / nfp, nz)

    # Get beta in each volume
    nvol = self.input.physics.Nvol

    vols = np.zeros((nvol,))
    betavol = np.zeros((nvol,))

    if (press==0).all(): return 0

    if nvol==1:
        sarr=np.linspace(-0.999,1, ns)
        vols = self.get_volume( 0 )
        _, _, sg, g = self.get_grid_and_jacobian_and_metric(
                        lvol=0, sarr=sarr, tarr=tarr, zarr=zarr
                    )
        Bcontrav = self.get_B(
            lvol=0, jacobian=sg, sarr=sarr, tarr=tarr, zarr=zarr
        )
        modB = self.get_modB( Bcontrav, g )

        betavol = 2 * nfp * press * integrate.simpson( 
            y=integrate.simpson( 
                y=integrate.simpson( 
                    y=sg / modB**2, x=zarr ), x=tarr ), x=sarr )
        
        return betavol / vols
    else:
        for ivol in range(0,nvol):
            if ivol==0: sarr=np.linspace(-0.999,1, ns)
            if ivol!=0: sarr=np.linspace(-1,    1, ns)

            vols[ivol] = self.get_volume( ivol )

            _, _, sg, g = self.get_grid_and_jacobian_and_metric(
                lvol=ivol, sarr=sarr, tarr=tarr, zarr=zarr
            )
            Bcontrav = self.get_B(
                lvol=ivol, jacobian=sg, sarr=sarr, tarr=tarr, zarr=zarr
            )
            modB = self.get_modB( Bcontrav, g )

            betavol[ivol] = 2 * nfp * press[ivol] * integrate.simpson( 
                y=integrate.simpson( 
                    y=integrate.simpson( 
                        y=sg / modB**2, x=zarr ), x=tarr ), x=sarr )

        return betavol.sum() / vols.sum()

def get_peak_beta(self, ns=64, nt=64, nz=64):
    if self.input.physics.pressure.size>1:
        press = self.input.physics.pressure[0] * self.input.physics.pscale
    else:
        press = self.input.physics.pressure * self.input.physics.pscale

    nfp = self.input.physics.Nfp
    tarr = np.linspace(0, 2*np.pi, nt)
    zarr = np.linspace(0, 2*np.pi / nfp, nz)
    sarr=np.linspace(-0.999,1, ns)

    vol = self.get_volume( 0 )
    _, _, sg, g = self.get_grid_and_jacobian_and_metric(
            lvol=0, sarr=sarr, tarr=tarr, zarr=zarr
        )
    Bcontrav = self.get_B(
            lvol=0, jacobian=sg, sarr=sarr, tarr=tarr, zarr=zarr
        )
    modB = self.get_modB( Bcontrav, g )

    return 2 * nfp * press * integrate.simpson( 
            y=integrate.simpson( 
                y=integrate.simpson( 
                    y=sg / modB**2, x=zarr ), x=tarr ), x=sarr ) / vol

def test_derivatives(self, lvol=0, s=0.3, t=0.4, z=0.5, delta=1e-6, tol=1e-6):
    ds = delta
    R, Z, j, g = self.get_grid_and_jacobian_and_metric(lvol, np.array([s-ds, s+ds]), np.array([t-ds, t+ds]), np.array([z-ds, z+ds]))
    Bcontra = self.get_B(lvol, j, np.array([s-ds, s+ds]), np.array([t-ds, t+ds]), np.array([z-ds, z+ds] ))
    modB = self.get_modB(Bcontra, g)
    B2 = modB ** 2
    R1, Z1, j1, g1, dj, dg = self.get_grid_and_jacobian_and_metric(lvol, np.array([s]), np.array([t]), np.array([z]), derivative=True)
    Bcontra1, dBcontra = self.get_B(lvol, j1, np.array([s]), np.array([t]), np.array([z] ), False, True, dj ) # 
    modB1, dB2 = self.get_modB(Bcontra, g, True, dBcontra, dg)

    print('Differences in dBcontra')
    print((Bcontra[1,0,0,:] - Bcontra[0,0,0,:])/ds/2 - dBcontra[0,0,0,0,:])
    print((Bcontra[0,1,0,:] - Bcontra[0,0,0,:])/ds/2 - dBcontra[0,0,0,1,:])
    print((Bcontra[0,0,1,:] - Bcontra[0,0,0,:])/ds/2 - dBcontra[0,0,0,2,:])
    print('Differences in Jacobian')
    print(np.array([j[1,0,0] - j[0,0,0], j[0,1,0] - j[0,0,0], j[0,0,1] - j[0,0,0]])/ds/2- dj[0,0,0,:])
    print('Differences in B**2')
    print(np.array([B2[1,0,0] - B2[0,0,0], B2[0,1,0] - B2[0,0,0], B2[0,0,1] - B2[0,0,0]])/ds/2- dB2[0,0,0,:])
    print('Differences in g')
    print((g[1,0,0,:,:] - g[0,0,0,:,:])/ds/2-dg[0,0,0,0,:,:])
    print((g[0,1,0,:,:] - g[0,0,0,:,:])/ds/2-dg[0,0,0,1,:,:])
    print((g[0,0,1,:,:] - g[0,0,0,:,:])/ds/2-dg[0,0,0,2,:,:])

def _validate_lsurf(lsurf:np.ndarray, mvol:int)->np.ndarray:
    """Check 
    
    Args:
        - lsurf: Interface number(s), between 1 and Mvol-1. default is np.arange(1, mvol)
        - mvol: Number of volumes
    Returns:
        - lsurf: 1d array of interface numbers
        
    Raises:
        - ValueError: if input is outside of range or wrong type (lsurf)
    """

    if lsurf is None:
        lsurf = np.arange(1,mvol)
    else:
        lsurf = np.atleast_1d(lsurf)
    if (lsurf<1).any() or (lsurf>mvol-1).any(): raise ValueError('lsurf should be in [1,mvol-1]')
    
    return lsurf

def get_surface_current_density(self, lsurf:np.ndarray, nt:int=64, nz:int=64)->typing.Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute j_surf.B on each side of the provided interfaces
    
    Args:
        - lsurf: Interface number(s), between 1 and Mvol-1. default is np.arange(1, mvol)
        - nt: Number of poloidal points
        - nz: Number of toroidal points
        
    Returns:
        - j_dot_B: mu0*j_surf.B evaluated on the grid. Shape (nsurf, nt, nz),
                   with nsurf the size of lsurf
        - tarr: theta array, size (nt,)
        - zarr: zeta array, size (nz,)
        
    Raises:
        - ValueError: if input is wrong (invalid lsurf, nt<=0, nz<=0)
    """

    mvol = self.output.Mvol
    nfp = self.input.physics.Nfp

    if mvol==1: raise ValueError('Mvol=1; no interface current!')
    lsurf = _validate_lsurf(lsurf, mvol)
    if nt<1: raise ValueError('nt should greater than zero')
    if nz<1: raise ValueError('nz should greater than zero')

    # Construct grid
    tarr = np.linspace(0, 2*np.pi, nt, endpoint=True)
    zarr = np.linspace(0, 2*np.pi/nfp, nz, endpoint=True)

    # Evaluate geometry elements
    nsurf = lsurf.size*2
    j_dot_B = np.zeros((mvol-1, 2, nt, nz))
    for s in lsurf:
        # Construct geometry elements - these are independent of the 
        # interface side
        R0, R1, R2, R3, Z0, Z1, Z2, Z3 = self.get_RZ_derivatives(
            lvol=int(s-1),
            sarr=np.asarray([1]),
            tarr=tarr,
            zarr=zarr
        )   
        et_x_ez = np.sqrt((R2*Z3)**2 + (R3*Z2)**2 + (R0*Z2)**2 + (R0*R2)**2 - 2*R2*R3*Z2*Z3)
        
        gtt = R2**2+Z2**2
        gzz = R0**2 + R3**2 + Z3**2
        gtz = R2*R3 + Z2*Z3
        g = gtt*gzz - gtz**2

        # project on each side of interface
        Bcontrav = np.zeros((2,nt,nz,3))
        for innout in [0,1]:
            # if innout=0, inner side of interface, thus vvol=s-1 and sarr=1
            # if innout=1, outer side of interface, thus vvol=s and sarr=-1
            lvol = s - np.mod(innout+1,2)
            sarr = np.asarray([-innout*2+1])

            # Get magnetic field
            Bcontrav[innout,:,:,:] = self.get_B(
                lvol=lvol,
                sarr=sarr,
                tarr=tarr,
                zarr=zarr
            )[0]
            
        Bcontrav_jump = Bcontrav[1]-Bcontrav[0]

        for innout in [0,1]:
            j_dot_B[s-1, innout] = g / et_x_ez * ( 
                Bcontrav[innout, :, :, 1]*Bcontrav_jump[:, :, 2]
              - Bcontrav[innout, :, :, 2]*Bcontrav_jump[:, :, 1]  
            )

    return j_dot_B, tarr, zarr

def get_surface_area(self, lsurf:np.ndarray=None, nt:int=64, nz:int=64):
    """Compute the surface area of a volume interface
    
    Args:
        - lsurf: Interface number(s), between 1 and Mvol-1. default is np.arange(1, mvol)
        - nt: Number of poloidal points for integration, default is 64
        - nz: Number of toroidal points for integration, default is 64
        
    Returns:
        - S: the surface area
        
    Raises:
        - ValueError: if input is wrong (invalid lsurf, nt<=0, nz<=0)
    """

    mvol = self.output.Mvol
    nfp = self.input.physics.Nfp

    if mvol==1: 
        raise ValueError('Mvol=1; no interface current!')
    lsurf = _validate_lsurf(lsurf, mvol)
    if nt<1: 
        raise ValueError('nt should greater than zero')
    if nz<1: 
        raise ValueError('nz should greater than zero')

    # Construct grid
    tarr = np.linspace(0, 2*np.pi, nt, endpoint=True)
    zarr = np.linspace(0, 2*np.pi/nfp, nz, endpoint=True)

    # Create variable for storing the surface area
    S = np.zeros((mvol-1,))

    # Loop on interfaces
    for s in lsurf:
        # Construct geometry elements
        R0, R1, R2, R3, Z0, Z1, Z2, Z3 = self.get_RZ_derivatives(
            lvol=int(s-1),
            sarr=np.asarray([1]),
            tarr=tarr,
            zarr=zarr
        )   
        #e theta x e phi
        et_x_ez = np.sqrt((R2*Z3)**2 + (R3*Z2)**2 + (R0*Z2)**2 + (R0*R2)**2 - 2*R2*R3*Z2*Z3)

        S = nfp*integrate.simps(integrate.simps(et_x_ez,zarr,axis=2),tarr,axis=1)

    return S

def get_flux_surface_average( self, lsurf, f, tarr, zarr ):
    """Returns the flux surface average of a function f. 

    For each surface lsurf, the average is made by computing the jacobian on the 
    inner side of the interface (i.e. lvol=lsurf-1, sarr=1)
    
    Args:
        - lsurf: Interface number(s), between 1 and Mvol-1. default is np.arange(1, mvol)
        - f (2D numpy array): function evaluated on a grid
        - tgrid (2D numpy array): theta grid
        - zgrid (2D numpy array): phi grid
    
    Returns>
        - fsavg (1D numpy array): The flux surface average of f on each surface
          given in lsurf
    Raises:
        - ValueError: if input is wrong (invalid lsurf)
    """
    
    lsurf = _validate_lsurf(lsurf, self.output.Mvol)

    # Get jacobian
    output = np.zeros(lsurf.shape)
    for ii, ll in enumerate(lsurf):
        sqrtg = self.get_jacobian(
            lvol=ll-1,
            sarr=np.array([1]),
            tarr=tarr,
            zarr=zarr
        )

        qrtg = np.squeeze(sqrtg)

        numerator   = integrate.simps(
            integrate.simps(np.multiply(f, qrtg), tarr, axis=0), zarr
            )
        denumerator = integrate.simps(
            integrate.simps(qrtg, tarr, axis=0), zarr
            )

        output[ii] = numerator / denumerator

    return output


## wishlist
# def get_polflux (see get_spec_polflux.m)
# def get_torcurr (see get_spec_volume_current.m)
# def plot_spec_boundary (see plot_spec_boundary.m)

def get_torflux(self,lvol=0,phi0=0,start=-0.999,send=1,ns=64,nt=64):

    Igeometry = self.input.physics.Igeometry
    if lvol==0 and Igeometry!=1 and start==-1.0:
        raise ValueError('InputError: start should be >1.0 in first volume')
    
    Mvol = self.output.Mvol
    if lvol<0 or lvol>Mvol:
        raise ValueError('InputError: Invalid lvol')


    if start<-1 or start>send:
        raise ValueError('InputError: invalid start')


    if send<start or send>1:
        raise ValueError('InputError: invalid send')


    if ns<1:
        raise ValueError('InputError: invalid ns')

    if nt<1:
        raise ValueError('InputError: invalid nt')


    # Prepare coordinate arrays
    sarr = np.linspace(start,send,ns)
    tarr = np.linspace(0,2*np.pi,nt,endpoint=True)
    zarr = np.linspace(phi0,phi0,1)
    jac = self.get_jacobian(lvol,sarr=sarr,tarr=tarr,zarr=zarr)
    Bcontrav = self.get_B(lvol,jacobian=jac,sarr=sarr,tarr=tarr,zarr=np.linspace(phi0,phi0,1))
    
    integrand = np.squeeze(Bcontrav[:,:,:,2]*jac)
    psitor = integrate.simpson( y=integrate.simpson( y=integrand,x=tarr ),x=sarr)
    
    return psitor







            



