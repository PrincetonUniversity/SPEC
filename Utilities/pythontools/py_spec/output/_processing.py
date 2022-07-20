import numpy as np


def get_grid_and_jacobian_and_metric(
    self,
    lvol=0,
    sarr=np.linspace(1, 1, 1),
    tarr=np.linspace(0, 0, 1),
    zarr=np.linspace(0, 0, 1),
    input1D=False,
    derivative=False,
):
    """!Compute the metric and Jacobian on a given grid

    @param lvol (int, optional): The SPEC volume of interest, starting from 0. Defaults to 0.
    @param sarr (1D numpy array, optional): The s grid. Defaults to np.linspace(1,1,1).
    @param tarr (1D numpy array, optional): The \f$\theta\f$ grid. Defaults to np.linspace(0,0,1).
    @param zarr (1D numpy array, optional): The \f$\zeta\f$ grid. Defaults to np.linspace(0,0,1).
    @param input1D If sarr, tarr and zarr should be treated as a grid in 3D or just 1D input. Defaults to False
    @param derivative (bool, optional): If the derivatives of jacobian and \f$g_{ij}\f$ is needed.

    @returns Rarr0, Zarr0, jacobian, g, [djacobian, dg]: \f$R, Z, J, g_{ij}\f$. If derivative==True also return the derivative of \f$J, g_{ij}\f$ (derivative is the first dimension).
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
    sbar = np.divide(np.add(sarr, 1.0), 2.0)
    fac = []

    Igeometry = self.input.physics.Igeometry
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

    # If the derivative of g and jacobian is needed
    if derivative:
        Rarr11 = np.sum(fac[2] * (Rbc - Rac) * cos, axis=0)
        Rarr12 = np.sum(-im * fac[1] * (Rbc - Rac) * sin, axis=0)
        Rarr13 = np.sum(in_ * fac[1] * (Rbc - Rac) * sin, axis=0)
        Rarr22 = np.sum(-(im ** 2) * dR1 * cos, axis=0)
        Rarr23 = np.sum(im * in_ * dR1 * cos, axis=0)
        Rarr33 = np.sum(-(in_ ** 2) * dR1 * cos, axis=0)

        dRarr = np.array(
            [
                [Rarr11, Rarr12, Rarr13],
                [Rarr12, Rarr22, Rarr23],
                [Rarr13, Rarr23, Rarr33],
            ]
        )

        if Igeometry == 3:
            Zarr11 = np.sum(fac[2] * (Zbs - Zas) * sin, axis=0)
            Zarr12 = np.sum(im * fac[1] * (Zbs - Zas) * cos, axis=0)
            Zarr13 = np.sum(-in_ * fac[1] * (Zbs - Zas) * cos, axis=0)
            Zarr22 = np.sum(-(im ** 2) * dZ1 * sin, axis=0)
            Zarr23 = np.sum(im * in_ * dZ1 * sin, axis=0)
            Zarr33 = np.sum(-(in_ ** 2) * dZ1 * sin, axis=0)

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


def grid(
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


def jacobian(
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


def metric(
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
    """!Compute the contravariant components of the magnetic field \f$(B^s, B^\theta, B^\zeta)\f$

    @param lvol (int, optional): The SPEC volume of interest, starting from 0. Defaults to 0.
    @param jacobian(numpy array, optional): if jacobian is already computed, provide it here
    @param sarr (1D numpy array, optional): The s grid. Defaults to np.linspace(1,1,1).
    @param tarr (1D numpy array, optional): The \f$\theta\f$ grid. Defaults to np.linspace(0,0,1).
    @param zarr (1D numpy array, optional): The \f$\zeta\f$ grid. Defaults to np.linspace(0,0,1).
    @param input1D If sarr, tarr and zarr should be treated as a grid in 3D or just 1D input. Defaults to False
    @param derivative (bool, optional): If the derivatives is needed.

    @returns Bcontrav, [dBcontrav]: \f$(B^s, B^\theta, B^\zeta)\f$. If derivative==True also return the derivative (derivative is the first dimension).
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


def get_B_covariant(self, Bcontrav, g, derivative=False):
    """Get covariant component of B"""
    Bco = np.einsum("...i,...ji->...j", Bcontrav, g)
    return Bco

def test_derivatives(self, lvol=0, s=0.3, t=0.4, z=0.5, delta=1e-6, tol=1e-6):
    ds = delta
    R, Z, j, g = self.get_grid_and_jacobian_and_metric(lvol, np.array([s-ds, s+ds]), np.array([t-ds, t+ds]), np.array([z-ds, z+ds]))
    Bcontra = self.get_B(lvol, j, np.array([s-ds, s+ds]), np.array([t-ds, t+ds]), np.array([z-ds, z+ds] ))
    modB = self.get_modB(Bcontra, g)
    B2 = modB ** 2
    R1, Z1, j1, g1, dj, dg = self.get_grid_and_jacobian_and_metric(lvol, np.array([s]), np.array([t]), np.array([z]), derivative=True)
    Bcontra1, dBcontra = self.get_B(lvol, j1, np.array([s]), np.array([t]), np.array([z] ), False, True, dj )
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
