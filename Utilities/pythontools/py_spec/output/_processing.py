import numpy as np

def get_grid_and_jacobian_and_metric(
    self,
    lvol=0,
    sarr=np.linspace(1, 1, 1),
    tarr=np.linspace(0, 0, 1),
    zarr=np.linspace(0, 0, 1),
):
    """[summary]

    Args:
        lvol (int, optional): [description]. Defaults to 0.
        sarr ([type], optional): [description]. Defaults to np.linspace(1,1,1).
        tarr ([type], optional): [description]. Defaults to np.linspace(0,0,1).
        zarr ([type], optional): [description]. Defaults to np.linspace(0,0,1).

    Returns:
        [type]: [description]
    """
    Rac, Rbc = self.output.Rbc[lvol : lvol + 2]
    Zas, Zbs = self.output.Zbs[lvol : lvol + 2]

    mn = Rac.size  # s.output.mn
    im = self.output.im
    in_ = self.output.in_

    sbar = (sarr + 1) / 2
    fac = []

    Igeometry = self.input.physics.Igeometry
    rpol = self.input.physics.rpol
    rtor = self.input.physics.rtor

    if Igeometry == 1:
        for j in range(mn):
            fac.append([sbar, 0.5 * np.ones(sarr.size)])
    elif Igeometry == 2:
        for j in range(mn):
            if lvol > 0 or im[j] == 0:
                fac.append([sbar, 0.5 * np.ones(sarr.size)])
            else:
                fac.append(
                    [sbar ** (im[j] + 1.0), (im[j] + 1.0) / 2.0 * sbar ** (im[j])]
                )
    elif Igeometry == 3:
        for j in range(mn):
            if lvol == 0 and im[j] == 0:
                fac.append([sbar ** 2, sbar])
            elif lvol == 0 and im[j] > 0:
                fac.append([sbar ** im[j], (im[j] / 2.0) * sbar ** (im[j] - 1.0)])
            else:
                fac.append([sbar, 0.5 * np.ones(sarr.size)])

    fac = np.array(fac)

    nax = np.newaxis
    ang_arg = im[:, nax, nax] * tarr[nax, :, nax] - in_[:, nax, nax] * zarr[nax, nax, :]
    cos = np.cos(ang_arg)
    sin = np.sin(ang_arg)
    dR1 = Rac[:, nax] + fac[:, 0, :] * (Rbc[:, nax] - Rac[:, nax])
    dZ1 = Zas[:, nax] + fac[:, 0, :] * (Zbs[:, nax] - Zas[:, nax])

    Rarr0 = np.sum(dR1[:, :, nax, nax] * cos[:, nax, :, :], axis=0)
    Zarr0 = np.sum(dZ1[:, :, nax, nax] * sin[:, nax, :, :], axis=0)

    Rarr1 = np.sum(
        fac[:, 1, :, nax, nax]
        * (Rbc[:, nax, nax, nax] - Rac[:, nax, nax, nax])
        * cos[:, nax, :, :],
        axis=0,
    )
    Zarr1 = np.sum(
        fac[:, 1, :, nax, nax]
        * (Zbs[:, nax, nax, nax] - Zas[:, nax, nax, nax])
        * sin[:, nax, :, :],
        axis=0,
    )

    Rarr2 = np.sum(
        -im[:, nax, nax, nax] * dR1[:, :, nax, nax] * sin[:, nax, :, :], axis=0
    )
    Zarr2 = np.sum(
        im[:, nax, nax, nax] * dZ1[:, :, nax, nax] * cos[:, nax, :, :], axis=0
    )

    Rarr3 = np.sum(
        in_[:, nax, nax, nax] * dR1[:, :, nax, nax] * sin[:, nax, :, :], axis=0
    )
    Zarr3 = np.sum(
        -in_[:, nax, nax, nax] * dZ1[:, :, nax, nax] * cos[:, nax, :, :], axis=0
    )

    if Igeometry == 1:
        jacobian = Rarr1 * rpol * rtor
        g11 = Rarr1 ** 2
        g22 = rpol ** 2 + Rarr2 ** 2
        g33 = rtor ** 2 + Rarr3 ** 2
        g12 = Rarr1 * Rarr2
        g13 = Rarr1 * Rarr3
        g23 = Rarr2 * Rarr3

    if Igeometry == 2:
        jacobian = Rarr1 * Rarr0
        g11 = Rarr1 ** 2
        g22 = Rarr2 ** 2 + Rarr0 ** 2
        g33 = Rarr3 ** 2 + 1.0
        g12 = Rarr1 * Rarr2
        g13 = Rarr1 * Rarr3
        g23 = Rarr2 * Rarr3
    elif Igeometry == 3:
        jacobian = Rarr0 * (Rarr2 * Zarr1 - Rarr1 * Zarr2)  # from matlab

        g11 = Rarr1 ** 2 + Zarr1 ** 2
        # gss
        g22 = Rarr2 ** 2 + Zarr2 ** 2
        # gtt
        g33 = Rarr0 ** 2 + Rarr3 ** 2 + Zarr3 ** 2
        # gzz
        g12 = Rarr1 * Rarr2 + Zarr1 * Zarr2
        # gst
        g13 = Rarr1 * Rarr3 + Zarr1 * Zarr3
        # gsz
        g23 = Rarr2 * Rarr3 + Zarr2 * Zarr3
        # gtz

    g = np.array([[g11, g12, g13], [g12, g22, g23], [g13, g23, g33]])

    return Rarr0, Zarr0, jacobian, g


def grid(
    self,
    lvol=0,
    sarr=np.linspace(1, 1, 1),
    tarr=np.linspace(0, 0, 1),
    zarr=np.linspace(0, 0, 1),
):

    Rarr0, Zarr0, _, _ = get_grid_and_jacobian_and_metric(self, lvol, sarr, tarr, zarr)
    return Rarr0, Zarr0


def jacobian(
    self,
    lvol=0,
    sarr=np.linspace(1, 1, 1),
    tarr=np.linspace(0, 0, 1),
    zarr=np.linspace(0, 0, 1),
):

    _, _, jacobian, _ = get_grid_and_jacobian_and_metric(self, lvol, sarr, tarr, zarr)
    return jacobian


def metric(
    self,
    lvol=0,
    sarr=np.linspace(1, 1, 1),
    tarr=np.linspace(0, 0, 1),
    zarr=np.linspace(0, 0, 1),
):

    _, _, _, g = get_grid_and_jacobian_and_metric(self, lvol, sarr, tarr, zarr)
    return g

def get_B(
    self,
    lvol=0,
    jacobian=None,
    sarr=np.linspace(1, 1, 1),
    tarr=np.linspace(0, 0, 1),
    zarr=np.linspace(0, 0, 1),
):
    Bcontra_and_der = get_B_and_der(self, lvol, jacobian, sarr, tarr, zarr)
    return Bcontra_and_der[0:3]

def get_s_der_B(
    self,
    lvol=0,
    jacobian=None,
    sarr=np.linspace(1, 1, 1),
    tarr=np.linspace(0, 0, 1),
    zarr=np.linspace(0, 0, 1),
):
    Bcontra_and_der = get_B_and_der(self, lvol, jacobian, sarr, tarr, zarr)
    return Bcontra_and_der[3:6]

def get_t_der_B(
    self,
    lvol=0,
    jacobian=None,
    sarr=np.linspace(1, 1, 1),
    tarr=np.linspace(0, 0, 1),
    zarr=np.linspace(0, 0, 1),
):
    Bcontra_and_der = get_B_and_der(self, lvol, jacobian, sarr, tarr, zarr)
    return Bcontra_and_der[6:9]

def get_z_der_B(
    self,
    lvol=0,
    jacobian=None,
    sarr=np.linspace(1, 1, 1),
    tarr=np.linspace(0, 0, 1),
    zarr=np.linspace(0, 0, 1),
):
    Bcontra_and_der = get_B_and_der(self, lvol, jacobian, sarr, tarr, zarr)
    return Bcontra_and_der[9:12]


def get_B_and_der(
    self,
    lvol=0,
    jacobian=None,
    sarr=np.linspace(1, 1, 1),
    tarr=np.linspace(0, 0, 1),
    zarr=np.linspace(0, 0, 1),
):

    if jacobian is None:
        R, Z, jacobian, g = get_grid_and_jacobian_and_metric(
            self, lvol=lvol, sarr=sarr, tarr=tarr, zarr=zarr
        )

    # Lrad = s.input.physics.Lrad[lvol]


    Ate = self.vector_potential.Ate[lvol]
    Aze = self.vector_potential.Aze[lvol]
    Ato = self.vector_potential.Ato[lvol]
    Azo = self.vector_potential.Azo[lvol]

    mn = Ate.shape[0]
    im = self.output.im
    in_ = self.output.in_
    Lrad = self.input.physics.Lrad[lvol]
    Mpol = self.input.physics.Mpol

    fac = []
    sbar = (sarr + 1) / 2

    nax = np.newaxis
    # [mn,it,iz]
    ang_arg = im[:, nax, nax] * tarr[nax, :, nax] - in_[:, nax, nax] * zarr[nax, nax, :]
    cosa = np.cos(ang_arg)
    sina = np.sin(ang_arg)

    LZernike = self.input.physics.Igeometry > 1 and lvol == 0

    if not LZernike:
        # make basis recombination for cheb
        lcoeff = np.arange(0, Lrad + 1) + 1

        Ate = Ate / lcoeff[None, :]
        Aze = Aze / lcoeff[None, :]
        Ato = Ato / lcoeff[None, :]
        Azo = Azo / lcoeff[None, :]

        Ate[:, 0] = np.sum(Ate * (-1.0) ** lcoeff[None, :], 1)
        Aze[:, 0] = np.sum(Aze * (-1.0) ** lcoeff[None, :], 1)
        Ato[:, 0] = np.sum(Ato * (-1.0) ** lcoeff[None, :], 1)
        Azo[:, 0] = np.sum(Azo * (-1.0) ** lcoeff[None, :], 1)

    c = (
          im[nax, :, nax, nax]  * Azo.T[:, :, nax, nax]
        + in_[nax, :, nax, nax] * Ato.T[:, :, nax, nax]
    ) * cosa[nax, :, :, :] - (
          im[nax, :, nax, nax]  * Aze.T[:, :, nax, nax]
        + in_[nax, :, nax, nax] * Ate.T[:, :, nax, nax]
    ) * sina[nax, :, :, :]

    ct = im[nax, :, nax, nax] * ( - (
          im[nax, :, nax, nax]  * Azo.T[:, :, nax, nax]
        + in_[nax, :, nax, nax] * Ato.T[:, :, nax, nax]
    ) * sina[nax, :, :, :] - (
          im[nax, :, nax, nax]  * Aze.T[:, :, nax, nax]
        + in_[nax, :, nax, nax] * Ate.T[:, :, nax, nax]
    ) * cosa[nax, :, :, :] )

    cz = in_[nax, :, nax, nax] * ( (
          im[nax, :, nax, nax]  * Azo.T[:, :, nax, nax]
        + in_[nax, :, nax, nax] * Ato.T[:, :, nax, nax]
    ) * sina[nax, :, :, :] + (
          im[nax, :, nax, nax]  * Aze.T[:, :, nax, nax]
        + in_[nax, :, nax, nax] * Ate.T[:, :, nax, nax]
    ) * cosa[nax, :, :, :] )


    c1 = (
          Aze.T[:, :, nax, nax] * cosa[nax, :, :, :]
        + Azo.T[:, :, nax, nax] * sina[nax, :, :, :])
    
    c1t = im[nax, :, nax, nax] * (
        - Aze.T[:, :, nax, nax] * sina[nax, :, :, :]
        + Azo.T[:, :, nax, nax] * cosa[nax, :, :, :])

    c1z = in_[nax, :, nax, nax] * (
        + Aze.T[:, :, nax, nax] * sina[nax, :, :, :]
        - Azo.T[:, :, nax, nax] * cosa[nax, :, :, :])


    c2 = (
          Ate.T[:, :, nax, nax] * cosa[nax, :, :, :]
        + Ato.T[:, :, nax, nax] * sina[nax, :, :, :])

    c2t = im[nax, :, nax, nax] * (
        - Ate.T[:, :, nax, nax] * sina[nax, :, :, :]
        + Ato.T[:, :, nax, nax] * cosa[nax, :, :, :])

    c2z = in_[nax, :, nax, nax] * (
        + Ate.T[:, :, nax, nax] * sina[nax, :, :, :]
        - Ato.T[:, :, nax, nax] * cosa[nax, :, :, :])


    if LZernike:
        # Zernike polynomial being used
        from ._processing import _get_zernike

        zernike, dzernike, d2zernike = _get_zernike(sbar, Lrad, Mpol)

        Bs   = np.sum( zernike[:, :, im, None, None] * c[None, :, :, :, :],  axis=(1, 2))
        dsBs = np.sum(dzernike[:, :, im, None, None] * c[None, :, :, :, :],  axis=(1, 2))
        dtBs = np.sum( zernike[:, :, im, None, None] * ct[None, :, :, :, :], axis=(1, 2))
        dzBs = np.sum( zernike[:, :, im, None, None] * cz[None, :, :, :, :], axis=(1, 2))

        Bt   = -np.sum( dzernike[:, :, im, None, None] * c1[None, :, :, :, :],  axis=(1, 2))
        dsBt = -np.sum(d2zernike[:, :, im, None, None] * c1[None, :, :, :, :],  axis=(1, 2))
        dtBt = -np.sum( dzernike[:, :, im, None, None] * c1t[None, :, :, :, :], axis=(1, 2))
        dzBt = -np.sum( dzernike[:, :, im, None, None] * c1z[None, :, :, :, :], axis=(1, 2))

        Bz   = np.sum( dzernike[:, :, im, None, None] * c2[None, :, :, :, :],  axis=(1, 2))
        dsBz = np.sum(d2zernike[:, :, im, None, None] * c2[None, :, :, :, :],  axis=(1, 2))
        dtBz = np.sum( dzernike[:, :, im, None, None] * c2t[None, :, :, :, :], axis=(1, 2))
        dzBz = np.sum( dzernike[:, :, im, None, None] * c2z[None, :, :, :, :], axis=(1, 2))

    else:
        # Chebyshev being used
        import numpy.polynomial.chebyshev as Cheb

        Bs   = np.rollaxis(np.sum(Cheb.chebval(sarr, c),               axis=0), 2)
        dsBs = np.rollaxis(np.sum(Cheb.chebval(sarr, Cheb.chebder(c)), axis=0), 2)
        dtBs = np.rollaxis(np.sum(Cheb.chebval(sarr, ct),              axis=0), 2)
        dzBs = np.rollaxis(np.sum(Cheb.chebval(sarr, cz),              axis=0), 2)

        Bt   = -np.rollaxis(np.sum(Cheb.chebval(sarr, Cheb.chebder(c1)),      axis=0),2)
        dsBt = -np.rollaxis(np.sum(Cheb.chebval(sarr, Cheb.chebder(c1,m=2)),  axis=0),2)
        dtBt = -np.rollaxis(np.sum(Cheb.chebval(sarr, Cheb.chebder(c1t)),     axis=0),2)
        dzBt = -np.rollaxis(np.sum(Cheb.chebval(sarr, Cheb.chebder(c1z)),     axis=0),2)

        Bz   = np.rollaxis(np.sum(Cheb.chebval(sarr, Cheb.chebder(c2)),      axis=0),2)
        dsBz = np.rollaxis(np.sum(Cheb.chebval(sarr, Cheb.chebder(c2,m=2)),  axis=0),2)
        dtBz = np.rollaxis(np.sum(Cheb.chebval(sarr, Cheb.chebder(c2t)),     axis=0),2)
        dzBz = np.rollaxis(np.sum(Cheb.chebval(sarr, Cheb.chebder(c2z)),     axis=0),2)


    Bcontrav_and_der = np.array([Bs, Bt, Bz, dsBs, dsBt, dsBz, dtBs, dtBt, dtBz, dzBs, dzBt, dzBz]) / jacobian
    return Bcontrav_and_der


def get_modB(self, Bcontrav, g):
    """Input - Bcontrav has to come from get_B function"""
    modB = np.sqrt(np.einsum("iabc,jiabc,jabc->abc", Bcontrav, g, Bcontrav))
    return modB

def get_B_covariant(self, Bcontrav, g):
    """Get covariant component of B"""
    Bco = np.einsum("iabc,jiabc->jabc", Bcontrav, g)
    return Bco    


def _get_zernike(sarr, lrad, mpol):
    """
    Get the value of the zernike polynomials, their first and second derivatives
    Adapted from basefn.f90
    """

    ns = sarr.size

    r = (sarr + 1.0) / 2
    rm = np.ones_like(r)  # r to the power of m'th
    rm1 = np.zeros_like(r)  # r to the power of m-1'th
    rm2 = np.zeros_like(r)  # r to the power of m-2'th

    zernike = np.zeros(shape=(ns, lrad + 1, mpol + 1), dtype=np.float64)
    dzernike = np.zeros_like(zernike)
    d2zernike = np.zeros_like(zernike)

    for m in range(mpol + 1):
        if lrad >= m:
            zernike[:, m, m]   = rm
            dzernike[:, m, m]  = m * rm1
            d2zernike[:, m, m] = m * (m-1) * rm2

        if lrad >= m + 2:
            zernike[:,   m + 2, m] = float(m + 2) * rm * r**2          - float(m + 1) * rm
            dzernike[:,  m + 2, m] = (float(m + 2)**2 * rm * r         - float((m + 1) * m) * rm1)
            d2zernike[:, m + 2, m] = (float((m + 2)**2 * (m + 1)) * rm - float((m + 1) * m * (m - 1)) * rm2)

        for n in range(m + 4, lrad + 1, 2):
            factor1 = float(n) / float(n ** 2 - m ** 2)
            factor2 = float(4 * (n - 1))
            factor3 = float((n - 2 + m) ** 2) / float(n - 2) + float((n - m) ** 2) / float(n)
            factor4 = float((n - 2) ** 2 - m ** 2) / float(n - 2)

            zernike[:, n, m] = factor1 * (
                (factor2 * r ** 2 - factor3) * zernike[:, n - 2, m]
                - factor4 * zernike[:, n - 4, m]
            )
            dzernike[:, n, m] = factor1 * (
                2 * factor2 * r * zernike[:, n - 2, m]
                + (factor2 * r ** 2 - factor3) * dzernike[:, n - 2, m]
                - factor4 * dzernike[:, n - 4, m]
            )
            d2zernike[:, n, m] = factor1 * (
                2 * factor2 * zernike[:, n - 2, m]
                + 4 * factor2 * r * dzernike[:, n - 2, m]
                + (factor2 * r ** 2 - factor3) * d2zernike[:, n - 2, m]
                - factor4 * d2zernike[:, n - 4, m]
            )

        rm2 = rm1
        rm1 = rm
        rm  = rm * r

    for n in range(2, lrad + 1, 2):
        zernike[:, n, 0] = zernike[:, n, 0] - (-1) ** (n / 2)

    if mpol >= 1:
        for n in range(3, lrad + 1, 2):
            zernike[:, n, 1] = (
                zernike[:, n, 1] - (-1) ** ((n - 1) / 2) * float((n + 1) / 2) * r
            )
            dzernike[:, n, 1] = dzernike[:, n, 1] - (-1) ** ((n - 1) / 2) * float(
                (n + 1) / 2
            )


    for m in range(mpol + 1):
        for n in range(m, lrad + 1, 2):
            zernike[:, n, m]   = zernike[:, n, m]   / float(n + 1)
            dzernike[:, n, m]  = dzernike[:, n, m]  / float(n + 1)
            d2zernike[:, n, m] = d2zernike[:, n, m] / float(n + 1)

    dzernike  = dzernike  * 0.5  # to account for the factor of half in sbar = (1+s)/2
    d2zernike = d2zernike * 0.25 # to account for the factor of half in sbar = (1+s)/2

    return zernike, dzernike, d2zernike
