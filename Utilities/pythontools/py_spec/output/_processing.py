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
    for j in range(mn):
        if lvol > 0 or im[j] == 0:
            fac.append([sbar, 0.5 * np.ones(sarr.size)])
        else:
            fac.append(
                [sbar ** (im[j] / 2.0), (im[j] / 4.0) * sbar ** (im[j] / 2.0 - 1.0)]
            )
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
    sarr=np.linspace(0, 0, 1),
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

    fac = []
    sbar = (sarr + 1) / 2
    for j in range(mn):
        if lvol > 0 or im[j] == 0:
            fac.append([np.ones(sarr.size), np.zeros(sarr.size)])
        else:
            fac.append(
                [sbar ** (im[j] / 2.0), (im[j] / 4.0) * sbar ** (im[j] / 2.0 - 1.0)]
            )
    fac = np.array(fac)

    import numpy.polynomial.chebyshev as Cheb

    nax = np.newaxis
    # [mn,it,iz]
    ang_arg = im[:, nax, nax] * tarr[nax, :, nax] - in_[:, nax, nax] * zarr[nax, nax, :]
    cosa = np.cos(ang_arg)
    sina = np.sin(ang_arg)

    # Ch ,mn,t ,z
    c = (
        im[nax, :, nax, nax] * Azo.T[:, :, nax, nax]
        + in_[nax, :, nax, nax] * Ato.T[:, :, nax, nax]
    ) * cosa[nax, :, :, :] - (
        im[nax, :, nax, nax] * Aze.T[:, :, nax, nax]
        + in_[nax, :, nax, nax] * Ate.T[:, :, nax, nax]
    ) * sina[
        nax, :, :, :
    ]

    Bs = np.rollaxis(np.sum(fac[:, 0, nax, nax, :] * Cheb.chebval(sarr, c), axis=0), 2)

    c1 = (
        Aze.T[:, :, nax, nax] * cosa[nax, :, :, :]
        + Azo.T[:, :, nax, nax] * sina[nax, :, :, :]
    )
    Bt = np.rollaxis(
        np.sum(
            -fac[:, 0, nax, nax, :] * Cheb.chebval(sarr, Cheb.chebder(c1))
            - fac[:, 1, nax, nax, :] * Cheb.chebval(sarr, c1),
            axis=0,
        ),
        2,
    )

    c2 = (
        Ate.T[:, :, nax, nax] * cosa[nax, :, :, :]
        + Ato.T[:, :, nax, nax] * sina[nax, :, :, :]
    )
    Bz = np.rollaxis(
        np.sum(
            fac[:, 0, nax, nax, :] * Cheb.chebval(sarr, Cheb.chebder(c2))
            + fac[:, 1, nax, nax, :] * Cheb.chebval(sarr, c2),
            axis=0,
        ),
        2,
    )

    Bcontrav = np.array([Bs, Bt, Bz]) / jacobian
    return Bcontrav


# Bcontrav = get_B(s,lvol=lvol,jacobian=jacobian,sarr=sarr,tarr=tarr,zarr=zarr)


def get_modB(Bcontrav, g):
    """Input - Bcontrav has to come from get_B function"""
    modB = np.sqrt(np.einsum("iabc,jiabc,jabc->abc", Bcontrav, g, Bcontrav))
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
