import numpy as np


def plot_modB(
    self,
    sarr=np.linspace(-0.999, 1, 4),
    tarr=np.linspace(0, 2 * np.pi, 5),
    zarr=np.linspace(0, 0, 6),
    lvol=0,
    ax=None,
    **kwargs
):
    """[summary]

    Args:
        sarr ([type], optional): [description]. Defaults to np.linspace(-0.999, 1, 4).
        tarr ([type], optional): [description]. Defaults to np.linspace(0, 2 * np.pi, 5).
        zarr ([type], optional): [description]. Defaults to np.linspace(0, 0, 6).
        lvol (int, optional): [description]. Defaults to 0.
        ax (Matplotlib axis, optional): Matplotlib axis to be plotted on. Defaults to None.
    """
    # parse data
    R, Z, jacobian, g = self.get_grid_and_jacobian_and_metric(
        lvol=lvol, sarr=sarr, tarr=tarr, zarr=zarr
    )
    Bcontrav = self.get_B(lvol=lvol, jacobian=jacobian, sarr=sarr, tarr=tarr, zarr=zarr)
    modB = self.get_modB(Bcontrav, g)
    # get axis data
    if ax is None:
        fig, ax = plt.subplots()
    plt.sca(ax)
    plot = ax.pcolormesh(R[:, :, 0], Z[:, :, 0], modB[:, :, 0], **kwargs)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$R$")
    ax.set_ylabel(r"$Z$")
    cbar = plt.colorbar(plot, ax=ax)
    cbar.set_label(r"$|\vec B|$", rotation=0)
    return