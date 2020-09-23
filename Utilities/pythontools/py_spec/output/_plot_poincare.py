def plot_poincare(self, toroidalIdx=0, prange="full", ax=None, **kwargs):
    """Poincare plots

    Args:
        toroidalIdx (int, optional): The index of toroidal cross-section to be plotted. Defaults to 0.
        prange (str, optional): Range of plotted points, one of ['full', 'upper', 'lower']. Defaults to 'full'.
        ax (Matplotlib axis, optional): Matplotlib axis to be plotted on. Defaults to None.
        kwargs (dict, optional): keyword arguments. Matplotlib.pyplot.scatter keyword arguments.
    Raises:
        ValueError: prange should be one of ['full', 'upper', 'lower']

    Returns:
        pyplot.scatter: Matplotlib.pyplot.scatter returns
    """
    import matplotlib.pyplot as plt
    import matplotlib.lines as mlines
    import numpy as np

    # extract slice corresponding to the given toroidal cutplane
    Igeometry = self.input.physics.Igeometry
    if Igeometry == 3:
        rr = self.poincare.R[:, :, toroidalIdx]
        zz = self.poincare.Z[:, :, toroidalIdx]
    elif Igeometry == 1:
        rr = np.mod(self.poincare.t[:, :, toroidalIdx], np.pi * 2)
        zz = self.poincare.R[:, :, toroidalIdx]
    elif Igeometry == 2:
        rr = self.poincare.R[:, :, toroidalIdx] * np.cos(
            self.poincare.t[:, :, toroidalIdx]
        )
        zz = self.poincare.R[:, :, toroidalIdx] * np.sin(
            self.poincare.t[:, :, toroidalIdx]
        )
    # get axix data
    if ax is None:
        fig, ax = plt.subplots()
    plt.sca(ax)
    # set default plotting parameters
    # use dots
    if kwargs.get("marker") == None:
        kwargs.update({"marker": "."})
    # use gray color
    if kwargs.get("c") == None:
        pass
    # size of marker
    if kwargs.get("s") == None:
        kwargs.update({"s": 0.3})
        # kwargs.update({"c": "gray"})
    # make plot depending on the 'range'
    if prange == "full":
        nptrj = rr.shape[0]
        for ii in range(nptrj):
            dots = ax.scatter(rr[ii, :], zz[ii, :], **kwargs)
    elif prange == "upper":
        dots = ax.scatter(rr[zz >= 0], zz[zz >= 0], **kwargs)
    elif prange == "lower":
        dots = ax.scatter(rr[zz <= 0], zz[zz <= 0], **kwargs)
    else:
        raise ValueError("prange should be one of ['full'(default), 'upper', 'lower'].")
    # adjust figure properties
    if self.input.physics.Igeometry == 3:
        plt.xlabel("R [m]", fontsize=20)
        plt.ylabel("Z [m]", fontsize=20)
        plt.axis("equal")
    if self.input.physics.Igeometry == 2:
        plt.xlabel("X [m]", fontsize=20)
        plt.ylabel("Y [m]", fontsize=20)
        plt.axis("equal")
    if self.input.physics.Igeometry == 1:
        plt.ylabel("R [m]", fontsize=20)
        plt.xlabel(r"$\theta$", fontsize=20)
        plt.xlim([0, np.pi])
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)

    return