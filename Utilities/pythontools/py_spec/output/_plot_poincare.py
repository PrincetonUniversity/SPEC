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

    # extract slice corresponding to the given toroidal cutplane
    rr = self.poincare.R[:, :, toroidalIdx]
    zz = self.poincare.Z[:, :, toroidalIdx]
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
        kwargs.update({"c": "gray"})
    # make plot depending on the 'range'
    if prange == "full":
        dots = ax.scatter(rr, zz, **kwargs)
    elif prange == "upper":
        dots = ax.scatter(rr[zz >= 0], zz[zz >= 0], **kwargs)
    elif prange == "lower":
        dots = ax.scatter(rr[zz <= 0], zz[zz <= 0], **kwargs)
    else:
        raise ValueError("prange should be one of ['full'(default), 'upper', 'lower'].")
    # adjust figure properties
    plt.xlabel("R [m]", fontsize=20)
    plt.ylabel("Z [m]", fontsize=20)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.axis("equal")
    return