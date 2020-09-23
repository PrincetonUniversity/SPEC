def plot_iota(self, xaxis="R", yaxis="i", ax=None, **kwargs):
    """Iota plots
    @param xaxis the choose of the xaxis, 'R'(default) or 's'
    @param yaxis the choose of the yaxis, 'i'(default) or 'q'
    @param ax (Matplotlib axis, optional): Matplotlib axis to be plotted on. Defaults to None.
    @param kwargs (dict, optional): keyword arguments. Matplotlib.pyplot.scatter keyword arguments.
    @raises ValueError: xaxis or yaxis illegal
    @returns pyplot.scatter: Matplotlib.pyplot.scatter returns
    """
    import matplotlib.pyplot as plt
    import matplotlib.lines as mlines
    import numpy as np

    if ax is None:
        fig, ax = plt.subplots()
    plt.sca(ax)
    # set default plotting parameters
    # use dots
    if kwargs.get("marker") == None:
        kwargs.update({"marker": "*"})
    # use gray color
    if kwargs.get("c") == None:
        pass
        # kwargs.update({"c": "gray"})

    if xaxis == "s":
        xdata = self.transform.fiota[0, :]
        ydata = self.transform.fiota[1, :]
        xlabel = r"s"
    elif xaxis == "R":
        xdata = self.poincare.R[:, 0, 0]
        ydata = self.transform.fiota[1, :]
        xlabel = r"R"
    else:
        raise ValueError("xaxis should be one of ['R', 's'].")

    dots = ax.scatter(xdata, ydata, **kwargs)

    if yaxis == "i":
        ylabel = r"$\iota$"
    elif yaxis == "q":
        ydata = 1.0 / ydata
        ylabel = r"q"
    else:
        raise ValueError("yaxis should be one of ['i', 'q'].")

    plt.xlabel(xlabel, fontsize=20)
    plt.ylabel(ylabel, fontsize=20)

    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)

    return