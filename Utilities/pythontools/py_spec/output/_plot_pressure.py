def plot_pressure(self, normalize=True, ax=None, **kwargs):
    """Plot stepped pressure profile

    Args:
        normalize (bool, optional): Normalize the pressure with mu_0. Defaults to True.
        ax (Matplotlib axis, optional): Matplotlib axis to be plotted on. Defaults to None.
        kwargs (dict, optional): Keyword arguments. Matplotlib.pyplot.plot keyword arguments.
    """
    import numpy as np
    import matplotlib.pyplot as plt

    pressure = self.input.physics.pressure * self.input.physics.pscale
    tflux = self.output.tflux[: len(pressure)]
    if not normalize:
        #  remove  mu_0
        pressure /= 4 * np.pi * 1.0e-7
    # get axis data
    if ax is None:
        fig, ax = plt.subplots()
    plt.sca(ax)
    # set default plotting parameters
    if kwargs.get("linewidth") == None:
        kwargs.update({"linewidth": 2.0})  # prefer thicker lines
    if kwargs.get("label") == None:
        kwargs.update({"label": "self_pressure"})  # default label
    # process data
    _tflux = np.insert(tflux, 0, 0)
    _pressure = np.append(pressure, 0)
    x_tflux = np.zeros(2 * len(tflux) + 1)
    x_tflux[0::2] = _tflux
    x_tflux[1::2] = tflux
    y_pressure = np.zeros(2 * len(pressure) + 1)
    y_pressure[0::2] = _pressure
    y_pressure[1::2] = pressure
    # plot
    ax.plot(x_tflux, y_pressure, **kwargs)
    # Figure properties
    plt.xlabel("Normalized flux", fontsize=20)
    plt.ylabel("Pressure", fontsize=20)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    return
