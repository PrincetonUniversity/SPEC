def plot_kam_surface(self, ns=[], ntheta=1000, zeta=0.0, ax=None, **kwargs):
    """Plot SPEC KAM surfaces

    Args:
        ns (list, optional): List of surface index to be plotted (0 for axis, -1 for the computational boundary if applied).
                             Defaults to [] (plot all).
        zeta (float, optional): The toroidal angle where the cross-sections are plotted. Defaults to 0.0.
        ax (Matplotlib axis, optional): Matplotlib axis to be plotted on. Defaults to None.
        kwargs (dict, optional): Keyword arguments. Matplotlib.pyplot.plot keyword arguments
    Returns:
        list : list of FourSurf classes
    """
    import numpy as np
    import matplotlib.pyplot as plt

    Igeometry = self.input.physics.Igeometry
    from coilpy import FourSurf

    surfs = []
    # check if plot all
    if len(ns) == 0:
        # 0 for the axis
        ns = np.arange(self.input.physics.Nvol + self.input.physics.Lfreebound + 1)
    else:
        ns = np.atleast_1d(ns)
    # get axix data
    if ax is None:
        fig, ax = plt.subplots()
    plt.sca(ax)
    # set default plotting parameters
    if kwargs.get("label") == None:
        kwargs.update({"label": "SPEC_KAM"})  # default label
    if kwargs.get("c") == None:
        kwargs.update({"c": "red"})
    # plot all the surfaces
    if Igeometry == 3:

        for i in ns:
            _surf = FourSurf.read_spec_output(self, i)
            if i == 0:
                # plot axis as a curve
                _r, _z = _surf.rz(0.0, zeta)
                plt.scatter(_r, _z, **kwargs)
            else:
                _surf.plot(zeta=zeta, **kwargs)
            surfs.append(_surf)
        return surfs
    elif Igeometry == 2:
        for i in ns:
            _surf = FourSurf.read_spec_output(self, i)
            if i == 0:
                pass  # don't do anything for the axis
            else:
                _theta = np.arange(
                    0, 2 * np.pi + 2 * np.pi / ntheta, 2 * np.pi / ntheta
                )
                _r, _z = _surf.rz(_theta, np.ones_like(_theta) * zeta)
                plt.scatter(_r * np.cos(_theta), _r * np.sin(_theta), **kwargs)
            surfs.append(_surf)
        return surfs
    elif Igeometry == 1:
        for i in ns:
            _surf = FourSurf.read_spec_output(self, i)
            # plot axis as a curve
            _theta = np.arange(0, 2 * np.pi + 2 * np.pi / ntheta, 2 * np.pi / ntheta)
            _r, _z = _surf.rz(_theta, np.ones_like(_theta) * zeta)
            plt.scatter(_theta, _r, **kwargs)
            surfs.append(_surf)
        return surfs
