########################################
# plot_kam_surface.py
# coded by @zhucaoxiang (czhu@pppl.gov)
# adapted by @smiet (csmiet@pppl.gov)


def plot_kam_surface(SPEC, ns=None, zeta=0.0, **kwargs):
    '''Plot SPEC KAM surfaces
    parameters:
        ns -- None (default) or integer list. list of surface index to be plotted.
             (0 for axis, -1 for the computational boundary if applied).
        zeta -- float, defaultr: 0.0. The toroidal angle where the cross-sections are plotted.
        **kwargs -- keyword arguments. fourier_surface.plot keyword arguments.
    returns:
        surfs -- list of fourier_surface classes
    '''
    # use fourier_surface, can be found at
    # https://github.com/zhucaoxiang/CoilPy/blob/master/surface.py
    from .. import math
    import matplotlib.pyplot as plt
    import numpy as np
    surfs = []
    # check if plot all
    if ns is None:
        # 0 for the axis
        ns = np.arange(SPEC.input.physics.Nvol+SPEC.input.physics.Lfreebound+1)
    else:
        ns = np.atleast_1d(ns)
    # set default plotting parameters
    if kwargs.get('label') is None:
        kwargs.update({'label': 'SPEC_KAM'})  # default label
    # plot all the surfaces
    for i in ns:
        _surf = math.fourier_surface.read_spec_output(SPEC, i)  #look into which spec fourier_surface needs.
        if i == 0:
            # plot axis as a curve
            _r, _z = _surf.rz(0.0, zeta)
            plt.scatter(_r, _z, **kwargs)
        else:
            _surf.plot(zeta=zeta, **kwargs)
        surfs.append(_surf)
    fig = plt.gcf()
    ax = plt.gca()
    return fig, ax
