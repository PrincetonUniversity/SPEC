########################################
# plot_poincare.py
# coded by @zhucaoxiang (czhu@pppl.gov)
# adapted by @smiet (csmiet@pppl.gov)


def plot_poincare(SPEC, toroidalIdx=0, prange='full', **kwargs):
    '''Poincare plots
    parameters:
        toroidalIdx -- int, default: 0. The index of toroidal cross-section to be plotted.
        prange -- str, ['full'(default), 'upper', 'lower']. Range of plotted points.
        **kwargs  -- keyword arguments. Matplotlib.pyplot scatter keyword arguments.

    return:
        dots -- Matplotlib.pyplot scatter class
    '''
    import matplotlib.pyplot as plt
    import matplotlib.lines as mlines
    import numpy as np
    # extract slice corresponding to the given toroidal cutplane
    if (SPEC.input.physics.Igeometry==3):
        rr = SPEC.poincare.R[:, :, toroidalIdx]
        zz = SPEC.poincare.Z[:, :, toroidalIdx]
    elif (SPEC.input.physics.Igeometry==1):
        rr = np.mod(SPEC.poincare.t[:, :, toroidalIdx],np.pi*2)
        zz = SPEC.poincare.R[:, :, toroidalIdx]
    # get current figure or build new one;
    if plt.get_fignums():
        fig = plt.gcf()
        ax = plt.gca()
    else:
        fig, ax = plt.subplots()
    # set default plotting parameters
    # use dots
    if kwargs.get('marker') is None:
        kwargs.update({'marker': '.'})
    # use gray color
    if kwargs.get('c') is None:
        kwargs.update({'c': 'gray'})
    # make plot depending on the 'range'
    if prange == 'full':
        dots = ax.scatter(rr, zz, **kwargs)
    elif prange == 'upper':
        dots = ax.scatter(rr[zz >= 0], zz[zz >= 0], **kwargs)
    elif prange == 'lower':
        dots = ax.scatter(rr[zz <= 0], zz[zz <= 0], **kwargs)
    else:
        raise ValueError("prange should be one of ['full'(default), 'upper', 'lower'].")
    # adjust figure properties
    plt.xlabel('R [m]', fontsize=20)
    plt.ylabel('Z [m]', fontsize=20)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    if (SPEC.input.physics.Igeometry==1):
        pass
    else:
        plt.axis('equal')
    return fig, ax
