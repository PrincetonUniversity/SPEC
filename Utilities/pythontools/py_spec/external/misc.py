#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import sys

def xy2rp(x, y):
    """Convert (x,y) to (R,phi) in polar coordinate
    """
    R = np.sqrt(x**2 + y**2)
    if   x >  0.0 and y >= 0.0 : # [0,pi/2)
        phi = np.arcsin(y/R)
    elif x <= 0.0 and y >  0.0 : # [pi/2, pi)
        phi = np.arccos(x/R)
    elif x <  0.0 and y <= 0.0 : # [pi, 3/2 pi)
        phi = np.arccos(-x/R) + np.pi
    elif x >= 0.0 and y <  0.0 : # [3/2 pi, 2pi)
        phi = np.arcsin(y/R) + 2*np.pi
    else:
        raise ValueError("Something wrong with your inputs ({:f}, {:f}).".format(x, y))
    return R, phi

def map_matrix(xx, first=True, second=True):
    """Map matrix to be complete (closed)
    Arguments:
      xx -- 2D numpy array
      first -- boolean, default: True, if increase the first dimension
      second -- boolean, default: True, if increase the second dimension

    Returns:
      new -- the new matrix with dimension increased
    """
    a, b = np.shape(xx)
    # only first
    if first and not second :
        new = np.zeros((a+1,b))
        new[0:a, 0:b] = xx[0:a, 0:b]
        new[  a, 0:b] = xx[0  , 0:b]
    # only second
    elif not first and second :
        new = np.zeros((a,b+1))
        new[0:a, 0:b] = xx[0:a, 0:b]
        new[0:a,   b] = xx[0:a, 0  ]
    # both direction
    elif first and second :
        new = np.zeros((a+1,b+1))
        new[0:a, 0:b] = xx[0:a, 0:b]
        new[  a, 0:b] = xx[0  , 0:b]
        new[0:a,   b] = xx[0:a, 0  ]
        new[  a,   b] = xx[0  , 0  ]
    # otherwise return the original matrix
    else :
        return xx
    return new

def toroidal_period(vec, nfp=1):
    """
    vec: [x,y,z] data
    Nfp: =1, toroidal number of periodicity
    """
    phi = 2*np.pi/nfp
    vec = np.atleast_2d(vec)
    new_vec = vec.copy()
    for ifp in range(nfp):
        if ifp==0:
            continue
        rotate = np.array([[  np.cos(ifp*phi), np.sin(ifp*phi), 0], \
                           [ -np.sin(ifp*phi), np.cos(ifp*phi), 0], \
                           [                0,               0, 1]])
        new_vec = np.concatenate((new_vec, np.matmul(vec, rotate)))
    return new_vec

# Print iterations progress
def print_progress(iteration, total, prefix='Progress', suffix='Complete', decimals=1, bar_length=60):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        bar_length  - Optional  : character length of bar (Int)
    """
    str_format = "{0:." + str(decimals) + "f}"
    percents = str_format.format(100 * (iteration / float(total)))
    filled_length = int(round(bar_length * iteration / float(total)))
    bar = '█' * filled_length + '-' * (bar_length - filled_length)

    sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix)),

    if iteration == total:
        sys.stdout.write('\n')
    sys.stdout.flush()
    return

# Smart way to check where to plot
def get_figure(axes=None, **kwargs):
    """
    Check where to plot
    Parameters:
        axes: matplotlib.pyplot axis, axis to plot (default: None)
        kwargs: keyword arguments
    Return:
        f, ax 
        f  : matplotlib.pyplot figure
        ax : matplotlib.pyplot axis 
    """
    import matplotlib.pyplot as plt
    if axes is None:
        # No axes provided
        f, axes = plt.subplots(**kwargs)
        '''
        f = plt.gcf()
        if len(f.axes):
            # normal situation in which existing figures should be respected and left alone
            f, axes = plt.subplots(**kwargs)
        else:
            #  made a empty figure for using
            axes = f.add_subplot(**kwargs)
         '''
    else:
        # axes = np.atleast_1d(axes)
        f = axes.get_figure()
    return f, axes

def colorbar(mappable, **kwargs):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    import matplotlib.pyplot as plt
    last_axes = plt.gca()
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(mappable, cax=cax, **kwargs)
    plt.sca(last_axes)
    return cbar

def kwargs2dict(**kwargs):
    return kwargs
