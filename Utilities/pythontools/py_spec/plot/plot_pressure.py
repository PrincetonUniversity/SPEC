########################################
# plot_pressure.py:
# coded by @zhucaoxiang (czhu@pppl.gov)
# adapted by @smiet (csmiet@pppl.gov)

def plot_pressure(SPEC, normalize=True, **kwargs):
    '''Plot stepped pressure profile
    Parameters:
       normalize -- Boolean, True (default). SPEC normalizes the pressure with mu_0.
                    If False, multiply mu_0 back to pressue.
    '''
    import numpy as np
    import matplotlib.pyplot as plt
    pressure = SPEC.input.physics.pressure * SPEC.input.physics.pscale
    tflux = SPEC.output.tflux[:len(pressure)]
    if not normalize :
        #  remove  mu_0
        pressure /= (4*np.pi*1.0E-7)
    # get figure and ax data
    if plt.get_fignums():
        fig = plt.gcf()
        ax = plt.gca()
    else :
        fig, ax = plt.subplots()
    # set default plotting parameters
    if kwargs.get('linewidth') == None:
        kwargs.update({'linewidth': 2.0}) # prefer thicker lines
    if kwargs.get('label') == None:
        kwargs.update({'label': 'SPEC_pressure'}) # default label
    # process data
    _tflux = np.insert(tflux, 0, 0)
    _pressure = np.append(pressure, 0)
    x_tflux = np.zeros(2*len(tflux)+1)
    x_tflux[0::2] = _tflux
    x_tflux[1::2] =  tflux
    y_pressure = np.zeros(2*len(pressure)+1)
    y_pressure[0::2] = _pressure
    y_pressure[1::2] =  pressure
    # plot
    ax.plot(x_tflux, y_pressure, **kwargs)
    # Figure properties
    plt.xlabel('Normalized flux',fontsize=20)
    plt.ylabel('Pressure',fontsize=20)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    return fig, ax


