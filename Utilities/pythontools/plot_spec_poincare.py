#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Jonathan Schilling (jonathan.schilling@ipp.mpg.de)
"""

# plot the Poincare data
def plot_spec_poincare(s, toroidalIdx=0, plt=None):
    
    # extract slice corresponding to the given toroidal cutplane
    pltR=s.poincare.R[:,:,toroidalIdx]
    pltZ=s.poincare.Z[:,:,toroidalIdx]
    
    # create a new figure if no object to plot into was given
    newFigure = False
    if plt == None:
        import matplotlib.pyplot as plt 
        plt.figure()
        newFigure = True
    
    # do the plot
    plt.plot(pltR, pltZ, 'k.', markersize=1)
    
    # basic decoration if newly created figure
    if newFigure:
        plt.axis("equal")
        plt.xlabel("R / m")
        plt.ylabel("Z / m")
    
    # return the target plot object for further decoration
    return plt
