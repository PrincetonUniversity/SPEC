#!/usr/bin/env python3
# coding: utf-8

# # Plot of |B| in (R,Z) cross-section by K. Aleynikova, ksenia.aleynikova@ipp.mpg.de, 2019

"""
@author Ksenia Aleynikova ksenia.aleynikov@ipp.mpg.de
"""

import numpy as np
import py_spec as p
import sys
import matplotlib.pyplot as plt

s = p.SPEC(sys.argv[1])

sarr    = np.linspace(-0.999,1,4)
tarr    = np.linspace(0,2*np.pi,5)
zarr    = np.linspace(0,0,6)
lvol    = 0

R, Z, jacobian, g = p.get_grid_and_jacobian_and_metric(s,lvol=lvol,sarr=sarr,tarr=tarr,zarr=zarr)
Bcontrav          = p.get_B(s,lvol=lvol,jacobian=jacobian,sarr=sarr,tarr=tarr,zarr=zarr)
modB              = p.get_modB(Bcontrav,g)

fig = plt.figure(figsize=(3*2,5*2))
ax = fig.add_subplot(111)
plot = ax.pcolormesh(R[:,:,0],Z[:,:,0],modB[:,:,0], vmin=2.8, vmax=3.5)#,linewidths=0.01/3,edgecolors='y')
ax.set_aspect('equal')
ax.set_xlabel(r"$R$")
ax.set_ylabel(r"$Z$")
cbar = fig.colorbar(plot, ax=ax)
cbar.set_label(r"$|\vec B|$", rotation=0)
plt.show()



