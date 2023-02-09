#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  9 13:32:45 2022

This file enables to test if the modB
Fourier modes were correctly implemented 
in Python, in order to pass them as a 
Booz_X_Form input

@author: S.Guinchard
"""
import py_spec as sp
import booz_xform as bx
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5 


d         = sp.SPECout('Filename.sp.h5')

verbose   = 1
asym      = False 
nfp       = d.input.physics.Nfp
mpol      = d.input.physics.Mpol
ntor      = d.input.physics.Ntor
nvol      = d.input.physics.Nvol
mnmax     = d.output.mn
mpol_nyq  = mpol
ntor_nyq  = ntor
mnmax_nyq = mnmax
xm        = d.output.im
xn        = d.output.in_
xm_nyq    = xm
xn_nyq    = xn
ns_in     = int(2*nvol-1) 
s_in      = np.ndarray(ns_in, dtype = np.float64)
s_in[0]   = d.output.tflux;
iota      = d.transform.fiota[1][:]
rmnc      = d.output.Rbc
rmns      = d.output.Rbs
zmnc      = d.output.Zbc
zmns      = d.output.Zbs
lmns      = d.output.lambdamn[1][:]

bsubumnc = d.output.Btemn
bsubvmnc = d.output.Bzemn 
mboz     = np.max(xm)
nboz     = int(1/nfp * np.max(xn))
aspect   = np.nan
toroidal_flux = d.input.physics.phiedge 

######## MOD B ########

Nt   = d.grid.Nt
Nz   = d.grid.Nz
sarr = np.linspace(0,1,2)
tarr = np.linspace(0,2*np.pi,Nt)
zarr = np.linspace(0,2*np.pi/nfp,Nz)  

Bcontrav       = d.get_B(lvol = 0, sarr = sarr, tarr = tarr, zarr = zarr)
[R, Z, jac, g] = d.get_grid_and_jacobian_and_metric(lvol = 0, sarr = sarr, tarr= tarr, zarr = zarr)
modB           = d.get_modB(Bcontrav, g)
modBsurf       = modB[1][:][:]
modBcos        = np.zeros(np.shape(modBsurf))
K              = 2*np.pi**2/nfp
Bmn_           = np.zeros((mpol+1,2*ntor+1))
X,Y            = np.meshgrid(tarr, zarr)
    


for m in range (0,mpol):
    for n in range (-ntor,ntor):
        if m == 0 and n<0:
            continue
        
        modBcos = modBsurf * np.cos(np.double(m)*X - np.double(n)*np.double(nfp)*Y)
        tmp     = np.trapz(modBcos, x=tarr, axis = 0)
        Bmn_[m,n+ntor] = 1/K*np.trapz(tmp, x=zarr) 
        
Bmn_[0,ntor] = 1/2*Bmn_[0,ntor]


Bmn       = np.zeros((mnmax_nyq,))
Bmn_trunc = np.zeros((mpol,2*ntor+1))

for jj in range(0,mpol):
	Bmn_trunc[jj][:] = Bmn_[jj+1][:]
Bmn_trunc = np.reshape(Bmn_trunc, np.size(Bmn_trunc))

for ii in range(-ntor,0):
	Bmn[ii+ntor] = Bmn_[0][ii+2*ntor]


for ll in range(0,np.size(Bmn_trunc)):
	Bmn[ll+ntor+1] = Bmn_trunc[ll]  

	
print('Bmn_ = ', Bmn_)

