#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  9 17:51:56 2022
This file initialises a booz x form instance 
without changing the SPEC output file
@authors: A.Baillod, S.Guinchard
"""

import py_spec as sp
import booz_xform as bx
import numpy as np

######################################################
############### EXTRACT ALL SPECOUT DATA #############
######################################################


def init_from_spec( filename ):
    d         = sp.SPECout( filename )
    
    compute_surfs = np.array([0])
    verbose   = 1
    asym      = False  #False if stellarator symmetric
    nfp       = d.input.physics.Nfp
    mpol      = d.input.physics.Mpol+1
    ntor      = d.input.physics.Ntor
    nvol      = d.input.physics.Nvol
    mpol_nyq  = mpol
    ntor_nyq  = ntor
    mnmax     = ntor    +1 + (mpol-1)    *(2*ntor    +1)
    mnmax_nyq = ntor_nyq+1 + (mpol_nyq-1)*(2*ntor_nyq+1)
    xm        = d.output.im
    xn        = d.output.in_
    xm_nyq    = xm
    xn_nyq    = xn
    ns_in     = int(2*nvol-1) 
    s_in      = np.ndarray(ns_in, dtype = np.float64)
    s_in[0]   = d.output.tflux;
    iota      = np.array([d.output.lambdamn[1][0][0]])
    rmnc      = np.array([list(d.output.Rbc[1,:])]).transpose()
    rmns      = np.array([])
    zmnc      = np.array([])
    zmns      = np.array([list(d.output.Zbs[1,:])]).transpose()
    
    # build lmns
    lambdamn  = d.output.lambdamn[1][:].transpose()
    xms       = d.output.ims
    xns       = d.output.ins
    mns       = d.output.mns
    lmns      = np.zeros([mnmax,ns_in])
    lmnc      = np.array([])
    
    for ii in range (0,mns):
        mm = xms[ii]
        nn = xns[ii]
        
        if mm==0 and nn==0: #mode (0,0) is zero
            continue
        if mm>mpol-1 or mm<0:
            continue
        if nn<-ntor*nfp or nn>ntor*nfp:
            continue
        
        for jj in range(0,mnmax):
            if mm==xm[jj] and nn==xn[jj]:
                lmns[jj] = lambdamn[ii]
                continue
                           
    
    bsubumnc = np.array([list(d.output.Btemn[1,:])]).transpose()
    bsubvmnc = np.array([list(d.output.Bzemn[1,:])]).transpose()
    bsubumns = np.array([])
    bsubvmns = np.array([])
    mboz     = np.max(xm)
    nboz     = int(1/nfp * np.max(xn))
    aspect   = np.nan
    toroidal_flux = d.input.physics.phiedge 
    
    # modB computation 
    Nt   = d.grid.Nt
    Nz   = d.grid.Nz
    sarr = np.linspace(0,1,2)
    tarr = np.linspace(0,2*np.pi,Nt)
    zarr = np.linspace(0,2*np.pi/nfp,Nz)  
    
    Bcontrav       = d.get_B(lvol = 0, sarr = sarr, tarr = tarr, zarr = zarr)
    [R, Z, jac, g] = d.get_grid_and_jacobian_and_metric(lvol = 0, sarr = sarr, tarr= tarr, zarr = zarr)
    
    modB     = d.get_modB(Bcontrav, g)
    modBsurf = modB[1][:][:]
    modBcos  = np.zeros(np.shape(modBsurf))
    K        = 2*np.pi**2/nfp
    Bmn_     = np.zeros((mpol,2*ntor+1))
      
    for m in range (0,mpol):
        for n in range (-ntor,ntor+1):
            if m == 0 and n<0:
                continue
            for line in range (0,Nt):
                for column in range (0,Nz):
                    modBcos[line,column] = modBsurf[line][column]* np.cos(np.double(m)*tarr[line] - np.double(n)*np.double(nfp)*zarr[column])

            tmp = np.trapz(modBcos, x = tarr, axis = 0)
            Bmn_[m,n+ntor] = 1/K*np.trapz(tmp, x=zarr)
            
    
    Bmn_[0,ntor] = 1/2*Bmn_[0,ntor]
    
    
    Bmnc       = np.zeros((mnmax,ns_in))
    jj = 0 #index of surface
    for ii in range(0,mnmax):
        m = xm[ii]
        n = int(xn[ii] / nfp)
        Bmnc[ii,jj] = Bmn_[m,n+ntor]
        
    Bmns = np.array([])


    ##########################################################
    # Initialisation of a Booz_xform instance b              # 
    ##########################################################
    
    b = bx.Booz_xform()
    
    if ns_in<1:
        raise ValueError('ns has to be larger than zero')
        
    if nfp<1:
        raise ValueError('Nfp has to be larger than zero')
        
    if iota.size!=ns_in:
        raise ValueError('Iota has not the size ns_in')
        
    if xm.size!=mnmax:
        raise ValueError('xm has not the size mnmax')
    if xn.size!=mnmax:
        raise ValueError('xn has not the size mnmax')
    if xm_nyq.size!=mnmax_nyq:
        raise ValueError('xm_nyq has not the size mnmax')
    if xn_nyq.size!=mnmax_nyq:
        raise ValueError('xn_nyq has not the size mnmax')
        
    if xm[0]!=0:
        raise ValueError('xm first element is not right')
    if xn[0]!=0:
        raise ValueError('xn first element is not right')
    if xm_nyq[0]!=0:
        raise ValueError('xm_nyq first element is not right')
    if xn_nyq[0]!=0:
        raise ValueError('xn_nyq first element is not right')
        
    if xm[-1]!=mpol-1:
        raise ValueError('xm last element is not right')
    if xn[-1]!=nfp*ntor:
        raise ValueError('xn last element is not right')
    if xm_nyq[-1]!=mpol_nyq-1:
        raise ValueError('xm_nyq last element is not right')
    if xn_nyq[-1]!=nfp*ntor_nyq:
        raise ValueError('xn_nyq last element is not right')
        
    if rmnc.shape[0]!=mnmax:
        raise ValueError('Rmnc should have mnmax rows')
    if zmns.shape[0]!=mnmax:
        raise ValueError('Zmns should have mnmax rows')
    if asym:
        if rmns.shape[0]!=mnmax:
            raise ValueError('Rmns should have mnmax rows')
        if zmnc.shape[0]!=mnmax:
            raise ValueError('Zmnc should have mnmax rows')
        
    if rmnc.shape[1]!=ns_in:
        raise ValueError('Rmnc should have ns_in columns')
    if zmns.shape[1]!=ns_in:
        raise ValueError('Zmns should have ns_in columns')
    if asym:
        if rmns.shape[1]!=ns_in:
            raise ValueError('Rmns should have ns_in columns')
        if zmnc.shape[1]!=ns_in:
            raise ValueError('Zmnc should have ns_in columns')
        
    if bsubumnc.shape[0]!=mnmax:
        raise ValueError('Rmnc should have mnmax rows')
    if bsubvmnc.shape[0]!=mnmax:
        raise ValueError('Rmns should have mnmax rows')
    if asym:
        if bsubumns.shape[0]!=mnmax:
            raise ValueError('Zmnc should have mnmax rows')
        if bsubvmns.shape[0]!=mnmax:
            raise ValueError('Zmns should have mnmax rows')
        
    if bsubumnc.shape[1]!=ns_in:
        raise ValueError('Rmnc should have ns_in columns')
    if bsubvmnc.shape[1]!=ns_in:
        raise ValueError('Rmns should have ns_in columns')
    if asym:
        if bsubumns.shape[1]!=ns_in:
            raise ValueError('Zmnc should have ns_in columns')
        if bsubvmns.shape[1]!=ns_in:
            raise ValueError('Zmns should have ns_in columns')
            
    if Bmnc.shape[0]!=mnmax:
        raise ValueError('Invalid number of modes for Bmnc')
    if Bmnc.shape[1]!=ns_in:
        raise ValueError('Invalid number of surfaces for Bmnc')
    if asym:
        if Bmns.shape[0]!=mnmax:
            raise ValueError('Invalid number of modes for Bmns')
        if Bmns.shape[1]!=ns_in:
            raise ValueError('Invalid number of surfaces for Bmns')
            
    if lmns.shape[0]!=mnmax:
        raise ValueError('Invalid number of modes for lmns')
    if lmns.shape[1]!=ns_in:
        raise ValueError('Invalid number of surfaces for lmns')
    if asym:
        if lmnc.shape[0]!=mnmax:
            raise ValueError('Invalid number of modes for lmnc')
        if lmnc.shape[1]!=ns_in:
            raise ValueError('Invalid number of surfaces for lmnc')
        
    if any(compute_surfs<0):
        raise ValueError('Compute_surfs should be zero or positive')
    if any(compute_surfs>=ns_in):
        raise ValueError('Compute_surfs should be smaller than ns_in')
    
    
    
    
    
    
    b.verbose = verbose
    b.asym = 0.0;
    b.nfp = nfp
    b.mpol = mpol
    b.ntor = ntor
    b.mnmax = mnmax
    b.mpol_nyq = mpol_nyq
    b.ntor_nyq = ntor_nyq
    b.mnmax_nyq = mnmax_nyq
    b.xn = xn
    b.xm=xm
    b.xm_nyq = xm_nyq
    b.xn_nyq = xn_nyq
    b.ns_in = ns_in
    b.s_in = s_in
    b.iota = iota
    b.rmnc = rmnc
    b.rmns = rmns
    b.zmnc = zmnc
    b.zmns = zmns
    b.lmns = lmns 
    b.lmnc = lmnc
    b.bmnc = Bmnc
    b.bmns = Bmns
    b.bsubumnc = bsubumnc
    b.bsubvmnc = bsubvmnc
    b.bsubumns = bsubumns
    b.bsubvmns = bsubvmns
    b.mboz = mboz
    b.nboz = nboz
    b.compute_surfs = compute_surfs
    b.aspect = aspect
    b.toroidal_flux = toroidal_flux 
    


    return b
