#!/usr/bin/env python
# coding: utf-8

# # SPEC postprocessing library by K. Aleynikova, ksenia.aleynikova@ipp.mpg.de, 2019

# Place holder for math introduction of the SPEC data: <br>
# $R_{ac} cos(m \theta - n \phi)$


"""
@author Ksenia Aleynikova ksenia.aleynikov@ipp.mpg.de
"""

import numpy as np

def get_grid_and_jacobian_and_metric(s,lvol=0,sarr=np.linspace(1,1,1),tarr=np.linspace(0,0,1),zarr=np.linspace(0,0,1)):
    
    Rac,Rbc = s.output.Rbc[lvol:lvol+2]
    Zas,Zbs = s.output.Zbs[lvol:lvol+2]
    
    mn = Rac.size #s.output.mn
    im = s.output.im
    in1 = s.output.in1
    
    sbar    = (sarr+1)/2;
    fac = []
    for j in range(mn):
        if (lvol>1 or im[j]==0):
            fac.append([sbar,0.5*np.ones(sarr.size)])
        else:
            fac.append([sbar**(im[j]/2.),(im[j]/4.)*sbar**(im[j]/2.-1.)])
    fac=np.array(fac)
    
    nax = np.newaxis
    ang_arg = im[:,nax,nax]*tarr[nax,:,nax]-in1[:,nax,nax]*zarr[nax,nax,:]
    cos = np.cos(ang_arg); sin = np.sin(ang_arg)
    dR1 = Rac[:,nax] + fac[:,0,:]*(Rbc[:,nax]-Rac[:,nax])
    dZ1 = Zas[:,nax] + fac[:,0,:]*(Zbs[:,nax]-Zas[:,nax])

    Rarr0 = np.sum(                                           dR1[:,:,nax,nax]*cos[:,nax,:,:],axis=0)
    Zarr0 = np.sum(                                           dZ1[:,:,nax,nax]*sin[:,nax,:,:],axis=0)

    Rarr1 = np.sum( fac[:,1,:,nax,nax]*(Rbc[:,nax,nax,nax]-Rac[:,nax,nax,nax])*cos[:,nax,:,:],axis=0)
    Zarr1 = np.sum( fac[:,1,:,nax,nax]*(Zbs[:,nax,nax,nax]-Zas[:,nax,nax,nax])*sin[:,nax,:,:],axis=0)

    Rarr2 = np.sum(                        -im[:,nax,nax,nax]*dR1[:,:,nax,nax]*sin[:,nax,:,:],axis=0)
    Zarr2 = np.sum(                         im[:,nax,nax,nax]*dZ1[:,:,nax,nax]*cos[:,nax,:,:],axis=0)

    Rarr3 = np.sum(                        in1[:,nax,nax,nax]*dR1[:,:,nax,nax]*sin[:,nax,:,:],axis=0)
    Zarr3 = np.sum(                       -in1[:,nax,nax,nax]*dZ1[:,:,nax,nax]*cos[:,nax,:,:],axis=0)
    
      
    jacobian = Rarr0*(Rarr2*Zarr1 - Rarr1*Zarr2) # from matlab
    
    g11 = Rarr1**2 + Zarr1**2;                #gss
    g22 = Rarr2**2 + Zarr2**2;                #gtt
    g33 = Rarr0**2 + Rarr3**2 + Zarr3**2;     #gzz
    g12 = Rarr1*Rarr2 + Zarr1*Zarr2;          #gst
    g13 = Rarr1*Rarr3 + Zarr1*Zarr3;          #gsz
    g23 = Rarr2*Rarr3 + Zarr2*Zarr3;          #gtz

    g = np.array([[g11,g12,g13],
                  [g12,g22,g23],
                  [g13,g23,g33]])
    
    #g_roll = np.rollaxis(np.rollaxis(g,0,5),0,5)
    #jacobian_from_metric = np.sqrt(np.linalg.det(g_roll))   #from Metric. Equal to np.abs of the one from matlab
    
    #WHY?
    #np.allclose(jacobian,jacobian_from_metric) -> False
    #np.allclose(np.abs(jacobian),jacobian_from_metric) -> True
    
    return Rarr0, Zarr0, jacobian, g
#R, Z, jacobian, g = get_jacobian_and_metric(s,lvol=lvol,sarr=sarr,tarr=tarr,zarr=zarr)

def get_grid_and_jacobian_and_metric_slow(s,lvol=0,sarr=np.linspace(0,0,1),tarr=np.linspace(0,0,1),zarr=np.linspace(0,0,1)):
    
    Rac,Rbc = s.output.Rbc[lvol:lvol+2]
    Zas,Zbs = s.output.Zbs[lvol:lvol+2]
    
    mn = Rac.size #s.output.mn
    im = s.output.im
    in1 = s.output.in1
    
    ns      = sarr.size
    nt      = tarr.size
    nz      = zarr.size
    
    sbar    = (sarr+1)/2;
    fac = []
    for j in range(mn):
        if (lvol>1 or im[j]==0):
            fac.append([sbar,0.5*np.ones(sarr.size)])
        else:
            fac.append([sbar**(im[j]/2.),(im[j]/4.)*sbar**(im[j]/2.-1.)])
    
    Rarr0 = np.zeros((ns,nt,nz))
    Zarr0 = np.zeros((ns,nt,nz))
    Rarr1 = np.zeros((ns,nt,nz))
    Zarr1 = np.zeros((ns,nt,nz))
    Rarr2 = np.zeros((ns,nt,nz))
    Zarr2 = np.zeros((ns,nt,nz))
    Rarr3 = np.zeros((ns,nt,nz))
    Zarr3 = np.zeros((ns,nt,nz))
    
    for j in range(mn):
        for it in range(nt):
            for iz in range(nz):
                cosa = np.cos(im[j]*tarr[it]-in1[j]*zarr[iz])
                sina = np.sin(im[j]*tarr[it]-in1[j]*zarr[iz])
                Rarr0[:,it,iz] = Rarr0[:,it,iz] +        (Rac[j] + fac[j][0]*(Rbc[j]-Rac[j]) )*cosa;
                Zarr0[:,it,iz] = Zarr0[:,it,iz] +        (Zas[j] + fac[j][0]*(Zbs[j]-Zas[j]) )*sina;

                Rarr1[:,it,iz] = Rarr1[:,it,iz] +                  fac[j][1]*(Rbc[j]-Rac[j])  *cosa;
                Zarr1[:,it,iz] = Zarr1[:,it,iz] +                  fac[j][1]*(Zbs[j]-Zas[j])  *sina;

                Rarr2[:,it,iz] = Rarr2[:,it,iz] -  im[j]*(Rac[j] + fac[j][0]*(Rbc[j]-Rac[j]) )*sina;
                Zarr2[:,it,iz] = Zarr2[:,it,iz] +  im[j]*(Zas[j] + fac[j][0]*(Zbs[j]-Zas[j]) )*cosa;

                Rarr3[:,it,iz] = Rarr3[:,it,iz] + in1[j]*(Rac[j] + fac[j][0]*(Rbc[j]-Rac[j]) )*sina;
                Zarr3[:,it,iz] = Zarr3[:,it,iz] - in1[j]*(Zas[j] + fac[j][0]*(Zbs[j]-Zas[j]) )*cosa;
    
    jacobian = Rarr0*(Rarr2*Zarr1 - Rarr1*Zarr2)
    
    g11 = Rarr1**2 + Zarr1**2;                #gss
    g22 = Rarr2**2 + Zarr2**2;                #gtt
    g33 = Rarr0**2 + Rarr3**2 + Zarr3**2;     #gzz
    g12 = Rarr1*Rarr2 + Zarr1*Zarr2;          #gst
    g13 = Rarr1*Rarr3 + Zarr1*Zarr3;          #gsz
    g23 = Rarr2*Rarr3 + Zarr2*Zarr3;          #gtz

    g = np.array([[g11,g12,g13],
                  [g12,g22,g23],
                  [g13,g23,g33]])
    
    return Rarr0, Zarr0, jacobian, g
#R_slow, Z_slow, jacobian_slow, g_slow  = get_jacobian_and_metric_slow(s,lvol=lvol,sarr=sarr,tarr=tarr,zarr=zarr)

def get_B_slow(s,lvol=0,jacobian=None,sarr=np.linspace(0,0,1),tarr=np.linspace(0,0,1),zarr=np.linspace(0,0,1)):
    
    if jacobian is None:
        R, Z, jacobian, g = get_grid_and_jacobian_and_metric_slow(s,lvol=lvol,sarr=sarr,tarr=tarr,zarr=zarr)

    #Lrad = s.input.physics.Lrad[lvol]
    Ate  = s.vector_potential.Ate[lvol]
    Aze  = s.vector_potential.Aze[lvol]
    Ato  = s.vector_potential.Ato[lvol]
    Azo  = s.vector_potential.Azo[lvol]
    
    mn = Ate.shape[0]
    im = s.output.im
    in1 = s.output.in1
    
    fac = []
    sbar    = (sarr+1)/2;
    for j in range(mn):
        if (lvol>1 or im[j]==0):
            fac.append([np.ones(sarr.size),np.zeros(sarr.size)])
        else:
            fac.append([sbar**(im[j]/2.),(im[j]/4.)*sbar**(im[j]/2.-1.)])

    import numpy.polynomial.chebyshev as Cheb
    Bs = np.zeros((sarr.size,tarr.size,zarr.size))
    Bt = np.zeros((sarr.size,tarr.size,zarr.size))
    Bz = np.zeros((sarr.size,tarr.size,zarr.size))

    for j in range(mn):
        for it in range(tarr.size):
            for iz in range(zarr.size):
                cosa = np.cos(im[j]*tarr[it]-in1[j]*zarr[iz])
                sina = np.sin(im[j]*tarr[it]-in1[j]*zarr[iz])

                c = ( (im[j]*Azo[j] + in1[j]*Ato[j])*cosa - (im[j]*Aze[j] + in1[j]*Ate[j])*sina  )
                Bs[:,it,iz] = Bs[:,it,iz] + fac[j][0]*Cheb.chebval(sarr, c);

                c1 = ( Aze[j]*cosa + Azo[j]*sina );
                Bt[:,it,iz] = Bt[:,it,iz] - ( fac[j][0]*Cheb.chebval(sarr, Cheb.chebder(c1))
                                             +fac[j][1]*Cheb.chebval(sarr, c1));

                c2 = ( Ate[j]*cosa + Ato[j]*sina );
                Bz[:,it,iz] = Bz[:,it,iz] + ( fac[j][0]*Cheb.chebval(sarr, Cheb.chebder(c2))
                                             +fac[j][1]*Cheb.chebval(sarr, c2));    
    
    Bcontrav = np.array([Bs,Bt,Bz])/jacobian
    return Bcontrav
#Bcontrav_slow = get_B_slow(s,lvol=lvol,jacobian=jacobian,sarr=sarr,tarr=tarr,zarr=zarr)

def get_B(s,lvol=0,jacobian=None,sarr=np.linspace(0,0,1),tarr=np.linspace(0,0,1),zarr=np.linspace(0,0,1)):
    
    if jacobian is None:
        R, Z, jacobian, g = get_grid_and_jacobian_and_metric(s,lvol=lvol,sarr=sarr,tarr=tarr,zarr=zarr)

    #Lrad = s.input.physics.Lrad[lvol]
    Ate  = s.vector_potential.Ate[lvol]
    Aze  = s.vector_potential.Aze[lvol]
    Ato  = s.vector_potential.Ato[lvol]
    Azo  = s.vector_potential.Azo[lvol]
    
    mn = Ate.shape[0]
    im = s.output.im
    in1 = s.output.in1
    
    fac = []
    sbar    = (sarr+1)/2;
    for j in range(mn):
        if (lvol>1 or im[j]==0):
            fac.append([np.ones(sarr.size),np.zeros(sarr.size)])
        else:
            fac.append([sbar**(im[j]/2.),(im[j]/4.)*sbar**(im[j]/2.-1.)])
    fac=np.array(fac)
    
    import numpy.polynomial.chebyshev as Cheb
    
    nax = np.newaxis
                #[mn,it,iz]
    ang_arg = im[:,nax,nax]*tarr[nax,:,nax]-in1[:,nax,nax]*zarr[nax,nax,:]
    cosa = np.cos(ang_arg); sina = np.sin(ang_arg)
        
            #Ch ,mn,t ,z                                  
    c = ((im[nax,:,nax,nax]*Azo.T[:,:,nax,nax] + in1[nax,:,nax,nax]*Ato.T[:,:,nax,nax])*cosa[nax,:,:,:] 
        -(im[nax,:,nax,nax]*Aze.T[:,:,nax,nax] + in1[nax,:,nax,nax]*Ate.T[:,:,nax,nax])*sina[nax,:,:,:])

    Bs = np.rollaxis(np.sum( fac[:,0,nax,nax,:]*Cheb.chebval(sarr, c) , axis=0),2)

    c1 = Aze.T[:,:,nax,nax]*cosa[nax,:,:,:] + Azo.T[:,:,nax,nax]*sina[nax,:,:,:]
    Bt = np.rollaxis(np.sum(-fac[:,0,nax,nax,:]*Cheb.chebval(sarr, Cheb.chebder(c1))
                            -fac[:,1,nax,nax,:]*Cheb.chebval(sarr, c1), axis=0),2)
    
    c2 = Ate.T[:,:,nax,nax]*cosa[nax,:,:,:] + Ato.T[:,:,nax,nax]*sina[nax,:,:,:]
    Bz = np.rollaxis(np.sum( fac[:,0,nax,nax,:]*Cheb.chebval(sarr, Cheb.chebder(c2))
                            +fac[:,1,nax,nax,:]*Cheb.chebval(sarr, c2), axis=0),2)
    
    Bcontrav = np.array([Bs,Bt,Bz])/jacobian
    return Bcontrav

#Bcontrav = get_B(s,lvol=lvol,jacobian=jacobian,sarr=sarr,tarr=tarr,zarr=zarr)

def get_modB_old(s,lvol=0,jacobian=None,sarr=np.linspace(0,0,1),tarr=np.linspace(0,0,1),zarr=np.linspace(0,0,1)):
    R, Z, jacobian, g = get_jacobian_and_metric(s,lvol=lvol,sarr=sarr,tarr=tarr,zarr=zarr)
    Bcontrav = get_B(s,lvol=lvol,jacobian=jacobian,sarr=sarr,tarr=tarr,zarr=zarr)
    modB = np.sqrt(np.einsum('iabc,jiabc,jabc->abc',Bcontrav,g,Bcontrav))
    return modB

def get_modB(Bcontrav,g):
    """Input - Bcontrav has to come from get_B function
    """
    modB = np.sqrt(np.einsum('iabc,jiabc,jabc->abc',Bcontrav,g,Bcontrav))
    return modB






