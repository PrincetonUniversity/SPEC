#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 17:07:08 2020

@author: jonathan
"""

import numpy as np
import matplotlib.pyplot as plt

Z = 0.0
m = 10
n = 20
r = 0.0

if (n-m)%2==1:
    raise ValueError("n-m is not even!")

# since n-m is even, it can be divided by 2
n_m_2 = int((n-m)/2)
print("(n-m)/2 = "+str(n_m_2))

# storage for recurrence relation
nRangeZ = np.zeros([n_m_2+1])

allR = np.linspace(0.0, 1.0, 1500)
allTheta = np.linspace(0.0, 2.0*np.pi, 1000)
allZ = np.zeros_like(allR)

allX = np.zeros([len(allR), len(allTheta)])
allY = np.zeros([len(allR), len(allTheta)])
Z2d = np.zeros([len(allR), len(allTheta)])

for idx,r in enumerate(allR):
        
    # j=0: (m,m) entry
    nRangeZ[0] = r**m
    
    # j=1: (m,m+2) entry
    nRangeZ[1] = (m+2) * r**(m+2) - (m+1) * r**m
    
    for j in range(2, int(n_m_2+1)):
        myN = m+2*j
        # print("compute for n="+str(myN)+" at j="+str(j))
        
        prefac = myN/(myN**2 - m**2)
        
        f1_s1 = 4*(myN-1) * r**2
        f1_s2 = (myN-2+m)**2/(myN-2)
        f1_s3 = (myN-m)**2/myN
        f1 = (f1_s1-f1_s2-f1_s3) * nRangeZ[j-1]
        
        f2 = ((myN-2)**2 - m**2)/(myN-2) * nRangeZ[j-2]
        
        nRangeZ[j] = prefac * (f1 - f2)

    
    allZ[idx] = nRangeZ[-1]
    
    allX[idx,:] = r*np.cos(allTheta)
    allY[idx,:] = r*np.sin(allTheta)
    Z2d[idx,:] = nRangeZ[-1] * np.cos(m*allTheta)

Z = nRangeZ[-1]

print("")
print("r = "+str(r))
print(" => Z = "+str(Z))

plt.figure()
plt.plot(allR, allZ)

plt.figure()
plt.semilogx(1.0-allR, allZ)

plt.figure()
plt.pcolormesh(allX, allY, Z2d)
plt.axis("equal")
plt.colorbar()