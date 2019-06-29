# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt

from spec import SPEC

filename = "/home/jonathan/Uni/04_PhD/00_programs/SPEC/SPEC/InputFiles/TestCases/G3V01L0Fi.002.h5"

s=SPEC(filename)

toroidalIdx = 40


pltR=s.poincare.R[:,:,toroidalIdx]
pltZ=s.poincare.Z[:,:,toroidalIdx]


plt.figure()
plt.plot(pltR, pltZ, 'k.', ms=1)
plt.axis("equal")
plt.xlabel("R / m")
plt.ylabel("Z / m")


