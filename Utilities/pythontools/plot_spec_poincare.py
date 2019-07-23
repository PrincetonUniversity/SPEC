# -*- coding: utf-8 -*-

# W7-X OP1.1 limiter, 1 volume
filename = "/home/jonathan/Uni/04_PhD/00_programs/SPEC/SPEC/InputFiles/TestCases/G3V01L0Fi.002.h5"

# read the output file
from read_spec import SPEC
s=SPEC(filename)

# which toroidal cutplane to plot
toroidalIdx = 40

# extract slice corresponding to the given toroidal cutplane
pltR=s.poincare.R[:,:,toroidalIdx]
pltZ=s.poincare.Z[:,:,toroidalIdx]

# plot the Poincare data
import matplotlib.pyplot as plt
plt.figure()
plt.plot(pltR, pltZ, 'k.', ms=1)
plt.axis("equal")
plt.xlabel("R / m")
plt.ylabel("Z / m")

