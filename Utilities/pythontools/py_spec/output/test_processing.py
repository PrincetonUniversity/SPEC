from py_spec import SPECout
import numpy as np

d = SPECout('/BIG_14TB/BetaScan/RotatingEllipse/Nfp5/BSmodel1_Pmodel1_Bz002/BetaScan_30/BetaScan_Bz002_run_5.sp.h5')

j_dot_B = d.get_surface_current_density(lvol=np.arange(0,7), nt=32, nz=16)

print(j_dot_B)