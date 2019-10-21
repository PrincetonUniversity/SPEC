The Hessian at fixed helicity can be verified against analytical predictions to machine precision. This is done in slab geometry by considering a 3-volume equilibrium (one stable current sheet, one unstable current sheet) at given radial resolution Lrad (given mpol=1,ntor=0) and comparing the SPEC output hessian with the exact analytical predictions given in [Loizu and Hudson, PoP 26, 030702 (2019)].

Attached input files for SPEC (G1V03L2Fi.001.sp is a stable current sheet, G1V03L2Fi.002.sp is an unstable current sheet). Varying Lrad from 2 to 16, one can show convergence of the error in the Hessian (maximum absolut difference between the coefficients) towards machine precision (attached figure showing such convergence). Attached also comparison of the entire matrices for the stable and unstable current sheet cases, for the SPEC calculation at Lrad=12 and for the analytical prediction.

Attached also matlab routines that can be used to generate the analytical Hessian. One has to first set, for the stable case,
mu2=0.1
tf2=0.2
pf3=0.04
R3=10
and run hessian.m, thereby obtaining Hmatrix from the theoretical prediction. For the unstable case, one has to set
mu2=0.1
tf2=0.1
pf3=0.0225
R3=10

The Hessian from SPEC can easily be obtained by using the matlab utilities,
hdata = read_spec_hessian(fname);
Hspec = hdata.Hmatrix;
and can be plotted using plot_spec_hessian(hdata).
