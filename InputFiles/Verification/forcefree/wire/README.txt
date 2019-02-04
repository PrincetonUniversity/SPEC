This is a test case designed to test and verify that the virtual casing routine in SPEC calculates correctly the normal component of the plasma-induced magnetic field on the computational boundary. The plasma volume is small compared to that enclosed by the computational boundary, thereby approximating a current toroidal wire (axisymmetry is imposed). Attached a Poincare plot illustrating the equilibrium model. To reproduce it, run any of the input files with mfreeits=0. In order to make SPEC use the virtual casing routine to calculate Bn from the plasma without iterating on the geometry, mfreeits=1 and Lfindzero=0 are set in all the attached input files.

The results from SPEC should be compared to the semi-analytical evaluation performed with the attached matlab routines (written by Zhisong Qu, slightly modified by Joaquim Loizu, and explained below). Attached is a plot showing how the difference in Bns harmonics between SPEC and the wire-model decreases with some power law as the size of the plasma is decreased. This is achieved by reducing phiedge, Rbc(m=1,n=0), Zbs(m=1,n=0) with the proportion given by the fact that Rbc(m=1,n=0)~Zbs(m=1,n=0)~phiedge^2. Also the total current in the wire is kept constant by increasing mu, i.e. mu_0*Iwire = mu*phiedge = constant.

Using the matlab routines:

First, run

[ theta, bn_raw ] = bn( Rbc, Zbs, nfft, nint, R0, mu0I )
 
which computes the normal field. One needs to provide Rbc and Zbs, number of points you want to compute in poloidal direction, the number of points for the integral evaluating the field at (R,Z), the major radius of the wire R0, the total current mu0*I.
 
Then, run
 
[ bns ] = fft_field( bn_raw, mpol)
 
where mpol is the number of poloidal harmonics you want to keep.
