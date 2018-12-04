This is a test case designed to test and verify that the virtual casing routine in SPEC calculates correctly the normal component of the plasma-induced magnetic field on the computational boundary. The plasma volume is small compared to that enclosed by the computational boundary, thereby approximating a current toroidal wire (axisymmetry is imposed).

The results from SPEC should be compared to the semi-analytical evaluation performed with the attached matlab routines, which were written by Zhisong Qu and explained below. Attached is a plot showing how the difference in Bns harmonics between SPEC and the wire-model decreases with some power law as the size of the plasma is decreased.

Using the matlab routines:

First, run

[ theta, bn_raw ] = bn( Rbc, Zbs, n, R0, mu0I )
 
which computes the normal field. One needs to provide Rbc and Zbs, number of points you want to compute in poloidal direction (typically 1024 or something), the major radius of the wire R0, the total current mu0*I (the unit your choice).
 
Then, run
 
[ bns ] = fft_field( bn_raw, mpol)
 
where mpol is the number of poloidal harmonics you want to keep.
