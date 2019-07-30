# SPEC-VMEC free boundary calculation benchmark

The attached files constitute a free-boundary verification between SPEC and VMEC.

---
## Comments from Stuart:

1) The SPEC input files are labelled with *.sp, and the output files are *.sp.end and *.sp.h5.

2) Caoxiang prepared the coils. The input for FOCUS is attached.

3) Caoxiang prepared the input for VMEC, namely the mgrid file. A high-resolution VMEC input and output is included.

4) The postscript file nvlimit.ps shows the convergence of two measures of the error. One error measures the root mean square difference between all the SPEC interfaces and the corresponding VMEC surfaces. The other error measure shows the difference between the i = Nv/2 SPEC interface and the corresponding VMEC interface. There is a small difference at the highest SPEC resolution, with Nvol = 128.

5) There is a Poincare plot of the Nv = 8 SPEC calculation that shows a comparison with VMEC. The normal field on the computational boundary is shown, both the normal field produced by the coils and the normal field produced by the plasma.

6) Caoxiang will provide a plot of the coils.

7) The SPEC calculation with Nvol = 256 did not succeed. I expect that this is because the ill-conditioning problem near the origin. The new version of SPEC prepared by Zhisong, which uses Zernike polynomials in the innermost volume, might work better. When I have this executable, I will try to get the Nv = 256 SPEC calculation to work.


In summary, the convergence looks good. I will begin work on the publication.

---

## Comments from Caoxiang

1. This branch includes the implementation of auto-initialization of Bns, which Stuart was not using when he ran free-boundary SPEC. It is optional because Bns is already close to the actual values. 

2. The coils consists of 12 PF coils and one central current (3.884526409876309MA). `modified_1r.focus` contains all the data in Fourier harmnics, while `coils.modified_1r` has the xyz points for all the PF coils.

3. The mgrid file 'mgrid.focus_modified_1r_vacuum_1601' is in binary format and the resolutions are `Np=1, Nr=1601, Nz=1601`.
