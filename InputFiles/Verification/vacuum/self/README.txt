An exact verification of the code can be made by calculating the error in the Beltrami field (curlB - muB) component by component and in each relaxed volume, and then showing that the error converges towards machine precision as resolution is increased. Usually, that happens as Mpol=Ntor is increased, so as long as Lrad is sufficiently large.

The input file here is the start of a sequence of runs for a 5-field-period rotating ellipse field. One can produce a scan like the one shown in MNscan.jpg, where Mpol=Ntor was scanned from 4 to 14 (Lrad=6 except for Mpol=Ntor=12,14 where Lrad=8 in order to keep the error descending at the same rate).

Runs were carried out with 2 processors and took a range of times: 2 minutes for M=N=4, 70 hours for M=N=12, 5 days for M=N=14.

The error in the beltrami field is written in the output HDF5 file and easily read, e.g. with the provided matlab functions.
