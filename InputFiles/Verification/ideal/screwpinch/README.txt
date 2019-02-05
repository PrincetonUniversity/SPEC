
SUMMARY
-------

Verification of SPEC in cylindrical geometry when the number of volumes grows to infinity. 

SPEC solution should converge to the ideal MHD solution of the screw pinch when the number of volumes increases. In mhd_limit_study.m, the number of volumes for multiple SPEC simulations can be set. The MHD solution is obtained with ScrewPinch_IotaSolver.m and used to prepare SPEC input file. Finally, SPEC input file is written in run_from_mhd.m and SPEC is run.

SPEC solution is read and transformed in the canonical components of B by get_full_field.m. These data can then be used to compute the rotational transform and the poloidal and toroidal flux. 

Six plots are generated. The pressure, the rotational transform, the poloidal and toroidal component of the magnetic field as well as the poloidal and toroidal flux are plotted.

-------------------------------------------------------------------------------------------------------

NOTES
-----

* Run in the same folder as SPEC execucutable (xspec).
* Don't forget to add the path to SPEC/Utilities/matlabtools/


A.Baillod (2019)

