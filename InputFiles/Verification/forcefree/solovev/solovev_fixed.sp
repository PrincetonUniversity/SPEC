&physicslist
 Igeometry   =         3
 Istellsym   =         1
 Lfreebound  =         0
 phiedge     =   1.000000000000000E+00
 curtor      =   0.5537492180127451
 curpol      =   0.000000000000000E+00
 gamma       =   0.000000000000000E+00
 Nfp         =         1
 Nvol        =         4
 Mpol        =         5
 Ntor        =         0
 Lrad        =       6     6       6        6
 tflux       =    0.25    0.5    0.75      1.000000000000000E+00    
 pflux       =    0.00  0.15   0.15     0.15
 helicity    =   1.559429589793997E-03  1.559429589793997E-03  1.559429589793997E-03  1.559429589793997E-03
 pscale      =   0.0 ! 0.125000000000000E+00
 Ladiabatic  =         0
 pressure    =   0.875   0.625   0.375   0.125  0.0
 adiabatic   =   1.0     1.0     1.0     1.0    0.0
 mu          =   0.1     0.1     0.1     0.1    0.1
 Lconstraint =         1
 pl          =                       0          0            0             0         0
 ql          =                       0          0            0             0         0
 pr          =                       0          0            0             0         0
 qr          =                       0          0            0             0         0
 iota        =   1.00000000000000E+00  1.0  1.00   1.00  1.00 
 lp          =                       0          0            0             0         0
 lq          =                       0          0            0             0         0
 rp          =                       0          0            0             0         0
 rq          =                       0          0            0             0         0
 oita        =   1.000000000000000E+00  1.00 1.00   1.00   1.00 
 mupftol     =   1.000000000000000E-12
 mupfits     =       128
 Rac         =   3.999
 Zas         =   0.000000000000000E+00 
 Ras         =   0.000000000000000E+00  
 Zac         =   0.000000000000000E+00  
  RBC( 0,0) =  3.999     ZBS( 0,0) =  0.000
  RBC( 0,1) =  1.026     ZBS( 0,1) = -1.580
  RBC( 0,2) = -0.068     ZBS( 0,2) = -0.010
/
&numericlist
 Linitialize =         1
 Ndiscrete   =         2
 Nquad       =        -1
 iMpol       =        -4
 iNtor       =        -4
 Lsparse     =         0
 Lsvdiota    =         0
 imethod     =         3
 iorder      =         2
 iprecon     =         1
 iotatol     =  -1.000000000000000E+00
/
&locallist
 LBeltrami   =         4
 Linitgues   =         1
/
&globallist
 Lfindzero   =         2
 escale      =   0.000000000000000E+00
 pcondense   =   4.000000000000000E+00
 forcetol    =   1.000000000000000E-12
 c05xtol     =   1.000000000000000E-12
 c05factor   =   1.000000000000000E-04
 LreadGF     =         F
 opsilon     =   1.000000000000000E+00
 epsilon     =   1.000000000000000E+00
 upsilon     =   1.000000000000000E+00
/
&diagnosticslist
 odetol      =   1.000000000000000E-07
 absreq      =   1.000000000000000E-08
 relreq      =   1.000000000000000E-08
 absacc      =   1.000000000000000E-04
 epsr        =   1.000000000000000E-08
 nPpts       =        2500
 nPtrj       =        10
 LHevalues   =         F
 LHevectors  =         F
/
&screenlist
 Wpp00aa = T
/
