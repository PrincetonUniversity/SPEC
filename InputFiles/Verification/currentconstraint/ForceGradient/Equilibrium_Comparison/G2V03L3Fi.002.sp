&physicslist
 Igeometry   =         2
 Istellsym   =         1
 Lfreebound  =         0
 phiedge     =   1.000000000000000E+00
 curtor      =   0.000000000000000E+00
 curpol      =   0.000000000000000E+00
 gamma       =   0.000000000000000E+00
 Nfp         =         1
 Nvol        =         3
 Mpol        =         0
 Ntor        =         0 
 Lrad        =         4 4 4
 tflux       =   0.5 1.0 1.7 
 pflux       =   0.0 0.1 0.2              
 helicity    =  -1.000000000000000E-01 
 pscale      =   0.000000000000000E+00
 Ladiabatic  =         0
 pressure    =   0.02  1.000000000000000E-02  0.0 
 adiabatic   =   0.000000000000000E+00 
 mu          =        0.8 1.2 1.4
 Ivolume     =   0.235294117647059  0.588235294117647  1.164705882352941
 Isurf       =   -0.106577555488019 -0.409209884152368 0.000000000000000
 Lconstraint =         3
 pl          =                       0                      0                      0
 ql          =                       0                      0                      0
 pr          =                       0                      0                      0
 qr          =                       0                      0                      0
 iota        =                       0                      0                      0
 lp          =                       0                      0                      0
 lq          =                       0                      0                      0
 rp          =                       0                      0                      0
 rq          =                       0                      0                      0
 oita        =                       0                      0                      0
 mupftol     =   1.000000000000000E-14
 mupfits     =         128
 Rac         =   0.000000000000000E+00
 Zas         =   0.000000000000000E+00
 Ras         =   0.000000000000000E+00
 Zac         =   0.000000000000000E+00
Rbc(0,0)    =  1.00000000000000E+00 Zbs(0,0)    =  0.000000000000000E+00 Rbs(0,0)    =  0.000000000000000E+00 Zbc(0,0)    =  0.000000000000000E+00
Rwc(0,0)    =  0.000000000000000E+00 Zws(0,0)    =  0.000000000000000E+00 Rws(0,0)    =  0.000000000000000E+00 Zwc(0,0)    =  0.000000000000000E+00
Vns(0,0)    =  0.000000000000000E+00 Bns(0,0)    =  0.000000000000000E+00 Vnc(0,0)    =  0.000000000000000E+00 Bnc(0,0)    =  0.000000000000000E+00
/
&numericlist
 Linitialize =         1
 Lzerovac    =         0
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
 Lextrap     =         0
 Mregular    =        -1
/
&locallist
 LBeltrami   =         4
 Linitgues   =         1
 Lposdef     =         0
/
&globallist
 Lfindzero   =         2
 escale      =   0.000000000000000E+00
 opsilon     =   1.000000000000000E+00
 pcondense   =   2.000000000000000E+00
 epsilon     =   0.000000000000000E+00
 wpoloidal   =   1.000000000000000E+00
 upsilon     =   1.000000000000000E+00
 forcetol    =   1.000000000000000E-14
 c05xmax     =   1.000000000000000E-06
 c05xtol     =   1.000000000000000E-14
 c05factor   =   1.000000000000000E-02
 LreadGF     =         F
 mfreeits    =         0
 gBntol      =   1.000000000000000E-06
 gBnbld      =   6.660000000000000E-01
 vcasingeps  =   1.000000000000000E-12
 vcasingtol  =   1.000000000000000E-08
 vcasingits  =         8
 vcasingper  =         1
/
&diagnosticslist
 odetol      =   1.000000000000000E-07
 nPpts       =         0
 nPtrj       =         10    10 10
 LHevalues   =         F
 LHevectors  =         F
 LHmatrix    =         T
 Lperturbed  =         0
 dpp         =        -1
 dqq         =        -1
 Lcheck      =       0
 dRZ         =        1E-5
 Ltiming     =         F
/
&screenlist
Wjo00aa      = F
Wpp00aa      = F
Wlbpol       = F
Wlforce      = F
/
     0     0  5.850250756131968E-01  0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00  1.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00
