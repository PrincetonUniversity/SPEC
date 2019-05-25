#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 25 18:01:55 2019

@author: Jonathan Schilling (jonathan.schilling@ipp.mpg.de)
"""

import h5py

# reader class for Stepped Pressure Equilibrium Code
# Hudson et al., Physics of Plasmas 19, 112502 (2012); doi: 10.1063/1.4765691
class SPEC:
    
    _file=None
    Igeometry=0
    
    def __init__(self, filename):
        self._file=h5py.File(filename, "r")
        
        # input/physics: /physicslist/ from input file
        self.Igeometry  =self._file['/input/physics/Igeometry'][0]
        self.Istellsym  =self._file['/input/physics/Istellsym'][0]
        self.Lfreebound =self._file['/input/physics/Lfreebound'][0]
        self.phiedge    =self._file['/input/physics/phiedge'][0]
        self.curtor     =self._file['/input/physics/curtor'][0]
        self.curpol     =self._file['/input/physics/curpol'][0]
        self.gamma      =self._file['/input/physics/gamma'][0]
        self.Nfp        =self._file['/input/physics/Nfp'][0]
        self.Nvol       =self._file['/input/physics/Nvol'][0]
        self.Mpol       =self._file['/input/physics/Mpol'][0]
        self.Ntor       =self._file['/input/physics/Ntor'][0]
        self.Lrad       =self._file['/input/physics/Lrad'][:]
        self.Lconstraint=self._file['/input/physics/Lconstraint'][0]
        self.tflux      =self._file['/input/physics/tflux'][:]
        self.pflux      =self._file['/input/physics/pflux'][:]
        self.helicity   =self._file['/input/physics/helicity'][:]
        self.pscale     =self._file['/input/physics/pscale'][0]
        self.pressure   =self._file['/input/physics/pressure'][:]
        self.Ladiabatic =self._file['/input/physics/Ladiabatic'][0]
        self.adiabatic  =self._file['/input/physics/adiabatic'][:]
        self.mu         =self._file['/input/physics/mu'][:]
        self.pl         =self._file['/input/physics/pl'][:]
        self.ql         =self._file['/input/physics/ql'][:]
        self.pr         =self._file['/input/physics/pr'][:]
        self.qr         =self._file['/input/physics/qr'][:]
        self.iota       =self._file['/input/physics/iota'][:]
        self.lp         =self._file['/input/physics/lp'][:]
        self.lq         =self._file['/input/physics/lq'][:]
        self.rp         =self._file['/input/physics/rp'][:]
        self.rq         =self._file['/input/physics/rq'][:]
        self.oita       =self._file['/input/physics/oita'][:]

        self.Rac        =self._file['/input/physics/Rac'][:] 
        self.Zas        =self._file['/input/physics/Zas'][:]
        self.Ras        =self._file['/input/physics/Ras'][:] 
        self.Zac        =self._file['/input/physics/Zac'][:]
        
        self.Rbc        =self._file['/input/physics/Rbc'][:,:]
        self.Zbs        =self._file['/input/physics/Zbs'][:,:]
        self.Rbs        =self._file['/input/physics/Rbs'][:,:]
        self.Zbc        =self._file['/input/physics/Zbc'][:,:]
        
        self.Rwc        =self._file['/input/physics/Rwc'][:,:]
        self.Zws        =self._file['/input/physics/Zws'][:,:]
        self.Rws        =self._file['/input/physics/Rws'][:,:]
        self.Zwc        =self._file['/input/physics/Zwc'][:,:]
        
        self.Vns        =self._file['/input/physics/Vns'][:,:]
        self.Bns        =self._file['/input/physics/Bns'][:,:]
        self.Vnc        =self._file['/input/physics/Vnc'][:,:]
        self.Bnc        =self._file['/input/physics/Bnc'][:,:]
        
        self.mupftol    =self._file['/input/physics/mupftol'][0]
        self.mupfits    =self._file['/input/physics/mupfits'][0]
        
        
if __name__=="__main__":
    s=SPEC("/home/jonathan/Uni/04_PhD/00_programs/SPEC/SPEC/InputFiles/Verification/forcefree/solovev/solovev_fb_vmec7vol_final.h5")