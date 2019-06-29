#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 25 18:01:55 2019

@author: Jonathan Schilling (jonathan.schilling@ipp.mpg.de)
"""

import h5py

# reader class for Stepped Pressure Equilibrium Code output file
# Hudson et al., Physics of Plasmas 19, 112502 (2012); doi: 10.1063/1.4765691
class SPEC:

    def __init__(self, *args, **kwargs):
        # args[0] should always be the name of a file or an item inside the root object
        # if args[0] is not a filename, kwagrs['content'] should be the content to be added
        #  as self.`args[0]`
        
        _content = None
        if kwargs.get('content') == None:
            # assume arg[0] is a filename
            _content = h5py.File(args[0], "r")
        elif isinstance(kwargs['content'], h5py.Group):
            _content = kwargs['content']
        
        if (_content != None):
            for key in _content:
                if isinstance(_content[key], h5py.Group):
#                    print("group '"+key+"'")
                    setattr(self, key, SPEC(content=_content[key]))
                elif isinstance(_content[key], h5py.Dataset):
#                    print("dataset '"+key+"' is '"+str(_content[key])+"'")
                    shape = _content[key].shape
                    nDim = len(shape)
                    if nDim==1:
                        if shape[0]==1:
                            setattr(self, key, _content[key][0])
                        else:
                            setattr(self, key, _content[key][:])
                    elif nDim==2:
                        setattr(self, key, _content[key][:,:])
                    else:
                        print("3 dim not supported: '"+key+"'")

        
if __name__=="__main__":
    
    filename = "/home/IPP-HGW/jons/04_PhD/00_programs/SPEC/InputFiles/TestCases/G3V02L1Fi.001.h5"

    s=SPEC(filename)
    print(s)
