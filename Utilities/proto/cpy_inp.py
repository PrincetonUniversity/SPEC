#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 14:41:05 2020

@author: jonathan
"""

import h5py

f_src = h5py.File("demo_hdf5.sp.h5", "r")
f_tgt = h5py.File("demo_hdf5_2.sp.h5", "a")

# copy group '/input' from f_src into f_tgt
f_src.copy('/input', f_tgt)

f_tgt.close()
f_src.close()