#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Caoxiang Zhu (czhu@ppp.gov)
For any help, type ./compare_spec.py -h
"""
import numpy as np
from py_spec import SPEC


def compare_spec(data, reference, tol=1e-12):
    match = True
    for key, value in vars(data).items():
        if isinstance(value, SPEC):  # recurse data
            print('------------------')
            print('Elements in '+key)
            compare_spec(value, reference.__dict__[key])
        else:
            if key in ['filename', 'version']:  # not compare filename and version
                continue
            elif key == 'iterations': # skip iteration data (might be revised)
                continue
            else:
                # print(key)
                diff = np.linalg.norm(np.abs(np.array(value) - np.array(reference.__dict__[key])))
                unmatch = diff > tol
                if unmatch:
                    match = False
                    print('UNMATCHED: '+key, ', diff={:12.5E}'.format(diff))
                else :
                    print('ok: '+key)
    return match