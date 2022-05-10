#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Caoxiang Zhu (czhu@ppp.gov)
For any help, type ./test.py -h
"""
import numpy as np
from py_spec.output import SPECout
import argparse

parser = argparse.ArgumentParser(description="Compare two SPEC HDF5 outputs")
parser.add_argument("filename", type=str, help="file name to be compared")
parser.add_argument("reference", type=str, help="reference data")
parser.add_argument("-t", "--tol", type=float, default=1E-12, help="difference tolerance")


def compare(data, reference, localtol=1e-6, action='ERR'):
    """
    compare all items in data to items in reference with the same keys.
    Throws an error or prints a warning when there is a diference.
    *action*: 'ERR' or 'WARN'. If ERR an error is thrown and the program
    fails on exit. if 'WARN' only an error is printed.
    """
    global match
    for key, value in vars(data).items():
        if isinstance(value, SPECout):  # recurse data (csmiet: I'm not the biggest fan of this recursion...)
            print('------------------')
            print('Elements in '+key)
            if key in ['poincare']:
                print('differences in ' + key + ' are not important to regression')
                compare(value, reference.__dict__[key], localtol, action='WARN')
            else:
                compare(value, reference.__dict__[key], localtol)
        else:
            if key in ['filename', 'version', 'iterations']:  # not compare filename and version and iterations
                continue
            else:
                if key in ['volume', 'fiota', 'Lmatsolver']:  # skip certain problematic variables. NOT A GOOD IDEA TO CHANGE (might be revised)
                    action = 'WARN'

                diff = 0.0
                if isinstance(value, list):  # compare each list
                    for ii, item in enumerate(value):
                        diff = np.max([diff, np.max(np.abs(np.array(item) - np.array(reference.__dict__[key][ii])))])
                else:
                    if np.shape(value) == np.shape(reference.__dict__[key]):
                        diff = np.max(np.abs(np.array(value) - np.array(reference.__dict__[key])))
                    else:
                        match = False
                        print('ERROR: '+key, ', dimensions mismatch: {} .ne. {}'.format(np.shape(value), np.shape(reference.__dict__[key])))
                        next # there is no point in computing differences if not even the dimensions match, so skip to next variable here
                unmatch = diff > localtol
                if unmatch:
                    if action == 'ERR':
                        match = False
                        print('ERROR: '+key, ', element average difference = {:12.5E}'.format(diff))
                    if action == 'WARN':
                        print('WARNING: '+key, ', element average difference = {:12.5E}'.format(diff))
                else:
                    print('ok: ',key, 'element average difference = {}'.format(diff))
    return

if __name__ == '__main__':
# parse command line arguments

    args = parser.parse_args()
    print('Compare SPEC outputs in {:s} and {:s} with tolerance {:12.5E}'.format(
            args.filename, args.reference, args.tol))
    data_A = SPECout(args.filename)
    data_B = SPECout(args.reference)
    tol = args.tol
    match = True

    compare(data_A, data_B, tol)
    print('===================')
    if match:
        print('All the terms are within tolerance.')
    else:
        print('Differences in some elements are larger than the tolerence.')
        raise SystemExit('Differences in some elements are larger than the tolerence.')

    exit

