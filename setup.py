#!/usr/bin/env python3

import os # environ
import sys
import json
import argparse
import setuptools

from os.path import basename, splitext
from glob import glob
import numpy

print("system.platform is {}".format(sys.platform))
if (sys.platform == "darwin"):
    from distutils import sysconfig
    vars = sysconfig.get_config_vars()
    vars['LDSHARED'] = vars['LDSHARED'].replace('-bundle', '-dynamiclib')

from setuptools import find_packages
from skbuild import setup

# Load machine-specific options from cmake_config.json.
# It should contain a dict with an array called cmake_args.
with open('cmake_config.json') as fp:
    d = json.load(fp)

# Include additional parameters from CMAKE_ARGS environment variable.
# This is the way Anaconda tells CMake its specific needs.
if 'CMAKE_ARGS' in os.environ:

  print("CMAKE_ARGS = '%s'"%(os.environ['CMAKE_ARGS']))
  for cmake_arg in os.environ['CMAKE_ARGS'].split(" "):
    d['cmake_args'].append(cmake_arg)

# Tell CMake about where to find numpy libraries
# see also: https://stackoverflow.com/a/14657667
d['cmake_args'].append("-DCMAKE_C_FLAGS=-I"+numpy.get_include())

setup(
    name="spec",
    version="0.0.3",
    #license="MIT",
    packages=['spec'],
    package_dir={'': 'Utilities/python_wrapper'},
    #py_modules=[splitext(basename(path))[0] for path in glob('src/vmec/*.py')],
    install_requires=['f90wrap', 'scikit-build'],
    classifiers=[
        "Development Status :: 1 - Alpha",
        "Intended Audience :: Nuclear Fusion Community",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Topic :: MHD Equilibrium Solver"],
    cmake_args=d['cmake_args'],
)
