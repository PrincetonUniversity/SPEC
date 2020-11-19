#!/usr/bin/env python3
import sys
import json
import argparse
import setuptools

from os.path import basename, splitext
from glob import glob

print("system.platform is {}".format(sys.platform))
if (sys.platform == "darwin"):
    from distutils import sysconfig
    vars = sysconfig.get_config_vars()
    vars['LDSHARED'] = vars['LDSHARED'].replace('-bundle', '-dynamiclib')

from setuptools import find_packages
from skbuild import setup

with open('cmake_config_file.json') as fp:
    d = json.load(fp)

setup(
    name="spec",
    version="0.0.1",
    #license="MIT",
    packages=['spec'],
    #package_dir={'': 'src'},
    #py_modules=[splitext(basename(path))[0] for path in glob('src/vmec/*.py')],
    install_requires=['f90wrap == v0.2.3'],
    classifiers=[
        "Development Status :: 1 - Alpha",
        "Intended Audience :: Nuclear Fusion Community",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Topic :: MHD Equilibrium Solver"],
    cmake_args=d['cmake_args'],
    cmake_source_dir="../.."
)
