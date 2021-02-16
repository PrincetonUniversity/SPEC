"""
SPEC python wrapper
Author: Caoxiang Zhu (caoxiangzhu@gmail.com)
"""
from __future__ import print_function, absolute_import, division
import numpy as np
import sys
import os
import logging
from mpi4py import MPI
import spec

logger = logging.getLogger("[{}]".format(MPI.COMM_WORLD.Get_rank()) + __name__)


class SPEC(object):
    def __init__(self, input_file="", verbose=False, comm=MPI.COMM_WORLD, **kwargs):
        """Initialization of SPEC runs

        Args:
            input_file (str) : Filename for SPEC input namelist. (default: '').
            verbose (bool): If wants scree outputs. (default: True).
            comm (MPI communicator): the comunicator assigned to SPEC. Default: MPI.COMM_WORLD
        Returns:
            None
        """
        # pass arguments and check
        assert isinstance(
            input_file, str
        ), "input_file should the input filename in str."
        if not input_file.endswith(
            ".sp"
        ):  # This causes problems if the filename starts with a drectory!
            input_file = input_file + ".sp"
        self.input_file = input_file
        self.comm = comm.py2f()
        assert isinstance(verbose, bool), "verbose is either True or False."
        self.verbose = verbose
        # Fortran libaries accessed via self.lib
        self.lib = spec

        # wrap around modules
        modules = [
            "constants",
            "numerical",
            "fileunits",
            "cputiming",
            "typedefns",
            "inputlist",
            "allglobal",
        ]
        for key in modules:
            setattr(self, key, getattr(spec, key))

        # assign ext and comm
        self.inputlist.ext = input_file[:-3]
        self.allglobal.mpi_comm_spec = self.comm
        self.initialized = False
        return

    def run(self):
        if not self.initialized:
            self.readin()
        self.lib.spec()
        return

    def readin(self):
        self.lib.allglobal.readin()
        return


if __name__ == "__main__":
    ext = sys.argv[1]
    if ".sp" in ext:
        ind = ext.index("sp")
        ext = ext[: ind - 1]
    print("Begin to run SPEC from python with input file at ", ext + ".sp")
    comm = MPI.COMM_WORLD
    test = SPEC(input_file=ext, comm=comm, verbose=True)
    test.run()
    print("SPEC called from python finished!")
