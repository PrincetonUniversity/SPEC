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
        self.comm = comm
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

        # assign ext and set MPI communicator
        self.allglobal.ext = input_file[:-3] # omit ".sp" at end

         # py2f converts the Python object to the Fortran integer identifying an MPI communicator.
        self.allglobal.set_mpi_comm(self.comm.py2f())

        self.initialized = False

        # mute screen output if necessary
        # TODO: relies on /dev/null being accessible (Windows!)
        if not self.verbose:
            self.fileunits.mute(1)
        return

    def run(self, save=True):
        if not self.initialized:
            self.read()
        self.lib.spec()
        if save:
            self.write()
        return

    def read(self, input_file=None):
        # initialize input quantities to a known state
        self.inputlist.initialize_inputs()

        if input_file is not None:
            print("Read SPEC input namelist from {:}.sp".format(input_file))
            self.allglobal.ext = input_file
        self.allglobal.readin()
        self.initialized = True
        return

    def write(self, output_file=None):
        if output_file is not None:
            print("Write SPEC output into {:d}.sp.h5".format(output_file))
            ext = self.inputlist.ext  # save original ext value
            self.inputlist.ext = output_file
            self.lib.write_hdf5()
            self.inputlist.ext = ext  # reset the ext value
        else:
            self.lib.write_hdf5()
        return


if __name__ == "__main__":
    ext = sys.argv[1]
    if ".sp" in ext:
        ind = ext.index("sp")
        ext = ext[: ind - 1]
    comm = MPI.COMM_WORLD
    rank = comm.rank
    if rank == 0:
        print("Begin to run SPEC from python with input file at ", ext + ".sp")
    test = SPEC(input_file=ext, comm=comm, verbose=True)
    test.run()
    if rank == 0:
        print("SPEC called from python finished!")
