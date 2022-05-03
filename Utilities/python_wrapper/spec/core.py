"""
SPEC python wrapper
Author: Caoxiang Zhu (caoxiangzhu@gmail.com)
        Jonatahn Schilling (jonathan.schilling@mail.de)
"""
from __future__ import print_function, absolute_import, division
import numpy as np
import sys
import os
import logging
from mpi4py import MPI
import spec.spec_f90wrapped as spec_lib


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
        self.lib = spec_lib

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
            setattr(self, key, getattr(self.lib, key))

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

    # def run(self, save=True):
    def run(self,
            save_output: bool = False):
        """
        Args:
            save_output: Whether or not to save the hdf5 and restart files.
        """
        if not self.initialized:
            self.read()

        if save_output:
            self.allglobal.skip_write = False
            self.lib.sphdf5.init_outfile()
            self.lib.sphdf5.mirror_input_to_outfile()
            if self.comm.rank == 0:
                self.lib.allglobal.wrtend()
            self.lib.sphdf5.init_convergence_output()
        else:
            self.allglobal.skip_write = True

        self.lib.spec()

        if save_output:
            self.lib.final_diagnostics()
            self.lib.sphdf5.write_grid()
            if self.comm.rank == 0:
                self.lib.allglobal.wrtend()
            self.lib.sphdf5.hdfint()
            self.lib.sphdf5.finish_outfile()
            self.lib.ending()

        return

    def read(self, input_file=None):
        # initialize input quantities to a known state
        if self.comm.rank == 0:
            self.inputlist.initialize_inputs()

            if input_file is not None:
                print("Read SPEC input namelist from {:}.sp".format(input_file))
                self.allglobal.ext = input_file

            self.allglobal.read_inputlists_from_file()
            self.allglobal.check_inputs()

        self.allglobal.broadcast_inputs()
        self.lib.preset()

        self.initialized = True
        return


if __name__ == "__main__":

    if len(sys.argv) < 2:
        # no command line arguments given
        print("usage: "+sys.argv[0]+" input_to_SPEC.sp")
        exit()

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
