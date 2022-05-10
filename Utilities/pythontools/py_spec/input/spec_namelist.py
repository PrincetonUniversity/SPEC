"""
SPECNamelist.py: contains the object that reads and writes the SPEC namelist
SPECNamelist extends the f90nml.Namelist object
Requires:
    numpy
    f90nml, install it using "pip install f90nml"
@author: Zhisong Qu (zhisong.qu@anu.edu.au)
"""
import f90nml
from f90nml import Namelist

import numpy as np


class SPECNamelist(Namelist):
    """The SPEC namelist class
    To get a content within the namelist, use:
        somevariable = spec_nml['whichlist']['whichitem']

    To change (or add) an item on the namelist, use:
        spec_nml['whichlist']['whichitem'] = somevalue
    e.g. spec_nml['physicslist']['Lconstraint'] = 1 sets Lconstraint in physicslist list to 1

    Please do not change Mpol, Ntor and Nvol directly. To change them, please use
        spec_nml.update_resolution(new_Mpol, new_Ntor)
        spec_nml.insert_volume(ivol, tflux)
        spec_nml.remove_volume(ivol)

    The guess of the interface Fourier harmonics, if exists, can be obtained by
        spec_nml.get_interface_guess(m_number, n_number, ivol, 'Rbc')
    and changed by
        spec_nml.set_interface_guess(m_number, n_number, ivol, 'Rbc')
        spec_nml.remove_interface_guess(m_number, n_number)
    Alternatively, one can directly modify spec_nml.interface_guess

    member functions:
        write, run, update_resolution, insert_volume, remove_volume,
        get_interface_guess, set_interface_guess, remove_interface_guess
    """

    # items in 'physicslist' that are arrays with index 0 being the quatities in the first volume
    zero_index_keys = [
        "tflux",
        "pflux",
        "mu",
        "helicity",
        "pressure",
        "adiabatic",
        "Lrad",
        "nPtrj",
        "Ivolume",
        "Isurf",
    ]

    # items in 'physicslist' that are arrays with index 1 being the quatities in the first volume
    one_index_keys = ["iota", "oita", "qr", "pr", "ql", "pl", "rq", "rp", "lq", "lp"]

    # items for the boundary
    boundary_keys = [
        "Rbc",
        "Rbs",
        "Zbc",
        "Zbs",
        "Vns",
        "Vnc",
        "Bns",
        "Bnc",
        "Rwc",
        "Rws",
        "Zwc",
        "Zws",
    ]

    # items for the axis
    axis_keys = ["Rac", "Zas", "Ras", "Zac"]

    def __init__(self, *args, **kwargs):
        """Initialization
        use one of the following
        1) spec_nml = SPECNamelist(filename), e.g. spec_nml=SPECNamelist("namelist.sp") or spec_nml=SPECNamelist("/path/to/namelist.sp")
        2) spec_nml = SPECNamelist(spec_output_object), where spec_output_object
        """

        from py_spec.output import SPECout

        if isinstance(args[0], str):
            # the first argument is a string, read from this namelist
            # first create the namelist using __init__ of the namelist object

            with open(args[0], "r") as file_object:

                super().__init__(f90nml.read(file_object))

                # now we should mind some variables that are important: we include them in the object itself and need to monitor them
                self._Mpol = self["physicslist"]["mpol"]
                self._Ntor = self["physicslist"]["ntor"]
                self._Nvol = self["physicslist"]["nvol"]

                # then read the part that specifies the guess of the interface geometry
                self._read_interface_guess(file_object)

            # we don't know the verson of SPEC from its namelist, so leave it as 'unknown
            self._spec_version = "unknown"

        elif isinstance(args[0], SPECout):
            # the first argumetn is a SPEC object, generate the quantities from the SPEC object
            # first initialize an empty namelist
            super().__init__()

            # then generate the namelist from the SPEC object
            self._generate_namelist_from_output(args[0])

            version_str = "{:5.2f}"
            self._spec_version = version_str.format(args[0].version)

        else:
            # the input is of unknown type, abort
            raise ValueError(
                "The first argument should contain either a path to the SPEC namelist, or a SPEC object generated from reading SPEC hdf5 output"
            )

        # correct the size/type of the namelist objects and so on
        self._rectify_namelist()

    def read(self, *args, **kwargs):
        raise NotImplementedError(
            "SPECNamelist does not support 'read' directly. Please call the constructor __init__."
        )

    def write(self, filename, force=False):
        """Write to a namelist file
        parameters:
            filename -- the filename to write
            force -- force to overwrite or not
        """

        import os

        if not isinstance(filename, str):
            raise ValueError(
                "The first argument should contain the path of the namelist to write to"
            )

        # check if some important quantities has changed. If yes, prompt an error
        if not self._Mpol == self["physicslist"]["mpol"]:
            raise RuntimeError(
                "Inconsistent Mpol. If one wishes to change Mpol or Ntor, please call member function update_resolution"
            )
        if not self._Ntor == self["physicslist"]["ntor"]:
            raise RuntimeError(
                "Inconsistent Ntor. If one wishes to change Mpol or Ntor, please call member function update_resolution"
            )
        if not self._Nvol == self["physicslist"]["nvol"]:
            raise RuntimeError(
                "Inconsistent Nvol. If one wishes to change Nvol, please call member function insert_volume, remove_volume"
            )

        if os.path.exists(filename):
            if not force:
                raise Exception(
                    filename
                    + " already exists! Please set force to True if you would like to force overwrite."
                )

        with open(filename, "w") as file_object:
            from datetime import datetime

            intro_str = "! auto generated by SPECNamelist.py "
            file_object.write(intro_str)

            # write the time for generating the namelist
            file_object.write(datetime.now().strftime("%Y-%m-%d %H:%M:%S "))

            # write the version of SPEC if known
            if not self._spec_version == "unknown":
                file_object.write("SPEC version " + self._spec_version)

            # conclude the first line
            file_object.write("\n")

            # convert all np.ndarray to list
            for key1 in self:
                for key2 in self[key1]:
                    if isinstance(self[key1][key2], np.ndarray):
                        self[key1][key2] = self[key1][key2].tolist()

            # write the main content of the namelist
            super().write(file_object)

            # write the interface guess
            self._write_interface_guess(file_object)

    def run(self, spec_command="./xspec", filename="spec.sp", force=False, quiet=False):
        """Run SPEC on the current namelist and obtain its output
        parameters:
            spec_command -- the command to call SPEC, usually it looks like '/path/to/spec/xspec'
                            or 'mpirun -np (ncpus) /path/to/spec/xspec', with (ncpus) replaced by the number of cpus
            filename -- write this namelist to the temporary file 'filename' in current folder
            force -- if file exists, force overwrite or not
            quiet -- True if you want some more screen output

        Returns:
            result -- after running SPEC, read the output hdf5 with py_spec and return the SPEC object
        """
        import os
        import subprocess
        from py_spec.output import SPECout

        self.write(filename, force=force)

        if not quiet:
            print("SPEC is running...")

        run_result = subprocess.run(spec_command + " " + filename, shell=True)

        if run_result.returncode == 0:  # the run is successful
            if not quiet:
                print("SPEC runs successfully.")
            return SPECout(filename + ".h5")
        else:
            print("SPEC runs unsuccessfully, check terminal output.")
            return None

    def insert_volume(self, ivol=0, tflux=None):
        """Insert a volume
        parameters:
            ivol - insert volume inside the ivol-th volume (Python index starts from 0)
                    if ivol == Nvol, insert a volume outside the plasma boundary
            tflux - the tflux of the volume inserted, default is the middle point of the volume
        """
        assert ivol >= 0 and ivol <= self._Nvol
        Nvol = self._Nvol

        # set default tflux
        if tflux is None:
            if ivol == 0:
                tflux = self["physicslist"]["tflux"][0] / 2.0
            elif ivol == Nvol and self["physicslist"]["Lfreebound"] == 0:
                tflux = self["physicslist"]["tflux"][ivol - 1] * 1.1
            else:
                tflux = (
                    self["physicslist"]["tflux"][ivol - 1]
                    + self["physicslist"]["tflux"][ivol]
                ) / 2.0

        # we interpolate the initial guess first
        self._interpolate_guess(ivol, tflux)

        for key in self.zero_index_keys:
            if key.lower() in self["physicslist"].keys():
                if key.lower == "tflux":
                    continue
                # make it a list if there is only a single item, and convert it to float
                if isinstance(self["physicslist"][key], int) or isinstance(
                    self["physicslist"][key], float
                ):
                    self["physicslist"][key] = [float(self["physicslist"][key])]
                self["physicslist"][key].insert(
                    ivol, self["physicslist"][key][min(ivol, Nvol - 1)]
                )

            if key.lower() in self["diagnosticslist"].keys():
                self["diagnosticslist"][key].insert(
                    ivol, self["diagnosticslist"][key][min(ivol, Nvol - 1)]
                )

        for key in self.one_index_keys:
            if key.lower() in self["physicslist"].keys():
                self["physicslist"][key].insert(
                    ivol + 1, self["physicslist"][key][min(ivol + 1, Nvol - 1)]
                )

        self["physicslist"]["tflux"].insert(ivol, tflux)

        # update Nvol
        self._Nvol = self._Nvol + 1
        self["physicslist"]["nvol"] = self["physicslist"]["nvol"] + 1

    def remove_volume(self, ivol=0):
        """Remove a volume
        parameters:
            ivol - Remove the ivol-th interface (Python index starts from 0)
        """
        assert ivol >= 0 and ivol < self._Nvol

        if self._Nvol == 1:
            raise RuntimeError("At least one volume should remain")
        if ivol == self._Nvol - 1 and self["physicslist"]["Lfreebound"] == 0:
            raise RuntimeError(
                "You cannot remove the last interface (plasma boundary)."
            )

        zero_index_keys = [
            "tflux",
            "pflux",
            "mu",
            "helicity",
            "pressure",
            "adiabatic",
            "Lrad",
            "nPtrj",
            "Ivolume",
            "Isurf",
        ]
        for key in zero_index_keys:
            if key.lower() in self["physicslist"].keys():
                self["physicslist"][key].remove(self["physicslist"][key][ivol])

            if key.lower() in self["diagnosticslist"].keys():
                self["diagnosticslist"][key].remove(self["diagnosticslist"][key][ivol])

        one_index_keys = [
            "iota",
            "oita",
            "qr",
            "pr",
            "ql",
            "pl",
            "rq",
            "rp",
            "lq",
            "lp",
        ]
        for key in one_index_keys:
            if key.lower() in self["physicslist"].keys():
                self["physicslist"][key].remove(self["physicslist"][key][ivol + 1])

        guess_index_keys = ["Rbc", "Zbs", "Rbs", "Zbc"]
        for harmonics, data in self.interface_guess.items():
            for key in guess_index_keys:
                data[key] = np.delete(data[key], ivol)

        # update Nvol
        self._Nvol = self._Nvol - 1
        self["physicslist"]["nvol"] = self["physicslist"]["nvol"] - 1

    def update_resolution(self, new_Mpol, new_Ntor):
        """Change the Fourier resolution of the SPEC namelist
        parameters:
            new_Mpol, new_Ntor -- the new Mpol and Ntor
        """
        if new_Mpol == self._Mpol and new_Ntor == self._Ntor:
            # nothing needs to be done
            return None
        if new_Mpol < 0 or new_Ntor < 0:
            raise ValueError("Mpol or Ntor >= 0")

        # We will need to update the size of Rbc, etc and their indexing
        changelist = [
            "Rbc",
            "Rbs",
            "Zbc",
            "Zbs",
            "Vns",
            "Vnc",
            "Bns",
            "Bnc",
            "Rwc",
            "Rws",
            "Zwc",
            "Zws",
        ]

        for key in changelist:
            # convert the list into numpy.ndarray
            data = np.array(self["physicslist"][key], dtype=np.float64)
            orgmlen, orgnlen = data.shape
            orgnmin = self["physicslist"].start_index[key.lower()][0]
            orgnmax = orgnmin + orgnlen - 1

            # the new data array with the new shape
            newnlen = new_Ntor * 2 + 1
            newmlen = new_Mpol + 1
            newnmin = -new_Ntor
            newdata = np.zeros([newmlen, newnlen])

            # copy the data over
            for ii in range(newnlen):
                newnid = newnmin + ii

                # if data is provided for this n number
                if newnid >= orgnmin and newnid <= orgnmax:
                    newdata[: min(orgmlen, newmlen), ii] = data[
                        : min(orgmlen, newmlen), newnid - orgnmin
                    ]

            self["physicslist"][key] = newdata.tolist()
            self["physicslist"].start_index[key.lower()][0] = -new_Ntor

        # We will need to update the size of Rac, etc and their indexing
        changelist = self.axis_keys

        for key in changelist:
            if not isinstance(self["physicslist"][key], list):
                data = np.array([self["physicslist"][key]], dtype=np.float64)
            else:
                data = np.array(self["physicslist"][key], dtype=np.float64)

            newnlen = new_Ntor + 1
            orgnlen = data.shape[0]
            newdata = np.zeros([newnlen])

            newdata[: min(orgnlen, newnlen)] = data[: min(orgnlen, newnlen)]
            self["physicslist"][key.lower()] = newdata.tolist()

        # update start index and self
        self._Ntor = new_Ntor
        self._Mpol = new_Mpol
        self["physicslist"]["mpol"] = new_Mpol
        self["physicslist"]["ntor"] = new_Ntor

    def get_interface_guess(self, m, n, ivol, key="Rbc"):
        """Get the guess of the interface Fourier harmonic
        parameters:
            m,n -- the m and n number of the guess, must be within the allowed Mpol and Ntor range
                   the n number is the one without multiplying by Nfp
            ivol -- which volume, Python convention, starting from 0
            key -- which item, can be 'Rbc', 'Zbs', 'Rbs', 'Zbc'
        Returns:
            guess -- the initial guess of the interface harmonic used in SPEC
        """
        if ivol >= self._Nvol or ivol < 0:
            raise ValueError("ivol must be between 0 and Nvol-1")
        if (m, n) not in self.interface_guess.keys():
            raise ValueError("unknown m or n")
        if key not in ["Rbc", "Rbs", "Zbc", "Zbs"]:
            raise ValueError("key must be in ['Rbc', 'Rbs', 'Zbc', 'Zbs']")

        return self.interface_guess[(m, n)][key][ivol]

    def set_interface_guess(self, value, m, n, ivol, key="Rbc"):
        """Set the guess of the interface Fourier harmonic
        parameters:
            value -- the value that one wants to set
            m,n -- the m and n number of the guess, must be within the allowed Mpol and Ntor range
                   the n number is the one without multiplying by Nfp
            ivol -- which volume, Python convention, starting from 0
            key -- which guess, can be 'Rbc', 'Zbs', 'Rbs', 'Zbc'
        """
        if ivol >= self._Nvol or ivol < 0:
            raise ValueError("ivol must be between 0 and Nvol-1")
        if m > self._Mpol or m < 0:
            raise ValueError("0 <= m <= Mpol")
        if n > self._Ntor or n < -self._Ntor:
            raise ValueError("-Ntor <= n <= Ntor")
        if key not in ["Rbc", "Rbs", "Zbc", "Zbs"]:
            raise ValueError("key must be in ['Rbc', 'Rbs', 'Zbc', 'Zbs']")

        if (m, n) not in self.interface_guess.keys():
            # add a new item
            self.interface_guess[(m, n)] = dict()
            for key in ["Rbc", "Rbs", "Zbc", "Zbs"]:
                self.interface_guess[(m, n)][key] = np.zeros([self._Nvol])

        self.interface_guess[(m, n)][key][ivol] = value

    def remove_interface_guess(self, m, n):
        """Remove the guess of the interface Fourier harmonic with some m,n
        parameters:
            m,n -- the m and n number of the guess, must be within the allowed Mpol and Ntor range
                   the n number is the one without multiplying by Nfp
        """

        if (m, n) not in self.interface_guess.keys():
            raise ValueError("unknown m or n")
        else:
            del self.interface_guess[(m, n)]

    def _rectify_namelist(self):
        """correct the size/type of the namelist objects and so on"""

        for key in self.zero_index_keys:
            if key.lower() in self["physicslist"].keys():
                if not isinstance(self["physicslist"][key], list):
                    self["physicslist"][key] = [self["physicslist"][key]]

            if key.lower() in self["diagnosticslist"].keys():
                if not isinstance(self["diagnosticslist"][key], list):
                    self["diagnosticslist"][key] = [self["diagnosticslist"][key]]

        for key in self.one_index_keys:
            if key.lower() in self["physicslist"].keys():
                if not isinstance(self["physicslist"][key], list):
                    self["physicslist"][key] = [self["physicslist"][key]]

        # check the boundary namelist, create non-existence items
        for key in self.boundary_keys:
            if key.lower() not in self["physicslist"].keys():
                self["physicslist"][key] = np.zeros(
                    [self._Mpol, self._Ntor * 2 + 1], dtype=np.float
                ).tolist()
                self["physicslist"].start_index[key.lower()] = [-self._Ntor, 0]
            else:
                # fill in all the none terms as zeros
                for item1 in self["physicslist"][key]:
                    for idx, item2 in enumerate(item1):
                        if item2 is None:
                            item1[idx] = 0.0

        # check the axis namelist, create non-existence items
        for key in self.axis_keys:
            if key.lower() not in self["physicslist"].keys():
                self["physicslist"][key] = np.zeros(
                    [self._Ntor + 1], dtype=np.float
                ).tolist()
                self["physicslist"].start_index[key.lower()] = [0]
            else:
                # check if the item is a list. if not, make it a list
                if not isinstance(self["physicslist"][key], list):
                    self["physicslist"][key] = [self["physicslist"][key]]
                # fill in all the none terms as zeros
                for idx, item1 in enumerate(self["physicslist"][key]):
                    if item1 is None:
                        self["physicslist"][key][idx] = 0.0

    def _dump_to_namelist(self, spec_hdf5_subgroup, target_namelist):
        """dump the properties in the SPEC hdf5 group to the namelist
        parameters:
            spec_hdf5_subgroup -- the SPEC property generated from hdf5 group
            target_namelist -- the target namelist object
        """

        for key in dir(spec_hdf5_subgroup):
            # add to the namelist if it is not starting with '_' (internal python functions)
            if not key.startswith("_") and not(key in self._not_to_dump_list):
                if key in self.boundary_keys:
                    # take care of all the boundary inputs
                    data = getattr(spec_hdf5_subgroup, key)
                    # we only need half of the mpol, since the data is -mpol:mpol
                    mdim = data.shape[0]
                    ndim = data.shape[1]
                    mdim = int((mdim - 1) / 2)
                    ndim = int((ndim - 1) / 2)

                    target_namelist[key] = data[mdim:, :].tolist()
                    target_namelist.start_index[key.lower()] = [-ndim, 0]

                elif key in [
                    "LreadGF",
                    "Ltiming",
                    "LHevalues",
                    "LHevectors",
                    "LHmatrix",
                ]:
                    # take care of all the bool inputs
                    target_namelist[key] = getattr(spec_hdf5_subgroup, key).item() == 1
                else:
                    if isinstance(getattr(spec_hdf5_subgroup, key), np.ndarray):
                        target_namelist[key] = getattr(spec_hdf5_subgroup, key).tolist()
                    else:
                        target_namelist[key] = getattr(spec_hdf5_subgroup, key)

    def _write_interface_guess(self, file_object):
        """Write the initialization of the interface to filename"""

        for harmonics, data in self.interface_guess.items():
            # write m,n
            output_str = "{:6d} {:6d} "
            output_str = output_str.format(harmonics[0], harmonics[1])
            file_object.write(output_str)

            # write the data
            for ii in range(self._Nvol):
                output_str = "{:23.15e} {:23.15e} {:23.15e} {:23.15e}"
                output_str = output_str.format(
                    data["Rbc"][ii], data["Zbs"][ii], data["Rbs"][ii], data["Zbc"][ii]
                )
                file_object.write(output_str)

            # end the line
            file_object.write("\n")

    def _read_interface_guess(self, file_object):
        """Read the initial guess of the interface
        parameters:
            file_object -- the file object, created using open(filename,'r')
        """

        # we can read all lines into the memory since the namelist file is usually not that big
        file_object.seek(0)
        lines = file_object.readlines()

        # first pass, we go through the file to locate the last '/' symbol, which indicates the end of the namelist
        slash_counter = 0
        for line_counter, line in enumerate(lines):
            if "/" in line:
                slash_counter = line_counter

        self.interface_guess = dict()

        # now go back and start from the next line exactly at the location of the last '/'
        for ii in range(slash_counter + 1, len(lines)):

            # split the line into objects
            line_split = lines[ii].split()

            # ignore empty lines
            if len(line_split) == 0:
                break

            # check if this line meet our expectation
            valid_line = True
            # first, the number of items should be equal or greater than nvol * 4 + 2
            if not len(line_split) == self._Nvol * 4 + 2:
                valid_line = False
            else:
                # the first item should be m, check if this is correct
                m = int(line_split[0])
                if m < 0:
                    valid_line = False

            # if something wrong is detected, report a warning message, jump this line
            if not valid_line:
                print("Initial guess of the interface geometry ignored: line ", ii + 1)
                break

            # extract the n number
            n = int(line_split[1])

            # now parse the line, put the data in a dictionary
            self.interface_guess[(m, n)] = dict()
            self.interface_guess[(m, n)]["Rbc"] = [
                float(line_split[2 + lvol * 4]) for lvol in range(self._Nvol)
            ]
            self.interface_guess[(m, n)]["Zbs"] = [
                float(line_split[2 + lvol * 4 + 1]) for lvol in range(self._Nvol)
            ]
            self.interface_guess[(m, n)]["Rbs"] = [
                float(line_split[2 + lvol * 4 + 2]) for lvol in range(self._Nvol)
            ]
            self.interface_guess[(m, n)]["Zbc"] = [
                float(line_split[2 + lvol * 4 + 3]) for lvol in range(self._Nvol)
            ]

    def _generate_namelist_from_output(self, spec_hdf5):
        """initialize the namelist from SPEC output:
        parameter:
            spec_hdf5 -- A SPEC object generated by reading SPEC hdf5 output
        """

        # create the namelists
        for key in [
            "physicslist",
            "numericlist",
            "locallist",
            "globallist",
            "diagnosticslist",
            "screenlist",
        ]:
            self[key] = Namelist()

        self._not_to_dump_list = dir(spec_hdf5)
        with spec_hdf5.input as i:
            self._dump_to_namelist(i.physics, self["physicslist"])
            self._dump_to_namelist(i.numerics, self["numericlist"])
            self._dump_to_namelist(i.local, self["locallist"])
            self._dump_to_namelist(i.global1, self["globallist"])
            self._dump_to_namelist(i.diagnostics, self["diagnosticslist"])

            # we don't dump screenlist since it is not saved in the hdf5
            # self._dump_to_namelist(i.screen, self['screenlist'])

        # now we should mind some variables that are important: we include them in the object itself and need to monitor them
        self._Mpol = self["physicslist"]["mpol"]
        self._Ntor = self["physicslist"]["ntor"]
        self._Nvol = self["physicslist"]["nvol"]

        # replace some namelist objects by those from the output
        with spec_hdf5.output as o, spec_hdf5.input.physics as p:
            # 1. replace guess of the geometry axis
            self["physicslist"]["Rac"] = o.Rbc[0, : p.Ntor + 1].tolist()
            self["physicslist"]["Zas"] = o.Zbs[0, : p.Ntor + 1].tolist()
            self["physicslist"]["Ras"] = o.Rbs[0, : p.Ntor + 1].tolist()
            self["physicslist"]["Zac"] = o.Zbc[0, : p.Ntor + 1].tolist()

            # 2. replace the boundary
            for ii in range(spec_hdf5.output.mn):
                mm = o.im[ii]
                nn = int((o.in_[ii]) / p.Nfp) + self._Ntor
                self["physicslist"]["Rbc"][mm][nn] = o.Rbc[p.Nvol, ii]
                self["physicslist"]["Zbs"][mm][nn] = o.Zbs[p.Nvol, ii]
                self["physicslist"]["Rbs"][mm][nn] = o.Rbs[p.Nvol, ii]
                self["physicslist"]["Zbc"][mm][nn] = o.Zbc[p.Nvol, ii]

            # 3. replace some physics quantities
            output_list = [
                "mu",
                "pflux",
                "helicity",
                "adabatic",
                "iota",
                "oita",
                "Isurf",
                "Ivolume",
            ]
            for key in output_list:
                if key in dir(o):
                    self["physicslist"][key] = getattr(o, key).tolist()

            # 4. generate the guess of the interface
            self.interface_guess = dict()
            for ii in range(o.mn):
                m = o.im[ii]
                n = int(o.in_[ii] / p.Nfp)
                self.interface_guess[(m, n)] = dict()
                self.interface_guess[(m, n)]["Rbc"] = o.Rbc[1:, ii]
                self.interface_guess[(m, n)]["Zbs"] = o.Zbs[1:, ii]
                self.interface_guess[(m, n)]["Rbs"] = o.Rbs[1:, ii]
                self.interface_guess[(m, n)]["Zbc"] = o.Zbc[1:, ii]

    def _interpolate_guess(self, ivol, tflux):
        """Interpolated interface harmonics guess
        parameters:
            ivol -- where is the new volume
            tflux -- the new tflux
        """

        # The way we interpolate
        Lcoordinatesingularity = False

        if self["physicslist"]["Igeometry"] >= 2:
            Lslab = False
            if ivol == 0 or ivol == self._Nvol:
                Lcoordinatesingularity = True
        else:
            Lslab = True

        if ivol == self._Nvol:
            if self["physicslist"]["Lfreebound"] == 0:
                Lextrapolate = True
            else:
                Lextrapolate = False
        else:
            Lextrapolate = False

        # if the inserted volume is outside the last volume, we need to extrapolate
        if Lextrapolate:
            raise RuntimeError[
                "Extrapolating outside the plasma boundary is not supported"
            ]
        else:
            if Lcoordinatesingularity:

                r_left = 0.0
                r_right = np.sqrt(self["physicslist"]["tflux"][0])
                r_int = np.sqrt(tflux)

                for key in ["Rbc", "Zbs", "Rbs", "Zbc"]:
                    self._interpolate_guess_singular_each(
                        key, ivol, r_left, r_right, r_int
                    )

            else:
                if ivol == 0:
                    tflux_left = 0.0
                else:
                    tflux_left = self["physicslist"]["tflux"][ivol - 1]

                tflux_right = self["physicslist"]["tflux"][ivol]

                if Lslab:  # for a slab, radius is the same as tflux
                    r_left = tflux_left
                    r_right = tflux_right
                    r_int = tflux
                else:  # for cylinder or toroidal, radius is the same as sqrt(tflux)
                    r_left = np.sqrt(tflux_left)
                    r_right = np.sqrt(tflux_right)
                    r_int = np.sqrt(tflux)

                # interpolate
                for key in ["Rbc", "Zbs", "Rbs", "Zbc"]:
                    self._interpolate_guess_normal_each(
                        key, ivol, r_left, r_right, r_int
                    )

    def _interpolate_guess_normal_each(self, key, ivol, r_left, r_right, r_int):
        """Interpolate the interface harmonic guess, normal way
        parameters:
            key -- which item? 'Rbc', etc
            ivol -- which volume?
            r_left, r_right -- the left and right "radius"
            r_int -- the interpolate "radius"
        """
        for mnkey, item in self.interface_guess.items():
            value_int = self._interpolate_normal(
                item[key], ivol, r_left, r_right, r_int
            )
            item[key] = np.insert(item[key], ivol, value_int)

    def _interpolate_guess_singular_each(self, key, ivol, r_left, r_right, r_int):
        """Interpolate the interface harmonic guess, normal way
        parameters:
            key -- which item? 'Rbc', etc
            ivol -- which volume?
            r_left -- dummy, not used
            r_right -- the right "radius"
            r_int -- the interpolate "radius"
        """
        # replacing 'b' by 'a' we will get the keys for Rac and so on
        key_axis = key.replace("b", "a")

        for mnkey, item in self.interface_guess.items():
            m = mnkey[0]
            n = mnkey[1]
            if m == 0:
                value_axis = self["physicslist"][key_axis][n]
            else:
                value_axis = 0.0
            value_int = self._interpolate_singular(
                item[key], ivol, m, r_right, r_int, value_axis
            )
            item[key] = np.insert(item[key], ivol, value_int)

    def _extrapolate_guess_singular_each(self, key, ivol, r_left, r_right, r_int):
        """Interpolate the interface harmonic guess, normal way
        parameters:
            key -- which item? 'Rbc', etc
            ivol -- which volume?
            r_left -- dummy, not used
            r_right -- the right "radius"
            r_int -- the interpolate "radius"
        """
        # replacing 'b' by 'a' we will get the keys for Rac and so on
        key_axis = key.replace("b", "a")

        for mnkey, item in self.interface_guess.items():
            m = mnkey[0]
            n = mnkey[1]
            if m == 0:
                value_axis = self["physicslist"][key_axis][n]
            else:
                value_axis = 0.0
            value_int = self._interpolate_singular(
                item[key], ivol - 1, m, r_right, r_int, value_axis
            )
            item[key] = np.insert(item[key], ivol, value_int)

    @staticmethod
    def _interpolate_normal(data, ivol, r_left, r_right, r_int):
        """Linear interpolation"""
        if ivol == 0:
            value_left = 0.0
        else:
            value_left = data[ivol - 1]
        value_right = data[ivol]
        # linear interpolate
        value_int = (value_right - value_left) / (r_right - r_left) * (
            r_int - r_left
        ) + value_left

        return value_int

    @staticmethod
    def _interpolate_singular(data, ivol, m, r_right, r_int, value_left):
        """Singular interpolation"""
        value_right = data[ivol]

        # interpolate
        s = float(r_int) / float(r_right)
        if (
            m == 0
        ):  # for m==0, we need to interpolate between axis value and boundary value
            value_int = value_right * s ** 2 + value_left * (1.0 - s ** 2)
        else:  # otherwise, we don't need the axis value
            value_int = value_right * s ** m

        return value_int
