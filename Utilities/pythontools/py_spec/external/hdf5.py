###############################################################################
# hdf5.py: adaptation of J. Schilling's SPEC script
#
import h5py
import numpy as np # for isscalar
import os          # for path.abspath
import keyword     # for getting python keywords

# Modified from J. Schilling (jonathan.schilling@ipp.mpg.de)'s script for SPEC

class HDF5:
    # use as s = HDF5(filename), e.g. s=HDF5("ext.h5")
    def __init__(self, *args, **kwargs):
        # args[0] should always be the name of a file or an item inside the root object
        # if args[0] is not a filename, kwargs['content'] should be the content to be added
        # as self.`args[0]`

        _content = None
        if kwargs.get('content') == None:
            # assume arg[0] is a filename
            _content = h5py.File(args[0], "r")

            # keep track of which file this object corresponds to
            self.filename = os.path.abspath(args[0])
        elif isinstance(kwargs['content'], h5py.Group):
            _content = kwargs['content']

        if (_content != None):
            for key in _content:
                if isinstance(_content[key], h5py.Group):
                    # recurse into group
                    setattr(self, key, HDF5(content=_content[key]))
                elif isinstance(_content[key], h5py.Dataset): # read dataset
                    if key in keyword.kwlist: # avoid assign python keywords
                        setattr(self, key+'1', _content[key][()])
                    else: # if just one element, use the value directly
                        if len(_content[key][()]) == 1:
                            setattr(self, key, _content[key][0])
                        else:
                            setattr(self, key, _content[key][()])
                    # glue string together
                    if isinstance(_content[key][0], bytes):
                        abc = ''
                        for i in _content[key]:
                            abc += i.decode('utf-8')
                        setattr(self, key, abc)

        if isinstance(_content, h5py.File):
            _content.close()

    # needed for iterating over the contents of the file
    def __iter__(self):
        return iter(self.__dict__)
    def __next__(self):
        return next(self.__dict__)

    # print a list of items contained in this object
    def inventory(self, prefix=""):
        _prefix = ""
        if prefix != "":
            _prefix = prefix+"/"

        for a in self:
            try:
                # recurse into member
                getattr(self, a).inventory(prefix=_prefix+a)
            except:
                # print item name
                print(_prefix+a)
    # delete
    def __del__(self):
        class_name = self.__class__.__name__

