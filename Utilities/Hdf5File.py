#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 22:48:32 2019

@author: Jonathan Schilling (jonathan.schilling@ipp.mpg.de)
"""

class Datatype:
    """A class to define a (custom) datatype in a HDF5 file."""
    
    parent = None
    name = None
    items = None
    
    def __init__(self, parent, name, dtype, rank=0):
        """Define a new Datatype.
        
        parent     -- Group or Hdf5File into which this Datatype belongs
        name       -- identifier for this Datatype
        """
        self.parent = parent
        self.name = name
        self.dtype = dtype
        self.rank = rank
        self.parent.add(self)
        
    def getFullName(self):
        """Get the full path inside the Hdf5File to this Dataset."""
        return self.parent.getFullName()+"/"+self.name
    
    def setDescription(self, description):
        """Provide a description of this Dataset to put into a HDF5 attribute or into the source code comments."""
        self.description = description

class Dataset:
    """A class to define a dataset in a HDF5 file."""
    
    parent = None
    name = None
    srcNames = None
    dtype = None
    dimensions = None
    description = None
    storageLocation_module = None
    storageLocation_name   = None
    routines = None
    indexMapping = None
    defaultValue = None
    
    def __init__(self, parent, name, dtype, rank=0):
        """Define a new Dataset.
        
        parent     -- Group or Hdf5File into which this Dataset belongs
        name       -- identifier for this Dataset
        dtype      -- data type in Java notation (currently one of int, double)
        dimensions -- list of matrix dimensions; number of elements determines rank of this object
        """
        self.parent = parent
        self.name = name
        self.dtype = dtype
        self.rank = rank
        self.parent.add(self)
        
    def getFullName(self):
        """Get the full path inside the Hdf5File to this Dataset."""
        return self.parent.getFullName()+"/"+self.name
    
    def getRank(self):
        """Get the rank of the Dataset, i.e. the number of matrix dimensions."""
        return self.rank
    
    def setSourceNames(self, srcNames):
        """Specify the name of this Dataset in the source code.
        
        srcNames['r'] -- name of this Dataset in the reading routine; defaults to Dataset's name
        srcNames['w'] -- name of this Dataset in the writing routine; defaults to Dataset's name
    """
    
    def setDescription(self, description):
        """Provide a description of this Dataset to put into a HDF5 attribute or into the source code comments."""
        self.description = description
    
    def setDefaultValue(self, defaultValue):
        """Provide a default value for documentation purposes."""
        self.defaultValue = defaultValue
    
    def setStorageLocation(self, module, aliasName=None):
        """Tell the source code generator that the variable corresponding to this dataset is stored somewhere else.
        
        module    -- the Fortran module in which the dataset's data is stored
                     can be a string --> same for reading and writing or a dict:
                     module['r'] is the module name for reading, module['w'] is the module name for writing
                     None or key not present default to the reading/writing module respectively.
        aliasName -- name of the variable this Dataset corresponds to in the remote module
        """
        self.storageLocation_module = module
        if aliasName is not None:
            self.storageLocation_name = aliasName
        else:
            self.storageLocation_name = self.name
            
    def setReadWriteRoutines(self, routines):
        """Specify specialized reading and writing routines that should be used to r/w this Dataset.
        
        routines['r'] -- name of routine for reading; defaults to one single routine for all variables
        routines['w'] -- name of routine for writing; defaults to one single routine for all variables
        """
        self.routines = routines
        
    def setIndexMapping(self, indexMapping):
        """Specify the dimensions of the array/matrix in the source code.
        All data is always written to the HDF5 file with index sets 0:(end-1); hence mapping to -n:n could be necessary.
        
        indexMapping -- list of index specifications in the form "minIndex:maxIndex"; one for each dimension of the data
        """
        self.indexMapping = indexMapping

class Group:
    """A class to define a group in a HDF5 file."""
    
    parent = None
    name = None
    items = None
    description = None
    
    def __init__(self, parent, name):
        """Define a new Group.
        
        parent -- Group or Hdf5File into which this Group belongs
        name   -- identifier for this Group
        """
        self.parent = parent
        self.name   = name
        self.items  = []
        if self.parent is not None:
            self.parent.add(self)
    
    def add(self, item):
        """Add an item to this Group, e.g. a Dataset or another Group."""
        self.items.append(item)
        
    def getFullName(self):
        """Get the full path inside the Hdf5File to this Group."""
        if self.parent is not None:
            return self.parent.getFullName()+"/"+self.name
        else:
            return self.name
    
    def setDescription(self, description):
        """Provide a description of this Dataset to put into a HDF5 attribute or into the source code comments."""
        self.description = description

class Hdf5File:
    """A class to define the contents of a HDF5 file."""
    
    name = None
    rootGroup = None
    
    def __init__(self, name):
        """Define a new HDF5 file.
        
        name -- name of the file
        """
        self.name = name
        self.rootGroup = Group(None, name)
        self.items = []
    
    def add(self, item):
        """Add a Group or a Dataset to the root group of this HDF5 file.
        
        item -- Group or Dataset to add
        """
        self.rootGroup.add(item)
    
    def getFullName(self):
        """Get the full path inside the Hdf5File to the root group, i.e. ''."""
        return ''
    
    def getDatatypes(self):
        """Get a list of all custom datatypes that occur in this HDF5 file."""
        return []
        
        
    def inventory(self):
        """Print a list of all Groups and Datasets in this HDF5 file definition."""
        # pre-order traversal of the content tree starting at the root Group
        grpStack=[self.rootGroup]
        while len(grpStack)>0:
            currentGroup = grpStack[-1]
            grpStack = grpStack[:-1]
            for item in currentGroup.items:
                if type(item)==Dataset:
                    print("D "+item.getFullName())
                elif type(item)==Group:
                    print("G "+item.getFullName())
                    grpStack.append(item)
        
