#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Algorithm Definition Framework

@author: Jonathan Schilling (jonathan.schilling@mail.de)
"""

#%% Class declarations
class Variable:
    
    name = None
    dtype = None
    defaultValue = None
    rank = 0 # scalar by default
    description = None
    isParameter = False
    unit = None
    maximumIndices = None
    startingIndices = None
    
    def __init__(self, name):
        self.name = name
        
    def setDescription(self, description):
        self.description = description
    
    def setType(self, dtype):
        self.dtype = dtype
    
    def setRank(self, rank):
        self.rank = rank
        
    def setDefaultValue(self, defaultValue):
        self.defaultValue = defaultValue
        
    def setUnit(self, unit):
        self.unit = unit
    
    def setMaximumIndices(self, maximumIndices):
        self.maximumIndices = maximumIndices
        
    def setStartingIndices(self, startingIndices):
        self.startingIndices = startingIndices
    
    def setIsParameter(self, isParameter):
        self.isParameter = isParameter
    
    def printDescription_HTML_lists(self):
        for d in self.description:
            if type(d) is list:
                print("<ul>")
                for dl in d:
                    print("<li> "+dl+" </li>")
                print("</ul>")
            else:
                print(d)
    