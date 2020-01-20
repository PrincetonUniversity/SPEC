#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 16:10:39 2020

@author: jonathan
"""

#%% Class declarations
class Variable:
    
    name = None
    dtype = None
    defaultValue = None
    rank = 0 # scalars by default
    description = None
    
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
        
    def printDescription_HTML_lists(self):
        for d in self.description:
            if type(d) is list:
                print("<ul>")
                for dl in d:
                    print("<li> "+dl+" </li>")
                print("</ul>")
            else:
                print(d)
    