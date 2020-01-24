#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Algorithm Definition Framework

@author: Jonathan Schilling (jonathan.schilling@mail.de)
"""

# Class declarations

# data container class for all information specifying a variable
class Variable:
    
    name = None
    description = None
    dtype = None
    defaultValue = None
    rank = 0 # scalar by default
    isParameter = False
    unit = None
    startingIndices = None
    maximumIndices = None
    
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
    
    def setStartingIndices(self, startingIndices):
        self.startingIndices = startingIndices
    
    def setMaximumIndices(self, maximumIndices):
        self.maximumIndices = maximumIndices
    
    def setIsParameter(self, isParameter):
        self.isParameter = isParameter

# data container class for all information specifying a namelist
class Namelist(object):
    
    name = None
    description = None
    variables = None
    
    def __init__(self, name):
        self.name = name
        self.variables = []
    
    def setDescription(self, description):
        self.description = description
    
    def addVariable(self, var):
        if type(var) is not Variable:
            raise TypeError("type of var is not Variable but '"+str(type(var))+"'")
        self.variables.append(var)

    def addVariables(self, listOfVars):
        if type(listOfVars) is not list:
            raise TypeError("type of listOfVars is not list but '"+str(type(listOfVars))+"'")
        for var in listOfVars:
            self.addVariable(var)







# source code generation utility methods

# indent a string (which might consist of multiple lines) by a given number of tabs or
# some other given character
def indented(tabs, strInput, indentationChar="\t"):
    indentation = ''
    for i in range(tabs):
        indentation += indentationChar
    indented = ''
    if "\n" in strInput:
        lines = strInput.split("\n")
        for line in lines[:-1]:
            indented += indentation+line+"\n"
        indented += indentation+lines[-1]
        if strInput[-1] == "\n":
            indented += "\n"
    else:
        indented = indentation+strInput
    return indented

# convert the description item from a Variable into the corresponding documentation
def toDoc(desc):
    if   type(desc) is str:
        return desc
    elif type(desc) is dict:
        return desc_dictToDoc(desc)
    elif type(desc) is list:
        return desc_listToDoc(desc)
    
    elif desc is not None:
        raise TypeError("what is this that you want to document of type "+str(type(desc))+"?")
    else:
        return ""

# convert a dict from a Variable's description into the corresponding documentation
def desc_dictToDoc(desc_dict):
    if type(desc_dict) is not dict:
        raise RuntimeError("desc_dictToDoc was called with "+str(type(desc_dict))+" instead of dict")
    result = ""
    iKey=0
    for key in desc_dict.keys():
        if type(key) is not str:
            raise RuntimeError("desc_dictToDoc was given a dict with key type "+str(type(desc_dict))+" instead of str keys")
        if iKey>0:
            result += "\n"
        result += key
        if desc_dict[key] is not None:
            result += "\n"+toDoc(desc_dict[key])
        iKey+=1
    return result

# convert a list from a Variable's description into the corresponding documentation
def desc_listToDoc(desc_list):
    htmlListsPutThere = "<ul>\n"
    for item in desc_list:
        itemStr = toDoc(item)
        # indent the item content by length of "<li> " so that it is nicely aligned => 5
        # first item shall not be indented => [5:]
        liIndented = "<li> "+indented(5, itemStr, " ")[5:]+" </li>"
        htmlListsPutThere += liIndented+"\n"
    htmlListsPutThere += "</ul>"
    return htmlListsPutThere

