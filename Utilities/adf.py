#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Algorithm Definition Framework

@author: Jonathan Schilling (jonathan.schilling@mail.de)
"""

# Class declarations
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

def toDoc(desc):
    if   type(desc) is dict:
        return desc_dictToDoc(desc)
    elif type(desc) is list:
        return desc_listToDoc(desc)
    elif type(desc) is str:
        return desc
    else:
        raise TypeError("what is this that you want to document of type "+str(type(desc)))

def desc_dictToDoc(desc_dict):
    if type(desc_dict) is not dict:
        raise RuntimeError("desc_dictToDoc was called with "+str(type(desc_dict))+" instead of dict")
    result = ""
    for key in desc_dict.keys():
        if type(key) is not str:
            raise RuntimeError("desc_dictToDoc was given a dict with key type "+str(type(desc_dict))+" instead of str keys")
        result += key+"\n"
        result += toDoc(desc_dict[key])
    return result

def desc_listToDoc(desc_list):
    htmlListsPutThere = "<ul>\n"
    for item in desc_list:
        itemStr = toDoc(item)
        liIndented = "<li> "+indented(5, itemStr, " ")[5:]+" </li>"
        htmlListsPutThere += liIndented+"\n"
    htmlListsPutThere += "</ul>"
    return htmlListsPutThere


#
#
#descNtor = {r"Control flag for solution of Beltrami equation":
#                                      [ r"if \c LBeltrami = 1,3,5 or 7, (SQP) then the Beltrami field in each volume is constructed"+"\n"
#                                       +r"by minimizing the magnetic energy with the constraint of fixed helicity;"+"\n"
#                                       +r"this is achieved by using sequential quadratic programming as provided by \c E04UFF ."+"\n"
#                                       +r"This approach has the benefit (in theory) of robustly constructing minimum energy solutions"+"\n"
#                                       +r"when multiple, i.e. bifurcated, solutions exist.",
#                                        r"if \c LBeltrami = 2,3,6 or 7, (Newton) then the Beltrami fields are constructed by employing a standard Newton method"+"\n"
#                                       +r"for locating an extremum of"+"\n"
#                                       +r"\f$F\equiv \int B^2 dv - \mu (\int {\bf A}\cdot{\bf B}dv-{\cal K})\f$,"+"\n"
#                                       +r"where \f$\mu\f$ is treated as an independent degree of freedom similar to the parameters describing the vector potential"+"\n"
#                                       +r"and \f${\cal K}\f$ is the required value of the helicity;"+"\n"
#                                       +r"this is the standard Lagrange multipler approach for locating the constrained minimum;"+"\n"
#                                       +r"this method cannot distinguish saddle-type extrema from minima, and which solution that will be obtained depends on the initial guess",
#                                        r"if \c LBeltrami = 4,5,6 or 7, (linear) it is assumed that the Beltrami fields are parameterized by \f$\mu\f$;"+"\n"
#                                       +r"in this case, it is only required to solve \f$\nabla \times {\bf B} = \mu {\bf B}\f$ which reduces to a system of linear equations;"+"\n"
#                                       +r"\f$\mu\f$ may or may not be adjusted iteratively, depending on \c Lconstraint,"+"\n"
#                                       +r"to satisfy either rotational-transform or helicity constraints",
#                                       {r"for flexibility and comparison, each of the above methods can be employed; for example:":
#                                        [r"if \c LBeltrami = 1, only the SQP    method will be employed",
#                                         r"if \c LBeltrami = 2, only the Newton method will be employed",
#                                         r"if \c LBeltrami = 4, only the linear method will be employed",
#                                         r"if \c LBeltrami = 3, the SQP and the Newton method are used",
#                                         r"if \c LBeltrami = 5, the SQP and the linear method are used",
#                                         r"if \c LBeltrami = 6, the Newton and the linear method are used",
#                                         r"if \c LBeltrami = 7, all three methods will be employed"]
#                                        }]
#                                      }
#
#print(toDoc(descNtor))