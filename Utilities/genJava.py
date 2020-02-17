#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 22:49:47 2019

@author: Jonathan Schilling (jonathan.schilling@ipp.mpg.de)
"""


#%% prepare for code generation

def indented(tabs, lines, indentationChar="\t"):
    indentation = ""
    for i in range(tabs):
        indentation += indentationChar
    indented = ''
    if '\n' in lines.strip():
        for line in lines.split('\n'):
            if line != '':
                indented += indentation+line+'\n'
    else:
        indented = indentation+lines#.strip()
    return indented

def indent(tabs, lines, indentationChar="\t"):
    return tabs+1, indented(tabs, lines, indentationChar)

def unindent(tabs, lines, indentationChar="\t"):
    return tabs-1, indented(tabs, lines, indentationChar)


#%% create a Java reading routine for a HDF5 file
from Hdf5File import Group, Dataset

def javaClassName(name):
    """Make a name like "asdf_adsf" into camel case at the locations of "_" and start with an uppercase letter."""
    while "_" in name:
        idx = name.index("_")
        name = name[:idx]+name[idx+1:idx+2].upper()+name[idx+2:]
    return name[0].upper()+name[1:]
        
def javaVarName(name):
    """Make a name like "asdf_adsf" into camel case at the locations of "_" and start with a lowercase letter."""
    while "_" in name:
        idx = name.index("_")
        name = name[:idx]+name[idx+1:idx+2].upper()+name[idx+2:]
    #return name[0].lower()+name[1:] # allow user to actually specify exact variable name
    return name

def javaDtype(dtype):
    """Translate the dtypes from the definition into valid Java primitive types or custom classes generated for compund datatypes."""
    if dtype=='int' or dtype=='double' or dtype=='boolean':
        return dtype
    else:
        return javaClassName(dtype)

def javaGenClassFromGroup(tabs, group, static=True):
    """Generate Java source code defining a corresponding class from the definition of a Group, recursing into the items.
    
    tabs   -- number of indentation marks to prepend to every line of the source code
    group  -- Group defininition of which a Java class should be generated
    static -- select if the generated class should be static or not (defaults to True)
    """
    classCode = "\n"
    readCodes = []
    if group.description is not None:
        classCode += indented(tabs, "/** "+group.description+" */\n")
    if static:
        tabs, decl = indent(tabs, "public static class "+javaClassName(group.getFullName().replace("/", "_"))+" {\n")
    else:
        tabs, decl = indent(tabs, "public class "+javaClassName(group.getFullName().replace("/", "_"))+" {\n")
    classCode += decl
    constructorPart=''
    memberPart = ''
    numComplexMembers = 0
    for item in group.items:
        if item.description is not None:
            memberPart += indented(tabs, "/** "+item.description+" */\n")
        memberPart += indented(tabs, 'public ')
        if type(item)==Dataset:
            memberPart += javaDtype(item.dtype)
            readCodes.append(javaRead(tabs, item))
        else:
            memberPart += javaClassName(item.getFullName().replace("/", "_"))
            constructorPart += indented(tabs+1, item.name)+" = new "+javaClassName(item.getFullName().replace("/", "_"))+"();\n"
            numComplexMembers+=1
        if type(item) == Dataset and item.getRank()>0:
            for i in range(item.getRank()):
                memberPart += '[]'
        memberPart += ' '+item.name+';\n'
    if numComplexMembers>0:
        classCode += indented(tabs, "/** initialize complex datatypes */\n")
        classCode += indented(tabs, 'public '+javaClassName(group.getFullName().replace("/", "_"))+'() {\n')
        classCode += constructorPart
        classCode += indented(tabs, '}\n\n')
    classCode += memberPart
    tabs -= 1
    classCode += indented(tabs, '} // end of '+javaClassName(group.getFullName().replace("/", "_")))
    return classCode, readCodes

def javaRead(tabs, dataset):
    """Generate Java code that reads the given Dataset from a NetcdfFile 'file'.
    
    dataset -- Dataset that should be read
    """
    varName = dataset.getFullName()
    javaName = varName[1:].replace("/", ".")
    rank = dataset.getRank()
    readCode = ''
    if dataset.dtype=='int':
        if rank==0:
            readCode = '{javaName} = file.findVariable("{varName}").readScalarInt()'
        else:
            readCode = '{javaName} = (int'
            for r in range(rank):
                readCode += '[]'
            readCode += ')file.findVariable("{varName}").read()'
            if rank==1:
                readCode += '.get1DJavaArray(DataType.INT)'
            else:
                readCode += '.copyToNDJavaArray()'
    elif dataset.dtype=='double':
        if rank==0:
            readCode = '{javaName} = file.findVariable("{varName}").readScalarDouble()'
        else:
            readCode = '{javaName} = (double'
            for r in range(rank):
                readCode += '[]'
            readCode += ')file.findVariable("{varName}").read()'
            if rank==1:
                readCode += '.get1DJavaArray(DataType.DOUBLE)'
            else:
                readCode += '.copyToNDJavaArray()'
    elif dataset.dtype=='boolean':
        if rank==0:
            readCode = '{javaName} = (file.findVariable("{varName}").readScalarInt() > 0 ? true : false)'
        else:
            print(dataset.getFullName()+" reading not implemented yet")
            readCode = '// read {varName} into {javaName}'
#            readCode = 'int '
#            dimBrackets = ''
#            firstElems = []
#            for r in range(rank):
#                dimBrackets += '[]'
#                if r>1:
#                    firstElems.append(firstElems[r-1]+'[0]')
#                else:
#                    firstElems.append('')
#            readCode += dimBrackets+' {javaName}_int = (int'+dimBrackets
#            readCode += ')file.findVariable("{varName}").read()'
#            if rank==1:
#                readCode += '.get1DJavaArray(DataType.INT)'
#            else:
#                readCode += '.copyToNDJavaArray()'
#            readCode += ';\n'
#            readCode += indented(tabs, '{javaName} = new boolean')
#            for r in range(rank):
#                readCode += '[{javaName}_int'+firstElems[r]+'.length];\n'
    else:
        # custom datatype
        print(dataset.getFullName()+" reading not implemented yet")
    return readCode.format(javaName=javaName, varName=varName)+';\n'


#%% document who created the reading routines when on which machine

from datetime import datetime
import getpass
import platform

# dd/mm/YY H:M:S in UTC
now_string = datetime.utcnow().strftime('%d/%m/%Y %H:%M:%S UTC')
username = getpass.getuser()
hostname = platform.node()

creation_tag = 'auto-created by a user called \''+username+'\' on a machine called \''+hostname+'\' at '+now_string

#%% actually generate Java class for reading SPEC output files
def genJavaReader(outdir, packageName, className, s):
    
    # we need to reverse the definition order so that types which are used inside other types
    # are already defined when used
    reverse_rootStack =  []
        
    rootStack = []
    rootStack.append(s.rootGroup)
    while len(rootStack)>0:
        currentItem = rootStack[-1]
        rootStack = rootStack[:-1]
        
        if currentItem is not s.rootGroup:
            reverse_rootStack.append(currentItem)
        if type(currentItem)==Group:
            for item in currentItem.items:
                rootStack.append(item)
    
    
    javaFilename = outdir+className+".java"
    print("creating Java reading class into '"+javaFilename+"'")
    
    # begin code for root group (== enclosing class)
    f=open(javaFilename, "w")
    tabs=0
    
    f.write("""package """+packageName+""";
// AUTO-GENERATED; DO NOT COMMIT CHANGES TO THIS FILE !
// """+creation_tag+"""

import java.io.IOException;
import java.util.Locale;
import ucar.ma2.DataType;
import ucar.nc2.NetcdfFile;

""")
    
    
    
    rootClassCode = ""
    if s.rootGroup.description is not None:
        rootClassCode += indented(tabs, "/** "+s.rootGroup.description+" */\n")
    tabs, decl = indent(tabs, "public class "+className+" {\n")
    rootClassCode += decl
    numComplexMembers = 0
    f.write(rootClassCode)
    
    readParts=[]
    
    # add nested groups
    while len(reverse_rootStack)>0:
        currentItem = reverse_rootStack[-1]
        reverse_rootStack = reverse_rootStack[:-1]
        
        if type(currentItem)==Group:
            defCode, readCodes = javaGenClassFromGroup(tabs, currentItem)
            f.write(defCode+'\n')
            for readCode in readCodes:
                readParts.append(readCode)
            numComplexMembers+=1
    
    # end code for root group (== enclosing class)
    constructorPart=''
    memberPart = ''
    rootClassCode = ""
    for item in s.rootGroup.items:
        if item.description is not None:
            memberPart += indented(tabs, "/** "+item.description+" */\n")
        memberPart += indented(tabs, "public ")
        if type(item)==Dataset:
            memberPart += javaDtype(item.dtype)
            readParts.append(javaRead(tabs, item))
        else:
            memberPart += javaClassName(item.getFullName().replace("/", "_"))
            constructorPart += indented(tabs+1, item.name+" = new "+javaClassName(item.getFullName().replace("/", "_"))+"();\n")
            numComplexMembers+=1
        if type(item) == Dataset and item.getRank()>0:
            for i in range(item.getRank()):
                memberPart += '[]'
        memberPart += ' '+item.name+';\n'
    rootClassCode += "\n"
    # constructor to initialize complex members
    if numComplexMembers>0:
        rootClassCode += indented(tabs, "/** Initialize complex datatypes. */\n")
        rootClassCode += indented(tabs, 'public '+className+'() {\n')
        rootClassCode += constructorPart
        rootClassCode += indented(tabs, '}\n')
        
    # constructors to load data from file
    rootClassCode += "\n"
    rootClassCode += indented(tabs, "/**\n")
    rootClassCode += indented(tabs, " * Initalize complex datatypes and load "+className+" contents from a HDF5 file identified by {@code filename}.\n")
    rootClassCode += indented(tabs, " * @param filename path to the HDF5 file to load\n")
    rootClassCode += indented(tabs, " */\n")
    tabs, line = indent(tabs, "public "+className+"(String filename) {\n")
    rootClassCode += line
    rootClassCode += indented(tabs, "this();\n")
    tabs, line = indent(tabs, "try {\n")
    rootClassCode += line
    rootClassCode += indented(tabs, "NetcdfFile file = NetcdfFile.open(filename);\n")
    rootClassCode += indented(tabs, "loadFrom(file);\n")
    rootClassCode += indented(tabs, "file.close();\n")
    rootClassCode += indented(tabs-1, "} catch (IOException e) {\n")
    rootClassCode += indented(tabs, "e.printStackTrace();\n")
    tabs -= 1
    rootClassCode +=   indented(tabs, "}\n")
    tabs -= 1
    rootClassCode +=   indented(tabs, "}\n")
    rootClassCode +=   "\n"
    rootClassCode +=   indented(tabs, "/**\n")
    rootClassCode +=   indented(tabs, " * Initalize complex datatypes and load "+className+" contents from an already-open NetCDF file identified by {@code file}.\n")
    rootClassCode +=   indented(tabs, " * @param file open file to load the data from\n")
    rootClassCode +=   indented(tabs, " */\n")
    tabs, line     =   indent  (tabs, "public "+className+"(NetcdfFile file) {\n") ; rootClassCode += line
    rootClassCode +=   indented(tabs, "this();\n")
    tabs, line     =   indent  (tabs, "try {\n")                     ; rootClassCode += line
    tabs, line     = unindent  (tabs, "loadFrom(file);\n")           ; rootClassCode += line
    tabs, line     =   indent  (tabs, "} catch (IOException e) {\n") ; rootClassCode += line
    tabs, line     = unindent  (tabs, "e.printStackTrace();\n")      ; rootClassCode += line
    tabs, line     = unindent  (tabs, "}\n")                         ; rootClassCode += line
    rootClassCode += indented  (tabs, "}\n")
    rootClassCode += "\n"
    rootClassCode += memberPart
    f.write(rootClassCode)
    
    # definitions part is done; now for the reading routines
    rootClassCode = "\n"
    rootClassCode += indented(tabs, "/**\n")
    rootClassCode += indented(tabs, " * Load "+className+" contents from an already-open NetCDF file identified by {@code file}.\n")
    rootClassCode += indented(tabs, " * @param file open file to load the data from\n")
    rootClassCode += indented(tabs, " * @return initialized "+className+" object\n")
    rootClassCode += indented(tabs, " */\n")
    tabs, line = indent(tabs, "public "+className+" loadFrom(NetcdfFile file) throws IOException {\n")
    rootClassCode += line
    
    # here goes the magic that actually loads the data from the file
    
    for readPart in readParts:
        rootClassCode += indented(tabs, readPart)
    
    
    
    
    
    tabs, line = unindent(tabs, "return this;\n")
    rootClassCode += line
    rootClassCode += indented(tabs, "}\n")
    
    f.write(rootClassCode)
    
    
    f.write("""
    public static void main(String[] args) {
        SpecOutput s = new SpecOutput("/home/jonathan/Uni/04_PhD/00_programs/SPEC/SPEC/InputFiles/TestCases/G3V02L1Fi.001.h5");
        System.out.printf(Locale.ENGLISH, "SPEC version: %.2f\\n", s.version);
    }
""")
    
    
    ## closing brace
    tabs -= 1
    f.write(indented(tabs, '} // end of '+className+"\n"))
    
    
    
    f.close()