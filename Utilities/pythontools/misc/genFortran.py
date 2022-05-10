#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 15:17:47 2019

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


#%% document who created the reading routines when on which machine

from datetime import datetime
import getpass
import platform

# dd/mm/YY H:M:S in UTC
now_string = datetime.utcnow().strftime('%d/%m/%Y %H:%M:%S UTC')
username = getpass.getuser()
hostname = platform.node()

creation_tag = 'auto-created by a user called \''+username+'\' on a machine called \''+hostname+'\' at '+now_string

#%% generate Fortran type declarations
from Hdf5File import Group, Dataset, Datatype

# datatype in Fortran from specification file
def fortran_dtype(dtype):
    if dtype=='int':
        return 'INTEGER'
    elif dtype=='double':
        return 'DOUBLE PRECISION'
    elif dtype=='boolean':
        return 'LOGICAL'
    else:
        return 'TYPE('+dtype.upper()+')'

# generate custom compound datatype declaration in Fortran
def fortran_genType(name, members):
    ret = 'TYPE '+name+'\n'
    for member in members:
        if type(member) == Group or type(member) == Datatype:
            ret += '    TYPE('+member.name+')'
        else:
            ret += '    '+fortran_dtype(member.dtype)
            if member.rank>0:
                ret += ', ALLOCATABLE'
        ret += ' :: '+member.name
        if type(member) != Group and member.rank>0:
            ret += '('
            for i in range(member.rank):
                if i>0:
                    ret += ',:'
                else:
                    ret += ':'
            ret += ')'
        ret += '\n'
    ret += 'END TYPE '+name
    return ret

# initial code of loading routine
def fortran_startLoader(f):
    f.write("""subroutine loadSpec(s, filename, ierr)
  use hdf5
  implicit none
  type(SpecOutput), intent(inout) :: s                 ! target datastructure
  character(len=*), intent(in)    :: filename          ! filename to load
  integer, intent(out), optional  :: ierr              ! error flag; .eq.0 if ok
  integer                         :: hdfier            ! error flag for HDF5 API calls
  integer(hid_t)                  :: file_id           ! identifier for current file
  integer(hid_t)                  :: dset_id           ! temporary dataset id
  integer(hid_t)                  :: dataspace         ! dataspace used to query dataset size
  integer(hsize_t)                :: dims_1(1)         ! current dimensions of rank-1 dataset
  integer(hsize_t)                :: dims_2(2)         ! current dimensions of rank-2 dataset
  integer(hsize_t)                :: dims_3(3)         ! current dimensions of rank-3 dataset
  integer(hsize_t)                :: max_dims_1(1)     ! maximum dimensions of rank-1 dataset
  integer(hsize_t)                :: max_dims_2(2)     ! maximum dimensions of rank-2 dataset
  integer(hsize_t)                :: max_dims_3(3)     ! maximum dimensions of rank-3 dataset
  integer                         :: logical_tmp       ! temporary integer used to read logicals
  
  call h5open_f(hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening HDF5 library" ; goto 9999 ; endif

  call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening HDF5 file '",filename,"'" ; goto 9998 ; endif
""")

# finalizing code of loading routine
def fortran_endLoader(f):
    f.write("""
9998 continue
  
  call h5fclose_f(file_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing HDF5 file '",filename,"'" ; ierr = hdfier ; endif

9999 continue

  call h5close_f(hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing HDF5 library" ; ierr = hdfier ; endif 
    
end subroutine loadSpec
""")

# write demo code
def fortran_demoLoader(f):
    f.write("""
program test_read_spec
  use read_spec
  implicit none
  type(SpecOutput) :: s
  character(*), parameter :: filename = "/home/jonathan/Uni/04_PhD/00_programs/SPEC/SPEC/InputFiles/TestCases/G3V02L1Fi.001.h5"
  
  write(*,*) "reading '",filename,"'..."
  call loadSpec(s, filename)
  write(*,*) "done"
  
  write(*,"(A,F4.2)") "SPEC version: ", s%version
  write(*,"(A,99I2)") "Lrad:", s%input%physics%Lrad
  
  call freeSpec(s)
end program test_read_spec
""")

# read a scalar (int or double) from HDF5 variable srcPath into the source code variable targetPath
def fortran_loadItem(f, item):
    
    srcName    = item.getFullName()
    
    targetName = "s"+srcName.replace("/","%")
    if item.rank>0:
        targetName += "("
        if item.indexMapping is not None:
            for dim,idxRange in enumerate(item.indexMapping):
                if dim==0:
                    targetName += idxRange
                else:
                    targetName += ", "+idxRange
        else:
            for dim in range(item.rank):
                if dim==0:
                    targetName += "1:dims_"+str(item.rank)+"(1)"
                else:
                    targetName += ", 1:dims_"+str(item.rank)+"("+str(dim+1)+")"
        targetName += ")"
    
    #print("read {} into {}".format(srcName, targetName))
    
    # translate dtype into HDF5 type
    h5type='ERROR'
    if item.dtype=='double':
        h5type='H5T_NATIVE_DOUBLE'
    elif item.dtype=='int' or item.dtype=='boolean':
        h5type='H5T_NATIVE_INTEGER'
    else:
        h5type='TYPE('+item.dtype.upper()+')'
    
        
    
    
    if item.rank==0:
        if (item.dtype=='boolean'):
            fmt="""
! {srcName} --> {targetName}; rank={rank}; h5type={h5type}
  call h5dopen_f(file_id, "{srcName}", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '{srcName}'" ; goto 9998 ; endif
  call h5dread_f(dset_id, {h5type}, logical_tmp, int((/1/), HSIZE_T), hdfier)
  {targetName} = merge(.TRUE., .FALSE., logical_tmp.ne.0)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '{srcName}'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '{srcName}'" ; goto 9998 ; endif
"""
        else:
            fmt="""
! {srcName} --> {targetName}; rank={rank}; h5type={h5type}
  call h5dopen_f(file_id, "{srcName}", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '{srcName}'" ; goto 9998 ; endif
  call h5dread_f(dset_id, {h5type}, {targetName}, int((/1/), HSIZE_T), hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '{srcName}'" ; goto 9998 ; endif
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '{srcName}'" ; goto 9998 ; endif
"""
    else:
        if (item.dtype=='boolean'):
            print("ERROR: cannot generate reader for logical array '"+srcName+"' yet!")
        fmt="""
! {srcName} --> {targetName}; rank={rank}
  call h5dopen_f(file_id, "{srcName}", dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error opening dataset '{srcName}'" ; goto 9998 ; endif
  
  ! open dataspace to get current state of dataset
  call h5dget_space_f(dset_id, dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error getting dataspace for dataset '{srcName}'" ; goto 9998 ; endif
  
  ! get current size of dataset
  call h5sget_simple_extent_dims_f(dataspace, dims_{rank}, max_dims_{rank}, hdfier)
  if (hdfier.ne.{rank}) then ; write(*,*) "unexpected rank of dataset '{srcName}': ",hdfier," .ne. {rank}" ; goto 9998 ; endif

  ! close dataspace after it has been used to query the size of the variable
  call h5sclose_f(dataspace, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataspace for dataset '{srcName}'" ; goto 9998 ; endif
  
  allocate({targetName})
  
  call h5dread_f(dset_id, {h5type}, {targetName}, dims_{rank}, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error reading dataset '{srcName}'" ; goto 9998 ; endif
  
  call h5dclose_f(dset_id, hdfier)
  if (hdfier.ne.0) then ; write(*,*) "error closing dataset '{srcName}'" ; goto 9998 ; endif
"""
    f.write(fmt.format(srcName=srcName, targetName=targetName, h5type=h5type, rank=item.rank))
    
# initial code of loading routine
def fortran_startFree(f):
    f.write("""subroutine freeSpec(s)
  implicit none
  type(SpecOutput), intent(inout) :: s ! datastructure to free
""")

# finalizing code of loading routine
def fortran_endFree(f):
    f.write("""end subroutine freeSpec
""")

# free an allocated item of rank .ge. 1
def fortran_freeItem(f, item):
    
    srcName    = item.getFullName()
    targetName = "s"+srcName.replace("/","%")
    
    if (item.rank > 0):
        print("free {}".format(targetName))
        f.write("  deallocate("+targetName+")\n")
    

#%% actually generate Fortran module for reading SPEC output files
def genFortranReader(outdir, moduleName, s):
    
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
    
    
    fortranFilename = outdir+moduleName+".f90"
    print("creating Fortran reading module into '"+fortranFilename+"'")
    
    # begin code for root group (== enclosing class)
    f=open(fortranFilename, "w")
    
    f.write("""! AUTO-GENERATED; DO NOT COMMIT CHANGES TO THIS FILE !
! """+creation_tag+"""
module """+moduleName+"\n")
    
    # custom datatypes come first
    for dtype in s.getDatatypes():
        f.write(fortran_genType(dtype.name, dtype.items)+'\n')
        
    # we need to reverse the definition order so that types which are used inside other types
    # are already defined when used
    reverse_groupStack =  []
    
    groupStack = []
    groupStack.append(s.rootGroup)
    while len(groupStack)>0:
        currentGroup = groupStack[-1]
        groupStack = groupStack[:-1]
        
        if type(currentGroup)==Group:
            reverse_groupStack.append(currentGroup)
    
        for item in currentGroup.items:
            if type(item)==Group:
                groupStack.append(item)
    
    # iterate in reverse order over the discovered variables to generate type definitions in correct order
    for currentGroup in reverse_groupStack[::-1]:
        f.write(fortran_genType(currentGroup.name, currentGroup.items)+'\n')
    
    f.write("contains\n")
    
    # initial code of loading routine
    fortran_startLoader(f)
    
    # loop over all variables again and put the loader code for each of them one after another
    for currentGroup in reverse_groupStack[::-1]:
        for item in currentGroup.items:
            if type(item)==Dataset:
                fortran_loadItem(f, item)
        
    # finalizing code of loading routine
    fortran_endLoader(f)
    
    # write the freeSpec subroutine to free the memory it occupied
    fortran_startFree(f)
    
    for currentGroup in reverse_groupStack[::-1]:
        for item in currentGroup.items:
            if type(item)==Dataset:
                fortran_freeItem(f, item)
    
    # finalizing code of freeing routine
    fortran_endFree(f)
    
    f.write("end module read_spec\n")

    # write demo code
    #fortran_demoLoader(f)

    f.close()