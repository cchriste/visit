###############################################################################
# Program: build_visitmodule.py
#
# Purpose: This program builds the VisIt Python module.
#
# Programmer: Brad Whitlock
# Creation: Wed Nov 22 12:41:39 PDT 2006
#
# Note: The "visitmodule" module is a self-contained (no library dependencies)
#       module that includes code for the other VisIt libraries on which it
#       would normally depend (comm, state, misc, plugins, ...). We use another
#       lightweight "visit" module to import into Python and that module
#       dynamically loads this module. This does give us the flexibility to
#       let us have a "visit" module that should not change much and allows us
#       to dynamically select at runtime which VisIt Python interface should
#       be loaded.
#
# Modifications:
#
###############################################################################
from distutils.core import setup, Extension
import string

def GetCXXFilesInDirectory(directory):
    import os
    files = os.listdir(directory)
    cxx = []
    for f in files:
       if f[-2:] == ".C":
          cxx = cxx + [directory + "/" + f]
    return cxx

def RemoveFileFromList(filelist, rm):
    ret = []
    for f in filelist:
        if f != rm:
           ret = ret + [f]
    return ret

# Get all of the VisIt sources that we need to build the VisIt module.
srcdirs = ("../common/comm",
"../common/state",
"../common/misc",
"../common/utility",
"../common/plugin",
"../common/Exceptions/Database",
"../common/Exceptions/Pipeline",
"../common/Exceptions/Plotter",
"../common/Exceptions/VisWindow",
"../common/expr",
"../common/parser",
"../components/DBAtts/MetaData",
"../components/DBAtts/SIL",
"../viewer/rpc",
"../viewer/proxy",
".")
visitmod_sources = []
for d in srcdirs:
   visitmod_sources = visitmod_sources + GetCXXFilesInDirectory(d)

# Remove some of the junk
rmlist = ("../common/comm/PipeConnection.C",
"../common/comm/testcomm.C",
"../common/comm/test1.C",
"../common/comm/win32commhelpers.C",
"../common/expr/test.C",
"../common/expr/ExprConfig.C",
"../common/parser/testparser.C",
"../common/state/TestUtil.C",
"../common/state/test1.C",
"../common/state/test2.C",
"../common/state/test3.C",
"../common/state/SimulationCommand.C",  # until I figure it out...
"./cli.C")
for rm in rmlist:
   visitmod_sources = RemoveFileFromList(visitmod_sources, rm)
#print visitmod_sources

# Determine the version
f = open("../VERSION")
lines = f.readlines()
f.close()
ver = string.replace(lines[0], "\n", "")

visitmod = Extension("visitmodule",
    language="c++",
    #extra_compile_args="-Wno-deprecated",
    include_dirs=['.', "../include", "../include/visit"],
    sources = visitmod_sources)
setup(name = "visit",
      version = ver,
      description = "This module contains functions to control VisIt.",
      ext_modules = [visitmod])

