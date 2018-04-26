#/Users/steve/Research/SCI/workspace/visit/deps/visit/cmake/3.8.1/i386-apple-darwin16_clang/bin/cmake
##
## ./build_visit generated host.cmake
## created: Thu Apr 26 09:32:11 MDT 2018
## system: Darwin dhcp219.sci.utah.edu 16.7.0 Darwin Kernel Version 16.7.0: Thu Jan 11 22:59:40 PST 2018; root:xnu-3789.73.8~1/RELEASE_X86_64 x86_64
## by: steve

##
## Setup VISITHOME & VISITARCH variables.
##
SET(VISITHOME /Users/steve/Research/SCI/workspace/visit/deps/visit)
SET(VISITARCH i386-apple-darwin16_clang)

## Compiler flags.
##
VISIT_OPTION_DEFAULT(VISIT_C_COMPILER clang TYPE FILEPATH)
VISIT_OPTION_DEFAULT(VISIT_CXX_COMPILER clang++ TYPE FILEPATH)
VISIT_OPTION_DEFAULT(VISIT_FORTRAN_COMPILER no TYPE FILEPATH)
VISIT_OPTION_DEFAULT(VISIT_C_FLAGS "-fno-common -fexceptions -O2" TYPE STRING)
VISIT_OPTION_DEFAULT(VISIT_CXX_FLAGS "-fno-common -fexceptions -O2" TYPE STRING)

##
## Parallel Build Setup.
##
VISIT_OPTION_DEFAULT(VISIT_PARALLEL ON TYPE BOOL)
## (configured w/ mpi compiler wrapper)
VISIT_OPTION_DEFAULT(VISIT_MPI_COMPILER /Users/steve/Research/SCI/workspace/visit/deps/visit/mpich/3.0.4/i386-apple-darwin16_clang/bin/mpicc TYPE FILEPATH)

##
## VisIt Thread Option
##
VISIT_OPTION_DEFAULT(VISIT_THREAD OFF TYPE BOOL)

##############################################################
##
## Database reader plugin support libraries
##
## The HDF4, HDF5 and NetCDF libraries must be first so that
## their libdeps are defined for any plugins that need them.
##
## For libraries with LIBDEP settings, order matters.
## Libraries with LIBDEP settings that depend on other
## Library's LIBDEP settings must come after them.
##############################################################
##

##
## ZLIB
##
VISIT_OPTION_DEFAULT(VISIT_ZLIB_DIR ${VISITHOME}/zlib/1.2.7/${VISITARCH})

##
## Python
##
VISIT_OPTION_DEFAULT(VISIT_PYTHON_DIR ${VISITHOME}/python/2.7.11/${VISITARCH})

##
## Qt
##
SETUP_APP_VERSION(QT 5.6.1)
VISIT_OPTION_DEFAULT(VISIT_QT5 ON TYPE BOOL)
VISIT_OPTION_DEFAULT(VISIT_QT_DIR ${VISITHOME}/qt/${QT_VERSION}/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_QT_BIN ${VISIT_QT_DIR}/bin)

##
## QWT
##
SETUP_APP_VERSION(QWT 6.1.2)
VISIT_OPTION_DEFAULT(VISIT_QWT_DIR ${VISITHOME}/qwt/${QWT_VERSION}/${VISITARCH})

##
## VTK
##
SETUP_APP_VERSION(VTK 6.1.0)
VISIT_OPTION_DEFAULT(VISIT_VTK_DIR ${VISITHOME}/vtk/${VTK_VERSION}/${VISITARCH})

##
## MPICH
##
# Give VisIt information so it can install MPI into the binary distribution.
VISIT_OPTION_DEFAULT(VISIT_MPICH_DIR ${VISITHOME}/mpich/3.0.4/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_MPICH_INSTALL ON TYPE BOOL)

# Tell VisIt the parallel compiler so it can deduce parallel flags
VISIT_OPTION_DEFAULT(VISIT_MPI_COMPILER ${VISIT_MPICH_DIR}/bin/mpicc TYPE FILEPATH)
VISIT_OPTION_DEFAULT(VISIT_PARALLEL ON TYPE BOOL)

##
## Ice-T
##
VISIT_OPTION_DEFAULT(VISIT_ICET_DIR ${VISITHOME}/icet/1.0.0/${VISITARCH})
##

##
## Uintah
##
SETUP_APP_VERSION(UINTAH 2.5.0)
VISIT_OPTION_DEFAULT(VISIT_UINTAH_DIR ${VISITHOME}/uintah/${UINTAH_VERSION}/${VISITARCH})

