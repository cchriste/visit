#/Projects/VisIt/Thirdparty/visit/cmake/3.5.2/i386-apple-darwin15_clang/bin/cmake
##
## ../trunk/src/svn_bin/build_visit generated host.cmake
## created: Fri Apr 15 16:03:11 MDT 2016
## system: Darwin lautertun.sci.utah.edu 15.3.0 Darwin Kernel Version 15.3.0: Thu Dec 10 18:40:58 PST 2015; root:xnu-3248.30.4~1/RELEASE_X86_64 x86_64
## by: allen

##
## Setup VISITHOME & VISITARCH variables.
##
SET(VISITHOME /Projects/VisIt/Thirdparty/visit)
SET(VISITARCH i386-apple-darwin15_clang)

## Compiler flags.
##
VISIT_OPTION_DEFAULT(VISIT_C_COMPILER clang TYPE FILEPATH)
VISIT_OPTION_DEFAULT(VISIT_CXX_COMPILER clang++ TYPE FILEPATH)
VISIT_OPTION_DEFAULT(VISIT_C_FLAGS "-fno-common -fexceptions" TYPE STRING)
VISIT_OPTION_DEFAULT(VISIT_CXX_FLAGS "-fno-common -fexceptions" TYPE STRING)

##
## Parallel Build Setup.
##
VISIT_OPTION_DEFAULT(VISIT_PARALLEL ON TYPE BOOL)
## (configured w/ mpi compiler wrapper)
VISIT_OPTION_DEFAULT(VISIT_MPI_COMPILER ${VISIT_MPICH_DIR}/bin/mpicc TYPE FILEPATH)

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
## Python
##
VISIT_OPTION_DEFAULT(VISIT_PYTHON_DIR ${VISITHOME}/python/2.7.11/${VISITARCH})

##
## Qt
##
#VISIT_OPTION_DEFAULT(VISIT_QT_DIR ${VISITHOME}/qt/4.8.6/${VISITARCH})
#VISIT_OPTION_DEFAULT(VISIT_QT_BIN ${VISITHOME}/qt/4.8.6/${VISITARCH}/bin)

VISIT_OPTION_DEFAULT(VISIT_QT5 ON TYPE BOOL)
VISIT_OPTION_DEFAULT(VISIT_QT_DIR ${VISITHOME}/qt/5.6.0-rc/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_QT_BIN ${VISITHOME}/qt/5.6.0-rc/${VISITARCH}/bin)

##
## VTK
##
SETUP_APP_VERSION(VTK 6.1.0)
VISIT_OPTION_DEFAULT(VISIT_VTK_DIR ${VISITHOME}/vtk/${VTK_VERSION}/${VISITARCH})
##

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
## ADIOS
## (configured w/ mpi compiler wrapper)
##
VISIT_OPTION_DEFAULT(VISIT_ADIOS_DIR ${VISITHOME}/adios/1.9.0/${VISITARCH})

##
## BOOST
##
SETUP_APP_VERSION(BOOST 1_60_0)
VISIT_OPTION_DEFAULT(VISIT_BOOST_DIR ${VISITHOME}/boost/1_60_0/${VISITARCH})

##
## Boxlib
##
VISIT_OPTION_DEFAULT(VISIT_BOXLIB_DIR ${VISITHOME}/boxlib/1.3.5/${VISITARCH})

##
## CCMIO
##
VISIT_OPTION_DEFAULT(VISIT_CCMIO_DIR ${VISITHOME}/ccmio/2.6.1/${VISITARCH})

##
## CFITSIO
##
VISIT_OPTION_DEFAULT(VISIT_CFITSIO_DIR ${VISITHOME}/cfitsio/3006/${VISITARCH})

##
## SZIP
##
VISIT_OPTION_DEFAULT(VISIT_SZIP_DIR ${VISITHOME}/szip/2.1/${VISITARCH})

##
## ZLIB
##
VISIT_OPTION_DEFAULT(VISIT_ZLIB_DIR ${VISITHOME}/zlib/1.2.7/${VISITARCH})

##
## HDF5
##
VISIT_OPTION_DEFAULT(VISIT_HDF5_DIR ${VISITHOME}/hdf5/1.8.14/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_HDF5_LIBDEP ${VISITHOME}/szip/2.1/${VISITARCH}/lib sz ${VISITHOME}/zlib/1.2.7/${VISITARCH}/lib z TYPE STRING)

##
## CGNS
##
VISIT_OPTION_DEFAULT(VISIT_CGNS_DIR ${VISITHOME}/cgns/3.2.1/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_CGNS_LIBDEP HDF5_LIBRARY_DIR hdf5 ${VISIT_HDF5_LIBDEP} TYPE STRING)

##
## FastBit
##
VISIT_OPTION_DEFAULT(VISIT_FASTBIT_DIR ${VISITHOME}/fastbit/1.2.0/${VISITARCH})

##
## GDAL
##
VISIT_OPTION_DEFAULT(VISIT_GDAL_DIR ${VISITHOME}/gdal/1.10.0/${VISITARCH})

##
## H5Part
##
VISIT_OPTION_DEFAULT(VISIT_H5PART_DIR ${VISITHOME}/h5part/1.6.6/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_H5PART_LIBDEP HDF5_LIBRARY_DIR hdf5 ${VISIT_HDF5_LIBDEP} TYPE STRING)

##
## Ice-T
##
VISIT_OPTION_DEFAULT(VISIT_ICET_DIR ${VISITHOME}/icet/1.0.0/${VISITARCH})

##
## NetCDF
##
VISIT_OPTION_DEFAULT(VISIT_NETCDF_DIR ${VISITHOME}/netcdf/4.1.1/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_NETCDF_LIBDEP HDF5_LIBRARY_DIR hdf5_hl HDF5_LIBRARY_DIR hdf5 ${VISIT_HDF5_LIBDEP} TYPE STRING)

##
## Nektar++
##
SETUP_APP_VERSION(NEKTAR++ 4.1.0)
VISIT_OPTION_DEFAULT(VISIT_NEKTAR++_DIR ${VISITHOME}/nektar++/${NEKTAR++_VERSION}/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_NEKTAR++_LIBDEP ${VISITHOME}/zlib/1.2.7/${VISITARCH}/lib z TYPE STRING)

##
## Silo
##
VISIT_OPTION_DEFAULT(VISIT_SILO_DIR ${VISITHOME}/silo/4.10.1/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_SILO_LIBDEP HDF5_LIBRARY_DIR hdf5 ${VISIT_HDF5_LIBDEP} TYPE STRING)

##
## Uintah
##
#SETUP_APP_VERSION(UINTAH 1.6.0)
#VISIT_OPTION_DEFAULT(VISIT_UINTAH_DIR ${VISITHOME}/uintah/1.6.0/${VISITARCH})

SETUP_APP_VERSION(UINTAH 1.7.0)
VISIT_OPTION_DEFAULT(VISIT_UINTAH_DIR /Projects/Uintah/trunk/opt)
SET(VISIT_USE_SYSTEM_UINTAH TRUE)

##
## VISUS 
##
VISIT_OPTION_DEFAULT(VISIT_VISUS_DIR ${VISITHOME}/visus/5f5fd6c/${VISITARCH})

##
## Xdmf
##
VISIT_OPTION_DEFAULT(VISIT_XDMF_DIR ${VISITHOME}/Xdmf/2.1.1/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_XDMF_LIBDEP HDF5_LIBRARY_DIR hdf5  VTK_LIBRARY_DIRS vtklibxml2-${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}  TYPE STRING)

