#/apps/visit/trunk/thirdparty/visit/cmake/3.8.1/linux-x86_64_gcc-5.4/bin/cmake
##
## ./build_visit generated host.cmake
## created: Wed Apr 11 13:39:22 EDT 2018
## system: Linux whoopingcough 4.4.0-119-generic #143-Ubuntu SMP Mon Apr 2 16:08:24 UTC 2018 x86_64 x86_64 x86_64 GNU/Linux
## by: pugmire

##
## Setup VISITHOME & VISITARCH variables.
##
SET(VISITHOME /apps/visit/trunk/thirdparty/visit)
SET(VISITARCH linux-x86_64_gcc-5.4)

## Compiler flags.
##
VISIT_OPTION_DEFAULT(VISIT_C_COMPILER gcc TYPE FILEPATH)
VISIT_OPTION_DEFAULT(VISIT_CXX_COMPILER g++ TYPE FILEPATH)
VISIT_OPTION_DEFAULT(VISIT_FORTRAN_COMPILER no TYPE FILEPATH)
VISIT_OPTION_DEFAULT(VISIT_C_FLAGS "-m64 -fPIC  -O2 -fvisibility=hidden" TYPE STRING)
VISIT_OPTION_DEFAULT(VISIT_CXX_FLAGS "-std=c++11 -m64 -fPIC  -O2 -fvisibility=hidden" TYPE STRING)

##
## Parallel Build Setup.
##
#VISIT_OPTION_DEFAULT(VISIT_PARALLEL ON TYPE BOOL)
## (configured w/ mpi compiler wrapper)
#VISIT_OPTION_DEFAULT(VISIT_MPI_COMPILER /usr/bin/mpicxx TYPE FILEPATH)

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
VISIT_OPTION_DEFAULT(VISIT_PYTHON_DIR /usr)
VISIT_OPTION_DEFAULT(PYTHON_INCLUDE_PATH /usr/include/python2.7 )
VISIT_OPTION_DEFAULT(PYTHON_LIBRARY /usr/lib/x86_64-linux-gnu/libpython2.7.so)
VISIT_OPTION_DEFAULT(PYTHON_LIBRARY_DIR /usr/lib/x86_64-linux-gnu)
VISIT_OPTION_DEFAULT(PYTHON_VERSION 2.7)
SET(VISIT_PYTHON_SKIP_INSTALL ON)

##
## Qt
##
SETUP_APP_VERSION(QT 5.6.1)
VISIT_OPTION_DEFAULT(VISIT_QT5 ON TYPE BOOL)
VISIT_OPTION_DEFAULT(VISIT_QT_DIR ${VISITHOME}/qt/${QT_VERSION}/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_QT_BIN ${VISIT_QT_DIR}/bin)

## system qt
#SETUP_APP_VERSION(QT 5.5.1)
#VISIT_OPTION_DEFAULT(QT_QTUITOOLS_INCLUDE_DIR /usr/include/x86_64-linux-gnu/qt5/QtUiTools)
#VISIT_OPTION_DEFAULT(VISIT_QT_BIN /usr/lib/x86_64-linux-gnu/qt5/bin)
#SET(VISIT_QT_SKIP_INSTALL ON)

##
## QWT
##
SETUP_APP_VERSION(QWT 6.1.2)
VISIT_OPTION_DEFAULT(VISIT_QWT_DIR ${VISITHOME}/qwt/${QWT_VERSION}/${VISITARCH})

##
## VTK
##
SETUP_APP_VERSION(VTK 8.1.0)
VISIT_OPTION_DEFAULT(VISIT_VTK_DIR ${VISITHOME}/vtk/${VTK_VERSION}/${VISITARCH})

##
## SZIP
##
VISIT_OPTION_DEFAULT(VISIT_SZIP_DIR ${VISITHOME}/szip/2.1/${VISITARCH})

##
## HDF5
##
VISIT_OPTION_DEFAULT(VISIT_HDF5_DIR ${VISITHOME}/hdf5/1.8.14/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_HDF5_LIBDEP ${VISITHOME}/szip/2.1/${VISITARCH}/lib sz ${VISITHOME}/zlib/1.2.7/${VISITARCH}/lib z TYPE STRING)

##
## ADIOS
## (configured w/ mpi compiler wrapper)
##
#SETUP_APP_VERSION(ADIOS 1.11.0)
#VISIT_OPTION_DEFAULT(VISIT_ADIOS_DIR ${VISITHOME}/adios/${ADIOS_VERSION}/${VISITARCH})

#SETUP_APP_VERSION(ADIOS 1.11.0)
VISIT_OPTION_DEFAULT(VISIT_ADIOS2_DIR /apps/adios2/install)

##
## BOOST
##
SETUP_APP_VERSION(BOOST 1_60_0)
VISIT_OPTION_DEFAULT(VISIT_BOOST_DIR ${VISITHOME}/boost/${BOOST_VERSION}/${VISITARCH})

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
## CGNS
##
VISIT_OPTION_DEFAULT(VISIT_CGNS_DIR ${VISITHOME}/cgns/3.2.1/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_CGNS_LIBDEP HDF5_LIBRARY_DIR hdf5 ${VISIT_HDF5_LIBDEP} TYPE STRING)

##
## GDAL
##
#VISIT_OPTION_DEFAULT(VISIT_GDAL_DIR ${VISITHOME}/gdal/2.2.4/${VISITARCH})

##
## H5Part
##
SETUP_APP_VERSION(H5PART 1.6.6)
VISIT_OPTION_DEFAULT(VISIT_H5PART_DIR ${VISITHOME}/h5part/${H5PART_VERSION}/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_H5PART_LIBDEP HDF5_LIBRARY_DIR hdf5 ${VISIT_HDF5_LIBDEP} TYPE STRING)

##
## Ice-T
##
VISIT_OPTION_DEFAULT(VISIT_ICET_DIR ${VISITHOME}/icet/77c708f9090236b576669b74c53e9f105eedbd7e/${VISITARCH})

##
## MFEM
##
VISIT_OPTION_DEFAULT(VISIT_MFEM_DIR ${VISITHOME}/mfem/3.3/${VISITARCH})
##

##
## Nektar++
##
#SETUP_APP_VERSION(NEKTAR++ 4.4.0)
#VISIT_OPTION_DEFAULT(VISIT_NEKTAR++_DIR ${VISITHOME}/nektar++/${NEKTAR++_VERSION}/${VISITARCH})
#VISIT_OPTION_DEFAULT(VISIT_NEKTAR++_LIBDEP ${VISITHOME}/zlib/1.2.7/${VISITARCH}/lib z TYPE STRING)

##
## NetCDF
##
VISIT_OPTION_DEFAULT(VISIT_NETCDF_DIR ${VISITHOME}/netcdf/4.1.1/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_NETCDF_LIBDEP HDF5_LIBRARY_DIR hdf5_hl HDF5_LIBRARY_DIR hdf5 ${VISIT_HDF5_LIBDEP} TYPE STRING)

##
## OpenEXR
##
VISIT_OPTION_DEFAULT(VISIT_OPENEXR_DIR ${VISITHOME}/openexr/2.2.0/${VISITARCH})

##
## PIDX
##
#SETUP_APP_VERSION(PIDX 0.9.1)
#VISIT_OPTION_DEFAULT(VISIT_PIDX_DIR ${VISITHOME}/pidx/${PIDX_VERSION}/${VISITARCH})

##
## Silo
##
VISIT_OPTION_DEFAULT(VISIT_SILO_DIR ${VISITHOME}/silo/4.10.1/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_SILO_LIBDEP HDF5_LIBRARY_DIR hdf5 ${VISIT_HDF5_LIBDEP} ZLIB_LIBRARY_DIR z TYPE STRING)

##
## VISUS
##
VISIT_OPTION_DEFAULT(VISIT_VISUS_DIR ${VISITHOME}/visus/5f5fd6c/${VISITARCH})

##
## Xdmf
##
VISIT_OPTION_DEFAULT(VISIT_XDMF_DIR ${VISITHOME}/Xdmf/2.1.1/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_XDMF_LIBDEP HDF5_LIBRARY_DIR hdf5  VTK_LIBRARY_DIRS vtklibxml2-${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}  TYPE STRING)

##
## ADIOS2
##
SETUP_APP_VERSION(ADIOS2 0.1.0)
VISIT_OPTION_DEFAULT(VISIT_ADIOS2_DIR ${VISITHOME}/adios2/${ADIOS2_VERSION}/${VISITARCH})
