#/projects/visit/thirdparty_shared_2_6_0/cmake/2.8.8/linux-x86_64_gcc-4.3/bin/cmake
##
## ./build_visit2_6_0 generated host.cmake
## created: Thu Nov 15 10:26:33 MST 2012
## system: Linux redsky-login3 2.6.18-108redsky_chaos #1 SMP Mon Dec 5 11:18:51 MST 2011 x86_64 x86_64 x86_64 GNU/Linux
## by: ebrugge

##
## Setup VISITHOME & VISITARCH variables.
##
SET(VISITHOME /projects/visit/thirdparty_shared_2_6_0)
SET(VISITARCH linux-x86_64_gcc-4.3)

## Compiler flags.
##
VISIT_OPTION_DEFAULT(VISIT_C_COMPILER gcc TYPE FILEPATH)
VISIT_OPTION_DEFAULT(VISIT_CXX_COMPILER g++ TYPE FILEPATH)
VISIT_OPTION_DEFAULT(VISIT_C_FLAGS " -m64 -fvisibility=hidden" TYPE STRING)
VISIT_OPTION_DEFAULT(VISIT_CXX_FLAGS " -m64 -fvisibility=hidden" TYPE STRING)

##
## Add parallel arguments.
##
VISIT_OPTION_DEFAULT(VISIT_PARALLEL ON)
VISIT_OPTION_DEFAULT(VISIT_MPI_CXX_FLAGS -I/apps/x86_64/mpi/openmpi/gcc-4.3.4/openmpi-1.4.3_oobpr/include)
VISIT_OPTION_DEFAULT(VISIT_MPI_C_FLAGS   -I/apps/x86_64/mpi/openmpi/gcc-4.3.4/openmpi-1.4.3_oobpr/include)
VISIT_OPTION_DEFAULT(VISIT_MPI_LD_FLAGS  "-L/apps/x86_64/mpi/openmpi/gcc-4.3.4/openmpi-1.4.3_oobpr/lib -Wl,-rpath=/apps/x86_64/mpi/openmpi/gcc-4.3.4/openmpi-1.4.3_oobpr/lib")
VISIT_OPTION_DEFAULT(VISIT_MPI_LIBS     mpi_cxx mpi)
VISIT_OPTION_DEFAULT(VISIT_PARALLEL_RPATH  "/apps/x86_64/mpi/openmpi/gcc-4.3.4/openmpi-1.4.3_oobpr/lib")

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
## Python
##
VISIT_OPTION_DEFAULT(VISIT_PYTHON_DIR /projects/visit/thirdparty_shared_2_6_0/python/2.6.4/linux-x86_64_gcc-4.3)

##
## Mesa
##
VISIT_OPTION_DEFAULT(VISIT_MESA_DIR ${VISITHOME}/mesa/7.8.2/${VISITARCH})

##
## VTK
##
VISIT_OPTION_DEFAULT(VISIT_VTK_DIR ${VISITHOME}/vtk/5.8.0.a/${VISITARCH})

##
## Qt
##
VISIT_OPTION_DEFAULT(VISIT_QT_BIN ${VISITHOME}/qt/4.8.3/${VISITARCH}/bin)

##
## SZIP
##
VISIT_OPTION_DEFAULT(VISIT_SZIP_DIR ${VISITHOME}/szip/2.1/${VISITARCH})

##
## HDF5
##
VISIT_OPTION_DEFAULT(VISIT_HDF5_DIR ${VISITHOME}/hdf5/1.8.7/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_HDF5_LIBDEP ${VISITHOME}/szip/2.1/${VISITARCH}/lib sz /usr/lib z TYPE STRING)

##
## Mili
##
VISIT_OPTION_DEFAULT(VISIT_MILI_DIR ${VISITHOME}/mili/1.11.1/${VISITARCH})

##
## Ice-T
##
VISIT_OPTION_DEFAULT(VISIT_ICET_DIR ${VISITHOME}/icet/1.0.0/${VISITARCH})

##
## HDF4
##
VISIT_OPTION_DEFAULT(VISIT_HDF4_DIR ${VISITHOME}/hdf4/4.2.5/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_HDF4_LIBDEP ${VISITHOME}/szip/2.1/${VISITARCH}/lib sz ${VISITHOME}/vtk/5.8.0.a/${VISITARCH}/lib vtkjpeg TYPE STRING)

##
## NetCDF
##
VISIT_OPTION_DEFAULT(VISIT_NETCDF_DIR ${VISITHOME}/netcdf/4.1.1/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_NETCDF_LIBDEP HDF5_LIBRARY_DIR hdf5_hl HDF5_LIBRARY_DIR hdf5 ${VISIT_HDF5_LIBDEP} TYPE STRING)

##
## CGNS
##
VISIT_OPTION_DEFAULT(VISIT_CGNS_DIR ${VISITHOME}/cgns/3.0.8/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_CGNS_LIBDEP HDF5_LIBRARY_DIR hdf5 ${VISIT_HDF5_LIBDEP} TYPE STRING)

##
## GDAL
##
VISIT_OPTION_DEFAULT(VISIT_GDAL_DIR ${VISITHOME}/gdal/1.7.1/${VISITARCH})

##
## Boxlib
##
VISIT_OPTION_DEFAULT(VISIT_BOXLIB_DIR ${VISITHOME}/boxlib/1.0.7/${VISITARCH})

##
## CFITSIO
##
VISIT_OPTION_DEFAULT(VISIT_CFITSIO_DIR ${VISITHOME}/cfitsio/3006/${VISITARCH})

##
## H5Part
##
VISIT_OPTION_DEFAULT(VISIT_H5PART_DIR ${VISITHOME}/h5part/1.6.6/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_H5PART_LIBDEP HDF5_LIBRARY_DIR hdf5 ${VISIT_HDF5_LIBDEP} TYPE STRING)

##
## FastBit
##
VISIT_OPTION_DEFAULT(VISIT_FASTBIT_DIR ${VISITHOME}/fastbit/1.2.0/${VISITARCH})

##
## CCMIO
##
VISIT_OPTION_DEFAULT(VISIT_CCMIO_DIR ${VISITHOME}/ccmio/2.6.1/${VISITARCH})

##
## Silo
##
VISIT_OPTION_DEFAULT(VISIT_SILO_DIR ${VISITHOME}/silo/4.9/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_SILO_LIBDEP HDF5_LIBRARY_DIR hdf5 ${VISIT_HDF5_LIBDEP} TYPE STRING)

##
## AdvIO
##
VISIT_OPTION_DEFAULT(VISIT_ADVIO_DIR ${VISITHOME}/AdvIO/1.2/${VISITARCH}/)

##
## Xdmf
##
VISIT_OPTION_DEFAULT(VISIT_XDMF_DIR ${VISITHOME}/Xdmf/2.1.1/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_XDMF_LIBDEP HDF5_LIBRARY_DIR hdf5   VTK_LIBRARY_DIRS vtklibxml2  TYPE STRING)

##
## PySide
##
VISIT_OPTION_DEFAULT(VISIT_PYSIDE_DIR ${VISITHOME}/pyside/1.1.1/${VISITARCH}/)

