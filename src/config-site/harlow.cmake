#/scratch/users/gweber/visit_3rdparty/cmake/2.8.0/linux-x86_64_gcc-4.3/bin/cmake
##
## ../visit_svn/trunk/src/svn_bin/build_visit generated host.cmake
## created: Thu Jan 21 17:56:19 PST 2010
## system: Linux harlow 2.6.27-16-generic #1 SMP Tue Dec 1 19:26:23 UTC 2009 x86_64 GNU/Linux
## by: gweber

##
## Setup VISITHOME & VISITARCH variables.
##
SET(VISITHOME /scratch/users/gweber/visit_3rdparty)
SET(VISITARCH linux-x86_64_gcc-4.3)
SET(VISIT_VERBOSE_MAKEFILE TRUE)
VISIT_OPTION_DEFAULT(CMAKE_INSTALL_PREFIX /usr/local/anag/pkg/VisIt-2.0.0b)
VISIT_OPTION_DEFAULT(VISIT_INSTALL_THIRD_PARTY ON)
VISIT_OPTION_DEFAULT(CMAKE_BUILD_TYPE Release)
VISIT_OPTION_DEFAULT(VISIT_RPATH_RELATIVE_TO_EXECUTABLE_PATH ON)

# FIXME: Python filters break static MPI. Long term, install shared library MPI version.
VISIT_OPTION_DEFAULT(VISIT_PYTHON_FILTERS OFF)

##
## Specify the location of the mesa.
##
VISIT_OPTION_DEFAULT(VISIT_MESA_DIR ${VISITHOME}/mesa/7.5/${VISITARCH})

##
## Specify the location of the vtk.
##
VISIT_OPTION_DEFAULT(VISIT_VTK_DIR ${VISITHOME}/vtk/5.0.0d/${VISITARCH}/lib/vtk-5.0/)

##
## Specify the Qt4 binary dir. 
## (qmake us used to locate & setup Qt4 dependencies)
##
VISIT_OPTION_DEFAULT(VISIT_QT_BIN ${VISITHOME}/qt/4.6.1/${VISITARCH}/bin)

##
## Specify the location of the python.
##
VISIT_OPTION_DEFAULT(VISIT_PYTHON_DIR ${VISITHOME}/python/2.6.4/${VISITARCH})

##
## Compiler flags.
##
VISIT_OPTION_DEFAULT(VISIT_C_COMPILER gcc)
VISIT_OPTION_DEFAULT(VISIT_CXX_COMPILER g++)
VISIT_OPTION_DEFAULT(VISIT_C_FLAGS "-m64 -fPIC -fvisibility=hidden")
VISIT_OPTION_DEFAULT(VISIT_MPI_C_FLAGS "-I/usr/local/anag/pkg/mpich-1.2.5.2/include -DUSE_STDARG -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_UNISTD_H=1 -DHAVE_STDARG_H=1 -DUSE_STDARG=1 -DMALLOC_RET_VOID=1")
VISIT_OPTION_DEFAULT(VISIT_CXX_FLAGS "-m64 -fPIC -fvisibility=hidden")
VISIT_OPTION_DEFAULT(VISIT_MPI_CXX_FLAGS "-I/usr/local/anag/pkg/mpich-1.2.5.2/include -DUSE_STDARG -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_UNISTD_H=1 -DHAVE_STDARG_H=1 -DUSE_STDARG=1 -DMALLOC_RET_VOID=1")
VISIT_OPTION_DEFAULT(VISIT_MPI_LIBS "-L/usr/local/anag/pkg/mpich-1.2.5.2/lib -lmpich")
VISIT_OPTION_DEFAULT(VISIT_PARALLEL ON)
VISIT_OPTION_DEFAULT(VISIT_NOLINK_MPI_WITH_LIBRARIES ON)


##############################################################
##
## Database reader plugin support libraries
##
##############################################################

##
## Boxlib
##
VISIT_OPTION_DEFAULT(VISIT_BOXLIB2D_DIR ${VISITHOME}/boxlib/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_BOXLIB3D_DIR ${VISITHOME}/boxlib/${VISITARCH})

##
## H5Part
##
VISIT_OPTION_DEFAULT(VISIT_H5PART_DIR ${VISITHOME}/h5part/1.6.0/${VISITARCH})

##
## FastBit
##
VISIT_OPTION_DEFAULT(VISIT_FASTBIT_DIR ${VISITHOME}/fastbit/1.0.9/${VISITARCH})

##
## HDF4
##
VISIT_OPTION_DEFAULT(VISIT_HDF4_DIR ${VISITHOME}/hdf4/4.2.1/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_HDF4_LIBDEP ${VISITHOME}/szip/2.1/${VISITARCH}/lib sz /usr/lib jpeg)

##
## HDF5
##
VISIT_OPTION_DEFAULT(VISIT_HDF5_DIR ${VISITHOME}/hdf5/1.8.4/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_HDF5_LIBDEP ${VISITHOME}/szip/2.1/${VISITARCH}/lib sz)

##
## NetCDF
##
VISIT_OPTION_DEFAULT(VISIT_NETCDF_DIR ${VISITHOME}/netcdf/3.6.3/${VISITARCH})

##
## SZIP
##
VISIT_OPTION_DEFAULT(VISIT_SZIP_DIR ${VISITHOME}/szip/2.1/${VISITARCH})

##
## Silo
##
VISIT_OPTION_DEFAULT(VISIT_SILO_DIR ${VISITHOME}/silo/4.7.2/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_SILO_LIBDEP HDF5_LIBRARY_DIR hdf5 ${VISIT_HDF5_LIBDEP})
