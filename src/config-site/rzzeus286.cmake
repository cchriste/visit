#/usr/gapps/visit/thirdparty_shared/2.8.0/cmake/2.8.12.2/linux-x86_64_gcc-4.4/bin/cmake
##
## ./build_visit2_8_0 generated host.cmake
## created: Fri Aug 15 14:17:07 PDT 2014
## system: Linux rzmerl156 2.6.32-431.23.3.1chaos.ch5.2.x86_64 #1 SMP Tue Jul 29 14:58:25 PDT 2014 x86_64 x86_64 x86_64 GNU/Linux
## by: brugger

##
## Setup VISITHOME & VISITARCH variables.
##
SET(VISITHOME /usr/gapps/visit/thirdparty_shared/2.8.0)
SET(VISITARCH linux-x86_64_gcc-4.4)
VISIT_OPTION_DEFAULT(VISIT_SLIVR TRUE)

## Compiler flags.
##
VISIT_OPTION_DEFAULT(VISIT_C_COMPILER gcc TYPE FILEPATH)
VISIT_OPTION_DEFAULT(VISIT_CXX_COMPILER g++ TYPE FILEPATH)
VISIT_OPTION_DEFAULT(VISIT_C_FLAGS " -m64 -fPIC -fvisibility=hidden" TYPE STRING)
VISIT_OPTION_DEFAULT(VISIT_CXX_FLAGS " -m64 -fPIC -fvisibility=hidden" TYPE STRING)

##
## VisIt Thread Option
##
VISIT_OPTION_DEFAULT(VISIT_THREAD OFF TYPE BOOL)

##
## Add parallel arguments.
##
VISIT_OPTION_DEFAULT(VISIT_PARALLEL ON TYPE BOOL)
VISIT_OPTION_DEFAULT(VISIT_MPI_CXX_FLAGS -I/usr/local/tools/mvapich-gnu/include TYPE STRING)
VISIT_OPTION_DEFAULT(VISIT_MPI_C_FLAGS   -I/usr/local/tools/mvapich-gnu/include TYPE STRING)
VISIT_OPTION_DEFAULT(VISIT_MPI_LD_FLAGS  "-L/usr/local/tools/mvapich-gnu/lib/shared -L/usr/local/tools/mvapich-gnu/lib -Wl,-rpath=/usr/local/tools/mvapich-gnu/lib/shared" TYPE STRING)
VISIT_OPTION_DEFAULT(VISIT_MPI_LIBS     mpich TYPE STRING)
VISIT_OPTION_DEFAULT(VISIT_PARALLEL_RPATH  "/usr/local/tools/mvapich-gnu/lib/shared")

##
## Turn on DDT support.
##
VISIT_OPTION_DEFAULT(VISIT_DDT ON TYPE BOOL)

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
## Mesa
##
VISIT_OPTION_DEFAULT(VISIT_MESA_DIR ${VISITHOME}/mesa/7.10.2/${VISITARCH})

##
## Python
##
VISIT_OPTION_DEFAULT(VISIT_PYTHON_DIR ${VISITHOME}/python/2.7.6/${VISITARCH})

##
## Qt
##
VISIT_OPTION_DEFAULT(VISIT_QT_BIN ${VISITHOME}/qt/4.8.3/${VISITARCH}/bin)

##
## VTK
##
SETUP_VTK_VERSION(6.1.0)
VISIT_OPTION_DEFAULT(VISIT_VTK_DIR ${VISITHOME}/vtk/${VTK_VERSION}/${VISITARCH})
##

##
## AdvIO
##
VISIT_OPTION_DEFAULT(VISIT_ADVIO_DIR ${VISITHOME}/AdvIO/1.2/${VISITARCH})

##
## Boxlib
##
VISIT_OPTION_DEFAULT(VISIT_BOXLIB_DIR ${VISITHOME}/boxlib/1.0.7/${VISITARCH})

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
## HDF5
##
VISIT_OPTION_DEFAULT(VISIT_HDF5_DIR ${VISITHOME}/hdf5/1.8.7/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_HDF5_LIBDEP ${VISITHOME}/szip/2.1/${VISITARCH}/lib sz /usr/lib z TYPE STRING)

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
## HDF4
##
VISIT_OPTION_DEFAULT(VISIT_HDF4_DIR ${VISITHOME}/hdf4/4.2.5/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_HDF4_LIBDEP ${VISITHOME}/szip/2.1/${VISITARCH}/lib sz ${VISITHOME}/vtk/${VTK_VERSION}/${VISITARCH}/lib vtkjpeg-${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION} TYPE STRING)

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
## ITAPS
##
## MOAB implementation
ITAPS_INCLUDE_DIRECTORIES(MOAB ${VISITHOME}/itaps/1.4/MOAB/4.6.0/${VISITARCH}/include)
ITAPS_FILE_PATTERNS(MOAB *.cub)
ITAPS_LINK_LIBRARIES(MOAB iMesh MOAB hdf5_hl hdf5 sz z netcdf_c++ netcdf )
ITAPS_LINK_DIRECTORIES(MOAB  ${VISITHOME}/itaps/1.4/MOAB/4.6.0/${VISITARCH}/lib  ${VISITHOME}/hdf5/1.8.7/${VISITARCH}/lib  ${VISITHOME}/szip/2.1/${VISITARCH}/lib  ${VISITHOME}/netcdf/4.1.1/${VISITARCH}/lib )
## FMDB implementation
ITAPS_INCLUDE_DIRECTORIES(FMDB ${VISITHOME}/itaps/1.4/FMDB/1.4.0/${VISITARCH}/include)
ITAPS_FILE_PATTERNS(FMDB *.sms)
ITAPS_LINK_LIBRARIES(FMDB FMDB SCORECModel SCORECUtil )
ITAPS_LINK_DIRECTORIES(FMDB ${VISITHOME}/itaps/1.4/FMDB/1.4.0/${VISITARCH}/lib)
## GRUMMP implementation
ITAPS_INCLUDE_DIRECTORIES(GRUMMP ${VISITHOME}/itaps/1.4/GRUMMP/0.6.5/${VISITARCH}/include)
ITAPS_FILE_PATTERNS(GRUMMP *.bdry *.smesh *.vmesh)
ITAPS_LINK_LIBRARIES(GRUMMP iMesh_GRUMMP GR_3D GR_surf GR_2D GR_geom GR_base SUMAAlog_lite OptMS cgm dl)
ITAPS_LINK_DIRECTORIES(GRUMMP ${VISITHOME}/itaps/1.4/GRUMMP/0.6.5/${VISITARCH}/lib)

##
## MFEM 
##
VISIT_OPTION_DEFAULT(VISIT_MFEM_DIR ${VISITHOME}/mfem/68e941f8fe/${VISITARCH})

##
## Mili
##
VISIT_OPTION_DEFAULT(VISIT_MILI_DIR ${VISITHOME}/mili/13.1.1-patch/${VISITARCH})

##
## PySide
##
VISIT_OPTION_DEFAULT(VISIT_PYSIDE_DIR ${VISITHOME}/pyside/1.2.1/${VISITARCH}/)

##
## Silo
##
VISIT_OPTION_DEFAULT(VISIT_SILO_DIR ${VISITHOME}/silo/4.10/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_SILO_LIBDEP HDF5_LIBRARY_DIR hdf5 ${VISIT_HDF5_LIBDEP} ${VISITHOME}/silo/4.10/${VISITARCH}/json/lib json TYPE STRING)

##
## Xdmf
##
VISIT_OPTION_DEFAULT(VISIT_XDMF_DIR ${VISITHOME}/Xdmf/2.1.1/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_XDMF_LIBDEP HDF5_LIBRARY_DIR hdf5  VTK_LIBRARY_DIRS vtklibxml2-${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION} TYPE STRING)

