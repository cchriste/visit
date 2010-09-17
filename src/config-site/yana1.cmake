#/usr/gapps/visit/cmake/2.8.0/linux-x86_64_gcc-4.1/bin/cmake
#
# Modifications:
#     Mark C. Miller, Fri Jul 16 08:57:10 PDT 2010
#     Updated to silo-4.8-pre3 for testing.
#
#     Mark C. Miller, Tue Jul 27 17:48:54 PDT 2010
#     Update to latest ITAPS implementation. Fixed link from
#     ITAPS_MOAB to 3.6.0 and not 3.6.3 of netcdf.
#
#     Mark C. Miller, Thu Sep 16 17:57:22 PDT 2010
#     Updated Silo to 4.8

##
## Set the VISITHOME environment variable.
##
SET(VISITHOME /usr/gapps/visit)
SET(VISITARCH linux-x86_64_gcc-4.1.2)
SET(VISITARCH2 linux-x86_64_gcc-4.1)

SET(VISIT_VERBOSE_MAKEFILE TRUE)

##
## MESA
##
VISIT_OPTION_DEFAULT(VISIT_MESA_DIR ${VISITHOME}/mesa/7.8.2/${VISITARCH2})

##
## VTK
##
VISIT_OPTION_DEFAULT(VISIT_VTK_DIR ${VISITHOME}/vtk/5.0.0g/${VISITARCH2}/lib/vtk-5.0)

##
## Qt
##
# Add Qt's qmake to the path so we can use CMake's Qt4 autodetection.
VISIT_OPTION_DEFAULT(VISIT_QT_BIN ${VISITHOME}/qt/4.6.1/${VISITARCH2}/bin)

##
## Use VisIt's Python.
##
VISIT_OPTION_DEFAULT(VISIT_PYTHON_DIR ${VISITHOME}/python/2.6.4/${VISITARCH2})

##
## Ice-T
##
VISIT_OPTION_DEFAULT(VISIT_ICET_DIR ${VISITHOME}/icet/1.0.0/${VISITARCH})

##
## Add parallel arguments.
##
VISIT_OPTION_DEFAULT(VISIT_PARALLEL ON)
VISIT_OPTION_DEFAULT(VISIT_MPI_CXX_FLAGS -I/usr/local/tools/mvapich-gnu/include)
VISIT_OPTION_DEFAULT(VISIT_MPI_C_FLAGS   -I/usr/local/tools/mvapich-gnu/include)
VISIT_OPTION_DEFAULT(VISIT_MPI_LD_FLAGS  "-L/usr/local/tools/mvapich-gnu/lib/shared -L/usr/local/tools/mvapich-gnu/lib -Wl,-rpath=/usr/local/tools/mvapich-gnu/lib/shared")
VISIT_OPTION_DEFAULT(VISIT_MPI_LIBS     mpich)
VISIT_OPTION_DEFAULT(VISIT_PARALLEL_RPATH  "/usr/local/tools/mvapich-gnu/lib/shared")

##
## Turn off warnings for deprecated features on g++
##
VISIT_OPTION_DEFAULT(VISIT_CXX_FLAGS "-Wno-deprecated -fPIC -fvisibility=hidden")
VISIT_OPTION_DEFAULT(VISIT_C_FLAGS "-Wno-deprecated -fPIC -fvisibility=hidden")

##
## Database reader plugin support libraries
##
###############################################################################
 
##
## Boxlib
##
VISIT_OPTION_DEFAULT(VISIT_BOXLIB2D_DIR ${VISITHOME}/boxlib/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_BOXLIB3D_DIR ${VISITHOME}/boxlib/${VISITARCH})

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
VISIT_OPTION_DEFAULT(VISIT_CGNS_DIR ${VISITHOME}/cgns/2.4/${VISITARCH})

##
## Exodus
##
VISIT_OPTION_DEFAULT(VISIT_EXODUSII_DIR ${VISITHOME}/exodus/4.46/${VISITARCH})

##
## GDAL
##
VISIT_OPTION_DEFAULT(VISIT_GDAL_DIR ${VISITHOME}/gdal/1.3.2/${VISITARCH})

##
## H5PART
##
VISIT_OPTION_DEFAULT(VISIT_H5PART_DIR ${VISITHOME}/h5part/1.6.0/${VISITARCH})

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
## ITAPS
##
## MOAB implementation
ITAPS_INCLUDE_DIRECTORIES(MOAB ${VISITHOME}/itaps/MOAB/3.99-10Jun10/${VISITARCH}/include)
ITAPS_FILE_PATTERNS(MOAB *.cub)
ITAPS_LINK_LIBRARIES(MOAB iMesh MOAB hdf5 sz z netcdf_c++ netcdf vtkGraphics)
ITAPS_LINK_DIRECTORIES(MOAB  ${VISITHOME}/itaps/MOAB/3.99-10Jun10/${VISITARCH}/lib  ${VISITHOME}/hdf5/1.8.4/${VISITARCH}/lib  ${VISITHOME}/szip/2.1/${VISITARCH}/lib  ${VISITHOME}/netcdf/3.6.0/${VISITARCH}/lib)
## FMDB implementation
ITAPS_INCLUDE_DIRECTORIES(FMDB ${VISITHOME}/itaps/FMDB/1.1-30Aug10/${VISITARCH}/include/FMDB)
ITAPS_FILE_PATTERNS(FMDB *.sms)
ITAPS_LINK_LIBRARIES(FMDB FMDB SCORECModel SCORECUtil vtkGraphics)
ITAPS_LINK_DIRECTORIES(FMDB ${VISITHOME}/itaps/FMDB/1.1-30Aug10/${VISITARCH}/lib)
## GRUMMP implementation
ITAPS_INCLUDE_DIRECTORIES(GRUMMP ${VISITHOME}/itaps/GRUMMP/0.5.0-03Jun10/${VISITARCH}/include)
ITAPS_FILE_PATTERNS(GRUMMP *.bdry *.smesh *.vmesh)
ITAPS_LINK_LIBRARIES(GRUMMP iMesh_GRUMMP GR_3D GR_surf GR_2D GR_base SUMAAlog_lite OptMS vtkGraphics)
ITAPS_LINK_DIRECTORIES(GRUMMP ${VISITHOME}/itaps/GRUMMP/0.5.0-03Jun10/${VISITARCH}/lib)

##
## Mili
##
VISIT_OPTION_DEFAULT(VISIT_MILI_DIR ${VISITHOME}/mili/1.10.0/${VISITARCH})

##
## netCDF
##
VISIT_OPTION_DEFAULT(VISIT_NETCDF_DIR ${VISITHOME}/netcdf/3.6.0/${VISITARCH})

##
## Silo
##
VISIT_OPTION_DEFAULT(VISIT_SILO_DIR ${VISITHOME}/silo/4.8/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_SILO_LIBDEP HDF5_LIBRARY_DIR hdf5 ${VISIT_HDF5_LIBDEP})

##
## Xdmf
##
VISIT_OPTION_DEFAULT(VISIT_XDMF_DIR ${VISITHOME}/Xdmf/2.1/${VISITARCH2})
VISIT_OPTION_DEFAULT(VISIT_XDMF_LIBDEP HDF5_LIBRARY_DIR hdf5 ${VISIT_HDF5_LIBDEP})

##
## Visus
##
VISIT_OPTION_DEFAULT(VISIT_VISUS_DIR ${VISITHOME}/visus/26Feb07/${VISITARCH})

