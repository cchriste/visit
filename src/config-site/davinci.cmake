#/project/projectdirs/visit/visit_3rdparty/cmake/2.8.0/linux-ia64_gcc-4.1/bin/cmake

##
## Modifications:
##   Mark C. Miller, Mon Apr 12 18:24:36 PDT 2010
##   Changed to use mpich-1.2.7p1 and an IceT compiled for mpich-1.2.7p1.
##

##
## Set the VISITHOME environment variable.
##
SET(VISITHOME /project/projectdirs/visit/visit_3rdparty)
SET(VISITARCH linux-ia64_gcc-4.1.2)
SET(VISITARCH2 linux-ia64_gcc-4.1)

SET(VISIT_VERBOSE_MAKEFILE TRUE)

##
## MESA
##
VISIT_OPTION_DEFAULT(VISIT_MESA_DIR ${VISITHOME}/mesa/7.5/${VISITARCH2})

##
## VTK
##
VISIT_OPTION_DEFAULT(VISIT_VTK_DIR ${VISITHOME}/vtk/5.0.0e/${VISITARCH2}/lib/vtk-5.0)

##
## Qt
##
# Add Qt's qmake to the path so we can use CMake's Qt4 autodetection.
VISIT_OPTION_DEFAULT(VISIT_QT_BIN ${VISITHOME}/qt/4.6.1/${VISITARCH2}/bin)

##
## Python
##
VISIT_OPTION_DEFAULT(VISIT_PYTHON_DIR ${VISITHOME}/python/2.6.4/${VISITARCH2})

##
## Ice-T
##
VISIT_OPTION_DEFAULT(VISIT_ICET_DIR ${VISITHOME}/icet/1.0.0-mpich-1.2.7p1/${VISITARCH})

##
## Turn off warnings for deprecated features and generate position independent code.
##
VISIT_OPTION_DEFAULT(VISIT_CXX_FLAGS "-Wno-deprecated -fPIC -fvisibility=hidden")
VISIT_OPTION_DEFAULT(VISIT_C_FLAGS "-fPIC -fvisibility=hidden")

##
## Add parallel arguments.
##
VISIT_OPTION_DEFAULT(VISIT_PARALLEL ON)
VISIT_OPTION_DEFAULT(VISIT_MPI_LIBRARY_DIR ${VISITHOME}/mpich/1.2.7p1/${VISITARCH}/lib)
VISIT_OPTION_DEFAULT(VISIT_MPI_CXX_FLAGS "-I${VISITHOME}/mpich/1.2.7p1/${VISITARCH}/include")
VISIT_OPTION_DEFAULT(VISIT_MPI_LD_FLAGS "-L${VISIT_MPI_LIBRARY_DIR}")
VISIT_OPTION_DEFAULT(VISIT_MPI_LIBS mpich)

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
## ExodusII
##
VISIT_OPTION_DEFAULT(VISIT_EXODUSII_DIR ${VISITHOME}/exodus/4.46/${VISITARCH})

##
## GDAL
##
VISIT_OPTION_DEFAULT(VISIT_GDAL_DIR ${VISITHOME}/gdal/1.3.2/${VISITARCH})

##
## FastBit
##
VISIT_OPTION_DEFAULT(VISIT_FASTBIT_DIR ${VISITHOME}/fastbit/1.0.8c/${VISITARCH})

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
## H5PART
##
VISIT_OPTION_DEFAULT(VISIT_H5PART_DIR ${VISITHOME}/h5part/1.6.0/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_H5PART_LIBDEP HDF5_LIBRARY_DIR hdf5 ${VISIT_HDF5_LIBDEP})

##
## ITAPS
##
# MOAB implementation
ITAPS_INCLUDE_DIRECTORIES(MOAB ${VISITHOME}/itaps/MOAB/3.99-10Jun10/linux-ia64_gcc-4.1/include)
ITAPS_FILE_PATTERNS(MOAB *.cub)
ITAPS_LINK_LIBRARIES(MOAB iMesh MOAB hdf5 sz z netcdf_c++ netcdf vtkGraphics)
ITAPS_LINK_DIRECTORIES(MOAB 
    ${VISITHOME}/itaps/MOAB/3.99-10Jun10/linux-ia64_gcc-4.1/lib
    ${VISITHOME}/hdf5/1.8.4/${VISITARCH}/lib
    ${VISITHOME}/szip/2.1/${VISITARCH}/lib
    ${VISITHOME}/netcdf/3.6.0/${VISITARCH}/lib)
# FMDB implementation
ITAPS_INCLUDE_DIRECTORIES(FMDB ${VISITHOME}/itaps/FMDB/1.1-03Jun10/linux-ia64_gcc-4.1/include/FMDB)
ITAPS_FILE_PATTERNS(FMDB *.sms)
ITAPS_LINK_LIBRARIES(FMDB FMDB SCORECModel SCORECUtil vtkGraphics)
ITAPS_LINK_DIRECTORIES(FMDB ${VISITHOME}/itaps/FMDB/1.1-03Jun10/linux-ia64_gcc-4.1/lib)
# GRUMMP implementation
ITAPS_INCLUDE_DIRECTORIES(GRUMMP ${VISITHOME}/itaps/GRUMMP/0.5.0-03Jun10/linux-ia64_gcc-4.1/include)
ITAPS_FILE_PATTERNS(GRUMMP *.bdry *.smesh *.vmesh)
ITAPS_LINK_LIBRARIES(GRUMMP iMesh_GRUMMP GR_3D GR_surf GR_2D GR_base SUMAAlog_lite OptMS vtkGraphics)
ITAPS_LINK_DIRECTORIES(GRUMMP ${VISITHOME}/itaps/GRUMMP/0.5.0-03Jun10/linux-ia64_gcc-4.1/lib)

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
VISIT_OPTION_DEFAULT(VISIT_SILO_DIR ${VISITHOME}/silo/4.7.2/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_SILO_LIBDEP HDF5_LIBRARY_DIR hdf5 ${VISIT_HDF5_LIBDEP})

##
## ViSUS
##
VISIT_OPTION_DEFAULT(VISIT_VISUS_DIR ${VISITHOME}/visus/26Feb07/${VISITARCH})

##
## Xdmf
##
VISIT_OPTION_DEFAULT(VISIT_XDMF_DIR ${VISITHOME}/Xdmf/2.1/${VISITARCH2})
VISIT_OPTION_DEFAULT(VISIT_XDMF_LIBDEP HDF5_LIBRARY_DIR hdf5 ${VISIT_HDF5_LIBDEP})
