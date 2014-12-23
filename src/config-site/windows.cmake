# this section should go in its own "prepare windows type cmake file"
IF (MSVC_VERSION)
  SET(ERRMSG "You have chosen a compiler version not currently supported.  You will need to build all dependent 3rd party libraries with this compiler, and provide the correct locations for them.")
  MESSAGE(STATUS "MSVC_VERSION=${MSVC_VERSION}")
  IF (MSVC11)
      SET(VISIT_MSVC_VERSION "MSVC2012")
  ELSEIF (MSVC10)
      SET(VISIT_MSVC_VERSION "MSVC2010")
  ELSE()
      MESSAGE(SEND_ERROR "${ERRMSG}")
  ENDIF()
ELSE (MSVC_VERSION)
  MESSAGE(SEND_ERROR "${ERRMSG}")
ENDIF(MSVC_VERSION)

# this can be left in generic_windows.cmake, but allow user to
# have their own config-site file

IF (CMAKE_CL_64)
    SET(VISIT_MSVC_VERSION "${VISIT_MSVC_VERSION}-x64")
ENDIF (CMAKE_CL_64)
MESSAGE(STATUS "Prebuilt library directory name: ${VISIT_MSVC_VERSION}")

IF(NOT DEFINED VISIT_WINDOWS_DIR)
    IF(EXISTS ${VISIT_SOURCE_DIR}/../windowsbuild)
        SET(VISIT_WINDOWS_DIR ${VISIT_SOURCE_DIR}/../windowsbuild)
    ELSE()
        message(FATAL_ERROR "Please set VISIT_WINDOWS_DIR to location of 'windowsbuild' files")
    ENDIF()
ENDIF()
get_filename_component(VISIT_WINDOWS_DIR ${VISIT_WINDOWS_DIR} ABSOLUTE)
SET(VISITHOME ${VISIT_WINDOWS_DIR}/${VISIT_MSVC_VERSION})


# Create the NSIS installer package (override default from root CMakeLists.txt
OPTION(VISIT_MAKE_NSIS_INSTALLER "Create an installer package using NSIS." ON)
SET(BOOST_INCLUDE_DIR ${VISIT_SOURCE_DIR}/third_party_builtin/bilib/boost_1_42_0 CACHE PATH "boost include")

##############################################################
##
## Required libraries
##
##############################################################


##
## VTK
##
IF(VISIT_QT5)
  VISIT_OPTION_DEFAULT(VISIT_VTK_DIR  ${VISITHOME}/vtk-qt5/${VTK_VERSION})
ELSE()
  VISIT_OPTION_DEFAULT(VISIT_VTK_DIR  ${VISITHOME}/vtk/${VTK_VERSION})
ENDIF()

##
## QT
##
if (VISIT_QT5)
  VISIT_OPTION_DEFAULT(VISIT_QT_DIR     ${VISITHOME}/Qt/5.2.1)
ELSE()
  VISIT_OPTION_DEFAULT(VISIT_QT_BIN     ${VISITHOME}/Qt/4.8.3/lib)
  VISIT_OPTION_DEFAULT(QT_DIR           ${VISITHOME}/Qt/4.8.3)
  VISIT_OPTION_DEFAULT(QT_QMAKE_EXECUTABLE ${VISIT_QT_BIN}/qmake.exe TYPE FILEPATH)
ENDIF()

##
## PYTHON
##
IF(MSVC11)
  VISIT_OPTION_DEFAULT(VISIT_PYTHON_DIR ${VISITHOME}/python/2.7.6)
ELSE()
  VISIT_OPTION_DEFAULT(VISIT_PYTHON_DIR ${VISITHOME}/python/2.7.5)
ENDIF()

##
## PYSIDE
##
VISIT_OPTION_DEFAULT(VISIT_PYSIDE_DIR ${VISITHOME}/PySide/1.2.2)

##############################################################
##
## Some misc libraries Database readers may depend upon
##
##############################################################

##
## JPEG
##
#IF(MSVC11)
#  IF (CMAKE_CL_64)
#     set(JPEG_LIBNAME_PREFIX_LIB false CACHE INTERNAL "Is 'lib' prefixed to jpeg libname")
#    VISIT_OPTION_DEFAULT(VISIT_JPEG_DIR     ${VISITHOME}/jpeg/9)
#  ENDIF()
#ELSE()
IF(MSVC10)
  set(JPEG_LIBNAME_PREFIX_LIB true CACHE INTERNAL "Is 'lib' prefixed to jpeg libname")
  VISIT_OPTION_DEFAULT(VISIT_JPEG_DIR     ${VISITHOME}/jpeg/8c)
      SET(JPEGVER "8c")
ENDIF()

##
## SZIP
##
IF(MSVC11)
  SET(SZIP_LIB_NAME "szip")
  VISIT_OPTION_DEFAULT(VISIT_SZIP_DIR     ${VISITHOME}/szip/2.1)
ELSE()
  SET(SZIP_LIB_NAME "szlibdll")
  VISIT_OPTION_DEFAULT(VISIT_SZIP_DIR     ${VISITHOME}/szip/2.0)
ENDIF()

##
## ZLIB
##
IF(MSVC11)
  SET(ZLIB_LIB_NAME "zlib")
  VISIT_OPTION_DEFAULT(VISIT_ZLIB_DIR     ${VISITHOME}/zlib/1.2.8)
ELSE()
  SET(ZLIB_LIB_NAME "zlib1")
  VISIT_OPTION_DEFAULT(VISIT_ZLIB_DIR     ${VISITHOME}/zlib/1.2.7)
ENDIF()

##############################################################
##
## Database reader plugin support libraries
##
## For libraries with LIBDEP settings, order matters.
## Libraries with LIBDEP settings that depend on other
## Library's LIBDEP settings must come after them.
##############################################################

##
## HDF4
##
IF(MSVC11)
  set(HDF4_LIBNAMES_AFFIX_DLL false CACHE INTERNAL "Is 'dll' appended to hdf4 libnames")
  VISIT_OPTION_DEFAULT(VISIT_HDF4_DIR     ${VISITHOME}/hdf4/4.2.10)
  VISIT_OPTION_DEFAULT(VISIT_HDF4_LIBDEP  
      VTK_LIBRARY_DIRS vtkjpeg-${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}
      SZIP_LIBRARY_DIR SZIP_LIB 
      ZLIB_LIBRARY_DIR ZLIB_LIBRARY
      TYPE STRING
  )
ELSE()
  set(HDF4_LIBNAMES_AFFIX_DLL true CACHE INTERNAL "Is 'dll' appended to hdf4 libnames")
  VISIT_OPTION_DEFAULT(VISIT_HDF4_DIR     ${VISITHOME}/hdf4/4.2.6)
  VISIT_OPTION_DEFAULT(VISIT_HDF4_LIBDEP  
    JPEG_LIBRARY_DIR JPEG_LIB 
    SZIP_LIBRARY_DIR SZIP_LIB 
    ZLIB_LIBRARY_DIR ZLIB_LIBRARY
    TYPE STRING
  )
ENDIF()

##
## HDF5
##
VISIT_OPTION_DEFAULT(VISIT_HDF5_DIR     ${VISITHOME}/hdf5/1.8.7)
VISIT_OPTION_DEFAULT(VISIT_HDF5_LIBDEP  
    SZIP_LIBRARY_DIR SZIP_LIB 
    ZLIB_LIBRARY_DIR ZLIB_LIBRARY
    TYPE STRING
)

##
## NETCDF
##
VISIT_OPTION_DEFAULT(VISIT_NETCDF_DIR   ${VISITHOME}/netcdf/4.1.1)
VISIT_OPTION_DEFAULT(VISIT_NETCDF_LIBDEP 
    HDF5_LIBRARY_DIR HDF5_LIB ${VISIT_HDF5_LIBDEP}
    TYPE STRING
)

##
## BOXLIB
##
IF(MSVC11)
  VISIT_OPTION_DEFAULT(VISIT_BOXLIB_DIR ${VISITHOME}/BoxLib/2011.04.28)
ELSE()
  VISIT_OPTION_DEFAULT(VISIT_BOXLIB_DIR ${VISITHOME}/BoxLib/0.1.8)
ENDIF()

##
## CCMIO
##
VISIT_OPTION_DEFAULT(VISIT_CCMIO_DIR    ${VISITHOME}/ccmio/2.6.1)

##
## CFITSIO
##
VISIT_OPTION_DEFAULT(VISIT_CFITSIO_DIR  ${VISITHOME}/cfitsio/3006)

##
## CGNS
##
IF(MSVC11)
  VISIT_OPTION_DEFAULT(VISIT_CGNS_DIR     ${VISITHOME}/cgns/3.2.1)
ELSE()
  VISIT_OPTION_DEFAULT(VISIT_CGNS_DIR     ${VISITHOME}/cgns/3.1.3)
ENDIF()
VISIT_OPTION_DEFAULT(VISIT_CGNS_LIBDEP 
    HDF5_LIBRARY_DIR HDF5_LIB ${VISIT_HDF5_LIBDEP}
    TYPE STRING
)

##
## GDAL
##
VISIT_OPTION_DEFAULT(VISIT_GDAL_DIR     ${VISITHOME}/gdal/1.10.0)

##
## H5PART
##
IF(MSVC11)
  VISIT_OPTION_DEFAULT(VISIT_H5PART_DIR   ${VISITHOME}/h5part/1.6.6)
ELSE()
  VISIT_OPTION_DEFAULT(VISIT_H5PART_DIR   ${VISITHOME}/h5part/1.6.5)
ENDIF()
VISIT_OPTION_DEFAULT(VISIT_H5PART_LIBDEP 
    HDF5_LIBRARY_DIR HDF5_LIB ${VISIT_HDF5_LIBDEP}
    TYPE STRING
)


##
## MILI
##

IF(CMAKE_CL_64 OR NOT MSVC11)
VISIT_OPTION_DEFAULT(VISIT_MILI_DIR     ${VISITHOME}/Mili/13.1.1-patch)
ENDIF()

##
## SILO
##
IF(MSVC11)
  VISIT_OPTION_DEFAULT(VISIT_SILO_DIR     ${VISITHOME}/silo/4.10.1)

  VISIT_OPTION_DEFAULT(VISIT_SILO_LIBDEP  
      HDF5_LIBRARY_DIR HDF5_LIB ${VISIT_HDF5_LIBDEP}
      TYPE STRING
  )
ENDIF()

##
## XDMF
##
if(VISIT_QT5)
  VISIT_OPTION_DEFAULT(VISIT_XDMF_DIR ${VISITHOME}/Xdmf-vtk-qt5/2.1.1)
else()
  VISIT_OPTION_DEFAULT(VISIT_XDMF_DIR ${VISITHOME}/Xdmf/2.1.1)
endif()
VISIT_OPTION_DEFAULT(VISIT_XDMF_LIBDEP 
    HDF5_LIBRARY_DIR HDF5_LIB ${VISIT_HDF5_LIBDEP}
    VTK_LIBRARY_DIRS vtklibxml2-${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}
    TYPE STRING
)
