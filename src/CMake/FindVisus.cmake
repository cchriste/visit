##**************************************************
## ViSUS Visualization Project                    **
## Copyright (c) 2010-2014 University of Utah     **
## Scientific Computing and Imaging Institute     **
## 72 S Central Campus Drive, Room 3750           **
## Salt Lake City, UT 84112                       **
##                                                **
## For information about this project see:        **
## http://www.pascucci.org/visus/                 **
##                                                **
##      or contact: pascucci@sci.utah.edu         **
##                                                **
##**************************************************/


# This module finds if ViSUS is installed, and sets the following variables
# indicating where it is:
#
# VISUS_FOUND               - system has ViSUS
# VISUS_INCLUDE_DIR         - path to Visus include directory
# VISUS_LIBRARIES           - all VISUS libraries
# VISUS_GLEW_INCLUDE_DIR    - path to glew included with visus
# VISUS_GLEW_LIBRARIES      - path to glew included with visus
# VISUS_TinyXML_INCLUDE_DIR - path to tinyxml included with visus
# VISUS_TinyXML_LIBRARIES   - path to tinyxml included with visus
#
# Execute cmake with "-DVISUS_DIR=/path/to/visus" to help find the library.
#


FIND_PATH(VISUS_INCLUDE_DIR  visuscpp    libs 
                                         ${VISUS_DIR}/include)

IF (VISUS_INCLUDE_DIR)

   ###########################################
   # Need to specify whether nvisusio was compiled against Qt or Juce.
   ###########################################

   SET(VISUS_GUI_LIBRARY "juce" CACHE STRING "one of juce, qt")
   SET_PROPERTY(CACHE VISUS_GUI_LIBRARY PROPERTY STRINGS juce qt)
   IF (${VISUS_GUI_LIBRARY} STREQUAL "qt")
     SET(VISUS_QT TRUE) 
   ELSEIF (${VISUS_GUI_LIBRARY} STREQUAL "juce") 
     SET(VISUS_JUCE TRUE)
   ENDIF()


   FIND_PATH(VISUS_GLEW_INCLUDE_DIR   glew.h 
                                ${VISUS_INCLUDE_DIR}/libs/glew/GL
   )

   FIND_PATH(VISUS_TinyXML_INCLUDE_DIR tinyxml.h 
                                ${VISUS_INCLUDE_DIR}/libs/tinyxml
   )

   FIND_LIBRARY(VISUS_KERNEL_LIB    visuskernel            ${VISUS_DIR}/lib)
   FIND_LIBRARY(VISUS_IDX_LIB       visusidx               ${VISUS_DIR}/lib)
   FIND_LIBRARY(VISUS_DB_LIB        visusdb                ${VISUS_DIR}/lib)
   FIND_LIBRARY(VISUS_DATAFLOW_LIB  visusdataflow          ${VISUS_DIR}/lib)
   #FIND_LIBRARY(VISUS_APPKIT_LIB    visusappkit            ${VISUS_DIR}/lib)
   FIND_LIBRARY(VISUS_SCENEGRAPH_LIB    visusscenegraph            ${VISUS_DIR}/lib)
   FIND_LIBRARY(VISUS_GUI_LIB       visusgui               ${VISUS_DIR}/lib)
   IF (VISUS_JUCE)
    FIND_LIBRARY(VISUS_GUI_IMPL_LIB Juce                   ${VISUS_DIR}/lib NO_DEFAULT_PATH)
   ENDIF()
   FIND_LIBRARY(VISUS_CURL_LIB      curl                   ${VISUS_DIR}/lib NO_DEFAULT_PATH)
   FIND_LIBRARY(VISUS_FREEIMAGE_LIB FreeImage              ${VISUS_DIR}/lib NO_DEFAULT_PATH)
   FIND_LIBRARY(VISUS_XML_LIB       tinyxml                ${VISUS_DIR}/lib NO_DEFAULT_PATH)
   FIND_LIBRARY(VISUS_LIBZ_LIB      libz                   ${VISUS_DIR}/lib NO_DEFAULT_PATH)
   FIND_LIBRARY(VISUS_GLEW_LIB      glew                   ${VISUS_DIR}/lib NO_DEFAULT_PATH)
   FIND_LIBRARY(VISUS_SSL_LIB       ssl                    ${VISUS_DIR}/lib NO_DEFAULT_PATH)
   FIND_LIBRARY(VISUS_CRYPTO_LIB    crypto                 ${VISUS_DIR}/lib NO_DEFAULT_PATH)

   SET(VISUS_CORE_LIBRARIES 
   #       ${VISUS_APPKIT_LIB}
       ${VISUS_SCENEGRAPH_LIB}
       ${VISUS_IDX_LIB}
       ${VISUS_GUI_LIB}
       ${VISUS_DATAFLOW_LIB}
       ${VISUS_DB_LIB}
       ${VISUS_KERNEL_LIB}
   )

   SET(VISUS_ADDL_LIBRARIES 
       ${VISUS_GUI_IMPL_LIB}
       ${VISUS_CURL_LIB}
       ${VISUS_FREEIMAGE_LIB}
       ${VISUS_XML_LIB}
       ${VISUS_LIBZ_LIB}
       ${VISUS_GLEW_LIB}
       ${VISUS_SSL_LIB}
       ${VISUS_CRYPTO_LIB}
   )

   SET(VISUS_GLEW_LIBRARIES     ${VISUS_GLEW_LIB})
   SET(VISUS_TinyXML_LIBRARIES  ${VISUS_XML_LIB})


   ###########################################
   # windows
   ###########################################

   IF (WIN32)
      SET(VISUS_WINDOWS 1)
      ADD_DEFINITIONS(-DVISUS_WINDOWS=1)

      #NOTE: These may cause problems in your build (e.g. for VisIt)
      #ADD_DEFINITIONS(-D_SCL_SECURE_NO_WARNINGS -D_CRT_SECURE_NO_WARNINGS -DWIN32_LEAN_AND_MEAN=1)
      #SET(CMAKE_CXX_FLAGS "/Oi ${CMAKE_CXX_FLAGS}") 

      SET (VISUS_LIBRARIES 
         ${VISUS_CORE_LIBRARIES}
         ${VISUS_ADDL_LIBRARIES}
         Vfw32.lib
         Version.lib
         Imm32.lib
         Winmm.lib;
         shlwapi.lib
         Wininet.lib
     )

   ###########################################
   # osx (desktop)
   ###########################################

   ELSEIF (APPLE)
     SET(VISUS_APPLE 1)
     SET(VISUS_OSX 1)
     ADD_DEFINITIONS(-DVISUS_APPLE=1 -DVISUS_OSX=1)  

     # On apple we need a bunch of other frameworks
     FIND_LIBRARY(COREMIDI_FRAMEWORK CoreMidi)
     FIND_LIBRARY(COREAUDIO_FRAMEWORK CoreAudio)
     FIND_LIBRARY(COREFOUNDATION_FRAMEWORK CoreFoundation)
     FIND_LIBRARY(QUARTZCORE_FRAMEWORK QuartzCore)
     FIND_LIBRARY(COCOA_FRAMEWORK Cocoa)
     FIND_LIBRARY(IOKIT_FRAMEWORK IoKit)
     IF (VISUS_QT)
        FIND_LIBRARY(QT_OPENGL_FRAMEWORK QtOpenGL)
        FIND_LIBRARY(QT_GUI_FRAMEWORK QtGui)
        FIND_LIBRARY(QT_CORE_FRAMEWORK QtCore)
     ENDIF()

     SET (VISUS_LIBRARIES 
         ${VISUS_CORE_LIBRARIES}
         ${VISUS_ADDL_LIBRARIES}
         ${COREMIDI_FRAMEWORK}
         ${COREAUDIO_FRAMEWORK}
         ${COREFOUNDATION_FRAMEWORK}
         ${QUARTZCORE_FRAMEWORK}
         ${COCOA_FRAMEWORK}
         ${IOKIT_FRAMEWORK}
         ${QT_OPENGL_FRAMEWORK}
         ${QT_GUI_FRAMEWORK}
         ${QT_CORE_FRAMEWORK}
     )

   ###########################################
   # linux
   ###########################################

   ELSEIF (UNIX)
     SET(VISUS_LINUX 1)
     ADD_DEFINITIONS(-DVISUS_LINUX=1)

     IF (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
       ADD_DEFINITIONS(-D_DEBUG=1)
     ENDIF()

     # On unix we need to ensure the libraries do not get pruned by the linker
     SET (VISUS_LIBRARIES 
       "-Wl,--whole-archive"
       ${VISUS_CORE_LIBRARIES}
       "-Wl,--no-whole-archive"
       ${VISUS_ADDL_LIBRARIES}
       freetype
       ssl
       rt
       pthread
       dl
     )
   ENDIF()


   # Set VISUS_FOUND and print summary messages
   SET(VISUS_FOUND "YES" CACHE BOOL "Visus library found" FORCE)
   IF (CMAKE_VERBOSE_MAKEFILE)
      MESSAGE(STATUS "Using VISUS_INCLUDE_DIR  = " ${VISUS_INCLUDE_DIR}) 
      MESSAGE(STATUS "Using VISUS_GLEW_INCLUDE_DIR   = " ${VISUS_GLEW_INCLUDE_DIR}) 
      MESSAGE(STATUS "Using VISUS_TinyXML_INCLUDE_DIR= " ${VISUS_TinyXML_INCLUDE_DIR}) 
      MESSAGE(STATUS "Found VISUS_LIBRARIES    = " ${VISUS_LIBRARIES}) 
      FOREACH(lib ${VISUS_LIBRARIES})
         MESSAGE(STATUS "is: " ${lib})
      ENDFOREACH()
      MESSAGE(STATUS "Found VISUS_GLEW_LIBRARIES     = " ${VISUS_GLEW_LIBRARIES}) 
      MESSAGE(STATUS "Found VISUS_TinyXML_LIBRARIES  = " ${VISUS_TinyXML_LIBRARIES}) 
   ENDIF (CMAKE_VERBOSE_MAKEFILE)

ELSE ()
   # change this to an error when logic is added such that we know if Visus was requested.
   MESSAGE(STATUS "Visus library not found on the system")
ENDIF ()
                         
