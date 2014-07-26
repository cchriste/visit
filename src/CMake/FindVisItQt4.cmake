#*****************************************************************************
#
# Copyright (c) 2000 - 2014, Lawrence Livermore National Security, LLC
# Produced at the Lawrence Livermore National Laboratory
# LLNL-CODE-442911
# All rights reserved.
#
# This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
# full copyright notice is contained in the file COPYRIGHT located at the root
# of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
#
# Redistribution  and  use  in  source  and  binary  forms,  with  or  without
# modification, are permitted provided that the following conditions are met:
#
#  - Redistributions of  source code must  retain the above  copyright notice,
#    this list of conditions and the disclaimer below.
#  - Redistributions in binary form must reproduce the above copyright notice,
#    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
#    documentation and/or other materials provided with the distribution.
#  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
#    be used to endorse or promote products derived from this software without
#    specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
# ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
# LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
# DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
# SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
# CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
# LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
# OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.
#
# Modifications:
#   Kathleen Bonnell, Thu Dec 10 17:51:27 MT 2009
#   Use QT_X_LIBRARY_RELEASE instead of simple QT_X_LIBRARY, which may 
#   list both optimized and debug versions of the library if present.
#
#   Tom Fogal, Mon May  3 10:22:15 MDT 2010
#   Remove QtScript requirement/installation.  We don't use it.
#
#   Kathleen Bonnell, Thu Dec 2 15:12:44 MST 2010
#   Install moc on Windows.  Add all include dirs on Windows.
# 
#   Kathleen Bonnell, Thu Feb 3 08:21:18 PST 2010
#   Allow for using installed QT on Windows (follow same path as on *nix).
#   Simplified code when using windowsbuild version of QT.
#
#   Cyrus Harrisond, Tue Oct  4 16:18:24 PDT 2011
#   Add back QtScript, we may want to use it in conjunction w/ 
#   PySide & QtDesigner
#
#   Brad Whitlock, Wed Apr 10 18:04:33 PDT 2013
#   Fix Qt dependent libraries to include some extra frameworks. These were
#   needed on a 10.8 machine.
#
#   Kathleen Biagas, Tues Oct 1 09:33:47 MST 2013
#   Removed VISIT_MSVC_VERSION from windows handling.
#
#****************************************************************************/

#
# Use the QT_BIN hint from the config-site cmake file along
# with the standard FindQt4 cmake module to set up Qt4.
#

IF(NOT "${QT_BIN}" MATCHES "OFF")
  IF(WIN32)
    IF(VISIT_MSVC_VERSION AND EXISTS ${QT_DIR}/lib)
      # using VisIt's windowsbuild Qt
      SET(USE_CMAKE_FIND OFF)
    ELSE()
      # using other Qt
      SET(USE_CMAKE_FIND ON)
    ENDIF()
  ELSE()
    SET(USE_CMAKE_FIND ON)
  ENDIF(WIN32)

  IF (USE_CMAKE_FIND)
    # Make sure the VISIT_QT_BIN path is valid & qmake exists.
    FIND_PROGRAM(VISIT_LOC_QMAKE_EXE NAMES qmake qmake4 qmake-qt4
                 PATHS
                 ${QT_BIN}
                 NO_DEFAULT_PATH
                 NO_CMAKE_ENVIRONMENT_PATH
                 NO_CMAKE_PATH
                 NO_SYSTEM_ENVIRONMENT_PATH)

    IF ( NOT VISIT_LOC_QMAKE_EXE)
      MESSAGE(FATAL_ERROR "Invalid Qt4 Binary path: ${QT_BIN}")
    ENDIF ( NOT VISIT_LOC_QMAKE_EXE)

    # add VISIT_QT_BIN to the env, so standard FindQt4 module can locate qmake.
    SET(ENV{PATH} "${QT_BIN}:$ENV{PATH}")

    # Invoke cmake's built in module for locating & setting up Qt4
    INCLUDE(${CMAKE_ROOT}/Modules/FindQt4.cmake)

    IF(NOT QT_FOUND)
      MESSAGE(FATAL_ERROR "Qt4 is required to build VisIt.")
    ENDIF(NOT QT_FOUND)
  ELSE (USE_CMAKE_FIND)
    # MESSAGE("QT_DIR = ${QT_DIR}")
    SET(QT_INCLUDE_DIR ${QT_DIR}/include)
    SET(QT_LIBRARY_DIR ${QT_DIR}/lib
        CACHE PATH "Qt library dir" FORCE )
    SET(QT_BINARY_DIR  ${QT_DIR}/lib
        CACHE INTERNAL "" FORCE )

    SET(QT_MOC_EXECUTABLE  ${QT_BINARY_DIR}/moc.exe)
    SET(QT_INCLUDES ${QT_INCLUDE_DIR})
    
    SET(QT_WIN_LIBS QtDesigner QtDesignerComponents QtSql QtSvg Qt QtTest 
                    QtMain QtAssistantClient QtHelp QtXMLPatterns QtUiTools
                    QtCore QtGui QtOpenGL QtNetwork QtXml)
    FOREACH(QTWINLIB ${QT_WIN_LIBS})
      STRING(TOUPPER ${QTWINLIB} upper_qtwinlib)
      SET(QT_${upper_qtwinlib}_FOUND 1)
      SET(QT_${upper_qtwinlib}_INCLUDE_DIR ${QT_INCLUDE_DIR}/${QTWINLIB} 
          CACHE PATH "The Qt ${QTWINLIB} include dir" FORCE)
      IF(EXISTS ${QT_${upper_qtwinlib}_INCLUDE_DIR})
        SET(QT_INCLUDES ${QT_INCLUDES} ${QT_${upper_qtwinlib}_INCLUDE_DIR})
      ENDIF()
      IF (EXISTS ${QT_LIBRARY_DIR}/${QTWINLIB}4.lib)
        SET(QT_${upper_qtwinlib}_LIBRARY
            ${QT_LIBRARY_DIR}/${QTWINLIB}4.lib CACHE STRING
            "The Qt ${QTWINLIB} library" FORCE)
        SET(QT_${upper_qtwinlib}_LIBRARY_RELEASE
            ${QT_LIBRARY_DIR}/${QTWINLIB}4.lib)
      ELSE ()
        SET(QT_${upper_qtwinlib}_LIBRARY
            ${QT_LIBRARY_DIR}/${QTWINLIB}.lib CACHE STRING
            "The Qt ${QTWINLIB} library" FORCE)
        SET(QT_${upper_qtwinlib}_LIBRARY_RELEASE
            ${QT_LIBRARY_DIR}/${QTWINLIB}.lib)
      ENDIF ()
    ENDFOREACH(QTWINLIB)
  ENDIF (USE_CMAKE_FIND)

  #
  # If we are using cocoa we need to define VISIT_MAC_NO_CARBON
  #
  IF(APPLE)
    IF(QT_MAC_USE_COCOA)
      ADD_DEFINITIONS(-DVISIT_MAC_NO_CARBON)

      IF(VISIT_STATIC)
         SET(QT_QTCORE_LIB_DEPENDENCIES ${QT_QTCORE_LIB_DEPENDENCIES} "-framework Security" "-framework AppKit" -lobjc)
         SET(QT_QTNETWORK_LIB_DEPENDENCIES ${QT_QTNETWORK_LIB_DEPENDENCIES} "-framework SystemConfiguration")
      ENDIF(VISIT_STATIC)
    ENDIF(QT_MAC_USE_COCOA)
  ENDIF(APPLE)

  IF(VISIT_QT_SKIP_INSTALL)
    MESSAGE(STATUS "Skipping installation of Qt headers and libraries..")
  ELSE (VISIT_QT_SKIP_INSTALL)
    # Since Qt was found, add install targets for its libraries.
    FOREACH(QTLIB
          QT_QT3SUPPORT
          QT_QTASSISTANT
          QT_QAXCONTAINER
          QT_QAXSERVER
          QT_QTCORE
          QT_QTDBUS
          QT_QTDESIGNER
          QT_QTDESIGNERCOMPONENTS
          QT_QTGUI
          QT_QTMOTIF
          QT_QTNETWORK
          QT_QTNSPLUGIN
          QT_QTOPENGL
          QT_QTSQL
          QT_QTXML
          QT_QTSVG
          QT_QTTEST
          QT_QTMAIN
          QT_QTUITOOLS
          QT_QTASSISTANTCLIENT
          QT_QTHELP
          QT_QTWEBKIT
          QT_QTXMLPATTERNS
          QT_PHONON
          QT_QTSCRIPT
    )
        IF(${${QTLIB}_FOUND})
            IF(EXISTS ${${QTLIB}_LIBRARY_RELEASE})
                THIRD_PARTY_INSTALL_LIBRARY(${${QTLIB}_LIBRARY_RELEASE})
            ENDIF(EXISTS ${${QTLIB}_LIBRARY_RELEASE})
        ENDIF(${${QTLIB}_FOUND})
    ENDFOREACH(QTLIB)

    # Add install targets for Qt headers too
    FOREACH(H ${QT_INCLUDES})
        IF(${H} MATCHES "/include/Qt")
        INSTALL(DIRECTORY ${H}
                DESTINATION ${VISIT_INSTALLED_VERSION_INCLUDE}/qt/include
                FILE_PERMISSIONS OWNER_WRITE OWNER_READ
                                   GROUP_WRITE GROUP_READ
                                   WORLD_READ
                DIRECTORY_PERMISSIONS OWNER_WRITE OWNER_READ OWNER_EXECUTE
                                        GROUP_WRITE GROUP_READ GROUP_EXECUTE
                                        WORLD_READ WORLD_EXECUTE
                PATTERN ".svn" EXCLUDE
        )
        ENDIF(${H} MATCHES "/include/Qt")
    ENDFOREACH(H)

    # Install moc, too
    IF(NOT WIN32)
        INSTALL(PROGRAMS ${QT_MOC_EXECUTABLE} #${QT_BIN}/moc
                DESTINATION ${VISIT_INSTALLED_VERSION_BIN}
                PERMISSIONS OWNER_WRITE OWNER_READ OWNER_EXECUTE
                            GROUP_WRITE GROUP_READ GROUP_EXECUTE
                            WORLD_READ WORLD_EXECUTE
        )
    ELSE(NOT WIN32)
        INSTALL(PROGRAMS ${QT_MOC_EXECUTABLE}
                DESTINATION ${VISIT_INSTALLED_VERSION_BIN}
                PERMISSIONS OWNER_WRITE OWNER_READ OWNER_EXECUTE
                            GROUP_WRITE GROUP_READ GROUP_EXECUTE
                            WORLD_READ WORLD_EXECUTE
        )
    ENDIF(NOT WIN32)
  ENDIF(VISIT_QT_SKIP_INSTALL)
ENDIF(NOT "${QT_BIN}" MATCHES "OFF")
