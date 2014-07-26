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
#   Kathleen Bonnell, Thu Dec  3 10:55:03 PST 2009
#   Wrap CMAKE_X_LIBS so that it won't parse on windows. Change ${MESA_FOUND}
#   to MESA_FOUND to remove cmake error.
#
#   Kathleen Bonnell, Wed Dec  9 15:13:27 MT 2009
#   Copy Mesa dlls to execution directory for OSMesa test on windows.
#
#   Kathleen Bonnell, Tue Jan  5 14:13:43 PST 2009
#   Use cmake 2.6.4 (rather than 2.8) compatible version of copying files.
#
#   Kathleen Bonnell, Tue Feb 16 14:00:02 MST 2010
#   Removed conditional check for OSMESA SIZE LIMIT, in case something wasn't
#   set up correctly during first configure pass (eg Mesa lib).
#
#   Kathleen Biagas, Tues Oct 1 09:33:47 MST 2013
#   Removed VISIT_MSVC_VERSION from windows handling.
#
#****************************************************************************/

# Use the VTK_DIR hint from the config-site .cmake file 
INCLUDE(${VISIT_SOURCE_DIR}/CMake/SetUpThirdParty.cmake)

IF (WIN32)
    SET_UP_THIRD_PARTY(MESA lib include MesaGL32 osmesa32)
ELSE (WIN32)
    SET_UP_THIRD_PARTY(MESA lib include OSMesa)

    # If we're on Apple, set up MesaGLU too. This is mostly to ensure that it gets installed.
    IF(APPLE)
        SET(MESAGLU_DIR ${MESA_DIR})
        SET_UP_THIRD_PARTY(MESAGLU lib include MesaGLU)
    ENDIF(APPLE)

    # Install Mesa headers
    IF(VISIT_MESA_SKIP_INSTALL)
        MESSAGE(STATUS "Skipping mesa installation")
    ELSE(VISIT_MESA_SKIP_INSTALL)
        IF(VISIT_HEADERS_SKIP_INSTALL)
            MESSAGE(STATUS "Skipping mesa headers installation")
        ELSE(VISIT_HEADERS_SKIP_INSTALL)
            INSTALL(DIRECTORY ${MESA_INCLUDE_DIR}
                DESTINATION ${VISIT_INSTALLED_VERSION_INCLUDE}/mesa
                FILE_PERMISSIONS OWNER_WRITE OWNER_READ GROUP_WRITE GROUP_READ WORLD_READ
                DIRECTORY_PERMISSIONS OWNER_WRITE OWNER_READ OWNER_EXECUTE GROUP_WRITE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
            )
        ENDIF(VISIT_HEADERS_SKIP_INSTALL)
    ENDIF(VISIT_MESA_SKIP_INSTALL)
ENDIF (WIN32)

IF(NOT WIN32)
  IF(NOT MESA_FOUND)
    MESSAGE(FATAL_ERROR "MESA is required to build VisIt")
  ENDIF(NOT MESA_FOUND)
# Need to have the mesa libs.
ELSE(NOT WIN32)
  IF(NOT MESA_FOUND)
    MESSAGE(WARNING "MESA not found.  Proceeding without.")
  ENDIF()
ENDIF(NOT WIN32)

SET(MY_LIBS ${MESA_LIB})

# Unix needs X_LIBS and THREAD_LIBS.
IF (NOT WIN32)
  IF (CMAKE_X_LIBS)
    SET(MY_LIBS ${MY_LIBS} ${CMAKE_X_LIBS})
  ENDIF (CMAKE_X_LIBS)
  MESSAGE(STATUS "Added unix libs.")
ENDIF (NOT WIN32)
IF (CMAKE_THREAD_LIBS)
    SET(MY_LIBS ${MY_LIBS} ${CMAKE_THREAD_LIBS})
ENDIF (CMAKE_THREAD_LIBS)


IF(MESA_FOUND)
  SET(MSG "Check for osmesa size limit")
  MESSAGE(STATUS ${MSG})
  SET(TRY_RUN_DIR ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}/CMakeTmp)

  IF (WIN32) 
    # Need these dlls to run the program
    IF(EXISTS ${MESA_LIBRARY_DIR}/MesaGL32.dll)
      EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy  ${MESA_LIBRARY_DIR}/MesaGL32.dll ${TRY_RUN_DIR}/CMakeFiles/CMakeTmp/debug/MesaGL32.dll)
      EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy  ${MESA_LIBRARY_DIR}/osmesa32.dll  ${TRY_RUN_DIR}/CMakeFiles/CMakeTmp/debug/osmesa32.dll)
    ENDIF(EXISTS ${MESA_LIBRARY_DIR}/MesaGL32.dll)
  ENDIF (WIN32) 
 
  TRY_RUN(TRY_RUN_RESULT HAVE_OSMESA_SIZE
    ${TRY_RUN_DIR}
    ${VISIT_SOURCE_DIR}/CMake/FindOSMesaSize.C
    CMAKE_FLAGS "-DINCLUDE_DIRECTORIES:STRING=${MESA_INCLUDE_DIR}"
                "-DLINK_DIRECTORIES:STRING=${MESA_LIBRARY_DIR}"
                "-DLINK_LIBRARIES:STRING=${MY_LIBS}"
    OUTPUT_VARIABLE OUTPUT
  )
  #MESSAGE(STATUS "${MSG} - OUTPUT_VARIABLE: ${OUTPUT}")

  IF (HAVE_OSMESA_SIZE)
    IF ("${TRY_RUN_RESULT}" MATCHES "FAILED_TO_RUN")
        MESSAGE(STATUS "${MSG} - failed to run, defaulting to 4096")
        SET(OSMESA_SIZE_LIMIT 4096)
    ELSE ("${TRY_RUN_RESULT}" MATCHES "FAILED_TO_RUN")
        IF (WIN32)
            SET(OSMESA_SIZE_LIMIT ${TRY_RUN_RESULT})
            MESSAGE(STATUS "${MSG} - found (${OSMESA_SIZE_LIMIT})")
            SET(HAVE_OSMESA_SIZE 1 CACHE INTERNAL "support for osmesa_size")
        ELSE (WIN32)
            IF (EXISTS ${CMAKE_BINARY_DIR}/junk.txt)
                FILE(STRINGS "${CMAKE_BINARY_DIR}/junk.txt" OSMESA_SIZE_LIMIT)
                FILE(REMOVE "${CMAKE_BINARY_DIR}/junk.txt")
                MESSAGE(STATUS "${MSG} - found (${OSMESA_SIZE_LIMIT})")
                SET(HAVE_OSMESA_SIZE 1 CACHE INTERNAL "support for osmesa_size")
            ELSE (EXISTS ${CMAKE_BINARY_DIR}/junk.txt)
                MESSAGE(STATUS "${MSG} - could not find junk.txt")
            ENDIF (EXISTS ${CMAKE_BINARY_DIR}/junk.txt)
        ENDIF (WIN32)
    ENDIF ("${TRY_RUN_RESULT}" MATCHES "FAILED_TO_RUN")
  ELSE(HAVE_OSMESA_SIZE)
    MESSAGE(STATUS "${MSG} - not found, defaulting to 4096")
    SET(HAVE_OSMESA_SIZE 0 CACHE INTERNAL "support for osmesa_size")
    SET(OSMESA_SIZE_LIMIT 4096)
  ENDIF()

    #
    # Create install symlinks so that we can use osmesa as libGL for an installed VisIt.
    #

    GET_FILENAME_COMPONENT(OSMESA_LIB_REAL ${MESA_LIBRARY_DIR}/${MESA_LIB} REALPATH)
    GET_FILENAME_COMPONENT(OSMESA_LIB_BASE ${OSMESA_LIB_REAL} NAME)

    INSTALL(CODE
            "EXECUTE_PROCESS(WORKING_DIRECTORY \$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}
            COMMAND ${CMAKE_COMMAND} -E make_directory ${VISIT_INSTALLED_VERSION_LIB}/osmesa/
            OUTPUT_VARIABLE OSMESA_DIR_OUT)
            MESSAGE(STATUS \"\${OSMESA_DIR_OUT}\")
            ")
    INSTALL(CODE
            "EXECUTE_PROCESS(WORKING_DIRECTORY \$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${VISIT_INSTALLED_VERSION_LIB}/osmesa/
            COMMAND ${CMAKE_COMMAND} -E remove -f libGL.so libGL.so.1
            OUTPUT_VARIABLE OSMESA_LINK_CLEAN_OUT)
            MESSAGE(STATUS \"\${OSMESA_LINK_CLEAN_OUT}\")
            ")
    INSTALL(CODE
            "EXECUTE_PROCESS(WORKING_DIRECTORY \$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${VISIT_INSTALLED_VERSION_LIB}/osmesa/
            COMMAND ${CMAKE_COMMAND} -E create_symlink ../${OSMESA_LIB_BASE} libGL.so
            OUTPUT_VARIABLE OSMESA_GL_SYMLINK)
            MESSAGE(STATUS \"\${OSMESA_GL_SYMLINK}\")
            ")

    INSTALL(CODE
            "EXECUTE_PROCESS(WORKING_DIRECTORY \$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${VISIT_INSTALLED_VERSION_LIB}/osmesa/
            COMMAND ${CMAKE_COMMAND} -E create_symlink ../${OSMESA_LIB_BASE} libGL.so.1
            OUTPUT_VARIABLE OSMESA_GL_SYMLINK)
            MESSAGE(STATUS \"\${OSMESA_GL_SYMLINK}\")
            ")
ELSE()
    MESSAGE(STATUS "Mesa not found, OSMESA_SIZE_LIMIT defaulting to 4096")
    SET(HAVE_OSMESA_SIZE 0 CACHE INTERNAL "support for osmesa_size")
    SET(OSMESA_SIZE_LIMIT 4096)
ENDIF ()

