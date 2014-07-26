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
#
#   Tom Fogal, Thu Mar 25 14:24:49 MDT 2010
#   Fix GDAL library naming convention on OS X.
#
#   Kathleen Bonnell, Tue Dec 21 14:45:51 MST 2010
#   Update gdal version to 1.7 on Windows.
#
#   Brad Whitlock, Fri Oct 14 10:56:28 PDT 2011
#   GDAL changed again on Mac.
#
#****************************************************************************/

# Use the GDAL_DIR hint from the config-site .cmake file 

INCLUDE(${VISIT_SOURCE_DIR}/CMake/SetUpThirdParty.cmake)

IF (WIN32)
    SET_UP_THIRD_PARTY(GDAL lib include gdal_i)
    # normally handled in InstallThirdParty.cmake, but gdal has a weird
    # naming convention on windows
    FOREACH(VER 17 19 110)
        IF(EXISTS ${GDAL_LIBRARY_DIR}/gdal${VER}.dll)
            EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy
                ${GDAL_LIBRARY_DIR}/gdal${VER}.dll
                ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/ThirdParty)
            INSTALL(FILES ${GDAL_LIBRARY_DIR}/gdal${VER}.dll
                DESTINATION ${VISIT_INSTALLED_VERSION_BIN}
                PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
                CONFIGURATIONS "" None Debug Release RelWithDebInfo MinSizeRel
                )
        ENDIF(EXISTS ${GDAL_LIBRARY_DIR}/gdal${VER}.dll)
    ENDFOREACH(VER)
ELSE (WIN32)
    SET_UP_THIRD_PARTY(GDAL lib include gdal)
ENDIF (WIN32)

