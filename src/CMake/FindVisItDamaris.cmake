#*****************************************************************************
#
# Copyright (c) 2000 - 2015, Lawrence Livermore National Security, LLC
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
#****************************************************************************/

# Use the DAMARIS_DIR hint from the config-site .cmake file 

INCLUDE(${VISIT_SOURCE_DIR}/CMake/SetUpThirdParty.cmake)

SET_UP_THIRD_PARTY(DAMARIS lib include damaris)
SET_UP_THIRD_PARTY(XERCESC lib include xerces-c)
SET_UP_THIRD_PARTY(XSD lib include NO_LIBS)

IF(DAMARIS_FOUND AND XERCESC_FOUND AND XSD_FOUND)

SET(BOOST_LIBS
        boost_date_time
        boost_filesystem
        boost_system)
SET_UP_THIRD_PARTY(BOOST lib include ${BOOST_LIBS})

INSTALL(DIRECTORY ${DAMARIS_INCLUDE_DIR}
                DESTINATION ${VISIT_INSTALLED_VERSION_INCLUDE}/../damaris
                FILE_PERMISSIONS OWNER_WRITE OWNER_READ GROUP_WRITE GROUP_READ WORLD_READ
                DIRECTORY_PERMISSIONS OWNER_WRITE OWNER_READ OWNER_EXECUTE GROUP_WRITE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
            )
INSTALL(FILES ${DAMARIS_LIBRARY_DIR}/libdamaris.a
                DESTINATION ${VISIT_INSTALLED_VERSION_LIB}/../damaris/lib
                PERMISSIONS OWNER_WRITE OWNER_READ GROUP_WRITE GROUP_READ WORLD_READ
            )
INSTALL(DIRECTORY ${XERCESC_INCLUDE_DIR}
                DESTINATION ${VISIT_INSTALLED_VERSION_INCLUDE}/../damaris
                FILE_PERMISSIONS OWNER_WRITE OWNER_READ GROUP_WRITE GROUP_READ WORLD_READ
                DIRECTORY_PERMISSIONS OWNER_WRITE OWNER_READ OWNER_EXECUTE GROUP_WRITE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
            )
INSTALL(FILES ${XERCESC_LIBRARY_DIR}/libxerces-c.a
                DESTINATION ${VISIT_INSTALLED_VERSION_LIB}/../damaris/lib
                PERMISSIONS OWNER_WRITE OWNER_READ GROUP_WRITE GROUP_READ WORLD_READ
            )
INSTALL(DIRECTORY ${XSD_INCLUDE_DIR}
                DESTINATION ${VISIT_INSTALLED_VERSION_INCLUDE}/../damaris
                FILE_PERMISSIONS OWNER_WRITE OWNER_READ GROUP_WRITE GROUP_READ WORLD_READ
                DIRECTORY_PERMISSIONS OWNER_WRITE OWNER_READ OWNER_EXECUTE GROUP_WRITE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
            )
INSTALL(DIRECTORY ${BOOST_INCLUDE_DIR}
                DESTINATION ${VISIT_INSTALLED_VERSION_INCLUDE}/../damaris
                FILE_PERMISSIONS OWNER_WRITE OWNER_READ GROUP_WRITE GROUP_READ WORLD_READ
                DIRECTORY_PERMISSIONS OWNER_WRITE OWNER_READ OWNER_EXECUTE GROUP_WRITE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
            )
INSTALL(FILES ${BOOST_LIBRARY_DIR}/libboost_date_time.a
	      ${BOOST_LIBRARY_DIR}/libboost_filesystem.a
	      ${BOOST_LIBRARY_DIR}/libboost_system.a
                DESTINATION ${VISIT_INSTALLED_VERSION_LIB}/../damaris/lib
                PERMISSIONS OWNER_WRITE OWNER_READ GROUP_WRITE GROUP_READ WORLD_READ
            )
ENDIF(DAMARIS_FOUND AND XERCESC_FOUND AND XSD_FOUND)
