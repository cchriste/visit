#*****************************************************************************
#
# Copyright (c) 2000 - 2009, Lawrence Livermore National Security, LLC
# Produced at the Lawrence Livermore National Laboratory
# LLNL-CODE-400142
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
#****************************************************************************/

#
# This function installs a library and any of its needed symlink variants.
#

FUNCTION(THIRD_PARTY_INSTALL_LIBRARY LIBFILE)
    IF(WIN32)
        IF(NOT EXISTS ${EXECUTABLE_OUTPUT_PATH}/ThirdParty)
            FILE(MAKE_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}/ThirdParty)
        ENDIF(NOT EXISTS ${EXECUTABLE_OUTPUT_PATH}/ThirdParty)
    ENDIF(WIN32)
    SET(tmpLIBFILE ${LIBFILE})
    GET_FILENAME_COMPONENT(LIBEXT ${tmpLIBFILE} EXT)
    IF(NOT ${LIBEXT} STREQUAL ".a")
        GET_FILENAME_COMPONENT(LIBREALPATH ${tmpLIBFILE} REALPATH)
#        MESSAGE("***tmpLIBFILE=${tmpLIBFILE}, LIBPATH=${LIBPATH}, LIBREALPATH=${LIBREALPATH}")
        IF(NOT ${tmpLIBFILE} STREQUAL ${LIBREALPATH})
            # We need to install a library and its symlinks
            GET_FILENAME_COMPONENT(curPATH ${LIBREALPATH} PATH)
            IF((NOT ${curPATH} STREQUAL "/usr/lib") AND (NOT ${curPATH} MATCHES "^\\/System\\/Library\\/Frameworks\\/.*"))
                GET_FILENAME_COMPONENT(curNAMEWE ${LIBREALPATH} NAME_WE)
                GET_FILENAME_COMPONENT(curEXT ${LIBREALPATH} EXT)
                STRING(REPLACE "." ";" extList ${curEXT})
                SET(curNAME "${curPATH}/${curNAMEWE}")
                FOREACH(X ${extList})
                    SET(curNAME "${curNAME}.${X}")
                    IF(EXISTS ${curNAME})
#                        MESSAGE("** Need to install ${curNAME}")
                        IF(IS_DIRECTORY ${curNAME})
                            INSTALL(DIRECTORY ${curNAME}
                                DESTINATION ${VISIT_INSTALLED_VERSION_LIB}
                                DIRECTORY_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
                                FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
                            )
                        ELSE(IS_DIRECTORY ${curNAME})
                            INSTALL(FILES ${curNAME}
                                DESTINATION ${VISIT_INSTALLED_VERSION_LIB}
                                PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
                            )
                            # On Windows, we also need to copy the file to the 
                            # binary dir so our out of source builds can run.
                            IF(WIN32)
                                EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy
                                        ${curPATH}/${curNAMEWE}.dll
                                        ${EXECUTABLE_OUTPUT_PATH}/ThirdParty)
                            ENDIF(WIN32)
                        ENDIF(IS_DIRECTORY ${curNAME})
                    ENDIF(EXISTS ${curNAME})
                ENDFOREACH(X)
            ENDIF((NOT ${curPATH} STREQUAL "/usr/lib") AND (NOT ${curPATH} MATCHES "^\\/System\\/Library\\/Frameworks\\/.*"))
        ELSE(NOT ${tmpLIBFILE} STREQUAL ${LIBREALPATH})
            # We need to install just the library
            IF(IS_DIRECTORY ${tmpLIBFILE})
                # It is a framework, install as a directory.
                INSTALL(DIRECTORY ${tmpLIBFILE}
                    DESTINATION ${VISIT_INSTALLED_VERSION_LIB}
                    DIRECTORY_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
                    FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
                )
            ELSE(IS_DIRECTORY ${tmpLIBFILE})
                # Create an install target for just the library file
                INSTALL(FILES ${tmpLIBFILE}
                    DESTINATION ${VISIT_INSTALLED_VERSION_LIB}
                    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
                )

                # On Windows, we also need to copy the file to the binary dir so
                # our out of source builds can run.
                IF(WIN32)
                    # determine the corresponding dll to the lib
                    GET_FILENAME_COMPONENT(curNAMEWE ${tmpLIBFILE} NAME_WE)
                    GET_FILENAME_COMPONENT(curPATH ${tmpLIBFILE} PATH)
                    EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy
                            ${curPATH}/${curNAMEWE}.dll
                            ${EXECUTABLE_OUTPUT_PATH}/ThirdParty)
                ENDIF(WIN32)
            ENDIF(IS_DIRECTORY ${tmpLIBFILE})
#            MESSAGE("**We need to install lib ${tmpLIBFILE}")
        ENDIF(NOT ${tmpLIBFILE} STREQUAL ${LIBREALPATH})
    ELSE(NOT ${LIBEXT} STREQUAL ".a")
        # We have a .a that we need to install to archives.
        IF(VISIT_INSTALL_THIRD_PARTY)
#            MESSAGE("***INSTALL ${LIBFILE} to ${VISIT_INSTALLED_VERSION_ARCHIVES}")
            INSTALL(FILES ${tmpLIBFILE}
                DESTINATION ${VISIT_INSTALLED_VERSION_ARCHIVES}
                PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ GROUP_WRITE WORLD_READ
            )

            # TODO: We could install windows import libraries here...

        ENDIF(VISIT_INSTALL_THIRD_PARTY)
    ENDIF(NOT ${LIBEXT} STREQUAL ".a")
ENDFUNCTION(THIRD_PARTY_INSTALL_LIBRARY)

#
# This function installs a library's includes.
#

FUNCTION(THIRD_PARTY_INSTALL_INCLUDE pkg incdir)
    IF(NOT WIN32)
        IF(VISIT_INSTALL_THIRD_PARTY)
            STRING(TOLOWER ${pkg} lcpkg)
#            MESSAGE("***INSTALL ${incdir} -> ${VISIT_INSTALLED_VERSION_INCLUDE}/${lcpkg}")
            INSTALL(DIRECTORY ${incdir}
                DESTINATION ${VISIT_INSTALLED_VERSION_INCLUDE}/${lcpkg}
                DIRECTORY_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
                FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
                FILES_MATCHING 
                PATTERN "*.h"
                PATTERN "*.H"
                PATTERN "*.hpp"
                PATTERN "*.HPP"
                PATTERN "*.inc"
                PATTERN "libccmio" EXCLUDE
            )
        ENDIF(VISIT_INSTALL_THIRD_PARTY)
    ENDIF(NOT WIN32)
ENDFUNCTION(THIRD_PARTY_INSTALL_INCLUDE)
