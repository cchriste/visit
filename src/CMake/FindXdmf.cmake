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
#   Kathleen Bonnell, Thu May 27 17:01:22 MST 2010
#   Windows builds can be handled the same as non-windows for XDMF.
#
#   Mark C. Miller, Fri Jul 30 22:06:06 PDT 2010
#   Removed logic setting XDMF_DIR (that is done to SET_UP_THIRD_PARTY)
#   as well as FIND_PACKAGE (also done by SET_UP_THIRD_PARTY).
#
#   Brad Whitlock, Fri Apr  6 11:00:10 PDT 2012
#   Also look for vtklibxml2 if we're building statically.
#
#   Cyrus Harrison, Fri Apr  6 11:00:10 PDT 2012
#   Static build: Only look for vtklibxml2 on OSX.
#
#   Cyrus Harrison, Tue Apr 10 13:07:08 PDT 2012
#   Revert to standard setup. Build_visit now handles vtk deps correctly
#   in the generated config-site.
#
#   Kathleen Biagas, Fri May 3 16:55:12 MST 2013
#   If our xdmf depends on vtlibxml2, ensure it exists.
#
#****************************************************************************/

# Use the XDMF_DIR hint from the config-site .cmake file 
#

IF(VISIT_XDMF_LIBDEP)
    LIST(FIND VISIT_XDMF_LIBDEP vtklibxml2-${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION} xdmf_needs_vtkxml)
    IF(${xdmf_needs_vtkxml} GREATER "-1" AND NOT TARGET vtklibxml2)
        MESSAGE(STATUS "Xdmf depends on vtklibxml2, but it doesn't exist")
        RETURN()
    ENDIF()
ENDIF()

INCLUDE(${VISIT_SOURCE_DIR}/CMake/SetUpThirdParty.cmake)

SET_UP_THIRD_PARTY(XDMF lib include Xdmf)

