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
#   Kathleen Bonnell, Thu Dec  3 10:30:15 PST 2009
#   Use only 1 if-def block, fix libarary names for windows.
#
#   Eric Brugger, Fri Jan  7 13:50:15 PST 2011
#   I replaced the BOXLIB2D and BOXLIB3D variables with just BOXLIB.
#
#   Kathleen Bonnell, Thu Jan 13 15:21:47 MST 2011
#   Restore separate vars for libraries (to handle different names on
#   different platforms).
#
#   Kathleen Bonnell, Mon Jan 17 17:24:44 MST 2011
#   Don't set BOXLIB_2D/3D_LIB unless BOXLIB_FOUND.
#
#   Kathleen Biagas, Tues Oct 1 09:33:47 MST 2013
#   Removed VISIT_MSVC_VER from Windows handling.
#
#****************************************************************************/

# Use the BOXLIB_DIR hint from the config-site .cmake file 

INCLUDE(${VISIT_SOURCE_DIR}/CMake/SetUpThirdParty.cmake)

IF (WIN32)
  SET_UP_THIRD_PARTY(BOXLIB lib include BoxLib2D BoxLib3D)
ELSE (WIN32)
  SET_UP_THIRD_PARTY(BOXLIB lib include box2D box3D)
ENDIF (WIN32)

IF(BOXLIB_FOUND)
  # place the 2D and 3D libraries into separate vars for plugin use.
  LIST(GET BOXLIB_LIB 0 tmp)
  SET(BOXLIB_2D_LIB ${tmp} CACHE STRING "2D boxlib" FORCE)

  LIST(GET BOXLIB_LIB 1 tmp)
  SET(BOXLIB_3D_LIB ${tmp} CACHE STRING "3D boxlib" FORCE)

  # unset unneeded vars.
  UNSET(tmp)
  UNSET(BOXLIB_LIB CACHE)
ENDIF(BOXLIB_FOUND)


