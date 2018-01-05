#*****************************************************************************
#
# Copyright (c) 2000 - 2017, Lawrence Livermore National Security, LLC
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
#   Kathleen Bonnell, Mon Dec 27, 17:52:39 MST 2010
#   Added high-level hdf5 lib to search on Windows. (hdf5_hldll).
#
#   Kathleen Biagas, Wed Oct 19 09:58:16 MST 2011
#   Remove ${VISIT_MSVC_VERSION} from lib location.
#
#   Kathleen Biagas, Thu Jan 9 18:47:21 PDT 2014
#   Add patch from John Cary for hdf5 without 'dll' suffix on name.
#
#****************************************************************************/

# Use the HDF5_DIR hint from the config-site .cmake file 

INCLUDE(${VISIT_SOURCE_DIR}/CMake/SetUpThirdParty.cmake)

OPTION(HDF5_LIBNAMES_AFFIX_DLL "Whether HDF5 library base names end with dll" ON)
IF(WIN32)
  if(HDF5_LIB_NAME)
    SET_UP_THIRD_PARTY(HDF5 lib include ${HDF5_LIB_NAME})
    IF(VISIT_PARALLEL)
        SET_UP_THIRD_PARTY(HDF5_MPI lib include ${HDF5_LIB_NAME})
    ENDIF(VISIT_PARALLEL)
  else()
    if(HDF5_LIBNAMES_AFFIX_DLL)
      SET_UP_THIRD_PARTY(HDF5 lib include hdf5dll hdf5_hldll)
      IF(VISIT_PARALLEL)
          SET_UP_THIRD_PARTY(HDF5_MPI lib include hdf5_mpidll hdf5_mpi_hldll)
      ENDIF(VISIT_PARALLEL)
    else()
      SET_UP_THIRD_PARTY(HDF5 lib include hdf5 hdf5_hl)
      IF(VISIT_PARALLEL)
          SET_UP_THIRD_PARTY(HDF5_MPI lib include hdf5_mpi hdf5_mpi_hl)
      ENDIF(VISIT_PARALLEL)
    endif()
  endif()
  if(MSVC_VERSION GREATER_EQUAL "1910")
  # using a version of hdf5 that requires this define:
      add_definitions(-DH5_BUILT_AS_DYNAMIC_LIB)
  endif()
ELSE()
  SET_UP_THIRD_PARTY(HDF5 lib include hdf5)
  IF(VISIT_PARALLEL)
      SET_UP_THIRD_PARTY(HDF5_MPI lib include hdf5_mpi)
  ENDIF(VISIT_PARALLEL)
ENDIF()
