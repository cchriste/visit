/*****************************************************************************
*
* Copyright (c) 2000 - 2014, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-442911
* All rights reserved.
*
* This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
* full copyright notice is contained in the file COPYRIGHT located at the root
* of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
*
* Redistribution  and  use  in  source  and  binary  forms,  with  or  without
* modification, are permitted provided that the following conditions are met:
*
*  - Redistributions of  source code must  retain the above  copyright notice,
*    this list of conditions and the disclaimer below.
*  - Redistributions in binary form must reproduce the above copyright notice,
*    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
*    documentation and/or other materials provided with the distribution.
*  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
*    be used to endorse or promote products derived from this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
* LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
* DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/

#ifndef VISIT_CSG_MESH_H
#define VISIT_CSG_MESH_H

#ifdef __cplusplus
extern "C"
{
#endif

int VisIt_CSGMesh_alloc(visit_handle *obj);
int VisIt_CSGMesh_free(visit_handle obj);
int VisIt_CSGMesh_setRegions(visit_handle obj, visit_handle boolops, visit_handle leftids, visit_handle rightids);
int VisIt_CSGMesh_setZonelist(visit_handle obj, visit_handle zl);
int VisIt_CSGMesh_setBoundaryTypes(visit_handle obj, visit_handle csgtypes);
int VisIt_CSGMesh_setBoundaryCoeffs(visit_handle obj, visit_handle coeffs);
int VisIt_CSGMesh_setExtents(visit_handle obj, double min[3], double max[3]);

int VisIt_CSGMesh_getRegions(visit_handle h, visit_handle *typeflags,
                             visit_handle *left, visit_handle *right);
int VisIt_CSGMesh_getZonelist(visit_handle h, visit_handle *zl);
int VisIt_CSGMesh_getBoundaryTypes(visit_handle h, visit_handle *cshtypes);
int VisIt_CSGMesh_getBoundaryCoeffs(visit_handle h, visit_handle *coeffs);
int VisIt_CSGMesh_getExtents(visit_handle h, double min[3], double max[3]);

#ifdef __cplusplus
}
#endif

#endif
