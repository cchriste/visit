/*****************************************************************************
*
* Copyright (c) 2000 - 2006, The Regents of the University of California
* Produced at the Lawrence Livermore National Laboratory
* All rights reserved.
*
* This file is part of VisIt. For details, see http://www.llnl.gov/visit/. The
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
*    documentation and/or materials provided with the distribution.
*  - Neither the name of the UC/LLNL nor  the names of its contributors may be
*    used to  endorse or  promote products derived from  this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED.  IN  NO  EVENT  SHALL  THE  REGENTS  OF  THE  UNIVERSITY OF
* CALIFORNIA, THE U.S.  DEPARTMENT  OF  ENERGY OR CONTRIBUTORS BE  LIABLE  FOR
* ANY  DIRECT,  INDIRECT,  INCIDENTAL,  SPECIAL,  EXEMPLARY,  OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/

// ************************************************************************* //
//                             vtkVisItUtility.h                             //
// ************************************************************************* //

#ifndef VTK_VISIT_UTILITY_H
#define VTK_VISIT_UTILITY_H

#include <visit_vtk_light_exports.h>
#include <vtkType.h>

class vtkCell;
class vtkDataSet;
class vtkPoints;
class vtkRectilinearGrid;

// ****************************************************************************
//  Modifications:
//
//    Hank Childs, Sat Mar 18 14:16:09 PST 2006
//    Added function CellContainsPoint.
//
// ****************************************************************************

namespace vtkVisItUtility
{
    VISIT_VTK_LIGHT_API vtkPoints  *GetPoints(vtkDataSet *);
    VISIT_VTK_LIGHT_API void        GetLogicalIndices(vtkDataSet *, const bool, 
                                                const int, int [3], 
                                                const bool = false,
                                                const bool = true);
    VISIT_VTK_LIGHT_API int         CalculateRealID(const int, const bool, 
                                              vtkDataSet *ds);
    VISIT_VTK_LIGHT_API int         ComputeStructuredCoordinates(
                                              vtkRectilinearGrid *, 
                                              double x[3], int ijk[3]);
    VISIT_VTK_LIGHT_API int         FindCell(vtkDataSet *, double pt[3]);
    VISIT_VTK_LIGHT_API void        GetDimensions(vtkDataSet *, int[3]);
    VISIT_VTK_LIGHT_API int         NodeGhostIdFromNonGhost(vtkDataSet *ds,
                                        const int);
    VISIT_VTK_LIGHT_API int         ZoneGhostIdFromNonGhost(vtkDataSet *ds,
                                        const int);
    VISIT_VTK_LIGHT_API int         CalculateGhostIdFromNonGhost(
                                        vtkDataSet *ds,
                                        const int cellId, 
                                        const bool forCell);
    VISIT_VTK_LIGHT_API int         GetLocalElementForGlobal(
                                        vtkDataSet *ds,
                                        const int elementId, 
                                        const bool forCell);

    VISIT_VTK_LIGHT_API void       GetCellCenter(vtkCell* cell, double center[3]);
    VISIT_VTK_LIGHT_API bool       ContainsMixedGhostZoneTypes(vtkDataSet *);
    VISIT_VTK_LIGHT_API bool       CellContainsPoint(vtkCell *, const double *);
    VISIT_VTK_LIGHT_API void       WriteDataSet(vtkDataSet*, char *);
    VISIT_VTK_LIGHT_API vtkRectilinearGrid *Create1DRGrid(int nXCoords,
                                                          int type = VTK_FLOAT); 
    VISIT_VTK_LIGHT_API bool       PointsEqual(double p1[3], double p2[3]);
}

#endif


