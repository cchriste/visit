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


#ifndef __vtkCellIntersections_h
#define __vtkCellIntersections_h

#include <visit_vtk_light_exports.h>
#include <vtkObject.h>

class vtkCell;
class vtkVertex;
class vtkPolyVertex;
class vtkLine;
class vtkPolyLine;
class vtkTriangle;
class vtkTriangleStrip;
class vtkPolygon;
class vtkPixel;
class vtkQuad;
class vtkTetra;
class vtkVoxel;
class vtkHexahedron;
class vtkWedge;
class vtkPyramid;
class vtkQuadraticHexahedron;

class VISIT_VTK_LIGHT_API vtkCellIntersections : public vtkObject
{
public:
  vtkTypeMacro(vtkCellIntersections,vtkObject);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  static vtkCellIntersections *New();

  // Description:
  // Boolean controls whether to test for Co-Planar condition.
  vtkSetMacro(TestCoPlanar,int);
  vtkGetMacro(TestCoPlanar,int);
  vtkBooleanMacro(TestCoPlanar,int);

  int CellIntersectWithLine(vtkCell *, double [3], double [3], 
                                double&, double [3]);

  static int IntersectBox(const double[6], const double [3],
                          const double[3], double [3]);
  static int LineIntersectBox(const double[6], const double [3],
                          const double[3], double [3]);

protected:
  vtkCellIntersections();
  ~vtkCellIntersections();

private:
  vtkCellIntersections(const vtkCellIntersections&);  // Not implemented.
  void operator=(const vtkCellIntersections&);  // Not implemented.


  int VertexIntersectWithLine(vtkVertex *, double [3], double [3], 
                                double&, double [3]);

  int PolyVertexIntersectWithLine(vtkPolyVertex *, double [3], double [3], 
                                double&, double [3]);

  int LineIntersectWithLine(vtkLine *, double [3], double [3], 
                                double&, double [3]);

  int PolyLineIntersectWithLine(vtkPolyLine *, double [3], double [3], 
                                double&, double [3]);

  int TriangleIntersectWithLine(vtkTriangle *, double [3], double [3], 
                                double&, double [3]);

  int TriStripIntersectWithLine(vtkTriangleStrip *, double [3], double [3], 
                                double&, double [3]);

  int PolygonIntersectWithLine(vtkPolygon *, double [3], double [3], 
                                double&, double [3]);

  int PixelIntersectWithLine(vtkPixel *, double [3], double [3], 
                                double&, double [3]);

  int QuadIntersectWithLine(vtkQuad *, double [3], double [3], 
                                double&, double [3]);

  int TetraIntersectWithLine(vtkTetra *, double [3], double [3], 
                                double&, double [3]);

  int VoxelIntersectWithLine(vtkVoxel *, double [3], double [3], 
                                double&, double [3]);

  int HexIntersectWithLine(vtkHexahedron *, double [3], double [3], 
                                double&, double [3]);

  int WedgeIntersectWithLine(vtkWedge *, double [3], double [3], 
                                double&, double [3]);

  int PyramidIntersectWithLine(vtkPyramid *, double [3], double [3], 
                                double&, double [3]);

  int QuadraticHexahedronIntersectWithLine(vtkQuadraticHexahedron *, 
                                double [3], double [3], double&, double [3]);

  vtkTriangle *triangle;
  vtkQuad *quad;
  int TestCoPlanar;
};

#endif


