/*****************************************************************************
*
* Copyright (c) 2000 - 2013, Lawrence Livermore National Security, LLC
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


//========================================================================
//
//  Program:   VisIt
//  Module:    $RCSfile: vtkSurfaceFilter.h,v $
//  Language:  C++
//  Date:      $Date: 2000/09/20 18:11:07 $
//  Version:   $Revision: 1.24 $
//
//
//
//=========================================================================
// .NAME vtkSurfaceFilter - 
// .SECTION Description
// vtkSurfaceFilter is a filter object that, 
//
// Output is DataSet, same concrete type as input except for inputs
// of rectilinear grids which are converted to unstructured grids for output.
//
// .SECTION See Also
//
// .SECTION Caveats
//

#ifndef __vtkSurfaceFilter_h
#define __vtkSurfaceFilter_h
#include <visit_vtk_exports.h>

#include <vtkUnstructuredGridAlgorithm.h>

class vtkDataAray;
class vtkPointSet;
class vtkRectilinearGrid;

//=======================================================================
// Modifications:
//   Kathleen Bonnell, Fri Feb  8 11:03:49 PST 2002
//   vtkScalars has been deprecated in VTK 4.0, use vtkDataArray instead.
//   Use ObjectMacro instead of plain Macro.
//=======================================================================

class VISIT_VTK_API  
vtkSurfaceFilter : public vtkUnstructuredGridAlgorithm
{
public:
  vtkTypeMacro(vtkSurfaceFilter,vtkUnstructuredGridAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  static vtkSurfaceFilter *New();

  // Description:
  // Set/Get the scalars to use for z-values in the output. 
  virtual void SetinScalars(vtkDataArray*); 
  vtkGetObjectMacro(inScalars, vtkDataArray); 
 
protected:
  vtkSurfaceFilter();
  ~vtkSurfaceFilter() ;

  virtual int RequestData(vtkInformation *,
                          vtkInformationVector **,
                          vtkInformationVector *);
  virtual int FillInputPortInformation(int port, vtkInformation *info);

  void ExecuteRectilinearGrid(vtkRectilinearGrid *, vtkUnstructuredGrid *);
  void ExecutePointSet(vtkPointSet *, vtkUnstructuredGrid *);

// Protected Data Members

  vtkDataArray *inScalars;

private:
  vtkSurfaceFilter(const vtkSurfaceFilter&);
  void operator=(const vtkSurfaceFilter&);
};

#endif
