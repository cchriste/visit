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

#ifndef VTK_VISIT_INTERPOLATED_VELOCITY_FIELD
#define VTK_VISIT_INTERPOLATED_VELOCITY_FIELD

#include <vtkObject.h>
#include <visit_vtk_exports.h>

class vtkDataSet;
class vtkVisItCellLocator;

class VISIT_VTK_API vtkVisItInterpolatedVelocityField  : public vtkObject
{
  public:
  vtkTypeMacro(vtkVisItInterpolatedVelocityField,vtkObject);
                  vtkVisItInterpolatedVelocityField();
    virtual      ~vtkVisItInterpolatedVelocityField();

    static vtkVisItInterpolatedVelocityField *New();
    void          SetDataSet(vtkDataSet *ds);
    void          SetLocator(vtkVisItCellLocator *loc);
    bool          Evaluate(double *pt, double *vel, double t = 0.);
 
    void          SetDoPathlines(bool b) { doPathlines = b; };
    void          SetCurrentTime(double t) { curTime = t; };
    void          SetNextTime(double t) { nextTime = t; };
    void          SetNextTimeName(const std::string &s) { nextTimeName = s; };

    vtkDataSet   *GetDataSet(void) { return ds; };
    double       *GetLastWeights(void) { return weights; };
    int           GetLastCell(void) { return lastCell; };
    double       *GetLastPCoords(void) { return pcoords; };

  private:
    vtkDataSet  *ds;
    vtkVisItCellLocator *locator;
    double       weights[1024];
    int          lastCell;
    double       pcoords[3];
    bool         doPathlines;
    std::string  nextTimeName;
    double       curTime;
    double       nextTime;
};

#endif

