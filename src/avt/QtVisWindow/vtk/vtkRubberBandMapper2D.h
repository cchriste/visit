/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkXPolyDataMapper2D.h,v $
  Language:  C++
  Date:      $Date: 2000/02/04 17:09:14 $
  Version:   $Revision: 1.11 $

Copyright (c) 1993-2000 Ken Martin, Will Schroeder, Bill Lorensen 
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

 * Neither name of Ken Martin, Will Schroeder, or Bill Lorensen nor the names
   of any contributors may be used to endorse or promote products derived
   from this software without specific prior written permission.

 * Modified source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ``AS IS''
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/
// .NAME vtkXPolyDataMapper2D - 2D PolyData support for windows
// .SECTION Description
// vtkXPolyDataMapper2D provides 2D PolyData annotation support for 
// vtk under windows.  Normally the user should use vtkPolyDataMapper2D 
// which in turn will use this class.

// .SECTION See Also
// vtkPolyDataMapper2D

#ifndef __vtkRubberBandMapper2D_h
#define __vtkRubberBandMapper2D_h
#include <qtviswindow_exports.h>

#include "vtkPolyDataMapper2D.h"

class QWidget;
struct vtkRubberBandMapper2DPrivate;

class QTVISWINDOW_API vtkRubberBandMapper2D : public vtkPolyDataMapper2D
{
public:
  vtkTypeMacro(vtkRubberBandMapper2D,vtkPolyDataMapper2D);
  static vtkRubberBandMapper2D *New();

  // Description:
  // Set the widget over which the drawing will happen.
  void SetWidget(QWidget *widget);

  // Description:
  // Actually draw the poly data.
  void RenderOverlay(vtkViewport* viewport, vtkActor2D* actor);

  // Description:
  // Release graphics resources.
  virtual void ReleaseGraphicsResources(vtkWindow *);

protected:
  vtkRubberBandMapper2D();
  ~vtkRubberBandMapper2D();

  void RenderOverlay_X11(vtkViewport* viewport, vtkActor2D* actor);
  void RenderOverlay_Qt(vtkViewport* viewport, vtkActor2D* actor);

  vtkRubberBandMapper2DPrivate *d;

private:
  vtkRubberBandMapper2D(const vtkRubberBandMapper2D&);
  void operator=(const vtkRubberBandMapper2D&);
  
};


#endif

