/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkStructuredGridMapper.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkStructuredGridMapper - map vtkStructuredGrid to graphics primitives
// .SECTION Description
// vtkStructuredGridMapper is a class that maps curvilinear data (i.e., vtkStructuredGrid)
// to graphics primitives. vtkStructuredGridMapper serves as a superclass for
// device-specific curvilinear data mappers, that actually do the mapping to the
// rendering/graphics hardware/software.

#ifndef __vtkStructuredGridMapper_h
#define __vtkStructuredGridMapper_h

#include "vtkMapper.h"

class vtkStructuredGrid;
class vtkRenderer;

class VTK_RENDERING_EXPORT vtkStructuredGridMapper : public vtkMapper 
{
public:
  static vtkStructuredGridMapper *New();
  vtkTypeRevisionMacro(vtkStructuredGridMapper,vtkMapper);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Specify the input data to map.
  void SetInput(vtkStructuredGrid *in);
  vtkStructuredGrid *GetInput();
  
  // Description:
  // Make a shallow copy of this mapper.
  void ShallowCopy(vtkAbstractMapper *m);

  // Description:
  // Return bounding box (array of six doubles) of data expressed as
  // (xmin,xmax, ymin,ymax, zmin,zmax).
  virtual double *GetBounds();
  virtual void GetBounds(double bounds[6])
    {this->Superclass::GetBounds(bounds);};

protected:  
  vtkStructuredGridMapper();
  ~vtkStructuredGridMapper() {};

  virtual int FillInputPortInformation(int, vtkInformation*);

private:
  vtkStructuredGridMapper(const vtkStructuredGridMapper&);  // Not implemented.
  void operator=(const vtkStructuredGridMapper&);  // Not implemented.
};

#endif
