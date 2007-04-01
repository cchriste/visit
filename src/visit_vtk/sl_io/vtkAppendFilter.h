/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkAppendFilter.h,v $
  Language:  C++
  Date:      $Date: 2002/11/03 15:57:42 $
  Version:   $Revision: 1.44 $

  Copyright (c) 1993-2002 Ken Martin, Will Schroeder, Bill Lorensen 
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkAppendFilter - appends one or more datasets together into a single unstructured grid
// .SECTION Description
// vtkAppendFilter is a filter that appends one of more datasets into a single
// unstructured grid. All geometry is extracted and appended, but point 
// attributes (i.e., scalars, vectors, normals, field data, etc.) are extracted 
// and appended only if all datasets have the point attributes available. 
// (For example, if one dataset has scalars but another does not, scalars will 
// not be appended.)

// .SECTION See Also
// vtkAppendPolyData

#ifndef __vtkAppendFilter_h
#define __vtkAppendFilter_h

#include "vtkDataSetToUnstructuredGridFilter.h"

class vtkDataSetCollection;

class VTK_GRAPHICS_EXPORT vtkAppendFilter : public vtkDataSetToUnstructuredGridFilter
{
public:
  static vtkAppendFilter *New();

  vtkTypeRevisionMacro(vtkAppendFilter,vtkDataSetToUnstructuredGridFilter);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Add a dataset to the list of data to append.
  void AddInput(vtkDataSet *in);

  // Description:
  // Get any input of this filter.
  vtkDataSet *GetInput(int idx);
  vtkDataSet *GetInput() 
    {return this->GetInput( 0 );}
  
  // Description:
  // Remove a dataset from the list of data to append.
  void RemoveInput(vtkDataSet *in);

  // Description:
  // Returns a copy of the input array.  Modifications to this list
  // will not be reflected in the actual inputs.
  vtkDataSetCollection *GetInputList();

protected:
  vtkAppendFilter();
  ~vtkAppendFilter();

  // Usual data generation method
  void Execute();

  // list of data sets to append together.
  // Here as a convenience.  It is a copy of the input array.
  vtkDataSetCollection *InputList;

private:
  // hide the superclass' AddInput() from the user and the compiler
  void AddInput(vtkDataObject *)
    { vtkErrorMacro( << "AddInput() must be called with a vtkDataSet not a vtkDataObject."); };
  void RemoveInput(vtkDataObject *input)
    { this->vtkProcessObject::RemoveInput(input); };
private:
  vtkAppendFilter(const vtkAppendFilter&);  // Not implemented.
  void operator=(const vtkAppendFilter&);  // Not implemented.
};


#endif


