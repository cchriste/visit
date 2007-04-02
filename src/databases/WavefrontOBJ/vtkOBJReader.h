/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkOBJReader.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkOBJReader - read Wavefront .obj files
// .SECTION Description
// vtkOBJReader is a source object that reads Wavefront .obj
// files. The output of this source object is polygonal data.
// .SECTION See Also
// vtkOBJImporter

#ifndef __vtkOBJReader_h
#define __vtkOBJReader_h

#include "vtkPolyDataAlgorithm.h"

class vtkOBJReader : public vtkPolyDataAlgorithm 
{
public:
  static vtkOBJReader *New();
  vtkTypeRevisionMacro(vtkOBJReader,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Specify file name of Wavefront .obj file.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

protected:
  vtkOBJReader();
  ~vtkOBJReader();
  
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  char *FileName;
private:
  vtkOBJReader(const vtkOBJReader&);  // Not implemented.
  void operator=(const vtkOBJReader&);  // Not implemented.
};

#endif
