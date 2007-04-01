/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkStructuredPointsWriter.h,v $
  Language:  C++
  Date:      $Date: 2002/05/31 23:12:41 $
  Version:   $Revision: 1.42 $

  Copyright (c) 1993-2002 Ken Martin, Will Schroeder, Bill Lorensen 
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkStructuredPointsWriter - write vtk structured points data file
// .SECTION Description
// vtkStructuredPointsWriter is a source object that writes ASCII or binary 
// structured points data in vtk file format. See text for format details.
// .SECTION Caveats
// Binary files written on one system may not be writable on other systems.

#ifndef __vtkStructuredPointsWriter_h
#define __vtkStructuredPointsWriter_h

#include "vtkDataWriter.h"
#include <vtk_sl_io_exports.h>

class vtkImageData;

class VTK_SL_IO_API vtkStructuredPointsWriter : public vtkDataWriter
{
public:
  static vtkStructuredPointsWriter *New();
  vtkTypeRevisionMacro(vtkStructuredPointsWriter,vtkDataWriter);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set / get  the input data or filter.
  void SetInput(vtkImageData *input);
  vtkImageData *GetInput();
                               
protected:
  vtkStructuredPointsWriter() {};
  ~vtkStructuredPointsWriter() {};

  void WriteData();

private:
  vtkStructuredPointsWriter(const vtkStructuredPointsWriter&);  // Not implemented.
  void operator=(const vtkStructuredPointsWriter&);  // Not implemented.
};

#endif


