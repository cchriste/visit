/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkRectilinearGridWriter.h,v $
  Language:  C++
  Date:      $Date: 2002/05/31 23:12:41 $
  Version:   $Revision: 1.22 $

  Copyright (c) 1993-2002 Ken Martin, Will Schroeder, Bill Lorensen 
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkRectilinearGridWriter - write vtk rectilinear grid data file
// .SECTION Description
// vtkRectilinearGridWriter is a source object that writes ASCII or binary 
// rectilinear grid data files in vtk format. See text for format details.

// .SECTION Caveats
// Binary files written on one system may not be readable on other systems.

#ifndef __vtkRectilinearGridWriter_h
#define __vtkRectilinearGridWriter_h

#include "vtkDataWriter.h"
#include <vtk_sl_io_exports.h>

class vtkRectilinearGrid;

class VTK_SL_IO_API vtkRectilinearGridWriter : public vtkDataWriter
{
public:
  static vtkRectilinearGridWriter *New();
  vtkTypeRevisionMacro(vtkRectilinearGridWriter,vtkDataWriter);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set / get the input data or filter.
  void SetInput(vtkRectilinearGrid *input);
  vtkRectilinearGrid *GetInput();
                               
protected:
  vtkRectilinearGridWriter() {};
  ~vtkRectilinearGridWriter() {};

  void WriteData();

private:
  vtkRectilinearGridWriter(const vtkRectilinearGridWriter&);  // Not implemented.
  void operator=(const vtkRectilinearGridWriter&);  // Not implemented.
};

#endif


