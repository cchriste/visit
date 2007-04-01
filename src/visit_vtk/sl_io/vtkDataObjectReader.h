/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkDataObjectReader.h,v $
  Language:  C++
  Date:      $Date: 2002/05/31 23:12:41 $
  Version:   $Revision: 1.20 $

  Copyright (c) 1993-2002 Ken Martin, Will Schroeder, Bill Lorensen 
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkDataObjectReader - read vtk field data file
// .SECTION Description
// vtkDataObjectReader is a source object that reads ASCII or binary field
// data files in vtk format. Fields are general matrix structures used
// represent complex data. (See text for format details).  The output of this
// reader is a single vtkDataObject.  The superclass of this class,
// vtkDataReader, provides many methods for controlling the reading of the
// data file, see vtkDataReader for more information.
// .SECTION Caveats
// Binary files written on one system may not be readable on other systems.
// .SECTION See Also
// vtkFieldData vtkDataObjectWriter

#ifndef __vtkDataObjectReader_h
#define __vtkDataObjectReader_h

#include "vtkDataReader.h"
#include <vtk_sl_io_exports.h>

class vtkDataObject;

class VTK_SL_IO_API vtkDataObjectReader : public vtkDataReader
{
public:
  static vtkDataObjectReader *New();
  vtkTypeRevisionMacro(vtkDataObjectReader,vtkDataReader);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Get the output field of this reader.
  vtkDataObject *GetOutput();
  vtkDataObject *GetOutput(int idx)
    {return this->vtkSource::GetOutput(idx); };
  void SetOutput(vtkDataObject *);
  
protected:
  vtkDataObjectReader();
  ~vtkDataObjectReader();

  void Execute();
private:
  vtkDataObjectReader(const vtkDataObjectReader&);  // Not implemented.
  void operator=(const vtkDataObjectReader&);  // Not implemented.
};

#endif


