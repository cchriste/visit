/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkVisItEnSightGoldBinaryReader.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkVisItEnSightGoldBinaryReader - class to read binary EnSight Gold files
// .SECTION Description
// vtkVisItEnSightGoldBinaryReader is a class to read EnSight Gold files into vtk.
// Because the different parts of the EnSight data can be of various data
// types, this reader produces multiple outputs, one per part in the input
// file.
// All variable information is being stored in field data.  The descriptions
// listed in the case file are used as the array names in the field data.
// For complex vector variables, the description is appended with _r (for the
// array of real values) and _i (for the array if imaginary values).  Complex
// scalar variables are stored as a single array with 2 components, real and
// imaginary, listed in that order.
// .SECTION Caveats
// You must manually call Update on this reader and then connect the rest
// of the pipeline because (due to the nature of the file format) it is
// not possible to know ahead of time how many outputs you will have or
// what types they will be.
// This reader can only handle static EnSight datasets (both static geometry
// and variables).

#ifndef __vtkVisItEnSightGoldBinaryReader_h
#define __vtkVisItEnSightGoldBinaryReader_h

#include "vtkVisItEnSightReader.h"

class vtkVisItEnSightGoldBinaryReader : public vtkVisItEnSightReader
{
public:
  static vtkVisItEnSightGoldBinaryReader *New();
  vtkTypeRevisionMacro(vtkVisItEnSightGoldBinaryReader, vtkVisItEnSightReader);
  virtual void PrintSelf(ostream& os, vtkIndent indent);
 
protected:
  vtkVisItEnSightGoldBinaryReader();
  ~vtkVisItEnSightGoldBinaryReader();
  
  // Returns 1 if successful.  Sets file size as a side action.
  int OpenFile(const char* filename);

  // Description:
  // Read the geometry file.  If an error occurred, 0 is returned; otherwise 1.
  virtual int ReadGeometryFile(const char* fileName, int timeStep);

  // Description:
  // Read the measured geometry file.  If an error occurred, 0 is returned;
  // otherwise 1.
  virtual int ReadMeasuredGeometryFile(const char* fileName, int timeStep);

  // Description:
  // Read scalars per node for this dataset.  If an error occurred, 0 is
  // returned; otherwise 1.  If there will be more than one component in
  // the data array, it is assumed that 0 is the first component added.
  virtual int ReadScalarsPerNode(const char* fileName, const char* description,
                                 int timeStep, int measured = 0,
                                 int numberOfComponents = 1,
                                 int component = 0);
  
  // Description:
  // Read vectors per node for this dataset.  If an error occurred, 0 is
  // returned; otherwise 1.
  virtual int ReadVectorsPerNode(const char* fileName, const char* description,
                                 int timeStep, int measured = 0);

  // Description:
  // Read tensors per node for this dataset.  If an error occurred, 0 is
  // returned; otherwise 1.
  virtual int ReadTensorsPerNode(const char* fileName, const char* description,
                                 int timeStep);

  // Description:
  // Read scalars per element for this dataset.  If an error occurred, 0 is
  // returned; otherwise 1.  If there will be more than one componenet in the
  // data array, it is assumed that 0 is the first component added.
  virtual int ReadScalarsPerElement(const char* fileName, const char* description,
                                    int timeStep, int numberOfComponents = 1,
                                    int component = 0);

  // Description:
  // Read vectors per element for this dataset.  If an error occurred, 0 is
  // returned; otherwise 1.
  virtual int ReadVectorsPerElement(const char* fileName, const char* description,
                                    int timeStep);

  // Description:
  // Read tensors per element for this dataset.  If an error occurred, 0 is
  // returned; otherwise 1.
  virtual int ReadTensorsPerElement(const char* fileName, const char* description,
                                    int timeStep);

  // Description:
  // Read an unstructured part (partId) from the geometry file and create a
  // vtkUnstructuredGrid output.  Return 0 if EOF reached. Return -1 if
  // an error occurred.
  virtual int CreateUnstructuredGridOutput(int partId, 
                                           char line[80],
                                           const char* name);
  
  // Description:
  // Read a structured part from the geometry file and create a
  // vtkStructuredGrid output.  Return 0 if EOF reached.
  virtual int CreateStructuredGridOutput(int partId, 
                                         char line[256],
                                         const char* name);
  
  // Description:
  // Read a structured part from the geometry file and create a
  // vtkRectilinearGrid output.  Return 0 if EOF reached.
  int CreateRectilinearGridOutput(int partId, char line[256], const char* name);
  
  // Description:
  // Read a structured part from the geometry file and create a
  // vtkImageData output.  Return 0 if EOF reached.
  int CreateImageDataOutput(int partId, char line[80], const char* name);
  
  // Description:
  // Internal function to read in a line up to 80 characters.
  // Returns zero if there was an error.
  int ReadLine(char result[80]);

  // Description:
  // Internal function to read in a single integer.
  // Returns zero if there was an error.
  int ReadInt(int *result);
  int ReadPartId(int *result);

  // Description:
  // Internal function to read in an integer array.
  // Returns zero if there was an error.
  int ReadIntArray(int *result, int numInts);

  // Description:
  // Internal function to read in a float array.
  // Returns zero if there was an error.
  int ReadFloatArray(float *result, int numFloats);

  // Description:
  // Read to the next time step in the geometry file.
  int SkipTimeStep();
  int SkipStructuredGrid(char line[256]);
  int SkipUnstructuredGrid(char line[256]);
  int SkipRectilinearGrid(char line[256]);
  int SkipImageData(char line[256]);
  
  int NodeIdsListed;
  int ElementIdsListed;
  
  ifstream *IFile;
  // The size of the file could be used to choose byte order.
  int FileSize;

private:
  vtkVisItEnSightGoldBinaryReader(const vtkVisItEnSightGoldBinaryReader&);  // Not implemented.
  void operator=(const vtkVisItEnSightGoldBinaryReader&);  // Not implemented.
};

#endif
