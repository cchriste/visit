// NOTE THAT THIS IS VTK's PLOT3D READER WITH OUR BUG FIXES AND ENHANCEMENT
// FOR READING IBLANKING.  THIS SHOULD ONLY BE INCLUDED IN VISIT'S
// REPOSITORY UNTIL THE KITWARE REPOSITORY HAS THESE CHANGES
/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkVisItPLOT3DReader.h,v $
  Language:  C++
  Date:      $Date: 2002/01/22 15:38:17 $
  Version:   $Revision: 1.48 $

  Copyright (c) 1993-2002 Ken Martin, Will Schroeder, Bill Lorensen 
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkVisItPLOT3DReader - read PLOT3D data files
// .SECTION Description
// vtkVisItPLOT3DReader is a reader object that reads PLOT3D formatted files and 
// generates a structured grid on output. PLOT3D is a computer graphics 
// program designed to visualize the grids and solutions of computational 
// fluid dynamics. Please see the "PLOT3D User's Manual" available from 
// NASA Ames Research Center, Moffett Field CA.
//
// PLOT3D files consist of a grid file (also known as XYZ file), an 
// optional solution file (also known as a Q file), and an optional function 
// file that contains user created data. The Q file contains solution 
// information as follows: the four parameters free stream mach number 
// (Fsmach), angle of attack (Alpha), Reynolds number (Re), and total 
// integration time (Time). In addition, the solution file contains 
// the flow density (scalar), flow momentum (vector), and flow energy (scalar).
//
// The reader can generate additional scalars and vectors (or "functions")
// from this information. To use vtkVisItPLOT3DReader, you must specify the 
// particular function number for the scalar and vector you want to visualize.
// This implementation of the reader provides the following functions. The
// scalar functions are:
//    -1  - don't read or compute any scalars
//    100 - density
//    110 - pressure
//    120 - temperature
//    130 - enthalpy
//    140 - internal energy
//    144 - kinetic energy
//    153 - velocity magnitude
//    163 - stagnation energy
//    170 - entropy
//    184 - swirl.
//
// The vector functions are:
//    -1  - don't read or compute any vectors
//    200 - velocity
//    201 - vorticity
//    202 - momentum
//    210 - pressure gradient.
//
// (Other functions are described in the PLOT3D spec, but only those listed are
// implemented here.) Note that by default, this reader creates the density 
// scalar (100) and momentum vector (202) as output. (These are just read in
// from the solution file.) Please note that the validity of computation is
// a function of this class's gas constants (R, Gamma) and the equations used.
// They may not be suitable for your computational domain.
//
// Additionally, you can read other data and associate it as a vtkDataArray
// into the output's point attribute data. Use the method AddFunction()
// to list all the functions that you'd like to read. AddFunction() accepts
// an integer parameter that defines the function number.
//
// The format of the function file is as follows. An integer indicating 
// number of grids, then an integer specifying number of functions per each 
// grid. This is followed by the (integer) dimensions of each grid in the 
// file. Finally, for each grid, and for each function, a float value per 
// each point in the current grid. Note: if both a function from the function
// file is specified, as well as a scalar from the solution file (or derived
// from the solution file), the function file takes precedence.

#ifndef __vtkVisItPLOT3DReader_h
#define __vtkVisItPLOT3DReader_h

#include <stdio.h>
#include <visit_vtk_exports.h>
#include "vtkStructuredGridSource.h"

class vtkIntArray;
class vtkFloatArray;
class vtkPointData;
class vtkPoints;
class vtkStructuredGrid;

// file formats
#define VTK_WHOLE_SINGLE_GRID_NO_IBLANKING 0
#define VTK_WHOLE_MULTI_GRID_NO_IBLANKING 2
#define VTK_WHOLE_MULTI_GRID_WITH_IBLANKING 3

class VISIT_VTK_API vtkVisItPLOT3DReader : public vtkStructuredGridSource 
{
public:
  static vtkVisItPLOT3DReader *New();
  vtkTypeRevisionMacro(vtkVisItPLOT3DReader,vtkStructuredGridSource);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Specify the PLOT3D file format to use
  vtkSetClampMacro(FileFormat,int,0,7);
  vtkGetMacro(FileFormat,int);

  // Description:
  // Set/Get the PLOT3D geometry FileName.
  vtkSetStringMacro(XYZFileName);
  vtkGetStringMacro(XYZFileName);

  // Description:
  // Set/Get the PLOT3D solution FileName.
  vtkSetStringMacro(QFileName);
  vtkGetStringMacro(QFileName);

  // Description:
  // Set/Get the PLOT3D function FileName.
  vtkSetStringMacro(FunctionFileName);
  vtkGetStringMacro(FunctionFileName);

  // Description:
  // Set/Get the PLOT3D vector FileName.
  vtkSetStringMacro(VectorFunctionFileName);
  vtkGetStringMacro(VectorFunctionFileName);

  // Description:
  // Specify the grid to read.
  vtkSetMacro(GridNumber,int);
  vtkGetMacro(GridNumber,int);

  // Description:
  // Specify the scalar function to extract. If ==(-1), then no scalar 
  // function is extracted.
  vtkSetMacro(ScalarFunctionNumber,int);
  vtkGetMacro(ScalarFunctionNumber,int);

  // Description:
  // Specify the vector function to extract. If ==(-1), then no vector
  // function is extracted.
  vtkSetMacro(VectorFunctionNumber,int);
  vtkGetMacro(VectorFunctionNumber,int);

  // Description:
  // Specify additional functions to read. These are placed into the
  // point data as data arrays. Later on they can be used by labeling
  // them as scalars, etc.
  void AddFunction(int functionNumber);
  void RemoveFunction(int);
  void RemoveAllFunctions();

  // these are read from PLOT3D file
  // Description:
  // Get the free-stream mach number.
  vtkGetMacro(Fsmach,float);

  // Description:
  // Get the angle of attack.
  vtkGetMacro(Alpha,float);

  // Description:
  // Get the Reynold's number.
  vtkGetMacro(Re,float);

  // Description:
  // Get the total integration time.
  vtkGetMacro(Time,float);

  // Description:
  // Set/Get the gas constant.
  vtkSetMacro(R,float);
  vtkGetMacro(R,float);

  // Description:
  // Set/Get the ratio of specific heats.
  vtkSetMacro(Gamma,float);
  vtkGetMacro(Gamma,float);

  // Description:
  // Set/Get the x-component of the free-stream velocity.
  vtkSetMacro(Uvinf,float);
  vtkGetMacro(Uvinf,float);

  // Description:
  // Set/Get the y-component of the free-stream velocity.
  vtkSetMacro(Vvinf,float);
  vtkGetMacro(Vvinf,float);

  // Description:
  // Set/Get the z-component of the free-stream velocity.
  vtkSetMacro(Wvinf,float);
  vtkGetMacro(Wvinf,float);

  // Description:
  // Get the number of grids. This is valid only after a
  // read has been performed.
  vtkGetMacro(NumberOfGrids, int);

protected:
  vtkVisItPLOT3DReader();
  ~vtkVisItPLOT3DReader();

  void ExecuteInformation();
  void Execute();
  int GetFileType(FILE *fp);

  //plot3d FileNames
  int FileFormat; //various PLOT3D formats
  char *XYZFileName;
  char *QFileName;
  char *FunctionFileName;
  char *VectorFunctionFileName;

  //flags describing data to be read
  int GridNumber; //for multi-grid files, the one we're interested in
  int ScalarFunctionNumber;
  int VectorFunctionNumber;
  int FunctionFileFunctionNumber;
  void MapFunction(int fNumber,vtkPointData *outputPD);

  //functions to read that are not scalars or vectors
  vtkIntArray *FunctionList;

  //temporary variables used during read
  float *TempStorage;
  int NumberOfPoints;
  int NumberOfGrids;

  //supplied in PLOT3D file
  float Fsmach;
  float Alpha;
  float Re;
  float Time;

  //parameters used in computing derived functions
  float R; 
  float Gamma;
  float Uvinf;
  float Vvinf;
  float Wvinf;
  
  //methods to read data
  int ReadBinaryGrid(FILE *fp,vtkStructuredGrid *output);
  int ReadBinaryGridDimensions(FILE *fp, vtkStructuredGrid *output);
  int ReadBinarySolution(FILE *fp, vtkStructuredGrid *output);
  int ReadBinaryFunctionFile(FILE *fp, vtkStructuredGrid *output);
  int ReadBinaryVectorFunctionFile(FILE *fp, vtkStructuredGrid *output);

  vtkPoints *Grid;
  vtkFloatArray *Density;
  vtkFloatArray *Energy;
  vtkFloatArray *Momentum;

  // derived functions from data in PLOT3D files
  void ComputeDensity(vtkPointData *outputPD);
  void ComputePressure(vtkPointData *outputPD);
  void ComputeTemperature(vtkPointData *outputPD);
  void ComputeEnthalpy(vtkPointData *outputPD);
  void ComputeInternalEnergy(vtkPointData *outputPD);
  void ComputeKineticEnergy(vtkPointData *outputPD);
  void ComputeVelocityMagnitude(vtkPointData *outputPD);
  void ComputeStagnationEnergy(vtkPointData *outputPD);
  void ComputeEntropy(vtkPointData *outputPD);
  void ComputeSwirl(vtkPointData *outputPD);

  void ComputeVelocity(vtkPointData *outputPD);
  void ComputeVorticity(vtkPointData *outputPD);
  void ComputeMomentum(vtkPointData *outputPD);
  void ComputePressureGradient(vtkPointData *outputPD);

private:
  vtkVisItPLOT3DReader(const vtkVisItPLOT3DReader&);  // Not implemented.
  void operator=(const vtkVisItPLOT3DReader&);  // Not implemented.
};

#endif


