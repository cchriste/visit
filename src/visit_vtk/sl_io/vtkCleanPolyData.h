/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkCleanPolyData.h,v $
  Language:  C++
  Date:      $Date: 2002/09/03 12:52:22 $
  Version:   $Revision: 1.52 $

  Copyright (c) 1993-2002 Ken Martin, Will Schroeder, Bill Lorensen 
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkCleanPolyData - merge duplicate points, and/or remove unused points and/or remove degenerate cells
// .SECTION Description
// vtkCleanPolyData is a filter that takes polygonal data as input and
// generates polygonal data as output. vtkCleanPolyData can merge duplicate
// points (within specified tolerance and if enabled), eliminate points 
// that are not used, and if enabled, transform degenerate cells into 
// appropriate forms (for example, a triangle is converted into a line 
// if two points of triangle are merged). 
//
// Conversion of degenerate cells is controlled by the flags
// ConvertLinesToPoints, ConvertPolysToLines, ConvertStripsToPolys which act
// cumulatively such that a degenerate strip may become a poly.
// The full set is
// Line with 1 points -> Vert (if ConvertLinesToPoints)
// Poly with 2 points -> Line (if ConvertPolysToLines)
// Poly with 1 points -> Vert (if ConvertPolysToLines && ConvertLinesToPoints)
// Strp with 3 points -> Poly (if ConvertStripsToPolys)
// Strp with 2 points -> Line (if ConvertStripsToPolys && ConvertPolysToLines)
// Strp with 1 points -> Vert (if ConvertStripsToPolys && ConvertPolysToLines
//   && ConvertLinesToPoints)
//
// If tolerance is specified precisely=0.0, then vtkCleanPolyData will use
// the vtkMergePoints object to merge points (which is faster). Otherwise the
// slower vtkPointLocator is used.  Before inserting points into the point
// locator, this class calls a function OperateOnPoint which can be used (in
// subclasses) to further refine the cleaning process. See
// vtkQuantizePolyDataPoints.
//
// Note that merging of points can be disabled. In this case, a point locator
// will not be used, and points that are not used by any cells will be 
// eliminated, but never merged.

// .SECTION Caveats
// Merging points can alter topology, including introducing non-manifold
// forms. The tolerance should be chosen carefully to avoid these problems.
// Subclasses should handle OperateOnBounds as well as OperateOnPoint
// to ensure that the locator is correctly initialized (i.e. all modified 
// points must lie inside modified bounds).
//
// .SECTION See Also
// vtkQuantizePolyDataPoints

#ifndef __vtkCleanPolyData_h
#define __vtkCleanPolyData_h

#include "vtkPolyDataToPolyDataFilter.h"

class vtkPointLocator;

class VTK_GRAPHICS_EXPORT vtkCleanPolyData : public vtkPolyDataToPolyDataFilter
{
public:
  static vtkCleanPolyData *New();
  void PrintSelf(ostream& os, vtkIndent indent);
  vtkTypeRevisionMacro(vtkCleanPolyData,vtkPolyDataToPolyDataFilter);

  // Description:
  // By default ToleranceIsAbsolute is false and Tolerance is
  // a fraction of Bounding box diagonal, if true, AbsoluteTolerance is
  // used when adding points to locator (merging)
  vtkSetMacro(ToleranceIsAbsolute,int);
  vtkBooleanMacro(ToleranceIsAbsolute,int);
  vtkGetMacro(ToleranceIsAbsolute,int);

  // Description:
  // Specify tolerance in terms of fraction of bounding box length.
  vtkSetClampMacro(Tolerance,float,0.0,1.0);
  vtkGetMacro(Tolerance,float);

  // Description:
  // Specify tolerance in absolute terms
  vtkSetClampMacro(AbsoluteTolerance,float,0.0,VTK_LARGE_FLOAT);
  vtkGetMacro(AbsoluteTolerance,float);

  // Description:
  // Turn on/off conversion of degenerate lines to points
  vtkSetMacro(ConvertLinesToPoints,int);
  vtkBooleanMacro(ConvertLinesToPoints,int);
  vtkGetMacro(ConvertLinesToPoints,int);

  // Description:
  // Turn on/off conversion of degenerate polys to lines
  vtkSetMacro(ConvertPolysToLines,int);
  vtkBooleanMacro(ConvertPolysToLines,int);
  vtkGetMacro(ConvertPolysToLines,int);

  // Description:
  // Turn on/off conversion of degenerate strips to polys
  vtkSetMacro(ConvertStripsToPolys,int);
  vtkBooleanMacro(ConvertStripsToPolys,int);
  vtkGetMacro(ConvertStripsToPolys,int);

  // Description:
  // Set/Get a boolean value that controls whether point merging is
  // performed. If on, a locator will be used, and points laying within 
  // the appropriate tolerance may be merged. If off, points are never
  // merged. By default, merging is on.
  vtkSetMacro(PointMerging,int);
  vtkGetMacro(PointMerging,int);
  vtkBooleanMacro(PointMerging,int);

  // Description:
  // Set/Get a spatial locator for speeding the search process. By
  // default an instance of vtkMergePoints is used.
  void SetLocator(vtkPointLocator *locator);
  vtkGetObjectMacro(Locator,vtkPointLocator);

  // Description:
  // Create default locator. Used to create one when none is specified.
  void CreateDefaultLocator(void);

  // Description:
  // Release locator
  void ReleaseLocator(void);

  // Description:
  // Get the MTime of this object also considering the locator.
  unsigned long int GetMTime();

  // Description:
  // Perform operation on a point
  virtual void OperateOnPoint(float in[3], float out[3]);

  // Description:
  // Perform operation on bounds
  virtual void OperateOnBounds(float in[6], float out[6]);

  // This filter is difficult to stream.
  // To get invariant results, the whole input must be processed at once.
  // This flag allows the user to select whether strict piece invariance
  // is required.  By default it is on.  When off, the filter can stream,
  // but results may change.
  vtkSetMacro(PieceInvariant, int);
  vtkGetMacro(PieceInvariant, int);
  vtkBooleanMacro(PieceInvariant, int);

protected:
  vtkCleanPolyData();
 ~vtkCleanPolyData();

  // Usual data generation method
  void Execute();
  void ExecuteInformation();
  virtual void ComputeInputUpdateExtents(vtkDataObject *output);

  int   PointMerging;
  float Tolerance;
  float AbsoluteTolerance;
  int ConvertLinesToPoints;
  int ConvertPolysToLines;
  int ConvertStripsToPolys;
  int ToleranceIsAbsolute;
  vtkPointLocator *Locator;

  int PieceInvariant;
private:
  vtkCleanPolyData(const vtkCleanPolyData&);  // Not implemented.
  void operator=(const vtkCleanPolyData&);  // Not implemented.
};

#endif


