/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkOpenGLStructuredGridMapper.h,v $
  Language:  C++
  Date:      $Date: 2002/08/22 18:39:31 $
  Version:   $Revision: 1.31 $

  Copyright (c) 1993-2002 Ken Martin, Will Schroeder, Bill Lorensen 
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkOpenGLStructuredGridMapper - a StructuredGridMapper for the OpenGL library
// .SECTION Description
// vtkOpenGLStructuredGridMapper is a subclass of vtkStructuredGridMapper.
// vtkOpenGLStructuredGridMapper is a geometric StructuredGridMapper for the OpenGL 
// rendering library.

#ifndef __vtkOpenGLStructuredGridMapper_h
#define __vtkOpenGLStructuredGridMapper_h

#include "vtkStructuredGridMapper.h"
#include <plotter_exports.h>

class vtkProperty;
class vtkRenderWindow;
class vtkOpenGLRenderer;

class PLOTTER_API vtkOpenGLStructuredGridMapper : public vtkStructuredGridMapper
{
public:
  static vtkOpenGLStructuredGridMapper *New();
  vtkTypeMacro(vtkOpenGLStructuredGridMapper,vtkStructuredGridMapper);
  virtual void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Implement superclass render method.
  virtual void Render(vtkRenderer *ren, vtkActor *a);

  // Description:
  // Release any graphics resources that are being consumed by this mapper.
  // The parameter window could be used to determine which graphic
  // resources to release.
  void ReleaseGraphicsResources(vtkWindow *);

  // Description:
  // Draw method for OpenGL.
  virtual int Draw(vtkRenderer *ren, vtkActor *a);

protected:
  vtkOpenGLStructuredGridMapper();
  ~vtkOpenGLStructuredGridMapper();

  int ListStart;
  int CurrentList;
  int nLists;
  bool doingDisplayLists;
  int  primsInCurrentList;

  bool          ColorTexturingAllowed;
  bool          ColorTextureLoaded;
  unsigned int  ColorTextureName;
  float        *ColorTexture;
  int           ColorTextureSize;
  bool          OpenGLSupportsVersion1_2;
  double        LastOpacity;

  bool LooksDiscrete() const;
  bool MapScalarsWithTextureSupport(double);
  void BeginColorTexturing();
  void EndColorTexturing();
  bool UsesPointData(vtkDataSet *input, int scalarMode,
                     int arrayAccessMode, int arrayId, const char *arrayName,
                     int& offset);

private:
  vtkOpenGLStructuredGridMapper(const vtkOpenGLStructuredGridMapper&);  // Not implemented.
  void operator=(const vtkOpenGLStructuredGridMapper&);  // Not implemented.
};

#endif
