/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtk_expat.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#ifndef __vtk_expat_h
#define __vtk_expat_h

/* Use the expat library configured for VTK.  */
#include "vtkToolkits.h"
#ifdef VTK_USE_SYSTEM_EXPAT
# include <expat.h>
#else
# include <vtkexpat/expat.h>
#endif

#endif
