/*****************************************************************************
*
* Copyright (c) 2000 - 2006, The Regents of the University of California
* Produced at the Lawrence Livermore National Laboratory
* All rights reserved.
*
* This file is part of VisIt. For details, see http://www.llnl.gov/visit/. The
* full copyright notice is contained in the file COPYRIGHT located at the root
* of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
*
* Redistribution  and  use  in  source  and  binary  forms,  with  or  without
* modification, are permitted provided that the following conditions are met:
*
*  - Redistributions of  source code must  retain the above  copyright notice,
*    this list of conditions and the disclaimer below.
*  - Redistributions in binary form must reproduce the above copyright notice,
*    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
*    documentation and/or materials provided with the distribution.
*  - Neither the name of the UC/LLNL nor  the names of its contributors may be
*    used to  endorse or  promote products derived from  this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED.  IN  NO  EVENT  SHALL  THE  REGENTS  OF  THE  UNIVERSITY OF
* CALIFORNIA, THE U.S.  DEPARTMENT  OF  ENERGY OR CONTRIBUTORS BE  LIABLE  FOR
* ANY  DIRECT,  INDIRECT,  INCIDENTAL,  SPECIAL,  EXEMPLARY,  OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/

// ************************************************************************* //
//                     avtOpenGLSplattingVolumeRenderer.h                    //
// ************************************************************************* //

#ifndef AVT_OPEN_GL_SPLATTING_VOLUME_RENDERER_H
#define AVT_OPEN_GL_SPLATTING_VOLUME_RENDERER_H

#include <avtVolumeRendererImplementation.h>

// ****************************************************************************
//  Class: avtOpenGLSplattingVolumeRenderer
//
//  Purpose:
//      An implementation of a gaussian splat volume renderer that
//      uses OpenGL calls.
//
//  Programmer: Hank Childs
//  Creation:   April 24, 2002
//
//  Modifications:
//    Jeremy Meredith, Tue Sep 30 11:52:27 PDT 2003
//    Renamed to include the word "splatting".  Moved alphatex from the base
//    class to here.  Added OpenGL texture object ID.  Added method to
//    release texture object when we are done with it.
//
//    Jeremy Meredith, Thu Oct  2 13:36:28 PDT 2003
//    Made this class not inherit from avtVolumeRenderer.  It now 
//    gets most of its state from the arguments to Render.
//
// ****************************************************************************

class avtOpenGLSplattingVolumeRenderer : public avtVolumeRendererImplementation
{
  public:
                            avtOpenGLSplattingVolumeRenderer();
    virtual                ~avtOpenGLSplattingVolumeRenderer();

  protected:
    virtual void            Render(vtkRectilinearGrid *grid,
                                   vtkDataArray *data,
                                   vtkDataArray *opac,
                                   const avtViewInfo &view,
                                   const VolumeAttributes&,
                                   float vmin, float vmax, float vsize,
                                   float omin, float omax, float osize,
                                   float *gx, float *gy, float *gz,
                                   float *gmn);

    float                  *alphatex;
    unsigned int            alphatexId;
};


#endif


