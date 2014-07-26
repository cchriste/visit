/*****************************************************************************
*
* Copyright (c) 2000 - 2014, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-442911
* All rights reserved.
*
* This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
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
*    documentation and/or other materials provided with the distribution.
*  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
*    be used to endorse or promote products derived from this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
* LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
* DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/

#ifndef AVT_STREAMLINE_RENDERER_IMPLEMENTATION_H
#define AVT_STREAMLINE_RENDERER_IMPLEMENTATION_H

#include <math.h>

class avtLookupTable;
class vtkPolyData;
class vtkDataArray;
class vtkCamera;
class StreamlineAttributes;
struct avtViewInfo;

// ****************************************************************************
//  Class:  avtStreamlineRendererImplementation
//
//  Purpose:
//    Implements the rendering-only portion of a molecule renderer in a
//    relatively stateless manner.  Meant to be instantiated at render
//    time by avtStreamlineRenderer::Render, though it can be kept around
//    across renderers while the implementation itself has not changed.
//
//  Programmer:  Jeremy Meredith
//  Creation:    February  3, 2006
//
//  Modifications:
//
//  Dave Pugmire, Fri Feb 12 14:02:57 EST 2010
//  Support for transparency sorting.
//
//  Hank Childs, Thu Sep 30 00:45:38 PDT 2010
//  Store the bounding box.
//
// ****************************************************************************
class avtStreamlineRendererImplementation
{
  public:
    avtStreamlineRendererImplementation(): varMin(0.0), varMax(0.0) {}
    virtual       ~avtStreamlineRendererImplementation() {}
    void           SetVarRange(const double &min, const double &max)
    {
        varMin = min;
        varMax = max;
    }
    void           SetBoundingBox(const double *b)
    {
        for (int i = 0 ; i < 6 ; i++)
            bbox[i] = b[i];
    }
    double         GetBBoxSize(void)
    {
        double vol = 1;
        int    numDims = 0;
        if (bbox[1] > bbox[0])
        {
            vol *= (bbox[1]-bbox[0]);
            numDims++;
        }
        if (bbox[3] > bbox[2])
        {
            vol *= (bbox[3]-bbox[2]);
            numDims++;
        }
        if (bbox[5] > bbox[4])
        {
            vol *= (bbox[5]-bbox[4]);
            numDims++;
        }

        double length = pow(vol, 1.0/numDims);
        return length;
    }
   
    virtual void   Render(vtkPolyData *data, const StreamlineAttributes&,
                          bool immediateModeRendering,
                          double vMin, double vMax,
                          vtkCamera *camera,
                          float ambient_coeff,
                          float spec_coeff, float spec_power,
                          float spec_r, float spec_g, float spec_b, 
                          const int *) = 0;
    virtual void   InvalidateColors() = 0;
    virtual void   SetLevelsLUT(avtLookupTable *) = 0;
    
protected:
    double varMin, varMax;
    double bbox[6];
};

#endif
