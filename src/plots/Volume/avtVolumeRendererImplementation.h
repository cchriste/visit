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

#ifndef AVT_VOLUME_RENDERER_IMPLEMENTATION_H
#define AVT_VOLUME_RENDERER_IMPLEMENTATION_H
#include <avtViewInfo.h>
#include <VolumeAttributes.h>
#include <vtkPointData.h>

class vtkDataSet;
class vtkDataArray;

// ****************************************************************************
//  Class:  avtVolumeRendererImplementation
//
//  Purpose:
//    Implements the rendering-only portion of a volume renderer in a
//    relatively stateless manner.  Meant to be instantiated at render
//    time by avtVolumeRenderer::Render, though it can be kept around
//    across renderers while the implementation itself has not changed.
//
//  Programmer:  Jeremy Meredith
//  Creation:    October  2, 2003
//
//  Modifications:
//    Brad Whitlock, Thu Jan 10 14:46:11 PST 2008
//    Added reducedDetail flag.
//
//    Brad Whitlock, Tue Apr 21 11:32:05 PDT 2009
//    I gathered properties into structs that we can pass to Render. I
//    also added the window size to the struct.
//
//    Jeremy Meredith, Tue Jan  5 15:53:30 EST 2010
//    Added storage of unnormalize gradient magnitude and its maximum value.
//
//    Allen Harvey, Thurs Nov 3 7:21:13 EST 2011
//    Added data items to support volume rendering without resampling
//
// ****************************************************************************

class avtVolumeRendererImplementation
{
  public:
    struct RenderProperties
    {
        RenderProperties() : view(), atts()
        {
            backgroundColor[0] = backgroundColor[1] = backgroundColor[2] = 0.;
            windowSize[0] = windowSize[1] = 0;
            reducedDetail = false;
        }

        avtViewInfo      view;
        VolumeAttributes atts;
        float            backgroundColor[3];
        int              windowSize[2];
        bool             reducedDetail;
    };

    struct VariableData
    {
        VariableData()
        {
            data = NULL;
            min = max = size = 0.f;
        }

        vtkDataArray *data;
        float         min;
        float         max;
        float         size;
    };

    struct VolumeData
    {
        VolumeData() : data(), opacity()
        {
            grid = NULL;
            gx = gy = gz = gm = gmn = hs = NULL;
            gm_max = 0;
        }

        VariableData        data;
        VariableData        opacity;
        vtkDataSet         *grid;
        float              *gx;
        float              *gy;
        float              *gz;
        float              *gm;
        float              *gmn;
        float               gm_max;
        float              *hs;
        float               hs_min;
    };

                   avtVolumeRendererImplementation() { }
    virtual       ~avtVolumeRendererImplementation() { }
    virtual void   Render(const RenderProperties &props, const VolumeData &volume) = 0;
};

#endif
