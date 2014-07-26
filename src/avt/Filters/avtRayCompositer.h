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

// ************************************************************************* //
//                              avtRayCompositer.h                           //
// ************************************************************************* //

#ifndef AVT_RAY_COMPOSITER_H
#define AVT_RAY_COMPOSITER_H

#include <filters_exports.h>

#include <avtSamplePointsToImageFilter.h>

#define BACKGROUND_SOLID           0
#define BACKGROUND_GRADIENT_TB     1
#define BACKGROUND_GRADIENT_BT     2
#define BACKGROUND_GRADIENT_LR     3
#define BACKGROUND_GRADIENT_RL     4
#define BACKGROUND_GRADIENT_RADIAL 5

class  avtRayFunction;
class  avtPixelizer;


// ****************************************************************************
//  Class: avtRayCompositer
//
//  Purpose:
//      Composites rays from sample points.  The output is an avtImage.
//
//  Programmer: Hank Childs
//  Creation:   December 4, 2000
//
//  Modifications:
//
//    Hank Childs, Sat Feb  3 20:17:20 PST 2001
//    Remove dependence on pixelizers.
//
//    Hank Childs, Tue Feb 13 15:15:50 PST 2001
//    Allowed for opaque images to be inserted in to volume rendering.
//
//    Brad Whitlock, Wed Dec 5 10:32:35 PDT 2001
//    Added methods to draw the background.
//
// ****************************************************************************

class AVTFILTERS_API avtRayCompositer : public avtSamplePointsToImageFilter
{
  public:
                          avtRayCompositer(avtRayFunction *);

    virtual const char   *GetType(void) { return "avtRayCompositer"; };
    virtual const char   *GetDescription(void) 
                                             { return "Compositing samples"; };
    void                  SetBackgroundColor(const unsigned char [3]);
    void                  SetBackgroundMode(int mode);
    void                  SetGradientBackgroundColors(const double [3],
                                                      const double [3]);
    void                  InsertOpaqueImage(avtImage_p);
    void                  UpdateCompositeProgress(int, int);

  protected:
    avtRayFunction       *rayfoo;
    int                   backgroundMode;
    unsigned char         background[3];
    double                gradBG1[3];
    double                gradBG2[3];
    avtImage_p            opaqueImage;

    virtual void          Execute(void);
    void                  FillBackground(unsigned char *, int, int);
    void                  DrawRadialGradient(unsigned char *, int, int);
};


#endif


