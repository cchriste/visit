/*****************************************************************************
*
* Copyright (c) 2000 - 2013, Lawrence Livermore National Security, LLC
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
//                              avtRayFunction.h                             //
// ************************************************************************* //

#ifndef AVT_RAY_FUNCTION_H
#define AVT_RAY_FUNCTION_H

#include <pipeline_exports.h>

#include <avtCellTypes.h>


class     avtGradients;
class     avtLightingModel;
class     avtRay;


// ****************************************************************************
//  Class: avtRayFunction
//
//  Purpose:
//      This is the base type for any ray function.  A ray function is a
//      routine that takes a series of sample points along a ray (avtRay) and
//      a lighting model and determines what the shading for the pixel should
//      be.
//
//  Programmer: Hank Childs
//  Creation:   December 1, 2000
//
//  Modifications:
//
//    Hank Childs, Sat Feb  3 20:14:19 PST 2001
//    Made the output of GetRayValue be a pixel instead of a value.
//
//    Hank Childs, Thu Feb  5 17:11:06 PST 2004
//    Moved inlined constructor and destructor definitions to .C files
//    because certain compilers have problems with them.
//
//    Hank Childs, Mon Sep 11 14:59:30 PDT 2006
//    Add method SetPrimaryIndex.  Also add methods for needing pixel
//    indices.
//
//    Hank Childs, Tue Dec 18 10:04:43 PST 2007
//    Define private copy constructor and assignment operator to prevent
//    accidental use of default, bitwise copy implementations.
//
//    Hank Childs, Sun Aug 31 08:04:42 PDT 2008
//    Remove infrastructure for gradients.  This is now done a different way.
//    Add support for lighting using gradients.
//
// ****************************************************************************

class PIPELINE_API avtRayFunction
{
  public:
                         avtRayFunction(avtLightingModel *);
    virtual             ~avtRayFunction();

    void                 SetPrimaryVariableIndex(int vi)
                                    { primaryVariableIndex = vi; };

    virtual void         GetRayValue(const avtRay *,
                                     unsigned char rgb[3], double) = 0;
    virtual bool         CanContributeToPicture(int,
                                          const double (*)[AVT_VARIABLE_LIMIT]);
    virtual bool         NeedPixelIndices(void) { return false; };

    void                 SetPixelIndex(int i, int j)
                             { pixelIndexI = i; pixelIndexJ = j; };

    void                 SetGradientVariableIndex(int gvi);

  protected:
    avtLightingModel    *lighting;
    int                  gradientVariableIndex;
    int                  primaryVariableIndex;
    int                  pixelIndexI, pixelIndexJ;

    inline int           IndexOfDepth(const double &, const int &);

  private:
    // These methods are defined to prevent accidental use of bitwise copy
    // implementations.  If you want to re-define them to do something
    // meaningful, that's fine.
                         avtRayFunction(const avtRayFunction &) {;};
    avtRayFunction      &operator=(const avtRayFunction &) { return *this; };
};


// ****************************************************************************
//  Method: avtRayFunction::IndexOfDepth
//
//  Purpose:
//      Determines the index of a depth in the z-buffer.  Assumes 0. is the
//      near plane, 1. is the far plane.
//
//  Programmer: Hank Childs
//  Creation:   February 13, 2001
//
// ****************************************************************************

inline int
avtRayFunction::IndexOfDepth(const double &depth, const int &numSamples)
{
    int rv = (int) (depth*numSamples);
    if (rv >= numSamples)
    {
        rv = numSamples-1;
    }
    else if (rv < 0)
    {
        rv = 0;
    }

    return rv;
}


#endif


