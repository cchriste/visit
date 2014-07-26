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
//                    avtWorldSpaceToImageSpaceTransform.h                   //
// ************************************************************************* //

#ifndef AVT_WORLD_SPACE_TO_IMAGE_SPACE_TRANSFORM_H
#define AVT_WORLD_SPACE_TO_IMAGE_SPACE_TRANSFORM_H

#include <filters_exports.h>

#include <avtTransform.h>

#include <avtViewInfo.h>

#include <vector>

class   avtIntervalTree;


// ****************************************************************************
//  Class: avtWorldSpaceToImageSpaceTransform
//
//  Purpose:
//      Transforms an avtDataset by a view matrix.
//
//  Note:     This class should probably redefine the CalcDomainList method and
//            have it cull away unused domains using a spatial extents interval
//            tree.
//
//  Programmer: Hank Childs
//  Creation:   November 27, 2000
//
//  Modifications:
//
//    Hank Childs, Mon Nov 26 18:33:16 PST 2001
//    Add support for aspect ratios.
//
//    Hank Childs, Fri Nov 19 13:38:21 PST 2004
//    Define ExecuteData so we can pass rectilinear grids through if specified.
//
//    Jeremy Meredith, Wed Jan 17 11:39:01 EST 2007
//    Added support for transformed rectilinear grids.
//
//    Hank Childs, Thu May 29 10:23:39 PDT 2008
//    Added argument to method for culling domains for the aspect ratio.
//
//    Hank Childs, Wed Dec 24 14:17:30 PST 2008
//    Make CalculateTransform be public and remove PreExecute, the majority
//    of whose logic went into avtRayTracer.
//
//    Hank Childs, Tue Sep 22 20:40:39 PDT 2009
//    Redefine virtual method in order to disable transformation of vectors.
//
//    Dave Pugmire, Fri May 14 08:04:43 EDT 2010
//    Move vector transform flag into base class.
//
//    Eric Brugger, Tue Jul 22 12:32:20 PDT 2014
//    Modified the class to work with avtDataRepresentation.
//
// ****************************************************************************

class AVTFILTERS_API avtWorldSpaceToImageSpaceTransform : public avtTransform
{
  public:
                       avtWorldSpaceToImageSpaceTransform(const avtViewInfo &,
                                                          double);
                       avtWorldSpaceToImageSpaceTransform(const avtViewInfo &,
                                                          const double *);
    virtual           ~avtWorldSpaceToImageSpaceTransform();

    virtual const char  *GetType(void)
                              { return "avtWorldSpaceToImageSpaceTransform"; };
    virtual const char  *GetDescription(void)
                              { return "Transforming data to image cube"; };

    static void        GetDomainsList(const avtViewInfo &, std::vector<int> &,
                                      const avtIntervalTree *, double aspect);

    void               TightenClippingPlanes(bool t)
                              { tightenClippingPlanes = t; };
    void               SetPassThruRectilinearGrids(bool t)
                              { passThruRectilinear = t; };
    static void        CalculateTransform(const avtViewInfo &,
                                       vtkMatrix4x4 *, const double *, double);

  protected:
    vtkMatrix4x4           *transform;
    avtViewInfo             view;
    double                  scale[3];
    double                  aspect;
    bool                    tightenClippingPlanes;
    bool                    passThruRectilinear;

    virtual vtkMatrix4x4   *GetTransform(void);
    virtual avtDataRepresentation *ExecuteData(avtDataRepresentation *);

    static void             CalculatePerspectiveTransform(const avtViewInfo &,
                                                          vtkMatrix4x4 *);
    static void             CalculateOrthographicTransform(const avtViewInfo &,
                                                           vtkMatrix4x4 *);
    virtual avtContract_p
                            ModifyContract(avtContract_p);

    virtual void            UpdateDataObjectInfo(void);
    virtual bool            FilterUnderstandsTransformedRectMesh();
};


#endif


