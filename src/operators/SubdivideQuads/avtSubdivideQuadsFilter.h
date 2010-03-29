/*****************************************************************************
*
* Copyright (c) 2000 - 2010, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-400124
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
//  File: avtSubdivideQuadsFilter.h
// ************************************************************************* //

#ifndef AVT_SubdivideQuads_FILTER_H
#define AVT_SubdivideQuads_FILTER_H


#include <avtPluginDataTreeIterator.h>
#include <SubdivideQuadsAttributes.h>


class vtkDataSet;


// ****************************************************************************
//  Class: avtSubdivideQuadsFilter
//
//  Purpose:
//      This operator will subdivide quadrilaterals based on a threshold
//      criteria.  The threshold is based on the maximum allowable change in 
//      variable value in a quad (so this operator really only makes sense for
//      nodal variables).  The point is to subdivide to the point that we 
//      aren't at the mercy of graphics hardware/rasterization techniques when
//      we are looking at only a handful of zones and those zones have a large
//      change in scalar value/color.
//
//  Note: the purpose of this operator is to subdivide quads using bilinear
//        interpolation because graphics rasterization does such a bad job.
//        Triangles are correctly interpolated, so this isn't a big deal.
//        That said, interpolation happens in color space, not data space.
//        So there is an option to subdivide triangles as well, because color
//        space interpolations can be deceiving.  Note that this option would
//        be totally worthless if we used textures instead.
//
//  Note: the edge connectivity of the surface after the application of this
//        operator will be poor (ie lots of T-intersections).
//        There is an option to take the points that are candidates to be
//        at T-intersections and move them away from the quad center, creating
//        overlapping quads.  This at least prevents artifacts where there are
//        "flashing holes" in your surface.  It can, however, lead to problems
//        where there is z-buffer fighting.  In practice, however, this isn't
//        too bad since the two quads fighting in z-buffer space typically
//        have very similar values.
//
//  Programmer: childs -- generated by xml2avt
//  Creation:   Tue Nov 2 05:54:25 PDT 2004
//
// ****************************************************************************

class avtSubdivideQuadsFilter : public avtPluginDataTreeIterator
{
  public:
                         avtSubdivideQuadsFilter();
    virtual             ~avtSubdivideQuadsFilter();

    static avtFilter    *Create();

    virtual const char  *GetType(void)  { return "avtSubdivideQuadsFilter"; };
    virtual const char  *GetDescription(void)
                             { return "SubdivideQuads"; };

    virtual void         SetAtts(const AttributeGroup*);
    virtual bool         Equivalent(const AttributeGroup*);

  protected:
    SubdivideQuadsAttributes   atts;
    bool                       haveIssuedWarning;

    virtual vtkDataSet   *ExecuteData(vtkDataSet *, int, std::string);
};


#endif
