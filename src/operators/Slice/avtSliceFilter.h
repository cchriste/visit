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
//                            avtSliceFilter.h                               //
// ************************************************************************* //

#ifndef AVT_SLICE_FILTER_H
#define AVT_SLICE_FILTER_H


#include <avtPluginStreamer.h>
#include <SliceAttributes.h>

class vtkDataSet;
class vtkRectilinearGrid;
class vtkTransformFilter;
class vtkMatrix4x4;
class vtkSlicer;

class avtPointAttribute;


// ****************************************************************************
//  Class: avtSliceFilter
//
//  Purpose:
//      A filter that takes a slice of domains of an avtDataSet.
//
//  Programmer: Hank Childs
//  Creation:   July 25, 2000
//
//  Modifications:
//
//    Hank Childs, Tue Aug  8 16:52:28 PDT 2000
//    Made constructor take a plane instead of a cut-filter, since we want
//    this filter to be driven by a cut-plane and the cut-filter allowed
//    for other implicit functions besides planes.
//
//    Jeremy Meredith, Tue Sep 19 22:30:05 PDT 2000
//    Made constructor take raw origin and normal, added origin and
//    normal data members, and added Equivalent method.
//
//    Jeremy Meredith, Thu Sep 28 12:50:55 PDT 2000
//    Removed CreateOutputDatasets.  Changed interface to ExecuteDomain.
//
//    Jeremy Meredith, Thu Mar  1 13:29:27 PST 2001
//    Made attributes be stored as an SliceAttributes class.
//
//    Jeremy Meredith, Sun Mar  4 16:59:57 PST 2001
//    Added a static Create method.
//
//    Kathleen Bonnell, Tue Apr 10 10:49:10 PDT 2001 
//    Changed ExecuteDomain to ExecuteData. 
//
//    Hank Childs, Wed Jun  6 08:44:38 PDT 2001
//    Renamed some methods to fit changes in base class.
//
//    Hank Childs, Fri Mar 15 19:33:24 PST 2002
//    Made use of dynamic resolves points.
//
//    Hank Childs, Tue Aug  6 10:30:25 PDT 2002
//    Calculate the cells that intersect a rectilinear slice before slicing.
//
//    Kathleen Bonnell, Thu Apr 10 11:25:01 PDT 2003   
//    Added PostExecute method. Store inverse transform matrix for possible
//    use later in the pipeline (project-2d only). 
//
//    Jeremy Meredith, Mon May  5 14:31:45 PDT 2003
//    Removed "point" for now.  The slice window has changed, and dynamically
//    resolved attributes will work differently soon.
//
//    Kathleen Bonnell, Wed Jun  2 09:11:01 PDT 2004
//    Store transform matrix for possible use later in the pipeline. 
//
//    Hank Childs, Thu Jan 20 10:36:10 PST 2005
//    Added extra argument to ProjectExtents.
//
//    Hank Childs, Fri Aug 19 08:57:27 PDT 2005
//    Use vtkTransformFilter instead of vtkTransformPolyDataFilter, since 
//    vtkTransformPolyDataFilter does not pass on names of vectors, which
//    can screw us up down stream. ['6471]
//
// ****************************************************************************

class avtSliceFilter : public avtPluginStreamer
{
  public:
                            avtSliceFilter();
    virtual                ~avtSliceFilter();

    static avtFilter       *Create();

    virtual const char     *GetType(void) { return "avtSliceFilter"; };
    virtual const char     *GetDescription(void) { return "Slicing"; };
    virtual void            ReleaseData(void);

    virtual void            SetAtts(const AttributeGroup*);
    virtual bool            Equivalent(const AttributeGroup*);
    void                    ProjectExtents(const double *, double *);

  protected:
    SliceAttributes               atts;
    float                         D;
    double                        cachedOrigin[3];

    vtkSlicer                    *slicer;
    vtkTransformFilter           *transform;
    int                          *celllist;
    vtkMatrix4x4                 *invTrans;
    vtkMatrix4x4                 *origTrans;

    virtual avtPipelineSpecification_p
                            PerformRestriction(avtPipelineSpecification_p);
    virtual vtkDataSet     *ExecuteData(vtkDataSet *, int, std::string);
    virtual void            PreExecute(void);
    virtual void            PostExecute(void);

    virtual void            RefashionDataObjectInfo(void);

    void                    CalculateRectilinearCells(vtkRectilinearGrid *);
    void                    SetPlaneOrientation(double *);

    void                    GetOrigin(double &, double &, double &);
    void                    SetUpProjection(void);
};


#endif


