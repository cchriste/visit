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
//                            avtResampleFilter.h                            //
// ************************************************************************* //

#ifndef AVT_RESAMPLE_FILTER_H
#define AVT_RESAMPLE_FILTER_H

#include <filters_exports.h>

#include <avtDatasetToDatasetFilter.h>

#include <ResampleAttributes.h>


// ****************************************************************************
//  Class: avtResampleFilter
//
//  Purpose:
//      Resamples a dataset onto a rectilinear grid.
//
//  Programmer: Hank Childs
//  Creation:   March 26, 2001
//
//  Modifications:
//
//    Hank Childs, Fri Apr  6 17:39:40 PDT 2001
//    Added ability to bypass filter with already valid rectilinear grids.
//
//    Mark C. Miller, Tue Sep 13 20:09:49 PDT 2005
//    Added selID to support data selections
//
//    Hank Childs, Sat Apr 29 15:53:13 PDT 2006
//    Add argument to GetDimensions.
//
// ****************************************************************************

class AVTFILTERS_API avtResampleFilter : public avtDatasetToDatasetFilter
{
  public:
                          avtResampleFilter(const AttributeGroup*);
    virtual              ~avtResampleFilter();

    static avtFilter     *Create(const AttributeGroup*);

    virtual const char   *GetType(void)  { return "avtResampleFilter"; };
    virtual const char   *GetDescription(void) { return "Resampling"; };

  protected:
    ResampleAttributes    atts;
    char                 *primaryVariable;
    int                   selID;

    virtual void          Execute(void);
    virtual void          RefashionDataObjectInfo(void);

    void                  GetDimensions(int &, int &, int &, const double *,
                                        bool);
    bool                  InputNeedsNoResampling(void);
    void                  ResampleInput(void);
    void                  BypassResample(void);

    virtual int           AdditionalPipelineFilters(void) { return 2; };

    virtual avtPipelineSpecification_p
                          PerformRestriction(avtPipelineSpecification_p);
};


#endif


