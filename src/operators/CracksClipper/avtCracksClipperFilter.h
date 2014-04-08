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
//  File: avtCracksClipperFilter.h
// ************************************************************************* //

#ifndef AVT_CracksClipper_FILTER_H
#define AVT_CracksClipper_FILTER_H


#include <avtPluginFilter.h>
#include <avtDatasetToDatasetFilter.h>
#include <avtDeferExpressionBaseFilter.h>
#include <CracksClipperAttributes.h>


class vtkDataSet;


// ****************************************************************************
//  Class: avtCracksClipperFilter
//
//  Purpose:
//    A plugin operator for clipping away Cracks.
//
//  Programmer: Kathleen Bonnell
//  Creation:   Thu Oct 13 08:17:36 PDT 2005
//
//  Modifications:
//    Kathleen Bonnell, Fri Oct 13 11:05:01 PDT 2006
//    Removed int arg from RemoveCracks method.
//
//    Kathleen Bonnell, Thu May  3 07:51:38 PDT 2007 
//    Changed inheritance for this filter, so that it can create a pipeline
//    with multiple filters.  Moved bulk of code to new avtRemoveCracksFilter.
//
//    Kathleen Bonnell, Wed Sep 29 08:58:10 PDT 2010
//    Add ivar 'calculateDensity'.
//
//    Kathleen Biagas, Tue Aug 14 15:05:23 MST 2012
//    Add ivar 'varname'.
//
// ****************************************************************************

class avtCracksClipperFilter : virtual public avtPluginFilter,
                               virtual public avtDatasetToDatasetFilter 
{
  public:
                         avtCracksClipperFilter();
    virtual             ~avtCracksClipperFilter();

    static avtFilter    *Create();

    virtual const char  *GetType(void)  { return "avtCracksClipperFilter"; };
    virtual const char  *GetDescription(void)
                             { return "CracksClipper"; };

    virtual void         SetAtts(const AttributeGroup*);
    virtual bool         Equivalent(const AttributeGroup*);

  protected:
    CracksClipperAttributes   atts;
    bool                  calculateDensity;
    std::string           varname;

    virtual void          Execute(void);
    virtual void          PostExecute(void);
    virtual void          UpdateDataObjectInfo(void);
    virtual int           AdditionalPipelineFilters(void);
    virtual avtContract_p ModifyContract(avtContract_p);
};


#endif
