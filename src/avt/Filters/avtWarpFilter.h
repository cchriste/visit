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
//                             avtWarpFilter.h                               //
// ************************************************************************* //

#ifndef AVT_WARP_FILTER_H
#define AVT_WARP_FILTER_H

#include <filters_exports.h>

#include <avtDataTreeIterator.h>


// ****************************************************************************
//  Class: avtWarpFilter
//
//  Purpose:
//    A filter that 'warps' (elevates) a dataset by the point scalars.
//
//  Programmer: Kathleen Bonnell 
//  Creation:   July 12, 2006 
//
//  Modifications:
//    Kathleen Bonnell, Fri Nov 17 13:05:53 PST 2006
//    Needs FILTERS api for windows build.
//
//    Eric Brugger, Tue Jul 22 12:24:30 PDT 2014
//    Modified the class to work with avtDataRepresentation.
//
// ****************************************************************************

class AVTFILTERS_API avtWarpFilter : public avtDataTreeIterator
{
  public:
                              avtWarpFilter();
    virtual                  ~avtWarpFilter();

    virtual const char       *GetType(void)   { return "avtWarpFilter"; };
    virtual const char       *GetDescription(void)
                                  { return "Warping dataset"; };

  protected:
    virtual avtDataRepresentation *ExecuteData(avtDataRepresentation *);
    virtual void              PostExecute(void);
    virtual void              UpdateDataObjectInfo(void);
    virtual avtContract_p
                              ModifyContract(avtContract_p);
};


#endif


