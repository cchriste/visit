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
//                         avtLabelSubsetsFilter.h                           //
// ************************************************************************* //

#ifndef AVT_LABEL_SUBSETS_FILTER_H
#define AVT_LABEL_SUBSETS_FILTER_H

#include <avtSIMODataTreeIterator.h>

#include <string>


// ****************************************************************************
//  Class: avtLabelSubsetsFilter
//
//  Purpose:  Ensures that the correct subset names are passed along
//            as labels.
//
//  Programmer: Kathleen Bonnell 
//  Creation:   October 16, 2001 
//
//  Modifications:
//    Brad Whitlock, Wed Aug 3 17:58:50 PST 2005
//    Copied from the Subset plot and made its only job to split up materials
//    and pass along other subset variables.
//
//    Eric Brugger, Tue Aug 19 10:51:44 PDT 2014
//    Modified the class to work with avtDataRepresentation.
//
// ****************************************************************************

class avtLabelSubsetsFilter : public avtSIMODataTreeIterator
{
  public:
                          avtLabelSubsetsFilter();
    virtual              ~avtLabelSubsetsFilter(){}; 

    virtual const char   *GetType(void) {return "avtLabelSubsetsFilter";};
    virtual const char   *GetDescription(void) 
                              { return "Setting subset names"; };

    void                  SetNeedMIR(bool val) { needMIR = val; };
  protected:
    virtual avtDataTree_p ExecuteDataTree(avtDataRepresentation *);
    virtual avtContract_p
                          ModifyContract(avtContract_p);

    virtual void          PostExecute(void);

    bool needMIR;
};


#endif
