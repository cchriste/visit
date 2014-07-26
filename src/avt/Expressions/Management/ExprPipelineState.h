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
//                           ExprPipelineState                               //
// ************************************************************************* //

#ifndef EXPR_PIPELINE_STATE_H
#define EXPR_PIPELINE_STATE_H

#include <avtDataObject.h>
#include <expression_exports.h>

class avtExpressionFilter;

// ****************************************************************************
//   Class: ExprPipelineState
//
//   Purpose:
//     Holds information about the pipeline state for expressions.
//
//  Programmer: Sean Ahern
//  Creation:   Thu Nov 21 15:15:07 PST 2002
//
//  Modifications:
//    Kathleen Bonnell, Thu Apr 22 14:42:38 PDT 2004
//    Moved code to new Source file.  Added ReleaseData method.
//
//    Hank Childs, Fri Dec 31 11:50:07 PST 2004
//    Add a Clear method.
//
//    Hank Childs, Fri Aug 25 17:26:59 PDT 2006
//    Add method GetNumNames.
//
// ****************************************************************************

class EXPRESSION_API ExprPipelineState
{
public:
                    ExprPipelineState();
                   ~ExprPipelineState();

    void            PushName(std::string s) {name_stack.push_back(s);} 
    std::string     PopName();
    int             GetNumNames(void) const { return static_cast<int>(name_stack.size()); };

    void            SetDataObject(avtDataObject_p d) {dataObject = d;}
    avtDataObject_p GetDataObject() {return dataObject;}
    void            AddFilter(avtExpressionFilter *f) {filters.push_back(f);}
    std::vector<avtExpressionFilter*>& GetFilters() {return filters;}

    void            ReleaseData(void);
    void            Clear();

protected:
    std::vector<std::string>    name_stack;
    avtDataObject_p             dataObject;
    std::vector<avtExpressionFilter*> filters;
};

#endif
