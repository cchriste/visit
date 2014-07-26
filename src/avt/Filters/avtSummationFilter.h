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
//                            avtSummationFilter.h                           //
// ************************************************************************* //

#ifndef AVT_SUMMATION_FILTER_H
#define AVT_SUMMATION_FILTER_H


#include <avtDataTreeIterator.h>
#include <filters_exports.h>
#include <string>


// ****************************************************************************
//  Class: avtSummationFilter
//
//  Purpose:
//      This should really be a query, not a filter.  It sums all of the values
//      for a variable.
//
//  Programmer: Hank Childs
//  Creation:   August 30, 2002
//
//  Modifications:
//
//    Hank Childs, Thu Feb  5 17:11:06 PST 2004
//    Moved inlined destructor definition to .C file because certain compilers
//    have problems with them.
//
// ****************************************************************************

class AVTFILTERS_API avtSummationFilter : public avtDataTreeIterator
{
  public:
                                    avtSummationFilter();
    virtual                        ~avtSummationFilter();

    void                            SetVariableName(std::string &);
    void                            SetSumType(std::string &);

    virtual const char             *GetType(void)
                                             { return "avtSummationFilter"; };
    virtual const char             *GetDescription(void)
                                             { return descriptionBuffer; };

    void                            PassData(bool);
    void                            SumGhostValues(bool);
    void                            IssueWarning(bool);

  protected:
    double                          sum;
    std::string                     variableName;
    std::string                     sumType;
    bool                            passData;
    bool                            sumGhostValues;
    bool                            issueWarning;
    char                            descriptionBuffer[1024];

    virtual void                    PreExecute(void);
    virtual void                    PostExecute(void);

    virtual vtkDataSet             *ExecuteData(vtkDataSet *, int,std::string);
};


#endif


