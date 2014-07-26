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
//                     avtOriginalDataMinMaxQuery.h                          //
// ************************************************************************* //

#ifndef AVT_ORIGINAL_MINMAX_QUERY_H
#define AVT_ORIGINAL_MINMAX_QUERY_H
#include <query_exports.h>

#include <avtMinMaxQuery.h>

class avtExpressionEvaluatorFilter;


// ****************************************************************************
//  Class: avtOriginalDataMinMaxQuery
//
//  Purpose:
//    A query that retrieves the min/max of a var from the original data. 
//
//  Programmer: Kathleen Bonnell
//  Creation:   February 10, 2004 
//
//  Modifications:
//    Kathleen Bonnell, Wed Mar 31 16:07:50 PST 2004
//    Added optional constructor args.
//
//    Kathleen Bonnell, Wed Apr 14 18:05:08 PDT 2004 
//    Added condense filter. 
//
//    Kathleen Bonnell, Tue Jun 29 08:14:35 PDT 2004 
//    Removed condense filter. 
//
//    Hank Childs, Wed Dec 22 15:14:33 PST 2010
//    Add QuerySupportsTimeParallelization.
//
// ****************************************************************************

class QUERY_API avtOriginalDataMinMaxQuery : public avtMinMaxQuery
{
  public:
                                  avtOriginalDataMinMaxQuery(
                                      bool m = true, bool x = true);
    virtual                      ~avtOriginalDataMinMaxQuery();

    virtual bool                  OriginalData(void) { return true; };
    virtual bool                  QuerySupportsTimeParallelization(void)
                                         { return true; };

  protected:
    virtual avtDataObject_p       ApplyFilters(avtDataObject_p);   
    virtual int                   GetNFilters() { return 1; };

  private:
    avtExpressionEvaluatorFilter *eef;
};

#endif


