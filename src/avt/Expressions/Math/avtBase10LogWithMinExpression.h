/*****************************************************************************
*
* Copyright (c) 2000 - 2017, Lawrence Livermore National Security, LLC
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
//                          avtBase10LogWithMinExpression.h                         //
// ************************************************************************* //

#ifndef AVT_BASE10LOG_WITH_MIN_FILTER_H
#define AVT_BASE10LOG_WITH_MIN_FILTER_H


#include <avtBinaryMathExpression.h>

class     vtkDataArray;
class     ArgsExpr;
class     ExprPipelineState;


// ****************************************************************************
//  Class: avtBase10LogWithMinExpression
//
//  Purpose:
//      A filter that calculates the base 10 logarithm of its input and also
//      take a minimum value.  This expression can be done with other 
//      expressions, but is more efficiently implemented.
//
//  Programmer: Hank Childs
//  Creation:   May 20, 2010
//
//  Modifications:
//
//    Hank Childs, Thu Sep 23 14:10:45 PDT 2010
//    Indicate this expression can deal with singletons (improves performance).
//
// ****************************************************************************

class EXPRESSION_API avtBase10LogWithMinExpression 
                                               : public avtBinaryMathExpression
{
  public:
                           avtBase10LogWithMinExpression();
    virtual               ~avtBase10LogWithMinExpression();

    virtual const char    *GetType(void)   { return "avtBase10LogWithMinExpression"; };
    virtual const char    *GetDescription(void) 
                                { return "Calculating base 10 logarithm"; };

  protected:
    virtual void           DoOperation(vtkDataArray *in1, vtkDataArray *in2,
                                       vtkDataArray *out, int ncomps,int ntups);
    virtual bool           CanHandleSingleton(void) { return true; };
};


#endif


