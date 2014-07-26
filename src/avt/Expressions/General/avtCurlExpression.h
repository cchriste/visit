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
//                           avtCurlExpression.h                             //
// ************************************************************************* //

#ifndef AVT_CURL_FILTER_H
#define AVT_CURL_FILTER_H


#include <avtMacroExpressionFilter.h>


// ****************************************************************************
//  Class: avtCurlExpression
//
//  Purpose:
//      A filter that calculates the curl.  The curl takes in a vector and
//      produces a vector.  This depends on several partial derivatives, 
//      which are accomplished using the gradient expression.
//
//      Because we need to use other expressions, this is a derived type of
//      the macro expression filter.
//
//      curl of vector {u,v,w} = { grad(w)[1]-grad(v)[2], 
//                                 grad(u)[2]-grad(w)[0],
//                                 grad(v)[0]-grad(u)[1] }
//
//      Curl has the following physical interpretation:
//      Imagine you have a pinwheel -- a disc with some squares placed along
//      the disc in a direction orthogonal to the disc (so that air can spin
//      the disc by pushing on the squares).
//      If at some point (X,Y,Z) the curl is (a,b,c), then placing the disc so
//      that it is normal to (a,b,c) will give the fastest possible rotational 
//      speed that is attainable by having the center of the pinwheel at 
//      (X,Y,Z).
//
//      Also: John Boyd felt that we should define curl for 2D variables as 
//      well.  In this case, only the third component of the vector will be
//      non-zero, so we return a scalar (instead of a vector) in this case.
//
//  Programmer: Hank Childs
//  Creation:   December 27, 2004
//
//  Modifications:
//
//    Hank Childs, Fri Aug 19 08:50:02 PDT 2005
//    Move definition of GetVariableDimension to the .C file.
//
// ****************************************************************************

class EXPRESSION_API avtCurlExpression : public avtMacroExpressionFilter
{
  public:
                              avtCurlExpression();
    virtual                  ~avtCurlExpression();

    virtual const char       *GetType(void)   { return "avtCurlExpression"; };
    virtual const char       *GetDescription(void)
                               { return "Calculating Curl"; };

  protected:
    virtual int               GetVariableDimension();
    virtual void              GetMacro(std::vector<std::string> &, 
                                       std::string &, Expression::ExprType &);
};


#endif

