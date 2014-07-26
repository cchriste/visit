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
//                         avtLaplacianExpression.h                          //
// ************************************************************************* //

#ifndef AVT_LAPLACIAN_FILTER_H
#define AVT_LAPLACIAN_FILTER_H


#include <avtMacroExpressionFilter.h>


// ****************************************************************************
//  Class: avtLaplacianExpression
//
//  Purpose:
//      A filter that calculates the Laplacian.  The Laplacian takes in a 
//      scalar and produces a scalar.  The Laplacian is the divergence of
//      a gradient, which is why we need to be a derived type of the macro
//      expression filter.
//
//      Laplacian of scalar S  = divergence(gradient(S))
//
//      If F describes temperature, Laplacian(F) gives information about
//      the gain or loss of heat in a region.
//
//  Programmer: Hank Childs
//  Creation:   December 28, 2004
//
// ****************************************************************************

class EXPRESSION_API avtLaplacianExpression : public avtMacroExpressionFilter
{
  public:
                              avtLaplacianExpression();
    virtual                  ~avtLaplacianExpression();

    virtual const char       *GetType(void)  { return "avtLaplacianExpression"; };
    virtual const char       *GetDescription(void)
                               { return "Calculating the Laplacian"; };

  protected:
    virtual int               GetVariableDimension() { return 1; }
    virtual void              GetMacro(std::vector<std::string> &,
                                  std::string &ne, Expression::ExprType &type);
};


#endif

