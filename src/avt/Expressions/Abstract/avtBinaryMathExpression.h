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
//                             avtBinaryMathExpression.h                         //
// ************************************************************************* //

#ifndef AVT_BINARY_MATH_FILTER_H
#define AVT_BINARY_MATH_FILTER_H

#include <expression_exports.h>
#include <avtMultipleInputExpressionFilter.h>

class     vtkDataArray;


// ****************************************************************************
//  Class: avtBinaryMathExpression
//
//  Purpose:
//      A filter that performs a calculation on two variables.
//
//  Programmer: Sean Ahern
//  Creation:   Tue Jun 25 19:40:54 PDT 2002
//
//  Modifications:
//      Sean Ahern, Fri Jun 13 11:16:09 PDT 2003
//      Added the NumVariableArguments function, specifying that all
//      subclasses of BinaryMathFilter process two variable arguments.
//
//      Hank Childs, Thu Feb  5 17:11:06 PST 2004
//      Moved inlined constructor and destructor definitions to .C files
//      because certain compilers have problems with them.
//
//      Hank Childs, Thu Jan 20 17:30:09 PST 2005
//      Add data member cur_mesh so derived types can access mesh when 
//      deriving variables.
//
//      Hank Childs, Fri Aug 19 14:04:25 PDT 2005
//      Do a better job of determining the variable dimension.
//
//      Jeremy Meredith, Thu Feb 15 11:44:28 EST 2007
//      Added support for rectilinear grids with an inherent transform.
//      Binary math filters can handle these with no modifications.
//
//      Kathleen Biagas, Thu Apr 5 08:55:32 PDT 2012
//      Added second array to CreateArray method.
//
// ****************************************************************************

class EXPRESSION_API avtBinaryMathExpression 
    : public avtMultipleInputExpressionFilter
{
  public:
                              avtBinaryMathExpression();
    virtual                  ~avtBinaryMathExpression();

    virtual const char       *GetType(void)   { return "avtBinaryMathExpression";};
    virtual const char       *GetDescription(void) = 0;
    virtual int               NumVariableArguments() { return 2; }

  protected:
    virtual vtkDataArray     *DeriveVariable(vtkDataSet *, int currentDomainsIndex);
    virtual vtkDataArray     *CreateArray(vtkDataArray *, vtkDataArray *);
    virtual void              DoOperation(vtkDataArray *in1, vtkDataArray *in2,
                                          vtkDataArray *out, int, int) = 0;
    virtual int               GetNumberOfComponentsInOutput(int ncompsIn1,
                                                            int ncompsIn2)
                                         { return (ncompsIn1 > ncompsIn2
                                                    ? ncompsIn1 
                                                    : ncompsIn2);
                                         };
    virtual int               GetVariableDimension(void);

    virtual bool              FilterUnderstandsTransformedRectMesh();

    avtCentering              centering;
    vtkDataSet               *cur_mesh;
};


#endif


