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
//                        avtNeighborExpression.h                            //
// ************************************************************************* //

#ifndef AVT_NEIGHBOR_FILTER_H
#define AVT_NEIGHBOR_FILTER_H

#include <expression_exports.h>

#include <avtSingleInputExpressionFilter.h>

class     vtkDataArray;


// ****************************************************************************
//  Class: avtNeighborExpression
//
//  Purpose:
//    This is a filter that takes a mesh and decomposes it into a mesh of
//    vertexes. The points are also assigned a variable based on the distance
//    to the next nearest node.
//
//  Programmer: Akira Haddox
//  Creation:   June 27, 2002
//
//  Modifications:
//
//    Hank Childs, Fri Feb 20 15:51:54 PST 2004
//    Re-define GetVariableDimension.
//
//    Eric Brugger, Wed Aug 20 16:18:18 PDT 2014
//    Modified the class to work with avtDataRepresentation.
//
// ****************************************************************************

class EXPRESSION_API avtNeighborExpression : public avtSingleInputExpressionFilter
{
  public:
                             avtNeighborExpression();
    virtual                 ~avtNeighborExpression();

  protected:
    static const char        *variableName;

    virtual const char       *GetType(void)   { return "avtNeighborExpression"; };
    virtual const char       *GetDescription(void)
                             { return "Create vertex mesh,"
                                      " find distance to nearest node"; };

    // Used to fullfill parent's requirement, but unused since
    // ExecuteData exists for this derived class.
    virtual vtkDataArray     *DeriveVariable(vtkDataSet *, int currentDomainsIndex) { return NULL; }

    virtual avtDataRepresentation *ExecuteData(avtDataRepresentation *);
    virtual void             UpdateDataObjectInfo(void);

    virtual bool             IsPointVariable()     { return true; };
    virtual int              GetVariableDimension()   { return 1; };
};

#endif
