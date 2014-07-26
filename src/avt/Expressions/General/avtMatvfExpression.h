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
//                          avtMatvfExpression.h                             //
// ************************************************************************* //

#ifndef AVT_MATVF_FILTER_H
#define AVT_MATVF_FILTER_H

#include <avtSingleInputExpressionFilter.h>

class     vtkDataArray;
class     ArgsExpr;
class     ExprPipelineState;
class     ConstExpr;

// ****************************************************************************
//  Class: avtMatvfExpression
//
//  Purpose:
//      Creates the material fraction at each point.
//          
//  Programmer: Sean Ahern
//  Creation:   March 18, 2003
//
//  Modifications:
//    Jeremy Meredith, Mon Sep 29 12:13:04 PDT 2003
//    Added support for integer material indices.
//
//    Hank Childs, Fri Oct 24 14:49:23 PDT 2003
//    Added ModifyContract.  This is because matvf does not work with
//    ghost zone communication.  It cannot get the avtMaterial object with
//    ghost information and it causes an exception.  This will tell the
//    database that it cannot communicate ghost zones until a better solution
//    comes along.
//
//    Hank Childs, Thu Feb  5 17:11:06 PST 2004
//    Moved inlined constructor and destructor definitions to .C files
//    because certain compilers have problems with them.
//
//    Hank Childs, Wed Feb 18 09:15:23 PST 2004
//    Issue a warning if we encounter a bad material name.
//
//    Cyrus Harrison, Mon Jun 18 10:52:45 PDT 2007
//    Ensure that the output is always scalar (added GetVariableDimension)
//
//    Cyrus Harrison, Tue Feb 19 13:20:50 PST 2008
//    Added doPostGhost memeber. 
//
// ****************************************************************************

class EXPRESSION_API avtMatvfExpression : public avtSingleInputExpressionFilter
{
  public:
                              avtMatvfExpression();
    virtual                  ~avtMatvfExpression();

    virtual const char       *GetType(void) { return "avtMatvfExpression"; };
    virtual const char       *GetDescription(void)
                                           {return "Calculating Material VF";};
    virtual void              ProcessArguments(ArgsExpr*, ExprPipelineState *);


  protected:
    virtual int               GetVariableDimension(void) { return 1; };

    virtual avtContract_p
                              ModifyContract(avtContract_p);

    virtual vtkDataArray     *DeriveVariable(vtkDataSet *, int currentDomainsIndex);
    virtual bool              IsPointVariable(void)  { return false; };
    virtual void              PreExecute(void);

    void                      AddMaterial(ConstExpr *);

    bool                      doPostGhost;
    bool                      issuedWarning;
    std::vector<std::string>  matNames;
    std::vector<int>          matIndices;
};


#endif


