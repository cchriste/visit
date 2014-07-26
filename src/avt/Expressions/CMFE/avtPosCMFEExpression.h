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
//                            avtPosCMFEExpression.h                         //
// ************************************************************************* //

#ifndef AVT_POS_CMFE_EXPRESSION_H
#define AVT_POS_CMFE_EXPRESSION_H

#include <avtCMFEExpression.h>

#include <avtIntervalTree.h>

class     vtkCell;
class     vtkDataArray;
class     ArgsExpr;
class     ExprPipelineState;
class     ConstExpr;


// ****************************************************************************
//  Class: avtPosCMFEExpression
//
//  Purpose:
//      Does a position based cross-mesh field evaluation.
//          
//  Programmer: Hank Childs
//  Creation:   October 10, 2005
//
//  Modifications:
//    Jeremy Meredith, Thu Jan 18 11:04:51 EST 2007
//    Report explicitly that this filter does NOT understand transformed
//    rectilinear grids.  This method should default to returning false
//    anyway, but there are specific reasons this filter cannot yet be
//    optimized in this fashion, so ensure that even if other CMFE's change
//    to default to true, this one remains false until it can be fixed.
//
//    Hank Childs, Tue Mar 13 08:26:09 PDT 2012
//    Define virtual method OnlyRequiresSpatiallyOverlappingData.
//
// ****************************************************************************

class EXPRESSION_API avtPosCMFEExpression : public avtCMFEExpression
{
  public:
                              avtPosCMFEExpression();
    virtual                  ~avtPosCMFEExpression();

    virtual const char       *GetType(void){ return "avtPosCMFEExpression"; };
    virtual const char       *GetDescription(void)
                                           {return "Evaluating field";};
  protected:
    virtual avtDataTree_p     PerformCMFE(avtDataTree_p, avtDataTree_p,
                                          const std::string &,
                                          const std::string &);
    virtual bool              UseIdenticalSIL(void) { return false; };
    virtual bool              HasDefaultVariable(void) { return true; };
    virtual bool              FilterUnderstandsTransformedRectMesh();
    virtual bool              OnlyRequiresSpatiallyOverlappingData()
                                          { return true; };
};


#endif


