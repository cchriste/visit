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
//                             avtGradientExpression.h                       //
// ************************************************************************* //

#ifndef AVT_GRADIENT_FILTER_H
#define AVT_GRADIENT_FILTER_H

#include <avtSingleInputExpressionFilter.h>


class     vtkCell;
class     vtkDataArray;
class     vtkDataSet;
class     vtkIdList;
class     vtkRectilinearGrid;
class     vtkStructuredGrid;


// ****************************************************************************
//  Class: avtGradientExpression
//
//  Purpose:
//      A filter that calculates the gradient at each node. Note uses simple
//      definition of looking in the x,y, and z directions. 
//
//  Programmer: Akira Haddox
//  Creation:   July 30, 2002
//
//  Modifications:
//
//    Hank Childs, Sat Dec 13 10:42:15 PST 2003
//    Added support for rectilinear meshes.  Also don't force cell data to be
//    point data in the output.  Added ReleaseData to avoid memory bloat.
//
//    Hank Childs, Fri Mar  4 08:21:04 PST 2005
//    Removed data centering conversion modules.
//
//    Hank Childs, Mon Feb 13 14:45:18 PST 2006
//    Add support for logical gradients.  Also add perform restriction, so we
//    can request ghost zones.
//
//    Cyrus Harrison, Wed Aug  8 11:17:51 PDT 2007
//    Add support for multiple gradient algorithms.
//
//    Cyrus Harrison, Tue Apr  1 11:06:28 PDT 2008
//    Added IsPointVariable() to deal with NZQH centering change.
//
//    Hank Childs, Fri Jan  9 17:56:00 CST 2009
//    Added the approximate gradient option, to be used with ray casted
//    volume rendering.
//
//    Hank Childs, Sat Dec  4 11:30:22 PST 2010
//    Added a static method for calculating expressions.  This allows for other
//    places in VisIt to access a gradient without instantiating the filter.
//
// ****************************************************************************

typedef enum
{
    SAMPLE  =  0,
    LOGICAL , /* 1 */
    NODAL_TO_ZONAL_QUAD_HEX, /* 2 */
    FAST /* 3 */
} GradientAlgorithmType;


class EXPRESSION_API avtGradientExpression : public avtSingleInputExpressionFilter
{
  public:
                              avtGradientExpression();
    virtual                  ~avtGradientExpression();

    void                      SetAlgorithm(GradientAlgorithmType algo)
                               {gradientAlgo = algo;}

    virtual const char       *GetType(void)   { return "avtGradientExpression"; };
    virtual const char       *GetDescription(void)
                               { return "Calculating Gradient"; };

    virtual void              PreExecute(void);
    virtual void              ProcessArguments(ArgsExpr*, ExprPipelineState *);

    static vtkDataArray      *CalculateGradient(vtkDataSet *, const char *,
                                                GradientAlgorithmType=SAMPLE);

  protected:
    GradientAlgorithmType     gradientAlgo;

    virtual vtkDataArray     *DeriveVariable(vtkDataSet *, int currentDomainsIndex);
    virtual int               GetVariableDimension() { return 3; }
    virtual bool              IsPointVariable(void);
    
    virtual avtContract_p     ModifyContract(avtContract_p);

    static vtkDataArray      *RectilinearGradient(vtkRectilinearGrid *, 
                                                  const char *);
    static vtkDataArray      *LogicalGradient(vtkStructuredGrid *, 
                                              const char *);
    static vtkDataArray      *NodalToZonalQuadHexGrad(vtkStructuredGrid *, 
                                                      const char *);
    static vtkDataArray      *FastGradient(vtkDataSet *, const char *);
    static void               CalculateNodalToZonalQuadGrad(vtkDataSet *,
                                                            vtkDataArray *,
                                                            int ,
                                                            double *);
    static void               CalculateNodalToZonalHexGrad(vtkDataSet *,
                                                           vtkDataArray *,
                                                           int ,
                                                           double *);
    static double             EvaluateComponent(double, double, double, double,
                                                double, double, double,
                                                vtkDataSet *, vtkDataArray *,
                                                vtkIdList *);
    static double             EvaluateValue(double, double, double, 
                                            vtkDataSet *,
                                            vtkDataArray *,vtkIdList *,bool &);
};


#endif

