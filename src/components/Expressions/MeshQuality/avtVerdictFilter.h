// ************************************************************************* //
//                             avtVerdictFilter.h                            //
// ************************************************************************* //

// Caveat: Verdict filters currently support triangles, but not triangle strips

#ifndef AVT_VERDICT_FILTER_H
#define AVT_VERDICT_FILTER_H

#include <expression_exports.h>

#include <avtSingleInputExpressionFilter.h>

class     vtkDataArray;


// ****************************************************************************
//  Class: avtVerdictFilter
//
//  Purpose:
//    This is a abstract base class for the verdict metric filters.
//
//  Programmer: Akira Haddox
//  Creation:   June 13, 2002
//
//  Modifications:
//
//    Hank Childs, Thu Feb  5 17:11:06 PST 2004
//    Moved inlined destructor definition to .C file because certain compilers
//    have problems with them.
//
//    Hank Childs, Thu May 19 10:44:29 PDT 2005
//    Add support for sub-types operating directly on the mesh.
//
// ****************************************************************************

class EXPRESSION_API avtVerdictFilter : public avtSingleInputExpressionFilter
{
  public:
                              avtVerdictFilter();
    virtual                  ~avtVerdictFilter();

    virtual void              PreExecute();

    virtual const char       *GetType(void)   { return "avtVerdictFilter"; };
    virtual const char       *GetDescription(void)
                                 { return "Calculating Verdict expression."; };

    virtual bool              OperateDirectlyOnMesh(vtkDataSet *)
                                       { return false; };
    virtual void              MetricForWholeMesh(vtkDataSet *, vtkDataArray *);

  protected:
    virtual vtkDataArray     *DeriveVariable(vtkDataSet *);

    virtual double            Metric(double coordinates[][3], int type) = 0;

    virtual bool              RequiresSizeCalculation() { return false; };

    virtual bool              IsPointVariable() {  return false; };
};

#endif
