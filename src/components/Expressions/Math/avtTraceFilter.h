// ************************************************************************* //
//                               avtTraceFilter.h                            //
// ************************************************************************* //

#ifndef AVT_TRACE_FILTER_H
#define AVT_TRACE_FILTER_H

#include <avtUnaryMathFilter.h>

class     vtkDataArray;


// ****************************************************************************
//  Class: avtTraceFilter
//
//  Purpose:
//      A filter that calculates trace (sum of diagonals of a tensor).
//
//  Programmer: Hank Childs
//  Creation:   September 19, 2003
//
// ****************************************************************************

class EXPRESSION_API avtTraceFilter : public avtUnaryMathFilter
{
  public:
                              avtTraceFilter() {;};
    virtual                  ~avtTraceFilter() {;};

    virtual const char       *GetType(void)   { return "avtTraceFilter"; };
    virtual const char       *GetDescription(void) 
                                              { return "Calculating trace"; };

  protected:
    virtual void              DoOperation(vtkDataArray *in, vtkDataArray *out,
                                          int ncomponents, int ntuples);
    virtual int               GetNumberOfComponentsInOutput(int)
                                          { return 1; };
};


#endif


