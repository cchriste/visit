// ************************************************************************* //
//                           avtEffectiveTensorFilter.h                      //
// ************************************************************************* //

#ifndef AVT_EFFECTIVE_TENSOR_FILTER_H
#define AVT_EFFECTIVE_TENSOR_FILTER_H


#include <avtUnaryMathFilter.h>


// ****************************************************************************
//  Class: avtEffectiveTensorFilter
//
//  Purpose:
//      Calculates the principals of a tensor.
//
//  Programmer: Hank Childs
//  Creation:   September 23, 2003
//
// ****************************************************************************

class EXPRESSION_API avtEffectiveTensorFilter : public avtUnaryMathFilter
{
  public:
                               avtEffectiveTensorFilter() {;};
    virtual                   ~avtEffectiveTensorFilter() {;};

    virtual const char       *GetType(void)  
                               { return "avtEffectiveTensorFilter"; };
    virtual const char       *GetDescription(void)
                               {return "Calculating the effective tensor";};

  protected:
    virtual void              DoOperation(vtkDataArray *in, vtkDataArray *out,
                                          int ncomponents, int ntuples);
    virtual int               GetNumberOfComponentsInOutput(int) { return 1; };
};

#endif


