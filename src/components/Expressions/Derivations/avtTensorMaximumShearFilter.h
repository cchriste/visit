// ************************************************************************* //
//                         avtTensorMaximumShearFilter.h                     //
// ************************************************************************* //

#ifndef AVT_TENSOR_MAXIMUM_SHEAR_FILTER_H
#define AVT_TENSOR_MAXIMUM_SHEAR_FILTER_H


#include <avtUnaryMathFilter.h>


// ****************************************************************************
//  Class: avtTensorMaximumShearFilter
//
//  Purpose:
//      Calculates the maximum shear for a tensor.
//
//  Programmer: Hank Childs
//  Creation:   September 23, 2003
//
// ****************************************************************************

class EXPRESSION_API avtTensorMaximumShearFilter : public avtUnaryMathFilter
{
  public:
                               avtTensorMaximumShearFilter() {;};
    virtual                   ~avtTensorMaximumShearFilter() {;};

    virtual const char       *GetType(void)  
                                     { return "avtTensorMaximumShearFilter"; };
    virtual const char       *GetDescription(void)
                                     { return "Calculating maximum shear"; };

  protected:
    virtual void              DoOperation(vtkDataArray *in, vtkDataArray *out,
                                          int ncomponents, int ntuples);
    virtual int               GetNumberOfComponentsInOutput(int) { return 1; };
};

#endif


