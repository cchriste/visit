// ************************************************************************* //
//                             avtAbsValFilter.h                             //
// ************************************************************************* //

#ifndef AVT_ABSVAL_FILTER_H
#define AVT_ABSVAL_FILTER_H

#include <avtUnaryMathFilter.h>

class     vtkDataArray;


// ****************************************************************************
//  Class: avtAbsValFilter
//
//  Purpose:
//      A filter that calculates the absolute value of its input.
//
//  Programmer: Sean Ahern
//  Creation:   Tue Jun 11 16:23:45 PDT 2002
//
//  Modifications:
//
//    Hank Childs, Fri Nov 15 15:23:54 PST 2002
//    Added support for vectors and non-float arrays.
//
//    Hank Childs, Thu Feb  5 17:11:06 PST 2004
//    Moved inlined constructor and destructor definitions to .C files 
//    because certain compilers have problems with them.
//
// ****************************************************************************

class EXPRESSION_API avtAbsValFilter : public avtUnaryMathFilter
{
  public:
                              avtAbsValFilter();
    virtual                  ~avtAbsValFilter();

    virtual const char       *GetType(void)   { return "avtAbsValFilter"; };
    virtual const char       *GetDescription(void) 
                                      { return "Calculating absolute value"; };

  protected:
    virtual void              DoOperation(vtkDataArray *in, vtkDataArray *out,
                                          int ncomponents, int ntuples);
};


#endif


