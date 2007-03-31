// ************************************************************************* //
//                            avtConstantCreatorFilter.h                     //
// ************************************************************************* //

#ifndef AVT_CONSTANT_CREATOR_FILTER_H
#define AVT_CONSTANT_CREATOR_FILTER_H

#include <avtUnaryMathFilter.h>

// ****************************************************************************
//  Class: avtConstantCreatorFilter
//
//  Purpose:
//      This class creates constant variables.
//
//  Programmer: Sean Ahern
//  Creation:   Fri Feb 21 23:28:13 America/Los_Angeles 2003
//
//  Modifications:
//
//    Hank Childs, Mon Nov  3 15:15:47 PST 2003
//    Override base class' definition of CreateArray.
//
// ****************************************************************************

class EXPRESSION_API avtConstantCreatorFilter : public avtUnaryMathFilter
{
  public:
                             avtConstantCreatorFilter() {value = 0;}
    virtual                 ~avtConstantCreatorFilter() {;}

    void                     SetValue(double v) {value = v;}
    virtual const char *     GetType(void) 
                                         { return "avtConstantCreatorFilter"; }
    virtual const char *     GetDescription(void) 
                                         { return "Generating constant"; }

  protected:
    virtual void             DoOperation(vtkDataArray *in, vtkDataArray *out,
                                         int ncomponents, int ntuples);
    virtual int              GetNumberOfComponentsInOutput(int) { return 1; };
    virtual vtkDataArray    *CreateArray(vtkDataArray *);

    double                   value;
};


#endif


