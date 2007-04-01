// ************************************************************************* //
//                            avtConditionalFilter.h                         //
// ************************************************************************* //

#ifndef AVT_CONDITIONAL_FILTER_H
#define AVT_CONDITIONAL_FILTER_H

#include <avtMultipleInputExpressionFilter.h>


// ****************************************************************************
//  Class: avtConditionalFilter
//
//  Purpose:
//      Performs an if-then-else test.  The first argument is the conditional
//      for the 'if' and the second two are the 'then' and 'else' clauses.
//
//  Programmer: Hank Childs
//  Creation:   August 20, 2003
//
//  Modifications:
//
//    Hank Childs, Thu Feb  5 17:11:06 PST 2004
//    Moved inlined constructor and destructor definitions to .C files
//    because certain compilers have problems with them.
//
// ****************************************************************************

class EXPRESSION_API avtConditionalFilter 
    : public avtMultipleInputExpressionFilter
{
  public:
                              avtConditionalFilter();
    virtual                  ~avtConditionalFilter();

    virtual const char       *GetType(void)  
                                    { return "avtConditionalFilter"; };
    virtual const char       *GetDescription(void)
                                    { return "Doing if-then-else"; };
    virtual int               NumVariableArguments() { return 3; }

  protected:
    virtual vtkDataArray     *DeriveVariable(vtkDataSet *);
};


#endif


