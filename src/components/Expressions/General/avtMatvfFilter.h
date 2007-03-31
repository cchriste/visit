// ************************************************************************* //
//                              avtMatvfFilter.h                             //
// ************************************************************************* //

#ifndef AVT_MATVF_FILTER_H
#define AVT_MATVF_FILTER_H

#include <avtSingleInputExpressionFilter.h>

class     vtkDataArray;
class     ArgsExpr;
class     ExprPipelineState;
class     ConstExpr;

// ****************************************************************************
//  Class: avtMatvfFilter
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
//    Added PerformRestriction.  This is because matvf does not work with
//    ghost zone communication.  It cannot get the avtMaterial object with
//    ghost information and it causes an exception.  This will tell the
//    database that it cannot communicate ghost zones until a better solution
//    comes along.
//
// ****************************************************************************

class EXPRESSION_API avtMatvfFilter : public avtSingleInputExpressionFilter
{
  public:
                              avtMatvfFilter() {;};
    virtual                  ~avtMatvfFilter() {;};

    virtual const char       *GetType(void) { return "avtMatvfFilter"; };
    virtual const char       *GetDescription(void)
                                           {return "Calculating Material VF";};
    virtual void              ProcessArguments(ArgsExpr*, ExprPipelineState *);

  protected:
    virtual avtPipelineSpecification_p
                              PerformRestriction(avtPipelineSpecification_p);

    virtual vtkDataArray     *DeriveVariable(vtkDataSet *);
    virtual bool              IsPointVariable(void)  { return false; };

    void                      AddMaterial(ConstExpr *);
    std::vector<std::string>  matNames;
    std::vector<int>          matIndices;
};


#endif


